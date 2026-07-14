#include "CentroidPostprocessorTransfer.h"

#include "FEProblemBase.h"
#include "MooseMesh.h"
#include "MooseVariableFieldBase.h"
#include "MultiApp.h"
#include "SystemBase.h"

#include "libmesh/mesh_base.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/system.h"

registerMooseObject("babblerApp", CentroidPostprocessorTransfer);

InputParameters
CentroidPostprocessorTransfer::validParams()
{
  InputParameters params = MultiAppTransfer::validParams();
  params.addClassDescription(
      "Writes each sub-app's postprocessor value into the parent element that "
      "owns the sub-app position (exact per-element map for CentroidMultiApp; "
      "target must be a CONSTANT MONOMIAL aux variable).");
  params.addRequiredParam<PostprocessorName>("postprocessor",
                                             "Postprocessor in the sub-apps to read");
  params.addRequiredParam<AuxVariableName>("variable",
                                           "CONSTANT MONOMIAL aux variable to write");
  params.addParam<bool>("aitken",
                        false,
                        "Apply scalar Aitken Delta-squared dynamic relaxation to the transferred "
                        "field, estimating the extrapolation factor from the last two residuals "
                        "instead of writing the incoming value directly.");
  return params;
}

CentroidPostprocessorTransfer::CentroidPostprocessorTransfer(const InputParameters & parameters)
  : MultiAppTransfer(parameters),
    _pp_name(getParam<PostprocessorName>("postprocessor")),
    _var_name(getParam<AuxVariableName>("variable")),
    _aitken(getParam<bool>("aitken")),
    _omega(1.0),
    _last_timestep(-1)
{
  if (_directions.size() != 1 || !_directions.contains(FROM_MULTIAPP))
    mooseError("CentroidPostprocessorTransfer only supports direction = from_multiapp");
}

void
CentroidPostprocessorTransfer::execute()
{
  FEProblemBase & problem = getFromMultiApp()->problemBase();
  MooseVariableFieldBase & var = problem.getVariable(0, _var_name);
  SystemBase & sys = var.sys();
  NumericVector<Number> & solution = sys.solution();

  MooseMesh & mesh = problem.mesh();
  std::unique_ptr<libMesh::PointLocatorBase> locator = mesh.getMesh().sub_point_locator();
  locator->enable_out_of_mesh_mode();

  const auto sys_num = sys.number();
  const auto var_num = var.number();

  if (_aitken)
  {
    const int current_step = problem.timeStep();
    if (current_step != _last_timestep)
    {
      _omega = 1.0;
      _prev_r.clear();
      _last_timestep = current_step;
    }
  }

  // Gather (position, value) pairs from ALL ranks: some MOOSE versions
  // distribute CentroidMultiApp sub-apps by index, not by element ownership,
  // so the app and its element can live on different ranks.
  std::vector<Real> packed;
  for (unsigned int i = 0; i < getFromMultiApp()->numGlobalApps(); ++i)
  {
    if (!getFromMultiApp()->hasLocalApp(i))
      continue;
    const Real value = getFromMultiApp()->appPostprocessorValue(i, _pp_name);
    const Point & pos = getFromMultiApp()->position(i);
    packed.push_back(pos(0));
    packed.push_back(pos(1));
    packed.push_back(pos(2));
    packed.push_back(value);
  }
  _communicator.allgather(packed, false);

  // First pass: locate owned elements and compute the raw update r = incoming - current.
  // For Aitken, also accumulate the dot products needed for the scalar relaxation factor
  // (a single global omega applied uniformly to the whole field, not a per-dof one, since
  // the transferred quantity is one coupled interface variable).
  std::unordered_map<dof_id_type, std::pair<Real, Real>> local_data; // dof -> (old_x, r)
  unsigned int n_nolocate = 0;
  Real local_r_dot_delta = 0.0;     // sum r_prev . (r - r_prev)
  Real local_delta_dot_delta = 0.0; // sum (r - r_prev) . (r - r_prev)

  for (std::size_t k = 0; k + 3 < packed.size(); k += 4)
  {
    const Point pos(packed[k], packed[k + 1], packed[k + 2]);
    const Real value = packed[k + 3];

    const Elem * elem = (*locator)(pos);
    if (!elem)
    {
      n_nolocate++;
      continue;
    }
    if (elem->processor_id() != processor_id())
      continue; // another rank owns it and will write it

    if (elem->n_dofs(sys_num, var_num) < 1)
      mooseError("CentroidPostprocessorTransfer: variable '", _var_name,
                 "' has no elemental dof — it must be CONSTANT MONOMIAL");

    const dof_id_type dof = elem->dof_number(sys_num, var_num, 0);
    const Real old_x = solution(dof);
    const Real r = value - old_x;
    local_data[dof] = std::make_pair(old_x, r);

    if (_aitken)
    {
      const auto it = _prev_r.find(dof);
      if (it != _prev_r.end())
      {
        const Real delta = r - it->second;
        local_r_dot_delta += it->second * delta;
        local_delta_dot_delta += delta * delta;
      }
    }
  }

  Real omega = 1.0;
  if (_aitken)
  {
    std::vector<Real> sums = {local_r_dot_delta, local_delta_dot_delta};
    _communicator.sum(sums);
    if (!_prev_r.empty() && sums[1] > std::numeric_limits<Real>::epsilon())
    {
      omega = -_omega * (sums[0] / sums[1]);
      omega = std::min(omega, 1e5);
      omega = std::max(omega, -10.0);
    }
    else
      omega = _omega;
    _omega = omega;
  }

  unsigned int n_written = 0;
  Real vmin = std::numeric_limits<Real>::max();
  Real vmax = -std::numeric_limits<Real>::max();

  for (auto & entry : local_data)
  {
    const dof_id_type dof = entry.first;
    const Real old_x = entry.second.first;
    const Real r = entry.second.second;
    const Real new_value = old_x + omega * r; // omega==1 reduces to writing the incoming value

    solution.set(dof, new_value);
    n_written++;
    vmin = std::min(vmin, new_value);
    vmax = std::max(vmax, new_value);

    if (_aitken)
      _prev_r[dof] = r;
  }

  solution.close();
  sys.update();

  std::cerr << "CentroidPostprocessorTransfer(" << name() << ") rank " << processor_id()
            << ": wrote " << n_written << " of " << packed.size() / 4 << " gathered, range ["
            << vmin << ", " << vmax << "], unlocated " << n_nolocate;
  if (_aitken)
    std::cerr << ", omega " << omega;
  std::cerr << std::endl;
}
