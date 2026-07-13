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
  return params;
}

CentroidPostprocessorTransfer::CentroidPostprocessorTransfer(const InputParameters & parameters)
  : MultiAppTransfer(parameters),
    _pp_name(getParam<PostprocessorName>("postprocessor")),
    _var_name(getParam<AuxVariableName>("variable"))
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

  unsigned int n_written = 0, n_nolocate = 0, n_offproc = 0;
  Real vmin = std::numeric_limits<Real>::max();
  Real vmax = -std::numeric_limits<Real>::max();

  for (unsigned int i = 0; i < getFromMultiApp()->numGlobalApps(); ++i)
  {
    if (!getFromMultiApp()->hasLocalApp(i))
      continue;

    const Real value = getFromMultiApp()->appPostprocessorValue(i, _pp_name);
    const Point & pos = getFromMultiApp()->position(i);

    const Elem * elem = (*locator)(pos);
    if (!elem)
    {
      n_nolocate++;
      continue;
    }
    if (elem->processor_id() != processor_id())
    {
      n_offproc++;
      continue;
    }

    if (elem->n_dofs(sys_num, var_num) < 1)
      mooseError("CentroidPostprocessorTransfer: variable '", _var_name,
                 "' has no elemental dof — it must be CONSTANT MONOMIAL");

    const dof_id_type dof = elem->dof_number(sys_num, var_num, 0);
    solution.set(dof, value);
    n_written++;
    vmin = std::min(vmin, value);
    vmax = std::max(vmax, value);
  }

  solution.close();
  sys.update();

  _console << "CentroidPostprocessorTransfer(" << name() << "): wrote " << n_written
           << " values [" << vmin << ", " << vmax << "], skipped " << n_nolocate
           << " unlocated / " << n_offproc << " off-proc (rank " << processor_id() << ")"
           << std::endl;
}
