#pragma once

#include "MultiAppTransfer.h"

#include <unordered_map>

// Exact backward map for CentroidMultiApp: writes each sub-app's
// postprocessor value into the element that owns the sub-app's position
// (one value per element -> use a CONSTANT MONOMIAL aux variable).
// Replaces MultiAppPostprocessorInterpolationTransfer, whose num_points=1
// mode produces corrupted fields in some MOOSE versions and whose
// num_points>1 mode blends neighboring sub-apps, making the macro and
// micro problems inconsistent (Picard cannot contract).
//
// Optional scalar Aitken Delta-squared relaxation (aitken=true): the plain
// fixed-point map on this macro<->micro coupling has spectral radius very
// close to 1 (flat LFP OCV plateau), so it contracts correctly but far too
// slowly for a handful of Picard sweeps per step. Aitken estimates the
// extrapolation factor a near-singular map actually needs from the last two
// residuals, instead of guessing a fixed relaxation/acceleration a priori
// (which is what MOOSE's built-in secant does, unsafeguarded, and why it
// overshoots into the BV exponential and diverges on this problem).
class CentroidPostprocessorTransfer : public MultiAppTransfer
{
public:
  static InputParameters validParams();
  CentroidPostprocessorTransfer(const InputParameters & parameters);

  virtual void execute() override;

protected:
  const PostprocessorName _pp_name;
  const AuxVariableName _var_name;
  const bool _aitken;

  Real _omega;
  int _last_timestep;
  std::unordered_map<dof_id_type, Real> _prev_r;
};
