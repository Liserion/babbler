#pragma once

#include "MultiAppTransfer.h"

// Exact backward map for CentroidMultiApp: writes each sub-app's
// postprocessor value into the element that owns the sub-app's position
// (one value per element -> use a CONSTANT MONOMIAL aux variable).
// Replaces MultiAppPostprocessorInterpolationTransfer, whose num_points=1
// mode produces corrupted fields in some MOOSE versions and whose
// num_points>1 mode blends neighboring sub-apps, making the macro and
// micro problems inconsistent (Picard cannot contract).
class CentroidPostprocessorTransfer : public MultiAppTransfer
{
public:
  static InputParameters validParams();
  CentroidPostprocessorTransfer(const InputParameters & parameters);

  virtual void execute() override;

protected:
  const PostprocessorName _pp_name;
  const AuxVariableName _var_name;
};
