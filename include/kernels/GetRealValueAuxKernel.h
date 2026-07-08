#pragma once

#include "AuxKernel.h"

class GetRealValueAuxKernel : public AuxKernel
{
public:
  static InputParameters validParams();
  GetRealValueAuxKernel(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _input_dof;
  const Real & _FactorValue;
};
