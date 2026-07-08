#pragma once

#include "Kernel.h"

class SeparatorCeKernel : public Kernel
{
public:
  static InputParameters validParams();
  SeparatorCeKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableGradient & _grad_couple_phi2;
  unsigned int _couple_phi2_var;

  const Real & _D;
  const Real & _K;
  const Real & _eps;

private:
  Real Deff, Keff;
};
