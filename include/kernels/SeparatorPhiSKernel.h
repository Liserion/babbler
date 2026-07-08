#pragma once

#include "Kernel.h"

class SeparatorPhiSKernel : public Kernel
{
public:
  static InputParameters validParams();
  SeparatorPhiSKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const Real & _Sigma;
};
