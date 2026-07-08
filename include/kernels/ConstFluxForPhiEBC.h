#pragma once

#include "IntegratedBC.h"

class ConstFluxForPhiEBC : public IntegratedBC
{
public:
  static InputParameters validParams();
  ConstFluxForPhiEBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override { return 0.0; }

  const Real & _I;
  const Real & _ChargeTime;
};
