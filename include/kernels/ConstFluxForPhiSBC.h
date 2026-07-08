#pragma once

#include "IntegratedBC.h"

class ConstFluxForPhiSBC : public IntegratedBC
{
public:
  static InputParameters validParams();
  ConstFluxForPhiSBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override { return 0.0; }

  const Real & _I;
  const Real & _ChargeTime;
};
