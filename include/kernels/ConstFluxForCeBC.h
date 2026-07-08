#pragma once

#include "IntegratedBC.h"

class ConstFluxForCeBC : public IntegratedBC
{
public:
  static InputParameters validParams();
  ConstFluxForCeBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const Real & _I;
  const Real & _ChargeTime;
};
