#pragma once

#include "Kernel.h"

// Solid-phase potential equation in the cathode for the Subramanian
// polynomial particle model: Sigma*laplace(phis) + a*(1-eps)*j_n, with the
// substituted parabolic-closure flux j_n = (cs_avg - cs) * 5*Ds/Rs.
class CathodePhiSKernelPoly : public Kernel
{
public:
  static InputParameters validParams();
  CathodePhiSKernelPoly(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _couple_cs;
  const VariableValue & _cs_avg;
  unsigned int _couple_cs_var, _cs_avg_var;

  const Real & _Sigma;
  const Real & _eps;
  const Real & _a;
  const Real & _Ds;
  const Real & _Rs;
};
