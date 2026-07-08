#pragma once

#include "Kernel.h"

// Electrolyte potential equation in the cathode for the Subramanian
// polynomial particle model. Identical transport terms to CathodePhiEKernel,
// with the substituted parabolic-closure flux j_n = (cs_avg - cs) * 5*Ds/Rs.
class CathodePhiEKernelPoly : public Kernel
{
public:
  static InputParameters validParams();
  CathodePhiEKernelPoly(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _couple_c;
  const VariableValue & _couple_cs;
  const VariableValue & _cs_avg;
  const VariableGradient & _grad_couple_c;
  unsigned int _couple_c_var, _couple_cs_var, _cs_avg_var;

  const Real & _K;
  const Real & _eps;
  const Real & _Cm;
  const Real & _a;
  const Real & _Ds;
  const Real & _Rs;
};
