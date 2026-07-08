#pragma once

#include "Kernel.h"

// Electrolyte concentration equation in the cathode for the Subramanian
// polynomial particle model. Identical transport terms to CathodeCeKernel,
// but the Butler-Volmer reaction rate is replaced by the substituted
// parabolic-closure flux  j_n = (cs_avg - cs) * 5*Ds/Rs  (exact at the
// solution of PolynomialCsClosureKernel), which is linear and bounded, so
// the exponential BV nonlinearity lives only in the closure equation.
class CathodeCeKernelPoly : public Kernel
{
public:
  static InputParameters validParams();
  CathodeCeKernelPoly(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _couple_cs;
  const VariableValue & _cs_avg;
  const VariableGradient & _grad_couple_phi2;
  unsigned int _couple_phi2_var, _couple_cs_var, _cs_avg_var;

  const Real & _D;
  const Real & _K;
  const Real & _eps;
  const Real & _Cm;
  const Real & _a;
  const Real & _Ds;
  const Real & _Rs;

  Real substitutedFlux() const; // a*(1-eps) * (cs_avg - cs) * 5*Ds/Rs
};
