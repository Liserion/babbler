#pragma once

#include "Kernel.h"

// Algebraic closure equation for the surface concentration (variable u = cs_surf)
// in the Subramanian (2005) two-parameter polynomial particle model:
//   (u - cs_avg)/beta + j_n(ce, phis, phie, u) = 0,   beta = Rs/(5*Ds)
// Solved monolithically with the macro variables (no lagged AuxVariable).
class PolynomialCsClosureKernel : public Kernel
{
public:
  static InputParameters validParams();
  PolynomialCsClosureKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _cs_avg;
  const VariableValue & _couple_c;
  const VariableValue & _couple_phi1;
  const VariableValue & _couple_phi2;
  unsigned int _cs_avg_var, _couple_c_var, _couple_phi1_var, _couple_phi2_var;

  const Real & _Cm;
  const Real & _K2;
  const int & _MateChoice;
  const Real _T;
  const Real & _Ds;
  const Real & _Rs;
  const Real & _Omega;

  const VariableValue & _couple_damage;
  const VariableValue & _couple_sigmaH;

private:
  void OpenCircuitV(const int & matechoice, Real x, Real & U, Real & dUdx);
  // per-area BV flux and its partial derivatives at the current qp
  void BV(const Real & c, const Real & phi1, const Real & phi2, const Real & cs,
          Real & J, Real & dJdcs, Real & dJdc, Real & dJdphi1, Real & dJdphi2);
};
