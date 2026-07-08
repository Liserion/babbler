#pragma once

#include "Kernel.h"

class CathodePhiSKernel : public Kernel
{
public:
  static InputParameters validParams();
  CathodePhiSKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _couple_c;
  const VariableValue & _couple_phi2;
  const VariableValue & _couple_cs;
  unsigned int _couple_c_var;
  unsigned int _couple_phi2_var;

  const Real & _Sigma;
  const Real & _eps;
  const Real & _K2;
  const Real & _Cm;
  const Real & _a;
  const int & _MateChoice;
  const Real & _T;

  const VariableValue & _couple_damage;
  const VariableValue & _couple_sigmaH;
  const Real & _Omega;

private:
  Real Jeff, dJdc, dJdphi1, dJdphi2;

  Real OpenCircuitV(const int & matechoice, Real x);
  void BV(const Real & c, const Real & phi1, const Real & phi2, const Real & cs,
          Real & JEFF, Real & DJDC, Real & DJDPHI1, Real & DJDPHI2);
};
