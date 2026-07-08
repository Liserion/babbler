#pragma once

#include "Kernel.h"

// Volume-averaged solid concentration ODE (variable u = cs_avg) in the
// Subramanian polynomial particle model:
//   d(cs_avg)/dt = -(3/Rs) * j_n,   j_n = (cs_avg - cs_surf) * 5*Ds/Rs
// This kernel adds the source term a*(u - cs)*5*Ds/Rs (a = 3/Rs); pair it
// with a TimeDerivative kernel on the same variable.
class CsAvgODEKernelPoly : public Kernel
{
public:
  static InputParameters validParams();
  CsAvgODEKernelPoly(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _couple_cs;
  unsigned int _couple_cs_var;

  const Real & _a;
  const Real & _Ds;
  const Real & _Rs;
};
