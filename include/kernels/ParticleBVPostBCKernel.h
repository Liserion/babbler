#pragma once

#include "IntegratedBC.h"

class ParticleBVPostBCKernel : public IntegratedBC
{
public:
  static InputParameters validParams();
  ParticleBVPostBCKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const PostprocessorName & _pps_c2;
  const PostprocessorName & _pps_phi1;
  const PostprocessorName & _pps_phi2;

  const Real & _K2;
  const Real & _Cm;
  const Real & _T;
  const int & _MateChoice;

  Real J, dJdc;

  Real Sech(const Real & x)
  {
    return 2.0 / (exp(x) + exp(-x));
  }

  void OpenCircuitV(const Real & x, Real & u, Real & dudx);
  void BV(const Real & c, const Real & phi1, const Real & phi2,
          const Real & cs, Real & J, Real & dJdc);
};
