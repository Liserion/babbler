#include "CsAvgODEKernelPoly.h"

registerMooseObject("babblerApp", CsAvgODEKernelPoly);

InputParameters
CsAvgODEKernelPoly::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addParam<Real>("a", 3.0, "Area/Volume (= 3/Rs)");
  params.addParam<Real>("Ds", 1.0e-5, "solid diffusion coefficient");
  params.addParam<Real>("Rs", 0.5, "particle radius");
  params.addRequiredCoupledVar("Cs", "surface concentration (closure variable)");
  return params;
}

CsAvgODEKernelPoly::CsAvgODEKernelPoly(const InputParameters & parameters)
  : Kernel(parameters),
    _couple_cs(coupledValue("Cs")),
    _couple_cs_var(coupled("Cs")),
    _a(getParam<Real>("a")),
    _Ds(getParam<Real>("Ds")),
    _Rs(getParam<Real>("Rs"))
{
}

Real
CsAvgODEKernelPoly::computeQpResidual()
{
  // residual: d(cs_avg)/dt + a * j_n = 0 with j_n = (u - cs) * 5*Ds/Rs
  return _a * (_u[_qp] - _couple_cs[_qp]) * 5.0 * _Ds / _Rs * _test[_i][_qp];
}

Real
CsAvgODEKernelPoly::computeQpJacobian()
{
  return _a * 5.0 * _Ds / _Rs * _phi[_j][_qp] * _test[_i][_qp];
}

Real
CsAvgODEKernelPoly::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _couple_cs_var)
    return -_a * 5.0 * _Ds / _Rs * _phi[_j][_qp] * _test[_i][_qp];

  return 0.0;
}
