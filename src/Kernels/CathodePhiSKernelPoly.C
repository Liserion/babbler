#include "CathodePhiSKernelPoly.h"

registerMooseObject("babblerApp", CathodePhiSKernelPoly);

InputParameters
CathodePhiSKernelPoly::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<Real>("Sigma", "conductivity of solid phase");
  params.addRequiredParam<Real>("eps", "porosity");
  params.addParam<Real>("a", 3.0, "Area/Volume");
  params.addParam<Real>("Ds", 1.0e-5, "solid diffusion coefficient");
  params.addParam<Real>("Rs", 0.5, "particle radius");
  params.addRequiredCoupledVar("Cs", "surface concentration (closure variable)");
  params.addRequiredCoupledVar("CsAvg", "volume-averaged solid concentration");
  return params;
}

CathodePhiSKernelPoly::CathodePhiSKernelPoly(const InputParameters & parameters)
  : Kernel(parameters),
    _couple_cs(coupledValue("Cs")),
    _cs_avg(coupledValue("CsAvg")),
    _couple_cs_var(coupled("Cs")),
    _cs_avg_var(coupled("CsAvg")),
    _Sigma(getParam<Real>("Sigma")),
    _eps(getParam<Real>("eps")),
    _a(getParam<Real>("a")),
    _Ds(getParam<Real>("Ds")),
    _Rs(getParam<Real>("Rs"))
{
}

Real
CathodePhiSKernelPoly::computeQpResidual()
{
  const Real Jeff = _a * (1 - _eps) * (_cs_avg[_qp] - _couple_cs[_qp]) * 5.0 * _Ds / _Rs;

  return _Sigma * _grad_u[_qp] * _grad_test[_i][_qp]
       + Jeff * _test[_i][_qp];
}

Real
CathodePhiSKernelPoly::computeQpJacobian()
{
  return _Sigma * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

Real
CathodePhiSKernelPoly::computeQpOffDiagJacobian(unsigned int jvar)
{
  const Real dJ = _a * (1 - _eps) * 5.0 * _Ds / _Rs;

  if (jvar == _couple_cs_var)
    return -dJ * _phi[_j][_qp] * _test[_i][_qp];
  if (jvar == _cs_avg_var)
    return dJ * _phi[_j][_qp] * _test[_i][_qp];

  return 0.0;
}
