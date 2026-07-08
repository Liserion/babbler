#include "CathodePhiEKernelPoly.h"

registerMooseObject("babblerApp", CathodePhiEKernelPoly);

InputParameters
CathodePhiEKernelPoly::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<Real>("K", "conductivity");
  params.addRequiredParam<Real>("Cm", "Max concentration of electrolyte");
  params.addRequiredParam<Real>("eps", "porosity");
  params.addParam<Real>("a", 3.0, "Area/Volume");
  params.addParam<Real>("Ds", 1.0e-5, "solid diffusion coefficient");
  params.addParam<Real>("Rs", 0.5, "particle radius");
  params.addRequiredCoupledVar("Ce", "concentration for electrolyte");
  params.addRequiredCoupledVar("Cs", "surface concentration (closure variable)");
  params.addRequiredCoupledVar("CsAvg", "volume-averaged solid concentration");
  return params;
}

CathodePhiEKernelPoly::CathodePhiEKernelPoly(const InputParameters & parameters)
  : Kernel(parameters),
    _couple_c(coupledValue("Ce")),
    _couple_cs(coupledValue("Cs")),
    _cs_avg(coupledValue("CsAvg")),
    _grad_couple_c(coupledGradient("Ce")),
    _couple_c_var(coupled("Ce")),
    _couple_cs_var(coupled("Cs")),
    _cs_avg_var(coupled("CsAvg")),
    _K(getParam<Real>("K")),
    _eps(getParam<Real>("eps")),
    _Cm(getParam<Real>("Cm")),
    _a(getParam<Real>("a")),
    _Ds(getParam<Real>("Ds")),
    _Rs(getParam<Real>("Rs"))
{
}

Real
CathodePhiEKernelPoly::computeQpResidual()
{
  const Real t0 = 0.0107907 + _couple_c[_qp] * 1.48837e-4;
  const Real Keff = _K * _eps * sqrt(_eps);
  const Real Jeff = _a * (1 - _eps) * (_cs_avg[_qp] - _couple_cs[_qp]) * 5.0 * _Ds / _Rs;

  // protect the 1/c migration term against local depletion (c -> 0)
  const Real cpos = std::max(_couple_c[_qp], 1.0e-3 * _Cm);

  return Keff * (_grad_u[_qp] - (1 - t0) * _grad_couple_c[_qp] / cpos) * _grad_test[_i][_qp]
       - Jeff * _test[_i][_qp];
}

Real
CathodePhiEKernelPoly::computeQpJacobian()
{
  const Real Keff = _K * _eps * sqrt(_eps);
  return Keff * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

Real
CathodePhiEKernelPoly::computeQpOffDiagJacobian(unsigned int jvar)
{
  const Real t0 = 0.0107907 + _couple_c[_qp] * 1.48837e-4;
  const Real dt0 = 1.48837e-4;
  const Real Keff = _K * _eps * sqrt(_eps);
  const Real dJ = _a * (1 - _eps) * 5.0 * _Ds / _Rs;

  if (jvar == _couple_c_var)
  {
    const Real cpos = std::max(_couple_c[_qp], 1.0e-3 * _Cm);
    return Keff * dt0 * (_grad_couple_c[_qp] / cpos) * _phi[_j][_qp] * _grad_test[_i][_qp]
         + Keff * (1 - t0) * (_grad_couple_c[_qp] / (cpos * cpos)) * _phi[_j][_qp] * _grad_test[_i][_qp]
         - Keff * ((1 - t0) / cpos) * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
  }
  if (jvar == _couple_cs_var)
    return dJ * _phi[_j][_qp] * _test[_i][_qp];
  if (jvar == _cs_avg_var)
    return -dJ * _phi[_j][_qp] * _test[_i][_qp];

  return 0.0;
}
