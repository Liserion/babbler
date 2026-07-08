#include "CathodeCeKernelPoly.h"

registerMooseObject("babblerApp", CathodeCeKernelPoly);

InputParameters
CathodeCeKernelPoly::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<Real>("D", "diffusivity");
  params.addRequiredParam<Real>("Cm", "Max concentration of electrolyte");
  params.addRequiredParam<Real>("eps", "porosity");
  params.addRequiredParam<Real>("K", "conductivity of electrolyte");
  params.addParam<Real>("a", 3.0, "Area/Volume");
  params.addParam<Real>("Ds", 1.0e-5, "solid diffusion coefficient");
  params.addParam<Real>("Rs", 0.5, "particle radius");
  params.addRequiredCoupledVar("PhiE", "potential for electrolyte phase");
  params.addRequiredCoupledVar("Cs", "surface concentration (closure variable)");
  params.addRequiredCoupledVar("CsAvg", "volume-averaged solid concentration");
  return params;
}

CathodeCeKernelPoly::CathodeCeKernelPoly(const InputParameters & parameters)
  : Kernel(parameters),
    _couple_cs(coupledValue("Cs")),
    _cs_avg(coupledValue("CsAvg")),
    _grad_couple_phi2(coupledGradient("PhiE")),
    _couple_phi2_var(coupled("PhiE")),
    _couple_cs_var(coupled("Cs")),
    _cs_avg_var(coupled("CsAvg")),
    _D(getParam<Real>("D")),
    _K(getParam<Real>("K")),
    _eps(getParam<Real>("eps")),
    _Cm(getParam<Real>("Cm")),
    _a(getParam<Real>("a")),
    _Ds(getParam<Real>("Ds")),
    _Rs(getParam<Real>("Rs"))
{
}

Real
CathodeCeKernelPoly::substitutedFlux() const
{
  // JEFF = a*(1-eps) * j_n with j_n = (cs_avg - cs_surf) * 5*Ds/Rs
  return _a * (1 - _eps) * (_cs_avg[_qp] - _couple_cs[_qp]) * 5.0 * _Ds / _Rs;
}

Real
CathodeCeKernelPoly::computeQpResidual()
{
  const Real t0 = 0.0107907 + _u[_qp] * 1.48837e-4;
  const Real dt0 = 1.48837e-4;

  const Real Deff = _D * _eps;
  const Real Keff = _K * _eps * sqrt(_eps);
  const Real Jeff = substitutedFlux();

  // protect the 1/c migration terms against local depletion (c -> 0)
  const Real cpos = std::max(_u[_qp], 1.0e-3 * _Cm);

  return Deff * _grad_u[_qp] * _grad_test[_i][_qp]
       - (1 - t0) * Jeff * _test[_i][_qp]
       - dt0 * Keff * (_grad_couple_phi2[_qp] - (1 - t0) * _grad_u[_qp] / cpos) * _grad_u[_qp] * _test[_i][_qp];
}

Real
CathodeCeKernelPoly::computeQpJacobian()
{
  const Real t0 = 0.0107907 + _u[_qp] * 1.48837e-4;
  const Real dt0 = 1.48837e-4;

  const Real Deff = _D * _eps;
  const Real Keff = _K * _eps * sqrt(_eps);
  const Real Jeff = substitutedFlux();

  const Real cpos = std::max(_u[_qp], 1.0e-3 * _Cm);

  return Deff * _grad_phi[_j][_qp] * _grad_test[_i][_qp]
       - dt0 * Keff * dt0 * (_grad_u[_qp] / cpos) * _grad_u[_qp] * _phi[_j][_qp] * _test[_i][_qp]
       - dt0 * Keff * (1 - t0) * (_grad_u[_qp] * _grad_u[_qp] / (cpos * cpos)) * _phi[_j][_qp] * _test[_i][_qp]
       + dt0 * Keff * (1 - t0) * (_grad_u[_qp] / cpos) * _grad_phi[_j][_qp] * _test[_i][_qp]
       - dt0 * Keff * (_grad_couple_phi2[_qp] - (1 - t0) * _grad_u[_qp] / cpos) * _grad_phi[_j][_qp] * _test[_i][_qp]
       + dt0 * Jeff * _phi[_j][_qp] * _test[_i][_qp];
}

Real
CathodeCeKernelPoly::computeQpOffDiagJacobian(unsigned int jvar)
{
  const Real t0 = 0.0107907 + _u[_qp] * 1.48837e-4;
  const Real dt0 = 1.48837e-4;

  const Real Keff = _K * _eps * sqrt(_eps);
  const Real dJ = _a * (1 - _eps) * 5.0 * _Ds / _Rs; // dJeff/dcs_avg = -dJeff/dcs

  if (jvar == _couple_phi2_var)
    return -dt0 * Keff * _grad_phi[_j][_qp] * _grad_u[_qp] * _test[_i][_qp];
  if (jvar == _couple_cs_var)
    return (1 - t0) * dJ * _phi[_j][_qp] * _test[_i][_qp];
  if (jvar == _cs_avg_var)
    return -(1 - t0) * dJ * _phi[_j][_qp] * _test[_i][_qp];

  return 0.0;
}
