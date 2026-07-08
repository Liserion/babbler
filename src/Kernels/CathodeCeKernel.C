#include "CathodeCeKernel.h"

registerMooseObject("babblerApp", CathodeCeKernel);

InputParameters
CathodeCeKernel::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<Real>("D", "diffusivity");
  params.addRequiredParam<Real>("Cm", "Max concentration of electrolyte");
  params.addRequiredParam<Real>("eps", "porosity");
  params.addRequiredParam<Real>("K", "conductivity of electrolyte");
  params.addRequiredParam<Real>("K2", "reaction rate");
  params.addParam<Real>("a", 3.0, "Area/Volume");
  params.addRequiredParam<int>("MateChoice",
                               "1->TiS2, 2->Mn2O4, 3->TiS2 new, "
                               "4->LiFePO4, 5->LiFePO4 from Safari, 6->V2O5");
  params.addParam<int>("IsDebug", 0, "0:->Not print, 1:->Print Js");
  params.addParam<Real>("T", 298.15, "temperature");
  params.addRequiredCoupledVar("PhiS", "potential for solid phase");
  params.addRequiredCoupledVar("PhiE", "potential for electrolyte phase");
  params.addRequiredCoupledVar("Cs", "surface concentration of solid particle");
  params.addRequiredCoupledVar("Damage", "damage from rve");
  params.addRequiredCoupledVar("SigmaH", "Hydrostatic stress from RVE homogenization");
  params.addParam<Real>("Omega", 0.0, "Partial molar volume of RVE material");
  return params;
}

CathodeCeKernel::CathodeCeKernel(const InputParameters & parameters)
  : Kernel(parameters),
    _couple_phi1(coupledValue("PhiS")),
    _couple_phi2(coupledValue("PhiE")),
    _couple_cs(coupledValue("Cs")),
    _grad_couple_phi2(coupledGradient("PhiE")),
    _couple_phi1_var(coupled("PhiS")),
    _couple_phi2_var(coupled("PhiE")),
    _D(getParam<Real>("D")),
    _K(getParam<Real>("K")),
    _eps(getParam<Real>("eps")),
    _K2(getParam<Real>("K2")),
    _Cm(getParam<Real>("Cm")),
    _a(getParam<Real>("a")),
    _T(getParam<Real>("T")),
    _MateChoice(getParam<int>("MateChoice")),
    _IsDebug(getParam<int>("IsDebug")),
    _couple_damage(coupledValue("Damage")),
    _couple_sigmaH(coupledValue("SigmaH")),
    _Omega(getParam<Real>("Omega"))
{
}

Real
CathodeCeKernel::OpenCircuitV(const int & matechoice, Real x)
{
  Real V = 0.0;
  const Real R = 8.3144598;
  const Real F = 96485.3329;

  if (matechoice == 1)
    V = 2.17 + (R * _T / F) * (log(fabs((1 - x) / x)) - 16.2 * x + 8.1);
  else if (matechoice == 2)
    V = 4.06279 + 0.0677504 * tanh(12.8268 - 21.8502 * x)
      - 0.105734 * (pow(1.00167 - x, -0.379571) - 1.575994)
      - 0.045 * exp(-71.69 * pow(x, 8))
      + 0.01 * exp(-200.0 * (x - 0.19));
  else if (matechoice == 3)
    V = 2.17 + R * _T * (-0.000558 * x + 8.10) / F;
  else if (matechoice == 4)
    V = 3.114559
      + 4.438792 * atan(-71.7352 * x + 70.85337)
      - 4.240252 * atan(-68.5605 * x + 67.730082);
  else if (matechoice == 5)
    V = 3.4324
      - 0.8428 * exp(-80.2493 * pow(1 - x, 1.3198))
      - (3.2474e-6) * exp(20.2645 * pow(1 - x, 3.8003))
      + (3.2482e-6) * exp(20.2646 * pow(1 - x, 3.7995));
  else if (matechoice == 6)
    V = 3.3059
      + 0.092769 * tanh(-14.362 * x + 6.6874)
      - 0.034252 * exp(100 * (x - 0.96))
      + 0.00724 * exp(80.0 * (0.01 - x));

  return V * F / (R * _T);
}

void
CathodeCeKernel::BV(const Real & c, const Real & phi1, const Real & phi2, const Real & cs,
                     Real & JEFF, Real & DJDC, Real & DJDPHI1, Real & DJDPHI2)
{
  Real a = _a * (1 - _eps);

  // Clamp c into its physical range (0, Cm) so the sqrt prefactor stays real
  // and differentiable when the electrolyte locally depletes or saturates;
  // outside the range J is extended flat (dJ/dc = 0). Bound eta so exp()
  // cannot overflow during Newton iterations.
  const Real cmin = 1.0e-3 * _Cm;
  const Real cb = std::min(std::max(c, cmin), _Cm - cmin);
  const bool c_clamped = (cb != c);

  Real eta = phi1 - phi2 - OpenCircuitV(_MateChoice, cs) - _Omega * _couple_sigmaH[_qp];
  eta = std::min(std::max(eta, -40.0), 40.0);

  JEFF = a * _K2 * sqrt(_Cm * cb - cb * cb) * (cs * exp(0.5 * eta) - (1 - cs) * exp(-0.5 * eta));
  DJDC = c_clamped ? 0.0 : a * 0.5 * _K2 * (_Cm - 2 * cb) / sqrt(_Cm * cb - cb * cb) * (cs * exp(0.5 * eta) - (1 - cs) * exp(-0.5 * eta));
  DJDPHI1 = a * 0.5 * _K2 * sqrt(_Cm * cb - cb * cb) * (cs * exp(0.5 * eta) + (1 - cs) * exp(-0.5 * eta));
  DJDPHI2 = -a * 0.5 * _K2 * sqrt(_Cm * cb - cb * cb) * (cs * exp(0.5 * eta) + (1 - cs) * exp(-0.5 * eta));

  // add damage influence
  JEFF = (1 - _couple_damage[_qp]) * JEFF;
  DJDC = (1 - _couple_damage[_qp]) * DJDC;
  DJDPHI1 = (1 - _couple_damage[_qp]) * DJDPHI1;
  DJDPHI2 = (1 - _couple_damage[_qp]) * DJDPHI2;
  if (_couple_damage[_qp] >= 0.9)
  {
    JEFF = 0.0;
    DJDC = 0.0;
    DJDPHI1 = 0.0;
    DJDPHI2 = 0.0;
  }
}

Real
CathodeCeKernel::computeQpResidual()
{
  Real t0 = 0.0107907 + _u[_qp] * 1.48837e-4;
  Real dt0 = 1.48837e-4;

  Deff = _D * _eps;
  Keff = _K * _eps * sqrt(_eps);

  BV(_u[_qp], _couple_phi1[_qp], _couple_phi2[_qp], _couple_cs[_qp], Jeff, dJdc, dJdphi1, dJdphi2);

  if (_IsDebug)
  {
    std::cout << "Js=" << Jeff
              << ",dJdc=" << dJdc
              << ",dJdphi1=" << dJdphi1
              << ",dJdphi2=" << dJdphi2 << std::endl;
    std::cout << "c=" << _u[_qp]
              << ",Phi1=" << _couple_phi1[_qp]
              << ",Phi2=" << _couple_phi2[_qp]
              << ",Cs=" << _couple_cs[_qp] << std::endl << std::endl;
  }

  // protect the 1/c migration terms against local depletion (c -> 0)
  const Real cpos = std::max(_u[_qp], 1.0e-3 * _Cm);

  return Deff * _grad_u[_qp] * _grad_test[_i][_qp]
       - (1 - t0) * Jeff * _test[_i][_qp]
       - dt0 * Keff * (_grad_couple_phi2[_qp] - (1 - t0) * _grad_u[_qp] / cpos) * _grad_u[_qp] * _test[_i][_qp];
}

Real
CathodeCeKernel::computeQpJacobian()
{
  Real t0 = 0.0107907 + _u[_qp] * 1.48837e-4;
  Real dt0 = 1.48837e-4;

  BV(_u[_qp], _couple_phi1[_qp], _couple_phi2[_qp], _couple_cs[_qp], Jeff, dJdc, dJdphi1, dJdphi2);

  Deff = _D * _eps;
  Keff = _K * _eps * sqrt(_eps);

  const Real cpos = std::max(_u[_qp], 1.0e-3 * _Cm);

  return Deff * _grad_phi[_j][_qp] * _grad_test[_i][_qp]
       - dt0 * Keff * dt0 * (_grad_u[_qp] / cpos) * _grad_u[_qp] * _phi[_j][_qp] * _test[_i][_qp]
       - dt0 * Keff * (1 - t0) * (_grad_u[_qp] * _grad_u[_qp] / (cpos * cpos)) * _phi[_j][_qp] * _test[_i][_qp]
       + dt0 * Keff * (1 - t0) * (_grad_u[_qp] / cpos) * _grad_phi[_j][_qp] * _test[_i][_qp]
       - dt0 * Keff * (_grad_couple_phi2[_qp] - (1 - t0) * _grad_u[_qp] / cpos) * _grad_phi[_j][_qp] * _test[_i][_qp]
       + dt0 * Jeff * _phi[_j][_qp] * _test[_i][_qp]
       - (1 - t0) * dJdc * _phi[_j][_qp] * _test[_i][_qp];
}

Real
CathodeCeKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real t0 = 0.0107907 + _u[_qp] * 1.48837e-4;
  Real dt0 = 1.48837e-4;

  Deff = _D * _eps;
  Keff = _K * _eps * sqrt(_eps);

  if (jvar == _couple_phi1_var)
  {
    BV(_u[_qp], _couple_phi1[_qp], _couple_phi2[_qp], _couple_cs[_qp], Jeff, dJdc, dJdphi1, dJdphi2);
    return -(1 - t0) * dJdphi1 * _phi[_j][_qp] * _test[_i][_qp];
  }
  else if (jvar == _couple_phi2_var)
  {
    BV(_u[_qp], _couple_phi1[_qp], _couple_phi2[_qp], _couple_cs[_qp], Jeff, dJdc, dJdphi1, dJdphi2);
    return -dt0 * Keff * _grad_phi[_j][_qp] * _grad_u[_qp] * _test[_i][_qp]
         - (1 - t0) * dJdphi2 * _phi[_j][_qp] * _test[_i][_qp];
  }

  return 0.0;
}
