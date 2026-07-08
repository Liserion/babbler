// Subramanian et al. (2005) two-parameter (parabolic profile) closure:
//   cs_surf = cs_avg - j_n * Rs / (5 * Ds)
// written as a residual for the nonlinear variable u = cs_surf, scaled by
// beta = Rs/(5*Ds) so its magnitude is comparable to the flux terms:
//   R = (u - cs_avg)/beta + j_n(ce, phis, phie, u)
// j_n is the per-particle-area Butler-Volmer flux (same expression as the
// micro-app surface BC in ParticleBVPostBCKernel, incl. damage factor).

#include "PolynomialCsClosureKernel.h"

registerMooseObject("babblerApp", PolynomialCsClosureKernel);

InputParameters
PolynomialCsClosureKernel::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredCoupledVar("cs_avg", "volume-averaged solid concentration");
  params.addRequiredCoupledVar("Ce", "electrolyte concentration");
  params.addRequiredCoupledVar("PhiS", "solid phase potential");
  params.addRequiredCoupledVar("PhiE", "electrolyte potential");
  params.addRequiredCoupledVar("Damage", "damage variable");
  params.addRequiredCoupledVar("SigmaH", "hydrostatic stress");

  params.addRequiredParam<Real>("Cm", "max electrolyte concentration");
  params.addRequiredParam<Real>("K2", "reaction rate");
  params.addRequiredParam<int>("MateChoice",
                               "1->TiS2, 2->Mn2O4, 3->TiS2 new, "
                               "4->LiFePO4, 5->LiFePO4 Safari, 6->V2O5");
  params.addParam<Real>("T", 298.15, "temperature");
  params.addParam<Real>("Ds", 1.0e-5, "solid diffusion coefficient");
  params.addParam<Real>("Rs", 0.5, "particle radius");
  params.addParam<Real>("Omega", 0.0, "partial molar volume");

  return params;
}

PolynomialCsClosureKernel::PolynomialCsClosureKernel(const InputParameters & parameters)
  : Kernel(parameters),
    _cs_avg(coupledValue("cs_avg")),
    _couple_c(coupledValue("Ce")),
    _couple_phi1(coupledValue("PhiS")),
    _couple_phi2(coupledValue("PhiE")),
    _cs_avg_var(coupled("cs_avg")),
    _couple_c_var(coupled("Ce")),
    _couple_phi1_var(coupled("PhiS")),
    _couple_phi2_var(coupled("PhiE")),
    _Cm(getParam<Real>("Cm")),
    _K2(getParam<Real>("K2")),
    _MateChoice(getParam<int>("MateChoice")),
    _T(getParam<Real>("T")),
    _Ds(getParam<Real>("Ds")),
    _Rs(getParam<Real>("Rs")),
    _Omega(getParam<Real>("Omega")),
    _couple_damage(coupledValue("Damage")),
    _couple_sigmaH(coupledValue("SigmaH"))
{
}

void
PolynomialCsClosureKernel::OpenCircuitV(const int & matechoice, Real x, Real & U, Real & dUdx)
{
  const Real R = 8.3144598;
  const Real F = 96485.3329;
  auto Sech = [](Real y) { return 1.0 / cosh(y); };

  U = 0.0;
  dUdx = 0.0;

  if (matechoice == 1)
  {
    U = 2.17 + (R * _T / F) * (log(fabs((1 - x) / x)) - 16.2 * x + 8.1);
    dUdx = (-16.2 * R * _T / F) * (x * x - x - 0.0617284) / (x * (x - 1));
  }
  else if (matechoice == 2)
  {
    U = 4.06279 + 0.0677504 * tanh(12.8268 - 21.8502 * x)
      - 0.105734 * (pow(1.00167 - x, -0.379571) - 1.575994)
      - 0.045 * exp(-71.69 * pow(x, 8))
      + 0.01 * exp(-200.0 * (x - 0.19));
    dUdx = -2.0 * exp(-200. * (x - 0.19))
         - 0.0401336 / pow(1.00167 - x, 1.37957)
         + 25.8084 * exp(-71.69 * pow(x, 8)) * pow(x, 7)
         - 1.48036 * Sech(12.8268 - 21.8502 * x) * Sech(12.8268 - 21.8502 * x);
  }
  else if (matechoice == 3)
  {
    U = 2.17 + R * _T * (-0.000558 * x + 8.10) / F;
    dUdx = R * _T * -0.000558 / F;
  }
  else if (matechoice == 4)
  {
    U = 3.114559
      + 4.438792 * atan(-71.7352 * x + 70.85337)
      - 4.240252 * atan(-68.5605 * x + 67.730082);
    dUdx = 4.438792 * (-71.7352) / (1 + pow(-71.7352 * x + 70.85337, 2))
         + (-4.240252) * (-68.5605) / (1 + pow(-68.5605 * x + 67.730082, 2));
  }
  else if (matechoice == 5)
  {
    U = 3.4324
      - 0.8428 * exp(-80.2493 * pow(1 - x, 1.3198))
      - (3.2474e-6) * exp(20.2645 * pow(1 - x, 3.8003))
      + (3.2482e-6) * exp(20.2646 * pow(1 - x, 3.7995));
    dUdx = -89.2635 * exp(-80.2493 * pow(1 - x, 1.3198)) * pow(1 - x, 0.3198)
         + 0.000250086 * exp(20.2645 * pow(1 - x, 3.8003)) * pow(1 - x, 2.8003)
         - 0.000250104 * exp(20.2646 * pow(1 - x, 3.7995)) * pow(1 - x, 2.7995);
  }
  else if (matechoice == 6)
  {
    U = 3.3059
      + 0.092769 * tanh(-14.362 * x + 6.6874)
      - 0.034252 * exp(100 * (x - 0.96))
      + 0.00724 * exp(80.0 * (0.01 - x));
    dUdx = 1.33235 * Sech(6.6874 - 14.362 * x) * Sech(6.6874 - 14.362 * x)
         - 3.4252 * exp(100 * (x - 0.96))
         - 0.5792 * exp(80 * (0.01 - x));
  }

  U = U * F / (R * _T);
  dUdx = dUdx * F / (R * _T);
}

void
PolynomialCsClosureKernel::BV(const Real & c, const Real & phi1, const Real & phi2,
                              const Real & cs, Real & J, Real & dJdcs, Real & dJdc,
                              Real & dJdphi1, Real & dJdphi2)
{
  // clamp cs into (0,1) for the OCV evaluation; extend flat outside
  const Real xb = std::min(std::max(cs, 1.0e-6), 1.0 - 1.0e-6);
  const bool cs_clamped = (xb != cs);

  Real U, dUdx;
  OpenCircuitV(_MateChoice, xb, U, dUdx);
  if (cs_clamped)
    dUdx = 0.0;

  // clamp c into (0, Cm) so the sqrt stays real; extend flat outside
  const Real cmin = 1.0e-3 * _Cm;
  const Real cb = std::min(std::max(c, cmin), _Cm - cmin);
  const bool c_clamped = (cb != c);

  Real eta = phi1 - phi2 - U - _Omega * _couple_sigmaH[_qp];
  Real deta_dcs = -dUdx;
  if (eta > 40.0) { eta = 40.0; deta_dcs = 0.0; }
  if (eta < -40.0) { eta = -40.0; deta_dcs = 0.0; }

  const Real sq = sqrt((_Cm - cb) * cb);
  const Real ep = exp(0.5 * eta);
  const Real em = exp(-0.5 * eta);
  const Real dmg = (_couple_damage[_qp] >= 0.9) ? 0.0 : (1.0 - _couple_damage[_qp]);

  J = dmg * _K2 * sq * (xb * ep - (1 - xb) * em);
  dJdcs = dmg * _K2 * sq * (ep + em + 0.5 * deta_dcs * (xb * ep + (1 - xb) * em));
  dJdc = c_clamped ? 0.0
                   : dmg * 0.5 * _K2 * (_Cm - 2 * cb) / sq * (xb * ep - (1 - xb) * em);
  dJdphi1 = dmg * 0.5 * _K2 * sq * (xb * ep + (1 - xb) * em);
  dJdphi2 = -dJdphi1;
}

Real
PolynomialCsClosureKernel::computeQpResidual()
{
  const Real beta = _Rs / (5.0 * _Ds);

  Real J, dJdcs, dJdc, dJdphi1, dJdphi2;
  BV(_couple_c[_qp], _couple_phi1[_qp], _couple_phi2[_qp], _u[_qp],
     J, dJdcs, dJdc, dJdphi1, dJdphi2);

  return ((_u[_qp] - _cs_avg[_qp]) / beta + J) * _test[_i][_qp];
}

Real
PolynomialCsClosureKernel::computeQpJacobian()
{
  const Real beta = _Rs / (5.0 * _Ds);

  Real J, dJdcs, dJdc, dJdphi1, dJdphi2;
  BV(_couple_c[_qp], _couple_phi1[_qp], _couple_phi2[_qp], _u[_qp],
     J, dJdcs, dJdc, dJdphi1, dJdphi2);

  return (1.0 / beta + dJdcs) * _phi[_j][_qp] * _test[_i][_qp];
}

Real
PolynomialCsClosureKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  const Real beta = _Rs / (5.0 * _Ds);

  Real J, dJdcs, dJdc, dJdphi1, dJdphi2;
  BV(_couple_c[_qp], _couple_phi1[_qp], _couple_phi2[_qp], _u[_qp],
     J, dJdcs, dJdc, dJdphi1, dJdphi2);

  if (jvar == _cs_avg_var)
    return (-1.0 / beta) * _phi[_j][_qp] * _test[_i][_qp];
  if (jvar == _couple_c_var)
    return dJdc * _phi[_j][_qp] * _test[_i][_qp];
  if (jvar == _couple_phi1_var)
    return dJdphi1 * _phi[_j][_qp] * _test[_i][_qp];
  if (jvar == _couple_phi2_var)
    return dJdphi2 * _phi[_j][_qp] * _test[_i][_qp];

  return 0.0;
}
