// Polynomial approximation kernel for volume-averaged solid concentration.
// Replaces the CentroidMultiApp micro solve with:
//   d(cs_avg)/dt = -a * K2 * sqrt(Cm*c - c^2) * BV_terms
// where BV_terms is the Butler-Volmer reaction rate.
// The coupled variable Cs (surface concentration) is computed by
// PolynomialCsSurfaceAux from cs_avg and the BV flux.

#include "CathodeCSAvgKernel.h"

registerMooseObject("babblerApp", CathodeCSAvgKernel);

InputParameters
CathodeCSAvgKernel::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredParam<Real>("Cm", "Max concentration of electrolyte");
  params.addRequiredParam<Real>("eps", "porosity");
  params.addRequiredParam<Real>("K2", "reaction rate");
  params.addParam<Real>("a", 3.0, "Area/Volume");

  params.addRequiredParam<int>("MateChoice",
                               "1->TiS2, 2->Mn2O4, 3->TiS2 new, "
                               "4->LiFePO4, 5->LiFePO4 Safari, 6->V2O5");

  params.addParam<Real>("T", 298.15, "temperature");
  params.addParam<Real>("Ds", 1.0e-5, "solid phase diffusion coefficient");
  params.addParam<Real>("Rs", 0.5, "particle radius");

  params.addRequiredCoupledVar("Ce", "electrolyte concentration");
  params.addRequiredCoupledVar("PhiS", "solid phase potential");
  params.addRequiredCoupledVar("PhiE", "electrolyte potential");
  params.addRequiredCoupledVar("Cs", "surface concentration (AuxVariable from polynomial approx)");
  params.addRequiredCoupledVar("Damage", "damage variable");
  params.addRequiredCoupledVar("SigmaH", "hydrostatic stress");

  params.addParam<Real>("Omega", 0.0, "partial molar volume");

  return params;
}

CathodeCSAvgKernel::CathodeCSAvgKernel(const InputParameters &parameters)
  : Kernel(parameters),
    _couple_c(coupledValue("Ce")),
    _couple_phi1(coupledValue("PhiS")),
    _couple_phi2(coupledValue("PhiE")),
    _couple_cs(coupledValue("Cs")),
    _couple_c_var(coupled("Ce")),
    _couple_phi1_var(coupled("PhiS")),
    _couple_phi2_var(coupled("PhiE")),
    _couple_cs_var(coupled("Cs")),
    _eps(getParam<Real>("eps")),
    _K2(getParam<Real>("K2")),
    _Cm(getParam<Real>("Cm")),
    _a(getParam<Real>("a")),
    _MateChoice(getParam<int>("MateChoice")),
    _T(getParam<Real>("T")),
    _Ds(getParam<Real>("Ds")),
    _Rs(getParam<Real>("Rs")),
    _couple_damage(coupledValue("Damage")),
    _couple_sigmaH(coupledValue("SigmaH")),
    _Omega(getParam<Real>("Omega"))
{
}

Real CathodeCSAvgKernel::OpenCircuitV(const int &matechoice, Real x)
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

void CathodeCSAvgKernel::BV(const Real &c, const Real &phi1, const Real &phi2, const Real &cs,
                             Real &JEFF, Real &DJDC, Real &DJDPHI1, Real &DJDPHI2)
{
    // Source for cs_avg: d(cs_avg)/dt = -JEFF/(1-eps) = -a*K2*sqrt(...)*BV
    // So the coefficient is a (not a*(1-eps) as in the cathode kernels)
    Real coeff = _a;

    // Same guards as ParticleBVPostBCKernel: keep the sqrt real and eta bounded
    const Real cmin = 1.0e-3 * _Cm;
    const Real cb = std::min(std::max(c, cmin), _Cm - cmin);

    Real eta = phi1 - phi2 - OpenCircuitV(_MateChoice, cs) - _Omega * _couple_sigmaH[_qp];
    eta = std::min(std::max(eta, -40.0), 40.0);

    const Real sq = sqrt((_Cm - cb) * cb);

    JEFF = coeff * _K2 * sq * (cs * exp(0.5 * eta) - (1 - cs) * exp(-0.5 * eta));
    DJDC = coeff * 0.5 * _K2 * (_Cm - 2 * cb) / sq * (cs * exp(0.5 * eta) - (1 - cs) * exp(-0.5 * eta));
    DJDPHI1 = coeff * 0.5 * _K2 * sq * (cs * exp(0.5 * eta) + (1 - cs) * exp(-0.5 * eta));
    DJDPHI2 = -coeff * 0.5 * _K2 * sq * (cs * exp(0.5 * eta) + (1 - cs) * exp(-0.5 * eta));

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

Real CathodeCSAvgKernel::computeQpResidual()
{
    BV(_couple_c[_qp], _couple_phi1[_qp], _couple_phi2[_qp], _couple_cs[_qp],
       Jeff, dJdc, dJdphi1, dJdphi2);

    // Residual: JEFF * test (source term for cs_avg)
    // The TimeDerivative is handled by a separate TimeDerivative kernel
    return Jeff * _test[_i][_qp];
}

Real CathodeCSAvgKernel::computeQpJacobian()
{
    // cs_avg (u) doesn't appear directly in JEFF — it appears through Cs (AuxVariable)
    // Since Cs is an AuxVariable (lagged), no Jacobian w.r.t. u
    return 0.0;
}

Real CathodeCSAvgKernel::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
    // Return 0 — let PJFNK handle the coupling through finite differences.
    // Providing explicit off-diagonal entries here makes the preconditioning
    // matrix ill-conditioned (huge off-diag vs small TimeDerivative diagonal).
    return 0.0;
}
