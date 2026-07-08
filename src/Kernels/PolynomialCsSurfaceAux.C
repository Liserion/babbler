// Computes surface concentration from volume-averaged concentration
// using the parabolic profile approximation (Subramanian et al. 2005):
//   cs_surface = cs_avg - j_n * Rs / (5 * Ds)
// where j_n is the Butler-Volmer surface flux per unit particle area.
//
// Since j_n itself depends on cs_surface (through the BV prefactor and the
// OCV inside eta), this is a scalar nonlinear equation in cs_surface:
//   f(x) = x - cs_avg + (Rs/(5*Ds)) * j_n(x) = 0
// j_n is strictly increasing in x for a physically valid OCV (dU/dx < 0),
// so f is monotone and the damped Newton below converges to the unique root.
// Execute this AuxKernel at NONLINEAR for tight coupling with the macro solve.

#include "PolynomialCsSurfaceAux.h"

registerMooseObject("babblerApp", PolynomialCsSurfaceAux);

InputParameters
PolynomialCsSurfaceAux::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("cs_avg", "volume-averaged solid concentration (Variable)");
  params.addRequiredCoupledVar("Ce", "electrolyte concentration");
  params.addRequiredCoupledVar("PhiS", "solid phase potential");
  params.addRequiredCoupledVar("PhiE", "electrolyte potential");
  params.addRequiredCoupledVar("Damage", "damage variable");
  params.addRequiredCoupledVar("SigmaH", "hydrostatic stress");

  params.addRequiredParam<Real>("Cm", "max electrolyte concentration");
  params.addRequiredParam<Real>("K2", "reaction rate");
  params.addRequiredParam<Real>("eps", "porosity");
  params.addParam<Real>("a", 3.0, "Area/Volume");
  params.addRequiredParam<int>("MateChoice", "material choice (4 = LiFePO4)");
  params.addParam<Real>("T", 298.15, "temperature");
  params.addParam<Real>("Ds", 1.0e-5, "solid diffusion coefficient");
  params.addParam<Real>("Rs", 0.5, "particle radius");
  params.addParam<Real>("Omega", 0.0, "partial molar volume");

  return params;
}

PolynomialCsSurfaceAux::PolynomialCsSurfaceAux(const InputParameters &parameters)
  : AuxKernel(parameters),
    _cs_avg(coupledValue("cs_avg")),
    _couple_c(coupledValue("Ce")),
    _couple_phi1(coupledValue("PhiS")),
    _couple_phi2(coupledValue("PhiE")),
    _Cm(getParam<Real>("Cm")),
    _K2(getParam<Real>("K2")),
    _eps(getParam<Real>("eps")),
    _a(getParam<Real>("a")),
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
PolynomialCsSurfaceAux::OpenCircuitV(const int &matechoice, Real x, Real &U, Real &dUdx)
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

Real
PolynomialCsSurfaceAux::computeValue()
{
    const Real xmin = 1.0e-6;
    const Real xmax = 1.0 - 1.0e-6;

    Real cs_avg = _cs_avg[_qp];
    if (cs_avg < xmin) cs_avg = xmin;
    if (cs_avg > xmax) cs_avg = xmax;

    const Real damage = _couple_damage[_qp];
    if (damage >= 0.9)
        return cs_avg; // no reaction flux -> uniform particle

    const Real beta = _Rs / (5.0 * _Ds);

    // Clamp the electrolyte concentration like ParticleBVPostBCKernel does,
    // so the sqrt stays real during hard transients.
    const Real cmin = 1.0e-3 * _Cm;
    Real c = _couple_c[_qp];
    if (c < cmin) c = cmin;
    if (c > _Cm - cmin) c = _Cm - cmin;
    const Real pref = (1.0 - damage) * _K2 * sqrt((_Cm - c) * c);

    const Real dphi = _couple_phi1[_qp] - _couple_phi2[_qp] - _Omega * _couple_sigmaH[_qp];

    // Initial guess: previous value of this AuxVariable, else cs_avg
    Real x = _u[_qp];
    if (!(x > xmin && x < xmax))
        x = cs_avg;

    for (unsigned int it = 0; it < 100; ++it)
    {
        Real U, dUdx;
        OpenCircuitV(_MateChoice, x, U, dUdx);

        Real eta = dphi - U;
        if (eta > 40.0) { eta = 40.0; dUdx = 0.0; }   // clamped -> eta no longer
        if (eta < -40.0) { eta = -40.0; dUdx = 0.0; } //   depends on x
        const Real ep = exp(0.5 * eta);
        const Real em = exp(-0.5 * eta);

        const Real J = pref * (x * ep - (1.0 - x) * em);
        const Real dJdx = pref * (ep + em - 0.5 * dUdx * (x * ep + (1.0 - x) * em));

        const Real f = x - cs_avg + beta * J;
        Real fp = 1.0 + beta * dJdx;
        if (fp < 1.0) fp = 1.0; // guard: f is monotone increasing for valid OCV

        Real dx = -f / fp;
        if (dx > 0.05) dx = 0.05;   // damp steps across the steep OCV ends
        if (dx < -0.05) dx = -0.05;

        x += dx;
        if (x < xmin) x = xmin;
        if (x > xmax) x = xmax;

        if (fabs(dx) < 1.0e-12)
            break;
    }

    return x;
}
