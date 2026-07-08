#pragma once

#include "AuxKernel.h"

class PolynomialCsSurfaceAux : public AuxKernel
{
public:
    static InputParameters validParams();
    explicit PolynomialCsSurfaceAux(const InputParameters &parameters);

protected:
    Real computeValue() override;

    const VariableValue &_cs_avg;
    const VariableValue &_couple_c, &_couple_phi1, &_couple_phi2;

    const Real &_Cm, &_K2, &_eps, &_a;
    const int &_MateChoice;
    const Real &_T, &_Ds, &_Rs;
    const Real &_Omega;

    const VariableValue &_couple_damage, &_couple_sigmaH;

    void OpenCircuitV(const int &matechoice, Real x, Real &U, Real &dUdx);
};
