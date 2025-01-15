#pragma once

#include "IntegratedBC.h"  // Ensure all dependencies are included here

class ConstFluxForCeBC : public IntegratedBC
{
public:
    ConstFluxForCeBC(const InputParameters &parameters);

    // New validParams format
    static InputParameters validParams();

protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;

    const Real &_I;
    const Real &_ChargeTime;
};
