/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    The Thermo class is used to define calculations of thermodynamic
    properties of chemical species.
*/

#ifndef GPC_THERMO_H
#define GPC_THERMO_H

#include "gpc_params.h"
#include "gpc_mech.h"
#include <vector>

namespace Sprog
{
namespace Thermo
{
class ThermoInterface
{
public:
    // Equation of State related functions.
    virtual real T(void) const = 0; // Returns temperature (K).
    virtual real P(void) const = 0; // Returns pressure in current units.
    virtual real Density(void) const = 0; // Returns molar density in current units.
    virtual real MassDensity(void) const = 0; // Returns mass density in current units.

    // Species concentrations/fractions.
    virtual void GetConcs(vector<real> &concs) const = 0; // Returns the molar concentrations of all species in current units.
    virtual void GetFracs(vector<real> &fracs) const = 0; // Returns the mole fractions of all species.
    virtual void GetMassFracs(vector<real> &fracs) const = 0; // Returns the mass fractions of all species.

    // Thermodynamic property calculation.
    virtual void CalcHs(vector<real> &Hs) const = 0;   // Calculates enthalpies of all species.
    virtual void CalcSs(vector<real> &Ss) const = 0;   // Calculates entropies of all species.
    virtual void CalcCps(vector<real> &Cps) const = 0; // Calculates molar heat capacity at const. P of all species.
    virtual void CalcCvs(vector<real> &Cvs) const = 0; // Calculates molar heat capacity at const. V of all species.
    virtual void CalcUs(vector<real> &Us) const = 0;   // Calculates molar internal energies of each species.
    virtual void CalcGs(vector<real> &Gs) const = 0;   // Calculates molar Gibbs free energies of each species.

    // Bulk thermodynamic property calculations.
    virtual real BulkH(void) const = 0;  // Calculates the bulk enthalpy in current units.
    virtual real BulkS(void) const = 0;  // Calculates the bulk entropy in current units.
    virtual real BulkCp(void) const = 0; // Calculates the mean molar heat capacity at const. P.
    virtual real BulkCv(void) const = 0; // Calculates the mean molar heat capacity at const. V.

    // Mechanism.
    const Sprog::Mechanism *const Mechanism(void) const; // Returns a pointer to the defining mechanism.
    void SetMechanism(const Mechanism *const mech);      // Sets the defining mechanism.

protected:
    const Sprog::Mechanism *m_mech; // Mechanism which defines species.
};
};
};
#endif