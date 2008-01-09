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
#include <vector>

namespace Sprog
{
//class Mechanism; // Forward declaration of mechanism.

namespace Thermo
{
class ThermoInterface
{
public:
    // Constructors
    ThermoInterface(void); // Default constructor.

    // Destructors.
    virtual ~ThermoInterface(void); // Default destructor.

    // Equation of State related functions.
    virtual real Pressure(void) const = 0; // Returns pressure in current units.

    // Thermodynamic property calculation.
    virtual void CalcHs(std::vector<real> &Hs) const = 0;   // Calculates enthalpies of all species.
    virtual void CalcSs(std::vector<real> &Ss) const = 0;   // Calculates entropies of all species.
    virtual void CalcCps(std::vector<real> &Cps) const = 0; // Calculates molar heat capacity at const. P of all species.
    virtual void CalcCvs(std::vector<real> &Cvs) const = 0; // Calculates molar heat capacity at const. V of all species.
    virtual void CalcUs(std::vector<real> &Us) const = 0;   // Calculates molar internal energies of each species.
    virtual void CalcGs(std::vector<real> &Gs) const = 0;   // Calculates molar Gibbs free energies of each species.

    // Bulk thermodynamic property calculations.
    virtual real BulkH(void) const = 0;  // Calculates the bulk enthalpy in current units.
    virtual real BulkS(void) const = 0;  // Calculates the bulk entropy in current units.
    virtual real BulkCp(void) const = 0; // Calculates the mean molar heat capacity at const. P.
    virtual real BulkCv(void) const = 0; // Calculates the mean molar heat capacity at const. V.

    /*
    // Mechanism.
    const Sprog::Mechanism *const Mechanism(void) const;   // Returns a pointer to the defining mechanism.
    void SetMechanism(const Sprog::Mechanism *const mech); // Sets the defining mechanism.
*/

protected:
//    const Sprog::Mechanism *m_mech; // Mechanism which defines species.
};
};
};

#endif