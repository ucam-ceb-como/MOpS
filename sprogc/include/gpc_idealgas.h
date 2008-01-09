/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    The IdealGas class defines a homogeneous ideal gas mixture.
*/

#ifndef GPC_IDEALGAS_H
#define GPC_IDEALGAS_H

#include "gpc_params.h"
#include "gpc_gasphase.h"
#include <vector>

namespace Sprog
{
namespace Thermo
{
class IdealGas : public GasPhase
{
public:
    // Constructors.
    IdealGas(void); // Default constructor.

    // Destructors.
    virtual ~IdealGas(void); // Default destructor.


    // Required function overrides for GasPhase:

    // Equation of State related functions.
    real Pressure(void) const; // Returns pressure in current units.

    // Thermodynamic property calculation.
    void CalcHs(std::vector<real> &Hs) const;   // Calculates enthalpies of all species.
    void CalcSs(std::vector<real> &Ss) const;   // Calculates entropies of all species.
    void CalcCps(std::vector<real> &Cps) const; // Calculates molar heat capacity at const. P of all species.
    void CalcCvs(std::vector<real> &Cvs) const; // Calculates molar heat capacity at const. V of all species.
    void CalcUs(std::vector<real> &Us) const;   // Calculates molar internal energies of each species.
    void CalcGs(std::vector<real> &Gs) const;   // Calculates molar Gibbs free energies of each species.

    void CalcHs_RT(std::vector<real> &Hs) const;  // Calculates enthalpies of all species / RT.
    void CalcSs_R(std::vector<real> &Ss) const;   // Calculates entropies of all species / R.
    void CalcCps_R(std::vector<real> &Cps) const; // Calculates molar heat capacity at const. P of all species / RT.
    void CalcCvs_R(std::vector<real> &Cvs) const; // Calculates molar heat capacity at const. V of all species / RT.
    void CalcUs_RT(std::vector<real> &Us) const;  // Calculates molar internal energies of each species / RT.
    void CalcGs_RT(std::vector<real> &Gs) const;  // Calculates molar Gibbs free energies of each species / RT.

    // Bulk thermodynamic property calculations.
    real BulkH(void) const;  // Calculates the bulk enthalpy in current units.
    real BulkS(void) const;  // Calculates the bulk entropy in current units.
    real BulkCp(void) const; // Calculates the mean molar heat capacity at const. P.
    real BulkCv(void) const; // Calculates the mean molar heat capacity at const. V.

private:
    // Calculates a polynomial fit of any thermo property given the
    // temperature terms.  The polynomial coefficients are found per
    // species.
    void sumTerms(real *t, int n, std::vector<real> &Xs) const;
};
};
};

#endif