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
#include "gpc_species.h"
#include <iostream>

namespace Sprog
{
namespace Thermo
{
class IdealGas : public GasPhase
{
public:
    // Constructors.
    IdealGas(const SpeciesPtrVector &sp); // Default constructor (requires species list).
    IdealGas(const IdealGas &copy);       // Copy constructor.
    IdealGas(                      // Stream-reading constructor.
        std::istream &in,          //  - Stream from which to read.
        const SpeciesPtrVector &sp //  - Species list.
        );     

    // Destructors.
    virtual ~IdealGas(void); // Default destructor.

    // Operator overloads.
    IdealGas &operator=(const Mixture &mix);  // Assignment operator (Mixture).
    IdealGas &operator=(const GasPhase &gas); // Assignment operator (GasPhase).
    IdealGas &operator=(const IdealGas &gas); // Assignment operator.


    // EQUATION OF STATE (overrides from GasPhase).

    // Returns pressure in current units.
    real Pressure(void) const;

    // Sets the pressure of the gas phase.
    void SetPressure(real p);


    // THERMODYNAMIC PROPERTIES (overrides from GasPhase::ThermoInterface):

    // INTERNAL ENERGY.
    using GasPhase::CalcBulkU;
    using GasPhase::CalcBulkU_RT;

    // Calculates molar internal energies of each species.
    void CalcUs(
        real T,    // Temperature (K).
        fvector &U // Return vector for the internal energies.
        ) const;

    // Calculates dimensionless molar internal energies of each species.
    void CalcUs_RT(
        real T,    // Temperature (K).
        fvector &U // Return vector for the internal energies.
        ) const;

    // Calculates the bulk internal energy given the temperature and the 
    // species mole fractions.  Also returns the species energies.
    real CalcBulkU(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions array.
        unsigned int n,      // Number of values in x.
        fvector &U           // Output vector for species internal energies.
        ) const;

    // Calculates the dimensionless bulk internal energy given the
    // temperature and the species mole fractions.  C array version.
    real CalcBulkU_RT(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &U           // Output vector for species internal energies.
        ) const;


    // ENTHALPY.
    using GasPhase::CalcBulkH;
    using GasPhase::CalcBulkH_RT;

    // Calculates enthalpies of all species using the 
    // given temperature.
    void CalcHs(
        real T,    // Temperature (K).
        fvector &H // Return vector for the enthalpies
        ) const;

    // Calculates dimensionless enthalpies of all species using the 
    // given temperature.
    void CalcHs_RT(
        real T,    // Temperature (K).
        fvector &H // Return vector for the enthalpies
        ) const;

    // Calculates the bulk enthalpy given the temperature and the 
    // species mole fractions.  Also returns the species enthalpies.
    real CalcBulkH(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &H           // Output vector for species enthalpies.
        ) const;

    // Calculates the dimensionless bulk enthalpy given the temperature and the 
    // species mole fractions.  Also returns the species enthalpies.
    real CalcBulkH_RT(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &H           // Output vector for species enthalpies.
        ) const;


    // ENTROPY.
    using GasPhase::CalcBulkS;
    using GasPhase::CalcBulkS_R;

    // Calculates entropies of all species using the given
    // temperature.
    void CalcSs(
        real T,     // Temperature (K).
        fvector &S  // Return vector for the entropies.
        ) const;

    // Calculates dimensionless entropies of all species using the given
    // temperature.
    void CalcSs_R(
        real T,     // Temperature (K).
        fvector &S  // Return vector for the entropies.
        ) const;

    // Calculates the bulk entropy given the temperature and the 
    // species mole fractions.  Also returns the species entropies.
    real CalcBulkS(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &S           // Output vector for species entropies.
        ) const;

    // Calculates the dimensionless bulk entropy given the temperature and the 
    // species mole fractions.  Also returns the species entropies.
    // C array version.
    real CalcBulkS_R(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &S           // Output vector for species entropies.
        ) const;


    // GIBBS FREE ENERGY.
    using GasPhase::CalcBulkG;
    using GasPhase::CalcBulkG_RT;

    // Calculates molar Gibbs free energies of each species
    // at the given temperature.
    void CalcGs(
        real T,     // Temperature (K).
        fvector &G  // Return vector for the Gibbs free energies.
        ) const;

    // Calculates dimensionless molar Gibbs free energies of each species
    // at the given temperature.
    void CalcGs_RT(
        real T,     // Temperature (K).
        fvector &G  // Return vector for the Gibbs free energies.
        ) const;

    // Calculates the species' Gibbs free energies given the temperature and
    // the species' enthalpies and entropies.
    void CalcGs(
        real T,           // Temperature (K).
        const fvector &H, // Species enthalpies.
        const fvector &S, // Species entropies.
        fvector &G        // Return vector for the Gibbs free energies.
        ) const;

    // Calculates the bulk Gibbs free energy given the temperature and the 
    // species mole fractions.  Also returns the species Gibbs free energies.
    real CalcBulkG(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &G           // Output vector for species Gibbs free energies.
        ) const;

    // Calculates the dimensionless species' Gibbs free energies given 
    // the temperature and the species' enthalpies and entropies.
    void CalcGs_RT(
        real T,              // Temperature (K).
        const fvector &H_RT, // Species enthalpies.
        const fvector &S_R,  // Species entropies.
        fvector &G_RT        // Return vector for the Gibbs free energies.
        ) const;

    // Calculates the dimensionless bulk Gibbs free energy given the temperature and the 
    // species mole fractions.  Also returns the species Gibbs free energies.
    real CalcBulkG_RT(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &G           // Output vector for species Gibbs free energies.
        ) const;


    // CONSTANT PRESSURE HEAT CAPACITY.
    using GasPhase::CalcBulkCp;
    using GasPhase::CalcBulkCp_R;

    // Calculates molar heat capacity at const. P of all species
    // using the given temperature.
    void CalcCps(
        real T,      // Temperature (K).
        fvector &Cp  // Returns vector for the heat capacities.
        ) const;

    // Calculates dimensionless molar heat capacity at const. P of 
    // all species using the given temperature.
    void CalcCps_R(
        real T,      // Temperature (K).
        fvector &Cp  // Returns vector for the heat capacities.
        ) const;

    // Calculates the bulk constant pressure heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.
    real CalcBulkCp(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &Cp          // Output vector for species heat capacities.
        ) const;

    // Calculates the dimensionless bulk constant pressure heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.  C array version.
    real CalcBulkCp_R(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &Cp_R        // Output vector for species heat capacities.
        ) const;


    // CONSTANT VOLUME HEAT CAPACITY.
    using GasPhase::CalcBulkCv;
    using GasPhase::CalcBulkCv_R;

    // Calculates molar heat capacity at const. V of all species
    // using the given temperature.
    void CalcCvs(
        real T,      // Temperature (K).
        fvector &Cv  // Returns vector for the heat capacities.
        ) const;

    // Calculates dimensionless molar heat capacity at const. V 
    // of all species using the given temperature.
    void CalcCvs_R(
        real T,      // Temperature (K).
        fvector &Cv  // Returns vector for the heat capacities.
        ) const;

    // Calculates the bulk constant volume heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.
    real CalcBulkCv(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &Cv          // Output vector for species heat capacities.
        ) const;

    // Calculates the dimensionless bulk constant volume heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.  C array version.
    real CalcBulkCv_R(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &Cvs_R       // Output vector for species heat capacities.
        ) const;


    // MULTIPLE PROPERTIES.

    // Calculates the species' enthalpies, entropies and constant pressure
    // heat capacities using the given temperature.  This function can be 
    // more efficient than calculating the properties individually.
    void CalcCpHSs(
        real T,      // Temperature (K)
        fvector &Cp, // Output vector for constant pressure heat capacities.
        fvector &H,  // Output vector for enthalpies.
        fvector &S   // Output vector for entropies.
        ) const;

    // Calculates the species' enthalpies (H/RT), entropies (S/R) and constant pressure
    // heat capacities (Cp/R) all in dimensionless units using the given 
    // temperature.  This function can be more efficient than calculating 
    // the properties individually.
    void CalcCpHSs_RT(
        real T,      // Temperature (K)
        fvector &Cp, // Output vector for constant pressure heat capacities (Cp/R).
        fvector &H,  // Output vector for enthalpies (H/RT).
        fvector &S   // Output vector for entropies (S/R).
        ) const;


    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the mixture object.
    virtual IdealGas *const Clone() const;

    /* NOT NEEDED UNTIL CLASS HAS SOME DATA MEMBERS!
    // Writes the mixture to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the mixture data from a binary data stream.
    void Deserialize(std::istream &in);
    */

    // Identifies the mixture type for serialisation.
    Serial_MixtureType SerialType() const;

protected:
    // Like in the GasPhase class it makes no sense to define an IdealGas
    // without knowledge of the defining species, therefore the default
    // constructor is declared protected.
    IdealGas(void);

private:
    // Calculates a polynomial fit of any thermo property given the
    // temperature terms.  The polynomial coefficients are found per
    // species.
    void sumTerms(
        real T,     // Temperature (K).
        real *t,    // Temperature terms in polynomial.
        int n,      // Number of terms in polynomial.
        fvector &Xs // Output vector for all species.
        ) const;
};
};
};

#endif
