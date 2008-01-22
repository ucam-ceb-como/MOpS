/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    The Thermo class is used to define calculations of thermodynamic
    properties of chemical species and mixtures.  This is a virtual
    class which provides a common interface for thermodynamic property
    calculator classes.

    It is expected that all derived classes follow the sprogc units
    methodology, such that temperatures are in Kelvin (K) and species
    concentrations are supplied as mole fractions and, where applicable,
    the mixture total concentration.
*/

#ifndef GPC_THERMO_H
#define GPC_THERMO_H

#include "gpc_params.h"
#include <vector>

namespace Sprog
{
namespace Thermo
{
class ThermoInterface
{
public:
    // Constructors
    ThermoInterface(void); // Default constructor.

    // Destructors.
    virtual ~ThermoInterface(void); // Default destructor.


    // INTERNAL ENERGY.

    // Calculates molar internal energies of each species.
    virtual void CalcUs(
        real T,    // Temperature (K).
        fvector &U // Return vector for the internal energies.
        ) const = 0;

    // Calculates dimensionless molar internal energies of each species.
    virtual void CalcUs_RT(
        real T,    // Temperature (K).
        fvector &U // Return vector for the internal energies.
        ) const = 0;

    // Calculates the bulk internal energy given the
    // temperature and the species mole fractions.
    real CalcBulkU(
        real T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;
    real CalcBulkU(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk internal energy given the temperature and the 
    // species mole fractions.  Also returns the species energies.
    real CalcBulkU(
        real T,           // Temperature (K).
        const fvector &x, // Species mole fractions array.
        fvector &U        // Output vector for species internal energies.
        ) const;
    virtual real CalcBulkU(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &U           // Output vector for species internal energies.
        ) const = 0;


    // ENTHALPY.

    // Calculates enthalpies of all species using the 
    // given temperature.
    virtual void CalcHs(
        real T,    // Temperature (K).
        fvector &H // Return vector for the enthalpies
        ) const = 0;

    // Calculates dimensionless enthalpies of all species using the 
    // given temperature.
    virtual void CalcHs_RT(
        real T,    // Temperature (K).
        fvector &H // Return vector for the enthalpies
        ) const = 0;

    // Calculates the bulk enthalpy given the
    // temperature and the species mole fractions.
    real CalcBulkH(
        real T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;
    real CalcBulkH(
        real T,             // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk enthalpy given the temperature and the 
    // species mole fractions.  Also returns the species enthalpies.
    real CalcBulkH(
        real T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &H        // Output vector for species enthalpies.
        ) const;
    virtual real CalcBulkH(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &H           // Output vector for species enthalpies.
        ) const = 0;


    // ENTROPY.

    // Calculates entropies of all species using the given
    // temperature.
    virtual void CalcSs(
        real T,     // Temperature (K).
        fvector &S  // Return vector for the entropies.
        ) const = 0;

    // Calculates dimensionless entropies of all species using the given
    // temperature.
    virtual void CalcSs_R(
        real T,     // Temperature (K).
        fvector &S  // Return vector for the entropies.
        ) const = 0;

    // Calculates the bulk entropy given the
    // temperature and the species mole fractions.
    real CalcBulkS(
        real T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;
    real CalcBulkS(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk entropy given the temperature and the 
    // species mole fractions.  Also returns the species entropies.
    real CalcBulkS(
        real T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &S        // Output vector for species entropies.
        ) const;
    virtual real CalcBulkS(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &S           // Output vector for species entropies.
        ) const = 0;


    // GIBBS FREE ENERGY.

    // Calculates molar Gibbs free energies of each species
    // at the given temperature.
    virtual void CalcGs(
        real T,     // Temperature (K).
        fvector &G  // Return vector for the Gibbs free energies.
        ) const = 0;

    // Calculates dimensionless molar Gibbs free energies of each species
    // at the given temperature.
    virtual void CalcGs_RT(
        real T,     // Temperature (K).
        fvector &G  // Return vector for the Gibbs free energies.
        ) const = 0;

    // Calculates the species' Gibbs free energies given the temperature and
    // the species' enthalpies and entropies.
    virtual void CalcGs(
        real T,           // Temperature (K).
        const fvector &H, // Species enthalpies.
        const fvector &S, // Species entropies.
        fvector &G        // Return vector for the Gibbs free energies.
        ) const = 0;

    // Calculates the bulk Gibbs free energy given the
    // temperature and the species mole fractions.
    real CalcBulkG(
        real T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;
    real CalcBulkG(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk Gibbs free energy given the temperature and the 
    // species mole fractions.  Also returns the species Gibbs free energies.
    real CalcBulkG(
        real T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &G        // Output vector for species Gibbs free energies.
        ) const;
    virtual real CalcBulkG(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &G           // Output vector for species Gibbs free energies.
        ) const = 0;


    // CONSTANT PRESSURE HEAT CAPACITY.

    // Calculates molar heat capacity at const. P of all species
    // using the given temperature.
    virtual void CalcCps(
        real T,      // Temperature (K).
        fvector &Cp  // Returns vector for the heat capacities.
        ) const = 0;

    // Calculates dimensionless molar heat capacity at const. P of 
    // all species using the given temperature.
    virtual void CalcCps_R(
        real T,      // Temperature (K).
        fvector &Cp  // Returns vector for the heat capacities.
        ) const = 0;

    // Calculates the bulk constant pressure heat capacity given the
    // temperature and the species mole fractions.
    real CalcBulkCp(
        real T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;
    real CalcBulkCp(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk constant pressure heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.
    real CalcBulkCp(
        real T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &Cp       // Output vector for species heat capacities.
        ) const;
    virtual real CalcBulkCp(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &Cp          // Output vector for species heat capacities.
        ) const = 0;


    // CONSTANT VOLUME HEAT CAPACITY.

    // Calculates molar heat capacity at const. V of all species
    // using the given temperature.
    virtual void CalcCvs(
        real T,      // Temperature (K).
        fvector &Cv  // Returns vector for the heat capacities.
        ) const = 0;

    // Calculates dimensionless molar heat capacity at const. V 
    // of all species using the given temperature.
    virtual void CalcCvs_R(
        real T,      // Temperature (K).
        fvector &Cv  // Returns vector for the heat capacities.
        ) const = 0;

    // Calculates the bulk constant volume heat capacity given the
    // temperature and the species mole fractions.
    real CalcBulkCv(
        real T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;
    real CalcBulkCv(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk constant volume heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.
    real CalcBulkCv(
        real T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &Cv       // Output vector for species heat capacities.
        ) const;
    virtual real CalcBulkCv(
        real T,              // Temperature (K).
        const real *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &Cv          // Output vector for species heat capacities.
        ) const = 0;


    // MULTIPLE PROPERTIES.

    // Calculates the species' enthalpies, entropies and constant pressure
    // heat capacities using the given temperature.  This function can be 
    // more efficient than calculating the properties individually.
    virtual void CalcCpHSs(
        real T,      // Temperature (K)
        fvector &Cp, // Output vector for constant pressure heat capacities.
        fvector &H,  // Output vector for enthalpies.
        fvector &S   // Output vector for entropies.
        ) const = 0;

    // Calculates the species' enthalpies (H/RT), entropies (S/R) and constant pressure
    // heat capacities (Cp/R) all in dimensionless units using the given 
    // temperature.  This function can be more efficient than calculating 
    // the properties individually.
    virtual void CalcCpHSs_RT(
        real T,      // Temperature (K)
        fvector &Cp, // Output vector for constant pressure heat capacities (Cp/R).
        fvector &H,  // Output vector for enthalpies (H/RT).
        fvector &S   // Output vector for entropies (S/R).
        ) const = 0;
};
};
};

#endif
