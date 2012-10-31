/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (chemical thermodynamics & kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ThermoInterface class is used to define calculations of thermodynamic
    properties of chemical species and mixtures.  This is a virtual
    class which provides a common interface for thermodynamic property
    calculator classes.

    It is expected that all derived classes follow the sprogc units
    methodology, such that temperatures are in Kelvin (K) and species
    concentrations are supplied as mole fractions and, where applicable,
    the mixture total concentration.

    The thermodynamic properties currently accounted for in this interface
    are:

    1.  Internal Energy (U).
    2.  Enthalpy (H).
    3.  Entropy (S).
    4.  Gibbs Free Energy (G).
    5.  Const. P Heat Capacity (Cp).
    6.  Const. V Heat Capacity (Cv).

    Functions are provided which calculate the individual species properties
    (partial properties) and bulk mixture properties (if the composition is
    provided).  When the composition is provided two functional forms are
    given:  STL vector form and C array form.  Only the C array form is
    actually required in the derived class, as the STL form simply calls
    the C array form.

    Functions are provided which calculate true thermodynamic property
    values, and function are provided which calculate dimensionless forms
    (either X/R or X/RT).

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#ifndef GPC_THERMO_H
#define GPC_THERMO_H

#include "gpc_params.h"

namespace Sprog
{
namespace Thermo
{
class ThermoInterface
{
public:
    // Destructors.
    virtual ~ThermoInterface(void); // Default destructor.



    // THE FOLLOWING FUNCTIONS MUST BE DEFINED BY ANY
    // DERIVED CLASS:


    // INTERNAL ENERGY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates molar internal energies of each species.
    virtual void CalcUs(
        double T,    // Temperature (K).
        fvector &U // Return vector for the internal energies.
        ) const = 0;

    // Calculates the bulk internal energy given the temperature and the 
    // species mole fractions.  Also returns the species energies.  
    // C array version.
    virtual double CalcBulkU(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &U           // Output vector for species internal energies.
        ) const = 0;


    // DIMENSIONLESS INTERNAL ENERGY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates dimensionless molar internal energies of each species.
    virtual void CalcUs_RT(
        double T,    // Temperature (K).
        fvector &U // Return vector for the internal energies.
        ) const = 0;

    // Calculates the dimensionless bulk internal energy given the
    // temperature and the species mole fractions.  C array version.
    virtual double CalcBulkU_RT(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &U           // Output vector for species internal energies.
        ) const = 0;


    // ENTHALPY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates enthalpies of all species using the 
    // given temperature.
    virtual void CalcHs(
        double T,    // Temperature (K).
        fvector &H // Return vector for the enthalpies
        ) const = 0;

    // Calculates the bulk enthalpy given the temperature and the 
    // species mole fractions.  Also returns the species enthalpies.
    // C array version.
    virtual double CalcBulkH(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &H           // Output vector for species enthalpies.
        ) const = 0;


    // DIMENSIONLESS ENTHALPY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates dimensionless enthalpies of all species using the 
    // given temperature.
    virtual void CalcHs_RT(
        double T,    // Temperature (K).
        fvector &H // Return vector for the enthalpies
        ) const = 0;

    // Calculates the dimensionless bulk enthalpy given the temperature and the 
    // species mole fractions.  Also returns the species enthalpies.
    // C array version.
    virtual double CalcBulkH_RT(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &H           // Output vector for species enthalpies.
        ) const = 0;


    // ENTROPY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates entropies of all species using the given
    // temperature.
    virtual void CalcSs(
        double T,     // Temperature (K).
        fvector &S  // Return vector for the entropies.
        ) const = 0;

    // Calculates the bulk entropy given the temperature and the 
    // species mole fractions.  Also returns the species entropies.
    // C array version.
    virtual double CalcBulkS(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &S           // Output vector for species entropies.
        ) const = 0;


    // DIMENSIONLESS ENTROPY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates dimensionless entropies of all species using the given
    // temperature.
    virtual void CalcSs_R(
        double T,     // Temperature (K).
        fvector &S  // Return vector for the entropies.
        ) const = 0;

    // Calculates the dimensionless bulk entropy given the temperature and the 
    // species mole fractions.  Also returns the species entropies.
    // C array version.
    virtual double CalcBulkS_R(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &S           // Output vector for species entropies.
        ) const = 0;


    // GIBBS FREE ENERGY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates molar Gibbs free energies of each species
    // at the given temperature.
    virtual void CalcGs(
        double T,     // Temperature (K).
        fvector &G  // Return vector for the Gibbs free energies.
        ) const = 0;

    // Calculates the species' Gibbs free energies given the temperature and
    // the species' enthalpies and entropies.
    virtual void CalcGs(
        double T,           // Temperature (K).
        const fvector &H, // Species enthalpies.
        const fvector &S, // Species entropies.
        fvector &G        // Return vector for the Gibbs free energies.
        ) const = 0;

    // Calculates the bulk Gibbs free energy given the temperature and the 
    // species mole fractions.  Also returns the species Gibbs free energies.
    // C array version.
    virtual double CalcBulkG(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &G           // Output vector for species Gibbs free energies.
        ) const = 0;


    // DIMENSIONLESS GIBBS FREE ENERGY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates dimensionless molar Gibbs free energies of each species
    // at the given temperature.
    virtual void CalcGs_RT(
        double T,     // Temperature (K).
        fvector &G  // Return vector for the Gibbs free energies.
        ) const = 0;

    // Calculates the dimensionless species' Gibbs free energies given 
    // the temperature and the species' enthalpies and entropies.
    virtual void CalcGs_RT(
        double T,              // Temperature (K).
        const fvector &H_RT, // Species enthalpies.
        const fvector &S_R,  // Species entropies.
        fvector &G           // Return vector for the Gibbs free energies.
        ) const = 0;

    // Calculates the dimensionless bulk Gibbs free energy given the temperature and the 
    // species mole fractions.  Also returns the species Gibbs free energies.
    // C array version.
    virtual double CalcBulkG_RT(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &G           // Output vector for species Gibbs free energies.
        ) const = 0;


    // CONSTANT PRESSURE HEAT CAPACITY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates molar heat capacity at const. P of all species
    // using the given temperature.
    virtual void CalcCps(
        double T,      // Temperature (K).
        fvector &Cp  // Returns vector for the heat capacities.
        ) const = 0;

    // Calculates the bulk constant pressure heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.  C array version.
    virtual double CalcBulkCp(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &Cp          // Output vector for species heat capacities.
        ) const = 0;


    // DIMENSIONLESS CONSTANT PRESSURE HEAT CAPACITY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates dimensionless molar heat capacity at const. P of 
    // all species using the given temperature.
    virtual void CalcCps_R(
        double T,       // Temperature (K).
        fvector &Cp_R // Returns vector for the heat capacities.
        ) const = 0;

    // Calculates the dimensionless bulk constant pressure heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.  C array version.
    virtual double CalcBulkCp_R(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &Cp_R        // Output vector for species heat capacities.
        ) const = 0;


    // CONSTANT VOLUME HEAT CAPACITY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates molar heat capacity at const. V of all species
    // using the given temperature.
    virtual void CalcCvs(
        double T,      // Temperature (K).
        fvector &Cv  // Returns vector for the heat capacities.
        ) const = 0;

    // Calculates the bulk constant volume heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.  C array version.
    virtual double CalcBulkCv(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &Cv          // Output vector for species heat capacities.
        ) const = 0;


    // DIMENSIONLESS CONSTANT VOLUME HEAT CAPACITY (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates dimensionless molar heat capacity at const. V 
    // of all species using the given temperature.
    virtual void CalcCvs_R(
        double T,      // Temperature (K).
        fvector &Cv  // Returns vector for the heat capacities.
        ) const = 0;

    // Calculates the dimensionless bulk constant volume heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.  C array version.
    virtual double CalcBulkCv_R(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n,      // Number of values in x.
        fvector &Cvs_R       // Output vector for species heat capacities.
        ) const = 0;


    // MULTIPLE PROPERTIES (REQUIRED DERIVED CLASS INTERFACE).

    // Calculates the species' enthalpies, entropies and constant pressure
    // heat capacities using the given temperature.  This function can be 
    // more efficient than calculating the properties individually.
    virtual void CalcCpHSs(
        double T,      // Temperature (K)
        fvector &Cp, // Output vector for constant pressure heat capacities.
        fvector &H,  // Output vector for enthalpies.
        fvector &S   // Output vector for entropies.
        ) const = 0;

    // Calculates the species' enthalpies (H/RT), entropies (S/R) and constant pressure
    // heat capacities (Cp/R) all in dimensionless units using the given 
    // temperature.  This function can be more efficient than calculating 
    // the properties individually.
    virtual void CalcCpHSs_RT(
        double T,      // Temperature (K)
        fvector &Cp, // Output vector for constant pressure heat capacities (Cp/R).
        fvector &H,  // Output vector for enthalpies (H/RT).
        fvector &S   // Output vector for entropies (S/R).
        ) const = 0;




    // THE FOLLOWING FUNCTION SIMPLY CALL THE ABOVE FUNCTIONS.  THEY
    // DO NOT NEED TO BE RE-DEFINED BY A DERIVED CLASS.


    // INTERNAL ENERGY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the bulk internal energy given the
    // temperature and the species mole fractions.  STL version.
    double CalcBulkU(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the bulk internal energy given the temperature
    // and the species mole fractions.  C array version.
    double CalcBulkU(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk internal energy given the temperature and the 
    // species mole fractions.  Also returns the species energies.  STL version.
    double CalcBulkU(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions array.
        fvector &U        // Output vector for species internal energies.
        ) const;


    // DIMENSIONLESS INTERNAL ENERGY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the dimensionless bulk internal energy given the
    // temperature and the species mole fractions.  STL version.
    double CalcBulkU_RT(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the dimensionless bulk internal energy given the
    // temperature and the species mole fractions.  C array version.
    double CalcBulkU_RT(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the dimensionless bulk internal energy given the temperature and the 
    // species mole fractions.  Also returns the species energies.  STL version.
    double CalcBulkU_RT(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions array.
        fvector &U        // Output vector for species internal energies.
        ) const;


    // ENTHALPY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the bulk enthalpy given the temperature
    // and the species mole fractions.  STL version.
    double CalcBulkH(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the bulk enthalpy given the temperature
    // and the species mole fractions.  C array version.
    double CalcBulkH(
        double T,             // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk enthalpy given the temperature and the 
    // species mole fractions.  Also returns the species enthalpies.
    // STL version.
    double CalcBulkH(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &H        // Output vector for species enthalpies.
        ) const;


    // DIMENSIONLESS ENTHALPY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the bulk enthalpy given the temperature
    // and the species mole fractions.  STL version.
    double CalcBulkH_RT(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the bulk enthalpy given the temperature
    // and the species mole fractions.  C array version.
    double CalcBulkH_RT(
        double T,             // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk enthalpy given the temperature and the 
    // species mole fractions.  Also returns the species enthalpies.
    // STL version.
    double CalcBulkH_RT(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &H        // Output vector for species enthalpies.
        ) const;


    // ENTROPY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the bulk entropy given the temperature
    // and the species mole fractions.  STL version.
    double CalcBulkS(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the bulk entropy given the temperature
    // and the species mole fractions.  C array version.
    double CalcBulkS(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk entropy given the temperature and the 
    // species mole fractions.  Also returns the species entropies.
    // STL version.
    double CalcBulkS(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &S        // Output vector for species entropies.
        ) const;


    // DIMENSIONLESS ENTROPY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the bulk entropy given the temperature
    // and the species mole fractions.  STL version.
    double CalcBulkS_R(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the bulk entropy given the temperature
    // and the species mole fractions.  C array version.
    double CalcBulkS_R(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk entropy given the temperature and the 
    // species mole fractions.  Also returns the species entropies.
    // STL version.
    double CalcBulkS_R(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &S        // Output vector for species entropies.
        ) const;


    // GIBBS FREE ENERGY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the bulk Gibbs free energy given the temperature
    // and the species mole fractions. STL version.
    double CalcBulkG(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the bulk Gibbs free energy given the temperature
    // and the species mole fractions. C array version.
    double CalcBulkG(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk Gibbs free energy given the temperature and the 
    // species mole fractions.  Also returns the species Gibbs free energies.
    // STL version.
    double CalcBulkG(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &G        // Output vector for species Gibbs free energies.
        ) const;


    // DIMENSIONLESS GIBBS FREE ENERGY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the dimensionless bulk Gibbs free energy given the temperature
    // and the species mole fractions. STL version.
    double CalcBulkG_RT(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the dimensionless bulk Gibbs free energy given the temperature
    // and the species mole fractions. C array version.
    double CalcBulkG_RT(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the dimensionless bulk Gibbs free energy given the temperature and the 
    // species mole fractions.  Also returns the species Gibbs free energies.
    // STL version.
    double CalcBulkG_RT(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &Gs_RT    // Output vector for species Gibbs free energies.
        ) const;


    // CONSTANT PRESSURE HEAT CAPACITY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the bulk constant pressure heat capacity given the
    // temperature and the species mole fractions.  STL version.
    double CalcBulkCp(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the bulk constant pressure heat capacity given the
    // temperature and the species mole fractions.  C array version.
    double CalcBulkCp(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk constant pressure heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.  STL version.
    double CalcBulkCp(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &Cp       // Output vector for species heat capacities.
        ) const;


    // DIMENSIONLESS CONSTANT PRESSURE HEAT CAPACITY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the dimensionless bulk constant pressure heat capacity given the
    // temperature and the species mole fractions.  STL version.
    double CalcBulkCp_R(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the dimensionless bulk constant pressure heat capacity given the
    // temperature and the species mole fractions.  C array version.
    double CalcBulkCp_R(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the dimensionless bulk constant pressure heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.  STL version.
    double CalcBulkCp_R(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &Cp       // Output vector for species heat capacities.
        ) const;


    // CONSTANT VOLUME HEAT CAPACITY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the bulk constant volume heat capacity given the
    // temperature and the species mole fractions.  STL version.
    double CalcBulkCv(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the bulk constant volume heat capacity given the
    // temperature and the species mole fractions.  C array version.
    double CalcBulkCv(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the bulk constant volume heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.  STL version.
    double CalcBulkCv(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &Cv       // Output vector for species heat capacities.
        ) const;


    // DIMENSIONLESS CONSTANT VOLUME HEAT CAPACITY (FUNCTIONS DERIVED FROM OTHERS).

    // Calculates the dimensionless bulk constant volume heat capacity given the
    // temperature and the species mole fractions.  STL version.
    double CalcBulkCv_R(
        double T,          // Temperature (K).
        const fvector &x // Species mole fractions.
        ) const;

    // Calculates the dimensionless bulk constant volume heat capacity given the
    // temperature and the species mole fractions.  C array version.
    double CalcBulkCv_R(
        double T,              // Temperature (K).
        const double *const x, // Species mole fractions.
        unsigned int n       // Number of values in x.
        ) const;

    // Calculates the dimensionless bulk constant volume heat capacity given 
    // the temperature and the species mole fractions.  Also 
    // returns the species heat capacities.  STL version.
    double CalcBulkCv_R(
        double T,           // Temperature (K).
        const fvector &x, // Species mole fractions.
        fvector &Cv       // Output vector for species heat capacities.
        ) const;


private:
    // Temporary vector of reals for doing internal calculations.  We define this
    // here in order to prevent continuous reallocation of the memory every time
    // a function which requires it is called.  It is defined as mutable as the
    // data it stores is not important to the ThermoInterface and this allows it
    // to be manipulated in const member functions.
    mutable fvector m_tmpvec;
};
};
};

#endif
