/*!
 * @file    swp_silicon_inception.cpp
 * @author  William J Menz
 * @brief   Defines inception of silicon particles by homogeneous nucleation.
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      Based on swp_dimer_inception by Robert Patterson, Markus Sander
 *      and Matthew Celnik.
 *
 *      This file uses the same rate calculation expressions as DimerInception
 *      to model the inception of particles. However, the silicon critical
 *      nucleus size is also dynamically calculated, and the incepting
 *      particle must have a diameter greater than this to yield a non-zero
 *      inception rate.
 *
 *      The idea of this is to model inception as a homogeneous nucleation
 *      process, where the critical diameter is calculated via the Kelvin
 *      equation.
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
 *      sweepc is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public License
 *      as published by the Free Software Foundation; either version 2
 *      of the License, or (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this program; if not, write to the Free Software
 *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *      02111-1307, USA.
 *
 *   Contact:
 *      Prof Markus Kraft
 *      Dept of Chemical Engineering
 *      University of Cambridge
 *      New Museums Site
 *      Pembroke Street
 *      Cambridge
 *      CB2 3RA, UK
 *
 *      Email:       mk306@cam.ac.uk
 *      Website:     http://como.cheng.cam.ac.uk
*/


#ifndef SWEEP_SILICON_INCEPTION_H
#define SWEEP_SILICON_INCEPTION_H

#include <iostream>
#include "swp_inception.h"

namespace Sweep {

// Forward declare the Mechanism class.
class Mechanism;

// Forward declare interface to gas phase data
class EnvironmentInterface;

namespace Transport
{
    // Forward declaration of unused argument type
    struct TransportOutflow;
}

namespace Processes
{

class SiliconInception: public Inception {
public:
    //! Initialising constructor.
    SiliconInception(const Sweep::Mechanism &mech);

    //! Copy constructor.
    SiliconInception(const SiliconInception &copy);

    //! Stream-reading constructor.
    SiliconInception(
        std::istream &in,                    //  - Input stream.
        const Sweep::Mechanism &mech         //  - Parent mechanism.
        );

    //! Destructor.
    ~SiliconInception(void);

    //! Operators.
    SiliconInception &operator=(const SiliconInception &rhs);


    // TOTAL RATE CALCULATIONS.

    //! Returns rate of the process for the given system.
    double Rate(
        double t,          // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d& local_geom // Information regarding surrounding cells and boundaries
        ) const;

    // PERFORMING THE PROCESS.

    //! Performs the process on the given system.
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng) const;


    // RATE TERM CALCULATIONS.

    //! Returns the number of rate terms for this process.
    unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a double vector. The
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all terms.
    double RateTerms(
        double t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        const Geometry::LocalGeometry1d &local_geom,                  // position information
        fvector::iterator &iterm // Iterator to the first term.
        ) const;


    // Sets the coagulation kernel constants given incepting species
    // masses and diameters.
    void SetInceptingSpecies(
        double m1, // Mass of species 1.
        double m2, // Mass of species 2.
        double d1, // Collision diameter of species 1.
        double d2  // Collision diameter of species 2.
        );

    //! Inception only according to free mol kernel with these sizes
    void SetInceptingSpeciesFreeMol(double m1, double m2, double d1, double d2);

    //! Sets the volume of an incepting particle
    void SetInceptingVolume(const Sweep::Mechanism &mech);

    //! Sets the diameter of an incepting particle
    void SetInceptingDiameter(const Sweep::Mechanism &mech);

    void GenerateSpeciesData(const Sweep::Mechanism &mech);

    //! Calculates the precursor concentration
    double GetPrecursorFraction(const Cell &sys) const;

    //! Returns the supersaturation
    double GetSupersaturation(const Sweep::Cell &sys) const;

    //! Calculates the critical nucleus size
    double GetCriticalNucleus(const Cell &sys) const;

    //! Checks whether inception should proceed, based on dcrit
    bool IsInceptionAllowed(const Cell &sys) const;

    //! Enum for inception equation type
    enum InceptionType {iCollisional, iVBDZ, iGirshick};

    //! Sets the inception mechanism
    void SetInceptionMechanism(InceptionType itype) {m_itype = itype;}


    //! Definition of the Species Data struct
    struct SiliconData
    {
        //! The index at which the silicon species is found in the gas-phase
        unsigned int _fracIndex;

        //! The name of the species
        std::string _name;

        //! The particle type-space if that species is selected
        fvector _track;

        //! The incepting diameter of the species
        double _diam;

        //! Default constructor
        SiliconData() :
            _fracIndex(0), _name(" "), _track(0), _diam(0.0) {}
    };

    // READ/WRITE/COPY.

    //! Copy function
    SiliconInception *const Clone(void) const;

    //! Returns the process type ID
    ProcessType ID(void) const;

    //! Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    //! Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    // Default constructor is protected to prevent an inception being
    // defined without knowledge of the parent mechanism.
    SiliconInception(void);

    // A faster rate calculation routine for Inception events only.  Requires all the
    // parameters that would otherwise be calculated by the routine to be passed as
    // arguments.
    double Rate(
        const EnvironmentInterface &gas,
        double sqrtT,           // Square-root of temperature.
        double T_mu,            // T / viscosity
        double MFP,             // Gas mean free path.
        double vol,             // Particle ensemble sample volume.
        const Cell &sys       // System for which to calculate rate terms.
        ) const;


    //! Calculates the gas-phase chemistry contribution to the rate expression
    double chemRatePart(
        const EnvironmentInterface &gas
        ) const;

private:

    //! A debugging function to print objects
    void printVector() const {
        for (unsigned int i(0); i!=m_sidata.size(); i++) {
            printSiliconData(m_sidata[i]);
        }
    }

    //! A debugging function to print the silicon data structure
    void printSiliconData(const SiliconData &data) const {
        std::cout << "nam: " << data._name;
        std::cout << " ind: " << data._fracIndex;
        std::cout << " dim: " << data._diam;
        for (unsigned int i(0); i!=data._track.size(); i++) {
            std::cout << " c" << i << ": " << data._track[i];
        }
        std::cout << std::endl;
    }

    //! Given a species name, indicate if it can incept
    bool IsCandidate(std::string name) const;

    //! Choses a species from the gas-phase for inception
    const SiliconData* ChooseData(const Sweep::Cell &sys, rng_type &rng) const;

    //! Adjusts the gas-phase based on the silicon species data
    void adjustGasPhase(Sweep::Cell &sys,
            const SiliconData &species, double wt) const;

    //! Returns the surface energy of silicon
    double GetSurfaceEnergy(double T) const;

    //! Returns the saturation vapour pressure of silicon
    double GetSatVapourPressure(double T) const;

    //! Returns the saturated monomer number concentration
    double GetMonomerConc(double T) const;

    // Rate parameters.
    //! Free-molecular kernel parameter.
    double m_kfm;
    double m_ksf1, m_ksf2; //! Slip-flow kernel parameters.

    //! Free-molecular enhancement factor.
    const double m_efm;

    //! Inception equation type
    InceptionType m_itype;

    //! Volume of a monomer (Si atom), hardcoded
    static const double m_v1;

    //! Mass of a monomer (Si atom), hardcoded
    static const double m_m1;

    //! Diameter of a monomer (Si atom), hardcoded
    static const double m_d1;

    //! Volume of an incepting particle
    double m_vi;

    //! Diameter of an incepting particle
    double m_di;

    //! A container for the silicon species data
    std::vector <SiliconData> m_sidata;

    //! Counter of number of reaction events for SiliconData index
    mutable std::vector<unsigned int> m_reacs;

    //! Counter of amount of concentration lost to an inception
    mutable std::vector<double> m_concs;
};

} // Processes

} // Sweep
#endif /* SILICONINCEPTION_H_ */
