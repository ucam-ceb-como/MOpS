/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    A specialised process for Sweep, which acts on Sweep::Cell
    objects.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

/*!
 *\file swp_process.h
 *\brief Interface declaration for process classes
 */

#ifndef SWEEP_PROCESS_H
#define SWEEP_PROCESS_H

#include "swp_params.h"
#include "swp_particle.h"
#include "swp_process_type.h"

#include "gpc_stoich.h"

#include <map>
#include <string>
#include <iostream>

namespace Sweep
{
// Forward declare parent mechanism.
class Mechanism;
// Forward declare the Cell class.
class Cell;

namespace Transport
{
    // Forward declaration of unused argument type
    struct TransportOutflow;
}

namespace Processes
{
class Process
/*!
 * \brief Define an interface for all processes
 *
 * Processes are the principal means for specifying simulation dynamics.
 */
{
public:
	/// Constructors.
    Process(const Sweep::Mechanism &mech); // Default Constructor.
    Process(const Process &copy);          // Copy-constructor.
    Process(                         // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    // Destructor.
    virtual ~Process(void);

    // Operators.
    Process &operator=(const Process &rhs);


    // PROCESS NAME/DESCRIPTION.

    // Returns the process name.
    const std::string &Name(void) const;

    // Sets the process name.
    void SetName(const std::string &name);


	// PARENT MECHANISM.

    // Returns reference to parent mechanism.
    const Sweep::Mechanism *const Mechanism() const;

    // Sets the parent mechanism
    void SetMechanism(const Sweep::Mechanism &mech);

    // VISCOSITY MODEL.

    //! Returns the viscosity model
    Sweep::ViscosityModel ViscosityModel() const;

    //! Sets the process' viscosity model
    void SetViscosityModel(Sweep::ViscosityModel vmodel);

    //! Returns the viscosity for the given model
    real GetViscosity(const Cell &sys) const;

    // REACTANTS.

    // Returns the number of reactants.
    unsigned int ReactantCount() const;

    // Returns the stoichiometric reactant coefficients.
    const Sprog::StoichMap &Reactants(void) const; 

    // Returns the stoichiometry of the ith reactant.
    int Reactants(unsigned int i) const;

    // Adds a reactant to the reaction.
    void AddReactant(
        unsigned int isp, // Species index.
        int mu            // Stoichiometry value.
        );

    // Adds a reactant given the species name.
    void AddReactant(
        const std::string &name, // Species name.
        int mu                   // Stoichiometry value.
        );

    // Removes a reactant, given by name, from the reaction.
    void RemoveReactant(const std::string &name);


    // PRODUCTS.

    // Returns the number of reactants.
    unsigned int ProductCount() const;

    // Returns the stoichiometric reactant coefficients.
    const Sprog::StoichMap &Products(void) const; 

    // Returns the stoichiometry of the ith reactant.
    int Products(unsigned int i) const;

    // Adds a reactant to the reaction.
    void AddProduct(
        unsigned int isp, // Species index.
        int mu            // Stoichiometry value.
        );

    // Adds a reactant given the species name.
    void AddProduct(
        const std::string &name, // Species name.
        int mu                   // Stoichiometry value.
        );

    // Removes a reactant, given by name, from the reaction.
    void RemoveProduct(const std::string &name);


	// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

    //! Rate of the process for the given system.
    virtual real Rate(
        real t,          // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d& local_geom
        ) const = 0;

	// RATE TERM CALCULATIONS.
    //   These routines return the individual rate terms for a 
    //   process, which may have multiple terms (e.g. condensation).

    //! Returns the number of rate terms for this process.
    virtual unsigned int TermCount(void) const = 0;

    // Calculates the rate terms given an iterator to a real vector. The 
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all rate terms.
    virtual real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        const Geometry::LocalGeometry1d& local_geom,
        fvector::iterator &iterm // Iterator to the first term.
        ) const = 0;

    // PERFORMING THE PROCESS.

    /*!
     * \brief Performs the process on the given system.
     *
     * \param[in]       t           Time
     * \param[in,out]   sys         System to update
     * \param[in]       local_geom  Details of local phsyical layout
     * \param[in]       iterm       Process term responsible for this event
     * \param[in,out]   rng         Random number generator
     *
     * \return      0 on success, otherwise negative.
     */
    virtual int Perform(
        real t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng
        ) const = 0;


    // FICTICIOUS EVENTS.

    //! See whether an event is fictitious
    static bool Fictitious(real majr, real truer, rng_type &rng);

    // READ/WRITE/COPY.

    // Returns a copy of the process
    virtual Process *const Clone(void) const = 0;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    virtual ProcessType ID(void) const = 0;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    // Process name/description.
    std::string m_name;

    // Process data.
    const Sweep::Mechanism *m_mech; // Pointer to the parent Mechanism.
    Sprog::StoichMap m_reac; // Reactant species stoichiometry.
    Sprog::StoichMap m_prod; // Product species stoichiometry.

    // Model for the process' viscosity calculation
    Sweep::ViscosityModel m_viscosity_model;

    // Default constructor is protected so that processes cannot
    // be defined without knowledge of the parent mechanism.
    Process(void);

    // Calculates the gas-phase chemistry contribution to the rate
    // expression.
    real chemRatePart(
        const fvector &fracs, // Species mole fractions in gas phase.
        real density          // Gas phase molar density.
        ) const;

    // Adjusts the gas-phase composition using the reactants and
    // products defined for this process.
    void adjustGas(
        Cell &sys,     // System to update.
        real wt,       // Stochastic weight of particle
        unsigned int n // Number of times to apply process.
         = 1           //  - Default is one time.
         ) const;
}; // class Process
} // namespace Processes
} // namespace Sweep

#endif
