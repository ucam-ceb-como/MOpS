/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Stochastic mechanism for Sweep.  The mechanism holds all the processes
    which can be enacted on a system with a particle ensemble.  It also holds
    auxilliary info which defines how those processes work.

    The mechanism class defines routines for calculating process rates and perform
    a process.  It also contains the Linear Process Deferment Algorithm (LPDA) for
    updating linear processes that have been removed from the stochastic jump process.

    The mechanism also defines particle components and additional particle values.  Particle
    values are particle properties which can be changed by processes, but which are not 
    components, therefore they do not contribute to particle mass.

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

#ifndef SWEEP_MECHANISM_H
#define SWEEP_MECHANISM_H

#include "swp_params.h"
#include "swp_particle_model.h"
#include "swp_cell.h"
#include "swp_ensemble.h"
#include "swp_process.h"
#include "swp_inception.h"
#include "swp_particle_process.h"
#include "swp_coagulation.h"
#include "swp_fragmentation.h"

#include <vector>
#include <string>
#include <iostream>

// Forward declaration
namespace Geometry {
    class LocalGeometry1d;
}

namespace Sweep
{

/*!
 *@brief The Mechanism collects together all of the processes which change the system
 *
 * It is not entirely clear why Mechanism must inherit from ParticleModel
 */
class Mechanism : public ParticleModel
{
public:
	// Constructors.
	Mechanism(void);                  // Default Constructor.
	Mechanism(const Mechanism &copy); // Copy-Constructor.

    // Destructor.
    ~Mechanism(void);

    // Operators.
    Mechanism &operator=(const Mechanism &rhs);

	// INCEPTIONS.

    // Returns the vector of inceptions.
    const Processes::IcnPtrVector &Inceptions(void) const;

    // Returns the inception with the given index.
    const Processes::Inception *const Inceptions(unsigned int i) const;

    // Adds an inception to the mechanism.
    void AddInception(Processes::Inception &icn);


    // PARTICLE PROCESSES.

    // Returns the vector of particle processes.
    const Processes::PartProcPtrVector &Processes(void) const;

    // Returns the process with the given index.
    const Processes::ParticleProcess *const Processes(unsigned int i) const;

    // Adds a process to the mechanism.
    void AddProcess(Processes::ParticleProcess &p);

    // COAGULATIONS.

    // Adds a coagulation process to the mechanism.
    void AddCoagulation(Processes::Coagulation &coag);
    
    //! Returns the vector of coagulations.
    const Processes::CoagPtrVector &Coagulations(void) const;

    // Adds a coagulation process to the mechanism.
    void AddFragmentation(Processes::Fragmentation &frag);
    
    //! Returns the vector of coagulations.
    const Processes::FragPtrVector &Fragmentations(void) const;

    // PROCESS INFORMATION.

    // Returns the number of processes (including 
    // inceptions) in the mechanism.
    unsigned int ProcessCount(void) const;

    // Returns the number of terms in all process rate expressions.
    unsigned int TermCount(void) const;

    // Returns true if the mechanism contains deferred (LPDA) 
    // processes otherwise false.
    bool AnyDeferred(void) const;

    // Checks all processes to see if any are deferred.
    void CheckDeferred(void) const;

    // Returns a vector containing the names of all processes.
    void GetProcessNames(
        std::vector<std::string> &names, // Output vector for names.
        unsigned int start=0             // Optional vector start index.
        ) const;

	// RATE CALCULATION.

    // Get total rates of all processes.  Returns the sum of
    // all rates.
    double CalcRates(
        double t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        const Geometry::LocalGeometry1d& local_geom, // Information regarding surrounding cells
        fvector &rates,  // Return vector for process rates.
        bool scale=false // Scale the rates to the ensemble (=true), leave per unit vol (=false).
        ) const;

    //! Get total number of jumps of all processes
    double CalcJumps(
        double t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        const Geometry::LocalGeometry1d& local_geom, // Information regarding surrounding cells
        fvector &jumps  // Return vector for process rates.
        ) const;

    //! Reset the jump number vectors
    void ResetJumpCount() const;

    // Get rates of all processes separated into different
    // terms.  Rate terms are useful for subsequent particle
    // selection by different properties for the same process.
    // In particular this is used for the condensation and
    // coagulation processes.  Returns the sum of all rates.
    double CalcRateTerms(
        double t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        const Geometry::LocalGeometry1d& local_geom, // Information regarding surrounding cells
        fvector &terms   // Return vector for process rates.
        ) const;

    // Get total rates of non-deferred processes.  Returns the sum
    // of all rates.
    double CalcJumpRateTerms(
        double t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        const Geometry::LocalGeometry1d& local_geom, // Information regarding surrounding cells
        fvector &rates   // Return vector for process rates.
        ) const;

    //! Rate of processes that are deferred
    double CalcDeferredRateTerms(
        double t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        const Geometry::LocalGeometry1d& local_geom, // Information regarding surrounding cells
        fvector &rates   // Return vector for process rates.
        ) const;

    // Calculates the rates-of-change of the chemical species fractions, 
    // gas-phase temperature and density due to particle processes.
    void CalcGasChangeRates(
        double t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        const Geometry::LocalGeometry1d& local_geom, // Information regarding surrounding cells
        fvector &rates   // Return vector for rates-of-change.
        ) const;


	// PERFORMING THE PROCESSES.

    // Performs the Process specified.  Process index could be
    // an inception, particle process or a coagulation event.
    void DoProcess(
        unsigned int i, // Index of process to perform.
        double t,         // Current time (s).
        Cell &sys,      // System to update (includes ensemble).
        const Geometry::LocalGeometry1d& local_geom, // Information regarding surrounding cells
        rng_type &rng
        ) const;

    //! Perform the transport processes in and out of the cell.
    void DoParticleFlow(
            double t,
            double dt,
            Cell &sys,
            const Geometry::LocalGeometry1d& local_geom,
            rng_type &rng) const;

    //! Transfer mass from gas phase to particle ensemble used for PAH-PP model.
    void MassTransfer(
        int i,                                         //!< The number of pyrene supposed in the ensemble.
        double t,                                      //!< Current time (s).
        Cell &sys,                                     //!< System to update (includes ensemble).
        rng_type &rng,                                 //!< Random number generator.
        const Geometry::LocalGeometry1d& local_geom    //!< Information regarding surrounding cells.
        ) const;


    // LINEAR PROCESS DEFERMENT ALGORITHM.

    //! LPDA for all particles
    void LPDA(
        double t,   // Time up to which to integrate.
        Cell &sys,// System to update.
        rng_type &rng
        ) const;


    //! LPDA for one particle
    void UpdateParticle(
        Particle &sp, // Particle to update.
        Cell &sys,    // System to which the particle belongs.
        double t,       // Time up to which to integrate.
		int ind,        // Index of particle in the emsemble
        rng_type &rng, 
		PartPtrVector &overflow
        ) const;

    // READ/WRITE/COPY.

    // Creates a copy of the mechanism.
    Mechanism *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

    //! Get the number of times each process has been performed
    std::vector<unsigned int> GetProcessUsageCounts() const {return m_proccount;}
	
	//! return a vector contain the information of particular primary particle with X molecules
	void Mass_pah(Ensemble &m_ensemble) const;

	//! write data in colunm for dimer and mononer
	//void writeMononer(fvector &out) const;
	//void writeDimer(fvector &out) const;
	
	//! write data in colum for particlar primary particle
	//void writeParimary(std::vector<fvector > &out) const;

private:
    // True if the mechanism contains deferred (LPDA)
    // processes, otherwise false. 
    mutable bool m_anydeferred;

    // Processes in mechanism.
    Processes::IcnPtrVector m_inceptions;     // Inception process list.
    Processes::PartProcPtrVector m_processes; // Particle process list.

    //! List of coagulation processes
    Processes::CoagPtrVector m_coags;

    //! List of coagulation processes
    Processes::FragPtrVector m_frags;

    // Auxilliary information about the processes.
    int m_icoag;                 // Index of first coagulation process in mechanism.
    unsigned int m_termcount;    // The rate term count of all processes.
    unsigned int m_processcount; // The process count.

    // Process counters.
    mutable std::vector<unsigned int> m_proccount, m_fictcount; 

    // Clears the mechanism from memory.
    void releaseMem(void);

};
} // namespace Sweep
#endif
