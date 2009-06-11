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
#include "swp_actsites_type.h"
#include "swp_process.h"
#include "swp_inception.h"
#include "swp_particle_process.h"
#include "swp_coagulation.h"
#include "sprog.h"
#include <vector>
#include <string>
#include <iostream>

#include "swp_ensemble_image.h"
#include "swp_particle_image.h"



namespace Sweep
{
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


    // ACTIVE-SITES MODELS.

    // Returns the set of particle model ID used by this mechanism
    const ActSites::ActSitesTypeSet &ActSiteModels(void) const;

    // Returns true if the mechanism include the given model.
    bool ContainsActSiteModel(ActSites::ActSitesType id) const;

    // Adds an active-sites model to the mechanism.
    void AddActSitesModel(ActSites::ActSitesType id);


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
    real CalcRates(
        real t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        fvector &rates,  // Return vector for process rates.
        bool scale=false // Scale the rates to the ensemble (=true), leave per unit vol (=false).
        ) const;

    // Get rates of all processes separated into different
    // terms.  Rate terms are useful for subsequent particle
    // selection by different properties for the same process.
    // In particular this is used for the condensation and
    // coagulation processes.  Returns the sum of all rates.
    real CalcRateTerms(
        real t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        fvector &terms   // Return vector for process rates.
        ) const;

/*
    // Get total rates of all processes.  Returns the sum of
    // all rates.  Uses supplied gas-phase conditions rather than
    // those in the given system.
    real CalcRates(
        real t,          // Time at which to get rates.
        const Sprog::Thermo::IdealGas &gas, // Use these gas-phase conditions instead.
        const Cell &sys, // System cell for which to get rates.
        fvector &rates   // Return vector for process rates.
        ) const;
*/

    // Get total rates of non-deferred processes.  Returns the sum
    // of all rates.
    real CalcJumpRateTerms(
        real t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        fvector &rates   // Return vector for process rates.
        ) const;

/*
    // Get total rates of non-deferred processes.  Returns the sum
    // of all rates.  Uses supplied gas-phase conditions rather than
    // those in the given system.
    real CalcJumpRates(
        real t,          // Time at which to get rates.
        const Sprog::Thermo::IdealGas &gas, // Use these gas-phase conditions instead.
        const Cell &sys, // System cell for which to get rates.
        fvector &rates   // Return vector for process rates.
        ) const;
*/

    // Calculates the rates-of-change of the chemical species fractions, 
    // gas-phase temperature and density due to particle processes.
    void CalcGasChangeRates(
        real t,          // Time at which to get rates.
        const Cell &sys, // System cell for which to get rates.
        fvector &rates   // Return vector for rates-of-change.
        ) const;


	// PERFORMING THE PROCESSES.

    // Performs the Process specified.  Process index could be
    // an inception, particle process or a coagulation event.
	void DoProcess(
        unsigned int i, // Index of process to perform.
        real t,         // Current time (s).
        Cell &sys       // System to update (includes ensemble).
        ) const;


    // LINEAR PROCESS DEFERMENT ALGORITHM.

    // Performs linear update algorithm on the 
    // given system up to given time.
    void LPDA(
        real t,   // Time up to which to integrate.
        Cell &sys // System to update.
        ) const;

/*
    // Performs linear update algorithm on the given system up to given time,
    // with the given chemical conditions rather than those in the 
    // given system.
    void LPDA(
        real t,    // Time up to which to integrate.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        Cell &sys // System to update.
        ) const;
*/

    // Performs linear process updates on a particle in the given system.
    void UpdateParticle(
        Particle &sp, // Particle to update.
        Cell &sys,    // System to which the particle belongs.
        real t        // Time up to which to integrate.
        ) const;

/*
    // Performs linear process updates on a particle in the given system,
    // with the current chemical conditions precalculated.
    void UpdateParticle(
        Particle &sp, // Particle to update.
        const Sprog::Thermo::IdealGas &gas, // Gas-phase conditions.
        Cell &sys,    // System to which the particle belongs.
        real t        // Time up to which to integrate.
        ) const;
*/

    // READ/WRITE/COPY.

    // Creates a copy of the mechanism.
    Mechanism *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);


	void  output(Cell &sys, real t) const;
	void  Mill(Cell &sys, real t) const;


private:
    // True if the mechanism contains deferred (LPDA)
    // processes, otherwise false. 
    mutable bool m_anydeferred;

    // Set of active-sites models which are implemented in this
    // mechanism.
    ActSites::ActSitesTypeSet m_actsites;

    // Processes in mechanism.
    Processes::IcnPtrVector m_inceptions;     // Inception process list.
    Processes::PartProcPtrVector m_processes; // Particle process list.
    Processes::Coagulation *m_coag;           // Coagulation process.

    // Auxilliary information about the processes.
    int m_icoag;                 // Index of first coagulation process in mechanism.
    unsigned int m_termcount;    // The rate term count of all processes.
    unsigned int m_processcount; // The process count.

    // Process counters.
    mutable std::vector<unsigned int> m_proccount, m_fictcount; 

    // Clears the mechanism from memory.
    void releaseMem(void);

    //time when the last particle emsenble has been written to a file
	mutable real last3dout;


	
};
};
#endif
