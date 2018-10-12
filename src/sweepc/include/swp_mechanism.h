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

	// COAGULATION PROCESS WEIGHTED.

	// aab64 Returns TRUE if coagulation process uses weighted transfer function.
	bool IsWeightedCoag(void) const;

	// aab64 Sets the coagulation process to be SWA or not. 
	// Note that this is only for original activation, it is not meant to change the state
	// during simulation and does not provide a means of doing so. 
	virtual void SetWeightedCoag(bool weightedCoag);

	// aab64 Returns TRUE if inception process uses variable weights.
	bool IsVariableWeightedInception(void) const;

	// aab64 Returns variable inception max weight.
	double GetMaxInceptionWeight(void) const;

	// aab64 Returns variable inception min weight.
	double GetMinInceptionWeight(void) const;

	// aab64 Returns minimum particles threshold to start
	// adjusting incepting weight
	double GetMinSPForAIWOnset(void) const;

	// aab64 Returns the type of inception weight scaling function
	void GetWeightScalingFn(std::string &weightfn) const;

	// aab64 Sets flag for the inception process to use variable weighting.
	// Weights fluctuate between wmax and wmin depending on number of 
	// particles in ensemble relative to ensemble capacity.
	virtual void SetVariableWeightedInception(bool isVarInceptWeight, double wmax, double wmin, double nmin, std::string &weightfn);

	// aab64 Set flag for heavy inceptions
	virtual void SetIsHeavy(bool heavyflag, double upperdlimval, double lowerdlimval);

	// aab64 Get flag for heavy inceptions
	bool GetIsHeavy(void) const;

	// aab64 Get onset value for heavy inceptions
	double GetHeavyValue(void) const;

	// aab64 Get cutoff value for heavy inceptions
	double GetHeavyCutoffValue(void) const;

	// aab64 Set flag and onset value for surface inceptions
	virtual void SetIsSurfInc(bool surfincflag, double upperdlimval, double lowerdlimval, std::string &psitype);

	// aab64 Get flag for surface inceptions
	bool GetIsSurfInc(void) const;

	// aab64 Get onset value for surface inceptions
	double GetSurfIncValue(void) const;

	// aab64 Get cutoff value for surface inceptions
	double GetSurfIncCutoffValue(void) const;

	// aab64 Get psi type
	void GetPSItype(std::string &psitype) const;

	// aab64 Returns threshold to adjust ensemble weights
	double GetWeightOnsetRatio(void) const;

	// aab64 Returns factor to adjust ensemble weights
	double GetWeightScalingFactor(void) const;

	// aab64 Returns the flag for adaptive ensemble weights
	bool GetWeightScalingFlag(void) const;

	// aab64 Sets flag for the adaptive ensemble weights and the onset ratio
	virtual void SetWeightScaling(bool isWeightScaling, double ratio, double factor);

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

    //aab64 Add version to return concentration and fraction rates
    // Calculates the rates-of-change of the chemical species fractions, 
    // gas-phase temperature and density due to particle processes.
    void CalcGasChangeRates(
	double t,          // Time at which to get rates.
	const Cell &sys, // System cell for which to get rates.
	const Geometry::LocalGeometry1d& local_geom, // Information regarding surrounding cells
	fvector &xrates,   // Return vector for rates-of-change.
	fvector &crates   // Return vector for rates-of-change.
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

	//! Moment update for the incepting class - use a hybrid 
    //! method to track surface growth, assuming a lognormal 
    //! diameter distribution for the single primaries 
	void UpdateSections(
		double t,   // Time up to which to integrate.
		double dt,
		Cell &sys,// System to update.
		rng_type &rng
		) const;

    //! Set properties of particle picked for coagulation/outflow using 
    //! distribution parameters
	unsigned int SetRandomParticle(
		bool isSP1,
		bool isSP2,
		Sweep::Ensemble &ens,
		double t,
		double random_number, 
		Sweep::PropID prop,
		rng_type &rng) const;

    //! LPDA for one particle
    void UpdateParticle(
        Particle &sp, // Particle to update.
        Cell &sys,    // System to which the particle belongs.
        double t,       // Time up to which to integrate.
		int ind,        // Index of particle in the emsemble
        rng_type &rng, 
		PartPtrVector &overflow
        ) const;


    /* aab64 these two functions could replace the UpdateParticle 
    function when it is called in LPDA, to allow OpenMP to be used
    in the sintering part of the update - this needs to be checked */ 
    //! non-sintering part of LPDA for one particle
    void UpdateParticleNS(
	Particle &sp, // Particle to update.
	Cell &sys,    // System to which the particle belongs.
	double t,       // Time up to which to integrate.
	rng_type &rng, 
	double &dtvec
	) const;

    //! sintering part of LPDA for one particle
    void UpdateParticleS(
	Particle &sp, // Particle to update.
	Cell &sys,    // System to which the particle belongs.
	double t,       // Time up to which to integrate.
	rng_type &rng,
	double dtvec
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
	


//////////////////////////////////////////// aab64 ////////////////////////////////////////////
    // Get the number of times each fictitious process has been performed 
    std::vector<unsigned int> GetFictitiousProcessCounts() const { return m_fictcount; }

    // Get the term count
    unsigned int GetTermCount() const { return m_termcount; }
	
    // Get the addition count
    unsigned int GetDeferredAddCount() const { return m_addcount; }
	
    // Get the inflow count
    unsigned int GetInflowCount() const { return m_inflowcount; }

    // Get the outflow count
    unsigned int GetOutflowCount() const { return m_outflowcount; }
    
	// aab64 particle number hybrid parameters
	void SetHybrid(bool hybrid_flag) const { m_hybrid = hybrid_flag; }
	bool IsHybrid() const { return m_hybrid; }
	void SetCriticalThreshold(unsigned int threshold) const { m_threshold = threshold; } 
	unsigned int GetCriticalThreshold() const { return m_threshold; }
//////////////////////////////////////////// aab64 ////////////////////////////////////////////
	


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



//////////////////////////////////////////// aab64 ////////////////////////////////////////////
	mutable unsigned int m_addcount;        // The addition count for deferred additions
	mutable unsigned int m_inflowcount;     // The inflow count for stochastic inflow
	mutable unsigned int m_outflowcount;    // The outflow count for stochastic inflow

	mutable bool m_weighted_coag;           // Is the coagulation process one of the weighted ones?
	mutable bool m_var_incept_weight;       // Is variable inception weighting turned on?
	mutable double m_minsp_for_aiw;         // Minimum particle threshold for which inception weight should be adapted
	mutable double m_min_incept_weight;     // Minimum incepting weight, corresponding to nmin
	mutable double m_max_incept_weight;     // Maximum incepting weight, corresponding to Nmax SPs
	mutable std::string m_incept_weight_fn; // The type of inception weight scaling to use

	mutable bool m_heavyallowed;            // Flag to allow heavier inception particles
	mutable double m_upp_dval_heavy;        // Onset for heavier inceptions
	mutable double m_low_dval_heavy;        // Cutoff for heavier inceptions
	mutable bool m_surfincflag;             // Flag to allow surface inceptions
	mutable double m_upp_dval_surfinc;      // Onset for surface inception
	mutable double m_low_dval_surfinc;      // Cutoff for surface inception
	mutable std::string m_psi_type;         // Type of particle surface inception to do

	mutable bool m_weightscaling_flag;      // Flag for adaptive ensemble weights
	mutable double m_weightscaling_onset;   // Onset ratio for adaptive ensemble weights
	mutable double m_weightscaling_factor;  // Factor multiplying N/sum(w) in weight scaling

	mutable bool m_hybrid;                  // identify hybrid particle model
	mutable unsigned int m_threshold;       // hybrid threshold value
//////////////////////////////////////////// aab64 ////////////////////////////////////////////



    // Clears the mechanism from memory.
    void releaseMem(void);

};
} // namespace Sweep
#endif
