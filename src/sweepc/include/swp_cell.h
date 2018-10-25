/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    A specialised ideal gas class for Sweep, which contains an ensemble
    with particles of type Particle.

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

#ifndef SWEEP_CELL_H
#define SWEEP_CELL_H

#include "swp_params.h"
#include "swp_particle_model.h"
#include "swp_ensemble.h"
#include "swp_ensemble_stats.h"
#include "swp_birth_process.h"
#include "swp_death_process.h"
#include "swp_environment_interface.h"

#include <string>
#include <iostream>

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

/*!
 *@brief Population particles in an ideal gas
 *
 * Also contains a pointer to a particle model, which is needed to interpret
 * the particles, but also brings in mechanism information, which should not
 * be in a cell, which capture state not dynamics.
 */
class Cell
{
public:
	// Constructors.
	Cell(const Sweep::ParticleModel &model, const bool const_gas = false); // Default constructor.
	Cell(const Cell &copy);                  // Copy constructor.
	Cell(                                 // Stream-reading constructor.
		std::istream &in,                 //   - Stream from which to read.
		const Sweep::ParticleModel &model //   - Model used to define particles.
		);

	// Destructor.
	virtual ~Cell(void);

	//! Overwrite contents
	Cell &operator=(const Cell &rhs);

	//! Overload of the << operator
	friend std::ostream& operator<<(
		std::ostream &os,
		const Sweep::Cell &c);

	// THE GAS-PHASE INTERFACE.

	//!Returns the description of the gas-phase mixture.
	const EnvironmentInterface &GasPhase(void) const { return *m_gas; };
	//!Returns the description of the gas-phase mixture.
	EnvironmentInterface &GasPhase(void) { return *m_gas; };

	// THE PARTICLE ENSEMBLE.

	// Returns the particle ensemble.
	Sweep::Ensemble &Particles(void);
	const Sweep::Ensemble &Particles(void) const;

	// Returns the number of particles in the ensemble.
	unsigned int ParticleCount(void) const;

	// Returns the number of particles in the ensemble.
	double ParticleWeightSum(void) const;

    // Returns particle statistics.
    void GetVitalStats(Stats::EnsembleStats &stats) const;

	//! Initialise with some particles, downsampling as required
	void SetParticles(
		std::list<Particle*>::iterator particle_list_begin,
		std::list<Particle*>::iterator particle_list_end,
		double statistical_weight, rng_type &rng);

	//! Initialise with some particles, which must not exceed the ensemble capacity
	void SetParticles(
		std::list<Particle*>::iterator particle_list_begin,
		std::list<Particle*>::iterator particle_list_end,
		double statistical_weight);

	// THE PARTICLE MODEL.

	// Returns the particle model used to define particles in this
	// cell.
	const Sweep::ParticleModel *const ParticleModel(void) const { return m_model; }


	// SCALING ROUTINES INCL. SAMPLE VOLUME.

	//! Returns the double system to stochastic system scaling factor.
	double SampleVolume() const;

	//! Multiply the sample volume by a scaling factor
	void AdjustSampleVolume(double scale_factor);

	//! Empty the cell and set the sample volume so a full particle ensemble would have the specified m0
	void Reset(double m0);

	// FIXED/VARIABLE CHEMISTRY.

	// Returns whether or not the chemical conditions are fixed.
	bool FixedChem() const;

	// Sets whether or not the chemical conditions are fixed.
	void SetFixedChem(bool fixed = true);

	// Set the chemical conditions to be variable.
	void SetVariableChem(bool vari = true);


	// PARTICLE INFLOW PROCESSES.

	//! Returns the number of inflow processes defined for this Cell.
	unsigned int InflowCount(void) const;

	//! Returns the vector of particle inflow processes.
	const Processes::BirthPtrVector &Inflows(void) const;

	//! Returns the ith particle inflow process.
	Processes::BirthProcess *const Inflows(unsigned int i) const;

	//! Add an inflow process to the Cell. The process is copied.
	void AddInflow(Processes::BirthProcess &inf);

	// PARTICLE OUTFLOW PROCESSES.

	//! Returns the number of outflow processes defined for this Cell
	unsigned int OutflowCount(void) const;

	//! Returns the vector of particle outflow processes.
	const Processes::DeathPtrVector &Outflows(void) const;

	//! Returns the ith particle outflow process.
	Processes::DeathProcess *const Outflows(unsigned int i) const;

	//! Add an outflow process to the Cell. The process is copied.
	void AddOutflow(Processes::DeathProcess &out);

	// Add an outflow process with the given rate to the Cell.
	void AddOutflow(
		double rate, // Rate constant for outflow.
		const Sweep::Mechanism &mech // Mechanism which defines LPDA for outflow.
		);

	// figure out how many starting species is supposed to be in the particle ensemble
	unsigned int NumOfStartingSpecies(const int index) const;

	// READ/WRITE/COPY.

	// Writes the object to a binary stream.
	virtual void Serialize(std::ostream &out) const;

	// Reads the object from a binary stream.
	virtual void Deserialize(
		std::istream &in,                 // Input stream.
		const Sweep::ParticleModel &model // Model used to define particles.
		);

	// aab64 set particle temperature
	void SetBulkParticleTemperature(double ptemp) {
		m_bulk_particle_temp = ptemp;
	}

	// aab64 get particle temperature
	double GetBulkParticleTemperature() const { return m_bulk_particle_temp; }

	// aab64 Set the process time interval
	void SetCurrentProcessTau(double tau) { m_proc_tau = tau; }

	// aab64 Get the process time interval 
	double GetCurrentProcessTau() const { return m_proc_tau; }

	// aab64 Set the process time interval
	void SetIsAdiabaticFlag(bool flag) { m_adiabatic_flag = flag; }

	// aab64 Get the process time interval 
	bool GetIsAdiabaticFlag() const { return m_adiabatic_flag; }

	// aab64 Set the incepting particle weight
	void SetInceptingWeight(double wt) { m_incepting_weight = wt; }

	// aab64 Get the incepting particle weight 
	double GetInceptingWeight() const { return m_incepting_weight; }

	// aab64 Set the inception factor to change how many units are present
	void SetInceptionFactor(double incfac) { m_incFactor = incfac; }

	// aab64 Get the inception factor that determines how many units are added
	double GetInceptionFactor() const { return m_incFactor; }

	// aab64 Set PSI flag used to tell surface reaction not to update gas-phase and 
	// temperature because this is handled using inception stoichiometry
	void SetNotPSIFlag(bool psiflag) { m_notpsiflag = psiflag; }

	// aab64 Get PSI flag used to tell surface reaction not to update gas-phase and 
	// temperature because this is handled using inception stoichiometry
	bool GetNotPSIFlag() const { return m_notpsiflag; }

	// aab64 Store and access coagulation properties used to select particles for PSI
	void SetCoagProps(PropID prop1, PropID prop2) { m_cprop1 = prop1; m_cprop2 = prop2; };
	PropID getCoagProp1() const { return m_cprop1; };
	PropID getCoagProp2() const { return m_cprop2; };


	// aab64 coagulation scaling for weighted events
	void SetRateFactor(double rateFac) { m_rateFactor = rateFac; }

	// aab64 coagulation scaling for weighted events
	int GetRateFactor() const { return m_rateFactor; }

	// aab64 Temporary functions for gas-phase properties
	void setGasPhaseProperties(double cp_bulk, double rhop, fvector enthalpies);
	double getBulkHeatCapacity() const { return m_bulk_heat_capacity; }
	double getParticleDensity() const { return m_particle_density; }
	void getEnthalpies(fvector &enthalpies) { enthalpies = m_enthalpies; }

protected:
    // Default constructor is protected as it makes no
    // sense to define a mixture without knowledge of the
    // defining species.  This trait is brought over from Sprog.
    Cell();
	
	// Check internal consistency 
	bool isValid() const;

private:
    //! Gas mixture
    EnvironmentInterface *m_gas;

    //! Particle ensemble.
    Ensemble m_ensemble;

    //! Particle model.
    const Sweep::ParticleModel *m_model;

    //! The volume in which the ensemble represents the complete double system.
    double m_smpvol;

    // Flag determining whether or not the chemistry in this system is fixed.
    // If the chemical conditions are fixed, then they cannot be altered by
    // any particle processes.  Default is false.
    bool m_fixed_chem;

    // Particle inflow processes.  These are not stored in the
    // mechanism, but are used by the Mechanism class when
    // calculating rates.
    Processes::BirthPtrVector m_inflow;

    // Particle outflow processes.  These are not stored in the
    // mechanism, but are used by the Mechanism class when
    // calculating rates.
    Processes::DeathPtrVector m_outflow;

    // aab64 new particle temp
    double m_bulk_particle_temp;

    // aab64 The process time interval
    double m_proc_tau; 

	// aab64 Flag for adiabatic operation
	bool m_adiabatic_flag;

	// aab64 Current incepting particle weight in Bintree primary case
	double m_incepting_weight;

	// aab64 Scale the inception process to increase primary size more rapidly
	double m_incFactor;

	// aab64 PSI flag used to tell surface reaction not to update gas-phase and 
	// temperature because this is handled using inception stoichiometry
	bool m_notpsiflag;

	// aab64 coagulation scaling for weighted events
	int m_rateFactor;

	// aab64 temp
	PropID m_cprop1;
	PropID m_cprop2;

	double m_bulk_heat_capacity;
	double m_particle_density;

	fvector m_enthalpies;

};

} //namespace Sweep

#endif
