/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Primary class defines the smallest sub-particle unit.  It uses 
    the spherical particle model, though is designed so that derived
    classes can relax that assumption.  It stores the fundamental
    data structures required to describe different particle models.

    In the context of the sub-particle tree, primaries can be thought of
    as the leaves of the tree.  The SubParticle class describes the branches
    and the Particle class describes the trunk.

    A Primary particle subscribes to a particular ParticleModel object.  The
    particle model defines what constitutes a primary and what sub-models are
    enabled.  A Primary is constructed with knowledge of the particle model,
    and the particle model cannot be changed.  It is expected that the
    data in the ParticleModel object does not change once Primaries have
    been created (e.g. by adding/removing components), and this class is
    defined uner this assumption.

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

#ifndef SWEEP_PRIMARY_H
#define SWEEP_PRIMARY_H

#include "swp_params.h"
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_particle_model.h"
#include "swp_aggmodel_type.h"
#include "swp_sintering_model.h"
#include "swp_titania_melting_model.h"
#include "swp_property_indices.h"

#include <iostream>

namespace Sweep
{
	enum Xmer{ MOMOMER=1,DIMER=2,TRIMER=3};
// Forward declaration
class Cell;

namespace AggModels {

/*!
 *  \brief Base class for primary particles
 *
 *  Primary particles are the place in sweep where particle structure and composition
 *  are stored.  Each computational particle contains one primary particle, which
 *  in some cases serves as the head of a tree of further primary particles.  When the
 *  computational particle needs to know about its own physical properties it queries
 *  its (top-level) primary particle.
 *
 *  There are some weaknesses in the current design: Primary and its subclasses try to
 *  satisfy two roles, one is to capture aggregate structure and the second is represent
 *  the contents of a (piece of a) particle.
 */
class Primary
{
public:
    // Constructors.
    //   Note:  Default constructor is protected to prevent a
    //          Primary being created without knowledge of the
    //          defining particle model.
    Primary(                              // Initialising constructor.
        double time,                        //  - Create time.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );
    Primary(const Primary &copy); // Copy constructor.
    Primary(                              // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );

    // Destructors.
    virtual ~Primary(void);

    // Operators.
    virtual Primary &operator=(const Primary &rhs);

    // DEFINING PARTICLE MODEL.
    
    // Returns the particle model used to create this primary.
    const Sweep::ParticleModel *const ParticleModel(void) const;


    // PRIMARY COMPOSITION.

	//Returns named component
	double GetComponent(std::string name) const;

    // Returns the composition vector.
    const fvector &Composition(void) const;

    // Returns the ith component value.  Returns 0.0 if i invalid.
    double Composition(unsigned int i) const;

    // Sets the composition vector.
    void SetComposition(const fvector &comp);


    // PRIMARY TRACKER VALUES.

    // Returns the values vector.
    const fvector &Values(void) const;

    // Returns the ith value.  Returns 0.0 if i invalid.
    double Values(unsigned int i) const;

    // Sets the values vector.
    void SetValues(const fvector &vals);

    // Sets the ith trackervalue.
    void SetValue(unsigned int i, double val);

	// PHASE

	// Get mass of a specific phase
	double GetPhaseMass(int i) const;

    // PRIMARY CREATE TIME.

    // Returns the particle create time.
    double CreateTime(void) const;


    // LAST UPDATE TIME.
    //   This is provided to assist numerics only.  It is not available
    //   to other programs/libraries outside sweep.

    // Returns the last update time of the particle.
    double LastUpdateTime(void) const;

    // Sets the last update time of the particle.
    virtual void SetTime(double t);


    // AGGREGATION MODEL.

    // Returns the aggregation model which this primary describes.
    virtual AggModels::AggModelType AggID(void) const;

    // BASIC DERIVED PROPERTIES.

    // Calculates the derived properties from the unique properties.
    virtual void UpdateCache(void);

    // Returns the equivalent-sphere diameter.
    double SphDiameter(void) const;

    // Returns the collision diameter.
    double CollDiameter(void) const;

    // Returns the mobility diameter.
    double MobDiameter(void) const;

    // Returns the surface area.
    double SurfaceArea(void) const;

    // Returns the equivalent sphere surface area, based
    // on the volume.
    double SphSurfaceArea(void) const;

    // Returns the volume.
    double Volume(void) const;

    // Returns the mass.
    double Mass(void) const;

    //! Returns the number of carbon atoms.
    int NumCarbon(void) const;

    //! Returns fragmentation flag.
    int Frag(void) const;
	//! Returns the number of rings.
	int NumRings(void) const;

    //! Returns the property with the given ID.
    double Property(const Sweep::PropID id) const;


    // BASIC DERIVED PROPERTY OVERWRITES.

    // Sets the spherical particle diameter
    void SetSphDiameter(double diam);

    // Sets the collision diameter of the particle.
    void SetCollDiameter(double dcol);

    // Sets the mobility diameter.
    void SetMobDiameter(double dmob);

    // Sets the surface area, subject to minimum spherical area condition.
    void SetSurfaceArea(double surf);

    // Sets the volume.
    void SetVolume(double vol);

    // Sets the mass.
    void SetMass(double m);

    //! Sets the number of carbon atoms.
    void SetNumCarbon(int numcarbon);

    //! Sets fragmentation flag.
    void SetFrag(int frag);
    //! Check particle still meets physical conditions for being a particle.
    bool IsValid() const;

    // OPERATIONS.

    // Adjusts the primary with the given composition and 
    // tracker values changes n times.  If the particle cannot be adjust
    // n times, then this function returns the number of times
    // it was adjusted.
    virtual unsigned int Adjust(
        const fvector &dcomp,   // Composition changes.
        const fvector &dvalues, // Tracker variable changes.
        rng_type &rng,			// Random number for leaf node
        unsigned int n=1        // Number of times to perform adjustment.
        );

    // As for Adjust, but applies to SilicaPrimary's IntParticle reaction
    virtual unsigned int AdjustIntPar(
        const fvector &dcomp,   // Composition changes.
        const fvector &dvalues, // Tracker variable changes.
        rng_type &rng,			// Random number for leaf node
        unsigned int n=1        // Number of times to perform adjustment.
        );

	//Adjusts the particle n times for the phase transformation process
	virtual unsigned int AdjustPhase(
		const fvector &dcomp,	// Composition changes.
        const fvector &dvalues,	// Tracker variable changes.
        rng_type &rng,			// Random number for leaf node
        unsigned int n			// Number of times to perform adjustment.
		);

	//Melting point dependent phase change
	virtual void Melt( rng_type &rng, Cell &sys	);

    // Combines this primary with another.  This is also the
    // implementation of the + and += operators.
    virtual Primary &Coagulate(const Primary &rhs,
                               rng_type &rng);

    // Combines this primary with another.  This is also the
    // implementation of the + and += operators.
    virtual Primary &Fragment(const Primary &rhs,
                               rng_type &rng);

    // This routine sinters the Primary for the given length of
    // time using the provided sintering model.
    virtual void Sinter(
        double dt, // Delta-t for sintering to occur.
        Cell &sys, // System which defines primary's environment.
        const Processes::SinteringModel &model, // Sintering model to use.
        rng_type &rng,  // Random number generator
        double wt     // Statistical weight
        );
	
	// PARTICLE TRACKING FOR VIDEOS

	//! Returns the frame position and orientation, and primary coordinates
	//	Used by bintree_primary for particle tracking
	virtual void GetFrameCoords(std::vector<fvector> &coords) const;

	//! Set primary tracking
	virtual void setTracking();

	//! Remove primary tracking
	virtual void removeTracking();

    // READ/WRITE/COPY.

    // Returns a copy of the primary.
    virtual Primary *const Clone(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );

	//! Gets the number of active sites (always component 0)
	virtual double GetSites() const { return m_comp[0]; }

	// Return the sintering rate for a SilicaPrimary
	virtual double GetSintRate() const { return 0.0; }

	virtual double GetCoverageFraction() const;

    //! Check whether the number of carbon atoms in the primary is equal to
    //! that of the inception species.
    int InceptedPAH() const;

	// Term for titania phase transformation
	virtual double GetPhaseTerm(void) const;

protected:
    // Particle model used to define the Primary.
    const Sweep::ParticleModel *m_pmodel;

    // Unique properties.
    fvector m_comp;   // Primary composition.
    fvector m_values; // Other primary values (defined at run time).
    double m_createt;   // Time at which primary was created.
    double m_time;      // Last time primary was updated.  Required for LPDA.

    //! Basic derived properties (calculated from above properties).
    double m_diam; // Equivalent spherical diameter.
    double m_dcol; // Collision diameter.
    double m_dmob; // Mobility diameter.
    double m_surf; // Surface area.
    double m_vol;  // Volume.
    double m_mass; // Mass.
    int m_numcarbon; //!< Number of carbon atoms.
    int m_frag;      //!< Fragmentation flag.
	int m_numOf6Rings;

	// Property for titania phase transformation model
	double m_phaseterm;

    // Primary class cannot be created without knowledge of the
    // particle model, therefore default constructor is protected.
    Primary(void);


    // DERIVED PROPERTIES.

    // Sets the particle cache to that of a spherical particle.
    void setSphereCache(void);

    //! Calculate the number of allowable adjustments for a LPDA process
    unsigned int CalculateMaxAdjustments(const fvector &dcomp, unsigned int n) const;

	//! Calculate the maxmimum number of allowable adjustments for a process
	unsigned int CalculateMaxAdjustments(const fvector &dcomp) const;

    // MEMORY MANAGEMENT.

    // Release all memory associated with object.
    void releaseMem(void);

    // Initialisation routine.
    void init(void);


};
} //namespace AggModels
} //namespace Sweep

#endif
