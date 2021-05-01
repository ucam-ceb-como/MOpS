/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Particle class represents the top node in the sub-particle tree.  This
    object is placed in the Ensemble and exhibits an interface which allows
    bulk particle properties to be summed and stored in the ensemble binary tree.

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

#ifndef SWEEP_PARTICLE_H
#define SWEEP_PARTICLE_H

#include "swp_params.h"
#include "swp_particle_model.h"
#include "swp_property_indices.h"
#include "swp_model_factory.h"

#include "camxml.h"

#include <vector>
#include <list>
#include <iostream>

namespace Sweep
{

namespace AggModels {
    // Forward declaration
    class Primary;
}

/*!
 * \brief Particle which can move around, coagulate and have internal structure
 * 
 * The interface that this class should offer is not entirely clear.  This class
 * may eventually become a template parameter so that simulations can be compiled
 * for different particle types.  Alternatively an interface may be defined and then
 * subclassed.  This will require some thought. 
 */
class Particle
{
public:
	// Constructors.
    Particle(                             // Initialising constructor.
        double time,                        // Create time.
        const Sweep::ParticleModel &model // Defining particle model.
        );

    //! Create a particle with a specified statistical weight
    Particle(
        double time,
        double weight,
        const Sweep::ParticleModel &model
        );

    Particle(Sweep::AggModels::Primary &pri);        // Initialising constructor (from primary).
    Particle(const Particle &copy);       // Copy constructor.
    
    //! Stream-reading constructor
    Particle(
        std::istream &in,                     // Input stream       
        const Sweep::ParticleModel &model,    // Model to which this particle subscribes
        void *duplicates                      // Information on duplicated PAH
        );

	// Destructor.
    virtual ~Particle(void);
    
    //! Create a new particle using the model according to the xml data
    static Particle* createFromXMLNode(const CamXML::Element& xml,
                                       const Sweep::ParticleModel& model);

    //! Create a new polydisperse particle using the model according to the xml data
    static Particle* createFromXMLNodeDetailed(const CamXML::Element& xml,
		const Sweep::ParticleModel& model, rng_type &rng);

    // Operators.
    Particle &operator=(const Particle &rhs);


    // PRIMARY PARTICLE CHILD.

    //! Pointer to the child primary particle
    Sweep::AggModels::Primary *const Primary();
    //! Pointer to the child primary particle
    const Sweep::AggModels::Primary *const Primary() const;

    // BASIC PROPERTIES.

    //! Returns the particle equivalent sphere diameter.
    double SphDiameter(void) const;

    //! Returns the collision diameter.
    double CollDiameter(void) const;

    //! Returns the mobility diameter.
    double MobDiameter(void) const;

    //! Returns the surface area.
    double SurfaceArea(void) const;

    //! Returns the equivalent sphere surface area, based on the volume.
    double SphSurfaceArea(void) const;

    //! Returns the volume.
    double Volume(void) const;

    //! Returns the mass.
    double Mass(void) const;

    //! Returns the number of carbons.
    int NumCarbon(void) const;

	//! Returns the number of rings.
	int NumRings(void) const;

    //! Returns fragmentation flag.
    int Frag(void) const;

    //! Returns the number of surface reaction sites.
    double GetSites(void) const;

    //! Returns the sintering rate (silica).
    double GetSintRate(void) const;

    //! Returns the coverage fraction
    double GetCoverageFraction(void) const;

    //! Geometric average diameter of aggregate sub-units
    double avgeomdiam(double) const;

	//! Phase composition term for phase transformation process
    double GetPhaseTerm(void) const;


    //! Returns the property with the given ID.
    double Property(Sweep::PropID id) const;


    // COMPOSITION.

    //! Returns the composition vector.
    const fvector &Composition(void) const;

    //! Returns the ith component value.  Returns 0.0 if i invalid.
    double Composition(unsigned int i) const;


    // TRACKER VALUES.

    //! Returns the tracked values vector.
    const fvector &Values(void) const;

    //! Returns the ith tracked value.  Returns 0.0 if i invalid.
    double Values(unsigned int i) const;

	// PHASE

	double PhaseMass(int i) const;

    // POSITION DATA
    
    //! Get spatial position
    double getPosition() const {return m_Position;}
    
    //! Get time at which spatial position was valid
    double getPositionTime() const {return m_PositionTime;}
    
    //! Set spatial position of particle and time at which it applied
    void setPositionAndTime(const double x, const double t);

    // STATISTICAL WEIGHT
    //! Get statistical weight
    double getStatisticalWeight() const {return m_StatWeight;}

    //! Set statistical weight
    void setStatisticalWeight(double wt) {m_StatWeight = wt;}

    // CREATE TIME
    //! Time particle was created.
    double CreateTime(void) const {return m_createt;}

    // LAST UPDATE TIME FOR LPDA
    //! Time at which LPDA was last performed on this particle.
    void SetTime(double t);

    //! Time at which LPDA was last performed on this particle.
    double LastUpdateTime(void) const {return mLPDAtime;}

    // COAG COUNT
    //! Number of coagulations since count was reset
    unsigned int getCoagCount() const;

    //! Increment count of coagulation events
    void incrementCoagCount() {++m_CoagCount;}

    //! Reset count of coagulation events
    void resetCoagCount() {m_CoagCount=0;}

    // COAG COUNT
    //! Number of coagulations since count was reset
    unsigned int getFragCount() const;

    //! Increment count of coagulation events
    void incrementFragCount() {++m_FragCount;}

    //! Reset count of coagulation events
    void resetFragCount() {m_FragCount=0;}

    // PARTICLE OPERATIONS.

    //! Adjust particle composition as a result of surface reactions and other processes
    unsigned int Adjust(
        const fvector &dcomp,             // Composition changes.
        const fvector &dvalues,           // Tracker variable changes.
        rng_type &rng,                    // Random number for leaf node
        unsigned int n                    // Number of times to perform adjustment.
        );

    //! Adjust particle composition as a result of IntParticle processes
    unsigned int AdjustIntPar(
        const fvector &dcomp,             // Composition changes.
        const fvector &dvalues,           // Tracker variable changes.
        rng_type &rng,                    // Random number for leaf node
        unsigned int n                    // Number of times to perform adjustment.
        );

	//! Adjust particle composition as a result of a phase transformation process
	unsigned int AdjustPhase(const fvector &dcomp,
                              const fvector &dvalues,
                              rng_type &rng,
                              unsigned int n);
	
	//! Melting point dependent phase change
	void Melt(rng_type &rng, Cell &sys);

    //! Combines this particle with another.
    Particle &Coagulate(const Particle &sp, rng_type &rng);

    //! Combines this particle with another.
    Particle &Fragment(const Particle &sp, rng_type &rng);

    //! Sinter over a given time step
    void Sinter(
        double dt,         // Delta-t for sintering.
        Cell &sys, // System which defines particle's environment.
        const Processes::SinteringModel &model, // Sintering model to use.
        rng_type &rng,   // Random number generator
        double wt     // Statistical weight
        );


    // Recalculate derived properties from the primary particle
    void UpdateCache();

	// PARTICLE TRACKING FOR VIDEOS

	//! Returns the frame position and orientation, and primary coordinates
	void getFrameCoords(std::vector<fvector> &coords) const;

	//! Initialise primary particle tracking for videos
	void setTracking();

	//! Remove primary tracking
	void removeTracking();

    // READ/WRITE/COPY.

    //! Clone the particle.
    Particle *const Clone() const;
    
    //! Internal consistency check
    bool IsValid() const;

    //! Write a POVRAY image of the particle
    void writeParticlePOVRAY(std::ofstream &out) const;

    //! Write the object to a binary stream.
	virtual void Serialize(std::ostream &out, void *duplicates) const;

    // Deserialisation is handled through the stream reading constructor

private:

    //! Spatial position of particle
    double m_Position;

    //! Time at which position was valid
    double m_PositionTime;

    //! Statistical weight of particle units: \f$ m^{-3}\f$
    double m_StatWeight;

    //! Primary particle containing physical details of this particle
    Sweep::AggModels::Primary *m_primary;

    //! Number of coagulations experienced by this particle
    unsigned int m_CoagCount;

    //! Number of coagulations experienced by this particle
    unsigned int m_FragCount;

    //! Time at which particle was created (earliest part).
    double m_createt;

    //! Last time particle was updated.  Required for LPDA.
    double mLPDAtime;

    // Can't create a particle without knowledge of the components
    // and the tracker variables.
    Particle(void);
};

typedef std::vector<Particle*> PartPtrVector;
typedef std::list<Particle*> PartPtrList;
}

#endif
