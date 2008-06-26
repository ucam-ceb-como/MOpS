/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The SurfVolPrimary class extends the Primary particle class to include
    the surface-volume model.
*/

#ifndef SWEEP_SURFVOL_PRIMARY_H
#define SWEEP_SURFVOL_PRIMARY_H

#include "swp_params.h"
#include "swp_primary.h"
#include "swp_particle_model.h"
#include "swp_sintering_model.h"
#include "swp_aggmodel_type.h"
#include "swp_surfvol_cache.h"
#include <iostream>

namespace Sweep
{
namespace AggModels 
{
class SurfVolPrimary : public Primary
{
public:
    // Constructors.
    //   Note:  Default constructor is protected to prevent a
    //          SurfVolPrimary being created without knowledge of the
    //          defining particle model.
    SurfVolPrimary(                       // Initialising constructor.
        real time,                        //  - Create time.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );
    SurfVolPrimary(const SurfVolPrimary &copy); // Copy constructor.
    SurfVolPrimary(                       // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );

    // Destructors.
    virtual ~SurfVolPrimary(void);

    // Operators.
    virtual SurfVolPrimary &operator=(const Primary &rhs);
    virtual SurfVolPrimary &operator=(const SurfVolPrimary &rhs);
    virtual SurfVolPrimary &operator+=(const Primary &rhs);


    // AGGREGATION MODEL.

    // Returns the aggregation model which this primary describes.
    virtual AggModels::AggModelType AggID(void) const;

    // Creates an aggregation data cache for this primary type.
    virtual AggModels::SurfVolCache *const CreateAggCache(ParticleCache &pcache) const;


    // BASIC DERIVED PROPERTIES.

    // Calculates the derived properties from the unique properties.
    virtual void UpdateCache(void);

    // Returns the equivalent spherical particle surface area.
    real SphSurfaceArea(void) const;

    // Returns the number of primary particles if the aggregate is assumed
    // to consist of mono-sized primaries.
    unsigned int PP_Count(void) const;

    // Returns the primary particle diameter if the aggregate is assumed
    // to consist of mono-sized primaries.
    real PP_Diameter(void) const;


    // OPERATIONS.

    // Adjusts the primary with the given composition and 
    // tracker values changes n times.  If the particle cannot be adjust
    // n times, then this function returns the number of times
    // it was adjusted.
    virtual unsigned int Adjust(
        const fvector &dcomp,   // Composition changes.
        const fvector &dvalues, // Tracker variable changes.
        unsigned int n=1        // Number of times to perform adjustment.
        );

    // Combines this primary with another.  This is also the
    // implementation of the + and += operators.
    virtual SurfVolPrimary &Coagulate(const Primary &rhs);

    
    // This routine sinters the Primary for the given length of
    // time using the provided sintering model.
    virtual void Sinter(
        real dt, // Delta-t for sintering to occur.
        const Cell &sys, // System which defines primary's environment.
        const Processes::SinteringModel &model // Sintering model to use.
        );


    // READ/WRITE/COPY.

    // Returns a copy of the primary.
    virtual SurfVolPrimary *const Clone(void) const;

    // Returns this object's instance.  This may seem rather circular, but
    // it has an important purpose for getting the correct object type reference
    // from a base class reference/pointer.
    virtual SurfVolPrimary &Instance();
    virtual const SurfVolPrimary &Instance() const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );

    // Returns a correctly typed reference to the Primary object.
//    virtual const SurfVolPrimary &TypedRef(void) const;

/*
    // Copies this surface-volume primary into a spherical
    // primary.  This implements the assignment operator for
    // the situation Primary = SurfVolPrimary.
    virtual Primary &CopyTo(Primary &lhs) const;

    // Copies this surface-volume primary into another surface-volume
    // primary.  This implements the assignment operator.
    virtual SurfVolPrimary &CopyTo(SurfVolPrimary &lhs) const;

    // Adds this surface-volume primary to a spherical 
    // primary.  This implements coagulation up inheritance chain.
    virtual Primary &AddTo(Primary &lhs) const;

    // Adds this surface-volume primary to another surface-volume 
    // primary.  This implements coagulation.
    virtual SurfVolPrimary &AddTo(SurfVolPrimary &lhs) const;
*/

protected:
    // The equivalent spherical-particle surface area.
    real m_sphsurf;

    // Primary class cannot be created without knowledge of the
    // particle model, therefore default constructor is protected.
    SurfVolPrimary(void);
};
};
};
#endif
