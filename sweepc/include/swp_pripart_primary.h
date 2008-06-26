/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The PriPartPrimary class extends the SurfVolPrimary particle class to include
    the primary-particle list as described in:

       West et al., Ind. Eng. Chem. Res., 46, 6147-6156, 2007
*/

#ifndef SWEEP_PRIPART_PRIMARY_H
#define SWEEP_PRIPART_PRIMARY_H

#include "swp_params.h"
#include "swp_surfvol_primary.h"
#include "swp_particle_model.h"
#include "swp_aggmodel_type.h"
#include "swp_pripart_cache.h"
#include <iostream>
#include <vector>

namespace Sweep
{
namespace AggModels 
{
class PriPartPrimary : public SurfVolPrimary
{
public:
    // Constructors.
    //   Note:  Default constructor is protected to prevent a
    //          SurfVolPrimary being created without knowledge of the
    //          defining particle model.
    PriPartPrimary(                       // Initialising constructor.
        real time,                        //  - Create time.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );
    PriPartPrimary(const PriPartPrimary &copy); // Copy constructor.
    PriPartPrimary(                       // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );

    // Destructors.
    virtual ~PriPartPrimary(void);

    // Operators.
    virtual PriPartPrimary &operator=(const PriPartPrimary &rhs);
    virtual PriPartPrimary &operator=(const SurfVolPrimary &rhs);
    virtual PriPartPrimary &operator=(const Primary &rhs);
//    virtual PriPartPrimary &operator+=(const PriPartPrimary &rhs);
//    virtual PriPartPrimary &operator+=(const SurfVolPrimary &rhs);
    virtual PriPartPrimary &operator+=(const Primary &rhs);


    // AGGREGATION MODEL.

    // Returns the aggregation model which this primary describes.
    virtual AggModels::AggModelType AggID(void) const;

    // Creates an aggregation data cache for this primary type.
    AggModels::PriPartCache *const CreateAggCache(ParticleCache &pcache) const;


    // PRIMARY PARTICLE LIST.

    // Returns the number of primary particles in the list.
    unsigned int PriCount(void) const;

    // Returns the number of monomers in the ith primary.
    unsigned int MonomerCount(unsigned int i) const;

    // Returns the diameter of the ith primary.
    real PriDiameter(unsigned int i) const;

    // Returns the average primary particle diameter.
    real AvgPriDiameter(void) const;

    // Returns the surface area of the ith primary.
    real PriSurface(unsigned int i) const;

    // Returns the total surface area of the primaries in the list.
    real PriSurface(void) const;

    // Returns the volume of the ith primary.
    real PriVolume(unsigned int i) const;

    // Returns the mass of the ith primary.
    real PriMass(unsigned int i) const;


    // DEFINING PARTICLE COMPONENT.

    // Returns the index of the particle model component which
    // defines the primary monomers.
    static unsigned int CompID(void);

    // Sets the index of the particle model component which
    // defines the primary monomers.
    static void SetCompID(unsigned int id);


    // BASIC DERIVED PROPERTIES.

    // Calculates the derived properties from the unique properties.
    void UpdateCache(void);


    // OPERATIONS.

    // Adjusts the primary with the given composition and 
    // tracker values changes n times.  If the particle cannot be adjust
    // n times, then this function returns the number of times
    // it was adjusted.
    unsigned int Adjust(
        const fvector &dcomp,   // Composition changes.
        const fvector &dvalues, // Tracker variable changes.
        unsigned int n=1        // Number of times to perform adjustment.
        );

    // Combines this primary with another.  This is also the
    // implementation of the + and += operators.
    PriPartPrimary &Coagulate(const Primary &rhs);

    // This routine sinters the Primary for the given length of
    // time using the provided sintering model.
    virtual void Sinter(
        real dt,         // Delta-t for sintering to occur.
        const Cell &sys, // System which defines primary's environment.
        const Processes::SinteringModel &model // Sintering model to use.
        );


    // READ/WRITE/COPY.

    // Returns a copy of the primary.
    virtual PriPartPrimary *const Clone(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );


protected:
    // Structure to describe a primary particle in the list.  Primary
    // particles are described only by the number of monomers they
    // contain.  They are assumed to be spherical.  Monomers are 
    // associated with a particular component in the particle model.  The
    // index of this component is given by the static variable m_icomp.  In
    // this way all pripartlist primaries use the same component.
    struct PriPart {
        unsigned int Monomers;
        real Surface;
    };

    // Vector of primary particles (the list).
    std::vector<PriPart> m_primaries;

    // Total surface area of primaries.
    real m_totprisurf;

    // Average primary particle diameter.
    real m_avgpridiam;

    // Index of the component in the particle model which defines the
    // primary particle monomer.
    static unsigned int m_icomp;

    // Primary class cannot be created without knowledge of the
    // particle model, therefore default constructor is protected.
    PriPartPrimary(void);


    // MEMORY MANAGEMENT.

    // Release all memory associated with object.
    void releaseMem(void);

    // Initialisation routine.
    void init(void);


    // PRIPART-LIST MANAGEMENT.

    // Calculates the surface area of a spherical primary with the
    // given monomer count.
    real calcSurf(unsigned int n) const;

    // Calculates the average primary-particle diameter for
    // the current list.  The primaries are assumed to be
    // spherical.
    real calcAvgDiam();

    // Randomly creates the pripart-list to match the number of monomers
    // and as-good-as match the surface area.
//    void createPriList(
//        unsigned int n, // Number of monomers.
//        real surf       // Surface area to match.
//        );

    // Adds primaries to the list of random size (where total size
    // sums to n), until the added surface area is as close to s
    // as reasonably possible without going over.
//    void addRandom(
//        unsigned int n, // Number of monomers to add.
//        real s          // Surface area to attempt to match.
//        );

    // Distributes a number of monomers among the pri-particles, weighted
    // by their surface areas.
    void distMonomers(unsigned int n);

    // Removes n monomers from the primary particle list.  Monomers
    // are chosen by weighting the particles by their surface area, 
    // and starting selection from the largest primary.
    void removeMonomers(unsigned int n);

    // Updates the primary particle list to attempt to match the current
    // primary surface area (given by surface-volume model).
    void updatePrimaries(void);

    // Merges a pri-particle list into the current list.
    void mergeInList(const std::vector<PriPart> &list);

    // Sorts the primary particle list between the two
    // given iterators in descending order of mass.
    void sortList(
        unsigned int i1, // Index to start of sort range.
        unsigned int i2  // Index to end of sort range.
        );
};
};
};
#endif
