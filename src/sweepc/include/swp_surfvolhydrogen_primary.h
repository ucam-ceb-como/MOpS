/*
  Author(s):      Robert I A Patterson
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2012 Robert I A Patterson.

  File purpose:
    The SurfVolHydrogenPrimary class implements the Primary particle class to include
    the surface-volume model and a count of the hydrogenised sites on the surface of
    soot particles.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
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

#ifndef SWEEP_SURFVOLHYDROGEN_PRIMARY_H
#define SWEEP_SURFVOLHYDROGEN_PRIMARY_H

#include "swp_params.h"
#include "swp_primary.h"
#include "swp_particle_model.h"
#include "swp_sintering_model.h"
#include "swp_aggmodel_type.h"

#include <iostream>

namespace Sweep
{
namespace AggModels 
{
class SurfVolHydrogenPrimary : public Primary
{
public:
    // Constructors.
    //   Note:  Default constructor is protected to prevent a
    //          SurfVolHydrogenPrimary being created without knowledge of the
    //          defining particle model.
    SurfVolHydrogenPrimary(                       // Initialising constructor.
        double time,                        //  - Create time.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );
    SurfVolHydrogenPrimary(const SurfVolHydrogenPrimary &copy); // Copy constructor.
    SurfVolHydrogenPrimary(                       // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );

    // Destructors.
    virtual ~SurfVolHydrogenPrimary(void);

    // Operators.
    virtual SurfVolHydrogenPrimary &operator=(const Primary &rhs);
    virtual SurfVolHydrogenPrimary &operator=(const SurfVolHydrogenPrimary &rhs);

    // AGGREGATION MODEL.

    //! Returns the aggregation model which this primary describes.
    virtual AggModels::AggModelType AggID(void) const;


    // BASIC DERIVED PROPERTIES.

    // Calculates the derived properties from the unique properties.
    virtual void UpdateCache(void);

    // Returns the equivalent spherical particle surface area.
    double SphSurfaceArea(void) const;

    // Returns the number of primary particles if the aggregate is assumed
    // to consist of mono-sized primaries.
    unsigned int PP_Count(void) const;

    // Returns the primary particle diameter if the aggregate is assumed
    // to consist of mono-sized primaries.
    double PP_Diameter(void) const;


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

    // Combines this primary with another.
    virtual SurfVolHydrogenPrimary &Coagulate(const Primary &rhs,
                                      rng_type &rng);

    // Combines this primary with another.
    virtual SurfVolHydrogenPrimary &Fragment(const Primary &rhs,
                                      rng_type &rng);

    // This routine sinters the Primary for the given length of
    // time using the provided sintering model.
    virtual void Sinter(
        double dt, // Delta-t for sintering to occur.
        Cell &sys, // System which defines primary's environment.
        const Processes::SinteringModel &model, // Sintering model to use.
        rng_type &rng,   // Random number generator
        double wt     // Statistical weight
        );


    // READ/WRITE/COPY.

    //! Returns a pointer to a new copy of the primary.
    virtual SurfVolHydrogenPrimary *const Clone(void) const;

    //! Number of active sites (only implemented for some particle models).
    virtual double GetSites() const;

    //! Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    //! Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Defining particle model.
        );

protected:

    // Primary class cannot be created without knowledge of the
    // particle model, therefore default constructor is protected.
    SurfVolHydrogenPrimary(void);

private:
    //! Index of H in the gas phase
    const unsigned iH; // Will need serialising if made non-const
};
};
}
#endif
