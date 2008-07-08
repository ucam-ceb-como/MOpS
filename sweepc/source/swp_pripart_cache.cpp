/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the PriPartCache class declared in the
    swp_pripart_cache.h header file.

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

#include "swp_pripart_cache.h"
#include "swp_pripart_primary.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
PriPartCache::PriPartCache(void)
: m_npri(0), m_prisurf(0.0), m_avgdiam(0.0)
{
}

// Default constructor (public).
PriPartCache::PriPartCache(ParticleCache &parent)
: SurfVolCache(parent), m_npri(0), m_prisurf(0.0), m_avgdiam(0.0)
{
}

// Copy constructor.
PriPartCache::PriPartCache(const PriPartCache &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
PriPartCache::PriPartCache(std::istream &in, ParticleCache &parent)
{
    Deserialize(in, parent);
}

// Default destructor.
PriPartCache::~PriPartCache()
{
    // Nothing special to destruct.
}


// ASSIGNMENT OPERATOR OVERLOADING.

// Assignment operator (PriPartCache RHS).
PriPartCache &PriPartCache::operator=(const PriPartCache &rhs)
{
    if (this != &rhs) {
        SurfVolCache::operator=(rhs);
        m_npri    = rhs.m_npri;
        m_prisurf = rhs.m_prisurf;
        m_avgdiam = rhs.m_avgdiam;
    }
    return *this;
}

// Assignment operator (PriPartPrimary RHS).
PriPartCache &PriPartCache::operator=(const PriPartPrimary &rhs)
{
    SurfVolCache::operator=(rhs);
    m_npri    = rhs.PriCount();
    m_prisurf = rhs.PriSurface();
    m_avgdiam = rhs.AvgPriDiameter();
    return *this;
}

// Assignment operator (SurfVolCache RHS).
PriPartCache &PriPartCache::operator=(const SurfVolCache &rhs)
{
    if (this != &rhs) {
        if (rhs.ID() == SurfVol_ID) {
            // The RHS actually is a SurfVolCache.
            SurfVolCache::operator=(rhs);
            m_npri    = rhs.PP_Count();
            m_prisurf = rhs.SurfaceArea();
            m_avgdiam = rhs.PP_Diameter();
        } else {
            // Attempt to cast RHS as a PriPartCache.  This will
            // throw an exception if the cast fails.
            *this = dynamic_cast<const PriPartCache&>(rhs);
        }
    }
    return *this;
}

// Assignment operator (SurfVolPrimary RHS).
PriPartCache &PriPartCache::operator=(const SurfVolPrimary &rhs)
{
    if (rhs.AggID() == SurfVol_ID) {
        // The RHS actually is a surface-volume primary.
        SurfVolCache::operator=(rhs);
        m_npri    = rhs.PP_Count();
        m_prisurf = rhs.SurfaceArea();
        m_avgdiam = rhs.PP_Diameter();
    } else {
        // Attempt to cast RHS as a PriPartPrimary.  This will
        // throw an exception if the cast fails.
        *this = static_cast<const PriPartPrimary&>(rhs);
    }
    return *this;
}

// Assignment operator (AggModelCache RHS).
PriPartCache &PriPartCache::operator=(const AggModelCache &rhs)
{
    // Attempt to cast the RHS as a PriPartCache.  This will throw
    // an exception if it isn't possible to cast.
    return operator=(dynamic_cast<const PriPartCache&>(rhs));
}

// Assignment operator (Primary RHS).
PriPartCache &PriPartCache::operator=(const Primary &rhs)
{
    // Continue based on the type of the primary.
    if (rhs.AggID() == Spherical_ID) {
        // Pass assignment down inheritance chain.
        SurfVolCache::operator=(rhs);
        m_npri    = 1;
        m_prisurf = rhs.SurfaceArea();
        m_avgdiam = rhs.SphDiameter();
    } else if (rhs.AggID() == SurfVol_ID) {
        // Pass assignment down inheritance chain.
        SurfVolCache::operator=(rhs);
        m_npri = static_cast<const SurfVolPrimary&>(rhs).PP_Count();
        m_prisurf = rhs.SurfaceArea();
        m_avgdiam = static_cast<const SurfVolPrimary&>(rhs).PP_Diameter();
    } else {
        // Attempt to cast the RHS as a PriPartPrimary.  This will throw
        // an exception if it isn't possible to cast.
        *this = dynamic_cast<const PriPartPrimary&>(rhs);
    }
    return *this;
}


// COMPOUND ASSIGNMENT OPERATOR OVERLOADING.

// Compound assignment (PriPartCache RHS).
PriPartCache &PriPartCache::operator+=(const PriPartCache &rhs)
{
    SurfVolCache::operator+=(rhs);
    m_npri    += rhs.m_npri;
    m_prisurf += rhs.m_prisurf;
    m_avgdiam += rhs.m_avgdiam;
    return *this;
}

// Compound assignment (PriPartPrimary RHS).
PriPartCache &PriPartCache::operator+=(const PriPartPrimary &rhs)
{
    SurfVolCache::operator+=(rhs);
    m_npri    += rhs.PriCount();
    m_prisurf += rhs.PriSurface();
    m_avgdiam += rhs.AvgPriDiameter();
    return *this;
}

// Compound assignment (SurfVolCache RHS).
PriPartCache &PriPartCache::operator+=(const SurfVolCache &rhs)
{
    if (rhs.ID() == SurfVol_ID) {
        // Pass operation down inheritance chain.
        SurfVolCache::operator+=(rhs);
        m_npri    += rhs.PP_Count();
        m_prisurf += rhs.SurfaceArea();
        m_avgdiam += rhs.PP_Diameter();
    } else {
        *this += dynamic_cast<const PriPartCache&>(rhs);
    }
    return *this;
}

// Compound assignment (SurfVolPrimary RHS).
PriPartCache &PriPartCache::operator+=(const SurfVolPrimary &rhs)
{
    if (rhs.AggID() == SurfVol_ID) {
        // Pass operation down inheritance chain.
        SurfVolCache::operator+=(rhs);
        m_npri    += rhs.PP_Count();
        m_prisurf += rhs.SurfaceArea();
        m_avgdiam += rhs.PP_Diameter();
    } else {
        *this += dynamic_cast<const PriPartPrimary&>(rhs);
    }
    return *this;
}

// Compound assignment (AggModelCache RHS).
PriPartCache &PriPartCache::operator+=(const AggModelCache &rhs)
{
    // Attempt to cast the RHS as a PriPartCache.  This will throw
    // an exception if it isn't possible to cast.
    return operator+=(dynamic_cast<const PriPartCache&>(rhs));
}

// Compound assignment (Primary RHS).
PriPartCache &PriPartCache::operator+=(const Primary &rhs)
{
    // Continue based on the type of the primary.
    if (rhs.AggID() == Spherical_ID) {
        // Pass assignment down inheritance chain.
        SurfVolCache::operator+=(rhs);
        m_npri    += 1;
        m_prisurf += rhs.SurfaceArea();
        m_avgdiam += rhs.SphDiameter();
    } else if (rhs.AggID() == SurfVol_ID) {
        // Pass assignment down inheritance chain.
        SurfVolCache::operator+=(rhs);
        m_npri    += static_cast<const SurfVolPrimary&>(rhs).PP_Count();
        m_prisurf += rhs.SurfaceArea();
        m_avgdiam    += static_cast<const SurfVolPrimary&>(rhs).PP_Diameter();
    } else {
        // Attempt to cast the RHS as a PriPartPrimary.  This will throw
        // an exception if it isn't possible to cast.
        *this += dynamic_cast<const PriPartPrimary&>(rhs);
    }
    return *this;
}


// DATA MANAGEMENT.

// Resets the model data to the default state.
void PriPartCache::Clear()
{
    m_npri    = 0;
    m_prisurf = 0.0;
    m_avgdiam = 0.0;
}


// PROPERTIES.

// Returns the number of primary particles.
unsigned int PriPartCache::Count(void) const {return m_npri;}

// Sets the number of primary particles.
void PriPartCache::SetCount(unsigned int n) {m_npri = n;}

// Returns the primary particle surface area.
real PriPartCache::PriSurfaceArea(void) const {return m_prisurf;}

// Sets the primary particle surface area.
void PriPartCache::SetPriSurfaceArea(real s) {m_prisurf = s;}

// Returns the average primary-particle diameter.
real PriPartCache::AvgPriDiameter(void) const {return m_avgdiam;}

// Sets the average primary-particle diameter.
void PriPartCache::SetAvgPriDiameter(real d) {m_avgdiam = d;}


// READ/WRITE/COPY.

// Returns a copy of the model data.
PriPartCache *const PriPartCache::Clone(void) const
{
    return new PriPartCache(*this);
}

// Returns the model data type.  Used to identify different models
// and for serialisation.
AggModelType PriPartCache::ID(void) const {return PriPartList_ID;}

// Writes the object to a binary stream.
void PriPartCache::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output parent class.
        SurfVolCache::Serialize(out);

        // Output primary count.
        unsigned int n = (unsigned int)m_npri;
        out.write((char*)&n, sizeof(n));

        // Output primary surface.
        double s = (double)m_prisurf;
        out.write((char*)&s, sizeof(s));

        // Output average primary diameter.
        double d = (double)m_avgdiam;
        out.write((char*)&d, sizeof(d));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PriPartCache::Serialize).");
    }
}

// Reads the object from a binary stream.
void PriPartCache::Deserialize(std::istream &in, ParticleCache &parent)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        double surf    = 0.0;
        double diam    = 0.0;

        switch (version) {
            case 0:
                // Read parent class.
                SurfVolCache::Deserialize(in, parent);

                // Read the number of primaries.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_npri = n;

                // Read the primary surface area.
                in.read(reinterpret_cast<char*>(&surf), sizeof(surf));
                m_prisurf = (real)surf;

                // Read the average primary diameter.
                in.read(reinterpret_cast<char*>(&diam), sizeof(diam));
                m_avgdiam = (real)diam;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, PriPartCache::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PriPartCache::Deserialize).");
    }
}
