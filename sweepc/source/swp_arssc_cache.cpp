/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ARSSC_Cache class declared in the
    swp_arssc_cache.h header file.

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

#include "swp_arssc_cache.h"
#include <stdexcept>

using namespace std;
using namespace Sweep;
using namespace Sweep::SubModels;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
SubModels::ARSSC_Cache::ARSSC_Cache(void)
: SubModelCache(), m_sites(ARSSC_Model::SiteTypeCount, 0.0)
{
}

// Default constructor (public).
ARSSC_Cache::ARSSC_Cache(Sweep::ParticleCache &parent)
: SubModelCache(parent), m_sites(ARSSC_Model::SiteTypeCount, 0.0)
{
}

// Copy constructor.
ARSSC_Cache::ARSSC_Cache(const ARSSC_Cache &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Default destructor.
ARSSC_Cache::~ARSSC_Cache()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator (ARSSC_Cache RHS).
ARSSC_Cache &ARSSC_Cache::operator=(const ARSSC_Cache &rhs)
{
    if (this != &rhs) {
        m_sites.assign(rhs.m_sites.begin(), rhs.m_sites.end());
    }
    return *this;
}

// Assignment operator (ARSSC_Model RHS).
ARSSC_Cache &ARSSC_Cache::operator=(const ARSSC_Model &rhs)
{
    for (unsigned int i=0; i!=ARSSC_Model::SiteTypeCount; ++i) {
        m_sites[i] = rhs.SiteCount((ARSSC_Model::SiteType)i);
    }
    return *this;
}

// Assignment operator (SubModelCache RHS).
ARSSC_Cache &ARSSC_Cache::operator=(const SubModelCache &rhs)
{
    // Attempt to cast the RHS as a ARSSC_Cache.  This will throw
    // an exception if it isn't possible to cast.
    return operator=(dynamic_cast<const ARSSC_Cache&>(rhs));
}

// Assignment operator (SubModel RHS).
ARSSC_Cache &ARSSC_Cache::operator=(const SubModel &rhs)
{
    // Attempt to cast the RHS as a ARSSC_Model.  This will throw
    // an exception if it isn't possible to cast.
    return operator=(dynamic_cast<const ARSSC_Model&>(rhs));
}

// Compound assignment operator (ARSSC_Cache RHS).
ARSSC_Cache &ARSSC_Cache::operator+=(const ARSSC_Cache &rhs)
{
    for (unsigned int i=0; i!=min(m_sites.size(), rhs.m_sites.size()); ++i) {
        m_sites[i] += rhs.m_sites[i];
    }
    return *this;
}

// Compound assignment operator (ARSSC_Model RHS).
ARSSC_Cache &ARSSC_Cache::operator+=(const ARSSC_Model &rhs)
{
    for (unsigned int i=0; i!=ARSSC_Model::SiteTypeCount; ++i) {
        m_sites[i] += rhs.SiteCount((ARSSC_Model::SiteType)i);
    }
    return *this;
}

// Compound assignment (SubModelCache RHS).
ARSSC_Cache &ARSSC_Cache::operator+=(const SubModelCache &rhs)
{
    // Attempt to cast the RHS as a ARSSC_Cache.  This will throw
    // an exception if it isn't possible to cast.
    return operator+=(dynamic_cast<const ARSSC_Cache&>(rhs));
}

// Compound assignment (SubModel RHS).
ARSSC_Cache &ARSSC_Cache::operator+=(const SubModel &rhs)
{
    // Attempt to cast the RHS as a ARSSC_Model.  This will throw
    // an exception if it isn't possible to cast.
    return operator+=(dynamic_cast<const ARSSC_Model&>(rhs));
}


// Resets the model data to the default state.
void ARSSC_Cache::Clear()
{
    m_sites.assign(ARSSC_Model::SiteTypeCount, 0.0);
}


// PROPERTIES.

// Returns the property with the given ID.
real ARSSC_Cache::Property(unsigned int id) const
{
    if (id < ARSSC_Model::SiteTypeCount) {
        return m_sites[id];
    } else {
        return 0.0;
    }
}


// READ/WRITE/COPY.

// Returns a copy of the model data.
ARSSC_Cache *const ARSSC_Cache::Clone(void) const
{
    return new ARSSC_Cache(*this);
}

// Returns the model data type.  Used to identify different models
// and for serialisation.
SubModelType ARSSC_Cache::ID(void) const
{
    return SubModels::ARSSC_Model_ID;
}

// Writes the object to a binary stream.
void ARSSC_Cache::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write number of sites.
        unsigned int n = (unsigned int)m_sites.size();
        out.write((char*)&n, sizeof(n));

        // Write site counts.
        double val = 0.0;
        for (unsigned int i=0; i!=n; ++i) {
            val = (double)m_sites[i];
            out.write((char*)&val, sizeof(val));
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ARSSC_Cache::Serialize).");
    }
}

// Reads the object from a binary stream.
void ARSSC_Cache::Deserialize(std::istream &in, ParticleCache &parent)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0, i=0;
        double       val = 0.0;

        switch (version) {
            case 0:
                // Read number of sites.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Resize sites vector.
                m_sites.resize(n, 0.0);

                // Read sites.
                for (i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_sites[i] = (real)val;
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ARSSC_Cache::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ARSSC_Cache::Deserialize).");
    }
}
