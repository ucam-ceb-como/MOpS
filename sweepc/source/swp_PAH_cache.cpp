/*
  Author(s):      Markus Sander (ms785)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Markus Sander.

  File purpose:


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

#include "swp_PAH_cache.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_particle_cache.h"
#include "swp_PAH_primary.h"
#include <stdexcept>
#include <iostream>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
PAHCache::PAHCache(void)
{
	m_numPAH=1;
	m_PAHDiameter=0.;
    m_numcarbon=0;
    m_numprimary=0;
    m_sqrtLW=0.0;
    m_LdivW=0.0;
    m_primarydiam=0.0;
    m_fdim=0.0;
    m_Rg=0.0;
    m_avg_coalesc=0;
}

// Copy constructor.
PAHCache::PAHCache(const PAHCache &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
PAHCache::PAHCache(std::istream &in)
{
    Deserialize(in);
}

// Default destructor.
PAHCache::~PAHCache()
{
    // Nothing special to destruct.
}


// ASSIGNMENT OPERATOR OVERLOADING.

// Assignment operator (PAHCache RHS).
PAHCache &PAHCache::operator=(const PAHCache &rhs)
{
    if (this != &rhs) {
        m_numPAH = rhs.m_numPAH;
		m_PAHDiameter = rhs.m_PAHDiameter;
		m_numcarbon=rhs.m_numcarbon;
        m_numprimary=rhs.m_numprimary;
        m_sqrtLW=rhs.m_sqrtLW;
        m_LdivW=rhs.m_LdivW;
        m_primarydiam=rhs.m_primarydiam;
        m_fdim=rhs.m_fdim;
        m_Rg=rhs.m_Rg;
        m_avg_coalesc=rhs.m_avg_coalesc;
    }
    return *this;
}

// Assignment operator (PAHPrimary RHS).
PAHCache &PAHCache::operator=(const PAHPrimary &rhs)
{
    m_numPAH = rhs.NumPAH();
	m_PAHDiameter = rhs.PAHCollDiameter();
	m_numcarbon=rhs.NumCarbon();
    m_numprimary=rhs.Numprimary();
    m_sqrtLW=rhs.sqrtLW();
    m_LdivW=rhs.LdivW();
    m_primarydiam=rhs.PrimaryDiam();
    m_fdim=rhs.Fdim();
    m_Rg=rhs.Rg();
    m_avg_coalesc=rhs.AvgCoalesc();
    return *this;
}

// Assignment operator (AggModelCache RHS).
PAHCache &PAHCache::operator=(const AggModelCache &rhs)
{
    // Attempt to cast the RHS as a PAHCache.  This will throw
    // an exception if it isn't possible to cast.
    return operator=(dynamic_cast<const PAHCache&>(rhs));
}

// Assignment operator (Primary RHS).
PAHCache &PAHCache::operator=(const Primary &rhs)
{
    if (rhs.AggID() == Spherical_ID) {
        // If the RHS is actually a spherical primary then copy
        // the required data.
       

    } else {
        // Attempt to cast the RHS as a PAHPrimary.  This will throw
        // an exception if it isn't possible to cast.
        operator=(dynamic_cast<const PAHPrimary&>(rhs));
    }
    return *this;
}


// COMPOUND ASSIGNMENT OPERATOR OVERLOADING.

// Compound assignment (PAHCache RHS).
PAHCache &PAHCache::operator+=(const PAHCache &rhs)
{
    m_numPAH += rhs.m_numPAH;
	m_PAHDiameter += rhs.m_PAHDiameter;
	m_numcarbon += rhs.m_numcarbon;
    m_numprimary+=rhs.m_numprimary;
    m_sqrtLW+=rhs.m_sqrtLW;
    m_LdivW+=rhs.m_LdivW;
    m_primarydiam+=rhs.m_primarydiam;
    m_fdim+=rhs.m_fdim;
    m_Rg+=rhs.m_Rg;
    m_avg_coalesc+=rhs.m_avg_coalesc;
    return *this;
}

// Compound assignment (PAHPrimary RHS).
PAHCache &PAHCache::operator+=(const PAHPrimary &rhs)
{
    m_numPAH += rhs.NumPAH();
	m_PAHDiameter += rhs.PAHCollDiameter();
	m_numcarbon+=rhs.NumCarbon();
    m_numprimary+=rhs.Numprimary();
    m_sqrtLW+=rhs.sqrtLW();
    m_LdivW+=rhs.LdivW();
    m_primarydiam+=rhs.PrimaryDiam();
    m_fdim+=rhs.Fdim();
    m_Rg+=rhs.Rg();
    m_avg_coalesc+=rhs.AvgCoalesc();

    return *this;
}

// Compound assignment (AggModelCache RHS).
PAHCache &PAHCache::operator+=(const AggModelCache &rhs)
{
    // Attempt to cast the RHS as a PAHCache.  This will throw
    // an exception if it isn't possible to cast.
    return operator+=(dynamic_cast<const PAHCache&>(rhs));
}

// Compound assignment (Primary RHS).
PAHCache &PAHCache::operator+=(const Primary &rhs)
{
    if (rhs.AggID() == Spherical_ID) {
       
    } else {
        // Attempt to cast the RHS as a PAHPrimary.  This will throw
        // an exception if it isn't possible to cast.
        operator+=(dynamic_cast<const PAHPrimary&>(rhs));
    }
    return *this;
}


real PAHCache::NumPAH() const
{
    return m_numPAH;
}


// DATA MANAGEMENT.

// Resets the model data to the default state.
void PAHCache::Clear()
{
    m_numPAH = 0;
	m_PAHDiameter = 0.0;
	m_numcarbon =0;
    m_numprimary=0;
    m_sqrtLW=0.0;
    m_LdivW=0.0;
    m_primarydiam=0.0;
    m_fdim=0.0;
    m_Rg=0.0;
    m_avg_coalesc=0.0;
}


// AGGREGATION MODEL PROPERTIES.


// READ/WRITE/COPY.

// Returns a copy of the data.
PAHCache *const PAHCache::Clone(void) const {return new PAHCache(*this);}

// Returns the model ID.
AggModelType PAHCache::ID(void) const {return PAH_ID;}

// Writes the object to a binary stream.
void PAHCache::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        double v = (int)m_numPAH;
        out.write((char*)&v, sizeof(v));

        v = (double)m_PAHDiameter;
        out.write((char*)&v, sizeof(v));

        v = (int)m_numcarbon;
        out.write((char*)&v, sizeof(v));
		
        v = (int)m_numprimary;
        out.write((char*)&v, sizeof(v));

        v = (double)m_sqrtLW;
        out.write((char*)&v, sizeof(v));

        v = (double)m_LdivW;
        out.write((char*)&v, sizeof(v));

        v = (double)m_primarydiam;
        out.write((char*)&v, sizeof(v));

        v = (double)m_fdim;
        out.write((char*)&v, sizeof(v));

        v = (double)m_Rg;
        out.write((char*)&v, sizeof(v));

        v = (double)m_avg_coalesc;
        out.write((char*)&v, sizeof(v));


    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PAHCache::Serialize).");
    }
}

// Reads the object from a binary stream.
void PAHCache::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val     = 0.0;

        switch (version) {
            case 0:
               //read the number of PAHs
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_numPAH = (int)val;

                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_PAHDiameter = (real)val;

				in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_numcarbon = (int)val;

				in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_numprimary = (int)val;

				in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_sqrtLW = (real)val;

				in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_LdivW = (real)val;

				in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_primarydiam = (real)val;

				in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_fdim = (real)val;

				in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_Rg = (real)val;	

				in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_avg_coalesc = (real)val;
                   
                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, PAHCache::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PAHCache::Deserialize).");
    }
}
