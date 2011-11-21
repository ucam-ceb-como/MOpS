/*
  Author(s):      Markus Sander (ms785)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Markus Sander.
  Extended by Shraddha Shekar (ss663) and William Menz (wjm34)

  File purpose:
	Implementation of SilicaCache class defined in swp_silica_cache.h


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

#include "swp_silica_cache.h"
#include "swp_aggmodel_type.h"
#include "swp_aggmodel_cache.h"
#include "swp_silica_primary.h"
#include <stdexcept>
#include <iostream>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS

// Default constructor (private).
SilicaCache::SilicaCache(void)
{
	m_numSi=2;
	m_numO=1;
	m_numOH=6;
	m_silicaDiameter=0.0;
    m_numprimary=0;
    m_sqrtLW=0.0;
    m_LdivW=0.0;
    m_primarydiam=0.0;
    m_fdim=0.0;
    m_Rg=0.0;
    m_avg_sinter=0.0;
    m_sintrate=0.0;
}

// Copy constructor.
SilicaCache::SilicaCache(const SilicaCache &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
SilicaCache::SilicaCache(std::istream &in)
{
    Deserialize(in);
}

// Default destructor.
SilicaCache::~SilicaCache()
{
    // Nothing special to destruct.
}

// ASSIGNMENT OPERATOR OVERLOADING.

// Assignment operator (silicaCache RHS).
SilicaCache &SilicaCache::operator=(const SilicaCache &rhs)
{
    if (this != &rhs) {
        m_numSi = rhs.m_numSi;
        m_numO = rhs.m_numO;
		m_numOH = rhs.m_numOH;
		m_silicaDiameter = rhs.m_silicaDiameter;
		m_numprimary=rhs.m_numprimary;
        m_sqrtLW=rhs.m_sqrtLW;
        m_LdivW=rhs.m_LdivW;
        m_primarydiam=rhs.m_primarydiam;
        m_fdim=rhs.m_fdim;
        m_Rg=rhs.m_Rg;
        m_avg_sinter=rhs.m_avg_sinter;
        m_sintrate=rhs.m_sintrate;
    }
    return *this;
}

// Assignment operator (silicaPrimary RHS).
SilicaCache &SilicaCache::operator=(const SilicaPrimary &rhs)
{
	m_numSi = rhs.NumSi();
    m_numO = rhs.NumO();
    m_numOH = rhs.NumOH();
	//m_silicaDiameter = rhs.silicaCollDiameter();
    m_numprimary=rhs.Numprimary();
    m_sqrtLW=rhs.sqrtLW();
    m_LdivW=rhs.LdivW();
    m_primarydiam=rhs.PrimaryDiam();
    m_fdim=rhs.Fdim();
    m_Rg=rhs.Rg();
    m_avg_sinter=rhs.AvgSinter();
    m_sintrate=rhs.GetSintRate();
    return *this;
}

// Assignment operator (AggModelCache RHS).
SilicaCache &SilicaCache::operator=(const AggModelCache &rhs)
{
    // Attempt to cast the RHS as a silicaCache.  This will throw
    // an exception if it isn't possible to cast.
    return operator=(dynamic_cast<const SilicaCache&>(rhs));
}

// Assignment operator (Primary RHS).
SilicaCache &SilicaCache::operator=(const Primary &rhs)
{
    if (rhs.AggID() == Spherical_ID) {
        // If the RHS is actually a spherical primary then copy
        // the required data.


    } else {
        // Attempt to cast the RHS as a silicaPrimary.  This will throw
        // an exception if it isn't possible to cast.
        operator=(dynamic_cast<const SilicaPrimary&>(rhs));
    }
    return *this;
}


// COMPOUND ASSIGNMENT OPERATOR OVERLOADING.

// Compound assignment (silicaCache RHS).
SilicaCache &SilicaCache::operator+=(const SilicaCache &rhs)
{
    m_numSi += rhs.m_numSi;
    m_numO += rhs.m_numO;
    m_numOH += rhs.m_numOH;
	m_silicaDiameter += rhs.m_silicaDiameter;
    m_numprimary+=rhs.m_numprimary;
    m_sqrtLW+=rhs.m_sqrtLW;
    m_LdivW+=rhs.m_LdivW;
    m_primarydiam+=rhs.m_primarydiam;
    m_fdim+=rhs.m_fdim;
    m_Rg+=rhs.m_Rg;
    m_avg_sinter+=rhs.m_avg_sinter;
    m_sintrate+=rhs.m_sintrate;
    return *this;
}

// Compound assignment (silicaPrimary RHS).
SilicaCache &SilicaCache::operator+=(const SilicaPrimary &rhs)
{
	m_numSi += rhs.NumSi();
    m_numO += rhs.NumO();
    m_numOH += rhs.NumOH();
	//m_silicaDiameter += rhs.silicaCollDiameter();
    m_numprimary+=rhs.Numprimary();
    m_sqrtLW+=rhs.sqrtLW();
    m_LdivW+=rhs.LdivW();
    m_primarydiam+=rhs.PrimaryDiam();
    m_fdim+=rhs.Fdim();
    m_Rg+=rhs.Rg();
    m_avg_sinter+=rhs.AvgSinter();
    m_sintrate+=rhs.GetSintRate();

    return *this;
}

// Compound assignment (AggModelCache RHS).
SilicaCache &SilicaCache::operator+=(const AggModelCache &rhs)
{
    // Attempt to cast the RHS as a silicaCache.  This will throw
    // an exception if it isn't possible to cast.
    return operator+=(dynamic_cast<const SilicaCache&>(rhs));
}

// Compound assignment (Primary RHS).
SilicaCache &SilicaCache::operator+=(const Primary &rhs)
{
    if (rhs.AggID() == Spherical_ID) {

    } else {
        // Attempt to cast the RHS as a silicaPrimary.  This will throw
        // an exception if it isn't possible to cast.
        operator+=(dynamic_cast<const SilicaPrimary&>(rhs));
    }
    return *this;
}

int SilicaCache::NumSi() const
{
    return m_numSi;
}

int SilicaCache::NumO() const
{
    return m_numO;
}

int SilicaCache::NumOH() const
{
    return m_numOH;
}

// DATA MANAGEMENT.

// Resets the model data to the default state.
void SilicaCache::Clear()
{
    m_numSi = 0;
	m_numO = 0;
	m_numOH = 0;
	m_silicaDiameter = 0.0;
    m_numprimary=0;
    m_sqrtLW=0.0;
    m_LdivW=0.0;
    m_primarydiam=0.0;
    m_fdim=0.0;
    m_Rg=0.0;
    m_avg_sinter=0.0;
    m_sintrate=0.0;
}


// AGGREGATION MODEL PROPERTIES.


// READ/WRITE/COPY.

// Returns a copy of the data.
SilicaCache *const SilicaCache::Clone(void) const {return new SilicaCache(*this);}

// Returns the model ID.
AggModelType SilicaCache::ID(void) const {return Silica_ID;}

// Writes the object to a binary stream.
void SilicaCache::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        double v = (int)m_numSi;
        out.write((char*)&v, sizeof(v));

		v = (int)m_numO;
        out.write((char*)&v, sizeof(v));

		v = (int)m_numOH;
        out.write((char*)&v, sizeof(v));

        v = (double)m_silicaDiameter;
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

        v = (double)m_avg_sinter;
        out.write((char*)&v, sizeof(v));

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, SilicaCache::Serialize).");
    }
}

// Reads the object from a binary stream.
void SilicaCache::Deserialize(std::istream &in)
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
               //read the number of silicas
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_numSi = (int)val;

				in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_numO = (int)val;

				in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_numOH = (int)val;

                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_silicaDiameter = (real)val;

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
                m_avg_sinter = (real)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, SilicaCache::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SilicaCache::Deserialize).");
    }
}
