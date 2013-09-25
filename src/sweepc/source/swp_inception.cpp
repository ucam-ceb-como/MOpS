/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Inception class declared in the
    swp_inception.h header file.

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

#include "swp_inception.h"

#include "swp_mechanism.h"

#include "local_geometry1d.h"

#include <cmath>
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Inception::Inception(void)
: Process()
{
}

// Initialising constructor.
Inception::Inception(const Sweep::Mechanism &mech)
: Process(mech)
{
}

// Copy constructor.
Inception::Inception(const Inception &copy)
{
    *this = copy;
}

// Stream-reading constructor.
Inception::Inception(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
Inception::~Inception(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
Inception &Inception::operator =(const Inception &rhs)
{
    if (this != &rhs) {
        Process::operator =(rhs);
        m_a    = rhs.m_a;
        m_newcomp = rhs.m_newcomp;
        m_newvals = rhs.m_newvals;
    }
    return *this;
}


// PROPERTIES OF INCEPTED PARTICLES.

// Returns the composition vector of the new particle.
const fvector &Inception::ParticleComp(void) const {return m_newcomp;}

// Returns the amount of the ith component of the new particle.
double Inception::ParticleComp(unsigned int i) const
{

    if (i < m_newcomp.size()) {
        return m_newcomp[i];
    } else {
        return 0.0;
    }
}

// Sets the amount of the ith component in the new particle.
void Inception::SetParticleComp(unsigned int i, double comp)
{
    if (i < m_mech->ComponentCount()) {
        // Ensure vector is sufficiently long.
        if (m_newcomp.size() < m_mech->ComponentCount()) {
            m_newcomp.resize(m_mech->ComponentCount(),0.0);
        }
        // Set value.
        m_newcomp[i] = comp;
    }
}

// Returns the tracker variable vector of the new particle.
const fvector &Inception::ParticleTrackers(void) const
{
    return m_newvals;
}

// Returns the value of the ith tracker variable of the
// new particle.
double Inception::ParticleTrackers(unsigned int i) const
{
    if (i < m_newvals.size()) {
        return m_newvals[i];
    } else {
        return 0.0;
    }
}

// Sets the new particle tracker variable vector.
void Inception::SetParticleTrackers(const fvector &track)
{
    m_newvals.assign(track.begin(), track.end());
}

// Sets the value of the ith tracker variable in the
// new particle.
void Inception::SetParticleTracker(unsigned int i, double track)
{
    if (i < m_mech->TrackerCount()) {
        // Ensure vector is sufficiently long.
        if (m_newvals.size() < m_mech->TrackerCount()) {
            m_newvals.resize(m_mech->TrackerCount(),0.0);
        }
        // Set value.
        m_newvals[i] = track;
    }
}


// TOTAL RATE CALCULATIONS.


// Calculates the rate of multiple inceptions given a
// vector of inceptions and an iterator to a vector of
// reals for output.
double Inception::CalcRates(double t, const Cell &sys,
                          const Geometry::LocalGeometry1d &local_geom,
                          const IcnPtrVector &icns,
                          fvector &rates, unsigned int start)
{
    IcnPtrVector::const_iterator p;
    fvector::iterator i = (rates.begin()+start);
    double sum = 0.0;
    for (p=icns.begin(); p!=icns.end(); ++p,++i) {
        *i = (*p)->Rate(t, sys, local_geom);
        sum += *i;
    }
    return sum;
}





// READ/WRITE/COPY.

// Writes the object to a binary stream.
void Inception::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        Process::Serialize(out);

        // Write new component count.
        unsigned int n = (unsigned int)m_newcomp.size();
        out.write((char*)&n, sizeof(n));

        // Write component values.
        double v(0.0);
        for (fvector::const_iterator i=m_newcomp.begin(); i!=m_newcomp.end(); ++i) {
            v = (double)*i;
            out.write((char*)&v, sizeof(v));
        }

        // Write new tracker values count.
        n = (unsigned int)m_newvals.size();
        out.write((char*)&n, sizeof(n));

        // Write tracker values.
        for (fvector::const_iterator i=m_newvals.begin(); i!=m_newvals.end(); ++i) {
            v = (double)*i;
            out.write((char*)&v, sizeof(v));
        }

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Inception::Serialize).");
    }
}

// Reads the object from a binary stream.
void Inception::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    m_newcomp.clear();
    m_newvals.clear();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        double val = 0.0;

        switch (version) {
            case 0:
                // Deserialize base class.
                Process::Deserialize(in, mech);

                // Read new component count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read component values.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_newcomp.push_back((double)val);
                }

                // Read new tracker values count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read tracker values.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_newvals.push_back((double)val);
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Inception::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Inception::Deserialize).");
    }
}
