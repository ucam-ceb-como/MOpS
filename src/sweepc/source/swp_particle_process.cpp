/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ParticleProcess class declared in the
    swp_particle_process.h header file.

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

#include "swp_particle_process.h"
#include "swp_mechanism.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ParticleProcess::ParticleProcess()
: Process(), m_defer(false)
{
}

// Initialising constructor.
ParticleProcess::ParticleProcess(const Sweep::Mechanism &mech)
: Process(mech), m_defer(false)
{
}

// Copy constructor.
ParticleProcess::ParticleProcess(const ParticleProcess &copy)
{
    *this = copy;
}

// Default destructor.
ParticleProcess::~ParticleProcess()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
ParticleProcess &ParticleProcess::operator=(const ParticleProcess &rhs)
{
    if (this != &rhs) {
        Process::operator=(rhs);
        m_defer = rhs.m_defer;
        m_dcomp = rhs.m_dcomp;
        m_dvals = rhs.m_dvals;
    }
    return *this;
}


// DEFERRED PROCESSES.

// Returns TRUE if process should be deferred, otherwise false.
bool ParticleProcess::IsDeferred(void) const
{
    return m_defer;
}

// Sets the process to be deferred or not.
void ParticleProcess::SetDeferred(bool defer)
{
    m_defer = defer; 
    if (m_mech!=NULL) m_mech->CheckDeferred();
}


// CHANGES TO PARTICLE ON PROCESS OCCURANCE.

// Returns the composition vector of the new particle.
const fvector &ParticleProcess::CompChange(void) const
{
    return m_dcomp;
}

// Returns the amount of the ith component of the new particle.
double ParticleProcess::CompChange(unsigned int i) const
{
    if (i < m_dcomp.size()) {
        return m_dcomp[i];
    } else {
        return 0.0;
    }
}

// Sets the particle composition vector.
void ParticleProcess::SetCompChange(const fvector &comp)
{
    m_dcomp.assign(comp.begin(), comp.end());
}

// Sets the amount of the ith component in the new particle.
void ParticleProcess::SetCompChange(unsigned int i, double comp)
{
    if (i < m_mech->ComponentCount()) {
        // Ensure vector is sufficiently long.
        if (m_dcomp.size() < m_mech->ComponentCount()) {
            m_dcomp.resize(m_mech->ComponentCount(),0.0);
        }
        // Set value.
        m_dcomp[i] = comp;
    }
}

// Returns the tracker variable vector of the new particle.
const fvector &ParticleProcess::TrackChange(void) const
{
    return m_dvals;
}

// Returns the value of the ith tracker variable of the
// new particle.
double ParticleProcess::TrackChange(unsigned int i) const
{
    if (i < m_dvals.size()) {
        return m_dvals[i];
    } else {
        return 0.0;
    }   
}

// Sets the new particle tracker variable vector.
void ParticleProcess::SetTrackChange(const fvector &track)
{
    m_dvals.assign(track.begin(), track.end());
}

// Sets the value of the ith tracker variable in the
// new particle.
void ParticleProcess::SetTrackChange(unsigned int i, double track)
{
    if (i < m_mech->TrackerCount()) {
        // Ensure vector is sufficiently long.
        if (m_dvals.size() < m_mech->TrackerCount()) {
            m_dvals.resize(m_mech->TrackerCount(),0.0);
        }
        // Set value.
        m_dvals[i] = track;
    }
}


// RATE CALCULATION.

// Calculates the rates of multiple particle processes.
double ParticleProcess::CalcRates(double t, const Cell &sys,
                                const Geometry::LocalGeometry1d &local_geom,
                                const PartProcPtrVector &proc,
                                fvector &rates, unsigned int start)
{
    PartProcPtrVector::const_iterator p;
    fvector::iterator i = (rates.begin()+start);
    double sum = 0.0;
    for (p=proc.begin(); p!=proc.end(); ++p,++i) {
        *i = (*p)->Rate(t, sys, local_geom);
        sum += *i;
    }
    return sum;
}

// Return rate constant and chemistry part for hybrid method
double ParticleProcess::Rate(double t, const Cell &sys) const
{
	std::cout << "Only used with surface growth for hybrid particle model\n";
	return -1;
}

// Do surface growth gas-phase adjustment for hybrid method
int ParticleProcess::Perform(double t, Cell &sys, rng_type &rng, unsigned int n) const
{
	std::cout << "Only used with surface growth for hybrid particle model\n";
	return -1;
}

// Do surface growth gas-phase adjustment for hybrid method
int ParticleProcess::Perform(double t, Cell &sys, Particle &sp, rng_type &rng, unsigned int n, bool isParticleNumberUpdate) const
{
	std::cout << "Only used with surface growth for hybrid particle model\n";
	return -1;
}

// READ/WRITE/COPY.

// Writes the object to a binary stream.
void ParticleProcess::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        Process::Serialize(out);

        // Write if the process is deferred.
        if (m_defer) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write component count.
        unsigned int n = (unsigned int)m_dcomp.size();
        out.write((char*)&n, sizeof(n));

        // Write component changes.
        double v = 0.0;
        for (fvector::const_iterator i=m_dcomp.begin(); i!=m_dcomp.end(); ++i) {
            v = (double)*i;
            out.write((char*)&v, sizeof(v));
        }

        // Write tracker values count.
        n = (unsigned int)m_dvals.size();
        out.write((char*)&n, sizeof(n));

        // Write tracker changes.
        for (fvector::const_iterator i=m_dvals.begin(); i!=m_dvals.end(); ++i) {
            v = (double)*i;
            out.write((char*)&v, sizeof(v));
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ParticleProcess::Serialize).");
    }
}

// Reads the object from a binary stream.
void ParticleProcess::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    m_dcomp.clear();
    m_dvals.clear();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;
        unsigned int n = 0;

        switch (version) {
            case 0:
                // Deserialize base class.
                Process::Deserialize(in, mech);

                // Read if the process is deferred.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_defer = (n==1);

                // Read component count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read component changes.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_dcomp.push_back(val);
                }

                // Read tracker values count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read tracker changes.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_dvals.push_back(val);
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ParticleProcess::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ParticleProcess::Deserialize).");
    }
}
