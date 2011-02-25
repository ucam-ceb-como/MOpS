/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ParticleStats class declared in the
    swp_particle_stats.h header file.

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

#include "swp_particle_stats.h"
#include "swp_submodel_type.h"
#include "swp_particle.h"
#include "swp_particle_model.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;
using namespace Sweep::Stats;

// STATIC CONST MEMBER VARIABLES.

const std::string ParticleStats::m_statnames[ParticleStats::STAT_COUNT] = {
    std::string("Particle Count"),
    std::string("M0 (cm-3)"),
    std::string("Equiv. Sphere Diameter (nm)"),
    std::string("Collision Diameter (nm)"),
    std::string("Mobility Diameter (nm)"),
    std::string("Surface Area (cm2/cm3)"),
    std::string("Avg. Surface Area (cm2)"),
    std::string("Fv"),
    std::string("Avg. Volume (cm3)"),
    std::string("Mass (g/cm3)"),
    std::string("Avg. Mass (g)"),
    std::string("Secondary Particle Count"),
    std::string("Secondary M0 (cm-3)")
};

const IModelStats::StatType ParticleStats::m_mask[ParticleStats::STAT_COUNT] = {
    IModelStats::None, // Particle count.
    IModelStats::Sum,  // M0.
    IModelStats::Avg,  // Diameter.
    IModelStats::Avg,  // Collision diameter.
    IModelStats::Avg,  // Mobility diameter.
    IModelStats::Sum,  // Surface area.
    IModelStats::Avg,  // Avg. Surface area.
    IModelStats::Sum,  // Volume.
    IModelStats::Avg,  // Avg. volume.
    IModelStats::Sum,  // Mass.
    IModelStats::Avg,  // Avg. mass.
    IModelStats::None, // Secondary particle count.
    IModelStats::Sum   // Secondary M0.
};

const std::string ParticleStats::m_const_pslnames[ParticleStats::PSL_COUNT] = {
    std::string("Equiv. Sphere Diameter (nm)"),
    std::string("Collision Diameter (nm)"),
    std::string("Mobility Diameter (nm)"),
    std::string("Surface Area (cm2)"),
    std::string("Volume (cm3)"),
    std::string("Mass (g)"),
    std::string("Age (s)")
};


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
ParticleStats::ParticleStats()
: m_ncomp(0), m_ntrack(0), m_stats(STAT_COUNT,0.0)
{
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        m_names.push_back(m_statnames[i]);
    }

    for (unsigned int i=0; i!=PSL_COUNT; ++i) {
        m_pslnames.push_back(m_const_pslnames[i]);
    }
}

// Default constructor (public).
ParticleStats::ParticleStats(const Sweep::ParticleModel &model)
{
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        m_names.push_back(m_statnames[i]);
    }

    for (unsigned int i=0; i!=PSL_COUNT; ++i) {
        m_pslnames.push_back(m_const_pslnames[i]);
    }

    m_ncomp  = model.ComponentCount();
    m_ntrack = model.TrackerCount();
    m_stats.resize(STAT_COUNT+2*(m_ncomp+m_ntrack), 0.0);

    // Add component and tracker names to m_names and m_pslnames vectors.
    for (unsigned int i=0; i!=m_ncomp; ++i) {
        m_names.push_back(model.Components(i)->Name());
        m_names.push_back(string("Avg. ").append(model.Components(i)->Name()));
        m_pslnames.push_back(model.Components(i)->Name());
    }
    for (unsigned int i=0; i!=m_ntrack; ++i) {
        m_names.push_back(model.Trackers(i)->Name());
        m_names.push_back(string("Avg. ").append(model.Trackers(i)->Name()));
        m_pslnames.push_back(model.Trackers(i)->Name());
    }
}

// Copy constructor.
ParticleStats::ParticleStats(const ParticleStats &copy)
{
    *this = copy;
}

// Stream-reading constructor.
ParticleStats::ParticleStats(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
ParticleStats::~ParticleStats()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

ParticleStats &ParticleStats::operator=(const ParticleStats &rhs)
{
    if (this != &rhs) {
        m_stats.assign(rhs.m_stats.begin(), rhs.m_stats.end());
        m_ncomp = rhs.m_ncomp;
        m_ntrack = rhs.m_ntrack;
    }
    return *this;
}


// IMPLEMENTATION.

// Returns the number of basic particle stats.
unsigned int ParticleStats::Count() const
{
    return STAT_COUNT + 2*(m_ncomp + m_ntrack);
}

// Calculates the model stats for a single particle.
void ParticleStats::Calculate(const Particle &data)
{
    m_stats[iNP]   = 1.0;
    m_stats[iM0]   = 1.0;
    m_stats[iD]    = data.SphDiameter() * 1.0e9;  // Convert from m to nm.
    m_stats[iDcol] = data.CollDiameter() * 1.0e9; // Convert from m to nm.
    m_stats[iDmob] = data.MobDiameter() * 1.0e9;  // Convert from m to nm.
    m_stats[iS]    = data.SurfaceArea() * 1.0e4;  // Convert from m2 to cm2.
    m_stats[iS+1]  = m_stats[iS];
    m_stats[iV]    = data.Volume() * 1.0e6;       // Convert from m3 to cm3.
    m_stats[iV+1]  = m_stats[iV];
    m_stats[iM]    = data.Mass() * 1.0e3;         // Convert from kg to g.
    m_stats[iM+1]  = m_stats[iM+1];
    m_stats[i2NP]  = 0.0; // Not used here, because particle assumed not to be secondary
    m_stats[i2M0]  = 0.0; // Not used here, because particle assumed not to be secondary

    // Add component and tracker values to stats.
    fvector::iterator i = m_stats.begin()+STAT_COUNT;
    for (unsigned int j=0; j!=m_ncomp; ++j,++i) {
        *i     = data.Composition(j);
        *(++i) = data.Composition(j);
    }
    for (unsigned int j=0; j!=m_ncomp; ++j,++i) {
        *i     = data.Values(j);
        *(++i) = data.Values(j);
    }
}

// Calculates the model stats for a particle ensemble.
void ParticleStats::Calculate(const Ensemble &e, real scale, real secondary_scale)
{
    m_stats.assign(m_stats.size(), 0.0);

    // Loop over all particles, getting the stats from each.
    Ensemble::const_iterator ip;
    unsigned int n = 0;
    for (ip=e.begin(); ip!=e.end(); ++ip) {
        real sz = (*ip)->Property(m_statbound.PID);
        // Check if the value of the property is within the stats bound
        if ((m_statbound.Lower < sz) && (sz < m_statbound.Upper) ) {
            // Sum stats from this particle.
            m_stats[iD]    += (*ip)->SphDiameter() * 1.0e9;  // Convert from m to nm.
            m_stats[iDcol] += (*ip)->CollDiameter() * 1.0e9; // Convert from m to nm.
            m_stats[iDmob] += (*ip)->MobDiameter() * 1.0e9;  // Convert from m to nm.
            m_stats[iS]    += (*ip)->SurfaceArea() * 1.0e4;  // Convert from m2 to cm2.
            m_stats[iS+1]  += (*ip)->SurfaceArea() * 1.0e4;  // Convert from m2 to cm2.
            m_stats[iV]    += (*ip)->Volume() * 1.0e6;       // Convert from m3 to cm3.
            m_stats[iV+1]  += (*ip)->Volume() * 1.0e6;       // Convert from m3 to cm3.
            m_stats[iM]    += (*ip)->Mass() * 1.0e3;         // Convert from kg to g.
            m_stats[iM+1]  += (*ip)->Mass() * 1.0e3;         // Convert from kg to g.

            // Sum component and tracker values.
            fvector::iterator i = m_stats.begin()+STAT_COUNT;
            for (unsigned int j=0; j!=m_ncomp; ++j,++i) {
                *i     += (*ip)->Composition(j);
                *(++i) += (*ip)->Composition(j);
            }
            for (unsigned int j=0; j!=m_ntrack; ++j,++i) {
                *i     += (*ip)->Values(j);
                *(++i) += (*ip)->Values(j);
            }
            ++n;
        }
    }

    // Get the particle count.
    m_stats[iNP] = (real)e.Count();
    m_stats[iM0] = (real)n;
    const real invNumPart = (n>0) ? 1.0 / n : 0.0;

    // Scale the summed stats and calculate the averages.
    for (unsigned int i=1; i!=STAT_COUNT - 2; ++i) {
        if (m_mask[i] == Sum) {
            m_stats[i] *= (scale * 1.0e-6); // Convert scale from 1/m3 to 1/cm3.
        } else {
            m_stats[i] *= invNumPart;
        }
    }

    // Secondary population quantities
    m_stats[i2NP] = e.SecondaryCount();
    m_stats[i2M0] = e.SecondaryCount() * (secondary_scale * 1.0e-6); // Convert scale from 1/m3 to 1/cm3.



    // Scale and calculate averages for components and tracker
    // variables.
    for (unsigned int i=STAT_COUNT; i!=Count(); ++i) {
        m_stats[i]   *= (scale * 1.0e-6); // Convert scale from 1/m3 to 1/cm3.
        m_stats[++i] *= invNumPart;
    }
}

// Returns a vector containing the stats.
const fvector &ParticleStats::Get(void) const
{
    return m_stats;
}

// Returns a vector containing the stats.
void ParticleStats::Get(fvector &stats, unsigned int start) const
{
    // Get an iterator to the first point of insertion in the
    // output stats array.
    fvector::iterator i;
    if (start < stats.size()) {
        i = stats.begin()+start;
    } else {
        i = stats.end();
    }

    // Add stats to output array until the end of that array is
    // reached or we have run out of stats to add.
    fvector::const_iterator j = m_stats.begin();
    for (; (i!=stats.end()) && (j!=m_stats.end()); ++i,++j) {
        *i = *j;
    }

    // If we still have stats to add, but the output array is full,
    // then we need to add elements to the end of the array.
    for (; j!=m_stats.end(); ++j) {
        stats.push_back(*j);
    }
}

// Returns a vector containing the stat names.
const std::vector<std::string> &ParticleStats::Names(void) const
{
    return m_names;
}

// Adds to a vector containing stat names.
void ParticleStats::Names(std::vector<std::string> &names,
                          unsigned int start) const
{
    // Get an iterator to the first point of insertion in the
    // output names array.
    std::vector<std::string>::iterator i;
    if (start < names.size()) {
        i = names.begin()+start;
    } else {
        i = names.end();
    }

    // Add stats to output array until the end of that array is
    // reached or we have run out of names to add.
    std::vector<std::string>::const_iterator j = m_names.begin();
    for (; (i!=names.end()) && (j!=m_names.end()); ++i,++j) {
        *i = *j;
    }

    // If we still have names to add, but the output array is full,
    // then we need to add elements to the end of the array.
    while(j!=m_names.end()) {
        names.push_back(*j);
        ++j;
    }
}


// AVAILABLE BASIC STATS.

// Returns the particle count.
real ParticleStats::PCount(void) const {return m_stats[iNP];}

//! Returns the secondary particle count.
real ParticleStats::SecondaryPCount(void) const {return m_stats[i2NP];}

// Returns the number density.
real ParticleStats::M0(void) const {return m_stats[iM0];}

// Returns the secondary number density.
real ParticleStats::SecondaryM0(void) const {return m_stats[i2M0];}

// Returns the avg. equiv. sphere diameter.
real ParticleStats::AvgDiam(void) const {return m_stats[iD];}

// Returns the avg. collision diameter.
real ParticleStats::AvgCollDiam(void) const {return m_stats[iDcol];}

// Returns the avg. mobility diameter.
real ParticleStats::AvgMobDiam(void) const {return m_stats[iDmob];}

// Returns the total surface area.
real ParticleStats::SurfaceArea(void) const {return m_stats[iS];}

// Returns the avg. surface area.
real ParticleStats::AvgSurfaceArea(void) const {return m_stats[iS+1];}

// Returns the total volume.
real ParticleStats::Fv(void) const {return m_stats[iV];}

// Returns the average volume.
real ParticleStats::AvgVolume(void) const {return m_stats[iV+1];}

// Returns the total mass.
real ParticleStats::Mass(void) const {return m_stats[iM];}

// Returns the average mass.
real ParticleStats::AvgMass(void) const {return m_stats[iM+1];}


// PARTICLE SIZE LISTS.

// Returns the number of PSL output variables.
unsigned int ParticleStats::PSL_Count(void) const
{
    return PSL_COUNT + m_ncomp + m_ntrack;
}

// Returns a vector of PSL variable names.
void ParticleStats::PSL_Names(std::vector<std::string> &names,
                              unsigned int start) const
{
    // Get an iterator to the first point of insertion in the
    // output names array.
    std::vector<std::string>::iterator i;
    if (start < names.size()) {
        i = names.begin()+start;
    } else {
        i = names.end();
    }

    // Add stats to output array until the end of that array is
    // reached or we have run out of names to add.
    std::vector<std::string>::const_iterator j = m_pslnames.begin();
    for (; (i!=names.end()) && (j!=m_pslnames.end()); ++i,++j) {
        *i = *j;
    }

    // If we still have names to add, but the output array is full,
    // then we need to add elements to the end of the array.
    while(j!=m_pslnames.end()) {
        names.push_back(*j);
        ++j;
    }
}

// Returns the PSL entry for the given particle.
void ParticleStats::PSL(const Sweep::Particle &sp, real time,
                        fvector &psl, unsigned int start) const
{
    // Resize vector if too small.
    if (start+PSL_Count()-1 >= psl.size()) {
        psl.resize(start+PSL_Count(), 0.0);
    }

    // Get an iterator to the first point of insertion in the
    // output stats array.
    fvector::iterator j = psl.begin()+start-1;

    // Get the PSL stats.
    *(++j) = sp.SphDiameter() * 1.0e9;  // m to nm.
    *(++j) = sp.CollDiameter() * 1.0e9; // m to nm.
    *(++j) = sp.MobDiameter() * 1.0e9;  // m to nm.
    *(++j) = sp.SurfaceArea() * 1.0e4;  // m2 to cm2.
    *(++j) = sp.Volume() * 1.0e6;       // m3 to cm3.
    *(++j) = sp.Mass() * 1.0e3;         // kg to g.
    *(++j) = time - sp.CreateTime();    // Particle age.

    // Get component and tracker stats.
    for (unsigned int k=0; k!=m_ncomp; ++k) {
        *(++j) = sp.Composition(k);
    }
    for (unsigned int k=0; k!=m_ntrack; ++k) {
        *(++j) = sp.Values(k);
    }
}


// READ/WRITE/COPY.

// Creates a copy of the object.
ParticleStats *const ParticleStats::Clone(void) const
{
    return new ParticleStats(*this);
}

// Returns the model data type.  Used to identify different models
// and for serialisation.
unsigned int ParticleStats::ID(void) const
{
    return (unsigned int)SubModels::BasicModel_ID;
}

// Writes the object to a binary stream.
void ParticleStats::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output the component count.
        unsigned int n = (unsigned int)m_ncomp;
        out.write((char*)&n, sizeof(n));

        // Output tracker variable count.
        n = (unsigned int)m_ntrack;
        out.write((char*)&n, sizeof(n));

        // Output number of stats in vector.
        n = (unsigned int)m_stats.size();
        out.write((char*)&n, sizeof(n));

        // Output stats.
        double v = 0.0;
        for (unsigned int i=0; i!=n; ++i) {
            v = m_stats[i];
            out.write((char*)&v, sizeof(v));
        }

        // Output number of stat names in vector.
        n = (unsigned int)m_names.size();
        out.write((char*)&n, sizeof(n));

        // Output names.
        for (unsigned int i=0; i!=n; ++i) {
            // Write the length of the stat name to the stream.
            unsigned int m = m_names[i].length();
            out.write((char*)&m, sizeof(m));

            // Write the stat name to the stream.
            out.write(m_names[i].c_str(), m);
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ParticleStats::Serialize).");
    }
}

// Reads the object from a binary stream.
void ParticleStats::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    // TODO:  Deserialize ParticleStats should reset to state with no components
    //        or tracker variables in the first instance.
    m_stats.clear();
    m_names.clear();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        double val = 0.0;
        char *name = NULL;

        switch (version) {
            case 0:

                // Read the component count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_ncomp = n;

                // Read tracker variable count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_ntrack = n;

                // Read number of stats in vector.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read stats.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_stats.push_back((real)val);
                }

                // Read number of stat names in vector.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read names.
                for (unsigned int i=0; i!=n; ++i) {
                    // Read the length of the element name.
                    unsigned int m = 0;
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));

                    // Read the element name.
                    name = new char[m];
                    in.read(name, m);
                    m_names.push_back(string(name, m));
                    delete [] name;
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ParticleStats::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ParticleStats::Deserialize).");
    }
}
