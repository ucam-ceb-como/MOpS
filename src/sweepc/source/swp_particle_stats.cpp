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
#include "swp_particle.h"
#include "swp_particle_model.h"

#include <stdexcept>

using namespace Sweep;
using namespace std;
using namespace Sweep::Stats;

// STATIC CONST MEMBER VARIABLES.

const std::string ParticleStats::m_statnames[ParticleStats::STAT_COUNT] = {
    std::string("Particle Count"),
    std::string("M0 (m-3)"),
    std::string("Equiv. Sphere Diameter (m)"),
    std::string("Collision Diameter (m)"),
    std::string("Mobility Diameter (m)"),
    std::string("Surface Area (m2/m3)"),
    std::string("Avg. Surface Area (m2)"),
    std::string("Fv (-)"),
    std::string("Avg. Volume (m3)"),
    std::string("Mass (kg/m3)"),
    std::string("Avg. Mass (kg)"),
    std::string("Mass2 (kg2/m6)"),
    std::string("Mass3 (kg3/m9)"),
    std::string("Avg num coags (-)"),
    std::string("Max num coags (-)")
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
    IModelStats::Sum,  // Mass2.
    IModelStats::Sum,  // Mass3.
    IModelStats::Avg,  // Average number of coagulations
    IModelStats::None  // Max number of coagulations
};

const std::string ParticleStats::m_const_pslnames[ParticleStats::PSL_COUNT] = {
    std::string("Equiv. Sphere Diameter (nm)"),
    std::string("Collision Diameter (nm)"),
    std::string("Mobility Diameter (nm)"),
    std::string("Surface Area (cm2)"),
    std::string("Volume (cm3)"),
    std::string("Mass (g)"),
    std::string("Age (s)"),
    std::string("Stat. Weight (-)"),
    std::string("Num coags (-)")
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

// Calculates the model stats for a particle ensemble.
void ParticleStats::Calculate(const Ensemble &e, double scale)
{
    m_stats.assign(m_stats.size(), 0.0);

    // Find out the maximum number of coags that has been experienced
    // by one of the particles in this ensemble.  Find the max
    // using this int before copying into the stats array
    unsigned int maxNumCoag = 0;


    // Loop over all particles, getting the stats from each.
    Ensemble::const_iterator ip;
    for (ip=e.begin(); ip!=e.end(); ++ip) {
        double sz = (*ip)->Property(m_statbound.PID);
        // Check if the value of the property is within the stats bound
        if ((m_statbound.Lower < sz) && (sz < m_statbound.Upper) ) {
            const double wt = (*ip)->getStatisticalWeight();
            // Sum stats from this particle.
            m_stats[iM0]   += wt;
            m_stats[iD]    += (*ip)->SphDiameter() * wt;
            m_stats[iDcol] += (*ip)->CollDiameter() * wt;
            m_stats[iDmob] += (*ip)->MobDiameter() * wt;
            m_stats[iS]    += (*ip)->SurfaceArea() * wt;
            m_stats[iS+1]  += (*ip)->SurfaceArea() * wt;
            m_stats[iV]    += (*ip)->Volume() * wt;
            m_stats[iV+1]  += (*ip)->Volume() * wt;
            const double m = (*ip)->Mass();
            m_stats[iM]    += m * wt;
            m_stats[iM+1]  += m * wt;
            m_stats[iM2]   += m * m * wt;
            m_stats[iM3]   += m * m * m * wt;

            // Coagulations experienced by this particle
            const unsigned int coagCount = (*ip)->getCoagCount();
            m_stats[iCoag] += coagCount * wt;

            // Check if this the highest count so far
            if(coagCount > maxNumCoag)
                maxNumCoag = coagCount;

            // Sum component and tracker values.
            fvector::iterator i = m_stats.begin()+STAT_COUNT;
            for (unsigned int j=0; j!=m_ncomp; ++j,++i) {
                *i     += (*ip)->Composition(j) * wt;
                *(++i) += (*ip)->Composition(j) * wt;
            }
            for (unsigned int j=0; j!=m_ntrack; ++j,++i) {
                *i     += (*ip)->Values(j) * wt;
                *(++i) += (*ip)->Values(j) * wt;
            }
        }
    }

	// Add contributions from hybrid particle-number/particle model
	if (e.GetTotalParticleNumber() != 0)
	{
		double wt = e.GetTotalParticleNumber();
		double dpri = e.GetTotalDiameter();
		double surf = PI * e.GetTotalDiameter2();
		double vol = (PI / 6.0) * e.GetTotalDiameter3();
		double mass = e.GetTotalMass();
		double mass2 = e.GetTotalMass2();
		double mass3 = e.GetTotalMass3();
		double comp = e.GetTotalComponent();

		// Sum component and tracker values.
		fvector::iterator i = m_stats.begin() + STAT_COUNT;
		for (unsigned int j = 0; j != m_ncomp; ++j, ++i) {
			*i += comp;
			*(++i) += comp;
		}

		m_stats[iM0] += wt;
		m_stats[iD] += dpri;
		m_stats[iDcol] += dpri;
		m_stats[iDmob] += dpri;
		m_stats[iS] += surf;
		m_stats[iS + 1] += surf;
		m_stats[iV] += vol;
		m_stats[iV + 1] += vol;
		m_stats[iM] += mass;
		m_stats[iM + 1] += mass;
		m_stats[iM2] += mass2;
		m_stats[iM3] += mass3;
	}


    // Get the particle count.
    m_stats[iNP] = (double)e.Count();

    // also the maximum number of coag events that has occurred to one particle
    m_stats[iMaxCoag] = static_cast<double>(maxNumCoag);

    // Now calculate the denominator for the average per particle quantities.
    // Note that m_stats[iM0] at this point does not in fact contain an M0 value,
    // since it has not yet been scaled by sample volume.  This is intentional
    // since here one should divide by the total statistical weight of all particles.
    const double invWeight = ((e.Count() + e.GetTotalParticleNumber())>0) ? 1.0 / m_stats[iM0] : 0.0;

    // Scale the summed stats and calculate the averages,
    for (unsigned int i=1; i!=STAT_COUNT; ++i) {
        if (m_mask[i] == Sum) {
            m_stats[i] *= scale;
        } else if(m_mask[i] == Avg){
            m_stats[i] *= invWeight;
        }
    }

    // Scale and calculate averages for components and tracker
    // variables.
    for (unsigned int i=STAT_COUNT; i!=Count(); ++i) {
        m_stats[i]   *= (scale * 1.0e-6); // Convert scale from 1/m3 to 1/cm3.
        m_stats[++i] *= invWeight;
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
double ParticleStats::PCount(void) const {return m_stats[iNP];}

// Returns the number density.
double ParticleStats::M0(void) const {return m_stats[iM0];}

// Returns the avg. equiv. sphere diameter.
double ParticleStats::AvgDiam(void) const {return m_stats[iD];}

// Returns the avg. collision diameter.
double ParticleStats::AvgCollDiam(void) const {return m_stats[iDcol];}

// Returns the avg. mobility diameter.
double ParticleStats::AvgMobDiam(void) const {return m_stats[iDmob];}

// Returns the total surface area.
double ParticleStats::SurfaceArea(void) const {return m_stats[iS];}

// Returns the avg. surface area.
double ParticleStats::AvgSurfaceArea(void) const {return m_stats[iS+1];}

// Returns the total volume.
double ParticleStats::Fv(void) const {return m_stats[iV];}

// Returns the average volume.
double ParticleStats::AvgVolume(void) const {return m_stats[iV+1];}

// Returns the total mass.
double ParticleStats::Mass(void) const {return m_stats[iM];}

// Returns the average mass.
double ParticleStats::AvgMass(void) const {return m_stats[iM+1];}


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
void ParticleStats::PSL(const Sweep::Particle &sp, double time,
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
    *(++j) = sp.getStatisticalWeight(); // Statistical weight
    *(++j) = sp.getCoagCount();         // Number of coag events for this particle in current cell

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
    return (unsigned int)AggModels::Spherical_ID;
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
                    m_stats.push_back((double)val);
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

void ParticleStats::PrintPrimary(const Sweep::Particle &sp, std::vector<fvector> &nodes, std::vector<fvector> &primaries, int k) const
{
}
