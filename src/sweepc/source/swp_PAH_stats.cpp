/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the PAHStats class declared in the
    swp_PAH_stats.h header file.

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

#include "swp_PAH_stats.h"
#include "swp_aggmodel_type.h"
#include "swp_particle.h"
#include "swp_PAH_primary.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Stats;
using namespace std;

// STATIC CONST MEMBER VARIABLES.

const std::string PAHStats::m_statnames[PAHStats::STAT_COUNT] = {
    std::string("Avg. Number of PAHs"),
    std::string("Number of PAHs cm-3"),
    std::string("Surface double Part"),
    std::string("Avg. Mass double Part"),
    std::string("Avg. PAH double Part"),
    std::string("Avg. PAH Collision Diameter"),
	std::string("Avg. Number of Carbon Atoms"),
	std::string("Avg. Number of Hydrogen Atoms"),
	std::string("Avg. Number of Edge Carbon Atoms"),
	std::string("Avg. Number of Rings"),
    //std::string("Avg. Coalesc Threshold"), //should be the same thing as Avg. Sintering Level
	std::string("Avg. Sint Level"),
    std::string("Num Primaries double Part"),
};

const IModelStats::StatType PAHStats::m_mask[PAHStats::STAT_COUNT] = {
    IModelStats::Avg,  // Avg. Number of PAHs.
    IModelStats::Sum,  // Number of PAHs.
    IModelStats::Sum,  // Surface double Part
    IModelStats::Avg,  // Avg. Mass double Part
    IModelStats::Avg,  // Avg. PAH double Part
    IModelStats::Avg,  // Avg. PAH Collision Diameter
	IModelStats::Avg,  // Avg. Number of Carbon atoms
	IModelStats::Avg,  // Avg. Number of Hydrogen atoms
	IModelStats::Avg,  // Avg. Number of Edge Carbon Atoms
	IModelStats::Avg,  // Avg. Number of Rings
    //IModelStats::Avg,  // Avg. Coalesc Threshold, should be the same thing as Avg. Sintering Level
	IModelStats::Avg,  // Avg. Sintering Level, if primary coordinates are tracked
    IModelStats::Avg,  // Num Primaries double Part
};

const std::string PAHStats::m_const_pslnames[PAHStats::PSL_COUNT] = {
    std::string("Number of PAHs"),
    std::string("PAH Diameter"),
	std::string("Number of Carbon atoms"),
    std::string("Number of Hydrogen atoms"),
	std::string("Number primaries"),
	std::string("Number of Edge Carbon Atoms"),
	std::string("Number of Rings"),
	std::string("sqrt(LW)"),
	std::string("LdivW"),
	std::string("Avg. primary diameter"),
    std::string("Radius of gyration"),
    std::string("fdim"),
};



// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
PAHStats::PAHStats()
: m_stats(STAT_COUNT,0.0)
{
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        m_names.push_back(m_statnames[i]);
    }

    for (unsigned int i=0; i!=PSL_COUNT; ++i) {
        m_pslnames.push_back(m_const_pslnames[i]);
    }


}

// Copy constructor.
PAHStats::PAHStats(const Sweep::Stats::PAHStats &copy)
{
    *this = copy;
}

// Stream-reading constructor.
PAHStats::PAHStats(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
PAHStats::~PAHStats()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

PAHStats &PAHStats::operator=(const Sweep::Stats::PAHStats &rhs)
{
    if (this != &rhs) {
        m_stats.assign(rhs.m_stats.begin(), rhs.m_stats.end());
    }
    return *this;
}


// IMPLEMENTATION.

// Returns the number of basic particle stats.
unsigned int PAHStats::Count() const
{
    return STAT_COUNT;
}

// Calculates the model stats for a particle ensemble.
void PAHStats::Calculate(const Ensemble &e, double scale)
{
    // Empty the stats array.
    fill(m_stats.begin(), m_stats.end(), 0.0);

    // Forward define sum of 'double particle' weights
    // used to count particles with more than one PAH
	// double part properties => particle only considered, sum of X / num of particle
    double wtreal = 0.0;

    // Loop over all particles, getting the stats from each.
    Ensemble::const_iterator ip;

    for (ip=e.begin(); ip!=e.end(); ++ip) {
		const AggModels::PAHPrimary *pah = NULL;
			pah = dynamic_cast<const AggModels::PAHPrimary*>((*(*ip)).Primary());
        double sz = (*ip)->Property(m_statbound.PID);
        double wt = (*ip)->getStatisticalWeight();

        // Check if the value of the property is within the stats bound
        if ((m_statbound.Lower < sz) && (sz < m_statbound.Upper) ) {
            // Sum stats from this particle.
			m_stats[iNPAH]    		+= pah->NumPAH() * wt;
			m_stats[iPAHD]    		+= pah->PAHCollDiameter()*1e9 * wt;
			m_stats[iNCARB]	  		+= pah->NumCarbon() * wt;
			m_stats[iNHYDROGEN]	  	+= pah->NumHydrogen() * wt;
            m_stats[iNEDGEC]	  	+= pah->NumEdgeC() * wt;
            m_stats[iNRINGS]	  	+= pah->NumRings() * wt;
            m_stats[iNPAH+1]    	+= pah->NumPAH() * wt; //used to calculate sum of Number of PAHs.
			//m_stats[iCOAL]          += pah->AvgCoalesc() * wt; //should be the same thing as AvgSinter()
            m_stats[iSINT]    		+= pah->AvgSinter() * wt; //AvgSinter() can return coalescence level (in the old model, coordinates are not tracked)
			                                                  //or return sintering level (in the new model, coordinates are not tracked)
            if (pah->NumPAH() > 1)
            {
                wtreal += wt;
                m_stats[iPARTSURF]+=(*ip)->SurfaceArea() * wt;
                m_stats[iNPRIM]+= pah->Numprimary() * wt;
                m_stats[iPARTMASS]+=(*ip)->Primary()->Mass() * wt;
                m_stats[iNAVGPAH]+= pah->NumPAH() * wt;  //used to calculate Avg. PAH double Part
            }
        }
    }

    // Calculate total weight
    double invTotalWeight = e.Count()>0 ? 1.0/e.GetSum(iW) : 0.0;
    double invRealWeight = wtreal==0.0 ? 1.0/wtreal : 0.0;

    // Scale the summed stats and calculate the averages.
    for (unsigned int i=0; i!=STAT_COUNT; ++i) {
        if (m_mask[i] == Sum) {
            m_stats[i] *= (scale * 1.0e-6); // Convert scale from 1/m3 to 1/cm3.
        } else {
            if (i==iNPRIM || i==iPARTMASS || i==iNAVGPAH)
            {
                if (wtreal>0)
                    m_stats[i] *= invRealWeight;
            }
            else
                m_stats[i] *= invTotalWeight;
        }
    }

}

// Returns a vector containing the stats.
const fvector &PAHStats::Get(void) const
{
    return m_stats;
}

// Returns a vector containing the stats.
void PAHStats::Get(fvector &stats, unsigned int start) const
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
const std::vector<std::string> &PAHStats::Names(void) const
{
    return m_names;
}

// Adds to a vector containing stat names.
void PAHStats::Names(std::vector<std::string> &names,
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




// PARTICLE SIZE LISTS.

// Returns the number of PSL output variables.
unsigned int PAHStats::PSL_Count(void) const
{
    return PSL_COUNT;
}

// Returns a vector of PSL variable names.
void PAHStats::PSL_Names(std::vector<std::string> &names,
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
void PAHStats::PSL(const Sweep::Particle &sp, double time,
                       fvector &psl, unsigned int start) const
{
    // Resize vector if too small.
    if (start+PSL_Count()-1 >= psl.size()) {
        psl.resize(start+PSL_Count(), 0.0);
    }

    // Get an iterator to the first point of insertion in the
    // output stats array.
    fvector::iterator j = psl.begin()+start-1;

    // This should succeed, there is no way the PAH stats should be called
    // unless we also have a PAHPrimary
    const AggModels::PAHPrimary* const pah =
            dynamic_cast<const AggModels::PAHPrimary*>(sp.Primary());

    // Get the PSL stats.
    if (pah != NULL) {
		*(++j) = (double)(pah->NumPAH());
		*(++j) = (double)(pah->PAHCollDiameter())*1e9;			//convert to nm
		*(++j) = (double) (pah->NumCarbon());
		*(++j) = (double) (pah->NumHydrogen());
		*(++j) = (double) (pah->Numprimary());
        *(++j) = (double) (pah->NumEdgeC());
        *(++j) = (double) (pah->NumRings());
        *(++j) = (double) (pah->sqrtLW());
		*(++j) = (double) (pah->LdivW());
		*(++j) = (double) (pah->PrimaryDiam())*1e9/(double)(pah->Numprimary());//convert to nm
        //*(++j) = (double) (pah->Rg())*1e9; //use function GetRadiusOfGyration() to calculate
		*(++j) = (double)(pah->GetRadiusOfGyration())*1e9;
        *(++j) = (double) (pah->Fdim());

    } else {
        fill (j+1, j+2, 0.0);
    }
}

// READ/WRITE/COPY.

// Creates a copy of the object.
PAHStats *const PAHStats::Clone(void) const
{
    return new PAHStats(*this);
}

// Returns the model data type.  Used to identify different models
// and for serialisation.
unsigned int PAHStats::ID(void) const
{
    return (unsigned int)AggModels::PAH_KMC_ID;
}

// Writes the object to a binary stream.
void PAHStats::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output number of stats in vector.
        unsigned int n = (unsigned int)m_stats.size();
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
                               "(Sweep, PAHStats::Serialize).");
    }
}

// Reads the object from a binary stream.
void PAHStats::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
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
                    // Read the length of the name.
                    unsigned int m = 0;
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));

                    // Read the name.
                    name = new char[m];
                    in.read(name, m);
                    m_names.push_back(string(name, m));
                    delete [] name;
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, PAHStats::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PAHStats::Deserialize).");
    }
}

// Get primary particle details and connectivity
void PAHStats::PrintPrimary(const Sweep::Particle &sp, std::vector<fvector> &nodes, std::vector<fvector> &primaries, int k) const
{
	const AggModels::PAHPrimary* const prim =
		dynamic_cast<const AggModels::PAHPrimary *>(sp.Primary());

	if (prim != NULL) {
		prim->PrintPrimary(nodes, primaries, k);
	}
}
