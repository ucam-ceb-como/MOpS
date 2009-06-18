/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ARSSC_Process class declared in the
    swp_arssc_process.h header file.

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

#include "swp_arssc_process.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ARSSC_Process::ARSSC_Process()
: m_sites(SubModels::ARSSC_Model::SiteTypeCount,0.0), m_upd_count(0), 
  m_upd_parent(SubModels::ARSSC_Model::InvalidSite), m_use_parent_wts(true),
  m_neigh_wts(4, 0.25)
{
}

// Copy constructor.
ARSSC_Process::ARSSC_Process(const ARSSC_Process &copy)
{
    *this = copy;
}

// Stream-reading constructor.
ARSSC_Process::ARSSC_Process(std::istream &in)
{
    Deserialize(in);
}

// Default destructor.
ARSSC_Process::~ARSSC_Process(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
ARSSC_Process &ARSSC_Process::operator=(const ARSSC_Process &rhs)
{
    if (this != &rhs) {
        m_sites.assign(rhs.m_sites.begin(), rhs.m_sites.end());
        m_upd_count      = rhs.m_upd_count;
        m_upd_parent     = rhs.m_upd_parent;
        m_use_parent_wts = rhs.m_use_parent_wts;
        m_neigh_wts.assign(rhs.m_neigh_wts.begin(), rhs.m_neigh_wts.end());
    }
    return *this;
}


// SITE COUNTS.

// Returns the change in the count of the given site type.
real ARSSC_Process::SiteCount(SubModels::ARSSC_Model::SiteType type) const
{
    if (type != SubModels::ARSSC_Model::InvalidSite) {
        return m_sites[(fvector::size_type)type];
    } else {
        return 0.0;
    }
}

// Returns the change in the number of free-edges.
real ARSSC_Process::FreeEdgeCount(void) const
{
    return m_sites[(fvector::size_type)SubModels::ARSSC_Model::FreeEdge];
}

// Returns the change in the number of armchairs.
real ARSSC_Process::ArmchairCount(void) const
{
    return m_sites[(fvector::size_type)SubModels::ARSSC_Model::Armchair];
}

// Returns the change in the number of zigzags.
real ARSSC_Process::ZigzagCount(void) const
{
    return m_sites[(fvector::size_type)SubModels::ARSSC_Model::Zigzag];
}

// Returns the change in the number of bays.
real ARSSC_Process::BayCount(void) const
{
    return m_sites[(fvector::size_type)SubModels::ARSSC_Model::Bay];
}

// Returns the change in the number of R5s.
real ARSSC_Process::R5Count(void) const
{
    return m_sites[(fvector::size_type)SubModels::ARSSC_Model::R5];
}

// Sets the change in the site count of the given site type.
void ARSSC_Process::SetSiteCount(SubModels::ARSSC_Model::SiteType type, real n)
{
    if (type != SubModels::ARSSC_Model::InvalidSite) {
        m_sites[(fvector::size_type)type] = n;
    }
}



// PARENT SITE FOR NEIGHBOUR UPDATES.

// Returns the parent site used for neighbour updates.
SubModels::ARSSC_Model::SiteType ARSSC_Process::ParentSite(void) const
{
    return m_upd_parent;
}

// Sets the parent site used for neighbour updates.
void ARSSC_Process::SetParentSite(SubModels::ARSSC_Model::SiteType parent)
{
    m_upd_parent = parent;
}


// Enables use of parent site for neighbour updates, as opposed
// to custom neighbour weights defined for this reaction.
void ARSSC_Process::EnableParentWts(void) {m_use_parent_wts = true;}


// NEIGHBOUR UDPATES.

// Returns the number of sites to be incremented.  Negative
// values indicate that sites are to be decremented.
int ARSSC_Process::UpdateCount(void) const {return m_upd_count;}

// Sets the number of sites to be incremented.  Negative
// values indicate that sites are to be decremented.
void ARSSC_Process::SetUpdateCount(int n) {m_upd_count = n;}


// CUSTOM NEIGHBOUR WEIGHTS.

// Returns the custom neighbour weight for the given site.
real ARSSC_Process::NeighWt(SubModels::ARSSC_Model::SiteType site) const
{
    if ((int)site < 4) {
        return m_neigh_wts[(fvector::size_type)site];
    } else {
        return 0.0;
    }
}

// Sets the custom neighbour weight of the given site.
void ARSSC_Process::SetNeighWt(SubModels::ARSSC_Model::SiteType site, real wt)
{
    if ((int)site < 4) {
        m_neigh_wts[(fvector::size_type)site] = wt;
    }
}

// Enable custom neighbour weights, as opposed to those
// specified for the parent site.
void ARSSC_Process::EnableCustomWts(void) {m_use_parent_wts = false;}

    
// READ/WRITE/COPY.

// Writes the object to a binary stream.
void ARSSC_Process::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;

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

        // Write update site count.
        int i = (unsigned int)m_upd_count;
        out.write((char*)&i, sizeof(i));

        // Write parent site type.
        n = (unsigned int)m_upd_parent;
        out.write((char*)&n, sizeof(n));

        // Write if parent weights are used or not.
        if (m_use_parent_wts) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write the custom neighbour weight count.
        n = (unsigned int)m_neigh_wts.size();
        out.write((char*)&n, sizeof(n));

        // Write the custom neighbour weights.
        for (unsigned int i=0; i!=n; ++i) {
            val = (double)m_neigh_wts[i];
            out.write((char*)&val, sizeof(val));
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ARSSC_Reaction::Serialize).");
    }
}

// Reads the object from a binary stream.
void ARSSC_Process::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n=0;
        int i=0;
        double val=0.0;

        switch (version) {
            case 0:
                // Read number of sites.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_sites.resize(n, 0.0);

                // Read sites.
                for (unsigned int j=0; j!=n; ++j) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_sites[j] = (real)val;
                }

                // Read update site count.
                in.read(reinterpret_cast<char*>(&i), sizeof(i));
                m_upd_count = i;

                // Read parent site type.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_upd_parent = (SubModels::ARSSC_Model::SiteType)n;

                // Read if parent weights are used or not.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_use_parent_wts = (n==1);

                // Read custom neighbour weight count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_neigh_wts.resize(n, 0.0);

                // Read custom neighbour weights.
                for (unsigned int j=0; j!=n; ++j) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_neigh_wts[j] = (real)val;
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ARSSC_Reaction::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ARSSC_Reaction::Deserialize).");
    }
}
