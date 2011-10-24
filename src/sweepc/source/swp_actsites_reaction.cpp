/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ActSiteReaction class declared in the
    swp_actsites_reaction.h header file.

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

#include "swp_actsites_reaction.h"
#include "swp_mechanism.h"
#include "swp_process_type.h"
#include "swp_model_factory.h"

#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor is protected to prevent reactions being
// defined without knowledge of the parent mechanism.
ActSiteReaction::ActSiteReaction(void)
: SurfaceReaction(), m_asmodel(NULL)
{
    m_name = "Active-site Reaction";
}

// Default constructor.
ActSiteReaction::ActSiteReaction(const Sweep::Mechanism &mech)
: SurfaceReaction(mech), m_asmodel(NULL)
{
    m_name = "Active-site Reaction";
}

// Copy constructor.
ActSiteReaction::ActSiteReaction(const ActSiteReaction &copy)
{
    *this = copy;
}

// Stream-reading constructor.
ActSiteReaction::ActSiteReaction(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
ActSiteReaction::~ActSiteReaction(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
ActSiteReaction &ActSiteReaction::operator =(const ActSiteReaction &rhs)
{
    if (this != &rhs) {
        SurfaceReaction::operator =(rhs);
        m_asmodel = rhs.m_asmodel;
    }
    return *this;
}


// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

/*!
 *@param[in]            t           Time at which rate is being calculated
 *@param[in]            sys         System for which rate is to be calculated
 *@param[in]            local_geom  Spatial configuration information (ignored)
 *
 *@return   Process rate
 */
real ActSiteReaction::Rate(real t, const Cell &sys,
                           const Geometry::LocalGeometry1d &local_geom) const
{
    return SurfaceReaction::Rate(t, sys, local_geom) *
           m_asmodel->SiteDensity(t, sys.GasPhase(), sys.Particles());
}


// SINGLE PARTICLE RATE CALCULATIONS.

// Returns the rate of the process for the given particle in
// the system. Process must be linear in particle number.
real ActSiteReaction::Rate(real t, const Cell &sys, const Particle &sp) const
{
    return SurfaceReaction::Rate(t, sys, sp) *
           m_asmodel->SiteDensity(t, sys.GasPhase(), sp);
}


// ACTIVE SITES MODEL.

// Sets the active sites model.
void ActSiteReaction::SetModel(ActSites::ABFModel &model)
{
    m_asmodel = &model;
}

// Returns the active sites model.
ActSites::ABFModel *const ActSiteReaction::Model(void) const
{
    return m_asmodel;
}


// READ/WRITE/COPY.

// Creates a copy of the particle process.
ActSiteReaction *const ActSiteReaction::Clone(void) const
{
    return new ActSiteReaction(*this);
}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType ActSiteReaction::ID(void) const {return ActSiteRxn_ID;}

// Writes the object to a binary stream.
void ActSiteReaction::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        SurfaceReaction::Serialize(out);

        // Output the active sites model ID.
        if (m_asmodel != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            unsigned int n = (unsigned int)m_asmodel->ID();
            out.write((char*)&n, sizeof(n));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ActSiteReaction::Serialize).");
    }
}

// Reads the object from a binary stream.
void ActSiteReaction::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0, id = 0;

        switch (version) {
            case 0:
                // Deserialize base class.
                SurfaceReaction::Deserialize(in, mech);

                // Read the active sites model ID.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n==1) {
                    in.read(reinterpret_cast<char*>(&id), sizeof(id));
                    m_asmodel = ModelFactory::GetActSitesModel((ActSites::ActSitesType)id);
                } else {
                    m_asmodel = NULL;
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ActSiteReaction::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ActSiteReaction::Deserialize).");
    }
}
