/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ARSSC_Reaction class declared in the
    swp_arssc_reaction.h header file.

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

#include "swp_arssc_reaction.h"
#include "swp_mechanism.h"
#include "swp_process_type.h"
#include "swp_model_factory.h"
#include "swp_actsites_model.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ARSSC_Reaction::ARSSC_Reaction(const Sweep::Mechanism &mech)
: ActSiteReaction(mech), ARSSC_Process()
{
    m_name = "ARS-SC Reaction";
}

// Copy constructor.
ARSSC_Reaction::ARSSC_Reaction(const ARSSC_Reaction &copy)
{
    *this = copy;
}

// Stream-reading constructor.
ARSSC_Reaction::ARSSC_Reaction(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
ARSSC_Reaction::~ARSSC_Reaction(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
ARSSC_Reaction &ARSSC_Reaction::operator=(const ARSSC_Reaction &rhs)
{
    if (this != &rhs) {
        ActSiteReaction::operator =(rhs);
        ARSSC_Process::operator =(rhs);
    }
    return *this;
}


// PERFORMING THE PROCESS.

/*!
 * 
 *
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       local_geom  Details of local phsyical layout
 * \param[in]       iterm       Process term responsible for this event
 * \param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 * \param[in,out]   rand_u01    Pointer to function that generates U[0,1] deviates
 * \param[out]      out         Details of any particle being transported out of system
 *
 * \return      0 on success, otherwise negative.
 */
int ARSSC_Reaction::Perform(Sweep::real t, Sweep::Cell &sys, 
                            const Geometry::LocalGeometry1d& local_geom,
                            unsigned int iterm,
                            int (*rand_int)(int, int), 
                            Sweep::real(*rand_u01)(), 
                            Sweep::Transport::TransportOutflow *out) const
{
    int i = sys.Particles().Select(static_cast<ParticleCache::PropID>(m_pid), rand_int, rand_u01);

    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);

        // Update particle with deferred processes.
        if (m_mech->AnyDeferred()) {
            // Calculate majorant rate then update the particle.
            real majr = MajorantRate(t, sys, *sp);
            m_mech->UpdateParticle(*sp, sys, t);

            // Check that the particle is still valid.
            if (sp->IsValid()) {
                real truer = Rate(t, sys, *sp);

                if (!Ficticious(majr, truer)) {
                    // Choose a primary particle to update.
                    SubParticle *sub = sp->SelectLeaf(m_modelid, m_pid);
                    Primary *pri = sub->Primary();

                    // Do ARS-SC primary update.
                    unsigned int n = adjustPri(*pri);

                    // Update the sub-particle tree above this primary.
                    pri->UpdateCache();
                    sub->UpdateTree();

                    // Update the particle ensemble.
                    sys.Particles().Update(i);

                    // Apply changes to gas-phase chemistry.
                    adjustGas(sys, n);
                }
            } else {
                // If not valid then remove the particle.
                sys.Particles().Remove(i);
            }
        } else {
            // No particle update required, just perform the surface
            // reaction.

            // Choose a primary particle to update.
            SubParticle *sub = sp->SelectLeaf(m_modelid, m_pid);
            Primary *pri = sub->Primary();

            // Do ARS-SC primary update.
            unsigned int n = adjustPri(*pri);

            // Update the sub-particle tree above this primary.
            pri->UpdateCache();
            sub->UpdateTree();

            // Update the particle ensemble.
            sys.Particles().Update(i);

            // Apply changes to gas-phase chemistry.
            adjustGas(sys, n);
        }
    } else {
        // Failed to select a particle.
        return -1;
    }

    return 0;
}

// Performs the process on a given particle in the system.  Particle
// is given by index.  The process is performed n times.
int ARSSC_Reaction::Perform(real t, Cell &sys, Particle &sp, unsigned int n) const
{
    // Choose a primary particle to update.
    SubParticle *sub = sp.SelectLeaf(m_modelid, m_pid);
    Primary *pri = sub->Primary();

    // Do ARS-SC primary update.
    unsigned int m = adjustPri(*pri, n);

    // Update the sub-particle tree above this primary.
    pri->UpdateCache();
    sub->UpdateTree();

    // Apply changes to gas-phase chemistry.
    adjustGas(sys, m);

    // Return number of times process was performed.
    return m;
}

// Adjusts a primary particle according to the rules of the reaction.
unsigned int ARSSC_Reaction::adjustPri(Sweep::Primary &pri, unsigned int n) const
{
    int m = n;

    // Get the ARS-SC sub-model from primary.
    SubModels::ARSSC_Model *ars =
        dynamic_cast<SubModels::ARSSC_Model*>(pri.SubModel(SubModels::ARSSC_Model_ID));

    // Remove those sites destroyed by the reaction.  First check that
    // there are sufficient sites to remove.
    for (unsigned int j=0; j!=SubModels::ARSSC_Model::SiteTypeCount; ++j) {
        if (m_sites[j] < 0.0) {
            m = min(m, (int)(ars->SiteCount((SubModels::ARSSC_Model::SiteType)j) -
                             ((real)m * m_sites[j])));
        }
    }
    m = max(0, m); // Avoid negative.

    if (m > 0) {
        // Remove sites destroyed by reaction.
        for (unsigned int j=0; j!=SubModels::ARSSC_Model::SiteTypeCount; ++j) {
            if (m_sites[j] < 0.0) {
                ars->RemoveSites((SubModels::ARSSC_Model::SiteType)j, -(real)m*m_sites[j]);
            }
        }

        // Adjust neighbour sites.
        if (m_use_parent_wts) {
            ars->AdjustNeighbourSites(m_upd_parent, m*m_upd_count);
        } else {
            ars->AdjustNeighbourSites(m_neigh_wts, m*m_upd_count);
        }

        // Add those sites created by the reaction.
        for (unsigned int j=0; j!=SubModels::ARSSC_Model::SiteTypeCount; ++j) {
            if (m_sites[j] > 0.0) {
                ars->AddSites((SubModels::ARSSC_Model::SiteType)j, (real)m*m_sites[j]);
            }
        }

        // Perform surface reaction update.
        SurfaceReaction::adjustPri(pri, m);
    }

    // Return number of events that were performed.
    return m;
}


// READ/WRITE/COPY.

// Creates a copy of the particle process.
ARSSC_Reaction *const ARSSC_Reaction::Clone(void) const
{
    return new ARSSC_Reaction(*this);
}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType ARSSC_Reaction::ID(void) const {return ActSiteRxn_ID;}

// Writes the object to a binary stream.
void ARSSC_Reaction::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        ActSiteReaction::Serialize(out);
        ARSSC_Process::Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ARSSC_Reaction::Serialize).");
    }
}

// Reads the object from a binary stream.
void ARSSC_Reaction::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                // Deserialize base class.
                ActSiteReaction::Deserialize(in, mech);
                ARSSC_Process::Deserialize(in);

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
