/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of a surface reaction process which includes a term
    for active sites.  Active site densities are calculated using an
    ActiveSitesModel, which must be set before the reaction is used.

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

#ifndef SWEEP_ARSSC_RXN_H
#define SWEEP_ARSSC_RXN_H

#include "swp_params.h"
#include "swp_process_type.h"
#include "swp_actsites_reaction.h"
#include "swp_arssc_model.h"
#include "swp_cell.h"
#include "swp_particle.h"
#include <iostream>

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Processes
{
class ARSSC_Reaction : public ActSiteReaction
{
public:
    // Constructors.
    ARSSC_Reaction(const Sweep::Mechanism &mech); // Default constructor.
    ARSSC_Reaction(const ARSSC_Reaction &copy);   // Copy constructor.
    ARSSC_Reaction(                  // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    // Destructor.
    virtual ~ARSSC_Reaction(void);

    // Operators.
    ARSSC_Reaction &operator=(const ARSSC_Reaction &rhs);


    // PERFORMING THE PROCESS.

    // Performs the process on the given system.  The responsible rate term is given
    // by index.  Returns 0 on success, otherwise negative.
    virtual int Perform(
        real t,              // Time.
        Cell &sys,           // System to update.
        unsigned int iterm=0 // The process term responsible for this event.
        ) const;

    // Performs the process on a given particle in the system.  Particle
    // is given by index.  The process is performed n times.
    virtual int Perform(
        real t,        // Current time (s).
        Cell &sys,     // System to which the particle belongs.
        Particle &sp,  // Particle for which to perform process.
        unsigned int n // Number of times to perform the process.
        ) const;


    // SITE COUNTS.

    // Returns the change in the count of the given site type.
    real SiteCount(SubModels::ARSSC_Model::SiteType type) const;

    // Returns the change in the number of free-edges.
    real FreeEdgeCount(void) const;

    // Returns the change in the number of armchairs.
    real ArmchairCount(void) const;

    // Returns the change in the number of zigzags.
    real ZigzagCount(void) const;

    // Returns the change in the number of bays.
    real BayCount(void) const;

    // Returns the change in the number of R5s.
    real R5Count(void) const;

    // Sets the change in the site count of the given site type.
    void SetSiteCount(SubModels::ARSSC_Model::SiteType type, real n);


    // PARENT SITE FOR NEIGHBOUR UPDATES.

    // Returns the parent site used for neighbour updates.
    SubModels::ARSSC_Model::SiteType ParentSite(void) const;

    // Sets the parent site used for neighbour updates.
    void SetParentSite(SubModels::ARSSC_Model::SiteType parent);


    // Enables use of parent site for neighbour updates, as opposed
    // to custom neighbour weights defined for this reaction.
    void EnableParentWts(void);


    // CUSTOM NEIGHBOUR WEIGHTS.

    // Returns the custom neighbour weight for the given site.
    real NeighWt(SubModels::ARSSC_Model::SiteType site) const;

    // Sets the custom neighbour weight of the given site.
    void SetNeighWt(
        SubModels::ARSSC_Model::SiteType site, // Set of which to set weight.
        real wt // Neighbour weight.
        );

    // Enable custom neighbour weights, as opposed to those
    // specified for the parent site.
    void EnableCustomWts(void);


    // READ/WRITE/COPY.
    
    // Creates a copy of the particle process.
    virtual ARSSC_Reaction *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    virtual ProcessType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    // Aromatic site count deviations due to reaction.
    fvector m_sites; 

    // Number of neighbour sites to update after reaction.  Negative
    // value means sites are decremented.
    int m_upd_count;

    // Parent site type used for selecting neighbours with correct
    // weights.
    SubModels::ARSSC_Model::SiteType m_upd_parent;

    // Flag indicating if the parent site neighbour weights are used
    // for selecting a neighbour.  If false then the custom weights
    // defined for the reaction are used.
    bool m_use_parent_wts;

    // Neighbour weights (if not using those for the parent site).
    fvector m_neigh_wts;

    // Default constructor is protected to prevent reactions being
    // defined without knowledge of the parent mechanism.
    ARSSC_Reaction(void);

    // Adjusts a primary particle according to the rules of the reaction.
    unsigned int adjustPri(
        Sweep::Primary &pri, // Primary to adjust.
        unsigned int n=1     // Number of times to perform adjustment.
        ) const;
};
};
};

#endif
