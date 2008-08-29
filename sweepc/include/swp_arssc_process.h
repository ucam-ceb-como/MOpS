/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of the interface for ARS-SC model enabled processes.
    Includes additional data storage for sites, which is common to
    all processes under the ARS-SC model.

    Some of the site information (in particular update info) is not
    required by all processes; the inception process does not require
    knowledge of site updates, for example.  For these processes, this
    additional data is ignored.

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

#ifndef SWEEP_ARSSC_PROCESS_H
#define SWEEP_ARSSC_PROCESS_H

#include "swp_params.h"
#include "swp_arssc_model.h"
#include <iostream>

namespace Sweep
{
namespace Processes
{
class ARSSC_Process
{
public:
    // Destructor.
    ~ARSSC_Process(void);


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


    // NEIGHBOUR UDPATES.

    // Returns the number of sites to be incremented.  Negative
    // values indicate that sites are to be decremented.
    int UpdateCount(void) const;

    // Sets the number of sites to be incremented.  Negative
    // values indicate that sites are to be decremented.
    void SetUpdateCount(int n);


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
    
    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

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

    // Constructors.  These are protected, as this object is to
    // be used as an interface for other objects, and should not
    // be constructed separately.
    ARSSC_Process();                          // Default constructor.
    ARSSC_Process(const ARSSC_Process &copy); // Copy constructor.
    ARSSC_Process(std::istream &in);          // Stream-reading constructor.

    // Assignment operator is also protected for the same reason
    // as the constructors.
    ARSSC_Process &operator=(const ARSSC_Process &rhs);
};
};
};

#endif
