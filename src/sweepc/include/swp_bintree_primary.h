/*!
 * @file    swp_bintree_primary.h
 * @author  William J Menz
 * @brief   Generic binary tree model for multicomponent particles
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      A generalised version of the binary tree models used in SilicaPrimary
 *      and PAHPrimary. This particle model uses the 'Component' interface
 *      of the Primary parent class to track different types of atoms, rather
 *      than hard-coding the atomic types.
 *
 *      The structure is very much based on SilicaPrimary and PAHPrimary,
 *      originally developed by ss663 and ms785. Ideally these classes would
 *      inherit from this class, but this remains to be implemented.
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
 *      sweepc is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public License
 *      as published by the Free Software Foundation; either version 2
 *      of the License, or (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this program; if not, write to the Free Software
 *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *      02111-1307, USA.
 *
 *   Contact:
 *      Prof Markus Kraft
 *      Dept of Chemical Engineering
 *      University of Cambridge
 *      New Museums Site
 *      Pembroke Street
 *      Cambridge
 *      CB2 3RA, UK
 *
 *      Email:       mk306@cam.ac.uk
 *      Website:     http://como.cheng.cam.ac.uk
*/

#ifndef SWEEP_BINTREE_PRIMARY_H
#define SWEEP_BINTREE_PRIMARY_H

#include "swp_primary.h"
#include "swp_particle_model.h"
#include "swp_bintree_serializer.h"

namespace Sweep {

namespace AggModels {

class BintreePrimary: public Sweep::Primary {
public:

    //! Build a new primary with one molecule
    BintreePrimary(const real time, const Sweep::ParticleModel &model);

    //! Build a new primary with one molecule
    BintreePrimary(const real time, const real position,
               const Sweep::ParticleModel &model);

    //! Stream-reading constructor
    BintreePrimary(std::istream &in, const Sweep::ParticleModel &model);

    //! Default destructor
    virtual ~BintreePrimary();

    // COPYING PARTICLES AND MAINTAINING TREE STRUCTURE
    //! Sets all the properties of two primaries to be equal
    BintreePrimary &operator=(const Primary &rhs);

    //! Returns a copy of the primary.
    virtual BintreePrimary *const Clone(void) const;

    //! Copy constructor
    BintreePrimary(const BintreePrimary &copy);


    // GENERAL PARTICLE MODEL PROPERTIES
    //! Returns the aggregation model type
    AggModelType AggID() const {return AggModels::Bintree_ID;}

    //! Coagulates this particle with rhs
    BintreePrimary &Coagulate(const Primary &rhs, rng_type &rng);

    //! Updates the particle cache using the particle details
    void UpdateCache();

    //! Updates sintering level
    real SinteringLevel();

    //! Prints the tree to a file that can be converted to a graph using graphviz
    void PrintTree(std::string filename) const;

    //! Adjusts a particle according to a surface reaction
    unsigned int Adjust(const fvector &dcomp,
            const fvector &dvalues, rng_type &rng, unsigned int n);

    //! Adjusts a particle according to an interparticle reaction
    unsigned int AdjustIntPar(const fvector &dcomp,
            const fvector &dvalues, rng_type &rng, unsigned int n);

    //! Sinters a particle for time dt
    void Sinter(real dt, Cell &sys, const Processes::SinteringModel &model,
            rng_type &rng,
            real wt);


    // GENERAL DATA ACCESS METHODS
    //! Overload of the Mobility Diameter
    real MobDiameter() const;

    //! Get the number of primaries in the particle
    int  GetNumPrimary() const {return m_numprimary;}

    //! Get the total primary diameter of the particle (dpri,1 + dpri,2..)
    real GetPrimaryDiam() const {return m_primarydiam;}

    //! Gets the average sintering level
    real GetAvgSinterLevel() const {return m_avg_sinter;}

    //! Gets the sintering rate
    real GetSintRate() const;

    //! Gets the number of active sites (always component 0)
    int GetSites() const {return m_comp[0];}

    //! Gets the sintering time
    real GetSintTime() const {return m_sint_time;}

    // SERIALISATION/DESERIALISATION
    // The binary tree serialiser needs full access to private attributes.
    friend class BinTreeSerializer<class BintreePrimary>;

    //! Serialise a BintreePrimary particle
    void Serialize(std::ostream &out) const;

    //! Deserialise a BintreePrimary particle
    void Deserialize(std::istream &in, const Sweep::ParticleModel &model);

private:
    //! Empty primary not meaningful
    BintreePrimary();

    // GENERAL PARTICLE MODEL PROPERTIES
    //! Helper function to update the particle
    void UpdateCache(BintreePrimary *root);

    //! Update the tree structure's surface area by increment dS
    void UpdateParents(real dS);


    // COPYING PARTICLES AND MAINTAINING TREE STRUCTURE
    //! Copies the subtree of a node
    void CopyTree( const BintreePrimary *source);

    //! Sets the pointers to the primary particles correct after a copy event
    void UpdateAllPointers( const BintreePrimary *source);

    //! Find the path through the tree from node top to node bottom
    static std::stack<bool> recordPath(const BintreePrimary* bottom,
                                       const BintreePrimary* const top);

    //! Follow a path down the tree
    static BintreePrimary* descendPath(BintreePrimary *here,
                                   std::stack<bool> &takeLeftBranch);

    //! Copy elements of just the node of interest
    void CopyParts(const BintreePrimary *source);

    //! Returns a uniformly chosen primary particle
    BintreePrimary *SelectRandomSubparticle(rng_type &rng);

    //! Helper function for SelectRandomSubparticle
    BintreePrimary *SelectRandomSubparticleLoop(int target);



    // SINTERING LEVEL THINGS
    //! Set the sintering time of a tree
    void SetSinteringTime(real time);

    //! Checks if the sintering level, merges particles if necessary
    bool CheckSintering();

    //! Merges the two children primaries together
    BintreePrimary &Merge();

    //! Updates the pointers after a merge event
    void ChangePointer(BintreePrimary *source, BintreePrimary *target);

    //! Sets the children properties to 0
    void ResetChildrenProperties();

    //! Updates the properties of a primary only, not the entire tree
    void UpdatePrimary(void);

    // PRINTING TREES
    //! Recursive loop function for print tree
    void PrintTreeLoop(std::ostream &out) const;

    //! Helper function for printing a tree node
    void PrintTreeNode(std::ostream &out) const;

    void PrintComponents() const;


    // SERIALISATION
    //! Serialise a BintreePrimary node
    void SerializePrimary(std::ostream &out) const;

    //! Deserialise a BintreePrimary node
    void DeserializePrimary(std::istream &in,
            const Sweep::ParticleModel &model);


    // STATE SPACE OF PARTICLE MODEL
    // The number of components are to be contained in Primary's m_comp
    // vector. This ensures standardised access of this information.

    // Derived quantities necessary for calculating particle trees

    //! Number of primaries below this node
    int m_numprimary;

    //! Sum of the diameter of the primaries under this treenode
    real m_primarydiam;

    //! Equivalent spherical radius of sum of childrens' volume
    real m_children_radius;

    //! Total volume of children under this node
    real m_children_vol;

    //! Common surface area between two connected children
    real m_children_surf;

    //! Sintering level of children connected by this node
    real m_children_sintering;

    //! Average sintering level of primaries under this node
    real m_avg_sinter;

    //! Sintering rate of particle
    real m_sint_rate;

    //! Absolute amount of time for which particles are sintered
    real m_sint_time;

    // TREE STRUCTURE PROPERTIES
    // The children are the next nodes in the binary tree and are used to
    // ascend/descend the tree in a standard manner.
    // The particles are the physical particles that the node is connecting,
    // thus always refer to leaf nodes.
    // The parent is the node above the present one in the tree (NULL) for the
    // root node.

    //! Left child node
    BintreePrimary *m_leftchild;

    //! Right child node
    BintreePrimary *m_rightchild;

    //! Parent node
    BintreePrimary *m_parent;

    //! Left particle node (always a leaf)
    BintreePrimary *m_leftparticle;

    //! Right particle node (always a leaf)
    BintreePrimary *m_rightparticle;
};

}

}

#endif /* SWEEP_BINTREE_PRIMARY_H */
