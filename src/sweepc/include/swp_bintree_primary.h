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
 *      originally developed by Shraddha Shekar (ss663) and Markus Sander
 *      (ms785). In its initial implementation, it was a separate class
 *      to SilicaPrimary. As there was a very large amount of code duplication,
 *      these classes were merged in Oct 2012.
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
#include "swp_particle_image.h"

namespace Sweep {

namespace AggModels {

class BinTreePrimary: public Primary {
public:

    //! Build a new primary with one molecule
    BinTreePrimary(const double time, const Sweep::ParticleModel &model);

    //! Build a new primary with one molecule
    BinTreePrimary(const double time, const double position,
               const Sweep::ParticleModel &model);

    //! Stream-reading constructor
    BinTreePrimary(std::istream &in, const Sweep::ParticleModel &model);

    //! Default destructor
    virtual ~BinTreePrimary();

    // COPYING PARTICLES AND MAINTAINING TREE STRUCTURE
    //! Primary assignment operator
    virtual BinTreePrimary &operator=(const Primary &rhs);

    //! Returns a copy of the primary.
    virtual BinTreePrimary *const Clone(void) const;

    //! Copy constructor
    BinTreePrimary(const BinTreePrimary &copy);


    // GENERAL PARTICLE MODEL PROPERTIES
    //! Returns the aggregation model type
    AggModelType AggID() const {return AggModels::BinTree_ID;}

    //! Coagulates this particle with rhs
    BinTreePrimary &Coagulate(const Primary &rhs, rng_type &rng);

    //! Updates the particle cache using the particle details
    void UpdateCache();

	//! Updates the particle cache from the root node
	//  calls UpdateCache()
	void UpdateCacheRoot();

    //! Updates sintering level
    double SinteringLevel();

    //! Prints the tree to a file that can be converted to a graph using graphviz
    void PrintTree(std::string filename) const;

    //! Adjusts a particle according to a surface reaction
    unsigned int Adjust(const fvector &dcomp,
            const fvector &dvalues, rng_type &rng, unsigned int n);

	//! Adjusts a particle according to a phase transformation reaction
	unsigned int AdjustPhase(const fvector &dcomp,
        const fvector &dvalues, rng_type &rng, unsigned int n);

    //! Adjusts a particle according to an interparticle reaction
    unsigned int AdjustIntPar(const fvector &dcomp,
            const fvector &dvalues, rng_type &rng, unsigned int n);

    //! Sinters a particle for time dt
    void Sinter(double dt, Cell &sys, const Processes::SinteringModel &model,
            rng_type &rng,
            double wt);


    // GENERAL DATA ACCESS METHODS
    //! Overload of the Mobility Diameter
    double MobDiameter() const;

    //! Get the number of primaries in the particle
    int  GetNumPrimary() const {return m_numprimary;}

    //! Get the total primary diameter of the particle (dpri,1 + dpri,2..)
    double GetPrimaryDiam() const {return m_primarydiam;}

    //! Gets the average sintering level
    double GetAvgSinterLevel() const {return m_avg_sinter;}

    //! Gets the sintering rate
    double GetSintRate() const;

    //! Gets the sintering time
    double GetSintTime() const {return m_sint_time;}

    //! Gets the arithmetic standard deviation of primary diameter
    double GetPrimaryAStdDev() const;

    //! Get the geometric mean primary diameter
    double GetPrimaryGMean() const;

    //! Gets the geometric standard deviation of primary diameter
    double GetPrimaryGStdDev() const;

    //! Gets the value of one of the chemical components
    double GetComponent(std::string name) const;

    //! Sets the value of one of the chemical components
    void SetComponent(std::string name, double val);

	//!Gets the distance between centres of primary particles
	double GetDistance() const {return m_distance_centreToCentre;}

    // SERIALISATION/DESERIALISATION
    // The binary tree serialiser needs full access to private attributes.
    friend class BinTreeSerializer<class BinTreePrimary>;

    ///////////////////////////////////////////////////////////////////////////
    /// The following ParticleImage functions have to be declared as friends to
    /// be able to access the private members of the BinTreePrimary class.
    ///////////////////////////////////////////////////////////////////////////
    template <class ParticleClass>
    friend void Sweep::Imaging::ParticleImage::ConstructTreeLoop(const ParticleClass *p);

    template <class ParticleClass>
    friend void Sweep::Imaging::ParticleImage::ConstructTree(const ParticleClass *p, Sweep::rng_type &rng, const bool trackPrimarySeparation);

    template <class ParticleClass>
    friend void Sweep::Imaging::ParticleImage::CopyTree(ImgNode &node, const ParticleClass *source);

    template <class ParticleClass>
    friend void Sweep::Imaging::ParticleImage::CopyParts(ImgNode &node, const ParticleClass *source);

    template <class ParticleClass>
    friend void Sweep::Imaging::ParticleImage::UpdateAllPointers(ImgNode &node, const ParticleClass *original);

    template <class ParticleClass>
    friend std::stack<bool> Sweep::Imaging::ParticleImage::recordPath(const ParticleClass* bottom, const ParticleClass* const top);

    template <class ParticleClass>
    friend ParticleClass* Sweep::Imaging::ParticleImage::descendPath(ParticleClass *here, std::stack<bool> &takeLeftBranch);

    //! Serialise a BinTreePrimary particle
    void Serialize(std::ostream &out) const;

    //! Deserialise a BinTreePrimary particle
    void Deserialize(std::istream &in, const Sweep::ParticleModel &model);

protected:
    //! Empty primary not meaningful
    BinTreePrimary();

    //! Print the state space of the particle to stdout
    void PrintComponents() const;

    //! Sets the children properties to 0
    void ResetChildrenProperties();

    //! Updates the properties of a primary only, not the entire tree
    void UpdatePrimary(void);

    //! Copies the subtree of a node
    void CopyTree( const BinTreePrimary *source);

    //! Sets the pointers to the primary particles correct after a copy event
    void UpdateAllPointers( const BinTreePrimary *source);

    //! Sinter a node for time dt
    void SinterNode(double dt,
            Cell &sys,
            const Processes::SinteringModel &model,
            rng_type &rng,
            double wt);

    //! Checks if the sintering level, merges particles if necessary
    bool CheckSintering();

	//! Checks if condition for merger is met
	bool MergeCondition();

    //! Set the sintering time of a tree
    void SetSinteringTime(double time);

    //! Set the last LPDA time throughout the particle tree
    void SetTime(double t);


    // STATE SPACE OF PARTICLE MODEL
    // The number of components are to be contained in Primary's m_comp
    // vector. This ensures standardised access of this information.

    // Derived quantities necessary for calculating particle trees

    //! Number of primaries below this node
    int m_numprimary;

    //! Sum of the diameter of the primaries under this treenode
	// This is usually the spherical equivalent diameter (= Primary::m_diam) 
	// unless centre to centre distance tracking is turned on
    double m_primarydiam;

	//! Sum of primary free surface areas under this node
	double m_free_surf;

    //! Equivalent spherical radius of sum of childrens' volume
    double m_children_radius;

    //! Total volume of children under this node
    double m_children_vol;

    //! Common surface area between two connected children
    double m_children_surf;

    //! Distance between the centres of primary particles.
    double m_distance_centreToCentre;

    //! Sintering level of children connected by this node
    double m_children_sintering;

    //! Average sintering level of primaries under this node
    double m_avg_sinter;

    //! Sintering rate of particle
    double m_sint_rate;

    //! Absolute amount of time for which particles are sintered
    double m_sint_time;

    // TREE STRUCTURE PROPERTIES
    // The children are the next nodes in the binary tree and are used to
    // ascend/descend the tree in a standard manner.
    // The particles are the physical particles that the node is connecting,
    // thus always refer to leaf nodes.
    // The parent is the node above the present one in the tree (NULL) for the
    // root node.

    //! Left child node
    BinTreePrimary *m_leftchild;

    //! Right child node
    BinTreePrimary *m_rightchild;

    //! Parent node
    BinTreePrimary *m_parent;

    //! Left particle node (always a leaf)
    BinTreePrimary *m_leftparticle;

    //! Right particle node (always a leaf)
    BinTreePrimary *m_rightparticle;

private:
    // GENERAL PARTICLE MODEL PROPERTIES
    //! Helper function to update the particle
    void UpdateCache(BinTreePrimary *root);

    //! Update the tree structure's surface area by increment dS
    void UpdateParents(double dS);

    //! Helper function to get a list of all primary diameters
    void GetAllPrimaryDiameters(fvector &diams) const;


    // COPYING PARTICLES AND MAINTAINING TREE STRUCTURE
    //! Find the path through the tree from node top to node bottom
    static std::stack<bool> recordPath(const BinTreePrimary* bottom,
                                       const BinTreePrimary* const top);

    //! Follow a path down the tree
    static BinTreePrimary* descendPath(BinTreePrimary *here,
                                   std::stack<bool> &takeLeftBranch);

    //! Copy elements of just the node of interest
    void CopyParts(const BinTreePrimary *source);

    //! Returns a uniformly chosen primary particle
    BinTreePrimary *SelectRandomSubparticle(rng_type &rng);

    //! Helper function for SelectRandomSubparticle
    BinTreePrimary *SelectRandomSubparticleLoop(int target);

    //! Merges the two children primaries together
    BinTreePrimary &Merge();

    //! Updates the pointers after a merge event
    void ChangePointer(BinTreePrimary *source, BinTreePrimary *target);

	//! Overloaded ChangePointer for centre to centre separation tracking model
	void ChangePointer(BinTreePrimary *source, BinTreePrimary *target, double d_ij, BinTreePrimary *small_prim);

	//! function to identify neighbours and sum their contribution to surface 
	void SumNeighbours(BinTreePrimary *prim, double &sumterm);

	//function to modify the centre to centre separations and returns free surface area
	void UpdateConnectivity(BinTreePrimary *prim, double delta_r, double &sumterm);
	
	//overload of function ignore update to neighbour
	void UpdateConnectivity(BinTreePrimary *prim, double delta_r, double &sumterm, BinTreePrimary *prim_ignore);

    // PRINTING TREES
    //! Recursive loop function for print tree
    void PrintTreeLoop(std::ostream &out) const;

    //! Helper function for printing a tree node
    void PrintTreeNode(std::ostream &out) const;


    // SERIALISATION
    //! Serialise a BinTreePrimary node
    virtual void SerializePrimary(std::ostream &out, void*) const;

    //! Deserialise a BinTreePrimary node
    virtual void DeserializePrimary(std::istream &in,
            const Sweep::ParticleModel &model, void*);

};

}

}

#endif /* SWEEP_BINTREE_PRIMARY_H */
