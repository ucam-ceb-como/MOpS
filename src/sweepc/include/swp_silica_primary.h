/*
  Author(s): Shraddha Shekar (ss663), Markus Sander (ms785)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2010 Shraddha Shekar.


  File purpose:
	The SilicaCache class is a specialisation of the AggModelCache
	class for holding cached data of silica primary particles.

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

#ifndef SWEEP_SILICA_PRIMARY_H
#define SWEEP_SILICA_PRIMARY_H

#include "swp_params.h"
#include "swp_primary.h"
#include "swp_particle_model.h"
#include "swp_sintering_model.h"
#include "swp_aggmodel_type.h"
#include "swp_surfvol_cache.h"
#include "swp_silica_cache.h"
#include "swp_cell.h"

#include <iostream>
#include <stack>

namespace Sweep
{
namespace AggModels
{
//! Silica Primary particle class
class SilicaPrimary : public Primary
{
public:
    //! Build a new primary with one molecule
    SilicaPrimary(const real time, const Sweep::ParticleModel &model);

    //! Build a new primary with one molecule
    SilicaPrimary(const real time, const real position,
               const Sweep::ParticleModel &model);

	//! Copy constructor
    SilicaPrimary(const SilicaPrimary &copy);
    
    //! Stream-reading constructor
    SilicaPrimary(                            // Stream-reading constructor.
        std::istream &in,                     //  - Input stream.
        const Sweep::ParticleModel &model     //  - Defining particle model.
        );


    SilicaPrimary(real time, const Sweep::ParticleModel &model, bool nosilica);

    //! Destructors.
    virtual ~SilicaPrimary(void);

    //! Sets all the properties of two primaries to be equal by using the = operator
	SilicaPrimary &operator=(const Primary &rhs);

    //! Returns a copy of the primary.
    virtual SilicaPrimary *const Clone(void) const;

    //! Sets the state space when initialising a primary from XML
    void SetStateSpace(const int numSi, const int numO, const int numOH);

    //! Coagulates this particle with rhs
    SilicaPrimary &Coagulate(const Primary &rhs, rng_type &rng);

    //! Prints the tree to a file that can be converted to a graph using graphviz
    void PrintTree(std::string filename);

    //! Updates the particle cache using the particle details
    void UpdateCache(void);

    //! Sinter particles
    virtual void Sinter(
            real dt, // Delta-t for sintering to occur.
            Cell &sys, // System which defines primary's environment.
            const Processes::SinteringModel &model, // Sintering model to use.
            rng_type &rng,  // Random number generator
            real wt     // Statistical weight
            );

	//! Updates Sintering level
	real SinteringLevel();

	//! Adjusts the number of primaries for a surface reaction
	unsigned int Adjust(
			const fvector &dcomp,
			const fvector &dvalues,
			rng_type &rng,			// Random number for leaf node
			unsigned int n
			);

	//! Adjusts the number of primaries for interparticle reaction
	unsigned int AdjustIntPar(
			const fvector &dcomp,
			const fvector &dvalues,
			rng_type &rng,
			unsigned int n
			);

	//! Gets the number of active sites for interparticle reaction
	int GetSites() const;

	//! Gets the sintering rate for interparticle reaction
	real GetSintRate() const;

	//! Gets the sintering time
	real GetSintTime() const;

    //! Updates the fractal dimension
    void CalcFractalDimension();

    AggModels::SilicaCache *const CreateAggCache() const;

    //! Deserialize
    void Deserialize(std::istream &in, const Sweep::ParticleModel &model);
    
    //! Serialize
    void Serialize(std::ostream &out) const;

    AggModels::AggModelType AggID(void) const;

    //! Returns the left child
    const SilicaPrimary *LeftChild() const;
    //! Returns the right child
    const SilicaPrimary *RightChild() const;

    // Functions used to gather data for statistics

    //! returns L divided by W
    real LdivW() const;
    //! Sum of the diameter of the primaries under this treenode needed for stats
    real PrimaryDiam() const;
    //! Returns the fractal dimension
    real Fdim() const;
    //! Returns the radius of gyration
    real Rg() const;
    //! Returns the number of primary particles
    int Numprimary() const;
    //! Returns the number of silicon atoms in the particle
    int NumSi() const;
	//! Returns the number of oxygen atoms in the particle
    int NumO() const;
	//! Returns the number of hydroxyl units in the particle
    int NumOH() const;
    //! Returns sqrt(L*W)
    real sqrtLW() const;
	//! Returns average coalescence level
    real AvgSinter() const;

private:
    //! Empty primary not meaningful
    SilicaPrimary();

    //! Checks if the sintering level is higher then the threshold and merges the primaries if necessary
    bool CheckSintering();

    //! Overload of Primary's SetTime function
    void SetTime(real t);

    //! Help function for printree
    void PrintTreeLoop(std::ostream &out);
    //! Help function for printree
    void PrintTreeNode(std::ostream &out);
    //! Sets the children properties to 0
    void ResetChildrenProperties();
    //! Updates the particle
    void UpdateCache(SilicaPrimary *root);
    //! Help function
    SilicaPrimary *SelectRandomSubparticleLoop(int target);
    //! Sets the pointers to the primary particles correct after a copy event
    void UpdateAllPointers( const SilicaPrimary *source);
    //! Updates the properties of a primary only, not the entire tree
    void UpdatePrimary(void);
    //! Sets some properties to 0
    void Reset();
    //! Merges the two children primaries together
    SilicaPrimary &Merge();

    //! Releases the memory associated with the object
    void ReleaseMem();

    //! Updates the pointers after a merge event
    void ChangePointer(SilicaPrimary *source, SilicaPrimary *target);
    //! Copies the node without the children
    void CopyParts( const SilicaPrimary *source);
    //! Copies the subtree of a node
    void CopyTree( const SilicaPrimary *source);
    //! Returns a uniformly chosen primary particle
    SilicaPrimary *SelectRandomSubparticle(rng_type &rng);
    
    //! Update the surface area and sintering level of all parents
    void UpdateParents(real dS);

    //! Find the path through the tree from node top to node bottom
    static std::stack<bool> recordPath(const SilicaPrimary* bottom,
                                       const SilicaPrimary* const top);

    //! Follow a path down the tree
    static SilicaPrimary* descendPath(SilicaPrimary *here,
                                   std::stack<bool> &takeLeftBranch);

    //! Set the sintering time of a tree
    void SetSinteringTime(real time);

    //! Number of silicon units in primary
    int m_numSi;
    
    //! Number of O units in primary
    int m_numO;
    
    //! Number of OH units in primary
    int m_numOH;

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

    /* Imaging properties
     * These are presently unused, however may be useful if 
     * one wishes to generate an image of the particle
     */
    
    //! Radius of gyration (currently unused)
    real m_Rg;
    
    //! Fractal dimension (currently unused)
    real m_fdim;
    
    //! Square-root of the length times width (currently unused)
    real m_sqrtLW;
    
    //! Length divided by width (currently unused)
    real m_LdivW;

    /*
     * Definition of the silica primaries
     * L/R child are the next nodes of the binary tree
     * L/R particle are the actual particles which *this* connects
     * Parent is the node above the current object in the bin. tree
     * 
     * Refer to Markus Sander's thesis for more detailed description
     */    
    SilicaPrimary *m_leftchild, *m_rightchild, *m_parent, *m_leftparticle, *m_rightparticle;

	//! Absolute amount of time for which particles are sintered
	real m_sint_time;

};
} //namespace AggModels
} //namespace Sweep
#endif
