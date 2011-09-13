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

    SilicaPrimary(const SilicaPrimary &copy); // Copy constructor.
    SilicaPrimary(                            // Stream-reading constructor.
        std::istream &in,                     //  - Input stream.
        const Sweep::ParticleModel &model     //  - Defining particle model.
        );


    SilicaPrimary(real time, const Sweep::ParticleModel &model, bool nosilica);

    // Destructors.
    virtual ~SilicaPrimary(void);

    // Sets all the properties of two primaries to be equal by using the = operator
	SilicaPrimary &operator=(const Primary &rhs);

    //! Returns a copy of the primary.
    virtual SilicaPrimary *const Clone(void) const;

    //! coagulates this particle with rhs
    SilicaPrimary &Coagulate(const Primary &rhs, int (*rand_int)(int, int),
                          real(*rand_u01)());

    //! prints the tree to a file that can be converted to a graph using graphviz
    void PrintTree(std::string filename);

    //! updates the particle cache using the particle details
    void UpdateCache(void);

    //! Sinter particles
    virtual void Sinter(
            real dt, // Delta-t for sintering to occur.
            Cell &sys, // System which defines primary's environment.
            const Processes::SinteringModel &model, // Sintering model to use.
            real (*rand_u01)() // Uniform [0,1] sample generator
            );

	// Updates Sintering level
	double SinteringLevel();

	// Adjusts the number of primaries for a surface reaction
	unsigned int Adjust(
			const fvector &dcomp,
			const fvector &dvalues,
			unsigned int n,
			real (*rand_u01)()
			);

	// Adjusts the number of primaries for interparticle reaction
	unsigned int AdjustIntPar(
			const fvector
			&dcomp,
			const fvector &dvalues,
			unsigned int n,
			real (*rand_u01)()
			);

	// Gets the number of active sites for interparticle reaction
	int GetSites() const;

	// Gets the number of active sites for interparticle reaction
	real GetSintRate() const;

    //! returns the left child
    const SilicaPrimary *LeftChild() const;

    //! returns the right child
    const SilicaPrimary *RightChild() const;

    //! Checks if the sintering level is higher then the treshold and merges the primaries if necessary
    bool CheckSintering();

    //! Updates the fractal dimension
    void CalcFractalDimension();

    AggModels::SilicaCache *const CreateAggCache() const;

    //serialize
    void Deserialize(std::istream &in, const Sweep::ParticleModel &model);
    void Serialize(std::ostream &out) const;

    AggModels::AggModelType AggID(void) const;

    //! returns L divided by W
    double LdivW() const;
    //! sum of the diameter of the primaries under this treenode needed for stats
    double PrimaryDiam() const;
    //! returns the fractal dimension
    double Fdim() const;
    //! returns the radius of gyration
    double Rg() const;
    //! returns the diameter of the largest silica
    //double silicaCollDiameter() const;
    //! returns the number of primary particles
    int Numprimary() const;
    //! returns the number of silicon atoms in the particle
    int NumSi() const;
	//! returns the number of oxygen atoms in the particle
    int NumO() const;
	//! returns the number of hydroxyl units in the particle
    int NumOH() const;
    //! returns sqrt(L*W)
    double sqrtLW() const;
	//! returns average coalescence level
    double AvgSinter() const;



protected:
    //! Empty primary not meaningful
    SilicaPrimary();

    //! help function for printree
    void PrintTreeLoop(std::ostream &out);
    //! sets the children properties to 0
    void ResetChildrenProperties();
    //! updates the particle
    void UpdateCache(SilicaPrimary *root);
    //! help function
    SilicaPrimary *SelectRandomSubparticleLoop(int target);
    //! sets the pointers to the primary particles correct after a copy event
    void UpdateAllPointers( const SilicaPrimary *source);
	//! Get All unique parents from source
    void GetAllParents(SilicaPrimary *source);
	//! Get All unique parents from source
    void FindAllParents(SilicaPrimary *source);
	//Delete target parent from m_allparents
	void DeleteParent(SilicaPrimary *target);
	//Adds source parent to m_allparents
	void AddParent(SilicaPrimary *source);
    //! updates the properties of a primary only, not the entire tree
    void UpdatePrimary(void);
    //! sets some properties to 0
    void Reset();
    //! merges the two children primaries together
    SilicaPrimary &Merge();
    //! updates the pointers after a merge event
    void ChangePointer(SilicaPrimary *source, SilicaPrimary *target);
    //! copies the node without the children
    void CopyParts( const SilicaPrimary *source);
    //! copies the subtree of a node
    void CopyTree( const SilicaPrimary *source);
    //! returns a uniformly chosen primary particle
    SilicaPrimary *SelectRandomSubparticle(Sweep::real(*rand_u01)());
	void ResetVol();
    void ReleaseMem();

	//time the two subparticles are connected
	real m_connect_time;





private:

    //! Find the path through the tree from node top to node bottom
    static std::stack<bool> recordPath(const SilicaPrimary* bottom,
                                       const SilicaPrimary* const top);

    //! Follow a path down the tree
    static SilicaPrimary* descendPath(SilicaPrimary *here,
                                   std::stack<bool> &takeLeftBranch);

    //some basic properties
    //derived from the silicas by UpdateCache()
    int m_numSi;
	int m_numO;
	int m_numOH;
    //double m_silicamass;
    //double m_silicaCollDiameter;


    //! Number of primaries below this node
    int m_numprimary;

    //sum of the diameter of the primaries under this treenode needed for stats
    double m_primarydiam;

    //properties needed to calculate the coalesence level
    //of the two primaries connected by this node
    double m_children_radius;

    double m_children_vol;

    //double m_leftparticle_vol_old;

    //double m_rightparticle_vol_old;

    double m_children_surf;

    double m_children_sintering;


    // radius of gyration and fractal dimension
    // the values are only updated in CalcFractaldimension()
    double m_Rg;
    double m_fdim;
    double m_sqrtLW;
    double m_LdivW;
    double m_avg_sinter;
	real m_sint_rate;


    SilicaPrimary *m_leftchild, *m_rightchild, *m_parent, *m_leftparticle, *m_rightparticle;
	std::vector<SilicaPrimary*> m_allparents;

};
} //namespace AggModels
} //namespace Sweep
#endif
