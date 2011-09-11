/*!
 * \file   swp_PAH_primary.h
 * \author Markus Sander
 *  Copyright (C) 2009 Markus Sander.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Defines a primary including detailed PAH information
  Author(s):      Markus Sander (ms785)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2009 Markus Sander.

  File purpose:


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

#ifndef SWEEP_PAH_PRIMARY_H
#define SWEEP_PAH_PRIMARY_H

#include "swp_params.h"
#include "swp_primary.h"
#include "swp_particle_model.h"
#include "swp_sintering_model.h"
#include "swp_aggmodel_type.h"
#include "swp_surfvol_cache.h"
#include "swp_PAH_cache.h"
#include "swp_PAH.h"
#include "swp_kmc_pah_structure.h"
#include "swp_kmc_simulator.h"
#include "swp_cell.h"
#include <boost/shared_ptr.hpp>

#include <iostream>
#include <stack>

namespace Sweep
{
namespace AggModels
{
//! PAH primary particle class
class PAHPrimary : public Primary
{
public:
    //! Build a new primary with one molecule
    PAHPrimary(const real time, const Sweep::ParticleModel &model, int (*rand_int)(int, int));

    //! Build a new primary with one molecule
    PAHPrimary(const real time, const real position,
               const Sweep::ParticleModel &model,
               int (*rand_int)(int, int));

    PAHPrimary(const PAHPrimary &copy); // Copy constructor.
    PAHPrimary(                       // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );


    PAHPrimary(real time, const Sweep::ParticleModel &model, bool noPAH);

    // Destructors.
    virtual ~PAHPrimary(void);

    PAHPrimary &operator=(const Primary &rhs);

    //! Returns a copy of the primary.
    virtual PAHPrimary *const Clone(void) const;
	

    //! coagulates this particle with rhs
    PAHPrimary &Coagulate(const Primary &rhs, int (*rand_int)(int, int),
                          real(*rand_u01)());

    //! prints the tree to a file that can be converted to a graph using graphviz
    void PrintTree(std::string filename);

    //! updates the particle cache using the particle details
    void UpdateCache(void);

    //! updates the evolution of the PAHs using the database and the current time
	void UpdatePAHs(double t, const Sweep::ParticleModel &model, Cell &sys);

    //! adds a PAH to a particle
    void AddPAH(real time, const Sweep::ParticleModel &model);

    //! returns the coalescence level
    double CoalescenceLevel();

    //! returns the left child
    const PAHPrimary *LeftChild() const;
    //! returns the right child
    const PAHPrimary *RightChild() const;

    //! Checks if the coalescence level is higher then the treshold and merges the primaries if necessary
    bool CheckCoalescence();

    //! Updates the fractal dimension
    void CalcFractalDimension();

    AggModels::PAHCache *const CreateAggCache() const;

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
    //! returns the diameter of the largest PAH
    double PAHCollDiameter() const;
    //! returns the number of primary particles
    int Numprimary() const;
    //! returns the number of carbon atoms in the particle
    int NumCarbon() const;
    //! returns the number of hydrogen atoms in the particle
    int NumHydrogen() const;
    //! returns the number of PAH in the particle
    int NumPAH() const;
    //! returns sqrt(L*W)
    double sqrtLW() const;
    double AvgCoalesc() const;

    //! Find Xmer, and store their information in a vector
    void FindXmer(std::vector<double> &out, Xmer m_xmer) const;
    //! Find Xmer which statisfy sort of condition, like target_num_PAH
    void FindXmer(std::vector<std::vector<double> > &out, int target_num_PAH) const;
    //! store PAH information in a vector
    void mass_PAH(std::vector<std::vector<double> > &out) const;
    //! return the mass of Xmer including C and H 
    double MassforXmer() const;
    //! set pah_structure=Null before destructor delete it
    //void ReleasePAH(Primary &rhs);

protected:
    //! Empty primary not meaningful
    PAHPrimary();

    //! help function for printree
    void PrintTreeLoop(std::ostream &out);
    //! sets the children properties to 0
    void ResetChildrenProperties();
    //! updates the particle
    void UpdateCache(PAHPrimary *root);
    //! help function
    PAHPrimary *SelectRandomSubparticleLoop(int target);
    //! sets the pointers to the primary particles correct after a copy event
    void UpdateAllPointers( const PAHPrimary *source);
    //! updates the properties of a primary only, not the entire tree
    void UpdatePrimary(void);
    //! sets some properties to 0
    void Reset();
    //! merges the two children primaries together
    PAHPrimary &Merge();
    //! updates the pointers after a merge event
    void ChangePointer(PAHPrimary *source, PAHPrimary *target);
    //! copies the node withoud the children
    void CopyParts( const PAHPrimary *source);
    //! copies the subtree of a node
    void CopyTree( const PAHPrimary *source);
    //! returns a uniformly chosen primary particle
    PAHPrimary *SelectRandomSubparticle(Sweep::real(*rand_u01)());
    void ReleaseMem();
	

  //  double pow(double a, double b);

private:

    //! Find the path through the tree from node top to node bottom
    static std::stack<bool> recordPath(const PAHPrimary* bottom,
                                       const PAHPrimary* const top);

    //! Follow a path down the tree
    static PAHPrimary* descendPath(PAHPrimary *here,
                                   std::stack<bool> &takeLeftBranch);


    // Vector of PAHs.
    // PAHStructure class now have proper copy constructor (under testing)
    // , but it is still not worthy to copy PAH from one vector to another
    // so we will use vector<std::tr1::shared_ptr<PAH>> instead of  vector <PAH>
    // Vector of std::tr1::shared_ptr<PAH>.
    std::vector<boost::shared_ptr<PAH> > m_PAH;
    //some basic properties
    //derived from the PAHs by UpdataCache()
    int m_numcarbon;
    int m_numH;
    double m_PAHmass;
    double m_PAHCollDiameter;
    int m_numPAH;
//! Number of primaries below this node
    int m_numprimary;

    //sum of the diameter of the primaries under this treenode needed for stats
    double m_primarydiam;

    //properties needed to calculate the coalesence level
    //of the two primaries connected by this node
    double m_children_radius;

    double m_children_vol;

    double m_leftparticle_vol_old;

    double m_rightparticle_vol_old;

    int m_rightparticle_numPAH;

    int m_leftparticle_numPAH;

    double m_children_surf;

    double m_children_coalescence;



    // radius of gyration and fractal dimension
    // the values are only update in CalcFractaldimension()
    double m_Rg;
    double m_fdim;
    double m_sqrtLW;
    double m_LdivW;
    double m_avg_coalesc;


    PAHPrimary *m_leftchild, *m_rightchild, *m_parent, *m_leftparticle, *m_rightparticle;
};
} //namespace AggModels
} //namespace Sweep
#endif
