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
#include "swp_particle_image.h"
#include "swp_sintering_model.h"
#include "swp_aggmodel_type.h"
#include "swp_PAH.h"
#include "swp_kmc_pah_structure.h"
#include "swp_kmc_simulator.h"
#include "swp_cell.h"
#include "swp_bintree_serializer.h"
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
    // The binary tree serialiser needs full access to private attributes.
    friend class BinTreeSerializer<class PAHPrimary>;

    // The image writer needs full access to private attributes
    template <class ParticleClass>
    friend void Sweep::Imaging::ParticleImage::ConstructTreeLoop(const ParticleClass *p);

    //! Build a new primary with one molecule
    PAHPrimary(const double time, const Sweep::ParticleModel &model);

    //! Build a new primary with one molecule
    PAHPrimary(const double time, const double position,
               const Sweep::ParticleModel &model);

    PAHPrimary(const PAHPrimary &copy); // Copy constructor.
    PAHPrimary(                       // Stream-reading constructor.
        std::istream &in,                 //  - Input stream.
        const Sweep::ParticleModel &model //  - Defining particle model.
        );


    PAHPrimary(double time, const Sweep::ParticleModel &model, bool noPAH);

    // Destructors.
    virtual ~PAHPrimary(void);

    PAHPrimary &operator=(const Primary &rhs);

    //! Returns a copy of the primary.
    virtual PAHPrimary *const Clone(void) const;
	

    //! coagulates this particle with rhs
    PAHPrimary &Coagulate(const Primary &rhs, rng_type &rng);

    //! investigate the sintering of PAH cluster
    void Sinter(double dt, Cell &sys,
                const Processes::SinteringModel &model,
                rng_type &rng,
                double wt);

    //! prints the tree to a file that can be converted to a graph using graphviz
    void PrintTree(std::string filename);

    //! updates the particle cache using the particle details
    void UpdateCache(void);

    //! updates the evolution of the PAHs using the database and the current time
	void UpdatePAHs(double t, const Sweep::ParticleModel &model, Cell &sys, rng_type &rng);

    //! adds a PAH to a particle
    void AddPAH(double time, const Sweep::ParticleModel &model);

    //! returns the rounding level due to mass addition
    double CoalescenceLevel();
    //! returns the Rounding Level according the Eq 6.3 on the markus sander's thesis
    double RoundingLevel();
    //! returns the left child
    const PAHPrimary *LeftChild() const;
    //! returns the right child
    const PAHPrimary *RightChild() const;

    //! Checks if the Rounding level is higher then the treshold (0.95) and merges the primaries if necessary
    bool CheckRounding();

    //! Updates the fractal dimension
    void CalcFractalDimension();

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
    //! returns the number of hydrogen atoms in the particle
    int NumEdgeC() const;
    //! returns the number of hydrogen atoms in the particle
    int NumRings() const;
    //! returns sqrt(L*W)
    double sqrtLW() const;
    double AvgCoalesc() const;

    //! Find Xmer, and store their information in a vector
    void FindXmer(std::vector<double> &out, int m_xmer) const;
    //! Find Xmer which statisfy sort of condition, like target_num_PAH
    void FindXmer(std::vector<std::vector<double> > &out, int target_num_PAH) const;
    //! store PAH information in a vector
    void mass_PAH(std::vector<std::vector<double> > &out) const;
    //! return the mass of Xmer including C and H 
    double MassforXmer() const;
    //! store the mass of individual PAH within this soot aggregate
    void mass_PAH(std::vector<double> &out) const;
    //! output PAH information in a vector of vector and then index the info
    void OutputPAHPSL(std::vector<std::vector<double> > &out, const int index, const double density) const;
    //! set pah_structure=Null before destructor delete it
    //void ReleasePAH(Primary &rhs);
    //! find soot particle with only one Incepted molecule (A1,A2 or A4)
    int InceptedPAH() const;
    //! check whether this PAH is invalid
    bool CheckInvalidPAHs(const boost::shared_ptr<PAH> & it) const;
    //! remove invalid PAHs under this primary particle
    void RemoveInvalidPAHs();
    //! test fragmentation assumption for the mass spectra
    void Fragtest(std::vector<double> &out, const int k, std::string mode, double threshold) const;

    double ReducedMass()const;

    //! Serialise a single PAHPrimary
    void SerializePrimary(std::ostream &out) const;

    //! Deserialise a single PAHPrimary
    void DeserializePrimary(std::istream &in, const Sweep::ParticleModel &model);

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
    //! return ture if it is a false rounding, false rounding is used to merge primary particle containing only one or no PAH after the InvalidPAHs are removed.
    bool FakeRounding();
    //! merges the two children primaries together
    void Merge();
    //! updates the pointers after a merge event
    void ChangePointer(PAHPrimary *source, PAHPrimary *target);
    //! copies the node withoud the children
    void CopyParts( const PAHPrimary *source);
    //! copies the subtree of a node
    void CopyTree( const PAHPrimary *source);
    //! returns a uniformly chosen primary particle
    PAHPrimary *SelectRandomSubparticle(rng_type &rng);
    void ReleaseMem();

  //  double pow(double a, double b);

private:

    //! Find the path through the tree from node top to node bottom
    static std::stack<bool> recordPath(const PAHPrimary* bottom,
                                       const PAHPrimary* const top);

    //! Follow a path down the tree
    static PAHPrimary* descendPath(PAHPrimary *here,
                                   std::stack<bool> &takeLeftBranch);

    void outputPAHs(std::ostream &out) const;
    void inputPAHs(std::istream &in, const Sweep::ParticleModel &model,  const int PAHcount);

    //! Set the sintering time of a tree
    void SetSinteringTime(double time);

    //some basic properties
    //derived from the PAHs by UpdataCache()
    int m_numcarbon;
    int m_numH;

    //! total num of edge C in this soot particle
    int m_numOfEdgeC;
    //! total num of rings (inculding 5, 6- menber rings) in this soot particle
    int m_numOfRings;
    //! Number of PAHs below this node
    int m_numPAH;
    //! Number of primaries below this node
    int m_numprimary;

    double m_PAHmass;
    double m_PAHCollDiameter;

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
    // store the RoundingLevel
    double m_children_roundingLevel;

    // radius of gyration and fractal dimension
    // the values are only update in CalcFractaldimension()
    double m_Rg;
    double m_fdim;
    double m_sqrtLW;
    double m_LdivW;
    double m_avg_coalesc;

    //! Absolute amount of time for which particles are sintered
    double m_sint_time;

    //! pointers to construct the connectivity and the tree stuctures
    PAHPrimary *m_leftchild, *m_rightchild, *m_parent, *m_leftparticle, *m_rightparticle;

    // Vector of PAHs.
    // PAHStructure class now have proper copy constructor
    // , but it is still not worthy to copy PAH from one vector to another
    // so we will use vector<std::tr1::shared_ptr<PAH>> instead of  vector <PAH>
    // Vector of std::tr1::shared_ptr<PAH>.
    std::vector<boost::shared_ptr<PAH> > m_PAH;


};
} //namespace AggModels
} //namespace Sweep
#endif
