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
#include "swp_coords.h"
#include <boost/shared_ptr.hpp>

#include "swp_coords.h"

#include <iostream>
#include <stack>
#include <map>
#include <set>

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

    ///////////////////////////////////////////////////////////////////////////
    /// The following ParticleImage functions have to be declared as friends to
    /// be able to access the private members of the PAHPrimary class.
    ///////////////////////////////////////////////////////////////////////////
    template <class ParticleClass>
    friend void Sweep::Imaging::ParticleImage::ConstructTree(const ParticleClass *p, Sweep::rng_type &rng, const bool trackPrimaryCoordinates);

    template <class ParticleClass>
    friend void Sweep::Imaging::ParticleImage::ConstructTreeLoop(const ParticleClass *p);

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

    //! For use while avoiding repeated deserialisation of a single PAH
    typedef std::map<void*, boost::shared_ptr<PAH> > PahDeserialisationMap;

    //! For use while avoiding repeated serialisation of a single PAH
    typedef std::set<void*> PahSerialisationMap;

    //! Build a new primary with one molecule
    PAHPrimary(const double time, const Sweep::ParticleModel &model);

    //! Build a new primary with one molecule
    PAHPrimary(const double time, const double position,
               const Sweep::ParticleModel &model);

    PAHPrimary(const PAHPrimary &copy); // Copy constructor.
    
    //! Stream-reading constructor
    PAHPrimary(                       
        std::istream &in,                        // Input stream
        const Sweep::ParticleModel &model,       // Defining particle model
        PahDeserialisationMap &pah_duplicates    // Information on duplicated PAH
        );


    PAHPrimary(double time, const Sweep::ParticleModel &model, bool noPAH);

    // Destructors.
    virtual ~PAHPrimary(void);

    PAHPrimary &operator=(const Primary &rhs);

    //! Returns a copy of the primary.
    virtual PAHPrimary *const Clone(void) const;
	
    //! coagulates this particle with rhs
    PAHPrimary &Coagulate(const Primary &rhs, rng_type &rng);

    //! coagulates this particle with rhs
    PAHPrimary &Fragment(const Primary &rhs, rng_type &rng);

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
	void UpdatePAHs(double t, double dt, const Sweep::ParticleModel &model, Cell &sys, int statweight, int ind, rng_type &rng,
		PartPtrVector &overflow);
	
	//! overload function
	//! updates the evolution of the PAHs using the database and the current time
	//! if coordinates of primary are tracked, the free surface area of the particle can be related to surface growth
	void UpdatePAHs(double t, double dt, const Sweep::ParticleModel &model, Cell &sys, int statweight, int ind, rng_type &rng,
		PartPtrVector &overflow, double fs);

	//! adjust the primary after surface growth in PAH_KMC model.
	void Adjust(const double old_vol);

    //! adds a PAH to a particle
    void AddPAH(double time, const Sweep::ParticleModel &model);

    //! returns the rounding level due to mass addition
    double CoalescenceLevel();

    //! Returns the distance between the centres of primary particles.
    double Distance() const;

	//! Calculates the radius of gyration.
	double GetRadiusOfGyration() const;

	//! Returns a vector of primary coordinates and radius (4D).
	void GetPriCoords(std::vector<fvector> &coords) const;

	//! Returns the free surface area of a particle.
	double GetFreeSurfArea() const;

    //! returns the Rounding Level according the Eq 6.3 on the markus sander's thesis
    double RoundingLevel();
    //! returns the left child
    const PAHPrimary *LeftChild() const;
    //! returns the right child
    const PAHPrimary *RightChild() const;

    //! Checks if the Rounding level is higher then the treshold (0.95) and merges the primaries if necessary
	//! Used if primary coordinates are not tracked.
    bool CheckRounding();

	//! Checks if the sintering level, merges particles if necessary
	//! Used if primary coordinates are tracked.
	bool CheckSintering();

    //! Updates the fractal dimension
    void CalcFractalDimension();

    //! Deserialize object from input binary stream
    void Deserialize(std::istream &in, const Sweep::ParticleModel &model, PahDeserialisationMap &pah_duplicates);

	//! Return primary particle details and connectivity
	void PrintPrimary(std::vector<fvector> &surface, std::vector<fvector> &primary_diameter, int k) const;
    
    //! Serialize object to output binary stream
    void Serialize(std::ostream &out, void *duplicates) const;

    //! Serialise a single PAHPrimary
    void SerializePrimary(std::ostream &out, void *duplicates) const;

    //! Deserialise a single PAHPrimary
    void DeserializePrimary(std::istream &in, const Sweep::ParticleModel &model, void *duplicates);

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
    //! returns the number of 6-member rings in the particle.
    int NumRings() const;
	//! returns the number of 5-member rings in the particle.
	int NumRings5() const;
    //! returns sqrt(L*W)
    double sqrtLW() const;

    double AvgCoalesc() const;
	double AvgSinter() const;

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
    void OutputPAHPSL(std::vector<std::vector<double> > &out, const int index, const double density, std::set<void*> &pah_duplicates, std::vector<std::string> &Mapping, const double timeStep) const;
	//! output primary particle information in a vector of vector and then index the info
	void OutputPPPSL(std::vector<std::vector<double> > &out, const int index, const double density,const double timeStep) const;
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

	std::vector<boost::shared_ptr<PAH> > GetPAHVector() const;

	//! Updates sintering level
	//! used if primary coordinates are tracked
	double SinteringLevel();

	//! Returns the cap volume of two connected primary particles
	double CalcChildrenSumCap();


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
	//void Merge();
    //! updates the pointers after a merge event
    void ChangePointer(PAHPrimary *source, PAHPrimary *target);

    //! copies the node withoud the children
    void CopyParts( const PAHPrimary *source);
    //! copies the subtree of a node
    void CopyTree( const PAHPrimary *source);
    //! returns a uniformly chosen primary particle
    PAHPrimary *SelectRandomSubparticle(rng_type &rng);
    void ReleaseMem();

	//! Distance between the centres of primary particles.
	double m_distance_centreToCentre;

	//! For tracking the coordinates of primary particles.
	Coords::Vector m_cen_bsph; //!< Bounding-sphere centre.
	Coords::Vector m_cen_mass; //!< Centre-of-mass coordinates.

	//! Radius of bounding sphere raised to powers of 1, 2 and 3.
	double m_r;  //!< Bounding sphere radius of aggregate/primary.
	double m_r2; //!< r squared (useful for efficient collision detection computation).
	double m_r3; //!< r cubed (useful for calculating centre-of-mass).

	//! if primary coordinates are tracked
	//! Sintering level of children connected by this node
	double m_children_sintering;

	///////////////////////////////////////////////////////////////////////////
	/// Functions for manipulating coordinates of primary particles in an
	/// aggregate.
	///////////////////////////////////////////////////////////////////////////

	//! Returns the bounding-sphere centre.
	const Coords::Vector &boundSphCentre(void) const;

	//! Calculates the bounding sphere position and radius using the left and
	//! right child node values.
	void calcBoundSph(void);

	//! Calculates the centre-of-mass using the left and right child node
	//! values.
	void calcCOM(void);

	//! Put the bounding-sphere at the origin.
	void centreBoundSph(void);

	//! Put the centre-of-mass at the origin.
	void centreCOM(void);

	//! Returns true if this node is a leaf (has no children).
	bool isLeaf(void) const;

	//! Returns the bounding sphere radius.
	double Radius(void) const;

	//! Rotates the aggregate node and child structure about its centre-of-mass
	//! by the given angles (spherical coordinates).
	void rotateCOM(double theta, fvector V);

	//! Translates (moves) the aggregate node and child structure by the given
	//! amounts along the cartesian axes.
	void Translate(double dx, double dy, double dz);

	//! Write the coordinates of the primaries belonging to the node pointed to
	//! by the this pointer.
	void writePrimaryCoordinatesRadius(void);

	//! From bintree model.
	//! Check for the overlap of primary particles.
	bool checkForOverlap(
		PAHPrimary &target, //!< Target node.
		PAHPrimary &bullet, //!< Bullet node.
		int &numberOfOverlaps,  //!< Number of overlaps.
		double &Separation      //!< Separation between the centres of the primary particles for use with the Newton bisection method.   
		);

	//! Determine whether the particles overlap.
	static bool particlesOverlap(
		const Coords::Vector &p1, //!< Positional vector of sphere 1.
		double r1,                //!< Radius of sphere 1.
		const Coords::Vector &p2, //!< Positional vector of sphere 2.
		double r2,                //!< Radius of sphere 2.
		double &Separation        //!< Separation between the centres of the primary particles for use with the Newton bisection method.   
		);

	//! Sets the radius of the bounding sphere.
	void setRadius(double r);

	//! Transforms the node coordinates using the given transformation matrix.
	void transform(const Coords::Matrix &mat);

	//! Sum of primary free surface areas under this node
	double m_free_surf;

	//! From bintree model.
	//! Sinter a node for time dt
	void SinterNode(double dt,
		Cell &sys,
		const Processes::SinteringModel &model,
		rng_type &rng,
		double wt);

	//! Analogous to bintree model.
	//! Checks if condition for merger is met
	bool MergeCondition();

	//! Sintering rate of particle
	double m_sint_rate;

	double m_avg_sinter;

	//primary volume -- different to m_vol if centre to centre seapration is tracked
	double m_primaryvol;

	//csl37-rewrite sum of neck radii * ri/xij
	//necessary for sintering process
	//cached
	double m_sum_necks;

	//! From bintree model
	//! Function to return the separation unit vector between two coordinates
	Coords::Vector UnitVector(Coords::Vector x_i, Coords::Vector x_j); 

	//! From bintree model
	//! Translates a primary particle
	void TranslatePrimary(Coords::Vector u, double delta_d);

	//! From bintree model
	//! Function to translate neighbours of a primary except prim_ignore
	void TranslateNeighbours(PAHPrimary *prim, Coords::Vector u, double delta_d, PAHPrimary *prim_ignore);

	//! From bintree model
	//! Calculates distance between two points
	double Separation(Coords::Vector x_i, Coords::Vector x_j);

  //  double pow(double a, double b);

private:

    //! Find the path through the tree from node top to node bottom
    static std::stack<bool> recordPath(const PAHPrimary* bottom,
                                       const PAHPrimary* const top);

    //! Follow a path down the tree
    static PAHPrimary* descendPath(PAHPrimary *here,
                                   std::stack<bool> &takeLeftBranch);

    //! Writes individual PAHs to a binary stream 
    void outputPAHs(std::ostream &out, PahSerialisationMap &pah_duplicates) const;

    //! Reads individual PAHs from the binary stream 
    void inputPAHs(std::istream &in, const Sweep::ParticleModel &model,  const int PAHcount, PahDeserialisationMap &pah_duplicates);

    //! Set the sintering time of a tree
    void SetSinteringTime(double time);

	////! From bintree model.
	////! function to identify neighbours and sum their contribution to surface
	//void SumNeighbours(PAHPrimary *prim, double &sumterm);

    //some basic properties
    //derived from the PAHs by UpdataCache()
    int m_numcarbon;
    int m_numH;

    //! total num of edge C in this soot particle
    int m_numOfEdgeC;
    //! Total number of 6-member rings in this soot particle
    int m_numOfRings;
	//! Total number of 6-member rings in this soot particle
	int m_numOfRings5;
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

	//used to calculate geometric volume of a particle
	double m_children_sumCap;
	double m_sum_cap;
	double m_sph_prim_vol;

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

	//! From bintree model
	//! Update primary free surface area and volume
	void UpdateOverlappingPrimary();

	//! From bintree model
	//! function to identify neighbours and sum their cap areas and volumes
	void SumCaps(PAHPrimary *prim, double &CapAreas, double &CapVolumes, double &SumNecks);

	//! From bintree model
	//! function to modify the centre to centre separations and coordinates and neighbours
	void UpdateConnectivity(PAHPrimary *prim, double delta_r, PAHPrimary *prim_ignore);

	//! From bintree model
	//Function to adjust primary properties
	void AdjustPrimary(double dV, double d_ij, PAHPrimary *prim_ignore);

	//! From bintree model
	//! Overloaded ChangePointer for centre to centre separation and coordinate tracking models
	void ChangePointer(PAHPrimary *source, PAHPrimary *target, PAHPrimary *small_prim, PAHPrimary *node);

	//! From bintree model
	//! Add new neighbours during a merger event
	double AddNeighbour(double A_n_k, PAHPrimary *small_prim, PAHPrimary *node);

	//! analogous to bintree model
	//! Merges the two children primaries together
	PAHPrimary &Merge();
	
};
} //namespace AggModels
} //namespace Sweep
#endif
