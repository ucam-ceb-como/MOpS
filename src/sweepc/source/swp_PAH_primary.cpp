/*!
 * \file   swp_PAH_primary.cpp
 * \author Markus Sander
 *
 * \brief  Particle Model that stores the individual PAHs of a soot particle
 */
/*
  Author(s):      Markus Sander (ms785)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Markus Sander.

  File purpose:
    Implementation of the PAHPrimary class declared in the
    swp_PAH_primary.h header file.

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

#include "swp_primary.h"
#include "swp_PAH_primary.h"
#include "swp_aggmodel_type.h"
#include "swp_model_factory.h"
#include "swp_particle_image.h"
#include "swp_cell.h"
#include "swp_kmc_pah_process.h"
#include "swp_kmc_pah_structure.h"
#include "swp_PAH.h"
#include "swp_ensemble.h"
#include "swp_particle_model.h"

#include <stdexcept>
#include <cassert>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/bind.hpp>
#include <boost/bind/placeholders.hpp>

#include "string_functions.h"

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace Sweep::KMC_ARS;

using namespace std;
using namespace Strings;

//used for debugging, testing clone function for PAHStructure.
static unsigned int ID=0; 
static bool m_clone=false;
/*
double PAHPrimary::pow(double a, double b) {
    int tmp = (*(1 + (int *)&a));
    int tmp2 = (int)(b * (tmp - 1072632447) + 1072632447);
    double p = 0.0;
    *(1 + (int * )&p) = tmp2;
    return p;
}
*/

// CONSTRUCTORS AND DESTRUCTORS.
PAHPrimary::PAHPrimary() : Primary(),
    m_numcarbon(0),
    m_numH(0),
    m_numOfEdgeC(0),
    m_numOfRings(0),
    m_numPAH(0),
    m_numprimary(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_roundingLevel(0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL)
{
}

/*!
 * @param[in]       time        Time at which particle is being created
 * @param[in]       model       Model which defines the meaning of the primary
 *
 */
PAHPrimary::PAHPrimary(const double time, const Sweep::ParticleModel &model)
: Primary(time, model),
    m_numcarbon(0),
    m_numH(0),
    m_numOfEdgeC(0),
    m_numOfRings(0),
    m_numPAH(0),
    m_numprimary(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_roundingLevel(0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL)
{
    // Other parts of the code check for a non-zero composition
    m_comp[0]=1;

    AddPAH(time, model);

    //Update the other properties
    UpdateCache();
}

/*!
 * @param[in]       time        Time at which particle is being created
 * @param[in]       position    Position at which particle is being created
 * @param[in]       model       Model which defines the meaning of the primary
 *
 */
PAHPrimary::PAHPrimary(const double time, const double position,
                       const Sweep::ParticleModel &model)
: Primary(time, model),
    m_numcarbon(0),
    m_numH(0),
    m_numOfEdgeC(0),
    m_numOfRings(0),
    m_numPAH(0),
    m_numprimary(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_roundingLevel(0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL)
{
    // Other parts of the code check for a non-zero composition
    m_comp[0]=1;
   
	AddPAH(time, model);

    //Update the other properties
    UpdateCache();
}


// Initialising constructor.
PAHPrimary::PAHPrimary(double time, const Sweep::ParticleModel &model, bool noPAH)
: Primary(time, model),
    m_numcarbon(0),
    m_numH(0),
    m_numOfEdgeC(0),
    m_numOfRings(0),
    m_numPAH(0),
    m_numprimary(0),
    m_PAHmass(0),
    m_PAHCollDiameter(0),
    m_primarydiam(0.0),
    m_children_radius(0),
    m_children_vol(0),
    m_leftparticle_vol_old(0),
    m_rightparticle_vol_old(0),
    m_rightparticle_numPAH(0),
    m_leftparticle_numPAH(0),
    m_children_surf(0),
    m_children_roundingLevel(0),
    m_Rg(0),
    m_fdim(0),
    m_sqrtLW(0),
    m_LdivW(0),
    m_avg_coalesc(0),
    m_sint_time(0.0),
    m_leftchild(NULL),
    m_rightchild(NULL),
    m_parent(NULL),
    m_leftparticle(NULL),
    m_rightparticle(NULL)
{
    m_comp[0]=1;
}


/*double PAHPrimary::Fdim() const
{
    return m_fdim;
}*/


/*!
 * Add a PAH to the primary particle
 *
 * @param[in]   time        create time of the PAH

 * @param[in]   model       Particle model containing molecule database
*/
void PAHPrimary::AddPAH(double time,const Sweep::ParticleModel &model)
{
    boost::shared_ptr<PAH> new_PAH (new PAH(time, model.InceptedPAH()));
    new_PAH->PAH_ID=ID;
    m_PAH.push_back(new_PAH);
    ID++;
    // Set the particle mass, diameter etc
    UpdatePrimary();
}

// Copy constructor.
PAHPrimary::PAHPrimary(const PAHPrimary &copy)
{
    *this = copy;
    if (copy.m_leftchild!=NULL)
    {
        CopyTree(&copy);
    }
	//m_clone=false;
}


// Default destructor.
PAHPrimary::~PAHPrimary()
{
    delete m_leftchild;
    delete m_rightchild;
    // it is not necessary to delete m_leftparticle because
    // this is also m_leftchild somewhere down the tree
    if (m_PAH.size()!=0) m_PAH.clear();
    // delete the PAH list
    releaseMem();
    m_clone=false;
}

/*!
 * Recursively copy the tree for non-leaf nodes.
 *
 * @param[in] source Pointer to the primary to be copied
*/
void PAHPrimary::CopyTree( const PAHPrimary *source)
{
    //create the new left and right children with nothing in them
    m_leftchild = new PAHPrimary(source->CreateTime(),*m_pmodel,false);
    m_rightchild = new PAHPrimary(source->CreateTime(),*m_pmodel,false);

    // copy the properties such as the volume, surface area and list of
	// constituent PAH molecules
	m_leftchild->CopyParts(source->m_leftchild);
	m_rightchild->CopyParts(source->m_rightchild);

    //set the pointers to the parent
	m_leftchild->m_parent=this;
	m_rightchild->m_parent=this;

    // the left and right particle are set further down in UpdateAllPointers
	// These are the pointers that specify which primary particles touch each
	// other in the aggregat structure.
	m_leftparticle=NULL;
	m_rightparticle=NULL;

    //recurse to copy the subtrees
	if (source->m_leftchild->m_leftchild!=NULL)
        m_leftchild->CopyTree(source->m_leftchild);

	if (source->m_rightchild->m_leftchild!=NULL)
        m_rightchild->CopyTree(source->m_rightchild);

    //set the leftparticle and rightparticle
	UpdateAllPointers(source);
}

PAHPrimary &PAHPrimary::operator=(const Primary &rhs)
{
    if (this != &rhs) {
        const AggModels::PAHPrimary *pahprimary = NULL;
        pahprimary = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);
        UpdateAllPointers(pahprimary);
    }
    return *this;
}

/*
 * @brief Stream-reading constructor
 *
 * @param[in,out]	 in		       Input binary stream  
 * @param[in]        model	       Particle model defining interpretation of particle data
 * @param[in,out]    duplicates    Addresses of PAHs for use when reading primary particles
 */ 
PAHPrimary::PAHPrimary(std::istream &in, const Sweep::ParticleModel &model, PahDeserialisationMap &pah_duplicates)
{
    Deserialize(in, model, pah_duplicates);
}


/*!
 * This function is like a limited assignment operator, except that the
 * children are not copied and the pointers to the particles may need
 * adjusting after this method has finished.
 *
 * @param[in] source Pointer to the primary to be copied
 */
void PAHPrimary::CopyParts( const PAHPrimary *source)
{
    SetCollDiameter(source->CollDiameter());
    SetMobDiameter(source->MobDiameter());
    SetSphDiameter(source->SphDiameter());
    m_PAHCollDiameter=source->m_PAHCollDiameter;
    SetSurfaceArea(source->SurfaceArea());
    m_time=source->m_time;
    m_PAHmass=source->m_PAHmass;
    m_leftchild=source->m_leftchild;
    m_rightchild=source->m_rightchild;
    m_leftparticle=source->m_leftparticle;
    m_rightparticle=source->m_rightparticle;
    m_parent=source->m_parent;
    SetMass(source->Mass());
    m_numPAH=source->m_numPAH;
    m_numprimary=source->m_numprimary;
    m_primarydiam=source->m_primarydiam;
    m_sqrtLW=source->m_sqrtLW;
    m_LdivW=source->m_LdivW;
    m_pmodel=source->m_pmodel;
    m_surf=source->m_surf;
    m_vol=source->m_vol;
    m_children_surf=source->m_children_surf;
    m_children_vol=source->m_children_vol;
    m_children_radius=source->m_children_radius;
    m_children_roundingLevel=source->m_children_roundingLevel;
    m_rightparticle_numPAH=source->m_rightparticle_numPAH;
    m_leftparticle_numPAH=source->m_leftparticle_numPAH;
    m_leftparticle_vol_old=source->m_leftparticle_vol_old;
    m_rightparticle_vol_old=source->m_rightparticle_vol_old;
    m_fdim=source->m_fdim;
    m_Rg=source->m_Rg;
    m_avg_coalesc=source->m_avg_coalesc;
    m_numcarbon = source->m_numcarbon;
    m_numH = source->m_numH;
    m_numOfEdgeC = source->m_numOfEdgeC;
    m_numOfRings = source->m_numOfRings;
    m_values=source->m_values;
    m_comp=source->m_comp;
    m_sint_time=source->m_sint_time;

    // Replace the PAHs with those from the source
    if (m_clone==true){
		    m_PAH.resize(source->m_PAH.size());
            m_PAH.clear();
    for (size_t i=0; i!=source->m_PAH.size();++i){
		
        //Create a copy of the shared pointers for PAHs in particles
		//if (source->m_PAH.size() > 1) {
  //          boost::shared_ptr<PAH> new_m_PAH = source->m_PAH[i];
  //          new_m_PAH->PAH_ID=source->m_PAH[i]->PAH_ID+100000;
  //          m_PAH.push_back(new_m_PAH);
  //      }
  //      //Always create new shared pointers for single PAHs
  //      else {
  //          boost::shared_ptr<PAH> new_m_PAH (source->m_PAH[i]->Clone());
  //          new_m_PAH->PAH_ID=source->m_PAH[i]->PAH_ID+100000;
  //          m_PAH.push_back(new_m_PAH);
  //      }
		
        //Create new shared pointers for single PAHs or PAHs in particles
        boost::shared_ptr<PAH> new_m_PAH (source->m_PAH[i]->Clone());
        //m_PAH[i]=source->m_PAH[i]->Clone();
        new_m_PAH->PAH_ID=source->m_PAH[i]->PAH_ID+100000;
        m_PAH.push_back(new_m_PAH);
//    //plus 100000 to tell which pah is cloned and also according to Id, we could easily calculate how many times this pah duplicate
        }
    }
    else m_PAH.assign(source->m_PAH.begin(),source->m_PAH.end());
}
	

/*!
 * Select a primary uniformly at random from this particle
 * and descend the aggregate tree to find the primary.
 * Note that most PAHPrimarys are nodes in a tree representing
 * connectivity within an aggregate, so it is necessary to
 * descend the tree to find a primary that really is a
 * primary.
 *
 * @param[in,out]   rng     Random number generator
 *
 * @return      Pointer to an object representing a phsyical primary
 */
PAHPrimary *PAHPrimary::SelectRandomSubparticle(rng_type &rng)
{
    // We want to choose an integer uniformly on the range [0, m_numprimary - 1]
    typedef boost::uniform_smallint<int> uniform_integer;
    uniform_integer uniformDistrib(0, m_numprimary - 1);

    // Now build an object that can generate a sample and use it
    boost::variate_generator<rng_type&, uniform_integer> uniformGenerator(rng, uniformDistrib);
    return SelectRandomSubparticleLoop(uniformGenerator());

}
/*!
 * @param[in] target The primary to be selected
 *
 * @return      Pointer to the next node down the tree on the path to target primary
*/
PAHPrimary *PAHPrimary::SelectRandomSubparticleLoop(int target)
{
	if (m_leftchild==NULL) return this;
	if (target<=m_leftchild->m_numprimary)
	{
		return m_leftchild->SelectRandomSubparticleLoop(target);
	}
	else
	{
		return m_rightchild->SelectRandomSubparticleLoop(target-(m_leftchild->m_numprimary));
	}
}


/*!
 * Each node contains two pointers (m_leftparticle and m_rightparticle)
 * to primary particles that are connected by this node
 * This function is used when the entire particle tree is duplicated.
 * It sets the pointers in the copied node (this), that the connectivity
 * of the primary particles in this node is the same as in the original node.
 *
 * A large number of asserts are provided to check that the pointer network
 * describing the aggregate structure is valid.  This can be re-enabled if
 * there is suspicion of a bug,
 *
 *
 *@todo give this method a more accurate name
 *
 * @param[in] original    Pointer to the primary to be copied
*/
void PAHPrimary::UpdateAllPointers( const PAHPrimary *original)
{
    //the primary has no children => there are no left and right particles
    if (original->m_leftchild == NULL) {
        // Check that both children are missing in the original and the copy
        //assert(original->m_rightchild == NULL);
        //assert(m_leftchild == NULL);
        //assert(m_rightchild == NULL);

        // Check the original really does not have left and right particles
        // since it should be representing a physical primary, not a connection
        // inside an aggregate unit.
        //assert(original->m_leftparticle == NULL);
        //assert(original->m_rightparticle == NULL);

        // Since this is not a connecting node it does not have left and
        // right particles.
		m_leftparticle=NULL;
		m_rightparticle=NULL;
    } else {
//        assert(original->m_rightchild != NULL);
//        assert(m_leftchild != NULL);
//        assert(m_rightchild != NULL);

        // Check the original really does have left and right particles
        // since it should be representing a connection
        // inside an aggregate unit.
//        assert(original->m_leftparticle != NULL);
//        assert(original->m_rightparticle != NULL);

        // m_leftparticle should be a leaf node so assert that it does not have children
//        assert(original->m_leftparticle->m_leftchild == NULL);
//        assert(original->m_leftparticle->m_rightchild == NULL);

        // Find the route to m_leftparticle in the original tree
        std::stack<bool> route = recordPath(original->m_leftparticle, original);

        // Now follow the same route down the new tree to find the new left particle
        m_leftparticle = descendPath(this, route);

        // the new m_leftparticle should be a leaf node so assert
        // that it does not have children
//        assert(m_leftparticle->m_leftchild == NULL);
//        assert(m_leftparticle->m_rightchild == NULL);

        // m_rightparticle should be a leaf node so assert that it does not have children
//        assert(original->m_rightparticle->m_leftchild == NULL);
//        assert(original->m_rightparticle->m_rightchild == NULL);

        // Find the route to m_rightparticle in the original tree
        route = recordPath(original->m_rightparticle, original);

        // Now follow the same route down the new tree to find the new right particle
        m_rightparticle = descendPath(this, route);

        // the new m_rightparticle should be a leaf node so assert
        // that it does not have children
//        assert(m_rightparticle->m_leftchild == NULL);
//        assert(m_rightparticle->m_rightchild == NULL);
	}
}
//void PAHPrimary::ReleasePAH(Primary &rhs)
//{
//    PAHPrimary *rhsparticle = NULL;
//    rhsparticle = dynamic_cast<AggModels::PAHPrimary*>(&rhs);
//    rhsparticle->m_PAH.clear();
//}


/*!
 * Combines this primary with another.
 *
 * \param[in]       rhs         Particle to add to current instance
 * \param[in,out]   rng         Random number generator
 *
 * \return      Reference to the current instance after rhs has been added
 */
PAHPrimary &PAHPrimary::Coagulate(const Primary &rhs, rng_type &rng)
{
    vector<PAH>::const_iterator j;
    const PAHPrimary *rhsparticle = NULL;
    rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);

    //only one PAH in rhs or this particle -> condensation or inception process.
    if ( (rhsparticle->m_numPAH==1) || (m_numPAH==1 ))
    {
        if (rhsparticle->Numprimary()>1)
        {
            // rhsparticle is a soot particle but this paricle repesents a PAH
            // this paricle will condense on this rhsparticle
            // Get a copy of the rhs ready to add to the current particle
            PAHPrimary copy_rhs(*rhsparticle);
            PAHPrimary *target = copy_rhs.SelectRandomSubparticle(rng);

            target->m_PAH.insert(target->m_PAH.end(),m_PAH.begin(),m_PAH.end());
			target->UpdatePrimary();
            CopyParts(&copy_rhs);
            CopyTree(&copy_rhs);
        }
        else
        {
            // this paricle is a soot particle but rhsparticle repesents a PAH
            // rhsparticle will condense on this particle
            // particle has more then one primary select the primary where
            // the PAH condenses to and add it to the list
            if (m_leftchild!=NULL)
            {
                PAHPrimary *target = SelectRandomSubparticle(rng);

                target->m_PAH.insert(target->m_PAH.end(),rhsparticle->m_PAH.begin(),rhsparticle->m_PAH.end());
                target->UpdatePrimary();

            }
            else
            {
                //this particle and rhsparticle are both PAHs, this process should be inception
                m_PAH.insert(m_PAH.end(),rhsparticle->m_PAH.begin(),rhsparticle->m_PAH.end());
                UpdatePrimary();
            }
        }
		UpdateCache();
        //Check the coalescence ratio
        CheckRounding();

	}

    else
    {
            //coagulation process
            PAHPrimary *newleft = new PAHPrimary;
		    PAHPrimary *newright = new PAHPrimary;
            PAHPrimary copy_rhs(*rhsparticle);
		    //rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);

           // bool print=false;
           // if (this->m_numprimary>2 && notconstpah->m_numprimary>2)
           // {
            //    PrintTree("before1");
           //     notconstpah->PrintTree("before2");
             //   print=true;
              //  cout << "printing tree"<<endl;
            //}

            //select where to add the second particle
            boost::bernoulli_distribution<> bernoulliDistrib;
            boost::variate_generator<rng_type &, boost::bernoulli_distribution<> > leftRightChooser(rng, bernoulliDistrib);
		    if (leftRightChooser())
		    {
			    newleft->CopyParts(this);
			    newright->CopyParts(&copy_rhs);
		    }
		    else
		    {
			    newright->CopyParts(this);
			    newleft->CopyParts(&copy_rhs);
		    }
            //set the pointers
		    m_leftchild=newleft;
		    m_rightchild=newright;
		    newright->m_parent=this;
		    newleft->m_parent=this;
            //set the pointers to the parent node
		    if (newleft->m_leftchild!=NULL)
		    {
			    newleft->m_leftchild->m_parent=newleft;
			    newleft->m_rightchild->m_parent=newleft;
		    }
		    if (newright->m_leftchild!=NULL)
		    {
			    newright->m_leftchild->m_parent=newright;
			    newright->m_rightchild->m_parent=newright;
		    }
            m_children_roundingLevel=0;
		    UpdateCache();
            //select the primaries that are touching
		    this->m_leftparticle=m_leftchild->SelectRandomSubparticle(rng);
		    this->m_rightparticle=m_rightchild->SelectRandomSubparticle(rng);

            //set the sintertime for the new created primary particle
            SetSinteringTime(std::max(this->m_sint_time, rhsparticle->m_sint_time));
            //initialise the variables used to calculate the coalesence ratio
            m_children_vol=m_leftparticle->m_vol+m_rightparticle->m_vol;
            m_children_surf=(m_leftparticle->m_surf+m_rightparticle->m_surf);
            m_leftparticle_vol_old=m_leftparticle->m_vol;
            m_rightparticle_vol_old=m_rightparticle->m_vol;
            m_leftparticle_numPAH=m_leftparticle->m_numPAH;
            m_rightparticle_numPAH=m_rightparticle->m_numPAH;
            m_children_radius=pow(3.0/(4.0*PI)*(m_children_vol),(ONE_THIRD));
            m_children_roundingLevel=CoalescenceLevel();
            CheckRounding();
            //must set all the pointer to NULL otherwise the delete function
            //will also delete the children
	        copy_rhs.m_leftchild=NULL;
	        copy_rhs.m_rightchild=NULL;
	        copy_rhs.m_parent=NULL;
            copy_rhs.m_leftparticle=NULL;
            copy_rhs.m_rightparticle=NULL;
		//	rhsparticle->Clear();
        //    if (fabs(surfacebeforerhs+surfacebefore-m_surf)/m_surf > 1e-6)
        //    {
        //        cout << "error" << surfacebeforerhs<<' ' <<surfacebefore << ' '<<m_surf<<endl;
        ////         PrintTree("after");
        //    }
        // if (print)
    //     PrintTree("after");

	}
    return *this;
}

//returns the CoalescenceLevel of the two primaries that are connected by this node
double PAHPrimary::CoalescenceLevel()
{
    if (m_leftparticle!=NULL)
    {

        double dV=0;


        //make sure that the volume growth does not arise from a former coalescence event
        //(number of PAHs must remain const, or increase by 1 in case of a condensation
        // event
        // if the volume changes due to a coalescence event the chilren surface is changed in
        // PAHPrimary::ChangePointer(PAHPrimary *source, PAHPrimary *target)

        // m_leftparticle->m_numPAH is the actual number of PAHs in the left one of the two
        // touching primary particles.

        // m_leftparticle_numPAH is the the number of PAHs in the left one of the two touching
        // particles last time this method was called.

        // The following test checks that at most one PAH has been added to the left one
        // of the two touching primaries since this method was last called.
        if (m_leftparticle->m_numPAH - m_leftparticle_numPAH<2)
        {
            //calculate the volume change of the left particle in the last timestep
            dV+= m_leftparticle->m_vol - m_leftparticle_vol_old;
        }

        //adjust the number of PAHs and vol for the next step
        m_leftparticle_numPAH=m_leftparticle->m_numPAH;
        m_leftparticle_vol_old=m_leftparticle->m_vol;

        //make sure that the volume growth does not arise from a coalescence event
        //(number of PAHs must remain const, or increase by1 in case of a condensation
        // event
        if (m_rightparticle->m_numPAH-m_rightparticle_numPAH<2)
        {
            //calculate the volume change of the right particle in the last timestep
            dV+= m_rightparticle->m_vol - m_rightparticle_vol_old;
        }

        //adjust the number of PAHs and vol for the next step
        m_rightparticle_numPAH=m_rightparticle->m_numPAH;
        m_rightparticle_vol_old=m_rightparticle->m_vol;

        //update the children volume, this value is used to calculate dV n the next timestep
        m_children_vol=m_rightparticle->m_vol+m_leftparticle->m_vol;
        //calculate dS, at the moment it is assumed that the particles always grow
        double ct=m_pmodel->Components(0)->CoalescThresh();
        const double dS=dV*ct/m_children_radius;
        //update the radius for the next event
        m_children_radius=pow(3.0/(4.0*PI)*(m_children_vol),(ONE_THIRD));
        m_children_surf+=dS;

       // const double spherical_surface=4*PI*m_children_radius*m_children_radius;
       //// double two_1_3=pow(2,-1*ONE_THIRD);
       // const double two_1_3=0.79370052231642452;
       // double clevel= ((spherical_surface/m_children_surf)-two_1_3)/(1-two_1_3);
       // if (clevel<0)
       //     return 0;
       // else
       //     return clevel;
       return RoundingLevel();
    }
    else
        return 0;
}

//calculates the fractal dimension of the particle and stores it in m_fdim
//void PAHPrimary::CalcFractalDimension()
//{
//	Sweep::Imaging::ParticleImage img;
//    // construct the particle by colliding the primary particles
//	img.constructSubParttree(this);
//	double L,W;
//    // calculate the length and the width of the particle
//    img.LengthWidth(L,W);
//    // calculate the radius of gyration
//    m_Rg=img.RadiusofGyration();
//	m_sqrtLW=sqrt(L*W);
//	m_LdivW=L/W;
//    m_Rg=m_Rg*1e-9;
//    m_fdim=log((double)m_numprimary)/log(2*m_Rg/(m_primarydiam/m_numprimary));
//  /*  if (m_fdim>0 && m_fdim<3)
//    {
//       string filename;
//       filename=cstr(m_fdim)+".3d";
//       ofstream out;
//       out.open(filename.c_str());
//       img.Write3dout(out,0,0,0);
//       out.close();
//    }*/
//
//}

// sets all the childrenproperties to zero, this function is used after the children are coalesced
void PAHPrimary::ResetChildrenProperties()
{
            m_children_roundingLevel=0.0;
            m_children_surf=0.0;
            m_children_vol=0.0;
            m_rightparticle_vol_old=0.0;
            m_leftparticle_vol_old=0.0;
            m_leftparticle_numPAH=0;
            m_rightparticle_numPAH=0;
            m_avg_coalesc=0;
}

//merges the left and the right primary particle to one primary particle
void PAHPrimary::Merge()
{
//    if(this->m_numprimary>5)
 //       cout<<"tset";
 //      PrintTree("before.inp");

    vector<PAH>::const_iterator j;
      //make sure this primary has children to merge
	  if(m_leftchild!=NULL)
           {
		if ( m_leftchild==m_leftparticle && m_rightchild==m_rightparticle)
		{
            //this node has only two primaries in its subtree
            //it is possible that this node is not the root node and belongs to a bigger particle
            //copy the PAHs of both children to the parent node

            //m_PAH is empty, therefore no need to append
            m_PAH=m_rightparticle->m_PAH;

            m_PAH.insert(this->m_PAH.end(),m_leftparticle->m_PAH.begin(),m_leftparticle->m_PAH.end());
            //update the pointers that pointed to the two former children
			ChangePointer(m_leftchild,this);
			ChangePointer(m_rightchild,this);

			//delete the children (destructor is recursive for this class)
            delete m_leftchild;
            delete m_rightchild;
			m_leftchild=NULL;
			m_rightchild=NULL;
			m_leftparticle=NULL;
			m_rightparticle=NULL;

            //set the children properties to zero, this node has no more children
            ResetChildrenProperties();
            UpdatePrimary();
		}


		else
		{
			if (m_leftchild->m_numprimary<m_rightchild->m_numprimary)
			{
				//append to left subtree because there are fewer primaries
                //this is only to keep the tree balanced
				PAHPrimary *oldleftparticle=m_leftparticle;

                //copy the PAHs

				//for (j=oldleftparticle->m_PAH.begin(); j!=oldleftparticle->m_PAH.end(); ++j) {
				//	m_rightparticle->m_PAH.insert(m_rightparticle->m_PAH.end(),PAH(*j));
				//}
                m_rightparticle->m_PAH.insert(m_rightparticle->m_PAH.end(),oldleftparticle->m_PAH.begin(),oldleftparticle->m_PAH.end());
                m_rightparticle->UpdatePrimary();
                //set the pointers from the leftprimary to the rightprimary
                //this will be the new bigger primary
				oldleftparticle->ChangePointer(oldleftparticle,m_rightparticle);
                m_rightparticle->ChangePointer(m_rightparticle,m_rightparticle);
                //set the pointer to the parent node
				if (oldleftparticle->m_parent->m_leftchild==oldleftparticle)
				{
					oldleftparticle->m_parent->m_leftchild=m_rightchild;
				}
				else
				{
					oldleftparticle->m_parent->m_rightchild=m_rightchild;
				}
				m_rightchild->m_parent=oldleftparticle->m_parent;


				PAHPrimary *oldleftchild=m_leftchild;
				PAHPrimary *oldparent=m_parent;

				//copy the properties of the former leftchild to this node
				// so that it can be removed from the aggregate tree structure
				CopyParts(oldleftchild);

				// Now break the links to the tree structure in oldleftchild and free it
				oldleftchild->m_leftchild = NULL;
				oldleftchild->m_rightchild = NULL;
				delete oldleftchild;

				m_parent=oldparent;
                if (m_leftchild!=NULL)
                {
				    m_rightchild->m_parent=this;
				    m_leftchild->m_parent=this;
                }

                delete oldleftparticle;

			}

			else
			{
				//append to right subtree
				PAHPrimary *oldrightparticle=m_rightparticle;

			//	for (j=oldrightparticle->m_PAH.begin(); j!=oldrightparticle->m_PAH.end(); ++j) {
			//		m_leftparticle->m_PAH.insert(m_leftparticle->m_PAH.end(),PAH(*j));
			//	}
                m_leftparticle->m_PAH.insert(m_leftparticle->m_PAH.end(),oldrightparticle->m_PAH.begin(),oldrightparticle->m_PAH.end());
                m_leftparticle->UpdatePrimary();
                oldrightparticle->ChangePointer(oldrightparticle,m_leftparticle);
                m_leftparticle->ChangePointer(m_leftparticle,m_leftparticle);
				if (oldrightparticle->m_parent->m_leftchild==oldrightparticle)
				{
					oldrightparticle->m_parent->m_leftchild=m_leftchild;
				}
				else
				{
					oldrightparticle->m_parent->m_rightchild=m_leftchild;
				}
				m_leftchild->m_parent=oldrightparticle->m_parent;

				//ReleaseMem(oldrightparticle);
				PAHPrimary *oldrightchild=m_rightchild;
				PAHPrimary *oldparent=m_parent;

                //copy the properties of the former leftchild to this node
                // so that it can be removed from the aggregate tree structure
				CopyParts(oldrightchild);

                // Now break the links to the tree structure in oldrightchild and free it
                oldrightchild->m_leftchild = NULL;
                oldrightchild->m_rightchild = NULL;
                delete oldrightchild;

				m_parent=oldparent;
                if (m_leftchild!=NULL)
                {
				    m_rightchild->m_parent=this;
				    m_leftchild->m_parent=this;
                }
                delete oldrightparticle;
			}
		}
        UpdateCache();
  //      PrintTree("after.inp");
        }
}

void PAHPrimary::ReleaseMem()
{ 
    m_PAH.clear();
}
/*!//add explanation
 * @param[in] source Pointer to the original particle
 * @param[in,out] target Pointer to the new particle
*/
void PAHPrimary::ChangePointer(PAHPrimary *source, PAHPrimary *target)
{
		if(m_rightparticle==source){
			m_rightparticle=target;
            double sphericalsurface=
                4*PI*pow(3*(m_leftparticle->Volume()+m_rightparticle->Volume())/(4*PI),TWO_THIRDS);
            m_children_surf=sphericalsurface/(m_children_roundingLevel*0.2063+0.7937);    //sphericalsurface/(m_children_roundingLevel*(1-2^(-1/3))+2^(-1/3))
		}
		if(m_leftparticle==source){
			m_leftparticle=target;
            double sphericalsurface=
                4*PI*pow(3*(m_leftparticle->Volume()+m_rightparticle->Volume())/(4*PI),TWO_THIRDS);
            m_children_surf=sphericalsurface/(m_children_roundingLevel*0.2063+0.7937);    //sphericalsurface/(m_children_roundingLevel*(1-2^(-1/3))+2^(-1/3))

		}

    // Update the tree above this sub-particle.
    if (m_parent != NULL) {
        m_parent->ChangePointer(source,target);
    }

}   

/*!
 *@param[in]        t       Time upto which to update
 *@param[in]        model   Particle model defining interpretation of particle data
 *@param[in]        sys     Cell containing particle and providing gas phase
 *@param[in,out]    rng     Random number generator
 *
 * The actual interval over which the update is carried out on a PAH is from
 * lastupdated to t.
 */
void PAHPrimary::UpdatePAHs(const double t, const Sweep::ParticleModel &model,Cell &sys, rng_type &rng)
{
    // Either the primary has two children or it is a leaf of the
    // tree
    if (m_leftchild!=NULL)
    {
        // Recurse down to the leaves
        m_leftchild->UpdatePAHs(t, model,sys, rng);
        m_rightchild->UpdatePAHs(t, model,sys, rng);
    }
    else
    {
        // There are PAHs in this primary so update them, if needed
        // Flag to show if any PAH has been changed
        // this flag store the condition of this cluster(primary particle)
        bool m_PAHclusterchanged = false;
        // this flag store the condition of a particular PAH within this cluster
        bool m_InvalidPAH = false;
        // check whether this updated primary particle is a inceptedPAH, used later to track the num of InceptedPAH in the ensemble
        const int m_InceptedPAH = InceptedPAH();

        // Loop over each PAH in this primary
        const std::vector<boost::shared_ptr<PAH> >::iterator itEnd = m_PAH.end();
        for (std::vector<boost::shared_ptr<PAH> >::iterator it = m_PAH.begin(); it != itEnd; ++it) {

            // This is a model parameter that defines when primary particles
            // contain too many PAHs for them all to have full access to the
            // surrounding gas phase.  Once a primary contains more than minPAH
            // PAHs the growth rate of each PAH is reduced according to the
            // growth factor
            bool m_PAHchanged = false;
            const double minPAH = model.Components(0)->MinPAH();

            double growthfact = 1.0;
            if (m_numPAH>=minPAH)
            {
                // concentration of gasphase species is reduced by this factor to 
                // reflect the slower growth of molecules that are closely surrounded
                // by many other PAHs.
                growthfact = model.Components(0)->GrowthFact();
            }

            // Time for one particular PAH to grow
            const double growtime = t - (*it)->lastupdated;
            assert(growtime >= 0.0);

            const int oldNumCarbon = (*it)->m_pahstruct->numofC(); 
            const int oldNumH = (*it)->m_pahstruct->numofH();

            // Here updatePAH function in KMC_ARS model is called.
            // waitingSteps is set to be 1 by dc516, details seeing KMCSimulator::updatePAH()
            sys.Particles().Simulator()->updatePAH((*it)->m_pahstruct, (*it)->lastupdated, growtime, 1,
                                                   rng, growthfact, (*it)->PAH_ID);

            (*it)->lastupdated=t;

            // See if anything changed, as this will required a call to UpdatePrimary() below
            if(oldNumCarbon != (*it)->m_pahstruct->numofC() || oldNumH != (*it)->m_pahstruct->numofH())
            {
                m_PAHclusterchanged = true; 
                m_PAHchanged = true;
            }
			// the second condition ensures that the m_InvalidPAH is modified in the correct way
			// consider 2 PAH, the first is invalid, the other is valid. If not using the 
			// second condition, the m_InvalidPAH will be false enventually, but it should be true.
            if (m_PAHchanged && !m_InvalidPAH)
                m_InvalidPAH = CheckInvalidPAHs(*it);
        }
        if (m_InvalidPAH)
        {
            RemoveInvalidPAHs();
        }

        //sys.Particles().Add();
        //this->ParticleModel()->Mode();
        // Calculate derived quantities such as collision diameter and surface
        // area by iterating through all the PAHs.  This call is rather expensive.
        if(m_PAHclusterchanged) {
            UpdatePrimary();
            // if (m_InceptedPAH!=0) {
            //        sys.Particles().SetNumOfInceptedPAH(-1);
            //        //sys.Particles().NumOfInceptedPAH();
            // }
            //else if (m_InceptedPAH == 0 && InceptedPAH()==1){
            //        sys.Particles().SetNumOfInceptedPAH(1);
            //        //sys.Particles().NumOfInceptedPAH();
            //}
        }
        // otherwise there is no need to update
    }
}

// currently, only A1,A2 and A4 can be specified as InceptedPAH due to the underlying KMC code
bool PAHPrimary::CheckInvalidPAHs(const boost::shared_ptr<PAH> & it) const
{
    ParticleModel::PostProcessStartingStr str = this->ParticleModel()->InceptedPAH();
    int m_control;
    switch (str){
    case ParticleModel::A1:
        m_control=Sweep::KMC_ARS::BENZENE_C;
        break;
    case ParticleModel::A2:
        m_control=10;
        break;
    case ParticleModel::A4:
        m_control=Sweep::KMC_ARS::PYRENE_C;
        break;
    case ParticleModel::A5:
        m_control=Sweep::KMC_ARS::BENZOPYRENE_C;
        break;
    default:
        throw std::runtime_error("no information about the incepted PAH is available (Sweep::PAHPrimary::CheckInvalidPAHs())");
    }
    // if the PAH in the cluster is as the same size as the incepted PAH, it will be released but the current implementation is directly removed which causes mass loss, for a fully coupled model this part should be redesigned.
    return (it->m_pahstruct->numofC() < m_control||(it->m_pahstruct->numofC() <= m_control && NumPAH()!=1));
}

//struct compare_class
//{
//    bool operator () (const boost::shared_ptr<PAH>  & it) const
//    {
//        if (m_pyreneInception) 
//        {
//            if (it->Structure()->numofC() < Sweep::KMC_ARS::PYRENE_C)
//                return true;
//        }
//        else 
//        {
//            if (it->Structure()->numofC() < Sweep::KMC_ARS::BENZENE_C)
//                return true;
//        }
//        return false;
//    }
//};

void PAHPrimary::RemoveInvalidPAHs()
{
    // if the num of PAH in this primary particle is smaller than 2 after removing the InvalidPAHs,
    // FakeRounding() will return true and it will be merged with other part of the aggregate.
    // TODO: thinking suitable solution for moving the InvalidPAH back to gasphase when coupling with the gasphase chemistry.
    // current implementation is only suitable for post-pocessing due to mass loss.
    //remove_if(m_PAH.begin(),m_PAH.end(),compare_class());
    std::vector<boost::shared_ptr<PAH> >::iterator NewEnd = remove_if(m_PAH.begin(),m_PAH.end(),boost::bind(&PAHPrimary::CheckInvalidPAHs, this,_1));
    m_PAH.resize(NewEnd-m_PAH.begin());
}

// Find Xmer (monomer, dimer or trimer) and record its mass which will be used to create mass spectra
// however, this function currently is only implenmented in PAHPrimary class
// this means it is limited to PAH-PP model. 
void PAHPrimary::FindXmer(std::vector<double> &out, int m_xmer) const
{
    if (m_leftchild!=NULL)
        m_leftchild->FindXmer(out, m_xmer);

	// dimer: a partilce with one primary containing 2 PAHs
	if (this->NumPAH()==m_xmer)
		out.push_back(this->MassforXmer());
}

// Find particular primary particle with target number of PAHs ( this can be a range, details please see the implementation) 
// and record num of C and H for each PAH, also push_back a divider (0,0,ID) at end to distinguish 
void PAHPrimary::FindXmer(std::vector<std::vector<double> > &out, int target_num_PAH) const
{
	if (m_leftchild != NULL)
		m_leftchild->FindXmer(out,target_num_PAH);
	if (m_rightchild != NULL)
		m_rightchild->FindXmer(out,target_num_PAH);
    if (target_num_PAH-10<=0) std::cout<<"Warning: (target_num_PAH-10) return 0 or negative value in PAHPrimary::FindXmer()"<<std::endl;
	if (m_PAH.size() <= size_t(target_num_PAH+10) && m_PAH.size() >= size_t(target_num_PAH-10))
		mass_PAH(out);
}

//calculate the mass of Xmer including C and H
double PAHPrimary::MassforXmer() const
{ 
	int sum = 0;
	for (size_t i = 0; i != m_PAH.size(); ++i)
	{
		sum+=m_PAH[i]->m_pahstruct->numofH();
	    sum+=12 * m_PAH[i]->m_pahstruct->numofC();
	}
	return sum;	
}

int PAHPrimary::InceptedPAH() const
{
    // m_parent == NULL is to check whether this primary particle is part of an aggregate.
    if (m_parent == NULL && Numprimary() == 1 && NumPAH() == 1){
        //currently only Num of C and H is used to identify the Pyrene, Naphthalene and benzene
        ParticleModel::PostProcessStartingStr str = ParticleModel()->InceptedPAH();
        switch (str){
        case ParticleModel::A1:
            if (NumCarbon() == BENZENE_C && NumHydrogen() == BENZENE_H)
                return 1;
            else return 0;
            break;
        case ParticleModel::A2:
            if (NumCarbon() == NAPHTHALENE_C && NumHydrogen() == NAPHTHALENE_H)
                return 1;
            else return 0;
            break;
        case ParticleModel::A4:
            if (NumCarbon() == PYRENE_C && NumHydrogen() == PYRENE_H)
                return 1;
            else return 0;
            break;
        case ParticleModel::A5:
            if (NumCarbon() == BENZOPYRENE_C && NumHydrogen() == BENZOPYRENE_H)
                return 1;
            else return 0;
            break;
        default:
            return 0;
        }
    }
    else return 0;
}

// dump information of this Xmer to a vector<vector<double> > 
 void PAHPrimary::mass_PAH(std::vector<std::vector<double> > &out) const
 {
    std::vector<double> temp;
    std::vector<double> divider(2,0);
    for (size_t i = 0; i != m_PAH.size(); ++i)
    {
        temp.push_back(m_PAH[i]->m_pahstruct->numofC());
        temp.push_back(m_PAH[i]->m_pahstruct->numofH());
        m_PAH[i]->saveDOTperLoop((int)ID,(int)i);
        out.push_back(temp);
        temp.clear();
    }
    divider.push_back(ID);
    out.push_back(divider);
    ID++;
}
 // only for num of Primary == 1, carbon only
double PAHPrimary::ReducedMass() const
{
    double val=0;
    for (size_t i = 0; i != m_PAH.size(); ++i)
    {
        int num_C=0;
        num_C=m_PAH[i]->m_pahstruct->numofC();
        val+=1/num_C;
    }
    if (val==0) return 1;
    else 
    return 1/val;
}
bool IsSticked(double val)
{
    return val<2*(32*12+14);
}
bool NotValid(double val)
{
    return val==0;
}

void PAHPrimary::Fragtest(std::vector<double> &out, const int k, std::string mode, double threshold) const
{
    //test8 start
    //int thres=0;
    //if (mode =="MIN"||mode=="MAX")
    //    thres = threshold/2;
    //else if (mode =="COMBINED")
    //    thres = threshold/4;
    //else if (mode =="REDUCED")
    //    thres = threshold;
    //if (m_leftchild!=NULL)
    //m_leftchild->Fragtest(out, k, mode, threshold);

    //std::vector<double> temp;

    ////if (k==0 && ReducedMass()<=thres) {
    ////    for (size_t i = 0; i != m_PAH.size(); i+=1)
    ////        {
    ////            int num_C=0;
    ////            int num_H=0;
    ////            int val=0;
    ////            int val1=0;
    ////            num_C=m_PAH[i]->m_pahstruct->numofC();
    ////            num_H=m_PAH[i]->m_pahstruct->numofH();
    ////            // PAH mass (u)
    ////            val = 12*num_C + num_H;
    ////            temp.push_back(val);
    ////        }
    ////}
    ////else if (k==1 && ReducedMass()>thres) {
    
    //if (k==1 && NumPAH()<=10) {
    //    for (size_t i = 0; i != m_PAH.size(); i+=2)
    //    {
    //        if (i+1>=m_PAH.size())
    //            break;
    //        int num_C=0;
    //        int num_H=0;
    //        int val=0;
    //        int val1=0;
    //        num_C=m_PAH[i]->m_pahstruct->numofC();
    //        num_H=m_PAH[i]->m_pahstruct->numofH();
    //        // PAH mass (u)
    //        val = 12*num_C + num_H;
    //        num_C=m_PAH[i+1]->m_pahstruct->numofC();
    //        num_H=m_PAH[i+1]->m_pahstruct->numofH();
    //        // PAH mass (u)
    //        val1 = 12*num_C + num_H;
    //        val+=val1;
    //        temp.push_back(val);
    //    }
    //}
    //if (temp.size()!=0) {
    //    std::vector<double>::iterator NewEnd = remove_if(temp.begin(),temp.end(),IsSticked);
    //    temp.resize(NewEnd-temp.begin());
    //    out.insert(out.end(),temp.begin(),temp.end());
    //}
    //test8 end
    //test9 start
    std::vector<double> temp;
    if (k==1 && NumPAH()<=5) {
        mass_PAH(temp);
    }
    if (temp.size()!=0)
    {
        for (size_t i = 0; i != static_cast<size_t>(NumPAH()); ++i){
            //c32h14
            if (temp[i]<=398)
                temp[i]=0;
        }
        std::vector<double>::iterator NewEnd = remove_if(temp.begin(),temp.end(),NotValid);
        temp.resize(NewEnd-temp.begin());
        for (size_t i = 0; i !=temp.size(); i+=2) {
            if (i+1>=temp.size()) {
                temp[i]=0;
                break;
            }
            temp[i]=temp[i]+temp[i+1];
            temp[i+1]=0;
        }
        NewEnd = remove_if(temp.begin(),temp.end(),NotValid);
        temp.resize(NewEnd-temp.begin());
        out.insert(out.end(),temp.begin(),temp.end());
    }
    //test9 end

}

// use to output detailed info about particular PAH
void PAHPrimary::OutputPAHPSL(std::vector<std::vector<double> > &out, const int index, const double density) const
{
    if (m_leftchild!=NULL)
        m_leftchild->OutputPAHPSL(out, index, density);
    if (m_rightchild!=NULL) m_rightchild->OutputPAHPSL(out, index, density);


    std::vector<double> temp;

    for (size_t i = 0; i != m_PAH.size(); ++i)
    {
        int num_C=0;
        int num_H=0;
        double val=0.0;
        double m_mass=0.0;
        double PAHCollDiameter=0.0;
        double diameter=0.0;
        temp.push_back(index);

        num_C=m_PAH[i]->m_pahstruct->numofC();
        temp.push_back(num_C);
        num_H=m_PAH[i]->m_pahstruct->numofH();
        temp.push_back(num_H);
        temp.push_back(m_PAH[i]->m_pahstruct->numofRings());
        temp.push_back(m_PAH[i]->m_pahstruct->numofRings5());
        temp.push_back(m_PAH[i]->m_pahstruct->numofEdgeC());

        // PAH mass (u)
        val = 12*num_C + num_H;
        temp.push_back(val);
        // PAH mass (kg)
        m_mass = num_C*1.9945e-26 + num_H*1.6621e-27;
        temp.push_back(m_mass);

        // PAHCollDiameter (Angstrom)
        PAHCollDiameter = sqrt(num_C*2.0/3.);
        PAHCollDiameter *= 2.4162*1e-10;         //convert from Angstrom to m
        temp.push_back(PAHCollDiameter);

        // PAH denbsity (kg/m3)
        temp.push_back(density);

        // PAH volume (m3)
        val = m_mass / density;
        temp.push_back(val);

        // diameter (m)
        diameter = pow(6.0 * val / PI, ONE_THIRD);
        temp.push_back(diameter);

        // collision diameter (m)
        val = max(diameter, PAHCollDiameter);
        temp.push_back(val);

        // time created (s)
        val = m_PAH[i]->time_created;
        temp.push_back(val);

        // index of PAH
        val = m_PAH[i]->PAH_ID;
        temp.push_back(val);

        out.push_back(temp);
        temp.clear();
    }
}

// this function is only used to create a vector containing all the mass of individual PAH within this soot particle
void PAHPrimary::mass_PAH(std::vector<double> &out) const
{
    if (m_leftchild != NULL)
        m_leftchild->mass_PAH(out);
    if (m_rightchild != NULL)
        m_rightchild->mass_PAH(out);

    double temp_mass=0.0;
    for (size_t i = 0; i != m_PAH.size(); ++i)
    {
        temp_mass = 12*m_PAH[i]->m_pahstruct->numofC() + m_PAH[i]->m_pahstruct->numofH();
        out.push_back(temp_mass);
    }
}
void PAHPrimary::UpdateCache(void)
{
    UpdateCache(this);
}

bool PAHPrimary::FakeRounding()
{
    // there are two conditions that it should perform FakeRounding, no PAH or only one PAH in this primary particle, thus it should be merged with other primary particle.
    if (m_leftparticle!=NULL)
    {
        if (m_leftparticle->m_numPAH==1||(m_leftparticle->m_numPAH==0&&m_leftparticle->m_numcarbon==0))
            return true;
        else if (m_rightparticle->m_numPAH==1||(m_rightparticle->m_numPAH==0&&m_rightparticle->m_numcarbon==0))
            return true;
        else return false;
    }
    else return false;
}

bool PAHPrimary::CheckRounding()
{
    bool hascoalesced=false;
    if ((m_children_roundingLevel> 0.95 && m_leftparticle!=NULL)||FakeRounding())
        {
           // PrintTree("before.inp");
           // cout <<"merging"<<m_children_roundingLevel<<endl;
             Merge();
             //  PrintTree("after.inp");
             hascoalesced=true;
             //check again because this node has changed
             CheckRounding();
        }
    if (m_leftchild!=NULL)
    {
        hascoalesced=m_leftchild->CheckRounding();
        hascoalesced=m_rightchild->CheckRounding();
    }
    UpdateCache();
    return hascoalesced;
}


//this function updates a primary particle (m_PAH.empty() is not true)
void PAHPrimary::UpdatePrimary(void)
{
    if (m_PAH.empty())
    {
        // in this case this primary particle is invalid, so set the cached values, all the m_??? could be set to 0
        // but currently only the mass is concerned in the Primary::IsValid()
        m_mass = 0;
        m_numcarbon=0;
        m_numH=0;
        m_numOfEdgeC=0;
        m_numOfRings=0;
        m_PAHmass=0;
        m_PAHCollDiameter=0;
        m_numPAH= m_PAH.size();
    }
    else 
    {
        m_numcarbon=0;
        m_numH=0;
        m_numOfEdgeC=0;
        m_numOfRings=0;
        m_PAHmass=0;
        m_PAHCollDiameter=0;
        m_numPAH= m_PAH.size();

        int maxcarbon=0;
        for (vector<boost::shared_ptr<PAH> >::iterator i=m_PAH.begin(); i!=m_PAH.end(); ++i) {
            m_numcarbon += (*i)->m_pahstruct->numofC();
            m_numH += (*i)->m_pahstruct->numofH();
            m_numOfEdgeC += (*i)->m_pahstruct->numofEdgeC();
            m_numOfRings += (*i)->m_pahstruct->numofRings();
            maxcarbon=max(maxcarbon, (*i)->m_pahstruct->numofC());    // search for the largest PAH in the PRimary, in Angstrom
        }
        m_PAHmass=m_numcarbon * 1.9945e-26 + m_numH * 1.6621e-27;     //convert to kg, hydrogen atoms are not considered
        m_PAHCollDiameter=sqrt(maxcarbon*2.0/3.);
        m_PAHCollDiameter*=2.4162*1e-10;         //convert from Angstrom to m

        if(m_pmodel->ComponentCount()!=1)        //at the moment we have only one component: soot
        {
            throw std::runtime_error("Model contains more then one component. Only soot is supported. (PAHPrimary::UpdatePrimary)");
        }
        m_vol  = m_PAHmass / m_pmodel->Components(0)->Density();        //in m^3
        m_mass = m_PAHmass;
        m_diam = pow(6.0 * m_vol / PI, ONE_THIRD);
        m_dmob = m_diam;
        m_dcol = max(m_diam,m_PAHCollDiameter);
        m_surf = PI * m_diam * m_diam;
        m_primarydiam = m_diam;
        m_avg_coalesc = 0;
    }
}



void PAHPrimary::Reset()
{
    m_numcarbon=0;
    m_numH=0;
    m_numOfEdgeC=0;
    m_numOfRings=0;
    m_primarydiam=0.0;
    m_surf=0;
    m_vol=0;
    m_PAH.clear();
    m_avg_coalesc=0;
}
/*!
 * @param[in] root The root node of this particle
*/
void PAHPrimary::UpdateCache(PAHPrimary *root)
{
    //Update the children
	if (m_leftchild!=NULL)
	{
		m_leftchild->UpdateCache(root);
		m_rightchild->UpdateCache(root);
		m_numprimary=m_leftchild->m_numprimary+m_rightchild->m_numprimary;
	}
    //this is a primary and the number of primaries below this node is one (this node and no children)
	else
	{
            m_numprimary=1;

            // Check that the primary is begin kept up to date
//            const double oldNumCarbons = m_numcarbon;
//            UpdatePrimary();
//            if(m_numcarbon != oldNumCarbons)
//                std::cerr << "UpdatePrimary has changed num carbons inside UpdateCache\n";

            m_avg_coalesc=0;
	}

    //this is not a primary, sum up the properties
	if (m_leftchild!=NULL)
    {
        // remove the PAHs from this node to free memory
        Reset();
		m_surf = m_leftchild->m_surf+m_rightchild->m_surf;
		m_primarydiam = (m_leftchild->m_primarydiam+m_rightchild->m_primarydiam);
        m_vol=m_leftchild->m_vol+m_rightchild->m_vol;
        m_numPAH = m_leftchild->m_numPAH+m_rightchild->m_numPAH;
        m_primarydiam = (m_leftchild->m_primarydiam+m_rightchild->m_primarydiam);
        m_mass=(m_leftchild->m_mass+m_rightchild->m_mass);
        m_PAHCollDiameter=max(m_leftchild->m_PAHCollDiameter,m_rightchild->m_PAHCollDiameter);

        m_numcarbon=m_leftchild->m_numcarbon + m_rightchild->m_numcarbon;
        m_numH=m_leftchild->m_numH + m_rightchild->m_numH;

        m_numOfEdgeC=m_leftchild->m_numOfEdgeC + m_rightchild->m_numOfEdgeC;
        m_numOfRings=m_leftchild->m_numOfRings + m_rightchild->m_numOfRings;
        // calculate the coalescence level of the two primaries connected by this node
        m_children_roundingLevel=CoalescenceLevel();
        //sum up the avg coal level
        m_avg_coalesc=m_children_roundingLevel+m_leftchild->m_avg_coalesc+m_rightchild->m_avg_coalesc;

        // calculate the different diameters only for the root node because this goes into the
        // particle tree and gets used by the coagulation kernel
        if (this==root)
        {
             //spherical eqiv radius
            double spherical_radius=pow(3*m_vol/(4*PI),ONE_THIRD);
            m_diam=2*spherical_radius;
            // there are m_numprimary-1 connections between the primary particles
            m_avg_coalesc=m_avg_coalesc/(m_numprimary-1);
            //approxmiate the surface of the particle
            const double numprim_1_3=pow(m_numprimary,-0.333333);
            m_surf=4*PI*spherical_radius*spherical_radius/
                (m_avg_coalesc*(1-numprim_1_3)+numprim_1_3);

            //calculate the surface equivalent radius
           // const double radius_surf=sqrt(m_surf/(4*PI));
            // the average between the surface and voluem equiv diameter
            //const double meandiam=spherical_radius+radius_surf;            //2*0.5*(radius(vol)+radius(sphere))
            const double aggcolldiam=(6*m_vol/m_surf)*
                                     pow(pow(m_surf,3)/(36*PI*m_vol*m_vol),(1.0/1.8));
            // the maximum of the largest PAH diameter and
            // the average between the surface and voluem equiv diameter
            const double cdiam=max(aggcolldiam,m_PAHCollDiameter);
            m_dmob = aggcolldiam;
            SetCollDiameter(cdiam);
        }
        else
        {
            m_diam=0;
            m_dmob=0;

        }
    }
}

/*!
 * @param[in] filename Output filename
*/
void PAHPrimary::PrintTree(string filename)
{
  ofstream out;
  out.open(filename.c_str());
  out << "digraph unix {"<<endl;
  out <<"graph [rankdir = \"LR\"];"<<endl;
  PrintTreeLoop(out);
  out << "}"<<endl;
  out.close();
}

/*!
 * @param[in] out Output stream
*/
void PAHPrimary::PrintTreeLoop(std::ostream &out)
{
  if (m_leftchild!=NULL)
  { //out<<"leftchild "<<10E8*m_leftchild->SphDiameter()<<endl;
	//out<<"rightchild "<<10E8*m_rightchild->SphDiameter()<<endl;
      out<<"\" "<<this<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_roundingLevel<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this->m_leftchild<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_roundingLevel<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this->m_rightchild<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_roundingLevel<<"|"<<this<<"\"];"<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftchild<<"\"; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightchild<<"\"; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftparticle<<"\"[label=\""<<this<<"\",color=\"blue\"]; "<<endl;
	out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightparticle<<"\"[label=\""<<this<<"\",color=\"blue\"]; "<<endl;
	m_leftchild->PrintTreeLoop(out);
    m_rightchild->PrintTreeLoop(out);
  }

  else
  {
      //out<<"\" "<<this<<"\" "<<" [label = \""<<2*pow((this->vol_sinter)*3/(4*PI),ONE_THIRD)<<"\"];"<<endl;
      out<<"\" "<<this<<"\" "<<" [shape = \"record\" label = \"surf="<<this->m_surf<<"|m_children_surf="<<this->m_children_surf<<"|m_vol="<<this->m_vol<<"|"<<this->m_children_roundingLevel<<"|"<<this<<"\"];"<<endl;
  }
}


const PAHPrimary *PAHPrimary::RightChild() const
{
    return m_rightchild;
}

 const PAHPrimary *PAHPrimary::LeftChild() const
{
    return m_leftchild;
}

double PAHPrimary::PAHCollDiameter() const
{
    return m_PAHCollDiameter;
}

double PAHPrimary::Rg() const
{
    return m_Rg;
}

double PAHPrimary::Fdim() const
{
    return m_fdim;
}

double PAHPrimary::PrimaryDiam() const
{
    return m_primarydiam;
}

double PAHPrimary::LdivW() const
{
    return m_LdivW;
}

int PAHPrimary::Numprimary() const
{
    return m_numprimary;
}

int PAHPrimary::NumCarbon() const
{
    return m_numcarbon;
}

int PAHPrimary::NumHydrogen() const
{
    return m_numH;
}

int PAHPrimary::NumEdgeC() const
{
    return m_numOfEdgeC;
}

int PAHPrimary::NumRings() const
{
    return m_numOfRings;
}

int PAHPrimary::NumPAH() const
{
    return m_numPAH;
}

double PAHPrimary::sqrtLW() const
{
    return m_sqrtLW;
}

double PAHPrimary::AvgCoalesc() const
{
    return m_avg_coalesc;
}

// READ/WRITE/COPY.

// Returns a copy of the model data.
PAHPrimary *const PAHPrimary::Clone(void) const
{
    m_clone=true;
	PAHPrimary* newPAHPrimary=new PAHPrimary();
	newPAHPrimary->CopyParts(this);
	if (this->m_leftchild!=NULL)
	newPAHPrimary->CopyTree(this);
	m_clone=false;
	return newPAHPrimary;
}


// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType PAHPrimary::AggID(void) const {return AggModels::PAH_KMC_ID;}

/*
 * @brief Writes a object to a binary stream
 *
 * @param[in,out]    out                 Output binary stream
 * @param[in,out]    duplicates          Addresses of PAHs that have already been serialised
 *
 * @exception		 invalid_argument    Stream not ready
 */
void PAHPrimary::Serialize(std::ostream &out, void *duplicates) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        if (m_pmodel->WriteBinaryTrees()) {
            // Call the binary tree serialiser...
            BinTreeSerializer <PAHPrimary> tree;
            tree.Serialize(out, this, duplicates);
        } else {
            // Just serialise the root node.
            SerializePrimary(out, duplicates);
        }

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PAHPrimary::Serialize).");
    }
}

/*!
 * @brief Writes an individual primary to a binary stream
 *
 * @param[in,out]    out                 Output binary stream
 * @param[in,out]    duplicates          Addresses of PAHs that have already been serialised
 *
 * @exception		 invalid_argument    Stream not ready
 */
void PAHPrimary::SerializePrimary(std::ostream &out, void *duplicates) const
{
    if (out.good()) {

        int  val_int(0);
        double val(0.0);
        // Serialise state space
        val_int = m_numcarbon;
        out.write((char*)&val_int, sizeof(val_int));

        val_int = m_numH;
        out.write((char*)&val_int, sizeof(val_int));

        val_int = m_numOfEdgeC;
        out.write((char*)&val_int, sizeof(val_int));

        val_int = m_numOfRings;
        out.write((char*)&val_int, sizeof(val_int));

        val_int =  m_numPAH;
        out.write((char*)&val_int, sizeof(val_int));

        val_int = m_numprimary;
        out.write((char*)&val_int, sizeof(val_int));

        val = m_PAHmass;
        out.write((char*)&val, sizeof(val));

        val = m_PAHCollDiameter;
        out.write((char*)&val, sizeof(val));

        val = m_primarydiam;
        out.write((char*)&val, sizeof(val));

        val = m_children_radius;
        out.write((char*)&val, sizeof(val));

        val = m_children_vol;
        out.write((char*)&val, sizeof(val));

        val = m_leftparticle_vol_old;
        out.write((char*)&val, sizeof(val));

        val = m_rightparticle_vol_old;
        out.write((char*)&val, sizeof(val));

        val_int = m_rightparticle_numPAH;
        out.write((char*)&val_int, sizeof(val_int));

        val_int =  m_leftparticle_numPAH;
        out.write((char*)&val_int, sizeof(val_int));

        val = m_children_surf;
        out.write((char*)&val, sizeof(val));

        val = m_children_roundingLevel;
        out.write((char*)&val, sizeof(val));

        /*Imaging properties
        val = (double)m_Rg;
        out.write((char*)&val, sizeof(val));

        val = (double)m_fdim;
        out.write((char*)&val, sizeof(val));

        val = (double)m_sqrtLW;
        out.write((char*)&val, sizeof(val));

        val = (double)m_LdivW;
        out.write((char*)&val, sizeof(val));*/

        val = m_avg_coalesc;
        out.write((char*)&val, sizeof(val));

        val = m_sint_time;
        out.write((char*)&val, sizeof(val));

        val_int = (int) m_PAH.size();
        out.write((char*)&val_int, sizeof(val_int));

        // write the PAH stack (m_PAH)
		PahSerialisationMap *pahDuplicates = reinterpret_cast<PahSerialisationMap*>(duplicates);
        outputPAHs(out, *pahDuplicates);

        // Output base class.
        Primary::Serialize(out);

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PAHPrimary::SerializePrimary).");
    }
}

/*
 * @brief Writes individual PAHs to a binary stream 
 *
 * @param[in]        out               Output binary stream
 * @param[in,out]    pah_duplicates    Addresses of PAHs that have already been serialised
 */
void PAHPrimary::outputPAHs(std::ostream &out, PahSerialisationMap &pah_duplicates) const
{
    if (m_PAH.size() != 0) {
        //count the number of PAH should be serialized
        unsigned count = 0;
		int PAHSize= m_PAH.size();
        
        // Keep a record of which PAHs have been serialised to avoid serialising the same
        // memory location twice.
        //PahSerialisationMap pahUniqueAdresses;

        while (count != m_PAH.size())
        {
            // Serialise the raw memory locations of the PAHs, so that we
            // can recognise PAHs that are referenced multiple times on
            // deserialisation.
            void *ptr = m_PAH[count].get();
            out.write(reinterpret_cast<char*>(&ptr), sizeof(ptr));
            //if(pahUniqueAdresses.find(m_PAH[count].get()) == pahUniqueAdresses.end()) {
            //    m_PAH[count]->Serialize(out);
            //    pahUniqueAdresses.insert(m_PAH[count].get());
            //}
            //else {
            //    std::cout << "PAH at " << m_PAH[count].get() << " has already been serialised\n";
            //}
			if(pah_duplicates.find(m_PAH[count].get()) == pah_duplicates.end()){
				m_PAH[count]->Serialize(out);
				pah_duplicates.insert(m_PAH[count].get());
			}
			else {
                //std::cout << "PAH at " << m_PAH[count].get() << " has already been serialised\n";
            }
            ++count;
        }
    }
}

/* 
 * @brief Reads the object from a binary stream.
 *
 * @param[in,out]	 in		             Stream from which to read
 * @param[in]        model	             Particle model defining interpretation of particle data
 * @param[in,out]    duplicates          Addresses of PAHs for use when reading primary particles
 *
 * @exception		 invalid_argument    Stream not ready
 */
void PAHPrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model, PahDeserialisationMap &pah_duplicates)
{
	if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        if (model.WriteBinaryTrees()) {
            // Call the binary tree serialiser...
            BinTreeSerializer <PAHPrimary> tree;
            tree.Deserialize(in, this, model, &pah_duplicates);
        } else {
            // Just deserialise the root node.
            DeserializePrimary(in, model, &pah_duplicates);
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PAHPrimary::Deserialize).");
    }
}

/*
 * @brief Reads an individual primary from the binary stream
 *
 * @param[in,out]	 in		             Input binary stream
 * @param[in]        model	             Particle model defining interpretation of particle data
 * @param[in,out]    duplicates          Addresses of PAHs for use when reading primary particles
 *
 * @exception		 invalid_argument    Stream not ready
 */
void PAHPrimary::DeserializePrimary(std::istream &in, const Sweep::ParticleModel &model, void *duplicates)
{
    if (in.good()) {

        int  val_int(0);
        double val(0.0);

        // Serialise state space
        in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
        m_numcarbon = val_int;

        in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
        m_numH = val_int;

        in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
        m_numOfEdgeC = val_int;

        in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
        m_numOfRings = val_int;

        in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
        m_numPAH = val_int;

        in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
        m_numprimary = val_int;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_PAHmass = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_PAHCollDiameter = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_primarydiam = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_children_radius = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_children_vol = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_leftparticle_vol_old = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_rightparticle_vol_old = val;

        in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
        m_rightparticle_numPAH = val_int;

        in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
        m_leftparticle_numPAH = val_int;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_children_surf = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_children_roundingLevel = val;
        /* imaging properies
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_Rg = val;
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_fdim = val;
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_sqrtLW = val;
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_LdivW = val;*/

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_avg_coalesc = val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_sint_time = val;

        in.read(reinterpret_cast<char*>(&val_int), sizeof(val_int));
        // the size of m_PAH
        const int PAHcount = val_int;
        // read the PAH stack (m_PAH)
        PahDeserialisationMap *pahDuplicates = reinterpret_cast<PahDeserialisationMap*>(duplicates);
        if (PAHcount != 0)
            inputPAHs(in, model, PAHcount, *pahDuplicates);

        m_leftchild = NULL;
        m_rightchild = NULL;
        m_leftparticle = NULL;
        m_rightparticle = NULL;
        m_parent = NULL;
        // Output base class.
        Primary::Deserialize(in, model);

    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PAHPrimary::DeserializePrimary).");
    }
}

/*
 * @brief Reads individual PAHs from the binary stream 
 *
 * @param[in]        in                Input binary stream
 * @param[in]        model	           Particle model defining interpretation of particle data
 * @param[in]        PAHcount          Number of PAHs in the particle
 * @param[in,out]    pah_duplicates    Addresses of PAHs for use when reading primary particles
 */
void PAHPrimary::inputPAHs(std::istream &in, const Sweep::ParticleModel &model, const int PAHcount, PahDeserialisationMap &pah_duplicates)
{
    int m_count=0;
    
    //PahDeserialisationMap pahPointerMappings;
    while (m_count != PAHcount)
    {
        // Raw memory location of PAH prior to serialisation.
        void* preSerializePointer;
        in.read(reinterpret_cast<char*>(&preSerializePointer), sizeof(preSerializePointer));
        //cout << preSerializePointer;

        // See if this PAH has already been deserialised
        //PahDeserialisationMap::iterator searchResult = pahPointerMappings.find(preSerializePointer);
        //if(searchResult != pahPointerMappings.end())
        //{
        //    // Found PAH in list that has already been deserialised
        //    m_PAH.push_back(searchResult->second);
        //}
        //else {
        //    // create a new PAH to load in its details
        //    boost::shared_ptr<PAH> new_PAH (new PAH());
        //    new_PAH->Deserialize(in);
        //    m_PAH.push_back(new_PAH);
        //    pahPointerMappings.insert(std::make_pair(preSerializePointer, new_PAH));
        //}
        PahDeserialisationMap::iterator searchResult = pah_duplicates.find(preSerializePointer);
        if(searchResult != pah_duplicates.end())
        {
            // Found PAH in list that has already been deserialised
            m_PAH.push_back(searchResult->second);
        }
        else {
            // create a new PAH to load in its details
            boost::shared_ptr<PAH> new_PAH (new PAH());
            new_PAH->Deserialize(in);
            m_PAH.push_back(new_PAH);
            pah_duplicates.insert(std::make_pair(preSerializePointer, new_PAH));
        }
        ++m_count;
    }
}

/*!
 * This is a helper function for ???????
 * It climbs up the tree from bottom to top recording a route
 * suitable for use in call to @see descendPath.
 *
 *@param[in]    bottom      Tree node from which to start climbing
 *@param[in]    top         Tree node at which to stop climbing
 *
 *@pre      top must be above bottom in a tree
 *
 *@return   Stack that can be used to descend the same path by moving to the
 *          left child each time the top of the stack is true.
 */
std::stack<bool> PAHPrimary::recordPath(const PAHPrimary* bottom,
                                        const PAHPrimary* const top) {
    std::stack<bool> wasLeftChild;

    while(bottom != top) {
        // check whether bottom was a left child of its parent
        wasLeftChild.push(bottom == bottom->m_parent->m_leftchild);

        // Climb one level up the tree
        bottom = bottom->m_parent;
    }
    return wasLeftChild;
}

/*!
 *@param[in]        here            Point in tree from which to start descent
 *@param[in,out]    takeLeftBranch  Instructions for which child to move to at each level
 *
 *@return   The node at the bottom of the path
 *
 *@pre  here must be a node of tree in which takeLeftBranch is a valid path
 *@post takeLeftBranch.empty() == true
 */
PAHPrimary* PAHPrimary::descendPath(PAHPrimary *here,
                                    std::stack<bool> &takeLeftBranch) {
    while(!takeLeftBranch.empty()) {
        // Move one step down the tree in the instructed direction
        if(takeLeftBranch.top())
            here = here->m_leftchild;
        else
            here = here->m_rightchild;

        // This instuction has now been processed
        takeLeftBranch.pop();
    }

    return here;
}

/*!
 * @brief       Sinters particles for time dt
 * 
 * This function only operates on non-leaf nodes. It begins at the root 
 * node, which sinters for time dt. It then descends the tree to sinter
 * nodes below the root. If the sintering level rises above 95%, Merge
 * is called and the particles are combined. 
 * 
 * 
 * @param[in]   dt      Time for which to sinter
 * @param[in]   sys     Environment for particles
 * @param[in]   model   Sintering model to apply
 * @param[in]   rng     Random number generator
 * @param[in]   wt      Statistical weight
 */
void PAHPrimary::Sinter(double dt, Cell &sys,
                            const Processes::SinteringModel &model,
                            rng_type &rng,
                            double wt)
{
    //Only update the time on the root node
    if (m_parent == NULL) {
        m_sint_time += dt;
        SetSinteringTime(m_sint_time);
    }

	//Do only if there is a particle to sinter
	if (m_leftparticle!=NULL)
    {
        // Store the old surface area of particles
        // double surf_old = m_children_surf;

        // Calculate the spherical surface
        const double spherical_surface=4*PI*m_children_radius*m_children_radius;
        // Declare time step variables.
        double t1=0.0, delt=0.0, tstop=dt;
        double r=0.0;

        // Define the maximum allowed change in surface
        // area in one internal time step (10% spherical surface).
        double dAmax = 0.1 * spherical_surface;

        // The scale parameter discretises the delta-S when using
        // the Poisson distribution.  This allows a smoother change
        // (smaller scale = higher precision).
        double scale = 0.01;

        // Perform integration loop.
        while (t1 < tstop)
        {
            // Calculate sintering rate.
            r = model.Rate(m_time+t1, sys, *this);

            if (r > 0) {
                // Calculate next time-step end point so that the
                // surface area changes by no more than dAmax.
                delt = dAmax / max(r, 1.0e-300);

                // Approximate sintering by a poisson process.  Calculate
                // number of poisson events.
                double mean;

                if (tstop > (t1+delt)) {
                    // A sub-step, we have changed surface by dAmax, on average
                    mean = 1.0 / scale;
                } else {
                    // Step until end.  Calculate degree of sintering explicitly.
                    mean = r * (tstop - t1) / (scale*dAmax);
                }
                boost::random::poisson_distribution<unsigned, double> repeatDistribution(mean);
                const unsigned n = repeatDistribution(rng);

                // Adjust the surface area.
                if (n > 0) {
                    m_children_surf -= (double)n * scale * dAmax;

                    // Check that primary is not completely sintered.
                    if (m_children_surf <= spherical_surface) {
                        m_children_surf = spherical_surface;
                        break;
                    }
                }
                // Set t1 for next time step.
                t1 += delt;
            }
        }
        m_children_roundingLevel=RoundingLevel();
        // one can specify a member in PAHPrimary class to store the rate of sintering, but now, it is not useful.
        //m_sint_rate = r;

        // Check if the sintering level is above the threshold, and merge if true
        if(m_children_roundingLevel>0.95)
          {
                CheckRounding();
                UpdateCache();
                if (m_leftchild!=NULL && m_rightchild!=NULL)
                {
                    m_leftchild->Sinter(dt,sys,model,rng,wt);
                    m_rightchild->Sinter(dt,sys,model,rng,wt);
                }
           }
         else
         {
             m_leftchild->Sinter(dt, sys, model,rng,wt);
             m_rightchild->Sinter(dt, sys, model,rng,wt);
         }

        UpdateCache();
        //// Adjust the gas-phase concentration
        //fvector dc(sys.GasPhase().Species()->size(), 0.0);

        //double n_NAvol_sint = wt * (double)num_H2O / (NA * sys.SampleVolume());
        //dc[Sprog::Species::Find(string("H2O"),*sys.GasPhase().Species())] += n_NAvol_sint;
        //sys.AdjustConcs(dc);
        m_children_roundingLevel=RoundingLevel();
    }  // endif m_leftparticle != NULL
}

void PAHPrimary::SetSinteringTime(double time) 
{
    m_sint_time = time;
    // Update children
    if (m_leftchild != NULL) {
        m_leftchild->SetSinteringTime(time);
        m_rightchild->SetSinteringTime(time);
    }
    // Update particles
    if (m_leftparticle != NULL) {
        m_leftparticle->SetSinteringTime(time);
        m_rightparticle->SetSinteringTime(time);
    }
}

//returns the CoalescenceLevel of the two primaries that are connected by this node
double PAHPrimary::RoundingLevel()
{
    if (m_leftparticle!=NULL) {
        // Calculate the spherical surface
        const double spherical_surface=4*PI*m_children_radius*m_children_radius;
        const double two_1_3=0.79370052231642452;
        double slevel;

        if (m_children_surf <= spherical_surface) {
            m_children_surf = spherical_surface;
            return 1.0;
        }

        if (m_children_surf == 0.0) {
            slevel = 0.0;
        } else {
            slevel= ((spherical_surface/m_children_surf)-two_1_3)/(1-two_1_3);
        }

        if (slevel < 0.0) {
            return 0.0;
        } else if (slevel > 1.0) {
            cout << "sweep: PAHPrimary::CoalescenceLevel() > 1.0";
            return 1.0;
        } else {
            return slevel;
        }
    } else {
        // Particle is a primary, should have 0 as default properties
        m_children_surf = 0.0;
        m_children_radius = 0.0;
        m_children_vol = 0.0;
        return 0.0;
    }
}
