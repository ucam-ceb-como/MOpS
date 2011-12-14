/*
  Author(s): Shraddha Shekar (ss663)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2010 Shraddha Shekar.


  File purpose:
	This file defines the state space of a primary in the silica system.


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
#include "swp_aggmodel_type.h"
#include "swp_model_factory.h"
#include "swp_silica_primary.h"
#include "swp_silica_cache.h"
#include "swp_particle_image.h"
#include "swp_cell.h"
#include "swp_kmc_pah_process.h"
#include "swp_kmc_pah_structure.h"
#include "swp_PAH.h"
#include "gpc_species.h"
#include <stdexcept>
#include <cassert>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_smallint.hpp>
#include "string_functions.h"

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace std;
using namespace Strings;


// CONSTRUCTORS AND DESTRUCTORS.
/*!
 * @brief       Initialising constructor with no arguments
 * 
 * Initialises a particle from no arguments. 
 *
 */
SilicaPrimary::SilicaPrimary() : Primary(),
    //State Space:: number of Si, O and OH units
    m_numSi(0),
    m_numO(0),
    m_numOH(0),
    m_numprimary(0),
    m_primarydiam(0.0),
    //Properties of children
    m_children_radius(0.0),
    m_children_vol(0.0),
    m_children_surf(0.0),
    m_children_sintering(0.0),
    //Sintering properties
    m_avg_sinter(0.0),
    m_sint_rate(0.0),
    //Imaging properties
    m_Rg(0.0),
    m_fdim(0.0),
    m_sqrtLW(0.0),
    m_LdivW(0.0),
    //Children are nodes holding pointers to other children and/or primary particles
    m_leftchild(NULL),
    m_rightchild(NULL),
    //Parent node is the top node of the tree
    m_parent(NULL),
    //Particles are leaf nodes containing primary particles
    m_leftparticle(NULL),
    m_rightparticle(NULL),
    m_sint_time(0.0)
{
}

/*!
 * @brief       Initialising constructor using time and model
 * 
 * Initialises a particle with the default chemical properties. Calls
 * UpdateCache() to calculate the derived properties after creation.
 * 
 * @param[in]   time        Time at which particle is being created
 * @param[in]   model       Model which defines the meaning of the primary
 *
 */
SilicaPrimary::SilicaPrimary(const real time, const Sweep::ParticleModel &model)
: Primary(time, model),
    //State Space:: number of Si, O and OH units
    m_numSi(0),
    m_numO(0),
    m_numOH(0),
    m_numprimary(0),
    m_primarydiam(0.0),
    //Properties of children
    m_children_radius(0.0),
    m_children_vol(0.0),
    m_children_surf(0.0),
    m_children_sintering(0.0),
    //Sintering properties
    m_avg_sinter(0.0),
    m_sint_rate(0.0),
    //Imaging properties
    m_Rg(0.0),
    m_fdim(0.0),
    m_sqrtLW(0.0),
    m_LdivW(0.0),
    //Children are nodes holding pointers to other children and/or primary particles
    m_leftchild(NULL),
    m_rightchild(NULL),
    //Parent node is the top node of the tree
    m_parent(NULL),
    //Particles are leaf nodes containing primary particles
    m_leftparticle(NULL),
    m_rightparticle(NULL),
    m_sint_time(0.0)
{
    // Other parts of the code check for a non-zero composition.
    m_comp[0]=1;
    // Create new particle with these specifications
	m_numSi = 2;
	m_numO = 1;
	m_numOH = 6;

	//Update the derived properties
    UpdateCache();
}


/*!
 * @brief       Initialising constructor using time, position and model
 * 
 * Initialises a particle with the default chemical properties. Calls
 * UpdateCache() to calculate the derived properties after creation.
 * 
 * @param[in]   time        Time at which particle is being created
 * @param[in]   position    Position at which particle is being created
 * @param[in]   model       Model which defines the meaning of the primary
 *
 */
SilicaPrimary::SilicaPrimary(const real time, const real position,
                       const Sweep::ParticleModel &model)
: Primary(time, model),
    //State Space:: number of Si, O and OH units
    m_numSi(0),
    m_numO(0),
    m_numOH(0),
    m_numprimary(0),
    m_primarydiam(0.0),
    //Properties of children
    m_children_radius(0.0),
    m_children_vol(0.0),
    m_children_surf(0.0),
    m_children_sintering(0.0),
    //Sintering properties
    m_avg_sinter(0.0),
    m_sint_rate(0.0),
    //Imaging properties
    m_Rg(0.0),
    m_fdim(0.0),
    m_sqrtLW(0.0),
    m_LdivW(0.0),
    //Children are nodes holding pointers to other children and/or primary particles
    m_leftchild(NULL),
    m_rightchild(NULL),
    //Parent node is the top node of the tree
    m_parent(NULL),
    //Particles are leaf nodes containing primary particles
    m_leftparticle(NULL),
    m_rightparticle(NULL),
    m_sint_time(0.0)
{
    // Other parts of the code check for a non-zero composition
    m_comp[0]=1;
	// Create new particle with these specifications
	m_numSi = 2;
	m_numO = 1;
	m_numOH = 6;

    //Update the other properties
    UpdateCache();
}


/*!
 * @brief       Initialising constructor using time, model and boolean
 * 
 * 
 * @param[in]   time        Time at which particle is being created
 * @param[in]   model       Model which defines the meaning of the primary
 * @param[in]   nosilica    Boolean.. unsure of purpose?
 *
 */
SilicaPrimary::SilicaPrimary(real time, const Sweep::ParticleModel &model, bool nosilica)
: Primary(time, model),
    //State Space:: number of Si, O and OH units
    m_numSi(0),
    m_numO(0),
    m_numOH(0),
    m_numprimary(0),
    m_primarydiam(0.0),
    //Properties of children
    m_children_radius(0.0),
    m_children_vol(0.0),
    m_children_surf(0.0),
    m_children_sintering(0.0),
    //Sintering properties
    m_avg_sinter(0.0),
    m_sint_rate(0.0),
    //Imaging properties
    m_Rg(0.0),
    m_fdim(0.0),
    m_sqrtLW(0.0),
    m_LdivW(0.0),
    //Children are nodes holding pointers to other children and/or primary particles
    m_leftchild(NULL),
    m_rightchild(NULL),
    //Parent node is the top node of the tree
    m_parent(NULL),
    //Particles are leaf nodes containing primary particles
    m_leftparticle(NULL),
    m_rightparticle(NULL),
    m_sint_time(0.0)
{
    m_comp[0]=1;

}


//! Copy constructor.
SilicaPrimary::SilicaPrimary(const SilicaPrimary &copy)
{
    *this = copy;
    if (copy.m_leftchild!=NULL)
    {
        CopyTree(&copy);
    }
}


//! Default destructor.
SilicaPrimary::~SilicaPrimary()
{
	delete m_leftchild;
	delete m_rightchild;

    releaseMem();

}


/*!
 * @brief       Recursively copy the tree for non-leaf nodes.
 *
 * @param[in] source Pointer to the primary to be copied
*/
void SilicaPrimary::CopyTree(const SilicaPrimary *source)
{
    //create the new left and right children with nothing in them
    m_leftchild = new SilicaPrimary(source->CreateTime(),*m_pmodel,false);
    m_rightchild = new SilicaPrimary(source->CreateTime(),*m_pmodel,false);

    // copy the properties such as the volume, surface area and list of constituent chemical units
	m_leftchild->CopyParts(source->m_leftchild);
	m_rightchild->CopyParts(source->m_rightchild);

    //set the pointers to the parent
	m_leftchild->m_parent=this;
	m_rightchild->m_parent=this;

	// the left and right particle are set further down in UpdateAllPointers
	// These are the pointers that specify which primary particles touch each
	// other in the aggregate structure.
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

SilicaPrimary &SilicaPrimary::operator=(const Primary &rhs)
{
    /*if (this != &rhs) {
       const AggModels::SilicaPrimary *Silicaprimary = NULL;
       Silicaprimary = dynamic_cast<const AggModels::SilicaPrimary*>(&rhs);
    //  UpdateAllPointers(silicaprimary,this);
    }*/

    return *this;
}



//! Stream-reading constructor.
SilicaPrimary::SilicaPrimary(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}


/*!
 * @brief       Copy state-space and derived properties from source
 * 
 * This function is like a limited assignment operator, except that the
 * children are not copied and the pointers to the particles may need
 * adjusting after this method has finished.
 *
 * @param[in] source Pointer to the primary to be copied
 */
void SilicaPrimary::CopyParts(const SilicaPrimary *source)
{
	SetCollDiameter(source->CollDiameter());
	SetMobDiameter(source->MobDiameter());
	SetSphDiameter(source->SphDiameter());
	//m_silicaCollDiameter=source->m_silicaCollDiameter;
	SetSurfaceArea(source->SurfaceArea());
	m_time=source->m_time;
	//m_silicamass=source->m_silicamass;
    m_leftchild=source->m_leftchild;
	m_rightchild=source->m_rightchild;
	m_leftparticle=source->m_leftparticle;
	m_rightparticle=source->m_rightparticle;
	m_parent=source->m_parent;
	SetMass(source->Mass());
	m_numSi=source->m_numSi;
	m_numO=source->m_numO;
	m_numOH=source->m_numOH;
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
    m_children_sintering=source->m_children_sintering;
    m_fdim=source->m_fdim;
    m_Rg=source->m_Rg;
    m_avg_sinter=source->m_avg_sinter;
    m_sint_rate=source->m_sint_rate;
    m_sint_time=source->m_sint_time;

}


/*!
 * @brief       Randomly selects a primary in the binary tree
 * 
 * Select a primary uniformly at random from this particle
 * and descend the aggregate tree to find the primary.
 * Note that most silicaPrimarys are nodes in a tree representing
 * connectivity within an aggregate, so it is necessary to
 * descend the tree to find a primary that really is a
 * primary.
 *
 * @param[in]   rng     Random number generator
 * 
 * @return      Pointer to an object representing a physical primary
 */
SilicaPrimary *SilicaPrimary::SelectRandomSubparticle(rng_type &rng)
{
    // We want to choose an integer uniformly on the range [0, m_numprimary - 1]
    boost::random::uniform_smallint<int> uniformDistrib(0, m_numprimary - 1);

    return SelectRandomSubparticleLoop(uniformDistrib(rng));

}
/*!
 * @brief       Helper function for SelectRandomSubparticle
 * 
 * @param[in] target The primary to be selected
 *
 * @return      Pointer to the next node down the tree on the path to target primary
*/
SilicaPrimary *SilicaPrimary::SelectRandomSubparticleLoop(int target)
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
 * @brief       Updates all pointers, ensuring the connectivity of the tree
 *                 is preserved during copying.
 * 
 * Each node contains two pointers (m_leftparticle and m_rightparticle)
 * to primary particles that are connected by this node
 * This function is used when the entire particle tree is duplicated.
 * It sets the pointers in the copied node (this), such that the connectivity
 * of the primary particles in this node is the same as in the original node.
 *
 * A large number of asserts are provided to check that the pointer network
 * describing the aggregate structure is valid.  This can be re-enabled if
 * there is suspicion of a bug,
 *
 * @param[in] original Pointer to the primary to be copied
*/
void SilicaPrimary::UpdateAllPointers(const SilicaPrimary *original)
{
    //the primary has no children => there are no left and right particles
    if (original->m_leftchild == NULL)
	{
        m_leftparticle=NULL;
		m_rightparticle=NULL;
     }
	else
	{
        // Find the route to m_leftparticle in the original tree
        std::stack<bool> route = recordPath(original->m_leftparticle, original);

        // Now follow the same route down the new tree to find the new left particle
        m_leftparticle = descendPath(this, route);

         // Find the route to m_rightparticle in the original tree
        route = recordPath(original->m_rightparticle, original);

        // Now follow the same route down the new tree to find the new right particle
        m_rightparticle = descendPath(this, route);


	}
}

/*!
 * @brief       Coagulates *this* and rhs particle
 * 
 * Called by coagulate, this function joins sp1 (this) and sp2 (rhs).
 * Basic properties are first calculated, followed by derived.
 * 
 * @param[in] rhs   Pointer to the particle to be coagulated with this
 * @param[in] rng   Random number generator
*/
SilicaPrimary &SilicaPrimary::Coagulate(const Primary &rhs, rng_type &rng)
{
    const SilicaPrimary *rhsparticle = NULL;
    rhsparticle = dynamic_cast<const AggModels::SilicaPrimary*>(&rhs);

    // Add the components.
    for (unsigned int i=0; i!=min(m_comp.size(),rhsparticle->m_comp.size()); ++i) {
        m_comp[i] += rhsparticle->m_comp[i];
    }

    // Add the tracker values.
    for (unsigned int i=0; i!=min(m_values.size(),rhsparticle->m_values.size()); ++i) {
        m_values[i] += rhsparticle->m_values[i];
    }

            //coagulation process
            SilicaPrimary *newleft = new SilicaPrimary;
		    SilicaPrimary *newright = new SilicaPrimary;
			SilicaPrimary copy_rhs(*rhsparticle);


            //Randomly select where to add the second particle
            boost::bernoulli_distribution<> bernoulliDistrib;
		    if (bernoulliDistrib(rng))
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
            m_children_sintering=0;
			UpdateCache();

			//select the primaries that are touching
		    this->m_leftparticle=m_leftchild->SelectRandomSubparticle(rng);
		    this->m_rightparticle=m_rightchild->SelectRandomSubparticle(rng);

		    // set the sintering times
		    this->SetSinteringTime(std::max(this->m_sint_time, rhsparticle->m_sint_time));
		    m_createt = max(this->m_createt, rhsparticle->m_createt);

			//initialise the variables used to calculate the sintering level
            m_children_vol=(m_leftparticle->m_vol+m_rightparticle->m_vol);
			m_children_surf=(m_leftparticle->m_surf+m_rightparticle->m_surf);
            m_children_radius=pow(3.0/(4.0*PI)*(m_children_vol),(ONE_THIRD));
            m_children_sintering=SinteringLevel();
			CheckSintering();

			//must set all the pointer to NULL otherwise the delete function
            //will also delete the children
	        copy_rhs.m_leftchild=NULL;
	        copy_rhs.m_rightchild=NULL;
	        copy_rhs.m_parent=NULL;
			copy_rhs.m_leftparticle=NULL;
            copy_rhs.m_rightparticle=NULL;

	return *this;
}


/*!
 * @brief       Sinters particles for time dt
 * 
 * This function only operates on non-leaf nodes. It begins at the root 
 * node, which sinters for time dt. It then descends the tree to sinter
 * nodes below the root. If the sintering level rises above 95%, Merge
 * is called and the particles are combined. 
 * 
 * It concludes by adjusting the gas-phase for the number of OH units
 * which are lost as H2O due to the change in surface area of the
 * sintering mechansim.
 * 
 * @param[in]   dt      Time for which to sinter
 * @param[in]   sys     Environment for particles
 * @param[in]   model   Sintering model to apply
 * @param[in]   rng     Random number generator
 * @param[in]   wt      Statistical weight
 */
void SilicaPrimary::Sinter(real dt, Cell &sys,
                            const Processes::SinteringModel &model,
                            rng_type &rng,
                            real wt)
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
        double surf_old = m_children_surf;
        int numOH_old = m_numOH;
        // First calculate the sintering rate

        // Calculate the spherical surface
        const double spherical_surface=4*PI*m_children_radius*m_children_radius;

        // Declare time step variables.
        real t1=0.0, delt=0.0, tstop=dt;
        real r=0.0;

        // Define the maximum allowed change in surface
        // area in one internal time step (10% spherical surface).
        real dAmax = 0.1 * spherical_surface;

        // The scale parameter discretises the delta-S when using
        // the Poisson distribution.  This allows a smoother change
        // (smaller scale = higher precision).
        real scale = 0.01;

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
                real mean;

                if (tstop > (t1+delt)) {
                    // A sub-step, we have changed surface by dAmax, on average
                    mean = 1.0 / scale;
                } else {
                    // Step until end.  Calculate degree of sintering explicitly.
                    mean = r * (tstop - t1) / (scale*dAmax);
                }
                boost::random::poisson_distribution<unsigned, real> repeatDistribution(mean);
                const unsigned n = repeatDistribution(rng);

                // Adjust the surface area.
                if (n > 0) {
                    m_children_surf -= (real)n * scale * dAmax;

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

        m_children_sintering=SinteringLevel();
        m_sint_rate = r;
        double rho_site = m_numOH/m_surf;

        // Adjust the units of OH and O due to release of water
        m_leftparticle->m_numOH -= int(0.5*rho_site*abs(m_children_surf - surf_old));
        m_rightparticle->m_numOH -= int(0.5*rho_site*abs(m_children_surf - surf_old));
        m_leftparticle->m_numO += int((0.5*rho_site*abs(m_children_surf - surf_old))/2);
        m_rightparticle->m_numO += int((0.5*rho_site*abs(m_children_surf - surf_old))/2);
        if(m_leftparticle->m_numOH < 0)
        {
            m_leftparticle->m_numOH = 0;
        }
        if(m_rightparticle->m_numOH < 0)
        {
            m_rightparticle->m_numOH = 0;
        }

        // Check if the sintering level is above the threshold, and merge if true
        if(m_children_sintering>0.95)
          {
                CheckSintering();
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

        // Adjust the gas-phase concentration
        fvector dc(sys.GasPhase().Species()->size(), 0.0);
        int num_H2O = int(abs(numOH_old - m_numOH)/2);

        real n_NAvol_sint = wt * (real)num_H2O / (NA * sys.SampleVolume());
        dc[Sprog::Species::Find(string("H2O"),*sys.GasPhase().Species())] += n_NAvol_sint;
        sys.AdjustConcs(dc);

        m_children_sintering=SinteringLevel();
    }  // endif m_leftparticle != NULL

}

void SilicaPrimary::SetSinteringTime(real time) {
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

/*!
 * @brief       Calculates the sintering level for particles connected
 *                  by this node
 * 
 * @return      Sintering level
 */
double SilicaPrimary::SinteringLevel()
{
    if (m_leftchild != NULL && m_rightchild != NULL) {
        // Calculate the spherical surface
        const double spherical_surface=4*PI*m_children_radius*m_children_radius;
        const double two_1_3=0.79370052231642452;
        double slevel;

        //Added by ss663
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
            cout << "sweep: SilicaPrimary::SinteringLevel() s.l. > 1.0";
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

/*!
 * @brief       Overload of the SetTime function for SilicaPrimary
 *
 * Sets the LDPA update time throughout the binary tree. This is important
 * as Merge can sometimes delete the originial root node, losing m_time
 * stored in m_primary. This potentially leads to longer dt LDPA action
 * times when UpdateParticle is called.
 *
 * @param t     LDPA update time
 */
void SilicaPrimary::SetTime(real t) {
    m_time = t;

    // Set LDPA time of children
    if (m_leftchild != NULL) {
        m_leftchild->SetTime(t);
        m_rightchild->SetTime(t);
    }

    // Set LDPA time of particles
    if (m_leftparticle != NULL) {
        m_leftparticle->SetTime(t);
        m_rightparticle->SetTime(t);
    }
}

//calculates the fractal dimension of the particle and stores it in m_fdim
//void SilicaPrimary::CalcFractalDimension()
//{
//	double create_time=this->CreateTime();
//	Sweep::Imaging::ParticleImage img;
//    // construct the particle by colliding the primary particles
//
//	img.constructSubParttree(this);
//
//	double L,W;
//    // calculate the length and the width of the particle
//    img.LengthWidth(L,W);
//    // calculate the radius of gyration
//    m_Rg=img.RadiusofGyration();
//	m_sqrtLW=sqrt(L*W);
//	m_LdivW=L/W;
//    m_Rg=m_Rg*1e-9;
//    m_fdim=log((double)m_numprimary)/log(2*m_Rg/(m_primarydiam/m_numprimary));
//    /*if (m_fdim>0 && m_fdim<3)
//	{
//       string filename;
//	   //cout<<"Creating 3d file";
//       filename=cstr(m_fdim)+".3d";
//       ofstream out;
//       out.open(filename.c_str());
//       img.Write3dout(out,0,0,0);
//       out.close();
//    }*/
//
//}


/*!
 * @brief       Sets all of the properties of the children to zero
 * 
 * Used after the children are coalesced
 */
void SilicaPrimary::ResetChildrenProperties()
{
            m_children_sintering=0.0;
            m_children_surf=0.0;
            m_children_vol=0.0;
            m_avg_sinter=0;
			m_sint_rate=0;
}


/*!
 * @brief       Merges the left and right particles
 * 
 * Considers two cases of particle structure. If the node only has two
 * primaries in its subtree, they are simply deleted and *this* becomes
 * the new primary. Otherwise, the subtrees are moved in the manner
 * explained in Markus Sander's thesis.
 * 
 * @return      Pointer to new merged particle
 */
SilicaPrimary &SilicaPrimary::Merge()
{

      //make sure this primary has children to merge
	  if(m_leftchild!=NULL)
      {
		if (m_leftchild==m_leftparticle && m_rightchild==m_rightparticle)
		{
            //this node has only two primaries in its subtree
            //it is possible that this node is not the root node and belongs to a bigger particle
            //copy the silicas of both children to the parent node
            // sum up the si, o and oh atoms
			m_numSi = m_leftparticle->m_numSi + m_rightparticle->m_numSi;
			m_numO = m_leftparticle->m_numO + m_rightparticle->m_numO;
			m_numOH = m_leftparticle->m_numOH + m_rightparticle->m_numOH;


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

            // Only update the cache on m_parent if the sintering level of m_parent if the sintering
            // level won't call a merge on the parent node. Otherwise, the *this* memory address could
            // be removed from the tree and segmentation faults will result!
			if(m_parent!=NULL) {
				if (m_parent->SinteringLevel() <= 0.95) {
					m_parent->UpdateCache();
				}
			}
		}


		else
		{
			if (m_leftchild->m_numprimary<m_rightchild->m_numprimary)
			{
				//append to left subtree because there are fewer primaries
                //this is only to keep the tree balanced
				SilicaPrimary *oldleftparticle=m_leftparticle;
                m_rightparticle->m_numSi = m_leftparticle->m_numSi + m_rightparticle->m_numSi;
			    m_rightparticle->m_numO = m_leftparticle->m_numO + m_rightparticle->m_numO;
		     	m_rightparticle->m_numOH = m_leftparticle->m_numOH + m_rightparticle->m_numOH;

				m_rightparticle->UpdatePrimary();

				//set the pointers from the leftprimary to the rightprimary
                //this will be the new bigger primary
				oldleftparticle->ChangePointer(oldleftparticle,m_rightparticle);
				m_rightparticle->ChangePointer(m_rightparticle,m_rightparticle);
				//m_rightparticle->GetAllParents(oldleftparticle);

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

				SilicaPrimary *oldleftchild=m_leftchild;
				SilicaPrimary *oldparent=m_parent;

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
				// Append to right subtree
				SilicaPrimary *oldrightparticle=m_rightparticle;
				//Add all mass to leftparticle
				m_leftparticle->m_numSi = m_leftparticle->m_numSi + m_rightparticle->m_numSi;
			    m_leftparticle->m_numO = m_leftparticle->m_numO + m_rightparticle->m_numO;
		     	m_leftparticle->m_numOH = m_leftparticle->m_numOH + m_rightparticle->m_numOH;

				m_leftparticle->UpdatePrimary();
                //All pointers to m_leftparticle now point to oldright particle
				oldrightparticle->ChangePointer(oldrightparticle,m_leftparticle);
				//oldrightparticle->DeleteParent(this);
                m_leftparticle->ChangePointer(m_leftparticle,m_leftparticle); //????
				//m_leftparticle->GetAllParents(oldrightparticle);
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
				SilicaPrimary *oldrightchild=m_rightchild;
				SilicaPrimary *oldparent=m_parent;

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


       }

		return *this;


}

/*!
 * @brief       Adjusts the particle after a surface rxn event
 * 
 * Analogous to the implementation in Primary. The function will 
 * however descend the tree to find a primary adjust. The surface area
 * and volume added to the particle is also rippled through the binary
 * tree.
 * 
 * @param[in]   dcomp   Vector storing changes in particle composition
 * @param[in]   dvalues Vector storing changes in gas-phase comp
 * @param[in]   rng     Random number generator
 * @param[in]   n       Number of times for adjustment
 */
unsigned int SilicaPrimary::Adjust(const fvector &dcomp,
		const fvector &dvalues, rng_type &rng, unsigned int n)

{
    //cout<<"Performing Surface Reaction\n";
	//Check if this is a leafnode (or a primary particle)
	//if(this->m_leftchild == NULL && this->m_rightchild == NULL)
	if(m_numprimary == 1)
	{
		unsigned int i = 0;

		double dV;
		double m_vol_old = m_vol;
		//double m_surf_old = m_surf;

		// Add the components.
		for (i=0; i!=min(m_comp.size(),dcomp.size()); ++i)
		{
			m_comp[i] += dcomp[i] * (real)n;
		}

		//Increase the number of Si, O and OH units
		m_numSi += 1 * n;
		m_numO += 1 * n;
		m_numOH += 2 * n;

		// Add the tracker values.
		for (i=0; i!=min(m_values.size(),dvalues.size()); ++i)
		{
			m_values[i] += dvalues[i] * (real)n;
		}

		// Update only the primary
		UpdatePrimary();

		dV = m_vol - m_vol_old;

		double ct=m_pmodel->Components(0)->CoalescThresh();

		// Surface change due to volume addition
		double dS=dV*ct/(m_diam/2.0);

		// Climb back-up the tree and update the surface area and
		// sintering of a particle
		this->UpdateParents(dS);

	}
	//Else this a non-leaf node (not a primary)
	else
	{
		// Generate random numbers
        boost::bernoulli_distribution<> leftRightChooser;
        // Select particle
		if(leftRightChooser(rng))
			return m_leftparticle->Adjust(dcomp, dvalues, rng, n);
		else
			return m_rightparticle->Adjust(dcomp, dvalues, rng, n);
	}

   	// Update property cache.
    UpdateCache();
	return n;

}

/*!
 * @brief       Updates the surface area and sintering level of all parents
 *
 * @param[in]   dS      Surface area increment to adjust area by
 */
void SilicaPrimary::UpdateParents(double dS) {
    if (m_parent != NULL) {
        m_parent->m_children_surf += dS;
        m_parent->m_children_sintering = m_parent->SinteringLevel();
        m_parent->UpdateCache();
    }
}

/*!
 * @brief       Adjusts the particle after an IntP event
 * 
 * Interparticle reactions have two OH units combining to release
 * a H2O molecule to the gas-phase. The leftover O atom is retained in
 * the particle phase.
 * 
 * @param[in]   dcomp   Vector storing changes in particle composition
 * @param[in]   dvalues Vector storing changes in gas-phase comp
 * @param[in]   rng     Random number generator
 * @param[in]   n       Number of times for adjustment
 */
unsigned int SilicaPrimary::AdjustIntPar(const fvector &dcomp,
		const fvector &dvalues, rng_type &rng, unsigned int n)

{
	//cout<<"Performing Inter-particle Reaction\n";
    if(m_numprimary == 1)
	{
		unsigned int i = 0;

		// Add the components.
		for (i=0; i!=min(m_comp.size(),dcomp.size()); ++i)
		{
			m_comp[i] += dcomp[i] * (real)n;
		}

		//Change the number of Si, O and OH units
		//if(m_numOH>0)
		//{
		m_numO += 1*n;
		m_numOH -= 2*n;
		//}
		if(m_numOH<0)
		{
			m_numOH = 0;
		}

		// Add the tracker values.
		for (i=0; i!=min(m_values.size(),dvalues.size()); ++i)
		{
			m_values[i] += dvalues[i] * (real)n;
		}
	}
	else
	{
		// Generate random numbers
        boost::bernoulli_distribution<> leftRightChooser;

        // Select particle
		if(leftRightChooser(rng))
			return m_leftparticle->AdjustIntPar(dcomp, dvalues, rng, n);
		else
			return m_rightparticle->AdjustIntPar(dcomp, dvalues, rng, n);
	}


		// Update property cache.
		UpdateCache();
		return n;

}

//! Returns the number of OH sites
int SilicaPrimary::GetSites() const

{
	return m_numOH;
}

//! Returns the sintering rate
real SilicaPrimary::GetSintRate() const
{
	real sint_rate = m_sint_rate;
	if(m_leftchild!=NULL)
	{
		sint_rate += (m_leftchild->GetSintRate() + m_rightchild->GetSintRate());
	}

	return sint_rate;
}

//! Returns the sintering time
real SilicaPrimary::GetSintTime() const {
    return m_sint_time;
}

//! Releases memory
void SilicaPrimary::ReleaseMem()
{
}

/*!
 * @brief       Changes pointer from source to target
 * 
 * @param[in] source Pointer to the original particle
 * @param[in] target Pointer to the new particle
*/
void SilicaPrimary::ChangePointer(SilicaPrimary *source, SilicaPrimary *target)
{
		if(m_rightparticle==source){
			m_rightparticle=target;
            double sphericalsurface=
                4*PI*pow(3*(m_leftparticle->Volume()+m_rightparticle->Volume())/(4*PI),TWO_THIRDS);
            m_children_surf=sphericalsurface/(m_children_sintering*0.2063+0.7937);    //sphericalsurface/(m_children_coalescence*(1-2^(-1/3))+2^(-1/3))
		}
		if(m_leftparticle==source){
			m_leftparticle=target;
            double sphericalsurface=
                4*PI*pow(3*(m_leftparticle->Volume()+m_rightparticle->Volume())/(4*PI),TWO_THIRDS);
            m_children_surf=sphericalsurface/(m_children_sintering*0.2063+0.7937);    //sphericalsurface/(m_children_coalescence*(1-2^(-1/3))+2^(-1/3))

		}

    // Update the tree above this sub-particle.
    if (m_parent != NULL) {
        m_parent->ChangePointer(source,target);
    }

}

//! UpdateCache helper function
void SilicaPrimary::UpdateCache(void)
{
    UpdateCache(this);
}

/*!
 * @brief       Checks the sintering level of the particle
 * 
 * If the sintering level is above 95%, Merge is called and the cache
 * is updated. 
 * 
 * @return      Boolean telling if the particle has sintered
 */
bool SilicaPrimary::CheckSintering()
{
	bool hassintered=false;
    if (m_children_sintering> 0.95 && m_leftparticle!=NULL)
        {

             Merge();
			 UpdateCache();
             hassintered=true;

             //check again because this node has changed
             CheckSintering();
        }
    if (m_leftchild!=NULL)
    {
        hassintered=m_leftchild->CheckSintering();
        hassintered=m_rightchild->CheckSintering();

    }

//	if(hassintered)
//		ResetVol();
    return hassintered;
}


/*!
 * @brief       Updates the properties of a primary
 */
void SilicaPrimary::UpdatePrimary(void)
{

    m_mass=(m_numSi*4.6621e-26 + m_numO*2.6565e-26 + m_numOH*2.8225e-26);  //convert to kg

    double silica_density = m_pmodel->Components(0)->Density(); // get density
    m_vol = m_mass / silica_density;							//in m^3
	m_diam = pow(6.0 * m_vol / PI, ONE_THIRD);
    m_dmob = m_diam;
	m_dcol=m_diam;
	m_surf = PI * m_diam * m_diam;
    m_primarydiam = m_diam;
	m_avg_sinter=0;
	m_numprimary = 1;
}


//! Set properties of particle to zero
void SilicaPrimary::Reset()
{
	m_numSi=0;
	m_numO=0;
	m_numOH=0;
	m_primarydiam=0.0;
	m_numprimary=0;
    m_surf=0;
    m_vol=0;
    m_avg_sinter=0;
}


/*!
 * @brief       Updates the SilicaPrimary cache
 * 
 * Updates the whole particle, from root to the lowest leaf-node.
 * Calculates the sintering level and merges particles (through Check
 * Sintering) if necessary. Only accurately calculates collision 
 * diameter and other properties used outside of SilicaPrimary for 
 * the root node.
 * 
 * @param[in] root The root node of this particle
*/
void SilicaPrimary::UpdateCache(SilicaPrimary *root)
{
    //Update the children
	if (m_leftchild!=NULL)
	{
		m_leftchild->UpdateCache(root);
		m_rightchild->UpdateCache(root);
	}
    //this is a primary and the number of primaries below this node is one (this node and no children)
	else
	{
            m_numprimary=1;
            m_avg_sinter=0;
			UpdatePrimary();
	}

    //this is not a primary, sum up the properties
	if (m_leftchild!=NULL)
    {

        Reset();
		m_numprimary=m_leftchild->m_numprimary+m_rightchild->m_numprimary;
		m_numSi=m_leftchild->m_numSi+m_rightchild->m_numSi;
		m_numO=m_leftchild->m_numO+m_rightchild->m_numO;
		m_numOH=m_leftchild->m_numOH+m_rightchild->m_numOH;
		m_surf = m_leftchild->m_surf+m_rightchild->m_surf;
		m_primarydiam = (m_leftchild->m_primarydiam+m_rightchild->m_primarydiam);
        m_vol=(m_leftchild->m_vol+m_rightchild->m_vol);
        m_mass=(m_leftchild->m_mass+m_rightchild->m_mass);


		// Calculate the sintering level of the two primaries connected by this node
		m_children_sintering=SinteringLevel();
		if (m_children_sintering>.95)
		{
			CheckSintering();
		}

		// Sum up the avg sintering level
		if((m_leftchild != NULL) && (m_rightchild != NULL))
		{
			m_avg_sinter=m_children_sintering+m_leftchild->m_avg_sinter+m_rightchild->m_avg_sinter;
		}
		else
		{
			m_avg_sinter=m_children_sintering;
		}

        // calculate the different diameters only for the root node because this goes into the
        // particle tree and gets used by the coagulation kernel

		if (this==root)
        {
             // Get spherical equivalent radius and diameter
            double spherical_radius=pow(3*m_vol/(4*PI),ONE_THIRD);
            m_diam=2*spherical_radius;

			// there are m_numprimary-1 connections between the primary particles
            m_avg_sinter=m_avg_sinter/(m_numprimary-1);

			/*if (m_avg_sinter == 0)
			{
				m_avg_sinter = 1e-5;
			}*/

			// Approxmiate the surface of the particle
            const real numprim_1_3=pow(m_numprimary,-0.333333);
            m_surf=4*PI*spherical_radius*spherical_radius/(m_avg_sinter*(1-numprim_1_3)+numprim_1_3);

            // Calculate dcol based-on formula given in Lavvas et al. (2011)
            // assume fractal dimension Df = 1.8
			const double aggcolldiam=(6*m_vol/m_surf)*pow(pow(m_surf,3)/(36*PI*m_vol*m_vol),(1.0/1.8));
			m_dmob = aggcolldiam;
            SetCollDiameter(aggcolldiam);

        }
        else
        {
            m_diam=0;
            m_dmob=0;

        }
    }
}

/*!
 * @brief       Prints a graphical output of the binary tree structure
 * 
 * Graph outputted in GraphViz format. Use command 'dot' to convert.
 * e.g. dot -Tpng -o name.png name.dat
 * for PNG terminal, plotting output file name.dat
 * 
 * @param[in] filename Output filename
*/
void SilicaPrimary::PrintTree(string filename)
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
 * @brief       Helper function for printing tree.
 * 
 * @param[in] out Output stream
*/
void SilicaPrimary::PrintTreeLoop(std::ostream &out)
{
    // Non leaf-node case
    if (m_leftchild!=NULL)
    {
        out << "\" " << this << "\" " << " [shape = \"record\" label = \"";
        this->PrintTreeNode(out);
        out <<"\"];"<<endl;

        out << "\" " << this->m_leftchild << "\" " << " [shape = \"record\" label = \"";
        this->PrintTreeNode(out);
        out << "|" << this << "\"];"<<endl;

        out << "\" " << this->m_rightchild << "\" " << " [shape = \"record\" label = \"";
        this->PrintTreeNode(out);
        out << "|" << this << "\"];"<<endl;

        out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftchild<<"\"; "<<endl;
        out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightchild<<"\"; "<<endl;
        out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_leftparticle<<"\"[label=\""<<this<<"\",color=\"blue\"]; "<<endl;
        out<<"\" "<<this<<"\" "<<"->"<<"\" "<<this->m_rightparticle<<"\"[label=\""<<this<<"\",color=\"blue\"]; "<<endl;
        m_leftchild->PrintTreeLoop(out);
        m_rightchild->PrintTreeLoop(out);
    }

    // Case when the node is a primary
    else
    {
        out << "\" " << this << "\" " << " [shape = \"record\" color=\"blue\" label = \"";
        this->PrintTreeNode(out);
        out <<"\"];"<<endl;
    }
}

/*!
 * @brief       Prints out each node of the tree
 *
 * @param[in] out Output stream
*/
void SilicaPrimary::PrintTreeNode(std::ostream &out) {
    out << "m_surf="             << this->m_surf
            << "|m_child_surf="      << this->m_children_surf
            << "|m_vol="             << this->m_vol
            << "|m_child_sint="      << this->m_children_sintering
            << "|m_numSi="           << this->m_numSi
            << "|m_numO="            << this->m_numO
            << "|m_numOH="           << this->m_numOH
            << "|m_time="            << this->m_time
            << "|m_numprimary="      << this->m_numprimary
            << "|m_child_rad="       << this->m_children_radius
            << "|m_parent="          << this->m_parent
            << "|this=" << this;
}


//! Creates an aggregation data cache for this primary type.
AggModels::SilicaCache *const SilicaPrimary::CreateAggCache() const
{
	SilicaCache *cache =
        static_cast<SilicaCache*>(ModelFactory::CreateAggCache(AggModels::Silica_ID));
    if (cache != NULL) *cache = *this;
    return cache;
}

/*!
 * @brief Returns the right child.
 * @return  Right child pointer
 */
const SilicaPrimary *SilicaPrimary::RightChild() const
{
   return m_rightchild;
}

/*!
 * @brief Returns the left child.
 * @return  Left child pointer
 */
const SilicaPrimary *SilicaPrimary::LeftChild() const
{
   return m_leftchild;
}

double SilicaPrimary::Rg() const
{
   return m_Rg;
}

double SilicaPrimary::Fdim() const
{
   return m_fdim;
}

double SilicaPrimary::PrimaryDiam() const
{
   return m_primarydiam;
}

double SilicaPrimary::LdivW() const
{
   return m_LdivW;
}

int SilicaPrimary::Numprimary() const
{
	return m_numprimary;
}

int SilicaPrimary::NumSi() const
{
   return m_numSi;
}

int SilicaPrimary::NumO() const
{
   return m_numO;
}

int SilicaPrimary::NumOH() const
{
   return m_numOH;
}


double SilicaPrimary::sqrtLW() const
{
   return m_sqrtLW;
}

double SilicaPrimary::AvgSinter() const
{
   return m_avg_sinter;
}

// READ/WRITE/COPY.

// Returns a copy of the model data.
SilicaPrimary *const SilicaPrimary::Clone(void) const
{
   return new SilicaPrimary(*this);
}


// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType SilicaPrimary::AggID(void) const {return AggModels::Silica_ID;}


// Writes the object to a binary stream.
void SilicaPrimary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

		double val = 0.0;
		// Write number of Si atoms.
        val = (int) m_numSi;
        out.write((char*)&val, sizeof(val));

		// Write number of O atoms.
        val = (int) m_numO;
        out.write((char*)&val, sizeof(val));

		// Write number of OH atoms.
        val = (int) m_numOH;
        out.write((char*)&val, sizeof(val));

	    val = (int)m_numprimary;
        out.write((char*)&val, sizeof(val));

	    val = (double)m_sqrtLW;
        out.write((char*)&val, sizeof(val));

	    val = (double)m_LdivW;
        out.write((char*)&val, sizeof(val));

	    val = (double)m_primarydiam;
        out.write((char*)&val, sizeof(val));

        val = (double)m_fdim;
        out.write((char*)&val, sizeof(val));

        val = (double)m_Rg;
        out.write((char*)&val, sizeof(val));


        // Output base class.
        Primary::Serialize(out);

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, SilicaPrimary::Serialize).");
    }
}

// Reads the object from a binary stream.
void SilicaPrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    //UpdateCache();
	if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

		double val = 0.0;
		// Read number of Si atoms.
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_numSi = (int)val;

		// Read number of O atoms.
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_numO = (int)val;

		// Read number of OH atoms.
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_numOH = (int)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_numprimary = (int)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_sqrtLW = (real)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_LdivW = (real)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_primarydiam = (real)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_fdim = (real)val;

		in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_Rg = (real)val;

		m_leftchild=NULL;
		m_rightchild=NULL;
		m_parent=NULL;
		m_leftparticle=NULL;
		m_rightparticle=NULL;
		//m_allparents.clear();

		switch (version) {
            case 0:
                // Read base class.
                Primary::Deserialize(in, model);

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                  "(Sweep, SilicaPrimary::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SilicaPrimary::Deserialize).");
    }
}


/*!
 *
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
std::stack<bool> SilicaPrimary::recordPath(const SilicaPrimary* bottom,
                                        const SilicaPrimary* const top)
{
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
SilicaPrimary* SilicaPrimary::descendPath(SilicaPrimary *here,
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
