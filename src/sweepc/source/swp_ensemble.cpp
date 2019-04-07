/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Ensemble class declared in the
    swp_ensemble.h header file.

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

#include "swp_ensemble.h"
#include "swp_particle_model.h"
#include "swp_kmc_simulator.h"
#include "swp_PAH_primary.h"
#include "swp_kmc_pah_structure.h"
#include "swp_kmc_structure_comp.h"
#include "swp_kmc_typedef.h"
#include "swp_kmc_pah_process.h"

#include <cmath>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <limits>
#include <functional>
#include <algorithm>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Sweep::Ensemble::Ensemble(void)
{
	m_kmcsimulator= NULL;
    init();
}

// Initialising constructor.
Sweep::Ensemble::Ensemble(unsigned int count)
: m_tree(count)
{
    // Call initialisation routine.
    //If there are no particles, do not initialise binary tree
    //Instead, initialise the ensemble to zero and skip the binary tree
    if (count ==0)
        init();

    else
    //Initialise binary tree, which also initialises the ensemble
        Initialise(count);
}

// Copy contructor.
Sweep::Ensemble::Ensemble(const Sweep::Ensemble &copy)
:m_tree(copy.m_tree)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
Sweep::Ensemble::Ensemble(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Destructor.
Sweep::Ensemble::~Ensemble(void)
{
	delete     m_kmcsimulator;
    // Clear the ensemble.
    Clear();
}


// OPERATOR OVERLOADING.

// Assignment operator.
Ensemble & Sweep::Ensemble::operator=(const Sweep::Ensemble &rhs)
{
    if (this != &rhs) {
        // Clear current particles.
        Clear();

        if (rhs.m_capacity > 0) {
            // Resize particle vector.
            m_particles.resize(rhs.m_capacity, NULL);

            // Capacity.
            m_levels   = rhs.m_levels;
            m_capacity = rhs.m_capacity;
            m_halfcap  = rhs.m_halfcap;
            m_count    = rhs.m_count;
            //m_numofInceptedPAH = rhs.m_numofInceptedPAH;
            // Scaling.
            m_ncont      = rhs.m_ncont;
            m_contfactor = rhs.m_contfactor;
            m_contwarn   = rhs.m_contwarn;
            // Doubling.
            m_maxcount   = rhs.m_maxcount;
            m_ndble      = rhs.m_ndble;
            m_dbleactive = rhs.m_dbleactive;
            m_dblecutoff = rhs.m_dblecutoff;
            m_dblelimit  = rhs.m_dblelimit;
            m_dbleslack  = rhs.m_dbleslack;
            m_dbleon     = rhs.m_dbleon;

            // Copy particle vector.
            for (unsigned int i=0; i!=rhs.Count(); ++i) {
                m_particles[i] = rhs.m_particles[i]->Clone();
            }

			//csl37
			//copy tracked particles
			m_tracked_number = rhs.m_tracked_number;
			m_tracked_particles.resize(rhs.m_tracked_number, NULL);
			for (unsigned int i=0; i!=rhs.m_tracked_number; ++i) {
                m_tracked_particles[i] = rhs.m_tracked_particles[i]->Clone();
            }

            m_tree.resize(m_capacity);
            rebuildTree();
        }
    }
    return *this;
}


/*!
 * @param os    Output stream
 * @param net   Ensemble object to print
 * @return      Output stream
 */
std::ostream& Sweep::operator<<(
        std::ostream &os,
        const Sweep::Ensemble &e)
{
  os << "[Ensemble],";
  os << " nmax=" << e.Capacity();
  os << " n(t)=" << e.Count();
  os << "\n";
  return os;
}

// INITIALISATION.

/*!
 * @param[in]   capacity    New size of ensemble
 *
 * Any particles in the ensemble will be destroyed
 */
void Sweep::Ensemble::Initialise(unsigned int capacity)
{
    // Clear current ensemble.
    Clear();

    //Check that there is no ensemble with zero capacity
    if (capacity == 0)
        throw std::invalid_argument("Cannot create a binary tree for an ensemble with zero capacity! (Sweep::Ensemble::Initialise)");

    // Calculate nearest power of 2 capacity.  Ensemble capacity must be a power
    // of 2.  This constraint is due to the binary tree implementation.
    double rl    = log((double)capacity) / log(2.0);
    m_levels   = static_cast<int>(rl + 0.5);

    // Capacity is 2^levels, which is most efficiently calculated with a bit shift.
    // Note that 1 has the binary representation (8 bits only for brevity) 00000001.
    // The left shift used here moves bits left the specified numer of places,
    // while filling the right hand places with 0.  Bits at the left hand end of
    // the number are discarded when shifted out of the range of the data type.
    // For example 5 << 3 or with 5 in binary notation 00000101 << 3 == 00101000 == 40 == 5 * 2^3.
    // In particular 1 << n == 2^n.
    m_capacity = 1 << m_levels;

    m_halfcap  = m_capacity / 2;

    m_count    = 0;
    //m_numofInceptedPAH = 0;

    // Reserve memory for ensemble.
    m_particles.resize(m_capacity, NULL);

    m_tree.resize(m_capacity);

    // Initialise scaling.
    m_ncont      = 0;
    m_contfactor = (double)(m_capacity) / (double)(m_capacity+1);
    m_contwarn   = false;

    // Initialise doubling.
    m_maxcount   = 0;
    m_ndble      = 0;
    m_dbleon     = true;
    m_dbleactive = false;
    //m_dblecutoff = (int)(3.0 * (double)m_capacity / 4.0 / 4.0);
	m_dblecutoff = (int)(3.0 * (double)m_capacity / 4.0 );
	//m_dblecutoff = m_capacity+10;

	m_dbleslack = (unsigned int)pow(2.0, (int)((m_levels - 5)>0 ? m_levels - 5 : 0));
	//m_dbleslack = (int) m_capacity*0.075;
	//m_dbleslack = m_dbleslack / 4.0;
	//m_dblelimit = m_halfcap/4.0 - m_dbleslack;
	m_dblelimit = m_halfcap - m_dbleslack;

	//csl37: for testing purposes 
	//in future initialise such that the number of tracked particles can be modified
	m_tracked_number = 0; //csl37-no tracking
	m_tracked_particles.resize(m_tracked_number, NULL);
}

/*!
 * Turn on (true) or off (false) doubling.
 *
 * @param val   Switch for doubling.
 */
void Sweep::Ensemble::SetDoubling(const bool val) {
    m_dbleon = val;
}

unsigned int Sweep::Ensemble::DoubleLimit() {
	return m_dblelimit;
}

bool Sweep::Ensemble::IsDoublingOn() {
	return m_dbleactive;
}

/**
 * Initialise the ensemble to hold particles of the type specified
 * by the model and containing the particular particles contained
 * in the range [first, last).  This is equivalent to multiple applications
 * of Add on an Ensemble instance after a call to Initialise(m_capacity).
 *
 *@param[in]        first      Iterator to first in range of particle pointers to insert
 *@param[in]        last       Iterator to one past end of range of particle pointers to insert
 *@param[in,out]    rng        Random number generator
 */
void Sweep::Ensemble::SetParticles(std::list<Particle*>::iterator first, std::list<Particle*>::iterator last,
                                   rng_type &rng)
{
    // Clear any existing particles
    for(iterator it = m_particles.begin(); it != m_particles.end(); ++it) {
        delete *it;
    }
    m_particles.assign(m_capacity, NULL);

    unsigned count = 0;
    // Read up to m_capacity particles straight into the array
    while((first != last) && (count < m_capacity)) {
        m_particles[count++] = (*first++);
    }

    // Now we have to decide whether or not to accept particles
    while(first != last) {
        // Possible index in which to store this particle
        boost::uniform_smallint<unsigned> indexGenerator(0, count);
        const unsigned possibleIndex = indexGenerator(rng);

        // Accept the index with probability m_capacity / count
        if(possibleIndex < m_capacity) {
            delete m_particles[possibleIndex];
            m_particles[possibleIndex] = *first;
        }
        else {
            delete *first;
        }
        ++count;
        ++first;
    }

    if(count > m_capacity) {
        // Some particles were thrown away and we must rescale
        m_count = m_capacity;
        m_ncont = 0;

        iterator it = begin();
        const iterator itEnd = end();
        while(it != itEnd) {
            (*it)->setStatisticalWeight((*it)->getStatisticalWeight() * static_cast<double>(count) / static_cast<double>(m_capacity));
            ++it;
        }
    }
    else {
        m_count = count;
        m_ncont = 0;
    }
    m_maxcount = m_count;

    //std::cout << m_count << " particles set on ensemble of capacity " << m_capacity << '\n';

    // Initialise scaling.

    m_contwarn   = false;

    // Initialise doubling.
    m_ndble      = 0;
    m_dbleon     = true;

    // Check for doubling activation.
    if (!m_dbleactive && (m_count >= m_dblecutoff-1)) {
        m_dbleactive = true;
    } else
        m_dbleactive = false;

	//csl37: add tracked particles -- we can remove this now that we a function to initialise tracking
	int i = 0;
	while(i < m_tracked_number && i<m_count){
		m_tracked_particles[i] = m_particles[i];
		i++;
	}

    // Build the tree with the weights for the new particles.
    rebuildTree();

    assert(m_tree.size() == m_count);
}

Sweep::KMC_ARS::KMCSimulator* Sweep::Ensemble::Simulator(void)
{   
	return m_kmcsimulator;
}

void Sweep::Ensemble::SetSimulator(Sweep::GasProfile& gp)
{   
    Sweep::KMC_ARS::KMCSimulator* kmc = new Sweep::KMC_ARS::KMCSimulator(gp);
    m_kmcsimulator= kmc;
    m_kmcsimulator->TestGP();
}



// PARTICLE ADDITION AND REMOVAL.

// Returns a pointer to the particle with the given index.
Particle *const Sweep::Ensemble::At(unsigned int i)
{
    // Check that the index in within range, then return the particle.
    if (i < m_count) {
        return m_particles[i];
    } else {
        return NULL;
    }
}

// Returns a pointer to the particle with the given index.
const Particle *const Sweep::Ensemble::At(unsigned int i) const
{
    // Check that the index in within range, then return the particle.
    if (i < m_count) {
        return m_particles[i];
    } else {
        return NULL;
    }
}

/*!
 * @param[in,out]   sp          Particle to add to the ensemble
 * @param[in,out]   rng         Random number generator
 *
 * @return      Index of added particle
 *
 * Particles must be heap allocated, because the ensemble takes the
 * address of sp and, by default, calls delete on the resulting pointer
 * when the particle is no longer needed.  It would probably make
 * sense to change the type of the first argument to Particle* to
 * reflect the fact that the ensemble mainly deals with the pointer.
 *
 */
int Sweep::Ensemble::Add(Particle &sp, rng_type &rng)
{
    // Check for doubling activation.
    if (!m_dbleactive && (m_count >= m_dblecutoff-1)) {
        m_dbleactive = true;
        printf("sweep: Particle doubling activated.\n");
    }

    // Check ensemble for space, if there is not enough space then need
    // to generate some by contracting the ensemble.
    int i = -1;
    if (m_count < m_capacity) {
        // There is space in the tree for a new particle.
        i = -1;
    } else {
        // We must contract the ensemble to accommodate a new particle.
        boost::uniform_smallint<int> indexDistrib(0, m_capacity);
        boost::variate_generator<Sweep::rng_type&, boost::uniform_smallint<int> > indexGenerator(rng, indexDistrib);
        i = indexGenerator();

        ++m_ncont;
        if (!m_contwarn && ((double)(m_ncont)/(double)m_capacity > 0.01)) {
            m_contwarn = true;
            printf("sweep: Ensemble contracting too often; "
                   "possible stiffness issue.\n");
        }
    }

    if (i < 0) {
        // We are adding a new particle.
        i=m_count++;
        m_particles[i] = &sp;
        m_tree.push_back(tree_type::value_type(sp, m_particles.begin() + i));
        //m_numofInceptedPAH++;

		//csl37:ToDo check if need to add to tracked particle here

    } else if ((unsigned)i < m_capacity) {
        // Replace an existing particle (if i=m_capacity) then
        // we are removing the new particle, so just ignore it.
        Replace(i, sp);

		//csl37:ToDo replace function will remove tracked particle
    } else {
        // The new particle is to be removed immediately
        assert(static_cast<unsigned int>(i) == m_capacity);
        delete &sp;
    }

    m_maxcount = std::max(m_maxcount, m_count);

    assert(m_tree.size() == m_count);

    return i;
}

int Sweep::Ensemble::CheckforPAH(Sweep::KMC_ARS::PAHStructure &m_PAH, double t, int ind)
{
	iterator it1;
	int count = 0;
	for (it1 = begin(); it1 != end(); it1++){
		//First, check if this particle is updated to the correct time
		if ((*it1)->LastUpdateTime() == t && count > ind){
			AggModels::PAHPrimary *pah =
				dynamic_cast<AggModels::PAHPrimary*>((*it1)->Primary());
			//Check if this particle contains a single primary with a single PAH with the same 
			//amount of hydrogens and carbons
			if (pah->NumPAH() == 1){
				if (pah->NumCarbon() == m_PAH.numofC() && pah->NumHydrogen() == m_PAH.numofH() 
					&& pah->NumRings() == m_PAH.numofRings()
					&& pah->NumRings5() == m_PAH.numofRings5()){

					std::map<KMC_ARS::kmcSiteType, KMC_ARS::svector> sitemapInput = m_PAH.GetSiteMap();
					std::map<KMC_ARS::kmcSiteType, KMC_ARS::svector> sitemapComp = (*(pah->GetPAHVector())[0]).GetPAHStruct()->GetSiteMap();
					if (sitemapInput[KMC_ARS::FE].size() == sitemapComp[KMC_ARS::FE].size() &&
						sitemapInput[KMC_ARS::ZZ].size() == sitemapComp[KMC_ARS::ZZ].size() &&
						sitemapInput[KMC_ARS::AC].size() == sitemapComp[KMC_ARS::AC].size() &&
						sitemapInput[KMC_ARS::BY6].size() == sitemapComp[KMC_ARS::BY6].size() &&
						sitemapInput[KMC_ARS::BY5].size() == sitemapComp[KMC_ARS::BY5].size() &&
						sitemapInput[KMC_ARS::R5].size() == sitemapComp[KMC_ARS::R5].size() &&
						sitemapInput[KMC_ARS::RFE].size() == sitemapComp[KMC_ARS::RFE].size() &&
						sitemapInput[KMC_ARS::RZZ].size() == sitemapComp[KMC_ARS::RZZ].size() &&
						sitemapInput[KMC_ARS::RAC].size() == sitemapComp[KMC_ARS::RAC].size() &&
						sitemapInput[KMC_ARS::RBY5].size() == sitemapComp[KMC_ARS::RBY5].size() &&
						sitemapInput[KMC_ARS::RFER].size() == sitemapComp[KMC_ARS::RFER].size() &&
						sitemapInput[KMC_ARS::RZZR].size() == sitemapComp[KMC_ARS::RZZR].size() &&
						sitemapInput[KMC_ARS::RACR].size() == sitemapComp[KMC_ARS::RACR].size() &&
						sitemapInput[KMC_ARS::FE3].size() == sitemapComp[KMC_ARS::FE3].size() &&
						sitemapInput[KMC_ARS::FE2].size() == sitemapComp[KMC_ARS::FE2].size() &&
						sitemapInput[KMC_ARS::AC_FE3].size() == sitemapComp[KMC_ARS::AC_FE3].size() &&
						sitemapInput[KMC_ARS::BY5_FE3].size() == sitemapComp[KMC_ARS::BY5_FE3].size() &&
						sitemapInput[KMC_ARS::FE_HACA].size() == sitemapComp[KMC_ARS::FE_HACA].size()){
						return count;
					}
				}
			}
		}
		count++;
	}
	return -1;
}

/*!
 * @param[in]   i       Index of particle to remove
 * @param[in]   fdel    True if delete should be called on removed particle
 *
 * Calls to remove may invalidate some or all iterators and indices
 * referring to particles in the ensemble.
 */
void Sweep::Ensemble::Remove(unsigned int i, bool fdel)
{
    //if (m_particles[i]->Primary()->AggID() ==AggModels::PAH_KMC_ID)
    //    {
    //        const Sweep::AggModels::PAHPrimary *rhsparticle = NULL;
    //        rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(m_particles[i]->Primary());
    //        if (rhsparticle->Pyrene()!=0)
    //        {
    //            --m_numofInceptedPAH;
    //        }
    //    }
    //SetNumOfInceptedPAH(-1,m_particles[i]->Primary());
    
	//csl37: for only one particle, set to NULL
	//expand to multiple
//	if(m_tracked_particles[0] == m_particles[i] && i<m_count) m_tracked_particles[0] = NULL;	//csl37-no tracking
//	if(m_tracked_particles[0] == m_particles[i] && i<m_count) m_tracked_particles[0] = NULL; //csl37-ttip

	// See if IWDSA is being used. If so, do not attempt doubling at the end of this routine.
	bool doubling = true;
	if (m_particles[i]->Primary()->AggID() == AggModels::PAH_KMC_ID){
		if (m_particles[i]->Primary()->ParticleModel()->Components(0)->WeightedPAHs()){
			doubling = false;
		}
	}
	
	// Check that particle index is valid.
    if (i<m_count-1) {
        // First delete particle from memory, then
        // overwrite it with the last particle, which
        // is subsequently removed from the vector.
        if (fdel) delete m_particles[i];
        --m_count;
        m_particles[i] = m_particles[m_count];
        m_particles[m_count] = NULL;

        // Iterator to the particle that is being removed
        iterator itPart = m_particles.begin() + i;
        m_tree.replace(m_tree.begin() + i, tree_type::value_type(**itPart, itPart));
        m_tree.pop_back();

    } else if (i==m_count-1) {
        // This is the last particle in the ensemble, we don't
        // need to swap it with another, just delete it.

        // Erase particle.
        if (fdel) delete m_particles[i];
        m_particles[i] = NULL;
        --m_count;

        m_tree.pop_back();
    }

    // Particle removal might reduce the particle count
    // sufficiently to require particle doubling.
	// But only do so if not using IWDSA as otherwise ensemble is likely to overflow during next LPDA
	if(doubling) dble();

    assert(m_tree.size() == m_count);
}

// Removes invalid particles from the ensemble.
void Sweep::Ensemble::RemoveInvalids(void)
{
    // This function loops forward through the list finding invalid
    // particles and backwards finding valid particles.  Once an invalid
    // and a valid particle are found they are swapped.  This results in
    // all the invalid particles being collected at the end of the vector,
    // from where they can be deleted.

    // Rearrange the particle list in m_particles so that all the valid particles
    // are in the range [m_particles.begin(), validEnd) and all the invalid
    // particles are in the range [validEnd, m_particles.end()).
    iterator validEnd = std::partition(m_particles.begin(),
                                       m_particles.begin() + m_count,
                                       std::mem_fun(&Particle::IsValid));

    // Update the number of particles in the tree
    m_count = validEnd - m_particles.begin();

    // Now delete the invalid particles and nullify the corresponding pointers
    while(validEnd != m_particles.end()) {

		//csl37: for only one particle, set to NULL
		//expand to multiple
		//if(m_tracked_particles[0] == *validEnd) m_tracked_particles[0] = NULL;	//csl37- no tracking
//		if(m_tracked_particles[0] == *validEnd) m_tracked_particles[0] = NULL;	//csl37-ttip

        delete *validEnd;
        *validEnd = NULL;
        ++validEnd;
    }


    // Rebuild the binary tree structure
    rebuildTree();

    // Stop doubling because the number of particles has dropped from above
    // m_dblelimit during this function, which means a rapid loss of particles
    // so doubling will make the sample volume needlessly large.
    if(m_count < m_capacity - m_dblecutoff) {
        m_dbleactive = false;
    }

    // If we removed too many invalid particles then we'll have to double.
    dble();
    assert(m_tree.size() == m_count);
}

/*!
 * @param[in]       i       Index of particle to replace
 * @param[in,out]   sp      New particle to use
 *
 * Particles must be heap allocated, because the ensemble takes the
 * address of sp and, by default, calls delete on the resulting pointer
 * when the particle is no longer needed.  It would probably make
 * sense to change the type of the second argument to Particle* to
 * reflect the fact that the ensemble mainly deals with the pointer.
 *
 */
void Sweep::Ensemble::Replace(unsigned int i, Particle &sp)
{
    // if (m_particles[i]->Primary()->AggID() ==AggModels::PAH_KMC_ID)
    //{
    //    const Sweep::AggModels::PAHPrimary *rhsparticle = NULL;
    //    rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(m_particles[i]->Primary());
    //    if (rhsparticle->Pyrene()!=0)
    //    {
    //        m_numofInceptedPAH--;
    //        std::cout<<"Warning: it removes a starting pah before adding a new one"<<std::endl;
    //    }
    //    m_numofInceptedPAH++;
    //}
    //SetNumOfInceptedPAH(-1, m_particles[i]->Primary());
    //SetNumOfInceptedPAH(1);
    // Check index is within range.
    if (i<m_count) {
		//csl37: for only one particle, set to NULL
		//expand to multiple
	//	if(m_tracked_particles[0] == m_particles[i]) m_tracked_particles[0] = NULL; //csl37- no tracking
//		if(m_tracked_particles[0] == m_particles[i]) m_tracked_particles[0] = NULL;		//csl37-ttip

        // First delete current particle, then
        // set pointer to new particle.
        delete m_particles[i];
        m_particles[i] = &sp;

        m_tree.replace(m_tree.begin() + i, tree_type::value_type(sp, m_particles.begin() + i));
    }
    assert(m_tree.size() == m_count);
}

/*!
 *  Clears all particles from the ensemble, deleting them to release the
 *  memory and reseting all details of the ensemble except its capacity
 */
void Sweep::Ensemble::Clear() {
    ClearMain();
}


/*!
 *  Clears all main particles from the ensemble, deleting them to release the
 *  memory
 */
void Sweep::Ensemble::ClearMain()
{
    // Delete particles from memory and delete vectors.
    for (PartPtrVector::size_type i = 0; i != m_particles.size(); ++i) {
        delete m_particles[i];
        m_particles[i] = NULL;
    }
    m_count = 0;
    //m_numofInceptedPAH = 0;

    m_ncont = 0; // No contractions any more.

    m_tree.clear();

    // Reset doubling.
    m_maxcount   = 0;
    m_ndble      = 0;
    m_dbleactive = false;

	//csl37: set tracked particle pointers to NULL
    for (PartPtrVector::size_type i = 0; i != m_tracked_particles.size(); ++i) {
        m_tracked_particles[i] = NULL;
    }

}

// SELECTING PARTICLES.

/*!
 * @param[in,out]   rng    Random number generator
 *
 * @return      Index of a uniformly selected particle from the ensemble
 */
int Sweep::Ensemble::Select(rng_type &rng) const
{
    assert(m_tree.size() == m_count);

    // Set up the rng sample
    boost::uniform_smallint<int> indexDistrib(0, m_count - 1);
    boost::variate_generator<Sweep::rng_type&, boost::uniform_smallint<int> > indexGenerator(rng, indexDistrib);

    // Take one sample and return it
    return indexGenerator();
}


/*!
 * @param[in]       id     Property by which to weight particle selection
 * @param[in,out]   rng    Random number generator
 *
 * @return      Index of selected particle
 *
 * id must refer to a basic property from the ParticleData class
 */
int Sweep::Ensemble::Select(Sweep::PropID id, rng_type &rng) const
{
    assert(m_tree.size() == m_count);

    // This routine uses the binary tree to select a particle weighted
    // by a given particle property (by index).

    // Do not try to use the tree for uniform selection
    if(id == Sweep::iUniform)
        return Select(rng);

    // Calculate random number weighted by sum of desired property (wtid).
    // Set up the rng sample
    boost::uniform_01<rng_type&, double> unifDistrib(rng);
    double r = unifDistrib() * m_tree.head().Property(id);

    WeightExtractor we(id);
    assert(abs((m_tree.head().Property(id) - we(m_tree.head()))/m_tree.head().Property(id)) < 1e-9);
    tree_type::const_iterator it2 = m_tree.select(r, we);

    return (it2 - m_tree.begin());
}

// SCALING AND PARTICLE DOUBLING.

// Returns the scaling factor due to internal ensemble processes.
double Sweep::Ensemble::Scaling() const
{
    // The scaling factor includes the contraction term and the doubling term.
    return pow(m_contfactor, (double)m_ncont) * pow(2.0,(double)m_ndble);
}

// Resets the ensemble scaling.
void Sweep::Ensemble::ResetScaling()
{
    m_ncont = 0;
    m_ndble = 0;
    m_contwarn = false;
}


// GET SUMS OF PROPERTIES.

// Returns a ParticleData object which contains property
// sums for all particles in the ensemble.
const Sweep::Ensemble::particle_cache_type & Sweep::Ensemble::GetSums(void) const
{
    return m_tree.head();
}

// Returns the sum of a property in the ParticleData class
// over all particles.
double Sweep::Ensemble::GetSum(Sweep::PropID id) const
{
    if(id != Sweep::iUniform)
        return m_tree.head().Property(id);
    else
        return m_count;
}

/*!
 * Returns the fitting factor 'alpha' for the particle ensemble. Reference is
 * Appel, J. et al (2000) Combustion & Flame 121, 122-136, also known as the
 * ABF soot model.
 *
 *  alpha = tanh(a/log10{mu1} + b)
 *
 *  a = 12.65 - 0.00563 * T
 *  b = -1.38 + 0.00068 * T
 *  mu1 = average number of C atoms per particle
 *
 * Note that mu1 is called the 'reduced first size moment' and is usually the
 * average mass of particles in the ensemble; however in the referenced
 * version the mass is represented by the number of carbon atoms per particle.
 * See Frenklach & Wang (1994) in Soot Formation in Combustion: Mechanisms
 * and Models pp 165.
 *
 * Here the average mass is first estimated using the Tree Cache, then the
 * number of C atoms using the the molecular weight of carbon (hard-coded).
 * Could automatically get MW through something like:
 * m_particles[0]->Primary()->ParticleModel()->Components()[0]->MolWt()
 *
 * @param[in]   T   Local temperature
 *
 * @return          Alpha for the ensemble (ABF model)
 */
double Sweep::Ensemble::Alpha(double T) const {
    double alpha(0.0), mu1(0.0);

    // Get mu1 (average mass per particle)
    // Use the cache for MUCH faster calculation.
    double mW = Count()>0 ? GetSum(iWM) : 0.0;
    double invTotalWeight = Count()>0 ? 1.0/GetSum(iW) : 0.0;
    mu1 = mW * invTotalWeight * Sweep::NA / 0.01201;    // Mw in kg/mol

    // Now find alpha
    if (mu1 > 0.0) {
        double a = 12.65 - 0.00563 * T;
        double b = -1.38 + 0.00068 * T;
        alpha = std::max(0.0, tanh((a / log10(mu1)) + b));
    }

    return alpha;
}

// UPDATE ENSEMBLE.

/*!
 * @param[in]   i       Index of particle which has been updated
 *
 * This method tells the ensemble that the particle at index i
 * has been changed by the calling code, so that the ensemble
 * can update its internal data structures to reflect the new
 * particle properties.  This allows for coagulation in place,
 * rather than particle copying during coagulation.
 */
void Sweep::Ensemble::Update(unsigned int i)
{
    m_tree.replace(m_tree.begin() + i, tree_type::value_type(*m_particles[i], m_particles.begin() + i));
}

/*!
 * Replace the contents of the weights tree
 */
void Ensemble::rebuildTree() {

    // Iterators to loop over all the particles
    iterator itPart = begin();
    const iterator itPartEnd = end();

    // Build up new data for binary tree
    std::vector<std::pair<tree_type::weight_type, tree_type::return_pointer_type> > newTreeValues;
    newTreeValues.reserve(m_count);
    while(itPart != itPartEnd) {
        newTreeValues.push_back(std::make_pair(static_cast<tree_type::weight_type>(**itPart), itPart));
        ++itPart;
    }

    // Put the data into the tree
    m_tree.assign(newTreeValues.begin(), newTreeValues.end());
}

// PRIVATE FUNCTIONS.

/*!
 * The doubling algorithm is activated if the number of particles
 * in the ensemble falls below half capacity.  It copies the whole particle
 * list and changes the scaling factor to keep it consistent.  Once the
 * ensemble is back above half full, the routine updates the binary tree.
 *
 *@exception    std::runtime_error  Attempt to double with 0 particles
 */
void Sweep::Ensemble::dble()
{
    // The doubling algorithm is activated if the number of particles
    // in the ensemble falls below half capacity.  It copies the whole particle
    // list and changes the scaling factor to keep it consistent.  Once the
    // ensemble is back above half full, the routine updates the binary tree.

    // Check that doubling is on and the activation condition has been met.
    if (m_dbleon && m_dbleactive && m_count > 0) {
        const unsigned originalCount = m_count;
		bool proceed = true;

        // Continue while there are too few particles in the ensemble.
        while (m_count < m_dblelimit && proceed) {
            if(m_count == 0) {
                throw std::runtime_error("Attempt to double particle ensemble with 0 particles");
            }
			std::cout << "Doubling!" <<std::endl;
			std::cout << m_count << std::endl;

            // Copy particles.
            const size_t prevCount = m_count;
			int ii = 0;
			bool IWDSA;
			if (m_particles[0]->Primary()->AggID() == AggModels::PAH_KMC_ID){
				IWDSA = m_particles[0]->Primary()->ParticleModel()->Components(0)->WeightedPAHs();
			}
			else{
				IWDSA = false;
			}
            for (size_t i = 0; i != prevCount; ++i) {

            //if (m_particles[i]->Primary()->AggID() ==AggModels::PAH_KMC_ID)
            //{
            //    const Sweep::AggModels::PAHPrimary *rhsparticle = NULL;
            //    rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(m_particles[i]->Primary());
            //    // if not 0, it is a pyrene
            //    if (rhsparticle->Pyrene()!=0)
            //        m_numofInceptedPAH++;
            //}
				int numberPAH = 0;
				if (IWDSA){
					const Sweep::AggModels::PAHPrimary *rhsparticle = NULL;
					if (m_particles[i]->Primary()->AggID() == AggModels::PAH_KMC_ID){

						rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(m_particles[i]->Primary());
						numberPAH = rhsparticle->NumPAH();
					}
				}
				if (numberPAH > 1 || !IWDSA){ //If this particle is not just a single PAH
					
					size_t iCopy = prevCount + ii;
					// Create a copy of a particle and add it to the ensemble.
					m_particles[iCopy] = m_particles[i]->Clone();

					// Keep count of the added particles
					++m_count;
					++ii;
				}
				else{ //If this particle is a single PAH, double its statistical weight
					double oldweight = m_particles[i]->getStatisticalWeight();
					m_particles[i]->setStatisticalWeight(2.0*oldweight);
					Update(i);
				}
            }

            // Update scaling.
            ++m_ndble;

			std::cout << "Doubling done" << std::endl;
			std::cout << m_count << std::endl;
			if (IWDSA) proceed = false;
        }

        m_maxcount = std::max(m_maxcount, m_count);

        // Reset the contents of the binary tree to match the new population, if it has been changed
        if(originalCount < m_count)
            rebuildTree();

		//csl37-tracking
		//unflag primaries in untracked particles that have been copies from tracked particles 
		/*	//csl37-no tracking
		/*
		for(int j = 0; j != m_count; j++) {
			//if particle is not tracked unflag primaries
			//currently assumes that only one particle is tracked
			if(m_particles[j] != m_tracked_particles[0]) {
				m_particles[j]->removeTracking();
			}
		}
		*/
    }
}

/*!
 * Empty the tree and pass of list of pointers to the particles in
 * the tree to the caller, which must take ownership of them.  This
 * clears all scaling information in the tree, but leaves all the
 * storage allocated.
 */
Sweep::PartPtrList Sweep::Ensemble::TakeParticles() {
    // Copy the pointers to particles
    PartPtrList listOfParticles(begin(), end());

    // Now set the pointers in m_particles to NULL so that cannot be used
    // to delete the particles that will now be owned by the caller
    iterator it = begin();
    const iterator itEnd = end();
    while(it != itEnd) {
        *it++ = NULL;
    }

    // Reset the tree and other data
    ClearMain();

    return listOfParticles;
}

/*
 * @brief Writes the object to a binary stream.
 *
 * @param        out                 Output binary stream
 *
 * @exception    invalid_argument    Stream not ready
 */
void Sweep::Ensemble::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output ensemble capacity.
        unsigned int n = (unsigned int)m_capacity;
        out.write((char*)&n, sizeof(n));

        // Nothing more to do if ensemble has 0 capacity
        if(n == 0)
            return;

        // Output ensemble particle count.
        n = (unsigned int)m_count;
        out.write((char*)&n, sizeof(n));

        // Output the particles.
		std::set<void*> uniquePAHAdresses;
        for (unsigned int i=0; i!=m_count; ++i) {
            m_particles[i]->Serialize(out, &uniquePAHAdresses);
        }

        // Output number of contractions.
        n = (unsigned int)m_ncont;
        out.write((char*)&n, sizeof(n));

        // Output number of doublings.
        n = (unsigned int)m_ndble;
        out.write((char*)&n, sizeof(n));

        // Output if doubling is active.
        if (m_dbleactive) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Output if doubling is turned on.
        if (m_dbleon) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Contraction warning flag.
        if (m_contwarn) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

		//csl37
		// Output the tracked particles.
		n = (unsigned int)m_tracked_number;
        out.write((char*)&n, sizeof(n));
        for (unsigned int i=0; i!=m_tracked_number; ++i) {
            m_tracked_particles[i]->Serialize(out, &uniquePAHAdresses);
        }
		

    } else {
        throw std::invalid_argument("Output stream not ready "
                               "(Sweep, Ensemble::Serialize).");
    }
}

/*
 * @brief Reads the object from the binary stream.
 *
 * @param[in,out]    in                  Input binary stream
 * @param[in]        model	             Particle model defining interpretation of particle data
 *
 * @exception        invalid_argument    Stream not ready
 * @exception        runtime_error       Invalid serialized version number
 */
void Sweep::Ensemble::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    Clear();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;

        switch (version) {
            case 0:
			{
                // Read the ensemble capacity.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // capacity of 0 should mean the ensemble was never initialised
                if (n == 0){
                    init();
                    return;
                }
                Initialise(n);

                // Read the particle count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_count = n;

                // Read the particles.
                // Provide a way to detect multiple instances of PAHs
                std::map<void*, boost::shared_ptr<AggModels::PAHPrimary> > duplicates;
                for (unsigned int i=0; i!=m_count; ++i) {
                    Particle *p = new Particle(in, model, &duplicates);
                    m_particles[i] = p;
                }

                // Read number of contractions.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_ncont = n;

                // Read number of doublings.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_ndble = n;

                // Read if doubling is active.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n==1) {
                    m_dbleactive = true;
                } else {
                    m_dbleactive = false;
                }

                // Read if doubling is turned on.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n==1) {
                    m_dbleon = true;
                } else {
                    m_dbleon = false;
                }

                // Read contraction warning flag.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n==1) {
                    m_contwarn = true;
                } else {
                    m_contwarn = false;
                }

				//csl37
				// Read the particles.
				in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_tracked_number = n;
                for (unsigned int i=0; i!=m_tracked_number; ++i) {
                    Particle *p = new Particle(in, model, &duplicates);
                    m_tracked_particles[i] = p;
                }
				
			
                // Calculate binary tree.
                rebuildTree();

                break;
			}
            default:
                throw std::runtime_error("Serialized version number is invalid "
                                    "(Sweep, Ensemble::Deserialize).");
        }
    } else {
        throw std::invalid_argument("Input stream not ready "
                               "(Sweep, Ensemble::Deserialize).");
    }
}

// MEMORY MANAGEMENT.

// Releases all memory resources used by the ensemble.
void Sweep::Ensemble::releaseMem(void)
{
    // Delete particles from memory and delete vectors.
    for (int i=0; i!=(int)m_particles.size(); ++i) {
        delete m_particles[i];
        m_particles[i] = NULL;
    }
    m_particles.clear();

	//csl37
	// Delete particles from memory and delete vectors.
    for (int i=0; i!=(int)m_tracked_particles.size(); ++i) {
        delete m_tracked_particles[i];
        m_tracked_particles[i] = NULL;
    }
    m_tracked_particles.clear();
}

// Sets the ensemble to its initial condition.  Used in constructors.
void Sweep::Ensemble::init(void)
{
    releaseMem();

    // Capacity.
    m_levels     = 0;
    m_capacity   = 0;
    m_halfcap    = 0;
    m_count      = 0;
    //m_numofInceptedPAH = 0;

    // Scaling.
    m_contfactor = 0;
    m_ncont      = 0;
    m_contwarn   = false;

    // Doubling algorithm.
    m_maxcount   = 0;
    m_ndble      = 0;
    m_dbleactive = false;
    m_dblecutoff = 0;
    m_dblelimit  = 0;
    m_dbleslack  = 0;
    m_dbleon     = true;

	//csl37
	m_tracked_number = 0;

}

//int Sweep::Ensemble::NumOfInceptedPAH() const
//{
    //int numofpyrene = 0;
    //for (int i =0;i<m_count;i++){
    //    const Sweep::AggModels::PAHPrimary *rhsparticle = NULL;
    //    rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(m_particles[i]->Primary());

    //    numofpyrene += rhsparticle->InceptedPAH();
    //    }
    //if (numofpyrene!= m_numofInceptedPAH)
    //    std::cout<<"something goes wrong, the number ofPAH in ensemble is not consistent"<<std::endl;
    //return numofpyrene;
//    return m_numofInceptedPAH;
//}

/*!
 * Iterate through all the stochastic particles and count the number of incepted PAHs
 *
 * An incepted PAH is a stochastic particle made up of a single primary and a single PAH matching the gas transfer species (PAH that bridges gas-phase profile and MOPS)
 */
int Sweep::Ensemble::NumOfInceptedPAH(int ID) const
{
    int numOfInceptedPAHs = 0;

    if (ID == AggModels::Spherical_ID || ID == AggModels::BinTree_ID) {
        for (int i = 0; i < m_count; i++){
            if (m_particles[i]->Primary()->InceptedPAH()) {
                numOfInceptedPAHs += 1;
            }
        }
    } else {
        for (int i = 0; i < m_count; i++){
            const Sweep::AggModels::PAHPrimary *rhsparticle = NULL;
            rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(m_particles[i]->Primary());

            numOfInceptedPAHs += rhsparticle->InceptedPAH();
        }
    }

    return numOfInceptedPAHs;
}

int Sweep::Ensemble::IndexOfInceptedPAH(int ID) const
{
    if (ID == AggModels::Spherical_ID) {
        for (int i =m_count-1;i>=0;--i){
            if (m_particles[i]->Primary()->InceptedPAH())
                return i;
	    }
    } else {
        for (int i =m_count-1;i>=0;--i){
            const Sweep::AggModels::PAHPrimary *rhsparticle = NULL;
            rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(m_particles[i]->Primary());
            if (rhsparticle->InceptedPAH() == 1)
                return i;
        }
    }

    return -1;
}

//void Sweep::Ensemble::SetNumOfInceptedPAH(int m_amount, Sweep::AggModels::Primary *m_primary)
//{
//    if (m_primary->AggID() ==AggModels::PAH_KMC_ID)
//    {
//        //m_primary->ParticleModel()->IsPyreneInception()
//        const Sweep::AggModels::PAHPrimary *rhsparticle = NULL;
//        rhsparticle = dynamic_cast<const Sweep::AggModels::PAHPrimary*>(m_primary);
//
//        if (rhsparticle->InceptedPAH()!=0)
//            // rhsparticle is pyrene
//            SetNumOfInceptedPAH(m_amount);
//    }
//
//}
//void Sweep::Ensemble::SetNumOfInceptedPAH(int m_amount)
//{
//    m_numofInceptedPAH += m_amount;
//}
/*!
 * @param[in]   id      Index of property which will be extracted
 */
Ensemble::WeightExtractor::WeightExtractor(const Sweep::PropID id)
: mId(id)
{}

/*!
 * @param[in]   cache   Particle cache from which to extract a weight
 *
 * @return      The weight extracted from the cache
 */
double Ensemble::WeightExtractor::operator()(const particle_cache_type& cache) const {
    return cache.Property(mId);
}

void Sweep::Ensemble::UpdateTracking(int p_old, int p_merged){

	//csl37
	//works for one particle at the moment
	if(m_tracked_particles[0] == m_particles[p_old])
		m_tracked_particles[0] = m_particles[p_merged];

}

// Returns a pointer to the given tracked.
Particle *const Sweep::Ensemble::TrackedAt(unsigned int i)
{
    // Check that the index in within range, then return the particle.
    if (i < m_tracked_number) {
        return m_tracked_particles[i];
    } else {
        return NULL;
    }
}

void Sweep::Ensemble::InitialiseTracking(){
	//csl37: add tracked particles
	m_tracked_number = 0;	//only done for one particle	//csl37-no tracking
	m_tracked_number = 0;	//only done for one particle	//csl37-TTIP
	int i = 0;
	while(i < m_tracked_number && i<m_count){
		m_tracked_particles[i] = m_particles[i];
		//initialise tracking of primary
		m_tracked_particles[i]->setTracking();
		i++;
	}
}