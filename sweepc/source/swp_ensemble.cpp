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
#include "swp_submodel_type.h"
#include "swp_submodel.h"
#include "swp_particle_cache.h"

#include <cmath>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <limits>
#include <functional>
#include <algorithm>

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Sweep::Ensemble::Ensemble(void)
{
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
            // Scaling.
            m_ncont      = rhs.m_ncont;
            m_scale      = rhs.m_scale;
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
        }
        m_tree = rhs.m_tree;
    }
    return *this;
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
    real rl    = log((real)capacity) / log(2.0);
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

    // Reserve memory for ensemble.
    m_particles.resize(m_capacity, NULL);

    m_tree.resize(m_capacity);

    // Initialise scaling.
    m_ncont      = 0;
    m_scale      = 1.0;
    m_contfactor = (real)(m_capacity-1) / (real)(m_capacity);
    m_contwarn   = false;

    // Initialise doubling.
    m_maxcount   = 0;
    m_ndble      = 0;
    m_dbleon     = true;
    m_dbleactive = false;
    m_dblecutoff = (int)(3.0 * (real)m_capacity / 4.0);
    m_dblelimit  = (m_halfcap - (unsigned int)pow(2.0, (int)((m_levels-5)>0 ? m_levels-5 : 0)));
    m_dbleslack  = (unsigned int)pow(2.0, (int)((m_levels-5)>0 ? m_levels-5 : 0));
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

// Adds the particle to the ensemble.  Takes control of particle
// destruction.
int Sweep::Ensemble::Add(Particle &sp, int (*rand_int)(int, int))
{
    // Check for doubling activation.
    if (!m_dbleactive && (m_count >= m_dblecutoff-1)) {
        m_dbleactive = true;
        //printf("sweep: Particle doubling activated.\n");
    }

    // Check ensemble for space, if there is not enough space then need
    // to generate some by contracting the ensemble.
    int i = -1;
    if (m_count < m_capacity) {
        // There is space in the tree for a new particle.
        i = -1;
    } else {
        // We must contract the ensemble to accomodate a new particle.
        i = rand_int(0, m_capacity);
		++m_ncont;
        if (!m_contwarn && ((real)(m_ncont)/(real)m_capacity > 0.01)) {
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
    } else if ((unsigned)i < m_capacity) {
        // Replace an existing particle (if i=m_capacity) then
        // we are removing the new particle, so just ignore it.
        Replace(i, sp);
    } else {
        // The new particle is to be removed immediately
        assert(static_cast<unsigned int>(i) == m_capacity);
        delete &sp;
    }

    m_maxcount = std::max(m_maxcount, m_count);

    assert(m_tree.size() == m_count);

    return i;
}

// Removes the particle at the given index, if the index is
// valid.
void Sweep::Ensemble::Remove(unsigned int i, bool fdel)
{
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
    dble();
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
    m_count = std::distance(m_particles.begin(), validEnd);

    // Now delete the invalid particles and nullify the corresponding pointers
    while(validEnd != m_particles.end()) {
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

// Replaces the particle at the given index with the given particle,
// as long as the index is valid.
void Sweep::Ensemble::Replace(unsigned int i, Particle &sp)
{
    // Check index is within range.
    if (i<m_count) {
        // First delete current particle, then
        // set pointer to new particle.
        delete m_particles[i];
        m_particles[i] = &sp;

        m_tree.replace(m_tree.begin() + i, tree_type::value_type(sp, m_particles.begin() + i));
    }
    assert(m_tree.size() == m_count);
}

// Clears all particles from the ensemble.
void Sweep::Ensemble::Clear()
{
    // Delete particles from memory and delete vectors.
    for (PartPtrVector::size_type i = 0; i != m_particles.size(); ++i) {
        delete m_particles[i];
        m_particles[i] = NULL;
    }
    m_count = 0;

    m_ncont = 0; // No contractions any more.

    m_tree.clear();

    // Reset doubling.
    m_maxcount   = 0;
    m_ndble      = 0;
    m_dbleactive = false;
}


// SELECTING PARTICLES.

/*!
 * @param[in]   rand_int    Pointer to function that can generate uniform integers on a range
 *
 * @return      Index of a uniformly selected particle from the ensemble
 */
int Sweep::Ensemble::Select(int (*rand_int)(int, int)) const
{
    assert(m_tree.size() == m_count);

    // Uniformly select a particle index.
    return rand_int(0, m_count-1);
}

/*!
 * @param[in]   id          Property by which to weight particle selection
 * @param[in]   rand_int    Pointer to function that generates uniform integers on a range
 * @param[in]   rand_u01    Pointer to function that generates U[0,1] deviates
 *
 * @return      Index of selected particle
 *
 * id must refer to a basic property from the ParticleData class
 */
int Sweep::Ensemble::Select(ParticleCache::PropID id,
                            int (*rand_int)(int, int),
                            real (*rand_u01)()) const
{
    assert(m_tree.size() == m_count);

    // This routine uses the binary tree to select a particle weighted
    // by a given particle property (by index).

    // Do not try to use the tree for uniform selection
    if(id == ParticleCache::iUniform)
        return Select(rand_int);

    // Calculate random number weighted by sum of desired property (wtid).
    real r = rand_u01() * m_tree.head().Property(id);

    WeightExtractor we(id);
    assert(abs((m_tree.head().Property(id) - we(m_tree.head()))/m_tree.head().Property(id)) < 1e-9);
    tree_type::const_iterator it2 = m_tree.select(r, we);

    return (it2 - m_tree.begin());
}

// SCALING AND PARTICLE DOUBLING.

// Returns the scaling factor due to internal ensemble processes.
real Sweep::Ensemble::Scaling() const
{
    // The scaling factor includes the contraction term and the doubling term.
    return m_scale * pow(m_contfactor, (double)m_ncont) * pow(2.0,(double)m_ndble);
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
const ParticleCache & Sweep::Ensemble::GetSums(void) const
{
    return m_tree.head();
}

// Returns the sum of a property in the ParticleData class
// over all particles.
real Sweep::Ensemble::GetSum(ParticleCache::PropID id) const
{
    if(id != ParticleCache::iUniform)
        return m_tree.head().Property(id);
    else
        return m_count;
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
    std::vector<tree_type::value_type> newTreeValues(m_count);
    std::vector<tree_type::value_type>::iterator newTreeValuesIt = newTreeValues.begin();
    while(itPart != itPartEnd) {
        *newTreeValuesIt++ = tree_type::value_type(**itPart, itPart);
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
    if (m_dbleon && m_dbleactive) {
       const unsigned originalCount = m_count;

        // Continue while there are too few particles in the ensemble.
        while (m_count < m_dblelimit) {
            if(m_count == 0) {
                throw std::runtime_error("Attempt to double particle ensemble with 0 particles");
            }

            // Copy particles.
            const size_t prevCount = m_count;
            for (size_t i = 0; i != prevCount; ++i) {

                size_t iCopy = prevCount + i;
                // Create a copy of a particle and add it to the ensemble.
                m_particles[iCopy] = m_particles[i]->Clone();

                // Keep count of the added particles
                ++m_count;
            }

            // Update scaling.
            ++m_ndble;
        }

        m_maxcount = std::max(m_maxcount, m_count);

        // Reset the contents of the binary tree to match the new population, if it has been changed
        if(originalCount < m_count);
            rebuildTree();
    }
}
/*!
 * Empty the tree and pass of list of pointers to the particles in
 * the tree to the caller, which must take ownership of them.
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
    Clear();

    return listOfParticles;
}

// READ/WRITE/COPY.

// Writes the object to a binary stream.
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
        for (unsigned int i=0; i!=m_count; ++i) {
            m_particles[i]->Serialize(out);
        }

        // Output the scaling factor.
        double val = (double)m_scale;
        out.write((char*)&val, sizeof(val));

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
    } else {
        throw std::invalid_argument("Output stream not ready "
                               "(Sweep, Ensemble::Serialize).");
    }
}

// Reads the object from a binary stream.
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
        double val     = 0.0;

        switch (version) {
            case 0:
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
                for (unsigned int i=0; i!=m_count; ++i) {
                    Particle *p = new Particle(in, model);
                    m_particles[i] = p;
                }

                // Read the scaling factor.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_scale = (real)val;

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

                // Calculate binary tree.
                rebuildTree();

                break;
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

    // Scaling.
    m_scale      = 1.0;
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
}

/*!
 * @param[in]   id      Index of property which will be extracted
 */
Ensemble::WeightExtractor::WeightExtractor(const ParticleCache::PropID id)
: mId(id)
{}

/*!
 * @param[in]   cache   Particle cache from which to extract a weight
 *
 * @return      The weight extracted from the cache
 */
real Ensemble::WeightExtractor::operator()(const ParticleCache& cache) const {
    return cache.Property(mId);
}
