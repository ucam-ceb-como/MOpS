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
#include "rng.h"
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

// Initialising constructor (no particles).
Sweep::Ensemble::Ensemble(const Sweep::ParticleModel &model)
: m_sums(0.0, model)
{
    init();
    m_model = &model;
}

// Initialising constructor.
Sweep::Ensemble::Ensemble(unsigned int count, const Sweep::ParticleModel &model)
{
    // Call initialisation routine.
    Initialise(count, model);
}

// Copy contructor.
Sweep::Ensemble::Ensemble(const Sweep::Ensemble &copy)
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

            // Resize binary tree vector.
            m_model = rhs.m_model;
            m_tree.resize(rhs.m_capacity-1, TreeNode(*m_model));

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

            // Set up the pointers to the related nodes
            unsigned int i, j;
            m_tree[0].Parent = NULL;
            for (i=0,j=1; i!=m_halfcap-1; ++i) {
                j = (2*i) + 1;
                m_tree[i].Left          = &m_tree.at(j);
                m_tree[i].Left->Parent  = &(m_tree[i]);
                m_tree[i].Right         = &m_tree.at(j+1);
                m_tree[i].Right->Parent = &(m_tree[i]);
            }
            //  Nodes in last half of tree refer to particles, and
            //  have no children.
            for (i=m_halfcap-1; i!=m_capacity-1; ++i) {
                m_tree[i].Left  = NULL;
                m_tree[i].Right = NULL;
            }

            // Now copy data from rhs tree to this (pointers are not affected)
            for (unsigned int i=0; i!=m_capacity-1; ++i) {
                m_tree[i] = rhs.m_tree[i];
            }

            m_sums = rhs.m_sums;
        }
    }
    return *this;
}


// INITIALISATION.

// Initialises the ensemble to the given size.  Any particles currently
// in the ensemble will be destroyed.
void Sweep::Ensemble::Initialise(unsigned int capacity, const Sweep::ParticleModel &model)
{
    // Clear current ensemble.
    Clear();
 
    // Store new particle model.
    m_model = &model;

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
    m_tree.resize(m_capacity-1, TreeNode(model));
    
    // Set all tree nodes to correct size and link up tree.
    unsigned int i, j;
    m_tree[0].Parent = NULL;
    for (i=0,j=1; i<m_halfcap-1; i++,j=(2*i)+1) {
        m_tree[i].Left = &m_tree.at(j);
        m_tree[i].Left->Parent = &(m_tree[i]);
        m_tree[i].Right = &m_tree.at(j+1);
        m_tree[i].Right->Parent = &(m_tree[i]);
    }
    for (i=m_halfcap-1; i!=m_capacity-1; ++i) {
        m_tree[i].Left  = NULL;
        m_tree[i].Right = NULL;
    }

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


// THE PARTICLE MODEL.

// Returns a pointer to the particle model to which this
// ensemble subscribes.
const ParticleModel *const Sweep::Ensemble::ParticleModel(void) const
{
    return m_model;
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
int Sweep::Ensemble::Add(Particle &sp)
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
        i = irnd(0, m_capacity);
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
    } else if ((unsigned)i < m_capacity) {
        // Replace an existing particle (if i=m_capacity) then
        // we are removing the new particle, so just ignore it.
        Replace(i, sp);
    } else {
        // The new particle is to be removed immediately
        assert(static_cast<unsigned int>(i) == m_capacity);
        delete &sp;
    }

    // Now we must recalculate the tree by inserting the particle properties
    // into the bottom row and calculating up.
    if ((unsigned)i < m_capacity) {
        int j = treeIndex((unsigned)i);
        if (isLeftBranch((unsigned)i)) {
            m_tree[j].LeftData = *m_particles[i];
        } else {
            m_tree[j].RightData = *m_particles[i];
        }
        ascendingRecalc(j);
    }

    m_maxcount = max(m_maxcount, m_count);

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

        // Recalculate the tree up from leaf.
        int k = treeIndex(i);
        if (isLeftBranch(i)) {
            m_tree[k].LeftData = *m_particles[i];
        } else {
            m_tree[k].RightData = *m_particles[i];
        }
        ascendingRecalc(k);

        // Recalculate tree up from last position (clear it first).
        k = treeIndex(m_count);
        if (isLeftBranch(m_count)) {
            m_tree[k].LeftData.Clear();
        } else {
            m_tree[k].RightData.Clear();
        }
        ascendingRecalc(k);
    } else if (i==m_count-1) {
        // This is the last particle in the ensemble, we don't
        // need to swap it with another, just delete it.

        // Erase particle.
        if (fdel) delete m_particles[i];
        m_particles[i] = NULL;
        --m_count;

        // Recalculate tree up from last position (clear it first).
        int k = treeIndex(i);
        if (isLeftBranch(i)) {
            m_tree[k].LeftData.Clear();
        } else {
            m_tree[k].RightData.Clear();
        }
        ascendingRecalc(k);
    }

    // Particle removal might reduce the particle count 
    // sufficiently to require particle doubling.
    dble();
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
    Update();

    // If we removed too many invalid particles then we'll have to double.
    dble();
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

        // Recalculate tree up from branch.
        int j = treeIndex(i);
        if (isLeftBranch(i)) {
            m_tree[j].LeftData = sp;
        } else {
            m_tree[j].RightData = sp;
        }
        ascendingRecalc(j);
    }
}

// Clears all particles from the ensemble.
void Sweep::Ensemble::Clear()
{
    // Delete particles from memory and delete vectors.
    for (int i=0; i!=(int)m_particles.size(); ++i) {
        delete m_particles[i];
        m_particles[i] = NULL;
    }
    m_count = 0;

    m_ncont = 0; // No contractions any more.

    vector<TreeNode>::iterator j;
    for(j=m_tree.begin(); j!=m_tree.end(); ++j) j->Clear();

    m_sums.Clear();

    // Reset doubling.
    m_maxcount   = 0;
    m_ndble      = 0;
    m_dbleactive = false;
}


// SELECTING PARTICLES.

// Returns the index of a uniformly selected particle
// from the ensemble.
int Sweep::Ensemble::Select(void) const
{
    // Uniformly select a particle index.
    return irnd(0, m_count-1);
}

// Returns the index of a particle selected using the property
// weight given which refers to a basic property in the
// ParticleData class.
int Sweep::Ensemble::Select(ParticleCache::PropID id) const
{
    // This routine uses the binary tree to select a particle weighted
    // by a given particle property (by index).

    // Do not try to use the tree for uniform selection
    if(id == ParticleCache::iUniform)
        return Select();

    int isp=-1;

    // Calculate random number weighted by sum of desired property (wtid).
    real r = rnd() * (m_tree[0].LeftData.Property(id) + 
                      m_tree[0].RightData.Property(id));

    // Fall down the binary tree until reaching a base node.
    const TreeNode *n = &m_tree[0];
    int j=0;
    while(n->Left!=NULL) {
        if (r <= n->LeftData.Property(id)) {
            n = n->Left;
            j = (2*j) + 1;
        } else {
            r -= n->LeftData.Property(id);
            n = n->Right;
            j = (2*j) + 2;
        }
    }
    // last level!
    if (r <=n->LeftData.Property(id)) {
        isp = (2*j) + 2 - m_capacity;
    } else {
        isp = (2*j) + 3 - m_capacity;
    }

    // One final check that the index we've chosen is valid.
    if ((unsigned)isp < m_count) {
        return isp;
    } else {
        // Chosen an non-existent particle.
        return -1;
    }
}

// Selects particle according to the particle property
// specified by the given model and the given property id
// of the model.
int Sweep::Ensemble::Select(SubModels::SubModelType model_id, unsigned int id) const
{
    // This routine uses the binary tree to select a particle weighted
    // by a given particle property (by index).

    if (model_id == SubModels::BasicModel_ID) return Select((ParticleCache::PropID)id);

    int isp=-1;

    // Calculate random number weighted by sum of desired property (wtid).
    real r = rnd();

    r *= (m_tree[0].LeftData.SubModel(model_id)->Property(id) + 
          m_tree[0].RightData.SubModel(model_id)->Property(id));

    // Fall down the binary tree until reaching a base node.
    const TreeNode *n = &m_tree[0];
    int j=0;
    while(n->Left!=NULL) {
        if (r <= n->LeftData.SubModel(model_id)->Property(id)) {
            n = n->Left;
            j = (2*j) + 1;
        } else {
            r -= n->LeftData.SubModel(model_id)->Property(id);
            n = n->Right;
            j = (2*j) + 2;
        }
    }
    // Last level!
    if (r <=n->LeftData.SubModel(model_id)->Property(id)) {
        isp = (2*j) + 2 - m_capacity;
    } else {
        isp = (2*j) + 3 - m_capacity;
    }

    // One final check that the index we've chosen is valid.
    if ((unsigned)isp < m_count) {
        return isp;
    } else {
        // Chosen an non-existent particle.
        return -1;
    }
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
    return m_sums;
}

// Returns the sum of a property in the ParticleData class
// over all particles.
real Sweep::Ensemble::GetSum(ParticleCache::PropID id) const
{
    if(id != ParticleCache::iUniform)
        return m_sums.Property(id);
    else
        return m_count;
}

// Returns the sum of one particle property with the given index 
// from the given model from the binary tree.
real Sweep::Ensemble::GetSum(SubModels::SubModelType model_id, unsigned int id) const
{
    if (model_id == SubModels::BasicModel_ID) {
        return GetSum(static_cast<ParticleCache::PropID>(id));
    } else {
        return m_sums.SubModel(model_id)->Property(id);
    }
}


// UPDATE ENSEMBLE.

// Updates the ensemble tree from the given particle index.
void Sweep::Ensemble::Update(unsigned int i)
{
    if (i < m_count) {
        // Get tree index of this particle.
        unsigned int j = treeIndex(i);

        // Update binary tree at this index.
        if (isLeftBranch(i)) {
            m_tree[j].LeftData = *m_particles[i];
        } else {
            m_tree[j].RightData = *m_particles[i];
        }
        ascendingRecalc(j);
    }
}

// Updates the ensemble tree completely.
void Sweep::Ensemble::Update()
{
    // Put the particle values into the bottom of the tree
    bool odd=true;
    unsigned int j=0;
    for (unsigned int i=0; i!=m_count; ++i) {
        j = treeIndex(i);
        if (odd) {
            m_tree[j].LeftData = *m_particles[i];
        } else {
            m_tree[j].RightData = *m_particles[i];
        }
        odd = !odd;
    }
    for (unsigned int i=m_count; i!=m_capacity; ++i) {
        j = treeIndex(i);
        if (odd) {
            m_tree[j].LeftData.Clear();
        } else {
            m_tree[j].RightData.Clear();
        }
        odd = !odd;
    }

    // Calculate the nodes further up the tree
    recalcAllNonLeaf();

    assert(m_tree.front().LeftData.Property(ParticleCache::iM) < std::numeric_limits<real>::max());
    assert(m_tree.front().LeftData.Property(ParticleCache::iM) >= 0.0);
}


// PRIVATE FUNCTIONS.

// Recalculates a branch of the tree from the given node upwards.
void Sweep::Ensemble::ascendingRecalc(unsigned int i)
{
    // This function starts at the bottom of the binary tree and works
    // upwards recalculating the sums of all particle properties under
    // the nodes.

    if (i<m_capacity) {
        // Get the node at the bottom of the branch.
        TreeNode *n = &m_tree[i];

        // Climb up the tree until the root node, 
        // summing up the properties.
        while (n->Parent != NULL) {
            assert(n->LeftData.Property(ParticleCache::iM) < std::numeric_limits<real>::max());
            assert(n->LeftData.Property(ParticleCache::iM) >= 0.0);
            assert(n->RightData.Property(ParticleCache::iM) < std::numeric_limits<real>::max());
            assert(n->RightData.Property(ParticleCache::iM) >= 0.0);

            n = n->Parent;
            (n->LeftData  = n->Left->LeftData)  += n->Left->RightData;
            (n->RightData = n->Right->LeftData) += n->Right->RightData;
        }
        m_sums = m_tree.front().LeftData + m_tree.front().RightData;

        assert(m_tree.front().LeftData.Property(ParticleCache::iM) < std::numeric_limits<real>::max());
        assert(m_tree.front().LeftData.Property(ParticleCache::iM) >= 0.0);
    }
}

void Sweep::Ensemble::dble()
{
    // The doubling algorithm is activated if the number of particles
    // in the ensemble falls below half capacity.  It copies the whole particle
    // list and changes the scaling factor to keep it consistent.  Once the
    // ensemble is back above half full, the routine updates the binary tree.

    // As an additional check, in case the maximum count did not reach 
    // the required threshold (75%), we check for a 20% reduction from
    // the maximum achieved particle count, assuming that this is statistically
    // significant.
    if (!m_dbleactive && (m_maxcount >= (unsigned int)(0.1*(real)m_capacity)) && 
        (m_count <= (unsigned int)(0.9*(real)m_maxcount))) {
        m_dbleactive = true;
        printf("sweep: Particle doubling activated.\n");
    }

    // Check that doubling is on and the activation condition has been met.
    if (m_dbleon && m_dbleactive) {     
        bool left;

        // Continue while there are too few particles in the ensemble.
        while (m_count < m_dblelimit) {
            printf("sweep: Doubling particle ensemble.\n");

            left = isLeftBranch(m_count);

            // Copy particles.
            const size_t prevCount = m_count;
            for (size_t i = 0; i != prevCount; ++i) {
                
                size_t iCopy = prevCount + i;
                // Create a copy of a particle and add it to the ensemble.
                m_particles[iCopy] = m_particles[i]->Clone();

                // Put properties into bottom of tree.
                size_t leafIndex = treeIndex(iCopy);
                if (left) {
                    m_tree[leafIndex].LeftData = *m_particles[iCopy];
                } else {
                    m_tree[leafIndex].RightData = *m_particles[iCopy];
                }
                
                // Keep count of the added particles
                ++m_count;
                // Next particle will be on the other side of a node
                left = !left;
            }
            
            // Update the tree from second last level up.
            recalcAllNonLeaf();

            // Update scaling.
            ++m_ndble;
        }

        m_maxcount = max(m_maxcount, m_count);
    }
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
        throw invalid_argument("Output stream not ready "
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
                Initialise(n, model);

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
                Update();

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Ensemble::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
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
    m_tree.clear();
}

// Sets the ensemble to its initial condition.  Used in constructors.
void Sweep::Ensemble::init(void)
{
    releaseMem();

    // Set particle model.
    m_model = NULL;

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

    m_sums.Clear();
}
