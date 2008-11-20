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
#include "rng.h"
#include <cmath>
#include <vector>
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Ensemble::Ensemble(void)
{
    init();
}

// Initialising constructor (no particles).
Ensemble::Ensemble(const Sweep::ParticleModel &model)
{
    init();
    m_model = &model;
}

// Initialising constructor.
Ensemble::Ensemble(unsigned int count, const Sweep::ParticleModel &model)
{
    // Call initialisation routine.
    Initialise(count, model);
}

// Copy contructor.
Ensemble::Ensemble(const Sweep::Ensemble &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
Ensemble::Ensemble(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Destructor.
Ensemble::~Ensemble(void)
{
    // Clear the ensemble.
    Clear();
}


// OPERATOR OVERLOADING.

// Assignment operator.
Ensemble &Ensemble::operator=(const Sweep::Ensemble &rhs)
{
    if (this != &rhs) {
        // Clear current particles.
        Clear();

        if (rhs.m_capacity > 0) {
            // Resize particle vector.
            m_particles.resize(rhs.m_capacity, NULL);

            // Resize binary tree vector.
            bool relink_tree = (m_capacity != rhs.m_capacity);
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

            if (relink_tree) {
                // Set all tree nodes to correct size and link up tree.
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
            }

            // Now copy data from rhs tree to this.
            for (unsigned int i=0; i!=m_capacity-1; ++i) {
                m_tree[i] = rhs.m_tree[i];
            }
        }
    }
    return *this;
}


// INITIALISATION.

// Initialises the ensemble to the given size.  Any particles currently
// in the ensemble will be destroyed.
void Ensemble::Initialise(unsigned int capacity, const Sweep::ParticleModel &model)
{
    // Clear current ensemble.
    Clear();
 
    // Store new particle model.
    m_model = &model;

    // Calculate nearest power of 2 capacity.  Ensemble capacity must be a power
    // of 2.  This constraint is due to the binary tree implementation.
    real rl    = log((real)capacity) / log(2.0);
    m_levels   = (int)(rl + 0.5);
    m_capacity = (int)pow(2.0, (int)m_levels);
    m_halfcap  = capacity / 2;
    m_count    = 0;

    // Reserve memory for ensemble.
    m_particles.resize(capacity, NULL);
    m_tree.resize(capacity-1, TreeNode(model));
    
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
    m_dblelimit  = m_halfcap - (unsigned int)pow(2.0, (int)((m_levels-5)>0 ? m_levels-5 : 0));
    m_dbleslack  = (unsigned int)pow(2.0, (int)((m_levels-5)>0 ? m_levels-5 : 0));
}

    
// THE PARTICLE MODEL.

// Returns a pointer to the particle model to which this
// ensemble subscribes.
const Sweep::ParticleModel *const Ensemble::ParticleModel(void) const
{
    return m_model;
}


// PARTICLE ADDITION AND REMOVAL.

// Returns a pointer to the particle with the given index.
Particle *const Ensemble::At(unsigned int i)
{
    // Check that the index in within range, then return the particle.
    if (i < m_count) {
        return m_particles[i];
    } else {
        return NULL;
    }
}

// Returns a pointer to the particle with the given index.
const Particle *const Ensemble::At(unsigned int i) const
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
int Ensemble::Add(Particle &sp)
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
        // We must contract the ensemble to accomodate a new particle.
        i = irnd(0, m_capacity);
        if (!m_contwarn && ((real)(++m_ncont)/(real)m_capacity > 0.01)) {
            m_contwarn = true;
            printf("sweep: Ensemble contracting too often; "
                   "possible stiffness issue.\n");
        }
    }

    if (i < 0) {
        // We are adding a new particle.
        i=m_count++;
        m_particles[i] = &sp;
        sp.SetEnsemble(*this);
    } else if ((unsigned)i < m_capacity) {
        // Replace an existing particle (if i=m_capacity) then
        // we are removing the new particle, so just ignore it.
        Replace(i, sp);
        sp.SetEnsemble(*this);
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
void Ensemble::Remove(unsigned int i, bool fdel)
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
        unsigned int j = m_count - 1;
        k = treeIndex(j);
        if (isLeftBranch(j)) {
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
void Ensemble::RemoveInvalids(void)
{
    // This function loops forward through the list finding invalid
    // particles and backwards finding valid particles.  Once an invalid
    // and a valid particle are found they are swapped.  This results in
    // all the invalid particles being collected at the end of the vector,
    // from where they can be deleted.

    unsigned int i=0, k=m_count-1;
    unsigned int ix=0, tx=0, n=0;
    bool found=false;

    if (m_particles[i] != NULL) {
        while (i<k) {
            found = false;

            // Locate next invalid particle in list.
            while ((i!=k) && m_particles[i]->IsValid()) {++i; ++ix;}

            // Locate next valid particle from end of list, remember
            // to delete invalid particles at same time.
            while ((k!=i) && (!m_particles[k]->IsValid())) {
                delete m_particles[k];
                m_particles[k] = NULL;
                --m_count;
                --k;
            }

            if (i<k) {            
                delete m_particles[i]; // Delete invalid particle from memory.
                m_particles[i] = m_particles[k]; // Put valid particle at index i.
                m_particles[k] = NULL;           // Clear duplicate pointer.
                ++n;       // Increment number of invalid particles.
                --m_count; // Decrement number of particles.

                // Update binary tree.
                tx = treeIndex(ix);
                if (isLeftBranch(ix)) {
                    m_tree[tx].LeftData = *m_particles[i];
                } else {
                    m_tree[tx].RightData = *m_particles[i];
                }
                ascendingRecalc(tx);

                // Move on indices.
                ++i; --k;
            } else {
                // All invalid particles have been deleted.
                break;
            }
        }
    }

    // If we removed too many invalid particles then we'll have to double.
    dble();
}

// Replaces the particle at the given index with the given particle,
// as long as the index is valid.
void Ensemble::Replace(unsigned int i, Particle &sp)
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
void Ensemble::Clear()
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

    // Reset doubling.
    m_maxcount   = 0;
    m_ndble      = 0;
    m_dbleactive = false;
}


// SELECTING PARTICLES.

// Returns the index of a uniformly selected particle
// from the ensemble.
int Ensemble::Select(void) const
{
    // Uniformly select a particle index.
    return irnd(0, m_count-1);
}

// Returns the index of a particle selected using the property
// weight given which refers to a basic property in the
// ParticleData class.
int Ensemble::Select(ParticleCache::PropID id) const
{
    // This routine uses the binary tree to select a particle weighted
    // by a given particle property (by index).

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
int Ensemble::Select(SubModels::SubModelType model_id, unsigned int id) const
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
real Ensemble::Scaling() const
{
    // The scaling factor includes the contraction term and the doubling term.
    return m_scale * pow(m_contfactor, (double)m_ncont) * pow(2.0,(double)m_ndble);
}

// Resets the ensemble scaling.
void Ensemble::ResetScaling()
{
    m_ncont = 0;
    m_ndble = 0;
    m_contwarn = false;
}


// GET SUMS OF PROPERTIES.

// Returns a ParticleData object which contains property
// sums for all particles in the ensemble.
const ParticleCache &Ensemble::GetSums(void) const
{
	m_sums = m_tree[0].LeftData + m_tree[0].RightData;
    return m_sums;
}

// Returns the sum of a property in the ParticleData class
// over all particles.
real Ensemble::GetSum(ParticleCache::PropID id) const
{
    return m_tree[0].LeftData.Property(id) + m_tree[0].RightData.Property(id);
}

// Returns the sum of one particle property with the given index 
// from the given model from the binary tree.
real Ensemble::GetSum(SubModels::SubModelType model_id, unsigned int id) const
{
    if (model_id == SubModels::BasicModel_ID) {
        return GetSum(static_cast<ParticleCache::PropID>(id));
    } else {
        return m_tree[0].LeftData.SubModel(model_id)->Property(id) + 
               m_tree[0].RightData.SubModel(model_id)->Property(id);
    }
}


// UPDATE ENSEMBLE.

// Updates the ensemble tree from the given particle index.
void Ensemble::Update(unsigned int i)
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
void Ensemble::Update()
{
    // This flavour updates the whole binary tree.

    bool odd=true;
    unsigned int j=0;
    for (unsigned int i=0; i!=m_count; ++i) {
        j = treeIndex(i);
        if (odd) {
            m_tree[j].LeftData = *m_particles[i];
        } else {
            m_tree[j].RightData = *m_particles[i];
            ascendingRecalc(j);
        }
        odd = !odd;
    }
    for (unsigned int i=m_count; i!=m_capacity; ++i) {
        j = treeIndex(i);
        if (odd) {
            m_tree[j].LeftData.Clear();
        } else {
            m_tree[j].RightData.Clear();
            ascendingRecalc(j);
        }
        odd = !odd;
    }

    // Need to do one last recalc if there are an odd number of particles.
    if (!odd) {
        ascendingRecalc(j);
    }
}


// PRIVATE FUNCTIONS.

// Recalculates a branch of the tree from the given node upwards.
void Ensemble::ascendingRecalc(unsigned int i)
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
            n = n->Parent;
            (n->LeftData  = n->Left->LeftData)  += n->Left->RightData;
            (n->RightData = n->Right->LeftData) += n->Right->RightData;
        }
    }
}

void Ensemble::dble()
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
/*	if (m_dbleactive)
	{
		m_dbleon=false;
		cout << " Double?";
		string ans;
		cin >> ans;
			if (ans=="y")
				m_dbleon=true;
	}*/
	/*
	while (m_count<=0.9*(real)m_maxcount) {
				unsigned int i,k;
				 Particle* sp;
			    vector<TreeNode>::iterator inode;
                sp = m_particles[1]->Clone();
                m_particles[i=m_count++] = sp;
				// Put properties into bottom of tree.
                k = treeIndex(m_count-1);
                if (left) {
                    m_tree[k].LeftData = *m_particles[i];
                } else {
                    m_tree[k].RightData = *m_particles[i];
                }
            

            // Update the tree from second last level up.
            for (inode=(m_tree.begin()+m_halfcap-2); inode!=m_tree.begin(); --inode) {
                inode->LeftData = inode->Left->LeftData + inode->Left->RightData;
                inode->RightData = inode->Right->LeftData + inode->Right->RightData;
			}
	}

*/
    // Check that doubling is on and the activation condition has been met.
    if (m_dbleon && m_dbleactive) {
        Particle* sp;
        bool left;
        unsigned int isp, i, j, k;
        unsigned int n=m_count-1;
        vector<TreeNode>::iterator inode;
        fvector::iterator iterp, iterl, iterr;

        // Continue while there are too few particles in the ensemble.
        while (n<=m_dblelimit) {
            printf("sweep: Doubling particle ensemble.\n");

            left = isLeftBranch(n);

            // Copy particles.
            for (isp=0,j=n; isp!=n; ++isp,++j,left=!left) {
                // Create a copy of a particle and add it to the ensemble.
                sp = m_particles[isp]->Clone();
                m_particles[i=m_count++] = sp;

                // Put properties into bottom of tree.
                k = treeIndex(j);
                if (left) {
                    m_tree[k].LeftData = *m_particles[i];
                } else {
                    m_tree[k].RightData = *m_particles[i];
                }
            }

            // Update the tree from second last level up.
            for (inode=(m_tree.begin()+m_halfcap-2); inode!=m_tree.begin(); --inode) {
                inode->LeftData = inode->Left->LeftData + inode->Left->RightData;
                inode->RightData = inode->Right->LeftData + inode->Right->RightData;
            }

            // Update scaling.
            ++m_ndble;
            n *= 2;
        }

        m_maxcount = max(m_maxcount, m_count);
    }
}

// READ/WRITE/COPY.

// Writes the object to a binary stream.
void Ensemble::Serialize(std::ostream &out) const
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
void Ensemble::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
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
                    p->SetEnsemble(*this);
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
void Ensemble::releaseMem(void)
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
void Ensemble::init(void)
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
}
