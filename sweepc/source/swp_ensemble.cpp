#include "swp_ensemble.h"
#include "rng.h"
#include <cmath>
#include <vector>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Ensemble::Ensemble(void)
{
    // Capacity.
    m_levels     = 0;
    m_capacity   = 0;
    m_halfcap    = 0;
    m_count      = 0;
    // Scaling.
    m_scale      = 1.0;
    m_contfactor = 0;
    m_ncont      = 0;
    // Doubling algorithm.
    m_ndble      = 0;
    m_dbleactive = false;
    m_dblecutoff = 0;
    m_dblelimit  = 0;
    m_dbleslack  = 0;
    m_dbleon     = true;
}

// Initialising constructor.
Ensemble::Ensemble(unsigned int count)
{
    // Call initialisation routine.
    Initialise(count);
}

// Copy contructor.
Ensemble::Ensemble(const Sweep::Ensemble &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Destructor.
Ensemble::~Ensemble(void)
{
    // Clear the ensemble.
    Clear();
}


// OPERATOR OVERLOADING.

// Assignment operator.
Ensemble &Ensemble::operator =(const Sweep::Ensemble &rhs)
{
    if (this != &rhs) {
        // Clear current particles.
        Clear();
        Initialise(rhs.m_capacity);

        // Capacity.
        m_levels   = rhs.m_levels;
        m_capacity = rhs.m_capacity;
        m_halfcap  = rhs.m_halfcap;
        m_count    = rhs.m_count;
        // Scaling.
        m_ncont      = rhs.m_ncont;
        m_scale      = rhs.m_scale;
        m_contfactor = rhs.m_contfactor;
        // Doubling.
        m_ndble      = rhs.m_ndble;
        m_dbleactive = rhs.m_dbleactive;
        m_dblecutoff = rhs.m_dblecutoff;
        m_dblelimit  = rhs.m_dblelimit;
        m_dbleslack  = rhs.m_dbleslack;
        m_dbleon     = rhs.m_dbleon;

        // Copy particle vector.
        for (const_iterator ip=rhs.begin(); ip!=rhs.end(); ++ip) {
            m_particles.push_back((*ip)->Clone());
        }

        // Copy binary tree.
        m_tree = rhs.m_tree;
        m_sums = rhs.m_sums;
    }
    return *this;
}


// INITIALISATION.

// Initialises the ensemble to the given size.  Any particles currently
// in the ensemble will be destroyed.
void Ensemble::Initialise(unsigned int capacity)
{
    Clear(); // Clear current ensemble.
    
    // Calculate nearest power of 2 capacity.  Ensemble capacity must be a power
    // of 2.  This constraint is due to the binary tree implementation.
    real rl    = log((real)capacity) / log(2.0);
    m_levels   = (int)(rl + 0.5);
    m_capacity = (int)pow(2.0, (int)m_levels);
    m_halfcap  = capacity / 2;
    m_count    = 0;

    // Reserve memory for ensemble.
    m_particles.assign(capacity, NULL);
    m_tree.resize(capacity);
    
    // Set all tree nodes to correct size and link up tree.
    unsigned int i, j;
    m_tree[0].Parent = NULL;
    for (i=0,j=1; i<m_halfcap-1; i++,j=(2*i)+1) {
        m_tree[i].Left = &m_tree.at(j);
        m_tree[i].Left->Parent = &(m_tree[i]);
        m_tree[i].Right = &m_tree.at(j+1);
        m_tree[i].Right->Parent = &(m_tree[i]);
    }
    for (i=m_halfcap-1; i<m_capacity; i++) {
        m_tree[i].Left  = NULL;
        m_tree[i].Right = NULL;
    }

    // Initialise scaling.
    m_ncont      = 0;
    m_scale      = 1.0;
    m_contfactor = (real)(m_capacity-1) / (real)(m_capacity);

    // Initialise doubling.
    m_ndble      = 0;
    m_dbleon     = true;
    m_dbleactive = false;
    m_dblecutoff = 3 * m_capacity / 4;
    m_dblelimit  = m_halfcap - (unsigned int)pow(2.0, (int)((m_levels-5)>0 ? m_levels-5 : 0));
    m_dbleslack  = (unsigned int)pow(2.0, (int)((m_levels-5)>0 ? m_levels-5 : 0));
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

// Adds the particle to the ensemble.  Takes control of particle
// destruction.
int Ensemble::Add(Particle &sp)
{
    // Check for doubling activation.
    m_dbleactive = m_dbleactive || (m_count >= m_dblecutoff-1);

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
    }

    if (i < 0) {
        // We are adding a new particle.
        m_particles[i=m_count++] = &sp;
    } else if (i < (int)m_capacity) {
        // Replace an existing particle (if i=m_capacity) then
        // we are removing the new particle, so just ignore it.
        ReplaceParticle(i, sp);
    }

    // Now we must recalculate the tree by inserting the particle properties
    // into the bottom row and calculating up.
    int j = treeIndex(i);
    if (isLeftBranch(i)) {
        m_tree[j].LeftData = *m_particles[i];
    } else {
        m_tree[j].RightData = *m_particles[i];
    }
    ascendingRecalc(j);

    return i;
}

// Removes the particle at the given index, if the index is
// valid.
void Ensemble::Remove(unsigned int i)
{
    // Check that particle index is valid.
    if (i<m_count-1) {
        // First delete particle from memory, then
        // overwrite it with the last particle, which
        // is subsequently removed from the vector.
        delete m_particles[i];
        m_particles[i] = m_particles[m_count];
        m_particles[m_count] = NULL;
        --m_count;

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
        delete m_particles[i];
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

    iterator i=begin(), k=end()-1;
    unsigned int ix=0, tx=0, n=0;
    bool found=false;

    while (i<k) {
        found = false;

        // Locate next invalid particle in list.
        while ((i!=k) && (*i)->IsValid()) {i++; ix++;}

        // Locate next valid particle from end of list.
        while ((k!=i) && (!(*k)->IsValid())) k--;

        if (i<k) {            
            delete *i; // Delete invalid particle from memory.
            *i = *k;   // Put valid particle at index i.
            n++;       // Increment number of invalid particles.
            *k = NULL; // Clear duplicate pointer.
            --m_count; // Decrrement number of particles.

            // Update binary tree.
            tx = treeIndex(ix);
            if (isLeftBranch(ix)) {
                m_tree[tx].LeftData = *(*i);
            } else {
                m_tree[tx].RightData = *(*i);
            }
            ascendingRecalc(tx);
        } else {
            // All invalid particles have been deleted.
            break;
        }
    }

    // If we removed too many invalid particles then we'll have to double.
    dble();
}

// Replaces the particle at the given index with the given particle,
// as long as the index is valid.
void Ensemble::ReplaceParticle(unsigned int i, Particle &sp)
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
    for (int i=0; i<(int)m_particles.size(); i++) {
        delete m_particles[i];
    }
    m_particles.clear();

    m_ncont = 0; // No contractions any more.

    vector<TreeNode>::iterator j;
    for(j=m_tree.begin(); j!=m_tree.end(); ++j) j->Clear();

    // Reset doubling.
    m_ndble      = 0;
    m_dbleactive = false;
}


// SELECTING PARTICLES.

// Returns the index of a uniformly selected particle
// from the ensemble.
int Ensemble::SelectParticle(void) const
{
    // Uniformly select a particle index.
    return irnd(0, m_count-1);
}

// Returns the index of a particle selected using the property
// weight given which refers to a basic property in the
// ParticleData class.
int Ensemble::SelectParticle(ParticleData::PropertyID id) const
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
    if (isp < m_count) {
        return isp;
    } else {
        // Chosen an non-existent particle.
        return -1;
    }
}

// Selects particle according to the particle property
// specified by the given model and the given property id
// of the model.
int Ensemble::SelectParticle(ModelType model_id, unsigned int id) const
{
    // This routine uses the binary tree to select a particle weighted
    // by a given particle property (by index).

    int isp=-1;

    // Calculate random number weighted by sum of desired property (wtid).
    real r = rnd();

    if (model_id == BasicModel_ID) {
        r *= (m_tree[0].LeftData.Property(static_cast<ParticleData::PropertyID>(id)) + 
              m_tree[0].RightData.Property(static_cast<ParticleData::PropertyID>(id)));
    } else {
        r *= (m_tree[0].LeftData.ModelCache(model_id)->Property(id) + 
              m_tree[0].RightData.ModelCache(model_id)->Property(id));
    }

    // Fall down the binary tree until reaching a base node.
    const TreeNode *n = &m_tree[0];
    int j=0;
    while(n->Left!=NULL) {
        if (r <= n->LeftData.ModelCache(model_id)->Property(id)) {
            n = n->Left;
            j = (2*j) + 1;
        } else {
            r -= n->LeftData.ModelCache(model_id)->Property(id);
            n = n->Right;
            j = (2*j) + 2;
        }
    }
    // Last level!
    if (r <=n->LeftData.ModelCache(model_id)->Property(id)) {
        isp = (2*j) + 2 - m_capacity;
    } else {
        isp = (2*j) + 3 - m_capacity;
    }

    // One final check that the index we've chosen is valid.
    if (isp < m_count) {
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
}


// GET SUMS OF PROPERTIES.

// Returns a ParticleData object which contains property
// sums for all particles in the ensemble.
const ParticleData &Ensemble::GetSums(void) const
{
	m_sums = m_tree[0].LeftData + m_tree[0].RightData;
    return m_sums;
}

// Returns the sum of a property in the ParticleData class
// over all particles.
real Ensemble::GetSum(ParticleData::PropertyID id) const
{
    return m_tree[0].LeftData.Property(id) + m_tree[0].RightData.Property(id);
}

// Returns the sum of one particle property with the given index 
// from the given model from the binary tree.
real Ensemble::GetSum(ModelType model_id, unsigned int id) const
{
    if (model_id == BasicModel_ID) {
        return GetSum(static_cast<ParticleData::PropertyID>(id));
    } else {
        return m_tree[0].LeftData.ModelCache(model_id)->Property(id) + 
               m_tree[0].RightData.ModelCache(model_id)->Property(id);
    }
}


// UPDATE ENSEMBLE.

// Updates the ensemble tree from the given particle index.
void Ensemble::Update(unsigned int i)
{
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

// Updates the ensemble tree completely.
void Ensemble::Update()
{
    // This flavour updates the whole binary tree.

    bool odd = true;
    iterator i;
    unsigned int j = treeIndex(0);
    for (i=begin(); i!=end(); i++) {
        if (odd) {
            m_tree[0].LeftData = *(*i);
        } else {
            m_tree[0].RightData = *(*i);
            ascendingRecalc(j);
            j++;
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

    if (i<m_count) {
        // Get the node at the bottom of the branch.
        TreeNode *n = &m_tree[i];

        // Climb up the tree until the root node, 
        // summing up the properties.
        while (n->Parent != NULL) {
            n = n->Parent;
            n->LeftData = n->Left->LeftData + n->Left->RightData;
            n->RightData = n->Right->LeftData + n->Right->RightData;
        }
    }
}

void Ensemble::dble()
{
    // The doubling algorithm is activated if the number of particles
    // in the ensemble falls below half capacity.  It copies the whole particle
    // list and changes the scaling factor to keep it consistent.  Once the
    // ensemble is back above half full, the routine updates the binary tree.

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
            left = isLeftBranch(n);

            // Copy particles.
            for (isp=0,j=n; isp!=n; isp++,j++,left=!left) {
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
            for (inode=(m_tree.begin()+m_halfcap-2); inode!=m_tree.begin(); inode--) {
                inode->LeftData = inode->Left->LeftData + inode->Left->RightData;
                inode->RightData = inode->Right->LeftData + inode->Right->RightData;
            }

            // Update scaling.
            m_ndble++;
            n *= 2;
        }
    }
}
