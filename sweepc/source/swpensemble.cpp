#include "swpensemble.h"
#include "rng.h"
#include <cmath>

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Ensemble::Ensemble(void)
{
    // On creation set all variables to default values.
    m_levels = 0;
    m_capacity = 0;
    m_halfcap = 0;
    m_ncont = 0;
    m_scale = 1.0;
    m_contfactor = 0;
    m_ndble = 0;
    m_dbleactive = false;
    m_dblecutoff = 0;
    m_dblelimit  = 0;
    m_dbleslack  = 0;
}

// Initialising constructor.
Ensemble::Ensemble(const unsigned int count, const unsigned int nprop)
{
    // Call initialisation routine.
    Initialise(count, nprop);
}

// Copy contructor.
Ensemble::Ensemble(const Sweep::Ensemble &ens)
{
    // Easy to copy things.
    m_levels = ens.m_levels;
    m_capacity = ens.m_capacity;
    m_halfcap = ens.m_halfcap;
    m_ncont = ens.m_ncont;
    m_scale = ens.m_scale;
    m_contfactor = ens.m_contfactor;
    m_ndble = ens.m_ndble;
    m_dbleactive = ens.m_dbleactive;
    m_dblecutoff = ens.m_dblecutoff;
    m_dblelimit  = ens.m_dblelimit;
    m_dbleslack  = ens.m_dbleslack;
    m_dbleon = ens.m_dbleon;

    // Particle vector.
    const_iterator ip;
    for (ip=ens.begin(); ip!=ens.end(); ip++) {
        m_particles.push_back(*(*ip));
    }

    // Binary tree.
    m_tree = ens.m_tree;
}

// Destructor.
Ensemble::~Ensemble(void)
{
    // Clear the ensemble.
    Destroy();
}

void Ensemble::Initialise(const unsigned int capacity, const unsigned int nprop)
{
    // Clear current ensemble.
    Destroy();
    
    // Calculate nearest power of 2 capacity.  Ensemble capacity must be a powe
    // of 2.  This constraint is due to the binary tree.
    real rl = log((real)capacity) / log(2.0);
    m_levels = (int)(rl + 0.5);
    m_capacity = (int)pow(2.0, (int)m_levels);
    m_halfcap = capacity / 2;

    // Reserve memory for ensemble.
    m_particles.reserve(capacity);
    m_tree.resize(capacity);
    
    // Set all tree nodes to correct size and link up tree.
    unsigned int i, j;
    m_tree[0].Parent = NULL;
    for (i=0,j=1; i<m_halfcap-1; i++,j=(2*i)+1) {
        m_tree[i].SetSize(nprop);
        m_tree[i].m_left = &m_tree.at(j);
        m_tree[i].m_left->Parent = &(m_tree[i]);
        m_tree[i].m_right = &m_tree.at(j+1);
        m_tree[i].m_right->Parent = &(m_tree[i]);
    }
    for (i=m_halfcap-1; i<m_capacity; i++) {
        m_tree[i].SetSize(nprop);
        m_tree[i].m_left = NULL;
        m_tree[i].m_right = NULL;
    }

    // Initialise scaling.
    m_ncont = 0;
    m_scale = 1.0;
    m_contfactor = real(m_capacity-1) / real(m_capacity);

    // Initialise doubling.
    m_ndble      = 0;
    m_dbleon     = true;
    m_dbleactive = false;
    m_dblecutoff = 3 * m_capacity / 4;
    m_dblelimit  = m_halfcap - (unsigned int)pow(2.0, (int)((m_levels-5)>0 ? m_levels-5 : 0));
    m_dbleslack  = (unsigned int)pow(2.0, (int)((m_levels-5)>0 ? m_levels-5 : 0));
}

void Ensemble::Destroy(void)
{
    // Delete particles from memory and delete vectors.
    vector<Particle*>::iterator i;
    for (i=m_particles.begin(); i!=m_particles.end(); i++) {
        delete (*i);
        (*i) = NULL;
    }
    m_particles.clear();
    m_tree.clear();

    // Reset variables.
    m_levels = 0;
    m_capacity = 0;
    m_halfcap = 0;
    m_ncont = 0;
    m_scale = 1.0;
    m_contfactor = 0;
    m_ndble = 0;
    m_dbleactive = false;
    m_dblecutoff = 0;
    m_dblelimit  = 0;
    m_dbleslack  = 0;
}

Particle *Ensemble::GetParticle(const unsigned int i) const
{
    // Check that the index in within range, then return the particle.
    if (i < (unsigned int)m_particles.size()) {
        return m_particles[i];
    } else {
        return NULL;
    }
}

int Ensemble::AddParticle(Particle &sp)
{
    // Check for doubling activation.
    m_dbleactive = m_dbleactive || ((unsigned int)m_particles.size() >= m_dblecutoff-1);

    // Check ensemble for space, if there is not enough space then need
    // to generate some by contracting the ensemble.
    int i=-1;
    if ((int)m_particles.size() < m_capacity) {
        // There is space in the tree for a new particle.
        i = -1;
    } else {
        // We must contract the ensemble to accomodate a new particle.
        i = irnd(0, m_capacity-1);
        m_ncont++;
    }

    if (i < 0) {
        // We are adding a new particle.
        m_particles.push_back(&sp);
        i = (int)m_particles.size() - 1;
    } else if (i < (int)m_capacity) {
        // Replace an existing particle (if i=m_capacity) then
        // we are removing the new particle, so just ignore it.
        ReplaceParticle(i,sp);
    }

    // Now we must recalculate the tree by inserting the particle properties
    // into the bottom row and calculating up.
    int j = TreeIndex(i);
    if (IsLeftBranch(i)) {
		//m_tree[j].Resize(ncache  SORT THIS OUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//
		//
		//
		//
		//
		//
		//
		//
         m_particles[i]->GetProperties(m_tree[j].m_leftsum);
    } else {
        m_particles[i]->GetProperties(m_tree[j].m_rightsum);
    }
    AscendingRecalc(j);

    return i;
}

void Ensemble::RemoveParticle(const unsigned int i)
{
    // Check that particle index is valid.
    if (i<(unsigned int)m_particles.size()-1) {
        // First delete particle from memory, then
        // overwrite it with the last particle, which
        // is subsequently removed from the vector.
        delete m_particles[i];
        m_particles[i] = *(m_particles.end()-1);
        m_particles.erase(m_particles.end()-1);

        // Recalculate the tree up from leaf.
        int k = TreeIndex(i);
        if (IsLeftBranch(i)) {
            m_particles[i]->GetProperties(m_tree[k].m_leftsum);
        } else {
            m_particles[i]->GetProperties(m_tree[k].m_rightsum);
        }
        AscendingRecalc(k);

        // Recalculate tree up from last position (clear it first).
        unsigned int j = (unsigned int)m_particles.size();
        k = TreeIndex(j);
        if (IsLeftBranch(j)) {
            m_tree[k].m_leftsum.Clear();
        } else {
            m_tree[k].m_rightsum.Clear();
        }
        AscendingRecalc(k);
    } else if (i==(unsigned int)m_particles.size()-1) {
        // This is the last particle in the ensemble, we don't
        // need to swap it with another, just delete it.

        // Erase particle.
        delete m_particles[i];
        m_particles.erase(m_particles.end()-1);

        // Recalculate tree up from last position (clear it first).
        int k = TreeIndex(i);
        if (IsLeftBranch(i)) {
            m_tree[k].m_leftsum.Clear();
        } else {
            m_tree[k].m_rightsum.Clear();
        }
        AscendingRecalc(k);
    }

    // Particle removal might reduce the particle count sufficiently to require
    // particle doubling.
    Double();
}

void Ensemble::RemoveInvalids(void)
{
    // This function loops forward through the list finding invalid
    // particles and backwards finding valid particles.  Once an invalid
    // and a valid particle are found they are swapped.  This results in
    // all the invalid particles being collected at the end of the vector,
    // from where they can be deleted.

    iterator i=begin(), k=end()-1;
    unsigned int ix=0, tx=0, n=0;
    bool found;

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

            // Update binary tree.
            tx = TreeIndex(ix);
            if (IsLeftBranch(ix)) {
                (*i)->GetProperties(m_tree[tx].m_leftsum);
            } else {
                (*i)->GetProperties(m_tree[tx].m_rightsum);
            }
            AscendingRecalc(tx);
        } else {
            // All invalid particles have been deleted.
            break;
        }
    }

    // Erase duplicate pointers at end of vector.
    if (n>0) m_particles.erase(end()-n, end());

    // If we removed too many invalid particles then we'll have to double.
    Double();
}

void Ensemble::ReplaceParticle(const unsigned int i, Particle &sp)
{
    // Check index is within range.
    if (i<(int)m_particles.size()) {
        // First delete current particle, then
        // set pointer to new particle.
        delete m_particles[i];
        m_particles[i] = &sp;

        // Recalculate tree up from branch.
        int j = TreeIndex(i);
        if (IsLeftBranch(i)) {
            m_particles[i]->GetProperties(m_tree[j].m_leftsum);
        } else {
            m_particles[i]->GetProperties(m_tree[j].m_rightsum);
        }
        AscendingRecalc(j);
    }
}

void Ensemble::Clear()
{
    // Delete particles from memory and delete vectors.
    for (int i=0; i<(int)m_particles.size(); i++) {
        delete m_particles[i];
    }
    m_particles.clear();

    m_ncont = 0; // No contractions any more.

    vector<TreeNode>::iterator i;
    for(i=m_tree.begin(); i!=m_tree.end(); i++) (*i).ClearNode();

    // Reset doubling.
    m_ndble = 0;
    m_dbleactive = false;
}

int Ensemble::SelectParticle(void) const
{
    // Uniformly select a particle index.
    return irnd(0, (int)m_particles.size()-1);
}

int Ensemble::SelectParticle(const int wtid) const
{
    // This routine uses the binary tree to select a particle weighted
    // by a given particle property (by index).

    int isp=-1;

    if ((wtid>=0) && (wtid<(int)m_tree[0].m_leftsum.Size())) {
        // Weight ID is valid, so choose using binary tree.

        // Calculate random number weighted by sum of desired property (wtid).
        real r = rnd() * (m_tree[0].m_leftsum[wtid] + m_tree[0].m_rightsum[wtid]);

        // Fall down the binary tree until reaching a base node.
        const TreeNode *n = &m_tree[0];
        int j=0;
        while(n->m_left!=NULL) {
            if (r <=n->m_leftsum[wtid]) {
                n = n->m_left;
                j = (2*j) + 1;
            } else {
                r -= n->m_leftsum[wtid];
                n = n->m_right;
                j = (2*j) + 2;
            }
        }
        // last level!
        if (r <=n->m_leftsum[wtid]) {
            isp = (2*j) + 2 - m_capacity;
        } else {
            isp = (2*j) + 3 - m_capacity;
        }
    } else {
        // Select uniformly if wtid is invalid.
        isp = SelectParticle();
    }

    // One final check that the inex we've chosen is valid.
    if (isp >= (int)m_particles.size()) {
        // Chosen an non-existent particle.
        return -1;
    } else {
        return isp;
    }
}

real Ensemble::Scaling() const
{
    // The scaling factor includes the contraction term and the doubling term.
    return m_scale * pow(m_contfactor, (double)m_ncont) * pow(2.0,(double)m_ndble);
}

void Ensemble::ResetScaling()
{
    m_ncont = 0;
    m_ndble = 0;
}

void Ensemble::GetSums(std::vector<real> &sums) const
{
	sums.assign(m_tree[0].m_leftsum.GetCache().begin(), m_tree[0].m_leftsum.GetCache().end());
	vector<real>::iterator i;
	vector<real>::const_iterator j;
	for(i=sums.begin(),j=m_tree[0].m_rightsum.GetCache().begin(); i!=sums.end(); i++,j++)
		*i += *j;
}

real Ensemble::GetSum(unsigned int i) const
{
    if (i<(unsigned int)m_tree[0].m_leftsum.Size()) {
        return m_tree[0].m_leftsum[i] + m_tree[0].m_rightsum[i];
    } else {
        return 0.0;
    }
}


void Ensemble::Update(const unsigned int i)
{
    // Get tree index of this particle.
    unsigned int j = TreeIndex(i);

    // Update binary tree at this index.
    if (IsLeftBranch(i)) {
        m_particles[i]->GetProperties(m_tree[j].m_leftsum);
    } else {
        m_particles[i]->GetProperties(m_tree[j].m_rightsum);
    }
    AscendingRecalc(j);
}

void Ensemble::Update()
{
    // This flavour updates the whole binary tree.

    bool odd = true;
    iterator i;
    unsigned int j = TreeIndex(0);
    for (i=begin(); i!=end(); i++) {
        if (odd) {
            (*i)->GetProperties(m_tree[j].m_leftsum);
        } else {
            (*i)->GetProperties(m_tree[j].m_rightsum);
            AscendingRecalc(j);
            j++;
        }
        odd = !odd;
    }

    // Need to do one last recalc if there are an odd number of particles.
    if (!odd) {
        AscendingRecalc(j);
    }
}

/* Protected functions. */

void Ensemble::AscendingRecalc(const unsigned int i)
{
    // This function starts at the bottom of the binary tree and works
    // upwards recalculating the sums of all particle properties under
    // the nodes.

    vector<real>::iterator j, k, l;
    if ((i>=0) && (i<(int)m_tree.size())) {
        // Get the node at the bottom of the branch.
        TreeNode *n = &m_tree[i];

        // Climb up the tree until the root node, summing up the properties.
        while (n->Parent != NULL) {
            n = n->Parent;
            for(j=n->m_leftsum.GetCache().begin(),k=n->m_left->m_leftsum.GetCache().begin(),l=n->m_left->m_rightsum.GetCache().begin(); 
                j!=n->m_leftsum.GetCache().end(); j++,k++,l++) {
                *j = *k + *l;
            }
            for(j=n->m_rightsum.GetCache().begin(),k=n->m_right->m_leftsum.GetCache().begin(),l=n->m_right->m_rightsum.GetCache().begin(); 
                j!=n->m_rightsum.GetCache().end(); j++,k++,l++) {
                *j = *k + *l;
            }
        }
    }
}

void Ensemble::Double()
{
    // The doubling algorithm is activated if the number of particles
    // in the ensemble falls below half capacity.  It copies the whole particle
    // list and changes the scaling factor to keep it consistent.  Once the
    // ensemble is back above half full, the routine updates the binary tree.

    // Check that doubling is on and the activation condition has been met.
    if (m_dbleon && m_dbleactive) {
        Particle* sp;
        bool left;
        unsigned int isp, j, k;
        unsigned int n=(int)m_particles.size();
        vector<TreeNode>::iterator inode;
        vector<real>::iterator iterp, iterl, iterr;

        // Continue while there are too few particles in the ensemble.
        while (n<=m_dblelimit) {
            left = IsLeftBranch(n);

            // Copy particles.
            for (isp=0,j=n; isp<n; isp++,j++,left=!left) {
                // Create a copy of a particle and add it to the ensemble.
                sp = &m_particles[isp]->CreateCopy();
                m_particles.push_back(sp);
                // Put properties into bottom of tree.
                k = TreeIndex(j);
                if (left) {
                    sp->GetProperties(m_tree[k].m_leftsum);
                } else {
                    sp->GetProperties(m_tree[k].m_rightsum);
                }
            }

            // Update the tree from second last level up.
            for (inode=(m_tree.begin()+m_halfcap-2); inode!=m_tree.begin(); inode--) {
                for(iterp=(*inode).m_leftsum.GetCache().begin(),
                    iterl=(*inode).m_left->m_leftsum.GetCache().begin(),
                    iterr=(*inode).m_left->m_rightsum.GetCache().begin(); 
                    iterp!=(*inode).m_leftsum.GetCache().end(); iterp++,iterl++,iterr++) {
                    *iterp = *iterl + *iterr;
                }
                for(iterp=(*inode).m_rightsum.GetCache().begin(),
                    iterl=(*inode).m_right->m_leftsum.GetCache().begin(),
                    iterr=(*inode).m_right->m_rightsum.GetCache().begin(); 
                    iterp!=(*inode).m_rightsum.GetCache().end(); iterp++,iterl++,iterr++) {
                    *iterp = *iterl + *iterr;
                }
            }

            // Update scaling.
            m_ndble++;
            n *= 2;
        }
    }
}
