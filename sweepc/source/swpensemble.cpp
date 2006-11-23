#include "swpensemble.h"
#include "rng.h"
#include <cmath>

using namespace Sweep;

Ensemble::Ensemble(void)
{
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

Ensemble::Ensemble(const unsigned int count, const unsigned int nprop)
{
    Initialise(count, nprop);
}

Ensemble::~Ensemble(void)
{
    Destroy();
}

void Ensemble::Initialise(const unsigned int capacity, const unsigned int nprop)
{
    // Clear current ensemble.
    Destroy();
    
    // Calculate nearest power of 2 capacity.
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
        m_tree[i].Left = &m_tree.at(j);
        m_tree[i].Left->Parent = &(m_tree[i]);
        m_tree[i].Right = &m_tree.at(j+1);
        m_tree[i].Right->Parent = &(m_tree[i]);
    }
    for (i=m_halfcap-1; i<m_capacity; i++) {
        m_tree[i].SetSize(nprop);
        m_tree[i].Left = NULL;
        m_tree[i].Right = NULL;
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
    vector<DefaultParticle*>::iterator i;
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

DefaultParticle *Ensemble::GetParticle(const unsigned int i) const
{
    if (i < (unsigned int)m_particles.size()) {
        return m_particles[i];
    } else {
        return NULL;
    }
}

int Ensemble::AddParticle(DefaultParticle &sp)
{
    m_dbleactive = m_dbleactive || ((unsigned int)m_particles.size() >= m_dblecutoff-1);
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

    // Now we must recalculate the tree.
    int j = TreeIndex(i);
    if (IsLeftBranch(i)) {
         m_particles[i]->GetProperties(m_tree[j].LeftSum);
    } else {
        m_particles[i]->GetProperties(m_tree[j].RightSum);
    }
    AscendingRecalc(j);

    return i;
}

void Ensemble::RemoveParticle(const unsigned int i)
{
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
            m_particles[i]->GetProperties(m_tree[k].LeftSum);
        } else {
            m_particles[i]->GetProperties(m_tree[k].RightSum);
        }
        AscendingRecalc(k);

        // Recalculate tree up from last position (clear it first).
        unsigned int j = (unsigned int)m_particles.size();
        k = TreeIndex(j);
        if (IsLeftBranch(j)) {
            m_tree[k].LeftSum.assign(m_tree[k].LeftSum.size(), 0.0);
        } else {
            m_tree[k].RightSum.assign(m_tree[k].RightSum.size(), 0.0);
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
            m_tree[k].LeftSum.assign(m_tree[k].LeftSum.size(), 0.0);
        } else {
            m_tree[k].RightSum.assign(m_tree[k].RightSum.size(), 0.0);
        }
        AscendingRecalc(k);
    }
    Double();
}

void Ensemble::RemoveInvalids(void)
{
    iterator i=begin(), k=end()-1;
    unsigned int ix=0, tx=0;
    bool found;

    while (i!=k) {
        found = false;

        // Keep replacing particle with one from the end of
        // the array until we find one that is valid.
        while ((i!=k) && (!(*i)->IsValid())) {
            found = true;
            delete *i;
            *i = *k;
            k--;
        }

        if (found && (i!=k)) {
            // Now *i is a valid particle so recalculate the tree.
            tx = TreeIndex(ix);
            if (IsLeftBranch(ix)) {
                (*i)->GetProperties(m_tree[tx].LeftSum);
            } else {
                (*i)->GetProperties(m_tree[tx].RightSum);
            }
            AscendingRecalc(tx);
        }

        // Next particle.
        i++; ix++;
    }
}

void Ensemble::ReplaceParticle(const unsigned int i, DefaultParticle &sp)
{
    if (i<(int)m_particles.size()) {
        // First delete current particle, then
        // set pointer to new particle.
        delete m_particles[i];
        m_particles[i] = &sp;

        // Recalculate tree up from branch.
        int j = TreeIndex(i);
        if (IsLeftBranch(i)) {
            m_particles[i]->GetProperties(m_tree[j].LeftSum);
        } else {
            m_particles[i]->GetProperties(m_tree[j].RightSum);
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

    m_ncont = 0;

    vector<NODE>::iterator i;
    for(i=m_tree.begin(); i!=m_tree.end(); i++) (*i).Clear();
    m_ndble = 0;
}

int Ensemble::SelectParticle(void) const
{
    return irnd(0, (int)m_particles.size()-1);
}

int Ensemble::SelectParticle(const int wtid) const
{
    int isp=-1;

    if ((wtid>=0) && (wtid<(int)m_tree[0].LeftSum.size())) {
        // Weight ID is valid, so choose using binary tree.

        // Calculate random number weighted by sum of desired property (wtid).
        real r = rnd() * (m_tree[0].LeftSum[wtid] + m_tree[0].RightSum[wtid]);

        // Fall down the binary tree until reaching a base node.
        const NODE *n = &m_tree[0];
        int j=0;
        while(n->Left!=NULL) {
            if (r <=n->LeftSum[wtid]) {
                n = n->Left;
                j = (2*j) + 1;
            } else {
                r -= n->LeftSum[wtid];
                n = n->Right;
                j = (2*j) + 2;
            }
        }
        // last level!
        if (r <=n->LeftSum[wtid]) {
            isp = (2*j) + 2 - m_capacity;
        } else {
            isp = (2*j) + 3 - m_capacity;
        }
    } else {
        // Select uniformly if wtid is invalid.
        isp = SelectParticle();
    }

    if (isp >= (int)m_particles.size()) {
        // Chosen an non-existent particle.
        return -1;
    } else {
        return isp;
    }
}

real Ensemble::Scaling() const
{
    return m_scale * pow(m_contfactor, (double)m_ncont) * pow(2.0,(double)m_ndble);
}

void Ensemble::GetSums(std::vector<real> &sums) const
{
    sums.assign(m_tree[0].LeftSum.begin(), m_tree[0].LeftSum.end());
    vector<real>::iterator i;
    vector<real>::const_iterator j;
    for(i=sums.begin(),j=m_tree[0].RightSum.begin(); i!=sums.end(); i++,j++)
        *i += *j;
}

real Ensemble::GetSum(unsigned int i) const
{
    if (i<(unsigned int)m_tree[0].LeftSum.size()) {
        return m_tree[0].LeftSum[i] + m_tree[0].RightSum[i];
    } else {
        return 0.0;
    }
}

/* Protected functions. */

void Ensemble::AscendingRecalc(const unsigned int i)
{
    vector<real>::iterator j, k, l;
    if ((i>=0) && (i<(int)m_tree.size())) {
        // Get the node at the bottom of the branch.
        NODE *n = &m_tree[i];

        // Climb up the until the root node, summing up the properties.
        while (n->Parent != NULL) {
            n = n->Parent;
            for(j=n->LeftSum.begin(),k=n->Left->LeftSum.begin(),l=n->Left->RightSum.begin(); 
                j!=n->LeftSum.end(); j++,k++,l++) {
                *j = *k + *l;
            }
            for(j=n->RightSum.begin(),k=n->Right->LeftSum.begin(),l=n->Right->RightSum.begin(); 
                j!=n->RightSum.end(); j++,k++,l++) {
                *j = *k + *l;
            }
        }
    }
}

void Ensemble::Double()
{
    if (m_dbleon && m_dbleactive) {
        DefaultParticle* sp;
        bool left;
        unsigned int isp, j, k;
        unsigned int n=(int)m_particles.size();
        vector<NODE>::iterator inode;
        vector<real>::iterator iterp, iterl, iterr;

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
                    sp->GetProperties(m_tree[k].LeftSum);
                } else {
                    sp->GetProperties(m_tree[k].RightSum);
                }
            }

            // Update the tree from second last level up.
            for (inode=(m_tree.begin()+m_halfcap-2); inode!=m_tree.begin(); inode--) {
                for(iterp=(*inode).LeftSum.begin(),
                    iterl=(*inode).Left->LeftSum.begin(),
                    iterr=(*inode).Left->RightSum.begin(); 
                    iterp!=(*inode).LeftSum.end(); iterp++,iterl++,iterr++) {
                    *iterp = *iterl + *iterr;
                }
                for(iterp=(*inode).RightSum.begin(),
                    iterl=(*inode).Right->LeftSum.begin(),
                    iterr=(*inode).Right->RightSum.begin(); 
                    iterp!=(*inode).RightSum.end(); iterp++,iterl++,iterr++) {
                    *iterp = *iterl + *iterr;
                }
            }

            // Update scaling.
            m_ndble++;
            n *= 2;
        }
    }
}

/* Binary tree functions. */

Ensemble::NODE::NODE(void)
{
    Left = NULL; Right = NULL;
}

Ensemble::NODE &Ensemble::NODE::operator=(const Ensemble::NODE &n)
{
    if (this == &n) return *this;
    LeftSum = n.LeftSum;
    RightSum = n.RightSum;
    return *this;
}

Ensemble::NODE &Ensemble::NODE::operator +=(const Ensemble::NODE &n)
{
    LeftSum.resize(max(LeftSum.size(),n.LeftSum.size()), 0.0);
    RightSum.resize(max(RightSum.size(),n.RightSum.size()), 0.0);
    vector<real>::iterator i;
    vector<real>::const_iterator j=n.LeftSum.begin();
    for (i=LeftSum.begin(); i!=LeftSum.end(); i++,j++) *i = *i + *j;
    for (i=RightSum.begin(),j = n.RightSum.begin(); i!=RightSum.end(); i++,j++) *i = *i + *j;
    return *this;
}

const Ensemble::NODE Ensemble::NODE::operator+(const Ensemble::NODE &n) const
{
    NODE lhs = *this;
    lhs += n;
    return lhs;
}

void Ensemble::NODE::Clear()
{
    LeftSum.assign(LeftSum.size(), 0.0);
    RightSum.assign(RightSum.size(), 0.0);
}

void Ensemble::NODE::SetSize(const unsigned int size)
{
    LeftSum.resize(size,0.0);
    RightSum.resize(size,0.0);
}