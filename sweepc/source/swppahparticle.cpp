#include "swppahparticle.h"
#include "swppahmodel.h"

using namespace Sweep;

PAHParticle::PAHParticle(void)
{
    m_ed = NULL;
    m_zz = NULL;
    m_r5 = NULL;
    m_ac = NULL;
    m_hl = NULL;
}

PAHParticle::~PAHParticle(void)
{
}

void PAHParticle::Initialise(vector<Component*> &components, const unsigned int nvals)
{
    // Initialise default particle components.
    DefaultParticle::Initialise(components, nvals);

    // Get references to PAH structure values.
    m_ed = &m_values.at(PAHModel::iED);
    m_zz = &m_values.at(PAHModel::iZZ);
    m_r5 = &m_values.at(PAHModel::iR5);
    m_ac = &m_values.at(PAHModel::iAC);
    m_hl = &m_values.at(PAHModel::iHL);
}

void PAHParticle::PerformEdgeRingGrowth(unsigned int &n)
{
    // Edge ring growth occurs when 2 C2H2s add to a PAH edge and close
    // to form a new R6 ring.  The rules for this have been hardcoded
    // here.  Note that this routine does not alter the composition of
    // the particle, so C and H counts have to be changed in the standard
    // way using the mechanism file.

    // Check that there are sufficient edges to perform n additions.
    n = min(n, (unsigned int)m_ed);

    // Edge growth removes one edge, but adds three more, therefore add
    // two edges to particle.
    m_ed += (2*n);

    // Additionally, two neighbouring sites must be incremented.
    IncrementSites(2*n);
}

void PAHParticle::PerformArmchairRingGrowth(unsigned int &n)
{
    // Armchair ring growth occurs when 1 C2H2 adds to a PAH armchair
    // to form a new R6 ring.  The rules for this have been hardcoded
    // here.  Note that this routine does not alter the composition of
    // the particle, so C and H counts have to be changed in the standard
    // way using the mechanism file.

    // Check that there are sufficient armchairs to perform n additions.
    n = min(n, (unsigned int)m_ac);

    // Armchair growth removes one armchair and adds one edge.
    m_ac -= n;
    m_ed += n;

    // Additionally, two neighbouring sites must be incremented.
    IncrementSites(2*n);
}

void PAHParticle::PerformR5Addition(unsigned int &n)
{
    // R5 addition occurs when 1 C2H2 adds to a PAH zig-zag to form a 
    // new R5 ring.  The rules for this have been hardcoded here.  Note
    // that this routine does not alter the composition of the particle, 
    // so C and H counts have to be changed in the standard way using 
    // the mechanism file.  No neighbouring sites are altered by R5 addition.

    // Check that there are sufficient zig-zags to perform n additions.
    n = min(n, (unsigned int)m_zz);

    // R5 addition replaces one zig-zag with an R5.
    m_zz -= n;
    m_r5 += n;
}

void PAHParticle::PerformR5EdgeConversion(unsigned int &n)
{
    // R5 edge conversion occurs when 1 C2H2 adds to a PAH edge next to an
    // R5 ring, which then convert to a new R6 ring.  The rules for this have
    // been hardcoded here.  Note that this routine does not alter the 
    // composition of the particle, so C and H counts have to be changed in 
    // the standard way using the mechanism file.

    // Check that there are sufficient R5s and edges to perform n additions.
    n = min(n, (unsigned int)m_ed);
    n = min(n, (unsigned int)m_r5);

    // R5 edge conversion removes one R5 and one edge, but adds one armchair and 
    // three edges.
    m_r5 -= n;
    m_ac += n;
    m_ed += (2*n);

    // Additionally, one neighbouring site must be incremented.
    IncrementSites(n);
}

void PAHParticle::PerformR5ArmchairConversion(unsigned int &n)
{
    // R5 armchair conversion occurs when a R5 ring collides with a PAH
    // armchair, which then convert to a new R6 ring.  The rules for this have
    // been hardcoded here.  Note that this routine does not alter the 
    // composition of the particle, so C and H counts have to be changed in 
    // the standard way using the mechanism file.

    // Check that there are sufficient R5s and edges to perform n additions.
    n = min(n, (unsigned int)m_ed);
    n = min(n, (unsigned int)m_r5);

    // R5 armchair conversion removes one R5 and one armchair, but adds one
    // armchair and one edge.
    m_r5 -= n;
    m_ed += n;

    // Additionally, one neighbouring site must be incremented.
    IncrementSites(n);
}

void PAHParticle::IncrementSites(unsigned int n)
{
    // PAH Sites are incremented when a growth process occurs or when a
    // R5 ring converts at an armchair site.  If sites are considered by
    // the number of carbon atoms which they contain, then site incrementing
    // increases this count by one.

    // Count the total number of sites which can be incremented.
    real sitetotal = *m_ed + *m_zz + *m_ac + *m_hl;

    unsigned int i;
    real rand;

    // Loop over number of sites to be incremened.
    for (i=0; i<n; i++)
    {
        // Generate a random deviate.
        rand = rnd() * sitetotal;

        // Use DIV algorithm to select a site to increment.
        rand = rand - *m_ed;
        if (rand <= 0.0) {
            // Increment an edge to a zig-zag.
            m_ed -= 1;
            m_zz += 1;
        } else {
            rand = rand - *m_zz;
            if (rand <= 0.0) {
                // Increment a zig-zag to an armchair.
                m_zz -= 1;
                m_ac += 1;
            } else {
                rand = rand - *m_ac;
                if (rand <= 0.0) {
                    // Increment an armchair to a hole.
                    m_ac -= 1;
                    m_hl += 1;
                } else {
                    // Increment a hole to nothing (assumption).
                    m_hl -= 1;
                    sitetotal -= 1;
                }
            }
        }
    }
}

void PAHParticle::DecrementSites(unsigned int n)
{
    // PAH Sites are decremented when a shrinking process occurs such as 
    // oxidation.  If sites are considered by the number of carbon atoms 
    // which they contain, then site decrementing decreases this count by one.

    // Count the total number of sites which can be decremented.
    real sitetotal = *m_zz + *m_ac + *m_hl;

    unsigned int i;
    real rand;

    // Loop over number of sites to be decremented.
    for (i=0; i<n; i++)
    {
        // Generate a random deviate.
        rand = rnd() * sitetotal;

        // Use DIV algorithm to select a site to decrement.
        rand = rand - *m_zz;
        if (rand <= 0.0) {
            // Decrement a zig-zag to an edge.
            m_zz -= 1;
            m_ed += 1;
        } else {
            rand = rand - *m_ac;
            if (rand <= 0.0) {
                // Decrement an armchair to a zig-zag.
                m_ac -= 1;
                m_zz += 1;
            } else {
                // Decrement a hole to an armchair.
                m_hl -= 1;
                m_ac += 1;
            }
        }
    }
}