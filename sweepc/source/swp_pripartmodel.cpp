#include "swp_pripartmodel.h"
#include "swp_pripartdata.h"
#include "swp_particledata.h"
#include "swp_particle.h"
#include "swp_modeltype.h"
#include <vector>
#include <algorithm>

using namespace Sweep;
using namespace std;

// SINGLETON IMPLEMENTATION.

// Default constructor (private).
PriPartModel::PriPartModel(void)
{
    // Nothing to do here.
}

// Copy constructor.
PriPartModel::PriPartModel(const Sweep::PriPartModel &copy)
{
    // Nothing to do here.
}

// Default destructor.
PriPartModel::~PriPartModel(void)
{
}

// Assignment operator.
const PriPartModel &PriPartModel::operator=(const Sweep::PriPartModel &rhs)
{
    // Nothing to do here.
    return *this;
}

// Returns the one and only instance of the CoagModel class.
// If the instance hasn't been created then this is done
// also.
PriPartModel &PriPartModel::Instance(void)
{
    static PriPartModel inst;
    return inst;
}


// PARTICLE UPDATES.

void PriPartModel::UpdateParticle(Sweep::Particle &p, 
                               const Sweep::fvector &dcomp, 
                               const Sweep::fvector &dval) const
{
    // Currently only the first component is added to the primary
    // particles.
    PriPartModelData &cache = dynamic_cast<PriPartModelData&>(*p.ModelCache(PriPartModel_ID));
    vector<Primary> &pri = cache.Primaries();

    // Distribute the mass, weighted by surface area, starting with the
    // largest primary particle.
    distMass(pri, dcomp[0], (*p.Components())[0]);

    // Get current primary surface area and calculate required change in
    // surface area.
    real surf = 0.0;
    for (vector<Primary>::iterator i=pri.begin(); i!=pri.end(); i++) {
        surf += i->SurfaceArea();
    }
    real ds = p.SurfaceArea() - surf;

    // Now the difficult bit, probably need to destroy primaries to 
    // reduce surface area.  The methodology is to remove the smallest
    // primaries, adding their mass to the larger, until the total surface
    // area is less than or equal to the aggregate surface area.
    if (ds < 0.0) {
        while ((surf > p.SurfaceArea()) && (pri.size() > 1)) {
            // Calculate change in surface area if smallest primary
            // is removed.
            surf -= (pri.end()-1)->SurfaceArea();

            // Remove the smallest primary, but save its mass first.
            real dm = (pri.end()-1)->Mass();
            pri.erase(pri.end()-1);

            // Distribute mass weighted by surface area, starting with 
            // largest particle.
            distMass(pri, dm, (*p.Components())[0]);
        }
    }
}

void PriPartModel::UpdateParticle(Sweep::Particle &p, 
                               const Sweep::fvector &dcomp, 
                               const Sweep::fvector &dval, 
                               unsigned int n) const
{
    PriPartModelData &cache = dynamic_cast<PriPartModelData&>(*p.ModelCache(PriPartModel_ID));
    cache.Primaries();
}


// PARTICLE-PARTICLE COAGULATION.

void PriPartModel::CoagParticles(Particle &p1, const Particle &p2) const
{
    // Add primary lists together.
    vector<Primary> &pri1 = dynamic_cast<PriPartModelData&>(*(p1.ModelCache(PriPartModel_ID))).Primaries();
    const vector<Primary> &pri2 = dynamic_cast<const PriPartModelData&>(*(p2.ModelCache(PriPartModel_ID))).Primaries();
    for (vector<Primary>::const_iterator i=pri2.begin(); i!=pri2.end(); ++i) {
        pri1.push_back(*i);
    }

    // Sort the resultant list.
    sort(pri1.begin(), pri1.end());
}


// PROPERTY CALCULATION.

// Recalculates those model properties which are function of
// other particle properties.
void PriPartModel::UpdateCache(Sweep::ParticleData &p) const
{
    PriPartModelData &cache = dynamic_cast<PriPartModelData&>(*p.ModelCache(PriPartModel_ID));
}


// PRIVATE IMPLEMENTATION.

// Distributes mass over a vector of primaries.
void PriPartModel::distMass(std::vector<Primary> &pri, 
                            real dmass, 
                            const Component *const comp)
{
    // TODO: Complete this function.
}
