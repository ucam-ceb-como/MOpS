#include "swp_primary.h"
#include "swp_surfvol_primary.h"
#include "swp_pripart_primary.h"
#include "swp_submodel_type.h"
#include "swp_aggmodel_type.h"
#include "swp_model_factory.h"
#include "swp_surfvol_cache.h"
#include "rng.h"
#include <stdexcept>
#include <vector>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace Sweep::SubModels;
using namespace std;

// Default component is the first component.
unsigned int PriPartPrimary::m_icomp = 0;


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
PriPartPrimary::PriPartPrimary(void)
: m_totprisurf(0.0), m_avgpridiam(0.0)
{
}

// Initialising constructor.
PriPartPrimary::PriPartPrimary(real time, const Sweep::ParticleModel &model)
: SurfVolPrimary(time, model), m_totprisurf(0.0), m_avgpridiam(0.0)
{
}

// Copy constructor.
PriPartPrimary::PriPartPrimary(const PriPartPrimary &copy) 
{
    *this = copy;
}

// Stream-reading constructor.
PriPartPrimary::PriPartPrimary(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
PriPartPrimary::~PriPartPrimary()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator (Primary RHS).
PriPartPrimary &PriPartPrimary::operator=(const Primary &rhs)
{
    // Now check the type of the RHS to see if additional data
    // copying is required.
    if (rhs.AggID() == Spherical_ID) {
        // Use the SurfVolPrimary version of the operator in order
        // to set correctly the surface-volume parameters.
        SurfVolPrimary::operator=(rhs);
        
        // Additionally we have to set the pripart-list correctly with
        // a single primary particle.
        m_primaries.clear();
        vector<PriPart>::iterator i = m_primaries.insert(m_primaries.end(), PriPart());
        i->Monomers = (unsigned int)rhs.Composition(m_icomp);
        i->Surface  = calcSurf(i->Monomers);

        // Set bulk primary properties.
        m_totprisurf = i->Surface;
        m_avgpridiam = PriDiameter(0);
    } else {
        // Attempt to cast the primary to a SurfVolPrimary.  This will
        // throw an exception if the cast fails.  Note, the subsequent
        // SurfVolPrimary version of the operator= function will attempt
        // to cast to a PriPartPrimary, if it isn't a SurfVolPrimary.
        operator=(dynamic_cast<const SurfVolPrimary&>(rhs));
    }

    return *this;
}

// Assignment operator (SurfVolPrimary RHS).
PriPartPrimary &PriPartPrimary::operator=(const SurfVolPrimary &rhs)
{
    if (rhs.AggID() == SurfVol_ID) {
        // This IS a surface-volume primary.
        SurfVolPrimary::operator=(rhs);

        // Set the primary particle list to have approx. the
        // same monomer count and surface area as the RHS.
//        createPriList((unsigned int)rhs.Composition(m_icomp), rhs.SurfaceArea());

        // Add a single primary to the list with the correct mass.
        m_primaries.clear();
        vector<PriPart>::iterator i = m_primaries.insert(m_primaries.end(), PriPart());
        i->Monomers = (unsigned int)rhs.Composition(m_icomp);
        i->Surface  = calcSurf(i->Monomers);

        // Set bulk primary properties.
        m_totprisurf = i->Surface;
        m_avgpridiam = PriDiameter(0);
    } else {
        // Attempt to cast the primary to a PriPartPrimary.  This will
        // throw an exception if the cast fails.
        operator=(dynamic_cast<const PriPartPrimary&>(rhs));
    }

    return *this;
}

// Assignment operator (PriPartPrimary RHS).
PriPartPrimary &PriPartPrimary::operator=(const PriPartPrimary &rhs)
{
    // Assign the surface-volume properties.
    SurfVolPrimary::operator=(rhs);

    // Copy the pripart-list.
    m_primaries.assign(rhs.m_primaries.begin(), rhs.m_primaries.end());

    // Set bulk properties.
    m_totprisurf = rhs.m_totprisurf;
    m_avgpridiam = rhs.m_avgpridiam;

    return *this;
}

// Compound assignment (coagulation).
PriPartPrimary &PriPartPrimary::operator+=(const Primary &rhs)
{
    return PriPartPrimary::Coagulate(rhs);
}


// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType PriPartPrimary::AggID(void) const {return AggModels::PriPartList_ID;}

// Creates an aggregation data cache for this primary type.
AggModels::PriPartCache *const PriPartPrimary::CreateAggCache(ParticleCache &pcache) const
{
    PriPartCache *cache = 
        static_cast<PriPartCache*>(ModelFactory::CreateAggCache(AggModels::PriPartList_ID, pcache));
    if (cache != NULL) *cache = *this;
    return cache;
}


// PRIMARY PARTICLE LIST.

// Returns the number of primary particles in the list.
unsigned int PriPartPrimary::PriCount(void) const {return m_primaries.size();}

// Returns the number of monomers in the ith primary.
unsigned int PriPartPrimary::MonomerCount(unsigned int i) const
{
    if (i < m_primaries.size()) {
        return m_primaries[i].Monomers;
    } else {
        return 0;
    }
}

// Returns the diameter of the ith primary.
real PriPartPrimary::PriDiameter(unsigned int i) const
{
    real vol = PriVolume(i);
    return pow(vol * 6.0 / PI, ONE_THIRD);
}

// Returns the average primary particle diameter.
real PriPartPrimary::AvgPriDiameter(void) const {return m_avgpridiam;}

// Returns the surface area of the ith primary.
real PriPartPrimary::PriSurface(unsigned int i) const
{
    if (i < m_primaries.size()) {
        return m_primaries[i].Surface;
    } else {
        return 0.0;
    }
}

// Returns the total surface area of the primaries in the list.
real PriPartPrimary::PriSurface(void) const {return m_totprisurf;}

// Returns the volume of the ith primary.
real PriPartPrimary::PriVolume(unsigned int i) const
{
    return PriMass(i) / m_pmodel->Components(m_icomp)->Density();
}

// Returns the mass of the ith primary.
real PriPartPrimary::PriMass(unsigned int i) const
{
    if (i < m_primaries.size()) {
        return (m_primaries[i].Monomers * m_pmodel->Components(m_icomp)->MolWt()) / NA;
    } else {
        return 0.0;
    }
}


// BASIC DERIVED PARTICLE PROPERTIES.

// Calculates the derived properties from the unique properties.  This
// function is broadly similar to the version in the spherical Primary
// class except that it uses the surface-volume model.  Therefore the
// surface area is not altered and the collision diameter is calculated
// using the arithmetic mean function.
void PriPartPrimary::UpdateCache(void)
{
    // Use the surface-volume model to update basic particle
    // properties.
    SurfVolPrimary::UpdateCache();

    // Now update the pripart-list.
    updatePrimaries();
}



// OPERATIONS.

// Adjusts the primary with the given composition and 
// tracker values changes n times.  If the particle cannot be adjust
// n times, then this function returns the number of times
// it was adjusted.
unsigned int PriPartPrimary::Adjust(const fvector &dcomp, const fvector &dvalues, 
                                    unsigned int n)
{
    // Use the surface-volume model to adjust basic
    // primary properties.
    SurfVolPrimary::Adjust(dcomp, dvalues, n);

    // Distribute the added monomers among the pripart-list.
    distMonomers((unsigned int)dcomp[m_icomp]*n);

    // Ensure the primaries are up-to-date.
    updatePrimaries();

    return n;
}

// Combines this primary with another.  This is also the
// implementation of the + and += operators.
PriPartPrimary &PriPartPrimary::Coagulate(const Primary &rhs)
{
    // Use surface-volume model to coagulate.
    SurfVolPrimary::Coagulate(rhs);

    // Now deal with the primary particle list.
    if (rhs.AggID() == Spherical_ID) {
        // Add a new pri-particle to the list.
        vector<PriPart>::iterator i = m_primaries.insert(m_primaries.end(), PriPart());
        i->Monomers = (unsigned int)rhs.Composition(m_icomp);
        i->Surface  = calcSurf(i->Monomers);
    } else if (rhs.AggID() == SurfVol_ID) {
        // Add a number of new pri-particles to the
        // list until the surface area is matched.
        vector<PriPart>::iterator i = m_primaries.insert(m_primaries.end(), PriPart());
        i->Monomers = (unsigned int)rhs.Composition(m_icomp);
        i->Surface  = calcSurf(i->Monomers);
//        addRandom((unsigned int)rhs.Composition(m_icomp), rhs.SurfaceArea());
    } else {
        // Attempt to cast the RHS as a PriPartPrimary.  This
        // will throw an exception if the cast fails.
        const PriPartPrimary &prirhs = dynamic_cast<const PriPartPrimary&>(rhs);

        // Merge the pri-particle lists.
        mergeInList(prirhs.m_primaries);
    }

    return *this;
}

// This routine sinters the Primary for the given length of
// time using the provided sintering model.
void PriPartPrimary::Sinter(real dt,const Cell &sys, 
                            const Processes::SinteringModel &model)
{
    // Use surface-volume model to adjust the primary m_surf
    // property.
    SurfVolPrimary::Sinter(dt, sys, model);

    // Now update the primary particle list to match the new 
    // surface area.
    updatePrimaries();
}


// READ/WRITE/COPY.

// Returns a copy of the model data.
PriPartPrimary *const PriPartPrimary::Clone(void) const
{
    return new PriPartPrimary(*this);
}

// Writes the object to a binary stream.
void PriPartPrimary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output base class.
        SurfVolPrimary::Serialize(out);

        // Write the number of primary particles.
        unsigned int n = (unsigned int)m_primaries.size();
        out.write((char*)&n, sizeof(n));

        // Write the primary particle list.
        double val = 0.0;
        for (vector<PriPart>::const_iterator i=m_primaries.begin(); 
             i!=m_primaries.end(); ++i) {
            n   = i->Monomers;
            val = i->Surface;
            out.write((char*)&n, sizeof(n));
            out.write((char*)&val, sizeof(val));
        }

        // Write the total surface area.
        val = (double)m_totprisurf;
        out.write((char*)&val, sizeof(val));

        // Write the average particle diameter.
        val = (double)m_avgpridiam;
        out.write((char*)&val, sizeof(val));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PriPartPrimary::Serialize).");
    }
}

// Reads the object from a binary stream.
void PriPartPrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0, num = 0;
        double val     = 0.0;

        switch (version) {
            case 0:
                // Read base class.
                SurfVolPrimary::Deserialize(in, model);

                // Read number of primary particles.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read primary particle list.
                for (unsigned int i=0; i!=n; ++i) {
                    // Read monomer count.
                    in.read(reinterpret_cast<char*>(&num), sizeof(num));
                    // Read surface area.
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    // Add pri-part to the list.
                    vector<PriPart>::iterator j = m_primaries.insert(m_primaries.end(), 
                                                                     PriPart());
                    j->Monomers = num;
                    j->Surface  = val;
                }

                // Read the total surface area.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_totprisurf = (real)val;

                // Read the average primary diameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_avgpridiam = (real)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, PriPartPrimary::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PriPartPrimary::Deserialize).");
    }
}


// MEMORY MANAGEMENT.

// Release all memory associated with object.
void PriPartPrimary::releaseMem(void)
{
    m_primaries.clear();
}

// Initialisation routine.
void PriPartPrimary::init(void)
{
    releaseMem();
    m_totprisurf = 0.0;
    m_avgpridiam = 0.0;
}


// PRIPART-LIST MANAGEMENT.

// Calculates the surface area of a spherical primary with the
// given monomer count.
real PriPartPrimary::calcSurf(unsigned int n) const
{
    // The primary is assumed to be spherical with n monomers.
    // In order to calculate the surface area, the volume must
    // first be calculated.
    real vol = (real)n * m_pmodel->Components(m_icomp)->MolWt() / 
               (NA * m_pmodel->Components(m_icomp)->Density());

    // Return the surface area.
    return PI * pow(6.0*vol/PI, TWO_THIRDS);
}

// Calculates the average primary-particle diameter for
// the current list.  The primaries are assumed to be
// spherical.
real PriPartPrimary::calcAvgDiam()
{
    real sumd = 0.0;
    for (vector<PriPart>::const_iterator i=m_primaries.begin(); 
         i!=m_primaries.end(); ++i) {
        real vol = (real)i->Monomers  * m_pmodel->Components(m_icomp)->MolWt() / 
                   (NA * m_pmodel->Components(m_icomp)->Density());
        sumd += pow(vol*6.0/PI, ONE_THIRD);
    }
    m_avgpridiam = sumd / (real)(m_primaries.size());
    return m_avgpridiam;
}

// Distributes a number of monomers among the pri-particles, weighted
// by their surface areas.
void PriPartPrimary::distMonomers(unsigned int n)
{
    // Big N = number of primaries.
    unsigned int N = m_primaries.size();

    // dm = n = Number of monomer to distribute.
    unsigned int dm = n;

    // Sum up the surface areas in order to calculate probabilities.
    // Note:  The surface areas are summed in reverse (from smallest particle),
    //        and each monomer uses a different sum.  This ensures that the 
    //        correct total surface area is used when calculating the probabilities.
    static vector<real> surfsums;
    surfsums.resize(N, 0.0);
    surfsums[N-1] = m_primaries[N-1].Surface;
    for (int i=N-2; i>=0; --i) {
        surfsums[i] = surfsums[i+1] + m_primaries[i].Surface;
    }

    // Now distribute mass weighted by surface area, starting with the largest
    // particle.
    for (unsigned int i=0; i!=N; ++i) {
        // Determine number of monomers to be added to this primary (m).
        unsigned int m = ignbin(dm, (float)(m_primaries[i].Surface/surfsums[i]));
        // Add monomers.
        m_primaries[i].Monomers += m;
        // Calculate surface area.
        m_totprisurf -= m_primaries[i].Surface;
        m_primaries[i].Surface   = calcSurf(m_primaries[i].Monomers);
        m_totprisurf += m_primaries[i].Surface;
        // Update number of remaining monomers.
        if ((dm-=m) == 0) break;
    }

    // Ensure the primary list is sorted.
    sortList(0, m_primaries.size()-1);

    // Calculate the average diameter.
    calcAvgDiam();
}

// Removes n monomers from the primary particle list.  Monomers
// are chosen by weighting the particles by their surface area, 
// and starting selection from the largest primary.
void PriPartPrimary::removeMonomers(unsigned int n)
{
    // Big N = number of primaries.
    unsigned int N = m_primaries.size();

    // dm = n = Number of monomer to remove.
    unsigned int dm = n;

    // Sum up the surface areas in order to calculate probabilities.
    // Note:  The surface areas are summed in reverse (from smallest particle),
    //        and each monomer uses a different sum.  This ensures that the 
    //        correct total surface area is used when calculating the probabilities.
    static vector<real> surfsums;
    surfsums.resize(N, 0.0);
    surfsums[N-1] = m_primaries[N-1].Surface;
    for (unsigned int i=N-2; i>=0; --i) {
        surfsums[i] = surfsums[i+1] + m_primaries[i].Surface;
    }

    // Now remove mass weighted by surface area, starting with the largest
    // particle.
    for (unsigned int i=0; i!=N; ++i) {
        // Determine number of monomers to be added to this primary (m).
        unsigned int m = ignbin(dm, (float)(m_primaries[i].Surface/surfsums[i]));
        // Remove monomers.
        m_primaries[i].Monomers -= m;
        // Calculate surface area.
        m_totprisurf -= m_primaries[i].Surface;
        m_primaries[i].Surface   = calcSurf(m_primaries[i].Monomers);
        m_totprisurf += m_primaries[i].Surface;
        // Update number of remaining monomers to remove.
        if ((dm-=m) == 0) break;
    }

    // Ensure the primary list is sorted.
    sortList(0, m_primaries.size()-1);

    // Calculate the average primary diameter.
    calcAvgDiam();
}

// Updates the primary particle list to attempt to match the current
// primary surface area (given by surface-volume model).
void PriPartPrimary::updatePrimaries(void)
{
    // Check that the primaries have got the correct number
    // of monomers before adjusting.
    unsigned int mass = 0;
    for (vector<PriPart>::iterator i=m_primaries.begin(); i!=m_primaries.end(); ++i) {
        mass += i->Monomers;
    }
    if (mass < m_comp[m_icomp]) {
        if (m_primaries.size() > 0) {
            // If primaries exist already then distribute the added
            // mass.
            distMonomers((unsigned int)m_comp[m_icomp] - mass);
        } else {
            // Add the first primary.
            m_primaries.push_back(PriPart());
            m_primaries[0].Monomers = (unsigned int)m_comp[m_icomp];
            m_primaries[0].Surface = calcSurf((unsigned int)m_comp[m_icomp]);
            m_totprisurf = m_primaries[0].Surface;
            m_avgpridiam = calcAvgDiam();
        }
    } else if (mass > m_comp[m_icomp]) {
        // Can't remove monomers if there are no primaries.  This
        // should never happen anyway.
        if (m_primaries.size() > 0) {
            removeMonomers(mass - (unsigned int)m_comp[m_icomp]);
        }
    }

    // Calculate required change in surface area to match
    // surface-volume model.
    real ds = m_surf - m_totprisurf;

    // Now the difficult bit, probably need to destroy primaries to 
    // reduce surface area.  The methodology is to remove the smallest
    // primaries, adding their mass to the larger, until the total surface
    // area is less than or equal to the aggregate surface area.
    if (ds < 0.0) {
        while ((m_totprisurf > m_surf) && (m_primaries.size() > 1)) {
            // Store the mass of the smallest primary, then remove it.
            vector<PriPart>::iterator i = (m_primaries.end()-1);
            unsigned int m = i->Monomers; // Store monomer count.
            m_totprisurf -= i->Surface;   // Remove contribution to total surface area.
            m_primaries.erase(i);         // Erase the primary.

            // Distribute this mass over the remaining primaries.
            distMonomers(m);
        }
    }

    // Ensure that the list is sorted.
    sortList(0, m_primaries.size()-1);

    // Calculate the average primary diameter.
    calcAvgDiam();
}

// Merges a pri-particle list into the current list.
void PriPartPrimary::mergeInList(const std::vector<PriPart> &list)
{
    vector<PriPart>::const_iterator j;
    unsigned int i;

    // Loop over list to be merged.
    for (j=list.begin(); j!=list.end(); ++j) {
        // Loop over this list and find the first primary smaller than 
        // that being inserted.
        for (i=0; i!=m_primaries.size(); ++i) {
            if (j->Monomers > m_primaries[i].Monomers) {
                m_primaries.insert(m_primaries.begin()+i, PriPart(*j));
                ++i;
                break;
            }
        }
        // If the primary being inserted is smaller than all primaries in
        // the list then add it at the end.
        if (i==m_primaries.size()) {
            m_primaries.insert(m_primaries.end(), PriPart(*j));
        }
    }

    // Calculate the average primary diameter.
    calcAvgDiam();
}

// Sorts the primary particle list in descending order of mass.
void PriPartPrimary::sortList(unsigned int i1, unsigned int i2)
{
    if (m_primaries.size() <= 1) {
        // No need to sort if the array only has one element.
        return;
    } else if (i2-i1 < 15) {
        // Too few elements for a quick-sort to be efficient.  Use a 
        // bubble sort instead.

        // Temporary storage variable.
        PriPart temp;

        // Perform bubble sort.
        for (unsigned int j=i1; j<i2; ++j) {
            for (unsigned int k=j+1; k<=i2; ++k) {
                if (m_primaries[j].Monomers < m_primaries[k].Monomers) {
                    // Store j in temp.
                    temp.Monomers = m_primaries[j].Monomers;
                    temp.Surface  = m_primaries[j].Surface;
                    // Put k into j.
                    m_primaries[j].Monomers = m_primaries[k].Monomers;
                    m_primaries[j].Surface  = m_primaries[k].Surface;
                    // Put temp into k.
                    m_primaries[k].Monomers = temp.Monomers;
                    m_primaries[k].Surface  = temp.Surface;
                }
            }
        }
    } else {
        // Use a quick-sort algorithm.

        // Choose a pivot at the mid-point.
        unsigned int ip = i1 + ((i2 - i1)/2);

        // Temporary storage variable.
        PriPart temp;

        // Iteration indices.
        unsigned int j=i1, k=i2;

        // Sort everything larger than the pivot to the front of 
        // the list, and everthing smaller to the end of the list.
        while (j < k) {
            // Find next element below pivot which is smaller, and next
            // element above pivot which is larger.
            while(m_primaries[ip].Monomers <= m_primaries[j].Monomers) {
                ++j;
                if (j >= k) break;
            }
            while(m_primaries[ip].Monomers > m_primaries[k].Monomers) {
                if (k == 0) break;
                --k;
            }

            // Swap j and k.
            if (j < k) {
                // Store j in temp.
                temp.Monomers = m_primaries[j].Monomers;
                temp.Surface  = m_primaries[j].Surface;
                // Put k into j.
                m_primaries[j].Monomers = m_primaries[k].Monomers;
                m_primaries[j].Surface  = m_primaries[k].Surface;
                // Put temp into k.
                m_primaries[k].Monomers = temp.Monomers;
                m_primaries[k].Surface  = temp.Surface;
            } else {
                break;
            }
        }

        // Now sort the two halves of the list either side of the
        // pivot.
        if (i1 < ip-1) sortList(i1, ip-1);
        if (ip < i2) sortList(ip, i2);
    }   
}
