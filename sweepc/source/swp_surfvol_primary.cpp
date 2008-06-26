#include "swp_primary.h"
#include "swp_surfvol_primary.h"
#include "swp_submodel_type.h"
#include "swp_aggmodel_type.h"
#include "swp_model_factory.h"
#include "swp_surfvol_cache.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace Sweep::SubModels;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
SurfVolPrimary::SurfVolPrimary(void)
: m_sphsurf(0.0)
{
}

// Initialising constructor.
SurfVolPrimary::SurfVolPrimary(real time, const Sweep::ParticleModel &model)
: Primary(time, model), m_sphsurf(0.0)
{
}

// Copy constructor.
SurfVolPrimary::SurfVolPrimary(const SurfVolPrimary &copy) 
{
    *this = copy;
}

// Stream-reading constructor.
SurfVolPrimary::SurfVolPrimary(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
SurfVolPrimary::~SurfVolPrimary()
{
    releaseMem();
}


// OPERATOR OVERLOADS.

// Assignment operator (Primary RHS).
SurfVolPrimary &SurfVolPrimary::operator=(const Primary &rhs)
{
    // Now check the type of the RHS to see if additional data
    // copying is required.
    if (rhs.AggID() == Spherical_ID) {
        Primary::operator=(rhs);
        // The RHS is a spherical primary, therefore the 
        // surface area is equal to the equiv. sphere
        // surface.  We must assign this separately because the
        // copyTo() implementation for a spherical primary is not
        // aware of this variable.
        m_sphsurf = rhs.SurfaceArea();
    } else {
        // Attempt to cast the primary to a SurfVolPrimary.  This will
        // throw an exception if the cast fails.
        operator=(dynamic_cast<const SurfVolPrimary&>(rhs));
    }

    return *this;
}

// Assignment operator (SurfVolPrimary RHS).
SurfVolPrimary &SurfVolPrimary::operator=(const SurfVolPrimary &rhs)
{
    // First copy everything for spherical primaries.
    Primary::operator=(rhs);

    // Now set surface-volume specific variables.
    m_sphsurf = rhs.m_sphsurf;

    return *this;
}

// Compound assignment (coagulation).
SurfVolPrimary &SurfVolPrimary::operator+=(const Primary &rhs)
{
    return SurfVolPrimary::Coagulate(rhs);
}


// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType SurfVolPrimary::AggID(void) const {return AggModels::SurfVol_ID;}

// Creates an aggregation data cache for this primary type.
AggModels::SurfVolCache *const SurfVolPrimary::CreateAggCache(ParticleCache &pcache) const
{
    SurfVolCache *cache = static_cast<SurfVolCache*>(ModelFactory::CreateAggCache(AggModels::SurfVol_ID, pcache));
    if (cache != NULL) *cache = *this;
    return cache;
}


// BASIC DERIVED PARTICLE PROPERTIES.

// Calculates the derived properties from the unique properties.  This
// function is broadly similar to the version in the spherical Primary
// class except that it uses the surface-volume model.  Therefore the
// surface area is not altered and the collision diameter is calculated
// using the arithmetic mean function.
void SurfVolPrimary::UpdateCache(void)
{
    // Store the correct surface area.
    real s = m_surf;

    // Pretend that the primary is spherical and set the cache
    // accordingly.
    Primary::UpdateCache();

    // The surface area is now set to that of a sphere. Set 
    // correct surface area (minimum is spherical surface area).
    m_sphsurf = m_surf;
    m_surf = max(s, m_sphsurf);

   
    // Calculate diameters.
    m_dcol = (m_diam + sqrt(m_surf / PI)) * 0.5;
    m_dmob = m_dcol; // TODO:  Correct expression for Dmob.
}

// Returns the equivalent spherical particle surface area.
real SurfVolPrimary::SphSurfaceArea(void) const {return m_sphsurf;}

// Returns the number of primary particles if the aggregate is assumed
// to consist of mono-sized primaries.
unsigned int SurfVolPrimary::PP_Count(void) const
{
    // Note the minimum number of primary particles must be 1.
    return max(1u, (unsigned int)((m_surf * m_surf * m_surf) / 
                                  (36.0 * PI * m_vol * m_vol)));
}

// Returns the primary particle diameter if the aggregate is assumed
// to consist of mono-sized primaries.
real SurfVolPrimary::PP_Diameter(void) const
{
    // This should always be <= equiv. sphere diameter.
    return 6.0 * m_vol / m_surf;
}


// OPERATIONS.

// Adjusts the primary with the given composition and 
// tracker values changes n times.  If the particle cannot be adjust
// n times, then this function returns the number of times
// it was adjusted.
unsigned int SurfVolPrimary::Adjust(const fvector &dcomp, const fvector &dvalues, 
                                    unsigned int n)
{
    // Calculate change in volume.
    real dvol = 0.0;
    for (unsigned int i=0; i!=dcomp.size(); ++i) {
        dvol += dcomp[i] * m_pmodel->Components(i)->MolWt() / 
                m_pmodel->Components(i)->Density();
    }
    dvol *= (real)n / NA;

    // Calculate change in surface area.
    real rad = 0.0;
    if (dvol > 0.0) {
        // Inverse growth radius.
        rad = sqrt(4.0 * PI / m_surf);
    } else {
        // Inverse oxidation radius.    
        rad = m_surf / (3.0 * m_vol);
    }

    // Save new surface area.
    real s = m_surf + (2.0 * dvol * rad);

    // Adjust the particle assuming that it is spherical.
    Primary::Adjust(dcomp, dvalues, n);

    // Set correct surface area, which was incorrectly set by
    // Primary::Adjust.
    m_sphsurf = m_surf;
    m_surf    = max(s, m_sphsurf);

    // This has a knock-on affect of changing the collision diameter.
    // Note, we can avoid recalling UpdateCache() here, because only
    // a couple of cached values will have changed.
    m_dcol = 0.5 * (m_dcol * sqrt(m_surf));
    m_dmob = m_dcol;

    return n;
}

// Combines this primary with another.  This is also the
// implementation of the + and += operators.
SurfVolPrimary &SurfVolPrimary::Coagulate(const Primary &rhs)
{
    // Store the resultant surface area.
    real s = m_surf + rhs.SurfaceArea();

    // Perform the coagulation.
    Primary::Coagulate(rhs);

    // The spherical particle Coagulate() routine has set the
    // surface area incorrectly.  We now replace the surface area
    // to the correct point-contact value.
    m_sphsurf = m_surf;
    m_surf    = max(m_sphsurf, s);

    // This has a knock-on affect of changing the collision diameter.
    // Note, we can avoid recalling UpdateCache() here, because only
    // a couple of cached values will have changed.
    m_dcol = 0.5 * (m_dcol * sqrt(m_surf));
    m_dmob = m_dcol;

    return *this;
}


// READ/WRITE/COPY.

// Returns a copy of the model data.
SurfVolPrimary *const SurfVolPrimary::Clone(void) const
{
    return new SurfVolPrimary(*this);
}

// Returns this object's instance.  This may seem rather circular, but
// it has an important purpose for getting the correct object type reference
// from a base class reference/pointer.
SurfVolPrimary &SurfVolPrimary::Instance() {return *this;}
const SurfVolPrimary &SurfVolPrimary::Instance() const {return *this;}

// Writes the object to a binary stream.
void SurfVolPrimary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output base class.
        Primary::Serialize(out);

        // Write spherical surface area.
        double val = (double)m_sphsurf;
        out.write((char*)&val, sizeof(val));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, SurfVolPrimary::Serialize).");
    }
}

// Reads the object from a binary stream.
void SurfVolPrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;

        switch (version) {
            case 0:
                // Read base class.
                Primary::Deserialize(in, model);

                // Read spherical surface area.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_sphsurf = (real)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, SurfVolPrimary::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SurfVolPrimary::Deserialize).");
    }
}


// Returns this object's instance.  This may seem rather circular, but
// it has an important purpose for getting the correct object type reference
// from a base class reference/pointer.
//const SurfVolPrimary &SurfVolPrimary::TypedRef() const {return *this;}

/*
// Copies this surface-volume primary into a spherical
// primary.  This implements the assignment operator for
// the situation Primary = SurfVolPrimary.  This results
// in a loss of information.
Primary &SurfVolPrimary::CopyTo(Primary &lhs) const
{
    // Copy spherical primary data.
    Primary::CopyTo(lhs);
    // Set the correct spherical surface area.
    lhs.SetSurfaceArea(m_sphsurf);
    return lhs;
}

// Copies this surface-volume primary into another surface-volume
// primary.  This implements the assignment operator.
SurfVolPrimary &SurfVolPrimary::CopyTo(SurfVolPrimary &lhs) const
{
    // Copy spherical primary data.
    Primary::CopyTo(lhs);
    // Set the correct spherical surface area.
    lhs.m_sphsurf = m_sphsurf;
    return lhs;
}

// Adds this surface-volume primary to a spherical 
// primary.  This implements coagulation up inheritance chain,
// and results in a loss of data (resultant particle is spherical):
//     Primary + SurfVolPrimary = Primary.
Primary &SurfVolPrimary::AddTo(Primary &lhs) const
{
    // Perform spherical primary coagulation.
    return Primary::AddTo(lhs);
}

// Adds this surface-volume primary to another surface-volume 
// primary.  This implements coagulation:
//     SurfVolPrimary + SurfVolPrimary = SurfVolPrimary.
SurfVolPrimary &SurfVolPrimary::AddTo(SurfVolPrimary &lhs) const
{
    // Save the LHS surface area.
    real s = lhs.m_surf;

    // Perform spherical primary coagulation.
    Primary::AddTo(lhs);

    // Store the equiv. sphere surface area.
    lhs.m_sphsurf = lhs.m_surf;

    // Set the surface area to be the sum of the areas before
    // coagulation occured (point contact).
    lhs.m_surf = m_surf + s;

    return lhs;
}
*/
