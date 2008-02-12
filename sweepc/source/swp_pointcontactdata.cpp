#include "swp_pointcontactdata.h"

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
PointContactData::PointContactData(void)
: m_sphsurf(0.0), m_surf(0.0)
{
}

// Default constructor (public).
PointContactData::PointContactData(Sweep::ParticleData &parent)
: CoagModelData(parent)
{
}

// Copy constructor.
PointContactData::PointContactData(const Sweep::PointContactData &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Default destructor.
PointContactData::~PointContactData()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator.
PointContactData &PointContactData::operator=(const Sweep::PointContactData &rhs)
{
    if (this != &rhs) {
        CoagModelData::operator=(rhs);
        m_sphsurf = rhs.m_sphsurf;
        m_surf    = rhs.m_surf;
    }
    return *this;
}

// Compound assignment.
PointContactData &PointContactData::operator +=(const Sweep::PointContactData &rhs)
{
    CoagModelData::operator +=(rhs);
    m_sphsurf += rhs.m_sphsurf;
    m_surf    += rhs.m_surf;
    return *this;
}

// Addition operator.
const PointContactData PointContactData::operator +(const Sweep::PointContactData &rhs) const
{
    return PointContactData(*this) += rhs;
}


// COAGULATION MODEL PARTICLE PROPERTIES.

// Returns the equivalent spherical surface area.
real PointContactData::SphSurfaceArea(void) const {return m_sphsurf;}

// Returns the actual surface area.
real PointContactData::SurfaceArea(void) const {return m_surf;}


// MODEL WHICH USES THIS DATA.

// Returns the PointContactModel which operates on this data.
const PointContactModel &PointContactData::Model(void) const
{
    return PointContactModel::Instance();
}


// READ/WRITE/COPY.

// Returns a copy of the data.
PointContactData *const PointContactData::Clone(void) const
{
    return new PointContactData(*this);
}

// Returns the model ID.
ModelType PointContactData::ID(void) const
{
    return SVModel_ID;
}
