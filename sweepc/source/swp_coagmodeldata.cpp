#include "swp_coagmodeldata.h"
#include "swp_coagmodel.h"

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
CoagModelData::CoagModelData(void)
{
    m_dcolsqr      = 0.0;
    m_inv_dcol     = 0.0;
    m_inv_dcolsqr  = 0.0;
    m_inv_sqrtmass = 0.0;
    m_d2_m_1_2     = 0.0;
}

// Default constructor (public).
CoagModelData::CoagModelData(Sweep::ParticleData &parent)
: IModelData(parent)
{
}

// Copy constructor.
CoagModelData::CoagModelData(const Sweep::CoagModelData &copy)
{
    // Use assignment operator.
    *this = copy;
}

// Default destructor.
CoagModelData::~CoagModelData()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
CoagModelData &CoagModelData::operator=(const Sweep::CoagModelData &rhs)
{
    if (this != &rhs) {
        IModelData::operator=(rhs);
        m_dcolsqr      = rhs.m_dcolsqr;
        m_inv_dcol     = rhs.m_inv_dcol;
        m_inv_dcolsqr  = rhs.m_inv_dcolsqr;
        m_inv_sqrtmass = rhs.m_inv_sqrtmass;
        m_d2_m_1_2     = rhs.m_d2_m_1_2;
    }
    return *this;
}

// Compound assignment operator.
CoagModelData &CoagModelData::operator+=(const Sweep::CoagModelData &rhs)
{
    IModelData::operator+=(rhs);
    m_dcolsqr      += rhs.m_dcolsqr;
    m_inv_dcol     += rhs.m_inv_dcol;
    m_inv_dcolsqr  += rhs.m_inv_dcolsqr;
    m_inv_sqrtmass += rhs.m_inv_sqrtmass;
    m_d2_m_1_2     += rhs.m_d2_m_1_2;
    return *this;
}

// Addition operator.
const CoagModelData CoagModelData::operator +(const Sweep::CoagModelData &rhs) const
{
    return CoagModelData(*this) += rhs;
}


// COAGULATION MODEL PARTICLE PROPERTIES.

// Collision diameter squared (cm2).
real CoagModelData::CollDiamSquared() const {return m_dcolsqr;}

// Inverse collision diameter (cm-1).
real CoagModelData::InvCollDiam() const {return m_inv_dcol;}

// Inverse squared collision diameter (cm-2).
real CoagModelData::InvCollDiamSquared() const {return m_inv_dcolsqr;}

// Inverse of square root of mass (g-1/2).
real CoagModelData::InvSqrtMass() const {return m_inv_sqrtmass;}

// Collision diameter squared times the inverse square root of mass.
real CoagModelData::CollDiamSqrdInvSqrtMass() const {return m_d2_m_1_2;}

// Returns the property with the given ID.
real CoagModelData::Property(Sweep::CoagModelData::PropertyID id) const
{
    switch (id) {
        case iD2:
            return m_dcolsqr;
        case iD_1:
            return m_inv_dcol;
        case iD_2:
            return m_inv_dcolsqr;
        case iM_1_2:
            return m_inv_sqrtmass;
        case iD2_M_1_2:
            return m_d2_m_1_2;
        default:
            return 0.0;
    }
}

// Returns the property with the given ID.
real CoagModelData::Property(unsigned int i) const
{
    return Property((PropertyID)i);
}


// MODEL.

const CoagModel &CoagModelData::Model() const
{
    return CoagModel::Instance();
}


// READ/WRITE/COPY.

// Returns a clone of the coag model data object.
CoagModelData *const CoagModelData::Clone(void) const
{
    return new CoagModelData(*this);
}

// Returns the model data type. 
ModelType CoagModelData::ID(void) const {return CoagModel_ID;}
