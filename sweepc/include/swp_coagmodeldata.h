/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The CoagModelData class is a specialisation of the ModelData
    class which defines the additional particle properties required
    to calculate coagulation rates.
*/

#ifndef SWEEP_COAGMODELDATA_H
#define SWEEP_COAGMODELDATA_H

#include "swp_params.h"
#include "swp_modeldata.h"
#include "swp_coagmodel.h"
#include <iostream>

namespace Sweep
{
class CoagModelData : public IModelData
{
friend class CoagModel;

public:
    enum PropertyID {iD2, iD_1, iD_2, iM_1_2, iD2_M_1_2};

    // Constructors.
    CoagModelData(ParticleData &parent);      // Default constructor.
    CoagModelData(const CoagModelData &copy); // Copy constructor.
    CoagModelData(            // Stream-reading constructor.
        std::istream &in,     //   - Input stream.
        ParticleData &parent  //   - Parent ParticleData object.
        );

    // Destructors.
    virtual ~CoagModelData(void);

    // Operators.
    virtual CoagModelData &operator=(const CoagModelData &rhs);
    virtual CoagModelData &operator+=(const CoagModelData &rhs);
    virtual const CoagModelData operator+(const CoagModelData &rhs) const;


    // COAGULATION MODEL PARTICLE PROPERTIES.

    // Collision diameter squared (cm2).
    real CollDiamSquared() const;

    // Inverse collision diameter (cm-1).
    real InvCollDiam() const;

    // Inverse squared collision diameter (cm-2).
    real InvCollDiamSquared() const;

    // Inverse of square root of mass (g-1/2).
    real InvSqrtMass() const;

    // Collision diameter squared times the inverse square root of mass.
    real CollDiamSqrdInvSqrtMass() const;

    // Returns the property with the given ID.
    real Property(PropertyID id) const;
    real Property(unsigned int id) const;

    // MODEL WHICH USES THIS DATA.

    virtual const CoagModel &Model(void) const;


    // READ/WRITE/COPY.

    // Returns a copy of the coagulation model data.
    virtual CoagModelData *const Clone(void) const;

    // Returns the model data type. 
    virtual ModelType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(std::istream &in);

protected:
    // Can't create a CoagModelData without knowledge
    // of the parent ParticleData.
    CoagModelData(void);

private:
    // Coagulation model particle properties.
    real m_dcolsqr;      // Collision diameter squared.
    real m_inv_dcol;     // Inverse collision diameter.
    real m_inv_dcolsqr;  // Inverse of the diameter squared.
    real m_inv_sqrtmass; // Inverse of the square-root of the mass.
    real m_d2_m_1_2;     // D^2 * M^-1/2.

};
};

#endif
