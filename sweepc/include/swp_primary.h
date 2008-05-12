/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The Primary class defines a spherical primary particle 
    used by the PriPartModel.
*/

#ifndef SWEEP_PRIMARY_H
#define SWEEP_PRIMARY_H

#include "swp_params.h"
#include <iostream>

namespace Sweep
{
class Primary
{
public:
    // Constructors.
    Primary(void);                // Default constructor.
    Primary(const Primary &copy); // Copy constructor.
    Primary(std::istream &in);    // Stream-reading constructor.

    // Destructors.
    ~Primary(void);

    // Operators.
    Primary &operator=(const Primary &rhs);
    Primary &operator+=(const Primary &rhs);
    const Primary operator+(const Primary &rhs) const;
    bool operator==(const Primary &rhs) const;
    bool operator!=(const Primary &rhs) const;
    bool operator<(const Primary &rhs) const;
    bool operator>(const Primary &rhs) const;
    bool operator<=(const Primary &rhs) const;
    bool operator>=(const Primary &rhs) const;

    // PRIMARY COMPOSITION.

    // Returns the composition of the primary (only one component).
    unsigned int Composition(void) const;

    // Sets the primary composition.
    void SetComposition(unsigned int comp);

    // Changes the composition by the given amount.
    void ChangeComposition(int dc);


    // DERIVED PROPERTIES.

    // Returns the primary volume.
    real Volume(void) const;

    // Returns the primary mass.
    real Mass(void) const;

    // Returns the primary diameter.
    real Diameter(void) const;

    // Returns the primary surface area.
    real SurfaceArea(void) const;


    // READ/WRITE/COPY.

    // Returns a copy of the model data.
    Primary *const Clone(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(std::istream &in);

private:
    // Primary composition (only one component).
    unsigned int m_comp;

    // Properties calculated from the composition.
    real m_vol;  // Primary volume.
    real m_mass; // Primary mass.
    real m_diam; // Diameter.
    real m_surf; // Surface area.
};
};

#endif
