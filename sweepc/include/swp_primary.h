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

namespace Sweep
{
class Primary
{
public:
    // Constructors.
    Primary(void);   // Default constructor.
    Primary(const Primary &copy); // Copy constructor.

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

    // PROPERTIES.

    // Returns the primary diameter.
    real Diameter(void) const;

    // Returns the primary volume.
    real Volume(void) const;

    // Sets the primary volume.
    void SetVolume(real vol);

    // Returns the primary mass.
    real Mass(void) const;

    // Sets the primary mass.
    void SetMass(real mass);

    // Adds/removes some mass to/from the primary.
    void ChangeMass(real dm);

    // Returns the primary surface area.
    real SurfaceArea(void) const;


    // READ/WRITE/COPY.

    // Returns a copy of the model data.
    Primary *const Clone(void) const;

private:
    real m_vol;  // Primary volume.
    real m_mass; // Primary mass.

    // Properties calculated from the volume.
    real m_diam; // Diameter.
    real m_surf; // Surface area.
};
};

#endif
