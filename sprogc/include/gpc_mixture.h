/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    The Mixture class is designed as a base class for all mixtures.  It maintains a group
    of species mole fractions, the temperature and the mixture density.
*/

#ifndef GPC_MIXTURE_H
#define GPC_MIXTURE_H

#include "gpc_params.h"
#include "gpc_species.h"
#include "gpc_unit_systems.h"
#include <vector>

namespace Sprog
{
namespace Thermo
{
class Mixture
{
public:
    // Constructors.
    Mixture(void); // Default constructor.

    // Destructors.
    virtual ~Mixture(void); // Default destructor.

    // Operator overloads.
    Mixture &operator=(const Mixture &mix); // Assignment operator.

    // Creates a copy of the mixture object.
    virtual Mixture* Clone() const;

    // Get/Set temperature.
    real T(void) const; // Returns temperature.
    void SetT(real T);  // Set the temperature.

    // Get/Set species concentrations/fractions.
    const std::vector<real> & MoleFractions() const;       // Returns the mole fractions of all species.
    void GetConcs(std::vector<real> &concs) const;         // Returns the molar concentrations of all species in current units.
    void GetMassFractions(std::vector<real> &fracs) const; // Returns the mass fractions of all species.
    real MoleFraction(unsigned int i) const; // Returns the mole fraction of species i.
    real MolarConc(unsigned int i) const;    // Returns the molar concentration of species i.
    real MassFraction(unsigned int i) const; // Returns the mass fraction of species i.
    void SetFracs(const std::vector<real> &fracs);         // Sets the mole fractions of all species.
    void SetConcs(const std::vector<real> &concs);         // Sets the molar concentrations of all species.
    void SetMassFracs(const std::vector<real> &fracs);     // Sets the mass fractions of all species.

    // Get/Set density.
    real Density(void) const;       // Returns molar density in current units.
    real MassDensity(void) const;   // Returns mass density in current units.
    void SetDensity(real dens);     // Sets the molar density.
    void SetMassDensity(real dens); // Sets the mass density.

    // Get/Set the species which define the mixture.
    const SpeciesPtrVector *const Species() const; // Returns the species for which this mixture is defined.
    void SetSpecies(const SpeciesPtrVector *const sp); // Sets the species for which this mixture is defined.

protected:
    // Mixture properties.
    std::vector<real> m_frac; // Species mole fractions.
    real m_T;    // Mixture temperature.
    real m_dens; // Mixture density.

    // Mixture context.
    const SpeciesPtrVector *m_species; // Vector of species for which this mixture is defined.
};
};
};

#endif