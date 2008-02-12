/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    A specialised ideal gas class for Sweep, which contains an ensemble
    with particles of type Particle.  The Cell class inherits the
    Sprog::IdealGas class and adds an additional layer of data to it
    to describe the particles within the cell.
*/

#ifndef SWEEP_CELL_H
#define SWEEP_CELL_H

#include "swp_params.h"
#include "swp_ensemble.h"
#include "sprog.h"
#include <string>
#include <iostream>

namespace Sweep
{
class Cell : public Sprog::Thermo::IdealGas
{
public:
    // Constructors.
    Cell(const Sprog::SpeciesPtrVector &sp); // Default constructor.
    Cell(const Cell &copy);                  // Copy constructor.
    Cell(                                    // Stream-reading constructor.
        std::istream &in,                    //   - Stream from which to read.
        const Sprog::SpeciesPtrVector &sp    //   - Definition of gas-phase species.
        );

    // Destructor.
    ~Cell(void);

    // Operators.
    Cell &operator=(const Cell &rhs);


    // THE GAS-PHASE INTERFACE.

    // Returns the description of the gas-phase mixture.
    const Sprog::Thermo::IdealGas &GasPhase(void) const;

    // Sets the gas-phase mixture.
    void SetGasPhase(const Sprog::Thermo::IdealGas &gas);

    // Adjusts the concentration of the ith species.
    void AdjustConc(unsigned int i, real dc);

    // Adjusts the concentration of all species.
    void AdjustConcs(const fvector &dc);


    // THE PARTICLE ENSEMBLE.

    // Returns the particle ensemble.
    Ensemble &Particles(void);
    const Ensemble &Particles(void) const;

    // Returns the number of particles in the ensemble.
    unsigned int ParticleCount(void) const;


    // SCALING ROUTINES INCL. SAMPLE VOLUME.

    // Returns the real system to stochastic system scaling factor.
    real SampleVolume() const;

    // Sets the number density which the ensemble represents.
    int SetM0(real m0);

    // Sets the number density which the full 
    // ensemble would represent.
    int SetMaxM0(real m0);

    void Reset(real m0);

private:
    // Particle ensemble.
    Sweep::Ensemble m_ensemble;

    // The volume in which the ensemble represents
    // the complete real system.
    real m_smpvol;

    // Default constructor is protected as it makes no
    // sense to define a mixture without knowledge of the
    // definin species.  This trait is brought over from Sprog.
    Cell(void);
};
};

#endif
