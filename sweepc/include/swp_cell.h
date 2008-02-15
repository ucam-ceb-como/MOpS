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
#include "swp_ensemblestats.h"
#include "sprog.h"
#include <string>
#include <iostream>

namespace Sweep
{
class Mechanism;

class Cell : public Sprog::Thermo::IdealGas
{
public:
    // Constructors.
    Cell(const Sprog::SpeciesPtrVector &sp); // Default constructor.
    Cell(const Cell &copy);                  // Copy constructor.
    Cell(                                    // Stream-reading constructor.
        std::istream &in,                    //   - Stream from which to read.
        const Mechanism &mech                //   - Mechanism used to define particles.
        );

    // Destructor.
    virtual ~Cell(void);

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

    // Returns particle statistics.
    void GetVitalStats(EnsembleStats &stats) const;


    // SCALING ROUTINES INCL. SAMPLE VOLUME.

    // Returns the real system to stochastic system scaling factor.
    real SampleVolume() const;

    // Sets the number density which the ensemble represents.
    int SetM0(real m0);

    // Sets the number density which the full 
    // ensemble would represent.
    int SetMaxM0(real m0);

    void Reset(real m0);


    // READ/WRITE/COPY.

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,     // Input stream.
        const Mechanism &mech //   - Mechanism used to define particles.
        );

protected:
    // Default constructor is protected as it makes no
    // sense to define a mixture without knowledge of the
    // definin species.  This trait is brought over from Sprog.
    Cell(void);

private:
    // Particle ensemble.
    Sweep::Ensemble m_ensemble;

    // The volume in which the ensemble represents
    // the complete real system.
    real m_smpvol;
};
};

#endif
