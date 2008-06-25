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
#include "swp_particle_model.h"
#include "swp_ensemble.h"
#include "swp_ensemble_stats.h"
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
    Cell(const Sweep::ParticleModel &model); // Default constructor.
    Cell(const Cell &copy);                  // Copy constructor.
    Cell(                                 // Stream-reading constructor.
        std::istream &in,                 //   - Stream from which to read.
        const Sweep::ParticleModel &model //   - Model used to define particles.
        );

    // Destructor.
    virtual ~Cell(void);

    // Operators.
    Cell &operator=(const Cell &rhs);
    Cell &operator=(const Sprog::Thermo::IdealGas &rhs);

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
    Sweep::Ensemble &Particles(void);
    const Sweep::Ensemble &Particles(void) const;

    // Returns the number of particles in the ensemble.
    unsigned int ParticleCount(void) const;

    // Returns particle statistics.
    void GetVitalStats(Stats::EnsembleStats &stats) const;


    // THE PARTICLE MODEL.

    // Returns the particle model used to define particles in this
    // cell.
    const Sweep::ParticleModel *const ParticleModel(void) const;


    // SCALING ROUTINES INCL. SAMPLE VOLUME.

    // Returns the real system to stochastic system scaling factor.
    real SampleVolume() const;

    // Sets the number density which the ensemble represents.
    int SetM0(real m0);

    // Sets the number density which the full 
    // ensemble would represent.
    int SetMaxM0(real m0);

    void Reset(real m0);


    // FIXED/VARIABLE CHEMISTRY.

    // Returns whether or not the chemical conditions are fixed.
    bool FixedChem() const;

    // Sets whether or not the chemical conditions are fixed.
    void SetFixedChem(bool fixed = true);

    // Set the chemical conditions to be variable.
    void SetVariableChem(bool vari = true);


    // READ/WRITE/COPY.

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,                 // Input stream.
        const Sweep::ParticleModel &model // Model used to define particles.
        );

protected:
    // Default constructor is protected as it makes no
    // sense to define a mixture without knowledge of the
    // definin species.  This trait is brought over from Sprog.
    Cell(void);

private:
    // Particle ensemble.
    Sweep::Ensemble m_ensemble;

    // Particle model.
    const Sweep::ParticleModel *m_model;

    // The volume in which the ensemble represents
    // the complete real system.
    real m_smpvol;

    // Flag determining whether or not the chemistry in this system is fixed.
    // If the chemical conditions are fixed, then they cannot be altered by
    // any particle processes.  Default is false.
    bool m_fixed_chem;
};
};

#endif
