/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    A specialised system for Sweep, which contains an ensemble
    with particles of type Particle.
*/

#ifndef SWEEP_SYSTEM_H
#define SWEEP_SYSTEM_H

#include "swpparams.h"
#include "swpensemble.h"
#include "gpcspecieslist.h"
#include <string>

namespace Sweep
{
class System
{
protected:
    Sweep::Ensemble m_ensemble;
    real m_smpvol; // The volume in which this represents the real system.
    SpeciesList *m_species;
public:
    System(void);
    ~System(void);
public: // Ensemble and particle related routines.
    Sweep::Ensemble &Ensemble(void);
    const Sweep::Ensemble &ConstEnsemble(void) const;
    /* Returns the number of particles in the ensemble. */
    unsigned int ParticleCount(void) const;
public: // Chemical environment get/sets.
    virtual void SetSpeciesList(SpeciesList &list);
    SpeciesList &GetSpeciesList(void);
    virtual real GetTemperature(const real t) const = 0;
    virtual real GetPressure(const real t) const = 0;
    virtual real GetSpeciesConc(const unsigned int i, const real t) const = 0;
    virtual real GetSpeciesConc(const string name, const real t) const = 0;
    /* Returns a vector of species concentrations at given time. */
    virtual void GetSpeciesConcs(const real t, vector<real> &chem) const = 0;
    /* Returns a vector of all chemical conditions at given time.  First com
       the species concentrations, then the temperature and then pressure. */
    virtual void GetConditions(const real time, vector<real> &chem, real &T, real &P) const = 0;
    virtual void SetTemperature(const real t) = 0;
    virtual void SetPressure(const real p) = 0;
    virtual void SetSpeciesConc(const unsigned int i, const real c) = 0;
    /* Adds the given value to the required species concentration. */
    virtual void AdjustSpeciesConc(const unsigned int i, const real delta) = 0;
    /* Sets the initial conditions. */
    virtual void SetInitConditions(const real time, const vector<real> &chem, const real T, const real P) = 0;    
public: // Scaling routines.
    /* Returns the real system to stochastic system scaling factor. */
    virtual real SampleVolume() const;
    /* Sets the number density which the ensemble represents. */
    virtual int SetM0(const real m0);
    /* Sets the number density which the full ensemble would represent. */
    virtual int SetMaxM0(const real m0);
public:
    virtual void Reset(const real m0);
};
};

#endif