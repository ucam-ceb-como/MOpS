/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a 2D surface-volume particle used by sweep.  Particles
    are described by their volume and surface area.  All processes lead to
    a rounding of particles except coagulation.  On coagulation surface area
    is conserved, i.e point contact.
*/

#ifndef SWEEP_SVPARTICLE_H
#define SWEEP_SVPARTICLE_H

#include <vector>
#include "swpparticle1d.h"

namespace Sweep
{
class SVParticle : public Sweep::DefaultParticle
{
public:
    SVParticle(void);
    ~SVParticle(void);
public:
    /* Adjusts the particle with the given changes in composition and values. */
    virtual void Adjust(const vector<real> &dcomp, const vector<real> &dvalues);
    /* Same as Adjust but applies it n times. */
    virtual void Adjust(const vector<real> &dcomp, const vector<real> &dvalues, const unsigned int n);
    /* Coagulates a particle with this one. */
    virtual DefaultParticle &Coagulate(const DefaultParticle &sp);
    virtual void CalcCache(void); // Calculates the cache of particle properties.
};
};

#endif