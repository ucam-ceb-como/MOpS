/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Class used to change particle composition and other properties.  Processes
    which change particles should inherit this.
*/

#ifndef SWEEP_PARTICLECHANGER_H
#define SWEEP_PARTICLECHANGER_H

#include "swpparams.h"
#include "swpcomponent.h"
#include "swpparticle1d.h"
#include <vector>

using namespace std;

namespace Sweep
{
class ParticleChanger
{
protected:
    vector<real> m_comp, m_values;    // Component counts and values for new particle.
    vector<Component*> *m_components; // Pointer to components for which this object is defined.
public:
    ParticleChanger(void);
    ~ParticleChanger(void);
public:
    /* Sets the reference to the components used to define the changer. */
    void SetComponents(vector<Component*> &components);
    /* Sets the component change if component i. */
    void SetCompChange(const unsigned int i, const real dcomp);
    /* Sets the value change for value i. */
    void SetValueChange(const unsigned int i, const real value);
public:
    /* Returns reference to the components used to define this changer. */
    inline vector<Component*> &GetComponents(void) const {return *m_components;}
    /* Returns the defined change in component i. */
    inline real GetCompChange(const unsigned int i) const;
    /* Returns the defined change in value i. */
    inline real GetValueChange(const unsigned int i) const;
public:
    /* Changes the particle composition and values by the amounts stored in
       the changer. The changes are applied n times. */
    void AdjustParticle(Particle &sp, const real t, const unsigned int n) const;
    /* Sets the particle composition and values to be the amounts stored in
       the changer. */
    void SetParticle(Particle &sp, const real t) const;
};

inline void ParticleChanger::AdjustParticle(Particle &sp, const Sweep::real t, 
                                            const unsigned int n) const
{
    sp.SetTime(t);
    sp.Adjust(m_comp, m_values, n);
    sp.CalcCache();
}

inline void ParticleChanger::SetParticle(Sweep::Particle &sp, const Sweep::real t) const
{
    sp.SetTime(t);
    sp.SetComposition(m_comp);
    sp.SetValues(m_values);
    sp.CalcCache();
}
};

#endif