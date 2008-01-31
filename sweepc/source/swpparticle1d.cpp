#include "swpparticle1d.h"
#include "swpcomponent.h"
#include "swpparams.h"
#include <vector>
#include <math.h>

using namespace Sweep;
using namespace std;

/* Constructors and destructors. */

Particle::Particle(void)
{
    m_comp.clear();
    m_values.clear();
    m_createt = 0.0;
    m_components = NULL;
}

Particle::Particle(std::vector<Component*> &components, const unsigned int nvals)
{
    Initialise(components, nvals);
}

Particle::~Particle(void)
{
    m_comp.clear();
    m_values.clear();
    m_components = NULL;
}

void Particle::Initialise(std::vector<Component*> &components, const unsigned int nvals)
{
    m_components = &components;
    m_comp.assign(m_components->size(),0.0);
    m_values.assign(nvals,0.0);
}


/* Overridden routines from base class. */

void Particle::GetProperties(ParticleCache & p) const
{
    vector<real>::iterator i;
    vector<real>::const_iterator j;
    p.Resize(NCACHE+(int)m_comp.size()+(int)m_values.size());
    i = p.GetCache().begin();
    for (j=m_cache.GetCache().begin(); j!=m_cache.GetCache().end(); i++, j++) {*i = *j;}
    for (j=m_comp.begin(); j!=m_comp.end(); i++, j++) {*i = *j;}
    for (j=m_values.begin(); j!=m_values.end(); i++, j++) {*i = *j;}
}

bool Particle::IsValid(void) const
{
    // Particle is considered valid if it has no negative
    // components.

    vector<real>::const_iterator i;
    for (i=m_comp.begin(); i!=m_comp.end(); i++) {
        if (*i < 0.0) return false;
    }
    return true;
}

Particle &Particle::CreateCopy()
{
    Particle *sp = new Particle(*m_components, (int)m_values.size());
    sp->SetComposition(m_comp);
    sp->SetValues(m_values);
    sp->SetCreateTime(m_createt);
    sp->CalcCache();
    return *sp;
}

/* Property get/set routines. */

real Particle::GetProperty(const unsigned int i) const
{
    if (i < m_cache.GetCache().size()) 
        return m_cache[i]; 
    else 
        return 0.0;
}

const vector<real> &Particle::Composition(void) const
{
    return m_comp;
}

void Particle::SetComposition(const std::vector<real> &comp)
{
    m_comp.assign(comp.begin(), comp.end());
}

void Particle::Adjust(const vector<real> &dcomp, const vector<real> &dvalues)
{
    // Adjusts a particle given a change in composition and values.

    vector<real>::iterator i;
    vector<real>::const_iterator j;

    // Apply changes to particle composition.
    for (i=m_comp.begin(), j=dcomp.begin(); 
         (i!=m_comp.end()) && (j!=dcomp.end()); i++, j++) {
        *i += *j;
    }

    // Apply changes to values.
    for (i=m_values.begin(), j=dvalues.begin(); 
         (i!=m_values.end()) && (j!=dvalues.end()); i++, j++) {
        *i += *j;
    }

    if (IsValid()) {
        CalcCache();
    }
};

void Particle::Adjust(const vector<real> &dcomp, const vector<real> &dvalues, const unsigned int n)
{
    // Adjusts a particle n times given a change in composition and values.

    vector<real>::iterator i;
    vector<real>::const_iterator j;

    // Apply changes to particle composition.
    for (i=m_comp.begin(), j=dcomp.begin(); 
         (i!=m_comp.end()) && (j!=dcomp.end()); i++, j++) {
        *i += *j * (real)n;
    }

    // Apply changes to values.
    for (i=m_values.begin(), j=dvalues.begin(); 
         (i!=m_values.end()) && (j!=dvalues.end()); i++, j++) {
        *i += (*j) * (real)n;
    }
};

const vector<real> &Particle::Values() const
{
    return m_values;
}

void Particle::SetValues(const std::vector<real> &values) 
{
    m_values.assign(values.begin(), values.end());
}

real Particle::CreateTime() const 
{
    return m_createt;
}

void Particle::SetCreateTime(Sweep::real t) 
{
    m_createt = t;
}

real Particle::Time(void) const 
{
    return m_time;
}

void Particle::SetTime(const real t) 
{
    m_time = t;
}

/* Particle actions and interactions. */

Particle &Particle::operator=(const Sweep::Particle &sp)
{
    // Definition of assignment (=) operator for particles.

    if (this == &sp) return *this;
    SetTime(sp.Time());
    m_createt = sp.CreateTime();
    m_comp.assign(const_cast<Particle&>(sp).Composition().begin(), const_cast<Particle&>(sp).Composition().end());
    m_values.assign(const_cast<Particle&>(sp).Values().begin(), const_cast<Particle&>(sp).Values().end());
    CalcCache();
    return *this;
}

Particle &Particle::operator +=(const Sweep::Particle &sp)
{
    // Definition of += operator for particles.  This operator is used to define 
    // coagulation.
    return Coagulate(sp);
}

const Particle Particle::operator+(const Sweep::Particle &sp) const
{
    Particle newp = *this;
    newp += sp;
    return newp;
}

Particle &Particle::Coagulate(const Particle &sp)
{
    // Add together particle components.
    vector<real>::iterator j = m_comp.begin();
    vector<real>::const_iterator k = sp.Composition().begin();
    for(j=m_comp.begin(); j!=m_comp.end(); j++,k++) *j += *k;

    // Add together particle values.
    j = m_values.begin();
    k = sp.Values().begin();
    for(j=m_values.begin(); j!=m_values.end(); j++,k++) *j += *k;

    // Particle create time is earliest time.
    m_createt = min(m_createt, sp.CreateTime());

    // Calculate particle's property cache, and return it.
    CalcCache();
    return *this;
}


/* Protected functions. */

void Particle::CalcCache(void)
{
    // Clear current cache.
	m_cache.Clear();

    real m;
    vector<Component*>::iterator ic;
    vector<real>::iterator j = m_comp.begin();

    // Loop over composition and calculate mass and volume.
    for (ic=m_components->begin(); ic!=m_components->end(); ic++, j++) {
        m = (*ic)->MolWt() * (*j) / NA;
        m_cache[m_cache.iM] += m;
        m_cache[m_cache.iV] += m / (*ic)->Density();
    }
    m_cache[m_cache.iD]       = pow(6.0 * m_cache[m_cache.iV] / PI, ONE_THIRD);
    m_cache[m_cache.iD2]      = m_cache[m_cache.iD] * m_cache[m_cache.iD];
    m_cache[m_cache.iS]       = PI * m_cache[m_cache.iD2];
    m_cache[m_cache.iD_1]     = 1.0 / m_cache[m_cache.iD];
    m_cache[m_cache.iD_2]     = 1.0 / m_cache[(m_cache.iD2];
    m_cache[m_cache.iM_1_2]   = pow(m_cache[m_cache.iM], -0.5);
    m_cache[m_cache.iD2M_1_2] = m_cache[m_cache.iD2] * m_cache[m_cache.iM_1_2];
}

