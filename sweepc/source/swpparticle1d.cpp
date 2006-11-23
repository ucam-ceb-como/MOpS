#include "swpparticle1d.h"
#include "swpcomponent.h"
#include "swpparams.h"
#include <vector>
#include <math.h>

using namespace Sweep;
using namespace std;

/* Constructors and destructors. */

DefaultParticle::DefaultParticle(void)
{
    m_comp.clear();
    m_values.clear();
    m_createt = 0.0;
    m_components = NULL;
    m_cache.resize(NCACHE,0.0);
}

DefaultParticle::DefaultParticle(std::vector<Component*> &components, const unsigned int nvals)
{
    Initialise(components, nvals);
}

DefaultParticle::~DefaultParticle(void)
{
    m_comp.clear();
    m_values.clear();
    m_components = NULL;
    m_cache.clear();
}

void DefaultParticle::Initialise(std::vector<Component*> &components, const unsigned int nvals)
{
    m_components = &components;
    m_comp.assign(m_components->size(),0.0);
    m_values.assign(nvals,0.0);
    m_cache.resize(NCACHE,0.0);
}

/* Overridden routines from base class. */

void DefaultParticle::GetProperties(std::vector<real> &p) const
{
    vector<real>::iterator i;
    vector<real>::const_iterator j;
    p.resize(NCACHE+(int)m_comp.size()+(int)m_values.size());
    i = p.begin();
    for (j=m_cache.begin(); j!=m_cache.end(); i++, j++) {*i = *j;}
    for (j=m_comp.begin(); j!=m_comp.end(); i++, j++) {*i = *j;}
    for (j=m_values.begin(); j!=m_values.end(); i++, j++) {*i = *j;}
}

bool DefaultParticle::IsValid(void) const
{
    vector<real>::const_iterator i;
    for (i=m_comp.begin(); i!=m_comp.end(); i++) {
        if (*i < 0.0) return false;
    }
    return true;
}

DefaultParticle &DefaultParticle::CreateCopy()
{
    DefaultParticle *sp = new DefaultParticle(*m_components, (int)m_values.size());
    sp->SetComposition(m_comp);
    sp->SetValues(m_values);
    sp->SetCreateTime(m_createt);
    sp->CalcCache();
    return *sp;
}

/* Property get/set routines. */

real DefaultParticle::GetProperty(const unsigned int i) const
{
    if (i < m_cache.size()) 
        return m_cache[i]; 
    else 
        return 0.0;
}

const vector<real> &DefaultParticle::Composition(void) const
{
    return m_comp;
}

void DefaultParticle::SetComposition(const std::vector<real> &comp)
{
    m_comp.assign(comp.begin(), comp.end());
}

void DefaultParticle::Adjust(const vector<real> &dcomp, const vector<real> &dvalues)
{
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

void DefaultParticle::Adjust(const vector<real> &dcomp, const vector<real> &dvalues, const unsigned int n)
{
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

const vector<real> &DefaultParticle::Values() const
{
    return m_values;
}

void DefaultParticle::SetValues(const std::vector<real> &values) 
{
    m_values.assign(values.begin(), values.end());
}

real DefaultParticle::CreateTime() const 
{
    return m_createt;
}

void DefaultParticle::SetCreateTime(Sweep::real t) 
{
    m_createt = t;
}

real DefaultParticle::Time(void) const 
{
    return m_time;
}

void DefaultParticle::SetTime(const real t) 
{
    m_time = t;
}

/* Particle actions and interactions. */

DefaultParticle &DefaultParticle::operator=(const Sweep::DefaultParticle &sp)
{
    if (this == &sp) return *this;
    SetTime(sp.Time());
    m_createt = sp.CreateTime();
    m_comp.assign(const_cast<DefaultParticle&>(sp).Composition().begin(), const_cast<DefaultParticle&>(sp).Composition().end());
    m_values.assign(const_cast<DefaultParticle&>(sp).Values().begin(), const_cast<DefaultParticle&>(sp).Values().end());
    CalcCache();
    return *this;
}

DefaultParticle &DefaultParticle::operator +=(const Sweep::DefaultParticle &sp)
{
    return Coagulate(sp);
}

const DefaultParticle DefaultParticle::operator+(const Sweep::DefaultParticle &sp) const
{
    DefaultParticle newp = *this;
    newp += sp;
    return newp;
}

DefaultParticle &DefaultParticle::Coagulate(const DefaultParticle &sp)
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

void DefaultParticle::CalcCache(void)
{
    // Clear current cache.
    m_cache.assign(m_cache.size(), 0.0);

    real m;
    vector<Component*>::iterator ic;
    vector<real>::iterator j = m_comp.begin();
    for (ic=m_components->begin(); ic!=m_components->end(); ic++, j++) {
        m = (*ic)->MolWt() * (*j) / NA;
        m_cache[iM] += m;
        m_cache[iV] += m / (*ic)->Density();
    }
    m_cache[iD]       = pow(6.0 * m_cache[iV] / PI, ONE_THIRD);
    m_cache[iD2]      = m_cache[iD] * m_cache[iD];
    m_cache[iS]       = PI * m_cache[iD2];
    m_cache[iD_1]     = 1.0 / m_cache[iD];
    m_cache[iD_2]     = 1.0 / m_cache[iD2];
    m_cache[iM_1_2]   = pow(m_cache[iM], -0.5);
    m_cache[iD2M_1_2] = m_cache[iD2] * m_cache[iM_1_2];
}

