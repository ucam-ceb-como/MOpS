/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a particle component.
*/

#ifndef SWEEP_COMPONENT_H
#define SWEEP_COMPONENT_H

#include "swpparams.h"
#include <string>

using namespace std;

namespace Sweep
{
class Component
{
protected:
    real m_molwt;   // Component molecular weight (g/mol).
    real m_density; // Component density (g/cm3).
    string m_name;
public:
    Component(void);
    Component(const real molwt, const real dens, const string &name);
    void Initialise(const real molwt, const real dens, const string &name);
    ~Component(void);
    inline real MolWt() const {return m_molwt;};     // Returns component molecular weight (g/mol).
    inline real Density() const {return m_density;}; // Returns component density (g/cm3).
    inline void SetMolWt(const real molwt) {m_molwt = molwt;};   // Sets the molecular weight (g/mol).
    inline void SetDensity(const real dens) {m_density = dens;}; // Sets the density (g/cm3).
    inline const string &Name() const {return m_name;};
    inline void SetName(const string &name) {m_name = name;};
};
};

#endif