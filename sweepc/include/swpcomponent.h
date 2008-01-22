/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a particle component.  A component is a chemical element
    or molecule which is the smallest repeating unit in a particle.  A 
    particle can consist of more than one type of component, though the assumption
    is generally that each component exists in discrete phases rather than being
    well mixed.  This makes calculating particle volume and mass much easier.
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
    string m_name;  // Component symbol or name.
public:
    /* Default constructor. */
    Component(void);
    /* Initialising constructor. */
    Component(const real molwt, const real dens, const string &name);
    /* Initialisation routine which takes as arguments all properties of a component. */
    void Initialise(const real molwt, const real dens, const string &name);
    /* Default destructor. */
    ~Component(void);
public: // Property gets.
    inline real MolWt() const {return m_molwt;};        // Returns component molecular weight (g/mol).
    inline real Density() const {return m_density;};    // Returns component density (g/cm3).
    inline const string &Name() const {return m_name;}; // Returns component symbol or name.
    inline void SetMolWt(const real molwt) {m_molwt = molwt;};   // Sets the molecular weight (g/mol).
    inline void SetDensity(const real dens) {m_density = dens;}; // Sets the density (g/cm3).
    inline void SetName(const string &name) {m_name = name;};    // Sets the symbol or name.
};
};

#endif