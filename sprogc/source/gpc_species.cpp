#include "gpc_species.h"

using namespace Sprog;
using namespace std;

// ** THE Species CLASS. **

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Species::Species(void)
{
    m_name = "";
    m_molwt = 0.0;
    m_T1 = 0.0;
}

// Copy constructor.
Species::Species(const Sprog::Species &sp)
{
    *this = sp;
}

// Destructor.
Species::~Species(void)
{
    m_name.clear();
    m_elcomp.clear();
    m_thermoparams.clear();
}


// OPERATOR OVERLOADING.

// Assignment operator.
Species &Species::operator=(const Sprog::Species &sp)
{
    // Remember to check for self-assignment!
    if (this!=&sp) {
        // Copy simple member data.
        m_name = sp.m_name;
        m_elcomp.assign(sp.m_elcomp.begin(), sp.m_elcomp.end());
        m_molwt = sp.m_molwt;
        m_T1 = sp.m_T1;
        
        // Copy thermo parameters.
        map<real,THERMO_PARAMS>::const_iterator i;
        for (i=sp.m_thermoparams.begin(); i!=sp.m_thermoparams.end(); i++) {
            m_thermoparams[i->first] = i->second;
        }
    }
    return *this;
}

// Comparison operator:  Compares two Species objects.  Returns true
// of they have the same name.
bool Species::operator==(const Sprog::Species &sp) const
{
    return m_name.compare(sp.m_name)==0;
}

// Comparison operator:  Compares a species object to a string.  Returns
// true if the species name is the same as the string.
bool Species::operator==(const std::string &name) const
{
    return m_name.compare(name)==0;
}

// Inequality operator:  Returns false of the two species objects have the
// same name.
bool Species::operator!=(const Sprog::Species &sp) const
{
    return !(*this==sp);
}

// Inequality operator:  Returns false if the species name matches
// the string.
bool Species::operator!=(const std::string &name) const
{
    return !(*this==name);
}


// SPECIES COMPOSITION.

// Adds an element to the species composition vector.
void Species::AddElement(const Sprog::ElComp &elcomp)
{
    m_elcomp.insert(m_elcomp.end(), elcomp);
}


// THERMODYNAMIC FITTING PARAMETERS.

// Returns the thermo parameters which are valid for the given temperature.
const THERMO_PARAMS &Species::ThermoParams(const Sprog::real T) const
{
    // The thermo params are indexed by temperature.  The lower_bound() function
    // returns the first point at which the key (Temperature) is greater than
    // the supplied value.
    ThermoMap::const_iterator i = m_thermoparams.lower_bound(T);
    return i->second;
}

// Adds a set of thermo parameters valid up to the given temperature to the
// species object.
void Species::AddParams(const Sprog::real T, const Sprog::THERMO_PARAMS &params)
{
    m_thermoparams[T] = params;
}

// Removes the thermo parameters associated with the given temperature, assuming
// it is in the list.
void Species::RemoveParams(const Sprog::real T)
{
    m_thermoparams.erase(T);
}


// ** THE IndexedSpecies STRUCTURE. **


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
IndexedSpecies::IndexedSpecies(void)
{
    Pointer = NULL;
    Index = 0;
}

// Copy constructor.
IndexedSpecies::IndexedSpecies(const Sprog::IndexedSpecies &sp)
{
    *this = sp;
}

// Destructor.
IndexedSpecies::~IndexedSpecies(void)
{
    Pointer = NULL;
    Index = 0;
}


// OPERATOR OVERLOADING.

// Assignment operator.
IndexedSpecies &IndexedSpecies::operator=(const Sprog::IndexedSpecies &sp)
{
    if (this!=&sp) {
        Pointer = sp.Pointer;
        Index = sp.Index;
    }
    return *this;
}

// Equality operator:  Returns true if both IndexedSpecies objects point
// to the same species object.
bool IndexedSpecies::operator==(const Sprog::IndexedSpecies &sp) const
{
    return (Pointer==sp.Pointer);
}

// Inequality operator:  Returns false if both IndexedSpecies objects point
// to the same species object.
bool IndexedSpecies::operator!=(const Sprog::IndexedSpecies &sp) const
{
    return !(*this==sp);
}