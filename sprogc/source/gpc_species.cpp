/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Species class declared in the
    gpc_species.h header file.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#include "gpc_species.h"
#include "gpc_el_comp.h"
#include "gpc_mech.h"
#include <stdexcept>
#include <string>

using namespace Sprog;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Species::Species(void)
{
    m_name     = "";
    m_elcomp.clear();
    m_molwt    = 0.0;
    m_mech     = NULL;
    m_thermoparams.clear();
    m_T1       = 0.0;
}

// Copy constructor.
Species::Species(const Sprog::Species &sp)
{
    m_mech = NULL;
    *this  = sp;
}

// Stream-reading constructor.
Species::Species(istream &in)
{
    Deserialize(in);
}

// Destructor.
Species::~Species(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator.
Species &Species::operator=(const Sprog::Species &sp)
{
    // Remember to check for self-assignment!
    if (this!=&sp) {
        // Copy simple member data.
        m_name  = sp.m_name;
        m_molwt = sp.m_molwt;
        m_mech  = sp.m_mech;
        m_T1    = sp.m_T1;

        // Copy composition.
        m_elcomp.assign(sp.m_elcomp.begin(), sp.m_elcomp.end());

        // Copy thermo parameters.
        Sprog::Thermo::ThermoMap::const_iterator i;
        for (i=sp.m_thermoparams.begin(); i!=sp.m_thermoparams.end(); i++) {
            m_thermoparams[i->first] = i->second;
        }
    }
    return *this;
}

// Comparison operator:  Compares two Species objects.  Returns true
// if they have the same name.
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


// SPECIES NAME.

// Sets the species name.
void Species::SetName(const std::string &name)
{
    // If this species is part of a mechanism then we must check
    // if a species with this name is already defined.
    if (m_mech != NULL) {
        if (m_mech->FindSpecies(name) >= 0) {
            // Oh dear:  Species with this name already defined.  We'd
            // better throw an error.
            throw invalid_argument("Cannot have two species with "
                                   "the same name (Sprog, Species::SetName).");
        }
    }

    // It is ok to set the name now.
    m_name = name;
}


// SPECIES COMPOSITION.

// Returns the total number of atoms in the species.
unsigned int Species::AtomCount(void) const
{
    // Sum up all ElComps.
    unsigned int n = 0;
    for (ElCompVector::const_iterator i=m_elcomp.begin(); i!=m_elcomp.end(); i++) {
        n += (*i).Count();
    }
    return n;
}

// Returns the number of the given element in the species.
unsigned int Species::AtomCount(unsigned int iel) const
{
    // Locate the ElComp with the given index.
    for (ElCompVector::const_iterator i=m_elcomp.begin(); i!=m_elcomp.end(); i++) {
        if (iel == (unsigned)i->Index()) {
            return i->Count();
        }
    }

    // Element not part of species.
    return 0;
}

// Adds an element to the species composition vector.
void Species::AddElement(const Sprog::ElComp &elcomp)
{
    if (m_mech != NULL) {
        if ((elcomp.Index() >= 0) && ((unsigned)elcomp.Index() < m_mech->ElementCount())) {
            // Must check if this element is already defined for this
            // species.
            ElCompVector::iterator i;
            bool found = false;
            for (i=m_elcomp.begin(); i!=m_elcomp.end(); i++) {
                if ((*i) == elcomp) {
                    // Found element:  Append value.
                    (*i) += elcomp.Count();
                    found = true;
                    break;
                }
            }

            // Element not currently defined for species, so add it.
            if (!found) m_elcomp.push_back(elcomp);

            // Update species molecular weight:  Sum of element weights.
            m_molwt += (m_mech->Elements(elcomp.Index())->MolWt() * (real)elcomp.Count());
        } else {
            // Element index in elcomp is out-of-range.
            throw out_of_range("Element index is out of range "
                               "(Sprog, Species::AddElement).");
        }
    } else {
        // Parent mechanism has not been assigned yet.
        throw logic_error("Attempt to add element before "
                          "mechanism is assigned (Sprog, Species::AddElement).");
    }
}

// Adds an element to the species given the index and count.
void Species::AddElement(unsigned int i, unsigned int n)
{
    if (m_mech != NULL) {
        if (i < m_mech->ElementCount()) {
            // Must check if this element is already defined for this
            // species.
            ElCompVector::iterator j;
            bool found = false;
            for (j=m_elcomp.begin(); j!=m_elcomp.end(); j++) {
                if ((unsigned)(*j).Index() == i) {
                    // Found element:  Append value.
                    (*j) += n;
                    found = true;
                    break;
                }
            }

            // Element not currently defined for species, so add it.
            if (!found) m_elcomp.push_back(ElComp(i, n));

            // Update species molecular weight:  Sum of element weights.
            m_molwt += (m_mech->Elements(i)->MolWt() * (real)n);
        } else {
            // Element index is out-of-range.
            throw out_of_range("Element index is out of range "
                               "(Sprog, Species::AddElement).");
        }
    } else {
        // Mechanism has not been assigned yet.
        throw logic_error("Attempt to add element before parent "
                          "mechanism is assigned (Sprog, Species::AddElement).");
    }
}

// Adds an element to the species given the element name and count.
void Species::AddElement(const std::string &name, unsigned int n)
{
    if (m_mech != NULL) {
        // Must find the index of the element in the list.
        int i = m_mech->FindElement(name);

        if (i >= 0) {
            // Found element in the list!

            // Must check if this element is already defined for this
            // species.
            ElCompVector::iterator j;
            bool found = false;
            for (j=m_elcomp.begin(); j!=m_elcomp.end(); j++) {
                if ((*j).Index() == i) {
                    // Found element:  Append value.
                    (*j) += n;
                    found = true;
                    break;
                }
            }

            // Element not currently defined for species, so add it.
            if (!found) m_elcomp.push_back(ElComp(i, n));

            // Update species molecular weight:  Sum of element weights.
            m_molwt += (m_mech->Elements(i)->MolWt() * (real)n);
        } else {
            // We have got here because the element wasn't found in the list.
            throw invalid_argument(string(name).append(" not found in "
                                   "element list (Sprog, Species::AddElement)."));
        }
    } else {
        // Element vector has not been assigned yet.
        throw logic_error("Attempt to add element before parent "
                          "mechanism is assigned (Sprog, Species::AddElement).");
    }
}

// Returns true if the species contains the element (given by index).
bool Species::ContainsElement(unsigned int i) const 
{
    // Loop over composition to find element.
    ElCompVector::const_iterator el;
    for (el=m_elcomp.begin(); el!=m_elcomp.end(); el++) {
        if ((unsigned)el->Index() == i) {
            return true;
        }
    }

    // We have arrived here because the element wasn't found.
    return false;
}

// Returns true if the species contains the element (given by reference).
bool Species::ContainsElement(const Sprog::Element &el) const
{
    if (m_mech != NULL) {
        // Loop over the elements to find the index.
        unsigned int i;
        for (i=0; i<m_mech->ElementCount(); i++) {
            if (*m_mech->Elements(i) == el) {
                return ContainsElement(i);
            }
        }

        // We have got here because the element isn't in the list, therefore
        // the species cannot contain it.
        return false;
    } else {
        // The parent mechanism is unassigned.
        throw logic_error("Parent mechanism is unassigned "
                          "(Sprog, Species::ContainsElement).");
    }
}

// Returns true if the species contains the element (given by name).
bool Species::ContainsElement(const std::string &name) const
{
    if (m_mech != NULL) {
        // Loop over the elements to find the index.
        unsigned int i;
        for (i=0; i<m_mech->ElementCount(); i++) {
            if (*m_mech->Elements(i) == name) {
                return ContainsElement(i);
            }
        }

        // We have got here because the element isn't in the list, therefore
        // the species cannot contain it.
        return false;
    } else {
        // The elements vector is unassigned.
        throw logic_error("Parent mechanism is unassigned "
                          "(Species::ContainsElement).");
    }
}


// MOLECULAR WEIGHT.

// Recalculates the species molecular weight using the elements.
real Species::CalcMolWt()
{
    m_molwt = 0.0; // Reset.

    if (m_mech != NULL) {
        // Loop over composition vector, summing up the molecular weight.
        ElCompVector::const_iterator el;
        for (el=m_elcomp.begin(); el!=m_elcomp.end(); el++) {
            m_molwt += m_mech->Elements(el->Index())->MolWt() * (real)el->Count();
        }
        return m_molwt;
    } else {
        // Parent mechanism not assigned.
        throw logic_error("Cannot calculate molecular weight before "
                          "assigning parent mechanism (Sprog, Species::CalcMolWt).");
    }
}


/*
// ELEMENTS VECTOR.

// Returns a pointer to the vector of elements used to define the species.
const ElementPtrVector *const Species::Elements()
{
    return m_elements;
}

// Sets the pointer to the vector of elements used to define the species.
void Species::SetElements(const Sprog::ElementPtrVector *const els)
{
    m_elements = els;
}
*/


// THERMODYNAMIC FITTING PARAMETERS.

// Returns the thermo parameters which are valid for the given temperature.
const Thermo::THERMO_PARAMS &Species::ThermoParams(const Sprog::real T) const
{
    if (T >= m_T1) {
        // The thermo params are indexed by temperature.  The lower_bound() function
        // returns the first point at which the key (Temperature) is greater than
        // the supplied value.
        Thermo::ThermoMap::const_iterator i = m_thermoparams.lower_bound(T);

        if (i != m_thermoparams.end()) {
            // Success.
            return i->second;
        } else {
            // T is to large!
            return (m_thermoparams.end()--)->second; // Return last set.
//            throw out_of_range("Temperature above acceptable range "
//                               "(Sprog, Species::ThermoParams)");
        }
    } else {
        // The given temperature is lower than the thermo range start temperature.
        return m_thermoparams.begin()->second; // Return first set.
//        throw out_of_range("Temperature under acceptable range "
//                           "(Sprog, Species::ThermoParams)");
    }
}

// Adds a set of thermo parameters valid up to the given temperature to the
// species object.
void Species::AddThermoParams(const Sprog::real T, 
                              const Sprog::Thermo::THERMO_PARAMS &params)
{
    m_thermoparams[T] = params;
}

// Removes the thermo parameters associated with the given temperature, assuming
// it is in the list.
void Species::RemoveThermoParams(const Sprog::real T)
{
    m_thermoparams.erase(T);
}


// PARENT MECHANISM.

// Sets the parent mechanism.
void Species::SetMechanism(Sprog::Mechanism &mech)
{
    m_mech = &mech;
    CalcMolWt();
}


// READ/WRITE FUNCTIONS.

// Writes the species to a binary data stream.
void Species::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Write the serialisation version number to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write the length of the species name to the stream.
        unsigned int n = m_name.length();
        out.write((char*)&n, sizeof(n));

        // Write the species name to the stream.
        out.write(m_name.c_str(), n);

        // Write the length of the composition vector to the stream.
        n = m_elcomp.size();
        out.write((char*)&n, sizeof(n));

        // Write the elemental composition to the stream.
        for (ElCompVector::const_iterator i=m_elcomp.begin(); i!=m_elcomp.end(); i++) {
            // Write the element index.
            int ix = (*i).Index();
            out.write((char*)&ix, sizeof(ix));
            
            // Write the element count.
            n = (*i).Count();
            out.write((char*)&n, sizeof(n));
        }

        // Write the molecular weight to the stream.
        double wt = (double)m_molwt;
        out.write((char*)&wt, sizeof(wt));

        // Write the number of thermo params to the stream.
        n = m_thermoparams.size();
        out.write((char*)&n, sizeof(n));

        // Write the thermo params to the stream.
        for (Thermo::ThermoMap::const_iterator i=m_thermoparams.begin(); i!=m_thermoparams.end(); i++) {
            // Write the temperature.
            double T = (double)(i->first);
            out.write((char*)&(T), sizeof(T));

            // Write the number of thermo params.
            n = i->second.Count;
            out.write((char*)&n, sizeof(n));

            // Write the thermo params.
            for (unsigned int j=0; j<n; j++) {
                double param = (double)(i->second.Params[j]);
                out.write((char*)&param, sizeof(param));
            }
        }

        // Output the start temperature for the thermo range.
        double T = (double)m_T1;
        out.write((char*)&T, sizeof(T));

    } else {
        throw invalid_argument("Output stream not ready (Sprog, Species::Serialize).");
    }
}

// Reads the species data from a binary data stream.
void Species::Deserialize(std::istream &in)
{
    // Clear the species of all current data.
    m_name  = "";
    m_molwt = 0.0;
    m_mech  = NULL;
    m_T1    = 0.0;
    m_elcomp.clear();
    m_thermoparams.clear();

    if (in.good()) {
        // Read the serialized species version number.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0; // Need for reading name length.
        char *name = NULL;
        double T =0.0, wt=0.0;

        switch (version) {
            case 0:
                // Read the length of the species name.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                
                // Read the species name.
                name = new char[n];
                in.read(name, n);
                m_name.assign(name, n);
                delete [] name;

                // Read length of composition vector and reserve memory.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_elcomp.reserve(n);

                // Read the elemental composition.
                for (unsigned int i=0; i<n; i++) {
                    // Read the element index.
                    int ix = -1;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Read the element count.
                    unsigned int m = 0;
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));

                    // Add new ElComp object to vector.
                    m_elcomp.push_back(ElComp(ix, m));
                }

                // Read the species mol. wt.
                in.read(reinterpret_cast<char*>(&wt), sizeof(wt));
                m_molwt = (real)wt;

                // Read the number of thermo parameters .
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read the thermo parameters.
                for (unsigned int i=0; i<n; i++) {
                    Thermo::THERMO_PARAMS params;

                    // Read the temperature.
                    in.read(reinterpret_cast<char*>(&T), sizeof(T));

                    // Read the number of thermo params.
                    unsigned int m = 0;
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));
                    params.Count = m;

                    // Read the parameters.
                    for (unsigned int j=0; j<m; j++) {
                        double param = 0.0;
                        in.read(reinterpret_cast<char*>(&param), sizeof(param));
                        params.Params[j] = (real)param;
                    }

                    // Add the parameters to the thermo map.
                    m_thermoparams[(real)T] = params;
                }

                // Read the start temperature for the thermo range.
                in.read(reinterpret_cast<char*>(&T), sizeof(T));
                m_T1 = (real)T;

                break;
            default:
                throw runtime_error("Species serialized version number "
                                    "is unsupported (Sprog, Species::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Sprog, Species::Deserialize).");
    }
}
