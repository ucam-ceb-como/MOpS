/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Element class declared in the
    gpc_element.h header file.

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

#include "gpc_element.h"
#include "gpc_mech.h"
#include <stdexcept>
#include <iostream>
#include <string>
#include "string_functions.h"

using namespace Sprog;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Element::Element(void)
{
    m_name  = "";
    m_molwt = 0.0;
    m_mech  = NULL;
}

// Stream-reading constructor.
Element::Element(std::istream &in)
{
    Deserialize(in);
}

// Copy constructor.
Element::Element(const Element &e)
{
    *this = e;
}

// Initialising constructor.
Element::Element(const std::string &name, const double molwt)
{
    m_name  = name;
    m_molwt = molwt;
    m_mech  = NULL;
}

// Destructor.
Element::~Element(void)
{
    // There is nothing special to destruct here.
}


// OPERATOR OVERLOADING.

// Assignment operator.
Element &Element::operator=(const Sprog::Element &el)
{
    if (this != &el) {
        m_name  = el.m_name;
        m_molwt = el.m_molwt;
        m_mech  = el.m_mech;
    }
    return *this;
}

// Comparison operator:  Compares two elements.  Returns true
// if the element names are the same.
bool Element::operator==(const Sprog::Element &el) const
{
    return m_name.compare(el.m_name)==0;
}

// Comparison operator: Compares an element to a string (name).  Returns
// true if the element name matches the string.
bool Element::operator==(const std::string &name) const
{
    return m_name.compare(name)==0;
}

// Inequality operator:  Compares two elements.  Returns false if the
// element names are equal.
bool Element::operator!=(const Sprog::Element &el) const
{
    return !(*this==el);
}

// Inequality operator:  Compares an element to a string (name).  Returns
// false of the element name matches the string.
bool Element::operator!=(const std::string &name) const
{
    return !(*this==name);
}


// NAME.

// Sets the element name.  Also checks the parent mechanism
// to see if an element with that name already exists.  If one does
// then this function throws an exception.
void Element::SetName(const std::string &name)
{
    // If the element is part of a mechanism then we must check
    // that the mechanism does not already contain an element with
    // that name.
    if (m_mech != NULL) {
        if (m_mech->FindElement(name) >= 0) {
            // Oh dear:  mechanism already contains an element with
            // this name.
            throw invalid_argument("Cannot have two elements with the "
                                   "same name (Sprog, Element::SetName).");
        }
    }

    // Element is not in a mechanism or name is not defined
    // in mechanism.  Can set name here in either case.
    m_name = name;
}


// MOLECULAR WEIGHT.

// Sets the element molecular weight.
void Element::SetMolWt(const double molwt)
{
    if (molwt > 0.0) {
        m_molwt = molwt;

        // Tell mechanism that element mol. wt. has changed:  It probably
        // wants to update its species.
        if (m_mech != NULL) {
            m_mech->CheckElementChanges(*this);
        }
    } else {
        // Attempting to set zero or negative mol. wt.  This
        // is unphysical.
        throw out_of_range("Molecular weight must be positive "
                           "and non-zero! (Sprog, Element::SetMolWt).");
    }
}

// Sets the element molecular weight by search the
// list of known elements.
bool Element::SetMolWtFromLibrary()
{
    for (unsigned int i=0; i<m_nlib; i++) {
        if (m_name.compare(m_lib[i].Name()) == 0) {
            // We have found a matching element in the library.
            m_molwt = m_lib[i].MolWt();
            return true;
        }
    }

    // We have not found an element in the library.
    return false;
}


// PARENT MECHANISM.

const Sprog::Mechanism *const Element::Mechanism() const
{
    return m_mech;
}

// Sets the parent mechanism.
void Element::SetMechanism(Sprog::Mechanism &mech)
{
    m_mech = &mech;
}


// READ/WRITE FUNCTIONS.

Element *const Element::Clone(void) const
{
    return new Element(*this);
}

// Prints a diagnostic output file containing all the
// element data.  This is used to debug.
void Element::WriteDiagnostics(std::ostream &out) const
{
    string data = "";

    if (out.good()) {
        // Name.
        out.write(string(m_name+ "\n").c_str(), m_name.length());
        // Mol. Wt.
        data = cstr(m_molwt) + "\n";
        out.write(data.c_str(), data.length());
    }
}

/*!
@param[in]      out         the ostream used to output the elements to the Chemkin mechanism file.
*/
void Element::WriteElements(std::ostream &out) const {

    if (out.good()) {
        // Name.
        out << m_name << " ";
        out << "\n";
    }
}

// Writes the element to a binary data stream.
void Element::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Write the serialisation version number to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write the length of the element name to the stream.
        unsigned int n = m_name.length();
        out.write((char*)&n, sizeof(n));

        // Write the element name to the stream.
        out.write(m_name.c_str(), n);

        // Write the molecular weight to the stream.
        double wt = (double)m_molwt;
        out.write((char*)&wt, sizeof(wt));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sprog, Element::Serialize).");
    }
}

// Reads the element data from a binary data stream.
void Element::Deserialize(std::istream &in)
{
    // Clear the element of its current data.
    m_name  = "";
    m_molwt = 0.0;
    m_mech  = NULL;

    if (in.good()) {
        // Read the serialized element version number.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0; // Need for reading name length.
        char *name = NULL;
        double wt = 0.0;

        switch (version) {
            case 0:
                // Read the length of the element name.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read the element name.
                name = new char[n];
                in.read(name, n);
                m_name.assign(name, n);
                delete [] name;

                // Read the element mol. wt.
                in.read(reinterpret_cast<char*>(&wt), sizeof(wt));
                m_molwt = (double)wt;

                break;
            default:
                throw runtime_error("Element serialized version "
                                    "number is unsupported "
                                    "(Sprog, Element::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sprog, Element::Deserialize).");
    }
}
