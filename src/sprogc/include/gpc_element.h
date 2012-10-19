/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Element class describes a chemical element.  The file also contains
    typedefs related to Element objects.  Chemical elements have a symbol/name
    and a molecular weight.  They also belong to a mechanism, which is responsible
    for creating, destroying and manipulating elements.  In particular the
    mechanism provides a routine for checking if an element is already defined.  This
    functionality is used when setting the element symbol/name to ensure duplicate
    elements are not defined.

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

#ifndef GPC_ELEMENT_H
#define GPC_ELEMENT_H

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "gpc_params.h"
#include <vector>
#include <string>
#include <iostream>

namespace Sprog
{
class Mechanism; // Forward declaration of mechanism class.

class Element
{
public:
    // Constructors.
    Element(void);             // Default constructor.
    Element(const Element &e); // Copy constructor.
    Element(                   // Initialising constructor.
        const std::string &name,  // - Element name.
        const double molwt          // - Molecular weight.
        );
    Element(std::istream &in);

    // Destructor.
    ~Element(void);

    // Operator overloads.
    Element &operator=(const Element &el);
    bool operator==(const Element &el) const;
    bool operator==(const std::string &name) const;
    bool operator!=(const Element &el) const;
    bool operator!=(const std::string &name) const;


    // ELEMENT NAME.

    // Returns the name of the element.
    const std::string Name(void) const;

    // Sets the name of the element.
    void SetName(const std::string &name);


    // MOLECULAR WEIGHT.

    // Returns molecular weight of the element.
    double MolWt() const;

    // Sets the molecular weight of the element.
    void SetMolWt(const double molwt);

    // Searches for the element in the library of known elements.
    bool SetMolWtFromLibrary();


    // PARENT MECHANISM.

    // Returns pointer to parent mechanism.
    const Sprog::Mechanism *const Mechanism(void) const;

    // Sets the parent mechanism.
    void SetMechanism(Sprog::Mechanism &mech);


    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the element object.
    Element *const Clone(void) const;

    // Writes the element to a binary data stream.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */)
    {
        ar & m_name & m_molwt & m_mech;
    }

    friend class boost::serialization::access;

    // Prints a diagnostic output file containing all the
    // element data.  This is used to debug.
    void WriteDiagnostics(std::ostream &out) const;

    //! Writes the elements to a Chemkin output file.
    void WriteElements(std::ostream &out) const;

    // Writes the element to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the element data from a binary data stream.
    void Deserialize(std::istream &in);

private:
    // Element data.
    std::string m_name;       // Element name/symbol.
    double m_molwt;             // Molecular weight (kg/mol).
    Sprog::Mechanism *m_mech; // Parent mechanism.

    // Library of known elements.
    const static unsigned int m_nlib = 74;
    const static Element m_lib[m_nlib];


};

// Inline function definitions.
#include "gpc_element_inl.h"

// A typedef for a STL vector of elements.
typedef std::vector<Element> ElementVector;

// A typedef for a STL vector of pointers to elements.
typedef std::vector<Element*> ElementPtrVector;
};

#endif // GPC_ELEMENT_H
