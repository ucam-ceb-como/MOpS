/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    The Element class describes a chemical element.  The file also contains
    typedefs related to Element objects.  Chemical elements have a symbol/name
    and a molecular weight.  They also belong to a mechanism, which is responsible
    for creating, destroying and manipulating elements.  In particular the
    mechanism provides a routine for checking if an element is already defined.  This
    functionality is used when setting the element symbol/name to ensure duplicate 
    elements are not defined.
*/

#ifndef GPC_ELEMENT_H
#define GPC_ELEMENT_H

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
    Element(std::istream &in); // Stream-reading constructor.
    Element(                   // Initialising constructor.
        const std::string &name,  // - Element name.
        const real molwt          // - Molecular weight.
        ); 

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
    real MolWt() const;

    // Sets the molecular weight of the element.
    void SetMolWt(const real molwt);

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
    void Serialize(std::ostream &out) const;

    // Reads the element data from a binary data stream.
    void Deserialize(std::istream &in);

protected:
    // Element data.
    std::string m_name;       // Element name/symbol.
    real m_molwt;             // Molecular weight (kg/mol).
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
