/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    The Element class describes a chemical element.  The file also contains
    typedefs related to Element objects.
*/

#ifndef GPC_ELEMENT_H
#define GPC_ELEMENT_H

#include "gpc_params.h"
#include <vector>
#include <string>

namespace Sprog
{
class Mechanism; // Forward declaration of mechanism class.

class Element
{
public:
    // Constructors.
    Element(void); // Default constructor.
    Element(const Element &e); // Copy constructor.    
    Element(const std::string &name, const real molwt); // Initialising constructor.
    
    // Destructor.
    virtual ~Element(void);

    // Operator overloads.
    Element &operator=(const Element &el);
    bool operator==(const Element &el) const;
    bool operator==(const std::string &name) const;
    bool operator!=(const Element &el) const;
    bool operator!=(const std::string &name) const;
    
    // Element name.
    const std::string Name(void) const;    // Returns the name of the element.
    void SetName(const std::string &name); // Sets the name of the element.

    // Element molecular weight.
    real MolWt() const;              // Returns molecular weight of the element.
    void SetMolWt(const real molwt); // Sets the molecular weight of the element.
    bool SetMolWtFromLibrary();      // Searches for the element in the library of known elements.

    // Parent mechanism.
    Sprog::Mechanism *const Mechanism(void);         // Returns pointer to parent mechanism.
    void SetMechanism(Sprog::Mechanism *const mech); // Sets the parent mechanism.

    // Cloning.
    virtual Element *const Clone(void) const; // Returns a pointer to a copy of the Element object.

protected:
    // Element data.
    std::string m_name;       // Element name/symbol.
    real m_molwt;             // Molecular weight (kg/mol).
    Sprog::Mechanism *m_mech; // Parent mechanism.

    // Library.
    const static unsigned int m_nlib = 5;
    const static Element m_lib[m_nlib];
};

// Inline function definitions.
#include "gpc_element_inl.h"

// A typedef for a STL vector of elements.
typedef std::vector<Element> ElementVector;
typedef std::vector<Element*> ElementPtrVector;
};

#endif // GPC_ELEMENT_H
