/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This file holds data structures to describe the elemental composition of
    chemical species.
*/

#ifndef GPC_EL_COMP_H
#define GPC_EL_COMP_H

#include <vector>
#include <map>

namespace Sprog
{
class ElComp
{
public:
    // Constructors.
    ElComp(void);             // Default constructor.
    ElComp(const ElComp &el); // Copy constructor.
    ElComp(                   // Initialising constructor.
        unsigned int i, 
        unsigned int n
        ); 

    // Destructor.
    ~ElComp(void);

    // Operators.
    ElComp &operator=(const ElComp &el);
    ElComp &operator+=(const unsigned int n);
    ElComp &operator-=(const unsigned int n);
    const ElComp operator+(const unsigned int n);
    const ElComp operator-(const unsigned int n);
    bool operator==(const ElComp &el) const;
    bool operator!=(const ElComp &el) const;


    // ELEMENT INDEX.
    
    // Returns the element index.  A negative number indicates an
    // invalid index, which probably means that the ElComp has
    // not been initialise properly.
    int Index(void) const;

    // Sets the related element index.  Note only positive values may
    // be set, even though the Index() function can return negative
    // values.  This ensures that invalid indices cannot be set.
    void SetIndex(unsigned int el);


    // ELEMENT COUNT.
    
    // Returns number of the element in this composition.
    unsigned int Count(void) const;
    
    // Sets the element count in this composition.
    void SetCount(unsigned int n);

private:
    // Composition data.
    int m_index;          // Index of element.
    unsigned int m_count; // Element count in this composition.
};

// A typedef for a STL vector of ElComp objects.
typedef std::vector<ElComp> ElCompVector;

// A typedef for a STL vector of pointers to ElComp objects.
typedef std::vector<ElComp*> ElCompPtrVector;


// Alternative method of defining element composition.
typedef std::map<unsigned int, unsigned int> ElementMap;
/*
typedef std::pair<unsigned int, unsigned int> ElComp;
*/
};

#endif
