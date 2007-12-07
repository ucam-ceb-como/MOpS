/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This file holds data structures to describe the elemental composition of
    chemical species.
*/

#ifndef GPC_EL_COMP_H
#define GPC_EL_COMP_H

#include "gpc_element.h"
#include <vector>

namespace Sprog
{
class ElComp
{
public:
    // Constructors.
    ElComp(void); // Default constructor.
    ElComp(const ElComp &el); // Copy constructor.
    ElComp(unsigned int i, unsigned int n); // Initialising constructor.

    // Destructor.
    virtual ~ElComp(void);

    // Operators.
    ElComp &operator=(const ElComp &el);
    ElComp &operator+=(const unsigned int n);
    ElComp &operator-=(const unsigned int n);
    const ElComp operator+(const unsigned int n);
    const ElComp operator-(const unsigned int n);
    bool operator==(const ElComp &el) const;
    bool operator!=(const ElComp &el) const;

    // Related element object.
    int Index(void) const; // Returns a pointer to the related element.
    void SetIndex(unsigned int el); // Sets the related element object.

    // Number of the element in this composition.
    unsigned int Count(void) const; // Returns number of the element in this composition.
    void SetCount(unsigned int n);  // Sets the element count in this composition.
private:
    // Data.
    int m_index;          // Index of element.
    unsigned int m_count; // Element count in this composition.
};

// A typedef for a STL vector of ElComp objects.
typedef std::vector<ElComp> ElCompVector;
typedef std::vector<ElComp*> ElCompPtrVector;
};

#endif