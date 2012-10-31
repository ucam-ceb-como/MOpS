/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This file holds data structures to describe the elemental composition of
    chemical species.

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

#ifndef GPC_EL_COMP_H
#define GPC_EL_COMP_H
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
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

    // Writes the element to a binary data stream.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */)
    {
      ar & m_index & m_count; // setting for m_elcomp in gpc_species.h
    }

    friend class boost::serialization::access;

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
