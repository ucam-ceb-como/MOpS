/*
  Author(s):      Martin Martin (mm864)
  Project:        sprog (gas-phase and surface chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Martin Martin.

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

#ifndef GPC_SP_COMP_H
#define GPC_SP_COMP_H
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <vector>
#include <map>

namespace Sprog
{
class SpComp
{
public:
    // Constructors.
    SpComp(void);             // Default constructor.
    SpComp(const SpComp &sp); // Copy constructor.
    SpComp(                   // Initialising constructor.
	   unsigned int i, // index
	   unsigned int n  // number of sp 
        );

    // Destructor.
    ~SpComp(void);

    // Operators.
    SpComp &operator=(const SpComp &sp);
    SpComp &operator+=(const unsigned int n);
    SpComp &operator-=(const unsigned int n);
    const SpComp operator+(const unsigned int n);
    const SpComp operator-(const unsigned int n);
    bool operator==(const SpComp &sp) const;
    bool operator!=(const SpComp &sp) const;


    // SPECIES INDEX.

    // Returns the species index.  A negative number indicates an
    // invalid index, which probably means that the SpComp has
    // not been initialise properly.
    int Index(void) const;

    // Sets the related species index.  Note only positive values may
    // be set, even though the Index() function can return negative
    // values.  This ensures that invalid indices cannot be set.
    void SetIndex(unsigned int sp);


    // SPECIES COUNT.

    // Returns number of the species in this phase composition.
    unsigned int Count(void) const;

    // Sets the species count in this phase composition.
    void SetCount(unsigned int n);

    // Writes the species to a binary data stream.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */)
    {
        ar & m_index & m_count;
    }

    friend class boost::serialization::access;

private:
    // Composition data.
    int m_index;          // Index of species.
    unsigned int m_count; // Species count in this composition.
};

// A typedef for a STL vector of SpComp objects.
typedef std::vector<SpComp> SpCompVector;

// A typedef for a STL vector of pointers to SpComp objects.
typedef std::vector<SpComp*> SpCompPtrVector;

// Alternative method of defining element composition.
typedef std::map<unsigned int, unsigned int> SpMap;
/*
typedef std::pair<unsigned int, unsigned int> ElComp;
*/
};

#endif
