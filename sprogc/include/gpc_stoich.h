/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This files contains data structures which describe the stoichiometry of chemical
    reactions.

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

#ifndef GPC_STOICH_H
#define GPC_STOICH_H

#include "gpc_species.h"
#include "gpc_params.h"
#include <map>

namespace Sprog
{
// Structure to hold species stoichiometry for a reaction.
template<class T>
class Stoichiometry
{
public:
    typedef T stoich_val; // Stoich value type.

    // Constructors.
    Stoichiometry(void); // Default constructor.
    Stoichiometry(const Stoichiometry<T> &s); // Copy constructor.
    Stoichiometry(unsigned int isp, const stoich_val &mu); // Initialising constructor.

    // Destructor.
    ~Stoichiometry(void);

    // Operator overloads.
    Stoichiometry<T> &operator=(const Stoichiometry<T> &s);

    // Species data.
    int Index(void) const; // Returns index of species.
    void SetSpecies(const unsigned int &sp); // Sets the species associated with this stoichiometry.

    // Stoichiometry value.
    const stoich_val &Mu(void) const; // Returns the stoichiometry value.
    void SetMu(const stoich_val &mu); // Sets the stoichiometry value.
    void IncMu(const stoich_val &mu); // Increments the stoichiometry value.
private:
    // Data.
    int m_species;
    stoich_val m_stoich;
};

// Typedefs for basic stoichiometry types.
typedef Stoichiometry<int> Stoich;   // Integer stoichiometry data.
typedef Stoichiometry<real> Stoichf; // Real stoichiometry data.

// Typedefs and structure for species/reaction stoichiometry cross-referencing.
typedef std::map<unsigned int, real> RxnStoichMap;
typedef std::pair<unsigned int, real> RxnStoichPair;
struct StoichXRef
{
    unsigned int Species;
    RxnStoichMap RxnStoich;
};
typedef std::vector<StoichXRef> StoichXRefVector;

// Alternative method of defining stoichiometry.
typedef std::map<unsigned int, int> StoichMap;
typedef std::map<unsigned int, real> StoichMapf;
/*
typedef std::pair<unsigned int, int> Stoich;
typedef std::pair<unsigned int, real> Stoichf;
*/
};

#endif
