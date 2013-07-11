/*
  Author(s):      Martin Martin (mm864)
  Project:        sprog (gas-phase and surface chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2012 Martin Martin.

  File purpose:
    Implementation of the Stoichiometry class declared in the
    gpc_stoich.h header file.

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
#ifndef GPC_DELTA_STOICH_H
#define GPC_DELTA_STOICH_H
#include <vector>
#include <string>
#include <iostream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include "gpc_params.h"

namespace Sprog
{
namespace Kinetics
{
 class Reaction; // Forward declaration of reaction. 
 class DeltaStoich
{
 public: 
  DeltaStoich(void);              // Default constructor.
  DeltaStoich(const DeltaStoich &delta_stoich); // Copy constructor.
  DeltaStoich(std::string &name); // Initialising constructor.
  DeltaStoich(std::istream &in);

  // Destructor
    ~DeltaStoich(void);

 // Operator overloads.
    DeltaStoich &operator=(const DeltaStoich &delta_stoich);

    bool operator==(const std::string &name) const;

   // Returns the species name. (inline)
    const std::string &Name(void) const;

    //Returns the total stoichiometry change of this species. (inline)
    const double &TotalStoich(void) const;

   //Returns the reactant stoichiometry change of this species. (inline)
    const double &ReacStoich(void) const;

    //Returns the product stoichiometry change of this species. (inline)
    const double &ProdStoich(void) const;

    // Sets the species name.
    void SetSpeciesName(const std::string &name);

     // Sets the phase name (type) and store in a vector.
    void SetPhaseName(const std::string &name);

    // Sets/change the total stoichiometry change of this species. 
    void IncrementTotalStoich(const double &value); // can be + or -
    
    // Sets/change the reactant stoichiometry change of this species. 
    void IncrementReacStoich(const double &value); // can be + or -
    
   // Sets/change the product stoichiometry change of this species. 
    void IncrementProdStoich(const double &value); // can be + or -

    // Sets the parent reaction.
    void SetReaction(Sprog::Kinetics::Reaction &rct);


    // Creates a copy of the delta stoich object. (inline)
    DeltaStoich *const Clone(void) const;

    // Prints a diagnostic output file containing all the
    // delta stoich data.  This is used to debug.
    void WriteDiagnostics(std::ostream &out) const;

  // Writes the delta stoich to a binary data stream.
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */)
    {
      ar & m_spName & total_stoich_sp & reac_stoich_sp & prod_stoich_sp /*& m_react*/;
    }
	
	friend class boost::serialization::access;
  // Prints a diagnostic output file containing all the stoichiometry change for each phase.  This is used to debug.
    

	// Writes the delta stoich to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the delta stoich data from a binary data stream.
    void Deserialize(std::istream &in);

 protected:

    Sprog::Kinetics::Reaction *m_react; 

private:
    // Data.
    std::string m_spName;
    double total_stoich_sp;
    double reac_stoich_sp;
    double prod_stoich_sp;
};
// Inline function definitions.
#include "gpc_delta_stoich_inl.h"
// A typedef for a STL vector of stoich changes.
typedef std::vector<DeltaStoich>  DeltaStoichVector;
// A typedef for a STL vector of pointers to species specific stoich changes.
typedef std::vector<DeltaStoich*>  DeltaStoichPtrVector;

};
};

#endif
