/*
  Author(s):      Martin Martin (mm864)
  Project:        sprog (gas-phase and surface chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2012 Martin Martin.

  File purpose:
    This files contains the definition of a structure for chemical reactions.  Also
    contains typedefs and other data structures related to chemical reactions.

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

#ifndef GPC_PHASE_H
#define GPC_PHASE_H
#include "gpc_species.h"
#include "gpc_species_comp.h"
#include "gpc_thermo_params.h"
#include <vector>
#include <string>
#include <iostream>

namespace Sprog
{
class Mechanism; // Forward declaration of mechanism.

// A typedef for a STL vector of pointers to phase.
class Phase;
typedef std::vector<Phase*> PhasePtrVector;

class Phase 
{
   

 public:
   // Constructors.
    Phase(void);              // Default constructor.
    Phase(const Phase &phase); // Copy constructor.
    /* DON'T REQUIRE 
    Phase(  
	  const std::string &name,  // - Element name.
	  const std::string &id,
	  const double siteden          // - Molecular weight.
        );
    */


    Phase(std::istream &in);

    // Destructor
    ~Phase(void);

     // Operator overloads.
    Phase &operator=(const Phase &phase);
    bool operator==(const Phase &phase) const;
    bool operator==(const std::string &name) const;
    bool operator!=(const Phase &phase) const;
    bool operator!=(const std::string &name) const;


   // PHASE NAME and ID.

    // Returns the phase name. (inline)
    const std::string &Name(void) const;

    //Returns the phase ID. (inline)
    const std::string &ID(void) const;

    // Sets the phase name.
    void SetName(const std::string &name);

	// Sets the phase ID: gas, and surface ONLY. 
    void SetID(const double site_density, const std::string &name);
	
    // Sets the phase ID: gas, surface or bulk.
    void SetID(const int site_occupancy, const double site_density, const std::string &name);
    
      // PHASE COMPOSITION.

    // Returns the species composition of the phase (inline).
    const SpCompVector &Composition(void) const;

    // Returns the number of species in the phase (inline).
    unsigned int ComponentCount(void) const;

    // Adds a species to the phase composition using a SpComp object.
    void AddSpecies(const SpComp &spcomp);

    // Adds a species given the index and count.
    void AddSpecies(unsigned int i, unsigned int n);

    // Adds a species given the name.  Species found using parent mechanism.
    void AddSpecies(const std::string &name, unsigned int n);

    // Returns true if phase contains the species (given by index).
    bool ContainsSpecies(unsigned int i) const;

    // Returns true if phase contains the element (given by name).
    bool ContainsSpecies(const std::string &name) const;

    // Returns true if phase contains the element (given by object).
    bool ContainsSpecies(const Species &sp) const;


    // SITE DENSITY.

    // Returns the phase site density () (inline).
    const double SiteDen(void) const;

	// Sets the phase site density. 
    void SetSiteDensity(const double site_density, const std::string &name);

    // PARENT MECHANISM.

    // Returns pointer to parent mechanism. (inline)
    const Sprog::Mechanism *const Mechanism(void) const;

    // Sets the parent mechanism.
    void SetMechanism(Sprog::Mechanism &mech);


    // PHASE LOOKUP. (inline)

    // Returns the index of a phase with the given name in the list.
    // Returns -1 if the phase is not in the list.
    static int Find(
        const std::string &name,     // Name of phase to find.
        const PhasePtrVector &list // List of phases to search.
        );


    // READ/WRITE/COPY FUNCTIONS. 

    // Creates a copy of the phase object. (inline)
    Phase *const Clone(void) const;

    // Writes the species to a binary data stream.
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */)
    {
        ar.template register_type<Sprog::SpComp>();
        // ar.register_type(static_cast< ::IO::Transport * >(NULL));
        ar & m_name & m_spcomp & m_id & m_siteden & m_mech;
    }

    // Writes the species to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the species data from a binary data stream.
    void Deserialize(std::istream &in);

    // Prints a diagnostic output file containing all the
    // phase data.  This is used to debug.
    void WriteDiagnostics(std::ostream &out) const;



protected:
    // Phase data.
    std::string m_name;       // Phase name/symbol.
    std::string m_id;         // Phase id (g or s or b).
    double m_siteden;           // Site density (mole/m^2).
    Sprog::Mechanism *m_mech; // Parent mechanism.
    SpCompVector m_spcomp;    // Species composition.
    
};
// Inline function definitions.
#include "gpc_phase_inl.h"

// A typedef for a STL vector of phase.
typedef std::vector<Phase> PhaseVector;

// A typedef for a STL vector of pointers to phase.
typedef std::vector<Phase*> PhasePtrVector;

};

#endif
