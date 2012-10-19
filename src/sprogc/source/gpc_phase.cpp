/*
 Author(s):      Martin Martin (mm864)
  Project:        sprog (gas-phase and surface chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2012 Martin Martin.

  File purpose:
    Implementation of the Species class declared in the
    gpc_species.h header file.

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

#include "gpc_phase.h"
#include "gpc_species_comp.h"
#include "gpc_mech.h"
#include "gpc_params.h"
#include <stdexcept>
#include <string>
#include <cmath>
#include "string_functions.h"

using namespace Sprog;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Phase::Phase(void)
    : m_name("")
    , m_id("")
    , m_siteden (0.0)
    , m_mech(NULL)
{}

// Stream-reading constructor.
Phase::Phase(istream &in)
{
    Deserialize(in);
}

// Copy constructor.
Phase::Phase(const Sprog::Phase &phase)
    : m_name("")
    , m_id("")
    , m_siteden (0.0)
    , m_mech(NULL)
{
    *this  = phase;
}

// Destructor.
Phase::~Phase(void)
{
   // There is nothing special to destruct here.
}


// OPERATOR OVERLOADING.

// Assignment operator.
Phase &Phase::operator=(const Sprog::Phase &phase)
{
    // Remember to check for self-assignment!
    if (this!=&phase) {
        // Copy simple member data.
        m_name  = phase.m_name;
        m_id    = phase.m_id;
	m_siteden = phase.m_siteden;
        m_mech  = phase.m_mech;

        // Copy composition.
        m_spcomp.assign(phase.m_spcomp.begin(), phase.m_spcomp.end());
    }
    return *this;
}

// Comparison operator:  Compares two Phase objects.  Returns true
// if they have the same name.
bool Phase::operator==(const Sprog::Phase &phase) const
{
    return m_name.compare(phase.m_name)==0;
}

// Comparison operator:  Compares a Phase object to a string.  Returns
// true if the phase name is the same as the string.
bool Phase::operator==(const std::string &name) const
{
    return m_name.compare(name)==0;
}

// Inequality operator:  Returns false of the two Phase objects have the
// same name.
bool Phase::operator!=(const Sprog::Phase &phase) const
{
    return !(*this==phase);
}

// Inequality operator:  Returns false if the phase name matches
// the string.
bool Phase::operator!=(const std::string &name) const
{
    return !(*this==name);
}


// PHASE NAME & ID.

// Sets the phase name.
void Phase::SetName(const std::string &name)
{
    // If this phase is part of a mechanism then we must check
    // if a phase with this name is already defined.
    if (m_mech != NULL) {
        if (m_mech->FindPhase(name) >= 0) {
            // Oh dear:  Phase with this name already defined.  We'd
            // better throw an error.
            throw invalid_argument("Cannot have two phase named " + name +
                                   "in the same mechanism (Sprog, Phase::SetName).");
        }
    }

    // It is ok to set the name now.
    m_name = name;
}


// Sets the phase site density.
void Phase::SetSiteDensity(const double site_density, const std::string &name)
{
    // If this phase is part of a mechanism then we must check
    // if a phase with this name is already defined.
    if (m_mech != NULL) {
        if ((name.compare("gas") != 0) && ((m_mech->FindSiteDensity(name)) == site_density)) {
            // If the species is gas phase, it won't have a phase name ("") but still have the phase id "g"
			// Oh dear:  Phase with this name already defined.  We'd
            // better throw an error.
            throw invalid_argument("Cannot have two phase site density (Sprog, Phase::SetSiteDensity).");
        }
    }

    // It is ok to set the site density now.
    m_siteden = site_density;
}


// Sets the phase ID: gas, surface ONLY. 
void Phase::SetID(const double site_density, const std::string &name)
{
    // If this phase is part of a mechanism then we must check
    // if a phase with this name is already defined.
    if (m_mech != NULL) {
      if ((m_mech->FindID(name)).length() > 0) {
            // Oh dear:  Phase with this name already defined.  We'd
            // better throw an error.
            throw invalid_argument("This phase name have already got an ID in the same mechanism (Sprog, Phase::SetID).");
        }
    }

    // It is ok to set the id now.
    if (site_density > 0.0){
       m_id = "s";
      }
	  else { m_id = "g";}
   
}


// Sets the phase ID: gas, surface or bulk.
void Phase::SetID(const int site_occupancy, const double site_density, const std::string &name)
{
    // If this phase is part of a mechanism then we must check
    // if a phase with this name is already defined.
    if (m_mech != NULL) {
      if ((m_mech->FindID(name)).length() > 0) {
            // Oh dear:  Phase with this name already defined.  We'd
            // better throw an error.
            throw invalid_argument("This phase have an ID in the same mechanism (Sprog, Phase::SetID).");
        }
    }

    // It is ok to set the id now.
    if (site_density > 0.0){
      if (site_occupancy > 0){ m_id = "s";}
      else { m_id = "b";}
    }
    else{ m_id = "g";
    }
}


// PHASE COMPOSITION.

// Adds a species to the phase composition vector.
void Phase::AddSpecies(const Sprog::SpComp &spcomp)
{
    if (m_mech != NULL) {
        if ((spcomp.Index() >= 0) && ((unsigned)spcomp.Index() < m_mech->SpeciesCount())) {
            // Must check if this species is already defined for this
            // phase.
            SpCompVector::iterator i;
            bool found = false;
            for (i=m_spcomp.begin(); i!=m_spcomp.end(); i++) {
                if ((*i) == spcomp) {
                    // Found species:  Append value.
                    (*i) += spcomp.Count();
                    found = true;
                    break;
                }
            }

            // Species not currently defined for species, so add it.
            if (!found) m_spcomp.push_back(spcomp);

        } else {
            // Species index in spcomp is out-of-range.
            throw out_of_range("Species index is out of range "
                               "(Sprog, Phase::AddSpecies).");
        }
    } else {
        // Parent mechanism has not been assigned yet.
        throw logic_error("Attempt to add species before "
                          "mechanism is assigned (Sprog, Phase::AddSpecies).");
    }
}

// Adds a species to the phase given the index and count.
void Phase::AddSpecies(unsigned int i, unsigned int n)
{
    if (m_mech != NULL) {
        if (i < m_mech->SpeciesCount()) {
            // Must check if this species is already defined for this
            // phase.
            SpCompVector::iterator j;
            bool found = false;
            for (j=m_spcomp.begin(); j!=m_spcomp.end(); j++) {
                if ((unsigned)(*j).Index() == i) {
                    // Found species:  Append value.
                    (*j) += n;
                    found = true;
                    break;
                }
            }

            // Species not currently defined for the phase, so add it.
            if (!found) m_spcomp.push_back(SpComp(i, n));

           
        } else {
            // Species index is out-of-range.
            throw out_of_range("Species index is out of range "
                               "(Sprog, Phase::AddSpecies).");
        }
    } else {
        // Mechanism has not been assigned yet.
        throw logic_error("Attempt to add species before parent "
                          "mechanism is assigned (Sprog, Phase::AddSpecies).");
    }
}

// Adds a species to the phase given the species name and count.
void Phase::AddSpecies(const std::string &name, unsigned int n)
{
    if (m_mech != NULL) {
        // Must find the index of the species in the list.
        int i = m_mech->FindSpecies(name);

        if (i >= 0) {
            // Found species in the list!

            // Must check if this species is already defined for this
            // phase.
            SpCompVector::iterator j;
            bool found = false;
            for (j=m_spcomp.begin(); j!=m_spcomp.end(); j++) {
                if ((*j).Index() == i) {
                    // Found species:  Append value.
                    (*j) += n;
                    found = true;
                    break;
                }
            }

            //Species not currently defined for phase, so add it.
            if (!found) m_spcomp.push_back(SpComp(i, n));

        } else {
            // We have got here because the species wasn't found in the list.
            throw invalid_argument(string(name).append(" not found in "
                                   "species list (Sprog, Phase::AddSpecies)."));
        }
    } else {
        // Species vector has not been assigned yet.
        throw logic_error("Attempt to add species before parent "
                          "mechanism is assigned (Sprog, Phase::AddSpecies).");
    }
}

// Returns true if the phase contains the species (given by index).
bool Phase::ContainsSpecies(unsigned int i) const
{
    // Loop over composition to find species.
    SpCompVector::const_iterator sp;
    for (sp=m_spcomp.begin(); sp!=m_spcomp.end(); sp++) {
        if ((unsigned)sp->Index() == i) {
            return true;
        }
    }

    // We have arrived here because the species wasn't found.
    return false;
}

// Returns true if the phase contains the species (given by reference).
bool Phase::ContainsSpecies(const Sprog::Species &sp) const
{
    if (m_mech != NULL) {
        // Loop over the species to find the index.
        unsigned int i;
        for (i=0; i<m_mech->SpeciesCount(); i++) {
            if (*m_mech->Species(i) == sp) {
                return ContainsSpecies(i);
            }
        }

        // We have got here because the species isn't in the list, therefore
        // the phase cannot contain it.
        return false;
    } else {
        // The parent mechanism is unassigned.
        throw logic_error("Parent mechanism is unassigned "
                          "(Sprog, Phase::ContainsSpecies).");
    }
}

// Returns true if the phase contains the species (given by name).
bool Phase::ContainsSpecies(const std::string &name) const
{
    if (m_mech != NULL) {
        // Loop over the species to find the index.
        unsigned int i;
        for (i=0; i<m_mech->SpeciesCount(); i++) {
            if (*m_mech->Species(i) == name) {
                return ContainsSpecies(i);
            }
        }

        // We have got here because the species isn't in the list, therefore
        // the phase cannot contain it.
        return false;
    } else {
        // The species vector is unassigned.
        throw logic_error("Parent mechanism is unassigned "
                          "(Phase::ContainsSpecies).");
    }
}



// PARENT MECHANISM.

// Sets the parent mechanism.
void Phase::SetMechanism(Sprog::Mechanism &mech)
{
    m_mech = &mech;
}

// Prints a diagnostic output file containing all the
// phase data.  This is used to debug.
void Phase::WriteDiagnostics(std::ostream &out) const
{
    string data = "";
    double val = 0.0;
    int ival = 0;

    if (out.good()) {
        // Name.
        data = m_name + " ";
        out.write(data.c_str(), data.length());
	// ID. 
	data = m_id + " ";
	out.write(data.c_str(), data.length());
        // Composition.
        for (unsigned int j=0; j!=m_spcomp.size(); ++j) {
            // Species name.
            data = m_mech->Species(m_spcomp[j].Index())->Name() + " ";
            out.write(data.c_str(), data.length());
            // Species count.
            ival = m_spcomp[j].Count();
            data = cstr(ival) + " ";
            out.write(data.c_str(), data.length());
        }

       // Site density 
        val = m_siteden;
        data = cstr(val) + " ";
        out.write(data.c_str(), data.length());

    }
}




// Writes the phase to a binary data stream.
void Phase::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Write the serialisation version number to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write the length of the phase name to the stream.
        unsigned int n = m_name.length();
        out.write((char*)&n, sizeof(n));

        // Write the phase name to the stream.
        out.write(m_name.c_str(), n);

		
		// Write the length of the phase id to the stream.
        n = m_id.length();
        out.write((char*)&n, sizeof(n));

		// Write the phase id to the stream.
        out.write(m_id.c_str(), n);
	

        // Write the length of the composition vector to the stream.
        n = m_spcomp.size();
        out.write((char*)&n, sizeof(n));

        // Write the species composition to the stream.
        for (SpCompVector::const_iterator i=m_spcomp.begin(); i!=m_spcomp.end(); i++) {
            // Write the species index.
            int ix = (*i).Index();
            out.write((char*)&ix, sizeof(ix));

            // Write the species count.
            n = (*i).Count();
            out.write((char*)&n, sizeof(n));
        }

        // Write the site density to the stream. - Added by mm864
        double den = (double)m_siteden;
        out.write((char*)&den, sizeof(den));


    } else {
        throw invalid_argument("Output stream not ready (Sprog, Phase::Serialize).");
    }
}

// Reads the phase data from a binary data stream.
void Phase::Deserialize(std::istream &in)
{
    // Clear the phase of all current data.
    m_name  = "";
	m_id = "";
    m_mech  = NULL; 
    m_siteden = 0.0;
    m_spcomp.clear();
   

    if (in.good()) {
        // Read the serialized phase version number.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0; // Need for reading name length.
        char *name = NULL;
        double den= 0.0;

        switch (version) {
            case 0:
                // Read the length of the phase name.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read the species name.
                name = new char[n];
                in.read(name, n);
                m_name.assign(name, n);
                delete [] name;

		// Read the length of the phase id.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read the phase id.
                name = new char[n];
                in.read(name, n);
                m_id.assign(name, n);
                delete [] name;

                // Read length of composition vector and reserve memory.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_spcomp.reserve(n);

                // Read the species composition.
                for (unsigned int i=0; i<n; i++) {
                    // Read the species index.
                    int ix = -1;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Read the species count.
                    unsigned int m = 0;
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));

                    // Add new ElComp object to vector.
                    m_spcomp.push_back(SpComp(ix, m));
                }


		// Read the species site density. - Added by mm864
                in.read(reinterpret_cast<char*>(&den), sizeof(den));
				m_siteden = (double)den;
                break;
            default:
                throw runtime_error("Phase serialized version number "
                                    "is unsupported (Sprog, Phase::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Sprog, Phase::Deserialize).");
    }
}
