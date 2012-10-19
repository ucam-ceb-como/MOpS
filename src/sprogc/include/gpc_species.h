/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This file contains the definition of a structure for a chemical species.  File
    also contains typdefs and structure definitions related to Species objects.

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

#ifndef GPC_SPECIES_H
#define GPC_SPECIES_H

#include "gpc_element.h"
#include "gpc_el_comp.h"
#include "gpc_thermo_params.h"
#include "transport.h"
#include <vector>
#include <string>
#include <iostream>

namespace Sprog
{
class Mechanism; // Forward declaration of mechanism.
class Phase; // Forward declaration of phase

// A typedef for a STL vector of pointers to species.
class Species;
typedef std::vector<Species*> SpeciesPtrVector;

class Species
{
public:
    // Constructors.
    Species(void);              // Default constructor.
    Species(const Species &sp); // Copy constructor.
    Species(std::istream &in);

    // Default destructor.
    ~Species(void);

    // Operator overloads.
    Species &operator=(const Species &sp);
    bool operator==(const Species &sp) const;
    bool operator==(const std::string &name) const;
    bool operator!=(const Species &sp) const;
    bool operator!=(const std::string &name) const;


    // SPECIES NAME.

    // Returns the species name. (inline)
    const std::string &Name(void) const;

    // Sets the species name.
    void SetName(const std::string &name);


    // SPECIES SITE OCCUPANCY

    // Return the species site occupancy. (inline)
    const int SiteOccupancy(void) const;

    // Set the species site occupancy. 
    void SetSiteOccupancy(const std::string &name, const int siteOccp); 

	// SPECIES PHASE NAME 
	
	// Returns the species phase name. (inline)
    const std::string &PhaseName(void) const;
	
	// Sets the species phase name.
    void SetPhaseName(const std::string &phaseName, const std::string &name);
	
    // ELEMENTAL COMPOSITION.

    // Returns the elemental composition of the species.
    const ElCompVector &Composition(void) const;

    // Returns the number of elements required to define the species.
    unsigned int ComponentCount(void) const;

    // Returns the total number of atoms in the species.
    unsigned int AtomCount(void) const;

    // Returns the number of the given element in the species.
    unsigned int AtomCount(unsigned int iel) const;

    // Returns the number of the given element (by name) in the species.
    unsigned int AtomCount(std::string name) const;

    // Adds an element to the species composition using an ElComp object.
    void AddElement(const ElComp &elcomp);

    // Adds an element given the index and count.
    void AddElement(unsigned int i, unsigned int n);

    // Adds an element given the name.  Element found using parent mechanism.
    void AddElement(const std::string &name, unsigned int n);

    // Returns true if species contains the element (given by index).
    bool ContainsElement(unsigned int i) const;

    // Returns true if species contains the element (given by name).
    bool ContainsElement(const std::string &name) const;

    // Returns true if species contains the element (given by object).
    bool ContainsElement(const Element &el) const;

    // SURFACE COVERAGE 
    //    void setSurfaceCoverages(const double* theta) const;


    // MOLECULAR WEIGHT.

    // Returns the species molecular weight.
    double MolWt(void) const;

    // Recalculates the molecular weight of the
    // species using the elements in the parent mechanism.
    double CalcMolWt(void);


    // Elements.
    /* Removed as the elements can be defined using the parent mechanism.
    const ElementPtrVector *const Elements(void); // Returns the vector of elements used to define species.
    void SetElements(const ElementPtrVector *const els); // Sets the vector of elements.
    */


    // THERMODYNAMIC FITTING PARAMETERS.

    // Returns the number of thermo parameter ranges.
    unsigned int ThermoRangeCount(void) const;

    // Sets the start temperature for the thermo parameter range.
    void SetThermoStartTemperature(const double T);

    // Returns the set of thermo parameters valid for the given temperature.
    const Sprog::Thermo::THERMO_PARAMS &ThermoParams(const double T) const;

    // Adds a set of thermo parameters with the given end point temperature.
    void AddThermoParams(
        const double T, // Maximum temperature for which parameters are valid.
        const Sprog::Thermo::THERMO_PARAMS &params // Thermo params to add.
        );

    // Removes the parameters from the list valid for the given temperature.
    void RemoveThermoParams(const double T);


    // PARENT MECHANISM.

    // Returns pointer to parent mechanism.
    const Sprog::Mechanism *const Mechanism(void) const;

    // Sets the parent mechanism.
    void SetMechanism(Sprog::Mechanism &mech);


    // SPECIES LOOKUP.

    // Returns the index of a species with the given name in the list.
    // Returns -1 if the species is not in the list.
    static int Find(
        const std::string &name,     // Name of species to find.
        const SpeciesPtrVector &list // List to search.
        );


    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the species object.
    Species *const Clone(void) const;

    // Writes the element to a binary data stream.
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */)
    {
        ar.template register_type<Sprog::ElComp>();
        ar.register_type(static_cast< ::IO::Transport * >(NULL));
        ar & m_name & m_elcomp & m_molwt & m_mech & m_thermoparams & m_T1 & m_phaseName & site_occupancy & m_phase & m_transport; 
    }

    // Writes the element to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the element data from a binary data stream.
    void Deserialize(std::istream &in);

    // Prints a diagnostic output file containing all the
    // species data.  This is used to debug.
    void WriteDiagnostics(std::ostream &out) const;

    //! Writes out the kept species as decided by LOI analysis to the reduced mechanism output file.
    void WriteSpecies(std::ostream &out) const;

	// TRANSPORT RELATED FUNCTIONS Added by vinod

	void setTransportData(const IO::Transport& transport);
	IO::Transport& getTransportData() const;

	//! Viscosity of the pure species in Kg/m-s
	double getViscosity(double T) const;

	//! Self diffusion coefficient m^2/s
	double getSelfDiffusion(double T, double p) const;

	//! Thermal conductivity in J/m-s-K
	double getThermalConductivity(double T, double p, double cp) const;

	//! Collision diameter in m
	double getCollisionDiameter() const;

	//! Check for transport data
	bool hasTransportData() const;

protected:
    // Species data.
    std::string m_name;              // Name/symbol.

    ElCompVector m_elcomp;           // Elemental composition.
    double m_molwt;                    // Molecular weight (kg/mol).
    Sprog::Mechanism *m_mech;        // Parent mechanism.
    
    //double activity;                 // Bulk species activity (NOT YET).  

   
    // Thermo parameters for different temperature ranges.  The map key is the
    // end point temperature up to which the parameters are valid.
    Sprog::Thermo::ThermoMap m_thermoparams;
    double m_T1; // Start temperature for range.
    
    int site_occupancy;              // Site occupancy for surface species (sigma).
   
    std::string m_phaseName; 	     // Phase Name for species
    Sprog::Phase *m_phase;			 // Phase pointer for species
private:
    //! Data for calculating transport properties
    IO::Transport *m_transport;
};

// Inline function definitions.
#include "gpc_species_inl.h"

// A typedef for a STL vector of species.
typedef std::vector<Species> SpeciesVector;

// A typedef for a STL vector of pointers to species.
typedef std::vector<Species*> SpeciesPtrVector;
};

#endif
