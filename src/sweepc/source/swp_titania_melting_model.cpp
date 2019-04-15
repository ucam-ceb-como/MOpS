/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the SinteringModel class declared in the
    swp_sintering_model.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

#include "swp_titania_melting_model.h"
#include "swp_primary.h"
#include "swp_particle.h"
#include "swp_cell.h"

#include <cmath>
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
MeltingModel::MeltingModel()
	: m_enable(false), m_liquidindex(0)
{
}

// Copy constructor.
MeltingModel::MeltingModel(const Sweep::Processes::MeltingModel &copy)
{
    *this = copy;
}

// Stream-reading constructor.
MeltingModel::MeltingModel(std::istream &in)
{
    Deserialize(in);
}

// Destructor.
MeltingModel::~MeltingModel()
{
    // Nothing special to destruct.
}

// OPERATOR OVERLOADING.

// Assignment operator.
MeltingModel &MeltingModel::operator=(const Sweep::Processes::MeltingModel &rhs)
{
    if (this != &rhs) {		// TODO change this
        m_enable = rhs.m_enable;
		m_liquidindex = rhs.m_liquidindex;

		// Copy phases.
		for (PhasePtrVector::const_iterator i = rhs.m_phases.begin();i != rhs.m_phases.end(); ++i) {
			m_phases.push_back((*i)->Clone());
		}
    }
    return *this;
}


// STRUCTURE DEFINITIONS

// Default constructor.
MeltingModel::PHASE::PHASE(void)
{
	name = "";
	T_bulk = 0.0;
	A = 0.0;
};

// Initialising constructor.
MeltingModel::PHASE::PHASE(std::string aname, double aA, double aT, fvector adcomp)
{
	name = aname;
	T_bulk = aT;
	A = aA;
	dcomp = adcomp;
};

// Copy constructor.
MeltingModel::PHASE::PHASE(const PHASE &copy)
{
	name = copy.name;
	T_bulk = copy.T_bulk;
	A = copy.A;
	dcomp = copy.dcomp;
};

// Clone
MeltingModel::PHASE* const MeltingModel::PHASE::Clone(void) const
{
	return new PHASE(*this);
}

// Returns melting point for given particle size
double MeltingModel::PHASE::MeltingPoint(double d){
	double T = 0.0;
	T = std::max(T_bulk*(1 - A / d), 0.0); //lower bound T=0;
	return T;
}


// MODEL ENABLE/DISABLE.

// Returns true is the model is enabled, otherwise false.
bool MeltingModel::IsEnabled(void) const { return m_enable; }

// Enables this sintering model.
void MeltingModel::Enable(void) { m_enable = true; }

// Disables this sintering model.
void MeltingModel::Disable(void) { m_enable = false; }

// Create a melting point transformation 
void MeltingModel::AddPhase(std::string name, double A, double T, fvector dcomp){
	PHASE *new_phase = NULL;
	new_phase = new PHASE(name, A, T, dcomp);
	m_phases.push_back(new_phase);
	return;
}


// PHASE FORMATION
void MeltingModel::CompositionChange(const Cell &sys, const AggModels::Primary &p, fvector &dcomp)
{
	//primary size
	double dp = 6.0 * p.Volume() / p.SurfaceArea();
	double T = 0.0;
	double Tmelt = 0.0;

	//loop over phases
	PhasePtrVector::const_iterator i;
	for (i = m_phases.begin(); i != m_phases.end(); ++i) {
		//get melting point
		Tmelt = (*i)->MeltingPoint(dp);
		//return composition change for highest melting point
		if (Tmelt > T){
			T = Tmelt;
			dcomp = (*i)->dcomp;
		}
	}
}

//return higher melting point
double MeltingModel::MeltingTemp(const AggModels::Primary &p) const
{
	//primary size
    double dp = 6.0 * p.Volume() / p.SurfaceArea();
	double T = 0.0;
	double Tmelt = 0.0;

	//loop over phases
	PhasePtrVector::const_iterator i;
	for (i = m_phases.begin(); i != m_phases.end(); ++i) {
		//get melting point
		Tmelt = (*i)->MeltingPoint(dp);
		//return composition change for highest melting point
		if (Tmelt > T){
			T = Tmelt;
		}
	}

	//return highest melting temperature
	return T;
}

bool MeltingModel::IsLiquid(const Cell &sys, const AggModels::Primary &p) const{
	
	bool liquid = false;

	double Tmelt = MeltingTemp(p);

	if (sys.GasPhase().Temperature() > Tmelt){ liquid = true; };

	return liquid;
}

//return index of liquid component
void MeltingModel::SetLiquidIndex(unsigned int index) {
	m_liquidindex = index;
}

//return index of liquid component
unsigned int MeltingModel::GetLiquidIndex(void) const{
	return m_liquidindex;
}

// READ/WRITE/COPY.

// Creates a copy of the sintering model.
MeltingModel *const MeltingModel::Clone(void) const
{
	return new MeltingModel(*this);
}

// Writes the object to a binary stream.
void MeltingModel::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;
	unsigned int n=0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write if enabled or disabled model.
        if (m_enable) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write liquid component index.
		unsigned int val = (unsigned int)m_liquidindex;
        out.write((char*)&val, sizeof(val));

		//write number of phases
		unsigned int n = (unsigned int) m_phases.size();
		out.write((char*)&n, sizeof(n));

		// Write phases
		for (PhasePtrVector::const_iterator i = m_phases.begin(); i != m_phases.end(); ++i) {
			
			// Write the length of the name to the stream.
			unsigned int n = (*i)->name.length();
			out.write((char*)&n, sizeof(n));

			// Write the name to the stream.
			out.write((*i)->name.c_str(), n);
			
			//write parameter
			double aA = (double)(*i)->A;
			out.write((char*)&aA, sizeof(aA));
			
			//write bulk temperature
			double aT_bulk = (double)(*i)->T_bulk;
			out.write((char*)&aT_bulk, sizeof(aT_bulk));

			// Write size of dcomp
			n = (unsigned int) (*i)->dcomp.size();
			out.write((char*)&n, sizeof(n));
			
			// write dcomp
			for (unsigned int j = 0; j != n; j++){
				double adcomp = (double)(*i)->dcomp[j];
				out.write((char*)&adcomp, sizeof(adcomp));
			}			
			
		}

	} else {
		throw std::invalid_argument("Output stream not ready "
                               "(Sweep, SinteringModel::Serialize).");
    }
}

// Reads the object from a binary stream.
void MeltingModel::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
		unsigned int m = 0;
        double val=0.0;
		std::string str="";
		char *name = NULL;

        switch (version) {
            case 0:
                // Read if enabled
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_enable = (n==1);

                // Read liquid component index.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_liquidindex = (unsigned int)val;

				// Read number of phases
				in.read(reinterpret_cast<char*>(&n), sizeof(n));

				// Read phases.
				for (unsigned int i = 0; i != n; ++i) {
					
					//read size of name
					in.read(reinterpret_cast<char*>(&m), sizeof(m));

					// Read the name.
					name = new char[m];
					in.read(name, m);
					
					//read parameter
					double aA;
					in.read(reinterpret_cast<char*>(&aA), sizeof(aA));

					//read bulk temperature
					double aT;
					in.read(reinterpret_cast<char*>(&aT), sizeof(aT));

					//read size of dcomp
					unsigned int an;
					in.read(reinterpret_cast<char*>(&an), sizeof(an));

					//read dcomp
					fvector dcomp;
					for (unsigned int j = 0; j != an; ++j) {
						in.read(reinterpret_cast<char*>(&val), sizeof(val));
						dcomp.push_back(val);
					}

					//create new phase
					PHASE new_phase(name, aA, aT, dcomp);
					m_phases.push_back(&new_phase);

					delete[] name;
				}

                break;
            default:
                throw std::runtime_error("Serialized version number is invalid "
                                    "(Sweep, SinteringModel::Deserialize).");
        }
    } else {
		throw std::invalid_argument("Input stream not ready "
                               "(Sweep, SinteringModel::Deserialize).");
    }
}
