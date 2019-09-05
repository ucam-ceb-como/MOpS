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
	: m_enable(false), m_crossover(false)
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
		m_crossover = rhs.m_crossover;

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
	dmin = 0.0;
	dmax = 0.0;
	id = 0;
};

// Initialising constructor.
MeltingModel::PHASE::PHASE(std::string aname, double aA, double aT, double admin,
	double admax, std::vector<fvector> adcomp, unsigned int iid)
{
	name = aname;
	T_bulk = aT;
	A = aA;
	dcomp = adcomp;
	dmin = admin;
	dmax = admax;
	id = iid;
};

// Copy constructor.
MeltingModel::PHASE::PHASE(const PHASE &copy)
{
	name = copy.name;
	T_bulk = copy.T_bulk;
	A = copy.A;
	dcomp = copy.dcomp;
	dmin = copy.dmin;
	dmax = copy.dmax;
	id = copy.id;
};

// Clone
MeltingModel::PHASE* const MeltingModel::PHASE::Clone(void) const
{
	return new PHASE(*this);
}

// Returns melting point for given particle size
double MeltingModel::PHASE::MeltingPoint(double d){
	double T = 0.0;
	if (d < dmax && d >= dmin){
		T = std::max(T_bulk*(1 - A / d), 0.0); //lower bound T=0;
	}
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
void MeltingModel::AddPhase(std::string name, double A, double T,
	double dmin, double dmax, std::vector<fvector> ddcomp, unsigned int iid){
	
	PHASE *new_phase = NULL;

	if (!m_crossover){	// If not using a fixed cross size
		dmin = 0.0;	
		dmax = 1.0e30;	// some large number
	}

	//add phase
	new_phase = new PHASE(name, A, T, dmin, dmax, ddcomp, iid);
	m_phases.push_back(new_phase);

	//create reverse melting process
	for (std::vector<fvector>::const_iterator it = ddcomp.begin(); it != ddcomp.end(); ++it){
		fvector dcomp = *it;
		fvector dmelt(dcomp.size(),0.0);

		for (unsigned int j = 0; j != dcomp.size(); j++){
			if (dcomp[j] != 0.0){ dmelt[j] = -1.0*dcomp[j]; }
		}
		m_dcompmelt.push_back(dmelt);
	}

	return;
}

// Return the composition change for melting
void MeltingModel::MeltingCompositionChange(std::vector<fvector> &dcomp) const{
	
	dcomp = m_dcompmelt;
	return;
}

// Return component transformation for a particular phase
void MeltingModel::CompositionChange(unsigned int iid, std::vector<fvector> &dcomp){

	//loop over phases
	PhasePtrVector::const_iterator i;
	for (i = m_phases.begin(); i != m_phases.end(); ++i) {
		//compare phase id
		if ((*i)->id == iid){
			// return phase transformation
			dcomp = (*i)->dcomp;
		}
	}
	return;
}

// PHASE FORMATION
void MeltingModel::CompositionChange(const Cell &sys, const AggModels::Primary &p, std::vector<fvector> &dcomp)
{
	//primary size
	double dp = 6.0 * p.Volume() / p.SurfaceArea();
	double T = 0.0;
	double Tmelt = 0.0;

	//loop over phases
	PhasePtrVector::const_iterator i;
	if (m_crossover){
		// crossover size specified
		for (i = m_phases.begin(); i != m_phases.end(); ++i) {
			if (dp > (*i)->dmin && dp <= (*i)->dmax)
			{
				dcomp = (*i)->dcomp;
			}
		}
	}else{
		//melting point dependent crossover size
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
	return;
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

//enable fixed crossover diameter
void MeltingModel::EnableCrossover(){
	m_crossover = true;
}

//return if fixed crossover diameter is enabled
bool MeltingModel::IsEnableCrossover() const{
	return m_crossover;
}

// READ/WRITE/COPY.

// Creates a copy of the melting model.
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
	unsigned int nn = 0;

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

		// Write if fixed crossover is enabled
		if (m_crossover) {
			out.write((char*)&trueval, sizeof(trueval));
		}
		else {
			out.write((char*)&falseval, sizeof(falseval));
		}

		//writing melting transformations
		unsigned int n = (unsigned int)m_dcompmelt.size();
		out.write((char*)&n, sizeof(n));

		for (unsigned int j = 0; j != n; j++){

			nn = (unsigned int)m_dcompmelt[j].size();
			out.write((char*)&n, sizeof(nn));

			for (unsigned int k = 0; k != nn; k++){
				double adcomp = (double)m_dcompmelt[j][k];
				out.write((char*)&adcomp, sizeof(adcomp));
			}
		}


		// Write number of phases
		n = (unsigned int) m_phases.size();
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

			//write dmin
			double admin = (double)(*i)->dmin;
			out.write((char*)&admin, sizeof(admin));

			//write dmax
			double admax = (double)(*i)->dmax;
			out.write((char*)&admax, sizeof(admax));

			//write phase id
			unsigned int iid = (unsigned int)(*i)->id;
			out.write((char*)&iid, sizeof(iid));

			// Write size of dcomp
			n = (unsigned int) (*i)->dcomp.size();
			out.write((char*)&n, sizeof(n));
			
			for (unsigned int j = 0; j != n; j++){
				// Write size of composition in dcomp
				nn = (unsigned int)(*i)->dcomp[j].size();
				out.write((char*)&n, sizeof(nn));

				// Write composition in dcomp
				for (unsigned int k = 0; k != nn; k++){
					double adcomp = (double)(*i)->dcomp[j][k];
					out.write((char*)&adcomp, sizeof(adcomp));
				}
				
			}		
			
		}

	} else {
		throw std::invalid_argument("Output stream not ready "
                               "(Sweep, MeltingModel::Serialize).");
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

				// Read if fixed crossover is enabled
				in.read(reinterpret_cast<char*>(&n), sizeof(n));
				m_crossover = (n == 1);

				//Read number of melting transformations
				unsigned int an;
				in.read(reinterpret_cast<char*>(&an), sizeof(an));

				//Read melting transformations
				for (unsigned int j = 0; j != an; ++j) {

					//read size of composition in dmelt
					unsigned int ann;
					in.read(reinterpret_cast<char*>(&ann), sizeof(ann));

					fvector ddmelt;

					//read composition in dcomp
					for (unsigned int k = 0; k != ann; ++k) {
						in.read(reinterpret_cast<char*>(&val), sizeof(val));
						ddmelt.push_back(val);
					}

					m_dcompmelt.push_back(ddmelt);
				}

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

					//read dmin
					double admin;
					in.read(reinterpret_cast<char*>(&admin), sizeof(admin));

					//read dmax
					double admax;
					in.read(reinterpret_cast<char*>(&admax), sizeof(admax));

					//read phase is
					unsigned int iid;
					in.read(reinterpret_cast<char*>(&iid), sizeof(iid));

					//read size of dcomp
					unsigned int an;
					in.read(reinterpret_cast<char*>(&an), sizeof(an));

					//read dcomp
					std::vector<fvector> dcomp;
					for (unsigned int j = 0; j != an; ++j) {
						
						//read size of composition in dcomp
						unsigned int ann;
						in.read(reinterpret_cast<char*>(&ann), sizeof(ann));
						
						fvector ddcomp;

						//read composition in dcomp
						for (unsigned int k = 0; k != ann; ++k) {
							in.read(reinterpret_cast<char*>(&val), sizeof(val));
							ddcomp.push_back(val);
						}

						dcomp.push_back(ddcomp);
					}

					//create new phase
					PHASE new_phase(name, aA, aT, admin, admax, dcomp, iid);
					m_phases.push_back(&new_phase);

					delete[] name;
				}

                break;
            default:
                throw std::runtime_error("Serialized version number is invalid "
                                    "(Sweep, MeltingModel::Deserialize).");
        }
    } else {
		throw std::invalid_argument("Input stream not ready "
                               "(Sweep, MeltingModel::Deserialize).");
    }
}
