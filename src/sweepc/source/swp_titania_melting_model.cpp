/*
	Author(s):      Casper Lindberg (csl37)
	Project:        sweep (population balance solver)

	Copyright (C) 2019 Casper Lindberg.

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
	Prof. Markus Kraft
	Dept of Chemical Engineering
	University of Cambridge
	Philippa Fawcett Drive
	Cambridge
	CB3 0AS
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

//! Default constructor.
MeltingModel::MeltingModel()
	: m_enable(false)
{
	m_crossover = NULL;
}

//! Copy constructor.
MeltingModel::MeltingModel(const Sweep::Processes::MeltingModel &copy)
{
    *this = copy;
}

//! Stream-reading constructor.
MeltingModel::MeltingModel(std::istream &in)
{
    Deserialize(in);
}

//! Destructor.
MeltingModel::~MeltingModel()
{
	// delete the crossover
	delete m_crossover;
	m_crossover = NULL;
	// delete the phase transformations
	for (PhasePtrVector::const_iterator i = m_phases.begin(); i != m_phases.end(); ++i){
		delete *i;
	}
	m_phases.clear();
}

// OPERATOR OVERLOADING.

//! Assignment operator.
MeltingModel &MeltingModel::operator=(const Sweep::Processes::MeltingModel &rhs)
{
    if (this != &rhs) {	
        m_enable = rhs.m_enable;
		m_type = rhs.m_type;
		m_OIndex = rhs.m_OIndex;
		m_TiIndex = rhs.m_TiIndex;
		m_dcompmelt = rhs.m_dcompmelt;
		
		// Copy phase transformations
		for (PhasePtrVector::const_iterator i = rhs.m_phases.begin();i != rhs.m_phases.end(); ++i) {
			m_phases.push_back((*i)->Clone());
		}

		// Copy crossover
		m_crossover = rhs.m_crossover->Clone();		
    }
    return *this;
}


// PHASE STRUCTURE DEFINITIONS

//! Default constructor.
MeltingModel::PHASE::PHASE(void)
{
	name = "";
	T_bulk = 0.0;
	A = 0.0;
	id = 0;
};

//! Initialising constructor.
MeltingModel::PHASE::PHASE(std::string aname, double aA, double aT, std::vector<fvector> adcomp, unsigned int iid)
{
	name = aname;
	T_bulk = aT;
	A = aA;
	dcomp = adcomp;
	id = iid;
};

//! Copy constructor.
MeltingModel::PHASE::PHASE(const PHASE &copy)
{
	name = copy.name;
	T_bulk = copy.T_bulk;
	A = copy.A;
	dcomp = copy.dcomp;
	id = copy.id;
};

//! Clone
MeltingModel::PHASE* const MeltingModel::PHASE::Clone(void) const
{
	return new PHASE(*this);
}

//! Returns the melting point for given particle size
double MeltingModel::PHASE::MeltingPoint(double d){
	double T = 0.0;

	T = std::max(T_bulk*(1 - A / d), 0.0); //lower bound T=0;

	return T;
}

//! Create a phase transformation 
void MeltingModel::AddPhase(std::string name, double A, double T, std::vector<fvector> ddcomp, unsigned int iid){

	PHASE *new_phase = NULL;

	//create and add new phase transformation
	new_phase = new PHASE(name, A, T, ddcomp, iid);
	m_phases.push_back(new_phase);

	//create reverse tranformation for melting process
	for (std::vector<fvector>::const_iterator it = ddcomp.begin(); it != ddcomp.end(); ++it){
		fvector dcomp = *it;
		fvector dmelt(dcomp.size(), 0.0);

		for (unsigned int j = 0; j != dcomp.size(); j++){
			if (dcomp[j] != 0.0){ dmelt[j] = -1.0*dcomp[j]; }
		}
		m_dcompmelt.push_back(dmelt);
	}

	return;
}

// CROSSOVER STRUCTURE DEFINITIONS

//! Default constructor.
MeltingModel::CROSSOVER::CROSSOVER(void)
{
	name = "";
	phase_l = NULL;
	phase_r = NULL;
	a = 0.0;
	b = 0.0;
	c = 0.0;
	d = 0.0;
};

//! Initialising constructor (Composition model).
MeltingModel::CROSSOVER::CROSSOVER(std::string aname, PHASE* aphase_l, PHASE* aphase_r, 
	double aa, double ab, double ac, double ad)
{
	name = aname;
	phase_l = aphase_l;
	phase_r = aphase_r;
	a = aa;
	b = ab;
	c = ac;
	d = ad;
};

//! Initialising constructor (Gibbs model).
MeltingModel::CROSSOVER::CROSSOVER(std::string aname, PHASE* aphase_l, PHASE* aphase_r, double aa)
{
	name = aname;
	phase_l = aphase_l;
	phase_r = aphase_r;
	a = aa;
	b = 0.0;
	c = 0.0;	
	d = 0.0;
};

//! Copy constructor.
MeltingModel::CROSSOVER::CROSSOVER(const CROSSOVER &copy)
{
	name = copy.name;
	phase_l = copy.phase_l;
	phase_r = copy.phase_r;
	a = copy.a;
	b = copy.b;
	c = copy.c;
	d = copy.d;
};

//! Clone
MeltingModel::CROSSOVER* const MeltingModel::CROSSOVER::Clone(void) const
{
	return new CROSSOVER(*this);
}

//! Returns the crossover diameter for the composition model c4e-238 Eq. (30). 
double MeltingModel::CROSSOVER::CrossoverDiameter(double sigma_o) const
{
	return a + b * erf(c*(sigma_o - d));
}

//! Returns the fixed crossover diameter for the Gibbs model c4e-238 Eq. (30).
double MeltingModel::CROSSOVER::CrossoverDiameter() const
{
	return a;
}

//! Return phase transformation (Composition model)
const MeltingModel::PHASE* const MeltingModel::CROSSOVER::GetTransformation(double dp, double sigma_o) const
{
	if (dp < CrossoverDiameter(sigma_o)){
		return phase_l;
	}
	else{
		return phase_r;
	}
}

//! Return phase transformation (Gibbs model)
const MeltingModel::PHASE* const MeltingModel::CROSSOVER::GetTransformation(double dp) const
{
	if (dp < CrossoverDiameter()){
		return phase_l;
	}
	else{
		return phase_r;
	}
}

//! Create a phase crossover diameter
void MeltingModel::AddCrossover(std::string name, std::string id_l, std::string id_r, double aa, double ab, double ac, double ad){

	CROSSOVER* new_crossover = NULL;
	PHASE* phase_l = NULL;
	PHASE* phase_r = NULL;

	//identify phases
	PhasePtrVector::const_iterator i;
	for (i = m_phases.begin(); i != m_phases.end(); ++i) {
		if (id_l.compare((*i)->name) == 0){
			phase_l = *i;
		}
		else if (id_r.compare((*i)->name) == 0){
			phase_r = *i;
		}
	}

	//create a crossover structure
	if (m_type == Gibbs){
		// for the Gibbs model (constant crossover diameter), the b,c,d parameters are set to 0
		new_crossover = new CROSSOVER(name, phase_l, phase_r, aa);
	}
	else{
		new_crossover = new CROSSOVER(name, phase_l, phase_r, aa, ab, ac, ad);
	}
	//add crossover
//	m_crossovers.push_back(new_crossover);
	m_crossover = new_crossover;
}

// MODEL ENABLE/DISABLE.

//! Returns true is the model is enabled, otherwise false.
bool MeltingModel::IsEnabled(void) const { return m_enable; }

//! Enables this model.
void MeltingModel::Enable(void) { m_enable = true; }

//! Disables this model.
void MeltingModel::Disable(void) { m_enable = false; }

// MODEL TYPE

//! Returns the model type.
MeltingModel::TransformationType MeltingModel::Type(void) const { return m_type; }

//! Sets the model type.
void MeltingModel::SetType(TransformationType t) { m_type = t; }

// DATA ACCESS

//! Set oxygen component indices
void MeltingModel::SetOxygenComp(std::vector<unsigned int> OIndex)
{
	m_OIndex = OIndex;
}

//! Return oxygen component indices
std::vector<unsigned int> MeltingModel::GetOxygenComp() const
{
	return m_OIndex;
}

//! Set titanium component indices
void MeltingModel::SetTitaniumComp(std::vector<unsigned int> TiIndex)
{
	m_TiIndex = TiIndex;
}

//! Return titanium component indices
std::vector<unsigned int> MeltingModel::GetTitaniumComp() const
{
	return m_TiIndex;
}

//! Return number of phase transformations
unsigned int MeltingModel::PhaseChangeCount(void) const
{
	return m_phases.size();
}

//! Check a phase transformation exists
bool MeltingModel::PhaseChangeValid(std::string id) const
{
	PhasePtrVector::const_iterator i;
	for (i = m_phases.begin(); i != m_phases.end(); ++i) {

		if (id.compare((*i)->name) == 0) return true;

	}

	return false;
}

// PERFORMING THE PROCESS

//! Return O:Ti ratio
double MeltingModel::OTiRatio(const AggModels::Primary &p) const
{
	double nO = 0.0;
	double nTi = 0.0;
	std::vector<unsigned int>::const_iterator it;

	// Sum O and Ti compositions
	for (it = m_OIndex.begin(); it != m_OIndex.end(); ++it)
	{
		nO += p.Composition((*it));
	}
	for (it = m_TiIndex.begin(); it != m_TiIndex.end(); ++it)
	{
		nTi += p.Composition((*it));
	}

	// Return ratio
	return nO / nTi;
}

//! Check if particle is liquid i.e. above transformation temperature
bool MeltingModel::IsLiquid(const Cell &sys, const AggModels::Primary &p) const{

	bool liquid = false;

	double Tmelt = MeltingTemp(p);

	if (sys.GasPhase().Temperature() > Tmelt){ liquid = true; };

	return liquid;
}

//! Return the melting temperature for a particle
//  This is the highest melting temperature
double MeltingModel::MeltingTemp(const AggModels::Primary &p) const
{
	double dp = p.SphDiameter(); //primary size
	double T = 0.0;
	double Tmelt = 0.0;

	//loop over phases
	PhasePtrVector::const_iterator i;
	for (i = m_phases.begin(); i != m_phases.end(); ++i) {
		//get (size dependent) phase transformation temperature
		Tmelt = (*i)->MeltingPoint(dp);
		if (Tmelt > T){
			T = Tmelt;
		}
	}

	//return highest melting temperature
	return T;
}

//! Returns composition change for melting (crystal -> non-crystalline)
void MeltingModel::MeltingCompositionChange(std::vector<fvector> &dcomp) const{
	
	dcomp = m_dcompmelt;
	return;
}

//! Return composition change for given non-crystalline -> crystal phase transformation
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

//! Return composition change for non-crystalline -> crystal transformation 
void MeltingModel::CompositionChange(const Cell &sys, const AggModels::Primary &p, std::vector<fvector> &dcomp)
{
	double dp = p.SphDiameter(); 	//primary size
	double T = 0.0;
	double Tmelt = 0.0;
	PhasePtrVector::const_iterator i;
	
	//Get transformation model
	switch (m_type){
		case Melting: //melting point dependent crossover size
			//loop over phases
			for (i = m_phases.begin(); i != m_phases.end(); ++i) {
				//get melting point
				Tmelt = (*i)->MeltingPoint(dp);
				//return composition change for highest melting point
				if (Tmelt > T){
					T = Tmelt;
					dcomp = (*i)->dcomp;
				}
			}
			break;
		case Gibbs:
			dcomp = m_crossover->GetTransformation(dp)->dcomp;
			break;
		case Composition:
			dcomp = m_crossover->GetTransformation(dp, OTiRatio(p))->dcomp;			
			break;
	}

	return;
}

// READ/WRITE/COPY.

//! Creates a copy of the melting model.
MeltingModel *const MeltingModel::Clone(void) const
{
	return new MeltingModel(*this);
}

//! Writes the object to a binary stream.
void MeltingModel::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;
	unsigned int n=0;
	unsigned int nn = 0;
	double val;

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

		// Write type.
		unsigned int t = (unsigned int)m_type;
		out.write((char*)&t, sizeof(t));

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

		// Write Ti indices
		n = (unsigned int)m_TiIndex.size();
		out.write((char*)&n, sizeof(n));
		for (unsigned int k = 0; k != n; k++){
			unsigned int TiIndex = (unsigned int)m_TiIndex[k];
			out.write((char*)&TiIndex, sizeof(TiIndex));
		}

		// Write O indices
		n = (unsigned int)m_OIndex.size();
		out.write((char*)&n, sizeof(n));
		for (unsigned int k = 0; k != n; k++){
			unsigned int OIndex = (unsigned int)m_OIndex[k];
			out.write((char*)&OIndex, sizeof(OIndex));
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

		// Write crossover
		if (m_type == Gibbs || m_type == Composition){
			// Write the length of the name to the stream.
			unsigned int n = m_crossover->name.length();
			out.write((char*)&n, sizeof(n));

			// Write the name to the stream.
			out.write(m_crossover->name.c_str(), n);

			//write parameters 
			val = (double)m_crossover->a;
			out.write((char*)&val, sizeof(val));

			val = (double)m_crossover->b;
			out.write((char*)&val, sizeof(val));

			val = (double)m_crossover->c;
			out.write((char*)&val, sizeof(val));

			val = (double)m_crossover->d;
			out.write((char*)&val, sizeof(val));

			//pointers to phases: write names
			//write length of name
			n = m_crossover->phase_l->name.length();
			out.write((char*)&n, sizeof(n));
			//write name
			out.write(m_crossover->phase_l->name.c_str(), n);

			//write length of name
			n = m_crossover->phase_r->name.length();
			out.write((char*)&n, sizeof(n));
			//write name
			out.write(m_crossover->phase_r->name.c_str(), n);
		}
	} else {
		throw std::invalid_argument("Output stream not ready "
                               "(Sweep, MeltingModel::Serialize).");
    }
}

//! Reads the object from a binary stream.
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
		unsigned int intval = 0;
		std::string str="";
		char *name = NULL;

        switch (version) {
            case 0:
                // Read if enabled
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_enable = (n==1);

				// Read type.
				in.read(reinterpret_cast<char*>(&n), sizeof(n));
				m_type = (TransformationType)n;

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

				//Read Ti Indices
				in.read(reinterpret_cast<char*>(&an), sizeof(an));
				for (unsigned int j = 0; j != an; ++j) {					
					in.read(reinterpret_cast<char*>(&intval), sizeof(intval));
					m_TiIndex.push_back(intval);
				}

				//Read O Indices
				in.read(reinterpret_cast<char*>(&an), sizeof(an));
				for (unsigned int j = 0; j != an; ++j) {
					in.read(reinterpret_cast<char*>(&intval), sizeof(intval));
					m_OIndex.push_back(intval);
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

					//read phase id
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
					PHASE new_phase(name, aA, aT, dcomp, iid);
					m_phases.push_back(&new_phase);

					delete[] name;
				}

				// read crossover
				if (m_type == Gibbs || m_type == Composition){
					//read size of name
					in.read(reinterpret_cast<char*>(&m), sizeof(m));

					// Read the name.
					char* cross_name = new char[m];
					in.read(cross_name, m);

					//read parameters
					double cross_a;
					in.read(reinterpret_cast<char*>(&cross_a), sizeof(cross_a));

					double cross_b;
					in.read(reinterpret_cast<char*>(&cross_b), sizeof(cross_b));

					double cross_c;
					in.read(reinterpret_cast<char*>(&cross_c), sizeof(cross_c));

					double cross_d;
					in.read(reinterpret_cast<char*>(&cross_d), sizeof(cross_d));

					//read below and above crossover phase names 
					//Read size of name
					in.read(reinterpret_cast<char*>(&m), sizeof(m));
					//Read the name.
					char*  cross_id_l = new char[m];
					in.read(cross_id_l, m);

					//Read size of name
					in.read(reinterpret_cast<char*>(&m), sizeof(m));
					//Read the name.
					char*  cross_id_r = new char[m];
					in.read(cross_id_r, m);

					//Create the crossover
					AddCrossover(cross_name, cross_id_l, cross_id_r, cross_a, cross_b, cross_c, cross_d);

					delete[] cross_name;
					delete[] cross_id_l;
					delete[] cross_id_r;
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
