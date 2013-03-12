/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the PSR class declared in the
    mops_psr.h header file.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
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

#include "mops_psr.h"
#include "mops_mixture.h"

#include <vector>
#include <cmath>
#include <stdexcept>


using namespace std;

namespace Mops {

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (public, requires mechanism).
PSR::PSR(const Mops::Mechanism &mech)
: Reactor(mech),
  m_restime(1.0),
  m_in(NULL),
  m_out(NULL),
  m_invrt(1.0),
  m_infH(0.0),
  m_default_birth(Sweep::Processes::BirthProcess::iStochastic),
  m_default_death(Sweep::Processes::DeathProcess::iContRescale)
{
    m_infHs.resize(m_neq);
}

// Copy constructor.
PSR::PSR(const Mops::PSR &copy)
: Reactor(copy),
  m_restime(copy.m_restime),
  m_in(NULL),
  m_out(copy.m_out),
  m_invrt(copy.m_invrt),
  m_infH(copy.m_infH),
  m_infHs(copy.m_infHs),
  m_default_birth(copy.m_default_birth),
  m_default_death(copy.m_default_death)
{
    if (copy.m_in) m_in = copy.m_in->Clone();
}

// Stream-reading constructor.
PSR::PSR(std::istream &in, const Mops::Mechanism &mech)
: Reactor(mech)
{
    Deserialize(in, mech);
}

// Default destructor.
PSR::~PSR(void)
{
    if (m_in != NULL) delete m_in;
}

// OPERATORS.

Mops::PSR &PSR::operator=(const Mops::PSR &rhs)
{
    if (this != &rhs) {
        // Copy base class.
        Reactor::operator=(rhs);

        // Residence time.
        m_restime = rhs.m_restime;

        // Inflow stream.
        delete m_in;
        m_in = NULL;
        if (rhs.m_in) m_in = rhs.m_in->Clone();

        // Precalculated terms.
        m_infH = rhs.m_infH;
        m_infHs.assign(rhs.m_infHs.begin(), rhs.m_infHs.end());

        // Birth and death types
        m_default_birth = rhs.m_default_birth;
        m_default_death = rhs.m_default_death;
    }
    return *this;
}

/*!
 * @param os    Output stream
 * @param net   PSR object to print
 * @return      Output stream
 */
std::ostream& operator<<(
        std::ostream &os,
        const Mops::PSR &r) {
    os << "[PSR]," <<
            " ConstP=" << r.IsConstP() <<
            " ConstV=" << r.IsConstV() <<
            "\n";
    os << "filled with: " << *(r.Mixture());
    os << "with inflow: " << *(r.Inflow()->Mixture());
    return os;
}


// RESIDENCE TIME.

// Returns the reactor residence time.
double PSR::ResidenceTime(void) const
{
    return m_restime;
}

// Sets the reactor residence time.
void PSR::SetResidenceTime(double t)
{
    if (t > 0.0) {
        m_restime = t;
        m_invrt = 1.0 / m_restime;
    } else {
        throw out_of_range("Residence time must be greater "
                           "than zero (Mops, PSR::SetResidenceTime).");
    }
}


// INFLOW CONDITIONS.

// Returns the inflow stream
Mops::FlowStream *const PSR::Inflow(void) const {return m_in;}

// Returns the outflow stream
Mops::FlowStream *const PSR::Outflow(void) const {return m_out;}

/*!
 * Intialise the inflow's birth process
 */
void PSR::InitialiseInflow() {

    // Calculate inflow enthalpies
    // TODO: ensure dynamic mixtures have proper enthalpies
    m_infH = m_in->Mixture()->GasPhase().BulkH() /
            (Sprog::R * m_in->Mixture()->GasPhase().Temperature());
    m_in->Mixture()->GasPhase().Hs_RT(m_infHs);

    // Create a new birth process
    Sweep::Processes::BirthProcess bp(Mech()->ParticleMech());
    bp.SetBirthType(m_default_birth);
    bp.SetCell(m_in->Mixture());
    bp.SetA(m_invrt);

    // If an upstream process is a move process, ensure the inflow is turned-off
    if (m_in->Mixture()->OutflowCount() > 0) {
        const Sweep::Processes::DeathPtrVector &outf = m_in->Mixture()->Outflows();
        for (Sweep::Processes::DeathPtrVector::const_iterator it = outf.begin();
                it != outf.end(); ++it) {
            if ((*it)->GetDeathType() == Sweep::Processes::DeathProcess::iContMove
                    || (*it)->GetDeathType() == Sweep::Processes::DeathProcess::iStochMove) {
                bp.SetProcessSwitch(false);
            }
        }
    }

    // Add the process to the Mixture, which clones it and takes ownership
    m_mix->AddInflow(bp);
}

/*!
 * Intialise the outflow's death process
 */
void PSR::InitialiseOutflow() {

    Sweep::Processes::DeathProcess dp(Mech()->ParticleMech());

    // Set the form of the death process.
    dp.SetDeathType(m_default_death);
    dp.SetA(m_invrt);

    // Set the outlet cell if a reactor is present for outflow
    if (m_out != NULL) {
        if (m_out->HasReacOutflow()) dp.SetCell(m_out->Outflow()->Mixture());
    }

    m_mix->AddOutflow(dp);
}

/*!
 * Sets the flowstream which points to the inflow mixture conditions. The
 * reactor takes ownership of the flowstream and therefore responsibility
 * for its deletion.
 *
 * @param inf   Inflow stream
 */
void PSR::SetInflow(Mops::FlowStream &inf)
{
    if (m_in != NULL) delete m_in;
    m_in = inf.Clone();

    // Initialise the inflow
    InitialiseInflow();
}

/*!
 * Sets the reactors out pointer to an outflow stream, then creates a new
 * death process in the reactor's mixture.
 *
 * @param out   Outflow stream
 */
void PSR::SetOutflow(Mops::FlowStream &out) {
    m_out = &out;

    InitialiseOutflow();
}

// Does the reactor have an outflow set?
bool PSR::HasOutflow() const {
    if (m_out!=NULL) return true;
    else return false;
}

// Does the reactor have an inflow set?
bool PSR::HasInflow() const {
    if (m_in!=NULL) return true;
    else return false;
}

// Set the inflow process type
void PSR::SetInflowType(Sweep::Processes::BirthProcess::BirthType btype) {
    m_default_birth = btype;
}

// Set the inflow process type
void PSR::SetOutflowType(Sweep::Processes::DeathProcess::DeathType dtype) {
    m_default_death = dtype;
}

// READ/WRITE/COPY FUNCTIONS.

// Creates a copy of the PSR object.
PSR* PSR::Clone() const
{
    return new PSR(*this);
}

// Writes the PSR to a binary data stream.
void PSR::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;
    
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output base class Reactor.
        this->Reactor::Serialize(out);

        // Output the residence time.
        double val = (double)m_restime;
        out.write((char*)&val, sizeof(val));

        // Output the inflow stream.
        if (m_in != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            m_in->Serialize(out);
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

    } else {
        throw invalid_argument("Output stream not ready (Mops, PSR::Serialize).");
    }
}

// Reads the Reactor data from a binary data stream.
void PSR::Deserialize(std::istream &in, const Mops::Mechanism &mech)
{

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;
        unsigned int n = 0;

        switch (version) {
            case 0:
                // Read the base class Reactor.
                this->Reactor::Deserialize(in, mech);

                // Read the residence time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_restime = (double)val;

                // Calculate inverse RT.
                m_invrt = m_restime;

                // Read the inflow mixture.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_in = new Mops::FlowStream(in, mech);
                }

                // Recalculate enthalpy.
                m_infH = m_in->Mixture()->GasPhase().BulkH();
                m_in->Mixture()->GasPhase().Hs(m_infHs);
                
                break;
            default:
                throw runtime_error("Reactor serialized version number "
                                    "is invalid (Mops, PSR::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Mops, PSR::Deserialize).");
    }
}

// Identifies the mixture type for serialisation.
Serial_ReactorType PSR::SerialType() const
{
    return Serial_PSR;
}


// GOVERNING EQUATIONS.

// Definition of RHS form for constant temperature energy equation.
void PSR::RHS_ConstT(double t, const double *const y,  double *ydot) const
{
    static fvector wdot, sdot;
    double wtot = 0.0, stot =0.0;

    // Calculate molar production rates.
    wtot = m_mech->GasMech().Reactions().GetMolarProdRates(y[m_iT], y[m_iDens], y,
                                                 m_nsp, m_mix->GasPhase(), wdot);
	
	stot = m_mech->GasMech().Reactions().GetSurfaceMolarProdRates(y[m_iT], y[m_iDens], y,
                                                 m_nsp, m_mix->GasPhase(), sdot);

    // Calculate mole fraction derivatives of gas-phase species.
	for (unsigned int i=0; i!=m_mech->GasMech().GasSpeciesCount(); ++i) {
	    ydot[i] = wdot[i] - y[i] * wtot;
	    // Deal with the flow terms
	    ydot[i] += (m_in->Mixture()->GasPhase().Density()
	            * m_in->Mixture()->GasPhase().MoleFraction(i) - y[m_iDens] * y[i])
	            / m_restime;
	    // Add any surface terms
	    if (Area > 0.0)
	        ydot[i] += (sdot[i] - y[i] * stot) * Area / Volume;
	    // Normalise to mole fractions
	    ydot[i] /= y[m_iDens];
	}

	// Calculate mole fraction derivatives of surface species.
	if (Area > 0.0) {
        for (unsigned int i=m_mech->GasMech().GasSpeciesCount(); i!=m_nsp; ++i) {
        ydot[i] = sdot[i]
                * double(m_mech->GasMech().FindSiteOccup(m_mech->GasMech().Species(i)->Name()))
                / m_mech->GasMech().FindSiteDensity(m_mech->GasMech().Species(i)->PhaseName());
        }
	}

    // Temperature derivative.
    if (m_Tfunc) {
        // Add imposed temperature gradient, if defined.
        ydot[m_iT] = m_Tfunc(t, y, ydot, *this);
    } else {
        // Constant temperature.
        ydot[m_iT] = 0.0;
    }

    // Density derivative.
    if (m_constv) {
        // Constant volume.
        ydot[m_iDens] = wtot + (m_invrt * (m_in->Mixture()->GasPhase().Density() - y[m_iDens]));
        if (Area > 0.0) ydot[m_iDens] += stot * Area / (Volume);

	} else {
        // Constant pressure.
        ydot[m_iDens] = 0.0; // NOT SURE
    }
}

/*!
 * Gets the bulk inflow enthalpy. Best to access via this method, as if the
 * mixture is dynamic the cached value could be incorrect.
 *
 * @return  Bulk enthalpy of the inflow
 */
double PSR::InflowBulkEnthalpy() const {
    double infH(0.0);
    if (m_in != NULL) {
        if (m_in->HasReacInflow()) infH = m_in->Mixture()->GasPhase().BulkH() /
                (Sprog::R * m_in->Mixture()->GasPhase().Temperature());
        else infH = m_infH;
    } else {
        throw std::runtime_error("Couldn't find an inflow to calculate H for."
                " Mops, (PSR::InflowBulkEnthalpy)");
    }
    return infH;
}

// Definition of RHS form for adiabatic energy equation.
void PSR::RHS_Adiabatic(double t, const double *const y,  double *ydot) const
{
    static fvector wdot, Hs, Cps, sdot;
    double wtot = 0.0, Cp = 0.0, H = 0.0, stot =0.0, avrMW = 0.0;

    // Calculate mixture thermodynamic properties.
    m_mix->GasPhase().CalcHs_RT(y[m_iT], Hs);
    H = m_mix->GasPhase().BulkH();
    Cp = m_mix->GasPhase().CalcBulkCp_R(y[m_iT], y, m_nsp); // suppose to take all species

    // Calculate molar production rates of species (mol/m3s).
    wtot = m_mech->GasMech().Reactions().GetMolarProdRates(y[m_iT], y[m_iDens], y,
                                                 m_nsp, m_mix->GasPhase(), wdot);

	
	stot = m_mech->GasMech().Reactions().GetSurfaceMolarProdRates(y[m_iT], y[m_iDens], y,
                                                 m_nsp, m_mix->GasPhase(), sdot);

	 for(unsigned int i=0; i!= m_mech->GasMech().GasSpeciesCount(); i++) {
        avrMW += y[i]*m_mech->GasMech().Species(i)->MolWt();

   }
   
    // Calculate mole fraction and temperature derivatives.
    ydot[m_iT] = 0.0;

	if (Area == 0){ 
    for (unsigned int i=0; i!=m_nsp; ++i) {
        // Mole fraction derivative.
        ydot[i] = ((wdot[i] - (y[i]*wtot)) +
                  // Inflow/Outflow term:
                  (m_in->Mixture()->GasPhase().Density() * m_invrt *
                   (m_in->Mixture()->GasPhase().MoleFraction(i) - y[i]))) / y[m_iDens];
	}
	} else{
		for (unsigned int i=0; i!=m_mech->GasMech().GasSpeciesCount(); ++i) {
		ydot[i] = ((wdot[i] - (y[i]*wtot)) / y[m_iDens]) + (sdot[i] * Area /(Volume * y[m_iDens]))   -  ( y[i] * stot * Area/( y[m_iDens] * Volume) ) +
                  // Inflow/Outflow term:
                  (m_in->Mixture()->GasPhase().Density() * m_invrt *
                   (m_in->Mixture()->GasPhase().MoleFraction(i) - y[i]) ) / y[m_iDens];
		}

		for (unsigned int i=m_mech->GasMech().GasSpeciesCount(); i!=m_nsp; ++i) {
			string speciesName = m_mech->GasMech().Species(i)->Name();
		   string phaseName =  m_mech->GasMech().Species(i)->PhaseName();
		   double siteDensity = m_mech->GasMech().FindSiteDensity(phaseName);
		   int siteOccupancy =  m_mech->GasMech().FindSiteOccup(speciesName);
		   ydot[i] = sdot[i]*siteOccupancy / siteDensity;
		}	

	}

	for (unsigned int i=0; i!=m_nsp; ++i) {
        // Temperature derivative.
        ydot[m_iT] += Volume * wdot[i] * Hs[i] + Area * sdot[i] * Hs[i] /* * m_mech->GasMech().Species(i)->MolWt()*/;
    }

    // Complete temperature derivative (including inflow/outflow term).
    ydot[m_iT] *= - y[m_iT] / (Cp * y[m_iDens] * Volume);
    ydot[m_iT] += (m_in->Mixture()->GasPhase().Density() / (y[m_iDens] * Cp * m_restime)) *
                  ((H*y[m_iT]) - (InflowBulkEnthalpy() * m_in->Mixture()->GasPhase().Temperature()));

    // Add imposed temperature gradient, if defined.
    if (m_Tfunc) ydot[m_iT] += m_Tfunc(t, y, ydot, *this);

    // Calculate density derivative.
    if (m_constv) {
		// Constant volume.
		if (Area != 0){
        ydot[m_iDens] = wtot  + stot * Area / (Volume) + (m_invrt * (m_in->Mixture()->GasPhase().Density() - y[m_iDens]));
		}
		else{
		 ydot[m_iDens] =  wtot  + (m_invrt * (m_in->Mixture()->GasPhase().Density() - y[m_iDens]));
		}
	} else {
        // Constant pressure (use EoS to evaluate).(THIS INCLUDES THE SURFACE SOURCE TERM in dT/dt) 
        ydot[m_iDens] = - y[m_iDens] * ydot[m_iT] / y[m_iT];
    }
}

} // Mops namespace
