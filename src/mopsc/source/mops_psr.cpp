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
  m_inflow_ptrs(),
  m_outflow_ptrs(),
  m_invrt(1.0),
  m_iscaling(1.0),
  m_oscaling(1.0),
  m_default_birth(Sweep::Processes::BirthProcess::iStochastic),
  m_default_death(Sweep::Processes::DeathProcess::iContRescale)
{
}

// Copy constructor.
PSR::PSR(const Mops::PSR &copy)
: Reactor(copy),
  m_restime(copy.m_restime),
  m_inflow_ptrs(copy.m_inflow_ptrs),
  m_outflow_ptrs(copy.m_outflow_ptrs),
  m_invrt(copy.m_invrt),
  m_iscaling(copy.m_iscaling),
  m_oscaling(copy.m_oscaling),
  m_default_birth(copy.m_default_birth),
  m_default_death(copy.m_default_death)
{
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
    // Parent destructor automatically called.
}

// OPERATORS.

Mops::PSR &PSR::operator=(const Mops::PSR &rhs)
{
    if (this != &rhs) {
        // Copy base class.
        Reactor::operator=(rhs);

        // Residence time.
        m_restime = rhs.m_restime;
        m_invrt = rhs.m_invrt;
        m_iscaling = rhs.m_iscaling;
        m_oscaling = rhs.m_oscaling;

        // Create links to IO streams
        for (Mops::FlowPtrVector::const_iterator it=rhs.m_inflow_ptrs.begin();
                it!=rhs.m_inflow_ptrs.end(); ++it) {
            m_inflow_ptrs.push_back(*it);
        }
        for (Mops::FlowPtrVector::const_iterator it=rhs.m_outflow_ptrs.begin();
                it!=rhs.m_outflow_ptrs.end(); ++it) {
            m_outflow_ptrs.push_back(*it);
        }

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
            " Include particles in energy model=" << r.IncludeParticles() <<
            "\n";
    os << "filled with: " << *(r.Mixture());
    for (Mops::FlowPtrVector::const_iterator it=r.m_inflow_ptrs.begin();
            it!=r.m_inflow_ptrs.end(); ++it)
        os << "with inflow: " << *((*it)->Mixture());
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

// Returns a inflow stream
Mops::FlowStream *const PSR::Inflow(unsigned int i) const {
    assert(i < m_inflow_ptrs.size());
    return m_inflow_ptrs[i];
}

// Returns the inflow stream pointers
// aab64: changed definition from 
// Mops::FlowPtrVector PSR::Inflows() const {
const Mops::FlowPtrVector & PSR::Inflows() const {
    return m_inflow_ptrs;
}

// Returns a outflow stream
Mops::FlowStream *const PSR::Outflow(unsigned int i) const {
    assert(i < m_outflow_ptrs.size());
    return m_outflow_ptrs[i];
}

// Returns the Outflow stream pointers
// aab64: changed definition from 
//Mops::FlowPtrVector PSR::Outflows() const {
const Mops::FlowPtrVector & PSR::Outflows() const {
    return m_outflow_ptrs;
}


/*!
 * Create the birth process for a given inflow
 *
 * @param inf   Inflow stream
 */
void PSR::InitialiseInflow(Mops::FlowStream& inf) {
    // Create a new birth process
    Sweep::Processes::BirthProcess bp(Mech()->ParticleMech());
    bp.SetBirthType(m_default_birth);
    bp.SetCell(inf.Mixture());
    bp.SetA(inf.GetFlowFraction() * m_invrt);

    // If an upstream process is a move process, ensure the inflow is turned-off
    if (inf.Mixture()->OutflowCount() > 0) {
        const Sweep::Processes::DeathPtrVector &outf = inf.Mixture()->Outflows();
        for (Sweep::Processes::DeathPtrVector::const_iterator i = outf.begin();
                i != outf.end(); ++i) {
            if ((*i)->GetDeathType() == Sweep::Processes::DeathProcess::iContMove
                    || (*i)->GetDeathType() == Sweep::Processes::DeathProcess::iStochMove) {
                bp.SetProcessSwitch(false);
            }
        }
    }

    // Add the process to the Mixture, which clones it and takes ownership
    m_mix->AddInflow(bp);
}

/*!
 * Create the death process for a given flowstream
 *
 * @param outf  Outflow stream
 */
void PSR::InitialiseOutflow(Mops::FlowStream& outf) {
    Sweep::Processes::DeathProcess dp(Mech()->ParticleMech());

    // Set the form of the death process.
    dp.SetDeathType(m_default_death);
    dp.SetA(outf.GetFlowFraction() * m_invrt);

    // Set the outlet cell if a reactor is present for outflow
    if (HasOutflow()) {
        if (outf.HasReacOutflow()) dp.SetCell(outf.Outflow()->Mixture());
    }

    m_mix->AddOutflow(dp);
}

/*!
 * Initialise the inflows' birth processes
 */
void PSR::InitialiseInflows() {

    // Loop over the inflows
    for (Mops::FlowPtrVector::const_iterator it=m_inflow_ptrs.begin();
            it != m_inflow_ptrs.end(); ++it) {
        InitialiseInflow(*(*it));
    }
}

/*!
 * Initialise the outflows' death processes
 */
void PSR::InitialiseOutflows() {

    // Loop over the outflows
    for (Mops::FlowPtrVector::const_iterator it=m_outflow_ptrs.begin();
            it != m_outflow_ptrs.end(); ++it) {
        InitialiseOutflow(*(*it));
    }
}

/*!
 * Sets the flowstream which points to the inflow mixture conditions.
 *
 * @param inf   Inflow stream
 */
void PSR::SetInflow(Mops::FlowStream &inf)
{
    m_inflow_ptrs.push_back(&inf);

    // Recalculate the scaling factor
    double sum(0.0);
    for (Mops::FlowPtrVector::const_iterator it=m_inflow_ptrs.begin();
            it!=m_inflow_ptrs.end(); ++it) {
        sum += (*it)->GetFlowFraction();
    }
    assert(sum > 0.0);
    m_iscaling = 1.0/sum;

    // Initialise the inflow
    InitialiseInflow(inf);
}

/*!
 * Sets the reactors out pointer to an outflow stream, then creates a new
 * death process in the reactor's mixture.
 *
 * @param out   Outflow stream
 */
void PSR::SetOutflow(Mops::FlowStream &out) {
    m_outflow_ptrs.push_back(&out);

    // Recalculate the scaling factor
    double sum(0.0);
    for (Mops::FlowPtrVector::const_iterator it=m_outflow_ptrs.begin();
            it!=m_outflow_ptrs.end(); ++it) {
        sum += (*it)->GetFlowFraction();
    }
    assert(sum > 0.0);
    m_oscaling = 1.0/sum;

    InitialiseOutflow(out);
}

// Does the reactor have an outflow set?
bool PSR::HasOutflow() const {
    if (m_outflow_ptrs.size() > 0) return true;
    else return false;
}

// Does the reactor have an inflow set?
bool PSR::HasInflow() const {
    if (m_inflow_ptrs.size() > 0) return true;
    else return false;
}

/*!
 * Ensure that all the 'fractional flowrates' specified for the FlowStreams
 * add up to 1.0, for inflows and outflows. This is necessary to check that
 * each reactor receives the correct volumetric flow.
 */
void PSR::NormaliseIOProcessRates() {
    // First loop over inflow processes
    Sweep::Processes::BirthPtrVector bps = m_mix->Inflows();
    for (Sweep::Processes::BirthPtrVector::iterator it=bps.begin();
            it!=bps.end(); ++it) {
        (*it)->SetA((*it)->A() * m_iscaling);
    }

    // Then do the outflow processes
    Sweep::Processes::DeathPtrVector dps = m_mix->Outflows();
    for (Sweep::Processes::DeathPtrVector::iterator it=dps.begin();
            it!=dps.end(); ++it) {
        (*it)->SetA((*it)->A() * m_oscaling);
    }
}

/*!
 * Clears the memory associated with any linked inflow or outflow streams. Only
 * to be used with caution. This is a result of transferring ownership of
 * flowstreams to objects outside of the PSR (e.g. network).
 */
void PSR::ClearStreamMemory() {
    std::cout << "mops: Clearing in/outflow streams and their mixtures." << std::endl;
    // Clear all streams connected to this reactor
    for (Mops::FlowPtrVector::iterator it=m_inflow_ptrs.begin();
            it!=m_inflow_ptrs.end(); ++it) {
        delete (*it);
    }
    for (Mops::FlowPtrVector::iterator it=m_outflow_ptrs.begin();
            it!=m_outflow_ptrs.end(); ++it) {
        delete (*it);
    }
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
    if (out.good()) {
        // Output base class Reactor.
        this->Reactor::Serialize(out);

        // Output the residence time.
        double val = (double)m_restime;
        out.write((char*)&val, sizeof(val));

    } else {
        throw invalid_argument("Output stream not ready (Mops, PSR::Serialize).");
    }
}

// Reads the Reactor data from a binary data stream.
void PSR::Deserialize(std::istream &in, const Mops::Mechanism &mech)
{

    if (in.good()) {

        double val = 0.0;
        // Read the base class Reactor.
        this->Reactor::Deserialize(in, mech);

        // Read the residence time.
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_restime = (double)val;

        // Calculate inverse RT.
        m_invrt = 1.0/m_restime;

    } else {
        throw invalid_argument("Input stream not ready (Mops, PSR::Deserialize).");
    }
}

// Identifies the mixture type for serialisation.
Serial_ReactorType PSR::SerialType() const
{
    return Serial_PSR;
}

/*!
 * Provides the RHS for the ODE solver in a PSR. Because of the complex
 * inter-dependency of many of the equations upon each other, all PSRs
 * have their ydots calculated here.
 *
 * 1. Get wdot, wtotal
 * 2. Get dT_dt (needed by gas-phase expansion factor)
 * 3. Calculate gamma, gas-phase expansion factor
 * 4. Get drho_dt
 * 6. Get dy_k_dt for each speices
 *
 * TODO: Verify this equation for surface chemistry.
 *
 * @param t     Time
 * @param y     Solution at the time
 * @param ydot  First derivative of the solution.
 */
void PSR::RHS_Complete(double t, const double *const y, double *ydot) const {
	
	// aab64: Debug loop for inflow pointers checking the pointer to the 
	// inflow mixture doesn't go NULL
	for (Mops::FlowPtrVector::const_iterator it = m_inflow_ptrs.begin();
		it != m_inflow_ptrs.end(); ++it) {
		assert(&(*it)->Mixture()->GasPhase());
	}

    if (m_sarea > 0.0)
        throw std::runtime_error("PSR untested for surface chemistry.");

    // Calculate molar production rates.
    fvector wdot, sdot;
    double wtot(0.0); double stot(0.0);
    wtot = m_mech->GasMech().Reactions().GetMolarProdRates(y[m_iT], y[m_iDens], y,
                                                 m_nsp, m_mix->GasPhase(), wdot);
    if (m_sarea > 0.0) {
        stot = m_mech->GasMech().Reactions().GetSurfaceMolarProdRates(y[m_iT], y[m_iDens], y,
                                                     m_nsp, m_mix->GasPhase(), sdot);
    }

    // Calculate mixture thermodynamic properties.
    ydot[m_iT] = 0.0;
    if (m_emodel == Reactor::Adiabatic) {
        double C(0.0);
        fvector H_reac, H_wdot;
        // We need to use enthalpy for the inflow terms.
        // If it's constant volume, we use internal energy for the production
        // terms, otherwise use enthalpy.
        m_mix->GasPhase().CalcHs(y[m_iT], H_reac);
        if (m_constv) {
            C = m_mix->GasPhase().ThermoInterface::CalcBulkCv_R(y[m_iT], y, m_nsp);
            m_mix->GasPhase().CalcUs(y[m_iT], H_wdot);
        } else {
            C = m_mix->GasPhase().ThermoInterface::CalcBulkCp_R(y[m_iT], y, m_nsp);
            m_mix->GasPhase().CalcHs(y[m_iT], H_wdot);
        }

        // Loop over incoming streams to get inflow energy change
        double hsum(0.0);
        for (Mops::FlowPtrVector::const_iterator it=m_inflow_ptrs.begin();
                it!=m_inflow_ptrs.end(); ++it) {
            double hval(0.0);
            fvector H_in;
            fvector fracs = (*it)->Mixture()->GasPhase().MoleFractions();

            // Get incoming enthalpies.
            (*it)->Mixture()->GasPhase().CalcHs((*it)->Mixture()->GasPhase().Temperature(), H_in);

            // Then get contribution to temperature change from this stream
            for (unsigned int i=0; i!=m_nsp; ++i) {
                hval += fracs[i] * (H_in[i] - H_reac[i]);
            }
            H_in.clear();

            hsum += hval * (*it)->Mixture()->GasPhase().Density()
                    * (*it)->GetFlowFraction();
        }
        // Account for the particles in thermal bulk term on denominator
        if (m_include_particle_terms) 
        {
            ydot[m_iT] += hsum * m_invrt * m_iscaling / (y[m_iDens] * C * Sprog::R + m_mix->getParticleDensity());
        }
        else
        {
            ydot[m_iT] += hsum * m_invrt * m_iscaling / (y[m_iDens] * C * Sprog::R);
        }
        // Loop over species in the reactor to get energy change due to mol change
        hsum = 0.0;
        for (unsigned int i=0; i!=m_nsp; ++i) {
            hsum -= wdot[i] * H_wdot[i];
            if (m_sarea > 0.0) hsum -= sdot[i] * H_wdot[i] * m_sarea / m_svol;
        }
        if (m_include_particle_terms)
        {
            ydot[m_iT] += hsum / (y[m_iDens] * C * Sprog::R + m_mix->getParticleDensity());
        }
        else
        {
            ydot[m_iT] += hsum / (y[m_iDens] * C * Sprog::R);
        }
    }
    // Add imposed temperature gradient, if defined
    if (m_Tfunc) ydot[m_iT] += m_Tfunc(t, y, ydot, *this);

    // Now we need to get dn_dt.
    // Here, dn_dt is used to refer to (1/V) * dn/dt
    double dn_dt(0.0);
    for (Mops::FlowPtrVector::const_iterator it=m_inflow_ptrs.begin();
            it!=m_inflow_ptrs.end(); ++it) {
        // Debug check for inflow pointer
        assert(&(*it)->Mixture()->GasPhase());
		
        dn_dt += (*it)->GetFlowFraction() * (*it)->Mixture()->GasPhase().Density();
    }
    dn_dt = (dn_dt * m_iscaling - y[m_iDens]) * m_invrt + wtot;

    // Can now get the gas-phase expansion factor
    // gamma = (1/V) * dV/dt
    double gamma(0.0);
    if (!m_constv) gamma = dn_dt / y[m_iDens] + ydot[m_iT] / y[m_iT];

    // So let's get the density derivative now
    ydot[m_iDens] = 0.0;
    if (m_constv) ydot[m_iDens] = dn_dt;
    else ydot[m_iDens] = -y[m_iDens] * ydot[m_iT] / y[m_iT];
    if (m_sarea > 0.0)  ydot[m_iDens] += stot * m_sarea / m_svol;

    // And finally, we can calculate the changes to species fractions
    for (unsigned int i=0; i!=m_nsp; ++i) {
        ydot[i] = (wdot[i] - y[i] * ydot[m_iDens]) / y[m_iDens];
        if (m_sarea > 0.0) ydot[i] += (sdot[i] - (y[i] * stot)) * m_sarea / (m_svol * y[m_iDens]);

        // Calculate inflow contribution
        double yval(0.0);
        for (Mops::FlowPtrVector::const_iterator it=m_inflow_ptrs.begin();
                it!=m_inflow_ptrs.end(); ++it) {
            yval += (*it)->Mixture()->GasPhase().Density()
                    * (*it)->GetFlowFraction()
                    * (*it)->Mixture()->GasPhase().MoleFraction(i);
        }
        yval *= m_iscaling;
        yval -= y[m_iDens] * y[i];
        ydot[i] += yval * m_invrt / y[m_iDens];

        // Subtract the expansion term
        ydot[i] -= y[i] * gamma;
    }

    // Get the ydot for some surface species.
    if (m_sarea > 0.0) {
        for (unsigned int i=m_mech->GasMech().GasSpeciesCount(); i!=m_neq-2; ++i) {
            ydot[i] = ((double) m_mech->GasMech().FindSiteOccup(m_mech->GasMech().Species(i)->Name()))
                   * sdot[i] / m_mech->GasMech().FindSiteDensity(m_mech->GasMech().Species(i)->PhaseName());
        }
    }
}

// GOVERNING EQUATIONS.

/*!
 * Supplies the RHS for the ODE solver for a constant temperature reactor.
 *
 * @param t     Time
 * @param y     Solution at the time
 * @param ydot  First derivative of the solution.
 */
void PSR::RHS_ConstT(double t, const double *const y,  double *ydot) const
{
    PSR::RHS_Complete(t, y, ydot);
}

/*!
 * Supplies the RHS for the ODE solver for a adiabatic reactor.
 *
 * @param t     Time
 * @param y     Solution at the time
 * @param ydot  First derivative of the solution.
 */
void PSR::RHS_Adiabatic(double t, const double *const y,  double *ydot) const
{
    PSR::RHS_Complete(t, y, ydot);
}

} // Mops namespace
