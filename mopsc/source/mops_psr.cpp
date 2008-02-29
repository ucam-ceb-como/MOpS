#include "mops_psr.h"
#include "mops_mixture.h"

#include <vector>
#include <cmath>
#include <stdexcept>

using namespace Mops;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
PSR::PSR(void)
{
    init();
}

// Default constructor (public, requires mechanism).
PSR::PSR(const Mops::Mechanism &mech)
: Reactor(mech)
{
    init();
    m_infHs.resize(m_neq);
}

// Copy constructor.
PSR::PSR(const Mops::PSR &copy)
{
    m_restime = copy.m_restime;
    m_inflow  = copy.m_inflow->Clone();
    m_infH    = copy.m_infH;
    m_infHs.assign(copy.m_infHs.begin(), copy.m_infHs.end());
}

// Stream-reading constructor.
PSR::PSR(std::istream &in, const Mops::Mechanism &mech)
: Reactor(mech)
{
    init();
    Deserialize(in, mech);
}

// Default destructor.
PSR::~PSR(void)
{
    releaseMemory();
}


// RESIDENCE TIME.

// Returns the reactor residence time.
real PSR::ResidenceTime(void) const
{
    return m_restime;
}

// Sets the reactor residence time.
void PSR::SetResidenceTime(real t)
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

// Returns the mixture which describes the inflow conditions.
const Mops::Mixture *const PSR::Inflow(void) const
{
    return m_inflow;
}

// Sets the mixture which describes the inflow conditions.
void PSR::SetInflow(Mops::Mixture &inf)
{
    // Delete current inflow mixture.
    if (m_inflow != NULL) delete m_inflow;

    m_inflow = &inf;

    // Calculate inflow enthalpies to speed up
    // later calculations.
    m_infH = m_inflow->BulkH() / (Sprog::R * m_inflow->Temperature());
    m_inflow->Hs_RT(m_infHs);
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

        // Output the inflow mixture.
        if (m_inflow != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            m_inflow->Serialize(out);
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
    // Clear current reactor data.
    releaseMemory();
    init();

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
                m_restime = (real)val;

                // Calculate inverse RT.
                m_invrt = m_restime;

                // Read the inflow mixture.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_inflow = new Mops::Mixture(in, mech);
                }

                // Recalculate enthalpy.
                m_infH = m_inflow->BulkH();
                m_inflow->Hs(m_infHs);
                
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
void PSR::RHS_ConstT(real t, const real *const y,  real *ydot)
{
    fvector wdot;
    real wtot = 0.0;

    // Calculate molar production rates.
    wtot = m_mech->Reactions().GetMolarProdRates(y[m_iT], y[m_iDens], y, 
                                                 m_nsp, *m_mix, wdot);

    // Calculate mole fraction derivatives.
    for (unsigned int i=0; i!=m_nsp; ++i) {
        ydot[i] = ((wdot[i] - (y[i]*wtot)) / y[m_iDens]) +
                  // Inflow/Outflow term:
                  (m_inflow->Density() * m_invrt * 
                   (m_inflow->MoleFraction(i) - y[i]));
    }

    // Temperature derivative.
    ydot[m_iT] = 0.0;

    // Density derivative.
    if (m_constv) {
        // Constant volume.
        ydot[m_iDens] = wtot + (m_invrt * (m_inflow->Density() - y[m_iDens]));
    } else {
        // Constant pressure.
        ydot[m_iDens] = 0.0;
    }
}

// Definition of RHS form for adiabatic energy equation.
void PSR::RHS_Adiabatic(real t, const real *const y,  real *ydot)
{
    fvector wdot, Hs, Cps;
    real wtot = 0.0, Cp = 0.0, H = 0.0;

    // Calculate mixture thermodynamic properties.
    m_mix->CalcHs_RT(y[m_iT], Hs);
    H = m_mix->BulkH();
    Cp = m_mix->CalcBulkCp(y[m_iT], y, m_nsp, Cps) / Sprog::R;

    // Calculate molar production rates of species (mol/m3s).
    wtot = m_mech->Reactions().GetMolarProdRates(y[m_iT], y[m_iDens], y, 
                                                 m_nsp, *m_mix, wdot);

    // Calculate mole fraction and temperature derivatives.
    ydot[m_iT] = 0.0;
    for (unsigned int i=0; i!=m_nsp; ++i) {
        // Mole fraction derivative.
        ydot[i] = ((wdot[i] - (y[i]*wtot)) / y[m_iDens]) +
                  // Inflow/outflow term:
                  (m_inflow->Density() * m_invrt * 
                   (m_inflow->MoleFraction(i) - y[i]));

        // Temperature derivative.
        ydot[m_iT] += Hs[i] * wdot[i];
    }

    // Complete temperature derivative (including inflow/outflow term).
    ydot[m_iT] *= - y[m_iT] / (Cp * y[m_iDens]);
    ydot[m_iT] += (m_inflow->Density() / (y[m_iDens] * Cp * m_restime)) * 
                  ((H*y[m_iT]) - (m_infH * m_inflow->Temperature()));

    // Calculate density derivative.
    if (m_constv) {
        // Constant volume.
        ydot[m_iDens] = wtot + (m_invrt * (m_inflow->Density() - y[m_iDens]));
    } else {
        // Constant pressure (use EoS to evaluate).
        ydot[m_iDens] = - y[m_iDens] * ydot[m_iT] / y[m_iT];
    }
}


// INITIALISATION.

// Initialises the PSR to a default state, this is called
// by the constructors.
void PSR::init(void)
{
    m_restime = 0.0;
    m_inflow  = NULL;
    m_invrt   = 0.0;
    m_infH    = 0.0;
    m_infHs.clear();
}

// Releases all object memory.
void PSR::releaseMemory(void)
{
    if (m_inflow != NULL) delete m_inflow;
    m_infHs.clear();
}
