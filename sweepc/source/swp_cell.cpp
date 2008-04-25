#include "swp_cell.h"
#include "swp_mechanism.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
Cell::Cell(void)
: m_smpvol(1.0), m_fixed_chem(false)
{
}

// Default constructor (public).
Cell::Cell(const Sprog::SpeciesPtrVector &sp)
: Sprog::Thermo::IdealGas(sp), m_smpvol(1.0), m_fixed_chem(false)
{
}

// Copy constructor.
Cell::Cell(const Cell &copy)
{
    *this = copy;
}

// Stream-reading constructor.
Cell::Cell(std::istream &in, const Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
Cell::~Cell(void)
{
}


// OPERATOR OVERLOADS.

// Assignment operator.
Cell &Cell::operator=(const Sweep::Cell &rhs)
{
    if (this != &rhs) {
        Sprog::Thermo::IdealGas::operator=(rhs);
        m_ensemble   = rhs.m_ensemble;
        m_smpvol     = rhs.m_smpvol;
        m_fixed_chem = rhs.m_fixed_chem;
    }
    return *this;
}

// Assignment operator - ideal gas.
Cell &Cell::operator=(const Sprog::Thermo::IdealGas &rhs)
{
    if (this != &rhs) {
        Sprog::Thermo::IdealGas::operator=(rhs);
    }
    return *this;
}

// THE GAS-PHASE INTERFACE.

// Returns the description of the gas-phase mixture.
const Sprog::Thermo::IdealGas &Cell::GasPhase(void) const
{
    return *this;
}

// Sets the gas-phase mixture.
void Cell::SetGasPhase(const Sprog::Thermo::IdealGas &gas)
{
    Sprog::Thermo::IdealGas::operator=(gas);
}

// Adjusts the concentration of the ith species.
void Cell::AdjustConc(unsigned int i, real dc)
{
    if (!m_fixed_chem) {
        unsigned int k;

        // Precalculate DC / density.
        real dc_rho = dc / *m_pdens;

        // Calculate change to all mole fractions k < i.
        for (k=0; k<i; ++k) {
            m_data[k] -= dc_rho * m_data[k];
        }

        // Calculate change for ith species.
        m_data[i] += dc_rho * (1.0 - m_data[i]);

        // Calculate change for all mole fractions k > i.
        for (k=i+1; k<m_species->size(); ++k) {
            m_data[k] -= dc_rho * m_data[k];
        }
    }
}

// Adjusts the concentration of all species.
void Cell::AdjustConcs(const fvector &dc)
{
    if (!m_fixed_chem) {
        // Calculate total change in density.
        real drho = 0.0;
        unsigned int k;
        for (k=0; k!=m_species->size(); ++k) {
            drho += dc[k];
        }

        // Calculate changes to the mole fractions.
        real invrho = 1.0 / *m_pdens;
        for (k=0; k!=m_species->size(); ++k) {
            m_data[k] += (invrho * dc[k]) - (invrho * m_data[k] * drho);
        }
    }
}


// THE PARTICLE ENSEMBLE.

// Returns the particle ensemble.
Ensemble &Cell::Particles(void) {return m_ensemble;}
const Ensemble &Cell::Particles(void) const {return m_ensemble;}

// Returns the particle count.
unsigned int Cell::ParticleCount(void) const 
{
    return m_ensemble.Count();
}

// Returns particle statistics.
void Cell::GetVitalStats(EnsembleStats &stats) const
{
    stats.Calculate(m_ensemble, 1.0/SampleVolume());
}


// SCALING ROUTINES INCL. SAMPLE VOLUME.

// Returns the real system to stochastic system scaling factor.
real Cell::SampleVolume() const
{
    return m_smpvol * m_ensemble.Scaling();
}

// Sets the number density which the ensemble represents.
int Cell::SetM0(const real m0)
{
    if ((m_ensemble.Count() > 0) && (m0 > 0.0)) {
        m_smpvol = (real)m_ensemble.Count() / m0;
        m_ensemble.ResetScaling();
        return 0;
    }// else {
    //    // The ensemble contains no particles, so assume this
    //    // is the maximum M0.        
    //    return SetMaxM0(m0);
    //}
    return 1;
}

// Sets the number density which the full 
// ensemble would represent.
int Cell::SetMaxM0(const real m0)
{
    if ((m_ensemble.Capacity() > 0) && (m0 > 0.0)) {
        m_smpvol = m_ensemble.Scaling() * (real)m_ensemble.Capacity() / m0;
        m_ensemble.ResetScaling();
        return 0;
    } else {
        // The ensemble has not yet been initialised, hence guess
        // unit volume and report an error.
        m_smpvol = 1.0;
        return -1;
    }
}

void Cell::Reset(const real m0)
{
    m_ensemble.Clear();
    SetMaxM0(m0);
}


// FIXED/VARIABLE CHEMISTRY.

// Returns whether or not the chemical conditions are fixed.
bool Cell::FixedChem() const {return m_fixed_chem;}

// Sets whether or not the chemical conditions are fixed.
void Cell::SetFixedChem(bool fixed) {m_fixed_chem = fixed;}

// Set the chemical conditions to be variable.
void Cell::SetVariableChem(bool vari) {m_fixed_chem = !vari;}


// READ/WRITE/COPY.

// Writes the object to a binary stream.
void Cell::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output the base class.
        Sprog::Thermo::IdealGas::Serialize(out);

        // Output the sample volume.
        double v = (double)m_smpvol;
        out.write((char*)&v, sizeof(v));

        // Output if fixed chem.
        out.write((char*)&m_fixed_chem, sizeof(m_fixed_chem));

        // Output the ensemble.
        m_ensemble.Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Cell::Serialize).");
    }
}

// Reads the object from a binary stream.
void Cell::Deserialize(std::istream &in, const Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;

        switch (version) {
            case 0:
                // Read the base class.
                Sprog::Thermo::IdealGas::Deserialize(in);

                // Read the sample volume.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_smpvol = (real)val;

                // Read if fixed chem.
                in.read(reinterpret_cast<char*>(&m_fixed_chem), sizeof(m_fixed_chem));

                // Read the ensemble.
                m_ensemble.Deserialize(in, mech);

                // Set the species.
                m_species = mech.Species();

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Cell::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Cell::Deserialize).");
    }
}
