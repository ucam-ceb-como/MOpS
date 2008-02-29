#include "swp_surfacereaction.h"
#include "swp_mechanism.h"
#include "swp_processtype.h"
#include <cmath>
#include <stdexcept>

using namespace Sweep;
using namespace std;

const real SurfaceReaction::m_majfactor = 2.0;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
SurfaceReaction::SurfaceReaction(void)
: m_arr(0.0,0.0,0.0), m_pid(0), m_modelid(BasicModel_ID)
{
    m_defer = true;
}

// Copy constructor.
SurfaceReaction::SurfaceReaction(const SurfaceReaction &copy)
{
    *this = copy;
}

// Stream-reading constructor.
SurfaceReaction::SurfaceReaction(std::istream &in)
{
    Deserialize(in);
}

// Default destructor.
SurfaceReaction::~SurfaceReaction(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
SurfaceReaction &SurfaceReaction::operator =(const Sweep::SurfaceReaction &rhs)
{
    if (this != &rhs) {
        ParticleProcess::operator =(rhs);
        m_arr     = rhs.m_arr;
        m_pid     = rhs.m_pid;
        m_modelid = rhs.m_modelid;
    }
    return *this;
}


// ARRHENIUS COEFFICIENTS.

// Returns the Arrhenius parameter.
Sprog::Kinetics::ARRHENIUS &SurfaceReaction::Arrhenius() {return m_arr;}
const Sprog::Kinetics::ARRHENIUS &SurfaceReaction::Arrhenius() const {return m_arr;}

// Sets the Arrhenius parameters.
void SurfaceReaction::SetArrhenius(Sprog::Kinetics::ARRHENIUS &arr) {m_arr = arr;}




// PARTICLE PROPERTY ID.

// Returns the ID number of the particle property to which
// the rate of this process is proportional.
unsigned int SurfaceReaction::PropertyID(void) const {return m_pid;}

// Returns the ID number of the particle number for which the
// PropertyID is valid.  The mechanism should check that this
// model is enabled.
ModelType SurfaceReaction::ModelID(void) const {return m_modelid;}

// Sets the ID number of the particle property to which
// the rate of this process is proportional.
void SurfaceReaction::SetPropertyID(unsigned int i, ModelType modelid)
{
    m_pid     = i;
    m_modelid = modelid;
}


// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

// Returns rate of the process for the given system.
real SurfaceReaction::Rate(real t, const Cell &sys) const
{
    return Rate(t, sys, sys);
}

// Calculates the process rate using the given 
// chemical conditions, rather than those conditions in the
// given system.
real SurfaceReaction::Rate(real t, const Sprog::Thermo::IdealGas &gas, 
                           const Cell &sys) const
{
    // Rate constant.
    real rate = m_arr.A;

    // Chemical species concentration dependence.
    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        real conc = gas.MolarConc(i->first);
        for (int j=0; j<i->second; ++j) {
            rate *= conc;
        }
    }

    // Tempearature dependance.
    real T = gas.Temperature();
    rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));

    // Particle dependence.
    rate *= sys.Particles().GetSum(m_modelid, m_pid);

    if (m_mech->AnyDeferred()) {
        return rate * m_majfactor;
    } else {
        return rate;
    }
}


// SINGLE PARTICLE RATE CALCULATIONS.

// Returns the rate of the process for the given particle in
// the system. Process must be linear in particle number.
real SurfaceReaction::Rate(real t, const Cell &sys, const Particle &sp) const
{
    return Rate(t,sys, sys, sp);
}

// Returns rate of the process for the given particle using the
// given chemical conditions rather than those conditions in the
// the given system.
real SurfaceReaction::Rate(real t, const Sprog::Thermo::IdealGas &gas, 
                           const Cell &sys, const Particle &sp) const
{
    // Rate constant.
    real rate = m_arr.A;

    // Chemical species concentration dependence.
    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        real conc =  gas.MolarConc(i->first);
        for (int j=0; j<i->second; ++j) {
            rate *= conc;
        }
    }

    // Tempearature dependance.
    real T = gas.Temperature();
    rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));

    // Paticle dependence.
    if (m_modelid == BasicModel_ID) {
        rate *= sp.Property(static_cast<ParticleData::PropertyID>(m_pid));
    } else {
        rate *= sp.ModelCache(m_modelid)->Property(m_pid);
    }

    return rate;
}

// Returns majorant rate of the process for the given system.
real SurfaceReaction::MajorantRate(real t, const Cell &sys, 
                                   const Particle &sp) const
{
    return Rate(t, sys, sys, sp) * m_majfactor;
}

// Calculates the majorant process rate using the given 
// chemical conditions, rather than those conditions in the
// given system.
real SurfaceReaction::MajorantRate(real t, const Sprog::Thermo::IdealGas &gas, 
                                   const Cell &sys, const Particle &sp) const
{
    return Rate(t, gas, sys, sp) * m_majfactor;
}


// RATE TERM CALCULATIONS.
//   These routines return the individual rate terms for a 
//   process, which may have multiple terms (e.g. condensation).

// Returns the number of rate terms for this process.
unsigned int SurfaceReaction::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.
real SurfaceReaction::RateTerms(real t, const Cell &sys, 
                                fvector::iterator &iterm) const
{
    return *(iterm++) = Rate(t, sys, sys);
}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  The given chemical conditions are used instead of those
// in the given system object.
real SurfaceReaction::RateTerms(real t, const Sprog::Thermo::IdealGas &gas, 
                                const Cell &sys, fvector::iterator &iterm) const
{
    return *(iterm++) = Rate(t, gas, sys);
}


// PERFORMING THE PROCESS.

// Performs the process on the given system.  The responsible rate term is given
// by index.  Returns 0 on success, otherwise negative.
int SurfaceReaction::Perform(real t, Cell &sys, unsigned int iterm) const
{
    int i = sys.Particles().Select(m_modelid, m_pid);

    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);


        // Update particle with deferred processes.
        if (m_mech->AnyDeferred()) {
            // Calculate majorant rate then update the particle.
            real majr = MajorantRate(t, sys, *sp);
            m_mech->UpdateParticle(*sp, sys, t);

            // Check that the particle is still valid.
            if (sp->IsValid()) {
                real truer = Rate(t, sys, *sp);

                if (!Ficticious(majr, truer)) {
                    // Adjust particle.
                    sp->Adjust(m_dcomp, m_dvals);
                    sys.Particles().Update(i);

                    // Apply changes to gas-phase chemistry.
                    adjustGas(sys);
                }
            } else {
                // If not valid then remove the particle.
                sys.Particles().Remove(i);
            }
        } else {
            // No particle update required, just perform the surface
            // reaction.
            sp->Adjust(m_dcomp, m_dvals);
            sys.Particles().Update(i);

            // Apply changes to gas-phase chemistry.
            adjustGas(sys);
        }
    } else {
        // Failed to select a particle.
        return -1;
    }

    return 0;
}

// Performs the process on a given particle in the system.  Particle
// is given by index.  The process is performed n times.
int SurfaceReaction::Perform(real t, Cell &sys, Particle &sp, 
                             unsigned int n) const
{
    unsigned int m = sp.Adjust(m_dcomp, m_dvals, n);
    adjustGas(sys, m);
    return 0;
}


// READ/WRITE/COPY.

// Creates a copy of the particle process.
SurfaceReaction *const SurfaceReaction::Clone(void) const
{
    return new SurfaceReaction(*this);
}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType SurfaceReaction::ID(void) const {return SurfaceReaction_ID;}

// Writes the object to a binary stream.
void SurfaceReaction::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        ParticleProcess::Serialize(out);

        // Write arrhenius coefficients.
        double A  = (double)m_arr.A;
        double nn = (double)m_arr.n;
        double E  = (double)m_arr.E;
        out.write((char*)&A, sizeof(A));
        out.write((char*)&nn, sizeof(nn));
        out.write((char*)&E, sizeof(E));

        // Write particle property ID.
        unsigned int n = (unsigned int)m_pid;
        out.write((char*)&n, sizeof(n));

        // Write the model ID.
        n = (unsigned int)m_modelid;
        out.write((char*)&n, sizeof(n));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, SurfaceReaction::Serialize).");
    }
}

// Reads the object from a binary stream.
void SurfaceReaction::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double A = 0.0, nn = 0.0, E = 0.0;
        unsigned int n = 0;

        switch (version) {
            case 0:
                // Deserialize base class.
                ParticleProcess::Deserialize(in);

                // Read arrhenius coefficients.
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
                in.read(reinterpret_cast<char*>(&E), sizeof(E));
                m_arr.A = (real)A;
                m_arr.n = (real)nn;
                m_arr.E = (real)E;

                // Read particle property ID.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_pid = n;

                // Read the model ID.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_modelid = (ModelType)n;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, SurfaceReaction::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SurfaceReaction::Deserialize).");
    }
}
