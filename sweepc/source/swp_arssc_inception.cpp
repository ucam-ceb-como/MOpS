#include "swp_arssc_inception.h"
#include "swp_arssc_model.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
ARSSC_Inception::ARSSC_Inception(void)
: m_sites(SubModels::ARSSC_Model::SiteTypeCount,0.0)
{
}

// Initialising constructor.
ARSSC_Inception::ARSSC_Inception(const Sweep::Mechanism &mech)
: Inception(mech), m_sites(SubModels::ARSSC_Model::SiteTypeCount,0.0)
{
}

// Copy constructor.
ARSSC_Inception::ARSSC_Inception(const ARSSC_Inception &copy)
{
    *this = copy;
}

// Stream-reading constructor.
ARSSC_Inception::ARSSC_Inception(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
ARSSC_Inception::~ARSSC_Inception(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
ARSSC_Inception &ARSSC_Inception::operator =(const ARSSC_Inception &rhs)
{
    if (this != &rhs) {
        Inception::operator=(rhs);
        m_sites.assign(rhs.m_sites.begin(), rhs.m_sites.end());
    }
    return *this;
}


// PERFORMING THE PROCESS.

// Performs the process on the given system.  The responsible rate term is given
// by index.  Returns 0 on success, otherwise negative.
int ARSSC_Inception::Perform(real t, Cell &sys, unsigned int iterm) const 
{
    // This routine performs the inception on the given chemical system.

    // Create a new particle of the type specified
    // by the system ensemble.
    Particle *sp = m_mech->CreateParticle(t);
    
    // Set the new particle's aromatic sites.
    SubModels::ARSSC_Model *ars = 
        static_cast<SubModels::ARSSC_Model*const>(sp->Primary()->SubModel(SubModels::ARSSC_Model_ID));
    for (unsigned int i=0; i!=m_sites.size(); ++i) {
        ars->SetSiteCount((SubModels::ARSSC_Model::SiteType)i, m_sites[i]);
    }

    // Initialise the new particle.
    sp->Primary()->SetComposition(m_newcomp);
    sp->Primary()->SetValues(m_newvals);
    sp->UpdateCache();

    // Add particle to system's ensemble.
    sys.Particles().Add(*sp);

    // Update gas-phase chemistry of system.
    adjustGas(sys);

    return 0;
}


// READ/WRITE/COPY.

// Creates a copy of the inception.
ARSSC_Inception *const ARSSC_Inception::Clone(void) const {return new ARSSC_Inception(*this);}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType ARSSC_Inception::ID(void) const {return ARSSC_Inception_ID;}

// Writes the object to a binary stream.
void ARSSC_Inception::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        Inception::Serialize(out);

        // Write number of sites
        unsigned int n = (unsigned int)m_sites.size();
        out.write((char*)&n, sizeof(n));

        // Write site values.
        for (fvector::const_iterator i=m_sites.begin(); i!=m_sites.end(); ++i) {
            double v = (double)*i;
            out.write((char*)&v, sizeof(v));
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ARSSC_Inception::Serialize).");
    }
}

// Reads the object from a binary stream.
void ARSSC_Inception::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        double val = 0.0;

        switch (version) {
            case 0:
                // Deserialize base class.
                Inception::Deserialize(in, mech);

                // Read new site count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read component values.
                m_sites.clear();
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_sites.push_back((real)val);
                }

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ARSSC_Inception::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ARSSC_Inception::Deserialize).");
    }
}
