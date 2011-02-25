#include "gpc_fall_off_reaction.h"
#include "gpc_third_body_reaction.h"
#include "gpc_rate_params.h"

using namespace Sprog;
using namespace Sprog::Kinetics;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Sprog::Kinetics::FallOffReaction::FallOffReaction()
{
    m_fothirdbody = -1;
}

// Copy constructor.
Sprog::Kinetics::FallOffReaction::FallOffReaction(const Sprog::Kinetics::FallOffReaction &rxn)
{
    *this = rxn;
}

// Destructor.
Sprog::Kinetics::FallOffReaction::~FallOffReaction()
{
}


// OPERATOR OVERLOADING.

// Assignment operator.
FallOffReaction &FallOffReaction::operator=(const Sprog::Kinetics::FallOffReaction &rxn)
{
    // Check for self assignment!
    if (this != &rxn) {
        // Remember to copy base class parts.
        ThirdBodyReaction::operator=(rxn);

        // Copy fall-off specific parts.
        m_lowplimit = rxn.m_lowplimit;
        m_fothirdbody = rxn.m_fothirdbody;
    }

    return *this;
}


// LOW PRESSURE LIMIT.

// Returns the Arrhenius parameters for the low pressure limit.
const ARRHENIUS &FallOffReaction::LowPressureLimit() const
{
    return m_lowplimit;
}

// Sets the low pressure limit Arrhenius parameters.
void FallOffReaction::SetLowPressureLimit(const Sprog::Kinetics::ARRHENIUS &lowp)
{
    m_lowplimit = lowp;
}


// FALL-OFF THIRD BODY.

// Returns a pointer to the species used as a third body in the fall-off reaction.
const Species &FallOffReaction::FallOffThirdBody() const
{
    return *(*m_species)[m_fothirdbody];
}

// Sets the species to use as a third body for fall-off calculations.
void FallOffReaction::SetFallOffThirdBody(unsigned int sp)
{
    m_fothirdbody = sp;
}

// Sets the species to use as a third body for fall-off calculations
// given the species name.
void FallOffReaction::SetFallOffThirdBody(const std::string &name)
{
    // Locate the species by name in the vector.
    unsigned int i;
    for (i=0; i<m_species->size(); i++) {
        if (*(*m_species)[i] == name) {
            // Found the species!
            m_fothirdbody = i;
            return;
        }
    }
}

// CLONING.

// Returns a pointer to a copy of the fall-off reaction.
FallOffReaction* FallOffReaction::Clone(void) const
{
    return new FallOffReaction(*this);
}