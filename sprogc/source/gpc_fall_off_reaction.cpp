#include "gpc_fall_off_reaction.h"
#include "gpc_third_body_reaction.h"
#include "gpc_rate_params.h"

using namespace Sprog;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
FallOffReaction::FallOffReaction()
{
}

// Copy constructor.
FallOffReaction::FallOffReaction(const Sprog::FallOffReaction &rxn)
{
    *this = rxn;
}

// Destructor.
FallOffReaction::~FallOffReaction()
{
}


// OPERATOR OVERLOADING.

// Assignment operator.
FallOffReaction &FallOffReaction::operator=(const Sprog::FallOffReaction &rxn)
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
void FallOffReaction::SetLowPressureLimit(const Sprog::ARRHENIUS &lowp)
{
    m_lowplimit = lowp;
}


// FALL-OFF THIRD BODY.

// Returns a pointer to the species used as a third body in the fall-off reaction.
const Species *const FallOffReaction::FallOffThirdBody() const
{
    return m_fothirdbody.Pointer;
}

// Sets the species to use as a third body for fall-off calculations.
void FallOffReaction::SetFallOffThirdBody(const Sprog::IndexedSpecies &sp)
{
    m_fothirdbody = sp;
}


// CLONING.

// Returns a pointer to a copy of the fall-off reaction.
FallOffReaction* FallOffReaction::Clone(void) const
{
    return new FallOffReaction(*this);
}