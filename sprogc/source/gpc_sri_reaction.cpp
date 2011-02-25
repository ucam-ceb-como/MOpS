#include "gpc_sri_reaction.h"

using namespace Sprog;
using namespace Sprog::Kinetics;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
SRIReaction::SRIReaction()
{
}

// Copy constructor.
SRIReaction::SRIReaction(const Sprog::Kinetics::SRIReaction &rxn)
{
    *this = rxn;
}

// Destructor.
SRIReaction::~SRIReaction()
{
}


// OPERATOR OVERLOADING.

// Assignment operator.
Kinetics::SRIReaction &SRIReaction::operator=(const Sprog::Kinetics::SRIReaction &rxn)
{
    // Check for self assignment!
    if (this != &rxn) {
        // Remember to assign base class.
        FallOffReaction::operator=(rxn);

        // Copy SRI parameters.
        m_params = rxn.m_params;
    }
    return *this;
}


// SRI PARAMETERS.

// Returns the SRI parameters.
const Kinetics::SRIReaction::SRI_PARAMS &SRIReaction::FallOffParams() const
{
    return m_params;
}

// Sets the SRI fall-off parameters.
void SRIReaction::SetFallOffParams(const Sprog::Kinetics::SRIReaction::SRI_PARAMS &params)
{
    m_params = params;
}


// CLONING.

// Returns a pointer to a copy of the SRI reaction.
Kinetics::SRIReaction* SRIReaction::Clone(void) const
{
    return new SRIReaction(*this);
}