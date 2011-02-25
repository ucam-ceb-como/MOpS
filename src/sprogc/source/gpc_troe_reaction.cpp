#include "gpc_troe_reaction.h"

using namespace Sprog;
using namespace Sprog::Kinetics;
using namespace std;

// THE TroeReaction CLASS.

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
TroeReaction::TroeReaction()
{
}

// Copy constructor.
TroeReaction::TroeReaction(const Sprog::Kinetics::TroeReaction &rxn)
{
    *this = rxn;
}

// Destructor.
TroeReaction::~TroeReaction()
{
}


// OPERATOR OVERLOADING.

Kinetics::TroeReaction &TroeReaction::operator=(const Sprog::Kinetics::TroeReaction &rxn)
{
    // Check for self assignment.
    if (this != &rxn) {
        // Remember to copy base class.
        FallOffReaction::operator=(rxn);

        // Copy Troe parameters.
        m_params = rxn.m_params;
    }
    return *this;
}


// TROE PARAMETERS.

// Returns the Troe fall-off parameters.
const TroeReaction::TroeParams &TroeReaction::FallOffParams() const
{
    return m_params;
}

// Sets the Troe fall-off parameters.
void TroeReaction::SetFallOffParams(const Sprog::Kinetics::TroeReaction::TroeParams &params)
{
    m_params = params;
}


// CLONING.

// Returns a pointer to a copy of the Troe reaction.
TroeReaction* TroeReaction::Clone(void) const
{
    return new TroeReaction(*this);
}


// THE TroeParams STRUCT.

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
TroeReaction::TroeParams::TroeParams()
{
    Alpha = Ts = Tsss = 0.0;
    Tss = NULL;
}

// Copy constructor.
TroeReaction::TroeParams::TroeParams(const Sprog::Kinetics::TroeReaction::TroeParams &troe)
{
    Alpha = troe.Alpha;
    Ts = troe.Ts;
    Tsss = troe.Tsss;

    // Copy Tss if it isn't NULL.
    if (troe.Tss != NULL) {
        Tss = new real(*troe.Tss);
    } else {
        Tss = NULL;
    }
}

// Destructor.
TroeReaction::TroeParams::~TroeParams()
{
    // Clear memory.
    if (Tss != NULL) delete Tss;
}


// OPERATOR OVERLOADING.

// Assignment operator.
Sprog::Kinetics::TroeReaction::TroeParams &TroeReaction::TroeParams::operator=(const Sprog::Kinetics::TroeReaction::TroeParams &troe)
{
    // Check for self assignment!
    if (this != &troe) {
        Alpha = troe.Alpha;
        Ts = troe.Ts;
        Tsss = troe.Tsss;

        // Copy Tss if it isn't NULL.
        if (troe.Tss != NULL) {
            Tss = new real(*troe.Tss);
        } else {
            Tss = NULL;
        }
    }
    return *this;
}