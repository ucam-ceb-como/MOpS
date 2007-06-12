#include "gpc_third_body_reaction.h"
#include "gpc_reaction.h"
#include "gpc_stoich.h"
#include "gpc_species.h"

using namespace Sprog;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ThirdBodyReaction::ThirdBodyReaction()
{
}

// Copy constructor.
ThirdBodyReaction::ThirdBodyReaction(const Sprog::ThirdBodyReaction &rxn)
{
    *this = rxn;
}

// Destructor.
ThirdBodyReaction::~ThirdBodyReaction()
{
}


// OPERATOR OVERLOADING.

// Assignment operator.
ThirdBodyReaction &ThirdBodyReaction::operator =(const Sprog::ThirdBodyReaction &rxn)
{
    // Check for self assignment!
    if (this != &rxn) {
        // Remember to copy base class.
        Reaction::operator=(rxn);

        // Copy third bodies.
        m_thirdbodies.assign(rxn.m_thirdbodies.begin(), rxn.m_thirdbodies.end());
    }
    return *this;
}


// THIRD BODIES.

// Returns a constant reference to the vector of third bodies with their
// coefficients.
const vector<Stoichf> &ThirdBodyReaction::ThirdBodies() const
{
    return m_thirdbodies;
}

// Adds a third body to the reaction given the stoichiometric structure.
void ThirdBodyReaction::AddThirdBody(const Sprog::Stoichf &tb)
{
    m_thirdbodies.push_back(tb);
}

// Adds a third body to the reaction given the species and the
// coefficient.
void ThirdBodyReaction::AddThirdBody(const Sprog::IndexedSpecies &sp, Sprog::real coeff)
{
    // Build a Stoichf type.
    Stoichf mu;
    mu.SetSpecies(sp);
    mu.SetMu(coeff);

    // Add it to the array of third bodies.
    m_thirdbodies.push_back(mu);
}

// Removes a third body, given by name, from the reaction.  If the third
// body is not defined for this reaction then does nothing.
void ThirdBodyReaction::RemoveThirdBody(const std::string &name)
{
    // Search through the list of third bodies to find that with
    // the given name.
    vector<Stoichf>::iterator i;
    for (i=m_thirdbodies.begin(); i!=m_thirdbodies.end(); i++) {
        if (name.compare((*i).Species()->Name()) == 0) {
            // We have found the species in the list.
            m_thirdbodies.erase(i);
        }
    }
}


// CLONING.

// Returns a pointer to a copy of the third-body reaction.
ThirdBodyReaction* ThirdBodyReaction::Clone(void) const
{
    return new ThirdBodyReaction(*this);
}