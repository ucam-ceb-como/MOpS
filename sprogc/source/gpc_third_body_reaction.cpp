#include "gpc_third_body_reaction.h"
#include "gpc_reaction.h"
#include "gpc_stoich.h"
#include "gpc_species.h"

using namespace Sprog;
using namespace Sprog::Kinetics;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ThirdBodyReaction::ThirdBodyReaction()
{
}

// Copy constructor.
ThirdBodyReaction::ThirdBodyReaction(const Sprog::Kinetics::ThirdBodyReaction &rxn)
{
    *this = rxn;
}

// Destructor.
ThirdBodyReaction::~ThirdBodyReaction()
{
}


// OPERATOR OVERLOADING.

// Assignment operator.
Kinetics::ThirdBodyReaction &ThirdBodyReaction::operator =(const Sprog::Kinetics::ThirdBodyReaction &rxn)
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
void ThirdBodyReaction::AddThirdBody(const unsigned int sp, Sprog::real coeff)
{
    // Add a new Stoichf to the array of third bodies.
    m_thirdbodies.push_back(Stoichf(sp, coeff));
}

// Adds a third body to the reaction given the species name.
void ThirdBodyReaction::AddThirdBody(const std::string &name, Sprog::real coeff)
{
    // Find the species in the vector with the given name.
    unsigned int i;
    for (i=0; i<m_species->size(); i++) {
        if (*(*m_species)[i] == name) {
            // Found the species!
            m_thirdbodies.push_back(Stoichf(i, coeff));
            return;
        }
    }
}

// Removes a third body, given by name, from the reaction.  If the third
// body is not defined for this reaction then does nothing.
void ThirdBodyReaction::RemoveThirdBody(const std::string &name)
{
    // Search through the list of species to find that with
    // the given name.
    unsigned int i;
    for (i=0; i<m_species->size(); i++) {
        if (*(*m_species)[i] == name) {
            // Found the species:  Loop though third bodies and find
            // that with this index.
            vector<Stoichf>::iterator j;
            for (j=m_thirdbodies.begin(); j!=m_thirdbodies.end(); j++) {
                if ((*j).Index() == i) {
                    // We have found the species in the list.
                    m_thirdbodies.erase(j);
                }
            }        
        }
    }
}


// CLONING.

// Returns a pointer to a copy of the third-body reaction.
ThirdBodyReaction* ThirdBodyReaction::Clone(void) const
{
    return new ThirdBodyReaction(*this);
}