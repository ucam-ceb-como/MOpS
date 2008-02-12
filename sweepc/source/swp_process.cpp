#include "swp_process.h"
#include "swp_mechanism.h"
#include "rng.h"

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Process::Process(void)
: m_mech(NULL)
{
}

// Copy constructor.
Process::Process(const Sweep::Process &copy)
{
    *this = copy;
}

// Default destructor.
Process::~Process(void)
{
}


// OPERATOR OVERLOADS.

// Assigment operator.
Process &Process::operator=(const Sweep::Process &rhs)
{
    if (this != &rhs) {
        m_mech = rhs.m_mech;

        // Copy reactants.
        Sprog::StoichMap::const_iterator i;
        for (i=rhs.m_reac.begin(); i!=rhs.m_reac.end(); ++i) {
            m_reac[i->first] = i->second;
        }

        // Copy products.
        for (i=rhs.m_prod.begin(); i!=rhs.m_prod.end(); ++i) {
            m_prod[i->first] = i->second;
        }
    }
    return *this;
}


// PARENT MECHANISM.

// Returns reference to parent mechanism.
const Sweep::Mechanism *const Process::Mechanism() const
{
    return m_mech;
}

// Sets the parent mechanism
void Process::SetMechanism(const Sweep::Mechanism &mech)
{
    m_mech = &mech;
}


// REACTANTS.

// Returns the reactant count.
unsigned int Process::ReactantCount() const
{
    return m_reac.size();
}

// Returns the stoichiometric reactant coefficients.
const Sprog::StoichMap &Process::Reactants(void) const
{
    return m_reac;
}

// Returns the stoichiometry of the ith reactant.
int Process::Reactants(unsigned int i) const
{
    Sprog::StoichMap::const_iterator j = m_reac.find(i);
    if (j != m_reac.end()) {
        return j->second;
    } else {
        return 0;
    }
}

// Adds a reactant to the reaction.
void Process::AddReactant(unsigned int isp, int mu)
{
    m_reac[isp] = mu;
}

// Adds a reactant given the species name.
void Process::AddReactant(const std::string &name, int mu)
{
    // Locate the species.
    unsigned int i;
    for (i=0; i!=m_mech->Species()->size(); ++i) {
        if (name.compare((*m_mech->Species())[i]->Name()) == 0) {
            // Set reactant.
            m_reac[i] = mu;
            return;
        }
    }
}

// Removes a reactant, given by name, from the reaction.
void Process::RemoveReactant(const std::string &name)
{
    // Locate the species.
    unsigned int i;
    for (i=0; i!=m_mech->Species()->size(); ++i) {
        if (name.compare((*m_mech->Species())[i]->Name()) == 0) {
            // Delete reactant.
            m_reac.erase(i);
            return;
        }
    }
}


// PRODUCTS.

// Returns the product count.
unsigned int Process::ProductCount() const
{
    return m_prod.size();
}

// Returns the stoichiometric product coefficients.
const Sprog::StoichMap &Process::Products(void) const
{
    return m_prod;
}

// Returns the stoichiometry of the ith product.
int Process::Products(unsigned int i) const
{
    Sprog::StoichMap::const_iterator j = m_prod.find(i);
    if (j != m_prod.end()) {
        return j->second;
    } else {
        return 0;
    }
}

// Adds a product to the reaction.
void Process::AddProduct(unsigned int isp, int mu)
{
    m_prod[isp] = mu;
}

// Adds a product given the species name.
void Process::AddProduct(const std::string &name, int mu)
{
    // Locate the species.
    unsigned int i;
    for (i=0; i!=m_mech->Species()->size(); ++i) {
        if (name.compare((*m_mech->Species())[i]->Name()) == 0) {
            // Set product.
            m_prod[i] = mu;
            return;
        }
    }
}

// Removes a product, given by name, from the reaction.
void Process::RemoveProduct(const std::string &name)
{
    // Locate the species.
    unsigned int i;
    for (i=0; i!=m_mech->Species()->size(); ++i) {
        if (name.compare((*m_mech->Species())[i]->Name()) == 0) {
            // Delete product.
            m_prod.erase(i);
            return;
        }
    }
}


// FICTICIOUS EVENTS.

// Determines whether a rate is ficticious given 
// the majorant and true values.
bool Process::Ficticious(real majk, real truek)
{
    return !((majk*rnd()) < truek);
}


// PROTECTED HELPER FUNCTIONS.

// Adjusts the gas-phase composition using the reactants and
// products defined for this process.
void Process::adjustGas(Cell &sys, unsigned int n) const
{
    fvector dc(m_mech->Species()->size(), 0.0);
    Sprog::StoichMap::const_iterator i;
    real n_NAvol = (real)n / (NA * sys.SampleVolume());
    for (i=m_reac.begin(); i!=m_reac.end(); ++i)
        dc[i->first] -= (real)(i->second) * n_NAvol;
    for (i=m_prod.begin(); i!=m_prod.end(); ++i)
        dc[i->first] +=(real)(i->second) * n_NAvol;
    sys.AdjustConcs(dc);
}
