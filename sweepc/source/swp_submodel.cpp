#include "swp_submodel.h"
#include "swp_primary.h"

using namespace Sweep;
using namespace Sweep::SubModels;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
SubModel::SubModel(void)
{
    m_parent = NULL;
}

// Initialising constructor.
SubModel::SubModel(Primary &parent)
{
    m_parent = &parent;
}

// Copy constructor (protected).
SubModel::SubModel(const SubModel &copy)
{
    *this = copy;
}

// Default destructor.
SubModel::~SubModel(void)
{
}


// OPERATOR OVERLOADS.

// Assignment operator.
SubModel &SubModel::operator=(const SubModel &rhs)
{
    if (this != &rhs) {
        // The parent primary particle is not set here.
    }
    return *this;
}


// PARENT PRIMARY PARTICLE.

// Returns the parent primary particle.
const Primary *const SubModel::Parent(void) const
{
    return m_parent;
}

// Sets the parent primary particle.
void SubModel::SetParent(Primary &parent)
{
    m_parent = &parent;
}
