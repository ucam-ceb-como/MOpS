#include "swp_model.h"

using namespace Sweep;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
IModel::IModel(void)
{
}

// Copy constructor (protected).
IModel::IModel(const IModel &copy)
{
}

// Default destructor.
IModel::~IModel(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
const IModel &IModel::operator=(const IModel &rhs)
{
    return *this;
}
