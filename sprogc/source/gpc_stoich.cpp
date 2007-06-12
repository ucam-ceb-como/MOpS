#include "gpc_stoich.h"

using namespace Sprog;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
template<class T>
Stoichiometry<T>::Stoichiometry(void)
{
}

// Copy constructor.
template<class T>
Stoichiometry<T>::Stoichiometry(const Sprog::Stoichiometry<T> &s)
{
    m_species = s.m_species;
    m_stoich = s.m_stoich;
}

// Destructor.
template<class T>
Stoichiometry<T>::~Stoichiometry()
{
}


// OPERATOR OVERLOADING.

// Assignment operator.
template<class T>
Stoichiometry<T> &Stoichiometry<T>::operator=(const Sprog::Stoichiometry<T> &s)
{
    // Check for self-assignment!
    if (this != &s) {
        m_species = s.m_species;
        m_stoich = s.m_stoich;
    }

    return *this;
}


// SPECIES DATA.

// Returns the index of the associated species.
template<class T>
int Stoichiometry<T>::Index() const
{
    return m_species.Index;
}

// Returns a const pointer to the associated species.
template<class T>
const Sprog::Species *const Stoichiometry<T>::Species() const
{
    return m_species.Species;
}

// Sets the associated species.
template<class T>
void Stoichiometry<T>::SetSpecies(const Sprog::IndexedSpecies &sp)
{
    m_species = sp;
}


// STOICHIOMETRY VALUE.

// Returns the stoichiometry value.
template<class T>
const T &Stoichiometry<T>::Mu() const
{
    return m_stoich;
}

// Sets the stoichiometry value.
template<class T>
void Stoichiometry<T>::SetMu(const T &mu)
{
    m_stoich = mu;
}

template<class T>
void Stoichiometry<T>::IncMu(const T &mu)
{
    m_stoich += mu;
}