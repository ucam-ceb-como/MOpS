#include "gpc_stoich.h"

using namespace Sprog;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
template<class T>
Stoichiometry<T>::Stoichiometry(void)
{
    m_species = -1;
}

// Copy constructor.
template<class T>
Stoichiometry<T>::Stoichiometry(const Sprog::Stoichiometry<T> &s)
{
    m_species = s.m_species;
    m_stoich = s.m_stoich;
}

// Initialising constructor.
template<class T>
Stoichiometry<T>::Stoichiometry(unsigned int isp, const typename Stoichiometry<T>::stoich_val &mu)
{
    m_species = isp;
    m_stoich = mu;
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
    return m_species;
}

// Sets the associated species.
template<class T>
void Stoichiometry<T>::SetSpecies(const unsigned int &sp)
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


template class Stoichiometry<int>;
template class Stoichiometry<real>;
