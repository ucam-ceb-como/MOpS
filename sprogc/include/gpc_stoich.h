/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).

  File purpose:
    This files contains data structures which describe the stoichiometry of chemical
    reactions.
*/

#ifndef GPC_STOICH_H
#define GPC_STOICH_H

#include "gpc_species.h"
#include "gpc_params.h"

namespace Sprog
{
    // Structure to hold species stoichiometry for a reaction.
    template<class T>
    class Stoichiometry
    {
    public:
        typedef T stoich_val; // Stoich value type.

        // Constructors.
        Stoichiometry(void); // Default constructor.
        Stoichiometry(const Stoichiometry<T> &s); // Copy constructor.
        Stoichiometry(unsigned int isp, const stoich_val &mu); // Initialising constructor.

        // Destructor.
        ~Stoichiometry(void);

        // Operator overloads.
        Stoichiometry<T> &operator=(const Stoichiometry<T> &s);

        // Species data.
        int Index(void) const; // Returns index of species.
        void SetSpecies(const unsigned int &sp); // Sets the species associated with this stoichiometry.

        // Stoichiometry value.
        const stoich_val &Mu(void) const; // Returns the stoichiometry value.
        void SetMu(const stoich_val &mu); // Sets the stoichiometry value.
        void IncMu(const stoich_val &mu); // Increments the stoichiometry value.
    private:
        // Data.
        int m_species;
        stoich_val m_stoich;
    };

    // typedefs for basic stoichiometry types.
    typedef Stoichiometry<int> Stoich;   // Integer stoichiometry data.
    typedef Stoichiometry<real> Stoichf; // Real stoichiometry data.
};

#endif