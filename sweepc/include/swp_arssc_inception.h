/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of an inception process which is also aware of the 
    ARS-SC model. This allows the inception to create particles with
    a number of aromatic surface sites.
*/

#ifndef SWEEP_ARSSC_INCEPTION_H
#define SWEEP_ARSSC_INCEPTION_H

#include "swp_params.h"
#include "swp_process_type.h"
#include "swp_inception.h"
#include "swp_cell.h"
#include "swp_mechanism.h"
#include <vector>
#include <iostream>

namespace Sweep
{
namespace Processes
{
class ARSSC_Inception : public Inception
{
public: 
    // Constructors.
    ARSSC_Inception(const Sweep::Mechanism &mech); // Initialising constructor.
    ARSSC_Inception(const ARSSC_Inception &copy);  // Copy constructor.
    ARSSC_Inception(                 // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    // Destructors.
    ~ARSSC_Inception(void);

    // Operators.
    ARSSC_Inception &operator=(const ARSSC_Inception &rhs);


	// PERFORMING THE PROCESS.

    // Performs the process on the given system.  The responsible rate term is given
    // by index.  Returns 0 on success, otherwise negative.
    int Perform(
        real t,                // Time.
        Cell &sys,             // System to update.
        unsigned int iterm = 0 // The process term responsible for this event.
        ) const;


    // READ/WRITE/COPY.

    // Creates a copy of the inception.
    ARSSC_Inception *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    ProcessType ID(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    // Aromatic site counts of newly incepted particle.
    fvector m_sites; 

    // Default constructor is protected to prevent an inception being
    // defined without knowledge of the parent mechanism.
    ARSSC_Inception(void);
};
};
};

#endif
