/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of a condensation process which includes terms
    for ARS-SC sites.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#ifndef SWEEP_ARSSC_CONDENSATION_H
#define SWEEP_ARSSC_CONDENSATION_H

#include "swp_params.h"
#include "swp_process_type.h"
#include "swp_condensation.h"
#include "swp_arssc_model.h"
#include "swp_arssc_process.h"
#include "swp_cell.h"
#include "swp_particle.h"
#include <iostream>

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Transport {
    struct TransportOutflow;
}

namespace Processes
{
class ARSSC_Condensation : public Condensation, public ARSSC_Process
{
public:
    // Constructors.
    ARSSC_Condensation(const Sweep::Mechanism &mech); // Default constructor.
    ARSSC_Condensation(const ARSSC_Condensation &copy);   // Copy constructor.
    ARSSC_Condensation(              // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    // Destructor.
    virtual ~ARSSC_Condensation(void);

    // Operators.
    ARSSC_Condensation &operator=(const ARSSC_Condensation &rhs);


    // PERFORMING THE PROCESS.

    //! Perform an arssc condensation event on a particle from the ensemble
    virtual int Perform(
        real t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        int (*rand_int)(int, int), 
        real(*rand_u01)(),
        Transport::TransportOutflow *out = 0
        ) const;

    // Performs the process on a given particle in the system.  Particle
    // is given by index.  The process is performed n times.
    virtual int Perform(
        real t,        // Current time (s).
        Cell &sys,     // System to which the particle belongs.
        Particle &sp,  // Particle for which to perform process.
        unsigned int n // Number of times to perform the process.
        ) const;


    // READ/WRITE/COPY.
    
    // Creates a copy of the particle process.
    virtual ARSSC_Condensation *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    virtual ProcessType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    // Default constructor is protected to prevent reactions being
    // defined without knowledge of the parent mechanism.
    ARSSC_Condensation(void);

    // Adjusts a primary particle according to the rules of the condensation.
    unsigned int adjustPri(
        Sweep::Primary &pri, // Primary to adjust.
        unsigned int n=1     // Number of times to perform adjustment.
        ) const;
};
};
};

#endif
