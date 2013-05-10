/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The ParticleSolver class is an extension of the Solver class which
    also contains routines common to all mops solvers which run Sweep.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
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

#ifndef MOPS_PARTICLE_SOLVER_H
#define MOPS_PARTICLE_SOLVER_H

#include "mops_params.h"
#include "mops_reactor.h"
#include "mops_timeinterval.h"
#include "mops_mechanism.h"
#include "mops_solver.h"
#include "console_io.h"
#include "swp_gas_profile.h"
#include <vector>
#include <string>
#include <fstream>
#include <time.h>

namespace Mops
{
class ParticleSolver : public Solver
{
public:
    // Constructors.
    ParticleSolver(void); // Default constructor.

    //! Copy constructor
    ParticleSolver(const ParticleSolver &sol);

    //! Clone the object
    ParticleSolver *const Clone() const;

    // Destructors.
    virtual ~ParticleSolver(void); // Default destructor.


    // COMPUTATION TIME.
    
    // Returns the number of CT time variables tracked by this
    // solver type.
    virtual unsigned int CT_Count(void) const;

    // Outputs internal computation time data to the given
    // binary stream.
    virtual void OutputCT(std::ostream &out) const;

    // Adds the CT descriptions to a vector of strings.
    virtual void CT_Names(
        std::vector<std::string> &names, // Vector of CT names.
        unsigned int start=0 // Optional start index in vector.
        ) const;

protected:
    // COMPUTATION TIME.

    // Sweep computation time.
    double m_swp_ctime;

};
};

#endif
