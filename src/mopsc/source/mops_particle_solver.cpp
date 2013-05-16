/*
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ParticleSolver class declared in the
    mops_particle_solver.h header file.

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

#include "mops_particle_solver.h"
#include "mops_reactor_factory.h"
#include "csv_io.h"
#include "string_functions.h"
#include "sweep.h"
#include <vector>
#include <string>
#include <time.h>
#include <stdexcept>

using namespace Mops;
using namespace std;
using namespace Strings;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ParticleSolver::ParticleSolver(void)
: Solver(), m_swp_ctime(0.0)
{
}

// Copy constructor
ParticleSolver::ParticleSolver(const ParticleSolver &sol)
: Solver(sol),
  m_swp_ctime(sol.m_swp_ctime) {}

// Clone the object
ParticleSolver *const ParticleSolver::Clone() const {
    return new ParticleSolver(*this);
}

// Default destructor.
ParticleSolver::~ParticleSolver(void)
{
}


// COMPUTATION TIME.

// Returns the number of CT time variables tracked by this
// solver type.
unsigned int Mops::ParticleSolver::CT_Count(void) const 
{
    return Solver::CT_Count() + 1;
}

// Outputs internal computation time data to the given
// binary stream.
void Mops::ParticleSolver::OutputCT(std::ostream &out) const
{
    Solver::OutputCT(out);
    out.write((char*)&m_swp_ctime, sizeof(m_swp_ctime));
}

// Adds the CT descriptions to a vector of strings.
void Mops::ParticleSolver::CT_Names(vector<string> &names, unsigned int start) const
{
    // Resize output vector to hold names, and get iterator
    // to first insertion point.
    if (start+CT_Count() > names.size()) names.resize(start+CT_Count(), "");

    // Get Solver names.
    Solver::CT_Names(names, start);

    // Add names to output array.
    names[start+Solver::CT_Count()] = "MC (particle) Comput. Time (s)";
}
