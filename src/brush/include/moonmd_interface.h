/*!
  \author      Robert I A Patterson
  \author      Sashikumaar Ganesan

  Copyright (C) 2011 Robert I A Patterson & Sashikumaar Ganesan

  \file moonmd_interface.h
  \brief Interface between continuum flow solver and stochastic particle
         population balance solver

  Licence:

    This file is free software; you can redistribute it and/or
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
    Dr. Robert Patterson
    Weierstrass Institute
    Mohrenstrasse 39
    10117 Berlin
    Germany
    riap@cantab.net
*/

#ifndef MOONMD_INTERFACE_H
#define MOONMD_INTERFACE_H

#include <vector>
#include <string>
#include <iostream>

namespace Brush {

// Forward declaration
class Reactor1d;

/*!
 * \brief Interface to provide particle stepping to the MooNMD solver
 *
 */
namespace MooNMDInterface {

//! Alias for the particle system (will eventually be Brush::Reactor1d)
typedef Brush::Reactor1d particle_reactor;

//! Pointer to a system of particles in 1 spatial dimension
typedef particle_reactor* particle_reactor_pointer;

//! Create a reactor ready for simulation
particle_reactor_pointer InitialiseBrush(const std::string& chemfile, const std::string& thermfile,
                                         const std::string& settfile, const std::string& swpfile,
                                         const std::string& partsolnfile, const std::string& chemsolnfile,
                                         const size_t num_grid_nodes, const double grid_nodes[]);

//! Run the particle phase and calculate source terms for the reacting flow solver
particle_reactor_pointer RunParticlePhase(particle_reactor& reac, const double t_stop,
    const size_t solution_length,
    const size_t num_species,
    const double solution_nodes[],
    const double temperature[],
    const double velocity[],
    const double mass_concs[],
    double energy_source[],
    double mass_conc_sources[],
    std::ostream &moment_output);


} //MooNMDInterface
} //Brush

#endif //MOONMD_INTERFACE_H
