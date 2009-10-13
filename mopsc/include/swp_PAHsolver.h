/*
  Author(s):      Markus Sander (ms785)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2009 Markus Sander.

  File purpose:


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

#ifndef SWEEP_PAH_SOLVER_H
#define SWEEP_PAH_SOLVER_H

#include "mops_particle_solver.h"
#include "swp_gas_profile.h"
#include "sweep.h"
#include "swp_PAH_primary.h"
#include "swp_PAH_trajectory.h"
#include "swp_PAH_database.h"
#include <map>

namespace Sweep
{
class PAHSolver : public Mops::ParticleSolver, public Sweep::Solver
{
public:
/*
    // A map of Time/Gas-Phase pairs which describes a gas-phase
    // chemistry profile.
    struct GasPoint {
        real Time;
        Sprog::Thermo::IdealGas Gas;
        GasPoint(const Sprog::SpeciesPtrVector &species);
    };
    typedef std::vector<GasPoint> GasProfile;
    typedef std::map<real, Sprog::Thermo::IdealGas> GasProfile;
    typedef std::pair<real, Sprog::Thermo::IdealGas> GasPoint;
*/

    // Constructors.
    PAHSolver(void); // Default constructor.

    // Destructors.
    virtual ~PAHSolver(void); // Default destructor.

    void Solve(
        Mops::Reactor &r,           // The reactor to solve.
        real tstop,                 // The end time for the step.
        int nsteps,                 // Number of internal steps to take.
        int niter,                  // Number of internal iterations to take.
        Mops::Solver::OutFnPtr out, // Output function pointer.
        void *data                  // Custom data object which will be passed as argument to out().
        );

	 void Output(Mops::Reactor &r);
	 void linInterpGas(Sweep::real t, 
                               const GasProfile &gasphase, 
                               Sprog::Thermo::IdealGas &gas) const;

	 void LoadGasProfile(const std::string &file, Mops::Mechanism &mech);
	 void LoadPAHProfile(const std::string &file);
	 void UpdateFractalDimension(Mops::Reactor &r);


private:
    // The gas-phase chemistry profile.
    GasProfile m_gasprof;
	// the id of the PAH. 
	//int ID;
/*	int maxID;
	struct trajectory {
		std::vector<double> time;
		std::vector<int> n_carbon_t;
    };
	std::vector<trajectory> alltrajectories;*/
    void UpdateNumberPAHs(Mops::Reactor &r,real time); 
    void UpdatePAHs(Mops::Reactor &r, real t);
	//std::vector<Trajectory>  m_trajectory;
   // Trajectory  m_trajectory;
    PAH_database m_PAHDatabase;
};
};

#endif
