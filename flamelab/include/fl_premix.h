/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This is the premix class. This class manages all functions related
	to premix problem solutions. This class contains implemention of 
	solver call for both steady state and transient problems. The residual
	functions are also defined in this class. 
  Licence:
    This file is part of "flameLab".

    flameLab is free software; you can redistribute it and/or
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
#ifndef FL_PREMIX_H
#define FL_PREMIX_H
#include "fl_params.h"
#include "fl_reactor.h"
#include "fl_io.h"
#include "fl_solver_manager.h"
#include "fl_solver_control.h"
#include "fl_single_cell.h"
#include <cmath>
#include <map>
namespace FlameLab{	
	class Premix : public SolverManager, SingleCell{
		int nEq, nVar, nSpecies;
		real aTol;
		Reactor *ptrToReactor;
		
		std::vector<SingleCell> cells;
		
	public:
		const int MFLX;
		const int TEMP;
		const int VEL;
		const int DENS;
		real heatConvection, heatSource;

		Premix(const Sprog::Mechanism &mech);
		~Premix(){}

		
		
		void solve(Sprog::Mechanism &mech,SolverControl &sc, Reactor &reac, FlameLabIO &io);

		void initSolver(SolverControl &sc, Reactor &reac);
		void reInitSolver(SolverControl &sc, Reactor &reac);
		void updateVariables(FlameLab::real *y, void *object);
		
		//residual function for premix/plug flow reactor
		static int residual(double time, // present integration time
							N_Vector y,  // solution vector
							N_Vector ydot, // time derivatibe of soln vector
							void *object); // some usefuel object

		static int trResidual(double time, // present integration time
							N_Vector y,  // solution vector
							N_Vector ydot, // time derivatibe of soln vector
							void *object); // some usefuel object

		//mass flux residual for single cell
		static void scMassFlux(real &time, real *y, real *ydot, void *object);
		// species residual definitions for single cell
		static void scSpeciesResidual(real &time, real *y, real *ydot, void *object);
		// energy residual for single cell
		static void scEnergyResidual(real &time, real *y, real *ydot, void *object);


		//residual definition for multi cell problem (fully transient or case with diffusion)
		static void mcResidual(real &time, real *y, real *ydot, void *object);
		// returns the number of equations
		static void boundary(real &time, real *y, real *ydot, void *object);
		int getNeq();


	private:
		// initialize
		void initVariables(Reactor &reac);
		//initialize transient case
		void initTransVector(Reactor &reac);
		// save the solution for multi cell problem
		void saveSolution(real *y, Reactor &reac);
		bool checkConvergance(real *dy);


		// returns the mass flx given the cell id
		//real getFlux(int &cellId);

		static Premix *p;

	};
}

#endif