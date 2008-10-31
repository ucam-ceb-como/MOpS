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
		//std::vector<Sprog::Thermo::Mixture> cstrMix;
		std::vector<SingleCell> cells;
		
	public:
		const int MFLX;
		const int TEMP;
		const int VEL;
		const int DENS;
		
		Premix(const Sprog::Mechanism &mech);
		~Premix(){}

		
		
		void solve(Sprog::Mechanism &mech,SolverControl &sc, Reactor &reac, FlameLabIO &io);

		void initSolver(SolverControl &sc, Reactor &reac);
		void reInitSolver(SolverControl &sc);
		void updateVariables(FlameLab::real *y, void *object);
		
		//residual function for premix/plug flow reactor
		static int residual(double time, // present integration time
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
		real getFlux(int &cellId);

	};
}

#endif