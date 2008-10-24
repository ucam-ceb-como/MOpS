#ifndef FL_PREMIX_H
#define FL_PREMIX_H
#include "fl_params.h"
#include "fl_reactor.h"
#include "fl_io.h"
#include "fl_solver_manager.h"
#include "fl_solver_control.h"
#include <cmath>
namespace FlameLab{	
	class Premix : public SolverManager{
		int nEq, nVar, nSpecies;
		real aTol;
		Reactor *ptrToReactor;
		
	public:
		const int MFLX;
		const int TEMP;
		const int VEL;
		const int DENS;
		
		Premix(const Sprog::Mechanism &mech);
		~Premix(){}

		void initVariables(Reactor &reac);
		
		void solve(Sprog::Mechanism &mech,SolverControl &sc, Reactor &reac, FlameLabIO &io);

		void initSolver(SolverControl &sc, Reactor &reac);
		void reInitSolver(SolverControl &sc);
		void updateVariables(FlameLab::real *y, void *object);

		static int residual(double time, // present integration time
							N_Vector y,  // solution vector
							N_Vector ydot, // time derivatibe of soln vector
							void *object); // some usefuel object
		static void scMassFlux(real &time, real *y, real *ydot, void *object);
		static void scSpeciesResidual(real &time, real *y, real *ydot, void *object);
		static void scEnergyResidual(real &time, real *y, real *ydot, void *object);

		// save the solution for multi cell problem
		void saveSolution(real *y, Reactor &reac);
		bool checkConvergance(real *dy);

		// returns the mass fractions for the given cellid
		void getMassFracs(int &cellId, vector<real> &fracs);
		// returns the temperature given the cell id
		real getCellTemp(int &cellId);
		//returns the density diven the cell id
		real getDensity(int &cellId);
		// returns the velocity given the cell id
		real getVelocity(int &cellId);
		// returns the mass flx given the cell id
		real getFlux(int &cellId);
		// returns the number of equations
		int getNeq();

	};
}

#endif