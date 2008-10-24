#ifndef FL_SOLVER_MANAGER_H
#define FL_SOLVER_MANAGER_H
#include "gpc.h"
#include "fl_params.h"
#include "fl_reactor.h"
#include "fl_solver_control.h"
//#include "fl_io.h"
#include "nvector/nvector_serial.h"
#include "cvodes_impl.h" // For CVodeMem.
#include "cvodes_band_impl.h" // For band mat.
#include "cvodes_dense_impl.h"
#include <vector>
namespace FlameLab{
	class FlameLabIO;
	class SolverManager{

	protected:
		void *cvode_mem; // memory space
		real *ptrToSlnVector; // pointer to solution vector
		real tMax,timeStep, currentTime;
		

		N_Vector solVect; // solution vector
		N_Vector derivative;
		vector<real> variables;
		vector<real> mcVariables;//variables for multi cell

	public:

		virtual void initSolver(SolverControl &sc,Reactor &reac)=0;
		virtual void solve(Sprog::Mechanism &mech, SolverControl &sc, Reactor &reac, FlameLabIO &io)=0;
		

	};
}
#endif