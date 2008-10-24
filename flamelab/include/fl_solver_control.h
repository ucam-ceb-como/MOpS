#ifndef FL_SOLVER_CONTROL_H
#define FL_SOLVER_CONTROL_H
#include "fl_params.h"
#include "fl_reactor.h"
//#include "gpc.h"

namespace FlameLab{
	class SolverControl{

		real aTol,rTol; // aobsolute and relative tolerances
		real iniStep, maxStep; // inital step and max allowed step size
		real maxTime, timeStep;
		


	public:
		enum SolutionMode{
			steadyState,
			transient
		};

		SolverControl(){}
		~SolverControl(){}
		void setATol(real aTol);
		void setRTol(real rTol);
		void setIniStep(real iniStep);
		void setMaxStep(real maxStep);
		void setMaxTime(real maxTime);
		void setTimeStep(real timeStep);
		void setSolMode(SolutionMode mode);

		real getATol()const;
		real getRTol() const;
		real getIniStep() const;
		real getMaxStep() const;
		real getMaxTime() const;
		real getTimeStep() const;
		bool operator ==(SolutionMode sm) ;



		//void initSolver(Reactor &reac, Sprog::Mechanism &mech, Sprog::Thermo::Mixture &mix);

	protected:
		SolutionMode solMode;

	
	};
}

#endif