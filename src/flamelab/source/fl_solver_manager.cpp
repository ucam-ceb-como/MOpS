#include "fl_solver_manager.h"
using namespace FlameLab;

// by default should initialize a premix reactor:
void FlameLab::SolverManager::initSolver(Sprog::Mechanism &mech, Sprog::Thermo::Mixture &mix, 
										 FlameLab::SolverControl &sc)
{
	
}

// By default will simulate a premix reactor
void SolverManager::solve(Sprog::Mechanism &mech, FlameLab::SolverControl &sc){
}
