#include "fl_solver_control.h"
using namespace FlameLab;
//sets the absolute tolerence
void SolverControl::setATol(real aTol){
	
	this->aTol = aTol;
}
//sets the relative tolerence
void SolverControl::setRTol(real rTol){
	this->rTol = rTol;
}
// sets the initial step size
void SolverControl::setIniStep(real iniStep){
	this->iniStep = iniStep;
}
// sets the maximum step size
void SolverControl::setMaxStep(real maxStep){
	this->maxStep = maxStep;
}
// sets the maximum integration time
void SolverControl::setMaxTime(real maxTime){
	this->maxTime = maxTime;
}
// sets the time step in (s)
void SolverControl::setTimeStep(real timeStep){
	this->timeStep = timeStep;
}
// returns the absolute tolerence
real SolverControl::getATol() const{
	return this->aTol;
}
// returns the relative tolerence
real SolverControl::getRTol() const{
	return this->rTol;
}
// returns the max allowed step 0->initinite
real SolverControl::getMaxStep() const{
	return this->maxStep;
}
// returns the initial step size
real SolverControl::getIniStep() const{
	return this->iniStep;
}
// retuns the maximum integration time
real SolverControl::getMaxTime() const{
	return this->maxTime;
}
//returns the step size for time integration
real SolverControl::getTimeStep() const{
	return this->timeStep;
}
//sets the solution mode such as steady or transient
void SolverControl::setSolMode(FlameLab::SolverControl::SolutionMode mode) {
	this->solMode = mode;
}
//returns the solution mode steady/transient

bool SolverControl::operator ==(SolverControl::SolutionMode sm){
	return this->solMode == sm;
}
