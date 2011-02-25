/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This class can be used to control the behavior of the solver
	based on the controls defined in the input file
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
//set the residual tolerance
void SolverControl::setResTol(FlameLab::real rsTol){
	this->resTol = rsTol;
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
//sets the mininum allowed step
void SolverControl::setMinStep(FlameLab::real minStep){
	this->minStep = minStep;
}

// sets the time step in (s)
void SolverControl::setTimeStep(real timeStep){
	this->timeStep = timeStep;
}
//set the report mode
void SolverControl::setReportMode(int mode){
	this->rptMode = mode;
}
// returns the absolute tolerence
real SolverControl::getATol() const{
	return this->aTol;
}
// returns the relative tolerence
real SolverControl::getRTol() const{
	return this->rTol;
}
// returns the residual tolerance
real SolverControl::getResTol() const{
	return this->resTol;
}
// returns the max allowed step 0->initinite
real SolverControl::getMaxStep() const{
	return this->maxStep;
}
//return the min step
FlameLab::real SolverControl::getMinStep() const{
	return this->minStep;
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
void SolverControl::setSolMode(int mode) {
	this->solMode = mode;
}
//returns the solution mode steady/transient

//bool SolverControl::operator ==(SolverControl::SolutionMode sm){
//	return this->solMode == sm;
//}

void SolverControl::setOutInterval(std::map<real,real> interval, real dt){
	outInterval.insert(pair<map<real,real>,real>(interval,dt));
}

const int SolverControl::getIntervalCount() const{
	return outInterval.size();
}

const map<map<real,real>,real>& SolverControl::getOutputInterval() const{
	return outInterval;
}

int SolverControl::getSolMode() const{
	return this->solMode;
}
//return the report mode
int SolverControl::getReportMode() const{
	return rptMode;
}

//set the number of sweep
void SolverControl::setNSweep(int n){
	this->nSweep = n;
}

// return the number of sweep
int SolverControl::getNSweep() const{
	return this->nSweep;
}

// set the mrthod
void SolverControl::setMethod(int n){
	method = n;
}
//return the method

int SolverControl::getMethod(){
	return method;
}