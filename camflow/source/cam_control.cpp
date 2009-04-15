/*
 * File:   cam_control.cpp
 * Author: vinod
 *
 * Created on January 17, 2009, 1:00 PM
 * Copyright (C) 2008 Vinod M Janardhanan.

 * File purpose:
 *  This class can be used to control the behavior of the solver
 *  based on the controls defined in the input file
 * Licence:
 *  This file is part of "Camflow".
 *
 *  Camflow is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * Contact:
 *  Dr Markus Kraft
 *  Dept of Chemical Engineering
 *  University of Cambridge
 *  New Museum Site
 *  Pembroke Street
 *  Cambridge
 *  CB2 3RA
 *  UK
 *
 *  Email   :   mk306@cam.ac.uk
 *  Website :   http://como.cheng.cam.ac.uk
 */

#include "cam_control.h"
#include "cam_configuration.h"
using namespace Camflow;

/*
 *
 */
//default constructor implementation
CamControl::CamControl(){
    setSpeciesAbsTol(1.e-04);
    setSpeciesRelTol(1.e-06);
    setTempAbsTol(1.e-02);
    setTempRelTol(1.e-03);
    setFlowAbsTol(1.e-03);
    setFlowRelTol(1.e-03);
    setResTol(1.e-05);
    setReportInterval(CamControl::FINAL);
    setIniStep(0.0);
    setMaxTime(1e5);
}

void CamControl::setSpeciesRelTol(doublereal tol){
    if(tol != 0)this->relTolSp = tol;
}

void CamControl::setSpeciesAbsTol(doublereal tol){
    if(tol != 0) this->absTolSp = tol;
}

void CamControl::setTempRelTol(doublereal tol){
    if(tol != 0) this->relTolT = tol;
}

void CamControl::setTempAbsTol(doublereal tol){
    if(tol != 0) this->absTolT = tol;
}

void CamControl::setFlowRelTol(doublereal tol){
    if(tol != 0) this->relTolFlow = tol;
}

void CamControl::setFlowAbsTol(doublereal tol){
    if(tol != 0) this->absTolFlow = tol;
}

void CamControl::setResTol(doublereal tol){
    if(tol != 0) this->resTol = tol;
}

void CamControl::setSolutionMode(int n){
    this->solMode = n;
}

void CamControl::setReportInterval(int n){
    this->repotMode = n;
}

void CamControl::setIniStep(doublereal step){
    this->iniStep = step;
}

void CamControl::setMinStep(doublereal step){
    this->minStep = step;
}

void CamControl::setMaxStep(doublereal step){
    this->maxStep = step;
}

void CamControl::setMaxTime(doublereal time){
    this->maxTime = time;
}


void CamControl::setSolver(int n){
    solver = n;
}

doublereal CamControl::getSpeciesRelTol() const{   
    return this->relTolSp;
}

doublereal CamControl::getSpeciesAbsTol() const{
    return this->absTolSp;
}

doublereal CamControl::getTempRelTol() const{
    return this->relTolT;
}

doublereal CamControl::getTempAbsTol() const{
    return this->absTolT;
}

doublereal CamControl::getFlowRelTol() const{
    return this->relTolFlow;
}

doublereal CamControl::getFlowAbsTol() const{
    return this->absTolFlow;
}

doublereal CamControl::getResTol() const{
    return this->resTol;
}

int CamControl::getSolutionMode() const{
    return this->solMode;
}

int CamControl::getReportInterval() const{
    return this->repotMode;
}

int CamControl::getSolver() const{
    return this->solver;
}

doublereal CamControl::getIniStep() const{
    return this->iniStep;
}

doublereal CamControl::getMinStep() const{
    return this->minStep;
}

doublereal CamControl::getMaxStep() const{
    return this->maxStep;
}

doublereal CamControl::getMaxTime() const{
    return this->maxTime;
}
