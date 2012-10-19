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
    setResTol(1.e-07);
    setReportInterval(CamControl::FINAL);
    setIniStep(0.0);
    setMaxTime(1e5);
    setNumIterations(1);
    setSpeciesUnderRelax(1.0);
}

void CamControl::setSpeciesRelTol(double tol){
    if(tol != 0)this->relTolSp = tol;
}

void CamControl::setSpeciesAbsTol(double tol){
    if(tol != 0) this->absTolSp = tol;
}

void CamControl::setTempRelTol(double tol){
    if(tol != 0) this->relTolT = tol;
}

void CamControl::setTempAbsTol(double tol){
    if(tol != 0) this->absTolT = tol;
}

void CamControl::setFlowRelTol(double tol){
    if(tol != 0) this->relTolFlow = tol;
}

void CamControl::setFlowAbsTol(double tol){
    if(tol != 0) this->absTolFlow = tol;
}

void CamControl::setResTol(double tol){
    if(tol != 0) this->resTol = tol;
}

void CamControl::setSolutionMode(int n){
    this->solMode = n;
}

void CamControl::setReportInterval(int n){
    this->repotMode = n;
}

void CamControl::setIniStep(double step){
    this->iniStep = step;
}

void CamControl::setMinStep(double step){
    this->minStep = step;
}

void CamControl::setMaxStep(double step){
    this->maxStep = step;
}

void CamControl::setMaxTime(double time){
    this->maxTime = time;
}


void CamControl::setSolver(int n){
    solver = n;
}

void CamControl::setResidualMonitor(bool lopt){
    resMonitor = lopt;
}

void CamControl::setNumIterations(int n){
    nIter = n;
}

void CamControl::setSpeciesUnderRelax(double ur){
    urSpecies = ur;
}

double CamControl::getSpeciesRelTol() const{   
    return this->relTolSp;
}

double CamControl::getSpeciesAbsTol() const{
    return this->absTolSp;
}

double CamControl::getTempRelTol() const{
    return this->relTolT;
}

double CamControl::getTempAbsTol() const{
    return this->absTolT;
}

double CamControl::getFlowRelTol() const{
    return this->relTolFlow;
}

double CamControl::getFlowAbsTol() const{
    return this->absTolFlow;
}

double CamControl::getResTol() const{
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

double CamControl::getIniStep() const{
    return this->iniStep;
}

double CamControl::getMinStep() const{
    return this->minStep;
}

double CamControl::getMaxStep() const{
    return this->maxStep;
}

double CamControl::getMaxTime() const{
    return this->maxTime;
}

bool CamControl::getResidualMonitor() const{
    return resMonitor;
}

int CamControl::getNumIterations() const{
    return nIter;
}

double CamControl::getSpeciesUnderRelax() const{
    return urSpecies;
}

