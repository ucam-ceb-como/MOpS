/*
 * File:   cam_admin.cpp
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2009 Vinod M Janardhanan
 *
 * File purpose:
 *  This class contains the implementation of process conditions
 *  and the boundary conditions
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
 *
 * Created on January 17, 2009, 4:36 PM
 */

#include <string>


#include "cam_admin.h"
#include <cmath>
using namespace Strings;
using namespace Camflow;


bool CamAdmin::radiation = false;

/*
 *set the reactor pressure
 */
void CamAdmin::setPressure(double p_){
    this->pre = p_;
}

/**
 *  Set the ignition step for temperature
 */
void CamAdmin::setIgnitionStep(double step){
    this->stepIgnition = step;
}

/*
 *set the wall temperature
 */

void CamAdmin::setWallTemp(double Tw){
    T_wall = Tw;
}

void CamAdmin::setFlameletEquationType(const std::string type)
{
    if (type == "SIMPLE")
    {
        std::cout << "Flamelet Equation Type is " << type << std::endl;
        flameletEquationType_ = SIMPLE;
    }
    else if (type == "COMPLETE")
    {
        std::cout << "Flamelet Equation Type is " << type << std::endl;
        flameletEquationType_ = COMPLETE;
    }
    else
    {
        throw std::invalid_argument(type+" is not a valid flameletEquation.");
    }
}

/*
 *set the energy model ISOTHERMAL, ADIABATIC, USERDEFINED
 */
void CamAdmin::setEnergyModel(int n){
    this->energyModel = n;
}

/*
 *set the energy model ISOTHERMAL, ADIABATIC, USERDEFINED
 */
void CamAdmin::setEnergyModel(std::string model){
    if(!convertToCaps(model).compare("ADIABATIC"))
        setEnergyModel(ADIABATIC);
    if(!convertToCaps(model).compare("USERDEFINED"))
        setEnergyModel(USERDEFINED);
    if(!convertToCaps(model).compare("ISOTHERMAL"))
        setEnergyModel(ISOTHERMAL);
    if(!convertToCaps(model).compare("NONISOTHERMAL"))
        setEnergyModel(NONISOTHERMAL);
}


double CamAdmin::getIgnitionStep() const{
    return stepIgnition;
}

void CamAdmin::setLeftBoundary(CamBoundary &cb){
    this->left = cb;
}

void CamAdmin::setRightBoundary(CamBoundary &cb){
    this->right = cb;
}

void CamAdmin::setRadiationModel(bool radiation_){
    radiation = radiation_;
}


void CamAdmin::setSpeciesOut(int n){
    speciesOut = n;
}

void CamAdmin::setRestartType(const std::string& restartType)
{
    if (!convertToCaps(restartType).compare("BINARY"))
        restartType_ = BINARY;
    else if (!convertToCaps(restartType).compare("TEXT"))
        restartType_ = TEXT;
    else if (!convertToCaps(restartType).compare("NONE"))
        restartType_ = NONE;
    else throw std::runtime_error("What restart file type?! Check:"
                                  "<restart file=""></restart>.");
}

void CamAdmin::setRestartFile(const std::string& restartFile)
{
    restartFile_ = restartFile;
}

void CamAdmin::setInputFile(std::string inputFileName)
{
    inputFileName_ = inputFileName;
}

const std::string& CamAdmin::getInputFile() const
{
    return inputFileName_;
}

int CamAdmin::getSpeciesOut() const{
    return speciesOut;
}

double CamAdmin::getPressure() const{
    return this->pre;
}

//double CamAdmin::getTemperature() const{
//    return this->T;
//}

double CamAdmin::getWallTemp() const{
    return this->T_wall;
}

int CamAdmin::getFlameletEquationType() const
{
    return flameletEquationType_;
}

int CamAdmin::getEnergyModel() const{
    return this->energyModel;
}


Camflow::CamBoundary& CamAdmin::getLeftBoundary(){
    return left;
}

Camflow::CamBoundary& CamAdmin::getRightBoundary(){
    return right;
}

bool CamAdmin::getRadiationModel() const{
    return radiation;
}

int CamAdmin::getRestartType() const
{
    return restartType_;
}

const std::string& CamAdmin::getRestartFile() const
{
    return restartFile_;
}

double CamAdmin::getNre(const double& hd,
                            const double& u,
                            const double& rho,
                            const double& eta){

    return (hd*u*rho/eta);
}

double CamAdmin::getPrandtl(const double& eta,
                                const double& lambda,
                                const double& cp){
    return (cp*eta/lambda);
}

double CamAdmin::getGraetzInv(const double& x,
                                    const double& dh,
                                    const double& Nre,
                                    const double& Pr){

    return ((x+1e-10)/(dh*Nre*Pr));
}

double CamAdmin::getNusselt(const double& gzInv){

    //this is strictly valid for non reacting cylindrical tubes
     return  3.657+ 8.827*pow(1000.0*gzInv,-0.545)*exp(-48.2*gzInv);
}

double CamAdmin::getHeatTransferCoeff(  const double& x,
                                            const double& vel,
                                            const double& hd,
                                            const double& rho,
                                            const double& eta,
                                            const double& lambda,
                                            const double& cp){

    double nre = getNre(hd,vel,rho,eta);
    double pr = getPrandtl(eta,lambda,cp);
    double greInv = getGraetzInv(x,hd,nre,pr);
    double nusselt = getNusselt(greInv);
    return (nusselt*lambda/hd);

}
