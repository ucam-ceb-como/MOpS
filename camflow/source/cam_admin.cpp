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

/*
 *set the reactor pressure
 */
void CamAdmin::setPressure(doublereal p_){
    this->pre = p_;
}



/*
 *set the wall temperature
 */

void CamAdmin::setWallTemp(doublereal Tw){
    T_wall = Tw;
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
void CamAdmin::setEnergyModel(string model){
    if(!convertToCaps(model).compare("ADIABATIC"))
        setEnergyModel(ADIABATIC);
    if(!convertToCaps(model).compare("USERDEFINED"))
        setEnergyModel(USERDEFINED);
    if(!convertToCaps(model).compare("ISOTHERMAL"))
        setEnergyModel(ISOTHERMAL);
    if(!convertToCaps(model).compare("NONISOTHERMAL"))
        setEnergyModel(NONISOTHERMAL);
}


void CamAdmin::setReportSchedule(string repo){
    if(!convertToCaps(repo).compare("FINAL"))
        setReportSchedule(FINAL);
    if(!convertToCaps(repo).compare("INTER"))
        setReportSchedule(INTER);
}


void CamAdmin::setReportSchedule(int n){
    repSch = n;
}


int CamAdmin::getReportSchedule(){
    return repSch;
}


void CamAdmin::setLeftBoundary(CamBoundary &cb){    
    this->left = cb;
}

void CamAdmin::setRightBoundary(CamBoundary &cb){
    this->right = cb;
}


void CamAdmin::setSpeciesOut(int n){
    speciesOut = n;
}

int CamAdmin::getSpeciesOut(){
    return speciesOut;
}

doublereal CamAdmin::getPressure() const{
    return this->pre;
}


//doublereal CamAdmin::getTemperature() const{
//    return this->T;
//}

doublereal CamAdmin::getWallTemp() const{
    return this->T_wall;
}

int CamAdmin::getEnergyModel() const{
    return this->energyModel;
}


void CamAdmin::getLeftBoundary(CamBoundary &cb){    
    cb = left;
}

void CamAdmin::getRightBoundary(CamBoundary &cb){
    cb = right;
}

doublereal CamAdmin::getNre(const doublereal& hd,
                            const doublereal& u,
                            const doublereal& rho,
                            const doublereal& eta){

    return (hd*u*rho/eta);
}

doublereal CamAdmin::getPrandtl(const doublereal& eta,
                                const doublereal& lambda,
                                const doublereal& cp){
    return (cp*eta/lambda);
}

doublereal CamAdmin::getGraetzInv(const doublereal& x,
                                    const doublereal& dh,
                                    const doublereal& Nre,
                                    const doublereal& Pr){

    return ((x+1e-10)/(dh*Nre*Pr));
}

doublereal CamAdmin::getNusselt(const doublereal& gzInv){

    //this is strictly valid for non reacting cylindrical tubes
     return  3.657+ 8.827*pow(1000.0*gzInv,-0.545)*exp(-48.2*gzInv);
}

doublereal CamAdmin::getHeatTransferCoeff(  const doublereal& x,
                                            const doublereal& vel,
                                            const doublereal& hd,
                                            const doublereal& rho,
                                            const doublereal& eta,
                                            const doublereal& lambda,
                                            const doublereal& cp){

    doublereal nre = getNre(hd,vel,rho,eta);
    doublereal pr = getPrandtl(eta,lambda,cp);
    doublereal greInv = getGraetzInv(x,hd,nre,pr);
    doublereal nusselt = getNusselt(greInv);
    return (nusselt*lambda/hd);

}