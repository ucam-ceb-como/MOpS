/*
 * File:   cam_boundary.cpp
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2009 Vinod M Janardhanan
 *
 * File purpose:
 *  This class contains the implementation boundary conditions
 *
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
 * Created on January 17, 2009, 4:50 PM
 */

#include "cam_boundary.h"
#include "comostrings.h"
#include "cam_error.h"
using namespace Camflow;
using namespace Strings;
//copy cnstructor
CamBoundary::CamBoundary(const CamBoundary &cb){

    /*std::cout << "Invoking copy constructor@@@@@@@@@@@@ \n" ;
    std::cout << cb.velocity << std::endl;
    std::cout << cb.flowRate << std::endl;
    std::cout << cb.fracType << getFracType() << std::endl;
    std::cout << cb.T << std::endl;*/

    velocity = cb.velocity;
    flowRate = cb.flowRate;
    fracType = cb.fracType;
    inletFracs = cb.inletFracs;
    species = cb.species;
    T = cb.T;
}

void CamBoundary::setVelocity(doublereal vel){
    this->velocity = vel;
}

void CamBoundary::setTemperature(doublereal T){
    this->T = T;
}


void CamBoundary::setFlowRate(doublereal flow){
    this->flowRate = flow;
}

void CamBoundary::setSpecies(std::map<std::string,doublereal> species){
    this->species = species;
}

void CamBoundary::setInletMassfracs(std::vector<doublereal> fracs){
    this->inletFracs = fracs;
}

std::vector<doublereal>& CamBoundary::setInletfracs(Mechanism& mech){
    int index;
    inletFracs.resize(mech.SpeciesCount(),0.0);
    std::map<std::string,doublereal>::iterator p;
    p = species.begin();
    while(p!=species.end()){
        index = mech.FindSpecies(convertToCaps(trim(p->first)));
        if(index < 0)
            throw CamError("Species "+p->first +" not found in species list\n");
        else
            inletFracs[index] = p->second;
        p++;
    }

    return inletFracs;
}

const doublereal& CamBoundary::getVelocity() const{
    return this->velocity;
}

const doublereal& CamBoundary::getTemperature() const{
    return this->T;
}

const doublereal& CamBoundary::getFlowRate() const{
    return this->flowRate;
}

const std::map<std::string,doublereal>& CamBoundary::getInletSpecies() const{
    return this->species;
}

const std::vector<doublereal>& CamBoundary::getInletMassfracs() const{
    return this->inletFracs;
}


