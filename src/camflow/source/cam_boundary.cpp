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

void CamBoundary::setVelocity(double vel){
    this->velocity = vel;
}

void CamBoundary::setTemperature(double T){
    this->T = T;
}


void CamBoundary::setFlowRate(double flow){
    this->flowRate = flow;
}

void CamBoundary::setSpecies(std::map<std::string,double> species){
    this->species = species;
}

void CamBoundary::setInletMassfracs(std::vector<double> fracs){
    this->inletFracs = fracs;
}

std::vector<double>& CamBoundary::setInletfracs(Mechanism& mech){
    int index;
    inletFracs.resize(mech.SpeciesCount(),0.0);
    std::map<std::string,double>::iterator p;
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

const double& CamBoundary::getVelocity() const{
    return this->velocity;
}

const double& CamBoundary::getTemperature() const{
    return this->T;
}

const double& CamBoundary::getFlowRate() const{
    return this->flowRate;
}

const std::map<std::string,double>& CamBoundary::getInletSpecies() const{
    return this->species;
}

const std::vector<double>& CamBoundary::getInletMassfracs() const{
    return this->inletFracs;
}


