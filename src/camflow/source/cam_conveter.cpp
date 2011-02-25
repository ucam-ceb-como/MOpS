/*
 * File:   cam_converter.cpp
 * Author: vinod
 *
 * Created on January 17, 2009, 2:40 PM
 ** Copyright (C) 2008 Vinod M Janardhanan.

 * File purpose:
 *  Unit conversion
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
#include "cam_converter.h"

#include "cam_conc.h"
#include <sstream>
/*
 *
 */
using namespace Camflow;
using namespace Strings;

Camflow::doublereal CamConverter::getConvertionFactor(std::string unit){
    std::string unitName = convertToCaps(unit);
    if(!unitName.compare("CM") ){
        factor = 0.01;
    }else if(!unitName.compare("IN")){
        factor = 0.0254;
    }else if(!unitName.compare("M")){
        factor = 1.0;
    }else if(!unitName.compare("K")){
        factor = 0.0;
    }else if(!unitName.compare("C")){
        factor = 273.15;
    }else if(!unitName.compare("ATM")){
        factor = 101325.0;
    }else if(!unitName.compare("PA")){
        factor = 1.0;
    }else if(!unitName.compare("BAR")){
        factor = 100000.0;
    }else if(!unitName.compare("CM/S")){
        factor = 0.01;
    }else if(!unitName.compare("M/S")){
        factor = 1.0;
    }else if(!unitName.compare("SI")){
        factor = 1.0;
    }else if(!unitName.compare("CGS")){
        factor = 10.0;
    }else{
        std::string err = "No conversion factor found for unit" + unit + "\n";
        throw CamError(err);
    }

    return factor;
    
}

std::string CamConverter::getString(int n){

    std::stringstream out;
    out << n;
    return out.str();
}

std::string CamConverter::getString(doublereal a){
    std::stringstream out;
    out << a;
    return out.str();
}

//convert from mole 2 mass
Camflow::doublereal CamConverter::mole2mass(std::vector<doublereal>& mole,
                std::vector<doublereal>& mass, Mechanism &mech ){

    const SpeciesPtrVector spv = mech.Species();
    doublereal avgMolWt = 0;
    int len = mole.size();
    for (int l = 0; l < len; l++) {
        avgMolWt += mole[l]* spv[l]->MolWt();
    }

    mass.resize(len,0.0);
    for (int l = 0; l < len; l++) {
        mass[l] = mole[l]* spv[l]->MolWt() /avgMolWt;
    }

    return avgMolWt;

}


Camflow::doublereal CamConverter::mass2mole(std::vector<doublereal>& mass,
                                            std::vector<doublereal>& mole,
                                            Mechanism& mech){
    const SpeciesPtrVector spv = mech.Species();
    doublereal avgMolWt = 0;
    int len = mass.size();
    for(int i=0; i<len; i++){
        avgMolWt += mass[i]/spv[i]->MolWt();
    }

    mole.resize(len,0.0);
    for(int i=0; i<len; i++){
        mole[i] = (mass[i]/spv[i]->MolWt())/avgMolWt;
    }
    avgMolWt = 1./avgMolWt;
    return avgMolWt;
}
