
#include <map>
#include <vector>

/*
 * File:   cam_profile.h
 * Author: vinod
 *
 * Created on January 18, 2009, 10:43 AM
 *File purpose:
 *  This class contains the implementation of initial profiles
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
#include "cam_profile.h"
#include "cam_plug.h"
#include "comostrings.h"
using namespace Camflow;
using namespace Strings;

void CamProfile::setSpecies(map<string,doublereal> spec){
    iniSpec = spec;
}

//return the species initial guesses
vector<doublereal>& CamProfile::getInitialSpeciesGuess(Mechanism& mech){
    int index;
    initialfrac.resize(mech.SpeciesCount(),0.0);
    vector<doublereal> temp;
    temp.resize(mech.SpeciesCount(),0.0);
    map<string,doublereal>::iterator p;
    p = iniSpec.begin();
    while(p!=iniSpec.end()){
        index = mech.FindSpecies(convertToCaps(trim(p->first)));
        if(index < 0)
           throw CamError("Species "+p->first +" not found in species list\n");
        else
            temp[index] = p->second;
        p++;
    }

    if(getFracType() == MASS){
        initialfrac = temp;
    }else{
        CamConverter cc;
        cc.mole2mass(temp,initialfrac,mech);
    }

    return initialfrac;

}

void CamProfile::setUserTemp(doublereal pos, doublereal temp){
    u_pos.push_back(pos);
    u_temp.push_back(temp);
    
   
}

//return the user defined temperature
doublereal CamProfile::getUserDefTemp(const doublereal& pos){


    doublereal tu, tl, xu, xl, temp=0.0;
    int len = u_pos.size();



    for(int i=0; i<len; i++){
        if(pos == u_pos[i]) {            
            temp= u_temp[i];
            break;
        }else if( i>0 && (pos > u_pos[i-1]) && (pos < u_pos[i]) ){
            tu = u_temp[i];
            xu = u_pos[i];
            //cout << "location " << i << " pos " << pos << endl;
            tl = u_temp[i-1];
            xl = u_pos[i-1];
            //cout << tu << "  " << xu << endl;
            //cout << tl << "  " << xl << endl;
            doublereal slope = (tu-tl)/(xu-xl);
            doublereal intersect = tu- (slope*xu);
            //cout << "slope " << slope << endl;
            //cout << "intersect " << intersect << endl;
            temp= slope*pos + intersect;
            break;
        }
    }

    return temp;
}

vector<doublereal>& CamProfile::getPosition(){
    return this->u_pos;
}
