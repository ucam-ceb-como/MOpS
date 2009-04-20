
#include <map>
#include <vector>

#include "cam_params.h"

#include "array.h"

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

/*
 *set products
 */
void CamProfile::setProductSpecies(map<string,doublereal> spec){
    list_prdt = spec;
}
/*
 *set intermediates
 */
void CamProfile::setIntermediateSpecies(map<string,doublereal> spec){
    list_intmd = spec;
}
/*
 *populate the product and intermediates
 */
void CamProfile::populateProducts(Mechanism& mech){
    if(list_prdt.size()>0)
        getmassFracs(list_prdt,mech,m_prdt);
}
void CamProfile::populateIntermdts(Mechanism& mech){
    if(list_intmd.size()>0)
        getmassFracs(list_intmd,mech,m_intmd);
}
/*
 *set the geometry object
 */
void CamProfile::setGeometryObj(CamGeometry& cg){
    geom = &cg;
}
/*
 *set mixing center
 */
void CamProfile::setMixingCenter(doublereal len){
    mCenter = len;
}
/*
 *set the mixing width
 */
void CamProfile::setMixingWidth(doublereal len){
    mWidth = len;
}
/*
 *set the start profile
 */
void CamProfile::setStarProfile(CamBoundary& cb, Mechanism& mech){

    vector<doublereal> m_in = cb.getInletMassfracs();
    vector<doublereal> position = geom->getAxpos();
    int len = position.size();

    start.resize(len,mech.SpeciesCount());

    populateIntermdts(mech);
    populateProducts(mech);

   
    if(mWidth != 0 && mCenter != 0 && m_prdt.size() != 0 && m_intmd.size() != 0){

        for(int i=0; i<len; i++){
            /*
             *sum the intermediates
             */
            doublereal sumInter=0;
            for(unsigned int l=0; l<mech.SpeciesCount();l++){
                sumInter += start(i,l);
            }
            doublereal factor = 1-sumInter;
            /*
             * scaling factor
             */
            doublereal f_prdt, f_reac;
            if(position[i] <= (mCenter-mWidth/2.0)) {
                f_prdt = 0.0;                                        
            }else{
                if(position[i] < (mCenter+mWidth/2.0)){
                    f_prdt = (1.0/mWidth)*(position[i]-mCenter)+0.5;
                }else{
                    f_prdt = 1.0;
                }
            }
            f_reac = 1-f_prdt;

            map<string, doublereal>::iterator p;
            map<string,doublereal> spec = cb.getInletSpecies();
            p = spec.begin();
            while(p!= spec.end()){
                int index = mech.FindSpecies(convertToCaps(trim(p->first)));
                start(i,index) = factor*(f_prdt*m_prdt[index]+f_reac*m_in[index]);
                p++;
            }
            p=list_prdt.begin();
            while(p!=list_prdt.end()){
                int index = mech.FindSpecies(convertToCaps(trim(p->first)));
                if(m_in[index]==0){
                    start(i,index) = factor*(f_prdt*m_prdt[index] + f_reac*m_in[index]);
                }
                p++;
            }
            
        }

    }else{
        for(int i=0; i<len; i++){
            for(unsigned int l=0; l<mech.SpeciesCount(); l++){
                start(i,l) = m_in[l];
            }
        }

    }

}
/*
 *set the gaussian
 */
void CamProfile::setGaussian(Mechanism& mech){
    vector<doublereal> position = geom->getAxpos();
    int len = position.size();
    map<string, doublereal>::iterator p;
    p = list_intmd.begin();
    while(p!=list_intmd.end()){
        int index = mech.FindSpecies(convertToCaps(trim(p->first)));
        doublereal gWidth = -log(0.15*m_intmd[index])/pow(mWidth/2.0,2);
        for(int i=0; i<len; i++){
            start(i,index) = m_intmd[index]*exp(-gWidth*pow(position[i]-mCenter,2));
        }

    }

}
/*
 *return the array
 */
Array2D& CamProfile::getStartProfile(){
    return start;
}
//return the species initial guesses
void CamProfile::getmassFracs(map<string,doublereal>& spec, Mechanism& mech, vector<doublereal>& frac){
    int index;
    frac.resize(mech.SpeciesCount(),0.0);
    vector<doublereal> temp;
    temp.resize(mech.SpeciesCount(),0.0);
    map<string,doublereal>::iterator p;
    p = spec.begin();
    while(p!=spec.end()){
        index = mech.FindSpecies(convertToCaps(trim(p->first)));
        if(index < 0)
           throw CamError("Species "+p->first +" not found in species list\n");
        else
            temp[index] = p->second;
        p++;
    }

    if(getFracType() == MASS){
        frac = temp;
    }else{
        CamConverter cc;
        cc.mole2mass(temp,frac,mech);
    }

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
