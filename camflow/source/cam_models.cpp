/*
 * File:   cam_models.cpp
 * Author: vinod
 *
 * Created on January 18, 2009, 6:50 PM
 * File purpose:
 *  This class decides which model to simulate based on the input
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
 */
#include "cam_models.h"
#include "cam_plug.h"
#include "batch.h"
#include "cam_premix.h"
#include "stagflow.h"
#include "flamelet.h"
#include "cam_error.h"
#include <iostream>
using namespace Camflow;
using namespace std;
/*
 *this member function chooses the reactor model based on the
 *configuration ID and calls the passes the control to the right
 * reactor object
 */
void CamModels::solve(CamAdmin& ca, 
        CamConfiguration& config,
        CamControl& cc,
        CamGeometry& cg,
        CamProfile& cp,
        CamSoot& cs,
        Mechanism& mech){

    int configID;
    configID = config.getConfiguration();


    CamResidual *rModel;
    if(configID == config.PLUG){
        rModel = new CamPlug();
    }else if(configID == config.PREMIX){
        rModel = new CamPremix();
    }else if(configID == config.BATCH_CV){
        rModel = new Batch();
    }else if(configID == config.STAGFLOW || configID == config.COUNTERFLOW){
        rModel = new StagFlow();
    }else if(configID==config.FLAMELET){
        rModel = new FlameLet();
    }else{
        throw CamError("Unknown reactor model\n");
    }
    try{
        rModel->solve(cc,ca,cg,cp,config,cs,mech);
    }catch(CamError &ce){
        cout << ce.errorMessge << endl;
    }

}
/*
 *this function is involked if the call initiates from the
 *external interface. Once the steady state is attained
 *species properties are set to the mixtures
 */

//void CamModels::solve(CamAdmin& ca,
//                                    CamConfiguration& config,
//                                    CamControl& cc,
//                                    CamGeometry& cg,
//                                    CamProfile& cp,
//                                    CamSoot& cs,
//                                    Mechanism& mech,
//                                    vector<Thermo::Mixture>& mixtures){
//
//    int configID;
//    configID = config.getConfiguration();
//
//
//    CamResidual *rModel;
//    if(configID == config.PLUG){
//        rModel = new CamPlug();
//    }else if(configID == config.PREMIX){
//        rModel = new CamPremix();
//    }else if(configID == config.BATCH_CV){
//        rModel = new Batch();
//    }else if(configID == config.STAGFLOW || configID == config.COUNTERFLOW){
//        rModel = new StagFlow();
//    }else if(configID==config.FLAMELET){
//        rModel = new FlameLet();
//    }else{
//        throw CamError("Unknown reactor model\n");
//    }
//    try{
//        rModel->solve(cc,ca,cg,cp,config,cs,mech);
//    }catch(CamError &ce){
//        cout << ce.errorMessge << endl;
//    }
//
//    //Get the species mass fractions
//    Array2D massFracs;
//    rModel->getSpeciesMassFracs(massFracs);
//
//    //storage for density, velocity, and temperature
//    vector<doublereal> density, vel, temp;
//
//    //Get the density
//    rModel->getDensityVector(density);
//    //Get the velocity
//    rModel->getVelocity(vel);
//    //Get the temperature
//    rModel->getTemperatureVector(temp);
//
//    /*
//     * Set the properties to the mixtures. If the mixture
//     *size is zero, then create the mixture and push back
//     *to the vector
//     */
//    int nCells = cg.getnCells();
//    if(mixtures.size() >0){
//        if(mixtures.size() != nCells-2)
//            throw("size of mixtures is not consistant with the grid\n");
//        int nSp = mech.SpeciesCount();
//        for(int i=0; i<nCells-2;i++){
//            vector<doublereal> mf;
//            for(int l=0; l<nSp; l++){
//                mf.push_back(massFracs(i+1,l));
//            }
//            mixtures[i].SetMassFracs(mf);
//            mixtures[i].SetMassDensity(density[i+1]);
//            mixtures[i].SetTemperature(temp[i+1]);
//            mixtures[i].SetVelocity(vel[i+1]);
//        }
//    }else{
//
//        Thermo::Mixture mix(mech.Species());
//        int nSp = mech.SpeciesCount();
//        for(int i=0; i<nCells-2;i++){
//            vector<doublereal> mf;
//            for(int l=0; l<nSp; l++){
//                mf.push_back(massFracs(i+1,l));
//            }
//            mix.SetMassFracs(mf);
//            mix.SetMassDensity(density[i+1]);
//            mix.SetTemperature(temp[i+1]);
//            mix.SetVelocity(vel[i+1]);
//            mixtures.push_back(mix);
//        }
//    }
//}
//










