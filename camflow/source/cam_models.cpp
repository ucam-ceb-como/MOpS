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
void CamModels::solve(CamAdmin& ca, CamBoundary& cb,
        CamConfiguration& config,
        CamControl& cc,
        CamGeometry& cg,
        CamProfile& cp,
        CamSoot& cs,
        Mechanism& mech){

    int configID;
    configID = config.getConfiguration();

    

    if(configID == config.PLUG){
        CamPlug cplug;
        try{
            cplug.solve(cc,ca,cg,cp,mech);
        }catch(CamError ce){
            cout << ce.errorMessge << endl;
        }
        
        
    }else if(configID == config.PREMIX){
        CamPremix cpremix;
        try{
            cpremix.solve(cc,ca,cg,cp,cs,mech);
        }catch(CamError ce){
            cout << ce.errorMessge << endl;
        }
        
    }else if(configID == config.STAGFLOW || configID == config.COUNTERFLOW){
        
        StagFlow stflow;

        try{
            stflow.solve(cc,ca,cg,cp,config, mech);
        }catch(CamError ce){
            cout << ce.errorMessge << endl;
        }
        
    }else if(configID == config.BATCH_CV){
        Batch batch;
        try{
            batch.solve(cc,ca,cg,cp,cs,mech);
        }catch(CamError ce){
            cout << ce.errorMessge << endl;
        }
        
    }else if(configID == config.FLAMELET){
        FlameLet flmlt;
        
        try{
            flmlt.solve(cc,ca,cg,cp,mech,false);
        }catch(CamError ce){
            cout << ce.errorMessge << endl;
        }
        
    } else{
        cout << "Unknown reactor model\n";
    }
    
}


//void CamModels::solve(CamAdmin& ca,
//                        CamBoundary& cb,
//                        CamConfiguration& config,
//                        CamControl& cc,
//                        CamGeometry& cg,
//                        CamProfile& cp,
//                        CamSoot& cs,
//                        Mechanism& mech,
//                        FlameLet& flmlt, doublereal sdr){
//
//        flmlt.setExternalScalarDissipationRate(sdr);
//        try{
//            flmlt.solve(cc,ca,cg,cp,mech,true);
//        }catch(CamError ce){
//            cout << ce.errorMessge << endl;
//        }
//}

