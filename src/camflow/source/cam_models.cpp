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

using namespace Camflow;
using namespace std;

CamModels::CamModels()
:
  rModel_(NULL)
{}

CamModels::~CamModels()
{
    if (rModel_ != NULL) delete rModel_;
}

/*
 *this member function chooses the reactor model based on the
 *configuration ID and calls the passes the control to the right
 * reactor object
 */
void CamModels::solve
(
    CamAdmin& ca,
    CamConfiguration& config,
    CamControl& cc,
    CamGeometry& cg,
    CamProfile& cp,
    CamSoot& cs,
    Mechanism& mech
)
{

    int configID;
    configID = config.getConfiguration();

    if(configID == config.PLUG){
        //rModel = new CamPlug();
    }else if(configID == config.PREMIX){
        //rModel = new CamPremix();
    }else if(configID == config.BATCH_CV){
        rModel_ = new Batch(ca, config, cc, cg, cp, cs, mech);
    }else if(configID == config.STAGFLOW || configID == config.COUNTERFLOW){
        rModel_ = new StagFlow(ca, config, cc, cg, cp, cs, mech);
    }else if(configID==config.FLAMELET){
        rModel_ = new FlameLet(ca, config, cc, cg, cp, cs, mech);
    }else{
        throw CamError("Unknown reactor model\n");
    }
    try{
        rModel_->solve();
    }catch(CamError &ce){
        cout << ce.errorMessage << endl;
    }

}

