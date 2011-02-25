/* 
 * File:   cam_models.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2008 Vinod M Janardhanan
 *
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
 *
 * Created on January 18, 2009, 6:50 PM
 */

#ifndef _CAM_MODELS_H
#define	_CAM_MODELS_H
#include "cam_admin.h"
#include "cam_boundary.h"
#include "cam_configuration.h"
#include "cam_control.h"
#include "cam_geometry.h"
#include "cam_profile.h"
#include "gpc.h"
#include "cam_soot.h"
#include "cam_residual.h"
#include "flamelet.h"
using namespace Sprog;
namespace Camflow{
    class CamModels{

    public:
        CamModels(){};
        ~CamModels(){};
        /*
         *this caling interface is for the kernel
         */
        void solve(CamAdmin &ca,                   
                   CamConfiguration &config,
                   CamControl &cc,
                   CamGeometry &cg,
                   CamProfile &cp,
                   CamSoot &cs,
                   Mechanism &mech
                   );

//        /*
//         *calling interface for external programs
//         */
//        void solve(CamAdmin &ca,
//                   CamConfiguration &config,
//                   CamControl &cc,
//                   CamGeometry &cg,
//                   CamProfile &cp,
//                   CamSoot &cs,
//                   Mechanism &mech,
//                   vector<Thermo::Mixture>& mixtures
//                   );
        /*
         *flamelet interface. This is provided since the handle to the
         *object call still remain with the interface object. That helps
         *querying the reactor object in a much easier fashion
         */
//        void solve(CamAdmin &ca,
//                   CamBoundary &cb,
//                   CamConfiguration &config,
//                   CamControl &cc,
//                   CamGeometry &cg,
//                   CamProfile &cp,
//                   CamSoot &cs,
//                   Mechanism &mech,
//                   FlameLet &flmlt,
//                   doublereal sdr
//                   );

    };
}


#endif	/* _CAM_MODELS_H */

