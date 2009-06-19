/* 
 * File:   cam_kernel.cpp
 * Author: vinod
 *
 * Created on January 17, 2009, 1:47 PM
 * Copyright (C) 2008 Vinod M Janardhanan.

 * File purpose:
 *  Main file controling the program execution
 * Licence:
 *  This file is part of "Camflow".
 *
 *  flameLab is free software; you can redistribute it and/or
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

#include <iostream>
#include <string>
#include "gpc.h"
#include "cam_read.h"
#include "cam_admin.h"
#include "cam_control.h"
#include "cam_geometry.h"
#include "cam_converter.h"
#include "cam_boundary.h"
#include "cam_profile.h"
#include "cam_configuration.h"
#include "cam_read.h"
#include "cam_models.h"
#include "cam_soot.h"
/*
 * 
 */
using namespace std;
using namespace Sprog;
using namespace Camflow;

int main() {

    string fChem("chem.inp");
    string fThermo("therm.dat");
    string fTrans("tran.dat");
    string fCamFlow("camflow.xml");

    //mechanism object
    static Mechanism mech;

    //read Camflow input file
    CamControl cc;
    CamGeometry cg;
    CamConverter convert;
    CamAdmin ca;
    CamBoundary cb;
    CamProfile cp;
    CamConfiguration config;
    CamRead cm;
    CamModels models;
    CamSoot cSoot;
    try{
        cm.readInput(fCamFlow,cc,cg,convert,ca,cb,cp,config,cSoot);
    }catch(CamError &ce){
        cout << ce.errorMessge;
    }
    //CamBoundary *cb2 = &cb;
    
    //read mechanism, thermo and trasnport data
    IO::MechanismParser::ReadChemkin(fChem,mech,fThermo,fTrans);
    
    try{
        models.solve(ca,cb,config,cc,cg,cp,cSoot,mech);
    }catch(CamError &ce){
        cout << ce.errorMessge;
    }
    cout << "\nCamflow end of execution..\n";
    return (EXIT_SUCCESS);
}

