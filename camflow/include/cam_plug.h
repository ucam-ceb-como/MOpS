/* 
 * File:   cam_plug.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class implements the plug flow model
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
 * Created on January 24, 2009, 12:08 PM
 */

#ifndef _CAM_PLUG_H
#define	_CAM_PLUG_H
#include "cam_residual.h"
#include "cam_control.h"
#include "cam_admin.h"
#include "cam_configuration.h"
#include "cam_reporter.h"
#include "gpc.h"
#include "cam_setup.h"

using namespace Sprog;
namespace Camflow{
    class CamPlug: public CamSetup{

    public:
        CamPlug(){};
        virtual ~CamPlug(){}
        
        /*
         *the following 3 functions are called by
         *the DAE solver
         */
        int eval(doublereal x, doublereal* y, doublereal* ydot, bool jacEval);
        //console output
        void report(doublereal x, doublereal* soln);
        //console output with residuals
        void report(doublereal x, doublereal* soln, doublereal& res);
        //prepare the data vector for output
        void vectorize(doublereal x, doublereal* soln, vector<doublereal>& data);
        //create the cummary file
        void createSummary();
        //write the summary
        void reportSummary(doublereal x, doublereal * soln);
        //mass matrix evaluation
        void massMatrix(doublereal **M);

        
        //solve
        void solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg,CamProfile &cp,
             CamConfiguration &config, CamSoot &cs,  Mechanism &mech );

        //return the initial solution vector
        void getInitial(vector<doublereal>& initial);
        
        //residual function definition----------------
        void residual(const doublereal& x, doublereal *y, doublereal *f);
        //species residual
        void speciesResidual(const doublereal& x, doublereal *y, doublereal *f);
        //mass flow
        void massFlowResidual(const doublereal& x, doublereal *y, doublereal *f);
        //temperature
        void energyResidual(const doublereal& x, doublereal *y, doublereal *f);
        //residence time
        void residenceTime(const doublereal& x, doublereal *y, doublereal *f);
        
        //update the mixture properties
        void updateMixture(const doublereal& x, doublereal *y);

        //header information
        void header();

    };
}
#endif	/* _CAM_PLUG_H */

