/* 
 * File:   cam_premix.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class implements the premix flame model
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
 * Created on January 18, 2009, 7:57 PM
 */

#ifndef _CAM_PREMIX_H
#define	_CAM_PREMIX_H
#include "gpc.h"
#include "cam_residual.h"
#include "cam_control.h"
#include "cam_control.h"
#include "cam_boundary.h"
#include "cam_profile.h"
#include "cam_setup.h"
#include "cam_profile.h"
#include "cam_soot.h"
using namespace Sprog;
namespace Camflow{
    class CamPremix : public CamSetup {

        
    public:

        CamPremix(){};
        //CamPremix(){}
        virtual ~CamPremix(){};


        /*
         *the following 3 functions are called by
         *the DAE solver
         */
        int eval(doublereal x, doublereal* y, doublereal* ydot,bool jacEval);
        //report the results
        void report(doublereal t, doublereal* solution);
        //console output with resisual monitoring
        void report(doublereal t, doublereal* solution, doublereal& res);
        //write the results to output file
        void reportToFile(doublereal t, doublereal* soln);
        //mass matrix evaluation
        //void massMatrix(doublereal **M);


        //solve
        void solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg,CamProfile &cp,
                CamSoot &cs,Mechanism &mech );
        void csolve(CamControl &cc);
        void ssolve(CamControl &cc);

        //return the initial solution vector
        void getInitial(vector<doublereal>& initial);

        //residual function definition----------------
        void residual(const doublereal& t, doublereal *y, doublereal *f);

        //boundary residual
        void massFlowBoundary(const doublereal& t, doublereal *y, doublereal *f);
        void speciesBoundary(const doublereal& t, doublereal *y, doublereal *f);
        void energyBoundary(const doublereal& t, doublereal *y, doublereal *f);
        void momentBoundary(const doublereal& t, doublereal *y, doublereal *f);

        //set up the solution vector
        void initSolutionVector(CamBoundary &cb, CamControl &cc);

        //header information
        void header();

        //store the inlet
        //void storeInlet(CamBoundary &cb);
        //update the diffusion fluxes
        void updateDiffusionFluxes();
        //update the thermal fluxes
        void updateThermo();
        //save the flow variables
        void saveFlowVariables(doublereal* y);

    protected:
        inletStruct ud_inlet;
        doublereal resNorm;
        

        
    };
}


#endif	/* _CAM_PREMIX_H */

