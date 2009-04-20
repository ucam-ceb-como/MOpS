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
#include "radau_wrapper.h"
#include "cam_profile.h"
using namespace Sprog;
namespace Camflow{
    class CamPremix : public CamSetup {

        
    public:

        typedef struct{
            vector<doublereal> Species;
            doublereal FlowRate;
            doublereal Vel;
            doublereal Dens;
            doublereal T;
            vector<doublereal> Dk;
            vector<doublereal> jk;
        } inletStruct;

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
        //write the results to output file
        void reportToFile(doublereal t, doublereal* soln);
        //mass matrix evaluation
        void massMatrix(doublereal **M);


        //solve
        void solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg,CamProfile &cp, Mechanism &mech );

        //return the initial solution vector
        void getInitial(vector<doublereal>& initial);

        //residual function definition----------------
        void residual(const doublereal& t, doublereal *y, doublereal *f);

        //boundary residual
        void boundaryCondition(const doublereal& t, doublereal *y, doublereal *f);

        //set up the solution vector
        void setupSolutionVector(CamBoundary &cb, CamControl &cc);

        //header information
        void header();

        //store the inlet
        void storeInlet(CamBoundary &cb);
        //update the diffusion fluxes
        void updateDiffusionFluxes();
        //update the thermal fluxes
        void updateThermo();

    protected:
        inletStruct ud_inlet;

        
    };
}


#endif	/* _CAM_PREMIX_H */

