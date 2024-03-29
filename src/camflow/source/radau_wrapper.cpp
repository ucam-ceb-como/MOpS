/*
 * File:   RadauWrapper.cpp
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2009 Vinod M Janardhanan.
 *
 * File purpose:
 *  This is a wrapper class for Radau C++ solver
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
 * Created on January 24, 2009, 6:03 PM
 */


#include "radau_wrapper.h"
#include "cam_error.h"
#include "cam_residual.h"
using namespace Camflow;

//RadauWrapper::~RadauWrapper(){
//
//}
void RadauWrapper::setControl(CamControl &cc){
    solverControl = &cc;
}

void RadauWrapper::setBandWidth(int n){
    bandWidth = n;
}

extern "C"{
    /*
     *function called by RADAU for evaluation of the residual function
     *@x - independent variable
     *@y - solution vector at x
     *@f - residual function f dy/dt or dy/dz
     */
    void residual(double x, double* y, double* f, void* udata, bool jacCall){
        CamResidual *residual = (CamResidual*)udata;
        residual->eval(x,y,f, jacCall);
    }

    void Jacobian(double x, double *y, double **J)
    {

    } // Jacobian

    /*
     *function called by radau for writing the output
     */
    double report(double x, double* y, void* udata){

        CamResidual *residual = (CamResidual*)udata;
        double resNorm = residual->getResidual();
        // Output time and residual to screen.
        std::cout << x << "    " << resNorm << std::endl;

        // If the residual is below a certain value, stop RADAU.
        //if(resNorm < 1e-5){return 1;}

	    return resNorm;

    }

    void massMatrix(double **M, void* udata){
        CamResidual *residual = (CamResidual*)udata;
        residual->massMatrix(M);
    }

}


//solver init
void RadauWrapper::initSolver(  int nEq,
                                const double tBeg,
                                const double tEnd,
                                std::vector<double>& y,
                                std::vector<double>& rTol,
                                std::vector<double>& aTol,
                                CamResidual &cr){
    /*
     *initialisation of the Radau solver
     *all fine tuning for the solver performance should be
     *done here
     */

   // initial value for x
   //double tBeg = 0.0;
   // final value for x
   if(tEnd <= 0) throw CamError("Invalid end time for integration\n");

   // interval of x for printing output
   double dx = 1e-3;
   // rtoler and atoler are vectors
   int itoler = 0;
   // use SolutionOutput routine
   const int iout = 1;
   // numerical Jacobian
   const int ijac = 0;

   // Mass matrix routine is identity
   const int imas = 0;
   int mlmas = nEq;
   int mumas = 0;

   // Use default values (see header files) for these parameters:
   double hinit = 0.0;//solverControl->getIniStep();   //initial step size
   double hmax(0.0);                             //max step size
   int nmax(100000);                                  //max number of allowed steps (100000)
   int nit(1000);                                   //max number of newton iteration in each step (7)
   double uround(1e-16), safe(0.0);                //safety factor for step size prediction
   double facl(0.0), facr(0.0);

   bool startn(true);
   int nind1(nEq), nind2(0), nind3(0), npred(1), m1(0), m2(0);
   bool hess(true);
   double fnewt(0.1), quot1(1), quot2(1.2);
   double thet(0.001);                             // Jac evaluation decision (0.001) increase of eval is expensive

   solver = new StiffIntegratorT(
                nEq,                    //system size
                &y[0],                  //solution vector
                tBeg,
                tEnd,
                dx,
                itoler,
                &rTol[0],
                &aTol[0],
                iout,
                hinit,
                hmax,
                nmax,
                uround,
                safe,
                facl,
                facr,
                ijac,
                bandWidth, bandWidth,
                imas,
                mlmas, mumas,
                nit, startn,
                nind1, nind2, nind3,
                npred,
                m1, m2,
                hess, fnewt,
                quot1, quot2,
                thet);

   /*
    * set the function that the solver should call for
    * evaluation of the right hand side function
    */

   solver->setRHSFcn(&residual);

   /*
    * An object that can be passed to the solver.
    * On calling the resiual function this will be
    * returned
    */

   solver->setUserData( (void*)&cr);

   /*
    *Report function
    */
   solver->setReportFcn(&report);

   /*
    *set the mass matrix function
    */
   solver->setMass(&massMatrix);
   solver->setJacobian(&Jacobian);
}


void RadauWrapper::Integrate(){
    /*
     *Actual integration call
     */
    solver->Integrate();
}
