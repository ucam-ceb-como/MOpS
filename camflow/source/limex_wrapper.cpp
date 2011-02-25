#include <iostream>
#include "limex_wrapper.h"
#include "cam_residual.h"

using namespace std;
using namespace Camflow;
using namespace Gadgets;









extern "C" {

    typedef void (*ResidFcn) (int* n, int* nz, doublereal* y,
            doublereal* f, doublereal* b, int* ir, int* ic, int* FcnInfo );
    
    typedef void (*JacFcn)();
    doublereal tInitial;

    /**
     *  Limex integration call
     */
//    extern void xlimex_(int* n , ResidFcn fcn, JacFcn jacobian, doublereal* 
//            t_begin,
//            doublereal* t_end, doublereal* y, doublereal* ys,
//            doublereal* rtol, doublereal* atol, doublereal* h,
//            int* iopt, doublereal* ropt, int* ipos, int* ifail);
    /**
     *  Limex integration call
     */
  

   
            extern void XLIMEX(int* n, ResidFcn fcn, JacFcn jacobian,
            doublereal* t_begin, doublereal* t_end, doublereal* y, doublereal* ys,
            doublereal* rtol, doublereal* atol, doublereal* h,
            int* iopt, doublereal* ropt, int* ipos, int* ifail);
   

    /**
     *  Function called by LIMEX for the evaluation of residuals
     */
    static void limex_res(int* n, int* nz, doublereal* y,
            doublereal* f, doublereal* b, int* ir, int* ic, int* FcnInfo){


    //Cast one of xlimex's arguments to the base class, in order to call the appropriate reactor model.
    //Check to see if this is the right way to do it.


    void **iddres_res = reinterpret_cast<void **> (&(FcnInfo[0]));
    void *hndl = *iddres_res; 
    Camflow::CamResidual* g = (Camflow::CamResidual*)hndl;
    g->CamResidual::speciesResidual(tInitial,y,f);

   //Question is whether we need to create a residual evaluation function, or whether we can use the one 
   //in CamResidual?

    }

    static void limex_jac(){}

   

}



/**
 *  Default constructor sets all the default
 *  setting, which may be modified by the
 *  public methods
 */
LimexWrapper::LimexWrapper(){

    iOpt.resize(30,0);

    //resuse the Jacobian
    iOpt[9] = 1;
    //vector tolerances
    iOpt[10] = 1;
    //return after each step
    iOpt[11] = 1;

    rOpt.resize(5,0);

    //error indicator
    iFail.resize(3,0);
}

/**
 *  Initialize LIMEX
 *  @nEq    : number of equations
 *  @y        : vector of solution variables
 *  @tStart: time for start of integration
 *  @tFinal : time for end of integration
 */
void LimexWrapper::initLimex(const int nEq,std::vector<doublereal> y, const doublereal tStart)
     {
 
    n = nEq;
    eqnSize = nEq;

    if(n <= 0)
        throw runtime_error("Number of equations zero less\n");
    v_y = y;
    if(n != (int)v_y.size())
        throw runtime_error("size of solution vector differs from number");
    
    //set the derivatives to zero
    v_yPrime.resize(n,0);

    //integration times
    tBegin = tStart;
    tInitial = tStart;
}

/**
 *  Set the tolerance parameters
 *  @rTols : vector of relatIve tolerances
 *  @aTols  : vector of absolute tolerances
 */


void LimexWrapper::setTolerances(const std::vector<doublereal> rTols, const std::vector<doublereal> aTols){
    int lrtol, latol;
    lrtol = (int)rTols.size();
    latol = (int)aTols.size();
    if(lrtol != n || latol !=n )
        throw runtime_error("Size of tolerance vectors differ from N equations\n");

    v_aTol = aTols;
    v_rTol = rTols;
}

/**
 *  Set the band widths for the Jacobian matrix
 *  @upper  : number of upper diagonals
 *  @lower  : number of lower diagonals
 */
void LimexWrapper::setBandWidth(const int upper, const int lower){
    iOpt[7] = lower;
    iOpt[8] = upper;
}


/**
 *  Limex integration call
 */
void LimexWrapper::integrate(const doublereal tFinal){

    tEnd = tFinal;
    do {

         XLIMEX(&n, limex_res, limex_jac, &tBegin, &tEnd,
                 &v_y[0], &v_yPrime[0], &v_rTol[0], &v_aTol[0],
                 &iniStep, &iOpt[0], &rOpt[0], &iPos[0], &iFail[0]);

         if(iFail[0] < 0){
             iniStep /= 2.0;
         }
         
    }while(tBegin < tEnd);

}



/**
 *  Limex DAE call
 */
void LimexWrapper::integrateDAE(const doublereal tFinal){

     static int iter = 0;
     do{
        if(iter%10 == 0){
        //solve the algebraic equation system
            reactorPtr->calcFlowField(tBegin, &v_y[0]);
        }
       
     //call Limex to solve the ODE
   
     tEnd = tFinal;
     do {
         XLIMEX(&n, limex_res, limex_jac, &tBegin, &tEnd,
                 &v_y[0], &v_yPrime[0], &v_rTol[0], &v_aTol[0],
                 &iniStep, &iOpt[0], &rOpt[0], &iPos[0], &iFail[0]);

         if(iFail[0] < 0){
             iniStep /= 2.0;
         }
        iter++;
    }while(tBegin < tEnd);
  }while (resNorm > resTol);
}

void LimexWrapper::calcResidualNorm(){
    resNorm = 0;
    doublereal *yp;
    yp = &v_yPrime[0];
    for(int i=0; i<eqnSize; i++)
        resNorm += yp[i]*yp[i];

    resNorm = sqrt(resNorm);
}





//Stuff that was in the radau wrapper that we may need here:
// Note that because the types are different we will need to use the reinterpretor cast operator.



/*
    * Set the function that the solver should call for evaluation of the right hand side function.
    * This sould give a default of fcn evaluated to 1, and the user would need to define an alternative
    * if it is needed.  Key is to determine what the function would normally need to be in the case of each 
    * of the solvers.
    */

   // solver->setRHSFcn(&residual);



   /*
    * An object that can be passed to the solver.
    * On calling the residual function this will be returned.
    *
    * The idea here would be to enable the user to set more 
    * parameters, and to make it possible to return data to 
    * the user.
    */



   /* solver->setUserData((void*)&cr);
    * This function could make it possible to select from other options in 
    * camxml.  For example, dense matrix solver, or banded matrix.  It would
    * also take data back.  
   */ 




   /*
    *Report function
    */

  // solver->setReportFcn(&report);
  // Here we would also need to change the typee of the data in order to get it back into
 //  Camflow.  


