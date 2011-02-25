
#include <stdlib.h>
#include <stdio.h>
#include "giant_wrapper.h"
using namespace Camflow;

extern "C"{
    #define AJDEL 1.0e-08
    #define AJMIN 1.0e-04
    void giantResidual(int *neq, doublereal* y, doublereal* f,void *uData ,
                                            int *fail){

        CamResidual *residual = (CamResidual*)uData;
        residual->eval(0,y,f,false);
        fail = 0;
    }

    void giantJac(int n ,doublereal* x, doublereal* fscale ,doublereal* fxk, int njac,
                            void *uData, int *fail){
        vector<doublereal> f, fdelta;
        f.resize(n);fdelta.resize(n);
        CamResidual *residual = (CamResidual*)uData;
        residual->eval(0,x,&f[0],false);
        

        int i,k,kn;
        for(k=0;k<n;k++){
            doublereal w=x[k];
            doublereal u = fabs(x[k]);
            u = MAX(u,AJMIN)*AJDEL*SIGN(x[k]);
            x[k] = w+u;
            
            residual->eval(0,x,&fdelta[0],false);
            
            x[k] = w;
            kn = k*n;
            for(i=0; i<n;i++){
                fxk[kn+i] = (fdelta[i]-f[i])/u;
            }
        }
        njac++;
        fail = 0;
    }
    
}


void GiantWrapper::solve(int neq,
        doublereal tol,vector<doublereal>& iniGuess,CamResidual& cr){
    /*
     *memory allocation for the structures
     */
//      resid = (GIANT_FUN*)malloc(sizeof(struct GIANT_FUN));
//      optional = (GIANT_OPT*)malloc(sizeof(struct GIANT_OPT));
//      info = (GIANT_INFO*)malloc(sizeof(struct GIANT_INFO));

      resid = (NLEQ_FUN*)malloc(sizeof(struct NLEQ_FUN));
      optional = (NLEQ_OPT*)malloc(sizeof(struct NLEQ_OPT));
      info = (NLEQ_INFO*)malloc(sizeof(struct NLEQ_INFO));

      /*
       *set the residual function pointer
       */
      resid->fun = &giantResidual;
      /*
       *set the jacobian function
       */
      //resid->jac = &giantJac;
      resid->jac = NULL;
      /*
       *optional parameter setting for giant
       */
//      optional->tol = tol;
//      optional->maxiter = 500;
//      optional->lin_maxiter = 5000;
//      optional->i_max = 10;
//      optional->datalevel = None;
//      optional->iterfile = NULL;
//      optional->resfile = NULL;
//      optional->miscfile = NULL;
//      optional->datafile = NULL;
//      optional->nonlin = Highly_Nonlinear;
//      optional->restricted = False;
//      optional->scale = NULL;
//      optional->scaleopt = GIANT_OPT::StandardScale;
//      optional->monitorlevel = None;
//      optional->errorlevel = None;
//      optional->monitorfile = NULL;
//      optional->linmonlevel = None;

      optional->uData = &cr;

      optional->datalevel = None;
      optional->errorlevel = None;
      optional->monitorlevel = None;
      optional->datafile = NULL;
      optional->errorfile = NULL;
      optional->iterfile = NULL;
      optional->miscfile = NULL;
      optional->monitorfile = NULL;
      optional->resfile = NULL;      
      optional->scaleopt = NLEQ_OPT::StandardScale;

      optional->tol = tol;
      optional->maxiter = 500;
      optional->nonlin = Extremely_Nonlinear;
      optional->restricted = False;
      optional->nleqcalled = False;





      /*
       *solver call
       */
      //giant_gmres(*resid,neq,&iniGuess[0],optional,info);
      nleq_res(*resid,neq,&iniGuess[0],optional,info);



}