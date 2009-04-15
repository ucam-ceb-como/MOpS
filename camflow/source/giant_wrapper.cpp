
#include <cmath>
#include <stdlib.h>


#include "cam_residual.h"
#include "giant_wrapper.h"
#include "cam_read.h"

using namespace Camflow;
extern "C"{
    void giantFun(int n, doublereal *u, doublereal *fu, int nfcn, int *fail,
            void* udata){
        CamResidual *residual = (CamResidual*)(udata);
        residual->eval(0.0,u,fu,true);
        *fail = 0;
        
    }

    void giantJac(int n, doublereal *x, doublereal *xscale, doublereal *fx,
            int nJac, int *fail, void* udata){
        CamResidual *residual = (CamResidual*)(udata);

        /*
         *this is a dense approximation. On call fx
         *contains the current function values
         */
        *fail = 0;
        int kn;
        doublereal w,u;
        vector<doublereal> v, temp;
        v.clear();
        v.resize(n,0.0);
        for(int i=0; i<n; i++){            
            temp.push_back(fx[i]);
        }
        for(int k=0; k<n; k++){
            w = x[k];
            if(xscale)
                u = MAX(fabs(x[k]),xscale[k]);
            else
                u = fabs(x[k]);
            u = MAX(u,1e-04)*1e-08*SIGN(x[k]);
            cout << u << endl;
            x[k] = w+u;
            residual->eval(0.0,x,&v[0],true);
            x[k] = w;

//            kn = k*n;
//            for(int i=0; i<n; i++){
//                fx[kn+i] = (v[i]-temp[i])/u;
//            }

                        
        }

    }

 
}

void GiantWrapper::init(int n, vector<doublereal>& solnVec, doublereal rtol,
                        CamResidual& cr){



    /*
     *prepeate giant fun
     */
    fun.fun = &giantFun;
    fun.jac = &giantJac;

    if(n!=solnVec.size())
        throw CamError("Number of Eqns and soln vec size do not match\n");
    nEq = n;
    y = &solnVec[0];

    /*
     *prepare opt
     */

    opt.tol = rtol;
    opt.maxiter = 1000;
    opt.nonlin = Highly_Nonlinear;
    opt.restricted = False;
    opt.errorlevel = Verbose;
    opt.monitorlevel = Verbose;
    opt.datalevel = None;
    opt.errorfile = NULL;
    opt.monitorfile = NULL;
    opt.datafile = NULL;
    opt.iterfile = NULL;
    opt.resfile =NULL;
    opt.miscfile = NULL;
    opt.scale = NULL;
    opt.scaleopt = opt.StandardScale;
    opt.udata = (void*)&cr;
    

    /*
     *prepare info
     */
    
}

void GiantWrapper::solve(){

    giant_gmres(fun,nEq,y,&opt,&info);
}