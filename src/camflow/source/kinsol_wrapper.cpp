
#include "kinsol_wrapper.h"

using namespace Camflow;

extern "C"{

    int func(N_Vector y, N_Vector fval, void *udata){
        CamResidual *residual = (CamResidual*)udata;
        residual->eval(NV_DATA_S(y),NV_DATA_S(fval));
		  return 0;
    }

}


void KinsolWrapper::init(int n, vector<doublereal>& solnVec, doublereal rtol, int band, CamResidual &cr){
    
  y = scale = NULL;
  kmem = NULL;

  y = N_VMake_Serial(n,&solnVec[0]);
  //if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  scale = N_VNew_Serial(n);
  //if (check_flag((void *)scale, "N_VNew_Serial", 0)) return(1);

    /* -----------------------------------------
   * Initialize and allocate memory for KINSOL
   * ----------------------------------------- */

  kmem = KINCreate();
  //if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

  /* y is used as a template */

  flag = KINMalloc(kmem, func, y);
  //if (check_flag(&flag, "KINMalloc", 1)) return(1);

  /* -------------------
   * Set optional inputs
   * ------------------- */

  /* Specify stopping tolerance based on residual */

  fnormtol  = rtol;
  flag = KINSetFuncNormTol(kmem, fnormtol);
  //if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);

  /* -------------------------
   * Attach band linear solver
   * ------------------------- */
  if(band < n)
      flag = KINBand(kmem,n,band,band);
  else
      flag = KINDense(kmem,n);

  flag = KINSetNumMaxIters(kmem,5000);

  flag = KINSetFdata(kmem,(void*)&cr);


  /* No scaling used */
  N_VConst_Serial(1.0,scale);

}

void KinsolWrapper::solve(){
    flag = KINSol(kmem,y,KIN_LINESEARCH,scale,scale);
    
    if(flag != KIN_SUCCESS){
        cout << "KINsol Error " << flag << endl;
        string errorStr = KINGetReturnFlagName(flag);
        errorStr += "\n";
        throw CamError(errorStr);
    }
    
    
}

void KinsolWrapper::destroy(){
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(scale);
    KINFree(&kmem);
}
