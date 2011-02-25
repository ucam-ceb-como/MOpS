#include "ida_wrapper.h"
#include "cam_read.h"
using namespace Camflow;


extern "C"{

    int idaRes(doublereal tt, N_Vector yy, N_Vector yp, N_Vector rr,
            void *udata){
        CamResidual *residual = (CamResidual*)(udata);
        int n = residual->getNEqn();        
        residual->eval(tt,NV_DATA_S(yy),NV_DATA_S(rr),true);
        doublereal *rsPtr = NV_DATA_S(rr);
        doublereal *ypPtr = NV_DATA_S(yp);
        for(int i=0; i<n; i++){
            rsPtr[i] -= ypPtr[i];
        }
        return 0;

    }
}




void IDAWrapper::init(int n,
                    vector<doublereal>& solnVec,
                    doublereal tol,
                    doublereal maxIntTime,
                    int band,
                    CamResidual& cr){
    reacPtr = &cr;
    ida_mem = NULL;
    y0 = NULL;
    yp0 = NULL;

    /*
     *create solution vector
     */
    y0 = N_VMake_Serial(n,&solnVec[0]);
    yp0 = N_VNew_Serial(n);
    N_VConst_Serial(0.0,yp0);

    /*
     *starting time
     */
    currentTime = 0.0;
    atol = tol;
    maxTime = maxIntTime;

    ida_mem = IDACreate();

    flag = IDAMalloc(ida_mem,idaRes,currentTime,y0,yp0,IDA_SS,tol,(void*)&atol);
    if(flag != IDA_SUCCESS) throw CamError("IDA memory alloc failure\n");

    /*
     *optional setting
     */
    flag = IDASetInitStep(ida_mem,1e-06);
    flag = IDASetRdata(ida_mem,(void*)&cr);
    if( n!= band)
        flag = IDABand(ida_mem,n,band,band);
    else
        flag = IDADense(ida_mem,n);
    
    flag = IDASetMaxNonlinIters(ida_mem,20);
    flag = IDASetMaxNumSteps(ida_mem,5000);

}

void IDAWrapper::solve(){
    do{
        flag = IDASolve(ida_mem,maxTime,&currentTime,y0,yp0,IDA_ONE_STEP);
        if(flag <0 ){
            cout << "Ida integration error\n";
        }else{
            cout << currentTime << endl;
        }
    }while(currentTime<maxTime);
}
