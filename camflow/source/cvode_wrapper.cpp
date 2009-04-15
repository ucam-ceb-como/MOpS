
#include "cvode_wrapper.h"
using namespace Camflow;

extern "C"{

    int cvodeResid(doublereal time, N_Vector y, N_Vector ydot, void *udata){
        CamResidual *residual = (CamResidual*)(udata);
        residual->eval(time,NV_DATA_S(y),NV_DATA_S(ydot), true);
        return 0;
    }
}

void CVodeWrapper::init(int n, vector<doublereal>& solnVec, doublereal tol,
               doublereal rtol,doublereal maxIntTime ,int band, CamResidual& cr){

    reacPtr = &cr;

    cvode_mem = NULL;
    y=NULL;

    /*
     *create the soln vector
     */
    y = N_VMake_Serial(n,&solnVec[0]);
    atol = tol;
    currentTime = 0.0;
    maxTime = maxIntTime;
    cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);
    
    CVodeMalloc(cvode_mem,cvodeResid,currentTime,y,CV_SS,rtol,(void*)&atol);


    if(band==n){
        CVDense(cvode_mem,n);
    }else{
        CVBand(cvode_mem,n,band,band);
    }
    CVodeSetFdata(cvode_mem,(void*)&cr);
    CVodeSetMaxNumSteps(cvode_mem,5000);

}

/*
 *additional solver control
 */
void CVodeWrapper::setIniStep(doublereal istep){
    CVodeSetInitStep(cvode_mem,istep);
}

void CVodeWrapper::solve(int stopMode){

    int flag;
    //CVodeSetStopTime(cvode_mem,maxTime);
    do{
        flag = CVode(cvode_mem,maxTime,y,&currentTime,stopMode);
        if(flag < 0){
            cout << "Cvode Integration error\n";            
        }else{
            //cout << currentTime << endl;
            reacPtr->report(currentTime,NV_DATA_S(y));
        }
    }while(currentTime < maxTime);

}

void CVodeWrapper::destroy(){
    CVodeFree(&cvode_mem);
    N_VDestroy(y);
}