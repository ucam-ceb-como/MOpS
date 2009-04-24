
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
    yPrime = NULL;
    /*
     *create the soln vector
     */
    y = N_VMake_Serial(n,&solnVec[0]);
    yPrime = N_VNew_Serial(n);
    eqnSize = n;
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

doublereal& CVodeWrapper::solve(int stopMode){

    int flag;    
    do{
        flag = CVode(cvode_mem,maxTime,y,&currentTime,stopMode);
        if(flag < 0){
            cout << "Cvode Integration error\n";            
        }else{
            reacPtr->report(currentTime,NV_DATA_S(y));
        }
    }while(currentTime < maxTime);
    
    CVodeGetDky(cvode_mem,currentTime,1,yPrime);
    calcResNorm();
    return resNorm;

}

void CVodeWrapper::solve(int stopMode, doublereal resTol){

    int flag;
    do{
        flag = CVode(cvode_mem,maxTime,y,&currentTime,stopMode);
        if(flag < 0){
            cout << "Cvode Integration error\n";
        }else{
            CVodeGetDky(cvode_mem,currentTime,1,yPrime);
            calcResNorm();
            reacPtr->report(currentTime,NV_DATA_S(y),resNorm);
        }
    }while(resNorm > resTol);

}

void CVodeWrapper::destroy(){
    CVodeFree(&cvode_mem);
    N_VDestroy(y);
    N_VDestroy(yPrime);
}

void CVodeWrapper::calcResNorm(){
    resNorm = 0;
    doublereal *yp;
    yp = NV_DATA_S(yPrime);
    for(int i=0; i<eqnSize; i++)
        resNorm += yp[i]*yp[i];

    resNorm = sqrt(resNorm);
}
