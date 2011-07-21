
#include "cvode_wrapper.h"
using namespace Camflow;

extern "C"{

    int cvodeResid(doublereal time, N_Vector y, N_Vector ydot, void *udata){
        ((CamResidual*)(udata))->eval(time,NV_DATA_S(y),NV_DATA_S(ydot), false);
        return 0;
    }
}

CVodeWrapper::CVodeWrapper()
:
cvode_mem(NULL),
reacPtr(NULL)
{}

CVodeWrapper::~CVodeWrapper()
{
    destroy();
}

void CVodeWrapper::init(int n, std::vector<doublereal>& solnVec, doublereal tol,
               doublereal rtol,doublereal maxIntTime ,int band, CamResidual& cr, doublereal iniTime){

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
    currentTime = iniTime;
    maxTime = maxIntTime;
    cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON); // ank25: original
    //cvode_mem = CVodeCreate(CV_BDF, CV_FUNCTIONAL); // ank25: trying something different
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
/*
 *set the max allowed step size
 */
void CVodeWrapper::setMaxStep(doublereal maxStep){
    CVodeSetMaxStep(cvode_mem,maxStep);
}

doublereal& CVodeWrapper::solve(int stopMode){

    int flag;
    do{
        flag = CVode(cvode_mem,maxTime,y,&currentTime,stopMode);
        if(flag < 0){
            std::cout << "Cvode Integration error\n";
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
        //std::cout << "Calling CVODE" << std::endl;
        flag = CVode(cvode_mem,maxTime,y,&currentTime,stopMode);
        if(flag < 0){
            std::cout << "Cvode Integration error\n";
        }else{
            CVodeGetDky(cvode_mem,currentTime,1,yPrime);
            calcResNorm();
            //std::cout << "CvodeWrapper: ResNorm: " << resNorm << std::endl;
            reacPtr->report(currentTime,NV_DATA_S(y),resNorm);
        }
        if(currentTime > maxTime)
            break;
    }while(resNorm > resTol );

}

void CVodeWrapper::solveDAE(int stopMode, doublereal resTol){

    int flag;
    static int iter = 0;
    do{
        if(iter%10 == 0){
        //solve the algebraic equation system
            reacPtr->calcFlowField(currentTime,NV_DATA_S(y));
        }
        //call Cvode to solve the ODE
        flag = CVode(cvode_mem,maxTime,y,&currentTime,stopMode);
        if(flag < 0){
            std::cout << "Cvode Integration error\n";
        }else{
            CVodeGetDky(cvode_mem,currentTime,1,yPrime);
            calcResNorm();
            reacPtr->report(currentTime,NV_DATA_S(y),resNorm);
        }
        iter++;
    }while(resNorm > resTol);
}

void CVodeWrapper::destroy(){
    CVodeFree(&cvode_mem);
    //N_VDestroy(y);
    //N_VDestroy(yPrime);
}

void CVodeWrapper::calcResNorm(){
    resNorm = 0;
    doublereal *yp;
    yp = NV_DATA_S(yPrime);
    for(int i=0; i<eqnSize; i++)
        resNorm += yp[i]*yp[i];

    resNorm = sqrt(resNorm);
}
