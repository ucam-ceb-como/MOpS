
#include "cvode_wrapper.h"
using namespace Camflow;


extern "C"{

    int cvodeResid(double time, N_Vector y, N_Vector ydot, void *udata){
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

void CVodeWrapper::init(int n, std::vector<double>& solnVec, double tol,
               double rtol,double maxIntTime ,int band, CamResidual& cr, double iniTime){

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

void CVodeWrapper::initVectorTol(int n, std::vector<double>& solnVec, double aTolTemp[],
               double rtol,double maxIntTime ,int band, CamResidual& cr, double iniTime){

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
    ///atol = tol;					////  ank25
    currentTime = iniTime;
    maxTime = maxIntTime;

    std::cout << "Calling CVodeCreate " << std::endl;   /// ank25 
    cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);

    aTolVector = NULL;
    aTolVector = N_VMake_Serial(n,aTolTemp);

    std::cout << "Calling CVodeMalloc " << std::endl; /// ank25  
    CVodeMalloc(cvode_mem,cvodeResid,currentTime,y,CV_SV,rtol,aTolVector);

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
void CVodeWrapper::setIniStep(double istep){
    CVodeSetInitStep(cvode_mem,istep);
}
/*
 *set the max allowed step size
 */
void CVodeWrapper::setMaxStep(double maxStep){
    CVodeSetMaxStep(cvode_mem,maxStep);
}

double& CVodeWrapper::solve(int stopMode){

    int flag;
    do{
        flag = CVode(cvode_mem,maxTime,y,&currentTime,stopMode);
        if(flag < 0){
            std::cout << "Cvode Integration error\n";
        }else{
            CVodeGetDky(cvode_mem,currentTime,1,yPrime);
            calcResNorm();
            reacPtr->report(currentTime,NV_DATA_S(y),resNorm);
        }
    }while(currentTime < maxTime);

    CVodeGetDky(cvode_mem,currentTime,1,yPrime);
    calcResNorm();
    reacPtr->report(currentTime,NV_DATA_S(y),resNorm);
    return resNorm;

}

void CVodeWrapper::solve(int stopMode, double resTol){

    int flag;
    do{
        flag = CVode(cvode_mem,maxTime,y,&currentTime,stopMode);
        if(flag < 0){
            std::cout << "Cvode Integration error\n";
        }else{
            CVodeGetDky(cvode_mem,currentTime,1,yPrime);
            calcResNorm();
            reacPtr->report(currentTime,NV_DATA_S(y),resNorm);
        }
        if(currentTime > maxTime)
            break;
    }while(resNorm > resTol );

}

void CVodeWrapper::solveDAE(int stopMode, double resTol){

    int flag;
    static int iter = 0;
    
    do
    {
        if (iter % 10 == 0)
        {
            //solve the algebraic equation system
            reacPtr->calcFlowField(currentTime, NV_DATA_S(y));
        }
        //call Cvode to solve the ODE
        flag = CVode(cvode_mem, maxTime, y, &currentTime, stopMode);
        
        if (flag < 0)
        {
            std::cout << "Cvode Integration error\n";
        }
        else
        {
            CVodeGetDky(cvode_mem, currentTime, 1, yPrime);
            calcResNorm();
            reacPtr->report(currentTime, NV_DATA_S(y), resNorm);
        }
        
        iter++;
        
    } while (resNorm > resTol);
}

void CVodeWrapper::destroy(){
    CVodeFree(&cvode_mem);
    N_VDestroy(y);
    N_VDestroy(yPrime);
}

void CVodeWrapper::calcResNorm(){
    resNorm = 0;
    double *yp;
    yp = NV_DATA_S(yPrime);
    for(int i=0; i<(eqnSize); i++)   
        resNorm += yp[i]*yp[i];

    resNorm = sqrt(resNorm);
}
