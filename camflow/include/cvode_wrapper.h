/* 
 * File:   cvode_wrapper.h
 * Author: vj231
 *
 * Created on 25 February 2009, 17:30
 */

#ifndef _CVODE_WRAPPER_H
#define	_CVODE_WRAPPER_H

#include "cam_params.h"
#include "cam_residual.h"
#include <cmath>
#include <vector>
#include "cvode/cvode.h"
#include "cvode/cvode_band.h"
#include "cvode/cvode_dense.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_types.h"

namespace Camflow{
    class CVodeWrapper{
        void *cvode_mem;
        N_Vector y,yPrime;
        doublereal atol, currentTime, maxTime;
        int eqnSize;
        CamResidual *reacPtr;
        doublereal resNorm;
    public:
        CVodeWrapper();
        ~CVodeWrapper();
        void init(int n, std::vector<doublereal> &solnVec, 
                                        doublereal atol,
                                        doublereal rtol,
                                        doublereal maxIntTime,
                                        int band,
                                        CamResidual &cr, doublereal iniTime=0);
        
        /*
         *additional solver control
         */
        void setIniStep(doublereal istep);
        void setMaxStep(doublereal maxStep);
        doublereal& solve(int stopMode);
        void solve(int stopMode, doublereal resTol);
        void solveDAE(int stopMode, doublereal resTol);
        void calcResNorm();
        void destroy();
    };
}


#endif	/* _CVODE_WRAPPER_H */

