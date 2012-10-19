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
        N_Vector y,yPrime,  aTolVector;
        double atol, currentTime, maxTime;
        int eqnSize;
        CamResidual *reacPtr;
        double resNorm;
    public:
        CVodeWrapper();
        ~CVodeWrapper();
        void init(int n, std::vector<double> &solnVec, 
                                        double atol,
                                        double rtol,
                                        double maxIntTime,
                                        int band,
                                        CamResidual &cr, double iniTime=0);
        
        void initVectorTol(int n, std::vector<double>& solnVec,
										double aTolTemp[],
										double rtol,
										double maxIntTime ,
										int band,
										CamResidual &cr, double iniTime=0);

        /*
         *additional solver control
         */
        void setIniStep(double istep);
        void setMaxStep(double maxStep);
        double& solve(int stopMode);
        void solve(int stopMode, double resTol);
        void solveDAE(int stopMode, double resTol);
        void calcResNorm();
        void destroy();
    };
}


#endif	/* _CVODE_WRAPPER_H */

