/* 
 * File:   ida_wrapper.h
 * Author: vj231
 *
 * Created on 03 March 2009, 10:30
 */

#ifndef _IDA_WRAPPER_H
#define	_IDA_WRAPPER_H
#include "cam_params.h"
#include "cam_residual.h"
#include "ida/ida.h"
#include "ida/ida_band.h"
#include "ida/ida_dense.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_types.h"


namespace Camflow{
    class IDAWrapper{
        void *ida_mem;
        N_Vector y0, yp0;
        double atol, currentTime, maxTime;
        int flag;
        CamResidual *reacPtr;
    public:
        IDAWrapper(){}
        ~IDAWrapper(){}
        void init(int n, vector<double> &solnVec, double tol,
                                        double maxIntTime,
                                        int band,
                                        CamResidual &cr);

        void solve();

    };
}

#endif	/* _IDA_WRAPPER_H */

