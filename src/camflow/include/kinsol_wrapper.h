/* 
 * File:   kinsol_wrapper.h
 * Author: vj231
 *
 * Created on 24 February 2009, 16:36
 */

#ifndef _KINSOL_WRAPPER_H
#define	_KINSOL_WRAPPER_H

#include "cam_params.h"
#include "cam_residual.h"
#include <kinsol/kinsol.h>
#include "kinsol/kinsol_band.h"
#include "kinsol/kinsol_dense.h"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
namespace Camflow{
    class KinsolWrapper{

    public:
        KinsolWrapper(){}
        ~KinsolWrapper(){}
        void init(int n, vector<double> &solnVec, double rtol, int band, CamResidual &cr);
        void solve();
        void destroy();

    private:
        N_Vector y, scale;
        double fnormtol, fnorm;
        void *kmem;
        int mset, msubset, flag;
    };
}


#endif	/* _KINSOL_WRAPPER_H */

