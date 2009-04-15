/* 
 * File:   giant_wrapper.h
 * Author: vj231
 *
 * Created on 24 February 2009, 18:33
 */

#ifndef _GIANT_WRAPPER_H
#define	_GIANT_WRAPPER_H
#include "cam_params.h"
#include "giant.h"
#include "cam_residual.h"
namespace Camflow{
    class GiantWrapper{
        int nEq;
        doublereal *y;
        GIANT_FUN fun;
        GIANT_OPT opt;
        GIANT_INFO info;
    public:
        GiantWrapper(){}
        ~GiantWrapper(){}
        void init(int n, vector<doublereal>& solnVec,
                doublereal rtol,
                CamResidual &cr);
        void solve();
    };
}


#endif	/* _GIANT_WRAPPER_H */

