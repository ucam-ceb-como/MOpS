/* 
 * File:   giant_wrapper.h
 * Author: vj231
 *
 * Created on 24 February 2009, 18:33
 */

#ifndef _GIANT_WRAPPER_H
#define	_GIANT_WRAPPER_H
#include <vector>
#include "cam_residual.h"
#include "cam_params.h"

extern "C"{
   // #include "giant.h"
    #include "nleq.h"
}

namespace Camflow{


    class GiantWrapper{
    public:
        GiantWrapper(){}
        virtual ~GiantWrapper(){};
        void solve(int neq,
                doublereal tol,std::vector<doublereal>& iniGuess,CamResidual& cr);

    private:
        //struct GIANT_FUN *resid;
        //struct GIANT_OPT *optional;
        //struct GIANT_INFO *info;
        struct NLEQ_FUN *resid;
        struct NLEQ_OPT *optional;
        struct NLEQ_INFO *info;
        
    };
}



#endif	/* _GIANT_WRAPPER_H */

