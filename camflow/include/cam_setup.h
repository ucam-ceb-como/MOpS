/* 
 * File:   cam_setup.h
 * Author: vj231
 *
 * Created on 17 February 2009, 10:25
 */

#ifndef _CAM_SETUP_H
#define	_CAM_SETUP_H
#include "cam_boundary.h"
#include "gpc.h"
#include "cam_params.h"
#include "cam_residual.h"
#include "cam_profile.h"
#include "cam_control.h"
#include <vector>
using namespace Sprog;
namespace Camflow{
    class CamSetup : public CamResidual{
    public:
        //return the inlet species mass fractions for the given boundary
        void getInletMassFrac(CamBoundary &cb, vector<doublereal>& fracs);
        //return the inlet temperature
        const doublereal getInletTemperature(CamBoundary &cb);
        //return the inlet flow rate
        const doublereal getInletFlowRate(CamBoundary &cb);
        //return the inlet velocity
        const doublereal getInletVelocity(CamBoundary &cb);
        //return the initial  guess for species mass fractions
        void getInitialGuess(vector<doublereal> &fracs);
        //create the solution vector at the specified boundary
        void createSolnVector(CamBoundary &cb, CamControl &cc, vector<doublereal> &soln);
        //create the solution vector for intermediate soln points
        void createSolnVector(CamBoundary &cb, CamControl &cc,int n1, int n2, vector<doublereal> &soln);
        
    protected:
        CamProfile *profile;
    };
}

#endif	/* _CAM_SETUP_H */

