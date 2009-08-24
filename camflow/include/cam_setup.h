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
        
        typedef struct{
            vector<doublereal> Species;
            doublereal FlowRate;
            doublereal Vel;
            doublereal rVelGrad;
            doublereal Dens;
            doublereal T;
            vector<doublereal> Dk;
            vector<doublereal> jk;
        } inletStruct;

        CamSetup(){}
        virtual ~CamSetup(){}
        
        //return the inlet species mass fractions for the given boundary
        void getInletMassFrac(CamBoundary &cb, vector<doublereal>& fracs);
        //return the inlet temperature
        const doublereal getInletTemperature(CamBoundary &cb);
        //return the inlet flow rate
        const doublereal getInletFlowRate(CamBoundary &cb);
        //return the inlet velocity
        const doublereal getInletVelocity(CamBoundary &cb);

        /*
         *init species
         */
        void initSpecies(CamBoundary &cb, CamControl &cc, vector<doublereal>& soln);
        /*
         *initialize the species vector for a counter flow flame
         */
        void initSpecies(CamBoundary &left, CamBoundary &right,
                                CamControl &cc, vector<doublereal>& soln);
        /*
         *init mass flow
         */
        void initMassFlow(CamBoundary &cb, CamControl &cc, vector<doublereal> &soln);
        /*
         *init temperature
         */
        void initTemperature(CamBoundary &cb, CamControl &cc, vector<doublereal> &soln);
        /*
         *init temperature based on a gauss profile
         */
        void initTempGauss(vector<doublereal> &soln);
        /*
         *store the inlet conditions
         */
        void storeInlet(CamBoundary &cb, inletStruct& ud_inlet);
        
    protected:
        CamProfile *profile;

    };
}

#endif	/* _CAM_SETUP_H */

