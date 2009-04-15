
#include <vector>


#include "cam_setup.h"
#include "cam_converter.h"
using namespace Camflow;

void CamSetup::getInletMassFrac(CamBoundary& cb, vector<doublereal>& fracs){

     /*
     *Function returns the inlet mass fraction based on the
     *boundary. Base class function called by all the reactor models
     */

    CamConverter converter;
    vector<doublereal> initialFracs = cb.setInletfracs(*camMech);
    if(cb.getFracType() == cb.MOLE){
        converter.mole2mass(initialFracs,fracs,*camMech);
        cb.setInletMassfracs(fracs);
    }else{
        fracs = initialFracs;
        cb.setInletMassfracs(fracs);
    }

}


const doublereal CamSetup::getInletTemperature(CamBoundary& cb){
    doublereal T=0;
    if(admin->getEnergyModel() == admin->ISOTHERMAL ){
        T = cb.getTemperature();
    }else{
        T = profile->getUserDefTemp(0.0);
    }
    return T;
}

const doublereal CamSetup::getInletFlowRate(CamBoundary& cb){

    doublereal flow;
    if(cb.getFlowRate() == 0){
        doublereal avgMolWt = camMixture->getAvgMolWt();
        rho = opPre*avgMolWt/(R*camMixture->Temperature());
        flow = rho*cb.getVelocity();
    }else{
        flow = cb.getFlowRate();
    }

    return flow;
}

const doublereal CamSetup::getInletVelocity(CamBoundary& cb){
    doublereal vel;
    doublereal avgMolWt = camMixture->getAvgMolWt();
    rho = opPre*avgMolWt/(R*camMixture->Temperature());
    vel= (getInletFlowRate(cb)/rho);
    return vel;
}


void CamSetup::getInitialGuess(vector<doublereal>& fracs){

    fracs = profile->getInitialSpeciesGuess(*camMech);
}

void CamSetup::createSolnVector(CamBoundary& cb, CamControl &cc, vector<doublereal>& soln){
    /*
     *This function generates the solution vector
     *at the inlet in the following order
     *Species, temperature, mass flow
     */
    getInletMassFrac(cb, soln);
    rTol.resize(soln.size(),cc.getSpeciesRelTol());
    aTol.resize(soln.size(),cc.getSpeciesAbsTol());

    doublereal flow = getInletFlowRate(cb);
    soln.push_back(flow);
    rTol.push_back(cc.getFlowRelTol());
    aTol.push_back(cc.getFlowAbsTol());

    doublereal T = getInletTemperature(cb);
    soln.push_back(T);
    rTol.push_back(cc.getTempRelTol());
    aTol.push_back(cc.getTempAbsTol());
}

void CamSetup::createSolnVector(CamBoundary &cb, CamControl &cc, int n1, int n2, vector<doublereal>& soln){
    /*
     *This function generates the solution
     *vector and the associated tolerances for the interior cells.
     *The vector elements are arranged in the following order
     *Species,massflow,Temperature
     */
    vector<doublereal> initial, position;
    doublereal flow = getInletFlowRate(cb);
    doublereal T = cb.getTemperature();

    doublereal sTol_r = cc.getSpeciesRelTol();
    doublereal sTol_a = cc.getSpeciesAbsTol();
    doublereal TTol_r = cc.getTempRelTol();
    doublereal TTol_a = cc.getTempAbsTol();
    doublereal cTol_r = cc.getFlowRelTol();
    doublereal cTol_a = cc.getFlowAbsTol();

    getInitialGuess(initial);
    position = reacGeom->getAxpos();


    for(int i = n1; i< n2; i++){
        for(int l=0; l<nSpc; l++ ){
            soln.push_back(initial[l]);
            rTol.push_back(sTol_r);
            aTol.push_back(sTol_a);
        }
 
        //mass flow assignment
        soln.push_back(flow);
        rTol.push_back(cTol_r);
        aTol.push_back(cTol_a);


        //temperature assignment
        if(admin->getEnergyModel() == admin->ISOTHERMAL){
            soln.push_back(T);
        }else{
            T = profile->getUserDefTemp(position[i]);
            soln.push_back(T);

        }
        rTol.push_back(TTol_r);
        aTol.push_back(TTol_a);
      
    }

}
