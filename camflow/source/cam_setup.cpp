
#include "cam_profile.h"


#include <vector>


#include "cam_setup.h"
#include "cam_converter.h"
using namespace Camflow;

CamSetup::CamSetup()
:
profile(NULL)
{}

CamSetup::~CamSetup()
{

//    if (profile != NULL) delete profile;

}


void CamSetup::getInletMassFrac(CamBoundary& cb, std::vector<doublereal>& fracs){

     /*
     *Function returns the inlet mass fraction based on the
     *boundary. Base class function called by all the reactor models
     */

    CamConverter converter;
    std::vector<doublereal> initialFracs = cb.setInletfracs(*camMech);
    if(cb.getFracType() == cb.MOLE){
        converter.mole2mass(initialFracs,fracs,*camMech);
        cb.setInletMassfracs(fracs);
    }else{
        fracs = initialFracs;
        cb.setInletMassfracs(fracs);
    }

}


doublereal CamSetup::getInletTemperature(CamBoundary& cb){
    doublereal T=0;
    if(admin->getEnergyModel() == admin->ISOTHERMAL ){
        T = cb.getTemperature();
    }else{
        T = profile->getUserDefTemp(0.0);
    }
    return T;
}

doublereal CamSetup::getInletFlowRate(CamBoundary& cb){

    doublereal flow;
    if(cb.getFlowRate() == 0){
        std::vector<doublereal> massfracs;
        getInletMassFrac(cb,massfracs);
        camMixture->SetMassFracs(massfracs);
        camMixture->SetTemperature(getInletTemperature(cb));
        doublereal avgMolWt = camMixture->getAvgMolWt();
        rho = opPre*avgMolWt/(R*camMixture->Temperature());
        flow = rho*cb.getVelocity();
    }else{
        flow = cb.getFlowRate();
    }

    return flow;
}

doublereal CamSetup::getInletVelocity(CamBoundary& cb){
    doublereal vel;
    doublereal avgMolWt = camMixture->getAvgMolWt();
    rho = opPre*avgMolWt/(R*camMixture->Temperature());
    vel= (getInletFlowRate(cb)/rho);
    return vel;
}

/*
 *init species
 */
void CamSetup::initSpecies(CamBoundary& cb, CamControl& cc,
                                    std::vector<doublereal>& soln){

    soln.resize(nSpc*cellEnd,0);
    std::vector<doublereal>  position;
    position = reacGeom->getAxpos();
    profile->setStartProfile(cb,*camMech);
    Array2D start = profile->getStartProfile();

    for(int i=cellBegin; i<cellEnd; i++){
        for(int l=0; l<nSpc; l++){
            soln[i*nSpc+l] = start(i,l);
        }
    }

}
/*
 *init species given 2 inlets
 */
void CamSetup::initSpecies(CamBoundary& left, CamBoundary& right,
                    CamControl& cc, std::vector<doublereal>& soln){

    soln.resize(nSpc*cellEnd,0);
    std::vector<doublereal>  position;
    position = reacGeom->getAxpos();
    profile->setStartprofile(left,right,*camMech);
    Array2D start = profile->getStartProfile();
    for(int i=cellBegin; i<cellEnd; i++){
        for(int l=0; l<nSpc; l++){
            soln[i*nSpc+l] = start(i,l);
        }
    }

}
/*
 *init mass
 */
void CamSetup::initMassFlow(CamBoundary& cb, CamControl& cc,
                                        std::vector<doublereal>& soln){

    soln.resize(cellEnd,0);
    std::vector<doublereal>  position;
    position = reacGeom->getAxpos();
    doublereal flow = getInletFlowRate(cb);
    for(int i=cellBegin; i<cellEnd; i++)
        soln[i] = flow;

}

/*
 *init temperature
 */
void CamSetup::initTemperature(CamBoundary& cb, CamControl& cc,
                        std::vector<doublereal>& soln){

    soln.resize(cellEnd,0);
    std::vector<doublereal>  position;
    position = reacGeom->getAxpos();
    doublereal T = getInletTemperature(cb);



    if(admin->getEnergyModel() == admin->ISOTHERMAL){
        for(int i=cellBegin; i<cellEnd; i++)
            soln[i] = T;
    }else{
        for(int i=cellBegin; i<cellEnd; i++){
            T = profile->getUserDefTemp(position[i]);
            soln[i] = T;

        }
    }

}

void CamSetup::initTempGauss(std::vector<doublereal>& soln){
    profile->setGaussTempProfile(soln);
}



void CamSetup::storeInlet(CamBoundary& cb, inletStruct& ud_inlet){
    /*
     *stote the inlet properties in the
     *structure to use with the inlet boundary
     */
    std::vector<doublereal> temp;
    getInletMassFrac(cb,temp);
    camMixture->SetMassFracs(temp);
    doublereal T = getInletTemperature(cb);
    camMixture->SetTemperature(T);
    ud_inlet.T = T;
    ud_inlet.FlowRate = getInletFlowRate(cb);
    ud_inlet.Vel = getInletVelocity(cb);
    ud_inlet.rVelGrad = 0.0;
    ud_inlet.Dens = ud_inlet.FlowRate/ud_inlet.Vel;
    ud_inlet.Species = temp;
    ud_inlet.Dk = camMixture->getMixtureDiffusionCoeff(opPre);

}


void CamSetup::storeObjects(CamControl& cc,
                            CamAdmin& ca,
                            CamGeometry& cg,
                            CamProfile& cp,
                            CamBoundary& cb,
                            Mechanism& mech){

    camMech = &mech;
    control = &cc;
    spv = camMixture->Species();
    profile = &cp;
    admin = &ca;
    reacGeom = &cg;
    CamBoundary cblocal;
    admin->getLeftBoundary(cblocal);
    cb = cblocal;
    opPre = ca.getPressure();
    profile->setGeometryObj(cg);
}

