#include "cam_setup.h"

using namespace Camflow;

CamSetup::CamSetup
(
    CamAdmin& ca,
    CamConfiguration& config,
    CamControl& cc,
    CamGeometry& cg,
    CamProfile& cp,
    CamSoot& cs,
    Mechanism& mech
)
:
  CamResidual(ca, config, cc, cg, cp, cs, mech),
  profile_(cp)
{}

CamSetup::~CamSetup()
{}


void CamSetup::getInletMassFrac(CamBoundary& cb, std::vector<doublereal>& fracs){

     /*
     *Function returns the inlet mass fraction based on the
     *boundary. Base class function called by all the reactor models
     */

    CamConverter converter;
    std::vector<doublereal> initialFracs = cb.setInletfracs(*camMech_);
    if(cb.getFracType() == cb.MOLE){
        converter.mole2mass(initialFracs,fracs,*camMech_);
        cb.setInletMassfracs(fracs);
    }else{
        fracs = initialFracs;
        cb.setInletMassfracs(fracs);
    }

}


doublereal CamSetup::getInletTemperature(CamBoundary& cb){
    doublereal T=0;
    if(admin_.getEnergyModel() == admin_.ISOTHERMAL ){
        T = cb.getTemperature();
    }else{
        T = profile_.getUserDefTemp(0.0);
    }
    return T;
}

doublereal CamSetup::getInletFlowRate(CamBoundary& cb){

    doublereal flow;
    if(cb.getFlowRate() == 0){
        std::vector<doublereal> massfracs;
        getInletMassFrac(cb,massfracs);
        camMixture_->SetMassFracs(massfracs);
        camMixture_->SetTemperature(getInletTemperature(cb));
        doublereal avgMolWt = camMixture_->getAvgMolWt();
        rho = opPre*avgMolWt/(R*camMixture_->Temperature());
        flow = rho*cb.getVelocity();
    }else{
        flow = cb.getFlowRate();
    }

    return flow;
}

doublereal CamSetup::getInletVelocity(CamBoundary& cb){
    doublereal vel;
    doublereal avgMolWt = camMixture_->getAvgMolWt();
    rho = opPre*avgMolWt/(R*camMixture_->Temperature());
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
    position = reacGeom_.getAxpos();
    profile_.setStartProfile(cb,*camMech_);
    Array2D start = profile_.getStartProfile();

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

    soln.resize(nSpc*cellEnd,0.0);
    std::vector<doublereal>  position;
    position = reacGeom_.getAxpos();
    profile_.setStartprofile(left,right,*camMech_);
    Array2D start = profile_.getStartProfile();
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
    position = reacGeom_.getAxpos();
    doublereal flow = getInletFlowRate(cb);
    for(int i=cellBegin; i<cellEnd; i++)
        soln[i] = flow;

}

/*
 *init temperature
 */
void CamSetup::initTemperature
(
    CamBoundary& cb,
    CamControl& cc,
    std::vector<doublereal>& soln
)
{

    soln.resize(cellEnd,0);
    std::vector<doublereal>  position;
    position = reacGeom_.getAxpos();
    doublereal T = getInletTemperature(cb);

    if(admin_.getEnergyModel() == admin_.ISOTHERMAL)
    {
        for (int i=cellBegin; i<cellEnd; i++)
            soln[i] = T;
    }
    else
    {
        for (int i=cellBegin; i<cellEnd; i++)
        {
            T = profile_.getUserDefTemp(position[i]);
            soln[i] = T;
        }
    }

}

void CamSetup::initTempGauss(std::vector<doublereal>& soln){
    profile_.setGaussTempProfile(soln);
}



void CamSetup::storeInlet(CamBoundary& cb, inletStruct& ud_inlet){
    /*
     *stote the inlet properties in the
     *structure to use with the inlet boundary
     */
    std::vector<doublereal> temp;
    getInletMassFrac(cb,temp);
    camMixture_->SetMassFracs(temp);
    doublereal T = getInletTemperature(cb);
    camMixture_->SetTemperature(T);
    ud_inlet.T = T;
    ud_inlet.FlowRate = getInletFlowRate(cb);
    ud_inlet.Vel = getInletVelocity(cb);
    ud_inlet.rVelGrad = 0.0;
    ud_inlet.Dens = ud_inlet.FlowRate/ud_inlet.Vel;
    ud_inlet.Species = temp;
    ud_inlet.Dk = camMixture_->getMixtureDiffusionCoeff(opPre);

}


void CamSetup::storeObjects(CamControl& cc,
                            CamAdmin& ca,
                            CamGeometry& cg,
                            CamProfile& cp,
                            CamBoundary& cb,
                            Mechanism& mech){

    //camMech = &mech;
    //control = &cc;
    //spv = camMixture->Species();
    //profile_ = &cp;
    //admin = &ca;
    //reacGeom = &cg;
    //CamBoundary cblocal;
    //admin_.getLeftBoundary(cblocal);
    //cb = cblocal;
    //opPre = ca.getPressure();
    //profile_.setGeometryObj(cg);
}

