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


std::vector<double> CamSetup::getInletMassFrac(CamBoundary& cb){

     /*
     *Function returns the inlet mass fraction based on the
     *boundary. Base class function called by all the reactor models
     */
    std::vector<double> fracs;
    CamConverter converter;
    std::vector<double> initialFracs = cb.setInletfracs(*camMech_);
    if(cb.getFracType() == cb.MOLE){
        converter.mole2mass(initialFracs,fracs,*camMech_);
        cb.setInletMassfracs(fracs);
    }else{
        fracs = initialFracs;
        cb.setInletMassfracs(fracs);
    }
    return fracs;

}


double CamSetup::getInletTemperature(CamBoundary& cb){
    double T=0;
    //if(admin_.getEnergyModel() == admin_.ISOTHERMAL ){
        T = cb.getTemperature();
    //}else{
    //    T = profile_.getUserDefTemp(0.0);
    //}
    return T;
}

double CamSetup::getInletFlowRate(CamBoundary& cb)
{
    double flow;
    if (cb.getFlowRate() == 0)
    {
        std::vector<double> massfracs = getInletMassFrac(cb);
        camMixture_->SetMassFracs(massfracs);
        camMixture_->SetTemperature(getInletTemperature(cb));
        double avgMolWt = camMixture_->getAvgMolWt();
        rho = opPre*avgMolWt/(R*camMixture_->Temperature());
        flow = rho*cb.getVelocity();
    }
    else
    {
        flow = cb.getFlowRate();
    }

    return flow;
}

double CamSetup::getInletVelocity(CamBoundary& cb){
    double vel;
    double avgMolWt = camMixture_->getAvgMolWt();
    rho = opPre*avgMolWt/(R*camMixture_->Temperature());
    vel= (getInletFlowRate(cb)/rho);
    return vel;
}

/*
 *init species
 */
std::vector<double> CamSetup::initSpecies(CamBoundary& cb)
{
    std::vector<double> soln(nSpc*cellEnd,0);
    std::vector<double> position = reacGeom_.getAxpos();

    profile_.setStartProfile(cb,*camMech_);
    Array2D start = profile_.getStartProfile();

    for(int i=cellBegin; i<cellEnd; i++){
        for(int l=0; l<nSpc; l++){
            soln[i*nSpc+l] = start(i,l);
        }
    }
    return soln;
}
/*
 *init species given 2 inlets
 */
std::vector<double> CamSetup::initSpecies
(
    CamBoundary& left,
    CamBoundary& right
)
{
    std::vector<double> soln(nSpc*cellEnd,0);
    std::vector<double> position = reacGeom_.getAxpos();

    profile_.setStartprofile(left,right,*camMech_);
    Array2D start = profile_.getStartProfile();
    for(int i=cellBegin; i<cellEnd; i++){
        for(int l=0; l<nSpc; l++){
            soln[i*nSpc+l] = start(i,l);
        }
    }
    return soln;
}


/*
 *init mass
 */
std::vector<double> CamSetup::initMassFlow(CamBoundary& cb)
{
    std::vector<double> soln(cellEnd, 0.0);

    std::vector<double> position = reacGeom_.getAxpos();

    double flow = getInletFlowRate(cb);

    for (int i=cellBegin; i<cellEnd; i++)
    {
        soln[i] = flow;
    }

    return soln;
}

/*
 *init temperature
 */
std::vector<double> CamSetup::initTemperature(CamBoundary& cb)
{
    std::vector<double> soln(cellEnd, 0.0);
    std::vector<double> position = reacGeom_.getAxpos();

    double T = getInletTemperature(cb);

    if (admin_.getEnergyModel() == admin_.ISOTHERMAL)
    {
        for (int i=cellBegin; i<cellEnd; i++)
        {
            soln[i] = T;
        }
    }
    else
    {
        for (int i=cellBegin; i<cellEnd; i++)
        {
            T = profile_.getUserDefTemp(position[i]);
            soln[i] = T;
        }
    }
    
    return soln;
}

void CamSetup::initTempGauss(std::vector<double>& soln){
    profile_.setGaussTempProfile(soln);
}


void CamSetup::storeInlet(CamBoundary& cb, inletStruct& ud_inlet)
{
    /*
     *store the inlet properties in the
     *structure to use with the inlet boundary
     */
    std::vector<double> temp = getInletMassFrac(cb);
    camMixture_->SetMassFracs(temp);
    double T = cb.getTemperature();
    camMixture_->SetTemperature(T);
    ud_inlet.T = T;
    ud_inlet.FlowRate = getInletFlowRate(cb);
    ud_inlet.Vel = getInletVelocity(cb);
    ud_inlet.rVelGrad = 0.0;
    ud_inlet.Dens = ud_inlet.FlowRate/ud_inlet.Vel;
    ud_inlet.Species = temp;
    ud_inlet.Dk = camMixture_->getMixtureDiffusionCoeff(opPre);


    std::cout<< "Inlet Boundary:" << "\n"
        << "    Temp      : " << ud_inlet.T << "\n"
        << "    Flow Rate : " << ud_inlet.FlowRate << "\n"
        << "    Velocity  : " << ud_inlet.Vel << "\n"
        << "    Density   : " << ud_inlet.Dens << "\n"
        << "    nDk       : " << ud_inlet.Dk.size() << "\n"
        << "    nSpecies  : " << ud_inlet.Species.size() << std::endl;
}

