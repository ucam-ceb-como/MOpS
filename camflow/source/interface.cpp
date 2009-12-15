
#include <vector>
#include "cam_residual.h"
#include "cam_configuration.h"
#include "cam_control.h"
#include "interface.h"
#include "flamelet.h"
using namespace Camflow;


/*
 *this interface is for external programs to call any of the models
 *within camflow. The interface constructor does read the input files.
 *all the files required for the stand alone program is also required
 *here. Once the solution is complete, the interface object can be
 *queried to get the information on dependent variables
 */
Interface::Interface(){
    
    std::string fChem("chem.inp");
    std::string fThermo("therm.dat");
    std::string fTrans("tran.dat");
    std::string fCamFlow("camflow.xml");

    try{
        cm.readInput(fCamFlow,cc,cg,convert,ca,cb,cp,config,cSoot);
    }catch(CamError &ce){
        std::cout << ce.errorMessge;
    }
    

    //read mechanism, thermo and trasnport data
    IO::MechanismParser::ReadChemkin(fChem,mech,fThermo,fTrans);

    //get the number of species
    nSpecies = mech.SpeciesCount();
    //create the mixture 
    Thermo::Mixture mix(mech.Species());
    const SpeciesPtrVector *spv = mix.Species();
    /*
     *populate the species names
     */
    for(int l=0; l<nSpecies; l++){
        speciesNames.push_back((*spv)[l]->Name());
    }

    flmlt = NULL;
}

/*
 *Interface for external drivers which passes a reference to the
 *mechanism object
 */
Interface::Interface(Mechanism& mech_in,
        std::vector<doublereal>& dz,
        std::vector<Thermo::Mixture>& cstrs,
        void* rModel, const doublereal sdr
                                                                ){

    std::string fCamFlow("camflow.xml");    
    try{
        cm.readInput(fCamFlow,cc,cg,convert,ca,cb,cp,config,cSoot);
    }catch(CamError &ce){
        throw;
    }
    /*
     *reset the reactor mesh to the passed in values
     */
    if(dz.size() != 0)cg.setGeometry(dz);
    
    CamSoot cs;

    model = (CamResidual*)rModel;
    if(sdr==0){
        model->solve(cc,ca,cg,cp,config,cs,mech_in);
        resetMixtures(cstrs);
    }else{
        /**
         * assumes that the call is made for flamelet
         * model if the scalar dissipation rate is
         * non-zero
         */
        model->setExternalScalarDissipationRate(sdr);
        model->solve(cc,ca,cg,cp,config,cs,mech_in);
    }


}

/*
 *reset the mixtures with the newly evaluated rties
 */
void Interface::resetMixtures(std::vector<Thermo::Mixture>& cstrs){
     //Get the species mass fractions
    Array2D massFracs;
    model->getSpeciesMassFracs(massFracs);

    //storage for density, velocity, and temperature
    std::vector<doublereal> density, vel, temp;

    //Get the density
    model->getDensityVector(density);
    //Get the velocity
    model->getVelocity(vel);
    //Get the temperature
    model->getTemperatureVector(temp);

    int nCells = cg.getnCells();
    if(cstrs.size() >0){
        if(cstrs.size() != nCells-2)
            throw("size of mixtures is not consistant with the grid\n");
        int nSp = mech.SpeciesCount();
        for(int i=0; i<nCells-2;i++){
            std::vector<doublereal> mf;
            for(int l=0; l<nSp; l++){
                mf.push_back(massFracs(i+1,l));
            }
            cstrs[i].SetMassFracs(mf);
            cstrs[i].SetMassDensity(density[i+1]);
            cstrs[i].SetTemperature(temp[i+1]);
            cstrs[i].SetVelocity(vel[i+1]);
        }
    }else{
        /*
         *create the cstrs
         */
        Thermo::Mixture mix(mech.Species());
        int nSp = mech.SpeciesCount();
        for(int i=0; i<nCells-2;i++){
            std::vector<doublereal> mf;
            for(int l=0; l<nSp; l++){
                mf.push_back(massFracs(i+1,l));
            }
            mix.SetMassFracs(mf);
            mix.SetMassDensity(density[i+1]);
            mix.SetTemperature(temp[i+1]);
            mix.SetVelocity(vel[i+1]);
            cstrs.push_back(mix);
        }
    }

}

/*
 *solve the reactor problem
 */
void Interface::solve(std::vector<Thermo::Mixture>& cstrs,
        const std::vector<doublereal>& dz,
        const std::vector< std::vector<doublereal> >& initalSource,
        const std::vector< std::vector<doublereal> >& finalSource,
        CamControl& ccObj,
        CamConfiguration& confObj,
        Mechanism& mech_in,
        void* reactorModel,
        const doublereal sdr){


    cc = ccObj;    
    config = confObj;
    if(cstrs.size() != dz.size()){
        throw CamError("Mismatch between the number of mixtures passed and the cell geometry\n");
    }else{
        //set the rector geometry passed in by the external code to the
        //geometry object
        cg.setGeometry(dz);
    }
    /*
     *solve the reactor problem
     */
    model = (CamResidual*)reactorModel;
    if(sdr==0){
        model->solve(cstrs,initalSource,finalSource,mech_in,ccObj,ca,cg,cp);
        resetMixtures(cstrs);
    }       

}

/*
 *return the number of species
 */
const int Interface::getNumberOfSpecies() const {
    return nSpecies;
}
/*
 *return the argument vector with the species names
 */
void Interface::getSpeciesNames(std::vector<std::string>& names){
    names.clear();
    names = speciesNames;
}

/**
 *  This function is called by the external code that
 *  passes a scalar dissipation rate with time history
 */
void Interface::flamelet(const std::vector<doublereal>& sdr, const std::vector<doublereal>& intTime, bool continuation, bool lnone){

    if(sdr.size() != intTime.size())
        throw CamError("Mismatch in the size of SDR and TIME vector\n");

    if(flmlt == NULL ) flmlt = new FlameLet();    
    //Set the time history of the scalar dissipation rate
    flmlt->setRestartTime(intTime[0]);
    flmlt->setExternalScalarDissipationRate(intTime,sdr);

    // Build up a vector of zero soot volume fractions
    std::vector<doublereal> zeroSoot(cg.getnCells(), 0.0);
    flmlt->setExternalSootVolumeFraction(zeroSoot);

    try{
        const size_t len = sdr.size();
        flamelet(sdr[len-1], intTime[len-1],continuation,lnone);
    }catch(CamError& ce){
        throw;
    }

}

/**
 *  This function is called by the external code that
 *  passes a scalar dissipation rate with time history
 *  and a spatial profile of soot volume fraction.
 */
void Interface::flameletWithSoot(const std::vector<doublereal>& soot_fv, const std::vector<doublereal>& sdr, 
                                 const std::vector<doublereal>& intTime, bool continuation, bool lnone){

    if(sdr.size() != intTime.size())
        throw CamError("Mismatch in the size of SDR and TIME vector\n");

    if(flmlt == NULL ) flmlt = new FlameLet();    
    //Set the time history of the scalar dissipation rate
    flmlt->setRestartTime(intTime[0]);
    flmlt->setExternalScalarDissipationRate(intTime,sdr);

    //@todo check length of the vector
    flmlt->setExternalSootVolumeFraction(soot_fv);

    try{
        const size_t len = sdr.size();
        flamelet(sdr[len-1], intTime[len-1],continuation,lnone);
    }catch(CamError& ce){
        throw;
    }

}


/*
 *this function is the interface for calling the flamelet code
 *for the implementation of interactive falmelets. The function
 *takes the scalar dissipation rate as input. The integration time
 *cam be passed as an optional argument. If its not passed, the
 *integration will be carried upto the time specified in the input
 *file/default value. If steady state is obtained before the specified
 *integration time, the program return with the converged solution.
 */
void Interface::flamelet(doublereal sdr, doublereal intTime, bool continuation, bool lnone){
    
    if(intTime!=0)cc.setMaxTime(intTime);
    if(flmlt == NULL ) flmlt = new FlameLet();
    if(!lnone) flmlt->setLewisNumber(FlameLet::LNNONE);
    flmlt->setExternalScalarDissipationRate(sdr);
    if(sdr==0){
        std::cout << "Passed in value of SDR is zero\n SDR automatic calcilation activated\n" ;
    }
    try{
        if(!continuation){
            flmlt->solve(cc,ca,cg,cp,mech,true);
        }else{
            flmlt->restart(cc);
        }
        /*
         *store the results for lookup
         */
        flmlt->getDensityVector(rhoVector);
        flmlt->getSpeciesMassFracs(spMassFracs);
        flmlt->getTemperatureVector(TVector);
        flmlt->getIndepedantVar(indVar);
        flmlt->getViscosityVector(muVector);
        flmlt->getSpecificHeat(spHeat);
        flmlt->getThermalConductivity(lambda);
        flmlt->getDiffusionCoefficient(mDiff);
        flmlt->getVelocity(mVelocity);
        stMixtureFrac = flmlt->stoichiometricMixtureFraction();
        std::cout << "Size of thermal conductivity " << lambda.size() << std::endl;
    }catch(CamError &ce){
        throw ;
    }
}
/*
 *return the stoichiometric mixture fraction
 */
const doublereal Interface::getStMixtureFrac(){
    return stMixtureFrac;
}
/*
 *
 *return the density
 */
const doublereal Interface::getDensity(const doublereal axpos){

    doublereal dens = getVariableAt(axpos,rhoVector);
    return dens;
}

/*!
 *@param[in]    spIndex     Index of species for which mass fractions requested
 *
 *@return   Vector of mass fractions at all the independent variable points
 */
const std::vector<doublereal> Interface::getMassFracsBySpecies(const int spIndex) const {
    // Find the length of the mass fraction profile and create an empty
    // vector with this much space
    const size_t len = indVar.size();
    std::vector<doublereal> mf(len);

    for(size_t i=0; i<len; i++){
        mf[i] = spMassFracs(i,spIndex);
    }

    return mf;
}

/*!
 *@param[in]    indVarIndex     Mass fractions requested at indVar[indVarIndex]
 *
 *@return   Vector of all mass fractions at one independent variable point
 */
const std::vector<doublereal> Interface::getMassFracsByPoint(const int indVarIndex) const {
    std::vector<doublereal> mf(nSpecies);

    for(int i = 0; i < nSpecies; i++){
        mf[i] = spMassFracs(indVarIndex, i);
    }

    return mf;
}
/*
 *return the mass fractions
 */
const doublereal Interface::getMassFrac(const int spIndex, const doublereal axpos){
    doublereal massfrac = getVariableAt(axpos, getMassFracsBySpecies(spIndex));
    return massfrac;
}
/*
 *return the temperature
 */
const doublereal Interface::getTemperature(const doublereal axpos){
    
    doublereal temp = getVariableAt(axpos,TVector);
    return temp;
}
/*
 *return the viscosity
 */
const doublereal Interface::getViscosity(const doublereal axpos){
    
    doublereal temp = getVariableAt(axpos,muVector);
    return temp;
}

/**
 *  Return the specific heat
 */
const doublereal Interface::getSpecificHeat(const doublereal axPos){
    doublereal temp = getVariableAt(axPos,spHeat);
    return temp;
}

/**
 *  Return the thermal conductivity
 */
const doublereal Interface::getThermalConductivity(const doublereal axPos){
    doublereal temp = getVariableAt(axPos,lambda);
    return temp;
}
/**
 *  Return a vector of diffusion coefficients
 */
const std::vector<doublereal> Interface::getDiffusionCoefficients(const doublereal axPos){
    std::vector<doublereal> diff, rmDiff;
    int len = indVar.size();
    rmDiff.clear();
    for(int k=0; k<nSpecies; k++){
        diff.clear();
        //get the diffusion coefficient of species k for all mixture fraction
        for(int i=0; i<len; i++){
            diff.push_back(mDiff(i,k));
        }
        doublereal spDiff = getVariableAt(axPos,diff);
        rmDiff.push_back(spDiff);
    }

    return rmDiff;

}
/*
 *private function to do the interpolation of the solution variables
 */
doublereal Interface::getVariableAt(const doublereal& pos,
                                    const std::vector<doublereal>& var) const{

    doublereal vu, vl, xu, xl, temp=0.0;
    int len = indVar.size();



    for(int i=0; i<len; i++){
        if(pos == indVar[i]) {
            temp= var[i];
            break;
        }else if( i>0 && (pos > indVar[i-1]) && (pos < indVar[i]) ){
            vu = var[i];
            xu = indVar[i];
            
            vl = var[i-1];
            xl = indVar[i-1];
            
            doublereal slope = (vu-vl)/(xu-xl);
            doublereal intersect = vu- (slope*xu);
            
            temp= slope*pos + intersect;
            break;
        }
    }

    return temp;

}

/*
 *return the handle to CamAdmin
 */
CamAdmin& Interface::getCamAdmin(){
    return ca;
}
/*
 *return the boundary object
 */
CamBoundary& Interface::getCamBoundary(){
    return cb;
}
/*
 *return the object for controlling solver
 */
CamControl&  Interface::getCamControl(){
    return cc;
}
/*
 *return the geometry object
 */
CamGeometry& Interface::getCamGeometry(){
    return cg;
}

/*
 *return the profile object which handles the
 *user defined profiles
 */

CamProfile& Interface::getCamProfile(){
    return cp;
}
/*
 *return the configuration object
 */
CamConfiguration& Interface::getCamConfiguration(){
    return config;
}

