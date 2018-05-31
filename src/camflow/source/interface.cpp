/*!
 * \file   interface.cpp
 * \author V. Janardhanan
 *
 * \brief Interface for coupling flamelet calculations to external codes.
 *
 *  Copyright (C) 2009 Vinod Janardhanan.
 *

 Licence:
    This file is part of "camflow".

    brush is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Prof Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
 */

#include <vector>
#include "cam_residual.h"
#include "cam_configuration.h"
#include "cam_control.h"
#include "interface.h"
#include "flamelet.h"
#include "linear_interpolator.hpp"

using namespace Camflow;

/*
 *this interface is for external programs to call any of the models
 *within camflow. The interface constructor does read the input files.
 *all the files required for the stand alone program is also required
 *here. Once the solution is complete, the interface object can be
 *queried to get the information on dependent variables
 */
Interface::Interface()
:
cp(cg),
model(NULL),
speciesPointerVector(NULL),
flmlt(NULL)
{

    std::string fChem("chem.inp");
    std::string fThermo("therm.dat");
    std::string fTrans("tran.dat");
    std::string fCamFlow("camflow.xml");

    try{
        cm.readInput(fCamFlow,cc,cg,convert,ca,cb,cp,config,cSoot);
    }catch(CamError &ce){
        std::cout << "ERROR! " << ce.errorMessage;
        std::exit(EXIT_FAILURE);
    }

    //read mechanism, thermo and trasnport data
    Sprog::IO::MechanismParser::ReadChemkin(fChem,mech,fThermo,1.0,fTrans);

    //get the number of species
    nSpecies = mech.SpeciesCount();

    //get the number of moments
    nMoments = cSoot.getNumMoments();

    //create the mixture
    Thermo::Mixture mix(mech.Species());
    speciesPointerVector = mix.Species();

    /*
     *populate the species names
     */
    for(int l=0; l<nSpecies; l++){
        speciesNames.push_back((*speciesPointerVector)[l]->Name());
    }

}

/*
 *Interface for external drivers which passes a reference to the
 *mechanism object
 */
Interface::Interface(Mechanism& mech_in,
        std::vector<double>& dz,
        std::vector<Thermo::Mixture>& cstrs,
        void* rModel, const double sdr
)
:
cp(cg),
model(NULL),
speciesPointerVector(NULL),
flmlt(NULL)
{

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
       // model->solve(cc,ca,cg,cp,config,cs,mech_in);
        resetMixtures(cstrs);
    }else{
        /**
         * assumes that the call is made for flamelet
         * model if the scalar dissipation rate is
         * non-zero
         */
        model->setExternalSDR(sdr);
       // model->solve(cc,ca,cg,cp,config,cs,mech_in);
    }
}


//! Destructor.
Interface::~Interface()
{

    if (flmlt != NULL) delete flmlt;
    if (model != NULL) delete model;
    flmlt = NULL;
    model = NULL;

}


/*
 *reset the mixtures with the newly evaluated rties
 */
void Interface::resetMixtures(std::vector<Thermo::Mixture>& cstrs){
     //Get the species mass fractions
    Array2D massFracs;
    model->getSpeciesMassFracs(massFracs);

    //storage for density, velocity, and temperature
    std::vector<double> density, vel, temp;

    //Get the density
    model->getDensityVector(density);
    //Get the velocity
    model->getVelocity(vel);
    //Get the temperature
    model->getTemperatureVector(temp);

    int nCells = cg.getnCells();
    
    if(cstrs.size() >0){
        size_t neededSize = nCells-2;
        if(cstrs.size() != neededSize)
            throw("size of mixtures is not consistant with the grid\n");
        int nSp = mech.SpeciesCount();
        for(int i=0; i<nCells-2;i++){
            std::vector<double> mf;
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
            std::vector<double> mf;
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
        const std::vector<double>& dz,
        const std::vector< std::vector<double> >& initalSource,
        const std::vector< std::vector<double> >& finalSource,
        CamControl& ccObj,
        CamConfiguration& confObj,
        Mechanism& mech_in,
        void* reactorModel,
        const double sdr){


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
int Interface::getNumberOfSpecies() const {
    return nSpecies;
}

/*
 *return the number of moments
 */
int Interface::getNumberOfMoments() const {
    return nMoments;
}



/*
 *return the number of reactions in the mechanism
 */
int Interface::getNumberOfReactions() const {
    return mech.ReactionCount();
}

/*
 *return the argument vector with the species names
 */
std::vector<std::string> Interface::getSpeciesNames(){
    return speciesNames;
}

/*!
 * Stores the results for lookup by the CFD program.
 */
void Interface::getFlameletVariables(FlameLet* const flmlt)
{

    flmlt->getDensityVector(rhoVector);
    flmlt->getSpeciesMassFracs(spMassFracs);
    flmlt->getMoments(sootMoments);		// ank25 added
    flmlt->getMomentsWdot(sootMomentsWdot);		// ank25 added for ELFM
    flmlt->getTemperatureVector(TVector);
    flmlt->getIndepedantVar(indVar);
    flmlt->getViscosityVector(muVector);
    flmlt->getSpecificHeat(spHeat);
    flmlt->getThermalConductivity(lambda);
    flmlt->getDiffusionCoefficient(mDiff);
    flmlt->getVelocity(mVelocity);
    flmlt->getAverageMolarWeight(avgMolWtVector);
    flmlt->getWdotA4interface(wdotA4);
    flmlt->getSootAverageDiameterVector(sootAverageDiameterVector);
    flmlt->getSootDispersionVector(sootDispersionVector);
    flmlt->getSootSurfaceAreaVector(sootSurfaceAreaVector);
    flmlt->getSootVolumeFractionVector(sootVolumeFractionVector);
}

/**
 *  This function is called by the external code that
 *  passes a scalar dissipation rate with time history
 */
void Interface::flamelet(const std::vector<double>& sdr, const std::vector<double>& intTime, bool continuation){

    if(sdr.size() != intTime.size())
        throw CamError("Mismatch in the size of SDR and TIME vector\n");

    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);

    //Set the time history of the scalar dissipation rate
    flmlt->setRestartTime(intTime[0]);
    //flmlt->setExternalScalarDissipationRate(intTime,sdr,true);
    flmlt->setExternalTimeSDR(intTime,sdr);

    // Build up a vector of zero soot volume fractions
    std::vector<double> zeroSoot(cg.getnCells(), 0.0);
    // \todo Check if can get rid of this.
    flmlt->setExternalSootVolumeFraction(zeroSoot);

    try{
        const size_t len = sdr.size();
        flamelet(sdr[len-1], intTime[len-1],continuation);
    }catch(CamError& ce){
        throw;
    }

}

/*!
 * This function is called externally to solve a flamelet for a given
 * strain rate.
 *
 *\param[in]    strainRate           Value of strain rate to solve for.
 *\param[in]    lnone				 If 'true', Le=1.
 *
 */
void Interface::flameletStrainRate(const double& strainRate) {

    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);

    flmlt->setExternalStrainRate(strainRate);

    try{
        flmlt->solve(false);
        getFlameletVariables(flmlt);
    }catch(CamError &ce){
        throw ;
    }

}

/*!
 * This function is called externally to solve a flamelet for a given SDR.
 *
 *\param[in]    SDR                 Value of SDR to solve for.
 *\param[in]    lnone                If 'true', Le=1.
 *
 */
void Interface::flameletSDR(const double& SDR) {

    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);

    flmlt->setExternalSDR(SDR);

    try{
        // Solve flamelet at base of flame with soot residual set to zero.
    	flmlt->solve(false,true);
    	//flmlt->solve(false,false);
    	//flmlt->solve(false);
        getFlameletVariables(flmlt);
    }catch(CamError &ce){
        throw ;
    }
}

/**
 * This is called when an sdr profile is required instead of just a constant one.
 */
void Interface::flameletSDRprofile(const std::vector< std::vector<double> >& sdr,
								   const std::vector< std::vector<double> >& Zcoords,
								   const std::vector<double>& intTime,
								   bool continuation) {

    if(sdr.size() != intTime.size())
        throw CamError("Mismatch in the size of SDR and TIME vector\n");

    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);

    flmlt->setRestartTime(intTime[0]);
    if(intTime[1]!=0)cc.setMaxTime(intTime[1]);

    //Set the time history of the scalar dissipation rate
    flmlt->setRestartTime(intTime[0]);

    // Set the scalar dissipation rate profile.
    //flmlt->setExternalScalarDissipationRate(intTime,sdr,Zcoords);


    // Build up a vector of zero soot volume fractions
    std::vector<double> zeroSoot(cg.getnCells(), 0.0);
    flmlt->setExternalSootVolumeFraction(zeroSoot);

    try{
        if (!continuation)
        {
            flmlt->solve(true);
        }
        else
        {
            flmlt->restart(intTime[0]);
        }
        getFlameletVariables(flmlt);
    }catch(CamError &ce){
        throw ;
    }

}

/**
 *  This function is called by the external code that
 *  passes a scalar dissipation rate with time history
 *  and a spatial profile of soot volume fraction.
 */
void Interface::flameletWithSoot(const std::vector<double>& soot_fv, const std::vector<double>& sdr,
                                 const std::vector<double>& intTime, bool continuation){

    if(sdr.size() != intTime.size())
        throw CamError("Mismatch in the size of SDR and TIME vector\n");

    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);
    //Set the time history of the scalar dissipation rate
    flmlt->setRestartTime(intTime[0]);
    //flmlt->setExternalScalarDissipationRate(intTime,sdr,true);

    //@todo check length of the vector
    flmlt->setExternalSootVolumeFraction(soot_fv);

    try{
        const size_t len = sdr.size();
        flamelet(sdr[len-1], intTime[len-1],continuation);
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

void Interface::flamelet(double sdr, double intTime, bool continuation){

    if(intTime!=0)cc.setMaxTime(intTime);
    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);

    if(sdr==0){
        throw std::logic_error("You passed a Scalar Dissipation Rate of zero! Wrong!");
    }
    try{
        if (!continuation)
        {
            // Solve flamelet at base of flame with soot residual set to zero.
        	flmlt->solve(false,true);
        }
        else
        {
            flmlt->restart(intTime);
        }
        getFlameletVariables(flmlt);
    }catch(CamError &ce){
        throw ;
    }
}


/*
 *return the stoichiometric mixture fraction
 */
double Interface::getStMixtureFrac()
{
    if (flmlt == NULL)
    {
        FlameLet* flmlt_temp = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);
        stMixtureFrac = flmlt_temp->stoichiometricMixtureFraction();
        delete flmlt_temp;
    }
    else
    {
        stMixtureFrac = flmlt->stoichiometricMixtureFraction();
    }
    return stMixtureFrac;
}

/*
 *
 *return the density
 */
double Interface::getDensity(const double axpos){

    double dens = getVariableAt(axpos,rhoVector);
    return dens;
}

/*!
 *@param[in]    spIndex     Index of species for which mass fractions requested
 *
 *@return   Vector of mass fractions at all the independent variable points
 */
std::vector<double> Interface::getMassFracsBySpecies(const int spIndex) const {
    // Find the length of the mass fraction profile and create an empty
    // vector with this much space
    const size_t len = indVar.size();
    std::vector<double> mf(len);

    for(size_t i=0; i<len; i++){
        mf[i] = spMassFracs(i,spIndex);
    }

    return mf;
}

std::vector<double> Interface::getMomentsByIndex(const int momentIndex) const {
    // Find the length of the moment profile and create an empty
    // vector with this much space
    const size_t len = indVar.size();
    std::vector<double> momentProfile(len);

    for(size_t i=0; i<len; i++){
        momentProfile[i] = sootMoments(i,momentIndex);
    }
    return momentProfile;
}

// ank25 added for ELFM
std::vector<double> Interface::getMomentsWdotByIndex(const int momentWdotIndex) const {
    // Find the length of the moment profile and create an empty
    // vector with this much space
    const size_t len = indVar.size();
    std::vector<double> momentWdotProfile(len);

    for(size_t i=0; i<len; i++){
        momentWdotProfile[i] = sootMomentsWdot(i,momentWdotIndex);
    }
    return momentWdotProfile;
}



/*!
 *@param[in]    indVarIndex     Mass fractions requested at indVar[indVarIndex]
 *
 *@return   Vector of all mass fractions at one independent variable point
 */
std::vector<double> Interface::getMassFracsByPoint(const int indVarIndex) const {
    std::vector<double> mf(nSpecies);

    for(int i = 0; i < nSpecies; i++){
        mf[i] = spMassFracs(indVarIndex, i);
    }

    return mf;
}
/*
 *return the mass fractions
 */
double Interface::getMassFrac(const int spIndex, const double axpos){
    double massfrac = getVariableAt(axpos, getMassFracsBySpecies(spIndex));
    return massfrac;
}

/*
 *return the moments
 */
double Interface::getMoment(const int momIndex, const double axpos){
    double moment = getVariableAt(axpos, getMomentsByIndex(momIndex));
    return moment;
}


/*
 *return the moments rates
 *ank25 added for ELFM
 */
double Interface::getMomentWdot(const int momIndex, const double axpos){
    double momentWdot = getVariableAt(axpos, getMomentsWdotByIndex(momIndex));
    return momentWdot;
}



/*
 *return the mole fractions
 */
double Interface::getMoleFrac(const int spIndex, const double axpos){

    const double speciesMolwt = (*speciesPointerVector)[spIndex] -> MolWt();
    const double averageMolwt = getVariableAt(axpos,avgMolWtVector);

    double molefrac = getVariableAt(axpos, getMassFracsBySpecies(spIndex))*(averageMolwt/speciesMolwt);

    return molefrac;
}
/*
 *return the temperature
 */
double Interface::getTemperature(const double axpos){

    double temp = getVariableAt(axpos,TVector);
    return temp;
}
/*
 *return the sootAverageDiameter
 */
double Interface::getSootAverageDiameter(const double axpos){

    double temp = getVariableAt(axpos,sootAverageDiameterVector);
    return temp;
}
/*
 *return the sootDispersion
 */
double Interface::getSootDispersion(const double axpos){

    double temp = getVariableAt(axpos,sootDispersionVector);
    return temp;
}
/*
 *return the sootSurfaceArea
 */
double Interface::getSootSurfaceArea(const double axpos){

    double temp = getVariableAt(axpos,sootSurfaceAreaVector);
    return temp;
}
/*
 *return the sootVolumeFraction
 */
double Interface::getSootVolumeFraction(const double axpos){

    double temp = getVariableAt(axpos,sootVolumeFractionVector);
    return temp;
}
/*
 *return the viscosity
 */
double Interface::getViscosity(const double axpos){

    double temp = getVariableAt(axpos,muVector);
    return temp;
}

/**
 *  Return the specific heat
 */
double Interface::getSpecificHeat(const double axPos){
    double temp = getVariableAt(axPos,spHeat);
    return temp;
}

/**
 *  Return the thermal conductivity
 */
double Interface::getThermalConductivity(const double axPos){
    double temp = getVariableAt(axPos,lambda);
    return temp;
}
/**
 *  Return a vector of diffusion coefficients
 */
std::vector<double> Interface::getDiffusionCoefficients(const double axPos){
    std::vector<double> diff, rmDiff;
    int len = indVar.size();
    rmDiff.clear();
    for(int k=0; k<nSpecies; k++){
        diff.clear();
        //get the diffusion coefficient of species k for all mixture fraction
        for(int i=0; i<len; i++){
            diff.push_back(mDiff(i,k));
        }
        double spDiff = getVariableAt(axPos,diff);
        rmDiff.push_back(spDiff);
    }

    return rmDiff;

}

/*!
 * Gets the rate of formation of pyrene. Not sure if this works properly!
 *
 *\param[in]    axpos           Value of independent variable.
 *
 *\return       Rate of formation of pyrene (\f$\mathrm{mol\,m^{-3}\,s^{-1}}\f$).
 */
double Interface::getWdotA4(const double axpos){

    double temp = getVariableAt(axpos,wdotA4);
    return temp;
}

/*!
 * Gets the value of var at a given pos by calling linear_interpolator.
 *
 *\param[in]    pos           Value of independent variable.
 *\param[in]    var           Variable to be interpolated.
 *
 *\return       Value of var at a given pos.
 */
double Interface::getVariableAt
(
    const double& pos,
    const std::vector<double>& var
)
const
{

    Utils::LinearInterpolator<double, double> interpolator
    (
        indVar,
        var
    );

    return interpolator.interpolate(pos);

}


CamAdmin& Interface::getCamAdmin() {
    return ca;
}

CamBoundary& Interface::getCamBoundary() {
    return cb;
}

CamControl&  Interface::getCamControl() {
    return cc;
}

CamGeometry& Interface::getCamGeometry() {
    return cg;
}

CamProfile& Interface::getCamProfile() {
    return cp;
}

CamConfiguration& Interface::getCamConfiguration() {
    return config;
}

