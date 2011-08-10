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
model(NULL),
speciesPointerVector(NULL),
flmlt(NULL),
cp(cg)
{

    std::string fChem("chem.inp");
    std::string fThermo("therm.dat");
    std::string fTrans("tran.dat");
    std::string fCamFlow("camflow.xml");

    try{
        cm.readInput(fCamFlow,cc,cg,convert,ca,cb,cp,config,cSoot);
    }catch(CamError &ce){
        std::cout << ce.errorMessage;
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

    /*
     * populate the moment names
     */
    std::stringstream tempMomentName;
    for(int l=0; l<nMoments; l++)
    {
    	tempMomentName << "M" << l;
    	momentNames.push_back(tempMomentName.str());
    }

}

/*
 *Interface for external drivers which passes a reference to the
 *mechanism object
 */
Interface::Interface(Mechanism& mech_in,
        std::vector<doublereal>& dz,
        std::vector<Thermo::Mixture>& cstrs,
        void* rModel, const doublereal sdr
)
:
model(NULL),
speciesPointerVector(NULL),
flmlt(NULL),
cp(cg)
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

/*
 *return the argument vector with the species names
 */
std::vector<std::string> Interface::getMomentNames(){
    return momentNames;
}


/*!
 * Stores the results for lookup by the CFD program.
 */
void Interface::getFlameletVariables(FlameLet* const flmlt)
{

    flmlt->getDensityVector(rhoVector);
    flmlt->getSpeciesMassFracs(spMassFracs);
    // ank25:  Same for moments.
    flmlt->getTemperatureVector(TVector);
    flmlt->getIndepedantVar(indVar);
    flmlt->getViscosityVector(muVector);
    flmlt->getSpecificHeat(spHeat);
    flmlt->getThermalConductivity(lambda);
    flmlt->getDiffusionCoefficient(mDiff);
    flmlt->getVelocity(mVelocity);
    flmlt->getAverageMolarWeight(avgMolWtVector);
    flmlt->getWdotA4(wdotA4);

}

/**
 *  This function is called by the external code that
 *  passes a scalar dissipation rate with time history
 */
void Interface::flamelet(const std::vector<doublereal>& sdr, const std::vector<doublereal>& intTime, bool continuation, bool lnone){

    if(sdr.size() != intTime.size())
        throw CamError("Mismatch in the size of SDR and TIME vector\n");

    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);

    //Set the time history of the scalar dissipation rate
    flmlt->setRestartTime(intTime[0]);
    //flmlt->setExternalScalarDissipationRate(intTime,sdr,true);
    flmlt->setExternalTimeSDR(intTime,sdr);

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

/*!
 * This function is called externally to solve a flamelet for a given
 * strain rate.
 *
 *\param[in]    strainRate           Value of strain rate to solve for.
 *\param[in]    lnone				 If 'true', Le=1.
 *
 */
void Interface::flameletStrainRate(const doublereal& strainRate, bool lnone) {

    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);
    if(!lnone) flmlt->setLewisNumber(FlameLet::LNNONE);

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
void Interface::flameletSDR(const doublereal& SDR, bool lnone) {

    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);
    if(!lnone) flmlt->setLewisNumber(FlameLet::LNNONE);

    flmlt->setExternalSDR(SDR);

    try{
        // Solve flamelet at base of flame with soot residual set to zero.
    	flmlt->solve(false,true);
    	//flmlt->solve(false);
        getFlameletVariables(flmlt);
    }catch(CamError &ce){
        throw ;
    }
}

/**
 * This is called when an sdr profile is required instead of just a constant one.
 */
void Interface::flameletSDRprofile(const std::vector< std::vector<doublereal> >& sdr,
								   const std::vector< std::vector<doublereal> >& Zcoords,
								   const std::vector<doublereal>& intTime,
								   bool continuation,
								   bool lnone) {

    if(sdr.size() != intTime.size())
        throw CamError("Mismatch in the size of SDR and TIME vector\n");

    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);

    flmlt->setRestartTime(intTime[0]);
    if(intTime[1]!=0)cc.setMaxTime(intTime[1]);
    if(!lnone) flmlt->setLewisNumber(FlameLet::LNNONE);

    //Set the time history of the scalar dissipation rate
    flmlt->setRestartTime(intTime[0]);

    // Set the scalar dissipation rate profile.
    //flmlt->setExternalScalarDissipationRate(intTime,sdr,Zcoords);


    // Build up a vector of zero soot volume fractions
    std::vector<doublereal> zeroSoot(cg.getnCells(), 0.0);
    flmlt->setExternalSootVolumeFraction(zeroSoot);

    try{
        if (!continuation)
        {
            flmlt->solve(true);
        }
        else
        {
            flmlt->restart();
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
void Interface::flameletWithSoot(const std::vector<doublereal>& soot_fv, const std::vector<doublereal>& sdr,
                                 const std::vector<doublereal>& intTime, bool continuation, bool lnone){

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
    if(flmlt == NULL ) flmlt = new FlameLet(ca, config, cc, cg, cp, cSoot, mech);
    if(!lnone) flmlt->setLewisNumber(FlameLet::LNNONE);

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
            flmlt->restart();
        }
        getFlameletVariables(flmlt);
    }catch(CamError &ce){
        throw ;
    }
}
/*
 *return the stoichiometric mixture fraction
 */
doublereal Interface::getStMixtureFrac()
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
doublereal Interface::getDensity(const doublereal axpos){

    doublereal dens = getVariableAt(axpos,rhoVector);
    return dens;
}

/*!
 *@param[in]    spIndex     Index of species for which mass fractions requested
 *
 *@return   Vector of mass fractions at all the independent variable points
 */
std::vector<doublereal> Interface::getMassFracsBySpecies(const int spIndex) const {
    // Find the length of the mass fraction profile and create an empty
    // vector with this much space
    const size_t len = indVar.size();
    std::vector<doublereal> mf(len);

    for(size_t i=0; i<len; i++){
        mf[i] = spMassFracs(i,spIndex);
    }

    return mf;
}

std::vector<doublereal> Interface::getMomentsByIndex(const int momentIndex) const {
    // Find the length of the moment profile and create an empty
    // vector with this much space
    const size_t len = indVar.size();
    std::vector<doublereal> momentProfile(len);

    for(size_t i=0; i<len; i++){
        momentProfile[i] = sootMoments(i,momentIndex);
    }
    return momentProfile;
}


/*!
 *@param[in]    indVarIndex     Mass fractions requested at indVar[indVarIndex]
 *
 *@return   Vector of all mass fractions at one independent variable point
 */
std::vector<doublereal> Interface::getMassFracsByPoint(const int indVarIndex) const {
    std::vector<doublereal> mf(nSpecies);

    for(int i = 0; i < nSpecies; i++){
        mf[i] = spMassFracs(indVarIndex, i);
    }

    return mf;
}
/*
 *return the mass fractions
 */
doublereal Interface::getMassFrac(const int spIndex, const doublereal axpos){
    doublereal massfrac = getVariableAt(axpos, getMassFracsBySpecies(spIndex));
    return massfrac;
}

/*
 *return the moments
 */
doublereal Interface::getMoment(const int momIndex, const doublereal axpos){
    doublereal moment = getVariableAt(axpos, getMomentsByIndex(momIndex));
    return moment;
}


/*
 *return the mole fractions
 */
doublereal Interface::getMoleFrac(const int spIndex, const doublereal axpos){

    const doublereal speciesMolwt = (*speciesPointerVector)[spIndex] -> MolWt();
    const doublereal averageMolwt = getVariableAt(axpos,avgMolWtVector);

    doublereal molefrac = getVariableAt(axpos, getMassFracsBySpecies(spIndex))*(averageMolwt/speciesMolwt);

    return molefrac;
}
/*
 *return the temperature
 */
doublereal Interface::getTemperature(const doublereal axpos){

    doublereal temp = getVariableAt(axpos,TVector);
    return temp;
}
/*
 *return the viscosity
 */
doublereal Interface::getViscosity(const doublereal axpos){

    doublereal temp = getVariableAt(axpos,muVector);
    return temp;
}

/**
 *  Return the specific heat
 */
doublereal Interface::getSpecificHeat(const doublereal axPos){
    doublereal temp = getVariableAt(axPos,spHeat);
    return temp;
}

/**
 *  Return the thermal conductivity
 */
doublereal Interface::getThermalConductivity(const doublereal axPos){
    doublereal temp = getVariableAt(axPos,lambda);
    return temp;
}
/**
 *  Return a vector of diffusion coefficients
 */
std::vector<doublereal> Interface::getDiffusionCoefficients(const doublereal axPos){
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

/*!
 * Gets the rate of formation of pyrene. Not sure if this works properly!
 *
 *\param[in]    axpos           Value of independent variable.
 *
 *\return       Rate of formation of pyrene (\f$\mathrm{mol\,m^{-3}\,s^{-1}}\f$).
 */
doublereal Interface::getWdotA4(const doublereal axpos){

    doublereal temp = getVariableAt(axpos,wdotA4);
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
doublereal Interface::getVariableAt
(
    const doublereal& pos,
    const std::vector<doublereal>& var
)
const
{

    Utils::LinearInterpolator<doublereal, doublereal> interpolator
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

