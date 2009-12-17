
#include "array.h"

#include <cmath>
#include "cam_setup.h"
#include "cam_residual.h"
#include <vector>
#include "cam_reporter.h"
#include <iostream>
#include <stdexcept>
#include "cam_profile.h"
#include "cam_admin.h"
#include "flamelet.h"
#include "cam_math.h"
#include "cvode_wrapper.h"


using namespace Camflow;
using namespace std;

/*
 *this is called by the model object. The boolean interface decides
 *if the call originates from the interface or from camflow kernel
 */
void FlameLet::solve(CamControl& cc, CamAdmin& ca, CamGeometry& cg, CamProfile& cp,
        CamConfiguration& config, CamSoot& cs, Mechanism& mech){

    solve(cc,ca,cg,cp,mech,false);
    
}

void FlameLet::solve(CamControl& cc, CamAdmin& ca, CamGeometry& cg,
                                CamProfile& cp, Mechanism& mech, bool interface){


    /*
     *check for mixture fraction bounds. In this case the input given in the
     * grid file is treated as mixture fraction coordinates. The reactor
     * length specified in the reactor section is simply ignored.
     */
   
    CamBoundary cb;
    Thermo::Mixture mix(mech.Species());
    camMixture = &mix;

    storeObjects(cc,ca,cg,cp,cb,mech);
   
    nSpc = camMech->SpeciesCount(); // number of species
    ptrT = nSpc;                    // temperature offset
    nVar = nSpc+1;                  // number of variables

    //if(reacGeom->getnCells()==0) reacGeom->discretize();

    reacGeom->addZeroWidthCells();

    mCord = reacGeom->getnCells();
   
    nEqn = nVar*mCord;              // number of equations
    iMesh_s = 1;                    // internal grid starting
    iMesh_e = mCord-1;              // internal grid ending

    cellBegin = 0;
    cellEnd = mCord;
    
    /*
     *get the strain rate
     */
    strain = ca.getStrainRate();
    //profile->setGeometryObj(cg);
    reporter = new CamReporter();
    

    /*
     *init the solution vector
     */
    initSolutionVector(cc);
    reporter->header("Flamelet");
    header();
    if(cc.getSolutionMode() == cc.COUPLED)
        csolve(cc,interface);
    else{
        ssolve(cc);
        csolve(cc,interface);
    }


}

/**
 */
void FlameLet::setRestartTime(doublereal t){
    rstartTime = t;
}

/**
 *contuation call from an external code that
 *solves for population balance
 */
void FlameLet::solve(vector<Thermo::Mixture>& cstrs,
        const vector<vector<doublereal> >& iniSource,
        const vector<vector<doublereal> >& fnlSource,
        Mechanism& mech,
        CamControl& cc,
        CamAdmin& ca,
        CamGeometry& cg,
        CamProfile& cp){

    CamBoundary cb;
    Thermo::Mixture mix(mech.Species());
    camMixture = &mix;

    storeObjects(cc,ca,cg,cp,cb,mech);

    nSpc = camMech->SpeciesCount(); // number of species
    ptrT = nSpc;                    // temperature offset
    nVar = nSpc+1;

    reacGeom->addZeroWidthCells();

    mCord = reacGeom->getnCells();

    nEqn = nVar*mCord;              // number of equations
    iMesh_s = 1;                    // internal grid starting
    iMesh_e = mCord-1;              // internal grid ending

    cellBegin = 0;
    cellEnd = mCord;
    
     /*
     *set the source terms for the particle process
     */
    setParticleSource(iniSource,fnlSource);
    
    /*
     *  reset the solution vector. cstrs contain only
     *  the interor cells. The inlet condisions need to
     *  be taken care of.
     */
    CamBoundary left, right;
    admin->getLeftBoundary(left);
    admin->getRightBoundary(right);
    storeInlet(left,fuel);    
    storeInlet(right,oxid);
    fuel.T = left.getTemperature();
    oxid.T = right.getTemperature();

    solvect.clear();
    solvect.resize(nEqn,0.0);

    /**!
     *  Inlet  boundary ie. z=0
     *  this is the oxidizer
     */
    //Species
    for(int l=0; l<nSpc; l++){
        solvect[l] = oxid.Species[l];
    }
    solvect[ptrT] = oxid.T;

    /**!
     *  Interior mesh points
     */
    for(int i=iMesh_s; i<iMesh_e;i++){
        vector<doublereal> massFrac;
        cstrs[i-1].GetMassFractions(massFrac);
        for(int l=0; l<nSpc; l++){
            solvect[i*nVar+l] = massFrac[l];
        }
        //Temperature
        solvect[i*nVar+ptrT] = cstrs[i-1].Temperature();
    }

    /**
     *  outlet boundary
     *  this is the fuel inlet
     */
    //Species
    for(int l=0; l<nSpc; l++){
        solvect[iMesh_e*nVar+l] = fuel.Species[l];
    }
    //Temperature
    solvect[iMesh_e*nVar+ptrT] = fuel.T;


}


void FlameLet::initSolutionVector(CamControl &cc){


    /*
     *initialize the geometry
     */
    dz = reacGeom->getGeometry();

    /*
     *left boundary is for the fuel and right boundary is
     *for oxidizer
     */
    CamBoundary left, right;
    admin->getLeftBoundary(left);
    admin->getRightBoundary(right);
    storeInlet(left,fuel);    
    storeInlet(right,oxid);
    fuel.T = left.getTemperature();
    oxid.T = right.getTemperature();

    stoichZ = stoichiometricMixtureFraction();
    profile->setMixingCenter(stoichZ);
    profile->setMixingWidth(0.5*stoichZ);
    /*
     *initialize the ODE vector
     */
    solvect.resize(nEqn,0);
    vector<doublereal> vSpec, vT;
    /*
     * actual signature follows (left,right,cc,vSpec)
     * but in the case of flamelets the species mass fractions
     * are solved for the mixture fraction coordinate, whose
     * direction is taken the same as physical space.
     * z=0 corresponds to oxidizer and z=1 corresponds to fuel
     * therefore the inlets are interchanged here to initialize
     * the species vector properly
     */
    initSpecies(right,left,cc,vSpec);
    /*
     *the following will initialize the temperature vecotor with
     *a linear profile
     */
    //initTempGauss(vT);
    doublereal inrsctOx, inrsctFl;
    doublereal slopeOx, slopeFl;
    inrsctOx = oxid.T;
    
    slopeOx = (2000.0-oxid.T)/stoichZ;
    slopeFl = (2000.0-fuel.T)/(stoichZ-1.0);
    inrsctFl = fuel.T - slopeFl;
    vector<doublereal> position = reacGeom->getAxpos();
    int len = position.size();
    vT.resize(len,0.0);
    for(unsigned int i=0; i<dz.size();i++){
        if(position[i] < stoichZ){
            vT[i] = slopeOx*position[i] + inrsctOx;
        }else{
            vT[i] = slopeFl*position[i] + inrsctFl;
        }
    }

    /*
     *fix the temperature for mixture fraction zero
     *i.e oxidizer inlet
     */
    vT[0] =oxid.T;
   
    /*
     *fix the temperature for mixture fraction 1.
     *i.e the fuel inlet
     */
    vT[iMesh_e] = fuel.T;
            
    /*
     *create the actual solution vector by merging the species
     *vector and the temperature vector
     */
    mergeSpeciesVector(&vSpec[0]);
    mergeEnergyVector(&vT[0]);
    

    stoichiometricMixtureFraction();
}

/*
 *coupled solver
 */
void FlameLet::csolve(CamControl& cc, bool interface){
    int solverID = cc.getSolver();
    if ( solverID == cc.CVODE){
        CVodeWrapper cvw;
        eqn_slvd = EQN_ALL;
        int band = nVar*2;
        
        cvw.init(nEqn,solvect,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(),
                            cc.getMaxTime(),band,*this);

        cvw.solve(CV_ONE_STEP,cc.getResTol());
        /*
         *write the output to file only if the call is not
         *from the interface
         */
        if(!interface) {
            reportToFile(cc.getMaxTime(),&solvect[0]);            
        }
        cvw.destroy();
        
    }else if( solverID == cc.NEWTON){
        cout << "Not implemented\n";
    }
}



/*
 *restart the solution. This is normally called from the interface routine
 *The solver is reinitialized each time with the previous solution.
 */
void FlameLet::restart(CamControl& cc){

    Thermo::Mixture mix(camMech->Species());
    camMixture = &mix;
    
    CVodeWrapper cvw;
    eqn_slvd = EQN_ALL;
    int band = nVar*2;
    
    cvw.init(nEqn,solvect,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(),
                            cc.getMaxTime(),band,*this,rstartTime);
    cvw.solve(CV_ONE_STEP,cc.getResTol());
}
/*
 *segregated solver
 */
void FlameLet::ssolve(CamControl& cc){

    CVodeWrapper cvw;
    int seg_eqn, band;
    vector<doublereal> seg_soln_vec;

    for(int i=0; i<cc.getNumIterations();i++){
        /*
         *solve species equation
         */
        cout << "Solving species equations  " << i << endl;
        int dd; cin >> dd;
        eqn_slvd = EQN_SPECIES;
        seg_eqn = nSpc*mCord;        
        band = nSpc*2;
        extractSpeciesVector(seg_soln_vec);
        cvw.init(seg_eqn,seg_soln_vec,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(),
                            cc.getMaxTime(),band,*this);
        cvw.solve(CV_ONE_STEP,1e-03);
        mergeSpeciesVector(&seg_soln_vec[0]);
        cvw.destroy();

        /*
         *solve energy equation
         */
        cout << "Solving energy equation  " << i << endl;
        cin >> dd;
        eqn_slvd = EQN_ENERGY;
        seg_eqn = mCord;
        band = 1;
        extractEnergyVector(seg_soln_vec);
        cvw.init(seg_eqn,seg_soln_vec,cc.getSpeciesAbsTol(),cc.getSpeciesRelTol(),
                                    cc.getMaxTime(),band,*this);
        cvw.solve(CV_ONE_STEP,1e-03);
        mergeEnergyVector(&seg_soln_vec[0]);
        cvw.destroy();
    }


    
}
/*
 *residual definitions
 */
void FlameLet::residual(const doublereal& t, doublereal* y, doublereal* f){

    resSp.resize(nSpc*mCord,0.0);
    resT.resize(mCord,0.0);

    if(eqn_slvd == EQN_ALL){                         // coupled solver part
        saveMixtureProp(y);        
        speciesResidual(t,y,&resSp[0]);
        energyResidual(t,y,&resT[0]);

        for(int i=0;i<mCord;i++){
            for(int l=0;l<nSpc;l++){
                f[i*nVar+l] = resSp[i*nSpc+l];
            }
            f[i*nVar+ptrT] = resT[i];
        }
        
    }else{                                          // segregated solver part
        if(eqn_slvd == EQN_SPECIES){
            mergeSpeciesVector(y);
            saveMixtureProp(&solvect[0]);
            speciesResidual(t,y,f);
        }else if(eqn_slvd==EQN_ENERGY){
            mergeEnergyVector(y);
            saveMixtureProp(&solvect[0]);
            energyResidual(t,y,f);
        }
    }
}

/*
 *species residual definitions
 */
void FlameLet::speciesResidual(const doublereal& t, doublereal* y,
                                                        doublereal* f){

    int i; //mixture fraction coordinate index
    doublereal grad_e, grad_w;
    doublereal zPE, zPW;
    doublereal source;

    /*
     *starting with mixture fraction zero: i.e oxidizer
     *inlet. The fuel composition is zero. Left inlet is
     *considered as the oxidizer inlet and the concentrationd
     *are help constant
     */


    /*
     *set the externally specified scalar dissipation rate
     *to sdr
     */
    if(timeHistory){
        sdr = getSDR(t);
       
    }else  if(sdr_ext!=0){
        sdr=sdr_ext;
    }
    i=0;
    //sdr = scalarDissipationRate(dz[i]);
    for(int l=0; l<nSpc; l++){
        f[l] = 0.0;
        
    }

    /**
     *  For non-unity Lewis numbers
     */
    
    Le.resize(mCord,nSpc,1.0);
    convection.resize(mCord,nSpc,0);
    doublereal onebLe;
    doublereal oneby16 = 1.0/16;
    if(Lewis == FlameLet::LNNONE){       
        
        //Iterating over the interior points
        for(int i=1; i<mCord-1; i++){
            for(int l=0; l<nSpc; l++){
                Le(i,l) = m_k[i]/(m_rho[i]*m_cp[i]*s_Diff(i,l));
                onebLe = 1/Le(i,l);
                        
                convection(i,l) = oneby16*(onebLe-1)*(s_mf(i+1,l)-s_mf(i-1,l))*(m_rho[i+1]-m_rho[i-1]);
            }
    
        }
        for(int l=0; l<nSpc; l++){
            // Computation at the i = 0 end of the grid
            Le(0,l) = m_k[0]/(m_rho[0]*m_cp[0]*s_Diff(0,l));
            onebLe = 1/Le(0,l);
            convection(0,l) = oneby16*(onebLe-1)*4*(s_mf(1,l)-s_mf(0,l))*(m_rho[1]-m_rho[0]);


            // Computation at the i = mCord - 1 end of the grid
            Le(mCord-1,l) = m_k[mCord-1]/(m_rho[mCord -1]*m_cp[mCord-1]*s_Diff(mCord-1,l));
            onebLe = 1/Le(mCord-1,l);
            convection(mCord-1,l) = oneby16*(onebLe-1)*4*(s_mf(mCord-1,l)-s_mf(mCord-2,l))*(m_rho[mCord-1]-m_rho[mCord-2]);
        }


    }
    /*
     *interior mixture fraction coordinates
     */
    
    for(i=iMesh_s; i<iMesh_e;i++){

        zPE = 0.5*(dz[i]+dz[i+1]);
        zPW = 0.5*(dz[i]+dz[i-1]);
        if(sdr_ext==0)sdr = scalarDissipationRate(dz[i]);
        for(int l=0; l<nSpc; l++){
            grad_e = (s_mf(i+1,l)-s_mf(i,l))/zPE;
            grad_w = (s_mf(i,l)-s_mf(i-1,l))/zPW;
            source = s_Wdot(i,l)*(*spv)[l]->MolWt()/m_rho[i];
            f[i*nSpc+l] = sdr*(grad_e-grad_w)/(2*Le(i,l)*dz[i]) + source + convection(i,l)*sdr/(m_rho[i]*dz[i]*dz[i]);
        }

    }
    /*
     *Mixture fraction 1. The fuel inlet. Concentrations are
     *held constant at the fuel composition
     */
    
    for(int l=0; l<nSpc; l++){
        f[iMesh_e*nSpc+l] = 0.0;
    }
    

}
/* This code computes the spectral and soot related components or radiative heat loss.  

    References:
    1.	Radiation Models, International Workshop on Measurement and Computation of Turbulent  Nonpremixed Flames,
        www.sandia.gov/TNF/radiation.html.  This reference covers the details of implementing the spectral part of 
        the radiation term.

    2.	Kim, S-K., Kim, Y., Assessment of the Eulerian particle flamelet model for nonpremixed turbulent jet flames, 
        Combustion and Flame 154 (2008) 232-247.  This article shows a modern implementation of the procedure used 
        in [1], for the spectral part of the radiation term.

    3.	Carbonnel, D., Oliva, A., Perez-Segarra, C.D., Implementation of two-equation soot flamelet models for laminar
        diffusion flames.  This paper includes a term, used here, for modelling soot radiative heat loss.

    4.  Grosshandler, W.L., RADCAL: A Narrow-Band Model For Radiation Calculations in a Combustion Environment, NIST 
        technical note 1402, 1993.
   
   */



 
/*!
    *Computes the Planck mean absorption coefficients
 
    *@param[in]      Temperature      Temperature at a point of the grid, as obtained by the energy residual function.
    *@param[out]     Absorption       An output array containing the computed Planck absorption coefficients.      
    *
    * The coefficients are produced for H2O, CO2, and CO, because these are the species that tend to produce the 
    * greatest spectral radiative effect in mixtures commonly associated with combustion. Coefficients are known
    * for CH4 as well, but they can lead to inaccurate results, and therefore are not included here.  The computation
    * of these coefficients are based on data from the RADCAL model.  This data was for temperatures between 300K 
    * and 2500K.
    
    */    
void FlameLet::PlanckAbsorption (const doublereal temperature, doublereal Absorption[3])const
{
        
    
    //This quantity is reused repeatedly in the equations below
    const doublereal beta = 1000/temperature;

    //This helps avoid repeated, costly calls to pow().  There is a modest tradeoff of speed vs. accuracy here.
    const doublereal beta2 = beta * beta;
    const doublereal beta3 = beta2 * beta;
    const doublereal beta4 = beta3 * beta;
    const doublereal beta5 = beta4 * beta;
 
    // This code computes the Planck mean absorption coefficient for H2O, using an equation in reference [1].  The specific 
    // equation is based on fitted results from the RADCAL program. See reference [4].  Multiplication by AtmToPascal 
    // converts the absorption coefficients from units of 1/m * 1/atm  to 1/m * 1/pa.

    const doublereal AtmToPascal = 101325.0;

    Absorption[0] = (-0.23093 - 1.12390 * beta + 9.41530 * beta2 - 2.99880 * beta3 + 0.51382 * beta4
                    - 1.86840e-5 * beta5)*AtmToPascal;
 

    // This code computes the Planck mean absorption coefficient for CO2, using an equation in reference [1].  The specific 
    // equation is based on fitted results from the RADCAL program. See reference [4]. Multiplication by AtmToPascal 
    // converts the absorption coefficients from units of 1/m * 1/atm  to 1/m * 1/pa.


    Absorption[1] = (-18.7410 - 121.3000 * beta + 273.5000 * beta2 - 194.0500 * beta3 + 56.3100 * beta4
                    - 5.8169 * beta5)*AtmToPascal;


    // This code computes the Planck mean absorption coefficient for CO, using two equations in reference [1].  It uses a 
    // conditional statement on temperature to determine which of the two equations to use. The specific  equation is based
    // on fitted results from the RADCAL program. See reference [4]. Multiplication by AtmToPascal converts the absorption 
    // coefficients from units of 1/m * 1/atm  to 1/m * 1/pa.


    if (temperature <= 750) {
        Absorption[2]  =  (4.7869 + temperature * (-0.06953 + temperature * (2.95775e-4 + temperature * 
                          (-4.25732e-7 + temperature * 2.202849e-10))))*AtmToPascal; 
    }
    else {
        Absorption[2]  = (10.0900 + temperature * (-0.01183 + temperature * (4.7753e-6 + temperature * 
                         (-5.87209e-10 + temperature * 2.5334e-14))))*AtmToPascal;
    }
    
    
    // After this function has been called the end result is a vector populated with absorption coefficients, to be used in 
    // the function RadiativeLoss.
}

//
/*! 
      *Computes  the radiative heat loss term for radiative heat dissipation model

      *@param[in]      rho                             Density of the mixture, as obtained by the energy residual function.
      *@param[in]      cp                              Specific heat of the mixture, as obtained by the energy residual function.
      *@param[in]      Temperature                     Temperature at a point of the grid, as obtained by the energy residual function.
      *@param[in]      SootVolFrac                     The soot value fraction, as provided by the user or obtained from an external source.  
      *@return                                         A term equivalent to the total radiative heat loss.
      *
      *                                                This function will be called from the energy residual code in Camflow's Flamelet class,
      *                                                which will provide such parameters as temperature, density, specific heat, mass fractions 
      *                                                and soot volume fractions.  
      *
      *
      */ 


doublereal FlameLet::RadiativeLoss(const doublereal rho, const doublereal cp, const doublereal temperature, 
                                   const doublereal soot_vol_frac, const doublereal mole_frac_H2O, 
                                   const doublereal mole_frac_CO2, const doublereal mole_frac_CO) const {               
    // RadiativeLoss requires a background setting temperature.  It is usually assumed to be 300K, unless experimental conditions 
    // suggest another temperature.
    const doublereal BackgroundTemp = 300;


    // Absorption coefficients are for H2O, CO2, and CO, in that order.
    doublereal AbsorptionVector[3] = {0,0,0};
    

    // The function calls PlanckAbsorption to produce the absorption coefficients.
    PlanckAbsorption(temperature, AbsorptionVector);


    // Using a single equation in reference [1], the function computes a Radiative Heat Dissipation term.  It uses the Planck mean 
    // absorption coefficients for H2O, CO2 abd CO, as well as the density of each of these species, the specific heat of each species,
    // the temperature, the partial pressure of each species, and the soot volume fraction.

             
    // Partial pressures of species k.  0 corresponds to H2O, 1 corresponds to CO2, 2 corresponds to CO. Operating pressure
    // must be expressed in Pascals in the camflow.xml file.
     
    doublereal partialPress[3];
    partialPress[0]  =  mole_frac_H2O * opPre;
    partialPress[1]  =  mole_frac_CO2 * opPre;
    partialPress[2]  =  mole_frac_CO * opPre;

    
    doublereal radiationScalar;  //To hold intermediate results
    doublereal spectralRadiation = 0;   //An intermediate result    
    
    //The following is used repeatedly in the loop below
    const doublereal temperaturePowers = 1/(rho * cp) * 4* 5.669e-8 * (pow(temperature, 4) - pow(BackgroundTemp, 4));


    for (unsigned int j = 0; j < 3; ++j){
        // 0 = H2O,  1 = CO2,  2 = CO


        //Spectral radiative heat loss due to each of H2O, CO2, and CO 
        radiationScalar =  temperaturePowers * partialPress[j] * AbsorptionVector[j];

        //Total spectral radiative heat loss.
        spectralRadiation += radiationScalar;
    }
     
    //Total spectral radiative heat loss + radiative heat loss due to soot
    return spectralRadiation + 3.337e-4 * soot_vol_frac * pow(temperature, 5);
    
    cout << "Spectral Radiation\n" << spectralRadiation;


}
       


/*
 *energy residual
 */
void FlameLet::energyResidual(const doublereal& t, doublereal* y, doublereal* f){
    int i; //mixture fraction coordinate index
    doublereal grad_e, grad_w;
    doublereal zPE, zPW;
    doublereal source;
    
    /*
     *starting with mixture fraction zero: i.e oxidizer
     *inlet. The temperature is fixed at the oxidizer
     *inlet temperature
     */
    i=0;    
    f[i] = 0.0;

    /*
     *set the externally specified scalar dissipation rate
     *to sdr
     */
//    if(timeHistory){
//        sdr = getSDR(t);
//
//    }else  if(sdr_ext!=0){
//        sdr=sdr_ext;
//    }
    
    /*
     *intermediate mixture fraction coordinates
     */
    
    

    for(int i=iMesh_s; i<iMesh_e; i++){
        zPE = 0.5*(dz[i]+dz[i+1]);
        zPW = 0.5*(dz[i]+dz[i-1]);
        source = 0.0;
        
        for(int l=0; l<nSpc; l++){
            source += s_Wdot(i,l)*s_H(i,l);
        }
        
        if(sdr_ext==0)sdr = scalarDissipationRate(dz[i]);
        grad_e = (m_T[i+1]-m_T[i])/zPE;
        grad_w = (m_T[i]-m_T[i-1])/zPW;

        f[i] = sdr*(grad_e-grad_w)/(2*dz[i])-(source/(m_rho[i]*m_cp[i]));

        /**
         *  Accounting for non-unity Lewis number
         */
        if(Lewis == FlameLet::LNNONE){
            doublereal conduction = 0;            
            conduction = sdr*(m_cp[i+1]-m_cp[i-1])*(m_T[i+1]-m_T[i-1])/(8*m_cp[i]*dz[i]*dz[i]);            
            doublereal enthFlux = 0;
            doublereal tGrad = (m_T[i+1]-m_T[i-1])/dz[i];
            doublereal cpterm=0;
            doublereal spGrad =0;
            for(int l=0; l<nSpc; l++){
                cpterm = (1-CpSpec(i,l)/m_cp[i]);
                spGrad = ((s_mf(i+1,l)-s_mf(i-1,l))/dz[i]) + (s_mf(i,l)*(avgMolWt[i+1]-avgMolWt[i-1])/(avgMolWt[i]*dz[i]));
                enthFlux += spGrad*cpterm/Le(i,l);
            }

            
            f[i] -= enthFlux*sdr*tGrad/8.0 ;

        }
    
        //======Radiative Heat Loss Term===============

               
        //Get species indexes corresponding to H20, CO2, CO
        const int iH2O = camMech->FindSpecies("H2O");
        const int iCO2 = camMech->FindSpecies("CO2");
        const int iCO  = camMech->FindSpecies("CO");

     
      

        //The next few lines access the molecular weights of H2O, CO2 and CO.
        const Sprog::SpeciesPtrVector *speciesDetails = camMixture->Species();
        const doublereal molwtH2O =   (*speciesDetails)[iH2O] -> MolWt();
        const doublereal molwtCO2 =   (*speciesDetails)[iCO2] -> MolWt();
        const doublereal molwtCO =    (*speciesDetails)[iCO] -> MolWt();


        //Computation of mole fractions as inputs to RadiativeLoss
        //Computed using the following equation: 
        //mole fraction (species) = mass fraction(species) * molecular mass (average) * (1/molecular mass (species)) 
        const doublereal mole_fracsH2O = s_mf(i,iH2O)*avgMolWt[i]/molwtH2O;
        const doublereal mole_fracsCO2 = s_mf(i,iCO2)*avgMolWt[i]/molwtCO2;                       
        const doublereal mole_fracsCO = s_mf(i,iCO)*avgMolWt[i]/molwtCO;

        //Soot Volume Fraction is set to zero here
        radiation.resize(mCord,0.0);
        
        //This radiation term is sentt to as output to profile.h 
        radiation[i] = RadiativeLoss(m_rho[i], m_cp[i], m_T[i], m_SootFv[i], mole_fracsH2O, mole_fracsCO2, mole_fracsCO);
        
        //This is the new energy residual term, accounting for radiation.
        f[i] += RadiativeLoss(m_rho[i], m_cp[i], m_T[i], m_SootFv[i], mole_fracsH2O, mole_fracsCO2, mole_fracsCO);

    }

    
    f[iMesh_e] = 0.0;
    
    
}
/*
 *save the mixture property
 */
void FlameLet::saveMixtureProp(doublereal* y){
    s_mf.resize(mCord,nSpc);
    s_Wdot.resize(mCord,nSpc);
    s_H.resize(mCord,nSpc);
    s_Diff.resize(mCord,nSpc);
    CpSpec.resize(mCord,nSpc);
    m_T.resize(mCord,0.0);
    m_rho.resize(mCord,0.0);
    m_cp.resize(mCord,0.0);
    m_mu.resize(mCord,0.0);
    m_u.resize(mCord,0.0);
    m_k.resize(mCord,0.0);
    
    avgMolWt.resize(mCord,0.0);
 

    vector<doublereal> mf, htemp, temp, cptemp;
    htemp.resize(nSpc,0.0);
    for(int i=0; i<mCord; i++){
        mf.clear();
        for(int l=0; l<nSpc; l++){
            mf.push_back(y[i*nVar+l]); 
        }
        m_T[i] = y[i*nVar+ptrT];                                   //temperature
        camMixture->SetMassFracs(mf);                              //mass fraction 
        camMixture->SetTemperature(m_T[i]);                        //temperature
       
        avgMolWt[i] = camMixture->getAvgMolWt();
        m_rho[i] = opPre*avgMolWt[i]/(R*m_T[i]);                   //density
        camMixture->SetMassDensity(m_rho[i]);                      //density
        camMech->Reactions().GetMolarProdRates(*camMixture,wdot);
        htemp = camMixture->getMolarEnthalpy();                    //enthalpy
        m_cp[i] = camMixture->getSpecificHeatCapacity();           //specific heat
        m_k[i] = camMixture->getThermalConductivity(opPre);        //thermal conductivity
        m_mu[i] = camMixture->getViscosity();                      //mixture viscosity
        temp = camMixture->getMixtureDiffusionCoeff(opPre);
        cptemp = camMixture->getMolarSpecificHeat();
        for(int l=0; l<nSpc; l++){
            s_mf(i,l) = mf[l];
            s_Wdot(i,l) = wdot[l];
            s_H(i,l) = htemp[l];
            s_Diff(i,l) = temp[l];
            //Specific heat capacity of species in J/Kg K
            CpSpec(i,l) =cptemp[l]/(*spv)[l]->MolWt();
        }
    }
}

doublereal FlameLet::stoichiometricMixtureFraction(){
    /*
     *check for C and H atoms
     */

    vector<Sprog::Species*> fspecies, ospecies;

    int indx_C = camMech->FindElement("C");
    int indx_H = camMech->FindElement("H");

    /*
     *fuel inlet
     */
    CamBoundary fuelInlet, oxInlet;
    map<string, doublereal> species;
    map<string, doublereal>::iterator sIterator;
    admin->getLeftBoundary(fuelInlet);
    species = fuelInlet.getInletSpecies();
    sIterator = species.begin();

    while(sIterator != species.end()){
        fspecies.push_back(camMech->GetSpecies(sIterator->first));
        sIterator++;
    }

    admin->getRightBoundary(oxInlet);
    species = oxInlet.getInletSpecies();
    sIterator = species.begin();
    while(sIterator != species.end()){
        ospecies.push_back(camMech->GetSpecies(sIterator->first));
        sIterator++;
    }

    int cAtoms=0;
    int hAtoms=0;
    doublereal avgMolWt=0;
    doublereal fuelMassFrac=0;
    unsigned int i;
    vector<doublereal> temp;
    getInletMassFrac(fuelInlet,temp);

    int iN2 = camMech->FindSpecies("N2");
    int iAR = camMech->FindSpecies("AR");
    int iHe = camMech->FindSpecies("HE");
    int iO2 = camMech->FindSpecies("O2");

    for(i=0; i<fspecies.size();i++){
        int icAtoms = fspecies[i]->AtomCount(indx_C);
        int ihAtoms = fspecies[i]->AtomCount(indx_H);
        if (icAtoms != 0 || ihAtoms !=0) avgMolWt += fspecies[i]->MolWt();
        cAtoms += icAtoms; hAtoms += ihAtoms;

        int spIndx = camMech->FindSpecies(fspecies[i]->Name());
        if(spIndx != iN2 && spIndx != iAR && spIndx != iHe){
            fuelMassFrac += temp[spIndx];
        }
    }

    getInletMassFrac(oxInlet,temp);
    doublereal o2MassFrac = temp[iO2];


    cout << "Number of H Atoms in the fuel species  " << hAtoms << endl;
    cout << "Number of C Atoms in the fuel species  " << cAtoms << endl;
    cout << "avg mol wt of fuel " << avgMolWt << endl;
    cout << "Total fuel mass fraction " << fuelMassFrac << endl;
    cout << "O2 mass frac " << temp[iO2] << endl;

    doublereal stO2 = cAtoms + hAtoms/4.0;
    
    /*
     *stoichiometric mass ratio
     */
    smr = stO2*0.032/avgMolWt;
    cout << "Avg mol wt " << avgMolWt << endl;
    /*
     *stoichiometric mixture fraction
     */
    stoichZ = 1.0/(1+ smr*fuelMassFrac/o2MassFrac);

    cout << "Stoichiometric mixture fraction " << stoichZ << endl;

    return stoichZ;

}
/*
 *calculate the scalar dissipation rate
 */
doublereal FlameLet::scalarDissipationRate(const doublereal m_frac){
    /*
     *Eq. 9.38 SummerSchool by N. Peters
     */
    CamMath cm;
    doublereal erterm = cm.inverfc(2*m_frac);    
    doublereal arg = -2*cm.SQR(erterm);
    strain = admin->getStrainRate();
    sdr = strain*exp(arg)/pi;
    return sdr;
    //return 1e-03;
}

/*
 *solver call for residual evaluation
 */
int FlameLet::eval(doublereal x, doublereal* y, doublereal* ydot, bool jacEval){


    //Sets soot volume fraction vector to zeros
    setExternalSootVolumeFraction(std::vector<doublereal>(mCord, 0.0));
    residual(x,y,ydot);
    return 0;
    
}

void FlameLet::report(doublereal x, doublereal* solution){
    cout << "Dummy function called\n";
}

/*
 *consol output function
 */
void FlameLet::report(doublereal x, doublereal* solution, doublereal& res){
    static int nStep=0;
    cout.width(5);
    cout.setf(ios::scientific);
    if(nStep%10==0) reporter->consoleHead("time(s) \t residual");
    cout << x <<"\t" << res << endl;
    nStep++;

//    cout << "---------------------------------\n";
//    for(int i=0; i<iMesh_e; i++){
//        cout << solution[i*nVar+ptrT] << endl;
//    }
//    cout << "---------------------------------\n";
    //int dd ; cin >> dd;
}
/*
 *output function for file output
 */
void FlameLet::reportToFile(doublereal t, doublereal* soln){
    
    doublereal sum;
    reporter->openFiles();
    reporter->writeCustomHeader(headerData);
    vector<doublereal> data, axpos;
    vector<doublereal> molfrac, massfrac;
    axpos = reacGeom->getAxpos();
    int len = axpos.size();
    for(int i=0; i<len; i++){
        data.clear();
        data.push_back(t);
        data.push_back(axpos[i]);
        data.push_back(scalarDissipationRate(axpos[i]));
        data.push_back(m_rho[i]);
        data.push_back(soln[i*nVar+ptrT]);
        data.push_back(radiation[i]); 

        massfrac.clear();
        molfrac.clear();
        for(int l=0; l<nSpc; l++){
            massfrac.push_back(soln[i*nVar+l]);
        }
        if(admin->getSpeciesOut() == admin->MASS){
            sum = 0;
            for(int l=0; l<nSpc; l++){
                data.push_back(fabs(massfrac[l]));
                sum += massfrac[l];
            }
        }else{
            CamConverter cc;
            cc.mass2mole(massfrac,molfrac,*camMech);
            sum = 0;
            for(int l=0; l<nSpc; l++){
                data.push_back(fabs(molfrac[l]));
                sum += molfrac[l];
            }
        }
        data.push_back(sum);
        reporter->writeCustomFileOut(data);

    }

    reporter->closeFiles();

}
/*
 *output file header
 */
void FlameLet::header(){
    headerData.clear();
    headerData.push_back("time");
    headerData.push_back("Z");
    headerData.push_back("SDR");
    headerData.push_back("rho");
    headerData.push_back("T");
    headerData.push_back("Radiation");
    for (int l = 0; l < nSpc; l++) {
        headerData.push_back( (*spv)[l]->Name() );
    }
    headerData.push_back("sumfracs");
}


/*
 *set the scalar dissipation rate provided by the external
 *calling program
 */
void FlameLet::setExternalScalarDissipationRate(const doublereal sr){
    sdr_ext = sr;
}

/**
 *  When the scalar dissipation rate has a time history
 *  use that during intergration
 */
void FlameLet::setExternalScalarDissipationRate(const vector<doublereal>& time, const vector<doublereal>& sdr){
    v_sdr = sdr;
    v_time = time;
    
    timeHistory = true;

}

/**
 *  Interpolate and return the scalar dissipation rate
 */
doublereal FlameLet::getSDR(const doublereal time) const {
    doublereal tsdr=0;
    doublereal vu, vl, xu, xl;
    for(size_t i =0; i < v_sdr.size(); i++){
        if(time == v_time[i]){
            tsdr = v_sdr[i];
            
            break;
        }else if( i>0 && (time > v_time[i-1]) && (time < v_time[i])) {
            vu = v_sdr[i];
            xu = v_time[i];

            vl = v_sdr[i-1];
            xl = v_time[i-1];

            doublereal slope = (vu-vl)/(xu-xl);
            doublereal intersect = vu- (slope*xu);

            tsdr =  slope*time + intersect;            
            break;
        }else{
            size_t l = v_sdr.size();
            vu = v_sdr[l-1];
            xu = v_time[l-1];

            vl = v_sdr[l-2];
            xl = v_time[l-2];

            doublereal slope = (vu-vl)/(xu-xl);
            doublereal intersect = vu- (slope*xu);

            tsdr =  slope*time + intersect;
            break;
        }
    }
    return tsdr;
}

/*!
 *@param[in]    soot_fv     Vector of soot volume fractions, one for each grid cell
 */
void FlameLet::setExternalSootVolumeFraction(const std::vector<doublereal>& soot_fv) {
    m_SootFv = soot_fv;

    if (m_T.size()!= m_SootFv.size()){
        throw std::runtime_error("soot volume fraction vector is of the wrong size");
     }
}
