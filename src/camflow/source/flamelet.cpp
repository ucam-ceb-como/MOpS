#include "flamelet.h"

using namespace Camflow;
using namespace std;
using namespace Gadgets;

FlameLet::FlameLet
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
  CamSetup(ca, config, cc, cg, cp, cs, mech),
  sdr_ext(0),
  timeHistory(false),
  sdrProfile(false),
  sdrAnalytic(false),
  scalarDissipationRate_(mCord, 1),
  CpSpec(mCord,nSpc)
{}

FlameLet::~FlameLet()
{
    //if (radiation != NULL) delete radiation;
}

/*
 *this is called by the model object. The boolean interface decides
 *if the call originates from the interface or from camflow kernel
 */
void FlameLet::setRestartTime(doublereal t){
    restartTime = t;
}

void FlameLet::solve()
{
    solve(false);
}

void FlameLet::solve
(
    bool interface
)
{

    reacGeom_.addZeroWidthCells();

    /*
     *init the solution vector
     */
    initSolutionVector();

    reporter_->header("Flamelet");
    header();

    if(!interface)
    {
        reportToFile("initialProfile.dat",control_.getMaxTime(),&solvect[0]);
    }

    if (control_.getSolutionMode() == control_.COUPLED)
    {
        csolve(interface);
    }
    else
    {
        ssolve();
        csolve(interface);
    }

}


/**
 *continuation call from an external code that
 *solves for population balance
 */
void FlameLet::solve
(
    vector<Thermo::Mixture>& cstrs,
    const vector<vector<doublereal> >& iniSource,
    const vector<vector<doublereal> >& fnlSource,
    Mechanism& mech,
    CamControl& cc,
    CamAdmin& ca,
    CamGeometry& cg,
    CamProfile& cp
)
{

    reacGeom_.addZeroWidthCells();

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
    admin_.getLeftBoundary(left);
    admin_.getRightBoundary(right);
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


void FlameLet::initSolutionVector()
{

    /*
     *left boundary is for the fuel and right boundary is
     *for oxidizer
     */
    CamBoundary left, right;
    admin_.getLeftBoundary(left);
    admin_.getRightBoundary(right);
    storeInlet(left,fuel);
    storeInlet(right,oxid);
    fuel.T = left.getTemperature();
    oxid.T = right.getTemperature();

    stoichiometricMixtureFraction();

    profile_.setMixingCenter(stoichZ);
    profile_.setMixingWidth(0.5*stoichZ);
    /*
     *initialize the ODE vector
     */
    solvect.resize(nEqn,0.0);
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
    initSpecies(right,left,control_,vSpec);
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
    vector<doublereal> position = reacGeom_.getAxpos();
    int len = position.size();
    vT.resize(len,0.0);
    for (unsigned int i=0; i<dz.size();i++)
    {
       // vT[i] = profile->getUserDefTemp(position[i]);

        if (position[i] < stoichZ)
        {
            vT[i] = slopeOx*position[i] + inrsctOx;
        }
        else
        {
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

}

/*
 *coupled solver
 */
void FlameLet::csolve
(
    bool interface
)
{

    if (solverID == control_.CVODE)
    {
        CVodeWrapper cvw;
        eqn_slvd = EQN_ALL;
        int band = nVar*2;

        cvw.init(nEqn,solvect,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
            control_.getMaxTime(),band,*this);

        cvw.solve(CV_ONE_STEP,control_.getResTol());
        /*
         *write the output to file only if the call is not
         *from the interface
         */
        if(!interface)
        {
            reportToFile("profile.dat",control_.getMaxTime(),&solvect[0]);
        }

    }
    else if (solverID == control_.NEWTON)
    {
        std::cout << "Not implemented\n";
    }
    else if (solverID == control_.RADAU)
    {

        RadauWrapper radauWrapper;

        eqn_slvd = EQN_ALL;

        std::vector<doublereal> relTolVector;
        std::vector<doublereal> absTolVector;

        relTolVector.push_back(control_.getSpeciesRelTol());
        absTolVector.push_back(control_.getSpeciesAbsTol());

        radauWrapper.setBandWidth(nVar);

    	radauWrapper.initSolver(nEqn,
                                0.0,
                                control_.getMaxTime(),
                                solvect,
                                relTolVector,
                                absTolVector,
                                *this);

        radauWrapper.Integrate();

        /*
         *write the output to file only if the call is not
         *from the interface
         */
        if(!interface) {
            reportToFile("profile.dat", control_.getMaxTime(),&solvect[0]);
        }
       // radauWrapper.destroy();

    }
    else if (solverID == control_.LIMEX)
    {
        throw std::logic_error("Error -- Limex is not yet supported");
    }

}

/*
 *mass matrix evaluation
 */
void FlameLet::massMatrix(doublereal** M)
{}

/*
 *restart the solution. This is normally called from the interface routine
 *The solver is reinitialized each time with the previous solution.
 */
void FlameLet::restart()
{

    if (solverID == control_.CVODE){
        CVodeWrapper cvw;
        eqn_slvd = EQN_ALL;
        int band = nVar*2;

        cvw.init(nEqn,solvect,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
            control_.getMaxTime(),band,*this,restartTime);
        cvw.solve(CV_ONE_STEP,control_.getResTol());
    }

     else if (solverID == control_.LIMEX) {
         throw std::logic_error("Error -- Limex is not yet supported");
     }

}

/*
 *segregated solver
 */
void FlameLet::ssolve()
{
    
    int seg_eqn, band;
    vector<doublereal> seg_soln_vec;

    if ( solverID == control_.CVODE){
    
       CVodeWrapper cvw;
       
       for (int i=0; i<control_.getNumIterations();i++){
        /*
         *solve species equation
         */
        cout << "Solving species equations  " << i << endl;
        int dd; cin >> dd;
        eqn_slvd = EQN_SPECIES;
        seg_eqn = nSpc*mCord;
        band = nSpc*2;
        extractSpeciesVector(seg_soln_vec);
        cvw.init(seg_eqn,seg_soln_vec,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
            control_.getMaxTime(),band,*this);
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
        cvw.init(seg_eqn,seg_soln_vec,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
            control_.getMaxTime(),band,*this);
        cvw.solve(CV_ONE_STEP,1e-03);
        mergeEnergyVector(&seg_soln_vec[0]);
        cvw.destroy();
      }
    }

    else if (solverID == control_.LIMEX) {
        throw std::logic_error("Error -- Limex is not yet supported");
    }
  }

  
/*
 *residual definitions
 */
void FlameLet::residual
(
    const doublereal& t,
    doublereal* y,
    doublereal* f
)
{

    if (eqn_slvd == EQN_ALL) // coupled solver part
    {
        saveMixtureProp(y);
        speciesResidual(t,y,&resSp[0]);
        energyResidual(t,y,&resT[0]);

        for (int i=0; i<mCord; ++i)
        {
            for (int l=0; l<nSpc; ++l)
            {
                f[i*nVar+l] = resSp[i*nSpc+l];
            }
            f[i*nVar+ptrT] = resT[i];
        }
    }
    else // segregated solver part
    {
        if (eqn_slvd == EQN_SPECIES)
        {
            mergeSpeciesVector(y);
            saveMixtureProp(&solvect[0]);
            speciesResidual(t,y,f);
        }
        else if (eqn_slvd==EQN_ENERGY)
        {
            mergeEnergyVector(y);
            saveMixtureProp(&solvect[0]);
            energyResidual(t,y,f);
        }
    }

    // Refine Grid
    //if(getResidual() < 1e-6)
    //{
    //    reacGeom->refine(y,nVar,nSpc,ptrT);
    //}

}

doublereal FlameLet::getResidual()
const
{

	doublereal resNorm=0;

    for (int i=0; i<nSpc*mCord; ++i)
    {
        resNorm += resSp[i]*resSp[i];
    }
    for (int i=0; i<mCord; ++i)
    {
        resNorm += resT[i]*resT[i];
    }

    return std::sqrt(resNorm);

}

/*
 *species residual definitions
 */
void FlameLet::speciesResidual
(
    const doublereal& t,
    doublereal* y,
    doublereal* f
)
{

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
    if (timeHistory) {

        sdr = getSDR(t);

    } else if (sdr_ext!=0) {

        sdr=sdr_ext;

    }

    //sdr = scalarDissipationRate(dz[i]);
    for (int l=0; l<nSpc; ++l)
    {
        f[l] = 0.0;
    }

    /**
     *  For non-unity Lewis numbers
     */

    Le.resize(mCord,nSpc,1.0);
    convection.resize(mCord,nSpc,0);
    doublereal onebLe;
    doublereal oneby16 = 1.0/16;
    if (Lewis == FlameLet::LNNONE)
    {
        //Iterating over the interior points
        for (int i=1; i<mCord-1; ++i)
        {
            for (int l=0; l<nSpc; ++l)
            {
                Le(i,l) = m_k[i]/(m_rho[i]*m_cp[i]*s_Diff(i,l));
                onebLe = 1/Le(i,l);

                convection(i,l) = oneby16*(onebLe-1)*(s_mf(i+1,l)-s_mf(i-1,l))*(m_rho[i+1]-m_rho[i-1]);
            }
        }
        for (int l=0; l<nSpc; ++l)
        {
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
    for (int i=iMesh_s; i<iMesh_e; ++i)
    {

        zPE = 0.5*(dz[i]+dz[i+1]);
        zPW = 0.5*(dz[i]+dz[i-1]);

        //if(sdr_ext==0)sdr = scalarDissipationRate(dz[i]);
        //if(sdrProfile)sdr = getSDRfromProfile(t,dz[i]);
        //if(sdrAnalytic)sdr = scalarDissipationRateProfile(dz[i],getSDR(t),i);
        sdr = scalarDissipationRateProfile(dz[i],2.0,i);

        doublereal diffusionConstant = sdr/(2.0*dz[i]);

        for (int l=0; l<nSpc; ++l)
        {
            grad_e = (s_mf(i+1,l)-s_mf(i,l))/zPE;
            grad_w = (s_mf(i,l)-s_mf(i-1,l))/zPW;
            source = s_Wdot(i,l)*(*spv_)[l]->MolWt()/m_rho[i];
            f[i*nSpc+l] = diffusionConstant*(grad_e-grad_w)/Le(i,l)
                          + source;
        }

    }

    /*
     *Mixture fraction 1. The fuel inlet. Concentrations are
     *held constant at the fuel composition
     */

    for (int l=0; l<nSpc; ++l)
    {
        f[iMesh_e*nSpc+l] = 0.0;
    }

}

/*!
 *energy residual
 *
 */
void FlameLet::energyResidual
(
    const doublereal& t,
    doublereal* y,
    doublereal* f
)
{

    doublereal grad_e=0, grad_w=0;
    doublereal zPE=0, zPW=0;
    doublereal source=0;

    if(admin_.getRadiationModel())
    {
        radiation = new Radiation(mCord);
    }

    /*
     *starting with mixture fraction zero: i.e oxidizer
     *inlet. The temperature is fixed at the oxidizer
     *inlet temperature
     */

    f[0] = 0.0;

    /*
     *set the externally specified scalar dissipation rate
     *to sdr
     */
    if(timeHistory){
        sdr = getSDR(t);

    } else if(sdr_ext!=0){
        sdr=sdr_ext;
    }

    /*
     *intermediate mixture fraction coordinates
     */

    for (int i=iMesh_s; i<iMesh_e; ++i)
    {

        zPE = 0.5*(dz[i]+dz[i+1]);
        zPW = 0.5*(dz[i]+dz[i-1]);
        source = 0.0;

        for (int l=0; l<nSpc; ++l)
        {
            source += s_Wdot(i,l)*s_H(i,l);
        }

        //if(sdr_ext==0)sdr = scalarDissipationRate(dz[i]);
        //if(sdrProfile)sdr = getSDRfromProfile(t,dz[i]);
        //if(sdrAnalytic)sdr = scalarDissipationRateProfile(dz[i],getSDR(t),i);
        sdr = scalarDissipationRateProfile(dz[i],2.0,i);

        grad_e = (m_T[i+1]-m_T[i])/zPE;
        grad_w = (m_T[i]-m_T[i-1])/zPW;

        f[i] = sdr*(grad_e-grad_w)/(2.0*dz[i])
               - source/(m_rho[i]*m_cp[i]);

        /**
         *  Accounting for non-unity Lewis number
         */
        if (Lewis == FlameLet::LNNONE)
        {

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

        if (admin_.getRadiationModel())
        {

            doublereal mole_fracsH2O;
            doublereal mole_fracsCO2;
            doublereal mole_fracsCO;

            //Get species indexes corresponding to H20, CO2, CO
            if(camMech_->FindSpecies("H2O") == -1){
                mole_fracsH2O = 0.0;
            } else {
                const int iH2O = camMech_->FindSpecies("H2O");
                const doublereal molwtH2O =   (*spv_)[iH2O] -> MolWt();
                mole_fracsH2O = s_mf(i,iH2O)*avgMolWt[i]/molwtH2O;
            }
            if(camMech_->FindSpecies("CO2") == -1){
                mole_fracsCO2 = 0.0;
            } else {
                const int iCO2 = camMech_->FindSpecies("CO2");
                const doublereal molwtCO2 =   (*spv_)[iCO2] -> MolWt();
                mole_fracsCO2 = s_mf(i,iCO2)*avgMolWt[i]/molwtCO2;
            }
            if(camMech_->FindSpecies("CO") == -1){
                mole_fracsCO = 0.0;
            } else {
                const int iCO  = camMech_->FindSpecies("CO");
                const doublereal molwtCO =    (*spv_)[iCO] -> MolWt();
                mole_fracsCO = s_mf(i,iCO)*avgMolWt[i]/molwtCO;
            }

            //This radiation term is sent to as output to profile.h
            radiation->RadiativeLoss
            (
               i,
               m_T[i],
               opPre,
               m_SootFv[i],
               mole_fracsH2O,
               mole_fracsCO2,
               mole_fracsCO
            );

            // This is the new energy residual term, accounting for radiation.
            // This is DEFINITELY NEGATIVE!
            f[i] -= radiation->getRadiation(i)/(m_rho[i]*m_cp[i]);

        }

    }

    f[iMesh_e] = 0.0;

}
/**
 *save the mixture property
 * \todo This should be trivially parallelisable. A lot of time is spent here.
 */
void FlameLet::saveMixtureProp(doublereal* y)
{

    vector<doublereal> mf;
    vector<doublereal> htemp(nSpc,0.0);
    vector<doublereal> temp(nSpc,0.0);
    vector<doublereal> cptemp(nSpc,0.0);

    for (int i=0; i<mCord; ++i)
    {

        mf.clear();
        for(int l=0; l<nSpc; l++)
        {
            mf.push_back(y[i*nVar+l]);
        }

        m_T[i] = y[i*nVar+ptrT];                                   //temperature
        camMixture_->SetMassFracs(mf);                              //mass fraction
        camMixture_->SetTemperature(m_T[i]);                        //temperature

        avgMolWt[i] = camMixture_->getAvgMolWt();
        m_rho[i] = opPre*avgMolWt[i]/(R*m_T[i]);                   //density

        camMixture_->SetMassDensity(m_rho[i]);                      //density

        camMech_->Reactions().GetMolarProdRates(*camMixture_,wdot);

        htemp = camMixture_->getMolarEnthalpy();                    //enthalpy
        m_cp[i] = camMixture_->getSpecificHeatCapacity();           //specific heat
        m_k[i] = camMixture_->getThermalConductivity(opPre);        //thermal conductivity

        // MOVE THIS OUTSIDE LOOP
        m_mu[i] = camMixture_->getViscosity();                      //mixture viscosity

        temp = camMixture_->getMixtureDiffusionCoeff(opPre);
        cptemp = camMixture_->getMolarSpecificHeat();
        for(int l=0; l<nSpc; l++)
        {
            s_mf(i,l) = mf[l];
            s_Wdot(i,l) = wdot[l];
            s_H(i,l) = htemp[l];
            s_Diff(i,l) = temp[l];
            //Specific heat capacity of species in J/Kg K
            CpSpec(i,l) =cptemp[l]/(*spv_)[l]->MolWt();
        }

    }
}

doublereal FlameLet::stoichiometricMixtureFraction()
{
    /*
     *check for C and H atoms
     */

    vector<Sprog::Species*> fspecies, ospecies;

    int indx_C = camMech_->FindElement("C");
    int indx_H = camMech_->FindElement("H");

    /*
     *fuel inlet
     */
    CamBoundary fuelInlet, oxInlet;
    map<string, doublereal> species;
    map<string, doublereal>::iterator sIterator;
    admin_.getLeftBoundary(fuelInlet);
    species = fuelInlet.getInletSpecies();
    sIterator = species.begin();

    while(sIterator != species.end()){
        fspecies.push_back(camMech_->GetSpecies(sIterator->first));
        sIterator++;
    }

    admin_.getRightBoundary(oxInlet);
    species = oxInlet.getInletSpecies();
    sIterator = species.begin();
    while(sIterator != species.end()){
        ospecies.push_back(camMech_->GetSpecies(sIterator->first));
        sIterator++;
    }

    int cAtoms=0;
    int hAtoms=0;
    doublereal avgMolWt=0;
    doublereal fuelMassFrac=0;
    unsigned int i;
    vector<doublereal> temp;
    getInletMassFrac(fuelInlet,temp);

    int iN2 = camMech_->FindSpecies("N2");
    int iAR = camMech_->FindSpecies("AR");
    int iHe = camMech_->FindSpecies("HE");
    int iO2 = camMech_->FindSpecies("O2");

    for(i=0; i<fspecies.size();i++){
        int icAtoms = fspecies[i]->AtomCount(indx_C);
        int ihAtoms = fspecies[i]->AtomCount(indx_H);
        if (icAtoms != 0 || ihAtoms !=0) avgMolWt += fspecies[i]->MolWt();
        cAtoms += icAtoms; hAtoms += ihAtoms;

        int spIndx = camMech_->FindSpecies(fspecies[i]->Name());
        if(spIndx != iN2 && spIndx != iAR && spIndx != iHe){
            fuelMassFrac += temp[spIndx];
        }
    }

    getInletMassFrac(oxInlet,temp);
    doublereal o2MassFrac = temp[iO2];


    cout << "Number of H Atoms in the fuel species  " << hAtoms << endl;
    cout << "Number of C Atoms in the fuel species  " << cAtoms << endl;
    cout << "avg mol wt of fuel " << avgMolWt << endl;
    cout << "nEqn " << nEqn << endl;
    cout << "Total fuel mass fraction " << fuelMassFrac << endl;
    cout << "O2 mass frac " << temp[iO2] << endl;

    doublereal stO2 = cAtoms + hAtoms/4.0;

    /*
     *stoichiometric mass ratio
     */
    doublereal smr = stO2*0.032/avgMolWt;
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
doublereal FlameLet::scalarDissipationRate(const doublereal m_frac)
{
    /*
     *Eq. 9.38 SummerSchool by N. Peters
     */
    CamMath cm;
    doublereal erterm = cm.inverfc(2*m_frac);
    doublereal arg = -2*cm.SQR(erterm);
    sdr = admin_.getStrainRate()*exp(arg)/PI;
    return sdr;
    //return 0.88;
}

/*
 *calculate the scalar dissipation rate profile. Method 1 in Carbonell(2009).
 */
doublereal FlameLet::scalarDissipationRateProfile
(
    const doublereal m_frac,
    const doublereal stoichSDR,
    const int cell
)
{

	CamMath cm;

	doublereal fZ = exp(-2*cm.SQR(cm.inverfc(2*m_frac)));
	doublereal fZst = exp(-2*cm.SQR(cm.inverfc(2*stoichZ)));

	Utils::LinearInterpolator<doublereal, doublereal> rhoInterpolate(reacGeom_.getAxpos(),m_rho);
	doublereal rhoStoich = rhoInterpolate.interpolate(stoichZ);

	doublereal phi = 0.75 *
	                 (
	                     cm.SQR(std::sqrt(m_rho[0]/m_rho[cell])+1.0)
	                   / (2.0*std::sqrt(m_rho[0]/m_rho[cell])+1.0)
	                 );
	doublereal phist = 0.75 *
	                   (
	                       cm.SQR(std::sqrt(m_rho[0]/rhoStoich)+1.0)
	                     / (2.0*std::sqrt(m_rho[0]/rhoStoich)+1.0)
	                   );

	return stoichSDR * (fZ/fZst);// * (phi/phist);

}

/*
 *solver call for residual evaluation
 */
int FlameLet::eval
(
    doublereal x,
    doublereal* y,
    doublereal* ydot,
    bool jacEval
)
{


    //Sets soot volume fraction vector to zeros
    setExternalSootVolumeFraction(std::vector<doublereal>(mCord, 0.0));
    residual(x,y,ydot);
    return 0;

}

/*
 *consol output function
 */
void FlameLet::report(doublereal x, doublereal* solution, doublereal& res)
{

    static int nStep=0;
    cout.width(5);
    cout.setf(ios::scientific);
    //if(nStep%10==0) reporter->consoleHead("time(s) \t residual");
    cout << x <<"\t" << res << endl;
    nStep++;

}
/*
 *output function for file output
 */
void FlameLet::reportToFile(std::string fileName, doublereal t, doublereal* soln)
{

    doublereal sum;
    reporter_->openFile(fileName,false);
    reporter_->writeCustomHeader(headerData);
    vector<doublereal> data, axpos;
    vector<doublereal> molfrac, massfrac, temperatureVec;
    axpos = reacGeom_.getAxpos();
    int len = axpos.size();

    for(int i=0; i<len; i++)
    {
        temperatureVec.push_back(soln[i*nVar+ptrT]);
    }
    reporter_->writeTempProfiletoXML("camflow.xml",temperatureVec);

    for (int i=0; i<len; i++) {

        data.clear();
        data.push_back(t);
        data.push_back(axpos[i]);
        data.push_back(scalarDissipationRateProfile(axpos[i],2.0,i));
        data.push_back(m_rho[i]);
        data.push_back(m_mu[i]);
        data.push_back(m_cp[i]);
        data.push_back(soln[i*nVar+ptrT]);
        if (admin_.getRadiationModel())
        {
            data.push_back(radiation->getRadiation(i));
        }
        else
        {
            data.push_back(0.0);
        }

        massfrac.clear();
        molfrac.clear();
        for(int l=0; l<nSpc; l++){
            massfrac.push_back(soln[i*nVar+l]);
        }
        if(admin_.getSpeciesOut() == admin_.MASS){
            sum = 0;
            for(int l=0; l<nSpc; l++){
                data.push_back(fabs(massfrac[l]));
                sum += massfrac[l];
            }
        }else{
            CamConverter cc;
            cc.mass2mole(massfrac,molfrac,*camMech_);
            sum = 0;
            for(int l=0; l<nSpc; l++){
                data.push_back(fabs(molfrac[l]));
                sum += molfrac[l];
            }
        }
        data.push_back(sum);
        reporter_->writeCustomFileOut(data);

    }

    reporter_->closeFile();

    //reacGeom->refine(soln,nVar,nSpc,ptrT);

}
/*
 *output file header
 */
void FlameLet::header()
{

    headerData.clear();
    headerData.push_back("time");
    headerData.push_back("Z");
    headerData.push_back("SDR");
    headerData.push_back("rho");
    headerData.push_back("mu");
    headerData.push_back("cp");
    headerData.push_back("T");
    headerData.push_back("Radiation");
    for (int l = 0; l < nSpc; l++) {
        headerData.push_back( (*spv_)[l]->Name() );
    }
    headerData.push_back("sumfracs");

}


/*
 *set the scalar dissipation rate provided by the external
 *calling program
 */
void FlameLet::setExternalScalarDissipationRate(const doublereal sr)
{
    sdr_ext = sr;
}

/**
 *  When the scalar dissipation rate has a time history
 *  use that during intergration
 */
void FlameLet::setExternalScalarDissipationRate
(
    const std::vector<doublereal>& time,
    const std::vector<doublereal>& sdr,
    const bool analytic
)
{

    v_sdr = sdr;
    v_time = time;

    timeHistory = true;
    sdrAnalytic = analytic;

}

/**
 *  When the scalar dissipation rate has a time history
 *  and has a profile with mixture fraction from the CFD.
 */
void FlameLet::setExternalScalarDissipationRate
(
    const std::vector<doublereal>& time,
	const std::vector< std::vector<doublereal> >& sdr,
	const std::vector< std::vector<doublereal> >& Zcoords
)
{

    profile_sdr = sdr;
    v_time = time;
    cfdMixFracCoords = Zcoords;

    sdrProfile = true;

}

/**
 *  Interpolate and return the scalar dissipation rate
 */
doublereal
FlameLet::getSDR(const doublereal time)
const
{

	Utils::LinearInterpolator<doublereal, doublereal> timeInterpolate(v_time, v_sdr);

	return timeInterpolate.interpolate(time);

}

/**
 *  Interpolate and return the scalar dissipation rate from a profile that varies through time.
 */
doublereal
FlameLet::getSDRfromProfile
(
    const doublereal time,
    const doublereal Z
)
const
{

	std::vector<doublereal> sdrTime, sdrInterpolated;
	std::vector<doublereal> cfdMixFracCoordsTime, cfdMixFracCoordsInterpolated;

	sdrInterpolated.clear();
	cfdMixFracCoordsInterpolated.clear();

	for (size_t i=0; i<cfdMixFracCoords[0].size(); ++i) {
		sdrTime.clear();
		sdrTime.push_back(profile_sdr[0][i]);
		sdrTime.push_back(profile_sdr[1][i]);

		cfdMixFracCoordsTime.clear();
		cfdMixFracCoordsTime.push_back(cfdMixFracCoords[0][i]);
		cfdMixFracCoordsTime.push_back(cfdMixFracCoords[1][i]);

		Utils::LinearInterpolator<doublereal, doublereal> timeInterpolate(v_time, sdrTime);
		doublereal sdrInterpolatedTime = timeInterpolate.interpolate(time);

		Utils::LinearInterpolator<doublereal, doublereal> time2Interpolate(v_time, cfdMixFracCoordsTime);
		doublereal cfdMixFracCoordsInterpolatedTime = time2Interpolate.interpolate(time);

		sdrInterpolated.push_back(sdrInterpolatedTime);
		cfdMixFracCoordsInterpolated.push_back(cfdMixFracCoordsInterpolatedTime);
	}

	Utils::LinearInterpolator<doublereal, doublereal> spaceInterpolate(cfdMixFracCoordsInterpolated, sdrInterpolated);

	return spaceInterpolate.interpolate(Z);

}

/*!
 *@param[in]    soot_fv     Vector of soot volume fractions, one for each grid cell
 */
void
FlameLet::setExternalSootVolumeFraction(const std::vector<doublereal>& soot_fv)
{
    m_SootFv = soot_fv;

    //Solver assumes one soot volume fraction for each cell; temperature is already initialized
    /*if (static_cast<size_t>(reacGeom->getnCells()) != m_SootFv.size()) {
        std::ostringstream msg ("new soot volume fraction vector length is ");
        msg << m_SootFv.size() << " but geometry length is " << reacGeom->getnCells();
        throw std::runtime_error (msg.str());
    }*/
}

/*
 *  Return the pyrene(A4) molar production rate term.
 */
void
FlameLet::getWdotA4(std::vector<doublereal>& wdotA4)
const
{

	wdotA4.clear();
	// Check the species exists first (returns -1 if it does not).
	if(camMech_->FindSpecies("A4") == -1){
		std::cout << "Species A4 not found." << std::endl;
		for(int i=iMesh_s; i<iMesh_e;i++){
			wdotA4.push_back(0.0);
		}
	} else {
		const int iA4 = camMech_->FindSpecies("A4");
		for(int i=iMesh_s; i<iMesh_e;i++){
			wdotA4.push_back(s_Wdot(i,iA4));
		}
	}

}
