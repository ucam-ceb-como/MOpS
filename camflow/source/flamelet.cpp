
#include "cam_setup.h"
#include "cam_residual.h"
#include <vector>
#include "cam_reporter.h"
#include <iostream>
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

    if(reacGeom->getnCells()==0) reacGeom->discretize();

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
                            cc.getMaxTime(),band,*this);
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
            f[i*nSpc+l] = sdr*(grad_e-grad_w)/(2*dz[i]) + source;
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
        

    }

    
    f[iMesh_e] = 0.0;
    
    
}
/*
 *save the mixture proerty
 */
void FlameLet::saveMixtureProp(doublereal* y){
    s_mf.resize(mCord,nSpc);
    s_Wdot.resize(mCord,nSpc);
    s_H.resize(mCord,nSpc);
    
    m_T.resize(mCord,0.0);
    m_rho.resize(mCord,0.0);
    m_cp.resize(mCord,0.0);
    m_mu.resize(mCord,0.0);
    m_u.resize(mCord,0.0);

    vector<doublereal> mf,htemp;
    htemp.resize(nSpc,0.0);
    for(int i=0; i<mCord; i++){
        mf.clear();
        for(int l=0; l<nSpc; l++){
            mf.push_back(y[i*nVar+l]);            
        }
        m_T[i] = y[i*nVar+ptrT];                                //temperature
        camMixture->SetMassFracs(mf);                           //mass fraction
        camMixture->SetTemperature(m_T[i]);                     //temperature
        m_rho[i] = opPre*camMixture->getAvgMolWt()/(R*m_T[i]);  //density
        camMixture->SetMassDensity(m_rho[i]);                   //density
        camMech->Reactions().GetMolarProdRates(*camMixture,wdot);
        htemp = camMixture->getMolarEnthalpy();                 //enthalpy
        m_cp[i] = camMixture->getSpecificHeatCapacity();        //specific heat
        m_mu[i] = camMixture->getViscosity();                   //mixture viscosity
        for(int l=0; l<nSpc; l++){
            s_mf(i,l) = mf[l];
            s_Wdot(i,l) = wdot[l];
            s_H(i,l) = htemp[l];
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
    reporter->writeHeader(headerData);
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
        reporter->writeStdFileOut(data);

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
    v_sdr.resize(time.size());
    v_time.resize(time.size());

    v_sdr = sdr;
    v_time = time;
    timeHistory = true;

}

/**
 *  Interpolate and return the scalar dissipation rate
 */
const doublereal FlameLet::getSDR(const doublereal time){
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
