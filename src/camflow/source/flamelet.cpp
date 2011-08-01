#include "flamelet.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/lexical_cast.hpp>

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
  stoichZ(stoichiometricMixtureFraction()),
  timeHistory(false),
  sdrProfile(false),
  sdrAnalytic(false),
  radiation(NULL),
  scalarDissipationRate_(admin_.getInputFile(), stoichZ, reacGeom_.getAxpos(), 1),
  CpSpec(mCord,nSpc)
{}

FlameLet::~FlameLet()
{
    if (radiation != NULL) delete radiation;
}

void FlameLet::checkSetup()
{

    if(reacGeom_.getAxpos()[mCord-1] != 1)
        throw std::invalid_argument("Mixture fraction does not go from 0 to 1. "
                                    "Check <length unit=\"m\">1.0</length>\n");

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

    // Check that the problem has been setup properly.
    checkSetup();

    /*
     *init the solution vector
     */
    initSolutionVector();

    reporter_->header("Flamelet");

    if(!interface)
    {
        reportToFile("initialProfile.dat",control_.getMaxTime(),solvect);
        //writeXMLFile(scalarDissipationRate_.getRefSDR(), solvect);
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

    if (admin_.getRestartType() == admin_.BINARY)
    {
        std::ofstream ofs(admin_.getRestartFile().c_str());
        boost::archive::binary_oarchive oa(ofs);
        oa << reacGeom_.getAxpos() << solvect;
        ofs.close();
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
    CamBoundary& left = admin_.getLeftBoundary();
    CamBoundary& right = admin_.getRightBoundary();
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

    // Initialise the radiation class if necessary.
    if(admin_.getRadiationModel())
    {
        radiation = new Radiation
        (
            admin_.getInputFile(),
            mCord,
            camMech_,
            avgMolWt,
            s_mf
        );
    }

    /*
     *left boundary is for the fuel and right boundary is
     *for oxidizer
     */
    CamBoundary& left = admin_.getLeftBoundary();
    CamBoundary& right = admin_.getRightBoundary();
    storeInlet(left,fuel);
    storeInlet(right,oxid);

    profile_.setMixingCenter(stoichZ);
    profile_.setMixingWidth(0.5*stoichZ);
    /*
     *initialize the ODE vector
     */
    solvect.resize(nEqn,0.0);
    vector<doublereal> vSpec, vT, vMom;
    /*
     * actual signature follows (left,right,cc,vSpec)
     * but in the case of flamelets the species mass fractions
     * are solved for the mixture fraction coordinate, whose
     * direction is taken the same as physical space.
     * z=0 corresponds to oxidizer and z=1 corresponds to fuel
     * therefore the inlets are interchanged here to initialize
     * the species vector properly
     */
    vSpec = initSpecies(right,left);
    /*
     *the following will initialize the temperature vector with
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
    for (size_t i=0; i<dz.size();i++)
    {
        if (position[i] < stoichZ)
        {
            vT[i] = slopeOx*position[i] + inrsctOx;
        }
        else
        {
            vT[i] = slopeFl*position[i] + inrsctFl;
        }
    }

    if(profile_.flagLoadTemp())
    {
        // Loop over all points, EXCLUDING boundaries
        for (size_t i=1; i<dz.size()-1; i++)
        {
            vT[i] = profile_.getUserDefTemp(position[i]);
        }
    }

    if(profile_.flagLoadFracs())
    {
        // Loop over all points, EXCLUDING boundaries
        for (size_t i=cellBegin+1; i<cellEnd-1; ++i)
        {
            for (size_t l=0; l<nSpc; ++l)
            {
                vSpec[i*nSpc+l] = profile_.getUserDefFracs(position[i],(*spv_)[l]->Name());
            }
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

    // Set the initial moment values (interior and boundary)
    // Also initial constants
    if (sootMom_.active())
    {
        vMom.resize(len*nMoments,0.0);
        for (size_t i=0; i<dz.size(); i++)
        {
        	cout << "First Moment " << sootMom_.getFirstMoment() << endl;
        	vMom[i*nMoments] = sootMom_.getFirstMoment();
        	cout << "vMom[i*nMoments]  " << i*nMoments << vMom[i*nMoments] << endl;
        	for (size_t l=1; l<nMoments; ++l)
            {
            	// ank25: Do we need to multiply by 1e6 here?
        		vMom[i*nMoments+l] = vMom[i*nMoments+l-1] + 1e6 * log(doublereal(sootMom_.getAtomsPerDiamer()));
            	cout << "vMom[i*nMoments+l]  " << i*nMoments+l << vMom[i*nMoments+l] << endl;
            }
        }
        // Call the soot constants function
        sootMom_.initMomentsConstants(*camMech_);
    }


    /*
     *create the actual solution vector by merging the species
     *vector, the temperature vector, and soot vector (if present)
     */

    mergeSpeciesVector(&vSpec[0]);
    mergeEnergyVector(&vT[0]);
    if (sootMom_.active())
    {
    	mergeSootMoments(&vMom[0]);
    }



    if (admin_.getRestartType() == admin_.BINARY)
    {
        std::vector<double> solvect_temp, mixFracCoords_temp;
        std::ifstream ifs(admin_.getRestartFile().c_str());
        if (ifs.good())
        {
            boost::archive::binary_iarchive oi(ifs);
            oi >> mixFracCoords_temp >> solvect_temp;

            if (mixFracCoords_temp.size() == reacGeom_.getAxpos().size())
                solvect = solvect_temp;
            else
            {
                throw std::runtime_error
                (
                    "The solution vector is not the same size "
                    "as the old one you are trying to read in. The old solution "
                    "should be interpolated onto the new grid but that function "
                    "has not been written yet.\n");
            }
        }
    }
    else if (admin_.getRestartType() == admin_.TEXT)
    {
        throw std::runtime_error("Text restart files not used yet.");
    }

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

        // Output the tolerances
        // control_.setResTol(1.e-04);
        cout << "Species Abs  Tol: " << control_.getSpeciesAbsTol() << endl;
        cout << "Species Rel Tol: " << control_.getSpeciesRelTol() << endl;
        cout << "Residual Tol: " << control_.getResTol() << endl;

        cvw.init(nEqn,solvect,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
            control_.getMaxTime(),band,*this);

        cvw.solve(CV_ONE_STEP,control_.getResTol());

        // Calculate the mixture viscosity.
        for (int i=0; i<mCord; ++i)
        {
            std::vector<doublereal> mf;
            for(int l=0; l<nSpc; l++)
            {
                mf.push_back(solvect[i*nVar+l]);
            }
            camMixture_->SetMassFracs(mf);
            camMixture_->SetTemperature(m_T[i]);
            camMixture_->SetMassDensity(m_rho[i]);                      //density
            m_mu[i] = camMixture_->getViscosity();                      //mixture viscosity
        }

        /*
         *write the output to file only if the call is not
         *from the interface
         */
        if(!interface)
        {
            reportToFile("profile.dat",control_.getMaxTime(), solvect);
            //writeXMLFile(scalarDissipationRate_.getStoichSDR(), solvect);
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
            reportToFile("profile.dat", control_.getMaxTime(), solvect);
        }

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

    if (solverID == control_.CVODE) {
        CVodeWrapper cvw;
        eqn_slvd = EQN_ALL;
        int band = nVar*2;

        cvw.init(nEqn,solvect,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
            control_.getMaxTime(),band,*this,restartTime);
        cvw.solve(CV_ONE_STEP,control_.getResTol());
        // Calculate the mixture viscosity.
        for (int i=0; i<mCord; ++i)
        {
            std::vector<doublereal> mf;
            for(int l=0; l<nSpc; l++)
            {
                mf.push_back(solvect[i*nVar+l]);
            }
            camMixture_->SetMassFracs(mf);
            camMixture_->SetTemperature(m_T[i]);
            camMixture_->SetMassDensity(m_rho[i]);                      //density
            m_mu[i] = camMixture_->getViscosity();                      //mixture viscosity
        }

        /*
         *write the output to file only if the call is not
         *from the interface
         */
        //if(!interface)
        //{
        string filename = "profile"+boost::lexical_cast<std::string>(restartTime)+".dat";
        reportToFile(filename,control_.getMaxTime(), solvect);
            //writeXMLFile(scalarDissipationRate_.getStoichSDR(), solvect);
        //}
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
        if (sootMom_.active())
        {
        	sootMomentResidual(t,y,&resMom[0]);
        }


        for (int i=0; i<mCord; ++i)
        {
            // Put species residual into master residual
        	for (int l=0; l<nSpc; ++l)
            {
                f[i*nVar+l] = resSp[i*nSpc+l];
            }

        	// Put temperature into master residual
            f[i*nVar+ptrT] = resT[i];

            // Put moments into master residual
        	for (int l=0; l<nMoments; ++l)
            {
                f[i*nVar+ptrT+1+l] = resMom[i*nMoments+l];
            }

        }
    }
    else // segregated solver part   // Soot not implemented in segregated solver
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
    for (int i=0; i<nMoments*mCord; ++i)
    {
    	resNorm += resMom[i]*resMom[i];
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
     *considered as the oxidizer inlet and the concentrations
     *are held constant
     */


    /*
     *set the externally specified scalar dissipation rate
     *to sdr
     */
    if (timeHistory) {

        //sdr = getSDR(t);

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
        //sdr = scalarDissipationRateProfile(reacGeom_.getAxpos()[i],5.0,i);
        sdr = scalarDissipationRate_(reacGeom_.getAxpos()[i],t);

        doublereal diffusionConstant = sdr/(2.0*dz[i]);

        for (int l=0; l<nSpc; ++l)
        {
            grad_e = (s_mf(i+1,l)-s_mf(i,l))/zPE;
            grad_w = (s_mf(i,l)-s_mf(i-1,l))/zPW;
            source = s_Wdot(i,l)/m_rho[i];
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

/*
 * soot residual definitions
 * This is written for the fully simplied soot flamelet equations
 * See Mauss 2006
 */
void FlameLet::sootMomentResidual
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
     *set the externally specified scalar dissipation rate
     *to sdr
     */
    //if (timeHistory) {
    //    sdr = getSDR(t);
    //} else if (sdr_ext!=0) {
    //    sdr=sdr_ext;
    //}

    //sdr = scalarDissipationRate(dz[i]);
    for (int l=0; l<nSpc; ++l)
    {
        f[l] = 0.0;
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
        //sdr = scalarDissipationRateProfile(reacGeom_.getAxpos()[i],2.0,i);
        sdr = scalarDissipationRate_(reacGeom_.getAxpos()[i],t);

        doublereal diffusionConstant = sdr/(2.0*dz[i]);

        for (int l=0; l<nMoments; ++l)
        {
            grad_e = (moments(i+1,l)-moments(i,l))/zPE;
            grad_w = (moments(i,l)-moments(i-1,l))/zPW;
            source = moments_dot(i,l)/m_rho[i];
            f[i*nMoments+l] = diffusionConstant*(grad_e-grad_w)
                          + source;								// LE = 1
        }

    }

    /*
     *Mixture fraction 1. The fuel inlet. Moments held constant here.
     */

    for (int l=0; l<nMoments; ++l)
    {
        f[iMesh_e*nMoments+l] = 0;
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
      //  sdr = getSDR(t);
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
        //sdr = scalarDissipationRateProfile(reacGeom_.getAxpos()[i],5.0,i);
        sdr = scalarDissipationRate_(reacGeom_.getAxpos()[i],t);

        grad_e = (m_T[i+1]-m_T[i])/zPE;
        grad_w = (m_T[i]-m_T[i-1])/zPW;

        f[i] = sdr*(grad_e-grad_w)/(2.0*dz[i])
               - source/(m_rho[i]*m_cp[i]);

        /**
         * Add some extra terms so that we agree with FlameMaster?
         */
        doublereal deltax = 0.5*(dz[i+1]+dz[i]) + 0.5*(dz[i]+dz[i-1]);
        doublereal tGrad = (m_T[i+1]-m_T[i-1])/(deltax);
        doublereal cpGrad = (m_cp[i+1]-m_cp[i-1])/(deltax);
        doublereal sumYGrad = 0.0;
        for (int l=0; l<nSpc; ++l)
        {
            sumYGrad += CpSpec(i,l) * (s_mf(i+1,l)-s_mf(i-1,l))/(deltax);
        }
        f[i] += (sdr/(2.0*m_cp[i])) * tGrad * (cpGrad + sumYGrad);

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
            //This radiation term is sent to as output to profile.h
            radiation->calculateRadiativeHeatLoss
            (
               i,
               m_T[i],
               opPre,
               m_SootFv[i]
            );

            // This is the new energy residual term, accounting for radiation.
            // This is DEFINITELY NEGATIVE!
            f[i] -= radiation->getRadiation(i)/(m_rho[i]*m_cp[i]);
        }

    }

    // Hold temperature constant at fuel inlet
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
    vector<doublereal> moments_dot_temp(nMoments,0.0);
    vector<doublereal> mom_temp(nMoments,0.0);
    vector<doublereal> conc(nSpc,0.0);
    vector<doublereal> wdotSootGasPhase(nMoments,0.0);

    for (int i=0; i<mCord; ++i)
    {

        // Extract the mass fractions from the solution vector
    	mf.clear();
        for(int l=0; l<nSpc; l++)
        {
            mf.push_back(y[i*nVar+l]);
        }

        // Extract temperature from the solution vector
        m_T[i] = y[i*nVar+ptrT];

        // Extract the moments from solution vector
        mom_temp.clear();
        for(int l=0; l<nSpc; l++)
        {
        	mom_temp.push_back(y[i*nVar+ptrT+1+l]);
        }

        camMixture_->SetMassFracs(mf);                              //mass fraction
        camMixture_->SetTemperature(m_T[i]);                        //temperature
        camMixture_->GetConcs(conc);								//molar conc used by soot

        avgMolWt[i] = camMixture_->getAvgMolWt();
        m_rho[i] = opPre*avgMolWt[i]/(R*m_T[i]);                    //density
        camMixture_->SetMassDensity(m_rho[i]);                      //density
        camMech_->Reactions().GetMolarProdRates(*camMixture_,wdot);
        htemp = camMixture_->getMolarEnthalpy();                    //enthalpy
        m_cp[i] = camMixture_->getSpecificHeatCapacity();           //specific heat
        m_k[i] = camMixture_->getThermalConductivity(opPre);        //thermal conductivity

        // MOVE THIS OUTSIDE LOOP TO CSOLVE
        //m_mu[i] = camMixture_->getViscosity();                      //mixture viscosity
        if (Lewis == FlameLet::LNNONE) temp = camMixture_->getMixtureDiffusionCoeff(opPre);
        cptemp = camMixture_->getMolarSpecificHeat();

        for(int l=0; l<nSpc; l++)
        {
            s_mf(i,l) = mf[l];
            s_Wdot(i,l) = wdot[l]*(*spv_)[l]->MolWt();
            s_H(i,l) = htemp[l]/(*spv_)[l]->MolWt();
            s_Diff(i,l) = temp[l];
            //Specific heat capacity of species in J/Kg K
            CpSpec(i,l) =cptemp[l]/(*spv_)[l]->MolWt();
        }
        if (sootMom_.active())
        {
          for(int l=0; l<nMoments; l++)
          {
        	  moments(i,l) = mom_temp[l];
          }

          moments_dot_temp = sootMom_.rateAll(conc, mom_temp, m_T[i], opPre, 1);
          for(int l=0; l<nMoments; l++)
          {
        	  moments_dot(i,l) =  moments_dot_temp[l];
          }

          // Now get the corresponding gas phase rates and add them to s_Wdot
          wdotSootGasPhase = sootMom_.showGasPhaseRates(nSpc);
    	  for (int l=0; l< nSpc; l++)
  	      {
             s_Wdot(i,l) = s_Wdot(i,l) + wdotSootGasPhase[l] * (*spv_)[l]->MolWt();
  	      }

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
    map<string, doublereal> species;
    map<string, doublereal>::iterator sIterator;
    CamBoundary& fuelInlet = admin_.getLeftBoundary();
    species = fuelInlet.getInletSpecies();
    sIterator = species.begin();

    while(sIterator != species.end()){
        fspecies.push_back(camMech_->GetSpecies(sIterator->first));
        sIterator++;
    }

    CamBoundary& oxInlet = admin_.getRightBoundary();
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

    vector<doublereal> temp = getInletMassFrac(fuelInlet);

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

    temp = getInletMassFrac(oxInlet);
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

void FlameLet::setExternalStrainRate(const doublereal strainRate)
{
    scalarDissipationRate_.setStrainRate(strainRate);
}

void FlameLet::setExternalSDR(const doublereal sdr)
{
    scalarDissipationRate_.setSDRRate(sdr);
}

void FlameLet::setExternalTimeSDR
(
    const std::vector<doublereal>& time,
    const std::vector<doublereal>& sdr
)
{
    scalarDissipationRate_.setExternalScalarDissipationRate(time,sdr);
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

    if(nStep%10==0) cout << "Time" <<"\t" << "Residual" << endl;
    cout << x <<"\t" << res << endl;
    nStep++;

}
/*
 *output function for file output
 */
void FlameLet::reportToFile(std::string fileName, doublereal t, std::vector<double>& soln)
{

    doublereal sum;
    reporter_->openFile(fileName,false);

    reporter_->writeCustomHeader(header());

    vector<doublereal> data, axpos;
    vector<doublereal> molfrac, massfrac, temperatureVec;
    axpos = reacGeom_.getAxpos();
    int len = axpos.size();

    for(int i=0; i<len; i++)
    {
        temperatureVec.push_back(soln[i*nVar+ptrT]);
    }

    for (int i=0; i<len; i++) {

        data.clear();
        data.push_back(t);
        data.push_back(axpos[i]);
        data.push_back(scalarDissipationRate_(axpos[i],0));
        data.push_back(m_rho[i]);
        data.push_back(m_mu[i]);
        data.push_back(m_cp[i]);
        data.push_back(soln[i*nVar+ptrT]);
        if (radiation != NULL)
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

        // Add the moments to to the data output
        if (sootMom_.active())
        {
            for(int l=0; l<nMoments; l++)
            {
            	data.push_back(soln[i*nVar+ptrT+1+l]);
            }
        }


        reporter_->writeCustomFileOut(data);

    }

    reporter_->closeFile();

    //reacGeom->refine(soln,nVar,nSpc,ptrT);

}

/*
 *output file header
 */
std::vector<std::string> FlameLet::header()
{

    std::vector<std::string> headerData;

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
    if (sootMom_.active())
    {
        headerData.push_back("M0");
        headerData.push_back("M1");
        headerData.push_back("M2");
        headerData.push_back("M3");
        headerData.push_back("M4");
        headerData.push_back("M5");
    }



    return headerData;

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
