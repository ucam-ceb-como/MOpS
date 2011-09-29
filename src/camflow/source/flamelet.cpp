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
  CpSpec(mCord,nSpc),
  sootResidualZeroed(false),
  Lewis(admin_.getInputFile(),camMech_,mCord,nSpc)
{}

FlameLet::~FlameLet()
{
    if (radiation != NULL) delete radiation;
}

void FlameLet::checkSetup()
{

    if (reacGeom_.getAxpos()[mCord-1] != 1)
        throw std::invalid_argument("Mixture fraction does not go from 0 to 1. "
                                    "Check <length unit=\"m\">1.0</length>\n");

    if (stoichZ <= 0.0)
        throw std::invalid_argument("The stoichiometric mixture fraction is not"
                                    " a positive real!");

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

/*
 * Use this call method when calling from openFoam.
 * if sootResidualZeroed is TRUE then, flamelet will
 * be solved to steady state and with soot residual set to zero.
 * (i.e. no soot present at base of flame)
 *
 * When calling the Lagrangian Flamelet (i.e. dynamic)
 * do this via restart() and set steadyStateAtFlameBase to FALSE
 */
void FlameLet::solve(bool interface, bool steadyStateNoSoot)
{
	sootResidualZeroed = steadyStateNoSoot;
	solve(interface);
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
        csolve(interface);   //ank25:temp while testing splitSolve
    	//splitSolve(interface);
    }
    else
    {
        ssolve(interface);
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
        	vMom[i*nMoments] = sootMom_.getFirstMoment();
        	//cout << "vMom[i*nMoments]  " << i*nMoments <<"  "  << vMom[i*nMoments] << endl;
        	for (size_t l=1; l<nMoments; ++l)
            {
            	// ank25: Do we need to multiply by 1e6 here?
        		vMom[i*nMoments+l] = vMom[i*nMoments+l-1] + 1e6 * log(doublereal(sootMom_.getAtomsPerDiamer()));
            	//cout << "vMom[i*nMoments+l]  " << i*nMoments+l <<"  " << vMom[i*nMoments+l] << endl;
            }

        	// Change of variable of the initial moments.
        	//M=exp(M_hat-1)
        	//M_hat = ln(M+1)
       // 	for (size_t l=0; l<nMoments; ++l)
       // 	{
       // 		vMom[i*nMoments+l] = log(vMom[i*nMoments+l] + 1.0);
       // 	}

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
            {
            	// If the size of solvect_temp is not equal to nVar*cellEnd
            	// then we assume that the binary file was previously
            	// generated with soot moments switched off.  In that case
            	// load species and temperature from binary file, but don't
            	// load up moments.

            	if 	(solvect_temp.size() == nVar*cellEnd)
                {
                	std::cout << "Loading species, temperature and moments from bin file"
										<< std::endl;
            		solvect = solvect_temp;
                }
            	else
                {
                	std::cout << "Loading species and temperature from binary file" << std::endl;
                	std::cout << "Not loading moments from binary file. " << std::endl;
            		for(int i=0; i<cellEnd; i++)
                    {
                	   for(int l=0; l<nSpc; l++)
                		   solvect[i*nVar+l] = solvect_temp[i*(nVar-nMoments)+l];
                       solvect[i*nVar+ptrT] = solvect_temp[i*(nVar-nMoments)+ptrT];
                    }
                }
            }
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

        //cvw.initVectorTol(nEqn,solvect,atolVector,control_.getSpeciesRelTol(),
        //		control_.getMaxTime(),band,*this);

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
void FlameLet::restart(doublereal flameTime)
{
	// Assumption is that a restart is always a Lagrangian flamelet (i.e. not steady state)
	//steadyStateAtFlameBase = false;

	/// Change restart so that it takes in flameHeight ( or  flameTime or gridnumber?
	// Then decide set steadyStateAtFlameBase to true if above a ceetain threshold.
	// i.e. we stop calculating soot when high enough out of the flame.
	// Rename "steadyStateAtFlameBase" to "sootResidualZeroed"

	if (flameTime < Lewis.sootFlameTimeThreshold)
	{
		// Still below the time at which we stop calculating soot residual
		sootResidualZeroed = false;
	    std::cout << "Soot residual is active " << std::endl;
	}
	else
	{
		// Past the time, beyond which we no longer calculate soot.
		sootResidualZeroed = true;
	    std::cout << "Soot residual is zeroed out " << std::endl;
	}

    if (solverID == control_.CVODE) {
        CVodeWrapper cvw;
        eqn_slvd = EQN_ALL;
        int band = nVar*2;

        cvw.init(nEqn,solvect,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
            control_.getMaxTime(),band,*this,restartTime);

        // When restarting don't call cvw.solve with a global tolerance
        // as a stop criteria.  We are solving dynamic flamelets and don't
        // want to stop at steady state.
        cvw.solve(CV_ONE_STEP);
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
            string filename = "interfaceProfiles/profile"+boost::lexical_cast<std::string>(restartTime)+".dat";
            reportToFile(filename,control_.getMaxTime(), solvect);
            if (sootMom_.active())
            {
                string filenameSoot = "interfaceSootRates/sootRatesProfile"+boost::lexical_cast<std::string>(restartTime)+".dat";
                reportSootRatesToFile(filenameSoot,control_.getMaxTime(), sootComponentRatesAllCells);
            }
        //}
    }
    else if (solverID == control_.LIMEX) {
        throw std::logic_error("Error -- Limex is not yet supported");
    }

}

/*
 *segregated solver
 */
void FlameLet::ssolve
(
    bool interface
)
{

    int seg_eqn, band;
    vector<doublereal> seg_soln_vec;

    if ( solverID == control_.CVODE){

       CVodeWrapper cvw;

       for (int i=0; i<control_.getNumIterations();i++){

           /*
            *solve soot moment equations
            */
           if (sootMom_.active())
           {
             cout << "Solving moment equations  " << i << endl;
             //int dd; cin >> dd;
             eqn_slvd = EQN_MOMENTS;
             seg_eqn = nMoments*mCord;
             band = nMoments*2;
             extractSootMoments(seg_soln_vec);
             // Might need to change tolerances for moments
             cvw.init(seg_eqn,seg_soln_vec,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
                 control_.getMaxTime(),band,*this,0.0);
             cvw.solve(CV_ONE_STEP,control_.getResTol());
             mergeSootMoments(&seg_soln_vec[0]);
           }

        /*
         *solve species equations
         */
        cout << "Solving species equations  " << i << endl;
        //int dd; cin >> dd;
        eqn_slvd = EQN_SPECIES;
        seg_eqn = nSpc*mCord;
        band = nSpc*2;
        extractSpeciesVector(seg_soln_vec);
        cvw.init(seg_eqn,seg_soln_vec,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
            control_.getMaxTime(),band,*this,0.0);
        cvw.solve(CV_ONE_STEP,control_.getResTol());
        mergeSpeciesVector(&seg_soln_vec[0]);

        /*
         *solve energy equation
         */
        cout << "Solving energy equation  " << i << endl;
        //cin >> dd;
        eqn_slvd = EQN_ENERGY;
        seg_eqn = mCord;
        band = 1;
        extractEnergyVector(seg_soln_vec);
        cvw.init(seg_eqn,seg_soln_vec,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
            control_.getMaxTime(),band,*this,0.0);
        cvw.solve(CV_ONE_STEP,control_.getResTol());
        mergeEnergyVector(&seg_soln_vec[0]);

      }

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

    else if (solverID == control_.LIMEX) {
        throw std::logic_error("Error -- Limex is not yet supported");
    }
  }
//----------------------

/*
 *splitting solver
 */
void FlameLet::splitSolve
(
    bool interface
)
{

    int seg_eqn, band;
    vector<doublereal> seg_soln_vec;

    // Get the time to which we are integrating.
    doublereal nextTime = control_.getMaxTime();

    // Break this time step into Nsteps small steps
    // todo: Get NSteps via function call and ultimately from xml file or similar
    // Note: The assumption here is that the current (or initial time) is t=0.
    // This seems to be always be the case when solving flamlets.
    // This is even true when calling via the interface as the integration occurs
    // over tau.

    // NOT TRUE:  Flamelet restart has a restart time.  Fix this !!

    int Nsteps = 100;
    doublereal deltaTime = nextTime / (doublereal)Nsteps;
    doublereal intermediateTime = 0.0 ;

    if ( solverID == control_.CVODE){

       CVodeWrapper cvw;

       //for (int i=0; i<control_.getNumIterations();i++){
       for (int i=0; i<Nsteps;i++){

    	// Set the integration time
    	nextTime = intermediateTime + deltaTime;


        // solve species equations
        cout << "Solving species equations  " << i << endl;
        //int dd; cin >> dd;
        eqn_slvd = EQN_SPECIES;
        seg_eqn = nSpc*mCord;
        band = nSpc*2;
        extractSpeciesVector(seg_soln_vec);

        cvw.init(seg_eqn,seg_soln_vec,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
        		nextTime,band,*this,intermediateTime);

        cvw.solve(CV_ONE_STEP); //,control_.getResTol());

        mergeSpeciesVector(&seg_soln_vec[0]);


         //solve energy equation
        cout << "Solving energy equation  " << i << endl;
        //cin >> dd;
        eqn_slvd = EQN_ENERGY;
        seg_eqn = mCord;
        band = 1;
        extractEnergyVector(seg_soln_vec);

        cvw.init(seg_eqn,seg_soln_vec,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
        		nextTime,band,*this,intermediateTime);

        cvw.solve(CV_ONE_STEP); //,control_.getResTol());

        mergeEnergyVector(&seg_soln_vec[0]);
/*

        // solve combined species and energy equations
        cout << "Solving species and energy equations  " << i << endl;
        //int dd; cin >> dd;
        eqn_slvd = EQN_SPECIES_ENERGY;
        seg_eqn = (nSpc+1)*mCord;
        band = (nSpc+1)*2;
        extractSpeciesAndEnergyVector(seg_soln_vec);

        cvw.init(seg_eqn,seg_soln_vec,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
        		nextTime,band,*this,intermediateTime);

        cvw.solve(CV_ONE_STEP,control_.getResTol());

        mergeSpeciesAndEnergyVector(&seg_soln_vec[0]);

*/
        /*
         *solve soot moment equations
         */
        if (sootMom_.active())
        {
          cout << "Solving moment equations  " << i << endl;
          //int dd; cin >> dd;
          eqn_slvd = EQN_MOMENTS;
          seg_eqn = nMoments*mCord;
          band = nMoments*2;
          extractSootMoments(seg_soln_vec);
          // Might need to change tolerances for moments
          cvw.init(seg_eqn,seg_soln_vec,control_.getSpeciesAbsTol(),control_.getSpeciesRelTol(),
        		  nextTime,band,*this,intermediateTime);

          cvw.solve(CV_ONE_STEP); //,control_.getResTol());

          mergeSootMoments(&seg_soln_vec[0]);
        }

    	// Increment time
    	intermediateTime = intermediateTime + deltaTime;

      }

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
    else if (solverID == control_.LIMEX) {
        throw std::logic_error("Error -- Limex is not yet supported");
    }
  }

//----------------------

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
        	if (sootResidualZeroed)
        	{
        		sootMomentResidualAtFlameBase(t,y,&resMom[0]);
        	}
        	else
        	{
        		sootMomentResidual(t,y,&resMom[0]);
        	}
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

            if (sootMom_.active())
            {
                // Put moments into master residual
                for (int l=0; l<nMoments; ++l)
                {
                    f[i*nVar+ptrT+1+l] = resMom[i*nMoments+l];
                }
            }

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
        else if (eqn_slvd==EQN_SPECIES_ENERGY)
        {
        	mergeSpeciesAndEnergyVector(y);
            saveMixtureProp(&solvect[0]);
            speciesResidual(t,y,&resSp[0]);
            energyResidual(t,y,&resT[0]);
            for (int i=0; i<mCord; ++i)
            {
                // Put species residual into master residual
            	for (int l=0; l<nSpc; ++l)
                {
                    f[i*nSpc+l] = resSp[i*nSpc+l];
                }
            	// Put temperature into master residual
                f[i*nSpc+ptrT] = resT[i];
            }
        }

        else if (eqn_slvd==EQN_MOMENTS)
        {
        	mergeSootMoments(y);
            saveMixtureProp(&solvect[0]);
        	if (sootResidualZeroed)
        	{
        		sootMomentResidualAtFlameBase(t,y,&resMom[0]);
        	}
        	else
        	{
        		sootMomentResidual(t,y,&resMom[0]);
        	}
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
    if (sootMom_.active())
    {
        for (int i=0; i<nMoments*mCord; ++i)
        {
            resNorm += resMom[i]*resMom[i];
        }
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
    doublereal sdr, sdrPE, sdrPW;
    doublereal source;
    doublereal deltax = 0;

    /*
     *starting with mixture fraction zero: i.e oxidizer
     *inlet. The fuel composition is zero. Left inlet is
     *considered as the oxidizer inlet and the concentrations
     *are held constant
     */

    for (int l=0; l<nSpc; ++l)
    {
        f[l] = 0.0;
    }

    /*for (int l=0; l<nSpc; ++l)
    {
        // Computation at the i = 0 end of the grid
        convection(0,l) = oneby16*(1.0/Le(0,l)-1)*4*(s_mf(1,l)-s_mf(0,l))*(m_rho[1]-m_rho[0]);

        // Computation at the i = mCord - 1 end of the grid
        convection(mCord-1,l) = oneby16*(1.0/Le(mCord-1,l)-1)*4*(s_mf(mCord-1,l)-s_mf(mCord-2,l))*(m_rho[mCord-1]-m_rho[mCord-2]);
    }*/

    /*
     *interior mixture fraction coordinates
     */
    for (int i=iMesh_s; i<iMesh_e; ++i)
    {

        zPE = 0.5*(dz[i]+dz[i+1]);
        zPW = 0.5*(dz[i]+dz[i-1]);
        deltax = zPE + zPW;

        sdr = scalarDissipationRate_(reacGeom_.getAxpos()[i],t);

        doublereal diffusionConstant = sdr/(2.0*dz[i]);
        for (int l=0; l<nSpc; ++l)
        {
            grad_e = (s_mf(i+1,l)-s_mf(i,l))/zPE;
            grad_w = (s_mf(i,l)-s_mf(i-1,l))/zPW;
            source = s_Wdot(i,l)/m_rho[i];
            f[i*nSpc+l] = diffusionConstant*(grad_e-grad_w)/Lewis(i,l)
                          + source;
        }

        if (   admin_.getFlameletEquationType() == admin_.COMPLETE
            && Lewis.type() != LewisNumber::UNITY)
        {
            sdrPE = scalarDissipationRate_(reacGeom_.getAxpos()[i+1],t);
            sdrPW = scalarDissipationRate_(reacGeom_.getAxpos()[i-1],t);

            doublereal convectionConstant
                       = 0.25/m_rho[i]
                         *(
                            (m_rho[i+1]*sdrPE - m_rho[i-1]*sdrPW)/deltax
                           +(m_rho[i]*sdr*m_cp[i]/m_k[i])
                            *(m_k[i+1]/m_cp[i+1] - m_k[i-1]/m_cp[i-1])/deltax
                         );

            for (int l=0; l<nSpc; ++l)
            {
                f[i*nSpc+l] +=
                  convectionConstant
                *((1.0/Lewis(i,l))-1.0)
                *(s_mf(i+1,l)-s_mf(i-1,l))/deltax;
            }
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
    doublereal sdr, sdrPE, sdrPW;

    for (int l=0; l<nMoments; ++l)
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

        // Get the scalar dissipation rate for the cell.
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

/*
        // Code below is for soot flamelet equations after change of variable
    	// M=exp(M_hat -1)
    	// M_hat = ln(M+1)
        doublereal diffusionConstant = sdr/(2.0);
        for (int l=0; l<nMoments; ++l)
        {
            grad_e = (moments(i+1,l)-moments(i,l))/zPE;
            grad_w = (moments(i,l)-moments(i-1,l))/zPW;
            source = moments_dot(i,l)/ ( m_rho[i] * exp(moments(i,l)) ) ;
            f[i*nMoments+l] = diffusionConstant
            		        * ( (grad_e-grad_w)/dz[i] +
            		        	pow((grad_e+grad_w)/2.0,2.0) )
                            + source;
       }
*/
    }

    /*
     *Mixture fraction 1. The fuel inlet. Moments held constant here.
     */

    for (int l=0; l<nMoments; ++l)
    {
        f[iMesh_e*nMoments+l] = 0.0;
    }

}

/*
 * If we are solving a flamelet at the base of a flame then
 * set the soot residual to zero.
 * i.e. we solve a steady state flamelet with no soot present.
 */

void FlameLet::sootMomentResidualAtFlameBase
(
    const doublereal& t,
    doublereal* y,
    doublereal* f
)
{
    for (int i=iMesh_s-1; i<iMesh_e+1; ++i)
    {
        for (int l=0; l<nMoments; ++l)
        {
        	f[i*nMoments+l] = 0.0;
        }
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
    doublereal deltax=0;
    doublereal sdr, sdrPE, sdrPW;

    /*
     *starting with mixture fraction zero: i.e oxidizer
     *inlet. The temperature is fixed at the oxidizer
     *inlet temperature
     */
    f[0] = 0.0;

    /*
     *intermediate mixture fraction coordinates
     */
    for (int i=iMesh_s; i<iMesh_e; ++i)
    {

        zPE = 0.5*(dz[i]+dz[i+1]);
        zPW = 0.5*(dz[i]+dz[i-1]);
        deltax = zPE + zPW;

        source = 0.0;
        for (int l=0; l<nSpc; ++l)
        {
            source += s_Wdot(i,l)*s_H(i,l);
        }

        // Get the scalar dissipation rate.
        sdr = scalarDissipationRate_(reacGeom_.getAxpos()[i],t);

        grad_e = (m_T[i+1]-m_T[i])/zPE;
        grad_w = (m_T[i]-m_T[i-1])/zPW;

        f[i] = 0.5*sdr*(grad_e-grad_w)/dz[i]
               - source/(m_rho[i]*m_cp[i]);

        /**
         * Add some extra terms so that we agree with FlameMaster?
         */
        if (admin_.getFlameletEquationType() == admin_.COMPLETE)
        {
            doublereal tGrad = (m_T[i+1]-m_T[i-1])/deltax;
            doublereal cpGrad = (m_cp[i+1]-m_cp[i-1])/deltax;
            doublereal sumYGrad = 0.0;
            for (int l=0; l<nSpc; ++l)
            {
                sumYGrad += (1.0/Lewis(i,l)) * CpSpec(i,l) * (s_mf(i+1,l)-s_mf(i-1,l))/deltax;
            }
            f[i] += (sdr/(2.0*m_cp[i])) * tGrad * (cpGrad + sumYGrad);
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
               sootVolumeFractionMaster[i]
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
    vector<doublereal> exp_mom_temp(nMoments,0.0);
    vector<doublereal> conc(nSpc,0.0);
    vector<doublereal> wdotSootGasPhase(nMoments,0.0);
    vector<doublereal> sootComponentRatesTemp(nMoments*4,0.0);

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
        if (Lewis.type() == LewisNumber::CALCULATED) temp = camMixture_->getMixtureDiffusionCoeff(opPre);
        cptemp = camMixture_->getMolarSpecificHeat();

        for(int l=0; l<nSpc; l++)
        {
            s_mf(i,l) = mf[l];
            s_Wdot(i,l) = wdot[l]*(*spv_)[l]->MolWt();
            s_H(i,l) = htemp[l]/(*spv_)[l]->MolWt();
            //Specific heat capacity of species in J/Kg K
            CpSpec(i,l) =cptemp[l]/(*spv_)[l]->MolWt();
            if (Lewis.type() == LewisNumber::CALCULATED)
            {
                s_Diff(i,l) = temp[l];
                Lewis.calcLewis(i,l) = m_k[i]/(m_rho[i]*m_cp[i]*temp[l]);
            }
        }

        if (sootMom_.active())
        {
          for(int l=0; l<nMoments; l++)
          {
        	  moments(i,l) = mom_temp[l];			// ank25:  This looks correct
        	//  exp_mom_temp[l] = exp(mom_temp[l])-1.0;
          //}

          // Change of variable before calling moment rates:
      	  // M=exp(M_hat -1)
      	  // M_hat = ln(M+1)

          // DEBUG:  Before calling rateAll check if any moments have gone negative.
          // If so then:
          // a) Output the cell number
          // b) Output the moment
          if (moments(i,l) <0.0)
          {
        	  std::cout << "Negative moment found before calling ratesAll" << std::endl;
        	  std::cout << "Cell index : " << i << std::endl;
        	  std::cout << "Moment index : " << l << std::endl;
        	  std::cout << "Moment value : " << moments(i,l) << std::endl;
          }
		  // end DEBUG

          }


          moments_dot_temp = sootMom_.rateAll(conc, mom_temp, m_T[i], opPre, 1);
          for(int l=0; l<nMoments; l++)
          {
        	  moments_dot(i,l) =  moments_dot_temp[l];
          }

          // Now get the corresponding gas phase rates and add them to s_Wdot
          // Only do this if we are solving a Lagrangian flamelet (not steady state)
          if (sootResidualZeroed == false)
          {
            wdotSootGasPhase = sootMom_.showGasPhaseRates(nSpc);
    	    for (int l=0; l< nSpc; l++)
  	        {
               s_Wdot(i,l) = s_Wdot(i,l) + wdotSootGasPhase[l] * (*spv_)[l]->MolWt();
  	        }

    	    // Also get the component soot rates for output.
            // Ideally we would only do this at the major output times rather than at each call
    	    // of saveMixtureProp. (Inefficient)
    	    // ToDo: Move this it a better place.
    	    sootComponentRatesTemp = sootMom_.showSootComponentRates(nMoments);
    	    for (int l=0; l< nMoments*4; l++)
  	        {
               sootComponentRatesAllCells(i,l) = sootComponentRatesTemp[l];
  	        }

    	    // Calculate soot properties at each Z point.
    	    // We need volume fraction for radiation.
    	    avgSootDiamMaster[i] = sootMom_.avgSootDiam();
    	    dispersionMaster[i] = sootMom_.dispersion();
    	    sootSurfaceAreaMaster[i] = sootMom_.sootSurfaceArea(moments(i,0));
    	    sootVolumeFractionMaster[i] = sootMom_.sootVolumeFraction(moments(i,0));
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

        // Add the moments and soot properties to the data output
        if (sootMom_.active())
        {
            data.push_back(avgSootDiamMaster[i]);
            data.push_back(dispersionMaster[i]);
            data.push_back(sootSurfaceAreaMaster[i]);
            data.push_back(sootVolumeFractionMaster[i]);
        	for(int l=0; l<nMoments; l++)
            {
            	data.push_back(soln[i*nVar+ptrT+1+l]);
            }
        }
        reporter_->writeCustomFileOut(data);
    }

    reporter_->closeFile();

    // Output Lewis Numbers to File.
    std::ofstream file("LewisNumbers");
    file << setw(5) << "Z" << " ";
    for (int l=0; l<nSpc; l++)
    {
        file << setw(8) << (*spv_)[l]->Name() << " ";
    }
    file<<std::endl;
    for (int i=0; i<mCord; ++i)
    {
        file << setw(5) << axpos[i] << " ";
        for (int l=0; l<nSpc; l++)
        {
            file << setprecision(5) << setw(8) <<  Lewis(i,l) << " ";
        }
        file << std::endl;
    }
    file.close();

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
        headerData.push_back("SootAvDiam");
        headerData.push_back("SootDisp");
        headerData.push_back("SootArea");
        headerData.push_back("SootVolFrac");
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


/*
 *output function for file output
 */
void FlameLet::reportSootRatesToFile(std::string fileName, doublereal t, Array2D& rates)
{

    doublereal sum;

    std::cout << "Reporting Component Soot Rates to File" << std::endl;

    reporter_->openFile(fileName,false);

    reporter_->writeCustomHeader(sootRatesHeader());

    vector<doublereal> data, axpos;

    //vector<doublereal> molfrac, massfrac, temperatureVec;
    axpos = reacGeom_.getAxpos();
    int len = axpos.size();

    for (int i=0; i<len; i++) {

        data.clear();
        data.push_back(t);
        data.push_back(axpos[i]);
        data.push_back(scalarDissipationRate_(axpos[i],0));

        for(int l=0; l<nMoments*4; l++){
            data.push_back(rates(i,l));
        }

        reporter_->writeCustomFileOut(data);

    }
    reporter_->closeFile();
}

/*
 *output file header
 */
std::vector<std::string> FlameLet::sootRatesHeader()
{
    std::vector<std::string> headerData;

    headerData.clear();
    headerData.push_back("time");
    headerData.push_back("Z");
    headerData.push_back("SDR");

    for(int l=0; l<nMoments; l++){
    	headerData.push_back("Nuc_M"+boost::lexical_cast<std::string>(l));
    }
    for(int l=0; l<nMoments; l++){
    	headerData.push_back("Coag_M"+boost::lexical_cast<std::string>(l));
    }
    for(int l=0; l<nMoments; l++){
    	headerData.push_back("Cond_M"+boost::lexical_cast<std::string>(l));
    }
    for(int l=0; l<nMoments; l++){
    	headerData.push_back("Surf_M"+boost::lexical_cast<std::string>(l));
    }
    return headerData;
}
