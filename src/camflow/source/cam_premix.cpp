/*
 * File:   cam_premix.cpp
 * Author: vinod
 *
 * File purpose:
 *  This class implements the premix model
 *
 * Licence:
 *  This file is part of "Camflow".
 *
 *  Camflow is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * Contact:
 *  Dr Markus Kraft
 *  Dept of Chemical Engineering
 *  University of Cambridge
 *  New Museum Site
 *  Pembroke Street
 *  Cambridge
 *  CB2 3RA
 *  UK
 *
 *  Email   :   mk306@cam.ac.uk
 *  Website
 *
 * Created on January 18, 2009, 7:57 PM
 */


#include <vector>
#include <map>
#include <math.h>
#include <cstring>
#include <stdexcept>

#include "array.h"
#include "cam_soot.h"
#include <iostream>
#include <sstream>
#include "cam_params.h"
#include "cam_admin.h"
#include "cam_profile.h"
#include "cam_residual.h"
#include "cam_setup.h"
#include "cam_premix.h"
#include "cam_reporter.h"
#include "cvode_wrapper.h"
#include "limex_wrapper.h"
//#include "ida_wrapper.h"


using namespace Camflow;
using namespace std;
using namespace Gadgets;


CamPremix::CamPremix
(
    CamAdmin& ca,
    CamConfiguration& config,
    CamControl& control_,
    CamGeometry& cg,
    CamProfile& cp,
    CamSoot& cs,
    Mechanism& mech
)
:
    CamSetup(ca, config, control_, cg, cp, cs, mech)
{}


double CamPremix::dydx(double nr1, double nr2, double dr) const
{
    return (nr1 - nr2)/dr;
}


double CamPremix::dydx(double nr1, double nr2, double nr3, double dr) const
{
    return (nr1 - 2.0*nr2 + nr3)/(dr*dr);
}


int CamPremix::eval(double t, double* y, double* ydot, bool jacEval)
{
    /*
     *this is called by the DAE wrapper object. Given y, ydot is returned
     */
    residual(t, y, ydot);

    return 0;
}

void CamPremix::speciesResidual
(
    const double& time, 
    double* y,
    double* f
)
{
    /*
     *prepare flux terms
     */
    double convection, diffusion, source;

    for (int i= iMesh_s; i< iMesh_e; i++ )
    {
        for (int l=0; l<nSpc; l++)
        {
            /*
             *convection term
             */
            convection = -m_u[i]*dydx(s_mf(i,l),s_mf(i-1,l),dz[i]);
            //convection = -m_u[i]*dydx(s_mf(i+1,l),s_mf(i,l),s_mf(i-1,l),dz[i]);
            /*
             *diffusion term
             */
            diffusion = dydx(s_jk(i,l),s_jk(i+1,l),dz[i])/m_rho[i];

            source = s_Wdot(i,l)   *(*spv_)[l]->MolWt()/m_rho[i];
            f[i*nSpc+l] = convection + diffusion + source ;

        }
    }
}


void CamPremix::energyResidual
(
    const double& time,
    double* y,
    double* f
)
{
    bool flxTrans = false;

    if (admin_.getEnergyModel() == admin_.ADIABATIC)
    {
        double sumHSource, flxk, sumFlx;;
        for (int i=iMesh_s; i<iMesh_e; i++)
        {
            /*
             *evaluate the summation terms
             */
            sumHSource = 0.0;
            sumFlx = 0.0;
            for (int l=0; l<nSpc; l++)
            {
                sumHSource += s_Wdot(i,l)*s_H(i,l);
                flxk = 0.5*(s_jk(i-1,l)+s_jk(i,l));
                sumFlx += flxk*s_cp(i,l)/(*spv_)[l]->MolWt();
            }

            /*
             *convection
             */
            double convection =  -m_u[i]*dydx(m_T[i],m_T[i-1],dz[i]);
            //double convection = -m_u[i]*dydx(m_T[i+1], m_T[i], m_T[i-1], dz[i]);
            /*
             *conduction
             */
            double conduction = dydx(m_q[i+1],m_q[i],dz[i]);
            /*
             *flux transport
             */
            sumFlx *= dydx(m_T[i],m_T[i-1],dz[i]);

            /*
             *residual
             */
            if (flxTrans)
            {
                f[i] = convection + (conduction - sumHSource - sumFlx)/(m_rho[i]*m_cp[i]);
            }
            else
            {
                f[i] = convection + (conduction - sumHSource )/(m_rho[i]*m_cp[i]);
            }
        }
    }
    else
    {
        for (int i=iMesh_s; i< iMesh_e; i++)
        {
            f[i] = 0;
            //f[ptrTP] = y[ptrTW]-y[ptrTP];
        }
    }
}


/**
 *@parameter[in]    y       current solution value
 *@parameter[out]   f       residuals
 */
void CamPremix::residual(const double& t, double* y, double* f)
{
    resSp.resize(cellEnd*nSpc,0);
    resT.resize(cellEnd,0);
    resFlow.resize(cellEnd,0);

    /*
     *residual evaluation function for premix. Internal
     *cell residuals are evaluated by calling the base class
     *function
     */
    if (eqn_slvd == EQN_ALL)
    {
        if (admin_.getEnergyModel() == admin_.ADIABATIC)
        {
            saveMixtureProp(t,y,true,false);
            updateThermo();
        } else
        {
            saveMixtureProp(t,y,false,false);
        }

        saveFlowVariables(y);
        updateDiffusionFluxes();

        //massFlowBoundary(t,y,&resM[0]);
        speciesBoundary(t,y,&resSp[0]);
        energyBoundary(t,y,&resT[0]);

        speciesResidual(t,y,&resSp[0]);
        energyResidual(t,y,&resT[0]);
        //massFlowResidual(t,y,&resM[0]);

        for (int i=0; i<cellEnd; i++)
        {
            for (int l=0; l<nSpc; l++)
            {
                f[i*nVar+l] = resSp[i*nSpc+l];
            }

            f[i*nVar+ptrF] = resFlow[i];
            f[i*nVar+ptrT] = resT[i];
        }
    }
    else
    {
        if (eqn_slvd == EQN_SPECIES)
        {
            mergeSpeciesVector(y);
            saveMixtureProp(t,&solvect[0],false,false);
            saveFlowVariables(&solvect[0]);
            updateDiffusionFluxes();
            speciesBoundary(t,y,f);
            speciesResidual(t,y,f);
        }
        else if (eqn_slvd == EQN_CONTINUITY)
        {
            mergeContinuity(y);
            saveMixtureProp(t,&solvect[0],false,false);
            saveFlowVariables(&solvect[0]);
            massFlowBoundary(t,y,f);
            massFlowResidual(t,y,f);
        }
        else if (eqn_slvd == EQN_ENERGY)
        {
            mergeEnergyVector(y);
            saveMixtureProp(t,&solvect[0],true,false);
            saveFlowVariables(&solvect[0]);
            updateDiffusionFluxes();
            updateThermo();
            energyResidual(t,y,f);
            energyBoundary(t,y,f);
        }
    }
}


void CamPremix::massFlowResidual
(
    const double& time,
    double* y,
    double* f
)
{
    for (int i=iMesh_s; i<iMesh_e; i++)
    {
        f[i] = -m_u[i]*(m_flow[i] - m_flow[i-1])/dz[i];
        //f[ptrCP] = y[ptrCW]-y[ptrCP];
    }
}


void CamPremix::massFlowBoundary
(
    const double& t,
    double* y,
    double* f)
{
    //f[0] = ud_inlet.FlowRate - y[0]
    //f[iMesh_e] = y[iMesh_e-1] - y[iMesh_e];

    f[0] = ud_inlet.Vel*(ud_inlet.FlowRate - m_flow[0])/dz[0];

    /*
     *Mass flow
     */

    f[iMesh_e]
        = -m_u[iMesh_e]*(m_flow[iMesh_e] - m_flow[iMesh_e-1])/dz[iMesh_e];
}


void CamPremix::speciesBoundary
(
    const double& t,
    double* y,
    double* f
)
{
    //-------------------------------------------
    //
    //  Left Boundary Settings
    //
    //------------------------------------------

    double convection, diffusion;

    for (int l=0; l<nSpc ; l++)
    {
        convection = m_u[0]*dydx(ud_inlet.Species[l],y[l],dz[0]);
        diffusion = dydx(0,s_jk(iMesh_s,l),dz[0])/m_rho[0];
        f[l] = convection + diffusion;
    }

    // remaining nSpc * num Cells entries are not touched by this loop

     //---------------------------------------------------
     //
     // Right Boundary Settings
     //
     //---------------------------------------------------

    for (int l=0; l<nSpc; l++)
    {
        // Only deals with the final cell (offset iMesh_e * nSpc
        f[iMesh_e*nSpc+l]
            = -m_u[iMesh_e]*dydx(s_mf(iMesh_e,l),s_mf(iMesh_e-1,l),dz[iMesh_e]);
    }
}


/*
 *soot moment boundary conditions
 */
void CamPremix::momentBoundary
(
    const double& t,
    double* y,
    double* f
)
{
    //-------------------------------------------
    //
    //  Left Boundary Settings
    //
    //------------------------------------------
    for (int l=0; l<nMoments; l++)
    {
        f[l] = 0;
    }

    //-------------------------------------------
    //
    //  Right Boundary Settings
    //
    //------------------------------------------
    double convection;
    for (int l=0; l<nMoments; l++)
    {
        double phi_e = y[iMesh_e*nVar+nSpc+l];
        double phi_w = y[(iMesh_e-1)*nVar+nSpc+l];
        convection = -m_u[iMesh_e]*(phi_e-phi_w)/dz[iMesh_e];
        f[iMesh_e*nMoments+l] = convection;
    }

}

void CamPremix::energyBoundary
(
    const double& t,
    double* y,
    double* f
)
{
    //-------------------------------------------
    //
    //  Left Boundary Settings
    //
    //------------------------------------------
    f[0] = 0.0;

    //---------------------------------------------------
    //
    // Right Boundary Settings
    //
    //---------------------------------------------------
    if (admin_.getEnergyModel() == admin_.ADIABATIC)
    {
        f[iMesh_e]
            = -m_u[iMesh_e]*(dydx(m_T[iMesh_e],m_T[iMesh_e-1],dz[iMesh_e]));
    }
    else
    {
        f[iMesh_e] = 0;
    }
}


void CamPremix::solve()
{
    /*
     *function to set up the solver.
     *Initialisation of the solution vector is
     *done here
     */

    CamBoundary cb;

    //storeObjects(control_,ca,cg,cp,cb,mech);

    // Hack!
    set_nVar(nSpc + 2);
    set_ptrF(nSpc + 1);

    //nVar = nSpc + 2;
    //ptrF = nSpc;
    //ptrT = ptrF + 1;

//    nEqn = nVar*reacGeom_.getnCells();

    /*
     *this loop is for residual evaluation by the base class
     */

    initSolutionVector(cb);

    reporter_->header("PREMIX");
    reporter_->problemDescription(cb, *this);
    header();

    /*
     *solution call
     */
    if (control_.getSolutionMode() == control_.COUPLED)
    {
        csolve();
    }
    else
    {
        ssolve();
        csolve();
    }
}

/**
 *  Save the inlet after integration
 */
void CamPremix::saveInlet()
{
    soln_inlet.Species.clear();

    for (int l=0; l<nSpc; l++)
    {
        soln_inlet.Species.push_back(s_mf(0,l));
    }

    soln_inlet.Dens = m_rho[0];
    soln_inlet.FlowRate = m_flow[0];
    soln_inlet.T = m_T[0];
    soln_inlet.Vel = m_u[0];
}


/*
 *solve the problem  for continuation calls from the external interface
 */
void CamPremix::solve
(
    vector<Thermo::Mixture>& cstrs,
    const vector< vector<double> >& iniSource,
    const vector< vector<double> >& fnlSource,
    Mechanism& mech,
    CamControl& control_,
    CamAdmin& ca,
    CamGeometry& cg,
    CamProfile& cp
)
{
    CamBoundary cb;

    //storeObjects(control_,ca,cg,cp,cb,mech);
    set_nVar(nSpc + 2);
    set_ptrF(nSpc + 1);
    //nVar = nSpc + 2;
    //ptrF = nSpc;
    //ptrT = ptrF + 1;

 //   reacGeom_.addZeroWidthCells();

    //nEqn = nVar *(reacGeom_.getnCells());

    /*
     *this loop is for residual evaluation by the base class
     */

    /*
     *set the source terms for the particle process
     */
    setParticleSource(iniSource,fnlSource);

    /*
     *initialise the solution vector. The cstrs contain only the interior
     *cells. The inlet concentration needs to be taken care of
     */
    solvect.clear();
    solvect.resize(nEqn,0.0);

    /*
     *-----------------------------------------
     *inlet cell
     *-----------------------------------------//
     */

    //species
    for (int l=0; l<nSpc; l++)
    {
        solvect[l] = soln_inlet.Species[l];
    }
    //Temperature
    solvect[ptrT] = ud_inlet.T;
    //mass flow
    solvect[ptrF] = ud_inlet.FlowRate;

    /*
     *interior cells
     */
    int len = cstrs.size();
    for (int i=iMesh_s; i<iMesh_e; i++)
    {
        vector<double> massFrac;
        cstrs[i-1].GetMassFractions(massFrac);

        for (int l=0; l<nSpc; l++)
        {
            solvect[i*nVar+l]=massFrac[l];
        }

        //Temperature
        solvect[i*nVar+ptrT] = cstrs[i-1].Temperature();
        //Mass flow
        solvect[i*nVar+ptrF] = cstrs[i-1].MassDensity()*cstrs[i-1].Velocity();
    }

    /*-------------------------------------
     *last cell
     *
     *-----------------------------------//
     */

    //species
    for (int l=0; l<nSpc; l++)
    {
        solvect[iMesh_e*nVar+l] = solvect[(iMesh_e-1)*nVar+l];
    }

    //Temperature
    solvect[iMesh_e*nVar+ptrT] = solvect[(iMesh_e-1)*nVar+ptrT];
    //Mass flow
    solvect[iMesh_e*nVar+ptrF] = solvect[(iMesh_e-1)*nVar+ptrF];


    /*
     *solve the problem
     */
    if (control_.getSolutionMode() == control_.COUPLED)
    {
        csolve();
    }
    else
    {
        ssolve();
        csolve();
    }


    /*
     *set the solution back into the vector of mixtures
     */
    for (int i=iMesh_s; i<iMesh_e; i++)
    {
        vector<double> massfracs;

        for (int l=0; l<nSpc; l++)
        {
            massfracs.push_back(solvect[i*nVar+l]);
        }

        cstrs[i-1].SetMassFracs(massfracs);
        cstrs[i-1].SetTemperature(m_T[i]);
        cstrs[i-1].SetMassDensity(m_rho[i]);
        cstrs[i-1].SetVelocity(m_u[i]);
    }

}

void CamPremix::initSolutionVector(CamBoundary& cb)
{
    //dz = reacGeom->getGeometry();

    cb = admin_.getLeftBoundary();

    storeInlet(cb, ud_inlet);

    /*
     *store the inlet species mass fraction, and flow rate for
     *use in the boundary residual definitions
     */

    //createSolnVector(cb,control_,solvect);
    solvect.resize(nEqn, 0);

    vector<double> vMass = initMassFlow(cb);
    vector<double> vSpec = initSpecies(cb);
    vector<double> vT = initTemperature(cb);

    cout<< "ptrF " << ptrF << " ptrT " << ptrT << endl;

    mergeEnergyVector(&vT[0]);
    mergeContinuity(&vMass[0]);
    mergeSpeciesVector(&vSpec[0]);
}


/*
 *coupled solver
 */
void CamPremix::csolve()
{
    cout<< "\nPremix coupled solver\n";
    eqn_slvd = EQN_ALL;
    int solver = control_.getSolver();
    int band = nVar*2;

    if (solver == control_.CVODE)
    {
        CVodeWrapper cvw;

        cvw.init
        (
            nEqn,
            solvect,
            control_.getSpeciesAbsTol(),
            control_.getSpeciesRelTol(),
            control_.getMaxTime(),
            band,
            *this
        );

        if (control_.getResidualMonitor())
        {
            cvw.solve(CV_ONE_STEP, control_.getResTol());
        }
        else
        {
            resNorm = cvw.solve(CV_ONE_STEP);
            cout << "residual at the end of maxTime " << resNorm << endl;
        }

        reportToFile(control_.getMaxTime(), &solvect[0]);
    }
    else if (solver == control_.LIMEX)
    {
        throw std::logic_error("Error -- Limex is not yet supported");
    }
}


/*
 *segregated solver
 */
void CamPremix::ssolve()
{
    /*
     *preparations
     */
    int seg_eqn, band;
    vector<double> seg_soln_vec;

    int solverID = control_.getSolver();

    if (solverID == control_.CVODE)
    {
        CVodeWrapper cvw;
        for (int i=0; i<control_.getNumIterations(); i++)
        {
            /*
             *integrate mass flow
             */
            cout << "solving mass flow: " << i << endl;
            eqn_slvd = EQN_CONTINUITY;
            band = 1;
            extractContinuity(seg_soln_vec);

            cvw.init
            (
                mCord,
                seg_soln_vec,
                control_.getSpeciesAbsTol(),
                control_.getSpeciesRelTol(),
                control_.getMaxTime(),
                band,
                *this
            );

            cvw.solve(CV_ONE_STEP,1e-03);
            mergeContinuity(&seg_soln_vec[0]);

            /*
             *integrate species
             */
            cout << "Solving species: " << i << endl;
            eqn_slvd = EQN_SPECIES;
            seg_eqn = nSpc*mCord;
            band = nSpc*2;
            extractSpeciesVector(seg_soln_vec);
            
            

            cvw.init
            (
                seg_eqn,
                seg_soln_vec,
                control_.getSpeciesAbsTol(),
                control_.getSpeciesRelTol(),
                control_.getMaxTime(),
                band,
                *this
            );

            cvw.solve(CV_ONE_STEP,1e-03);
            
            
            
            mergeSpeciesVector(&seg_soln_vec[0]);

            if (admin_.getEnergyModel() == admin_.ADIABATIC)
            {
                /*
                 *integrate energy
                 */
                cout << "Solving energy: " << i << endl;
                eqn_slvd = EQN_ENERGY;
                band = 1;
                extractEnergyVector(seg_soln_vec);

                cvw.init
                (
                    mCord,
                    seg_soln_vec,
                    control_.getSpeciesAbsTol(),
                    control_.getSpeciesRelTol(),
                    control_.getMaxTime(),
                    band,
                    *this
                );

                cvw.solve(CV_ONE_STEP, 1e-03);
                mergeEnergyVector(&seg_soln_vec[0]);
            }
        }
    }
    else if (solverID == control_.LIMEX)
    {
        throw std::logic_error("Error -- Limex is not yet supported");
    }
}


/*
 *console output
 */
void CamPremix::report(double t, double* soln)
{
    static int nSteps=0;
    cout.width(5);
    cout.setf(ios::scientific);
    cout << t << endl;
    if(nSteps%10==0)reporter_->consoleHead("time(s)");
    nSteps++;
}


/*
 *console output with residual monitoring
 */
void CamPremix::report(double t, double* soln, double& res)
{
    static int nStep=0;
    cout.width(5);
    cout.setf(ios::scientific);
    if(nStep%10==0) reporter_->consoleHead("time(s) \t residual");
    cout << t <<"\t" << res << endl;
    nStep++;
}


void CamPremix::reportToFile(double t, double* soln)
{
    saveMixtureProp(t,soln,false,false);
    double sum =0;

    reporter_->openFiles();
    reporter_->writeHeader(headerData);
    vector<double> data, axpos;
    vector<double> molfrac, massfrac;
    axpos = reacGeom_.getAxpos();
    int len = axpos.size();

    for (int i=0; i<len; i++)
    {
        data.clear();
        data.push_back(t);
        data.push_back(axpos[i]);
        data.push_back(axpos[i]/m_u[i]);
        data.push_back(m_rho[i]);
        data.push_back(m_u[i]);
        data.push_back(soln[i*nVar+ptrF]);
        data.push_back(soln[i*nVar+ptrT]);

        massfrac.clear();
        molfrac.clear();
        for (int l=0; l<nSpc; l++)
        {
            massfrac.push_back(soln[i*nVar+l]);
        }

        if (admin_.getSpeciesOut() == admin_.MASS)
        {
            sum = 0;
            for (int l=0; l<nSpc; l++)
            {
                data.push_back(fabs(massfrac[l]));
                sum += massfrac[l];
            }
        }
        else
        {
            CamConverter control_;
            control_.mass2mole(massfrac,molfrac,*camMech_);
            sum = 0;
            for (int l=0; l<nSpc; l++)
            {
                data.push_back(fabs(molfrac[l]));
                sum += molfrac[l];
            }
        }

        data.push_back(sum);

        reporter_->writeStdFileOut(data);
    }

    reporter_->closeFiles();
}


//return the initial solution vector
void CamPremix::getInitial(vector<double>& initial)
{
    initial = solvect;
}


/*
 *update the diffusion fluxes. The east side of
 *last cell flux need to be set as zero
 */
void CamPremix::updateDiffusionFluxes()
{
    /*
     *call the base class function and then
     *pad the species outlet boundary condition
     *at the exit
     */
    CamResidual::updateDiffusionFluxes();
    for (int l=0; l<nSpc; l++)
    {
        s_jk(iMesh_e + 1,l) = 0.0;
    }
}


/*
 *update the thermal fluxes by calling
 *the base class function
 */
void CamPremix::updateThermo()
{
    CamResidual::updateThermo();
}


/*
 *save flow variables
 */
void CamPremix::saveFlowVariables(double* y)
{
    m_u.clear();
    m_flow.clear();

    for (int i=cellBegin; i<cellEnd; i++)
    {
        m_u.push_back(y[i*nVar+ptrF]/m_rho[i]);
        m_flow.push_back(y[i*nVar+ptrF]);
    }
}

void CamPremix::header()
{
    headerData.clear();
    headerData.push_back("int_time");
    headerData.push_back("x");
    headerData.push_back("tau");
    headerData.push_back("rho");
    headerData.push_back("u");
    headerData.push_back("mdot");
    headerData.push_back("T");
    for (int l = 0; l < nSpc; l++)
    {
        headerData.push_back( (*spv_)[l]->Name() );
    }
    headerData.push_back("sumfracs");
}


// Mass matrix evaluation
void CamPremix::massMatrix(double** M)
{
//    /*
//     *inlet boundary: all governing equations are algebraic except energy
//     */

//    /*
//     * on all interior mesh points except the mass flow
//     * all other governing equations are ODE. Two additional
//     * imaginary cells are padded to handle the inlet and exhaust
//     */
//    int nCells = reacGeom->getnCells()+2;

//    for(int i=0; i< nEqn; i++)
//        M[0][i] = 1.0;
//    /*
//     *at the inlet species conservation eqauations are algebraic
//     */
//    for(int l=0; l<nSpc; l++)
//        M[0][l] = 0.0;
//    /*
//     *mass flow
//     */
//    for(int i=0; i< nCells; i++){
//        int ptr = nVar*i + ptrF;
//        M[0][ptr] = 0.0;
//    }
//    /*
//     *at the exit both species and energy eqn are algebraic
//     */
//    int ptr;
//    for(int l=0; l<nSpc; l++){
//        ptr = nVar*iMesh_e+l;
//        M[0][ptr] = 0.0;
//    }
//    if(admin_->getEnergyModel() == admin_->ADIABATIC){

//        M[0][ptrT] = 0.0;

//        ptr = nVar*iMesh_e + ptrT;
//        M[0][ptr] = 0.0;
//    }
}

