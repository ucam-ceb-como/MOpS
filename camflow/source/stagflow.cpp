
#include <stdlib.h>
#include <vector>

#include "array.h"

#include "cam_geometry.h"


#include "cam_setup.h"

/*
 * File:   stagflow.cpp
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class implements the stagnation flow and twin flame model
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
 *  Website :   http://como.cheng.cam.ac.uk
 *
 * Created on 30 April 2009, 15:03
 */
#include "cam_reporter.h"
#include "cam_residual.h"
#include "cvode_wrapper.h"
#include "kinsol_wrapper.h"
#include "stagflow.h"

using namespace Camflow;
/*
 *solve the stagnation/twinflame model
 */
void StagFlow::solve(CamControl& cc, CamAdmin& ca, CamGeometry& cg,
                            CamProfile& cp, Mechanism& mech){
    /*
     *function to set up the solver.
     *Initialisation of the solution vector is
     *done here
     */

    camMech = &mech;
    Thermo::Mixture mix(mech.Species());
    camMixture = &mix;
    spv = camMixture->Species();
    profile = &cp;

    
    admin = &ca;
    reacGeom = &cg;
    reacGeom->descretize();
    /*
     * 2 additional cells are padded to consider the
     * inlet and the exhaust
     */

    reacGeom->addZeroWidthCells();

    opPre = ca.getPressure();
    profile->setGeometryObj(cg);
    reporter = new CamReporter();
    nSpc = camMech->SpeciesCount();

    /*
     *array offsets for newton solver
     */
    ptrF = 0;
    ptrG = ptrF+1;
    ptrH = ptrG+1;
    /*
     *array offsets for ODE
     */
    ptrT = nSpc;

    

    /*
     *number of finite volume cells and
     *loop initializations
     */
    nCells = reacGeom->getnCells();
    loopBegin = 1;
    loopEnd = nCells-1;

    cellBegin = 0;
    cellEnd = nCells;


    /*
     *number of variables and the number of equations
     *to be solved by the newton solver
     */
    alg_nVar = 3; // F(x), G(x), H
    alg_nEqn =  alg_nVar*nCells;
    alg_band = alg_nVar;
    

    /*
     *number of equations and the number of variables to
     *be slved by the ODE
     */
    nVar = nSpc+1;
    nEqn = nVar*nCells;

    /*
     *init the solution vector
     */
    initSolutionVector(cc);

}
/*
 *coupled solver
 */
void StagFlow::csolve(CamControl& cc){

}

/*
 *segregated solver
 */
void StagFlow::ssolve(CamControl& cc){
    
}
/*
 *function called by the solver (DAEs and ODEs)
 */
int StagFlow::eval(doublereal t, doublereal* y, doublereal* ydot, bool jacEval){

    return 0;
}
/*
 *function called by newton solver
 */
int StagFlow::eval(doublereal* y, doublereal* ydot){
    /*
     *this calls the residual functions for
     *mass flow, momentum, and radial pressure gradient
     *which are solved by the newton solver
     */
    return 0;
}
/*
 *mass flow residual
 */
void StagFlow::massFlowResidual(const doublereal& time, doublereal* y,
                                                doublereal* f){

    

}
/*
 *initialize the solution vector
 */
void StagFlow::initSolutionVector(CamControl& cc){
    
    /*
     *initialize the geometry
     */
    dz = reacGeom->getGeometry();

    /*
     *set the pressure
     */
    opPre = admin->getPressure();



    CamBoundary left, right;
    admin->getLeftBoundary(left);
    admin->getRightBoundary(right);
    
    /*
     *store the inlets and change the
     *sign of flow for the oxidizer
     */
    storeInlet(left,fuel);
    storeInlet(right,oxid);
    
    oxid.FlowRate *= -1;
    oxid.Vel *= -1;
    right.setVelocity(-right.getVelocity());
    right.setFlowRate(-right.getFlowRate());
    

    /*
     *initialize the ODE vector
     */
    solvect.resize(nEqn,0);
    vector<doublereal> vSpec, vT;
    initSpecies(left,right,cc,vSpec);
    initTemperature(left,cc,vT);

    mergeEnergyVector(&vT[0]);
    mergeSpeciesVector(&vSpec[0]);

    /*
     *initialize the algebraic equations vector
     */
    alg_solvect.resize(alg_nEqn,0.0);
    vector<doublereal> vFlow, vMom, vPGrad;
    initMassFlow(fuel.FlowRate, oxid.FlowRate,vFlow);
    initMomentum(fuel.FlowRate,oxid.FlowRate,vMom);
    initPressureGrad(vPGrad);
    /*
     *merge the dependent vaiables into
     *algebraic solution vector
     */
    for(int i=0; i<nCells; i++){
        alg_solvect[i*alg_nVar+ptrF] = vFlow[i];
        alg_solvect[i*alg_nVar+ptrG] = vMom[i];
        alg_solvect[i*alg_nVar+ptrH] = vPGrad[i];
    }
    

}
/*
 *init pressure gradianet
 */
void StagFlow::initPressureGrad(vector<doublereal>& soln){
    soln.resize(nCells,-1);
    for(int i=0; i<nCells/2; i++)
        soln[i] = -100;

}
/*
 *initialize the mass flow
 */
void StagFlow::initMassFlow(const doublereal fLeft, const doublereal fRight,
                            vector<doublereal>& soln){
    soln.resize(nCells,0);
    vector<doublereal> position = reacGeom->getAxpos();    
    doublereal flowGrad = (fRight-fLeft)/reacGeom->getLenth();
    
    for(int i=0; i<nCells; i++){

        soln[i] = flowGrad*(position[i]-position[0]) + fLeft;
    }

}
/*
 *initialize momentum
 */
void StagFlow::initMomentum(const doublereal fLeft, const doublereal fRight,
                            vector<doublereal>& soln){

    vector<doublereal> position = reacGeom->getAxpos();
    doublereal flowGrad = (fRight-fLeft)/reacGeom->getLenth();
    soln.resize(nCells,flowGrad);

}

/*
 *save the flow variables
 */
void StagFlow::saveFlowVariables(doublereal* y){
    
    m_flow.resize(nCells,0);
    m_G.resize(nCells,0);
    for(int i=cellBegin; i<cellEnd; i++){
        
        m_flow[i] = y[i*nVar+ptrF];
        m_G[i] = y[i*nVar+ptrG];

    }
}
/*
 *report functions
 */
void StagFlow::report(doublereal t, doublereal* solution){

}

void StagFlow::report(doublereal t, doublereal* solutio, doublereal& res){

}
