/*
 * File:   stagflow.h
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

#ifndef _STAGFLOW_H
#define	_STAGFLOW_H
#include "gpc.h"
#include "cam_residual.h"
#include "cam_control.h"
#include "cam_control.h"
#include "cam_boundary.h"
#include "cam_profile.h"
#include "cam_setup.h"
#include "cam_profile.h"

namespace Camflow{
    class StagFlow : public CamSetup {

    public:
        StagFlow(){}
        virtual ~StagFlow(){}

        /*
         *solve function
         */
        void solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg,CamProfile &cp, Mechanism &mech );
        /*
         *coupled solver
         */
        void csolve(CamControl &cc);
        /*
         *segregated solver
         */
        void ssolve(CamControl &cc);
        /*
         *function called by the solver (DAE and ODEs)
         */
        int eval(doublereal t, doublereal* y, doublereal *ydot, bool jacEval);
        /*
         *function called by newton solver
         */
        int eval(doublereal* y, doublereal* ydot);
        /*
         *mass flow residual
         */
        void massFlowResidual(const doublereal& time, doublereal* y, doublereal *f);
        /*
         *initialize the solution vector
         */
        void initSolutionVector(CamControl &cc);
        /*
         *initialize pressure grad
         */
        void initPressureGrad(vector<doublereal> &soln);
        /*
         *initialize the mass flow
         */
        void initMassFlow(const doublereal fLeft, const doublereal fRight,
                    vector<doublereal>& soln);
        /*
         *initialize momentum
         */
        void initMomentum(const doublereal fLeft, const doublereal fRight,
                    vector<doublereal>& soln);
        /*
         *save mixture properties
         */
        void saveFlowVariables(doublereal* y);
        /*
         *report functions
         */
        void report(doublereal t, doublereal* solution);
        void report(doublereal t, doublereal* solutio, doublereal& res);
    private:
        inletStruct fuel, oxid;
        int nCells;
        /*
         *Newton solver variables
         */
        int alg_nEqn, alg_nVar, alg_band;
        vector<doublereal> alg_solvect;
    };
}

#endif	/* _STAGFLOW_H */

