/*!
* \file   stagflow.h
* \author V. Janardhanan
*
* \brief This class implements the stagnation flow and twin flame model.
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

#ifndef _STAGFLOW_H
#define _STAGFLOW_H

#include "gpc.h"
#include "cam_residual.h"
#include "cam_control.h"
#include "cam_control.h"
#include "cam_boundary.h"
#include "cam_profile.h"
#include "cam_setup.h"
#include "cam_profile.h"
#include "cam_configuration.h"

#include <vector>

namespace Camflow
{
    /*!
    *@brief    Stagnation flow / Counterflow equation solver class.
    *
    * Include a more detailed description here.
    */
    class StagFlow
    :
        public CamSetup
    {

        public:

            typedef struct
            {
                std::vector<double> a, b, c;
            } tdma_coeff;

            StagFlow
            (
                CamAdmin& ca,
                CamConfiguration& config,
                CamControl& cc,
                CamGeometry& cg,
                CamProfile& cp,
                CamSoot& cs,
                Mechanism& mech
            );

            virtual ~StagFlow(){}

            /*
            *solve function
            */
            void solve();

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
            int eval
            (
                double t,
                double* y,
                double *ydot,
                bool jacEval
            );

            void residual(const double& t, double* y, double* f);

            double getResidual() const;

            void speciesResidual
            (
                const double& time,
                double* y,
                double* f
            );

            void energyResidual
            (
                const double& time,
                double* y,
                double* f
            );

            void massMatrix(double** M);

            double dydx(double nr1, double nr2, double dr) const;

            double dydx
            (
                const int i,
                double nr1,
                double nr2,
                double nr3,
                double dr
            ) const;

            /*
            *calculate flow field
            */
            void calcFlowField(const double& time, double* y);

            void calcVelocity(std::vector<double>& flow);

            void calcMomentum(std::vector<double>& mom);

            /*
            *species boundary condition
            */
            void speciesBoundary
            (
                const double& t,
                double* y,
                double* f
            );

            /*
            *energy boundary
            */
            void energyBoundary
            (
                const double& t,
                double* y,
                double* f
            );

            /*
            *initialize the solution vector
            */
            void initSolutionVector(CamControl &cc);

            /*
            *initialize the mass flow
            */
            void initMassFlow();

            /*
            *initialize momentum
            */
            void initMomentum();

            /*
            *update diffusion fluxes
            */
            void updateDiffusionFluxes();

            /*
            *report functions
            */
            void report(double t, double* solution);

            void report(double t, double* solutio, double& res);

            void header(std::vector<std::string>& headerData);

            void reportToFile(double t, double* soln);

            //! Calculate the mixture fraction variable using Bilger's formula.
            double getBilgerMixFrac(const int& cell);


        private:

            inletStruct fuel, oxid;
            tdma_coeff tdmaFlow;
            int configID;

            double strainRate;

    }; // End StagFlow class declaration.

} // End Camflow namespace.

#endif /* _STAGFLOW_H */
