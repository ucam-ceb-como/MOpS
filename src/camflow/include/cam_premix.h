/*
 * File:   cam_premix.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class implements the premix flame model
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
 * Created on January 18, 2009, 7:57 PM
 */

#ifndef _CAM_PREMIX_H
#define	_CAM_PREMIX_H

#include "gpc.h"
#include "cam_residual.h"
#include "cam_control.h"
#include "cam_control.h"
#include "cam_boundary.h"
#include "cam_profile.h"
#include "cam_setup.h"
#include "cam_profile.h"
#include "cam_configuration.h"
#include "cam_soot.h"

using namespace Sprog;

namespace Camflow
{

    class CamPremix
    :
        public CamSetup
    {

        inletStruct ud_inlet;

        inletStruct soln_inlet;

        double resNorm;

        std::vector<std::string> headerData;

        double dydx(double nr1, double nr2, double dr) const;

        double dydx(double nr1, double nr2, double nr3, double dr) const;


    public:

        CamPremix
        (
            CamAdmin& ca,
            CamConfiguration& config,
            CamControl& cc,
            CamGeometry& cg,
            CamProfile& cp,
            CamSoot& cs,
            Mechanism& mech
        );

        virtual ~CamPremix(){};

        /*
         *the following 3 functions are called by
         *the DAE solver
         */
        int eval(double x, double* y, double* ydot, bool jacEval);

        /*
         * report the results without residual monitoring
         */
        void report(double t, double* solution);

        /*
         * console output with resisual monitoring
         */
        void report(double t, double* solution, double& res);

        /*
         * write the results to output file
         */
        void reportToFile(double t, double* soln);

        /*
         *solve the premix reactor for the stand alone case
         */
        void solve();

        /*
         *solve the premix reactor problem for continuation calls from the
         *external interface
         */
        void solve
        (
            std::vector<Thermo::Mixture>& cstrs,
            const std::vector< std::vector<double> >& iniSource,
            const std::vector< std::vector<double> >& fnlSource,
            Mechanism& mech,
            CamControl& cc,
            CamAdmin& ca,
            CamGeometry& cg,
            CamProfile& cp
        );

        /*
         *  Save the inlet after solution for the
         *  next call from population balance solver
         */
        void saveInlet();

        /*
         *coupled solver. In this case all the equations are solved
         *simultaneousl
         */
        void csolve();

        /*
         *segregated solver. It is found that CVode failes to
         *handle the energy equation for large problems. In this case
         *segregated solvers are used to solve the equation system by
         *decoupling the energy equation and the species transport equation.
         *Once the specified number of iterations are completed, the control
         *is given back to the coupled solver to complete the intergration
         */
        void ssolve();

        /*
         * return the initial solution vector
         */
        void getInitial(std::vector<double>& initial);

        /*
         *residual function definition.
         */
        void residual(const double& t, double *y, double *f);

        /*
         *definition of boundary condition for the total continuity
         */
        void massFlowBoundary
        (
            const double& t,
            double *y,
            double *f
        );

        /*
         * base class definition for mass flow. This function can be used for
         * any boundary value problems. The boundary condition has to be implemented
         * in the respective reactor models
         */
        void massFlowResidual
        (
            const double& time,
            double* y,
            double* f
        );

        /*
         *definition of species boundary conditions
         */
        void speciesBoundary(const double& t, double *y, double *f);

        /*
         *definition of energy equation boundary conditions
         */
        void energyBoundary(const double& t, double *y, double *f);

        /*
         *definition of moment conditions for soot model
         */
        void momentBoundary(const double& t, double *y, double *f);

        /*
         * set up the solution vector of dependent variables
         */
        void initSolutionVector(CamBoundary& cb);

        /*
         * header information to write the file output
         */
        void header();

        /*
         * update the species diffusion fluxes to solve the species
         * transport equations
         */
        void updateDiffusionFluxes();

        /*
         * update the thermal fluxes for solving the energy transport
         */
        void updateThermo();

        /*
         * save the flow variables. i.e. extract the bulk velocity
         * and store it for later use
         */
        void saveFlowVariables(double* y);

        double getResidual() const { return 0.0; };

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
    };

} // End Camflow Namespace

#endif  /* _CAM_PREMIX_H */

