/*!
 * \file   flamelet.h
 * \author V. Janardhanan
 *
 * \brief Solver class for flamelets.
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

#ifndef _FLAMELET_H
#define	_FLAMELET_H

#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "cam_control.h"
#include "cam_radiation.h"
#include "gpc.h"
#include "array.h"
#include "cam_setup.h"
#include "cam_residual.h"
#include "cam_reporter.h"
#include "cam_profile.h"
#include "cam_admin.h"
#include "flamelet.h"
#include "cam_math.h"
#include "cvode_wrapper.h"
#include "radau_wrapper.h"
#include "limex_wrapper.h"
#include "linear_interpolator.hpp"
#include "cam_sdr.h"

namespace Camflow
{

    /*!
     *@brief    Flamelet equation solver class.
     *
     * Include a more detailed description here.
     */
    class FlameLet
    :
      public CamSetup
    {

        public:

            //! Enumerator for Unity or Fixed Lewis number.
            enum
            {
                LNONE,
                LNNONE
            };

            //! Default constructor. Give boolean flags some default values.
            FlameLet
            (
                CamAdmin& ca,
                CamConfiguration& config,
                CamControl& cc,
                CamGeometry& cg,
                CamProfile& cp,
                CamSoot& cs,
                Sprog::Mechanism& mech
            );

            //! Destructor.
            virtual ~FlameLet();

            // Run some tests to check the setup.
            void checkSetup();

            //! Solve the flamelet. Call for coupling without solving population balance.
            void solve(bool interface);

            //! Stand alone call as well as first call from an external code that solves the population balance
            void solve();

            // Call from openFoam with flag to indicate if we are solving a steady state
            // flamelet (at the base of the flame) with no soot, or if we are solving a
            // Lagrangian flamelet with soot (if soot mech is present)
            void solve(bool interface, bool steadyStateNoSoot);

            //! Continulation call from an external code that solves the population balance.
            void solve(std::vector<Thermo::Mixture>& cstrs,
                       const std::vector< std::vector<doublereal> >& iniSource,
                       const std::vector< std::vector<doublereal> >& fnlSource,
                       Sprog::Mechanism& mech,
                       CamControl& cc,
                       CamAdmin& ca,
                       CamGeometry& cg,
                       CamProfile& cp);

            //! Initialize the solution vector.
            void initSolutionVector();

            //! Coupled solver.
            void csolve(bool interface=false);

            //! Segregated solver.
            void ssolve(bool interface=false);

            //! Split solver.
            void splitSolve(bool interface=false);

            //! Restart the solution with the converged solution.
            void restart();

            //! Save the solution vector.
            void saveMixtureProp(doublereal* y);

            //! Called by ODE/DAE solver to evaluate residuals.
            int eval(doublereal x, doublereal* y, doublereal* ydot, bool jacEval);

            //! Residual definitions.
            void residual(const doublereal& t, doublereal* y, doublereal* f);

            //! Species residual calculation.
            void speciesResidual(const doublereal& t, doublereal* y, doublereal* f);

            //! Soot moments residual calculation.
            void sootMomentResidual(const doublereal& t,doublereal* y,doublereal* f);

            //! Soot moments residual calculation. (When SS at base of flame)
            void sootMomentResidualAtFlameBase(const doublereal& t,doublereal* y,doublereal* f);

            //! Energy residual calculation.
            void energyResidual(const doublereal& t, doublereal* y, doublereal* f);

            //! Calculate the stoichiometric mixture fraction.
            doublereal stoichiometricMixtureFraction();

            void setExternalStrainRate(const doublereal strainRate);

            void setExternalSDR(const doublereal sdr);

            void setExternalTimeSDR
            (
                const std::vector<doublereal>& time,
                const std::vector<doublereal>& sdr
            );

            //! Provide a soot volume fraction from an external calculation.
            void setExternalSootVolumeFraction(const std::vector<doublereal>& soot_fv);

            //! Get the residual for use in radauWrapper.
            doublereal getResidual() const;

            //! Mass matrix evaluation. Called by radauWrapper.
            void massMatrix(doublereal **M);

            //! Set the restart time for the solver.
            void setRestartTime(doublereal t);

            //! Get a vector for the rate of production of pyrene.
            void getWdotA4(std::vector<doublereal>& wdotA4) const;

            //! Create header for file output.
            std::vector<std::string> header();

            //! Create header for soot rates file output.
            std::vector<std::string> sootRatesHeader();

            //! Console output with residuals.
            void report(doublereal x, doublereal* solution, doublereal& res);

            //! File output.
            void reportToFile(std::string fileName, doublereal x, std::vector<double>& solution);

            // Output soot compoent rates to file
            void reportSootRatesToFile(std::string fileName, doublereal t, Array2D& rates);



        private:

            doublereal stoichZ; //stoichiometric mixture fraction
            doublereal sdr;     // scalar dissipation rate
            doublereal restartTime;
            bool steadyStateAtFlameBase;
            std::vector< std::vector<doublereal> > cfdMixFracCoords; // Mixture Fraction coords from CFD.

            //! Spatial profile of soot volume fraction
            std::vector<doublereal> m_SootFv;


            bool timeHistory, sdrProfile, sdrAnalytic;
            inletStruct fuel, oxid;

            Radiation *radiation;
            ScalarDissipationRate scalarDissipationRate_;

            Array2D Le, convection, CpSpec; //Lewis numbers

    }; // End FlameLet class declaration.

} // End Camflow namespace.

#endif	/* _FLAMELET_H */
