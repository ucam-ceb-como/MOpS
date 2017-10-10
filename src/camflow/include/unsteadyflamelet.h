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

#ifndef _UNSTEADYFLAMELET_H
#define	_UNSTEADYFLAMELET_H

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
#include "unsteadyflamelet.h"
#include "cam_math.h"
#include "cvode_wrapper.h"
#include "radau_wrapper.h"
#include "limex_wrapper.h"
#include "linear_interpolator.hpp"
#include "cam_sdr.h"
#include "cam_lewis.h"

namespace Camflow
{

    /*!
     *@brief    UnsteadyFlameLet equation solver class.
     *
     * eb656:
     * 
     * Test class to solve MOMIC equations for soot as unsteady
     * Idea: Replicate camflow side of momicFOAM coupling:
     *      Solve steady state flamelet calculation at base of the flame
     *          ie. No soot and no radiation
     *      With SS solution as IC, solve unsteady flamelet from 0 to $\tau$.
     *
     * This class is a near duplicate of the flamelet class, difference should be how it solves/restarts (see idea above)
     * 
     */
    class UnsteadyFlameLet
    :
      public CamSetup
    {

        public:


            //! Default constructor. Give boolean flags some default values.
            UnsteadyFlameLet
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
            virtual ~UnsteadyFlameLet();

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
                       const std::vector< std::vector<double> >& iniSource,
                       const std::vector< std::vector<double> >& fnlSource,
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
            void restart(double flameTime);

            //! Save the solution vector.
            void saveMixtureProp(double* y);

            //! Called by ODE/DAE solver to evaluate residuals.
            int eval(double x, double* y, double* ydot, bool jacEval);

            //! Residual definitions.
            void residual(const double& t, double* y, double* f);

            //! Species residual calculation.
            void speciesResidual(const double& t, double* y, double* f);

            //! Soot moments residual calculation.
            void sootMomentResidual(const double& t,double* y,double* f);

            //! Soot moments residual calculation. (When SS at base of flame)
            void sootMomentResidualZeroedOut(const double& t,double* y,double* f);

            //! Energy residual calculation.
            void energyResidual(const double& t, double* y, double* f);

            //! Calculate the stoichiometric mixture fraction.
            double stoichiometricMixtureFraction();

            void setExternalStrainRate(const double strainRate);

            void setExternalSDR(const double sdr);

            void setExternalTimeSDR
            (
                const std::vector<double>& time,
                const std::vector<double>& sdr
            );

            //! Provide a soot volume fraction from an external calculation.
            void setExternalSootVolumeFraction(const std::vector<double>& soot_fv);

            //! Get the residual for use in radauWrapper.
            double getResidual() const;

            //! Mass matrix evaluation. Called by radauWrapper.
            void massMatrix(double **M);

            //! Set the restart time for the solver.
            void setRestartTime(double t);

            //! Get a vector for the rate of production of pyrene.
            void getWdotA4(std::vector<double>& wdotA4) const;

            //! Create header for file output.
            std::vector<std::string> header();

            //! Create header for soot rates file output.
            std::vector<std::string> sootRatesHeader();

            //! Console output with residuals.
            void report(double x, double* solution, double& res);

            //! File output.
            void reportToFile(std::string fileName, double x, std::vector<double>& solution);

            // Output soot compoent rates to file
            void reportSootRatesToFile(std::string fileName, double t, Array2D& rates);



        private:

            double stoichZ; //stoichiometric mixture fraction
            double sdr;     // scalar dissipation rate
            double restartTime;
            bool sootResidualZeroed;
            //double sootFlameTimeThreshold;  // Flametime beyond which we switch off soot.
            std::vector< std::vector<double> > cfdMixFracCoords; // Mixture Fraction coords from CFD.

            //! Spatial profile of soot volume fraction
            std::vector<double> m_SootFv;


            bool timeHistory, sdrProfile, sdrAnalytic;
            inletStruct fuel, oxid;

            Radiation *radiation;
            ScalarDissipationRate scalarDissipationRate_;
            LewisNumber Lewis;

            Array2D convection, CpSpec;

    }; // End UnsteadyFlameLet class declaration.

} // End Camflow namespace.

#endif	/* _UNSTEADYFLAMELET_H */
