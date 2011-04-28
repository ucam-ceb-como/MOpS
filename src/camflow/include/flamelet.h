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

#include "cam_residual.h"
#include "cam_control.h"
#include "cam_admin.h"
#include "cam_reporter.h"
#include "cam_radiation.h"
#include "gpc.h"
#include "cam_setup.h"

using namespace Sprog;

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
            enum{
                LNONE,
                LNNONE
            };

            //! Default constructor. Give boolean flags some default values.
            FlameLet();

            //! Destructor.
            virtual ~FlameLet();

            //! Solve the flamelet. Call for coupling without solving population balance.
            void solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg, CamProfile&cp,
                                Mechanism &mech, bool interface=false);

            //! Stand alone call as well as first call from an external code that solves the population balance
            void solve(CamControl& cc, CamAdmin& ca, CamGeometry& cg, CamProfile& cp,
                    CamConfiguration& config, CamSoot& cs, Mechanism& mech);

            //! Continulation call from an external code that solves the population balance.
            void solve(std::vector<Thermo::Mixture>& cstrs,
                       const std::vector< std::vector<doublereal> >& iniSource,
                       const std::vector< std::vector<doublereal> >& fnlSource,
                       Mechanism& mech,
                       CamControl& cc,
                       CamAdmin& ca,
                       CamGeometry& cg,
                       CamProfile& cp);

            //! Initialize the solution vector.
            void initSolutionVector(CamControl &cc);

            //! Coupled solver.
            void csolve(CamControl &cc, bool interface=false);

            //! Segregated solver.
            void ssolve(CamControl &cc);

            //! Restart the solution with the converged solution.
            void restart(CamControl &cc);

            //! Save the solution vector.
            void saveMixtureProp(doublereal* y);

            //! Called by ODE/DAE solver to evaluate residuals.
            int eval(doublereal x, doublereal* y, doublereal* ydot, bool jacEval);

            //! Residual definitions.
            void residual(const doublereal& t, doublereal* y, doublereal* f);

            //! Species residual calculation.
            void speciesResidual(const doublereal& t, doublereal* y, doublereal* f);

            //! Energy residual calculation.
            void energyResidual(const doublereal& t, doublereal* y, doublereal* f);

            //! Set the external scalar dissipation rate.
            void setExternalScalarDissipationRate(const doublereal sr);

            //! Set the time history of scalar dissipation rate.
            void setExternalScalarDissipationRate(const std::vector<doublereal>& time,
                                                  const std::vector<doublereal>& sdr,
                                                  const bool analytic);

            //! Set the time history of scalar dissipation rate with a profile from the CFD.
            void setExternalScalarDissipationRate(const std::vector<doublereal>& time,
                                                  const std::vector< std::vector<doublereal> >& sdr,
                                                  const std::vector< std::vector<doublereal> >& Zcoords);

            //! Calculate the stoichiometric mixture fraction.
            doublereal stoichiometricMixtureFraction();

            //! Calculate the scalar dissipation rate.
            doublereal scalarDissipationRate(const doublereal m_frac);

            //! Get the SDR at a given time.
            doublereal getSDR(const doublereal time) const;

            //! Get the SDR at a given time and mixture fraction coordinate.
            doublereal getSDRfromProfile(const doublereal time, const doublereal Z) const;

            //! Get the SDR at a given mixture fraction coordinate using an analytical profile.
            doublereal scalarDissipationRateProfile(const doublereal m_frac, const doublereal stoichSDR, const int cell);

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
            void header();

            //! Console output.
            void report(doublereal x, doublereal* solution);

            //! Console output with residuals.
            void report(doublereal x, doublereal* solution, doublereal& res);

            //! File output.
            void reportToFile(doublereal x, doublereal* solution);

        private:

            doublereal stoichZ; //stoichiometric mixture fraction
            doublereal smr;     //stoichiometric mass ratio
            doublereal sdr;     // scalar dissipation rate
            doublereal sdr_ext; // scalar dissipation rate passed by exteranl program
            doublereal rstartTime;
            std::vector<doublereal> v_sdr;   //scalar dissipation rate that has a time history
            std::vector<doublereal> v_time; //time profile of scalar dissipation rates
            std::vector< std::vector<doublereal> > profile_sdr; // SDR profile with a time history.
            std::vector< std::vector<doublereal> > cfdMixFracCoords; // Mixture Fraction coords from CFD.

            //! Spatial profile of soot volume fraction
            std::vector<doublereal> m_SootFv;


            bool timeHistory, sdrProfile, sdrAnalytic;
            doublereal strain;  // strain rate
            int mCord;          // this is the mixture fraction coordinates
            inletStruct fuel, oxid;

            Array2D Le, convection, CpSpec; //Lewis numbers

            Radiation *radiation;

    }; // End FlameLet class declaration.

} // End Camflow namespace.

#endif	/* _FLAMELET_H */
