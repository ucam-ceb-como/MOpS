/*!
 * \file   cam_setup.h
 * \author V. Janardhanan
 *
 * \brief Class to setup initial conditions.
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

#ifndef _CAM_SETUP_H
#define	_CAM_SETUP_H

#include "cam_boundary.h"
#include "gpc.h"
#include "cam_params.h"
#include "cam_residual.h"
#include "cam_profile.h"
#include "cam_control.h"
#include <vector>

using namespace Sprog;

namespace Camflow{
    /*!
     *@brief    Class to setup initial conditions.
     *
     * Include a more detailed description here.
     */
    class CamSetup
    : public CamResidual {

        public:

            //! Structure containing inlet conditions.
            typedef struct{
                std::vector<doublereal> Species;
                doublereal FlowRate;
                doublereal Vel;
                doublereal rVelGrad;
                doublereal Dens;
                doublereal T;
                std::vector<doublereal> Dk;
                std::vector<doublereal> jk;
            } inletStruct;

            CamSetup(){}
            virtual ~CamSetup(){}

            //! Return the inlet species mass fractions for the given boundary.
            void getInletMassFrac(CamBoundary &cb, std::vector<doublereal>& fracs);

            //! Return the inlet temperature.
            const doublereal getInletTemperature(CamBoundary &cb);

            //! Return the inlet flow rate.
            const doublereal getInletFlowRate(CamBoundary &cb);

            //! Return the inlet velocity.
            const doublereal getInletVelocity(CamBoundary &cb);

            //! Initialise the species vector.
            void initSpecies(CamBoundary &cb, CamControl &cc, std::vector<doublereal>& soln);

            //! Initialize the species vector for a counter flow flame.
            void initSpecies(CamBoundary &left, CamBoundary &right,
                                    CamControl &cc, std::vector<doublereal>& soln);

            //! Initialise the mass flow rate.
            void initMassFlow(CamBoundary &cb, CamControl &cc, std::vector<doublereal> &soln);

            //! Initialise the temperature.
            void initTemperature(CamBoundary &cb, CamControl &cc, std::vector<doublereal> &soln);

            //! Initialise the temperature field with a Gaussian profile.
            void initTempGauss(std::vector<doublereal> &soln);

            //! Store the inlet conditions.
            void storeInlet(CamBoundary &cb, inletStruct& ud_inlet);

            //! Store the objects that are common across the reactor models.
            virtual void storeObjects(CamControl& cc,   CamAdmin& ca,    CamGeometry& cg,
                        CamProfile& cp,   CamBoundary& cb,    Mechanism& mech);


        protected:

            CamProfile *profile;

        }; // End CamSetup class declaration.

} // End Camflow namespace.

#endif	/* _CAM_SETUP_H */

