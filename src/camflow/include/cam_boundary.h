/*
 * File:   cam_boundary.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2009 Vinod M Janardhanan
 *
 * File purpose:
 *  This class contains the implementation boundary conditions
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
 *  Website :   http://como.cheng.cam.ac.uk
 *
 * Created on January 17, 2009, 4:50 PM
 */

#ifndef _CAM_BOUNDARY_H
#define	_CAM_BOUNDARY_H

#include "cam_params.h"
#include "cam_conc.h"
#include "gpc.h"
#include <string>
#include<map>
#include<vector>

using namespace Sprog;

namespace Camflow{
    class CamBoundary : public CamConc{
        doublereal velocity;
        doublereal flowRate;
        doublereal T;
        std::map<std::string, doublereal> species;
        std::vector<doublereal> inletFracs;
        int fracType;


        // Disallow copy constructor
        CamBoundary(const CamBoundary &cb);

    public:

        CamBoundary(){}
        ~CamBoundary(){}


        //set the inlet velocity
        void setVelocity(doublereal vel);

        //set the temperature
        void setTemperature(doublereal T);

        //set fraction type
        //void setFractionType(int n);

        //set flow rate
        void setFlowRate(doublereal flow);

        //set the inlet species
        void setSpecies(std::map<std::string, doublereal> species);

        //set the inlet mass fraction
        void setInletMassfracs(std::vector<doublereal> fracs);

        /*set the inlet mass fraction vector based on the
         *species map and the mechanism object. This function
         *may be called only after setting the inlet species.
         *This will do a valiity check for the species listed
         *in the inlet against the species available in the
         *mechanism. If any species listed in the inlet do not
         *belong to the mechanism the function will throw
         *an instance of CamError. If the validity check is
         *successful the function will return the inlet mass/mole vector
         */
        std::vector<doublereal>& setInletfracs( Mechanism &mech);

        //return the velocity
        const doublereal& getVelocity() const;

        //return the temperature
        const doublereal& getTemperature() const;

        //return the fraction type
        //const int& getFractionType() const;

        //return the flow rate
        const doublereal& getFlowRate() const;

        //return the inlet species
        const std::map<std::string,doublereal>& getInletSpecies() const;


        //return the inlet mass fractions
        const std::vector<doublereal>& getInletMassfracs() const;

    };
}
#endif	/* _CAM_BOUNDARY_H */

