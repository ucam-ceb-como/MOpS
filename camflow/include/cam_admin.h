/* 
 * File:   cam_admin.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2009 Vinod M Janardhanan
 *
 * File purpose:
 *  This class contains the implementation of process conditions
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
 * Created on January 17, 2009, 4:36 PM
 */

#ifndef _CAM_ADMIN_H
#define	_CAM_ADMIN_H

#include "cam_params.h"
#include "cam_boundary.h"
#include "comostrings.h"

namespace Camflow{
    class CamAdmin : public CamConc {
        doublereal T, T_wall;     //temperature in case of isothermal operation
        doublereal pre;   //constant pressure system
        doublereal strain; // strain rate
        doublereal stepIgnition;    //step change of temperature for ignition
        int energyModel;
        int speciesOut;  //identifier for species output mass/mole
        int repSch;     //report schedule
        CamBoundary left; //fuel 
        CamBoundary right;//oxidizer
    public:
        enum EnergyModel{
            ADIABATIC,
            ISOTHERMAL,
            USERDEFINED,
            NONISOTHERMAL
        };

        enum Report{
            FINAL,
            INTER
        };
        //set the temperature
        //void setTemperature(doublereal T_);

        //set the operating pressure
        void setPressure(doublereal p_);

        //set strain rate
        void setStrainRate(doublereal sr);

        //set the step for ignition temp
        void setIgnitionStep(doublereal step);

        //set the wall temperature
        void setWallTemp(doublereal Tw);

        //set the energy model
        void setEnergyModel(int n);

        //set the energymodel
        void setEnergyModel(string model);

        //set left side boundary (normally for fuel inlet)
        void setLeftBoundary(CamBoundary &cb);

        //set right boundary (normally for oxidizer inlet)
        void setRightBoundary(CamBoundary &cb);

        //set the species output
        void setSpeciesOut(int n);

        int getSpeciesOut() const;

        //return the temperature
        //doublereal getTemperature() const;

        //return the wall temperature
        doublereal getWallTemp() const;
        
        //return the pressure
        doublereal getPressure()const;

        //return the strain rate
        doublereal getStrainRate() const;

        //return the ignition step
        doublereal getIgnitionStep() const;

        //return the energy model
        int getEnergyModel() const;

        //return the monitor option
        int getMonitor() const;

        //return fuel boundary
        void getRightBoundary(CamBoundary &cb);

        //return oxidizer boundary
        void getLeftBoundary(CamBoundary &cb) ;

        //set the report schedule
        void setReportSchedule(string repo);
        void setReportSchedule(int n);

        //return the report schedule
        int getReportSchedule();

        //Evaluation of dimensionaless numbers
        doublereal getNre(const doublereal& hd, //hydraulic dia
                const doublereal& u, //velocity
                const doublereal& rho,//density
                const doublereal& eta // viscosity
                );

        //return the prandtl number
        doublereal getPrandtl(const doublereal& eta, //viscosity
                const doublereal& lambda,//thermal conductivity
                const doublereal& cp //specific heat
                );

        //return the Graetz number
        doublereal getGraetzInv(const doublereal& x,//axial position
                    const doublereal& dh,   //hydraulic dia
                    const doublereal& Nre,  //Reynolds number
                    const doublereal& Pr    //Prandtl number
                    );

        //return the nussult number
        doublereal getNusselt(const doublereal& gzInv);

        //return the heat transfer coefficient
        doublereal getHeatTransferCoeff(const doublereal& x,
                                        const doublereal& vel,
                                        const doublereal& hd,
                                        const doublereal& rho,
                                        const doublereal& eta,
                                        const doublereal& lambda,
                                        const doublereal& cp);

    };
}


#endif	/* _CAM_ADMIN_H */

