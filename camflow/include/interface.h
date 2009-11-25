/* 
 * File:   interface.h
 *Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2009 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class contains the definition of interface for coupling
 *  with other programs, which handles the  particle processes.
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
 * Created on 20 August 2009, 17:20
 */

#ifndef _INTERFACE_H
#define	_INTERFACE_H

#include <iostream>
#include <string>
#include "gpc.h"
#include "cam_read.h"
#include "cam_admin.h"
#include "cam_control.h"
#include "cam_geometry.h"
#include "cam_converter.h"
#include "cam_boundary.h"
#include "cam_profile.h"
#include "cam_configuration.h"
#include "cam_read.h"
#include "cam_models.h"
#include "cam_soot.h"
using namespace std;
using namespace Sprog;

/*--------------------------------------------------------------------------------------------------
 * How to use Code Snippet:
 *
 *  Creating the interface object
 * -----------------------------------------------------------------
 * The interface object can be created passing in the
 * mechanism object, grid sizes, and a vector of
 * mixtures to hold the mixture object for each
 * finite volume cells. The following code snippet
 * demonstrates inistantiating the interface a
 * premix reactor
 *
 * CamPremix cp;
 * Interface camFlowInterface(mech_in, dz, mixtures,&cp)
 *
 *
 * @mech_in        : Sprog::Mechanism object
 * @dz                 : Grid points vector of type double
 *                           This vector contains the width of
 *                            finite volume cells
 * @cstrs              : vector of type
 *                           Sprog::Thermo::Mixture
 * &cp                 : pointer to premix reactor
 *
 * dz, and mixtures can be empty vectors. If these
 * arguments are provided with the corresponding
 * contents then the size of  cstrs must be same
 * as dz.  if the size of dz is foubd to be zero then
 * the reactor geometry will be read in from the
 * grid.inp file, which is the grid specification file for
 * camflow. If the size of cstrs is found to be
 * zero, then it will be populated with n number of
 * mixture objects, where n is the number of cells
 *
 * Interface instantiation is followed by the
 * integration and the control will be, returned to
 * the calling program once steady state is obtained.
 * Once the control is returned to the calling
 * program, handles for all objects, which controls
 * execution of camflow can be obtained by
 * calling the get function. The most important
 * oncis the CamControl object
 *
 *  Getting the cam control object
 * -----------------------------------------------------------
 *  camFlowInterface.getCamControl(CamControl& ccObj);
 *
 *  Getting configuration object
 *---------------------------------------------------------------------------------
 *  camFlowInterface.getCamConfiguration(CamConfiguration& configObj);
 *
 *
 *  Making successive calls
 *--------------------------------------------------------------------------------------
 *  camFlowInterface.solve(vector<Thermo::Mixture>& cstrs,
                const vector<doublereal>& dz,
                const vector< vector<doublereal> >& initalSource,
                const vector< vector<doublereal> >& finalSource,
                CamControl& ccObj,
                CamConfig& confObj,
                Mechanism& mech_in,
                void* reactorModel
                );
 *
 *@cstrs        : vector of Sprog::Thermo:;Mixture
 *@dz           : vector of type double containing
 *                   cell widths
 *@initialSource: vector of source terms for
 *                   particle processes at time 0
 *@finalSource: vector of source term for particle
 *                  processes at the end of time t+dt
 *@ccObj        : object of the type CamControl. This
 *                    object must be set with the maximum
 *                    integration time by calling function
 *  setMaxTIme(double time). A number of
 * other optional parameters can also be set.
 * Refer to the respective header file.
 *
 *@mech_in  : Sprog::Mechanism object
 *@*reactorModel : This is a pointer to any of the
 * following reactor type (see instantiation of
 * Interface for how to pass the pointer)
 *
 *      CamPremix
 *      CamPlug
 *      Batch
 *      StagFlow
 *
 */

namespace Camflow{
    class Interface{
    public:
        /*
         *default interface; for external calling programs which
         *do not read the mechanism file.
         */
        Interface();
        /*
         *Interface for external calling program which reads in the
         *mechanism file. This interface will also perform the
         *integration until steady state, and return the species massfractions,
         *velocity, density, and temprature the mixture objects.
         *The mechanism object and the pointer to the reactor
         *model are mandatory input parameters, while the vector
         *containing the width of cells (dz) and the vector of
         *mixture objects, which corresponds to each finite volume cells
         *are optional.
         *
         */
        Interface(Mechanism& mech_in, 
                vector<doublereal>& dz,
                vector<Thermo::Mixture>& cstrs,
                void* rModel,
                const doublereal sdr=0
                );
        /*
         *Default destructor
         */
        ~Interface(){}
        /*
         *calling interface to solve the flameletes. continuation true means the
         *call is a continuation call to the same problem
         */
        void flamelet(doublereal sdr, doublereal intTime=0, bool continuation=false, bool lnone=true);
        /*
         *calling interface to solve flamelets with time-history-scalar dissipation rates
         */
        void flamelet(const vector<doublereal>& sdr, const vector<doublereal>& intTime, bool continuation=false, bool lnone=true);
        /*
         *assigen the species mass fractions, temperature, density, and
         *velocity into the vector of mixture objects
         */
        void resetMixtures(vector<Thermo::Mixture>& cstrs);
        /*
         *return the vector of species names
         */
        void getSpeciesNames(vector<string>& names);
        /*
         *return the number of species
         */
        const int getNumberOfSpecies() const ;
        /*
         *return the species mass frac given the indenpendant variable
         */
        const doublereal getMassFrac(const int spIndex, const doublereal axpos);
        /*
         *return the temperature given the independant variable
         */
        const doublereal getTemperature(const doublereal axpos);
        /*
         *return the density given the independant variable
         */
        const doublereal getDensity(const doublereal axpos);
        /*
         *return the viscosity
         */
        const doublereal getViscosity(const doublereal axpos);

        /*
         *  return the specific heat at a given mixture fraction /axpos
         */
        const doublereal getSpecificHeat(const doublereal axPos);
        /*
         *  return the thermal conductivity at a given mixture fraction / axpos
         */
        const doublereal getThermalConductivity(const doublereal axPos);
        /*
         *  Return a vector of diffusion coefficents (size number of gasphase species)
         */
        const vector<doublereal> getDiffusionCoefficients(const doublereal axPos);
        /*
         *return the stoichiometric mixture fraction
         */
        const doublereal getStMixtureFrac();

        /*
         *return the controller
         */
        CamControl& getCamControl();
        /*
         *return the geometry object
         */
        CamGeometry& getCamGeometry();
        /*
         *return CamAdmin
         */
        CamAdmin& getCamAdmin();
        /*
         *return the boundary object which stores the boundary condition
         * information
         */
        CamBoundary& getCamBoundary();
        /*
         *return the profile object which contains the information
         *of user defined temperature profile
         */
        CamProfile& getCamProfile();
        /*
         *return the configuration that stores the reactor type
         */
        CamConfiguration& getCamConfiguration();

        /*
         *solve the reactor for the successive calls
         */
        void solve(vector<Thermo::Mixture>& cstrs,                      //Sprog mixture for each cell
                const vector<doublereal>& dz,                                   //geomtry cell widths
                const vector< vector<doublereal> >& initalSource,       //initial source terms ( under relaxed )
                const vector< vector<doublereal> >& finalSource,         //final source terms (under relaxed)
                CamControl& ccObj,                                                    //Solver control
                CamConfiguration& confObj,                                   //configuration object
                Mechanism& mech_in,                                             //mechanism object
                void* reactorModel,                                                  //reactor model to solve
                const doublereal sdr = 0                                        //scalar dissipation rate in case of flamelets
                );
        

    private:
        CamControl cc;
        CamGeometry cg;
        CamConverter convert;
        CamAdmin ca;
        CamBoundary cb;
        CamProfile cp;
        CamConfiguration config;
        CamRead cm;
        CamResidual *model;

        CamSoot cSoot;        
        Mechanism mech;
        int nSpecies;
        vector<string> speciesNames;
        /*
         *memberes to hold dependent variabes
         */
        Array2D spMassFracs;            //species mass fractions
        vector<doublereal> TVector;     //temperature
        vector<doublereal> rhoVector;   //density
        vector<doublereal> muVector;    //viscosity
        vector<doublereal> indVar;      //independant variable
        vector<doublereal> spHeat;              //specific heats
        vector<doublereal> lambda;      //thermal conductivity
        Array2D mDiff;                      //Diffusion coefficients

        doublereal stMixtureFrac;       //stoichiometric mixture fraction
        /*
         *reactor pinters
         */
        FlameLet* flmlt;

        doublereal getVariableAt(const doublereal& pos, vector<doublereal>& var);

    };
}


#endif	/* _INTERFACE_H */

