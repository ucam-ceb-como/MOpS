/*!
 * \file   interface.h
 * \author V. Janardhanan
 *
 * \brief Interface for coupling flamelet calculations to external codes.
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
#include "cam_params.h"

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
                const vector<double>& dz,
                const vector< vector<double> >& initalSource,
                const vector< vector<double> >& finalSource,
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

namespace Camflow {
/*!
 *@brief    Wrapper to combine flamelet calculations with external codes (e.g. OpenFOAM).
 *
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
                const vector<double>& dz,
                const vector< vector<double> >& initalSource,
                const vector< vector<double> >& finalSource,
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
 */
    class Interface{

        public:

            //! Default interface.
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
                    std::vector<double>& dz,
                    std::vector<Thermo::Mixture>& cstrs,
                    void* rModel,
                    const double sdr=0
                    );

             //! Default destructor.
            ~Interface();

            //! Calling interface to solve flamelets with time-history scalar dissipation rates.
            void flamelet(const std::vector<double>& sdr, const std::vector<double>& intTime, bool continuation=false);

            //! Call to solve a flamelet for a given strain rate.
            void flameletStrainRate(const double& strainRate);

            //! Call to solve a flamelet for a given SDR rate.
            void flameletSDR(const double& SDR);

            //! Calling interface to solve flamelets with a profile of SDRs with a time history.
            void flameletSDRprofile(const std::vector< std::vector<double> >& sdr,
                                    const std::vector< std::vector<double> >& Zcoords,
                                    const std::vector<double>& intTime,
                                    bool continuation=false);

            //! Interface for use when soot volume fraction is available.
            void flameletWithSoot(const std::vector<double>& soot_fv, const std::vector<double>& sdr,
                                  const std::vector<double>& intTime, bool continuation=false);

            //! Assign the species mass fractions, temperature, density, and velocity into the vector of mixture objects.
            void resetMixtures(std::vector<Thermo::Mixture>& cstrs);

            //! Solve the reactor for the successive calls.
            void solve(std::vector<Thermo::Mixture>& cstrs,                      //Sprog mixture for each cell
                       const std::vector<double>& dz,                                   //geomtry cell widths
                       const std::vector< std::vector<double> >& initalSource,       //initial source terms ( under relaxed )
                       const std::vector< std::vector<double> >& finalSource,         //final source terms (under relaxed)
                       CamControl& ccObj,                                                    //Solver control
                       CamConfiguration& confObj,                                   //configuration object
                       Mechanism& mech_in,                                             //mechanism object
                       void* reactorModel,                                                  //reactor model to solve
                       const double sdr = 0);                                        //scalar dissipation rate in case of flamelets

            //! Return a vector of species names.
            std::vector<std::string> getSpeciesNames();

            //! Return the number of species.
            int getNumberOfSpecies() const;

            //! Return the number of moments.
            int getNumberOfMoments() const;

            //! Return the number of reactions.
            int getNumberOfReactions() const;

            //! Return the species mass fraction given the independent variable.
            double getMassFrac(const int spIndex, const double axpos);

            //! Return the soot moment given the independent variable.
            double getMoment(const int momIndex, const double axpos);

            //! Return the soot moment rates given the independent variable.
            // ank25 added for ELFM
            double getMomentWdot(const int momIndex, const double axpos);

            //! Return the species mole fraction given the independent variable.
            double getMoleFrac(const int spIndex, const double axpos);

            //! Get spatial profile of one species.
            std::vector<double> getMassFracsBySpecies(const int spIndex) const;

            //! Get spatial profile of one moment
            std::vector<double> getMomentsByIndex(const int momentIndex) const;

            //! Get spatial profile of rate of change of one moment
            // ank25 added for ELFM
            std::vector<double> getMomentsWdotByIndex(const int momentWdotIndex) const;

            //! Get mass fractions for all species at one point.
            std::vector<double> getMassFracsByPoint(const int indVarIndex) const;

            //! Return the temperature given the independent variable.
            double getTemperature(const double axpos);

            //! Return a vector of temperatures at all the independent variable points.
            const std::vector<double>& getTemperatures() const {return TVector;}

            //! Return the density given the independent variable.
            double getDensity(const double axpos);

            //! Return the viscosity given the independent variable.
            double getViscosity(const double axpos);

            //! Return the specific heat given the independent variable.
            double getSpecificHeat(const double axPos);

            //! return the sootAverageDiameter
            double getSootAverageDiameter(const double axpos);

            //!return the sootDispersion
            double getSootDispersion(const double axpos);

             //!return the sootSurfaceArea
            double getSootSurfaceArea(const double axpos);

             //!return the sootVolumeFraction
            double getSootVolumeFraction(const double axpos);

            //! Return the thermal conductivity given the independent variable.
            double getThermalConductivity(const double axPos);

            //! Return a vector of diffusion coefficents (size equals number of gas phase species).
            std::vector<double> getDiffusionCoefficients(const double axPos);

            //! Return pyrene production rate.
            double getWdotA4(const double axPos);

            //! Return the stoichiometric mixture fraction.
            double getStMixtureFrac();

            //! Return a vector of densities at all the independent variable points.
            const std::vector<double>& getDensities() const {return rhoVector;}

            //! Return a vector of the independent variable values.
            const std::vector<double>& getIndepVars() const {return indVar;}

            //! Return a vector of velocities at all the independent variable points.
            const std::vector<double>& getVelocities() const {return mVelocity;}


            //! Return the controller object.
            CamControl& getCamControl();

            //! Return the geometry object.
            CamGeometry& getCamGeometry();

            //! Return CamAdmin object.
            CamAdmin& getCamAdmin();

            //! Return the boundary object which stores the boundary condition information.
            CamBoundary& getCamBoundary();

            //! Return the profile object which contains the information of user defined temperature profile.
            CamProfile& getCamProfile();

            //! Return the configuration that stores the reactor type.
            CamConfiguration& getCamConfiguration();


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
            Sprog::Mechanism mech;
            int nSpecies;
            int nMoments;
            std::vector<std::string> speciesNames;
            std::vector<std::string> momentNames;
            const SpeciesPtrVector *speciesPointerVector;
            /*
             *memberes to hold dependent variabes
             */
            Array2D spMassFracs;            //species mass fractions
            Array2D sootMoments;			// Soot moments
            Array2D sootMomentsWdot;		// ank25 added for ELFM
            std::vector<double> TVector;     //temperature
            std::vector<double> rhoVector;   //density
            std::vector<double> muVector;    //viscosity
            std::vector<double> indVar;      //independant variable
            std::vector<double> spHeat;      //specific heats
            std::vector<double> lambda;      //thermal conductivity
            std::vector<double> mVelocity;   //velocity
            std::vector<double> avgMolWtVector;   // Average Molar Weight of the mixture
            std::vector<double> wdotA4;   // rate of production of pyrene
            std::vector<double> sootAverageDiameterVector;
            std::vector<double> sootDispersionVector;
            std::vector<double> sootSurfaceAreaVector;
            std::vector<double> sootVolumeFractionVector;

            Array2D mDiff;                  //Diffusion coefficients

            double stMixtureFrac;       //stoichiometric mixture fraction

            //! Pointer to flamelet class.
            FlameLet* flmlt;

            //! Calling interface to solve a flamelet with a given scalar dissipation rate.
            void flamelet(double sdr, double intTime=0, bool continuation=false);

            //! Function to return a variable's value given the independent variable.
            double getVariableAt(const double& pos, const std::vector<double>& var) const;

            void getFlameletVariables(FlameLet* const flmlt);

        }; // End Interface class declaration.

} // End Camflow namespace.


#endif	/* _INTERFACE_H */

