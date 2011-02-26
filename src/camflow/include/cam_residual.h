/*!
 * \file   cam_residual.h
 * \author V. Janardhanan
 *
 * \brief This is the base class for residual function definitions.
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

#ifndef _CAM_RESIDUAL_H
#define	_CAM_RESIDUAL_H

#include "cam_params.h"
#include "cam_converter.h"
#include "cam_boundary.h"
#include "cam_admin.h"
#include "cam_geometry.h"
#include "cam_control.h"
#include "cam_configuration.h"
#include "array.h"
#include "gpc.h"
#include "cam_soot.h"
#include "cam_profile.h"
#include <map>
#include <vector>
#include <stdexcept>

using namespace Sprog;

namespace Camflow{

    class CamReporter;

    /*!
     *@brief    Base class from which all of the solver classes are derived.
     *
     * Include a more detailed description here.
     */
    class CamResidual{

        public:

            enum EqnType{
                EQN_ALL,
                EQN_SPECIES,
                EQN_ENERGY,
                EQN_MOMENTUM,
                EQN_CONTINUITY,
                EQN_MOMENTS,
                EQN_FLOW //solve momentum and continuity coupled
            };

            //! Default constructor.
            CamResidual();

            //! Destructor.
            virtual ~CamResidual();



            virtual void solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg,CamProfile &cp,
                 CamConfiguration &config, CamSoot &cs,  Mechanism &mech )=0;

            /*
             * residual evaluation function for ODE and DAEs
             */
            virtual int eval(doublereal t, doublereal* y, doublereal* ydot, bool jacEval)=0;

            //mass matrix evaluation
            virtual void massMatrix(doublereal **M);

            /*
             *console output for monitoring the integration process
             */
            virtual void report(doublereal t, doublereal* solution)=0;
            /*
             *console output for monitoring the integration process
             *with residual output
             */
            virtual void report(doublereal t, doublereal* solutio, doublereal& res)=0;

            /*
             *call to solve the reactor from the external interface
             */
            virtual  void solve(std::vector<Thermo::Mixture>& cstrs,
                    const std::vector< std::vector<doublereal> >& iniSource,
                    const std::vector< std::vector<doublereal> >& fnlSource,
                    Mechanism& mech,
                    CamControl &cc,
                    CamAdmin &ca,
                    CamGeometry &cg,
                    CamProfile& cp);


            /*
             * residual evaluation for newton solvers
             */
            virtual int eval(doublereal* y, doublereal* ydot);

            /*
             * base class definition for species residuals. This function can be used for
             * any boundary value problems. The boundary condition has to be implemented
             * in the respective reactor models
             */
            virtual void speciesResidual(const doublereal& time, doublereal* y, doublereal* f);

            /*
             * base class definition for energy residuals. This function can be used for
             * any boundary value problems. The boundary condition has to be implemented
             * in the respective reactor models
             */
            virtual void energyResidual(const doublereal& time, doublereal* y, doublereal* f, bool flxTrans=false);

            /*
             * base class definition for mass flow. This function can be used for
             * any boundary value problems. The boundary condition has to be implemented
             * in the respective reactor models
             */
            virtual void massFlowResidual(const doublereal& time, doublereal* y, doublereal* f);

            /*
             * flow field evaluation (algebraic equation system)
             */
            virtual void calcFlowField(const doublereal& time, doublereal* y);

            /**
             *set the scalar dissipation rate in case of flamelet odel
             */
            virtual void setExternalScalarDissipationRate(const doublereal sr);

            /*
             * stores the mixture properties for the calculation of fluxes. This
             * is normally done for each function call from the solver as the
             * species transport properties and thermal properties are
             * dependent on the species compositions
             */
            virtual void saveMixtureProp(const doublereal time, doublereal* y, bool thermo, bool mom);




            /*
             * update the thermal flux. The flux is stored on the west cell
             * face for any given cell index. For eg the flux
             * returned for cell index 2 is the interface flux between 1 and 2
             */
            virtual void updateThermo();

            /* Calculates the diffusion fluxes and stores in cr_jk.
             * The fluxes are stored at for the left cell face (w), i.e
             * the flux stored at cell index 2 is the flux at the interface
             * between cell index 1 and 2. The fluxes are calculated based
             * on mixture averaged approach
             */
            void updateDiffusionFluxes();

            /*
             * return the number of species
             */
            const int& getNSpecies() const ;
            /*
             * return the number of variables
             */
            const int& getNVar() const;
            /*
             * return the total number of equations
             */
            const int& getNEqn() const;

            /*
             *return the species mass fractions to the calling program
             */
            virtual void getSpeciesMassFracs(Array2D& mf);
            /*
             *return the average molar weight of the mixture to the calling program
             */
            virtual void getAverageMolarWeight(std::vector<doublereal>& avgMolWt_);
            /*
             *return the temperature
             */
            virtual void getTemperatureVector(std::vector<doublereal>& temp);
            /*
             *return the density vecot
             */
            virtual void getDensityVector(std::vector<doublereal>& density);
            /*
             *return the viscosity vector
             */
            virtual void getViscosityVector(std::vector<doublereal>& viscosity);
            /*
             *return point velocities
             */
            virtual void getVelocity(std::vector<doublereal>& vel);

            /*
             *  return the specific heat
             */
            virtual void getSpecificHeat(std::vector<doublereal>& spHeat);
            /*
             * return a 2D array for Diffusion coefficient
             */
            virtual void getDiffusionCoefficient(Array2D& dCoeff);
            /*
             * Return the thermal conductivity
             */
            virtual void getThermalConductivity(std::vector<doublereal>& lambda);
            /*
             *return the independant variable
             */
            virtual void getIndepedantVar(std::vector<doublereal>& indVar);

            //Get the species residual
            virtual doublereal getResidual() const;

            /*
             *store the particle source terms for use in species residual
             *by passing the initial and final source terms
             */
            virtual void setParticleSource(const std::vector< std::vector<doublereal> >& initial,
                   const  std::vector< std::vector<doublereal> >& final);


            //derivative calculation
            doublereal dydx(doublereal nr1, doublereal nr2, doublereal dr){
                return ((nr1-nr2)/dr);
            }

            doublereal dydx(const int i, doublereal phi_E, doublereal phi_P, doublereal phi_W,
                                                doublereal dr){
                doublereal val = m_u[i] > 0 ? (phi_P-phi_W)/dr : (phi_E-phi_P)/dr;
                return val;

            }

            /*
             *extract various dependent variabes from the master
             *solution vector
             */
            void extractContinuity(std::vector<doublereal>& vec);
            void extractSpeciesVector(std::vector<doublereal>& vec);
            void extractEnergyVector(std::vector<doublereal>& vec);
            void extractMomentum(std::vector<doublereal>& vec);
            void extractSootMoments(std::vector<doublereal>& vec);

            /*
             *merge the denepdent variables to the
             *master solution vector
             */
            void mergeContinuity(doublereal* vec);
            void mergeSpeciesVector(doublereal* vec);
            void mergeEnergyVector(doublereal* vec);
            void mergeMomentum(doublereal* vec);
            void mergeSootMoments(doublereal* vec);


            //The following functions are specifically for flamelets

            /**
             *  Set the lewis number to be used
             */
            void setLewisNumber(int n){
                Lewis = n;
            }

            /**
             *  Return the Lewis number flag
             */
            int getLewisNumber() const{
                return Lewis;
            }


        protected:


            Array2D s_H;
            Array2D s_mf;
            Array2D s_Diff;
            Array2D s_Wdot;
            Array2D s_ParticleBegin, s_ParticleEnd;
            Array2D s_jk;
            Array2D s_cp;

            std::vector<doublereal> m_T;                  //mixture temperature
            std::vector<doublereal> m_cp;                 //mixture specific heat
            std::vector<doublereal> m_k;                  //mixture thermal conductivity
            std::vector<doublereal> m_rho;                //mixture density
            std::vector<doublereal> m_u;                  //axial velocity
            std::vector<doublereal> m_v_grd;              //radial velocity gradient
            std::vector<doublereal> m_q;                  //thermal conduction flux
            std::vector<doublereal> m_flow;               //mass flow
            std::vector<doublereal> m_G;                  //radial velocity gradient
            std::vector<doublereal> m_mu;                 //viscosity
            std::vector<doublereal> m_shear;              //shear rate
            std::vector<doublereal> m_eigen;              //pressure gradient eigen value
            std::vector<doublereal> wdot;                 //rate of production
            std::vector<doublereal> dz;                   //grid spacting
            std::vector<doublereal> axpos;                //axial position
            std::vector<doublereal> avgMolWt;             //average molecular weight
            std::vector<doublereal> slopes;               //slopes of piece-wise linear particle sources
            std::vector<doublereal> radiation;            //radiative heat loss term for output to profile.h

            std::vector<doublereal> resSp, resT, resFlow, resMoment, resAxVel;


            doublereal opPre;                        //operating pressure
            doublereal Tignition;                   //Ignition temperature





            int nEqn;    //number of equations
            int nVar;    //number of variables
            int nMoments;//number of moments
            int nSpc;    //number of species
            int ptrT;    //array offset to temperature
            int ptrF;    //array offset to continuity
            int ptrG;    //array offset to momentum
            int ptrH;    // 1/r (dp/dr)
            int ptrS;    //array offset to species
            int ptrR;    //aray offset for residence time
            int iMesh_s, cellBegin;//first cell
            int iMesh_e, cellEnd;//last cell
            int eqn_slvd;

            int Lewis;


            Sprog::Thermo::Mixture *camMixture;
            Sprog::Mechanism *camMech;
            const SpeciesPtrVector *spv;
           // doublereal* yData;

            //members for the reactor models
            std::vector<doublereal> solvect, rTol, aTol;
            doublereal vel,rho, Ac, As;
            std::vector<std::string> headerData;
            CamAdmin *admin;
            CamReporter *reporter;
            CamGeometry *reacGeom;
            CamSoot *sootMom;
            CamControl *control;


    }; // End CamResidual class declaration.

} // End Camflow namespace.

#endif	/* _CAM_RESIDUAL_H */
