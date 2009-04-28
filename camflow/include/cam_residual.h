/* 
 * File:   cam_residual.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2009 Vinod M Janardhanan.
 *
 * File purpose:
 *  This is the base class for residual function definitions
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
 * Created on January 24, 2009, 7:57 PM
 */

#ifndef _CAM_RESIDUAL_H
#define	_CAM_RESIDUAL_H
#include "cam_params.h"
#include "cam_converter.h"
#include "cam_boundary.h"
#include "cam_admin.h"
#include "cam_geometry.h"
#include "array.h"
#include "gpc.h"
#include <map>
#include <vector>
using namespace Sprog;
using namespace std;
namespace Camflow{
    class CamReporter;
    class CamResidual{

    public:

        enum EqnType{
            EQN_ALL,
            EQN_SPECIES,
            EQN_ENERGY,
            EQN_MASSFLOW
        };

        CamResidual(){};
        virtual ~CamResidual(){}

        //residual evaluation function
        virtual int eval(doublereal t, doublereal* y, doublereal* ydot, bool jacEval)=0;
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
         *mass matrix definition for Radau
         */
        virtual void massMatrix(doublereal **M)=0;
        //get intitial
        virtual void getInitial(vector<doublereal> &initial)=0;
        
        //base class definition for species residuals
        virtual void speciesResidual(const doublereal& time, doublereal* y, doublereal* f);

        // base class definition for energy residuals
        virtual void energyResidual(const doublereal& time, doublereal* y, doublereal* f);

        //base class definition for mass flow
        virtual void massFlowResidual(const doublereal& time, doublereal* y, doublereal* f);
        
        //stores the mixture properties for the calculation of fluxes
        virtual void saveMixtureProp(doublereal* y, bool thermo);
        

        //update the thermal flux. The flux is stored on the west cell
        //face for any given cell index. For eg the flux
        //returned for cell index 2 is the interface flux between 1 and 2
        virtual void updateThermo();

        //Calculates the diffusion fluxes and stores in cr_jk.
        //The fluxes are stored at for the left cell face (w), i.e
        //the flux stored at cell index 2 is the flux at the interface
        //between cell index 1 and 2. The fluxes are calculated based
        //on mixture averaged approach
        void updateDiffusionFluxes();

        //return the number of species
        const int& getNSpecies() const ;
        //return the number of variables
        const int& getNVar() const;
        //return the total number of equations
        const int& getNEqn() const;
        //derivative calculation
        doublereal dydx(doublereal nr1, doublereal nr2, doublereal dr){return ((nr1-nr2)/dr);}

        /*
         *extract various dependent variabes
         */
        void extractMassFlowVector(vector<doublereal>& vec);
        void extractSpeciesVector(vector<doublereal>& vec);
        void extractEnergyVector(vector<doublereal>& vec);

        /*
         *merge
         */
        void mergeMassFlowVector(doublereal* vec);
        void mergeSpeciesVector(doublereal* vec);
        void mergeEnergyVector(doublereal* vec);


        

    protected:
        
        
        Array2D s_H;
        Array2D s_mf;
        Array2D s_Diff;
        Array2D s_Wdot;        
        Array2D s_jk;
        Array2D s_cp;
        
        vector<doublereal> m_T;                  //mixture temperature
        vector<doublereal> m_cp;                 //mixture specific heat
        vector<doublereal> m_k;                  //mixture thermal conductivity
        vector<doublereal> m_rho;                //mixture density
        vector<doublereal> m_u;                  //fluid velocity
        vector<doublereal> m_q;                  //thermal conduction flux
        vector<doublereal> m_flow;               //mass flow
        vector<doublereal> wdot;                 //rate of production
        vector<doublereal> dz;                   //grid spacting
        vector<doublereal> axpos;                //axial position

        vector<doublereal> resSp, resT, resM;
       
        doublereal opPre;                        //operating pressure

        int nEqn;    //number of equations
        int nVar;    //number of variables
        int nSpc;    //number of species
        int ptrT;    //array offset to temperature
        int ptrC;    //array offset to continuity
        int ptrM;    //array offset to momentum
        int ptrS;    //array offset to species
        int ptrR;    //aray offset for residence time
        int loopBegin, cellBegin;//first cell
        int loopEnd, cellEnd;//last cell
        int eqn_slvd;

        Sprog::Thermo::Mixture *camMixture;
        Sprog::Mechanism *camMech;
        const SpeciesPtrVector *spv;
       // doublereal* yData;

        //members for the reactor models
        vector<doublereal> solvect, rTol, aTol;
        doublereal vel,rho, Ac, As;
        vector<string> headerData;
        CamAdmin *admin;
        CamReporter *reporter;
        CamGeometry *reacGeom;

    };

}

#endif	/* _CAM_RESIDUAL_H */

