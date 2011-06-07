/*
 * File:   cam_soot.h
 * Author: vj231
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class implements the soot moment (MOMIC) based on
 *  chemkin implementation. In this file everything is
 *  evaluated in cgs unit
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
 * Created on 19 May 2009, 10:38
 */

#ifndef _CAM_SOOT_H
#define	_CAM_SOOT_H

#include <algorithm>
#include <vector>
#include <cmath>
#include <map>

#include "array.h"

#include "gpc.h"
#include "cam_math.h"
#include "comostrings.h"
#include "cam_params.h"
#include "cam_converter.h"
#include "camxml.h"

namespace Camflow{
    class CamSoot{
    public:

        enum regime{
            FM,
            CT,
            TR
        };
        CamSoot();
        ~CamSoot(){}
        void readSoot(CamConverter& convert,
                                            const CamXML::Element& node);
        void setSootSpecies(std::vector<std::string> species);
        void setNucCutoff(doublereal cutOff);
        void setPAHDia(doublereal dPAH);
        void setNumCAtomInception(int nCAtom);
        void setSootMomentActive();
        void setInceptionSpecies(std::string species);
        void setRegime(int n);
        void setFirstMom(doublereal m0);
        int getRegime();
        void initialize(int nCells, Sprog::Mechanism &mech, std::vector<doublereal> &mSolnVec);
        bool active();
        void setNumMoments(int n);
        int getNumMoments() const;
        doublereal nucleationConst(doublereal& T){
            return constNucleation*std::sqrt(T);
        }
        doublereal const_coag_freeMol(doublereal& T){
            return Kf*std::sqrt(T);
        }
        /*
         *grid function
         */
        doublereal grid(int k, int n, int m);
        doublereal betaC1(int i, int j);
        doublereal betaC2(int i, int j);
        void linear(int n, std::vector<doublereal>& y, doublereal& a, doublereal& b,
                    doublereal& rsq);
        void setSizeMoments(std::vector<doublereal>& conc);
        void sums(int hMoment, doublereal massAdded, doublereal coeff,
                            std::vector<doublereal>& rates);
        /*
         *soot reaction processes
         */
        void sootReactions(int cell, std::vector<doublereal>& conc, std::vector<doublereal>& mom,
                    int nSpec, doublereal T,   doublereal p);
        /*
         *nucleation rate
         */
        void nucleation(std::vector<doublereal> &conc,
                                    doublereal T,
                                    doublereal p);
        void coagulation(doublereal T,
                         doublereal p,
                        doublereal M0);

        void condensation(doublereal T,            //temperature
                    doublereal M0,           //zeroth moment
                    doublereal concPAH,      //concentration of PAH
                    doublereal ratePAH);     //rate of PAH consumption due to condensation
        void surface(int hMoment, int cell,
                    doublereal T,
                    doublereal M0,
                    std::vector<doublereal>& conc,
                    std::vector<doublereal>& totalRates);


        void momentResidual(const doublereal& time, int iMesh_s, int iMesh_e,
                                        int nVar, int nSpec,
                                        std::vector<doublereal>& dz,
                                        std::vector<doublereal>& u,
                                        std::vector<doublereal>& rho,
                                        doublereal* y, doublereal* f);

        void momentResidual(const doublereal& time, int iMesh_s, int iMesh_e,
                                        std::vector<doublereal>& dz,
                                        std::vector<doublereal>& u,
                                        std::vector<doublereal>& rho,
                                        doublereal* y, doublereal* f);

        void clearRates(int ncells);

        /*
         *adjust the species production rates for soot formation
         */
        void addRates(int nCells, Array2D& wdot);

        /*
         *report the results
         */
        void report(int nCells);

        /*-------------------------------------------                           
         *
         * NEW FUNCTIONS ADDED
         *                  
         *---------------------------------------------/
         */
        
        /*
         *Initialize the soot moments
         */
        void initMoments(Sprog::Mechanism &mech, std::vector<doublereal>& soln,int nCells=0);
        /*
         *calculate the residual functions for moments
         */
        void residual(const doublereal& time, std::vector<doublereal>& wdot,
                                doublereal* y, doublereal* f );
        /*
         *cauclate the all rates
         */
        void rateAll(std::vector<doublereal>& conc,
                     std::vector<doublereal>& moments,
                     doublereal& T,
                     doublereal& p,
                     std::vector<doublereal>& rate,
                     int cellID=0);
        /*
         *calculate nucleation rate
         */
        void rateNucleation(doublereal& concPAH,         //conc of PAH 4 nucln
                            doublereal& T,               //temperature
                            std::vector<doublereal>& nucRates //rates returned
                            );
        /*
         *calculate coagulation rate
         */
        void rateCoagulation(std::vector<doublereal>& mom,     //vector of moments
                            doublereal& T,                //temperature
                            std::vector<doublereal>& coagRates //rates returned
                            );
        /*
         *calculate condensation rate
         */
        doublereal rateCondensation(std::vector<doublereal>& mom, //vector of moments
                            doublereal& T,                   //temperature
                            doublereal& conc,                //PAH concentration
                            std::vector<doublereal>& cdRates      //rates returned
                            );
        /*
         *cauculate surface reaction rates
         */
        void rateSurface(std::vector<doublereal>& conc,          //concentration
                            doublereal T,                   //temperature
                            std::vector<doublereal>& mom,        //moment vector
                            std::vector<doublereal>& prodRates,  //surf prod rates
                            std::vector<doublereal>& totalRates);//rates returned

    private:
        const doublereal cMass;
        const doublereal rhoSoot;
        const doublereal ohMass;
        const doublereal lambda;
        const int atomsPerDimer;
        int numCAtomInception;
        int nMoments;
        int lowFrac, highFrac;
        int tRegime;
        /*
         * species index for surface production rates
         */
        int iCO, iH, iH2, iH2O, iO2, iOH, iC2H2, iInception;
        bool momentON;
        std::string iSpecies;
        doublereal particleDia;
        doublereal constNucleation;
        doublereal betaFree;
        doublereal betaContinuum;
        doublereal Kf,Kc_;
        doublereal kOH;
        doublereal kSurf;
        doublereal kPAH;
        doublereal dia_PAH;
        doublereal nucRateCutOff;
        doublereal D1_PAH;
        doublereal firstMom;

        std::vector<std::string> sootSpecies;
        std::vector<int> sootSpecIndex;
        std::vector<doublereal> base;

        Array2D bnCoeff,powPAH, prime;
        Array2D wdot,surfProdRate;
        Array2D conc_received, momReceived;

        static Array1D sizeMoments, reducedMoments;


        /*
         *members for storing the rates
         */
        //rates for nucleation, coagulation and condensation
        std::vector<doublereal> nucRate, cgRate, cdRate;
        //map<string,doublereal> surfProdRate; //surface production rate
        std::map<std::string, std::vector<doublereal> > smRates; //moment rates due to surface reactions

        //CamReporter creport;


        /*
         *New variables
         */
        doublereal Beta_nucl, Beta_fm, Beta_cd, Beta_surf;
        //int mmn2Calc;

        void interpolateReducedMoments(std::vector<doublereal>& wom);
        doublereal gridFunction(int k, int n, int m);


    };
}

#endif	/* _CAM_SOOT_H */

