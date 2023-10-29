/*
 * File:   cam_soot.h
 * Author: vj231
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *
 *  This class implements the soot moment (MOMIC) based on
 *  chemkin implementation. In this file everything is
 *  evaluated in cgs unit
 *
 * Licence:
 *
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

namespace Camflow {

/*!
 * CamSoot
 *
 * \brief Solves the MoMIC equations (Frenklach, 2002).
 *
 * \todo Get rid of all these Set functions, probably by incorporating into
 *       readSoot(). PARTIALLY DONE.
 * \todo Put physical processes such as coagulation in their own classes.
 * \todo Use abstract factory method instead of the enum regime.
 */
class CamSoot {

    typedef std::vector<double> realVector;
    typedef std::vector<std::string> stringVector;
    typedef std::vector<int> intVector;

    // Should get these from Sprog::Mechanism
    static const double cMass;
    static const double ohMass;

    public:

        /*!
         * Free Molecular (FM), Continuum (CT) or Transition (TR) Regime
         */
        enum regime {FM, CT, TR};

        /*!
         * Constructor.
         */
        CamSoot();

        /*!
         * Destructor.
         */
        ~CamSoot(){}

        /*!
         * Read in all the constants and settings.
         */
        void readSoot(CamConverter& convert, const CamXML::Element& node);

        //! Get the Regime (Free Molecular / Continuum / Transition).
        int getRegime() const;

        //! Get the number of moments.
        int getNumMoments() const;

        //! Get the user specified initial value of M0
        double getFirstMoment() const;

        //! Get number of atoms per diamer
        int getAtomsPerDiamer() const;

        //! Returns true if soot is going to be solved.
        bool active() const;

        void linear(int n, realVector& y,
                    double& a, double& b, double& rsq);

        void setSizeMoments(realVector& conc);

        void sums(int hMoment, double massAdded, double coeff,
                  realVector& rates);

        //! These are probably to be called from the flamelet code.
        // No longer needed. Flamelt code has its own residual
        void momentResidual(const double& time,
                            const int iMesh_s, const int iMesh_e,
                            const realVector& dz,
                            const realVector& u,
                            const realVector& rho,
                            const double* y, double* f,
                            const int nVar, const int nSpec);

        //! This just calls the other momentResidual function.
        void momentResidual(const double& time,
                            const int iMesh_s, const int iMesh_e,
                            const realVector& dz,
                            const realVector& u,
                            const realVector& rho,
                            const double* y, double* f);

        //void clearRates(int ncells);
        //void clearRates();

        /*!
         *adjust the species production rates for soot formation
         */
        void addRates(int nCells, Array2D& wdot);

        /*!
         *report the results
         */
        void report(int nCells);

        /*-------------------------------------------
         *
         * NEW FUNCTIONS ADDED
         *
         *---------------------------------------------/
         */

        /*!
         * Initialize the soot moments
         */
        //void initMoments(Sprog::Mechanism &mech, realVector& soln,int nCells=0);
        void initMoments(realVector& soln,int nCells);

        void initMomentsConstants(Mechanism &mech);

        /*!
         * Calculate the residual functions for moments. There does not seem to
         * be any point to this -> it just assigns wdot to f.
         */
        void residual(const double& time, realVector& wdot,
                                double* y, double* f );

        /*!
         * This calculates all of the moment rates and returns a vector of them.
         */
        realVector rateAll(const realVector& conc,
                const realVector& moments,
                const double& T,
                const double& p,
                const int cellID=0);

        realVector showGasPhaseRates(int nSpecies);

        realVector showSootComponentRates(int nMoments);

        double avgSootDiam();
        double dispersion();
        double sootSurfaceArea(double M0);
        double sootVolumeFraction(double M0);

    private:

        const double rhoSoot;
        const double lambda;
        int atomsPerDimer;
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
        double particleDia;
        double constNucleation;
        double betaFree;
        double betaContinuum;
        double Kf,Kc_;
        double kOH;
        double kSurf;
        double kPAH;
        double dia_PAH;
        double m0Threshold;
        double D1_PAH;
        double firstMom;
        double CD1;

        stringVector sootSpecies;

        realVector nucRates, coagRates, cdRates, sRates;	// Rate of change of moments

        Array2D bnCoeff,powPAH, prime;
        Array2D wdot; //,surfProdRate;
        realVector surfProdRate;
        Array2D conc_received, momReceived;

        static Array1D sizeMoments, reducedMoments;
        /*
         *members for storing the rates
         */
        //rates for nucleation, coagulation and condensation
        // Better to have these as local to rates All.
        // realVector nucRate, cgRate, cdRate;

        //map<string,double> surfProdRate; //surface production rate
        std::map<std::string, realVector > smRates; //moment rates due to surface reactions

        /*
         *New variables
         */
        // Beta_surf in CGS units. (as surface chem done in CGS)
        double Beta_nucl, Beta_fm, Beta_cd, Beta_surf_CGS, CBOHCGS;

        /*
         *calculate nucleation rate
         */
        realVector rateNucleation(const double& concPAH,  // conc of PAH 4 nucln
                                  const double& T);       // temperature

        /*
         * Calculate coagulation rate.
         * Seems to follow Eq 19. in Frenklach (2002).
         */
        realVector rateCoagulation(const realVector& mom,     // vector of moments
                                   const double& T);      // temperature

        /*
         *calculate condensation rate
         */
        realVector rateCondensation(const realVector& mom, //vector of moments
        		const double& T,           //temperature
        		const double& concPAH     //PAH concentration
                            );

        /*
         *calculate surface reaction rates
         */
        void rateSurface(const realVector& conc,       //concentration
                            const double& T,           //temperature
                            const realVector& mom,     //moment vector
                            realVector& prodRates,     //surf prod rates
                            realVector& sRates);       //rates returned

        /*
         *calculate surface reaction rates using HACARC 
         */
        void rateSurfaceHACARC(const realVector& conc, //concentration
                            const double& T,           //temperature
                            const realVector& mom,     //moment vector
                            realVector& prodRates,     //surf prod rates
                            realVector& sRates);       //rates returned



        void interpolateReducedMoments(realVector& wom);

        /*!
         * The grid function.
         * k = , n = , m =
         */
        double gridFunction(const int k, const int n, const int m);


}; // End Class CamSoot

} // End namespace Camflow

#endif	/* _CAM_SOOT_H */

