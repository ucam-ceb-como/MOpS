/* 
 * File:   flamelet.h
 * Author: vj231
 *
 * Created on 06 July 2009, 12:13
 */

#ifndef _FLAMELET_H
#define	_FLAMELET_H

#include "cam_residual.h"
#include "cam_control.h"
#include "cam_admin.h"
#include "cam_reporter.h"
#include "gpc.h"
#include "cam_setup.h"

using namespace Sprog;

namespace Camflow{
    class FlameLet: public CamSetup{
    public:

        enum{
            LNONE,
            LNNONE
        };
        
        FlameLet(){sdr_ext=0;timeHistory=false;}
        virtual ~FlameLet(){}

        int eval(doublereal x, doublereal* y, doublereal* ydot, bool jacEval);
        //console output
        void report(doublereal x, doublereal* solution);
        //console output with residuals
        void report(doublereal x, doublereal* solution, doublereal& res);
        //file output
        void reportToFile(doublereal x, doublereal* solution);

        /**
         * solve the flamelet: call for coupling without
         * solving population balance
         */
        void solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg, CamProfile&cp,
                            Mechanism &mech, bool interface=false);

        /**
         * stand alone call as well as first call from an external
         * code that solves the population balance
         */
        void solve(CamControl& cc, CamAdmin& ca, CamGeometry& cg, CamProfile& cp,
                CamConfiguration& config, CamSoot& cs, Mechanism& mech);

        /**
         * continulation call from an external code that
         * solves the population balance
         */
        void solve(std::vector<Thermo::Mixture>& cstrs,
                const std::vector< std::vector<doublereal> >& iniSource,
                const std::vector< std::vector<doublereal> >& fnlSource,
                Mechanism& mech,
                CamControl &cc,
                CamAdmin &ca,
                CamGeometry &cg,
                CamProfile& cp);


        
        //initialize the solution vector
        void initSolutionVector(CamControl &cc);
        //calculate the stoichiometric mixture fraction
        doublereal stoichiometricMixtureFraction();
        //calculate the scalar dissipation rate
        doublereal scalarDissipationRate(const doublereal m_frac);
        /*
         *create header for file output
         */
        void header();
        /*
         *coupled solver
         */
        void csolve(CamControl &cc, bool interface=false);
        /*
         *segregated solver
         */
        void ssolve(CamControl &cc);
        /*
         *restart the solution with the converged solution
         */
        void restart(CamControl &cc);
        /*
         *save the solution vector
         */
        void saveMixtureProp(doublereal* y);
        /*
         *residual defiitions
         */
        void residual(const doublereal& t, doublereal* y, doublereal* f);
        /*
         *species residual
         */
        void speciesResidual(const doublereal& t, doublereal* y, doublereal* f);
        /*
         *Energy residual
         */
        void energyResidual(const doublereal& t, doublereal* y, doublereal* f);
        /*
         *set the external scalar dissipation rate
         */
        void setExternalScalarDissipationRate(const doublereal sr);

        /*
         *set the time history of scalar dissipation rate
         */
        void setExternalScalarDissipationRate(const std::vector<doublereal>& time, const std::vector<doublereal>& sdr);

        /**
         *  get scalar dissipation rate
         */
        
        doublereal getSDR(const doublereal time) const;

        //! Provide a soot volume fraction from an external calculation
        void setExternalSootVolumeFraction(const std::vector<doublereal>& soot_fv);


        /**
         *  Set restart time
         */
        void setRestartTime(doublereal t);
 
        
        
    //! Computes the Planck mean absorption constants, as part of the radiative heat loss dissipation model.
    void PlanckAbsorption (const doublereal Temperature, doublereal Absorption[3])const;
    
    //! Computes  the radiative heat loss term for radiative heat dissipation model
   doublereal RadiativeLoss(const doublereal Temperature, 
                            const doublereal soot_vol_frac, const doublereal mole_frac_H2O, 
                            const doublereal mole_frac_CO2, const doublereal mole_frac_CO) const; 

    private:
        doublereal stoichZ; //stoichiometric mixture fraction
        doublereal smr;     //stoichiometric mass ratio
        doublereal sdr;     // scalar dissipation rate
        doublereal sdr_ext; // scalar dissipation rate passed by exteranl program
        doublereal rstartTime;
        std::vector<doublereal> v_sdr;   //scalar dissipation rate that has a time history
        std::vector<doublereal> v_time; //time profile of scalar dissipation rates

        //! Spatial profile of soot volume fraction
        std::vector<doublereal> m_SootFv;
        

        bool timeHistory;
        doublereal strain;  // strain rate
        int mCord;          // this is the mixture fraction coordinates
        inletStruct fuel, oxid;

        Array2D Le, convection, CpSpec; //Lewis numbers
        
    };
}

#endif	/* _FLAMELET_H */

