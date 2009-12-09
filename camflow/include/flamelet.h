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
        void solve(vector<Thermo::Mixture>& cstrs,
                const vector< vector<doublereal> >& iniSource,
                const vector< vector<doublereal> >& fnlSource,
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
        void setExternalScalarDissipationRate(const vector<doublereal>& time, const vector<doublereal>& sdr);

        /**
         *  get scalar dissipation rate
         */
        
        const doublereal getSDR(const doublereal time);

        /**
         *  Set restart time
         */
        void setRestartTime(doublereal t);
 
        
        
    //! Computes the Planck mean absorption constants, as part of the radiative heat loss dissipation model.
    void PlanckAbsorption (const doublereal Temperature, doublereal Absorption[3])const;
    
    //! Computes  the radiative heat loss term for radiative heat dissipation model
    doublereal RadiativeLoss (const doublereal rho, const doublereal cp, const doublereal Temperature,
                              const doublereal SootVolFrac)const;


    private:
        doublereal stoichZ; //stoichiometric mixture fraction
        doublereal smr;     //stoichiometric mass ratio
        doublereal sdr;     // scalar dissipation rate
        doublereal sdr_ext; // scalar dissipation rate passed by exteranl program
        doublereal rstartTime;
        vector<doublereal> v_sdr;   //scalar dissipation rate that has a time history
        vector<doublereal> v_time; //time profile of scalar dissipation rates
        

        bool timeHistory;
        doublereal strain;  // strain rate
        int mCord;          // this is the mixture fraction coordinates
        inletStruct fuel, oxid;

        Array2D Le, convection, CpSpec; //Lewis numbers
        
    };
}

#endif	/* _FLAMELET_H */

