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
        FlameLet(){sdr_ext=0;}
        virtual ~FlameLet(){}

        int eval(doublereal x, doublereal* y, doublereal* ydot, bool jacEval);
        //console output
        void report(doublereal x, doublereal* solution);
        //console output with residuals
        void report(doublereal x, doublereal* solution, doublereal& res);
        //file output
        void reportToFile(doublereal x, doublereal* solution);
        //solve the flamelet
        void solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg, CamProfile&cp,
                            Mechanism &mech, bool interface=false);
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
        void setExternalScalarDissipationRate(doublereal sr);
        


    private:
        doublereal stoichZ; //stoichiometric mixture fraction
        doublereal smr;     //stoichiometric mass ratio
        doublereal sdr;     // scalar dissipation rate
        doublereal sdr_ext; // scalar dissipation rate passed by exteranl program
        doublereal strain;  // strain rate
        int mCord;          // this is the mixture fraction coordinates
        inletStruct fuel, oxid;
    };
}

#endif	/* _FLAMELET_H */

