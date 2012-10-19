/* 
 * File:   cam_control.h
 * Author: vinod (vj231@cam.ac.uk)
 * 
 * Copyright (C) 2008 Vinod M Janardhanan.

 * File purpose:
 *  This class can be used to control the behavior of the solver
 *  based on the controls defined in the input file
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
 * Created on January 17, 2009, 1:00 PM
 */

#ifndef _CAM_CONTROL_H
#define	_CAM_CONTROL_H
#include "cam_params.h"
#include "cam_conc.h"
/*
 *
 */
namespace Camflow{
    class CamControl{
        double relTolSp, absTolSp;        //species rel and abs tols
        double relTolT, absTolT;          //temperature rel and abs tols
        double relTolFlow, absTolFlow;    //flow rel and abs tols
        double resTol;                    //residual tolerence
        double iniStep, maxStep,minStep;  //solver stepsize control
        double maxTime;                   //max integration time
        double urSpecies;                 //under relaxation for species sources
        bool resMonitor;
    protected:
        int solMode;                    //solution mode steady or trans
        int repotMode;                  //repot mode intermediate or final
        int solver;
        int nIter;                      //number of iterations for segregated solver

    public:
        enum Solution{                  //steady state or transient solution
            COUPLED,                    // coupled solver
            SEGREGATED                  // segregated solver
        };

        enum Report{                    //output control
            INTERMEDIATE,
            FINAL
        };
        /*
         *only CVode is supported at the moment. The wrapper classes
         *and their implementations may be found in the distribution.
         *However they are not included in the makefile. You may try other
         *solvers at own risk.
         */
        enum Solvers{
            RADAU,
            CVODE,
            KINSOL,
            IDA,
            NEWTON,
            LIMEX
        };

        //constructor
        CamControl();
        //distructor
        ~CamControl(){}

        //set the species rel tol
        void setSpeciesRelTol(double tol);

        //set the species abs tol
        void setSpeciesAbsTol(double tol);

        //set the temperature rel tol
        void setTempRelTol(double tol);

        //set temperature abs tol
        void setTempAbsTol(double tol);

        //set the flow rel tol
        void setFlowRelTol(double tol);

        //set the flow abs tol
        void setFlowAbsTol(double tol);

        //set the residual tolerence
        void setResTol(double tol);

        //set the solution option
        void setSolutionMode(int n);

        //set the report mode
        void setReportInterval(int n);

        //set initial step size
        void setIniStep(double step);

        //set the maximum allowed step size
        void setMaxStep(double step);

        //set the minimum allowed step
        void setMinStep(double step);

        //set max integration time
        void setMaxTime(double time);

        //set the solver
        void setSolver(int n);

        void setResidualMonitor(bool lopt);

        //set the under relaxation for the species
        void setSpeciesUnderRelax(double ur);

        //number of iterations
        void setNumIterations(int n);
        
        //return the species rel tol
        double getSpeciesRelTol() const;

        //return the species abs tol
        double getSpeciesAbsTol() const;

        //return the temperature rel tol
        double getTempRelTol() const;

        //return the temperature abs tol
        double getTempAbsTol() const;

        //return the flow rel tol
        double getFlowRelTol() const;

        //return the flow abs tol
        double getFlowAbsTol() const;

        //return the residual tolerence
        double getResTol() const;

        //return the solution option
        int getSolutionMode() const;

        //return the reporting mode
        int getReportInterval() const;

        //return the solver
        int getSolver() const;

        //return the initial step
        double getIniStep() const;

        //return the maximum step size
        double getMaxStep() const;

        //return the minimum step size
        double getMinStep() const;

        //return the max integration time
        double getMaxTime() const;

        bool getResidualMonitor() const;

        //return the number of iterations
        int getNumIterations() const;

        //return the species under relaxation factor
        double getSpeciesUnderRelax() const;

    };
}

#endif	/* _CAM_CONTROL_H */

