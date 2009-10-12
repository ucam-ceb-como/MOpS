/* 
 * File:   batch.h
 * Author: vj231
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class implements the batch reactor model
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
 * Created on 04 June 2009, 11:43
 */

#ifndef _BATCH_H
#define	_BATCH_H

#include "cam_residual.h"
#include "cam_control.h"
#include "cam_admin.h"
#include "cam_reporter.h"
#include "gpc.h"
#include "cam_setup.h"
#include "cam_configuration.h"
using namespace Sprog;
namespace Camflow{
    class Batch : public CamSetup{

    public:
        enum type{
            C_VOLUME,
            C_PRES
        };
        /*
         *default constructor
         */
        Batch(){};
        /*
         *virtual destructor
         */
        virtual ~Batch(){};
        /*
         *set the reactor type
         */
        void setType(int n);
        /*
         *get the reactor type
         */
        int getType();
        /*
         *function called by solver
         */
        int eval(doublereal x, doublereal* y, doublereal* ydot, bool jacEval);
        //console output
        void report(doublereal x, doublereal* solution);
        //console output with residuals
        void report(doublereal x, doublereal* solution, doublereal& res);
        void reportToFile(doublereal time, doublereal* solution);


        //solve
//        void solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg, CamProfile&cp,
//                     CamSoot &cs,Mechanism &mech);
        void solve(CamControl &cc, CamAdmin &ca, CamGeometry &cg,CamProfile &cp,
             CamConfiguration &config, CamSoot &cs,  Mechanism &mech );

        //return the initial solution vector
        void getInitial(vector<doublereal>& initial);

        //residual function definition----------------
        void residual(const doublereal& time, doublereal *y, doublereal *f);
        //species residual
        void speciesResidual(const doublereal& time, doublereal *y, doublereal *f);
        
        //temperature
        void energyResidual(const doublereal& time, doublereal *y, doublereal *f);
        
        //update the mixture properties
        void updateMixture(const doublereal& time, doublereal *y);

        //header information
        void header();
        
    private:
        

        int batchType;
        vector<doublereal> momRates;

    };
}

#endif	/* _BATCH_H */

