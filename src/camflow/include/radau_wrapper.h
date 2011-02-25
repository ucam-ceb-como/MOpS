/*
 * File:   RadauWrapper.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2009 Vinod M Janardhanan.
 *
 * File purpose:
 *  This is a wrapper class for Radau C++ solver
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
 * Created on January 24, 2009, 6:03 PM
 */

#ifndef _RADAUWRAPPER_H
#define	_RADAUWRAPPER_H
#include <vector>
#include "cam_control.h"
#include "StiffIntegratorT.h"
#include "cam_residual.h"
namespace Camflow{
    class RadauWrapper {

        CamControl *solverControl;
        StiffIntegratorT *solver;

        int bandWidth;

    public:

        RadauWrapper(){}
        ~RadauWrapper(){}

        //set the control parameters
        void setControl(CamControl &cc);

        //set the band width of Jacobian
        void setBandWidth(int n);


        //initialise the solver parameters
        void initSolver(int nEq,
                        const doublereal tBeg,
                        const doublereal tEnd,
                        std::vector<doublereal>& y,
                        std::vector<doublereal>& rTol,
                        std::vector<doublereal>& aTol,
                        CamResidual &cr);

        //integrate the DAE system
        void Integrate();
  
    };
}

#endif	/* _RADAUWRAPPER_H */

