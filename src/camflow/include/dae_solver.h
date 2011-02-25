/* 
 * File:   dae_solver.h
 * Author: vj231
 *
 * Copyright (C) 2009 Vinod M Janardhanan.
 *
 * File purpose:
 *  base class for DAE solvers
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
 * Created on 09 February 2009, 12:56
 */

#ifndef _DAE_SOLVER_H
#define	_DAE_SOLVER_H
#include "cam_residual.h"
//#include "cam_control.h"
namespace Camflow{

    class DAESolver{
    protected:
        //CamResidual& mResidual;
    public:
        //DAESolver(CamResidual& cr):mResidual(cr){}
        DAESolver(){}
        ~DAESolver(){}
        virtual void setBand(){
            cout << "Base class: nothing to do\n";
        }
        virtual void setDens(){
            cout << "Base class: nothing to do\n";
        }
    };

}


#endif	/* _DAE_SOLVER_H */

