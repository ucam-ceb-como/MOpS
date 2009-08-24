/* 
 * File:   cam_conc.h
 * Author: vinod(vj231@cam.ac.uk)
 *
 * Copyright (C) 2009 Vinod M Janardhanan
 *
 * File purpose:
 *  This class contains the implementation of inlet concentration
 *  definitions
 *
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
 * Created on January 18, 2009, 11:18 AM
 */

#ifndef _CAM_CONC_H
#define	_CAM_CONC_H

namespace Camflow{
    class CamConc{
    protected:
        int fracType;

    public:
        
        enum Conc{
            MASS,
            MOLE
        };
        CamConc(){}
        ~CamConc(){}
        CamConc(CamConc &cc);
        //set the fraction type mass or mole
        void setFracType(int n);

        //return the frac type such as mass or mole
        int getFracType() const;
    };
}

#endif	/* _CAM_CONC_H */

