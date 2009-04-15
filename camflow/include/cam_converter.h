/* 
 * File:   cam_converter.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  Unit conversion
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
 * Created on January 17, 2009, 2:40 PM
 */

#ifndef _CAM_CONVERTER_H
#define	_CAM_CONVERTER_H
#include "cam_params.h"
#include "comostrings.h"
#include "cam_error.h"
#include "gpc.h"
#include <string>
#include <vector>
using namespace Sprog;
namespace Camflow{
    class CamConverter{
        doublereal factor;
    public:
        CamConverter(){}
        ~CamConverter(){}
        //return the conversion factor based in the string argument
        doublereal getConvertionFactor(string unit);

        //return a string given an integer
        string getString(int n);

        //return a string given a doublereal
        string getString(doublereal a);

        //convert from mole to mass and returns tha average molecular weight
        doublereal mole2mass(vector<doublereal>& mole,
                             vector<doublereal>& mass,
                             Mechanism &mech );

        doublereal mass2mole(vector<doublereal>& mass,
                             vector<doublereal>& mole,
                             Mechanism &mech);
    };
}


#endif	/* _CAM_CONVERTER_H */

