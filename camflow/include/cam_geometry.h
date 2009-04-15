/* 
 * File:   cam_geometry.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2008 Vinod M Janardhanan
 *
 *  File purpose:
 *  This class contains the implementation of geometry details
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
 * Created on January 17, 2009, 3:29 PM
 */

#ifndef _CAM_GEOMETRY_H
#define	_CAM_GEOMETRY_H
#include "cam_params.h"
#include "cam_error.h"
#include <vector>
#include <map>
namespace Camflow{
    class CamGeometry{
        doublereal length;                    //length of the model geometry
        doublereal dia;                 //rector diameter
        
        int nCell;                      //total number of FV cells
        vector<doublereal> dz, axPos;                //cell width
        string gridFile;

    public:
        CamGeometry(){
            length = 0;
            nCell = 0;
        }
        ~CamGeometry(){}


        void setGridFile(string name);
        //descretise the geometry. This will set the length 
        //of the reactor the total number of cells
        void descretize();

        //set the reactor length
        void setLength(doublereal len);
        //set the reactor diameter
        void setDia(doublereal dia);

        //return the reactor dia
        doublereal getDia() const;

        //return the cross sectional area
        doublereal getArea() const;
        
        //return the surface area
        doublereal getSurfArea();

        //return the surface area per unit length
        doublereal getSurfAres_l() const;

        //return the length of the model geometry
        doublereal getLenth() const;

        //return the total number of cells
        int getnCells() const;

        //reuturn the geometry info
        vector<doublereal>& getGeometry();

        //return the axial position vector
        vector<doublereal>& getAxpos();
        //add cells with zero width
        void addZeroWidthCells();
    };
}

#endif	/* _CAM_GEOMETRY_H */

