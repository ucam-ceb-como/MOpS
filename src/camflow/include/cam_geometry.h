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

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "cam_params.h"
#include "cam_error.h"
#include "comostrings.h"

namespace Camflow{

    class CamGeometry
    {

        double length;                   //length of the model geometry
        double dia;                      //rector diameter

        int nCell;                           //total number of FV cells
        std::vector<double> dz, axPos;        //cell width
        std::string gridFile;

        /*
         *grid refinement, not well implemented in the present version
         */
        double a_slope;
        double curve;
        double minRange;
        double mPrune;
        std::map<int, int> z_loc, z_keep;

    public:

        CamGeometry();

        ~CamGeometry(){}


        void setGridFile(std::string name);
        //discretise the geometry. This will set the length
        //of the reactor the total number of cells
        //void discretize();

        //set the reactor length
        void setLength(double len);
        //set the reactor diameter
        void setDia(double dia);

        //get the name of the grid input file.
        std::string getGridFileName() const;

        //return the reactor dia
        double getDia() const;

        //return the cross sectional area
        double getArea() const;

        //return the surface area
        double getSurfArea();

        //return the surface area per unit length
        double getSurfAres_l() const;

        //return the length of the model geometry
        double getLenth() const;

        //return the total number of cells
        int getnCells() const;

        //reuturn the geometry info
        const std::vector<double>& getGeometry() const ;

        //return the axial position vector
        const std::vector<double>& getAxpos() const ;
        //add cells with zero width
        void addZeroWidthCells();

        /*
         * set the geometry information. this will set the vector of
         *cell widths into dz vctor
         */
        void setGeometry(const std::vector<double>& dz_);

        //refine the grid
        void refine(double* y, const int nVar, const int nSpec, int ptrT=0);

    };
}

#endif	/* _CAM_GEOMETRY_H */

