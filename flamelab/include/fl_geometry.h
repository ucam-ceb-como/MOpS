/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	Geometry object contains information on the reactor geometry and the
	descretization information. Each reactor object inherits from Geometry 
	class
  Licence:
    This file is part of "flameLab".

    flameLab is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/
#ifndef FL_GEOMETRY_H
#define FL_GEOMETRY_H
#include "fl_params.h"
#include "fl_error_handler.h"
#include<map>
namespace FlameLab{

	class Geometry{

		int nCell; // number of max finite volume cells
		int axialPosition;
		real rLength; // length of the reactor
		//real aspectRatio. First map contains from-to data and the
		//second map contains number of cells and aspect ratio
		map< map<real,real>, real > aspectRatioMap;
		vector<int> numCells;
		
	public:
		Geometry():nCell(10){// default number of max u,ber of cells			
			dz.resize(nCell);
		}

		Geometry(int n);		
		~Geometry();
		void setnCells(int n);
		void setLength(real len);
		void setAspectRatio(map< map<real,real>,real> fromTo, 
							vector<int> numCells);
		void setAspectRatio(vector<real> from,
							vector<real> to,
							vector<int> numCells,
							vector<real> ar);
		//return the number of cells
		int getnCells() const;

		//return the length of the reactor
		real getLength() const;

		//do the descretization
		void descretize();

		// return the descretizes geometry
		vector<real> getGeometry() const;

		//real getAspectRatio() const;
		void setAxialPosition(int n);

		// return the cell id corresponding to the axial position
		int getAxialPosition()const;

		//return the axial position in m
		real getAxialPosition(int n) const;

	protected:
		vector<real> dz; // vector holding descretization info
		typedef struct{
			real from;
			real to;
			int numCells;
			real ar;
		}aspectRatioStruct;
		
		vector<aspectRatioStruct> data_AspectRatio;

				


	};
};

#endif
