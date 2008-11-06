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

#include "fl_geometry.h"
using namespace FlameLab;


// constructor
Geometry::Geometry(int n){
	nCell = n;
	dz.resize(nCell);
}
// Geometry desctructor
Geometry::~Geometry(){
	
}

// sets the number of finite volume cells
void Geometry::setnCells(int n){
	nCell = n;	
	descretize();
}
// set the reactor length in case of premix and 
//sets the nozzle gaps in case of counter flow flames
void Geometry::setLength(real len){
	rLength = len;
}
// sets the grid aspect ratio
void Geometry::setAspectRatio(real ar){
	aspectRatio = ar;
}

// descretize the geometry in to nCells
void Geometry::descretize(){
	dz.resize(nCell);
	for(int i=0; i<nCell; i++)
		dz[i] = rLength/nCell;
}

// returns the number of finite volume cells
int Geometry::getnCells() const{
	return nCell;
}
// returns the length of the reactor in (m)
real Geometry::getLength() const{
	return this->rLength;
}
// returns the grid aspect ratio
real Geometry::getAspectRatio() const{
	return this->aspectRatio;
}

// returns a vector containig the width of finite volume cells
vector<real> Geometry::getGeometry() const{
	return this->dz;
}
// set the axial position
void Geometry::setAxialPosition(int n){
	this->axialPosition = n;
}

// returns the axial position

int Geometry::getAxialPosition() const{
	return this->axialPosition;
}