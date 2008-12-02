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
	//descretize();
}
// set the reactor length in case of premix and 
//sets the nozzle gaps in case of counter flow flames
void Geometry::setLength(real len){
	rLength = len;
}
// sets the grid aspect ratio
void Geometry::setAspectRatio(map< map<real,real>,real> fromTo, 							  
							  vector<int>ncells)
{
	aspectRatioMap=fromTo;
	numCells = ncells;
	descretize();


}

// descretize the geometry in to nCells.
// This function descretizes the geometry based on aspect ratios
// Aspect ratios can be > or < 1. If the AR is > 1, then  fine 
// grids are generated at the beginning of the reactor. If AR<1
// fine grifs are generated at the exit of the reactor. For counter
// flow diffusion flames it is ideal to define 2 ARs with one AR>1 and
// the other with AR<1
void Geometry::descretize(){
	//get the size of aspect ratio map
	int arSize = aspectRatioMap.size();
	if(arSize > 2 )
		throw ErrorHandler("More than 2 section not supported\n",4004);
	map< map<real,real>, real>::iterator p; //iterator for from-to value paies
	real totalLength=0;
	vector<real> length, aspectRatio;
	nCell = 0;
	int counter = 0;
	for(p=aspectRatioMap.begin(); p!= aspectRatioMap.end(); ++p){
		map<real,real> fromTo = p->first; // get the from-to key from the map
		map<real,real>::iterator ftp= fromTo.begin(); //iterate from-key to strip from and to
		real from = ftp->first;
		real to = ftp->second;
		totalLength += (to-from); // sum of aspect ratio segments should not exceed the defined reactor length
		length.push_back(to-from);
		aspectRatio.push_back(p->second);
		// this over writes any prescribed finite volume cells
		nCell += numCells[counter]; 
		counter ++;
		
	}

	if(totalLength != rLength) throw ErrorHandler("Aspect ratio in-consistancy\n",606);
	dz.resize(nCell,0.0);
	for(unsigned int j=0; j<length.size(); j++){

		real sum = 0;
		real d =1.0;
		real aspect;
		/**
		** implementation is same for ARs >&< 1. The directions are chosen
		** while assignment to dz
		**/
		if(aspectRatio[j] <1) 
			aspect = 1+aspectRatio[j];
		else
			aspect = aspectRatio[j];

		for(int i=0; i<numCells[j]; i++){
			sum += d;
			d *= aspect;
		}
		d = length[j]/sum;
		int pos,cellBegin,cellEnd;
		if(j==0){
			cellBegin = 0; 
			cellEnd = numCells[j];
		}else{
			cellBegin = numCells[j-1];
			cellEnd = cellBegin+numCells[j];
		}
		for(int i=cellBegin; i<cellEnd; i++){
			pos = i;
			if(aspectRatio[j]<1) // reverse the direction
				pos=numCells[j]-1-i;
			dz[pos] = d;
			d *= aspect;
		}
	
	}
				
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
//real Geometry::getAspectRatio() const{
//	return this->aspectRatio;
//}

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

real Geometry::getAxialPosition(int n) const{
	real pos=0.0;
	int i;
	if(n==0) return dz[0];
	for(i=0; i<n; i++)
		pos += dz[i];
	return pos;
}