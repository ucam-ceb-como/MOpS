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