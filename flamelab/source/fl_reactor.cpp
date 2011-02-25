/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        flameLab (premix solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
	This is the reactor class. It contains information
	regarding the reactor under simulation. It holds information
	on the inlet condition, geometry of the reactor, and
	operating conditions defined for the reactor.
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

#include "fl_reactor.h"
#include "fl_premix.h"
#include <iostream>
using namespace FlameLab;
using namespace Strings;
using namespace Sprog;
using namespace std;

int Reactor::nSpecies;
const Sprog::Mechanism* Reactor::mech;

Reactor::Reactor(){
	fuel = NULL;
	oxidizer = NULL;
}

Reactor::~Reactor(){
	if(fuel != NULL)
		delete fuel;
	if(oxidizer != NULL)
		delete oxidizer;
}

Reactor::Reactor(const Reactor &copy){
	*this = copy;
}
		
// sets the reactor with premix or cdflame model
void Reactor::setReactorModel(int modelID){
	rModel = modelID;
}

// returns the reactor model PremixFlame or CDflame
int Reactor::getReactorModel(){
	return rModel;
}

// returns the fuel initial conditions object
FlameLab::InitialConditions& Reactor::getFuelInletConditions(){
	if(fuel == NULL) 
		fuel = new InitialConditions();
	return *fuel;
}
//returns the oxidizer initial conditions object
FlameLab::InitialConditions& Reactor::getOxidizerInletConditions(){
	if(oxidizer == NULL) 
		oxidizer = new InitialConditions();
	return *oxidizer;
}
// operator overloaded function for checking the reactor model (premix or cdflame)
//bool Reactor::operator ==(Reactor::ReactorModel rm){
//	return this->rModel == rm;
//}

// fill the reactor with the mixture
void Reactor::setMixture(Sprog::Thermo::Mixture &mix){
	this->mix = &mix;
}
// return the mixture in the reactor
Sprog::Thermo::Mixture& Reactor::getMixture() const{
	return *mix;
}

// set the current mechsnism
void Reactor::setMechanism(const Sprog::Mechanism &mech){
	this->mech = &mech;
}
// returns the current mechanism
const Sprog::Mechanism *const Reactor::getMechanism() const{
	return mech;
}


// set the number of species in the mechanism
void Reactor::setSpeciesCount(int n){
	this->nSpecies = n;
}
// get the number of species
int Reactor::getSpeciesCount(){
	return this->nSpecies;
}

// sets the guess values for species mass/mole fraction
void Reactor::setFraction(const std::string &name, FlameLab::real fraction){
	guessSpecies.insert(pair<std::string, FlameLab::real>(name,fraction));
}
//returns the initial guess conditions
map<string,FlameLab::real> Reactor::getInitialGuess() const{
	return guessSpecies;
}

void Reactor::setUserTemperature(FlameLab::real x, FlameLab::real T){	
	fitProperty fP;
	fP.T = T;
	//set 0 to the rest
	fP.a = 0.0;
	fP.b = 0.0;
	fP.c = 0.0;
	userTemperature.insert(pair<real,fitProperty>(x,fP));
}
/* Natural cubic spline fitting function
Over n intervals the routine fits n equations subject to the
boundary conditions of n+1 data points. The spline variable is zero
for the frist segment and the last segment for natural spline

         z1       z2       z3       z4       z5  
|--------|--------|--------|--------|--------|-------|
0  h0    1   h1   2   h2   3   h3   4   h4   5  h5   6
 
*/

void Reactor::naturalCubicSplineFit(){

	setTempFit(SPLINE);
	vector<FlameLab::real> h,x,T;	
	vector<real> a,b,c,d;
	fitProperty s;

	int nData = userTemperature.size();
	map<real,fitProperty>::iterator p;

	for(p=userTemperature.begin(); p!= userTemperature.end(); p++){
		x.push_back(p->first);
		s = p->second;
		T.push_back(s.T);
	}
	// get the distance between adjuscent points
	for(int i=0;i<nData-1;i++)
		h.push_back(x[i+1]-x[i]);


	//Formulate the matrix
	a.push_back(0);
	b.push_back(0);
	c.push_back(0);
	for(int i=1; i<nData-1; i++){
		//diagonal elements
		b.push_back((h[i]+h[i-1])*2);
		if(i == nData-2)
			c.push_back(0);
		else
			c.push_back(h[i]);
		if(i==1)
			a.push_back(0);
		else
			a.push_back(h[i-1]);
	}	



	d.push_back(0);
	for(int i=1; i<nData-1; i++){
		real dVal = ( (T[i+1]-T[i])/h[i] -  (T[i]-T[i-1])/h[i-1] );
		d.push_back(6*dVal);
	}



	// TDMA implementation
	c[1] /= b[1];
	d[1] /= b[1];
	
	double *z = new double[d.size()];

	for(int i=2; i<nData-1; i++){

		FlameLab::real divide = b[i]- c[i-1]*a[i];
		c[i] /= divide;
		d[i] = (d[i]-d[i-1]*a[i])/divide;
	}
	// back substitution
	z[nData-2] = d[nData-2];
	for(int i = nData-3; i >0; i--)
		z[i] = d[i] - z[i+1]*c[i];

	//store the spline coefficients
	splineCoeff.push_back(0);
	for(int i=1;i<nData-1; i++)
		splineCoeff.push_back(z[i]);

	splineCoeff.push_back(0);
	delete []z;
	
	p = userTemperature.begin();
	for(int i=0;i<nData-1;i++){
		
		p->second.a = (splineCoeff[i+1]-splineCoeff[i])/(6*h[i]);
		p->second.b = splineCoeff[i]/2;
		p->second.c = (T[i+1]-T[i])/h[i] - (h[i]*(2*splineCoeff[i]+splineCoeff[i+1])/6.0);
		//p->second.d = T[i];

		p++;
	}
	

}
//linear fit
void Reactor::linearFit(){

	setTempFit(LINEAR);
	vector<FlameLab::real> x,T;
	//int nData = userTemperature.size();
	map<real,fitProperty>::iterator p;
	fitProperty linear;

	for(p=userTemperature.begin(); p!= userTemperature.end(); p++){
		x.push_back(p->first);
		linear = p->second;
		T.push_back(linear.T);
	}
	p = userTemperature.begin();
	for(unsigned int i=0; i< x.size()-1; i++){
		real m =  (T[i]-T[i+1])/(x[i]-x[i+1]);
		p->second.a = m;      //slope
		p->second.b = T[i]-(m*x[i]);	// intersect	
		p++;
	}

}
// calculates and returns the user defined temperature for an interpolated point
FlameLab::real Reactor::getUserTemperature(const FlameLab::real position){

	map<real,fitProperty>::iterator bound, lowerBound;
	fitProperty fp;
	//check for out of bound exception
	bound = userTemperature.end();
	bound--;
	if(bound->first < position)
		throw ErrorHandler("User defined temperature position out of bounds exception\n",304);
	bound = userTemperature.find(position);
	if(bound == userTemperature.end()){
		bound = userTemperature.upper_bound(position);
		real x = bound->first;
		lowerBound = --bound;
		real xLower = lowerBound->first;
		bound = fabs(position-x) > fabs(position-xLower)? lowerBound : ++bound;
	}

	real temp;
	if (getTempFit() == SPLINE){
		fp = bound->second;
		real diff = position- bound->first;
		temp = (fp.a*pow(diff,3) + fp.b*pow(diff,2) + fp.c*diff+ fp.T);	
	}else{		
		bound = userTemperature.upper_bound(position);
		--bound;
		real slope = bound->second.a;
		real intersect = bound->second.b;
		temp = slope*position + intersect;
	}
	return  temp;
	
}
		
void Reactor::setInitialGuessCondition(int n){
	guessCond = n;
}

int Reactor::getInitialGuessCondition() const{
	return guessCond;
}
//set the strain rate 1/s
void Reactor::setStrainRate(FlameLab::real srate){
	strainRate = srate;
}
//return the strain rate in 1/s
FlameLab::real Reactor::getStrainRate() const{
	return strainRate;
}
//set the temperature fit
void Reactor::setTempFit(int n){
	tempFit = n;
}
// return the temperature fit
int Reactor::getTempFit() const{
	return tempFit;
}

//set the energy solution option
void Reactor::setSolveEnergy(bool option){
	solveEnergy = option;
}

//return the energy solution option
bool Reactor::getSolveEnergy(){
	return solveEnergy;
}


		
