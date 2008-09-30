/*
  Author(s):      Vinod Janardhanan (vj231)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Vinod M Janardhanan.

  File purpose:
    This file contains the definition of a structure for a chemical species.  File
    also contains typdefs and structure definitions related to Species objects.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
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


#include "gpc_transport_data.h"
#include "gpc_species.h"
#include "gpc_mech.h"
#include <string>


using namespace Sprog::Transport;
using namespace std;

TransportData::TransportData(){
	molIndex = 0;
	LJwellDepth = 0.;
	LJcollisionDia = 0.;
	dipol = 0.;
	polarity = 0.;
	rotRelaxNum = 0.;
}

TransportData::TransportData(int mi, 
							 double LJDepth, 
							 double LJcd, 
							 double dpol, 
							 double polar, 
							 double rotRelax){

	
	TransportData::molIndex = mi;	
	TransportData::LJwellDepth = LJDepth;
	TransportData::LJcollisionDia = LJcd;
	TransportData::dipol = dpol;
	TransportData::polarity = polar;
	TransportData::rotRelaxNum = rotRelax;
}

//TransportData::TransportData(const Sprog::Transport::TransportData &td){
//	*this = td;
//}

TransportData TransportData::operator =(Sprog::Transport::TransportData td){
	molIndex = td.molIndex;
	LJwellDepth = td.LJwellDepth;
	LJcollisionDia = td.LJcollisionDia;
	dipol = td.dipol;
	polarity = td.polarity;
	rotRelaxNum = td.rotRelaxNum;

	return *this;

}


void TransportData::setMolIndex(int mi){
	this->molIndex = mi;
}

int TransportData::getMolIndex() const{
	return this->molIndex;
}

void TransportData::setWellDepth(double wd){
	this->LJwellDepth = wd;
}

double TransportData::getWellDepth() const{
	return this->LJwellDepth;
}

void TransportData::setCollisionDia(double cd){
	this->LJcollisionDia = cd;
}

double TransportData::getCollisionDia() const{
	return this->LJcollisionDia;
}


void TransportData::setDipole(double dipol){
	this->dipol = dipol;
}

double TransportData::getDipole() const{
	return this->dipol;
}



void TransportData::setPolarity(double p){
	this->polarity = p;
}

double TransportData::getPolarity() const{
	return this->polarity;
}

void TransportData::setRotRelaxNum(double rrn){
	this->rotRelaxNum = rrn;
}

double TransportData::getRotRelaxNum() const{
	return this->rotRelaxNum;
}

void TransportData::validateTransport(std::map<string,vector<string> > &trMap,
									  Sprog::Mechanism &mech)
{
	vector<Species*> vSpecies = mech.Species();
	map<std::string,vector<string>>::iterator p;
	
	int vsize = vSpecies.size();
	
	for(int i=0; i< vsize; i++){
		p = trMap.find(vSpecies[i]->Name());
		if( p== trMap.end()){
			cout << " Species :" << vSpecies[i]->Name() << "not found in transport data!!" << endl;
			exit(1);
		
		}

	}

}




	
	






