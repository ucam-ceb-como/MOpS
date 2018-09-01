 /*
  Author(s):      Dongping Chen(dc516)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2011 Dongping Chen

  File purpose:


  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

#include "swp_kmc_pah_structure.h"
#include "swp_kmc_typedef.h"
#include "swp_PAH.h"

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace Sweep::KMC_ARS;

using namespace std;


PAH::PAH():time_created(0),lastupdated(0),PAH_ID(0),m_pahstruct(new PAHStructure()){
}

PAH::PAH(double time, Sweep::ParticleModel::PostProcessStartingStr str):
time_created(time),
lastupdated(time),
PAH_ID(0),
m_pahstruct(new PAHStructure())
{
    switch (str){
    case ParticleModel::A1:
        m_pahstruct->initialise(BENZENE_C);
        break;
    case ParticleModel::A2:
        m_pahstruct->initialise(NAPHTHALENE_C);
        break;
    case ParticleModel::A4:
        m_pahstruct->initialise(PYRENE_C);
        break;
    case ParticleModel::A5:
        m_pahstruct->initialise(BENZOPYRENE_C);
        break;
    }
}

PAH::PAH(const PAH &copy){
	*this=copy;
}

PAH::~PAH() {
	delete m_pahstruct;
}

PAH &PAH::operator=(const PAH &rhs){
	if (this != &rhs){
		time_created=rhs.time_created;
	    lastupdated=rhs.lastupdated;
		PAH_ID=rhs.PAH_ID;
		m_pahstruct=rhs.m_pahstruct->Clone();
		m_pahstruct->setParent(this);
	}
	return *this;
}

PAH *const PAH::Clone() const{
	return new PAH(*this);
}

int PAH::ID() const {
	return PAH_ID;
}

void PAH::saveDOTperLoop(int ID, int i) const{
	m_pahstruct->saveDOTperLoop(ID,i);
}

Sweep::KMC_ARS::PAHStructure* PAH::Structure()
{
   return m_pahstruct;
}


void PAH::Serialize(std::ostream &out) const
{
    double val=0.0;

    m_pahstruct->Serialize(out);

    val=time_created;
    out.write((char*)&(val), sizeof(val));
    val=lastupdated;
    out.write((char*)&(val), sizeof(val));
    val=PAH_ID;
    out.write((char*)&(val), sizeof(val));
}
void PAH::Deserialize(std::istream &in)
{
    if (in.good()) {
    // Read the output version.  Currently there is only one
    // output version, so we don't do anything with this variable.
    // Still needs to be read though.

	double val = 0.0;

    m_pahstruct=new PAHStructure();
    m_pahstruct->Deserialize(in);

    in.read(reinterpret_cast<char*>(&val), sizeof(val));
    time_created = (double)val;

    in.read(reinterpret_cast<char*>(&val), sizeof(val));
    lastupdated = (double)val;

    in.read(reinterpret_cast<char*>(&val), sizeof(val));
    PAH_ID = (int)val;
    }
}

Sweep::KMC_ARS::PAHStructure* PAH::GetPAHStruct(){
	return m_pahstruct;
}
