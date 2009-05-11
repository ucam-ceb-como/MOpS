/*
  Author(s):      Markus Sander (ms785)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Markus Sander.

  File purpose:
    Implementation of the PAHPrimary class declared in the
    swp_PAH_primary.h header file.

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

#include "swp_primary.h"
#include "swp_PAH_primary.h"
#include "swp_submodel_type.h"
#include "swp_aggmodel_type.h"
#include "swp_model_factory.h"
#include "swp_PAH_cache.h"
#include <math.h>
#include "rng.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace Sweep::SubModels;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.



// Initialising constructor.
PAHPrimary::PAHPrimary(real time, const Sweep::ParticleModel &model)
: Primary(time, model)
{
	m_numcarbon=0;
	m_comp[0]=1;
	PAH *new_PAH=new PAH;
	new_PAH->ID=0;
	new_PAH->m_numcarbon=0;
	new_PAH->time_created=time;
	m_PAH.insert(m_PAH.begin(),*new_PAH);
	this->UpdateCache();
}

// Copy constructor.
PAHPrimary::PAHPrimary(const PAHPrimary &copy) 
{
    *this = copy;
}


PAHPrimary &PAHPrimary::operator=(const Primary &rhs)
{
    if (this != &rhs) {
    }
    return *this;
}



// Returns the collision diameter.
real PAHPrimary::PAHCollDiameter(void) const {return m_dcol;}


// Stream-reading constructor.
PAHPrimary::PAHPrimary(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
PAHPrimary::~PAHPrimary()
{
    releaseMem();
}



PAHPrimary &PAHPrimary::Coagulate(const Primary &rhs)
{
    // Use spherical model to coagulate.
    //Primary::Coagulate(rhs);

    // Now deal with the PAHs.
    vector<PAH>::const_iterator j;
	const AggModels::PAHPrimary *pah = NULL;
	pah = dynamic_cast<const AggModels::PAHPrimary*>(&rhs);
	for (j=pah->m_PAH.begin(); j!=pah->m_PAH.end(); ++j) {
		this->m_PAH.insert(m_PAH.end(),PAH(*j));
	}
	UpdateCache();
    return *this;
}


void PAHPrimary::UpdateTime(double t)
{
	unsigned int j=0;
	for (vector<PAH>::iterator i=m_PAH.begin(); i!=m_PAH.end(); ++i) {
		for (j=0;j<i->time.size();j++)
		{	
			i->m_numcarbon=i->n_carbon_t.at(j);
			//double timeatpos=i->time.at(j);
			if (i->time.at(j)>=t-i->time_created)
			{
				break;
			}
		}
	}
}

void PAHPrimary::UpdateCache(void)
{	m_numcarbon=0;
	m_PAHmass=0;
	m_PAHCollDiameter=0;
	m_numPAH=0;
	 for (vector<PAH>::iterator i=m_PAH.begin(); i!=m_PAH.end(); ++i) {
        m_numcarbon += i->m_numcarbon;
		m_PAHCollDiameter=max(m_PAHCollDiameter,2.4162*sqrt((i->m_numcarbon)*2.0/3.));    // in Angstrom
		m_numPAH++;
    }
	 m_PAHmass=m_numcarbon*1.9945e-26;        //convert to kg
	m_PAHCollDiameter*=1e-10;//convert from Angstrom to m

	 double V=m_PAHmass/1.6e3;
	 double cdiamsphere=pow((3*V/(4*PI)),(1.0/3.0));
	 double alpha=1.0/(pow(m_numPAH,0.1));
	// double cdiam=alpha*m_PAHCollDiameter+(1-alpha)*cdiamsphere;
	double cdiam=max(cdiamsphere,m_PAHCollDiameter);
	 SetCollDiameter(cdiam);	
	 SetMass(m_PAHmass);
	

	

}

// Creates an aggregation data cache for this primary type.
AggModels::PAHCache *const PAHPrimary::CreateAggCache(ParticleCache &pcache) const
{
    PAHCache *cache = 
        static_cast<PAHCache*>(ModelFactory::CreateAggCache(AggModels::PAH_ID, pcache));
    if (cache != NULL) *cache = *this;
    return cache;
}





// READ/WRITE/COPY.

// Returns a copy of the model data.
PAHPrimary *const PAHPrimary::Clone(void) const
{
    return new PAHPrimary(*this);
}

// Returns this object's instance.  This may seem rather circular, but
// it has an important purpose for getting the correct object type reference
// from a base class reference/pointer.
PAHPrimary &PAHPrimary::Instance() {return *this;}
const PAHPrimary &PAHPrimary::Instance() const {return *this;}


