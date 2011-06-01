 /*\author Dongping Chen
 *  Copyright (C) 2011 Dongping Chen.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Defines a primary including detailed PAH information
  Author(s):      Dongping Chen (dc516)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2009 Markus Sander.

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
#include "swp_PAH.h"

using namespace Sweep;
using namespace Sweep::AggModels;


using namespace std;

PAH::PAH(){
}

PAH::PAH(const PAH &copy){
	*this=copy;
}

PAH::~PAH() {
	//delete m_pahstruct;
}

PAH &PAH::operator=(const PAH &rhs){
	if (this != &rhs){
		m_numcarbon=rhs.m_numcarbon;
		time_created=rhs.time_created;
	    lastupdated=rhs.lastupdated;
		PAH_ID=rhs.PAH_ID;//should add a value to diff those two PAHs after cloning
		//PAHStructure* m_pahstruct=new PAHStructure();
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