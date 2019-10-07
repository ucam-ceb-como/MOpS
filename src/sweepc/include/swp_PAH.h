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
#include "swp_particle_model.h"

#ifndef SWEEP_PAH_H
#define SWEEP_PAH_H
	
namespace Sweep
{
//forward declaration 
namespace KMC_ARS{
class PAHStructure;}

namespace AggModels
{
//! PAH primary particle class
class PAH 
{
	friend class PAHPrimary;
public:
    //! Build a new PAH with time
    PAH(double time, Sweep::ParticleModel::PostProcessStartingStr str);

    PAH(const PAH &copy); // Copy constructor.
    ~PAH();
    PAH &operator=(const PAH &rhs);
    PAH *const Clone() const;
    //returen ID of PAH
    int ID() const;
    void saveDOTperLoop(int ID, int i) const;
    Sweep::KMC_ARS::PAHStructure* Structure();

    //! this default conscustor is only used for serialization.
    PAH();
    void Serialize(std::ostream &out) const;
    void Deserialize(std::istream &in);

	Sweep::KMC_ARS::PAHStructure* GetPAHStruct();
private:

    //! Simulated time at point when molecule created
    double time_created;

    //! Simulated time when molecule size was last updated
    double lastupdated;

    //! Index of the PAH, used to track specified PAH
    int PAH_ID;

    //! PAH structure provide by KMC PAH growth part
    Sweep::KMC_ARS::PAHStructure* m_pahstruct;

};
} //namespace AggModels
} //namespace Sweep

#endif

