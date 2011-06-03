/*!
  * \Author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_pah_structure.h
  *
  * \brief      Defines the data structure which holds information of PAH molecules
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Holds the data structure which contains information on the molecular structure
    of PAH accounted by the kMC model. The structure is a directed graph connecting
    C atoms on the PAH perimeter, with the direction going clockwise along the edges.

  Licence:
    This file is part of "sweep".

    Sweep is free software; you can redistribute it and/or
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

#ifndef SWP_KMC_PAH_STRUCTURE_H
#define SWP_KMC_PAH_STRUCTURE_H

#include "swp_kmc_gasph.h"
#include "swp_kmc_structure_comp.h"
#include "swp_kmc_typedef.h"
#include "swp_kmc_jump_process.h" 
#include "swp_PAH.h"
#include <list>
#include <vector>
#include <map>
#include <set>

namespace Sweep{
    namespace KMC_ARS{
        class JumpProcess;
        typedef std::vector<Spointer> svector;
        class PAHStructure{
        public:
            friend class PAHProcess;
            // Constructors and destructors
            //! Default Constructor
            PAHStructure();
            //! Default Destructor
            ~PAHStructure();
			
			void setParent(Sweep::AggModels::PAH* parent);

			//! Overloaded Operators
			//PAHStructure &operator=(const PAHStructure &rhs);
			bool operator==(PAHStructure &rhs) const;
			bool operator!=(PAHStructure &rhs) const;

			//void CopyLoop(const PAHStructure* source);
            //! Stores coordinates of all Carbon atoms (not according to order)
            std::set<cpair> m_cpositions;
							//add by dongping 12.04. however, we should keep the PAHStructure
							//object to be small in case memory problem, so we should not put these member
							//functions in this class, however, it seems that member functions
							//in specified class do not take any extra memory
			void initialise(StartingStructure ss);
			PAHStructure* Clone() ;
							//void clearStructure();
							//Cpointer addC(Cpointer C_1, angletype angle1, angletype angle2);
							//Cpointer addC();
							//void connectToC(Cpointer C_1, Cpointer C_2);
							//void addCount(int C_in, int H_in);
							//void setCount(int CCount, int HCount);
							//void updateA(Cpointer C, char sp);
							//void updateA(Cpointer C_1, Cpointer C_2, char spc);
							//cpair jumpToPos(const cpair& starting, const angletype& direction) const;
							//void updateSites();
							//void updateSites(Spointer& st, // site to be updated
				   //                            Cpointer Carb1, Cpointer Carb2, // new C members
				   //                            int bulkCchange);
							//void updateCombinedSites();
							//void updateCombinedSites(Spointer& st);
							//void delSiteFromMap(const std::vector<kmcSiteType>& v,const Spointer& si);
							//void delSiteFromMap(const kmcSiteType& stype, const Spointer& si);
				   //         Spointer moveIt(Spointer i, int x);
							//Spointer addSite(kmcSiteType stype, Cpointer C_1, Cpointer C_2, Spointer& b4site);
				   //         void printSites();
							//void printSites(Spointer& stt);
							//Cpointer moveCPointer(Cpointer &previous, Cpointer &current) const;
			int numofC();
			bool havebridgeC();
							////Carbon* Carbon::Clone() const;//added by dongping14.04
			Sweep::AggModels::PAH* m_parent; // pointer to parent PAH
			
        protected:
            //! First and last Carbon atom in list
            Cpointer m_cfirst;
            Cpointer m_clast;
            //! Stores all principal PAH sites in order from m_cfirst-m_clast.
            std::list<Site> m_siteList;
			//! Stores iterators to the PAH sites according to their site type
            std::map<kmcSiteType, svector> m_siteMap;
			//! Stores total counts of carbon and hydrogen
            intpair m_counts;
            //! Stores number of rings
            int m_rings;
		private:
			//! Copy Constructor
			PAHStructure(const PAHStructure& copy);
        };
    }
}

#endif
