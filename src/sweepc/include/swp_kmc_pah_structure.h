/*!
  * \author     Zakwan Zainuddin (zz260)
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

#include "swp_kmc_mech.h"
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
            
            //! Remove all data and free any memory
            void clear();

            void setParent(Sweep::AggModels::PAH* parent);

            //! Overloaded Operators
            PAHStructure &operator=(const PAHStructure &rhs);
            bool operator==(PAHStructure &rhs) const;
            bool operator!=(PAHStructure &rhs) const;

            //! Stores coordinates of all Carbon atoms (not according to order)
            std::set<cpair> m_cpositions;
            //! Initialise pah with pyrene (currently) or benzene
            void initialise(StartingStructure ss);
            PAHStructure* Clone() ;
            //! return number of carbon and hydrogen for particular PAH
            int numofC() const;
            int numofH() const;
            //! return num of 6-membered rings
            int numofRings() const;
            //! return num of 5-membered rings
            int numofRings5() const;
            //! return num of edge carbon
            int numofEdgeC() const;
            //! return num of site
            int numofSite() const;
            //! set number of carbon and hydrogen for particular PAH
            void setnumofC(int val);
            void setnumofH(int val);

            //! set number of rings for particular PAH
            void setnumofRings(int val); // 6-membered
            void setnumofRings5(int val); // 5-membered

            //! check PAH have bridge or not
            bool havebridgeC();

            Sweep::AggModels::PAH* m_parent; // pointer to parent PAH
            
            void saveDOTperLoop(int PAH_ID, int i);

            //! serialization (incomplete)
            void Serialize(std::ostream &out) const;
            void Deserialize(std::istream &in);

			std::list<Site> GetSiteList() const;

			std::map<kmcSiteType, svector> GetSiteMap() const;

        private:
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
            int m_rings; // 6-membered rings
            int m_rings5; // 5-membered rings
        private:
            //! write m_cpositions
            void WriteCposition(std::ostream &out) const;
            //! read in m_cpositions
            void ReadCposition(std::istream &in, const int size);
            //! Copy Constructor
            PAHStructure(const PAHStructure& copy);
            //! Set storing carbon objects
            Ccontainer m_carbonList;

            //Cpointer NULLC;
        };
    }
}

#endif
