/*!
  * \author     Zakwan Zainuddin (zz260) && Gustavo Leon (gl413)
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
			//! Initialise pah from file (debugging purposes)
            void initialise_fromfile();
            PAHStructure* Clone() ;
            //! return number of carbon and hydrogen for particular PAH
            int numofC() const;
            int numofH() const;
			int numofCH3() const;
            //! return num of 6-membered rings
            int numofRings() const;
			//! return num of lone 5-membered rings
			int numofLoneRings5() const;
			//! return num of embedded 5-membered rings
			int numofEmbeddedRings7() const;
			//! return num of lone 7-membered rings
			int numofLoneRings7() const;
			//! return num of embedded 7-membered rings
			int numofEmbeddedRings5() const;
            //! return num of edge carbon
            int numofEdgeC() const;
            //! return num of site
            int numofSite() const;
			//! return num of bridges
            int numofBridges() const;
            //! set number of carbon and hydrogen for particular PAH
            void setnumofC(int val);
            void setnumofH(int val);
			void setnumofCH3(int val);
			

            //! set number of rings for particular PAH
			void setnumofRings(int val); // 6-membered
			void setnumofLoneRings5(int val); // 5-membered
			void setnumofEmbeddedRings5(int val); // 5-membered
			void setnumofLoneRings7(int val); // 7-membered
			void setnumofEmbeddedRings7(int val); // 7-membered

            //! check PAH have bridge or not
            bool havebridgeC();
			
			//! check if the PAH has embedded R5s, R7s or partially embedded R7s that may give curvature or if PAH is large enough.
			bool hasCurvedSubunits();

            Sweep::AggModels::PAH* m_parent; // pointer to parent PAH
            
            void saveDOTperLoop(int PAH_ID, int i);

            //! Save PAH structure as XYZ files. 
            void saveXYZ(const std::string &filename, bool optimise);

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
			//! Stores number of methyl groups.
			int m_methyl_counts;
            //! Stores number of rings
			int m_rings; // 6-membered rings
			int m_rings5_Lone; // 5-membered rings
			int m_rings5_Embedded; // 5-membered rings
			int m_rings7_Lone; // 7-membered rings
			int m_rings7_Embedded; // 7-membered rings
			std::list<cpair> m_InternalCarbons; //List of internal Carbon atoms coordinates
			std::list<cpair> m_R5loc; //List of R5s center coordinates
			std::list<cpair> m_R7loc; //List of R7s center coordinates
            std::vector<std::tuple<Spointer,Spointer,int>> m_R5walker_sites; //Vector of R5 random walkers.
			bool m_optimised = false; //Flag to know if the PAH has been optimised or not.
			
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
