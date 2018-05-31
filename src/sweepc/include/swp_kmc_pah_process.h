/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_pah_process.h
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

#ifndef SWP_KMC_STRUCTURE_PROCESSES_H
#define SWP_KMC_STRUCTURE_PROCESSES_H

//#include "swp_kmc_mech.h"
#include "swp_kmc_structure_comp.h"
#include "swp_kmc_typedef.h"
//#include "swp_kmc_jump_process.h"
#include "swp_kmc_pah_structure.h"
#include <list>
#include <vector>
#include <map>
#include <set>

namespace Sweep{
namespace KMC_ARS{

typedef std::vector<Spointer> svector;
class JumpProcess;

class PAHProcess {
public:
    // Constructors and destructors
    //! Default Constructor
    PAHProcess();
    //! Overloaded constructor; specifying PAH structure to perform processes on
    PAHProcess(PAHStructure& pah);
    //! Copy Constructor
    PAHProcess(const PAHProcess &pahp);
    //! Default Destructor
    virtual ~PAHProcess();

    PAHStructure* returnPAH();
    bool havebridgeC();

    // Check to validate if coordinates of C matches bond angles
    bool checkCoordinates() const;
    // Check to see if all sites are connected to each other
    bool checkSiteContinuity() const;
    // Check to see if site neighbours has a valid combined site type
    bool checkCombinedSiteType(Spointer& stt);

    //! Sets target PAH structure to perform processes on
    void setPAH(PAHStructure& pah);
    //! Returns a copy of PAH structure
    PAHStructure* clonePAH() const;
    //PAHStructure* clonePAH_new() const; // --> new method

    size_t SiteListSize() const;
    size_t CarbonListSize() const;
    std::list<Site>& SiteList() const;


    // Structure change processes
    //! Initialisation of structure given a starting structure
    virtual PAHStructure& initialise(StartingStructure ss);
    //! Initialisation of structure given a starting structure (new method)
    virtual PAHStructure& initialise_new(StartingStructure ss);
    //! Initialisation of structure given a string of site types (separated by ',')
    virtual PAHStructure& initialise(std::string siteList_str, int R6_num, int R5_num);
    //! Create Structure from vector of site types and number of rings
    void createPAH(std::vector<kmcSiteType>& vec, int R6, int R5);
    //! Structure processes: returns success or failure
    bool performProcess(const JumpProcess& jp, rng_type &rng, int PAH_ID);

    // Read Processes
    //! Get Counts
    intpair getCHCount() const;
    //! Get Site Counts
    unsigned int getSiteCount(const kmcSiteType& st) const;
    //! Get Ring Counts
    intpair getRingsCount() const;
    //! Get number of bridges
    int numberOfBridges() const;
    //! Print structure in console
    void printStruct() const;
    //! Print Structure in console, with arrow pointing at current C
    void printStruct(Cpointer c) const;
    //! Print sites in console
    void printSites() const;
    //! Print sites & site members in console, with arrow pointing at site stt
    void printSitesMemb(Spointer& stt) const;
    void printSitesMemb() const;
    //! Print sites in console, with an arrow pointing at site stt
    void printSites(Spointer& stt) const;
    //! Store structure results in external file, returns success/failure
    bool saveDOT(const std::string &filename) const;
    bool saveDOT(const std::string &filename, const std::string &title) const;
    //! obtains a vector of the PAH site list
    std::vector<kmcSiteType> SiteVector() const;
    //! obtains a string containing the PAH site list
    std::string SiteString(char delimiter) const;
    
    // Update Processes

    //! Jump processes: defined in swp_kmc_processes_list.cpp
    //! Growth Processes
    // name structure: proc_processtype_sitetype.
    // refer StructureProc enum in swp_kmc_typedef.h
    void proc_G6R_FE(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_G6R_AC(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_L6_BY6(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_PH_benz(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);
    void proc_D6R_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_O6R_FE3_O2(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_O6R_FE3_OH(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_O6R_FE_HACA_O2(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_O6R_FE_HACA_OH(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_G5R_ZZ(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_D5R_R5(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_C6R_AC_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);
    void proc_C5R_RFE(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_C5R_RAC(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_M5R_RZZ(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_C6R_BY5_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);
    void proc_C6R_BY5_FE3violi(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);

    void proc_L5R_BY5(Spointer& stt, Cpointer C_1, Cpointer C_2);
    void proc_M6R_BY5_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);

    void proc_O6R_FE2(Spointer& stt, Cpointer C_1, Cpointer C_2);
    //void proc_M5R_eZZ(Spointer& stt, Cpointer C_1, Cpointer C_2);
    
    // true: saves rates only, returns all site count as 1
    // false: doesn't save rates, returns actual site counts
    bool m_rates_save;

private:
    // Read Process
    //! Get other member of the site a particular C atom is a member of
    Cpointer getPair(
        const Cpointer Carb, 
        // Is position of other member preceeding or after the C atom?
        // 0:prev, 1:next
        bool after 
        ) const;
    //! Move current and previous Carbon iterator to the next carbon
    Cpointer moveCPointer(Cpointer &previous, Cpointer &current) const;
    //! Check if process is allowed
    bool allowed(const Spointer& st, StructureProc proc) const;
    //! Choose a random site of site type st
    Spointer chooseRandomSite(kmcSiteType st, rng_type &rng);
    //! Choose a random site of any site types in vtype
    Spointer chooseRandomSite(std::vector<kmcSiteType> vtype, rng_type &rng);
    //! Jump to a position coordinate given starting position and angle towards new position
    cpair jumpToPos(const cpair& starting, const angletype& direction) const;
    //! Search a particular site (si) from svector associated with stype and erases it from sitemap
    void delSiteFromMap(const kmcSiteType& stype, const Spointer& si);
    //! Overload: search and erase from svectors associated with all site types in vector v
    void delSiteFromMap(const std::vector<kmcSiteType>& v, const Spointer& si);
    //! Check for steric hindrance for a site growth process
    bool checkHindrance(const Spointer& st) const;
    //! Check steric hindrance for phenyl addition reactions
    bool checkHindrancePhenyl(const Cpointer C_1) const;
    //! Returns site iterator x steps after i
    Spointer moveIt(Spointer i, int x);
    //! Finds C atom with specific coordinates
    Cpointer findC(cpair coordinates);

    // Write Process
    //! Creates a lone carbon atom
    Cpointer addC();
    //! Creates a new carbon atom attached next to C_1.
    Cpointer addC(Cpointer C_1, angletype angle1, angletype angle2, bool bulk=false);
    //! Creates a carbon atom bridging next to C_1. 
    Cpointer bridgeC(Cpointer C_1);
    /*//! Creates a bulk carbon atom connected to C_1
    Cpointer addBC(Cpointer C_1);
    //! Connects an edge carbon atom to another edge or bulk carbon atom
    void connectToC(Cpointer C_1, Cpointer C_2, bool);*/
    //! Connects a carbon atom to a carbon (to close loop)
    void connectToC(Cpointer C_1, Cpointer C_2);
    /*! Removes a carbon atom from perimeter taking into account if
      ! the atom becomes bulk carbon within the PAH
    */ 
    void removeC(Cpointer C_1, bool bulk);
    //! Adds a site with members C_1 & C_2 before site sIt
    Spointer addSite(kmcSiteType stype, 
        Cpointer C_1, Cpointer C_2, Spointer& sIt);
    //! Adds a site with members C_1 & C_2 at end of SiteList
    Spointer addSite(kmcSiteType stype, 
        Cpointer C_1, Cpointer C_2);
    //! Removes a site
    //void removeSite(Spointer st);
    //! Changes site type into combined site with R5 (e.g. FE -> RFE)
    void addR5toSite(Spointer& st, Cpointer Carb1, Cpointer Carb2);
    //! Changes site type into combined site without R5 (e.g. RFE -> FE)
    void remR5fromSite(Spointer& st, Cpointer Carb1, Cpointer Carb2);
    //! Changes site type into another site type
    void convSiteType(Spointer& st, Cpointer Carb1, Cpointer Carb2, kmcSiteType t);
    //! Remove site
    void removeSite(Spointer& stt);
    //! Sets the number of counts of C and H
    void setCount(int CCount, int HCount);
    //! Add counts
    void addCount(int C_in, int H_in);
    //! For createPAH function: drawing type 0 sites
    Cpointer drawType0Site(Cpointer Cnow, int bulkC);
    //! For createPAH function: drawing type 1 sites
    Cpointer drawType1Site(Cpointer Cnow, int bulkC, kmcSiteType prevType);
    //! For createPAH function: drawing type 2 sites
    Cpointer drawType2Site(Cpointer Cnow, int bulkC);

    //! Update Sites and its members in structure
    //! All principal sites
    void updateSites();
    //! Updates particular site
    void updateSites(Spointer& st, // site to be updated
        Cpointer Carb1, Cpointer Carb2, // new C members
        int bulkCchange); // addition to number of bulk C in site
    //! Combined sites for all sites
    void updateCombinedSites();
    //! Combined site for a particular site
    void updateCombinedSites(Spointer& st);
    //! Sets third species bonded to C to a species sp if it is a reactive surface carbon
    void updateA(Cpointer C, char sp);
    //! Overload function, updateA for all C from C_1 to C_2 inclusive
    void updateA(Cpointer C_1, Cpointer C_2, char spc);
    
    // PAH data structure to perform processes on
   PAHStructure* m_pah;

    //Cpointer NULLC;
};

}
}


#endif
