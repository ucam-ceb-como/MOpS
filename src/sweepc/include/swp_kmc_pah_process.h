/*!
  * \author     Zakwan Zainuddin (zz260) && Gustavo Leon (gl413)
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

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

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
	// Check if site stt is a valid site type.
	bool checkSiteValid(int type);
	// Check if site stt is a valid site type.
	bool checkSiteValid(Spointer& stt);
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
	//! Initialisation of structure given an existing PAH (cloning)
	virtual PAHStructure& initialise(std::string siteList_str, int R6_num, int R5_num_Lone, int R5_num_Embedded, int R7_num_Lone, int R7_num_Embedded, Ccontainer edgeCarbons, std::list<cpair> internalCarbons, std::list<cpair> R5_locs, std::list<cpair> R7_locs);
    //! Create Structure from given an existing PAH (cloning)
	void createPAH(std::vector<kmcSiteType>& vec, std::vector<int>& carb_vec, int R6, int R5_Lone, int R5_Embedded, int R7_Lone, int R7_Embedded, Ccontainer edCarbons, std::list<cpair> inCarbs, std::list<cpair> R5loc, std::list<cpair> R7loc, cpair first_carbon_coords);
    //! Initialisation of structure given a string of site types (separated by ',')
	virtual PAHStructure& initialise(std::string siteList_str, int R6_num, int R5_num_Lone, int R5_num_Embedded, int R7_num_Lone, int R7_num_Embedded, int num_C, int num_H, std::list<cpair> internalCarbons);
    //! Create Structure from vector of site types and number of rings
	void createPAH(std::vector<kmcSiteType>& vec, int R6, int R5_Lone, int R5_Embedded, int R7_Lone, int R7_Embedded, int num_C, int num_H);
    //! Structure processes: returns success or failure
    bool performProcess(const JumpProcess& jp, rng_type &rng, int PAH_ID);

    // Read Processes
    //! Get Counts
    intpair getCHCount() const;
    //! Get Site Counts
    unsigned int getSiteCount(const kmcSiteType& st) const;
    //! Get Ring Counts
    std::tuple <int, int, int> getRingsCount() const;
	//! Get R5 embedded count
    int getR5EmbeddedCount() const;
	//! Get R7 embedded count
    int getR7EmbeddedCount() const;
    //! Get number of bridges
    int numberOfBridges() const;
	//! Get number of methyl moeities
	int numberOfMethyl() const;
	//! Is probably curved (has embedded or partially embedded R5s, R7s, more than 8 rings)
	bool HasCurvedMoeities() const;
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
	//! Create a copy of the sites type
	std::list<std::string> copySites() const;
	//! Create a copy of the sites type
	std::list<std::string> copySites(Spointer& stt) const;
	//! Print copied sites before JP in console
	void printBeforeSites(std::list<std::string>& before_site_list) const;
    //! Store structure results in external file, returns success/failure
    bool saveDOT(const std::string &filename) const;
    bool saveDOT(const std::string &filename, const std::string &title) const;
	void saveXYZ(const std::string &filename, bool optimise=false);
	void save_trajectory_xyz(const double timer, const std::string &filename, bool optimise=false);
	bool saveDOT3D(const std::string &filename) const;
	bool saveDOT3D(const std::string &filename, const std::string &title) const;
    //! obtains a vector of the PAH site list
    std::vector<kmcSiteType> SiteVector() const;
	//! obtains a vector of the carbons per site in PAH site list
	std::vector<int> SiteIntVector() const;
    //! obtains a string containing the PAH site list
    std::string SiteString(char delimiter) const;
	//! Passes a PAH from MOpS to OpenBabel. Returns a mol object.
	OpenBabel::OBMol passPAH(bool detectBonds=true);
	//! Connects the atoms in a PAH using OpenBabel routines. Equivalent to OpenBabel::OBMol::ConnectTheDots();
	void connectPAH(OpenBabel::OBMol my_mol);
	//! Passes a PAH from OpenBabel to MOpS.
	void passbackPAH(OpenBabel::OBMol mol);
	//! Runs optimisation of a PAHusing Openbabel
	OpenBabel::OBMol optimisePAH(OpenBabel::OBMol mol, int nsteps=4000, std::string forcefield="mmff94") ;
	//! Includes curvature in a PAH after a pentagon is integrated in the structure.
	//OpenBabel::OBMol includeCurvature(OpenBabel::OBMol mol, cpair CR5_1, cpair CR5_2) ;
    
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
    
	void proc_D6R_FE_AC(Spointer& stt, Cpointer C_1, Cpointer C_2);
	void proc_B6R_ACR5(Spointer& stt, Cpointer C_1, Cpointer C_2);                          //!< ID22.
	void proc_M5R_ACR5_ZZ(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);        //!< ID23.
	void proc_G6R_RZZ(Spointer& stt, Cpointer C_1, Cpointer C_2);                           //!< ID24.
	void proc_G6R_RFER(Spointer& stt, Cpointer C_1, Cpointer C_2);                          //!< ID25.
	void proc_G6R_R5(Spointer& stt, Cpointer C_1, Cpointer C_2);                            //!< ID26.
	void proc_L6_RBY5(Spointer& stt, Cpointer C_1, Cpointer C_2);                           //!< ID27.
	void proc_L6_RACR(Spointer& stt, Cpointer C_1, Cpointer C_2);                           //!< ID28.
	void proc_G5R_RFE(Spointer& stt, Cpointer C_1, Cpointer C_2);                           //!< ID29.
	void proc_C6R_RAC_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);        //!< ID30.
	void proc_C6R_RAC_FE3violi(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);   //!< ID31.
	void proc_M6R_RAC_FE3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);        //!< ID32.
	void proc_MR5_R6(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);             //!< ID34.
	void proc_GR7_R5R6AC(Spointer& stt, Cpointer C_1, Cpointer C_2);             			//!< ID35.
	void proc_GR7_FEACR5(Spointer& stt, Cpointer C_1, Cpointer C_2);             			//!< ID36.
	void proc_G6R_R5R6ZZ(Spointer& stt, Cpointer C_1, Cpointer C_2);                        //!< ID37.
	void proc_L7_ACACR5(Spointer& stt, Cpointer C_1, Cpointer C_2);                         //!< ID38.
	void proc_G6R_R5R6FER(Spointer& stt, Cpointer C_1, Cpointer C_2);                       //!< ID39.
	void proc_G6R_R5R6FER5R6(Spointer& stt, Cpointer C_1, Cpointer C_2);                    //!< ID40.
	void proc_L7_FEZZACR5(Spointer& stt, Cpointer C_1, Cpointer C_2);                       //!< ID41
	void proc_C5R_RZZR(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);           //!< ID42
	void proc_C5R_R5R6ZZR(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);        //!< ID43
	void proc_L6_R5R6BY5(Spointer& stt, Cpointer C_1, Cpointer C_2);                        //!< ID44
	void proc_L6_R5R6ACR(Spointer& stt, Cpointer C_1, Cpointer C_2);                        //!< ID45
	void proc_L6_R5R6ACR5R6(Spointer& stt, Cpointer C_1, Cpointer C_2);                     //!< ID46
	void proc_L6_ZZACR5(Spointer& stt, Cpointer C_1, Cpointer C_2);                         //!< ID47
	void proc_L6_R5FEACR5(Spointer& stt, Cpointer C_1, Cpointer C_2);                       //!< ID48
	void proc_L6_FEACR5FE(Spointer& stt, Cpointer C_1, Cpointer C_2);                       //!< ID49
	void proc_L6_R5ACR5R5(Spointer& stt, Cpointer C_1, Cpointer C_2);                       //!< ID50
	void proc_L7_R5ZZACR5(Spointer& stt, Cpointer C_1, Cpointer C_2);                       //!< ID51
	void proc_L6_ACR5R5R6(Spointer& stt, Cpointer C_1, Cpointer C_2);                       //!< ID52
	void proc_L7_ACR5R5R6ZZ(Spointer& stt, Cpointer C_1, Cpointer C_2);                     //!< ID53
	void proc_A_CH3(Spointer& stt, Cpointer C_1, Cpointer C_2, rng_type &rng);              //!< ID54
	void proc_D_CH3(Spointer& stt, Cpointer C_1, Cpointer C_2);                     		//!< ID55

    // true: saves rates only, returns all site count as 1
    // false: doesn't save rates, returns actual site counts
    bool m_rates_save;
	// true, adds intermediate save points before each jump process to debug.
	bool debug_pah=false;

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
	cpair jumpToPos(const cpair& starting, const angletype& direction, const angletype& direction2, const bondlength& distance) const;
	//! Jump to a position coordinate given starting position and a vector. Handles 3D growth.
	cpair jumpToPos(const cpair& starting, const cpair& direction, const bondlength& distance) const;
	//! Gets a vector between two points.
	cpair get_vector(cpair p1, cpair p2)const;
	//! Rescales a vector.
	cpair scale_vector(cpair vec)const;
	//! Changes vector direction.
	cpair invert_vector(cpair vec)const;
	//! Returns the resultant of two vectors.
	cpair add_vector(cpair vec1, cpair vec2)const;
	//! Returns the normal vector to the face of a ring.
	cpair norm_vector(cpair p1, cpair p2, cpair p3)const;
	//! Returns the cross product of two vectors. If v1 is the surface normal vector and v2 goes from C->newC v1xv2 redefines the next growth vector.
	cpair cross_vector(cpair vec1, cpair vec2)const;
    //! Search a particular site (si) from svector associated with stype and erases it from sitemap
    void delSiteFromMap(const kmcSiteType& stype, const Spointer& si);
    //! Overload: search and erase from svectors associated with all site types in vector v
    void delSiteFromMap(const std::vector<kmcSiteType>& v, const Spointer& si);
    //! Check for steric hindrance for a site growth process
    bool checkHindrance(const Spointer& st) const;
    //! Check steric hindrance for phenyl addition reactions
    bool checkHindrancePhenyl(const Cpointer C_1) const;
	//! Check steric hindrance between a carbon and all other carbons in a PAH.
	bool checkHindrance_C_PAH(cpair coords) const;
	//! Check steric hindrance between a new carbon position and all other carbons in a PAH.
	bool checkHindrance_newC(Cpointer C_1) const;
	//! Check position between a carbon and internal carbons in a PAH.
	cpair checkHindrance_C_intPAH(cpair coords) const;
	//! Check steric hindrance between two carbons.
	bool checkHindrance_twoC(const Cpointer C_1, const Cpointer C_2) const;
	//!Gets distance between two carbons.
	double getDistance_twoC(const Cpointer C_1, const Cpointer C_2) const;
	//!Overload. Gets distance between two coordinates.
	double getDistance_twoC(const cpair C_1, const cpair C_2) const;
    //! Returns site iterator x steps after i
    Spointer moveIt(Spointer i, int x);
    //! Finds C atom with specific coordinates
    Cpointer findC(cpair coordinates);
	//! Find Site to the other side of a bridge
	Spointer findSite(Cpointer C_1);
	//! Finds C with distance less than 1.8 of current ZZ carbon
	Cpointer findThirdC(Cpointer C_1);

    // Write Process
    //! Creates a lone carbon atom
    Cpointer addC();
    //! Creates a new carbon atom attached next to C_1.
    Cpointer addC(Cpointer C_1, angletype angle1, angletype angle2, bondlength length, bool bulk=false);
	//! Creates a new carbon atom attached next to C1 using a vector to tell the position of the next carbon atom.
	Cpointer addC(Cpointer C_1, cpair direction, bondlength length, bool bulk=false);
	//! Overload routine. Passes a OpenBabel molecule so it is modified. 
	Cpointer addC(Cpointer C_1, cpair direction, bondlength length, OpenBabel::OBMol mol, bool bulk=false);
	//! Adds bond to OpenBabel structure between two carbons. 
	void addOBbond(Cpointer C_1, Cpointer C_2, OpenBabel::OBMol mol);
	//! Moves a carbon position .
	void moveC(Cpointer C_1, cpair newpos);
	//! Moves a carbon position .
	void moveC(Cpointer C_1, Cpointer Cprev, double new_distance);
	//! Move a carbon in z direction.
	void moveC_z(Cpointer C_1, double z_distance);
	//! Adds an R5 to the list of R5s and R7s
	void addR5internal(Cpointer C_1, Cpointer C_2, bool invert_dir=false);
	//! Removes an R5 from the list of R5s and R7s
	void removeR5internal(Cpointer C_1, Cpointer C_2);
	//! Removes an R7 from the list of R5s and R7s
	void removeR7internal(Cpointer C_1, Cpointer C_2);
	//! Return internal R5 associated to two carbons
	cpair findR5internal(Cpointer C_1, Cpointer C_2);
	//! Are the two carbon atoms members of an R5 with coordinates in R5Internal??
	bool isR5internal(Cpointer C_1, Cpointer C_2, bool invert_dir=false);
	//! Return coords of final position of an internal R5 based on two carbons
	cpair endposR5internal(Cpointer C_1, Cpointer C_2, bool invert_dir=false);
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
	//! Redraws 5 member rings
	void redrawR5(Spointer& st, Cpointer Carb1, Cpointer Carb2);
    //! Changes site type into another site type
    void convSiteType(Spointer& st, Cpointer Carb1, Cpointer Carb2, kmcSiteType t);
    //! Remove site
    void removeSite(Spointer& stt);
    //! Sets the number of counts of C and H
    void setCount(int CCount, int HCount);
    //! Add counts
    void addCount(int C_in, int H_in);
    //! For createPAH function: drawing type 0 sites
    Cpointer drawType0Site(Cpointer Cnow, int bulkC, double angle);
    //! For createPAH function: drawing type 1 sites
    Cpointer drawType1Site(Cpointer Cnow, int bulkC, kmcSiteType prevType, double angle);
    //! For createPAH function: drawing type 2 sites
    Cpointer drawType2Site(Cpointer Cnow, int bulkC, double angle);
	//! For createPAH function: drawing type 3 sites
    Cpointer drawType3Site(Cpointer Cnow, int bulkC, double angle);

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
	//! Overload function, forces char to C atom.
	void updateA(Cpointer C, char sp, bool flag);
	//! Overload function, forces char to C atom and defines growth vector.
	void updateA(Cpointer C, char sp, cpair gro_vec);
	//! Overload function, forces char to C atom and defines growth vector.
	void updateA(Cpointer C, char sp, cpair gro_vec, OpenBabel::OBMol mol);
    //! Overload function, updateA for all C from C_1 to C_2 inclusive
    void updateA(Cpointer C_1, Cpointer C_2, char spc);
    
    // PAH data structure to perform processes on
   PAHStructure* m_pah;

    //Cpointer NULLC;
};

}
}


#endif
