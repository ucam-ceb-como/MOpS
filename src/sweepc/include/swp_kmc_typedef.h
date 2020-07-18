/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_typedef.h
  *
  * \brief        Defines typedefs and types used by the kmc model.
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Defines typedefs and types used by the kmc model.

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

#ifndef SWP_KMC_TYPEDEF_H
#define SWP_KMC_TYPEDEF_H

#include <vector>
#include <string>
#include <tuple>
#include "swp_params.h"

namespace Sweep {
    namespace KMC_ARS {
        //! Enumeration of site types.
        enum kmcSiteType {
			// Principal Sites, uncombined basic types
			FE = 0, ZZ = 1, AC = 2, BY5 = 3, BY6 = 4, R5 = 100,
			// Principal sites with R5 at one side
			RFE = 101, RZZ = 102, RAC = 103, RBY5 = 104,
			// Principal sites with R5 at both sides
			RFER = 202, RZZR = 203, RACR = 204,
			// Combined Sites: FE3:FE with 2 FEs on both sides, AC_FE3:AC next to FE3
			// FE_HACA:FE with non-FEs on both sides
			FE3 = 13, AC_FE3 = 14, FE_HACA = 15, BY5_FE3 = 16,
			// Combined Sites: FE2:FE with FE at one side, but not an FE3
			FE2 = 17, R5R7 = 18,
			// Armchair site formed by the sides of a 6-member aromatic ring, followed by a 5-member ring, then a 6-member ring
			//ACR5 = 18,
			//
			RAC_FE3 = 19, MIGR = 28, MIGR2 = 38, R5R6_MIGR = 48,
			// R5 shoulder with shoulder with an R6
			R5R6 = 501, R5R6ZZ = 502, R5R6AC = 503, R5R6BY5 = 504,
			// Double R5, one shouldered with R6
			R5R6FER = 602, R5R6ZZR = 603, R5R6ACR = 604,
			// Principal sites with Embedded R5R6 to both sides
			R5R6FER5R6 = 1002, R5R6ZZR5R6 = 1003, R5R6ACR5R6 = 1004,
			//eBY5 = 20,
			// Armchair site formed by the sides of a 6-member aromatic ring, followed by a 5-member ring, then a 6-member ring
			ACR5 = 2002, FEACR5 = 2003, ZZACR5 = 2004, FEACR5FE = 2014, ACACR5 = 2005, FEZZACR5 = 2015,
			// Armchair with R5 and R5s
			R5ACR5 = 2103, R5FEACR5 = 2104, ACR5R5R6 = 2114, R5ZZACR5 = 2105, ACR5R5R6ZZ = 2115, R5ACR5R5 = 2204, ACR5RFER = 2205, SPIRAL = 9999,
			// Combined Sites
			None,
			// Error ID
			Inv = -1,
			// Any Site:
			any = 1000,
			// Benzene addition sites:
			benz = 1001,
			// Methyl sites:
			Methyl = 6666
			/*
			// Principal Sites, uncombined basic types
            FE=0,ZZ=1,AC=2,BY5=3,BY6=4,R5=5,
            // Principal sites with R5 at one side
            RFE=6, RZZ=7, RAC=8, RBY5=9,
            // Principal sites with R5 at both sides
            RFER=10, RZZR=11, RACR=12,
            // Combined Sites: FE3:FE with 2 FEs on both sides, AC_FE3:AC next to FE3
            // FE_HACA:FE with non-FEs on both sides
            FE3=13, AC_FE3=14, FE_HACA=15, BY5_FE3=16,
            // Combined Sites: FE2:FE with FE at one side, but not an FE3
            FE2=17,
			// Armchair site formed by the sides of a 6-member aromatic ring, followed by a 5-member ring, then a 6-member ring
			//ACR5 = 18,
			//
			RAC_FE3 = 19,
			// R5 shoulder with shoulder with an R6
			R5R6 = 20, R5R6_emb = 21, BY5R5 = 22, BY6R5 = 23,
			// Double R5 shoulder with R6
			R5R6R5 = 25, R5R6_embR5 = 26, R5BY5R5 = 27,
			//eBY5 = 20,
			// Armchair site formed by the sides of a 6-member aromatic ring, followed by a 5-member ring, then a 6-member ring
			ACR5 = 30, R6ACR5 = 31, R6ACR5R6 = 32, R5ACR5 = 35, R6ACR5R5 = 36, R5ACR5R5 = 40,
			// Principal sites with Embedded R5R6 at one side
			R5R6FE = 56, R5R6ZZ = 57, R5R6AC = 58, R5R6BY5 = 59,
			// Principal sites with Embedded R5R6 at one side and R5 to the other
			R5R6FER = 62, R5R6ZZR = 63, R5R6ACR = 64,
			// Combined Sites
            None,
            // Error ID
            Inv=-1,
            // Any Site:
            any=100, 
            // Benzene addition sites:
            benz=101*/
        };
        
        //! Enumeration of processes on data structure available.
        enum StructureProc {
            // Growth
            G6R=0, G5R=1, GR7=13,
            // Desorption
            D6R=2, D5R=3,
            // Migration
            M6R=4, M5R=5,
            // Oxidation
            O6R=6, O5R=7,
            // Conversion
            C6R=8, C5R=9,
            // Closure
            L6=10, L5=11,
            // Phenyl Addition
            PH=12,
            // Error
            err=-1
        };
        
        //! Enumeration of starting structures
        enum StartingStructure {
            BENZENE_C=6, BENZENE_H=6, TOLUENE_C=7, TOLUENE_H=8, NAPHTHALENE_C=10, NAPHTHALENE_H=8, PYRENE_C=16, PYRENE_H=10, 
            METHYLPYRENE_C=17, METHYLPYRENE_H=12, MPHENANTHRENER_C=15, MPHENANTHRENER_H=9, BENZOPYRENE_C=20, BENZOPYRENE_H=12, CORONENE_C=24, CORONENE_H=12,
            TEST_STRUCT, NONE};

        //! typedefs used in the model
        typedef std::vector<double> rvector;
        typedef std::vector<int> intvector;
        typedef double angletype; //Changed to allow changing the side length of a R5 (gl413)
		typedef double bondlength; //Changed to allow changing the side length of a R5 (gl413)
        typedef double coordtype; //Changed to allow changing the side length of a R5 (gl413)
        typedef std::tuple<coordtype, coordtype, coordtype> cpair;
        typedef std::pair<int, int> intpair;
        
        
        //! Change site type to sitename
        std::string inline kmcSiteName(kmcSiteType site_type) {
            switch(site_type) {
                case FE: return "FE";
                case ZZ: return "ZZ";
                case AC: return "AC";
                case BY5: return "BY5";
                case BY6: return "BY6";
                case R5: return "R5";
                case RFE: return "RFE";
                case RZZ: return "RZZ";
                case RAC: return "RAC";
                case RBY5: return "RBY5";
                case RFER: return "RFER";
                case RZZR: return "RZZR";
                case RACR: return "RACR";
                case FE3: return "FE3";
                case AC_FE3: return "AC_FE3";
                case FE_HACA: return "FE_HACA";
                case BY5_FE3: return "BY5_FE3";
                case FE2: return "FE2";
                case R5R7: return "R5R7";
				case R5R6: return "R5R6";
				case R5R6ZZ: return "R5R6ZZ";
				case R5R6AC: return "R5R6AC";
				case R5R6BY5: return "R5R6BY5";
				case R5R6FER: return "R5R6FER";
				case R5R6ZZR: return "R5R6ZZR";
				case R5R6ACR: return "R5R6ACR";
				case R5R6FER5R6: return "R5R6FER5R6";
				case R5R6ZZR5R6: return "R5R6ZZR5R6";
				case R5R6ACR5R6: return "R5R6ACR5R6";
				case ACR5: return "ACR5";
				case FEACR5: return "FEACR5";
				case ZZACR5: return "ZZACR5";
				case FEACR5FE: return "FEACR5FE";
				case ACACR5: return "ACACR5";
				case FEZZACR5: return "FEZZACR5";
				case R5ACR5: return "R5ACR5";
				case R5FEACR5: return "R5FEACR5";
				case ACR5R5R6: return "ACR5R5R6";
				case R5ZZACR5: return "R5ZZACR5";
				case ACR5R5R6ZZ: return "ACR5R5R6ZZ";
				case R5ACR5R5: return "R5ACR5R5";
				case ACR5RFER: return "ACR5RFER";
				case RAC_FE3: return "RAC_FE3";
                case MIGR: return "MIGR";
                case MIGR2: return "MIGR2";
                case R5R6_MIGR: return "R5R6_MIGR";
				case SPIRAL: return "SPIRAL";
                case None: return "None";
                case Inv: return "Invalid";
                case any: return "any";
                case benz: return "benz";
				case Methyl: return "Methyl";
            }
            return "ERROR";
        }

        //! Change site name to site type (only for principal site types).
        //! Returns Inv if site type is invalid
        kmcSiteType inline kmcSiteType_str(std::string str) {
            if(str == "FE") return FE;
            else if(str == "ZZ") return ZZ;
            else if(str == "AC") return AC;
            else if(str == "BY5") return BY5;
            else if(str == "BY6") return BY6;
            else if(str == "R5") return R5;
            else if(str == "RFE") return RFE;
            else if(str == "RZZ") return RZZ;
            else if(str == "RAC") return RAC;
            else if(str == "RBY5") return RBY5;
            else if(str == "RFER") return RFER;
            else if(str == "RZZR") return RZZR;
            else if(str == "RACR") return RACR;
			else if (str == "R5R6") return R5R6;
			else if (str == "R5R6ZZ") return R5R6ZZ;
			else if (str == "R5R6AC") return R5R6AC;
			else if (str == "R5R6BY5") return R5R6BY5;
			else if (str == "R5R6FER") return R5R6FER;
			else if (str == "R5R6ZZR") return R5R6ZZR;
			else if (str == "R5R6ACR") return R5R6ACR;
			else if (str == "R5R6FER5R6") return R5R6FER5R6;
			else if (str == "R5R6ZZR5R6") return R5R6ZZR5R6;
			else if (str == "R5R6ACR5R6") return R5R6ACR5R6;
			else if (str == "ACR5") return ACR5;
			else if (str == "FEACR5") return FEACR5;
			else if (str == "ZZACR5") return ZZACR5;
			else if (str == "FEACR5FE") return FEACR5FE;
			else if (str == "ACACR5") return ACACR5;
			else if (str == "FEZZACR5") return FEZZACR5;
			else if (str == "R5ACR5") return R5ACR5;
			else if (str == "R5FEACR5") return R5FEACR5;
			else if (str == "ACR5R5R6") return ACR5R5R6;
			else if (str == "R5ZZACR5") return R5ZZACR5;
			else if (str == "ACR5R5R6ZZ") return ACR5R5R6ZZ;
			else if (str == "R5ACR5R5") return R5ACR5R5;
			else if (str == "ACR5RFER") return ACR5RFER;
			else if (str == "SPIRAL") return SPIRAL;
			else if (str == "Methyl") return Methyl;
            return Inv;
        }
        //! Get a vector of all site types
        std::vector<kmcSiteType> inline vectSiteType() {
            std::vector<kmcSiteType> temp;
            temp.push_back(FE);
            temp.push_back(ZZ);
            temp.push_back(AC);
            temp.push_back(BY5);
            temp.push_back(BY6);
            temp.push_back(R5);
            temp.push_back(RFE);
            temp.push_back(RZZ);
            temp.push_back(RAC);
            temp.push_back(RBY5);
            temp.push_back(RFER);
            temp.push_back(RZZR);
            temp.push_back(RACR);
            temp.push_back(FE2);
            temp.push_back(FE3);
            temp.push_back(AC_FE3);
            temp.push_back(RAC_FE3);
            temp.push_back(MIGR);
            temp.push_back(MIGR2);
            temp.push_back(R5R6_MIGR);
            temp.push_back(FE_HACA);
            temp.push_back(BY5_FE3);
            temp.push_back(R5R7);
			temp.push_back(R5R6);
			temp.push_back(R5R6ZZ);
			temp.push_back(R5R6AC);
			temp.push_back(R5R6BY5);
			temp.push_back(R5R6FER);
			temp.push_back(R5R6ZZR);
			temp.push_back(R5R6ACR);
			temp.push_back(R5R6FER);
			temp.push_back(R5R6ZZR);
			temp.push_back(R5R6ACR);
			temp.push_back(ACR5);
			temp.push_back(FEACR5);
			temp.push_back(ZZACR5);
			temp.push_back(FEACR5FE);
			temp.push_back(ACACR5);
			temp.push_back(FEZZACR5);
			temp.push_back(R5ACR5);
			temp.push_back(R5FEACR5);
			temp.push_back(ACR5R5R6);
			temp.push_back(R5ZZACR5);
			temp.push_back(ACR5R5R6ZZ);
			temp.push_back(R5ACR5R5);
			temp.push_back(ACR5RFER);
			temp.push_back(SPIRAL);
			temp.push_back(Methyl);
            return temp;
        }
        //! Get a vector of all site types for phenyl addition.
		//Assumes only additions to R6 sites. Discourages additions to restricted sites.
        std::vector<kmcSiteType> inline vectPHsites() {
            std::vector<kmcSiteType> temp;
            temp.push_back(FE);
            temp.push_back(ZZ);
            temp.push_back(AC);
            temp.push_back(RFE);
            temp.push_back(RZZ);
			temp.push_back(R5R6ZZ);
            return temp;
        }
    };
};
#endif
