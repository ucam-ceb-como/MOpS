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
#include "swp_params.h"

namespace Sweep {
    namespace KMC_ARS {
        //! Enumeration of site types.
        enum kmcSiteType {
			//Reserve negatives for non-reactive sites
			//Preserve the number conventions between Principle and Combine site types
            // Principal Sites involving 6-member rings
			FE = 1, ZZ = 2, AC = 3, BY5 = 4, BY6 = 5, ACBL = 6, ACBR = 7,
			// Principal Sites involving 6-member rings that are non-reactive (due to hinderances)
			NFE = -1, NZZ = -2, NAC = -3, NBY5 = -4, NBY6 = -5,
			// Principal Sites involving 5-member rings. HER5 = One side embedded R5, ER5 = embedded R5
			R5 = 11, HER5 =12, ER5 =13,
            // Principal sites combining R6 and R5 sites
			// Site with R5 at one side
            RFE=21, RZZ=22, RAC=23, RBY5=24,
			// Site with ER5 at one side
			ERFE = 31, ERZZ = 32, ERAC = 33, ERBY5 = 34,
            // Sites with R5 at both sides
			RFER = 41, RZZR = 42, RACR = 43, RBY5R = 44,
			// Sites with ER5 at both sides
			ERFEER = 51, ERZZER = 52, ERACER = 53,
			// Sites with one R5 and one ER5
			RFEER = 61, RZZER = 62, RACER = 63,
			//End principal site types
			//Start combined site types
			//First, basic combined sites with R6s. FE3:FE with 2 FEs on both sides. FE_HACA: FE with non FE's on each side
			FE2 =1001, FE3 = 1002, FE_HACA = 1003, 
			//Combined sites with only R5
			//None at the moment, reserve 2000's for this
			//Combined sites with R5 and R6s
			ACR5 = 3001, CAC = 3002, CBY5 = 3003, CBY6 = 3004,
			//Second, combined sites with other principle and combined sites, AC_FE3:AC next to FE3
			AC_FE3 = 4001, BY5_FE3 = 4002 , RAC_FE3 = 4003, ER5_FE3 = 4004,
			//none
            None = 0,
            // Error ID
            Inv=-100000,
            // Any Site:
            any=100000,
            // Benzene addition sites:
            benz=100001
        };
        
        //! Enumeration of processes on data structure available.
        enum StructureProc {
            // Growth
            G6R=0, G5R=1,
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
            BENZENE_C=6, BENZENE_H=6, NAPHTHALENE_C=10, NAPHTHALENE_H=8, PYRENE_C=16, PYRENE_H=10, 
            BENZOPYRENE_C=20, BENZOPYRENE_H=12, CORONENE_C=24, CORONENE_H=12, BY5_C=18, BY5_H=12,
            TEST_STRUCT, NONE};

        //! typedefs used in the model
        typedef std::vector<double> rvector;
        typedef std::vector<int> intvector;
        //typedef int angletype;
        //typedef int coordtype;
        //typedef std::pair<coordtype, coordtype> cpair;
        typedef std::pair<int, int> intpair;
        
        
        //! Change site type to sitename
        std::string inline kmcSiteName(kmcSiteType site_type) {
            switch(site_type) {
                case FE: return "FE";
                case ZZ: return "ZZ";
                case AC: return "AC";
                case BY5: return "BY5";
                case BY6: return "BY6";
				case ACBL: return "ACBL";
				case ACBR: return "ACBR";

				case NFE: return "NFE";
				case NZZ: return "NZZ";
				case NAC: return "NAC";
				case NBY5: return "NBY5";
				case NBY6: return "NBY6";

                case R5: return "R5";
				case HER5: return "HER5";
				case ER5: return "ER5";

                case RFE: return "RFE";
                case RZZ: return "RZZ";
                case RAC: return "RAC";
                case RBY5: return "RBY5";

				case ERFE: return "ERFE";
				case ERZZ: return "ERZZ";
				case ERAC: return "ERAC";
				case ERBY5: return "ERBY5";

                case RFER: return "RFER";
                case RZZR: return "RZZR";
                case RACR: return "RACR";
				case RBY5R: return "RBY5R";

				case ERFEER: return "ERFEER";
				case ERZZER: return "ERZZER";
				case ERACER: return "ERACER";

				case RFEER: return "RFEER";
				case RZZER: return "RZZER";
				case RACER: return "RACER";

				case FE2: return "FE2";
                case FE3: return "FE3";
				case FE_HACA: return "FE_HACA";

				case ACR5: return "ACR5";
				case CAC: return "CAC";
				case CBY5: return "CBY5";
				case CBY6: return "CBY6";

                case AC_FE3: return "AC_FE3";

                case BY5_FE3: return "BY5_FE3";
                case RAC_FE3: return "RAC_FE3";
				case ER5_FE3: return "ER5_FE3";

                case None: return "None";
                case Inv: return "Invalid";
                case any: return "any";
                case benz: return "benz";
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
			else if (str == "ACBL") return ACBL;
			else if (str == "ACBR") return ACBR;

			if (str == "NFE") return NFE;
			else if (str == "NZZ") return NZZ;
			else if (str == "NAC") return NAC;
			else if (str == "NBY5") return NBY5;
			else if (str == "NBY6") return NBY6;

            else if(str == "R5") return R5;
			else if (str == "HER5") return HER5;
			else if (str == "ER5") return ER5;

            else if(str == "RFE") return RFE;
            else if(str == "RZZ") return RZZ;
            else if(str == "RAC") return RAC;
            else if(str == "RBY5") return RBY5;

			else if (str == "ERFE") return ERFE;
			else if (str == "ERZZ") return ERZZ;
			else if (str == "ERAC") return ERAC;
			else if (str == "ERBY5") return ERBY5;

            else if(str == "RFER") return RFER;
            else if(str == "RZZR") return RZZR;
            else if(str == "RACR") return RACR;
			else if (str == "RBY5R") return RBY5R;

			else if (str == "ERFEER") return ERFEER;
			else if (str == "ERZZER") return ERZZER;
			else if (str == "ERACER") return ERACER;

			else if (str == "RFEER") return RFEER;
			else if (str == "RZZER") return RZZER;
			else if (str == "RACER") return RACER;

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
			temp.push_back(ACBL);
			temp.push_back(ACBR);

			temp.push_back(NFE);
			temp.push_back(NZZ);
			temp.push_back(NAC);
			temp.push_back(NBY5);
			temp.push_back(NBY6);

            temp.push_back(R5);
			temp.push_back(HER5);
			temp.push_back(ER5);

            temp.push_back(RFE);
            temp.push_back(RZZ);
            temp.push_back(RAC);
            temp.push_back(RBY5);

			temp.push_back(ERFE);
			temp.push_back(ERZZ);
			temp.push_back(ERAC);
			temp.push_back(ERBY5);

            temp.push_back(RFER);
            temp.push_back(RZZR);
            temp.push_back(RACR);
			temp.push_back(RBY5R);

			temp.push_back(ERFEER);
			temp.push_back(ERZZER);
			temp.push_back(ERACER);

			temp.push_back(RFEER);
			temp.push_back(RZZER);
			temp.push_back(RACER);

			temp.push_back(FE2);
            temp.push_back(FE3);
			temp.push_back(FE_HACA);

			temp.push_back(ACR5);
			temp.push_back(CAC);
			temp.push_back(CBY5);
			temp.push_back(CBY6);

            temp.push_back(AC_FE3);
            temp.push_back(BY5_FE3);
            temp.push_back(RAC_FE3);
			temp.push_back(ER5_FE3);

            return temp;
        }
        //! Get a vector of all site types for phenyl addition
        std::vector<kmcSiteType> inline vectPHsites() {
            std::vector<kmcSiteType> temp;
            temp.push_back(FE);
            temp.push_back(ZZ);
            temp.push_back(AC);
            temp.push_back(RFE);
            temp.push_back(RZZ);
            return temp;
        }
    };
};
#endif
