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
            // 
            // Combined Sites
            None,
            // Error ID
            Inv=-1,
            // Any Site:
            any=100,
            // Benzene addition sites:
            benz=101
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
            BENZOPYRENE_C=20, BENZOPYRENE_H=12, CORONENE_C=24, CORONENE_H=12,
            TEST_STRUCT, NONE};

        //! typedefs used in the model
        typedef std::vector<double> rvector;
        typedef std::vector<int> intvector;
        typedef int angletype;
        typedef int coordtype;
        typedef std::pair<coordtype, coordtype> cpair;
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
            else if(str == "R5") return R5;
            else if(str == "RFE") return RFE;
            else if(str == "RZZ") return RZZ;
            else if(str == "RAC") return RAC;
            else if(str == "RBY5") return RBY5;
            else if(str == "RFER") return RFER;
            else if(str == "RZZR") return RZZR;
            else if(str == "RACR") return RACR;
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
            temp.push_back(FE_HACA);
            temp.push_back(BY5_FE3);
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
