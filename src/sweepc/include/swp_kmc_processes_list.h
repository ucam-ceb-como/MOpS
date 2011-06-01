/*!
  * \Author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_processes_list.h
  *
  * \brief      defines namespace which contains all jump processes and related functions
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Defines namespace which contains all jump processes and related functions

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
#ifndef SWP_KMC_PROCESSES_LIST_H
#define SWP_KMC_PROCESSES_LIST_H

#include "swp_kmc_gasph.h"
#include "swp_kmc_gaspoint.h"
#include "swp_kmc_jump_process.h"
#include "swp_kmc_pah_structure.h"
#include "swp_kmc_pah_process.h"
#include "swp_kmc_structure_comp.h"

#include <vector>

/* To add a new jump process, additions should be made at FIVE places:
    -swp_kmc_pah_structure.h: under protected: "Jump processes:", in the form of
        proc_StructureProc_kmcSiteType
    -swp_kmc_pah_structure.cpp: under protected function performProcess(..). Add
        the new jump process accordingly.
    -this file: in namespace JumpProcessList, under "Process list:", in the form of
        StructureProc_kmcSiteType : public JumpProcess
    -swp_kmc_processes_list.cpp: under "Process list (energy units in kcal)"
        defines the functions NEW_PROCESS::initialise() [elementary reactions,
        site type, structure change process and name], NEW_PROCESS::setRate(..)
        [jump rate calculation] and PAHProcess::proc_NEW_PROCESS [information
        for structure change].
    -swp_kmc_processes_list.cpp: in JumpProcessList::obtainJumpProcess(), under
        "Initialise all jump processes". The process can then be chosen to be included
        or not by specifying it under "Jump Processes included in the model:" in the
        same function
*/
namespace Sweep {
namespace KMC_ARS {
namespace JumpProcessList {
    //! Returns a vector of jump processes
    std::vector<JumpProcess*> obtainJumpProcess(const KMCGasPoint& gp);
    //! Calculates jump rate for each jump process, returns total rate
    real calculateRates(const KMCGasPoint& gp, 
        PAHProcess& st, 
        const real& t, 
        std::vector<JumpProcess*>& jp,
        rvector& rateV);

    //! Process list:
    class G6R_FE : public JumpProcess { //R6 growth on FE
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class G6R_AC : public JumpProcess { //R6 growth on AC
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class L6_BY6 : public JumpProcess { //BY6 closure
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class PH_benz : public JumpProcess { //phenyl addition
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class D6R_FE3 : public JumpProcess { //R6 desorption
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class O6R_FE3_O2 : public JumpProcess { //R6 oxidation at FE by O2
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class O6R_FE3_OH : public JumpProcess { //R6 oxidation at FE by OH
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class O6R_FE_HACA_O2 : public JumpProcess { //R6 oxidation at AC by O2
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class O6R_FE_HACA_OH : public JumpProcess { //R6 oxidation at AC by OH
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class G5R_ZZ : public JumpProcess { //R5 growth on ZZ
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class D5R_R5 : public JumpProcess { //R5 desorption
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class C6R_AC_FE3 : public JumpProcess { //R6 conversion to R5
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class C5R_RFE : public JumpProcess { //R5 conversion to R6 on FE
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class C5R_RAC : public JumpProcess { //R5 conversion to R6 on AC
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class M5R_RZZ : public JumpProcess { //R5 migration to neighbouring ZZ
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class C6R_BY5_FE3 : public JumpProcess { //R6 migration & conversion to R5 at BY5 (pathway 1)
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class C6R_BY5_FE3violi : public JumpProcess { //R6 migration & conversion to R5 at BY5 (pathway 2; violi)
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class L5R_BY5 : public JumpProcess { //BY5 closure
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class M6R_BY5_FE3 : public JumpProcess { //R6 desorption at bay -> pyrene
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
}
}
}

#endif
