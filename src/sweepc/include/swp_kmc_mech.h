/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_mech.h
  *
  * \brief        Defines mechanism for the KMC simulator
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Defines the mechanism used in the KMC simulator. Includes
    all elementary reactions and compound rate equation resulting
    from steady-state assumptions.

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

#ifndef SWEEP_KMC_MECH_H
#define SWEEP_KMC_MECH_H

//#include "swp_kmc_reaction.h"
#include "swp_kmc_jump_process.h"
//#include "swp_kmc_pah_structure.h"
#include "swp_kmc_pah_process.h"
#include "swp_kmc_typedef.h"
//#include "swp_kmc_structure_comp.h"
#include "swp_kmc_gaspoint.h"
#include "csv_io.h"

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cmath>


namespace Sweep {
namespace KMC_ARS {
    //! Forward declaration of classes
    // JumpProcess: Class which contains information on a jump process
    class JumpProcess;
    class PAHProcess;

    // PAHStructure: Class which contains the data structure of the model
    //class PAHStructure;
        
    typedef std::pair<JumpProcess*, int> ChosenProcess;

    //! The KMC model
    /*!
     * Timestep and reaction chosen calculated using a kinetic monte-carlo
     * algorithm by Gillespie (1977).
    */
    class KMCMechanism {
    public:
        // Constructors
        // Default
        KMCMechanism();

        //! Copy Constructor
        KMCMechanism(KMCMechanism& m);

        //! Destructor
        virtual ~KMCMechanism();

        // WRITE PROCESSES

        //! Load processes from process list
        void loadProcesses(std::vector<JumpProcess*> (*jp)());

        //! Choosing a reaction to be taken place, returns pointer to jump process
        //! and index of process in m_jplist
        ChosenProcess chooseReaction(rng_type &rng) const;

        //! Calculates jump rate for each jump process
        void calculateRates(const KMCGasPoint& gp, 
            PAHProcess& st, 
            const real& t);
        
        // DATA ACCESS

        //! Returns vector of jump processes
        std::vector<JumpProcess*> JPList() const;
        
        //! Returns vector of jump rates
        std::vector<real> Rates() const;

        //! Returns total rates
        real TotalRate() const;
    private:
        //! Vector of jump processes
        std::vector<JumpProcess*> m_jplist;
        //! Returns a vector of jump processes implemented in model
        std::vector<JumpProcess*> obtainJumpProcess();
        //! Checks if this mechanism object is a copy
        bool isACopy;
        //! Vector of jump rates
        std::vector<real> m_rates;
        //! Total rate
        real m_totalrate;
    };
        
    //! Process list:
    class G6R_FE : public Sweep::KMC_ARS::JumpProcess { //R6 growth on FE
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class G6R_AC : public Sweep::KMC_ARS::JumpProcess { //R6 growth on AC
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class L6_BY6 : public Sweep::KMC_ARS::JumpProcess { //BY6 closure
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class PH_benz : public Sweep::KMC_ARS::JumpProcess { //phenyl addition
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class D6R_FE3 : public Sweep::KMC_ARS::JumpProcess { //R6 desorption
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class O6R_FE3_O2 : public Sweep::KMC_ARS::JumpProcess { //R6 oxidation at FE by O2
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class O6R_FE3_OH : public Sweep::KMC_ARS::JumpProcess { //R6 oxidation at FE by OH
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class O6R_FE_HACA_O2 : public Sweep::KMC_ARS::JumpProcess { //R6 oxidation at AC by O2
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class O6R_FE_HACA_OH : public Sweep::KMC_ARS::JumpProcess { //R6 oxidation at AC by OH
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class G5R_ZZ : public Sweep::KMC_ARS::JumpProcess { //R5 growth on ZZ
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class D5R_R5 : public Sweep::KMC_ARS::JumpProcess { //R5 desorption
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class C6R_AC_FE3 : public Sweep::KMC_ARS::JumpProcess { //R6 conversion to R5
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class C5R_RFE : public Sweep::KMC_ARS::JumpProcess { //R5 conversion to R6 on FE
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class C5R_RAC : public Sweep::KMC_ARS::JumpProcess { //R5 conversion to R6 on AC
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class M5R_RZZ : public Sweep::KMC_ARS::JumpProcess { //R5 migration to neighbouring ZZ
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class C6R_BY5_FE3 : public Sweep::KMC_ARS::JumpProcess { //R6 migration & conversion to R5 at BY5 (pathway 1)
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class C6R_BY5_FE3violi : public Sweep::KMC_ARS::JumpProcess { //R6 migration & conversion to R5 at BY5 (pathway 2; violi)
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class L5R_BY5 : public Sweep::KMC_ARS::JumpProcess { //BY5 closure
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class M6R_BY5_FE3 : public Sweep::KMC_ARS::JumpProcess { //R6 desorption at bay -> pyrene
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
    class O6R_FE2_OH : public Sweep::KMC_ARS::JumpProcess { //R6 desorption at bay -> pyrene
    public:
        real setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        real setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const real& time_now*/);
        void initialise();
    };
}

}

#endif
