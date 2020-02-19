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
            const double& t);
        
        // DATA ACCESS

        //! Returns vector of jump processes
        const std::vector<JumpProcess*>& JPList() const;
        
        //! Returns vector of jump rates
        const std::vector<double>& Rates() const;

        //! Returns total rates
        double TotalRate() const;
    private:
        //! Vector of jump processes
        std::vector<JumpProcess*> m_jplist;
        //! Returns a vector of jump processes implemented in model
        std::vector<JumpProcess*> obtainJumpProcess();
        //! Checks if this mechanism object is a copy
        bool isACopy;
        //! Vector of jump rates
        std::vector<double> m_rates;
        //! Total rate
        double m_totalrate;
    };
        
    //! Process list:

    //! ID1.
    class G6R_FE : public Sweep::KMC_ARS::JumpProcess { //R6 growth on FE
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID2.
    class G6R_AC : public Sweep::KMC_ARS::JumpProcess { //R6 growth on AC
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID3.
    class L6_BY6 : public Sweep::KMC_ARS::JumpProcess { //BY6 closure
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };
    class PH_benz : public Sweep::KMC_ARS::JumpProcess { //phenyl addition
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID5.
    class D6R_FE3 : public Sweep::KMC_ARS::JumpProcess { //R6 desorption
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID6.
    class O6R_FE3_O2 : public Sweep::KMC_ARS::JumpProcess { //R6 oxidation at FE by O2
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID7.
    class O6R_FE3_OH : public Sweep::KMC_ARS::JumpProcess { //R6 oxidation at FE by OH
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID8.
    class O6R_FE_HACA : public Sweep::KMC_ARS::JumpProcess { //R6 oxidation at AC by O2
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID9.
    class O6R_FE_HACA_OH : public Sweep::KMC_ARS::JumpProcess { //R6 oxidation at AC by OH
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID10.
    class G5R_ZZ : public Sweep::KMC_ARS::JumpProcess { //R5 growth on ZZ
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID11.
    class D5R_R5 : public Sweep::KMC_ARS::JumpProcess { //R5 desorption
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID12.
    class C6R_AC_FE3 : public Sweep::KMC_ARS::JumpProcess { //R6 conversion to R5
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID13.
    class C5R_RFE : public Sweep::KMC_ARS::JumpProcess { //R5 conversion to R6 on FE
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID14.
    class C5R_RAC : public Sweep::KMC_ARS::JumpProcess { //R5 conversion to R6 on AC
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };
    //! ID15.
    class M5R_RZZ : public Sweep::KMC_ARS::JumpProcess { //R5 migration to neighbouring ZZ
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID16.
    class C6R_BY5_FE3 : public Sweep::KMC_ARS::JumpProcess { //R6 migration & conversion to R5 at BY5 (pathway 1)
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID17.
    class C6R_BY5_FE3violi : public Sweep::KMC_ARS::JumpProcess { //R6 migration & conversion to R5 at BY5 (pathway 2; violi)
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID18.
    class L5R_BY5 : public Sweep::KMC_ARS::JumpProcess { //BY5 closure
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID19.
    class M6R_BY5_FE3 : public Sweep::KMC_ARS::JumpProcess { //R6 desorption at bay -> pyrene
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID20.
    class O6R_FE2_side : public Sweep::KMC_ARS::JumpProcess { //R6 desorption at bay -> pyrene
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

    //! ID21.
    class O6R_FE2_top : public Sweep::KMC_ARS::JumpProcess { //R6 desorption at bay -> pyrene
    public:
        double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
        void initialise();
    };

	//! ID22.
	class D6R_FE_AC : public Sweep::KMC_ARS::JumpProcess { //R6 Desorption from FE to AC site.
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID22.
	class B6R_ACR5 : public Sweep::KMC_ARS::JumpProcess { //R6 desorption at bay -> pyrene
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID23.
	class M5R_ACR5_ZZ : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID24.
	class G6R_RZZ : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID25.
	class G6R_RFER : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID26.
	class G6R_R5 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID27.
	class L6_RBY5 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID28.
	class L6_RACR : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID29.
	class G5R_RFE : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID30.
	class C6R_RAC_FE3 : public Sweep::KMC_ARS::JumpProcess { //R6 migration & conversion to R5 at BY5 (pathway 1)
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID31.
	class C6R_RAC_FE3violi : public Sweep::KMC_ARS::JumpProcess { //R6 migration & conversion to R5 at BY5 (pathway 2; violi)
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID32.
	class M6R_RAC_FE3 : public Sweep::KMC_ARS::JumpProcess { //R6 desorption at bay -> pyrene
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID34.
	class MR5_R6 : public Sweep::KMC_ARS::JumpProcess { //R5 flip with R6
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};

	//! ID35.
	class GR7_R5R6AC : public Sweep::KMC_ARS::JumpProcess { //R7 growth at embedded-obstructed R5
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID36.
	class GR7_FEACR5 : public Sweep::KMC_ARS::JumpProcess { //R7 growth at embedded-obstructed R5-2
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID37.
	class G6R_R5R6ZZ : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID38.
	class L7_ACACR5 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID39.
	class G6R_R5R6FER : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID40.
	class G6R_R5R6FER5R6 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID41.
	class L7_FEZZACR5 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID42.
	class C5R_RZZR : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID43.
	class C5R_R5R6ZZR : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID44.
	class L6_R5R6BY5 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID45.
	class L6_R5R6ACR : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID46.
	class L6_R5R6ACR5R6 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID47.
	class L6_ZZACR5 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID48.
	class L6_R5FEACR5 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID49.
	class L6_FEACR5FE : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID50.
	class L6_R5ACR5R5 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID51.
	class L7_R5ZZACR5 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	//! ID52.
	class L6_ACR5R5R6 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID53.
	class L7_ACR5R5R6ZZ : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID54.
	class A_CH3 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID55.
	class D_CH3 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID56.
	class O5R_R5R6 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID57.
	class O5R_R5R6ZZ : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID58.
	class O5R_R5R6AC : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID59.
	class O5R_R5R6BY5 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID60.
	class O5R_R5R6FER : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID61.
	class O5R_R5R6ZZR : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID62.
	class O5R_R5R6ACR : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID63.
	class O5R_R5R6FER5R6 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID64.
	class O5R_R5R6ZZR5R6 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
	
	//! ID65.
	class O5R_R5R6ACR5R6 : public Sweep::KMC_ARS::JumpProcess {
	public:
		double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
		void initialise();
	};
}

}

#endif
