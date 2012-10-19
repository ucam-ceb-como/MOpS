/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_jump_process.h
  *
  * \brief        Defines the JumpProcess class used by the kMC model.
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    The JumpProcess class is a base class which contains the rate constants for a jump process 
    related to several Reaction classes as described by Celnik (2008) and Raj (2009) for the 
    calculation of reaction rates of PAH processes with intermediates assumed as steady 
    states.

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

#ifndef SWP_KMC_JUMP_PROCESS_H
#define SWP_KMC_JUMP_PROCESS_H

#include "swp_kmc_typedef.h"
#include "swp_kmc_reaction.h"
//#include "swp_kmc_mech.h"
#include "swp_kmc_gaspoint.h"
//#include "swp_kmc_structure_comp.h"
//#include "swp_kmc_pah_structure.h"
//#include "swp_kmc_pah_process.h"
#include <vector>
#include <string>

namespace Sweep {
    namespace KMC_ARS {
        class PAHProcess;
        class Reaction;
        // Pointer to member function in PAHProcess
        //typedef void(Sweep::KMC_ARS::PAHProcess::*pt2Proc)(Spointer&,Cpointer,Cpointer);
        class JumpProcess {
        public:
            //! Constructor
            JumpProcess();
            //! Copy contructor
            JumpProcess(const JumpProcess& p);
            //! Destructors
            virtual ~JumpProcess(); //default
            

            // Set Processes
            //! Loads elementary reaction details (defined in derived classes), kmcSiteType and StructureProc related to it
            virtual void initialise();
            //! Adds reaction
            void addReaction(std::vector<Sweep::KMC_ARS::Reaction>& rxnv, const Sweep::KMC_ARS::Reaction& rxn);
            //! Calculate rates of each elementary reaction
            void calculateElemRxnRate(std::vector<Sweep::KMC_ARS::Reaction>& rxnv, const KMCGasPoint& gp/*, const double t_now*/);
            //! Calculates jump process rates and store (for Pressures 0.0267, 0.12 & 1 atm; defined in derived classes)
            virtual double setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
            virtual double setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);
            virtual double setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/);

            // Read Processes
            //! Jump Rate
            double getRate() const;
            //! Gets site type associated with the jump process
            kmcSiteType getSiteType() const;
            //! Returns name of process
            std::string getName() const;
            //! Returns process ID
            int getID() const;
            //! Returns references to elementary reactions vectors
            std::vector<Sweep::KMC_ARS::Reaction>& getVec0p0267();
            std::vector<Sweep::KMC_ARS::Reaction>& getVec0p12();
            std::vector<Sweep::KMC_ARS::Reaction>& getVec1();
        protected:
            //! Site type associated to the jump process
            kmcSiteType m_sType;
            //! List of elementary reactions making up the jump process
            //! for pressures P=0.027, 0.12 and 1 atm
            std::vector<Sweep::KMC_ARS::Reaction> m_rxnvector0p0267;
            std::vector<Sweep::KMC_ARS::Reaction> m_rxnvector0p12;
            std::vector<Sweep::KMC_ARS::Reaction> m_rxnvector1;
            //! Vector which stores rates of elementary reactions
            std::vector<double> m_r;
            //! Rate of reaction
            double m_rate;
            //! Name of process
            std::string m_name;
            //! ID of process
            int m_ID;
        };
    }
}

#endif
