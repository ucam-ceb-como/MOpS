/*!
  * \Author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_gasph.h
  *
  * \brief        Defines the gas phase for the kmc simulator
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Defines the gas phase for the kmc simulator, 
    stores the gas phase profiles (time, pressure, temperature, species concentrations)

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

#ifndef SWEEP_KMC_GASPH_H
#define SWEEP_KMC_GASPH_H

#include "swp_kmc_reaction.h"
#include "swp_kmc_jump_process.h"
#include "swp_kmc_pah_structure.h"
#include "swp_kmc_pah_process.h"
#include "swp_kmc_typedef.h"
#include "swp_kmc_structure_comp.h"
#include "swp_kmc_gaspoint.h"
#include "linear_interpolator.hpp"
#include "rng.h"
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
    // PAHStructure: Class which contains the data structure of the model
    class PAHStructure;

    const std::string Tlabel = "T[K]"; // Label for temperature in profiles
    const std::string Timelabel = "Time"; // Label for time in profiles
    const std::string Plabel = "P"; // Label for pressure in profiles

    // Standard vector for profile values
    //typedef Utils::LinearInterpolator<real, real> rInterp;
    typedef std::map<std::string, rvector> pmap;

    //! The KMC model
    /*!
     * Timestep and reaction chosen calculated using a kinetic monte-carlo
     * algorithm by Gillespie (1977).
    */
    class KMCGasph {
    public:
        // Constructors
        // Default
        KMCGasph();
        //! Copy Constructor
        KMCGasph(KMCGasph& m);
        //! From external file containing profiles and reactions
        KMCGasph(const std::string& filename); 
        //! Destructor
        virtual ~KMCGasph();

        // Write Processes
        //! Load profiles from an external file.
        void loadProfileCSV();
        //! Load processes from process list
        void loadProcesses(std::vector<JumpProcess*> (*jp)(const KMCGasPoint&));
        //! Set csv input file name
        void setCSVname(const std::string& filename);
        
        // KMC Algorithm and Read Processes
        //! Calculates each jump rates and stores in a vector
        //rvector calculateRates(std::vector<JumpProcess*>& p) const;
        //! Calculates total jump rates
        //real getTotalRate(const rvector& jrate);
        //! Choosing a reaction to be taken place, returns vector index of process
        int chooseReaction(const rvector& jrate, real (*rand_u01)()) const;
        //! Gets t_max
        real gett_stop();
        //! Creates profile interpolators
        void setInterpolators(std::vector<rvector>& values, // vector of values
            std::vector<std::string>& labels, // vector of profile labels
            std::string& Xp); // profile describing position of data (independent profile)
        //! Interpolates all profiles and calculates gas conc. (if interpolate_gas = true) at time tnow
        void interpolateProfiles(const real& tnow, bool interpolate_gas, real r_factor);
        //! Vector of all the jump processes associated with the model
        std::vector<JumpProcess*> jplist;
        
        KMCGasPoint m_gpoint;
        //! 
        void getIt();
        void initIt();
    private:
        //! CSV file name
        std::string m_csvfilename;
        //! CSV file reader class
        CSV_IO m_csvfile;
        //! Labels for profiles
        std::vector<std::string> m_profile_labels;
        //! values for profiles
        std::vector<rvector> m_profile_val;
        //! interpolators for profiles
        std::vector<Utils::LinearInterpolator<real, real> > m_gasProfile;
        //! interpolated profiles and calculated concentrations of species
        rvector m_interpolated;
        
        //! Stores all indexes for profiles in a vector
        //! order: (Time, T, H2, H, O2, OH, C2H2, C2H6, C6H6, H2O, CH4, CO, CO2, P)
        intvector m_profile_it;
    };
}
}

#endif
