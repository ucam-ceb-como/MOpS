/*!
  * \Author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_gaspoint.h
  *
  * \brief        Defines the gas phase species data point for the kmc simulator
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Defines the gas phase species for the kmc simulator, 
    stores current datapoint of the gas phase profiles (time, pressure, temperature, species concentrations)

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

#ifndef SWEEP_KMC_GASPOINT_H
#define SWEEP_KMC_GASPOINT_H

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "swp_kmc_typedef.h"
#include "swp_gas_profile.h"

namespace Sweep {
namespace KMC_ARS {
        class KMCGasPoint {
        public:
            //! Default Constructor
            KMCGasPoint();
            //! Constructor from a GasProfile object
            KMCGasPoint(Sweep::GasProfile& gasprof,
                const Sprog::SpeciesPtrVector& sptrv);
            //! Copy Constructor
            KMCGasPoint(const KMCGasPoint& gp);
            //! Default Destructor
            virtual ~KMCGasPoint();
            //! Initialise data point
            void initData();
            //! Interpolate data
            void Interpolate(real t, real fact=1);
            //! Convert Mole frac to Conc
            void ConvertMoleFrac();
            
            //! Profile number and column index
            const int    Time;
            const int    T;
            const int    H2;
            const int    H;
            const int    O2;
            const int    OH;
            const int    C2H2;
            const int    C2H6;
            const int    C6H6;
            const int    H2O;
            const int    CH4;
            const int    CO;
            const int    CO2;
            const int    P;
            const int    None;

            const int    m_total;
            //! Accessing data
            real operator[](const int n) const;
            KMCGasPoint& operator=(const KMCGasPoint& gp);
            //! Get species names
            std::vector<std::string> SpNames() const;

        private:
            //! Datapoint for each variable (arranged according to above const int order)
            std::vector<real> m_data;
            Sweep::GasProfile* m_gasprof;
            std::vector<std::string> m_spnames;

            //! Map profile number with column index
            std::map<int, size_t> m_prof_in;
        };
}
}

#endif
