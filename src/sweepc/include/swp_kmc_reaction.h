/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_reaction.h
  *
  * \brief        Defines the Reaction class used by the kMC model.
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    The Reaction class contains the rate constants for individual PAH growth reactions
    as described by Celnik (2008) and Raj (2009) for the calculation of individual
    reaction rates.

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

#ifndef SWP_KMC_REACTION_H
#define SWP_KMC_REACTION_H

#include <string>
#include <vector>
#include <map>
#include "swp_kmc_gaspoint.h"
#include "swp_kmc_typedef.h"
#include "swp_params.h"

namespace Sweep {
    namespace KMC_ARS {
        //class KMCGasph;
        class PAHStructure;
        class Reaction {
        public:
            //typedef Utils::LinearInterpolator<double, double> Interpolator;
            //! Default Constructor
            Reaction();
            //! Constructor from Arrhenius constants
            Reaction(double A_c, double n_c, double E_c, int species_c);
            //! Destructors
            virtual ~Reaction();//default

            // Read Processes
            //! Calculate rate of equation using gas profiles in model at time t
            double getRate(const KMCGasPoint& gp/*, const double& t*/) const;

        private:
            //! Arrhenius Constants arranged in order: A, n, E
            double m_A;
            double m_n;
            double m_E;
            // Species name
            int m_r_species;
        };
    }
}
#endif
