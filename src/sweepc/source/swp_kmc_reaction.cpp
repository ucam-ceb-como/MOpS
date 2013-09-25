/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_reaction.cpp
  *
  * \brief        Implementation file for swp_kmc_reaction.h
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

#include "swp_kmc_reaction.h"
#include <cmath>
#include <string>
#include "swp_params.h"

using namespace std;
using namespace Sweep;
using namespace Sweep::KMC_ARS;


//! Default Constructor
Reaction::Reaction() {
    // set all Arrhenius parameters into 0
    m_A = 0;
    m_n = 0;
    m_E = 0;
    m_r_species = 0;
}
//! Constructor from Arrhenius constants
Reaction::Reaction(double A_c, double n_c, double E_c, int species_c) {
    m_A = A_c;
    m_n = n_c;
    m_E = E_c;
    m_r_species = species_c;
}

//! Destructor
Reaction::~Reaction() {
}

//! Get Rate
double Reaction::getRate(const KMCGasPoint& gp/*, const double& t*/) const {
	double T(gp[gp.T]);

	// Putting the nlog(T) into the exp is faster than using pow separately
	// (for n != int, n < 1 or n > 3)
	return m_A * gp[m_r_species] * exp(m_n * log(T) - m_E / (RCAL * T));
}
