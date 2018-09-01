/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_gaspoint.cpp
  *
  * \brief        Implementation file for swp_kmc_gaspoint.h.
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Holds the data structure which contains information on the molecular structure
    of PAH accounted by the kMC model.

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

#include "swp_kmc_gaspoint.h"
#include "gpc_species.h"
#include "comostrings.h"

using namespace Sweep;
using namespace Sweep::KMC_ARS;

// Constructors and Destructor
//! Default Constructor
KMCGasPoint::KMCGasPoint():
		m_data(),
		m_gasprof(NULL),
		m_spnames(),
		m_prof_in()
{}

//! Constructor from a GasProfile object
KMCGasPoint::KMCGasPoint(Sweep::GasProfile& gasprof,
    const Sprog::SpeciesPtrVector& sptrv):
				m_data(),
				m_gasprof(NULL),
				m_spnames(),
				m_prof_in()
{
    m_gasprof = &gasprof;
    initData();
    std::vector<std::string> spname;
    for(size_t i=0; i<sptrv.size(); i++)
        spname.push_back(sptrv[i]->Name());
    for(int i=H2; i<=CO2; i++) {
        m_prof_in[i] = Strings::findinlist(m_spnames[i], spname);
    }
}

//! Copy Constructor
KMCGasPoint::KMCGasPoint(const KMCGasPoint &gp):
				m_data(),
				m_gasprof(NULL),
				m_spnames(),
				m_prof_in()
{
    *this = gp;
}
// Default Destructor
KMCGasPoint::~KMCGasPoint() {
}
//! Initialise data point
void KMCGasPoint::initData() {
    m_data.clear();
    m_spnames.clear();
    m_spnames.push_back("Time");
    m_spnames.push_back("T");
    m_spnames.push_back("H2");
    m_spnames.push_back("H");
    m_spnames.push_back("O2");
    m_spnames.push_back("OH");
    m_spnames.push_back("C2H2");
    m_spnames.push_back("C2H6");
    m_spnames.push_back("A1");
    m_spnames.push_back("H2O");
    m_spnames.push_back("CH4");
    m_spnames.push_back("CO");
    m_spnames.push_back("CO2");
    m_spnames.push_back("P");
    m_spnames.push_back("None");
    // set all to zero
    for(int i=0; i!=m_total-1; i++) {
        m_data.push_back(0);
    }
    // set 1 for no variables
    m_data.push_back(1);
}

//! Interpolate data
void KMCGasPoint::Interpolate(double t, double fact) {
    // get time point after t
    GasProfile::const_iterator j = LocateGasPoint(*m_gasprof, t);
    if(j == m_gasprof->begin() || j == m_gasprof->end()-1) {
        m_data[Time] = j->Time;
        m_data[T] = j->Gas.Temperature();
        m_data[P] = j->Gas.Pressure();
        for(int i=H2; i<(m_total-2); i++) { // exclude P & None
            m_data[i] = j->Gas.RawData()[m_prof_in[i]];
            m_data[i] *= fact;
        }
    }else {
        GasProfile::const_iterator i = j; --i;
        double dt_ij = j->Time - i->Time;
        double dt = t - i->Time;
        double wx = dt/dt_ij;
        double wy = 1-wx;
        m_data[Time] = t;
        m_data[T] = i->Gas.Temperature()*wy + j->Gas.Temperature()*wx;
        m_data[P] = i->Gas.Pressure()*wy + j->Gas.Pressure()*wx;
        for(int k=H2; k<(m_total-2); k++) { // exclude P & None
            m_data[k] = i->Gas.RawData()[m_prof_in[k]]*wy + j->Gas.RawData()[m_prof_in[k]]*wx;
            m_data[k] *= fact;
        }
    }
    ConvertMoleFrac();
}

//! Convert Mole frac to Conc
void KMCGasPoint::ConvertMoleFrac() {
    double factor = m_data[P]/(R*m_data[T]*1e6); // convert to mol/cm^3
    for(int k=H2; k<(m_total-2); k++)
        m_data[k] *= factor;
}

//! Accessing data
double KMCGasPoint::operator[](const int n) const {
    return m_data[n];
}
KMCGasPoint& KMCGasPoint::operator=(const KMCGasPoint& gp) {
    if(this != &gp) {
    m_data = gp.m_data;
    m_gasprof = gp.m_gasprof;
    m_spnames = gp.m_spnames;
    m_prof_in = gp.m_prof_in;
    return *this;
    } else return *this;
}
//! Get species names
std::vector<std::string> KMCGasPoint::SpNames() const {
    return m_spnames;
}
