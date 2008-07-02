/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the GasPhase class declared in the
    gpc_gasphase.h header file.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
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

#include "gpc_gasphase.h"
#include<vector>

using namespace Sprog;
using namespace Sprog::Thermo;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
GasPhase::GasPhase(void)
{
    /*
    m_precalc = false;
    m_Us_valid = false;
    m_Hs_valid = false;
    m_Ss_valid = false;
    m_Gs_valid = false;
    m_Cps_valid = false;
    m_Cvs_valid = false;
    m_Us_RT_valid = false;
    m_Hs_RT_valid = false;
    m_Ss_R_valid = false;
    m_Gs_RT_valid = false;
    m_Cps_R_valid = false;
    m_Cvs_R_valid = false;
    */
}

// Default constructor (public, requires species vector).
GasPhase::GasPhase(const SpeciesPtrVector &sp) : Mixture(sp)
{
}

// Default destructor.
GasPhase::~GasPhase(void)
{
}


// THERMODYNAMIC PROPERTIES.

void GasPhase::Us(Sprog::fvector &U) const
{
    CalcUs(*m_pT, U);
}

void GasPhase::Hs(Sprog::fvector &H) const
{
    CalcHs(*m_pT, H);
}

void GasPhase::Ss(Sprog::fvector &S) const
{
    CalcSs(*m_pT, S);
}

void GasPhase::Gs(Sprog::fvector &G) const
{
    CalcGs(*m_pT, G);
}

void GasPhase::Cps(Sprog::fvector &Cp) const
{
    CalcCps(*m_pT, Cp);
}

void GasPhase::Cvs(Sprog::fvector &Cv) const
{
    CalcCvs(*m_pT, Cv);
}


// THERMODYNAMIC PROPERTIES (DIMENSIONLESS).

void GasPhase::Us_RT(Sprog::fvector &U) const
{
    CalcUs_RT(*m_pT, U);
}

void GasPhase::Hs_RT(Sprog::fvector &H) const
{
    CalcHs_RT(*m_pT, H);
}

void GasPhase::Ss_R(Sprog::fvector &S) const
{
    CalcSs_R(*m_pT, S);
}

void GasPhase::Gs_RT(Sprog::fvector &G) const
{
    CalcGs_RT(*m_pT, G);
}

void GasPhase::Cps_R(Sprog::fvector &Cp) const
{
    CalcCps_R(*m_pT, Cp);
}

void GasPhase::Cvs_R(Sprog::fvector &Cv) const
{
    CalcCvs_R(*m_pT, Cv);
}

// BULK MIXTURE PROPERTIES.

// Calculates the bulk internal energies in current units.
real GasPhase::BulkU() const
{
    return CalcBulkU(*m_pT, &m_data[0], m_species->size());
}

// Calculates the bulk enthalpy in current units.
real GasPhase::BulkH() const
{
    return CalcBulkH(*m_pT, &m_data[0], m_species->size());
}

// Calculates the bulk entropy in current units.
real GasPhase::BulkS() const
{
    return CalcBulkS(*m_pT, &m_data[0], m_species->size());
}

// Calculates the bulk Gibbs free energies in current units.
real GasPhase::BulkG() const
{
    return CalcBulkG(*m_pT, &m_data[0], m_species->size());
}

// Calculates the mean molar heat capacity at const. P.
real GasPhase::BulkCp() const
{
    return CalcBulkCp(*m_pT, &m_data[0], m_species->size());
}

// Calculates the mean molar heat capacity at const. V.
real GasPhase::BulkCv() const
{
    return CalcBulkCv(*m_pT, &m_data[0], m_species->size());
}
