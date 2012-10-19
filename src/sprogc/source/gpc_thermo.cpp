/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ThermoInterface class declared in the
    gpc_thermo.h header file.

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

#include "gpc_thermo.h"
#include "gpc_params.h"

using namespace Sprog;
using namespace Sprog::Thermo;


// CONSTRUCTORS AND DESTRUCTORS.

ThermoInterface::~ThermoInterface(void)
{
}


// INTERNAL ENERGY.

double ThermoInterface::CalcBulkU(double T, const Sprog::fvector &x) const
{
    return CalcBulkU(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkU(double T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &U) const
{
    return CalcBulkU(T, &x[0], x.size(), U);
}

double ThermoInterface::CalcBulkU(double T, 
                                const double *const x, 
                                unsigned int n) const
{
    return CalcBulkU(T, x, n, m_tmpvec);
}


// DIMENSIONLESS INTERNAL ENERGY.

double ThermoInterface::CalcBulkU_RT(double T, const Sprog::fvector &x) const
{
    return CalcBulkU_RT(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkU_RT(double T, 
                                   const Sprog::fvector &x, 
                                   Sprog::fvector &U) const
{
    return CalcBulkU_RT(T, &x[0], x.size(), U);
}

double ThermoInterface::CalcBulkU_RT(double T, 
                                   const double *const x, 
                                   unsigned int n) const
{
    return CalcBulkU_RT(T, x, n, m_tmpvec);
}


// ENTHALPY.

double ThermoInterface::CalcBulkH(double T, const Sprog::fvector &x) const
{
    return CalcBulkH(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkH(double T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &H) const
{
    return CalcBulkH(T, &x[0], x.size(), H);
}

double ThermoInterface::CalcBulkH(double T, 
                                const double *const x, 
                                unsigned int n) const
{
    return CalcBulkH(T, x, n, m_tmpvec);
}


// DIMENSIONLESS ENTHALPY.

double ThermoInterface::CalcBulkH_RT(double T, const Sprog::fvector &x) const
{
    return CalcBulkH_RT(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkH_RT(double T, const Sprog::fvector &x, 
                                   Sprog::fvector &H) const
{
    return CalcBulkH_RT(T, &x[0], x.size(), H);
}

double ThermoInterface::CalcBulkH_RT(double T, 
                                   const double *const x, 
                                   unsigned int n) const
{
    return CalcBulkH_RT(T, x, n, m_tmpvec);
}


// ENTROPY.

double ThermoInterface::CalcBulkS(double T, const Sprog::fvector &x) const
{
    return CalcBulkS(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkS(double T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &S) const
{
    return CalcBulkS(T, &x[0], x.size(), S);
}

double ThermoInterface::CalcBulkS(double T, 
                                const double *const x, 
                                unsigned int n) const
{
    return CalcBulkS(T, x, n, m_tmpvec);
}


// DIMENSIONLESS ENTROPY.

double ThermoInterface::CalcBulkS_R(double T, const Sprog::fvector &x) const
{
    return CalcBulkS_R(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkS_R(double T, 
                                  const Sprog::fvector &x, 
                                  Sprog::fvector &S) const
{
    return CalcBulkS_R(T, &x[0], x.size(), S);
}

double ThermoInterface::CalcBulkS_R(double T, 
                                  const double *const x, 
                                  unsigned int n) const
{
    return CalcBulkS_R(T, x, n, m_tmpvec);
}


// GIBBS FREE ENERGY.

double ThermoInterface::CalcBulkG(double T, const Sprog::fvector &x) const
{
    return CalcBulkG(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkG(double T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &G) const
{
    return CalcBulkG(T, &x[0], x.size(), G);
}

double ThermoInterface::CalcBulkG(double T, 
                                const double *const x, 
                                unsigned int n) const
{
    return CalcBulkG(T, x, n, m_tmpvec);
}


// DIMENSIONLESS GIBBS FREE ENERGY.

double ThermoInterface::CalcBulkG_RT(double T, const Sprog::fvector &x) const
{
    return CalcBulkG_RT(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkG_RT(double T, 
                                   const Sprog::fvector &x, 
                                   Sprog::fvector &G) const
{
    return CalcBulkG_RT(T, &x[0], x.size(), G);
}

double ThermoInterface::CalcBulkG_RT(double T, 
                                   const double *const x, 
                                   unsigned int n) const
{
    return CalcBulkG_RT(T, x, n, m_tmpvec);
}


// CONSTANT PRESSURE HEAT CAPACITY.

double ThermoInterface::CalcBulkCp(double T, const Sprog::fvector &x) const
{
    return CalcBulkCp(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkCp(double T, 
                                 const Sprog::fvector &x, 
                                 Sprog::fvector &Cp) const
{
    return CalcBulkCp(T, &x[0], x.size(), Cp);
}

double ThermoInterface::CalcBulkCp(double T, 
                                 const double *const x, 
                                 unsigned int n) const
{
    return CalcBulkCp(T, x, n, m_tmpvec);
}


// DIMENSIONLESS CONSTANT PRESSURE HEAT CAPACITY.

double ThermoInterface::CalcBulkCp_R(double T, const Sprog::fvector &x) const
{
    return CalcBulkCp_R(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkCp_R(double T, 
                                   const Sprog::fvector &x, 
                                   Sprog::fvector &Cp) const
{
    return CalcBulkCp_R(T, &x[0], x.size(), Cp);
}

double ThermoInterface::CalcBulkCp_R(double T, 
                                   const double *const x, 
                                   unsigned int n) const
{
    return CalcBulkCp_R(T, x, n, m_tmpvec);
}


// CONSTANT VOLUME HEAT CAPACITY.

double ThermoInterface::CalcBulkCv(double T, const Sprog::fvector &x) const
{
    return CalcBulkCv(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkCv(double T, 
                                 const Sprog::fvector &x, 
                                 Sprog::fvector &Cv) const
{
    return CalcBulkCv(T, &x[0], x.size(), Cv);
}

double ThermoInterface::CalcBulkCv(double T, 
                                 const double *const x, 
                                 unsigned int n) const
{
    return CalcBulkCv(T, x, n, m_tmpvec);
}


// DIMENSIONLESS CONSTANT VOLUME HEAT CAPACITY.

double ThermoInterface::CalcBulkCv_R(double T, const Sprog::fvector &x) const
{
    return CalcBulkCv_R(T, &x[0], x.size());
}

double ThermoInterface::CalcBulkCv_R(double T, 
                                   const Sprog::fvector &x, 
                                   Sprog::fvector &Cv) const
{
    return CalcBulkCv_R(T, &x[0], x.size(), Cv);
}

double ThermoInterface::CalcBulkCv_R(double T, 
                                   const double *const x, 
                                   unsigned int n) const
{
    return CalcBulkCv_R(T, x, n, m_tmpvec);
}
