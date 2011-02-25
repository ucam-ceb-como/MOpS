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

real ThermoInterface::CalcBulkU(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkU(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkU(Sprog::real T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &U) const
{
    return CalcBulkU(T, &x[0], x.size(), U);
}

real ThermoInterface::CalcBulkU(Sprog::real T, 
                                const Sprog::real *const x, 
                                unsigned int n) const
{
    return CalcBulkU(T, x, n, m_tmpvec);
}


// DIMENSIONLESS INTERNAL ENERGY.

real ThermoInterface::CalcBulkU_RT(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkU_RT(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkU_RT(Sprog::real T, 
                                   const Sprog::fvector &x, 
                                   Sprog::fvector &U) const
{
    return CalcBulkU_RT(T, &x[0], x.size(), U);
}

real ThermoInterface::CalcBulkU_RT(Sprog::real T, 
                                   const Sprog::real *const x, 
                                   unsigned int n) const
{
    return CalcBulkU_RT(T, x, n, m_tmpvec);
}


// ENTHALPY.

real ThermoInterface::CalcBulkH(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkH(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkH(Sprog::real T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &H) const
{
    return CalcBulkH(T, &x[0], x.size(), H);
}

real ThermoInterface::CalcBulkH(Sprog::real T, 
                                const Sprog::real *const x, 
                                unsigned int n) const
{
    return CalcBulkH(T, x, n, m_tmpvec);
}


// DIMENSIONLESS ENTHALPY.

real ThermoInterface::CalcBulkH_RT(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkH_RT(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkH_RT(Sprog::real T, const Sprog::fvector &x, 
                                   Sprog::fvector &H) const
{
    return CalcBulkH_RT(T, &x[0], x.size(), H);
}

real ThermoInterface::CalcBulkH_RT(Sprog::real T, 
                                   const Sprog::real *const x, 
                                   unsigned int n) const
{
    return CalcBulkH_RT(T, x, n, m_tmpvec);
}


// ENTROPY.

real ThermoInterface::CalcBulkS(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkS(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkS(Sprog::real T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &S) const
{
    return CalcBulkS(T, &x[0], x.size(), S);
}

real ThermoInterface::CalcBulkS(Sprog::real T, 
                                const Sprog::real *const x, 
                                unsigned int n) const
{
    return CalcBulkS(T, x, n, m_tmpvec);
}


// DIMENSIONLESS ENTROPY.

real ThermoInterface::CalcBulkS_R(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkS_R(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkS_R(Sprog::real T, 
                                  const Sprog::fvector &x, 
                                  Sprog::fvector &S) const
{
    return CalcBulkS_R(T, &x[0], x.size(), S);
}

real ThermoInterface::CalcBulkS_R(Sprog::real T, 
                                  const Sprog::real *const x, 
                                  unsigned int n) const
{
    return CalcBulkS_R(T, x, n, m_tmpvec);
}


// GIBBS FREE ENERGY.

real ThermoInterface::CalcBulkG(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkG(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkG(Sprog::real T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &G) const
{
    return CalcBulkG(T, &x[0], x.size(), G);
}

real ThermoInterface::CalcBulkG(Sprog::real T, 
                                const Sprog::real *const x, 
                                unsigned int n) const
{
    return CalcBulkG(T, x, n, m_tmpvec);
}


// DIMENSIONLESS GIBBS FREE ENERGY.

real ThermoInterface::CalcBulkG_RT(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkG_RT(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkG_RT(Sprog::real T, 
                                   const Sprog::fvector &x, 
                                   Sprog::fvector &G) const
{
    return CalcBulkG_RT(T, &x[0], x.size(), G);
}

real ThermoInterface::CalcBulkG_RT(Sprog::real T, 
                                   const Sprog::real *const x, 
                                   unsigned int n) const
{
    return CalcBulkG_RT(T, x, n, m_tmpvec);
}


// CONSTANT PRESSURE HEAT CAPACITY.

real ThermoInterface::CalcBulkCp(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkCp(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkCp(Sprog::real T, 
                                 const Sprog::fvector &x, 
                                 Sprog::fvector &Cp) const
{
    return CalcBulkCp(T, &x[0], x.size(), Cp);
}

real ThermoInterface::CalcBulkCp(Sprog::real T, 
                                 const Sprog::real *const x, 
                                 unsigned int n) const
{
    return CalcBulkCp(T, x, n, m_tmpvec);
}


// DIMENSIONLESS CONSTANT PRESSURE HEAT CAPACITY.

real ThermoInterface::CalcBulkCp_R(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkCp_R(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkCp_R(Sprog::real T, 
                                   const Sprog::fvector &x, 
                                   Sprog::fvector &Cp) const
{
    return CalcBulkCp_R(T, &x[0], x.size(), Cp);
}

real ThermoInterface::CalcBulkCp_R(Sprog::real T, 
                                   const Sprog::real *const x, 
                                   unsigned int n) const
{
    return CalcBulkCp_R(T, x, n, m_tmpvec);
}


// CONSTANT VOLUME HEAT CAPACITY.

real ThermoInterface::CalcBulkCv(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkCv(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkCv(Sprog::real T, 
                                 const Sprog::fvector &x, 
                                 Sprog::fvector &Cv) const
{
    return CalcBulkCv(T, &x[0], x.size(), Cv);
}

real ThermoInterface::CalcBulkCv(Sprog::real T, 
                                 const Sprog::real *const x, 
                                 unsigned int n) const
{
    return CalcBulkCv(T, x, n, m_tmpvec);
}


// DIMENSIONLESS CONSTANT VOLUME HEAT CAPACITY.

real ThermoInterface::CalcBulkCv_R(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkCv_R(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkCv_R(Sprog::real T, 
                                   const Sprog::fvector &x, 
                                   Sprog::fvector &Cv) const
{
    return CalcBulkCv_R(T, &x[0], x.size(), Cv);
}

real ThermoInterface::CalcBulkCv_R(Sprog::real T, 
                                   const Sprog::real *const x, 
                                   unsigned int n) const
{
    return CalcBulkCv_R(T, x, n, m_tmpvec);
}
