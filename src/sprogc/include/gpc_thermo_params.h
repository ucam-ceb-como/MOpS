/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This files contains the definition of a structure for chemical species thermodynamic
    fitting parameters.  Thermo parameters are split into ranges defined by the temperatures
    for which the ranges are valid.

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

#ifndef GPC_THERMO_PARAMS_H
#define GPC_THERMO_PARAMS_H

#include <map>
#include "gpc_params.h"

namespace Sprog
{
namespace Thermo
{
    // Maximum lengths of thermo parameter arrays (S_PARAM_COUNT).
    const unsigned int CP_PARAM_COUNT    = 5;
    const unsigned int H_PARAM_COUNT     = 6;
    const unsigned int S_PARAM_COUNT     = 7;
    const unsigned int MAX_THERMO_PARAMS = 10;

    struct THERMO_PARAMS
    {
        unsigned int Count;
        double Params[S_PARAM_COUNT];

        // Writes the element to a binary data stream.
        template<class Archive>
        void serialize(Archive & ar, const unsigned int /* file_version */)
        {
            ar & Count & Params;
        }
        friend class boost::serialization::access;

    };

    typedef std::map<double, THERMO_PARAMS> ThermoMap;
};
};
#endif
