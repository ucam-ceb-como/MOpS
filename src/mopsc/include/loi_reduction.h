/*!
  \author       Rebecca Riehl (rcr33)
  \file         loi_reduction.h
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2010 Rebecca Riehl.

  \brief Namespace for mechanism reduction using the Level of Importance (LOI) technique.
  File purpose:
    Level of Importance (LOI) mechanism reduction namespace.
    This method uses both the Jacobian and the sensitivity matrix for
    a reaction to assign a LOI parameter to each species. Species with
    a LOI beneath user tolerance specifications are deemed redundant.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
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

#ifndef LOI_REDUCTION_H
#define LOI_REDUCTION_H

#include <vector>

#include "mops_ode_solver.h"
#include "mops_simulator.h"
#include "mops_gpc_sensitivity.h"


namespace Mops
{
//! Calculates level of importance reaction mechanism reduction
namespace LOIReduction
{

    //!Calculates a Level of Importance for each species in a mechanism by multiplying the Jacobian and sensitivity matrices.
    std::vector<fvector> CalcLOI(double** J, double** Sensi, std::vector<fvector> LOI, int n_species, int n_sensi);

    //! Compares LOI value to a given cutoff then makes lists of kept or rejected species
    void RejectSpecies(std::vector<fvector> LOI, double LOICompVal, const Mechanism *const mech,
        std::vector<std::string>& RejectSpecies, std::vector<std::string> Kept_Spec);

    //! Creates an output file for each species' LOI at each timepoint.
    void CreateLOIFile(std::ofstream &LOIFile, const Mechanism *const mech);

    //! Build the name for the LOI file.
    std::string buildLOIFileName(const std::string &n);

    //! Save the data from one run in the file created by the function above.
    void SaveLOI(std::vector<fvector> LOI, double t, std::ofstream &out, const Mechanism *const mech);


} //namespace LOIReduction
} //namespace Mops

#endif
