/*!
  \Author       Rebecca Riehl (rcr33)
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

#include "loi_reduction.h"
using namespace std;
using namespace Mops;

/*!
@param[in]        J           a real matrix
@param[in]        Sensi       a real N_Vector
@param[in]        LOI         a double vector array, stores LOI
@param[in]        n_species   number of species in the reaction mechanism
@param[in]        n_sensi     number of sensitivities calculated in the reaction mechanism
@return           LOI         fvector with LOI results
*/
vector<fvector> Mops::LOIReduction::CalcLOI(double** J, double** Sensi, vector<fvector> LOI, int n_species, int n_sensi)
{ 
    for (int i = 0; i < n_sensi; i++){
        for (int j = 0; j < n_species; j++){ 
            if (Sensi[i][j] < 0){
                Sensi[i][j] = -Sensi[i][j];
            }
            LOI[i][j] = LOI[i][j] - Sensi[i][j] / J[j][j]; 
            cout << LOI[i][j] << "\n";
       }
    }
    return LOI;
}

