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

#include "loi_reduction.h"
#include "gpc_species.h"
using namespace std;
using namespace Mops;

/*!
@param[in]        J           a double matrix
@param[in]        Sensi       a double N_Vector
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
            if (J[j][j] == 0.0){
                LOI[i][j] = LOI[i][j] + 0.0;
            }
            else{
                LOI[i][j] = LOI[i][j] - Sensi[i][j] / J[j][j]; 
            }
       }
    }
    return LOI;
}


/*!
@param[in]      LOI             fvector with LOI results
@param[in]      LOICompVal      double against which to compare an LOI value
@param[in]      mech            Reaction mechanism as defined in Sprog.
@param[in]      Kept_Spec       User-defined species that must be kept in the reduced mechanism.
@param[in, out] RejectSpecies   A string vector to store rejected species' names
*/
void Mops::LOIReduction::RejectSpecies(vector<fvector> LOI, double LOICompVal, const Mechanism *const mech,
                                       std::vector<std::string>& RejectSpecies, std::vector<std::string> Kept_Spec)
{
    for (unsigned int i = 0; i < mech->GasMech().SpeciesCount(); i++){
        if (LOI[0][i] < LOICompVal){
            for (unsigned int j = 0; j < Kept_Spec.size(); j++){
                if (mech->GasMech().Species(i)->Name() == Kept_Spec[j]){
                    RejectSpecies.push_back("NULL");
                break;
                }
                else if (j == (Kept_Spec.size()-1) && Kept_Spec[j] != mech->GasMech().Species(i)->Name())
                    RejectSpecies.push_back(mech->GasMech().Species(i)->Name());
            }
        }
        else {
            RejectSpecies.push_back("NULL");
        }
    }
}

/*!
@param[in]      LOIFile     Ofstream used to write data to LOI csv output file
@param[in]      mech        Reaction mechanism as defined in Sprog.
*/
void Mops::LOIReduction::CreateLOIFile(std::ofstream &LOIFile, const Mechanism *const mech)
{
    typedef std::vector<std::string> string_vector;

    string_vector LOIHeadings; 
    LOIHeadings.push_back("Time/ s");

    for (unsigned int isp=0; isp<mech->GasMech().SpeciesCount(); ++isp) {
    LOIHeadings.push_back(mech->GasMech().Species(isp)->Name() + " (LOI)");
    }
    
    
    // Write the column headings to the file, with commas between names, but not after the last name
    string_vector::const_iterator it = LOIHeadings.begin();
    const string_vector::const_iterator itStop = LOIHeadings.end() - 1;
    while(it != itStop) {
        LOIFile << *it++ << ',';
    }
    // Final heading
    LOIFile << LOIHeadings.back() << '\n';

}

/*!
return      name        Returns the name of the LOI file
*/
std::string Mops::LOIReduction::buildLOIFileName(const std::string& n)
{
    std::stringstream name;
    name << n << std::string("-loi.csv");
    return name.str();

}

/*!
@param[in]      LOI     fvector of LOI values
@param[in]      t       The current time.
@param[in, out] out     File that is being written to.
@param[in]      mech    Reaction mechanism as defined in Sprog.
*/
void Mops::LOIReduction::SaveLOI(std::vector<fvector> LOI, double t, std::ofstream &out, const Mechanism *const mech)
{
    //Write out the time point
    out << t;

    //Write out the LOI for each species
    for (unsigned int i = 0; i < mech->GasMech().SpeciesCount(); ++i){
        out << ',' << LOI[0][i];
    }

    //Go to a new line
    out << '\n';

    // Flush the buffer, so that a future crash does not lose data before it actually
    // reaches its destination.
    out.flush();

}

