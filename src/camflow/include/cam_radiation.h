/*!
 * \file   cam_radiation.h
 * \author L. R. McGlashan
 *
 * \brief Radiation.
 *
 *  Copyright (C) 2011 L. R. McGlashan.
 *

 Licence:
    This file is part of "camflow".

    brush is free software; you can redistribute it and/or
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
    Prof Markus Kraft
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

#ifndef _RADIATION_H
#define	_RADIATION_H

#include "cam_params.h"
#include "gpc_mech.h"
#include "array.h"

#include <map>

namespace Camflow
{

/*! @brief Radiation equation solver class.
 *
 *  This code computes the spectral and soot related components or radiative heat loss.
 *
 *   References:
 *   1.  Radiation Models, International Workshop on Measurement and Computation of Turbulent  Nonpremixed Flames,
 *       www.sandia.gov/TNF/radiation.html.  This reference covers the details of implementing the spectral part of
 *       the radiation term.
 *
 *   2.  Kim, S-K., Kim, Y., Assessment of the Eulerian particle flamelet model for nonpremixed turbulent jet flames,
 *       Combustion and Flame 154 (2008) 232-247.  This article shows a modern implementation of the procedure used
 *       in [1], for the spectral part of the radiation term.
 *
 *   3.  Carbonnel, D., Oliva, A., Perez-Segarra, C.D., Implementation of two-equation soot flamelet models for laminar
 *       diffusion flames.  This paper includes a term, used here, for modelling soot radiative heat loss.
 *
 *   4.  Grosshandler, W.L., RADCAL: A Narrow-Band Model For Radiation Calculations in a Combustion Environment, NIST
 *       technical note 1402, 1993.
*/
class Radiation
{

    //! Conversion factor.
    static const doublereal AtmToPascal;

    //! Names of the radiative species.
    std::vector<std::string> radiativeSpecies_;
    //! Indices of the species for looking up in Sprog::Mechanism.
    std::vector<doublereal> speciesIndex_;
    //! Molecular weights of each species.
    std::vector<doublereal> speciesMolWt_;

    //! Absorption coefficients of each species.
    std::vector<doublereal> absorption;
    //! Partial pressures for each species.
    std::vector<doublereal> partialPress;

    //! Stores the radiation sources.
    std::vector<doublereal> radiation;

    const Sprog::Mechanism *const mech_;
    const std::vector<doublereal>& avgMolWt_;
    const Array2D& speciesMassFracs_;

    //! Computes the Planck mean absorption constants,
    //! as input to the radiative heat loss dissipation model.
    void PlanckAbsorption(const doublereal Temperature);

public:

    //! Default constructor.
    Radiation
    (
        const std::string& inputFileName,
        const int totalCells,
        const Sprog::Mechanism *const mech,
        const std::vector<doublereal>& avgMolWt,
        const Array2D& s_mf
    );

    //! Destructor.
    ~Radiation();

    //! Computes the radiative heat loss term
    //! for the radiative heat dissipation model.
    void calculateRadiativeHeatLoss
    (
        const int i,
        const doublereal& Temperature,
        const doublereal& opPre,
        const doublereal& soot_vol_frac
    );

    inline const doublereal& getRadiation(const int i)
    {
        return radiation[i];
    }

}; // End Radiation class declaration.

} // End Camflow namespace.

#endif	/* _RADIATION_H */
