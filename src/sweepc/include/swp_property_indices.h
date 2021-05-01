/*!
 * \file   swp_property_indices.h
 * \author Robert I A Patterson
 *  Copyright (C) 2010 Robert I A Patterson.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 *
 * \brief  Symbolic indices for the properties of particles
 *
 Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

#ifndef SWEEP_PROPERTY_INDICES_H
#define SWEEP_PROPERTY_INDICES_H


namespace Sweep
{
    //! Symbolic indices for particle properties.
    enum PropID {
        iUniform=-1, /**< Special Case:  Always returns 1.0.  Used to select particles uniformly. */
        iDsph,       /**< Equivalent sphere diameter. */
        iDcol,       /**< Collision diameter. */
        iDmob,       /**< Mobility diameter. */
        iS,          /**< Surface area. */
        iV,          /**< Volume. */
        iM,          /**< Mass. */
		iNumCarbon,  /**< Number of carbon atoms. */
        iFrag,       /**< Fragmentation flag. */

		// The next properties are provided for calculation of
		// collision rates.
		iD2,      // Collision diameter squared.
		iD_1,     // Inverse collision diameter.
		iD_2,     // Inverse of the diameter squared.
		iM_1_2,   // Inverse of the square-root of the mass.

		//! Statistical weight
		iW,

		//! Statistical weight time physical mass
		iWM,
		iDW,		// dcol * weight
		iD2W,		// dcol * dcol * weight
		iD_1W,		// weight / dcol
		iD_2W,		// weight / dcol ^ 2
		iM_1_2W,	// mass ^ -1/2 * weight
		iD2_M_1_2W, // dcol * dcol * mass ^ -1/2 * weight

		iD2_M_1_2, // D^2 * M^-1/2.
		iFS,		// the free surface available for other particles to sinter

		// Silica model properties
		iASN, // Number of active (OH) sites available
		iSintRate, // Sintering rate of a particle

		// Silicon properties
		iCoverage,     // Ratio of component 0 to component 1

		iUniform1,


		//Titania model properties
		iAn_2_3_comp	//Anatase fraction ^ (2/3) * total composition
    };
}

#endif
