/*!
 * \file   cam_params.h
 * \author V. Janardhanan
 *
 * \brief Definition of Camflow real type variables.
 *
 *  Copyright (C) 2009 Vinod Janardhanan.
 *
 * Licence:
 *  This file is part of "Camflow".
 *
 *  flameLab is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * Contact:
 *  Dr Markus Kraft
 *  Dept of Chemical Engineering
 *  University of Cambridge
 *  New Museum Site
 *  Pembroke Street
 *  Cambridge
 *  CB2 3RA
 *  UK
 *
 *  Email   :   mk306@cam.ac.uk
 *  Website :   http://como.cheng.cam.ac.uk
 *
 * Created on January 17, 2009, 1:10 PM
 */

#ifndef _CAM_PARAMS_H
#define	_CAM_PARAMS_H

#include  "sprog.h"

/*!
 *@brief    Camflow namespace.
 *
 * Include a more detailed description here.
 */
namespace Camflow {

    //! Typedef for doublereal.
    typedef Sprog::real doublereal;

    //! cgs units used in the soot model.
    const doublereal kB_cgs = 1.3806504e-16;
    
}

#endif	/* _CAM_PARAMS_H */
