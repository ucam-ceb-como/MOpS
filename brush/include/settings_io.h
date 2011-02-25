/*!
 * \file   settings_io.h
 * \author Robert I A Patterson
 *
 * \brief  Methods for parsing XML settings fragments
 *
 *  Copyright (C) 2009 Robert I A Patterson.
 *

 Licence:
    This file is part of "brush".

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
#ifndef BRUSH_SETTINGS_IO_H
#define BRUSH_SETTINGS_IO_H

#include <vector>
#include <string>

#include "brush_params.h"

#include "camxml.h"

namespace Brush {
namespace Settings_IO {

    //! Read a 1d spatial profile of elements with specified name
    std::vector<std::pair<real, real> > readProfile(const CamXML::Element * const document,
                                                    const std::string &name);

} // namespace Settings_IO
} // namespace Brush


#endif // BRUSH_SETTINGS_IO_H

