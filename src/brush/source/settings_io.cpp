/*!
 * \file   settings_io.cpp
 * \author Robert I A Patterson
 *
 * \brief Methods for parsing XML settings fragments
 *
 *  Copyright (C) 2009 Robert I A Patterson.

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

#include "settings_io.h"

#include <cstdlib>
#include <stdexcept>

/*!
 * Read in a 1d spatial profile of values of `blah' from an XML tree containing
 * one or more top level children of the form <blah x=1.23>5.2e7</blah>.
 * Alternatively the XML tree may contain exactly one global value in a top
 * level child of the form <blah>5.2e7</blah>.  Numbers can be specified in any
 * text format that can be converted by atof.
 *
 *@param[in]    document    root node of XML document containing profile
 *@param[in]    name        Name of XML element for which a profile is to be built
 *
 *@return       Profile of element values and positions as a vector of pairs
 *
 *@exception    std::runtime_error  No data
 *@exception    std::runtime_error  Mixture of global and local data
 */
std::vector<std::pair<double, double> >
  Brush::Settings_IO::readProfile(const CamXML::Element * const document,
                                  const std::string &name) {
    // Use to build up the return value
    std::vector<std::pair<double, double> > profile;

    // Extract the relevant XML elements
    std::vector<CamXML::Element*> profileElements;
    document->GetChildren(name, profileElements);

    // Check some data has been provided
    if (profileElements.empty()) {
        std::string msg("At least one <");
        msg += name;
        msg += "> element must be supplied to specify a profile";
        throw std::runtime_error(msg);
    }
    else {
        // Read the XML one element at a time
        for(std::vector<CamXML::Element*>::const_iterator it = profileElements.begin();
            it != profileElements.end(); ++it) {
            // Loop over the pcount entries
            const double value = std::atof((*it)->Data().c_str());

            // Look for position information
            const CamXML::Attribute* const  posnAttrib = (*it)->GetAttribute("x");
            double x = 0;
            if(posnAttrib) {
                //Found a position so store the full setting
                x = std::atof(posnAttrib->GetValue().c_str());
            }
            else if(profileElements.size() == 1) {
                // Allow one element to be specified without a position, if no other
                // elements are specified so a global setting may be given
                x = 0.0;
            }
            else {
                std::string msg("x (position) attribute must be supplied for all <");
                msg += name;
                msg += "> elements\n";
                msg += "or exactly one global <";
                msg += name;
                msg += "> must be supplied without an x attribute";
                throw std::runtime_error(msg);
            }
            profile.push_back(std::make_pair<double, double>(x, value));
        }
    }
    return profile;
}
