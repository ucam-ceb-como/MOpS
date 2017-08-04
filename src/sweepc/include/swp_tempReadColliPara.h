/*
  Author(s):      Dongping Chen (dc516)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2011 Dongping Chen.

  File purpose:
    temporarily used to read collision efficiency parameters, including mode, NONE, MAX, MIN, COMBINED

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

#ifndef SWEEP_TEMPREADCOLLIPARA_H
#define SWEEP_TEMPREADCOLLIPARA_H

#include "swp_mechanism.h"
#include "swp_maths_functional.h"

#include "camxml.h"
#include "string_functions.h"
#include "csv_io.h"
#include <stdexcept>
#include <string>
#include <cassert>
#include <memory>
#include <cstdlib>


namespace Sweep
{
void readColliPara(std::vector<CamXML::Element*>::iterator i, Sweep::Mechanism &mech)
{
    CamXML::Element *el;
    CamXML::Element *el0;
    CamXML::Element *el1;
    CamXML::Element *el2;

    el = (*i)->GetFirstChild("mode");
    if (el!=NULL) {
        std::string str = el->Data();
        if (str == "NONE") {
                //! Get component's parameter for collision efficiency model.
            el0 = (*i)->GetFirstChild("ColliParaA");
            el1 = (*i)->GetFirstChild("ColliParaB");
            el2 = (*i)->GetFirstChild("ColliParaC");
            if (el0!=NULL && el1!=NULL && el2!=NULL) {
                std::string str0 = el0->Data();
                std::string str1 = el1->Data();
                std::string str2 = el2->Data();
                if (str0 != "" && str1 != "" && str2 != "") { 
                    mech.SetCollisionEffPara(Strings::cdble(str0), Strings::cdble(str1), Strings::cdble(str2));
                } else {
                    //! Collision parameter contains no data.
                    std::string msg("one of collision parameter is left to be empty (Sweep::readColliPara).");

                    throw std::runtime_error(msg);
                }
            } else  {
                std::string msg("one of collision parameter (ColliParaA, ColliParaB,ColliParaC) is missing (Sweep::readColliPara).");

                throw std::runtime_error(msg);
            }
        }

        else if (str == "MAX" || str == "MIN" || str == "COMBINED" || str == "REDUCED") {
            mech.SetMode(str);
            el0 = (*i)->GetFirstChild("inceptionThreshold");
            if (el0!=NULL) {
            str = el0->Data();
                if (str != "") {
                    mech.setInceptionThreshold(Strings::cdble(str)); //!< Number of 6-member rings.
                } else {
                    //! inceptionThreshold contains no data.
                    std::string msg("threshold contains no data (Sweep::readColliPara).");

                    throw std::runtime_error(msg);
                }
            } else {
            std::string msg("inceptionThreshold is required for collision efficiency model (Sweep::readColliPara).");

            throw std::runtime_error(msg);
            }

            el0 = (*i)->GetFirstChild("condensationThreshold");
            if (el0!=NULL) {
            str = el0->Data();
                if (str != "") {
                    mech.setCondensationThreshold(Strings::cdble(str)); //!< Number of 6-member rings.
                } else {
                    //! condensationThreshold contains no data.
                    std::string msg("condensationThreshold contains no data (Sweep::readColliPara).");

                    throw std::runtime_error(msg);
                }
            } else {
            std::string msg("condensationThreshold is required for collision efficiency model (Sweep::readColliPara).");

            throw std::runtime_error(msg);
            }
        }
        else {
            std::string msg(" This type of collision efficiency model is not supported ");
            msg += str;
            msg += " at this moment (Sweep::readColliPara).";

            throw std::runtime_error(msg);
        }
    }
    else  mech.SetCollisionEffPara(2, 1100, 5); //!< Set to the default one abhjeet's choice, this is use to pass the regression tests. 
}
};

#endif
