/*!
  * \author     Zakwan Zainuddin (zz260)
  * \file       swp_kmc_structure_comp.h
  *
  * \brief        Defines the component struct and classes for the data structure
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Defines the component struct and classes for the data structure

  Licence:
    This file is part of "sweep".

    Sweep is free software; you can redistribute it and/or
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

#ifndef SWP_KMC_STRUCTURE_COMP_H
#define SWP_KMC_STRUCTURE_COMP_H

#include "swp_kmc_typedef.h"
#include "boost/shared_ptr.hpp"
#include <list>
#include <cmath>
#include <cstdlib>
#include <set>

namespace Sweep {
    namespace KMC_ARS {
        // Forward Declarations
        struct Site;
        class Carbon;

        // Pointer to a Carbon atom and a Site.
        typedef std::set< Carbon* > Ccontainer;
        typedef Carbon* Cpointer;
        typedef std::list< Site >::iterator Spointer;
        

        //! Carbon atom on PAH edges
        class Carbon {
        public:
            Carbon();

            //! Default Destructor
            ~Carbon();
            //! Neighbouring carbon atoms
            Cpointer C1; // prev
            Cpointer C2; // next
            //! 3rd carbon attached
            Cpointer C3;
            //! Bulk or edge C atom
            bool edge;
             //! For bridging carbon atoms
            bool bridge;
            //! Third species bonded to it, 0 if none [radical]
            char A;
            //! Bond angle made with C2
            angletype bondAngle1;
            //! Bond angle at bridge
            angletype bondAngle2;
            //! Coordinates of atom
            cpair coords;
        };

        //static Ccontainer NULLSET(1, Carbon());
        static Cpointer NULLC(new Carbon);//SET.begin();

        //! Active site on PAH edges
        struct Site {
            //! Type of site
            kmcSiteType type;
            //! Type of combined site
            kmcSiteType comb;
            //! Reactive carbon members
            Cpointer C1;
            Cpointer C2;
        };
        

        // To calculate increments in x coordinates and y coordinates given angle direction
        coordtype inline x_inc(const angletype &angle) {
            //coordtype scale = 2.0f;
            if(abs(angle)==0) return (2);
            else if(std::abs(angle)==60) return (1); // for 6R
            else if(std::abs(angle)==120) return (-1);
            else if(std::abs(angle)==180) return (-2);
            else if(std::abs(angle)==30) return (3); // for 5R
            else if(std::abs(angle)==90) return (0);
            else if(std::abs(angle)==150) return(-3);
            return -10;
        }
        coordtype inline y_inc(const angletype &angle) {
            //coordtype scale = 2.0f;
            if(angle==0 || std::abs(angle)==180) return 0; // for 6R
            else if(angle==60 || angle==120) return (2);
            else if(angle==-60 || angle==-120) return (-2);
            else if(angle==90) return (4); // for 5R
            else if(angle==-90) return (-4);
            else if(angle==30 || angle==150) return (2);
            else if(angle==-30 || angle==-150) return (-2);
            return -10;
        }
        //! Changes reflex angle into +ve and -ve angle format
        angletype inline normAngle(angletype angle) {
            angletype temp = angle;
            while(std::abs(temp)>180) {
                if (temp > 180) {
                    temp-= 360;
                }else if (temp < -180) {
                    temp+= 360;
                }
            }
            return temp;
        }
        
    }
}

#endif
