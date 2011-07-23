/*!
 * File:   cam_read.h
 * Author: vinod (vj231@cam.ac.uk)
 *
 * Copyright (C) 2008 Vinod M Janardhanan
 *
 * File purpose:
 *  This class contains the implementation IO functions
 *
 * Licence:
 *  This file is part of "Camflow".
 *
 *  Camflow is free software; you can redistribute it and/or
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
 * Created on January 17, 2009, 5:38 PM
 */

#ifndef _CAM_READ_H
#define	_CAM_READ_H

#include <string>
#include <vector>

#include "cam_geometry.h"
#include "cam_admin.h"
#include "cam_configuration.h"
#include "cam_profile.h"
#include "camxml.h"
#include "cam_converter.h"
#include "cam_soot.h"
#include "cam_control.h"
#include "cam_boundary.h"
#include "string_functions.h"
#include "cam_error.h"

namespace Camflow{

	/*!
	 *@brief    This class contains the implementation of IO functions.
	 *
	 * Include a more detailed description here.
	 */

    class CamRead{

    public:
        CamRead(){}
        ~CamRead(){}

        //main call to read input data set
        void readInput(const std::string fileName,
                       CamControl& cc,
                       CamGeometry& cg,
                       CamConverter& convert,
                       CamAdmin& ca,
                       CamBoundary& cb,
                       CamProfile& cp,
                       CamConfiguration& config, CamSoot &cSoot);
        //read the geometry information
        void readGeometry(CamGeometry &cg,
                          CamConfiguration& config,
                          CamConverter& convert,
                          const CamXML::Element &node);
        //read process conditions
        void readProcessConditions(CamConverter& convert,
                                   CamAdmin& ca,
                                   const CamXML::Element &node);
        //read boundary condition
        void readBoundary(CamAdmin& ca,
                            CamBoundary& cb,
                            CamConverter& convert,
                            const CamXML::Element &node);

        //common function for inlet reading
        void readNozzle(CamBoundary& cb, CamConverter convert,
                        const CamXML::Element &node);

        //read control parameters
        void readControl(CamControl& cc, const CamXML::Element &node);

        //read initial guess
        void readInitialGuess(
                CamAdmin& ca,
                CamProfile& cp,
                CamConverter& convert,
                const CamXML::Element &node);

        //read the report mode
        void readReport(CamAdmin& ca, const CamXML::Element &node);

        //read the grid file
        void readGrid(CamGeometry& cg, const CamXML::Element &node);

    private:

        //read tolerences
        void readTol(const CamXML::Element &node, doublereal& atol, doublereal& rtol);

        //read mass or mole fractions
        void readFrac(std::string& member, std::map<std::string,doublereal>& fracs, const CamXML::Element &subnode);


    };
}

#endif	/* _CAM_READ_H */

