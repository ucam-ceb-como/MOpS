/* 
 * File:   cam_reporter.h
 * Author: vinod
 *
 * Copyright (C) 2009 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class implements the file output and the screen output
 * for all the reactor models
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
 * Created on January 24, 2009, 7:54 PM
 */

#ifndef _CAM_REPORTER_H
#define	_CAM_REPORTER_H

#include "cam_boundary.h"
#include "cam_residual.h"
#include "gpc.h"
#include "comostrings.h"
#include <fstream>

namespace Camflow{

    class CamReporter{

        DataIO *standard;
        DataIO *rates;
        DataIO *transport;
        DataIO *custom;

    public:

        CamReporter();

        ~CamReporter();

        void header(std::string prog);
        void problemDescription(CamBoundary &cb, CamResidual &cr);
        //title for the consol output
        void consoleHead(std::string head);
        //write the data to file
        void openFile(std::string fileName, bool old);
        void closeFile();
        void openFiles(bool stdrd = true, bool ratesOut = false, bool transOut = false);
        void closeFiles(bool stdrd = true, bool ratesOut = false, bool transOut = false);
        void writeHeader(std::vector<std::string>& stdHeader);
        void writeHeader(std::vector<std::string>& stdHeader, std::vector<std::string>& ratesHeader);
        void wrteHeader(std::vector<std::string>& stdHeader, std::vector<std::string>& ratesHeader, std::vector<std::string>& transHeader);
        void writeToFile(std::string fname,CamResidual &resid);
        void writeStdFileOut(std::vector<doublereal>& data);
        void writeCustomHeader(std::vector<std::string>& header);
        void writeCustomFileOut(std::vector<doublereal>& data);
    };
}


#endif	/* _CAM_REPORTER_H */

