/*
 * File:   cam_reporter.h
 * Author: vinod
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

#include "cam_reporter.h"

using namespace Camflow;
using namespace std;


CamReporter::CamReporter()
:
    standard(NULL),
    rates(NULL),
    transport(NULL),
    custom(NULL)
{}


CamReporter::~CamReporter()
{
    if (standard != NULL) delete standard;
    if (rates != NULL) delete rates;
    if (transport != NULL) delete transport;
    if (custom != NULL) delete custom;
}


void CamReporter::header(string prog)
{
    cout << endl;
    cout << "!----------------------------------------------------------!\n";
    cout << "!               CamFlow " << prog << " Version 1.0               !\n";
    cout << "!        This program is distributed without any warranty  !\n";
    cout << "!            Author: Vinod M. J. (vj231@cam.ac.uk)         !\n";
    cout << "!----------------------------------------------------------!\n";
}


void CamReporter::problemDescription(CamBoundary& cb, CamResidual& cr)
{
    //screen out put for simulation case
    cout << "\t Problem description.....\n";
    cout << "\t Number of species: " << cr.getNSpecies() << endl;
    cout << "\t Total number of equations solving: " << cr.getNEqn() << endl;
    cout << "\t Total number of variables solving: " << cr.getNVar() << endl;

    //initializing solution vector
    cout << "\t Initializing solution vector...\n";
    cout << "\t Inlet velocity: "<<cb.getVelocity() <<  " m/s" << endl;
    cout << "\t Inlet temperature "<< cb.getTemperature() << " K"<< endl;
    map<string, doublereal> species = cb.getInletSpecies();
    map<string, doublereal>::iterator p = species.begin();
    //inlet species mass/mole fraction
    cout << "\n\t Inlet species mass/mole fraction...\n";
    cout << "\t----------------------------------------\n";

    while (p!= species.end())
    {
        cout <<"\t "<< p->first << "\t"<< p->second << endl;
        p++;
    }
}


void CamReporter::consoleHead(string head)
{
    cout << "\n " << head << endl;
    int len = head.length();

    for (int i = 0; i < len + 2; i++)
    {
        cout << "-";
    }

    cout << endl;
}


void CamReporter::openFile(string fileName, bool old)
{
    custom = new DataIO(fileName, old);
}


void CamReporter::closeFile()
{
    custom->close();

    delete custom;
}


void CamReporter::openFiles(bool stdrd, bool ratesOut, bool transOut)
{
    if (stdrd)
    {
        standard  = new DataIO("profile.dat");
    }
    if (ratesOut)
    {
        rates = new DataIO("rates.dat");
    }
}


void CamReporter::closeFiles(bool stdrd, bool ratesOut, bool transOut)
{
    if (stdrd)
    {
        standard->close();
    }

    if (ratesOut)
    {
        rates->close();
    }

    delete standard;
    delete rates;
}


//write a custom header
void CamReporter::writeCustomHeader(std::vector<std::string> header)
{
    if (header.empty())
    {
        throw CamError("header info missing\n");
    }

    custom->write(header);
}

//write the header for standard output data
void CamReporter::writeHeader(vector<string> stdHeader)
{
    if (stdHeader.empty())
    {
        throw CamError("header info missing\n");
    }

    standard->write(stdHeader);
}


//write header for standard data as well as the rates header
void CamReporter::writeHeader
(
    vector<string>& stdHeader,
    vector<string>& ratesHeader
)
{
    standard->write(stdHeader);
    rates->write(ratesHeader);
}


void CamReporter::writeStdFileOut(vector<doublereal>& data)
{
    standard->write(data);
}


void CamReporter::writeCustomFileOut(vector<doublereal>& data)
{
    custom->write(data);
}
