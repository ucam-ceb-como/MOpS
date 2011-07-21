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

void CamReporter::header(std::string prog){
    std::cout << std::endl;
    std::cout << "!----------------------------------------------------------!\n";
    std::cout << "!               CamFlow  " << prog << " Version 1.0                !\n";
    std::cout << "!        This program is distributed without any warranty  !\n";
    std::cout << "!            Author: Vinod M. J. (vj231@cam.ac.uk)         !\n";
    std::cout << "!----------------------------------------------------------!\n";

}
void CamReporter::problemDescription(CamBoundary& cb, CamResidual& cr){
        /* screen out put for simulation case
     */
    std::cout << "\t Problem description.....\n";
    std::cout << "\t Number of species: " << cr.getNSpecies() << std::endl;
    std::cout << "\t Total number of equations solving: " << cr.getNEqn() << std::endl;
    std::cout << "\t Total number of variables solving: " << cr.getNVar() << std::endl;

    //initializing solution vector
    std::cout << "\t Initializing solution vector...\n";
    std::cout << "\t Inlet velocity: "<<cb.getVelocity() <<  " m/s" << std::endl;
    std::cout << "\t Inlet temperature "<< cb.getTemperature() << " K"<< std::endl;
    std::map<std::string, doublereal> species = cb.getInletSpecies();   
    std::map<std::string, doublereal>::iterator p = species.begin();
    //inlet species mass/mole fraction
    std::cout << "\n\t Inlet species mass/mole fraction...\n";
    std::cout << "\t----------------------------------------\n";
    while (p!= species.end()){
        std::cout <<"\t "<< p->first << "\t"<< p->second << std::endl;
        p++;

    }

}

void CamReporter::consoleHead(std::string head){
    std::cout << "\n";
    std::cout << " " << head << std::endl;
    int len = head.length();
    for (int i = 0; i < len+2; i++) {
        std::cout << "-";
    }
    std::cout << std::endl;
}
void CamReporter::openFile(std::string fileName, bool old){
    custom = new DataIO(fileName, old);
}
void CamReporter::closeFile(){
    custom->close();
}

void CamReporter::openFiles(bool stdrd, bool ratesOut, bool transOut){
    if(stdrd) standard  = new DataIO("profile.dat");
    if(ratesOut) rates = new DataIO("rates.dat");
}

void CamReporter::closeFiles(bool stdrd, bool ratesOut, bool transOut){
    if(stdrd) standard->close();
    if(ratesOut) rates->close();
}
//write a custom header
void CamReporter::writeCustomHeader(std::vector<std::string> header){
    if(header.size() == 0)throw CamError("header info missing\n");
    custom->write(header);
}
//write the header for standard output data
void CamReporter::writeHeader(std::vector<std::string> stdHeader){
    if(stdHeader.size() == 0)throw CamError("header info missing\n");
    standard->write(stdHeader);
}

//write header for standard data as well as the rates header
void CamReporter::writeHeader(std::vector<std::string>& stdHeader, std::vector<std::string>& ratesHeader){
    standard->write(stdHeader);
    rates->write(ratesHeader);
}

//write header for standard data, rates as well as transport data
void CamReporter::wrteHeader(std::vector<std::string>& stdHeader, std::vector<std::string>& ratesHeader, std::vector<std::string>& transHeader){
    standard->write(stdHeader);
    rates->write(ratesHeader);
    transport->write(transHeader);
}

void CamReporter::writeStdFileOut(std::vector<doublereal>& data){
    standard->write(data);
}
void CamReporter::writeCustomFileOut(std::vector<doublereal>& data){
    custom->write(data);
}

void CamReporter::writeTempProfiletoXML
(
    const std::string fileName,
    const doublereal sdr,
    const std::vector<doublereal>& axPos,
    const std::vector<doublereal>& temperature
)
{

    CamXML::Document doc;
    const CamXML::Element* node = NULL;

    if(doc.Load(fileName + ".xml") == 0){
        node = doc.Root();
    }

    Utils::LinearInterpolator<doublereal, doublereal> mTempInterpolator(axPos, temperature);

    CamXML::Element *initialize;
    CamXML::Element* subsubnode;
    std::vector<CamXML::Element*> subsubnodes;
    std::vector<CamXML::Element*>::const_iterator p;
    initialize = node->GetFirstChild("initialize");
    if(initialize != NULL)
    {
        //temperature profile
        subsubnode = initialize->GetFirstChild("Tprofile");
        if(subsubnode != NULL)
        {
            subsubnode->GetChildren("position",subsubnodes);

            for(p=subsubnodes.begin(); p<subsubnodes.end(); ++p)
            {
                doublereal temp = mTempInterpolator.interpolate(Strings::cdble((*p)->GetAttributeValue("x")));
                (*p)->SetData(Strings::cstr(temp));
            }
        }
        else
        {
            initialize->AddChild("Tprofile");
            subsubnode = initialize->GetFirstChild("Tprofile");
            subsubnode->SetAttribute("unit_L","m");
            subsubnode->SetAttribute("unit_T","K");
            for(size_t q=0; q<axPos.size(); ++q)
            {
                subsubnode->AddChild("position");
                subsubnode->GetChildren("position",subsubnodes);
                subsubnodes[q]->SetAttribute("x",Strings::cstr(axPos[q]));
                doublereal temp = mTempInterpolator.interpolate(Strings::cdble(subsubnodes[q]->GetAttributeValue("x")));
                subsubnodes[q]->SetData(Strings::cstr(temp));
            }
        }
    }

    std::string outputFile = "sdr_" + Strings::cstr(sdr) + ".xml";

    doc.Save("camflow.xml");

}

void CamReporter::writeMassFracProfiletoXML
(
    const std::string fileName,
    const doublereal sdr,
    const std::vector<doublereal>& axPos,
    const std::vector<doublereal>& massFrac,
    const std::string speciesName
)
{
    CamXML::Document doc;
    CamConverter convert;
    const CamXML::Element* node = NULL;

    std::string outputFile = "sdr_" + Strings::cstr(sdr) + ".xml";

    if(doc.Load(fileName + ".xml") == 0){
        node = doc.Root();
    }

    Utils::LinearInterpolator<doublereal, doublereal> mMassFracInterpolator(axPos, massFrac);

    CamXML::Element *initialize;
    std::vector<CamXML::Element*> subsubnode;
    std::vector<CamXML::Element*> subsubnodes;
    std::vector<CamXML::Element*>::const_iterator p, q;
    initialize = node->GetFirstChild("initialize");
    if(initialize != NULL)
    {

        initialize->GetFirstChild("massfracs");
        initialize->GetChildren("massfracs",subsubnode);

        bool foundSpecies = false;

        for(p=subsubnode.begin(); p<subsubnode.end(); ++p)
        {
            const CamXML::Attribute *species;
            species = (*p)->GetAttribute("species");
            if(species->GetValue() == speciesName)
            {
                foundSpecies = true;

                (*p)->GetChildren("position",subsubnodes);
                for (q=subsubnodes.begin(); q<subsubnodes.end(); ++q)
                {
                    doublereal temp = mMassFracInterpolator.interpolate(Strings::cdble((*q)->GetAttributeValue("x")));
                    (*q)->SetData(Strings::cstr(temp));
                }

                doc.Save("camflow.xml");

            }

        }

        if(!foundSpecies)
        {
            initialize->AddChild("massfracs");
            initialize->GetChildren("massfracs",subsubnode);

            size_t r = subsubnode.size()-1;

            subsubnode[r]->SetAttribute("species",speciesName);
            for(size_t s=0; s<axPos.size(); ++s)
            {
                subsubnode[r]->AddChild("position");
                subsubnode[r]->GetChildren("position",subsubnodes);
                subsubnodes[s]->SetAttribute("x",Strings::cstr(axPos[s]));
                doublereal temp = mMassFracInterpolator.interpolate(Strings::cdble(subsubnodes[s]->GetAttributeValue("x")));
                subsubnodes[s]->SetData(Strings::cstr(temp));
            }

            doc.Save("camflow.xml");

        }

    }

}
