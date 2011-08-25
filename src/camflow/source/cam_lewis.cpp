#include "cam_lewis.h"
#include "camxml.h"
#include <boost/algorithm/string.hpp>

using namespace Camflow;

LewisNumber::LewisNumber
(
    const std::string& inputFileName,
    const int& mCord,
    const int& nSpc
)
:
    Le(mCord,nSpc,1.0)
{

    loadSettings(inputFileName);
    if (lewisType_ == FIXEDFROMFILE) readFromFile();
    if (lewisType_ == CALCULATED) calculateLe();

}


void LewisNumber::loadSettings(const std::string& inputFileName)
{
    CamXML::Document doc;
    doc.Load(inputFileName);
    const CamXML::Element* root = doc.Root();
    CamXML::Element *subnode, *opNode;
    std::vector<CamXML::Element*> radiativeSpecies;
    std::vector<CamXML::Element*>::const_iterator p;

    opNode = root->GetFirstChild("op_condition");
    subnode = opNode->GetFirstChild("Lewis");
    if (subnode != NULL)
    {
        std::string lewisOption = subnode->Data();
        boost::to_upper(lewisOption);
        if (lewisOption == "UNITY") lewisType_ = UNITY;
        else if (lewisOption == "FIXEDFROMFILE") lewisType_ = FIXEDFROMFILE;
        else if (lewisOption == "CALCULATED") lewisType_ = CALCULATED;
        else
        {
            throw std::runtime_error
            (
                lewisOption+" does not exist."
                " Specify unity / fixedFromFile / calculated."
            );
        }
    }
    else
    {
        throw std::runtime_error
        (
            "No Lewis number setting. Specify either\n"
            " <Lewis>unity</Lewis> or "
            " <Lewis>fixedFromFile</Lewis> or "
            " <Lewis>calculated</Lewis> in <op_condition>."
        );
    }
}

void LewisNumber::readFromFile()
{

}

void LewisNumber::calculateLe()
{

}

const doublereal&
LewisNumber::operator()
(
    const int& Z,
    const int& species
)
const
{

    return Le(Z,species);

}

doublereal&
LewisNumber::operator()
(
    const int& Z,
    const int& species
)
{

    return Le(Z,species);

}
