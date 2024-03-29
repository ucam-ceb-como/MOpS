#include "cam_lewis.h"
#include "camxml.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace Camflow;

LewisNumber::LewisNumber
(
    const std::string& inputFileName,
    const Sprog::Mechanism *const mech,
    const int& mCord,
    const int& nSpc
)
:
    mCord_(mCord),
    mech_(mech),
    Le(mCord,nSpc,1.0)
{
    loadSettings(inputFileName);
    if (lewisType_ == FIXEDFROMFILE) readFromFile("LewisNumbersInput.xml");
    //if (lewisType_ == CALCULATED) calculateLe();
}


void LewisNumber::loadSettings(const std::string& inputFileName)
{
    CamXML::Document doc;
    doc.Load(inputFileName);
    const CamXML::Element* root = doc.Root();
    CamXML::Element *subnode, *opNode;

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

    // Read in which soot Flamlet model to use
    subnode = opNode->GetFirstChild("SootFlameletType");
    if (subnode != NULL)
    {
        std::string sootFlameletOption = subnode->Data();
        boost::to_upper(sootFlameletOption);
        if (sootFlameletOption == "MAUSS06") sootFlameletType_ = MAUSS06;
        else if (sootFlameletOption == "PITSCH00DD") sootFlameletType_ = PITSCH00DD;
        else if (sootFlameletOption == "CARBONELL09") sootFlameletType_ = CARBONELL09;
        else if (sootFlameletOption == "EXTENDEDLAGRANGIAN") sootFlameletType_ = EXTENDEDLAGRANGIAN;
        else
                {
                    throw std::runtime_error
                    (
                    	sootFlameletOption+" does not exist."
                        " Specify MAUSS06 / PITSCH00DD / CARBONELL09 / EXTENDEDLAGRANGIAN"
                    );
                }

    }
    else
    {
        throw std::runtime_error
        (
            "No soot flamelet model setting. Specify either\n"
            " <SootFlameletType>MAUSS06</SootFlameletType> or "
            " <SootFlameletType>PITSCH00DD</SootFlameletType> or "
            " <SootFlameletType>EXTENDEDLAGRANGIAN</SootFlameletType> or "
            " <SootFlameletType>CARBONELL09</SootFlameletType> in <op_condition>."
        );
    }


    // This might not be the best place for this.
    // Put it here for now as Lewis number is tied in with Flamelet calcs.
    subnode = opNode->GetFirstChild("sootFlameTimeCutoff");
    if (subnode != NULL)
    {
    	sootFlameTimeThreshold = Strings::cdble(subnode->Data());
    }
    else
    {
        throw std::runtime_error
        (
            "Soot flametime cutoff not specified"
            "e.g. <sootFlameTimeCutoff>1e-3</sootFlameTimeCutoff>"
        	"Place in <op_condition>"
        );
    }

}

void LewisNumber::readFromFile(const std::string& fixedLewisFile)
{
    CamXML::Document doc;
    doc.Load(fixedLewisFile);
    const CamXML::Element* root = doc.Root();
    CamXML::Element *subnode;
    std::vector<CamXML::Element*> subsubnodes;
    std::vector<CamXML::Element*>::const_iterator p;

    subnode = root->GetFirstChild("Lewis");

    std::cout << "Lewis Numbers read from " << fixedLewisFile << std::endl;
    if (subnode != NULL)
    {

        subnode->GetChildren("species",subsubnodes);
        for (p=subsubnodes.begin(); p<subsubnodes.end(); ++p)
        {
            std::string speciesName = (*p)->GetAttributeValue("name");
            const int speciesIndex = mech_->FindSpecies(speciesName);
            if (speciesIndex == -1)
            {
                throw std::runtime_error(
                    "You have specified radiative species " + speciesName + " that "
                    "does not exist in the mechanism!");
            }
            else
            {
                double LeNumber = boost::lexical_cast<double>((*p)->Data());
                std::cout << "    " << speciesName << " " << LeNumber << std::endl;
                for (size_t i=0; i<mCord_; ++i) Le(i,speciesIndex) = LeNumber;
            }
        }
    }

}

//!Not used yet
void LewisNumber::calculateLe()
{

}

double
LewisNumber::operator()
(
    const int& Z,
    const int& species
)
const
{

    return Le(Z,species);

}

double&
LewisNumber::calcLewis
(
    const int& Z,
    const int& species
)
{

    return Le(Z,species);

}
