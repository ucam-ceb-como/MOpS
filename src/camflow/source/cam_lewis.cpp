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
    if (lewisType_ == CONSTANTCALCULATED) calculateConstantLe();
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
        else if (lewisOption == "CONSTANTCALCULATED") lewisType_ = CONSTANTCALCULATED;
        else if (lewisOption == "CALCULATED") lewisType_ = CALCULATED;
        else
        {
            throw std::runtime_error
            (
                lewisOption+" does not exist."
                " Specify unity / fixedFromFile / calculated / constantcalculated (vs N2 at 1800K)."
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

    
    //Moved SootFLameletType to cam_soot.cpp
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

void LewisNumber::calculateConstantLe()
{
    // <Lewis>CONSTANTCALCULATED</Lewis>

    int numSpecies = mech_->SpeciesCount();

    // Create mixture of hot nitrogen to get Lewis number for species
    Sprog::Thermo::IdealGas* tmpMixture = new Sprog::Thermo::IdealGas(mech_->Species());

    // Set mixture properties: NB Set T before P
    tmpMixture->SetTemperature(1800.0); // K
    tmpMixture->SetPressure(101325.0); // Pa (Actually sets Density using ideal gas law)
    
    // set composition of Pure N2
    const int iN2 = mech_->FindSpecies("N2");
    std::vector<double> tmpMixComposition; 
    int j=0;
    for (j=0; j<numSpecies; ++j)
    {
        if (j==iN2)
        {
            tmpMixComposition.push_back(1.0);
        }
        else   
        {
            tmpMixComposition.push_back(0.0);
        }
    }
    tmpMixture -> SetMassFracs(tmpMixComposition);
    tmpMixture -> SetFracs(tmpMixComposition);
    std::cout << "Calculating Lewis Numbers in a N2 mixture at 101325 Pa and 1800K..."<< std::endl;
    std::cout << "Pressure    :" << tmpMixture->Pressure() << std::endl;
    std::cout << "Temperature :" << tmpMixture->Temperature() << std::endl;
    std::cout << "Denisty     : " << tmpMixture -> Density() << "  mol / m3" <<std::endl;
    std::cout << "k           : " << tmpMixture -> getThermalConductivity(tmpMixture->Pressure()) << std::endl;
    std::cout << "Cp          : " << tmpMixture -> getSpecificHeatCapacity() << std::endl; 
    //Matches expected NIST data after unit conversion.


    Sprog::Mechanism::const_sp_iterator spIt = mech_->SpBegin();
    const Sprog::Mechanism::const_sp_iterator spItEnd = mech_->SpEnd();
    
    //N2 gas constants 
    double k=tmpMixture->getThermalConductivity(tmpMixture->Pressure()); // Pa
    double rho=tmpMixture -> Density(); //mol/m3
    double cp=tmpMixture->getSpecificHeatCapacity() // J/kg/K
        * tmpMixture->getAvgMolWt(tmpMixComposition); // kg/mol
    double k_rho_cp=k/rho/cp;
    

    std::cout << "Calculating and generating LewisNumbersInput.xml text to be used with\n" << "<Lewis>fixedFromFile</Lewis>" <<std::endl;
    std::cout << "<Camflow>\n    <Lewis>" <<std::endl;
    while(spIt != spItEnd) {
        // Dereference unincremented value
        std::string speciesName = (*spIt)->Name();
        const int speciesIndex = mech_->FindSpecies(speciesName);
        
        //std::cout << "SpeciesIndex = " << speciesIndex << std::endl;
       
        // declare Lewis Number 
        double LeNum=0.0;
          
        //Binary Diffusion Coefficient into Hot N2 Gas
        if(speciesIndex == iN2)
        {
            LeNum=1;
        }
        else
        {  
            Sprog::Transport::MixtureTransport* mt;
            double Dab= mt->binaryDiffusionCoeff(
                speciesIndex,
                iN2,
                tmpMixture->Temperature(),
                tmpMixture->Pressure(),
                (*tmpMixture) //dereferenced
            );
            //Calculate Lewis Number
            LeNum = (k_rho_cp)/(Dab);
        }
        for (size_t i=0; i<mCord_; ++i) Le(i,speciesIndex) = LeNum;
        std::cout << "        <species name=\""<<speciesName <<"\">" << LeNum << "</species> " <<std::endl;

        *spIt++;
    }
    std::cout << "    </Lewis>\n</Camflow>" <<std::endl;
    //~tmpMixture();
    std::cout << "Finished loading Lewis Number" <<std::endl;
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
