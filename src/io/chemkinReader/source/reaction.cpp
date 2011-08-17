/*
 * reaction.cpp
 *
 *  Created on: Jun 25, 2011
 *      Author: lrm29
 */

#include "reaction.h"
#include <iomanip>
#include <stdexcept>

using namespace std;

IO::Reaction::Reaction()
:flagReversible_(true)
,flagDuplicate_(false)
,reactants_()
,products_()
,forwardArrhenius_()
,reverseArrhenius_()
,flagHasREV_(false)
,flagThirdBody_(false)
,flagLOW_(false)
,flagTROE_(false)
,flagSRI_(false)
,flagPressureDependent_(false)
,thirdBodies_()
,fallOffBody_()
{}

void IO::Reaction::setReversible(const bool flag)
{
    flagReversible_ = flag;
}
const bool& IO::Reaction::isReversible() const
{
    return flagReversible_;
}

const bool& IO::Reaction::hasREV() const
{
    return flagHasREV_;
}

void IO::Reaction::setArrhenius
(
    double A,
    double n,
    double E,
    bool reverse
)
{
    if (reverse)
    {
        flagHasREV_ = true;
        reverseArrhenius_.A = A;
        reverseArrhenius_.n = n;
        reverseArrhenius_.E = E;
    }
    else
    {
        forwardArrhenius_.A = A;
        forwardArrhenius_.n = n;
        forwardArrhenius_.E = E;
    }
}
const IO::Arrhenius& IO::Reaction::getArrhenius(bool reverse) const
{
    if (reverse)
    {
        return reverseArrhenius_;
    }
    return forwardArrhenius_;
}

void IO::Reaction::setReactants(multimap<string, double> reactants)
{
    checkForThirdBody(reactants);
    reactants_ = reactants;
}

const multimap<string, double>& IO::Reaction::getReactants() const
{
    return reactants_;
}

void IO::Reaction::setProducts(multimap<string, double> products)
{
    checkForThirdBody(products);
    products_ = products;
}

const multimap<string, double>& IO::Reaction::getProducts() const
{
    return products_;
}

void IO::Reaction::checkForThirdBody(multimap<string, double>& species)
{
    multimap<string,double>::iterator iter;
    for (iter = species.begin(); iter != species.end(); ++iter)
    {
        if(iter->first == "M")
        {
            flagThirdBody_ = true;
            species.erase(iter);
            break;
        }
    }
}

void IO::Reaction::setThirdBodies(const multimap<string, double>& thirdBodies)
{
    thirdBodies_ = thirdBodies;
}

const multimap<string, double>& IO::Reaction::getThirdBodies() const
{
    return thirdBodies_;
}

void IO::Reaction::setFallOffBody(const string& fallOffBody)
{
    fallOffBody_ = fallOffBody;
}

const string& IO::Reaction::getFallOffBody() const
{
    return fallOffBody_;
}

void IO::Reaction::setLOW(const vector<double>& LOW)
{
    flagLOW_ = true;
    LOW_ = LOW;
}

const vector<double>& IO::Reaction::getLOW() const
{
    return LOW_;
}

void IO::Reaction::setTROE(const vector<double>& TROE)
{
    if ((TROE.size() != 3) && (TROE.size() != 4))
    {
        throw runtime_error("TROE is not of size 3 or 4.");
    }

    TROE_ = TROE;
    if (TROE.size() == 3) TROE_.push_back(0);
    flagTROE_ = true;
}

const vector<double>& IO::Reaction::getTROE() const
{
    return TROE_;
}

void IO::Reaction::setSRI(const vector<double>& SRI)
{

    if ((SRI.size() != 3) && (SRI.size() != 5))
    {
        throw runtime_error("SRI is not of size 3 or 5.");
    }

    SRI_ = SRI;
    if (SRI.size() == 3)
    {
        SRI_.push_back(1.0);
        SRI_.push_back(0);
    }
    flagSRI_ = true;
}

const bool& IO::Reaction::hasSRI() const
{
    return flagSRI_;
}

const bool& IO::Reaction::hasTROE() const
{
    return flagTROE_;
}

const bool& IO::Reaction::hasLOW() const
{
    return flagLOW_;
}

const vector<double>& IO::Reaction::getSRI() const
{
    return SRI_;
}

void IO::Reaction::setPressureDependent()
{
    flagPressureDependent_ = true;
}

const bool& IO::Reaction::isPressureDependent() const
{
    return flagPressureDependent_;
}

void IO::Reaction::setDuplicate()
{
    flagDuplicate_ = true;
}

const bool& IO::Reaction::hasDuplicate() const
{
    return flagDuplicate_;
}

namespace IO
{

    ostream& operator<<(ostream& output, const Reaction& reaction)
    {
        multimap<string,double>::const_iterator iter;
        vector<double>::const_iterator iterVec;

        output << "    Reaction Data\n"
               << "    (\n"
               << "        Reversible  : " << reaction.flagReversible_ << "\n"
               << "        P Dependent : " << reaction.flagPressureDependent_ << "\n"
               << "        Duplicate : " << reaction.flagDuplicate_ << "\n"
               << "        Reactants\n"
               << "        (\n";
       for (iter = reaction.reactants_.begin(); iter != reaction.reactants_.end(); ++iter)
       {
        output << "            Name : " << setw(10) <<iter->first << " | Stoich : " << iter->second << "\n";
       }
        output << "        )\n"
               << "        Products\n"
               << "        (\n";
       for (iter = reaction.products_.begin(); iter != reaction.products_.end(); ++iter)
       {
        output << "            Name : " << setw(10) <<iter->first << " | Stoich : " << iter->second << "\n";
       }
        output << "        )\n"
               << "        Arrhenius\n"
               << "        (\n"
               << "            Forward A = " << reaction.forwardArrhenius_.A << "\n"
               << "                    n = " << reaction.forwardArrhenius_.n << "\n"
               << "                    E = " << reaction.forwardArrhenius_.E << "\n"
               << "            Reverse A = " << reaction.reverseArrhenius_.A << "\n"
               << "                    n = " << reaction.reverseArrhenius_.n << "\n"
               << "                    E = " << reaction.reverseArrhenius_.E << "\n"
               << "        )\n"
               << "        Third Bodies\n"
               << "        (\n"
               << "            Third Body : " << reaction.flagThirdBody_ << "\n";
       for (iter = reaction.thirdBodies_.begin(); iter != reaction.thirdBodies_.end(); ++iter)
       {
        output << "            Name : " << setw(10) <<iter->first << " | Efficiency Factor : " << iter->second << "\n";
       }
       for (iterVec = reaction.LOW_.begin(); iterVec != reaction.LOW_.end(); ++iterVec)
       {
        output << "            LOW : " << *iterVec << "\n";
       }
       for (iterVec = reaction.TROE_.begin(); iterVec != reaction.TROE_.end(); ++iterVec)
       {
        output << "            TROE : " << *iterVec << "\n";
       }
       for (iterVec = reaction.SRI_.begin(); iterVec != reaction.SRI_.end(); ++iterVec)
       {
        output << "            SRI : " << *iterVec << "\n";
       }
        output << "        )\n"
               << "    )";
        return output;
    }

}
