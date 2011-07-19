/*
 * element.cpp
 *
 *  Created on: Jun 22, 2011
 *      Author: lrm29
 */

#include "element.h"

using namespace std;

IO::Element::Element
(
    const string name,
    const double atomicWeight
)
:
    name_(name),
    atomicWeight_(atomicWeight)
{}

const string& IO::Element::getName() const
{
    return this->name_;
}

const double& IO::Element::getAtomicWeight() const
{
    return this->atomicWeight_;
}

namespace IO
{

    ostream& operator<<(ostream& output, const Element& element)
    {
        output << "(" << element.name_ << ", " << element.atomicWeight_ << ")";
        return output;
    }

}
