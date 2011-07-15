/*
 * transport.cpp
 *
 *  Created on: Jun 24, 2011
 *      Author: Laurence R. McGlashan
 */

#include "transport.h"
#include "gpc_params.h"
#include "fast_math_functions.hpp"

using namespace std;

IO::Transport::Transport()
:
    moleculeIndex_(-1),
    potentialWellDepth_(-1.0),
    collisionDiameter_(-1.0),
    dipoleMoment_(-1.0),
    polarizability_(-1.0),
    rotRelaxationNumber_(-1.0),
    reducedDipoleMoment_(-1.0)
{}

void IO::Transport::setMoleculeIndex(const int moleculeIndex)
{
   this->moleculeIndex_ = moleculeIndex;
}

const int& IO::Transport::getMoleculeIndex() const
{
    return this->moleculeIndex_;
}

void IO::Transport::setCollisionDiameter(const double collisionDiameter)
{
   this->collisionDiameter_ = collisionDiameter*Sprog::Angstroem__;
}

const double& IO::Transport::getCollisionDiameter() const
{
    return this->collisionDiameter_;
}

void IO::Transport::setDipoleMoment(const double dipoleMoment)
{
   this->dipoleMoment_ = dipoleMoment*Sprog::Debye__;
}

const double& IO::Transport::getDipoleMoment() const
{
    return this->dipoleMoment_;
}

void IO::Transport::setPolarizability(const double polarizability)
{
   this->polarizability_ = polarizability*fastMath::pow3(Sprog::Angstroem__);
}

const double& IO::Transport::getPolarizability() const
{
    return this->polarizability_;
}

void IO::Transport::setPotentialWellDepth(const double potentialWellDepth)
{
   this->potentialWellDepth_ = potentialWellDepth*Sprog::kB;
}

const double& IO::Transport::getPotentialWellDepth() const
{
    return this->potentialWellDepth_;
}

void IO::Transport::setRotRelaxationNumber(const double rotRelaxationNumber)
{
   this->rotRelaxationNumber_ = rotRelaxationNumber;
}

const double& IO::Transport::getRotRelaxationNumber() const
{
    return this->rotRelaxationNumber_;
}

void IO::Transport::setReducedDipoleMoment()
{
   this->reducedDipoleMoment_
   =
   fastMath::pow2(this->dipoleMoment_)
   /
   (
       8
       *Sprog::PI
       *Sprog::EPSILON0
       *this->potentialWellDepth_
       *fastMath::pow3(this->collisionDiameter_)
   );
}

void IO::Transport::setReducedDipoleMoment(const double reducedDipoleMoment)
{
   this->reducedDipoleMoment_ = reducedDipoleMoment;
}

const double& IO::Transport::getReducedDipoleMoment()
{
    return this->reducedDipoleMoment_;
}

namespace IO {

    ostream& operator<<(ostream& output, const Transport& transport) {
        output << "    Transport Data:\n"
               << "    (\n"
               << "        Mol. Index      : " << transport.moleculeIndex_ << "\n"
               << "        Pot. Well Depth : " << transport.potentialWellDepth_ << "\n"
               << "        Collision Diam. : " << transport.collisionDiameter_ << "\n"
               << "        Dipole Moment   : " << transport.dipoleMoment_ << "\n"
               << "        Polarizability  : " << transport.polarizability_ << "\n"
               << "        Rot. Relax. No. : " << transport.rotRelaxationNumber_ << "\n"
               << "    )";
        return output;
    }
}
