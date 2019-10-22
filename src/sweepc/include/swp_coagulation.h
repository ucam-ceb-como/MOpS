/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of the coagulation process.  This coagulation process uses a transition
    kernel (REF).  It calculates the rates for the free molecular and slip-flow regimes.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#ifndef SWEEP_COAGULATION_H
#define SWEEP_COAGULATION_H

#include "swp_process.h"
#include "swp_coag_weight_rules.h"

#include <iostream>
    
namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Processes
{
    //Forward declaration to make typedef possible
    class Coagulation;

    //! Vector of polymorphic coagulation processes
    typedef std::vector<Coagulation*> CoagPtrVector;

/*!
 * \brief Processes that stick two particles together to form one
 */
class Coagulation : public Process
{
public:

    //! Ordinary method of construction just passes argument through.
    Coagulation(const Sweep::Mechanism &mech);

    //! Returns a copy of the coagulation process
    virtual Coagulation *const Clone(void) const = 0;

    //! Calculate rates for a sequence of coagulation processes
    static double CalcRates(double t, const Cell &sys,
                          const Geometry::LocalGeometry1d &local_geom,
                          const CoagPtrVector &coags,
                          fvector &rates, unsigned int start = 0);

    //! Calculate rate terms for a sequence of coagulation processes
    static double CalcRateTerms(double t, const Cell &sys,
                              const Geometry::LocalGeometry1d &local_geom,
                              const CoagPtrVector &coags,
                              fvector::iterator &iterm);

    //! Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    //! Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

    //! Instructions on how to choose position of coagulated particle in spatially resovled simulations
    enum ParticlePositionChoice {
        //! Do whatever was most convenient in the homogeneous case (may introduce bias)
        NoPositionChoice,

        //! Choose randomly with equal probabilities between the positions of the two particles
        UniformPositionChoice,

        //! Choose randomly between the positions of the incoming particles with probabilities proportional to their mass
        MassPositionChoice,

        //! Always take the position of the largest particle
        LargestMassPositionChoice,

        //! Take the midpoint position
        MidpointPositionChoice,

        //! Centre of mass
        CentreOfMassPositionChoice,
    };

    //! Rule for choosing post coagulation position (not relevant to homogeneous sims)
    ParticlePositionChoice PositionChoiceRule() const {return mPositionChoice;}

    //! Set rule for choosing post coagulation position (not relevant to homogeneous sims)
    void SetPositionChoiceRule(const ParticlePositionChoice rule) {mPositionChoice = rule;}

    // Majorant types are important for the transition regime kernel,
    // but additional kernels are free to add addition values to this
    // enum.  Most kernels will ignore most enum values.
    enum MajorantType {
        Default, // Place holder value
        FreeMol, // Free-molecular majorant.
        SlipFlow // Slip-flow majorant.
    };

protected:

    /*!
     * @brief  Default constructor is protected to prevent coagulations being
     *         defined without knowledge of the parent mechanism.  An architecture
     *         that removed this restriction would be nice, so that there are
     *         fewer pointers connecting apparently separate objects.
     */
    Coagulation(void) {};
    
    //! Actually stick two particles together (intended for DSA use)
    int JoinParticles(const double t, const int ip1, Particle *sp1,
                      const int ip2, Particle *sp2,
                      Cell &sys, rng_type &rng) const;

    //! Select two particles and stick them together in a weighted particle event
    int WeightedPerform(const double t, const Sweep::PropID prop1,
                        const Sweep::PropID prop2,
                        const Sweep::Processes::CoagWeightRule weight_rule,
                        Cell &sys, rng_type &rng,
                        Sweep::Processes::Coagulation::MajorantType maj) const;

    //! Calculate kernel between two particles
    virtual double CoagKernel(const Particle &sp1, const Particle &sp2,
                    const Cell& sys) const = 0;



    //! Calculate majorant kernel between two particles
    virtual double MajorantKernel(const Particle &sp1, const Particle &sp2,
                                const Cell& sys, const MajorantType maj) const = 0;
private:

    //! Rule for determining position of particle after coagulation (only relevant for spatial sims)
    ParticlePositionChoice mPositionChoice;
};

} //namespace Processes
} //namespace Sweep

#endif
