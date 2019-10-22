/*!
 * \file   swp_weighted_transcoag.h
 * \author William J Menz
 *  Copyright (C) 2010 William J Menz
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 * 
 * \brief  Definition of the weighted transition regime kernel.
 *

  File purpose:
    Definition of the coagulation process.  This coagulation process uses a transition
    kernel, calculating the rates for the free-molecular and slip-flow kernels. The
    kernel is modified to include statstical weights.

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

#ifndef SWP_WEIGHTED_TRANSCOAG_H
#define	SWP_WEIGHTED_TRANSCOAG_H

#include "swp_coagulation.h"
#include "swp_process_type.h"

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

// Forward declare class used for sums in the binary tree
class TreeTransCoagWeightedCache;

namespace Processes
{

class WeightedTransitionCoagulation : public Coagulation
{
public:
    //! Default constructors.
    WeightedTransitionCoagulation(const Sweep::Mechanism &mech, const CoagWeightRule weight_rule);

    //! Deserialisation
    WeightedTransitionCoagulation(
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    //! Virtual destructor
    virtual ~WeightedTransitionCoagulation() {};

    //! Clone object
    virtual WeightedTransitionCoagulation* const Clone() const;

    //! Returns the process type for identification during serialisation
    virtual ProcessType ID(void) const {return Weighted_Transition_Coagulation_ID;};

    // TOTAL RATE CALCULATION.

    //! Returns the rate of the process for the given system.
    virtual double Rate(double t,         // Time.
                      const Cell &sys, // System for which to calculate rate.
                      const Geometry::LocalGeometry1d &local_geom
                      ) const;

    // RATE TERM CALCULATION.

    //! Returns the number of rate terms for this process.
    virtual unsigned int TermCount(void) const;

    //! Calculate rate terms given an iterator to a double vector, return sum of all terms.
    virtual double RateTerms(
        double t,                  // Time.
        const Cell &sys,       // Indicates true kernel (not majorant).
        const Geometry::LocalGeometry1d &local_geom,
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

    //! Perform a coagulation with particles chosen according to the transition kernel
    virtual int Perform(
        double t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng) const;

    //! Write the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

protected:
    //! Transition coagulation kernel between two particles
    virtual double CoagKernel(
        const Particle &sp1, // First particle.
        const Particle &sp2, // Second particle.
        const Cell &sys
        ) const;

    //! Majorant coagulation kernel between two particles
    virtual double MajorantKernel(
        const Particle &sp1, // First particle.
        const Particle &sp2, // Second particle.
        const Cell &sys,
        const MajorantType maj) const;

private:

    //! Calculate the individual rate terms for the weighted transition kernel
    double RateTerms(
    const TreeTransCoagWeightedCache &data, // Particle model data.
    double n,     // Number of particles.
    double sqrtT, // Square root of the temperature
    double T_mu,  // T / viscosity
    double MFP,   // Gas mean-free path.
    double vol,   // System sample volume.
    fvector::iterator &iterm // Iterator to first coagulation term.
    ) const;

    // COAGULATION KERNEL ROUTINES.

    //! Returns the free-molecular kernel for two given particles.
    double FreeMolKernel(
    const Particle &sp1, // First particle.
    const Particle &sp2, // Second particle.
    double T,              // Temperature.
    double P,              // Pressure.
    const bool maj       // true=majorant kernel, false=true kernel.
    ) const;

    //! Returns the slip-flow kernel for two given particles.
    double SlipFlowKernel(
    const Particle &sp1, // First particle.
    const Particle &sp2, // Second particle.
    double T,              // Temperature.
    double P,              // Pressure.
    double mu,             // Viscosity
    const bool maj       // true=majorant kernel, false=true kernel.
    ) const;

    //! Total number of terms for coagulation property choices
    static const unsigned int TYPE_COUNT = 11;
    
    //! Enum of terms defining possible coagulation particle propety choices
    enum TermType {
        FreeMol1,
        FreeMol2,
        FreeMol3,
        FreeMol4,
        SlipFlow1,
        SlipFlow2,
        SlipFlow3,
        SlipFlow4,
        SlipFlow5,
        SlipFlow6,
        SlipFlow7
    };

    //! Free-molecular enhancement factor.
    const double m_efm;

    //! Specify what to do with weights on coagulation
    CoagWeightRule m_CoagWeightRule;


}; // class WeightedTransitionCoagulation

} // namespace Processes

} // namespace Sweep

#endif	/* SWP_WEIGHTED_TRANSCOAG_H */
