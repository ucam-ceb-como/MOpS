/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of a surface reaction process which includes a term
    for active sites.  Active site densities are calculated using an
    ActiveSitesModel, which must be set before the reaction is used.

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

#ifndef SWEEP_ACTSITE_RXN_H
#define SWEEP_ACTSITE_RXN_H

#include "swp_params.h"
#include "swp_surface_reaction.h"
#include "swp_environment_interface.h"

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

namespace Processes
{
/*!
 * \brief Surface reactions that require an additional model for surface activity
 */
class ActSiteReaction : public SurfaceReaction
{
public:
    //! Rule for calculating what fraction of surface Hydrogens are radicals
    enum RadicalSiteFractionModel {
        ABFRadicalSiteModel, //! ABF see radicalSiteFractionABF
        BPRadicalSiteModel,  //! Blanquart & Pitsch see radicalSiteFractionBP
        ABFConstant          //! ABF with alpha = 1.0
    };

    // Constructors.
    ActSiteReaction(const Sweep::Mechanism &mech,
                    const RadicalSiteFractionModel rad_site_model,
                    const EnvironmentInterface::PropertyIndex alpha_index); // Default constructor.

    ActSiteReaction(                 // Stream-reading constructor.
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    // TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

    // Returns rate of the process for the given system.
    virtual double Rate(
        double t,         // Time.
        const Cell &sys, // System for which to calculate rate.
        const Geometry::LocalGeometry1d &local_geom // cell location information
        ) const;


	// SINGLE PARTICLE RATE CALCULATIONS.

    // Returns the rate of the process for the given particle in
    // the system. Process must be linear in particle number.
    virtual double Rate(
        double t,             // Current time (s).
        const Cell &sys,    // System to which the particle belongs.
        const Particle &sp  // Particle for which to calculate rate.
        ) const;

    //! Concentration of surface sites participating in this reaction (?units)
    double SiteDensity(const Cell &sys) const;

    //! Index for looking up an active sites factor in the gas phase
    void setActiveSitesIndex(const EnvironmentInterface::PropertyIndex index) {mAlphaIndex = index;}

    // READ/WRITE/COPY.
    
    // Creates a copy of the particle process.
    virtual ActSiteReaction *const Clone(void) const;

    // Returns the process type.  Used to identify different
    // processes and for serialisation.
    virtual ProcessType ID(void) const;

    // Writes the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    virtual void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    // Default constructor is protected to prevent reactions being
    // defined without knowledge of the parent mechanism.
    ActSiteReaction(void);

    //! Fraction of surface hydrogen sites that are radicals in the ABF model
    double radicalSiteFractionABF(const EnvironmentInterface &gas) const;

    //! Fraction of surface hydrogen sites that are radicals according to Blanquart & Pitsch
    double radicalSiteFractionBP(const EnvironmentInterface &gas) const;

private:
    //! Index of acetylene in mechanism and gas phase mixture data
    int iC2H2;

    //! Index of O2 in mechanism and gas phase mixture data
    int iO2;

    //! Index of OH in mechanism and gas phase mixture data
    int iOH;

    //! Index of CO in mechanism and gas phase mixture data
    int iCO;

    //! Index of H in mechanism and gas phase mixture data
    int iH;

    //! Index of H2 in mechanism and gas phase mixture data
    int iH2;

    //! Index of water in mechanism and gas phase mixture data
    int iH2O;

    //! Index of active sites fraction in gas phase mixture data
    EnvironmentInterface::PropertyIndex mAlphaIndex;

    //! Specify how to calculate radical site fraction
    RadicalSiteFractionModel mRadicalSiteModel;
};
}
}

#endif
