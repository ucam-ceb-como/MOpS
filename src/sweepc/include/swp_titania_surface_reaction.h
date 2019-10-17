/*!
 * @file    swp_titania_surface_reaction.h
 * @author  William J Menz
 * @brief   Different types of surface reactions for titania particles
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      There has been considerable effort dedicated to characterising the
 *      surface reaction of TiCl4 on TiO2 surfaces. This file provides support
 *      for several titania surface reaction mechanisms.
 *
 *      1) First-order
 *      This is a mechanism first-order in O2 concentration (as defined in
 *      Preprint 99) and is essentially no different from a normal
 *      surface reaction process.
 *
 *      2) Ghoshtagore
 *      This assumes the rate is proportional to [02]^0.5, although
 *      stoichiometrically 1 O2 units is removed from the gas phase per Ti unit
 *
 *      3) EleyRidealAdsorption & EleyRidealDesorption
 *      Preprint 100 highlights an Eley-Rideal mechansim for the growth of
 *      TiCl4 on TiO2. EleyRidealAdsorption refers to the first reaction:
 *      TiCl4 + S* -> TiCl4*
 *
 *      EleyRidealDesorption refers to the second:
 *      TiCl4* + O2 -> TiO2 + 2Cl2 + S*
 *      The rate expressions for these processes are calculated assuming
 *      pseudo steady state of active sites (same as ABF model).
 *      These are intended to be used in the same sweep.xml file.
 *
 *      4) Multivariate
 *      Assuming the number of sites on a particle are represented by the
 *      number of Cl atoms, the rate of the first reaction above is
 *      estimated using:
 *      R = k1 (rho_s * S - \eta_{Cl}/4) * [TiCl4]
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
 *      sweepc is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public License
 *      as published by the Free Software Foundation; either version 2
 *      of the License, or (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this program; if not, write to the Free Software
 *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *      02111-1307, USA.
 *
 *   Contact:
 *      Prof Markus Kraft
 *      Dept of Chemical Engineering
 *      University of Cambridge
 *      New Museums Site
 *      Pembroke Street
 *      Cambridge
 *      CB2 3RA, UK
 *
 *      Email:       mk306@cam.ac.uk
 *      Website:     http://como.cheng.cam.ac.uk
*/
#ifndef SWP_TITANIA_SURFACE_REACTION_H_
#define SWP_TITANIA_SURFACE_REACTION_H_

#include "swp_surface_reaction.h"

namespace Sweep {

namespace Processes {

class TitaniaSurfaceReaction: public Sweep::Processes::SurfaceReaction {
public:
    //! Definitions of reactions types
    enum TitaniaSRForm
        {iFirstOrder, iGhoshtagore,
         iEleyRidealAdsorption, iEleyRidealDesorption,
         iMultivariate};

    //! Mechanism constructor
    TitaniaSurfaceReaction(const Sweep::Mechanism &mech, TitaniaSRForm form);

    //! Copy constructor
    TitaniaSurfaceReaction(const TitaniaSurfaceReaction &copy);

    //! Stream reading constructor
    TitaniaSurfaceReaction(
        std::istream &in,
        const Sweep::Mechanism &mech
        );

    //! Set reaction form
    void SetReactionForm(TitaniaSRForm form) {m_sr_type = form;}

    //! Get reaction form
    TitaniaSRForm GetRateForm(void) const {return m_sr_type;}

    //! Creates a copy of the particle process
    TitaniaSurfaceReaction *const Clone(void) const;

    //! Assignment operator
    TitaniaSurfaceReaction &operator=(const TitaniaSurfaceReaction &rhs);

    //! Returns the process type
    ProcessType ID(void) const;

    //! Returns rate of the process for the given system
    double Rate(
        double t,
        const Cell &sys,
        const Geometry::LocalGeometry1d &local_geom
        ) const;

    //! Returns single-particle rate of the process
    double Rate(
        double t,
        const Cell &sys,
        const Particle &sp
        ) const;

    // Return just the rate constant and chemistry part for hybrid method
    double Rate(
        double t,
        const Cell &sys) const;

    //! Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    //! Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,
        const Sweep::Mechanism &mech
        );

private:
    //! Private default constructor
    TitaniaSurfaceReaction();

    //! Returns the fraction of occupied sites
    double OccupiedSiteFraction(const EnvironmentInterface &gas) const;

    //! Initialise indices for gas-phase species
    void init(const Sprog::SpeciesPtrVector &sp);

    //! Reaction type selected
    TitaniaSRForm m_sr_type;

    //! Index of TiCl4 in the gas-phase
    unsigned int m_i_ticl4;

    //! Index of O2 in the gas-phase
    unsigned int m_i_o2;

};

}

}

#endif /* SWP_TITANIA_SURFACE_REACTION_H_ */
