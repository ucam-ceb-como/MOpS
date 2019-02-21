/*!
 * @file    swp_titania_surface_reaction.cpp
 * @author  William J Menz
 * @brief   Implementation of the titania surface reactions
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      Implementation of the titania surface reactions
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
#include "swp_titania_surface_reaction.h"
#include "swp_mechanism.h"

// Anonymous namespace just for use in this file
namespace {

// Surface site density of TiO2 sites (#/m2)
static const double rho_s(5.677e18);
}

namespace Sweep {

namespace Processes {

// Default constructor (private)
TitaniaSurfaceReaction::TitaniaSurfaceReaction()
: SurfaceReaction(),
  m_sr_type(iFirstOrder),
  m_i_ticl4(0u),
  m_i_o2(0u)
{}

/*!
 * Mechanism constructor
 *
 * @param mech  Mechanism to create Process with
 * @return      Initialised Process
 */
TitaniaSurfaceReaction::TitaniaSurfaceReaction(
        const Sweep::Mechanism &mech,
        TitaniaSRForm form)
: SurfaceReaction(mech),
  m_sr_type(form),
  m_i_ticl4(0u),
  m_i_o2(0u)
{
    // Call initialise to assign the indices
    init(*mech.Species());
}

/*!
 * Assigns gas-phase indices to object
 *
 * @param Vector of pointers to gas-phase species
 */
void TitaniaSurfaceReaction::init(const Sprog::SpeciesPtrVector &sp)
{
    // Loop over species, find O2 and TiCl4
    for (unsigned int i = 0; i != sp.size(); ++i) {
        if (sp[i]->Name() == "O2") m_i_o2 = i;
        if ((sp[i]->Name() == "TiCl4") ||
                (sp[i]->Name() == "TICl4") ||
                (sp[i]->Name() == "TICl4"))
                    m_i_ticl4 = i;
    }

    // Check they have been assigned
    if ((m_i_ticl4 == 0u) && (m_i_o2 == 0u))
        throw std::runtime_error("Could not find TiCl4 & O2 in gas-phase "
                "in TitaniaSurfaceReaction::init");
}

/*!
 * Copy constructor
 *
 * @param copy  The object to be copied
 * @return      A copy of the object
 */
TitaniaSurfaceReaction::TitaniaSurfaceReaction(const TitaniaSurfaceReaction &copy)
{
    *this = copy;
}

/*!
 * Stream-reading constructor
 *
 * @param in    Input binary stream
 * @param mech  Mechanism to create Process with
 * @return      Initialised Process
 */
TitaniaSurfaceReaction::TitaniaSurfaceReaction(
        std::istream &in,
        const Sweep::Mechanism &mech
        )
{
    Deserialize(in, mech);
}

// Assignment operator
TitaniaSurfaceReaction &TitaniaSurfaceReaction::operator =(const TitaniaSurfaceReaction &rhs)
{
    if (this != &rhs) {
        SurfaceReaction::operator =(rhs);
        m_sr_type = rhs.m_sr_type;
        m_i_ticl4 = rhs.m_i_ticl4;
        m_i_o2    = rhs.m_i_o2;
    }
    return *this;
}

// Creates a copy of the particle process
TitaniaSurfaceReaction *const TitaniaSurfaceReaction::Clone() const
{
    return new TitaniaSurfaceReaction(*this);
}

// Returns the process type
ProcessType TitaniaSurfaceReaction::ID() const {return TitaniaSR_ID;}

/*!
 * Calculates the fraction of occupied sites, assuming quasi steady state
 * This is done for the recommended model of Preprint 100 (Model C)
 *
 * f = k1 [TiCl4] / (k2 [O2] + k1 [TiCl4])
 * where k1: A = 2.1e12 cm3/mols, Ea = 50 kJ/mol
 *       k2: A = 6.8e11 cm3/mols, Ea = 25 kJ/mol, n = 0.5
 *
 * @param gas   Gas to evaluate the site fraction for
 * @return      Occupied site fraction
 */
double TitaniaSurfaceReaction::OccupiedSiteFraction(
        const EnvironmentInterface &gas
        ) const
{
    double T  = gas.Temperature();
    double k1 = 2.10e12 * exp(50000.0 / (Sweep::R * T));
    k1 *= gas.SpeciesConcentration(m_i_ticl4);
    double k2 = 6.80e11 * exp(25000.0 / (Sweep::R * T)) * sqrt(T);
    k2 *= gas.SpeciesConcentration(m_i_o2);

    return k1 / (k1 + k2);
}

/*!
 * Returns the rate for the whole system
 *
 * @param t             Time
 * @param sys           The cell to calculate the rate for
 * @param local_geom    Location information
 * @return              Value of the rate
 */
double TitaniaSurfaceReaction::Rate(
        double t,
        const Cell &sys,
        const Geometry::LocalGeometry1d &local_geom
        ) const
{
    double rate(1.0);

    if (GetRateForm() == iMultivariate) {
        // Can't use usual infrastructure (due to PropID dependence)
        // Rate = k (S * \rho_S - \eta_{Cl}) * [TiCl4]

        rate *= m_arr.A;
        rate *= chemRatePart(sys.GasPhase());

        double T = sys.GasPhase().Temperature();
        rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));
        if (m_mech->AnyDeferred()) rate *= m_majfactor;

        // Now calculate the number of active sites [S*]
        //  = [S*]total - [TiCl4*]
        //  = S * \rho - \eta_{Cl} / 4
        rate *= std::max(sys.Particles().GetSum(Sweep::iS) * rho_s -
                         sys.Particles().GetSum(Sweep::iASN)/4.0, 0.0);

    } else {
        // Just use parent class infrastructure, hoping the correct parameters
        // have been chosen in sweep.xml
        rate *= SurfaceReaction::Rate(t, sys, local_geom);

        if (GetRateForm() == iEleyRidealAdsorption) {
            // Proportional to available site fraction (1-occupied sites)
            rate *= (1.0 - OccupiedSiteFraction(sys.GasPhase()));

        } else if (GetRateForm() == iEleyRidealDesorption) {
            // Proportional to available site fraction (1-occupied sites)
            rate *= OccupiedSiteFraction(sys.GasPhase());

        } else if (GetRateForm() == iGhoshtagore) {
            // Proportional to [O2]^0.5
            // Assume no O2 conc dependence in parent class rate function
            rate *= sqrt(sys.GasPhase().SpeciesConcentration(m_i_o2));
        }
    }

    return rate;
}

/*!
 * Calculates the single-particle rate
 *
 * @param t             Time
 * @param sys           The cell to calculate the rate for
 * @param sp            Particle for which to calculate rate
 * @return              Value of the rate
 */
double TitaniaSurfaceReaction::Rate(
        double t,
        const Cell &sys,
        const Particle &sp
        ) const
{
    double rate(1.0);

    if (GetRateForm() == iMultivariate) {
        // Can't use usual infrastructure (due to PropID dependence)
        // Rate = k (S * \rho_S - \eta_{Cl}) * [TiCl4]

        rate *= m_arr.A;
        rate *= chemRatePart(sys.GasPhase());

        double T = sys.GasPhase().Temperature();
        rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));
        if (m_mech->AnyDeferred()) rate *= m_majfactor;

        // Now calculate the number of active sites [S*]
        //  = [S*]total - [TiCl4*]
        //  = S * \rho - \eta_{Cl} / 4
        rate *= std::max(sp.Property(Sweep::iS) * rho_s -
                sp.Property(Sweep::iASN)/4.0, 0.0);

    } else {
        // Just use parent class infrastructure, hoping the correct parameters
        // have been chosen in sweep.xml
        rate *= SurfaceReaction::Rate(t, sys, sp);

        if (GetRateForm() == iEleyRidealAdsorption) {
            // Proportional to available site fraction (1-occupied sites)
            rate *= (1.0 - OccupiedSiteFraction(sys.GasPhase()));

        } else if (GetRateForm() == iEleyRidealDesorption) {
            // Proportional to available site fraction (1-occupied sites)
            rate *= OccupiedSiteFraction(sys.GasPhase());

        } else if (GetRateForm() == iGhoshtagore) {
            // Proportional to [O2]^0.5
            // Assume no O2 conc dependence in parent class rate function
            rate *= sqrt(sys.GasPhase().SpeciesConcentration(m_i_o2));
        }
    }
    return rate;
}

// aab64 Return rate constant and chemistry part for hybrid method
double TitaniaSurfaceReaction::Rate(double t, const Cell &sys) const
{
	// Rate constant.
	double rate = SurfaceReaction::Rate(t, sys);
	return rate;
}


/*!
 * Writes the object to a binary stream.
 *
 * @param out   Output binary stream
 */
void TitaniaSurfaceReaction::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Serialize base class.
        SurfaceReaction::Serialize(out);

        // Write derived class data
        unsigned int val(0);

        val = (unsigned int) m_sr_type;
        out.write((char*)&val, sizeof(val));

        val = m_i_ticl4;
        out.write((char*)&val, sizeof(val));

        val = m_i_o2;
        out.write((char*)&val, sizeof(val));

    } else {
        throw invalid_argument("Output stream not ready "
                               "in TitaniaSurfaceReaction::Serialize");
    }
}

/*!
 * Reads object from a binary stream
 *
 * @param in    Input binary stream
 * @param mech  Mechanism for process
 */
void TitaniaSurfaceReaction::Deserialize(
        std::istream &in,
        const Sweep::Mechanism &mech
        )
{
    if (in.good()) {
        // Deserialize base class.
        SurfaceReaction::Deserialize(in, mech);

        // Read derived class data
        unsigned int val(0);

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_sr_type = (TitaniaSRForm) val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_i_ticl4 = (unsigned int) val;

        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        m_i_o2 = (unsigned int) val;

    } else {
        throw invalid_argument("Input stream not ready "
                               "in TitaniaSurfaceReaction::Deserialize");
    }
}

} // Processes

} // Sweep
