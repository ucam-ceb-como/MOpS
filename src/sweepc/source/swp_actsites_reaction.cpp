/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ActSiteReaction class declared in the
    swp_actsites_reaction.h header file.

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

#include "swp_actsites_reaction.h"
#include "swp_mechanism.h"
#include "swp_process_type.h"
#include "swp_model_factory.h"

#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor is protected to prevent reactions being
// defined without knowledge of the parent mechanism.
ActSiteReaction::ActSiteReaction(void)
: SurfaceReaction()
, iC2H2(-1)
, iO2(-1)
, iOH(-1)
, iCO(-1)
, iH(-1)
, iH2(-1)
, iH2O(-1)
, mAlphaIndex(-1)
, mRadicalSiteModel(static_cast<RadicalSiteFractionModel>(-1))
{
    m_name = "Active-site Reaction";
}

// Default constructor.
ActSiteReaction::ActSiteReaction(const Sweep::Mechanism &mech,
                                 const RadicalSiteFractionModel rad_site_model,
                                 const EnvironmentInterface::PropertyIndex alpha_index)
: SurfaceReaction(mech)
, iC2H2(-1)
, iO2(-1)
, iOH(-1)
, iCO(-1)
, iH(-1)
, iH2(-1)
, iH2O(-1)
, mAlphaIndex(alpha_index)
, mRadicalSiteModel(rad_site_model)
{
    m_name = "Active-site Reaction";

    // Set up the species indices that will be needed to calculate active site density
    Sprog::SpeciesPtrVector::const_iterator i;
    int j = 0;
    for (i=mech.Species()->begin(); i!=mech.Species()->end(); ++i, ++j) {
        if  ((*i)->Name().compare("C2H2")==0) {
            iC2H2 = j;
        } else if  ((*i)->Name().compare("O2")==0) {
            iO2 = j;
        } else if  ((*i)->Name().compare("OH")==0) {
            iOH = j;
        } else if  ((*i)->Name().compare("CO")==0) {
            iCO = j;
        } else if  ((*i)->Name().compare("H")==0) {
            iH = j;
        } else if  ((*i)->Name().compare("H2")==0) {
            iH2 = j;
        } else if  ((*i)->Name().compare("H2O")==0) {
            iH2O = j;
        }
    }

    // Should really check that all the indices are now >= 0
}

// Stream-reading constructor.
ActSiteReaction::ActSiteReaction(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

/*!
 *@param[in]            t           Time at which rate is being calculated
 *@param[in]            sys         System for which rate is to be calculated
 *@param[in]            local_geom  Spatial configuration information (ignored)
 *
 *@return   Process rate
 */
double ActSiteReaction::Rate(double t, const Cell &sys,
                           const Geometry::LocalGeometry1d &local_geom) const
{
    return SurfaceReaction::Rate(t, sys, local_geom) * SiteDensity(sys);
}


// SINGLE PARTICLE RATE CALCULATIONS.

// Returns the rate of the process for the given particle in
// the system. Process must be linear in particle number.
double ActSiteReaction::Rate(double t, const Cell &sys, const Particle &sp) const
{
    return SurfaceReaction::Rate(t, sys, sp) * SiteDensity(sys);
}

/*!
 *  This calculation is part of the ABF Hydrogen Abstraction - C2H2 Addtion (HACA) model
 *  for soot particles as discussed by Appel et al., Proc. Combust. Inst. (2000).
 *
 *\param[in]    gas     Gas mixture containing the soot particles
 *
 *\return       Fraction of surface sites which are radicals.
 *
 */
double ActSiteReaction::radicalSiteFractionABF(const EnvironmentInterface &gas) const
{
    double r1f, r1b, r2f, r2b, r3f, r4f, r5f, rdenom;
    double T  = gas.Temperature();
    double RT = RCAL * T;

    // Calculate the forward and back reaction rates.
    r1f = 4.2e+07 * exp(-13.0/RT)                * gas.SpeciesConcentration(iH);
    r1b = 3.9e+06 * exp(-11.0/RT)                * gas.SpeciesConcentration(iH2);
    r2f = 1.0e+04 * exp(-1.43/RT) * pow(T,0.734) * gas.SpeciesConcentration(iOH);
    r2b = 3.68e+2 * exp(-17.1/RT) * pow(T,1.139) * gas.SpeciesConcentration(iH2O);
    r3f = 2.0e+07                                * gas.SpeciesConcentration(iH);
    r4f = 8.0e+01 * exp( -3.8/RT) * pow(T,1.56)  * gas.SpeciesConcentration(iC2H2);
    r5f = 2.1e+06 * exp( -7.47/RT)               * gas.SpeciesConcentration(iO2);
    rdenom = r1b+r2b+r3f+r4f+r5f;

    if (rdenom > 0.0) {
        double f = (r1f+r2f) / rdenom;
        return f / (f + 1.0);
    } else {
        return 0.0;
    }
}

/*!
 *  This calculation is an alternative to the ABF Hydrogen Abstraction - C2H2 Addtion (HACA) model
 *  for soot particles.  Constants are taken from Blanquart & Pitsch, Combustion & Flame 156 (2009)
 *  1614-1626.
 *
 *\param[in]    gas     Gas mixture containing the soot particles
 *
 *\return       Fraction of surface hydrogens which are radicals.
 *
 */
double ActSiteReaction::radicalSiteFractionBP(const EnvironmentInterface &gas) const
{
    const double T  = gas.Temperature();
    const double RT = R * T / 1000;

    // Calculate the forward and back reaction rates.
    const double r1f = 1.0e+08  * std::exp(-68.42/RT)  * std::pow(T,1.8)   * gas.SpeciesConcentration(iH);
    const double r1b = 8.68e+04 * std::exp(-25.46/RT)  * std::pow(T,2.36)  * gas.SpeciesConcentration(iH2);
    const double r2f = 6.72e+01 * std::exp(-6.09/RT)   * std::pow(T,3.33)  * gas.SpeciesConcentration(iOH);
    const double r2b = 6.44e-01 * std::exp(-27.96/RT)  * std::pow(T,3.79)  * gas.SpeciesConcentration(iH2O);
    const double r3f = 1.13e+16 * std::exp(-476.05/RT) * std::pow(T,-0.06);
    const double r3b = 4.17e+13 *                        std::pow(T,0.15)  * gas.SpeciesConcentration(iH);
    const double r4f = 2.52e+09 * std::exp(-17.13/RT)  * std::pow(T,1.10)  * gas.SpeciesConcentration(iC2H2);

    // This is a quasi-steady state calculation, the number of radical active sites is
    // calculated by assuming that there is a steady state between radical and non-radical
    // surface hydrogens.
    const double radicalDepletion = r1b + r2b + r3b + r4f;
    const double radicalProduction  = r1f + r2f + r3f;

    if (radicalProduction > 0.0) {
        return radicalProduction / (radicalProduction + radicalDepletion);
    } else {
        return 0.0;
    }
}

/*!
 *  This calculation is part of the ABF Hydrogen Abstraction - C2H2 Addtion (HACA) model
 *  for soot particles as discussed by Appel et al., Proc. Combust. Inst. (2000) or of
 *  the alternative taken from  Blanquart & Pitsch, Combustion & Flame 156 (2009)
 *  1614-1626, which requires the modelling of the number of surface hydrogen atoms
 *  on a particle (surf vol hydrogen model).
 *
 *\param[in]    sys     System for which rate is to be calculated
 *
 *\return       Concentration of surface sites available for reaction
 *
 */
double ActSiteReaction::SiteDensity(const Cell &sys) const
{
    double alpha(0.0);

    if (!sys.FixedChem())
        alpha = sys.Particles().Alpha(sys.GasPhase().Temperature());
    else
        alpha = sys.GasPhase().PropertyValue(mAlphaIndex);

    switch(mRadicalSiteModel) {
    case ABFRadicalSiteModel:
        // This value of 2.3e19m^-2 is the density of H sites (whether active or not)
        // on the surface of a soot particle.  It can be found on p179 of the article
        // by Frenklach & Wang in "Soot Formation in Combustion", ed Bockhorn,
        // pub Springer 1994.
        return 2.3e19 * radicalSiteFractionABF(sys.GasPhase()) * alpha;
    case BPRadicalSiteModel:
        return radicalSiteFractionBP(sys.GasPhase());
    case ABFConstant:
        // Same as ABFRadicalSiteModel except here alpha = 1.0.
        return 2.3e19 * radicalSiteFractionABF(sys.GasPhase());
    default:
        throw std::runtime_error("Unrecognised radical site model in ActSiteReaction::radicalSiteFraction");
        return 0.0;
    }
}

// READ/WRITE/COPY.

// Creates a copy of the particle process.
ActSiteReaction *const ActSiteReaction::Clone(void) const
{
    return new ActSiteReaction(*this);
}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType ActSiteReaction::ID(void) const {return ActSiteRxn_ID;}

// It is not obvious why the serialisation/deserialisation does
// not include the species index members. (Rob P 05.04.2012)

// Writes the object to a binary stream.
void ActSiteReaction::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        SurfaceReaction::Serialize(out);

        out.write(reinterpret_cast<const char*>(&mAlphaIndex), sizeof(mAlphaIndex));
        out.write(reinterpret_cast<const char*>(&mRadicalSiteModel), sizeof(mRadicalSiteModel));

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ActSiteReaction::Serialize).");
    }
}

// Reads the object from a binary stream.
void ActSiteReaction::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                // Deserialize base class.
                SurfaceReaction::Deserialize(in, mech);

                in.read(reinterpret_cast<char*>(&mAlphaIndex), sizeof(mAlphaIndex));
                in.read(reinterpret_cast<char*>(&mRadicalSiteModel), sizeof(mRadicalSiteModel));

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ActSiteReaction::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ActSiteReaction::Deserialize).");
    }
}
