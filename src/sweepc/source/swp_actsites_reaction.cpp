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
{
    m_name = "Active-site Reaction";
}

// Default constructor.
ActSiteReaction::ActSiteReaction(const Sweep::Mechanism &mech)
: SurfaceReaction(mech)
, iC2H2(-1)
, iO2(-1)
, iOH(-1)
, iCO(-1)
, iH(-1)
, iH2(-1)
, iH2O(-1)
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
real ActSiteReaction::Rate(real t, const Cell &sys,
                           const Geometry::LocalGeometry1d &local_geom) const
{
    return SurfaceReaction::Rate(t, sys, local_geom) * SiteDensity(sys.GasPhase());
}


// SINGLE PARTICLE RATE CALCULATIONS.

// Returns the rate of the process for the given particle in
// the system. Process must be linear in particle number.
real ActSiteReaction::Rate(real t, const Cell &sys, const Particle &sp) const
{
    return SurfaceReaction::Rate(t, sys, sp) * SiteDensity(sys.GasPhase());
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
real ActSiteReaction::radicalSiteFraction(const Sprog::Thermo::IdealGas &gas) const
{
    real r1f, r1b, r2f, r2b, r3f, r4f, r5f, rdenom;
    real T  = gas.Temperature();
    real RT = RCAL * T;

    // Calculate the forward and back reaction rates.
    r1f = 4.2e+07 * exp(-13.0/RT)                * gas.MolarConc(iH);
    r1b = 3.9e+06 * exp(-11.0/RT)                * gas.MolarConc(iH2);
    r2f = 1.0e+04 * exp(-1.43/RT) * pow(T,0.734) * gas.MolarConc(iOH);
    r2b = 3.68e+2 * exp(-17.1/RT) * pow(T,1.139) * gas.MolarConc(iH2O);
    r3f = 2.0e+07                                * gas.MolarConc(iH);
    r4f = 8.0e+01 * exp( -3.8/RT) * pow(T,1.56)  * gas.MolarConc(iC2H2);
    r5f = 2.1e+06 * exp( -7.47/RT)               * gas.MolarConc(iO2);
    rdenom = r1b+r2b+r3f+r4f+r5f;

    if (rdenom > 0.0) {
        real f = (r1f+r2f) / rdenom;
        return f / (f + 1.0);
    } else {
        return 0.0;
    }
}

/*!
 *  This calculation is part of the ABF Hydrogen Abstraction - C2H2 Addtion (HACA) model
 *  for soot particles as discussed by Appel et al., Proc. Combust. Inst. (2000).
 *
 *\param[in]    gas     Gas mixture containing the soot particles
 *
 *\return       Concentration of surface sites available for reaction
 *
 */
real ActSiteReaction::SiteDensity(const Sprog::Thermo::IdealGas &gas) const
{
    return 2.3e19 * radicalSiteFraction(gas) * gas.Alpha();
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

// Writes the object to a binary stream.
void ActSiteReaction::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        SurfaceReaction::Serialize(out);

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

        unsigned int n = 0, id = 0;

        switch (version) {
            case 0:
                // Deserialize base class.
                SurfaceReaction::Deserialize(in, mech);

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
