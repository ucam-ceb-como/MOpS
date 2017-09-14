/*
  Author(s):	  Shraddha Shekar (ss663)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2010 Shraddha Shekar.

  File purpose:
    Implementation of the InterParticle class declared in the
    swp_InterParticle.h header file.

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

#include "swp_silica_interparticle.h"
#include "swp_sintering_model.h"
#include "swp_particle_model.h"
#include "swp_mechanism.h"
#include "swp_process_type.h"
#include "gpc_species.h"

#include <cmath>
#include <stdexcept>
#include <cassert>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

const double InterParticle::m_majfactor = 2.0;

// CONSTRUCTORS AND DESTRUCTORS.

//! Default constructor (protected).
InterParticle::InterParticle(void)
: ParticleProcess(), m_arr(0.0,0.0,0.0)
{
    m_defer = true;
    m_name = "InterParticle";
}

/*!
 * Initialising constructor
 *
 *@param[in]    mech        Mechanism of which this process is a part
 *@param[in]
 */
InterParticle::InterParticle(const Sweep::Mechanism &mech,
                             const EnvironmentInterface::SpeciesIndex h4o4si)
: ParticleProcess(mech)
, m_arr(0.0,0.0,0.0)
, m_H4O4SI_Index(h4o4si)
{
    // Assume the InterParticle is simulated as a deferred process (LPDA).
    m_defer = true;
    m_name = "InterParticle";
}

//! Copy constructor.
InterParticle::InterParticle(const InterParticle &copy)
{
    *this = copy;
}

//! Stream-reading constructor.
InterParticle::InterParticle(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

//! Default destructor.
InterParticle::~InterParticle(void)
{
    // Nothing to destruct.
}


// OPERATOR OVERLOADS.

//! Assignment operator.
InterParticle &InterParticle::operator=(const InterParticle &rhs)
{
    if (this != &rhs) {
        ParticleProcess::operator =(rhs);
        m_arr    = rhs.m_arr;
		m_pid     = rhs.m_pid;
		m_H4O4SI_Index = rhs.m_H4O4SI_Index;
    }
    return *this;
}


// RATE CONSTANT AND PARAMETERS.

//! Returns the Arrhenius parameter.
Sprog::Kinetics::ARRHENIUS &InterParticle::Arrhenius() {return m_arr;}
const Sprog::Kinetics::ARRHENIUS &InterParticle::Arrhenius() const {return m_arr;}

//! Sets the fixed rate constant.
void InterParticle::SetArrhenius(Sprog::Kinetics::ARRHENIUS &arr) {m_arr = arr;}


// PARTICLE PROPERTY ID.

/*!
 * @brief       Returns the PropID to which the rate is proportional
 * @return      ID of particle property
 */
unsigned int InterParticle::PropertyID(void) const {return m_pid;}

/*!
 * @brief       Sets the PropID to which the rate is proportional
 * @param[in]   pid    ID of particle property
 */
void InterParticle::SetPropertyID(PropID pid)
{
    m_pid = pid;
}


// TOTAL RATE CALCULATIONS (ALL PARTICLES IN A SYSTEM).

/*!
 * @brief       Returns the rate of the process in a given system
 * 
 * The interparticle reaction rate is given by the difference between
 * the surface reaction rate and the sintering rate. This function
 * calcualtes the present value of SR rate, and gets the sintrate from 
 * the particle cache.
 * 
 * @param[in]   t             Time at which process occurs
 * @param[in]   sys           System for rate calculation
 * @param[in]   local_geom    Local geometry
 * 
 * @return      Rate of process
 */
double InterParticle::Rate(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom) const
{

	// From theoretical calculations: R_int = R_surf - sintering_contribution + R_inc
    double rate = 0.0;

	// First calculate surface reaction contribution:
	double T = sys.GasPhase().Temperature();

	// Get the total number of OH sites from cache
	double numOH = sys.Particles().GetSum(static_cast<Sweep::PropID>(m_pid));

	// Calculate the concentration of SiOH4
	// Note that the usual ParticleProcess interface can NOT be used, as this
	// would erroneously remove H4O4SI from the gas-phase.
	double conc = sys.GasPhase().SpeciesConcentration(m_H4O4SI_Index);

	// Rate of surface reaction
	double R_surf = m_arr.A * conc *pow(T, m_arr.n)
			* exp(-m_arr.E / (R * T)) * numOH;

	// Forward-declare the total sintering rate
	double total_sint_rate = 0;

	// Get total surface area from cache
	double surface = sys.Particles().GetSum(static_cast<Sweep::PropID>(Sweep::iS));

	// Check if particle surface area exists. If not, set the sintering rate to zero.
	if (surface == 0) {
		surface = 1;
		total_sint_rate = 0;
	} else {
		// Get total sintering rate from cache.
		total_sint_rate = sys.Particles().GetSum(static_cast<Sweep::PropID>(Sweep::iSintRate));
	}

	// Calculate the total surface density of sites
	double rho_s =  numOH/surface;
	//Rate of sintering
	double R_sint = (rho_s*total_sint_rate)/(2.0);

	rate += (R_surf - R_sint);

	if(rate < 0)
		rate = 0;

	if (m_mech->AnyDeferred())
	{
        return rate * m_majfactor;
    }
	else
	{
        return rate;
    }
}



// SINGLE PARTICLE RATE CALCULATIONS.

/*!
 * @brief       Returns the rate of the process for a particle
 * 
 * The interparticle reaction rate is given by the difference between
 * the surface reaction rate and the sintering rate. This function
 * calculates the present value of SR rate, and gets the sintrate from
 * the particle cache.
 * 
 * @param[in]   t       Time at which process occurs
 * @param[in]   sys     System for rate calculation
 * @param[in]   sp      Particle for rate calculation
 * 
 * @return      Rate of process
 */
double InterParticle::Rate(double t, const Cell &sys, const Particle &sp) const
{

	// Rate constant.
    double rate = m_arr.A;

    // Calculate the concentration of SiOH4
    // Note that the usual ParticleProcess interface can NOT be used, as this
    // would erroneously remove H4O4SI from the gas-phase.
    double conc = sys.GasPhase().SpeciesConcentration(m_H4O4SI_Index);

    // Do calculation for surface-reaction part of rate:
    // Chemical species concentration dependence.
    rate *= conc;

    // Temperature dependance.
    double T = sys.GasPhase().Temperature();
    rate *= pow(T, m_arr.n) * exp(-m_arr.E / (R * T));

	// Get the number of OH sites from cache
	double numOH = sp.Property(static_cast<Sweep::PropID>(m_pid));
	rate *= numOH;

	// Do calculation for sintering part of rate:
    // Forward-declare some parameters
    double sint_rate = 0.0;
    double rho_s = 0.0;

	// Get surface area from cache
    double surface = sp.Property(static_cast<Sweep::PropID>(Sweep::iS));

	if(surface == 0) {
		rho_s = 1;
		sint_rate = 0;
	} else {
		rho_s = numOH / surface;
		//Calculate sintering contribution
		sint_rate = (sp.Property(static_cast<Sweep::PropID>(Sweep::iSintRate)));
	}


	double sint_comp = (sint_rate*rho_s)/(2.0);

	// Sintering rate dependency
	rate -= sint_comp;

	if(rate<0)
		rate = 0;

    return rate; //trm[0] + trm[1] + trm[2];
}

/*!
 * @brief       Returns the majorant rate.
 * 
 * @param[in]   t       Time at which process occurs
 * @param[in]   sys     System for rate calculation
 * @param[in]   sp      Particle for majorant rate calculation
 * 
 * @return      Rate of process
 */
double InterParticle::MajorantRate(double t, const Cell &sys, const Particle &sp) const
{
    // Return the single particle rate multiplied by the
    // majorant factor.
    return Rate(t, sys, sp) * m_majfactor;
}


// RATE TERM CALCULATIONS.
//   These routines return the individual rate terms for a
//   process, which may have multiple terms (e.g. InterParticle).

//! Returns the number of rate terms for this process.
unsigned int InterParticle::TermCount(void) const {return 1;}

/*!
 * @brief       Passes the system rate to an iterator
 * 
 * Calculates the rate terms given an iterator to a double vector. The
 * iterator is advanced to the position after the last term for this
 * process.
 * 
 * @return      Rate of process
 */
double InterParticle::RateTerms(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom,
                             fvector::iterator &iterm) const
{
    return *(iterm++) = Rate(t, sys, local_geom);
}

// PERFORMING THE PROCESS.

/*!
 * @brief       Performs the interparticle process system-wide
 * 
 * Selects a particle, calculates the rate, and updates the system
 * accordingly. Note that the gas-phase is actually adjusted in the
 * Sinter() function of SilicaPrimary.
 *
 * @param[in]       t           Time
 * @param[in,out]   sys         System to update
 * @param[in]       local_geom  Details of local physical layout
 * @param[in]       iterm       Process term responsible for this event
 * @param[in,out]   rng         Random number generator
 *
 * @return      0 on success, otherwise negative.
 */
int InterParticle::Perform(double t, Sweep::Cell &sys,
                             const Geometry::LocalGeometry1d& local_geom,
                             unsigned int iterm,
                             rng_type &rng) const
{
	PartPtrVector dummy;

    int i = sys.Particles().Select(static_cast<Sweep::PropID>(m_pid), rng);

    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);


        // Update particle with deferred processes.
        if (m_mech->AnyDeferred()) {
            // Calculate majorant rate then update the particle.
            double majr = MajorantRate(t, sys, *sp);
            m_mech->UpdateParticle(*sp, sys, t, i, rng, dummy);

            // Check that the particle is still valid.
            if (sp->IsValid()) {
                double truer = Rate(t, sys, *sp);

                if (!Fictitious(majr, truer, rng)) {
                    // Adjust particle.
                    sp->AdjustIntPar(m_dcomp, m_dvals, rng, 1);
                    sys.Particles().Update(i);

                    // Apply changes to gas-phase chemistry.
                    adjustGas(sys, sp->getStatisticalWeight());
                }
            } else {
                // If not valid then remove the particle.
                sys.Particles().Remove(i);
            }
        } else {
            // No particle update required, just perform the surface
            // reaction.
            sp->AdjustIntPar(m_dcomp, m_dvals, rng, 1);

            if (sp->IsValid()) {
                // Tell the binary tree to recalculate
                sys.Particles().Update(i);
            }
            else {
                // Particle has been removed due to oxidation
                sys.Particles().Remove(i);
            }

            // Apply changes to gas-phase chemistry.
            adjustGas(sys, sp->getStatisticalWeight());
        }
    } else {
        // Failed to select a particle.
        return -1;
    }

    return 0;
}

/*!
 * @brief       Performs process on a given particle n times
 * 
 * @param[in]   t   Time for process
 * @param[in]   sys System in which to act
 * @param[in]   sp  Particle to adjust
 * @param[in]   rng Random number generator
 * @param[in]   n   Number of times to do process
 */
int InterParticle::Perform(double t, Cell &sys, Particle &sp, rng_type &rng,
                          unsigned int n) const
{
    unsigned int m = sp.AdjustIntPar(m_dcomp, m_dvals, rng, n);
    adjustGas(sys, sp.getStatisticalWeight(), m);
    return m;
}

//! Returns the process type.
ProcessType InterParticle::ID(void) const {return InterParticle_ID;}

//! Creates a copy of the particle process.
InterParticle *const InterParticle::Clone(void) const
{
    return new InterParticle(*this);
}

//! Writes the object to a binary stream.
void InterParticle::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        ParticleProcess::Serialize(out);

        // Write arrhenius coefficients.
        double A  = (double)m_arr.A;
        double nn = (double)m_arr.n;
        double E  = (double)m_arr.E;
        out.write((char*)&A, sizeof(A));
        out.write((char*)&nn, sizeof(nn));
        out.write((char*)&E, sizeof(E));

        // Write particle property ID.
        unsigned int n = (unsigned int)m_pid;
        out.write((char*)&n, sizeof(n));

        out.write(reinterpret_cast<const char*>(&m_H4O4SI_Index), sizeof(m_H4O4SI_Index));

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, InterParticle::Serialize).");
    }
}

//! Reads the object from a binary stream.
void InterParticle::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double A = 0.0, nn = 0.0, E = 0.0;

        switch (version) {
            case 0:
                // Deserialize base class.
                ParticleProcess::Deserialize(in, mech);

                // Read arrhenius coefficients.
                in.read(reinterpret_cast<char*>(&A), sizeof(A));
                in.read(reinterpret_cast<char*>(&nn), sizeof(nn));
                in.read(reinterpret_cast<char*>(&E), sizeof(E));
                m_arr.A = (double)A;
                m_arr.n = (double)nn;
                m_arr.E = (double)E;

                // Read particle property ID.
                in.read(reinterpret_cast<char*>(&m_pid), sizeof(m_pid));

                in.read(reinterpret_cast<char*>(&m_H4O4SI_Index), sizeof(m_H4O4SI_Index));

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, InterParticle::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, InterParticle::Deserialize).");
    }
}
