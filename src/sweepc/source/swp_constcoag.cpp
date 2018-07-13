/*
  Author(s):      Robert I A Patterson
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2011 Robert I A Patterson.

  File purpose:
    Implementation of constant additive coagulation kernel

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
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

#include "swp_constcoag.h"

#include "swp_params.h"
#include "swp_cell.h"
#include "swp_mechanism.h"
#include <boost/random/uniform_01.hpp>

using namespace Sweep::Processes;

const double Sweep::Processes::ConstantCoagulation::s_MajorantFactor = 1.5;

/**
 * Main way of building a new coagulation object
 * @param[in] mech      Mechanism to which coagulation will belong
 *
 */
Sweep::Processes::ConstantCoagulation::ConstantCoagulation(const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "ConstantCoagulation";
}

// Stream-reading constructor.
Sweep::Processes::ConstantCoagulation::ConstantCoagulation(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "ConstantCoagulation";
    Deserialize(in, mech);
}

// Returns the rate of the process for the given system.
double Sweep::Processes::ConstantCoagulation::Rate(double t, const Cell &sys,
                                                        const Geometry::LocalGeometry1d &local_geom) const
{
    // Get the number of particles in the system.
    unsigned int n = sys.ParticleCount(); // aab64 Removed const type to be compatible with the below. 

	// aab64 for hybrid particle model
	if (m_mech->IsHybrid())
		n += (unsigned int)sys.GetIncepted();

    return A() * n * (n - 1) * s_MajorantFactor / sys.SampleVolume() / 2;
}

/**
 * Number of terms in the expression for the sum of the majorant
 * kernel over all particle pairs.
 */
unsigned int Sweep::Processes::ConstantCoagulation::TermCount() const {return TYPE_COUNT;}


/**
 * Calculate the terms in the sum of the majorant kernel over all particle
 * pairs, placing each term in successive positions of the sequence
 * beginning at iterm and return the sum of the terms added to that
 * vector.
 *
 * @param[in]     t            Time for which rates are requested
 * @param[in]     sys          Details of the particle population and environment
 * @param[in]     local_geom   Position information
 * @param[in,out] iterm        Pointer to start of sequence to hold the rate terms, returned as one past the end.
 */
double Sweep::Processes::ConstantCoagulation::RateTerms(double t, const Cell &sys,
                                                             const Geometry::LocalGeometry1d &local_geom,
                                                             fvector::iterator &iterm) const
{
    return *(iterm++) = Rate(t, sys, local_geom);
}

/*!
 * 
 *
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       local_geom  Details of local phsyical layout
 * \param[in]       iterm       Process term responsible for this event
 * \param[in,out]   rng         Random number generator
 *
 * \return      0 on success, otherwise negative.
 */
int ConstantCoagulation::Perform(double t, Sweep::Cell &sys, 
                             const Geometry::LocalGeometry1d& local_geom,
                             unsigned int iterm,
                             Sweep::rng_type &rng) const
{
    // Select properties by which to choose particles.
    // Note we need to choose 2 particles.  One particle must be chosen
    // uniformly and one with probability proportional
    // to particle mass.
	int ip1 = -1, ip2 = -1;

	// aab64 hybrid particle model flags
	bool hybrid_flag = m_mech->IsHybrid();
	bool ip1_flag = false;
	bool ip2_flag = false;

	double n_incep = sys.GetIncepted();
	double n_other = sys.ParticleCount();
	double n_total = n_incep + n_other;

	boost::uniform_01<rng_type&, double> unifDistrib(rng);
	double frac = 0;
	if (n_total > 0)
		frac = n_other / n_total;
	else
		return -1;

    if (n_total < 2)  // if there are < 2 SPs but incepting class has weight >= 2, we can still act
        return 1;
	else 
	{
		if (hybrid_flag)
		{
			// Particle 1 is picked uniformly. Here, the
			// incepting class has multiple weight 1 particles
			// Account for this by selecting this by default and 
			// switching with probability n_other/n_total
			ip1 = -2;
			ip2 = -2;
			if (frac >= unifDistrib() && frac != 0)
			{
				ip1 = sys.Particles().Select(rng);
			}
			if (frac >= unifDistrib() && frac != 0)
			{
				ip2 = sys.Particles().Select(rng);
			}
		}
		else
		{
			ip1 = sys.Particles().Select(rng);
			ip2 = sys.Particles().Select(rng);
		}
	}

	// Choose and get first particle, then update it.
	Particle *sp1 = NULL;

	// Is this an incepting class particle?
	if (hybrid_flag && ip1 == -2)
	{
		m_mech->SetRandomParticle(true, sys, t, unifDistrib(), true, sys.GetMuLN(), sys.GetSigmaLN(),
			sys.Particles().GetInceptedSP_youngest().CollDiameter(), sys.Particles().GetInceptedSP_oldest().CollDiameter(), 1.0, rng);
		
		sp1 = sys.Particles().GetInceptedSP_tmp_d_1().Clone();
		ip1_flag = true;                                                             // Flag sp1 as an incepting class particle
		sys.AdjustIncepted(-1.0);                                                    // Reduce the incepting class count
		double dinc = sp1->CollDiameter();                                           // Update moments for removal of a particle from the space
		sys.SetMomentsk(sys.GetMomentsk_0() + 1.0,
			sys.GetMomentsk_1() + dinc,
			sys.GetMomentsk_2() + dinc * dinc,
			sys.GetMomentsk_3() + dinc * dinc * dinc);

		// If incepting class is now empty, pick another particle before adding sp1 to the ensemble
		if (ip2 == -2 && sys.GetIncepted() == 0)
			ip2 = sys.Particles().Select(rng);

		ip1 = sys.Particles().Add(*sp1, rng);                                        // Add the particle to the ensemble
	}
	else
	{
		if (ip1 >= 0) {
			sp1 = sys.Particles().At(ip1);
		}
		else {
			// Failed to choose a particle.
			return -1;
		}
	}

    // Choose and get unique second particle, then update it.  Note, we are allowed to do
    // this even if the first particle was invalidated.
    unsigned int guard = 0;
	while ((ip2 == ip1) && (++guard < 1000))
	{
		if (hybrid_flag)
		{
			ip2 = -2;
			if (ip1_flag)
			{
				++n_other;
				frac = n_other / n_total;
			}
			if (frac >= unifDistrib() && frac != 0)
			{
				ip2 = sys.Particles().Select(rng);
			}
		}
		else
			ip2 = sys.Particles().Select(rng);
	}

	// Choose and get second particle, then update it.
    Particle *sp2 = NULL;
	double dsp2 = 0.0;

	// Is this an incepting class particle?
	if (hybrid_flag && ip2 == -2)
	{
		m_mech->SetRandomParticle(false, sys, t, unifDistrib(), true, sys.GetMuLN(), sys.GetSigmaLN(),
			sys.Particles().GetInceptedSP_youngest().CollDiameter(), sys.Particles().GetInceptedSP_oldest().CollDiameter(), 1.0, rng);
		// Note don't need to add it to the ensemble unless coagulation is successful
		sp2 = sys.Particles().GetInceptedSP_tmp_d_2().Clone();
		ip2_flag = true;                                                             // Flag sp2 as an incepting class particle
		dsp2 = sp2->CollDiameter();
	}
	else
	{
		if ((ip2 >= 0) && (ip2 != ip1)) {
			sp2 = sys.Particles().At(ip2);
		}
		else {
			// Failed to select a unique particle.
			return -1;
		}
	}

    //Calculate the majorant rate before updating the particles
    const double majk = MajorantKernel(*sp1, *sp2, sys, Default);
	
    //Update the particles
	m_mech->UpdateParticle(*sp1, sys, t, rng);

    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp1->IsValid()) {
        // Must remove first particle now.
		sys.Particles().Remove(ip1);
		
        // Invalidating the index tells this routine not to perform coagulation.
        ip1 = -1;
        return 0;
    }

	m_mech->UpdateParticle(*sp2, sys, t, rng);
	
	// Check validity of particles after update.
    if (!sp2->IsValid()) {
        // Tell the ensemble to update particle one before we confuse things
        // by removing particle 2
        sys.Particles().Update(ip1);

		if (!ip2_flag)
			// Must remove second particle now.
			sys.Particles().Remove(ip2);
		else
		{
			// Particle sp2 is not in the ensemble, must manually delete it
			delete sp2;
			sp2 = NULL;
		}

        // Invalidating the index tells this routine not to perform coagulation.
        ip2 = -1;

        return 0;
    }

    // Check that both the particles are still valid.
    if ((ip1 != -1) && (ip2 != -1)) {
        // Must check for ficticious event now by comparing the original
        // majorant rate and the current (after updates) true rate.
		
        double truek = CoagKernel(*sp1, *sp2, sys);
		
        if (!Fictitious(majk, truek, rng)) {
			// If particle sp2 is used, we now need to remove it from the incepting class
			if (ip2_flag)
			{
				sys.AdjustIncepted(-1.0);                                                    // Reduce the incepting class count
				sys.SetMomentsk(sys.GetMomentsk_0() + 1.0,                                   // Update moments for removal of a particle from the space
					sys.GetMomentsk_1() + dsp2,
					sys.GetMomentsk_2() + dsp2 * dsp2,
					sys.GetMomentsk_3() + dsp2 * dsp2 * dsp2);
			}
			JoinParticles(t, ip1, sp1, ip2, sp2, sys, rng);
			if (ip2_flag && sp2 != NULL)
			{
				// Particle sp2 is not in the ensemble, must manually delete it
				delete sp2;
				sp2 = NULL;
			}

        } else {
			sys.Particles().Update(ip1);
			if (!ip2_flag)
				sys.Particles().Update(ip2);
			else if (sp2 != NULL)
			{
				// Particle sp2 is not in the ensemble, must manually delete it
				delete sp2;
				sp2 = NULL;
			}
            return 1; // Ficticious event.
        }
    } else {
        // One or both particles were invalidated on update,
        // but that's not a problem.  Information on the update
        // of valid particles must be propagated into the binary
        // tree
		if (ip1 != -1)
            sys.Particles().Update(ip1);

		if (ip2 != -1 && !ip2_flag)
			sys.Particles().Update(ip2);

		if (ip2_flag && sp2 != NULL)
		{
			// Particle sp2 is not in the ensemble, must manually delete it
			delete sp2;
			sp2 = NULL;
		}
    }

	if (ip2_flag && sp2 != NULL)
	{
		// Particle sp2 is not in the ensemble, must manually delete it
		delete sp2;
		sp2 = NULL;
	}

    return 0;
}

/**
 * Calculate the coagulation kernel between two particles in a given environment.
 * Note that the details of the environment are not currently used.
 *
 *@param[in]    sp1         First particle
 *@param[in]    sp2         Second particle
 *@param[in]    sys         Details of the environment
 *
 *@return       Value of kernel
 */
double Sweep::Processes::ConstantCoagulation::CoagKernel(const Particle &sp1,
                                                              const Particle &sp2,
                                                              const Cell& sys) const
{
    return A();
}


/**
 * Calculate the majorant kernel between two particles in a given environment.
 * Note that the details of the environment are not currently used.
 *
 *@param[in]    sp1         First particle
 *@param[in]    sp2         Second particle
 *@param[in]    sys         Details of the environment
 *@param[in]    maj         Unused flag to indicate which majorant kernel is required
 *
 *@return       Value of majorant kernel
 */
double Sweep::Processes::ConstantCoagulation::MajorantKernel(const Particle &sp1,
                                                                  const Particle &sp2,
                                                                  const Cell& sys,
                                                                  const MajorantType maj) const
{
    return CoagKernel(sp1, sp2, sys) * s_MajorantFactor;
}

