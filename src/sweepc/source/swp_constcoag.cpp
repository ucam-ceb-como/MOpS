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
    //const unsigned int n = sys.ParticleCount();

	// aab64 temp
	bool hybrid_flag = true;    
	unsigned int n = sys.ParticleCount();
	unsigned int n2 = n;
	if (hybrid_flag)
		n2 = n2 + sys.GetIncepted();
	if (n > 1)
		n2 -= 1;

	//return A() * (n2 * n2 - n) * s_MajorantFactor / sys.SampleVolume() / 2;
    return A() * n2 * (n2 - 1) * s_MajorantFactor / sys.SampleVolume() / 2;
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

	// aab64 temp
	bool hybrid_flag = true;
	bool ip1_flag = false;
	bool ip2_flag = false;
	boost::uniform_01<rng_type&, double> unifDistrib(rng);
	unsigned int n_others = 0;
	unsigned int n_total  = sys.GetIncepted();
	if (sys.ParticleCount() > 1)
	{
		n_others += sys.ParticleCount() - 1;
		n_total  += n_others;
	}
	double frac = 0;
	if (n_total > 0)
		frac = (double(n_others)) / (double(n_total));

	int ip1 = -1, ip2 = -1;

    if (sys.ParticleCount() < 2) 
	{
		if (!hybrid_flag || sys.GetIncepted() < 2)
            return 1;
		else
		{
			ip1 = 0;
			ip2 = 0; // should not be required
		}
    }
	else 
	{
		ip1 = 0;
		if (frac > unifDistrib())
		{
			while (ip1 == 0)
				ip1 = sys.Particles().Select(rng); // aab64 temporary
		}

		ip2 = 0;
		if (frac > unifDistrib())
		{
			while (ip2 == 0)
				ip2 = sys.Particles().Select(rng); // aab64 temporary
		}

		//ip1 = sys.Particles().Select(iW, rng);
		//ip2 = sys.Particles().Select(iW, rng);
	}

	// Choose and get first particle, then update it.
	Particle *sp1 = NULL;

	if (hybrid_flag && ip1 == 0)
	{
		// incept particle
		//ip1 = sys.ParticleCount();
		sp1 = m_mech->Inceptions()[0]->Perform_incepted(t, sys, local_geom, 0, rng);
		ip1 = sys.Particles().Add(*sp1, rng); // must add here to avoid contraction removing chosen sp2
		sys.AdjustIncepted(-1);
		sys.Particles().At(0)->setStatisticalWeight(sys.GetIncepted());
		m_mech->UpdateParticle(*(sys.Particles().At(0)), sys, t, rng);
		ip1_flag = true;
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
		//if (!(hybrid_flag && ip1 == 0))
		//ip2 = sys.Particles().Select(iW, rng);
		ip2 = 0;
		if (frac > unifDistrib())
		{
			while (ip2 == 0)
				ip2 = sys.Particles().Select(rng); // aab64 temporary
		}

	}

    Particle *sp2 = NULL;

	if (hybrid_flag && ip2 == 0)
	{
		// incept particle
		sp2 = m_mech->Inceptions()[0]->Perform_incepted(t, sys, local_geom, 0, rng);
		ip2_flag = true;
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
			delete sp2;
			sp2 = NULL;
		}

        // Invalidating the index tells this routine not to perform coagulation.
        ip2 = -1;

        return 0;
    }

    // Check that both the particles are still valid.
    if ((ip1>=0) && (ip2>=0)) {
        // Must check for ficticious event now by comparing the original
        // majorant rate and the current (after updates) true rate.
		
        double truek = CoagKernel(*sp1, *sp2, sys);
		
        if (!Fictitious(majk, truek, rng)) {
			// Add particle to system's ensemble.
			//if (ip1_flag)
			//{
			//	ip1 = sys.Particles().Add(*sp1, rng);
			//}

			if (ip2_flag)
			{
				sys.AdjustIncepted(-1);
				sys.Particles().At(0)->setStatisticalWeight(sys.GetIncepted());
				sys.Particles().At(0)->UpdateCache();
				m_mech->UpdateParticle(*(sys.Particles().At(0)), sys, t, rng);
			}
            JoinParticles(t, ip1, sp1, ip2, sp2, sys, rng);

        } else {
			sys.Particles().Update(ip1);
			if (!ip2_flag)
				sys.Particles().Update(ip2);
            return 1; // Ficticious event.
        }
    } else {
        // One or both particles were invalidated on update,
        // but that's not a problem.  Information on the update
        // of valid particles must be propagated into the binary
        // tree
        if(ip1 >= 0)
            sys.Particles().Update(ip1);

        if(ip2 >= 0 && !ip2_flag)
            sys.Particles().Update(ip2);
    }

	if (ip2_flag)
	{
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

