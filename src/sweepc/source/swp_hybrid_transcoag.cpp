/*
  Author(s):      Robert I A Patterson
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2009 Robert I A Patterson.

  File purpose:
    Implementation of transition regime coagulation kernel
	for hybrid particle-number/particle model

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
#include "swp_hybrid_transcoag.h"

#include "swp_params.h"
#include "swp_cell.h"
#include "swp_mechanism.h"
#include "swp_PAH_primary.h"
#include <boost/random/uniform_01.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Sweep::Processes;


// Default constructor.
Sweep::Processes::HybridTransitionCoagulation::HybridTransitionCoagulation(const Sweep::Mechanism &mech)
: Coagulation(mech), m_efm(mech.GetEnhancementFM())
{
    m_name = "HybridTransitionRegimeCoagulation";
}

Sweep::Processes::HybridTransitionCoagulation* const Sweep::Processes::HybridTransitionCoagulation::Clone() const
{
	return new HybridTransitionCoagulation(*this);
}

// Stream-reading constructor.
Sweep::Processes::HybridTransitionCoagulation::HybridTransitionCoagulation(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech), m_efm(mech.GetEnhancementFM())
{
    m_name = "HybridTransitionRegimeCoagulation";
    Deserialize(in, mech);
}

// TOTAL RATE CALCULATION.

// Returns the rate of the process for the given system.
double Sweep::Processes::HybridTransitionCoagulation::Rate(double t, const Cell &sys,
                                                          const Geometry::LocalGeometry1d &local_geom) const
{
    // Get the number of particles in the system.
    unsigned int n = sys.ParticleCount() + sys.Particles().GetTotalParticleNumber();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        double T = sys.GasPhase().Temperature();
        double P = sys.GasPhase().Pressure();

        // Get the particle-number contributions
        fvector props(6, 0);
        props[0] = sys.Particles().GetTotalDiameter();
        props[1] = sys.Particles().GetTotalDiameter2();
        props[2] = sys.Particles().GetTotalDiameter_1();
        props[3] = sys.Particles().GetTotalDiameter_2();
        props[4] = sys.Particles().GetTotalMass_1_2();
        props[5] = sys.Particles().GetTotalDiameter2_mass_1_2();
	
        // Calculate the rate.
        double rate = Rate(sys.Particles().GetSums(), (double)n, sqrt(T),
                    T/sys.GasPhase().Viscosity(), MeanFreePathAir(T,P),
                    sys.SampleVolume(), props);

        props.clear();
        fvector().swap(props);

        return rate;
    } else {
        return 0.0;
    }
}

// More efficient rate routine for coagulation only.
// All parameters required to calculate rate passed
// as arguments.
double Sweep::Processes::HybridTransitionCoagulation::Rate(const Ensemble::particle_cache_type &data, double n, double sqrtT,
                       double T_mu, double MFP, double vol, fvector & props) const
{
    // Some prerequisites.
    double n_1 = n - 1.0;
    double a = CSF * T_mu * A();
    double b = a * MFP * 1.257 * A();
    double c = CFMMAJ * m_efm * CFM * sqrtT * A();

    // Summed particle properties required for coagulation rate.
    double d       = data.Property(Sweep::iDcol);
    double d2      = data.Property(Sweep::iD2);
    double d_1     = data.Property(Sweep::iD_1);
    double d_2     = data.Property(Sweep::iD_2);
    double m_1_2   = data.Property(Sweep::iM_1_2);
    double d2m_1_2 = data.Property(Sweep::iD2_M_1_2);

    // Add particle-number contributions
    d += props[0];
    d2 += props[1];
    d_1 += props[2];
    d_2 += props[3];
    m_1_2 += props[4];
    d2m_1_2 += props[5];

    // Get individual terms.
    double terms[TYPE_COUNT];
    // Slip-flow.
    terms[0] = n * n_1 * a / vol;
    terms[1] = ((d * d_1) - n) * a / vol;
    terms[2] = d_1 * n_1 * b / vol;
    terms[3] = ((d * d_2) - d_1) * b / vol;
    // Free-molecular.
    terms[4] = n_1 * d2m_1_2  * c / vol;
    terms[5] = (m_1_2 * d2 - d2m_1_2) * c / vol;

    // Sum up total coagulation rates for different regimes.
    double sf = terms[0] + terms[1] + terms[2] + terms[3];
    double fm = terms[4] + terms[5];

    if ((sf>0.0) || (fm>0.0)) {
        // There is some coagulation.
        if (sf > fm) {
            // Use free-mol majorant.
            return fm;
        } else {
            // Use slip-flow majorant.
            return sf;
        }
    }
    return 0.0;
}


// RATE TERM CALCULATION.

// Returns the number of rate terms for this process.
unsigned int Sweep::Processes::HybridTransitionCoagulation::TermCount(void) const { return TYPE_COUNT; }

/**
 * Calculate the terms in the sum of the majorant kernel over all particle
 * pairs, placing each term in successive positions of the sequence
 * beginning at iterm and return the sum of the terms added to that
 * vector.
 *
 * @param[in]     t          Time for which rates are requested
 * @param[in]     sys        Details of the particle population and environment
 * @param[in]     local_geom Details of spatial position and boundaries
 * @param[in,out] iterm      Pointer to start of sequence to hold the rate terms, returned as one past the end.
 *
 * @return      Sum of all rate terms for this process
 */
double Sweep::Processes::HybridTransitionCoagulation::RateTerms(double t, const Cell &sys,
                            const Geometry::LocalGeometry1d &local_geom,
                            fvector::iterator &iterm) const
{
    // Get the number of particles in the system.
    unsigned int n = sys.ParticleCount() + sys.Particles().GetTotalParticleNumber();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        double T = sys.GasPhase().Temperature();
        double P = sys.GasPhase().Pressure();

        // Add contributions from particle-number list
        fvector props(6, 0);
        props[0] = sys.Particles().GetTotalDiameter();
        props[1] = sys.Particles().GetTotalDiameter2();
        props[2] = sys.Particles().GetTotalDiameter_1();
        props[3] = sys.Particles().GetTotalDiameter_2();
        props[4] = sys.Particles().GetTotalMass_1_2();
        props[5] = sys.Particles().GetTotalDiameter2_mass_1_2();
		
        // Calculate the rate terms.
        double rate = RateTerms(sys.Particles().GetSums(), (double)n, sqrt(T), T/sys.GasPhase().Viscosity(),
                         MeanFreePathAir(T,P), sys.SampleVolume(), iterm, props);
        props.clear();
        fvector().swap(props);
        return rate;
    } else {
        // No coagulation as there are too few particles.
        for(unsigned int i=0; i<TYPE_COUNT; i++) *(iterm++) = 0.0;
        return 0.0;
    }
}

// More efficient rate routine for coagulation only.
// All parameters required to calculate rate terms
// passed as arguments.
double Sweep::Processes::HybridTransitionCoagulation::RateTerms(const Ensemble::particle_cache_type &data, double n, double sqrtT,
                            double T_mu, double MFP, double vol,
                            fvector::iterator &iterm, fvector & props) const
{
    // Some prerequisites.
    double n_1 = n - 1.0;
    double a   = CSF * T_mu * A();
    double b   = a * MFP * 1.257 * 2.0;
    double c   = CFMMAJ * m_efm * CFM * sqrtT * A();

    // Summed particle properties required for coagulation rate.
    double d       = data.Property(Sweep::iDcol);
    double d2      = data.Property(Sweep::iD2);
    double d_1     = data.Property(Sweep::iD_1);
    double d_2     = data.Property(Sweep::iD_2);
    double m_1_2   = data.Property(Sweep::iM_1_2);
    double d2m_1_2 = data.Property(Sweep::iD2_M_1_2);

    // Add particle-number contributions
    d       += props[0];
    d2      += props[1];
    d_1     += props[2];
    d_2     += props[3];
    m_1_2   += props[4];
    d2m_1_2 += props[5];

    fvector::iterator isf = iterm;
    fvector::iterator ifm = iterm+4;

    // Slip-flow.
    *(iterm)   = n * n_1 * a / vol;
    *(++iterm) = ((d * d_1) - n) * a / vol;
    *(++iterm) = d_1 * n_1 * b / vol;
    *(++iterm) = ((d * d_2) - d_1) * b / vol;
    // Free-molecular.
    *(++iterm) = n_1 * d2m_1_2  * c / vol;
    *(++iterm) = (m_1_2 * d2 - d2m_1_2) * c / vol;

    // Return iterator to next term after the coagulation terms.
    ++iterm;

    // Sum up total coagulation rates for different regimes.
    double sf = *(isf) + *(isf+1) + *(isf+2) + *(isf+3);
    double fm = *(ifm) + *(ifm+1);

    if ((sf>0.0) || (fm>0.0)) {
        // There is some coagulation.
        if (sf > fm) {
            // Use free-mol majorant.
            *(isf) = 0.0;
            *(isf+1) = 0.0;
            *(isf+2) = 0.0;
            *(isf+3) = 0.0;
            return fm;
        } else {
            // Use slip-flow majorant.
            *(ifm) = 0.0;
            *(ifm+1) = 0.0;
            return sf;
        }
    } else {
        // Something went wrong with the rate calculation.
        *(isf)   = 0.0;
        *(isf+1) = 0.0;
        *(isf+2) = 0.0;
        *(isf+3) = 0.0;
        *(ifm)   = 0.0;
        *(ifm+1) = 0.0;
        return 0.0;
    }
}


// PERFORMING THE PROCESS.

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
int HybridTransitionCoagulation::Perform(double t, Sweep::Cell &sys, 
                                   const Geometry::LocalGeometry1d& local_geom,
                                   unsigned int iterm,
                                   Sweep::rng_type &rng) const
{
	PartPtrVector dummy;

    int ip1=-1, ip2=-1;
    MajorantType maj;
    TermType term = (TermType)iterm;

    // Store number of particles in each particle model
    double n_incep = sys.Particles().GetTotalParticleNumber();
    double n_other = sys.ParticleCount();
    double n_total = n_incep + n_other;
    unsigned int index1 = 0, index2 = 0, n_index1 = 0; 

    // Hybrid particle model flags
    bool hybrid_flag = m_mech->IsHybrid() && n_incep > 0.0; // Check if particle-number list is relevant
    bool ip1_flag = false;      // Flag particle 1 as a particle-number particle
    bool ip2_flag = false;      // Flag particle 2 as a particle-number particle
    bool coag_in_place = false; // Flag for coagulation in particle-number list

    // Choose and get first particle.
    Particle *sp1 = NULL;
    Particle *sp2 = NULL;

    // Hybrid model: if < 2 SPs but incepting class >= 2, we can still act
    if (n_total < 2)
        return 1;
    else
    {
        // Select properties by which to choose particles. 
        // Note we need to choose 2 particles.  There are six possible 
        // rate terms to choose from; 4 slip-flow and 2 free molecular.
        boost::uniform_01<rng_type&, double> unifDistrib(rng);
	
        // Get property sums for number list	
        double dc_incep = sys.Particles().GetTotalDiameter();
        double dc2_incep = sys.Particles().GetTotalDiameter2();
        double dc_1_incep = sys.Particles().GetTotalDiameter_1();
        double dc_2_incep = sys.Particles().GetTotalDiameter_2();
        double m_1_2_incep = sys.Particles().GetTotalMass_1_2();
        double dc2_m_1_2_incep = sys.Particles().GetTotalDiameter2_mass_1_2();

        // Get property sums for ensemble
        double dc_other = sys.Particles().GetSum(iDcol);
        double dc2_other = sys.Particles().GetSum(iD2);
        double dc_1_other = sys.Particles().GetSum(iD_1);
        double dc_2_other = sys.Particles().GetSum(iD_2);
        double m_1_2_other = sys.Particles().GetSum(iM_1_2);
        double dc2_m_1_2_other = sys.Particles().GetSum(iD2_M_1_2);

        // Select the first particle and note the majorant type
        double alpha1 = 0.0, alpha2 = 0.0;
        switch (term) {
            case SlipFlow1:
                alpha1 = unifDistrib() * n_total;
                if (alpha1 <= n_incep)
                {
                    index1 = m_mech->SetRandomParticle(sys.Particles(), t, alpha1, iUniform, rng);
                    ip1 = -2;
                } else {
                    ip1 = sys.Particles().Select_usingGivenRand(iUniform, alpha1 - n_incep, rng);
                }
                maj = SlipFlow;
                break;
            case SlipFlow2:
                alpha1 = unifDistrib() * (dc_incep + dc_other);
                if (alpha1 <= dc_incep)
                {
                    index1 = m_mech->SetRandomParticle(sys.Particles(), t, alpha1, iDcol, rng);
                    ip1 = -2;
                } else {
                    ip1 = sys.Particles().Select_usingGivenRand(iDcol, alpha1 - dc_incep, rng);
                }
                maj = SlipFlow;
                break;
            case SlipFlow3:
                alpha1 = unifDistrib() * n_total;
                if (alpha1 <= n_incep)
                {
                    index1 = m_mech->SetRandomParticle(sys.Particles(), t, alpha1, iUniform, rng);
                    ip1 = -2;
                } else {
                    ip1 = sys.Particles().Select_usingGivenRand(iUniform, alpha1 - n_incep, rng);
                }
                maj = SlipFlow;
                break;
            case SlipFlow4:
                alpha1 = unifDistrib() * (dc_incep + dc_other);
                if (alpha1 <= dc_incep)
                {
                    index1 = m_mech->SetRandomParticle(sys.Particles(), t, alpha1, iDcol, rng);
                    ip1 = -2;
                } else {
                    ip1 = sys.Particles().Select_usingGivenRand(iDcol, alpha1 - dc_incep, rng);
                }
                maj = SlipFlow;
                break;
            case FreeMol1:
                alpha1 = unifDistrib() * n_total;
                if (alpha1 <= n_incep)
                {
                    index1 = m_mech->SetRandomParticle(sys.Particles(), t, alpha1, iUniform, rng);
                    ip1 = -2;
                } else {
                    ip1 = sys.Particles().Select_usingGivenRand(iUniform, alpha1 - n_incep, rng);
                }
                maj = FreeMol;
                break;
            case FreeMol2:
                alpha1 = unifDistrib() * (dc2_incep + dc2_other);
                if (alpha1 <= dc2_incep)
                {
                    index1 = m_mech->SetRandomParticle(sys.Particles(), t, alpha1, iD2, rng);
                    ip1 = -2;
                } else {
                    ip1 = sys.Particles().Select_usingGivenRand(iD2, alpha1 - dc2_incep, rng);
                }
                maj = FreeMol;
                break;
            default:
                alpha1 = unifDistrib() * n_total;
                if (alpha1 <= n_incep)
                {
                    index1 = m_mech->SetRandomParticle(sys.Particles(), t, alpha1, iUniform, rng);
                    ip1 = -2;
                } else {
                    ip1 = sys.Particles().Select_usingGivenRand(iUniform, alpha1 - n_incep, rng);
                }
                maj = SlipFlow;
                break;
        }

        // Is this a particle-number particle?
        if (ip1 == -2)
        {
            if (index1 == 0)
            {
		// Property sums update triggered - round-off was too severe and
                // no suitable particle could be found. The easiest option
                // is not to perform this coagulation event.
		return -1;
            }
	    n_index1 = sys.Particles().NumberAtIndex(index1);
            // Flag sp1 as a particle-number particle
            ip1_flag = true; 
            // Clone a template particle of the correct size for sp1
            sp1 = sys.Particles().GetPNParticleAt(index1)->Clone();
            sp1->SetTime(t); 
        } else {
            // An ensemble particle
            if (ip1 >= 0) {
                sp1 = sys.Particles().At(ip1);
            } else {
            // Failed to choose a particle.
            return -1;
            }
        }

        // Choose and get unique second particle.  Note, we are allowed to do
        // this even if the first particle was invalidated.
        ip2 = ip1;
        index2 = index1;
        unsigned int guard = 0;
        bool mustSwitch = (ip1 == -2) && (n_incep == 1); // sp1 was the only list particle
        bool unsuitableChoice = true;

        switch (term) {
           case SlipFlow1:
               while (unsuitableChoice && (++guard < 1000))
                {
                    alpha2 = unifDistrib() * n_total;
                    if (alpha2 <= n_incep && !mustSwitch)
                    {
                        index2 = m_mech->SetRandomParticle(sys.Particles(), t, alpha2, iUniform, rng);
                        if (!(index2 == index1 && n_index1 == 1))
                        {
                            unsuitableChoice = false;
                            ip2 = -2;
                        }
                    } else {
                        ip2 = sys.Particles().Select_usingGivenRand(iUniform, alpha2 - n_incep, rng);
                        if (!(ip2 == ip1))
                            unsuitableChoice = false;
                    }
                }
                break;
            case SlipFlow2:
                while (unsuitableChoice && (++guard < 1000))
                {
                    alpha2 = unifDistrib() * (dc_1_incep + dc_1_other);
                    if (alpha2 <= dc_1_incep && !mustSwitch)
                    {
                        index2 = m_mech->SetRandomParticle(sys.Particles(), t, alpha2, iD_1, rng);
                        if (!(index2 == index1 && n_index1 == 1))
                        {
                            unsuitableChoice = false;
                            ip2 = -2;
                        }
                    } else {
                        ip2 = sys.Particles().Select_usingGivenRand(iD_1, alpha2 - dc_1_incep, rng);
                        if (!(ip2 == ip1))
                            unsuitableChoice = false;
                    }
                }
                break;
            case SlipFlow3:
                while (unsuitableChoice && (++guard < 1000))
                {
                    alpha2 = unifDistrib() * (dc_1_incep + dc_1_other);
                    if (alpha2 <= dc_1_incep && !mustSwitch)
                    {
                        index2 = m_mech->SetRandomParticle(sys.Particles(), t, alpha2, iD_1, rng);
                        if (!(index2 == index1 && n_index1 == 1))
                        {
                            unsuitableChoice = false;
                            ip2 = -2;
                        }
                    } else {
                        ip2 = sys.Particles().Select_usingGivenRand(iD_1, alpha2 - dc_1_incep, rng);
                        if (!(ip2 == ip1))
                            unsuitableChoice = false;
                    }
                }
                break;
            case SlipFlow4:
                while (unsuitableChoice && (++guard < 1000))
                {
                    alpha2 = unifDistrib() * (dc_2_incep + dc_2_other);
                    if (alpha2 <= dc_2_incep && !mustSwitch)
                    {
                        index2 = m_mech->SetRandomParticle(sys.Particles(), t, alpha2, iD_2, rng);
                        if (!(index2 == index1 && n_index1 == 1))
                        {
                            unsuitableChoice = false;
                            ip2 = -2;
                        }
                    } else {
                        ip2 = sys.Particles().Select_usingGivenRand(iD_2, alpha2 - dc_2_incep, rng);
                        if (!(ip2 == ip1))
                            unsuitableChoice = false;
                    }
                }
                break;
            case FreeMol1:
                while (unsuitableChoice && (++guard < 1000))
                {
                    alpha2 = unifDistrib() * (dc2_m_1_2_incep + dc2_m_1_2_other);
                    if (alpha2 <= dc2_m_1_2_incep && !mustSwitch)
                    {
                        index2 = m_mech->SetRandomParticle(sys.Particles(), t, alpha2, iD2_M_1_2, rng);
                        if (!(index2 == index1 && n_index1 == 1))
                        {
                            unsuitableChoice = false;
                            ip2 = -2;
                        }
                    } else {
                        ip2 = sys.Particles().Select_usingGivenRand(iD2_M_1_2, alpha2 - dc2_m_1_2_incep, rng);
                        if (!(ip2 == ip1))
                            unsuitableChoice = false;
                    }
                }
                break;
            case FreeMol2:
                while (unsuitableChoice && (++guard < 1000))
                {
                    alpha2 = unifDistrib() * (m_1_2_incep + m_1_2_other);
                    if (alpha2 <= m_1_2_incep && !mustSwitch)
                    {
                        index2 = m_mech->SetRandomParticle(sys.Particles(), t, alpha2, iM_1_2, rng);
                        if (!(index2 == index1 && n_index1 == 1))
                        {
                            unsuitableChoice = false;
                            ip2 = -2;
                        }
                    } else {
                        ip2 = sys.Particles().Select_usingGivenRand(iM_1_2, alpha2 - m_1_2_incep, rng);
                        if (!(ip2 == ip1))
                            unsuitableChoice = false;
                        }
                }
                break;
            default:
                while (unsuitableChoice && (++guard < 1000))
                {
                    alpha2 = unifDistrib() * n_total;
                    if (alpha2 <= n_incep && !mustSwitch)
                    {
                        index2 = m_mech->SetRandomParticle(sys.Particles(), t, alpha2, iUniform, rng);
                        if (!(index2 == index1 && n_index1 == 1))
                        {
                            unsuitableChoice = false;
                            ip2 = -2;
                        }
                    } else {
                        ip2 = sys.Particles().Select_usingGivenRand(iUniform, alpha2 - n_incep, rng);
                        if (!(ip2 == ip1))
                            unsuitableChoice = false;
                    }
                }
                break;
            }
        }

        // Is this a particle-number particle?
        if (ip2 == -2)
        {
            // Note this could be factored into the loops above
            // but that would require constantly checking for something
            // that should very infrequently be a problem.
            // Easiest option is to exit without performing the coagulation event. 
            if (index2 == 0)
            {
                if (ip1_flag)
                {
                    delete sp1;
                    sp1 = NULL;
                }
                // Property sums updated
                return -1;
            }
            ip2_flag = true; // Flag sp2 as a particle-number particle
            sp2 = sys.Particles().GetPNParticleAt(index2)->Clone();
            sp2->SetTime(t); 
        }
        else
        {
            // An ensemble particle
            if ((ip2 >= 0) && (ip2 != ip1)) {
                sp2 = sys.Particles().At(ip2);
            } else {
                // Failed to select a unique particle.
                return -1;
            }
        }

        //Calculate the majorant rate before updating the particles
        double majk = MajorantKernel(*sp1, *sp2, sys, maj);

        //Update the particles
        if (t > sp1->LastUpdateTime())
            m_mech->UpdateParticle(*sp1, sys, t, ip1, rng, dummy);
        // Check that particle is still valid.  If not,
        // remove it and cease coagulating.
        if (!sp1->IsValid()) {
            if (!ip1_flag) {
                // Must remove first particle now.
                sys.Particles().Remove(ip1);
            } else {
                delete sp1;
                sp1 = NULL;
            }
            // Invalidating the index tells this routine not to perform coagulation.
            ip1 = -1;
            return 0;
        }

        if (t > sp2->LastUpdateTime())
            m_mech->UpdateParticle(*sp2, sys, t, ip2, rng, dummy);
        // Check validity of particles after update.
        if (!sp2->IsValid()) {
            // Tell the ensemble to update particle one before we confuse things
            // by removing particle 2
            if (!ip1_flag)
                sys.Particles().Update(ip1);
            // Must remove second particle now.
            if (!ip2_flag)
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
        if ((ip1 != -1) && (ip2 != -1)) 
        {
            // Must check for ficticious event now by comparing the original
            // majorant rate and the current (after updates) true rate.
            double truek = CoagKernel(*sp1, *sp2, sys);
            double ceff = 0;
            if (majk<truek)
                std::cout << "maj< true" << std::endl;

            //added by ms785 to include the collision efficiency in the calculation of the rate
            if (sys.ParticleModel()->AggModel() == AggModels::PAH_KMC_ID)
            {
                ceff = sys.ParticleModel()->CollisionEff(sp1, sp2);
                truek *= ceff;
            }

            //! For the purpose of checking consistency with the spherical soot
            //! model solved using the method of moments with interpolative closure
            //! which assumes that only pyrene (A4) is able to incept and condense.
            else if (sys.ParticleModel()->AggModel() == AggModels::BinTree_ID || sys.ParticleModel()->AggModel() == AggModels::Spherical_ID) {
                if (sp1->Primary()->InceptedPAH() && sp2->Primary()->InceptedPAH()) {
                    ceff = 1;
                }
                else if (sp1->Primary()->InceptedPAH() && sp2->NumCarbon() > 16 || sp1->NumCarbon() > 16 && sp2->Primary()->InceptedPAH()) {
                    ceff = 1;
                } else {
                    ceff = 1;
                }
                truek *= ceff;
            }

            if (!Fictitious(majk, truek, rng)) {
                if (ip1_flag)
                {
                    // We are removing the particle from the PN model and
                    // adding it to the ensemble
                    sys.Particles().UpdateTotalsWithIndex(index1, -1.0);
                    sys.Particles().UpdateNumberAtIndex(index1, -1);
                    sys.Particles().UpdateTotalParticleNumber(-1);
                    unsigned int index12 = index1 + index2;
                    // Allow for coagulation in place if the combined particle is small enough
                    if ((m_mech->CoagulateInList()) && ip2_flag && (index12 < sys.Particles().GetHybridThreshold()))
                    {
                        coag_in_place = true;
                        sys.Particles().UpdateTotalsWithIndex(index12, 1.0);
                        sys.Particles().UpdateNumberAtIndex(index12, 1);
                        sys.Particles().UpdateTotalParticleNumber(1);
                        if (sp1 != NULL) 
                        {
                            delete sp1;
                            sp1 = NULL;
                        }
                    } else 
                        ip1 = sys.Particles().Add(*sp1, rng, ip2, true);    
                }
                if (ip2_flag)
                {
                    // We are removing the particle from the PN model and
                    // adding it to the ensemble
                    sys.Particles().UpdateTotalsWithIndex(index2, -1.0);
                    sys.Particles().UpdateNumberAtIndex(index2, -1);
                    sys.Particles().UpdateTotalParticleNumber(-1);
                }
                if (!coag_in_place)
                    JoinParticles(t, ip1, sp1, ip2, sp2, sys, rng);
			
                if (ip2_flag && sp2 != NULL)                                                     // Particle sp2 is not in the ensemble, must manually delete it
                {
                    delete sp2;
                    sp2 = NULL;
                }
            } else {
                if (!ip1_flag)
                    sys.Particles().Update(ip1);
                else if (sp1 != NULL)                                                            // Particle sp1 is not in the ensemble, must manually delete it
                {
                    delete sp1;
                    sp1 = NULL;
                }
                if (!ip2_flag)
                    sys.Particles().Update(ip2);
                else if (sp2 != NULL)                                                            // Particle sp2 is not in the ensemble, must manually delete it
                {
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
                if (!ip1_flag)
                    sys.Particles().Update(ip1);
            if (ip2 != -1)
            {
                if (!ip2_flag)
                    sys.Particles().Update(ip2);
            }
            if (ip1_flag && sp1 != NULL)                                                         // Particle sp2 is not in the ensemble, must manually delete it
            {
                delete sp1;
                sp1 = NULL;
            }
            if (ip2_flag && sp2 != NULL)                                                         // Particle sp2 is not in the ensemble, must manually delete it
            {
                delete sp2;
                sp2 = NULL;
            }
    	}
    return 0;
}


// COAGULATION KERNELS.

/**
 * Calculate the coagulation kernel between two particles in the given environment.
 *
 *@param[in]    sp1         First particle
 *@param[in]    sp2         Second particle
 *@param[in]    sys         Details of the environment, including temperature and pressure
 *
 *@return       Value of kernel
 */
double Sweep::Processes::HybridTransitionCoagulation::CoagKernel(const Particle &sp1,
                                                                const Particle &sp2,
                                                                const Cell &sys) const
{
    const double T = sys.GasPhase().Temperature();
    const double P = sys.GasPhase().Pressure();
    const double fm = FreeMolKernel(sp1, sp2, T, P, false);
    const double sf = SlipFlowKernel(sp1, sp2, T, P, sys.GasPhase().Viscosity(), false);
    return (fm*sf)/(fm+sf);
}

/**
 * Calculate the majorant kernel between two particles in the given environment.
 *
 *@param[in]    sp1         First particle
 *@param[in]    sp2         Second particle
 *@param[in]    sys         Details of the environment, including temperature and pressure
 *@param[in]    maj         Flag to indicate which majorant kernel is required
 *
 *@return       Value of kernel
 */
double Sweep::Processes::HybridTransitionCoagulation::MajorantKernel(const Particle &sp1,
                                                                    const Particle &sp2,
                                                                    const Cell &sys,
                                                                    const MajorantType maj) const
{
    // This routine calculates the coagulation kernel for two particles.  The kernel
    // type is chosen by the majorant type requested.
    switch (maj) {
        case Default:
            // This should never happen for the transition coagulation kernel
            assert(maj != Default);
            break;
        case FreeMol:
            // Free molecular majorant.
            return FreeMolKernel(sp1, sp2, sys.GasPhase().Temperature(),
                    sys.GasPhase().Pressure(), true);
        case SlipFlow:
            // Slip-flow majorant.
            return SlipFlowKernel(sp1, sp2, sys.GasPhase().Temperature(),
                    sys.GasPhase().Pressure(), sys.GasPhase().Viscosity(), true);
    }

    // Invalid majorant, return zero.
    return 0.0;
}
// Returns the free-molecular coagulation kernel value for the
// two given particles.  Can return either the majorant or
// true kernel.
double Sweep::Processes::HybridTransitionCoagulation::FreeMolKernel(const Particle &sp1, const Particle &sp2,
                                double T, double P, bool maj) const
{
    // This routine calculate the free molecular coagulation kernel for two particles.
    // There are two forms of kernel; a majorant form and a non-majorant form.

    // Collect the particle properties
    const double d1 = sp1.CollDiameter();
    const double d2 = sp2.CollDiameter();
    const double invm1 = 1.0 / sp1.Mass();
    const double invm2 = 1.0 / sp2.Mass();

    if (maj) {
        // The majorant form is always >= the non-majorant form.
        return CFMMAJ * m_efm * CFM * sqrt(T) * A() *
               (std::sqrt(invm1) + std::sqrt(invm2)) *
               (d1 * d1 + d2 * d2);
    } else {
        const double dterm = d1 + d2;
        return m_efm * CFM * A() *
               sqrt(T * (invm1 + invm2)) *
               dterm * dterm;
    }
}

// Returns the slip-flow coagulation kernel value for the
// two given particles.  Can return either the majorant or
// true kernel.
double Sweep::Processes::HybridTransitionCoagulation::SlipFlowKernel(const Particle &sp1, const Particle &sp2,
                                 double T, double P, double mu, bool maj) const
{
    // Collect the particle properties
    const double d1 = sp1.CollDiameter();
    const double d2 = sp2.CollDiameter();

    // For the slip-flow kernel the majorant and non-majorant forms are identical.
    return ((1.257 * 2.0 * MeanFreePathAir(T,P) *
             (1.0 / d1 / d1 + 1.0 / d2 / d2)) +
            (1.0 / d1 + 1.0 / d2)) *
           CSF * T * (d1 + d2)
           * A() / mu;
}

