/*
  Author(s):      Robert I A Patterson
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2009 Robert I A Patterson.

  File purpose:
    Implementation of transition regime coagulation kernel

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
#include "swp_transcoag_weighted_PAHs.h"

#include "swp_params.h"
#include "swp_cell.h"
#include "swp_mechanism.h"
#include "swp_PAH_primary.h"

using namespace Sweep::Processes;


// Default constructor.
Sweep::Processes::TransitionCoagulationWeightedPAHs::TransitionCoagulationWeightedPAHs(const Sweep::Mechanism &mech)
: Coagulation(mech), m_efm(mech.GetEnhancementFM())
{
    m_name = "TransitionRegimeCoagulationwithWeightedPAHs";
}

Sweep::Processes::TransitionCoagulationWeightedPAHs* const Sweep::Processes::TransitionCoagulationWeightedPAHs::Clone() const
{
	return new TransitionCoagulationWeightedPAHs(*this);
}

// Stream-reading constructor.
Sweep::Processes::TransitionCoagulationWeightedPAHs::TransitionCoagulationWeightedPAHs(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech), m_efm(mech.GetEnhancementFM())
{
    m_name = "TransitionRegimeCoagulation";
    Deserialize(in, mech);
}

// TOTAL RATE CALCULATION.

// Returns the rate of the process for the given system.
double Sweep::Processes::TransitionCoagulationWeightedPAHs::Rate(double t, const Cell &sys,
                                                          const Geometry::LocalGeometry1d &local_geom) const
{
    // Get the number of particles in the system.
	unsigned int n = sys.Particles().GetSum(Sweep::iW);

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        double T = sys.GasPhase().Temperature();
        double P = sys.GasPhase().Pressure();

        // Calculate the rate.
        return Rate(sys.Particles().GetSums(), (double)n, sqrt(T),
                    T/sys.GasPhase().Viscosity(), MeanFreePathAir(T,P),
                    sys.SampleVolume());
    } else {
        return 0.0;
    }
}

// More efficient rate routine for coagulation only.
// All parameters required to calculate rate passed
// as arguments.
double Sweep::Processes::TransitionCoagulationWeightedPAHs::Rate(const Ensemble::particle_cache_type &data, double n, double sqrtT,
                       double T_mu, double MFP, double vol) const
{
    // Some prerequisites.
    double n_1 = n - 1.0;
    double a = CSF * T_mu * A();
    double b = a * MFP * 1.257 * A();
    double c = CFMMAJ * m_efm * CFM * sqrtT * A();

    // Summed particle properties required for coagulation rate.
    const double d       = data.Property(Sweep::iDW);
    const double d2      = data.Property(Sweep::iD2W);
    const double d_1     = data.Property(Sweep::iD_1W);
    const double d_2     = data.Property(Sweep::iD_2W);
    const double m_1_2   = data.Property(Sweep::iM_1_2W);
    const double d2m_1_2 = data.Property(Sweep::iD2_M_1_2W);

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
unsigned int Sweep::Processes::TransitionCoagulationWeightedPAHs::TermCount(void) const { return TYPE_COUNT; }

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
double Sweep::Processes::TransitionCoagulationWeightedPAHs::RateTerms(double t, const Cell &sys,
                            const Geometry::LocalGeometry1d &local_geom,
                            fvector::iterator &iterm) const
{
    // Get the number of particles in the system.
	unsigned int n = sys.Particles().GetSum(Sweep::iW);

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        double T = sys.GasPhase().Temperature();
        double P = sys.GasPhase().Pressure();

        // Calculate the rate terms.
        return RateTerms(sys.Particles().GetSums(), (double)n, sqrt(T), T/sys.GasPhase().Viscosity(),
                         MeanFreePathAir(T,P), sys.SampleVolume(), iterm);
    } else {
        // No coagulation as there are too few particles.
        for(unsigned int i=0; i<TYPE_COUNT; i++) *(iterm++) = 0.0;
        return 0.0;
    }
}

// More efficient rate routine for coagulation only.
// All parameters required to calculate rate terms
// passed as arguments.
double Sweep::Processes::TransitionCoagulationWeightedPAHs::RateTerms(const Ensemble::particle_cache_type &data, double n, double sqrtT,
                            double T_mu, double MFP, double vol,
                            fvector::iterator &iterm) const
{
    // Some prerequisites.
    double n_1 = n - 1.0;
    double a   = CSF * T_mu * A();
    double b   = a * MFP * 1.257 * 2.0;
    double c   = CFMMAJ * m_efm * CFM * sqrtT * A();

    // Summed particle properties required for coagulation rate.
	const double d = data.Property(Sweep::iDW);
	const double d2 = data.Property(Sweep::iD2W);
	const double d_1 = data.Property(Sweep::iD_1W);
	const double d_2 = data.Property(Sweep::iD_2W);
	const double m_1_2 = data.Property(Sweep::iM_1_2W);
	const double d2m_1_2 = data.Property(Sweep::iD2_M_1_2W);

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
int TransitionCoagulationWeightedPAHs::Perform(double t, Sweep::Cell &sys,
                                   const Geometry::LocalGeometry1d& local_geom,
                                   unsigned int iterm,
                                   Sweep::rng_type &rng) const
{
	PartPtrVector dummy;

	bool new1 = false, new2 = false;
    // Select properties by which to choose particles (-1 means
    // choose uniformly).  Note we need to choose 2 particles.  There
    // are six possible rate terms to choose from; 4 slip-flow and 2
    // free molecular.
	if (sys.Particles().GetSum(Sweep::iW) < 2) {
        return 1;
    }

    int ip1=-1, ip2=-1;
    MajorantType maj;
    TermType term = (TermType)iterm;

    // Select the first particle and note the majorant type.
    switch (term) {
        case SlipFlow1:
			ip1 = sys.Particles().Select(Sweep::iUniform1, rng);
            maj = SlipFlow;
            break;
        case SlipFlow2:
            ip1 = sys.Particles().Select(Sweep::iDW, rng);
            maj = SlipFlow;
            break;
        case SlipFlow3:
			ip1 = sys.Particles().Select(Sweep::iUniform1, rng);
            maj = SlipFlow;
            break;
        case SlipFlow4:
            ip1 = sys.Particles().Select(Sweep::iDW, rng);
            maj = SlipFlow;
            break;
        case FreeMol1:
			ip1 = sys.Particles().Select(Sweep::iUniform1, rng);
            maj = FreeMol;
            break;
        case FreeMol2:
            ip1 = sys.Particles().Select(Sweep::iD2W, rng);
            maj = FreeMol;
            break;
        default :
			ip1 = sys.Particles().Select(Sweep::iUniform1, rng);
            maj = SlipFlow;
            break;
    }

    // Choose and get first particle.
    Particle *sp1=NULL;
    if (ip1 >= 0) {
        sp1 = sys.Particles().At(ip1);
    } else {
        // Failed to choose a particle.
        return -1;
    }

    // Choose and get unique second particle.  Note, we are allowed to do
    // this even if the first particle was invalidated.
    ip2 = ip1;
	unsigned int max;
	if (sys.Particles().At(ip1)->getStatisticalWeight() > 1.0){
		max = 1;
	}
	else{
		max = 1000;
	}
    unsigned int guard = -1;
    switch (term) {
        case SlipFlow1:
            while ((ip2 == ip1) && (++guard<max))
				ip2 = sys.Particles().Select(Sweep::iUniform1, rng);
            break;
        case SlipFlow2:
			while ((ip2 == ip1) && (++guard<max))
                ip2 = sys.Particles().Select(Sweep::iD_1W, rng);
            break;
        case SlipFlow3:
			while ((ip2 == ip1) && (++guard<max))
                ip2 = sys.Particles().Select(Sweep::iD_1W, rng);
            break;
        case SlipFlow4:
			while ((ip2 == ip1) && (++guard<max))
                ip2 = sys.Particles().Select(Sweep::iD_2W, rng);
            break;
        case FreeMol1:
			while ((ip2 == ip1) && (++guard<max))
                ip2 = sys.Particles().Select(Sweep::iD2_M_1_2W, rng);
            break;
        case FreeMol2:
			while ((ip2 == ip1) && (++guard<max))
                ip2 = sys.Particles().Select(Sweep::iM_1_2W, rng);
            break;
        default :
			while ((ip2 == ip1) && (++guard<max))
				ip2 = sys.Particles().Select(Sweep::iUniform1, rng);
            break;
    }

	Particle *sp2 = NULL;
	if ((ip2 >= 0) && ((ip2 != ip1) || sys.Particles().At(ip1)->getStatisticalWeight() > 1)) {
		sp2 = sys.Particles().At(ip2);
	}
	else {
		// Failed to select a unique particle.
		return -1;
	}

	Particle part1 = Particle(t, (*sys.ParticleModel()));
	Particle part2 = Particle(t, (*sys.ParticleModel()));
	if (sp1->getStatisticalWeight() > 1){
		part1 = *sp1;
		part1.setStatisticalWeight(1);
		sp1 = &part1;
		new1 = true;
	}
	if (sp2->getStatisticalWeight() > 1){
		part2 = *sp2;
		part2.setStatisticalWeight(1);
		sp2 = &part2;
		new2 = true;
	}

    //Calculate the majorant rate before updating the particles
    double majk = MajorantKernel(*sp1, *sp2, sys, maj);

    //Update the particles
    m_mech->UpdateParticle(*sp1, sys, t, ip1, rng, dummy);
    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp1->IsValid()) {
        // Must remove first particle now.
		if (!new1){
			sys.Particles().Remove(ip1);
		}
		else{
			int oldweight = sys.Particles().At(ip1)->getStatisticalWeight();
			int newweight = oldweight - 1;
			if (newweight > 0){
				sys.Particles().At(ip1)->setStatisticalWeight(newweight);
				sys.Particles().Update(ip1);
			}
			else{
				sys.Particles().Remove(ip1);
			}
		}

        // Invalidating the index tells this routine not to perform coagulation.
        ip1 = -1;
        return 0;
    }

    m_mech->UpdateParticle(*sp2, sys, t, ip2, rng, dummy);
    // Check validity of particles after update.
    if (!sp2->IsValid()) {
		// Tell the ensemble to update particle one before we confuse things
		// by removing particle 2 (if sp1 was not just a copy of the particle at ip1)
		if (!new1){
			sys.Particles().Update(ip1);
		}

		// Must remove second particle now.
		if (!new2){
			sys.Particles().Remove(ip2);
		}
		else{
			int oldweight = sys.Particles().At(ip2)->getStatisticalWeight();
			int newweight = oldweight - 1;
			if (newweight > 0){
				sys.Particles().At(ip2)->setStatisticalWeight(newweight);
				sys.Particles().Update(ip2);
			}
			else{
				sys.Particles().Remove(ip2);
			}
		}

		if (new1){
			Particle* adder = new Particle(*sp1);
			int newip1;
			if (sys.ParticleCount() < sys.Particles().Capacity()){
				newip1 = sys.Particles().Add(*adder, rng);
				sys.Particles().Update(newip1);
			}
			else
			{
				std::cout << "No room in ensemble after coagulation" << std::endl;
			}
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
        double ceff=0;
        if (majk<truek)
            std::cout << "maj< true"<< std::endl;

        //added by ms785 to include the collision efficiency in the calculation of the rate
        if (sys.ParticleModel()->AggModel() == AggModels::PAH_KMC_ID)
        {
            ceff=sys.ParticleModel()->CollisionEff(sp1,sp2);
            truek*=ceff;
        }

        //! For the purpose of checking consistency with the spherical soot
        //! model solved using the method of moments with interpolative closure
        //! which assumes that only pyrene (A4) is able to incept and condense.
        else if (sys.ParticleModel()->AggModel() == AggModels::BinTree_ID || sys.ParticleModel()->AggModel() == AggModels::Spherical_ID) {
            if (sp1->Primary()->InceptedPAH() && sp2->Primary()->InceptedPAH()) {
                ceff = 1;
            } else if (sp1->Primary()->InceptedPAH() && sp2->NumCarbon() > 16 || sp1->NumCarbon() > 16 && sp2->Primary()->InceptedPAH()) {
                ceff = 1;
            } else {
                ceff = 1;
            }
            truek*=ceff;
        }
        
        if (!Fictitious(majk, truek, rng)) {
			if (new1 && new2){ //This coagulation was between two weighted monomer PAHs
				if (ip1 == ip2){ //This was a coagulation between the same weight monomer
					//First, reduce the weight of the original particle or remove it
					int oldweight = (*sys.Particles().At(ip1)).getStatisticalWeight();
					int newweight = oldweight - 2;
					if (newweight > 0){
						(*sys.Particles().At(ip1)).setStatisticalWeight(newweight);
						sys.Particles().Update(ip1);
					}
					else{
						sys.Particles().Remove(ip1);
					}

					//Now coagulate the particles and add sp1 to the ensemble
					Particle* adder = new Particle(*sp1);
					int ipnew1;
					if (sys.ParticleCount() < sys.Particles().Capacity()){
						ipnew1 = sys.Particles().Add(*adder, rng);
						sys.Particles().Update(ipnew1);
					}
					else
					{
						std::cout << "No room in ensemble after coagulation" << std::endl;
					}

					adder->Coagulate(*sp2, rng);
					adder->SetTime(t);
					adder->incrementCoagCount();
					sys.Particles().Update(ipnew1);

				}
				else { //Coagulation between two different weighted monomers
					//First, reduce the weight of the original particles or remove them
					int oldweight = (*sys.Particles().At(ip1)).getStatisticalWeight();
					int newweight = oldweight - 1;
					if (newweight > 0){
						(*sys.Particles().At(ip1)).setStatisticalWeight(newweight);
						sys.Particles().Update(ip1);
					}

					oldweight = (*sys.Particles().At(ip2)).getStatisticalWeight();
					int newweight2 = oldweight - 1;
					if (newweight2 > 0){
						(*sys.Particles().At(ip2)).setStatisticalWeight(newweight2);
						sys.Particles().Update(ip2);
					}
					else{
						sys.Particles().Remove(ip2);
					}

					if(newweight <= 0){
						if (ip2 > ip1){
							sys.Particles().Remove(ip1);
						}
						else{
							if (newweight2 > 0){
								sys.Particles().Remove(ip1);
							}
							else{
								sys.Particles().Remove(ip1 - 1);
							}
						}
					}

					//Now coagulate the particles and add sp1 to the ensemble
					Particle* adder = new Particle(*sp1);
					int ipnew1;
					if (sys.ParticleCount() < sys.Particles().Capacity()){
						ipnew1 = sys.Particles().Add(*adder, rng);
						sys.Particles().Update(ipnew1);
					}
					else
					{
						std::cout << "No room in ensemble after coagulation" << std::endl;
					}

					adder->Coagulate(*sp2, rng);
					adder->SetTime(t);
					adder->incrementCoagCount();
					sys.Particles().Update(ipnew1);
				}
			}
			else if (new1 || new2){ //if one of the particles was a weighted monomer PAH
				int index;
				bool first;
				if (new1){ //First one was the weighted monomer PAH
					index = ip1;
					first = true;
				}
				else{ //It was the second one
					index = ip2;
					first = false;
				}
				//reduce the weight of the original particle or remove it
				int oldweight = (*sys.Particles().At(index)).getStatisticalWeight();
				int newweight = oldweight - 1;
				if (newweight > 0){
					(*sys.Particles().At(index)).setStatisticalWeight(newweight);
					sys.Particles().Update(index);
				}
				if (first){//if the first was the weighted monomer, add particle 1 to particle 2
					sp2->Coagulate(*sp1, rng);
					sp2->SetTime(t);
					sp2->incrementCoagCount();
					sys.Particles().Update(ip2);
				}
				else{ //else, add particle 2 to particle 1
					sp1->Coagulate(*sp2, rng);
					sp1->SetTime(t);
					sp1->incrementCoagCount();
					sys.Particles().Update(ip1);
				}

				if (newweight <= 0){
					sys.Particles().Remove(index);
				}
			}
			else{ //neither particle was a weight monomer. perform coagulation as normal with a non-weighted algorithm
				JoinParticles(t, ip1, sp1, ip2, sp2, sys, rng);
			}
        } else {
			sys.Particles().Update(ip1);
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

        if(ip2 >= 0)
            sys.Particles().Update(ip2);
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
double Sweep::Processes::TransitionCoagulationWeightedPAHs::CoagKernel(const Particle &sp1,
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
double Sweep::Processes::TransitionCoagulationWeightedPAHs::MajorantKernel(const Particle &sp1,
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
double Sweep::Processes::TransitionCoagulationWeightedPAHs::FreeMolKernel(const Particle &sp1, const Particle &sp2,
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
double Sweep::Processes::TransitionCoagulationWeightedPAHs::SlipFlowKernel(const Particle &sp1, const Particle &sp2,
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

