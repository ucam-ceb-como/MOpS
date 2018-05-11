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
#include "swp_transcoag.h"

#include "swp_params.h"
#include "swp_cell.h"
#include "swp_mechanism.h"
#include "swp_PAH_primary.h"
#include <boost/random/uniform_01.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Sweep::Processes;


// Default constructor.
Sweep::Processes::TransitionCoagulation::TransitionCoagulation(const Sweep::Mechanism &mech)
: Coagulation(mech), m_efm(mech.GetEnhancementFM())
{
	m_name = "TransitionRegimeCoagulation";
}

Sweep::Processes::TransitionCoagulation* const Sweep::Processes::TransitionCoagulation::Clone() const
{
	return new TransitionCoagulation(*this);
}

// Stream-reading constructor.
Sweep::Processes::TransitionCoagulation::TransitionCoagulation(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech), m_efm(mech.GetEnhancementFM())
{
	m_name = "TransitionRegimeCoagulation";
	Deserialize(in, mech);
}

// TOTAL RATE CALCULATION.

// Returns the rate of the process for the given system.
double Sweep::Processes::TransitionCoagulation::Rate(double t, const Cell &sys,
	const Geometry::LocalGeometry1d &local_geom) const
{
	// Get the number of particles in the system.
	unsigned int n = sys.ParticleCount();

	// aab64 for hybrid particle model
	if (m_mech->IsHybrid())
		n += sys.GetIncepted();

	// Check that there are at least 2 particles before calculating rate.
	if (n > 1) {
		// Get system properties required to calculate coagulation rate.
		double T = sys.GasPhase().Temperature();
		double P = sys.GasPhase().Pressure();

		// Calculate the rate.
		if (m_mech->IsHybrid() && sys.GetIncepted() > 0)
		{
			Particle * spInc_tmp1 = sys.Particles().GetInceptedSP_tmp().Clone();
			Particle * spInc_inc = sys.Particles().GetInceptedSP().Clone();

			spInc_inc->setStatisticalWeight(sys.GetInceptions_tmp());
			spInc_inc->SetTime(t);
			double tmp_age = sys.GetIncepted() - sys.GetInceptions_tmp();
			tmp_age = (tmp_age > 0) ? tmp_age : 0.0;
			spInc_tmp1->setStatisticalWeight(tmp_age);

			double rate = Rate(sys.Particles().GetSums(), (double)n, sqrt(T),
				T / sys.GasPhase().Viscosity(), MeanFreePathAir(T, P),
				sys.SampleVolume(), spInc_tmp1, spInc_inc); //sys.Particles().GetInceptedSP_tmp()
			//std::cout << "Coag: " << t << " , " << rate << std::endl;

			delete spInc_inc;
			spInc_inc = NULL;
			delete spInc_tmp1;
			spInc_tmp1 = NULL;

			return rate;
		}
		else
		{
			double rate = Rate(sys.Particles().GetSums(), (double)n, sqrt(T),
				T / sys.GasPhase().Viscosity(), MeanFreePathAir(T, P),
				sys.SampleVolume(), nullptr, nullptr);
			//std::cout << "Coag: " << t << " , " << rate << std::endl;
			return rate;
		}
	}
	else {
		return 0.0;
	}
}

// More efficient rate routine for coagulation only.
// All parameters required to calculate rate passed
// as arguments.
double Sweep::Processes::TransitionCoagulation::Rate(const Ensemble::particle_cache_type &data, double n, double sqrtT,
	double T_mu, double MFP, double vol, Particle * spInc_tmp1, Particle * spInc_inc) const
{
	// Some prerequisites.
	double n_1 = n - 1.0;
	double a = CSF * T_mu * A();
	double b = a * MFP * 1.257 * A();
	double c = CFMMAJ * m_efm * CFM * sqrtT * A();

	// Summed particle properties required for coagulation rate.
	double d = data.Property(iDcol); // aab64 removed const to adapt as below
	double d2 = data.Property(iD2);
	double d_1 = data.Property(iD_1);
	double d_2 = data.Property(iD_2);
	double m_1_2 = data.Property(iM_1_2);
	double d2m_1_2 = data.Property(iD2_M_1_2);

	if (m_mech->IsHybrid())
	{
		/*double d_inc = 2.6766e-08;
		double d2_inc = 2.6827e-15;
		double d_1_inc = 5.0186e+08;
		double d_2_inc = 6.8816e+17;
		double m_1_2_inc = 3.7800e+11;
		double d2m_1_2_inc = 1.0808e-23;
		double wt = (spInc_tmp1 != nullptr) ? spInc_tmp1->getStatisticalWeight() : 0.0;
		if (spInc_inc != nullptr)
			wt += spInc_inc->getStatisticalWeight();
		d += (wt * d_inc);
		d2 += (wt * d2_inc);
		d_1 += (wt * d_1_inc);
		d_2 += (wt * d_2_inc);
		m_1_2 += (wt * m_1_2_inc);
		d2m_1_2 += (wt * d2m_1_2_inc);*/
		if (spInc_tmp1 != nullptr)
		{
			double wt = spInc_tmp1->getStatisticalWeight();
			double dc = spInc_tmp1->CollDiameter();
			double inv_sqrt_m = 1.0 / std::sqrt(spInc_tmp1->Mass());
			double dc_sqrd = dc * dc;
			d += (wt * dc);
			d2 += (wt * dc_sqrd);
			d_1 += (wt / dc);
			d_2 += (wt / dc_sqrd);
			m_1_2 += (wt * inv_sqrt_m);
			d2m_1_2 += (wt * dc_sqrd * inv_sqrt_m);
		}
		if (m_mech->IsHybrid() && spInc_inc != nullptr)
		{
			double wt = spInc_inc->getStatisticalWeight();
			double dc = spInc_inc->CollDiameter();
			double inv_sqrt_m = 1.0 / std::sqrt(spInc_inc->Mass());
			double dc_sqrd = dc * dc;
			d += (wt * dc);
			d2 += (wt * dc_sqrd);
			d_1 += (wt / dc);
			d_2 += (wt / dc_sqrd);
			m_1_2 += (wt * inv_sqrt_m);
			d2m_1_2 += (wt * dc_sqrd * inv_sqrt_m);
		}
	}

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
		}
		else {
			// Use slip-flow majorant.
			return sf;
		}
	}
	return 0.0;
}


// RATE TERM CALCULATION.

// Returns the number of rate terms for this process.
unsigned int Sweep::Processes::TransitionCoagulation::TermCount(void) const { return TYPE_COUNT; }

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
double Sweep::Processes::TransitionCoagulation::RateTerms(double t, const Cell &sys,
	const Geometry::LocalGeometry1d &local_geom,
	fvector::iterator &iterm) const
{
	// Get the number of particles in the system.
	unsigned int n = sys.ParticleCount();

	// aab64 for hybrid particle model
	if (m_mech->IsHybrid())
		n += sys.GetIncepted();

	// Check that there are at least 2 particles before calculating rate.
	if (n > 1) {
		// Get system properties required to calculate coagulation rate.
		double T = sys.GasPhase().Temperature();
		double P = sys.GasPhase().Pressure();

		// Calculate the rate terms.
		if (m_mech->IsHybrid() && sys.GetIncepted() > 0)
		{
			Particle * spInc_tmp1 = sys.Particles().GetInceptedSP_tmp().Clone();
			Particle * spInc_inc = sys.Particles().GetInceptedSP().Clone();

			spInc_inc->setStatisticalWeight(sys.GetInceptions_tmp());
			spInc_inc->SetTime(t);
			double tmp_age = sys.GetIncepted() - sys.GetInceptions_tmp();
			tmp_age = (tmp_age > 0) ? tmp_age : 0.0;
			spInc_tmp1->setStatisticalWeight(tmp_age);

			double rate = RateTerms(sys.Particles().GetSums(), (double)n, sqrt(T), T / sys.GasPhase().Viscosity(),
				MeanFreePathAir(T, P), sys.SampleVolume(), iterm, spInc_tmp1, spInc_inc); //&sys.Particles().GetInceptedSP_tmp()

			delete spInc_inc;
			spInc_inc = NULL;
			delete spInc_tmp1;
			spInc_tmp1 = NULL;

			return rate;
		}
		else
		{
			double rate = RateTerms(sys.Particles().GetSums(), (double)n, sqrt(T), T / sys.GasPhase().Viscosity(),
				MeanFreePathAir(T, P), sys.SampleVolume(), iterm, nullptr, nullptr);
			//std::cout << "Coag: " << t << " , " << rate << std::endl;
			return rate;
		}
	}
	else {
		// No coagulation as there are too few particles.
		for (unsigned int i = 0; i<TYPE_COUNT; i++) *(iterm++) = 0.0;
		return 0.0;
	}
}

// More efficient rate routine for coagulation only.
// All parameters required to calculate rate terms
// passed as arguments.
double Sweep::Processes::TransitionCoagulation::RateTerms(const Ensemble::particle_cache_type &data, double n, double sqrtT,
	double T_mu, double MFP, double vol,
	fvector::iterator &iterm, Particle * spInc_tmp1, Particle * spInc_inc) const
{
	// Some prerequisites.
	double n_1 = n - 1.0;
	double a = CSF * T_mu * A();
	double b = a * MFP * 1.257 * 2.0;
	double c = CFMMAJ * m_efm * CFM * sqrtT * A();

	// Summed particle properties required for coagulation rate.
	double d = data.Property(iDcol); // aab64 removed const to adapt as below
	double d2 = data.Property(iD2);
	double d_1 = data.Property(iD_1);
	double d_2 = data.Property(iD_2);
	double m_1_2 = data.Property(iM_1_2);
	double d2m_1_2 = data.Property(iD2_M_1_2);
	
	if (m_mech->IsHybrid())
	{
		/*double d_inc = 2.7081e-08;
		double d2_inc = 2.7183e-15;
		double d_1_inc = 4.9635e+08;
		double d_2_inc = 6.8089e+17;
		double m_1_2_inc = 3.7392e+11;
		double d2m_1_2_inc = 2.6527e-06;*/
		/*double d_inc = 7.4509e-08;
		double d2_inc = 4.2153e-13;// 4.0256e - 14;
		double d_1_inc = 1.0191e+09;
		double d_2_inc = 7.8855e+19;// 1.4092e + 19;
		double m_1_2_inc = 3.4935e+12;// 1.4064e + 12;
		double d2m_1_2_inc = 3.3639e-06;// 2.7389e - 06;
		double wt = (spInc_tmp1 != nullptr) ? spInc_tmp1->getStatisticalWeight() : 0.0;
		if (spInc_inc != nullptr)
			wt += spInc_inc->getStatisticalWeight();
		d += (wt * d_inc);
		d2 += (wt * d2_inc);
		d_1 += (wt * d_1_inc);
		d_2 += (wt * d_2_inc);
		m_1_2 += (wt * m_1_2_inc);
		d2m_1_2 += (wt * d2m_1_2_inc);*/

		if (spInc_tmp1 != nullptr)
		{
			double wt = spInc_tmp1->getStatisticalWeight();

			double dc = spInc_tmp1->CollDiameter();
			double inv_sqrt_m = 1.0 / std::sqrt(spInc_tmp1->Mass());

			double dc_sqrd = dc * dc;
			d += (wt * dc);
			d2 += (wt * dc_sqrd);
			d_1 += (wt / dc);
			d_2 += (wt / dc_sqrd);
			m_1_2 += (wt * inv_sqrt_m);
			d2m_1_2 += (wt * dc_sqrd * inv_sqrt_m);
		}
		if (spInc_inc != nullptr)
		{
			double wt = spInc_inc->getStatisticalWeight();

			double dc = spInc_inc->CollDiameter();
			double inv_sqrt_m = 1.0 / std::sqrt(spInc_inc->Mass());

			double dc_sqrd = dc * dc;
			d += (wt * dc);
			d2 += (wt * dc_sqrd);
			d_1 += (wt / dc);
			d_2 += (wt / dc_sqrd);
			m_1_2 += (wt * inv_sqrt_m);
			d2m_1_2 += (wt * dc_sqrd * inv_sqrt_m);
		}
	}

	fvector::iterator isf = iterm;
	fvector::iterator ifm = iterm + 4;

	// Slip-flow.
	*(iterm) = n * n_1 * a / vol;
	*(++iterm) = ((d * d_1) - n) * a / vol;
	*(++iterm) = d_1 * n_1 * b / vol;
	*(++iterm) = ((d * d_2) - d_1) * b / vol;
	// Free-molecular.
	*(++iterm) = n_1 * d2m_1_2  * c / vol;
	*(++iterm) = (m_1_2 * d2 - d2m_1_2) * c / vol;

	// Return iterator to next term after the coagulation terms.
	++iterm;

	// Sum up total coagulation rates for different regimes.
	double sf = *(isf)+*(isf + 1) + *(isf + 2) + *(isf + 3);
	double fm = *(ifm)+*(ifm + 1);

	if ((sf>0.0) || (fm>0.0)) {
		// There is some coagulation.
		if (sf > fm) {
			// Use free-mol majorant.
			*(isf) = 0.0;
			*(isf + 1) = 0.0;
			*(isf + 2) = 0.0;
			*(isf + 3) = 0.0;
			return fm;
		}
		else {
			// Use slip-flow majorant.
			*(ifm) = 0.0;
			*(ifm + 1) = 0.0;
			return sf;
		}
	}
	else {
		// Something went wrong with the rate calculation.
		*(isf) = 0.0;
		*(isf + 1) = 0.0;
		*(isf + 2) = 0.0;
		*(isf + 3) = 0.0;
		*(ifm) = 0.0;
		*(ifm + 1) = 0.0;
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
int TransitionCoagulation::Perform(double t, Sweep::Cell &sys,
	const Geometry::LocalGeometry1d& local_geom,
	unsigned int iterm,
	Sweep::rng_type &rng) const
{
	// aab64 hybrid particle model flags
	bool hybrid_flag = m_mech->IsHybrid();
	bool ip1_flag = false;
	bool ip2_flag = false;

	int ip1 = -1, ip2 = -1;

	//double create_t = hybrid_flag ? sys.Particles().GetInceptedSP().CreateTime() : 0;
	//double exist_t = t - create_t;
	//double inceptingcoagulationchange = (exist_t > 0) ? (sys.GetInceptionCoagulationChange() / exist_t) : 0;
	//double sp_age = 0;

	//if (!hybrid_flag)
	//double inceptingcoagulationchange = 1; // distribution can't handle zero
	//boost::exponential_distribution<double> waitDistrib(inceptingcoagulationchange);
	//boost::variate_generator<Sweep::rng_type&, boost::exponential_distribution<double> > waitGenerator(rng, waitDistrib);

	MajorantType maj;
	TermType term = (TermType)iterm;

	// Select properties by which to choose particles (-1 means
	// choose uniformly).  Note we need to choose 2 particles.  There
	// are six possible rate terms to choose from; 4 slip-flow and 2
	// free molecular.
	if ((sys.GetIncepted() + sys.ParticleCount()) < 2)  // if there are < 2 SPs but incepting class has weight >= 2, we can still act
		return 1;
	else
		Select_ip12(t, sys, iterm, rng, maj, ip1, ip2);                              // particle selection moved to its own function for neatness

	// Choose and get first particle.
	Particle *sp1 = NULL;

	// Is this an incepting class particle?
	if (hybrid_flag && ip1 <= -2)
	{
		sp1 = (ip1 == -2) ? sys.Particles().GetInceptedSP().Clone() : sys.Particles().GetInceptedSP_tmp().Clone();//m_mech->CreateParticle(t);
		sp1->setStatisticalWeight(1);
		//m_mech->Inceptions()[0]->Perform_incepted(t, sys, local_geom, 0, rng, *sp1); // Incept a new particle from SP[0]
		sys.AdjustIncepted(-(sp1->getStatisticalWeight()));                          // Reduce the incepting class count
		ip1_flag = true;                                                             // Flag sp1 as an incepting class particle
		sys.AdjustInceptingCoagulations();                                           // Increment number of times particles have left the incepting class

		//sp_age = waitGenerator();                                                    // Choose an LPDA last update time from the class residence time distribution
		//if (sp_age > t - create_t)                                                   // t \in [tcreate,t], age \in [0,t-tcreate]
		//	sp_age = t - create_t;
		//sp1->SetTime(t - sp_age);

		unsigned int m = sp1->Composition()[0] - 2; //added
		if (m > 0)                                                                   // Adjust the gas phase to account for the surface growth performed
		{
			adjustGas(sys, sp1->getStatisticalWeight(), m);
			adjustParticleTemperature(sys, sp1->getStatisticalWeight(), m, sys.GetIsAdiabaticFlag(), 1, 2);
		}

		if (ip1 == -2)
			sp1->SetTime(t);

		ip1 = sys.Particles().Add(*sp1, rng);                                        // Add the particle to the ensemble
	}
	else
	{
		if (ip1 >= 0) {
			sp1 = sys.Particles().At(ip1);
			//double just = waitGenerator(); // temporary
			//double just2 = waitGenerator(); // temporary
		}
		else {
			// Failed to choose a particle.
			return -1;
		}
	}

	Particle *sp2 = NULL;

	// Is this an incepting class particle?
	if (hybrid_flag && ip2 <= -2)
	{
		sp2 = (ip2 == -2) ? sys.Particles().GetInceptedSP().Clone() : sys.Particles().GetInceptedSP_tmp().Clone();// m_mech->CreateParticle(t);
		sp2->setStatisticalWeight(1);
		//m_mech->Inceptions()[0]->Perform_incepted(t, sys, local_geom, 0, rng, *sp2); // Incept a new particle from SP[0]
		// Note we do not need to add it to the ensemble
		ip2_flag = true;                                                             // Flag sp2 as an incepting class particle

		//sp_age = waitGenerator();                                                    // Choose an LPDA last update time from the class residence time distribution
		//if (sp_age > t - create_t)                                                   // t \in [tcreate,t], age \in [0,t-tcreate]
		//	sp_age = t - create_t;
		//sp2->SetTime(t - sp_age);

		if (ip2 == -2)
			sp2->SetTime(t);

		ip2 = -2;
	}
	else
	{
		if ((ip2 >= 0) && (ip2 != ip1)) {
			sp2 = sys.Particles().At(ip2);
			//double just = waitGenerator(); // temporary
			//double just2 = waitGenerator(); // temporary
		}
		else {
			// Failed to select a unique particle.
			return -1;
		}
	}


	/*std::ofstream pscFile;
	std::string pscfname;
	pscfname = "sp_update_times.csv";
	if (t > 0.0015 && t <= 0.001505)
	{
		unsigned int counter = 0;
		for (unsigned int i = 0; i < sys.ParticleCount(); ++i)
		{
			if (sys.Particles().At(i)->CollDiameter() < 5e-10)
				++counter;
		}
		counter += sys.GetIncepted();

		pscFile.open(pscfname.c_str(), std::ios::app);
		// t maj ip1 ip1flag ip2 ip2flag sp1_wt sp2_wt sp1_d sp2_d sp1_age sp2_age sp1_lut sp2_lut coag successful
		pscFile << t << " , " << term << " , "
			<< sys.ParticleCount() << " , " << sys.GetIncepted() << " , " << counter << " , "
			<< ip1_flag << " , " << ip1 << " , " << ip2_flag << " , " << ip2 << " , "
			<< sp1->getStatisticalWeight() << " , " << sp2->getStatisticalWeight() << " , "
			<< sp1->Primary()->CollDiameter() << " , " << sp2->Primary()->CollDiameter() << " , "
			<< sp1->CreateTime() << " , " << sp2->CreateTime() << " , "
			<< sp1->LastUpdateTime() << " , " << sp2->LastUpdateTime() << " , ";
		pscFile.close();
	}*/

	//Calculate the majorant rate before updating the particles
	double majk = MajorantKernel(*sp1, *sp2, sys, maj);

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

	// If sp2 does not belong to the ensemble, do not adjust the gas phase in surface growth
	// so sp2 can be returned to the bin unchanged if the coagulation is unsuccessful
	//if (ip2_flag)
	//	sys.SetNotPSIFlag(false);
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
			}
			else {
				ceff = 1;
			}
			truek *= ceff;
		}

		if (!Fictitious(majk, truek, rng)) {
			// If particle sp2 is used, we need to remove formally it from the incepting class
			if (ip2_flag)
			{
				sys.AdjustIncepted(-(sp2->getStatisticalWeight()));
				sys.AdjustInceptingCoagulations();                                           // Increment number of times particles have left the incepting class
			}
			JoinParticles(t, ip1, sp1, ip2, sp2, sys, rng);
			if (ip2_flag && sp2 != NULL)
			{
				unsigned int m = sp2->Composition()[0] - 2;
				if (m > 0)                                                                   // Adjust the gas phase to account for the surface growth performed
				{
					adjustGas(sys, sp2->getStatisticalWeight(), m);
					adjustParticleTemperature(sys, sp2->getStatisticalWeight(), m, sys.GetIsAdiabaticFlag(), 1, 2);
				}
				//sys.SetNotPSIFlag(true);
				// Particle sp2 is not in the ensemble, must manually delete it
				delete sp2;
				sp2 = NULL;
			}
			/*if (t > 0.0015 && t <= 0.001505)
			{
				pscFile.open(pscfname.c_str(), std::ios::app);
				pscFile << 1 << "\n";
				pscFile.close();
			}*/
		}
		else
		{
			sys.Particles().Update(ip1);
			if (!ip2_flag)
				sys.Particles().Update(ip2);
			else if (sp2 != NULL)
			{
				//double ncomp = sp2->Composition()[0];
				//double adjust = (ncomp - 2) / 2;
				//sys.AdjustIncepted(adjust);
				//sys.AdjustInceptions(adjust);
				// Particle sp2 is not in the ensemble, must manually delete it
				delete sp2;
				sp2 = NULL;
				//sys.SetNotPSIFlag(true);
			}
			/*if (t > 0.0015 && t <= 0.001505)
			{
				pscFile.open(pscfname.c_str(), std::ios::app);
				pscFile << 0 << "\n";
				pscFile.close();
			}*/
			return 1; // Ficticious event.
		}
	}
	else
	{
		// One or both particles were invalidated on update,
		// but that's not a problem.  Information on the update
		// of valid particles must be propagated into the binary
		// tree
		if (ip1 != -1)
			sys.Particles().Update(ip1);

		if (ip2 != -1)
		{
			if (!ip2_flag)
				sys.Particles().Update(ip2);
			//else
			//{
			//double ncomp = sp2->Composition()[0];
			//double adjust = (ncomp - 2) / 2;
			//sys.AdjustIncepted(adjust);
			//sys.AdjustInceptions(adjust);
			//}
		}

		if (ip2_flag && sp2 != NULL)
		{
			// Particle sp2 is not in the ensemble, must manually delete it
			delete sp2;
			sp2 = NULL;
		}
		/*if (t > 0.0039 && t < 0.0052)
		{
			pscFile << -1 << "\n";
			pscFile.close();
		}*/
	}

	if (ip2_flag && sp2 != NULL)
	{
		// Particle sp2 is not in the ensemble, must manually delete it
		delete sp2;
		sp2 = NULL;
	}

	// make sure surface growth is updating the gas phase again
	/*if (ip2_flag)
		sys.SetNotPSIFlag(true);*/

	return 0;
}

void TransitionCoagulation::ChooseProps(Sweep::Cell &sys, unsigned int iterm) const
{
	Sweep::PropID prop1, prop2;

	// Properties to which the probabilities of particle selection will be proportional
	switch (static_cast<TermType>(iterm)) {
	case FreeMol1:
		prop1 = iUniform;
		prop2 = iD2_M_1_2;
		break;
	case FreeMol2:
		prop1 = iD2;
		prop2 = iM_1_2;
		break;
	case SlipFlow1:
		prop1 = iUniform;
		prop2 = iUniform;
		break;
	case SlipFlow2:
		prop1 = iDcol;
		prop2 = iD_1;
		break;
	case SlipFlow3:
		prop1 = iUniform;
		prop2 = iD_1;
		break;
	case SlipFlow4:
		prop1 = iDcol;
		prop2 = iD_2;
		break;
	default:
		throw std::logic_error("Unrecognised term, (Sweep, TransitionCoagulation::Perform)");
	}

	sys.SetCoagProps(prop1, prop2);
}

// Select properties by which to choose particles (-1 means
// choose uniformly).  Note we need to choose 2 particles.  There
// are six possible rate terms to choose from; 4 slip-flow and 2
// free molecular.
void Sweep::Processes::TransitionCoagulation::Select_ip12(double t,
	Sweep::Cell &sys,
	unsigned int iterm,
	Sweep::rng_type &rng,
	mutable Sweep::Processes::TransitionCoagulation::MajorantType &maj,
	mutable int &ip1,
	mutable int &ip2) const
{
	boost::uniform_01<rng_type&, double> unifDistrib(rng); // Select the first particle and note the majorant type.

	bool hybrid_flag = m_mech->IsHybrid() && sys.GetIncepted() > 0 ;

	double wt_incfrac = 0, d_incfrac = 0, d2_incfrac = 0, d_1_incfrac = 0, d_2_incfrac = 0, m_1_2_incfrac = 0, d2m_1_2_incfrac = 0;
	double wt_tmp1frac = 0, d_tmp1frac = 0, d2_tmp1frac = 0, d_1_tmp1frac = 0, d_2_tmp1frac = 0, m_1_2_tmp1frac = 0, d2m_1_2_tmp1frac = 0;

	TermType term = (TermType)iterm;

	double wt_incep = 0, dc_incep = 0, sqrt_m_incep = 0;
	double wt_tmp1 = 0, dc_tmp1 = 0, sqrt_m_tmp1 = 0;//

	if (hybrid_flag)
	{
		ip1 = -2;
		//wt_incep = sys.Particles().GetInceptedSP_tmp().getStatisticalWeight();
		//dc_incep = sys.Particles().GetInceptedSP_tmp().CollDiameter();
		//sqrt_m_incep = std::sqrt(sys.Particles().GetInceptedSP_tmp().Mass());
		wt_incep = sys.GetInceptions_tmp();
		dc_incep = sys.Particles().GetInceptedSP().CollDiameter();
		sqrt_m_incep = std::sqrt(sys.Particles().GetInceptedSP().Mass());
		wt_tmp1 = sys.GetIncepted() - wt_incep;
		wt_tmp1 = (wt_tmp1 > 0) ? wt_tmp1 : 0.0;
		dc_tmp1 = sys.Particles().GetInceptedSP_tmp().CollDiameter();
		sqrt_m_tmp1 = std::sqrt(sys.Particles().GetInceptedSP_tmp().Mass());
	}
	double test = unifDistrib();

	switch (term) {
	case SlipFlow1:
		if (hybrid_flag)
		{
			double wt_other = sys.Particles().Count();
			wt_incfrac = wt_incep / (wt_other + wt_incep + wt_tmp1);
			wt_tmp1frac = wt_tmp1 / (wt_other + wt_tmp1);
		}
		if (test > wt_incfrac) //!hybrid_flag || 
		{
			test = unifDistrib();
			if (test > wt_tmp1)
				ip1 = sys.Particles().Select(rng);
			else
				ip1 = -3;
		}
		maj = SlipFlow;
		break;
	case SlipFlow2:
		if (hybrid_flag)
		{
			double d_incep = wt_incep * dc_incep;
			double d_tmp1 = wt_tmp1 * dc_tmp1;
			double d_other = sys.Particles().GetSum(iDcol);
			d_incfrac = d_incep / (d_other + d_incep + d_tmp1);
			d_tmp1frac = d_tmp1 / (d_other + d_tmp1);
		}
		if (test > d_incfrac)
		{
			test = unifDistrib();
			if (test > d_tmp1frac)
				ip1 = sys.Particles().Select(iDcol, rng);
			else
				ip1 = -3;
		}
		maj = SlipFlow;
		break;
	case SlipFlow3:
		if (hybrid_flag)
		{
			double wt_other = sys.Particles().Count();
			wt_incfrac = wt_incep / (wt_other + wt_incep + wt_tmp1);
			wt_tmp1frac = wt_tmp1 / (wt_other + wt_tmp1);
		}
		if (test > wt_incfrac)
		{
			test = unifDistrib();
			if (test > wt_tmp1frac)
				ip1 = sys.Particles().Select(rng);
			else
				ip1 = -3;
		}
		maj = SlipFlow;
		break;
	case SlipFlow4:
		if (hybrid_flag)
		{
			double d_incep = wt_incep * dc_incep;
			double d_tmp1 = wt_tmp1 * dc_tmp1;
			double d_other = sys.Particles().GetSum(iDcol);
			d_incfrac = d_incep / (d_other + d_incep + d_tmp1);
			d_tmp1frac = d_tmp1 / (d_other + d_tmp1);
		}
		if (test > d_incfrac)
		{
			test = unifDistrib();
			if (test > d_tmp1frac)
				ip1 = sys.Particles().Select(iDcol, rng);
			else
				ip1 = -3;
		}
		maj = SlipFlow;
		break;
	case FreeMol1:
		if (hybrid_flag)
		{
			double wt_other = sys.Particles().Count();
			wt_incfrac = wt_incep / (wt_other + wt_incep + wt_tmp1);
			wt_tmp1frac = wt_tmp1 / (wt_other + wt_tmp1);
		}
		if (test > wt_incfrac)
		{
			test = unifDistrib();
			if (test > wt_tmp1frac)
				ip1 = sys.Particles().Select(rng);
			else
				ip1 = -3;
		}
		maj = FreeMol;
		break;
	case FreeMol2:
		if (hybrid_flag)
		{
			double d2_incep = wt_incep * dc_incep * dc_incep;
			double d2_tmp1 = wt_tmp1 * dc_tmp1 * dc_tmp1;
			double d2_other = sys.Particles().GetSum(iD2);
			d2_incfrac = d2_incep / (d2_other + d2_incep + d2_tmp1);
			d2_tmp1frac = d2_tmp1 / (d2_other + d2_tmp1);
		}
		if (test > d2_incfrac)
		{
			test = unifDistrib();
			if (test > d2_tmp1frac)
				ip1 = sys.Particles().Select(iD2, rng);
			else
				ip1 = -3;
		}
		maj = FreeMol;
		break;
	default:
		if (hybrid_flag)
		{
			double wt_other = sys.Particles().Count();
			wt_incfrac = wt_incep / (wt_other + wt_incep + wt_tmp1);
			wt_tmp1frac = wt_tmp1 / (wt_other + wt_tmp1);
		}
		if (test > wt_incfrac)
		{
			test = unifDistrib();
			if (test > wt_tmp1frac)
				ip1 = sys.Particles().Select(rng);
			else
				ip1 = -3;
		}
		maj = SlipFlow;
		break;
	}

	// Choose and get unique second particle.  Note, we are allowed to do
	// this even if the first particle was invalidated.
	ip2 = ip1;
	unsigned int guard = 0;
	test = unifDistrib();
	bool mustSwitch = (ip1 < -1) && (sys.GetIncepted() <= 1);

	switch (term) {
	case SlipFlow1:
		if (hybrid_flag)
		{
			double wt_other = sys.Particles().Count();
			wt_incfrac = wt_incep / (wt_other + wt_incep + wt_tmp1);
			wt_tmp1frac = wt_tmp1 / (wt_other + wt_tmp1);
		}
		if (test > wt_incfrac || mustSwitch)
		{
			test = unifDistrib();
			if (test > wt_tmp1frac || mustSwitch)
			{
				while ((ip2 == ip1) && (++guard < 1000))
					ip2 = sys.Particles().Select(rng);
			}
			else
				ip2 = -3;
		}
		else
		{
			ip2 = -2;
		}
		break;
	case SlipFlow2:
		if (hybrid_flag)
		{
			double d_1_incep = wt_incep / dc_incep;
			double d_1_tmp1 = wt_tmp1 / dc_tmp1;
			double d_1_other = sys.Particles().GetSum(iD_1);
			d_1_incfrac = d_1_incep / (d_1_other + d_1_incep + d_1_tmp1);
			d_1_tmp1frac = d_1_tmp1 / (d_1_other + d_1_tmp1);
		}
		if (test > d_1_incfrac || mustSwitch)
		{
			test = unifDistrib();
			if (test > d_1_tmp1frac || mustSwitch)
			{
				while ((ip2 == ip1) && (++guard < 1000))
					ip2 = sys.Particles().Select(iD_1, rng);
			}
			else
				ip2 = -3;
		}
		else
		{
			ip2 = -2;
		}
		break;
	case SlipFlow3:
		if (hybrid_flag)
		{
			double d_1_incep = wt_incep / dc_incep;
			double d_1_tmp1 = wt_tmp1 / dc_tmp1;
			double d_1_other = sys.Particles().GetSum(iD_1);
			d_1_incfrac = d_1_incep / (d_1_other + d_1_incep + d_1_tmp1);
			d_1_tmp1frac = d_1_tmp1 / (d_1_other + d_1_tmp1);
		}
		if (test > d_1_incfrac || mustSwitch)
		{
			test = unifDistrib();
			if (test > d_1_tmp1frac || mustSwitch)
			{
				while ((ip2 == ip1) && (++guard < 1000))
					ip2 = sys.Particles().Select(iD_1, rng);
			}
			else
				ip2 = -3;
		}
		else
		{
			ip2 = -2;
		}
		break;
	case SlipFlow4:
		if (hybrid_flag)
		{
			double d_2_incep = wt_incep / (dc_incep * dc_incep);
			double d_2_tmp1 = wt_tmp1 / (dc_tmp1 * dc_tmp1);
			double d_2_other = sys.Particles().GetSum(iD_2);
			d_2_incfrac = d_2_incep / (d_2_other + d_2_incep + d_2_tmp1);
			d_2_tmp1frac = d_2_tmp1 / (d_2_other + d_2_tmp1);
		}
		if (test > d_2_incfrac || mustSwitch)
		{
			test = unifDistrib();
			if (test > d_2_tmp1frac || mustSwitch)
			{
				while ((ip2 == ip1) && (++guard < 1000))
					ip2 = sys.Particles().Select(iD_2, rng);
			}
			else
				ip2 = -3;
		}
		else
		{
			ip2 = -2;
		}
		break;
	case FreeMol1:
		if (hybrid_flag)
		{
			double d2m_1_2_incep = wt_incep * dc_incep * dc_incep / sqrt_m_incep;
			double d2m_1_2_tmp1 = wt_tmp1 * dc_tmp1 * dc_tmp1 / sqrt_m_tmp1;
			double d2m_1_2_other = sys.Particles().GetSum(iD2_M_1_2);
			d2m_1_2_incfrac = d2m_1_2_incep / (d2m_1_2_other + d2m_1_2_incep + d2m_1_2_tmp1);
			d2m_1_2_tmp1frac = d2m_1_2_tmp1 / (d2m_1_2_other + d2m_1_2_tmp1);
		}
		if (test > d2m_1_2_incfrac || mustSwitch)
		{
			test = unifDistrib();
			if (test > d2m_1_2_tmp1frac || mustSwitch)
			{
				while ((ip2 == ip1) && (++guard < 1000))
					ip2 = sys.Particles().Select(iD2_M_1_2, rng);
			}
			else
				ip2 = -3;
		}
		else
		{
			ip2 = -2;
		}
		break;
	case FreeMol2:
		if (hybrid_flag)
		{
			double m_1_2_incep = wt_incep / sqrt_m_incep;
			double m_1_2_tmp1 = wt_tmp1 / sqrt_m_tmp1;
			double m_1_2_other = sys.Particles().GetSum(iM_1_2);
			m_1_2_incfrac = m_1_2_incep / (m_1_2_other + m_1_2_incep + m_1_2_tmp1);
			m_1_2_tmp1frac = m_1_2_tmp1 / (m_1_2_other + m_1_2_tmp1);
		}
		if (test > m_1_2_incfrac || mustSwitch)
		{
			test = unifDistrib();
			if (test > m_1_2_tmp1frac || mustSwitch)
			{
				while ((ip2 == ip1) && (++guard < 1000))
					ip2 = sys.Particles().Select(iM_1_2, rng);
			}
			else
				ip2 = -3;
		}
		else
		{
			ip2 = -2;
		}
		break;
	default:
		if (hybrid_flag)
		{
			double wt_other = sys.Particles().Count();
			wt_incfrac = wt_incep / (wt_other + wt_incep + wt_tmp1);
			wt_tmp1frac = wt_tmp1 / (wt_other + wt_tmp1);
		}
		if (test > wt_incfrac || mustSwitch)
		{
			test = unifDistrib();
			if (test > wt_tmp1frac || mustSwitch)
			{
				while ((ip2 == ip1) && (++guard < 1000))
					ip2 = sys.Particles().Select(rng);
			}
			else
				ip2 = -3;
		}
		else
		{
			ip2 = -2;
		}
		break;
	}
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
double Sweep::Processes::TransitionCoagulation::CoagKernel(const Particle &sp1,
	const Particle &sp2,
	const Cell &sys) const
{
	const double T = sys.GasPhase().Temperature();
	const double P = sys.GasPhase().Pressure();
	const double fm = FreeMolKernel(sp1, sp2, T, P, false);
	const double sf = SlipFlowKernel(sp1, sp2, T, P, sys.GasPhase().Viscosity(), false);
	return (fm*sf) / (fm + sf);
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
double Sweep::Processes::TransitionCoagulation::MajorantKernel(const Particle &sp1,
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
double Sweep::Processes::TransitionCoagulation::FreeMolKernel(const Particle &sp1, const Particle &sp2,
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
	}
	else {
		const double dterm = d1 + d2;
		return m_efm * CFM * A() *
			sqrt(T * (invm1 + invm2)) *
			dterm * dterm;
	}
}

// Returns the slip-flow coagulation kernel value for the
// two given particles.  Can return either the majorant or
// true kernel.
double Sweep::Processes::TransitionCoagulation::SlipFlowKernel(const Particle &sp1, const Particle &sp2,
	double T, double P, double mu, bool maj) const
{
	// Collect the particle properties
	const double d1 = sp1.CollDiameter();
	const double d2 = sp2.CollDiameter();

	// For the slip-flow kernel the majorant and non-majorant forms are identical.
	return ((1.257 * 2.0 * MeanFreePathAir(T, P) *
		(1.0 / d1 / d1 + 1.0 / d2 / d2)) +
		(1.0 / d1 + 1.0 / d2)) *
		CSF * T * (d1 + d2)
		* A() / mu;
}

