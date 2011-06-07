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

using namespace Sweep::Processes;

//* Free-molecular enhancement factor.
const Sweep::real Sweep::Processes::TransitionCoagulation::m_efm = 2.2; // 2.2 is for soot.


// Default constructor.
Sweep::Processes::TransitionCoagulation::TransitionCoagulation(const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "TransitionRegimeCoagulation";
}

Sweep::Processes::TransitionCoagulation* const Sweep::Processes::TransitionCoagulation::Clone() const
{
    return new TransitionCoagulation(*this);
}

// Stream-reading constructor.
Sweep::Processes::TransitionCoagulation::TransitionCoagulation(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech)
{
    m_name = "TransitionRegimeCoagulation";
    Deserialize(in, mech);
}

// TOTAL RATE CALCULATION.

// Returns the rate of the process for the given system.
Sweep::real Sweep::Processes::TransitionCoagulation::Rate(real t, const Cell &sys) const
{
    // Get the number of particles in the system.
    unsigned int n = sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        real T = sys.Temperature();
        real P = sys.Pressure();

        // Calculate the rate.
        return Rate(sys.Particles().GetSums(), (real)n, sqrt(T),
                    T/ViscosityAir(T), MeanFreePathAir(T,P),
                    sys.SampleVolume());
    } else {
        return 0.0;
    }
}

// More efficient rate routine for coagulation only.
// All parameters required to calculate rate passed
// as arguments.
Sweep::real Sweep::Processes::TransitionCoagulation::Rate(const TreeWeightedCache &data, real n, real sqrtT,
                       real T_mu, real MFP, real vol) const
{
    // Some prerequisites.
    real n_1 = n - 1.0;
    real a = CSF * T_mu * A();
    real b = a * MFP * 1.257 * A();
    real c = CFMMAJ * m_efm * CFM * sqrtT * A();

    // Summed particle properties required for coagulation rate.
    const real d       = data.Property(Sweep::iDcol);
    const real d2      = data.Property(Sweep::iD2);
    const real d_1     = data.Property(Sweep::iD_1);
    const real d_2     = data.Property(Sweep::iD_2);
    const real m_1_2   = data.Property(Sweep::iM_1_2);
    const real d2m_1_2 = data.Property(Sweep::iD2_M_1_2);

    // Get individual terms.
    real terms[TYPE_COUNT];
    // Slip-flow.
    terms[0] = n * n_1 * a / vol;
    terms[1] = ((d * d_1) - n) * a / vol;
    terms[2] = d_1 * n_1 * b / vol;
    terms[3] = ((d * d_2) - d_1) * b / vol;
    // Free-molecular.
    terms[4] = n_1 * d2m_1_2  * c / vol;
    terms[5] = (m_1_2 * d2 - d2m_1_2) * c / vol;

    // Sum up total coagulation rates for different regimes.
    real sf = terms[0] + terms[1] + terms[2] + terms[3];
    real fm = terms[4] + terms[5];

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
unsigned int Sweep::Processes::TransitionCoagulation::TermCount(void) const {return TYPE_COUNT;}

// Calculates the rate terms given an iterator to a real vector. The
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all rate terms.
Sweep::real Sweep::Processes::TransitionCoagulation::RateTerms(real t, const Cell &sys,
                            fvector::iterator &iterm) const
{
    // Get the number of particles in the system.
    unsigned int n = sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        real T = sys.Temperature();
        real P = sys.Pressure();

        // Calculate the rate terms.
        return RateTerms(sys.Particles().GetSums(), (real)n, sqrt(T), T/ViscosityAir(T),
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
Sweep::real Sweep::Processes::TransitionCoagulation::RateTerms(const TreeWeightedCache &data, real n, real sqrtT,
                            real T_mu, real MFP, real vol,
                            fvector::iterator &iterm) const
{
    // Some prerequisites.
    real n_1 = n - 1.0;
    real a   = CSF * T_mu * A();
    real b   = a * MFP * 1.257 * 2.0 * A();
    real c   = CFMMAJ * m_efm * CFM * sqrtT * A();

    // Summed particle properties required for coagulation rate.
    const real d       = data.Property(Sweep::iDcol);
    const real d2      = data.Property(Sweep::iD2);
    const real d_1     = data.Property(Sweep::iD_1);
    const real d_2     = data.Property(Sweep::iD_2);
    const real m_1_2   = data.Property(Sweep::iM_1_2);
    const real d2m_1_2 = data.Property(Sweep::iD2_M_1_2);

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
    real sf = *(isf) + *(isf+1) + *(isf+2) + *(isf+3);
    real fm = *(ifm) + *(ifm+1);

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
 * \param[in,out]   rand_int    Pointer to function that generates uniform integers on a range
 * \param[in,out]   rand_u01    Pointer to function that generates U[0,1] deviates
 * \param[out]      out         Details of any particle being transported out of system
 *
 * \return      0 on success, otherwise negative.
 */
int TransitionCoagulation::Perform(Sweep::real t, Sweep::Cell &sys, 
                                   const Geometry::LocalGeometry1d& local_geom,
                                   unsigned int iterm,
                                   int (*rand_int)(int, int), 
                                   Sweep::real(*rand_u01)(), 
                                   Sweep::Transport::TransportOutflow *out) const
{
    // Select properties by which to choose particles (-1 means
    // choose uniformly).  Note we need to choose 2 particles.  There
    // are six possible rate terms to choose from; 4 slip-flow and 2
    // free molecular.
    if (sys.ParticleCount() < 2) {
        return 1;
    }

    int ip1=-1, ip2=-1;
    MajorantType maj;
    TermType term = (TermType)iterm;

    // Select the first particle and note the majorant type.
    switch (term) {
        case SlipFlow1:
            ip1 = sys.Particles().Select(rand_int);
            maj = SlipFlow;
            break;
        case SlipFlow2:
            ip1 = sys.Particles().Select(Sweep::iDcol, rand_int, rand_u01);
            maj = SlipFlow;
            break;
        case SlipFlow3:
            ip1 = sys.Particles().Select(rand_int);
            maj = SlipFlow;
            break;
        case SlipFlow4:
            ip1 = sys.Particles().Select(Sweep::iDcol, rand_int, rand_u01);
            maj = SlipFlow;
            break;
        case FreeMol1:
            ip1 = sys.Particles().Select(rand_int);
            maj = FreeMol;
            break;
        case FreeMol2:
            ip1 = sys.Particles().Select(Sweep::iD2, rand_int, rand_u01);
            maj = FreeMol;
            break;
        default :
            ip1 = sys.Particles().Select(rand_int);
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
    unsigned int guard = 0;
    switch (term) {
        case SlipFlow1:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(rand_int);
            break;
        case SlipFlow2:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(Sweep::iD_1, rand_int, rand_u01);
            break;
        case SlipFlow3:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(Sweep::iD_1, rand_int, rand_u01);
            break;
        case SlipFlow4:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(Sweep::iD_2, rand_int, rand_u01);
            break;
        case FreeMol1:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(Sweep::iD2_M_1_2, rand_int, rand_u01);
            break;
        case FreeMol2:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(Sweep::iM_1_2, rand_int, rand_u01);
            break;
        default :
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(rand_int);
            break;
    }


    Particle *sp2=NULL;
    if ((ip2>=0) && (ip2!=ip1)) {
        sp2 = sys.Particles().At(ip2);
    } else {
        // Failed to select a unique particle.
        return -1;
    }

    //Calculate the majorant rate before updating the particles
    real majk = MajorantKernel(*sp1, *sp2, sys, maj);

    //Update the particles
    m_mech->UpdateParticle(*sp1, sys, t, rand_u01);
    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp1->IsValid()) {
        // Must remove first particle now.
        sys.Particles().Remove(ip1);

        // Invalidating the index tells this routine not to perform coagulation.
        ip1 = -1;
        return 0;
    }

    m_mech->UpdateParticle(*sp2, sys, t, rand_u01);
    // Check validity of particles after update.
    if (!sp2->IsValid()) {
        // Tell the ensemble to update particle one before we confuse things
        // by removing particle 2
        sys.Particles().Update(ip1);

        // Must remove second particle now.
        sys.Particles().Remove(ip2);

        // Invalidating the index tells this routine not to perform coagulation.
        ip2 = -1;

        return 0;
    }


    // Check that both the particles are still valid.
    if ((ip1>=0) && (ip2>=0)) {
        // Must check for ficticious event now by comparing the original
        // majorant rate and the current (after updates) true rate.

        real truek = CoagKernel(*sp1, *sp2, sys);
		double ceff=0;
		if (majk<truek)
			std::cout << "maj< true"<< std::endl;

		//added by ms785 to include the collision efficiency in the calculation of the rate
		if (sys.ParticleModel()->AggModel()==AggModels::PAH_ID||sys.ParticleModel()->AggModel()==AggModels::PAH_KMC_ID)
		{
			 ceff=sys.ParticleModel()->CollisionEff(sp1,sp2);
			 truek*=ceff;
		}



        if (!Fictitious(majk, truek, rand_u01)) {
            JoinParticles(t, ip1, sp1, ip2, sp2, sys, rand_int, rand_u01);
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
Sweep::real Sweep::Processes::TransitionCoagulation::CoagKernel(const Particle &sp1,
                                                                const Particle &sp2,
                                                                const Cell &sys) const
{
    const real T = sys.Temperature();
    const real P = sys.Pressure();
    const real fm = FreeMolKernel(sp1, sp2, T, P, false);
    const real sf = SlipFlowKernel(sp1, sp2, T, P, false);
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
Sweep::real Sweep::Processes::TransitionCoagulation::MajorantKernel(const Particle &sp1,
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
            return FreeMolKernel(sp1, sp2, sys.Temperature(), sys.Pressure(), true);
        case SlipFlow:
            // Slip-flow majorant.
            return SlipFlowKernel(sp1, sp2, sys.Temperature(), sys.Pressure(), true);
    }

    // Invalid majorant, return zero.
    return 0.0;
}
// Returns the free-molecular coagulation kernel value for the
// two given particles.  Can return either the majorant or
// true kernel.
Sweep::real Sweep::Processes::TransitionCoagulation::FreeMolKernel(const Particle &sp1, const Particle &sp2,
                                real T, real P, bool maj) const
{
    // This routine calculate the free molecular coagulation kernel for two particles.
    // There are two forms of kernel; a majorant form and a non-majorant form.

    // Collect the particle properties
    const real d1 = sp1.CollDiameter();
    const real d2 = sp2.CollDiameter();
    const real invm1 = 1.0 / sp1.Mass();
    const real invm2 = 1.0 / sp2.Mass();

    if (maj) {
        // The majorant form is always >= the non-majorant form.
        return CFMMAJ * m_efm * CFM * sqrt(T) * A() *
               (std::sqrt(invm1) + std::sqrt(invm2)) *
               (d1 * d1 + d2 * d2);
    } else {
        const real dterm = d1 + d2;
        return m_efm * CFM * A() *
               sqrt(T * (invm1 + invm2)) *
               dterm * dterm;
    }
}

// Returns the slip-flow coagulation kernel value for the
// two given particles.  Can return either the majorant or
// true kernel.  This kernel is currently evaluated assuming
// that the surrounding gas is air, i.e. principally N2.
Sweep::real Sweep::Processes::TransitionCoagulation::SlipFlowKernel(const Particle &sp1, const Particle &sp2,
                                 real T, real P, bool maj) const
{
    // Collect the particle properties
    const real d1 = sp1.CollDiameter();
    const real d2 = sp2.CollDiameter();

    // For the slip-flow kernel the majorant and non-majorant forms are identical.
    return ((1.257 * 2.0 * MeanFreePathAir(T,P) *
             (1.0 / d1 / d1 + 1.0 / d2 / d2)) +
            (1.0 / d1 + 1.0 / d2)) *
           CSF * T * (d1 + d2)
           * A() / ViscosityAir(T);
}

