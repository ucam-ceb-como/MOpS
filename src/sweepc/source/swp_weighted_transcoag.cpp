/*!
 * \file   swp_weighted_transcoag.cpp
 * \author William J Menz 
 *  Copyright (C) 2011 William J Menz.
 *
 *  Project:        sweepc (population balance solver)
 *  Sourceforge:    http://sourceforge.net/projects/mopssuite
 * 
 * \brief Implementation of weighted transition regime coagulation kernel
 * 
  Author(s):      William J Menz
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2009 Robert I A Patterson.
   updated from DSA transition regime to weighted by William J Menz 2011

  File purpose:
    Implementation of weighted transition regime coagulation kernel

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
#include "swp_weighted_transcoag.h"
#include "swp_params.h"
#include "swp_cell.h"
#include "swp_mechanism.h"
#include <numeric>
using namespace std;
using namespace Sweep::Processes;


/**
 * Main way of building a new coagulation object.
 * 
 * \param[in]   mech            Mechanism to which coagulation will belong
 * \param[in]   weight_rule     Specify how to calculate statistical weight of newly coagulation particles
 *
 */
Sweep::Processes::WeightedTransitionCoagulation::WeightedTransitionCoagulation(
        const Sweep::Mechanism &mech,
        const CoagWeightRule weight_rule)
: Coagulation(mech)
, m_efm(mech.GetEnhancementFM())
, m_CoagWeightRule(weight_rule)
{
    m_name = "WeightedTransitionRegimeCoagulation";
}

// Clone object
Sweep::Processes::WeightedTransitionCoagulation* const Sweep::Processes::WeightedTransitionCoagulation::Clone() const
{
    return new WeightedTransitionCoagulation(*this);
}

/**
 * Load an instance of this process from a binary stream
 *
 * \param[in,out]   in      Input stream
 * \param[in]       mech    Mechanism to which process will belong
 *
 * \exception   runtime_error   Input stream not ready
 */
Sweep::Processes::WeightedTransitionCoagulation::WeightedTransitionCoagulation(std::istream &in, const Sweep::Mechanism &mech)
: Coagulation(mech)
, m_efm(mech.GetEnhancementFM())
{
    m_name = "WeightedTransitionCoagulation";
    Deserialize(in, mech);

    if(in.good())
        in.read(reinterpret_cast<char*>(&m_CoagWeightRule), sizeof(m_CoagWeightRule));
    else
        throw std::runtime_error("Input stream not ready (WeightedTransitionCoagulation::WeightedTransitionCoagulation)");
}

// TOTAL RATE CALCULATION.

/**
 * Calculate the total coagulation rate. Public function which calls the
 * helper function RateTerms to calculate the individual terms. Returns
 * the sum of all rate terms.
 *
 * \param[in] t             Time for which rates are requested
 * \param[in] sys           Details of the particle population and environment
 * \param[in] local_geom    Cell physical size and position information
 *
 * \return      Sum of all rate terms for this process
 */
double Sweep::Processes::WeightedTransitionCoagulation::Rate(double t, const Cell &sys,
                                                                  const Geometry::LocalGeometry1d &local_geom) const
{
    // Get the number of particles in the system.
    unsigned int n = sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        double T = sys.GasPhase().Temperature();
        double P = sys.GasPhase().Pressure();

        // Create a vector so we can call through to RateTerms
        Sweep::fvector vec(TYPE_COUNT);
        fvector::iterator it = vec.begin();

        return RateTerms(sys.Particles().GetSums(), (double)n, sqrt(T), T/sys.GasPhase().Viscosity(),
                MeanFreePathAir(T,P), sys.SampleVolume(), it);
    } else {
        // Not enough particles so return 0
        return 0.0;
    }
}

/**
 * Number of terms in the expression for the sum of the majorant
 * kernel over all particle pairs.
 */
unsigned int Sweep::Processes::WeightedTransitionCoagulation::TermCount() const {return TYPE_COUNT;}

/**
 * Publicly available calculation of weighted transition kernel rate terms.
 * 
 * Calculate the terms in the sum of the majorant kernel over all particle
 * pairs, placing each term in successive positions of the sequence
 * beginning at iterm and return the sum of the terms added to that
 * vector.
 *
 * \param[in]       t           Time for which rates are requested
 * \param[in]       sys         Details of the particle population and environment
 * \param[in]       local_geom  Cell physical size and position information
 * \param[in,out]   iterm       Pointer to start of sequence to hold the rate terms, returned as one past the end.
 *
 * \return      Sum of all rate terms for this process
 */
double Sweep::Processes::WeightedTransitionCoagulation::RateTerms(double t, const Cell &sys,
                                                                       const Geometry::LocalGeometry1d &local_geom,
                                                                       fvector::iterator &iterm) const
{
    // Get the number of particles in the system.
    unsigned int n = sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        double T = sys.GasPhase().Temperature();
        double P = sys.GasPhase().Pressure();

        double r = RateTerms(sys.Particles().GetSums(), (double)n, sqrt(T), T/sys.GasPhase().Viscosity(),
                MeanFreePathAir(T,P), sys.SampleVolume(), iterm);
        return r;

    } else {
        // Free-molecular.
    	fill(iterm, iterm+4, 0.0);
    	iterm = iterm+4;

        // Slip-flow.
        fill(iterm, iterm+7, 0.0);
    	iterm = iterm+7;
        return 0.0;
    }
}

/**
 * Calculate the individual rate terms for weighted transition coagulation.
 * 
 * Calculate the terms in the sum of the majorant kernel over all particle
 * pairs, placing each term in successive positions of the sequence
 * beginning at iterm and return the sum of the terms added to that
 * vector.
 *
 * \param[in] data      Particle data cache pointer
 * \param[in] n       	Number of particles
 * \param[in] sqrtT		Square-root of temperature
 * \param[in] T_mu		Temperature divided by viscosity
 * \param[in] MFP		Mean free path in air
 * \param[in] vol		Sample volume used for scaling
 * \param[in,out] iterm  Pointer to start of sequence to hold the rate terms, returned as one past the end.
 *
 * \return      Sum of all rate terms for this process
 */
double Sweep::Processes::WeightedTransitionCoagulation::RateTerms(
		const Ensemble::particle_cache_type &data, double n, double sqrtT,
        double T_mu, double MFP, double vol,
        fvector::iterator &iterm) const
{
    // Some prerequisites.
    double n_1 = n - 1.0;
    double a   = CSF * T_mu * A();
    double b   = a * MFP * 1.257 * 2.0;
    double c   = CFMMAJ * m_efm * CFM * sqrtT * A();

    // Summed particle properties required for coagulation rate.
    const double d        = data.Property(Sweep::iDcol);
    const double d2       = data.Property(Sweep::iD2);
    const double d_1      = data.Property(Sweep::iD_1);
    const double d_2      = data.Property(Sweep::iD_2);
    const double m_1_2    = data.Property(Sweep::iM_1_2);
    const double d2m_1_2  = data.Property(Sweep::iD2_M_1_2);

    // Weighted particle specific properties
    const double w        = data.Property(Sweep::iW);
    const double dw       = data.Property(Sweep::iDW);
    const double d2w	    = data.Property(Sweep::iD2W);
    const double d_1w     = data.Property(Sweep::iD_1W);
    const double d_2w     = data.Property(Sweep::iD_2W);
    const double m_1_2w   = data.Property(Sweep::iM_1_2W);
    const double d2m_1_2w = data.Property(Sweep::iD2_M_1_2W);

    fvector::iterator ifm = iterm;
    fvector::iterator isf = iterm+4;

    // Get individual terms

    // Free-molecular.
    *(iterm) =  n_1 * d2m_1_2w * c / vol;
    *(++iterm) = (d2 * m_1_2w - d2m_1_2w) * c / vol;
    *(++iterm) = (d2w * m_1_2 - d2m_1_2w) * c / vol;
    *(++iterm) = (d2m_1_2 * w - d2m_1_2w) * c / vol;

    // Slip-flow.
    *(++iterm) = 2 * n_1 * w * a / vol;
    *(++iterm) = (d * d_1w - w) * a / vol;
    *(++iterm) = (dw * d_1 - w) * a / vol;
    *(++iterm) = n_1 * d_1w * b / vol;
    *(++iterm) = (d * d_2w - d_1w) * b / vol;
    *(++iterm) = (dw * d_2 - d_1w) * b / vol;
    *(++iterm) = (d_1 * w - d_1w) * b / vol;

    // Return iterator to next term after the coagulation terms.
    ++iterm;

    // Sum up total coagulation rates for different regimes.
    double fm = *(ifm) + *(ifm+1) + *(ifm+2) + *(ifm+3);
    double sf = *(isf) + *(isf+1) + *(isf+2) + *(isf+3) + *(isf+4) + *(isf+5) + *(isf+6);
    //cout << "fm =" << fm << "sf = " << sf << endl;
    if ((sf>0.0) || (fm>0.0)) {
        // There is some coagulation.
        if (sf > fm) {
            // Use free-mol majorant.
        	fill(isf, isf+7, 0.0);
            return fm;
        } else {
            // Use slip-flow majorant.
            fill(ifm, ifm+4, 0.0);
            return sf;
        }
    } else {
        // Something went wrong with the rate calculation.
    	fill(isf, isf+7, 0.0);
    	fill(ifm, ifm+4, 0.0);
        return 0.0;
    }
}

/*!
 * Perform coagulation process.
 * Choose the properties of two particles for coagulation based on the index
 * passed to the function. Then call WeightedPerform to conduct the process.
 *
 * \param[in]       t           Time
 * \param[in,out]   sys         System to update
 * \param[in]       local_geom  Details of local phsyical layout
 * \param[in]       iterm       Process term responsible for this event
 * \param[in,out]   rng         Random number generator
 *
 * \return      0 on success, otherwise negative.
 * \exception   logic_error     Unrecognised weight rule
 * \exception   logic_error     Unrecognised rate term index (ie iterm value)
 */
int Sweep::Processes::WeightedTransitionCoagulation::Perform(
        double t, Sweep::Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng) const
{
    assert(iterm < TYPE_COUNT);

    // Select properties by which to choose particles.
    // Note we need to choose 2 particles.  One particle must be chosen
    // uniformly and one with probability proportional
    // to particle mass.

    if (sys.ParticleCount() < 2) {
        return 1;
    }

    MajorantType maj;
    Sweep::PropID prop1, prop2;

    // Properties to which the probabilities of particle selection will be proportional
    switch(static_cast<TermType>(iterm)) {
        case FreeMol1:
            prop1 = iUniform;
            prop2 = iD2_M_1_2W;
            maj = FreeMol;
            break;
        case FreeMol2:
            prop1 = iD2;
            prop2 = iM_1_2W;
            maj = FreeMol;
            break;
        case FreeMol3:
            prop1 = iM_1_2;
            prop2 = iD2W;
            maj = FreeMol;
            break;
        case FreeMol4:
            prop1 = iD2_M_1_2;
            prop2 = iW;
            maj = FreeMol;
            break;
        case SlipFlow1:
            prop1 = iUniform;
            prop2 = iW;
            maj = SlipFlow;
            break;
        case SlipFlow2:
            prop1 = iDcol;
            prop2 = iD_1W;
            maj = SlipFlow;
            break;
        case SlipFlow3:
            prop1 = iD_1;
            prop2 = iDW;
            maj = SlipFlow;
            break;
        case SlipFlow4:
            prop1 = iUniform;
            prop2 = iD_1W;
            maj = SlipFlow;
            break;
        case SlipFlow5:
            prop1 = iDcol;
            prop2 = iD_2W;
            maj = SlipFlow;
            break;
        case SlipFlow6:
            prop1 = iD_2;
            prop2 = iDW;
            maj = SlipFlow;
            break;
        case SlipFlow7:
            prop1 = iD_1;
            prop2 = iW;
            maj = SlipFlow;
            break;
        default :
        	throw std::logic_error("Unrecognised term, (Sweep, WeightedTransitionCoagulation::Perform)");
    }

    return WeightedPerform(t, prop1, prop2, m_CoagWeightRule, sys, rng, maj);
}


// COAGULATION KERNELS.

/**
 * Calculate the coagulation kernel between two particles in the given environment.
 * Note that the weights are included in FreeMolKernel and SlipFlowKernel rather than here.
 *
 *\param[in]    sp1         First particle
 *\param[in]    sp2         Second particle
 *\param[in]    sys         Details of the environment, including temperature and pressure
 *
 *\return       Value of kernel
 */
double Sweep::Processes::WeightedTransitionCoagulation::CoagKernel(const Particle &sp1,
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
 *\param[in]    sp1         First particle
 *\param[in]    sp2         Second particle
 *\param[in]    sys         Details of the environment, including temperature and pressure
 *\param[in]    maj         Flag to indicate which majorant kernel is required
 *
 *\return       Value of kernel
 */
double Sweep::Processes::WeightedTransitionCoagulation::MajorantKernel(const Particle &sp1,
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

/**
 * Calculate the free-molecular kernel between two particles in the given environment.
 * Either the majorant or non-majorant (true) kernel can be calculated.
 * The kernel is evaluated assuming that air is the surrounding gas (principally N2)
 *
 *\param[in]    sp1         First particle
 *\param[in]    sp2         Second particle
 *\param[in]    T        	Temperature
 *\param[in]    P        	Pressure
 *\param[in]    maj			Flag to indicate which majorant kernel is required
 *
 *\return       Value of kernel
 */
double Sweep::Processes::WeightedTransitionCoagulation::FreeMolKernel(const Particle &sp1, const Particle &sp2,
                                double T, double P, bool maj) const
{
    // This routine calculate the free molecular coagulation kernel for two particles.
    // There are two forms of kernel; a majorant form and a non-majorant form.

    // Collect the particle properties
    const double d1 = sp1.CollDiameter();
    const double d2 = sp2.CollDiameter();
    const double invm1 = 1.0 / sp1.Mass();
    const double invm2 = 1.0 / sp2.Mass();
    const double w2 = sp2.getStatisticalWeight();

    if (maj) {
        // The majorant form is always >= the non-majorant form.
        return CFMMAJ * m_efm * CFM * sqrt(T) * A() * w2 *
               (std::sqrt(invm1) + std::sqrt(invm2)) *
               (d1 * d1 + d2 * d2);
    } else {
        const double dterm = d1 + d2;
        return m_efm * CFM * A() * w2 *
               sqrt(T * (invm1 + invm2)) *
               dterm * dterm;
    }
}


/**
 * Calculate the slip-flow kernel between two particles in the given environment.
 * The kernel can be evaluated using air or Chapman-Enskog theory for viscosity
 *
 *\param[in]    sp1         First particle
 *\param[in]    sp2         Second particle
 *\param[in]    T        	Temperature
 *\param[in]    P        	Pressure
 *\param[in]    mu          Viscosity
 *\param[in]    maj			Flag to indicate which majorant kernel is required
 *
 *\return       Value of kernel
 */
double Sweep::Processes::WeightedTransitionCoagulation::SlipFlowKernel(const Particle &sp1, const Particle &sp2,
                                 double T, double P, double mu, bool maj) const
{
    // Collect the particle properties
    const double d1 = sp1.CollDiameter();
    const double d2 = sp2.CollDiameter();
    const double w2 = sp2.getStatisticalWeight();

    // For the slip-flow kernel the majorant and non-majorant forms are identical.
    return ((1.257 * 2.0 * MeanFreePathAir(T,P) *
             (1.0 / d1 / d1 + 1.0 / d2 / d2)) +
            (1.0 / d1 + 1.0 / d2)) *
           CSF * T * (d1 + d2) * w2
           * A() / mu;
}


/**
 * Write output to binary stream.
 * 
 * \param[in,out]   out     Binary output stream
 *
 * \exception   runtime_error   Output stream not ready
 */
void Sweep::Processes::WeightedTransitionCoagulation::Serialize(std::ostream &out) const
{
    // Serialise the parent class
    Coagulation::Serialize(out);

    // Note that m_efm is not written out in serialize as it is always directly
    // taken from the ParitcleModel. (allows const initialisation)

    // Now the data in this class
    if(out.good())
        out.write(reinterpret_cast<const char*>(&m_CoagWeightRule), sizeof(m_CoagWeightRule));
    else
        throw std::runtime_error("Output stream not ready (WeightedTransitionCoagulation::Serialize).");
}
