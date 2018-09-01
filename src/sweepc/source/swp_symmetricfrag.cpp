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

#include "swp_symmetricfrag.h"

#include "swp_params.h"
#include "swp_cell.h"
#include "swp_mechanism.h"

using namespace Sweep::Processes;

const double Sweep::Processes::SymmetricFragmentation::s_MajorantFactor = 2.0;

/**
 * Main way of building a new coagulation object
 * @param[in] mech      Mechanism to which coagulation will belong
 *
 */
Sweep::Processes::SymmetricFragmentation::SymmetricFragmentation(const Sweep::Mechanism &mech)
: Fragmentation(mech)
{
    m_name = "SymmetricFragmentation";
}

// Stream-reading constructor.
Sweep::Processes::SymmetricFragmentation::SymmetricFragmentation(std::istream &in, const Sweep::Mechanism &mech)
: Fragmentation(mech)
{
    m_name = "SymmetricFragmentation";
    Deserialize(in, mech);
}

// Returns the rate of the process for the given system.
double Sweep::Processes::SymmetricFragmentation::Rate(double t, const Cell &sys,
                                                        const Geometry::LocalGeometry1d &local_geom) const
{
    // Get system properties required to calculate coagulation rate.
    double rate = sys.Particles().GetSum(m_pid);

    rate *= A();

    if (m_mech->AnyDeferred()) {
        return rate * s_MajorantFactor;
    } else {
        return rate;
    }
}

/**
 * Number of terms in the expression for the sum of the majorant
 * kernel over all particle pairs.
 */
unsigned int Sweep::Processes::SymmetricFragmentation::TermCount() const {return TYPE_COUNT;}


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
double Sweep::Processes::SymmetricFragmentation::RateTerms(double t, const Cell &sys,
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
int SymmetricFragmentation::Perform(double t, Sweep::Cell &sys, 
                             const Geometry::LocalGeometry1d& local_geom,
                             unsigned int iterm,
                             Sweep::rng_type &rng) const
{
	PartPtrVector dummy;

    int i = -1;
    i = sys.Particles().Select(m_pid, rng);
    if (i >= 0) {
        Particle *sp = sys.Particles().At(i);
        std::vector<double> m_dcomp(1);
        m_dcomp[0] = - sp->NumCarbon() / 2;
        std::vector<double> m_dvals(1);
        m_dvals[0] = 0.0;
        std::vector<double> newComposition(1);
        newComposition[0] = - m_dcomp[0];
        if (m_mech->AnyDeferred()) {
            m_mech->UpdateParticle(*sp, sys, t, i, rng, dummy);
            if (sp->IsValid()) {
                double majk = MajorantKernel(*sp, sys, Default);
                double truek = FragKernel(*sp, sys);
                if (!Fictitious(majk, truek, rng)) {
                    sp->Adjust(m_dcomp, m_dvals, rng, 1);
                    sys.Particles().Update(i);
                    Particle * const sp = m_mech->CreateParticle(t);
                    sp->setPositionAndTime(-1.0, t);
                    sp->Primary()->SetComposition(newComposition);
                    sp->UpdateCache();
                    sys.Particles().Add(*sp, rng);
                }
            } else {
                sys.Particles().Remove(i);
            }
        } else {
            sp->Adjust(m_dcomp, m_dvals, rng, 1);     
            sys.Particles().Update(i);
            Particle * const sp = m_mech->CreateParticle(t);
            sp->setPositionAndTime(-1.0, t);
            sp->Primary()->SetComposition(newComposition);
            sp->UpdateCache();
            sys.Particles().Add(*sp, rng);
        }
    } else {
        return -1;
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
double Sweep::Processes::SymmetricFragmentation::FragKernel(const Particle &sp, const Cell& sys) const
{
    return sp.Property(m_pid) * A();
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
double Sweep::Processes::SymmetricFragmentation::MajorantKernel(const Particle &sp,
                                                                const Cell& sys,
                                                                const MajorantType maj) const
{
    return FragKernel(sp, sys) * s_MajorantFactor;
}

