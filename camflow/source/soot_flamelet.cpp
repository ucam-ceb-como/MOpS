/*!
 * \file   soot_flamelet.cpp
 * \author Robert I A Patterson
 *
 * \brief  Interface for flamelet calculations including soot.
 *
 *  Copyright (C) 2009 Robert I A Patterson.
 *

 Licence:
    This file is part of "camflow".

    brush is free software; you can redistribute it and/or
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
    Prof Markus Kraft
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
#include <vector>
#include <numeric>

#include "soot_flamelet.h"

#include "pred_corr_solver.h"

/*!
 *
 *
 */
SootFlamelet::SootFlamelet()
    : mChemistry()
    , mParticles(NULL)
    , mCalcStarted(false)
{
    // Build a geometry object that matches the camflow geometry
    Geometry::Geometry1d geom = buildBrushGeometry();

    //========= Load chemical mechanism ==========================
    // This will load a separate copy of the mechanism to camflow,
    // but from the same files.
    Mops::Mechanism mech;

    std::cout << "Reading chemical species...\n";
    Sprog::IO::MechanismParser::ReadChemkin("chem.inp", mech, "therm.dat", "tran.dat");

    //========= Load particle mechanism ==========================
    std::cout << "Setting species on particle mechanism...\n";
    mech.ParticleMech().SetSpecies(mech.Species());
    std::cout << "Reading particle mechanism...\n";

    Sweep::MechParser::Read("sweep.xml", mech.ParticleMech());
    std::cout << "Read particle mechanism with " << mech.ParticleMech().ProcessCount()
              << " processes\n";

    // Max particle count and max m0 estimates
    std::vector<std::pair<real, real> > maxPCounts, maxM0s;
    maxPCounts.push_back(std::make_pair<real, real>(0, 1024));
    maxPCounts.push_back(std::make_pair<real, real>(1, 1024));
    maxM0s.push_back(std::make_pair<real, real>(0, 1e11));
    maxM0s.push_back(std::make_pair<real, real>(1, 1e11));

    // Finally create the object that will contain the particle populations.
    mParticles = new Brush::Reactor1d(geom, mech, maxPCounts, maxM0s);
}

/*!
 * The particles are stored in a dynamically allocoated Reactor1d which
 * must be deleted to avoid leaks.
 */
SootFlamelet::~SootFlamelet() {
    delete mParticles;
    mParticles = NULL;
}

/*!
 *@param[in]    Instance to copy
 */
/*SootFlamelet::SootFlamelet(const SootFlamelet& rhs)
    : mChemistry(rhs.mChemistry)
    , mParticles(NULL)
    , mCalcStarted(rhs.mCalcStarted)
{
    mParticles = new Brush::Reactor1d(*rhs.mParticles);
}*/

/*!
 *@param[in]    Value to assign
 *@return       The instance with its new value
 */
/*SootFlamelet& SootFlamelet::operator=(const SootFlamelet& rhs) {
    // avoid self-assignment corrupting the object
    if(this != &rhs) {
        mChemistry = rhs.mChemistry;
        mCalcStarted = rhs.mCalcStarted;

        // delete any particles and replace with a copy of those in rhs
        delete mParticles;
        mParticles = new Brush::Reactor1d(*rhs.mParticles);
    }
    return *this;
}*/



/*!
 *@return    Geometry1d with cell boundaries taken from Camflow
 */
Geometry::Geometry1d SootFlamelet::buildBrushGeometry() {
    std::vector<real> cellLengths = mChemistry.getCamGeometry().getGeometry();

    // There will be one more vertex than there are cells
    std::vector<real> vertices(cellLengths.size() + 1u);

    // First vertex is assumed to be at 0
    vertices.front() = 0.0;

    // Other vertices are found by adding up the cell lengths (note that the first
    // vertex is set to 0 above and skipped here.
    std::partial_sum(cellLengths.begin(), cellLengths.end(), vertices.begin() + 1);

    // @todo Remove hardcoding of boundary condition types.
    return Geometry::Geometry1d(vertices, Geometry::dirichlet, Geometry::dirichlet);
}

/*!
 *@param[in]    t_stop                  Time upto which to advance the flamelet
 *@param[in]    data_times              Vector of times at which data points apply
 *@param[in]    mix_frac_diff           Mixture fraction diffusion coefficient
 *@param[in]    grad_mix_frac           Gradient of mixture fraction
 *@param[in]    lapl_mix_frac           Laplacian of mixture fraction
 *@param[in]    grad_rho_mix_frac_diff  Gradient of gas mass density times mixture fraction diffusion coefficient
 *
 */
void SootFlamelet::run(const real t_stop, const vector<real>& data_times,
                       const std::vector<real>& mix_frac_diff,
                       const std::vector<real>& grad_mix_frac,
                       const std::vector<real>& lapl_mix_frac,
                       const std::vector<real>& grad_rho_mix_frac_diff) {
    //run the chemistry
    mChemistry.flamelet(sdr_data, data_times, mCalcStarted);
    mCalcStarted = true;

    // Extract the chemistry from mChemistry ready to use for the particle calculations
    Brush::ResetChemistry newChem = buildResetChemistry(data_times,
                                                        mix_frac_diff,
                                                        grad_mix_frac,
                                                        lapl_mix_frac,
                                                        grad_rho_mix_frac_diff);

    // Create and run the particle solver
    Brush::PredCorrSolver particleSolver(newChem, 0, 0.0, 0.0, true, true);
    particleSolver.solve(*mParticles, t_stop, 1, 0);
}

/*!
 *@param[in]    data_times              Vector of times at which data points apply
 *@param[in]    mix_frac_diff           Mixture fraction diffusion coefficient
 *@param[in]    grad_mix_frac           Gradient of mixture fraction
 *@param[in]    lapl_mix_frac           Laplacian of mixture fraction
 *@param[in]    grad_T                  Gradient of temperature
 *
 *@return    Object that can be used to set the gas phase mixture details on the 1d reactor
 *
 *@exception    std::invalid_argument   Input vectors must be of length 1 or 2
 *@excpetion    std::invalid_argument   Input vectors must all have same length
 */
Brush::ResetChemistry SootFlamelet::buildResetChemistry(const vector<real>& data_times,
                                                        const std::vector<real>& mix_frac_diff,
                                                        const std::vector<real>& grad_mix_frac,
                                                        const std::vector<real>& lapl_mix_frac,
                                                        const std::vector<real>& grad_T) {
    //========= Check the lengths of the input vectors ===============
    const size_t len = data_times.size();
    if((len != 1) && (len != 2)) {
        std::ostringstream msg("List of times to which data applies must have length 2, not ");
        msg << len << " (SootFlamelet::buildResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != mix_frac_diff.size()) {
        std::ostringstream msg;
        msg << "Length of mixture fraction diffusion coefficient vector is "
            << mix_frac_diff.size()
            << ", but it must match the time vector length of " << len
            << " (SootFlamelet::buildResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != grad_mix_frac.size()) {
        std::ostringstream msg;
        msg << "Length of mixture fraction gradient vector is "
            << grad_mix_frac.size()
            << ", but it must match the time vector length of " << len
            << " (SootFlamelet::buildResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != lapl_mix_frac.size()) {
        std::ostringstream msg;
        msg << "Length of mixture fraction Laplacian vector is "
            << lapl_mix_frac.size()
            << ", but it must match the time vector length of " << len
            << " (SootFlamelet::buildResetChemistry)";
        throw std::invalid_argument(msg.str());
    }

    if(len != grad_T.size()) {
        std::ostringstream msg;
        msg << "Length of gradient of temperature vector is "
            << grad_rho_mix_frac_diff.size()
            << ", but it must match the time vector length of " << len
            << " (SootFlamelet::buildResetChemistry)";
        throw std::invalid_argument(msg.str());
    }
    //========================================================================
    // Input vectors now guaranteed to all have length 1 or all to have
    // length 2.

    real mixFracDiffusion = mix_frac_diff.front();
    real gradMixFrac = grad_mix_frac.front();
    real laplMixFrac = lapl_mix_frac.front();
    real gradT = grad_T.front();

    // If there is only one value for each quantity use that value unchanged,
    // otherwise use a midpoint approach
    if(len == 2) {
        mixFracDiffusion += mix_frac_diff[1];
        mixFracDiffusion /= 2.0;

        gradMixFrac += grad_mix_frac[1];
        gradMixFrac /= 2.0;

        laplMixFrac += lapl_mix_frac[1];
        laplMixFrac /= 2.0;

        gradTn += grad_T[1];
        gradT /= 2.0;
    }

    const size_t numDataPoints = mChemistry.getIndepVars().size();

    // Assume that these values hold for all values of the indep variable,
    // so create vectors with an appropriate number of repeats
    std::vector<real> mixFracDiffusions(numDataPoints, mixFracDiffusion);
    std::vector<real> gradMixFracs(numDataPoints, gradMixFrac);
    std::vector<real> laplMixFracs(numDataPoints, laplMixFrac);
    std::vector<real> gradTs(numDataPoints, gradT);

    // Now collect the species mass fraction data
    std::vector<std::vector<real> > chemData(numDataPoints);
    for(size_t i = 0; i != numDataPoints; ++i) {
        chemData[i] = mChemistry.getMassFracsByPoint(i);
    }

    //PAH formation rate assumed 0 for now
    std::vector<real> pahFormations(numDataPoints, 0.0);

    return Brush::ResetChemistry(mChemistry.getIndepVars(),
                                 mChemistry.getTemperatures(),
                                 mChemistry.getDensities(), 
                                 mChemistry.getVelocities(),
                                 pahFormations,
                                 mixFracDiffuions,
                                 gradMixFracs,
                                 laplMixFracs,
                                 gradT,
                                 chemData);
}

