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
#include "swp_particle_cache.h"
#include "settings_io.h"

#include "swp_ensemble_stats.h"


/*!
 * Load all settings from files with standard names.
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


    //========= Load numerical parameters for brush ==============
    CamXML::Document settingsXML;
    settingsXML.Load("brush.xml");
    const CamXML::Element * const root = settingsXML.Root();

    // Maximum number of computational particles as a function
    // of position.
    std::vector<std::pair<real, real> > maxPCounts =
        Brush::Settings_IO::readProfile(root, "pcount");


    // Note that the max m0 estimates are interpreted as having
    // units of m^-3, in contrast to Mops which uses cm^-3
    // for inputs.  Max m0 estimates are used to estimate an
    // initial sample volume (scaling factor)
    std::vector<std::pair<real, real> > maxM0s =
        Brush::Settings_IO::readProfile(root, "maxm0");

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

// These two methods cannot be implemented until the corresponding methods are
// implemented for Interface.
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
 * Run the flamelet upto the time specified by the last element of data_times.
 *
 *@param[in]    data_times              Vector of times at which data points apply
 *@param[in]    mix_frac_diff           Mixture fraction diffusion coefficient
 *@param[in]    grad_mix_frac           Gradient of mixture fraction
 *@param[in]    lapl_mix_frac           Laplacian of mixture fraction
 *@param[in]    grad_rho_mix_frac_diff  Gradient of gas mass density times mixture fraction diffusion coefficient
 *
 *@exception    std::invalid_argument   Input vectors must be of length 1 or 2
 *@excpetion    std::invalid_argument   Input vectors must all have same length
 *
 */
void SootFlamelet::run(const std::vector<real>& data_times,
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
            << grad_T.size()
            << ", but it must match the time vector length of " << len
            << " (SootFlamelet::buildResetChemistry)";
        throw std::invalid_argument(msg.str());
    }
    //========================================================================
    // Input vectors now guaranteed to all have length 1 or all to have
    // length 2.

    // make a vector of scalar dissipation rates
    std::vector<real> scalarDissipRate(grad_mix_frac.size());
    for(size_t i = 0; i != grad_mix_frac.size(); ++i) {
        scalarDissipRate[i] = grad_mix_frac[i] * grad_mix_frac[i] * mix_frac_diff[i] * 2;
    }

    //=========================================================================
    // get a vector of soot volume fractions
    std::vector<real> sootFv(mParticles->getNumCells());
    for(size_t i = 0; i != mParticles->getNumCells(); ++i) {
        sootFv[i] = mParticles->getCell(i).Mixture()->Particles().GetSum(Sweep::ParticleCache::iV) /
                    mParticles->getCell(i).Mixture()->SampleVolume();
    }

    // run the chemistry
    mChemistry.flameletWithSoot(sootFv, scalarDissipRate, data_times, mCalcStarted);
    mCalcStarted = true;

    //=========================================================================
    // Particles
    // Extract the chemistry from mChemistry ready to use for the particle calculations
    Brush::ResetChemistry newChem = buildResetChemistry(data_times,
                                                        mix_frac_diff,
                                                        grad_mix_frac,
                                                        lapl_mix_frac,
                                                        grad_T);

    // Create and run the particle solver
    Brush::PredCorrSolver particleSolver(newChem, 0, 0.0, 0.0, true, true);
    particleSolver.solve(*mParticles, data_times.back(), 1, 0);
}

/*!
 * Load in an input file for the chemistry instead of calling camflow to solve
 * a flamelet. 
 */
void SootFlamelet::run(const real stop_time, const std::string fileName) {

    //=========================================================================
    // Particles
    // Extract the chemistry from mChemistry ready to use for the particle calculations

    // Declare input variables for Brush::ResetChemistry
    int verbosity = 1;
    Sprog::Mechanism mech;
    //std::string inputFileName = "camflowFlamelet.dat";   
    Brush::ResetChemistry::InputFileType file_type = Brush::ResetChemistry::CamflowFlamelet;

    // Set the stop time. If calling from OpenFOAM, it will be whatever has been 
    // set for 'stopTime' in constant/flowToFlameletDict.
    std::vector<real> data_times(1);
    data_times[0] = stop_time;

    // Read the chemical species into Sprog::Mechanism mech.
    std::cout << "Reading chemical species...\n";
    Sprog::IO::MechanismParser::ReadChemkin("chem.inp", mech, "therm.dat", "tran.dat");

    Brush::ResetChemistry newChem = Brush::ResetChemistry(fileName,
                                                          file_type,
                                                          mech,
                                                          verbosity);

    // Create and run the particle solver
    Brush::PredCorrSolver particleSolver(newChem, 0, 0.0, 0.0, true, true);

    mParticles->ReplaceChemistry(newChem,false);

    particleSolver.solve(*mParticles, data_times.back(), 1, 0);

}

/*!
 *@param[in]    data_times              Vector of times at which data points apply
 *@param[in]    mix_frac_diff           Mixture fraction diffusion coefficient
 *@param[in]    grad_mix_frac           Gradient of mixture fraction
 *@param[in]    lapl_mix_frac           Laplacian of mixture fraction
 *@param[in]    grad_T                  Gradient of temperature
 *
 *@pre  All of the agument vectors must have the same length and this length
 *      must be 1 or 2.
 *
 *@return    Object that can be used to set the gas phase mixture details on the 1d reactor
 */
Brush::ResetChemistry SootFlamelet::buildResetChemistry(const std::vector<real>& data_times,
                                                        const std::vector<real>& mix_frac_diff,
                                                        const std::vector<real>& grad_mix_frac,
                                                        const std::vector<real>& lapl_mix_frac,
                                                        const std::vector<real>& grad_T) {
    real mixFracDiffusion = mix_frac_diff.front();
    real gradMixFrac = grad_mix_frac.front();
    real laplMixFrac = lapl_mix_frac.front();
    real gradT = grad_T.front();

    // If there is only one value for each quantity use that value unchanged,
    // otherwise use a midpoint approach
    if(data_times.size() == 2) {
        mixFracDiffusion += mix_frac_diff[1];
        mixFracDiffusion /= 2.0;

        gradMixFrac += grad_mix_frac[1];
        gradMixFrac /= 2.0;

        laplMixFrac += lapl_mix_frac[1];
        laplMixFrac /= 2.0;

        gradT += grad_T[1];
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
                                 mixFracDiffusions,
                                 gradMixFracs,
                                 laplMixFracs,
                                 gradTs,
                                 chemData);
}

/*!
 * Get a vector of the moments and other statistics of the soot particle
 * population.  The result will have length 0 if statistics are requested for
 * a position outside the range of the simulation grid.
 *
 *@param[in]    position    Mixture fraction value at which statistics required
 *
 *@return       Vector of statistics
 */
std::vector<real> SootFlamelet::getSootStatistics(const real position) const {
    const int cellIndex = mParticles->getGeometry().containingCell(position);

    if(cellIndex < 0) {
        // Position is not inside any cell so return an empty vector
        return std::vector<real>(0u);
    }

    // Create a stats object
    Sweep::Stats::EnsembleStats stats(mParticles->getMechanism().ParticleMech());

    // Reject particles with a collision diameter outside the range (0, 1.0e30)
    Sweep::Stats::IModelStats::StatBound statsBound;
    statsBound.Lower = 0;
    statsBound.Upper = 1.0e30;
    statsBound.PID = Sweep::ParticleCache::iDcol;
    stats.SetStatBoundary(statsBound);

    // Collect the statistics
    stats.Calculate(mParticles->getCell(cellIndex).Mixture()->Particles(),
                    1.0 / mParticles->getCell(cellIndex).Mixture()->SampleVolume());

    return stats.Get();

}

/*!
 * Get a list of vectors each describing one computational particle from the
 * population.  The list will have length 0 if requested for a position outside
 * the range of the simulation grid.
 *
 *@param[in]    position    Mixture fraction value at which statistics required
 *
 *@return       List of particle description vectors
 */
std::list<std::vector<real> > SootFlamelet::getSootParticleList(const real position) const {
    // Find index of cell for which the particle list is requested
    const int cellIndex = mParticles->getGeometry().containingCell(position);

    // Collect result in this variable (default constructs with length 0)
    std::list<std::vector<real> > particleList;

    if(cellIndex < 0) {
        // Position is not inside any cell so return an empty list
        return particleList;
    }

    // Convenient way of referring to the mixture in the cell of interest
    const Mops::Mixture &mix = *(mParticles->getCell(cellIndex).Mixture());
    const Brush::real time = mParticles->getTime();

    // Create a stats object
    Sweep::Stats::EnsembleStats stats(mParticles->getMechanism().ParticleMech());

    // Loop over the particles
    for(unsigned int particleIndex = 0; particleIndex < mix.ParticleCount(); ++particleIndex) {
        // Get a vector containing the particle details
        fvector particleListEntry;
        stats.PSL(mix.Particles(), particleIndex, time,
                  particleListEntry, 1.0 / mix.SampleVolume());

        particleList.push_back(particleListEntry);
    }

    return particleList;
}
