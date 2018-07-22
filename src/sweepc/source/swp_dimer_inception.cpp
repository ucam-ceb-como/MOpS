/*
  Author(s):      Robert Patterson and Markus Sander
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the DimerInception class declared in the
    swp_DimerInception.h header file.

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

#include "swp_dimer_inception.h"

#include "swp_mechanism.h"
#include "swp_primary.h"

#include <boost/random/uniform_01.hpp>
#include <boost/random/lognormal_distribution.hpp>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
DimerInception::DimerInception(void)
: Inception(), m_kfm(0.0), m_ksf1(0.0), m_ksf2(0.0), m_efm(2.2)
{
    m_name = "DimerInception";
}

// Initialising constructor.
DimerInception::DimerInception(const Sweep::Mechanism &mech)
: Inception(mech), m_kfm(0.0), m_ksf1(0.0), m_ksf2(0.0),
  m_efm(mech.GetEnhancementFM())
{
    m_name = "DimerInception";
}

// Copy constructor.
DimerInception::DimerInception(const DimerInception &copy)
: m_efm(copy.m_efm)
{
    *this = copy;
}

// Stream-reading constructor.
DimerInception::DimerInception(std::istream &in, const Sweep::Mechanism &mech)
: m_efm(mech.GetEnhancementFM())
{
    Deserialize(in, mech);
}

// Default destructor.
DimerInception::~DimerInception(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
DimerInception &DimerInception::operator =(const DimerInception &rhs)
{
    if (this != &rhs) {
        Inception::operator =(rhs);
        m_kfm  = rhs.m_kfm;
        m_ksf1 = rhs.m_ksf1;
        m_ksf2 = rhs.m_ksf2;
    }
    return *this;
}

/*!
 * Create a new particle and add it to the ensemble with position uniformly
 * distributed over the grid cell.
 *
 * The iterm parameter is included because it will be needed for many process
 * types and this function is meant to have a general signature.
 *
 * \param[in]       t               Time
 * \param[in]       local_geom      Details of geometry around current location
 * \param[in,out]   sys             System to update
 * \param[in]       iterm           Process term responsible for this event
 * \param[in,out]   rng             Random number generator
 *
 * \return      0 on success, otherwise negative.
 */
int DimerInception::Perform(const double t, Cell &sys,
                            const Geometry::LocalGeometry1d &local_geom,
                            const unsigned int iterm,
                            rng_type &rng) const {

    // This routine performs the inception on the given chemical system.

	// aab64 Inception approximation ideas to reduce ensemble filling rate
	// (1) PSI - particle surface inception: do surface growth on existing large particle
	// (2) HCI - heavy cluster inception: incept larger primary particle
	// Caution - may introduce large errors in the primary size distribution.

	// aab64 hybrid particle model
	// If hybrid_flag is active, track the number of incepting particles
	int iprng = -1;
	int iprng2 = -1;
	unsigned int nsp = sys.ParticleCount();
	double dcol_1, dcol_2, sw_0, sw_1;
	Particle *sprng = NULL;
	Particle *sprng2 = NULL;
	// Get surface inception settings
	bool surfincflag = m_mech->GetIsSurfInc();
	double dcol_switch = m_mech->GetSurfIncValue();
	double dcol_switch_min = m_mech->GetSurfIncCutoffValue();
	bool sizeflag = false;
	std::string PSItype;
	m_mech->GetPSItype(PSItype); // event, both, weight
	Sweep::PropID prop1 = sys.getCoagProp1();
	Sweep::PropID prop2 = sys.getCoagProp2();
	if (nsp > 1 && surfincflag) 
	{
		// Get average particle collision diameter
		//dcol_ave = sys.Particles().GetSum(Sweep::iDW) / sys.Particles().GetSum(Sweep::iW);

		// Select a particle at random, weighted by collision diameter sqrd
		//Sweep::PropID prop1 = iD2;
		iprng = sys.Particles().Select(prop1, rng);
		if (iprng >= 0)
			sprng = sys.Particles().At(iprng);
		else
			return -1; // Failed to choose a particle.

		// Select a second particle at random, weighted by inverse sqrt mass times weight
		//Sweep::PropID prop2 = iM_1_2W;
		iprng2 = sys.Particles().Select(prop2, rng);
		if (iprng2 >= 0) 
			sprng2 = sys.Particles().At(iprng2);
		else 
			return -1; // Failed to choose a particle.

		// Toggle size flag if selected particle has collision diameter > switch collision diameter
		dcol_1 = sprng->CollDiameter();
		dcol_2 = sprng2->CollDiameter();
		if (dcol_2 > dcol_1)
			sprng = sprng2;

		sw_0 = sprng->getStatisticalWeight();
		sizeflag = (sys.ParticleCount() == sys.Particles().Capacity()); 
		// aab64 temp change (dcol_1 > dcol_switch) && (dcol_2 <= dcol_switch_min);
	}

		// If surface inception is active AND average particle size is large enough
		// THEN do surface inception on a random particle.
		if (surfincflag && sizeflag)
		{
			// Inception stoichiometry determines particle composition
			unsigned int nInceptingParticle_ui = (unsigned int)ParticleComp()[0];
			double nInceptingParticle_d = ParticleComp()[0];

			// PSItype = event: Update the number of rutiles in the particle, leave the weight
			if (PSItype == "E")
			{
				// Use inception composition and weight ratio to determine number of units to add to particle
				unsigned int weightRatio = (unsigned int)(sys.GetInceptingWeight() / sprng->getStatisticalWeight());
				unsigned int nChosenParticle = weightRatio * nInceptingParticle_ui;

				// Avoid doing extra update of gas-phase and temperature during surface growth
				// Store a flag in the sys object that can be checked inside Perform(.)
				sys.SetNotPSIFlag(false);

				// 2. Call perform for surface growth process, and do one surface event.
				// ParticleComp()[0] is the number of TiO2 units added, dx, in this case...
				// How does this generalise for other systems with different # processes and comp?
				int m = m_mech->Processes(0)->Perform(t, sys, *sprng, rng, nChosenParticle);

				// Reset flag in the sys object for outside this function
				sys.SetNotPSIFlag(true);

				// 3. Update the cached properties (note could this be left off to make 
				// it more efficient, assuming that the small change is insignificant?).
				sprng->UpdateCache();
			}
			// PSItype = weight: Update the particle weight and leave the number of rutiles unchanged
			else if (PSItype == "W")
			{
				double nChosenParticle = sprng->Composition()[0];
				double newWeight = (nChosenParticle * sprng->getStatisticalWeight()) + (nInceptingParticle_d * sys.GetInceptingWeight());
				newWeight *= (1.0 / (nChosenParticle));
				sprng->setStatisticalWeight(newWeight);
			}
			// PSItype = both: Update both the weight and the number of rutiles of the particle
			else
			{
				// Avoid doing extra update of gas-phase and temperature during surface growth
				// Store a flag in the sys object that can be checked inside Perform(.)
				sys.SetNotPSIFlag(false);

				// 2. Call perform for surface growth process, and do one surface event.
				// ParticleComp()[0] is the number of TiO2 units added, dx, in this case...
				// How does this generalise for other systems with different # processes and comp?
				int m = m_mech->Processes(0)->Perform(t, sys, *sprng, rng, nInceptingParticle_d);

				// Reset flag in the sys object for outside this function
				sys.SetNotPSIFlag(true);

				// Update the weight 
				double nChosenParticle = sprng->Composition()[0];
				double newWeight = (nChosenParticle * sprng->getStatisticalWeight());
				newWeight += (nInceptingParticle_d * sys.GetInceptingWeight());
				newWeight *= (1.0 / (nChosenParticle + nInceptingParticle_d));
				sprng->setStatisticalWeight(newWeight);

				// 3. Update the cached properties (note could this be left off to make 
				// it more efficient, assuming that the small change is insignificant?).
				sprng->UpdateCache();
			}

			// Update gas-phase chemistry of system 
			double particleWt = sys.GetInceptingWeight();
			adjustGas(sys, particleWt, 1, sys.GetInceptionFactor()); // nInceptingParticle_ui
			adjustParticleTemperature(sys, particleWt, 1, sys.GetIsAdiabaticFlag(), nInceptingParticle_d, 1, sys.GetInceptionFactor());
		}
		else
		{

			// Get the cell vertices
			fvector vertices = local_geom.cellVertices();

			// Sample a uniformly distributed position, note that this method
			// works whether the vertices come in increasing or decreasing order,
			// but 1d is assumed for now.
			double posn = vertices.front();

			const double width = vertices.back() - posn;
			boost::uniform_01<rng_type&, double> uniformGenerator(rng);
			posn += width * uniformGenerator();

			// Add particle to system's ensemble.
			if (!m_mech->IsHybrid())
			{
				// Create a new particle of the type specified
				// by the system ensemble.
				Particle *sp = m_mech->CreateParticle(t);

				// aab64 Get incepting particle weight for cases where wt != 1.0:
				// Check if SWA is in play and variable inception weighting is active.
				// Update newly incepted particle weight if necessary.
				if (m_mech->IsWeightedCoag() && m_mech->IsVariableWeightedInception()) {
					sp->setStatisticalWeight(sys.GetInceptingWeight());
				}

				sp->setPositionAndTime(posn, t);

				// Initialise the new particle.
				sp->Primary()->SetComposition(ParticleComp());
				sp->Primary()->SetValues(ParticleTrackers());

				// aab64 Adjust composition to allow inception of heavier particles
				// Perform adjustment before updating the cache and computing properties.
				bool heavyAllowed = m_mech->GetIsHeavy();
				if (heavyAllowed)
					sp->Primary()->AdjustForInception(sys.GetInceptionFactor());

				sp->UpdateCache();

				sys.Particles().Add(*sp, rng);

				// Update gas-phase chemistry of system.
				adjustGas(sys, sp->getStatisticalWeight(), 1, sys.GetInceptionFactor());
				adjustParticleTemperature(sys, sp->getStatisticalWeight(), 1, sys.GetIsAdiabaticFlag(), ParticleComp()[0], 1, sys.GetInceptionFactor());
			}
			else
			{
				if (!sys.Particles().IsFirstSP())
				{
					// Flag that register of particle properties is set up
					sys.Particles().SetInceptedSP();

					// Initialise lookup of particles below critical size
					for (unsigned int i = 0; i < sys.Particles().GetCritialNumber(); i++)
					{
						Particle * sp_pn = m_mech->CreateParticle(t);
						std::vector<double> newComposition(1);
						newComposition[0] = i;
						sp_pn->setPositionAndTime(posn, t);
						sp_pn->Primary()->SetComposition(newComposition);
						sp_pn->Primary()->SetValues(ParticleTrackers());
						sp_pn->UpdateCache();
						sys.Particles().SetPNParticle(*sp_pn, rng, i);
					}
					sys.Particles().InitialiseDiameters(sys.ParticleModel()->Components()[0]->MolWt(),
						sys.ParticleModel()->Components()[0]->Density()); // Works for current TiO2 -> Need to generalise
				}

				// Adjust particle number properties
				sys.Particles().UpdateNumberAtIndex(ParticleComp()[0], 1);
				sys.Particles().UpdateTotalParticleNumber(1);
				sys.Particles().UpdateTotalsWithIndex(ParticleComp()[0], 1.0);

				// Update gas-phase chemistry of system.
				adjustGas(sys, sys.GetInceptingWeight(), 1, sys.GetInceptionFactor());
				adjustParticleTemperature(sys, sys.GetInceptingWeight(), 1, sys.GetIsAdiabaticFlag(), ParticleComp()[0], 1, sys.GetInceptionFactor());
			}
		}

    return 0;
}

// PERFORMING THE PROCESS.


// INCEPTION KERNEL.

/*!
 * The inception rate is based on a transition coagulation kernel that
 * is calculated using a molecule diameter and mass.  These values do not
 * have to correspond the the physical properties of the molecule, but they
 * are fed into the transition regime kernel.
 *
 * @param[in]    m1    mass of first molecule
 * @param[in]    m2    mass of second molecule
 * @param[in]    d1    diameter of first molecule
 * @param[in]    d2    diameter of second molecule
 */
void DimerInception::SetInceptingSpecies(double m1, double m2, double d1, double d2)
{
    // The free mol part can be handled by the free mol specific method
    SetInceptingSpeciesFreeMol(m1, m2, d1, d2);

    // Now the slip flow part
    double invd1=1.0/d1, invd2=1.0/d2;
    m_ksf1 = CSF * (d1+d2);
    m_ksf2 = 2.0 * 1.257 * m_ksf1 * ((invd1*invd1) + (invd2*invd2));
    m_ksf1 = m_ksf1 * (invd1+invd2);
}

/*!
 * The inception rate is based on a free molecular coagulation kernel only,
 * contrast \ref SetInceptingSpecies that
 * is calculated using a molecule diameter and mass.  These values do not
 * have to correspond the the physical properties of the molecule, but they
 * are fed into the transition regime kernel.
 *
 * @param[in]    m1    mass of first molecule
 * @param[in]    m2    mass of second molecule
 * @param[in]    d1    diameter of first molecule
 * @param[in]    d2    diameter of second molecule
 */
void DimerInception::SetInceptingSpeciesFreeMol(double m1, double m2, double d1, double d2)
{
    // This routine sets the free-mol and slip flow kernel parameters given
    // the mass and diameter of the incepting species.
    m_kfm  = m_efm * CFM * sqrt((1.0/m1) + (1.0/m2)) * (d1+d2) * (d1+d2);
    m_ksf1 = 0.0;
    m_ksf2 = 0.0;
}

// TOTAL RATE CALCULATIONS.

// Returns rate of the process for the given system.
double DimerInception::Rate(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom) const
{
    // Get the current chemical conditions.
    double T = sys.GasPhase().Temperature();
    double P = sys.GasPhase().Pressure();

	// aab64 Divide the inception rate by the current particle weight
	// aab64 Divide by current incepting composition scale factor for cases with heavier inception
	// This really belongs in the rate fn below for consistency. 
	// To do: ensure that it is always correctly implemented when rates 
	// are computed. 
	double scaleFac = sys.GetInceptionFactor() * sys.GetInceptingWeight();
	double scaledRate = Rate(sys.GasPhase(), sqrt(T),
		MeanFreePathAir(T, P),
		sys.SampleVolume());
	scaledRate *= (1.0 / scaleFac); 

    // Calculate the rate.
    return scaledRate;
}



/*!
 * Calculate inception rate using a the transition coagulation kernel
 * with the values provided by the user.  The result is this value
 * multiplied by the square of the number concentration of the gas
 * phase species.
 *
 * Requires all the parameters that would otherwise be calculated by
 * the routine to be passed as arguments.
 *
 * @param[in]    gas      Gas phase mixture
 * @param[in]    sqrtT    square root of temperature
 * @param[in]    MFP      mean free path in gas
 * @param[in]    vol      sample volume
 *
 * @return    Inception rate for a cell of size vol. (\f$ \mathrm{s}^{-1}\f$)
 */
double DimerInception::Rate(const EnvironmentInterface &gas, double sqrtT,
                     double MFP, double vol) const
{
    double rate = A() * vol * chemRatePart(gas);

    const double fm   = sqrtT * m_kfm;
    if((m_ksf1 > 0) || (m_ksf2 > 0))  {
        const double Temperature = gas.Temperature();

        // Temperature divided by viscosity
        double T_viscosity = Temperature / gas.Viscosity();

        // Transition regime
        double sf   = T_viscosity  * (m_ksf1 + (MFP*m_ksf2));
        rate *= ((fm*sf) / (fm+sf));
    }
    else {
        // Free mol regime only
        rate *= fm;
    }

    return rate;
}

/*!
 * Calculates the gas-phase chemistry contribution to the rate
 * expression.  This is overloaded as Avogadro's number must be
 * included in the terms for inception processes.
 *
 * @param[in]    gas      Gas phase mixture
 */
double DimerInception::chemRatePart(const EnvironmentInterface &gas) const
{
    // Factor of 0.5 adjusts for doubling counting of pairs of molecules in the number
    // of possible collisions.
    double rate = 0.5;

    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        //std::cerr << "Mole frac to use " << fracs[i->first] << std::endl;
        double conc = gas.SpeciesConcentration(i->first);
        for (int j=0; j!=i->second; ++j) {
            rate *= (NA * conc);
        }
    }
	
    return rate;
}


// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int DimerInception::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a double vector. The
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
double DimerInception::RateTerms(const double t, const Cell &sys,
                               const Geometry::LocalGeometry1d &local_geom,
                               fvector::iterator &iterm) const
{
    // Get the current chemical conditions.
    double T = sys.GasPhase().Temperature();
    double P = sys.GasPhase().Pressure();

    // Calculate the single rate term and advance iterator.
    
    if (sys.ParticleModel()->Postprocessing() == ParticleModel::wdotA4) {
        double Rate = NA * sys.GasPhase().PropertyValue(1007) * sys.SampleVolume();

        if (Rate < 0.0)
            Rate = 0.0;

        *iterm = Rate; 
    } else {
		// aab64 Divide by current incepting particle weight for cases when wt != 1.0. 
		// aab64 Divide by current incepting composition scale factor for cases with heavier inception
		double scaleFac = 1.0 / (sys.GetInceptingWeight() * sys.GetInceptionFactor());
        *iterm = Rate(sys.GasPhase(), sqrt(T),
                      MeanFreePathAir(T,P),
					  sys.SampleVolume()) * scaleFac; 
    }

    return *(iterm++);
}


// READ/WRITE/COPY.

// Creates a copy of the inception.
DimerInception *const DimerInception::Clone(void) const {return new DimerInception(*this);}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType DimerInception::ID(void) const {return Dimer_Inception_ID;}

// Writes the object to a binary stream.
void DimerInception::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        Inception::Serialize(out);

        // Write free-mol parameter.
        double v = (double)m_kfm;
        out.write((char*)&v, sizeof(v));

        // Write slip-flow parameters.
        v = (double)m_ksf1;
        out.write((char*)&v, sizeof(v));
        v = (double)m_ksf2;
        out.write((char*)&v, sizeof(v));

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, DimerInception::Serialize).");
    }
}

// Reads the object from a binary stream.
void DimerInception::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;

        switch (version) {
            case 0:
                // Deserialize base class.
                Inception::Deserialize(in, mech);

                // Read free-mol parameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_kfm = (double)val;

                // Read slip-flow parameters.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_ksf1 = (double)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_ksf2 = (double)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, DimerInception::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, DimerInception::Deserialize).");
    }
}
