/*!
 * @file    swp_silicon_inception.cpp
 * @author  William J Menz
 * @brief   Implementation of the SiliconInception class.
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      Based on swp_dimer_inception by Robert Patterson, Markus Sander
 *      and Matthew Celnik.
 *
 *      Implementation of the SiliconInception class.
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
 *      sweepc is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public License
 *      as published by the Free Software Foundation; either version 2
 *      of the License, or (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this program; if not, write to the Free Software
 *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *      02111-1307, USA.
 *
 *   Contact:
 *      Prof Markus Kraft
 *      Dept of Chemical Engineering
 *      University of Cambridge
 *      New Museums Site
 *      Pembroke Street
 *      Cambridge
 *      CB2 3RA, UK
 *
 *      Email:       mk306@cam.ac.uk
 *      Website:     http://como.cheng.cam.ac.uk
*/

#include "swp_silicon_inception.h"
#include "swp_mechanism.h"
#include "swp_params.h"
#include "swp_primary.h"
#include "swp_sprog_idealgas_wrapper.h"

#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// Mass of a silicon monomer, kg
const double SiliconInception::m_m1  = 4.664e-26;
// Vol of a silicon atom, m3
const double SiliconInception::m_v1  = 2.003e-29;
// Spherical diameter of a bulk silicon atom, m
const double SiliconInception::m_d1  = 3.369e-10;


// CONSTRUCTORS AND DESTRUCTORS.

//! Default constructor (protected).
SiliconInception::SiliconInception(void)
: Inception(), m_kfm(0.0), m_ksf1(0.0), m_ksf2(0.0),
  m_efm(2.2),
  m_itype(iCollisional), m_vi(0.0), m_di(0.0),
  m_sidata(0), m_reacs(0), m_concs(0)
{
    m_name = "SiliconInception";
}

//! Initialising constructor.
SiliconInception::SiliconInception(const Sweep::Mechanism &mech)
: Inception(mech), m_kfm(0.0), m_ksf1(0.0), m_ksf2(0.0),
  m_efm(mech.GetEnhancementFM()),
  m_itype(iCollisional), m_vi(0.0), m_di(0.0),
  m_sidata(0), m_reacs(0), m_concs(0)
{
    m_name = "SiliconInception";
}

//! Copy constructor.
SiliconInception::SiliconInception(const SiliconInception &copy)
: m_efm(copy.m_efm)
{
    *this = copy;
}

//! Stream-reading constructor.
SiliconInception::SiliconInception(std::istream &in, const Sweep::Mechanism &mech)
: m_efm(mech.GetEnhancementFM())
{
    Deserialize(in, mech);
}

//! Default destructor.
SiliconInception::~SiliconInception(void)
{
}

// OPERATOR OVERLOADS.

/*!
 * @brief           Assignment operator
 *
 * @param rhs       Pointer to right hand side object
 * @return          Pointer to object
 */
SiliconInception &SiliconInception::operator =(const SiliconInception &rhs)
{
    if (this != &rhs) {
        Inception::operator =(rhs);
        m_kfm  = rhs.m_kfm;
        m_ksf1 = rhs.m_ksf1;
        m_ksf2 = rhs.m_ksf2;
        m_itype = rhs.m_itype;
        m_vi   = rhs.m_vi;
        m_di   = rhs.m_di;

        m_sidata = rhs.m_sidata;
        m_reacs  = rhs.m_reacs;
        m_concs  = rhs.m_concs;
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
int SiliconInception::Perform(const double t, Cell &sys,
                            const Geometry::LocalGeometry1d &local_geom,
                            const unsigned int iterm,
                            rng_type &rng) const {

    // This routine performs the inception on the given chemical system.

    // Create a new particle of the type specified
    // by the system ensemble.
    Particle *sp = m_mech->CreateParticle(t);

    // Get the cell vertices
    fvector vertices = local_geom.cellVertices();

    // Sample a uniformly distributed position, note that this method
    // works whether the vertices come in increasing or decreasing order,
    // but 1d is assumed for now.
    double posn = vertices.front();

    const double width = vertices.back() - posn;
    boost::uniform_01<rng_type&, double> uniformGenerator(rng);
    posn += width * uniformGenerator();

    sp->setPositionAndTime(posn, t);

    if (m_itype == iCollisional) {

        // If the normal collisional-type rate is used, create the particle
        // and adjust the gas-phase as in normal Inception

        // Initialise the new particle.
        sp->Primary()->SetComposition(ParticleComp());
        sp->Primary()->SetValues(ParticleTrackers());
        sp->UpdateCache();

        // Add particle to system's ensemble.
        sys.Particles().Add(*sp, rng);

        // Update gas-phase chemistry of system.
        adjustGas(sys, sp->getStatisticalWeight());

    } else {

        // Otherwise, use information from m_sidata to select a gas-phase
        // species to transform into a new particle
        const SiliconData* species = ChooseData(sys, rng);
        if (species != NULL) {
            sp->Primary()->SetComposition(species->_track);
            sp->Primary()->SetValues(ParticleTrackers());
            sp->UpdateCache();

            // Add particle to system's ensemble.
            sys.Particles().Add(*sp, rng);

            // Update gas-phase chemistry of system.
            adjustGasPhase(sys, (*species), sp->getStatisticalWeight());
        } else {
            std::cout << "Warning: silicon inception chosen with no precursor."
                    << endl;
        }
    }

    return 0;
}


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
void SiliconInception::SetInceptingSpecies(double m1, double m2, double d1, double d2)
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
void SiliconInception::SetInceptingSpeciesFreeMol(double m1, double m2, double d1, double d2)
{
    // This routine sets the free-mol and slip flow kernel parameters given
    // the mass and diameter of the incepting species.
    m_kfm  = m_efm * CFM * sqrt((1.0/m1) + (1.0/m2)) * (d1+d2) * (d1+d2);
    m_ksf1 = 0.0;
    m_ksf2 = 0.0;
}

/*!
 * @brief           Sets the incepting volume
 *
 * The volume of the incepting particle is determined upon initialisation,
 * based on the density of molwt of each component.
 *
 * @param mech      Particle mechanism
 */
void SiliconInception::SetInceptingVolume(const Sweep::Mechanism &mech)
{
    double v1(0.0);
    CompPtrVector comp = mech.Components();
    for (unsigned int i=0; i!=comp.size(); i++) {
        // returns in m3, due to Density()
        v1 += 1.0 * mech.Components(i)->MolWt() /
                (Sweep::NA * mech.Components(i)->Density());
    }

    m_vi = v1;
}

/*!
 * @brief       Calculates the diameter of an incepting particle
 *
 * Pre-calculates the incepting diameter of a particle (when a collisional
 * inception process is used) for this inception reaction. This is done
 * using the molwt, density and initial number of each component in
 * the new particle.
 *
 * @param mech  Particle mechanism
 */
void SiliconInception::SetInceptingDiameter(const Sweep::Mechanism &mech)
{
    if (m_itype != iCollisional) {
        double dmin(1.0);
        for (unsigned int i=0; i!=m_sidata.size(); i++) {
            if (m_sidata[i]._diam < dmin) dmin = m_sidata[i]._diam;
        }
        m_di = dmin;
    } else {
        //Pre-define volume and diameters
        double v(0.0), d(0.0);

        // Loop over list of components in incepting particle to find initial vol
        for (unsigned int i=0; i!=ParticleComp().size(); i++) {

            // returns in m3, due to Density()
            v += ParticleComp()[i] * mech.Components(i)->MolWt() /
                    (Sweep::NA * mech.Components(i)->Density());
        }

        // Now calculate the equivalent spherical diameter
        d = pow(6 * v / Sweep::PI, Sweep::ONE_THIRD);
        m_di = d;
    }
}

/*!
 * @brief           Generates the silicon species data
 *
 * The silicon species data hold information on each species such as name,
 * diameter and the incepting particle composition.
 *
 * @param mech      Particle mechanism
 */
void SiliconInception::GenerateSpeciesData(const Sweep::Mechanism &mech)
{
    // First scan the species to determine possible candidates for inception
    // Loop over all gas-phase species
    const Sprog::SpeciesPtrVector *spec = mech.Species();
    for (unsigned int i=0; i!=spec->size(); i++) {

        string name = spec->at(i)->Name();
        CompPtrVector comp = mech.Components();

        if (IsCandidate(name)) {
            SiliconData *data;
            data = new SiliconData();
            data->_fracIndex = i;
            data->_name      = name;

            // Create a new component vector for storing new particle comp
            data->_track.resize(comp.size(), 0.0);

            // Initialise variable for volume
            double v(0.0);

            // Loop over particle components to work out how many the
            // gas-phase species of interest will  contribute
            for (unsigned int j=0; j!=comp.size(); j++) {

                if (comp[j]->Name() == "silicon") {
                    data->_track[j] = spec->at(i)->AtomCount("SI");
                    v += data->_track[j] * comp[j]->MolWt() /
                            (NA * comp[j]->Density());
                } else if (comp[j]->Name() == "hydrogen") {
                    data->_track[j] = spec->at(i)->AtomCount("H");
                    v += data->_track[j] * comp[j]->MolWt() /
                            (NA * comp[j]->Density());
                } else {
                    std::cout << "Warning! Unrecognised component name!" << std::endl;
                }
            }
            data->_diam = pow(6.0 * v / Sweep::PI, Sweep::ONE_THIRD);
            printSiliconData(*data);
            m_sidata.push_back(*data);
        }

    }


    // Now generate other data
    SetInceptingDiameter(mech);
    SetInceptingVolume(mech);

    // Resize the counter vectors
    m_reacs.resize(spec->size(), 0);
    m_concs.resize(spec->size(), 0.0);
}

const SiliconInception::SiliconData* SiliconInception::ChooseData(const Sweep::Cell &sys,
        rng_type &rng) const
{
    // Choose a species which is ABOVE THE CRITICAL DIAMETER and
    // weight the choice based on mole fractions.
    const SiliconData* ans(NULL);

    double dcrit = GetCriticalNucleus(sys);

    // Hold fractions of available species in vector
    fvector availFracs;
    std::vector<const SiliconData*> availSpecies;
    double sum(0.0);
    double frac(0.0);

    // Loop over candidate species
    for (unsigned int i=0; i != m_sidata.size(); i++) {
        if (m_sidata[i]._diam >= dcrit) {
            frac = sys.GasPhase().SpeciesConcentration(m_sidata[i]._fracIndex) /
                    sys.GasPhase().MolarDensity();
            sum += frac;
            availFracs.push_back(frac);
            availSpecies.push_back(&(m_sidata[i]));
        }
    }

    // Scale these to be probabilities
    if (sum > 0.0) {
        for (unsigned int i=0; i != availFracs.size(); i++) {
            availFracs[i] /= sum;
        }

        // Select an index based on its mole fraction
        boost::random::discrete_distribution<> dist(availFracs);
        unsigned int j = dist(rng);
        ans = availSpecies[j];
    }
    return ans;
}

/*!
 * Adjust the gas phase for the effects of this process.
 *
 *@param[in]    sys         Cell in which the process is taking place
 *@param[in]    species     Silicon hydride species to be incepted
 *@param[in]    wt          Statistical weight of particle
 *
 *@pre      The gas phase in sys must be of type SprogIdealGasWrapper
 *
 *@exception    std::runtime_error      Could not cast gas phase to SprogIdealGasWrapper
 */
void SiliconInception::adjustGasPhase(Sweep::Cell &sys,
        const SiliconInception::SiliconData &species,
        double wt) const
{
    if(!sys.FixedChem()) {
    // This method requires write access to the gas phase, which is not
    // standard in sweep.  This means it cannot use the generic gas
    // phase interface
    SprogIdealGasWrapper *gasWrapper = dynamic_cast<SprogIdealGasWrapper*>(&sys.GasPhase());
    if(gasWrapper == NULL)
        throw std::runtime_error("Could not cast gas phase to SprogIdealGasWrapper in SiliconInception::adjustGasPhase");

    // If excecution reaches here, the cast must have been successful
    Sprog::Thermo::IdealGas *gas = gasWrapper->Implementation();

    // Get the existing concentrations
    fvector newConcs;
    gas->GetConcs(newConcs);

    double n_NAvol = wt / (NA * sys.SampleVolume());

    // Update the internal counters
    m_reacs[species._fracIndex] += 1;
    m_concs[species._fracIndex] += n_NAvol;

    // Use the fracindex to calculate the amount of precursor to be
    // removed.
    newConcs[species._fracIndex] -= 1.0 * n_NAvol;

    // Now adjust the gas-phase!
    gas->SetConcs(newConcs);
    }
}

/*!
 * @brief       Does this species contribute to inception?
 *
 * @param name  Species name
 * @return      True/false
 */
bool SiliconInception::IsCandidate(std::string name) const
{
    bool ans(false);
    // Permit silylenes to incept
    if (std::string::npos != name.find("B")
            && std::string::npos != name.find("SI")) {
        ans = true;
    } else if (std::string::npos != name.find("A")
            && std::string::npos != name.find("SI")) {
        ans = true;
    } else if (name == "SI2H2") {
        ans = true;
    } else if (name == "SIH2") {
        ans = true;
    } else if (name == "SIH") {
        ans = true;
    } else if (name == "SI") {
        ans = true;
    }
    return ans;
}


/*!
 * @brief       Calculates the mole fraction of precursor
 *
 * @param sys   System for analysis
 * @return      Mole fraction of precursor
 */
double SiliconInception::GetPrecursorFraction(const Sweep::Cell &sys) const
{
    // Initialise the fraction
    double frac(0.0);

    // Loop over all gas-phase species
    for (unsigned int i=0; i!=m_sidata.size(); i++) {
        frac += sys.GasPhase().SpeciesConcentration(m_sidata[i]._fracIndex) /
                sys.GasPhase().MolarDensity();
    }

    if (frac > 1.0) {
        throw runtime_error("Mole fraction greater than 1.0! "
                            "(Sweep, SiliconInception::GetPrecursorFraction).");
    }
    return frac;
}

/*!
 * @brief       Returns the supersaturation
 *
 * Supersaturation for this system is defined by the ratio of the
 * monomer pressure to the saturation vapour pressure.
 *
 * @param sys   System for which to calculate
 * @return      Supersaturation (-)
 */
double SiliconInception::GetSupersaturation(const Sweep::Cell &sys) const
{
    double s(1.0);

    // First get the silicon pressure
    s *= (GetPrecursorFraction(sys) * sys.GasPhase().Pressure());

    // Now divide by the saturation vapour pressure
    s /= GetSatVapourPressure(sys.GasPhase().Temperature());

    return s;

}

/*!
 * @brief       Returns the surface energy of silicon
 *
 * Koermer et al (2010) J. Aerosol Sci, 41, 1007-
 *
 * @param T     Temperature (K)
 * @return      Surface energy (N/m)
 */
double SiliconInception::GetSurfaceEnergy(double T) const
{
    return 1.152 - 1.574e-4 * T;
}

/*!
 * @brief       Returns the saturation vapour pressure of silicon
 *
 * Koermer et al (2010) J. Aerosol Sci, 41, 1007-
 *
 * @param T     Temperature (K)
 * @return      Saturation vapour pressure (Pa)
 */
double SiliconInception::GetSatVapourPressure(double T) const
{
    return 101325.0 * pow(10, 7.534 - (23399 / T));
}


/*!
 * @brief       Returns the saturated monomer number concentration
 *
 * Estimates the number concentration from the saturation vapour pressure
 *
 * @param T     Temperature (K)
 * @return      Sat. monomer number concentration (#/m3)
 */
double SiliconInception::GetMonomerConc(double T) const
{
    double psat = GetSatVapourPressure(T);
    psat *= (NA / (T * R));
    return psat;
}


/*!
 * @brief       Calculates the critical nucleus size
 *
 * Uses the Kelvin equation to calculate the critical nucleus size. This
 * is given by:
 * d* = 4 gamma v1 / (kB T ln(S))
 * where:
 *  gamma is the surface energy (N/m)
 *  v1 is the volume of a monomer (m)
 *  kB is Boltzmann's constant
 *  S is the supersaturation
 *
 * @param sys   System for analysis
 * @return      Critical diameter, in metres
 */
double SiliconInception::GetCriticalNucleus(const Sweep::Cell &sys) const
{
    double d = 4.0/Sweep::KB;
    double s = GetSupersaturation(sys);

    // Get chemical conditions
    double T = sys.GasPhase().Temperature();      // in K

    // Set d to an arbitrarily high value if psi = 0 (i.e. no precursor)
    if (s <= 1.0) {
        d = 1.0;                            // in m
    } else {
        d *= (GetSurfaceEnergy(T) * m_v1 / (T * log(s)));      // in m
    }

    return d;
}


bool SiliconInception::IsInceptionAllowed(const Sweep::Cell &sys) const
{
    bool val(false);
    // Only allow inception if d_incep >= d_crit
    if (m_di >= GetCriticalNucleus(sys)) {
        val = true;
    } else{
        val = false;
    }
    return val;
}

// TOTAL RATE CALCULATIONS.

// Returns rate of the process for the given system.
double SiliconInception::Rate(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom) const
{
    // Get the current chemical conditions.
    double T = sys.GasPhase().Temperature();
    double P = sys.GasPhase().Pressure();

    // Calculate the rate.
    return Rate(sys.GasPhase(), sqrt(T),
                T/sys.GasPhase().Viscosity(), MeanFreePathAir(T,P),
                sys.SampleVolume(), sys);
}



/*!
 * Calculate inception rate using a the transition coagulation kernel
 * with the values provided by the user.  The result is this value
 * multiplied by the square of the number concentration of the gas
 * phase species.
 *
 * @param[in]    gas      interface to gas mixture
 * @param[in]    sqrtT    square root of temperature
 * @param[in]    T_mu     temperature divided by air viscosity
 * @param[in]    MFP      mean free path in gas
 * @param[in]    vol      sample volume
 *
 * @return    Inception rate for a cell of size vol. (\f$ \mathrm{s}^{-1}\f$)
 */
double SiliconInception::Rate(const EnvironmentInterface &gas, double sqrtT,
                     double T_mu, double MFP, double vol, const Cell &sys) const
{
    double rate(1.0);
    double fm(0.0);
    double sf(0.0);
    double n(0.0);
    double s(0.0);
    double y(0.0);
    double Theta(0.0);
    double T = sys.GasPhase().Temperature();

    // Check if the rate is > 0 as found by homogeneous nucleation
    if (IsInceptionAllowed(sys)) {
        switch (m_itype) {
        case iCollisional:
            // Use the 'normal' MOPS collisional rate
            rate = A() * vol * chemRatePart(gas);

            fm   = sqrtT * m_kfm;
            if((m_ksf1 > 0) || (m_ksf2 > 0))  {
                // Transition regime
                sf   = T_mu  * (m_ksf1 + (MFP*m_ksf2));
                rate *= ((fm*sf) / (fm+sf));
            }
            else {
                // Free mol regime only
                rate *= fm;
            }

            break;
        case iVBDZ:
            // Wu et al (1987) Langmuir 3:266-271
            rate = A() * vol;

            s = GetSupersaturation(sys);
            y = GetSurfaceEnergy(T);
            n = GetMonomerConc(T);

            rate *= (s * s * m_v1 * n * n);
            rate *= sqrt(2.0 * y / (PI * m_m1));
            rate *= exp(-16.0 * PI * y * y * y * m_v1 * m_v1 * m_v1
                    / (3.0 * KB * KB * KB * T * T * T * log(s) * log(s)));

            break;
        case iGirshick:
            // Girshick & Chiu (1990) J. Chem. Phys. 93:1273-1277
            rate = A() * vol;

            s = GetSupersaturation(sys);
            y = GetSurfaceEnergy(T);
            n = GetMonomerConc(T);

            Theta = PI * m_d1 * m_d1 * y / (KB * T);

            rate *= (s  * m_v1 * n * n);
            rate *= sqrt(2.0 * y / (PI * m_m1));
            rate *= exp(Theta - (4.0 * Theta * Theta * Theta / (27.0 * log(s) * log(s))));

            break;
        }
    } else {
        // We have d < d_crit, return 0.
        rate = 0.0;
    }

    return rate;
}

/*!
 * Calculates the gas-phase chemistry contribution to the rate
 * expression.  This is overloaded as Avogadro's number must be
 * included in the terms for inception processes.
 *
 * @param[in]    gas      interface to gas mixture
 */
double SiliconInception::chemRatePart(const EnvironmentInterface &gas) const
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
unsigned int SiliconInception::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a double vector. The
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
double SiliconInception::RateTerms(const double t, const Cell &sys,
                               const Geometry::LocalGeometry1d &local_geom,
                               fvector::iterator &iterm) const
{
    // Get the current chemical conditions.
    double T = sys.GasPhase().Temperature();
    double P = sys.GasPhase().Pressure();

    // Calculate the single rate term and advance iterator.
    *iterm = Rate(sys.GasPhase(), sqrt(T),
                  T/sys.GasPhase().Viscosity(), MeanFreePathAir(T,P),
                  sys.SampleVolume(), sys);
    return *(iterm++);
}


// READ/WRITE/COPY.

// Creates a copy of the inception.
SiliconInception *const SiliconInception::Clone(void) const {return new SiliconInception(*this);}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType SiliconInception::ID(void) const {return Silicon_Inception_ID;}

// Writes the object to a binary stream.
void SiliconInception::Serialize(std::ostream &out) const
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

        // Write HNT parameters
        unsigned int num = (unsigned int)m_itype;
        out.write((char*)&num, sizeof(num));
        v = (double)m_vi;
        out.write((char*)&v, sizeof(v));
        v = (double)m_di;
        out.write((char*)&v, sizeof(v));

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, SiliconInception::Serialize).");
    }
}

// Reads the object from a binary stream.
void SiliconInception::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;
        unsigned int num(0);

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

                // Read HNT parameters
                in.read(reinterpret_cast<char*>(&num), sizeof(num));
                m_itype = (InceptionType)num;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_vi = (double)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_di = (double)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, SiliconInception::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SiliconInception::Deserialize).");
    }
}
