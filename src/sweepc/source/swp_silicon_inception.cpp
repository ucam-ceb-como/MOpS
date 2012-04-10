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
#include <boost/random/uniform_01.hpp>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

//! Free-molecular enhancement factor.
const real SiliconInception::m_efm = 2.2; // 2.2 is for soot.

// CONSTRUCTORS AND DESTRUCTORS.

//! Default constructor (protected).
SiliconInception::SiliconInception(void)
: Inception(), m_kfm(0.0), m_ksf1(0.0), m_ksf2(0.0), m_v1(0.0), m_di(0.0)
{
    m_name = "SiliconInception";
}

//! Initialising constructor.
SiliconInception::SiliconInception(const Sweep::Mechanism &mech)
: Inception(mech), m_kfm(0.0), m_ksf1(0.0), m_ksf2(0.0), m_v1(0.0), m_di(0.0)
{
    m_name = "SiliconInception";
}

//! Copy constructor.
SiliconInception::SiliconInception(const SiliconInception &copy)
{
    *this = copy;
}

//! Stream-reading constructor.
SiliconInception::SiliconInception(std::istream &in, const Sweep::Mechanism &mech)
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
        m_v1   = rhs.m_v1;
        m_di   = rhs.m_di;
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
int SiliconInception::Perform(const real t, Cell &sys,
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
    real posn = vertices.front();

    const real width = vertices.back() - posn;
    boost::uniform_01<rng_type&, real> uniformGenerator(rng);
    posn += width * uniformGenerator();

    sp->setPositionAndTime(posn, t);


    // Initialise the new particle.
    sp->Primary()->SetComposition(ParticleComp());
    sp->Primary()->SetValues(ParticleTrackers());
    sp->UpdateCache();

    // Add particle to system's ensemble.
    sys.Particles().Add(*sp, rng);

    // Update gas-phase chemistry of system.
    adjustGas(sys, sp->getStatisticalWeight());

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
void SiliconInception::SetInceptingSpecies(real m1, real m2, real d1, real d2)
{
    // The free mol part can be handled by the free mol specific method
    SetInceptingSpeciesFreeMol(m1, m2, d1, d2);

    // Now the slip flow part
    real invd1=1.0/d1, invd2=1.0/d2;
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
void SiliconInception::SetInceptingSpeciesFreeMol(real m1, real m2, real d1, real d2)
{
    // This routine sets the free-mol and slip flow kernel parameters given
    // the mass and diameter of the incepting species.
    m_kfm  = m_efm * CFM * sqrt((1.0/m1) + (1.0/m2)) * (d1+d2) * (d1+d2);
    m_ksf1 = 0.0;
    m_ksf2 = 0.0;
}

/*!
 * @brief           Sets the volume of a monomer particle
 *
 * The volume of a monomer particle (1 of each component) is needed for the
 * critical nucleus calculation. This is determined upon initialisation, based
 * on the density of molwt of each component.
 *
 * @param mech      Particle mechanism
 */
void SiliconInception::SetMonomerVolume(const Sweep::Mechanism &mech)
{
    real v1(0.0);
    CompPtrVector comp = mech.Components();
    for (unsigned int i=0; i!=comp.size(); i++) {
        // returns in cm3, due to Density()
        v1 += 1.0 * mech.Components(i)->MolWt() /
                (Sweep::NA * mech.Components(i)->Density());
    }

    m_v1 = v1;
}

/*!
 * @brief       Calculates the diameter of an incepting particle
 *
 * Pre-calculates the incepting diameter of a particle for this inception
 * reaction. This is done using the molwt, density and inital number of
 * each component in the new particle.
 *
 * @param mech  Particle mechanism
 */
void SiliconInception::SetInceptingDiameter(const Sweep::Mechanism &mech)
{
    //Pre-define volume and diameters
    real v(0.0), d(0.0);

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


/*!
 * @brief       Calculates the mole fraction of precursor
 *
 * @param sys   System for analysis
 * @return      Mole fraction of precursor
 */
real SiliconInception::GetPrecursorFraction(const Sweep::Cell &sys) const
{
    // Initialise the fraction
    real frac(0.0);

    // Loop over all gas-phase species
    for (unsigned int i=0; i!=sys.GasPhase().Species()->size(); i++) {
        if (sys.GasPhase().Species()->at(i)->ContainsElement("SI")) {
            string name = sys.GasPhase().Species()->at(i)->Name();
            // First check for 'B' string (silylene)
            if (std::string::npos != name.find("B")) {
                frac += sys.GasPhase().MoleFraction(i);
            } else if (name == "SIH2") {
                frac += sys.GasPhase().MoleFraction(i);
            } else {
                frac += 0.0;
            }

        }
    }

    if (frac > 1.0) {
        throw runtime_error("Mole fraction greater than 1.0! "
                            "(Sweep, SiliconInception::GetPrecursorFraction).");
    }
    return frac;
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
real SiliconInception::GetCriticalNucleus(const Sweep::Cell &sys) const
{
    real d = 4.0/Sweep::KB;

    // Get chemical conditions
    real T = sys.GasPhase().Temperature();      // in K
    real P = sys.GasPhase().Pressure();         // in Pa

    // Calculate the surface energy of silicon
    // Koermer et al (2010) J. Aerosol Sci, 41, 1007-
    real gamma = 1.152 - 1.574e-4 * T;

    // Calculate the supersaturation factor, S
    // first get the saturation vapour pressure of silicon
    real psat = pow(10, 7.534 - (23399 / T));       // in atm
    psat *= 101325;                                 // in Pa
    // now get the pressure of the contributing silicon hydride species
    real psi = GetPrecursorFraction(sys) * P;

    // Set d to an arbitrarily high value if psi = 0 (i.e. no precursor)
    // avoid getting log(0)
    if (psi == 0.0) {
        d = 1.0;                            // in m
    } else {
        // now calculate the supersaturation
        real s = log(P / psat);

        d *= (gamma * m_v1 / (T * s));      // in m
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
real SiliconInception::Rate(real t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom) const
{
    // Get the current chemical conditions.
    real T = sys.GasPhase().Temperature();
    real P = sys.GasPhase().Pressure();

    // Calculate the rate.
    return Rate(sys.GasPhase().MoleFractions(), sys.GasPhase().Density(), sqrt(T),
                T/ViscosityAir(T), MeanFreePathAir(T,P),
                sys.SampleVolume(), sys);
}



/*!
 * Calculate inception rate using a the transition coagulation kernel
 * with the values provided by the user.  The result is this value
 * multiplied by the square of the number concentration of the gas
 * phase species.
 *
 * @param[in]    fracs    species molefractions
 * @param[in]    density  gas number density in \f$ \mathrm{mol}\ \mathrm{m}^{-3}\f$
 * @param[in]    sqrtT    square root of temperature
 * @param[in]    T_mu     temperature divided by air viscosity
 * @param[in]    MFP      mean free path in gas
 * @param[in]    vol      sample volume
 *
 * @return    Inception rate for a cell of size vol. (\f$ \mathrm{s}^{-1}\f$)
 */
real SiliconInception::Rate(const fvector &fracs, real density, real sqrtT,
                     real T_mu, real MFP, real vol, const Cell &sys) const
{
    real rate(1.0);

    // Check if the rate is > 0 as found by homogeneous nucleation
    if (IsInceptionAllowed(sys)) {
        real rate = A() * vol * chemRatePart(fracs, density);

        const real fm   = sqrtT * m_kfm;
        if((m_ksf1 > 0) || (m_ksf2 > 0))  {
            // Transition regime
            real sf   = T_mu  * (m_ksf1 + (MFP*m_ksf2));
            rate *= ((fm*sf) / (fm+sf));
        }
        else {
            // Free mol regime only
            rate *= fm;
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
 * @param[in]    fracs    species molefractions
 * @param[in]    density  gas number density in \f$ \mathrm{mol}\ \mathrm{m}^{-3}\f$
 */
real SiliconInception::chemRatePart(const fvector &fracs, real density) const
{
    // Factor of 0.5 adjusts for doubling counting of pairs of molecules in the number
    // of possible collisions.
    real rate = 0.5;

    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        //std::cerr << "Mole frac to use " << fracs[i->first] << std::endl;
        real conc = density * fracs[i->first];
        for (int j=0; j!=i->second; ++j) {
            rate *= (NA * conc);
        }
    }

    return rate;
}


// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int SiliconInception::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a real vector. The
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
real SiliconInception::RateTerms(const real t, const Cell &sys,
                               const Geometry::LocalGeometry1d &local_geom,
                               fvector::iterator &iterm) const
{
    // Get the current chemical conditions.
    real T = sys.GasPhase().Temperature();
    real P = sys.GasPhase().Pressure();

    // Calculate the single rate term and advance iterator.
    *iterm = Rate(sys.GasPhase().MoleFractions(), sys.GasPhase().Density(), sqrt(T),
                  T/ViscosityAir(T), MeanFreePathAir(T,P),
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
        v = (double)m_v1;
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

        switch (version) {
            case 0:
                // Deserialize base class.
                Inception::Deserialize(in, mech);

                // Read free-mol parameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_kfm = (real)val;

                // Read slip-flow parameters.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_ksf1 = (real)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_ksf2 = (real)val;

                // Read HNT parameters
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_v1 = (real)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_di = (real)val;

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
