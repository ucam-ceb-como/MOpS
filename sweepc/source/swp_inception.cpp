/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Inception class declared in the
    swp_inception.h header file.

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

#include "swp_inception.h"
#include "swp_mechanism.h"
#include <cmath>
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// Free-molecular enhancement factor.
const real Inception::m_efm = 2.2; // 2.2 is for soot.

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
Inception::Inception(void)
: Process(), m_a(0.5), m_kfm(0.0), m_ksf1(0.0), m_ksf2(0.0)
{
    m_name = "Inception";
}

// Initialising constructor.
Inception::Inception(const Sweep::Mechanism &mech)
: Process(mech), m_a(0.5), m_kfm(0.0), m_ksf1(0.0), m_ksf2(0.0)
{
    m_name = "Inception";
}

// Copy constructor.
Inception::Inception(const Inception &copy)
{
    *this = copy;
}

// Stream-reading constructor.
Inception::Inception(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
Inception::~Inception(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
Inception &Inception::operator =(const Inception &rhs)
{
    if (this != &rhs) {
        Process::operator =(rhs);
        m_a    = rhs.m_a;
        m_kfm  = rhs.m_kfm;
        m_ksf1 = rhs.m_ksf1;
        m_ksf2 = rhs.m_ksf2;
        m_newcomp = rhs.m_newcomp;
        m_newvals = rhs.m_newvals;
    }
    return *this;
}


// RATE CONSTANT.

// Returns the rate constant.
real Inception::A(void) const {return m_a;}

// Sets the rate constant.
void Inception::SetA(real a) {m_a = a;}


// INCEPTION KERNEL.

// Sets the coagulation kernel constants given incepting species
// masses and diameters.
void Inception::SetInceptingSpecies(real m1, real m2, real d1, real d2)
{
    // This routine sets the free-mol and slip flow kernel parameters given
    // the mass and diameter of the incepting species.
    real invd1=1.0/d1, invd2=1.0/d2;
    m_kfm  = m_efm * CFM * sqrt((1.0/m1) + (1.0/m2)) * (d1+d2) * (d1+d2);
    m_ksf1 = CSF * (d1+d2);
    m_ksf2 = 2.0 * 1.257 * m_ksf1 * ((invd1*invd1) + (invd2*invd2));
    m_ksf1 = m_ksf1 * (invd1+invd2);
}


// PROPERTIES OF INCEPTED PARTICLES.

// Returns the composition vector of the new particle.
const fvector &Inception::ParticleComp(void) const {return m_newcomp;}

// Returns the amount of the ith component of the new particle.
real Inception::ParticleComp(unsigned int i) const
{

    if (i < m_newcomp.size()) {
        return m_newcomp[i];
    } else {
        return 0.0;
    }
}

// Sets the particle composition vector.
void Inception::SetParticleComp(const fvector &comp)
{
    m_newcomp.assign(comp.begin(), comp.end());
}

// Sets the amount of the ith component in the new particle.
void Inception::SetParticleComp(unsigned int i, real comp)
{
    if (i < m_mech->ComponentCount()) {
        // Ensure vector is sufficiently long.
        if (m_newcomp.size() < m_mech->ComponentCount()) {
            m_newcomp.resize(m_mech->ComponentCount(),0.0);
        }
        // Set value.
        m_newcomp[i] = comp;
    }
}

// Returns the tracker variable vector of the new particle.
const fvector &Inception::ParticleTrackers(void) const
{
    return m_newvals;
}

// Returns the value of the ith tracker variable of the
// new particle.
real Inception::ParticleTrackers(unsigned int i) const
{
    if (i < m_newvals.size()) {
        return m_newvals[i];
    } else {
        return 0.0;
    }
}

// Sets the new particle tracker variable vector.
void Inception::SetParticleTrackers(const fvector &track)
{
    m_newvals.assign(track.begin(), track.end());
}

// Sets the value of the ith tracker variable in the
// new particle.
void Inception::SetParticleTracker(unsigned int i, real track)
{
    if (i < m_mech->TrackerCount()) {
        // Ensure vector is sufficiently long.
        if (m_newvals.size() < m_mech->TrackerCount()) {
            m_newvals.resize(m_mech->TrackerCount(),0.0);
        }
        // Set value.
        m_newvals[i] = track;
    }
}


// TOTAL RATE CALCULATIONS.

// Returns rate of the process for the given system.
real Inception::Rate(real t, const Cell &sys) const 
{
    // Get the current chemical conditions.
    real T = sys.Temperature();
    real P = sys.Pressure();

    // Calculate the rate.
    return Rate(sys.MoleFractions(), sys.Density(), sqrt(T), 
                T/ViscosityAir(T), MeanFreePathAir(T,P), 
                sys.SampleVolume());
}

/*
// Calculates the process rate using the given 
// chemical conditions, rather than those conditions in the
// given system.
real Inception::Rate(const real t, const Sprog::Thermo::IdealGas &gas, 
                     const Cell &sys) const
{
    // Get the current chemical conditions.
    real T = gas.Temperature();
    real P = gas.Pressure();

    // Calculate the rate.
    return Rate(gas.MoleFractions(), gas.Density(), sqrt(T), 
                T/ViscosityAir(T), MeanFreePathAir(T,P), 
                sys.SampleVolume());
}
*/

// Calculates the rate of multiple inceptions given a
// vector of inceptions and an iterator to a vector of
// reals for output.
real Inception::CalcRates(real t, const Cell &sys, const IcnPtrVector &icns, 
                          fvector &rates, unsigned int start)
{
    // Precalculate some values.
    real T     = sys.Temperature();
    real P     = sys.Pressure();
    real sqrtT = sqrt(T);
    real T_mu  = T / ViscosityAir(T);
    real MFP   = MeanFreePathAir(T,P);
    real vol   = sys.SampleVolume();

    IcnPtrVector::const_iterator p;
    fvector::iterator i = (rates.begin()+start);
    real sum = 0.0;
    for (p=icns.begin(); p!=icns.end(); ++p,++i) {
        *i = (*p)->Rate(sys.MoleFractions(), sys.Density(), sqrtT, 
                        T_mu, MFP, vol);
        sum += *i;
    }
    return sum;
}

// A faster rate calculation routine for Inception events only.  Requires all the
// parameters that would otherwise be calculated by the routine to be passed as
// arguments.
real Inception::Rate(const fvector &fracs, real density, real sqrtT, 
                     real T_mu, real MFP, real vol) const
{
    // Temperature and pressure dependence.
    real fm   = sqrtT * m_kfm;
    real sf   = T_mu  * (m_ksf1 + (MFP*m_ksf2));
    real rate = m_a   * ((fm*sf) / (fm+sf)) * vol;

    // Chemical species concentration dependence.
    rate *= chemRatePart(fracs, density);

    return rate;
}

// Calculates the gas-phase chemistry contribution to the rate
// expression.  This is overloaded as Avogadro's number must be
// included in the terms for inception processes.
real Inception::chemRatePart(const fvector &fracs, real density) const
{
    real rate = 1.0;

    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        real conc = density * fracs[i->first];
        for (int j=0; j!=i->second; ++j) {
            rate *= (NA * conc);
        }
    }

    return rate;
}


// RATE TERM CALCULATIONS.

// Returns the number of rate terms for this process (one).
unsigned int Inception::TermCount(void) const {return 1;}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all terms.
real Inception::RateTerms(const real t, const Cell &sys, 
                          fvector::iterator &iterm) const
{
    // Get the current chemical conditions.
    real T = sys.Temperature();
    real P = sys.Pressure();

    // Calculate the single rate term and advance iterator.
    *iterm = Rate(sys.MoleFractions(), sys.Density(), sqrt(T), 
                  T/ViscosityAir(T), MeanFreePathAir(T,P), 
                  sys.SampleVolume());
    return *(iterm++);
}

/*
// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  The given chemical conditions are used instead of those
// in the given system object.
real Inception::RateTerms(const real t, const Sprog::Thermo::IdealGas &gas,
                          const Cell &sys, fvector::iterator &iterm) const
{
    // Get the current chemical conditions.
    real T = gas.Temperature();
    real P = gas.Pressure();

    // Calculate rate term and advance iterator.
    *iterm = Rate(gas.MoleFractions(), gas.Density(), sqrt(T), 
                  T/ViscosityAir(T), MeanFreePathAir(T,P), 
                  sys.SampleVolume());
    return *(iterm++);
}
*/

// PERFORMING THE PROCESS.

// Performs the process on the given system.  The responsible rate term is given
// by index.  Returns 0 on success, otherwise negative.
int Inception::Perform(real t, Cell &sys, unsigned int iterm, TransportOutflow*) const
{

    // This routine performs the inception on the given chemical system.

    // Create a new particle of the type specified
    // by the system ensemble.
    Particle *sp = m_mech->CreateParticle(t);
    
    // Initialise the new particle.
    sp->Primary()->SetComposition(m_newcomp);
    sp->Primary()->SetValues(m_newvals);
    sp->UpdateCache();

    // Add particle to system's ensemble.
    sys.Particles().Add(*sp);

    // Update gas-phase chemistry of system.
    adjustGas(sys);

    return 0;
}


// READ/WRITE/COPY.

// Creates a copy of the inception.
Inception *const Inception::Clone(void) const {return new Inception(*this);}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType Inception::ID(void) const {return Inception_ID;}

// Writes the object to a binary stream.
void Inception::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        Process::Serialize(out);

        // Write rate constant.
        double v = (double)m_a;
        out.write((char*)&v, sizeof(v));

        // Write free-mol parameter.
        v = (double)m_kfm;
        out.write((char*)&v, sizeof(v));

        // Write slip-flow parameters.
        v = (double)m_ksf1;
        out.write((char*)&v, sizeof(v));
        v = (double)m_ksf2;
        out.write((char*)&v, sizeof(v));

        // Write new component count.
        unsigned int n = (unsigned int)m_newcomp.size();
        out.write((char*)&n, sizeof(n));

        // Write component values.
        for (fvector::const_iterator i=m_newcomp.begin(); i!=m_newcomp.end(); ++i) {
            v = (double)*i;
            out.write((char*)&v, sizeof(v));
        }

        // Write new tracker values count.
        n = (unsigned int)m_newvals.size();
        out.write((char*)&n, sizeof(n));

        // Write tracker values.
        for (fvector::const_iterator i=m_newvals.begin(); i!=m_newvals.end(); ++i) {
            v = (double)*i;
            out.write((char*)&v, sizeof(v));
        }
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Inception::Serialize).");
    }
}

// Reads the object from a binary stream.
void Inception::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    m_newcomp.clear();
    m_newvals.clear();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;
        double val = 0.0;

        switch (version) {
            case 0:
                // Deserialize base class.
                Process::Deserialize(in, mech);

                // Read rate constant.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_a = (real)val;

                // Read free-mol parameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_kfm = (real)val;

                // Read slip-flow parameters.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_ksf1 = (real)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_ksf2 = (real)val;

                // Read new component count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read component values.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_newcomp.push_back((real)val);
                }

                // Read new tracker values count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read tracker values.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_newvals.push_back((real)val);
                }
                
                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Inception::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Inception::Deserialize).");
    }
}
