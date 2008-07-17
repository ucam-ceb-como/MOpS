/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Coagulation class declared in the
    swp_coagulation.h header file.

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

#include "swp_coagulation.h"
#include "swp_mechanism.h"
#include "swp_ensemble.h"
#include "swp_particle_cache.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// Free-molecular enhancement factor.
const real Coagulation::m_efm = 2.2; // 2.2 is for soot.

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Coagulation::Coagulation(const Sweep::Mechanism &mech)
: Process(mech)
{
    m_name = "Coagulation";
}

// Copy contructor.
Coagulation::Coagulation(const Coagulation &copy)
{
    *this = copy;
}

// Stream-reading constructor.
Coagulation::Coagulation(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
Coagulation::~Coagulation(void)
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator.
Coagulation &Coagulation::operator=(const Coagulation &rhs)
{
    if (this != &rhs) {
        Process::operator=(rhs);
    }
    return *this;
}


// TOTAL RATE CALCULATION.

// Returns the rate of the process for the given system.
real Coagulation::Rate(real t, const Cell &sys) const 
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

/*
// Calculates the process rate using the given 
// chemical conditions, rather than those conditions in the
// given system.
real Coagulation::Rate(real t, const Sprog::Thermo::IdealGas &gas, 
                       const Cell &sys) const
{
    // Get the number of particles in the system.
    unsigned int n = sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        real T = gas.Temperature();
        real P = gas.Pressure();

        // Calculate the rate.
        return Rate(sys.Particles().GetSums(), (real)n, sqrt(T), 
                    T/ViscosityAir(T), MeanFreePathAir(T,P),
                    sys.SampleVolume());
    } else {
        return 0.0;
    }
}
*/

// More efficient rate routine for coagulation only.  
// All parameters required to calculate rate passed 
// as arguments.
real Coagulation::Rate(const ParticleCache &data, real n, real sqrtT, 
                       real T_mu, real MFP, real vol) const
{
    // Some prerequisites.
    real n_1 = n - 1.0;
    real a = CSF * T_mu;
    real b = a * MFP * 1.257;
    real c = CFMMAJ * m_efm * CFM * sqrtT;

    // Summed particle properties required for coagulation rate.
    real d       = data.CollDiameter();
    real d2      = data.CollDiamSquared();
    real d_1     = data.InvCollDiam();
    real d_2     = data.InvCollDiamSquared();
    real d2m_1_2 = data.CollDiamSqrdInvSqrtMass();
    real m_1_2   = data.InvSqrtMass();

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
unsigned int Coagulation::TermCount(void) const {return TYPE_COUNT;}

// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  Returns the sum of all rate terms.
real Coagulation::RateTerms(real t, const Cell &sys, 
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

/*
// Calculates the rate terms given an iterator to a real vector. The 
// iterator is advanced to the position after the last term for this
// process.  The given chemical conditions are used instead of those
// in the given system object.  Returns the sum of all rate terms.
real Coagulation::RateTerms(real t, const Sprog::Thermo::IdealGas &gas, 
                            const Cell &sys, fvector::iterator &iterm) const
{
    // Get the number of particles in the system.
    unsigned int n = sys.ParticleCount();

    // Check that there are at least 2 particles before calculating rate.
    if (n > 1) {
        // Get system properties required to calculate coagulation rate.
        real T = gas.Temperature();
        real P = gas.Pressure();

        // Calculate the rate terms.
        return RateTerms(sys.Particles().GetSums(), (real)n, sqrt(T), T/ViscosityAir(T), 
                         MeanFreePathAir(T,P), sys.SampleVolume(), iterm);
    } else {
        // No coagulation as there are too few particles.
        for(unsigned int i=0; i!=TYPE_COUNT; ++i) *(iterm++) = 0.0;
        return 0.0;
    }
}
*/

// More efficient rate routine for coagulation only.  
// All parameters required to calculate rate terms
// passed as arguments.
real Coagulation::RateTerms(const ParticleCache &data, real n, real sqrtT, 
                            real T_mu, real MFP, real vol,
                            fvector::iterator &iterm) const
{
    // Some prerequisites.
    real n_1 = n - 1.0;
    real a   = CSF * T_mu;
    real b   = a * MFP * 1.257 * 2.0;
    real c   = CFMMAJ * m_efm * CFM * sqrtT;

    // Summed particle properties required for coagulation rate.
    real d       = data.CollDiameter();
    real d2      = data.CollDiamSquared();
    real d_1     = data.InvCollDiam();
    real d_2     = data.InvCollDiamSquared();
    real d2m_1_2 = data.CollDiamSqrdInvSqrtMass();
    real m_1_2   = data.InvSqrtMass();

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

// Performs the process on the given system. Must return 0
// on success, otherwise negative.
int Coagulation::Perform(real t, Cell &sys, unsigned int iterm) const
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
            ip1 = sys.Particles().Select();
            maj = SlipFlow;
            break;
        case SlipFlow2:
            ip1 = sys.Particles().Select(ParticleCache::iD);
            maj = SlipFlow;
            break;
        case SlipFlow3:
            ip1 = sys.Particles().Select();
            maj = SlipFlow;
            break;
        case SlipFlow4:
            ip1 = sys.Particles().Select(ParticleCache::iD);
            maj = SlipFlow;
            break;
        case FreeMol1:
            ip1 = sys.Particles().Select();
            maj = FreeMol;
            break;
        case FreeMol2:
            ip1 = sys.Particles().Select(ParticleCache::iD2);
            maj = FreeMol;
            break;
        default :
            ip1 = sys.Particles().Select();
            maj = SlipFlow;
            break;
    }

    // Save current T and P.
    real T = sys.Temperature();
    real P = sys.Pressure();

    // Choose and get first particle, then update it.
    Particle *sp1=NULL, *sp1old=NULL;
    if (ip1 >= 0) {
        sp1 = sys.Particles().At(ip1);

        // Create a copy of the particle before updating.
        sp1old = sp1->Clone();

        m_mech->UpdateParticle(*sp1, sys, t);
    } else {
        // Failed to choose a particle.
        return -1;
    }

    // Check that particle is still valid.  If not,
    // remove it and cease coagulating.
    if (!sp1->IsValid()) {
        // Must remove first particle now.
        sys.Particles().Remove(ip1);

        // Invalidating the index tells this routine not to perform coagulation.
        ip1 = -1;
    }

    // Choose and get unique second particle, then update it.  Note, we are allowed to do
    // this even if the first particle was invalidated.
    ip2 = ip1;
    unsigned int guard = 0;
    switch (term) {
        case SlipFlow1:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select();            
            break;
        case SlipFlow2:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(ParticleCache::iD_1);
            break;
        case SlipFlow3:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(ParticleCache::iD_1);
            break;
        case SlipFlow4:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(ParticleCache::iD_2);
            break;
        case FreeMol1:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(ParticleCache::iD2_M_1_2);
            break;
        case FreeMol2:
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select(ParticleCache::iM_1_2);
            break;
        default :
            while ((ip2 == ip1) && (++guard<1000))
                ip2 = sys.Particles().Select();
            break;
    }

    Particle *sp2=NULL, *sp2old=NULL;
    if ((ip2>=0) && (ip2!=ip1)) {
        sp2 = sys.Particles().At(ip2);

        // Create a copy of the particle before updating.
        sp2old = sp2->Clone();

        m_mech->UpdateParticle(*sp2, sys, t);
    } else {
        // Failed to select a unique particle.
        delete sp1old;
        return -1;
    }

    // Check validity of particles after update.
    if (!sp2->IsValid()) {
        // Must remove second particle now.
        sys.Particles().Remove(ip2);

        // Invalidating the index tells this routine not to perform coagulation.
        ip2 = -1;
    }

    // Check that both the particles are still valid.
    if ((ip1>=0) && (ip2>=0)) {
        // Must check for ficticious event now by calculate the original
        // majorant rate and the current (after updates) true rate.
        real majk = CoagKernel(*sp1old, *sp2old, T, P, maj);
        real truek = CoagKernel(*sp1, *sp2, T, P, None);
        if (!Ficticious(majk, truek)) {
            // We can now coagulate the particles, remember to
            // remove second particle afterwards.
            if (ip1 < ip2) {
                *sp1 += *sp2;
                sp1->SetTime(t);
                sys.Particles().Update(ip1);
                sys.Particles().Remove(ip2, !m_mech->UseSubPartTree());
            } else {
                *sp2 += *sp1;
                sp2->SetTime(t);
                sys.Particles().Update(ip2);
                sys.Particles().Remove(ip1, !m_mech->UseSubPartTree());
            }
        } else {
            delete sp1old; delete sp2old;
            return 1; // Ficticious event.
        }
    } else {
        // One or both particles were invalidated on update,
        // but that's not a problem.
    }

    delete sp1old; delete sp2old;
    return 0;
}


// COAGULATION KERNELS.

// Returns the transition coagulation kernel value for the 
// two given particles.
real Coagulation::CoagKernel(const Particle &sp1, const Particle &sp2, 
                             real T, real P, MajorantType maj) const
{
    // This routine calculates the coagulation kernel for two particles.  The kernel
    // type is chosen by the majorant type requested.
    real fm=0.0, sf=0.0;
    switch (maj) {
        case None:
            // No majorant, so return half harmonic mean of free molecular and
            // slip-flow kernels (non-majorant).
            fm = FreeMolKernel(sp1, sp2, T, P, false);
            sf = SlipFlowKernel(sp1, sp2, T, P, false);
            return (fm*sf)/(fm+sf);
        case FreeMol:
            // Free molecular majorant.
            return FreeMolKernel(sp1, sp2, T, P, true);
        case SlipFlow:
            // Slip-flow majorant.
            return SlipFlowKernel(sp1, sp2, T, P, true);
    }

    // Invalid majorant, return zero.
    return 0.0;
}

// Returns the free-molecular coagulation kernel value for the 
// two given particles.  Can return either the majorant or
// true kernel.
real Coagulation::FreeMolKernel(const Particle &sp1, const Particle &sp2, 
                                real T, real P, bool maj) const
{
    // This routine calculate the free molecular coagulation kernel for two particles.
    // There are two forms of kernel; a majorant form and a non-majorant form.

    if (maj) {
        // The majorant form is always >= the non-majorant form.
        return CFMMAJ * m_efm * CFM * sqrt(T) * 
               (sp1.InvSqrtMass() + sp2.InvSqrtMass()) *
               (sp1.CollDiamSquared() + sp2.CollDiamSquared());
    } else {
        real dterm = sp1.CollDiameter()+sp2.CollDiameter();
        return m_efm * CFM * 
               sqrt(T * ((1.0/sp1.Mass())+(1.0/sp2.Mass()))) * 
               dterm * dterm;
    }
}

// Returns the slip-flow coagulation kernel value for the 
// two given particles.  Can return either the majorant or
// true kernel.  This kernel is currently evaluated assuming 
// that the surrounding gas is air, i.e. principally N2.
real Coagulation::SlipFlowKernel(const Particle &sp1, const Particle &sp2, 
                                 real T, real P, bool maj) const
{
    // For the slip-flow kernel the majorant and non-majorant forms are identical.
    return ((1.257 * 2.0 * MeanFreePathAir(T,P) * 
             (sp1.InvCollDiamSquared() + sp2.InvCollDiamSquared())) +
            (sp1.InvCollDiam() + sp2.InvCollDiam())) * 
           CSF * T * (sp1.CollDiameter()+sp2.CollDiameter()) / ViscosityAir(T);
}


// READ/WRITE/COPY.

// Creates a copy of the coagulation process.
Coagulation *const Coagulation::Clone(void) const
{
    return new Coagulation(*this);
}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType Coagulation::ID(void) const {return Coagulation_ID;}

// Writes the object to a binary stream.
void Coagulation::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        Process::Serialize(out);

    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Coagulation::Serialize).");
    }
}

// Reads the object from a binary stream.
void Coagulation::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                // Deserialize base class.
                Process::Deserialize(in, mech);

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Coagulation::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Coagulation::Deserialize).");
    }
}
