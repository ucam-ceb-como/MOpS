/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the SinteringModel class declared in the
    swp_sintering_model.h header file.

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

#include "swp_sintering_model.h"
#include "swp_primary.h"
#include "swp_particle.h"
#include "swp_cell.h"
#include "swp_subparticle.h"
#include <cmath>
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
SinteringModel::SinteringModel()
: m_enable(false), m_A(0.0), m_E(0.0), m_dpmin(0.0), m_type(GBD)
{
}

// Copy constructor.
SinteringModel::SinteringModel(const Sweep::Processes::SinteringModel &copy)
{
    *this = copy;
}

// Stream-reading constructor.
SinteringModel::SinteringModel(std::istream &in)
{
    Deserialize(in);
}

// Destructor.
SinteringModel::~SinteringModel()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADING.

// Assignment operator.
SinteringModel &SinteringModel::operator=(const Sweep::Processes::SinteringModel &rhs)
{
    if (this != &rhs) {
        m_enable = rhs.m_enable;
        m_A      = rhs.m_A;
        m_E      = rhs.m_E;
        m_type   = rhs.m_type;
    }
    return *this;
}


// MODEL ENABLE/DISABLE.

// Returns true is the model is enabled, otherwise false.
bool SinteringModel::IsEnabled(void) const {return m_enable;}

// Enables this sintering model.
void SinteringModel::Enable(void) {m_enable = true;}

// Disables this sintering model.
void SinteringModel::Disable(void) {m_enable = false;}


// PRE-EXPONENTIAL CONSTANT.

// Returns the pre-exponential constant.
real SinteringModel::A(void) const {return m_A;}

// Sets the pre-exponential constant.
void SinteringModel::SetA(real a) {m_A = a;}


// ENERGY PARAMETER.

// Returns the characteristic temperature (E) in Kelvin.
real SinteringModel::E(void) const {return m_E;}

// Sets the characteristic temperature (E) in Kelvin.
void SinteringModel::SetE(real e) {m_E = e;}

// MINIMUM PRIMARY PARTICLE DIAMETER PARAMETER.

// Returns the minimum primary particle diameter in the system (Dpmin) in m.
real SinteringModel::Dpmin(void) const {return m_dpmin;}

// Sets the minimum primary particle diameter in the system (Dpmin) in m.
void SinteringModel::SetDpmin(real dpmin) {m_dpmin = dpmin;}


// SINTERING MODEL TYPE.

// Returns the sintering model type.
SinteringModel::SintType SinteringModel::Type(void) const {return m_type;}

// Sets the sintering model type.
void SinteringModel::SetType(SinteringModel::SintType t) {m_type = t;}


// CHARACTERISTIC SINTERING TIME.

// Returns the characteristic sintering time for the
// given particle.
real SinteringModel::SintTime(const Cell &sys, const Particle &p) const
{
    if (p.Primary() != NULL) {
        return SintTime(sys, *p.Primary());
    } else {
        // TODO:  Complete SinteringModel::SintTime() function for sub-particle tree.    
		// This is included in the function SinteringModel::SintTime(const Cell &sys, const SubParticle &p) const
        return 0.0;
    }
}

// Returns the characteristic sintering time for the
// given primary.
real SinteringModel::SintTime(const Cell &sys,const Primary &p) const
{
    real dp = 6.0 * p.Volume() / p.SurfaceArea();
    switch (m_type) {
        case ViscousFlow:
            return m_A * dp * 
                   exp((m_E*(1-(m_dpmin/dp)))/sys.Temperature());
            break;
        case GBD:
        default:
            return m_A * dp * dp * dp * dp * sys.Temperature() *  
                   exp((m_E*(1-(m_dpmin/dp)))/sys.Temperature());
        case Rutile:
        	// Buesser et al., J. Phys. Chem. C, 2011, 115, 11030-11035
        	// SintTime function from MD calculations
        	// default: m_A = 3.7E16 s/m4K,
        	//			m_E = 258 kJ/mol = 31032 K
        	//			m_dpmin = 3.4 nm
        	return m_A * dp * dp * dp * dp * sys.Temperature() *
				   exp((m_E* (1 - pow( (m_dpmin/dp) - (sys.Temperature()/4100.0) , 3.76))
				     /sys.Temperature()));
    }
}

// RATE CALCULATION.

// Returns the rate of the process for the given particle.
// wjm34: include max() catch on tau to prevent getting stuck in infinite Sinter() loops
real SinteringModel::Rate(real t, const Cell &sys, const Particle &p) const
{
    real tau = max(1.0e-30, SintTime(sys, p));
    return (p.SurfaceArea() - p.SphSurfaceArea()) / tau;
}

// Returns the rate of the process for the given primary.
real SinteringModel::Rate(real t, const Cell &sys, const Primary &p) const
{
    real tau = max(1.0e-30, SintTime(sys, p));
    return max((p.SurfaceArea() - p.SphSurfaceArea()) / tau, 0.0);
}

// READ/WRITE/COPY.

// Creates a copy of the sintering model.
SinteringModel *const SinteringModel::Clone(void) const
{
    return new SinteringModel(*this);
}

// Writes the object to a binary stream.
void SinteringModel::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write if enabled or disabled model.
        if (m_enable) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write pre-exponential factor.
        double val = (double)m_A;
        out.write((char*)&val, sizeof(val));

        // Write characteristic temperature.
        val = (double)m_E;
        out.write((char*)&val, sizeof(val));

        // Write type.
        unsigned int t = (unsigned int)m_type;
        out.write((char*)&t, sizeof(t));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, SinteringModel::Serialize).");
    }
}

// Reads the object from a binary stream.
void SinteringModel::Deserialize(std::istream &in)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n=0;
        double val=0.0;

        switch (version) {
            case 0:
                // Read if enabled
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_enable = (n==1);

                // Read pre-exponential factor.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_A = (real)val;

                // Read characteristic temperature.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_E = (real)val;

                // Read type.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_type = (SintType)n;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, SinteringModel::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, SinteringModel::Deserialize).");
    }
}
