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
double SinteringModel::A(void) const {return m_A;}

// Sets the pre-exponential constant.
void SinteringModel::SetA(double a) {m_A = a;}


// ENERGY PARAMETER.

// Returns the characteristic temperature (E) in Kelvin.
double SinteringModel::E(void) const {return m_E;}

// Sets the characteristic temperature (E) in Kelvin.
void SinteringModel::SetE(double e) {m_E = e;}

// MINIMUM PRIMARY PARTICLE DIAMETER PARAMETER.

// Returns the minimum primary particle diameter in the system (Dpmin) in m.
double SinteringModel::Dpmin(void) const {return m_dpmin;}

// Sets the minimum primary particle diameter in the system (Dpmin) in m.
void SinteringModel::SetDpmin(double dpmin) {m_dpmin = dpmin;}


// SINTERING MODEL TYPE.

// Returns the sintering model type.
SinteringModel::SintType SinteringModel::Type(void) const {return m_type;}

// Sets the sintering model type.
void SinteringModel::SetType(SinteringModel::SintType t) {m_type = t;}


// CHARACTERISTIC SINTERING TIME.

/*!
 *
 * @param temp    Temperature to evaluate at (K)
 * @param diam    Diameter of particle (m)
 * @return        Characteristic sintering time (s)
 */
double SinteringModel::SintTime(double temp, double diam) const {
	double tau_s(m_A);
	switch (m_type) {
	case ViscousFlow:
		tau_s *= diam * exp((m_E*(1-(m_dpmin/diam)))/temp);
		break;
	case GBD:
		tau_s *= diam * diam * diam * diam * temp * exp((m_E*(1-(m_dpmin/diam)))/temp);
		break;
	case SSD:
		tau_s *= diam * diam * diam * temp * exp((m_E*(1-(m_dpmin/diam)))/temp);
		break;
	case Silicon:
        // Zachariah & Carrier, J. Aerosol Sci., 1999, 30, 1139-1151
        // implementation of the SSD silicon sintering kinetic
        // form: tau = A * d^3 * T / (gamma * diff)
        //       gamma = surface energy = 1.152 - 1.574e-4*T(K)    N/m
        //       diff. = SSD diffusivity = 4.69e-7 exp(-m_E / T)   m2/s
        // default: m_A = 5396 J/K.m3
        //          m_E = 7562 K
        //          m_dpmin = 0 nm
		tau_s *= diam * diam * diam * temp / (
                (1.152 - 1.574e-4 * temp) *
                4.69e-7 * exp((m_E*(1-(m_dpmin/diam)))/temp)
                );
		break;
	case SilicaKirchoff:
        // Kirchoff et al., J. Aerosol Sci., 2012, 45, 26-39
        // http://dx.doi.org/10.1016/j.jaerosci.2011.10.006
        // form: tau = A * r_sph^2 * exp (Ea/RT)
        //  where: A = (1-exp(-1)) / (2^(2/3) * k_S)
        //         r_sph = equiv. spherical radius
        //         k_S = 3.5e-4 m2/s
        //         Ea = 3.8 eV
        //  defaults:
        //         A = 284.44 s/m2
        //         Ea = 44080 K
        // EXPECTING diam AS SPHERICAL DIAMETER
		tau_s *= diam * diam * exp(m_E / temp);
		break;
	case Rutile:
    	// Buesser et al., J. Phys. Chem. C, 2011, 115, 11030-11035
    	// SintTime function from MD calculations
    	// default: m_A = 3.7E16 s/m4K,
    	//			m_E = 258 kJ/mol = 31032 K
    	//			m_dpmin = 3.4 nm
    	tau_s *= diam * diam * diam * diam * temp *
			   exp((m_E* (1 - pow( (m_dpmin/diam) - (temp/4100.0) , 3.76))/temp));
		break;
	case Constant:
		// Do nothing.
		break;
	default:
		throw std::runtime_error("Unknown sintering model requested, "
				"Sweep::Processes::SinteringModel::SintTime.");
	}
	return tau_s;
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
                m_A = (double)val;

                // Read characteristic temperature.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_E = (double)val;

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
