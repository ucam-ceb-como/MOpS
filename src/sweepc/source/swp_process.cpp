/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Process class declared in the
    swp_process.h header file.

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

#include "swp_process.h"

#include "swp_params.h"
#include "swp_mechanism.h"
#include "swp_sprog_idealgas_wrapper.h"

#include <stdexcept>
#include <cmath>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Process::Process(void)
: m_name(""), m_mech(NULL), m_a(1.0)
{
}

// Initialising constructor.
Process::Process(const Sweep::Mechanism &mech)
: m_name(""), m_mech(&mech), m_a(1.0)
{
}

// Copy constructor.
Process::Process(const Process &copy)
{
    *this = copy;
}

// Default destructor.
Process::~Process(void)
{
}


// OPERATOR OVERLOADS.

// Assigment operator.
Process &Process::operator=(const Process &rhs)
{
    if (this != &rhs) {
        m_name = rhs.m_name;
        m_mech = rhs.m_mech;

        m_a = rhs.m_a;

        // Copy reactants.
        Sprog::StoichMap::const_iterator i;
        for (i=rhs.m_reac.begin(); i!=rhs.m_reac.end(); ++i) {
            m_reac[i->first] = i->second;
        }

        // Copy products.
        for (i=rhs.m_prod.begin(); i!=rhs.m_prod.end(); ++i) {
            m_prod[i->first] = i->second;
        }
    }
    return *this;
}


// PROCESS NAME/DESCRIPTION.

// Returns the process name.
const std::string &Process::Name(void) const {return m_name;}

// Sets the process name.
void Process::SetName(const std::string &name) {m_name = name;}


// PARENT MECHANISM.

// Returns reference to parent mechanism.
const Sweep::Mechanism *const Process::Mechanism() const
{
    return m_mech;
}

// Sets the parent mechanism
void Process::SetMechanism(const Sweep::Mechanism &mech)
{
    m_mech = &mech;
}

// REACTANTS.

// Returns the reactant count.
unsigned int Process::ReactantCount() const
{
    return m_reac.size();
}

// Returns the stoichiometric reactant coefficients.
const Sprog::StoichMap &Process::Reactants(void) const
{
    return m_reac;
}

// Returns the stoichiometry of the ith reactant.
int Process::Reactants(unsigned int i) const
{
    Sprog::StoichMap::const_iterator j = m_reac.find(i);
    if (j != m_reac.end()) {
        return j->second;
    } else {
        return 0;
    }
}

// Adds a reactant to the reaction.
void Process::AddReactant(unsigned int isp, int mu)
{
    m_reac[isp] = mu;
}

// Adds a reactant given the species name.
void Process::AddReactant(const std::string &name, int mu)
{
    // Locate the species.
    unsigned int i;
    for (i=0; i!=m_mech->Species()->size(); ++i) {
        if (name.compare((*m_mech->Species())[i]->Name()) == 0) {
            // Set reactant.
            m_reac[i] = mu;
            return;
        }
    }
}

// Removes a reactant, given by name, from the reaction.
void Process::RemoveReactant(const std::string &name)
{
    // Locate the species.
    unsigned int i;
    for (i=0; i!=m_mech->Species()->size(); ++i) {
        if (name.compare((*m_mech->Species())[i]->Name()) == 0) {
            // Delete reactant.
            m_reac.erase(i);
            return;
        }
    }
}


// PRODUCTS.

// Returns the product count.
unsigned int Process::ProductCount() const
{
    return m_prod.size();
}

// Returns the stoichiometric product coefficients.
const Sprog::StoichMap &Process::Products(void) const
{
    return m_prod;
}

// Returns the stoichiometry of the ith product.
int Process::Products(unsigned int i) const
{
    Sprog::StoichMap::const_iterator j = m_prod.find(i);
    if (j != m_prod.end()) {
        return j->second;
    } else {
        return 0;
    }
}

// Adds a product to the reaction.
void Process::AddProduct(unsigned int isp, int mu)
{
    m_prod[isp] = mu;
}

/*!
 * Virtual parent class function.
 *
 * @param t     Current time of the system
 * @param dt    Time to remove particles over
 * @param sys   The system to do transport for
 * @param local_geom    Geometry of the process
 * @param rng   Random number generator
 */
void Process::PerformDT (
        const double t,
        const double dt,
        Sweep::Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        rng_type &rng) const {
    assert(dt > 0.0);
}

// FICTICIOUS EVENTS.

/*!
 * @param[in]       majr        Majorant rate for event
 * @param[in]       truer       True rate for event
 * @param[in,out]   rng         Random number generator
 *
 * @pre     truer <= majr, if truer > majr, cerr a warning and return false
 *
 * @return      True with probability 1-truer/majr, otherwise false
 */
bool Process::Fictitious(double majr, double truer, rng_type &rng)
{
    // ensure truer <= maj, otherwise it will crash the program.
    if (truer>majr)
    {
        // if the following warning msg show frequently, a better solution ensuring truer <= maj should be considered.
        std::cerr<<"maj is still smaller than true"<<std::endl;
        return false;
    }
    typedef boost::bernoulli_distribution<double> bernoulli_distrib;
    bernoulli_distrib fictitiousDistrib(1.0 - truer/majr);
    boost::variate_generator<rng_type&, bernoulli_distrib> fictitiousGenerator(rng, fictitiousDistrib);
    const bool isFictitious = fictitiousGenerator();
    return isFictitious;
}


// READ/WRITE/COPY.

// Writes the object to a binary stream.
void Process::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write the length of the name to the stream.
        unsigned int n = m_name.length();
        out.write((char*)&n, sizeof(n));

        // Write the name to the stream.
        out.write(m_name.c_str(), n);


        // Write reactant count.
        n = (unsigned int)m_reac.size();
        out.write((char*)&n, sizeof(n));

        // Write reactant stoichiometry.
        int m = 0;
        for (Sprog::StoichMap::const_iterator i=m_reac.begin(); i!=m_reac.end(); ++i) {
            // Write species ID.
            n = (unsigned int)i->first;
            out.write((char*)&n, sizeof(n));

            // Write stoichiometry.
            m = (int)i->second;
            out.write((char*)&m, sizeof(m));
        }

        // Write product count.
        n = (unsigned int)m_prod.size();
        out.write((char*)&n, sizeof(n));

        // Write product stoichiometry.
        for (Sprog::StoichMap::const_iterator i=m_prod.begin(); i!=m_prod.end(); ++i) {
            // Write species ID.
            n = (unsigned int)i->first;
            out.write((char*)&n, sizeof(n));

            // Write stoichiometry.
            m = (int)i->second;
            out.write((char*)&m, sizeof(m));
        }

        // Write scaling factor
        out.write(reinterpret_cast<const char*>(&m_a), sizeof(m_a));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Process::Serialize).");
    }
}

// Reads the object from a binary stream.
void Process::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    m_mech = &mech;
    m_prod.clear();
    m_reac.clear();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        int m = 0;
        unsigned int n = 0, id = 0;
        char *name = NULL;

        switch (version) {
            case 0:
                // Read the length of the process name.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read the process name.
                name = new char[n];
                in.read(name, n);
                m_name.assign(name, n);
                delete [] name;

                // Read reactant count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read reactant stoichiometry.
                for (unsigned int i=0; i!=n; ++i) {
                    // Read species ID.
                    in.read(reinterpret_cast<char*>(&id), sizeof(id));
                    // Read stoichiometry.
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));
                    // Create reactant.
                    m_reac[id] = m;
                }

                // Read product count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read product stoichiometry.
                for (unsigned int i=0; i!=n; ++i) {
                    // Read species ID.
                    in.read(reinterpret_cast<char*>(&id), sizeof(id));
                    // Read stoichiometry.
                    in.read(reinterpret_cast<char*>(&m), sizeof(m));
                    // Create product.
                    m_prod[id] = m;
                }

                // Read scaling factor
                double a;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                SetA(a);

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Process::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Process::Deserialize).");
    }
}


// PROTECTED HELPER FUNCTIONS.

/*!
 *@param[in]    gas     Gas phase mixture
 *
 *@return       Gas phase part of reaction rate
 */
double Process::chemRatePart(const EnvironmentInterface &gas) const
{
    double rate = 1.0;

    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        double conc = gas.SpeciesConcentration(i->first) ;
        for (int j=0; j!=i->second; ++j) {
            rate *= conc;
        }
    }

    return rate;
}

/*!
 * Adjusts the gas-phase composition using the reactants and
 * products defined for this process.
 *
 * Think carefully about whether calling this function is a
 * good idea:  The use of this method was found to be problematic
 * in "A predictor-corrector algorithm for the coupling of stiff
 * ODEs to a particle population balance", Celnik et al, Preprint 58.
 * The requirement that the gas phase mixture be implemented via a
 * SprogIdealGasWrapper will also restrict future applications of the
 * code.  In situations where updating the gas phase during stochastic
 * particle reactions is definitely required, this is the method
 * to use.
 *
 * @param[in,out]   sys     System in which the gas phase is changing
 * @param[in]       wt      Statistical weight of reaction
 * @param[in]       n       Repeat count of reaction
 *
 * @pre      The gas phase in sys must be of type SprogIdealGasWrapper
 *
 * @exception   std::runtime_error      Could not cast gas phase to SprogIdealGasWrapper
 */
void Process::adjustGas(Cell &sys, double wt, unsigned int n) const
{
    if(!sys.FixedChem()) {
        // This method requires write access to the gas phase, which is not
        // standard in sweep.  This means it cannot use the generic gas
        // phase interface
        SprogIdealGasWrapper *gasWrapper = dynamic_cast<SprogIdealGasWrapper*>(&sys.GasPhase());
        if(gasWrapper == NULL)
            throw std::runtime_error("Could not cast gas phase to SprogIdealGasWrapper in Process::adjustGas");

        // If excecution reaches here, the cast must have been successful
        Sprog::Thermo::IdealGas *gas = gasWrapper->Implementation();

        // Get the existing concentrations
        fvector newConcs;
        gas->GetConcs(newConcs);

        // Now adjust the concentrations
        Sprog::StoichMap::const_iterator i;
        double n_NAvol = wt * (double)n / (NA * sys.SampleVolume());
        for (i=m_reac.begin(); i!=m_reac.end(); ++i)
            newConcs[i->first] -= (double)(i->second) * n_NAvol;
        for (i=m_prod.begin(); i!=m_prod.end(); ++i)
            newConcs[i->first] +=(double)(i->second) * n_NAvol;

        // Set the new data
        gas->SetConcs(newConcs);
    }
}

// aab64 adjusts the particle-phase temperature using change in composition of the particle
/*
* Same warning as in adjustGas above applies
* @param[in, out]   sys      System in which the gas phase is changing
* @param[in]        wt       Statistical weight of reaction
* @param[in]        n        Repeat count of reaction
* @param[in]        adjustT  Update temperature for adiabatic case
*
* @pre      The gas phase in sys must be of type SprogIdealGasWrapper
*
* @exception   std::runtime_error      Could not cast gas phase to SprogIdealGasWrapper
*/
void Process::adjustParticleTemperature(Cell &sys, double wt, unsigned int n, bool adjustT, double dcomp, int processID) const {
	if (adjustT) {
		// Function will update the temperature oldTp, oldTg to temperature newTp, newTg for the particles and gas (K)
		double newTp, newTg, newRho;
		double oldTp = sys.GetBulkParticleTemperature();
		double oldTg = sys.GasPhase().Temperature();
		
		// This method requires write access to the gas phase, which is not
		// standard in sweep.  This means it cannot use the generic gas
		// phase interface
		SprogIdealGasWrapper *gasWrapper = dynamic_cast<SprogIdealGasWrapper*>(&sys.GasPhase());
		
		if (gasWrapper == NULL)
			throw std::runtime_error("Could not cast gas phase to SprogIdealGasWrapper in Process::adjustGas");

		// If excecution reaches here, the cast must have been successful
		Sprog::Thermo::IdealGas *gas = gasWrapper->Implementation();

		// Get particle-phase volume and area per m3 reactor
		double Vp = (sys.Particles().GetSum(iV)) / (sys.SampleVolume());
		double Sp = (sys.Particles().GetSum(iS)) / (sys.SampleVolume());

		// Overall heat transfer coefficient
		double h = 0.0;
		double rad = 0.0;
		int radSwitch = 1;   // 1 for PBB, 2 for RWR, 3 for RNR
		double radCoeff = 4; // Exponent of the temperature term in radiation

		if (processID == 0) {
			// From Tian et al. (2017) for titania and Michelsen et al. (2007)
			// Conduction
			double lambda = MeanFreePathAir(oldTg, gas->Pressure());                     // Mean free path, P.m.[K]^-1
			double ka = (5.83 * 1.0e7) * pow((oldTg / 273), 0.82);                       // Thermal conductivity of the gas, W.[m.K]^-1
			double alphaT = 1.0;                                                         // Thermal accomodation coefficient, -
			double gammaAir = 1.3;                                                       // Heat capacity ratio of air from Michelsen et al. (2007)
			double f = (9.0 * gammaAir - 5.0) / 4.0;                                     // Eucken correction
			double G = (8.0 * f) / (alphaT * (gammaAir + 1.0));                          // Factor in denominator of conductivity term
			double Dp = (sys.Particles().GetSum(iDcol)) / (sys.SampleVolume());          // Average particle diameter, m 
			h = (2.0 * ka) / (Dp + (G * lambda));                                        // Overall heat transfer coefficient, J.[m2.K.s]^-1

			// Radiation
			double Em = 0.251;                                                           // Function of the complex refractive index, m=2.51-1.7i => Em=Im[(m^2-1)/(m^2+2)]
			double hP = 6.62607004e-34;                                                  // Planck constant, m^2.kg.[s]^-1
			double cs = 299792458.0;                                                     // Speed of light, m.[s]^-1
			double kSB = 5.6704e-8;                                                      // Stefan-Boltzmann constant, W.[m2.K4]^-1
			if (radSwitch == 2 || radSwitch == 3) {
				rad = 1194.0 * pow(PI, 2) * Vp * Em * pow(KB, 5) / (pow(hP*cs, 3) * hP); // Build up radiation term for Rayleigh assumption
				radCoeff = 5;                                                            // T^5
			}
			else {
				rad = Sp * kSB;                                                          // Build up radiation term for perfect black body assumption
				radCoeff = 4;                                                            // T^4
			}
		}

		// Set enthalpy, heat capacity and density for both phases
		fvector Hs = gas->getMolarEnthalpy(oldTp);
		fvector Cs;
		gas->CalcCps(oldTp, Cs);
		double Cg = gas->BulkCp();          // bulk gp heat capacity at oldTg
		double Cp = Cs[28];                 // heat capacity of titania crystals at oldTp
		double Hp = Hs[28];                 // enthalpy of titania crystals at oldTp
		Hs = gas->getMolarEnthalpy(oldTg);  // enthalpy of gp species at oldTg
		double rhog = gas->Density();       // molar density of gp at oldTg
		double rhop = 53337;                // mol.[m3]^-1, molar density of titania (rutile)

		// Concentration change in system due to new particle(s)
		double n_NAvol = wt * (double)n / (NA * sys.SampleVolume());

		// Time step parameters
		double t0 = 0.0;
		double tf = sys.GetCurrentProcessTau();

		// Integration parameters a, b for gas (g) and particle phase (p)
		double ag = 0.0, ap = 0.0;
		double bg = 0.0, bp = 0.0;
		double gg = 0.0, gp = 0.0;

		// Distribute the formation energy between the phases
		double xg = 1.0; 
		double xp = 1.0 - xg; 

		if (tf != 0.0) {

			if (processID == 1 || processID == 2) { // inception or surface growth
				// Contributions of gp species
				Sprog::StoichMap::const_iterator i;
				for (i = m_reac.begin(); i != m_reac.end(); ++i) {
					ag -= (double)(i->second) * n_NAvol * Hs[i->first];
				}
				for (i = m_prod.begin(); i != m_prod.end(); ++i) {
					ag += (double)(i->second) * n_NAvol * Hs[i->first];
				}
				// Contribution of particle formation
				ag += (dcomp * n_NAvol * Hp * xg / (tf - t0)); // Divide by the time interval, multiply by arbitarily assigned gp component xg
				ag *= (-1.0);

				// Contributions of particles
				ap += (-1.0 * n_NAvol * dcomp * Hp * xp / (tf - t0));
			}
			else if (processID == 4) { // inflow
				ag += (-1.0 * n_NAvol * xg * (Hp - Hp) / (tf - t0)); // (temporary) Enthalpy at Tin - how to get Tin?
				ap += (-1.0 * n_NAvol * xp * (Hp - Hp) / (tf - t0)); // (temporary) Enthalpy at Tin - how to get Tin?
			}
			
			if (processID == 0) {
				double dt = tf - t0;

				// Calculate new temperatures using implicit difference method
				/*newTp = oldTp;
				newTp += dt * ap * oldTg / (1 + dt * ag);
				newTp *= 1.0 / (1 + (dt * ap) - ((dt * ap * dt * ag) / (1 + dt * ag)));

				newTg = oldTg;
				newTg += dt * ag * newTp;
				newTg *= 1.0 / (1 + dt * ag);*/

				// Calculate new temperatures using implicit difference method and adding nominal radiation from oldTp
				/*newTp = oldTp;
				newTp += (dt * ap / (1 + dt * ag)) * (oldTg + (dt * gg * rad));
				newTp += (dt * gp * rad);
				newTp *= 1.0 / (1 + (dt * ap) - ((dt * ap * dt * ag) / (1 + dt * ag)));

				newTg = oldTg;
				newTg += dt * ag * newTp;
				newTg += dt * gg * rad;
				newTg *= 1.0 / (1 + dt * ag);*/

				//rad = Sp * kSB * (pow(newTp, 4) - pow(newTg, 4));

				// Heat transfer pseudo-constants
				/*ag = (h * Sp) / (Cg * rhog);
				ap = (h * Sp) / (Cp * rhop * Vp);
				gg = rad / (Cg * rhog);
				gp = rad / (Cp * rhop * Vp);*/

				ag = (h * Sp) / (Cg);               // Divided by rhog. Technically also h is a function of tg
				ap = (h * Sp) / (Cp * rhop * Vp);   // Technically also h is a function of tg
				gg = rad / (Cg);                    // Divided by rhog
				gp = rad / (Cp * rhop * Vp);

				int newt_its = 1000;
				double tol = 1e-2;
				double c_a, c_b, c_c, c_d, c_e, c_f, c_g, c_h, c_i;
				double x1 = oldTp, x2 = oldTg, x3 = rhog;
				double dx1, dx2, dx3;
				double f1, f2, f3, det, cof11, cof12, cof21, cof22, x13, x14, x23, x24, x22;
				double j11, j12, j13, j21, j22, j23, j31, j32, j33;

				c_a = 1.0 + dt * ap;
				c_b = radCoeff * dt * gp;
				c_c = -1.0 * dt * ap;
				c_d = -1.0 * radCoeff * dt * gp;
				/*c_e = -1.0 * dt * ag;
				c_f = -1.0 * radCoeff * dt * gg;
				c_g = 1.0 + dt * ag;
				c_h = radCoeff * dt * gg;*/

				for (int i = 0; i < newt_its; i++) {

					c_e = -1.0 * dt * ag / x3;
					c_f = -1.0 * radCoeff * dt * gg / x3;
					c_g = 1.0 + (dt * ag / x3);
					c_h = radCoeff * dt * gg / x3;
					c_i = -1.0 * dt * gg / x3;

					x13 = pow(x1, radCoeff - 1);
					x14 = pow(x1, radCoeff);
					x22 = pow(x2, 2);
					if (radSwitch != 3){
						x23 = pow(x2, radCoeff - 1);
						x24 = pow(x2, radCoeff);
					}
					else {
						x23 = 0.0;
						x24 = 0.0;
					}

					f1 = (c_a * x1) + ((c_b / radCoeff) * x14) + (c_c * x2) + ((c_d / radCoeff) * x24) - oldTp;
					f2 = (c_g * x2) + ((c_h / radCoeff) * x24) + (c_e * x1) + ((c_f / radCoeff) * x14) - oldTg;
					f3 = (2.0 * x3) - (oldTg * (x3 / x2)) - rhog;
					
					/* Adding change in density makes this 3x3 (see below)
					det = ((c_a + (c_b * x13)) * (c_g + (c_h * x23))) - 
						  ((c_c + (c_d * x23)) * (c_e + (c_f * x13)));
					det = 1.0 / det;
					
					cof11 = c_g + (c_h * x23);
					cof12 = (-1.0 * c_c) + (-1.0 * c_d * x23);
					cof21 = (-1.0 * c_e) + (-1.0 * c_f * x13);
					cof22 = c_a + (c_b * x13);

					dx1 = det * ((cof11 * f1) + (cof12 * f2));
					dx2 = det * ((cof21 * f1) + (cof22 * f2));*/

					j11 = c_a + (c_b * x13);
					j12 = c_c + (c_d * x23);
					j13 = 0.0;
					j21 = c_e + (c_f * x13);
					j22 = c_g + (c_h * x23);
					j23 = (c_e / x3) * (x2 - x1) + (c_i / x3) * (x24 - x14); // Ignoring rhog: 0.0;
					j31 = 0.0;
					j32 = oldTg * (x3 / x22);
					j33 = 2.0 - (oldTg / x2);

					det = j11 * (j22 * j33 - j23 * j32) - j12 * (j21 * j33 - j23 * j31) + j13 * (j21 * j32 - j22 * j31);
					det = 1.0 / det;

					dx1 = det * (((j22 * j33 - j23 * j32) * f1) + ((j13 * j32 - j12 * j33) * f2) + ((j12 * j23 - j13 * j22) * f3));
					dx2 = det * (((j23 * j31 - j21 * j33) * f1) + ((j11 * j33 - j13 * j31) * f2) + ((j13 * j21 - j11 * j23) * f3));
					dx3 = det * (((j21 * j32 - j22 * j31) * f1) + ((j12 * j31 - j11 * j32) * f2) + ((j11 * j22 - j12 * j21) * f3));

					x1 = x1 - dx1;
					x2 = x2 - dx2;
					x3 = x3 - dx3;

					if (abs(dx1) < tol && abs(dx2) < tol && abs(dx3) < tol) {
						break;
					}
				}
				if (abs(dx1) > tol || abs(dx2) > tol || abs(dx3) > tol) {
					cout << "-------------------------------------\n";
					cout << "Max its reached without finding root!\n";
					cout << "Tp0=" << oldTp << ", Tg0=" << oldTg << ", rhog0=" << rhog << "\n";
					cout << "|dx1|=" << abs(dx1) << ", |dx2|=" << abs(dx2) << ", |dx3|=" << abs(dx3) << "\n";
					cout << "Tpf=" << x1 << ", Tgf=" << x2 << ", rhogf=" << x3 << "\n";
					cout << "-------------------------------------\n";
				}
				newTp = x1;
				newTg = x2;
				newRho = x3;
			}
			else {
				// Add denominator
				ag *= (1.0 / (Cg * rhog));
				ap *= (1.0 / (Cp * rhop * Vp));

				// Solve for new particle and gp temperatures
				newTp = oldTp + ap * (tf - t0);
				newTg = oldTg + ag * (tf - t0);

				//Solve for new gas density
				newRho = rhog;
				newRho *= 1.0 / (1.0 + (1.0 / newTg) * (newTg - oldTg));
			}
			
			// Update particle temperature and gas density
			sys.SetBulkParticleTemperature(newTp);
			gas->SetTemperature(newTg);
			gas->SetDensity(newRho);

			// Update density and sample volume using IG relationship
			//double oldRho = gas->Density();
			//double oldRhom = gas->MassDensity();
			//double newRho = oldRho - ((oldRho / oldTg) * (newTg - oldTg)); // constP
			//double newRho = rhog;
			//newRho *= 1.0 / (1.0 + (1.0/newTg) * (newTg - oldTg));
			//gas->SetDensity(newRho);
			//sys.AdjustSampleVolume(oldRhom / gas->MassDensity());
			//sys.AdjustSampleVolume(newRho * newTg / (oldRho * oldTg));
		}
	}
}

