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

// Adjusts the gas-phase composition and temperature 
// using the change in composition of the particle
// Currently only implemented for titania!
/*
* Same warning as in adjustGas above applies
* @param[in, out]   sys        System in which the gas phase is changing
* @param[in]        wt         Statistical weight of reaction
* @param[in]        n          Repeat count of reaction
* @param[in]        dcomp      Composition change in particle
* @param[in]        processID  Identify process as reaction (1, 2) or flow (3, 4)
*
* @pre      The gas phase in sys must be of type SprogIdealGasWrapper
*
* @exception   std::runtime_error      Could not cast gas phase to SprogIdealGasWrapper
*/
void Process::adjustParticleTemperature(Cell &sys, double wt, unsigned int n, double dcomp, int processID) const 
{
	if (!sys.FixedChem())
	{
		// This method requires write access to the gas phase, which is not
		// standard in sweep.  This means it cannot use the generic gas
		// phase interface
		SprogIdealGasWrapper *gasWrapper = dynamic_cast<SprogIdealGasWrapper*>(&sys.GasPhase());

		if (gasWrapper == NULL)
			throw std::runtime_error("Could not cast gas phase to SprogIdealGasWrapper in Process::adjustGas");

		// If excecution reaches here, the cast must have been successful
		Sprog::Thermo::IdealGas *gas = gasWrapper->Implementation();

		// Function will update the temperature oldTg to temperature newTg
		double oldTg = sys.GasPhase().Temperature();
		double newTg = oldTg;

		// Thermal parameters 
		// Cg, Cp, Cp * rhop and Hs/Us are stored on sys object for
		// faster computation. They are updated less frequently.
		fvector Hs, Cs;
		sys.getEnthalpies(Hs);
		double Cg = sys.getBulkHeatCapacity();
		double Cp = sys.getParticleHeatCapacity();
		bool constv = sys.IsConstV();
		int pindex = m_mech->GetParticleSpeciesIndex();
		// Set enthalpy, heat capacity and density for both phases
		// Uncomment to update properties each time this function is called. 
		/*if (constv)
		{
			gas->Us(Hs);
			gas->CalcCvs(oldTg, Cs);
			Cg = gas->BulkCv();
			Cp = Cs[pindex]; 
		}
		else
		{
			gas->Hs(Hs); 
			gas->CalcCps(oldTg, Cs);
			Cg = gas->BulkCp();
			Cp = Cs[pindex];
		}*/
		double mw = sys.ParticleModel()->Components()[0]->MolWt(); 
		double rhop = (sys.Particles().GetSum(iM) + sys.Particles().GetTotalMass()) / (mw * sys.SampleVolume());
		double rhog = gas->Density();
		double Cprhop = Cp * rhop;
		double Cgrhog = Cg * rhog;
		double rhog_mass = gas->MassDensity();
		double P_R = gas->Pressure() / R;
		double Hp = Hs[pindex];

		// Get the existing concentrations
		fvector newConcs;
		gas->GetConcs(newConcs);

		// Concentration change in system due to new particle(s)
		double n_NAvol = wt * (double)n / (NA * sys.SampleVolume());

		// Compute the change in concentrations and enthalpy
		double deltaHr = 0.0;
		double gdot = 0.0;
		double change = 0.0;
		if (processID == 1 || processID == 2) {                // Inception or surface growth
			Sprog::StoichMap::const_iterator i;
			for (i = m_reac.begin(); i != m_reac.end(); ++i)
			{
				change = (double)(i->second) * n_NAvol;
				newConcs[i->first] -= change;
				deltaHr -= change * Hs[i->first];
				gdot -= change;
			}
			for (i = m_prod.begin(); i != m_prod.end(); ++i)
			{
				change = (double)(i->second) * n_NAvol;
				newConcs[i->first] += change;
				deltaHr += change * Hs[i->first];
				gdot += change;
			}
			deltaHr += (dcomp * n_NAvol * Hp);                 // Contribution of particle formation
		} else if (processID == 4) {                           // Inflow - need to get (f/tau and?) Hp_in
			double Hp_in = Hp;
			deltaHr += n_NAvol * (Hp_in - Hp);
		}

		// Start to adjust the temperature by subtracting deltaHr/(rho*C)
		deltaHr *= 1.0 / (Cgrhog + Cprhop);
		newTg -= deltaHr;
		
		if (constv)
		{
			// a) dCk/dt = gdot_k
			// b) drho/dt = dC/dt (constant volume)
			// c) dT/dt * rho * Cv = - deltaHr_total + deltaHflow - (drho/dt * R * T) 
		
			// Adjust the concentrations as in a) and b) above
			gas->SetConcs(newConcs);

			// Add denominator to gdot
			gdot *= 1.0 / (Cgrhog + Cprhop);
			
			// Finish adjusting the temperature as in c) above
			// newTg *= 1.0 / (1.0 + R * gdot);                // This way assumes T at tf
			newTg -= (R * oldTg * gdot);                       // This way assumes T at t0
			gas->SetTemperature(newTg);
		}
		else
		{
			// a) dT/dt * rho * Cv = - deltaHr_total + deltaHflow
			// b) gamma = (drho/dt / rho) + (dT/dt / T)
			// c) dCk/dt = gdot_k - gamma * Ck
			// d) rho = P / (R * T)

			// Adjust the temperature as in a) above
			// Taking rho at tf:
			//double a = Cprhop;
			//double b = (Cg * P_R) - (oldTg * Cprhop) + (deltaHr);
			//double c = -oldTg * Cg * P_R;
			//newTg = (-b + sqrt((b * b) - (4 * a * c))) / (2 * a); 
			// Taking rho at t0, temperature update is already complete.
			gas->SetTemperature(newTg);

			// Compute gas-phase expansion as in b) above
			double gamma = (gdot / rhog) + ((newTg - oldTg) / oldTg);	

			// Adjust the concentrations as in c) above
			fvector::const_iterator j;					
			for (int j = 0; j != newConcs.size(); ++j)
				newConcs[j] -= (gamma * newConcs[j]);
			gas->SetConcs(newConcs);
			
			// Account for gas phase expansion (via gamma)
			sys.AdjustSampleVolume(rhog_mass / gas->MassDensity());
		}
	}
}
