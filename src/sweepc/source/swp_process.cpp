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

#include "swp_params.h"
#include "swp_process.h"
#include "swp_mechanism.h"
#include <stdexcept>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Process::Process(void)
: m_name(""), m_mech(NULL), m_viscosity_model(Sweep::iAir)
{
}

// Initialising constructor.
Process::Process(const Sweep::Mechanism &mech)
: m_name(""), m_mech(&mech), m_viscosity_model(Sweep::iAir)
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
        m_viscosity_model = rhs.m_viscosity_model;

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

//! Returns the viscosity model
Sweep::ViscosityModel Process::ViscosityModel() const
{
    return m_viscosity_model;
}

/*!
 * @brief           Sets the process' viscosity model
 * @param vmodel    Model type
 */
void Process::SetViscosityModel(Sweep::ViscosityModel vmodel)
{
    m_viscosity_model = vmodel;
}

/*!
 * @brief       Returns the viscosity.
 *
 * @param sys   System for which to calculate viscosity
 * @return      Viscosity (kg/ms)
 */
real Process::GetViscosity(const Cell &sys) const
{
    real mu(0.0);
    // Return the viscosity for air
    if (m_viscosity_model == Sweep::iAir) {
        real T = sys.GasPhase().Temperature();
        mu = Sweep::ViscosityAir(T);
    } else if (m_viscosity_model == Sweep::iChampanEnskog) {
        mu = sys.GasPhase().getViscosity();
    } else {
        throw runtime_error("Wrong viscosity identifier. "
                                "(Sweep, Process::GetViscosity).");
    }
    return mu;
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

// Adds a product given the species name.
void Process::AddProduct(const std::string &name, int mu)
{
    // Locate the species.
    unsigned int i;
    for (i=0; i!=m_mech->Species()->size(); ++i) {
        if (name.compare((*m_mech->Species())[i]->Name()) == 0) {
            // Set product.
            m_prod[i] = mu;
            return;
        }
    }
}

// Removes a product, given by name, from the reaction.
void Process::RemoveProduct(const std::string &name)
{
    // Locate the species.
    unsigned int i;
    for (i=0; i!=m_mech->Species()->size(); ++i) {
        if (name.compare((*m_mech->Species())[i]->Name()) == 0) {
            // Delete product.
            m_prod.erase(i);
            return;
        }
    }
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
bool Process::Fictitious(real majr, real truer, rng_type &rng)
{
    // ensure truer <= maj, otherwise it will crash the program.
    if (truer>majr)
    {
        // if the following warning msg show frequently, a better solution ensuring truer <= maj should be considered.
        std::cerr<<"maj is still smaller than true"<<std::endl;
        return false;
    }
    typedef boost::bernoulli_distribution<real> bernoulli_distrib;
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

// Calculates the gas-phase chemistry contribution to the rate
// expression.
real Process::chemRatePart(const fvector &fracs, real density) const
{
    real rate = 1.0;

    Sprog::StoichMap::const_iterator i;
    for (i=m_reac.begin(); i!=m_reac.end(); ++i) {
        real conc = density * fracs[i->first];
        for (int j=0; j!=i->second; ++j) {
            rate *= conc;
        }
    }

    return rate;
}

// Adjusts the gas-phase composition using the reactants and
// products defined for this process.
void Process::adjustGas(Cell &sys, real wt, unsigned int n) const
{
    if(!sys.FixedChem()) {
        //riap 7 Sep 2009
        // Only adjust the gas if the chemistry is not fixed.
        fvector dc(m_mech->Species()->size(), 0.0);
        Sprog::StoichMap::const_iterator i;
        real n_NAvol = wt * (real)n / (NA * sys.SampleVolume());
        for (i=m_reac.begin(); i!=m_reac.end(); ++i)
            dc[i->first] -= (real)(i->second) * n_NAvol;
        for (i=m_prod.begin(); i!=m_prod.end(); ++i)
            dc[i->first] +=(real)(i->second) * n_NAvol;
        sys.AdjustConcs(dc);
    }
}
