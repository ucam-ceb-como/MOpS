/*
  Author(s):      Matthew Celnik (msc37) and Markus Sander
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the Mechanism class declared in the
    swp_mechanism.h header file.

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

#include "swp_mechanism.h"
#include "swp_model_factory.h"
#include "swp_process_factory.h"
#include "swp_tempwriteXmer.h"
#include "swp_pah_inception.h"

#include "geometry1d.h"

#include <stdexcept>
#include <cassert>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>

#include "string_functions.h"
#include "swp_particle.h"
#include "swp_PAH_primary.h"


using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;
using namespace Strings;


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Mechanism::Mechanism(void)
: m_anydeferred(false), m_icoag(-1), m_termcount(0), m_processcount(0)
{
}

// Copy constructor.
Mechanism::Mechanism(const Mechanism &copy)
{
	*this = copy;
}

// Default destructor.
Mechanism::~Mechanism(void)
{
    releaseMem();
}


// OPERATOR OVERLOADS.

// Assignment operator.
Mechanism &Mechanism::operator=(const Mechanism &rhs)
{
    if (this != &rhs) {
        // Clear current mechanism from memory.
        releaseMem();

        // Copy base class.
        ParticleModel::operator=(rhs);

        // Easily copied data.
        m_anydeferred  = rhs.m_anydeferred;
        m_icoag        = rhs.m_icoag;
        m_termcount    = rhs.m_termcount;
        m_processcount = rhs.m_processcount;

        // Copy inceptions.
        for (IcnPtrVector::const_iterator i=rhs.m_inceptions.begin();
            i!=rhs.m_inceptions.end(); ++i) {
            m_inceptions.push_back((*i)->Clone());

            // Need to update the parent mechanism
            m_inceptions.back()->SetMechanism(*this);
        }

        // Copy particle processes.
        for (PartProcPtrVector::const_iterator i=rhs.m_processes.begin();
            i!=rhs.m_processes.end(); ++i) {
            m_processes.push_back((*i)->Clone());

            // Need to update the parent mechanism
            m_processes.back()->SetMechanism(*this);
        }

        // Copy coagulation processes.
        for (CoagPtrVector::const_iterator i=rhs.m_coags.begin();
            i!=rhs.m_coags.end(); ++i) {
            m_coags.push_back((*i)->Clone());

            // Need to update the parent mechanism
            m_coags.back()->SetMechanism(*this);
        }


        // Copy process counters.
        m_proccount.assign(rhs.m_proccount.begin(), rhs.m_proccount.end());
        m_fictcount.assign(rhs.m_fictcount.begin(), rhs.m_fictcount.end());
    }
    return *this;
}


// INCEPTIONS.

// Returns the vector of inceptions.
const IcnPtrVector &Mechanism::Inceptions(void) const
{
    return m_inceptions;
}

// Returns the inception with the given index.
const Inception *const Mechanism::Inceptions(unsigned int i) const
{
    if (i < m_inceptions.size()) {
        return m_inceptions[i];
    } else {
        return NULL;
    }
}



// Adds an inception to the mechanism.
void Mechanism::AddInception(Inception &icn)
{

    m_termcount += icn.TermCount();
    ++m_processcount;
    m_proccount.resize(m_termcount, 0);
    m_fictcount.resize(m_termcount, 0);

    // Set the inception to belong to this mechanism.
    icn.SetMechanism(*this);

    // Add the inception to the mechanism last as a step towards exception safety.
    m_inceptions.push_back(&icn);
}


// PARTICLE PROCESSES.

// Returns the vector of particle processes.
const PartProcPtrVector &Mechanism::Processes(void) const
{
    return m_processes;
}

// Returns the process with the given index.
const ParticleProcess *const Mechanism::Processes(unsigned int i) const
{
    if (i < m_processes.size()) {
        return m_processes[i];
    } else {
        return NULL;
    }
}

// Adds a process to the mechanism.
void Mechanism::AddProcess(ParticleProcess &p)
{
    // Add the process to the mechanism.
    m_processes.push_back(&p);
    m_termcount += p.TermCount();
    ++m_processcount;
    m_proccount.resize(m_termcount, 0);
    m_fictcount.resize(m_termcount, 0);

    // Check for any deferred.
    m_anydeferred = m_anydeferred || p.IsDeferred();

    // Set the process to belong to this mechanism.
    p.SetMechanism(*this);
}

// COAGULATIONS.

/*!
 * @param[in,out]   coag    New coagulation process
 *
 * Ownership of the process will be taken by the mechanism.  The
 * process must be heap allocated so that delete can be called on
 * it.
 */
void Mechanism::AddCoagulation(Coagulation& coag)
{

    m_coags.push_back(&coag);
    m_termcount += coag.TermCount();
    ++m_processcount;

    m_termcount += coag.TermCount();
    m_proccount.resize(m_termcount, 0);
    m_fictcount.resize(m_termcount, 0);

    // Set the coagulation to belong to this mechanism.
    coag.SetMechanism(*this);
}

const CoagPtrVector &Mechanism::Coagulations(void) const
{
    return m_coags;
}

// PROCESS INFORMATION.

// Returns the number of processes (including
// inceptions) in the mechanism.
unsigned int Mechanism::ProcessCount(void) const
{
    return m_processcount;
}

// Returns the number of terms in all process rate expressions.
unsigned int Mechanism::TermCount(void) const {return m_termcount;}

// Returns true if the mechanism contains deferred (LPDA)
// processes otherwise false.
bool Mechanism::AnyDeferred(void) const {return m_anydeferred;}

// Checks all processes to see if any are deferred.
void Mechanism::CheckDeferred(void) const
{
    // Loop though all processes checking if any are deferred.
    m_anydeferred = false;
    PartProcPtrVector::const_iterator i;
    for (i=m_processes.begin(); i!=m_processes.end(); ++i) {
        if ((*i)->IsDeferred()) {
            // Set anydeferred flag is true.
            m_anydeferred = true;
            return;
        }
    }
}


// Returns a vector containing the names of all processes.
void Mechanism::GetProcessNames(std::vector<std::string> &names,
                                unsigned int start) const
{
    // Resize the output vector to hold all process names from
    // start index onwards.
    unsigned int N = start + m_processcount; //m_inceptions.size() + m_processes.size();
//    if (m_coag) ++N;
    names.resize(N, "");

    // Get iterator to first insertion point.
    vector<string>::iterator i = names.begin()+start;

    // Add inception names.
    for (unsigned int j=0; j!=m_inceptions.size(); ++j) {
        *i = m_inceptions[j]->Name(); ++i;
    }

    // Add particle process names.
    for (unsigned int j=0; j!=m_processes.size(); ++j) {
        *i = m_processes[j]->Name(); ++i;
    }

    // Add coagulation name.
    for (CoagPtrVector::const_iterator it = m_coags.begin();
         it != m_coags.end(); ++it) {
        *i++ = (*it)->Name();
    }
}


// RATE CALCULATION.

// Get total rates of all processes.  Returns the sum of
// all rates.
double Mechanism::CalcRates(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom, fvector &rates, bool scale) const
{
    // Ensure rates vector is the correct length, then set to zero.
    rates.resize(m_processcount+sys.InflowCount()+sys.OutflowCount(), 0.0);
    fill(rates.begin(), rates.end(), 0.0);

    double sum = 0.0;

    // Get rates of inception processes.
    sum += Inception::CalcRates(t, sys, local_geom, m_inceptions, rates);

    // Query other processes for their rates.
    sum += ParticleProcess::CalcRates(t, sys, local_geom, m_processes, rates, m_inceptions.size());

    // Get coagulation rate.
    sum += Coagulation::CalcRates(t, sys, local_geom, m_coags, rates, m_inceptions.size() + m_processes.size());

    // Get birth rates from the Cell.
    fvector::iterator i = rates.begin() + m_inceptions.size() + m_processes.size() + m_coags.size();
    const BirthPtrVector &inf = sys.Inflows();
    for (BirthPtrVector::const_iterator j=inf.begin(); j!=inf.end(); ++j) {
        *i = (*j)->Rate(t, sys, local_geom);
        sum += *i;
        ++i;
    }

    // Get death rates from the Cell.
    const DeathPtrVector &outf = sys.Outflows();
    for (DeathPtrVector::const_iterator j=outf.begin(); j!=outf.end(); ++j) {
        *i = (*j)->Rate(t, sys, local_geom);
        sum += *i;
        ++i;
    }

    if (!scale) {
        // Need to return the rates to per unit vol.
        double invvol = 1.0 / sys.SampleVolume();
        for(i=rates.begin(); i!=rates.end(); ++i) {
            *i *= invvol;
        }
        sum *= invvol;
    }

    return sum;
}

/*!
 * @brief               Calculates the number of jump events for each process
 *
 * Calculates the absolute number of jump events for each inception, particle
 * process and coagulation event. Returns the sum of the jump events. Ignores
 * transport processes.
 *
 * @param t             Time
 * @param sys           Particle population
 * @param local_geom    Pointer to local geometry
 * @param jumps         Vector containing the number of jumps
 * @return              Sum of jump events
 */
double Mechanism::CalcJumps(double t, const Cell &sys, const Geometry::LocalGeometry1d &local_geom, fvector &jumps) const
{
    // Ensure jumps vector is the correct length, then set to zero.
    jumps.resize(m_processcount, 0.0);
    fill(jumps.begin(), jumps.end(), 0.0);

    // Iterator for filling jumps vector
    fvector::iterator iterm = jumps.begin();

    double sum = 0.0;

    // Get number of inception jumps
    for (unsigned int j=0; j!=m_inceptions.size(); ++j) {
        (*iterm++) = m_proccount[j];
        sum += m_proccount[j];
    }

    // Get number of particle process jumps
    for (unsigned int j=0; j!=m_processes.size(); ++j) {
        (*iterm++) = m_proccount[j+m_inceptions.size()];
        sum += m_proccount[j+m_inceptions.size()];
    }

    // Get number of coagulation jumps.
    unsigned int coagterms(0);       // Number of terms already used
    for (unsigned int j=0; j!=m_coags.size(); ++j) {
        unsigned int coagsum(0);     // Sum of double and fictitious jumps
        // Sum up all terms of this process
        for (unsigned int k=0; k!=m_coags[j]->TermCount(); ++k) {
            coagsum += m_proccount[k+m_inceptions.size()+m_processes.size()+coagterms];
            coagsum += m_fictcount[k+m_inceptions.size()+m_processes.size()+coagterms];
        }
        (*iterm++) = coagsum;
        sum += coagsum;
        coagterms += m_coags[j]->TermCount();
    }

    return sum;
}

/*!
 * @brief       Resets the counters which track the number of jumps
 *
 * This function is to only be called at the end of a run, to ensure
 * accurate capturing of information from CalcJumps
 */
void Mechanism::ResetJumpCount() const {
    // Do for number of double jumps
    fill(m_proccount.begin(), m_proccount.end(), 0.0);
    // Do for number of fictitious jumps
    fill(m_fictcount.begin(), m_fictcount.end(), 0.0);
}


// Get rates of all processes separated into different
// terms.  Rate terms are useful for subsequent particle
// selection by different properties for the same process.
// In particular this is used for the condensation and
// coagulation processes.  Returns the sum of all rates.
double Mechanism::CalcRateTerms(double t, const Cell &sys, const Geometry::LocalGeometry1d& local_geom, fvector &terms) const
{
    // Ensure rates vector is the correct length.
    terms.resize(m_termcount+sys.InflowCount()+sys.OutflowCount(), 0.0);
    fvector::iterator iterm = terms.begin();

    double sum = 0.0;

    // Get rates of inception processes.
    IcnPtrVector::const_iterator ii;
    for (ii=m_inceptions.begin(); ii!=m_inceptions.end(); ++ii) {
        sum += (*ii)->RateTerms(t, sys, local_geom, iterm);
    }

    // Query other processes for their rates.
    if (sys.ParticleCount() > 0) {
        for(PartProcPtrVector::const_iterator i=m_processes.begin();
            (i!=m_processes.end()) && (iterm!=terms.end()); ++i) {
            sum += (*i)->RateTerms(t, sys, local_geom, iterm);
        }
    } else {
        // Fill vector with zeros.
        for(PartProcPtrVector::const_iterator i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=terms.end()); ++i) {
            fill(iterm, iterm+(*i)->TermCount(), 0.0);
            iterm += (*i)->TermCount();
        }
    }

    // Coagulation
    sum += Coagulation::CalcRateTerms(t, sys, local_geom, m_coags, iterm);

    // Birth processes
    const BirthPtrVector &inf = sys.Inflows();
    for (BirthPtrVector::const_iterator j=inf.begin(); j!=inf.end(); ++j) {
        sum += (*j)->RateTerms(t, sys, local_geom, iterm);
    }

    // Death processes
    const DeathPtrVector &outf = sys.Outflows();
    for (DeathPtrVector::const_iterator j=outf.begin(); j!=outf.end(); ++j) {
        sum += (*j)->RateTerms(t, sys, local_geom, iterm);
    }

    return sum;
}

// Get total rates of non-deferred processes.  Returns the sum
// of all rates.
double Mechanism::CalcJumpRateTerms(double t, const Cell &sys, const Geometry::LocalGeometry1d& local_geom, fvector &terms) const
{
    // This routine only calculates the rates of those processes which are
    // not deferred.  The rate terms of deferred processes are returned
    // as zero.

    // Ensure vector is the correct length, then set to zero.
    terms.resize(m_termcount+sys.InflowCount()+sys.OutflowCount(), 0.0);
    fvector::iterator iterm = terms.begin();

    double sum = 0.0;

    // Get rates of inception processes.
    IcnPtrVector::const_iterator ii;
    for (ii=m_inceptions.begin(); ii!=m_inceptions.end(); ++ii) {
        sum += (*ii)->RateTerms(t, sys, local_geom, iterm);
    }

    // Query other processes for their rates.
    if (sys.ParticleCount() > 0) {
        for(PartProcPtrVector::const_iterator i=m_processes.begin();
            (i!=m_processes.end()) && (iterm!=terms.end()); ++i) {
            if (!(*i)->IsDeferred()) {
                // Calculate rate if not deferred.
                sum += (*i)->RateTerms(t, sys, local_geom, iterm);
            } else {
                // If process is deferred, then set rate to zero.
                for (unsigned int j=0; j!=(*i)->TermCount(); ++j) {*(iterm++)=0.0;}
            }
        }
    } else {
        // Fill vector with zeros.
        for(PartProcPtrVector::const_iterator i=m_processes.begin();
            (i!=m_processes.end()) && (iterm!=terms.end()); ++i) {
            fill(iterm, iterm+(*i)->TermCount(), 0.0);
            iterm += (*i)->TermCount();
        }
    }

    // Get coagulation rate.
    sum += Coagulation::CalcRateTerms(t, sys, local_geom, m_coags, iterm);

    // Get birth rates from the Cell.
    const BirthPtrVector &inf = sys.Inflows();
    for (BirthPtrVector::const_iterator j=inf.begin(); j!=inf.end(); ++j) {
        sum += (*j)->RateTerms(t, sys, local_geom, iterm);
    }

    const DeathPtrVector &outf = sys.Outflows();
    for (DeathPtrVector::const_iterator j=outf.begin(); j!=outf.end(); ++j) {
        sum += (*j)->RateTerms(t, sys, local_geom, iterm);
    }

    return sum;
}


/*!
 * LPDA allows some processes to be deferred and removed from the main simulation
 * loop in a kind of splitting.  This method calculates the combined rate of all
 * the processes that are deferred in this way.  The result is the expected number
 * of events per second in sys, the second argument.  The result will therefore
 * scale linearly with the sample volume of sys, if particle concentrations are
 * held constant.
 *
 *@param[in]    t           Time at which to calculate the rates
 *@param[in]    sys         System for which rates are to be calculated
 *@param[in]    local_geom  Information on location of surrounding cells
 *@param[out]   terms       Vector to fill with the terms contributing to the summed rate
 *
 *@return   The total rate of all deferred processes
 */
double Mechanism::CalcDeferredRateTerms(double t, const Cell &sys, const Geometry::LocalGeometry1d& local_geom, fvector &terms) const
{
    // This routine only calculates the rates of those processes which are
    // deferred.  The rate terms of non-deferred processes are returned
    // as zero.

    // Ensure vector is the correct length and full of zeros.
    fill(terms.begin(), terms.end(), 0.0);
    terms.resize(m_termcount+sys.InflowCount()+sys.OutflowCount(), 0.0);
    fvector::iterator iterm = terms.begin();

    double sum = 0.0;

    // Query other processes for their rates.
    if (sys.ParticleCount() > 0) {
        for(PartProcPtrVector::const_iterator i=m_processes.begin();
            (i!=m_processes.end()) && (iterm!=terms.end()); ++i) {
            if ((*i)->IsDeferred()) {
                // Calculate rate if not deferred.
                sum += (*i)->RateTerms(t, sys, local_geom, iterm);
            }
        }
    }

    return sum;
}

/*!
 * Calculates the rates-of-change of the chemical species fractions,
 * gas-phase temperature and density due to particle processes.
 *
 *@param[in]    t      Time at which rates are to be calculated
 *@param[in]    sys    System for which rates are to be calculated
 *@param[in]    local_geom  Position information
 *@param[out]   rates  Vector of time rates of change of mole fractions, temperature and (?number) density
 *
 *@post  rates.size() == m_species.size() + 2
 * There is no precondition on rates.size().
 */
void Mechanism::CalcGasChangeRates(double t, const Cell &sys,
                                   const Geometry::LocalGeometry1d &local_geom,
                                   fvector &rates) const
{
    // Resize vector to hold all species and set all rates to zero.
    rates.resize(m_species->size()+2, 0.0);
    fill(rates.begin(), rates.end(), 0.0);

    // Rate of change of total concentration
    double idrho(0.0);

    // Precalculate parameters.
    double invVolNA = 1.0 / (sys.SampleVolume() * NA);

    // Inceptions and surface processes can affect the gas-phase chemistry
    // at the moment.

    // Loop over the contributions of all inception processes.
    for (IcnPtrVector::const_iterator i=m_inceptions.begin();
         i!=m_inceptions.end(); ++i) {

        // Calculate the inception rate.
        double rate = (*i)->Rate(t, sys, local_geom);

        // Loop over all reactants, subtracting their contributions.
        for (Sprog::StoichMap::const_iterator j=(*i)->Reactants().begin();
             j!=(*i)->Reactants().end(); ++j) {
            double dc = rate * (double)j->second * invVolNA;
            rates[j->first] -= dc;
            idrho -= dc;
        }

        // Loop over all products, adding their contributions.
        for (Sprog::StoichMap::const_iterator j=(*i)->Products().begin();
             j!=(*i)->Products().end(); ++j) {
            double dc = rate * (double)j->second * invVolNA;
            rates[j->first] += dc;
            idrho += dc;
        }
    }

    // Loop over the contributions of all other processes (except coagulation and transport).
    for (PartProcPtrVector::const_iterator i=m_processes.begin();
         i!=m_processes.end(); ++i) {

        // Calculate the process rate.
        double rate = (*i)->Rate(t, sys, local_geom);

        // Loop over all reactants, subtracting their contributions.
        for (Sprog::StoichMap::const_iterator j=(*i)->Reactants().begin();
             j!=(*i)->Reactants().end(); ++j) {
            double dc = rate * (double)j->second * invVolNA;
            rates[j->first] -= dc;
            idrho -= dc;
        }

        // Loop over all products, adding their contributions.
        for (Sprog::StoichMap::const_iterator j=(*i)->Products().begin();
             j!=(*i)->Products().end(); ++j) {
            double dc = rate * (double)j->second * invVolNA;
            rates[j->first] += dc;
            idrho += dc;
        }
    }

    // Now convert to changes in mole fractions.
    double invrho = 1.0 / sys.GasPhase().MolarDensity();
    for (unsigned int k=0; k!=m_species->size(); ++k) {
        // Quotient rule for dXk/dt = d(Ck/CT)/dt
        rates[k] = (invrho * rates[k]) - (invrho * invrho * sys.GasPhase().SpeciesConcentration(k) * idrho);
    }
}


// PERFORMING THE PROCESSES.

/*!
 * Performs the Process specified.  Process index could be
 * an inception, particle process or a coagulation event.
 *
 * \param[in]       i           Index of process to perform
 * \param[in]       t           Time at which event is to take place
 * \param[in,out]   sys         System in which event is to take place
 * \param[in]       local_geom  Information on surrounding cells for use with transport processes
 * \param[in,out]   rng         Random number generator
 *
 * The support for transport processes may well no longer be needed, in that it is
 * rarely efficient to simulate such phenomena with stochastic jumps.
 */
void Mechanism::DoProcess(unsigned int i, double t, Cell &sys,
                          const Geometry::LocalGeometry1d& local_geom,
                          rng_type &rng) const
{
    // Test for now
    assert(sys.ParticleModel() != NULL);

    // Work out to which process this term belongs.
    int j = i - m_inceptions.size();

    if (j < 0) {
        // This is an inception process.
        m_inceptions[i]->Perform(t, sys, local_geom, 0, rng);
        m_proccount[i] += 1;
    } else {
        // This is another process.
        for(PartProcPtrVector::const_iterator ip=m_processes.begin(); ip!=m_processes.end(); ++ip) {
            if (j < (int)(*ip)->TermCount()) {
                // Do the process.
                if ((*ip)->Perform(t, sys, local_geom, j, rng) == 0) {
                    m_proccount[i] += 1;
                } else {
                    m_fictcount[i] += 1;
                }
                return;
            } else {
                j -= (*ip)->TermCount();
            }
        }

        // We are here because the process was neither an inception
        // nor a single particle process.  It is therefore either a
        // coagulation or a birth/death process.
        for(CoagPtrVector::const_iterator it = m_coags.begin(); it != m_coags.end(); ++it) {
            // Check if coagulation process.
            if (j < static_cast<int>((*it)->TermCount())) {
                // This is the coagulation process.
                if ((*it)->Perform(t, sys, local_geom, j, rng) == 0) {
                    m_proccount[i] += 1;
                } else {
                    m_fictcount[i] += 1;
                }
                return;
            } else {
                // This must be the birth/death process.
                j -= (*it)->TermCount();
            }
        }

        if ((j < (int)sys.InflowCount()) && (j>=0)) {
            // An inflow process.
            sys.Inflows(j)->Perform(t, sys, local_geom, 0, rng);
            return;
        } else {
            // Hopefully a death process then!
            j -= sys.InflowCount();;
        }

        if ((j < (int)sys.OutflowCount()) && (j>=0)) {
            // An outflow process.
            sys.Outflows(j)->Perform(t, sys, local_geom, 0, rng);
        } else {
            throw std::runtime_error("Unknown index of process, couldn't Perform."
                    " (Sweep, Mechanism::DoProcess)");
        }

    }
}


/*!
 * The equivalent of the DoProcess function, but for particle transport
 * due to reactor flow and over time dt. The processes are done in a
 * pseudorandom order, as preliminary investigations indicated that a
 * set order (e.g. inflow first, then outflow) led to a statistically
 * significant difference in moments.
 *
 * @param t             Current system time
 * @param dt            Time over which to do process
 * @param sys           The system to perform the process on
 * @param local_geom    Geometry information
 * @param rng           Random number generator
 */
void Mechanism::DoParticleFlow(
        double t,
        double dt,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        rng_type &rng) const {
    Processes::ProcessPtrVector flows;
    std::vector<double> rates;

    // Get the processes and their 'rates' first
    for (Processes::BirthPtrVector::const_iterator it = sys.Inflows().begin();
            it != sys.Inflows().end(); ++it) {
        flows.push_back(*it);
        rates.push_back((*it)->A());
    }
    for (Processes::DeathPtrVector::const_iterator it = sys.Outflows().begin();
            it != sys.Outflows().end(); ++it) {
        flows.push_back(*it);
        rates.push_back((*it)->A());
    }

    // Using a discrete distribution weighted by the relative rates of each
    // process, progressively PerformDT each one.
    unsigned int i(0), j(0), nprocs(flows.size());
    while (i != nprocs) {
        // Get the vector index of the process to perform.
        boost::random::discrete_distribution<> dist(rates);
        j = dist(rng);

        // Do the process
        flows.at(j)->PerformDT(t, dt, sys, local_geom, rng);

        // Remove the process and its rates from the distribution
        rates.erase(rates.begin() + j);
        flows.erase(flows.begin() + j);

        // Increment the iterator
        i++;
    }

    // Now check if the death process needs adapting.
    for (Processes::DeathPtrVector::const_iterator it = sys.Outflows().begin();
            it != sys.Outflows().end(); ++it) {
        (*it)->Adapt(sys);
    }

}

/*!
 * Transfer InceptedPAH molecule (i.e. A1, A2 or A4), from the gas phase to the particle ensemble and vice versa used for PAH-PP model only
 * 
 * @param[in]       i           the number of pyrene supposed in the emsemble
 * @param[in]       t           Time at which event is to take place
 * @param[in,out]   sys         System in which event is to take place
 * @param[in,out]   rng         Random number generator
 *
 */
void Mechanism::MassTransfer(int i, double t, Cell &sys, rng_type &rng) const
{
        // Test for now
        assert(sys.ParticleModel() != NULL);
        // This is an inception process.
        const Sweep::Processes::PAHInception *m_pahinception = NULL;
        m_pahinception = dynamic_cast<const Sweep::Processes::PAHInception*>(m_inceptions[0]);
        m_pahinception->AddInceptedPAH(i, t, sys, rng);
}

// LINEAR PROCESS DEFERMENT ALGORITHM.

/*!
 * Performs linear process updates on all particles in a system.
 *
 *@param[in,out]    sys         System containing particles to update
 *@param[in]        t           Time upto which particles to be updated
 *@param[in,out]    rng         Random number generator
 */
void Mechanism::LPDA(double t, Cell &sys, rng_type &rng) const
{
    // Check that there are particles to update and that there are
    // deferred processes to perform.
    if ((sys.ParticleCount() > 0) &&
        (m_anydeferred ||
                (AggModel() == AggModels::PAH_KMC_ID) ||
                (AggModel() == AggModels::BinTree_ID) ||
                (AggModel() == AggModels::BinTreeSilica_ID) ||
                (AggModel() == AggModels::SurfVolSilica_ID) ||
                (AggModel() == AggModels::SurfVol_ID))) {
        // Stop ensemble from doubling while updating particles.
        sys.Particles().FreezeDoubling();

        // Perform deferred processes on all particles individually.
        Ensemble::iterator i;
        for (i=sys.Particles().begin(); i!=sys.Particles().end(); ++i) {
            UpdateParticle(*(*i), sys, t, rng);
        }

        // Now remove any invalid particles and update the ensemble.
        sys.Particles().RemoveInvalids();

        // Start particle doubling again.  This will also double the ensemble
        // if too many particles have been removed.
        sys.Particles().UnfreezeDoubling();
    }
}


/*!
 * Performs linear process updates on a particle in the given system.
 *
 *@param[in,out]    sp          Particle to update
 *@param[in,out]    sys         System containing particle to update
 *@param[in]        t           Time upto which particle to be updated
 *@param[in,out]    rng         Random number generator
 */
void Mechanism::UpdateParticle(Particle &sp, Cell &sys, double t, rng_type &rng) const
{
    // Deal with the growth of the PAHs
    if (AggModel() == AggModels::PAH_KMC_ID)
    {
        // Calculate delta-t and update particle time.
        double dt;
        dt = t - sp.LastUpdateTime();
        sp.SetTime(t);

        // If the agg model is PAH_KMC_ID then all the primary
        // particles must be PAHPrimary.
        AggModels::PAHPrimary *pah =
                dynamic_cast<AggModels::PAHPrimary*>(sp.Primary());

        // Update individual PAHs within this particle by using KMC code
        // sys has been inserted as an argument, since we would like use Update() Fuction to call KMC code
        pah->UpdatePAHs(t, *this, sys, rng);

        pah->UpdateCache();
        pah->CheckRounding();
        if (sp.IsValid()) {
            sp.UpdateCache();

            // Sinter the particles for the soot model (as no deferred process)
            if (m_sint_model.IsEnabled()) {
                pah->Sinter(dt, sys, m_sint_model, rng, sp.getStatisticalWeight());
            }
            sp.UpdateCache();
        }
    }

    // Sinter the particles if no deferred processes
    if (m_sint_model.IsEnabled() && !m_anydeferred
            && AggModel() != AggModels::PAH_KMC_ID) {

        // Calculate delta-t and update particle time.
        double dt;
        dt = t - sp.LastUpdateTime();
        sp.SetTime(t);

        sp.Sinter(dt, sys, m_sint_model, rng, sp.getStatisticalWeight());

        // Check particle is valid and recalculate cache.
        if (sp.IsValid()) {
            sp.UpdateCache();
        }
    }

    // If there are no deferred processes then stop right now.
    if (m_anydeferred) {
        PartProcPtrVector::const_iterator i;
        double rate, dt;

        while ((sp.LastUpdateTime() < t) && sp.IsValid()) {
            // Calculate delta-t and update particle time.
            dt = t - sp.LastUpdateTime();
            sp.SetTime(t);

            // Loop through all processes, performing those
            // which are deferred.
            for (i=m_processes.begin(); i!=m_processes.end(); ++i) {
                if ((*i)->IsDeferred()) {
                    // Get the process rate x the time interval.
                    rate = (*i)->Rate(t, sys, sp) * dt;

                    // Use a Poission deviate to calculate number of
                    // times to perform the process.  If the rate is
                    // 0 then the count is guaranteed to be 0
                    if(rate > 0) {
                         boost::random::poisson_distribution<unsigned, double> repeatDistrib(rate);
                         unsigned num = repeatDistrib(rng);
                         if (num > 0) {
                             // Do the process to the particle.
                             (*i)->Perform(t, sys, sp, rng, num);
                         }
                    }
                }
            }

            // Perform sintering update.
            if (m_sint_model.IsEnabled()) {
			    sp.Sinter(dt, sys, m_sint_model, rng, sp.getStatisticalWeight());
            }
        }

        // Check that the particle is still valid, only calculate
        // cache if it is.
        if (sp.IsValid())
            sp.UpdateCache();
    }
}

void Mechanism::Mass_pah(Ensemble &m_ensemble) const
{
    if (AggModel() == AggModels::PAH_KMC_ID) {
        std::vector<fvector> pah_vector;  // used for storing primary particle with specified mass

        for (size_t i=0;i!=m_ensemble.Capacity();++i){
            if (m_ensemble.At(i)!=NULL)
            {
                const Sweep::AggModels::PAHPrimary *rhsparticle = NULL;
                rhsparticle = dynamic_cast<const AggModels::PAHPrimary*>(m_ensemble.At(i)->Primary());
                rhsparticle->FindXmer(pah_vector,20);//2nd arguement is target_c, it can not be less than 10, otherwise it will be meaningless. 
            }
        }
        // create csv file for target primary particles
        writeParimary(pah_vector);
     }
}

// READ/WRITE/COPY.

// Creates a copy of the mechanism.
Mechanism *const Mechanism::Clone(void) const
{
    return new Mechanism(*this);
}

// Writes the object to a binary stream.
void Mechanism::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;

    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write particle model base class.
        ParticleModel::Serialize(out);

        // Write if any processes are deferred.
        if (m_anydeferred) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write number of inceptions.
        unsigned int n = (unsigned int)m_inceptions.size();
        out.write((char*)&n, sizeof(n));

        // Write inceptions.
        for (IcnPtrVector::const_iterator i=m_inceptions.begin();
             i!=m_inceptions.end(); ++i) {
            ProcessFactory::Write(*(*i), out);
        }

        // Write number of particle processes.
        n = (unsigned int)m_processes.size();
        out.write((char*)&n, sizeof(n));

        // Write particle processes.
        for (PartProcPtrVector::const_iterator i=m_processes.begin();
             i!=m_processes.end(); ++i) {
            ProcessFactory::Write(*(*i), out);
        }

        // Coagulation
        n = m_coags.size();
        out.write(reinterpret_cast<const char*>(&n), sizeof(n));
        for (CoagPtrVector::const_iterator i = m_coags.begin();
             i != m_coags.end(); ++i) {
            ProcessFactory::Write(*(*i), out);
        }

        // Write index of first coag process.
        int m = (int)m_icoag;
        out.write((char*)&m, sizeof(m));

        // Write term count.
        n = (unsigned int)m_termcount;
        out.write((char*)&n, sizeof(n));

        // Write process count.
        n = (unsigned int)m_processcount;
        out.write((char*)&n, sizeof(n));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, Mechanism::Serialize).");
    }
}

// Reads the object from a binary stream.
void Mechanism::Deserialize(std::istream &in)
{
    releaseMem();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        int m=0;
        unsigned int n=0;

        switch (version) {
            case 0:
                // Read ParticleModel base class.
                ParticleModel::Deserialize(in);

                // Read if any processes are deferred.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_anydeferred = (n==1);

                // Read number of inceptions.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read inceptions.
                for (unsigned int i=0; i!=n; ++i) {
                    Inception *icn = ProcessFactory::ReadInception(in, *this);
                    icn->SetMechanism(*this);
                    m_inceptions.push_back(icn);
                }

                // Read number of particle processes.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read particle processes.
                for (unsigned int i=0; i!=n; ++i) {
                    ParticleProcess *p = ProcessFactory::ReadPartProcess(in, *this);
                    p->SetMechanism(*this);
                    m_processes.push_back(p);
                }

                // Read coagulation process.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                for(unsigned int i = 0; i != n; ++i) {
                    Coagulation *pCoag = ProcessFactory::ReadCoag(in, *this);
                    pCoag->SetMechanism(*this);
                    m_coags.push_back(pCoag);
                }

                // Read index of first coag process.
                in.read(reinterpret_cast<char*>(&m), sizeof(m));
                m_icoag = m;

                // Read term count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_termcount = n;

                // Read process count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_processcount = n;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, Mechanism::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, Mechanism::Deserialize).");
    }
}


// MEMORY MANAGEMENT.

// Clears the current mechanism from memory.
void Mechanism::releaseMem(void)
{
    // Clear base class.
    ParticleModel::releaseMem();

    // Delete inceptions.
    for (IcnPtrVector::iterator i=m_inceptions.begin();
         i!=m_inceptions.end(); ++i) {
        delete *i;
    }
    m_inceptions.clear();

    // Delete processes.
    for (PartProcPtrVector::iterator i=m_processes.begin();
         i!=m_processes.end(); ++i) {
        delete *i;
    }
    m_processes.clear();

    // Delete coagulation processes.
    for (CoagPtrVector::iterator i = m_coags.begin();
         i != m_coags.end(); ++i) {
        delete *i;
    }
    m_coags.clear();
    m_icoag = -1;

    m_anydeferred = false;
    m_termcount = 0;
    m_processcount = 0;
    m_proccount.clear();
    m_fictcount.clear();
}
