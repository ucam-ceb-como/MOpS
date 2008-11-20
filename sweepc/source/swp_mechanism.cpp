/*
  Author(s):      Matthew Celnik (msc37)
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
#include "swp_abf_model.h"
#include <stdexcept>

//added to test 3-d output
//#include "swp_imgnode.h"
#include "swp_particle_image.h"
#include "swp_ensemble_image.h"
#include "string_functions.h"


using namespace Sweep;
using namespace Sweep::Processes;
using namespace Sweep::Imaging;
using namespace std;
using namespace Strings;


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Mechanism::Mechanism(void)
: m_anydeferred(false), m_coag(NULL), m_icoag(-1), m_termcount(0), m_processcount(0)
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
        }

        // Copy particle processes.
        for (PartProcPtrVector::const_iterator i=rhs.m_processes.begin();
             i!=rhs.m_processes.end(); ++i) {
            m_processes.push_back((*i)->Clone());
        }

        // Copy coagulation process.
        if (rhs.m_coag) m_coag = rhs.m_coag->Clone();

        // Copy process counters.
        m_proccount.assign(rhs.m_proccount.begin(), rhs.m_proccount.end());
        m_fictcount.assign(rhs.m_fictcount.begin(), rhs.m_fictcount.end());
    }
    return *this;
}


// ACTIVE-SITES MODELS.

// Returns the set of particle model ID used by this mechanism
const ActSites::ActSitesTypeSet &Mechanism::ActSiteModels(void) const
{
    return m_actsites;
}

// Returns true if the mechanism include the given model.
bool Mechanism::ContainsActSiteModel(ActSites::ActSitesType id) const
{
    return m_actsites.find(id) != m_actsites.end();
}

// Adds an active-sites model to the mechanism.
void Mechanism::AddActSitesModel(ActSites::ActSitesType id)
{
    m_actsites.insert(id);
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
    // Add the inception to the mechanism.
    m_inceptions.push_back(&icn);
    m_termcount += icn.TermCount();
    ++m_processcount;
    m_proccount.resize(m_termcount, 0);
    m_fictcount.resize(m_termcount, 0);

    // Set the inception to belong to this mechanism.
    icn.SetMechanism(*this);
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

// Adds a coagulation process to the mechanism.
void Mechanism::AddCoagulation()
{
    if (m_coag != NULL) {
        m_termcount -= m_coag->TermCount();
    } else {
        ++m_processcount;
    }
    delete m_coag;
    m_coag = new Coagulation(*this);
    m_termcount += m_coag->TermCount();
    m_proccount.resize(m_termcount, 0);
    m_fictcount.resize(m_termcount, 0);
}


// PROCESS INFORMATION.

// Returns the number of processes (including 
// inceptions) in the mechanism.
unsigned int Mechanism::ProcessCount(void) const
{
/*
    unsigned int n = m_inceptions.size() + m_processes.size();
    if (m_coag != NULL) ++n;
    return n;
*/
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
    if (m_coag) {
        *i = m_coag->Name(); ++i;
    }
}


// RATE CALCULATION.

// Get total rates of all processes.  Returns the sum of
// all rates.
real Mechanism::CalcRates(real t, const Cell &sys, fvector &rates, bool scale) const
{
    // Ensure rates vector is the correct length, then set to zero.
    rates.resize(m_processcount+sys.InflowCount()+sys.OutflowCount(), 0.0);
    fill(rates.begin(), rates.end(), 0.0);

    real sum = 0.0;

    // Get rates of inception processes.
    sum += Inception::CalcRates(t, sys, m_inceptions, rates);

    // Query other processes for their rates.
    sum += ParticleProcess::CalcRates(t, sys, m_processes, rates, m_inceptions.size());

    // Get coagulation rate.
    fvector::iterator i = rates.begin()+m_inceptions.size()+m_processes.size();
    if (m_coag != NULL) {
        *i = m_coag->Rate(t, sys);
        sum += *i;
        ++i;
    }

    // Get birth and death rates from the Cell.
    const BirthPtrVector &inf = sys.Inflows();
    for (BirthPtrVector::const_iterator j=inf.begin(); j!=inf.end(); ++j) {
        *i = (*j)->Rate(t, sys);
        sum += *i;
        ++i;
    }
    const DeathPtrVector &out = sys.Outflows();
    for (DeathPtrVector::const_iterator j=out.begin(); j!=out.end(); ++j) {
        *i = (*j)->Rate(t, sys);
        sum += *i;
        ++i;
    }

    if (!scale) {
        // Need to return the rates to per unit vol.
        real invvol = 1.0 / sys.SampleVolume();
        for(i=rates.begin(); i!=rates.end(); ++i) {
            *i *= invvol;
        }
        sum *= invvol;
    }

    return sum;
}

// Get rates of all processes separated into different
// terms.  Rate terms are useful for subsequent particle
// selection by different properties for the same process.
// In particular this is used for the condensation and
// coagulation processes.  Returns the sum of all rates.
real Mechanism::CalcRateTerms(real t, const Cell &sys, fvector &terms) const
{
    // Ensure rates vector is the correct length.
    terms.resize(m_termcount+sys.InflowCount()+sys.OutflowCount(), 0.0);
    fvector::iterator iterm = terms.begin();

    real sum = 0.0;

    // Get rates of inception processes.
    IcnPtrVector::const_iterator ii;
    for (ii=m_inceptions.begin(); ii!=m_inceptions.end(); ++ii) {
        sum += (*ii)->RateTerms(t, sys, iterm);
    }

    // Query other processes for their rates.
    if (sys.ParticleCount() > 0) {
        PartProcPtrVector::const_iterator i;
        for(i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=terms.end()); ++i) {
            sum += (*i)->RateTerms(t, sys, iterm);
        }
    } else {
        // Fill vector with zeros.
        PartProcPtrVector::const_iterator i;
        for(i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=terms.end()); ++i) {
            fill(iterm, iterm+(*i)->TermCount(), 0.0);
            iterm += (*i)->TermCount();
        }
    }

    // Get coagulation rate.
    if (m_coag != NULL) {
        sum += m_coag->RateTerms(t, sys, iterm);
    }

    // Get birth and death rates from the Cell.
    const BirthPtrVector &inf = sys.Inflows();
    for (BirthPtrVector::const_iterator j=inf.begin(); j!=inf.end(); ++j) {
        sum += (*j)->RateTerms(t, sys, iterm);
    }
    const DeathPtrVector &out = sys.Outflows();
    for (DeathPtrVector::const_iterator j=out.begin(); j!=out.end(); ++j) {
        sum += (*j)->RateTerms(t, sys, iterm);
    }

    return sum;
}

// Get total rates of non-deferred processes.  Returns the sum
// of all rates.
real Mechanism::CalcJumpRateTerms(real t, const Cell &sys, fvector &terms) const
{
    // This routine only calculates the rates of those processes which are
    // not deferred.  The rate terms of deferred processes are returned
    // as zero.

    // Ensure vector is the correct length, then set to zero.
    terms.resize(m_termcount+sys.InflowCount()+sys.OutflowCount(), 0.0);
    fill(terms.begin(), terms.end(), 0.0);
    fvector::iterator iterm = terms.begin();

    real sum = 0.0;

    // Get rates of inception processes.
    IcnPtrVector::const_iterator ii;
    for (ii=m_inceptions.begin(); ii!=m_inceptions.end(); ++ii) {
        sum += (*ii)->RateTerms(t, sys, iterm);
    }

    // Query other processes for their rates.
    if (sys.ParticleCount() > 0) {
        PartProcPtrVector::const_iterator i;
        for(i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=terms.end()); ++i) {
            if (!(*i)->IsDeferred()) {
                // Calculate rate if not deferred.
                sum += (*i)->RateTerms(t, sys, iterm);
            } else {
                // If process is deferred, then set rate to zero.
                for (unsigned int j=0; j!=(*i)->TermCount(); ++j) {*(iterm++)=0.0;}
            }
        }
    } else {
        // Fill vector with zeros.
        PartProcPtrVector::const_iterator i;
        for(i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=terms.end()); ++i) {
            fill(iterm, iterm+(*i)->TermCount(), 0.0);
            iterm += (*i)->TermCount();
        }
    }

    // Get coagulation rate.
    if (m_coag != NULL) {
        sum += m_coag->RateTerms(t, sys, iterm);
    } 

    // Get birth and death rates from the Cell.
    const BirthPtrVector &inf = sys.Inflows();
    for (BirthPtrVector::const_iterator j=inf.begin(); j!=inf.end(); ++j) {
        sum += (*j)->RateTerms(t, sys, iterm);
    }
    const DeathPtrVector &out = sys.Outflows();
    for (DeathPtrVector::const_iterator j=out.begin(); j!=out.end(); ++j) {
        sum += (*j)->RateTerms(t, sys, iterm);
    }

    return sum;
}

/*
// Get total rates of non-deferred processes.  Returns the sum
// of all rates.  Uses supplied gas-phase conditions rather than
// those in the given system.
real Mechanism::CalcJumpRates(real t, const Sprog::Thermo::IdealGas &gas, 
                              const Cell &sys, fvector &rates) const
{
    // This routine only calculates the rates of those processes which are
    // not deferred.  The rate terms of deferred processes are returned
    // as zero.

    // Ensure vector is the correct length.
    rates.resize(m_termcount);
    fvector::iterator iterm = rates.begin()+m_inceptions.size();

    real sum = 0.0;

    // Get rates of inception processes.
    sum += Inception::CalcRates(t, gas, sys, m_inceptions, rates, 0);

    //IcnPtrVector::const_iterator ii;
    //for (ii=m_inceptions.begin(); ii!=m_inceptions.end(); ++ii) {
    //    sum += (*ii)->RateTerms(t, gas, sys, iterm);
    //}

    // Query other processes for their rates.
    if (sys.ParticleCount() > 0) {
        PartProcPtrVector::const_iterator i;
        for(i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=rates.end()); ++i) {
            if (!(*i)->IsDeferred()) {
                // Calculate rate if not deferred.
                sum += (*i)->RateTerms(t, gas, sys, iterm);
            } else {
                // If process is deferred, then set rate to zero.
                for (unsigned int j=0; j!=(*i)->TermCount(); ++j) {*(iterm++)=0.0;}
            }
        }
    }

    // Get coagulation rate.
    if (m_coag != NULL) {
        sum += m_coag->RateTerms(t, gas, sys, iterm);
    }

    // Get death process rate.
    if (m_death != NULL) {
        sum += m_death->RateTerms(t, sys, iterm);
    }

    return sum;
}
*/

// Calculates the rates-of-change of the chemical species fractions, 
// gas-phase temperature and density due to particle processes.
void Mechanism::CalcGasChangeRates(real t, const Cell &sys, fvector &rates) const
{
    // Resize vector to hold all species and set all rates to zero.
    rates.resize(m_species->size()+2, 0.0);
    fill(rates.begin(), rates.end(), 0.0);
    fvector::iterator idrho = rates.begin() + (rates.size()-1);

    // Precalculate parameters.
    real invVolNA = 1.0 / (sys.SampleVolume() * NA);

    // Inceptions and surface processes can affect the gas-phase chemistry
    // at the moment.

    // Loop over the contributions of all inception processes.
    for (IcnPtrVector::const_iterator i=m_inceptions.begin(); 
         i!=m_inceptions.end(); ++i) {

        // Calculate the inception rate.
        real rate = (*i)->Rate(t, sys);

        // Loop over all reactants, subtracting their contributions.
        for (Sprog::StoichMap::const_iterator j=(*i)->Reactants().begin(); 
             j!=(*i)->Reactants().end(); ++j) {
            real dc = rate * (real)j->second * invVolNA;
            rates[j->first] -= dc;
            *idrho -= dc;
        }
        
        // Loop over all products, adding their contributions.
        for (Sprog::StoichMap::const_iterator j=(*i)->Products().begin(); 
             j!=(*i)->Products().end(); ++j) {
            real dc = rate * (real)j->second * invVolNA;
            rates[j->first] += dc;
            *idrho += dc;
        }
    }

    // Loop over the contributions of all other processes (except coagulation).
    for (PartProcPtrVector::const_iterator i=m_processes.begin(); 
         i!=m_processes.end(); ++i) {

        // Calculate the process rate.
        real rate = (*i)->Rate(t, sys);

        // Loop over all reactants, subtracting their contributions.
        for (Sprog::StoichMap::const_iterator j=(*i)->Reactants().begin(); 
             j!=(*i)->Reactants().end(); ++j) {
            real dc = rate * (real)j->second * invVolNA;
            rates[j->first] -= dc;
            *idrho -= dc;
        }
        
        // Loop over all products, adding their contributions.
        for (Sprog::StoichMap::const_iterator j=(*i)->Products().begin(); 
             j!=(*i)->Products().end(); ++j) {
            real dc = rate * (real)j->second * invVolNA;
            rates[j->first] += dc;
            *idrho += dc;
        }
    }

    // Now convert to changes in mole fractions.
    real invrho = 1.0 / sys.Density();
    for (unsigned int k=0; k!=m_species->size(); ++k) {
        rates[k] = (invrho * rates[k]) - (invrho * sys.MoleFraction(k) * (*idrho));
    }
}


// PERFORMING THE PROCESSES.

// Performs the Process specified.  Process index could be
// an inception, particle process or a coagulation event.
void Mechanism::DoProcess(unsigned int i, real t, Cell &sys) const
{
    // Work out to which process this term belongs.
    int j = i - m_inceptions.size();

    if (j < 0) {
        // This is an inception process.
        m_inceptions[i]->Perform(t, sys, 0);
        m_proccount[i] += 1;
    } else {
        // This is another process. 
        PartProcPtrVector::const_iterator ip;
        for(ip=m_processes.begin(); ip!=m_processes.end(); ++ip) {
            if (j < (int)(*ip)->TermCount()) {
                // Do the process.
                if ((*ip)->Perform(t, sys, j) == 0) {
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
        if (m_coag) {
            // Check if coagulation process.
            if (j < (int)m_coag->TermCount()) {
                // This is the coagulation process.
                if (m_coag->Perform(t, sys, j) == 0) {
                    m_proccount[i] += 1;
                } else {
                    m_fictcount[i] += 1;
                }
                return;
            } else {
                // This must be the birth/death process.
                j -= m_coag->TermCount();
            }
        }

        // The coagulation process is undefined, so it 
        // must be either a birth or death process from
        // the cell.
        if ((j < (int)sys.InflowCount()) && (j>=0)) {
            // An inflow process.
            sys.Inflows(j)->SetMechanism(*this);
            sys.Inflows(j)->Perform(t, sys);
        } else {
            j -= sys.InflowCount();
            if ((j < (int)sys.OutflowCount()) && (j>=0)) {
                // An outflow process.
                sys.Outflows(j)->SetMechanism(*this);
                sys.Outflows(j)->Perform(t, sys);
            }
        }
    }
}


// LINEAR PROCESS DEFERMENT ALGORITHM.

// Performs linear update algorithm on the 
// given system up to given time.

void Mechanism::LPDA(real t, Cell &sys) const
{	
    // Check that there are particles to update and that there are
    // deferred processes to perform.
    if ((sys.ParticleCount() > 0) && (m_anydeferred)) {
        // Stop ensemble from doubling while updating particles.
        sys.Particles().FreezeDoubling();

        // Perform deferred processes on all particles individually.
        Ensemble::iterator i;
        unsigned int k = 0;
        for (i=sys.Particles().begin(); i!=sys.Particles().end(); ++i) {
            UpdateParticle(*(*i), sys, t);
            ++k;
        }

        
        // Now remove any invalid particles and update the ensemble.
        sys.Particles().RemoveInvalids();
        sys.Particles().Update();

        // Start particle doubling again.  This will also double the ensemble
        // if too many particles have been removed.
        sys.Particles().UnfreezeDoubling();
    }
}

/*
// Performs linear update algorithm on the given system up to given time,
// with the given chemical conditions rather than those in the 
// given system.
void Mechanism::LPDA(real t, const Sprog::Thermo::IdealGas &gas, 
                     Cell &sys) const
{
    // Check that there are particles to update and that there are
    // deferred processes to perform.
    if ((sys.ParticleCount() > 0) && (m_anydeferred)) {
        // Stop ensemble from doubling while updating particles.
        sys.Particles().FreezeDoubling();

        // Perform deferred processes on all particles individually.
        Ensemble::iterator i;
        unsigned int k = 0;
        for (i=sys.Particles().begin(); i!=sys.Particles().end(); ++i) {
            UpdateParticle(*(*i), gas, sys, t);
            ++k;
        }
        
        // Now remove any invalid particles and update the ensemble.
        sys.Particles().RemoveInvalids();
        sys.Particles().Update();

        // Start particle doubling again.  This will also double the ensemble
        // if too many particles have been removed.
        sys.Particles().UnfreezeDoubling();
    }
}
*/

// Performs linear process updates on a particle in the given system.
void Mechanism::UpdateParticle(Particle &sp, Cell &sys, real t) const
{
    // If there are no deferred processes then stop right now.
    if (m_anydeferred) {
        PartProcPtrVector::const_iterator i;
        unsigned int num;
        real rate, dt;

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
                    // times to perform the process.
                    num = ignpoi(rate);

                    if (num > 0) {
                        // Do the process to the particle.
                        (*i)->Perform(t, sys, sp, num);
                    }
                }
            }

            // Perform sintering update.
            if (m_sint_model.IsEnabled()) {
				sp.UpdateFreeSurface();

				//	 sp.FindRoot()->CheckTree();
			    // cout << "check before sinter passed\n";
			    sp.Sinter(dt, sys, m_sint_model);
    			//	sp.FindRoot()->CheckTree();
				// cout << "check after sinter passed\n";
				//added to test the 3-d output		
				  /*   if(sp.Property(ParticleCache::iM)>4E-23)   
					  {
						    Sweep::Imaging::ParticleImage img;
							img.Construct(sp);
							ofstream filepart; filepart.open("part.3d");
							img.Write3dout(filepart,0,0,0);
					   }*/
						 if (sys.ParticleCount() > 10 && t-last3dout>0.001)
						{	ofstream file;
						    string fname;
							Ensemble::iterator i;
							cout << "creating subpart image...";
							last3dout=t;
							fname = "subpart1" + cstr(t) + ".3d";
							 file.open(fname.c_str());
							i=sys.Particles().begin();
							Sweep::Imaging::ParticleImage img;
							img.Construct(*(*i));
							img.Write3dout(file,0,0,0);
							file.close();
							EnsembleImage *ensimg = new EnsembleImage;
							fname = "subparttree" + cstr(t) + ".3d";
							 file.open(fname.c_str());
							ensimg->PrintEnsemble(sys,file);
							file.close();
							delete ensimg;
							
				
							
							int numpart=0;
							int numsubpart=0;
							real avcoldiam=0.;
							numpart=0;
							numsubpart=0;
							int coldiamdistr[500]={0};
							fname = "avcoldiam.txt";
							file.open(fname.c_str(),ios::app);
					//		for (i=sys.Particles().begin(); i!=sys.Particles().end(); ++i) {
							for (i=sys.Particles().begin(); i!=sys.Particles().end() ; ++i) {
								avcoldiam=(*(*i)).CollDiameter()+avcoldiam;
								numpart++;
								numsubpart+=(*(*i)).NumSubPart();
								int intdiam=(int)((*(*i)).CollDiameter()*1e9);
								coldiamdistr[intdiam]++;
							}
							if (numpart>0)
							{	
								avcoldiam=avcoldiam/numpart;
								file << t << "     " << avcoldiam << endl;
							}
							file.close();
							fname = "coldiamdistr" + cstr(t) + ".txt";
							file.open(fname.c_str());
							for (int j=0;j<500;j++)
							{
								file<<j<< "    "<<coldiamdistr[j]<<endl;
							}
							file.close();
						


							fname = "numsubpart.txt";
							file.open(fname.c_str(),ios::app);
							if (numpart>0)
							{	
								file << t << "     " << numsubpart << endl;
							}

							file.close();


							real avsubpartdiam=0.;
							real avsubpartdiam2=1.;
							int avdiamdistr[500]={0};
							fname = "avdiam.txt";
							file.open(fname.c_str(),ios::app);
					//		for (i=sys.Particles().begin(); i!=sys.Particles().end(); ++i) {
							for (i=sys.Particles().begin(); i!=sys.Particles().end() ; ++i) {
								avsubpartdiam2=(*(*i)).avgeomdiam(1.0/numsubpart)*avsubpartdiam2;
								avsubpartdiam=(*(*i)).Volume()/(*(*i)).SurfaceArea()+avsubpartdiam;
							
							}
							//avsubpartdiam=pow(6*avsubpartdiam/(PI*numsubpart),ONE_THIRD);
							avsubpartdiam=6*avsubpartdiam/numpart;
							if (numpart>0)
							{	
								file << t << "     " << avsubpartdiam << endl;
							}
							file.close();





							
							fname = "temperature.txt";
							file.open(fname.c_str(),ios::app);
					//		for (i=sys.Particles().begin(); i!=sys.Particles().end(); ++i) {

							file << t << "     " << sys.Temperature() << endl;
							
							file.close();

							cout << "done"<<endl;
	
						}

			
            }
        }

        // Check that the particle is still valid, only calculate 
        // cache if it is.
        if (sp.IsValid()) sp.UpdateCache();
    }
}

/*
// Performs linear process updates on a particle in the given system,
// with the current chemical conditions precalculated.
void Mechanism::UpdateParticle(Particle &sp, const Sprog::Thermo::IdealGas &gas, 
                               Cell &sys, real t) const
{
    // If there are no deferred processes then stop right now.
    if (m_anydeferred) {
        PartProcPtrVector::const_iterator i;
        unsigned int num;
        real rate, dt;

        while ((sp.LastUpdateTime() < t) && sp.IsValid()) {
            // Calculate delta-t and update particle time.
            dt = t - sp.LastUpdateTime();
            sp.SetTime(t);

            // Loop through all processes, performing those
            // which are deferred.
            for (i=m_processes.begin(); i!=m_processes.end(); ++i) {
                if ((*i)->IsDeferred()) {
                    // Get the process rate x the time interval.
                    rate = (*i)->Rate(t, gas, sys, sp) * dt;

                    // Use a Poission deviate to calculate number of 
                    // times to perform the process.
                    num = ignpoi(rate);

                    if (num > 0) {
                        // Do the process to the particle.
                        (*i)->Perform(t, sys, sp, num);
                    }
                }
            }

            // Perform sintering update.
            if (m_sint_model.IsEnabled()) {
                sp.Sinter(dt, sys, m_sint_model);
            }
        }

        // Check that the particle is still valid, only calculate 
        // cache if it is.
        if (sp.IsValid()) sp.UpdateCache();
    }
}
*/


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

        // Write coagulation process.
        if (m_coag != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            ProcessFactory::Write(*m_coag, out);
        } else {
            out.write((char*)&falseval, sizeof(falseval));
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
                if (n == 1) {
                    m_coag = ProcessFactory::ReadCoag(in, *this);
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

    // Delete coagulation process.
    delete m_coag; m_coag = NULL;
    m_icoag = -1;

    m_anydeferred = false;
    m_termcount = 0;
    m_processcount = 0;
    m_proccount.clear();
    m_fictcount.clear();
}
