#include "swp_mechanism.h"
#include "swp_modelfactory.h"
#include "swp_processfactory.h"
#include "swp_abfmodel.h"
#include <stdexcept>

using namespace Sweep;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Mechanism::Mechanism(void)
: m_anydeferred(false), m_icoag(-1), m_termcount(0), 
  m_coag(NULL), m_species(NULL)
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

        // Easily copied data.
        m_anydeferred = rhs.m_anydeferred;
        m_species     = rhs.m_species;
        m_icoag       = rhs.m_icoag;
        m_termcount   = rhs.m_termcount;

        // Copy components.
        for (CompPtrVector::const_iterator i=rhs.m_components.begin();
             i!=rhs.m_components.end(); ++i) {
            m_components.push_back((*i)->Clone());
        }

        // Copy trackers.
        for (TrackPtrVector::const_iterator i=rhs.m_trackers.begin();
             i!=rhs.m_trackers.end(); ++i) {
            m_trackers.push_back((*i)->Clone());
        }

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
        m_coag = rhs.m_coag->Clone();

        // Copy particle model info.
        for (ModelTypeSet::const_iterator i=rhs.m_models.begin();
            i!=rhs.m_models.end(); ++i) {
            m_models.insert(*i);
        }

        // Copy process counters.
        m_proccount.assign(rhs.m_proccount.begin(), rhs.m_proccount.end());
        m_fictcount.assign(rhs.m_fictcount.begin(), rhs.m_fictcount.end());
    }
    return *this;
}


// CHEMICAL SPECIES.

// Returns the chemical species vector.
const Sprog::SpeciesPtrVector *const Mechanism::Species(void) const
{
    return m_species;
}

// Sets the chemical species vector.
void Mechanism::SetSpecies(const Sprog::SpeciesPtrVector &sp)
{
    m_species = &sp;
}


// COMPONENT DEFINITIONS.

// Returns the number of components in the mechanism.
unsigned int Mechanism::ComponentCount(void) const
{
    return m_components.size();
}

// Returns the vector of particle components.
const CompPtrVector &Mechanism::Components() const
{
    return m_components;
}

// Returns the component with the given index.
const Component *const Mechanism::Components(unsigned int i) const
{
    if (i < m_components.size()) {
        return m_components[i];
    } else {
        return NULL;
    }
}

// Returns the index of the component with the 
// given name in the mechanism if found, otherwise 
// return negative.
int Mechanism::ComponentIndex(const std::string &name) const
{
    for (unsigned int i=0; i!=m_components.size(); ++i) {
        if (name.compare(m_components[i]->Name())==0) {
            return i;
        }
    }
    return -1;
}

// Adds a component to the mechanism and returns the index
// of the component.
unsigned int Mechanism::AddComponent(Component &comp)
{
    m_components.push_back(&comp);
    return m_components.size()-1;
}

// Overwrites the ith component with that given.  Previous
// component is deleted from memory.
void Mechanism::ReplaceComponent(unsigned int i, Component &comp)
{
    if (i < m_components.size()) {
        delete m_components[i];
        m_components[i] = &comp;
    }
}

// Sets the particle components vector.
void Mechanism::SetComponents(const CompPtrVector &comps)
{
    // Delete current components.
    CompPtrVector::iterator i;
    for (i=m_components.begin(); i!=m_components.end(); ++i) {
        delete (*i);
    }

    // Resize component vector and copy components.
    m_components.resize(comps.size());
    CompPtrVector::const_iterator ic;
    for (ic=comps.begin(); ic!=comps.end(); ++ic) {
        m_components.push_back((*ic)->Clone());
    }
}


// TRACKER VARIABLES.

// Returns the number of tracker variables.
unsigned int Mechanism::TrackerCount(void) const 
{
    return m_trackers.size();
}

// Returns the vector of tracker variables.
const TrackPtrVector &Mechanism::Trackers(void) const
{
    return m_trackers;
}

// Returns the ith tracker variable.
const Tracker *const Mechanism::Trackers(unsigned int i) const
{
    if (i < m_trackers.size()) {
        return m_trackers[i];
    } else {
        return NULL;
    }
}

// Returns the index of the tracker variable with the given name 
// on success, otherwise returns negative.
int Mechanism::GetTrackerIndex(const std::string &name) const
{
    for (unsigned int i=0; i!=m_trackers.size(); ++i) {
        if (name.compare(m_trackers[i]->Name())==0) {
            return i;
        }
    }
    return -1;
}

// Adds a tracker variable to the mechanism.
void Mechanism::AddTracker(Tracker &track)
{
    m_trackers.push_back(&track);
}

// Replaces the tracker at the given index with the given tracker
// object.
void Mechanism::ReplaceTracker(unsigned int i, Tracker &track)
{
    if (i < m_trackers.size()) {
        delete m_trackers[i];
        m_trackers[i] = &track;
    }
}

// Sets the vector of tracker variables.
void Mechanism::SetTrackers(const TrackPtrVector &track)
{
    // Delete current trackers.
    TrackPtrVector::iterator i;
    for (i=m_trackers.begin(); i!=m_trackers.end(); ++i) {
        delete (*i);
    }

    // Resize tracker vector and copy tracker variables.
    m_trackers.resize(track.size());
    TrackPtrVector::const_iterator ic;
    for (ic=track.begin(); ic!=track.end(); ++i) {
        m_trackers.push_back((*ic)->Clone());
    }
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
    m_proccount.resize(m_termcount, 0.0);
    m_fictcount.resize(m_termcount, 0.0);

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
    m_proccount.resize(m_termcount, 0.0);
    m_fictcount.resize(m_termcount, 0.0);

    // Check for any deferred.
    m_anydeferred = m_anydeferred || p.IsDeferred();

    // Set the process to belong to this mechanism.
    p.SetMechanism(*this);
}


// COAGULATIONS.

// Adds a coagulation process to the mechanism.
void Mechanism::AddCoagulation()
{
    if (m_coag != NULL) m_termcount -= m_coag->TermCount();
    delete m_coag;
    m_coag = new Coagulation();
    m_coag->SetMechanism(*this);
    m_termcount += m_coag->TermCount();
    m_proccount.resize(m_termcount, 0.0);
    m_fictcount.resize(m_termcount, 0.0);
}


// PROCESS INFORMATION.

// Returns the number of processes (including 
// inceptions) in the mechanism.
unsigned int Mechanism::ProcessCount(void) const
{
    unsigned int n = m_inceptions.size() + m_processes.size();
    if (m_coag != NULL) ++n;
    return n;
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


// PARTICLE MODELS.

// Returns the set of particle model ID used by this mechanism
const ModelTypeSet &Mechanism::Models(void) const
{
    return m_models;
}

// Returns true if the mechanism include the given model.
bool Mechanism::ContainsModel(ModelType id) const
{
    return m_models.find(id) != m_models.end();
}

// Adds a model to the mechanism.  Any subsequent particles
// created with this mechanism will use this model.
void Mechanism::AddModel(ModelType id)
{
    m_models.insert(id);

    // Initialise the model as necessary.
    switch(id) {
        case ABFSites_ID:
            ABFModel::Instance().Initialise(*this);
            break;
    }
}


// RATE CALCULATION.

// Get total rates of all processes.  Returns the sum of
// all rates.
real Mechanism::CalcRates(real t, const Cell &sys, fvector &rates) const
{
    // Ensure rates vector is the correct length.
    rates.resize(m_termcount);
    fvector::iterator iterm = rates.begin();

    real sum = 0.0;

    // Get rates of inception processes.
    IcnPtrVector::const_iterator ii;
    for (ii=m_inceptions.begin(); ii!=m_inceptions.end(); ++ii) {
        sum += (*ii)->RateTerms(t, sys, iterm);
    }

    // Query other processes for their rates.
    if (sys.ParticleCount() > 0) {
        PartProcPtrVector::const_iterator i;
        for(i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=rates.end()); ++i) {
            sum += (*i)->RateTerms(t, sys, iterm);
        }
    }

    // Get coagulation rate.
    if (m_coag != NULL) {
        sum += m_coag->RateTerms(t, sys, iterm);
    }

    return sum;
}

// Get total rates of non-deferred processes.  Returns the sum
// of all rates.
real Mechanism::CalcJumpRates(real t, const Cell &sys, fvector &rates) const
{
    // This routine only calculates the rates of those processes which are
    // not deferred.  The rate terms of deferred processes are returned
    // as zero.

    // Ensure vector is the correct length.
    rates.resize(m_termcount);
    fvector::iterator iterm = rates.begin();

    real sum = 0.0;

    // Get rates of inception processes.
    IcnPtrVector::const_iterator ii;
    for (ii=m_inceptions.begin(); ii!=m_inceptions.end(); ++ii) {
        sum += (*ii)->RateTerms(t, sys, iterm);
    }

    // Query other processes for their rates.
    if (sys.ParticleCount() > 0) {
        PartProcPtrVector::const_iterator i;
        for(i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=rates.end()); ++i) {
            if (!(*i)->IsDeferred()) {
                // Calculate rate if not deferred.
                sum += (*i)->RateTerms(t, sys, iterm);
            } else {
                // If process is deferred, then set rate to zero.
                for (unsigned int j=0; j!=(*i)->TermCount(); ++j) {*(iterm++)=0.0;}
            }
        }
    }

    // Get coagulation rate.
    if (m_coag != NULL) {
        sum += m_coag->RateTerms(t, sys, iterm);
    }

    return sum;
}

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

    return sum;
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
            if (j < (*ip)->TermCount()) {
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
        // nor a single particle process.  It is therefore a 
        // coagulation process.
        if (m_coag->Perform(t, sys, j) == 0) {
            m_proccount[i] += 1;
        } else {
            m_fictcount[i] += 1;
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

// Performs linear process updates on a particle in the given system.
void Mechanism::UpdateParticle(Particle &sp, Cell &sys, real t) const
{
    UpdateParticle(sp, sys, sys, t);
}

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
                    // Get the process rate times the time interval.
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
        }

        // Check that the particle is still valid, only calculate 
        // cache if it is.
        if (sp.IsValid()) sp.UpdateCache();
    }
}


// PARTICLE FUNCTIONS.

// Creates a new particle and sets it up with all the models
// required by the mechanism.
Particle *const Mechanism::CreateParticle(void) const
{
    // Create new particle using this mechanism's components
    // and tracker variables.
    Particle *p = new Particle(m_components, m_trackers);

    // Add particle models.
    for (ModelTypeSet::const_iterator id=m_models.begin(); id!=m_models.end(); ++id) {
        IModelData *model = ModelFactory::CreateData(*id, *p);
        if (model != NULL) p->AddModel(*model);
    }

    // Returns particle.
    return p;
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

        // Write if any processes are deferred.
        if (m_anydeferred) {
            out.write((char*)&trueval, sizeof(trueval));
        } else {
            out.write((char*)&falseval, sizeof(falseval));
        }

        // Write number of components.
        unsigned int n = (unsigned int)m_components.size();
        out.write((char*)&n, sizeof(n));

        // Write components.
        for (CompPtrVector::const_iterator i=m_components.begin(); 
             i!=m_components.end(); ++i) {
            (*i)->Serialize(out);
        }

        // Write number of trackers.
        n = (unsigned int)m_trackers.size();
        out.write((char*)&n, sizeof(n));

        // Write trackers.
        for (TrackPtrVector::const_iterator i=m_trackers.begin(); 
             i!=m_trackers.end(); ++i) {
            (*i)->Serialize(out);
        }

        // Write number of inceptions.
        n = (unsigned int)m_inceptions.size();
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

        // Write model count.
        n = (unsigned int)m_models.size();
        out.write((char*)&n, sizeof(n));

        // Write model set.
        for (ModelTypeSet::const_iterator i=m_models.begin(); i!=m_models.end(); ++i) {
            n = (unsigned int)(*i);
            out.write((char*)&n, sizeof(n));
        }
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

        int m = 0;
        unsigned int n =0, id = 0;

        switch (version) {
            case 0:
                // Read if any processes are deferred.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_anydeferred = (n==1);

                // Read number of components.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read components.
                for (unsigned int i=0; i!=n; ++i) {
                    m_components.push_back(new Component(in));
                }

                // Read number of trackers.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read trackers.
                for (unsigned int i=0; i!=n; ++i) {
                    m_trackers.push_back(new Tracker(in));
                }

                // Read number of inceptions.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read inceptions.
                for (unsigned int i=0; i!=n; ++i) {
                    Inception *icn = ProcessFactory::ReadInception(in);
                    icn->SetMechanism(*this);
                    m_inceptions.push_back(icn);
                }

                // Read number of particle processes.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read particle processes.
                for (unsigned int i=0; i!=n; ++i) {
                    ParticleProcess *p = ProcessFactory::ReadPartProcess(in);
                    p->SetMechanism(*this);
                    m_processes.push_back(p);
                }

                // Read coagulation process.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_coag = ProcessFactory::ReadCoag(in);
                }

                // Read index of first coag process.
                in.read(reinterpret_cast<char*>(&m), sizeof(m));
                m_icoag = m;

                // Read term count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                m_termcount = n;

                // Read model count.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read model set.
                for (unsigned int i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&id), sizeof(id));
                    m_models.insert((ModelType)id);
                }

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
    // Delete components.
    for (CompPtrVector::iterator i=m_components.begin(); 
         i!=m_components.end(); ++i) {
        delete *i;
    }
    m_components.clear();

    // Delete trackers.
    for (TrackPtrVector::iterator i=m_trackers.begin(); 
         i!=m_trackers.end();++i) {
        delete *i;
    }
    m_trackers.clear();

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
    delete m_coag;
}
