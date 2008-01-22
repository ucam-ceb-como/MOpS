#include "swpmechanism.h"
#include "swpcoagulation.h"
#include "swpensemble.h"
#include "rng.h"
#include "swpabf.h"
#include "swppahmodel.h"

using namespace Sweep;

Sweep::Mechanism::Mechanism(void)
{
    m_anydeferred = false;
    m_icoag = -1;
    m_components.clear();
    m_termcount = 0;
    m_pmodel = SphericalParticle;
}

Sweep::Mechanism::~Mechanism(void)
{
    // Delete components.
    vector<Component*>::iterator icmp;
    for (icmp=m_components.begin(); icmp!=m_components.end(); icmp++) {
        delete *icmp;
    }
    m_components.clear();

    // Delete inceptions.
    vector<Inception*>::iterator ii;
    for (ii=m_inceptions.begin(); ii!=m_inceptions.end(); ii++) {
        delete *ii;
    }
    m_inceptions.clear();

    // Delete processes.
    vector<Process*>::iterator ip;
    for (ip=m_processes.begin(); ip!=m_processes.end(); ip++) {
        delete *ip;
    }
    m_processes.clear();

    // Clear value names.
    m_valuenames.clear();
}

int Sweep::Mechanism::GetRates(std::vector<real> &rates, const real t, 
                               const System &sys) const
{
    // Ensure rates vector is the correct length.
    rates.resize(m_termcount);
    vector<real>::iterator iterm = rates.begin();

    // Get chemical conditions.
    vector<real> chem, sums;
    real T=0.0, P=0.0;
    sys.GetConditions(t, chem, T, P);
    sys.ConstEnsemble().GetSums(sums);

    // Get rates of inception processes.
    vector<Inception*>::const_iterator ii;
    for (ii=m_inceptions.begin(); ii!=m_inceptions.end(); ii++) {
        (*ii)->RateTerms(t, chem, T, P, sums, sys, iterm);
    }

    // Query other processes for their rates.
    vector<Process*>::const_iterator i;
    for(i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=rates.end()); i++) {
        (*i)->RateTerms(t, chem, T, P, sums, sys, iterm);
    }

    return 0;
}

int Mechanism::GetStochasticRates(std::vector<real> &rates, const real t, 
                                  const System &sys) const
{
    // This routine only calculates the rates of those processes which ar
    // not deferred.

    // Ensure vector is the correct length.
    rates.resize(m_termcount);
    vector<real>::iterator iterm = rates.begin();

    // Get chemical conditions.
    vector<real> chem, sums;
    real T=0.0, P=0.0;
    sys.GetConditions(t, chem, T, P);
    sys.ConstEnsemble().GetSums(sums);

    // Get rates of inception processes.
    vector<Inception*>::const_iterator ii;
    for (ii=m_inceptions.begin(); ii!=m_inceptions.end(); ii++) {
        (*ii)->RateTerms(t, chem, T, P, sums, sys, iterm);
    }

    // Query other processes for their rates.
    vector<Process*>::const_iterator i;
    unsigned int j = 0;
    for(i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=rates.end()); i++) {
        if (!(*i)->IsDeferred()) {
            // Calculate rate if not deferred.
            (*i)->RateTerms(t, chem, T, P, sums, sys, iterm);
        } else {
            // If process is deferred, then set rate to zero.
            for (j=0; j<(*i)->TermCount(); j++) {*iterm=0.0; iterm++;}
        }
    }
    return 0;
}

real Mechanism::JumpRate(const std::vector<real> &rates) const
{
    // This function returns the combined rate for all non-deferred
    // processes.

    vector<Process*>::const_iterator i;
    vector<real>::const_iterator j=rates.begin();
    unsigned int k;
    real jrate = 0.0;

    // Sum up rates of inception processes.
    for (k=0; k<(unsigned int)m_inceptions.size(); k++,j++) {
        jrate += *j;
    }

    // Sum up rates for all processes which are not deferred.
    for (i=m_processes.begin(); (i!=m_processes.end()) && (j!=rates.end()); i++) {
        if (!(*i)->IsDeferred()) {
            for (k=0; k<(*i)->TermCount(); k++,j++) {
                jrate += *j;
            }
        }
    }

    return jrate;
}

void Mechanism::SetSpeciesList(SpeciesList &list)
{
    m_species = &list;
}

SpeciesList &Mechanism::GetSpeciesList(void) const
{
    return *m_species;
}

int Mechanism::AddInception(Inception *icn)
{
    // Add the inception to the mechanism.
    m_inceptions.push_back(icn);
    m_termcount += icn->TermCount();

    // Set the inception to belong to this mechanism.
    icn->SetMechanism(*const_cast<Mechanism*>(this));
    return 0;
}

int Mechanism::AddProcess(Process *p)
{
    // Add the process to the mechanism.
    m_processes.push_back(p);
    m_termcount += p->TermCount();

    // Check for any deferred.
    m_anydeferred = m_anydeferred || p->IsDeferred();

    // Set the process to belong to this mechanism.
    p->SetMechanism(*const_cast<Mechanism*>(this));
    return 0;
}

int Sweep::Mechanism::AddCoagulation()
{
    return AddProcess(new Coagulation());
}

void Mechanism::CheckDeferred(void)
{
    // Loop though all processes checking if any are deferred.
    m_anydeferred = false;
    vector<Process*>::iterator i;
    for (i=m_processes.begin(); i!=m_processes.end(); i++) {
        if ((*i)->IsDeferred()) {
            // Set anydeferred flag is true.
            m_anydeferred = true;
            return;
        }
    }
}

int Mechanism::DoProcess(const unsigned int i, const real t, System &sys) const
{
    // Work out to which process this term belongs.
    vector<Process*>::const_iterator ip;
    int j=(int)m_inceptions.size();

    if (i<(unsigned)j) {
        // This is an inception process.
        return m_inceptions[i]->Perform(t, sys, i);
    } else {
        // This is another process. 
        j -= 1;
        for(ip=m_processes.begin(); ip!=m_processes.end(); ip++) {
            j += (*ip)->TermCount();
            if ((unsigned)j>=i) {
                // Do the process.
                return (*ip)->Perform(t, sys, i+(*ip)->TermCount()-j-1);
            }
        }
    }

    // Not a valid process index.
    return -1;
}

int Mechanism::LPDA(const real t, System &sys)
{
    // Check that there are particles to update and that there are
    // deferred processes to perform.
    if (sys.ParticleCount() <= 0) return 0;
    if (!m_anydeferred) return 0;

    // Get conditions from the system.
    vector<real> chem, sums;
    real T, P;
    sys.GetConditions(t, chem, T, P);
    sys.Ensemble().GetSums(sums);

    // Stop ensemble from doubling while updating particles.
    sys.Ensemble().FreezeDoubling();

    // Perform deferred processes on all particles individually.
    Ensemble::iterator i;
    unsigned int k;
    int err;
    for (i=sys.Ensemble().begin(),k=0; i!=sys.Ensemble().end(); i++,k++) {
        err = UpdateParticle(*(*i), sys, t, chem, T, P, sums);
        if (err < 0) {
            return err;
        }
    }
    
    // Now remove any invalid particles and update the ensemble.
    sys.Ensemble().RemoveInvalids();
    sys.Ensemble().Update();

    // Start particle doubling again.  This will also double the ensemble
    // if too many particles have been removed.
    sys.Ensemble().UnfreezeDoubling();

    return 0;
}

int Mechanism::UpdateParticle(DefaultParticle &sp, System &sys, const real t)
{
    // If there are no deferred processes then stop right now.
    if (!m_anydeferred) return 0;

    // Get conditions from the system.
    vector<real> chem, sums;
    real T, P;
    sys.GetConditions(t, chem, T, P);
    sys.Ensemble().GetSums(sums);

    // Update the particle.
    return UpdateParticle(sp, sys, t, chem, T, P, sums);
}

int Mechanism::UpdateParticle(DefaultParticle &sp, System &sys, const real t,
                              const vector<real> &chem, const real T, const real P, 
                              const vector<real> &sums)
{
    // This function updates a single particle with deferred processes.  It is 
    // very important for LPDA.

    // If there are no deferred processes then stop right now.
    if (!m_anydeferred) return 0;

    vector<Process*>::const_iterator i;
    int err;
    unsigned int num;
    real rate, dt;

    while ((sp.Time() < t) && sp.IsValid()) {
        // Calculate delta-t and update particle time.
        dt = t - sp.Time();
        sp.SetTime(t);

        // Loop through all processes, performing those
        // which are deferred.
        for (i=m_processes.begin(); i!=m_processes.end(); i++) {
            if ((*i)->IsDeferred()) {
                // Get the process rate.
                rate = (*i)->Rate(t, chem, T, P, sums, sys, sp) * dt;

                // Use a Poission deviate to calculate number of times to perform
                // the process.
                num  = ignpoi(rate);

                if (num > 0) {
                    // Do the process to the particle.
                    err = (*i)->Perform(t, sys, sp, num);
                    if (err != 0) {
                        return err;
                    }
                }
            }
        }
    }

    // Check that the particle is still valid, only calculate cache if it is.
    if (sp.IsValid()) {
        sp.CalcCache();
        return 0;
    } else {
        return 1;
    }
}

Component const &Sweep::Mechanism::GetComponent(const unsigned int i) const
{
    // Check for valid index before returning component.  If index is invalid,
    // then just return the first component.
    if (i<(int)m_components.size()) {
        return *m_components[i];
    } else {
        return *m_components[0];
    }
}

int Sweep::Mechanism::GetComponentIndex(const std::string &name) const
{
    // Loop over all components checking the names.  When found the
    // desired component stop the loop.
    vector<Component*>::const_iterator i;
    unsigned int k;
    for (i=m_components.begin(),k=0; i!=m_components.end(); i++,k++) {
        if ((*i)->Name().compare(name)==0) {
            return k;
        }
    }

    // If the component name was not found then return an invalid index.
    return -1;
}

unsigned int Sweep::Mechanism::AddComponent(Sweep::Component &comp)
{
    m_components.push_back(&comp);
    return (unsigned int)m_components.size() - 1;
}

void Sweep::Mechanism::SetComponent(const unsigned int i, Sweep::Component &comp)
{
    m_components.reserve(i+1);
    m_components[i] = &comp;
}

void Sweep::Mechanism::SetComponents(const std::vector<Component*> &comps)
{
    m_components.assign(comps.begin(), comps.end());
}

int Sweep::Mechanism::GetValueIndex(const std::string &name) const
{
    // Loop over all value names checking the names.  When found the
    // desired value stop the loop.
    vector<string>::const_iterator i;
    unsigned int k;
    for (i=m_valuenames.begin(),k=0; i!=m_valuenames.end(); i++,k++) {
        if (i->compare(name)==0) {
            return k;
        }
    }

    // If the value name was not found then return an invalid index.
    return -1;
}

int Mechanism::InitReqdModels()
{
    set<Model>::const_iterator i;

    // Loop over all required models for the mechanism and initialise
    // them all.
    for (i=m_reqmodels.begin(); i!=m_reqmodels.end(); i++) {
        switch (*i) {
            case HACA :
                // Initialise the ABF HACA active sites model.
                ABF::ABFMech::InitHACA(*m_species);
            case PAH :
                // Initialise the PAH sites model.
                PAHModel::Initialise(*this);
//          case CNT :
                // Initialise the simple CNT model.  This isn't coded yet!
        }
    }

    return 0;
}