#include "swpmechanism.h"
#include "swpcoagulation.h"
#include "swpensemble.h"
#include "rng.h"

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
    for(i=m_processes.begin(); (i!=m_processes.end()) && (iterm!=rates.end()); i++) {
        (*i)->RateTerms(t, chem, T, P, sums, sys, iterm);
    }

    return 0;
}

int Mechanism::GetStochasticRates(std::vector<real> &rates, const real t, 
                                  const System &sys) const
{
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
    m_inceptions.push_back(icn);
    icn->SetMechanism(*const_cast<Mechanism*>(this));
    m_termcount += icn->TermCount();
    return 0;
}

int Mechanism::AddProcess(Process *p)
{
    m_processes.push_back(p);
    m_anydeferred = m_anydeferred || p->IsDeferred();
    p->SetMechanism(*const_cast<Mechanism*>(this));
    m_termcount += p->TermCount();
    return 0;
}

int Sweep::Mechanism::AddCoagulation()
{
    return AddProcess(new Coagulation());
}

void Mechanism::CheckDeferred(void)
{
    m_anydeferred = false;
    vector<Process*>::iterator i;
    for (i=m_processes.begin(); i!=m_processes.end(); i++) {
        if ((*i)->IsDeferred()) m_anydeferred = true;
        return;
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
    
    // Now remove any invalid particles.
    sys.Ensemble().RemoveInvalids();
    sys.Ensemble().Update();

    // Start particle doubling again.  This will also double the ensemble
    // if any particles have been removed.
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

    return UpdateParticle(sp, sys, t, chem, T, P, sums);
}

int Mechanism::UpdateParticle(DefaultParticle &sp, System &sys, const real t,
                              const vector<real> &chem, const real T, const real P, 
                              const vector<real> &sums)
{
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
                rate = (*i)->Rate(t, chem, T, P, sums, sys, sp) * dt;
                num  = ignpoi(rate);
                if (num > 0) {
                    err = (*i)->Perform(t, sys, sp, num);
                    if (err != 0) {
                        return err;
                    }
                }
            }
        }
    }

    if (sp.IsValid()) {
        sp.CalcCache();
        return 0;
    } else {
        return 1;
    }
}

Component const &Sweep::Mechanism::GetComponent(const unsigned int i) const
{
    if (i<(int)m_components.size()) {
        return *m_components[i];
    } else {
        return *m_components[0];
    }
}

int Sweep::Mechanism::GetComponentIndex(const std::string &name) const
{
    vector<Component*>::const_iterator i;
    unsigned int k;
    for (i=m_components.begin(),k=0; i!=m_components.end(); i++,k++) {
        if ((*i)->Name().compare(name)==0) {
            return k;
        }
    }
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
    vector<string>::const_iterator i;
    unsigned int k;
    for (i=m_valuenames.begin(),k=0; i!=m_valuenames.end(); i++,k++) {
        if (i->compare(name)==0) {
            return k;
        }
    }
    return -1;
}
