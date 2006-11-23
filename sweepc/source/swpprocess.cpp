#include "swpprocess.h"
#include "swpmechanism.h"

using namespace Sweep;

Process::Process(void)
{
    m_mech = NULL;
    m_defer = false;
}

Process::~Process(void)
{
    m_reac.clear();
    m_prod.clear();
}

void Process::SetDeferred(const bool defer) 
{
    m_defer = defer; 
    if (m_mech!=NULL) m_mech->CheckDeferred();
}

void Process::SetReactant(const unsigned int id, const int stoich)
{
    map<unsigned int,int>::iterator i = m_reac.find(id);
    if (i == m_reac.end()) {
        // This reactant is not currently present, so add it.
        m_reac.insert(pair<unsigned int,int>(id,stoich));
    } else {
        // This reactant is already there, so set its value.
        (*i).second = stoich;
    }
}

void Process::SetProduct(const unsigned int id, const int stoich)
{
    map<unsigned int,int>::iterator i = m_prod.find(id);
    if (i == m_prod.end()) {
        // This product is not currently present, so add it.
        m_prod.insert(pair<unsigned int,int>(id,stoich));
    } else {
        // This product is already there, so set its value.
        (*i).second = stoich;
    }
}
