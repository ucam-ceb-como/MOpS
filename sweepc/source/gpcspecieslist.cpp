#include "gpcspecieslist.h"

using namespace std;

SpeciesList::SpeciesList(void)
{
}

SpeciesList::~SpeciesList(void)
{
}

unsigned int SpeciesList::Count() const
{
    return (unsigned int)m_names.size();
}

string &SpeciesList::GetName(const unsigned int i)
{
    return m_names[i];
}

vector<string> &SpeciesList::GetNames()
{
    return m_names;
}

void SpeciesList::SetNames(const std::vector<string> &names)
{
    m_names.assign(names.begin(), names.end());
}

unsigned int SpeciesList::Add(const std::string &name)
{
    unsigned int i;
    for(i=0; i<(unsigned int)m_names.size(); i++) {
        if (m_names[i].compare(name) == 0) return i;
    }

    m_names.push_back(name);
    return (unsigned int)m_names.size() - 1;
}

void SpeciesList::Remove(const std::string &name)
{
    vector<string>::iterator i = m_names.begin();
    while((i!=m_names.end()) && ((*i).compare(name) != 0)) i++;
    if (i!=m_names.end()) m_names.erase(i);
}

void SpeciesList::Remove(const unsigned int i)
{
    if (i < (unsigned int)m_names.size()) {
        m_names.erase(m_names.begin()+i);
    }
}

int SpeciesList::GetIndex(const string &name) const
{
    int i;
    for (i=0; i<(int)m_names.size(); i++) {
        if (m_names[i].compare(name) == 0) return i;
    }
    return -1;
}