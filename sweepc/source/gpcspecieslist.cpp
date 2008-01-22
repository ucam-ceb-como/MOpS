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
    // Simply return the size of the names vector.
    return (unsigned int)m_names.size();
}

const string SpeciesList::GetName(const unsigned int i) const
{
    // If the index is within range then return the name, otherwise
    // return a blank string.
    if (i<(unsigned int)m_names.size()) {
        return m_names[i];
    } else {
        return "";
    }
}

const vector<string> &SpeciesList::GetNames() const
{
    // Return the names vector.
    return m_names;
}

void SpeciesList::SetNames(const std::vector<string> &names)
{
    // Overwrite the current list with the new names.
    m_names.assign(names.begin(), names.end());
}

unsigned int SpeciesList::Add(const std::string &name)
{
    // First check the list to see if the species is already
    // in it.  If it is, then return its index in the list.
    unsigned int i;
    for(i=0; i<(unsigned int)m_names.size(); i++) {
        if (m_names[i].compare(name) == 0) return i;
    }

    // The species is not in the list, so add it to the end and
    // return the index.
    m_names.push_back(name);
    return (unsigned int)m_names.size() - 1;
}

void SpeciesList::Remove(const std::string &name)
{
    // Locate the species with this name in the list.
    vector<string>::iterator i = m_names.begin();
    while((i!=m_names.end()) && ((*i).compare(name) != 0)) i++;

    // If the iterator is valid (species is in list) then erase it.
    if (i!=m_names.end()) m_names.erase(i);
}

void SpeciesList::Remove(const unsigned int i)
{
    // Erase the species name as long as the index is valid, i.e.
    // within range.
    if (i < (unsigned int)m_names.size()) {
        m_names.erase(m_names.begin()+i);
    }
}

const int SpeciesList::GetIndex(const string &name) const
{
    int i;

    // Loop through all species names by index until we find one
    // that matches the given name.
    for (i=0; i<(int)m_names.size(); i++) {
        if (m_names[i].compare(name) == 0) return i;
    }
    return -1;
}