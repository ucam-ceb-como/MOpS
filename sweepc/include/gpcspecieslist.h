/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    
*/

#ifndef SPECIESLIST_H
#define SPECIESLIST_H

#include <string>
#include <vector>
#include <set>

using namespace std;

class SpeciesList
{
private:
    vector<string> m_names;
public:
    SpeciesList(void);
    ~SpeciesList(void);
public:
    unsigned int Count(void) const;
public: // Property gets/sets.
    string &GetName(unsigned int i);
    vector<string> &GetNames(void);
    void SetNames(const vector<string> &names);
public: // Add and remove routines.
    unsigned int Add(const string &name);
    void Remove(const string &name);
    void Remove(const unsigned int i);
public: // By index.
    /* Returns the index of the species (given by name) in the
       list on success, otherwise return less than zero. */
    int GetIndex(const string &name) const;
};

#endif