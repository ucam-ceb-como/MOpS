/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver) - later sprogc.

  File purpose:
    This class holds a list of chemical species.  At the moment it only
    holds their names, but in future the class should be made part of the 
    c++ rewrite of the F90 sprog code.  sprog is a competitor to CHEMKIN and
    is designed to be less bloated than CHEMKIN or cantera, therefore easier
    to use.

    The class is designed to not store a list of "species" objects but rather
    just store arrays of species information.  This decision was taken because
    the majority of the time codes seem to use these arrays, rather than access
    information species-wise.  Therefore, it is faster to access the data in this
    format.
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
    vector<string> m_names; // List of species names.
public: // Default constructor and destructor.
    SpeciesList(void);
    ~SpeciesList(void);
public:
    unsigned int Count(void) const;
public: // Routines for getting and setting species names.
    /* Returns the name of the species with index i in the list. */
    const string  GetName(unsigned int i) const;
    /* Returns the vector of all species names in the list. */
    const vector<string> &GetNames(void) const;
    /* Set the names of all species in the list. */
    void SetNames(const vector<string> &names);
    /* Adds a species to the list given only the name. */
    unsigned int Add(const string &name);
    /* Removes the species with the given name from the list. */
    void Remove(const string &name);
    /* Removes the species at the given index from the list. */
    void Remove(const unsigned int i);
public: // Routines for getting species' indices.
    /* Returns the index of the species (given by name) in the
       list on success, otherwise return less than zero. */
    const int GetIndex(const string &name) const;
};

#endif