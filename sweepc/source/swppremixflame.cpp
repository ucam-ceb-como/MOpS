#include "swppremixflame.h"
#include <string>
#include "morestrings.h"
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace Sweep;

PremixFlame::PremixFlame()
{
    m_iT = 0;
    m_iP = 1;
    m_profile.clear();
}

PremixFlame::~PremixFlame()
{
    m_profile.clear();
}

void PremixFlame::SetSpeciesList(SpeciesList &list)
{
    System::SetSpeciesList(list);
    m_iT = list.Count();
    m_iP = m_iT + 1;
};

real PremixFlame::GetTemperature(const real t) const 
{
    map<real,vector<real>>::const_iterator iter = FindTime(t);
    real s = iter->first;
    real T = iter->second[m_iT];
    if (++iter != m_profile.end()) {
        return T + ((iter->second[m_iT] - T) * (t - s)/ (iter->first - s));
    } else {
        return T;
    }
}

real PremixFlame::GetPressure(const real t) const 
{
    map<real,vector<real>>::const_iterator iter = FindTime(t);
    real s = iter->first;
    real P = iter->second[m_iP];
    if (++iter != m_profile.end()) {
        return P + ((iter->second[m_iP] - P) * (t - s)/ (iter->first - s));
    } else {
        return P;
    }
}

real PremixFlame::GetSpeciesConc(const unsigned int i, const real t) const 
{
    if (i<m_species->Count()) {
        map<real,vector<real>>::const_iterator iter = FindTime(t);
        real s = iter->first;
        real C = iter->second[i];
        if (++iter != m_profile.end()) {
            return C + ((iter->second[i] - C) * (t - s)/ (iter->first - s));
        } else {
            return C;
        }
    } else {
        return 0.0;
    }
}

real PremixFlame::GetSpeciesConc(const string name, const real t) const
{
    int i = m_species->GetIndex(name);
    if (i>=0) 
        return GetSpeciesConc(i, t);
    else
        return 0.0;
}

void PremixFlame::GetSpeciesConcs(const real t, vector<real> &chem) const 
{
    // Find nearest point before required and copy it to the output
    // vector.
    map<real,vector<real>>::const_iterator iter = FindTime(t);
    chem.assign(iter->second.begin(), iter->second.end());
    real s = iter->first;

    if (++iter != m_profile.end()) {
        // We are between two points, so we must use linear
        // interpolation.
        vector<real>::iterator ic;
        vector<real>::const_iterator ip;
        for (ic=chem.begin(),ip=iter->second.begin(); 
             ic!=chem.end(); 
             ic++,ip++) {
            *ic += ((*ip - *ic) * (t - s) / (iter->first - s));
        }
    }
}

void PremixFlame::GetConditions(const real time, vector<real> &chem, real &T, real &P) const 
{
    // Find nearest point before required and copy it to the output
    // vector.
    map<real,vector<real>>::const_iterator iter = FindTime(time);
    chem.assign(iter->second.begin(), iter->second.end());
    real s = iter->first;

    if (++iter != m_profile.end()) {
        // We are between two points, so we must use linear
        // interpolation.
        vector<real>::iterator ic;
        vector<real>::const_iterator ip;
        for (ic=chem.begin(),ip=iter->second.begin(); 
             ic!=chem.end(); 
             ic++,ip++) {
            *ic += ((*ip - *ic) * (time - s) / (iter->first - s));
        }
    }

    // Note temperature and pressure.
    T = chem[m_iT];
    P = chem[m_iP];
}

int PremixFlame::SetProfile(const std::vector<Sweep::real> &times, 
                            const std::vector<vector<Sweep::real>> &profile)
{
    vector<real>::size_type sz = profile[0].size();
    vector<real>::const_iterator titer;
    vector<vector<real>>::const_iterator piter;

    // Loop through the profile and check that all vectors are the same length.
    for(piter=profile.begin(); piter!=profile.end(); piter++) {
        if (piter->size() != sz) return -1;
    }

    // Clear the current profile.
    m_profile.clear();

    // Loop over times or profiles (whichever is shortest), and add each
    // entry into the profile.
    for(titer=times.begin(), piter=profile.begin();
        (titer!=times.end()) && (piter!=profile.end());
        titer++, piter++)
    {
            m_profile.insert(pair<real,vector<real>>(*titer,*piter));    
    }

    return 0;
}

int PremixFlame::AddTimePoint(const Sweep::real t, const std::vector<Sweep::real> conditions)
{
    // Check that time point is of correct size.
    if ((unsigned int)conditions.size() != (m_species->Count()+2)) {
        return -1;
    } else {
        // Insert new time point.
        m_profile.insert(pair<real,vector<real>>(t,conditions));
        return 0;
    }
}

int PremixFlame::ReadProfile(const string &filename)
{
    // Clear the current profile and species list.
    m_profile.clear();
    m_iT = 0;
    m_iP = 0;
    if (m_species != NULL) delete m_species;
    m_species = new SpeciesList();


    // Open the file to read.
    ifstream file(filename.c_str(), ios::in);
    if (!file.good()) return -1;

    // Variables read from file.
    vector<string> subs;
    string delim = ",\t ";
    string line;

    // Get the first line (should be the header defining the columns).
    if (!std::getline(file, line).eof()) {

         // Split the line to get the column headings.
        split(line, subs, delim);

        // Get important column indices (time, temperature and pressure).
        unsigned int tcol=0, Tcol=0, Pcol=0;
        tcol = findinlist(string("Time"), subs);
        Tcol = findinlist(string("T"), subs);
        Pcol = findinlist(string("P"), subs);

        // All other columns are chemical species.  Add them, and note
        // their columns.
        unsigned int i;
        map<unsigned int,unsigned int> spcols;
        for (i=0; i<(unsigned int)subs.size(); i++) {
            if ((i!=tcol) && (i!=Tcol) && (i!=Pcol)) {
                m_species->Add(subs[i]);
                spcols.insert(pair<unsigned int, unsigned int>(i, m_species->GetIndex(subs[i])));
            }
        }

        // Set up indices in profile and provide somewhere to read the
        // data.
        m_iT = m_species->Count();
        m_iP = m_iT + 1;
        real t;
        vector<real> prof(m_species->Count()+2, 0.0);
        map<unsigned int,unsigned int>::iterator isp;

        // Now we can read the profile.
        while(!std::getline(file, line).eof()) {
            split(line, subs, delim);
            for (i=0; i<(unsigned int)subs.size(); i++) {
                if (i==tcol) {
                    t = (real)strtod(subs[i].c_str(), 0);
                } else if (i==Tcol) {
                    prof[m_iT] = (real)strtod(subs[i].c_str(), 0);
                } else if (i==Pcol) {
                    prof[m_iP] = (real)strtod(subs[i].c_str(), 0);
                } else {
                    isp = spcols.find(i);
                    prof[isp->second] = (real)strtod(subs[i].c_str(), 0);
                }
            }
            m_profile.insert(pair<real,vector<real>>(t,prof));
        }

        file.close();
        return 0;
    } else {
        file.close();
        return -1;
    }
}

map<real,vector<real>>::const_iterator PremixFlame::FindTime(const real t) const 
{
    map<real,vector<real>>::const_iterator iter = m_profile.upper_bound(t);
    return --iter;
};

