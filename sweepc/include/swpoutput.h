/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Output routines for Sweep.
*/

#ifndef SWEEP_OUTPUT_H
#define SWEEP_OUTPUT_H

#include "swpparams.h"
#include "swpsystem.h"
#include "swpmechanism.h"
#include "morestrings.h"

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

namespace Sweep
{
class SweepOutput
{
protected:
    static const int NSTATS = 10;
    static const int NPSL   = 4;
    static const real CONFA;
protected:
    ofstream m_file;
    string m_filename;
public:
    SweepOutput(void);
    ~SweepOutput(void);
public:
    inline bool IsOpen(void) const {return m_file.is_open();};
public:
    /* Compiles useful bulk and per-particle statistics about an ensemble. */
    static void CompileStatistics(const System &sys, const Mechanism &mech, vector<real> &stats);
    /* Compiles a particle size list for the given system, including all components and values. */
    static void CompilePSL(const System &sys, vector<vector<real>> &psl);
public:
    /* Opens an output stream to write binary data. The given
       number is appended to the file name as a way of allowing
       a series of files to be written. */
    int Open(const string &filename, const unsigned int n);
    /* Closes the output stream. */
    void Close();
    /* Writes current system to the output file. Returns 0 on success,
       otherwise negative. */
    int Write(const real t, const System &sys, const Mechanism &mech);
    /* Post processes a series of file containing a number of runs and
       writes the averages and errors to a CSV file. */
    int PostProcess(const string &filename, const unsigned int i1, const unsigned int i2, 
                    const unsigned int npoints, const Mechanism &mech) const;
public:
    /*  Compiles the Particle Size List from the given system and writes it
        to binary file.  The PSL is numbered with i and is written using the 
        filename stored in the output object. */
    int WritePSL(const real t, const System &sys, const unsigned int i);
    /*  Compiles the Particle Size List from the given system and writes it
        to binary file.  The PSL is numbered with i and is written using the 
        given filename. */
    static int WritePSL(const string &filename, const real t, const System &sys, const unsigned int i);
    /* Post-processes the particle size lists given their numbers. The PSLs are
       put consecutively into a single CSV file and the weights are recalculated
       based on the number of runs. */
    int PostProcessPSL(const string &filename, // Base file name for processing.
                       const unsigned int r1,  // ID of the first run.
                       const unsigned int r2,  // ID of last run (all runs in between are iterated over).
                       const vector<unsigned int> psln, // Numbers appended to PSL file names, also gives number of PSLS.
                       const Mechanism &mech); // Mechanism used to solve system.
    /* Post-processes a set of PSLs when there is only one PSL per run. */
    int PostProcessPSL(const string &filename, // Base file name for processing.
                       const unsigned int r1,  // ID of the first run.
                       const unsigned int r2,  // ID of last run (all runs in between are iterated over).
                       const unsigned int psln, // Number appended to PSL file names.
                       const Mechanism &mech); // Mechanism used to solve system.
public:
    /* Prints some data to the console. */
    static void PrintToConsole(const real t, const System &sys, const Mechanism &mech);
    /* Prints a string messageto the console. */
    inline static void PrintConsoleMsg(const char *const msg) {cout << msg << endl;};
public:
    /* Returns the column names in a post-processed output file. */
    static void GetColumnNames(const Mechanism &mech, vector<string> &names);
    /* Returns the column names in a post-processed PSL file. */
    static void GetPSLColumnNames(const Mechanism &mech, vector<string> &names);
protected:
    /* Returns the file name with the number appended in a standard way. */
    static string BuildFilename(const string &basename, const unsigned int i);
    static string BuildPSLFilename(const string &basename, const unsigned int i);
};

inline string SweepOutput::BuildFilename(const string &basename, const unsigned int i)
{
    string name;
    string::size_type iext = basename.find_last_of(".",basename.length());
    name = basename.substr(0, iext);
    name.append("(");
    name.append(cstr<const unsigned int>(i));
    name.append(")");
    name.append(basename.substr(iext, basename.length()-iext));
    return name;
}

inline string SweepOutput::BuildPSLFilename(const string &basename, const unsigned int i)
{
    string name;
    string::size_type iext = basename.find_last_of(".",basename.length());
    name = basename.substr(0, iext);
    name.append("-psl(");
    name.append(cstr<const unsigned int>(i));
    name.append(")");
    name.append(basename.substr(iext, basename.length()-iext));
    return name;
}

inline void SweepOutput::GetColumnNames(const Sweep::Mechanism &mech, std::vector<string> &names)
{
    names.resize(NSTATS+mech.ComponentCount()+mech.ValueCount());
    unsigned int i, j;

    // Stat names.
    names[0] = "Stoch. Particle Count";
    names[1] = "Number Density (cm-3)";
    names[2] = "Mass Concentration (g/cm3)";
    names[3] = "Volume Fraction";
    names[4] = "Surface Area (cm2/cm3)";
    names[5] = "Average Mass (g)";
    names[6] = "Average Volume (cm3)";
    names[7] = "Average Surface Area (cm2)";
    names[8] = "Average Diameter (nm)";
    names[9] = "Sample Volume (cm-3)";
    // Component names.
    for (i=NSTATS,j=0; j<mech.ComponentCount(); i++,j++) {
        names[i] = mech.GetComponent(j).Name();
    }
    // value names.
    for(j=0; j<mech.ValueCount(); i++,j++) {
        names[i] = mech.GetValueName(j);
    }
}

inline void SweepOutput::GetPSLColumnNames(const Sweep::Mechanism &mech, std::vector<string> &names)
{
    names.resize(NPSL+mech.ComponentCount()+mech.ValueCount());
    unsigned int i, j;

    // Stat names.
    names[0] = "Weight";
    names[1] = "Volume";
    names[2] = "Surface Area";
    names[3] = "Diameter";
    // Component names.
    for (i=NPSL,j=0; j<mech.ComponentCount(); i++,j++) {
        names[i] = mech.GetComponent(j).Name();
    }
    // Value names.
    for(j=0; j<mech.ValueCount(); i++,j++) {
        names[i] = mech.GetValueName(j);
    }
}

};

#endif