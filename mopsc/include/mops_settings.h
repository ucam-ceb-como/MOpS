/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).

  File purpose:
    The Settings class holds simulation settings for mops.
*/

#ifndef MOPS_SETTINGS_H
#define MOPS_SETTINGS_H

#include "mops_params.h"
#include <vector>
#include <string>

namespace Mops
{
class Settings
{
public:
    // Constructors.
    Settings(void); // Default constructor.

    // Destructors.
    ~Settings(void); // Default destructor.

    // Get/Set error tolerances.
    real ATOL() const;
    void SetATOL(real atol);
    real RTOL() const;
    void SetRTOL(real rtol);

    // Get/set console interval.
    unsigned int ConsoleInterval() const;
    void SetConsoleInterval(unsigned int cint);

    // Get/Add/Remove console variable name(s).
    const std::vector<std::string> &ConsoleVariables() const;
    const std::string &ConsoleVariable(unsigned int i) const;
    void AddConsoleVariable(const std::string &var);
    void RemoveConsoleVariable(const std::string &var);
    void RemoveConsoleVariable(unsigned int i);

    // Get/set console messages flag.
    bool UseConsoleMsgs() const;
    void SetUseConsoleMsgs(bool msgs);

    // Get/set output file name.
    const std::string &OutputFile() const;
    void SetOutputFile(const std::string &name);

private:
    // Default error tolerances for the ODE solver.
    real m_atol, m_rtol;

    // Program output parameters:
    
    // Interval of console output data (in terms of time steps).
    unsigned int m_console_interval;
    
    // Console column variables.
    std::vector<std::string> m_console_vars;
    
    bool m_console_msgs; // Set to true if the console is to print messages.

    // Name of output file.
    std::string m_output_filename;
};
};

#endif
