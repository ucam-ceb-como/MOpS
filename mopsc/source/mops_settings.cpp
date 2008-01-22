#include "mops_settings.h"
#include <vector>
#include <string>

using namespace Mops;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
Settings::Settings(void)
{
    m_atol = 1.0e-6;
    m_rtol = 1.0e-3;
    m_console_interval = 1;
    m_console_vars.clear();
    m_console_msgs = true;
    m_output_filename = "mops-out";
}

// Default destructor.
Settings::~Settings(void)
{
}


// ERROR TOLERANCES.

real Settings::ATOL() const
{
    return m_atol;
}

void Settings::SetATOL(real atol)
{
    m_atol = atol;
}

real Settings::RTOL() const
{
    return m_rtol;
}

void Settings::SetRTOL(real rtol)
{
    m_rtol = rtol;
}


// CONSOLE INTERVAL.

unsigned int Settings::ConsoleInterval() const
{
    return m_console_interval;
}

void Settings::SetConsoleInterval(unsigned int cint)
{
    m_console_interval = cint;
}


// CONSOLE VARIABLE NAMES.

const std::vector<std::string> &Settings::ConsoleVariables() const
{
    return m_console_vars;
}

const std::string &Settings::ConsoleVariable(unsigned int i) const
{
    if (i < m_console_vars.size()) {
        return m_console_vars[i];
    } else {
        // Returns the first variable name if index is invalid.     
        return m_console_vars[0];
    }
}

void Settings::AddConsoleVariable(const std::string &var)
{
    m_console_vars.push_back(var);
}

void Settings::RemoveConsoleVariable(const std::string &var)
{
    vector<string>::const_iterator i;
    for (i=m_console_vars.begin(); i!=m_console_vars.end(); i++) {
        if ((*i).compare(var) == 0) {
            m_console_vars.erase(i);
            return;
        }
    }
}

void Settings::RemoveConsoleVariable(unsigned int i)
{
    if (i < m_console_vars.size()) {
        m_console_vars.erase(m_console_vars.begin()+i);
    }
}


// CONSOLE MESSAGES.

bool Settings::UseConsoleMsgs() const
{
    return m_console_msgs;
}

void Settings::SetUseConsoleMsgs(bool msgs)
{
    m_console_msgs = msgs;
}


// OUTPUT FILE NAME.

const std::string &Settings::OutputFile() const
{
    return m_output_filename;
}

void Settings::SetOutputFile(const std::string &name)
{
    m_output_filename = name;
}
