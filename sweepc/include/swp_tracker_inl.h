/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of inline member functions of the Tracker class.
*/

#ifndef SWEEP_TRACKER_INL_H
#define SWEEP_TRACKER_INL_H

#include "swp_params.h"
#include "swp_tracker.h"
#include <string>

// TRACKER NAME.

// Returns component symbol or name.
inline const std::string &Sweep::Tracker::Name() const {return m_name;};

// Sets the symbol or name.
inline void Sweep::Tracker::SetName(const std::string &name) {m_name = name;};

#endif
