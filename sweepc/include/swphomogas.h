/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    A homogeneous gas system with a coupling between the
    gas and the particles.
*/

#ifndef SWEEP_HOMOGAS_H
#define SWEEP_HOMOGAS_H

#include <vector>
#include "swpparams.h"
#include "swpsystem.h"

namespace Sweep
{
class HomoGas : public System
{
protected:
    vector<real> m_initchem; // Initial gas-phase species concentrations.
    real m_initT, m_initP; //
    real m_inittime; // Time for which the initial conditions are valid.
    vector<real> m_chem; // Gas-phase species concentrations.
    real m_t, m_p; // Gas-phase temperature and pressure.
public:
    HomoGas(void);
    ~HomoGas(void);
public: // Definition of pure virtuals.
    inline real GetTemperature(const real t) const {return m_t;};
    inline real GetPressure(const real t) const {return m_p;};
    inline real GetSpeciesConc(const unsigned int i, const real t) const {
        if (i<(int)m_chem.size()) {
            return m_chem[i];
        } else {
            return 0.0;
        }
    };
    inline void SetTemperature(const real t) {m_t = t;};
    inline void SetPressure(const real p) {m_p = p;};
    inline void SetSpeciesConc(const unsigned int i, const real c) {
        if (i<m_chem.size()) m_chem[i] = c;
    }
    /* Adds the given value to the required species concentration. */
    inline void AdjustSpeciesConc(const unsigned int i, const real delta) {
        if (i<m_chem.size()) {
            m_chem[i] += delta;
            if (m_chem[i] < 0.0) m_chem[i] = 0.0;
        }
    }
    /* Sets the initial conditions and sets the current conditions to those. */
    void SetInitConditions(const real time, const vector<real> &chem, const real T, const real P);
};
};

#endif