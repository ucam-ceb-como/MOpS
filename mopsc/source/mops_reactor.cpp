#include "mops_reactor.h"
#include "mops_mixture.h"

// CVODE includes.
#include "cvode\cvode.h"
#include "cvode\cvode_dense.h"
#include "nvector\nvector_serial.h"

#include <vector>
#include <math.h>

using namespace Mops;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

Reactor::Reactor(void)
{
    m_time = 0.0;
    m_mix  = NULL;
    m_mech = NULL;
    m_deriv = NULL;
    m_atol = 1.0e-6;
    m_rtol = 1.0e-3;
    m_neq   = 0;
    m_iT    = -1;
    m_iDens = -1;

    // Init CVODE.
    m_odewk = CVodeCreate(CV_BDF, CV_NEWTON);
}

Reactor::~Reactor(void)
{
    CVodeFree(&m_odewk);
    if (m_deriv != NULL) delete [] m_deriv;
}


// REACTOR SOLUTION.

// Returns the current reactor time.
real Reactor::Time() const
{
    return m_time;
}

// Initialises the reactor to the given time.
void Reactor::Initialise(real time)
{
    m_time = time;
    
    // Fill solution and derivative vectors.
//    buildSoln();
    for (int i=0; i<m_neq; i++) {
        m_deriv[i] = 0.0;
    }

    // Allocate CVODE stuff.
    CVodeMalloc(m_odewk, &rhsFn, time, N_VMake_Serial(m_neq, m_mix->RawData()), 
                CV_SS, m_rtol, (void*)&m_atol);

    // Set up internal solution array.
    m_solvec = N_VNewEmpty_Serial(m_neq);

    // Set the f_data pointer to this Reactor object.
    CVodeSetFdata(m_odewk, (void*)this);

    // Set other parameters.
    CVodeSetMaxNumSteps(m_odewk, 2000);

    // Set CVDense as the linear system solver.
    CVDense(m_odewk, m_neq);
}

// Solves the reactor up to the given time.
void Reactor::Solve(real time)
{
    // Put the solution into an N_Vector data pointer.  This does not
    // involve copying the data.
    NV_DATA_S(m_solvec) = m_mix->RawData();
    CVodeSetStopTime(m_odewk, time);

    // Solve over time step.
    while (m_time < time) {
        CVode(m_odewk, time, m_solvec, &m_time, CV_NORMAL);

        m_mix->Normalise();
    }

    // Calculate derivatives at end point.
    if (m_emodel == ConstT) {
        RHS_ConstT(time, m_mix->RawData(), m_deriv);
    } else {
        RHS_Adiabatic(time, m_mix->RawData(), m_deriv);
    }
}


// REACTOR CONTENTS.

// Returns a pointer to the current reactor contents.
Mops::Mixture *const Reactor::Mixture() const
{
    return m_mix;
}

// Fills the reactor with the given mixture.  If the clearfirst flag
// is set then the current mixture is first deleted from memory.
void Reactor::Fill(Mops::Mixture *const mix, bool clearfirst)
{
    if (clearfirst && (m_mix!=NULL)) {
        // First delete current mixture from memory.
        delete m_mix;
        m_mix = NULL;
    }

    // Set new mixture.
    m_mix = mix;

    // The number of equations to solve is the number of chemical species
    // plus the temperature and the density.
    m_neq   = m_mix->Species()->size() + 2;
    m_iT    = m_neq - 2; // Second to last element is temperature.
    m_iDens = m_neq - 1; // Last element is density.

    // Allocate the derivative array.
    if (m_deriv != NULL) delete [] m_deriv;
    m_deriv = new real[m_neq];
}


// REACTOR MECHANISM.

// Returns the current mechanism.
const Mops::Mechanism *const Reactor::Mechanism() const
{
    return m_mech;
}

// Sets the reactor mechanism.
void Reactor::SetMechanism(const Mops::Mechanism *const mech)
{
    m_mech = mech;
}


// ENERGY MODEL.

// Returns the current energy model.
Reactor::EnergyModel Reactor::EnergyEquation() const
{
    return m_emodel;
}

// Sets the energy model.
void Reactor::SetEnergyEquation(Reactor::EnergyModel model)
{
    m_emodel = model;
}


// RHS FUNCTION.

// The right-hand side evaluator.  This function calculates the RHS of
// the reactor differential equations.  CVODE uses a void* pointer to
// allow the calling code to pass whatever information it wants to
// the RHS function.  In this case the void* pointer should be cast
// into a Reactor object.
int Reactor::rhsFn(double t,      // Independent variable.
                   N_Vector y,    // Solution array.
                   N_Vector ydot, // Derivatives of y wrt t.
                   void* reactor) // Pointer to reactor object.
{
    // Cast the Reactor object.
    Reactor *r = static_cast<Reactor*>(reactor);

    /*
    // Check values of y.
    real ytot = 0.0;
    for (int i=0; i<r->m_neq-2; i++) {
        ytot += NV_DATA_S(y)[i];
        if (NV_DATA_S(y)[i] < -r->m_atol) {
            return 1;
        }
        if (NV_DATA_S(y)[i] > 1.0+r->m_atol) {
            return 1;
        }
    }
    if (ytot > 1.0+r->m_atol) {
        return 2;
    }
    if (NV_DATA_S(y)[r->m_neq-2] < -r->m_atol) {
        return 3;
    }
    if (NV_DATA_S(y)[r->m_neq-1] < -r->m_atol) {
        return 4;
    }
    */

    // Get the RHS from the system model.
    if (r->m_emodel == ConstT) {
        r->RHS_ConstT(t, NV_DATA_S(y), NV_DATA_S(ydot));
    } else {
        r->RHS_Adiabatic(t, NV_DATA_S(y), NV_DATA_S(ydot));
    }

    return 0;
};

// Definition of RHS form for constant temperature energy equation.
void Reactor::RHS_ConstT(real t, const real *const y,  real *ydot)
{
    int i;
    vector<real> wdot, Gs, rop, kf, kr;
    real wtot = 0.0, ydot0;

    m_mix->CalcGs_RT(y[m_iT], Gs);
    m_mech->Reactions().GetRateConstants(y[m_iT], y[m_iDens], y, 
                                         m_mech->Species().size(), Gs, kf, kr);
    m_mech->Reactions().GetRatesOfProgress(y[m_iDens], y, m_mech->Species().size(), 
                                           kf, kr, rop);
    m_mech->Reactions().GetMolarProdRates(rop, wdot);

    // Calculate mole fraction derivatives.
    for (i=0; i<m_neq-2; i++) {
        wtot += wdot[i];
    }
    //ydot[0];
    ydot0 = 0.0;
    for (i=0; i<m_neq-2; i++) {
        ydot[i] = (wdot[i] - (y[i]*wtot)) / y[m_iDens];
        //ydot[0] -= ydot[i];
        if (i>0) ydot0 -= ydot[i];
    }

    ydot[m_neq-2] = 0.0; // Temperature derivative.
    ydot[m_neq-1] = 0.0; // Density derivative.
}

// Definition of RHS form for adiabatic energy equation.
void Reactor::RHS_Adiabatic(real t, const real *const y,  real *ydot)
{
    int i;
    vector<real> wdot, Gs, rop, kf, kr;
    real wtot = 0.0, ydot0;

    m_mix->CalcGs_RT(y[m_iT], Gs);
    m_mech->Reactions().GetRateConstants(y[m_iT], y[m_iDens], y, 
                                         m_mech->Species().size(), Gs, kf, kr);
    m_mech->Reactions().GetRatesOfProgress(y[m_iDens], y, m_mech->Species().size(), 
                                           kf, kr, rop);
    m_mech->Reactions().GetMolarProdRates(rop, wdot);

    // Calculate mole fraction derivatives.
    for (i=0; i<m_neq-2; i++) {
        wtot += wdot[i];
    }
    //ydot[0];
    ydot0 = 0.0;
    for (i=0; i<m_neq-2; i++) {
        ydot[i] = (wdot[i] - (y[i]*wtot)) / y[m_iDens];
        //ydot[0] -= ydot[i];
        if (i>0) ydot0 -= ydot[i];
    }

    // TODO:  Include adiabatic energy equation here.

    ydot[m_neq-2] = 0.0; // Temperature derivative.
    ydot[m_neq-1] = 0.0; // Density derivative.
}
