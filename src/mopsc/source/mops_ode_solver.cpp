/*
  Author(s):      Matthew Celnik (msc37), Weerapong Phadungsukanan (wp214)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ODE_Solver class declared in the
    mops_ode_solver.h header file.

  Future work:
    Copy constructure is not fully copy CVODES ODE workspace. This
    is need to be finished.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#include "mops_ode_solver.h"
#include "mops_reactor.h"
#include "mops_rhs_func.h"
#include "cvodes_utils.h"

// CVODE includes.
#include "cvodes/cvodes.h"
#include "cvodes/cvodes_dense.h"

#include <vector>
#include <cmath>
#include <stdexcept>

using namespace Mops;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ODE_Solver::ODE_Solver(void)
{
    init();
}

// Copy constructor.
ODE_Solver::ODE_Solver(const Mops::ODE_Solver &copy)
{
    init();
    // Use assignment operator.
    *this = copy;
}

// Stream-reading constructor.
ODE_Solver::ODE_Solver(std::istream &in)
{
    Deserialize(in);
}

// Default destructor.
ODE_Solver::~ODE_Solver(void)
{
    releaseMemory();
}


// OPERATORS.

// Assignment operator. (This function cannot be used with sensitivity problem.)
ODE_Solver &ODE_Solver::operator=(const Mops::ODE_Solver &rhs)
{
    if (this != &rhs) {
        m_time     = rhs.m_time;
        m_soln     = rhs.m_soln;
        m_reactor  = rhs.m_reactor;
        m_rtol     = rhs.m_rtol;
        m_atol     = rhs.m_atol;
        m_neq      = rhs.m_neq;
        m_srcterms = rhs.m_srcterms;
        _srcTerms  = rhs._srcTerms;

        // Delete memories
        if (m_yS != NULL) N_VDestroyVectorArray_Serial(m_yS, m_sensi.NParams());
        // Reset pointers to NULL in case of rsh is zero.
        m_yS = NULL;
        m_sensi = rhs.m_sensi;
        // No memory allocation require if m_NS is zero.
        if (m_sensi.NParams() != 0) {
            m_yS = N_VCloneVectorArrayEmpty_Serial(m_sensi.NParams(), rhs.m_yS[0]);
        }
        // Set values of varibles to rhs variables.
        for (unsigned int i = 0; i < m_sensi.NParams(); i++) {
            m_yS[i] = CVODES::N_VExactClone_Serial(rhs.m_yS[i]);
        }
        if (m_yvec != NULL) N_VDestroy_Serial(m_yvec);
        m_yvec = CVODES::N_VExactClone_Serial(rhs.m_yvec);
        //if (m_solvec != NULL) N_VDestroy_Serial(m_solvec);
        //VCopy_Serial(rhs.m_solvec, m_solvec);

        // Copy ODE workspace.
        if (rhs.m_odewk != NULL) {
            if (m_odewk == NULL) InitCVode();
            CVODES::CVodeMemRecCopy_Serial(*((CVodeMem)m_odewk), *((CVodeMem)rhs.m_odewk));
        } else {
            if (m_odewk != NULL) CVodeFree(&m_odewk);
            m_odewk = NULL;
        }

        // Copy ODE workspace, incl. derivatives.
        memcpy(m_deriv, rhs.m_deriv, sizeof(double)*m_neq);
    }
    return *this;
}

// SOLVER SET UP.

// Initialises the reactor to the given time.
void ODE_Solver::Initialise(Reactor &reac)
{
    // Store useful variables.
    m_time = reac.Time();
    m_neq  = reac.ODE_Count();
    
    // Store solution vector.
    m_soln = reac.Mixture()->GasPhase().RawData();

    // Fill derivative vector.
    delete [] m_deriv;
    m_deriv = new double[m_neq];
    for (unsigned int i=0; i!=m_neq; ++i) {
        m_deriv[i] = 0.0;
    }

    InitCVode();    
}

// Initialises the CVode ODE solver assuming that the
// remainder of the the solver has been correctly set up.
void ODE_Solver::InitCVode(void)
{
    // Create ODE workspace.
    if (m_odewk != NULL) CVodeFree(&m_odewk);
    m_odewk = CVodeCreate(CV_BDF, CV_NEWTON);

    // Allocate CVODE stuff.
    if (m_yvec != NULL) N_VDestroy_Serial(m_yvec);
    m_yvec = N_VMake_Serial(m_neq, m_soln);

    if (m_sensi.isEnable()) {
        CVodeInit(m_odewk, &rhsFn_CVODES, m_time, m_yvec);
    } else {
        CVodeInit(m_odewk, &rhsFn_CVODE, m_time, m_yvec);
    }

    CVodeSStolerances(m_odewk, m_rtol, m_atol);

    // Set up internal solution array.
    m_solvec = N_VNewEmpty_Serial(m_neq);
    NV_DATA_S(m_solvec) = m_soln;

    // Set the f_data pointer to this Reactor object.
    CVodeSetUserData(m_odewk, (void*)this);

    // Set other parameters.
    CVodeSetMaxNumSteps(m_odewk, 2000);

    // Set CVDense as the linear system solver.
    CVDense(m_odewk, m_neq);

    // Set the Jacobian function in CVODE.
    // Benchmark results (CVODES 2.1) from certain test case :
    // - No sensitivity analysis : External Jacobian is faster than CVODE internal jacobain.
    // - Sensitivity analyis (Rate parameters) : CVODE internal jacobian is fastest (114 s). jacFn_CVODES is slightly slower (120 s)
    //   and jacFn_CVODE is very slow (146 s). Thus, Jacobian function will be set according to this test.
    if (m_sensi.isEnable()) {
        if (m_sensi.ProblemType() == SensitivityAnalyzer::Reaction_Rates) {
            // Internal one is fastest so don't set jacobian function.
            // CVDenseSetJacFn(m_odewk, jacFn_CVODES, (void*)this);
        } else if (m_sensi.ProblemType() == SensitivityAnalyzer::Init_Conditions) {
            // CVDenseSetJacFn(m_odewk, jacFn_CVODES, (void*)this); //<== need test
        }
    } else {
        // CVDlsSetDenseJacFn(m_odewk, &jacFn_CVODE);
    }
    int flag = 0;
    if (m_sensi.isEnable()) {
        m_yS = N_VCloneVectorArray_Serial(m_sensi.NParams(), m_yvec);
        m_sensi.InitSensMatrix(m_yS);
        flag = CVodeSensInit(m_odewk, m_sensi.NParams(), m_sensi.GetMethod(), NULL, m_yS);

        flag = CVodeSensEEtolerances(m_odewk);
        //flag = CVodeSensSStolerances(m_odewk, m_rtol, &m_atol);
        flag = CVodeSetSensErrCon(m_odewk, m_sensi.isEnableErrorControl());

        flag = CVodeSetSensParams(m_odewk, m_sensi.ParamsPtr(), m_sensi.ParamBarsPtr(), NULL);

    } else {
        //printf("Sensitivity: NO ");
    }
}

// Reset the solver.  Need to do this if the the reactor
// contents has been changed between calls to Solve().
void ODE_Solver::ResetSolver(void)
{
    if (m_yvec != NULL) N_VDestroy_Serial(m_yvec);
    m_yvec = N_VMake_Serial(m_neq, m_soln);
    // m_yS cannot be reset since it need to know the previous values
    // in order to continue solving the next step.
    if (m_sensi.isEnable()) {
        CVodeReInit(m_odewk, m_time, 
                    m_yvec);
        CVodeSensReInit(m_odewk, m_sensi.GetMethod(), m_yS);
    } else {
        CVodeReInit(m_odewk, m_time, 
                    m_yvec);
    }
    
} 

// Reset the solver.  Need to do this if the the reactor
// contents has been changed between calls to Solve().
void ODE_Solver::ResetSolver(Reactor &reac)
{
    // Check that this reactor has the same problem size
    // as the last reactor.  If not then we have to resize the
    // workspace, which is done by the Initialise() routine.
    if (reac.ODE_Count() == m_neq) {
        m_time = reac.Time();
        m_soln = reac.Mixture()->GasPhase().RawData();
        if (m_yvec != NULL) N_VDestroy_Serial(m_yvec);
        m_yvec = N_VMake_Serial(m_neq, m_soln);
        // m_yS cannot be reset since it need to know the previous values.
        if (m_sensi.isEnable()) {
            CVodeReInit(m_odewk, m_time, m_yvec);
            CVodeSensReInit(m_odewk, m_sensi.GetMethod(), m_yS);
        } else {
            CVodeReInit(m_odewk, m_time, m_yvec);
        }
    } else {
        Initialise(reac);
    }
}


// SOLVING THE REACTOR.

// Solves the reactor up to the given time.
void ODE_Solver::Solve(Reactor &reac, double stop_time)
{
    // Check that the reactor has the same problem size
    // as set in the solver.
    if (m_neq != reac.ODE_Count()) {
        // Re-initialise the solver for a new problem size.
        Initialise(reac);
    } else {
        // Get the time and solution pointer from the reactor.
        m_time = reac.Time();
        m_soln = reac.Mixture()->GasPhase().RawData();
    }

    // Put the solution into an N_Vector data pointer.  This does not
    // involve copying the data.
    m_reactor = &reac;
    NV_DATA_S(m_solvec) = m_soln;
    CVodeSetStopTime(m_odewk, stop_time);

    // Solve over time step.
    while (m_time < stop_time) {
        int CVode_error = 0;
        CVode_error = CVode(m_odewk, stop_time, m_solvec, &m_time, CV_NORMAL);
        // If no error then this line will quickly skip throwing exception which help to speed up checking
        if (CVode_error < 0) {
            switch (CVode_error) {
                case -4 :
                    throw runtime_error("Convergence test failures occurred too many time or occurred with |h| = hmin (Mops, ODE_Solver::Solve).");
                    break;
                default :
                    throw runtime_error("CVode returns error, please find out the error code and add to report error (Mops, ODE_Solver::Solve).");
            }
            //throw invalid_argument(" (Mops, ODE_Solver::Solve).");
        }
        reac.Mixture()->GasPhase().Normalise(); // This should not be required if CVODE solves correctly.
    }
    

    // Calculate derivatives at end point.
    if (reac.EnergyEquation() == Reactor::ConstT) {
        reac.RHS_ConstT(stop_time, m_soln, m_deriv);
    } else {
        reac.RHS_Adiabatic(stop_time, m_soln, m_deriv);
    }
    // Add the source terms to derivatives, if defined.
    if (_srcTerms != NULL) 
        _srcTerms(m_deriv, m_neq, stop_time, *m_srcterms);

    // Set a pointer in sensitivity object to result for outputting.
    if (m_sensi.isEnable()) {
        CVodeGetSens(m_odewk, &m_time, m_yS);
        m_sensi.SetSensResult(m_yS);
    }

}

/*!
Method creates 2D array for sensitivity values and initialises it to zero
@param[in]      n_species        number of species
@param[in]      n_sensi          number of sensitivities computed
*/
void ODE_Solver::InitialiseSensArray(int n_sensi, int n_species)
{
    m_sensitivity = new double*[n_sensi];
    for (int i = 0; i < n_sensi; i++){
        m_sensitivity[i] = new double[n_species];
        for (int j = 0; j < n_species; j++){
            m_sensitivity[i][j] = 1.0;
        }
    }
}

/*!
Method deletes the memory allocated to double** m_sensitivity
@param[in]      n_species        number of species in sensitivity matrix
*/
void ODE_Solver::DestroySensArray(int n_species)
{
    for (int i = 0; i < n_species; i++){
        delete [] m_sensitivity[i];
    }
    delete m_sensitivity;
}

/*!
Method loads sensitivity values from N_Vector* into it
@param[in]      n_species        number of species
@param[in]      n_sensi          number of sensitivities computed
@return         m_sensitivity    Sensitivity Array
*/
double** ODE_Solver::GetSensSolution(int n_sensi, int n_species)
{
    for (int i = 0; i < (n_sensi); i++){
        //double* test; 
		double* test = new double[n_sensi];
        test[i] = 0.0;
        test = NV_DATA_S(m_yS[i]);
        for (int j = 0; j < (n_species); j++){
            m_sensitivity[i][j] = test[j];
        }
    }
    return m_sensitivity;
}


/*!
Method retrieves the pointer to the sensitivity solution
@return     m_NS        Number of sensitivities computed
*/
unsigned int ODE_Solver::GetNSensitivities() const
{
    return m_sensi.NParams();
}

// ERROR TOLERANCES.

double ODE_Solver::ATOL() const
{
    return m_atol;
}

void ODE_Solver::SetATOL(double atol)
{
    m_atol = atol;
}

double ODE_Solver::RTOL() const
{
    return m_rtol;
}

void ODE_Solver::SetRTOL(double rtol)
{
    m_rtol = rtol;
}


// EXTERNAL SOURCE TERMS.

// Returns the vector of external source terms (const version).
const SrcProfile *const ODE_Solver::ExtSrcTerms(void) const
{
    return m_srcterms;
}

// Sets the external source terms.
void ODE_Solver::SetExtSrcTerms(const SrcProfile &src)
{
    m_srcterms = &src;
}

// Returns the external source term function.
SrcTermFnPtr ODE_Solver::ExtSrcTermFn(void) const
{
    return _srcTerms;
}

// Sets the source term function pointer.
void ODE_Solver::SetExtSrcTermFn(SrcTermFnPtr fn)
{
    _srcTerms = fn;
}


// READ/WRITE/COPY FUNCTIONS.

// Creates a copy of the reactor object.
ODE_Solver* ODE_Solver::Clone() const
{
    return new ODE_Solver(*this);
}

// Writes the reactor to a binary data stream.
void ODE_Solver::Serialize(std::ostream &out) const
{
    const unsigned int trueval  = 1;
    const unsigned int falseval = 0;
    
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output the time.
        double val = (double)m_time;
        out.write((char*)&val, sizeof(val));

        // Output the error tolerances.
        val = (double)m_atol;
        out.write((char*)&val, sizeof(val));
        val = (double)m_rtol;
        out.write((char*)&val, sizeof(val));

        // Output equation count.
        unsigned int n = (unsigned int)m_neq;
        out.write((char*)&n, sizeof(n));

        // Output derivatives array.
        if (m_deriv != NULL) {
            out.write((char*)&trueval, sizeof(trueval));
            for (unsigned int i=0; i<m_neq; i++) {
                val = (double)m_deriv[i];
                out.write((char*)&val, sizeof(val));
            }
        } else {
            out.write((char*)falseval, sizeof(falseval));
        }
    } else {
        throw invalid_argument("Output stream not ready (Mops, Reactor::Serialize).");
    }
}

// Reads the Reactor data from a binary data stream.
void ODE_Solver::Deserialize(std::istream &in)
{
    // Clear current reactor data.
    releaseMemory();

    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        double val = 0.0;
        unsigned int n = 0;

        switch (version) {
            case 0:
                // Read the time.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_time = (double)val;

                // Read the error tolerances.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_atol = (double)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_rtol = (double)val;

                // Read equation count + special indices.
                in.read(reinterpret_cast<char*>(&m_neq), sizeof(m_neq));

                // Initialise CVODE.
                InitCVode();

                // Read derivatives array.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_deriv = new double[m_neq];
                    for (unsigned int i=0; i<m_neq; i++) {
                        in.read(reinterpret_cast<char*>(&val), sizeof(val));
                        m_deriv[i] = (double)val;
                    }
                }

                break;
            default:
                throw runtime_error("Solver serialized version number "
                                    "is invalid (Mops, ODE_Solver::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Mops, ODE_Solver::Deserialize).");
    }
}




// INITIALISATION AND MEMORY RELEASE.

// Initialises the reactor to a default state, this is called
// by the constructors.
void ODE_Solver::init(void)
{
    m_time     = 0.0;
    m_reactor  = NULL;
    m_soln     = NULL;
    m_deriv    = NULL;
    m_solvec   = NULL;
    m_atol     = 1.0e-6;
    m_rtol     = 1.0e-3;
    m_neq      = 0;
    m_srcterms = NULL;
    _srcTerms  = NULL;
    m_yvec     = NULL;
    // Space required by rate parameter sensitivity.
    m_yS       = NULL;
    m_sensitivity   = NULL;

    // Init CVODE.
    m_odewk = NULL;
}

// Releases all object memory.
void ODE_Solver::releaseMemory(void)
{
    delete [] m_deriv;
    if (m_odewk != NULL) CVodeFree(&m_odewk);

    // Free extra space
    if (m_yvec != NULL) N_VDestroy_Serial(m_yvec);
    if (m_yS != NULL) N_VDestroyVectorArray_Serial(m_yS, m_sensi.NParams());
    // Set values to defaults.
    init();
}


