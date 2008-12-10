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

// CVODE includes.
#include "nvector/nvector_serial.h"
#include "cvodes/cvodes.h"
#include "cvodes_impl.h"
#include "cvodes/cvodes_dense.h"
#include "cvodes_dense_impl.h"

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
        //m_abstolQ  = rhs.m_abstolQ;
        //m_reltolQ  = rhs.m_reltolQ;
        //m_abstolB  = rhs.m_abstolB;
        //m_reltolB  = rhs.m_reltolB;
        //m_abstolQB = rhs.m_abstolQB;
        //m_isInitB  = rhs.m_isInitB;

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
            m_yS[i] = N_VClone_Serial(rhs.m_yS[i]);
        }
        if (m_yvec != NULL) N_VDestroy_Serial(m_yvec);
        m_yvec = N_VClone_Serial(rhs.m_yvec);
        //if (m_solvec != NULL) N_VDestroy_Serial(m_solvec);
        //VCopy_Serial(rhs.m_solvec, m_solvec);

        // Copy ODE workspace.
        if (rhs.m_odewk != NULL) {
            if (m_odewk == NULL) InitCVode();
            assignCVMem(*((CVodeMem)rhs.m_odewk));
        } else {
            if (m_odewk != NULL) CVodeFree(&m_odewk);
            m_odewk = NULL;
        }

        // Copy ODE workspace, incl. derivatives.
        memcpy(m_deriv, rhs.m_deriv, sizeof(real)*m_neq);
    }
    return *this;
}

// SOLVER SET UP.

// Initialises the reactor to the given time.
void ODE_Solver::Initialise(const Reactor &reac)
{
    // Store useful variables.
    m_time = reac.Time();
    m_neq  = reac.ODE_Count();
    
    // Store solution vector.
    m_soln = reac.Mixture()->RawData();

    // Fill derivative vector.
    delete [] m_deriv;
    m_deriv = new real[m_neq];
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
        CVodeMalloc(m_odewk, &rhsFn_CVODES, m_time, 
                    m_yvec, 
                    CV_SS, m_rtol, (void*)&m_atol);
    } else {
        CVodeMalloc(m_odewk, &rhsFn_CVODE, m_time, 
                    m_yvec, 
                    CV_SS, m_rtol, (void*)&m_atol);
    }

    // Set up internal solution array.
    m_solvec = N_VNewEmpty_Serial(m_neq);
    NV_DATA_S(m_solvec) = m_soln;

    // Set the f_data pointer to this Reactor object.
    CVodeSetFdata(m_odewk, (void*)this);

    // Set other parameters.
    CVodeSetMaxNumSteps(m_odewk, 2000);

    // Set CVDense as the linear system solver.
    CVDense(m_odewk, m_neq);

    // Set the Jacobian function in CVODE.
    // Benchmark results from certain test case :
    // - No sensitivity analysis : External Jacobian is faster than CVODE internal jacobain.
    // - Sensitivity analyis (Rate parameters) : CVODE internal jacobian is fastest (114 s). jacFn_CVODES is slightly slower (120 s)
    //   and jacFn_CVODE is very slow (146 s). Thus, Jacobian function will be set according to this test.
    if (m_sensi.isEnable()) {
        if (m_sensi.ProblemType() == SensitivityAnalyzer::Reaction_Rates) {
            // Internal one is fastest so don't set jacobian function.
            // CVDenseSetJacFn(m_odewk, jacFn_CVODES, (void*)this);
        } else if (m_sensi.ProblemType() == SensitivityAnalyzer::Init_Conditions) {
            // CVDenseSetJacFn(m_odewk, jacFn_CVODES, (void*)this); <== need test
        }
    } else {
        CVDenseSetJacFn(m_odewk, jacFn_CVODE, (void*)this);
    }
    int flag = 0;
    if (m_sensi.isEnable()) {
        m_yS = N_VCloneVectorArray_Serial(m_sensi.NParams(), m_yvec);
        m_sensi.InitSensMatrix(m_yS);

        flag = CVodeSensMalloc(m_odewk, m_sensi.NParams(), m_sensi.GetMethod(), m_yS);

        //flag = CVodeSetSensRhs1Fn(cvode_mem, fS, data);

        flag = CVodeSetSensErrCon(m_odewk, m_sensi.isEnableErrorControl());

        flag = CVodeSetSensParams(m_odewk, m_sensi.ParamsPtr(), m_sensi.ParamBarsPtr(), NULL);

    } else {
        //printf("Sensitivity: NO ");
    }
}

// Initialises CVode for backward problem.
//void ODE_Solver::ReInitCVodeB(real t_final)
//{
//    if (!m_isInitB) {
//        // Allocate global memory
//
//        // Space required by initial condition sensitivity.
//        int steps = 10;
//        m_odeadj = CVadjMalloc(m_odewk, steps, CV_HERMITE);
//        // Initialize m_yB
//        m_yB = N_VNew_Serial(m_neq);
//        N_VConst(ZERO, m_yB);
//
//        // Initialize m_qB
//        m_qB = N_VNew_Serial(m_neq);
//        N_VConst(ZERO, m_qB);
//
//        // Setup problem
//        CVodeCreateB(m_odeadj, CV_BDF, CV_NEWTON);
//
//        CVodeMallocB(m_odeadj, &rhsLamdaFn_CVODES, t_final, m_yB, CV_SS, m_reltolB, &m_abstolB);
//
//        CVodeSetFdataB(m_odeadj, (void*)this);
//
//        CVDenseB(m_odeadj, m_neq);
//
//        //CVDenseSetJacFnB(m_odeadj, &rhs_, data);
//
//        CVodeQuadMallocB(m_odeadj, &rhsQuadBFn_CVODES, m_qB);
//
//        CVodeSetQuadFdataB(m_odeadj, (void*)this);
//
//        CVodeSetQuadErrConB(m_odeadj, TRUE, CV_SS, m_reltolB, &m_abstolQB);
//
//        m_isInitB = true;
//    } else {
//        // Backward problem has been initialised before.
//        N_VConst(ZERO, m_yB);
//
//        N_VConst(ZERO, m_qB);
//
//        CVodeReInitB(m_odeadj, &rhsLamdaFn_CVODES, t_final, m_yB, CV_SS, m_reltolB, &m_abstolB);
//
//        CVodeQuadReInitB(m_odeadj, &rhsQuadBFn_CVODES, m_qB); 
//   }
//}

// Reset the solver.  Need to do this if the the reactor
// contents has been changed between calls to Solve().
void ODE_Solver::ResetSolver(void)
{
    if (m_yvec != NULL) N_VDestroy_Serial(m_yvec);
    m_yvec = N_VMake_Serial(m_neq, m_soln);
    // m_yS cannot be reset since it need to know the previous values
    // in order to continue solving the next step.
    //if (m_yS != NULL) N_VDestroyVectorArray_Serial(m_yS, m_sensi.NParams());
    //m_yS = N_VCloneVectorArray_Serial(m_sensi.NParams(), m_yvec);
    for (unsigned int is=0;is<m_sensi.NParams();is++) N_VConst(ZERO, m_yS[is]);
    if (m_sensi.isEnable()) {
        CVodeReInit(m_odewk, &rhsFn_CVODES, m_time, 
                    m_yvec,
                    CV_SS, m_rtol, (void*)&m_atol);
        CVodeSensReInit(m_odewk, m_sensi.GetMethod(), m_yS);
        // Reintialise sensitivity related workspace.
        //if (m_sensi.ProblemType() == SensitivityAnalyzer::Reaction_Rates) {
        //    CVodeSensReInit(m_odewk, m_sensi.GetMethod(), m_yS);
        //} else if (m_sensi.ProblemType() == SensitivityAnalyzer::Init_Conditions) {
        //    CVodeQuadReInit(m_odewk, &rhsQuadFn_CVODES, m_q);
        //}
    } else {
        CVodeReInit(m_odewk, &rhsFn_CVODE, m_time, 
                    m_yvec,
                    CV_SS, m_rtol, (void*)&m_atol);
    }
    
} 

// Reset the solver.  Need to do this if the the reactor
// contents has been changed between calls to Solve().
void ODE_Solver::ResetSolver(const Reactor &reac)
{
    // Check that this reactor has the same problem size
    // as the last reactor.  If not then we have to resize the
    // workspace, which is done by the Initialise() routine.
    if (reac.ODE_Count() == m_neq) {
        m_time = reac.Time();
        m_soln = reac.Mixture()->RawData();
        if (m_yvec != NULL) N_VDestroy_Serial(m_yvec);
        m_yvec = N_VMake_Serial(m_neq, m_soln);
        // m_yS cannot be reset since it need to know the previous values.
        //if (m_yS != NULL) N_VDestroyVectorArray_Serial(m_yS, m_sensi.NParams());
        //m_yS = N_VCloneVectorArray_Serial(m_sensi.NParams(), m_yvec);
        for (unsigned int is=0;is<m_sensi.NParams();is++) N_VConst(ZERO, m_yS[is]);
        if (m_sensi.isEnable()) {
            CVodeReInit(m_odewk, &rhsFn_CVODES, m_time, 
                        m_yvec,
                        CV_SS, m_rtol, (void*)&m_atol);
            CVodeSensReInit(m_odewk, m_sensi.GetMethod(), m_yS);
            // Reintialise sensitivity related workspace.
            //if (m_sensi.ProblemType() == SensitivityAnalyzer::Reaction_Rates) {
            //    CVodeSensReInit(m_odewk, m_sensi.GetMethod(), m_yS);
            //} else if (m_sensi.ProblemType() == SensitivityAnalyzer::Init_Conditions) {
            //    CVodeQuadReInit(m_odewk, &rhsQuadFn_CVODES, m_q);
            //}
        } else {
            CVodeReInit(m_odewk, &rhsFn_CVODE, m_time, 
                        m_yvec,
                        CV_SS, m_rtol, (void*)&m_atol);
        }
    } else {
        Initialise(reac);
    }
}


// SOLVING THE REACTOR.

// Solves the reactor up to the given time.
void ODE_Solver::Solve(Reactor &reac, real stop_time)
{
    // Check that the reactor has the same problem size
    // as set in the solver.
    if (m_neq != reac.ODE_Count()) {
        // Re-initialise the solver for a new problem size.
        Initialise(reac);
    } else {
        // Get the time and solution pointer from the reactor.
        m_time = reac.Time();
        m_soln = reac.Mixture()->RawData();
    }

    // Put the solution into an N_Vector data pointer.  This does not
    // involve copying the data.
    m_reactor = &reac;
    NV_DATA_S(m_solvec) = m_soln;
    CVodeSetStopTime(m_odewk, stop_time);

    // Solve over time step.
    while (m_time < stop_time) {
        int CVode_error = 0;
        CVode_error = CVode(m_odewk, stop_time, m_solvec, &m_time, CV_NORMAL_TSTOP);
        //if (m_sensi.isEnable() && (m_sensi.ProblemType() == SensitivityAnalyzer::Init_Conditions)) {
        //    // Reinitialing is always needed for initial concentration problem due to change from previous backward solver.
        //    CVodeReInit(m_odewk, &rhsFn_CVODES, m_time, m_solvec, CV_SS, m_rtol, (void*)&m_atol);

        //    CVodeQuadReInit(m_odewk, &rhsQuadFn_CVODES, m_q); 

        //    int ncheck = 0;  
        //    CVode_error = CVodeF(m_odeadj, stop_time, m_solvec, &m_time, CV_NORMAL_TSTOP, &ncheck);

        //    CVodeGetQuad(m_odewk, stop_time, m_q);

        //    ReInitCVodeB(stop_time);

        //    CVode_error = CVodeB(m_odeadj, m_sensi.StartTime(), m_yB, &m_time, CV_NORMAL_TSTOP);

        //    CVodeGetQuadB(m_odeadj, m_qB);

        //} else {
        //    CVode_error = CVode(m_odewk, stop_time, m_solvec, &m_time, CV_NORMAL_TSTOP);
        //}
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
        reac.Mixture()->Normalise(); // This should not be required if CVODE solves correctly.
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
        int flag = CVodeGetSens(m_odewk, m_time, m_yS);
        //if (m_sensi.check_flag(&flag, "CVodeGetSens", 1)) break;
        m_sensi.SetSensResult(m_yS);
    }

}


// ERROR TOLERANCES.

real ODE_Solver::ATOL() const
{
    return m_atol;
}

void ODE_Solver::SetATOL(real atol)
{
    m_atol = atol;
}

real ODE_Solver::RTOL() const
{
    return m_rtol;
}

void ODE_Solver::SetRTOL(real rtol)
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
                m_time = (real)val;

                // Read the error tolerances.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_atol = (real)val;
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_rtol = (real)val;

                // Read equation count + special indices.
                in.read(reinterpret_cast<char*>(&m_neq), sizeof(m_neq));

                // Initialise CVODE.
                InitCVode();

                // Read derivatives array.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));
                if (n == 1) {
                    m_deriv = new real[m_neq];
                    for (unsigned int i=0; i<m_neq; i++) {
                        in.read(reinterpret_cast<char*>(&val), sizeof(val));
                        m_deriv[i] = (real)val;
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
    // Space required by initial condition sensitivity.
    //m_q        = NULL;  // Quadrature solution of the concentration integration.
    //m_odeadj   = NULL;  // CVODES worksapce for adjoint sensitivity.
    //m_yB       = NULL;  // Backwards solution of lamda for sensitivity.
    //m_qB       = NULL;  // Backwards integration solution of concentration.
    //m_abstolQ  = 1.0e-6;
    //m_reltolQ  = 1.0e-3;
    //m_abstolB  = 1.0e-6;
    //m_reltolB  = 1.0e-3;
    //m_abstolQB = 1.0e-6;
    //m_isInitB  = false;

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

// Copies the given CVode workspace into this ODE_Solver object. 
// (This function cannot be used with sensitivity problem.)
void ODE_Solver::assignCVMem(const CVodeMemRec &copy)
{
    if (m_odewk != NULL) {
        // Get structure from ODE workspace which belongs to this
        // ODE_Solver object.
        CVodeMemRec &mem = *((CVodeMem)m_odewk);

        // Machine unit round-off.
        mem.cv_uround = copy.cv_uround;

        // PROBLEM SPECIFICATION DATA.

        mem.cv_f      = copy.cv_f;      // RHS evaluator function.
        mem.cv_f_data = copy.cv_f_data; // User data pointer, passed to f.
        mem.cv_lmm    = copy.cv_lmm;    // CV_ADAMS or CV_BDF.
        mem.cv_iter   = copy.cv_iter;   // CV_FUNCTIONAL or CV_NEWTON.
        mem.cv_itol   = copy.cv_itol;   // CV_SS or CV_SV.

        mem.cv_reltol  = copy.cv_reltol;   // Relative tolerance.
        mem.cv_Sabstol = copy.cv_Sabstol;  // Scalar absolute tolerance.

        // Clear and re-assign vector abs. tolerance.
        if (mem.cv_lrw1 == copy.cv_lrw1) {
            // Workspaces are same size.  Just perform memory copy.
            if (copy.cv_VabstolMallocDone) {
                if (mem.cv_VabstolMallocDone) {
                    // Vector exists in both workspaces.
                    memcpy(NV_DATA_S(mem.cv_Vabstol), NV_DATA_S(copy.cv_Vabstol), 
                           copy.cv_lrw1*sizeof(realtype));
                } else {
                    // Vector only exists in copy workspace.
                    mem.cv_Vabstol = N_VClone_Serial(copy.cv_Vabstol);
                }
            } else {
                if (mem.cv_VabstolMallocDone) {
                    // Vector exists in mem workspace only.
                    N_VDestroy_Serial(mem.cv_Vabstol);
                    mem.cv_Vabstol = NULL;
                }
            }
        } else {
            // Workspaces are diferent sizes, need to resize vectors.
            if (mem.cv_VabstolMallocDone) N_VDestroy_Serial(mem.cv_Vabstol);
            if (copy.cv_VabstolMallocDone) {
                // Vector exists in copy workspace.
                mem.cv_Vabstol = N_VClone_Serial(copy.cv_Vabstol);
            } else {
                // Vector does not exist in copy workspace.
                mem.cv_Vabstol = NULL;
            }
        }

        mem.cv_efun    = copy.cv_efun;   // Function to set ewt.
        mem.cv_e_data  = copy.cv_e_data; // User data pointer passed to efun.
        
        //QUADRATURE RELATED DATA
        
        mem.cv_quadr    = copy.cv_quadr;    // TRUE if integrating quadratures

        mem.cv_fQ       = copy.cv_fQ;
        mem.cv_fQ_data  = copy.cv_fQ_data;  // user pointer passed to fQ
        mem.cv_itolQ    = copy.cv_itolQ;
        mem.cv_errconQ  = copy.cv_errconQ;

        mem.cv_reltolQ  = copy.cv_reltolQ;  // relative tolerance for quadratures
        mem.cv_SabstolQ = copy.cv_SabstolQ; // scalar absolute tolerance for quadratures

        // Clear and re-assign vector abs. tolerance.
        if (mem.cv_lrw1Q == copy.cv_lrw1Q) {
            // Workspaces are same size.  Just perform memory copy.
            if (copy.cv_VabstolQMallocDone) {
                if (mem.cv_VabstolQMallocDone) {
                    // Vector exists in both workspaces.
                    memcpy(NV_DATA_S(mem.cv_VabstolQ), NV_DATA_S(copy.cv_VabstolQ), 
                           copy.cv_lrw1Q*sizeof(realtype));
                } else {
                    // Vector only exists in copy workspace.
                    mem.cv_VabstolQ = N_VClone_Serial(copy.cv_VabstolQ);
                }
            } else {
                if (mem.cv_VabstolQMallocDone) {
                    // Vector exists in mem workspace only.
                    N_VDestroy_Serial(mem.cv_VabstolQ);
                    mem.cv_VabstolQ = NULL;
                }
            }
        } else {
            // Workspaces are diferent sizes, need to resize vectors.
            if (mem.cv_VabstolQMallocDone) N_VDestroy_Serial(mem.cv_VabstolQ);
            if (copy.cv_VabstolQMallocDone) {
                // Vector exists in copy workspace.
                mem.cv_VabstolQ = N_VClone_Serial(copy.cv_VabstolQ);
            } else {
                // Vector does not exist in copy workspace.
                mem.cv_VabstolQ = NULL;
            }
        }

        // SENSITIVITY RELATED DATA 

        mem.cv_sensi        = copy.cv_sensi;        // TRUE if computing sensitivities

        mem.cv_Ns           = copy.cv_Ns;           // Number of sensitivities

        mem.cv_ism          = copy.cv_ism;          // ism = SIMULTANEOUS or STAGGERED

        mem.cv_fS           = copy.cv_fS;           // fS = (df/dy)*yS + (df/dp)
        mem.cv_fS1          = copy.cv_fS1;          // fS1 = (df/dy)*yS_i + (df/dp)
        mem.cv_user_fS_data = copy.cv_user_fS_data; // user data pointer for fS
        mem.cv_fS_data      = copy.cv_fS_data;      // actual data pointer passed to fS
        mem.cv_fSDQ         = copy.cv_fSDQ;
        mem.cv_ifS          = copy.cv_ifS;          // ifS = ALLSENS or ONESENS

        mem.cv_p            = copy.cv_p;            // parameters in f(t,y,p), no copy here since cvode do not free this memory
        // Scale factors for parameters
        delete [] mem.cv_pbar; mem.cv_pbar = NULL;
        if (mem.cv_Ns > 0) {
            mem.cv_pbar = new realtype[mem.cv_Ns];
            memcpy(mem.cv_pbar, copy.cv_pbar, mem.cv_Ns*sizeof(realtype));
        }
        // List of sensitivities
        delete [] mem.cv_plist; mem.cv_plist = NULL;
        if (mem.cv_Ns > 0) {
            mem.cv_plist = new int[mem.cv_Ns];
            memcpy(mem.cv_plist, copy.cv_plist, mem.cv_Ns*sizeof(int));
        }
        mem.cv_DQtype   = copy.cv_DQtype;            // central/forward finite differences
        mem.cv_DQrhomax = copy.cv_DQrhomax;         // cut-off value for separate/simultaneous FD

        mem.cv_errconS  = copy.cv_errconS;          // TRUE if sensitivities are in err. control

        mem.cv_itolS  = copy.cv_itolS;
        mem.cv_reltolS  = copy.cv_reltolS;          // relative tolerance for sensitivities
        // Scalar absolute tolerances for sensi.
        delete [] mem.cv_SabstolS; mem.cv_SabstolS = NULL;
        if (mem.cv_Ns > 0) {
            mem.cv_SabstolS = new realtype[mem.cv_Ns];
            memcpy(mem.cv_SabstolS, copy.cv_SabstolS, mem.cv_Ns*sizeof(realtype));
        }
    //N_Vector *cv_VabstolS;   // vector absolute tolerances for sensi.

        // VECTORS.

        if (mem.cv_lrw1 == copy.cv_lrw1) {
            // Workspaces are same size.  Just need to perform
            // memory copies of vectors.

            // NORDSIECK HISTORY ARRAY.

            // Copy Nordsieck array vectors.
            for (int i=0; i!=mem.cv_qmax+1; ++i) {
                memcpy(NV_DATA_S(mem.cv_zn[i]), NV_DATA_S(copy.cv_zn[i]), 
                       NV_LENGTH_S(copy.cv_zn[i])*sizeof(realtype));
            }

            // OTHER VECTORS OF LENGTH N.

            // Error weight vector.
            if (copy.cv_ewt != NULL) {
                // Vector exists in copy workspace.  Copy.
                if (mem.cv_ewt != NULL) {
                    // Vector exists in mem workspace.  memcpy.
                    memcpy(NV_DATA_S(mem.cv_ewt), NV_DATA_S(copy.cv_ewt), 
                           NV_LENGTH_S(copy.cv_ewt)*sizeof(realtype));
                } else {
                    // Vector does not exist in mem workspace.  Clone.
                    mem.cv_ewt = N_VClone_Serial(copy.cv_ewt);
                }
            } else {
                // Vector does not exist in copy workspace.  Delete.
                if (mem.cv_ewt != NULL) N_VDestroy_Serial(mem.cv_ewt);
            }

            // Temporary storage vector.
            //if (copy.cv_y != NULL) {
            //    // Vector exists in copy workspace.  Copy.
            //    if (mem.cv_y != NULL) {
            //        // Vector exists in mem workspace.  memcpy.
            //        memcpy(NV_DATA_S(mem.cv_y), NV_DATA_S(copy.cv_y), 
            //               NV_LENGTH_S(copy.cv_y)*sizeof(realtype));
            //    } else {
            //        // Vector does not exist in mem workspace.  Clone.
            //        mem.cv_y = N_VClone_Serial(copy.cv_y);
            //    }
            //} else {
            //    // Vector does not exist in copy workspace.  Delete.
            //    if (mem.cv_y != NULL) N_VDestroy_Serial(mem.cv_y);
            //}

            // Estimated local error.
            if (copy.cv_acor != NULL) {
                // Vector exists in copy workspace.  Copy.
                if (mem.cv_acor != NULL) {
                    // Vector exists in mem workspace.  memcpy.
                    memcpy(NV_DATA_S(mem.cv_acor), NV_DATA_S(copy.cv_acor), 
                           NV_LENGTH_S(copy.cv_acor)*sizeof(realtype));
                } else {
                    // Vector does not exist in mem workspace.  Clone.
                    mem.cv_acor = N_VClone_Serial(copy.cv_acor);
                }
            } else {
                // Vector does not exist in copy workspace.  Delete.
                if (mem.cv_acor != NULL) N_VDestroy_Serial(mem.cv_acor);
            }

            // Temporary storage vector.
            if (copy.cv_tempv != NULL) {
                // Vector exists in copy workspace.  Copy.
                if (mem.cv_tempv != NULL) {
                    // Vector exists in mem workspace.  memcpy.
                    memcpy(NV_DATA_S(mem.cv_tempv), NV_DATA_S(copy.cv_tempv), 
                           NV_LENGTH_S(copy.cv_tempv)*sizeof(realtype));
                } else {
                    // Vector does not exist in mem workspace.  Clone.
                    mem.cv_tempv = N_VClone_Serial(copy.cv_tempv);
                }
            } else {
                // Vector does not exist in copy workspace.  Delete.
                if (mem.cv_tempv != NULL) N_VDestroy_Serial(mem.cv_tempv);
            }

            // Temporary storage vector.
            if (copy.cv_ftemp != NULL) {
                // Vector exists in copy workspace.  Copy.
                if (mem.cv_ftemp != NULL) {
                    // Vector exists in mem workspace.  memcpy.
                    memcpy(NV_DATA_S(mem.cv_ftemp), NV_DATA_S(copy.cv_ftemp), 
                           NV_LENGTH_S(copy.cv_ftemp)*sizeof(realtype));
                } else {
                    // Vector does not exist in mem workspace.  Clone.
                    mem.cv_ftemp = N_VClone_Serial(copy.cv_ftemp);
                }
            } else {
                // Vector does not exist in copy workspace.  Delete.
                if (mem.cv_ftemp != NULL) N_VDestroy_Serial(mem.cv_ftemp);
            }
        } else {
            // Workspaces are different sizes.  Need to reallocate
            // all vectors.  Use Clone facility of N_Vector.

            // NORDSIECK HISTORY ARRAY.

            // Resize and re-assign Nordsieck array vectors.
            for (int i=0; i!=mem.cv_qmax+1; ++i) {
                N_VDestroy_Serial(mem.cv_zn[i]);
                mem.cv_zn[i] = N_VClone_Serial(copy.cv_zn[i]);
            }

            // OTHER VECTORS OF LENGTH N.

            // Error weight vector.
            if (mem.cv_ewt != NULL) N_VDestroy_Serial(mem.cv_ewt);
            if (copy.cv_ewt != NULL) mem.cv_ewt = N_VClone_Serial(copy.cv_ewt);
            else mem.cv_ewt = NULL;

            // Temporary storage vector.
            if (mem.cv_y != NULL) N_VDestroy_Serial(mem.cv_y);
            if (copy.cv_y != NULL) mem.cv_y = N_VClone_Serial(copy.cv_y);
            else mem.cv_y = NULL;

            // Estimated local vector.
            if (mem.cv_acor != NULL) N_VDestroy_Serial(mem.cv_acor);
            if (copy.cv_acor != NULL) mem.cv_acor = N_VClone_Serial(copy.cv_acor);
            else mem.cv_acor = NULL;

            // Temporary storage vector.
            if (mem.cv_tempv != NULL) N_VDestroy_Serial(mem.cv_tempv);
            if (copy.cv_tempv != NULL) mem.cv_tempv = N_VClone_Serial(copy.cv_tempv);
            else mem.cv_tempv = NULL;

            // Temporary storage vector.
            if (mem.cv_ftemp != NULL) N_VDestroy_Serial(mem.cv_ftemp);
            if (copy.cv_ftemp != NULL) mem.cv_ftemp = N_VClone_Serial(copy.cv_ftemp);
            else mem.cv_ftemp = NULL;
        }

        
        // QUADRATURE RELATED VECTORS 

        // Assumming Workspaces are different sizes.  Need to reallocate
        // all vectors.  Use Clone facility of N_Vector.

        // NORDSIECK HISTORY ARRAY FOR QUADRATURES.

        // Resize and re-assign Nordsieck array vectors.
        for (int i=0; i!=mem.cv_qmax+1; ++i) {
            N_VDestroy_Serial(mem.cv_znQ[i]);
            mem.cv_znQ[i] = N_VClone_Serial(copy.cv_znQ[i]);
        }

        // OTHER VECTORS OF LENGTH N.

        // Error weight vector for quadratures
        if (mem.cv_ewtQ != NULL) N_VDestroy_Serial(mem.cv_ewtQ);
        if (copy.cv_ewtQ != NULL) mem.cv_ewtQ = N_VClone_Serial(copy.cv_ewtQ);
        else mem.cv_ewtQ = NULL;

        // Temporary storage vector. Unlike y, yQ is not allocated by the user.
        if (mem.cv_yQ != NULL) N_VDestroy_Serial(mem.cv_yQ);
        if (copy.cv_yQ != NULL) mem.cv_yQ = N_VClone_Serial(copy.cv_yQ);
        else mem.cv_yQ = NULL;

        // Estimated local vector. acorQ = yQ_n(m) - yQ_n(0).
        if (mem.cv_acorQ != NULL) N_VDestroy_Serial(mem.cv_acorQ);
        if (copy.cv_acorQ != NULL) mem.cv_acorQ = N_VClone_Serial(copy.cv_acorQ);
        else mem.cv_acorQ = NULL;

        // Temporary storage vector.
        if (mem.cv_tempvQ != NULL) N_VDestroy_Serial(mem.cv_tempvQ);
        if (copy.cv_tempvQ != NULL) mem.cv_tempvQ = N_VClone_Serial(copy.cv_tempvQ);
        else mem.cv_tempvQ = NULL;

        // DOES CVODESENSMALLOC ALLOCATE ADDITIONAL SPACE?

        mem.cv_stgr1alloc = copy.cv_stgr1alloc; // Did we allocate ncfS1, ncfnS1, and nniS1?

        // TSTOP INFORMATION.

        mem.cv_istop    = copy.cv_istop;
        mem.cv_tstopset = copy.cv_tstopset;
        mem.cv_tstop    = copy.cv_tstop;

        // STEP DATA.

        mem.cv_q      = copy.cv_q;      // Current order.
        mem.cv_qprime = copy.cv_qprime; // Order to be used on next step.
        mem.cv_next_q = copy.cv_next_q; // Order to be used on next step.
        mem.cv_qwait  = copy.cv_qwait;  // Number of internal steps before considering change in q.
        mem.cv_L      = copy.cv_L;      // = q + 1.

        mem.cv_hin      = copy.cv_hin;      // Initial step size.
        mem.cv_h        = copy.cv_h;        // Current step size.
        mem.cv_hprime   = copy.cv_hprime;   // Step size to be used on next step.
        mem.cv_next_h   = copy.cv_next_h;   // Step size to be used on next step.
        mem.cv_eta      = copy.cv_eta;      // = hprime / h.
        mem.cv_hscale   = copy.cv_hscale;   // Value of h used in zn.
        mem.cv_tn       = copy.cv_tn;       // Current internal value of t.
        mem.cv_tretlast = copy.cv_tretlast; // Value of tret last returned by CVode.

        // Array of previous q+1 successfull step sizes.
        for (unsigned int i=0; i!=L_MAX+1; ++i) {
            mem.cv_tau[i] = copy.cv_tau[i];
        }

        // Array of test quantities.
        for (unsigned int i=0; i!=NUM_TESTS+1; ++i) {
            mem.cv_tq[i] = copy.cv_tq[i];
        }

        // Coefficients of l(x) (degree q poly).
        for (unsigned int i=0; i!=L_MAX; ++i) {
            mem.cv_l[i] = copy.cv_l[i];
        }

        mem.cv_rl1    = copy.cv_rl1;    // The scalar 1/l[1].
        mem.cv_gamma  = copy.cv_gamma;  // = h * rl1.
        mem.cv_gammap = copy.cv_gammap; // gamma at last setup call.
        mem.cv_gamrat = copy.cv_gamrat; // gamma / gammap.

        mem.cv_crate   = copy.cv_crate;   // Estimated corrector convergence rate.
        mem.cv_crateS  = copy.cv_crateS;  // est. corrector conv. rate in NlsStgr
        mem.cv_acnrm   = copy.cv_acnrm;   // | acor | nrms.
        mem.cv_acnrmS  = copy.cv_acnrmS;  // | acorS |
        mem.cv_acnrmQ  = copy.cv_acnrmQ;  // | acorQ |
        mem.cv_nlscoef = copy.cv_nlscoef; // Coefficient in non-linear convergence test.
        mem.cv_mnewt   = copy.cv_mnewt;   // Newton iteration counter.
        // Array of Ns local counters for conv. failures (used in CVStep for STAGGERED1)
        delete [] mem.cv_ncfS1; mem.cv_ncfS1 = NULL;
        if (mem.cv_Ns > 0) {
            mem.cv_ncfS1 = new int[mem.cv_Ns];
            memcpy(mem.cv_ncfS1, copy.cv_ncfS1, mem.cv_Ns*sizeof(int));
        }

        // LIMITS.

        mem.cv_qmax   = copy.cv_qmax;     // q <= qmax.
        mem.cv_mxstep = copy.cv_mxstep;   // Max. number of internal steps for one user call.
        mem.cv_maxcor = copy.cv_maxcor;   // Max number of corrector iterations for soln. of non-linear eqn.

        mem.cv_maxcorS = copy.cv_maxcorS;
        mem.cv_mxhnil  = copy.cv_mxhnil;   // Max. warning msg. count.
        mem.cv_maxnef  = copy.cv_maxnef;   // Max. error test failure count.
        mem.cv_maxncf  = copy.cv_maxncf;   // Max. non-linear convergence failure count.

        mem.cv_hmin     = copy.cv_hmin;     // |h| >= hmin.
        mem.cv_hmax_inv = copy.cv_hmax_inv; // |h| <= 1/hmax_inv.
        mem.cv_etamax   = copy.cv_etamax;   // eta <= etamax.

        // COUNTERS.

        mem.cv_nst      = copy.cv_nst;      // Internal steps taken count.
        mem.cv_nfe      = copy.cv_nfe;      // F call count.
        mem.cv_nfSe      = copy.cv_nfSe;    // number of fS calls
        mem.cv_nfQe      = copy.cv_nfQe;    // number of fQ calls
        mem.cv_nfeS      = copy.cv_nfeS;    // number of f calls from sensi DQ

        mem.cv_ncfn     = copy.cv_ncfn;     // Corrected convergence failure count.
        mem.cv_ncfnS    = copy.cv_ncfnS;   // Number of total sensi. corr. conv. failures
        // Number of sensi. corrector conv. failures
        delete [] mem.cv_ncfnS1; mem.cv_ncfnS1 = NULL;
        if (mem.cv_Ns > 0) {
            mem.cv_ncfnS1 = new long int[mem.cv_Ns];
            memcpy(mem.cv_ncfnS1, copy.cv_ncfnS1, mem.cv_Ns*sizeof(long int));
        }

        mem.cv_nni      = copy.cv_nni;      // Newton iterations performed count.
        mem.cv_nniS     = copy.cv_nniS;     // Number of total sensi. nonlinear iterations
        // Number of sensi. nonlinear iterations
        delete [] mem.cv_nniS1; mem.cv_nniS1 = NULL;
        if (mem.cv_Ns > 0) {
            mem.cv_nniS1 = new long int[mem.cv_Ns];
            memcpy(mem.cv_nniS1, copy.cv_nniS1, mem.cv_Ns*sizeof(long int));
        }

        mem.cv_netf     = copy.cv_netf;     // Error test failure count.
        mem.cv_netfS    = copy.cv_netfS;    // Number of sensi. error test failures
        mem.cv_netfQ    = copy.cv_netfQ;    // Number of quadr. error test failures

        mem.cv_nsetups  = copy.cv_nsetups;  // Setup call count.
        mem.cv_nsetupsS = copy.cv_nsetupsS; // Number of setup calls due to sensitivities
        
        mem.cv_nhnil    = copy.cv_nhnil;   // Msgs. issued to user count.

        // SPACE REQUIREMENTS FOR CVODE.

        mem.cv_lrw1  = copy.cv_lrw1;  // No. of realtype words in 1 N_Vector.
        mem.cv_liw1  = copy.cv_liw1;  // No. of integer words in 1 N_Vector.
        mem.cv_lrw1Q = copy.cv_lrw1Q; // no. of realtype words in 1 N_Vector yQ
        mem.cv_liw1Q = copy.cv_liw1Q; // no. of integer words in 1 N_Vector yQ
        mem.cv_lrw   = copy.cv_lrw;   // No. of realtype words in CVODE work vectors.
        mem.cv_liw   = copy.cv_liw;   // No. of integer words in CVODE work vectors.

        // LINEAR SOLVER FUNCTIONS.

        mem.cv_linit  = copy.cv_linit;
        mem.cv_lsetup = copy.cv_lsetup;
        mem.cv_lsolve = copy.cv_lsolve;
        mem.cv_lfree  = copy.cv_lfree;
        
        //STEP SIZE RATIOS

        mem.cv_etaqm1  = copy.cv_etaqm1;      // Ratio of new to old h for order q-1
        mem.cv_etaq    = copy.cv_etaq;        // Ratio of new to old h for order q
        mem.cv_etaqp1  = copy.cv_etaqp1;      // Ratio of new to old h for order q+1

        // LINEAR SOLVER SPECIFIC MEMORY.

        if (copy.cv_lmem != NULL) {
            // Get references to the linear solver data.
            CVDenseMemRec &lmem  = *((CVDenseMem)mem.cv_lmem);
            CVDenseMemRec &lcopy = *((CVDenseMem)copy.cv_lmem);
            
            // Copy linear solver data.
            long N     = lmem.d_n;
            lmem.d_n   = lcopy.d_n;   // Problem dimension.
            lmem.d_jac = lcopy.d_jac; // Jacobian routine.

            // M = I - gamma J, gamma = h / l1.
            if (lcopy.d_M != NULL) {
                if ((lcopy.d_M->M == lmem.d_M->M) && (lcopy.d_M->N == lmem.d_M->N)) {
                    // Just copy the memory, no need to reallocate.
                    DenseCopy(lcopy.d_M, lmem.d_M);
                } else {
                    // Need to reallocate data.  First delete current data.
                    DenseFreeMat(lmem.d_M);
                    // Allocate memory to matrix..
                    lmem.d_M = DenseAllocMat(lcopy.d_M->M, lcopy.d_M->N);
                    // Copy data.
                    DenseCopy(lcopy.d_M, lmem.d_M);
                }
            } else {
                // Clear matrix.
                DenseFreeMat(lmem.d_M);
            }

            // Pivot array for PM = LU.
            if (lcopy.d_pivots != NULL) {
                if (lcopy.d_n != N) {
                    // Need to reallocate memory.
                    DenseFreePiv(lmem.d_pivots);
                    lmem.d_pivots = DenseAllocPiv(lcopy.d_n);
                }
                memcpy(lmem.d_pivots, lcopy.d_pivots, lcopy.d_n*sizeof(long int));
            } else {
                // Free pivots memory.
                DenseFreePiv(lmem.d_pivots);
            }

            // Old Jacobian.
            if (lcopy.d_savedJ != NULL) {
                if ((lcopy.d_savedJ->M == lmem.d_savedJ->M) && (lcopy.d_savedJ->N == lmem.d_savedJ->N)) {
                    // Just copy the memory, no need to reallocate.
                    DenseCopy(lcopy.d_savedJ, lmem.d_savedJ);
                } else {
                    // Need to reallocate data.  First delete current data.
                    DenseFreeMat(lmem.d_savedJ);
                    // Allocate memory to matrix..
                    lmem.d_savedJ = DenseAllocMat(lcopy.d_savedJ->M, lcopy.d_savedJ->N);
                    // Copy data.
                    DenseCopy(lcopy.d_savedJ, lmem.d_savedJ);
                }
            } else {
                // Clear matrix.
                DenseFreeMat(lmem.d_savedJ);
            }

            lmem.d_nstlj     = lcopy.d_nstlj;     // = nst at last Jacobian evaluation.
            lmem.d_nje       = lcopy.d_nje;       // Number of calls to jac.
            lmem.d_nfeD      = lcopy.d_nfeD;      // Number of calls to f due to difference quotient approximation of J.
            if (lcopy.d_J_data != &copy) {
                lmem.d_J_data    = lcopy.d_J_data; // J_data is passed to jac.
            } else {
                lmem.d_J_data    = (void*)&mem; // J_data is passed to jac.
            }
            lmem.d_last_flag = lcopy.d_last_flag; // Last error return flag.
        } else {
            mem.cv_lfree(&mem);
        }

        mem.cv_forceSetup = copy.cv_forceSetup;     // Flag to request a call to the setup routine


        // SAVED VALUES.

        mem.cv_qu           = copy.cv_qu;           // Last successful q value used.
        mem.cv_nstlp        = copy.cv_nstlp;        // Step number of last setup call.
        mem.cv_h0u          = copy.cv_h0u;          // Actual initial step size.
        mem.cv_hu           = copy.cv_hu;           // Last successfull h value used.
        mem.cv_saved_tq5    = copy.cv_saved_tq5;    // Saved value of tq[5].
        mem.cv_jcur         = copy.cv_jcur;         // Is Jacobian info used by linear solver current?
        mem.cv_tolsf        = copy.cv_tolsf;        // tolerance scale factor.
        mem.cv_qmax_alloc   = copy.cv_qmax_alloc;   // Value of qmax used when allocating memory.
        mem.cv_qmax_allocQ  = copy.cv_qmax_allocQ;  // value of qmax used when allocating quad. memory
        mem.cv_qmax_allocS  = copy.cv_qmax_allocS;  // value of qmax used when allocating sensi. memory
        mem.cv_indx_acor    = copy.cv_indx_acor;    // Index of xn vector in which acor is stored.
        mem.cv_setupNonNull = copy.cv_setupNonNull; // Does setup do something?

        mem.cv_VabstolMallocDone     = copy.cv_VabstolMallocDone;
        mem.cv_MallocDone            = copy.cv_MallocDone;

        mem.cv_VabstolQMallocDone    = copy.cv_VabstolQMallocDone;
        mem.cv_quadMallocDone        = copy.cv_quadMallocDone;

        mem.cv_VabstolSMallocDone    = copy.cv_VabstolSMallocDone;
        mem.cv_SabstolSMallocDone    = copy.cv_SabstolSMallocDone;
        mem.cv_sensMallocDone        = copy.cv_sensMallocDone;

        // ERROR HANDLER FUNCTION AND OUTPUT FILE.

        mem.cv_ehfun   = copy.cv_ehfun; // Error message handler.
        mem.cv_eh_data = copy.cv_eh_data; // User data pointer passed to ehfun.
        mem.cv_errfp   = mem.cv_errfp; // File buffer to which CVODE error messages are sent.
        
        // STABILITY LIMIT DETECTION.

        mem.cv_sldeton = copy.cv_sldeton; // Is stability limit detection on?

        // Scaled data array for STALD.
        for (unsigned int i=0; i!=6; ++i) {
            for (unsigned int j=0; j!=4; ++j) {
                mem.cv_ssdat[i][j] = copy.cv_ssdat[i][j];
            }
        }
        
        mem.cv_nscon = copy.cv_nscon; // Counter for STALD method.
        mem.cv_nor = copy.cv_nor; // Counter for order reduction count.

        // ROOT-FINDING DATA.

        mem.cv_gfun   = copy.cv_gfun;   // Function g for roots sought.
        int nrt       = mem.cv_nrtfn;
        mem.cv_nrtfn  = copy.cv_nrtfn;  // Number of components of g.
        mem.cv_g_data = copy.cv_g_data; // Pointer to user data for g.
        
        // int array for root information.
        if (copy.cv_nrtfn != nrt) {
            // Need to reallocate memory first.
            delete [] mem.cv_iroots; mem.cv_iroots = NULL;
            if (mem.cv_nrtfn > 0) {
                mem.cv_iroots = new int[mem.cv_nrtfn];
                memcpy(mem.cv_iroots, copy.cv_iroots, mem.cv_nrtfn*sizeof(int));
            }
        } else {
            if (nrt > 0) {
                memcpy(mem.cv_iroots, copy.cv_iroots, mem.cv_nrtfn*sizeof(int));
            }
        }

        mem.cv_tlo   = copy.cv_tlo;   // Nearest endpoint of iterval in root search.
        mem.cv_thi   = copy.cv_thi;   // Farthest endpoint of iterval in root search.
        mem.cv_trout = copy.cv_trout; // t value returned by rootfinding routine.

        // Saved array of g values at t = tlo.
        delete [] mem.cv_glo; mem.cv_glo = NULL;
        if (mem.cv_nrtfn > 0) {
            mem.cv_glo = new realtype[mem.cv_nrtfn];
            memcpy(mem.cv_glo, copy.cv_glo, mem.cv_nrtfn*sizeof(realtype));
        }

        // Saved array of g values at t= thi.
        delete [] mem.cv_ghi; mem.cv_ghi = NULL;
        if (mem.cv_nrtfn > 0) {
            mem.cv_ghi = new realtype[mem.cv_nrtfn];
            memcpy(mem.cv_ghi, copy.cv_ghi, mem.cv_nrtfn*sizeof(realtype));
        }

        // Array of g values at t = trout.
        delete [] mem.cv_grout; mem.cv_grout = NULL;
        if (mem.cv_nrtfn > 0) {
            mem.cv_grout = new realtype[mem.cv_nrtfn];
            memcpy(mem.cv_grout, copy.cv_grout, mem.cv_nrtfn*sizeof(realtype));
        }

        mem.cv_toutc = copy.cv_toutc; // Copy of tout (if NORMAL mode).
        mem.cv_ttol  = copy.cv_ttol;  // Tolerance on root location.
        mem.cv_taskc = copy.cv_taskc; // Copy of parameter task.
        mem.cv_irfnd = copy.cv_irfnd; // Flag showing whether last step had a root.
        mem.cv_nge   = copy.cv_nge;   // Counter for g evaluations.
    }
}
