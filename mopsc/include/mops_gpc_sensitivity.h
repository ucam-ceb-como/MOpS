/*
  Author(s):      Weerapong Phadungsukana (wp214)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    Gas phase sensitivity analyzer class. It contains setup function for
    mops_gpc_sensitivity_solver.h and mops_gpc_sensitivity_solver.cpp.

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

#ifndef MOPS_GPC_SENSITIVITY_H
#define MOPS_GPC_SENSITIVITY_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <camxml.h>
#include <mops_mechanism.h>

#include <cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDENSE fcts. and constants */
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

 /* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i=1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */

 /* Problem Constants */

#define NEQ   3             /* number of equations  */
#define Y1    RCONST(1.0)   /* initial y components */
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)
#define RTOLK  RCONST(1e-4)  /* scalar relative tolerance */
#define ATOL1 RCONST(1e-8)  /* vector absolute tolerance components */
#define ATOL2 RCONST(1e-14)
#define ATOL3 RCONST(1e-6)
#define T0    RCONST(0.0)   /* initial time */
#define T1    RCONST(0.4)   /* first output time */
#define TMULT RCONST(10.0)  /* output time factor */
#define NOUT  12            /* number of output times */

#define NP    3             /* number of problem parameters */
#define NS    3             /* number of sensitivities computed */


namespace Mops
{
class SensitivityAnalyzer
{
public:
    // Enumeration of parameter type.
    enum ARRHENIUS_TYPE {ARR_A,  // Pre-exponential factor.
                         ARR_n,  // Temperature exponent.
                         ARR_E}; // Activation energy.
    // Arrhenius parameters for sensitivity.
    struct ARRHENIUS_PARAMS
    {
        unsigned int Rxnth; // Reaction number in the mechanism.
        ARRHENIUS_TYPE Type;

        // Constructors.
        ARRHENIUS_PARAMS(void) // Default constructor.
        {Type = ARR_A; Rxnth=0;}; 
        ARRHENIUS_PARAMS(unsigned int aRxnth, ARRHENIUS_TYPE aType) // Initialising constructor.
        {Type=aType; Rxnth = aRxnth;};
    };

    // Enumeration of Sensitivity Types
    enum SensitivityType {Reaction_Rates, Init_Concentrations};

    // Constructor.
    SensitivityAnalyzer();

    // Destructor.
    ~SensitivityAnalyzer();
    
    // Copy constructor.
    SensitivityAnalyzer(const SensitivityAnalyzer &copy);
    
    // Operators.
    SensitivityAnalyzer &operator=(const SensitivityAnalyzer &rhs);

    // Reset all setting.
    void Clear();
    
    // Enable/Disable sensitivity analyzer.
    void Enable(bool enable);

    // Get enable status.
    bool isEnable() const;

    // Set sensitivity method. This is for CVODES.
    void SetMethod(int sensi_meth);

    // Return the sensitivity method used in CVODES.
    int GetMethod();

    // Enable error control in CVODES.
    void EnableErrorControl(booleantype err_con);

    // Get enable error control status in CVODES.
    booleantype isEnableErrorControl();

    // Define the mechanism and parameters.
    void SetupProblem(Mops::Mechanism &mech, const string &sfile);

    // Return number of parameters.
    unsigned int NParams();

    // Return a non-constant Mops mechanism.
    Mops::Mechanism &GetMech();

    // Change mechanism parameters to the values in m_params.
    void ChangeMechParams();

    // Reset mechanism parameters to its original values.
    // This might not be needed.
    void ResetMechParams();

    // File output.
    void outputTo() {};

    // Set output.
    void SetOutputFile(const string &ofile);

    // Parameter pointer to array of real. This is needed by CVODES.
    // CVODES will access and change these parameters.
    real *ParamsPtr();

    // Parameter scaling factors. approx by original params.
    real *ParamBarsPtr();

    // Test solve.
    int Solve();
    static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
    static void PrintOutputS(N_Vector *uS);
    static void PrintFinalStats(void *cvode_mem, booleantype sensi);
    static int check_flag(void *flagvalue, char *funcname, int opt);
private:
    typedef struct {
      realtype p[3];           /* problem parameters */
    } *UserData;

    // Define problem type.
    SensitivityType m_probType;

    // Enable sensitivity analyzer status.
    bool m_enable;

    // Error control method.
    booleantype m_err_con;

    // Sensitivity method.
    int m_sensi_meth;

    // Non-constant mechanism which will be used by rshFn_CVODES in
    // mops_rhs_func.cpp to access reaction set variables.
    // This can only be set via SetMechanism and where Mechanism
    // object is not made constant, ex from mops.cpp.
    Mops::Mechanism *m_mech;

    // PARAMETERS' MEMORIES
    // Number of sensitivity parameters
    unsigned int m_NS;
    real *m_org_params;
    real *m_params;
    real *m_parambars;
    vector<ARRHENIUS_PARAMS> m_arr_params;

    bool AddParam(const ARRHENIUS_PARAMS &arr);
     /* Prototypes of functions by CVODES */

    // Setup sensitivity problem from given setting xml version 1.0.
    // Clear function must be called before using this function, otherwise 
    // program might behave abnormally. m_mech pointer must be set after 
    // use clear function.
    void ReadSettingV1(const CamXML::Element &elemSA);

    static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);
    //static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

    //static int Jac(long int N, DenseMat J, realtype t,
    //               N_Vector y, N_Vector fy, void *jac_data, 
    //               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    //static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
    //              int iS, N_Vector yS, N_Vector ySdot, 
    //              void *fS_data, N_Vector tmp1, N_Vector tmp2);

    static int ewt(N_Vector y, N_Vector w, void *e_data);

     /* Prototypes of private functions */


};
};

#endif
