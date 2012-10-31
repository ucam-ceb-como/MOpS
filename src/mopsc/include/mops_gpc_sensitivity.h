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
#include <mops_reactor.h>

#include "cvodes/cvodes.h"           /* prototypes for CVODES fcts. and consts. */
#include "cvodes/cvodes_dense.h"     /* prototype for CVDENSE fcts. and constants */
#include "nvector/nvector_serial.h"  /* defs. of serial NVECTOR fcts. and macros  */
#include "sundials/sundials_types.h" /* def. of type realtype */
#include "sundials/sundials_math.h"  /* definition of ABS */

namespace Mops
{
class SensitivityAnalyzer
{
public:
    // Enumeration of parameter type.
    enum SENS_TYPE {ARR_A,  // Pre-exponential factor.
                    ARR_n,  // Temperature exponent.
                    ARR_E,  // Activation energy.
                    INIT_C,
                    INIT_T,
                    INIT_D
    }; 
    // Arrhenius parameters for sensitivity.
    struct SENS_PARAM
    {
        unsigned int Index; // Reaction number in the mechanism.
        SENS_TYPE Type;

        // Constructors.
        SENS_PARAM(void) // Default constructor.
        {Type = ARR_A; Index=0;}; 
        SENS_PARAM(unsigned int aIndex, SENS_TYPE aType) // Initialising constructor.
        {Type=aType; Index = aIndex;};
    };

    // Enumeration of Sensitivity Types
    enum SensitivityType {Reaction_Rates, Init_Conditions};

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
    booleantype isEnableErrorControl() const;

    // Define the mechanism and parameters.
    void SetupProblem(Mops::Mechanism &mech, Mops::Reactor &reactor, const std::string &sfile);

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
    void OutputSens(std::fstream &fout, const Mops::Reactor &r, void *sim);

    // File postprocessing.
    static void PostProcess(const std::string &filename);

    // Parameter pointer to array of double. This is needed by CVODES.
    // CVODES will access and change these parameters.
    double *ParamsPtr();

    // Parameter scaling factors. approx by original params.
    double *ParamBarsPtr();

    // Set a pointer to last sensitivity output result.
    void SetSensResult(N_Vector *sens_matrix);

    // Initialise sensitivity matrix.
    void InitSensMatrix(N_Vector *sens_matrix);

    // Return sensitivity problem type.
    SensitivityType ProblemType();

private:

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

    // Non-constant mechanism which will be used for initial temperature
    // sensitivity.
    Mops::Reactor *m_reactor;

    // PARAMETERS' MEMORIES
    // Number of sensitivity parameters
    unsigned int m_NS;
    double *m_org_params;
    double *m_params;
    double *m_parambars;
    std::vector<SENS_PARAM> m_sens_params;

    // SENSITIVITY TEMP RESULT
    N_Vector *m_sens_matrix;

    // Add parameters function. m_NS should not be touch anywhere in
    // the code except here.
    bool AddParam(const SENS_PARAM &arr);

    // Read sensitivity setting from file.
    void ReadSettingV1(const CamXML::Element &elemSA);

    // Read Sensitivity Matrix block.
    // Read block of n x m from fin to matrix and simulation time.
    static void ReadSensMatrix(std::ifstream &fin, const unsigned int n, const unsigned int m, double &time, double ** matrix, double ** matrix_sqr);


};
};

#endif
