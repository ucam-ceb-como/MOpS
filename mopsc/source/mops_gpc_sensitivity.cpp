 /*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    Implementation of the SensitivitySolver class declared in the
    mops_gpc_sensitivity_solver.h header file.

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
#include "mops_gpc_sensitivity.h"
#include "mops_simulator.h"

#include "camxml.h"
#include "string_functions.h"

#include <stdexcept>
#include <math.h>

using namespace Mops;

// Constructor.
SensitivityAnalyzer::SensitivityAnalyzer() :
    m_mech(NULL),
    m_org_params(NULL),
    m_params(NULL),
    m_parambars(NULL),
    m_sens_matrix(NULL)
{
    Clear();
}

// Destructor.
SensitivityAnalyzer::~SensitivityAnalyzer()
{
    if (m_org_params != NULL) delete [] m_org_params;
    if (m_params     != NULL) delete [] m_params;
    if (m_parambars  != NULL) delete [] m_parambars;
}

// Copy constructor.
SensitivityAnalyzer::SensitivityAnalyzer(const SensitivityAnalyzer &copy)
{
    *this = copy;
}

// Operators =.
SensitivityAnalyzer &SensitivityAnalyzer::operator=(const SensitivityAnalyzer &rhs)
{
    if (this != &rhs) {
        m_probType = rhs.m_probType;
        m_enable = rhs.m_enable;
        m_err_con = rhs.m_err_con;
        m_sensi_meth = rhs.m_sensi_meth;
        m_mech = rhs.m_mech;
        // Copy pointer array
        if(m_NS != rhs.m_NS) {
            // Delete memories
            if (m_org_params != NULL) delete [] m_org_params;
            if (m_params     != NULL) delete [] m_params;
            if (m_parambars  != NULL) delete [] m_parambars;
            // Reset pointers to NULL in case of rsh.m_NS is zero.
            m_org_params = NULL;
            m_params     = NULL;
            m_parambars  = NULL;
            m_NS = rhs.m_NS;
            // No memory allocation require if m_NS is zero.
            if (m_NS != 0) {
                m_org_params = new real[m_NS];
                m_params = new real[m_NS];
                m_parambars = new real[m_NS];
            }
        }
        m_arr_params.clear();
        // Set values of varibles to rhs variables.
        for (unsigned int i = 0; i < m_NS; i++) {
            m_arr_params.push_back(rhs.m_arr_params.at(i));
            m_org_params[i] = rhs.m_org_params[i];
            m_params[i] = rhs.m_params[i];
            m_parambars[i] = rhs.m_parambars[i];
        }
        m_sens_matrix = NULL;
    }
    return *this;
}

void SensitivityAnalyzer::Clear()
{ 
    m_probType = Reaction_Rates;
    m_enable  = false;
    m_err_con = FALSE;
    m_sensi_meth = CV_SIMULTANEOUS;
    m_mech = NULL;
    m_NS = 0;
    if (m_org_params != NULL) delete [] m_org_params;
    if (m_params     != NULL) delete [] m_params;
    if (m_parambars  != NULL) delete [] m_parambars;
    m_arr_params.clear();
    m_org_params  = NULL;
    m_params      = NULL;
    m_parambars   = NULL;
    m_sens_matrix = NULL;
}

// Enable/Disable sensitivity analyzer.
void SensitivityAnalyzer::Enable(bool enable)
{
    m_enable = enable;
}

// Get enable status.
bool SensitivityAnalyzer::isEnable() const
{
    return m_enable;
}

// Set sensitivity method. This is for CVODES.
void SensitivityAnalyzer::SetMethod(int sensi_meth)
{
    m_sensi_meth = sensi_meth;
}

// Return the sensitivity method used in CVODES.
int SensitivityAnalyzer::GetMethod()
{
    return m_sensi_meth;
}

// Enable error control in CVODES.
void SensitivityAnalyzer::EnableErrorControl(booleantype err_con)
{
    m_err_con = err_con;
}

// Get enable error control status in CVODES.
booleantype SensitivityAnalyzer::isEnableErrorControl() const
{
    return m_err_con;
}

// Define the mechanism and parameters.
void SensitivityAnalyzer::SetupProblem(Mops::Mechanism &mech, const string &sfile)
{
    Clear();
    m_mech = &mech;
    CamXML::Document xmlSA;
    try {
        xmlSA.Load(sfile);
    } catch (std::exception ex) {
        // do nothing
    }
    if (xmlSA.Root()->GetChildren().size() > 0) {
        if (xmlSA.Root()->GetAttributeValue("version").compare("1.0") == 0) {
            ReadSettingV1(*xmlSA.Root());
        } else {
            throw std::runtime_error("Semsitivity XML version not supported "
                                     "(Mops, SensitivityAnalyzer::SetupProblem).");
        }
    }
}

void SensitivityAnalyzer::ReadSettingV1(const CamXML::Element &elemSA)
{
    // General pre-settings
    CamXML::Element *settingElem = NULL;
    CamXML::Element *paramsElem = NULL;
    settingElem = elemSA.GetFirstChild("settings");
    if (settingElem != NULL)
    {
        CamXML::Element *sensiElem = NULL;
        CamXML::Element *errConElem = NULL;
        CamXML::Element *probTypeElem = NULL;
        // Read Sensitivity Method
        // Default value is CV_SIMULTANEOUS.
        sensiElem = settingElem->GetFirstChild("sensitivity");
        if (sensiElem != NULL) {
            if (sensiElem->GetAttributeValue("enable").compare("true") == 0) {
                m_enable = true;
            } else {
                m_enable = false;
            }
            if (sensiElem->GetAttributeValue("method").compare("1") == 0) {
                m_sensi_meth = CV_SIMULTANEOUS;
            } else if (sensiElem->GetAttributeValue("method").compare("2") == 0) {
                m_sensi_meth = CV_STAGGERED;
            } else if (sensiElem->GetAttributeValue("method").compare("3") == 0) {
                m_sensi_meth = CV_STAGGERED1;
            } else {
                m_sensi_meth = CV_SIMULTANEOUS;
            }
        }
        // Read Error Control.
        // Default value is TRUE.
        errConElem = settingElem->GetFirstChild("errorControl");
        if (errConElem == NULL) {
            m_err_con = 1;
        } else {
            int r = errConElem->GetAttributeValue("enable").compare("true");
            m_err_con = (r != 0) ? (0) : (1);
        }
        // Read Problem Type.
        // Default value of problemtype is Reaction_Rate.
        probTypeElem = settingElem->GetFirstChild("problemType");
        if (probTypeElem != NULL) {
            if (probTypeElem->GetAttributeValue("type").compare("Reaction_Rates") == 0) {
                m_probType = Reaction_Rates;
            } else if (probTypeElem->GetAttributeValue("type").compare("Init_Concentrations") == 0) {
                m_probType = Init_Concentrations;
            } else {
                m_probType = Reaction_Rates;
            }
        }
    } else {
        throw std::runtime_error("'settings' element not found in Sensitivity XML "
                                 "(Mops, SensitivityAnalyzer::ReadSettingV1).");
    }
    paramsElem = elemSA.GetFirstChild("parameters");
    if (paramsElem != NULL)
    {
        if (m_probType == Reaction_Rates) {
            CamXML::Element *AElem = NULL;
            CamXML::Element *nElem = NULL;
            CamXML::Element *EElem = NULL;

            CamXML::Element *AddiParams = NULL;
            vector<CamXML::Element *> AParams, nParams, EParams;

            bool enableA = false;
            bool enablen = false;
            bool enableE = false;
            bool allA = false;
            bool alln = false;
            bool allE = false;

            const int nrxn = m_mech->ReactionCount();

            AElem = paramsElem->GetFirstChild("A");
            nElem = paramsElem->GetFirstChild("n");
            EElem = paramsElem->GetFirstChild("E");
            if (AElem != NULL) {
                enableA = (AElem->GetAttributeValue("enable").compare("true") == 0);
                allA    = (AElem->GetAttributeValue("all").compare("true") == 0);
            }
            if (nElem != NULL) {
                enablen = (nElem->GetAttributeValue("enable").compare("true") == 0);
                alln    = (nElem->GetAttributeValue("all").compare("true") == 0);
            }
            if (EElem != NULL) {
                enableE = (EElem->GetAttributeValue("enable").compare("true") == 0);
                allE    = (EElem->GetAttributeValue("all").compare("true") == 0);
            }
            // m_enable = m_enable && (enableA || enablen || enableE);

            // Get additional user define parameters. This xml node is only in use if allX is false.
            AddiParams = paramsElem->GetFirstChild("additional");
            if (AddiParams != NULL) {
                AddiParams->GetChildren("paramA", AParams);
                AddiParams->GetChildren("paramn", nParams);
                AddiParams->GetChildren("paramE", EParams);
            }

            // Add parameters to parameter list.
            // Add A parameters
            if (enableA) {
                unsigned int nend = (allA) ? nrxn : AParams.size();
                for (unsigned int i = 0; i < nend; i++) {
                    int rxn_index = (allA) ? (i) : ((int)Strings::cdble(AParams.at(i)->GetAttributeValue("rxnth")));
                    ARRHENIUS_PARAMS arrp(rxn_index, ARR_A);
                    AddParam(arrp);
                }
            }
            // Add n parameters
            if (enablen) {
                unsigned int nend = (alln) ? nrxn : nParams.size();
                for (unsigned int i = 0; i < nend; i++) {
                    int rxn_index = (alln) ? (i) : ((int)Strings::cdble(nParams.at(i)->GetAttributeValue("rxnth")));
                    ARRHENIUS_PARAMS arrp(rxn_index, ARR_n);
                    AddParam(arrp);
                }
            }
            // Add E parameters
            if (enableE) {
                unsigned int nend = (allE) ? nrxn : EParams.size();
                for (unsigned int i = 0; i < nend; i++) {
                    int rxn_index = (allE) ? (i) : ((int)Strings::cdble(EParams.at(i)->GetAttributeValue("rxnth")));
                    ARRHENIUS_PARAMS arrp(rxn_index, ARR_n);
                    AddParam(arrp);
                }
            }
            // Sensitivity is also disable if there is no parameter defined.
            m_enable = m_enable && (m_NS > 0);
            // Allocating memories.
            m_org_params = new Mops::real[m_NS];
            m_params     = new Mops::real[m_NS];
            m_parambars  = new Mops::real[m_NS];
            // Initializing allocated memories.
            for (unsigned int i = 0; i < m_arr_params.size(); i++) {
                Mops::real val = 0.0;
                if (m_arr_params.at(i).Type == ARR_A) {
                    val = m_mech->Reactions(m_arr_params.at(i).Rxnth)->Arrhenius().A;
                } else if (m_arr_params.at(i).Type == ARR_n) {
                    val = m_mech->Reactions(m_arr_params.at(i).Rxnth)->Arrhenius().n;
                } else {
                    val = m_mech->Reactions(m_arr_params.at(i).Rxnth)->Arrhenius().E;
                }
                m_org_params[i] = m_params[i] = m_parambars[i] = val;
            }
        } else {
            throw std::runtime_error("Unknown problem type "
                                     "(Mops, SensitivityAnalyzer::ReadSettingV1).");
        }
    } else {
        throw std::runtime_error("'parameters' element not found in Sensitivity XML "
                                 "(Mops, SensitivityAnalyzer::ReadSettingV1).");
    }
}

// Return a non-constant Mops mechanism.
Mops::Mechanism &SensitivityAnalyzer::GetMech()
{
    return *m_mech;
}

// Change mechanism parameters to the values in m_params.
void SensitivityAnalyzer::ChangeMechParams()
{
    for(unsigned int i = 0; i < m_NS; i++) {
        Sprog::Kinetics::ARRHENIUS arr = m_mech->Reactions(m_arr_params.at(i).Rxnth)->Arrhenius();
        switch ( m_arr_params.at(i).Type ) {
            case ARR_A :
                arr.A = m_params[i];
                break;
            case ARR_n :
                arr.n = m_params[i];
                break;
            case ARR_E :
                arr.E = m_params[i];
                break;
            default :
                break;
        }
        m_mech->GetReactions(i)->SetArrhenius(arr);
    }
}

// Reset mechanism parameters to its original values.
// This might not be needed.
void SensitivityAnalyzer::ResetMechParams()
{
    for(unsigned int i = 0; i < m_NS; i++) {
        Sprog::Kinetics::ARRHENIUS arr = m_mech->Reactions(m_arr_params.at(i).Rxnth)->Arrhenius();
        switch ( m_arr_params.at(i).Type ) {
            case ARR_A :
                arr.A = m_org_params[i];
                break;
            case ARR_n :
                arr.n = m_org_params[i];
                break;
            case ARR_E :
                arr.E = m_org_params[i];
                break;
            default :
                break;
        }
        m_mech->GetReactions(i)->SetArrhenius(arr);
    }
}

// File output.
void SensitivityAnalyzer::OutputSens(std::fstream &fout, const Mops::Reactor &r, void *sim)
{
    // Cast the void pointer to a Simulator object.
    Simulator *me = static_cast<Simulator*>(sim);
    if (m_NS - m_arr_params.size()) {
        throw std::runtime_error("Number of sensitivity parameters miss matched "
                             "(Mops, SensitivityAnalyzer::OutputSens).");
    }
    // number of variables (number of species + 2). // T and P
    unsigned int n_vars = r.Mech()->SpeciesCount() + 2;

    // If file is empty then output the sensitivity information
    if (fout.tellg() <= 0) {
        // enable status of sensitivity.
        fout.write((char*)&m_enable, sizeof(m_enable));

        if (m_enable) {
            // output number of simulation runs.
            unsigned int n_runs = me->RunCount();
            fout.write((char*)&n_runs, sizeof(n_runs));

            // output number of simulation iters.
            unsigned int n_iters = me->IterCount();
            fout.write((char*)&n_iters, sizeof(n_iters));

            // output number of simulation iters.
            unsigned int n_timesteps = me->TimeStepCount();
            fout.write((char*)&n_timesteps, sizeof(n_timesteps));

            // output sensitivity problemtype.
            fout.write((char*)&m_probType, sizeof(m_probType));

            // output number of sensitivity parameters.
            fout.write((char*)&m_NS, sizeof(m_NS));

            // output number of variables (number of species + 2). // T and P
            fout.write((char*)&n_vars, sizeof(n_vars));

            // output variable name list.
            for (unsigned int i = 0; i < n_vars - 2; ++i) {
                string sname = r.Mech()->GetSpecies(i)->Name();
                unsigned int len = sname.length();
                fout.write((char*)&len, sizeof(len));
                fout.write(sname.c_str(), len * sizeof(char));
            }

            // output parameter list.
            for (unsigned int i = 0; i < m_arr_params.size(); ++i) {
                fout.write((char*)&m_arr_params.at(i).Type, sizeof(m_arr_params.at(i).Type));
                fout.write((char*)&m_arr_params.at(i).Rxnth, sizeof(m_arr_params.at(i).Rxnth));
            }
        }

        // Flush unwritten buffer to file output.
        fout.flush();
    } else {
        if (m_sens_matrix != NULL) {
            // Write simulation time.
            real time = r.Time();
            fout.write((char*)&time, sizeof(time));
            //cout << "Time : " << r.Time() << endl;
            // Write main sensitivity matrix.
            for (unsigned int i = 0; i < m_NS; ++i) {
                real *sdata;
                sdata = NV_DATA_S(m_sens_matrix[i]);
                for (unsigned int j = 0; j < n_vars; ++j) {
                    real val = sdata[j];
                    fout.write((char*)&val, sizeof(val));
                    //cout << val << ", " ;
                }
                //cout << endl;
            }
            fout.flush();
        } else {
            // This should not happen at all. If so, contact your software programmers.
            //throw std::runtime_error("No sensitivity result for outputting. "
            //                     "(Mops, SensitivityAnalyzer::OutputSens).");
        }
    }
}

// File postprocessing.
void SensitivityAnalyzer::PostProcess(const std::string &filename)
{
    // Build the sensitivity binary file name.
    string finname = filename + ".sen";
    string foutname = filename + "-sensi.csv";

    // Open the sensitivity binary file.
    ifstream fin;
    fin.open(finname.c_str(), ios_base::in | ios_base::binary);
    //// Read the gas-phase conditions.
    //fvector y(N, 0.0);
    //real T=0.0, D=0.0, P=0.0;

    //in.read(reinterpret_cast<char*>(&y[0]), sizeof(y[0])*N);
    //in.read(reinterpret_cast<char*>(&T), sizeof(T));
    //in.read(reinterpret_cast<char*>(&D), sizeof(D));
    //D *= 1.0e-6; // Convert density from m^3 to cm^3.
    //in.read(reinterpret_cast<char*>(&P), sizeof(P));


    // Enable status of sensitivity.
    bool enable = false;
    // Number of simulation runs.
    unsigned int n_runs = 0;
    // Number of simulation iters.
    unsigned int n_iters = 0;
    // Number of simulation iters.
    unsigned int n_timesteps = 0;
    // Sensitivity problemtype.
    SensitivityType probType = Reaction_Rates;
    // Number of sensitivity parameters.
    unsigned int NS = 0;
    // Number of variables (number of species + 2). // T and P
    unsigned int n_vars = 2;
    // Variable name list.
    vector<string> var_names;
    // Parameter list.
    vector<ARRHENIUS_PARAMS> arr_params;

    fin.read(reinterpret_cast<char*>(&enable), sizeof(enable));

    if (enable) {

        fin.read(reinterpret_cast<char*>(&n_runs), sizeof(n_runs));

        fin.read(reinterpret_cast<char*>(&n_iters), sizeof(n_iters));

        fin.read(reinterpret_cast<char*>(&n_timesteps), sizeof(n_timesteps));

        fin.read(reinterpret_cast<char*>(&probType), sizeof(probType));

        fin.read(reinterpret_cast<char*>(&NS), sizeof(NS));

        fin.read(reinterpret_cast<char*>(&n_vars), sizeof(n_vars));

        for (unsigned int i = 0; i < n_vars - 2; ++i) {
            unsigned int len = 0;
            fin.read(reinterpret_cast<char*>(&len), sizeof(len));
            char *csname = new char[len];
            fin.read(csname, len);
            string sname;
            sname.assign(csname, len);
            var_names.push_back(sname);
            //if (len == 1) delete csname; else delete [] csname;
            delete [] csname;
        }
        var_names.push_back(string("T"));
        var_names.push_back(string("P"));


        for (unsigned int i = 0; i < NS; ++i) {
            ARRHENIUS_PARAMS arr;
            fin.read(reinterpret_cast<char*>(&arr.Type), sizeof(arr.Type));
            fin.read(reinterpret_cast<char*>(&arr.Rxnth), sizeof(arr.Rxnth));
            arr_params.push_back(arr);
        }
        // Allocating memories for calculated results
        // Total memories require = 2 * n_timesteps * NS * n_vars * sizeof(real);
        //unsigned int mem_space = 2 * n_timesteps * NS * n_vars * sizeof(real);
        //if (mem_space < 209715200) {
        //  // if require mem less than 200MB then will calculate it on RAM.
        real *** avg;
        real *** err_sd;
        real **  temp_sum;
        real **  temp_sum_sqr;
        real *times = new real [n_timesteps];
        
        // Allocate memories for overall result
        avg     = new real ** [n_timesteps];
        err_sd  = new real ** [n_timesteps];
        for (unsigned int i = 0; i < n_timesteps; ++i) {
            avg[i]     = new real * [NS];
            err_sd[i]  = new real * [NS];
            for (unsigned int j = 0; j < NS; ++j) {
                avg[i][j]     = new real [n_vars];
                err_sd[i][j] = new real [n_vars];
                for (unsigned int k = 0; k < n_vars; ++k) {
                    avg[i][j][k]     = 0.0;
                    err_sd[i][j][k]  = 0.0;
                }
            }
        }
        // Allocate memories for temporarily calculations.
        temp_sum     = new real * [NS];
        temp_sum_sqr = new real * [NS];
        for (unsigned int i = 0; i < NS; ++i) {
            temp_sum[i]     = new real [n_vars];
            temp_sum_sqr[i] = new real [n_vars];
            for (unsigned int j = 0; j < n_vars; ++j) {
                temp_sum[i][j]     = 0.0;
                temp_sum_sqr[i][j] = 0.0;
            }
        }

        // Calculating Sum and Sum square for n_runs.
        for (unsigned int r = 0; r < n_runs; ++r) {
            for (unsigned int i = 0; i < n_timesteps; ++i) {
                ReadSensMatrix(fin, NS, n_vars, times[i], temp_sum, temp_sum_sqr);
                for (unsigned int j = 0; j < NS; ++j) {
                    for (unsigned int k = 0; k < n_vars; ++k) {
                        avg[i][j][k]     += temp_sum[j][k];
                        err_sd[i][j][k]  += temp_sum_sqr[j][k];
                    }
                }
            }
        }
        // Calculate Average and SD in to sum and sum_sqr.
        fstream fout;
        fout.open(foutname.c_str(), ios_base::in | ios_base::out | ios_base::trunc);
        if (fout.good()) {
            // Output header line.
            fout << "Time (s), Parameter Type, Parameter Index" ;
            for (unsigned int i = 0; i < var_names.size(); ++i) {
                fout << ", " << var_names.at(i);
                fout << ", Err";
            }
            fout << endl;
            for (unsigned int i = 0; i < n_timesteps; ++i) {
                for (unsigned int j = 0; j < NS; ++j) {
                    fout << times[i] << "," << arr_params.at(j).Type << "," << arr_params.at(j).Rxnth;
                    for (unsigned int k = 0; k < n_vars; ++k) {
                        avg[i][j][k]     /= n_runs;
                        err_sd[i][j][k]  /= n_runs;
                        err_sd[i][j][k]  -= avg[i][j][k] * avg[i][j][k];
                        err_sd[i][j][k]  /= n_runs;
                        // fabs might be needed before taking sqrt.
                        err_sd[i][j][k]  = sqrt(err_sd[i][j][k]);
                        fout << "," << avg[i][j][k] << "," << err_sd[i][j][k];
                    }
                    fout << endl;
                }
            }
        }
        fout.close();
        // Free memories
        for (unsigned int i = 0; i < n_timesteps; ++i) {
            for (unsigned int j = 0; j < NS; ++j) {
                delete [] avg[i][j];
                delete [] err_sd[i][j];
            }
            delete [] avg[i];
            delete [] err_sd[i];
        }
        delete [] avg;
        delete [] err_sd;
        delete [] times;
        // Free memories for temporarily calculations.
        for (unsigned int i = 0; i < NS; ++i) {
            delete [] temp_sum[i];
            delete [] temp_sum_sqr[i];
        }
        delete [] temp_sum;
        delete [] temp_sum_sqr;
        //} else {
        //}
    }
    fin.close();

    //if (m_sens_matrix != NULL) {
    //    // Write simulation time.
    //    real time = r.Time();
    //    fout.write((char*)&time, sizeof(time));

    //    // Write main sensitivity matrix.
    //    for (unsigned int i = 0; i < m_NS; ++i) {
    //        real *sdata;
    //        sdata = NV_DATA_S(m_sens_matrix[i]);
    //        for (unsigned int j = 0; j < n_vars; ++j) {
    //            real val = sdata[j];
    //            fout.write((char*)&val, sizeof(val));
    //        }
    //    }
    //}
}

// Read Sensitivity Matrix block.
// Read block of n x m from fin to matrix and simulation time.
void SensitivityAnalyzer::ReadSensMatrix(ifstream &fin, const unsigned int n, const unsigned int m, real &time, real ** matrix, real ** matrix_sql)
{
    fin.read(reinterpret_cast<char*>(&time), sizeof(time));

    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < m; ++j) {
            fin.read(reinterpret_cast<char*>(&matrix[i][j]), sizeof(real));
            matrix_sql[i][j] = matrix[i][j] * matrix[i][j];
        }
    }
}

// Parameter pointer to array of real. This is needed by CVODES.
// CVODES will access and change these parameters.
real * SensitivityAnalyzer::ParamsPtr()
{
    return m_params;
}

// Set a pointer to last sensitivity output result.
void SensitivityAnalyzer::SetSensResult(N_Vector *sens_matrix)
{
    m_sens_matrix = sens_matrix;
}

// Parameter scaling factors. approx by original params.
real * SensitivityAnalyzer::ParamBarsPtr()
{
    return m_parambars;
}

bool SensitivityAnalyzer::AddParam(const ARRHENIUS_PARAMS &arrp)
{
    bool paramExist = false;
    if (m_NS == m_arr_params.size()) {
        for (unsigned int i = 0; i < m_NS; i++) {
            if ((m_arr_params.at(i).Rxnth == arrp.Rxnth) &&
                (m_arr_params.at(i).Type == arrp.Type)) {
                paramExist = true;
                break;
            }
        }
        if (!paramExist) {
            m_arr_params.push_back(arrp);
            ++m_NS;
        }
    } else {
        throw runtime_error("Number of sensitivity parameters in array and in sensitivity "
                            "object are not the same. (Mops, SensitivityAnalyzer::AddParam).");
    }
    return !paramExist;
}

unsigned int SensitivityAnalyzer::NParams()
{
    return m_NS;
}

//void SensitivityAnalyzer::ReadSensiSettings(const CamXML::Document &doc) {
//    //const CamXML::Element * elem = doc.Root();
//    //std::string probStr = elem->GetFirstChild("problemType")->GetAttributeValue("name");
//    //if (probStr.compare("Init_Concentrations")==0) {
//    //    m_probType = Init_Concentrations;
//    //} else if (probStr.compare("Reaction_Rates")==0) {
//    //    m_probType = Reaction_Rates;
//
//    //}
//}

//int SensitivityAnalyzer::Solve()
//{
//    void *cvode_mem;
//    UserData data;
//    realtype t, tout;
//    N_Vector y;
//    int iout, flag;
//
//    realtype pbar[NS];
//    //int is; 
//    N_Vector *yS;
//
//    cvode_mem = NULL;
//    data      = NULL;
//    y         =  NULL;
//    yS        = NULL;
//
//    /* Process arguments */
//    //ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);
//    Enable(TRUE);
//    EnableErrorControl(TRUE);
//    
//    /* User data structure */
//    data = (UserData) malloc(sizeof *data);
//    if (check_flag((void *)data, "malloc", 2)) return(1);
//    data->p[0] = RCONST(0.04);
//    data->p[1] = RCONST(1.0e4);
//    data->p[2] = RCONST(3.0e7);
//
//    /* Initial conditions */
//    y = N_VNew_Serial(NEQ);
//    if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
//
//    Ith(y,1) = Y1;
//    Ith(y,2) = Y2;
//    Ith(y,3) = Y3;
//
//    /* Create CVODES object */
//    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
//    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
//
//    /* Allocate space for CVODES */
//    flag = CVodeMalloc(cvode_mem, f, T0, y, CV_WF, 0.0, NULL);
//    if (check_flag(&flag, "CVodeMalloc", 1)) return(1);
//
//    /* Use private function to compute error weights */
//    flag = CVodeSetEwtFn(cvode_mem, ewt, NULL);
//    if (check_flag(&flag, "CVodeSetEwtFn", 1)) return(1);
//
//    /* Attach user data */
//    flag = CVodeSetFdata(cvode_mem, data);
//    if (check_flag(&flag, "CVodeSetFdata", 1)) return(1);
//
//    /* Attach linear solver */
//    flag = CVDense(cvode_mem, NEQ);
//    if (check_flag(&flag, "CVDense", 1)) return(1);
//
//    //flag = CVDenseSetJacFn(cvode_mem, Jac, data);
//    //if (check_flag(&flag, "CVDenseSetJacFn", 1)) return(1);
//
//    printf("\n3-species chemical kinetics problem\n");
//
//    /* Sensitivity-related settings */
//    if (isEnable()) {
//        pbar[0] = data->p[2];
//        pbar[1] = data->p[2];
//        pbar[2] = data->p[2];
//
//        yS = N_VCloneVectorArray_Serial(NS, y);
//        if (check_flag((void *)yS, "N_VCloneVectorArray_Serial", 0)) return(1);
//        //for (is=0;is<NS;is++) N_VConst(ZERO, yS[is]);
//
//        flag = CVodeSensMalloc(cvode_mem, NS, GetMethod(), yS);
//        if(check_flag(&flag, "CVodeSensMalloc", 1)) return(1);
//
//        //flag = CVodeSetSensRhs1Fn(cvode_mem, fS, data);
//        //if (check_flag(&flag, "CVodeSetSensRhs1Fn", 1)) return(1);
//        flag = CVodeSetSensErrCon(cvode_mem, isEnableErrorControl());
//        if (check_flag(&flag, "CVodeSetSensErrCon", 1)) return(1);
//        flag = CVodeSetSensParams(cvode_mem, data->p, pbar, NULL);
//        if (check_flag(&flag, "CVodeSetSensParams", 1)) return(1);
//
//        printf("Sensitivity: YES ");
//        if(GetMethod() == CV_SIMULTANEOUS)   
//            printf("( SIMULTANEOUS +");
//        else 
//            if(GetMethod() == CV_STAGGERED) printf("( STAGGERED +");
//            else                           printf("( STAGGERED1 +");   
//        if(isEnableErrorControl()) printf(" FULL ERROR CONTROL )");
//        else        printf(" PARTIAL ERROR CONTROL )");
//    } else {
//        printf("Sensitivity: NO ");
//    }
//
//    /* In loop over output points, call CVode, print results, test for error */
//
//    printf("\n\n");
//    printf("===========================================");
//    printf("============================\n");
//    printf("     T     Q       H      NST           y1");
//    printf("           y2           y3    \n");
//    printf("===========================================");
//    printf("============================\n");
//
//    for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {
//
//        flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
//        if (check_flag(&flag, "CVode", 1)) break;
//
//        PrintOutput(cvode_mem, t, y);
//
//        if (isEnable()) {
//            flag = CVodeGetSens(cvode_mem, t, yS);
//            if (check_flag(&flag, "CVodeGetSens", 1)) break;
//            PrintOutputS(yS);
//        } 
//        printf("-----------------------------------------");
//        printf("------------------------------\n");
//
//    }
//
//    /* Print final statistics */
//    PrintFinalStats(cvode_mem, isEnable());
//
//    /* Free memory */
//
//    N_VDestroy_Serial(y);                    /* Free y vector */
//    if (isEnable()) {
//        N_VDestroyVectorArray_Serial(yS, NS);  /* Free yS vector */
//    }
//    free(data);                              /* Free user data */
//    CVodeFree(&cvode_mem);                   /* Free CVODES memory */
//
//    return(0);
//}

 /*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

 /*
 * f routine. Compute f(t,y). 
 */

//int SensitivityAnalyzer::f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
//{
//    realtype y1, y2, y3, yd1, yd3;
//    UserData data;
//    realtype k1, k2, k3;
//
//    y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
//    data = (UserData) f_data;
//    k1 = data->p[0]; k2 = data->p[1]; k3 = data->p[2];
//
//    yd1 = Ith(ydot,1) = -k1*y1 + k2*y2*y3;
//    yd3 = Ith(ydot,3) = k3*y2*y2;
//        Ith(ydot,2) = -yd1 - yd3;
//
//    return(0);
//}


 /* 
 * Jacobian routine. Compute J(t,y). 
 */

//int SensitivityAnalyzer::Jac(long int N, DenseMat J, realtype t,
//               N_Vector y, N_Vector fy, void *jac_data, 
//               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
//{
//    realtype y1, y2, y3;
//    UserData data;
//    realtype p1, p2, p3;
//
//    y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
//    data = (UserData) jac_data;
//    p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
//
//    IJth(J,1,1) = -p1;  IJth(J,1,2) = p2*y3;          IJth(J,1,3) = p2*y2;
//    IJth(J,2,1) =  p1;  IJth(J,2,2) = -p2*y3-2*p3*y2; IJth(J,2,3) = -p2*y2;
//                      IJth(J,3,2) = 2*p3*y2;
//
//    return(0);
//}

/* 
 * fS routine. Compute sensitivity r.h.s. 
 */

//int SensitivityAnalyzer::fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
//              int iS, N_Vector yS, N_Vector ySdot, 
//              void *fS_data, N_Vector tmp1, N_Vector tmp2)
//{
//    UserData data;
//    realtype p1, p2, p3;
//    realtype y1, y2, y3;
//    realtype s1, s2, s3;
//    realtype sd1, sd2, sd3;
//
//    data = (UserData) fS_data;
//    p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
//
//    y1 = Ith(y,1);  y2 = Ith(y,2);  y3 = Ith(y,3);
//    s1 = Ith(yS,1); s2 = Ith(yS,2); s3 = Ith(yS,3);
//
//    sd1 = -p1*s1 + p2*y3*s2 + p2*y2*s3;
//    sd3 = 2*p3*y2*s2;
//    sd2 = -sd1-sd3;
//
//    switch (iS) {
//    case 0:
//    sd1 += -y1;
//    sd2 +=  y1;
//    break;
//    case 1:
//    sd1 +=  y2*y3;
//    sd2 += -y2*y3;
//    break;
//    case 2:
//    sd2 += -y2*y2;
//    sd3 +=  y2*y2;
//    break;
//    }
//
//    Ith(ySdot,1) = sd1;
//    Ith(ySdot,2) = sd2;
//    Ith(ySdot,3) = sd3;
//
//    return(0);
//}
//
 /*
 * EwtSet function. Computes the error weights at the current solution.
 */

//int SensitivityAnalyzer::ewt(N_Vector y, N_Vector w, void *e_data)
//{
//    int i;
//    realtype yy, ww, rtol, atol[3];
//
//    rtol    = RTOLK;
//    atol[0] = ATOL1;
//    atol[1] = ATOL2;
//    atol[2] = ATOL3;
//
//    for (i=1; i<=3; i++) {
//        yy = Ith(y,i);
//        ww = rtol * ABS(yy) + atol[i-1];  
//        if (ww <= 0.0) return (-1);
//        Ith(w,i) = 1.0/ww;
//    }
//
//    return(0);
//}

 /*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

 /*
 * Print current t, step count, order, stepsize, and solution.
 */

//void SensitivityAnalyzer::PrintOutput(void *cvode_mem, realtype t, N_Vector u)
//{
//    long int nst;
//    int qu, flag;
//    realtype hu, *udata;
//
//    udata = NV_DATA_S(u);
//
//    flag = CVodeGetNumSteps(cvode_mem, &nst);
//    check_flag(&flag, "CVodeGetNumSteps", 1);
//    flag = CVodeGetLastOrder(cvode_mem, &qu);
//    check_flag(&flag, "CVodeGetLastOrder", 1);
//    flag = CVodeGetLastStep(cvode_mem, &hu);
//    check_flag(&flag, "CVodeGetLastStep", 1);
//
//#if defined(SUNDIALS_EXTENDED_PRECISION)
//    printf("%8.3Le %2d  %8.3Le %5ld\n", t, qu, hu, nst);
//#elif defined(SUNDIALS_DOUBLE_PRECISION)
//    printf("%8.3le %2d  %8.3le %5ld\n", t, qu, hu, nst);
//#else
//    printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
//#endif
//
//    printf("                  Solution       ");
//
//#if defined(SUNDIALS_EXTENDED_PRECISION)
//    printf("%12.4Le %12.4Le %12.4Le \n", udata[0], udata[1], udata[2]);
//#elif defined(SUNDIALS_DOUBLE_PRECISION)
//    printf("%12.4le %12.4le %12.4le \n", udata[0], udata[1], udata[2]);
//#else
//    printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
//#endif
//
//}
//
// /* 
// * Print sensitivities.
// */
//
//void SensitivityAnalyzer::PrintOutputS(N_Vector *uS)
//{
//    realtype *sdata;
//
//    sdata = NV_DATA_S(uS[0]);
//    printf("                  Sensitivity 1  ");
//
//#if defined(SUNDIALS_EXTENDED_PRECISION)
//    printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
//#elif defined(SUNDIALS_DOUBLE_PRECISION)
//    printf("%12.4le %12.4le %12.4le %12.4le %12.4le \n", sdata[0], sdata[1], sdata[2], sdata[3], sdata[4]);
//#else
//    printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
//#endif
//  
//    sdata = NV_DATA_S(uS[1]);
//    printf("                  Sensitivity 2  ");
//
//#if defined(SUNDIALS_EXTENDED_PRECISION)
//    printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
//#elif defined(SUNDIALS_DOUBLE_PRECISION)
//    printf("%12.4le %12.4le %12.4le %12.4le %12.4le \n", sdata[0], sdata[1], sdata[2], sdata[3], sdata[4]);
//#else
//    printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
//#endif
//
//    sdata = NV_DATA_S(uS[2]);
//    printf("                  Sensitivity 3  ");
//
//#if defined(SUNDIALS_EXTENDED_PRECISION)
//    printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
//#elif defined(SUNDIALS_DOUBLE_PRECISION)
//    printf("%12.4le %12.4le %12.4le %12.4le %12.4le \n", sdata[0], sdata[1], sdata[2], sdata[3], sdata[4]);
//#else
//    printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
//#endif
//}
//
// /* 
// * Print some final statistics from the CVODES memory.
// */
//
//void SensitivityAnalyzer::PrintFinalStats(void *cvode_mem, booleantype sensi)
//{
//    long int nst;
//    long int nfe, nsetups, nni, ncfn, netf;
//    long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
//    long int nje, nfeLS;
//    int flag;
//
//    flag = CVodeGetNumSteps(cvode_mem, &nst);
//    check_flag(&flag, "CVodeGetNumSteps", 1);
//    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
//    check_flag(&flag, "CVodeGetNumRhsEvals", 1);
//    flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
//    check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
//    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
//    check_flag(&flag, "CVodeGetNumErrTestFails", 1);
//    flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
//    check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
//    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
//    check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
//
//    if (sensi) {
//        flag = CVodeGetNumSensRhsEvals(cvode_mem, &nfSe);
//        check_flag(&flag, "CVodeGetNumSensRhsEvals", 1);
//        flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
//        check_flag(&flag, "CVodeGetNumRhsEvalsSens", 1);
//        flag = CVodeGetNumSensLinSolvSetups(cvode_mem, &nsetupsS);
//        check_flag(&flag, "CVodeGetNumSensLinSolvSetups", 1);
//        flag = CVodeGetNumSensErrTestFails(cvode_mem, &netfS);
//        check_flag(&flag, "CVodeGetNumSensErrTestFails", 1);
//        flag = CVodeGetNumSensNonlinSolvIters(cvode_mem, &nniS);
//        check_flag(&flag, "CVodeGetNumSensNonlinSolvIters", 1);
//        flag = CVodeGetNumSensNonlinSolvConvFails(cvode_mem, &ncfnS);
//        check_flag(&flag, "CVodeGetNumSensNonlinSolvConvFails", 1);
//    }
//
//    flag = CVDenseGetNumJacEvals(cvode_mem, &nje);
//    check_flag(&flag, "CVDenseGetNumJacEvals", 1);
//    flag = CVDenseGetNumRhsEvals(cvode_mem, &nfeLS);
//    check_flag(&flag, "CVDenseGetNumRhsEvals", 1);
//
//    printf("\nFinal Statistics\n\n");
//    printf("nst     = %5ld\n\n", nst);
//    printf("nfe     = %5ld\n",   nfe);
//    printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
//    printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);
//
//    if(sensi) {
//        printf("\n");
//        printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
//        printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
//        printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
//    }
//
//    printf("\n");
//    printf("nje    = %5ld    nfeLS     = %5ld\n", nje, nfeLS);
//
//}
//
// /* 
// * Check function return value.
// *    opt == 0 means SUNDIALS function allocates memory so check if
// *             returned NULL pointer
// *    opt == 1 means SUNDIALS function returns a flag so check if
// *             flag >= 0
// *    opt == 2 means function allocates memory so check if returned
// *             NULL pointer 
// */
//
//int SensitivityAnalyzer::check_flag(void *flagvalue, char *funcname, int opt)
//{
//    int *errflag;
//
//    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
//    if (opt == 0 && flagvalue == NULL) {
//        fprintf(stderr, 
//            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
//        funcname);
//        return(1); }
//
//    /* Check if flag < 0 */
//    else if (opt == 1) {
//        errflag = (int *) flagvalue;
//        if (*errflag < 0) {
//            fprintf(stderr, 
//                "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
//            funcname, *errflag);
//        return(1); }}
//
//    /* Check if function returned NULL pointer - no memory allocated */
//    else if (opt == 2 && flagvalue == NULL) {
//        fprintf(stderr, 
//            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
//        funcname);
//        return(1); }
//
//    return(0);
//}
