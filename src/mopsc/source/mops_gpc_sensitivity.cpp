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
        m_probType   = rhs.m_probType;
        m_enable     = rhs.m_enable;
        m_err_con    = rhs.m_err_con;
        m_sensi_meth = rhs.m_sensi_meth;
        m_mech       = rhs.m_mech;
        m_reactor    = rhs.m_reactor;
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
                m_org_params = new double[m_NS];
                m_params = new double[m_NS];
                m_parambars = new double[m_NS];
            }
        }
        m_sens_params.clear();
        // Set values of varibles to rhs variables.
        for (unsigned int i = 0; i < m_NS; i++) {
            m_sens_params.push_back(rhs.m_sens_params.at(i));
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
    m_probType      = Reaction_Rates;
    m_enable        = false;
    m_err_con       = FALSE;
    m_sensi_meth    = CV_SIMULTANEOUS;
    m_mech          = NULL;
    m_reactor       = NULL;
    m_NS            = 0;
    if (m_org_params != NULL) delete [] m_org_params;
    if (m_params     != NULL) delete [] m_params;
    if (m_parambars  != NULL) delete [] m_parambars;
    m_sens_params.clear();
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
void SensitivityAnalyzer::SetupProblem(Mops::Mechanism &mech, Mops::Reactor &reactor, const std::string &sfile)
{
    Clear();
    m_mech = &mech;
    m_reactor = &reactor;
    CamXML::Document xmlSA;
    xmlSA.Load(sfile);

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
        //CamXML::Element *errTolElem = NULL;
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
            } else if (probTypeElem->GetAttributeValue("type").compare("Init_Conditions") == 0) {
                m_probType = Init_Conditions;
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
            std::vector<CamXML::Element *> AParams, nParams, EParams;

            bool enableA = false;
            bool enablen = false;
            bool enableE = false;
            bool allA = false;
            bool alln = false;
            bool allE = false;

            const int nrxn = m_mech->GasMech().ReactionCount();

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
                    SENS_PARAM arrp(rxn_index, ARR_A);
                    AddParam(arrp);
                }
            }
            // Add n parameters
            if (enablen) {
                unsigned int nend = (alln) ? nrxn : nParams.size();
                for (unsigned int i = 0; i < nend; i++) {
                    int rxn_index = (alln) ? (i) : ((int)Strings::cdble(nParams.at(i)->GetAttributeValue("rxnth")));
                    SENS_PARAM arrp(rxn_index, ARR_n);
                    AddParam(arrp);
                }
            }
            // Add E parameters
            if (enableE) {
                unsigned int nend = (allE) ? nrxn : EParams.size();
                for (unsigned int i = 0; i < nend; i++) {
                    int rxn_index = (allE) ? (i) : ((int)Strings::cdble(EParams.at(i)->GetAttributeValue("rxnth")));
                    SENS_PARAM arrp(rxn_index, ARR_n);
                    AddParam(arrp);
                }
            }
            // Sensitivity is also disable if there is no parameter defined.
            m_enable = m_enable && (m_NS > 0);
            // Allocating memories.
            m_org_params = new double[m_NS];
            m_params     = new double[m_NS];
            m_parambars  = new double[m_NS];
            // Initializing allocated memories.
            for (unsigned int i = 0; i < m_sens_params.size(); i++) {
                double val = 0.0;
                if (m_sens_params.at(i).Type == ARR_A) {
                    val = m_mech->GasMech().Reactions(m_sens_params.at(i).Index)->Arrhenius().A;
                } else if (m_sens_params.at(i).Type == ARR_n) {
                    val = m_mech->GasMech().Reactions(m_sens_params.at(i).Index)->Arrhenius().n;
                } else {
                    val = m_mech->GasMech().Reactions(m_sens_params.at(i).Index)->Arrhenius().E;
                }
                m_org_params[i] = m_params[i] = m_parambars[i] = val;
                if (val ==0.0){
                m_parambars[i] = 0.9; //!This varies at the user's discretion.
                }
            }
        } else if (m_probType == Init_Conditions) {
            std::vector<CamXML::Element *> initParams;

            // Add initial condition parameter list
            paramsElem->GetChildren("init", initParams );
            for (unsigned int i = 0; i < initParams.size(); i++) {
                std::string id = Strings::trim(initParams.at(i)->GetAttributeValue("id"));
                if (id.compare("T") == 0) {
                    // initial temperature parameter, only work if a non-constant temperature is applied.
                    SENS_PARAM arrp(0, INIT_T);
                    AddParam(arrp);
                } else if (id.compare("D") == 0) {
                    // initial density parameter
                    SENS_PARAM arrp(0, INIT_D);
                    AddParam(arrp);
                } else {
                    // initial concentration parameters.
                    int sp_index = m_mech->GasMech().FindSpecies(id);
                    if (sp_index > -1) {
                        SENS_PARAM arrp(sp_index, INIT_C);
                        AddParam(arrp);
                    } else {
                        // ignore non-existing species.
                    }
                }
            }

            // Sensitivity is also disable if there is no parameter defined.
            m_enable = m_enable && (m_NS > 0);
            // Allocating memories.
            m_org_params = new double[m_NS];
            m_params     = new double[m_NS];
            m_parambars  = new double[m_NS];
            // Initializing allocated memories.
            unsigned int i_temp = m_reactor->Mixture()->GasPhase().temperatureIndex();
            unsigned int i_dens = m_reactor->Mixture()->GasPhase().densityIndex();
            for (unsigned int i = 0; i < m_sens_params.size(); i++) {
                double val = 0.0;
                if (m_sens_params.at(i).Type == INIT_T) {
                    val = m_reactor->Mixture()->GasPhase().RawData()[i_temp];
                } else if (m_sens_params.at(i).Type == INIT_D) {
                    val = m_reactor->Mixture()->GasPhase().RawData()[i_dens];
                } else {
                    val = m_reactor->Mixture()->GasPhase().RawData()[m_sens_params.at(i).Index];
                }
                m_org_params[i] = m_params[i] = m_parambars[i] = val;
                if (val == 0.0) {
                    m_parambars[i] = 0.9; //!This value varies at the user's discretion.
                }
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
    if (m_probType == Reaction_Rates) {
        for(unsigned int i = 0; i < m_NS; i++) {
            Sprog::Kinetics::ARRHENIUS arr = m_mech->GasMech().Reactions(m_sens_params.at(i).Index)->Arrhenius();
            switch ( m_sens_params.at(i).Type ) {
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
            m_mech->GasMech().GetReactions(i)->SetArrhenius(arr);
        }
    } else if (m_probType == Init_Conditions) {
        unsigned int i_temp = m_reactor->Mixture()->GasPhase().temperatureIndex();
        unsigned int i_dens = m_reactor->Mixture()->GasPhase().densityIndex();
        for(unsigned int i = 0; i < m_NS; i++) {
            if (m_sens_params.at(i).Type == INIT_T) {
                m_reactor->Mixture()->GasPhase().RawData()[i_temp] = m_params[i];
            } else if (m_sens_params.at(i).Type == INIT_D) {
                m_reactor->Mixture()->GasPhase().RawData()[i_dens] = m_params[i];
            } else {
                m_reactor->Mixture()->GasPhase().RawData()[m_sens_params.at(i).Index] = m_params[i];
            }
        }
    } else {
        throw std::runtime_error("Undefined Type: Impossible case. Contact you programmer."
                                 "(Mops, SensitivityAnalyzer::ChangeMechParams).");
    }
}

// Reset mechanism parameters to its original values.
// This might not be needed.
void SensitivityAnalyzer::ResetMechParams()
{
    if (m_probType == Reaction_Rates) {
        for(unsigned int i = 0; i < m_NS; i++) {
            Sprog::Kinetics::ARRHENIUS arr = m_mech->GasMech().Reactions(m_sens_params.at(i).Index)->Arrhenius();
            switch ( m_sens_params.at(i).Type ) {
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
            m_mech->GasMech().GetReactions(i)->SetArrhenius(arr);
        }
    } else if (m_probType == Init_Conditions) {
        unsigned int i_temp = m_reactor->Mixture()->GasPhase().temperatureIndex();
        unsigned int i_dens = m_reactor->Mixture()->GasPhase().densityIndex();
        for(unsigned int i = 0; i < m_NS; i++) {
            if (m_sens_params.at(i).Type == INIT_T) {
                m_reactor->Mixture()->GasPhase().RawData()[i_temp] = m_org_params[i];
            } else if (m_sens_params.at(i).Type == INIT_D) {
                m_reactor->Mixture()->GasPhase().RawData()[i_dens] = m_org_params[i];
            } else {
                m_reactor->Mixture()->GasPhase().RawData()[m_sens_params.at(i).Index] = m_org_params[i];
            }
        }
    } else {
        throw std::runtime_error("Undefined Type: Impossible case. Contact you programmer."
                                 "(Mops, SensitivityAnalyzer::ResetMechParams).");
    }
}

// File output.
void SensitivityAnalyzer::OutputSens(std::fstream &fout, const Mops::Reactor &r, void *sim)
{
    // Cast the void pointer to a Simulator object.
    Simulator *me = static_cast<Simulator*>(sim);
    if (m_NS - m_sens_params.size()) {
        throw std::runtime_error("Number of sensitivity parameters miss matched "
                             "(Mops, SensitivityAnalyzer::OutputSens).");
    }
    // number of variables (number of species + 2). // T and P
    unsigned int n_vars = r.Mech()->GasMech().SpeciesCount() + 2;

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
                std::string sname = r.Mech()->GasMech().GetSpecies(i)->Name();
                unsigned int len = sname.length();
                fout.write((char*)&len, sizeof(len));
                fout.write(sname.c_str(), len * sizeof(char));
            }

            // output parameter list.
            for (unsigned int i = 0; i < m_sens_params.size(); ++i) {
                fout.write((char*)&m_sens_params.at(i).Type, sizeof(m_sens_params.at(i).Type));
                fout.write((char*)&m_sens_params.at(i).Index, sizeof(m_sens_params.at(i).Index));
            }
        }

        // Flush unwritten buffer to file output.
        fout.flush();
    } else {
        if (m_sens_matrix != NULL) {
            // Write simulation time.
            double time = r.Time();
            fout.write((char*)&time, sizeof(time));
            //cout << "Time : " << r.Time() << std::endl;
            // Write main sensitivity matrix.
            for (unsigned int i = 0; i < m_NS; ++i) {
                double *sdata;
                sdata = NV_DATA_S(m_sens_matrix[i]);
                for (unsigned int j = 0; j < n_vars; ++j) {
                    double val = sdata[j];
                    fout.write((char*)&val, sizeof(val));
                    //cout << val << ", " ;
                }
                //cout << std::endl;
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
    std::string finname = filename + ".sen";
    std::string foutname = filename + "-sensi.csv";

    // Open the sensitivity binary file.
    std::ifstream fin;
    fin.open(finname.c_str(), std::ios_base::in | std::ios_base::binary);

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
    std::vector<std::string> var_names;
    // Parameter list.
    std::vector<SENS_PARAM> arr_params;

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
            std::string sname;
            sname.assign(csname, len);
            var_names.push_back(sname);
            //if (len == 1) delete csname; else delete [] csname;
            delete [] csname;
        }
        var_names.push_back(std::string("T"));
        var_names.push_back(std::string("P"));


        for (unsigned int i = 0; i < NS; ++i) {
            SENS_PARAM arr;
            fin.read(reinterpret_cast<char*>(&arr.Type), sizeof(arr.Type));
            fin.read(reinterpret_cast<char*>(&arr.Index), sizeof(arr.Index));
            arr_params.push_back(arr);
        }
        // Allocating memories for calculated results
        // Total memories require = 2 * n_timesteps * NS * n_vars * sizeof(double);
        //unsigned int mem_space = 2 * n_timesteps * NS * n_vars * sizeof(double);
        //if (mem_space < 209715200) {
        //  // if require mem less than 200MB then will calculate it on RAM.
        double *** avg;
        double *** err_sd;
        double **  temp_sum;
        double **  temp_sum_sqr;
        double *times = new double [n_timesteps];
        
        // Allocate memories for overall result
        avg     = new double ** [n_timesteps];
        err_sd  = new double ** [n_timesteps];
        for (unsigned int i = 0; i < n_timesteps; ++i) {
            avg[i]     = new double * [NS];
            err_sd[i]  = new double * [NS];
            for (unsigned int j = 0; j < NS; ++j) {
                avg[i][j]     = new double [n_vars];
                err_sd[i][j] = new double [n_vars];
                for (unsigned int k = 0; k < n_vars; ++k) {
                    avg[i][j][k]     = 0.0;
                    err_sd[i][j][k]  = 0.0;
                }
            }
        }
        // Allocate memories for temporarily calculations.
        temp_sum     = new double * [NS];
        temp_sum_sqr = new double * [NS];
        for (unsigned int i = 0; i < NS; ++i) {
            temp_sum[i]     = new double [n_vars];
            temp_sum_sqr[i] = new double [n_vars];
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
        std::fstream fout;
        fout.open(foutname.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::trunc);
        if (fout.good()) {
            // Output header line.
            fout << "Time (s), Parameter Type, Parameter Index" ;
            for (unsigned int i = 0; i < var_names.size(); ++i) {
                fout << ", " << var_names.at(i);
                fout << ", Err";
            }
            fout << std::endl;
            for (unsigned int i = 0; i < n_timesteps; ++i) {
                for (unsigned int j = 0; j < NS; ++j) {
                    fout << times[i] << "," << arr_params.at(j).Type << "," << arr_params.at(j).Index;
                    for (unsigned int k = 0; k < n_vars; ++k) {
                        avg[i][j][k]     /= n_runs;
                        err_sd[i][j][k]  /= n_runs;
                        err_sd[i][j][k]  -= avg[i][j][k] * avg[i][j][k];
                        err_sd[i][j][k]  /= n_runs;
                        // fabs might be needed before taking sqrt.
                        err_sd[i][j][k]  = sqrt(err_sd[i][j][k]);
                        fout << "," << avg[i][j][k] << "," << err_sd[i][j][k];
                    }
                    fout << std::endl;
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
}

// Read Sensitivity Matrix block.
// Read block of n x m from fin to matrix and simulation time.
void SensitivityAnalyzer::ReadSensMatrix(std::ifstream &fin, const unsigned int n, const unsigned int m, double &time, double ** matrix, double ** matrix_sql)
{
    fin.read(reinterpret_cast<char*>(&time), sizeof(time));

    for (unsigned int i = 0; i < n; ++i) {
        for (unsigned int j = 0; j < m; ++j) {
            fin.read(reinterpret_cast<char*>(&matrix[i][j]), sizeof(double));
            matrix_sql[i][j] = matrix[i][j] * matrix[i][j];
        }
    }
}

// Parameter pointer to array of double. This is needed by CVODES.
// CVODES will access and change these parameters.
double * SensitivityAnalyzer::ParamsPtr()
{
    return m_params;
}

// Set a pointer to last sensitivity output result.
void SensitivityAnalyzer::SetSensResult(N_Vector *sens_matrix)
{
    m_sens_matrix = sens_matrix;
}


// Set a pointer to last sensitivity output result.
void SensitivityAnalyzer::InitSensMatrix(N_Vector *sens_matrix)
{
    for (unsigned int i = 0; i < m_NS; i++) {
        N_VConst(ZERO, sens_matrix[i]);
    }
    // initialise with 1 if parameter is the same as variable, else zero.
    if (m_probType == Init_Conditions) {
        // number of variables (number of species + 2). // T and P
        unsigned int n_sp = m_mech->GasMech().SpeciesCount();
        unsigned int i_temp = m_reactor->Mixture()->GasPhase().temperatureIndex();
        unsigned int i_dens = m_reactor->Mixture()->GasPhase().densityIndex();
        for (unsigned int i = 0; i < m_sens_params.size(); i++) {
            for (unsigned int j = 0; j < n_sp + 2; j++) {
                if ((m_sens_params.at(i).Type == INIT_T) &&
                    (j == i_temp)) {
                    NV_Ith_S(sens_matrix[i],j) = RCONST(1.0);
                    break;
                } else if ((m_sens_params.at(i).Type == INIT_D) &&
                    (j == i_dens)) {
                    NV_Ith_S(sens_matrix[i],j) = RCONST(1.0);
                    break;
                } else if ((m_sens_params.at(i).Type == INIT_C) &&
                           (m_sens_params.at(i).Index == j)) {
                    NV_Ith_S(sens_matrix[i],j) = RCONST(1.0);
                    break;
                }
            }
        }
    }
}
// Return sensitivity problem type.
SensitivityAnalyzer::SensitivityType SensitivityAnalyzer::ProblemType()
{
    return m_probType;
}

// Parameter scaling factors. approx by original params.
double * SensitivityAnalyzer::ParamBarsPtr()
{
    return m_parambars;
}

bool SensitivityAnalyzer::AddParam(const SENS_PARAM &arrp)
{
    bool paramExist = false;
    if (m_NS == m_sens_params.size()) {
        for (unsigned int i = 0; i < m_NS; i++) {
            if ((m_sens_params.at(i).Index == arrp.Index) &&
                (m_sens_params.at(i).Type == arrp.Type)) {
                paramExist = true;
                break;
            }
        }
        if (!paramExist) {
            m_sens_params.push_back(arrp);
            ++m_NS;
        }
    } else {
        throw std::runtime_error("Number of sensitivity parameters in array and in sensitivity "
                            "object are not the same. (Mops, SensitivityAnalyzer::AddParam).");
    }
    return !paramExist;
}

unsigned int SensitivityAnalyzer::NParams()
{
    return m_NS;
}

