/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the MechParser class declared in the
    swp_mech_parser.h header file.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
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

#include "swp_mech_parser.h"
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_inception.h"
#include "swp_surface_reaction.h"
#include "swp_actsites_reaction.h"
#include "swp_condensation.h"
#include "swp_transcoag.h"
#include "swp_addcoag.h"
#include "swp_abf_model.h"
#include "swp_arssc_reaction.h"
#include "swp_arssc_condensation.h"
#include "swp_diffusion_process.h"
#include "swp_advection_process.h"
#include "swp_molecule_evolution.h"

#include "camxml.h"
#include "string_functions.h"
#include "swp_arssc_model.h"
#include "csv_io.h"
#include <stdexcept>
#include <string>
#include <cassert>
#include <memory>
#include <cstdlib>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;
using namespace Strings;

void MechParser::Read(const std::string &filename, Sweep::Mechanism &mech)
{
    CamXML::Document xml;

    // Try to load XML file into data structure.
    xml.Load(filename);

    // Read the mechanism file version from document, then parse
    // the correct version.
    string version = xml.Root()->GetAttributeValue("version");
    if ((version.compare("1")==0) || (version=="")) {
        return readV1(xml, mech);
    } else {
        // Unknown version.
        throw runtime_error("Unrecognised version attribute in XML file "
                            "(Sweep, MechParser::Read).");
    }
}

void MechParser::readV1(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i, j;

    // Read the particle components and the trackers.
    readComponents(xml, mech);
    readTrackers(xml, mech);

    // READ DEFINED MODELS.

    xml.Root()->GetChildren("model", items);

    for (i=items.begin(); i!=items.end(); ++i) {
        // Check the model type.
        string str = (*i)->GetAttributeValue("type");

        if (str.compare("activesites")==0) {
            // Read active sites model ID.
            str = (*i)->GetAttributeValue("id");

            if (str.compare("haca")==0) {
                // Hydrogen-abstraction, Carbon-addition model.
                mech.AddActSitesModel(ActSites::ABFSites_ID);

                // Initialise the HACA model.
                ActSites::ABFModel::Instance().Initialise(mech);

                // Read the form of alpha.
                str = (*i)->GetAttributeValue("alpha");

                if (str.compare("abf")==0) {
                    // The ABF correlation for alpha should be used.
                    ActSites::ABFModel::Instance().UseAlphaCorrelation();
                } else if (str.compare("profile")==0) {
                    // The profile for alpha which should have been
                    // loaded in the solver should be used.  (Nb, not
                    // sure why alpha profiles are loaded in the solver,
                    // but FlameSolver is definitely where it happens)
                    ActSites::ABFModel::Instance().UseAlphaProfile();
                } else {
                    // Hopefully a numeric value has been supplied, which
                    // shall be used as a constant value for alpha.
                    real alpha = cdble(str);
                    ActSites::ABFModel::Instance().SetAlphaConstant(alpha);
                }
            }
        } else if (str == "particle") {
            // Read particle model ID.
            str = (*i)->GetAttributeValue("id");

            if (str == "surfvol") {
                // Surface-volume model.
                mech.SetAggModel(AggModels::SurfVol_ID);
            } else if (str == "pripartlist") {
                // Basic primary particle model.
                mech.SetAggModel(AggModels::PriPartList_ID);
            } else if (str == "arssc") {
                // ARS-SC model.  Read from XML file.
                CamXML::Element *const ars = xml.Root()->GetFirstChild("arssc");
                if (ars != NULL) {
                    readARSSC_Model(*ars, mech);
                } else {
                    // ARS-SC model parameters not defined!
                    throw runtime_error("ARS-SC model not defined in XML file "
                                        "(Sweep, MechParser::readV1).");
                }

                // Set PAH count tracker.  First search for it in
                // the already defined trackers.  If it doesn't yet
                // exist then add it.
                unsigned int ipah = 0;
                for (ipah=0; ipah!=mech.TrackerCount(); ++ipah) {
                    if (mech.Trackers()[ipah]->Name() == "pahs") {
                        break;
                    }
                }
                if (ipah >= mech.TrackerCount()) {
                    ipah = mech.TrackerCount();
                    mech.AddTracker(*(new Tracker("pahs")));
                }
                SubModels::ARSSC_Model::SetPAH_Tracker(ipah);
            }
        } else if (str == "drag") {
            // Read drag model ID.
            str = (*i)->GetAttributeValue("id");

            if(str == "StokesCunningham") {
                // Parameters of the drag expression
                real A, B, E;

                CamXML::Element* numberXML = (*i)->GetFirstChild("A");
                if(numberXML != NULL) {
                    A = atof(numberXML->Data().c_str());
                }
                else {
                    throw std::runtime_error("Parameter A value must be given for StokesCunningham drag (Sweep, MechParser::readV1)");
                }

                numberXML = (*i)->GetFirstChild("B");
                if(numberXML != NULL) {
                    B = atof(numberXML->Data().c_str());
                }
                else {
                    throw std::runtime_error("Parameter B value must be given for StokesCunningham drag (Sweep, MechParser::readV1)");
                }

                numberXML = (*i)->GetFirstChild("E");
                if(numberXML != NULL) {
                    E = atof(numberXML->Data().c_str());
                }
                else {
                    throw std::runtime_error("Parameter E value must be given for StokesCunningham drag (Sweep, MechParser::readV1)");
                }

                // Set the drag parameters on the mechanism
                mech.SetKnudsenDragConstants(A, B, E);
                mech.SetDragType(Sweep::ParticleModel::KnudsenDrag);

            }
            else if(str == "FreeMol") {
                mech.SetDragType(Sweep::ParticleModel::FreeMolDrag);
            } else if(str == "Temperature") {
                mech.SetDragType(Sweep::ParticleModel::TemperatureDrag);

                // Get constant of proportionality between drag and temperature
                CamXML::Element* numberXML = (*i)->GetFirstChild("A");
                real coeff;
                if(numberXML != NULL) {
                    coeff = atof(numberXML->Data().c_str());
                }
                else {
                    throw std::runtime_error("Parameter A value must be given for drag proprtional to temperature (Sweep, MechParser::readV1)");
                }

                // This option is just for testing so just reuse the support
                // for the Knudsen regime
                mech.SetKnudsenDragConstants(coeff, 0.0, 0.0);


            } else {
                throw std::runtime_error("Unrecognised drag model id (Sweep, MechParser::readV1).");
            }

        } else if (str == "diffusion") {
            // Read diffusion model ID.
            str = (*i)->GetAttributeValue("id");

            if(str == "flamelet") {
                // Diffusion according to soot flamelet equation in flamelet space
                mech.setDiffusionType(Sweep::ParticleModel::FlameletDiffusion);
            }
            else if(str == "einstein") {
                // Diffusion to represent transport in physical space
                // according to Einstein's work
                mech.setDiffusionType(Sweep::ParticleModel::EinsteinDiffusion);
            }
            else {
                throw std::runtime_error("Unrecognised diffusion model id (Sweep, MechParser::readV1).");
            }

        } else if (str == "advection") {
            // Read advection model ID.
            str = (*i)->GetAttributeValue("id");
            if(str == "flamelet") {
                // Advection according to soot flamelet equation in flamelet space
                mech.setAdvectionType(Sweep::ParticleModel::FlameletAdvection);
            }
            else if(str == "physical") {
                // Advection at bulk gas velocity in physical space
                mech.setAdvectionType(Sweep::ParticleModel::BulkAdvection);
            }
            else {
                throw std::runtime_error("Unrecognised advection model id (Sweep, MechParser::readV1).");
            }
        }else if (str == "thermophoresis") {
            // Read advection model ID.
            str = (*i)->GetAttributeValue("id");
            if(str == "waldmann") {
                // Advection according to soot flamelet equation in flamelet space
                mech.setThermophoresisType(Sweep::ParticleModel::WaldmannThermophoresis);
            }
	    else if(str == "none") {
                // No thermophoresis
                mech.setThermophoresisType(Sweep::ParticleModel::NoThermophoresis);
            }
            else {
                throw std::runtime_error("Unrecognised thermophoresis model id (Sweep, MechParser::readV1).");
            }
        } else {
            // An invalid model type has been supplied.
            throw runtime_error("Invalid model type (" + str + ") in XML file "
                                "(Sweep, MechParser::readV1).");
        }
    }

    // READ PARTICLE MODEL.

    xml.Root()->GetChildren("particle", items);
    i = items.begin();

    // Check if the sub-particle tree is active.
    string str = (*i)->GetAttributeValue("subtree");
    if (str == "true") {
        mech.EnableSubPartTree();
    } else {
        mech.DisableSubPartTree();
    }

    // Check the aggregation model.
    str = (*i)->GetAttributeValue("model");
    if (str == "spherical") {
        mech.SetAggModel(AggModels::Spherical_ID);
    } else if (str == "surfvol") {
        mech.SetAggModel(AggModels::SurfVol_ID);
    } else if (str == "pripartlist") {
        mech.SetAggModel(AggModels::PriPartList_ID);
    } else if (str == "PAH") {
        mech.SetAggModel(AggModels::PAH_ID);
        loadPAHStories(*(xml.Root()), mech);
    } else {
        mech.SetAggModel(AggModels::Spherical_ID);
    }

    // Get the sintering model.
    (*i)->GetChildren("sintering", items);
    if (items.size() > 0) {
        i = items.begin();

        // Check if sintering model enabled.
        str = (*i)->GetAttributeValue("enable");
        if (str == "true") {
            mech.SintModel().Enable();
        } else {
            mech.SintModel().Disable();
        }

        // Get sintering model type.
        str = (*i)->GetAttributeValue("model");
        if (str == "viscous_flow") {
            mech.SintModel().SetType(SinteringModel::ViscousFlow);
        } else {
            // Grain-boundary diffusion is the default.
            mech.SintModel().SetType(SinteringModel::GBD);
        }

        // Get pre-exponential factor.
        str = (*i)->GetFirstChild("A")->Data();
        mech.SintModel().SetA(cdble(str));

        // Get characteristic temperature.
        str = (*i)->GetFirstChild("E")->Data();
        mech.SintModel().SetE(cdble(str));

        // Get minimum primary particle diameter.
        // S. Tsantilis, H. Briesen and S.E. Pratsinis,
        // Sintering Time for Silica Particle Growth
        // Aerosol Science and Technology 34: 237�246 (2001).
        if ((*i)->GetFirstChild("Dpmin") != NULL) {
            str = (*i)->GetFirstChild("Dpmin")->Data();
            mech.SintModel().SetDpmin(cdble(str));
        }
    } else {
        // No sintering model defined.
        mech.SintModel().Disable();
    }

    // Read the processes (inceptions, surface reactions and condensations).
    if (mech.AggModel()==AggModels::PAH_ID)
    {
        readPAHInceptions(xml, mech);
    }
    else
    {
        readInceptions(xml, mech);
        readSurfRxns(xml, mech);
        readCondensations(xml, mech);
    }
    readDiffusionProcs(xml, mech);
    readAdvectionProcs(xml, mech);

    // Read which coagulation kernel to use
    readCoagulation(xml, mech);

    // Read ARS-SC model processes.
    if (mech.ContainsSubModel(SubModels::ARSSC_Model_ID)) {
        readARSSC_Inceptions(xml, mech);
        readARSSC_SurfRxns(xml, mech);
        readARSSC_Condensations(xml, mech);
    }
}

// Reads components from a sweep mechanism XML file.
void MechParser::readComponents(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;
    CamXML::Element *el;

    // Get all component specifications from the XML data.
    xml.Root()->GetChildren("component", items);

    // Loop over the components and add them to the mechanism.
    for (i=items.begin(); i!=items.end(); ++i) {
        Component *comp = new Component();

        // Get component ID.
        string str = (*i)->GetAttributeValue("id");
        if (str != "") {
            comp->SetName(str);
        } else {
            // Component must have an id!
            delete comp;
            throw runtime_error("Component statement found with no ""id"" attribute "
                                "(Sweep, MechParser::readComponents).");
        }

        // Get component density.
        el = (*i)->GetFirstChild("density");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetDensity(cdble(str) * 1.0e3); // Convert from g/cm3 to kg/m3.
            } else {
                // Density contains no data.
                std::string msg("Component ");
                msg += comp->Name();
                msg += " density contains no data (Sweep, MechParser::readComponents).";

                delete comp;
                throw runtime_error(msg);
            }
        } else {
            // Component does not have density specified.
            std::string msg("Density required for component ");
            msg += comp->Name();
            msg += " specification (Sweep, MechParser::readComponents).";

            delete comp;
            throw runtime_error(msg);
        }


        // Get coalesc threshold.
        el = (*i)->GetFirstChild("coalthresh");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetCoalescThresh(cdble(str));
            } else {
                // coalthresh contains no data.
                std::string msg("Component ");
                msg += comp->Name();
                msg += " coalthresh contains no data (Sweep, MechParser::readComponents).";

                delete comp;
                throw runtime_error(msg);
            }
        } else {
            comp->SetCoalescThresh(1.0);
        }


        // Get m_growthfact.
        el = (*i)->GetFirstChild("growthfact");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetGrowthFact(cdble(str));
            } else {
                // coalthresh contains no data.
                std::string msg("Component ");
                msg += comp->Name();
                msg += " growthfact contains no data (Sweep, MechParser::readComponents).";

                delete comp;
                throw runtime_error(msg);
            }
        } else {
            comp->SetGrowthFact(1.0);
        }


        // Get minPAH.
        el = (*i)->GetFirstChild("minPAH");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetMinPAH(int(cdble(str)));
            } else {
                // coalthresh contains no data.
                std::string msg("Component ");
                msg += comp->Name();
                msg += " minPAH contains no data (Sweep, MechParser::readComponents).";

                delete comp;
                throw runtime_error(msg);
            }
        } else {
            comp->SetMinPAH(0);
        }


        // Get component mol. wt.
        el = (*i)->GetFirstChild("molwt");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetMolWt(cdble(str)*1.0e-3); // Convert from g/mol to kg/mol.
            } else {
                // Mol. wt. contains no data.
                std::string msg("Component ");
                msg += comp->Name();
                msg += " mol. wt. contains no data (Sweep, MechParser::readComponents).";

                delete comp;
                throw runtime_error(msg);
            }
        } else {
            // Component does not have molecular weight specified.
            std::string msg("Mol. wt. required for component ");
            msg += comp->Name();
            msg += " specification (Sweep, MechParser::readComponents).";

            delete comp;
            throw runtime_error(msg);
        }

		// Get minimum valid value for component
        el = (*i)->GetFirstChild("min");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetMinValid(atof(str.c_str()));
            } else {
                // Mol. wt. contains no data.
                std::string msg("Component ");
                msg += comp->Name();
                msg += " minimum valid value contains no data (Sweep, MechParser::readComponents).";

                delete comp;
                throw runtime_error(msg);
            }
        }

        // Add component to mechanism.
        mech.AddComponent(*comp);
    }
}

// Reads tracker variables from a sweep mechanism XML file.
void MechParser::readTrackers(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;

    // Get all component specifications from the XML data.
    xml.Root()->GetChildren("track", items);

    for (i=items.begin(); i!=items.end(); ++i) {
        Tracker *t = new Tracker();

        string str = (*i)->GetAttributeValue("id");
        if (str != "") {
            t->SetName(str);
        } else {
            // Variable must have a name.
            throw runtime_error("Tracker variable must have an ""id"" attribute "
                                "(Sweep, MechParser::readTrackers).");
        }

        mech.AddTracker(*t);
    }
}

// Reads inception processes from a sweep mechanism XML file.
void MechParser::readInceptions(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    string str;
    unsigned int k = 0;

    // Get list of inceptions from XML data.
    xml.Root()->GetChildren("inception", items);

    for (i=items.begin(),k=0; i!=items.end(); ++i,++k) {
        // Create new inception.
        DimerInception *icn = new DimerInception(mech);
        icn->SetMechanism(mech);
        icn->SetName("Inception " + cstr(k));

        try {
            readInception(*(*i), *icn);
        }
        catch (std::exception &e) {
            delete icn;
            throw;
        }

        // Add inception to mechanism.  Once entered into mechanism, the mechanism
        // takes control of the inception object for memory management.
        mech.AddInception(*icn);
    }
}

// Reads an inception process from an XML element.
void MechParser::readInception(CamXML::Element &xml, Processes::DimerInception &icn)
{
    string str;
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator j;
    int id  = -1;
    real dx = 0.0;

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "") icn.SetName(str);

    // Read reactants.
    readReactants(xml, icn);

    // Get reactant masses and diameters, and set inception
    // parameters.
    fvector mass, diam;
    readReactantMDs(xml, mass, diam);
    if (mass.size() > 1) {
        icn.SetInceptingSpecies(mass[0], mass[1], diam[0], diam[1]);
    } else if (mass.size() > 0) {
        // Dimer inception.
        icn.SetInceptingSpecies(mass[0], mass[0], diam[0], diam[0]);
    }

    // Read products.
    readProducts(xml, icn);

    // Read initial particle composition.
    xml.GetChildren("component", subitems);
    for (j=subitems.begin(); j!=subitems.end(); ++j) {
        // Get component ID.
        str = (*j)->GetAttributeValue("id");
        id = icn.Mechanism()->ComponentIndex(str);

        if (id >= 0) {
            // Get component change.
            str = (*j)->GetAttributeValue("dx");
            dx = cdble(str);
            // Set component change.
            icn.SetParticleComp(id, dx);
        } else {
            // Unknown component in mechanism.
            throw runtime_error(str + ": Component not found in mechanism "
                                "(Sweep, MechParser::readInception).");
        }
    }

    // Read initial tracker variable values.
    xml.GetChildren("track", subitems);
    for (j=subitems.begin(); j!=subitems.end(); j++) {
        // Get tracker ID.
        str = (*j)->GetAttributeValue("id");
        id = icn.Mechanism()->GetTrackerIndex(str);

        if (id >= 0) {
            // Get tracker change.
            str = (*j)->GetAttributeValue("dx");
            dx = cdble(str);
            // Set tracker change.
            icn.SetParticleTracker(id, dx);
        } else {
            // Unknown tracker variable in mechanism.
            throw runtime_error(str + ": Tracker variable not found in mechanism. "
                                "(Sweep, MechParser::readInception).");
        }
    }
}


// Reads inception processes from a sweep mechanism XML file.
void MechParser::readPAHInceptions(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    string str;
    unsigned int k = 0;

    // Get list of inceptions from XML data.
    xml.Root()->GetChildren("pahinception", items);

    for (i=items.begin(),k=0; i!=items.end(); ++i,++k) {
        // Create new inception.
        PAHInception *icn = new PAHInception(mech);
        icn->SetMechanism(mech);
        // This inception only involves one molecule and not two as in the usual inception process
        icn->SetA(1.0);
        icn->SetName("PAHInception " + cstr(k));

        try {
            readPAHInception(*(*i), *icn);
        }
        catch (std::exception &e) {
            delete icn;
            throw;
        }

        // Add inception to mechanism.  Once entered into mechanism, the mechanism
        // takes control of the inception object for memory management.
        mech.AddInception(*icn);
    }
}

// Reads an inception process from an XML element.
void MechParser::readPAHInception(CamXML::Element &xml, Processes::PAHInception &icn)
{
    string str;
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator j;
    int id  = -1;
    real dx = 0.0;

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "") icn.SetName(str);

    // Read reactants.
    //readReactants(xml, icn);

    // Read products.
    //readProducts(xml, icn);

    // Read initial tracker variable values.
    xml.GetChildren("track", subitems);
    for (j=subitems.begin(); j!=subitems.end(); j++) {
        // Get tracker ID.
        str = (*j)->GetAttributeValue("id");
        id = icn.Mechanism()->GetTrackerIndex(str);

        if (id >= 0) {
            // Get tracker change.
            str = (*j)->GetAttributeValue("dx");
            dx = cdble(str);
            // Set tracker change.
            icn.SetParticleTracker(id, dx);
        } else {
            // Unknown tracker variable in mechanism.
            throw runtime_error(str + ": Tracker variable not found in mechanism. "
                                "(Sweep, MechParser::readPAHInception).");
        }
    }
}




// SURFACE REACTIONS.

// Reads surface reactions from a sweep mechanism XML file.
void MechParser::readSurfRxns(CamXML::Document &xml, Mechanism &mech)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;
    string str;
    unsigned int k = 0;

    // Get all surface reactions.
    xml.Root()->GetChildren("reaction", items);

    for (i=items.begin(),k=0; i!=items.end(); ++i,++k) {
        // Get the type of reaction.
        str = (*i)->GetAttributeValue("type");

        SurfaceReaction *rxn = NULL;

        // The type of reaction defined in the XML dictates what type of
        // process class to use.  We now go through the options:
        if (str.compare("surface")==0) {
            // This is a bog-standard surface reaction.
            rxn = new SurfaceReaction(mech);
        } else if (str.compare("abf")==0) {
            // This is an ABF active-sites enabled reaction.
            ActSiteReaction *asrxn = new ActSiteReaction(mech);
            rxn = asrxn;

            // Must also set the reaction to use to ABF model.
            asrxn->SetModel(ActSites::ABFModel::Instance());
        } else {
            // Unrecognised reaction type.
            throw runtime_error("Unrecognised reaction type: " + str +
                                " (Sweep, MechParser::readSurfRxns).");
        }

        // Set default name.
        rxn->SetName("Reaction " + cstr(k));

        // Read the reaction properties.
        try {
            readSurfRxn(*(*i), *rxn);
        } catch (std::exception &e) {
            delete rxn;
            throw;
        }

        // Add process to mechanism.
        mech.AddProcess(*rxn);
    }
}

/*!
 * @param[in]   xml    XML node of type reaction
 * @param[out]  rxn    Reaction on which to set details from XML
 *
 * Rate expressions are of the form
 * \f[
 *    A \left[Y\right] T^n e^\frac{E}{RT} \times particle property,
 * \f]
 * where \f$\left[Y\right]\f$ gives the number of moles of a gas
 * phase species \f$Y\f$ per cubic centimetre and the \f$T^n\f$ term is
 * treated as dimensionless.  This means that for a dimensionless particle
 * property, such as the number of active sites the dimensions of A must
 * be \f$\mathrm{cm}^3 \mathrm{mol}^{-1} \mathrm{s}^{-1}\f$.  If the
 * particle property has dimensions, for example collision diameter
 * raised to the power \f$k\f$, then \f$A\f$ has different units:
 * \f$\mathrm{cm}^{3-n} \mathrm{mol}^{-1} \mathrm{s}^{-1}\f$.  If the
 * particle property is surface or active surface then any power
 * specified for the property in the input is ignored.
 * Active surface is also treated like a dimensionless quantity,
 * because active surface is always multiplied by an active site
 * density to produce a plain number.
 *
 * For historical reasons user inputs of these reaction rate
 * constants must be in units of cm, kcal and K.  However, when
 * the data is read into the program it is immediately converted
 * to SI units as all storage and calculations are performed with
 * SI quantities.
 * @todo Add support for SI inputs.
 *
 * @exception   std::runtime_error  More than one reactant
 * @exception   std::runtime_error  No rate constant
 * @exception   std::runtime_error  Power of particle diameter not specified
 * @exception   std::runtime_error  Particle property on which rate depends not specified
 */
void MechParser::readSurfRxn(CamXML::Element &xml, Processes::SurfaceReaction &rxn)
{
    string str;
    CamXML::Element *el = NULL;

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "") rxn.SetName(str);

    // Is reaction deferred.
    str = xml.GetAttributeValue("defer");
    if (str=="true") {
        rxn.SetDeferred(true);
    } else {
        rxn.SetDeferred(false);
    }

    readReactants(xml, rxn);
    if(rxn.ReactantCount() > 1) {
        std::ostringstream msg;
        msg << "Soot surface reactions may not have more than one reactant, but "
            << rxn.Name() << " has " << rxn.ReactantCount()
            << "MechParser::readSurfRxn";
        throw std::runtime_error(msg.str());
    }

    // Read products.
    readProducts(xml, rxn);

    // Read particle composition change.
    readCompChanges(xml, rxn);

    // Read tracker variable changes.
    readTrackChanges(xml, rxn);

    //========== Read Arrhenius rate parameters ======================
    Sprog::Kinetics::ARRHENIUS arr;
    el = xml.GetFirstChild("A");
    if (el != NULL) {
        arr.A = cdble(el->Data());
    } else {
        // Reaction must have constant.
        throw runtime_error("Surface reaction found with no rate constant "
                            "defined (Sweep, MechParser::readSurfRxns).");
    }
    el = xml.GetFirstChild("n");
    if (el!=NULL) {
        arr.n = cdble(el->Data());
    } else {
        // Default temperature power is 0.
        arr.n = 0.0;
    }
    el = xml.GetFirstChild("E");
    if (el!=NULL) {
        arr.E = cdble(el->Data()) * R / RCAL;
    } else {
        // Default activation energy is zero.
        arr.E = 0.0;
    }

    //========= Particle dependency ==================================
    el = xml.GetFirstChild("particleterm");
    if (el!=NULL) {
        // Get property ID.
        str = el->GetAttributeValue("id");

        // Get power.
        int power = (int)cdble(el->GetAttributeValue("power"));

        if (str.compare("as")==0) {
            // This reaction depends on surface area.  Ignore power,
            // they must have meant 1.
            rxn.SetPropertyID(ParticleCache::iS);

            // Must scale rate constant from cm3 to m3, surface area
            // is multiplied by the site density so that it is a dimensionless
            // quantity hence A has units cm^3 s^-1.
            arr.A *= (1.0e-6);
        }
        else if (str.compare("s")==0) {
            // This reaction depends on surface area.  Ignore power,
            // they must have meant 1.
            rxn.SetPropertyID(ParticleCache::iS);
            arr.A *= (1.0e-2);
        }
        else if (str.compare("d")==0) {
            // This reaction depends on some power of the diameter.
            switch (power) {
                case 1:
                    rxn.SetPropertyID(ParticleCache::iDcol);
                    arr.A *= (1.0e-4);
                    break;
                case 2:
                    rxn.SetPropertyID(ParticleCache::iD2);
                    arr.A *= (1.0e-2);
                    break;
                case -1:
                    rxn.SetPropertyID(ParticleCache::iD_1);
                    arr.A *= (1.0e2);
                    break;
                case -2:
                    rxn.SetPropertyID(ParticleCache::iD_2);
                    arr.A *= (1.0e4);
                    break;
                default:
                    // Oh dear, can't have a zero power.
                    throw runtime_error("""particleterm"" tag found with "
                                        "invalid or zero power attribute "
                                        "(Sweep, MechParser::readSurfRxn)");
            }
        } else {
            // Maybe it is an ARS-SC model property.
            if (rxn.Mechanism()->ContainsSubModel(SubModels::ARSSC_Model_ID)) {
                SubModels::ARSSC_Model::SiteType site =
                    SubModels::ARSSC_Model::IdentifySite(str);
                if (site != SubModels::ARSSC_Model::InvalidSite) {
                    rxn.SetPropertyID((unsigned int)site, SubModels::ARSSC_Model_ID);
                }
            }
            arr.A *= (1.0e-6);
        }
    } else {
        throw runtime_error("Surface process defined without ""particleterm"" "
                            "element (Sweep, MechParser::readSurfRxn).");
    }

    rxn.SetArrhenius(arr);

}


// CONDENSATIONS.

// Reads condensation processes from a sweep mechanism XML file.
void MechParser::readCondensations(CamXML::Document &xml, Mechanism &mech)
{
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    string str;
    unsigned int k = 0;

    // Get all surface reactions.
    xml.Root()->GetChildren("condensation", items);

    for (i=items.begin(),k=0; i!=items.end(); ++i,++k) {
        // Create a new condensation object.
        Condensation *cond = new Condensation(mech);

        // Set default name.
        cond->SetName("Condensation " + cstr(k));

        // Read the condensation properties.
        try {
            readCondensation(*(*i), *cond);
        } catch (std::exception &e) {
            delete cond;
            throw;
        }

        // Add condensation to mechanism.
        mech.AddProcess(*cond);
    }
}

// Reads a condensation process from a sweep mechanism XML file.
void MechParser::readCondensation(CamXML::Element &xml, Processes::Condensation &cond)
{
    string str;

    // Is condensation deferred.
    str = xml.GetAttributeValue("defer");
    if (str=="true") {
        cond.SetDeferred(true);
    } else {
        cond.SetDeferred(false);
    }

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "") {
        cond.SetName(str);
    }

    // Read reactants.
    try {
        // Get reactant stoichiometry.
        readReactants(xml, cond);

        // Get reactant masses and diameters, and set inception
        // parameters.
        fvector mass, diam;
        readReactantMDs(xml, mass, diam);
        if (mass.size() > 0) {
            cond.SetCondensingSpecies(mass[0], diam[0]);
        }
    } catch (std::exception &e) {
        throw;
    }

    // Read products.
    readProducts(xml, cond);

    // Read particle composition change.
    readCompChanges(xml, cond);

    // Read tracker variable changes.
    readTrackChanges(xml, cond);

    // The condensation should change something, otherwise what is the point of having it?
    assert((cond.CompChange().size() > 0) || (cond.TrackChange().size() > 0));

    // Read Arrhenius rate parameters.
    real A = 0.0;
    CamXML::Element *el = xml.GetFirstChild("A");
    if (el != NULL) {
        A = cdble(el->Data());
    } else {
        A = 1.0;
    }
    cond.SetA(A);
}

//COAGULATION

/**
 *
 */
void MechParser::readCoagulation(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items;

    // Get list of inceptions from XML data.
    xml.Root()->GetChildren("coagulation", items);

    if(items.size() == 0)
    {
        // Use the default (transition regime) kernel
        mech.AddCoagulation(*(new Processes::TransitionCoagulation(mech)));
    }
    else if(items.size() == 1)
    {
        // Read the user choice of kernel.  If no kernel is specified
        // a default will be used
        const CamXML::Element* const kernel = (items.front())->GetFirstChild("kernel");
        if(kernel ==  NULL)
        {
            throw std::runtime_error("No kernel given for coagulation \
                                     (Sweep, MechParser::readCoagulation)");
        }
        else
        {
            // Work out which kernel to use and create it
            const string kernelName = kernel->Data();

            std::auto_ptr<Processes::Coagulation> coag;
            if(kernelName == "transition")
                coag.reset(new Processes::TransitionCoagulation(mech));
            else if(kernelName == "additive")
                coag.reset(new Processes::AdditiveCoagulation(mech));
            else
                // Unrecognised option
                throw std::logic_error("Coagulation kernel " + kernelName + "not yet available \
                                        (Sweep, MechParser::readCoagulation)");

            // Rate scaling now that a process has been created
            real A = 0.0;
            CamXML::Element *el = (items.front())->GetFirstChild("A");
            if (el != NULL) {
                A = cdble(el->Data());
            } else {
                A = 1.0;
            }
            coag->SetA(A);

            mech.AddCoagulation(*coag);

            // Get rid of the auto_ptr without deleting the coagulation object
            coag.release();
        }
    }
    else
    {
        throw std::runtime_error("More than one coagulation process specified \
                                 (Sweep, MechParser::readCoagulation)");
    }
}

// TRANSPORT PROCESSES

/*!
 *\param[in]        xml         XML document containing a mechanism description
 *\param[inout]     mech        Mechanism to which transport processes will be added
 */
void MechParser::readDiffusionProcs(CamXML::Document &xml, Mechanism &mech)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;
    string str;
    unsigned int k = 0;

    // Get all transport processes.
    xml.Root()->GetChildren("diffusion", items);

    for (i=items.begin(),k=0; i!=items.end(); ++i,++k) {
        DiffusionProcess *tran = new DiffusionProcess();
        // Set default name.
        tran->SetName("Diffusion " + cstr(k));

        // Read the reaction properties.
        try {
            readDiffusionProc(*(*i), *tran);
        } catch (std::exception &e) {
            delete tran;
            throw;
        }

        // Add process to mechanism.
        mech.AddTransport(*tran);
    }
}

/*!
 *\param[in]        xml     XML node of type diffusion
 *\param[in]        rxn     mechanism to which to add the process
 */
void MechParser::readDiffusionProc(CamXML::Element &xml, Processes::DiffusionProcess &tran)
{
    string str;
    CamXML::Element *el = NULL;

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "")
        tran.SetName(str);

    // Particle dependency.
    el = xml.GetFirstChild("particleterm");
    if (el!=NULL) {
        // Get property ID.
        str = el->GetAttributeValue("id");

        // Get power.
        int id = atoi(el->GetAttributeValue("power").c_str());

        if (str.compare("d")==0) {
            // This reaction depends on some power of the collision diameter.
            switch (id) {
                case 0:
                    tran.SetPropertyID(ParticleCache::iUniform);
                    break;
                case 1:
                    tran.SetPropertyID(ParticleCache::iDcol);
                    break;
                case 2:
                    tran.SetPropertyID(ParticleCache::iD2);
                    break;
                case -1:
                    tran.SetPropertyID(ParticleCache::iD_1);
                    break;
                case -2:
                    tran.SetPropertyID(ParticleCache::iD_2);
                    break;
                default:
                    // Oh dear, can't have a zero power.
                    throw runtime_error("particleterm"" tag found with \
                                         invalid power attribute \
                                        (Sweep, MechParser::readTransportProc)");
            }
        }
        else {
            throw runtime_error("Surface process defined without particleterm \
                                 element (Sweep, MechParser::readTransportProc).");
        }
    }

    // Temperature dependency.
    el = xml.GetFirstChild("temperature");
    real power = 0;
    if (el!=NULL) {
        // Get power.
        power = atof(el->GetAttributeValue("power").c_str());
    }
    tran.SetTemperatureExponent(power);

    // Read scaling factor
    real A = 0.0;
    el = xml.GetFirstChild("A");
    if (el != NULL) {
        A = atof(el->Data().c_str());
    } else {
        A = 1.0;
    }
    tran.SetA(A);
}

/*!
 *\param[in]        xml         XML document containing a mechanism description
 *\param[inout]     mech        Mechanism to which transport processes will be added
 */
void MechParser::readAdvectionProcs(CamXML::Document &xml, Mechanism &mech)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;
    string str;
    unsigned int k = 0;

    // Get all transport processes.
    xml.Root()->GetChildren("advection", items);

    for (i=items.begin(),k=0; i!=items.end(); ++i,++k) {
        AdvectionProcess *tran = new AdvectionProcess();
        // Set default name.
        tran->SetName("Advection " + cstr(k));

        // Read the reaction properties.
        try {
            readAdvectionProc(*(*i), *tran);
        } catch (std::exception &e) {
            delete tran;
            throw;
        }

        // Add process to mechanism.
        mech.AddTransport(*tran);
    }
}

/*!
 *\param[in]        xml     XML node of type advection
 *\param[in]        rxn     mechanism to which to add the process
 */
void MechParser::readAdvectionProc(CamXML::Element &xml, Processes::AdvectionProcess &tran)
{
    string str;
    CamXML::Element *el = NULL;

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "")
        tran.SetName(str);

    // Read scaling factor
    real A = 0.0;
    el = xml.GetFirstChild("A");
    if (el != NULL) {
        A = atof(el->Data().c_str());
    } else {
        A = 1.0;
    }
    tran.SetA(A);
}


// REACTION SHARED COMPONENTS.

// Reads reactants into a process.
void MechParser::readReactants(CamXML::Element &xml, Process &proc)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;

    // Get the reactant.
    xml.GetChildren("reactant", items);

    for (i=items.begin(); i!=items.end(); ++i) {
        // Get the species ID.
        string str = (*i)->GetAttributeValue("id");
        int isp = Sprog::Species::Find(str, *proc.Mechanism()->Species());

        if (isp < 0) {
            // Reactant does not exist in list of species.
             throw runtime_error(str + ": Reactant species not "
                                 "found in mechanism (Sweep, "
                                 "MechParser::readReactants).");
        } else {
            // Get reactant stoichiometry.
            str = (*i)->GetAttributeValue("stoich");
            int mu = (int)cdble(str);

            // Add reactant to inception.
            proc.AddReactant(isp, mu);
        }
    }
}

// Reads products into a process.
void MechParser::readProducts(CamXML::Element &xml, Process &proc)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;

    // Get the products.
    xml.GetChildren("product", items);

    for (i=items.begin(); i!=items.end(); ++i) {
        // Get the species ID.
        string str = (*i)->GetAttributeValue("id");
        int isp = Sprog::Species::Find(str, *proc.Mechanism()->Species());

        if (isp < 0) {
            // Product does not exist in list of species.
            throw runtime_error(str + ": Product species not "
                                "found in mechanism (Sweep, "
                                "MechParser::readProducts).");
        } else {
            // Get product stoichiometry.
            str = (*i)->GetAttributeValue("stoich");
            int mu = (int)cdble(str);

            // Add product to inception.
            proc.AddProduct(isp, mu);
        }
    }
}

// Reads reactant masses and diameters.  This is required by inceptions
// and condensation to calculate collision rates.
void MechParser::readReactantMDs(CamXML::Element &xml,
                                 fvector &mass, fvector &diam)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;

    // Get reactants.
    xml.GetChildren("reactant", items);

    // Clear output vectors.
    mass.clear();
    diam.clear();

    for (i=items.begin(); i!=items.end(); ++i) {
        // Get species mass
        string str = (*i)->GetAttributeValue("m");
        real m = cdble(str);
        if (m > 0.0) {
            mass.push_back(m*1.0e-3); // Convert from g to kg.
        } else {
            // Species can't have zero mass!
            throw runtime_error("Species can't have zero mass "
                                "(Sweep, MechParser::readReactantMDs).");
        }

        // Get species diameter
        str = (*i)->GetAttributeValue("d");
        real d = cdble(str);
        if (d > 0.0) {
            diam.push_back(d*1.0e-2); // Convert from cm to m.
        } else {
            // Species can't have zero diameter!
            throw runtime_error("Species can't have zero diameter "
                                "(Sweep, MechParser::readReactantMDs).");
        }
    }
}

// Reads composition changes into a particle process.
void MechParser::readCompChanges(CamXML::Element &xml, ParticleProcess &proc)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;

    // Get list of component changes from XML.
    xml.GetChildren("component", items);

    for (i=items.begin(); i!=items.end(); ++i) {
        // Get component ID.
        string str = (*i)->GetAttributeValue("id");
        int id = proc.Mechanism()->ComponentIndex(str);

        if (id >= 0) {
            // Get component change.
            str = (*i)->GetAttributeValue("dx");
            real dx = cdble(str);
            // Set component change.
            proc.SetCompChange(id, dx);
        } else {
            throw runtime_error(str + ": Component not found in mechanism "
                                "Sweep, MechParser::readCompChanges).");
        }
    }
}

// Reads tracker variable changes into a particle process.
void MechParser::readTrackChanges(CamXML::Element &xml, ParticleProcess &proc)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;

    // Get list of component changes from XML.
    xml.GetChildren("track", items);

    for (i=items.begin(); i!=items.end(); ++i) {
        // Get component ID.
        string str = (*i)->GetAttributeValue("id");
        int id = proc.Mechanism()->GetTrackerIndex(str);

        if (id >= 0) {
            // Get component change.
            str = (*i)->GetAttributeValue("dx");
            real dx = cdble(str);
            // Set component change.
            proc.SetTrackChange(id, dx);
        } else {
            throw runtime_error(str + ": Tracker variable not found in mechanism "
                                "Sweep, MechParser::readTrackChanges).");
        }
    }
}

// Reads a maths functional from XML and creates an object to contain it.
Sweep::Maths::Functional *const MechParser::readFunctional(CamXML::Element &xml)
{
    Sweep::Maths::Functional *fun = NULL;
    CamXML::Element *el = NULL;

    // Read type.
    string str = xml.GetAttributeValue("type");
    if (str == "linear") {
        // A linear function.
        fun = new Sweep::Maths::Linear();
        // Read two parameters.
        el = xml.GetFirstChild("A");
        if (el) {
            str = el->Data();
            fun->SetParam(0, cdble(str));
        }
        el = xml.GetFirstChild("B");
        if (el) {
            str = el->Data();
            fun->SetParam(1, cdble(str));
        }
    } else {
        throw runtime_error(str + ": Invalid functional type "
                            "(Sweep, MechParser::readFunctional).");
    }

    // Return the new functional object.
    return fun;
}


// ARS-SC MODEL.

// Reads ARS-SC model static parameters into a mechanism.
void MechParser::readARSSC_Model(CamXML::Element &xml, Mechanism &mech)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;
    CamXML::Element *el = NULL;
    string str;
    SubModels::ARSSC_Model::SiteType parent, var;
    Sweep::Maths::Functional *fun = NULL;
    real wt;

    // Add ARS-SC model to mechanism.
    mech.AddSubModel(SubModels::ARSSC_Model_ID);

    // Read combined-site correlations.
    xml.GetChildren("combinedsite", items);
    for (i=items.begin(); i!=items.end(); ++i) {
        // Get combined-site ID and type.
        str = (*i)->GetAttributeValue("id");
        if (str == "r6") {
            parent = SubModels::ARSSC_Model::R6;
        } else if (str == "ac_r6") {
            parent = SubModels::ARSSC_Model::ACR6;
        } else if (str == "r5_ed") {
            parent = SubModels::ARSSC_Model::R5ED;
        } else if (str == "r5_ac") {
            parent = SubModels::ARSSC_Model::R5AC;
        } else {
            // Invalid site ID.
            throw runtime_error(str + ": Invalid ARS-SC combined-site ID "
                                "(Sweep, MechParser::readARSSC_Model).");
        }

        // Get dependent variable for correlation.
        el = (*i)->GetFirstChild("var");
        if (el != NULL) {
            str = el->Data();
            if (str == "edges") {
                var = SubModels::ARSSC_Model::FreeEdge;
            } else if (str == "armchairs") {
                parent = SubModels::ARSSC_Model::Armchair;
            } else if (str == "zigzags") {
                parent = SubModels::ARSSC_Model::Zigzag;
            } else if (str == "zigR5s") {
                parent = SubModels::ARSSC_Model::R5;
            } else if (str == "bays") {
                parent = SubModels::ARSSC_Model::Bay;
            } else {
                // Invalid site ID.
                throw runtime_error(str + ": Invalid ARS-SC site ID for dependent "
                                    "variable of combined-site correlation "
                                    "(Sweep, MechParser::readARSSC_Model).");
            }
        } else {
            // Dependent variable must be defined..
            throw runtime_error(str + ": No dependent variable defined for this "
                                "combined-site ID "
                                "(Sweep, MechParser::readARSSC_Model).");
        }

        // Get functional form.
        el  = (*i)->GetFirstChild("functional");
        fun = readFunctional(*el);
        if (fun == NULL) {
            // Functional form undefined.
            throw runtime_error("Functional form of combined-site correlation undefined "
                                "(Sweep, MechParser::readARSSC_Model).");
        }

        // Set combined-site correlation in model.
        switch(parent) {
            case SubModels::ARSSC_Model::R6:
                SubModels::ARSSC_Model::SetR6_Variable(var);
                SubModels::ARSSC_Model::SetR6_Functional(*fun);
                break;
            case SubModels::ARSSC_Model::ACR6:
                SubModels::ARSSC_Model::SetACR6_Variable(var);
                SubModels::ARSSC_Model::SetACR6_Functional(*fun);
                break;
            case SubModels::ARSSC_Model::R5ED:
                SubModels::ARSSC_Model::SetR5ED_Variable(var);
                SubModels::ARSSC_Model::SetR5ED_Functional(*fun);
                break;
            case SubModels::ARSSC_Model::R5AC:
                SubModels::ARSSC_Model::SetR5AC_Variable(var);
                SubModels::ARSSC_Model::SetR5AC_Functional(*fun);
                break;
            // Do nothing in all other cases
            default:
                throw std::runtime_error("Unhandled case in switch statement (MechParser::readARSSC_Model)");
                break;
        }
    }

    // Read neighbour weights.
    xml.GetChildren("neighbours", items);
    for (i=items.begin(); i!=items.end(); ++i) {
        // Get parent site ID and type.
        str = (*i)->GetAttributeValue("site");
        if (str == "edges") {
            parent = SubModels::ARSSC_Model::FreeEdge;
        } else if (str == "armchairs") {
            parent = SubModels::ARSSC_Model::Armchair;
        } else if (str == "zigzags") {
            parent = SubModels::ARSSC_Model::Zigzag;
        } else if (str == "bays") {
            parent = SubModels::ARSSC_Model::Bay;
        } else {
            // Invalid site ID.
            throw runtime_error(str + ": Invalid ARS-SC parent site ID for neighbour "
                                "weights (Sweep, MechParser::readARSSC_Model).");
        }

        // Get probability for free-edges.
        el = (*i)->GetFirstChild("edges");
        if (el != NULL) {
            wt = cdble(el->Data());
            SubModels::ARSSC_Model::SetNeighbour(
                parent, SubModels::ARSSC_Model::FreeEdge, wt);
        } else {
            // Error!  Weight for free-edges must be defined
            throw runtime_error(str + ": No neighbour weight defined for free-edges"
                                "for this site (Sweep, MechParser::readARSSC_Model).");
        }

        // Get probability for armchairs.
        el = (*i)->GetFirstChild("armchairs");
        if (el != NULL) {
            wt = cdble(el->Data());
            SubModels::ARSSC_Model::SetNeighbour(
                parent, SubModels::ARSSC_Model::Armchair, wt);
        } else {
            // Error!  Weight for armchairs must be defined
            throw runtime_error(str + ": No neighbour weight defined for armchairs"
                                "for this site (Sweep, MechParser::readARSSC_Model).");
        }

        // Get probability for zig-zags.
        el = (*i)->GetFirstChild("zigzags");
        if (el != NULL) {
            wt = cdble(el->Data());
            SubModels::ARSSC_Model::SetNeighbour(
                parent, SubModels::ARSSC_Model::Zigzag, wt);
        } else {
            // Error!  Weight for zig-zags must be defined
            throw runtime_error(str + ": No neighbour weight defined for zig-zags"
                                "for this site (Sweep, MechParser::readARSSC_Model).");
        }

        // Get probability for bays.
        el = (*i)->GetFirstChild("bays");
        if (el != NULL) {
            wt = cdble(el->Data());
            SubModels::ARSSC_Model::SetNeighbour(
                parent, SubModels::ARSSC_Model::Bay, wt);
        } else {
            // Error!  Weight for bays must be defined
            throw runtime_error(str + ": No neighbour weight defined for bays "
                                "for this site (Sweep, MechParser::readARSSC_Model).");
        }
    }
}

// Reads ARS-SC inception processes from a sweep mechanism XML file.
void MechParser::readARSSC_Inceptions(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    string str;

    // Get list of inceptions from XML data.
    xml.Root()->GetChildren("arsinception", items);

    unsigned int k = 0;
    for (i=items.begin(); i!=items.end(); ++i,++k) {
        // Create new inception.
        ARSSC_Inception *icn = new ARSSC_Inception(mech);
        icn->SetMechanism(mech);
        icn->SetName("ARS-SC Inception " + cstr(k));

        // Read basic inception properties.
        try {
            readInception(*(*i), *icn);
        } catch (std::exception &e) {
            delete icn;
            throw;
        }

        // Read ARS-SC inception parameters.
        try {
            readARSSC_Sites(*(*i), *icn);
        } catch (std::exception &e) {
            delete icn;
            throw;
        }

        // Add inception to mechanism.  Once entered into mechanism, the mechanism
        // takes control of the inception object for memory management.
        mech.AddInception(*icn);
    }
}

// Reads ARS-SC surface reactions from a sweep mechanism XML file.
void MechParser::readARSSC_SurfRxns(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    string str;
    unsigned int k = 0;

    // Get list of reactions from XML data.
    xml.Root()->GetChildren("arsreaction", items);

    for (i=items.begin(),k=0; i!=items.end(); ++i,++k) {
        // Create new reaction.
        ARSSC_Reaction *rxn = new ARSSC_Reaction(mech);
        rxn->SetMechanism(mech);
        rxn->SetName("ARS-SC Reaction " + cstr(k));

        // Read basic reaction properties.
        try {
            readSurfRxn(*(*i), *rxn);
        } catch (std::exception &e) {
            delete rxn;
            throw;
        }

        // Read ARS-SC parameters.
        try {
            readARSSC_Sites(*(*i), *rxn);
        } catch (std::exception &e) {
            delete rxn;
            throw;
        }

        // Add reaction to mechanism.  Once entered into mechanism, the mechanism
        // takes control of the object for memory management.
        mech.AddProcess(*rxn);
    }
}

// Reads ARS-SC condensations from a sweep mechanism XML file.
void MechParser::readARSSC_Condensations(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    string str;
    unsigned int k = 0;

    // Get list of reactions from XML data.
    xml.Root()->GetChildren("arscondensation", items);

    for (i=items.begin(),k=0; i!=items.end(); ++i,++k) {
        // Create new reaction.
        ARSSC_Condensation *cond = new ARSSC_Condensation(mech);
        cond->SetMechanism(mech);
        cond->SetName("ARS-SC Condensation " + cstr(k));

        // Read basic reaction properties.
        try {
            readCondensation(*(*i), *cond);
        } catch (std::exception &e) {
            delete cond;
            throw;
        }

        // Read ARS-SC parameters.
        try {
            readARSSC_Sites(*(*i), *cond);
        } catch (std::exception &e) {
            delete cond;
            throw;
        }

        // Add reaction to mechanism.  Once entered into mechanism, the mechanism
        // takes control of the object for memory management.
        mech.AddProcess(*cond);
    }
}

// Reads ARS-SC site parameters for a process into
// the given ARSSC_Process object.
void MechParser::readARSSC_Sites(CamXML::Element &xml, Processes::ARSSC_Process &proc)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;
    string str;
    SubModels::ARSSC_Model::SiteType id;
    real dx = 0.0;

    // Read site count changes.
    xml.GetChildren("sites", items);
    for (i=items.begin(); i!=items.end(); ++i) {
        // Get site type.
        str = (*i)->GetAttributeValue("id");
        id = SubModels::ARSSC_Model::IdentifySite(str);

        if (id != SubModels::ARSSC_Model::InvalidSite) {
            // Get site count.
            str = (*i)->GetAttributeValue("dx");
            dx = cdble(str);

            // Set site count in inception object.
            proc.SetSiteCount(id, dx);
        } else {
            // Error! Invalid site ID.
            throw runtime_error(str + ": Invalid site ID in ARS process definition "
                                "(Sweep, MechParser::readARSSC_Sites).");
        }
    }

    // Read site updates.  There can be only one tag with this
    // ID defined for each process.
    CamXML::Element *el = xml.GetFirstChild("updatesites");
    if (el != NULL) {
        // Get parent site type.
        str = el->GetAttributeValue("parent");
        id = SubModels::ARSSC_Model::IdentifySite(str);

        // Get update quantity.
        str = el->GetAttributeValue("count");
        dx = cdble(str);

        // Set up the ARS-SC process object.
        proc.SetParentSite(id);
        proc.SetUpdateCount((int)(dx+0.5));
        if (id == SubModels::ARSSC_Model::InvalidSite) {
            proc.EnableCustomWts();
        } else {
            proc.EnableParentWts();
        }
    }
}

/*!
 * @param[in]   xml     XML node with one or more <pahfile> children
 *
 * @exception   std::runtime_error  No pahfile child elements found
 * @exception   std::runtime_error  Neither time nor position specified for a file
 */
void MechParser::loadPAHStories(CamXML::Element &xml, Mechanism &mech) {
    // Extract the details of the pah story files
    std::vector<CamXML::Element*> items;
    xml.GetChildren("pahfile", items);

    std::vector<CamXML::Element*>::iterator it = items.begin();
    const std::vector<CamXML::Element*>::iterator itEnd = items.end();

    if(it == itEnd)
        throw std::runtime_error("No pahfile elements found in mechanism (MechParser::loadPAHStories)");

    MoleculeEvolution::Database db;

    while(it != itEnd) {
        std::string fileName = (*it)->GetAttributeValue("path");

        // See if a creation time has been specified for the molecules in this file
        real time = 0.0;
        const CamXML::Attribute *timeAttr = (*it)->GetAttribute("time");

        if(timeAttr != NULL)
            time = std::atof(timeAttr->GetValue().c_str());

        // See if a creation position has been specified for the molecules in this file
        real position = 0.0;
        const CamXML::Attribute *posnAttr = (*it)->GetAttribute("position");

        if(posnAttr != NULL)
            position = std::atof(posnAttr->GetValue().c_str());

        if((timeAttr == NULL) && (posnAttr == NULL)) {
            throw std::runtime_error("At least one of time and position must be specified for " +
                                      fileName + " (MechParser::loadPAHStories)");
        }

        db.addStoriesFromFile(fileName, time, position);

        // Move on to next file
        ++it;
    }

    mech.setMoleculeStories(db);
}
