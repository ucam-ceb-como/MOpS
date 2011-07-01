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
#include "swp_constant_inception.h"
#include "swp_surface_reaction.h"
#include "swp_actsites_reaction.h"
#include "swp_condensation.h"
#include "swp_transcoag.h"
#include "swp_addcoag.h"
#include "swp_constcoag.h"
#include "swp_weighted_addcoag.h"
#include "swp_weighted_constcoag.h"
#include "swp_coag_weight_rules.h"
#include "swp_abf_model.h"
#include "swp_diffusion_process.h"
#include "swp_advection_process.h"

#include "camxml.h"
#include "string_functions.h"
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

//! Anonymous namespace for functions only used within this file
namespace {

/*!
 * Read and set initial particle composition for an inception
 *
 *@param[in]		xml		XML node with component children
 *@param[in,out]	icn		Inception instance for which to set composition of newly incepted particles
 */
void readInceptedComposition(const CamXML::Element &xml, Sweep::Processes::Inception &icn) {
    // Get the component (= composition) XML
	vector<CamXML::Element*> subitems;
    xml.GetChildren("component", subitems);

	for (vector<CamXML::Element*>::iterator j=subitems.begin(); j!=subitems.end(); ++j) {
		// Get component ID.
		std::string str = (*j)->GetAttributeValue("id");
		int id = icn.Mechanism()->ComponentIndex(str);

		if (id >= 0) {
			// Get component change.
			str = (*j)->GetAttributeValue("dx");
			Sweep::real dx = cdble(str);
			// Set component change.
			icn.SetParticleComp(id, dx);
		} else {
			// Unknown component in mechanism.
			throw runtime_error(str + ": Component not found in mechanism (loadInceptedComposition)");
		}
	}
}

/*!
 * Read and set initial particle tracker values for an inception
 *
 *@param[in]		xml		XML node with track children
 *@param[in,out]	icn		Inception instance for which to set tracker values of newly incepted particles
 */
void readInceptedTrackers(const CamXML::Element &xml, Sweep::Processes::Inception &icn) {
    // Get the component (= composition) XML
	vector<CamXML::Element*> subitems;
	xml.GetChildren("track", subitems);

	for (vector<CamXML::Element*>::iterator j=subitems.begin(); j!=subitems.end(); j++) {
		// Get tracker ID.
		std::string str = (*j)->GetAttributeValue("id");
		int id = icn.Mechanism()->GetTrackerIndex(str);

		if (id >= 0) {
			// Get tracker change.
			str = (*j)->GetAttributeValue("dx");
			Sweep::real dx = cdble(str);
			// Set tracker change.
			icn.SetParticleTracker(id, dx);
		} else {
			// Unknown tracker variable in mechanism.
			throw runtime_error(str + ": Tracker variable not found in mechanism (loadInceptedTrackers)");
		}
	}
}
} // anonymous namespace

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
            throw std::runtime_error("<model id=\"particle\" is not a valid tag (Sweep::MechParser::readV1");
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

    const CamXML::Element* particleXML = xml.Root()->GetFirstChild("particle");

    // Check if the sub-particle tree is active.
    string str = particleXML->GetAttributeValue("subtree");
    if (str == "true") {
        throw std::runtime_error("Subtrees are no longer supported (Sweep::MechParser::readV1)");
    }

    // Check the aggregation model.
    str = particleXML->GetAttributeValue("model");
    if (str == "spherical") {
        mech.SetAggModel(AggModels::Spherical_ID);
    } else if (str == "surfvol") {
        mech.SetAggModel(AggModels::SurfVol_ID);
    } else if (str == "PAH") {
	// Reject all old style input files
		throw std::runtime_error("PAH-PP MODEL are no longer supported (Sweep::MechParser::readV1), you can use NEW PAH_KMC model");
    } else if (str == "PAH_KMC") {
        mech.SetAggModel(AggModels::PAH_KMC_ID);
	}else {
        mech.SetAggModel(AggModels::Spherical_ID);
    }

    // See if there are any secondary particle criteria
    const CamXML::Element* secondaryXML = particleXML->GetFirstChild("secondaryparticle");
    if(secondaryXML != NULL) {
        throw std::runtime_error("Secondary particles are no longer supported (Sweep::MechParser::readV1)");
    }

    // Get the sintering model.
    particleXML->GetChildren("sintering", items);
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
    if (mech.AggModel()==AggModels::PAH_KMC_ID)
    {
        readPAHInceptions(xml, mech);
    }
    else
    {
        readInceptions(xml, mech);
        readConstantInceptions(xml, mech);
        readSurfRxns(xml, mech);
        readCondensations(xml, mech);
    }
    readDiffusionProcs(xml, mech);
    readAdvectionProcs(xml, mech);

    // Read which coagulation kernel to use
    readCoagulation(xml, mech);
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

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "") icn.SetName(str);

    // Read reactants.
    readReactants(xml, icn);

    // Find the rate calculation method (a coagulation kernel)
    // Currently the only possibilities are free molecular and transition
    // regime kernels.  These are handled by a boolean flag.  An enum
    // will be needed if more cases are introduced (or possibly a pointer
    // to an appropriate member function of DimerInception.)
    bool useFreeMolRegime = false;
    str = xml.GetAttributeValue("rate");
    if(!str.empty()) {
        if(str ==  "freemolecular")
            useFreeMolRegime = true;
        else if(str == "transition")
            useFreeMolRegime = false;
        else
            throw std::runtime_error("Unrecognised rate type " + str + " in Sweep::MechParser::readInception");
    }

    // Get reactant masses and diameters, and set inception
    // parameters.
    fvector mass, diam;
    readReactantMDs(xml, mass, diam);
    if (mass.size() == 2) {
        if(useFreeMolRegime)
            icn.SetInceptingSpeciesFreeMol(mass[0], mass[1], diam[0], diam[1]);
        else
            icn.SetInceptingSpecies(mass[0], mass[1], diam[0], diam[1]);
    }
    else if (mass.size() == 1) {
        if(useFreeMolRegime)
            icn.SetInceptingSpeciesFreeMol(mass[0], mass[0], diam[0], diam[0]);
        else
            icn.SetInceptingSpecies(mass[0], mass[0], diam[0], diam[0]);
    }
    else
        throw std::runtime_error("One or two inception species must be specified in Sweep::MechParser::readInception");


    // Rate scaling now that a process has been created
    real A = 0.0;
    CamXML::Element *el = xml.GetFirstChild("A");
    if (el != NULL) {
        A = cdble(el->Data());
        icn.SetA(A);
    }

    // Read products.
    readProducts(xml, icn);

    // Read initial particle composition.
    readInceptedComposition(xml, icn);

    // Read initial tracker variable values.
    readInceptedTrackers(xml, icn);


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
        // icn->SetA(1.0); this is now done in the constructors, the old practice of defaulting
        // to A=0.5 has now been removed. riap 28Apr2011

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

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "") icn.SetName(str);

    // Read reactants.
    //readReactants(xml, icn);

    // Read products.
    //readProducts(xml, icn);

    // Read initial tracker values
    readInceptedTrackers(xml, icn);

    CamXML::Element* secondary = xml.GetFirstChild("usesecondary");
    if(secondary != NULL)
        throw std::runtime_error("Secondary particles are no longer supported (MechParser::readPAHInception)");
}

/*!
 * Read constant rate inception processes from a sweep mechanism XML file.
 *
 *@param[in]		xml		XML node with constantinception children
 */
void MechParser::readConstantInceptions(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    string str;
    unsigned int k = 0;

    // Get list of inceptions from XML data.
    xml.Root()->GetChildren("constantinception", items);

    for (i=items.begin(),k=0; i!=items.end(); ++i,++k) {
        // Create new inception.
        ConstantInception *icn = new ConstantInception(mech);
        icn->SetMechanism(mech);
        icn->SetName("ConstantInception " + cstr(k));

        try {
            readConstantInception(*(*i), *icn);
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

/*!
 *@param[in]		xml		XML node of type constantinception
 *@param[in,out]	icn		Instance of ConstantInception on which to set data
 *
 *@exception	std::runtime_error		No rate specified
 */
void MechParser::readConstantInception(CamXML::Element &xml, Processes::ConstantInception &icn)
{
    std::string str;
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator j;

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "") icn.SetName(str);

     // Read initial particle composition.
    readInceptedComposition(xml, icn);

    // Read initial tracker variable values.
    readInceptedTrackers(xml, icn);

    // Read the fixed rate (particles per cm^3 per s)
    CamXML::Element *el = xml.GetFirstChild("rate");
    if (el != NULL) {
    	real rate = std::atof(el->Data().c_str());
    	// But store with units of particles per m^3 per s
    	icn.setConstantVolumetricInceptionRate(rate * 1.0e6);
    }
    else {
    	throw std::runtime_error("No rate provided for constant inception (Sweep::MechParser::readConstantInception)");
    }

    // Rate scaling now that a process has been created
    real A = 0.0;
    el = xml.GetFirstChild("A");
    if (el != NULL) {
        A = cdble(el->Data());
        icn.SetA(A);
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
            rxn.SetPropertyID(Sweep::iS);

            // Must scale rate constant from cm3 to m3, surface area
            // is multiplied by the site density so that it is a dimensionless
            // quantity hence A has units cm^3 s^-1.
            arr.A *= (1.0e-6);
        }
        else if (str.compare("s")==0) {
            // This reaction depends on surface area.  Ignore power,
            // they must have meant 1.
            rxn.SetPropertyID(Sweep::iS);
            arr.A *= (1.0e-2);
        }
        else if (str.compare("d")==0) {
            // This reaction depends on some power of the diameter.
            switch (power) {
                case 1:
                    rxn.SetPropertyID(Sweep::iDcol);
                    arr.A *= (1.0e-4);
                    break;
                case 2:
                    rxn.SetPropertyID(Sweep::iD2);
                    arr.A *= (1.0e-2);
                    break;
                case -1:
                    rxn.SetPropertyID(Sweep::iD_1);
                    arr.A *= (1.0e2);
                    break;
                case -2:
                    rxn.SetPropertyID(Sweep::iD_2);
                    arr.A *= (1.0e4);
                    break;
                default:
                    // Oh dear, can't have a zero power.
                    throw runtime_error("""particleterm"" tag found with "
                                        "invalid or zero power attribute "
                                        "(Sweep, MechParser::readSurfRxn)");
            }
        } else {
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
 * @param[in]       xml     XML document containing zero or more top level <coagulation> nodes
 * @param[in,out]   mech    Mechanism to which to add coagulation processes
 *
 * @exception   runtime_error   Kernel no specified
 * @exception   runtime_error   Unrecongised/unsupported kernel
 * @exception   runtime_error   Unrecongised/unsupported weight rule
 *
 * By default a transition regime coagulation process is created and added to the mechanism,
 * this will happen if no coagulation nodes are found immediately below the root node of the
 * XML tree.
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
    else
    {
        for(vector<CamXML::Element*>::const_iterator it = items.begin();
            it != items.end(); ++it) {

            // Read the user choice of kernel.  If no kernel is specified
            // a default will be used
            const CamXML::Element* const kernel = (*it)->GetFirstChild("kernel");
            if(kernel ==  NULL)
            {
                throw std::runtime_error("No kernel given for coagulation \
                                         (Sweep, MechParser::readCoagulation)");
            }
            else
            {
                // Work out which kernel to use and create it
                const string kernelName = kernel->Data();

                // See if there is a weight rule, which would indicate the use
                // of individual statistical weights for the particles
                CamXML::Element *weightXML = (*it)->GetFirstChild("weightrule");

                // Create a process of the appropriate type, but wrap it with an auto_ptr so
                // that it gets deleted if an exception is thrown when reading in the value
                // of A.
                std::auto_ptr<Processes::Coagulation> coag;

                if(weightXML == NULL) {
                    // Unweighted case
                    if(kernelName == "transition")
                        coag.reset(new Processes::TransitionCoagulation(mech));
                    else if(kernelName == "additive")
                        coag.reset(new Processes::AdditiveCoagulation(mech));
                    else if(kernelName == "constant")
                        coag.reset(new Processes::ConstantCoagulation(mech));
                    else
                        // Unrecognised option
                        throw std::runtime_error("Coagulation kernel " + kernelName + " not yet available in DSA \
                                                (Sweep, MechParser::readCoagulation)");

                    // Choice of position of newly coagulated particle
                    const CamXML::Element *positionChoiceXML = (*it)->GetFirstChild("positionchoice");

                    // This is an optional input
                    if (positionChoiceXML != NULL) {
                        const std::string choice = positionChoiceXML->Data();
                        if(choice == "none")
                            coag->SetPositionChoiceRule(Processes::Coagulation::NoPositionChoice);
                        else if (choice == "uniform")
                            coag->SetPositionChoiceRule(Processes::Coagulation::UniformPositionChoice);
                        else if (choice == "mass")
                            coag->SetPositionChoiceRule(Processes::Coagulation::MassPositionChoice);
                        else if (choice == "largestmass")
                            coag->SetPositionChoiceRule(Processes::Coagulation::LargestMassPositionChoice);
                        else
                            // Unrecognised option
                            throw std::runtime_error("Position choice rule " + choice + " not yet available \
                                                     (Sweep, MechParser::readCoagulation)");
                    }
                }
                else {
                    // weightXML != NULL so must have a weighted kernel
                    // need to find out what weight rule to use
                    const std::string weightRuleName = weightXML->Data();
                    Processes::CoagWeightRule weightRule;

                    if(weightRuleName == "w1" || weightRuleName == "harmonic") {
                        weightRule = Processes::CoagWeightHarmonic;
                        std::cerr << "Found harmonic weight rule" << std::endl;
                    }
                    else if(weightRuleName == "w2" || weightRuleName == "half")
                        weightRule = Processes::CoagWeightHalf;
                    else if(weightRuleName == "w3" || weightRuleName == "mass") {
                        weightRule = Processes::CoagWeightMass;
                        std::cerr << "Found mass weight rule" << std::endl;
                    }
                    else if(weightRuleName == "w4")
                        weightRule = Processes::CoagWeightRule4;
                    else
                        throw std::runtime_error("Coagulation weight rule " + weightRuleName + " not supported \
                                                 (Sweep, MechParser::readCoagulation)");

                    if(NULL != (*it)->GetFirstChild("positionchoice"))
                        throw std::runtime_error("Position choice rule does not apply with weighted coagulation \
                                                 (Sweep, MechParser::readCoagulation)");

                    // Now create the process
                    if(kernelName == "weightedadditive")
                         coag.reset(new Processes::WeightedAdditiveCoagulation(mech, weightRule));
                    else if(kernelName == "weightedconstant")
                         coag.reset(new Processes::WeightedConstantCoagulation(mech, weightRule));
                    else
                        // Unrecognised option
                        throw std::runtime_error("Coagulation kernel " + kernelName + " not yet available with weights \
                                                (Sweep, MechParser::readCoagulation)");
                }

                // Rate scaling now that a process has been created
                real A = 0.0;
                CamXML::Element *el = (*it)->GetFirstChild("A");
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
    }
}

// TRANSPORT PROCESSES

/*!
 *\param[in]        xml         XML document containing a mechanism description
 *\param[in,out]     mech        Mechanism to which transport processes will be added
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
                    tran.SetPropertyID(Sweep::iUniform);
                    break;
                case 1:
                    tran.SetPropertyID(Sweep::iDcol);
                    break;
                case 2:
                    tran.SetPropertyID(Sweep::iD2);
                    break;
                case -1:
                    tran.SetPropertyID(Sweep::iD_1);
                    break;
                case -2:
                    tran.SetPropertyID(Sweep::iD_2);
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
