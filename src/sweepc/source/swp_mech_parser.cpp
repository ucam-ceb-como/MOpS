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
#include "swp_phase.h"
#include "swp_inception.h"
#include "swp_constant_inception.h"
#include "swp_surface_reaction.h"
#include "swp_titania_surface_reaction.h"
#include "swp_titania_phase_transformation.h"
#include "swp_actsites_reaction.h"
#include "swp_condensation.h"
#include "swp_transcoag.h"
#include "swp_hybrid_transcoag.h"
#include "swp_addcoag.h"
#include "swp_constcoag.h"
#include "swp_hybrid_constcoag.h"
#include "swp_erosionfrag.h"
#include "swp_symmetricfrag.h"
#include "swp_weighted_addcoag.h"
#include "swp_weighted_constcoag.h"
#include "swp_weighted_transcoag.h"
#include "swp_transcoag_weighted_PAHs.h"
#include "swp_weighted_erosionfrag.h"
#include "swp_weighted_symmetricfrag.h"
#include "swp_coag_weight_rules.h"
#include "swp_frag_weight_rules.h"
#include "swp_silica_interparticle.h"
#include "swp_sprog_idealgas_wrapper.h"

#include "swp_tempReadColliPara.h" //temporarily used to read collision efficiency parameters, including mode, NONE, MAX, MIN, COMBINED

#include "gpc_species.h"

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
 *@param[in]        xml     XML node with component children
 *@param[in,out]    icn     Inception instance for which to set composition of newly incepted particles
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
            double dx = cdble(str);
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
 *@param[in]        xml     XML node with track children
 *@param[in,out]    icn     Inception instance for which to set tracker values of newly incepted particles
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
            double dx = cdble(str);
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

	// Assign phases
	assignPhases(mech);

    // READ DEFINED MODELS.

    xml.Root()->GetChildren("model", items);

    for (i=items.begin(); i!=items.end(); ++i) {
        // Check the model type.
        string str = (*i)->GetAttributeValue("type");

        if (str == "particle") {
            throw std::runtime_error("<model id=\"particle\" is not a valid tag (Sweep::MechParser::readV1");
        } else if (str == "drag") {
            // Read drag model ID.
            str = (*i)->GetAttributeValue("id");

            if(str == "StokesCunningham") {
                // Parameters of the drag expression
                double A, B, E;

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
                double coeff;
                if(numberXML != NULL) {
                    coeff = atof(numberXML->Data().c_str());
                }
                else {
                    throw std::runtime_error("Parameter A value must be given for drag proportional to temperature (Sweep, MechParser::readV1)");
                }

                // This option is just for testing so just reuse the support
                // for the Knudsen regime
                mech.SetKnudsenDragConstants(coeff, 0.0, 0.0);


            } else if (str == "LiWang") {
                mech.SetDragType(Sweep::ParticleModel::LiWangDrag);
            } else if (str == "LiWangPat") {
                mech.SetDragType(Sweep::ParticleModel::LiWangPatDrag);
            } else {
                throw std::runtime_error("Unrecognised drag model id (Sweep, MechParser::readV1).");
            }

        } else if (str == "diffusion") {
            // Read diffusion model ID.
            str = (*i)->GetAttributeValue("id");

            mech.setAvgMolWtIndex(SprogIdealGasWrapper::sAvgMolWtIndex);

            if((str == "Flamelet") || (str == "flamelet")) {
                // Diffusion according to soot flamelet equation in flamelet space
                mech.setDiffusionType(Sweep::ParticleModel::FlameletDiffusion);
                mech.setMixFracGradIndex(SprogIdealGasWrapper::sMixFracGradientIndex);
            }
            else if((str == "Einstein") || (str == "einstein")) {
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
            if((str == "Flamelet") || (str == "flamelet")) {
                // Advection according to soot flamelet equation in flamelet space
                mech.setAdvectionType(Sweep::ParticleModel::FlameletAdvection);

                // Indices to enable the model to find the quantities needed to calculate
                // the advection in the flamelet model
                mech.setMixFracDiffIndex(SprogIdealGasWrapper::sMixFracDiffusionIndex);
                mech.setMixFracGradIndex(SprogIdealGasWrapper::sMixFracGradientIndex);
                mech.setMixFracLaplIndex(SprogIdealGasWrapper::sMixFracLaplacianIndex);
            }
            else if((str == "Physical") || (str == "physical")) {
                // Advection at bulk gas velocity in physical space
                mech.setAdvectionType(Sweep::ParticleModel::BulkAdvection);
            }
            else {
                throw std::runtime_error("Unrecognised advection model id (Sweep, MechParser::readV1).");
            }
        }else if (str == "thermophoresis") {
            // Read advection model ID.
            str = (*i)->GetAttributeValue("id");
            if((str == "Waldmann") || (str == "waldmann")) {
                // Thermophoretic velocity that is the same for all particles
                mech.setThermophoresisType(Sweep::ParticleModel::WaldmannThermophoresis);
                mech.setThermalConductIndex(SprogIdealGasWrapper::sThermalConductivityIndex);
                mech.setTGradIndex(SprogIdealGasWrapper::sTemperatureGradientIndex);
            }
            else if(str == "LiWang") {
                // Particle dependence according to model of Li and Wang
                mech.setThermophoresisType(Sweep::ParticleModel::LiWangThermophoresis);
                mech.setThermalConductIndex(SprogIdealGasWrapper::sThermalConductivityIndex);
                mech.setTGradIndex(SprogIdealGasWrapper::sTemperatureGradientIndex);
            }
            else if((str == "None") || (str == "none")) {
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
    } else if (str == "surfvolhydrogen") {
        mech.SetAggModel(AggModels::SurfVolHydrogen_ID);
    } else if (str == "surfvolsilica") {
        mech.SetAggModel(AggModels::SurfVolSilica_ID);
    } else if (str == "surfvolcubic") {
        mech.SetAggModel(AggModels::SurfVolCubic_ID);
    } else if (str == "PAH") {
    // Reject all old style input files
        throw std::runtime_error("PAH-PP MODEL are no longer supported (Sweep::MechParser::readV1), you can use NEW PAH_KMC model");
    } else if (str == "PAH_KMC") {
        mech.SetAggModel(AggModels::PAH_KMC_ID);
    } else if (str == "silica") {
        throw std::runtime_error("Old silica model is deprecated. Use bintreesilica or surfvolsilica."
                " in Sweep::MechParser::readV1");
    } else if (str == "bintree") {
        mech.SetAggModel(AggModels::BinTree_ID);
    } else if (str == "bintreesilica") {
        mech.SetAggModel(AggModels::BinTreeSilica_ID);

    } else {
        mech.SetAggModel(AggModels::Spherical_ID);
    }

    string strh = particleXML->GetAttributeValue("hybrid");
    if (strh == "true")
    {
        mech.SetHybrid(true);
        string strn = particleXML->GetAttributeValue("threshold");
        unsigned int threshold = 100;
        if (strn != "")
            threshold = (unsigned int)(cdble(strn));
        if (threshold < 1)
            throw std::runtime_error("Hybrid threshold must be positive. (Sweep::MechParser::readV1)");
		mech.SetHybridThreshold(threshold);
		string strc = particleXML->GetAttributeValue("coagulate-in-list");
		if (strc == "true")
			mech.SetCoagulateInList(true);
		else
			mech.SetCoagulateInList(false);
    }
    else
    {
        mech.SetHybrid(false);
		mech.SetHybridThreshold(0);
		mech.SetCoagulateInList(false);
    }

    //! Check whether to track the distance betweeen the centres of primary
    //! particles or the coordinates of the primary particles, but this only
    //! applies to the binary tree and PAH-KMC (and maybe surface-volume)
    //! models.
    //! By default, in terms of type space, primary particles are assumed to be
    //! in point contact but the partial rounding/sintering of particles is
    //! implicitly taken into account where mass addition (surface growth and
    //! condensation) and sintering events modify the common surface area
    //! between pairs of neighbouring primaries. By turning on this flag, mass
    //! addition events lead to an increase in primary diameter at a rate
    //! dependent on its neighbour while sintering leads to a simultaneous
    //! decrease in the distance betweeen the centres of neighbouring primaries
    //! and an increase in their diameters. Note that if the coordinates of the
    //! primary particles are tracked their separation is known.
    if (mech.AggModel() == AggModels::PAH_KMC_ID || mech.AggModel() == AggModels::BinTree_ID) {
        str = particleXML->GetAttributeValue("track");
        if (str == "primarySeparation") {
            mech.setTrackPrimarySeparation(true);
        } else if (str == "primaryCoordinates") {
            mech.setTrackPrimaryCoordinates(true);
        } else {
            mech.setTrackPrimarySeparation(false);
            mech.setTrackPrimaryCoordinates(false);
        }
    } else {
        mech.setTrackPrimarySeparation(false);
        mech.setTrackPrimaryCoordinates(false);
    }

    // See if there are any secondary particle criteria
    const CamXML::Element* secondaryXML = particleXML->GetFirstChild("secondaryparticle");
    if(secondaryXML != NULL) {
        throw std::runtime_error("Secondary particles are no longer supported (Sweep::MechParser::readV1)");
    }

    // Get the coalescence threshold for a multicomponent binary tree model
    const CamXML::Element* el = particleXML->GetFirstChild("coalthresh");
    if (mech.AggModel() == AggModels::BinTree_ID) {
        if (el != NULL) {
            double ct = cdble(el->Data());
            if (ct < 0.0 || ct > 2.0) {
                throw std::runtime_error("Coalescence threshold must be 0<ct<2.0. (Sweep::MechParser::readV1)");
            } else {
            mech.SetBinTreeCoalThresh(ct);
            }
        } else {
            throw std::runtime_error("Must specify coalescence threshold in <particle>"
                    "block with tag <coalthresh> for bintree particle model. (Sweep::MechParser::readV1)");
        }
    }

    el = particleXML->GetFirstChild("fractdim");
    // Get the fractal dimension
    if (el != NULL) {
        double df = cdble(el->Data());
        if (df < 1.0 || df > 3.0) {
            throw std::runtime_error("Fract. dim must be 1<df<3.0. (Sweep::MechParser::readV1)");
        } else {
            mech.SetFractDim(df);
        }
    } else {
        // Default fractal dimension is 1.8
        // Schaefer & Hurd (1990) Aer. Sci. Tech. 12:876-890
        // for synthesis of fumed silica nanoparticles.
        mech.SetFractDim(1.8);
    }

    el = particleXML->GetFirstChild("efm");
    // Get the free molecular enhancement factor
    if (el != NULL) {
        double efm = cdble(el->Data());
        if (efm < 0.0) {
            throw std::runtime_error("EFM must be >0. (Sweep::MechParser::readV1)");
        } else {
            mech.SetEnhancementFM(efm);
        }
    } else {
        // Default enhancement factor is 2.2
        // Kazakov & Frenklach, Combustion & Flame, 1998 113:484-501
        // for coagulation of spherical non-polar soot particles
        mech.SetEnhancementFM(2.2);
    }

    // Get the sintering model.
    particleXML->GetChildren("sintering", items);
    if (items.size() > 0) {
        i = items.begin();

        // Check if sintering model enabled.
        str = (*i)->GetAttributeValue("enable");
        if (str == "true") {
            mech.SintModel().Enable();

        // Get sintering model type.
        str = (*i)->GetAttributeValue("model");
        if (str == "viscous_flow") {
            mech.SintModel().SetType(SinteringModel::ViscousFlow);
        } else if (str == "rutile") {
            // Special MD fit for GBD sintering of rutile
            mech.SintModel().SetType(SinteringModel::Rutile);
        } else if (str == "ssd") {
            // Solid-state diffusion (d^3)
            mech.SintModel().SetType(SinteringModel::SSD);
        } else if (str == "silicon") {
            // Solid-state diffusion (d^3)
            mech.SintModel().SetType(SinteringModel::Silicon);
        } else if (str == "kirchoff") {
            // Constant characterisitic sintering time
            mech.SintModel().SetType(SinteringModel::SilicaKirchoff);
        } else if (str == "constant") {
            // Constant characterisitic sintering time
            mech.SintModel().SetType(SinteringModel::Constant);
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
            mech.SintModel().Disable();
        }

		// Get critical exponent alpha (default is alpha = 1)
		// TODO: add paper reference
		if ((*i)->GetFirstChild("alpha") != NULL) {
			str = (*i)->GetFirstChild("alpha")->Data();
			mech.SintModel().Setalpha(cdble(str));
		}
		else{
			mech.SintModel().Setalpha(1.0);
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
        readInterParticles(xml, mech);
        readCondensations(xml, mech);
		readPhaseTransformation(xml, mech);
    }
    readDiffusionProcs(xml, mech);
    readAdvectionProcs(xml, mech);

    // Read which coagulation kernel to use
    readCoagulation(xml, mech);

    // Read which coagulation kernel to use
    readFragmentation(xml, mech);
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
			
			// Applies to flamepp:
			// If any component is soot and the postprocessingtype is not wdotA4 
			// then set postprocessingtype to XA4 (the default type is concentration).
			// Note: this assumes that all soot simulations have a component id "soot",
			// otherwise the generic mechanism is used. 
			if (str == "soot" && mech.Postprocessing() != Sweep::ParticleModel::wdotA4){
				mech.setPostprocessingType(Sweep::ParticleModel::XA4);
			}

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

        readColliPara(i,mech);

        // Get InceptedPAH, only pyrene and benzene are supported.
        el = (*i)->GetFirstChild("InceptedPAH");
        if (el!=NULL) {
            str = el->Data();
            if (str != "" && str != "Multiple" && str != "FromFile") {
                mech.SetInceptedPAH(str);
            } else if (str == "Multiple"){
				cout << "Incepting multiple PAHs.\n";
				vector<std::string> incepted_pah_list_str;
				vector<CamXML::Element*> inception_items;
				vector<CamXML::Element*>::iterator it_inc;
				CamXML::Element *my_el;
				my_el = (*i)->GetFirstChild("InceptedPAH_items");
				if (my_el != NULL){
					//Read the list
					my_el->GetChildren("PAH", inception_items);
				}else {
					throw std::runtime_error("No inception PAH list found. (InceptedPAH_items)");
				}
				if(inception_items.size() == 0) {
					std::string msg("Component ");
					msg += comp->Name();
					msg += " InceptedPAH_items contains no data (Sweep, MechParser::readComponents).";
					delete comp;
					throw runtime_error(msg);
				}
				for (it_inc=inception_items.begin(); it_inc!=inception_items.end(); ++it_inc){
					std::string inc_pah_str = (*it_inc)->Data();
					cout << "Incepting PAH: " << inc_pah_str << "\n";
					mech.SetInceptedPAH(inc_pah_str);
				}
			}
			else if (str == "FromFile"){
				cout << "Incepting PAH from file. Reading file InceptedPAH.inx.\n";
				std::ifstream src("InceptedPAH.inx");
				if (src.good()){
					mech.SetInceptedPAH("FromFile");
				}
				else {
					cout << "Unable to open file InceptedPAH.inx. Setting incepted PAH to A4.\n";
					mech.SetInceptedPAH("A4");
				}
				src.close();
			}else {
                // coalthresh contains no data.
                std::string msg("Component ");
                msg += comp->Name();
                msg += " InceptedPAH contains no data (Sweep, MechParser::readComponents).";

                delete comp;
                throw runtime_error(msg);
            }
        } else {
            mech.SetInceptedPAH("A4");
        }
		
		// Get flag for test rates.
        el = (*i)->GetFirstChild("DebugPAH");
        if (el!=NULL) {
            str = el->Data();
			if (str != "") {
				if (std::stoi(str)==1){
					cout << "DebugPAH set to True. Saving test rates for each jump process.";
					comp->SetDebugPAH(true);
				}else{
					cout << "DebugPAH defined but not set to True. Not saving test rates for each jump process.";
					comp->SetDebugPAH(false);
				}
			}
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

        //! Model parameter. 
        /*!
         * The number of 6-member rings in a PAH in a particle below which it is removed.
         */
        el = (*i)->GetFirstChild("thresholdOxidation");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetThresholdOxidation(int(cdble(str)));
            } else {
                std::string msg("Component ");
                msg += comp->Name();
                msg += " thresholdOxidation contains no data (Sweep, MechParser::readComponents).";

                delete comp;
                throw runtime_error(msg);
            }
        } else {
            comp->SetThresholdOxidation(0);
        }

        //! Numerical parameter
        /*!
         * Allow PAHs in soot particles to point to the same memory location after a doubling event.
         */
        el = (*i)->GetFirstChild("sharedPointers");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetSharedPointers(int(cdble(str)));
            } else {
                std::string msg("Component ");
                msg += comp->Name();
                msg += " sharedPointers contains no data (Sweep, MechParser::readComponents).";

                delete comp;
                throw runtime_error(msg);
            }
        } else {
            comp->SetSharedPointers(0);
        }
		
		//! Numerical parameter
		/*!
		* Allow PAHs in soot particles to point to the same memory location after a doubling event.
		*/
		el = (*i)->GetFirstChild("weightedPAHs");
		if (el != NULL) {
			str = el->Data();
			if (str != "") {
				comp->SetWeightedPAHs(int(cdble(str)));
			}
			else {
				std::string msg("Component ");
				msg += comp->Name();
				msg += " weightedPAHs contains no data (Sweep, MechParser::readComponents).";

				delete comp;
				throw runtime_error(msg);
			}
		}
		else {
			comp->SetWeightedPAHs(0);
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

		//Phase

		//get the phase the component belongs to
		el = (*i)->GetFirstChild("phase");
		if (el != NULL) {
			str = el->Data();
			if (str != "") {
				comp->SetPhase(str);
			}
		}

		//get the element the component belongs to
		el = (*i)->GetFirstChild("element");
		if (el != NULL) {
			str = el->Data();
			if (str != "") {
				comp->SetElement(str);
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

// Assign components to phases
void MechParser::assignPhases(Sweep::Mechanism &mech)
{

	// Loop over components and assign to phases
	for (unsigned int i = 0; i < mech.ComponentCount(); ++i) {

		//Get phase of component
		std::string phaseName = mech.Components(i)->Phase();

		//Does phase already exist?
		int i_phase = mech.PhaseIndex(phaseName);

		if (i_phase >= 0) //phase exists
		{
			Phase *phase = mech.Phases(i_phase);
			
			// Add component to existing phase
			phase->AddComponent(i);
		}
		else{
			// Phase does not exist: create a new phase
			Phase *phase = new Phase();

			// Phase name/symbol
			phase->SetName(phaseName);

			// Add component to phase
			phase->AddComponent(i);

			// Add phase to mechanism.
			mech.AddPhase(*phase);
		}
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

        string inceptype = (*(*i)).GetAttributeValue("type");

        if (inceptype == "silicon") {
            // Create new silicon inception.
            SiliconInception *icn = new SiliconInception(mech);
            icn->SetMechanism(mech);
            icn->SetName("Inception " + cstr(k));

            try {
                readSiliconInception(*(*i), *icn);
            }
            catch (std::exception &e) {
                delete icn;
                throw;
            }

            // Set silicon-specific constants
            icn->GenerateSpeciesData(mech);

            // Add inception to mechanism.  Once entered into mechanism, the mechanism
            // takes control of the inception object for memory management.
            mech.AddInception(*icn);

        } else {

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
    double A = 0.0;
    const CamXML::Element *el = xml.GetFirstChild("A");
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


/*!
 * @brief           Reads a silicon inception process from an XML element.
 *
 * @param xml       XML node to read in
 * @param icn       Pointer to new inception process
 */
void MechParser::readSiliconInception(CamXML::Element &xml, Processes::SiliconInception &icn)
{
    string str;
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator j;

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "") icn.SetName(str);

    // Read inception mechanism to use
    str = xml.GetAttributeValue("mech");
    if (str == "vbdz") {
        icn.SetInceptionMechanism(Sweep::Processes::SiliconInception::iVBDZ);
    } else if (str == "girshick") {
        icn.SetInceptionMechanism(Sweep::Processes::SiliconInception::iGirshick);
    } else if (str == "collisional") {
        icn.SetInceptionMechanism(Sweep::Processes::SiliconInception::iCollisional);
    } else {
        throw std::runtime_error("Unrecognised mechanism " +
                str + " in Sweep::MechParser::readSiliconInception");
    }

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
            throw std::runtime_error("Unrecognised rate type " + str + " in Sweep::MechParser::readSiliconInception");
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
        throw std::runtime_error("One or two inception species must be specified in Sweep::MechParser::readSiliconInception");


    // Rate scaling now that a process has been created
    double A = 0.0;
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
        PAHInception *icn = new PAHInception(mech, SprogIdealGasWrapper::sPAHFormationIndex);
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
 * Read constant inception processes from a sweep mechanism XML file.
 *
 *@param[in]        xml     XML node with constantinception children
 *@param[in]        mech    Mechanism of which the inceptions are to be part
 *
 *@exception    std::runtime_error  Missing parameters for the lognormal distributions
 */
void MechParser::readConstantInceptions(CamXML::Document &xml, Sweep::Mechanism &mech)
{

    // Get list of inceptions from XML data.
    vector<CamXML::Element*> inceptions;
    xml.Root()->GetChildren("constantinception", inceptions);

    vector<CamXML::Element*>::iterator i = inceptions.begin();
    for (unsigned k=0; i!=inceptions.end(); ++i,++k) {
        // Read the data before creating the object
        vector<CamXML::Element*> components;
        (*i)->GetChildren("component", components);

        // Parameters specifying the log normal distributions for each component
        std::vector<double> locations, scales;
        // Component ids in the mechanism
        std::vector<int> ids;

        for (vector<CamXML::Element*>::iterator j=components.begin(); j!=components.end(); ++j) {
            // Get component ID.
            std::string str = (*j)->GetAttributeValue("id");
            ids.push_back(mech.ComponentIndex(str));

            if (ids.back() >= 0) {
                std::string str = (*j)->GetAttributeValue("dx");
                if(str != "") {
                    // Deterministic size increment specified
                    double val = atof(str.c_str());
                    if(val <= 0)
                        throw std::runtime_error("Amount of component in a newly incepted particle must be strictly positive, not "
                                                 + str + ". MechParser::readConstantInceptions");
                    locations.push_back(log(cdble(str)));
                    // Set the standard deviation parameter to 0
                    scales.push_back(0.0);

                }
                else {
                    str = (*j)->GetAttributeValue("distribution");
                    if(str == "lognormal") {
                        // Lognormal component distribution
                        // Mean of the logarithm of the amount of this component
                        str = (*j)->GetAttributeValue("location");
                        double val = atof(str.c_str());
                        if(val <= 0)
                            throw std::runtime_error("Location of lognormal size distribution must be positive, not "
                                                     + str + ". MechParser::readConstantInceptions");
                        locations.push_back(cdble(str));

                        // This is the std dev of the amount of this component
                        str = (*j)->GetAttributeValue("scale");
                        val =  atof(str.c_str());
                        if(val < 0)
                            throw std::runtime_error("Scale of lognormal size distribution must be non-negative, not "
                                                     + str + ". MechParser::readConstantInceptions");
                        scales.push_back(cdble(str));
                    }
                    else
                        throw std::runtime_error("Unrecognised component distribution " + str
                                                 + " in MechParser::readConstantInceptions");
                }
            }
            else {
                // Unknown component in mechanism.
                throw std::runtime_error(str + ": Component not found in mechanism (readConstantInceptions)");
            }
        }

        if(ids.size() != mech.ComponentCount())
            throw std::runtime_error("Missing component in (readConstantInceptions)");

        // Sort the component data by component id
        std::vector<double> locationsSorted(mech.ComponentCount()), scalesSorted(mech.ComponentCount());
        std::cerr << mech.ComponentCount() << ',' << ids.size() << ',' << locations.size() << ',' << scales.size() << std::endl;
        for(unsigned l = 0; l < mech.ComponentCount(); ++l) {
            std::cerr << l << ',' << ids[l] << std::endl;
            locationsSorted[ids[l]] = locations[l];
            scalesSorted[ids[l]]    = scales[l];
        }

        // Rate of inception events
        double A = 0.0;
        CamXML::Element* el = (*i)->GetFirstChild("A");
        if (el != NULL) {
            A = cdble(el->Data());
            // But store with units of particles per m^3 per s
            A *= 1.0e6;
        }

        std::auto_ptr<ConstantInception> icn(new ConstantInception(mech, A, locationsSorted, scalesSorted));

        // See if the inception is at a fixed position
        el = (*i)->GetFirstChild("fixedposition");
        if (el != NULL) {
            // Read position in cm (for consistency with rest of input file ...
            double posn = std::atof(el->Data().c_str());
            // ..., but store position in m
            icn->setFixedPosition(posn/100);
            icn->useFixedPosition(true);
        }

        std::string str = (*i)->GetAttributeValue("name");
        if (str != "")
            icn->SetName(str);
        else
            icn->SetName("ConstantInception" + cstr(k));

        // Read initial tracker variable values.
        readInceptedTrackers(**i, *icn);

        // Add inception to mechanism.  Once entered into mechanism, the mechanism
        // takes control of the inception object for memory management.
        mech.AddInception(*icn);
        icn.release();
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
            rxn = new ActSiteReaction(mech, ActSiteReaction::ABFRadicalSiteModel, SprogIdealGasWrapper::sAlphaIndex);
        } else if (str.compare("bp")==0) {
            // This is a Blanquart Pitsch active-sites enabled reaction.
            rxn = new ActSiteReaction(mech, ActSiteReaction::BPRadicalSiteModel, SprogIdealGasWrapper::sAlphaIndex);
        } else if (str.compare("abfconst")==0) {
            // This is an ABF active-sites reaction with alpha constant (=1).
            rxn = new ActSiteReaction(mech, ActSiteReaction::ABFConstant, SprogIdealGasWrapper::sAlphaIndex);
        } else if (str.compare("titania")==0) {
            // This is a titania surface reaction
            // Create a new reaction of a certain form
            string str2 = (*i)->GetAttributeValue("form");
            if (str2.compare("firstorder")==0) {
                rxn = new TitaniaSurfaceReaction(mech, TitaniaSurfaceReaction::iFirstOrder);
            } else if (str2.compare("ghoshtagore")==0) {
                rxn = new TitaniaSurfaceReaction(mech, TitaniaSurfaceReaction::iGhoshtagore);
            } else if (str2.compare("adsorption")==0) {
                rxn = new TitaniaSurfaceReaction(mech, TitaniaSurfaceReaction::iEleyRidealAdsorption);
            } else if (str2.compare("desorption")==0) {
                rxn = new TitaniaSurfaceReaction(mech, TitaniaSurfaceReaction::iEleyRidealDesorption);
            } else if (str2.compare("multivariate")==0) {
                rxn = new TitaniaSurfaceReaction(mech, TitaniaSurfaceReaction::iMultivariate);
            } else {
                throw runtime_error("Unrecognised titania reaction form" + str2 +
                        "in MechParser::readSurfRxns");
            }
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
 *    A \left[Y\right] T^n e^\frac{E}{RT} \times \mathrm{particle\ property},
 * \f]
 * where \f$\left[Y\right]\f$ gives the number of moles of a gas
 * phase species \f$Y\f$ per cubic centimetre and the \f$T^n\f$ term is
 * treated as dimensionless.  This means that for a dimensionless particle
 * property, such as the number of active sites the dimensions of A must
 * be \f$\mathrm{cm}^3 \mathrm{mol}^{-1} \mathrm{s}^{-1}\f$.  If the
 * particle property has dimensions, for example collision diameter
 * raised to the power \f$k\f$, then \f$A\f$ has different units:
 * \f$\mathrm{cm}^{3-k} \mathrm{mol}^{-1} \mathrm{s}^{-1}\f$.  If the
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

    // Read reactants
    readReactants(xml, rxn);

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
        } else if (str.compare("u")==0) {
            // Reaction is not dependent on any particle property and is a
            // constant rate.
            rxn.SetPropertyID(Sweep::iUniform);
        } else if (str.compare("s")==0) {
            // This reaction depends on surface area.  Ignore power,
            // they must have meant 1.
            rxn.SetPropertyID(Sweep::iS);
            arr.A *= (1.0e-2);
        }
        else if (str.compare("coverage")==0) {
            // This reaction depends on fractional coverage (dimensionless),
            // the constant does not require scaling as its units are s^-1.
            rxn.SetPropertyID(Sweep::iCoverage);
        }
        else if (str.compare("asn")==0) {
            rxn.SetPropertyID(Sweep::iASN);

            // Must scale rate constant from cm3 to m3, surface area
            // is multiplied by the site density so that it is a dimensionless
            // quantity hence A has units cm^3 s^-1.
            arr.A *= (1.0e-6);
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
                    arr.A *= (1.0e-8);
                    break;
                case -2:
                    rxn.SetPropertyID(Sweep::iD_2);
                    arr.A *= (1.0e-10);
                    break;
                default:
                    // Oh dear, can't have a zero power.
                    throw runtime_error("""particleterm"" tag found with "
                                        "invalid or zero power attribute "
                                        "(Sweep, MechParser::readSurfRxn)");
            }
        } else {
            throw std::runtime_error("Surface process defined with unrecognised ""particleterm"" "
                                        "(Sweep, MechParser::readSurfRxn).");
        }
    } else {
        throw runtime_error("Surface process defined without ""particleterm"" "
                            "element (Sweep, MechParser::readSurfRxn).");
    }

    rxn.SetArrhenius(arr);
    // Ensure the Process object has the same pre-exponential in case it
    // is mistakenly used instead of the Arrhenius object.
    //rxn.SetA(arr.A);

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
    double A = 0.0;
    CamXML::Element *el = xml.GetFirstChild("A");
    if (el != NULL) {
        A = cdble(el->Data());
    } else {
        A = 1.0;
    }
    cond.SetA(A);
}


// INTER PARTICLE REACTIONS.

/*!
 *  Read interparticle reaction processes from a sweep mechanism XML file.
 *
 *@exception    std::runtime_error      Failed to find \f$\mathrm{H}_4\mathrm{O}_4\mathrm{Si}\f$ in mechanism
 */
void MechParser::readInterParticles(CamXML::Document &xml, Mechanism &mech)
{
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    string str;
    unsigned int k = 0;

    // Get all surface reactions.
    xml.Root()->GetChildren("interparticle", items);

    for (i=items.begin(),k=0; i!=items.end(); ++i,++k) {
        // Look up the index for H4O4SI
        const int isp = Sprog::Species::Find("H4O4SI", *mech.Species());
        if(isp < 0)
            throw std::runtime_error("Could not find species H4O4SI in MechParser::ReadInterParticles");

        // Create a new interparticle object.
        InterParticle *intpar = new InterParticle(mech, isp);

        // Set default name.
        intpar->SetName("Inter Particle " + cstr(k));

        // Read the inter particle properties.
        try {
            readInterParticle(*(*i), *intpar);
        } catch (std::exception &e) {
            delete intpar;
            throw;
        }

        // Add Inter particle to mechanism.
        mech.AddProcess(*intpar);
    }
}

// Reads a interparticle process from a sweep mechanism XML file.
void MechParser::readInterParticle(CamXML::Element &xml, Processes::InterParticle &intpar)
{
    string str;
    CamXML::Element *el = NULL;

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "") intpar.SetName(str);

    // Is reaction deferred.
    str = xml.GetAttributeValue("defer");
    if (str=="true") {
        intpar.SetDeferred(true);
    } else {
        intpar.SetDeferred(false);
    }

    readReactants(xml, intpar);

    // Read products.
    readProducts(xml, intpar);

    // Read particle composition change.
    readCompChanges(xml, intpar);

    // Read tracker variable changes.
    readTrackChanges(xml, intpar);

    //========== Read Arrhenius rate parameters ======================
    Sprog::Kinetics::ARRHENIUS arr;
    el = xml.GetFirstChild("A");
    if (el != NULL) {
        arr.A = cdble(el->Data());
    } else {
        // Reaction must have constant.
        throw runtime_error("InterParticle reaction found with no rate constant "
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
            intpar.SetPropertyID(Sweep::iS);

            // Must scale rate constant from cm3 to m3, surface area
            // is multiplied by the site density so that it is a dimensionless
            // quantity hence A has units cm^3 s^-1.
            arr.A *= (1.0e-6);
        }
        else if (str.compare("s")==0) {
            // This reaction depends on surface area.  Ignore power,
            // they must have meant 1.
            intpar.SetPropertyID(Sweep::iS);
            arr.A *= (1.0e-2);
        }
        else if (str.compare("asn")==0) {
            intpar.SetPropertyID(Sweep::iASN);

            // Must scale rate constant from cm3 to m3, surface area
            // is multiplied by the site density so that it is a dimensionless
            // quantity hence A has units cm^3 s^-1.
            //arr.A *= (1.0e-6);
            arr.A *= (1.0e-6);
        }
        else if (str.compare("d")==0) {
            // This reaction depends on some power of the diameter.
            switch (power) {
                case 1:
                    intpar.SetPropertyID(Sweep::iDcol);
                    arr.A *= (1.0e-4);
                    break;
                case 2:
                    intpar.SetPropertyID(Sweep::iD2);
                    arr.A *= (1.0e-2);
                    break;
                case -1:
                    intpar.SetPropertyID(Sweep::iD_1);
                    arr.A *= (1.0e2);
                    break;
                case -2:
                    intpar.SetPropertyID(Sweep::iD_2);
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


    intpar.SetArrhenius(arr);
    // Ensure the Process object has the same pre-exponential in case it
    // is mistakenly used instead of the Arrhenius object.
    intpar.SetA(arr.A);

}

// PHASE TRANSFORMATION REACTION.

/*!
 *  Read phase transformation process from a sweep mechanism XML file.
 *
 *@exception    std::runtime_error   
 */
void MechParser::readPhaseTransformation(CamXML::Document &xml, Mechanism &mech)
{
	vector<CamXML::Element*> items, subitems, subsubitems;
	vector<CamXML::Element*>::iterator i, j, l, m;
	CamXML::Element *el = NULL;
	std::string str;
	std::string type;
	unsigned int k = 0;

	// Get all phase transformations
	xml.Root()->GetChildren("transformation", items);

	for (i = items.begin(), k = 0; i != items.end(); ++i, ++k) {

		// transformation type
		type = (*i)->GetAttributeValue("type");

		//----------------------------------------------------
		//--------------- titania_melting_model --------------
		//----------------------------------------------------
		if (type.compare("melting") == 0 || type.compare("gibbs") == 0 || type.compare("composition") == 0) {

			// Check that the particle model is valid
			if (mech.AggModel() == AggModels::Spherical_ID || mech.AggModel() == AggModels::BinTree_ID) {
			
				// Enable model
				mech.MeltModel().Enable();

				// Set model type
				if (type.compare("melting") == 0){
					mech.MeltModel().SetType(MeltingModel::Melting);
				}
				else if (type.compare("gibbs") == 0) {
					mech.MeltModel().SetType(MeltingModel::Gibbs);
				}
				else if (type.compare("composition") == 0){
					mech.MeltModel().SetType(MeltingModel::Composition);
				}

				// Get the name of liquid phase
				el = (*i)->GetFirstChild("liquid");
				if (el != NULL) {
					str = el->GetAttributeValue("id");
				}
				else{ // Otherwise assume the liquid phase is named "Liquid"
					str = "Liquid";
				}
				// Set liquid phase
				int i_phase = mech.PhaseIndex(str);
				if (i_phase >= 0){// Check that a liquid phase exists
					Phase *phase = mech.Phases(i_phase);
					phase->SetLiquid();
				}
				else{
					throw runtime_error(str+" Liquid phase not found in mechanism "
						"Sweep, MechParser::readPhaseTransformation).");
				}			

				// Loop over phase transformations
				(*i)->GetChildren("phase", subitems);
				for (j = subitems.begin(); j != subitems.end(); ++j) {

					double A = 0.0;
					double T = 0.0;

					// crystal phase name or symbol
					string phasename = (*j)->GetAttributeValue("id");

					// Get component ID.
					int iid = mech.PhaseIndex(phasename);
					if (iid < 0) {
						throw runtime_error(str + ": Phase not found in mechanism "
							"Sweep, MechParser::readPhaseTransformation).");
					}

					//Temperature parameters
					//parameter
					el = (*j)->GetFirstChild("A");
					if (el != NULL) A = cdble(el->Data());
					//bulk melting temperature
					el = (*j)->GetFirstChild("Tbulk");
					if (el != NULL) {
						T = cdble(el->Data());
					}
					else{
						throw runtime_error(phasename +
							"Temperature not specified (Sweep, MechParser::readPhaseTransformation).");
					}

					vector<fvector> ddcomp;		//vector of component changes
					vector<std::string> elems;	//vector of elements

					//get components and sort/separate by element
					(*j)->GetChildren("component", subsubitems);
					for (m = subsubitems.begin(); m != subsubitems.end(); ++m) {
						
						fvector dcomp(mech.ComponentCount(), 0.0);
						
						// Get component ID.
						string str = (*m)->GetAttributeValue("id");
						int id = mech.ComponentIndex(str);
						if (id >= 0) {
							// Get component element
							std::string elem = mech.Components()[id]->Element();
							// Get component change.
							str = (*m)->GetAttributeValue("dx");
							double dx = cdble(str);
							// Check if element already exists in element vector
							int i_elem = -1;							
							for (unsigned int ii = 0; ii < elems.size(); ++ii){
								if ((elems[ii]).compare(elem) == 0) i_elem = ii;
							}
							if (i_elem >= 0){ //element exists
								// Set component change.
								ddcomp[i_elem][id] += dx;
							}
							else{ // Add new element to vector
								elems.push_back(elem);
								// Add composition change
								dcomp[id] += dx;
								ddcomp.push_back(dcomp);
							}
						}
						else {
							throw runtime_error(str + ": Component not found in mechanism "
								"Sweep, MechParser::readPhaseTransformation).");
						}
					}

					// Check mass conservation of component change
					vector<fvector>::const_iterator jj;
					for (jj = ddcomp.begin(); jj != ddcomp.end(); ++jj){

						double mass = 0.0;
						double avg_mass = 0.0;
						for (unsigned int ii = 0; ii != mech.ComponentCount(); ii++){
							mass += mech.Components(ii)->MolWt() * (*jj)[ii] / NA;
							avg_mass += mech.Components(ii)->MolWt() / NA;
						}
						avg_mass /= mech.ComponentCount();

						if (mass > 0.001*avg_mass){ // tolerance of 0.1% average component mass
							throw runtime_error(phasename + "Phase transformation does not conserve mass "
								"Sweep, MechParser::readPhaseTransformation).");
						}

					}

					// Add phase change to melting model
					mech.MeltModel().AddPhase(phasename, A, T, ddcomp, iid);
				}

				// Get crossover diameter (Gibbs and Composition models). Only one crossover is supported
				if (type.compare("gibbs") == 0 || type.compare("composition") == 0) {
					(*i)->GetChildren("crossover", subitems);
					// Check that there are at least 2 phase changes
					if (mech.MeltModel().PhaseChangeCount() < 2){
						throw runtime_error("2 phases expected (Sweep, MechParser::readPhaseTransformation).");
					}
					// Check that there is one crossover
					if (subitems.size() == 1){

						for (m = subitems.begin(); m != subitems.end(); ++m){
							
							std::string name, id_below, id_above;
							name = (*m)->GetAttributeValue("name");
							
							// Get crossover diameter parameters
							double aa(0.0), ab(0.0), ac(0.0), ad(0.0);
							el = (*m)->GetFirstChild("a");
							if (el != NULL) aa = cdble(el->Data());
							el = (*m)->GetFirstChild("b");
							if (el != NULL) ab = cdble(el->Data());
							el = (*m)->GetFirstChild("c");
							if (el != NULL) ac = cdble(el->Data());
							el = (*m)->GetFirstChild("d");
							if (el != NULL) ad = cdble(el->Data());

							// Get phase changes 
							// These point to the phase changes above and below crossover diameter 
							id_above = (*m)->GetAttributeValue("id_above");
							//check phase exists
							if (!mech.MeltModel().PhaseChangeValid(id_above)){
								throw runtime_error(id_above + "Phase change not found "
									"Sweep, MechParser::readPhaseTransformation).");
							}

							id_below = (*m)->GetAttributeValue("id_below");
							//check phase exists
							if (!mech.MeltModel().PhaseChangeValid(id_below)){
								throw runtime_error(id_below + "Phase change not found "
									"Sweep, MechParser::readPhaseTransformation).");
							}

							mech.MeltModel().AddCrossover(name, id_below, id_above, aa, ab, ac, ad);
						}
					}
					else{
						throw runtime_error("1 crossover expected (Sweep, MechParser::readPhaseTransformation).");
					}
				}
				
				// Get oxygen and titanium components for Composition model
				if (type.compare("composition") == 0) {
					std::vector<unsigned int> oxygenindex;
					std::vector<unsigned int> titaniumindex;

					for (unsigned int ii = 0; ii < mech.ComponentCount(); ++ii){
						if (mech.Components()[ii]->Element().compare("O") == 0){
							oxygenindex.push_back(ii);
						}
						else if (mech.Components()[ii]->Element().compare("Ti") == 0){
							titaniumindex.push_back(ii);
						}

					}
					// Set Ti and O indices
					if (oxygenindex.size() > 0){
						mech.MeltModel().SetOxygenComp(oxygenindex);
					}
					else{
						throw runtime_error("Components of element O not found (Sweep, MechParser::readPhaseTransformation).");
					}
					if (titaniumindex.size() > 0){
						mech.MeltModel().SetTitaniumComp(titaniumindex);
					}
					else{
						throw runtime_error("Components of element Ti not found (Sweep, MechParser::readPhaseTransformation).");
					}
				}
			}
			else {
				throw runtime_error("Only Spherical and BinTree models supported for phase transformation (Sweep, MechParser::readPhaseTransformation).");
			}
		}
		//----------------------------------------------------
		//----------- titania_phase_transformation -----------
		//----------------------------------------------------
		// This is a kinetic transformation, which is read like an ordinary process 
		else if (str.compare("kinetic") == 0){

			TitaniaPhaseTransformation *phasetransform = NULL;

			//check that the particle model is valid
			if (mech.AggModel() == AggModels::Spherical_ID || mech.AggModel() == AggModels::BinTree_ID) {
				// Create a new Titania phase transformation object.
				phasetransform = new TitaniaPhaseTransformation(mech);
			}
			else {
				throw runtime_error("Only Spherical and BinTree models supported for phase transformation (Sweep, MechParser::readPhaseTransformation).");
			}

			// Set default name.
			phasetransform->SetName("Phase Transformation " + cstr(k));

			// Read the phase transformation properties.
			try {
				readPhaseTransformation(*(*i), *phasetransform);
			}
			catch (std::exception &e) {
				delete phasetransform;
				throw;
			}

			// Add phase transformation to mechanism.
			mech.AddProcess(*phasetransform);
		}
		else {
			// Unrecognised reaction type.
			throw runtime_error("Unrecognised phase transformation type: " + str +
				" (Sweep, MechParser::readPhaseTransformation).");
		}

		
	}
}

void MechParser::readPhaseTransformation(CamXML::Element &xml, Processes::TitaniaPhaseTransformation &phasetransform)
{
    string str;
    CamXML::Element *el = NULL;

    // Read name.
    str = xml.GetAttributeValue("name");
    if (str != "") phasetransform.SetName(str);

    // Is reaction deferred.
    str = xml.GetAttributeValue("defer");
    if (str=="true") {
    	phasetransform.SetDeferred(true);
    } else {
    	phasetransform.SetDeferred(false);
    }

	/////////////////////////////////////////////////////////////////
	//check that components are correct for phase transformation
	vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;
	bool Ru(false);
	bool An(false);

	// Get list of component changes from XML.
    xml.GetChildren("component", items);

    for (i=items.begin(); i!=items.end(); ++i) {
        // Get component ID.
        string str = (*i)->GetAttributeValue("id");
		if(str.compare("Ru") == 0)	Ru = true;
		if(str.compare("An") == 0)	An = true;
    }
	
	if( Ru != true || An != true) {
        throw runtime_error(" 'Ru' and 'An' components are required for phase transformation "
                                "Sweep, MechParser::readPhaseTransformation).");
    }
	/////////////////////////////////////////////////////////////////

    // Read particle composition change.
    readCompChanges(xml, phasetransform);

    // Read tracker variable changes.
    readTrackChanges(xml, phasetransform);

    //========== Read Arrhenius rate parameters ======================
    Sprog::Kinetics::ARRHENIUS arr;
    el = xml.GetFirstChild("A");
    if (el != NULL) {
        arr.A = cdble(el->Data());
    } else {
        // Reaction must have constant.
        throw runtime_error("Phase transformation found with no rate constant "
                            "defined (Sweep, MechParser::readPhaseTransformation).");
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

        if (str.compare("an")==0) {
            // This reaction depends on surface area.
            phasetransform.SetPropertyID(Sweep::iAn_2_3_comp);
		}
    } else {
        throw runtime_error("Surface process defined without ""particleterm"" "
                            "element (Sweep, MechParser::readPhaseTransformation).");
    }

    phasetransform.SetArrhenius(arr);

}

//COAGULATION

/**
 * @param[in]       xml     XML document containing zero or more top level \<coagulation\> nodes
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
						if (mech.Components()[0]->WeightedPAHs()){
							coag.reset(new Processes::TransitionCoagulationWeightedPAHs(mech));
						}
						else if (mech.IsHybrid())
						{
							coag.reset(new Processes::HybridTransitionCoagulation(mech));							
						}
						else{
							coag.reset(new Processes::TransitionCoagulation(mech));
						}
                    else if(kernelName == "additive")
                        coag.reset(new Processes::AdditiveCoagulation(mech));
                    else if(kernelName == "constant")
						if (mech.IsHybrid())
						{
							coag.reset(new Processes::HybridConstantCoagulation(mech));
						}
						else{
							coag.reset(new Processes::ConstantCoagulation(mech));
						}
					else if (kernelName == "transitionweightedPAHs")
						coag.reset(new Processes::TransitionCoagulationWeightedPAHs(mech));
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
                        else if (choice == "midpoint")
                            coag->SetPositionChoiceRule(Processes::Coagulation::MidpointPositionChoice);
                        else if (choice == "centreofmass")
                            coag->SetPositionChoiceRule(Processes::Coagulation::CentreOfMassPositionChoice);
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
                    else if(kernelName == "weightedtransition")
                         coag.reset(new Processes::WeightedTransitionCoagulation(mech, weightRule));
                    else
                        // Unrecognised option
                        throw std::runtime_error("Coagulation kernel " + kernelName + " not yet available with weights \
                                                (Sweep, MechParser::readCoagulation)");
                }

                // Rate scaling now that a process has been created
                double A = 0.0;
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

//COAGULATION

/**
 * @param[in]       xml     XML document containing zero or more top level \<coagulation\> nodes
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
void MechParser::readFragmentation(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    vector<CamXML::Element*> items;

    // Get list of inceptions from XML data.
    xml.Root()->GetChildren("fragmentation", items);

    for(vector<CamXML::Element*>::const_iterator it = items.begin();
        it != items.end(); ++it) {

        // Read the user choice of kernel.  If no kernel is specified
        // a default will be used
        const CamXML::Element* const kernel = (*it)->GetFirstChild("kernel");
        if(kernel ==  NULL)
        {
            throw std::runtime_error("No kernel given for fragmentation \
                                        (Sweep, MechParser::readFragmentation)");
        }
        else
        {
            Processes::FragWeightRule weightRule;

            // Work out which kernel to use and create it
            const string kernelName = kernel->Data();

            // Create a process of the appropriate type, but wrap it with an auto_ptr so
            // that it gets deleted if an exception is thrown when reading in the value
            // of A.
            std::auto_ptr<Processes::Fragmentation> frag;

            // Unweighted case
            if(kernelName == "erosion") {
                frag.reset(new Processes::ErosionFragmentation(mech));
            } else if(kernelName == "symmetric") {
                frag.reset(new Processes::SymmetricFragmentation(mech));
            } else if(kernelName == "weightederosion") {
                CamXML::Element *weightXML = (*it)->GetFirstChild("weightrule");
                if(weightXML == NULL) {
                    throw std::runtime_error("Fragulation weight rule must be specified for erosion \
                                                (Sweep, MechParser::readFragulation)");
                } else {
                    const std::string weightRuleName = weightXML->Data();
                    if(weightRuleName == "w1" || weightRuleName == "number") {
                        weightRule = Processes::FragWeightNumber;
                        std::cerr << "Found number weight rule" << std::endl;
                    }
                    else if(weightRuleName == "w2" || weightRuleName == "mass") {
                        weightRule = Processes::FragWeightMass;
                        std::cerr << "Found mass weight rule" << std::endl;
                    }
                    else
                        throw std::runtime_error("Fragmentation weight rule " + weightRuleName + " not supported \
                                                    (Sweep, MechParser::readFragmentation)");
                    frag.reset(new Processes::WeightedErosionFragmentation(mech, weightRule));
                }
            } else if(kernelName == "weightedsymmetric") {
                weightRule = Processes::FragWeightSymmetric;
                frag.reset(new Processes::WeightedSymmetricFragmentation(mech, weightRule));
            } else {
                // Unrecognised option
                throw std::runtime_error("Fragmentation kernel " + kernelName + " not yet available with weights \
                                        (Sweep, MechParser::readFragmentation)");
            }

            // Choice of position of newly coagulated particle
            const CamXML::Element *positionChoiceXML = (*it)->GetFirstChild("positionchoice");

            // This is an optional input
            if (positionChoiceXML != NULL) {
                const std::string choice = positionChoiceXML->Data();
                if(choice == "none")
                    frag->SetPositionChoiceRule(Processes::Fragmentation::NoPositionChoice);
                else if (choice == "uniform")
                    frag->SetPositionChoiceRule(Processes::Fragmentation::UniformPositionChoice);
                else if (choice == "mass")
                    frag->SetPositionChoiceRule(Processes::Fragmentation::MassPositionChoice);
                else if (choice == "largestmass")
                    frag->SetPositionChoiceRule(Processes::Fragmentation::LargestMassPositionChoice);
                else if (choice == "midpoint")
                    frag->SetPositionChoiceRule(Processes::Fragmentation::MidpointPositionChoice);
                else if (choice == "centreofmass")
                    frag->SetPositionChoiceRule(Processes::Fragmentation::CentreOfMassPositionChoice);
                else
                    // Unrecognised option
                    throw std::runtime_error("Position choice rule " + choice + " not yet available \
                                                (Sweep, MechParser::readFragmentation)");
            }

            // Rate scaling now that a process has been created
            double A = 0.0;
            CamXML::Element *el = (*it)->GetFirstChild("A");
            if (el != NULL) {
                A = cdble(el->Data());
            } else {
                A = 1.0;
            }
            frag->SetA(A);

            el = (*it)->GetFirstChild("particleterm");
            if (el!=NULL) {
                string str = el->GetAttributeValue("id");
                if (str.compare("size")==0) {
                    frag->SetPropertyID(Sweep::iNumCarbon);
                } else if (str.compare("uniform")==0) {
                    frag->SetPropertyID(Sweep::iFrag);
                } else {
                    throw std::runtime_error("Surface process defined with unrecognised ""particleterm"" "
                                            "(Sweep, MechParser::readFragmentation).");
                }
            } else {
                throw runtime_error("Surface process defined without ""particleterm"" "
                            "element (Sweep, MechParser::readFragmentation).");
            }

            mech.AddFragmentation(*frag);

            // Get rid of the auto_ptr without deleting the coagulation object
            frag.release();
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

    // Get all transport processes.
    xml.Root()->GetChildren("diffusion", items);

    if(!items.empty())
        throw std::runtime_error("Diffusion processes are no longer supported");
}

/*!
 *\param[in]        xml         XML document containing a mechanism description
 *\param[in,out]    mech        Mechanism to which transport processes will be added
 */
void MechParser::readAdvectionProcs(CamXML::Document &xml, Mechanism &mech)
{
    vector<CamXML::Element*> items;
    vector<CamXML::Element*>::iterator i;
    string str;

    // Get all transport processes.
    xml.Root()->GetChildren("advection", items);

    if(!items.empty())
        throw std::runtime_error("Advection processes are no longer supported");
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
        double m = cdble(str);
        if (m > 0.0) {
            mass.push_back(m*1.0e-3); // Convert from g to kg.
        } else {
            // Species can't have zero mass!
            throw runtime_error("Species can't have zero mass "
                                "(Sweep, MechParser::readReactantMDs).");
        }

        // Get species diameter
        str = (*i)->GetAttributeValue("d");
        double d = cdble(str);
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
            double dx = cdble(str);
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
            double dx = cdble(str);
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
