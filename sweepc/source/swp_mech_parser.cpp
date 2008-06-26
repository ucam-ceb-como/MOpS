#include "swp_mech_parser.h"
#include "swp_component.h"
#include "swp_tracker.h"
#include "swp_inception.h"
#include "swp_surface_reaction.h"
#include "swp_actsites_reaction.h"
#include "swp_condensation.h"
#include "swp_abf_model.h"
#include "camxml.h"
#include "string_functions.h"
#include <stdexcept>
#include <string>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;
using namespace Strings;

void MechParser::Read(const std::string &filename, Sweep::Mechanism &mech)
{
    CamXML::Document xml;

    // Try to load XML file into data structure.
    try {
        xml.Load(filename);
    } catch (exception &e) {;
        throw e;
    }

    // Read the mechanism file version from document, then parse
    // the correct version.
    string version = xml.Root()->GetAttributeValue("version");
    if ((version.compare("1")==0) || (version=="")) {
        try {
            return readV1(xml, mech);
        } catch (exception &e) {
            throw e;
        }
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

    // Read the processes (inceptions, surface reactions and condensations).
    readInceptions(xml, mech);
    readSurfRxns(xml, mech);
    readCondensations(xml, mech);

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
    } else {
        // No sintering model defined.
        mech.SintModel().Disable();
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
                delete comp;
                throw runtime_error("Component " + comp->Name() + " density "
                                    "contains no data (Sweep, MechParser::readComponents).");
            }
        } else {
            // Component does not have density specified.
            delete comp;
            throw runtime_error("Density required for component " + comp->Name() + 
                                " specification (Sweep, MechParser::readComponents).");
        }

        // Get component mol. wt.
        el = (*i)->GetFirstChild("molwt");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetMolWt(cdble(str)*1.0e-3); // Convert from g/mol to kg/mol.
            } else {
                // Mol. wt. contains no data.
                delete comp;
                throw runtime_error("Component " + comp->Name() + " mol. wt. "
                                    "contains no data (Sweep, MechParser::readComponents).");
            }
        } else {
            // Component does not have molecular weight specified.
            delete comp;
            throw runtime_error("Mol. wt. required for component " + comp->Name() + 
                                " specification (Sweep, MechParser::readComponents).");
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

    int id  = -1;
    real dx = 0.0;
    
    // Get list of inceptions from XML data.
    xml.Root()->GetChildren("inception", items);

    for (i=items.begin(); i!=items.end(); ++i) {
        // Create new inception.
        Inception *icn = new Inception(mech);
        icn->SetMechanism(mech);

        // Read reactants.
        try {
            // Get reactant stoichiometry.
            readReactants(*(*i), *icn);

            // Get reactant masses and diameters, and set inception
            // parameters.
            fvector mass, diam;
            readReactantMDs(*(*i), mass, diam);
            if (mass.size() > 1) {
                icn->SetInceptingSpecies(mass[0], mass[1], diam[0], diam[1]);
            } else if (mass.size() > 0) {
                // Dimer inception.
                icn->SetInceptingSpecies(mass[0], mass[0], diam[0], diam[0]);
            }
        } catch (exception &e) {
            delete icn;
            throw e;
        }

        // Read products.
        try {
            readProducts(*(*i), *icn);
        } catch (exception &e) {
            delete icn;
            throw e;
        }

        // Read initial particle composition.
        (*i)->GetChildren("component", subitems);
        for (j=subitems.begin(); j!=subitems.end(); ++j) {
            // Get component ID.
            str = (*j)->GetAttributeValue("id");
            id = mech.ComponentIndex(str);

            if (id >= 0) {
                // Get component change.
                str = (*j)->GetAttributeValue("dx");
                dx = cdble(str);
                // Set component change.
                icn->SetParticleComp(id, dx);
            } else {
                // Unknown component in mechanism.
                delete icn;
                throw runtime_error(str + ": Component not found in mechanism "
                                    "(Sweep, MechParser::readInceptions).");
            }
        }

        // Read initial tracker variable values.
        (*i)->GetChildren("track", subitems);
        for (j=subitems.begin(); j!=subitems.end(); j++) {
            // Get tracker ID.
            str = (*j)->GetAttributeValue("id");
            id = mech.GetTrackerIndex(str);

            if (id >= 0) {
                // Get tracker change.
                str = (*j)->GetAttributeValue("dx");
                dx = cdble(str);
                // Set tracker change.
                icn->SetParticleTracker(id, dx);
            } else {
                // Unknown tracker variable in mechanism.
                delete icn;
                throw runtime_error(str + ": Tracker variable not found in mechanism. "
                                    "(Sweep, MechParser::readInceptions).");
            }
        }

        // Add inception to mechanism.  Once entered into mechanism, the mechanism
        // takes control of the inception object for memory management.
        mech.AddInception(*icn);
    }
}

// Reads surface reactions from a sweep mechanism XML file.
void MechParser::readSurfRxns(CamXML::Document &xml, Mechanism &mech)
{
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    CamXML::Element *el;
    string str;

    // Get all surface reactions.
    xml.Root()->GetChildren("reaction", items);

    for (i=items.begin(); i!=items.end(); ++i) {
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

        // Is reaction deferred.
        str = (*i)->GetAttributeValue("defer");
        if (str=="true") {
            rxn->SetDeferred(true);
        } else {
            rxn->SetDeferred(false);
        }

        // Read reactants.
        try {
            // Get reactant stoichiometry.
            readReactants(*(*i), *rxn);
        } catch (exception &e) {
            delete rxn;
            throw e;
        }

        // Read products.
        try {
            readProducts(*(*i), *rxn);
        } catch (exception &e) {
            delete rxn;
            throw e;
        }

        // Read particle composition change.
        try {
            readCompChanges(*(*i), *rxn);
        } catch (exception &e) {
            delete rxn;
            throw e;
        }

        // Read tracker variable changes.
        try {
            readTrackChanges(*(*i), *rxn);
        } catch (exception &e) {
            delete rxn;
            throw e;
        }

        // Read Arrhenius rate parameters.
        Sprog::Kinetics::ARRHENIUS arr;
        el = (*i)->GetFirstChild("A");
        if (el != NULL) {
            arr.A = cdble(el->Data());
            // Must scale rate constant from cm3 to m3.
            arr.A *= (1.0e-6 * (real)rxn->ReactantCount());
        } else {
            // Reaction must have constant.
            throw runtime_error("Surface reaction found with no rate constant "
                                "defined (Sweep, MechParser::readSurfRxns).");
        }
        el = (*i)->GetFirstChild("n");
        if (el!=NULL) {
            arr.n = cdble(el->Data());
        } else {
            // Default temperature power is 0.
            arr.n = 0.0;
        }
        el = (*i)->GetFirstChild("E");
        if (el!=NULL) {
            arr.E = cdble(el->Data()) * R / RCAL;
        } else {
            // Default activation energy is zero.
            arr.E = 0.0;
        }
        rxn->SetArrhenius(arr);
 
        // Particle dependency.
        el = (*i)->GetFirstChild("particleterm");
        if (el!=NULL) {
            // Get property ID.
            str = el->GetAttributeValue("id");

            // Get power.
            int id = (int)cdble(el->GetAttributeValue("power"));

            if ((str.compare("s")==0) || (str.compare("as")==0)) {
                // This reaction depends on surface area.  Ignore power,
                // they must have meant 1.
                rxn->SetPropertyID(ParticleCache::iS);
            } else if (str.compare("d")==0) {
                // This reaction depends on some power of the diameter.
                switch (id) {
                    case 1:
                        rxn->SetPropertyID(ParticleCache::iD);
                        break;
                    case 2:
                        rxn->SetPropertyID(ParticleCache::iD2);
                        break;
                    case -1:
                        rxn->SetPropertyID(ParticleCache::iD_1);
                        break;
                    case -2:
                        rxn->SetPropertyID(ParticleCache::iD_2);
                        break;
                    default:
                        // Oh dear, can't have a zero power.
                        throw runtime_error("""particleterm"" tag found with "
                                            "invalid or zero power attribute "
                                            "(Sweep, MechParser::readSurfRxns)");
                }
            }
        } else {
            throw runtime_error("Surface process defined with ""particleterm"" "
                                "element (Sweep, MechParser::readSurfRxns).");
        }

        // Add process to mechanism.
        mech.AddProcess(*rxn);
    }
}

// Reads condensation processes from a sweep mechanism XML file.
void MechParser::readCondensations(CamXML::Document &xml, Mechanism &mech)
{
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    CamXML::Element *el;
    string str;

    // Get all surface reactions.
    xml.Root()->GetChildren("condensation", items);

    for (i=items.begin(); i!=items.end(); ++i) {
        // Create a new condensation object.
        Condensation *cond = new Condensation(mech);

        // Is condensation deferred.
        str = (*i)->GetAttributeValue("defer");
        if (str=="true") {
            cond->SetDeferred(true);
        } else {
            cond->SetDeferred(false);
        }

        // Read reactants.
        try {
            // Get reactant stoichiometry.
            readReactants(*(*i), *cond);

            // Get reactant masses and diameters, and set inception
            // parameters.
            fvector mass, diam;
            readReactantMDs(*(*i), mass, diam);
            if (mass.size() > 0) {
                cond->SetCondensingSpecies(mass[0], diam[0]);
            }
        } catch (exception &e) {
            delete cond;
            throw e;
        }

        // Read products.
        try {
            readProducts(*(*i), *cond);
        } catch (exception &e) {
            delete cond;
            throw e;
        }

        // Read particle composition change.
        try {
            readCompChanges(*(*i), *cond);
        } catch (exception &e) {
            delete cond;
            throw e;
        }

        // Read tracker variable changes.
        try {
            readTrackChanges(*(*i), *cond);
        } catch (exception &e) {
            delete cond;
            throw e;
        }

        // Read Arrhenius rate parameters.
        real A = 0.0;
        el = (*i)->GetFirstChild("A");
        if (el != NULL) {
            A = cdble(el->Data());
        } else {
            A = 1.0;
        }
        cond->SetA(A);

        // Add condensation to mechanism.
        mech.AddProcess(*cond);
    }
}


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
