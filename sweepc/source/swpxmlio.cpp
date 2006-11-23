#include "camxml.h"
#include "swpxmlio.h"
#include "swpcomponent.h"
#include "swpinception.h"
#include "swpsurfacereaction.h"
#include "swpactivesitesreaction.h"
#include "swpcondensation.h"
#include "swpparams.h"
#include "swpabf.h"
#include <exception>

using namespace Sweep;

XMLIO::XMLIO(void)
{
}

XMLIO::~XMLIO(void)
{
}

int XMLIO::ReadMechanism(const std::string &filename, Sweep::Mechanism &mech)
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
            return readVersion1(xml, mech);
        } catch (exception &e) {
            throw e;
        }
    } else {
        // Unknown version.
        throw range_error("Unrecognised version attribute in XML file.");
    }
    return 0;
}

int XMLIO::readVersion1(CamXML::Document &xml, Sweep::Mechanism &mech)
{
    CamXML::Element *root = xml.Root();
    vector<CamXML::Element*> items, subitems;
    vector<CamXML::Element*>::iterator i, j;
    string str;
    CamXML::Element *el=NULL;
    unsigned int uint=0;
    int sint=0;
    real num=0.0;
    ACTSURFMODEL acmodel = Const;

    // Read defined models.
    root->GetChildren("model", items);
    for (i=items.begin(); i!=items.end(); i++) {
        str = (*i)->GetAttributeValue("id");
        if (str.compare("actsurf")==0) {
            // Read active surface model type.
            str = (*i)->GetAttributeValue("type");                
            if (str.compare("const")==0) {
                // Constant active surface.
                acmodel = Const;
            } else if (str.compare("abf")==0) {
                // Active sites by ABF model.
                acmodel = ABF;
            } else if (str.compare("profile")==0) {
                // Active sites given in chemistry profile.
                acmodel = Profile;
            }
        }
    }

    // Read components.
    Component *comp;
    root->GetChildren("component", items);
    for (i=items.begin(); i!=items.end(); i++) {
        comp = new Component();

        // Get component id.
        str = (*i)->GetAttributeValue("id");
        if (str != "") {
            comp->SetName(str);
        } else {
            // Component must have an id!
            delete comp;
            throw exception("Component statement found with no ""id"" attribute.");
        }
        
        // Get component density.
        el = (*i)->GetFirstChild("density");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetDensity(atof(str.c_str()));
            } else {
                // Density contains no data.
                delete comp;
                throw exception("Component density contains no data.");
            }
        } else {
            // Component does not have density specified.
            delete comp;
            throw exception("Density required for component specification.");
        }

        // Get component mol. wt.
        el = (*i)->GetFirstChild("molwt");
        if (el!=NULL) {
            str = el->Data();
            if (str != "") {
                comp->SetMolWt(atof(str.c_str()));
            } else {
                // Mol. wt. contains no data.
                delete comp;
                throw exception("Component molecular weight contains no data.");
            }
        } else {
            // Component does not have molecular weight specified.
            delete comp;
            throw exception("Molecular weight required for component specification");
        }

        // Add component to mechanism.
        mech.AddComponent(*comp);
    }

    // Read tracker variable names.
    root->GetChildren("track", items);
    for (i=items.begin(); i!=items.end(); i++) {
        str = (*i)->GetAttributeValue("id");
        if (str != "") {
            mech.AddValueName(str);
        } else {
            // Variable must have a name.
            throw exception("Tracker variable must have an ""id"" attribute.");
        }
    }

    // Read inception processes.
    Inception *icn = NULL;
    real m[2], d[2];
    root->GetChildren("inception", items);
    for (i=items.begin(); i!=items.end(); i++) {
        // Create new inception.
        icn = new Inception();

        // Read reactants.
        (*i)->GetChildren("reactant", subitems);
        for (j=subitems.begin(); j!=subitems.end(); j++) {
            str = (*j)->GetAttributeValue("id");
            sint = mech.GetSpeciesList().GetIndex(str);
            if (sint < 0) {
                // Reactant does not exist in list of species.
                delete icn;
                throw exception(str.append(": Reactant species not found in mechanism.").c_str());
            } else {
                uint = sint;

                // Get reactant stoichiometry.
                str = (*j)->GetAttributeValue("stoich");
                sint = atoi(str.c_str());
                if (sint==0) sint = 1;

                // Get species mass and diameter, required to calculate
                // inception rate.
                str = (*j)->GetAttributeValue("m");
                m[icn->ReactantCount()] = atof(str.c_str());
                if (m[icn->ReactantCount()] == 0.0) {
                    // Species can't have zero mass!
                    delete icn;
                    throw exception("Incepting species can't have zero mass");
                }
                str = (*j)->GetAttributeValue("d");
                d[icn->ReactantCount()] = atof(str.c_str());
                if (d[icn->ReactantCount()] == 0.0) {
                    // Species can't have zero diameter!
                    delete icn;
                    throw exception("Incepting species can't have zero diameter.");
                }

                // Add reactant to inception.
                icn->SetReactant(uint, sint);

                if ((icn->ReactantCount()==1) && (sint==2)) {
                    // Dimer inception.
                    m[1] = m[0];
                    d[1] = d[0];
                    break;
                } else if (icn->ReactantCount()==2) {
                    // Only allowed 2 reactant molecules for inception.
                    break;
                }
            }            
        }
        
        // Finished reading reactants, need to validate and set reacting
        // species' masses and diameters.
        icn->SetInceptingSpecies(m[0], m[1], d[0], d[1]);

        // Read products.
        (*i)->GetChildren("product", subitems);
        for (j=subitems.begin(); j!=subitems.end(); j++) {
            str = (*j)->GetAttributeValue("id");
            sint = mech.GetSpeciesList().GetIndex(str);
            if (sint < 0) {
                // Product does not exist in list of species.
                delete icn;
                throw exception(str.append(": Product species does not exist in mechanism.").c_str());
            } else {
                uint = sint;

                // Get product stoichiometry.
                str = (*j)->GetAttributeValue("stoich");
                sint = atoi(str.c_str());
                if (sint==0) sint = 1;

                // Add product to inception.
                icn->SetProduct(uint, sint);
            }            
        }

        // Read initial particle composition.
        icn->SetComponents(mech.GetComponents());
        (*i)->GetChildren("component", subitems);
        for (j=subitems.begin(); j!=subitems.end(); j++) {
            str = (*j)->GetAttributeValue("id");
            sint = mech.GetComponentIndex(str);
            if (sint >= 0) {
                // Get component change.
                str = (*j)->GetAttributeValue("dx");
                num = atof(str.c_str());
                // Set component change.
                icn->SetCompChange(sint, num);
            } else {
                // Unknown component in mechanism.
                delete icn;
                throw exception(str.append(": Component not found in mechanism.").c_str());
            }
        }

        // Read initial tracker variable values.
        (*i)->GetChildren("track", subitems);
        for (j=subitems.begin(); j!=subitems.end(); j++) {
            str = (*j)->GetAttributeValue("id");
            sint = mech.GetValueIndex(str);
            if (sint >= 0) {
                // Get component change.
                str = (*j)->GetAttributeValue("dx");
                num = atof(str.c_str());
                // Set component change.
                icn->SetValueChange(sint, num);
            } else {
                // Unknown tracker variable in mechanism.
                delete icn;
                throw exception(str.append(": Tracker variable not found in mechanism.").c_str());
            }
        }

        // Add inception to mechanism.  Once entered into mechanism, the mechanism
        // takes control of the inception object for memory management.
        mech.AddInception(icn);
    }

    // Read surface processes.
    SurfaceReaction *surf = NULL;
    ActiveSitesReaction *act = NULL;
    Condensation *cond = NULL;
    Process *rxn = NULL;
    ParticleChanger *pc = NULL;
    RXNTYPE rxntype;

    root->GetChildren("reaction", items);
    for (i=items.begin(); i!=items.end(); i++) {
        // Create correct type of surface reaction.
        str = (*i)->GetAttributeValue("type");
        if (str.compare("surface")==0) {
            el = (*i)->GetFirstChild("particleterm");
            if (el!=NULL) {
                if (el->GetAttributeValue("id").compare("as")==0) {
                    // This is an active sites reaction.
                    act = new ActiveSitesReaction();
                    rxn = act;
                    pc = act;
                    surf = act;
                    rxntype = ActiveSites;

                    // Set active sites function.
                    if (acmodel==ABF) {
                        act->SetActiveSitesFn(Sweep::ABF::ABFMech::HACASites);
                    }
                } else {
                    // This is a normal surface reaction.
                    surf = new SurfaceReaction();
                    rxn = surf;
                    pc = surf;
                    rxntype = Surf;
                }
            } else {
                // Reaction contains no "particleterm" element.  This is an
                // error.
                throw exception("Reaction found with no ""particleterm"" attribute.");
            }
        } else if (str.compare("condensation")==0) {
            // This is a condensation reaction.
            cond = new Condensation();
            rxn = cond;
            pc = cond;
            rxntype = Cond;
        } else {
            // Unrecognised reaction type.
            throw exception(string("Unrecognised reaction type: ").append(str).c_str());
        }

        // Is reaction deferred.
        str = (*i)->GetAttributeValue("defer");
        if (str=="true") {
            rxn->SetDeferred(true);
        } else {
            rxn->SetDeferred(false);
        }

        // Read reactants.
        (*i)->GetChildren("reactant", subitems);
        for (j=subitems.begin(); j!=subitems.end(); j++) {
            str = (*j)->GetAttributeValue("id");
            sint = mech.GetSpeciesList().GetIndex(str);
            if (sint < 0) {
                // Reactant does not exist in list of species.
                if (rxntype==Surf) delete surf;
                if (rxntype==ActiveSites) delete act;
                if (rxntype==Cond) delete cond;
                throw exception(str.append(": Reactant species not found in mechanism.").c_str());
            } else {
                uint = sint;

                // Get reactant stoichiometry.
                str = (*j)->GetAttributeValue("stoich");
                sint = atoi(str.c_str());
                if (sint==0) sint = 1;

                // Add reactant to reaction.
                rxn->SetReactant(uint, sint);

                // If this is a condensation then we also need to get specie's
                // mass and diameter to calculate collision kernel.
                if (rxntype==Cond) {
                    str = (*j)->GetAttributeValue("m");
                    m[0] = atof(str.c_str());
                    if (m[0] == 0.0) {
                        // Species can't have zero mass!
                        delete cond;
                        throw exception("Condensing species can't have zero mass");
                    }
                    str = (*j)->GetAttributeValue("d");
                    d[0] = atof(str.c_str());
                    if (d[0] == 0.0) {
                        // Species can't have zero diameter!
                        delete cond;
                        throw exception("Condensing species can't have zero diameter.");
                    }

                    // Condensation may only have one molecule reactant.
                    if (cond->ReactantCount()==1) {break;}
                }
            }            
        }
        
        // Finished reading reactants, need to validate and set reacting
        // species' masses and diameters if condensation.
        if (rxntype==Cond) {
            cond->SetCondensingSpecies(m[0], d[0]);
        }

        // Read products.
        (*i)->GetChildren("product", subitems);
        for (j=subitems.begin(); j!=subitems.end(); j++) {
            str = (*j)->GetAttributeValue("id");
            sint = mech.GetSpeciesList().GetIndex(str);
            if (sint < 0) {
                // Product does not exist in list of species.
                if (rxntype==Surf) delete surf;
                if (rxntype==ActiveSites) delete act;
                if (rxntype==Cond) delete cond;
                throw exception(str.append(": Product species does not exist in mechanism.").c_str());
            } else {
                uint = sint;

                // Get product stoichiometry.
                str = (*j)->GetAttributeValue("stoich");
                sint = atoi(str.c_str());
                if (sint==0) sint = 1;

                // Add product to reaction.
                rxn->SetProduct(uint, sint);
            }            
        }

        // Read particle composition change.
        icn->SetComponents(mech.GetComponents());
        (*i)->GetChildren("component", subitems);
        for (j=subitems.begin(); j!=subitems.end(); j++) {
            str = (*j)->GetAttributeValue("id");
            sint = mech.GetComponentIndex(str);
            if (sint >= 0) {
                // Get component change.
                str = (*j)->GetAttributeValue("dx");
                num = atof(str.c_str());
                // Set component change.
                pc->SetCompChange(sint, num);
            } else {
                // Unknown component in mechanism.
                if (rxntype==Surf) delete surf;
                if (rxntype==ActiveSites) delete act;
                if (rxntype==Cond) delete cond;
                throw exception(str.append(": Component not found in mechanism.").c_str());
            }
        }

        // Read tracker variable changes.
        (*i)->GetChildren("track", subitems);
        for (j=subitems.begin(); j!=subitems.end(); j++) {
            str = (*j)->GetAttributeValue("id");
            sint = mech.GetValueIndex(str);
            if (sint >= 0) {
                // Get variable change.
                str = (*j)->GetAttributeValue("dx");
                num = atof(str.c_str());
                // Set variable change.
                pc->SetValueChange(sint, num);
            } else {
                // Unknown tracker variable in mechanism.
                if (rxntype==Surf) delete surf;
                if (rxntype==ActiveSites) delete act;
                if (rxntype==Cond) delete cond;
                throw exception(str.append(": Tracker variable not found in mechanism.").c_str());
            }
        }

        // If not a condensation process then also need to get Arrhenius parameters and
        // particle dependency.
        if (rxntype!=Cond) {
            // Arrhenius parameters.
            el = (*i)->GetFirstChild("A");
            if (el!=NULL) {
                num = atof(el->Data().c_str());
            } else {
                // Reaction must have constant.
                throw exception("Surface reaction found with no rate constant defined.");
            }
            el = (*i)->GetFirstChild("n");
            if (el!=NULL) {
                m[0] = atof(el->Data().c_str());
            } else {
                // Default temperature power is 0.
                m[0] = 0.0;
            }
            el = (*i)->GetFirstChild("E");
            if (el!=NULL) {
                m[1] = atof(el->Data().c_str()) * R / RCAL;
            } else {
                // Default activation energy is zero.
                m[1] = 0.0;
            }
            surf->SetArrhenius(num, m[0], m[1]);
 
            // Particle dependency.
            el = (*i)->GetFirstChild("particleterm");
            if (el!=NULL) {
                str = el->GetAttributeValue("id");
                sint = atoi(el->GetAttributeValue("power").c_str());
                if ((str.compare("s")==0) || (str.compare("as")==0)) {
                    // This reaction depends on surface area.  Ignore power,
                    // they must have meant 1.
                    surf->SetPropertyID(DefaultParticle::iS);
                } else if (str.compare("d")==0) {
                    // This reaction depends on some power of the diameter.
                    switch (sint) {
                        case 1:
                            surf->SetPropertyID(DefaultParticle::iD);
                            break;
                        case 2:
                            surf->SetPropertyID(DefaultParticle::iD2);
                            break;
                        case -1:
                            surf->SetPropertyID(DefaultParticle::iD_1);
                        case -2:
                            surf->SetPropertyID(DefaultParticle::iD_2);
                        case 0:
                            // Oh dear, can't have a zero power.
                            throw exception("""particleterm"" tag found with invalid or zero power attribute");
                    }
                }
            } else {
                throw exception("Surface process defined with ""particleterm"" element.");
            }
        } else {
            // Read only rate constant for condensations.
            el = (*i)->GetFirstChild("A");
            if (el!=NULL) {
                cond->SetA(atof(el->Data().c_str()));
            } else {
                // Reaction must have constant.
                throw exception("Surface reaction found with no rate constant defined.");
            }
        }

         // Add process to mechanism.
        if (rxntype==Surf) {
            mech.AddProcess(surf);
        } else if (rxntype==ActiveSites) {
            mech.AddProcess(act);
        } else if (rxntype==Cond) {
            mech.AddProcess(cond);
        }
    }
    return 0;
}