/*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    The element flux analysis class post-processes the forward and reverse reaction rates
    for FluxViewer.

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

#include "mops_flux_postprocessor.h"
#include "string_functions.h"

FluxAnalyser::FluxAnalyser(const Mechanism &mech, 
                           const timevector &times, 
                           const std::vector<fvector> &agpfwdrates,
                           const std::vector<fvector> &agprevrates,
                           const fvector &atemperatures) {
    m_mech = &mech;
    m_agpfwdrates = &agpfwdrates;
    m_agprevrates = &agprevrates;
    m_atemperatures = &atemperatures;
    // Make all time points.
    for (unsigned int i = 0; i < times.size(); i++) {
        for (unsigned int j = 0; j <= times.at(i).StepCount(); j++) {
            real t = times.at(i).StartTime();
            t += ((real) j) * times.at(i).StepSize();
            m_times.push_back(t);
        }
    }
}

void FluxAnalyser::addFluxElement(const Sprog::Element &elem) {
    int iel = m_mech->FindElement(elem);
    if (iel > -1) {
        m_ElementIndexes.push_back((unsigned int)iel);
    } else {
        throw runtime_error("Attemp to add unknow element in Flux analysis "
                            "(Mops, FluxAnalyser::addFluxElement).");
    }
}

int FluxAnalyser::getFluxPathIndex(Mops::FluxAnalyser::FluxPath &fpath) {
    for (unsigned int i = 0; i < m_flux_network.size(); i++) {
        if ((m_flux_network.at(i).SourceSpecies == fpath.SourceSpecies) &&
            (m_flux_network.at(i).TargetSpecies == fpath.TargetSpecies)) {
            return (int) i;
        }
    }
    return -1;
}

void FluxAnalyser::addToFluxPathRate(Mops::FluxAnalyser::FluxPath &fpath) {
    int index = getFluxPathIndex(fpath);
    if (index < 0) {
        m_flux_network.push_back(fpath);
    } else {
        m_flux_network.at(index).Rate += fpath.Rate;
    }
}

real FluxAnalyser::getTotalElementStoi(const Sprog::Kinetics::Reaction &rxn, unsigned int iel) {
    real n = 0.0;
    for (int i = 0; i < rxn.ReactantCount(); i++) {
        n += getNumberOfElementAtom(rxn.Reactant(i), iel);
    }
    for (int i = 0; i < rxn.FReactantCount(); i++) {
        n += getNumberOfElementAtom(rxn.FReactant(i), iel);
    }
    return n;
}

real FluxAnalyser::getNumberOfElementAtom(const Sprog::Stoich &sc, unsigned int iel) {
    if (m_mech->GetSpecies(sc.Index()) != NULL) {
        real n = (real) m_mech->GetSpecies(sc.Index())->AtomCount(iel);
        return n * ((real) sc.Mu());
    } else {
        return 0.0;
    }
}

real FluxAnalyser::getNumberOfElementAtom(const Sprog::Stoichf &sc, unsigned int iel) {
    if (m_mech->GetSpecies(sc.Index()) != NULL) {
        real n = (real) m_mech->GetSpecies(sc.Index())->AtomCount(iel);
        return n * ((real) sc.Mu());
    } else {
        return 0.0;
    }
}
void FluxAnalyser::writeFluxAt(unsigned int index, unsigned int iel, std::ofstream &fout) {
    // If index of time point is greater than total number of time point then do not do flux analysis
    if (index >= m_times.size()) {
        return;
    }
    if (fout.is_open()) {
        calculateFluxAt(index, iel);
        // Write flux network to file.
        if (m_flux_network.size() > 0) {
            fout << "Fluxes of " << m_mech->Elements(iel)->Name() << endl;
            for (unsigned int i = 0; i < m_flux_network.size(); i++) {
                fout << formatWhiteSpace(ComoString::int2string(i+1),5,false) << "  "
                     << formatWhiteSpace(m_mech->GetSpecies(m_flux_network.at(i).SourceSpecies)->Name(), 17)
                     << "=> "
                     << formatWhiteSpace(m_mech->GetSpecies(m_flux_network.at(i).TargetSpecies)->Name(), 19)
                     << formatWhiteSpace(Strings::cstr(m_flux_network.at(i).Rate), 12, false)
                     << endl;
            }
            fout << endl;
        }
    } else {
        throw runtime_error("File stream is not opened for writting flux analysis "
                            "(Mops, FluxAnalyser::writeFluxAt).");
    }
}

void FluxAnalyser::writeFluxes(const std::string &filename) {
    //  write output to file
    if (m_ElementIndexes.size() > 0) {
        ofstream fout;
        fout.open(filename.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);
        cout << "mops: Writting output from FluxAnalyser..." << endl;
        if (fout.is_open()) {
            // Write header
            fout << formatWhiteSpace(ComoString::int2string(m_mech->SpeciesCount()),5,false) << "/ Number of species" << endl;
            fout << formatWhiteSpace(ComoString::int2string(m_mech->ReactionCount()),5,false) << "/ Number of reactions" << endl;
            fout << formatWhiteSpace(ComoString::int2string(m_times.size() - 1),5,false) << "/ Number of points" << endl;
            for (unsigned int i = 0; i < m_ElementIndexes.size(); i++) {
                fout << formatWhiteSpace(m_mech->Elements(m_ElementIndexes.at(i))->Name(), 4);
            }
            fout << endl;
            for (unsigned int i = 0; i < m_mech->SpeciesCount(); i++) {
                fout << formatWhiteSpace(m_mech->GetSpecies(i)->Name(), 16);
            }
            fout << endl;
            fout << endl;
            for (unsigned int i = 1; i < m_times.size(); i++) {
                fout << formatWhiteSpace(ComoString::int2string(i),5,false) << "  "
                     << m_times.at(i) << " t/s  " << m_atemperatures->at(i) << " T/K " << endl;
                for (unsigned int j = 0; j < m_ElementIndexes.size(); j++) {
                    writeFluxAt(i, m_ElementIndexes.at(j), fout);
                }
            }
            fout.close();
        } else {
            // Throw error if the output file failed to open.
            throw runtime_error("Failed to open file for write element flux "
                                "output (Mops, Simulator::FluxAnalyser::writeFluxes).");
        }
    }
}

std::string FluxAnalyser::formatWhiteSpace(std::string str, unsigned int len, bool isright) {
    string wstr = "";
    if (str.length() <= len) {
        for (unsigned int i = 0; i < len - str.length(); i++) {
            wstr += " ";
        }
    } else {
        throw runtime_error("String is too long to output in flux analysis "
                            "(Mops, FluxAnalyser::formatWhiteSpace).");
    }
    if (isright) {
        wstr = str.append(wstr);
    } else {
        wstr = wstr.append(str);
    }
    return wstr;
}

void FluxAnalyser::calculateFluxAt(unsigned int index, unsigned int iel) {
    m_flux_network.clear();
    // If index of time point is greater than total number of time point then do not do flux analysis
    //if (index >= m_times-) {
    //    return;
    //}
    for (unsigned int i = 0; i < m_mech->ReactionCount(); i++) {
        if (m_mech->Reactions(i)->ReactantCount() > 0) {
            // Integer Stoichiometry
            for (int j = 0; j < m_mech->Reactions(i)->ReactantCount(); j++) {
                for (int k = 0; k < m_mech->Reactions(i)->ProductCount(); k++) {
                    real n_total_stoi = getTotalElementStoi(*m_mech->Reactions(i), iel);
                    real n_A_elem     = getNumberOfElementAtom(m_mech->Reactions(i)->Reactant(j), iel);
                    real n_B_elem     = getNumberOfElementAtom(m_mech->Reactions(i)->Product(k), iel);
                    real flux_fraction = 0.0;
                    if (n_total_stoi > 0) {
                        flux_fraction = n_A_elem * n_B_elem / n_total_stoi;
                    } else {
                        flux_fraction = 0.0;
                    }
                    // A => B Flux
                    FluxPath temp_fpath;
                    temp_fpath.Rate = (m_agpfwdrates->at(index)).at(i) * flux_fraction;
                    if (temp_fpath.Rate > 0.0) {
                        temp_fpath.SourceSpecies = m_mech->Reactions(i)->Reactant(j).Index();
                        temp_fpath.TargetSpecies = m_mech->Reactions(i)->Product(k).Index();
                        addToFluxPathRate(temp_fpath);
                    }

                    // B => A Flux, this is only if reaction is reversible.
                    if (m_mech->Reactions(i)->IsReversible()) {
                        temp_fpath.Rate = (m_agprevrates->at(index)).at(i) * flux_fraction;
                        if (temp_fpath.Rate > 0.0) {
                            temp_fpath.SourceSpecies = m_mech->Reactions(i)->Product(k).Index();
                            temp_fpath.TargetSpecies = m_mech->Reactions(i)->Reactant(j).Index();
                            addToFluxPathRate(temp_fpath);
                        }
                    }
                }
            }
        } else {
            // Real Stoichiometry
            for (int j = 0; j < m_mech->Reactions(i)->FReactantCount(); j++) {
                for (int k = 0; k < m_mech->Reactions(i)->FProductCount(); k++) {
                    real n_total_stoi = getTotalElementStoi(*m_mech->Reactions(i), iel);
                    real n_A_elem     = getNumberOfElementAtom(m_mech->Reactions(i)->FReactant(j), iel);
                    real n_B_elem     = getNumberOfElementAtom(m_mech->Reactions(i)->FProduct(k), iel);
                    real flux_fraction = n_A_elem * n_B_elem / n_total_stoi;
                    // A => B Flux
                    FluxPath temp_fpath;
                    temp_fpath.SourceSpecies = m_mech->Reactions(i)->FReactant(j).Index();
                    temp_fpath.TargetSpecies = m_mech->Reactions(i)->FProduct(k).Index();
                    temp_fpath.Rate          = (m_agpfwdrates->at(index)).at(i) * flux_fraction;
                    addToFluxPathRate(temp_fpath);

                    // B => A Flux, this is only if reaction is reversible.
                    if (m_mech->Reactions(i)->IsReversible()) {
                        temp_fpath.SourceSpecies = m_mech->Reactions(i)->FProduct(k).Index();
                        temp_fpath.TargetSpecies = m_mech->Reactions(i)->FReactant(j).Index();
                        temp_fpath.Rate          = (m_agprevrates->at(index)).at(i) * flux_fraction;
                        addToFluxPathRate(temp_fpath);
                    }
                }
            }
        }
    }

    // Sort current Flux Network
    sortDecendingFluxNetwork();
}

void FluxAnalyser::sortDecendingFluxNetwork() {
    unsigned int n = m_flux_network.size();
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = i; j > 0; j--) {
            if (m_flux_network[j-1].Rate < m_flux_network[j].Rate) {
                FluxPath temp;
                // Store j in temp.
                temp.Rate = m_flux_network[j].Rate;
                temp.SourceSpecies  = m_flux_network[j].SourceSpecies;
                temp.TargetSpecies  = m_flux_network[j].TargetSpecies;
                // Put j+1 into j.
                m_flux_network[j].Rate = m_flux_network[j-1].Rate;
                m_flux_network[j].SourceSpecies  = m_flux_network[j-1].SourceSpecies;
                m_flux_network[j].TargetSpecies  = m_flux_network[j-1].TargetSpecies;
                // Put temp into j+1.
                // Put temp into j+1.
                m_flux_network[j-1].Rate = temp.Rate;
                m_flux_network[j-1].SourceSpecies  = temp.SourceSpecies;
                m_flux_network[j-1].TargetSpecies  = temp.TargetSpecies;
            } else {
                break;
            }
        }
    }

}

