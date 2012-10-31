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
#include <stdexcept>
FluxAnalyser::FluxAnalyser(const Mechanism &mech,
                           const timevector &times,
                           const std::vector<fvector> &agpfwdrates,
                           const std::vector<fvector> &agprevrates,
                           const fvector &atemperatures) {
    m_mech = &mech;
    m_agpfwdrates = &agpfwdrates;
    m_agprevrates = &agprevrates;
    m_atemperatures = &atemperatures;
    // Make all time points and time stop points.
    double t = times.at(0).StartTime();
    m_times.push_back(t);
    for (unsigned int i = 0; i < times.size(); i++) {
        for (unsigned int j = 1; j <= times.at(i).StepCount(); j++) {
            t += times.at(i).StepSize();
            m_times.push_back(t);
        }
        m_times_stop.push_back(m_times.size()-1);
    }
}

void FluxAnalyser::addElement(const Sprog::Element &elem) {
    int iel = m_mech->GasMech().FindElement(elem);
    if (iel > -1) {
        // Check whether index already in the index list
        bool ExistIndex = false;
        for (unsigned int i = 0; i < m_ElementIndexes.size(); i++) {
            if (m_ElementIndexes.at(i) == static_cast<unsigned int>(iel)) {
                ExistIndex = true;
                break;
            }
        }
        if (!ExistIndex) {
            m_ElementIndexes.push_back((unsigned int)iel);
            FluxNetwork fnw;
            m_flux_networks.push_back(fnw);
            m_int_flux_networks.push_back(fnw);
        }
    } else {
        std::cout << "Element " << elem.Name() << " Not Found. Ignoring "
                  << elem.Name() << " from flux analyzer." << std::endl;
        //throw std::runtime_error("Attemp to add unknow element in Flux analysis "
        //                    "(Mops, FluxAnalyser::addFluxElement).");
    }
}

void FluxAnalyser::writeFluxes(const std::string &filenameprefix, bool doIntFluxes) {
    //  write output to file
    if (m_ElementIndexes.size() > 0) {
        std::cout << "mops: Writting output from FluxAnalyser..." << std::endl;

        std::ofstream fflux;
        std::ofstream fintflux;
        std::string file_flux = filenameprefix + "-elem-flux.fvr";
        std::string file_int_flux = filenameprefix + "-elem-int-flux.fvr";
        fflux.open(file_flux.c_str(), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
        writeHeader(fflux, m_times.size()-1);
        if (doIntFluxes) {
            fintflux.open(file_int_flux.c_str(), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
            writeHeader(fintflux, m_times_stop.size());
        }
        // std::ofstream of fflux and fintflux has been checked in writeHeaders
        int int_flux_time_stop_index = 0;
        for (unsigned int i = 1; i < m_times.size(); i++) {
            unsigned int i_times = m_times_stop.at(int_flux_time_stop_index);
            bool isAtTimeStop = (i_times <= i);

            fflux << formatWhiteSpace(ComoString::int2string(i),5,false) << "  "
                 << m_times.at(i) << " t/s  " << m_atemperatures->at(i) << " T/K " << std::endl;

            if (doIntFluxes && isAtTimeStop) {
                fintflux << formatWhiteSpace(ComoString::int2string(int_flux_time_stop_index+1),5,false) << "  "
                    << m_times.at(i_times) << " t/s  " << m_atemperatures->at(i_times) << " T/K " << std::endl;
            }

            for (unsigned int j = 0; j < m_ElementIndexes.size(); j++) {
                calculateFluxAt(i, m_ElementIndexes.at(j), m_flux_networks.at(j));
                writeFluxAt(m_ElementIndexes.at(j), fflux, m_flux_networks.at(j));
                if (doIntFluxes) {
                    calculateIntFluxAt(i, m_flux_networks.at(j), m_int_flux_networks.at(j));
                    if (isAtTimeStop) {
                        writeFluxAt(m_ElementIndexes.at(j), fintflux, m_int_flux_networks.at(j));
                    }
                }
            }
            // increase the index of the time stop if time point pass the current index.
            if (isAtTimeStop) int_flux_time_stop_index++;
        }
        fflux.close();
        if (doIntFluxes) fintflux.close();
    }
}

void FluxAnalyser::calculateFluxAt(unsigned int index, unsigned int iel, Mops::FluxAnalyser::FluxNetwork &flux_network) {
    flux_network.clear();

    for (unsigned int i = 0; i < m_mech->GasMech().ReactionCount(); i++) {
            // Integer Stoichiometry
            for (int j = 0; j < m_mech->GasMech().Reactions(i)->ReactantCount(); j++) {
                for (int k = 0; k < m_mech->GasMech().Reactions(i)->ProductCount(); k++) {
                    double n_total_stoi = getTotalElementStoi(*m_mech->GasMech().Reactions(i), iel);
                    double n_A_elem     = getNumberOfElementAtom(m_mech->GasMech().Reactions(i)->Reactant(j), iel);
                    double n_B_elem     = getNumberOfElementAtom(m_mech->GasMech().Reactions(i)->Product(k), iel);
                    double flux_fraction = 0.0;
                    if (n_total_stoi > 0) {
                        flux_fraction = n_A_elem * n_B_elem / n_total_stoi;
                    } else {
                        flux_fraction = 0.0;
                    }
                    // A => B Flux
                    FluxPath temp_fpath;
                    temp_fpath.Rate = (m_agpfwdrates->at(index)).at(i) * flux_fraction;
                    if (temp_fpath.Rate > 0.0) {
                        temp_fpath.SourceSpecies = m_mech->GasMech().Reactions(i)->Reactant(j).Index();
                        temp_fpath.TargetSpecies = m_mech->GasMech().Reactions(i)->Product(k).Index();
                        addToFluxPathRate(temp_fpath, flux_network);
                    }

                    // B => A Flux, this is only if reaction is reversible.
                    if (m_mech->GasMech().Reactions(i)->IsReversible()) {
                        temp_fpath.Rate = (m_agprevrates->at(index)).at(i) * flux_fraction;
                        if (temp_fpath.Rate > 0.0) {
                            temp_fpath.SourceSpecies = m_mech->GasMech().Reactions(i)->Product(k).Index();
                            temp_fpath.TargetSpecies = m_mech->GasMech().Reactions(i)->Reactant(j).Index();
                            addToFluxPathRate(temp_fpath, flux_network);
                        }
                    }
                }
            }
    }

    // Sort current Flux Network
    sortDecendingFluxNetwork(flux_network);
}

void FluxAnalyser::calculateIntFluxAt(unsigned int index, Mops::FluxAnalyser::FluxNetwork &flux_network, Mops::FluxAnalyser::FluxNetwork &int_flux_network) {
    for (unsigned int j = 0; j < flux_network.size(); j++) {
        // This is over estimated fluxes. Need modification for proper integration.
        flux_network.at(j).Rate = flux_network.at(j).Rate * (m_times.at(index) - m_times.at(index - 1));
        addToFluxPathRate(flux_network.at(j), int_flux_network);
    }
    sortDecendingFluxNetwork(int_flux_network);
}

int FluxAnalyser::getFluxPathIndex(Mops::FluxAnalyser::FluxPath &fpath, Mops::FluxAnalyser::FluxNetwork &flux_network) {
    for (unsigned int i = 0; i < flux_network.size(); i++) {
        if ((flux_network.at(i).SourceSpecies == fpath.SourceSpecies) &&
            (flux_network.at(i).TargetSpecies == fpath.TargetSpecies)) {
            return (int) i;
        }
    }
    return -1;
}

void FluxAnalyser::addToFluxPathRate(Mops::FluxAnalyser::FluxPath &fpath, Mops::FluxAnalyser::FluxNetwork &flux_network) {
    int index = getFluxPathIndex(fpath, flux_network);
    if (index < 0) {
        flux_network.push_back(fpath);
    } else {
        flux_network.at(index).Rate += fpath.Rate;
    }
}

void FluxAnalyser::sortDecendingFluxNetwork(Mops::FluxAnalyser::FluxNetwork &flux_network) {
    unsigned int n = flux_network.size();
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = i; j > 0; j--) {
            if (flux_network[j-1].Rate < flux_network[j].Rate) {
                FluxPath temp;
                // Store j in temp.
                temp.Rate = flux_network[j].Rate;
                temp.SourceSpecies  = flux_network[j].SourceSpecies;
                temp.TargetSpecies  = flux_network[j].TargetSpecies;
                // Put j+1 into j.
                flux_network[j].Rate = flux_network[j-1].Rate;
                flux_network[j].SourceSpecies  = flux_network[j-1].SourceSpecies;
                flux_network[j].TargetSpecies  = flux_network[j-1].TargetSpecies;
                // Put temp into j+1.
                // Put temp into j+1.
                flux_network[j-1].Rate = temp.Rate;
                flux_network[j-1].SourceSpecies  = temp.SourceSpecies;
                flux_network[j-1].TargetSpecies  = temp.TargetSpecies;
            } else {
                break;
            }
        }
    }

}

void FluxAnalyser::writeFluxAt(unsigned int iel, std::ofstream &fout, Mops::FluxAnalyser::FluxNetwork &flux_network) {
    if (fout.is_open()) {
        // Write flux network to file.
        if (flux_network.size() > 0) {
            fout << "Fluxes of " << m_mech->GasMech().Elements(iel)->Name() << std::endl;
            for (unsigned int i = 0; i < flux_network.size(); i++) {
                fout << formatWhiteSpace(ComoString::int2string(i+1),5,false) << "  "
                     << formatWhiteSpace(m_mech->GasMech().GetSpecies(flux_network.at(i).SourceSpecies)->Name(), 17)
                     << "=> "
                     << formatWhiteSpace(m_mech->GasMech().GetSpecies(flux_network.at(i).TargetSpecies)->Name(), 19)
                     << formatWhiteSpace(Strings::cstr(flux_network.at(i).Rate), 12, false)
                     << std::endl;
            }
            fout << std::endl;
        }
    } else {
        throw std::runtime_error("File stream is not opened for writting flux analysis "
                            "(Mops, FluxAnalyser::writeFluxAt).");
    }
}

void FluxAnalyser::writeHeader(std::ofstream &fout, unsigned int npoints) {
    if (fout.is_open()) {
        // Write header
        fout << formatWhiteSpace(ComoString::int2string(m_mech->GasMech().SpeciesCount()),5,false) << "/ Number of species" << std::endl;
        fout << formatWhiteSpace(ComoString::int2string(m_mech->GasMech().ReactionCount()),5,false) << "/ Number of reactions" << std::endl;
        fout << formatWhiteSpace(ComoString::int2string(npoints),5,false) << "/ Number of points" << std::endl;
        for (unsigned int i = 0; i < m_ElementIndexes.size(); i++) {
            fout << formatWhiteSpace(m_mech->GasMech().Elements(m_ElementIndexes.at(i))->Name(), 4);
        }
        fout << std::endl;
        for (unsigned int i = 0; i < m_mech->GasMech().SpeciesCount(); i++) {
            fout << formatWhiteSpace(m_mech->GasMech().GetSpecies(i)->Name(), 16);
        }
        fout << std::endl;
        fout << std::endl;
    } else {
        // Throw error if the output file failed to open.
        throw std::runtime_error("Failed to open file for write element flux "
                            "output (Mops, Simulator::FluxAnalyser::writeHeader).");
    }

}

double FluxAnalyser::getTotalElementStoi(const Sprog::Kinetics::Reaction &rxn, unsigned int iel) {
    double n = 0.0;
    for (int i = 0; i < rxn.ReactantCount(); i++) {
        n += getNumberOfElementAtom(rxn.Reactant(i), iel);
    }
    return n;
}

double FluxAnalyser::getNumberOfElementAtom(const Sprog::Stoich &sc, unsigned int iel) {
    Sprog::Species * sp = m_mech->GasMech().GetSpecies(sc.Index());
    if (sp != NULL) {
        double n = (double) sp->AtomCount(iel);
        return n * ((double) sc.Mu());
    } else {
        return 0.0;
    }
}

std::string FluxAnalyser::formatWhiteSpace(std::string str, unsigned int len, bool isright) {
    std::string wstr = "";
    if (str.length() <= len) {
        for (unsigned int i = 0; i < len - str.length(); i++) {
            wstr += " ";
        }
    } else {
        throw std::runtime_error("String is too long to output in flux analysis "
                            "(Mops, FluxAnalyser::formatWhiteSpace).");
    }
    if (isright) {
        wstr = str.append(wstr);
    } else {
        wstr = wstr.append(str);
    }
    return wstr;
}
