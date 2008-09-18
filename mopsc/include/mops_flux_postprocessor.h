/*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

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

#ifndef MOPS_FLUX_POSTPROCESSOR_H
#define MOPS_FLUX_POSTPROCESSOR_H

#include <vector>
#include <string>
#include "mops_timeinterval.h"
#include "mops_mechanism.h"
#include "sprog.h"
using namespace Mops;

namespace Mops
{
class FluxAnalyser
{
public:
    struct FluxPath {
        int SourceSpecies;
        int TargetSpecies;
        real Rate;
        FluxPath() {
            SourceSpecies = -1;
            TargetSpecies = -1;
            Rate = 0.0;
        }
    };
    FluxAnalyser(const Mechanism &mech, 
                 const timevector &times, 
                 const std::vector<fvector> &agpfwdrates,
                 const std::vector<fvector> &agprevrates,
                 const fvector &atemperatures);
    void addFluxElement(const Sprog::Element &elem);
    void writeFluxes(const std::string &filename);

private:
    vector<unsigned int> m_ElementIndexes;
    const Mops::Mechanism *m_mech; // The mechanism which defines reaction set.
    fvector m_times; 
    const std::vector<fvector> *m_agpfwdrates;
    const std::vector<fvector> *m_agprevrates;
    const fvector *m_atemperatures;

    std::vector<FluxPath> m_flux_network;

    // Private methods

    int getFluxPathIndex(Mops::FluxAnalyser::FluxPath &fpath);
    void addToFluxPathRate(Mops::FluxAnalyser::FluxPath &fpath);
    void calculateFluxAt(unsigned int index, unsigned int iel);
    real getTotalElementStoi(const Sprog::Kinetics::Reaction &rxn, unsigned int iel);
    real getNumberOfElementAtom(const Sprog::Stoich &sc, unsigned int iel);
    real getNumberOfElementAtom(const Sprog::Stoichf &sc, unsigned int iel);
    void sortDecendingFluxNetwork();
    void writeFluxAt(unsigned int index, unsigned int iel, std::ofstream &fout);

    string formatWhiteSpace(std::string str, unsigned int len, bool isright = true);
};
}

#endif
