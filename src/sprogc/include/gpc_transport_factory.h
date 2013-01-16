/*
 Author(s):      Vinod Janardhanan (vj231)
 Project:        sprog (gas-phase chemical kinetics).
 Sourceforge:    http://sourceforge.net/projects/mopssuite

 Copyright (C) 2008 Vinod M Janardhanan.

 File purpose:
 This file contains the definition of a structure for a chemical species.  File
 also contains typdefs and structure definitions related to Species objects.

 Licence:
 This file is part of "sprog".

 sprog is free software; you can redistribute it and/or
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

#ifndef GPC_TRANSPORT_FACTORY_H
#define GPC_TRANSPORT_FACTORY_H

#include "gpc_species.h"
#include "gpc_mixture.h"
#include "gpc_params.h"
#include "fast_math_functions.hpp"

namespace Sprog
{
namespace Transport
{

class TransportFactory
{

private:

    static double Omega22Matrix[37 * 5];
    static double Omega11Matrix[37 * 5];
    static double TStar[37];

public:
    TransportFactory (){}
    virtual ~TransportFactory (){}

    //helper function for doing lookup and curve fit
    double polyFitOmega (double delta, double *matrixPtr) const;
    double QuadInterPol (int imin, double *oPtr, double Tr) const;
    void getTempInterval (const double rT, int &iMin, int &iMax) const;

    double getReducedTemp (const double T, const Sprog::Species &sp) const;
    double getOmega11 (double rT, const double deltaStar) const;
    double getOmega22 (const double T, const Sprog::Species &sp) const;
};

class PureSpeciesTransport
:
    public TransportFactory
{
public:
    PureSpeciesTransport (){}
    virtual ~PureSpeciesTransport (){}

    //returns the viscosity in Kg/m-s
    double getViscosity (const double T, const Sprog::Species &sp) const;

    //returns the diffusion coefficeint in m^2/s
    double getSlefDiffusionCoeff
    (
        const double T,
        double p,
        const Sprog::Species &sp
    ) const;

    //returns the thermal conductivity in J/m-s-K
    double getThermalConductivity
    (
        const double T,
        const double p,
        const double cp,
        const Sprog::Species &sp
    ) const;
};

class MixtureTransport
:
    public TransportFactory
{

    static const double oneByRootEight;
    static const double deltaStarCoeff;
    static const double binaryDiffCoeff;

public:
    MixtureTransport (){}
    virtual ~MixtureTransport (){}

    //returns the viscosity in Kg/m-s
    double getViscosity (const double T, const Sprog::Thermo::Mixture &mix) const;

    //returns the thermal conductivity in J/m-s-K
    double getThermalConductivity
    (
        const double T,
        double p,
        const Sprog::Thermo::Mixture &mix
    ) const;

    //returns the diffusion coefficeint in m^2/s
    double binaryDiffusionCoeff
    (
        const int j,
        const int k,
        const double T,
        const double p,
        const Sprog::Thermo::Mixture &mix
    ) const;

    //returns the diffusion coefficeint in m^2/s
    std::vector<double> getMixtureDiffusionCoeff
    (
        const double T,
        const double p,
        const Sprog::Thermo::Mixture &mix
    ) const;

};

} // End namespace Transport
} // End namespace Sprog

#endif
