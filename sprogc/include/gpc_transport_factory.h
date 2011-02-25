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


namespace Sprog{

	namespace Transport{



		class TransportFactory{

		private:

			static real Omega22Matrix[37*5];
			static real Omega11Matrix[37*5];	
			static real TStar[37];


		public:
			TransportFactory(){}
			~TransportFactory(){}

			//helper function for doing lookup and curve fit
			real polyFitOmega(real delta,real *matrixPtr) const;
			real QuadInterPol(int imin, real *oPtr, real Tr) const;
			void getTempInterval(const real rT, int &iMin, int &iMax) const;




			real getReducedTemp(const real T,  const Sprog::Species &sp) const;
			real getReducedDipole(const Sprog::Species &sp) const;
			real getOmega11( real rT, const real deltaStar) const;
			real getOmega22(const real T, const Sprog::Species &sp) const;
		};

		class PureSpeciesTransport : public TransportFactory {
		public:
			PureSpeciesTransport(){}
			~PureSpeciesTransport(){}
			//returns the viscosity in Kg/m-s
			real getViscosity(const real T, const Sprog::Species &sp) const;
			//returns the diffusion coefficeint in m^2/s
			real getSlefDiffusionCoeff(const real T, real p, const Sprog::Species &sp) const;
			//returns the thermal conductivity in J/m-s-K
			real getThermalConductivity(const real T, const real p, const real cp,
				const Sprog::Species &sp) const;
		};

		class MixtureTransport : public TransportFactory{

            static const real oneByRootEight;

		public:
			MixtureTransport(){}
			~MixtureTransport(){}
			//returns the viscosity in Kg/m-s
			real getViscosity(const real T, const Sprog::Thermo::Mixture &mix) const;
			//returns the thermal conductivity in J/m-s-K
			real getThermalConductivity(const real T, real p, const Sprog::Thermo::Mixture &mix) const;
			//returns the diffusion coefficeint in m^2/s
			real binaryDiffusionCoeff(const int j, const int k,
				const real T, const real p, const Sprog::Thermo::Mixture &mix) const;
			//returns the diffusion coefficeint in m^2/s
			std::vector<double> getMixtureDiffusionCoeff(const real T, const real p, const Sprog::Thermo::Mixture &mix) const;
			
		};


	};
}
#endif
		
