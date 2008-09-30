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


#ifndef TRANSPORT_DATA_H
#define TRANSPORT_DATA_H
#include<string>
#include<map>
#include<vector>

using namespace std;
namespace Sprog{
	class Mechanism;

	namespace Transport{

		class TransportData{				
			int	   molIndex; // molecule index:- 0=atom, 1=linear molec, 2= non-linear molec
			double LJwellDepth; // L-J potential well depth e/kb (K)
			double LJcollisionDia; // L-J collistion dia
			double dipol; // Dipole moment Debye
			double polarity; // Polarizability
			double rotRelaxNum; // Rotational relaxation number

		public:
			TransportData();
			TransportData(int molIndex,
				double LJwellDepth,
				double LJcollisionDia,
				double dipol,
				double polarity,
				double rotRelaxNum);
			//TransportData(const TransportData &td); // Copy constructor

			
			~TransportData(){};

			TransportData operator=(TransportData td);// overloaded operator

			void setMolIndex(int molIndex);
			int getMolIndex() const;

			void setWellDepth(double LJwellDepth);
			double getWellDepth() const;
		
			void setCollisionDia(double LJcollisionDia);
			double getCollisionDia() const;

			void setDipole(double dipol);
			double getDipole()const;

			void setPolarity(double polarity);
			double getPolarity()const;

			void setRotRelaxNum(double rotRelaxNum);
			double getRotRelaxNum()const;

			void validateTransport(map<string,vector<string>> &trMap,Sprog::Mechanism &mech);
		};
	};
};

#endif


		


