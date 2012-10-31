/*
  Author(s):      Matthew Celnik (msc37) && Martin Martin (mm864)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This files contains structures which define reaction rate parameter sets and those for surface chemistry

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

#ifndef GPC_RATE_PARAMS_H
#define GPC_RATE_PARAMS_H
#include "gpc_params.h"
#include <string>
using namespace std; 

namespace Sprog
{
namespace Kinetics
{
    class Reaction; // Forward declare reaction class.

    // Arrhenius parameters.
    struct ARRHENIUS
    {
        double A; // Pre-exponential factor.
        double n; // Temperature exponent.
        double E; // Activation energy.

        // Constructors.
        ARRHENIUS(void) // Default constructor.
        {A=n=E=0.0;};
        ARRHENIUS(double aA, double an, double aE) // Initialising constructor.
        {A=aA; n=an; E=aE;};

        // Writes the element to a binary data stream.
        template<class Archive>
        void serialize(Archive & ar, const unsigned int /* file_version */)
        {
            ar & A & n & E;
        }
        friend class boost::serialization::access;

    };

	// Surface Coverage parameters.
    struct COVERAGE
    {
        double Eta; // Eta.
        double Miu; // Miu.
        double Epsilon; // e.
		std::string spName; // species name associated with the coverage param. 
		
        // Constructors.
        COVERAGE(void) // Default constructor.
        {Eta=Miu=Epsilon=0.0;
		spName = "";};
      COVERAGE(double e, double m, double eps, std::string &name) // Initialising constructor.
        {Eta=e; Miu=m; Epsilon=eps; spName=name;};

        // Writes the element to a binary data stream.
        template<class Archive>
        void serialize(Archive & ar, const unsigned int /* file_version */)
        {
            ar & Eta & Miu & Epsilon & spName;
        }
        friend class boost::serialization::access;
    };
	
	//  FORD parameters.
    struct FORD
    {
        double F_k; // Forward reaction order.
		std::string spName; // species name associated with the reaction order. 
		
        // Constructors.
        FORD(void) // Default constructor.
        {F_k=0.0;
		spName = "";};
      FORD(double f, std::string &name) // Initialising constructor.
        {F_k = f; spName=name;};

        // Writes the element to a binary data stream.
        template<class Archive>
        void serialize(Archive & ar, const unsigned int /* file_version */)
        {
            ar & F_k & spName;
        }
        friend class boost::serialization::access;
    };
	
	
    // Landau Teller reaction parameters.
    struct LTCOEFFS
    {
        double B, C;

        // Constructors.
        LTCOEFFS(void) {B=C=0.0;}; // Default constructor.
        LTCOEFFS(double aB, double aC) {B=aB; C=aC;}; // Initialising constructor.

        // Writes the element to a binary data stream.
        template<class Archive>
        void serialize(Archive & ar, const unsigned int /* file_version */)
        {
            ar & B & C;
        }
        friend class boost::serialization::access;

    };

    // Enumeration of possible fall-off reaction types.
    enum FALLOFF_FORM {
        None,      // No fall-off behaviour.
        Lindemann, // Lindemann form, only use low pressure limit.
        Troe3,     // 3-parameter TROE form.
        Troe4,     // 4-parameter TROE form.
        SRI,       // SRI form.
        Custom     // Custom form.  Forward + Reverse rates are adjusted by a user supplied routine.
    };

    // Fall-off parameters.
    struct FALLOFF_PARAMS
    {
        // Maximum number of fall-off parameters stored in this structure. Remember to
        // change constructor functions below if changed!
        static const int MAX_FALLOFF_PARAMS = 5;

        ARRHENIUS LowP_Limit; // Low pressure limit Arrhenius parameters.
        int ThirdBody;        // Index of the species to use as the third body.  -1 for all.
        double Params[MAX_FALLOFF_PARAMS]; // Fall-off parameters specific to whichever form is used:
                                         //   Lindemann:  No other parameters required.
                                         //   TROE3(4):  0=alpha, 1=T***, 2=T*, (3=T**).
                                         //   SRI:  0=a, 1=b, 2=c, 3=d, 4=e.
        // Constructors.
        FALLOFF_PARAMS(void) {ThirdBody=-1,Params[0]=Params[1]=Params[2]=Params[3]=Params[4]=0.0;};
        FALLOFF_PARAMS(ARRHENIUS lowp, int tb, double a=0.0, double b=0.0, double c=0.0, double d=0.0, double e=0.0) {
            LowP_Limit=lowp; ThirdBody=tb; Params[0]=a; Params[1]=b; Params[2]=c; Params[3]=d; Params[4]=e;
        };
        FALLOFF_PARAMS(double lpa, double lpn, double lpe, int tb, double a=0.0, double b=0.0, double c=0.0,
                       double d=0.0, double e=0.0) {
            LowP_Limit=ARRHENIUS(lpa,lpn,lpe); ThirdBody=tb; Params[0]=a; Params[1]=b;
            Params[2]=c; Params[3]=d; Params[4]=e;
        };

        // Writes the element to a binary data stream.
        template<class Archive>
        void serialize(Archive & ar, const unsigned int /* file_version */)
        {
            ar & LowP_Limit & ThirdBody & Params;
        }
        friend class boost::serialization::access;

    };

    // A function pointer type for a custom fall-off reaction.
    /*typedef void (*FallOffFnPtr)(
        const Reaction &rxn, // The reaction for which the fall-off is defined.
        double lowk, // The pre-calculated low-pressure limit rate constant.
        double tbc,  // The pre-calculated third-body concentration.
        double T,    // Temperature.
        double &kf,  // Forward rate constant (will be adjusted).
        double &kr   // Reverse rate constant (will be adjusted).
        );*/
};
};

#endif
