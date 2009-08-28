/* 
 * File:   interface.h
 * Author: vj231
 *
 * Created on 20 August 2009, 17:20
 */

#ifndef _INTERFACE_H
#define	_INTERFACE_H

#include <iostream>
#include <string>
#include "gpc.h"
#include "cam_read.h"
#include "cam_admin.h"
#include "cam_control.h"
#include "cam_geometry.h"
#include "cam_converter.h"
#include "cam_boundary.h"
#include "cam_profile.h"
#include "cam_configuration.h"
#include "cam_read.h"
#include "cam_models.h"
#include "cam_soot.h"
using namespace std;
using namespace Sprog;

namespace Camflow{
    class Interface{
    public:
        Interface();
        ~Interface(){}
        /*
         *calling interface to solve the flameletes. continuation true means the
         *call is a continuation call to the same problem
         */
        void flamelet(doublereal sdr, doublereal intTime=0, bool continuation=false);
        /*
         *return the vector of species names
         */
        void getSpeciesNames(vector<string>& names);
        /*
         *return the number of species
         */
        const int getNumberOfSpecies() const ;
        /*
         *return the species mass frac given the indenpendant variable
         */
        doublereal getMassFrac(int spIndex, doublereal axpos);
        /*
         *return the temperature given the independant variable
         */
        doublereal getTemperature(doublereal axpos);
        /*
         *return the density given the independant variable
         */
        doublereal getDensity(doublereal axpos);
        /*
         *return the viscosity
         */
        doublereal getViscosity(doublereal axpos);
        /*
         *return the stoichiometric mixture fraction
         */
        doublereal getStMixtureFrac();
    private:
        CamControl cc;
        CamGeometry cg;
        CamConverter convert;
        CamAdmin ca;
        CamBoundary cb;
        CamProfile cp;
        CamConfiguration config;
        CamRead cm;
        CamModels models;
        CamSoot cSoot;
        static Mechanism mech;
        int nSpecies;
        vector<string> speciesNames;
        /*
         *memberes to hold dependent variabes
         */
        Array2D spMassFracs;            //species mass fractions
        vector<doublereal> TVector;     //temperature
        vector<doublereal> rhoVector;   //density
        vector<doublereal> muVector;    //viscosity
        vector<doublereal> indVar;      //independant variable

        doublereal stMixtureFrac;       //stoichiometric mixture fraction
        /*
         *reactor pinters
         */
        FlameLet *flmlt;

        doublereal getVariableAt(const doublereal& pos, vector<doublereal>& var);

    };
}


#endif	/* _INTERFACE_H */

