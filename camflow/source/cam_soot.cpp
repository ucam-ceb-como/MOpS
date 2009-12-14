/*
 * File:   cam_soot.cpp
 * Author: vj231
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class implements the soot moment (MOMIC) based on
 *  chemkin implementation
 * Licence:
 *  This file is part of "Camflow".
 *
 *  Camflow is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * Contact:
 *  Dr Markus Kraft
 *  Dept of Chemical Engineering
 *  University of Cambridge
 *  New musseum Site
 *  Pembroke Street
 *  Cambridge
 *  CB2 3RA
 *  UK
 *
 *  Email   :   mk306@cam.ac.uk
 *  Website :   http://como.cheng.cam.ac.uk
 *
 * Created on 19 May 2009, 10:38
 */
#include <algorithm>
#include <vector>
#include <cmath>
#include <map>

#include "array.h"
#include "cam_soot.h"
#include "cam_reporter.h"
#include "gpc.h"
#include "cam_math.h"
#include "cam_setup.h"
#include "comostrings.h"
using namespace Camflow;
using namespace Sprog;
using namespace Strings;
Array1D CamSoot::sizeMoments;
Array1D CamSoot::reducedMoments;

CamSoot::CamSoot():cMass(0.012), //kg/mol
        rhoSoot(1.8e+03),        //kg/m^3
        ohMass(0.017),           //kg/mol
        lambda(2.3e15),          //cm^-2
        atomsPerDimer(32),
        D1_PAH(2.42e-10)         //m
{ 

    momentON = false;    
    lowFrac = -4;
    
}

void CamSoot::setNumMoments(int n){
    mmn2Calc = n;
    nMoments = 6;
    highFrac = 6*nMoments+1; //highet fractional moment
}
void CamSoot::setNucCutoff(doublereal cutOff){
    nucRateCutOff = cutOff;
}

void CamSoot::setNumCAtomInception(int nCAtom){
    numCAtomInception = nCAtom;

}

void CamSoot::setPAHDia(doublereal dPAH){
    dia_PAH = dPAH;

}

void CamSoot::setSootSpecies(std::vector<std::string> species){
    sootSpecies = species;

}

void CamSoot::setSootMomentActive(){
    momentON = true;

}

void CamSoot::setRegime(int n){
    tRegime = n;

}
/*
 *set the user defined first moment
 */
void CamSoot::setFirstMom(doublereal m0){
    firstMom = m0;
}
int CamSoot::getRegime(){
    return tRegime;
}
bool CamSoot::active(){
    return momentON;
}

int CamSoot::getNumMoments() const{
    return mmn2Calc;
}

void CamSoot::setInceptionSpecies(std::string species){
    iSpecies = species;
}

void CamSoot::initialize(int nCells, Mechanism &mech,
                                        std::vector<doublereal> &mSolnVec){


    /*
     *reduced mass for nucleation for calculating
     *constant part of nucleation. Needs to be sizeMomentsltiplied by sqrt(T)
     *before using (final unit = m^3/(s mol^2). mol^2 is required since
     *this is sizeMomentsltiplied by the concentraction of PAHs
     */
    doublereal rMass = numCAtomInception*cMass/NA;
    constNucleation = 2.2*sqrt(4*pi*kB_cgs/rMass)*dia_PAH*1.0e2*dia_PAH*1.0e2*NA*NA;


    /*
     *Beta for free molecular coagulation, needs to be multiplied by
     *sqrt(T)
     */
    Kf = 2.2*pow(((3*cMass/NA)/(4*pi*rhoSoot)),1.0/6.0)*sqrt(6*kB_cgs/rhoSoot);


    /*
     *constant for continuum coagulation, needs to be multiplied by
     *free path
     */
    doublereal mass = cMass/NA;
    doublereal CD1 = pow((6.0*mass/pi/rhoSoot),1.0/3.0); //unit [=] cm


    Kc_ = 2.514/CD1;


    /*
     *soot hydroxyl oxidation constant [=] cm^3/mol-s
     */
    doublereal oh = ohMass/NA;
    kOH = CD1*CD1* sqrt(pi*kB_cgs/2.0/oh)*NA;


    /*
     *constant for surface growth [=] --
     */
    kSurf = CD1*CD1*lambda*pi;


    /*
     *PAH condensation constant [=] cm/s-mol
     */
    kPAH = 2.2*sqrt(pi*kB_cgs/2.0/mass)*NA;



    doublereal CD1_PAH = D1_PAH*sqrt(2.0/3.0); //[=] cm
    doublereal CH2 = kPAH*CD1_PAH*CD1_PAH;           //[=]cm^3/mol-s
    doublereal CHS = kPAH*2*CD1_PAH*CD1;             //[=]cm^3/mol-s
    doublereal CS2 = kPAH*CD1*CD1;                   //[=]cm^3/mol-s


    powPAH.resize(4,nMoments);


    for(int i=1; i<nMoments; i++){
        doublereal x = (2.0*i+1)/2.0;
        powPAH(1,i) = CH2 * pow(doublereal(numCAtomInception),x);
        powPAH(2,i) = CHS * pow(doublereal(numCAtomInception),i);
        x  = (2.0*i-1)/2.0;
        powPAH(3,i) = CS2 * pow(numCAtomInception,x);

    }

    /*
     *populate the solution vector
     */

    mSolnVec.resize(nMoments,firstMom);
    for(int i=1; i<nMoments; i++){
        mSolnVec[i] = mSolnVec[i-1] + log(doublereal(atomsPerDimer));
    }
    std::vector<doublereal> temp = mSolnVec;
    mSolnVec.resize(nMoments*nCells,0.0);
    for(int i=0; i<nCells; i++){
        for(int l=0; l<nMoments; l++){
            mSolnVec[i*nMoments+l] = temp[l];
        }
    }


    CamMath cm;
    cm.binomCoeff(nMoments,bnCoeff);
    cm.prime(nMoments,prime);

    /*
     *species index initialization
     */
    iInception = mech.FindSpecies(convertToCaps(trim(iSpecies)));
    iC2H2 = mech.FindSpecies("C2H2");
    iCO = mech.FindSpecies("CO");
    iH = mech.FindSpecies("H");
    iH2 = mech.FindSpecies("H2");
    iH2O = mech.FindSpecies("H2O");
    iO2 = mech.FindSpecies("O2");
    iOH = mech.FindSpecies("OH");
    //rates
    wdot.resize(nCells,nMoments);
    surfProdRate.resize(nCells,mech.SpeciesCount());
    conc_received.resize(nCells,mech.SpeciesCount());
    momReceived.resize(nCells,nMoments);

}


/*
 *soot reactions
 */
void CamSoot::sootReactions(int cell, std::vector<doublereal>& conc, std::vector<doublereal>& mom,
                                    int nSpec, doublereal T, doublereal p){

    /*
     *convert all concentration to mol/cm^3
     */
    std::vector<doublereal> conc_cgs;
    conc_cgs.resize(conc.size(),0.0);
    for(unsigned int i=0; i<conc.size(); i++)
        conc_cgs[i] = conc[i]*1e-06;

    /*
     *store the concentrations for output
     */
    conc_received(cell,iC2H2) = conc_cgs[iC2H2];
    conc_received(cell,iCO) = conc_cgs[iCO];
    conc_received(cell,iH) = conc_cgs[iH];
    conc_received(cell,iH2) = conc_cgs[iH2];
    conc_received(cell,iH2O) = conc_cgs[iH2O];
    conc_received(cell,iInception) = conc_cgs[iInception];
    conc_received(cell,iO2) = conc_cgs[iO2];
    conc_received(cell,iOH) = conc_cgs[iOH];

    /*
     *store the moment concentrations received
     */
    for(int i=0; i<nMoments; i++)
        momReceived(cell,i) = mom[i];

   
    /*
     *nucleation rats
     */
    nucleation(conc_cgs,T,p);
    /*
     *update size moments
     */
    std::vector<doublereal> fp;
    fp.resize(nMoments,1.0);
    for(int i=1; i<nMoments; i++){
        fp[i] = mom[i]/mom[0];;
    }
    setSizeMoments(fp);
    /*
     *coagulation
     */
    coagulation(T,p/1e5,mom[0]);

    /*
     *condensation
     */
    doublereal ratePAH;
    condensation(T,mom[0],conc_cgs[iInception], ratePAH);


    /*
     *surface reactions
     */
    doublereal a,b,rsq;
    std::vector<doublereal> surfRates;
    surface(nMoments,cell,T,mom[0],conc_cgs,surfRates);
    surfRates[0] = 0.0;
    if(surfRates[1] < 0.0 ){
        for(int i=1; i<nMoments; i++){
            fp[i]=log(sizeMoments(6*i));

        }

        linear(nMoments,fp,a,b,rsq);

        for(int i=4; i<(6*nMoments-2);i+=6){
            sizeMoments(i) = exp(a+b*i/6.0);

        }

        surface(nMoments,cell,T,mom[0],conc_cgs,surfRates);

        for(int i=1; i<nMoments; i++){
            surfRates[i] *= (mom[i]/exp(a+b*i)/mom[0]);
        }
        surfRates[0]=0.0;
    }
    for(int i=0; i<nMoments;i++){
        surfRates[i] += cdRate[i];
    }

    /*
     *rates of moments
     */

    for(int i=0; i<nMoments; i++){
        wdot(cell,i) = (nucRate[i] + cgRate[i] + surfRates[i])/mom[i];
       
    }



}
/*
 *calculate the nucleation rate
 */
void CamSoot::nucleation(std::vector<doublereal>& conc, doublereal T, doublereal p){

    nucRate.resize(nMoments,0);

    /*
     *Nucleation rate for zeroth moment
     */
    doublereal concPAH = std::max(0.0,conc[iInception]);
    nucRate[0] = constNucleation*sqrt(T)*concPAH*concPAH;

    /*
     *nucleation rate for higher order moments
     */
    for(int i=1; i<nMoments; i++)
        nucRate[i] = 2*numCAtomInception*nucRate[i-1];


}

/*
 *coagulation
 */
void CamSoot::coagulation(doublereal T, doublereal p, doublereal M0){


    doublereal kCoag = Kf*sqrt(T);


    std::vector<doublereal> f;
    CamMath cm;
    f.clear();
    for(int i=0; i<4; i++){
        f.push_back(grid(i,0,0));//    int dd; cin >> dd;
        
    }

    doublereal crk = kCoag*cm.interpolateLG(0.5,4,prime,f);
    
    
    f.clear();
    for(int i=0; i<4; i++)
        f.push_back(grid(i,1,1));
    doublereal crk2 = kCoag*cm.interpolateLG(0.5,4,prime,f);


    f.clear();
    for(int i=0; i<4; i++)
        f.push_back(grid(i,1,2));
    doublereal crk3 = kCoag*cm.interpolateLG(0.5,4,prime,f);


    f.clear();
    for(int i=0; i<3; i++)
        f.push_back(grid(i,1,3));
    doublereal crk4a = kCoag*cm.interpolateLG(0.5,3,prime,f);


    f.clear();
    for(int i=0; i<4; i++)
        f.push_back(grid(i,2,2));
    doublereal crk4b = kCoag*cm.interpolateLG(0.5,4,prime,f);


    f.clear();
    for(int i=0; i<2; i++)
        f.push_back(grid(i,1,4));
    doublereal crk5a = kCoag*cm.interpolateLG(0.5,2,prime,f);


    f.clear();
    for(int i=0; i<3; i++)
        f.push_back(grid(i,2,3));
    doublereal crk5b = kCoag*cm.interpolateLG(0.5,3,prime,f);


    if(tRegime != FM){
        /*
         *viscosity for air
         */
        doublereal viscosity = 14.58e-06*pow(T,1.5)/(T+100.4);

        /*
         *mean free path
         */
        doublereal freePath = 2.37e-03*T/p/1.103e05;


        doublereal ckt = (2*kB_cgs/3.0)*T/viscosity;
        doublereal cat = Kc_*freePath;


        doublereal crkc     = ckt*( betaC1(0,0) + cat*betaC2(0,0) );
        doublereal crk2c    = ckt*( betaC1(1,1) + cat*betaC2(1,1) );
        doublereal crk3c    = ckt*( betaC1(1,2) + cat*betaC2(1,2) );
        doublereal crk4Ac   = ckt*( betaC1(1,3) + cat*betaC2(1,3) );
        doublereal crk4Bc   = ckt*( betaC1(2,2) + cat*betaC2(2,2) );
        doublereal crk5Ac   = ckt*( betaC1(1,4) + cat*betaC2(1,4) );
        doublereal crk5Bc   = ckt*( betaC1(2,3) + cat*betaC2(2,3) );

        /*
         *Harmonic mean
         */
        crk  *= crkc/(crk + crkc);
        crk2 *= crk2c /(crk2 + crk2c);
        crk3 *= crk3c /(crk3 + crk3c);
        crk4a *= crk4Ac/(crk4a + crk4Ac);
        crk4b *= crk4Bc/(crk4b + crk4Bc);
        crk5a *= crk5Ac/(crk5a + crk5Ac);
        crk5b *= crk5Bc/(crk5b + crk5Bc);
    }



    doublereal m02 = M0*M0;

    cgRate.resize(6,0.0);

    cgRate[0] = -0.5*crk * m02;
    cgRate[1] = 0;
    cgRate[2] = crk2 * m02;
    cgRate[3] = 3*crk3 * m02;
    cgRate[4] = (4*crk4a + 3*crk4b) * m02;
    cgRate[5] = (5* crk5a + 10*crk5b) * m02;



}


/*
 *condensation/
 */
void CamSoot::condensation( doublereal T,
                            doublereal M0,
                            doublereal concPAH,
                            doublereal ratePAH){
    cdRate.resize(nMoments,0.0);

    doublereal coeff = sqrt(T)*concPAH*M0;

    for (int r = 1; r < nMoments; r++) {
        for (int l = 0; l < r; l++) {
            int l6 = 6*l;
            cdRate[r] += bnCoeff(l,r)*(
                    powPAH(1,r-l)*sizeMoments(l6) +
                    powPAH(2,r-l)*sizeMoments(l6+2) +
                    powPAH(3,r-l)*sizeMoments(l6+4)
                    );
        }

        cdRate[r] *= coeff;

    }



    //PAH consumption due to condensation
    ratePAH = -cdRate[1]/numCAtomInception/NA;

}

/*
 *surface reaction rates
 */
void CamSoot::surface(int hMoment, int cell, doublereal T, doublereal M0,
                        std::vector<doublereal>& conc,
                        std::vector<doublereal>& totalRates){


    doublereal RT = 1.987e-3*T;
    doublereal fr1 = 4.2e13*exp(-13.0/RT)*conc[iH];
    doublereal rr1 = 3.9e12*exp(-11.0/RT)*conc[iH2];
    doublereal fr2 = 1.0e10 * pow(T,0.734) * exp(-1.43/RT)*conc[iOH];
    doublereal rr2 = 3.68e08 * pow(T,1.139) * exp(-17.1/RT)*conc[iH2O];
    doublereal fr3 = 2e13 * conc[iH];
    doublereal fr4 = 8.0e7*pow(T,1.56)*exp(-3.8/RT)*conc[iC2H2];
    doublereal fr5 = 2.2e12*exp(-7.5/RT)*conc[iO2];
    doublereal fr6 = 0.13*conc[iOH];

    doublereal par_a = 12.65 - 5.63e-03*T;
    doublereal par_b = -1.38 + 6.80e-04*T;

    doublereal alpha = tanh(par_a/log10(sizeMoments(6)) + par_b );

    
    doublereal denom = rr1 + rr2 + fr3 + fr4 + fr5;
    smRates.clear();
    doublereal p_rate;
    std::vector<doublereal> rateC2H2, rateO2, rateOH;
    if(denom != 0.){
        doublereal ssRatio = (fr1+fr2)/denom;
        doublereal cArea = alpha * kSurf * M0;
        doublereal cRad = cArea * ssRatio;
        /*
         *C2H2
         */
        
        
        rateC2H2.resize(hMoment,0.0);
        doublereal coef = fr4*cRad;
        sums(hMoment,2,coef,rateC2H2);
        smRates.insert(make_pair("C2H2",rateC2H2));

        /*
         *O2
         */
        rateO2.resize(hMoment,0.0);
        coef = fr5*cRad;
        sums(hMoment,-2,coef,rateO2);
        smRates.insert(make_pair("O2",rateO2));

        /*
         *OH the unit in 1/s
         */
        rateOH.resize(hMoment,0.0);
        coef = fr6*kOH*sqrt(T)*M0;
        sums(hMoment,-1,coef,rateOH);
        smRates.insert(make_pair("OH",rateOH));


        /*
         *production rates C2H2
         */
        p_rate = -rateC2H2[1]/2/NA;
        //surfProdRate.insert(make_pair("C2H2",p_rate));
        surfProdRate(cell,iC2H2) = p_rate;

        /*
         *production rate O2
         */
        p_rate = rateO2[1]/2/NA;
        //surfProdRate.insert(make_pair("O2",p_rate));
        surfProdRate(cell,iO2) = p_rate;


        /*
         *production rate OH
         */
        p_rate = ((rr2*cRad - fr2*cArea)*sizeMoments(4) + rateOH[1])/NA;
        //surfProdRate.insert(make_pair("OH",p_rate));
        surfProdRate(cell,iOH) = p_rate;


        /*
         *production rate H
         */
        p_rate =  ((rr1+fr4)*cRad - fr1*cArea - fr3*cRad)*sizeMoments(4)/NA;
        //surfProdRate.insert(make_pair("H",p_rate));
        surfProdRate(cell,iH) = p_rate;

        /*
         *production rate H2
         */
        p_rate = (fr1*cArea - rr1*cRad)*sizeMoments(4)/NA;
        //surfProdRate.insert(make_pair("H2",p_rate));
        surfProdRate(cell,iH2) = p_rate;

        /*
         *production rate H2O
         */
        p_rate = (fr2*cArea - rr2*cRad)*sizeMoments(4)/NA;
        //surfProdRate.insert(make_pair("H2O",p_rate));
        surfProdRate(cell,iH2O) = p_rate;

        /*
         *production rate CO
         */
        p_rate = -(rateO2[1]+rateOH[1])/NA;
        //surfProdRate.insert(make_pair("CO",p_rate));
        surfProdRate(cell,iCO) = p_rate;



    }else{
        rateC2H2.resize(hMoment,0.0);
        rateO2.resize(hMoment,0.0);
        rateOH.resize(hMoment,0.0);
        smRates.insert(make_pair("C2H2",rateC2H2));
        smRates.insert(make_pair("O2",rateC2H2));
        smRates.insert(make_pair("OH",rateC2H2));

        p_rate = 0.0;
//        surfProdRate.insert(make_pair("C2H2",p_rate));
//        surfProdRate.insert(make_pair("CO",p_rate));
//        surfProdRate.insert(make_pair("O2",p_rate));
//        surfProdRate.insert(make_pair("OH",p_rate));
//        surfProdRate.insert(make_pair("H",p_rate));
//        surfProdRate.insert(make_pair("H2",p_rate));
        surfProdRate(cell,iC2H2) = p_rate;
        surfProdRate(cell,iCO) = p_rate;
        surfProdRate(cell,iH) = p_rate;
        surfProdRate(cell,iH2) = p_rate;
        surfProdRate(cell,iH2O) = p_rate;
        surfProdRate(cell,iO2) = p_rate;
        surfProdRate(cell,iOH) = p_rate;


    }

    totalRates.resize(hMoment,0.0);
    for (int i = 0; i < hMoment; i++) {
        totalRates[i] = rateC2H2[i] + rateO2[i] + rateOH[i];

    }



}

/*
 *calculate the sum of surface growth term
 */
void CamSoot::sums(int hMoment, doublereal massAdded, doublereal coeff,
                        std::vector<doublereal>& rates){
    for (int r = 1; r < hMoment; r++) {
        rates[r] = 0.0;
        for (int l = 0; l < r; l++) {
            rates[r] += bnCoeff(l,r)*pow(massAdded,r-l)*reducedMoments(6*l+4);            
        }        
        rates[r] *= coeff;
        
    }

}
/*
 *grid function
 */
doublereal CamSoot::grid(int k, int n, int m){

    doublereal gfun = 0;
    for(int l=0; l<=k; l++){
        int i = 6*(k-l+n);
        int j = 6*(l+m);
        gfun += bnCoeff(l,k)*(
                sizeMoments(i+1)*sizeMoments(j-3)+
                2*sizeMoments(i-1)*sizeMoments(j-1)+
                sizeMoments(i-3)*sizeMoments(j+1)
                );
    }
    return gfun;
}

doublereal CamSoot::gridFunction(int k, int n, int m){

    doublereal gfun = 0;
    for(int l=0; l<=k; l++){
        int i = 6*(k-l+n);
        int j = 6*(l+m);
        gfun += bnCoeff(l,k)*(
                reducedMoments(i+1)*reducedMoments(j-3)+
                2*reducedMoments(i-1)*reducedMoments(j-1)+
                reducedMoments(i-3)*reducedMoments(j+1)
                );
    }
    return gfun;
}


doublereal CamSoot::betaC1(int i, int j){
    int i6,j6;
    i6 = 6*i;
    j6 = 6*j;

    doublereal x = sizeMoments(-2+i6)*sizeMoments(2+j6) +
                2*sizeMoments(i6)*sizeMoments(j6) +
            sizeMoments(2+i6)*sizeMoments(-2+j6);

    return x;
}

doublereal CamSoot::betaC2(int i, int j){

    int i6, j6;
    i6 = 6*i;
    j6 = 6*j;

    doublereal x = sizeMoments(-4+i6)*sizeMoments(2+j6)+
            sizeMoments(-2+i6)*sizeMoments(j6)+
            sizeMoments(i6)*sizeMoments(-2+j6)+
            sizeMoments(2+i6)*sizeMoments(-4+j6);

    return x;
}

void CamSoot::linear(int n, std::vector<doublereal>& y, doublereal& a, doublereal& b,
                        doublereal& rsq){
    doublereal sum = 1.0*n;
    doublereal x1 = 0.0;
    doublereal x2 = 0.0;
    doublereal y1 = 0.0;
    doublereal y2 = 0.0;

    doublereal yMean, ess, yss,yPred;

    for(int i=1; i<n;i++){
        x1 += doublereal(i);
        x2 += doublereal(i)*doublereal(i);

        y1 += y[i];
        y2 += y[i]*doublereal(i);

    }
    doublereal d = sum*x2 - x1*x1;
    a = (y1*x2 - x1*y2)/d;
    b = (sum*y2 - y1*x1)/d;
    yMean = y1/n;
    ess = 0.0;
    yss = 0.0;

    for(int i=1; i<n; i++){
        yPred = a+b*i;
        ess += (yMean-yPred)*(yMean-yPred);
        yss += (y[i]-yMean)*(y[i]-yMean);
    }

    rsq = ess/yss;

}

void CamSoot::setSizeMoments(std::vector<doublereal>& fp){
    //std::vector<doublereal> fp;
    //fp.resize(nMoments,1.0);
    CamMath cm;

    //for(int i=1; i<nMoments; i++){
    //    fp[i] = conc[i]/conc[0];;
    //}

    sizeMoments.resize(lowFrac,highFrac+1);

    for(int i=lowFrac; i<=-1; i++){
        sizeMoments(i) = cm.interpolateLG(i/6.0,3,prime,fp);

    }
//
//
    for(int i=0; i<=(highFrac-nMoments); i++){
        if((i%6) != 0.0){
            sizeMoments(i) = cm.interpolateLG(i/6.0,nMoments,prime,fp);

        }else{
            sizeMoments(i) = fp[i/6];

        }
    }

    for(int i=(highFrac-nMoments)+1; i<=highFrac; i++){
        sizeMoments(i) = 0;

    }
}


void CamSoot::clearRates(int nCells){
    surfProdRate.resize(nCells,nMoments);
}

void CamSoot::report(int nCells){
    std::vector<std::string> header;
    header.push_back("C2H2");
    header.push_back("CO");
    header.push_back("H");
    header.push_back("H2");
    header.push_back("H2O");
    header.push_back("incept");
    header.push_back("O2");
    header.push_back("OH");

    header.push_back("w0");
    header.push_back("w1");
    header.push_back("w2");
    header.push_back("w3");
    header.push_back("w4");
    header.push_back("w5");


    header.push_back("M0");
    header.push_back("M1");
    header.push_back("M2");
    header.push_back("M3");
    header.push_back("M4");
    header.push_back("M5");

    CamReporter cr;
    cr.openFile("moments.dat",false);
    cr.writeCustomHeader(header);

    for(int i=0; i<nCells; i++){
        std::vector<doublereal> data;
        data.push_back(conc_received(i,iC2H2));
        data.push_back(conc_received(i,iCO));
        data.push_back(conc_received(i,iH));
        data.push_back(conc_received(i,iH2));
        data.push_back(conc_received(i,iH2O));
        data.push_back(conc_received(i,iInception));
        data.push_back(conc_received(i,iO2));
        data.push_back(conc_received(i,iOH));

        for(int k=0; k<nMoments; k++)
            data.push_back(wdot(i,k));

        for(int k=0; k<nMoments; k++)
            data.push_back((momReceived(i,k)));

        cr.writeCustomFileOut(data);

    }

    cr.closeFile();

}

/*
 *Initialize the soot moments based on user defined
 *zeroth moments.
 */
void CamSoot::initMoments(Mechanism &mech, std::vector<doublereal>& soln,int nCells){
   
    int st = soln.size();
    /*
     *zeroth moment based on user input
     */
    soln.push_back(firstMom);
    /*
     *rest of the moments
     */
    for(int i=1; i<nMoments; i++){
       
        soln.push_back(soln[i-1+st]+ log(doublereal(atomsPerDimer)));

    }

    /*
     *preperation for interpolation
     */

    CamMath cm;
    cm.prime(nMoments,prime);
    cm.binomCoeff(nMoments,bnCoeff);

    if(nCells>0){
        wdot.resize(nCells,nMoments);
        surfProdRate.resize(nCells,mech.SpeciesCount());
    }

    iInception = mech.FindSpecies(convertToCaps(trim(iSpecies)));
    iC2H2 = mech.FindSpecies("C2H2");
    iCO = mech.FindSpecies("CO");
    iH = mech.FindSpecies("H");
    iH2 = mech.FindSpecies("H2");
    iH2O = mech.FindSpecies("H2O");
    iO2 = mech.FindSpecies("O2");
    iOH = mech.FindSpecies("OH");


    /*--------------------------------------------------------------------------
     *constant part of nucleation rate kernel
     *-------------------------------------------------------------------------/
     */

     //molecular weight if PAH
    doublereal mwPAH = numCAtomInception*cMass;

     //mass of PAH
    doublereal mPAH = mwPAH/NA;
    /*
     *      nucleation kernel  Balthase thesis
     *               |-----------
     *         e_ij *| 8*pi*kB*T
     *               |-----------  (r_i+r_j)^2
     *               | mu_ij
     *              \|
     */
    Beta_nucl = 2.2*sqrt(16*pi*kB/mPAH)*dia_PAH*dia_PAH;
    /*
     * Multiplication by NA is carried out to avoid it while
     * evaluating the rates
     */
    Beta_nucl *= NA*NA;


    /*--------------------------------------------------------------------------
     *                constant for free molecular coagulation
     *-------------------------------------------------------------------------/
     */
    Beta_fm = 2.2*sqrt(6*kB/rhoSoot)* pow(3*(cMass/NA)/(4*pi*rhoSoot),(1.0/6.0));
    
    /*--------------------------------------------------------------------------
     *               constant fro condensation
     *-------------------------------------------------------------------------/
     */
    Beta_cd = 2.2*sqrt(pi*kB/(2*cMass/NA))*NA;
    //diameter of one particle
    doublereal CD1 = pow((6*(cMass/NA)/(pi*rhoSoot)),1.0/3.0);
    doublereal CD1_PAH = D1_PAH*sqrt(2.0/3.0);
    doublereal CH2 = Beta_cd*CD1_PAH*CD1_PAH;
    doublereal CHS = Beta_cd*CD1_PAH*CD1;
    doublereal CS2 = Beta_cd*CD1*CD1;
    powPAH.resize(5,nMoments);
    for(int i=1; i<nMoments; i++){
        doublereal x = (2.0*i+1)/2.0;
        powPAH(1,i) = CH2 * pow(doublereal(numCAtomInception),x);
        powPAH(2,i) = CHS * pow(doublereal(numCAtomInception),i);
        x  = (2.0*i-1)/2.0;
        powPAH(3,i) = CS2 * pow(numCAtomInception,x);

    }
    /*--------------------------------------------------------------------------
     *                  constants for surface growth
     *-------------------------------------------------------------------------/
     */
    Beta_surf = lambda*1.0e+04*pi*CD1*CD1;

}
/*
 *residual evaluation
 */
void CamSoot::residual(const doublereal& time, std::vector<doublereal>& wdot,
                        doublereal* y, doublereal* f){
    for(int r=0; r<nMoments; r++){
        f[r] = wdot[r];        
    }
    
}

/*
 *calculate all rates
 */
void CamSoot::rateAll(std::vector<doublereal>& conc,     //species concentration
                      std::vector<doublereal>& moments,  //moments
                      doublereal& T,                //temperature
                      doublereal& p,                //pressure
                      std::vector<doublereal>& rates,    //return rates
                      int cellID){  
    
        
    std::vector<doublereal> nucRates, coagRates,cdRates, prodRates, sRates;
    //Calculate nucleation rate
    rateNucleation(conc[iInception],T,nucRates);
    
    //Calculate coagulation rates
    rateCoagulation(moments,T,coagRates);

//    //condensation rates
//    doublereal rA4 = rateCondensation(moments,T,conc[iInception],cdRates);
//    surfProdRate(cellID,iInception) = rA4;

    //surface rates : Always calculate this after coagulation
    rateSurface(conc,T,moments,prodRates,sRates);
    sRates[0] = 0.0;

    if(sRates[1] <= 0){
        std::vector<doublereal> f;
        f.resize(nMoments,0.0);
        for(int r = 1; r< nMoments; r++){
            f[r] = log(reducedMoments(6*r));
        }


        doublereal a, b, rsq;
        linear(nMoments,f,a,b,rsq);

        for(int n = 4; n<(6*nMoments-2); n += 6){
            reducedMoments(n) = exp(a+b*n/6.0);
        }

        rateSurface(conc,T,moments,prodRates,sRates);
        for(int n=1; n<nMoments; n++){
            sRates[n] *= (moments[n]/(exp(a+b*n)*moments[0]));
        }
        sRates[0] = 0.0;
    }
    

//    for(int i=0; i<nMoments; i++){
//        sRates[i] += cdRates[i];
//    }

    //store the surface production rate for use later
    surfProdRate(cellID,iC2H2) = prodRates[iC2H2];
    surfProdRate(cellID,iO2) = prodRates[iO2];
    surfProdRate(cellID,iH) = prodRates[iH];
    surfProdRate(cellID,iH2) = prodRates[iH2];
    surfProdRate(cellID,iH2O) = prodRates[iH2O];

    for(int m=0; m<nMoments;m++){
        rates[m] = (nucRates[m]+coagRates[m])/moments[m];
        //std::cout << nucRates[m] << "  " << coagRates[m] << std::endl;
        wdot(cellID,m) = rates[m];
        
    }
    

}

/*
 *nucleation rate
 */
void CamSoot::rateNucleation(doublereal& concPAH, doublereal& T,
                                std::vector<doublereal>& nucRates){

    doublereal kNucl = Beta_nucl*sqrt(T);
    nucRates.resize(nMoments,0.0);
    /*
     *nucleation rate for the zeroth moment
     */
    nucRates[0] = kNucl* std::max(concPAH,0.0)*concPAH;
    
    /*
     *number of carbon atoms per dimer
     */
    doublereal cDimer = 2*numCAtomInception;
    /*
     *nucleation rate for higher moments
     */
    for(int m=1; m<nMoments; m++)
        nucRates[m] = nucRates[m-1] * cDimer;

}
/*
 *coagulation rate
 */
void CamSoot::rateCoagulation(std::vector<doublereal>& mom, doublereal& T,
                            std::vector<doublereal>& coagRates){

    doublereal kCoag = Beta_fm * sqrt(T);
    CamMath cm;
    /*
     *since the lagrangian interpolations are
     *carried out on the log of reduced moments
     *store the exponents of the moments
     */
    std::vector<doublereal> wholeOrderRedMom; //whole order reduced moments
    wholeOrderRedMom.resize(mom.size(),0.0);
    for(int r=0; r<nMoments; r++){
        wholeOrderRedMom[r] = mom[r]/mom[0];

    }
    interpolateReducedMoments(wholeOrderRedMom);
//    /*
//     *evaluate the grid functions
//     */
    std::vector<doublereal> f;
    for(int i=0; i<4; i++){
        f.push_back(gridFunction(i,0,0));
        
    }
            
    doublereal crk = kCoag*cm.interpolateLG(0.5,4,prime,f);    
    

    f.clear();
    for(int i=0; i<4; i++)
        f.push_back(gridFunction(i,1,1));
    doublereal crk2 = kCoag*cm.interpolateLG(0.5,4,prime,f);


    f.clear();
    for(int i=0; i<4; i++)
        f.push_back(gridFunction(i,1,2));
    doublereal crk3 = kCoag*cm.interpolateLG(0.5,4,prime,f);


    f.clear();
    for(int i=0; i<3; i++)
        f.push_back(gridFunction(i,1,3));
    doublereal crk4a = kCoag*cm.interpolateLG(0.5,3,prime,f);


    f.clear();
    for(int i=0; i<4; i++)
        f.push_back(gridFunction(i,2,2));
    doublereal crk4b = kCoag*cm.interpolateLG(0.5,4,prime,f);


    f.clear();
    for(int i=0; i<2; i++)
        f.push_back(gridFunction(i,1,4));
    doublereal crk5a = kCoag*cm.interpolateLG(0.5,2,prime,f);


    f.clear();
    for(int i=0; i<3; i++)
        f.push_back(gridFunction(i,2,3));
    doublereal crk5b = kCoag*cm.interpolateLG(0.5,3,prime,f);

    doublereal M02 = mom[0]*mom[0];

    coagRates.resize(nMoments,0.0);
    coagRates[0] = -0.5*crk*M02;
    coagRates[1] = 0.0;
    coagRates[2] = crk2 * M02;
    coagRates[3] = 3*crk3 * M02;
    coagRates[4] = (4*crk4a + 3*crk4b) * M02;
    coagRates[5] = (5* crk5a + 10*crk5b) * M02;

    
}

/*
 *condensation rates
 */
doublereal CamSoot::rateCondensation(std::vector<doublereal>& mom,
                                    doublereal& T,
                                    doublereal& conc,
                                    std::vector<doublereal>& cdRates){

    doublereal k = sqrt(T)*conc*mom[0];
    cdRates.resize(nMoments,0.0);
    for(int r=1; r<nMoments; r++){
        for(int l=0; l<r;l++){
            int l6 = 6*l;
            cdRates[r] += bnCoeff(l,r)*(
                    powPAH(1,r-l)*reducedMoments(l6)+
                    powPAH(2,r-l)*reducedMoments(l6+2)+
                    powPAH(3,r-l)*reducedMoments(l6+4));
        }
        cdRates[r] *= k;
    }
    //return the PAH consumption rate due to condensation
    doublereal prodRate = -cdRates[1]/numCAtomInception/NA;
    return prodRate;

}

/*
 *surface reaction rates
 */
void CamSoot::rateSurface(std::vector<doublereal>& conc,
                            doublereal T,
                            std::vector<doublereal>& mom,
                            std::vector<doublereal>& prodRates,
                            std::vector<doublereal>& sRates){

 
    
    doublereal RT = 1.987e-3*T;
    doublereal fr1 = 4.2e13*exp(-13.0/RT)*conc[iH];
    doublereal rr1 = 3.9e12*exp(-11.0/RT)*conc[iH2];
    doublereal fr2 = 1.0e10 * pow(T,0.734) * exp(-1.43/RT)*conc[iOH];
    doublereal rr2 = 3.68e08 * pow(T,1.139) * exp(-17.1/RT)*conc[iH2O];
    doublereal fr3 = 2e13 * conc[iH];
    doublereal fr4 = 8.0e13*pow(T,1.56)*exp(-3.8/RT)*conc[iC2H2];
    doublereal fr5 = 2.2e12*exp(-7.5/RT)*conc[iO2];
    //doublereal fr6 = 0.13*conc[iOH];

    


    doublereal par_a = 12.65 - 5.63e-03*T;
    doublereal par_b = -1.38 + 6.80e-04*T;

    doublereal alpha = tanh(par_a/log10(reducedMoments(6)) + par_b );

    doublereal denom = rr1+rr2+fr3+fr4+fr5;
    std::vector<doublereal> rateC2H2, rateO2, rateOH;
    doublereal coef;

    prodRates.resize(conc.size(),0.0);

    if(denom != 0){
        doublereal ssRatio = (fr1+fr2)/denom;
        doublereal cArea = alpha*Beta_surf*mom[0];
        doublereal cRad = cArea * ssRatio;
       

        /*
         *c2h2
         */
        rateC2H2.resize(nMoments,0.0);
        coef = fr4*cRad;        
        sums(nMoments,2,coef,rateC2H2);        
        /*
         *O2
         */
        rateO2.resize(nMoments,0.0);
        coef = fr5 * cRad;
        sums(nMoments,-2,coef,rateO2);

        prodRates[iC2H2] = -rateC2H2[1]/(2*NA);
        prodRates[iO2] = rateO2[1]/(2*NA);
        prodRates[iH] = ((rr1+fr4)*cRad - fr1*cArea - fr3*cRad)*reducedMoments(4)/NA;
        prodRates[iH2] = (fr1*cArea - rr1*cRad)*reducedMoments(4)/NA;
        prodRates[iH2O] = (fr2*cArea - rr2*cRad)*reducedMoments(4)/NA;


    }else{
        rateC2H2.resize(nMoments,0.0);
        rateO2.resize(nMoments,0.0);
        rateOH.resize(nMoments,0.0);        
    }

    sRates.resize(nMoments,0.0);
    for(int r=0; r<nMoments; r++){
        sRates[r] = rateC2H2[r]+rateO2[r];
    }

}

/*
 *interpolate the whole order reduced moments to evaluate the
 *fractional order reduced moments
 */
void CamSoot::interpolateReducedMoments(std::vector<doublereal>& wom){

    /*
     *order of the highest fractional moment is nMoments+1/6
     */
    CamMath cm;
    int hMoments = nMoments-1;

    reducedMoments.resize(lowFrac,highFrac+1);

    /*
     *negative order moments
     */
    for(int i=lowFrac; i<=-1; i++){
        reducedMoments(i) = cm.interpolateLG(i/6.0,3,prime,wom);
    }

    /*
     *positive order moments
     */
    for(int i=0; i<= 6*hMoments+1; i++){
        if((i%6) != 0.0){
            reducedMoments(i) = cm.interpolateLG(i/6.0,nMoments,prime,wom);
        }else{            
            reducedMoments(i) = wom[i/6];
            
        }
    }

    for(int i=6*hMoments+2; i<=highFrac; i++){
        reducedMoments(i) = 0;
    }
    
}


/*
 *residual definitions
 */
void CamSoot::momentResidual(const doublereal& time, int iMesh_s,  int iMesh_e,
                            int nVar, int nSpec,
                            std::vector<doublereal>& u, std::vector<doublereal>& rho,
                           std::vector<doublereal>& dz,doublereal* y, doublereal* f){
    doublereal convection;
    for(int i=iMesh_s; i<iMesh_e; i++){
        for(int l=0; l<nMoments; l++){
            doublereal phi_e = y[i*nVar+l+nSpec];
            doublereal phi_w = y[(i-1)*nVar+l+nSpec];
            convection = -u[i]*(phi_e - phi_w)/dz[i];
            f[i*nMoments+l] =  convection + wdot(i,l);
            //std::cout << "moment resid " << phi_e << "  " << phi_w << "  " << wdot(i,l) << std::endl;

        }
        //int dd; cin >> dd;
    }

}
/*
 *residual definitions
 */
void CamSoot::momentResidual(const doublereal& time, int iMesh_s, int iMesh_e,
                                std::vector<doublereal>& dz,
                                std::vector<doublereal>& u,
                                std::vector<doublereal>& rho,
                                doublereal* y,
                                doublereal* f){
    doublereal convection;
    for(int i=iMesh_s;i<iMesh_e;i++){
        for(int l=0; l<nMoments;l++){
            doublereal phi_e = y[i*nMoments+l];
            doublereal phi_w = y[(i-1)*nMoments+l];
            convection = -u[i]*(phi_e - phi_w)/dz[i];
            f[i*nMoments+l] = convection + wdot(i,l);
        }
    }

}
void CamSoot::addRates(int nCells, Array2D& rates){
//    for(int i=0; i<nCells; i++){
//        rates(i,iC2H2) += surfProdRate(i,iC2H2)*1e+06;
//        rates(i,iCO) += surfProdRate(i,iCO)*1e+06;
//        rates(i,iH) += surfProdRate(i,iH)*1e+06;
//        rates(i,iH2) += surfProdRate(i,iH2)*1e+06;
//        rates(i,iH2O) += surfProdRate(i,iH2O)*1e+06;
//        rates(i,iO2) += surfProdRate(i,iO2)*1e+06;
//        rates(i,iOH) += surfProdRate(i,iOH)*1e+06;
//    }
    for(int i=0; i<nCells; i++){
        //rates(i,iInception) += surfProdRate(i,iInception);
        //std::cout << i << "  " << rates(i,iC2H2) << "  " << surfProdRate(i,iC2H2) << std::endl;
        rates(i,iC2H2) += surfProdRate(i,iC2H2);
        rates(i,iH) += surfProdRate(i,iH);
        rates(i,iH2) += surfProdRate(i,iH2);
        rates(i,iH2O) += surfProdRate(i,iH2O);
        rates(i,iO2) += surfProdRate(i,iO2);
    }
    //int dd; cin >> dd;
}
