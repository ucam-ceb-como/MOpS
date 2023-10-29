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
 *  New museum Site
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

#include "cam_soot.h"
#include "cam_reporter.h"
#include "cam_setup.h"

using namespace Camflow;
using namespace Sprog;
using namespace Strings;

Array1D CamSoot::sizeMoments;
Array1D CamSoot::reducedMoments;

const double CamSoot::cMass = 0.012; //kg/mol
const double CamSoot::ohMass = 0.017; //kg/mol

CamSoot::CamSoot()
:
    rhoSoot(1.8e+03),        //kg/m^3
    lambda(2.3e15),          //cm^-2 (Keep this in cgs units. We do all surface chem calcs in cgs)
    atomsPerDimer(32),    	 // Override this later as 2 x numCAtomInception
    D1_PAH(2.42e-10),        //m  (This is the size of a benzene ring, not the PAH diameter)
							 // 2.42 = 1.395*sqrt(3), 1.395 A is the C-C in aromatics
    momentON(false),
    lowFrac(-4),
    CD1(pow((6*(cMass/NA)/(PI*rhoSoot)),1.0/3.0))
{

    CamConverter convert;
    std::string fileName("camflow.xml");

    CamXML::Document doc;
    const CamXML::Element* root;
    if(doc.Load(fileName) == 0){
        root = doc.Root();
        readSoot(convert,*root);
    }

}

void CamSoot::readSoot
(
    CamConverter& convert,
    const CamXML::Element& node
)
{

    CamXML::Element *soot;
    CamXML::Element *subnode;
    std::vector<CamXML::Element*> sootSpeciesElement;
    std::vector<CamXML::Element*>::const_iterator p;
    const CamXML::Attribute *atr;
    soot = node.GetFirstChild("soot");

    if (soot != NULL)
    {
        momentON = true;
        soot->GetChildren("species",sootSpeciesElement);
        sootSpecies.clear();
        for (p=sootSpeciesElement.begin(); p<sootSpeciesElement.end(); ++p){
            atr = (*p)->GetAttribute("name");
            if(atr!=NULL){
                sootSpecies.push_back(atr->GetValue());
            }
        }

        subnode = soot->GetFirstChild("inception_species");
        if (subnode!=NULL){
            iSpecies = subnode->GetAttributeValue("name");
        }else{
            throw CamError("Inception species not specified\n");
        }

        subnode = soot->GetFirstChild("nMoments");
        if (subnode!=NULL){
            nMoments = (int(cdble(subnode->Data())));
            highFrac = 6*nMoments+1; //highest fractional moment
        }else{
            throw CamError("Number of moments not specified\n");
        }

        subnode = soot->GetFirstChild("M0threshold");
        if (subnode!=NULL){
            m0Threshold = cdble(subnode->Data());
        }else{
            throw CamError("M0 Threshold not specified\n");
        }

        subnode = soot->GetFirstChild("M0");
        if (subnode != NULL){
            firstMom = cdble(subnode->Data());
            std::cout << "M0 read in " << firstMom << std::endl;
        }else{
            throw CamError("First moment not specified\n");
        }

        subnode = soot->GetFirstChild("dPAH");
        if (subnode!=NULL){
            std::string atrVal = subnode->GetAttributeValue("unit");
            double convertL = convert.getConvertionFactor(atrVal);
            dia_PAH = cdble(subnode->Data())*convertL;
        }

        subnode = soot->GetFirstChild("cPAH");
        if (subnode != NULL){
            numCAtomInception = int(cdble(subnode->Data()));
            atomsPerDimer = 2.0 * numCAtomInception;

        }

        std::string atrVal = soot->GetAttributeValue("regime");

        if (!convertToCaps(trim(atrVal)).compare("FREEMOL"))
            tRegime = FM;
        else if (!convertToCaps(trim(atrVal)).compare("TRANSITION"))
            tRegime = TR;
        else if (!convertToCaps(trim(atrVal)).compare("CONTINUUM"))
            tRegime = CT;
        else
            throw CamError("Unknown transport regime\n");

    }
    else
    {
    	nMoments = 0;
        highFrac = 1;
    }

}

int CamSoot::getRegime() const {
    return tRegime;
}

int CamSoot::getNumMoments() const {
	return nMoments;
}

bool CamSoot::active() const {
    return momentON;
}

double CamSoot::getFirstMoment() const {
	return firstMom;
}

int CamSoot::getAtomsPerDiamer() const {
	return atomsPerDimer;
}


/*
 *calculate the sum of surface growth term
 */
void CamSoot::sums(int hMoment, double massAdded, double coeff,
                        realVector& rates){
    for (int r = 1; r < hMoment; r++) {
        rates[r] = 0.0;
        for (int l = 0; l < r; l++) {
            rates[r] += bnCoeff(l,r)*pow(massAdded,r-l)*reducedMoments(6*l+4);
        }
        rates[r] *= coeff;

    }

}

//------ OLD SOOT FUNCTIONS CULLED FROM HERE ----


double CamSoot::gridFunction(const int k, const int n, const int m)
{
    double gfun = 0;
    for (int l=0; l<=k; ++l)
    {
        int i = 6*(k-l+n);
        int j = 6*(l+m);
        gfun += bnCoeff(l,k)
                *(
                    reducedMoments(i+1)*reducedMoments(j-3)+
                    2*reducedMoments(i-1)*reducedMoments(j-1)+
                    reducedMoments(i-3)*reducedMoments(j+1)
                );
    }
    return gfun;
}



void CamSoot::linear(int n, realVector& y, double& a, double& b,
                        double& rsq)
{
    double sum = 1.0*n;
    double x1 = 0.0;
    double x2 = 0.0;
    double y1 = 0.0;
    double y2 = 0.0;

    double yMean, ess, yss,yPred;

    for (int i=1; i<n; ++i)
    {
        x1 += double(i);
        x2 += double(i)*double(i);

        y1 += y[i];
        y2 += y[i]*double(i);
    }
    double d = sum*x2 - x1*x1;
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


//------------

void CamSoot::setSizeMoments(realVector& fp){
    //realVector fp;
    //fp.resize(nMoments,1.0);
    CamMath cm;

    //for(int i=1; i<nMoments; i++){
    //    fp[i] = conc[i]/conc[0];;
    //}

    sizeMoments.resize(lowFrac,highFrac+1);

    for(int i=lowFrac; i<=-1; i++){
        sizeMoments(i) = cm.interpolateLG(i/6.0,3,prime,fp);
    }


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


//void CamSoot::clearRates(int nCells){
//    surfProdRate.resize(nCells,nMoments);
//void CamSoot::clearRates(){
//    surfProdRate.resize(nMoments);
//}

void CamSoot::report(int nCells)
{
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
        realVector data;
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
//void CamSoot::initMoments(Mechanism &mech, realVector& soln,int nCells){
void CamSoot::initMoments(realVector& soln,int nCells){

    int st = soln.size();
    /*
     *zeroth moment based on user input
     */
    soln.push_back(firstMom);
    /*
     *rest of the moments
     */
    for(int i=1; i<nMoments; i++){

    	// ank25: Do we need to multiply by 1e6 here?
        soln.push_back(soln[i-1+st] * log(double(atomsPerDimer)));

    }
}

void CamSoot::initMomentsConstants(Mechanism &mech){
    /*
     *preparation for interpolation
     */

    CamMath cm;
    cm.prime(nMoments,prime);
    cm.binomCoeff(nMoments,bnCoeff);

   // if(nCells>0){
   //     wdot.resize(nCells,nMoments);
   //     surfProdRate.resize(nCells,mech.SpeciesCount());
   // }
    surfProdRate.resize(mech.SpeciesCount(),0);

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
    double mwPAH = numCAtomInception*cMass;

     //mass of PAH
    double mPAH = mwPAH/NA;
    /*
     *      nucleation kernel  Balthase thesis
     *               |-----------
     *         e_ij *| 8*pi*kB*T
     *               |-----------  (r_i+r_j)^2
     *               | mu_ij
     *              \|
     */

     // ank25:  This looks correct.
     // The 16 (compared to Bathaser's 8) results from 1/2 factor
     // when calculating reduced mass.
     // Units of mPAH are kg per molecule.
     // Units kB are m^2 kg s^-2 K^-1
     Beta_nucl = 2.2*sqrt(4*PI*kB/mPAH)*dia_PAH*dia_PAH;

    /*·
     * Multiplication by NA is carried out to avoid it while
     * evaluating the rates
     */
    Beta_nucl *= NA*NA;
    //std::cout << "beta_nuc  " << Beta_nucl << std::endl;


    /*--------------------------------------------------------------------------
     *                constant for free molecular coagulation
     *-------------------------------------------------------------------------/
     */

    // ank25:  This almost matches eqn 13 in Frenklach 2002.
    // However here we have the additional term cMass/NA
    // Also need to multiply by sqrt(T) later on.
    Beta_fm = 2.2*sqrt(6*kB/rhoSoot)* pow(3*(cMass/NA)/(4*PI*rhoSoot),(1.0/6.0));

    /*--------------------------------------------------------------------------
     *               constant for condensation
     *-------------------------------------------------------------------------/
     */
    Beta_cd = 2.2*sqrt(PI*kB/(2*cMass/NA))*NA;
    //diameter of one particle
    double CD1 = pow((6*(cMass/NA)/(PI*rhoSoot)),1.0/3.0);
    //cout << "D1_PAH after unit conv: " << D1_PAH << endl;

    double CD1_PAH = D1_PAH*sqrt(2.0/3.0);
    double CH2 = Beta_cd*CD1_PAH*CD1_PAH;
    double CHS = Beta_cd* 2* CD1_PAH*CD1;
    double CS2 = Beta_cd*CD1*CD1;
    powPAH.resize(5,nMoments);
    for(int i=1; i<nMoments; i++){
        double x = (2.0*i+1)/2.0;
        powPAH(1,i) = CH2 * pow(double(numCAtomInception),x);
        powPAH(2,i) = CHS * pow(double(numCAtomInception),i);
        x  = (2.0*i-1)/2.0;
        powPAH(3,i) = CS2 * pow(numCAtomInception,x);
    }
    /*--------------------------------------------------------------------------
     *                  constants for surface growth
     *-------------------------------------------------------------------------/
     */
    // Be careful with the units of Beta_surf
    // In this code all surface reactions are done in cgs units
    // and then converted to m3 units at end.
    double CD1_CGS = pow((6*(cMass*1e3/NA)/(PI*rhoSoot*1e-3)),1.0/3.0) ;
    Beta_surf_CGS = lambda*PI*CD1_CGS*CD1_CGS;

    // CBOH is OH+SOOT oxidation constant (again keep in cgs)
    // CBOH = CD1^2 * SQRT(PI * BOLTZMANN / 2 / OH_MASS) * NA ;

    double kB_CGS = kB * 1e7;   // 1.3807D-16;
    double OH_MASS = ohMass*1e3/NA;  // mass per single molecule
    double CBOHCGS = pow(CD1_CGS,2.0) * sqrt(PI * kB_CGS / (2.0 * OH_MASS)) * NA ;
}

void CamSoot::residual
(
    const double& time,
    realVector& wdot,
    double* y,
    double* f
)
{
    for (int r=0; r<nMoments; ++r)
    {
        f[r] = wdot[r];
    }
}

CamSoot::realVector CamSoot::rateAll
(
    const realVector& conc,     //species concentration
    const realVector& moments,  //moments
    const double& T,        //temperature
    const double& p,        //pressure
    const int cellID
)
{

    realVector rates;  // Total rate of change of moments
    //realVector nucRates, coagRates, cdRates, sRates;	// Rate of change of moments
    realVector prodRates; 	// Rate of change of gas phase species due to surface chem
    double prodRatePAHCond; 	// Rate of change of PAH due to condensation

    rates.resize(moments.size(),0.0);
    nucRates.resize(moments.size(),0.0);
    coagRates.resize(moments.size(),0.0);
    cdRates.resize(moments.size(),0.0);
    sRates.resize(moments.size(),0.0);
    prodRates.resize(conc.size(),0.0);

    // Clear the surface production rates and set to zero.
    surfProdRate.clear();
    surfProdRate.resize(conc.size(),0.0);

    // Calculate nucleation rate
    // nucRates is rate of change of moments due to nucleation
    nucRates = rateNucleation(conc[iInception],T);
    // Remove two PAH.
    surfProdRate[iInception]= -2.0 * nucRates[0] / NA;

	// Check the nucleation rates
    //std::cout << "nucRates[0]  " << nucRates[0] << std::endl;
    //std::cout << "surfProdRate[iInception]  " << surfProdRate[iInception] << std::endl;


    // Calculate the interpolated reduced moments
    realVector wholeOrderRedMom;
    wholeOrderRedMom.resize(moments.size(),0.0);
    for(int r=0; r<nMoments; r++){
    	wholeOrderRedMom[r] = moments[r]/moments[0];
    }
    interpolateReducedMoments(wholeOrderRedMom);


    // Calculate coagulation rates
    // coagRates is rate of change of moments due to coagulation
    if (moments[0] > m0Threshold) 		// Threshold check
    	coagRates = rateCoagulation(moments,T);
    else // This else statement is essential.  Ommiting this can result in coagRates being undefined, which leads to negative moments. 
    {
        for(int r=0; r<nMoments; ++r)
          coagRates[r] = 0.0;
    }
	// Check the coagulation rates
    //std::cout << "coagRates[0]  " << coagRates[0] << std::endl;


    // Condensation rates
    // Call this after coagulation, as the reduced moments are calculated in the Coagulation function
    // cdRates is rate of change of moments due to condensation\
    // prodRatePAHCond is rate of change in PAH due to condensation.
    // double tempPAHConc = conc[iInception];

    if (moments[0] > m0Threshold) 		// Threshold check
    {
    	cdRates = rateCondensation(moments, T,conc[iInception]);
    	prodRatePAHCond = -cdRates[1]/numCAtomInception/NA;
    	surfProdRate[iInception] += prodRatePAHCond;
    }
    else  // This else statement is essential.  Ommiting this can result in cdRates being undefined, which leads to negative moments. 
    {
        for(int r=0; r<nMoments; ++r)
          cdRates[r] = 0.0;
    }

	// Check the condensation rates
    //std::cout << "cdRates[0]  " << cdRates[0] << std::endl;
    //std::cout << "prodRatePAHCond  " << prodRatePAHCond << std::endl;

    // Surface rates
    // Call this after coagulation as the reduced moments are calculated in the Coagulation function
    // prodRatesSurface is rate of change of gas phase species due to surface chemistry
    // sRates is rate of change of moments due to surface chemistry

    if (moments[0] > m0Threshold)		// Threshold check
    {
    	rateSurface(conc,T,moments,prodRates,sRates);
    	//rateSurfaceHACARC(conc,T,moments,prodRates,sRates);
        //sRates[0] = 0.0; // Surface growth should not affect M0. However oxidation should decrease the particle count (burnout) therefore commented out on 4/6/21


    // Oxidation regime.
    // TODO: It would be better to include all this in rateSurface function
    // 20 Jan 2022 -- commented out the oxidation sub-code below 
      if(sRates[1] <= 0){

        //std::cout << "In oxidation regime."  << std::endl;

        realVector f;
        f.resize(nMoments,0.0);
        for(int r = 1; r< nMoments; r++){
            f[r] = log(reducedMoments(6*r));
        }

        double a, b, rsq;
        linear(nMoments,f,a,b,rsq);

        for(int n = 4; n<(6*nMoments-2); n += 6){
            reducedMoments(n) = exp(a+b*n/6.0);
        }

        //rateSurface(conc,T,moments,prodRates,sRates);
    	rateSurfaceHACARC(conc,T,moments,prodRates,sRates);

        for(int n=1; n<nMoments; n++){
            sRates[n] *= (moments[n]/(exp(a+b*n)*moments[0]));
        }
        //sRates[0] = 0.0;   Commented out on 4/6/21
      } 
    }   
    else  // This else statement is essential.  Ommiting this can result in sRates being undefined, which leads to negative moments. 
    {
        for(int r=0; r<nMoments; ++r)
          sRates[r] = 0.0;
    }



	// Check the surface rates
    //std::cout << "sRates[0]  " << sRates[0] << std::endl;

    //std::cout << "prodRates[iInception]  " << prodRates[iInception] << std::endl;
    //std::cout << "prodRates[iC2H2]  " << prodRates[iC2H2] << std::endl;
    //std::cout << "prodRates[iCO]  " << prodRates[iCO] << std::endl;
    //std::cout << "prodRates[iO2]  " << prodRates[iO2] << std::endl;
    //std::cout << "prodRates[iH]  " << prodRates[iH] << std::endl;
    //std::cout << "prodRates[iH2]  " << prodRates[iH2] << std::endl;
    //std::cout << "prodRates[iH2O]  " << prodRates[iH2O] << std::endl;
    //std::cout << "prodRates[iOH]  " << prodRates[iOH] << std::endl;

    // Store the surface production rates
    // To be returned to batch reactor via CamSoot::showGasPhaseRates
    // Or used in calculating the residual for flamelets.
    surfProdRate[iInception] += prodRates[iInception];    // RHS is calculated above as part of nucleation step
    surfProdRate[iC2H2] += prodRates[iC2H2];
    surfProdRate[iCO] += prodRates[iCO];
    surfProdRate[iO2] += prodRates[iO2];
    surfProdRate[iH] += prodRates[iH];
    surfProdRate[iH2] += prodRates[iH2];
    surfProdRate[iH2O] += prodRates[iH2O];
    surfProdRate[iOH] += prodRates[iOH];

    rates.resize(nMoments,0.0);
    for (int m=0; m<nMoments; ++m)
    {
    	//rates[m] = (nucRates[m]+coagRates[m]+sRates[m]);
    	rates[m] = (nucRates[m]+coagRates[m]+sRates[m]+cdRates[m]);    /// Base case        
    }
    // Note: This only returns the rate of change of moments.
    // Does not return the rate of change of gas phase species.
    return rates;
}



CamSoot::realVector CamSoot::rateNucleation
(
    const double& concPAH,
    const double& T
)
// Note: This function does not return rate of change of PAH.
// Do this outside this function using:
// 	prodPAH = -2 * nucRate(0) / NA
{

    realVector nucRates;

    double kNucl = Beta_nucl*sqrt(T);
    nucRates.resize(nMoments,0.0);
    /*
     *nucleation rate for the zeroth moment
     */
    nucRates[0] = kNucl* std::max(concPAH,0.0)*concPAH;

    /*
     *number of carbon atoms per dimer
     */
    double cDimer = 2*numCAtomInception;
    /*
     *nucleation rate for higher moments
     */
    for (int m=1; m<nMoments; ++m)
    {
    	nucRates[m] = nucRates[m-1] * cDimer;
    }

    return nucRates;

}

CamSoot::realVector CamSoot::rateCoagulation
(
    const realVector& mom,
    const double& T
)
{

    realVector coagRates;
    double kCoag = Beta_fm * sqrt(T);
    CamMath cm;

    /*
    *ank25:  Move the calculation of interpolated reduced moments
    *out of coagulation and into ratesAll
    realVector wholeOrderRedMom;
    wholeOrderRedMom.resize(mom.size(),0.0);

    for(int r=0; r<nMoments; r++){
    	wholeOrderRedMom[r] = mom[r]/mom[0];
    }
    interpolateReducedMoments(wholeOrderRedMom);
     */

     /*
     *evaluate the grid functions
     */

    realVector f;

    if (nMoments == 6){

    for(int i=0; i<4; i++)
        f.push_back(gridFunction(i,0,0));
    double crk = kCoag*cm.interpolateLG(0.5,4,prime,f);

    f.clear();
    for(int i=0; i<4; i++)
        f.push_back(gridFunction(i,1,1));
    double crk2 = kCoag*cm.interpolateLG(0.5,4,prime,f);

    f.clear();
    for(int i=0; i<4; i++)
        f.push_back(gridFunction(i,1,2));
    double crk3 = kCoag*cm.interpolateLG(0.5,4,prime,f);

    f.clear();
    for(int i=0; i<3; i++)
        f.push_back(gridFunction(i,1,3));
    double crk4a = kCoag*cm.interpolateLG(0.5,3,prime,f);

    f.clear();
    for(int i=0; i<4; i++)
        f.push_back(gridFunction(i,2,2));
    double crk4b = kCoag*cm.interpolateLG(0.5,4,prime,f);

    f.clear();
    for(int i=0; i<2; i++)
        f.push_back(gridFunction(i,1,4));
    double crk5a = kCoag*cm.interpolateLG(0.5,2,prime,f);

    f.clear();
    for(int i=0; i<3; i++)
        f.push_back(gridFunction(i,2,3));
    double crk5b = kCoag*cm.interpolateLG(0.5,3,prime,f);

    double M02 = mom[0]*mom[0];

    coagRates.resize(nMoments,0.0);

    coagRates[0] = -0.5*crk*M02;
    coagRates[1] = 0.0;
    coagRates[2] = crk2 * M02;
    coagRates[3] = 3*crk3 * M02;
    coagRates[4] = (4*crk4a + 3*crk4b) * M02;
    coagRates[5] = (5* crk5a + 10*crk5b) * M02;
    }

    ///-----

    if (nMoments == 5){

    for(int i=0; i<4; i++)
        f.push_back(gridFunction(i,0,0));
    double crk = kCoag*cm.interpolateLG(0.5,4,prime,f);

    f.clear();
    for(int i=0; i<4; i++)
        f.push_back(gridFunction(i,1,1));
    double crk2 = kCoag*cm.interpolateLG(0.5,4,prime,f);

    f.clear();
    for(int i=0; i<3; i++)
        f.push_back(gridFunction(i,1,2));
    double crk3 = kCoag*cm.interpolateLG(0.5,3,prime,f);

    f.clear();
    for(int i=0; i<2; i++)
        f.push_back(gridFunction(i,1,3));
    double crk4a = kCoag*cm.interpolateLG(0.5,2,prime,f);

    f.clear();
    for(int i=0; i<3; i++)
        f.push_back(gridFunction(i,2,2));
    double crk4b = kCoag*cm.interpolateLG(0.5,3,prime,f);


    double M02 = mom[0]*mom[0];

    coagRates.resize(nMoments,0.0);

    coagRates[0] = -0.5*crk*M02;
    coagRates[1] = 0.0;
    coagRates[2] = crk2 * M02;
    coagRates[3] = 3*crk3 * M02;
    coagRates[4] = (4*crk4a + 3*crk4b) * M02;
    }

    ///-----

      if (nMoments == 4){

      for(int i=0; i<4; i++)
          f.push_back(gridFunction(i,0,0));
      double crk = kCoag*cm.interpolateLG(0.5,4,prime,f);

      f.clear();
      for(int i=0; i<3; i++)
          f.push_back(gridFunction(i,1,1));
      double crk2 = kCoag*cm.interpolateLG(0.5,3,prime,f);

      f.clear();
      for(int i=0; i<2; i++)
          f.push_back(gridFunction(i,1,2));
      double crk3 = kCoag*cm.interpolateLG(0.5,2,prime,f);


      double M02 = mom[0]*mom[0];

      coagRates.resize(nMoments,0.0);

      coagRates[0] = -0.5*crk*M02;
      coagRates[1] = 0.0;
      coagRates[2] = crk2 * M02;
      coagRates[3] = 3*crk3 * M02;
      }

      ///-----

        if (nMoments == 3){

        for(int i=0; i<3; i++)
            f.push_back(gridFunction(i,0,0));
        double crk = kCoag*cm.interpolateLG(0.5,3,prime,f);

        f.clear();
        for(int i=0; i<2; i++)
            f.push_back(gridFunction(i,1,1));
        double crk2 = kCoag*cm.interpolateLG(0.5,2,prime,f);

        double M02 = mom[0]*mom[0];

        coagRates.resize(nMoments,0.0);

        coagRates[0] = -0.5*crk*M02;
        coagRates[1] = 0.0;
        coagRates[2] = crk2 * M02;
        }

    return coagRates;

}

/*!
 * Condensation rates
 */
CamSoot::realVector CamSoot::rateCondensation(const realVector& mom,
									const double& T,
                                    const double& concPAH){

    double k_coeff = sqrt(T)*concPAH*mom[0];
    realVector cdRates;
    cdRates.resize(nMoments,0.0);
    for(int r=1; r<nMoments; r++){
        for(int l=0; l<r;l++){
            int l6 = 6*l;
            cdRates[r] += bnCoeff(l,r)*(
                    powPAH(1,r-l)*reducedMoments(l6)+
                    powPAH(2,r-l)*reducedMoments(l6+2)+
                    powPAH(3,r-l)*reducedMoments(l6+4));
        }
        cdRates[r] *= k_coeff;
    }

    // Do this calculation in the calling function
    //prodRatePAHCond = -cdRates[1]/numCAtomInception/NA;

    // Return moment rates
    return cdRates;
}

void CamSoot::rateSurfaceHACARC(const realVector& conc,
                            const double& T,
                            const realVector& mom,
                            realVector& prodRates,
                            realVector& sRates){

    // Balthasar = p20 of his thesis 
    // Appel = p130 of Appel, Henning, Bockhurn, 2000     
    // Mauss = Mauss, Trilkne, Breithbach. p325 of Bockurn (ed)
    // This is coded up as the HACARC mech but with the assumption that f3a = 1.
    // This means that the HACARC reduces to the HACA                    


    double RT = 1.987e-3*T;

    double A_1a_f = 4.2e13;               // Appel (S1 fwd).  Balthasar (1a fwd)
    double A_1a_r = 3.9e12;               // Appel (S1 rev).  Balthasar (1a rev)  
    double A_1b_f = 1.0e10;               // Appel (S2 fwd).  Balthasar (1b fwd)  
    double A_1b_r = 3.68e08;              // Appel (S2 rev).  Balthasar (1b rev) 
    double A_2 = 2e13;                    // Appel (S3). Balthasar (2)
    double A_3a_f = 8.0e7;                // Appel (S4). Balthasar (3a fwd)  
    double A_4a = 2.2e12;                 // Appel (S5). Balthasar (4a and/or 4b)  
    double A_5 = 0.13;    // Reac. prob.  // Appel (S6). Balthasar (5) 

    double k_1a_f = A_1a_f * exp(-13.0/RT);		           // Appel (S1 fwd).  Balthasar (1a fwd)
    double k_1a_r = A_1a_r * exp(-11.0/RT);                // Appel (S1 rev).  Balthasar (1a rev)  
    double k_1b_f = A_1b_f * pow(T,0.734) * exp(-1.43/RT); // Appel (S2 fwd).  Balthasar (1b fwd)  
    double k_1b_r = A_1b_r * pow(T,1.139) * exp(-17.1/RT); // Appel (S2 rev).  Balthasar (1b rev) 
    double k_2 = A_2;                                      // Appel (S3). Balthasar (2)
    double k_3a_f = A_3a_f * pow(T,1.56) * exp(-3.8/RT);   // Appel (S4). Balthasar (3a fwd)  
    double k_4a = A_4a * exp(-7.5/RT);                     // Appel (S5). Balthasar (4a and/or 4b)  
    double k_5 = A_5;   // React. Prob.                    // Appel (S6). Balthasar (5) 

    double r_1a_f = k_1a_f*conc[iH]/1e6;	         // Appel (S1 fwd).  Balthasar (1a fwd)
    double r_1a_r = k_1a_r*conc[iH2]/1e6;            // Appel (S1 rev).  Balthasar (1a rev)  
    double r_1b_f = k_1b_f*conc[iOH]/1e6;            // Appel (S2 fwd).  Balthasar (1b fwd)  
    double r_1b_r = k_1b_r*conc[iH2O]/1e6;           // Appel (S2 rev).  Balthasar (1b rev) 
    double r_2 = k_2 * conc[iH]/1e6;                 // Appel (S3). Balthasar (2)
    double r_3a_f = k_3a_f*conc[iC2H2]/1e6;          // Appel (S4). Balthasar (3a fwd)  
    double r_4a = k_4a*conc[iO2]/1e6;                // Appel (S5). Balthasar (4a and/or 4b)  
    double r_5 = k_5*conc[iOH]/1e6;                  // Appel (S6). Balthasar (5) 


    // The temperature dependent forulations for a and b can result in negatve alpha 
    // so instead use fixed paramters from table 4 of Appel 
    //double par_a = 12.65 - 5.63e-03*T;
    //double par_b = -1.38 + 6.80e-04*T;
    double par_a = 1.1;
    double par_b = -0.06;
    double alpha = tanh(par_a/log10(reducedMoments(6)) + par_b );         // reducedMoments(6) is first reduced whole moment.  

   
    // Eqn 2.51 in Balthasar
    // This formualtion assumes f3a = 1. See eqn 21.36 in Mauss p340 for full eqn. 
    double A_numer = r_1a_f + r_1b_f + r_5 ;
    double A_denom = r_1a_r +r_1b_r + r_2 + r_3a_f  + r_4a + 1e-10;  // Add small epsilon to avoid divide by zero
    double par_A = A_numer / A_denom;
 
    double par_Beta = 0.01;  //  Mauss p340 -- Alternativley can use Beta_surf_CGS as defined above. 
    //double f3a = k3b_f /  (k3b_f +k_3a_r + k_4a* conc[iO2]/1e6)  ;       
    double f3a = 1.0;    // f3a is between 0 and 1.  Setting f3a reduces the mech to HACA
                         // See p22 of Balthasar and p339 of Mauss


    realVector rateC2H2, rateO2, rateOH;		// In cgs units.
    double coef;

    // Rate of change of moments due to C2H2  (eqn 2.47 & 2.48 in Balthasar)
    // Note:  This formulation assumes that f3a  = 1 
    // If f3a !=0 then see p240 in Mauss for formualtion 
    rateC2H2.resize(nMoments,0.0);
    coef = alpha * r_3a_f * f3a * par_A; 

    sums(nMoments,2,coef,rateC2H2);
    rateC2H2[0] = 0;    // Eqn 2.47 in Balthasar

    // Rate of change of moments due to O2  (part of 2.50 & 2.52 in Balthasar)
    rateO2.resize(nMoments,0.0);
    coef = alpha * r_4a * par_A; 
    sums(nMoments,-2,coef,rateO2);
    // Eqn 2.50 in Balthasar and 21.38 p340 in Mauss
    // Note: Baltahasar has omitted Beta (typo?). Mauss has Beta = 0.01
    // reducedMoment(-2) is the negative 1/3 moment 
    rateO2[0] = -r_4a * par_A * par_Beta * reducedMoments(-2);

   // Rate of change of moments due to OH  (part of 2.50 & 2.52 in Balthasar)
    rateOH.resize(nMoments,0.0); 
    coef = alpha * r_5; 
    sums(nMoments,-2,coef,rateOH);
    // Eqn 2.50 in Balthasar and 21.38 p340 in Mauss
    // Note: Baltahasar has omitted Beta (typo?). Mauss has Beta = 0.01
    // reducedMoment(-2) is the negative 1/3 moment 
    rateOH[0] = -r_5 * par_A * par_Beta * reducedMoments(-2);


    // And finally calculate the contribution of each surface reaction
    // to the moments (and convert from cgs back to SI)
    sRates.resize(nMoments,0.0);
    for(int r=0; r<nMoments; r++){
        sRates[r] = (rateC2H2[r]+rateO2[r]+rateOH[r]) * 1e6;
     
    }


    prodRates.resize(conc.size(),0.0);	
 
    // Calcualte the rates of production and consumption of gas phase species 
    double cArea = alpha*par_Beta;    
    double cRad = par_A * alpha * par_Beta;   

    prodRates.resize(conc.size(),0.0);		
    prodRates[iC2H2] = -rateC2H2[1]/(2*NA);
    prodRates[iO2] = rateO2[1]/(2*NA);
    prodRates[iOH] =((r_1b_r*cRad - r_1b_f*cArea) *reducedMoments(4) + rateOH[1]) / NA;     
    prodRates[iH] = ((r_1a_r + r_3a_f )*cRad - r_1a_f *cArea - r_2*cRad)*reducedMoments(4)/NA;
    prodRates[iH2] = (r_1a_f*cArea - r_1a_r*cRad)*reducedMoments(4)/NA;                                  
    prodRates[iH2O] = (r_1b_f*cArea - r_1b_r*cRad)*reducedMoments(4)/NA; 	
    prodRates[iCO]= -(rateO2[1] + rateOH[1]) / NA;

    // Now convert prodRates from cgs to SI
    prodRates[iC2H2] = prodRates[iC2H2] * 1e6;
    prodRates[iO2] = prodRates[iO2] * 1e6;
    prodRates[iOH] = prodRates[iOH] * 1e6;
    prodRates[iH] = prodRates[iH] * 1e6;
    prodRates[iH2] = prodRates[iH2] * 1e6;
    prodRates[iH2O] = prodRates[iH2O] * 1e6;
    prodRates[iCO] = prodRates[iCO] * 1e6;
}

/*!
 * Surface reaction rates
 * This function modifies both prodRates and sRates.
 * prodRates is the rate of change of gas phase species due to surface chem
 * sRates is the rate of change of moments due to surface chem
 */
void CamSoot::rateSurface(const realVector& conc,
                            const double& T,
                            const realVector& mom,
                            realVector& prodRates,
                            realVector& sRates){

    // ank25:  This appears to be roughly the same as the Fortran version
	// Units below are cgs.  Rather than attempting to convert all constants
	// its easier to to just convert the inputs (moments and concentrations)
	// and outputs (rates of moments and rates of species).

	// Divide all concentration by 1e6 below. (SI --> cgs)
    double RT = 1.987e-3*T;
    double fr1 = 4.2e13*exp(-13.0/RT)*conc[iH]/1e6;		 // ank25: was 4.2e12.  Appel has 4.2e13 !!
    double rr1 = 3.9e12*exp(-11.0/RT)*conc[iH2]/1e6;
    double fr2 = 1.0e10 * pow(T,0.734) * exp(-1.43/RT)*conc[iOH]/1e6;
    double rr2 = 3.68e08 * pow(T,1.139) * exp(-17.1/RT)*conc[iH2O]/1e6;
    double fr3 = 2e13 * conc[iH]/1e6;
    double fr4 = 8.0e7*pow(T,1.56)*exp(-3.8/RT)*conc[iC2H2]/1e6;     // ank25: Was 8.0e13 (should be 8e7)
    double fr5 = 2.2e12*exp(-7.5/RT)*conc[iO2]/1e6;
    double fr6 = 0.13*conc[iOH]/1e6;

    // ank25: No unit change here
    double par_a = 12.65 - 5.63e-03*T;
    double par_b = -1.38 + 6.80e-04*T;

    // reducedMoments(6) is first reduced whole moment.
    double alpha = tanh(par_a/log10(reducedMoments(6)) + par_b );

    double denom = rr1+rr2+fr3+fr4+fr5;
    realVector rateC2H2, rateO2, rateOH;		// In cgs units.
    double coef;
    prodRates.resize(conc.size(),0.0);			// In cgs.

    // ank25: rateC2H2, rateO2, rateOH are vectors describing
    // the rate of change to the moments.
    // prodRates is a vector describing the rate of change to gas phase species.

    if(denom != 0){
        double ssRatio = (fr1+fr2)/denom;
        double cArea = alpha*Beta_surf_CGS*mom[0]/1e6;
        // Moment unit conv.  And beta_surf in units of CGS
        double cRad = cArea * ssRatio;

        // C2H2
        rateC2H2.resize(nMoments,0.0);
        coef = fr4*cRad;
        sums(nMoments,2,coef,rateC2H2);

        // O2
        rateO2.resize(nMoments,0.0);
        coef = fr5 * cRad;
        sums(nMoments,-2,coef,rateO2);

        // OH
        rateOH.resize(nMoments,0.0);
        coef = fr6 * CBOHCGS * sqrt(T) * mom[0]/1e6;
        sums(nMoments,-1,coef,rateOH);
        // Moment unit conversion and CBOH in CGS units

        prodRates[iC2H2] = -rateC2H2[1]/(2*NA);
        prodRates[iO2] = rateO2[1]/(2*NA);
        prodRates[iOH] =((rr2*cRad - fr2*cArea) *reducedMoments(4) + rateOH[1]) / NA;
        prodRates[iH] = ((rr1+fr4)*cRad - fr1*cArea - fr3*cRad)*reducedMoments(4)/NA;
        prodRates[iH2] = (fr1*cArea - rr1*cRad)*reducedMoments(4)/NA;
        prodRates[iH2O] = (fr2*cArea - rr2*cRad)*reducedMoments(4)/NA;
  	    prodRates[iCO]= -(rateO2[1] + rateOH[1]) / NA;
    }else{
        rateC2H2.resize(nMoments,0.0);
        rateO2.resize(nMoments,0.0);
        rateOH.resize(nMoments,0.0);
    }

    // Now convert prodRates from cgs to SI
    prodRates[iC2H2] = prodRates[iC2H2] * 1e6;
    prodRates[iO2] = prodRates[iO2] * 1e6;
    prodRates[iOH] = prodRates[iOH] * 1e6;
    prodRates[iH] = prodRates[iH] * 1e6;
    prodRates[iH2] = prodRates[iH2] * 1e6;
    prodRates[iH2O] = prodRates[iH2O] * 1e6;
    prodRates[iCO] = prodRates[iCO] * 1e6;


    // And finally calculate the contribution of each surface reaction
    // to the moments (and convert from cgs to SI)
    sRates.resize(nMoments,0.0);
    for(int r=0; r<nMoments; r++){
        sRates[r] = (rateC2H2[r]+rateO2[r]+rateOH[r]) * 1e6;
    }
}

/*!
 *interpolate the whole order reduced moments to evaluate the
 *fractional order reduced moments
 */
void CamSoot::interpolateReducedMoments(realVector& wom){

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
        //std::cout << "reducedMoments " << i << " : " << reducedMoments(i) << std::endl;
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
 * Don't use this moment residual for Flamelets.
 * Flamlet code provides its own moment residual.
 */
void CamSoot::momentResidual
(
    const double& time,
    const int iMesh_s,
    const int iMesh_e,
    const realVector& dz,
    const realVector& u,
    const realVector& rho,
    const double* y,
    double* f,
    const int nVar,
    const int nSpec
)
{
    double convection;
    for (int i=iMesh_s; i<iMesh_e; ++i)
    {
        for (int l=0; l<nMoments; ++l)
        {
            double phi_e = y[i*nVar+l+nSpec];
            double phi_w = y[(i-1)*nVar+l+nSpec];
            convection = -u[i]*(phi_e - phi_w)/dz[i];
            f[i*nMoments+l] =  convection + wdot(i,l);
        }
    }
}

/*
 * Don't use this moment residual for Flamelets.
 * Flamlet code provides its own moment residual.
 */
void CamSoot::momentResidual
(
    const double& time,
    const int iMesh_s,
    const int iMesh_e,
    const realVector& dz,
    const realVector& u,
    const realVector& rho,
    const double* y,
    double* f
)
{
    momentResidual(time,iMesh_s,iMesh_e,dz,u,rho,y,f,nMoments,0);
}

/*
void CamSoot::addRates(int nCells, Array2D& rates)
{
    for (int i=0; i<nCells; ++i)
    {
        rates(i,iInception) += surfProdRate(i,iInception);
        rates(i,iC2H2) += surfProdRate(i,iC2H2);
        rates(i,iCO) += surfProdRate(i,iCO);
        rates(i,iH) += surfProdRate(i,iH);
        rates(i,iH2) += surfProdRate(i,iH2);
        rates(i,iH2O) += surfProdRate(i,iH2O);
        rates(i,iO2) += surfProdRate(i,iO2);
        rates(i,iOH) += surfProdRate(i,iOH);
//std::cout << i << "  " << rates(i,iC2H2) << "  " << surfProdRate(i,iC2H2) << std::endl;

    }
}
*/

CamSoot::realVector CamSoot::showGasPhaseRates(int nSpecies)
{
    // set rates to zero
	realVector rates;
	rates.resize(nSpecies,0.0);
	rates[iInception] = surfProdRate[iInception];
	rates[iC2H2] = surfProdRate[iC2H2];
	rates[iCO] = surfProdRate[iCO];
	rates[iH] = surfProdRate[iH];
	rates[iH2] = surfProdRate[iH2];
	rates[iH2O] = surfProdRate[iH2O];
	rates[iO2] = surfProdRate[iO2];
	rates[iOH] = surfProdRate[iOH];
    //std::cout << "rates[iInception]  " << rates[iInception] << std::endl;
    //std::cout << "rates[iC2H2]  " << rates[iC2H2] << std::endl;
    //std::cout << "rates[iCO]  " << rates[iCO] << std::endl;
    //std::cout << "rates[iH]  " << rates[iH] << std::endl;
    //std::cout << "rates[iH2]  " << rates[iH2] << std::endl;
    //std::cout << "rates[iH2O]  " << rates[iH2O] << std::endl;
    //std::cout << "rates[iO2]  " << rates[iO2] << std::endl;
    //std::cout << "rates[iOH]  " << rates[iOH] << std::endl;

	return rates;
}


CamSoot::realVector CamSoot::showSootComponentRates(int nMoments)
{
    // set rates to zero
	realVector rates;
	rates.resize(nMoments*4,0.0);

    for (int l=0; l<nMoments; ++l){
	rates[l] = 	nucRates[l];
	rates[l+nMoments] = coagRates[l];
	rates[l+2*nMoments] = cdRates[l];
	rates[l+3*nMoments] = sRates[l];
    }
	return rates;
}

double CamSoot::avgSootDiam()
{
	double sootDiam = CD1 * reducedMoments(2);
	return sootDiam;
}

double CamSoot::dispersion()
{
	double disp = reducedMoments(12) * reducedMoments(6) * reducedMoments(6);
	// double disp = moments(2) * moments(2) * moment(1);
	return disp;
}

double CamSoot::sootSurfaceArea(double M0)
{
    double surfaceArea = PI * CD1 * CD1 * CD1 * reducedMoments(4) * M0;
    return surfaceArea;
}

double CamSoot::sootVolumeFraction(double M0)
{
    double volumeFraction = PI * CD1 * CD1 * CD1 * reducedMoments(6) * M0 / 6.0;
    return volumeFraction;
}



