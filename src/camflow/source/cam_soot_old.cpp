/*
 * gridFunction() is is common to both old and new soot functions
 * linear() seems to be used by both old and new. (See commented out part of new)
 * sums() is used by both
 *
 */

void CamSoot::initialize
(
    int nCells,
    Mechanism &mech,
    realVector &mSolnVec
)
{
    /*
     *reduced mass for nucleation for calculating
     *constant part of nucleation. Needs to be sizeMomentsltiplied by sqrt(T)
     *before using (final unit = m^3/(s mol^2). mol^2 is required since
     *this is sizeMomentsltiplied by the concentraction of PAHs
     */
    doublereal rMass = numCAtomInception*cMass/NA;
    constNucleation = 2.2*sqrt(4*PI*kB_cgs/rMass)*dia_PAH*1.0e2*dia_PAH*1.0e2*NA*NA;


    /*
     *Beta for free molecular coagulation, needs to be multiplied by
     *sqrt(T)
     */
    Kf = 2.2*pow(((3*cMass/NA)/(4*PI*rhoSoot)),1.0/6.0)*sqrt(6*kB_cgs/rhoSoot);


    /*
     *constant for continuum coagulation, needs to be multiplied by
     *free path
     */
    doublereal mass = cMass/NA;
    doublereal CD1 = pow((6.0*mass/PI/rhoSoot),1.0/3.0); //unit [=] cm


    Kc_ = 2.514/CD1;


    /*
     *soot hydroxyl oxidation constant [=] cm^3/mol-s
     */
    doublereal oh = ohMass/NA;
    kOH = CD1*CD1* sqrt(PI*kB_cgs/2.0/oh)*NA;


    /*
     *constant for surface growth [=] --
     */
    kSurf = CD1*CD1*lambda*PI;


    /*
     *PAH condensation constant [=] cm/s-mol
     */
    kPAH = 2.2*sqrt(PI)*NA;



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
    realVector temp = mSolnVec;
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


void CamSoot::sootReactions
(
    int cell,
    realVector& conc,
    realVector& mom,
    int nSpec,
    doublereal T,
    doublereal p
)
{
    /*
     *convert all concentration to mol/cm^3
     */
    realVector conc_cgs;
    conc_cgs.resize(conc.size(),0.0);
    for(unsigned int i=0; i<conc.size(); i++)
        conc_cgs[i] = conc[i]*1e-06;

    /*
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
     *nucleation rates
     */
    nucleation(conc_cgs,T,p);
    /*
     *update size moments
     */
    realVector fp;
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
    doublereal ratePAH = 0.0;
    condensation(T,mom[0],conc_cgs[iInception], ratePAH);


    /*
     *surface reactions
     */
    doublereal a,b,rsq;
    realVector surfRates;
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
    	// Need to include surface rate here
        wdot(cell,i) = 0.0; // (nucRate[i] + cgRate[i] + surfRates[i])/mom[i];

    }
}

//----------------------------------------------------//


/*
 *calculate the nucleation rate
 */
void CamSoot::nucleation(realVector& conc, doublereal T, doublereal p){

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

//----------------------------------------------------//

/*
 *coagulation
 */
void CamSoot::coagulation(doublereal T, doublereal p, doublereal M0){


    doublereal kCoag = Kf*sqrt(T);

    realVector f;
    CamMath cm;

    f.clear();
    for(int i=0; i<4; i++){
        f.push_back(gridFunction(i,0,0));//    int dd; cin >> dd;
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

//----------------------------------------------------//



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

//----------------------------------------------------//


/*
 *surface reaction rates
 */
void CamSoot::surface(int hMoment, int cell, doublereal T, doublereal M0,
                        realVector& conc,
                        realVector& totalRates){


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
    realVector rateC2H2, rateO2, rateOH;
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


//----------------------------------------------------//


//----------------------------------------------------//


doublereal CamSoot::betaC1(const int i, const int j)
{

    int i6 = 6*i;
    int j6 = 6*j;

    return sizeMoments(-2+i6)*sizeMoments(2+j6)
           + 2*sizeMoments(i6)*sizeMoments(j6)
           + sizeMoments(2+i6)*sizeMoments(-2+j6);

}

//----------------------------------------------------//

doublereal CamSoot::betaC2(const int i, const int j)
{

    int i6 = 6*i;
    int j6 = 6*j;

    return sizeMoments(-4+i6)*sizeMoments(2+j6)
           + sizeMoments(-2+i6)*sizeMoments(j6)
           + sizeMoments(i6)*sizeMoments(-2+j6)
           + sizeMoments(2+i6)*sizeMoments(-4+j6);

}

//----------------------------------------------------//
