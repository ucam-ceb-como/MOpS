/*!
  * \author     Zakwan Zainuddin (zz260) && Gustavo Leon (gl413)
  * \file       swp_kmc_mech.cpp
  *
  * \brief        Implementation of the KMCMechanism class declared in swp_kmc_gasph.h
  *
  Project:      sweep (gas-phase chemistry solver).
  Sourceforge:  http://sourceforge.net/projects/mopssuite
  Copyright (C) 2010 Zakwan Zainuddin.

  File purpose:
    Stores the gas phase profiles (time, pressure, temperature, species concentrations)

  Licence:
    This file is part of "sweep".

    Sweep is free software; you can redistribute it and/or
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

#include "swp_kmc_mech.h"
#include "swp_params.h"

#include "choose_index.hpp"
#include "string_functions.h"
#define NDEBUG
#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/random/uniform_01.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace Sweep;
using namespace Sweep::KMC_ARS;
using namespace boost::numeric::ublas;

typedef std::vector<Sweep::KMC_ARS::Reaction> rxnvector;
typedef Sweep::KMC_ARS::KMCGasPoint ggg;

// Constructors
// Default
KMCMechanism::KMCMechanism() {
    m_jplist = obtainJumpProcess();
    m_rates = std::vector<double>(m_jplist.size(),0);
    m_totalrate = 0;
    isACopy = false;
}
//! Copy Constructor
KMCMechanism::KMCMechanism(KMCMechanism& m) {
    m_jplist = m.m_jplist;
    m_rates = m.m_rates;
    m_totalrate = m.m_totalrate;
    isACopy = true;
}
//! Destructor
KMCMechanism::~KMCMechanism() {
    if (m_jplist.size()!=0 && !isACopy){
        {
            for (size_t i=0;i!=m_jplist.size();++i)
            delete m_jplist[i];
        }
        m_jplist.clear();
    }
}


//! WRITE PROCESSES
//! Load processes from process list
void KMCMechanism::loadProcesses(std::vector<JumpProcess*> (*jp)()) {//how##
    m_jplist = jp();
}
// KMC Algorithm and Read Processes
////! Calculates each jump rates and stores in a vector
//rvector KMCMechanism::calculateRates(std::vector<JumpProcess*>& p) const {
//    // temporary vector to store rates of each JumpProcess
//    rvector temp; 
//    for(int i=0; i!= (int)p.size(); i++) {
//        double a = p[i]->getRate(); // gets rate
//        temp.push_back(a); // store in temporary vector
//    }
//    return temp;
//}
//! Calculates total jump rates
//double KMCMechanism::getTotalRate(const rvector& jrate) {
//    // calculates sum of all elements of jrate
//    double sum=0;
//    for(int i=0; i!=(int) jrate.size(); i++) {
//        sum += jrate[i];
//    }
//    return sum;
//}

//! Choosing a reaction to be taken place, returns pointer to jump process
ChosenProcess KMCMechanism::chooseReaction(rng_type &rng) const {
    // chooses index from a vector of weights (double number in this case) randomly
    boost::uniform_01<rng_type &, double> uniformGenerator(rng);
    size_t ind = chooseIndex<double>(m_rates, uniformGenerator);
    return ChosenProcess(m_jplist[ind], ind);
}

typedef Sweep::KMC_ARS::KMCGasPoint sp;

//! Returns a vector of jump processes implemented in model.
std::vector<JumpProcess*> KMCMechanism::obtainJumpProcess(){
    std::vector<JumpProcess*> temp;
    //! Initialise all jump processes.
    JumpProcess* j_G6R_AC = new G6R_AC; j_G6R_AC->initialise();									//!  1- R6 Growth on AC [AR1].
    JumpProcess* j_G6R_FE = new G6R_FE; j_G6R_FE->initialise();									//!  2- R6 Growth on FE [AR2].
    JumpProcess* j_L6_BY6 = new L6_BY6; j_L6_BY6->initialise();									//!  3- BY6 closure [AR14].
    JumpProcess* j_PH_benz = new PH_benz; j_PH_benz->initialise();								//!  4- phenyl addition [AR15].
    JumpProcess* j_D6R_FE3 = new D6R_FE3; j_D6R_FE3->initialise();								//!  5- R6 Desorption at FE [AR8].
    //JumpProcess* j_O6R_FE3_O2 = new O6R_FE3_O2; j_O6R_FE3_O2->initialise();						//!  6- R6 Oxidation at FE by O2 [AR10].
    //JumpProcess* j_O6R_FE3_OH = new O6R_FE3_OH; j_O6R_FE3_OH->initialise();						//!  7- R6 Oxidation at FE by OH [AR11].
    JumpProcess* j_O6R_FE_HACA = new O6R_FE_HACA; j_O6R_FE_HACA->initialise();					//!  8- R6 Oxidation at FE-HACA [AR12].
    //JumpProcess* j_O6R_FE_HACA_OH = new O6R_FE_HACA_OH; j_O6R_FE_HACA_OH->initialise();			//!  9- R6 Oxidation at AC by OH [AR13].
    JumpProcess* j_G5R_ZZ = new G5R_ZZ; j_G5R_ZZ->initialise();									//! 10- R5 growth at ZZ [AR3].
    JumpProcess* j_D5R_R5 = new D5R_R5; j_D5R_R5->initialise();									//! 11- R5 desorption [AR7].
    JumpProcess* j_C6R_AC_FE3 = new C6R_AC_FE3; j_C6R_AC_FE3->initialise();						//! 12- R6 conversion to R5 [AR9].
    JumpProcess* j_C5R_RFE = new C5R_RFE; j_C5R_RFE->initialise();                              //! 13- R5 conversion to R6 on RFE [AR5].
    JumpProcess* j_C5R_RAC = new C5R_RAC; j_C5R_RAC->initialise();                              //! 14- R5 conversion to R6 on AC [AR4].
    JumpProcess* j_M5R_RZZ = new M5R_RZZ; j_M5R_RZZ->initialise();                              //! 15- R5 migration to neighbouring ZZ [AR6]. 
    JumpProcess* j_C6R_BY5_FE3 = new C6R_BY5_FE3; j_C6R_BY5_FE3->initialise();                  //! 16- R6 migration & conversion to R5 at BY5 [AR22].
    JumpProcess* j_C6R_BY5_FE3violi = new C6R_BY5_FE3violi; j_C6R_BY5_FE3violi->initialise();   //! 17- R6 migration & conversion to R5 at BY5 (violi) [AR24].
    JumpProcess* j_L5R_BY5 = new L5R_BY5; j_L5R_BY5->initialise();                              //! 18- BY5 closure [AR16].
    JumpProcess* j_M6R_BY5_FE3 = new M6R_BY5_FE3; j_M6R_BY5_FE3->initialise();                  //! 19- R6 desorption at bay -> pyrene [AR21].
    JumpProcess* j_O6R_FE2_side = new O6R_FE2_side; j_O6R_FE2_side->initialise();               //! 20- R6 Oxidation at FE2 by side reaction.
    JumpProcess* j_O6R_FE2_top = new O6R_FE2_top; j_O6R_FE2_top->initialise();                  //! 21- R6 Oxidation at FE2 by top reaction.
	JumpProcess* j_D6R_FE_AC = new D6R_FE_AC; j_D6R_FE_AC->initialise();						//! 22- R6 Desorption from FE to AC site.
	JumpProcess* j_B6R_ACR5 = new B6R_ACR5; j_B6R_ACR5->initialise();                           //!< 23 - Bay-capping.
	JumpProcess* j_M5R_ACR5_ZZ = new M5R_ACR5_ZZ; j_M5R_ACR5_ZZ->initialise();                  //!< 24 - Embedded 5-member ring migration to ZZ. //Moved to instant jump process list.
	JumpProcess* j_G6R_RZZ = new G6R_RZZ; j_G6R_RZZ->initialise();                              //!< 25 - R6 growth on RZZ.
	JumpProcess* j_G6R_RFER = new G6R_RFER; j_G6R_RFER->initialise();                           //!< 26 - R6 growth on RFER.
	//JumpProcess* j_G6R_R5 = new G6R_R5; j_G6R_R5->initialise();                                 //!< 27 - R6 growth on R5.
	JumpProcess* j_L6_RBY5 = new L6_RBY5; j_L6_RBY5->initialise();                              //!< 28 - RBY5 closure reaction.
	JumpProcess* j_L6_RACR = new L6_RACR; j_L6_RACR->initialise();                              //!< 29 - RACR closure reaction.
	//JumpProcess* j_G5R_RFE = new G5R_RFE; j_G5R_RFE->initialise();                              //!< 30 - R5 growth on RFE.
	JumpProcess* j_C6R_RAC_FE3 = new C6R_RAC_FE3; j_C6R_RAC_FE3->initialise();                  //!< 31 - R6 migration & conversion to R5 at RAC.
	JumpProcess* j_C6R_RAC_FE3violi = new C6R_RAC_FE3violi; j_C6R_RAC_FE3violi->initialise();   //!< 32 - R6 migration & conversion to R5 at RAC.
	JumpProcess* j_M6R_RAC_FE3 = new M6R_RAC_FE3; j_M6R_RAC_FE3->initialise();                  //!< 33 - R6 desorption at RAC -> pyrene.
	JumpProcess* j_MR5_R6 = new MR5_R6; j_MR5_R6->initialise();                                 //!< 34 - R5 exchange with R6. //Moved to instant jump process list.
	JumpProcess* j_GR7_R5R6AC = new GR7_R5R6AC; j_GR7_R5R6AC->initialise();                        //!< 35 - R7 growth on R5R6AC.
	JumpProcess* j_GR7_FEACR5 = new GR7_FEACR5; j_GR7_FEACR5->initialise();                        //!< 36 - R7 growth on FEACR5.
	JumpProcess* j_G6R_R5R6ZZ = new G6R_R5R6ZZ; j_G6R_R5R6ZZ->initialise();                        //!< 37 - R6 growth on R5R6ZZ.
	JumpProcess* j_L7_ACACR5 = new L7_ACACR5; j_L7_ACACR5->initialise();                        //!< 38 - R7 bay closure on ACACR5.
	JumpProcess* j_G6R_R5R6FER = new G6R_RFER; j_G6R_R5R6FER->initialise();                       //!< 39 - R6 growth on RFER.
	JumpProcess* j_G6R_R5R6FER5R6 = new G6R_R5R6FER5R6; j_G6R_R5R6FER5R6->initialise();           //!< 40 - R6 growth on RFER.
	JumpProcess* j_L7_FEZZACR5 = new L7_FEZZACR5; j_L7_FEZZACR5->initialise();                    //!< 41 - R7 bay closure on FEZZACR5.
	JumpProcess* j_C5R_RZZR = new C5R_RZZR; j_C5R_RZZR->initialise();                    		 //!< 42 - R5 conversion to R6 on RZZR. (2R5 collision)
	JumpProcess* j_C5R_R5R6ZZR = new C5R_R5R6ZZR; j_C5R_R5R6ZZR->initialise();                   //!< 43 - R5 conversion to R6 on R5R6ZZR. (2R5 collision)
	JumpProcess* j_L6_R5R6BY5 = new L6_R5R6BY5; j_L6_R5R6BY5->initialise();                      //!< 44 - R5R6BY5 closure reaction.
	JumpProcess* j_L6_R5R6ACR = new L6_R5R6ACR; j_L6_R5R6ACR->initialise();                      //!< 45 - R5R6ACR closure reaction.
	JumpProcess* j_L6_R5R6ACR5R6 = new L6_R5R6ACR5R6; j_L6_R5R6ACR5R6->initialise();             //!< 46 - R5R6ACR5R6 closure reaction.
	JumpProcess* j_L6_ZZACR5 = new L6_ZZACR5; j_L6_ZZACR5->initialise();             			 //!< 47 - ZZACR5 closure reaction.
	JumpProcess* j_L6_R5FEACR5 = new L6_R5FEACR5; j_L6_R5FEACR5->initialise();             		 //!< 48 - R5FEACR5 closure reaction.
	JumpProcess* j_L6_FEACR5FE = new L6_FEACR5FE; j_L6_FEACR5FE->initialise();             		 //!< 49 - FEACR5FE closure reaction.
	JumpProcess* j_L6_R5ACR5R5 = new L6_R5ACR5R5; j_L6_R5ACR5R5->initialise();             		 //!< 50 - R5ACR5R5 closure reaction.
	JumpProcess* j_L7_R5ZZACR5 = new L7_R5ZZACR5; j_L7_R5ZZACR5->initialise();                   //!< 51 - R7 bay closure on R5ZZACR5.
	JumpProcess* j_L6_ACR5R5R6 = new L6_ACR5R5R6; j_L6_ACR5R5R6->initialise();             		 //!< 52 - ACR5R5R6 closure reaction.
	JumpProcess* j_L7_ACR5R5R6ZZ = new L7_ACR5R5R6ZZ; j_L7_ACR5R5R6ZZ->initialise();             //!< 53 - R7 bay closure on ACR5R5R6ZZ.
	JumpProcess* j_A_CH3 = new A_CH3; j_A_CH3->initialise();             						 //!< 54 - CH3 addition.
	JumpProcess* j_D_CH3 = new D_CH3; j_D_CH3->initialise();             						 //!< 55 - CH3 desorption.*/
	JumpProcess* j_O5R_R5R6 = new O5R_R5R6; j_O5R_R5R6->initialise();             				//!< 56 - Oxidation of R5R6 site.
	JumpProcess* j_O5R_R5R6ZZ = new O5R_R5R6ZZ; j_O5R_R5R6ZZ->initialise();             		//!< 57 - Oxidation of R5R6ZZ site.
	JumpProcess* j_O5R_R5R6AC = new O5R_R5R6AC; j_O5R_R5R6AC->initialise();             		//!< 58 - Oxidation of R5R6AC site.
	JumpProcess* j_O5R_R5R6BY5 = new O5R_R5R6BY5; j_O5R_R5R6BY5->initialise();             		//!< 59 - Oxidation of R5R6BY5 site.
	JumpProcess* j_O5R_R5R6FER = new O5R_R5R6FER; j_O5R_R5R6FER->initialise();             		//!< 60 - Oxidation of R5R6FER site.
	JumpProcess* j_O5R_R5R6ZZR = new O5R_R5R6ZZR; j_O5R_R5R6ZZR->initialise();             		//!< 61 - Oxidation of R5R6ZZR site.
	JumpProcess* j_O5R_R5R6ACR = new O5R_R5R6ACR; j_O5R_R5R6ACR->initialise();             		//!< 62 - Oxidation of R5R6ACR site.
	JumpProcess* j_O5R_R5R6FER5R6 = new O5R_R5R6FER5R6; j_O5R_R5R6FER5R6->initialise();         //!< 63 - Oxidation of R5R6FER5R6 site.
	JumpProcess* j_O5R_R5R6ZZR5R6 = new O5R_R5R6ZZR5R6; j_O5R_R5R6ZZR5R6->initialise();         //!< 64 - Oxidation of R5R6ZZR5R6 site.
	JumpProcess* j_O5R_R5R6ACR5R6 = new O5R_R5R6ACR5R6; j_O5R_R5R6ACR5R6->initialise();         //!< 65 - Oxidation of R5R6ACR5R6 site.
	JumpProcess* j_M5R_ACR5_ZZ_ZZ = new M5R_ACR5_ZZ_ZZ; j_M5R_ACR5_ZZ_ZZ->initialise();         //!< 66 - ACR5 migration to both sides.
	JumpProcess* j_MR5R7_edge = new MR5R7_edge; j_MR5R7_edge->initialise();        				//!< 67 - R5R7 pair edge healing.
    
       
	//! Jump processes included in the model (Comment out any process to be omitted).
    temp.push_back(j_G6R_AC);            //!  1- R6 Growth on AC [AR1].
    temp.push_back(j_G6R_FE);            //!  2- R6 Growth on FE [AR2].
    temp.push_back(j_L6_BY6);            //!  3- BY6 closure [AR14].
    temp.push_back(j_PH_benz);           //!  4- phenyl addition [AR15].
    temp.push_back(j_D6R_FE3);           //!  5- R6 Desorption at FE [AR8].
    //temp.push_back(j_O6R_FE3_O2);        //!  6- R6 Oxidation at FE by O2 [AR10].
    //temp.push_back(j_O6R_FE3_OH);        //!  7- R6 Oxidation at FE by OH [AR11].
    temp.push_back(j_O6R_FE_HACA);    //!  R6 Oxidation at FE-HACA [AR12].
    //temp.push_back(j_O6R_FE_HACA_OH);    //!  9- R6 Oxidation at AC by OH [AR13].
    temp.push_back(j_G5R_ZZ);            //! 10- R5 growth at ZZ [AR3].
    temp.push_back(j_D5R_R5);            //! 11- R5 desorption [AR7].
    temp.push_back(j_C6R_AC_FE3);        //! 12- R6 conversion to R5 [AR9].
    temp.push_back(j_C5R_RFE);           //! 13- R5 conversion to R6 on RFE [AR5].
    temp.push_back(j_C5R_RAC);           //! 14- R5 conversion to R6 on AC [AR4].
    temp.push_back(j_M5R_RZZ);           //! 15- R5 migration to neighbouring ZZ [AR6]. 
    temp.push_back(j_C6R_BY5_FE3);       //! 16- R6 migration & conversion to R5 at BY5 [AR22].
    temp.push_back(j_C6R_BY5_FE3violi);  //! 17- R6 migration & conversion to R5 at BY5 (violi) [AR24].
    temp.push_back(j_L5R_BY5);           //! 18- BY5 closure [AR16].
    temp.push_back(j_M6R_BY5_FE3);       //! 19- R6 desorption at bay -> pyrene [AR21].
    temp.push_back(j_O6R_FE2_side);        //! 20- R6 Oxidation at FE2 by side reaction.
	temp.push_back(j_O6R_FE2_top);        //! 21- R6 Oxidation at FE2 by top reaction.
	temp.push_back(j_D6R_FE_AC);         //! 22- R6 Desorption from FE to AC site.
	temp.push_back(j_B6R_ACR5);         //!< 23 - Bay-capping.
	temp.push_back(j_M5R_ACR5_ZZ);      //!< 24 - Embedded 5-member ring migration to ZZ. //Moved to instant jump process list.
	temp.push_back(j_G6R_RZZ);          //!< 25 - R6 growth on RZZ.
	temp.push_back(j_G6R_RFER);         //!< 26 - R6 growth on RFER.
	//temp.push_back(j_G6R_R5);           //!< 27 - R6 growth on R5.
	temp.push_back(j_L6_RBY5);          //!< 28 - RBY5 closure reaction.
	temp.push_back(j_L6_RACR);          //!< 29 - RACR closure reaction.
	//temp.push_back(j_G5R_RFE);          //!< 30 - R5 growth on RFE.
	temp.push_back(j_C6R_RAC_FE3);      //!< 31 - R6 migration & conversion to R5 at RAC.
	temp.push_back(j_C6R_RAC_FE3violi); //!< 32 - R6 migration & conversion to R5 at RAC.
	temp.push_back(j_M6R_RAC_FE3);      //!< 33 - R6 desorption at RAC -> pyrene.
	temp.push_back(j_MR5_R6);           //!< 34 - R5 exchange with R6. //Moved to instant jump process list.
	temp.push_back(j_GR7_R5R6AC);           //!< 35 - R7 growth on R5R6AC.
	temp.push_back(j_GR7_FEACR5);           //!< 36 - R7 growth on FEACR5.
	temp.push_back(j_G6R_R5R6ZZ);          //!< 37 - R6 growth on R5R6ZZ.
	temp.push_back(j_L7_ACACR5);          //!< 38 - R7 bay closure on ACACR5.
	temp.push_back(j_G6R_R5R6FER);         //!< 39 - R6 growth on RFER.
	temp.push_back(j_G6R_R5R6FER5R6);      //!< 40 - R6 growth on RFER5R6.
	temp.push_back(j_L7_FEZZACR5);          //!< 41 - R7 bay closure on FEZZACR5.
	temp.push_back(j_C5R_RZZR);          //!< 42 - R5 conversion to R6 on RZZR. (2R5 collision)
	temp.push_back(j_C5R_R5R6ZZR);        //!< 43 - R5 conversion to R6 on R5R6ZZR. (2R5 collision)
	temp.push_back(j_L6_R5R6BY5); 		  //!< 44 - R5R6BY5 closure reaction.
	temp.push_back(j_L6_R5R6ACR);         //!< 45 - R5R6ACR closure reaction.
	temp.push_back(j_L6_R5R6ACR5R6);      //!< 46 - R5R6ACR5R6 closure reaction.
	temp.push_back(j_L6_ZZACR5);          //!< 47 - ZZACR5 closure reaction.
	temp.push_back(j_L6_R5FEACR5);        //!< 48 - R5FEACR5 closure reaction.
	temp.push_back(j_L6_FEACR5FE);        //!< 49 - FEFEACR5FE closure reaction.
	temp.push_back(j_L6_R5ACR5R5);        //!< 50 - R5ACR5R5 closure reaction.
	temp.push_back(j_L7_R5ZZACR5);          //!< 51 - R7 bay closure on R5ZZACR5.
	temp.push_back(j_L6_ACR5R5R6);        //!< 52 - ACR5R5R6 closure reaction.
	temp.push_back(j_L7_ACR5R5R6ZZ);          //!< 53 - R7 bay closure on ACR5R5R6ZZ.
	temp.push_back(j_A_CH3);          		//!< 54 - CH3 addition.
	temp.push_back(j_D_CH3);          		//!< 55 - CH3 desorption.*/
	temp.push_back(j_O5R_R5R6);          		//!< 56 - Oxidation of R5R6 site.
	temp.push_back(j_O5R_R5R6ZZ);               //!< 57 - Oxidation of R5R6ZZ site.
	temp.push_back(j_O5R_R5R6AC);               //!< 58 - Oxidation of R5R6AC site.
	temp.push_back(j_O5R_R5R6BY5);              //!< 59 - Oxidation of R5R6BY5 site.
	temp.push_back(j_O5R_R5R6FER);          		//!< 60 - Oxidation of R5R6FER site.
	temp.push_back(j_O5R_R5R6ZZR);          		//!< 61 - Oxidation of R5R6ZZR site.
	temp.push_back(j_O5R_R5R6ACR);          		//!< 62 - Oxidation of R5R6ACR site.
	temp.push_back(j_O5R_R5R6FER5R6);          		//!< 63 - Oxidation of R5R6FER5R6 site.
	temp.push_back(j_O5R_R5R6ZZR5R6);          		//!< 64 - Oxidation of R5R6ZZR5R6 site.
	temp.push_back(j_O5R_R5R6ACR5R6);          		//!< 65 - Oxidation of R5R6ACR5R6 site.
	temp.push_back(j_M5R_ACR5_ZZ_ZZ);      			//!< 66 - ACR5 migration to both sides.
	temp.push_back(j_MR5R7_edge);          		//!< 69 - Oxidation of R5R6ACR5R6 site.

    return temp;
}

//! Calculates jump rate for each jump process
void KMCMechanism::calculateRates(const KMCGasPoint& gp, 
                    PAHProcess& st, 
                    const double& t) {
    double temp=0;
    double pressure = gp[gp.P]/1e5;
    // Choose suitable mechanism according to P
    if(pressure > 0.5 && pressure <= 5) { // mechanism at 1 atm 
        for(int i = 0; i!= (int) m_jplist.size() ; i++) {
			(m_jplist[i])->calculateElemRxnRate((m_jplist[i])->getVec1(), gp);
			m_rates[i] = (m_jplist[i])->setRate1(gp, st/*, t*/);
			temp += m_rates[i];
        }
    }else if(pressure > 0.01 && pressure <= 0.07) { // mechanism at 0.0267atm
        for(int i = 0; i!= (int) m_jplist.size() ; i++) {
			(m_jplist[i])->calculateElemRxnRate((m_jplist[i])->getVec0p0267(), gp);
			m_rates[i] = (m_jplist[i])->setRate0p0267(gp, st/*, t*/);
			temp += m_rates[i];
        }
    }else if(pressure > 0.07 && pressure <= 0.5) { // mechanism at 0.12atm
        for(int i = 0; i!= (int) m_jplist.size() ; i++) {
			(m_jplist[i])->calculateElemRxnRate((m_jplist[i])->getVec1(), gp); // As a test, decided to use the mechanism for 1atm that has been debugged.
			//(m_jplist[i])->calculateElemRxnRate((m_jplist[i])->getVec0p12(), gp);
			m_rates[i] = (m_jplist[i])->setRate1(gp, st/*, t*/);
			//m_rates[i] = (m_jplist[i])->setRate0p12(gp, st/*, t*/);
			temp += m_rates[i];
        }
    }else std::cout<<"ERROR: No reaction mechanism for this pressure condition.\n";
    // update total rates
    if (temp < 1e-20) temp = 1e-20;
    m_totalrate = temp;
}

//! Returns vector of jump processes
const std::vector<JumpProcess*>& KMCMechanism::JPList() const {
    return m_jplist;
}

//! Returns vector of jump rates
const std::vector<double>& KMCMechanism::Rates() const {
    return m_rates;
}

//! Returns total rates
double KMCMechanism::TotalRate() const {
    return m_totalrate;
}

//! Process list (rate calculations, energy units in kcal)
// For C++ jump process ID X, search for IDX
// For Matlab jump process X, search for ARX
// TEST: 4.563046 -> 4.56937799
//        11.51512 -> 11.53110048
//        13.354675 -> 13.37320574
//        38.463376 -> 38.51674641
// ************************************************************
// ID1- R6 growth on AC (AR1 on Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void G6R_AC::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(2.50e14, 0, 16.00, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.40e9, 0.88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.10e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.35468, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(1.87e7, 1.787, 3.262, sp::C2H2));   //5 - r4f
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(2.50e14, 0, 16.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.40e9, 0.88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.10e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(1.87e7, 1.787, 3.262, sp::C2H2));   //5 - r4f
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	//------------Reactions for A3-4----------------------------
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4              - Forward
	//Using rates by On the low-temperature limit of HACA. Frenklach et al. 2018
	addReaction(rxnV3, Reaction(1.190E+22, -2.450E+00, 1.889E+01, sp::C2H2));         // A3* + C2H2 -> A4 + H              	- 5              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(1.060E+14, -4.900E-01, 8.204E+00, sp::C2H2));         // A3* + H -> A4 + H              	- 6              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 7              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 8              - Frenklach et al. 2018
	//Method using QSSA from ABF:
	/*addReaction(rxnV3, Reaction(1.000E+14, 0.000E+00, 0.000E+00, sp::None));       // A3-4 + H (+M) <=> A3 (+M)         	- 4              - Forward - High rate
	addReaction(rxnV3, Reaction(2.100E+139, -3.480E+01, 1.838E+01, sp::None));       // A3-4 + H (+M) <=> A3 (+M)           - 5              - Forward - High rate*/
	//Parameters = [  1.00000000e-03   1.71400000e+02   1.71400000e+02   4.99280000e+03]
	/*addReaction(rxnV3, Reaction(2.315E+31, -4.589E+00, 1.185E+02, sp::None));                // A3-4 + H (+M) <=> A3 (+M)           - 6              - Backward
	addReaction(rxnV3, Reaction(8.000E+61, -1.450E+01, 3.480E+01, sp::C2H2));                // A3-4 + C2H2 <=> A3C2H2              - 7              - Forward
	addReaction(rxnV3, Reaction(4.107E+69, -1.626E+01, 6.773E+01, sp::None));                // A3-4 + C2H2 <=> A3C2H2              - 8              - Backward
	addReaction(rxnV3, Reaction(1.200E+26, -3.440E+00, 3.020E+01, sp::C2H2));                // A3-4 + C2H2 <=> A3C2H + H           - 9              - Forward
	addReaction(rxnV3, Reaction(6.319E+33, -5.219E+00, 3.169E+01, sp::H));           // A3-4 + C2H2 <=> A3C2H + H           - 10             - Backward
	addReaction(rxnV3, Reaction(6.600E+24, -3.360E+00, 1.780E+01, sp::C2H2));        // A3-4 + C2H2 <=> A4 + H              - 11             - Forward
	addReaction(rxnV3, Reaction(1.815E+34, -4.875E+00, 7.521E+01, sp::H));           // A3-4 + C2H2 <=> A4 + H              - 12             - Backward
	////------------Reactions for A3C2H----------------------------
	addReaction(rxnV3, Reaction(5.000E+13, 0.000E+00, 0.000E+00, sp::C2H));          // A3 + C2H <=> A3C2H + H              - 13             - Forward
	addReaction(rxnV3, Reaction(9.332E+20, -1.526E+00, 2.401E+01, sp::H));           // A3 + C2H <=> A3C2H + H              - 14             - Backward
	addReaction(rxnV3, Reaction(1.900E+64, -1.512E+01, 2.930E+01, sp::H));           // A3C2H + H <=> A3C2H2                - 15             - Forward
	addReaction(rxnV3, Reaction(1.852E+64, -1.510E+01, 6.074E+01, sp::None));        // A3C2H + H <=> A3C2H2                - 16             - Backward
	////------------Reactions for A3C2H2----------------------------
	addReaction(rxnV3, Reaction(2.000E+63, -1.528E+01, 4.320E+01, sp::None));                // A3C2H2 <=> A4 + H           - 17             - Forward
	addReaction(rxnV3, Reaction(1.071E+65, -1.504E+01, 6.768E+01, sp::H));           // A3C2H2 <=> A4 + H           - 18             - Backward
	////------------Reactions for A4----------------------------
	addReaction(rxnV3, Reaction(9.000E+38, -7.390E+00, 2.070E+01, sp::H));           // A3C2H + H <=> A4 + H                - 19             - Forward
	addReaction(rxnV3, Reaction(4.700E+40, -7.125E+00, 7.662E+01, sp::H));           // A3C2H + H <=> A4 + H                - 20             - Backward
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));            // A4 + H <=> A4- + H2         - 21             - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));           // A4 + H <=> A4- + H2         - 22             - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A4 + OH <=> A4- + H2O               - 23             - Forward
	addReaction(rxnV3, Reaction(9.600E+04, 2.162E+00, 9.780E+00, sp::H2O));          // A4 + OH <=> A4- + H2O               - 24             - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A4* + H -> A4              			- 25              - Forward*/
	//addReaction(rxnV3, Reaction(1.000E+14, 0.000E+00, 0.000E+00, sp::H));            // A4- + H <=> A4              - 25             - Forward
	//addReaction(rxnV3, Reaction(4.635E+18, -7.927E-01, 1.113E+02, sp::None));                // A4- + H <=> A4              - 26             - Backward
	////------------Reactions for A4-----------------------------	
	/*
    addReaction(rxnV3, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));   //1 - r1b
    addReaction(rxnV3, Reaction(1.0e10, 0.734, 1.43, sp::OH));  //2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV3, Reaction(2.0e13, 0, 0, sp::H));  //4 - r3f
    addReaction(rxnV3, Reaction(8.0e7, 1.56, 3.8, sp::C2H2));   //5 - r4f
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
	*/

    m_sType = AC; // sitetype
    m_name = "G6R at AC"; // name of process
    m_ID = 1;
}
// Jump rate calculation
double G6R_AC::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = 2*m_r[5]*r_f* site_count; // Rate Equation
}
double G6R_AC::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double G6R_AC::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]+m_r[8]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = (m_r[5]+m_r[6])*r_f* site_count; // Rate Equation
	
	/*//ABF. Previous version.
	// Calculating R3b3
	// Parameters alpha = 0.001, T3 = 117.4, T1 = 171.4, T2 = 4992.8
	double T = gp[gp.T];
	double P = gp[gp.P]; //Pressure in Pa
	double M_eff1 = P / 8.314E+06 / T;
	double M = M_eff1 + (1 * gp[gp.H2]) + (5 * gp[gp.H2O]) + (1 * gp[gp.CH4]) + (0.5*gp[gp.CO]) + (1 * gp[gp.CO2]) + (2 * gp[gp.C2H6]);
	double k_o = m_r[5];
	double k_inf = m_r[4];
	double F_cent = (1 - 0.001) * exp(-T / 117.4) + 0.001*exp(-T / 171.4) + exp(-4992.8 / T);
	double logF_cent = log10(F_cent);
	double r_c = k_inf / (1 + (k_inf / (k_o*M)));
	double N0 = 0.75 - 1.27*logF_cent;
	double C0 = -0.4 - 0.67*logF_cent;
	double logP_r = log10(k_o*M / k_inf);
	double F = pow(10,(logF_cent / (1 + pow(((logP_r + C0) / (N0 - 0.14*(logP_r + C0))), 2))));
	r_c *= F;
	double R3b3 = r_c* gp[gp.H];

	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_f = 0; // radical fraction
	double r_f2 = 0;
	//double r_f3 = 0;
	//double r_f4 = 0;
	//double r_f5 = 0;

	if (T < 600) {
		r_f = 0;
		r_f2 = 0;
	}
	else{
		double Anrad = (m_r[0] + m_r[2]) / (m_r[1] + m_r[3] + R3b3);
		double AnC2H = (m_r[13]) / (m_r[14] + m_r[19]);
		double AnC2H_2 = (Anrad*m_r[7] + AnC2H*m_r[11]) / (m_r[16] + m_r[8] + m_r[17]);
		r_f = m_r[11] * Anrad + m_r[19] * AnC2H + m_r[17] * AnC2H_2;

		boost::numeric::ublas::matrix<double> arr1(5, 5);
		boost::numeric::ublas::vector<double> arr2(5);
		for (unsigned k = 0; k < arr1.size1(); ++k)
			for (unsigned l = 0; l < arr1.size2(); ++l)
				arr1(k, l) = 0.0;
		for (unsigned k = 0; k < arr2.size(); ++k)
			arr2(k) = 0.0;

		arr1(0, 0) = m_r[1] + m_r[3] + R3b3 + m_r[7] + m_r[9] + m_r[11];
		arr1(0, 1) = -m_r[10];
		arr1(0, 2) = -m_r[8];
		arr1(0, 3) = -m_r[12];
		arr1(1, 0) = -m_r[9];
		arr1(1, 1) = m_r[10] + m_r[14] + m_r[15] + m_r[19];
		arr1(1, 2) = -m_r[16];
		arr1(1, 3) = -m_r[20];
		arr1(2, 0) = -m_r[7];
		arr1(2, 1) = -m_r[15];
		arr1(2, 2) = m_r[8] + m_r[16] + m_r[17];
		arr1(2, 3) = -m_r[18];
		arr1(3, 0) = -m_r[11];
		arr1(3, 1) = -m_r[19];
		arr1(3, 2) = -m_r[17];
		arr1(3, 3) = m_r[12] + m_r[18] + m_r[20] + m_r[21] + m_r[23] + m_r[26];
		arr1(3, 4) = -(m_r[22] + m_r[24] + m_r[25]);
		arr1(4, 3) = -(m_r[21] + m_r[23] + m_r[26]);
		arr1(4, 4) = m_r[22] + m_r[24] + m_r[25];

		arr2(0) = m_r[0] + m_r[2] + m_r[6];
		arr2(1) = m_r[13];

		//std::cout << arr1 << std::endl;
		//std::cout << arr2 << std::endl;

		permutation_matrix<size_t> pm(arr1.size1());
		//std::cout << arr1 << std::endl;
		//std::cout << arr2 << std::endl;
		lu_factorize(arr1, pm);

		//std::cout << arr1 << std::endl;
		//std::cout << arr2 << std::endl;

		lu_substitute(arr1, pm, arr2);

		//std::cout << arr1 << std::endl;
		//std::cout << arr2 << std::endl;;
		r_f2 = arr2(0)*m_r[11] + arr2(1)*m_r[19] + arr2(2)*m_r[17]; //Lumped production of A4
		//r_f3 = arr2(0)*m_r[11] + arr2(1)*m_r[19] + arr2(2)*m_r[17] + arr2(4)*(m_r[22] + m_r[24] + m_r[25]); //Production of A4
		//r_f4 = arr2(3)*(m_r[12] + m_r[20] + m_r[18]); //Lumped consumption of A4
		//r_f5 = arr2(3)*(m_r[12] + m_r[20] + m_r[18] + m_r[21] + m_r[23] + m_r[26]); //Consumption of A4
		double A3rad = arr2(0);
		double A3C2H = arr2(1);
		double A3C2H2 = arr2(2);
		double A4 = arr2(3);
		double A4rad = arr2(4);
	}

	return m_rate = r_f2;*/

    /*return setRate0p0267(gp, pah_st);*/
}
// 
// ************************************************************
// ID2- R6 growth on FE (AR2 on Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void G6R_FE::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(2.50e14, 0, 16.00, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.40e9, 0.88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.10e13, 0, 4.569378, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(8.02e19, -2.011, 1.968, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(7.50e26, -3.96, 17.10, sp::C2H2));   //5 - r4f1
    addReaction(rxnV, Reaction(4.40e49, -11.6, 19.30, sp::C2H2));   //6 - r4f2 (r4f = r4f1+r4f2)
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //7 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.40e9, 0.88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.10e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(9.1e24, -3.39, 20.4, sp::C2H2));   //5 - r4f1
    addReaction(rxnV2, Reaction(9.9e41, -9.26, 15.7, sp::C2H2));   //6 - r4f2 (r4f = r4f1+r4f2)
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //7 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//Reactions from Narayanaswamy2010 and Mebel2016 for 1 atm.
	//rxns of A1-
	addReaction(rxnV3, Reaction(1.910e+61, -1.460e+01, 2.861e+01, sp::C2H2));                // A1- + C2H2 <=> A1C2H2                - 0             - Forward
	addReaction(rxnV3, Reaction(2.499e+69, -1.643e+01, 7.129e+01, sp::None));                // A1- + C2H2 <=> A1C2H2                - 1             - Backward
	addReaction(rxnV3, Reaction(1.100e+31, -4.830e+00, 2.662e+01, sp::C2H2));                // A1- + C2H2 <=> A1C2H + H             - 2             - Forward
	addReaction(rxnV3, Reaction(2.542e+37, -6.213e+00, 3.761e+01, sp::H));           		 // A1- + C2H2 <=> A1C2H + H             - 3             - Backward
	addReaction(rxnV3, Reaction(1.360e+75, -1.840e+01, 4.088e+01, sp::C2H2));                // A1- + C2H2 <=> A1C2H3*               - 4             - Forward
	addReaction(rxnV3, Reaction(4.055e+82, -2.012e+01, 7.940e+01, sp::None));                // A1- + C2H2 <=> A1C2H3*               - 5             - Backward
	addReaction(rxnV3, Reaction(6.000e+12, 0.000e+00, 0.000e+00, sp::C2H3));                 // A1- + C2H3 <=> A1C2H3                - 6             - Forward
	addReaction(rxnV3, Reaction(8.216e+23, -2.162e+00, 1.191e+02, sp::None));                // A1- + C2H3 <=> A1C2H3                - 7             - Backward
	addReaction(rxnV3, Reaction(9.450e-03, 4.470e+00, 4.472e+00, sp::C2H4));                 // A1- + C2H4 <=> A1 + C2H3             - 8             - Forward
	addReaction(rxnV3, Reaction(2.316e-02, 4.416e+00, 6.709e+00, sp::C2H3));                 // A1- + C2H4 <=> A1 + C2H3             - 9             - Backward
	addReaction(rxnV3, Reaction(1.260e+04, 2.610e+00, 1.434e+00, sp::C4H4));                 // A1- + C4H4 <=> A2 + H                - 10            - Forward
	addReaction(rxnV3, Reaction(1.130e+16, 7.538e-01, 6.694e+01, sp::H));            		 // A1- + C4H4 <=> A2 + H                - 11            - Backward
	addReaction(rxnV3, Reaction(4.570e+08, 1.880e+00, 1.484e+01, sp::H));            		 // A1 + H <=> A1- + H2          		 - 12            - Forward
	addReaction(rxnV3, Reaction(1.985e+04, 2.594e+00, 4.292e+00, sp::H2));           		 // A1 + H <=> A1- + H2          		 - 13            - Backward
	addReaction(rxnV3, Reaction(7.665e+04, 2.105e+00, 9.394e+00, sp::OH));           		 // A1 + OH <=> A1- + H2O                - 14            - Forward
	addReaction(rxnV3, Reaction(8.624e+01, 2.606e+00, 1.411e+01, sp::H2O));          		 // A1 + OH <=> A1- + H2O                - 15            - Backward
	addReaction(rxnV3, Reaction(4.170e+13, 1.500e-01, 0.000e+00, sp::H));            		 // A1- + H <=> A1               		 - 16            - Forward
	addReaction(rxnV3, Reaction(5.230e+17, -3.204e-01, 1.142e+02, sp::None));          		 // A1- + H <=> A1               		 - 17            - Backward
	//rxns of A1C2H3
	addReaction(rxnV3, Reaction(1.870e+07, 1.470e+00, 5.533e+00, sp::C2H3));                 // A1 + C2H3 <=> A1C2H3 + H             - 18            - Forward
	addReaction(rxnV3, Reaction(2.042e+14, -2.214e-01, 1.041e+01, sp::H));           		 // A1 + C2H3 <=> A1C2H3 + H             - 19            - Backward
	addReaction(rxnV3, Reaction(4.170e+13, 1.500e-01, 0.000e+00, sp::H));            		 // A1C2H3* + H <=> A1C2H3               - 20            - Forward
	addReaction(rxnV3, Reaction(2.523e+17, -3.422e-01, 1.161e+02, sp::None));        		 // A1C2H3* + H <=> A1C2H3               - 21            - Backward
	addReaction(rxnV3, Reaction(4.570e+08, 1.880e+00, 1.484e+01, sp::H));            		 // A1C2H3 + H <=> A1C2H3* + H2          - 22            - Forward
	addReaction(rxnV3, Reaction(4.115e+04, 2.615e+00, 2.435e+00, sp::H2));           		 // A1C2H3 + H <=> A1C2H3* + H2          - 23            - Backward
	addReaction(rxnV3, Reaction(7.665e+04, 2.105e+00, 9.394e+00, sp::OH));           		 // A1C2H3 + OH <=> A1C2H3* + H2O        - 24            - Forward
	addReaction(rxnV3, Reaction(1.788e+02, 2.628e+00, 1.226e+01, sp::H2O));          		 // A1C2H3 + OH <=> A1C2H3* + H2O        - 25            - Backward
	addReaction(rxnV3, Reaction(3.010e+14, 3.400e-01, 1.113e+02, sp::None));         		 // A1C2H3 <=> A1C2H2 + H                - 26            - Forward
	addReaction(rxnV3, Reaction(2.184e+11, 7.219e-01, -6.994e-01, sp::H));           		 // A1C2H3 <=> A1C2H2 + H                - 27            - Backward
	addReaction(rxnV3, Reaction(6.350e+04, 2.750e+00, 1.165e+01, sp::H));            		 // A1C2H3 + H <=> A1C2H2 + H2           - 28            - Forward
	addReaction(rxnV3, Reaction(2.509e+01, 3.375e+00, 3.404e+00, sp::H2));           		 // A1C2H3 + H <=> A1C2H2 + H2           - 29            - Backward
	addReaction(rxnV3, Reaction(6.550e-02, 4.200e+00, -8.604e-01, sp::OH));          		 // A1C2H3 + OH <=> A1C2H2 + H2O         - 30            - Forward
	addReaction(rxnV3, Reaction(6.705e-04, 4.613e+00, 6.162e+00, sp::H2O));          		 // A1C2H3 + OH <=> A1C2H2 + H2O         - 31            - Backward
	//rxns of A1C2H3*
	addReaction(rxnV3, Reaction(2.440e+30, -5.730e+00, 3.207e+01, sp::None));                // A1C2H2 <=> A1C2H3*           		 - 32            - Forward
	addReaction(rxnV3, Reaction(5.560e+29, -5.620e+00, 2.791e+01, sp::None));                // A1C2H2 <=> A1C2H3*           		 - 33            - Backward
	addReaction(rxnV3, Reaction(3.020e+10, 7.020e-01, 5.530e+00, sp::C2H2));                 // A1C2H3* + C2H2 <=> A2 + H    		 - 34            - Forward
	addReaction(rxnV3, Reaction(1.387e+21, -7.978e-01, 7.245e+01, sp::H));           		 // A1C2H3* + C2H2 <=> A2 + H    		 - 35            - Backward
	//rxns of A1C2H
	addReaction(rxnV3, Reaction(4.170e+13, 1.500e-01, 0.000e+00, sp::H));            		 // A1C2H* + H <=> A1C2H                 - 36            - Forward
	addReaction(rxnV3, Reaction(2.075e+17, -3.346e-01, 1.168e+02, sp::None));                // A1C2H* + H <=> A1C2H                 - 37            - Backward
	addReaction(rxnV3, Reaction(4.570e+08, 1.880e+00, 1.484e+01, sp::H));            		 // A1C2H + H <=> A1C2H* + H2            - 38            - Forward
	addReaction(rxnV3, Reaction(5.003e+04, 2.608e+00, 1.761e+00, sp::H2));           		 // A1C2H + H <=> A1C2H* + H2            - 39            - Backward
	addReaction(rxnV3, Reaction(7.665e+04, 2.105e+00, 9.394e+00, sp::OH));           		 // A1C2H + OH <=> A1C2H* + H2O          - 40            - Forward
	addReaction(rxnV3, Reaction(2.174e+02, 2.621e+00, 1.158e+01, sp::H2O));          		 // A1C2H + OH <=> A1C2H* + H2O          - 41            - Backward
	addReaction(rxnV3, Reaction(1.590e+62, -1.450e+01, 3.176e+01, sp::H));           		 // A1C2H + H <=> A1C2H2                 - 42            - Forward
	addReaction(rxnV3, Reaction(9.003e+63, -1.495e+01, 6.344e+01, sp::None));        		 // A1C2H + H <=> A1C2H2                 - 43            - Backward
	addReaction(rxnV3, Reaction(1.650e+11, 4.900e-01, 1.063e+01, sp::H));            		 // A1C2H2 + H <=> A1C2H + H2            - 44            - Forward
	addReaction(rxnV3, Reaction(1.587e+09, 1.184e+00, 8.265e+01, sp::H2));           		 // A1C2H2 + H <=> A1C2H + H2            - 45            - Backward
	addReaction(rxnV3, Reaction(2.500e+12, 0.000e+00, 0.000e+00, sp::OH));           		 // A1C2H2 + OH <=> A1C2H + H2O          - 46            - Forward
	addReaction(rxnV3, Reaction(6.230e+11, 4.815e-01, 8.728e+01, sp::H2O));          		 // A1C2H2 + OH <=> A1C2H + H2O          - 47            - Backward
	addReaction(rxnV3, Reaction(3.600e+17, -1.440e+00, 1.576e+01, sp::C2H3));        		 // A1C2H + C2H3 <=> A2 + H              - 48            - Forward
	addReaction(rxnV3, Reaction(1.619e+29, -3.226e+00, 7.470e+01, sp::H));           		 // A1C2H + C2H3 <=> A2 + H              - 49            - Backward
	//rxns of A1C2H*
	addReaction(rxnV3, Reaction(4.490e+82, -2.000e+01, 5.183e+01, sp::C2H2));                // A1C2H* + C2H2 <=> A2-                - 50            - Forward
	addReaction(rxnV3, Reaction(1.338e+94, -2.184e+01, 1.435e+02, sp::None));                // A1C2H* + C2H2 <=> A2-                - 51            - Backward
	addReaction(rxnV3, Reaction(1.180e+104, -2.570e+01, 7.682e+01, sp::C2H2));               // A1C2H* + C2H2 <=> A2*                - 52            - Forward
	addReaction(rxnV3, Reaction(3.917e+115, -2.755e+01, 1.688e+02, sp::None));               // A1C2H* + C2H2 <=> A2*                - 53            - Backward
	addReaction(rxnV3, Reaction(1.760e+40, -7.040e+00, 4.821e+01, sp::C2H2));                // A1C2H* + C2H2 <=> A1C2HX2 + H        - 54            - Forward
	addReaction(rxnV3, Reaction(9.718e+46, -8.438e+00, 6.084e+01, sp::H));           		 // A1C2H* + C2H2 <=> A1C2HX2 + H        - 55            - Backward
	addReaction(rxnV3, Reaction(3.620e+28, -4.240e+00, 2.386e+01, sp::C2H4));                // A1C2H* + C2H4 <=> A2 + H             - 56            - Forward
	addReaction(rxnV3, Reaction(1.583e+40, -6.094e+00, 8.758e+01, sp::H));           		 // A1C2H* + C2H4 <=> A2 + H             - 57            - Backward
	//rxns of A1C2H2
	addReaction(rxnV3, Reaction(3.570e+22, -2.720e+00, 1.447e+01, sp::C2H2));                // A1C2H2 + C2H2 <=> A2 + H             - 58            - Forward
	addReaction(rxnV3, Reaction(3.736e+32, -4.109e+00, 7.723e+01, sp::H));           		 // A1C2H2 + C2H2 <=> A2 + H             - 59            - Backward
	//rxns of A1C2HX2
	addReaction(rxnV3, Reaction(1.010e+86, -2.060e+01, 5.670e+01, sp::H));           		 // A1C2HX2 + H <=> A2-          		 - 60            - Forward
	addReaction(rxnV3, Reaction(5.450e+90, -2.104e+01, 1.358e+02, sp::None));                // A1C2HX2 + H <=> A2-          		 - 61            - Backward
	addReaction(rxnV3, Reaction(6.000e+108, -2.660e+01, 8.359e+01, sp::H));          		 // A1C2HX2 + H <=> A2*          		 - 62            - Forward
	addReaction(rxnV3, Reaction(3.607e+113, -2.705e+01, 1.629e+02, sp::None));               // A1C2HX2 + H <=> A2*          		 - 63            - Backward
	//rxns of A2-
	addReaction(rxnV3, Reaction(4.170e+13, 1.500e-01, 0.000e+00, sp::H));            		 // A2- + H <=> A2               		 - 64            - Forward
	addReaction(rxnV3, Reaction(4.126e+17, -3.390e-01, 1.195e+02, sp::None));                // A2- + H <=> A2               		 - 65            - Backward
	addReaction(rxnV3, Reaction(4.570e+08, 1.880e+00, 1.484e+01, sp::H));            		 // A2 + H <=> A2- + H2          		 - 66            - Forward
	addReaction(rxnV3, Reaction(2.516e+04, 2.612e+00, -9.970e-01, sp::H2));          		 // A2 + H <=> A2- + H2          		 - 67            - Backward
	addReaction(rxnV3, Reaction(7.665e+04, 2.105e+00, 9.394e+00, sp::OH));           		 // A2 + OH <=> A2- + H2O        		 - 68            - Forward
	addReaction(rxnV3, Reaction(1.093e+02, 2.625e+00, 8.825e+00, sp::H2O));          		 // A2 + OH <=> A2- + H2O        		 - 69            - Backward
	//rxns of A2*
	addReaction(rxnV3, Reaction(4.170e+13, 1.500e-01, 0.000e+00, sp::H));            		 // A2* + H <=> A2               		 - 70            - Forward
	addReaction(rxnV3, Reaction(3.703e+17, -3.255e-01, 1.192e+02, sp::None));                // A2* + H <=> A2               		 - 71            - Backward
	addReaction(rxnV3, Reaction(4.570e+08, 1.880e+00, 1.484e+01, sp::H));            		 // A2 + H <=> A2* + H2          		 - 72            - Forward
	addReaction(rxnV3, Reaction(2.803e+04, 2.599e+00, -6.926e-01, sp::H2));          		 // A2 + H <=> A2* + H2          		 - 73            - Backward
	addReaction(rxnV3, Reaction(7.665e+04, 2.105e+00, 9.394e+00, sp::OH));           		 // A2 + OH <=> A2* + H2O                - 74            - Forward
	addReaction(rxnV3, Reaction(1.218e+02, 2.611e+00, 9.129e+00, sp::H2O));          		 // A2 + OH <=> A2* + H2O                - 75            - Backward
	//rxns of A2

	//Reactions for the ABF mechanism
	/*//Reactions of A1 -
	addReaction(rxnV3, Reaction(2.510e+12, 0.000e+00, 6.190e+00, sp::C2H4));			// A1- + C2H4 <=> A1C2H3 + H            - 0             - Forward
	addReaction(rxnV3, Reaction(1.190e+18, -1.324e+00, 1.005e+01, sp::H));				// A1- + C2H4 <=> A1C2H3 + H            - 1             - Backward
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));				// A1 + H <=> A1- + H2					- 2             - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));				// A1 + H <=> A1- + H2					- 3             - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));               // A1- + H -> A1              			- 4              - Forward
	//addReaction(rxnV3, Reaction(1.000e+14, 0.000e+00, 0.000e+00, sp::None));            // A1- + H (+M) <=> A1 (+M)			- 4
	//addReaction(rxnV3, Reaction(6.600e+75, -1.630e+01, 7.000e+00, sp::None));           // A1- + H (+M) <=> A1 (+M)			- 5
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));              // A1 + OH <=> A1- + H2O                - 6             - Forward
	addReaction(rxnV3, Reaction(7.665e+04, 2.105e+00, 9.394e+00, sp::H2O));             // A1 + OH <=> A1- + H2O                - 7             - Backward
	addReaction(rxnV3, Reaction(3.300e+33, -5.700e+00, 2.550e+01, sp::C2H2));           // A1- + C2H2 <=> A1C2H + H             - 8             - Forward
	addReaction(rxnV3, Reaction(5.077e+41, -7.589e+00, 3.365e+01, sp::H));              // A1- + C2H2 <=> A1C2H + H             - 9             - Backward
	addReaction(rxnV3, Reaction(1.100e+07, 1.710e+00, 3.894e+00, sp::C2H2));           // A1- + C2H2 <=> n-A1C2H2              - 10            - Forward
	addReaction(rxnV3, Reaction(1.300e+04, 0.000e+00, 4.1778e+01, sp::None));           // A1- + C2H2 <=> n-A1C2H2              - 11            - Backward
	addReaction(rxnV3, Reaction(1.200e+27, -4.220e+00, 7.235e+00, sp::C2H3));           // A1- + C2H3 <=> A1C2H3                - 12            - Forward
	addReaction(rxnV3, Reaction(1.017e+37, -6.132e+00, 1.231e+02, sp::None));           // A1- + C2H3 <=> A1C2H3                - 13            - Backward
	//Reactions of n - A1C2H2
	addReaction(rxnV3, Reaction(3.100e+06, 2.000e+00, 3.430e+00, sp::OH));              // A1C2H3 + OH <=> H2O + n-A1C2H2       - 14            - Forward
	addReaction(rxnV3, Reaction(5.677e+06, 1.711e+00, 1.208e+01, sp::H2O));             // A1C2H3 + OH <=> H2O + n-A1C2H2       - 15            - Backward
	addReaction(rxnV3, Reaction(6.650e+06, 2.530e+00, 1.224e+01, sp::H));               // A1C2H3 + H <=> H2 + n-A1C2H2         - 16            - Forward
	addReaction(rxnV3, Reaction(3.776e+05, 2.486e+00, 5.248e+00, sp::H2));              // A1C2H3 + H <=> H2 + n-A1C2H2         - 17            - Backward
	addReaction(rxnV3, Reaction(1.500e+10, 0.850e+00, 1.194e+01, sp::H));              // A1C2H + H <=> n-A1C2H2               - 18            - Forward
	addReaction(rxnV3, Reaction(6.800e+11, 0.000e+00, 2.2025e+01, sp::None));           // A1C2H + H <=> n-A1C2H2               - 19            - Backward
	//Reactions of A1C2H
	addReaction(rxnV3, Reaction(5.000e+13, 0.000e+00, 0.000e+00, sp::C2H));             // A1C2H + C2H <=> A1C2H)2 + H          - 20            - Forward
	addReaction(rxnV3, Reaction(1.105e+22, -1.680e+00, 2.918e+01, sp::H));              // A1C2H + C2H <=> A1C2H)2 + H          - 21            - Backward
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));               // A1C2H + H <=> A1C2H* + H2            - 22            - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));              // A1C2H + H <=> A1C2H* + H2            - 23            - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));              // A1C2H + OH <=> A1C2H* + H2O          - 24            - Forward
	addReaction(rxnV3, Reaction(1.359e+06, 1.948e+00, 9.131e+00, sp::H2O));             // A1C2H + OH <=> A1C2H* + H2O          - 25            - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));               // A1C2H* + H -> A1C2H         			- 26              - Forward
	//addReaction(rxnV3, Reaction(1.000e+14, 0.000e+00, 0.000e+00, sp::None));            // A1C2H* + H (+M) <=> A1C2H (+M)		- 26
	//addReaction(rxnV3, Reaction(6.600e+75, -1.630e+01, 7.000e+00, sp::None));           // A1C2H* + H (+M) <=> A1C2H (+M)		- 27
	addReaction(rxnV3, Reaction(5.000e+13, 0.000e+00, 0.000e+00, sp::C2H));             // A1 + C2H <=> A1C2H + H               - 28            - Forward
	addReaction(rxnV3, Reaction(3.885e+20, -1.578e+00, 2.991e+01, sp::H));              // A1 + C2H <=> A1C2H + H               - 29            - Backward
	//Reactions of A1C2H3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));              // A1C2H3 + OH <=> A1C2H3* + H2O        - 30            - Forward
	addReaction(rxnV3, Reaction(1.336e+06, 1.973e+00, 1.121e+01, sp::H2O));             // A1C2H3 + OH <=> A1C2H3* + H2O        - 31            - Backward
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));               // A1C2H3 + H <=> A1C2H3* + H2          - 32            - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));              // A1C2H3 + H <=> A1C2H3* + H2          - 33            - Backward
	addReaction(rxnV3, Reaction(7.900e+11, 0.000e+00, 6.400e+00, sp::C2H3));            // A1 + C2H3 <=> A1C2H3 + H             - 34            - Forward
	addReaction(rxnV3, Reaction(4.561e+16, -1.046e+00, 1.067e+01, sp::H));              // A1 + C2H3 <=> A1C2H3 + H             - 35            - Backward
	//Reactions of A1C2H3*
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));                // A1C2H3* + H -> A1C2H3              	- 4              - Forward
	//addReaction(rxnV3, Reaction(1.000e+14, 0.000e+00, 0.000e+00, sp::None));            // A1C2H3* + H (+M) <=> A1C2H3 (+M)		- 36
	//addReaction(rxnV3, Reaction(6.600e+75, -1.630e+01, 7.000e+00, sp::None));           // A1C2H3* + H (+M) <=> A1C2H3 (+M)		- 37
	//Reactions of A1C2H*
	addReaction(rxnV3, Reaction(1.800e+19, -1.670e+00, 1.880e+01, sp::C2H2));           // A1C2H* + C2H2 <=> A1C2H)2 + H        - 38            - Forward
	addReaction(rxnV3, Reaction(4.444e+27, -3.504e+00, 2.648e+01, sp::H));              // A1C2H* + C2H2 <=> A1C2H)2 + H        - 39            - Backward
	addReaction(rxnV3, Reaction(2.200e+62, -1.456e+01, 3.310e+01, sp::C2H2));           // A1C2H* + C2H2 <=> A2-1               - 40            - Forward
	addReaction(rxnV3, Reaction(3.933e+73, -1.646e+01, 1.288e+02, sp::None));           // A1C2H* + C2H2 <=> A2-1               - 41            - Backward
	//Reactions of A1C2H)2
	addReaction(rxnV3, Reaction(1.400e+64, -1.457e+01, 2.990e+01, sp::H));              // A1C2H)2 + H <=> A2-1                 - 42            - Forward
	addReaction(rxnV3, Reaction(1.014e+67, -1.463e+01, 1.179e+02, sp::None));           // A1C2H)2 + H <=> A2-1                 - 43            - Backward
	//Reactions of A2 production
	addReaction(rxnV3, Reaction(1.600e+16, -1.330e+00, 6.660e+00, sp::C2H2));           // A1C2H3* + C2H2 <=> A2 + H            - 44            - Forward
	addReaction(rxnV3, Reaction(1.600e+16, -1.330e+00, 5.4000e+00, sp::C2H2));          // n-A1C2H2 + C2H2 <=> A2 + H           - 45            - Forward
	*/
	
    /*addReaction(rxnV3, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));   //1 - r1b
    addReaction(rxnV3, Reaction(1.0e10, 0.734, 1.43, sp::OH));  //2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV3, Reaction(2.0e13, 0, 0, sp::H));  //4 - r3f
    addReaction(rxnV3, Reaction(1.1e7, 1.71, 3.9, sp::C2H2));   //5 - r4f1
    addReaction(rxnV3, Reaction(3.3e33, -5.7, 25.5, sp::C2H2));   //6 - r4f2 (r4f = r4f1+r4f2)
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //7 - r5f*/
	
    m_sType = FE; // sitetype
    m_name = "G6R at FE"; // name of process
    m_ID = 2;
}
// Jump rate calculation
double G6R_FE::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = (m_r[5]+m_r[6])*r_f* site_count; // Rate Equation
}
double G6R_FE::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double G6R_FE::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {

	// check if site count is zero
	double T = gp[gp.T];
	double r_f;
	double r_f2_ss, r_f2_peq;
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	//int ring_count = ((int)pah_st.getRingsCount().first); // Ring count
	//if (ring_count > 2) return m_rate = 0.0;
	if (site_count == 0) return m_rate = 0;
	
	/*//For ABF mechanism with 3-body reactions
	// Calculating 3 body reactions
	//Parameters alpha = 0.2, T3 = 122.8, T1 = 478.4, T2 = 5411.9
	//Parameters = [1.00000000e+00   1.00000000e-01   5.84900000e+02   6.11300000e+03]
	double P = gp[gp.P]; //Pressure in Pa
	double M_eff1 = P / 8.314E+06 / T;
	double M = M_eff1 + (1 * gp[gp.H2]) + (5 * gp[gp.H2O]) + (1 * gp[gp.CH4]) + (0.5*gp[gp.CO]) + (1 * gp[gp.CO2]) + (2 * gp[gp.C2H6]);
	double k_o1 = m_r[5];
	double k_inf1 = m_r[4];
	double F_cent1 = 0.00 * exp(-T / 0.1) + 1.00*exp(-T / 584.9) + exp(-6113 / T);
	double logF_cent1 = log10(F_cent1);
	double r_c1 = k_inf1 / (1 + (k_inf1 / (k_o1*M)));
	double N01 = 0.75 - 1.27*logF_cent1;
	double C01 = -0.4 - 0.67*logF_cent1;
	double logP_r1 = log10(k_o1*M / k_inf1);
	double F1 = pow(10,(logF_cent1 / (1 + pow(((logP_r1 + C01) / (N01 - 0.14*(logP_r1 + C01))), 2))));
	r_c1 *= F1;
	double R3b1 = r_c1* gp[gp.H];
	//Parameters = [1.00000000e+00   1.00000000e-01   5.84900000e+02   6.11300000e+03]
	double k_o2 = m_r[27];
	double k_inf2 = m_r[26];
	double F_cent2 = 0.00 * exp(-T / 0.1) + 1.00*exp(-T / 584.9) + exp(-6113 / T);
	double logF_cent2 = log10(F_cent2);
	double r_c2 = k_inf2 / (1 + (k_inf2 / (k_o2*M)));
	double N02 = 0.75 - 1.27*logF_cent2;
	double C02 = -0.4 - 0.67*logF_cent2;
	double logP_r2 = log10(k_o2*M / k_inf2);
	double F2 = pow(10,(logF_cent2 / (1 + pow(((logP_r2 + C02) / (N02 - 0.14*(logP_r2 + C02))), 2))));
	r_c2 *= F2;
	double R3b2 = r_c2* gp[gp.H];
	//Parameters = [1.00000000e+00   1.00000000e-01   5.84900000e+02   6.11300000e+03]
	double k_o3 = m_r[37];
	double k_inf3 = m_r[36];
	double F_cent3 = 0.00 * exp(-T / 0.1) + 1.00*exp(-T / 584.9) + exp(-6113 / T);
	double logF_cent3 = log10(F_cent3);
	double r_c3 = k_inf3 / (1 + (k_inf3 / (k_o3*M)));
	double N03 = 0.75 - 1.27*logF_cent3;
	double C03 = -0.4 - 0.67*logF_cent3;
	double logP_r3 = log10(k_o3*M / k_inf3);
	double F3 = pow(10,(logF_cent3 / (1 + pow(((logP_r3 + C03) / (N03 - 0.14*(logP_r3 + C03))), 2))));
	r_c3 *= F3;
	double R3b3 = r_c3* gp[gp.H];

	//int ring_count = ((int)pah_st.getRingsCount().first); // Ring count
	//if (ring_count > 2) return m_rate = 0.0;
	//if (site_count == 0) return m_rate = 0;
	//double A1_conc = gp[gp.C6H6];
	//double C2H2_conc = gp[gp.C2H2];*/

	try{
		//Reactions from Mebel2016 and Narayanaswamy2010 and Blanquart2009
		matrix<double> arr1(10, 10);
		boost::numeric::ublas::vector<double> arr2(10);
		matrix<double> arr1_peq(7, 7);
		boost::numeric::ublas::vector<double> arr2_peq(7);
		for (unsigned k = 0; k < arr1.size1(); ++k)
			for (unsigned l = 0; l < arr1.size2(); ++l)
				arr1(k, l) = 0.0;
		for (unsigned k = 0; k < arr2.size(); ++k)
			arr2(k) = 0.0;
		for (unsigned k = 0; k < arr1_peq.size1(); ++k)
			for (unsigned l = 0; l < arr1_peq.size2(); ++l)
				arr1_peq(k, l) = 0.0;
		for (unsigned k = 0; k < arr2_peq.size(); ++k)
			arr2_peq(k) = 0.0;

		//Species_list_SS = ['A1-','A1C2H3','A1C2H3*','A1C2H','A1C2H*','A1C2H2','A1C2HX2','A2-','A2*','A2']
		arr1(0, 0) =  (m_r[0] + m_r[2] + m_r[4] + m_r[6] + m_r[8] + m_r[10] + m_r[13] + m_r[15] + m_r[16]);
		arr1(0, 1) = -(m_r[7]);
		arr1(0, 2) = -(m_r[5]);
		arr1(0, 3) = -(m_r[3]);
		arr1(0, 5) = -(m_r[1]);
		arr1(0, 9) = -(m_r[11]);
		arr1(1, 0) = -(m_r[6]);
		arr1(1, 1) =  (m_r[7] + m_r[19] + m_r[21] + m_r[22] + m_r[24] + m_r[26] + m_r[28] + m_r[30]);
		arr1(1, 2) = -(m_r[20] + m_r[23] + m_r[25]);
		arr1(1, 5) = -(m_r[27] + m_r[29] + m_r[31]);
		arr1(2, 0) = -(m_r[4]);
		arr1(2, 1) = -(m_r[21] + m_r[22] + m_r[24]);
		arr1(2, 2) =  (m_r[5] + m_r[33] + m_r[20] + m_r[23] + m_r[25] + m_r[34]);
		arr1(2, 5) = -(m_r[32]);
		arr1(2, 9) = -(m_r[35]);
		arr1(3, 0) = -(m_r[2]);
		arr1(3, 3) =  (m_r[3] + m_r[37] + m_r[38] + m_r[40] + m_r[42] + m_r[45] + m_r[47] + m_r[48]);
		arr1(3, 4) = -(m_r[36] + m_r[39] + m_r[41]);
		arr1(3, 5) = -(m_r[43] + m_r[44] + m_r[46]);
		arr1(3, 9) = -(m_r[49]);
		arr1(4, 3) = -(m_r[37] + m_r[38] + m_r[40]);
		arr1(4, 4) =  (m_r[36] + m_r[39] + m_r[41] + m_r[50] + m_r[52] + m_r[54] + m_r[56]);
		arr1(4, 6) = -(m_r[55]);
		arr1(4, 7) = -(m_r[51]);
		arr1(4, 8) = -(m_r[53]);
		arr1(4, 9) = -(m_r[57]);
		arr1(5, 0) = -(m_r[0]);
		arr1(5, 1) = -(m_r[26] + m_r[28] + m_r[30]);
		arr1(5, 2) = -(m_r[33]);
		arr1(5, 3) = -(m_r[42] + m_r[45] + m_r[47]);
		arr1(5, 5) =  (m_r[1] + m_r[32] + m_r[43] + m_r[44] + m_r[46] + m_r[27] + m_r[29] + m_r[31] + m_r[58]);
		arr1(5, 9) = -(m_r[59]);
		arr1(6, 4) = -(m_r[54]);
		arr1(6, 6) =  (m_r[55] + m_r[60] + m_r[62]);
		arr1(6, 7) = -(m_r[61]);
		arr1(6, 8) = -(m_r[63]);
		arr1(7, 4) = -(m_r[50]);
		arr1(7, 6) = -(m_r[60]);
		arr1(7, 7) =  (m_r[51] + m_r[61] + m_r[64] + m_r[67] + m_r[69]);
		arr1(7, 9) = -(m_r[65] + m_r[66] + m_r[68]);
		arr1(8, 4) = -(m_r[52]);
		arr1(8, 6) = -(m_r[62]);
		arr1(8, 8) =  (m_r[53] + m_r[63] + m_r[70] + m_r[73] + m_r[75]);
		arr1(8, 9) = -(m_r[71] + m_r[72] + m_r[74]);
		arr1(9, 0) = -(m_r[10]);
		arr1(9, 2) = -(m_r[34]);
		arr1(9, 3) = -(m_r[48]);
		arr1(9, 4) = -(m_r[56]);
		arr1(9, 5) = -(m_r[58]);
		arr1(9, 7) = -(m_r[64] + m_r[67] + m_r[69]);
		arr1(9, 8) = -(m_r[70] + m_r[73] + m_r[75]);
		arr1(9, 9) =  (m_r[35] + m_r[59] + m_r[49] + m_r[57] + m_r[11] + m_r[65] + m_r[66] + m_r[68] + m_r[71] + m_r[72] + m_r[74]);

		arr2(0) =  (m_r[9] + m_r[12] + m_r[14] + m_r[17]);
		arr2(1) =  (m_r[18]);
		
		//Species_list_PEQ = ['A1-','A1C2H3','A1C2H3*','A1C2H','A1C2H*','A1C2H2','A1C2HX2']
		arr1_peq(0, 0) =  (m_r[0] + m_r[2] + m_r[4] + m_r[6] + m_r[8] + m_r[10] + m_r[13] + m_r[15] + m_r[16]);
		arr1_peq(0, 1) = -(m_r[7]);
		arr1_peq(0, 2) = -(m_r[5]);
		arr1_peq(0, 3) = -(m_r[3]);
		arr1_peq(0, 5) = -(m_r[1]);
		arr1_peq(1, 0) = -(m_r[6]);
		arr1_peq(1, 1) =  (m_r[7] + m_r[19] + m_r[21] + m_r[22] + m_r[24] + m_r[26] + m_r[28] + m_r[30]);
		arr1_peq(1, 2) = -(m_r[20] + m_r[23] + m_r[25]);
		arr1_peq(1, 5) = -(m_r[27] + m_r[29] + m_r[31]);
		arr1_peq(2, 0) = -(m_r[4]);
		arr1_peq(2, 1) = -(m_r[21] + m_r[22] + m_r[24]);
		arr1_peq(2, 2) =  (m_r[5] + m_r[33] + m_r[20] + m_r[23] + m_r[25] + 1.0*m_r[34]);
		arr1_peq(2, 5) = -(m_r[32]);
		arr1_peq(3, 0) = -(m_r[2]);
		arr1_peq(3, 3) =  (m_r[3] + m_r[37] + m_r[38] + m_r[40] + m_r[42] + m_r[45] + m_r[47] + 1.0*m_r[48]);
		arr1_peq(3, 4) = -(m_r[36] + m_r[39] + m_r[41]);
		arr1_peq(3, 5) = -(m_r[43] + m_r[44] + m_r[46]);
		arr1_peq(4, 3) = -(m_r[37] + m_r[38] + m_r[40]);
		arr1_peq(4, 4) =  (m_r[36] + m_r[39] + m_r[41] + 1.0*m_r[50] + 0.0*m_r[52] + m_r[54] + 1.0*m_r[56]);
		arr1_peq(4, 6) = -(m_r[55]);
		arr1_peq(5, 0) = -(m_r[0]);
		arr1_peq(5, 1) = -(m_r[26] + m_r[28] + m_r[30]);
		arr1_peq(5, 2) = -(m_r[33]);
		arr1_peq(5, 3) = -(m_r[42] + m_r[45] + m_r[47]);
		arr1_peq(5, 5) =  (m_r[1] + m_r[32] + m_r[43] + m_r[44] + m_r[46] + m_r[27] + m_r[29] + m_r[31] + 1.0*m_r[58]);
		arr1_peq(6, 4) = -(m_r[54]);
		arr1_peq(6, 6) =  (m_r[55] + m_r[60] + m_r[62]);

		arr2_peq(0) =  (m_r[9] + m_r[12] + m_r[14] + m_r[17]);
		arr2_peq(1) =  (m_r[18]);
		
		//std::cout << "Steady-state matrix:\n"; 
		//std::cout << arr1 << std::endl;
		//std::cout << arr2 << std::endl;
		
		//std::cout << "Partial-equilibrium matrix:\n";
		//std::cout << arr1_peq << std::endl;
		//std::cout << arr2_peq << std::endl;
		
		permutation_matrix<size_t> pm(arr1.size1());
		//std::cout << pm << std::endl;
		lu_factorize(arr1, pm);
		//std::cout << pm << std::endl;
		lu_substitute(arr1, pm, arr2);
		//std::cout << arr2 << std::endl;

		permutation_matrix<size_t> pm_peq(arr1_peq.size1());
		//std::cout << pm_peq << std::endl;
		lu_factorize(arr1_peq, pm_peq);
		//std::cout << pm_peq << std::endl;
		lu_substitute(arr1_peq, pm_peq, arr2_peq);
		//std::cout << arr2_peq << std::endl;
		
		r_f2_ss =  arr2(4)*m_r[50] + arr2(6)*m_r[60] + arr2(4)*m_r[52] + arr2(6)*m_r[62] + arr2(2)*m_r[34] + arr2(5)*m_r[58] + arr2(3)*m_r[48] + arr2(4)*m_r[56] + arr2(0)*m_r[10]; //Production of the sum of A2* and A2*/
		
		r_f2_peq = arr2_peq(4)*m_r[50] + arr2_peq(6)*m_r[60] + arr2_peq(4)*m_r[52] + arr2_peq(6)*m_r[62] + arr2_peq(2)*m_r[34] + arr2_peq(5)*m_r[58] + arr2_peq(3)*m_r[48] + arr2_peq(4)*m_r[56] + arr2_peq(0)*m_r[10]; //Production of the sum of A2* and A2*/
		
		//cout << "Steady-state rate = " << r_f2_ss;
		//cout << "Partial-equilibrium rate = " << r_f2_peq; 
		
		//For ABF REACTIONS:
		/*matrix<double> arr1(7, 7);
		boost::numeric::ublas::vector<double> arr2(7);
		for (unsigned k = 0; k < arr1.size1(); ++k)
			for (unsigned l = 0; l < arr1.size2(); ++l)
				arr1(k, l) = 0.0;
		for (unsigned k = 0; k < arr2.size(); ++k)
			arr2(k) = 0.0;

		//Speclist = ['A1-', 'n-A1C2H2', 'A1C2H', 'A1C2H3', 'A1C2H3*', 'A1C2H*','A1C2H)2']
		arr1(0, 0) = m_r[0] + m_r[3] + m_r[7] + m_r[8] + m_r[10] + m_r[12] + R3b1; //A1- consumption
		arr1(0, 1) = -(m_r[11]);
		arr1(0, 2) = -(m_r[9]);
		arr1(0, 3) = -(m_r[1] + m_r[13]);
		arr1(1, 0) = -(m_r[10]);
		arr1(1, 1) = m_r[11] + m_r[15] + m_r[17] + m_r[19];
		arr1(1, 2) = -(m_r[18]);
		arr1(1, 3) = -(m_r[14] + m_r[16]); 
		arr1(2, 0) = -(m_r[8]);
		arr1(2, 1) = -(m_r[19]);
		arr1(2, 2) = m_r[9] + m_r[18] + m_r[20] + m_r[22] + m_r[24] + m_r[29];
		arr1(2, 5) = -(m_r[23] + m_r[25] + R3b2);
		arr1(2, 6) = -(m_r[21]);
		arr1(3, 0) = -(m_r[0] + m_r[12]);
		arr1(3, 1) = -(m_r[15]+ m_r[17]);
		arr1(3, 3) = m_r[1] + m_r[13] + m_r[14] + m_r[16] + m_r[30] + m_r[32] + m_r[35];
		arr1(3, 4) = -(m_r[31] + m_r[33] + R3b3);
		arr1(4, 3) = -(m_r[30] + m_r[32]);
		arr1(4, 4) = m_r[31] + m_r[33] + R3b3;
		arr1(5, 2) = -(m_r[22] + m_r[24]);
		arr1(5, 5) = m_r[23] + m_r[25] + m_r[38] + m_r[40] + R3b2;
		arr1(5, 6) = -(m_r[39]);
		arr1(6, 2) = -(m_r[20]);
		arr1(6, 5) = -(m_r[38]);
		arr1(6, 6) = m_r[21] + m_r[39] + m_r[42];

		arr2(0) = m_r[2] + m_r[6];
		arr2(2) = m_r[28];
		arr2(3) = m_r[34];

		//std::cout << arr1 << std::endl;
		//std::cout << arr2 << std::endl;

		permutation_matrix<size_t> pm(arr1.size1());
		lu_factorize(arr1, pm);
		lu_substitute(arr1, pm, arr2);

		double A1rad = arr2(0);
		double nA1C2H2 = arr2(1);
		double A1C2H = arr2(2);
		double A1C2H3 = arr2(3);
		double A1C2H3rad = arr2(4);
		double A1C2Hrad = arr2(5);
		double A1C2H2 = arr2(6);

		r_f2 = arr2(5)*m_r[40] + arr2(6)*m_r[42] + arr2(1)*m_r[45] + arr2(4)*m_r[44];*/ //Production of the sum of A2* and A2
		//double r_f3 = arr2(8)*(m_r[44] + m_r[46]);
		//double A2radtoA1ratio = arr2(9);
		//double A2toA1ratio = arr2(8);

		//double test = 1.0;
	}
	catch ( const std::runtime_error& (s)){
		return m_rate = 0.0; 
	}

	double r_denom = (m_r[47]+m_r[49]+m_r[50]+m_r[51]);

	if(r_denom>0) {
	 r_f = (m_r[46]+m_r[48])/r_denom; 
	}
	else r_f=0;
	
	
	/*JumpProcess* j_D6R_FE3_contained = new D6R_FE3; j_D6R_FE3_contained->initialise();
	(j_D6R_FE3_contained)->calculateElemRxnRate((j_D6R_FE3_contained)->getVec1(), gp);
	double desorpt_rate = (j_D6R_FE3_contained)->setRate1(gp, pah_st);
	*/
	
	//This is the reversible rate. Checks how far away from steady state are 6 FE sites and 2 FE3 sites  to decide if PEQ or SS approximation should be valid.
	double FE3_site_count = site_count / 6.0 * 2.0;
	double desorpt_rate;
	double r_denom_desorpt = (m_r[64] + m_r[67] + m_r[69]);
	if (r_denom_desorpt <=0) desorpt_rate = 0.0;
	else{ 
		double r_num_desorpt = (m_r[65] + m_r[66] + m_r[68]);
		desorpt_rate = ((m_r[11] + m_r[35] + m_r[49] + m_r[57] + m_r[59]) + r_num_desorpt/r_denom_desorpt*(m_r[51] + m_r[53] + m_r[61] + m_r[63])) * FE3_site_count / 2.0; 
	}
	
	double m_rate_old = m_r[51]*r_f* site_count; // Rate Equation
	//cout << "Temp = " << T; 
	double lambda = 1.0;
	double peq_rate = lambda * r_f2_peq * site_count / 6.0;
	if (peq_rate > desorpt_rate) m_rate = peq_rate;
	else m_rate = r_f2_ss * site_count / 6.0; // Rate Equation
	
		
	/*
	if (T<700.0) m_rate = r_f2_peq * site_count / 6.0; //1400K is a typical cutoff.
	else m_rate = r_f2_ss * site_count / 6.0; // Rate Equation*/
	
	return m_rate ; // Rate Equation
    //return setRate0p0267(gp, pah_st);
}
// 
// ************************************************************
// ID3- BY6 closure reaction (AR14 on Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void L6_BY6::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(9.24e7, 1.5, 9.646, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(9.6e4, 1.96, 9.021, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.35468, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(1.11e11, .658, 23.99, sp::None));   //5 - r4f
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(9.24e7, 1.5, 9.646, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(9.6e4, 1.96, 9.021, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(8.02e19, -2.011, 1.968, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(2.22e11, .658, 23.99, sp::None));   //5 - r4f
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//From Raj2009 OLD METHOD
    /*addReaction(rxnV3, Reaction(9.24e7, 1.5, 9.646, sp::H));      				//0 - r1f
    addReaction(rxnV3, Reaction(9.6e4, 1.96, 9.021, sp::H2));   				//1 - r1b
    addReaction(rxnV3, Reaction(2.1e13, 0, 4.56937799, sp::OH));  				//2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); 				//3 - r2b
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));    	//4 - r3f           // Harding et al. 2005 Recommended by Frenklach2018
    //addReaction(rxnV3, Reaction(8.02e19, -2.011, 1.968, sp::H));  			//4 - r3f
    addReaction(rxnV3, Reaction(1.11e11, .658, 23.99, sp::None));   			//5 - r4f
	addReaction(rxnV3, Reaction(3.49e12, -0.39, 2.44, sp::None));   			//6 - r4f			//Second route Raj2009
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

    m_sType = BY6; // sitetype
    m_name = "BY6 closure"; // name of process
    m_ID = 3;
}
// Jump rate calculation
double L6_BY6::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = 2*m_r[5]*r_f* site_count; // Rate Equation
}
double L6_BY6::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double L6_BY6::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//OLD METHOD
	/*// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = 2*(m_r[5]+m_r[6])*r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}
// 
// ************************************************************
// ID4- phenyl addition (AR15 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void PH_benz::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(2.5e14, 0, 16.00, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(2.08e30, -4.98, 5.290, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(2e76, -18.4456, 46.93, sp::C6H6));   //5 - r4f
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(2.2e36, -8.21, 9.92, sp::C6H6));   //5 - r4f
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4              - Forward
	addReaction(rxnV3, Reaction(2.220E+83,-2.079E+01, 4.689E+01, sp::C6H6));   		 // A2* + A1 -> A2--A1 + H				- 5 */
    /*addReaction(rxnV3, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));   //1 - r1b
    addReaction(rxnV3, Reaction(1.0e10, .734, 1.43, sp::OH));  //2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV3, Reaction(2.0e13, 0, 0, sp::H));  //4 - r3f
    addReaction(rxnV3, Reaction(1.9e76, -18.4043, 47.87, sp::C6H6));   //5 - r4f*/
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f

    m_sType = benz; // sitetype
    m_name = "Phenyl addition"; // name of process
    m_ID = 4;
}
// Jump rate calculation
double PH_benz::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    if((double)pah_st.getSiteCount(R5)!=0) site_count++;
    return m_rate = 2*m_r[5]*r_f* site_count; // Rate Equation
}
double PH_benz::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double PH_benz::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
// 
// ************************************************************
// ID5- R6 desorption at FE (AR8 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void D6R_FE3::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(3.23e7, 2.095, 15.84, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(1.5e11, 0.2, 42.57, sp::None));   //5 - r4f
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(2e11, 0.2, 42.57, sp::None));   //5 - r4f
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//Commented out by GLC
	//addReaction(rxnV3, Reaction(4.20e13, 0, 13.00, sp::H));     // 0 - r1f
	//addReaction(rxnV3, Reaction(3.90e12, 0, 11.00, sp::H2));    // 1 - r1b
	//addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	//addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	//addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	//addReaction(rxnV3, Reaction(1.30e11, 1.080, 70.40, sp::None));  //5 - r4f


	//Added by GLC
	//Reactions from Narayanaswamy2010 and Mebel2016 for 1 atm.
	//rxns of A1-
	//rxns of A2-		 
	addReaction(rxnV3, Reaction(4.170e+13, 1.500e-01, 0.000e+00, sp::H));            // A2- + H <=> A2               		 - 0            - Forward
	addReaction(rxnV3, Reaction(4.126e+17, -3.390e-01, 1.195e+02, sp::None));        // A2- + H <=> A2               		 - 1            - Backward
	addReaction(rxnV3, Reaction(4.570e+08, 1.880e+00, 1.484e+01, sp::H));            // A2 + H <=> A2- + H2          		 - 2            - Forward
	addReaction(rxnV3, Reaction(2.516e+04, 2.612e+00, -9.970e-01, sp::H2));          // A2 + H <=> A2- + H2          		 - 3            - Backward
	addReaction(rxnV3, Reaction(7.665e+04, 2.105e+00, 9.394e+00, sp::OH));           // A2 + OH <=> A2- + H2O                - 4            - Forward
	addReaction(rxnV3, Reaction(1.093e+02, 2.625e+00, 8.825e+00, sp::H2O));          // A2 + OH <=> A2- + H2O                - 5            - Backward
	addReaction(rxnV3, Reaction(1.130e+16, 7.538e-01, 6.694e+01, sp::H));            // A1- + C4H4 <=> A2 + H                - 6            - Backward
	addReaction(rxnV3, Reaction(1.387e+21, -7.978e-01, 7.245e+01, sp::H));           // A1C2H3* + C2H2 <=> A2 + H            - 7            - Backward
	addReaction(rxnV3, Reaction(1.619e+29, -3.226e+00, 7.470e+01, sp::H));           // A1C2H + C2H3 <=> A2 + H              - 8            - Backward
	addReaction(rxnV3, Reaction(1.338e+94, -2.184e+01, 1.435e+02, sp::None));        // A1C2H* + C2H2 <=> A2-                - 9            - Backward
	addReaction(rxnV3, Reaction(3.917e+115, -2.755e+01, 1.688e+02, sp::None));       // A1C2H* + C2H2 <=> A2*                - 10            - Backward
	addReaction(rxnV3, Reaction(1.583e+40, -6.094e+00, 8.758e+01, sp::H));           // A1C2H* + C2H4 <=> A2 + H             - 11            - Backward
	addReaction(rxnV3, Reaction(3.736e+32, -4.109e+00, 7.723e+01, sp::H));           // A1C2H2 + C2H2 <=> A2 + H             - 12            - Backward
	addReaction(rxnV3, Reaction(5.450e+90, -2.104e+01, 1.358e+02, sp::None));        // A1C2HX2 + H <=> A2-          		 - 13            - Backward
	addReaction(rxnV3, Reaction(3.607e+113, -2.705e+01, 1.629e+02, sp::None));       // A1C2HX2 + H <=> A2*          		 - 14            - Backward
	
	//Reactions for the ABF mechanism
	/*//------------Reactions for A1-----------------------------
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));            // A1 + H <=> A1- + H2					- 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));           // A1 + H <=> A1- + H2					- 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A1 + OH <=> A1- + H2O               - 2              - Forward
	addReaction(rxnV3, Reaction(1.136E+05, 2.055E+00, 9.453E+00, sp::H2O));          // A1 + OH <=> A1- + H2O               - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A1- + H -> A1              			- 4              - Forward
	addReaction(rxnV3, Reaction(1.100e+07, 1.710e+00, 3.894e+00, sp::C2H2));           // A1- + C2H2 <=> n-A1C2H2              - 10            - Forward
	addReaction(rxnV3, Reaction(1.300e+04, 0.000e+00, 4.1778e+01, sp::None));        // A1- + C2H2 <=> n-A1C2H2             - 5              - Backward
	addReaction(rxnV3, Reaction(3.300E+33, -5.700E+00, 2.550E+01, sp::C2H2));        // A1- + C2H2 <=> A1C2H + H            - 6              - Forward
	addReaction(rxnV3, Reaction(3.933E+41, -7.568E+00, 3.344E+01, sp::H));           // A1- + C2H2 <=> A1C2H + H            - 7              - Backward
	addReaction(rxnV3, Reaction(8.500E-02, 4.710E+00, 1.842E+01, sp::C2H3));         // A1- + C2H3 <=> H + i-A1C2H2         - 8              - Forward
	addReaction(rxnV3, Reaction(1.148E+08, 2.721E+00, 3.097E+01, sp::H));            // A1- + C2H3 <=> H + i-A1C2H2         - 9              - Backward
	addReaction(rxnV3, Reaction(2.510E+12, 0.000E+00, 6.190E+00, sp::C2H4));         // A1- + C2H4 <=> A1C2H3 + H           - 10             - Forward
	addReaction(rxnV3, Reaction(1.035E+19, -1.601E+00, 1.040E+01, sp::H));           // A1- + C2H4 <=> A1C2H3 + H           - 11             - Backward
	addReaction(rxnV3, Reaction(3.300E+33, -5.700E+00, 2.550E+01, sp::C4H4));        // A1- + C4H4 <=> A2 + H               - 12             - Forward
	addReaction(rxnV3, Reaction(4.337E+44, -7.353E+00, 8.577E+01, sp::H));           // A1- + C4H4 <=> A2 + H               - 13             - Backward
	//------------Reactions for n-A1C2H2----------------------------
	addReaction(rxnV3, Reaction(1.500e+10, 0.850e+00, 1.194e+01, sp::H));              // A1C2H + H <=> n-A1C2H2               - 18            - Forward
	addReaction(rxnV3, Reaction(6.800e+11, 0.000e+00, 2.2025e+01, sp::None));        // A1C2H + H <=> n-A1C2H2              - 15             - Backward
	addReaction(rxnV3, Reaction(6.650E+06, 2.530E+00, 1.224E+01, sp::H));            // A1C2H3 + H <=> H2 + n-A1C2H2        - 16             - Forward
	addReaction(rxnV3, Reaction(4.616E+02, 3.309E+00, 3.511E+00, sp::H2));           // A1C2H3 + H <=> H2 + n-A1C2H2        - 17             - Backward
	addReaction(rxnV3, Reaction(1.600E+16, -1.330E+00, 5.400E+00, sp::C2H2));        // C2H2 + n-A1C2H2 <=> A2 + H          - 18             - Forward
	addReaction(rxnV3, Reaction(2.594E+27, -2.986E+00, 6.880E+01, sp::H));           // C2H2 + n-A1C2H2 <=> A2 + H          - 19             - Backward
	//------------Reactions for A1C2H----------------------------
	addReaction(rxnV3, Reaction(5.000E+13, 0.000E+00, 0.000E+00, sp::C2H));          // A1 + C2H <=> A1C2H + H              - 20             - Forward
	addReaction(rxnV3, Reaction(5.220E+20, -1.629E+00, 2.974E+01, sp::H));           // A1 + C2H <=> A1C2H + H              - 21             - Backward
	addReaction(rxnV3, Reaction(3.000E+43, -9.220E+00, 1.527E+01, sp::H));           // A1C2H + H <=> i-A1C2H2              - 22             - Forward
	addReaction(rxnV3, Reaction(9.521E+44, -9.486E+00, 5.563E+01, sp::None));        // A1C2H + H <=> i-A1C2H2              - 23             - Backward
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));            // A1C2H + H <=> A1C2H* + H2           - 24             - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));           // A1C2H + H <=> A1C2H* + H2           - 25             - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A1C2H + OH <=> A1C2H* + H2O         - 26             - Forward
	addReaction(rxnV3, Reaction(1.263E+05, 2.245E+00, 8.610E+00, sp::H2O));          // A1C2H + OH <=> A1C2H* + H2O         - 27             - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A1- + H -> A1              			- 4              - Forward
	addReaction(rxnV3, Reaction(5.000E+13, 0.000E+00, 0.000E+00, sp::C2H));          // A1C2H + C2H <=> A1C2H)2 + H         - 28             - Forward
	addReaction(rxnV3, Reaction(4.296E+20, -1.286E+00, 2.828E+01, sp::H));           // A1C2H + C2H <=> A1C2H)2 + H         - 29             - Backward
	//------------Reactions for i-A1C2H2----------------------------
	addReaction(rxnV3, Reaction(3.330E+05, 2.530E+00, 9.240E+00, sp::H));            // A1C2H3 + H <=> H2 + i-A1C2H2        - 30             - Forward
	addReaction(rxnV3, Reaction(6.407E+03, 2.874E+00, 9.069E+00, sp::H2));           // A1C2H3 + H <=> H2 + i-A1C2H2        - 31             - Backward
	addReaction(rxnV3, Reaction(1.550E+06, 2.000E+00, 4.300E-01, sp::OH));           // A1C2H3 + OH <=> H2O + i-A1C2H2      - 32             - Forward
	addReaction(rxnV3, Reaction(1.095E+06, 2.090E+00, 1.607E+01, sp::H2O));          // A1C2H3 + OH <=> H2O + i-A1C2H2      - 33             - Backward
	//------------Reactions for A1C2H3----------------------------
	//------------Reactions for A1C2H*----------------------------
	//       A1C2H* + H (+M) <=> A1C2H (+M) Reaction ommitted
	addReaction(rxnV3, Reaction(2.200E+62, -1.456E+01, 3.310E+01, sp::C2H2));        // A1C2H* + C2H2 <=> A2-1              - 34             - Forward
	addReaction(rxnV3, Reaction(5.164E+72, -1.618E+01, 1.289E+02, sp::None));        // A1C2H* + C2H2 <=> A2-1              - 35             - Backward
	addReaction(rxnV3, Reaction(1.800E+19, -1.670E+00, 1.880E+01, sp::C2H2));        // A1C2H* + C2H2 <=> A1C2H)2 + H       - 36             - Forward
	addReaction(rxnV3, Reaction(1.588E+27, -3.384E+00, 2.613E+01, sp::H));           // A1C2H* + C2H2 <=> A1C2H)2 + H       - 37             - Backward
	addReaction(rxnV3, Reaction(5.700E+64, -1.441E+01, 5.700E+01, sp::C2H2));        // A1C2H* + C2H2 <=> H + naphthyne     - 38             - Forward
	addReaction(rxnV3, Reaction(6.479E+71, -1.544E+01, 7.452E+01, sp::H));           // A1C2H* + C2H2 <=> H + naphthyne     - 39             - Backward
	//------------Reactions for A1C2H)2----------------------------
	addReaction(rxnV3, Reaction(1.400E+64, -1.457E+01, 2.990E+01, sp::H));           // A1C2H)2 + H <=> A2-1                - 40             - Forward
	addReaction(rxnV3, Reaction(3.725E+66, -1.447E+01, 1.183E+02, sp::None));        // A1C2H)2 + H <=> A2-1                - 41             - Backward
	addReaction(rxnV3, Reaction(1.900E+73, -1.630E+01, 6.090E+01, sp::H));           // A1C2H)2 + H <=> H + naphthyne       - 42             - Forward
	addReaction(rxnV3, Reaction(2.448E+72, -1.561E+01, 7.109E+01, sp::H));           // A1C2H)2 + H <=> H + naphthyne       - 43             - Backward
	//------------Reactions for naphthyne----------------------------
	//------------Reactions for A2----------------------------
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));            // A2 + H <=> A2-1 + H2                - 44             - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));           // A2 + H <=> A2-1 + H2                - 45             - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A2 + OH <=> A2-1 + H2O              - 46             - Forward
	addReaction(rxnV3, Reaction(1.597E+05, 2.116E+00, 9.183E+00, sp::H2O));          // A2 + OH <=> A2-1 + H2O              - 47             - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A1- + H -> A1              			- 4              - Forward
	//I STARTED THE SECOND PART HERE!!!
	//------------3 body reactions------------------------------
	addReaction(rxnV3, Reaction(1.000E+14, 0.000E+00, 0.000E+00, sp::None));            // A1- + H (+M) <=> A1 (+M)            - 48              - Forward - High rate
	addReaction(rxnV3, Reaction(6.600E+72, -1.630E+01, 7.000E+00, sp::None));           // A1- + H (+M) <=> A1 (+M)            - 49              - Forward - High rate
	addReaction(rxnV3, Reaction(1.392E+34, -5.237E+00, 1.200E+02, sp::None));           // A1- + H (+M) <=> A1 (+M)            - 50              - Backward
	//Parameters = [1.00000000e+00   1.00000000e-01   5.84900000e+02   6.11300000e+03]
	addReaction(rxnV3, Reaction(1.000E+14, 0.000E+00, 0.000E+00, sp::None));            // A1C2H* + H (+M) <=> A1C2H (+M)      - 51              - Forward - High rate
	addReaction(rxnV3, Reaction(6.600E+72, -1.630E+01, 7.000E+00, sp::None));           // A1C2H* + H (+M) <=> A1C2H (+M)      - 52              - Forward - High rate
	addReaction(rxnV3, Reaction(1.252E+34, -5.427E+00, 1.209E+02, sp::None));           // A1C2H* + H (+M) <=> A1C2H (+M)      - 53              - Backward
	//Parameters = [1.00000000e+00   1.00000000e-01   5.84900000e+02   6.11300000e+03]
	addReaction(rxnV3, Reaction(1.000E+14, 0.000E+00, 0.000E+00, sp::None));            // A2-1 + H (+M) <=> A2 (+M)           - 54              - Forward - High rate
	addReaction(rxnV3, Reaction(3.800E+124, -3.143E+01, 1.868E+01, sp::None));          // A2-1 + H (+M) <=> A2 (+M)           - 55              - Forward - High rate
	addReaction(rxnV3, Reaction(2.662E+37, -6.271E+00, 1.225E+02, sp::None));           // A2-1 + H (+M) <=> A2 (+M)           - 56              - Backward
	//Parameters = [2.00000000e-01   1.22800000e+02   4.78400000e+02   5.41190000e+03]*/


	m_sType = FE3; // sitetype
	m_name = "R6 (FE3) desorption"; // name of process
	m_ID = 5;
}
// Jump rate calculation
double D6R_FE3::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
		r_f = r_f / (r_f + 1.0);
	}
	else r_f = 0;
	return m_rate = m_r[5] * r_f* site_count; // Rate Equation
}
double D6R_FE3::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate0p0267(gp, pah_st);
}
double D6R_FE3::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	//int ring_count = ((int)pah_st.getRingsCount().first); // Ring count
	//auto ring_count = pah_st.getRingsCount();
	//int ring_count_6 = std::get<0>(ring_count);
	//if (ring_count_6 <= 1) return m_rate = 0;
	if (site_count == 0) return m_rate = 0;
	

	// Calculating 3 body reactions for ABF
	/*//Parameters alpha = 0.2, T3 = 122.8, T1 = 478.4, T2 = 5411.9
	//Parameters = [1.00000000e+00   1.00000000e-01   5.84900000e+02   6.11300000e+03]
	double T = gp[gp.T];
	double P = gp[gp.P]; //Pressure in Pa
	double M_eff1 = P / 8.314E+06 / T;
	double M = M_eff1 + (1 * gp[gp.H2]) + (5 * gp[gp.H2O]) + (1 * gp[gp.CH4]) + (0.5*gp[gp.CO]) + (1 * gp[gp.CO2]) + (2 * gp[gp.C2H6]);
	double k_o1 = m_r[49];
	double k_inf1 = m_r[48];
	double F_cent1 = 0.00 * exp(-T / 0.1) + 1.00*exp(-T / 584.9) + exp(-6113 / T);
	double logF_cent1 = log10(F_cent1);
	double r_c1 = k_inf1 / (1 + (k_inf1 / (k_o1*M)));
	double N01 = 0.75 - 1.27*logF_cent1;
	double C01 = -0.4 - 0.67*logF_cent1;
	double logP_r1 = log10(k_o1*M / k_inf1);
	double F1 = exp(logF_cent1 / (1 + pow(((logP_r1 + C01) / (N01 - 0.14*(logP_r1 + C01))), 2)));
	r_c1 *= F1;
	double R3b1 = r_c1* gp[gp.H];
	//Parameters = [1.00000000e+00   1.00000000e-01   5.84900000e+02   6.11300000e+03]
	double k_o2 = m_r[52];
	double k_inf2 = m_r[51];
	double F_cent2 = 0.00 * exp(-T / 0.1) + 1.00*exp(-T / 584.9) + exp(-6113 / T);
	double logF_cent2 = log10(F_cent2);
	double r_c2 = k_inf2 / (1 + (k_inf2 / (k_o2*M)));
	double N02 = 0.75 - 1.27*logF_cent2;
	double C02 = -0.4 - 0.67*logF_cent2;
	double logP_r2 = log10(k_o2*M / k_inf2);
	double F2 = exp(logF_cent2 / (1 + pow(((logP_r2 + C02) / (N02 - 0.14*(logP_r2 + C02))), 2)));
	r_c2 *= F2;
	double R3b2 = r_c2* gp[gp.H];
	//Parameters = [2.00000000e-01   1.22800000e+02   4.78400000e+02   5.41190000e+03]
	double k_o3 = m_r[55];
	double k_inf3 = m_r[54];
	double F_cent3 = 0.20 * exp(-T / 122.8) + 0.80*exp(-T / 478.4) + exp(-5411.9 / T);
	double logF_cent3 = log10(F_cent3);
	double r_c3 = k_inf3 / (1 + (k_inf3 / (k_o3*M)));
	double N03 = 0.75 - 1.27*logF_cent3;
	double C03 = -0.4 - 0.67*logF_cent3;
	double logP_r3 = log10(k_o3*M / k_inf3);
	double F3 = exp(logF_cent3 / (1 + pow(((logP_r3 + C03) / (N03 - 0.14*(logP_r3 + C03))), 2)));
	r_c3 *= F3;
	double R3b3 = r_c3* gp[gp.H];

	//double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	//if (r_denom>0) {
	//r_f = (m_r[0] + m_r[2]) / r_denom;
	//}
	//else r_f = 0;

	//Used for ABF
	matrix<double> arr1(10, 10);
	boost::numeric::ublas::vector<double> arr2(10);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;

	//Speclist = ['A2*', 'n-A2C2H2', 'A2C2H', 'i-A2C2H2', 'A2C2H3','A2C2H*','A2C2H)2','naphthyneA2','A3','A3*']
	arr1(0, 0) = m_r[1] + m_r[3] + m_r[4] + m_r[6] + m_r[8] + m_r[10] + m_r[12] + R3b1; //A1- consumption
	arr1(0, 1) = -(m_r[5]);
	arr1(0, 2) = -(m_r[7]);
	arr1(0, 3) = -(m_r[9]);
	arr1(0, 4) = -(m_r[11]);
	arr1(0, 8) = -(m_r[13]);
	arr1(1, 0) = -(m_r[4]);
	arr1(1, 1) = m_r[5] + m_r[15] + m_r[17] + m_r[5] + m_r[18];
	arr1(1, 2) = -(m_r[14]);
	arr1(1, 4) = -(m_r[16]);
	arr1(1, 8) = -(m_r[19]);
	arr1(2, 0) = -(m_r[6]);
	arr1(2, 1) = -(m_r[15]);
	arr1(2, 2) = m_r[7] + m_r[14] + m_r[21] + m_r[22] + m_r[24] + m_r[26] + m_r[28] + m_r[53];
	arr1(2, 3) = -(m_r[23]);
	arr1(2, 5) = -(m_r[25] + m_r[27] + R3b2);
	arr1(2, 6) = -(m_r[29]);
	arr1(3, 0) = -(m_r[8]);
	arr1(3, 2) = -(m_r[22]);
	arr1(3, 3) = m_r[9] + m_r[23] + m_r[31] + m_r[33];
	arr1(3, 4) = -(m_r[30] + m_r[32]);
	arr1(4, 0) = -(m_r[10]);
	arr1(4, 1) = -(m_r[17]);
	arr1(4, 3) = -(m_r[31] + m_r[33]);
	arr1(4, 4) = m_r[11] + m_r[16] + m_r[30] + m_r[32];
	arr1(5, 2) = -(m_r[24] + m_r[26] + m_r[53]);
	arr1(5, 5) = m_r[25] + m_r[27] + m_r[34] + m_r[36] + m_r[38] + R3b2;
	arr1(5, 6) = -(m_r[37]);
	arr1(5, 7) = -(m_r[39]);
	arr1(5, 9) = -(m_r[35]);
	arr1(6, 2) = -(m_r[28]);
	arr1(6, 5) = -(m_r[36]);
	arr1(6, 6) = m_r[29] + m_r[37] + m_r[40] + m_r[42];
	arr1(6, 7) = -(m_r[39]);
	arr1(6, 9) = -(m_r[41]);
	arr1(7, 5) = -(m_r[38]);
	arr1(7, 6) = -(m_r[42]);
	arr1(7, 7) = m_r[39] + m_r[43];
	arr1(8, 0) = -(m_r[12]);
	arr1(8, 1) = -(m_r[18]);
	arr1(8, 8) = m_r[13] + m_r[19] + m_r[44] + m_r[46] + m_r[56];
	arr1(8, 9) = -(m_r[45] + m_r[47] + R3b3);
	arr1(9, 5) = -(m_r[34]);
	arr1(9, 6) = -(m_r[40]);
	arr1(9, 8) = -(m_r[44] + m_r[46] + m_r[56]);
	arr1(9, 9) = m_r[35] + m_r[41] + m_r[45] + m_r[47] + R3b3;

	arr2(0) = m_r[0] + m_r[2] + m_r[50];
	arr2(2) = m_r[20];

	//std::cout << arr1 << std::endl;
	//std::cout << arr2 << std::endl;



	permutation_matrix<size_t> pm(arr1.size1());

	//std::cout << arr1 << std::endl;
	//std::cout << arr2 << std::endl;

	lu_factorize(arr1, pm);

	//std::cout << arr1 << std::endl;
	//std::cout << arr2 << std::endl;

	lu_substitute(arr1, pm, arr2);

	//std::cout << arr1 << std::endl;
	//std::cout << arr2 << std::endl;;


	//double PA2toA3 = arr2(0)*(m_r[45] + m_r[47] + R3b3);
	//double PA2toA1 = arr2(9)*(m_r[45] + m_r[47] + R3b3);

	
	r_f = (m_r[44] + m_r[46] + m_r[56]) / (m_r[45] + m_r[47] + R3b3);
	double r_f_check = arr2(0);
	//double tot_A2_cons = m_r[44] + m_r[46];
	//double PA2toA1 = tot_A2_cons - PA2toA3;
	//double test = m_r[44] + m_r[46]; //This rate grows well back to A1.
	//double prod_A2 = (m_r[45] + m_r[47] + R3b3) * r_f;
	//double test2 = 1.0;
	m_rate = (m_r[41] + m_r[35]) * r_f* site_count / 2.0; // Rate Equation*/
	//m_rate = PA2toA1;
	//New rates
	double r_denom = (m_r[0] + m_r[3] + m_r[5]);
	//double r_denom = (m_r[0] + m_r[3] + m_r[5] + m_r[9] + m_r[13]);
	if (r_denom <=0) return m_rate = 0.0;
	double r_num = (m_r[1] + m_r[2] + m_r[4]);
	m_rate = ((m_r[6] + m_r[7] + m_r[8] + m_r[11] + m_r[12]) + r_num/r_denom*(m_r[9] + m_r[10] + m_r[13] + m_r[14])) * site_count / 2.0; 
	
	return m_rate;
	


    //return setRate0p0267(gp, pah_st);
}
// ************************************************************
// ID6- R6 oxidation at FE by O2 (AR10 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void O6R_FE3_O2::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(3.23e7, 2.095, 15.84, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));   //5 - r5f(1)
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f(2)
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
	addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));   //5 - r5f(1)
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //6 - r5f(2)
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    addReaction(rxnV3, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));   //1 - r1b
    addReaction(rxnV3, Reaction(1.0e10, 0.734, 1.43, sp::OH));  //2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV3, Reaction(2.0e13, 0, 0, sp::H));  //4 - r3f
    addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));   //5 - r4f
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
    m_sType = FE3; // sitetype
    m_name = "R6 (FE3) oxidation by O2"; // name of process
    m_ID = 6;
}
// Jump rate calculation
double O6R_FE3_O2::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation
}
double O6R_FE3_O2::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation
}
double O6R_FE3_O2::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p12(gp, pah_st);
}
// ************************************************************
// ID7- R6 oxidation at FE by OH (AR11 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void O6R_FE3_OH::initialise() {
    // Adding elementary reactions
    // 0.0267atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(3.23e7, 2.095, 15.84, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(1.3e11, 1.08, 70.42, sp::OH));   //5 - r4f
    addReaction(rxnV, Reaction(1.3e13, 0, 10.62, sp::OH));          //6 - r5f
    // 0.12atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    //addReaction(rxnV2, Reaction(4.2e13, 0, 13, sp::H));      //0 - r1f
    //addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    //addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    //addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    //addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(1.3e13, 0, 10.62, sp::OH));   //5 - r4f
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //6 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    //addReaction(rxnV3, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    //addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));   //1 - r1b
    //addReaction(rxnV3, Reaction(1.0e10, 0.734, 1.43, sp::OH));  //2 - r2f
    //addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    //addReaction(rxnV3, Reaction(2.0e13, 0, 0, sp::H));  //4 - r3f
    addReaction(rxnV3, Reaction(1.3e13, 0, 10.62, sp::OH));   //5 - r4f
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
    m_sType = FE3; // sitetype
    m_name = "R6 (FE3) oxidation by OH"; // name of process
    m_ID = 7;
}
// Jump rate calculation
double O6R_FE3_OH::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[6]*r_f* site_count; // Rate Equation
}
double O6R_FE3_OH::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    //double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    //double r_f; // radical fraction 
    //if(r_denom>0) {
    //    r_f = (m_r[0]+m_r[2])/r_denom; 
    //    r_f = r_f/(r_f+1.0);
    //}
    //else r_f=0;
	return m_rate = m_r[0] * site_count; // Rate Equation
}
double O6R_FE3_OH::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p12(gp, pah_st);
}
// ************************************************************
// ID8- R6 oxidation at AC by O2 (AR12 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void O6R_FE_HACA::initialise() {
    // Adding elementary reactions
    // 0.0267atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(3.23e7, 2.095, 15.84, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));   //5 - r5f(1)
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f(2)
    // 0.12atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
	addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));   //5 - r5f(1)
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //6 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4
	addReaction(rxnV3, Reaction(3.91e13, 0.00, 1496.0*8.314/4.184/1000.0,sp::O2));  						// S48 // S49	- 5
    addReaction(rxnV3, Reaction(1.76e23, -3.681, 24155.0*8.314/4.184/1000.0,sp::None));  // S53 S54                         - 6	//Using reaction from He Lin et al @1atm
    addReaction(rxnV3, Reaction(1.47e14, 0.00, 632.0*8.314/4.184/1000.0,sp::OH)); // S57 //S58                              - 7
    addReaction(rxnV3, Reaction(2.13e15, 0.00, 42305.0*8.314/4.184/1000.0,sp::None));  // S62 // S63                        - 8
    addReaction(rxnV3, Reaction(4.34e14, 0.00, 984.0*8.314/4.184/1000.0,sp::H));  // S67 // S68                             - 9
    addReaction(rxnV3, Reaction(2.14e16, 0.00, 52988.0*8.314/4.184/1000.0,sp::None));  // S72 // S73                        - 10
    addReaction(rxnV3, Reaction(5476782.795, 2.147, 5603.0*8.314/4.184/1000.0,sp::H));  // S77 // S78                        - 11
    addReaction(rxnV3, Reaction(1.70e14, 0.00, 9635.0*8.314/4.184/1000.0,sp::H2));  // S82 // S83                           - 12
    addReaction(rxnV3, Reaction(2.00e14, 0.00, 2670.0*8.314/4.184/1000.0,sp::H));  // S87 // S88                            - 13
    addReaction(rxnV3, Reaction(4.00e12, 0.00, 2328.0*8.314/4.184/1000.0,sp::O));  // S94                                   - 14
    addReaction(rxnV3, Reaction(0.0*1.00e14, 0.00, 0.0*8.314/4.184/1000.0,sp::OH));  // S97 // S98                              - 15
    addReaction(rxnV3, Reaction(0.0*1.00e10, 0.00, 0.0*8.314/4.184/1000.0,sp::H2O));  // S101 // S102                       - 16 //It is 1.00e12 in Frenkalch's paper but that's too high. Lowered here.
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f                                       - 
	addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 17              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 18              - Frenklach et al. 2018
    m_sType = FE_HACA; // sitetype
    m_name = "R6 (FE_HACA) oxidation"; // name of process
    m_ID = 8;
}
// Jump rate calculation
double O6R_FE_HACA::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
	/*// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation*/
}
double O6R_FE_HACA::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
	/*// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation*/
}
double O6R_FE_HACA::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, C-O, C-OH
	arr1(0,0) = m_r[1] + m_r[3] + m_r[4] + m_r[5] + m_r[7]+m_r[15]+m_r[16] + m_r[17] + m_r[18];
    arr1(0,2) = -m_r[10] - m_r[13];
    arr1(1,0) = -m_r[5];
    arr1(1,1) = m_r[6] + m_r[9] + m_r[12];
    arr1(1,2) = - m_r[8] - m_r[11];
	arr1(2,0) = -m_r[7];
    arr1(2,1) = - m_r[9] - m_r[12];
    arr1(2,2) = m_r[8] + m_r[10] + m_r[11] + m_r[13];

    arr2(0) = +m_r[0] + m_r[2];
    arr2(1) = +m_r[14];
    arr2(2) = 0.0;

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	double Crad_fraction = arr2(0);
	double CO_fraction = arr2(1);
	double COH_fraction = arr2(2);
	if (CO_fraction >= 1.0 || COH_fraction >= 1.0) return m_rate = 0.0;
	else return m_rate = (m_r[6] * CO_fraction + (m_r[15] + m_r[16]) * Crad_fraction)*site_count*2.0;
	/*if (arr2(1) < arr2(0) && arr2(2) < arr2(0)) {
		return m_rate = (m_r[6] * arr2(1) + m_r[16] * arr2(0) + m_r[6] * arr2(4) + m_r[16] * arr2(3))*site_count;
	}
	else return m_rate = 0.0; // The SS approximation breaks so the rate approximation is wrong.*/
}

// ************************************************************
// ID9- R6 oxidation at AC by OH (AR13 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void O6R_FE_HACA_OH::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(3.23e7, 2.095, 15.84, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(1.3e11, 1.08, 70.42, sp::OH));   //5 - r4f
    addReaction(rxnV, Reaction(1.3e13, 0, 10.62, sp::OH));          //6 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    //addReaction(rxnV2, Reaction(4.2e13, 0, 13, sp::H));      //0 - r1f
    //addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    //addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    //addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    //addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(1.3e13, 0, 10.62, sp::OH));   //5 - r4f
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::OH));          //6 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    //addReaction(rxnV3, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    //addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));   //1 - r1b
    //addReaction(rxnV3, Reaction(1.0e10, 0.734, 1.43, sp::OH));  //2 - r2f
    //addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    //addReaction(rxnV3, Reaction(2.0e13, 0, 0, sp::H));  //4 - r3f
    addReaction(rxnV3, Reaction(1.3e13, 0, 10.62, sp::OH));   //5 - r4f
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
    m_sType = FE_HACA; // sitetype
    m_name = "R6 (FE_HACA) oxidation by OH"; // name of process
    m_ID = 9;
}
// Jump rate calculation
double O6R_FE_HACA_OH::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;

    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[6]*r_f* site_count; // Rate Equation
}
double O6R_FE_HACA_OH::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;

    //double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    //double r_f; // radical fraction 
    //if(r_denom>0) {
    //    r_f = (m_r[0]+m_r[2])/r_denom; 
    //    r_f = r_f/(r_f+1.0);
    //}
    //else r_f=0;
	return m_rate = m_r[0] * site_count; // Rate Equation
}
double O6R_FE_HACA_OH::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p12(gp, pah_st);
}

// ************************************************************
// ID10- R5 growth on ZZ (AR3 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void G5R_ZZ::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(2.50e14, 0, 16.00, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(3.57e24, -3.176, 14.86, sp::C2H2));   //5 - r4f
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f(2)
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(1.8e33, -5.91, 19.7, sp::C2H2));   //5 - r4f
    addReaction(rxnV2, Reaction(3.4e43, -9.56, 18.2, sp::C2H2));   //6 - r4f (T<800)
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //7 - r5f(2)
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//Frenklach2018 and Blanquart
    addReaction(rxnV3, Reaction(4.170e+13, 1.500e-01, 0.000e+00, sp::H));            // A2- + H <=> A2               		 - 0            - Forward
	addReaction(rxnV3, Reaction(4.126e+17, -3.390e-01, 1.195e+02, sp::None));        // A2- + H <=> A2               		 - 1            - Backward
	addReaction(rxnV3, Reaction(4.570e+08, 1.880e+00, 1.484e+01, sp::H));            // A2 + H <=> A2- + H2          		 - 2            - Forward
	addReaction(rxnV3, Reaction(2.516e+04, 2.612e+00, -9.970e-01, sp::H2));          // A2 + H <=> A2- + H2          		 - 3            - Backward
	addReaction(rxnV3, Reaction(7.665e+04, 2.105e+00, 9.394e+00, sp::OH));           // A2 + OH <=> A2- + H2O                - 4            - Forward
	addReaction(rxnV3, Reaction(1.093e+02, 2.625e+00, 8.825e+00, sp::H2O));          // A2 + OH <=> A2- + H2O                - 5            - Backward
	addReaction(rxnV3, Reaction(1.250e+27, -3.950e+00, 1.6779e+01, sp::C2H2));       // A2- + C2H2 <=> A2R5 + H              - 6            - Forward
	addReaction(rxnV3, Reaction(3.090e+20, -2.780e+00, 8.8890e+00, sp::C2H2));       // A2- + C2H2 <=> A2R5 + H              - 7            - Forward
	addReaction(rxnV3, Reaction(3.090e+25, -3.110e+00, 3.1586e+01, sp::C2H2));       // A2- + C2H2 <=> A2C2H + H             - 8            - Forward
	addReaction(rxnV3, Reaction(2.850e+07,  1.520e+00, 1.3190e+01, sp::C2H2));       // A2- + C2H2 <=> A2C2H + H             - 9            - Forward
	/*addReaction(rxnV3, Reaction(1.250e+27, -3.950e+00, 1.6495e+01, sp::C2H2));     // A5- + C2H2 <=> A5R5 + H              - 10            - Forward
	addReaction(rxnV3, Reaction(3.090e+20, -2.780e+00, 8.8890e+00, sp::C2H2));       // A5- + C2H2 <=> A5R5 + H              - 11            - Forward
	addReaction(rxnV3, Reaction(3.090e+25, -3.110e+00, 3.1586e+01, sp::C2H2));       // A5- + C2H2 <=> A5C2H + H             - 12            - Forward
	addReaction(rxnV3, Reaction(2.850e+07,  1.520e+00, 1.3190e+01, sp::C2H2));       // A5- + C2H2 <=> A5C2H + H             - 13            - Forward*/
	
	//ABF
	/*addReaction(rxnV3, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));   //1 - r1b
    addReaction(rxnV3, Reaction(1.0e10, 0.734, 1.43, sp::OH));  //2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV3, Reaction(2.0e13, 0, 0, sp::H));  //4 - r3f
    addReaction(rxnV3, Reaction(6.8e11, 0, 22.0, sp::C2H2));   //5 - r4f
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f*/
    m_sType = ZZ; // sitetype
    m_name = "R5 growth on ZZ"; // name of process
    m_ID = 10;
}
// Jump rate calculation
double G5R_ZZ::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation
}
double G5R_ZZ::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    double r4f = m_r[5];
    if(gp[gp.T] < 800) r4f = m_r[6];
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+r4f);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r4f*r_f* site_count; // Rate Equation
}
double G5R_ZZ::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    double r_denom = (m_r[0]+m_r[3]+m_r[5]+m_r[6]+m_r[7]+m_r[8]+m_r[9]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[2]+m_r[4])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = (m_r[6] + m_r[7]) * r_f * site_count/2.0; // Rate Equation
}


// ************************************************************
// ID11- R5 desorption (AR7 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void D5R_R5::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(5.07e7, 1.93, 12.95, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(9.45e3, 2.56, 5.007, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    //r3f to be calculated in setRate0p0267
    addReaction(rxnV, Reaction(1.60e14, 0, 42.42, sp::None));   //4 - r4f
    addReaction(rxnV, Reaction(1.4e30, -3.86, 3.320, sp::None)); // 5 - k_o: rate calculation
    addReaction(rxnV, Reaction(6.08e12, .27, 0.280, sp::None)); // 6 - k_inf: rate calculation
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //7 - r5f(2)
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    rxnV2 = rxnV;
    //addReaction(rxnV2, Reaction(5.07e7, 1.93, 12.95, sp::H));      //0 - r1f
    //addReaction(rxnV2, Reaction(9.45e3, 2.56, 5.007, sp::H2));   //1 - r1b
    //addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    //addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    ////r3f to be calculated in setRate0p0267
    //addReaction(rxnV2, Reaction(1.60e14, 0, 42.42, sp::None));   //4 - r4f
    //addReaction(rxnV2, Reaction(1.4e30, -3.86, 3.320, sp::None)); // 5 - k_o: rate calculation
    //addReaction(rxnV2, Reaction(6.08e12, .27, 0.280, sp::None)); // 6 - k_inf: rate calculation
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //7 - r5f(2)
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    //Frenklach2018 and Blanquart -- Could not use because of lack of A2R5- thermal data (Most mechanisms assume A2R5- is the radical from any hydrogen abstractions except from R5 site)
	/*addReaction(rxnV3, Reaction(6.080e+12, 2.700e-01, 0.000e+00, sp::H));          // A2R5- + H -> A2R5            - 0 Frenklach2005
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));           // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(2.516e+04, 2.612e+00, -9.970e-01, sp::H2));          // A2 + H <=> A2- + H2          - 2 Fitted
	addReaction(rxnV3, Reaction(7.665e+04, 2.105e+00, 9.394e+00, sp::OH));           // A2 + OH <=> A2- + H2O        - 3            - Forward
	addReaction(rxnV3, Reaction(1.093e+02, 2.625e+00, 8.825e+00, sp::H2O));          // A2 + OH <=> A2- + H2O                - 4            - Backward
	
	addReaction(rxnV3, Reaction(1.736e+35, -5.235e+00, 6.138e+01, sp::H));     // A2- + C2H2 <=> A2R5 + H              - 6             - Backward
	addReaction(rxnV3, Reaction(4.292e+28, -4.065e+00, 5.359e+01, sp::H));     // A2- + C2H2 <=> A2R5 + H              - 7             - Backward*/
	
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	//addReaction(rxnV3, Reaction(5.070e+07, 1.930e+00, 1.2951e+01, sp::H));        //0 - r1f //R10 in Whitesides2010
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017       //R11 in Whitesides2010
	addReaction(rxnV3, Reaction(6.080e+12, 2.700e-01, 0.000e+00, sp::H));           // A2R5- + H -> A2R5            - 2 Frenklach2005  //R11 in Whitesides2010
	addReaction(rxnV3, Reaction(1.0e10, 0.734, 1.43, sp::OH));  					//3 - r2f 
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); 					//4 - r2b
	addReaction(rxnV3, Reaction(1.6e14, 0, 42.42, sp::None));   					//5 - r4f Frenklach1998 DFT with very old methods. I would not recommend it: gl413
	
	//H addition & migration
	addReaction(rxnV3, Reaction(5.400e+11, 4.500e-01, 1.820e+00, sp::H));          //6  //R12 in Whitesides2010
	addReaction(rxnV3, Reaction(1.791e+00, 0.000e+00, -3.5187e+01, sp::None));      //7  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(2.000e+12, 0.000e+00, 0.000e+00, sp::H));          //8  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(3.100e+11, 8.700e-01, 7.4323e+01, sp::None));      //9  //R7 in Whitesides2010
	addReaction(rxnV3, Reaction(6.700e+11, 8.400e-01, 7.0790e+01, sp::None));      //10 //R8 in Whitesides2010
	
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
	//addReaction(rxnV3, Reaction(1.736e+35, -5.235e+00, 6.138e+01, sp::H));     // A2- + C2H2 <=> A2R5 + H              - 6             - Backward
	//addReaction(rxnV3, Reaction(4.292e+28, -4.065e+00, 5.359e+01, sp::H));     // A2- + C2H2 <=> A2R5 + H              - 7             - Backward
	//route by hydrogen addition
	
    m_sType = R5; // sitetype
    m_name = "R5 desorption"; // name of process
    m_ID = 11;
}
// Jump rate calculation
double D5R_R5::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;

    // Calculating r3f
    double T = gp[gp.T];
	double P = gp[gp.P]; //Pressure in Pa
	double M_eff1 = P / 8.314E+06 / T;
	double M = M_eff1 + (1 * gp[gp.H2]) + (5 * gp[gp.H2O]) + (1 * gp[gp.CH4]) + (0.5*gp[gp.CO]) + (1 * gp[gp.CO2]) + (2 * gp[gp.C2H6]);
    double k_o = m_r[5];
    double k_inf = m_r[6];
    double F_cent = 0.218 * exp(-T/207.5) + 0.782*exp(-T/2663) + exp(-6095/T);
    double logF_cent = log10(F_cent);
    double r_c = k_inf/(1+ (k_inf/(k_o*M)));
    double N0 = 0.75- 1.27*logF_cent;
    double C0 = -0.4 - 0.67*logF_cent;
    double logP_r = log10(k_o*M/k_inf);
    double F = pow(10,(logF_cent/(1+ pow(((logP_r+C0)/(N0-0.14*(logP_r+C0))),2))));
    r_c *= F;
    double r3f = r_c* gp[gp.H];

    double r_denom = (m_r[1]+m_r[3]+r3f+m_r[4]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[4]*r_f* site_count; // Rate Equation
}
double D5R_R5::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double D5R_R5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    //Frenklach and Blanquart
    double rf_rad_denom = m_r[1] + m_r[2] + m_r[4] + m_r[5];
	double rf_rad;
	if (rf_rad_denom >0){
		rf_rad = (m_r[0] + m_r[3]) / rf_rad_denom;
	}
	else rf_rad = 0;
	double rf_H2_denom = m_r[6]/m_r[7] + m_r[8] + m_r[9] + m_r[10];
	double rf_H2;
	if (rf_H2_denom >0){
		rf_H2 = (m_r[6]) / rf_H2_denom;
	}
	else rf_H2 = 0;
	return m_rate = rf_rad*m_r[5] + rf_H2*(m_r[9] + m_r[10]);
	
	//ABF
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation*/
}
// ************************************************************
// ID12- R6 conversion to R5 (AR9 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void C6R_AC_FE3::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(3.23e7, 2.095, 15.84, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(1.5e11, 0.2, 42.57, sp::None));   //5 - r4f
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f(2)
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(4.2e13, 00, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(2e11, 0.2, 42.57, sp::None));   //5 - r4f
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //6 - r5f(2)
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // An + H <=> An* + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // An + H <=> An* + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // An + OH <=> An* + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // An + OH <=> An* + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // An* + H -> An              			- 4              - Forward
	//ABF
    /*addReaction(rxnV3, Reaction(4.2e13, 00, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));   //1 - r1b
    addReaction(rxnV3, Reaction(1.0e10, 0.734, 1.43, sp::OH));  //2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV3, Reaction(2.0e13, 0, 0, sp::H));  //4 - r3f*/
    addReaction(rxnV3, Reaction(1.3e11, 1.08, 70.42, sp::None));   //5 - r4f  //Frenklach2005
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //6 - r5f(2)
    m_sType = AC_FE3; // sitetype
    m_name = "R6 conversion to R5 on AC_FE3"; // name of process
    m_ID = 12;
}
// Jump rate calculation
double C6R_AC_FE3::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation
}
double C6R_AC_FE3::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double C6R_AC_FE3::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID13- R5 conversion to R6 on RFE (AR5 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void C5R_RFE::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(2.50e14, 0, 16.00, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(7.50e26, -3.96, 17.10, sp::C2H2));   //5 - r4f1
    addReaction(rxnV, Reaction(4.40e49, -11.6, 19.30, sp::C2H2));   //6 - r4f2 (r4f = r4f1+r4f2)
    addReaction(rxnV, Reaction(8.42e8, 1.49, 0.990, sp::H));         //7 - for rate calc
    addReaction(rxnV, Reaction(6.28e37, -8.24, 44.67, sp::None));       //8 - for rate calc
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //9 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(9.1e24, -3.39, 20.4, sp::C2H2));   //5 - r4f1
    addReaction(rxnV2, Reaction(9.90e41, -9.26, 15.7, sp::C2H2));   //6 - r4f2 (r4f = r4f1+r4f2)
    addReaction(rxnV2, Reaction(8.42e8, 1.49, 0.990, sp::H));         //7 - for rate calc
    addReaction(rxnV2, Reaction(6.28e37, -8.24, 44.67, sp::None));       //8 - for rate calc
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //9 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//H abstraction
	//addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // An + H <=> An* + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // An + H <=> An* + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // An + OH <=> An* + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // An + OH <=> An* + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // An* + H -> An              			- 4              - Forward
	//addReaction(rxnV3, Reaction(1.1e7, 1.61, 3.896, sp::C2H2));   //5 - r4f1						//Needs to be updated. This comes from 1994.
    //addReaction(rxnV3, Reaction(3.3e33, -5.7, 25.5, sp::C2H2));   //6 - r4f2 (r4f = r4f1+r4f2) 	//Needs to be updated. This comes from 1994.
	addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 5              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 6              - Frenklach et al. 2018
	
	//H addition & migration
	addReaction(rxnV3, Reaction(5.400e+11, 4.500e-01, 1.820e+00, sp::H));          //7  //R12 in Whitesides2010
	addReaction(rxnV3, Reaction(1.791e+00, 0.000e+00, -3.5187e+01, sp::None));      //8  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(2.000e+12, 0.000e+00, 0.000e+00, sp::H));          //9  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(3.100e+11, 8.700e-01, 7.4323e+01, sp::None));      //10  //R7 in Whitesides2010
	addReaction(rxnV3, Reaction(6.700e+11, 8.400e-01, 7.0790e+01, sp::None));      //11 //R8 in Whitesides2010
	addReaction(rxnV3, Reaction(1.300e+11, 1.600e-01, 4.5900e+01, sp::None));      //12 //R9 in Whitesides2010
	
	
	// These are the previous reactions. They do not agree with what the preprints mention.
    /*addReaction(rxnV3, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));   //1 - r1b
    addReaction(rxnV3, Reaction(1.0e10, 0.734, 1.43, sp::OH));  //2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV3, Reaction(2.0e13, 0, 0, sp::H));  //4 - r3f
    addReaction(rxnV3, Reaction(1.1e7, 1.71, 3.9, sp::C2H2));   //5 - r4f1
    addReaction(rxnV3, Reaction(3.3e33, -5.7, 25.5, sp::C2H2));   //6 - r4f2 (r4f = r4f1+r4f2)*/
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //7 - r5f
    m_sType = RFE; // sitetype
    m_name = "R5 conversion to R6 on RFE"; // name of process
    m_ID = 13;
}
// Jump rate calculation
double C5R_RFE::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        double f_R5H = m_r[7]/(m_r[8] + (m_r[0]/r_denom * (m_r[5]+m_r[6]))); 
        r_f *= f_R5H;
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = (m_r[5]+m_r[6])*r_f* site_count; // Rate Equation
}
double C5R_RFE::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double C5R_RFE::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r4f = m_r[5]+m_r[6];
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+r4f);
    double r_f; // radical fraction 
	double f_R5H; //fraction of 5-membered rings with an H2 moeity.
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
		//r4f * r_f gives the rate of acetylene addition to a C on a free edge
		double f_R5H_denom = (m_r[7]/m_r[8] + m_r[9] + m_r[10] + m_r[11] + m_r[12] + r_f * r4f);
		if (f_R5H_denom > 0.0){
			f_R5H = m_r[7] / f_R5H_denom;
		}
		else f_R5H=0;
        //double f_R5H = m_r[0]/(m_r[1] + ((m_r[0]/r_denom) * r4f)); //For previous method. The m_r[0] makes no sense.
        r_f = r_f*f_R5H;
        //r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r4f*r_f* site_count; // Rate Equation
}


// ************************************************************
// ID14- R5 conversion to R6 on AC (AR4 in Matlab)
// ************************************************************
// Elementary rate constants, site type, process type and name
void C5R_RAC::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(8.42e8, 1.49, 0.990, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(6.28e37, -8.24, 44.67, sp::None));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    // r3f defined in setRate0p0267, r4f = 0
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //4 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(8.42e8, 1.49, 0.990, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(1.25e38, -8.20, 45.31, sp::None));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    // r3f defined in setRate0p12, r4f = 0
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //4 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//H addition & migration
	addReaction(rxnV3, Reaction(5.400e+11, 4.500e-01, 1.820e+00, sp::H));          //0  //R12 in Whitesides2010
	addReaction(rxnV3, Reaction(1.791e+00, 0.000e+00, -3.5187e+01, sp::None));      //1  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(2.000e+12, 0.000e+00, 0.000e+00, sp::H));          //2  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(3.100e+11, 8.700e-01, 7.4323e+01, sp::None));      //3  //R7 in Whitesides2010
	addReaction(rxnV3, Reaction(6.700e+11, 8.400e-01, 7.0790e+01, sp::None));      //4  //R8 in Whitesides2010
	addReaction(rxnV3, Reaction(1.300e+11, 1.600e-01, 4.5900e+01, sp::None));      //5  //R9 in Whitesides2010
	
	/*Previous version which does not make sense.
    addReaction(rxnV3, Reaction(8.41e8, 1.49, 0.980, sp::H));      //0 - r1f
    addReaction(rxnV3, Reaction(3.81e11, 0.49, 59.05, sp::None));   //1 - r1b
    addReaction(rxnV3, Reaction(1.0e10, .734, 1.43, sp::OH));  //2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    // r3f defined in setRate1, r4f = 0*/
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //4 - r5f
    m_sType = RAC; // sitetype
    m_name = "R5 conversion to R6 on RAC"; // name of process
    m_ID = 14;
}
// Jump rate calculation
double C5R_RAC::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r3f = 7.297e8 * gp[gp.T] - 5.0641e11; 
    double r_denom = (m_r[1]+m_r[3]+r3f);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = 1.34e12*r_f* site_count; // Rate Equation
}
double C5R_RAC::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r3f = 7.297e8 * gp[gp.T] - 5.0641e11; 
    double r_denom = (m_r[1]+m_r[3]+r3f);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r3f*r_f* site_count; // Rate Equation
}
double C5R_RAC::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
	double r_denom = (m_r[0]/m_r[1] + m_r[2] + m_r[3] + m_r[4] + m_r[5]);
	double r_f;
	if (r_denom > 0) r_f = m_r[0] / r_denom;
	else r_f = 0.0;
	return m_rate = m_r[5] * r_f;
    
	/*Previous version
	// calculate rate
    double r3f = 1.34e12; 
    double r_denom = (m_r[1]+m_r[3]+r3f);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r3f*r_f* site_count; // Rate Equation*/
}


// ************************************************************
// ID15- R5 migration to neighbouring ZZ (AR6 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void M5R_RZZ::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(8.42e8, 1.49, 0.990, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(6.28e37, -8.24, 44.67, sp::None));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    // r3f defined in setRate0p0267, r4f = 0
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //4 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(8.42e8, 1.49, 0.990, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(1.25e38, -8.20, 45.31, sp::None));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    // r3f defined in setRate0p12, r4f = 0
    //addReaction(rxnV2, Reaction(2.20e12, 0, 7.5, sp::O2));          //4 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//H addition & migration
	addReaction(rxnV3, Reaction(5.400e+11, 4.500e-01, 1.820e+00, sp::H));          //0  //R12 in Whitesides2010
	addReaction(rxnV3, Reaction(1.791e+00, 0.000e+00, -3.5187e+01, sp::None));      //1  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(2.000e+12, 0.000e+00, 0.000e+00, sp::H));          //2  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(3.100e+11, 8.700e-01, 7.4323e+01, sp::None));      //3  //R7 in Whitesides2010
	addReaction(rxnV3, Reaction(6.700e+11, 8.400e-01, 7.0790e+01, sp::None));      //4  //R8 in Whitesides2010
	addReaction(rxnV3, Reaction(1.300e+11, 1.600e-01, 4.5900e+01, sp::None));      //5  //R9 in Whitesides2010
	
	/*Old version
    addReaction(rxnV3, Reaction(8.41e8, 1.49, 0.980, sp::H));      //0 - r1f
    addReaction(rxnV3, Reaction(3.81e11, .49, 59.05, sp::None));   //1 - r1b
    addReaction(rxnV3, Reaction(1.0e10, 0.734, 1.43, sp::OH));  //2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    // r3f defined in setRate1, r4f = 0
    //addReaction(rxnV3, Reaction(2.20e12, 0, 7.5, sp::O2));          //4 - r5f*/
    m_sType = RZZ; // sitetype
    m_name = "R5 migration to neighbouring ZZ"; // name of process
    m_ID = 15;
}
// Jump rate calculation
double M5R_RZZ::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r3f = 7.297e8 * gp[gp.T] - 5.0641e11; 
    double r_denom = (m_r[1]+m_r[3]+r3f);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
	return m_rate = 1.34e12*r_f* (site_count); // Rate Equation
}
double M5R_RZZ::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r3f = 7.297e8 * gp[gp.T] - 5.0641e11; 
    double r_denom = (m_r[1]+m_r[3]+r3f);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
	return m_rate = r3f*r_f* (site_count); // Rate Equation
}
double M5R_RZZ::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
	//New version 2019
	double r_denom = (m_r[0]/m_r[1] + m_r[2] + m_r[3] + m_r[4] + m_r[5]);
	double r_f;
	if (r_denom > 0) r_f = m_r[0] / r_denom;
	else r_f = 0.0;
	return m_rate = m_r[5] * r_f * site_count;
	
	/*Old version
    // calculate rate
    double r3f = 1.34e12;
    double r_denom = (m_r[1]+m_r[3]+r3f);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
	if (site_count <= 0){
		return m_rate = 0;
	}
	else{
		return m_rate = r3f*r_f* (site_count); // Rate Equation
	}*/
	
}


// ************************************************************
// ID16- R6 migration & conversion to R5 at BY5 (pyrene+R5; pathway 1; AR22 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void C6R_BY5_FE3::initialise() {
    // Adding elementary reactions
    // 0.0267atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(7.81e7, 1.772, 10.333, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(1.85e6, 1.86, 11.04, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H)); //4 - r3f
    addReaction(rxnV, Reaction(2.3e9, 1.6031, 61.85, sp::None)); //5 - r4f
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    rxnV2 = rxnV;
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    rxnV3 = rxnV;

    m_sType = BY5_FE3; // sitetype
    m_name = "R6 migration & conversion to R5 at BY5 (pathway 1)"; // name of process
    m_ID = 16;
}
// Jump rate calculation
double C6R_BY5_FE3::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation
}
double C6R_BY5_FE3::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double C6R_BY5_FE3::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID17- R6 migration & conversion to R5 at BY5 (pyrene+R5; pathway 2-violi; AR24 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void C6R_BY5_FE3violi::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(7.25e7, 1.76, 9.69, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(1.85e6, 1.86, 11.04, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H)); //4 - r3f
    addReaction(rxnV, Reaction(1.23e10, 1.41, 85.2, sp::None)); //5 - r4f
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    rxnV2 = rxnV;
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    rxnV3 = rxnV;

    m_sType = BY5_FE3; // sitetype
    m_name = "R6 migration & conversion to R5 at BY5 (pathway 2; violi)"; // name of process
    m_ID = 17;
}
// Jump rate calculation
double C6R_BY5_FE3violi::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation
}
double C6R_BY5_FE3violi::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double C6R_BY5_FE3violi::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
// ************************************************************
// ID18- BY5 closure (AR16 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
bool violi = false;
void L5R_BY5::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(1.73e7, 1.7453, 8.732, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.40e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H)); //4 - r3f
    if(violi) {
        addReaction(rxnV, Reaction(3.86e11, .21, 17.7, sp::None)); //5 - r4f - violi(2005) mechanism, include reverse.
    }
    else {
        addReaction(rxnV, Reaction(1.6451e12, .3239, 15.088227, sp::None)); // 5 - r4f - AR's work
        addReaction(rxnV, Reaction(9.9866e13, 0, 22.57, sp::None)); //6
        addReaction(rxnV, Reaction(5.75e10, .93, 30.4, sp::None));//7 - for rate calculation
    }
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(1.73e7, 1.7453, 8.732, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.40e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H)); //4 - r3f
    if(violi)
        addReaction(rxnV2, Reaction(7.72e11, .21, 17.7, sp::None)); //5 - r4f - violi(2005) mechanism, include reverse.
    else {
        addReaction(rxnV2, Reaction(3.2902e12, .3239, 15.088227, sp::None)); // 5 - r4f - AR's work
        addReaction(rxnV2, Reaction(9.9866e13, 0, 22.560985, sp::None)); //6
        addReaction(rxnV2, Reaction(5.75e10, .93, 30.4, sp::None));//7 - for rate calculation
    }
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //6
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
    /*addReaction(rxnV3, Reaction(7.25e7, 1.76, 9.69, sp::H));      //0 - r1f
    addReaction(rxnV3, Reaction(3.40e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV3, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV3, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV3, Reaction(3.49e39, -7.77, 13.37320574, sp::H)); //4 - r3f
    if(violi)
        addReaction(rxnV3, Reaction(3.86e11, .21, 17.7, sp::None)); //5 - r4f - violi(2005) mechanism, include reverse.
    else {
        addReaction(rxnV3, Reaction(1.6451e12, .3239, 15.088227, sp::None)); // 5 - r4f - AR's work
        addReaction(rxnV3, Reaction(9.9866e13, 0, 22.560985, sp::None)); //6
        addReaction(rxnV3, Reaction(5.75e10, .93, 30.4, sp::None));//7 - for rate calculation
    }
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //6*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.534E+07,1.854E+00,8.443E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.311E+05,1.977E+00,1.220E+01, sp::H2)); //-1
	addReaction(rxnV3, Reaction(1.485E+11,3.012E-01,1.732E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(8.560E+11,4.399E-01,2.480E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(7.973E+09,1.184E+00,2.938E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(1.611E+08,1.536E+00,1.781E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward

    m_sType = BY5; // sitetype
    m_name = "BY5 closure"; // name of process
    m_ID = 18;
}
// Jump rate calculation
double L5R_BY5::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        if(violi)
            r_f = (m_r[0]+m_r[2])/r_denom; 
        else
            r_f = (m_r[0]+m_r[2])/(r_denom+((m_r[6]/m_r[7])*(m_r[1]+m_r[4]))); 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = 2*m_r[5]*r_f* site_count; // Rate Equation
}
double L5R_BY5::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double L5R_BY5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    //return setRate0p0267(gp, pah_st); // OLD METHOD
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(2, 2);
	boost::numeric::ublas::vector<double> arr2(2);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR5-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = m_r[4] * arr2(1) * site_count;
}
// ************************************************************
// ID19- R6 desorption at bay -> pyrene (AR21 in Matlab)
// ************************************************************

// Elementary rate constants, site type, process type and name
void M6R_BY5_FE3::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(7.81e7, 1.772, 10.333, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(1.85e6, 1.86, 11.04, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H)); //4 - r3f
    addReaction(rxnV, Reaction(2.3e9, 1.6031, 61.85, sp::None)); //5 - r4f
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    rxnV2 = rxnV;
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    rxnV3 = rxnV;

    m_sType = BY5_FE3; // sitetype
    m_name = "R6 desorption at bay BY5_FE3"; // name of process
    m_ID = 19;
}
// Jump rate calculation
double M6R_BY5_FE3::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation
}
double M6R_BY5_FE3::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double M6R_BY5_FE3::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID20- R6 oxidation at FE2 (side)
// ************************************************************

// Elementary rate constants, site type, process type and name
void O6R_FE2_side::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(3.23e7, 2.095, 15.84, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(1.3e11, 1.08, 70.42, sp::OH));   //5 - r4f
    addReaction(rxnV, Reaction(1.3e13, 0, 10.62, sp::OH));          //6 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    //addReaction(rxnV2, Reaction(4.2e13, 0, 13, sp::H));      //0 - r1f
    //addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    //addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    //addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    //addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(1.3e13, 0, 10.62, sp::OH));   //5 - r4f
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::OH));          //6 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4
	addReaction(rxnV3, Reaction(3.17e13, 0.00, 2.021,sp::O2));  						// S46                              - 5
    addReaction(rxnV3, Reaction(8.82e9, -0.063, 16612.2*8.314/4.184/1000.0,sp::None));  // S51                                - 6 //Using reaction from He Lin et al @1atm
    addReaction(rxnV3, Reaction(1.47e14, 0.00, 632.0*8.314/4.184/1000.0,sp::OH)); // S55                                    - 7
    addReaction(rxnV3, Reaction(2.13e15, 0.00, 42305.0*8.314/4.184/1000.0,sp::None));  // S60                               - 8
    addReaction(rxnV3, Reaction(4.34e14, 0.00, 984.0*8.314/4.184/1000.0,sp::H));  // S65                                    - 9
    addReaction(rxnV3, Reaction(2.14e16, 0.00, 52988.0*8.314/4.184/1000.0,sp::None));  // S70                               - 10
    addReaction(rxnV3, Reaction(5476782.795, 2.147, 5603.0*8.314/4.184/1000.0,sp::H));  // S75                               - 11 // Correct n value is 2.147
    addReaction(rxnV3, Reaction(1.70e14, 0.00, 9635.0*8.314/4.184/1000.0,sp::H2));  // S80                                  - 12
    addReaction(rxnV3, Reaction(2.00e14, 0.00, 2670.0*8.314/4.184/1000.0,sp::H));  // S85                                   - 13
    addReaction(rxnV3, Reaction(4.00e12, 0.00, 2328.0*8.314/4.184/1000.0,sp::O));  // S93                                   - 14
    addReaction(rxnV3, Reaction(0.0*1.00e14, 0.00, 0.0*8.314/4.184/1000.0,sp::OH));  // S96                                     - 15
    addReaction(rxnV3, Reaction(0.0*1.00e10, 0.00, 0.0*8.314/4.184/1000.0,sp::H2O));  // S100    								- 16	//It is 1.00e12 in Frenkalch's paper but that's too high. Lowered here.
	addReaction(rxnV3, Reaction(1.76e23, -3.681, 24155.0*8.314/4.184/1000.0,sp::None));  // S52                               - 17	//Using reaction from He Lin et al @1atm
	addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 18              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 19              - Frenklach et al. 2018
	
	//He Lin et al reference :
	//New Insights into Thermal Decomposition of Polycyclic Aromatic Hydrocarbon Oxyradicals. Peng Liu, He Lin,* Yang Yang, Can Shao, Chen Gu, and Zhen Huang 2014
    //addReaction(rxnV3, Reaction(1.3e13, 0, 10.62, sp::OH));   //5 - r4f
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
    m_sType = FE2; // sitetype
    m_name = "R6 (FE2) Oxidation side"; // name of process
    m_ID = 20;
}
// Jump rate calculation
double O6R_FE2_side::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
	/*// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;

    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[6]*r_f* site_count; // Rate Equation*/
}
double O6R_FE2_side::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
    /*// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;

    //double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    //double r_f; // radical fraction 
    //if(r_denom>0) {
    //    r_f = (m_r[0]+m_r[2])/r_denom; 
    //    r_f = r_f/(r_f+1.0);
    //}
    //else r_f=0;
	return m_rate = m_r[0] * site_count; // Rate Equation*/
}
double O6R_FE2_side::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	if ((int)site_count % 2 != 0) std::cout << "Error. FE2 number of sites is " << site_count << " and is not even.\n";
    site_count = site_count / 2.0;
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	//Intermediaries are PAH*, PAH-O, PAH-OH
	arr1(0,0) = m_r[1] + m_r[3] + m_r[4] + m_r[5] + m_r[7]+m_r[15]+m_r[16]+m_r[18]+m_r[19];
    arr1(0,2) = -m_r[10] - m_r[13];
    arr1(1,0) = -m_r[5];
    arr1(1,1) = m_r[17] + m_r[9] + m_r[12];
    arr1(1,2) = - m_r[8] - m_r[11];
	arr1(2,0) = -m_r[7];
    arr1(2,1) = - m_r[9] - m_r[12];
    arr1(2,2) = m_r[8] + m_r[10] + m_r[11] + m_r[13];

    arr2(0) = +m_r[0] + m_r[2];
    arr2(1) = +m_r[14];
    arr2(2) = 0.0;

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	double Crad_fraction = arr2(0);
	double CO_fraction = arr2(1);
	double COH_fraction = arr2(2);
	if (CO_fraction >= 1.0 || COH_fraction >= 1.0) return m_rate = 0.0;
	else return m_rate = (m_r[17] * CO_fraction + (m_r[15] + m_r[16]) * Crad_fraction)*site_count;
	/*if (arr2(1) < arr2(0) && arr2(2) < arr2(0)) return m_rate = (m_r[17] * arr2(1) + m_r[16] * arr2(0))*site_count;
	else return m_rate = 0.0; // The SS approximation breaks so the rate approximation is wrong.*/
}

// ************************************************************
// ID21- R6 oxidation at FE2 (top)
// ************************************************************
// Elementary rate constants, site type, process type and name
void O6R_FE2_top::initialise() {
    // Adding elementary reactions
    // 0.0267atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(3.23e7, 2.095, 15.84, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));   //5 - r5f(1)
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f(2)
    // 0.12atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
	addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));   //5 - r5f(1)
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.51674641, sp::O2));          //6 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4
	addReaction(rxnV3, Reaction(3.17e13, 0.00, 2.021,sp::O2));  						// S46                              - 5
    addReaction(rxnV3, Reaction(8.82e9, -0.063, 16612.2*8.314/4.184/1000.0,sp::None));  // S51                              - 6 //Using reaction from He Lin et al @1atm
    addReaction(rxnV3, Reaction(1.47e14, 0.00, 632.0*8.314/4.184/1000.0,sp::OH)); // S55                                    - 7
    addReaction(rxnV3, Reaction(2.13e15, 0.00, 42305.0*8.314/4.184/1000.0,sp::None));  // S60                               - 8
    addReaction(rxnV3, Reaction(4.34e14, 0.00, 984.0*8.314/4.184/1000.0,sp::H));  // S65                                    - 9
    addReaction(rxnV3, Reaction(2.14e16, 0.00, 52988.0*8.314/4.184/1000.0,sp::None));  // S70                               - 10
    addReaction(rxnV3, Reaction(5476782.795, 2.147, 5603.0*8.314/4.184/1000.0,sp::H));  // S75                              - 11 // Correct n value is 2.147
    addReaction(rxnV3, Reaction(1.70e14, 0.00, 9635.0*8.314/4.184/1000.0,sp::H2));  // S80                                  - 12
    addReaction(rxnV3, Reaction(2.00e14, 0.00, 2670.0*8.314/4.184/1000.0,sp::H));  // S85                                   - 13
    addReaction(rxnV3, Reaction(4.00e12, 0.00, 2328.0*8.314/4.184/1000.0,sp::O));  // S93                                   - 14
    addReaction(rxnV3, Reaction(0.0*1.00e14, 0.00, 0.0*8.314/4.184/1000.0,sp::OH));  // S96                                     - 15
    addReaction(rxnV3, Reaction(0.0*1.00e10, 0.00, 0.0*8.314/4.184/1000.0,sp::H2O));  // S100									- 16	//It is 1.00e12 in Frenkalch's paper but that's too high. Lowered here.
	addReaction(rxnV3, Reaction(1.76e23, -3.681, 24155.0*8.314/4.184/1000.0,sp::None));  // S52                               - 17	//Using reaction from He Lin et al @1atm
	addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 18              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 19              - Frenklach et al. 2018
	
	//He Lin et al reference :
	//New Insights into Thermal Decomposition of Polycyclic Aromatic Hydrocarbon Oxyradicals. Peng Liu, He Lin,* Yang Yang, Can Shao, Chen Gu, and Zhen Huang 2014
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));   //5 - r4f
    //addReaction(rxnV3, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
    m_sType = FE2; // sitetype
    m_name = "R6 (FE2) Oxidation top"; // name of process
    m_ID = 21;
}
// Jump rate calculation
double O6R_FE2_top::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
    /*// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[5]*r_f* site_count; // Rate Equation*/
}
double O6R_FE2_top::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
	/*// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
	return m_rate = m_r[5] * r_f* site_count; // Rate Equation*/
}
double O6R_FE2_top::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	if ((int)site_count % 2 != 0) std::cout << "Error. FE2 number of sites is " << site_count << " and is not even.\n";
    site_count = site_count / 2.0;
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	//Intermediaries are PAH*, PAH-O, PAH-OH
	arr1(0,0) = m_r[1] + m_r[3] + m_r[4] + m_r[5] + m_r[7]+m_r[15]+m_r[16]+m_r[18]+m_r[19];
    arr1(0,2) = -m_r[10] - m_r[13];
    arr1(1,0) = -m_r[5];
    arr1(1,1) = m_r[6] + m_r[9] + m_r[12];
    arr1(1,2) = - m_r[8] - m_r[11];
	arr1(2,0) = -m_r[7];
    arr1(2,1) = - m_r[9] - m_r[12];
    arr1(2,2) = m_r[8] + m_r[10] + m_r[11] + m_r[13];

    arr2(0) = +m_r[0] + m_r[2];
    arr2(1) = +m_r[14];
    arr2(2) = 0.0;

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	double Crad_fraction = arr2(0);
	double CO_fraction = arr2(1);
	double COH_fraction = arr2(2);
	if (CO_fraction >= 1.0 || COH_fraction >= 1.0) return m_rate = 0.0;
	else return m_rate = (m_r[6] * CO_fraction + (m_r[15] + m_r[16]) * Crad_fraction)*site_count;
	/*if (arr2(1) < arr2(0) && arr2(2) < arr2(0)) return m_rate = (m_r[6] * arr2(1) + m_r[16] * arr2(0))*site_count;
	else return m_rate = 0.0; // The SS approximation breaks so the rate approximation is wrong.*/
}

// Elementary rate constants, site type, process type and name
void D6R_FE_AC::initialise() {
	// Adding elementary reactions
	// 0.0267 atm
	rxnvector& rxnV = m_rxnvector0p0267;
	addReaction(rxnV, Reaction(3.23e7, 2.095, 15.84, sp::H));      //0 - r1f
	addReaction(rxnV, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
	addReaction(rxnV, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
	addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
	addReaction(rxnV, Reaction(3.49e39, -7.77, 13.37320574, sp::H));  //4 - r3f
	addReaction(rxnV, Reaction(1.5e11, 0.2, 42.57, sp::None));   //5 - r4f
	//addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f
	// 0.12 atm
	rxnvector& rxnV2 = m_rxnvector0p12;
	addReaction(rxnV2, Reaction(4.2e13, 0, 13.00, sp::H));      //0 - r1f
	addReaction(rxnV2, Reaction(3.4e9, .88, 7.870, sp::H2));   //1 - r1b
	addReaction(rxnV2, Reaction(2.1e13, 0, 4.56937799, sp::OH));  //2 - r2f
	addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
	addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
	addReaction(rxnV2, Reaction(2e11, 0.2, 42.57, sp::None));   //5 - r4f
	//addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//Reverse reactions fitted from Blanquart thermal data with Frenklach rates from 2018
	addReaction(rxnV3, Reaction(5.465e+30, -3.657e+00, 8.624e+01, sp::H));          			// A3- + C2H2 <=> A4 + H                - 0             - Backward
	addReaction(rxnV3, Reaction(4.868e+22, -1.697e+00, 7.555e+01, sp::H));          			// A3- + C2H2 <=> A4 + H  
	addReaction(rxnV3, Reaction(1.815E+34, -4.875E+00, 7.521E+01, sp::H));           			// A3-4 + C2H2 <=> A4 + H              - 12             - Backward

	//Reactions for ABF
	/*//------------Reactions for A3-4----------------------------
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));            	// A3 + H <=> A3-4 + H2              	- 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));           	// A3 + H <=> A3-4 + H2              	- 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           	// A3 + OH <=> A3-4 + H2O              	- 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          	// A3 + OH <=> A3-4 + H2O              	- 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            	// A3- + H -> A3              			- 4              - Forward
	//addReaction(rxnV3, Reaction(1.000E+14, 0.000E+00, 0.000E+00, sp::None));         	// A3-4 + H (+M) <=> A3 (+M)           	- 4              - Forward - High rate
	//addReaction(rxnV3, Reaction(2.100E+139, -3.480E+01, 1.838E+01, sp::None));    	// A3-4 + H (+M) <=> A3 (+M)          	- 5              - Forward - High rate
	//Parameters = [  1.00000000e-03   1.71400000e+02   1.71400000e+02   4.99280000e+03]
	addReaction(rxnV3, Reaction(2.315E+31, -4.589E+00, 1.185E+02, sp::None));      		// A3-4 + H (+M) <=> A3 (+M)           - 6              - Backward
	addReaction(rxnV3, Reaction(8.000E+61, -1.450E+01, 3.480E+01, sp::C2H2));      		// A3-4 + C2H2 <=> A3C2H2              - 7              - Forward
	addReaction(rxnV3, Reaction(4.107E+69, -1.626E+01, 6.773E+01, sp::None));        	// A3-4 + C2H2 <=> A3C2H2              - 8              - Backward
	addReaction(rxnV3, Reaction(1.200E+26, -3.440E+00, 3.020E+01, sp::C2H2));         	// A3-4 + C2H2 <=> A3C2H + H           - 9              - Forward
	addReaction(rxnV3, Reaction(6.319E+33, -5.219E+00, 3.169E+01, sp::H));           	// A3-4 + C2H2 <=> A3C2H + H           - 10             - Backward
	addReaction(rxnV3, Reaction(6.600E+24, -3.360E+00, 1.780E+01, sp::C2H2));      		// A3-4 + C2H2 <=> A4 + H              - 11             - Forward
	addReaction(rxnV3, Reaction(1.815E+34, -4.875E+00, 7.521E+01, sp::H));        		// A3-4 + C2H2 <=> A4 + H              - 12             - Backward
	////------------Reactions for A3C2H----------------------------
	addReaction(rxnV3, Reaction(5.000E+13, 0.000E+00, 0.000E+00, sp::C2H));          	// A3 + C2H <=> A3C2H + H              - 13             - Forward
	addReaction(rxnV3, Reaction(9.332E+20, -1.526E+00, 2.401E+01, sp::H));           	// A3 + C2H <=> A3C2H + H              - 14             - Backward
	addReaction(rxnV3, Reaction(1.900E+64, -1.512E+01, 2.930E+01, sp::H));           	// A3C2H + H <=> A3C2H2                - 15             - Forward
	addReaction(rxnV3, Reaction(1.852E+64, -1.510E+01, 6.074E+01, sp::None));        	// A3C2H + H <=> A3C2H2                - 16             - Backward
	////------------Reactions for A3C2H2----------------------------
	addReaction(rxnV3, Reaction(2.000E+63, -1.528E+01, 4.320E+01, sp::None));           // A3C2H2 <=> A4 + H           - 17             - Forward
	addReaction(rxnV3, Reaction(1.071E+65, -1.504E+01, 6.768E+01, sp::H));           // A3C2H2 <=> A4 + H           - 18             - Backward
	////------------Reactions for A4----------------------------
	addReaction(rxnV3, Reaction(9.000E+38, -7.390E+00, 2.070E+01, sp::H));           // A3C2H + H <=> A4 + H                - 19             - Forward
	addReaction(rxnV3, Reaction(4.700E+40, -7.125E+00, 7.662E+01, sp::H));           // A3C2H + H <=> A4 + H                - 20             - Backward
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A4 + H <=> A4- + H2         - 21             - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A4 + H <=> A4- + H2         - 22             - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A4 + OH <=> A4- + H2O               - 23             - Forward
	addReaction(rxnV3, Reaction(9.600E+04, 2.162E+00, 9.780E+00, sp::H2O));          // A4 + OH <=> A4- + H2O               - 24             - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A4- + H -> A4              			- 25              - Forward
	//addReaction(rxnV3, Reaction(1.000E+14, 0.000E+00, 0.000E+00, sp::H));          // A4- + H <=> A4              - 25             - Forward
	//addReaction(rxnV3, Reaction(4.635E+18, -7.927E-01, 1.113E+02, sp::None));                // A4- + H <=> A4              - 26             - Backward
	////------------Reactions for A4-----------------------------	*/

	m_sType = FE_HACA; // sitetype
	m_name = "R6 (FE_HACA) desorption"; // name of process
	m_ID = 22;
}
// Jump rate calculation
double D6R_FE_AC::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
		r_f = r_f / (r_f + 1.0);
	}
	else r_f = 0;
	return m_rate = m_r[5] * r_f* site_count; // Rate Equation
}
double D6R_FE_AC::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate0p0267(gp, pah_st);
}
double D6R_FE_AC::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	//int ring_count = ((int)pah_st.getRingsCount().first); // Ring count
	//auto ring_count = pah_st.getRingsCount();
	//int ring_count_6 = std::get<0>(ring_count);
	//if (ring_count_6 <= 4) return m_rate = 0.0;
	double T = gp[gp.T];
	double r_f = 0.0;
	if (T < 600) {
		return m_rate = 0.0;
	}
	else{
		//For Frenklach2018
		return m_rate = (m_r[0] + m_r[1] + m_r[2])* site_count / 2.0; //Lumped consumption of A4
		//For ABF
		//return m_rate = (m_r[12] + m_r[20] + m_r[18]) * site_count / 2.0; //Lumped consumption of A4
	}
}

// ************************************************************
// ID22 - Bay-capping
// ************************************************************

void B6R_ACR5::initialise() {
	// Adding elementary reactions
	// 0.0267 atm
	rxnvector& rxnV = m_rxnvector0p0267;

	/**
	* Reaction mechanism: The units are mole, centimeter, second, Kelbins, and kilocalorie.
	*
	* C6H6 + H
	*
	* A. M. Mebel, M. C. Lin, T. Yu, and K. Morokuma. Theoretical study of potential
	* energy surface and thermal rate constants for the C6H5 + H2 and C6H6 + H reactions.
	* J. Phys. Chem. A, 101:3189�3196, 1997. doi:10.1021/jp9702356.
	*/
	addReaction(rxnV, Reaction(3.23e7, 2.095, 15.842, sp::H));

	/**
	* Cs* + H2 => Cs-H + H
	*
	* M. Frenklach. On surface growth mechanism of soot particles. Proc. Combust. Inst.,
	* 26:2285�2293, 1996. doi:10.1016/S0082-0784(96)80056-7.
	*/
	addReaction(rxnV, Reaction(3.4e9, 0.88, 7.86, sp::H2));

	/**
	* A1- + H = A1
	* A1- + C2H2 = A1C2H + H
	*
	* M. Frenklach and H. Wang. Detailed surface and gas-phase chemi-
	* cal kinetics of diamond deposition. Phys. Rev. B, 43:1520�1545, 1991.
	* doi:10.1103/PhysRevB.43.1520.
	*/
	addReaction(rxnV, Reaction(3.48e39, -7.77, 13.36, sp::H));
	addReaction(rxnV, Reaction(4.00e13, 0, 10.11, sp::C2H2));

	//
	// 0.12 atm
	rxnvector& rxnV2 = m_rxnvector0p12;
	rxnV2 = rxnV;

	rxnvector& rxnV3 = m_rxnvector1;
	//Rates taken from Raj 2019. Combustion and Flame 204 (2019) 331–340
	addReaction(rxnV3, Reaction(2.540E+11,9.310E-01,1.644E+01, sp::H)); //1
	addReaction(rxnV3, Reaction(1.830E+12,3.970E-01,8.815E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(1.630E+12,4.090E-01,5.975E+00, sp::C2H2)); //2
	addReaction(rxnV3, Reaction(9.130E+11,9.910E-01,4.599E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(6.320E+11,1.660E-01,1.805E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(9.750E+10,4.580E-01,1.583E+01, sp::None)); //-3
	addReaction(rxnV3, Reaction(9.580E+11,-6.400E-02,1.631E+01, sp::None)); //4
	addReaction(rxnV3, Reaction(9.650E+11,5.010E-01,4.150E+01, sp::None)); //-4
	addReaction(rxnV3, Reaction(3.160E+12,7.870E-01,3.651E+01, sp::None)); //5
	addReaction(rxnV3, Reaction(9.710E+11,5.070E-01,4.695E+00, sp::H)); //-5
	addReaction(rxnV3, Reaction(2.780E+11,6.300E-02,2.387E+01, sp::None)); //6
	addReaction(rxnV3, Reaction(5.470E+11,6.450E-01,3.277E+01, sp::None)); //-6
	addReaction(rxnV3, Reaction(8.150E+11,5.630E-01,2.486E+01, sp::None)); //7
	addReaction(rxnV3, Reaction(9.060E+11,4.560E-01,7.286E+00, sp::H)); //-7
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4              - Forward
	addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 7              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 8              - Frenklach et al. 2018


	//OLD
	/*addReaction(rxnV3, Reaction(4.20e13, 0, 13.00, sp::H));         // 0 - r1f
	addReaction(rxnV3, Reaction(5.53e12, 0, 11.83, sp::H2));        // 1 - rlb
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));             // 4 - r3b
	addReaction(rxnV3, Reaction(1.00e10, 0.955, 10.58, sp::C2H2));  // 5 - r4*/

	m_sType = ACR5; // sitetype
	m_name = "Bay-capping"; // name of process
	m_ID = 23;
}
// Jump rate calculation
double B6R_ACR5::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[2] + m_r[3]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = m_r[0] / r_denom;
	}
	else r_f = 0;
	return m_rate = m_r[3] * r_f*site_count; // Rate Equation
}
double B6R_ACR5::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate0p0267(gp, pah_st);
}
double B6R_ACR5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	matrix<double> arr1(5, 5);
	boost::numeric::ublas::vector<double> arr2(5);
	//matrix<double> arr1_peq(7, 7);
	//boost::numeric::ublas::vector<double> arr2_peq(7);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	/*for (unsigned k = 0; k < arr1_peq.size1(); ++k)
		for (unsigned l = 0; l < arr1_peq.size2(); ++l)
			arr1_peq(k, l) = 0.0;
	for (unsigned k = 0; k < arr2_peq.size(); ++k)
		arr2_peq(k) = 0.0;*/
	arr1(0,0) = m_r[1] + m_r[2] + m_r[14] + m_r[15] + m_r[16];
	arr1(0,1) = -m_r[3];
	arr1(1,0) = -m_r[2];
	arr1(1,1) = m_r[3] + m_r[4] + m_r[10];
	arr1(1,2) = -m_r[5];
	arr1(1,4) = -m_r[11];
	arr1(2,1) = -m_r[4];
	arr1(2,2) = m_r[5] + m_r[6];
	arr1(2,3) = -m_r[7];
	arr1(3,2) = -m_r[6];
	arr1(3,3) = m_r[8] + m_r[7];
	arr1(4,1) = -m_r[10];
	arr1(4,4) = m_r[11] + m_r[12];
	
	arr2(0) = m_r[0];
	
	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	return m_rate = (m_r[8] * arr2(3) + m_r[12] * arr2(4))*site_count;
	
	//OLD
	/*double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f*site_count*0.1; // Rate Equation slowed by 10*/
}

// ************************************************************
// ID23 - Embedded 5-member ring migration to ZZ
// ************************************************************

void M5R_ACR5_ZZ::initialise() {
	// Adding elementary reactions
	// 0.0267 atm
	rxnvector& rxnV = m_rxnvector0p0267;
	//------------Reactions for A3-4----------------------------
	addReaction(rxnV, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4   					- 4            - Forward
	addReaction(rxnV, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 5              - Frenklach et al. 2018
	addReaction(rxnV, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));        // A3* + C2H2 -> A3C2H + H        	- 6              - Frenklach et al. 2018
	addReaction(rxnV, Reaction(4.960e+11, 7.550e-01, 5.000e+01, sp::None));  		// Violi2005. 							- 7
	addReaction(rxnV, Reaction(5.188e+13, 0.000E+00, 3.340e+01, sp::None));  		// Whitesides2009 system b. 			- 8
	//Assumes the migration of embedded rings is similar to migration to the edge	
	
	//ABF and old
	/*addReaction(rxnV, Reaction(1.74e08, 1.740, 9.370, sp::H));     // 0 - r1f
	addReaction(rxnV, Reaction(3.90e12, 0, 11.00, sp::H2));    // 1 - r1b
	addReaction(rxnV, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2
	addReaction(rxnV, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV, Reaction(4.96e11, 0.755, 50, sp::None));  // 5*/
	//
	// 0.12 atm
	rxnvector& rxnV2 = m_rxnvector0p12;
	rxnV2 = rxnV;
	//
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	rxnV3 = rxnV;
	//
	m_sType = MIGR; // sitetype
	m_name = "ACR5 migration"; // name of process
	m_ID = 24;
}
// Jump rate calculation
double M5R_ACR5_ZZ::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5] + m_r[6] + m_r[7]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	//return m_rate = m_r[7] * r_f*site_count; // Rate Equation //Violi
	return m_rate = m_r[8] * r_f*site_count; // Rate Equation //Whitesides
}
double M5R_ACR5_ZZ::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate0p0267(gp, pah_st);
}
double M5R_ACR5_ZZ::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID24 - R6 growth on RZZ
// ************************************************************
// Elementary rate constants, site type, process type and name
void G6R_RZZ::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//Checked 2019
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 6 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 7 Hou2017
	
	//Previous. No idea where they got the rates from.
	/*addReaction(rxnV3, Reaction(4.20e13, 0, 13.00, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.90e12, 0, 11.00, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(8.00e07, 1.560, 3.800, sp::C2H2));  // 5 - r4f*/

	m_sType = RZZ; // sitetype
	m_name = "G6R at RZZ"; // name of process
	m_ID = 25;
}

double G6R_RZZ::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_denom2 = (m_r[7] + m_r[3] + m_r[4] + m_r[5]);
	double r_f, r_f2; // radical fraction
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	if (r_denom2>0) {
		r_f2 = (m_r[6] + m_r[2]) / r_denom2;
	}
	else r_f2 = 0;
	return m_rate = m_r[5] * (r_f + r_f2) * site_count; // Rate Equation
}

// ************************************************************
// ID25 - R6 growth on RFER
// ************************************************************
// Elementary rate constants, site type, process type and name
void G6R_RFER::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//Checked 2019
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010. Combination of S17 and S39
	
	//Previous. No idea about the C2H2 rate
	/*addReaction(rxnV3, Reaction(4.20e13, 0, 13.00, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.90e12, 0, 11.00, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(8.00e07, 1.560, 3.800, sp::C2H2));  // 5 - r4f*/

	m_sType = RFER; // sitetype
	m_name = "G6R at RFER"; // name of process
	m_ID = 26;
}

double G6R_RFER::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation
}

// ************************************************************
// ID26 - R6 growth on R5
// ************************************************************
// Elementary rate constants, site type, process type and name
void G6R_R5::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	addReaction(rxnV3, Reaction(4.20e13, 0, 13.00, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.90e12, 0, 11.00, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.10e07, 1.610, 3.896, sp::C2H2));  // 5 - r4f1

	m_sType = R5; // sitetype
	m_name = "G6R at R5"; // name of process
	m_ID = 27;
}

double G6R_R5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = m_r[5] * r_f* site_count; // Rate Equation
}

// ************************************************************
// ID27- RBY5 closure reaction
// ************************************************************
// Elementary rate constants, site type, process type and name
void L6_RBY5::initialise() {
	// Adding elementary reactions
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
	/*addReaction(rxnV3, Reaction(9.24e07, 1.500, 9.646, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(9.60e04, 1.960, 9.021, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.11e11, 0.658, 23.99, sp::None));  // 5 - r4f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

	m_sType = RBY5; // sitetype
	m_name = "RBY5 closure"; // name of process
	m_ID = 28;
}
double L6_RBY5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	/*// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}

// ************************************************************
// ID28- RACR closure reaction
// ************************************************************
// Elementary rate constants, site type, process type and name
void L6_RACR::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
	/*addReaction(rxnV3, Reaction(9.24e07, 1.500, 9.646, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(9.60e04, 1.960, 9.021, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.11e11, 0.658, 23.99, sp::None));  // 5 - r4f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

	m_sType = RACR; // sitetype
	m_name = "RACR closure"; // name of process
	m_ID = 29;
}
double L6_RACR::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	/*// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}

// ************************************************************
// ID29 - R5 growth on RFE
// ************************************************************
// Elementary rate constants, site type, process type and name
void G5R_RFE::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	addReaction(rxnV3, Reaction(4.20e13, 0, 13.00, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.90e12, 0, 11.00, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(6.80e11, 0, 22.02, sp::C2H2));  // 5 - r4f

	m_sType = RFE; // sitetype
	m_name = "R5 growth on RFE"; // name of process
	m_ID = 30;
}
double G5R_RFE::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = m_r[5] * r_f* site_count; // Rate Equation
}

// ************************************************************
// ID30 - R6 migration & conversion to R5 at RAC
// ************************************************************
// Elementary rate constants, site type, process type and name
void C6R_RAC_FE3::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	addReaction(rxnV3, Reaction(7.81e07, 1.772, 10.33, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));         // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(2.30e09, 1.603, 61.85, sp::None));  // 5 - r4f

	m_sType = RAC_FE3; // sitetype
	m_name = "R6 migration & conversion to R5 at RAC_FE3"; // name of process
	m_ID = 31;
}

double C6R_RAC_FE3::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = m_r[5] * r_f* site_count; // Rate Equation
}

// ************************************************************
// ID31 - R6 migration & conversion to R5 at RAC
// ************************************************************
// Elementary rate constants, site type, process type and name
void C6R_RAC_FE3violi::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	addReaction(rxnV3, Reaction(7.25e07, 1.772, 10.33, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));         // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.23e10, 1.410, 85.20, sp::None));  // 5 - r4f

	m_sType = RAC_FE3; // sitetype
	m_name = "R6 migration & conversion to R5 at RAC_FE3"; // name of process
	m_ID = 32;
}

double C6R_RAC_FE3violi::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = m_r[5] * r_f* site_count; // Rate Equation
}

// ************************************************************
// ID32 - R6 desorption at RAC -> pyrene
// ************************************************************

// Elementary rate constants, site type, process type and name
void M6R_RAC_FE3::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	addReaction(rxnV3, Reaction(7.81e07, 1.772, 10.33, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));         // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(2.30e09, 1.603, 61.85, sp::None));  // 5 - r4f

	m_sType = RAC_FE3; // sitetype
	m_name = "R6 desorption at RAC_FE3"; // name of process
	m_ID = 33;
}

double M6R_RAC_FE3::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = m_r[5] * r_f* site_count; // Rate Equation
}
// ************************************************************
// ID34 - R5 exchange with R6
// ************************************************************

// Elementary rate constants, site type, process type and name
void MR5_R6::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//------------Reactions for A3-4----------------------------
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4              - Forward
	//Using rates by On the low-temperature limit of HACA. Frenklach et al. 2018
	addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 5              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));        // A3* + C2H2 -> A3C2H + H        	- 6              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(4.960e+11, 7.550e-01, 5.000e+01, sp::None));  		// Violi2005. 							- 7
	addReaction(rxnV3, Reaction(5.188e+13, 0.000E+00, 3.340e+01, sp::None));  		// Whitesides2009 system b. 			- 8
	//Assumes the migration of embedded rings is similar to migration to the edge	
	
	//H addition & migration
	/*addReaction(rxnV3, Reaction(5.400e+11, 4.500e-01, 1.820e+00, sp::H));          //0  //R12 in Whitesides2010
	addReaction(rxnV3, Reaction(1.791e+00, 0.000e+00, -3.5187e+01, sp::None));      //1  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(2.000e+12, 0.000e+00, 0.000e+00, sp::H));          //2  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(3.100e+11, 8.700e-01, 7.4323e+01, sp::None));      //3  //R7 in Whitesides2010
	addReaction(rxnV3, Reaction(6.700e+11, 8.400e-01, 7.0790e+01, sp::None));      //4  //R8 in Whitesides2010
	addReaction(rxnV3, Reaction(1.300e+11, 1.600e-01, 4.5900e+01, sp::None));      //5  //R9 in Whitesides2010
	addReaction(rxnV3, Reaction(8.900e+05, 2.280e+00, 6.1489e+01, sp::None));      //6  //R14 in Whitesides2010*/
	
	//Old
	/*addReaction(rxnV3, Reaction(7.81e07, 1.772, 10.33, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.9e12, 0, 11.00, sp::H2));         // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f*/

	m_sType = R5R6_MIGR; // sitetype
	m_name = "R5R6 migration"; // name of process
	m_ID = 34;
}

double MR5_R6::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = m_r[1] + m_r[3] + m_r[4] + m_r[5] + m_r[6] + m_r[7];
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	//return m_rate = m_r[7] * r_f* site_count / 2.0; // Rate Equation // Violi
	return m_rate = m_r[8] * r_f* site_count * 2.0; // Rate Equation // Whitesides
}

// ************************************************************
// ID35 - R7 growth on embedded-obstructed R5
// ************************************************************

// Elementary rate constants, site type, process type and name
void GR7_R5R6AC::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//Reactions provided by Menon2019 at M06-2X level.
	//Pathway A
	addReaction(rxnV3, Reaction(6.586e+07, 1.766e+00, 1.477e+01, sp::H));           //7-HACAsite2-1 + H => 7-HACAsite2-2a + H2      //0
	addReaction(rxnV3, Reaction(1.155e+05, 2.310e+00, 8.819e+00, sp::H2));          //7-HACAsite2-2a + H2 => 7-HACAsite2-1 + H      //1
	addReaction(rxnV3, Reaction(3.886e+03, 2.592e+00, 4.012e+00, sp::C2H2));                //7-HACAsite2-2a + Acety => 7-HACAsite2-3a      //2
	addReaction(rxnV3, Reaction(6.507e+12, 7.101e-01, 4.505e+01, sp::None));                //7-HACAsite2-3a => 7-HACAsite2-2a + Acety      //3
	addReaction(rxnV3, Reaction(5.755e+11, 7.034e-02, 2.983e+00, sp::None));                //7-HACAsite2-3a => 7-HACAsite2-4a      //4
	addReaction(rxnV3, Reaction(1.742e+12, 4.186e-01, 2.904e+01, sp::None));                //7-HACAsite2-4a => 7-HACAsite2-3a      //5
	addReaction(rxnV3, Reaction(3.207e+10, 9.581e-01, 2.313e+01, sp::None));                //7-HACAsite2-4a => 7-HACAsite2-5 + H   //6
	addReaction(rxnV3, Reaction(1.293e+08, 1.505e+00, 7.425e+00, sp::H));           //7-HACAsite2-5 + H => 7-HACAsite2-4a   //7
	//Pathway B
	addReaction(rxnV3, Reaction(1.479e+07, 1.854e+00, 1.707e+01, sp::H));           //7-HACAsite2-1 + H => 7-HACAsite2-2b + H2      //8
	addReaction(rxnV3, Reaction(5.914e+04, 2.234e+00, 1.187e+01, sp::H2));          //7-HACAsite2-2b + H2 => 7-HACAsite2-1 + H      //9
	addReaction(rxnV3, Reaction(1.098e+03, 2.581e+00, 7.651e+00, sp::C2H2));                //7-HACAsite2-2b + Acety => 7-HACAsite2-3b      //10
	addReaction(rxnV3, Reaction(2.894e+12, 7.099e-01, 3.830e+01, sp::None));                //7-HACAsite2-3b => 7-HACAsite2-2b + Acety      //11
	addReaction(rxnV3, Reaction(5.097e+11, 1.378e-01, 1.974e+01, sp::None));                //7-HACAsite2-3b => 7-HACAsite2-4b      //12
	addReaction(rxnV3, Reaction(9.936e+11, 4.098e-01, 4.526e+01, sp::None));                //7-HACAsite2-4b => 7-HACAsite2-3b      //13
	addReaction(rxnV3, Reaction(3.590e+11, 6.035e-01, 3.005e+01, sp::None));                //7-HACAsite2-4b => 7-HACAsite2-5 + H   //14
	addReaction(rxnV3, Reaction(6.258e+08, 1.380e+00, 2.451e+01, sp::H));           //7-HACAsite2-5 + H => 7-HACAsite2-4b   //15
	
	/*//Reactions provided by Menon2019 at B3LYP level. Might need to be changed.
	//Pathway A
	addReaction(rxnV3, Reaction(6.58594E+07,1.76634E+00,9.33825E+00, sp::H));			//H + 7-memb-1 <=> 7-memb-2a + H2		//0
	addReaction(rxnV3, Reaction(1.15481E+05,2.30995E+00,5.60390E+00, sp::H2));         	//H + 7-memb-1 <=> 7-memb-2a + H2		//1
	addReaction(rxnV3, Reaction(3.88643E+03,2.59182E+00,5.48872E+00, sp::C2H2));		//Acety + 7-memb-2a <=> 7-memb-3a		//2
	addReaction(rxnV3, Reaction(6.50688E+12,7.10131E-01,4.39670E+01, sp::None));       	//Acety + 7-memb-2a <=> 7-memb-3a		//3
	addReaction(rxnV3, Reaction(5.75482E+11,7.03403E-02,1.55667E+00, sp::None));		//7-memb-3a <=> 7-memb-4a				//4
	addReaction(rxnV3, Reaction(1.74248E+12,4.18589E-01,2.88623E+01, sp::None));       	//7-memb-3a <=> 7-memb-4a				//5
	addReaction(rxnV3, Reaction(3.20705E+10,9.58088E-01,2.25574E+01, sp::None));		//7-memb-4a <=> 7-memb-5 + H			//6
	addReaction(rxnV3, Reaction(1.29267E+08,1.50516E+00,3.23074E+00, sp::H));          	//7-memb-4a <=> 7-memb-5 + H			//7
	//Pathway B                 
	addReaction(rxnV3, Reaction(1.4788E+07, 1.85386E+00,9.82012E+00, sp::H));			//H + 7-memb-1 <=> 7-memb-2b + H2		//8
	addReaction(rxnV3, Reaction(5.91414E+04,2.23375E+00,8.96025E+00, sp::H2));         	//H + 7-memb-1 <=> 7-memb-2b + H2		//9
	addReaction(rxnV3, Reaction(1.09789E+03,2.58074E+00,9.42378E+00, sp::C2H2));       	//Acety + 7-memb-2b <=> 7-memb-3b		//10
	addReaction(rxnV3, Reaction(2.89436E+12,7.09941E-01,3.62228E+01, sp::None));       	//Acety + 7-memb-2b <=> 7-memb-3b		//11
	addReaction(rxnV3, Reaction(5.09705E+11,1.37833E-01,1.86806E+01, sp::None));       	//7-memb-3b <=> 7-memb-4b				//12
	addReaction(rxnV3, Reaction(9.93601E+11,4.09784E-01,4.63933E+01, sp::None));       	//7-memb-3b <=> 7-memb-4b				//13
	addReaction(rxnV3, Reaction(3.58968E+11,6.03489E-01,2.80765E+01, sp::None));       	//7-memb-4b <=> 7-memb-5 + H			//14
	addReaction(rxnV3, Reaction(6.25772E+08,1.37967E+00,1.71007E+01, sp::H));          	//7-memb-4b <=> 7-memb-5 + H			//15*/
	
	//OLD
	/*addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4              - Forward
	//Menon2019 DO NOT TRUST.
	addReaction(rxnV3, Reaction(9.720E+02, 3.10E+00, 6.10E+00, sp::C2H2));         // A3* + C2H2 -> A4 + H              	- 5             
	addReaction(rxnV3, Reaction(6.33E+03, 2.61E+00, 4.71E+00, sp::C2H2));         // A3* + H -> A4 + H              	- 6      */        
	//Using rates by On the low-temperature limit of HACA. Frenklach et al. 2018
	
	//addReaction(rxnV3, Reaction(1.190E+22, -2.450E+00, 1.889E+01, sp::C2H2));         // A3* + C2H2 -> A4 + H              	- 5              - Frenklach et al. 2018
	//addReaction(rxnV3, Reaction(1.060E+14, -4.900E-01, 8.204E+00, sp::C2H2));         // A3* + H -> A4 + H              	- 6              - Frenklach et al. 2018
	//addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 7              - Frenklach et al. 2018
	//addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 8              - Frenklach et al. 2018

	m_sType = R5R6AC; // sitetype
	m_name = "R7 growth on R5R6AC"; // name of process
	m_ID = 35;
}

double GR7_R5R6AC::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(6, 6);
	boost::numeric::ublas::vector<double> arr2(6);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	arr1(0,0) = m_r[1] + m_r[2];
    arr1(0,1) = - m_r[3];
    arr1(1,0) = -m_r[2]; 
    arr1(1,1) = m_r[3] + m_r[4];
    arr1(1,2) = - m_r[5];
    arr1(2,1) = - m_r[4];
    arr1(2,2) = m_r[5] + m_r[6];
    arr1(3,3) = m_r[9] + m_r[10];
    arr1(3,4) = - m_r[11];
    arr1(4,3) = - m_r[10];
    arr1(4,4) = + m_r[11] + m_r[12];
    arr1(4,5) = - m_r[13];
    arr1(5,4) = - m_r[12];
    arr1(5,5) = + m_r[13] + m_r[14];
    arr2(0) = +m_r[0];
    arr2(2) = 0.0*(+1.0*m_r[7]);
    arr2(3) = +m_r[8];
    arr2(5) = 0.0*(+1.0*m_r[15]);
	
	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	return m_rate = (m_r[6] * arr2(2) + m_r[14] * arr2(5))*site_count;
	//OLD
	/*double r_denom = m_r[1] + m_r[3] + m_r[4] + m_r[5] + m_r[6];
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = (m_r[5] + m_r[6]) * r_f* site_count; // Rate Equation*/
}

// ************************************************************
// ID36 - R7 growth on embedded-obstructed R5-2
// ************************************************************

// Elementary rate constants, site type, process type and name
void GR7_FEACR5::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//Reactions provided by Menon2019 at M06-2X level.
	//Pathway A
	addReaction(rxnV3, Reaction(5.897e+07, 1.847e+00, 1.712e+01, sp::H));           		//7-memb-1 + H => 7-memb-2a + H2        //0
	addReaction(rxnV3, Reaction(1.215e+05, 2.229e+00, 7.720e+00, sp::H2));          		//7-memb-2a + H2 => 7-memb-1 + H        //1
	addReaction(rxnV3, Reaction(1.348e+03, 2.573e+00, 4.935e+00, sp::C2H2));                //7-memb-2a + Acety => 7-memb-3a        //2
	addReaction(rxnV3, Reaction(2.366e+12, 7.052e-01, 3.967e+01, sp::None));                //7-memb-3a => 7-memb-2a + Acety        //3
	addReaction(rxnV3, Reaction(1.958e+11, 1.107e-01, 2.533e+01, sp::None));                //7-memb-3a => 7-memb-4a        //4
	addReaction(rxnV3, Reaction(3.412e+11, 6.253e-01, 5.337e+01, sp::None));                //7-memb-4a => 7-memb-3a        //5
	addReaction(rxnV3, Reaction(1.770e+10, 1.094e+00, 2.715e+01, sp::None));                //7-memb-4a => 7-memb-5 + H     //6
	addReaction(rxnV3, Reaction(5.321e+07, 1.515e+00, 7.095e+00, sp::H));           		//7-memb-5 + H => 7-memb-4a     //7
	//Pathway B
	addReaction(rxnV3, Reaction(5.315e+07, 1.858e+00, 1.612e+01, sp::H));           		//7-memb-1 + H => 7-memb-2b + H2        //8
	addReaction(rxnV3, Reaction(9.106e+04, 2.277e+00, 7.007e+00, sp::H2));          		//7-memb-2b + H2 => 7-memb-1 + H        //9
	addReaction(rxnV3, Reaction(3.521e+03, 2.598e+00, 3.998e+00, sp::C2H2));                //7-memb-2b + Acety => 7-memb-3b        //10
	addReaction(rxnV3, Reaction(4.736e+12, 7.020e-01, 4.080e+01, sp::None));                //7-memb-3b => 7-memb-2b + Acety        //11
	addReaction(rxnV3, Reaction(1.125e+11, 1.277e-01, 3.051e+01, sp::None));                //7-memb-3b => 7-memb-4b        //12
	addReaction(rxnV3, Reaction(2.383e+11, 5.964e-01, 5.790e+01, sp::None));                //7-memb-4b => 7-memb-3b        //13
	addReaction(rxnV3, Reaction(1.505e+10, 1.076e+00, 2.884e+01, sp::None));                //7-memb-4b => 7-memb-5 + H     //14
	addReaction(rxnV3, Reaction(5.841e+07, 1.533e+00, 7.084e+00, sp::H));           		//7-memb-5 + H => 7-memb-4b     //15
		
	/*//Reactions provided by Menon2019 at B3LYP level. Might need to be changed.
	//Pathway A
	addReaction(rxnV3, Reaction(5.8971E+07,1.8466E+00,9.3299E+00, sp::H));			//H + 7-memb-1 <=> 7-memb-2a + H2		//0
	addReaction(rxnV3, Reaction(1.2148E+05,2.2285E+00,5.0102E+00, sp::H2));         //H + 7-memb-1 <=> 7-memb-2a + H2		//1
	addReaction(rxnV3, Reaction(1.3478E+03,2.5734E+00,6.1505E+00, sp::C2H2));		//Acety + 7-memb-2a <=> 7-memb-3a		//2
	addReaction(rxnV3, Reaction(2.3662E+12,7.0522E-01,3.5866E+01, sp::None));       //Acety + 7-memb-2a <=> 7-memb-3a		//3
	addReaction(rxnV3, Reaction(1.9579E+11,1.1067E-01,2.3440E+01, sp::None));		//7-memb-3a <=> 7-memb-4a				//4
	addReaction(rxnV3, Reaction(3.4121E+11,6.2533E-01,5.1779E+01, sp::None));       //7-memb-3a <=> 7-memb-4a				//5
	addReaction(rxnV3, Reaction(1.7700E+10,1.0938E+00,2.6150E+01, sp::None));		//7-memb-4a <=> 7-memb-5 + H			//6
	addReaction(rxnV3, Reaction(5.3209E+07,1.5146E+00,2.8370E+00, sp::H));          //7-memb-4a <=> 7-memb-5 + H			//7
	//Pathway B
	addReaction(rxnV3, Reaction(5.3149E+07,1.8580E+00,8.3759E+00, sp::H));			//H + 7-memb-1 <=> 7-memb-2b + H2		//8
	addReaction(rxnV3, Reaction(9.1063E+04,2.2769E+00,4.2404E+00, sp::H2));         //H + 7-memb-1 <=> 7-memb-2b + H2		//9
	addReaction(rxnV3, Reaction(3.5208E+03,2.5977E+00,4.7459E+00, sp::C2H2));       //Acety + 7-memb-2b <=> 7-memb-3b		//10
	addReaction(rxnV3, Reaction(4.7358E+12,7.0201E-01,3.7437E+01, sp::None));       //Acety + 7-memb-2b <=> 7-memb-3b		//11
	addReaction(rxnV3, Reaction(1.1252E+11,1.2767E-01,2.8620E+01, sp::None));       //7-memb-3b <=> 7-memb-4b				//12
	addReaction(rxnV3, Reaction(2.3828E+11,5.9635E-01,5.5326E+01, sp::None));       //7-memb-3b <=> 7-memb-4b				//13
	addReaction(rxnV3, Reaction(1.5046E+10,1.0755E+00,2.7612E+01, sp::None));       //7-memb-4b <=> 7-memb-5 + H			//14
	addReaction(rxnV3, Reaction(5.8411E+07,1.5329E+00,2.7734E+00, sp::H));          //7-memb-4b <=> 7-memb-5 + H			//15*/
	
	/*Previous.
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4              - Forward
	//Menon2019 DO NOT TRUST.
	addReaction(rxnV3, Reaction(9.720E+02, 3.10E+00, 6.10E+00, sp::C2H2));         // A3* + C2H2 -> A4 + H              	- 5             
	addReaction(rxnV3, Reaction(6.33E+03, 2.61E+00, 4.71E+00, sp::C2H2));         // A3* + H -> A4 + H              	- 6    */
	//Using rates by On the low-temperature limit of HACA. Frenklach et al. 2018
	//addReaction(rxnV3, Reaction(1.190E+22, -2.450E+00, 1.889E+01, sp::C2H2));         // A3* + C2H2 -> A4 + H              	- 5              - Frenklach et al. 2018
	//addReaction(rxnV3, Reaction(1.060E+14, -4.900E-01, 8.204E+00, sp::C2H2));         // A3* + H -> A4 + H              	- 6              - Frenklach et al. 2018
	//addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 7              - Frenklach et al. 2018
	//addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 8              - Frenklach et al. 2018
	m_sType = FEACR5; // sitetype
	m_name = "R7 growth on FEACR5"; // name of process
	m_ID = 36;
}

double GR7_FEACR5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(6, 6);
	boost::numeric::ublas::vector<double> arr2(6);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	arr1(0,0) = m_r[1] + m_r[2];
    arr1(0,1) = - m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
    arr1(1,2) = - m_r[5];
    arr1(2,1) = - m_r[4];
    arr1(2,2) = m_r[5] + m_r[6];
    arr1(3,3) = m_r[9] + m_r[10];
    arr1(3,4) = - m_r[11];
    arr1(4,3) = - m_r[10];
    arr1(4,4) = + m_r[11] + m_r[12];
    arr1(4,5) = - m_r[13];
    arr1(5,4) = - m_r[12];
    arr1(5,5) = + m_r[13] + m_r[14];
    arr2(0) = +m_r[0];
    arr2(2) = 0.0*(+1.0*m_r[7]);
    arr2(3) = +m_r[8];
    arr2(5) = 0.0*(+1.0*m_r[15]);
	
	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	return m_rate = (m_r[6] * arr2(2) + m_r[14] * arr2(5))*site_count;
	
	/*OLD
	double r_denom = m_r[1] + m_r[3] + m_r[4] + m_r[5] + m_r[6];
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = (m_r[5] + m_r[6]) * r_f* site_count; // Rate Equation*/
}

// ************************************************************
// ID37 - R6 growth on R5R6ZZ
// ************************************************************
// Elementary rate constants, site type, process type and name
void G6R_R5R6ZZ::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//Checked 2019
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 6 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 7 Hou2017
	
	//Previous. No idea about the C2H2 rate
	/*addReaction(rxnV3, Reaction(4.20e13, 0, 13.00, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.90e12, 0, 11.00, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(8.00e07, 1.560, 3.800, sp::C2H2));  // 5 - r4f*/

	m_sType = R5R6ZZ; // sitetype
	m_name = "G6R at R5R6ZZ"; // name of process
	m_ID = 37;
}

double G6R_R5R6ZZ::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_denom2 = (m_r[7] + m_r[3] + m_r[4] + m_r[5]);
	double r_f, r_f2; // radical fraction
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	if (r_denom2>0) {
		r_f2 = (m_r[6] + m_r[2]) / r_denom2;
	}
	else r_f2 = 0;
	return m_rate = m_r[5] * (r_f + r_f2) * site_count; // Rate Equation
}

// ************************************************************
// ID38 - R7-bay closure at ACACR5
// ************************************************************
// Elementary rate constants, site type, process type and name
void L7_ACACR5::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//M06-2X level.
	//Pathway A
	addReaction(rxnV3, Reaction(3.915e+07, 1.876e+00, 9.421e+00, sp::H));           		//7-memb2-start + H => 7-memb2-habs-1a + H2     //0
	addReaction(rxnV3, Reaction(5.369e+04, 2.275e+00, 5.583e+00, sp::H2));          		//7-memb2-habs-1a + H2 => 7-memb2-start + H     //1
	addReaction(rxnV3, Reaction(8.513e+11, 1.358e-01, 4.510e+00, sp::None));                //7-memb2-habs-1a => 7-memb2-habs-2a    //2
	addReaction(rxnV3, Reaction(3.523e+12, 2.926e-01, 2.567e+01, sp::None));                //7-memb2-habs-2a => 7-memb2-habs-1a    //3
	addReaction(rxnV3, Reaction(2.033e+10, 1.067e+00, 3.160e+01, sp::None));                //7-memb2-habs-2a => 7-memb2-end + H    //4
	addReaction(rxnV3, Reaction(1.033e+08, 1.495e+00, 2.895e+00, sp::H));           		//7-memb2-end + H => 7-memb2-habs-2a    //5
	//Pathway B
	addReaction(rxnV3, Reaction(3.091e+07, 1.891e+00, 9.308e+00, sp::H));           		//7-memb2-start + H => 7-memb2-habs-1b + H2     //6
	addReaction(rxnV3, Reaction(5.144e+04, 2.267e+00, 7.132e+00, sp::H2));          		//7-memb2-habs-1b + H2 => 7-memb2-start + H     //7
	addReaction(rxnV3, Reaction(7.041e+11, 1.839e-01, 1.034e+01, sp::None));                //7-memb2-habs-1b => 7-memb2-habs-2b    //8
	addReaction(rxnV3, Reaction(2.944e+12, 4.130e-01, 2.862e+01, sp::None));                //7-memb2-habs-2b => 7-memb2-habs-1b    //9
	addReaction(rxnV3, Reaction(1.861e+10, 1.136e+00, 2.957e+01, sp::None));                //7-memb2-habs-2b => 7-memb2-end + H    //10
	addReaction(rxnV3, Reaction(7.712e+07, 1.514e+00, 2.067e+00, sp::H));           		//7-memb2-end + H => 7-memb2-habs-2b    //11
	
	//Carbene route
	//Pathway A
	addReaction(rxnV3, Reaction(8.031e+10, 8.996e-01, 9.583e+01, sp::None));                //7-memb2-start => 7-memb2-carb-1a      //12
	addReaction(rxnV3, Reaction(4.398e+11, 3.592e-01, 3.385e+00, sp::None));                //7-memb2-carb-1a => 7-memb2-start      //13
	addReaction(rxnV3, Reaction(8.031e+11, 1.032e-02, 8.456e+00, sp::None));                //7-memb2-carb-1a => 7-memb2-carb-2a    //14
	addReaction(rxnV3, Reaction(1.897e+12, 2.234e-01, 1.779e+01, sp::None));                //7-memb2-carb-2a => 7-memb2-carb-1a    //15
	addReaction(rxnV3, Reaction(5.759e+11, 3.930e-01, 8.715e-01, sp::None));                //7-memb2-carb-2a => 7-memb2-carb-4a    //16
	addReaction(rxnV3, Reaction(1.052e+11, 9.049e-01, 5.350e+01, sp::None));                //7-memb2-carb-4a => 7-memb2-carb-2a    //17
	addReaction(rxnV3, Reaction(8.873e+10, 6.388e-01, 3.131e+01, sp::None));                //7-memb2-carb-4a => 7-memb2-carb-5a    //18
	addReaction(rxnV3, Reaction(1.728e+10, 7.118e-01, 6.065e+01, sp::None));                //7-memb2-carb-5a => 7-memb2-carb-4a    //19
	addReaction(rxnV3, Reaction(3.907e+09, 1.273e+00, 9.705e+01, sp::None));                //7-memb2-carb-5a => 7-memb2-end + H2   //20
	addReaction(rxnV3, Reaction(2.448e+05, 1.999e+00, 8.640e+01, sp::H2));          		//7-memb2-end + H2 => 7-memb2-carb-5a   //21
	//Pathway B
	addReaction(rxnV3, Reaction(1.061e+11, 7.990e-01, 8.426e+01, sp::None));                //7-memb2-start => 7-memb2-carb-1b      //22
	addReaction(rxnV3, Reaction(5.486e+11, 3.348e-01, 2.012e+00, sp::None));                //7-memb2-carb-1b => 7-memb2-start      //23
	addReaction(rxnV3, Reaction(1.000e+12, -1.397e-02, 3.568e+00, sp::None));               //7-memb2-carb-1b => 7-memb2-carb-2b    //24
	addReaction(rxnV3, Reaction(2.167e+12, 5.562e-01, 5.928e+01, sp::None));                //7-memb2-carb-2b => 7-memb2-carb-1b    //25
	addReaction(rxnV3, Reaction(3.063e+11, 8.244e-01, 6.356e+01, sp::None));                //7-memb2-carb-2b => 7-memb2-carb-4b    //26
	addReaction(rxnV3, Reaction(1.514e+10, 6.739e-01, 5.031e+01, sp::None));                //7-memb2-carb-4b => 7-memb2-carb-2b    //27
	addReaction(rxnV3, Reaction(4.479e+09, 7.144e-01, 2.710e+01, sp::None));                //7-memb2-carb-4b => 7-memb2-end + H2   //28
	addReaction(rxnV3, Reaction(2.332e+05, 1.742e+00, 5.510e+01, sp::H2));          		//7-memb2-end + H2 => 7-memb2-carb-4b   //29
	
	/*//B3LYP level.
	//Pathway A
	addReaction(rxnV3, Reaction(2.76732E+07,1.91279E+00,9.54190E+00, sp::H));		//H + 7-membHabs-1 <=> 7-membHabs-2a + H2		//0
	addReaction(rxnV3, Reaction(4.21185E+04,2.26366E+00,6.87815E+00, sp::H2));      //H + 7-membHabs-1 <=> 7-membHabs-2a + H2		//1
	addReaction(rxnV3, Reaction(4.70295E+11,1.42940E-01,4.72228E+00, sp::None));	//7-membHabs-2a <=> 7-membHabs-3a				//2
	addReaction(rxnV3, Reaction(1.47581E+12,3.66538E-01,2.73731E+01, sp::None));    //7-membHabs-2a <=> 7-membHabs-3a				//3
	addReaction(rxnV3, Reaction(6.42375E+09,1.09295E+00,3.21573E+01, sp::None));	//7-membHabs-3a <=> 7-membHabs-4 + H			//4
	addReaction(rxnV3, Reaction(1.69261E+08,1.52168E+00,1.63656E+00, sp::H));       //7-membHabs-3a <=> 7-membHabs-4 + H			//5
	//Pathway B
	addReaction(rxnV3, Reaction(2.84300E+07,1.90572E+00,9.53262E+00, sp::H));		//H + 7-membHabs-1 <=> 7-membHabs-2b + H2		//6
	addReaction(rxnV3, Reaction(5.33824E+04,2.26128E+00,7.52478E+00, sp::H2));      //H + 7-membHabs-1 <=> 7-membHabs-2b + H2		//7
	addReaction(rxnV3, Reaction(6.59893E+11,8.15588E-02,2.62478E+00, sp::None));    //7-membHabs-2b <=> 7-membHabs-3b				//8
	addReaction(rxnV3, Reaction(5.13229E+12,3.39661E-01,2.54015E+01, sp::None));    //7-membHabs-2b <=> 7-membHabs-3b				//9
	addReaction(rxnV3, Reaction(2.00638E+10,1.09912E+00,3.28066E+01, sp::None));    //7-membHabs-3b <=> 7-membHabs-4 + H			//10
	addReaction(rxnV3, Reaction(1.72899E+08,1.48867E+00,1.50412E+00, sp::H));       //7-membHabs-3b <=> 7-membHabs-4 + H			//11*/
	
    m_sType = ACACR5; // sitetype
    m_name = "BY7 closure on ACACR5"; // name of process
    m_ID = 38;
}
// Jump rate calculation
double L7_ACACR5::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double L7_ACACR5::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double L7_ACACR5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
	double T = gp[gp.T];
    if (T <= 1000) return m_rate = 0.0;
	//Rate assuming PEQ approximation
	matrix<double> arr1(4, 4);
	boost::numeric::ublas::vector<double> arr2(4);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	arr1(0,0) = m_r[1] + m_r[2];
	arr1(0,1) = - m_r[3];
	arr1(1,0) = -m_r[2];
	arr1(1,1) = + m_r[3] + m_r[4];
	arr1(2,2) = +m_r[7] + m_r[9];
	arr1(2,3) = - m_r[8];
	arr1(3,2) = - m_r[9];
	arr1(3,3) = m_r[8] + m_r[10];
	arr2(0) = +m_r[0];
	arr2(1) = 0.0*(m_r[5]);
	arr2(2) = +m_r[6];
	arr2(3) = 0.0*(m_r[11]);
	
	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	double rate = (m_r[4] * arr2(1) + m_r[10] * arr2(3))*site_count;
	
	if (T>=2200) {
		matrix<double> arr1_carbene(7, 7);
		boost::numeric::ublas::vector<double> arr2_carbene(7);
		for (unsigned k = 0; k < arr1_carbene.size1(); ++k)
			for (unsigned l = 0; l < arr1_carbene.size2(); ++l)
				arr1_carbene(k, l) = 0.0;
		for (unsigned k = 0; k < arr2_carbene.size(); ++k)
			arr2_carbene(k) = 0.0;
		
		arr1_carbene(0,0) = m_r[13] + m_r[14];
		arr1_carbene(0,1) = -m_r[15];
		arr1_carbene(1,0) = -m_r[14];
		arr1_carbene(1,1) = m_r[15] + m_r[16];
		arr1_carbene(1,2) = -m_r[17];
		arr1_carbene(2,1) = -m_r[16];
		arr1_carbene(2,2) = m_r[17] + m_r[18];
		arr1_carbene(2,3) = -m_r[19];
		arr1_carbene(3,2) = -m_r[18];
		arr1_carbene(3,3) = m_r[19] + m_r[20];
		arr1_carbene(4,4) = m_r[23] + m_r[24];
		arr1_carbene(4,5) = -m_r[25];
		arr1_carbene(5,4) = -m_r[24];
		arr1_carbene(5,5) = m_r[25] + m_r[26];
		arr1_carbene(5,6) = -m_r[27];
		arr1_carbene(6,5) = -m_r[26];
		arr1_carbene(6,6) = m_r[27] + m_r[28];
		arr2_carbene(0) = +m_r[12];
		arr2_carbene(3) = 0.0*(m_r[21]);
		arr2_carbene(4) = +m_r[22];
		arr2_carbene(6) = 0.0*(m_r[29]);
		
		permutation_matrix<size_t> pm_carbene(arr1_carbene.size1());
		lu_factorize(arr1_carbene, pm_carbene);
		lu_substitute(arr1_carbene, pm_carbene, arr2_carbene);
		
		rate += (m_r[20] * arr2_carbene(3) + m_r[28] * arr2_carbene(6)) * site_count;
	}
		
	return m_rate = rate;
}

// ************************************************************
// ID39 - R6 growth on R5R6FER
// ************************************************************
// Elementary rate constants, site type, process type and name
void G6R_R5R6FER::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//Checked 2019
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	
	//Previous. No idea about the C2H2 rate
	/*addReaction(rxnV3, Reaction(4.20e13, 0, 13.00, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.90e12, 0, 11.00, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(8.00e07, 1.560, 3.800, sp::C2H2));  // 5 - r4f*/

	m_sType = R5R6FER; // sitetype
	m_name = "G6R at R5R6FER"; // name of process
	m_ID = 39;
}

double G6R_R5R6FER::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation
}

// ************************************************************
// ID40 - R6 growth on R5R6FER5R6
// ************************************************************
// Elementary rate constants, site type, process type and name
void G6R_R5R6FER5R6::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//Checked 2019
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	
	//Previous. No idea about the C2H2 rate
	/*addReaction(rxnV3, Reaction(4.20e13, 0, 13.00, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(3.90e12, 0, 11.00, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(8.00e07, 1.560, 3.800, sp::C2H2));  // 5 - r4f*/

	m_sType = R5R6FER5R6; // sitetype
	m_name = "G6R at R5R6FER5R6"; // name of process
	m_ID = 40;
}

double G6R_R5R6FER5R6::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation
}

// ************************************************************
// ID41 - R7-bay closure at FEZZACR5
// ************************************************************
// Elementary rate constants, site type, process type and name
void L7_FEZZACR5::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//M06-2X level.
	//Pathway A
	addReaction(rxnV3, Reaction(2.767e+07, 1.913e+00, 9.542e+00, sp::H));           		//7-membHabs-1 + H => 7-membHabs-2a + H2        //0
	addReaction(rxnV3, Reaction(4.212e+04, 2.264e+00, 6.878e+00, sp::H2));          		//7-membHabs-2a + H2 => 7-membHabs-1 + H        //1
	addReaction(rxnV3, Reaction(4.703e+11, 1.429e-01, 4.722e+00, sp::None));                //7-membHabs-2a => 7-membHabs-3a        //2
	addReaction(rxnV3, Reaction(1.476e+12, 3.665e-01, 2.737e+01, sp::None));                //7-membHabs-3a => 7-membHabs-2a        //3
	addReaction(rxnV3, Reaction(6.424e+09, 1.093e+00, 3.216e+01, sp::None));                //7-membHabs-3a => 7-membHabs-4 + H     //4
	addReaction(rxnV3, Reaction(1.693e+08, 1.522e+00, 1.637e+00, sp::H));           		//7-membHabs-4 + H => 7-membHabs-3a     //5
	//Pathway B
	addReaction(rxnV3, Reaction(2.843e+07, 1.906e+00, 9.533e+00, sp::H));           		//7-membHabs-1 + H => 7-membHabs-2b + H2        //6
	addReaction(rxnV3, Reaction(5.338e+04, 2.261e+00, 7.525e+00, sp::H2));          		//7-membHabs-2b + H2 => 7-membHabs-1 + H        //7
	addReaction(rxnV3, Reaction(6.599e+11, 8.156e-02, 2.625e+00, sp::None));                //7-membHabs-2b => 7-membHabs-3b        //8
	addReaction(rxnV3, Reaction(5.132e+12, 3.397e-01, 2.540e+01, sp::None));                //7-membHabs-3b => 7-membHabs-2b        //9
	addReaction(rxnV3, Reaction(2.006e+10, 1.099e+00, 3.281e+01, sp::None));                //7-membHabs-3b => 7-membHabs-4 + H     //10
	addReaction(rxnV3, Reaction(1.729e+08, 1.489e+00, 1.504e+00, sp::H));           		//7-membHabs-4 + H => 7-membHabs-3b     //11
	//Pathway A - Carbene
	addReaction(rxnV3, Reaction(1.603e+11, 7.771e-01, 8.323e+01, sp::None));                //7-membHabs-1 => 7-memb-carb-1a        //12
	addReaction(rxnV3, Reaction(3.050e+11, 2.944e-01, 8.342e-02, sp::None));                //7-memb-carb-1a => 7-membHabs-1        //13
	addReaction(rxnV3, Reaction(2.488e+11, 1.200e-01, 1.162e+01, sp::None));                //7-memb-carb-1a => 7-memb-carb-2a      //14
	addReaction(rxnV3, Reaction(1.472e+12, 6.757e-01, 4.509e+01, sp::None));                //7-memb-carb-2a => 7-memb-carb-1a      //15
	addReaction(rxnV3, Reaction(1.397e+11, 5.806e-01, 2.701e+01, sp::None));                //7-memb-carb-2a => 7-memb-carb-4a      //16
	addReaction(rxnV3, Reaction(3.264e+10, 7.336e-01, 5.681e+01, sp::None));                //7-memb-carb-4a => 7-memb-carb-2a      //17
	addReaction(rxnV3, Reaction(5.515e+10, 8.485e-01, 5.932e+01, sp::None));                //7-memb-carb-4a => 7-memb-carb-5a      //18
	addReaction(rxnV3, Reaction(7.566e+10, 6.747e-01, 3.856e+01, sp::None));                //7-memb-carb-5a => 7-memb-carb-4a      //19
	addReaction(rxnV3, Reaction(6.206e+09, 8.471e-01, 2.861e+01, sp::None));                //7-memb-carb-5a => 7-membHabs-4 + H2   //20
	addReaction(rxnV3, Reaction(2.164e+05, 1.798e+00, 5.870e+01, sp::H2));          		//7-membHabs-4 + H2 => 7-memb-carb-5a   //21
	//Pathway B - Carbene
	addReaction(rxnV3, Reaction(1.162e+11, 8.365e-01, 8.601e+01, sp::None));                //7-membHabs-1 => 7-memb-carb-1b        //22
	addReaction(rxnV3, Reaction(4.948e+11, 3.313e-01, 1.476e+00, sp::None));                //7-memb-carb-1b => 7-membHabs-1        //23
	addReaction(rxnV3, Reaction(5.744e+11, 3.921e-02, 8.721e+00, sp::None));                //7-memb-carb-1b => 7-memb-carb-2b      //24
	addReaction(rxnV3, Reaction(2.164e+12, 2.924e-01, 1.861e+01, sp::None));                //7-memb-carb-2b => 7-memb-carb-1b      //25
	addReaction(rxnV3, Reaction(9.900e+11, 3.314e-01, -2.850e+00, sp::None));               //7-memb-carb-2b => 7-memb-carb-4b      //26
	addReaction(rxnV3, Reaction(1.025e+11, 8.748e-01, 5.711e+01, sp::None));                //7-memb-carb-4b => 7-memb-carb-2b      //27
	addReaction(rxnV3, Reaction(6.156e+10, 7.822e-01, 3.925e+01, sp::None));                //7-memb-carb-4b => 7-memb-carb-5b      //28
	addReaction(rxnV3, Reaction(9.567e+10, 6.963e-01, 3.520e+01, sp::None));                //7-memb-carb-5b => 7-memb-carb-4b      //29
	addReaction(rxnV3, Reaction(3.478e+09, 1.288e+00, 8.743e+01, sp::None));                //7-memb-carb-5b => 7-membHabs-4 + H2   //30
	addReaction(rxnV3, Reaction(1.696e+05, 2.086e+00, 9.563e+01, sp::H2));          		//7-membHabs-4 + H2 => 7-memb-carb-5b   //31
	
	//B3LYP level. Old.
	/*//Pathway A
	addReaction(rxnV3, Reaction(2.76732E+07,1.91279E+00,9.54190E+00, sp::H));		//H + 7-membHabs-1 <=> 7-membHabs-2a + H2		//0
	addReaction(rxnV3, Reaction(4.21185E+04,2.26366E+00,6.87815E+00, sp::H2));      //H + 7-membHabs-1 <=> 7-membHabs-2a + H2		//1
	addReaction(rxnV3, Reaction(4.70295E+11,1.42940E-01,4.72228E+00, sp::None));	//7-membHabs-2a <=> 7-membHabs-3a				//2
	addReaction(rxnV3, Reaction(1.47581E+12,3.66538E-01,2.73731E+01, sp::None));    //7-membHabs-2a <=> 7-membHabs-3a				//3
	addReaction(rxnV3, Reaction(6.42375E+09,1.09295E+00,3.21573E+01, sp::None));	//7-membHabs-3a <=> 7-membHabs-4 + H			//4
	addReaction(rxnV3, Reaction(1.69261E+08,1.52168E+00,1.63656E+00, sp::H));       //7-membHabs-3a <=> 7-membHabs-4 + H			//5
	//Pathway B
	addReaction(rxnV3, Reaction(2.84300E+07,1.90572E+00,9.53262E+00, sp::H));		//H + 7-membHabs-1 <=> 7-membHabs-2b + H2		//6
	addReaction(rxnV3, Reaction(5.33824E+04,2.26128E+00,7.52478E+00, sp::H2));      //H + 7-membHabs-1 <=> 7-membHabs-2b + H2		//7
	addReaction(rxnV3, Reaction(6.59893E+11,8.15588E-02,2.62478E+00, sp::None));    //7-membHabs-2b <=> 7-membHabs-3b				//8
	addReaction(rxnV3, Reaction(5.13229E+12,3.39661E-01,2.54015E+01, sp::None));    //7-membHabs-2b <=> 7-membHabs-3b				//9
	addReaction(rxnV3, Reaction(2.00638E+10,1.09912E+00,3.28066E+01, sp::None));    //7-membHabs-3b <=> 7-membHabs-4 + H			//10
	addReaction(rxnV3, Reaction(1.72899E+08,1.48867E+00,1.50412E+00, sp::H));       //7-membHabs-3b <=> 7-membHabs-4 + H			//11*/
	
    m_sType = FEZZACR5; // sitetype
    m_name = "BY7 closure on FEZZACR5"; // name of process
    m_ID = 41;
}
// Jump rate calculation
double L7_FEZZACR5::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double L7_FEZZACR5::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double L7_FEZZACR5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
	double T = gp[gp.T];
    if (T <= 800) return m_rate = 0.0;
	
	//Rate assuming PEQ approximation
	matrix<double> arr1(4, 4);
	boost::numeric::ublas::vector<double> arr2(4);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	arr1(0,0) = m_r[1] + m_r[2];
    arr1(0,1) = - m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = + m_r[3] + m_r[4];
    arr1(2,2) = +m_r[7] + m_r[9];
    arr1(2,3) = - m_r[8];
    arr1(3,2) = - m_r[9];
    arr1(3,3) = m_r[8] + m_r[10];
    arr2(0) = +m_r[0];
    arr2(1) = 0.0*(m_r[5]);
    arr2(2) = +m_r[6];
    arr2(3) = 0.0*(m_r[11]);
	
	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	double rate = (m_r[4] * arr2(1) + m_r[10] * arr2(3))*site_count;
	
	if (T>=2000) {
		matrix<double> arr1_carbene(8, 8);
		boost::numeric::ublas::vector<double> arr2_carbene(8);
		for (unsigned k = 0; k < arr1_carbene.size1(); ++k)
			for (unsigned l = 0; l < arr1_carbene.size2(); ++l)
				arr1_carbene(k, l) = 0.0;
		for (unsigned k = 0; k < arr2_carbene.size(); ++k)
			arr2_carbene(k) = 0.0;
		
		arr1_carbene(0,0) = m_r[13] + m_r[14];
		arr1_carbene(0,1) = -m_r[15];
		arr1_carbene(1,0) = -m_r[14];
		arr1_carbene(1,1) = m_r[15] + m_r[16];
		arr1_carbene(1,2) = -m_r[17];
		arr1_carbene(2,1) = -m_r[16];
		arr1_carbene(2,2) = m_r[17] + m_r[18];
		arr1_carbene(2,3) = -m_r[19];
		arr1_carbene(3,2) = -m_r[18];
		arr1_carbene(3,3) = m_r[19] + m_r[20];
		arr1_carbene(4,4) = m_r[23] + m_r[24];
		arr1_carbene(4,5) = -m_r[25];
		arr1_carbene(5,4) = -m_r[24];
		arr1_carbene(5,5) = m_r[25] + m_r[26];
		arr1_carbene(5,6) = -m_r[27];
		arr1_carbene(6,5) = -m_r[26];
		arr1_carbene(6,6) = m_r[27] + m_r[28];
		arr1_carbene(6,7) = -m_r[29];
		arr1_carbene(7,6) = -m_r[28];
		arr1_carbene(7,7) = m_r[29] + m_r[30];
		arr2_carbene(0) = +m_r[12];
		arr2_carbene(3) = 0.0*(m_r[21]);
		arr2_carbene(4) = +m_r[22];
		arr2_carbene(7) = 0.0*(m_r[31]);
		
		permutation_matrix<size_t> pm_carbene(arr1_carbene.size1());
		lu_factorize(arr1_carbene, pm_carbene);
		lu_substitute(arr1_carbene, pm_carbene, arr2_carbene);
		
		rate += (m_r[20] * arr2_carbene(3) + m_r[30] * arr2_carbene(7)) * site_count;
	}
		
	return m_rate = rate;
}

// ************************************************************
// ID42 - R5 conversion to R6 on RZZR (2R5 collision)
// ************************************************************
// Elementary rate constants, site type, process type and name
void C5R_RZZR::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    //H addition & migration
	addReaction(rxnV3, Reaction(5.400e+11, 4.500e-01, 1.820e+00, sp::H));          //0  //R12 in Whitesides2010
	addReaction(rxnV3, Reaction(1.791e+00, 0.000e+00, -3.5187e+01, sp::None));      //1  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(2.000e+12, 0.000e+00, 0.000e+00, sp::H));          //2  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(3.100e+11, 8.700e-01, 7.4323e+01, sp::None));      //3  //R7 in Whitesides2010
	addReaction(rxnV3, Reaction(6.700e+11, 8.400e-01, 7.0790e+01, sp::None));      //4  //R8 in Whitesides2010
	addReaction(rxnV3, Reaction(1.300e+11, 1.600e-01, 4.5900e+01, sp::None));      //5  //R9 in Whitesides2010
	addReaction(rxnV3, Reaction(8.900e+05, 2.280e+00, 6.1489e+01, sp::None));      //6  //R14 in Whitesides2010
    m_sType = RZZR; // sitetype
    m_name = "R5 conversion to R6 on RZZR (2R5 collision)"; // name of process
    m_ID = 42;
}
// Jump rate calculation
double C5R_RZZR::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    
    double r_denom = ( m_r[0]/m_r[1] + m_r[2]+m_r[3]+m_r[4] + m_r[6]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = 2*m_r[6]*r_f* site_count; // Rate Equation
}

// ************************************************************
// ID43 - R5 conversion to R6 on R5R6ZZR (2R5 collision)
// ************************************************************
// Elementary rate constants, site type, process type and name
void C5R_R5R6ZZR::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
    //H addition & migration
	addReaction(rxnV3, Reaction(5.400e+11, 4.500e-01, 1.820e+00, sp::H));          //0  //R12 in Whitesides2010
	addReaction(rxnV3, Reaction(1.791e+00, 0.000e+00, -3.5187e+01, sp::None));      //1  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(2.000e+12, 0.000e+00, 0.000e+00, sp::H));          //2  //R13 in Whitesides2010
	addReaction(rxnV3, Reaction(3.100e+11, 8.700e-01, 7.4323e+01, sp::None));      //3  //R7 in Whitesides2010
	addReaction(rxnV3, Reaction(6.700e+11, 8.400e-01, 7.0790e+01, sp::None));      //4  //R8 in Whitesides2010
	addReaction(rxnV3, Reaction(1.300e+11, 1.600e-01, 4.5900e+01, sp::None));      //5  //R9 in Whitesides2010
	addReaction(rxnV3, Reaction(8.900e+05, 2.280e+00, 6.1489e+01, sp::None));      //6  //R14 in Whitesides2010
    m_sType = R5R6ZZR; // sitetype
    m_name = "R5 conversion to R6 on R5R6ZZR (2R5 collision)"; // name of process
    m_ID = 43;
}
// Jump rate calculation
double C5R_R5R6ZZR::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    
    
    double r_denom = ( m_r[0]/m_r[1] + m_r[2]+m_r[3]+m_r[4] + m_r[6]);
    double r_f; // radical fraction 
    if(r_denom>0) {
        r_f = (m_r[0])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = 1*m_r[6]*r_f* site_count; // Rate Equation
}

// ************************************************************
// ID44- R5R6BY5 closure reaction
// ************************************************************
// Elementary rate constants, site type, process type and name
void L6_R5R6BY5::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
	/*addReaction(rxnV3, Reaction(9.24e07, 1.500, 9.646, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(9.60e04, 1.960, 9.021, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.11e11, 0.658, 23.99, sp::None));  // 5 - r4f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

	m_sType = R5R6BY5; // sitetype
	m_name = "R5R6BY5 closure"; // name of process
	m_ID = 44;
}
double L6_R5R6BY5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	/*// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}

// ************************************************************
// ID45- R5R6ACR closure reaction
// ************************************************************
// Elementary rate constants, site type, process type and name
void L6_R5R6ACR::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
	/*addReaction(rxnV3, Reaction(9.24e07, 1.500, 9.646, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(9.60e04, 1.960, 9.021, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.11e11, 0.658, 23.99, sp::None));  // 5 - r4f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

	m_sType = R5R6ACR; // sitetype
	m_name = "R5R6ACR closure"; // name of process
	m_ID = 45;
}
double L6_R5R6ACR::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	/*// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}

// ************************************************************
// ID46- R5R6ACR5R6 closure reaction
// ************************************************************
// Elementary rate constants, site type, process type and name
void L6_R5R6ACR5R6::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
	/*addReaction(rxnV3, Reaction(9.24e07, 1.500, 9.646, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(9.60e04, 1.960, 9.021, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.11e11, 0.658, 23.99, sp::None));  // 5 - r4f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

	m_sType = R5R6ACR5R6; // sitetype
	m_name = "R5R6ACR5R6 closure"; // name of process
	m_ID = 46;
}
double L6_R5R6ACR5R6::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	/*// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}

// ************************************************************
// ID47- ZZACR5 closure reaction
// ************************************************************
// Elementary rate constants, site type, process type and name
void L6_ZZACR5::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
	/*addReaction(rxnV3, Reaction(9.24e07, 1.500, 9.646, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(9.60e04, 1.960, 9.021, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.11e11, 0.658, 23.99, sp::None));  // 5 - r4f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

	m_sType = ZZACR5; // sitetype
	m_name = "ZZACR5 closure"; // name of process
	m_ID = 47;
}
double L6_ZZACR5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	/*// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}

// ************************************************************
// ID48- R5FEACR5 closure reaction
// ************************************************************
// Elementary rate constants, site type, process type and name
void L6_R5FEACR5::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
	/*addReaction(rxnV3, Reaction(9.24e07, 1.500, 9.646, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(9.60e04, 1.960, 9.021, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.11e11, 0.658, 23.99, sp::None));  // 5 - r4f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

	m_sType = R5FEACR5; // sitetype
	m_name = "R5FEACR5 closure"; // name of process
	m_ID = 48;
}
double L6_R5FEACR5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	/*// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}

// ************************************************************
// ID49- FEACR5FE closure reaction
// ************************************************************
// Elementary rate constants, site type, process type and name
void L6_FEACR5FE::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
	/*addReaction(rxnV3, Reaction(9.24e07, 1.500, 9.646, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(9.60e04, 1.960, 9.021, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.11e11, 0.658, 23.99, sp::None));  // 5 - r4f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

	m_sType = FEACR5FE; // sitetype
	m_name = "FEACR5FE closure"; // name of process
	m_ID = 49;
}
double L6_FEACR5FE::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	/*// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}

// ************************************************************
// ID50- R5ACR5R5 closure reaction
// ************************************************************
// Elementary rate constants, site type, process type and name
void L6_R5ACR5R5::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
	/*addReaction(rxnV3, Reaction(9.24e07, 1.500, 9.646, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(9.60e04, 1.960, 9.021, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.11e11, 0.658, 23.99, sp::None));  // 5 - r4f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

	m_sType = R5ACR5R5; // sitetype
	m_name = "R5ACR5R5 closure"; // name of process
	m_ID = 50;
}
double L6_R5ACR5R5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	/*// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}

// ************************************************************
// ID51 - R7-bay closure at R5ZZACR5
// ************************************************************
// Elementary rate constants, site type, process type and name
void L7_R5ZZACR5::initialise() {
    // Adding elementary reactions. The rate for this process was assumed to be similar to that of L&_ACACR5
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//M06-2X level.
	//Pathway A
	addReaction(rxnV3, Reaction(3.915e+07, 1.876e+00, 9.421e+00, sp::H));           		//7-memb2-start + H => 7-memb2-habs-1a + H2     //0
	addReaction(rxnV3, Reaction(5.369e+04, 2.275e+00, 5.583e+00, sp::H2));          		//7-memb2-habs-1a + H2 => 7-memb2-start + H     //1
	addReaction(rxnV3, Reaction(8.513e+11, 1.358e-01, 4.510e+00, sp::None));                //7-memb2-habs-1a => 7-memb2-habs-2a    //2
	addReaction(rxnV3, Reaction(3.523e+12, 2.926e-01, 2.567e+01, sp::None));                //7-memb2-habs-2a => 7-memb2-habs-1a    //3
	addReaction(rxnV3, Reaction(2.033e+10, 1.067e+00, 3.160e+01, sp::None));                //7-memb2-habs-2a => 7-memb2-end + H    //4
	addReaction(rxnV3, Reaction(1.033e+08, 1.495e+00, 2.895e+00, sp::H));           		//7-memb2-end + H => 7-memb2-habs-2a    //5
	//Pathway B
	addReaction(rxnV3, Reaction(3.091e+07, 1.891e+00, 9.308e+00, sp::H));           		//7-memb2-start + H => 7-memb2-habs-1b + H2     //6
	addReaction(rxnV3, Reaction(5.144e+04, 2.267e+00, 7.132e+00, sp::H2));          		//7-memb2-habs-1b + H2 => 7-memb2-start + H     //7
	addReaction(rxnV3, Reaction(7.041e+11, 1.839e-01, 1.034e+01, sp::None));                //7-memb2-habs-1b => 7-memb2-habs-2b    //8
	addReaction(rxnV3, Reaction(2.944e+12, 4.130e-01, 2.862e+01, sp::None));                //7-memb2-habs-2b => 7-memb2-habs-1b    //9
	addReaction(rxnV3, Reaction(1.861e+10, 1.136e+00, 2.957e+01, sp::None));                //7-memb2-habs-2b => 7-memb2-end + H    //10
	addReaction(rxnV3, Reaction(7.712e+07, 1.514e+00, 2.067e+00, sp::H));           		//7-memb2-end + H => 7-memb2-habs-2b    //11
	
	//Carbene route
	//Pathway A
	addReaction(rxnV3, Reaction(8.031e+10, 8.996e-01, 9.583e+01, sp::None));                //7-memb2-start => 7-memb2-carb-1a      //12
	addReaction(rxnV3, Reaction(4.398e+11, 3.592e-01, 3.385e+00, sp::None));                //7-memb2-carb-1a => 7-memb2-start      //13
	addReaction(rxnV3, Reaction(8.031e+11, 1.032e-02, 8.456e+00, sp::None));                //7-memb2-carb-1a => 7-memb2-carb-2a    //14
	addReaction(rxnV3, Reaction(1.897e+12, 2.234e-01, 1.779e+01, sp::None));                //7-memb2-carb-2a => 7-memb2-carb-1a    //15
	addReaction(rxnV3, Reaction(5.759e+11, 3.930e-01, 8.715e-01, sp::None));                //7-memb2-carb-2a => 7-memb2-carb-4a    //16
	addReaction(rxnV3, Reaction(1.052e+11, 9.049e-01, 5.350e+01, sp::None));                //7-memb2-carb-4a => 7-memb2-carb-2a    //17
	addReaction(rxnV3, Reaction(8.873e+10, 6.388e-01, 3.131e+01, sp::None));                //7-memb2-carb-4a => 7-memb2-carb-5a    //18
	addReaction(rxnV3, Reaction(1.728e+10, 7.118e-01, 6.065e+01, sp::None));                //7-memb2-carb-5a => 7-memb2-carb-4a    //19
	addReaction(rxnV3, Reaction(3.907e+09, 1.273e+00, 9.705e+01, sp::None));                //7-memb2-carb-5a => 7-memb2-end + H2   //20
	addReaction(rxnV3, Reaction(2.448e+05, 1.999e+00, 8.640e+01, sp::H2));          		//7-memb2-end + H2 => 7-memb2-carb-5a   //21
	//Pathway B
	addReaction(rxnV3, Reaction(1.061e+11, 7.990e-01, 8.426e+01, sp::None));                //7-memb2-start => 7-memb2-carb-1b      //22
	addReaction(rxnV3, Reaction(5.486e+11, 3.348e-01, 2.012e+00, sp::None));                //7-memb2-carb-1b => 7-memb2-start      //23
	addReaction(rxnV3, Reaction(1.000e+12, -1.397e-02, 3.568e+00, sp::None));               //7-memb2-carb-1b => 7-memb2-carb-2b    //24
	addReaction(rxnV3, Reaction(2.167e+12, 5.562e-01, 5.928e+01, sp::None));                //7-memb2-carb-2b => 7-memb2-carb-1b    //25
	addReaction(rxnV3, Reaction(3.063e+11, 8.244e-01, 6.356e+01, sp::None));                //7-memb2-carb-2b => 7-memb2-carb-4b    //26
	addReaction(rxnV3, Reaction(1.514e+10, 6.739e-01, 5.031e+01, sp::None));                //7-memb2-carb-4b => 7-memb2-carb-2b    //27
	addReaction(rxnV3, Reaction(4.479e+09, 7.144e-01, 2.710e+01, sp::None));                //7-memb2-carb-4b => 7-memb2-end + H2   //28
	addReaction(rxnV3, Reaction(2.332e+05, 1.742e+00, 5.510e+01, sp::H2));          		//7-memb2-end + H2 => 7-memb2-carb-4b   //29
	
	/*//B3LYP level.
	//Pathway A
	addReaction(rxnV3, Reaction(2.76732E+07,1.91279E+00,9.54190E+00, sp::H));		//H + 7-membHabs-1 <=> 7-membHabs-2a + H2		//0
	addReaction(rxnV3, Reaction(4.21185E+04,2.26366E+00,6.87815E+00, sp::H2));      //H + 7-membHabs-1 <=> 7-membHabs-2a + H2		//1
	addReaction(rxnV3, Reaction(4.70295E+11,1.42940E-01,4.72228E+00, sp::None));	//7-membHabs-2a <=> 7-membHabs-3a				//2
	addReaction(rxnV3, Reaction(1.47581E+12,3.66538E-01,2.73731E+01, sp::None));    //7-membHabs-2a <=> 7-membHabs-3a				//3
	addReaction(rxnV3, Reaction(6.42375E+09,1.09295E+00,3.21573E+01, sp::None));	//7-membHabs-3a <=> 7-membHabs-4 + H			//4
	addReaction(rxnV3, Reaction(1.69261E+08,1.52168E+00,1.63656E+00, sp::H));       //7-membHabs-3a <=> 7-membHabs-4 + H			//5
	//Pathway B
	addReaction(rxnV3, Reaction(2.84300E+07,1.90572E+00,9.53262E+00, sp::H));		//H + 7-membHabs-1 <=> 7-membHabs-2b + H2		//6
	addReaction(rxnV3, Reaction(5.33824E+04,2.26128E+00,7.52478E+00, sp::H2));      //H + 7-membHabs-1 <=> 7-membHabs-2b + H2		//7
	addReaction(rxnV3, Reaction(6.59893E+11,8.15588E-02,2.62478E+00, sp::None));    //7-membHabs-2b <=> 7-membHabs-3b				//8
	addReaction(rxnV3, Reaction(5.13229E+12,3.39661E-01,2.54015E+01, sp::None));    //7-membHabs-2b <=> 7-membHabs-3b				//9
	addReaction(rxnV3, Reaction(2.00638E+10,1.09912E+00,3.28066E+01, sp::None));    //7-membHabs-3b <=> 7-membHabs-4 + H			//10
	addReaction(rxnV3, Reaction(1.72899E+08,1.48867E+00,1.50412E+00, sp::H));       //7-membHabs-3b <=> 7-membHabs-4 + H			//11*/
	
    m_sType = R5ZZACR5; // sitetype
    m_name = "BY7 closure on R5ZZACR5"; // name of process
    m_ID = 51;
}
// Jump rate calculation
double L7_R5ZZACR5::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double L7_R5ZZACR5::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double L7_R5ZZACR5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
	double T = gp[gp.T];
    if (T <= 1000) return m_rate = 0.0;
	//Rate assuming PEQ approximation
	matrix<double> arr1(4, 4);
	boost::numeric::ublas::vector<double> arr2(4);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	arr1(0,0) = m_r[1] + m_r[2];
	arr1(0,1) = - m_r[3];
	arr1(1,0) = -m_r[2];
	arr1(1,1) = + m_r[3] + m_r[4];
	arr1(2,2) = +m_r[7] + m_r[9];
	arr1(2,3) = - m_r[8];
	arr1(3,2) = - m_r[9];
	arr1(3,3) = m_r[8] + m_r[10];
	arr2(0) = +m_r[0];
	arr2(1) = 0.0*(m_r[5]);
	arr2(2) = +m_r[6];
	arr2(3) = 0.0*(m_r[11]);
	
	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	double rate = (m_r[4] * arr2(1) + m_r[10] * arr2(3))*site_count;
	
	if (T>=2200) {
		matrix<double> arr1_carbene(7, 7);
		boost::numeric::ublas::vector<double> arr2_carbene(7);
		for (unsigned k = 0; k < arr1_carbene.size1(); ++k)
			for (unsigned l = 0; l < arr1_carbene.size2(); ++l)
				arr1_carbene(k, l) = 0.0;
		for (unsigned k = 0; k < arr2_carbene.size(); ++k)
			arr2_carbene(k) = 0.0;
		
		arr1_carbene(0,0) = m_r[13] + m_r[14];
		arr1_carbene(0,1) = -m_r[15];
		arr1_carbene(1,0) = -m_r[14];
		arr1_carbene(1,1) = m_r[15] + m_r[16];
		arr1_carbene(1,2) = -m_r[17];
		arr1_carbene(2,1) = -m_r[16];
		arr1_carbene(2,2) = m_r[17] + m_r[18];
		arr1_carbene(2,3) = -m_r[19];
		arr1_carbene(3,2) = -m_r[18];
		arr1_carbene(3,3) = m_r[19] + m_r[20];
		arr1_carbene(4,4) = m_r[23] + m_r[24];
		arr1_carbene(4,5) = -m_r[25];
		arr1_carbene(5,4) = -m_r[24];
		arr1_carbene(5,5) = m_r[25] + m_r[26];
		arr1_carbene(5,6) = -m_r[27];
		arr1_carbene(6,5) = -m_r[26];
		arr1_carbene(6,6) = m_r[27] + m_r[28];
		arr2_carbene(0) = +m_r[12];
		arr2_carbene(3) = 0.0*(m_r[21]);
		arr2_carbene(4) = +m_r[22];
		arr2_carbene(6) = 0.0*(m_r[29]);
		
		permutation_matrix<size_t> pm_carbene(arr1_carbene.size1());
		lu_factorize(arr1_carbene, pm_carbene);
		lu_substitute(arr1_carbene, pm_carbene, arr2_carbene);
		
		rate += (m_r[20] * arr2_carbene(3) + m_r[28] * arr2_carbene(6)) * site_count;
	}
		
	return m_rate = rate;
}

// ************************************************************
// ID52- ACR5R5R6 closure reaction
// ************************************************************
// Elementary rate constants, site type, process type and name
void L6_ACR5R5R6::initialise() {
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	//OLD METHOD
	/*addReaction(rxnV3, Reaction(9.24e07, 1.500, 9.646, sp::H));     // 0 - r1f
	addReaction(rxnV3, Reaction(9.60e04, 1.960, 9.021, sp::H2));    // 1 - r1b
	addReaction(rxnV3, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV3, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2b
	addReaction(rxnV3, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV3, Reaction(1.11e11, 0.658, 23.99, sp::None));  // 5 - r4f*/

	//Rates calculated by Angiras (unpublished as July2020)
	addReaction(rxnV3, Reaction(3.434E+07,1.879E+00,8.801E+00, sp::H)); //1
	addReaction(rxnV3, Reaction(2.644E+05,2.048E+00,7.507E+00, sp::H2)); //-1
	addReaction(rxnV3, Reaction(9.266E+11,1.131E-01,1.034E+00, sp::None)); //2
	addReaction(rxnV3, Reaction(1.037E+12,3.770E-01,3.743E+01, sp::None)); //-2
	addReaction(rxnV3, Reaction(4.695E+09,1.110E+00,2.157E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(2.793E+08,1.420E+00,4.145E+00, sp::H)); //-3
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	//Direct cyclisation
	addReaction(rxnV3, Reaction(1.998E+11,3.603E-01,6.303E+01, sp::None)); //2
	addReaction(rxnV3, Reaction(1.428E+12,3.177E-01,4.889E+00, sp::None)); //-2
	addReaction(rxnV3, Reaction(5.766E+09,1.110E+00,2.532E+01, sp::None)); //3
	addReaction(rxnV3, Reaction(4.133E+05,1.895E+00,1.011E+02, sp::H)); //-3

	m_sType = ACR5R5R6; // sitetype
	m_name = "ACR5R5R6 closure"; // name of process
	m_ID = 52;
}
double L6_ACR5R5R6::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	/*// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;

	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	return m_rate = 2 * m_r[5] * r_f* site_count; // Rate Equation*/

	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
	//Rate assuming PEQ approximation
	matrix<double> arr1(3, 3);
	boost::numeric::ublas::vector<double> arr2(3);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	// SS species are C*, ACR6-H, H-ACR6-H
	arr1(0,0) = m_r[1] + m_r[2] + m_r[7];
    arr1(0,1) = -m_r[3];
    arr1(1,0) = -m_r[2];
    arr1(1,1) = m_r[3] + m_r[4];
	arr1(2,2) = m_r[9] + m_r[10];

    arr2(0) = +m_r[0] + m_r[6];
    arr2(1) = 0.0;
	arr2(2) = +m_r[8];

	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	//std::cout << arr2 << std::endl;
	return m_rate = (m_r[4] * arr2(1) + m_r[10] * arr2(2) ) * site_count;
}

// ************************************************************
// ID53 - R7-bay closure at ACR5R5R6ZZ
// ************************************************************
// Elementary rate constants, site type, process type and name
void L7_ACR5R5R6ZZ::initialise() {
    // Adding elementary reactions. The rate for this jump process was assumed to be the same as L7_ACACR5.
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//M06-2X level.
	//Pathway A
	addReaction(rxnV3, Reaction(3.915e+07, 1.876e+00, 9.421e+00, sp::H));           		//7-memb2-start + H => 7-memb2-habs-1a + H2     //0
	addReaction(rxnV3, Reaction(5.369e+04, 2.275e+00, 5.583e+00, sp::H2));          		//7-memb2-habs-1a + H2 => 7-memb2-start + H     //1
	addReaction(rxnV3, Reaction(8.513e+11, 1.358e-01, 4.510e+00, sp::None));                //7-memb2-habs-1a => 7-memb2-habs-2a    //2
	addReaction(rxnV3, Reaction(3.523e+12, 2.926e-01, 2.567e+01, sp::None));                //7-memb2-habs-2a => 7-memb2-habs-1a    //3
	addReaction(rxnV3, Reaction(2.033e+10, 1.067e+00, 3.160e+01, sp::None));                //7-memb2-habs-2a => 7-memb2-end + H    //4
	addReaction(rxnV3, Reaction(1.033e+08, 1.495e+00, 2.895e+00, sp::H));           		//7-memb2-end + H => 7-memb2-habs-2a    //5
	//Pathway B
	addReaction(rxnV3, Reaction(3.091e+07, 1.891e+00, 9.308e+00, sp::H));           		//7-memb2-start + H => 7-memb2-habs-1b + H2     //6
	addReaction(rxnV3, Reaction(5.144e+04, 2.267e+00, 7.132e+00, sp::H2));          		//7-memb2-habs-1b + H2 => 7-memb2-start + H     //7
	addReaction(rxnV3, Reaction(7.041e+11, 1.839e-01, 1.034e+01, sp::None));                //7-memb2-habs-1b => 7-memb2-habs-2b    //8
	addReaction(rxnV3, Reaction(2.944e+12, 4.130e-01, 2.862e+01, sp::None));                //7-memb2-habs-2b => 7-memb2-habs-1b    //9
	addReaction(rxnV3, Reaction(1.861e+10, 1.136e+00, 2.957e+01, sp::None));                //7-memb2-habs-2b => 7-memb2-end + H    //10
	addReaction(rxnV3, Reaction(7.712e+07, 1.514e+00, 2.067e+00, sp::H));           		//7-memb2-end + H => 7-memb2-habs-2b    //11
	
	//Carbene route
	//Pathway A
	addReaction(rxnV3, Reaction(8.031e+10, 8.996e-01, 9.583e+01, sp::None));                //7-memb2-start => 7-memb2-carb-1a      //12
	addReaction(rxnV3, Reaction(4.398e+11, 3.592e-01, 3.385e+00, sp::None));                //7-memb2-carb-1a => 7-memb2-start      //13
	addReaction(rxnV3, Reaction(8.031e+11, 1.032e-02, 8.456e+00, sp::None));                //7-memb2-carb-1a => 7-memb2-carb-2a    //14
	addReaction(rxnV3, Reaction(1.897e+12, 2.234e-01, 1.779e+01, sp::None));                //7-memb2-carb-2a => 7-memb2-carb-1a    //15
	addReaction(rxnV3, Reaction(5.759e+11, 3.930e-01, 8.715e-01, sp::None));                //7-memb2-carb-2a => 7-memb2-carb-4a    //16
	addReaction(rxnV3, Reaction(1.052e+11, 9.049e-01, 5.350e+01, sp::None));                //7-memb2-carb-4a => 7-memb2-carb-2a    //17
	addReaction(rxnV3, Reaction(8.873e+10, 6.388e-01, 3.131e+01, sp::None));                //7-memb2-carb-4a => 7-memb2-carb-5a    //18
	addReaction(rxnV3, Reaction(1.728e+10, 7.118e-01, 6.065e+01, sp::None));                //7-memb2-carb-5a => 7-memb2-carb-4a    //19
	addReaction(rxnV3, Reaction(3.907e+09, 1.273e+00, 9.705e+01, sp::None));                //7-memb2-carb-5a => 7-memb2-end + H2   //20
	addReaction(rxnV3, Reaction(2.448e+05, 1.999e+00, 8.640e+01, sp::H2));          		//7-memb2-end + H2 => 7-memb2-carb-5a   //21
	//Pathway B
	addReaction(rxnV3, Reaction(1.061e+11, 7.990e-01, 8.426e+01, sp::None));                //7-memb2-start => 7-memb2-carb-1b      //22
	addReaction(rxnV3, Reaction(5.486e+11, 3.348e-01, 2.012e+00, sp::None));                //7-memb2-carb-1b => 7-memb2-start      //23
	addReaction(rxnV3, Reaction(1.000e+12, -1.397e-02, 3.568e+00, sp::None));               //7-memb2-carb-1b => 7-memb2-carb-2b    //24
	addReaction(rxnV3, Reaction(2.167e+12, 5.562e-01, 5.928e+01, sp::None));                //7-memb2-carb-2b => 7-memb2-carb-1b    //25
	addReaction(rxnV3, Reaction(3.063e+11, 8.244e-01, 6.356e+01, sp::None));                //7-memb2-carb-2b => 7-memb2-carb-4b    //26
	addReaction(rxnV3, Reaction(1.514e+10, 6.739e-01, 5.031e+01, sp::None));                //7-memb2-carb-4b => 7-memb2-carb-2b    //27
	addReaction(rxnV3, Reaction(4.479e+09, 7.144e-01, 2.710e+01, sp::None));                //7-memb2-carb-4b => 7-memb2-end + H2   //28
	addReaction(rxnV3, Reaction(2.332e+05, 1.742e+00, 5.510e+01, sp::H2));          		//7-memb2-end + H2 => 7-memb2-carb-4b   //29
	
	/*//B3LYP level.
	//Pathway A
	addReaction(rxnV3, Reaction(2.76732E+07,1.91279E+00,9.54190E+00, sp::H));		//H + 7-membHabs-1 <=> 7-membHabs-2a + H2		//0
	addReaction(rxnV3, Reaction(4.21185E+04,2.26366E+00,6.87815E+00, sp::H2));      //H + 7-membHabs-1 <=> 7-membHabs-2a + H2		//1
	addReaction(rxnV3, Reaction(4.70295E+11,1.42940E-01,4.72228E+00, sp::None));	//7-membHabs-2a <=> 7-membHabs-3a				//2
	addReaction(rxnV3, Reaction(1.47581E+12,3.66538E-01,2.73731E+01, sp::None));    //7-membHabs-2a <=> 7-membHabs-3a				//3
	addReaction(rxnV3, Reaction(6.42375E+09,1.09295E+00,3.21573E+01, sp::None));	//7-membHabs-3a <=> 7-membHabs-4 + H			//4
	addReaction(rxnV3, Reaction(1.69261E+08,1.52168E+00,1.63656E+00, sp::H));       //7-membHabs-3a <=> 7-membHabs-4 + H			//5
	//Pathway B
	addReaction(rxnV3, Reaction(2.84300E+07,1.90572E+00,9.53262E+00, sp::H));		//H + 7-membHabs-1 <=> 7-membHabs-2b + H2		//6
	addReaction(rxnV3, Reaction(5.33824E+04,2.26128E+00,7.52478E+00, sp::H2));      //H + 7-membHabs-1 <=> 7-membHabs-2b + H2		//7
	addReaction(rxnV3, Reaction(6.59893E+11,8.15588E-02,2.62478E+00, sp::None));    //7-membHabs-2b <=> 7-membHabs-3b				//8
	addReaction(rxnV3, Reaction(5.13229E+12,3.39661E-01,2.54015E+01, sp::None));    //7-membHabs-2b <=> 7-membHabs-3b				//9
	addReaction(rxnV3, Reaction(2.00638E+10,1.09912E+00,3.28066E+01, sp::None));    //7-membHabs-3b <=> 7-membHabs-4 + H			//10
	addReaction(rxnV3, Reaction(1.72899E+08,1.48867E+00,1.50412E+00, sp::H));       //7-membHabs-3b <=> 7-membHabs-4 + H			//11*/
	
    m_sType = ACR5R5R6ZZ; // sitetype
    m_name = "BY7 closure on ACR5R5R6ZZ"; // name of process
    m_ID = 53;
}
// Jump rate calculation
double L7_ACR5R5R6ZZ::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double L7_ACR5R5R6ZZ::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double L7_ACR5R5R6ZZ::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    // check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
	double T = gp[gp.T];
    if (T <= 1000) return m_rate = 0.0;
	//Rate assuming PEQ approximation
	matrix<double> arr1(4, 4);
	boost::numeric::ublas::vector<double> arr2(4);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	arr1(0,0) = m_r[1] + m_r[2];
	arr1(0,1) = - m_r[3];
	arr1(1,0) = -m_r[2];
	arr1(1,1) = + m_r[3] + m_r[4];
	arr1(2,2) = +m_r[7] + m_r[9];
	arr1(2,3) = - m_r[8];
	arr1(3,2) = - m_r[9];
	arr1(3,3) = m_r[8] + m_r[10];
	arr2(0) = +m_r[0];
	arr2(1) = 0.0*(m_r[5]);
	arr2(2) = +m_r[6];
	arr2(3) = 0.0*(m_r[11]);
	
	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	double rate = (m_r[4] * arr2(1) + m_r[10] * arr2(3))*site_count;
	
	if (T>=2200) {
		matrix<double> arr1_carbene(7, 7);
		boost::numeric::ublas::vector<double> arr2_carbene(7);
		for (unsigned k = 0; k < arr1_carbene.size1(); ++k)
			for (unsigned l = 0; l < arr1_carbene.size2(); ++l)
				arr1_carbene(k, l) = 0.0;
		for (unsigned k = 0; k < arr2_carbene.size(); ++k)
			arr2_carbene(k) = 0.0;
		
		arr1_carbene(0,0) = m_r[13] + m_r[14];
		arr1_carbene(0,1) = -m_r[15];
		arr1_carbene(1,0) = -m_r[14];
		arr1_carbene(1,1) = m_r[15] + m_r[16];
		arr1_carbene(1,2) = -m_r[17];
		arr1_carbene(2,1) = -m_r[16];
		arr1_carbene(2,2) = m_r[17] + m_r[18];
		arr1_carbene(2,3) = -m_r[19];
		arr1_carbene(3,2) = -m_r[18];
		arr1_carbene(3,3) = m_r[19] + m_r[20];
		arr1_carbene(4,4) = m_r[23] + m_r[24];
		arr1_carbene(4,5) = -m_r[25];
		arr1_carbene(5,4) = -m_r[24];
		arr1_carbene(5,5) = m_r[25] + m_r[26];
		arr1_carbene(5,6) = -m_r[27];
		arr1_carbene(6,5) = -m_r[26];
		arr1_carbene(6,6) = m_r[27] + m_r[28];
		arr2_carbene(0) = +m_r[12];
		arr2_carbene(3) = 0.0*(m_r[21]);
		arr2_carbene(4) = +m_r[22];
		arr2_carbene(6) = 0.0*(m_r[29]);
		
		permutation_matrix<size_t> pm_carbene(arr1_carbene.size1());
		lu_factorize(arr1_carbene, pm_carbene);
		lu_substitute(arr1_carbene, pm_carbene, arr2_carbene);
		
		rate += (m_r[20] * arr2_carbene(3) + m_r[28] * arr2_carbene(6)) * site_count;
	}
		
	return m_rate = rate;
}

// ************************************************************
// ID54- CH3 addition
// ************************************************************
// Elementary rate constants, site type, process type and name
void A_CH3::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    addReaction(rxnV, Reaction(2.50e14, 0, 16.00, sp::H));      //0 - r1f
    addReaction(rxnV, Reaction(3.40e9, 0.88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV, Reaction(2.10e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV, Reaction(3.49e39, -7.77, 13.35468, sp::H));  //4 - r3f
    addReaction(rxnV, Reaction(1.87e7, 1.787, 3.262, sp::C2H2));   //5 - r4f
	//addReaction(rxnV, Reaction(1.070E+65, -1.564E+01, 2.2720E+01, sp::CH3));         // A1- + CH3 -> A1CH3          	    - 9              - Richter, Howard 2005
	//addReaction(rxnV, Reaction(5.354E-09, 1.45, -104.4, sp::None));  //6				Keq: A1- + CH3 -> A1CH3             - Klippenstein 2007
	addReaction(rxnV, Reaction(2.331E+14, -0.283, -0.191, sp::CH3));  //6				kinf: A1- + CH3 -> A1CH3            - Klippenstein 2007
	//addReaction(rxnV, Reaction(6.528E+102, -20.22, 14.98, sp::CH3));  //8				k0: A1- + CH3 -> A1CH3              - Klippenstein 2007
	//addReaction(rxnV, Reaction(1.30E-03, 14.98, -3.11, sp::None));  //9				    Fcent: A1- + CH3 -> A1CH3           - Klippenstein 2007
	
    //addReaction(rxnV, Reaction(2.20e12, 0, 7.5, sp::O2));          //6 - r5f
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    addReaction(rxnV2, Reaction(2.50e14, 0, 16.00, sp::H));      //0 - r1f
    addReaction(rxnV2, Reaction(3.40e9, 0.88, 7.870, sp::H2));   //1 - r1b
    addReaction(rxnV2, Reaction(2.10e13, 0, 4.56937799, sp::OH));  //2 - r2f
    addReaction(rxnV2, Reaction(3.68e8, 1.139, 17.10, sp::H2O)); //3 - r2b
    addReaction(rxnV2, Reaction(2.18e35, -6.51, 11.53110048, sp::H));  //4 - r3f
    addReaction(rxnV2, Reaction(1.87e7, 1.787, 3.262, sp::C2H2));   //5 - r4f
	//addReaction(rxnV2, Reaction(1.070E+65, -1.564E+01, 2.2720E+01, sp::CH3));         // A1- + CH3 -> A1CH3          	    - 9              - Richter, Howard 2005
	//addReaction(rxnV, Reaction(5.354E-09, 1.45, -104.4, sp::None));  //6				Keq: A1- + CH3 -> A1CH3             - Klippenstein 2007
	addReaction(rxnV, Reaction(2.331E+14, -0.283, -0.191, sp::CH3));  //6				kinf: A1- + CH3 -> A1CH3            - Klippenstein 2007
	//addReaction(rxnV, Reaction(6.528E+102, -20.22, 14.98, sp::CH3));  //8				k0: A1- + CH3 -> A1CH3              - Klippenstein 2007
	//addReaction(rxnV, Reaction(1.30E-03, 14.98, -3.11, sp::None));  //9				    Fcent: A1- + CH3 -> A1CH3           - Klippenstein 2007
    //addReaction(rxnV2, Reaction(9.7e3, 2.42, 38.46338, sp::O2));          //6 - r5f
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	//------------Reactions for A3-4----------------------------
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4              - Forward
	//Using rates by On the low-temperature limit of HACA. Frenklach et al. 2018
	addReaction(rxnV3, Reaction(1.190E+22, -2.450E+00, 1.889E+01, sp::C2H2));         // A3* + C2H2 -> A4 + H              	- 5              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(1.060E+14, -4.900E-01, 8.204E+00, sp::C2H2));         // A3* + H -> A4 + H              	- 6              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 7              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 8              - Frenklach et al. 2018
	//addReaction(rxnV3, Reaction(1.070E+65, -1.564E+01, 2.2720E+01, sp::CH3));         // A1- + CH3 -> A1CH3          	    - 9              - Richter, Howard 2005
	//addReaction(rxnV, Reaction(5.354E-09, 1.45, -104.4, sp::None));  //9				Keq: A1- + CH3 -> A1CH3             - Klippenstein 2007
	addReaction(rxnV, Reaction(2.331E+14, -0.283, -0.191, sp::CH3));  //9				kinf: A1- + CH3 -> A1CH3            - Klippenstein 2007
	//addReaction(rxnV, Reaction(6.528E+102, -20.22, 14.98, sp::CH3));  //11				k0: A1- + CH3 -> A1CH3              - Klippenstein 2007
	//addReaction(rxnV, Reaction(1.30E-03, 14.98, -3.11, sp::None));  //12				Fcent: A1- + CH3 -> A1CH3           - Klippenstein 2007

    m_sType = benz; // sitetype
    m_name = "CH3 addition"; // name of process
    m_ID = 54;
}
// Jump rate calculation
double A_CH3::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[6]*r_f* site_count/6.0; // Rate Equation
}
double A_CH3::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate0p0267(gp, pah_st);
}
double A_CH3::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]+m_r[8]+m_r[9]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = m_r[9]*r_f* site_count/6.0; // Rate Equation
}

// ************************************************************
// ID55- CH3 desorption
// ************************************************************
// Elementary rate constants, site type, process type and name
void D_CH3::initialise() {
    // Adding elementary reactions
    // 0.0267 atm
    rxnvector& rxnV = m_rxnvector0p0267;
    //addReaction(rxnV, Reaction(3.52751561e79, -1.71557169e+01, 105.3773, sp::None));      // A1- + CH3 -> A1CH3         	    - 9              - Richter, Howard 2005 OLD
	//------------Reactions for A3-4----------------------------
	addReaction(rxnV, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4              - Forward
	//Using rates by On the low-temperature limit of HACA. Frenklach et al. 2018
	addReaction(rxnV, Reaction(1.190E+22, -2.450E+00, 1.889E+01, sp::C2H2));         // A3* + C2H2 -> A4 + H              	- 5              - Frenklach et al. 2018
	addReaction(rxnV, Reaction(1.060E+14, -4.900E-01, 8.204E+00, sp::C2H2));         // A3* + H -> A4 + H              	- 6              - Frenklach et al. 2018
	addReaction(rxnV, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 7              - Frenklach et al. 2018
	addReaction(rxnV, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 8              - Frenklach et al. 2018
	//addReaction(rxnV, Reaction(1.070E+65, -1.564E+01, 2.2720E+01, sp::CH3));         // A1- + CH3 -> A1CH3          	    - 9              - Richter, Howard 2005
	addReaction(rxnV, Reaction(1.950E+27, -3.160, 107.447, sp::None));  //9				kinf: A1CH3 -> A1- + CH3            - Klippenstein 2007
	addReaction(rxnV, Reaction(7.054562E-01, -22.966, 122.08, sp::None));  //10			k0: A1CH3 -> A1- + CH3              - Klippenstein 2007
    // 0.12 atm
    rxnvector& rxnV2 = m_rxnvector0p12;
    //addReaction(rxnV2, Reaction(3.52751561e79, -1.71557169e+01, 105.3773, sp::None));     // A1- + CH3 -> A1CH3         	    - 9              - Richter, Howard 2005 OLD
	//------------Reactions for A3-4----------------------------
	addReaction(rxnV2, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV2, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV2, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV2, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV2, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4              - Forward
	//Using rates by On the low-temperature limit of HACA. Frenklach et al. 2018
	addReaction(rxnV2, Reaction(1.190E+22, -2.450E+00, 1.889E+01, sp::C2H2));         // A3* + C2H2 -> A4 + H              	- 5              - Frenklach et al. 2018
	addReaction(rxnV2, Reaction(1.060E+14, -4.900E-01, 8.204E+00, sp::C2H2));         // A3* + H -> A4 + H              	- 6              - Frenklach et al. 2018
	addReaction(rxnV2, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 7              - Frenklach et al. 2018
	addReaction(rxnV2, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 8              - Frenklach et al. 2018
	//addReaction(rxnV2, Reaction(1.070E+65, -1.564E+01, 2.2720E+01, sp::CH3));         // A1- + CH3 -> A1CH3          	    - 9              - Richter, Howard 2005
	addReaction(rxnV2, Reaction(1.950E+27, -3.160, 107.447, sp::None));  //9				kinf: A1CH3 -> A1- + CH3            - Klippenstein 2007
	addReaction(rxnV2, Reaction(7.054562E-01, -22.966, 122.08, sp::None));  //10			k0: A1CH3 -> A1- + CH3              - Klippenstein 2007
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	//------------Reactions for A3-4----------------------------
	addReaction(rxnV3, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV3, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4              - Forward
	//Using rates by On the low-temperature limit of HACA. Frenklach et al. 2018
	addReaction(rxnV3, Reaction(1.190E+22, -2.450E+00, 1.889E+01, sp::C2H2));         // A3* + C2H2 -> A4 + H              	- 5              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(1.060E+14, -4.900E-01, 8.204E+00, sp::C2H2));         // A3* + H -> A4 + H              	- 6              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 7              - Frenklach et al. 2018
	addReaction(rxnV3, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H        	- 8              - Frenklach et al. 2018
	//addReaction(rxnV3, Reaction(1.070E+65, -1.564E+01, 2.2720E+01, sp::CH3));         // A1- + CH3 -> A1CH3          	    - 9              - Richter, Howard 2005
	addReaction(rxnV3, Reaction(1.950E+27, -3.160, 107.447, sp::None));  //9				kinf: A1CH3 -> A1- + CH3            - Klippenstein 2007
	addReaction(rxnV3, Reaction(7.054562E-01, -22.966, 122.08, sp::None));  //10			k0: A1CH3 -> A1- + CH3              - Klippenstein 2007

    m_sType = Methyl; // sitetype
    m_name = "CH3 desorption"; // name of process
    m_ID = 55;
}
// Jump rate calculation
double D_CH3::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st); 
}
double D_CH3::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double D_CH3::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper:Klippenstein2007, reverse notation retrieved from Blanquart2008
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]+m_r[8]+m_r[9]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
	double T = gp[gp.T];
	double P = gp[gp.P]; //Pressure in Pa
	double M_eff1 = P / 8.314E+06 / T;
	double M = M_eff1; // + (1 * gp[gp.H2]) + (5 * gp[gp.H2O]) + (1 * gp[gp.CH4]) + (0.5*gp[gp.CO]) + (1 * gp[gp.CO2]) + (2 * gp[gp.C2H6]);
	double k_o1 = m_r[10];
	double k_inf1 = m_r[9];
	double F_cent1 = (1.0-0.7054562) * exp(-T / 9.999989E+09) + 0.7054562*exp(-T / 459.9180) + exp(-8.213938E+09 / T);
	double logF_cent1 = log10(F_cent1);
	double k_eff = k_inf1 / (1 + (k_inf1 / (k_o1*M)));
	double N01 = 0.75 - 1.27*logF_cent1;
	double C01 = -0.4 - 0.67*logF_cent1;
	double logP_r1 = log10(k_o1*M / k_inf1);
	double F1 = exp(logF_cent1 / (1 + pow(((logP_r1 + C01) / (N01 - 0.14*(logP_r1 + C01))), 2)));
	k_eff *= F1;                                        
	
	//double for_rate = m_r[9]*r_f;
	//double rev_rate = for_rate * 100.0;
    return m_rate = k_eff* site_count; // Rate Equation
	//return setRate0p0267(gp, pah_st);
}


// ************************************************************
// ID56- Oxidation of R5R6 site.
// ************************************************************
// Elementary rate constants, site type, process type and name
void O5R_R5R6::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	//Using rates by Frenklach 2018
	//"Detailed, sterically-resolved modeling of soot oxidation: Role of O atoms, interplay with particle nanostructure, and emergence of inner particle burning"
	addReaction(rxnV3, Reaction(3.450E+11, 0.000E+00, 13422.0*8.314/4.184/1000.0, sp::O2));         // S90
	addReaction(rxnV3, Reaction(6.290E+12/17.76148, 0.505E+00, 306.0*8.314/4.184/1000.0, sp::O));         // S104

    m_sType = R5R6; // sitetype
    m_name = "R5R6 oxidation"; // name of process
    m_ID = 56;
}
// Jump rate calculation
double O5R_R5R6::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double O5R_R5R6::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double O5R_R5R6::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r_f*(m_r[6] + m_r[7]) * site_count/2.0; // Rate Equation*/
	return m_rate = (m_r[6] + m_r[7]) * site_count/2.0; // Rate Equation
	//return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID57- Oxidation of R5R6ZZ site.
// ************************************************************
// Elementary rate constants, site type, process type and name
void O5R_R5R6ZZ::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	//Using rates by Frenklach 2018
	//"Detailed, sterically-resolved modeling of soot oxidation: Role of O atoms, interplay with particle nanostructure, and emergence of inner particle burning"
	addReaction(rxnV3, Reaction(3.450E+11, 0.000E+00, 13422.0*8.314/4.184/1000.0, sp::O2));         // S90
	addReaction(rxnV3, Reaction(6.290E+12/17.76148, 0.505E+00, 306.0*8.314/4.184/1000.0, sp::O));         // S104

    m_sType = R5R6ZZ; // sitetype
    m_name = "Oxidation of R5R6ZZ"; // name of process
    m_ID = 57;
}
// Jump rate calculation
double O5R_R5R6ZZ::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double O5R_R5R6ZZ::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double O5R_R5R6ZZ::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r_f*(m_r[6] + m_r[7]) * site_count/2.0; // Rate Equation*/
	return m_rate = (m_r[6] + m_r[7]) * site_count/2.0; // Rate Equation
	//return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID58- Oxidation of R5R6AC site.
// ************************************************************
// Elementary rate constants, site type, process type and name
void O5R_R5R6AC::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	//Using rates by Frenklach 2018
	//"Detailed, sterically-resolved modeling of soot oxidation: Role of O atoms, interplay with particle nanostructure, and emergence of inner particle burning"
	addReaction(rxnV3, Reaction(3.450E+11, 0.000E+00, 13422.0*8.314/4.184/1000.0, sp::O2));         // S90
	addReaction(rxnV3, Reaction(6.290E+12/17.76148, 0.505E+00, 306.0*8.314/4.184/1000.0, sp::O));         // S104

    m_sType = R5R6AC; // sitetype
    m_name = "Oxidation of R5R6AC"; // name of process
    m_ID = 58;
}
// Jump rate calculation
double O5R_R5R6AC::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double O5R_R5R6AC::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double O5R_R5R6AC::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r_f*(m_r[6] + m_r[7]) * site_count/2.0; // Rate Equation*/
	return m_rate = (m_r[6] + m_r[7]) * site_count/2.0;
	//return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID59- Oxidation of R5R6BY5 site.
// ************************************************************
// Elementary rate constants, site type, process type and name
void O5R_R5R6BY5::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	//Using rates by Frenklach 2018
	//"Detailed, sterically-resolved modeling of soot oxidation: Role of O atoms, interplay with particle nanostructure, and emergence of inner particle burning"
	addReaction(rxnV3, Reaction(3.450E+11, 0.000E+00, 13422.0*8.314/4.184/1000.0, sp::O2));         // S90
	addReaction(rxnV3, Reaction(6.290E+12/17.76148, 0.505E+00, 306.0*8.314/4.184/1000.0, sp::O));         // S104

    m_sType = R5R6BY5; // sitetype
    m_name = "Oxidation of R5R6BY5"; // name of process
    m_ID = 59;
}
// Jump rate calculation
double O5R_R5R6BY5::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double O5R_R5R6BY5::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double O5R_R5R6BY5::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r_f*(m_r[6] + m_r[7]) * site_count/2.0; // Rate Equation*/
	return m_rate = (m_r[6] + m_r[7]) * site_count/2.0;
	//return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID60- Oxidation of R5R6FER site.
// ************************************************************
// Elementary rate constants, site type, process type and name
void O5R_R5R6FER::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	//Using rates by Frenklach 2018
	//"Detailed, sterically-resolved modeling of soot oxidation: Role of O atoms, interplay with particle nanostructure, and emergence of inner particle burning"
	addReaction(rxnV3, Reaction(3.450E+11, 0.000E+00, 13422.0*8.314/4.184/1000.0, sp::O2));         // S90
	addReaction(rxnV3, Reaction(6.290E+12/17.76148, 0.505E+00, 306.0*8.314/4.184/1000.0, sp::O));         // S104

    m_sType = R5R6FER; // sitetype
    m_name = "Oxidation of R5R6FER"; // name of process
    m_ID = 60;
}
// Jump rate calculation
double O5R_R5R6FER::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double O5R_R5R6FER::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double O5R_R5R6FER::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r_f*(m_r[6] + m_r[7]) * site_count/2.0; // Rate Equation*/
	return m_rate = (m_r[6] + m_r[7]) * site_count/2.0;
	//return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID61- Oxidation of R5R6ZZR site.
// ************************************************************
// Elementary rate constants, site type, process type and name
void O5R_R5R6ZZR::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	//Using rates by Frenklach 2018
	//"Detailed, sterically-resolved modeling of soot oxidation: Role of O atoms, interplay with particle nanostructure, and emergence of inner particle burning"
	addReaction(rxnV3, Reaction(3.450E+11, 0.000E+00, 13422.0*8.314/4.184/1000.0, sp::O2));         // S90
	addReaction(rxnV3, Reaction(6.290E+12/17.76148, 0.505E+00, 306.0*8.314/4.184/1000.0, sp::O));         // S104

    m_sType = R5R6ZZR; // sitetype
    m_name = "Oxidation of R5R6ZZR"; // name of process
    m_ID = 61;
}
// Jump rate calculation
double O5R_R5R6ZZR::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double O5R_R5R6ZZR::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double O5R_R5R6ZZR::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r_f*(m_r[6] + m_r[7]) * site_count/2.0; // Rate Equation*/
	return m_rate = (m_r[6] + m_r[7]) * site_count/2.0;
	//return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID62- Oxidation of R5R6ACR site.
// ************************************************************
// Elementary rate constants, site type, process type and name
void O5R_R5R6ACR::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	//Using rates by Frenklach 2018
	//"Detailed, sterically-resolved modeling of soot oxidation: Role of O atoms, interplay with particle nanostructure, and emergence of inner particle burning"
	addReaction(rxnV3, Reaction(3.450E+11, 0.000E+00, 13422.0*8.314/4.184/1000.0, sp::O2));         // S90
	addReaction(rxnV3, Reaction(6.290E+12/17.76148, 0.505E+00, 306.0*8.314/4.184/1000.0, sp::O));         // S104

    m_sType = R5R6ACR; // sitetype
    m_name = "Oxidation of R5R6ACR"; // name of process
    m_ID = 62;
}
// Jump rate calculation
double O5R_R5R6ACR::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double O5R_R5R6ACR::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double O5R_R5R6ACR::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r_f*(m_r[6] + m_r[7]) * site_count/2.0; // Rate Equation*/
	return m_rate = (m_r[6] + m_r[7]) * site_count/2.0;
	//return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID63- Oxidation of R5R6FER5R6 site.
// ************************************************************
// Elementary rate constants, site type, process type and name
void O5R_R5R6FER5R6::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	//Using rates by Frenklach 2018
	//"Detailed, sterically-resolved modeling of soot oxidation: Role of O atoms, interplay with particle nanostructure, and emergence of inner particle burning"
	addReaction(rxnV3, Reaction(3.450E+11, 0.000E+00, 13422.0*8.314/4.184/1000.0, sp::O2));         // S90
	addReaction(rxnV3, Reaction(6.290E+12/17.76148, 0.505E+00, 306.0*8.314/4.184/1000.0, sp::O));         // S104

    m_sType = R5R6FER5R6; // sitetype
    m_name = "Oxidation of R5R6FER5R6"; // name of process
    m_ID = 63;
}
// Jump rate calculation
double O5R_R5R6FER5R6::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double O5R_R5R6FER5R6::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double O5R_R5R6FER5R6::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r_f*(m_r[6] + m_r[7]) * site_count; // Rate Equation*/
	return m_rate = (m_r[6] + m_r[7]) * site_count;
	//return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID64- Oxidation of R5R6ZZR5R6 site.
// ************************************************************
// Elementary rate constants, site type, process type and name
void O5R_R5R6ZZR5R6::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	//Using rates by Frenklach 2018
	//"Detailed, sterically-resolved modeling of soot oxidation: Role of O atoms, interplay with particle nanostructure, and emergence of inner particle burning"
	addReaction(rxnV3, Reaction(3.450E+11, 0.000E+00, 13422.0*8.314/4.184/1000.0, sp::O2));         // S90
	addReaction(rxnV3, Reaction(6.290E+12/17.76148, 0.505E+00, 306.0*8.314/4.184/1000.0, sp::O));         // S104

    m_sType = R5R6ZZR5R6; // sitetype
    m_name = "Oxidation of R5R6ZZR5R6"; // name of process
    m_ID = 64;
}
// Jump rate calculation
double O5R_R5R6ZZR5R6::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double O5R_R5R6ZZR5R6::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double O5R_R5R6ZZR5R6::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r_f*(m_r[6] + m_r[7]) * site_count; // Rate Equation*/
	return m_rate = (m_r[6] + m_r[7]) * site_count;
	//return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID65- Oxidation of R5R6ACR5R6 site.
// ************************************************************
// Elementary rate constants, site type, process type and name
void O5R_R5R6ACR5R6::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(4.890e+09, 1.508e+00, 1.9862e+01, sp::H));          // A2R5 + H <=> A2R5- + H2      - 0 Hou2017
	addReaction(rxnV3, Reaction(5.0677e+04, 2.4449e+00, 4.5197e+00, sp::H2));       // A2R5 + H <=> A2R5- + H2      - 1 Hou2017
	addReaction(rxnV3, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // - 2              - Forward
	addReaction(rxnV3, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // - 3              - Backward
	addReaction(rxnV3, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // - 4              - Forward
	addReaction(rxnV3, Reaction(1.235E+07, 1.530E+00, 9.311E+00, sp::C2H2));  		 // - 5  Whitesides and Frenklach2010
	//Using rates by Frenklach 2018
	//"Detailed, sterically-resolved modeling of soot oxidation: Role of O atoms, interplay with particle nanostructure, and emergence of inner particle burning"
	addReaction(rxnV3, Reaction(3.450E+11, 0.000E+00, 13422.0*8.314/4.184/1000.0, sp::O2));         // S90
	addReaction(rxnV3, Reaction(6.290E+12/17.76148, 0.505E+00, 306.0*8.314/4.184/1000.0, sp::O));         // S104

    m_sType = R5R6ACR5R6; // sitetype
    m_name = "Oxidation of R5R6ACR5R6"; // name of process
    m_ID = 65;
}
// Jump rate calculation
double O5R_R5R6ACR5R6::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double O5R_R5R6ACR5R6::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double O5R_R5R6ACR5R6::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: On the low-temperature limit of HACA. Frenklach et al. 2018
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    /*double r_denom = (m_r[1]+m_r[3]+m_r[4]+m_r[5]+m_r[6]+m_r[7]);
    double r_f; // radical fraction
    if(r_denom>0) {
        r_f = (m_r[0]+m_r[2])/r_denom; 
        r_f = r_f/(r_f+1.0);
    }
    else r_f=0;
    return m_rate = r_f*(m_r[6] + m_r[7]) * site_count; // Rate Equation*/
	return m_rate = (m_r[6] + m_r[7]) * site_count;
	//return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID66 - Embedded 5-member ring migration to ZZ (both sites accesible)
// ************************************************************

void M5R_ACR5_ZZ_ZZ::initialise() {
	// Adding elementary reactions
	// 0.0267 atm
	rxnvector& rxnV = m_rxnvector0p0267;
	//------------Reactions for A3-4----------------------------
	addReaction(rxnV, Reaction(4.570E+08, 1.880E+00, 1.4839E+01, sp::H));           // A3 + H <=> A3-4 + H2                - 0              - Forward
	addReaction(rxnV, Reaction(1.690E+04, 2.6200E+00, 4.559E+00, sp::H2));          // A3 + H <=> A3-4 + H2                - 1              - Backward
	addReaction(rxnV, Reaction(5.190E+03, 3.040E+00, 3.675E+00, sp::OH));           // A3 + OH <=> A3-4 + H2O              - 2              - Forward
	addReaction(rxnV, Reaction(5.590E+00, 3.573E+00, 8.659E+00, sp::H2O));          // A3 + OH <=> A3-4 + H2O              - 3              - Backward
	addReaction(rxnV, Reaction(4.170E+13, 1.500E-01, 0.000E+00, sp::H));            // A3* + H -> A3              			- 4   
	addReaction(rxnV, Reaction(4.240E+14,  2.500E-02, 3.308E+01, sp::C2H2));         // A3* + C2H2 -> A3C2H + H            - 5              - Frenklach et al. 2018
	addReaction(rxnV, Reaction(7.640E-02,  3.950E+00, 1.6495E+01, sp::C2H2));        // A3* + C2H2 -> A3C2H + H        	- 6              - Frenklach et al. 2018
	addReaction(rxnV, Reaction(4.960e+11, 7.550e-01, 5.000e+01, sp::None));  		// Violi2005. 							- 7
	addReaction(rxnV, Reaction(5.188e+13, 0.000E+00, 3.340e+01, sp::None));  		// Whitesides2009 system b. 			- 8
	//Assumes the migration of embedded rings is similar to migration to the edge	
	
	//ABF and old
	/*addReaction(rxnV, Reaction(1.74e08, 1.740, 9.370, sp::H));     // 0 - r1f
	addReaction(rxnV, Reaction(3.90e12, 0, 11.00, sp::H2));    // 1 - r1b
	addReaction(rxnV, Reaction(1.00e10, 0.734, 1.430, sp::OH));    // 2 - r2f
	addReaction(rxnV, Reaction(3.68e08, 1.139, 17.10, sp::H2O));   // 3 - r2
	addReaction(rxnV, Reaction(2.00e13, 0, 0, sp::H));     // 4 - r3f
	addReaction(rxnV, Reaction(4.96e11, 0.755, 50, sp::None));  // 5*/
	//
	// 0.12 atm
	rxnvector& rxnV2 = m_rxnvector0p12;
	rxnV2 = rxnV;
	//
	// 1 atm
	rxnvector& rxnV3 = m_rxnvector1;
	rxnV3 = rxnV;
	//
	m_sType = MIGR2; // sitetype
	m_name = "ACR5 migration both sides"; // name of process
	m_ID = 66;
}
// Jump rate calculation
double M5R_ACR5_ZZ_ZZ::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	// check if site count is zero
	double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
	if (site_count == 0) return m_rate = 0;
	// calculate rate
	double r_denom = (m_r[1] + m_r[3] + m_r[4] + m_r[5] + m_r[6] + m_r[7]);
	double r_f; // radical fraction 
	if (r_denom>0) {
		r_f = (m_r[0] + m_r[2]) / r_denom;
	}
	else r_f = 0;
	//return m_rate = m_r[7] * r_f*site_count; // Rate Equation //Violi
	return m_rate = m_r[8] * r_f*site_count * 2.0; // Rate Equation //Whitesides
}
double M5R_ACR5_ZZ_ZZ::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate0p0267(gp, pah_st);
}
double M5R_ACR5_ZZ_ZZ::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate0p0267(gp, pah_st);
}

// ************************************************************
// ID69- R5R7 pair edge healing.
// ************************************************************
// Elementary rate constants, site type, process type and name
void MR5R7_edge::initialise() {
    // Adding elementary reactions
    // 1 atm
    rxnvector& rxnV3 = m_rxnvector1;
	//ABF - Frenklach and Whitesides2010
	//H abstraction
	addReaction(rxnV3, Reaction(3.48e+08, 1.49, 3.44, sp::H));          			// BP + H -> BP1		0 
	addReaction(rxnV3, Reaction(2.34e+09, 1.18, 30.77, sp::None));       			// BP1 -> BP + H		1 
	addReaction(rxnV3, Reaction(2.12e+10, 0.67, 40.32, sp::None));          		// BP1 -> BP2			2 
	addReaction(rxnV3, Reaction(8.39e+11, 0.29, 18.61, sp::None));          		// BP2 -> BP1			3 
	addReaction(rxnV3, Reaction(1.03e+12, 0.15, 9.25, sp::None));          			// BP2 -> BP3			4
	addReaction(rxnV3, Reaction(1.87e+12, 0.32, 15.04, sp::None));           		// BP3 -> BP2			5
	addReaction(rxnV3, Reaction(6.65e+12, -0.13, 4.09, sp::None));          		// BP3 -> BP4			6 
	addReaction(rxnV3, Reaction(1.52e+12, 0.00, 34.08, sp::None));           		// BP4 -> BP3			7
	addReaction(rxnV3, Reaction(1.13e+10, 1.00, 15.65, sp::None));          		// BP4 -> BF + H		8 
	addReaction(rxnV3, Reaction(2.82e+08, 1.36, 35.44, sp::H));           			// BF + H -> BP4		9


    m_sType = R5R7; // sitetype
    m_name = "R5R7 pair edge healing"; // name of process
    m_ID = 67;
}
// Jump rate calculation
double MR5R7_edge::setRate0p0267(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	return setRate1(gp, pah_st);
}
double MR5R7_edge::setRate0p12(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
    return setRate1(gp, pah_st);
}
double MR5R7_edge::setRate1(const KMCGasPoint& gp, PAHProcess& pah_st/*, const double& time_now*/) {
	//Rates updated according to paper: A density functional theory study on the kinetics of seven-member ring formation in polyaromatic hydrocarbons 2020
	// check if site count is zero
    double site_count = ((double)pah_st.getSiteCount(m_sType));
	//double site_count = 1; // Site count
    if(site_count==0) return m_rate=0;
    // calculate rate
    //Rate assuming PEQ approximation
	matrix<double> arr1(4, 4);
	boost::numeric::ublas::vector<double> arr2(4);
	for (unsigned k = 0; k < arr1.size1(); ++k)
		for (unsigned l = 0; l < arr1.size2(); ++l)
			arr1(k, l) = 0.0;
	for (unsigned k = 0; k < arr2.size(); ++k)
		arr2(k) = 0.0;
	
	arr1(0,0) = m_r[1] + m_r[2];
	arr1(0,1) = - m_r[3];
	arr1(1,0) = -m_r[2];
	arr1(1,1) = + m_r[3] + m_r[4];
	arr1(1,2) = -m_r[5];
	arr1(2,1) = -m_r[4];
	arr1(2,2) = + m_r[5] + m_r[6];
	arr1(2,3) = -m_r[7];
	arr1(3,2) = -m_r[6];
	arr1(3,3) = + m_r[7] + m_r[8];

	arr2(0) = +m_r[0];
	arr2(1) = 0.0;
	arr2(2) = 0.0;
	arr2(3) = 0.0;
	
	permutation_matrix<size_t> pm(arr1.size1());
	lu_factorize(arr1, pm);
	lu_substitute(arr1, pm, arr2);
	
	double rate = (m_r[8] * arr2(3))*site_count;
	return m_rate = rate;
}