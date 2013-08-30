/*!
 * @author     William Menz (wjm34)
 * @file       kmc_pah.cpp
 * @brief      Basic unit testing of the KMC model
 */


#include "swp_kmc_jump_process.h"
#include "swp_kmc_mech.h"
#include "swp_kmc_pah_structure.h"
#include "swp_kmc_pah_process.h"
#include "swp_kmc_reaction.h"
#include "swp_kmc_simulator.h"
#include "swp_kmc_structure_comp.h"
#include "swp_kmc_typedef.h"
#include "swp_kmc_gaspoint.h"

namespace KMC = Sweep::KMC_ARS;
using std::cerr;
using std::endl;
using std::cout;

/*!
 * Compare properties of two PAHs, to see if we get the same thing
 *
 * @param lhs    PAH Structure 1
 * @param rhs    PAH Structure 2
 * @return       Are they the same?
 */
bool compare(const KMC::PAHStructure &lhs, const KMC::PAHStructure &rhs) {
	bool ans(true);
	if (lhs.numofC() != rhs.numofC()) ans = false;
	if (lhs.numofH() != rhs.numofH()) ans = false;
	if (lhs.numofRings() != rhs.numofRings()) ans = false;
	if (lhs.numofRings5() != rhs.numofRings5()) ans = false;
	if (lhs.numofEdgeC() != rhs.numofEdgeC()) ans = false;
	if (lhs.numofSite() != rhs.numofSite()) ans = false;
	return ans;
}

int main(int argc, char *argv[])
{
	// Introduce ourselves
	cout << "KMC unit testing." << endl;

	// Create some basic objects
	Sweep::rng_type rng(1337);

	// Create some PAHs!
	KMC::PAHStructure pah1, pah2, pah3;
	pah1.initialise(KMC::BENZENE_C);
	pah2.initialise(KMC::NAPHTHALENE_C);
	pah3.initialise(KMC::PYRENE_C);

	// Now create a mechanism with some processes to adjust these PAHs
	KMC::PAHProcess pah_process;
	KMC::G6R_FE proc_g6r_fe; proc_g6r_fe.initialise();
	// Create some oxdiation processes
	KMC::O6R_FE3_OH proc_06r_fe3_oh; proc_06r_fe3_oh.initialise();
	KMC::O6R_FE3_O2 proc_06r_fe3_o2; proc_06r_fe3_o2.initialise();
	KMC::O6R_FE_HACA_OH proc_06r_fe_haca_oh; proc_06r_fe_haca_oh.initialise();
	KMC::O6R_FE_HACA_O2 proc_06r_fe_haca_o2; proc_06r_fe_haca_o2.initialise();

	// Do G6R_FE on A1, should get napthalene
	pah_process.setPAH(pah1);
	pah_process.performProcess(proc_g6r_fe, rng);
	if (! compare(pah1, pah2)) {
		cerr << "Failed comparing A1 post G6R_FE with A2" << endl;
		return 1;
	}

	// Test oxidation of pyrene
	// Do O6R_FE_HACA_OH then O6R_FE3_OH, should get napthalene
	pah_process.setPAH(pah3);
	pah_process.performProcess(proc_06r_fe_haca_oh, rng);
	pah_process.performProcess(proc_06r_fe3_oh, rng);
	if (! compare(pah2, pah3)) {
		cerr << "Failed comparing A4 post O6R_FE_HACA_OH, O6R_FE3_OH with A2" << endl;
		return 1;
	}

	// Do O6R_FE_HACA_O2 then O6R_FE3_O2, should get napthalene
	pah3.initialise(KMC::PYRENE_C);
	pah_process.performProcess(proc_06r_fe_haca_o2, rng);
	pah_process.performProcess(proc_06r_fe3_o2, rng);
	if (! compare(pah2, pah3)) {
		cerr << "Failed comparing A4 post O6R_FE_HACA_O2, O6R_FE3_O2 with A2" << endl;
		return 1;
	}

	// If we got to here, hopefully everything is working
	cout << "All tests passed! :D" << endl;

	return 0;
}
