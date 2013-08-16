//#include <boost/program_options.hpp>
//namespace po = boost::program_options;

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// Sweep includes for this file
#include "swp_particle_model.h"
#include "swp_bintree_primary.h"
#include "swp_pri_list.hpp"
#include "swp_pri_list_node.h"
#include "swp_params.h"
#include "swp_sintering_model.h"

template <typename T>
bool within_tolerance(T check, T obtained, double tol) {
	bool ans(false);
	if (std::abs(check - obtained) / check < tol) ans = true;
	return ans;
}

void fail() {
	std::cerr << "Test failed." << std::endl;
}

using namespace Sweep;

int main(int argc, char* argv[])
{
	std::cout << "sweep: Particle Models test" << std::endl;

	/*
	 * ###################################################################
	 * Set up of basic objects needed for the testing
	 */

	// First create particle components and component vector
	// molwt, density, min_value, name
	Component* c = new Component(Sweep::NA, 1.0, 1.0, "dummy");

	// Now create the particle model
	ParticleModel pm;
	pm.AddComponent(*c); // Takes ownership and clears mem

	// Create a sintering model
	Processes::SinteringModel sm;
	sm.SetA(1.0);        // tau_s = 1.0 s
	sm.SetType(Processes::SinteringModel::Constant);

	// Create a Cell (needed for sintering)
	const Sprog::SpeciesPtrVector species;
	pm.SetSpecies(species);
	Cell cell(pm, true);

	// Set basic properties of the model, used by particle models
	pm.SetBinTreeCoalThresh(1.0);	// Coalescence threshold
	pm.SetFractDim(1.8);			// Fractal dimension of aggregates

	// Typedef for easier use
	typedef AggModels::PrimaryList<AggModels::PrimaryListNode> DefaultParticle;
	typedef boost::mt19937 rng_type;

	// Create random number generator
	rng_type rng(1);

	// Set the relative tolerance for testing our doubles
	double tolerance(1.0e-5);

	/*
	 * ###################################################################
	 * Main testing of particle models
	 */
	// Create a basic computational particle
	fvector fv(1, 1.0);
	DefaultParticle part_1(0.0, pm);
	part_1.SetComposition(fv);

	// Create a clone of part_1 and adjust it
	DefaultParticle part_2(part_1);
	fvector dvalues(0);
	part_1.Adjust(fv, dvalues, rng, 1);

	// Check the adjusted particle
	std::cout << "Checking adjusted particle..." << std::endl;
	if (part_1.NumberOfPrimaries() != 1 || part_1.Volume() != 2.0) {
		fail();
		return 1;
	}

	part_1.Coagulate(part_2, rng);

	// Create a clone of part_1 again
	DefaultParticle part_3(part_1);

	// And let's go crazy and coagulate again!
	part_3.Coagulate(part_1, rng);
	//std::cout << "PART 3 AFTER COAGULATION" << std::endl;
	//std::cout << part_3 << "\n\n";

	// Check the coagulated particle
	std::cout << "Checking coagulated particle..." << std::endl;
	if (part_3.NumberOfPrimaries() != 4
			|| ! within_tolerance(part_3.CollDiameter(), 3.02403, tolerance)
			|| ! within_tolerance(part_3.SurfaceArea(), 24.8059, tolerance)) {
		fail();
		return 1;
	}

	part_3.Sinter(5.0, cell, sm, rng, 1.0);
	//std::cout << "PART 3 AFTER COMPLETE SINTERING" << std::endl;
	//std::cout << part_3 << "\n\n";

	// Check the number of primaries (1) and the collision diameter
	std::cout << "Checking sintered particle..." << std::endl;
	if (part_3.NumberOfPrimaries() != 1
			|| ! within_tolerance(part_3.CollDiameter(), 2.2545, tolerance)
			|| ! within_tolerance(part_3.SurfaceArea(), 15.968, tolerance)) {
		fail();
		return 1;
	}
/*
	// Okay, time to check serialisation
	std::ofstream ofs("prilist.dat");
	part_1.Serialize(ofs);
	ofs.close();

	// Create a new particle
	std::ifstream ifs("prilist.dat");
	DefaultParticle part_4(0.0, pm);
	part_4.Deserialize(ifs, pm);
	ifs.close();*/

	std::cout << "All tests passed successfully.\n" << std::endl;

	return 0;
}
