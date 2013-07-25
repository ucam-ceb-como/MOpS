//#include <boost/program_options.hpp>
//namespace po = boost::program_options;

#include <iostream>
#include <string>
#include <vector>

// Sweep includes for this file
#include "swp_particle_model.h"
#include "swp_bintree_primary.h"
#include "swp_pri_list.hpp"
#include "swp_pri_list_node.h"
#include "swp_params.h"

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
	Component* c = new Component(1.0, 1.0, 1.0, "dummy");

	// Now create the particle model
	ParticleModel pm;
	pm.AddComponent(*c); // Takes ownership and clears mem

	// Set basic properties of the model, used by particle models
	pm.SetBinTreeCoalThresh(1.0);	// Coalescence threshold
	pm.SetFractDim(1.8);			// Fractal dimension of aggregates

	// Typedef for easier use
	typedef AggModels::PrimaryList<
			AggModels::PrimaryListNode,
			AggModels::PrimaryListCache> DefaultParticle;

	/*
	 * ###################################################################
	 * Main testing of particle models
	 * 	prim = A primary particle (e.g. ListNode class)
	 * 	part = An actual computational particle
	 */
	// Create a basic node primary to insert into a particle
	fvector fv(1, 1.0);
	AggModels::PrimaryListNode prim_1(pm, fv);

	// Create some computational particles
	DefaultParticle part_1(0.0, pm);
	part_1.AddNode(prim_1); // Clones prim_1 and takes ownership
	std::cout << part_1;

	return 0;
}
