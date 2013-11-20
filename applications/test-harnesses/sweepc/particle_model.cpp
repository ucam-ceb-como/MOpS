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
#include "sprog.h"

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
    std::cout << "PRIMARY LIST UNIT TESTING" << std::endl;

    //######################################################################
    // Set up basic objects
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
    pm.SetBinTreeCoalThresh(1.0);   // Coalescence threshold
    pm.SetFractDim(1.8);            // Fractal dimension of aggregates

    // Typedef for easier use
    typedef AggModels::PrimaryList<AggModels::PrimaryListNode> DefaultParticle;
    typedef boost::mt19937 rng_type;

    // Create random number generator
    rng_type rng(1);

    // Set the relative tolerance for testing our doubles
    double tolerance(1.0e-5);

    //######################################################################
    // Create particles and check adjusting
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

    //######################################################################
    // Check coagulation
    part_1.Coagulate(part_2, rng);

    // Create a clone of part_1 again
    DefaultParticle part_3(part_1);

    // And let's go crazy and coagulate again!
    part_3.Coagulate(part_1, rng);

    // Check the coagulated particle
    std::cout << "Checking coagulated particle..." << std::endl;
    if (part_3.NumberOfPrimaries() != 4
            || ! within_tolerance(part_3.CollDiameter(), 3.00976, tolerance)
            || ! within_tolerance(part_3.SurfaceArea(), 24.6304, tolerance)) {
        fail();
        return 1;
    }

    //######################################################################
    // Okay, time to check serialisation
    std::ofstream ofs("prilist.binary");
    part_3.Serialize(ofs);
    ofs.close();

    // Create a new particle
    std::ifstream ifs("prilist.binary");
    DefaultParticle part_4(0.0, pm);
    part_4.Deserialize(ifs, pm);
    ifs.close();

    // Check the coagulated particle
    std::cout << "Checking serialisation/deserialisation..." << std::endl;
    if (part_4.NumberOfPrimaries() != 4
            || ! within_tolerance(part_3.CollDiameter(), 3.00976, tolerance)
            || ! within_tolerance(part_3.SurfaceArea(), 24.6304, tolerance)) {
        fail();
        return 1;
    }

    //######################################################################
    // Check sintering
    part_3.Sinter(5.0, cell, sm, rng, 1.0);

    // Check the number of primaries (1) and the collision diameter
    std::cout << "Checking sintered particle..." << std::endl;
    if (part_3.NumberOfPrimaries() != 1
            || ! within_tolerance(part_3.CollDiameter(), 2.2545, tolerance)
            || ! within_tolerance(part_3.SurfaceArea(), 15.968, tolerance)) {
        fail();
        return 1;
    }

    //######################################################################
    // Check combinations of sintering and adjusting
    part_2.Coagulate(part_2, rng);
    part_2.Adjust(fv, dvalues, rng, 1);
    part_2.Sinter(0.5, cell, sm, rng, 1.0);

    std::cout << "Checking partially sintered particles..." << std::endl;
    if (part_2.NumberOfPrimaries() != 2
            || ! within_tolerance(part_2.CollDiameter(), 1.83879, tolerance)
            || ! within_tolerance(part_2.SurfaceArea(), 10.4786, tolerance)) {
        fail();
        return 1;
    }

    // Check adjusting particles negatively and sintering
    part_3.Coagulate(part_3, rng);
    part_3.Sinter(0.3, cell, sm, rng, 1.0);
    fvector fv_neg(1, -4.0);
    part_3.Adjust(fv_neg, dvalues, rng, 1);
    part_3.Sinter(0.3, cell, sm, rng, 1.0);

    std::cout << "Checking negatively adjusted sintered particles..." << std::endl;
    if (part_3.NumberOfPrimaries() != 2
            || ! within_tolerance(part_3.CollDiameter(), 2.7274, tolerance)
            || ! within_tolerance(part_3.SurfaceArea(), 22.2906, tolerance)) {
        fail();
        return 1;
    }


    std::cout << "All tests passed successfully.\n" << std::endl;

    return 0;
}
