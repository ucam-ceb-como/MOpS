/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    The Imager class draws particle images and saves them to file.  Currently
    only POVRAY output files are generated, though it is expected that
    other formats will be added later.

    Particles are drawn as fractal aggregates using the primary-particle
    information.  Different particle models generate different primary-particle
    data structures, which are all catered for in this class.  The Draw
    functions automatically detect what form of data is available and draw
    the aggregate images accordingly.
*/

#ifndef SWEEP_IMAGER_H
#define SWEEP_IMAGER_H

#include "swp_params.h"
#include "swp_particle.h"
#include "swp_surfvol_primary.h"
#include "swp_pripart_primary.h"
#include "swp_imgnode.h"
#include <string>
#include <fstream>
#include <vector>

namespace Sweep
{
namespace Imaging
{
class ParticleImage
{
public:
    // Enumeration of possible collision regimes.  Different regimes
    // use different collision algorithms, so may produce differently
    // shaped aggregates.
    enum CollRegime{FreeMol=0};

    // Constructors.
    ParticleImage(void);               // Default constructor.
    ParticleImage(const Particle &sp); // Initialising constructor.

    // Destructors.
    ~ParticleImage(void);


    // COLLISION REGIME.

    // Returns the collision regime used to calculate the
    // particle image.
    CollRegime CollisionRegime(void) const;

    // Sets the collision regime used to calculate the
    // particle image.  If the collision regime is being 
    // changed then this function also forces a recalculation
    // of the image aggregate structure, assuming that it 
    // has already been initialised.
    void SetCollRegime(CollRegime creg);


    // PARTICLE IMAGE DATA CONSTRUCTION.

    // Constructs the particle image from the given particle.
    void Construct(const Particle &sp);


    // RENDERING FUNCTIONS.

    // Draws the particle image to a POVRAY file.
    void WritePOVRAY(
        std::ofstream &file // Output file stream.
        );

private:
    // The collision regime used to calculate the image
    // aggregate structure.  The default is free-molecular.
    CollRegime m_creg;

    // The image aggregate structure root node.
    ImgNode m_root;

    // Amount of necking between particles in output.
    static const real m_necking;

    // AGGREGATE SPHERE-TREE CONSTRUCTORS (FREE-MOLECULAR).

    // Constructs a PNode sphere-tree aggregate from the given 
    // particle using free-molecular collision dynamics.
    void constructAgg_FM(const Particle &sp);

    // Constructs a PNode sphere-tree aggregate from the given 
    // pri-part list primary using free-molecular collision dynamics.
    void constructAgg_FM(const AggModels::PriPartPrimary &pri);

    // Constructs a PNode sphere-tree aggregate with uniform 
    // primaries (equal diameter).  The diameter and primary
    // count are passed as arguments.
    void uniformAgg_FM(
        unsigned int n, // Primary particle count.
        real d          // Primary particle diameter.
        );

    // Calculates the aggregate structure down from the given
    // node.  Assumes that the tree leaves have been initialised
    // with the correct radii, and recalculates their positions.
    static void calc_FM(ImgNode &node);

    // Calculates the z-displacement of a bullet sphere for a +ve
    // collision with a target sphere.  Returns true if the
    // spheres collide, otherwise false.
    static bool calcCollZ(
        real sumrsqr,           // Sum of the target and bullet radii squared.
        real dxsqr, real dysqr, // The X and Y plane displacements of the bullet squared.
        real &dz                // The output z-axis displacement of the bullet (+ve).
        );

    // OUTPUT FUNCTIONS.

};
};
};

#endif
