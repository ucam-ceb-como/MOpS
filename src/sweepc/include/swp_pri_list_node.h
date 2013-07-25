// This class has the actual physical functionality of the tree.
// That is, a Primary object added to the tree is converted into one of these nodes
// The idea is that if a new particle model is to be added, this can be subclassed and inserted into the 
// binary tree template without having to worry about the caching and connectivity rubbish
// all other properties stored in the cache?

#ifndef SWEEP_PRI_LIST_NODE
#define SWEEP_PRI_LIST_NODE

#include "swp_particle_model.h"
#include "swp_primary.h"

namespace Sweep {

namespace AggModels {

class PrimaryListNode
{
public:
    //! Create empty primary with just particle model
    PrimaryListNode(const Sweep::ParticleModel &model);

    //! Create with a composition vector
    PrimaryListNode(
        const Sweep::ParticleModel &model,
        const fvector &comp);

    //! Copy constructor
    PrimaryListNode(const PrimaryListNode &copy);

    //! Destructor
    virtual ~PrimaryListNode(void);

    //! Overload of the = operator
    virtual PrimaryListNode &operator=(const PrimaryListNode &rhs);

    //! Overload of the << operator
    friend std::ostream& operator<<(
        std::ostream &os,
        const Sweep::AggModels::PrimaryListNode &pri);

    // DATA ACCESS METHODS
    //! Return the volume of the node
    double Volume() const;

    //! Return the surface area
    double SurfaceArea() const;

    //! Return the mass
    double Mass() const;

    //! Return the diameter
    double Diameter() const;

    // STATE SPACE ADJUSTMENT
    //! Add dcomp to the state space
    unsigned int Adjust(
        const fvector &dcomp,
        unsigned int n);

protected:
    //! Particle model used to define the node
    const Sweep::ParticleModel* mPModel;

    //! The composition vector of the node
    fvector mComp;

private:
    //! Default constructor is meaningless
    PrimaryListNode(void);
 
};

} // AggModels
} // Sweep

#endif
