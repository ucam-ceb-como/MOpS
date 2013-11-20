
#ifndef SWEEP_PRI_LIST_NODE
#define SWEEP_PRI_LIST_NODE

#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
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
    //! Return the value of component i
    double Component(unsigned int i) const;

    //! Return the volume of the node
    double Volume() const;

    //! Return the surface area
    double SurfaceArea() const;

    //! Return the surface area (static access)
    static double SurfaceArea(double diam);

    //! Return the mass
    double Mass() const;

    //! Return the diameter
    double Diameter() const;

    //! Return the diameter (static access)
    static double Diameter(double vol);

    // STATE SPACE ADJUSTMENT
    //! Set the particle model pointer
    void SetParticleModel(const Sweep::ParticleModel &model);

    //! Add dcomp to the state space
    unsigned int Adjust(
        const fvector &dcomp,
        unsigned int n);

    //! Merge this particle with rhs
    void Merge(PrimaryListNode &rhs);

protected:
    //! Particle model used to define the node
    const Sweep::ParticleModel* m_pmodel;

    //! The composition vector of the node
    fvector m_comp;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int /* version */) {
        ar & m_comp;
        // Prevent a clone of the original particle model from being serialized.
        //ar & m_pmodel;
    }

    //! Default constructor is meaningless
    PrimaryListNode(void);

};

} // AggModels
} // Sweep

#endif
