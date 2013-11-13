#ifndef SWP_PRI_LIST_CONNECTOR_H_
#define SWP_PRI_LIST_CONNECTOR_H_

#include "swp_cell.h"
#include "swp_params.h"
#include <boost/random/poisson_distribution.hpp>
#include <algorithm>

namespace Sweep {
namespace AggModels {

// Forward-declare the list class
template <class NodeT> class PrimaryList;

template <class NodeT>
class PrimaryListConnector {
protected:
    //! A pointer to the left node
    NodeT* m_left;

    //! A pointer to the right node
    NodeT* m_right;

    //! The common surface area shared between nodes
    double m_common_surface;

    //! The sintering level (cached)
    double m_sint_level;

    //! Flag indicating whether the left and right nodes should be merged, and
    //! this object deleted.
    bool m_merge_me;

    //! Calculate the sintering level
    double CalcSinteringLevel() const {
        double sl(0.0);

        if (m_common_surface > 0.0) {
            sl = (SphSurfaceArea() / m_common_surface - Sweep::TWO_ONE_THIRD) /
                    (1.0 - Sweep::TWO_ONE_THIRD);
        }
        // Ensure the parameter is between 0 and 1.
        if (sl > 1.0) sl = 1.0;
        if (sl < 0.0) sl = 0.0;
        return sl;
    }

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) const {
        ar & m_left;
        ar & m_right;
        ar & m_common_surface;
        ar & m_sint_level;
        ar & m_merge_me;
    }

public:
    // The PrimaryList must be a friend class
    friend class PrimaryList<NodeT>;


    //! Default constructor
    PrimaryListConnector():
        m_left(NULL),
        m_right(NULL),
        m_common_surface(0.0),
        m_sint_level(1.0),
        m_merge_me(false) {}

    //! Copy constructor
    PrimaryListConnector(const PrimaryListConnector &copy):
        m_left(copy.m_left),
        m_right(copy.m_right),
        m_common_surface(copy.m_common_surface),
        m_sint_level(copy.m_sint_level),
        m_merge_me(copy.m_merge_me)  {}

    //! Main constructor to connect two nodes
    PrimaryListConnector(NodeT *lhs, NodeT *rhs):
        m_left(NULL),
        m_right(NULL),
        m_common_surface(0.0),
        m_sint_level(1.0),
        m_merge_me(false) {

        // Set the pointers
        m_left = lhs;
        m_right = rhs;

        // Calculate the common surface as the sum of the spherical surface
        // areas of the nodes
        m_common_surface = lhs->SurfaceArea() + rhs->SurfaceArea();

        // Update other properties
        Update();
    }

    //! Overload of the << operator
    friend std::ostream& operator<<(
        std::ostream &os,
        const Sweep::AggModels::PrimaryListConnector<NodeT> &conn) {
            os << "[Primary List Connector],";
            os <<  " csurf=" << conn.m_common_surface;
            os <<  ", add=" << static_cast<void const *>(&conn) << std::endl;
            os << "  LHS: " << *(conn.m_left);
            os << "  RHS: " << *(conn.m_right);

            return os;
    }


    //! Default destructor
    virtual ~PrimaryListConnector() {
        // Nothing happens here, as the PrimaryList class should handle
        // the memory associated with nodes
        m_left = NULL;
        m_right = NULL;
    }

    //! Return the volume
    double Volume() const {
        return m_left->Volume() + m_right->Volume();
    }

    //! Return the common surface area
    double SurfaceArea() const {return m_common_surface;}

    //! Return the spherical surface area of the two nodes
    double SphSurfaceArea() const {
        // This uses the basic ListNode's static declarations of surface area and
        // diameter, but this could be changed if this is not possible to do
        // when more complex node classes are added
        return NodeT::SurfaceArea(NodeT::Diameter(Volume()));
    }

    //! Return the sintering level
    double SinteringLevel() const {return m_sint_level;}

    //! Should this connector be merged?
    bool Merge() const {return m_merge_me;}

    //! Update the cached properties of the connector
    void Update() {
        // Just need to update the sintering level at this stage
        m_sint_level = CalcSinteringLevel();
    }

    // Sinter the connection
    void Sinter(
            double dt,
            double t,
            Cell &sys,
            const Processes::SinteringModel &model,
            rng_type &rng) {

        const double surf_sph = SphSurfaceArea();

        // Declare time step variables.
        double t1=0.0, delt=0.0, tstop=dt;
        double r=0.0;

        // Define the maximum allowed change in surface
        // area in one internal time step (10% spherical surface).
        double dAmax = 0.1 * surf_sph;

        // The scale parameter discretises the delta-S when using
        // the Poisson distribution.  This allows a smoother change
        // (smaller scale = higher precision).
        double scale = 0.01;

        // Perform integration loop.
        while (t1 < tstop)
        {
            // Calculate sintering rate.
            r = model.Rate(t+t1, sys.GasPhase().Temperature(), *this);

            if (r > 0) {
                // Calculate next time-step end point so that the
                // surface area changes by no more than dAmax.
                delt = dAmax / std::max(r, 1.0e-300);

                // Approximate sintering by a poisson process.  Calculate
                // number of poisson events.
                double mean;

                if (tstop > (t1+delt)) {
                    // A sub-step, we have changed surface by dAmax, on average
                    mean = 1.0 / scale;
                } else {
                    // Step until end.  Calculate degree of sintering explicitly.
                    mean = r * (tstop - t1) / (scale*dAmax);
                }
                boost::random::poisson_distribution<unsigned, double> repeatDistribution(mean);
                const unsigned n = repeatDistribution(rng);

                // Adjust the surface area.
                if (n > 0) {
                    m_common_surface -= (double)n * scale * dAmax;

                    // Check that primary is not completely sintered.
                    if (m_common_surface <= surf_sph) {
                        m_common_surface = surf_sph;
                        // Activate the flag indicating that the nodes should be merged
                        m_merge_me = true;
                        break;
                    }
                }

                // Set t1 for next time step.
                t1 += delt;
            }

        }

        Update();
    }



};

} /* namespace AggModels */
} /* namespace Sweep */
#endif /* SWP_PRI_LIST_CONNECTOR_H_ */
