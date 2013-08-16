#ifndef SWP_PRI_LIST_CACHE_HPP
#define SWP_PRI_LIST_CACHE_HPP

#include "swp_pri_list.hpp"
#include "swp_params.h"

namespace Sweep {

namespace AggModels {

// Forward declare the list class
template <class NodeT> class PrimaryList;

template <class NodeT>
class PrimaryListCache
{
public:
    //! A typedef for functions used for data collection
    typedef std::vector<double (NodeT::*)()const> DataFunctionPtrVector;

    //! Default constructor
    PrimaryListCache(void) {};

    //! Copy constructor
    PrimaryListCache(const PrimaryListCache &copy) {
        // Nothing to copy in this particular cache
    }

    //! Default destructor
    ~PrimaryListCache(void) {};

    /*!
     * Given a particle and a node's data function, return the sum across
     * all nodes of that function.
     *
     * @param pri           Particle object
     * @param fn_pointer    Function to calculate sum with
     * @return              Sum across all nodes, using the supplied function
     */
    double GetSum(
            PrimaryList<NodeT> &pri,
            double (NodeT::*fn_pointer)()const) const {

        double sum(0.0);

        // Use the node's function pointer to get the desired value
        for (typename PrimaryList<NodeT>::PrimaryPtrVector::const_iterator i  = pri.begin();
                i != pri.end(); ++i) {
            sum += ((*i)->*fn_pointer)();
        }

        return sum;
    }

    /*!
     * Given a particle and a vector of functions, return a vector containing
     * the sums across all nodes for each function
     *
     * @param pri           Particle object
     * @param fn_ptrs       Vector of function pointers
     * @return              Vector of sums, using the functions supplied
     */
    fvector GetSums(
            PrimaryList<NodeT> &pri,
            DataFunctionPtrVector fn_ptrs) const {

        fvector sums(fn_ptrs.size(), 0.0);

        // Loop over the primaries in the list
        for (typename PrimaryList<NodeT>::PrimaryPtrVector::const_iterator i  = pri.begin();
                i != pri.end(); ++i) {
            // Now loop over the data collection function pointers
            for (typename DataFunctionPtrVector::const_iterator f = fn_ptrs.begin();
                    f != fn_ptrs.end(); ++f) {
                sums += ((*i)->*(*f))();
            }
        }

        return sums;
    }

    /*!
     * Directly update the properties of a PrimaryList particle. This function
     * is written to minimise double-looping over the PrimaryPtrVector and
     * recalculation of properties in the Node class.
     *
     * @param pri         Primary to update the cache for
     */
    void FastUpdate(PrimaryList<NodeT> &pri) {

        // Initialise some working values
        double vol(0.0);
        pri.m_mass = 0.0;
        pri.m_vol = 0.0;
        pri.m_avg_dpri = 0.0;


        // Loop over primaries - get the mass and the volume
        for (typename PrimaryList<NodeT>::PrimaryPtrVector::const_iterator i  = pri.begin();
                i != pri.end(); ++i) {
            // Loop over the components and sum them up
            for (size_t j(0); j != pri.m_comp.size(); ++j) {
                double c = (*i)->Component(j);
                pri.m_comp[j] += c;
                pri.m_mass += c * pri.m_pmodel->Components(j)->MolWt() / Sweep::NA;
            }

            // Get volume
            vol = (*i)->Volume();
            pri.m_vol += vol;
            pri.m_avg_dpri += (*i)->Diameter(vol);
        }
        pri.m_avg_dpri /= (double) pri.size();
        pri.m_diam = NodeT::Diameter(pri.m_vol);

        // Get the surface area
        // Note that spherical surface is handled by Primary::SphSurfaceArea
        pri.m_surf = pri.CalcSurfaceArea();

        // Now calculate collision and mobility diameter
        pri.m_dcol = pri.CalcCollisionDiameter();
        pri.m_dmob = pri.CalcMobilityDiameter();
    }

};

} // AggModels
} // Sweep

#endif
