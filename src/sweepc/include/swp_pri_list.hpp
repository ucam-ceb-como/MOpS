#ifndef SWP_PRI_LIST_HPP
#define SWP_PRI_LIST_HPP

#include <boost/random/uniform_int_distribution.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "swp_primary.h"
#include "swp_pri_list_cache.hpp"
#include "swp_pri_list_connector.hpp"


namespace Sweep {

namespace AggModels {

template <class NodeT>
class PrimaryList: public Primary
{
protected:

    // BASIC DEFINITIONS
    //! Default constructor is private (meaningless)
    PrimaryList(void);

    //! Typedef the primary particle list
    typedef std::vector<NodeT*> PrimaryPtrVector;

    //! Typedef the primary connector list
    typedef std::vector<PrimaryListConnector<NodeT>* > ConnectorPtrVector;

    //! The list of primary particles
    PrimaryPtrVector m_primaries;

    //! The list of connectors of primaries
    ConnectorPtrVector m_connectors;

    //! A pointer to the particle data cache being used
    PrimaryListCache<NodeT>* m_cache;

    //! The average primary diameter of the particle
    double m_avg_dpri;

    // BASIC STD-LIBRARY-LIKE MANIPLULATION OF THE CLASS
    //! Return a pointer to the first primary
    typename PrimaryPtrVector::iterator begin()
        {return m_primaries.begin();}
    typename PrimaryPtrVector::const_iterator begin() const
        {return m_primaries.begin();}

    //! Return a pointer to the last primary
    typename PrimaryPtrVector::iterator end()
        {return m_primaries.end();}
    typename PrimaryPtrVector::const_iterator end() const
        {return m_primaries.end();}

    //! Return the size of the primaries vector
    size_t size() const {return m_primaries.size();}

    //! Clear the PrimaryPtrVector and memory associated with them
    void ClearPrimaries() {
        for (size_t i(0); i != size(); i++) delete m_primaries[i];
        for (size_t i(0); i != m_connectors.size(); i++) delete m_connectors[i];
    }

    /*!
     * Given a pointer to a primary, return the numerical index of the node
     *
     * @param pri     Pointer to node
     * @return        Index in the primaries vector
     */
    unsigned int PrimaryIndex(NodeT* pri) const {
        unsigned int index(-1);
        for (size_t i(0); i != size(); ++i) {
            if (pri == m_primaries[i]) index = i;
        }
        if (index < 0 || index >= size())
            throw std::runtime_error("Error locating index of node, in"
                    "Sweep::AggModels::PrimaryList::PrimaryIndex.");
        return index;
    }

    // PHYSICAL AND CHEMICAL MANIPULATION OF THE CLASS
    //! Select a node in the list
    NodeT* SelectNode(rng_type &rng) const {
        // TODO: Expand in the future to weight selection
        // e.g., Select bigger nodes instead of smaller ones
        boost::random::uniform_int_distribution<unsigned> dist(0, size()-1);
        return m_primaries[dist(rng)];
    }

    //! Clear the cached properties of the particle
    void ResetCachedProperties() {
        m_mass = 0.0;
        m_diam = 0.0;
        m_dcol = 0.0;
        m_dmob = 0.0;
        m_surf = 0.0;
        m_vol  = 0.0;
    }

    //! Calculates the collision diameter from cached particle properties
    virtual double CalcCollisionDiameter() const {
        return (6.0 * m_vol / m_surf) * pow(
                    pow(m_surf, 3.0) / (36.0 * Sweep::PI * m_vol * m_vol),
                    (1.0/m_pmodel->GetFractDim()
                    )
                );
    }

    //! Calculates the mobility diameter from cached particle properties
    virtual double CalcMobilityDiameter() const {
        double d_mob(m_diam);
        double n_pri = (double) size();
        // Presently hard-coded as T, P are abstracted from this
        // function's view

        if (false) {
            // SF regime mobility diameter
            d_mob *= 0.9 * m_avg_dpri / n_pri;
            d_mob *= sqrt(m_pmodel->GetFractDim() / (m_pmodel->GetFractDim() + 2.0));
            d_mob *= pow(n_pri, (1.0/m_pmodel->GetFractDim()));
        } else {
            // FM regime mobility diameter
            d_mob *= m_avg_dpri / (double)n_pri;
            d_mob *= sqrt(0.802 * (n_pri - 1.0) + 1.0);
        }

        if (d_mob < m_diam) d_mob = m_diam;
        return d_mob;
    }

    //! Calculates the surface area from cached particle properties
    virtual double CalcSurfaceArea() const {
        // If there's just one primary, it's easy..
        if (size() == (size_t)1) return (m_primaries[0])->SurfaceArea();

        double sl(0.0);
        double n_1_3 = pow((double)size(), - Sweep::ONE_THIRD);

        // Loop over the connectors to get the sum of the sintering level
        for (typename PrimaryList<NodeT>::ConnectorPtrVector::const_iterator i = m_connectors.begin();
                i != m_connectors.end(); ++i) {
            sl += (*i)->SinteringLevel();
        }
        sl /= (double) size();

        return SphSurfaceArea() / (sl * (1.0 - n_1_3) + n_1_3);
    }

    /*!
     * Create a copy of this particle's structure, putting the object pointers
     * into the supplied vectors
     *
     * @param cloned_primaries
     * @param cloned_connectors
     */
    void CloneStructure(
            PrimaryPtrVector &cloned_primaries,
            ConnectorPtrVector &cloned_connectors) const {

        // Ensure these vectors are empty first
        cloned_primaries.clear();
        cloned_connectors.clear();

        // Loop over the primaries, creating a copy of each one
        for (size_t i(0); i != size(); ++i) {
            NodeT* pri = new NodeT(*(m_primaries[i]));
            cloned_primaries.push_back(pri);
        }

        // Now loop over the *this* connector vector, new connectors
        std::vector<std::pair<unsigned, unsigned> > index_pairs;
        for (typename PrimaryList<NodeT>::ConnectorPtrVector::const_iterator it = m_connectors.begin();
                it != m_connectors.end(); ++it) {
            // Copy the connector
            PrimaryListConnector<NodeT>* conn = new PrimaryListConnector<NodeT>(*(*it));

            // Get the pointers right
            conn->m_left = cloned_primaries[PrimaryIndex((*it)->m_left)];
            conn->m_right = cloned_primaries[PrimaryIndex((*it)->m_right)];

            cloned_connectors.push_back(conn);
        }
    }


    /*!
     * Delete an object and its entry in the list
     *
     * @param list  Vector to remove object from
     * @param id    ID of object to remove
     */
    template <typename T>
    void DeleteFrom(T &list, unsigned int id) {
        // Delete the node from memory
        delete list[id];
        // And erase its pointer in the list
        list.erase(list.begin() + id);
    }

    /*!
     * Merge the connector object, passed by reference
     *
     * @param conn   Connector object to merge
     */
    void MergeConnection(PrimaryListConnector<NodeT> &conn) {
        assert(conn.Merge());

        // First work out the id of the RHS
        NodeT* rhs = conn.m_right;
        const unsigned int rhs_id = PrimaryIndex(rhs);

        // Now loop over the connectors, finding which have pointers to the RHS
        for (typename PrimaryList<NodeT>::ConnectorPtrVector::iterator it = m_connectors.begin();
                        it != m_connectors.end(); ++it) {
            // Only an issue for connectors other than the one being merged
            if ((*it) != &conn) {
                // Check if our RHS pointer is the same as one of those in
                // another connector
                if (rhs == (*it)->m_left) {
                    // If it is, set the *it connector's LHS pointer to this
                    // connector's LHS
                    (*it)->m_left = conn.m_left;
                } else if (rhs ==(*it)->m_right) {
                    // Same as above, but for RHS
                    (*it)->m_right = conn.m_left;
                }
            }
        }

        // Add the state space of the RHS to the LHS particle
        conn.m_left->Merge(*(conn.m_right));

        // Now delete the RHS particle
        DeleteFrom(m_primaries, rhs_id);
    }


    /*!
     * Use Boost's inbuilt serialization to write out this object. We serialize
     * the following data (in this order):
     *    1. number of nodes
     *    2. each node
     *    3. the connections between the nodes
     *
     * @param ar   A boost archive object
     */
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) const {
        ar & m_primaries;
        ar & m_connectors;
    }


public:
    // BASIC CREATION AND DELETION METHODS
    // The cache must be a friend of the list class
    friend class PrimaryListCache<NodeT>;
    friend class boost::serialization::access;

    //! Construct a primary at time
    PrimaryList(const double time, const Sweep::ParticleModel &model):
        Primary(time, model),
        m_primaries(),
        m_connectors(),
        m_cache(NULL),
        m_avg_dpri(0.0) {

        // Create a new cache for use with this particle.
        m_cache = new PrimaryListCache<NodeT>;
    }

    //! Copy constructor
    PrimaryList(const PrimaryList &copy):
        Primary(copy),
        m_primaries(),
        m_connectors(),
        m_cache(NULL),
        m_avg_dpri(copy.m_avg_dpri) {

        // Clone the cache
        m_cache = new PrimaryListCache<NodeT>(*(copy.m_cache));

        // Clone the particle structure from the copy into *this*'s primaries
        // and connectors vectors
        copy.CloneStructure(m_primaries, m_connectors);
    }

    //! Stream-reading constructor
    PrimaryList(std::istream &in, const Sweep::ParticleModel &model):
        Primary(in, model),
        m_primaries(),
        m_connectors(),
        m_cache(NULL),
        m_avg_dpri (0.0) {}

    //! Constructor including position
    //PrimaryList(const double time, const double position,
    //    const Sweep::ParticleModel &model);

    //! Destructor
    virtual ~PrimaryList(void) {
        // Delete memory associated with primaries and cache
        ClearPrimaries();
        delete m_cache;
    }

    //! Equals operator
    virtual PrimaryList &operator=(const Primary &rhs) {
        operator=(dynamic_cast<const PrimaryList&>(rhs));
        return *this;
    }

    //! Get a clone of the particle
    virtual PrimaryList *const Clone(void) const {
        return new PrimaryList(*this);
    }

    //! Overload of the << operator
    friend std::ostream& operator<<(
        std::ostream &os,
        const Sweep::AggModels::PrimaryList<NodeT> &pri) {
            os << "[Primary List], n=" << pri.size();
            os << ", vol=" << pri.m_vol << ", surf=" << pri.m_surf;
            os << ", dcol=" << pri.m_dcol << "\n";
            os << " List of primary nodes\n";
            for (typename PrimaryList<NodeT>::PrimaryPtrVector::const_iterator it = pri.m_primaries.begin();
                it != pri.m_primaries.end(); ++it) {
                os << "  " << *(*it);
            }
            for (typename PrimaryList<NodeT>::ConnectorPtrVector::const_iterator it = pri.m_connectors.begin();
                it != pri.m_connectors.end(); ++it) {
                os << " " << *(*it);
            }
            return os;
    }


    // SWEEP MECHANISM METHODS
    //! Return the aggregation model index
    AggModelType AggID() const {return Sweep::AggModels::PrimaryList_ID;};

    //! Coagulation of this particle with rhs
    PrimaryList &Coagulate(const Primary &rhs, rng_type &rng) {

        // Cast the RHS as a PrimaryList
        const PrimaryList<NodeT>* rhs_pl
            = dynamic_cast<const AggModels::PrimaryList<NodeT>* >(&rhs);

        // First we clone the rhs particle's structure
        PrimaryPtrVector new_primaries(0);
        ConnectorPtrVector new_connectors(0);
        rhs_pl->CloneStructure(new_primaries, new_connectors);

        // Now get a random number to select a primary index of the rhs particle
        boost::random::uniform_int_distribution<unsigned> dist(0, new_primaries.size()-1);

        // Select a particle to join the lists with, and create a connector
        PrimaryListConnector<NodeT>* conn
            = new PrimaryListConnector<NodeT>(SelectNode(rng), new_primaries[dist(rng)]);
        m_connectors.push_back(conn);

        // Finally, we add the lists
        m_primaries.insert(m_primaries.end(), new_primaries.begin(), new_primaries.end());
        m_connectors.insert(m_connectors.begin(), new_connectors.begin(), new_connectors.end());

        UpdateCache();
        return *this;
    }

    //! Updates the particle cache using the particle details
    void UpdateCache() {
        ResetCachedProperties();

        // Get sums of various quantities from the cache
        //m_vol = mCache->GetSum(*this, &NodeT::Volume);
        m_cache->FastUpdate(*this);
    }

    /*!
     * Sets the composition vector of the particle. This is to be used
     * with EXTREME caution. Only really for initialisation of a particle.
     *
     * @param comp    Composition vector of node to set
     */
    void SetComposition(const fvector &comp) {
        // Create a new primary node and insert it
        NodeT* pri = new NodeT(*m_pmodel, comp);
        m_primaries.push_back(pri);
    }

    /*!
     * Overload of the Primary Adjust function. Used for performing a surface
     * reaction process on a particle.
     *
     * @param dcomp     Composition change vector
     * @param dvalues   Tracker change vector
     * @param rng       Random number generator
     * @param n         Number of times to adjust
     * @return          Number of times actually adjusted
     */
    unsigned int Adjust(
        const fvector &dcomp,
        const fvector &dvalues,
        rng_type &rng,
        unsigned int n) {

        // Select a primary randomly and adjust it
        NodeT* pri = SelectNode(rng);
        n = pri->Adjust(dcomp, n);

        // Add the tracker values (only implemented for the 'master' particle)
        for (size_t i(0); i != std::min(m_values.size(), dvalues.size()); ++i) {
            m_values[i] += dvalues[i] * (double)n;
        }

        // Update the cache
        UpdateCache();

        return n;

    }

    //! Adjusts a particle according to an interparticle reaction
    unsigned int AdjustIntPar(
        const fvector &dcomp,
        const fvector &dvalues,
        rng_type &rng,
        unsigned int n) {
        std::cout << "Adjusting for Interparticle!" << std::endl;
        return n;
    }

    /*!
     * Overload of the Primary::Sinter function. Sinters the list structure.
     * The sintering is first done on the connector level. Those which need
     * to be merged are flagged at this stage. Then we loop over the connectors,
     * merging those which need to be (this adds the nodes and maintains the
     * pointer structure). The leftover connector is then deleted.
     *
     * @param dt     Time to sinter over
     * @param sys    Cell to do it in
     * @param model  Sintering model to use
     * @param rng    Random number generator
     * @param wt     Statistical weight of particle
     */
    void Sinter(
        double dt,
        Cell &sys,
        const Processes::SinteringModel &model,
        rng_type &rng,
        double wt) {

        // Loop over the connectors, calculating the sintering of each
        // connection between nodes
        for (typename PrimaryList<NodeT>::ConnectorPtrVector::iterator it = m_connectors.begin();
                        it != m_connectors.end(); ++it) {
            (*it)->Sinter(dt, m_time, sys, model, rng);
        }

        // Now loop over the connectors again, merging those with their merge
        // activated
        std::vector<unsigned int> ids_to_merge;
        for (size_t i(0); i != m_connectors.size(); ++i) {
            if (m_connectors[i]->Merge()) {
                // MergeConnection adds the state space of the nodes and keeps
                // the pointer structure intact
                MergeConnection(*(m_connectors[i]));
                // Also need to flag this connector for deletion
                ids_to_merge.push_back(i);
            }
        }

        // Delete any now-redundant connectors (reversed to avoid segfault)
        for (size_t i = ids_to_merge.size(); i-- > 0; ) {
            DeleteFrom(m_connectors, i);
        }

        // And finally, update the particle cache
        UpdateCache();

    }


    /*!
     * Writes the object to a binary stream. Using boost_serialize.
     * Oh yeah!
     *
     * @param out
     */
    void Serialize(std::ostream &out) const {
        boost::archive::binary_oarchive oa(out);
        serialize(oa, 0);
    }

    /*!
     * Reads the object from a stream.
     *
     * @param in      Input stream
     * @param model   Particle model
     */
    void Deserialize(
        std::istream &in,
        const Sweep::ParticleModel &model
        ) {

        // Get boost to do the hard work for us
        boost::archive::binary_iarchive ia(in);
        //serialize(ia, 0);

        // Now need to make sure the particle model pointers are correct,
        // and update the particle
        m_pmodel = &model;
        UpdateCache();
    }

    //! Gets the number of active sites (always component 0)
    double GetSites() const { return m_comp[0]; }

    //! Return the sintering rate
    double GetSintRate() const {return 0.0;}

    //! Return the fractional coverage
    double GetCoverageFraction() const {return 0.0;}

    //! Return the number of primary particles
    unsigned int NumberOfPrimaries() const {return (unsigned int) size();}




};

} // AggModels

} // Sweep

#endif
