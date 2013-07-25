#pragma once
#include "swp_primary.h"
#include "swp_pri_list_cache.h"
#include <boost/random/uniform_int_distribution.hpp>

namespace Sweep {

namespace AggModels {

template <class NodeT, class CacheT>
class PrimaryList: public Primary
{
protected:

	// BASIC DEFINITIONS
    //! Default constructor is private (meaningless)
    PrimaryList(void);

    //! Typedef the primary particle list
    typedef std::vector<NodeT*> PrimaryPtrVector;

    //! The list of primary particles
    PrimaryPtrVector mPrimaries;

    // A pointer to the particle data cache being used
    CacheT* mCache;

    // BASIC STD-LIBRARY-LIKE MANIPLULATION OF THE CLASS
    //! Return a pointer to the first primary
    typename PrimaryPtrVector::iterator begin()
    	{return mPrimaries.begin();}
    typename PrimaryPtrVector::const_iterator begin() const
    	{return mPrimaries.begin();}

    //! Return a pointer to the last primary
    typename PrimaryPtrVector::iterator end()
    	{return mPrimaries.end();}
    typename PrimaryPtrVector::const_iterator end() const
    	{return mPrimaries.end();}

    //! Return the size of the primaries vector
    size_t size() const {return mPrimaries.size();}

    //! Clear the PrimaryPtrVector and memory associated with them
    void ClearPrimaries() {
    	for (size_t i(0); i != size(); i++) delete mPrimaries[i];
    }

    // PHYSICAL AND CHEMICAL MANIPULATION OF THE CLASS
    //! Select a node in the list
    NodeT* SelectNode(rng_type &rng) const {
    	// TODO: Expand in the future to weight selection
    	// e.g., Select bigger nodes instead of smaller ones
    	boost::random::uniform_int_distribution<unsigned> dist(0, size());
    	return mPrimaries[dist(rng)];
    }

public:
    // BASIC CREATION AND DELETION METHODS
    // The cache must be a friend of the list class
    friend class PrimaryListCache;

	//! Construct a primary at time
	PrimaryList(const double time, const Sweep::ParticleModel &model):
		Primary(time, model),
		mPrimaries(),
		mCache(NULL) {}

    //! Copy constructor
    PrimaryList(const PrimaryList &copy):
    	mPrimaries(copy.mPrimaries),
    	mCache(copy.mCache) {
    	//TODO: fill in copy constructor
    }

    //! Stream-reading constructor
    PrimaryList(std::istream &in, const Sweep::ParticleModel &model):
		Primary(in, model),
		mPrimaries(),
		mCache(NULL) {}

    //! Constructor including position
    //PrimaryList(const double time, const double position,
    //    const Sweep::ParticleModel &model);

	//! Destructor
	virtual ~PrimaryList(void) {
		// Delete memory associated with primaries and cache
		ClearPrimaries();
		delete mCache;
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
        const Sweep::AggModels::PrimaryList<NodeT, CacheT> &pri) {
            os << "[Primary List], n=" << pri.size() << "\n";
            for (typename PrimaryPtrVector::const_iterator i = pri.begin(); i != pri.end(); ++i) {
            	os << " " << *(*i);
            }
            return os;
    }


    // SWEEP MECHANISM METHODS
    //! Return the aggregation model index
    AggModelType AggID() const {return Sweep::AggModels::PrimaryList_ID;};

    //! Coagulation of this particle with rhs
    PrimaryList &Coagulate(const Primary &rhs, rng_type &rng) {
    	std::cout << "Coagulating!" << std::endl;
    	return *this;
    }

    //! Updates the particle cache using the particle details
    void UpdateCache() {
    	std::cout << "Updating ze cache!" << std::endl;
    }

    //! Adjusts a particle according to a surface reaction
    unsigned int Adjust(
        const fvector &dcomp,
        const fvector &dvalues,
        rng_type &rng,
        unsigned int n) {
    	std::cout << "Adjusting!" << std::endl;
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

     //! Sinters a particle for time dt
    void Sinter(
        double dt,
        Cell &sys, 
        const Processes::SinteringModel &model,
        rng_type &rng,
        double wt) {
    	std::cout << "Sintering!" << std::endl;
    }

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const {
    	std::cout << "Serialising!" << std::endl;
    }

    // Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,
        const Sweep::ParticleModel &model
        ) {
    	std::cout << "Deserialising!" << std::endl;
    }

	//! Gets the number of active sites (always component 0)
	double GetSites() const { return m_comp[0]; }

	// Return the sintering rate
	double GetSintRate() const {return 0.0;}

	//! Return the fractional coverage
	double GetCoverageFraction() const {return 0.0;}




	// NOT SURE IF THESE SHOULD BE PUBLIC.
	void AddNode(NodeT pri) {
		// Create a clone of the node, and add it to the list
		NodeT* new_pri = new NodeT(pri);
		mPrimaries.push_back(new_pri);
	}

};

} // AggModels

} // Sweep
