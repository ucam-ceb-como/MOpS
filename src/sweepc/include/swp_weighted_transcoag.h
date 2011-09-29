#ifndef SWP_WEIGHTED_TRANSCOAG_H
#define	SWP_WEIGHTED_TRANSCOAG_H

#include "swp_coagulation.h"
#include "swp_process_type.h"

namespace Sweep
{
// Forward declare Mechanism class.
class Mechanism;

// Forward declare class used for sums in the binary tree
class TreeTransCoagWeightedCache;

namespace Transport
{
    // Forward declaration of unused argument type
    struct TransportOutflow;
}

namespace Processes
{

class WeightedTransitionCoagulation : public Coagulation
{
public:
    //! Default constructors.
    WeightedTransitionCoagulation(const Sweep::Mechanism &mech, const CoagWeightRule weight_rule);

    //! Deserialisation
    WeightedTransitionCoagulation(
        std::istream &in,            //  - Input stream.
        const Sweep::Mechanism &mech //  - Parent mechanism.
        );

    //! Virtual destructor
    virtual ~WeightedTransitionCoagulation() {};

    //! Clone object
    virtual WeightedTransitionCoagulation* const Clone() const;

    //! Returns the process type for identification during serialisation
    virtual ProcessType ID(void) const {return Weighted_Transition_Coagulation_ID;};

    // TOTAL RATE CALCULATION.

    // Returns the rate of the process for the given system.
    virtual real Rate(real t,         // Time.
                      const Cell &sys // System for which to calculate rate.
                      ) const;

    // RATE TERM CALCULATION.

    // Returns the number of rate terms for this process.
    virtual unsigned int TermCount(void) const;

    // Calculates the rate terms given an iterator to a real vector. The
    // iterator is advanced to the position after the last term for this
    // process.  Returns the sum of all rate terms.
    virtual real RateTerms(
        real t,                  // Time.
        const Cell &sys,       // Indicates true kernel (not majorant).
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

    //! Perform a coagulation with particles chosen according to the transition kernel
    virtual int Perform(
        real t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        unsigned int iterm,
        rng_type &rng) const;

    //! Write the object to a binary stream.
    virtual void Serialize(std::ostream &out) const;

protected:
    //! Transition coagulation kernel between two particles
    virtual real CoagKernel(
        const Particle &sp1, // First particle.
        const Particle &sp2, // Second particle.
        const Cell &sys
        ) const;

    //! Majorant coagulation kernel between two particles
    virtual real MajorantKernel(
        const Particle &sp1, // First particle.
        const Particle &sp2, // Second particle.
        const Cell &sys,
        const MajorantType maj) const;

private:

    // More efficient rate routine for coagulation only.
	// All parameters required to calculate rate terms
    // passed as arguments.
    real RateTerms(
    const TreeTransCoagWeightedCache &data, // Particle model data.
    real n,     // Number of particles.
    real sqrtT, // Square root of the temperature
    real T_mu,  // T / viscosity of air.
    real MFP,   // Gas mean-free path.
    real vol,   // System sample volume.
    fvector::iterator &iterm // Iterator to first coagulation term.
    ) const;

    // COAGULATION KERNEL ROUTINES.

    // Returns the free-molecular coagulation kernel value for the
    // two given particles.  Can return either the majorant or
    // true kernel.
    real FreeMolKernel(
    const Particle &sp1, // First particle.
    const Particle &sp2, // Second particle.
    real T,              // Temperature.
    real P,              // Pressure.
    const bool maj       // true=majorant kernel, false=true kernel.
    ) const;

    // Returns the slip-flow coagulation kernel value for the
    // two given particles.
    real SlipFlowKernel(
    const Particle &sp1, // First particle.
    const Particle &sp2, // Second particle.
    real T,              // Temperature.
    real P,              // Pressure.
    const bool maj       // true=majorant kernel, false=true kernel.
    ) const;

	/* Coagulation rate types. these define how the rate is
	 * calculated and how the particles are chosen.
	 */
	static const unsigned int TYPE_COUNT = 11;
	enum TermType {
        FreeMol1,
        FreeMol2,
        FreeMol3,
        FreeMol4,
        SlipFlow1,
        SlipFlow2,
        SlipFlow3,
        SlipFlow4,
        SlipFlow5,
        SlipFlow6,
        SlipFlow7
	};

    /* Free-molecular enhancement factor.  Currently hardcoded
     *  for soot particles (m_efm = 2.2).
     */
    static const real m_efm;

    //! Specify what to do with weights on coagulation
    CoagWeightRule m_CoagWeightRule;


}; // class WeightedTransitionCoagulation

} // namespace Processes

} // namespace Sweep

#endif	/* SWP_WEIGHTED_TRANSCOAG_H */
