/*!
 * \file   swp_transport_process.h
 * \author Robert I A Patterson
 *
 * \brief  Process for stochastic transport
 */
#ifndef SWP_TRANSPORT_PROCESS_H
#define	SWP_TRANSPORT_PROCESS_H

#include "swp_process.h"

#include "local_geometry1d.h"

#include <vector>

namespace Sweep {

namespace Transport
{
    // Forward declaration of type to carry particles for onward transport
    struct TransportOutflow;
}

namespace Processes {

// Forward declaration
class TransportProcess;
//! Vector of pointers to transport processes
typedef std::vector<TransportProcess*> TransportPtrVector;


//! Transport of particles by jumps
class TransportProcess : public Process {
public:
    //! Create a copy of the transport process.
    virtual TransportProcess *const Clone(void) const = 0;

    //! Set the rate scaling factor
    void SetA(const real a) {m_a = a;}

    //! Get the rate scaling factor
    real A() const {return m_a;}

    //! Backwards compatiliblity method, always returns 0
    virtual real Rate(
        real t,
        const Cell &sys
        ) const;

    //! Rate of the process
    virtual real Rate(
        real t,
        const Cell &sys,
        const Geometry::LocalGeometry1d& local_geom
        ) const = 0;

    // SINGLE PARTICLE RATE CALCULATIONS.

    /*!
     * Calculate the particle dependent part of the process rate
     * for a single particle.  Scaling factors depending on the
     * chemical environment such as temperature, pressure and
     * species concentrations will not be included.
     *
     * This method (and its concrete implementations) do not seem
     * to need geometry information, because particle properties
     * are used in the same way for all possible directions.
     *
     *@param[in]        t       Time at which rate is to be calculated
     *@param[in]        sys     System containing particle
     *@param[in]        sp      Particle for which to calculate rate
     *
     *@return       Particle dependent part of process rate
     */
    virtual real Rate(
        real t,
        const Cell &sys,
        const Particle &sp
        ) const = 0;


    /*!
     * Calculate the particle dependent part of the process majorant rate
     * for a single particle.  Scaling factors depending on the
     * chemical environment such as temperature, pressure and
     * species concentrations will not be included.
     *
     *@param[in]        t       Time at which rate is to be calculated
     *@param[in]        sys     System containing particle
     *@param[in]        sp      Particle for which to calculate rate
     *
     *@return       Particle dependent part of process majorant rate
     */
    virtual real MajorantRate(
        real t,
        const Cell &sys,
        const Particle &sp
        ) const = 0;

    // RATE TERM CALCULATIONS.
    //   These routines return the individual rate terms for a
    //   process, which may have multiple terms (e.g. condensation).

    //! Puts rate terms into a vector assuming all transport blocked
    virtual real RateTerms(
        real t,
        const Cell &sys,
        fvector::iterator &iterm 
        ) const;

    //! Puts rate terms into a vector for transport in permitted directions
    virtual real RateTerms(
        real t,
        const Cell &sys,
        //const Sweep::Transport::DirectionMask &transport_mask,
        const Geometry::LocalGeometry1d& local_geom,
        fvector::iterator &iterm
        ) const = 0;

    //! Calculate rates for a vector of transport processes
    static real CalcRates(
        real t,                   // Time.
        const Cell &sys,          // System for which to calculate rates.
        //const Sweep::Transport::DirectionMask &transport_mask,
        const Geometry::LocalGeometry1d& local_geom,
        const TransportPtrVector &trans, // Vector of transport processes.
        fvector &rates,           // Output rates vector.
        unsigned int start = 0    // Vector position to start at in vector rates.
        );

    //! Remove a particle and put it into the TransportOutflow object for onward transport
    int Outflow(
        real t,
        Cell &sys,
        const Geometry::LocalGeometry1d& local_geom,
        const int particle_index,
        const Geometry::Direction &direction,
        real(*rand_u01)(),
        Transport::TransportOutflow *out
        ) const;

    /*!
     *\brief    Transport processes are not deferred (for now)
     *
     *\return   True if the process is deferred
     *
     * This method may eventually be moved into concrete derived classes.
     */
    bool IsDeferred() const {return false;}

private:
    //! Constant factor multiplying the rate (must not be negative)
    real m_a;
};

} // namespace Processes
} // namespace Sweep

#endif	/* SWP_TRANSPORT_PROCESS_H */

