/*!
 * \file   swp_transport_outflow.h
 * \author Robert I A Patterson
 *
 * \brief Types to carry details of a transported particle
 */

#ifndef SWP_TRANSPORT_OUTFLOW_H
#define	SWP_TRANSPORT_OUTFLOW_H

#include "swp_params.h"

namespace Sweep {
    // Forward declaration
    class Particle;

namespace Transport {

/*!
 *Collect together all the information needed when a particle is transported
 *between two ensembles.  The structure will never own the particle it may
 *point to, releasing the memory is always the responsibility of the calling
 *code.
 */
struct TransportOutflow {
    TransportOutflow()
        : particle(0), weight(0.0), destination(-1)
    {}

    //! Pointer to particle that is being moved
    Particle* particle;

    //! Statistical weight of particle begin moved
    real weight;

    //! Index of destination cell
    int destination;
};

} // namespace Transport
} // namespace Sweep

#endif	/* SWP_TRANSPORT_OUTFLOW_H */

