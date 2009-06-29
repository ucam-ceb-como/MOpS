/*!
 * \file   swp_transport_process.h
 * \author Robert I A Patterson
 *
 * \brief  Process for stochastic transport
 */
#ifndef SWP_TRANSPORT_PROCESS_H
#define	SWP_TRANSPORT_PROCESS_H

#include "swp_process.h"
#include "swp_transport_outflow.h"

namespace Sweep {
namespace Processes {

//! Transport of particles by jumps
class TransportProcess : public Process {
public:

    //! Set the ID number of the particle property used for rate calculation
    void SetPropertyID(
        unsigned int id,
        SubModels::SubModelType modelid = SubModels::BasicModel_ID
        );


    //! Rate of the process for the given system.
    virtual real Rate(
        real t,
        const Cell &sys
        ) const;

    // RATE TERM CALCULATIONS.
    //   These routines return the individual rate terms for a
    //   process, which may have multiple terms (e.g. condensation).

    //! Number of rate terms for this process.
    virtual unsigned int TermCount(void) const;

    //! Puts rates terms into a vector
    virtual real RateTerms(
        real t,                  // Time.
        const Cell &sys,         // System for which to calculate rate terms.
        fvector::iterator &iterm // Iterator to the first term.
        ) const;

    // PERFORMING THE PROCESS.

    //! Performs the process on the given system.
    /*!
     * \param       t       Time
     * \param       sys     System to update
     * \param       iterm   Process term responsible for this event
     * \param       out     Details of any particle being transported out of system
     *
     * \return      0 on success, otherwise negative.
     */
    virtual int Perform(
        real t,
        Cell &sys,
        unsigned int iterm = 0,
        TransportOutflow *out = 0
        ) const;

private:
    //! Direction of the transport
    TransportDirection m_direction;

    //! Particle property to which the rate of the process is proportional.
    unsigned int m_pid;

    //! Particle model for which the above particle property ID is valid.
    SubModels::SubModelType m_modelid;



};

} // namespace Processes
} // namespace Sweep

#endif	/* SWP_TRANSPORT_PROCESS_H */

