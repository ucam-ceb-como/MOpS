/*!
 * \file   swp_transport_process.h
 * \author Robert I A Patterson
 *
 * \brief  Process for stochastic transport
 */

#include "swp_transport_process.h"

#include "swp_cell.h"
#include "swp_mechanism.h"
#include "swp_transport_outflow.h"


/*!
 * Calculate the rates of all the transport processes from a vector of pointers.
 * The rates for each process put into the rates vector at consecutive locations
 * starting at index start.
 *
 *\param[in]        t           Time at which rates should be calculated
 *\param[in]        sys         System for which to calculate the rates
 *\param[in]        local_geom  Details of immediately surrounding grid
 *\param[in]        icns        Vector of pointers to tranport process for which rates are to be calculated
 *\param[in,out]    rates       Vector to store individual process rates
 *\param[in]        start       Index in vector of location to place first rate
 *
 *\return       Combined total rate of all supplied processes
 */
Sweep::real Sweep::Processes::TransportProcess::CalcRates(real t, const Cell &sys,
                                                          const Geometry::LocalGeometry1d& local_geom,
                                                          const TransportPtrVector &icns,
                                                          fvector &rates, unsigned int start)
{
    TransportPtrVector::const_iterator p;
    fvector::iterator i = (rates.begin()+start);
    real sum = 0.0;
    for (p=icns.begin(); p!=icns.end(); ++p,++i) {
        *i = (*p)->Rate(t, sys, local_geom);
        sum += *i;
    }
    return sum;
}

/*!
 * Default the process rate to 0, because no Transport::DirectionMask
 * has been provided.  To get a non-zero rate call the overload of this
 * function, which takes a Transport::DirectionMask as its third argument.
 *
 *\param[in]        t       Time at which rate should be calculated
 *\param[in]        sys     System for which to calculate the process rate
 *
 *\return           0
 */
Sweep::real Sweep::Processes::TransportProcess::Rate(real t, const Cell &sys) const {
    return 0;
}

/*!
 * Append the rate to a vector, assuming transport in all directions is blocked,
 * which should guarantee a rate of 0.  This method is provided for compatibility
 * with other processes.
 *
 *\param[in]        t       Time at which rate should be calculated
 *\param[in]        sys     System for which to calculate the rate
 *\param[in,out]    iterm   Position in vector at which to insert the rate
 *
 *\return       Total rate of process
 */
 Sweep::real Sweep::Processes::TransportProcess::RateTerms(real t, const Cell &sys,
                                                           fvector::iterator &iterm) const {
     return RateTerms(t, sys, Geometry::LocalGeometry1d(), iterm);
 }

/*!
 * \param[in]       t               Time
 * \param[in]       local_geom      Details of geometry around current location
 * \param[in,out]   sys             System from which to remove particle
 * \param[in]       particle_index  Index of particle to remove
 * \param[out]      out             Details of particle being transported out of system
 *
 * \return      0 on success, otherwise negative.
 */
int Sweep::Processes::TransportProcess::Outflow(const real t, Cell &sys,
                                                const Geometry::LocalGeometry1d& local_geom,
                                                const int particle_index,
                                                const Geometry::Direction &direction,
                                                Transport::TransportOutflow *out) const {
    // Check for a valid particle (particle_index>=0) before proceeding
    if (particle_index >= 0) {
        out->particle = sys.Particles().At(particle_index);

        // Get the majorant rate, which is proportional to the contribution this
        // particle made to the overall simulation jump rate for this process.
        const real majorantRate = MajorantRate(t, sys, *(out->particle));

        if (m_mech->AnyDeferred()) {
            // Update particle with deferred processes.
            m_mech->UpdateParticle(*(out->particle), sys, t);
        }

        // Check that the particle is still valid.
        if (out->particle->IsValid()) {
            // Get the true process rate (after updates).
            const real trueRate = Rate(t, sys, *(out->particle));

            // Check whether the event is  ficticious
            if (Ficticious(majorantRate, trueRate)) {
                // Fictitious event so update the particle and leave it in
                // its original location
                sys.Particles().Update(particle_index);
                out->particle = NULL;
            }
            else {
                // Remove particle from the cell, but do not delete the particle
                sys.Particles().Remove(particle_index, false);

                // Pass details of the transport back to the caller to handle onward routing
                out->weight = 1.0 / sys.SampleVolume();
                out->destination = local_geom.calcDestination(direction);
            }
        }
        else {
            // Particle was removed during the deferred update, so it cannot
            // be transported on.
            out->particle = NULL;
            sys.Particles().Remove(particle_index, true);
        }

        return 0;
    }

    // This will be < 0, because if particle_index >= 0 function will be exited in above if block
    return particle_index;
}

/*!
 * \param       t       Time
 * \param       sys     System to update
 * \param       iterm   Process term responsible for this event
 * \param       out     Details of any particle being transported out of system
 *
 * \return      0 on success, otherwise negative.
 *
 * This method is provided to implement a pure virtual method in the parent
 * class.  It calls through to an non-virtual overload with default geometry
 * information.
 */
 int Sweep::Processes::TransportProcess::Perform(real t, Cell &sys,
                                                 unsigned int iterm,
                                                 Transport::TransportOutflow *out) const {
      return Perform(t, sys, Geometry::LocalGeometry1d(), iterm, out);
 }