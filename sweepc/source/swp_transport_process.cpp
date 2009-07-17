/*!
 * \file   swp_transport_process.h
 * \author Robert I A Patterson
 *
 * \brief  Process for stochastic transport
 */

#include "swp_transport_process.h"

#include "swp_cell.h"

/*!
 * Set the ID number of the particle property to which
 * the rate of this process is proportional.
 * Copied from the surface reaction.
// Calculates the rate of multiple inceptions given a
// vector of inceptions and an iterator to a vector of
// reals for output.
 *
 *\param[in]            id          Particle property to which rate is proportional
 *\param[in,optional]   modelid     Model in which property applies
 */
void Sweep::Processes::TransportProcess::SetPropertyID(unsigned int id, SubModels::SubModelType modelid) {
    m_pid     = id;
    m_modelid = modelid;
}

/*!
 * Calculate the rates of all the transport processes from a vector of pointers.
 * The rates for each process put into the rates vector at consecutive locations
 * starting at index start.
 *
 *\param[in]        t       Time at which rates should be calculated
 *\param[in]        sys     System for which to calculate the rates
 *\param[in]        icns    Vector of pointers to tranport process for which rates are to be calculated
 *\param[in,out]    rates   Vector to store individual process rates
 *\param[in]        start   Index in vector of location to place first rate
 *
 *\return       Combined total rate of all supplied processes
 */
Sweep::real Sweep::Processes::TransportProcess::CalcRates(real t, const Cell &sys, const TransportPtrVector &icns,
                          fvector &rates, unsigned int start)
{
    TransportPtrVector::const_iterator p;
    fvector::iterator i = (rates.begin()+start);
    real sum = 0.0;
    for (p=icns.begin(); p!=icns.end(); ++p,++i) {
        *i = (*p)->Rate(t, sys);
        sum += *i;
    }
    return sum;
}
