/*!
 * @file    swp_surfvol_silica_primary.h
 * @author  William J Menz
 * @brief   Adaptation of the silica model to surfvol type space
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      A derived class of SurfVolPrimary, SurfVolSilicaPrimaryPrimary provides
 *      an adaptation of the silica model. The type space is given by:
 *      n_Si, n_O, n_OH, Surface area.
 *
 *      It is almost identical to the SurfVolPrimary model, except
 *      the sintering rate is cached to allow calculation of the
 *      interparticle reaction.
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
 *      sweepc is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public License
 *      as published by the Free Software Foundation; either version 2
 *      of the License, or (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this program; if not, write to the Free Software
 *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *      02111-1307, USA.
 *
 *   Contact:
 *      Prof Markus Kraft
 *      Dept of Chemical Engineering
 *      University of Cambridge
 *      New Museums Site
 *      Pembroke Street
 *      Cambridge
 *      CB2 3RA, UK
 *
 *      Email:       mk306@cam.ac.uk
 *      Website:     http://como.cheng.cam.ac.uk
*/

#ifndef SWP_SURFVOLSILICA_PRIMARY_H_
#define SWP_SURFVOLSILICA_PRIMARY_H_

#include "swp_surfvol_primary.h"

namespace Sweep {

namespace AggModels {

class SurfVolSilicaPrimary: public Sweep::AggModels::SurfVolPrimary {
public:

    //! Initialising constructor
    SurfVolSilicaPrimary(
        double time,
        const Sweep::ParticleModel &model
    );

    //! Copy constructor
    SurfVolSilicaPrimary(const SurfVolSilicaPrimary &copy);

    //! Stream-reading constructor
    SurfVolSilicaPrimary(
        std::istream &in,
        const Sweep::ParticleModel &model
    );

    //! Destructor
    ~SurfVolSilicaPrimary();

    //! Assignment operator
    SurfVolSilicaPrimary &operator=(const SurfVolSilicaPrimary &rhs);

    //! Clone particle
    SurfVolSilicaPrimary* const Clone() const;

    //! Returns the aggregation model which this primary describes.
    AggModels::AggModelType AggID() const;

    //! Sinters a primary for a given time
    void Sinter(
        double dt,
        Cell &sys,
        const Processes::SinteringModel &model,
        rng_type &rng,
        double wt
        );

    //! Updates properties of the particle
    void UpdateCache();

    //! Gets the number of active sites (hydroxides)
    double GetSites() const { return GetComponent("hydroxide"); }

    //! Gets the sintering rate of the particle
    double GetSintRate() const {return m_sinter_rate; }

    //! Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    //! Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,
        const Sweep::ParticleModel &model
        );

private:
    // Protect default constructor
    SurfVolSilicaPrimary();

    //! Gets the value of one of the chemical components
    double GetComponent(std::string name) const;

    //! Sets the value of one of the chemical components
    void SetComponent(std::string name, double val);

    // Sintering rate of particles
    double m_sinter_rate;
};

}

}

#endif /* SWP_SURFVOLSILICA_PRIMARY_H_ */
