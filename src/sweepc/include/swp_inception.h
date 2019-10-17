/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Definition of an inception process.  Inceptions are modelled as the coagulation
    of 2 gas-phase species.  The rate is given by the collision kernel of these species,
    therefore one must provide their masses and diameters.  Kernel parameters are
    calculated internally.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#ifndef SWEEP_INCEPTION_H
#define SWEEP_INCEPTION_H

#include "swp_params.h"
#include "swp_process.h"
#include "swp_process_type.h"
#include "swp_particle.h"
#include "swp_cell.h"
#include <vector>

namespace Geometry
{
    // Forward declaration
    class LocalGeometry1d;
}

namespace Sweep
{
// Forward declare the Mechanism class.
class Mechanism;

namespace Processes
{
class Inception;
typedef std::vector<Inception*> IcnPtrVector;

/*!
 * \brief Base class for inception processes.
 *
 */
class Inception : public Process
{
public: 
    // Constructors.
    Inception(const Sweep::Mechanism &mech); // Initialising constructor.
    Inception(const Inception &copy);        // Copy constructor.
    Inception(                               // Stream-reading constructor.
        std::istream &in,                    //  - Input stream.
        const Sweep::Mechanism &mech         //  - Parent mechanism.
        );

    // Destructors.
    ~Inception(void);

    // Operators.
    Inception &operator=(const Inception &rhs);


    // INCEPTION KERNEL.


    // PROPERTIES OF INCEPTED PARTICLES.

    // Returns the composition vector of the new particle.
    const fvector &ParticleComp(void) const;

    // Returns the amount of the ith component of the new particle.
    double ParticleComp(unsigned int i) const;

    // Sets the amount of the ith component in the new particle.
    void SetParticleComp(unsigned int i, double comp);

    // Returns the tracker variable vector of the new particle.
    const fvector &ParticleTrackers(void) const;

    // Returns the value of the ith tracker variable of the
    // new particle.
    double ParticleTrackers(unsigned int i) const;

    // Sets the new particle tracker variable vector.
    void SetParticleTrackers(const fvector &track);

    // Sets the value of the ith tracker variable in the
    // new particle.
    void SetParticleTracker(unsigned int i, double track);

	// TOTAL RATE CALCULATIONS.

    // Calculates the rate of multiple inceptions given a
    // vector of inceptions and an iterator to a vector of
    // reals for output.
    static double CalcRates(
        double t,                   // Time.
        const Cell &sys,          // System for which to calculate rates.
        const Geometry::LocalGeometry1d &local_geom,
        const IcnPtrVector &icns, // Vector of inception processes.
        fvector &rates,           // Output rates vector.
        unsigned int start = 0    // Vector position to start at in vector rates.
        );


    // Creates a copy of the inception.
    Inception *const Clone(void) const=0;

    //! Perform an event to transfer mass from gasphase to particle pahse. This function is only used for PAH-PP model
    // so it will be implemented in the PAHInception class
    int AddPyrene(
    int i,
    double t,
    Cell &sys,
    rng_type &rng) const;

    // READ/WRITE/COPY.

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(
        std::istream &in,            // Input stream.
        const Sweep::Mechanism &mech // Parent mechanism.
        );

protected:
    /*! Default constructor is protected to prevent an inception being
     * defined without knowledge of the parent mechanism.*/
    Inception(void);

private:

    //! Composition of newly incepted particles.
    fvector m_newcomp;

    //! Initial tracker values for newly incepted particles
    fvector m_newvals;
};
} // namespace Processes
} // namespace Sweep

#endif
