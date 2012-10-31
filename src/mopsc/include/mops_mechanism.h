/*
  Author(s):      Matthew Celnik (msc37)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Mechanism class combines the gas-phase mechanism of Sprog with
    the particle mechanism of Sweep.
    
    Additionally the Mixture class includes a description of a particle
    population which can be solved with Sweep.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
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

#ifndef MOPS_MECHANISM_H
#define MOPS_MECHANISM_H

#include "mops_params.h"
#include "gpc_mech.h"
#include "swp_mechanism.h"

namespace Mops
{
class Mechanism
{
public:
    // Constructors.
    Mechanism(void); // Default constructor.

    // Destructors.
    ~Mechanism(void); // Default destructors.

    // PARTICLE MECHANISM.

    // Returns a reference (non-const) to the particle mechanism.
    Sweep::Mechanism &ParticleMech(void) {return m_pmech;}
    const Sweep::Mechanism &ParticleMech(void) const {return m_pmech;}

	/*
	// Returns a reference (non-const) to the surface mechanism. (Added by mm864)
    Sprog::SurfaceMechanism &SurfaceMech(void) {return m_smech;}
    const Sprog::SurfaceMechanism &SurfaceMech(void) const {return m_smech;}
	*/
	
    //! Access the gas phase mechanism
    Sprog::Mechanism &GasMech() {return m_gmech;}
    //! Access the gas phase mechanism
    const Sprog::Mechanism &GasMech() const {return m_gmech;}

    // READ/WRITE/COPY FUNCTIONS.

    // Writes the mechanism to a binary data stream.
    void Serialize(std::ostream &out) const;

    // Reads the mechanism data from a binary data stream.
    void Deserialize(std::istream &in);

private:
    //! Gas phase mechanism
    Sprog::Mechanism m_gmech;

	//! Surface mechanism (Added by mm864)
	// Sprog::SurfaceMechanism m_smech;
	
    //! The particle mechanism.
    Sweep::Mechanism m_pmech;
};
}

#endif
