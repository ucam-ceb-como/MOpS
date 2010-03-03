/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ARSSC_Inception class declared in the
    swp_arssc_inception.h header file.

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

#include "swp_arssc_inception.h"
#include "swp_arssc_model.h"
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::Processes;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
ARSSC_Inception::ARSSC_Inception(void)
: DimerInception(), ARSSC_Process()
{
    m_name = "ARS-SC Inception";
}

// Initialising constructor.
ARSSC_Inception::ARSSC_Inception(const Sweep::Mechanism &mech)
: DimerInception(mech), ARSSC_Process()
{
    m_name = "ARS-SC Inception";
}

// Copy constructor.
ARSSC_Inception::ARSSC_Inception(const ARSSC_Inception &copy)
{
    *this = copy;
}

// Stream-reading constructor.
ARSSC_Inception::ARSSC_Inception(std::istream &in, const Sweep::Mechanism &mech)
{
    Deserialize(in, mech);
}

// Default destructor.
ARSSC_Inception::~ARSSC_Inception(void)
{
}

// OPERATOR OVERLOADS.

// Assignment operator.
ARSSC_Inception &ARSSC_Inception::operator =(const ARSSC_Inception &rhs)
{
    if (this != &rhs) {
        DimerInception::operator=(rhs);
        ARSSC_Process::operator =(rhs);
        m_sites.assign(rhs.m_sites.begin(), rhs.m_sites.end());
    }
    return *this;
}


// PERFORMING THE PROCESS.

// Performs the process on the given system.  The responsible rate term is given
// by index.  Returns 0 on success, otherwise negative.
int ARSSC_Inception::Perform(real t, Cell &sys, unsigned int iterm, Transport::TransportOutflow*) const  
{
    // This routine performs the inception on the given chemical system.

    // Create a new particle of the type specified
    // by the system ensemble.
    Particle *sp = m_mech->CreateParticle(t);
    
    // Set the new particle's aromatic sites.
    SubModels::ARSSC_Model *ars = 
        static_cast<SubModels::ARSSC_Model*const>(sp->Primary()->SubModel(SubModels::ARSSC_Model_ID));
    for (unsigned int i=0; i!=m_sites.size(); ++i) {
        ars->SetSiteCount((SubModels::ARSSC_Model::SiteType)i, m_sites[i]);
    }

    // Initialise the new particle.
    sp->Primary()->SetComposition(m_newcomp);
    sp->Primary()->SetValues(m_newvals);
    sp->UpdateCache();

    // Add particle to system's ensemble.
    sys.Particles().Add(*sp, Sweep::irnd);

    // Update gas-phase chemistry of system.
    adjustGas(sys);

    return 0;
}


// READ/WRITE/COPY.

// Creates a copy of the inception.
ARSSC_Inception *const ARSSC_Inception::Clone(void) const
{
    return new ARSSC_Inception(*this);
}

// Returns the process type.  Used to identify different
// processes and for serialisation.
ProcessType ARSSC_Inception::ID(void) const {return ARSSC_Inception_ID;}

// Writes the object to a binary stream.
void ARSSC_Inception::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Serialize base class.
        DimerInception::Serialize(out);
        ARSSC_Process::Serialize(out);
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ARSSC_Inception::Serialize).");
    }
}

// Reads the object from a binary stream.
void ARSSC_Inception::Deserialize(std::istream &in, const Sweep::Mechanism &mech)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        switch (version) {
            case 0:
                // Deserialize base class.
                DimerInception::Deserialize(in, mech);
                ARSSC_Process::Deserialize(in);

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ARSSC_Inception::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ARSSC_Inception::Deserialize).");
    }
}
