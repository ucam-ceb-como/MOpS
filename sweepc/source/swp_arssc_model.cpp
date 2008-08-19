/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the ARSSC_Model class declared in the
    swp_arssc_model.h header file.

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

#include "swp_arssc_model.h"
//#include "swp_arssc_cache.h"
#include "swp_primary.h"
#include "rng.h"
#include <cmath>
#include <stdexcept>

using namespace Sweep;
using namespace Sweep::SubModels;
using namespace std;

// STATIC VARIABLES.

unsigned int ARSSC_Model::m_ipah = 0;

// Neighbouring site probability weights.  These are static
// because they are common to all instances of the ARSSC model
// data.  Weights are stored in order Edge/Zigzag/Armchair/Bay.
fvector ARSSC_Model::m_fe_wt(ARSSC_Model::SiteTypeCount, 0.25);
fvector ARSSC_Model::m_zz_wt(ARSSC_Model::SiteTypeCount, 0.25);
fvector ARSSC_Model::m_ac_wt(ARSSC_Model::SiteTypeCount, 0.25);
fvector ARSSC_Model::m_bay_wt(ARSSC_Model::SiteTypeCount, 0.25);

// Combined-site correlation functions.  There are currently
// four hardcoded combined sites.
Maths::Functional *ARSSC_Model::m_r6_func    = NULL;
Maths::Functional *ARSSC_Model::m_ac_r6_func = NULL;
Maths::Functional *ARSSC_Model::m_r5_ed_func = NULL;
Maths::Functional *ARSSC_Model::m_r5_ac_func = NULL;

// Combine-site correlation independent variables (site types).
ARSSC_Model::SiteType ARSSC_Model::m_r6_var    = ARSSC_Model::FreeEdge;
ARSSC_Model::SiteType ARSSC_Model::m_ac_r6_var = ARSSC_Model::Armchair;
ARSSC_Model::SiteType ARSSC_Model::m_r5_ed_var = ARSSC_Model::FreeEdge;
ARSSC_Model::SiteType ARSSC_Model::m_r5_ac_var = ARSSC_Model::Armchair;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (private).
ARSSC_Model::ARSSC_Model()
: m_sites(ARSSC_Model::SiteTypeCount, 0.0)
{
}

// Initialising constructor (default).
ARSSC_Model::ARSSC_Model(Primary &parent)
: m_sites(ARSSC_Model::SiteTypeCount, 0.0)
{
    m_parent = &parent;
}

// Copy constructor.
ARSSC_Model::ARSSC_Model(const ARSSC_Model &copy)
{
    *this = copy;
}

// Stream-reading constructor.
ARSSC_Model::ARSSC_Model(std::istream &in, Primary &parent)
: m_sites(ARSSC_Model::SiteTypeCount, 0.0)
{
    Deserialize(in, parent);
}

// Destructor.
ARSSC_Model::~ARSSC_Model(void)
{
    // Nothing special to destruct.
}


// OPERATORS.

// Assignment operator.
ARSSC_Model &ARSSC_Model::operator=(const ARSSC_Model &rhs)
{
    if (this != &rhs) {
        m_sites.assign(rhs.m_sites.begin(), rhs.m_sites.end());
    }
    return *this;
}


// SITE COUNTS.

// Returns the count for the given site type.
real ARSSC_Model::SiteCount(SiteType type) const {return m_sites[type];}

// Returns the number of free-edges.
real ARSSC_Model::FreeEdgeCount(void) const {return m_sites[FreeEdge];}

// Returns the number of armchairs.
real ARSSC_Model::ArmchairCount(void) const {return m_sites[Armchair];}

// Returns the number of zigzags.
real ARSSC_Model::ZigzagCount(void) const {return m_sites[Zigzag];}

// Returns the number of bays.
real ARSSC_Model::BayCount(void) const {return m_sites[Bay];}

// Returns the number of R5s.
real ARSSC_Model::R5Count(void) const {return m_sites[R5];}

// Sets the site count of the given site type.
void ARSSC_Model::SetSiteCount(SiteType type, real n) {m_sites[type] = n;}

// Adds n to the site count of the given site type.
void ARSSC_Model::AddSites(SiteType type, real n) {
    m_sites[type] = max(m_sites[type]+=n, 0.0);
}

// Subtracts n from the site count of the given site type.
void ARSSC_Model::RemoveSites(SiteType type, real n) {
    m_sites[type] = max(m_sites[type]-=n, 0.0);
}


// Updates the site counts using the rules for neighbouring
// sites.  A negative count for the number of sites to update
// indicates site decrementing.
void ARSSC_Model::AdjustNeighbourSites(SiteType parent, int n)
{
    if (n >= 0) {
        for (int i=0; i!=n; ++i) {
            SiteType k = SelectNeighbour(parent);
            switch (k) {
                // For free-edges, zigzags and armchairs site
                // incrementing occurs.
                case FreeEdge:
                case Zigzag:
                case Armchair:
                    if (m_sites[k] >= 1.0) {
                        m_sites[k]   -= 1.0;
                        m_sites[k+1] += 1.0;
                    }
                    break;
                case Bay:
                    // For bay sites we do nothing.
                    break;
                default:
                    // Something has gone wrong to bring us here.
                    break;
            }

        }
    } else {
        for (int i=0; i!=n; --i) {
            SiteType k = SelectNeighbour(parent, true);
            switch (k) {
                // For zigzags and armchairs site
                // decrementing occurs.
                case Zigzag:
                case Armchair:
                    if (m_sites[k] >= 1.0) {
                        m_sites[k]   -= 1.0;
                        m_sites[k-1] += 1.0;
                    }
                    break;
                case Bay:
                    // For bay sites we do nothing.
                    break;
                default:
                    // Something has gone wrong to bring us here.
                    // This includes if a free-edge was chosen, which 
                    // should be impossible.
                    break;
            }

        }
    }
}

// Updates the site counts using the rules for neighbouring
// sites.  The neighbour sites weights are provided as a vector,
// rather than using the weights for a given parent site.  A negative
// count for the number of sites to update indicates site decrementing.
void ARSSC_Model::AdjustNeighbourSites(fvector wts, int n)
{
    if (n >= 0) {
        for (int i=0; i!=n; ++i) {
            SiteType k = SelectNeighbour(wts);
            switch (k) {
                // For free-edges, zigzags and armchairs site
                // incrementing occurs.
                case FreeEdge:
                case Zigzag:
                case Armchair:
                    if (m_sites[k] >= 1.0) {
                        m_sites[k]   -= 1.0;
                        m_sites[k+1] += 1.0;
                    }
                    break;
                case Bay:
                    // For bay sites we do nothing.
                    break;
                default:
                    // Something has gone wrong to bring us here.
                    break;
            }

        }
    } else {
        for (int i=0; i!=n; --i) {
            SiteType k = SelectNeighbour(wts, true);
            switch (k) {
                // For zigzags and armchairs site
                // decrementing occurs.
                case Zigzag:
                case Armchair:
                    if (m_sites[k] >= 1.0) {
                        m_sites[k]   -= 1.0;
                        m_sites[k-1] += 1.0;
                    }
                    break;
                case Bay:
                    // For bay sites we do nothing.
                    break;
                default:
                    // Something has gone wrong to bring us here.
                    // This includes if a free-edge was chosen, which 
                    // should be impossible.
                    break;
            }

        }
    }
}


// PAH COUNT.

// Returns the number of PAHs.
real ARSSC_Model::PAH_Count(void) const {return m_parent->Values(m_ipah);}

// Sets the number of PAHs.
void ARSSC_Model::SetPAH_Count(real n) {m_parent->SetValue(m_ipah, n);}

// Adds n to the PAH count.
void ARSSC_Model::AddPAHs(real n)
{
    m_parent->SetValue(m_ipah, max(m_parent->Values(m_ipah)+n, 0.0));
}

// Subtracts n from the PAH count.
void ARSSC_Model::RemovePAHs(real n)
{
    m_parent->SetValue(m_ipah, max(m_parent->Values(m_ipah)-n, 0.0));
}


// NEIGHBOURING SITE WEIGHTS.

// Returns the neighbouring site weights vector for
// free-edge sites.
const fvector &ARSSC_Model::FE_Neighbours(void) {return m_fe_wt;}

// Returns the neighbouring site weights vector for
// zigzag sites.
const fvector &ARSSC_Model::ZZ_Neighbours(void) {return m_zz_wt;}

// Returns the neighbouring site weights vector for
// armchair sites.
const fvector &ARSSC_Model::AC_Neighbours(void) {return m_ac_wt;}

// Returns the neighbouring site weights vector for
// bay sites.
const fvector &ARSSC_Model::BayNeighbours(void) {return m_bay_wt;}

// Sets the neighbour site (type y) weight for the
// given site type (type x).
void ARSSC_Model::SetNeighbour(SiteType x, SiteType y, real weight)
{
    switch (x) {
        case FreeEdge:
            m_fe_wt[(fvector::size_type)y] = weight;
            break;
        case Zigzag:
            m_zz_wt[(fvector::size_type)y] = weight;
            break;
        case Armchair:
            m_ac_wt[(fvector::size_type)y] = weight;
            break;
        case Bay:
            m_bay_wt[(fvector::size_type)y] = weight;
            break;
        default:
            //Do nothing.
            break;
    }
}

// Randomly selects a neighbour site for the given
// parent using the correct weights.  If the optional
// dec flag is set to true then a neighbour for the
// decrementing process is selected (free-edges are
// ignored).
ARSSC_Model::SiteType ARSSC_Model::SelectNeighbour(SiteType parent, bool dec) const
{
    // Select corrent site weights for the given parent
    // site type.
    static fvector wts(4, 0.0);
    switch (parent) {
        case FreeEdge:
            return SelectNeighbour(m_fe_wt, dec);
        case Zigzag:
            return SelectNeighbour(m_zz_wt, dec);
        case Armchair:
            return SelectNeighbour(m_ac_wt, dec);
        case Bay:
            return SelectNeighbour(m_bay_wt, dec);
        default:
            // Invalid parent site, chosen uniformly.
            wts[FreeEdge] = 0.25;
            wts[Zigzag]   = 0.25;
            wts[Armchair] = 0.25;
            wts[Bay]      = 0.25;
            return SelectNeighbour(wts, dec);
    }
}

// Randomly selects a neighbour site using the given custom
// weights.  If the optional
// dec flag is set to true then a neighbour for the
// decrementing process is selected (free-edges are
// ignored).
ARSSC_Model::SiteType ARSSC_Model::SelectNeighbour(const fvector &wts, bool dec) const
{
    static fvector probs(4, 0.0);
    probs[FreeEdge] = wts[FreeEdge];
    probs[Zigzag]   = wts[Zigzag];
    probs[Armchair] = wts[Armchair];
    probs[Bay]      = wts[Bay];

    // Multiply the site weights by the current site
    // counts.
    if (dec) {
        probs[FreeEdge] = 0.0;
    } else {
        probs[FreeEdge] *= m_sites[FreeEdge];
    }
    probs[Zigzag]   *= m_sites[Zigzag];
    probs[Armchair] *= m_sites[Armchair];
    probs[Bay]      *= m_sites[Bay];

    // Get the sum of the weighted counts.
    real probsum = probs[FreeEdge] + probs[Zigzag] + probs[Armchair] + probs[Bay];

    if (probsum > 0.0) {
        // Generate a random deviate.
        real r = rnd() * probsum;

        // Use a DIV algorithm to select a neighbour.
        if ((r-=probs[FreeEdge]) <= 0.0) {
            return FreeEdge;
        } else {
            if ((r-=probs[Zigzag]) <= 0.0) {
                return Zigzag;
            } else {
                if ((r-=probs[Armchair]) <= 0.0) {
                    return Armchair;
                } else {
                    if ((r-=probs[Bay]) <= 0.0) {
                        return Bay;
                    } else {
                        // An error has occured to be here.
                        throw std::runtime_error("Invalid site type chosen"
                                                 " (Sweep, ARSSC_Model::SelectNeighbour).");
                    }
                }
            }
        }
    } else {
        // There appear to be no sites.  This is an error.
        throw std::runtime_error("No sites are available for selection"
                                 " (Sweep, ARSSC_Model::SelectNeighbour).");
    }
}


// COMBINED-SITE CORRELATION FUNCTIONS.

// Returns the combined-site correlation functional
// for R6 rings.
Maths::Functional &ARSSC_Model::R6_Functional(void) {return *m_r6_func;}

// Returns the combined-site correlation functional
// for AC next to R6 rings.
Maths::Functional &ARSSC_Model::ACR6_Functional(void) {return *m_ac_r6_func;}

// Returns the combined-site correlation functional
// for R5 next to ED.
Maths::Functional &ARSSC_Model::R5ED_Functional(void) {return *m_r5_ed_func;}

// Returns the combined-site correlation functional
// for R5 next to AC.
Maths::Functional &ARSSC_Model::R5AC_Functional(void) {return *m_r5_ac_func;}

// Sets the combined-site correlation functional
// for R6 rings.  The given functional is stored
// by reference, so it would be unwise to delete
// the reference after calling this routine.
void ARSSC_Model::SetR6_Functional(Maths::Functional &func) {m_r6_func = &func;}

// Sets the combined-site correlation functional
// for AC next to R6 rings.  The given functional is stored
// by reference, so it would be unwise to delete
// the reference after calling this routine.
void ARSSC_Model::SetACR6_Functional(Maths::Functional &func) {m_ac_r6_func = &func;}

// Sets the combined-site correlation functional
// for R5 next to ED.  The given functional is stored
// by reference, so it would be unwise to delete
// the reference after calling this routine.
void ARSSC_Model::SetR5ED_Functional(Maths::Functional &func) {m_r5_ed_func = &func;}

// Sets the combined-site correlation functional
// for R5 next to AC.  The given functional is stored
// by reference, so it would be unwise to delete
// the reference after calling this routine.
void ARSSC_Model::SetR5AC_Functional(Maths::Functional &func) {m_r5_ac_func = &func;}


// COMBINED-SITE CORRELATION INDEPENDENT VARIABLES.

// Returns the independent variable (site type) used to 
// calculate R6 ring counts.
ARSSC_Model::SiteType ARSSC_Model::R6_Variable(void) {return m_r6_var;}

// Returns the independent variable (site type) used to 
// calculate AC next to R6 ring counts.
ARSSC_Model::SiteType ARSSC_Model::ACR6_Variable(void) {return m_ac_r6_var;}

// Returns the independent variable (site type) used to 
// calculate R5 next to ED counts.
ARSSC_Model::SiteType ARSSC_Model::R5ED_Variable(void) {return m_r5_ed_var;}

// Returns the independent variable (site type) used to 
// calculate R5 next to AC counts.
ARSSC_Model::SiteType ARSSC_Model::R5AC_Variable(void) {return m_r5_ac_var;}

// Sets the independent variable (site type) used to
// calculate R6 ring counts.
void ARSSC_Model::SetR6_Variable(ARSSC_Model::SiteType type) {m_r6_var = type;}

// Sets the independent variable (site type) used to
// calculate AC next to R6 ring counts.
void ARSSC_Model::SetACR6_Variable(ARSSC_Model::SiteType type) {m_ac_r6_var = type;}

// Sets the independent variable (site type) used to
// calculate R5 next to ED counts.
void ARSSC_Model::SetR5ED_Variable(ARSSC_Model::SiteType type) {m_r5_ed_var = type;}

// Sets the independent variable (site type) used to
// calculate R5 next to AC counts.
void ARSSC_Model::SetR5AC_Variable(ARSSC_Model::SiteType type) {m_r5_ac_var = type;}


// PARTICLE UPDATES.

// Updates the model data according to the rules of the model
// given the changes to particle composition and values.
// Performs the update n times.
void ARSSC_Model::Update(const fvector &dcomp, const fvector &dval, 
                         unsigned int n)
{
    // Do nothing, updates to sites are performed seperately.
}

// Updates the data according to the model rules for
// coagulation.
ARSSC_Model &ARSSC_Model::Coagulate(const SubModel &rhs)
{
    // Attempt to cast sub-model as ARS-SC model.
    const ARSSC_Model *ars = dynamic_cast<const ARSSC_Model*>(&rhs);

    if (ars) {
        // Sum up the site counts.
        for (unsigned int j=0; j!=SiteTypeCount; ++j) {
            m_sites[j] += ars->m_sites[j];
        }
    }
    return *this;
}


// CACHE UPDATES.

// Updates the sub-model cached variables.
void ARSSC_Model::UpdateCache(void)
{
    // Determine the combined-site counts.
    m_sites[R6] = m_r6_func->Evaluate(m_sites[m_r6_var]);
    m_sites[ACR6] = m_ac_r6_func->Evaluate(m_sites[m_ac_r6_var]);
    m_sites[R5ED] = m_r5_ed_func->Evaluate(m_sites[m_r5_ed_var]);
    m_sites[R5AC] = m_r5_ac_func->Evaluate(m_sites[m_r5_ac_var]);
}

// Creates a model cache which subscribes to this particle model.
SubModelCache *const ARSSC_Model::CreateCache(ParticleCache &pcache) const
{
//    return new ARSSC_Cache(pcache);
    return NULL;
}


// READ/WRITE/COPY.

// Returns a copy of the model data.
ARSSC_Model *const ARSSC_Model::Clone(void) const
{
    return new ARSSC_Model(*this);
}

// Returns the model data type.  Used to identify different models
// and for serialisation.
SubModelType ARSSC_Model::ID(void) const
{
    return ARSSC_Model_ID;
}

// Writes the object to a binary stream.
void ARSSC_Model::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Write number of sites.
        unsigned int n = (unsigned int)m_sites.size();
        out.write((char*)&n, sizeof(n));

        // Write site counts.
        double val = 0.0;
        for (unsigned int i=0; i!=n; ++i) {
            val = (double)m_sites[i];
            out.write((char*)&val, sizeof(val));
        }

        // NEXT VARIABLES ARE STATIC.  NEED SOME METHOD OF
        // OUTPUTTING THEM WITHOUT COPYING THEM FOR EACH
        // OUTPUTTED PRIMARY.

        // Write index of PAH tracker.
        unsigned int i = (unsigned int)m_ipah;
        out.write((char*)&i, sizeof(i));

        // Write neighbour weights.
        for (unsigned int i=0; i!=n; ++i) {
            // Free-edges.
            val = (double)m_fe_wt[i];
            out.write((char*)&val, sizeof(val));
            // Zigzags.
            val = (double)m_zz_wt[i];
            out.write((char*)&val, sizeof(val));
            // Armchairs.
            val = (double)m_ac_wt[i];
            out.write((char*)&val, sizeof(val));
            // Bays.
            val = (double)m_bay_wt[i];
            out.write((char*)&val, sizeof(val));
        }

        // Write combined-site functionals.
        Maths::Functional::Write(out, *m_r6_func);
        Maths::Functional::Write(out, *m_ac_r6_func);
        Maths::Functional::Write(out, *m_r5_ed_func);
        Maths::Functional::Write(out, *m_r5_ac_func);

        // Write combined-site independent variables.
        i = (unsigned int)m_r6_var;
        out.write((char*)&i, sizeof(i));
        i = (unsigned int)m_ac_r6_var;
        out.write((char*)&i, sizeof(i));
        i = (unsigned int)m_r5_ed_var;
        out.write((char*)&i, sizeof(i));
        i = (unsigned int)m_r5_ac_var;
        out.write((char*)&i, sizeof(i));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, ARSSC_Model::Serialize).");
    }
}

// Reads the object from a binary stream.
void ARSSC_Model::Deserialize(std::istream &in, Primary &parent)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0, i=0;
        double       val = 0.0;

        switch (version) {
            case 0:
                // Read number of sites.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Resize sites vectors.
                m_sites.resize(n, 0.0);
                m_fe_wt.resize(n, 0.0);
                m_zz_wt.resize(n, 0.0);
                m_ac_wt.resize(n, 0.0);
                m_bay_wt.resize(n, 0.0);

                // Read sites.
                for (i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_sites[i] = (real)val;
                }

                // Read PAH tracker index.
                in.read(reinterpret_cast<char*>(&i), sizeof(i));
                m_ipah = (unsigned int)i;

                // Read neighbour weights.
                for (i=0; i!=n; ++i) {
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_fe_wt[i] = (real)val;
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_zz_wt[i] = (real)val;
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_ac_wt[i] = (real)val;
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    m_bay_wt[i] = (real)val;
                }
                
                // Delete current combined-site functionals.
                delete m_r6_func;
                delete m_ac_r6_func;
                delete m_r5_ed_func;
                delete m_r5_ac_func;

                //  Read combined-site functionals.
                m_r6_func    = Maths::Functional::Read(in);
                m_ac_r6_func = Maths::Functional::Read(in);
                m_r5_ed_func = Maths::Functional::Read(in);
                m_r5_ac_func = Maths::Functional::Read(in);

                // Read combined-site independent variables.
                in.read(reinterpret_cast<char*>(&i), sizeof(i));
                m_r6_var = (SiteType)i;
                in.read(reinterpret_cast<char*>(&i), sizeof(i));
                m_ac_r6_var = (SiteType)i;
                in.read(reinterpret_cast<char*>(&i), sizeof(i));
                m_r5_ed_var = (SiteType)i;
                in.read(reinterpret_cast<char*>(&i), sizeof(i));
                m_r5_ac_var = (SiteType)i;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, ARSSC_Model::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, ARSSC_Model::Deserialize).");
    }
}
