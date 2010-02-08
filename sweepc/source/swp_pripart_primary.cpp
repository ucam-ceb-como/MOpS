/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweepc (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Implementation of the PriPartPrimary class declared in the
    swp_pripart_primary.h header file.

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

#include "swp_primary.h"
#include "swp_surfvol_primary.h"
#include "swp_pripart_primary.h"
#include "swp_submodel_type.h"
#include "swp_aggmodel_type.h"
#include "swp_model_factory.h"
#include "swp_surfvol_cache.h"
#include "rng.h"
#include <stdexcept>
#include <vector>

using namespace Sweep;
using namespace Sweep::AggModels;
using namespace Sweep::SubModels;
using namespace std;

// Default component is the first component.
unsigned int PriPartPrimary::m_icomp = 0;


// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor (protected).
PriPartPrimary::PriPartPrimary(void)
: m_totprisurf(0.0), m_avgpridiam(0.0)
{
}

// Initialising constructor.
PriPartPrimary::PriPartPrimary(real time, const Sweep::ParticleModel &model)
: SurfVolPrimary(time, model), m_totprisurf(0.0), m_avgpridiam(0.0)
{
}

// Copy constructor.
PriPartPrimary::PriPartPrimary(const PriPartPrimary &copy) 
{
    *this = copy;
}

// Stream-reading constructor.
PriPartPrimary::PriPartPrimary(std::istream &in, const Sweep::ParticleModel &model)
{
    Deserialize(in, model);
}

// Default destructor.
PriPartPrimary::~PriPartPrimary()
{
    // Nothing special to destruct.
}


// OPERATOR OVERLOADS.

// Assignment operator (Primary RHS).
PriPartPrimary &PriPartPrimary::operator=(const Primary &rhs)
{
    // Now check the type of the RHS to see if additional data
    // copying is required.
    if (rhs.AggID() == Spherical_ID) {
        // Use the SurfVolPrimary version of the operator in order
        // to set correctly the surface-volume parameters.
        SurfVolPrimary::operator=(rhs);
        
        // Additionally we have to set the pripart-list correctly with
        // a single primary particle.
        m_primaries.clear();
        vector<PriPart>::iterator i = m_primaries.insert(m_primaries.end(), PriPart());
        updatePriPartProperties(*i, (unsigned int)rhs.Composition(m_icomp));

        // Set bulk primary properties.
        m_totprisurf = i->Surface;
        m_avgpridiam = PriDiameter(0);
    } else {
        // Attempt to cast the primary to a SurfVolPrimary.  This will
        // throw an exception if the cast fails.  Note, the subsequent
        // SurfVolPrimary version of the operator= function will attempt
        // to cast to a PriPartPrimary, if it isn't a SurfVolPrimary.
        operator=(dynamic_cast<const SurfVolPrimary&>(rhs));
    }

    return *this;
}

// Assignment operator (SurfVolPrimary RHS).
PriPartPrimary &PriPartPrimary::operator=(const SurfVolPrimary &rhs)
{
    if (rhs.AggID() == SurfVol_ID) {
        // This IS a surface-volume primary.
        SurfVolPrimary::operator=(rhs);

        // Set the primary particle list to have approx. the
        // same monomer count and surface area as the RHS.
//        createPriList((unsigned int)rhs.Composition(m_icomp), rhs.SurfaceArea());

        // Add a single primary to the list with the correct mass.
        m_primaries.clear();
        vector<PriPart>::iterator i = m_primaries.insert(m_primaries.end(), PriPart());
        updatePriPartProperties(*i, (unsigned int)rhs.Composition(m_icomp));

        // Set bulk primary properties.
        m_totprisurf = i->Surface;
        m_avgpridiam = PriDiameter(0);
    } else {
        // Attempt to cast the primary to a PriPartPrimary.  This will
        // throw an exception if the cast fails.
        operator=(dynamic_cast<const PriPartPrimary&>(rhs));
    }

    return *this;
}

// Assignment operator (PriPartPrimary RHS).
PriPartPrimary &PriPartPrimary::operator=(const PriPartPrimary &rhs)
{
    // Assign the surface-volume properties.
    SurfVolPrimary::operator=(rhs);

    // Copy the pripart-list.
    m_primaries.assign(rhs.m_primaries.begin(), rhs.m_primaries.end());

    // Set bulk properties.
    m_totprisurf = rhs.m_totprisurf;
    m_avgpridiam = rhs.m_avgpridiam;

    return *this;
}

// Compound assignment (coagulation).
PriPartPrimary &PriPartPrimary::operator+=(const Primary &rhs)
{
    return PriPartPrimary::Coagulate(rhs);
}


// AGGREGATION MODEL.

// Returns the aggregation model which this primary describes.
AggModels::AggModelType PriPartPrimary::AggID(void) const {return AggModels::PriPartList_ID;}

// Creates an aggregation data cache for this primary type.
AggModels::PriPartCache *const PriPartPrimary::CreateAggCache(ParticleCache &pcache) const
{
    PriPartCache *cache = 
        static_cast<PriPartCache*>(ModelFactory::CreateAggCache(AggModels::PriPartList_ID, pcache));
    if (cache != NULL) *cache = *this;
    return cache;
}


// PRIMARY PARTICLE LIST.

// Returns the number of primary particles in the list.
unsigned int PriPartPrimary::PriCount(void) const {return m_primaries.size();}

// Returns the number of monomers in the ith primary.
unsigned int PriPartPrimary::MonomerCount(unsigned int i) const
{
    if (i < m_primaries.size()) {
        return m_primaries[i].Monomers;
    } else {
        return 0;
    }
}

// Returns the diameter of the ith primary.
real PriPartPrimary::PriDiameter(unsigned int i) const
{
    real vol = PriVolume(i);
    return pow(vol * 6.0 / PI, ONE_THIRD);
}

// Returns the average primary particle diameter.
real PriPartPrimary::AvgPriDiameter(void) const {return m_avgpridiam;}

// Returns the surface area of the ith primary.
real PriPartPrimary::PriSurface(unsigned int i) const
{
    if (i < m_primaries.size()) {
        return m_primaries[i].Surface;
    } else {
        return 0.0;
    }
}

// Returns the total surface area of the primaries in the list.
real PriPartPrimary::PriSurface(void) const {return m_totprisurf;}

// Returns the volume of the ith primary.
real PriPartPrimary::PriVolume(unsigned int i) const
{
    return PriMass(i) / m_pmodel->Components(m_icomp)->Density();
}

// Returns the mass of the ith primary.
real PriPartPrimary::PriMass(unsigned int i) const
{
    if (i < m_primaries.size()) {
        return (m_primaries[i].Monomers * m_pmodel->Components(m_icomp)->MolWt()) / NA;
    } else {
        return 0.0;
    }
}


// BASIC DERIVED PARTICLE PROPERTIES.

// Calculates the derived properties from the unique properties.  This
// function is broadly similar to the version in the spherical Primary
// class except that it uses the surface-volume model.  Therefore the
// surface area is not altered and the collision diameter is calculated
// using the arithmetic mean function.
void PriPartPrimary::UpdateCache(void)
{
    // Use the surface-volume model to update basic particle
    // properties.
    SurfVolPrimary::UpdateCache();

    // Now update the pripart-list.
    updatePrimaries();
}



// OPERATIONS.

// Adjusts the primary with the given composition and 
// tracker values changes n times.  If the particle cannot be adjust
// n times, then this function returns the number of times
// it was adjusted.
unsigned int PriPartPrimary::Adjust(const fvector &dcomp, const fvector &dvalues, 
                                    unsigned int n)
{
    // Use the surface-volume model to adjust basic
    // primary properties.
    SurfVolPrimary::Adjust(dcomp, dvalues, n);

    // Distribute the added monomers among the pripart-list.
    distMonomers((unsigned int)dcomp[m_icomp]*n);

    // Ensure the primaries are up-to-date.
    updatePrimaries();

    return n;
}

// Combines this primary with another.  This is also the
// implementation of the + and += operators.
PriPartPrimary &PriPartPrimary::Coagulate(const Primary &rhs)
{
    // Use surface-volume model to coagulate.
    SurfVolPrimary::Coagulate(rhs);

    // Now deal with the primary particle list.
    if (rhs.AggID() == Spherical_ID) {
        // Add a new pri-particle to the list.
        vector<PriPart>::iterator i = m_primaries.insert(m_primaries.end(), PriPart());
        updatePriPartProperties(*i, (unsigned int)rhs.Composition(m_icomp));
    } else if (rhs.AggID() == SurfVol_ID) {
        // Add a number of new pri-particles to the
        // list until the surface area is matched.
        vector<PriPart>::iterator i = m_primaries.insert(m_primaries.end(), PriPart());
        updatePriPartProperties(*i, (unsigned int)rhs.Composition(m_icomp));
//        addRandom((unsigned int)rhs.Composition(m_icomp), rhs.SurfaceArea());
    } else {
        // Attempt to cast the RHS as a PriPartPrimary.  This
        // will throw an exception if the cast fails.
        const PriPartPrimary &prirhs = dynamic_cast<const PriPartPrimary&>(rhs);

        // Merge the pri-particle lists.
        mergeInList(prirhs.m_primaries);
    }

    return *this;
}

// This routine sinters the Primary for the given length of
// time using the provided sintering model.
void PriPartPrimary::Sinter(real dt,const Cell &sys, 
                            const Processes::SinteringModel &model)
{
    // Use surface-volume model to adjust the primary m_surf
    // property.
    SurfVolPrimary::Sinter(dt, sys, model);

    // Begin Sintering Process (Currently, it is the same as SurfVolPrimary::Sinter).
    // Modification required for using primaries list infomation
    //// Perform a first order integration method to sinter
    //// the primary for the given time.
    //
    //// Declare time step variables.
    //real t1=0.0, delt=0.0, tstop=dt;
    //real r=0.0;

    //// Define the maximum allowed change in surface
    //// area in one internal time step (10% spherical surface).
    //real dAmax = 0.1 * m_sphsurf;

    //// The scale parameter discretises the delta-S when using
    //// the Poisson distribution.  This allows a smoother change
    //// (smaller scale = higher precision).
    //real scale = 0.01;

    //// Perform integration loop.
    //while (t1 < tstop) {
    //    // Calculate sintering rate.
    //    r = model.Rate(m_time+t1, sys, *this);

    //    // Calculate next time-step end point so that the
    //    // surface area changes by no more than dAmax.
    //    delt = dAmax / max(r, 1.0e-300);

    //    // Approximate sintering by a poisson process.  Calculate
    //    // number of poisson events.
    //    int n;
    //    if (tstop > (t1+delt)) {
    //        // A sub-step, we have changed surface by dAmax, on average.
    //        n = ignpoi(1.0 / scale);
    //    } else {
    //        // Step until end.  Calculate degree of sintering explicitly.
    //        n = ignpoi(r * (tstop - t1) / (scale*dAmax));
    //    }

    //    // Adjust the surface area.
    //    if (n > 0) {
    //        m_surf -= (real)n * scale * dAmax;
    //        // Check that primary is not completely sintered.
    //        if (m_surf <= m_sphsurf) {
    //            m_surf = m_sphsurf;
    //            break;
    //        }
    //    }

    //    // Set t1 for next time step.
    //    t1 += delt;
    //}

    // Now update the primary particle list to match the new 
    // surface area.
    updatePrimaries();
}


// READ/WRITE/COPY.

// Returns a copy of the model data.
PriPartPrimary *const PriPartPrimary::Clone(void) const
{
    return new PriPartPrimary(*this);
}

// Writes the object to a binary stream.
void PriPartPrimary::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Output the version ID (=0 at the moment).
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));

        // Output base class.
        SurfVolPrimary::Serialize(out);

        // Write the number of primary particles.
        unsigned int n = (unsigned int)m_primaries.size();
        out.write((char*)&n, sizeof(n));

        // Write the primary particle list.
        double val = 0.0;
        for (vector<PriPart>::const_iterator i=m_primaries.begin(); 
             i!=m_primaries.end(); ++i) {
            n   = i->Monomers;
            val = i->Surface;
            out.write((char*)&n, sizeof(n));
            out.write((char*)&val, sizeof(val));
        }

        // Write the total surface area.
        val = (double)m_totprisurf;
        out.write((char*)&val, sizeof(val));

        // Write the average particle diameter.
        val = (double)m_avgpridiam;
        out.write((char*)&val, sizeof(val));
    } else {
        throw invalid_argument("Output stream not ready "
                               "(Sweep, PriPartPrimary::Serialize).");
    }
}

// Reads the object from a binary stream.
void PriPartPrimary::Deserialize(std::istream &in, const Sweep::ParticleModel &model)
{
    if (in.good()) {
        // Read the output version.  Currently there is only one
        // output version, so we don't do anything with this variable.
        // Still needs to be read though.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0, num = 0;
        double val     = 0.0;

        switch (version) {
            case 0:
                // Read base class.
                SurfVolPrimary::Deserialize(in, model);

                // Read number of primary particles.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read primary particle list.
                for (unsigned int i=0; i!=n; ++i) {
                    // Read monomer count.
                    in.read(reinterpret_cast<char*>(&num), sizeof(num));
                    // Read surface area.
                    in.read(reinterpret_cast<char*>(&val), sizeof(val));
                    // Add pri-part to the list.
                    vector<PriPart>::iterator j = m_primaries.insert(m_primaries.end(), 
                                                                     PriPart());
                    j->Monomers = num;
                    j->Surface  = val;
                }

                // Read the total surface area.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_totprisurf = (real)val;

                // Read the average primary diameter.
                in.read(reinterpret_cast<char*>(&val), sizeof(val));
                m_avgpridiam = (real)val;

                break;
            default:
                throw runtime_error("Serialized version number is invalid "
                                    "(Sweep, PriPartPrimary::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready "
                               "(Sweep, PriPartPrimary::Deserialize).");
    }
}


// MEMORY MANAGEMENT.

// Release all memory associated with object.
void PriPartPrimary::releaseMem(void)
{
    m_primaries.clear();
}

// Initialisation routine.
void PriPartPrimary::init(void)
{
    releaseMem();
    m_totprisurf = 0.0;
    m_avgpridiam = 0.0;
}


// PRIPART-LIST MANAGEMENT.

// Calculates the average primary-particle diameter for
// the current list.  The primaries are assumed to be
// spherical.
real PriPartPrimary::calcAvgDiam()
{
    real sumd = 0.0;
    for (vector<PriPart>::const_iterator i=m_primaries.begin(); 
         i!=m_primaries.end(); ++i) {
        //real vol = (real)i->Monomers  * m_pmodel->Components(m_icomp)->MolWt() / 
        //           (NA * m_pmodel->Components(m_icomp)->Density());
        sumd += i->Diameter;
    }
    m_avgpridiam = sumd / (real)(m_primaries.size());
    return m_avgpridiam;
}

void PriPartPrimary::updatePriPartProperties(PriPart &pripart, unsigned int n) {
    pripart.Monomers = n;
    pripart.Volume = (real)n * m_pmodel->Components(m_icomp)->MolWt() / 
                     (NA * m_pmodel->Components(m_icomp)->Density());
    pripart.Surface = PI * pow(6.0*pripart.Volume/PI, TWO_THIRDS);
    pripart.Diameter = pow(pripart.Volume*6.0/PI, ONE_THIRD);
}

// Distributes a number of monomers among the pri-particles, weighted
// by their surface areas.
void PriPartPrimary::distMonomers(unsigned int n)
{
    // Big N = number of primaries.
    unsigned int N = m_primaries.size();

    // dm = n = Number of monomer to distribute.
    unsigned int dm = n;

    // Sum up the surface areas in order to calculate probabilities.
    // Note:  The surface areas are summed in reverse (from smallest particle),
    //        and each monomer uses a different sum.  This ensures that the 
    //        correct total surface area is used when calculating the probabilities.
    static vector<real> surfsums;
    surfsums.resize(N, 0.0);
    surfsums[N-1] = m_primaries[N-1].Surface;
    for (int i=N-2; i>=0; --i) {
        surfsums[i] = surfsums[i+1] + m_primaries[i].Surface;
    }

    // Now distribute mass weighted by surface area, starting with the largest
    // particle, until the N-1 particle.
    for (unsigned int i=0; i<N-1; ++i) {
        // Determine number of monomers to be added to this primary (m).
        unsigned int m = ignbin(dm, (float)(m_primaries[i].Surface/surfsums[i]));
        
        // If m is zero then don't update surface area and primary particle properties
        // This conditon is added to reduce computational time.
        if (m > 0) {
            // Calculate surface area.
            m_totprisurf -= m_primaries[i].Surface;
            // Calculate primary particle properties
            updatePriPartProperties(m_primaries[i], m_primaries[i].Monomers + m);
            m_totprisurf += m_primaries[i].Surface;
            // Update number of remaining monomers.
            if ((dm-=m) == 0) break;
        }
    }

    // Whatever mass is left shall be put on the smallest (last)
    // primary particle.
    if (dm > 0) {
        m_totprisurf -= m_primaries[N-1].Surface;
        updatePriPartProperties(m_primaries[N-1], m_primaries[N-1].Monomers + dm);
        m_totprisurf += m_primaries[N-1].Surface;
    }

    // Ensure the primary list is sorted.
    sortList(0, m_primaries.size()-1);

    // Calculate the average diameter.
    calcAvgDiam();
}

// Removes n monomers from the primary particle list.  Monomers
// are chosen by weighting the particles by their surface area, 
// and starting selection from the largest primary.
void PriPartPrimary::removeMonomers(unsigned int n)
{
    // Big N = number of primaries.
    unsigned int N = m_primaries.size();

    // dm = n = Number of monomer to remove.
    unsigned int dm = n;

    // Sum up the surface areas in order to calculate probabilities.
    // Note:  The surface areas are summed in reverse (from smallest particle),
    //        and each monomer uses a different sum.  This ensures that the 
    //        correct total surface area is used when calculating the probabilities.
    static vector<real> surfsums;
    surfsums.resize(N, 0.0);
    surfsums[N-1] = m_primaries[N-1].Surface;
    for (unsigned int i=N-2; i>=0; --i) {
        surfsums[i] = surfsums[i+1] + m_primaries[i].Surface;
    }

    // Now remove mass weighted by surface area, starting with the largest
    // particle.
    for (unsigned int i=0; i!=N; ++i) {
        // Determine number of monomers to be added to this primary (m).
        unsigned int m = ignbin(dm, (float)(m_primaries[i].Surface/surfsums[i]));
        // Remove monomers.
        if (m > 0) {
            // Calculate surface area.
            m_totprisurf -= m_primaries[i].Surface;
            // Calculate primary particle properties
            updatePriPartProperties(m_primaries[i], m_primaries[i].Monomers - m);
            m_totprisurf += m_primaries[i].Surface;
            // Update number of remaining monomers to remove.
            if ((dm-=m) == 0) break;
        }
    }

    // Ensure the primary list is sorted.
    sortList(0, m_primaries.size()-1);

    // Calculate the average primary diameter.
    calcAvgDiam();
}

// Updates the primary particle list to attempt to match the current
// primary surface area (given by surface-volume model).
void PriPartPrimary::updatePrimaries(void)
{
    // Check that the primaries have got the correct number
    // of monomers before adjusting.
    unsigned int mass = 0;
    for (vector<PriPart>::iterator i=m_primaries.begin(); i!=m_primaries.end(); ++i) {
        mass += i->Monomers;
    }
    if (mass < m_comp[m_icomp]) {
        if (m_primaries.size() > 0) {
            // If primaries exist already then distribute the added
            // mass.
            distMonomers((unsigned int)m_comp[m_icomp] - mass);
        } else {
            // Add the first primary.
            m_primaries.push_back(PriPart());
            updatePriPartProperties(m_primaries[0], (unsigned int)m_comp[m_icomp]);
            m_totprisurf = m_primaries[0].Surface;
            m_avgpridiam = calcAvgDiam();
        }
    } else if (mass > m_comp[m_icomp]) {
        // Can't remove monomers if there are no primaries.  This
        // should never happen anyway.
        if (m_primaries.size() > 0) {
            removeMonomers(mass - (unsigned int)m_comp[m_icomp]);
        }
    }

    real ds = 0.0;
    real surfremoved = 0.0;
    unsigned int m = 0;

    // Now the difficult bit, probably need to destroy primaries to 
    // reduce surface area.  The methodology is to remove the smallest
    // primaries, adding their mass to the larger, until the total surface
    // area is less than or equal to the aggregate surface area.
    while ((m_totprisurf > m_surf) && (m_primaries.size() > 1)) {
        // Calculate required change in surface area to match
        // surface-volume model.
        ds = m_totprisurf - m_surf;

        // Reset counters.
        m = 0;
        surfremoved = 0.0;

        // Get iterator to last (smallest) primary.
        vector<PriPart>::iterator i = (m_primaries.end()-1);
        
        // Calculate how many primaries must be removed to reduced the
        // surface area by ds.  Store total monomers in these primaries.
        while ((surfremoved < ds) && (m_primaries.size() > 1)) {
            // Store the mass of the smallest primary, then remove it.
            m += i->Monomers;             // Store monomer count.
            m_totprisurf -= i->Surface;   // Remove contribution to total surface area.
            surfremoved += i->Surface;    // Sum up removed surface area.
            m_primaries.erase(i);         // Erase the primary.
            i = (m_primaries.end()-1);    // Get pointer to new last (smallest) primary.
        }

        // Distribute this mass over the remaining primaries and update
        // the surface area which must be removed.
        distMonomers(m);
    }
    
    // These two functions are already sorted in any of above sections.
    // However, if there is any change to this function, you may need to
    // enable these two functions to esure that the list is sorted.
    // sortList(0, m_primaries.size()-1);

    // Calculate the average primary diameter.
    // calcAvgDiam();
}

// Merges a pri-particle list into the current list.
void PriPartPrimary::mergeInList(const std::vector<PriPart> &particle_list)
{
    vector<PriPart>::const_iterator j;
    unsigned int i;

    // Loop over list to be merged.
    for (j=particle_list.begin(); j!=particle_list.end(); ++j) {
        // Loop over this list and find the first primary smaller than 
        // that being inserted.
        for (i=0; i!=m_primaries.size(); ++i) {
            if (j->Monomers > m_primaries[i].Monomers) {
                m_primaries.insert(m_primaries.begin()+i, PriPart(*j));
                m_totprisurf += j->Surface;
                ++i;
                break;
            }
        }
        // If the primary being inserted is smaller than all primaries in
        // the list then add it at the end.
        if (i==m_primaries.size()) {
            m_primaries.insert(m_primaries.end(), PriPart(*j));
            m_totprisurf += j->Surface;
        }
    }

    // Calculate the average primary diameter.
    calcAvgDiam();
}

// Sorts the primary particle list in descending order of mass.
void PriPartPrimary::sortList(unsigned int i1, unsigned int i2)
{
    // Insertion sort - this sorting algorithm is very good for 
    // sorted or almosted list. Previous combination sorting algorithms
    // of quick sort and bubble sort take about 36% of CPU time.
    // This is reduced to 6% with insertion sort and total computational
    // time is reduced by 10%. An optimization to this sorting problem
    // is to make a primary particle list class to keep track if the list
    // should be sorted.
    unsigned int n = m_primaries.size();
    if (n <= 1) {
        // No need to sort if the array only has one element.
        return;
    } else {
        PriPart temp;
        for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = i; j > 0; j--) {
                if (m_primaries[j-1].Monomers < m_primaries[j].Monomers) {
                    // Store k in temp.
                    temp.Monomers = m_primaries[j].Monomers;
                    temp.Surface  = m_primaries[j].Surface;
                    temp.Volume  = m_primaries[j].Volume;
                    temp.Diameter  = m_primaries[j].Diameter;
                    // Put k+1 into k.
                    m_primaries[j].Monomers = m_primaries[j-1].Monomers;
                    m_primaries[j].Surface  = m_primaries[j-1].Surface;
                    m_primaries[j].Volume  = m_primaries[j-1].Volume;
                    m_primaries[j].Diameter  = m_primaries[j-1].Diameter;
                    // Put temp into k+1.
                    // Put temp into k+1.
                    m_primaries[j-1].Monomers = temp.Monomers;
                    m_primaries[j-1].Surface  = temp.Surface;
                    m_primaries[j-1].Volume  = temp.Volume;
                    m_primaries[j-1].Diameter  = temp.Diameter;
                } else {
                    break;
                }
            }
        }
    }

}
