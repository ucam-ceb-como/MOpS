/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    This file declares the ARomatic Site - Site Counting (ARS-SC) model.  This
    model is fully described in the following publication:

    Celnik et al., Combust. Flame, in press, 2008.


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

#ifndef SWEEP_ARSSC_MODEL_H
#define SWEEP_ARSSC_MODEL_H

#include "swp_params.h"
#include "swp_submodel_type.h"
#include "swp_submodel.h"
#include "swp_maths_functional.h"
#include <vector>
#include <iostream>

namespace Sweep
{
// Forward declare Primary, ParticleCache classes.
class Primary;
class ParticleCache;

namespace SubModels
{
// Forward declare ARSSC_Cache class.
class ARSSC_Cache;

class ARSSC_Model : public SubModel
{
public:
    // Number of aromatic site types defined in this model.
    static const unsigned int SiteTypeCount = 9;

    // Enumeration of different aromatic site types.  These are
    // given explicit values so that they can be used as 
    // vector indices.
    enum SiteType {
        // Principal sites.
        FreeEdge=0, Zigzag=1, Armchair=2, Bay=3, R5=4,
        // Combined sites.
        R6=5, ACR6=6, R5ED=7, R5AC=8,
        // Error site ID.
        InvalidSite=-1
    };

    // Constructors.
    //   Default constructor is protected to prevent a sub-model
    //   being created without knowledge of the parent primary.
    ARSSC_Model(Primary &parent);          // Initialising constructor.
    ARSSC_Model(const ARSSC_Model &copy);  // Copy constructor.
    ARSSC_Model(          // Stream-reading constructor.
        std::istream &in, //  - Input stream.
        Primary &parent   //  - Parent primary particle.
        );

    // Destructor.
    ~ARSSC_Model(void);

    // Operators.
    ARSSC_Model &operator=(const ARSSC_Model &rhs);


    // SITE INFO.

    // Returns the site type enum value given the site ID as
    // a string.  It is quite forgiving, but will return -1 if
    // the site name is genuinely not known.
    static SiteType IdentifySite(const std::string &id);


    // SITE COUNTS.

    // Returns the count for the given site type.
    real SiteCount(SiteType type) const;

    // Returns the number of free-edges.
    real FreeEdgeCount(void) const;

    // Returns the number of armchairs.
    real ArmchairCount(void) const;

    // Returns the number of zigzags.
    real ZigzagCount(void) const;

    // Returns the number of bays.
    real BayCount(void) const;

    // Returns the number of R5s.
    real R5Count(void) const;

    // Sets the site count of the given site type.
    void SetSiteCount(SiteType type, real n);

    // Adds n to the site count of the given site type.
    void AddSites(SiteType type, real n);

    // Subtracts n from the site count of the given site type.
    void RemoveSites(SiteType type, real n);

    // Updates the site counts using the rules for neighbouring
    // sites.  A negative count for the number of sites to update
    // indicates site decrementing.
    void AdjustNeighbourSites(
        SiteType parent, // Parent site type (to get correct neighbour weights).
        int n            // Number of sites to increment (+ve) / decrement (-ve).
        );

    // Updates the site counts using the rules for neighbouring
    // sites.  The neighbour sites weights are provided as a vector,
    // rather than using the weights for a given parent site.  A negative
    // count for the number of sites to update indicates site decrementing.
    void AdjustNeighbourSites(
        fvector wts, // Custom site weights.
        int n        // Number of sites to increment (+ve) / decrement (-ve).
        );

    // PAH COUNT.

    // Sets the index of the PAH count tracker variable.
    static void SetPAH_Tracker(unsigned int i);

    // Returns the number of PAHs.
    real PAH_Count(void) const;

    // Sets the number of PAHs.
    void SetPAH_Count(real n);

    // Adds n to the PAH count.
    void AddPAHs(real n);

    // Subtracts n from the PAH count.
    void RemovePAHs(real n);


    // NEIGHBOURING SITE WEIGHTS.

    // Returns the neighbouring site weights vector for
    // free-edge sites.
    static const fvector &FE_Neighbours(void);

    // Returns the neighbouring site weights vector for
    // zigzag sites.
    static const fvector &ZZ_Neighbours(void);

    // Returns the neighbouring site weights vector for
    // armchair sites.
    static const fvector &AC_Neighbours(void);

    // Returns the neighbouring site weights vector for
    // bay sites.
    static const fvector &BayNeighbours(void);

    // Sets the neighbour site (type y) weight for the
    // given site type (type x).
    static void SetNeighbour(
        SiteType x, // Parent site type.
        SiteType y, // Neighbouring site type.
        real weight // Probability weight.
        );

    // Randomly selects a neighbour site for the given
    // parent using the correct weights.  If the optional
    // dec flag is set to true then a neighbour for the
    // decrementing process is selected (free-edges are
    // ignored).
    SiteType SelectNeighbour(
        SiteType parent, // Parent site type.
        bool dec=false // Choose for decrementing process?
        ) const;

    // Randomly selects a neighbour site using the given custom
    // weights.  If the optional
    // dec flag is set to true then a neighbour for the
    // decrementing process is selected (free-edges are
    // ignored).
    SiteType SelectNeighbour(
        const fvector &wts, // Custom neighbour weights.
        bool dec=false      // Choose for decrementing process?
        ) const;


    // COMBINED-SITE CORRELATION FUNCTIONS.

    // Returns the combined-site correlation functional
    // for R6 rings.
    static Maths::Functional &R6_Functional(void);

    // Returns the combined-site correlation functional
    // for AC next to R6 rings.
    static Maths::Functional &ACR6_Functional(void);

    // Returns the combined-site correlation functional
    // for R5 next to ED.
    static Maths::Functional &R5ED_Functional(void);

    // Returns the combined-site correlation functional
    // for R5 next to AC.
    static Maths::Functional &R5AC_Functional(void);

    // Sets the combined-site correlation functional
    // for R6 rings.  The given functional is stored
    // by reference, so it would be unwise to delete
    // the reference after calling this routine.
    static void SetR6_Functional(Maths::Functional &func);

    // Sets the combined-site correlation functional
    // for AC next to R6 rings.  The given functional is stored
    // by reference, so it would be unwise to delete
    // the reference after calling this routine.
    static void SetACR6_Functional(Maths::Functional &func);

    // Sets the combined-site correlation functional
    // for R5 next to ED.  The given functional is stored
    // by reference, so it would be unwise to delete
    // the reference after calling this routine.
    static void SetR5ED_Functional(Maths::Functional &func);

    // Sets the combined-site correlation functional
    // for R5 next to AC.  The given functional is stored
    // by reference, so it would be unwise to delete
    // the reference after calling this routine.
    static void SetR5AC_Functional(Maths::Functional &func);


    // COMBINED-SITE CORRELATION INDEPENDENT VARIABLES.

    // Returns the independent variable (site type) used to 
    // calculate R6 ring counts.
    static SiteType R6_Variable(void);

    // Returns the independent variable (site type) used to 
    // calculate AC next to R6 ring counts.
    static SiteType ACR6_Variable(void);

    // Returns the independent variable (site type) used to 
    // calculate R5 next to ED counts.
    static SiteType R5ED_Variable(void);

    // Returns the independent variable (site type) used to 
    // calculate R5 next to AC counts.
    static SiteType R5AC_Variable(void);

    // Sets the independent variable (site type) used to
    // calculate R6 ring counts.
    static void SetR6_Variable(SiteType type);

    // Sets the independent variable (site type) used to
    // calculate AC next to R6 ring counts.
    static void SetACR6_Variable(SiteType type);

    // Sets the independent variable (site type) used to
    // calculate R5 next to ED counts.
    static void SetR5ED_Variable(SiteType type);

    // Sets the independent variable (site type) used to
    // calculate R5 next to AC counts.
    static void SetR5AC_Variable(SiteType type);


    // PARTICLE UPDATES.

    // Updates the model data according to the rules of the model
    // given the changes to particle composition and values.
    // Performs the update n times.
    void Update(
        const fvector &dcomp, // Changes to the composition.
        const fvector &dval,  // Changes to the tracker variables.
        unsigned int n=1      // Number of times to perform update.
        );

    // Updates the data according to the model rules for
    // coagulation.
    ARSSC_Model &Coagulate(const SubModel &rhs);


    // CACHE UPDATES.

    // Updates the sub-model cached variables.
    void UpdateCache(void);

    // Creates a model cache which subscribes to this particle model.
    SubModelCache *const CreateCache(ParticleCache &pcache) const;


    // READ/WRITE/COPY.

    // Returns a copy of the model data.
    ARSSC_Model *const Clone(void) const;

    // Returns the model data type.  Used to identify different models
    // and for serialisation.
    SubModelType ID(void) const;

    // Writes the object to a binary stream.
    void Serialize(std::ostream &out) const;

    // Reads the object from a binary stream.
    void Deserialize(
        std::istream &in, // Input stream.
        Primary &parent   // Parent object.
        );

private:
    // Principal and combined site counts.
    fvector m_sites;

    // Index of the PAH count tracker variable in the
    // parent primary.
    static unsigned int m_ipah;

    // Neighbouring site probability weights.  These are static
    // because they are common to all instances of the ARSSC model
    // data.  Weights are stored in order Edge/Zigzag/Armchair/Bay.
    static fvector m_fe_wt;  // Free-edge neighbour weights.
    static fvector m_zz_wt;  // Zig-zag neighbour weights.
    static fvector m_ac_wt;  // Armchair neighbour weights.
    static fvector m_bay_wt; // Bay neighbour weights.

    // Combined-site correlation functions.  There are currently
    // four hardcoded combined sites.
    static Maths::Functional *m_r6_func;    // R6 rings functional.
    static Maths::Functional *m_ac_r6_func; // AC next to R6 functional.
    static Maths::Functional *m_r5_ed_func; // R5 next to ED functional.
    static Maths::Functional *m_r5_ac_func; // R5 next to AC functional.

    // Combine-site correlation independent variables (site types).
    static SiteType m_r6_var;    // Site type used as independent variable for R6 rings.
    static SiteType m_ac_r6_var; // Site type used as independent variable for AC-R6.
    static SiteType m_r5_ed_var; // Site type used as independent variable for R5-ED.
    static SiteType m_r5_ac_var; // Site type used as independent variable for R5-AC.

    // Sub-models cannot be created without knowledge of the parent
    // primary particle.
    ARSSC_Model(void);
};
};
};

#endif
