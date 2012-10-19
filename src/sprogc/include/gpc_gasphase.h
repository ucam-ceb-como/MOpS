/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The GasPhase class is a child of the Mixture class that also has thermodynamic
    functionality.  In addition to maintaining a description of the mixture composition
    it also allows thermodynamic properties of the gas to be calculated.  This is intended
    as a base class for other classes which will define the equation of state used
    to define the mixture, and hence how the properties should be calculated.

    In addition to the routines defined by the ThermoInterface class, the GasPhase
    class defines routines to return the thermodynamic properties of the mixture,
    as defined by the Mixture class.  These routines return a constant reference
    to a vector of properties, which are stored in the GasPhase class.  In this
    way, speed improvements can be gained if the mixture is unchanged between
    calls to calculate the properties.  If the mixture is changed then the 
    properties are recalculated the next time the relevant routines are called.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
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

#ifndef GPC_GASPHASE_H
#define GPC_GASPHASE_H

#include "gpc_params.h"
#include "gpc_thermo.h"
#include "gpc_mixture.h"
#include "gpc_species.h"

namespace Sprog
{
namespace Thermo
{

/*!
 *@brief Gas with thermodynamic properties, but no equation of state.
 *
 *  The GasPhase class is a child of the Mixture class that also has thermodynamic
 *  functionality.  In addition to maintaining a description of the mixture composition
 *  it also allows thermodynamic properties of the gas to be calculated.  This is intended
 *  as a base class for other classes which will define the equation of state used
 *  to define the mixture, and hence how the properties should be calculated.

 *  In addition to the routines defined by the ThermoInterface class, the GasPhase
 *  class defines routines to return the thermodynamic properties of the mixture,
 *  as defined by the Mixture class.  These routines return a constant reference
 *  to a vector of properties, which are stored in the GasPhase class.  In this
 *  way, speed improvements can be gained if the mixture is unchanged between
 *  calls to calculate the properties.  If the mixture is changed then the 
 *  properties are recalculated the next time the relevant routines are called.
 */
class GasPhase : public ThermoInterface, public Mixture
{
public:
    // Constructors.
    GasPhase(const SpeciesPtrVector &sp); // Default constructor.

    // Destructors.
    virtual ~GasPhase(void); // Default destructor.


    // EQUATION OF STATE.
    //   The temperature and density routines exist in the Mixture class
    //   and do not need to be redefined here.

    // Returns pressure in current units.
    virtual double Pressure(void) const = 0;

    // Sets the pressure of the gas phase.
    virtual void SetPressure(double p) = 0;


    // MIXTURE THERMODYNAMIC PROPERTIES.

    // Returns vector of internal energies of all species.
//    virtual const fvector *const Us() const = 0;
    void Us(fvector &U) const;

    // Returns vector of enthalpies of all species.
//    virtual const fvector *const Hs() const = 0;
    void Hs(fvector &H) const;

    // Returns vector of entropies of all species.
//    virtual const fvector *const Ss() const = 0;
    void Ss(fvector &S) const;

    // Returns vector of Gibbs free energies of all species.
//    virtual const fvector *const Gs() const = 0;
    void Gs(fvector &G) const;

   // Returns vector of constant pressure heat capacities 
    // of all species.
//    virtual const fvector *const Cps() const = 0;
    void Cps(fvector &Cp) const;

    // Returns vector of constant volume heat capacities 
    // of all species.
//    virtual const fvector *const Cvs() const = 0;
    void Cvs(fvector &Cv) const;


    // MIXTURE THERMODYNAMIC PROPERTIES (DIMENSIONLESS).

    // Returns vector of dimensionless internal energies of all species.
//    virtual const fvector *const Us_RT() const = 0;
    void Us_RT(fvector &U) const;

    // Returns vector of dimensionless enthalpies of all species.
//    virtual const fvector *const Hs_RT() const = 0;
    void Hs_RT(fvector &H) const;

    // Returns vector of dimensionless entropies of all species.
//    virtual const fvector *const Ss_R() const = 0;
    void Ss_R(fvector &S) const;

    // Returns vector of dimensionless Gibbs free energies of all species.
//    virtual const fvector *const Gs_RT() const = 0;
    void Gs_RT(fvector &G) const;

    // Returns vector of dimensionless constant pressure heat capacities 
    // of all species.
//    virtual const fvector *const Cps_R() const = 0;
    void Cps_R(fvector &Cp) const;

    // Returns vector of dimensionless constant volume heat capacities 
    // of all species.
//    virtual const fvector *const Cvs_R() const = 0;
    void Cvs_R(fvector &Cv) const;


    // BULK MIXTURE PROPERTIES.

    // Calculates the bulk internal energies in current units.
    double BulkU() const;

    // Calculates the bulk enthalpy in current units.
    double BulkH() const;

    // Calculates the bulk entropy in current units.
    double BulkS() const;

    // Calculates the bulk Gibbs free energies in current units.
    double BulkG() const;

    // Calculates the mean molar heat capacity at const. P.
    double BulkCp() const;

    // Calculates the mean molar heat capacity at const. V.
    double BulkCv() const;


protected:
    // As in the Mixture class, the default constructor is protected to prevent
    // mixtures being created with no knowledge of the defining species set.
    GasPhase(void);
	Sprog::Mechanism *m_mech;
    /*
    // Flag to tell GasPhase class whether or not it should
    // store pre-calculated values or not.  This option may be
    // useful if the calling code has a small memory requirement.
    bool m_precalc;

    // Vectors to store pre-calculated thermodynamic properties.
    fvector m_Us, m_Hs, m_Ss, m_Gs, m_Cps, m_Cvs;
    fvector m_Us_RT, m_Hs_RT, m_Ss_R, m_Gs_RT, m_Cps_R, m_Cvs_R;

    // Flags indicated whether or not the pre-calculated thermodynamic
    // properties are valid.
    bool m_Us_valid, m_Hs_valid, m_Ss_valid, m_Gs_valid,
         m_Cps_valid, m_Cvs_valid;
    bool m_Us_RT_valid, m_Hs_RT_valid, m_Ss_R_valid, m_Gs_RT_valid,
         m_Cps_R_valid, m_Cvs_R_valid;
    */
};
};
};

#endif
