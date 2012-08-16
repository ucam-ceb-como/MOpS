/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    The Mixture class is designed as a base class for all mixtures.  It maintains a group
    of species mole fractions, the temperature and the mixture density.

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

#ifndef GPC_MIXTURE_H
#define GPC_MIXTURE_H

#include "gpc_params.h"
#include "gpc_species.h"
#include "gpc_unit_systems.h"
#include "gpc_mixture_type.h"
#include <vector>
#include <iostream>

namespace Sprog
{
namespace Thermo
{

/*!
 * The Mixture class is designed as a base class for all mixtures.  It maintains a group of species mole fractions, the temperature and the mixture density.
 *
 */
class Mixture
{
public:
    // Constructors.  It makes no sense to have a mixture without knowledge
    // of the species list, therefore the default constructor is protected.
    Mixture(const SpeciesPtrVector &sp); // Default constructor (requires species list).
    Mixture(const Mixture &copy);        // Copy constructor.
    Mixture(                       // Stream-reading constructor.
        std::istream &in,          //   - Stream from which to read.
        const SpeciesPtrVector &sp //   - Species list.
        );

    // Destructors.
    virtual ~Mixture(void); // Default destructor.

    // Operator overloads.
    Mixture &operator=(const Mixture &mix); // Assignment operator.

    //! Set viscosity model
    void SetViscosityModel(Sprog::ViscosityModel vmodel)
        {m_vmodel = vmodel;}

    //! Get viscosity model
    Sprog::ViscosityModel GetViscosityModel(void) const
        {return m_vmodel;}

    // TEMPERATURE.

    // Returns temperature.
    real Temperature(void) const;

    // Set the temperature.
    void SetTemperature(const real T);


    // Returns the PAH formation rate.
    real PAHFormationRate() const;

    // Set the PAH formation rate.
    void SetPAHFormationRate(const real r);

    //! Return bulk gas velocity
    real Velocity() const {return m_data[velocityIndex()];}

    //! Set bulk gas velocity (species must be defined already)
    void SetVelocity(const real v) {m_data[velocityIndex()] = v;}



    //PARTICLE RELATED DATA

    //! Return mixture fraction diffusion coefficient
    real MixFracDiffCoeff() const {return m_data[MixFracDiffCoeffIndex()];}

    //! Set mixture fraction diffusion coefficient
    void SetMixFracDiffCoeff (const real d) {m_data[MixFracDiffCoeffIndex()]= d;}

    //! Return derivative of mixture fraction
    real GradientMixFrac() const {return m_data[GradientMixFracIndex()];}

    //! Set derivative of mixture fraction
    void SetGradientMixFrac (const real g) {m_data[GradientMixFracIndex()] = g;}

    //! Return second derivative of mixture fraction
    real LaplacianMixFrac() const {return m_data[LaplacianMixFracIndex()];}

    //! Set second derivative of mixture fraction
    void SetLaplacianMixFrac (const real l) {m_data[LaplacianMixFracIndex()] = l;}

    //! Return derivative of density times mixture fraction diffusion coefficient
    real GradientTemperature() const {return m_data[GradientTemperatureIndex()];}

    //! Set derivative of density times mixture fraction diffusion coefficient
    void SetGradientTemperature (const real g) {m_data[GradientTemperatureIndex()] = g;}

    //! Alpha factor (used in ABF model to scale certain reaction rates)
    real Alpha() const {return m_data[AlphaIndex()];}

    //! Set alpha factor (used in ABF model to scale certain reaction rates)
    void SetAlpha(const real a) {m_data[AlphaIndex()] = a;}


    // SPECIES CONCENTRATIONS/FRACTIONS.

    // Returns the mole fractions of all species.
    const fvector &MoleFractions() const;

    // Returns the molar concentrations of all species in current units.
    void GetConcs(fvector &concs) const;

    // Returns the mass fractions of all species.
    void GetMassFractions(fvector &fracs) const;

    // Returns the mole fraction of species i.
    real MoleFraction(unsigned int i) const;

    // Returns the molar concentration of species i.
    real MolarConc(unsigned int i) const;

    // Returns the mass fraction of species i.
    real MassFraction(unsigned int i) const;

    // Sets the mole fractions of all species.
    void SetFracs(const fvector &fracs);

    // Sets the mole fractions of all species from an array of values.
    void SetFracs(
        const real fracs[], // The array of mole fractions.
        int n               // The length of the array.
        );

    // Sets the molar concentrations of all species.
    void SetConcs(const fvector &concs);

    // Sets the mass fractions of all species.
    void SetMassFracs(const fvector &fracs);

    // Checks the vector of mole fractions for validity by settings all negative
    // values to zero, and by normalising the values so that they sum
    // to one.
    void Normalise();


    // MIXTURE DENSITY.

    // Returns molar density in current units.
    real Density(void) const;

    // Returns mass density in current units.
    real MassDensity(void) const;

    // Sets the molar density.
    void SetDensity(real dens);

    // Sets the mass density.
    void SetMassDensity(real dens);


    // SPECIES WHICH DEFINE THE MIXTURE.

    // Returns the species for which this mixture is defined.
    const SpeciesPtrVector *const Species() const;

    // Sets the species for which this mixture is defined.
    void SetSpecies(const SpeciesPtrVector &sp);


    // RAW DATA ACCESS.

    // Returns a pointer to the raw data (mole fractions, temperature
    // and density)  This function is provided in order to allow
    // numerical operations to be performed on the mixture, for example
    // integration of ODEs.
    real *const RawData();
    const real *const RawData() const;


    // READ/WRITE/COPY FUNCTIONS.

    // Creates a copy of the mixture object.
    virtual Mixture *const Clone() const;

    // Writes the element to a binary data stream.
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */)
    {
        ar & m_data & m_species & m_vmodel;
    }


    // Identifies the mixture type for serialisation.
    virtual Serial_MixtureType SerialType() const;

    //FUNCTIONS added by Vinod


    // return the molar enthalpy of all species as a vector J/mol.
    const std::vector<real> getMolarEnthalpy(real T);
        //return the molar enthalpy for th given mixture
        const std::vector<real> getMolarEnthalpy();

    // routine to return the avg moleculat weight in Kg/mol.
    real getAvgMolWt(fvector &massFrac) const;

    //! Mixture fraction weighted average molecular weight
    real getAvgMolWt() const;

    //! RMS molecule collision diameter
    real getMeanCollisionSection() const;

    //TRANSPORT RELATED

    // returns the mixture viscosity in Kg/m-s.
    real getViscosity() const;

    //! Checks the mixture has sufficient transport data.
    void checkForTransportData() const;

    // returns the mixture thermal conductivity in J/m-s-K    .
    real getThermalConductivity(real pre) const;
    // returns the mixture specific heat capacity in J/Kg K.
    real getSpecificHeatCapacity(real T);
        //return the specifi heat capacity for a given temp in J/Kg K.
        real getSpecificHeatCapacity();
        //return the molar specific heats
        const std::vector<real> getMolarSpecificHeat();
    // returns the vector of mixture diffusion coefficient in m^2/s.
    const std::vector<real> getMixtureDiffusionCoeff(const real pre)const;

    //! Index of temperature in m_data
    size_t temperatureIndex() const {return m_species->size();}

    //! Index of density in m_data
    size_t densityIndex() const {return m_species->size() + 1;}

    //! Index of velocity in m_data
    size_t velocityIndex() const {return m_species->size() + 2;}

    //! Index of PAH formation rate in m_data
    size_t PAHFormationIndex() const {return m_species->size() + 3;}

    //! Index of mixture fraction diffusion coefficient
    size_t MixFracDiffCoeffIndex() const {return m_species->size() + 4;}

    //!  Index of derivative of mixture fraction
    size_t GradientMixFracIndex() const {return m_species->size() + 5;}

    //!  Index of second derivative of mixture fraction
    size_t LaplacianMixFracIndex() const {return m_species->size() + 6;}

    //!  Index of derivative of density times mixture fraction diffusion coefficient
    size_t GradientTemperatureIndex() const {return m_species->size() + 7;}

    //!  Index of alpha (used in ABF model to scale certain reaction rates)
    size_t AlphaIndex() const {return m_species->size() + 8;}

    //! Number of items of data that are not species concentrations that are stored at the end of m_data
    static const size_t sNumNonSpeciesData = 9;

    // Writes the mixture to a binary data stream.
    virtual void Serialize(std::ostream &out) const;

    // Reads the mixture data from a binary data stream.
    virtual void Deserialize(std::istream &in);

protected:
    // The default constructor is protected to prevent mixture objects being
    // generated without knowledge of the defining species list.
    Mixture(void);

private:
    // The data vector contains, in order, the species mole fractions,
    // the mixture temperature and the mixture density.
    fvector m_data;

    // Vector of species for which this mixture is defined.
    const SpeciesPtrVector *m_species;

    //! The viscosity model used
    Sprog::ViscosityModel m_vmodel;

};
} //namespace Thermo
} //namespace Sprog

#endif
