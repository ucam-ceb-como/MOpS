#include "cam_radiation.h"

using namespace Camflow;

const doublereal Radiation::AtmToPascal = 101325.0;

Radiation::Radiation(const int totalCells)
:
  absorption(3,0.0),
  radiation(totalCells,0.0)
{}

Radiation::~Radiation()
{}

/*!
    *Computes the Planck mean absorption coefficients
    *
    *@param[in]      Temperature      Temperature at a point of the grid, as obtained by the energy residual function.
    *@param[out]     Absorption       An output array containing the computed Planck absorption coefficients.
    *
    * The coefficients are produced for H2O, CO2, and CO, because these are the species that tend to produce the
    * greatest spectral radiative effect in mixtures commonly associated with combustion. Coefficients are known
    * for CH4 as well, but they can lead to inaccurate results, and therefore are not included here.  The computation
    * of these coefficients are based on data from the RADCAL model.  This data was for temperatures between 300K
    * and 2500K.
    */
void Radiation::PlanckAbsorption
(
    const doublereal temperature
)
{

    //This quantity is reused repeatedly in the equations below
    const doublereal beta = 1000/temperature;

    //This helps avoid repeated, costly calls to pow().  There is a modest tradeoff of speed vs. accuracy here.
    const doublereal beta2 = beta * beta;
    const doublereal beta3 = beta2 * beta;
    const doublereal beta4 = beta3 * beta;
    const doublereal beta5 = beta4 * beta;

    // This code computes the Planck mean absorption coefficient for H2O, using an equation in reference [1].  The specific
    // equation is based on fitted results from the RADCAL program. See reference [4].  Multiplication by AtmToPascal
    // converts the absorption coefficients from units of 1/m * 1/atm  to 1/m * 1/pa.

    // CHECK SECOND TERM COEFFICIENT!
    absorption[0] = (-0.23093 - 1.12390 * beta + 9.41530 * beta2 - 2.99880 * beta3 + 0.51382 * beta4
                    - 1.86840e-5 * beta5)/AtmToPascal;


    // This code computes the Planck mean absorption coefficient for CO2, using an equation in reference [1].  The specific
    // equation is based on fitted results from the RADCAL program. See reference [4]. Multiplication by AtmToPascal
    // converts the absorption coefficients from units of 1/m * 1/atm  to 1/m * 1/pa.


    absorption[1] = (18.7410 - 121.3100 * beta + 273.5000 * beta2 - 194.0500 * beta3 + 56.3100 * beta4
                    - 5.8169 * beta5)/AtmToPascal;


    // This code computes the Planck mean absorption coefficient for CO, using two equations in reference [1].  It uses a
    // conditional statement on temperature to determine which of the two equations to use. The specific  equation is based
    // on fitted results from the RADCAL program. See reference [4]. Multiplication by AtmToPascal converts the absorption
    // coefficients from units of 1/m * 1/atm  to 1/m * 1/pa.


    if (temperature <= 750)
    {
        absorption[2]  =  (4.7869 + temperature * (-0.06953 + temperature * (2.95775e-4 + temperature *
                          (-4.25732e-7 + temperature * 2.02894e-10))))/AtmToPascal;
    }
    else
    {
        absorption[2]  = (10.0900 + temperature * (-0.01183 + temperature * (4.7753e-6 + temperature *
                         (-5.87209e-10 + temperature * -2.5334e-14))))/AtmToPascal;
    }


    // After this function has been called the end result is a vector populated with absorption coefficients, to be used in
    // the function RadiativeLoss.
}

/*!
 * Computes  the radiative heat loss term for radiative heat dissipation model
 *
 *@param[in]    rho            Density of the mixture, as obtained by the energy residual function.
 *@param[in]    cp             Specific heat of the mixture, as obtained by the energy residual function.
 *@param[in]    Temperature    Temperature at a point of the grid, as obtained by the energy residual function.
 *@param[in]    SootVolFrac    The soot value fraction, as provided by the user or obtained from an external source.
 *@return                      A term equivalent to the total radiative heat loss.
 *
 * This function will be called from the energy residual code in Camflow's Flamelet class,
 * which will provide such parameters as temperature, density, specific heat, mass fractions
 * and soot volume fractions.
*/
void
Radiation::RadiativeLoss
(
    const int i,
    const doublereal temperature,
    const doublereal opPre,
    const doublereal soot_vol_frac,
    const doublereal mole_frac_H2O,
    const doublereal mole_frac_CO2,
    const doublereal mole_frac_CO
)
{

    // RadiativeLoss requires a background setting temperature.  It is usually assumed to be 300K, unless experimental conditions
    // suggest another temperature. Used 300^4 which is 81e8 in temperaturePowers below.
    //const doublereal BackgroundTemp = 300;

    doublereal temperature2 = temperature*temperature;
	doublereal temperature4 = temperature2*temperature2;
	doublereal temperature5 = temperature4*temperature;

    // The function calls PlanckAbsorption to produce the absorption coefficients.
    PlanckAbsorption(temperature);

    // Using a single equation in reference [1], the function computes a Radiative Heat Dissipation term.  It uses the Planck mean
    // absorption coefficients for H2O, CO2 abd CO, as well as the density of each of these species, the specific heat of each species,
    // the temperature, the partial pressure of each species, and the soot volume fraction.


    // Partial pressures of species k.  0 corresponds to H2O, 1 corresponds to CO2, 2 corresponds to CO. Operating pressure
    // must be expressed in Pascals in the camflow.xml file.

    doublereal partialPress[3];
    partialPress[0]  =  mole_frac_H2O * opPre;
    partialPress[1]  =  mole_frac_CO2 * opPre;
    partialPress[2]  =  mole_frac_CO * opPre;

    //The following is used repeatedly in the loop below
    const doublereal temperaturePowers = 4 * 5.669e-8 * (temperature4 - 81e8);

    doublereal spectralRadiation = 0;   //An intermediate result
    for (unsigned int j = 0; j < 3; ++j)
    {
        // 0 = H2O,  1 = CO2,  2 = CO
        //Total spectral radiative heat loss.
        spectralRadiation += partialPress[j] * absorption[j];
    }

    //Total spectral radiative heat loss + radiative heat loss due to soot
    radiation[i] = temperaturePowers * spectralRadiation + 3.3337e-4 * soot_vol_frac * temperature5;

}
