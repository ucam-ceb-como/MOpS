#include "cam_sdr.h"

using namespace Camflow;

ScalarDissipationRate::ScalarDissipationRate
(
    const int n_MixtureFractionPoints,
    const int n_TimePoints
)
:
  scalarDissipationRate_(n_MixtureFractionPoints, n_TimePoints, 0.0),
  scalarDissipationRateRef_(1.0)
{}

//! Destructor.
ScalarDissipationRate::~ScalarDissipationRate()
{}


const doublereal ScalarDissipationRate::calculate
(
    const doublereal mixtureFraction,
    const doublereal strainRate
)
{

    if (sdrType_ == ScalarDissipationRate::constant_fromStrainRate)
    {
        //scalarDissipationRateRef()
    }
    return 0.0;
}

/*!
 * Eq. 9.38 SummerSchool by N. Peters
 */
/*doublereal FlameLet::scalarDissipationRateRef(const doublereal m_frac)
{
    CamMath cm;
    doublereal erterm = cm.inverfc(2*m_frac);
    doublereal arg = -2*cm.SQR(erterm);
    scalarDissipationRateRef_ = admin_.getStrainRate()*exp(arg)/PI;
    return scalarDissipationRateRef_;
}*/
