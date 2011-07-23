#include "cam_sdr.h"
#include "camxml.h"
#include "stringFunctions.h"

using namespace Camflow;

ScalarDissipationRate::ScalarDissipationRate
(
    const std::string& inputFileName,
    const doublereal& stoichZ,
    const std::vector<doublereal>& mixFracCoords,
    const int n_TimePoints
)
:
  sdrType_(0),
  stoichZ_(stoichZ),
  mixFracCoords_(mixFracCoords),
  scalarDissipationRate_(mixFracCoords.size(), n_TimePoints, 0.0)
{

    readStrainRate(inputFileName);
    stoichSDR_ = scalarDissipationRate(stoichZ_);
    for (size_t i=0; i<mixFracCoords.size(); ++i)
    {
        scalarDissipationRate_(i,0) = calculate(mixFracCoords[i]);
    }

}

//! Destructor.
ScalarDissipationRate::~ScalarDissipationRate()
{}

void
ScalarDissipationRate::readStrainRate(const std::string& inputFileName)
{
    CamXML::Document doc;
    doc.Load(inputFileName);
    const CamXML::Element* root = doc.Root();
    CamXML::Element *subnode, *opNode;
    std::vector<CamXML::Element*> radiativeSpecies;
    std::vector<CamXML::Element*>::const_iterator p;

    opNode = root->GetFirstChild("op_condition");
    subnode = opNode->GetFirstChild("strain");
    strainRate_ = IO::from_string<double>(subnode->Data());
}

doublereal
ScalarDissipationRate::calculate(const doublereal& mixtureFraction)
const
{
    CamMath cm;
    doublereal fZ = exp(-2*cm.SQR(cm.inverfc(2*mixtureFraction)));
    doublereal fZst = exp(-2*cm.SQR(cm.inverfc(2*stoichZ_)));
    return scalarDissipationRate(stoichZ_) * (fZ/fZst);// * (phi/phist);
}

/*!
 * Eq. 9.38 SummerSchool by N. Peters
 */
doublereal
ScalarDissipationRate::scalarDissipationRate(const doublereal& mixtureFraction)
const
{
    CamMath cm;
    doublereal fZ = exp(-2*cm.SQR(cm.inverfc(2*mixtureFraction)));
    return strainRate_*fZ/Sprog::PI;
}

/*
 *calculate the scalar dissipation rate
 */
//doublereal ScalarDissipationRate::scalarDissipationRate(const doublereal m_frac)
//{
//    /*
//     *Eq. 9.38 SummerSchool by N. Peters
//     */
//    CamMath cm;
//    doublereal erterm = cm.inverfc(2*m_frac);
//    doublereal arg = -2*cm.SQR(erterm);
//    sdr = admin_.getStrainRate()*exp(arg)/PI;
//    return sdr;
//    //return 0.88;
//}
//
///*
// *calculate the scalar dissipation rate profile. Method 1 in Carbonell(2009).
// */
//doublereal ScalarDissipationRate::scalarDissipationRateProfile
//(
//    const doublereal m_frac,
//    const doublereal stoichSDR,
//    const int cell
//)
//{
//
//    CamMath cm;
//
//    doublereal fZ = exp(-2*cm.SQR(cm.inverfc(2*m_frac)));
//    doublereal fZst = exp(-2*cm.SQR(cm.inverfc(2*stoichZ)));
//
//    /*Utils::LinearInterpolator<doublereal, doublereal> rhoInterpolate(reacGeom_.getAxpos(),m_rho);
//    doublereal rhoStoich = rhoInterpolate.interpolate(stoichZ);
//
//    doublereal phi = 0.75 *
//                     (
//                         cm.SQR(std::sqrt(m_rho[0]/m_rho[cell])+1.0)
//                       / (2.0*std::sqrt(m_rho[0]/m_rho[cell])+1.0)
//                     );
//    doublereal phist = 0.75 *
//                       (
//                           cm.SQR(std::sqrt(m_rho[0]/rhoStoich)+1.0)
//                         / (2.0*std::sqrt(m_rho[0]/rhoStoich)+1.0)
//                       );*/
//
//    return stoichSDR * (fZ/fZst);// * (phi/phist);
//
//}
//
//
//
///*
// *set the scalar dissipation rate provided by the external
// *calling program
// */
//void ScalarDissipationRate::setExternalScalarDissipationRate(const doublereal sr)
//{
//    sdr_ext = sr;
//}
//
///**
// *  When the scalar dissipation rate has a time history
// *  use that during intergration
// */
//void ScalarDissipationRate::setExternalScalarDissipationRate
//(
//    const std::vector<doublereal>& time,
//    const std::vector<doublereal>& sdr,
//    const bool analytic
//)
//{
//
//    v_sdr = sdr;
//    v_time = time;
//
//    timeHistory = true;
//    sdrAnalytic = analytic;
//
//}
//
///**
// *  When the scalar dissipation rate has a time history
// *  and has a profile with mixture fraction from the CFD.
// */
//void ScalarDissipationRate::setExternalScalarDissipationRate
//(
//    const std::vector<doublereal>& time,
//    const std::vector< std::vector<doublereal> >& sdr,
//    const std::vector< std::vector<doublereal> >& Zcoords
//)
//{
//
//    profile_sdr = sdr;
//    v_time = time;
//    cfdMixFracCoords = Zcoords;
//
//    sdrProfile = true;
//
//}
//
///**
// *  Interpolate and return the scalar dissipation rate
// */
//doublereal
//ScalarDissipationRate::getSDR(const doublereal time)
//const
//{
//
//    Utils::LinearInterpolator<doublereal, doublereal> timeInterpolate(v_time, v_sdr);
//
//    return timeInterpolate.interpolate(time);
//
//}
//
///**
// *  Interpolate and return the scalar dissipation rate from a profile that varies through time.
// */
//doublereal
//ScalarDissipationRate::getSDRfromProfile
//(
//    const doublereal time,
//    const doublereal Z
//)
//const
//{
//
//    std::vector<doublereal> sdrTime, sdrInterpolated;
//    std::vector<doublereal> cfdMixFracCoordsTime, cfdMixFracCoordsInterpolated;
//
//    sdrInterpolated.clear();
//    cfdMixFracCoordsInterpolated.clear();
//
//    for (size_t i=0; i<cfdMixFracCoords[0].size(); ++i)
//    {
//
//        sdrTime.clear();
//        sdrTime.push_back(profile_sdr[0][i]);
//        sdrTime.push_back(profile_sdr[1][i]);
//
//        cfdMixFracCoordsTime.clear();
//        cfdMixFracCoordsTime.push_back(cfdMixFracCoords[0][i]);
//        cfdMixFracCoordsTime.push_back(cfdMixFracCoords[1][i]);
//
//        Utils::LinearInterpolator<doublereal, doublereal> timeInterpolate(v_time, sdrTime);
//        doublereal sdrInterpolatedTime = timeInterpolate.interpolate(time);
//
//        Utils::LinearInterpolator<doublereal, doublereal> time2Interpolate(v_time, cfdMixFracCoordsTime);
//        doublereal cfdMixFracCoordsInterpolatedTime = time2Interpolate.interpolate(time);
//
//        sdrInterpolated.push_back(sdrInterpolatedTime);
//        cfdMixFracCoordsInterpolated.push_back(cfdMixFracCoordsInterpolatedTime);
//    }
//
//    Utils::LinearInterpolator<doublereal, doublereal> spaceInterpolate(cfdMixFracCoordsInterpolated, sdrInterpolated);
//
//    return spaceInterpolate.interpolate(Z);
//
//}
