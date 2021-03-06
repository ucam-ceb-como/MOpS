#include "cam_sdr.h"
#include "camxml.h"
#include "stringFunctions.h"
#include "linear_interpolator.hpp"

using namespace Camflow;

ScalarDissipationRate::ScalarDissipationRate
(
    const std::string& inputFileName,
    const double& stoichZ,
    const std::vector<double>& mixFracCoords,
    const int n_TimePoints
)
:
  sdrType_(NONE),
  stoichZ_(stoichZ),
  mixFracCoords_(mixFracCoords),
  interpolatorZeroTime_(NULL),
  interpolatorNextTime_(NULL)
{

    scalarDissipationRate_.resize(n_TimePoints);
    for (size_t t=0; t<n_TimePoints; ++t)
    {
        scalarDissipationRate_[t].resize(mixFracCoords.size());
    }

    sdrType_ = notFromCFD;

    readStrainRate(inputFileName);

    for (size_t i=0; i<mixFracCoords.size(); ++i)
    {
        scalarDissipationRate_[0][i] = calculate(mixFracCoords[i]);
    }

    interpolator_ = new Utils::LinearInterpolator<double, double>
                        (
                            mixFracCoords_,
                            scalarDissipationRate_[0]
                        );

}

//! Destructor.
ScalarDissipationRate::~ScalarDissipationRate()
{
    if (interpolator_ != NULL) delete interpolator_;
    if (interpolatorZeroTime_ != NULL) delete interpolatorZeroTime_;
    if (interpolatorNextTime_ != NULL) delete interpolatorNextTime_;
}

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
    if (subnode != NULL)
    {
        strainRate_ = IO::from_string<double>(subnode->Data());
        stoichSDR_ = scalarDissipationRate(stoichZ_);
    }
    else
    {
        subnode = opNode->GetFirstChild("sdr");
        if (subnode != NULL)
        {
            stoichSDR_ = IO::from_string<double>(subnode->Data());
            strainRate_ = strainRate(stoichZ_);
        }
        else
        {
            throw std::runtime_error
            (
                "No strain rate or stoich SDR read. Specify either\n"
                " <sdr>VALUE</sdr> or "
                " <strain>VALUE</strain> in <op_condition>."
            );
        }
    }

    std::cout << "Strain Rate = " << strainRate_
              << " | Stoich SDR = " << stoichSDR_
              << " | Stoich Mix. Frac. = " << stoichZ_
              << std::endl;
}

/*
 * Calculate the scalar dissipation rate profile. Method 1 in Carbonell(2009).
 */
double
ScalarDissipationRate::calculate(const double& mixtureFraction)
const
{
    CamMath cm;
    double fZ = exp(-2*cm.SQR(cm.inverfc(2*mixtureFraction)));
    double fZst = exp(-2*cm.SQR(cm.inverfc(2*stoichZ_)));

    /*Utils::LinearInterpolator<double, double> rhoInterpolate(reacGeom_.getAxpos(),m_rho);
    double rhoStoich = rhoInterpolate.interpolate(stoichZ);

    double phi = 0.75 *
                     (
                         cm.SQR(std::sqrt(m_rho[0]/m_rho[cell])+1.0)
                       / (2.0*std::sqrt(m_rho[0]/m_rho[cell])+1.0)
                     );
    double phist = 0.75 *
                       (
                           cm.SQR(std::sqrt(m_rho[0]/rhoStoich)+1.0)
                         / (2.0*std::sqrt(m_rho[0]/rhoStoich)+1.0)
                       );*/

    return scalarDissipationRate(stoichZ_) * (fZ/fZst);// * (phi/phist);
}

/*!
 * Eq. 9.38 SummerSchool by N. Peters
 */
double
ScalarDissipationRate::scalarDissipationRate(const double& mixtureFraction)
const
{
    CamMath cm;
    double fZ = exp(-2*cm.SQR(cm.inverfc(2*mixtureFraction)));
    return strainRate_*fZ/Sprog::PI;
}

double
ScalarDissipationRate::strainRate(const double& mixtureFraction)
const
{
    CamMath cm;
    double fZ = exp(-2*cm.SQR(cm.inverfc(2*mixtureFraction)));
    return stoichSDR_*Sprog::PI/fZ;
}

void ScalarDissipationRate::setStrainRate(const double strainRate)
{

    sdrType_ = notFromCFD;
    strainRate_ = strainRate;
    stoichSDR_ = scalarDissipationRate(stoichZ_);
    for (size_t i=0; i<mixFracCoords_.size(); ++i)
    {
        scalarDissipationRate_[0][i] = calculate(mixFracCoords_[i]);
    }
    delete interpolator_;
    interpolator_ = new Utils::LinearInterpolator<double, double>
                        (
                            mixFracCoords_,
                            scalarDissipationRate_[0]
                        );
}

void ScalarDissipationRate::setSDRRate(const double sdr)
{

    sdrType_ = notFromCFD;
    stoichSDR_ = sdr;
    strainRate_ = strainRate(stoichZ_);
    for (size_t i=0; i<mixFracCoords_.size(); ++i)
    {
        scalarDissipationRate_[0][i] = calculate(mixFracCoords_[i]);
    }
    delete interpolator_;
    interpolator_ = new Utils::LinearInterpolator<double, double>
                        (
                            mixFracCoords_,
                            scalarDissipationRate_[0]
                        );
}


/*!
 *  When the scalar dissipation rate has a time history
 *  use that during integration. The spatial profile will
 *  be set assuming that the sdr values are the stoichiometric
 *  ones.
 */
void ScalarDissipationRate::setExternalScalarDissipationRate
(
    const std::vector<double>& time,
    const std::vector<double>& sdr
)
{

    sdrType_ = constant_fromCFD;
    v_stoichSDR = sdr;
    v_time = time;

    scalarDissipationRate_.resize(time.size());
    for (size_t t=0; t<time.size(); ++t)
    {
        scalarDissipationRate_[t].resize(mixFracCoords_.size());
    }

    // Now calculate the analytic profile of sdr
    for (size_t t=0; t<sdr.size(); ++t)
    {
        stoichSDR_ = v_stoichSDR[t];
        strainRate_ = strainRate(stoichZ_);
        for (size_t z=0; z<mixFracCoords_.size(); ++z)
        {
            scalarDissipationRate_[t][z] = calculate(mixFracCoords_[z]);
        }
    }

    interpolatorZeroTime_ = new Utils::LinearInterpolator<double, double>
                            (
                                mixFracCoords_,
                                scalarDissipationRate_[0]
                            );
    interpolatorNextTime_ = new Utils::LinearInterpolator<double, double>
                            (
                                mixFracCoords_,
                                scalarDissipationRate_[1]
                            );

}

///**
// *  When the scalar dissipation rate has a time history
// *  and has a profile with mixture fraction from the CFD.
// */
/*void ScalarDissipationRate::setExternalScalarDissipationRate
(
    const std::vector<double>& time,
    const std::vector< std::vector<double> >& sdr,
    const std::vector< std::vector<double> >& Zcoords
)
{

    sdrType_ =  profile_fromCFD;
    profile_sdr = sdr;
    v_time = time;
    //mixFracCoords_ = Zcoords;


}*/

/**
 *  Interpolate and return the scalar dissipation rate
 */
/*double
ScalarDissipationRate::getSDR(const double time)
const
{

    Utils::LinearInterpolator<double, double> timeInterpolate(v_time, v_sdr);

    return timeInterpolate.interpolate(time);

}*/

/*!
 *  Interpolate and return the scalar dissipation rate from a profile that varies through time.
 */
double
ScalarDissipationRate::operator()
(
    const double& Z,
    const double& time
)
const
{

    if (sdrType_ == notFromCFD)
    {
        return interpolator_->interpolate(Z);
    }
    else if (sdrType_ == constant_fromCFD)
    {
        std::vector<double> sdrAtTime_i;
        sdrAtTime_i.push_back(interpolatorZeroTime_->interpolate(Z));
        sdrAtTime_i.push_back(interpolatorNextTime_->interpolate(Z));

        Utils::LinearInterpolator<double, double> interpolateOverTime
        (v_time, sdrAtTime_i);

        return interpolateOverTime.interpolate(Z);
    }
    else
    {
        throw std::logic_error("sdrType_ is not set to a valid value.");
    }

}

