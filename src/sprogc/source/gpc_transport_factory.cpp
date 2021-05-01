/*
 Author(s):      Vinod Janardhanan (vj231)
 Project:        sprog (gas-phase chemical kinetics).
 Sourceforge:    http://sourceforge.net/projects/mopssuite

 Copyright (C) 2008 Vinod M Janardhanan.

 File purpose:
 This file contains the definition of a structure for a chemical species.  File
 also contains typdefs and structure definitions related to Species objects.

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

#include "gpc_transport_factory.h"
#include "gpc_idealgas.h"
#include "gpc_params.h"
#include <stdexcept> 
#include <cmath>

using namespace Sprog;
using namespace Sprog::Thermo;
using namespace Sprog::Transport;

// Collision integral values for the evaluation of diffusion coeffcient
// Adapted from DETCHEM
double TransportFactory::Omega11Matrix[37 * 5] = {
    4.008e+000, -1.0380e+00, 5.9659e+00, -2.9977e+00, 4.9812e-01, 3.1300e+00,
    -4.6244e-01, 2.2622e+00, -8.4707e-01, 1.0961e-01, 2.6490e+00, -3.6366e-01,
    1.4543e+00, -4.7309e-01, 5.2330e-02, 2.3140e+00, -2.6379e-01, 1.0326e+00,
    -2.9782e-01, 2.7717e-02, 2.0660e+00, -1.9750e-01, 8.0301e-01, -2.1612e-01,
    1.7991e-02, 1.8770e+00, -2.5097e-01, 9.7452e-01, -4.5890e-01, 9.1504e-02,
    1.7290e+00, -1.1931e-01, 5.5652e-01, -1.4021e-01, 1.0704e-02, 1.6122e+00,
    -8.8902e-02, 4.6422e-01, -1.0763e-01, 6.7026e-03, 1.5170e+00, -6.6921e-02,
    3.9820e-01, -8.7569e-02, 4.7731e-03, 1.4400e+00, -5.6468e-02, 3.5399e-01,
    -7.6371e-02, 4.0889e-03, 1.3204e+00, -4.0079e-02, 2.8619e-01, -5.9986e-02,
    3.2544e-03, 1.2340e+00, -3.4407e-02, 2.4399e-01, -5.2429e-02, 3.4858e-03,
    1.1680e+00, -1.4465e-02, 1.7825e-01, -2.5006e-02, -7.2435e-04, 1.1166e+00,
    -1.3643e-02, 1.5968e-01, -2.5965e-02, 6.4166e-04, 1.0750e+00, -6.3727e-03,
    1.2820e-01, -1.4746e-02, -9.7281e-04, 1.0006e+00, -9.1095e-03, 1.0503e-01,
    -1.6931e-02, 6.9808e-04, 9.5000e-01, -4.6060e-03, 7.7117e-02, -1.0881e-02,
    3.7158e-04, 9.1310e-01, -1.9581e-03, 5.6961e-02, -5.2681e-03, -3.0468e-04,
    8.8450e-01, -2.9286e-03, 5.0573e-02, -7.3139e-03, 4.8588e-04, 8.4280e-01,
    -1.4158e-03, 3.3049e-02, -2.7314e-03, -1.2165e-04, 8.1300e-01, -1.3882e-03,
    2.4887e-02, -2.0834e-03, -2.3635e-05, 7.8980e-01, -1.9075e-05, 1.7423e-02,
    -4.1225e-04, -2.2042e-04, 7.7110e-01, 2.0888e-05, 1.3897e-02, -2.4013e-04,
    -1.7556e-04, 7.5550e-01, -4.5415e-05, 1.1903e-02, -6.0156e-04, -3.7943e-05,
    7.4220e-01, 5.4903e-04, 8.6391e-03, 2.7416e-04, -1.5004e-04, 7.2022e-01,
    -2.7033e-04, 7.7439e-03, -6.1756e-04, 4.3613e-05, 7.0250e-01, 8.3428e-04,
    3.8508e-03, 8.2584e-04, -2.0902e-04, 6.8776e-01, 3.9894e-05, 3.9815e-03,
    2.0106e-04, -8.6618e-05, 6.7510e-01, 3.8786e-05, 3.3483e-03, 5.6978e-05,
    -3.8795e-05, 6.6400e-01, 3.2925e-04, 2.3460e-03, 3.2880e-04, -9.2696e-05,
    6.4140e-01, -1.6252e-04, 2.1257e-03, -4.1127e-05, -1.6239e-05, 6.2350e-01,
    3.8876e-05, 1.4220e-03, -2.7050e-05, -1.7042e-06, 6.0882e-01, -1.2964e-04,
    1.5607e-03, -4.0253e-04, 8.8076e-05, 5.9640e-01, 2.6913e-04, 8.8982e-05,
    6.1738e-04, -1.4247e-04, 5.7630e-01, -1.8613e-04, 9.0887e-04, -2.7653e-04,
    7.1825e-05, 5.4150e-01, -5.7840e-05, 4.9071e-04, -2.3284e-04, 5.6760e-05,
    5.1800e-01, 1.5100e-04, 7.2383e-04, -5.6918e-04, 1.2030e-04 };

// Collision integral values for the evaluation of viscosities
// Adapted from DETCHEM
double TransportFactory::Omega22Matrix[37 * 5] = {
    4.1000e+00, -4.9400e-01, 5.1705e+00, -2.4579e+00, 3.8700e-01, 3.2630e+00,
    -4.8580e-01, 2.4444e+00, -9.0334e-01, 1.1276e-01, 2.8400e+00, -4.7159e-01,
    1.5619e+00, -4.5765e-01, 4.0936e-02, 2.5310e+00, -3.6736e-01, 1.0905e+00,
    -2.6097e-01, 1.3879e-02, 2.2840e+00, -2.8365e-01, 8.4186e-01, -1.8024e-01,
    5.8124e-03, 2.0840e+00, -2.1794e-01, 6.9371e-01, -1.4468e-01, 4.3636e-03,
    1.9220e+00, -2.0893e-01, 6.9031e-01, -1.8193e-01, 1.4363e-02, 1.7902e+00,
    -1.2932e-01, 5.1435e-01, -1.0742e-01, 3.8897e-03, 1.6820e+00, -1.0196e-01,
    4.5847e-01, -9.7322e-02, 4.1072e-03, 1.5930e+00, -8.5797e-02, 4.2360e-01,
    -9.6560e-02, 5.8448e-03, 1.4550e+00, -4.8607e-02, 3.2642e-01, -6.4005e-02,
    1.9493e-03, 1.3550e+00, -3.2972e-02, 2.7273e-01, -5.1800e-02, 1.5361e-03,
    1.2800e+00, -2.4970e-02, 2.3581e-01, -4.5617e-02, 1.8995e-03, 1.2220e+00,
    -1.3736e-02, 1.9366e-01, -3.2066e-02, 3.8064e-04, 1.1760e+00, -1.2000e-02,
    1.6956e-01, -2.7242e-02, 3.7995e-04, 1.0933e+00, -6.3451e-03, 1.2403e-01,
    -1.7897e-02, 1.2523e-04, 1.0390e+00, -2.1846e-03, 9.0757e-02, -1.0363e-02,
    -2.7932e-04, 9.9960e-01, -1.8453e-03, 7.2675e-02, -8.3026e-03, 1.1244e-05,
    9.6990e-01, -2.0468e-03, 5.9254e-02, -5.6424e-03, -2.8151e-04, 9.2680e-01,
    -1.5662e-03, 4.1725e-02, -3.9558e-03, -4.8331e-05, 8.9620e-01, -5.9689e-04,
    2.9341e-02, -1.8027e-03, -1.5170e-04, 8.7270e-01, 1.2046e-04, 2.1443e-02,
    -6.0037e-04, -2.1344e-04, 8.5380e-01, 2.3847e-04, 1.6539e-02, -7.9724e-05,
    -2.1286e-04, 8.3790e-01, -2.8887e-04, 1.4856e-02, -8.8946e-04, -3.7426e-06,
    8.2430e-01, 5.6032e-04, 1.0655e-02, 3.2629e-04, -1.8930e-04, 8.0180e-01,
    5.2094e-04, 7.2131e-03, 6.8005e-04, -2.1524e-04, 7.8360e-01, -1.0792e-04,
    6.8273e-03, -4.4591e-04, 5.5406e-05, 7.6830e-01, 5.7900e-04, 3.4869e-03,
    1.1370e-03, -2.8172e-04, 7.5520e-01, 1.0191e-04, 3.7993e-03, 8.7340e-05,
    -2.1863e-05, 7.4360e-01, 1.1106e-04, 3.1512e-03, 1.9598e-04, -6.2322e-05,
    7.1982e-01, 5.9318e-05, 2.1699e-03, 9.4738e-05, -2.8537e-05, 7.0100e-01,
    1.2696e-04, 1.3406e-03, 2.2449e-04, -5.1240e-05, 6.8545e-01, 5.2502e-04,
    -3.2893e-05, 9.4361e-04, -2.1258e-04, 6.7230e-01, 3.0299e-05, 1.4509e-03,
    -5.6073e-04, 1.4627e-04, 6.5100e-01, -3.5478e-05, 6.2081e-04, 4.7129e-05,
    -2.3793e-05, 6.1400e-01, -5.1551e-05, 1.6994e-03, -1.2009e-03, 2.2986e-04,
    5.8870e-01, 4.7810e-04, 3.2386e-03, -2.7109e-03, 5.3094e-04 };

// Reduced temperature values against with the values of Omega11 and Omega22 are plotted
double TransportFactory::TStar[37] = {
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
    2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0,
    20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0 };

const double MixtureTransport::oneByRootEight =
                0.353553390593273762200422181052424;

// 1/(8 * PI * EPSILON0)
const double MixtureTransport::deltaStarCoeff = 4493775895;

// (3.0/16.0) * sqrt(2 * NA * pow(kB,3)) / sqrt(PI)
const double MixtureTransport::binaryDiffCoeff = 5.955831556e-24;

// F(T) function for Z_rot calculation Eq. 33 SAND86-8246
double rotFunctionFT (double T, double LJdepth)
{
    double factor = LJdepth / (kB * T);
    return (1 + 0.5 * sqrt(PI) * PI * sqrt(factor) + (0.25 * fastMath::pow2(PI)
                    + 2) * factor + sqrt(PI) * PI * sqrt(factor) * factor);
}

// correction factor for polar non polar interaction
double getChi
(
    const int polar,
    const int nonPolar,
    const Sprog::Thermo::Mixture &mix
)
{
    // routine returns polar non-polar interaction corrrection factor
    // @polar : index of polar species
    // @nonPolar : index of nonPolar species

    const SpeciesPtrVector *spv = mix.Species();
    IO::Transport& tdPolar = (*spv)[polar]->getTransportData();
    IO::Transport& tdNonPolar = (*spv)[nonPolar]->getTransportData();

    double alphaNonPolar = tdNonPolar.getPolarizability(); // polarity
    double sigmaNonPolar = tdNonPolar.getCollisionDiameter(); // LJ- collision dia
    double sigmaPolar = tdPolar.getCollisionDiameter(); // LJ-collision dia
    double muPolar = tdPolar.getDipoleMoment(); // dipole moment
    double epsilonPolar = tdPolar.getPotentialWellDepth(); // LJ well depth
    double epsilonNonPolar = tdNonPolar.getPotentialWellDepth(); // LJ well depth


    double reducedAlpha = alphaNonPolar / (sigmaNonPolar * sigmaNonPolar
                    * sigmaNonPolar);
    double reducedMu = muPolar / (sqrt(epsilonPolar * (sigmaPolar * sigmaPolar
                    * sigmaPolar)));

    return (
                1.0 + 0.25 * reducedAlpha * reducedMu
              * sqrt(epsilonPolar / epsilonNonPolar)
           );
}

double TransportFactory::polyFitOmega (double delta, double *matrixPtr) const
{
    int i;
    double omega = *(matrixPtr + 4);
    for (i = 0; i < 4; i++)
    {
        matrixPtr--;
        omega = omega * delta + *(matrixPtr + 4);
    }

    return omega;
}

double TransportFactory::QuadInterPol (int imin, double *oPtr, double Tr) const
{
    // Quadratic interpolation routine
    double x12, m12, x23, m23, x2x, m2x;
    x12 = 0.5 * (TStar[imin] + TStar[imin - 1]);

    m12 = (*(oPtr + 1) - *(oPtr)) / (TStar[imin] - TStar[imin - 1]);

    x23 = 0.5 * (TStar[imin] + TStar[imin + 1]);
    m23 = (*(oPtr + 2) - *(oPtr + 1)) / (TStar[imin + 1] - TStar[imin]);

    x2x = 0.5 * (TStar[imin] + Tr);

    m2x = (m12 * (x23 - x2x) + m23 * (x2x - x12)) / (x23 - x12);

    return (*(oPtr + 1) + m2x * (Tr - TStar[imin]));
}

// routine for finding the  interval in which the reduced temperature lays
void TransportFactory::getTempInterval (const double rT, int &iMin, int &iMax) const
{
    int i;

    for (i = 0; i < 37; i++)
        if (rT < TStar[i])
            break;

    iMin = i - 1;
    if (iMin < 0)
        iMin = 0;
    iMax = iMin + 3;
    if (iMax > 36)
    {
        iMax = 36;
        iMin = iMax - 3;
    }

}

//! Calculate the reduced temperature T* = kB*T/epsilon
double TransportFactory::getReducedTemp (const double T, const Sprog::Species &sp) const
{

    return (kB * T / sp.getTransportData().getPotentialWellDepth());

}

double TransportFactory::getOmega11 (double rT, const double deltaStar) const
{

    int iMin, iMax;
    double Omega11;
    double o11[3];

    // find the temperarature interval
    getTempInterval(rT, iMin, iMax);

    for (size_t i = 0; i < 3; ++i)
    {
        int j = iMin + i - 1;
        if (deltaStar < 1.0e-02)
        {
            o11[i] = Omega11Matrix[5 * j];
        }
        else
        {
            o11[i] = polyFitOmega(deltaStar, &Omega11Matrix[5 * j]);
        }
    }

    // Temperature denpendancy
    if (rT < TStar[1])
        rT = TStar[1];
    if (rT > 500.e0)
        rT = 500.0;

    if (rT <= 100)
        Omega11 = QuadInterPol(iMin, o11, rT); // quadratic interpolation
    else
        Omega11 = 0.623e0 + rT * (-1.36e-3 + rT * (3.46e-6 + rT * (-3.43e-9)));

    return Omega11;

}

double TransportFactory::getOmega22 (const double T, const Sprog::Species &sp) const
{

    int i;
    int iMin, iMax;
    double Omega22;
    double o22[3];

    double rT = getReducedTemp(T, sp);

    //get reduced Dipole moment and check the dependancy
    double deltaStar = sp.getTransportData().getReducedDipoleMoment();

    // find the temperarature interval
    getTempInterval(rT, iMin, iMax);

    // Temperature denpendancy
    if (rT < TStar[0])
        rT = TStar[1];
    if (rT > 500.e0)
        rT = 500.0;

    if (rT <= 100)
    {
        for (i = 0; i < 3; i++)
        {
            int j = iMin + i - 1; // row index
            if (deltaStar < 1.0e-02)
            {
                o22[i] = Omega22Matrix[5 * j];

            }
            else
            {
                o22[i] = polyFitOmega(deltaStar, &Omega22Matrix[5 * j]);
            }
        }
        Omega22 = QuadInterPol(iMin, o22, rT); // quadratic interpolation
    }
    else
    {
        Omega22 = 0.703e0 + rT * (-1.46e-3 + rT * (3.57e-6 + rT * (-3.43e-9)));
    }

    return Omega22;

}

double PureSpeciesTransport::getViscosity (const double T, const Sprog::Species &sp) const
{

    double eta, omega22;
    omega22 = getOmega22(T, sp);
    ::IO::Transport& td = sp.getTransportData();

    eta = (5.0 / 16.0) * sqrt(sp.MolWt() * kB * T / NA)
          / ( sqrt(Sprog::PI)* fastMath::pow2(td.getCollisionDiameter()) * omega22 );

    return eta;

}

double PureSpeciesTransport::getSlefDiffusionCoeff (
    const double T,
    double p,
    const Sprog::Species &sp) const
{

    IO::Transport& td = sp.getTransportData();

    double omega11;
    double rT = getReducedTemp(T, sp);
    double deltaStar = td.getReducedDipoleMoment();
    omega11 = getOmega11(rT, deltaStar);

    // return selfDiff.
    return (3.0 / 8.0) * sqrt(PI * NA * fastMath::pow3(kB * T) / sp.MolWt())
                    / (p * PI * fastMath::pow2(td.getCollisionDiameter())
                                    * omega11);

}

double PureSpeciesTransport::getThermalConductivity
(
    const double T,
    const double p,
    const double cp,
    const Sprog::Species &sp
) const
{

    // pure species thermal conductivities are used only for the calculation of mixture conductivity
    int molIndex;
    double CvTrans, CvRot, CvVib, Cv;
    double fTrans, fRot, fVib, Zrot;
    double A, B, T0 = 298.0;
    double selfDiff, eta_k, rho;
    double lambda;
    IO::Transport& td = sp.getTransportData();
    molIndex = td.getMoleculeIndex();

    if (molIndex == 0)
    {

        CvTrans = 1.5 * R;
        fTrans = 2.5;
        lambda = getViscosity(T, sp) * fTrans * CvTrans / sp.MolWt();

    }
    else
    {

        rho = p * sp.MolWt() / (R * T);

        Cv = cp - R;
        if (molIndex == 1)
        {
            CvTrans = 1.5 * R;
            CvRot = R;
            CvVib = Cv - 2.5 * R;
        }
        else
        {
            CvTrans = 1.5 * R;
            CvRot = 1.5 * R;
            CvVib = Cv - 3 * R;
        }

        selfDiff = getSlefDiffusionCoeff(T, p, sp);
        eta_k = getViscosity(T, sp);

        fVib = rho * selfDiff / eta_k;

        A = 2.5 - fVib;
        Zrot = td.getRotRelaxationNumber()
              +(
                   rotFunctionFT(T0, td.getPotentialWellDepth())
                 / rotFunctionFT(T, td.getPotentialWellDepth())
               );
        B = Zrot + (2 / PI) * ((5 * CvRot / (3 * R)) + fVib);

        fTrans = 2.5 * (1 - (2 * CvRot * A / (PI * CvTrans * B)));
        fRot = fVib * (1 + (2 * A / (PI * B)));

        lambda = eta_k * (fTrans * CvTrans + fRot * CvRot + fVib * CvVib)
                        / sp.MolWt();

    }

    return lambda;

}

//----------------------------------------------------------------------------------------//
// The following function are for mixture transport -------------------------------------//


/*!
 * \f[
 *      \eta = \sum^{N_{sp}}_{k=1} \left( \frac{x_k \eta_k}{\displaystyle\sum^{N_{sp}}_{j=1} x_j \Phi_{kj}} \right)
 * \f]
 *
 * \f[
 *      \Phi_{kj} = \frac{1}{\sqrt{8}}
                    \left( 1+ \frac{M_{w,k}}{M_{w,j}} \right)^{-1/2}
                    \left( 1+\left( \frac{\eta_k}{\eta_j} \right)^{1/2} \left( \frac{M_{w,j}}{M_{w,k}} \right)^{1/4}\right)^2
 * \f]
 *
 */
double MixtureTransport::getViscosity
(
    const double T,
    const Sprog::Thermo::Mixture &mix
) const
{
    double eta = 0.0;

    if (mix.GetViscosityModel() == Sprog::iChapmanEnskog) {
        // Use the Chapman Enskog viscosity model
        // requires a 'tran.dat' with the relevant parameters for
        // each species present in the 'chem.inp'

        const SpeciesPtrVector *spv = mix.Species();
        const std::vector<double>& moleFrac = mix.MoleFractions();

        std::vector<double> nkVec;
        double xTimesEta, xTimesPhi;
        double m_kj, m_jk, eta_kj, phi_kj, nk;

        // Precompute the viscosity of each species before the for loop
        // to save time.
        for (size_t k = 0; k != spv->size(); ++k)
        {
            nkVec.push_back((*spv)[k]->getViscosity(T));
        }

        for (size_t k = 0; k != spv->size(); ++k)
        {

            nk = nkVec[k];
            xTimesEta = moleFrac[k] * nk;

            xTimesPhi = 0.0;

            for (size_t j = 0; j != spv->size(); ++j)
            {

                m_kj = (*spv)[k]->MolWt() / (*spv)[j]->MolWt();
                m_jk = 1.0 / m_kj;
                eta_kj = nk / nkVec[j];
                phi_kj = oneByRootEight
                         / sqrt((1 + m_kj))
                         * fastMath::pow2(1.0 + sqrt(eta_kj) * sqrt(sqrt(m_jk)) );
                xTimesPhi += moleFrac[j] * phi_kj;

            }

            eta += xTimesEta / xTimesPhi;
        }

    } else if (mix.GetViscosityModel() == Sprog::iAir) {
        // Use viscosity of air
        // Ref: Kazakov & Frenklach (1998), Combust. Flame 114, 484-501.
        eta = 1.458e-6 * T * sqrt(T) / (T + 110.4);

    } else if (mix.GetViscosityModel() == Sprog::iArgon) {
        // Use viscosity of argon
        // Ref: Perry's Chemical Engineers' Handbook, 8e
        eta = 9.2121e-7 * pow(T, 0.60529) / (1 + (83.24 / T));

    } else if (mix.GetViscosityModel() == Sprog::iHydrogen) {
        // Use viscosity of hydrogen
        // Ref: Perry's Chemical Engineers' Handbook, 8e
        eta = 1.797e-7 * pow(T, 0.685) / (1 - (0.59 / T) + (140/fastMath::pow2(T)));
    } else {
        throw std::runtime_error("Unrecognised viscosity model"
                "Sprog::MixtureTransport::getViscosity");
    }
    return eta;
}


/*!
 * \f[
 *      \lambda = \frac{1}{2}
 *                \left(
 *                      \sum^{N_{sp}}_{k=1} x_k \lambda_k
 *                      + \frac{1}{\displaystyle\sum^{N_{sp}}_{k=1} x_k/\lambda_k}
 *                \right)
 * \f]
 */
double MixtureTransport::getThermalConductivity
(
    const double T,
    double p,
    const Sprog::Thermo::Mixture &mix
) const
{

    double lambdaProd = 0;
    double lambdaFrac = 0;
    double tc = 0.0;

    std::vector<double> cp;
    const std::vector<double>& moleFrac = mix.MoleFractions();

    const SpeciesPtrVector *spv = mix.Species();

    IdealGas ig(*spv);
    ig.CalcCps(T, cp);

    for (size_t k = 0; k != spv->size(); ++k)
    {
        tc = (*spv)[k]->getThermalConductivity(T, p, cp[k]);
        lambdaProd += moleFrac[k] * tc;
        lambdaFrac += moleFrac[k] / tc;
    }

    return 0.5 * (lambdaProd + (1.0 / lambdaFrac));
}

double MixtureTransport::binaryDiffusionCoeff
(
    const int j,
    const int k,
    const double T,
    const double p,
    const Sprog::Thermo::Mixture &mix
) const
{

    const SpeciesPtrVector *spv = mix.Species();
    IO::Transport& td_j = (*spv)[j]->getTransportData();
    IO::Transport& td_k = (*spv)[k]->getTransportData();

    double Chi = 1.0;

    double m_j = (*spv)[j]->MolWt();
    double m_k = (*spv)[k]->MolWt();
    double onebym_jk = (m_j + m_k) / (m_j * m_k);

    // polar polar interaction
    double mu_j = td_j.getDipoleMoment();
    double mu_k = td_k.getDipoleMoment();

    double mu_jk = 0.0;

    double sigma_j = td_j.getCollisionDiameter();
    double sigma_k = td_k.getCollisionDiameter();

    double epsilon_j = td_j.getPotentialWellDepth();
    double epsilon_k = td_k.getPotentialWellDepth();

    if (mu_j > 1.0e-50)
    { // j is polar
        if (mu_k > 1.0e-50)
        {
            mu_jk = sqrt(mu_j * mu_k);
        }
        else
        { // k is non -polar
            Chi = getChi(j, k, mix);
        }
    }
    else
    { // j is non polar
        if (mu_k > 1.0e-50)
        { // k is polar
            Chi = getChi(k, j, mix);
        }
        else
        {
            mu_jk = sqrt(mu_j * mu_k);
        }
    }

    double epsilon_jk = Chi * Chi * sqrt(epsilon_j * epsilon_k);

    double sigma_jk = 0.5 * (sigma_j + sigma_k);
    if (Chi != 1.0) sigma_jk *= pow(Chi, -(1.0 / 6.0));

    double rT = kB * T / epsilon_jk;
    double deltaStar = deltaStarCoeff * mu_j * mu_k
                     / ( epsilon_jk * fastMath::pow3(sigma_jk) );
    double omega11 = getOmega11(rT, deltaStar);

    return binaryDiffCoeff * (sqrt(fastMath::pow3(T) * onebym_jk)
                    / (p * fastMath::pow2(sigma_jk) * omega11));

}

std::vector<double> MixtureTransport::getMixtureDiffusionCoeff
(
    const double T,
    const double p,
    const Sprog::Thermo::Mixture &mix
) const
{

    const SpeciesPtrVector *spv = mix.Species();
    const int size = spv->size();
    const std::vector<double>& moleFracs = mix.MoleFractions();

    double bDiff;
    std::vector<double> Dmix(size);

    for (int k = 0; k != size; ++k)
    {
        for (int j = k + 1; j != size; ++j)
        {
            bDiff = binaryDiffusionCoeff(j, k, T, p, mix);
            Dmix[k] += moleFracs[j] / bDiff;
            Dmix[j] += moleFracs[k] / bDiff;
        }
    }

    for (int k = 0; k != size; ++k)
    {
        if (Dmix[k] != 0)
            Dmix[k] = (1 - mix.MassFraction(k)) / Dmix[k];
        else
            Dmix[k] = binaryDiffusionCoeff(k, k, T, p, mix);
    }

    return Dmix;

}

