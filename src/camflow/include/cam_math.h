/*!
 * File:   cam_help.h
 * Author: vj231
 *
 * Created on 13 February 2009, 11:35
 */

#ifndef _CAM_HELP_H
#define	_CAM_HELP_H

#include <vector>
#include "cam_error.h"
#include "cam_params.h"
#include "array.h"
#include <cmath>

namespace Camflow{

	/*!
	 *@brief    Math solver class.
	 *
	 * Include a more detailed description here.
	 */
    class CamMath{
    public:
        template <class T> T sum(const std::vector<T>& data);
        template <class T> T sum(const std::vector<T>& vec1, std::vector<double>& vec2);
        double sumVector(std::vector<double>& vec1);
        double sumVector(std::vector<double>& vec1, std::vector<double>& vec2);
        //double dydx(double nr1, double nr2, double dr);
        //lagrange interpolation
        double interpolateLG(double at, int size, Array2D &prime,
                                 const std::vector<double>& val);

        void binomCoeff(int n, Array2D &bCoeff);
        void prime(int size, Array2D &prime);
        void TDMA(std::vector<double>& a, std::vector<double>&b,
                  std::vector<double>&c, std::vector<double>& r,
                  std::vector<double>& u);

        double erf(double x);
        double erfc(double x);
        double erfccheb(double z);
        double inverfc(double p);
        double inverf(double p);

        template<class T>
        inline T SQR(const T a){
            return a*a;
        };

    private:
        static const int ncof=28;
        static const double cof[28];


    };
}

#endif	/* _CAM_HELP_H */

