/* 
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
    class CamMath{
    public:
        template <class T> T sum(const vector<T>& data);
        template <class T> T sum(const vector<T>& vec1, vector<double>& vec2);
        double sumVector(vector<double>& vec1);
        double sumVector(vector<double>& vec1, vector<double>& vec2);
        //double dydx(double nr1, double nr2, double dr);
        //lagrange interpolation
        doublereal interpolateLG(doublereal at, int size, Array2D &prime, const vector<doublereal>& val);

        void binomCoeff(int n, Array2D &bCoeff);
        void prime(int size, Array2D &prime);
        void TDMA(vector<doublereal>& a, vector<doublereal>&b, vector<doublereal>&c,
                vector<doublereal>& r, vector<doublereal>& u);

        doublereal erf(doublereal x);
        doublereal erfc(doublereal x);
        doublereal erfccheb(doublereal z);
        doublereal inverfc(doublereal p);
        doublereal inverf(doublereal p);

        template<class T>
        inline T SQR(const T a){
            return a*a;
        };

    private:
        static const int ncof=28;
        static const doublereal cof[28];
        

    };
}

#endif	/* _CAM_HELP_H */

