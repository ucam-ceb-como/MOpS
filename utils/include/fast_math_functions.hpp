#ifndef fastMathFunctions_H
#define fastMathFunctions_H

// Namespace for fast math functions.
namespace fastMath {

    typedef double real;

    inline real pow2(const real& x) 
    { 
        return x*x; 
    } 

    inline real pow3(const real& x) 
    { 
        return x*pow2(x); 
    } 

    inline real pow4(const real& x) 
    { 
        return pow2(x)*pow2(x); 
    } 

    // iterative cube root approximation using Halley's method (double)
    inline double cbrta_halleyd(const double a, const double R)
    {
	    const double a3 = a*a*a;
        const double b= a * (a3 + R + R) / (a3 + a3 + R);
	    return b;
    }

    // cube root approximation using bit hack for 64-bit float 
    // adapted from Kahan's cbrt
    inline double cbrt_5d(double d)
    {
	    const unsigned int B1 = 715094163;
	    double t = 0.0;
	    unsigned int* pt = (unsigned int*) &t;
	    unsigned int* px = (unsigned int*) &d;
	    pt[1]=px[1]/3+B1;
	    return t;
    }

    // cube root approximation using 2 iterations of Halley's method (double)
    inline double halley_cbrt2d(double d)
    {
	    double a = cbrt_5d(d);
	    a = cbrta_halleyd(a, d);
	    return cbrta_halleyd(a, d);
    }


} // End namespace fastMathFunctions

#endif
