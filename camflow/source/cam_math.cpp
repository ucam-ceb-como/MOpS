
#include "array.h"

#include "cam_math.h"
using namespace Camflow;

template <class T>
T CamMath::sum(const std::vector<T>& data){
    double ss =0;
    int len = data.size();
    for (int i = 0; i < len; i++) {
        ss += data[i];
    }
    return ss;

};


template <class T>
T CamMath::sum(const std::vector<T>& vec1, std::vector<double>& vec2){
    double ss =0;
    int l1 = vec1.size();
    int l2 = vec2.size();
    int size = (l1<l2) ? l1 : l2;
    for (int i = 0; i < size; i++) {
        ss += (vec1[i]*vec2[i]);
    }

    return ss;

};


double CamMath::sumVector(std::vector<double>& vec1){
    return sum<double>(vec1);
}

double CamMath::sumVector(std::vector<double>& vec1, std::vector<double>& vec2){
    return sum<double>(vec1,vec2);
}

doublereal CamMath::interpolateLG(doublereal at, int size, Array2D& dPrime,
                                      const std::vector<doublereal>& val){
    doublereal prod = 1;
    for( int i=0; i< size; i++ )
        prod *= (at-i);

    doublereal retVal = 1;
    for( int i=0; i<size; i++){        
        doublereal expnt = prod/((at-i)*dPrime(size,i+1));
        retVal *= pow(val[i],expnt);

    }

    return retVal;
}

void CamMath::binomCoeff(int n, Array2D& bCoeff){
    bCoeff.resize(n+1,n+1);
    for(int i=0; i<=n; i++){
        bCoeff(0,i) = 1.0;
        bCoeff(i,i) = 1.0;
        for(int j=1; j<=i; j++){
            bCoeff(j,i) = bCoeff(j-1,i-1) + bCoeff(j,i-1);
            
        }
    }
}

void CamMath::prime(int size, Array2D& prime){
    
    
    prime.resize(size+1,size+1);
    
    for (int n = 1; n < size+1; n++) {
        for(int m=1;m<=n;m++){
            prime(n,m) =1.0;
            for(int l=1; l<=n;l++){
                if(l!=m) prime(n,m) *= (m-l);
            }
            
        }


    }

}

/*
 *This is a Tridiagonal Matrix Algorith implementation to
 *solve equations of the following form
 *
 * ---                         ----- --
 * | b0 c0 0 ....                  | | u0  |   |r0  |
 * | a1 b1 c1 0 ....               | | u1  |   |r1  |
 * | 0  a2 b2 c2 0 ....            | | .   |   |    |
 * |                               | | .   | = |    |
 * |                               | |     |   |    |
 * |             an-2   bn-2  cn-2 | |     |   |    |
 * |__                  an-1  bn-1_| | un-1|   |rn-1|
 *
 */
void CamMath::TDMA(     std::vector<doublereal>& a,
                        std::vector<doublereal>& b,
                        std::vector<doublereal>& c,
                        std::vector<doublereal>& r,
                        std::vector<doublereal>& u){

    int n = a.size();
    std::vector<doublereal> gam(n);
    doublereal bet;
    if(b[0] == 0.0){
        throw CamError("Error in TDMA\n");
    }
    u[0] = r[0]/(bet=b[0]);
    for(int j=1; j<n; j++){
        gam[j] = c[j-1]/bet;
        bet = b[j]-a[j]*gam[j];
        if(bet == 0.0) throw CamError("Error in TDMA\n");
        u[j] = (r[j]-a[j]*u[j-1])/bet;
    }
    //Back substitution
    for(int j=(n-2); j>=0;j--){
        u[j] -= gam[j+1]*u[j+1];
    }

}


const doublereal CamMath::cof[28] = {-1.3026537197817094, 6.4196979235649026e-1,
	1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
	3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
	-1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
	6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
	9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
	-1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};


doublereal CamMath::erfccheb(doublereal z){
    int j;
    doublereal t,ty,tmp,d=0.,dd=0.;
    if (z < 0.) throw CamError("erfccheb requires nonnegative argument");
    t = 2./(2.+z);
    ty = 4.*t - 2.;
    for (j=ncof-1;j>0;j--) {
        tmp = d;
        d = ty*d - dd + cof[j];
        dd = tmp;
    }
    return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
}

doublereal CamMath::inverfc(doublereal p){
    doublereal x,err,t,pp;
    if (p >= 2.0) return -100.;
    if (p <= 0.0) return 100.;
    pp = (p < 1.0)? p : 2. - p;
    t = sqrt(-2.*log(pp/2.));
    x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
    for (int j=0;j<2;j++) {
        err = erfc(x) - pp;
        x += err/(1.12837916709551257*exp(-SQR(x))-x*err);
    }
    return (p < 1.0? x : -x);
}



doublereal CamMath::erf(doublereal x){
    if (x >=0.) return 1.0 - erfccheb(x);
    else return erfccheb(-x) - 1.0;
}

doublereal CamMath::erfc(doublereal x){
    if (x >= 0.) return erfccheb(x);
    else return 2.0 - erfccheb(-x);
}



