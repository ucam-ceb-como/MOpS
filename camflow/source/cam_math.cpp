
#include "array.h"

#include "cam_math.h"
using namespace Camflow;

template <class T>
T CamMath::sum(const vector<T>& data){
    double ss =0;
    int len = data.size();
    for (int i = 0; i < len; i++) {
        ss += data[i];
    }
    return ss;

};


template <class T>
T CamMath::sum(const vector<T>& vec1, vector<double>& vec2){
    double ss =0;
    int l1 = vec1.size();
    int l2 = vec2.size();
    int size = (l1<l2) ? l1 : l2;
    for (int i = 0; i < size; i++) {
        ss += (vec1[i]*vec2[i]);
    }

    return ss;

};


double CamMath::sumVector(vector<double>& vec1){
    return sum<double>(vec1);
}

double CamMath::sumVector(vector<double>& vec1, vector<double>& vec2){
    return sum<double>(vec1,vec2);
}

doublereal CamMath::interpolateLG(doublereal at, int size, Array2D& dPrime,
                                      const vector<doublereal>& val){
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
void CamMath::TDMA(     vector<doublereal>& a,
                        vector<doublereal>& b,
                        vector<doublereal>& c,
                        vector<doublereal>& r,
                        vector<doublereal>& u){

    int n = a.size();
    vector<doublereal> gam(n);
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
