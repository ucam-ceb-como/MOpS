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

//inline double CamMath::dydx(double nr1, double nr2, double dr){
//    return ((nr2-nr1)/dr);
//}

