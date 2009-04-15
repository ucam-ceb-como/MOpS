/* 
 * File:   array.h
 * Author: vj231
 *
 * Created on 14 April 2009, 10:25
 */

#ifndef _ARRAY_H
#define	_ARRAY_H
#include "cam_params.h"
#include <vector>
namespace Camflow{
    /*
     * This is costructed in a Column major format
     */
    class Array2D{
    public:
        /*
         *Default constructor
         */
        Array2D(){
            nRows = 0;
            nCols = 0;
            aData.clear();
        }
        /*
         *create an mXn Array and initialize
         *with zeros
         */
        Array2D(int m, int n, doublereal v=0.0){
            nRows = n;
            nCols = m;
            aData.resize(m*n,v);
        }
        /*
         *copy constructor/
         *
         */
        Array2D(const Array2D& a){

            nRows = a.nRows;
            nCols = a.nCols;
            aData.resize(nRows*nCols);
            aData = a.aData;
        }
        /*
         *Assign
         */
        Array2D& operator=(const Array2D& a){
            nRows = a.nRows;
            nCols = a.nCols;
            aData.resize(nRows*nCols);
            aData = a.aData;
            return *this;
        }
        /*
         *resize the array
         */
        virtual ~Array2D(){}

        void resize(int m, int n, doublereal v=0.0){
            nRows = n;
            nCols = m;
            aData.resize(m*n,v);
        }

        /*
         *Element setting column major
         */
        doublereal& operator()(int i, int j){
            return value(i,j);
        }

        /*
         *Retrive the elements column major
         */
        doublereal operator()(int i, int j) const{
            return value(i,j);
        }

        doublereal& value(int i, int j){
            return aData[nRows*i + j];
        }


        doublereal value(int i, int j) const{
            return aData[nRows*i + j];
        }
        
    protected:
        vector<doublereal> aData;
        int nRows, nCols;
    };
}

#endif	/* _ARRAY_H */

