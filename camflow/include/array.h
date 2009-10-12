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

        /*
         *return the size mxn
         */
        int size() const{
            return aData.size();
        }
        
    protected:
        vector<doublereal> aData;
        int nRows, nCols;
    };


    class Array1D{
    public:
        //Dfault constructor
        Array1D(){
            minIndex = 0;
            maxIndex = 0;
            aData.clear();
        }

        Array1D(int m, int n, doublereal v=0.0){
            int size = n-m;
            minIndex = m;
            maxIndex = n;
            aData.resize(size,v);
        }
        /*
         *copy constructor
         */
        Array1D(const Array1D &a){
            minIndex = a.minIndex;
            maxIndex = a.maxIndex;
            aData = a.aData;
        }

        Array1D& operator=(const Array1D& a){
            minIndex = a.minIndex;
            maxIndex = a.maxIndex;
            aData = a.aData;
            return *this;
        }

        void resize(int m, int n, doublereal v=0.0){
            int size = n-m;
            minIndex = m;
            maxIndex = n;
            aData.resize(size,v);
            
        }

        doublereal& operator ()(int i){
            return value(i);
        }

        doublereal operator ()(int i) const{
            return value(i);
        }
        
        doublereal value(int i) const{
            
            return aData[i-minIndex];
        }

        doublereal& value(int i){
            return aData[i-minIndex];
        }


    protected:
        vector<doublereal> aData;
        int minIndex, maxIndex;
    };
}

#endif	/* _ARRAY_H */

