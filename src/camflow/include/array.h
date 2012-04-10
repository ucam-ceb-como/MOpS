/*!
 * \file   array.h
 * \author V. Janardhanan
 *
 * \brief Class for creating and storing arrays.
 *
 *  Copyright (C) 2009 Vinod Janardhanan.
 *

 Licence:
    This file is part of "camflow".

    brush is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Prof Markus Kraft
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

#ifndef _ARRAY_H
#define _ARRAY_H

#include "cam_params.h"
#include <vector>
#include <iomanip>

namespace Camflow{
    /*!
     *@brief    Class to create and store data for a 2D array.
     *
     * Include a more detailed description here.
     */
    class Array2D{

            friend std::ostream& operator<<(std::ostream& output, const Array2D& array)
            {
                output << "Array ("
                       << array.nRows << " rows x " << array.nCols
                       << " columns):"
                       <<  std::endl;
                for (size_t i=0; i<array.nRows; ++i)
                {
                    for (size_t j=0; j<array.nCols; ++j)
                    {
                       output << std::scientific << std::setw(14) << array(i,j);
                    }
                    output << std::endl;
                }
                output << std::endl;

                return output;
            };

        public:

            //! Default constructor.
            Array2D(){
                nRows = 0;
                nCols = 0;
                aData.clear();
            }

            //! Create an m X n Array and initialize with zeros.
            Array2D(int m, int n, doublereal v=0.0){
                nRows = m;
                nCols = n;
                aData.resize(m*n,v);
            }

            //! Copy constructor.
            Array2D(const Array2D& a){
                nRows = a.nRows;
                nCols = a.nCols;
                aData.resize(nRows*nCols);
                aData = a.aData;
            }

            //! Assignment operator overload.
            Array2D& operator=(const Array2D& a){
                nRows = a.nRows;
                nCols = a.nCols;
                aData.resize(nRows*nCols);
                aData = a.aData;
                return *this;
            }

            //! Destructor.
            virtual ~Array2D(){}

            //! Resize the array.
            void resize(int m, int n, doublereal v=0.0){
                nRows = m;
                nCols = n;
                aData.resize(m*n,v);
            }

            //! Element setting column major.
            doublereal& operator()(int i, int j){
                return value(i,j);
            }

            //! Retrive the elements column major.
            doublereal operator()(int i, int j) const{
                return value(i,j);
            }

            //! Access an element in the array.
            doublereal& value(int i, int j){
                return aData[nCols*i + j];
            }

            //! Access an element in the array.
            doublereal value(int i, int j) const{
                return aData[nCols*i + j];
            }

            //! Return the size m X n.
            int size() const {
                return aData.size();
            }

        protected:

            std::vector<doublereal> aData;
            size_t nRows, nCols;

    };

    /*!
     *@brief    Class to create and store data for a 1D array.
     *
     * Include a more detailed description here.
     */
    class Array1D{

        public:

            //! Default constructor
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

            std::vector<doublereal> aData;
            int minIndex, maxIndex;

    };

} // End Camflow namespace.

#endif  /* _ARRAY_H */
