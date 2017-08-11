/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Defines coordinate structures and transform functions.

  Licence:
    This file is part of "sweepc".

    sweepc is free software; you can redistribute it and/or
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

#ifndef SWEEP_COORDS_H
#define SWEEP_COORDS_H

#include "swp_params.h"
#include <cmath>
#include <vector>
#include <memory.h>

namespace Sweep
{
namespace Coords
{
    // A 3D cartesian position vector.
    class Vector
    {
    private:
        double A[3]; // Vector data.
    public:
        // Constructor.
        Vector() {memset(&A[0], 0, sizeof(double)*3);}
        // Cartesian coordinates.
        inline double &X(void) {return A[0];}
        inline double &Y(void) {return A[1];}
        inline double &Z(void) {return A[2];}
        // Operators to return vector elements.
        double &operator[](int i) {return A[i];}
        const double &operator[](int i) const {return A[i];}

        // VECTOR TRANSLATION.

        // Translates the vector the given deviations.
        inline void Translate(double dx, double dy, double dz) {A[0]+=dx; A[1]+=dy; A[2]+=dz;}
        inline void Translate(Vector D) {A[0]+=D[0]; A[1]+=D[1]; A[2]+=D[2];}

    };

    // A 3D coordinate transform matrix.
    class Matrix 
    {
    private:
        double A[3][3]; // The matrix data [row][col].
    public:
        // Constructor.
        Matrix() {memset(&A[0][0], 0, sizeof(double)*9);}

        // OPERATORS.

        // Operators to return matrix elements.
        double* operator[](int i) {return &A[i][0];}
        const double *const operator[](int i) const {return &A[i][0];}

        // MATRIX MULTIPLICATION.

        // Performs matrix multiplication.  Returns
        // the resultant matrix.
        inline Matrix Mult(const Matrix &B) const
        {
            // C = A X B.
            Matrix C;
            C[0][0] = (A[0][0]*B[0][0]) + (A[0][1]*B[1][0]) + (A[0][2]*B[2][0]);
            C[0][1] = (A[0][0]*B[0][1]) + (A[0][1]*B[1][1]) + (A[0][2]*B[2][1]);
            C[0][2] = (A[0][0]*B[0][2]) + (A[0][1]*B[1][2]) + (A[0][2]*B[2][2]);
            C[1][0] = (A[1][0]*B[0][0]) + (A[1][1]*B[1][0]) + (A[1][2]*B[2][0]);
            C[1][1] = (A[1][0]*B[0][1]) + (A[1][1]*B[1][1]) + (A[1][2]*B[2][1]);
            C[1][2] = (A[1][0]*B[0][2]) + (A[1][1]*B[1][2]) + (A[1][2]*B[2][2]);
            C[2][0] = (A[2][0]*B[0][0]) + (A[2][1]*B[1][0]) + (A[2][2]*B[2][0]);
            C[2][1] = (A[2][0]*B[0][1]) + (A[2][1]*B[1][1]) + (A[2][2]*B[2][1]);
            C[2][2] = (A[2][0]*B[0][2]) + (A[2][1]*B[1][2]) + (A[2][2]*B[2][2]);
            return C;
        }

        // Performs matrix multiplication on a vector.
        inline Vector Mult(const Vector &B) const
        {
            // C = A X B.
            Vector C;
            C[0] = (A[0][0]*B[0]) + (A[0][1]*B[1]) + (A[0][2]*B[2]);
            C[1] = (A[1][0]*B[0]) + (A[1][1]*B[1]) + (A[1][2]*B[2]);
            C[2] = (A[2][0]*B[0]) + (A[2][1]*B[1]) + (A[2][2]*B[2]);
            return C;
        }
        
        // IDENTITY MATRIX.

        // Returns the identity matrix.
        static inline Matrix Identity(void)
        {
            Matrix I;
            I[0][0] = I[1][1] = I[2][2] = 1.0;
            return I;
        }

        // Sets this matrix to the identity matrix.
        inline void SetIdentity(void)
        {
            A[0][0] = A[1][1] = A[2][2] = 1.0;
            A[0][1] = A[0][2] = A[1][0] = 0.0;
            A[1][2] = A[2][0] = A[2][1] = 0.0;
        }

        // COORDINATE ROTATION.

        // Adds a rotation to the transform matrix around the Z-axis.
        inline void RotateZ(double phi)
        {
            // Precalculate cos and sin.
            double cosp = cos(phi);
            double sinp = sin(phi);

            // This function requires one temporary storage variable (a0)
            // to account for a change in the left column before it is used.

            // Top row.
            double a0 = A[0][0];
            A[0][0] = + (a0*cosp) + (A[0][1]*sinp);
            A[0][1] = - (a0*sinp) + (A[0][1]*cosp);
            //A[0][2] = A[0][2];
            // Middle row.
            a0 = A[1][0];
            A[1][0] = + (a0*cosp) + (A[1][1]*sinp);
            A[1][1] = - (a0*sinp) + (A[1][1]*cosp);
            //A[1][2] = A[1][2];
            // Bottom row.
            a0 = A[2][0];
            A[2][0] = + (a0*cosp) + (A[2][1]*sinp);
            A[2][1] = - (a0*sinp) + (A[2][1]*cosp);
            //A[2][2] = A[2][2];
        }

        // Sets the matrix to be a transform matrix for rotation around
        // the z-axis.
        inline void SetRotZ(double phi)
        {
            A[0][2] =   A[1][2] = A[2][0] = A[2][1] = 0.0;
            A[0][0] =   (A[1][1] = cos(phi));
            A[0][1] = - (A[1][0] = sin(phi));
            A[2][2] =   1.0;
        }

        // Adds a rotation to the transform matrix around the x-axis
        inline void RotateX(double theta)
        {
            // Precalculate cos and sin.
            double cosp = cos(theta);
            double sinp = sin(theta);

            // This function requires one temporary storage variable (a1)
            // to account for a change in the middle column before it is used.

            // Top row.
            double a1 = A[0][1];
            //A[0][0] = A[0][0];
            A[0][1] = + (a1*cosp) + (A[0][2]*sinp);
            A[0][2] = - (a1*sinp) + (A[0][2]*cosp);
            // Middle row.
            a1 = A[1][1];
            //A[1][0] = A[1][0];
            A[1][1] = + (a1*cosp) + (A[1][2]*sinp);
            A[1][2] = - (a1*sinp) + (A[1][2]*cosp);
            // Bottom row.
            a1 = A[2][1];
            //A[2][0] = A[2][0];
            A[2][1] = + (a1*cosp) + (A[2][2]*sinp);
            A[2][2] = - (a1*sinp) + (A[2][2]*cosp);
       }

        // Sets the matrix to be a transform matrix for rotation
        // about the x-axis.
        inline void SetRotX(double theta)
        {
            A[0][0] =   1.0;
            A[0][1] =   A[0][2] = A[1][0] = A[2][0] = 0.0;
            A[1][1] =   (A[2][2] = cos(theta));
            A[1][2] = - (A[2][1] = sin(theta));
        }

        //! Sets the matrix to be a transform matrix for rotations about the
        //! x-axis and the z-axis.
        inline void Rotate(
            double theta, //!< X-axis rotation (radians).
            double phi    //!< Z-axis rotation (radians).
            )
        {
            // M = Z x X.

            //! Precalculate trig terms.
            double sinp = sin(phi);
            double cosp = cos(phi);
            double sint = sin(theta);
            double cost = cos(theta);

            //! Set matrix.
            A[0][0] =   cosp;
            A[0][1] = - sinp * cost;
            A[0][2] =   sinp * sint;
            A[1][0] =   sinp;
            A[1][1] =   cosp * cost;
            A[1][2] = - cosp * sint;
            A[2][0] =   0.0;
            A[2][1] =   sint;
            A[2][2] =   cost;
        }

        //! Arvo's random rotation matrix.
        inline void rotateArvo(
            double theta, //!< Rotation about the pole. 
            fvector V     //!< Vector for performing the reflection.
            )
        {
            //! Household matrix.
            double minusH[3][3];

            //! Initial rotation matrix.
            double R[3][3];

            minusH[0][0] = 2 * V[0] * V[0] - 1;
            minusH[0][1] = 2 * V[0] * V[1];
            minusH[0][2] = 2 * V[0] * V[2];
            minusH[1][0] = 2 * V[1] * V[0];
            minusH[1][1] = 2 * V[1] * V[1] - 1;
            minusH[1][2] = 2 * V[1] * V[2];
            minusH[2][0] = 2 * V[2] * V[0];
            minusH[2][1] = 2 * V[2] * V[1];
            minusH[2][2] = 2 * V[2] * V[2] - 1;

            //! The rotation matrix is premultiplied by a rotation of pi around
            //! the world z axis as the algorithm leads to incorrect results
            //! when used to sample a perturbation. This flips the signs of the
            //! first two columns of the rotation matrix however this should
            //! not matter for the present simulations:
            //! http://demonstrations.wolfram.com/SamplingAUniformlyRandomRotation/
            R[0][0] = -cos(theta);
            R[0][1] = -sin(theta);
            R[0][2] = 0;
            R[1][0] = sin(theta);
            R[1][1] = -cos(theta);
            R[1][2] = 0;
            R[2][0] = 0;
            R[2][1] = 0;
            R[2][2] = 1;

            //! Construct the final rotation matrix by combining two simple
            //! rotations: first rotate about the the Z axis, then rotate the
            //! Z axis to a random orientation.
            A[0][0] = minusH[0][0] * R[0][0] + minusH[0][1] * R[1][0] + minusH[0][2]* R[2][0];
            A[0][1] = minusH[0][0] * R[0][1] + minusH[0][1] * R[1][1] + minusH[0][2]* R[2][1];
            A[0][2] = minusH[0][0] * R[0][2] + minusH[0][1] * R[1][2] + minusH[0][2]* R[2][2];

            A[1][0] = minusH[1][0] * R[0][0] + minusH[1][1] * R[1][0] + minusH[1][2]* R[2][0];
            A[1][1] = minusH[1][0] * R[0][1] + minusH[1][1] * R[1][1] + minusH[1][2]* R[2][1];
            A[1][2] = minusH[1][0] * R[0][2] + minusH[1][1] * R[1][2] + minusH[1][2]* R[2][2];

            A[2][0] = minusH[2][0] * R[0][0] + minusH[2][1] * R[1][0] + minusH[2][2]* R[2][0];
            A[2][1] = minusH[2][0] * R[0][1] + minusH[2][1] * R[1][1] + minusH[2][2]* R[2][1];
            A[2][2] = minusH[2][0] * R[0][2] + minusH[2][1] * R[1][2] + minusH[2][2]* R[2][2];
        }
    };
};
};
#endif
