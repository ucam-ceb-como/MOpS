/*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2009 Weerapong Phadungsukanan.

  File purpose:
    Header of of the extension of CVODES package.

  Future work:
    Copy constructure is not fully copy CVODES ODE workspace. This
    is need to be finished. Optimisation is needed.

  Licence:
    This file is part of "mops".

    mops is free software; you can redistribute it and/or
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
#ifndef CVODES_UTILS_H
#define CVODES_UTILS_H


#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include "cvodes_impl.h"
#include "cvodes_direct_impl.h"
namespace Mops
{
class CVODES
{
public:
    /**
     * Assign CVODE memory data functions with values copied/cloned.
     */
    static void CVodeMemRecCopy_Serial(CVodeMemRec &mem_dsc, CVodeMemRec &mem_src);
    /**
     * Clone functions with values copied.
     */
    static N_Vector N_VExactClone_Serial(N_Vector w);
private:
    // CVODES auxilary functions
    static N_Vector *N_VExactCloneVectorArray_Serial(int size, N_Vector *w);
    /**
     * Assign functions with values cloned.
     */
    static void N_VectorExactCopy_Serial(N_Vector &dsc, N_Vector &src, const int dsc_malloc = 0, const int src_malloc = 0);
    static void N_VectorArrayExactCopy_Serial(int dsc_size, N_Vector * &dsc, int src_size, N_Vector * &src);

    /**
     * Assign and clone array functions with values copied.
     */
    static void realtypeArrayExactCopy(int dsc_size, realtype * &dsc, int src_size, realtype * &src);
    static void intArrayExactCopy(int dsc_size, int * &dsc, int src_size, int *&src);
    //void longintArrayExactCopy(int dsc_size, long int * &dsc, int src_size, long int *&src);
    static void booleantypeArrayExactCopy(int dsc_size, booleantype * &dsc, int src_size, booleantype *&src);

    /**
     * Assign CVODE memory data functions with values copied/cloned.
     */
    static void CVDlsMemRecCopy_Serial(CVDlsMemRec &lmem_dsc, CVDlsMemRec &lmem_src);
};
};

#endif
