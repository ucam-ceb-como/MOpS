/*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2009 Weerapong Phadungsukanan.

  File purpose:
    Implementation of the extension of CVODES package.

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
#include "cvodes_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


using namespace Mops;
/**
 * Clone an N_Vector with the values copied from w. The normal cvodes N_VClone_Serial
 * does not copy the values.
 */
N_Vector CVODES::N_VExactClone_Serial(N_Vector w) {
    N_Vector v = NULL;
        
    v = N_VClone_Serial(w);
    if (v != NULL) {
        memcpy(NV_DATA_S(v), NV_DATA_S(w), NV_LENGTH_S(w)*sizeof(realtype));
    }
    return v;
}

/**
 * Clone an array of N_Vector with the values copied from *w. The normal cvodes
 * N_VCloneVectorArray_Serial does not copy the values.
 */
N_Vector *CVODES::N_VExactCloneVectorArray_Serial(int size, N_Vector *w) {
    N_Vector *vs = NULL;
    int j;

    if (size <= 0) return(NULL);

    vs = NULL;
    vs = (N_Vector *) malloc(size * sizeof(N_Vector));
    if(vs == NULL) return(NULL);

    for (j = 0; j < size; j++) {
        vs[j] = NULL;
        vs[j] = N_VExactClone_Serial(w[j]);
        if (vs[j] == NULL) {
            N_VDestroyVectorArray_Serial(vs, j-1);
            return(NULL);
        }
    }

    return(vs);
}

/**
 * Make an copy of N_Vector values from src to dsc.
 * TODO: need optimised.
 */
void CVODES::N_VectorExactCopy_Serial(N_Vector &dsc, N_Vector &src, const int dsc_malloc, const int src_malloc) {
    if (src_malloc) {
        if (dsc_malloc) {
            // TODO : optimise with [check the lenght] whether it should be destroy;
            // or destroy and copy
            N_VDestroy_Serial(dsc);
        }
        dsc = N_VExactClone_Serial(src);
    } else {
        if (dsc_malloc) {
            // delete dsc
            N_VDestroy_Serial(dsc);
            dsc = NULL;
        }
        // Nothing to do, because neither vector is initialised
    }
}

/**
 * Make an copy of an array of N_Vector with exact value cloned from src to dsc.
 * TODO: need optimised.
 */
void CVODES::N_VectorArrayExactCopy_Serial(int dsc_size, N_Vector * &dsc, int src_size, N_Vector * &src) {
    N_VDestroyVectorArray_Serial(dsc, dsc_size);
    dsc = NULL;
    dsc = N_VExactCloneVectorArray_Serial(src_size, src);
}

/**
 * Copy realtype array from src to dsc. Not sure if it is used and correctly implemented.
 */
void CVODES::realtypeArrayExactCopy(int dsc_size, realtype * &dsc, int src_size, realtype *&src) {
    if (src != NULL) {
        if ((dsc_size == src_size) && (src_size > 0)) {
            memcpy(dsc, src, src_size*sizeof(realtype));
        } else {
            if (dsc != NULL) {
                delete [] dsc;
                dsc = NULL;
            }
            if (src_size > 0) {
                dsc = new realtype[src_size];
                memcpy(dsc, src, src_size*sizeof(realtype));
            }
        }
    } else {
        if (dsc != NULL) {
            delete [] dsc;
            dsc = NULL;
        }
    }
}

/**
 * Copy int array from src to dsc. Not sure if it is used and correctly implemented.
 */
void CVODES::intArrayExactCopy(int dsc_size, int * &dsc, int src_size, int *&src) {
    if (src != NULL) {
        if ((dsc_size == src_size) && (src_size > 0)) {
            memcpy(dsc, src, src_size*sizeof(int));
        } else {
            if (dsc != NULL) {
                delete [] dsc;
                dsc = NULL;
            }
            if (src_size > 0) {
                dsc = new int[src_size];
                memcpy(dsc, src, src_size*sizeof(int));
            }
        }
    } else {
        if (dsc != NULL) {
            delete [] dsc;
            dsc = NULL;
        }
    }
}

///**
// * Copy long int array from src to dsc. Not sure if it is used and correctly implemented.
// */
//void CVODES::longintArrayExactCopy(long int dsc_size, long int * &dsc, long int src_size, long int *&src) {
//    if (src != NULL) {
//        if ((dsc_size == src_size) && (src_size > 0)) {
//            memcpy(dsc, src, src_size*sizeof(long int));
//        } else {
//            if (dsc != NULL) {
//                delete [] dsc;
//                dsc = NULL;
//            }
//            if (src_size > 0) {
//                dsc = new long int[src_size];
//                memcpy(dsc, src, src_size*sizeof(long int));
//            }
//        }
//    } else {
//        if (dsc != NULL) {
//            delete [] dsc;
//            dsc = NULL;
//        }
//    }
//}

/**
 * Copy booleantype array from src to dsc. Not sure if it is used and correctly implemented.
 */
void CVODES::booleantypeArrayExactCopy(int dsc_size, booleantype * &dsc, int src_size, booleantype *&src) {
    if (src != NULL) {
        if ((dsc_size == src_size) && (src_size > 0)) {
            memcpy(dsc, src, src_size*sizeof(booleantype));
        } else {
            if (dsc != NULL) {
                delete [] dsc;
                dsc = NULL;
            }
            if (src_size > 0) {
                dsc = new booleantype[src_size];
                memcpy(dsc, src, src_size*sizeof(booleantype));
            }
        }
    } else {
        if (dsc != NULL) {
            delete [] dsc;
            dsc = NULL;
        }
    }
}
/**
 * Copy CVDlsMemRec workspace from lmem_src to lmem_dsc. This function is used
 * in CVodeMemRecCopy_Serial to copy linear solver of CVDlsMemRec type.
 */
void CVODES::CVDlsMemRecCopy_Serial(CVDlsMemRec &lmem_dsc, CVDlsMemRec &lmem_src) {

  //int d_type;             /* SUNDIALS_DENSE or SUNDIALS_BAND              */
    lmem_dsc.d_type = lmem_src.d_type;

  //int d_n;                /* problem dimension                            */
    int d_n = lmem_dsc.d_n;
    lmem_dsc.d_n = lmem_src.d_n;

  //int d_ml;               /* lower bandwidth of Jacobian                  */
    lmem_dsc.d_ml = lmem_src.d_ml;
  //int d_mu;               /* upper bandwidth of Jacobian                  */ 
    lmem_dsc.d_mu = lmem_src.d_mu;
  //int d_smu;              /* upper bandwith of M = MIN(N-1,d_mu+d_ml)     */
    lmem_dsc.d_smu = lmem_src.d_smu;

  //booleantype d_jacDQ;    /* TRUE if using internal DQ Jacobian approx.   */
    lmem_dsc.d_jacDQ = lmem_src.d_jacDQ;
  //CVDlsDenseJacFn d_djac; /* dense Jacobian routine to be called          */
    lmem_dsc.d_djac = lmem_src.d_djac;
  //CVDlsBandJacFn d_bjac;  /* band Jacobian routine to be called           */
    lmem_dsc.d_bjac = lmem_src.d_bjac;
  //void *d_J_data;         /* data pointer passed to djac or bjac          */
    lmem_dsc.d_J_data = lmem_src.d_J_data;

  //DlsMat d_M;             /* M = I - gamma * df/dy                        */
    if (lmem_src.d_M != NULL) {
        if ((lmem_dsc.d_M->M == lmem_src.d_M->M) && (lmem_dsc.d_M->N == lmem_src.d_M->N)) {
            // Just copy the memory, no need to reallocate.
            DenseCopy(lmem_src.d_M, lmem_dsc.d_M);
        } else {
            // Need to reallocate data.  First delete current data.
            DestroyMat(lmem_dsc.d_M);
            // Allocate memory to matrix..
            lmem_dsc.d_M = NewDenseMat(lmem_src.d_M->M, lmem_src.d_M->N);
            // Copy data.
            DenseCopy(lmem_src.d_M, lmem_dsc.d_M);
        }
    } else {
        // Clear matrix.
        DestroyMat(lmem_dsc.d_M);
        lmem_dsc.d_M = NULL;
    }
  //DlsMat d_savedJ;        /* savedJ = old Jacobian                        */
    if (lmem_src.d_savedJ != NULL) {
        if ((lmem_src.d_savedJ->M == lmem_dsc.d_savedJ->M) && (lmem_src.d_savedJ->N == lmem_dsc.d_savedJ->N)) {
            // Just copy the memory, no need to reallocate.
            DenseCopy(lmem_src.d_savedJ, lmem_dsc.d_savedJ);
        } else {
            // Need to reallocate data.  First delete current data.
            DestroyMat(lmem_dsc.d_savedJ);
            // Allocate memory to matrix..
            lmem_dsc.d_savedJ = NewDenseMat(lmem_src.d_savedJ->M, lmem_src.d_savedJ->N);
            // Copy data.
            DenseCopy(lmem_src.d_savedJ, lmem_dsc.d_savedJ);
        }
    } else {
        // Clear matrix.
        DestroyMat(lmem_dsc.d_savedJ);
        lmem_dsc.d_savedJ = NULL;
    }

  //int *d_pivots;          /* pivots = pivot array for PM = LU             */
    intArrayExactCopy(d_n, lmem_dsc.d_pivots, lmem_src.d_n, lmem_src.d_pivots);
  //
  //long int  d_nstlj;      /* nstlj = nst at last Jacobian eval.           */
    lmem_dsc.d_nstlj = lmem_src.d_nstlj;

  //long int d_nje;         /* nje = no. of calls to jac                    */
    lmem_dsc.d_nje = lmem_src.d_nje;

  //long int d_nfeDQ;       /* no. of calls to f due to DQ Jacobian approx. */
    lmem_dsc.d_nfeDQ = lmem_src.d_nfeDQ;

  //int d_last_flag;        /* last error return flag                       */
    lmem_dsc.d_last_flag = lmem_src.d_last_flag;
}

/**
 * Copy CVodeMemRec workspace from mem_src to mem_dsc.
 * Currently, some variables for sensitivity and quadrature are ignored as ther are not
 * used anywhere (Sensitivity cannot be run in any other mode except GPC and There is no
 * solver which has to use quadrature integration.
 */
void CVODES::CVodeMemRecCopy_Serial(CVodeMemRec &mem_dsc, CVodeMemRec &mem_src) {
    //  realtype cv_uround;         /* machine unit roundoff                        */   
        mem_dsc.cv_uround           = mem_src.cv_uround;

    //  /*-------------------------- 
    //    Problem Specification Data 
    //    --------------------------*/
    //  CVRhsFn cv_f;               /* y' = f(t,y(t))                               */
        mem_dsc.cv_f            = mem_src.cv_f;
    //  void *cv_user_data;         /* user pointer passed to f                     */
        //mem_dsc.cv_user_data    = mem_src.cv_user_data;
    //  int cv_lmm;                 /* lmm = ADAMS or BDF                           */
        mem_dsc.cv_lmm          = mem_src.cv_lmm;
    //  int cv_iter;                /* iter = FUNCTIONAL or NEWTON                  */
        mem_dsc.cv_iter         = mem_src.cv_iter;   // CV_FUNCTIONAL or CV_NEWTON.
    //  int cv_itol;                /* itol = CV_SS, CV_SV, or CV_WF, or CV_NN      */
        mem_dsc.cv_itol         = mem_src.cv_itol;
    //  realtype cv_reltol;         /* relative tolerance                           */
        mem_dsc.cv_reltol       = mem_src.cv_reltol;
    //  realtype cv_Sabstol;        /* scalar absolute tolerance                    */
        mem_dsc.cv_Sabstol      = mem_src.cv_Sabstol;
    //  N_Vector cv_Vabstol;        /* vector absolute tolerance                    */
        N_VectorExactCopy_Serial(mem_dsc.cv_Vabstol, mem_src.cv_Vabstol, mem_dsc.cv_VabstolMallocDone, mem_src.cv_VabstolMallocDone);
    //  booleantype cv_user_efun;   /* TRUE if user sets efun                       */
        mem_dsc.cv_user_efun    = mem_src.cv_user_efun;
    //  CVEwtFn cv_efun;            /* function to set ewt                          */
        mem_dsc.cv_efun         = mem_src.cv_efun;
    //  void *cv_e_data;            /* user pointer passed to efun                  */
        mem_dsc.cv_e_data       = mem_src.cv_e_data;

    //  /*-----------------------
    //    Quadrature Related Data 
    //    -----------------------*/
    //  booleantype cv_quadr;       /* TRUE if integrating quadratures              */
        mem_dsc.cv_quadr        = mem_src.cv_quadr;
    //  CVQuadRhsFn cv_fQ;          /* q' = fQ(t, y(t))                             */
        mem_dsc.cv_fQ           = mem_src.cv_fQ;
    //  booleantype cv_errconQ;     /* TRUE if quadrs. are included in error test   */
        mem_dsc.cv_errconQ      = mem_src.cv_errconQ;
    //  int cv_itolQ;               /* itolQ = CV_SS or CV_SV                       */
        mem_dsc.cv_itolQ        = mem_src.cv_itolQ;
    //  realtype cv_reltolQ;        /* relative tolerance for quadratures           */
        mem_dsc.cv_reltolQ      = mem_src.cv_reltolQ;
    //  realtype cv_SabstolQ;       /* scalar absolute tolerance for quadratures    */
        mem_dsc.cv_SabstolQ     = mem_src.cv_SabstolQ;
    //  N_Vector cv_VabstolQ;       /* vector absolute tolerance for quadratures    */
        N_VectorExactCopy_Serial(mem_dsc.cv_VabstolQ, mem_src.cv_VabstolQ, mem_dsc.cv_VabstolQMallocDone, mem_src.cv_VabstolQMallocDone);

    //  /*------------------------
    //    Sensitivity Related Data 
    //    ------------------------*/
    //  booleantype cv_sensi;       /* TRUE if computing sensitivities              */
        mem_dsc.cv_sensi        = mem_src.cv_sensi;
    //  int cv_Ns;                  /* Number of sensitivities                      */
        // For N_Vector copying.
        mem_dsc.cv_Ns           = mem_src.cv_Ns;
    //  int cv_ism;                 /* ism = SIMULTANEOUS or STAGGERED              */
        mem_dsc.cv_ism          = mem_src.cv_ism;
    //  CVSensRhsFn cv_fS;          /* fS = (df/dy)*yS + (df/dp)                    */
        mem_dsc.cv_fS           = mem_src.cv_fS;
    //  CVSensRhs1Fn cv_fS1;        /* fS1 = (df/dy)*yS_i + (df/dp)                 */
        mem_dsc.cv_fS1          = mem_src.cv_fS1;
    //  void *cv_fS_data;           /* data pointer passed to fS                    */
        mem_dsc.cv_fS_data      = mem_src.cv_fS_data;
    //  booleantype cv_fSDQ;        /* TRUE if using internal DQ functions          */
        mem_dsc.cv_fSDQ         = mem_src.cv_fSDQ;
    //  int cv_ifS;                 /* ifS = ALLSENS or ONESENS                     */
        mem_dsc.cv_ifS          = mem_src.cv_ifS;
    //  realtype *cv_p;             /* parameters in f(t,y,p)                       */
        mem_dsc.cv_p            = mem_src.cv_p;        // no copy here since cvode do not free this memory
//TODO:Not needed  realtype *cv_pbar;          /* scale factors for parameters                 */
//TODO:Not needed        realtypeArrayExactCopy(old_Ns, mem_dsc.cv_pbar, mem_src.cv_Ns, mem_src.cv_pbar);
//TODO:Not needed  int *cv_plist;              /* list of sensitivities                        */
//TODO:Not needed        intArrayExactCopy(old_Ns, mem_dsc.cv_plist, mem_src.cv_Ns, mem_src.cv_plist);
    //  int cv_DQtype;              /* central/forward finite differences           */
        mem_dsc.cv_DQtype       = mem_src.cv_DQtype;
    //  realtype cv_DQrhomax;       /* cut-off value for separate/simultaneous FD   */
        mem_dsc.cv_DQrhomax     = mem_src.cv_DQrhomax;
    //  booleantype cv_errconS;     /* TRUE if yS are considered in err. control    */
        mem_dsc.cv_errconS      = mem_src.cv_errconS;
    //  int cv_itolS;
        mem_dsc.cv_itolS        = mem_src.cv_itolS;
    //  realtype cv_reltolS;        /* relative tolerance for sensitivities         */
        mem_dsc.cv_reltolS      = mem_src.cv_reltolS;
//TODO:Not needed  realtype *cv_SabstolS;      /* scalar absolute tolerances for sensi.        */
//TODO:Not needed        realtypeArrayExactCopy(old_Ns, mem_dsc.cv_SabstolS, mem_src.cv_Ns, mem_src.cv_SabstolS);
//TODO:Not needed  N_Vector *cv_VabstolS;      /* vector absolute tolerances for sensi.        */
//TODO:Not needed        N_VectorArrayExactCopy_Serial(old_Ns, mem_dsc.cv_VabstolS, mem_src.cv_Ns, mem_src.cv_VabstolS);

    //  /*-----------------------------------
    //    Quadrature Sensitivity Related Data 
    //    -----------------------------------*/

    //  booleantype cv_quadr_sensi; /* TRUE if computing sensitivties of quadrs.    */
        mem_dsc.cv_quadr_sensi  = mem_src.cv_quadr_sensi;

    //  CVQuadSensRhsFn cv_fQS;     /* fQS = (dfQ/dy)*yS + (dfQ/dp)                 */
        mem_dsc.cv_fQS          = mem_src.cv_fQS;
    //  void *cv_fQS_data;          /* data pointer passed to fQS                   */
        mem_dsc.cv_fQS_data     = mem_src.cv_fQS_data;
    //  booleantype cv_fQSDQ;       /* TRUE if using internal DQ functions          */
        mem_dsc.cv_fQSDQ        = mem_src.cv_fQSDQ;
    //  booleantype cv_errconQS;    /* TRUE if yQS are considered in err. con.      */
        mem_dsc.cv_errconQS     = mem_src.cv_errconQS;
    //  int cv_itolQS;
        mem_dsc.cv_itolQS       = mem_src.cv_itolQS;
    //  realtype cv_reltolQS;       /* relative tolerance for yQS                   */
        mem_dsc.cv_reltolQS     = mem_src.cv_reltolQS;
//TODO:Not needed  realtype *cv_SabstolQS;     /* scalar absolute tolerances for yQS           */
//TODO:Not needed        realtypeArrayExactCopy(old_Ns, mem_dsc.cv_SabstolQS, mem_src.cv_Ns, mem_src.cv_SabstolQS);
//TODO:Not needed  N_Vector *cv_VabstolQS;     /* vector absolute tolerances for yQS           */
//TODO:Not needed        N_VectorArrayExactCopy_Serial(old_Ns, mem_dsc.cv_VabstolQS, mem_src.cv_Ns, mem_src.cv_VabstolQS);

    //  /*-----------------------
    //    Nordsieck History Array 
    //    -----------------------*/
    //  N_Vector cv_zn[L_MAX];      /* Nordsieck array, of size N x (q+1).
    //                                 zn[j] is a vector of length N (j=0,...,q)
    //                                 zn[j] = [1/factorial(j)] * h^j * 
    //                                 (jth derivative of the interpolating poly.)  */
    //

        for (int i=0; i < mem_src.cv_q; i++) {
            N_VDestroy_Serial(mem_dsc.cv_zn[i]);
            mem_dsc.cv_zn[i] = N_VExactClone_Serial(mem_src.cv_zn[i]);
        }
    //  /*-------------------
    //    Vectors of length N 
    //    -------------------*/

    //  N_Vector cv_ewt;            /* error weight vector                          */
        N_VectorExactCopy_Serial(mem_dsc.cv_ewt, mem_src.cv_ewt);
    //  N_Vector cv_y;              /* y is used as temporary storage by the solver.
    //                                 The memory is provided by the user to CVode 
    //                                 where the vector is named yout.              */
        N_VectorExactCopy_Serial(mem_dsc.cv_y, mem_src.cv_y);
    //  N_Vector cv_acor;           /* In the context of the solution of the
    //                                 nonlinear equation, acor = y_n(m) - y_n(0).
    //                                 On return, this vector is scaled to give
    //                                 the estimated local error in y.              */
        N_VectorExactCopy_Serial(mem_dsc.cv_acor, mem_src.cv_acor);
    //  N_Vector cv_tempv;          /* temporary storage vector                     */
        N_VectorExactCopy_Serial(mem_dsc.cv_tempv, mem_src.cv_tempv);
    //  N_Vector cv_ftemp;          /* temporary storage vector                     */
        N_VectorExactCopy_Serial(mem_dsc.cv_ftemp, mem_src.cv_ftemp);

    //  /*--------------------------
    //    Quadrature Related Vectors 
    //    --------------------------*/
    //  N_Vector cv_znQ[L_MAX];     /* Nordsieck arrays for quadratures             */
//TODO:Not needed        for (int i=0; i < mem_src.cv_q; i++) {
//TODO:Not needed            N_VDestroy_Serial(mem_src.cv_znQ[i]);
//TODO:Not needed            mem_src.cv_znQ[i] = N_VExactClone_Serial(mem_dsc.cv_znQ[i]);
//TODO:Not needed        }
    //  N_Vector cv_ewtQ;           /* error weight vector for quadratures          */
        N_VectorExactCopy_Serial(mem_src.cv_ewtQ, mem_dsc.cv_ewtQ);
    //  N_Vector cv_yQ;             /* Unlike y, yQ is not allocated by the user    */
        N_VectorExactCopy_Serial(mem_src.cv_yQ, mem_dsc.cv_yQ);
    //  N_Vector cv_acorQ;          /* acorQ = yQ_n(m) - yQ_n(0)                    */
        N_VectorExactCopy_Serial(mem_src.cv_acorQ, mem_dsc.cv_acorQ);
    //  N_Vector cv_tempvQ;         /* temporary storage vector (~ tempv)           */
        N_VectorExactCopy_Serial(mem_src.cv_tempvQ, mem_dsc.cv_tempvQ);

//TODO:Not needed  /*---------------------------
//TODO:Not needed    Sensitivity Related Vectors 
//TODO:Not needed    ---------------------------*/
//TODO:Not needed  N_Vector *cv_znS[L_MAX];    /* Nordsieck arrays for sensitivities           */
//TODO:Not needed  N_Vector *cv_ewtS;          /* error weight vectors for sensitivities       */
//TODO:Not needed  N_Vector *cv_yS;            /* yS=yS0 (allocated by the user)               */
//TODO:Not needed  N_Vector *cv_acorS;         /* acorS = yS_n(m) - yS_n(0)                    */
//TODO:Not needed  N_Vector *cv_tempvS;        /* temporary storage vector (~ tempv)           */
//TODO:Not needed  N_Vector *cv_ftempS;        /* temporary storage vector (~ ftemp)           */

//TODO:Not needed  booleantype cv_stgr1alloc;  /* Did we allocate ncfS1, ncfnS1, and nniS1?    */

//TODO:Not needed  /*--------------------------------------
//TODO:Not needed    Quadrature Sensitivity Related Vectors 
//TODO:Not needed    --------------------------------------*/
//TODO:Not needed
//TODO:Not needed  N_Vector *cv_znQS[L_MAX];   /* Nordsieck arrays for quadr. sensitivities    */
//TODO:Not needed  N_Vector *cv_ewtQS;         /* error weight vectors for sensitivities       */
//TODO:Not needed  N_Vector *cv_yQS;           /* Unlike yS, yQS is not allocated by the user  */
//TODO:Not needed  N_Vector *cv_acorQS;        /* acorQS = yQS_n(m) - yQS_n(0)                 */
//TODO:Not needed  N_Vector *cv_tempvQS;       /* temporary storage vector (~ tempv)           */
//TODO:Not needed  N_Vector cv_ftempQ;         /* temporary storage vector (~ ftemp)           */
  
    //  /*-----------------
    //    Tstop information
    //    -----------------*/
    //  booleantype cv_tstopset;
        mem_dsc.cv_tstopset     = mem_src.cv_tstopset;
    //  realtype cv_tstop;
        mem_dsc.cv_tstop        = mem_src.cv_tstop;
    //
    //  /*---------
    //    Step Data 
    //    ---------*/
    //  int cv_q;                    /* current order                               */
        mem_dsc.cv_q            = mem_src.cv_q;
    //  int cv_qprime;               /* order to be used on the next step
        mem_dsc.cv_qprime       = mem_src.cv_qprime;
    //                                * qprime = q-1, q, or q+1                     */
    //  int cv_next_q;               /* order to be used on the next step           */
        mem_dsc.cv_next_q       = mem_src.cv_next_q;
    //  int cv_qwait;                /* number of internal steps to wait before
    //                                * considering a change in q                   */
        mem_dsc.cv_qwait        = mem_src.cv_qwait;
    //  int cv_L;                    /* L = q + 1                                   */
        mem_dsc.cv_L            = mem_src.cv_L;
    //  realtype cv_hin;
        mem_dsc.cv_hin          = mem_src.cv_hin;
    //  realtype cv_h;               /* current step size                           */
        mem_dsc.cv_h            = mem_src.cv_h;
    //  realtype cv_hprime;          /* step size to be used on the next step       */ 
        mem_dsc.cv_hprime       = mem_src.cv_hprime;
    //  realtype cv_next_h;          /* step size to be used on the next step       */ 
        mem_dsc.cv_next_h       = mem_src.cv_next_h;
    //  realtype cv_eta;             /* eta = hprime / h                            */
        mem_dsc.cv_eta          = mem_src.cv_eta;
    //  realtype cv_hscale;          /* value of h used in zn                       */
        mem_dsc.cv_hscale       = mem_src.cv_hscale;
    //  realtype cv_tn;              /* current internal value of t                 */
        mem_dsc.cv_tn           = mem_src.cv_tn;
    //  realtype cv_tretlast;        /* last value of t returned                    */
        mem_dsc.cv_tretlast     = mem_src.cv_tretlast;
    //
    //  realtype cv_tau[L_MAX+1];    /* array of previous q+1 successful step
    //                                * sizes indexed from 1 to q+1                 */
        for (int i=0; i < mem_src.cv_q+1; i++) {
            mem_dsc.cv_tau[i] = mem_src.cv_tau[i];
        }
    //  realtype cv_tq[NUM_TESTS+1]; /* array of test quantities indexed from
    //                                * 1 to NUM_TESTS(=5)                          */
        for (int i=0; i < NUM_TESTS+1; i++) {
            mem_dsc.cv_tq[i] = mem_src.cv_tq[i];
        }
    //  realtype cv_l[L_MAX];        /* coefficients of l(x) (degree q poly)        */
    //
        for (int i=0; i < mem_src.cv_q; i++) {
            mem_dsc.cv_l[i] = mem_src.cv_l[i];
        }
    //  realtype cv_rl1;             /* the scalar 1/l[1]                           */
        mem_dsc.cv_rl1          = mem_src.cv_rl1;
    //  realtype cv_gamma;           /* gamma = h * rl1                             */
        mem_dsc.cv_gamma        = mem_src.cv_gamma;
    //  realtype cv_gammap;          /* gamma at the last setup call                */
        mem_dsc.cv_gammap       = mem_src.cv_gammap;
    //  realtype cv_gamrat;          /* gamma / gammap                              */
        mem_dsc.cv_gamrat       = mem_src.cv_gamrat;
    //  realtype cv_crate;           /* est. corrector conv. rate in Nls            */
        mem_dsc.cv_crate        = mem_src.cv_crate;
    //  realtype cv_crateS;          /* est. corrector conv. rate in NlsStgr        */
        mem_dsc.cv_crateS       = mem_src.cv_crateS;
    //  realtype cv_acnrm;           /* | acor |                                    */
        mem_dsc.cv_acnrm        = mem_src.cv_acnrm;
    //  realtype cv_acnrmQ;          /* | acorQ |                                   */
        mem_dsc.cv_acnrmQ       = mem_src.cv_acnrmQ;
    //  realtype cv_acnrmS;          /* | acorS |                                   */
        mem_dsc.cv_acnrmS       = mem_src.cv_acnrmS;
    //  realtype cv_acnrmQS;         /* | acorQS |                                  */
        mem_dsc.cv_acnrmQS      = mem_src.cv_acnrmQS;
    //  realtype cv_nlscoef;         /* coeficient in nonlinear convergence test    */
        mem_dsc.cv_nlscoef      = mem_src.cv_nlscoef;
    //  int  cv_mnewt;               /* Newton iteration counter                    */
        mem_dsc.cv_mnewt        = mem_src.cv_mnewt;
//TODO:Not needed  int  *cv_ncfS1;              /* Array of Ns local counters for conv.  
//TODO:Not needed                                * failures (used in CVStep for STAGGERED1)    */
//TODO:Not needed        intArrayExactCopy(old_Ns, mem_dsc.cv_ncfS1, mem_src.cv_Ns, mem_src.cv_ncfS1);

    //  /*------
    //    Limits 
    //    ------*/
    //
    //  int cv_qmax;             /* q <= qmax                                       */
        mem_dsc.cv_qmax         = mem_src.cv_qmax;
    //  long int cv_mxstep;      /* maximum number of internal steps for one 
    //                  user call                                       */
        mem_dsc.cv_mxstep       = mem_src.cv_mxstep;
    //  int cv_maxcor;           /* maximum number of corrector iterations for 
    //                  the solution of the nonlinear equation          */
        mem_dsc.cv_maxcor       = mem_src.cv_maxcor;
    //  int cv_maxcorS;
        mem_dsc.cv_maxcorS      = mem_src.cv_maxcorS;
    //  int cv_mxhnil;           /* max. number of warning messages issued to the
    //                  user that t + h == t for the next internal step */
        mem_dsc.cv_mxhnil       = mem_src.cv_mxhnil;
    //  int cv_maxnef;           /* maximum number of error test failures           */
        mem_dsc.cv_maxnef       = mem_src.cv_maxnef;
    //  int cv_maxncf;           /* maximum number of nonlinear conv. failures      */
        mem_dsc.cv_maxncf       = mem_src.cv_maxncf;
    //
    //  realtype cv_hmin;        /* |h| >= hmin                                     */
        mem_dsc.cv_hmin         = mem_src.cv_hmin;
    //  realtype cv_hmax_inv;    /* |h| <= 1/hmax_inv                               */
        mem_dsc.cv_hmax_inv     = mem_src.cv_hmax_inv;
    //  realtype cv_etamax;      /* eta <= etamax                                   */
        mem_dsc.cv_etamax       = mem_src.cv_etamax;
    //
    //  /*----------
    //    Counters 
    //    ----------*/
    //
    //  long int cv_nst;         /* number of internal steps taken                  */
        mem_dsc.cv_nst          = mem_src.cv_nst;
     //  long int cv_nfe;         /* number of f calls                               */
        mem_dsc.cv_nfe          = mem_src.cv_nfe;
    //  long int cv_nfQe;        /* number of fQ calls                              */
        mem_dsc.cv_nfQe         = mem_src.cv_nfQe;
    //  long int cv_nfSe;        /* number of fS calls                              */
        mem_dsc.cv_nfSe         = mem_src.cv_nfSe;
    //  long int cv_nfeS;        /* number of f calls from sensi DQ                 */
        mem_dsc.cv_nfeS         = mem_src.cv_nfeS;
    //  long int cv_nfQSe;       /* number of fQS calls                             */
        mem_dsc.cv_nfQSe        = mem_src.cv_nfQSe;
    //  long int cv_nfQeS;       /* number of fQ calls from sensi DQ                */
        mem_dsc.cv_nfQeS        = mem_src.cv_nfQeS;
    //  long int cv_ncfn;        /* number of corrector convergence failures        */
        mem_dsc.cv_ncfn         = mem_src.cv_ncfn;
    //  long int cv_ncfnS;       /* number of total sensi. corr. conv. failures     */
        mem_dsc.cv_ncfnS = mem_src.cv_ncfnS;
//TODO:Not needed  long int *cv_ncfnS1;     /* number of sensi. corrector conv. failures       */
//TODO:Not needed       longintArrayExactCopy(old_ncfn, mem_dsc.cv_ncfnS1, mem_src.cv_ncfn, mem_src.cv_ncfnS1);
    //  long int cv_nni;         /* number of nonlinear iterations performed        */
        mem_dsc.cv_nni = mem_src.cv_nni;
    //  long int cv_nniS;        /* number of total sensi. nonlinear iterations     */
        mem_dsc.cv_nniS = mem_src.cv_nniS;
//TODO:Not needed   long int *cv_nniS1;      /* number of sensi. nonlinear iterations           */
//TODO:Not needed       longintArrayExactCopy(old_nniS, mem_dsc.cv_nniS1, mem_src.cv_nniS, mem_src.cv_nniS1);
    //  long int cv_netf;        /* number of error test failures                   */
        mem_dsc.cv_netf         = mem_src.cv_netf;
    //  long int cv_netfQ;       /* number of quadr. error test failures            */
        mem_dsc.cv_netfQ        = mem_src.cv_netfQ;
    //  long int cv_netfS;       /* number of sensi. error test failures            */
        mem_dsc.cv_netfS        = mem_src.cv_netfS;
    //  long int cv_netfQS;      /* number of quadr. sensi. error test failures     */
        mem_dsc.cv_netfQS       = mem_src.cv_netfQS;
    //  long int cv_nsetups;     /* number of setup calls                           */
        mem_dsc.cv_nsetups      = mem_src.cv_nsetups;
    //  long int cv_nsetupsS;    /* number of setup calls due to sensitivities      */
        mem_dsc.cv_nsetupsS     = mem_src.cv_nsetupsS;
    //  int cv_nhnil;            /* number of messages issued to the user that
    //                  t + h == t for the next iternal step            */
        mem_dsc.cv_nhnil        = mem_src.cv_nhnil;

    //  /*-----------------------------
    //    Space requirements for CVODES 
    //    -----------------------------*/
    //  long int cv_lrw1;        /* no. of realtype words in 1 N_Vector y           */ 
        mem_dsc.cv_lrw1         = mem_src.cv_lrw1;
    //  long int cv_liw1;        /* no. of integer words in 1 N_Vector y            */ 
        mem_dsc.cv_liw1         = mem_src.cv_liw1;
    //  long int cv_lrw1Q;       /* no. of realtype words in 1 N_Vector yQ          */ 
        mem_dsc.cv_lrw1Q        = mem_src.cv_lrw1Q;
    //  long int cv_liw1Q;       /* no. of integer words in 1 N_Vector yQ           */ 
        mem_dsc.cv_liw1Q        = mem_src.cv_liw1Q;
    //  long int cv_lrw;         /* no. of realtype words in CVODES work vectors    */
        mem_dsc.cv_lrw          = mem_src.cv_lrw;
    //  long int cv_liw;         /* no. of integer words in CVODES work vectors     */
        mem_dsc.cv_liw          = mem_src.cv_liw;

    //  /*----------------
    //    Step size ratios
    //    ----------------*/
    //  realtype cv_etaqm1;      /* ratio of new to old h for order q-1             */
        mem_dsc.cv_etaqm1       = mem_src.cv_etaqm1;
    //  realtype cv_etaq;        /* ratio of new to old h for order q               */
        mem_dsc.cv_etaq         = mem_src.cv_etaq;
    //  realtype cv_etaqp1;      /* ratio of new to old h for order q+1             */
        mem_dsc.cv_etaqp1       = mem_src.cv_etaqp1;

    // These are functions. No need to copy.
    //  /*------------------
    //    Linear Solver Data 
    //    ------------------*/
    //
    //  /* Linear Solver functions to be called */
    //
    //  int (*cv_linit)(struct CVodeMemRec *cv_mem);
    //
    //  int (*cv_lsetup)(struct CVodeMemRec *cv_mem, int convfail, 
    //           N_Vector ypred, N_Vector fpred, booleantype *jcurPtr, 
    //           N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3); 
    //
    //  int (*cv_lsolve)(struct CVodeMemRec *cv_mem, N_Vector b, N_Vector weight,
    //           N_Vector ycur, N_Vector fcur);
    //
    //  void (*cv_lfree)(struct CVodeMemRec *cv_mem);

    //  /* Linear Solver specific memory */
    //  void *cv_lmem;
        // Assuming a particular linear solver. 
        // TODO:It will need to be looked at if other linear solver are used (they can be the same, it has not been investigated.).
        if (mem_src.cv_lmem != NULL) {
            CVDlsMemRec &lmem_dsc  = *((CVDlsMem)mem_dsc.cv_lmem);
            CVDlsMemRec &lmem_src = *((CVDlsMem)mem_src.cv_lmem);
            CVDlsMemRecCopy_Serial(lmem_dsc, lmem_src);
        } else {
            mem_dsc.cv_lfree(&mem_dsc);
        }

    //  /* Flag to request a call to the setup routine */
    //  booleantype cv_forceSetup;
        mem_dsc.cv_forceSetup   = mem_src.cv_forceSetup;

    //  /*------------
    //    Saved Values
    //    ------------*/
    //  int cv_qu;                   /* last successful q value used                */
        mem_dsc.cv_qu           = mem_src.cv_qu;
    //  long int cv_nstlp;           /* step number of last setup call              */
        mem_dsc.cv_nstlp        = mem_src.cv_nstlp;
    //  realtype cv_h0u;             /* actual initial stepsize                     */
        mem_dsc.cv_h0u          = mem_src.cv_h0u;
    //  realtype cv_hu;              /* last successful h value used                */
        mem_dsc.cv_hu           = mem_src.cv_hu;
    //  realtype cv_saved_tq5;       /* saved value of tq[5]                        */
        mem_dsc.cv_saved_tq5    = mem_src.cv_saved_tq5;
    //  booleantype cv_jcur;         /* is Jacobian info for linear solver current? */
        mem_dsc.cv_jcur         = mem_src.cv_jcur;
    //  realtype cv_tolsf;           /* tolerance scale factor                      */
        mem_dsc.cv_tolsf        = mem_src.cv_tolsf;
    //  int cv_qmax_alloc;           /* qmax used when allocating mem               */
        mem_dsc.cv_qmax_alloc   = mem_src.cv_qmax_alloc;
    //  int cv_qmax_allocQ;          /* qmax used when allocating quad. mem         */
        mem_dsc.cv_qmax_allocQ  = mem_src.cv_qmax_allocQ;
    //  int cv_qmax_allocS;          /* qmax used when allocating sensi. mem        */
        mem_dsc.cv_qmax_allocS  = mem_src.cv_qmax_allocS;
    //  int cv_qmax_allocQS;         /* qmax used when allocating quad. sensi. mem  */
        mem_dsc.cv_qmax_allocQS = mem_src.cv_qmax_allocQS;
    //  int cv_indx_acor;            /* index of zn vector in which acor is saved   */
        mem_dsc.cv_indx_acor    = mem_src.cv_indx_acor;
    //  booleantype cv_setupNonNull; /* Does setup do something?                    */
        mem_dsc.cv_setupNonNull = mem_src.cv_setupNonNull;

    //  /*--------------------------------------------------------------------
    //    Flags turned ON by CVodeInit, CVodeSensMalloc, and CVodeQuadMalloc 
    //    and read by CVodeReInit, CVodeSensReInit, and CVodeQuadReInit
    //    --------------------------------------------------------------------*/
    //  booleantype cv_VabstolMallocDone;
        mem_dsc.cv_VabstolMallocDone    = mem_src.cv_VabstolMallocDone;
    //  booleantype cv_MallocDone;
        mem_dsc.cv_MallocDone           = mem_src.cv_MallocDone;
    //  booleantype cv_VabstolQMallocDone;
        mem_dsc.cv_VabstolQMallocDone   = mem_src.cv_VabstolQMallocDone;
    //  booleantype cv_QuadMallocDone;
        mem_dsc.cv_QuadMallocDone       = mem_src.cv_QuadMallocDone;
    //  booleantype cv_VabstolSMallocDone;
        mem_dsc.cv_VabstolSMallocDone   = mem_src.cv_VabstolSMallocDone;
    //  booleantype cv_SabstolSMallocDone;
        mem_dsc.cv_SabstolSMallocDone   = mem_src.cv_SabstolSMallocDone;
    //  booleantype cv_SensMallocDone;
        mem_dsc.cv_SensMallocDone       = mem_src.cv_SensMallocDone;
    //  booleantype cv_VabstolQSMallocDone;
        mem_dsc.cv_VabstolQSMallocDone  = mem_src.cv_VabstolQSMallocDone;
    //  booleantype cv_SabstolQSMallocDone;
        mem_dsc.cv_SabstolQSMallocDone  = mem_src.cv_SabstolQSMallocDone;
    //  booleantype cv_QuadSensMallocDone;
        mem_dsc.cv_QuadSensMallocDone   = mem_src.cv_QuadSensMallocDone;

    //  /*-------------------------------------------
    //    Error handler function and error ouput file 
    //    -------------------------------------------*/
    //  CVErrHandlerFn cv_ehfun;    /* Error messages are handled by ehfun          */
        mem_dsc.cv_ehfun            = mem_src.cv_ehfun;
    //  void *cv_eh_data;           /* dats pointer passed to ehfun                 */
        mem_dsc.cv_eh_data          = mem_src.cv_eh_data;
    //  FILE *cv_errfp;             /* CVODES error messages are sent to errfp      */    
        mem_dsc.cv_errfp            = mem_src.cv_errfp;

    //  /*-------------------------
    //    Stability Limit Detection
    //    -------------------------*/
    //  booleantype cv_sldeton;     /* Is Stability Limit Detection on?             */
        mem_dsc.cv_sldeton          = mem_src.cv_sldeton;
//TODO:Not needed  realtype cv_ssdat[6][4];    /* scaled data array for STALD                  */
//TODO:Not needed        mem_dsc.                    = mem_src.;
    //  int cv_nscon;               /* counter for STALD method                     */
        mem_dsc.cv_nscon            = mem_src.cv_nscon;
    //  long int cv_nor;            /* counter for number of order reductions       */
        mem_dsc.cv_nor              = mem_src.cv_nor;

    //  /*----------------
    //    Rootfinding Data
    //    ----------------*/
    //  CVRootFn cv_gfun;        /* Function g for roots sought                     */
        mem_dsc.cv_gfun             = mem_src.cv_gfun;
    //  int cv_nrtfn;            /* number of components of g                       */
        int old_nrtfn               = mem_dsc.cv_nrtfn;
        mem_dsc.cv_nrtfn            = mem_src.cv_nrtfn;
    //  int *cv_iroots;          /* array for root information                      */
        intArrayExactCopy(old_nrtfn, mem_dsc.cv_iroots, mem_src.cv_nrtfn, mem_src.cv_iroots);
    //  int *cv_rootdir;         /* array specifying direction of zero-crossing     */
        intArrayExactCopy(old_nrtfn, mem_dsc.cv_rootdir, mem_src.cv_nrtfn, mem_src.cv_rootdir);
    //  realtype cv_tlo;         /* nearest endpoint of interval in root search     */
        mem_dsc.cv_tlo              = mem_src.cv_tlo;
    //  realtype cv_thi;         /* farthest endpoint of interval in root search    */
        mem_dsc.cv_tlo              = mem_src.cv_tlo;
    //  realtype cv_trout;       /* t value returned by rootfinding routine         */
        mem_dsc.cv_trout            = mem_src.cv_trout;
    //  realtype *cv_glo;        /* saved array of g values at t = tlo              */
        realtypeArrayExactCopy(old_nrtfn, mem_dsc.cv_glo, mem_src.cv_nrtfn, mem_src.cv_glo);
    //  realtype *cv_ghi;        /* saved array of g values at t = thi              */
        realtypeArrayExactCopy(old_nrtfn, mem_dsc.cv_ghi, mem_src.cv_nrtfn, mem_src.cv_ghi);
    //  realtype *cv_grout;      /* array of g values at t = trout                  */
        realtypeArrayExactCopy(old_nrtfn, mem_dsc.cv_grout, mem_src.cv_nrtfn, mem_src.cv_grout);
    //  realtype cv_toutc;       /* copy of tout (if NORMAL mode)                   */
        mem_dsc.cv_toutc            = mem_src.cv_toutc;
    //  realtype cv_ttol;        /* tolerance on root location trout                */
        mem_dsc.cv_ttol             = mem_src.cv_ttol;
    //  int cv_taskc;            /* copy of parameter itask                         */
        mem_dsc.cv_taskc            = mem_src.cv_taskc;
    //  int cv_irfnd;            /* flag showing whether last step had a root       */
        mem_dsc.cv_irfnd            = mem_src.cv_irfnd;
    //  long int cv_nge;         /* counter for g evaluations                       */
        mem_dsc.cv_nge              = mem_src.cv_nge;
    //  booleantype *cv_gactive; /* array with active/inactive event functions      */
        booleantypeArrayExactCopy(old_nrtfn, mem_dsc.cv_gactive, mem_src.cv_nrtfn, mem_src.cv_gactive);
    //  int cv_mxgnull;          /* number of warning messages about possible g==0  */
        mem_dsc.cv_mxgnull          = mem_src.cv_mxgnull;
    //
    //  /*------------------------
    //    Adjoint sensitivity data
    //    ------------------------*/
    //
    //  booleantype cv_adj;             /* TRUE if performing ASA                   */
        mem_dsc.cv_adj              = mem_src.cv_adj;
//
////  struct CVadjMemRec *cv_adj_mem; /* Pointer to adjoint memory structure      */
//TODO:Not needed and difficult    mem_dsc.cv_adj_mem = mem_src.cv_adj_mem;
    //
    //  booleantype cv_adjMallocDone;
        mem_dsc.cv_adjMallocDone    = mem_src.cv_adjMallocDone;
}

