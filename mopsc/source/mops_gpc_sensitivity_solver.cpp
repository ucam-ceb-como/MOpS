/*
  Author(s):      Weerapong Phadungsukana (wp214)
  Project:        mopsc (gas-phase chemistry solver).
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    Implementation of the SensitivityAnalyzer class declared in the
    mops_gpc_sensitivity.h header file.

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
#include "mops_gpc_sensitivity_solver.h"
    
    //realtype cv_uround;      /* machine unit roundoff                        */   

    /*-------------------------- 
      Problem Specification Data 
      --------------------------*/

    //CVRhsFn cv_f;            /* y' = f(t,y(t))                               */
    //void *cv_f_data;         /* user pointer passed to f                     */
    //int cv_lmm;              /* lmm = ADAMS or BDF                           */
    //int cv_iter;             /* iter = FUNCTIONAL or NEWTON                  */
    //int cv_itol;             /* itol = SS or SV                              */

    //realtype cv_reltol;      /* relative tolerance                           */
    //realtype cv_Sabstol;     /* scalar absolute tolerance                    */
    //N_Vector cv_Vabstol;     /* vector absolute tolerance                    */
    //CVEwtFn cv_efun;         /* function to set ewt                          */
    //void *cv_e_data;         /* user pointer passed to efun                  */

    /*-----------------------
      Quadrature Related Data 
      -----------------------*/

    //booleantype cv_quadr;    /* TRUE if integrating quadratures              */

    //CVQuadRhsFn cv_fQ;
    //void *cv_fQ_data;        /* user pointer passed to fQ                    */
    //int cv_itolQ;
    //booleantype cv_errconQ;

    //realtype cv_reltolQ;     /* relative tolerance for quadratures           */
    //realtype cv_SabstolQ;    /* scalar absolute tolerance for quadratures    */
    //N_Vector cv_VabstolQ;    /* vector absolute tolerance for quadratures    */

    /*------------------------
      Sensitivity Related Data 
      ------------------------*/

    //booleantype cv_sensi;    /* TRUE if computing sensitivities              */

    //int cv_Ns;               /* Number of sensitivities                      */

    //int cv_ism;              /* ism = SIMULTANEOUS or STAGGERED              */

    //CVSensRhsFn cv_fS;       /* fS = (df/dy)*yS + (df/dp)                    */
    //CVSensRhs1Fn cv_fS1;     /* fS1 = (df/dy)*yS_i + (df/dp)                 */
    //void *cv_user_fS_data;   /* user data pointer for fS                     */
    //void *cv_fS_data;        /* actual data pointer passed to fS             */
    //booleantype cv_fSDQ;
    //int cv_ifS;              /* ifS = ALLSENS or ONESENS                     */

    //realtype *cv_p;          /* parameters in f(t,y,p)                       */
    //realtype *cv_pbar;       /* scale factors for parameters                 */
    //int *cv_plist;           /* list of sensitivities                        */
    //int cv_DQtype;           /* central/forward finite differences           */
    //realtype cv_DQrhomax;    /* cut-off value for separate/simultaneous FD   */

    //booleantype cv_errconS;  /* TRUE if sensitivities are in err. control    */

    //int cv_itolS;
    //realtype cv_reltolS;     /* relative tolerance for sensitivities         */
    //realtype *cv_SabstolS;   /* scalar absolute tolerances for sensi.        */
//=>    N_Vector *cv_VabstolS;   /* vector absolute tolerances for sensi.        */

    /*-----------------------
      Nordsieck History Array 
      -----------------------*/

    //N_Vector cv_zn[L_MAX];   /* Nordsieck array, of size N x (q+1).
    //                            zn[j] is a vector of length N (j=0,...,q)
    //                            zn[j] = [1/factorial(j)] * h^j * 
    //                            (jth derivative of the interpolating 
    //                            polynomial                                   */

    /*-------------------
      Vectors of length N 
      -------------------*/

    //N_Vector cv_ewt;         /* error weight vector                          */
    //N_Vector cv_y;           /* y is used as temporary storage by the solver.
    //                            The memory is provided by the user to CVode 
    //                            where the vector is named yout.              */
    //N_Vector cv_acor;        /* In the context of the solution of the
    //                            nonlinear equation, acor = y_n(m) - y_n(0).
    //                            On return, this vector is scaled to give
    //                            the estimated local error in y.              */
    //N_Vector cv_tempv;       /* temporary storage vector                     */
    //N_Vector cv_ftemp;       /* temporary storage vector                     */

    /*--------------------------
      Quadrature Related Vectors 
      --------------------------*/

//=>    N_Vector cv_znQ[L_MAX];  /* Nordsieck arrays for quadratures             */
//=>    N_Vector cv_ewtQ;        /* error weight vector for quadratures          */
//=>    N_Vector cv_yQ;          /* Unlike y, yQ is not allocated by the user    */
//=>    N_Vector cv_acorQ;       /* acorQ = yQ_n(m) - yQ_n(0)                    */
//=>    N_Vector cv_tempvQ;      /* temporary storage vector (~ tempv)           */

    /*---------------------------
      Sensitivity Related Vectors 
      ---------------------------*/

//=>    N_Vector *cv_znS[L_MAX]; /* Nordsieck arrays for sensitivities           */
//=>    N_Vector *cv_ewtS;       /* error weight vectors for sensitivities       */
//=>    N_Vector *cv_yS;         /* yS=yS0 (allocated by the user)               */
//=>    N_Vector *cv_acorS;      /* acorS = yS_n(m) - yS_n(0)                    */
//=>    N_Vector *cv_tempvS;     /* temporary storage vector (~ tempv)           */
//=>    N_Vector *cv_ftempS;     /* temporary storage vector (~ ftemp)           */

    /*-----------------------------------------------
      Does CVodeSensMalloc allocate additional space?
      -----------------------------------------------*/  

    //booleantype cv_stgr1alloc; /* Did we allocate ncfS1, ncfnS1, and nniS1?  */

    /*-----------------
      Tstop information
      -----------------*/
    //booleantype cv_istop;
    //booleantype cv_tstopset;
    //realtype cv_tstop;

    /*---------
      Step Data 
      ---------*/

    //int cv_q;                    /* current order                            */
    //int cv_qprime;               /* order to be used on the next step        */ 
    ///* = q-1, q, or q+1                         */
    //int cv_next_q;               /* order to be used on the next step        */
    //int cv_qwait;                /* number of internal steps to wait before  */
    ///* considering a change in q                */
    //int cv_L;                    /* L = q + 1                                */

    //realtype cv_hin;
    //realtype cv_h;               /* current step size                        */
    //realtype cv_hprime;          /* step size to be used on the next step    */ 
    //realtype cv_next_h;          /* step size to be used on the next step    */ 
    //realtype cv_eta;             /* eta = hprime / h                         */
    //realtype cv_hscale;          /* value of h used in zn                    */
    //realtype cv_tn;              /* current internal value of t              */
    //realtype cv_tretlast;        /* last value of t returned                 */

    //realtype cv_tau[L_MAX+1];    /* array of previous q+1 successful step
    //                                sizes indexed from 1 to q+1              */
    //realtype cv_tq[NUM_TESTS+1]; /* array of test quantities indexed from
    //                                1 to NUM_TESTS(=5)                       */
    //realtype cv_l[L_MAX];        /* coefficients of l(x) (degree q poly)     */

    //realtype cv_rl1;             /* the scalar 1/l[1]                        */
    //realtype cv_gamma;           /* gamma = h * rl1                          */
    //realtype cv_gammap;          /* gamma at the last setup call             */
    //realtype cv_gamrat;          /* gamma / gammap                           */

    //realtype cv_crate;           /* est. corrector conv. rate in Nls         */
    //realtype cv_crateS;          /* est. corrector conv. rate in NlsStgr     */
    //realtype cv_acnrm;           /* | acor |                                 */
    //realtype cv_acnrmS;          /* | acorS |                                */
    //realtype cv_acnrmQ;          /* | acorQ |                                */
    //realtype cv_nlscoef;         /* coeficient in nonlinear convergence test */
    //int  cv_mnewt;               /* Newton iteration counter                 */
    //int  *cv_ncfS1;              /* Array of Ns local counters for conv.  
    //                                failures (used in CVStep for STAGGERED1) */

    /*------
      Limits 
      ------*/

    //int cv_qmax;             /* q <= qmax                                    */
    //long int cv_mxstep;      /* maximum number of internal steps for one 
    //                            user call                                    */
    //int cv_maxcor;           /* maximum number of corrector iterations for 
    //                            the solution of the nonlinear equation       */
    //int cv_maxcorS;
    //int cv_mxhnil;           /* maximum number of warning messages issued to 
    //                            the user that t + h == t for the next 
    //                            internal step                                */
    //int cv_maxnef;           /* maximum number of error test failures        */
    //int cv_maxncf;           /* maximum number of nonlinear conv. failures   */

    //realtype cv_hmin;        /* |h| >= hmin                                  */
    //realtype cv_hmax_inv;    /* |h| <= 1/hmax_inv                            */
    //realtype cv_etamax;      /* eta <= etamax                                */

    /*----------
      Counters 
      ----------*/

    //long int cv_nst;         /* number of internal steps taken               */
    //long int cv_nfe;         /* number of f calls                            */
    //long int cv_nfSe;        /* number of fS calls                           */
    //long int cv_nfQe;        /* number of fQ calls                           */
    //long int cv_nfeS;        /* number of f calls from sensi DQ              */

    //long int cv_ncfn;        /* number of corrector convergence failures     */
    //long int cv_ncfnS;       /* number of total sensi. corr. conv. failures  */
    //long int *cv_ncfnS1;     /* number of sensi. corrector conv. failures    */

    //long int cv_nni;         /* number of nonlinear iterations performed     */
    //long int cv_nniS;        /* number of total sensi. nonlinear iterations  */
    //long int *cv_nniS1;      /* number of sensi. nonlinear iterations        */

    //long int cv_netf;        /* number of error test failures                */
    //long int cv_netfS;       /* number of sensi. error test failures         */
    //long int cv_netfQ;       /* number of quadr. error test failures         */

    //long int cv_nsetups;     /* number of setup calls                        */
    //long int cv_nsetupsS;    /* number of setup calls due to sensitivities   */

    //int cv_nhnil;            /* number of messages issued to the user that
    //                            t + h == t for the next iternal step    */

    /*-----------------------------
      Space requirements for CVODES 
      -----------------------------*/

    //long int cv_lrw1;        /* no. of realtype words in 1 N_Vector y        */ 
    //long int cv_liw1;        /* no. of integer words in 1 N_Vector y         */ 
    //long int cv_lrw1Q;       /* no. of realtype words in 1 N_Vector yQ       */ 
    //long int cv_liw1Q;       /* no. of integer words in 1 N_Vector yQ        */ 
    //long int cv_lrw;         /* no. of realtype words in CVODES work vectors */
    //long int cv_liw;         /* no. of integer words in CVODES work vectors  */

    /*----------------
      Step size ratios
      ----------------*/

    //realtype cv_etaqm1;      /* ratio of new to old h for order q-1          */
    //realtype cv_etaq;        /* ratio of new to old h for order q            */
    //realtype cv_etaqp1;      /* ratio of new to old h for order q+1          */

    /*------------------
      Linear Solver Data 
      ------------------*/

    /* Linear Solver functions to be called */

    //int (*cv_linit)(struct CVodeMemRec *cv_mem);

    //int (*cv_lsetup)(struct CVodeMemRec *cv_mem, int convfail, 
    //                 N_Vector ypred, N_Vector fpred, booleantype *jcurPtr, 
    //                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3); 

    //int (*cv_lsolve)(struct CVodeMemRec *cv_mem, N_Vector b, N_Vector weight,
    //                 N_Vector ycur, N_Vector fcur);

    //void (*cv_lfree)(struct CVodeMemRec *cv_mem);

    /* Linear Solver specific memory */

    //void *cv_lmem;           

    /* Flag to request a call to the setup routine */

    //booleantype cv_forceSetup;

    /*------------
      Saved Values
      ------------*/

    //int cv_qu;             /* last successful q value used                        */
    //long int cv_nstlp;     /* step number of last setup call                      */
    //realtype cv_h0u;       /* actual initial stepsize                             */
    //realtype cv_hu;        /* last successful h value used                        */
    //realtype cv_saved_tq5; /* saved value of tq[5]                                */
    //booleantype cv_jcur;   /* Is the Jacobian info used by linear solver current? */
    //realtype cv_tolsf;     /* tolerance scale factor                              */
    //int cv_qmax_alloc;     /* value of qmax used when allocating memory           */
    //int cv_qmax_allocQ;    /* value of qmax used when allocating quad. memory     */
    //int cv_qmax_allocS;    /* value of qmax used when allocating sensi. memory    */
    //int cv_indx_acor;      /* index of the zn vector in which acor is saved       */
    //booleantype cv_setupNonNull; /* Does setup do something?                      */

    /*--------------------------------------------------------------------
      Flags turned ON by CVodeMalloc, CVodeSensMalloc, and CVodeQuadMalloc 
      and read by CVodeReInit, CVodeSensReInit, and CVodeQuadReInit
      --------------------------------------------------------------------*/

    //booleantype cv_VabstolMallocDone;
    //booleantype cv_MallocDone;

    //booleantype cv_VabstolQMallocDone;
    //booleantype cv_quadMallocDone;

    //booleantype cv_VabstolSMallocDone;
    //booleantype cv_SabstolSMallocDone;
    //booleantype cv_sensMallocDone;

    /*-------------------------------------------
      Error handler function and error ouput file 
      -------------------------------------------*/

    //CVErrHandlerFn cv_ehfun;    /* Error messages are handled by ehfun       */
    //void *cv_eh_data;           /* user pointer passed to ehfun              */
    //FILE *cv_errfp;             /* CVODES error messages are sent to errfp   */    

    /*-------------------------
      Stability Limit Detection
      -------------------------*/

    //booleantype cv_sldeton;     /* Is Stability Limit Detection on?          */
    //realtype cv_ssdat[6][4];    /* scaled data array for STALD               */
    //int cv_nscon;               /* counter for STALD method                  */
    //long int cv_nor;            /* counter for number of order reductions    */

    /*----------------
      Rootfinding Data
      ----------------*/

    //CVRootFn cv_gfun;     /* Function g for roots sought                     */
    //int cv_nrtfn;         /* number of components of g                       */
    //void *cv_g_data;      /* pointer to user data for g                      */
    //int *cv_iroots;       /* int array for root information                  */
    //realtype cv_tlo;      /* nearest endpoint of interval in root search     */
    //realtype cv_thi;      /* farthest endpoint of interval in root search    */
    //realtype cv_trout;    /* t value returned by rootfinding routine         */
    //realtype *cv_glo;     /* saved array of g values at t = tlo              */
    //realtype *cv_ghi;     /* saved array of g values at t = thi              */
    //realtype *cv_grout;   /* array of g values at t = trout                  */
    //realtype cv_toutc;    /* copy of tout (if NORMAL mode)                   */
    //realtype cv_ttol;     /* tolerance on root location trout                */
    //int cv_taskc;         /* copy of parameter task                          */
    //int cv_irfnd;         /* flag showing whether last step had a root       */
    //long int cv_nge;      /* counter for g evaluations                       */

