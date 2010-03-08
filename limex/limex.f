

      program BETTERSOLVER

      call XLIMEX(A,B,C,D,E,F,G,H,I,J,K,L,M,N)

	end


      
      subroutine XLIMEX ( n, Fcn, Jacobian, t_Begin, t_End, y, ys, rTol,
	2                    aTol, h, Iopt, Ropt, IPos, IFail )

      implicit none

      integer MAXNEQ,LRW,LIW
      parameter(MAXNEQ=500)
      parameter(LRW=50000,LIW=50000)


c
   
c
c-----------------------------------------------------------------------
c
c     Extrapolation integrator  for the solution of  linearly-implicit
c     differential-algebraic systems of the form 
c
c          B (t,y) * y' (t) = f (t,y) 
c                                
c     with B a (n,n)-matrix of rank less or equal n.
c
c-----------------------------------------------------------------------
c  
c     Copyright (C) 2000, Konrad-Zuse-Zentrum fuer Informationstechnik
c     Berlin (ZIB)
c     ALL RIGHTS RESERVED
c
c     Written by:
c
c     R. Ehrig, U. Nowak,
c     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
c     Takustrasse 7 
c     D-14195 Berlin-Dahlem
c
c     phone : +49-30-84185-0
c     fax   : +49-30-84185-125
c     e-mail: ehrig@zib.de, nowak@zib.de
c     URL   : http://www.zib.de
c
c-----------------------------------------------------------------------
c
c            **************************************************
c            **                                              **
c            **   This is version 4.2A1 of March, 17, 2000   **
c            **                                              **
c            **************************************************
c
c-----------------------------------------------------------------------
c
c     Overview of current versions:
c
c     4.2A1  Non-sparse dense  or  banded  Jacobians. Direct  solvers:
c            LAPACK routines. Based on BLAS and LAPACK routines.
c
c     4.2A2  Non-sparse dense or banded Jacobians. Direct solvers: NAG
c            routines. Based on NAG routines (Mark 16).
c
c     4.2B1  Sparse Jacobians. Direct  solvers: MA28 routines from the
c            Harwell subroutine library, iterative  solvers: GMRES and
c            BICGSTAB  with a  variable ILU  preconditioner. Based  on
c            BLAS and LAPACK routines.
c
c     4.2B2  Sparse Jacobians. Direct solvers: NAG routines, iterative
c            solvers: GMRES  and  BICGSTAB  with a  variable ILU  pre-
c            conditioner. Based on NAG routines (Mark 16).
c
c     Versions with other solver routines are available on request.
c                             
c-----------------------------------------------------------------------
c
c     NOTICE: "The LIMEX  program may be used SOLELY  for educational, 
c     research, and benchmarking  purposes by no-profit organizations.
c     Commercial and other organizations  may make use of LIMEX SOLELY
c     for benchmarking  purposes only. LIMEX may  be modified by or on
c     behalf of the user for  such use but  at no time  shall LIMEX or 
c     any such  modified version  of LIMEX become the  property of the
c     user. LIMEX  is provided  without warranty  of any  kind, either
c     expressed  or implied.  Neither the  authors nor their employers
c     shall be liable  for any direct or consequential  loss or damage
c     whatsoever  arising out  of the  use or  misuse of  LIMEX by the
c     user. LIMEX  must not  be sold. You  may make  copies  of LIMEX, 
c     but  this NOTICE  and the  Copyright notice  must appear  in all 
c     copies. Any  other  use  of LIMEX  requires written  permission.
c     Your use of LIMEX is an implicit agreement to these conditions."
c
c     LIMEX is also  available for commercial use. Contact the authors 
c     who will provide details of price and condition of use.
c
c-----------------------------------------------------------------------
c
c     Description
c
c     LIMEX  solves linearly-implicit  differential-algebraic  systems 
c     (DAEs) of the form 
c
c        B (t,y) * y' (t) = f (t,y).
c 
c     If n is the size of the system, B is a (n,n)-matrix of rank less 
c     or equal n. General  conditions for  the applicability  of LIMEX
c     are a regular  matrix pencil B + h A with A the Jacobian of  the 
c     residual of the DAE and an index of the DAE less or equal 1. The
c     discretization  of LIMEX  is based  on the  elementary  linearly 
c     implicit Euler discretization
c
c        ( B(t,y(k)) - h J ) ( y(k+1) - y(k) ) = h f(t(k+1),y(k))
c
c     where J is the (approximate) Jacobian of the residual
c
c         d              |
c        -- ( f - B y' ) |      .
c        dy              |t=t(0)
c
c     Combined  with  extrapolation  this one-step  method permits  an
c     adaptive control of stepsize and order. Within the extrapolation
c     process  one  computes  for  a  basic stepsize H  approximations
c     T(j,1) for y(t(0)+H) using  the discretization  above with step-
c     sizes  h(j) = H / j, j = 1, ..., j(max). Then the  extrapolation
c     tableau recursively defines higher order approximations T(j,k):
c
c                            T(j,k-1) - T(j-1,k-1)
c        T(j,k) = T(j,k-1) + ---------------------, k = 1, ..., j.
c                            j / ( j - k + 1 ) - 1
c
c     As error estimates the subdiagonal differences T(j,j) - T(j,j-1)
c     are taken. 
c
c     The efficiency of LIMEX mainly depends on the performance of the 
c     evaluation of the Jacobian  and in particular on the solution of
c     the  linear systems. Therefore  different versions  of LIMEX are
c     available with  different methods to create the  Jacobian and to
c     solve the linear systems.
c
c     Jacobian generation: numerical or analytical (user supplied)
c
c     Linear systems solvers: direct for dense matrices, LIMEX version
c     4.2A1 and 4.2A2, direct or iterative methods (GMRES or BICGSTAB) 
c     for  sparse  matrices, LIMEX  version 4.2B1  and 4.2B2. The code
c     for dense matrices determines the optimal matrix representation,
c     full or banded, from the  ratio between the system size  and the 
c     number of lower and upper diagonals in the Jacobian.
c
c     Throughout  the code a local error control is implemented, which
c     requires that the  local error of a component y(i) is  less than
c     rTol * abs ( y(i) ) + aTol(i). This approach enables  to specify
c     more or less sensitive components of the solution vector.
c
c     The code has an additional option to compute  consistent initial
c     values (CIVs). For this  task an  extrapolated Newton  iteration
c     is used to  determine correct  start values at t = t(0) for  the
c     algebraic variables (which are known only  approximately in many
c     applications) and the derivatives  of the differential variables 
c     (which are unknown in most applications), which  satisfy exactly
c     the  DAE. The algebraic  variables are these  components of  the 
c     solution  vector, for  which  the  corresponding columns  in the
c     matrix B(t,y) contain zeros only. All other  components will  be
c     accounted  as  differential  variables,  whose  values  are  not 
c     changed during the CIV computation. However, in many cases where 
c     only the  right derivatives at  the start of the integration are 
c     unknown, LIMEX solves the DAE even without CIV determination.
c
c     Furthermore a dense output option is availabe. With this  option
c     one may compute solutions of the DAE on a dense set of points in
c     the integration interval with nearly the same accuracy as at the
c     automatically  selected integration  points. The method bases on
c     Hermitian interpolation and considers all intermediate solutions
c     computed in  the extrapolation process. The dense  output option
c     may be combined with the  one-step mode, i.e. the return  to the 
c     calling  program after  every  integration  step, and  also with 
c     continuation  calls. These options  enable to obtain efficiently 
c     solutions  on an  arbitrary user defined  set of values  for the
c     independent variable.
c
c     If information  about the  structure of the  Jacobian  matrix is  
c     needed, the  matrix in a specific  integration step is available
c     in a PostScript plot.
c
c     Since in many  applications some components of the solution have
c     to be nonnegative, due to their physical interpretation, one may
c     set for these  components that only positive values ( => 0 ) are 
c     allowed. If in the extrapolation negative values  are occurring,
c     the stepsize is reduced by a heuristic factor.
c
c-----------------------------------------------------------------------
c
c     References:
c
c     1) P. Deuflhard, U. Nowak:
c        Extrapolation integrators for quasilinear implicit ODEs.
c        In: P. Deuflhard, B. Engquist (eds.): Large  scale scientific
c        computing, Birkhaeuser, Prog. Sci. Comp. 7, pp. 37-50 (1987)
c
c     2) P. Deuflhard, E. Hairer, J. Zugck:
c        One step and extrapolation methods for differential-algebraic
c        systems.
c        Num. Math. 51, pp. 501-516 (1987)
c
c     3) P. Deuflhard:
c        Recent progress in extrapolation methods for ODEs.
c        SIAM Rev. 27, pp. 505-535 (1985)
c
c     4) R. Ehrig, U. Nowak, L. Oeverdieck, P. Deuflhard:
c        Advanced extrapolation  methods for large  scale differential
c        algebraic problems.
c        In: High  Performance Scientific  and Engineering  Computing,
c        H.-J. Bungartz,  F. Durst, and  Chr. Zenger  (eds.),  Lecture
c        Notes  in  Computational  Science  and  Engineering, Springer
c        Vol 8, pp. 233-244 (1999)
c
c-----------------------------------------------------------------------
c
c     Installation notes:
c
c     Requires the BLAS (Basic Linear Algebra Subprograms) library and
c     the LAPACK  library. Ideally, one  should  use  vendor-optimized
c     BLAS and LAPACK routines. If these are not available, one should
c     use the routines collected in source form in the file 
c
c        'LIMEX4_2_Auxiliaries.f'
c
c     This  file is  part  of the  LIMEX4_2A1 distribution. The  whole
c     Fortran BLAS and LAPACK libraries are  available from the netlib
c     repository  at http://cm.bell-labs.com/netlib/master/readme.html
c     or one of its mirrors.
c
c     This code is written almost  purely in ANSI Fortran 77 standard,
c     only the following non-standard features are used:
c
c        do ...  end do  constructions (without statement labels)
c
c        implicit none   to use only explicitly declared variables
c
c        include '....'  for inserting code  
c
c        symbolic names may be longer than 6 characters
c
c        lower case and underscore characters used
c
c     All dimension statements in LIMEX  refer to the size definitions
c     in the include file 'LIMEX4_2_Size_Definitions.h'. This file has
c     to consist of four parameter statements of the following form:
c
c-----------------------------------------------------------------------
c
c     Define the  maximum number of  equations. This statement sets an 
c     upper limit for the number of equations defining the DAE. 
c
c     parameter ( Max_Nr_of_Equations    =
c
c-----------------------------------------------------------------------
c
c     Define the maximum number  of non-zero entries in the  Jacobian.
c     This statement defines an upper limit for the number of non-zero
c     entries  in  the  Jacobian of  the  residual  of  the  DAE. This 
c     parameter is not used in  the LIMEX versions 4.2A1 and 4.2A2, it
c     is  included here  only for consistency  of the include files of
c     the different LIMEX versions.
c
c     parameter ( Max_Non_Zeros_Jacobian =
c
c-----------------------------------------------------------------------
c
c     Define the  maximum  number of non-zero entries in the left-hand
c     side matrix B.
c
c     parameter ( Max_Non_Zeros_B       =
c
c-----------------------------------------------------------------------
c
c     Define the  maximum number of  lower diagonals in  the Jacobian,
c     i.e. the lower  half-bandwith of  this matrix. The value of this
c     parameter  must  be  between 0 (no lower  diagonals) and n (full
c     lower triangular submatrix possible). The actual number of upper
c     diagonals is defined  in LIMEX by the control parameter Iopt(8).
c     This  parameter is  not  used  in  the LIMEX  versions 4.2B1 and
c     4.2B2.
c
c     parameter ( Max_Lower_Diagonals    =
c
c-----------------------------------------------------------------------
c
c     Define the  maximum number of  upper diagonals in  the Jacobian,
c     i.e. the upper  half-bandwith of  this matrix. The value of this
c     parameter  must  be  between 0 (no upper  diagonals) and n (full
c     upper triangular submatrix possible). The actual number of upper
c     diagonals is defined in LIMEX by the control  parameter Iopt(9).
c     This  parameter is  not  used  in  the LIMEX  versions 4.2B1 and
c     4.2B2.
c
c     parameter ( Max_Upper_Diagonals    =
c
c-----------------------------------------------------------------------
c
c     Define the maximum number of vectors which are  available within 
c     GMRES or  BICGSTAB. For GMRES  this statement  implies  an upper 
c     limit for the dimension of the Krylov  subspaces. This dimension 
c     is  defined by  Iopt(20), see  below. GMRES  then requires  that
c     Max_It_Vectors => Iopt(20) + 4 (this  relation  is  checked). If  
c     the Krylov subspace is exhausted without  attaining the required
c     accuracy, a  restart is done. Reasonable values for Iopt(20) are
c     in many cases 10 or 20, then Max_It_Vectors should  be set to 15
c     resp. 25. The greater  the dimension  of the Krylov subspace is, 
c     the  faster GMRES  converges  (in theory), but  needs a  rapidly 
c     increasing number of floating point  operations and is also less
c     stable. If  BICGSTAB is  used, Max_It_Vectors must be 6 at least 
c     (this  condition is  checked). If  no iterative  method will  be
c     used, Max_It_Vectors may  be set equal 1. This parameter  is not
c     used in the LIMEX versions 4.2A1 and 4.2A2, it is  included here 
c     only for consistency of the include files of the different LIMEX 
c     versions.
c
c     parameter ( Max_It_Vectors         =
c
c-----------------------------------------------------------------------
c
c     Arguments
c
c-----------------------------------------------------------------------
c
c     All restrictions  are  checked, but  no further  control of  the
c     consistency of the input data will be performed. 
c
c     n          Integer variable, must be set by caller on input, not 
c                modified. Size of the  differential algebraic system.
c                Restriction: Max_Nr_of_Equations => n => 1.
c 
c     Fcn        Name of an  external subroutine computing  the values 
c                of  the  right-hand  side  function  f(t,y)  and  the 
c                entries in the matrix B (see below).
c
c     Jacobian   Name of an  external subroutine computing  analytical
c                derivatives of the residual f(t,y) - B(t,y) * y' (see
c                below).
c
c     t_Begin    Real variable, must be set by caller on  input. Inde-
c                pendent variable, starting point  of the integration.
c                On exit, the  value  of the  independent value, up to
c                which LIMEX has solved the DAE successfully. 
c
c     t_End      Real variable, must  be set  by caller on  input, not
c                modified. Independent  variable, final  point  of the 
c                integration.
c
c     y          Real  array  of  size  n, must  be set  by caller  on 
c                input. The  initial  values  of the  solution at  the 
c                starting  point  t_Begin. On  exit,  y  contains  the 
c                solution  at  the  final  point t_End, resp. at  this  
c                value of the independent variable, up  to which LIMEX
c                has solved the DAE successfully. 
c
c     ys         Real array of size n. Derivatives  of the solution at
c                t_Begin. Only the values  of ys for  the differential
c                variables  are needed. If the  values for ys are  not 
c                known, set ys to  zero. If only estimates are  known, 
c                these  should be  used. However, the  option  Iopt(6) 
c                for the  computation of consistent initial values can
c                be used to  start always with  correct values for ys. 
c                On exit, ys contains  the derivatives at t_End, resp. 
c                at this  value  of  the  independent  variable, up to 
c                which LIMEX has solved the DAE successfully.
c
c     rTol       Real  array  of size 1 at least for Iopt(11) = 0 or n
c                for Iopt(11) = 1. Values  must be  set by  caller  on 
c                input,  not  modified.  Relative  error   tolerances, 
c                interpreted as a scalar or a vector  depending on the
c                value of Iopt(11). The code keeps the local  error of
c                y(i)  below  rTol(i)*abs(y(i))+aTol(i).  
c
c                Restriction: rTol(*) > 0.
c
c     aTol       Real array  of size 1 at least  for Iopt(11) = 0 or n
c                for  Iopt(11) = 1. Values  must be set by  caller  on
c                input,  not  modified.  Absolute  error   tolerances, 
c                interpreted as a scalar or a vector  depending on the
c                value of Iopt(11). 
c
c                Restriction: aTol(*) => 0. 
c
c                Remark: For any  component y(i), which could  be zero  
c                during the integration  one has to set aTol(i) > 0 to 
c                prevent  divisions  by  zero.  In  many  well  scaled 
c                problems a reasonable setting is aTol(i) = rTol(i).
c
c     h          Real  variable, must  be  set  by  caller  on  input.
c                Initial  stepsize  guess, if h  is  set  to zero, the
c                program puts  h = h_Min. This  variable is  currently 
c                set  in a  parameter statement to 1.d0-4, but  may be
c                changed. On exit h is the  estimated optimal stepsize
c                for the next  integration step. If an  error occurred
c                h is set to zero.
c
c     Iopt       Integer  array  of  size 30, values  from  Iopt(1) to
c                Iopt(18) must be set by caller on  input. Integration 
c                control parameters.
c
c                Iopt(1)  Integration monitoring, not modified 
c                         = 0 : no output
c                         = 1 : standard output
c                         = 2 : additional integration monitor
c
c                         Restriction: 0 <= Iopt(1) <= 2.
c
c                Iopt(2)  The  unit  number  for  monitor  output, not
c                         modified. If Iopt(2) = 0, the  default value
c                         for the unit number is 6.
c
c                         Restriction: 0 <= Iopt(2) if Iopt(1) > 0.
c
c                Iopt(3)  Solution output, not modified
c                         = 0 : no solution output
c                         = 1 : initial and final solution values
c                         = 2 : additional  solution  values at inter-
c                               mediate points
c
c                         Restriction: 0 <= Iopt(3) <= 2.
c
c                Iopt(4)  The  unit  number for  solution  output, not
c                         modified. If Iopt(4) = 0, the  default value
c                         for the unit number is 6.
c
c                         Restriction: 0 <= Iopt(4) if Iopt(3) > 0.
c
c                Iopt(5)  Singular or nonsingular  matrix B, not modi-
c                         fied
c                         = 0 : matrix B may be singular
c                         = 1 : matrix B is nonsingular
c                         For theoretical  reasons Iopt(5) = 0 reduces 
c                         the  maximum  order  to 5, even  if  in most 
c                         cases  LIMEX  also  runs  successfully  with
c                         Iopt(5) = 1  and  singular  matrices B(t,y).
c
c                         Restriction: 0 <= Iopt(5) <= 1.
c
c                Iopt(6)  Determination  of consistent  initial values 
c                         (CIVs), not modified 
c                         = 0 : no determination of CIVs
c                         = 1 : determination of CIVs
c                         If Iopt(6) = 1, the code computes before the
c                         integration  consistent initial  values. For 
c                         this  reason, new  values for  the algebraic 
c                         variables  and  for the  derivatives  of the 
c                         differential  variables are  computed, which
c                         satisfy the DAE up to the required accuracy.
c
c                         Restriction: 0 <= Iopt(6) <= 1.
c
c                Iopt(7)  Numerical  or analytical  generation  of the 
c                         Jacobian matrix, not modified
c                         = 0 : Numerical  difference approximation of
c                               Jacobian of the residual
c                         = 1 : Analytic Jacobian supplied by the user
c                               subroutine Jacobian
c
c                         Restriction: 0 <= Iopt(7) <= 1.
c
c                Iopt(8)  Lower bandwidth of the  Jacobian matrix, not 
c                         modified. If  the Jacobian is a band matrix,  
c                         the computing time needed for the evaluation
c                         of this matrix can be substantially reduced.
c                         If Iopt(8) => n, the whole lower  triangular 
c                         submatrix will be computed. If Iopt(8) (and/
c                         or Iopt(9)) is less  than n, the Jacobian is
c                         a band  matrix. Then LIMEX uses  specialized
c                         subroutines  for the factorization  and  the 
c                         solution of linear systems, which are faster
c                         than  the corresponding  algorithms for full 
c                         matrices. Also the storage needed is reduced
c                         if  the parameters  Max_Upper_Diagonals  and 
c                         Max_Lower_Diagonals are  set properly in the
c                         include file LIMEX4_2_Size_Definitions.h. 
c
c                         Restriction: 
c                         0 <= Iopt(8) <= Max_Lower_Diagonals. 
c
c                Iopt(9)  Upper bandwidth of the  Jacobian matrix, not 
c                         modified. If  the Jacobian is a band matrix, 
c                         the computing time needed for the evaluation
c                         of this matrix can be substantially reduced. 
c                         If Iopt(9) => n, the whole upper  triangular 
c                         submatrix will be computed. If Iopt(9) (and/
c                         or Iopt(8)) is less  than n, the Jacobian is
c                         a band matrix. Read then the hints above for
c                         Iopt(8).
c
c                         Restriction:
c                         0 <= Iopt(9) <= Max_Upper_Diagonals. 
c
c                Iopt(10) Determines whether the code tries to use the
c                         Jacobians in more than one integration step, 
c                         not modified. This can lead to a significant 
c                         increase of the overall performance.
c                         = 0 : no reuse of the Jacobian
c                         = 1 : reuse of the  Jacobian in the follwing
c                               integration steps
c                         If Iopt(10) = 1, the  code  computes the  so
c                         called contractivity factor, an estimate for
c                         the current linearity of  the DAE. Depending 
c                         from an  upper  bound for  this  factor, the
c                         current Jacobian  is used also  in the  next 
c                         integration step. The upper bound is defined
c                         by ThMin. This  variable  is set  within the
c                         parameter statement  below and should always
c                         selected between 0 (no  reuse the  Jacobian)
c                         and 1 (try almost in every step to reuse the
c                         Jacobian). A recommended value for the first
c                         computations  is 2.0d-2. If  in the  current 
c                         step a  new Jacobian  is computed, the  step
c                         number is marked in the monitor  output with
c                         an asterisk '*'. If  the stepsize is reduced
c                         due to  an unsatisfying  convergence  within 
c                         the extrapolation, the true current Jacobian  
c                         is recomputed if necessary.
c
c                         Restriction: 0 <= Iopt(10) <= 1.
c
c                Iopt(11) Switch  for error  tolerances, not modified.
c                         = 0 : both rTol and aTol are scalars
c                         = 1 : both rTol and aTol are vectors 
c
c                         Restriction: 0 <= Iopt(11) <= 1.
c
c                Iopt(12) Switch for the one step mode, not  modified.
c                         If Iopt(12) is equal 1, LIMEX returns during
c                         the integration to the calling program after
c                         each integration  step. If  Iopt(12) = 2 and
c                         the  dense  output  option is  active , i.e. 
c                         Iopt(13) > 0, the code  returns only  on the
c                         specified dense  output points (see  below). 
c                         If LIMEX returns during the integration, the
c                         following  parameters of  LIMEX contain  the 
c                         current  values and  may be used  within the 
c                         calling program:
c
c                         t_Begin: the current value of t (independent
c                                  variable)
c                         y      : the current solution 
c
c                         To  assure  a correct  continuation  of  the
c                         integration, no parameter of LIMEX should be
c                         changed  before the  next call  of LIMEX. If
c                         Iopt(12) = 0, LIMEX does  not return  to the
c                         calling program  until the whole integration
c                         is finished.
c
c                         Restriction: 0 <= Iopt(12) <= 2.
c
c                Iopt(13) Dense output option, not modified. 
c                         = 0 : no dense output
c                         = 1 : dense  output  on  equidistant  points
c                               within the whole integration interval.
c                               The number of such points is specified
c                               by Iopt(14).
c                         = 2 : dense  output  on  equidistant  points
c                               within  every  integration  step.  The 
c                               number of such points  is specified by
c                               Iopt(14).
c                         = 3 : output  of the solution  at the end of
c                               every integration  interval and  dense
c                               output on some additional points, such
c                               the distance  of two  output points is
c                               always  less or  equal t_Max. tMax  is
c                               specified by Ropt(2).
c
c                         If  one sets  Iopt(13) > 0 and Iopt(12) = 2, 
c                         LIMEX returns to  the calling program at all
c                         dense output points.
c
c                         The solution at such points is  computed via 
c                         Hermite interpolation. The size of the error 
c                         of the interpolated solution may be slightly 
c                         larger than  specified by  the rTol and aTol
c                         values. 
c
c                         Restriction: 0 <= Iopt(13) <= 3.
c
c                Iopt(14) The number of  equidistant output points, if 
c                         the dense output  option  Iopt(13) = 1 or 2, 
c                         not  modified. If  Iopt(13) = 1, Iopt(14) is
c                         the number of output points within the whole
c                         integration  interval. If Iopt(14) = 0 or 1,
c                         no dense output will be provided.
c
c                         Example: if  t_Begin = 0 and  t_End = 1 with
c                         Iopt(13) = 1 and  Iopt(14) = 6  the solution 
c                         at  the  points 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 
c                         will be provided.
c
c                         If Iopt(13) = 2, Iopt(14) is  the number  of
c                         output points within every integration step.
c                         If  Iopt(14) = 0, no  dense  output will  be
c                         provided.  If  Iopt(14) = 1 or 2, the  dense 
c                         output consists of the solutions at the ends
c                         of the integration intervals.
c
c                         Example: if the current integration interval 
c                         is [0,1], with Iopt(13) = 2 and Iopt(14) = 6
c                         the solution at 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
c                         will be computed in this interval.
c
c                         Restriction: 0 <= Iopt(14)  if  Iopt(13) = 1 
c                                      or 2.    
c
c                Iopt(15) The unit number  for dense solution  output,
c                         not  modified.  Iopt(15) = 0  supresses  the
c                         dense output, which is then available in the
c                         calling program by setting Iopt(12) = 2.
c
c                         Restriction: 0 <= Iopt(15) if Iopt(13) > 0.
c
c                Iopt(16) Type  of call, may be modified. Iopt(16) = 0
c                         specifies  the  current  call  as an initial
c                         call, i.e. the  first  call  of LIMEX  for a 
c                         given  problem. Thus  LIMEX performs  then a 
c                         full data check and initialization. Further-
c                         more the solution  output  and, if required, 
c                         the dense  output will  be printed  also  at 
c                         t_Begin. Iopt(16) = 1 specifies a successive
c                         call, i.e. the integration of the last LIMEX
c                         call will  be continued. Then  the  solution
c                         and/or the dense  output will not be printed
c                         at t_Begin. To assure a correct continuation
c                         of the  further integration, no parameter of 
c                         LIMEX except t_End should be changed between
c                         the LIMEX calls. If the code returns with no
c                         error, Iopt(16) on return will be always set
c                         to 1, thus enabling the continuation  of the 
c                         integration.
c
c                         Restriction: 0 <= Iopt(16) <= 1.
c
c                Iopt(17) This option determines the behavior of LIMEX
c                         at  t_End, not  modified. With  Iopt(17) = 0
c                         the code integrates  exactly up to t_End. If
c                         Iopt(17) = 1, then  LIMEX may  have  already
c                         computed internally solutions for t > t_End.
c                         This option enables to  get most efficiently
c                         solutions  on an arbitrary user defined  set
c                         of values of the independent  variable t. If
c                         f(t,y)  or B(t,y) is  defined only  until an
c                         upper  limit, Ropt(3) may  be used to assure
c                         that t is always less or equal this limit.
c
c                         Restriction: 0 <= Iopt(17) <= 1.
c
c                Iopt(18) By  means of  this option one may generate a 
c                         PostScript plot of the  Jacobian  matrix. If
c                         Iopt(18) = i > 0, in the  integration step i
c                         the  file 'Jacobian.ps'  will  be  generated
c                         representing the  non-zero structure  of the
c                         current  Jacobian. If in  this  step no  new 
c                         Jacobian is computed, since Iopt(10) = 1, no
c                         PostScript plot will be  generated. The sub-
c                         routine JacPlot internally  uses the logical 
c                         unit-number 99. If Iopt(18) = 0 no plot will
c                         be generated, if Iopt(18) = -1 the structure
c                         of the initial Jacobian (step 0) is plotted.
c
c                         Restriction: -1 <= Iopt(18).
c           
c                Iopt(19) - Iopt(23)
c                         Parameters used  only in the  LIMEX versions
c                         4.2B1 and 4.2B2 (sparse Jacobians).
c
c                Iopt(24) On return the number of function evaluations
c                         during the integration without those for the
c                         computation of the Jacobian.
c
c                Iopt(25) On return the number of function evaluations
c                         needed for the  numerical computation of the 
c                         Jacobian.
c
c                Iopt(26) On return the number of LU decompositions by
c                         the LAPACK routines dgetrf or dgbtrf.
c
c                Iopt(27) On return  the number of  back-substitutions
c                         by the LAPACK routines dgetrs or dgbtrs.
c
c                Iopt(28) On return the number of integration steps.
c
c                Iopt(29) On  return  the  number  of Jacobian  matrix 
c                         evaluations.
c
c                Iopt(30) Parameter  used  only in the  LIMEX versions
c                         4.2B1 and 4.2B2 (sparse Jacobians).
c
c     Ropt       Real array of size 5, Ropt(1) to Ropt(3) must  be set 
c                by caller on input. Integration control parameter.
c
c                Ropt(1)  The maximum allowed  stepsize, not modified.
c                         If Ropt(1) = 0, the maximum allowed stepsize
c                         is set to the default value t_End - t_Begin.
c
c                         Restriction: 0 <= Ropt(1).
c
c                Ropt(2)  The  maximal  distance of  two dense  output 
c                         points, not  modified. Used only if Iopt(13) 
c                         is equal 3.
c
c                         Restriction: 0 < Ropt(2) if Iopt(13) = 3.
c
c                Ropt(3)  The upper limit for the independent variable
c                         t, not modified. Used only  if Iopt(17) = 1. 
c                         If Ropt(3)  is less  than t_End the value of
c                         Ropt(3) will be ignored  and no upper  limit 
c                         will be considered.
c                         
c                Ropt(4), Ropt(5)
c                         Parameters used  only in the  LIMEX versions
c                         4.2B1 and 4.2B2 (sparse Jacobians).
c
c     IPos       Integer array of size n, values must be set by caller
c                on input, not modified. If IPos(i) = 1, the values of
c                the corresponding  solution component will be checked 
c                to be  nonnegative. If  during  the  extrapolation  a
c                negative value occurs, the stepsize is reduced by the
c                heuristic factor RedPos, which  is set in a parameter
c                statement to 0.5, but  may  be changed. If IPos(i) is 
c                not equal 1, no check will be performed.
c
c     IFail      Integer array of size 3, need not to be set by caller
c                on input. Error indicator field.
c
c                IFail(1) Internal LIMEX error code. Zero if  no error
c                         occurred. If  Iopt(1) < 0, for each error an
c                         explanatory message is output on the current
c                         monitor  unit. The error  codes -1 until -31  
c                         indicate  wrong or  inconsistent input data,
c                         the other error codes indicate errors during
c                         the integration.
c
c                Error-code IFail(1)            Description
c
c                         -1            Max_Nr_of_Equations  <  n   or  
c                                       n < 1 
c
c                         -2            rTol(i) <= 0 for some index i
c
c                         -3            aTol(i) < 0 for some index i
c
c                         -4            Iopt(1) < 0 or Iopt(1) > 2
c
c                         -5            Iopt(2) < 0 and Iopt(1) > 0
c
c                         -6            Iopt(3) < 0 or Iopt(3) > 2
c
c                         -7            Iopt(4) < 0 and Iopt(3) > 0
c
c                         -8            Iopt(5) < 0 or Iopt(5) > 1
c
c                         -9            Iopt(6) < 0 or Iopt(6) > 1
c
c                        -10            Iopt(7) < 0 or Iopt(7) > 1
c
c                        -11            Max_Lower_Diagonals <  Iopt(8)
c                                       or Iopt(8) < 0
c
c                        -12            Max_Upper_Diagonals <  Iopt(9)
c                                       or Iopt(9) < 0
c
c                        -13            Iopt(10) < 0 or Iopt(10) > 1
c
c                        -14            Iopt(11) < 0 or Iopt(11) > 1
c
c                        -15            Iopt(12) < 0 or Iopt(12) > 1
c
c                        -16            Iopt(13) < 0 or Iopt(13) > 3
c
c                        -17            Iopt(14) < 0 and  Iopt(13) = 1
c                                       or 2
c
c                        -18            Iopt(15) < 0 and Iopt(13) > 0
c
c                        -19            Iopt(16) < 0  or  Iopt(16) > 1 
c
c                        -20            Iopt(17) < 0  or  Iopt(17) > 1 
c
c                        -21            Iopt(18) < - 1
c
c                        -27            Ropt(1) < 0
c
c                        -28            Ropt(2) <= 0 and Iopt(13) = 3
c
c                        -32            An  initial value of y is less
c                                       than zero, but  by the setting 
c                                       of IPos it should be => 0.
c
c                        -33            The routine Fcn  returns  with
c                                       an   error ( FcnInfo < 0 )  in
c                                       the first call  of Fcn in  the
c                                       current integration interval.
c                
c                        -34            More than Max_Non_Zeros_B non-
c                                       zero entries in  matrix B(t,y) 
c                                       defined.
c
c                        -35            The routine Fcn  returns  with
c                                       an  error (FcnInfo < 0) in the
c                                       numerical  evaluation  of  the
c                                       Jacobian.
c
c                        -36            The  routine Jacobian  returns 
c                                       with an error (JacInfo < 0).
c
c                        -39            Internal error  in the  LAPACK
c                                       routines dgetrf or dgbtrf, see  
c                                       IFail(2) and the LAPACK User's 
c                                       Guide.
c
c                        -43            Problem not solvable with this
c                                       LIMEX version,  most  probable
c                                       reason: index of the DAE > 1.
c
c                        -44            Problem not solvable with this
c                                       LIMEX version,  most  probable
c                                       reason:  initial  values   not
c                                       consistent  or  index  of  the 
c                                       DAE > 1.
c
c                        -45            A value of the solution vector
c                                       y within  the CIV  computation
c                                       is negative, but  should be by
c                                       => 0 by the settings of IPos.
c
c                        -46            More  stepsize reductions than
c                                       allowed  due  to failing  con-
c                                       vergence in the extrapolation. 
c                                       This  is  controlled  by   the 
c                                       parameter Max_Step_Red_Ex.
c
C                        -47            Singular matrix pencil B - hA,
c                                       problem   not  solvable   with
c                                       LIMEX.
c
c                        -48            More  integration  steps  than
c                                       allowed. This is controlled by
c                                       the parameter Max_Int_Steps. 
c
c                        -49            More internal Newton steps for
c                                       the computation  of consistent
c                                       initial  values than  allowed. 
c                                       This  is  controlled  by   the 
c                                       parameter Max_Step_CIV.
c
c                        -50            Stepsize  too  small, in  most
c                                       cases  due  to  many  stepsize
c                                       reductions.
c
c                IFail(2) Internal error code of the subroutine dgetrf
c                         or  dgbtrf (IFail(1) = - 39). See the LAPACK
c                         User's Guide for  an explanation. Zero if no 
c                         error occurred.
c
c                IFail(3) This  error code  signals whether  the error 
c                         occurred  during the  computation of CIVs or 
c                         not. IFail(3) = 0: error during integration,
c                         IFail(3) = 1: error  during CIV computation. 
c                         For  both IFail(3) = 0 or 1 the error  codes 
c                         IFail(1)  and  IFail(2) specify the  type of
c                         the error.
c
c     Remark: In  addition to  these errors  LIMEX may  output on  the 
c             current monitor unit the following warning messages:
c
c             1. Warning: component xxx does  not have  an  asymptotic 
c                         expansion in the first integration step
c
c                This messsage  occurrs if for  some components of the
c                solution  no asymptotic expansion in  powers of h can 
c                be found. The reason  for this can be either an index 
c                of the DAE which is greater than 1 or (in most cases) 
c                the lack of consistent initial  values. Thus, if this
c                message occurs  very frequently, the  computation  of 
c                consistent initial  values by setting Iopt(6) = 1 may 
c                prevent the  warning. However, in  many  cases  where 
c                this  warning is  printed, the  integration  proceeds 
c                without problems.
c                
c             2. Warning: switch  to direct  solution of  higher order 
c                         approximations
c
c                This warning  occurrs, if during  the integration the 
c                higher order approximations can not be computed using
c                a Richardson iteration, if the direct solution of the
c                linear  systems is  selected. If an  iterative method
c                is  selected, this message  indicates that the higher
c                order  solutions can  not be obtained with  the first 
c                order matrix ILU factorization as  preconditioner. In
c                both  cases, an appropriate  direct approach  will be 
c                tried, but this can  be more time  consuming. If  the 
c                problem is well posed, then this warning  indicates a
c                strongly varying matrix B(t,y). If the linear systems
c                are solved with direct methods, the maximum number of
c                iterations in the  Richardson  process  is defined by 
c                the variable Max_Iter_Rich. The value of the limit is
c                defined in a parameter  statement and may  be changed 
c                if necessary.
c 
c-----------------------------------------------------------------------
c
c     Subroutines and functions called
c
c-----------------------------------------------------------------------
c
c     Subroutines to be supplied by the user as externals:
c
c     Fcn ( n, nz, t, y, f, b, ir, ic, FcnInfo )
c     ------------------------------------------
c
c     This subroutine  computes the right-hand side of the DAE and the
c     non-zero entries in the left-hand side matrix B(t,y).
c
c     Parameters:
c
c     n          Integer variable, size  of the DAE system, should not
c                be changed.
c
c     nz         Integer  variable, must  be  set  in the  subroutine.
c                Number of non-zero entries in the matrix B(t,y).
c
c     t          Real  variable, should not be  changed. Current value
c                of the independent variable.
c
c     y          Real array of size n, should not  be changed. Current
c                values of  the solution  vector, for  which f(t,y) is
c                required.
c
c     f          Real array of  size n, must be set in the subroutine.
c                Values of f(t,y) at the current values of t and y.
c
c     b          Real array of size nz, must be set in the subroutine.
c                Values of the non-zero entries in  the matrix B(t,y).
c
c     ir         Integer  array of  size nz, must  be set  in the sub-
c                routine. Row  indices of  the non-zero entries in the 
c                matrix B(t,y).
c
c     ic         Integer  array of  size nz, must  be set  in the sub-
c                routine. Column  indices  of the  non-zero entries in 
c                the matrix B(t,y).
c
c     FcnInfo    Integer variable. FcnInfo indicates  the type  of the
c                call  of Fcn. If FcnInfo = 0, the  call is  the first
c                call of  Fcn  in  the  current  integration  step. If  
c                FcnInfo = 1 this  call is  within  the  extrapolation
c                process. FcnInfo = 2 indicates  a call of  Fcn within
c                the  numerical  evaluation of  the Jacobian. In  some
c                application  it may be  reasonable for efficiency, to 
c                compute f and/or B during the numerical evaluation of
c                the Jacobian only with  limited accuracy. Furthermore 
c                FcnInfo on return is  used as error  indicator. If an
c                error in the subroutine  Fcn occurs, set FcnInfo < 0.
c                Then LIMEX reduces the stepsize by the factor RedFcn. 
c                This  variable   is  currently  set  in  a  parameter 
c                statement  to 0.5, but  may  be  changed. If  such an
c                error  occurs during  the numerical evaluation of the
c                Jacobian, the program exits with an error message. If 
c                no error occurs, set FcnInfo = 0.
c
c     Remarks:   The order  of the  entries in  matrix B is arbitrary.
c                If the indices  ir and ic are constant, they  must be
c                defined only once in the first call of Fcn.
c                
c
c     Jacobian ( n, t, y, ys, Jac, LDJac, ml, mu, Full_or_Band, 
c     ---------------------------------------------------------
c                JacInfo )
c                ---------
c
c     This  subroutine computes  the Jacobian  matrix, i.e. analytical
c     derivatives of the residual f(t,y) - B(t,y) * y'. 
c
c     Parameters:
c
c     n          Integer variable, size  of the DAE system, should not
c                be changed.
c
c     t          Real  variable, should not be  changed. Current value
c                of the independent variable.
c
c     y          Real array of size n, should not  be changed. Current
c                values of  the solution  vector, for  which f(t,y) is
c                required.
c
c     ys         Real array of size n, should not  be changed. Current
c                values of the derivatives y'.
c
c     Jac        Real matrix of  dimension at least (LDJac,n), must be  
c                set in the subroutine. The derivatives at the current 
c                values of y, ys and t. 
c
c     LDJac      Integer variable, the leading dimension of the matrix 
c                Jac as declared in LIMEX, should not be changed.
c
c     ml         Integer  variable, the  number of  lower diagonals in
c                the Jacobian, should not be changed. 
c
c     mu         Integer  variable, the  number of  upper diagonals in
c                the Jacobian, should not be changed. 
c
c     Full_or_Band
c                Integer variable, should not be changed.
c                = 0 : Jacobian matrix is a full matrix.
c                = 1 : Jacobian matrix is a band matrix.
c
c     JacInfo    Integer variable, error indicator. If an error in the
c                subroutine  Jacobian  occurs,  set  JacInfo < 0. Then 
c                LIMEX aborts with an error  message. If JacInfo => 0,
c                no error is assumed.
c
c     Remarks:   The general form of  this subroutine should look like
c                the  following  code  segment. This  works  for  both
c                matrix representations. The storage scheme for banded
c                matrices is the usual LINPACK/LAPACK format.
c
c                do j = 1, n
c
c                   iUpper = max ( 1, j - mu )
c                   iLower = min ( n, j + ml )
c
c                   if ( Full_or_Band .eq. 0 ) then
c                      k = 0
c                   else
c                      k = mu + 1 - j
c                   end if
c
c                   do i = iUpper, iLower
c
c                      Now define matrix element with indices i, j:
c
c                      Jac(i+k,j) = ...
c
c                   end do
c
c                end do
c
c                If the  numerical  approximation  of the  Jacobian is
c                used, one may supply a dummy routine. If the Jacobian
c                is a  sparse  matrix, it  can be  more  efficient, to 
c                compute  the numerical approximation  of the Jacobian
c                also  by this  routine, even  if  for  such  problems 
c                LIMEX4.2B should be preferable.
c
c     Other subroutines
c
c     From LIMEX_Dense:
c     -----------------
c
c     Comp_Herm  Computes  right  derivatives  approximations and  the
c                coefficients of the Hermite interpolation.
c
c     Eval_Herm  Evaluates the Hermite interpolation polynomial.
c
c     JacPlot    Generates a PostScript plot of the Jacobian matrix.
c
c     From the BLAS library:
c     ----------------------
c
c     daxpy      Adds a scalar multiple of a  vector to another vector
c                (used only in LIMEX4_2A1_Dense.f).
c
c     dcopy      Copies real vectors
c
c     dnrm2      Computes the Euclidean norm of a real vector
c
c     dscal      Scales a real vector by a constant
c
c     From the LAPACK Library:
c     ------------------------
c
c     dgbtrf     Factorizes a banded real matrix
c
c     dgbtrs     Solves a banded real linear system
c
c     dgetrf     Factorizes a general real matrix
c
c     dgetrs     Solves a general real linear system
c
c-----------------------------------------------------------------------
c
c     Declaration of the arguments
c
c-----------------------------------------------------------------------
c
      integer            n, Iopt ( 30 ), IPos ( * ), IFail ( 3 )
c
      double precision   t_Begin, t_End, y ( * ), ys ( * ),  rTol ( * ),
     2                   aTol ( * ), h, Ropt ( 5 )
c
      external           Fcn, Jacobian
c
c-----------------------------------------------------------------------
c
c     Define size of problem via the include file
c
c        'xlimex_def.h'. 
c
c     See the installation notes for a detailed description.
c
c-----------------------------------------------------------------------
c
      integer            Max_Nr_of_Equations, Max_Non_Zeros_Jacobian, 
     2                   Max_Non_Zeros_B,
     3                   Max_Lower_Diagonals, Max_Upper_Diagonals,
     4                   Max_It_Vectors
c
c-----------------------------------------------------------------------
c
      include 'xlimex_def.h'
c
c-----------------------------------------------------------------------
c
c     Derived size definitions
c
c-----------------------------------------------------------------------
c
      integer            Max_Size
c
      parameter        ( Max_Size = Max_Nr_of_Equations )
c
c-----------------------------------------------------------------------
c
c     Control parameters. May be changed by skillful users.
c
c-----------------------------------------------------------------------
c
c     Max_Int_Steps      The maximum  number of integration  steps per 
c                        interval permitted
c
c     Max_Step_Red_Ex    The  maximum  number of  stepsize  reductions 
c                        permitted due to  failing  convergence in the 
c                        extrapolation tableau or negative  components
c                        of the solution vector y, which by the values
c                        of IPos should be positive
c
c     Max_Step_Red_Dc    The  maximum  number  of stepsize  reductions 
c                        permitted  due to zero pivot  elements in the 
c                        LU decomposition
c
c     Max_Row_Tab        The maximum row  number in the  extrapolation 
c                        tableau
c
c     Max_Step_CIV       The maximum  number of Newton steps permitted
c                        for  the  computation of  consistent  initial 
c                        values
c
      integer            Max_Int_Steps, Max_Step_Red_Ex, 
     2                   Max_Step_Red_Dc, Max_Row_Tab, Max_Step_CIV
c
      parameter        ( Max_Int_Steps      = 1000000,
     2                   Max_Step_Red_Ex    =  5000,
     3                   Max_Step_Red_Dc    =     5,
     4                   Max_Row_Tab        =     7,
     5                   Max_Step_CIV       =    100)
c
c-----------------------------------------------------------------------
c     
c     Local scalars
c
c-----------------------------------------------------------------------
c
      logical            CIV_Convergence, Convergence, Direct
c
      integer            DensOut, FcnInfo, Full_or_Band, HoldJac, i, 
     2                   icc, iCIV, iCont, iLower, iRed, irh, irr, 
     3                   iSing, Iter_Rich, iUpper, j, j_Cons, JacAct, 
     4                   JacInfo, JacP, jCIV, jk, jo, joh, jmax, jmh, 
     5                   jmm, k, k_Cons, k_Dense, k1, kc1, kc2, kl, km, 
     6                   ko, koh, mb1, mb2, ml, MonOut, mu, n_Dense, 
     7                   nDec, NewJac, nFcn, nFcnE, nFcnJ, nJac, njj, 
     8                   NoDense, NonPosComp, nSol, nStep, nz, SolOut, 
     9                   Status
c
      double precision   abs_del, CErr, Cons_Err, d1, d2, delm,
     2                   delta_Dense, delta_y, Delta_y0,delta_y_rec, 
     3                   del_quot, eph, eps_Mach, eps_Mach10, EtaDif,
     4                   fc, fcm, fco, fn, g, g_inv, h_proposed, h_True,
     5                   hMax, hMax_Int, hSave, Mean_rTol, omj, omjo, 
     6                   red, rTol2, t, t_Dense, tEnd_t, th, tmp, 
     7                   tn_Dense, TolQuot, WarnCont, y_orig
c
      double precision   dnrm2
c
c-----------------------------------------------------------------------
c
c     Arrays for the LIMEX core algorithm
c
c-----------------------------------------------------------------------
c
      double precision   del   ( Max_Size ), 
     2                   dz    ( Max_Size ), 
     3                   Scal  ( Max_Size ),
     4                   Temp1 ( Max_Size ), 
     5                   Temp2 ( Max_Size ),
     6                   Temp3 ( Max_Size ),
     7                   Temp4 ( Max_Size ),
     8                   TolQu ( Max_Size ), 
     9                   yk    ( Max_Size ) 
c
c-----------------------------------------------------------------------
c
c     Arrays for the extrapolation and order control
c
c-----------------------------------------------------------------------
c
      integer            nj    ( Max_Row_Tab ), 
     2                   incr  ( Max_Row_Tab )
c
      double precision   aj    ( Max_Row_Tab )
c
      double precision   al    ( Max_Row_Tab, Max_Row_Tab ),
     2                   d     ( Max_Row_Tab, Max_Row_Tab )
c
      double precision   dt    ( Max_Size, Max_Row_Tab ),
     2                   dtp   ( Max_Size, Max_Row_Tab )
c
c-----------------------------------------------------------------------
c
c     Arrays for the left-hand side matrix B 
c
c-----------------------------------------------------------------------
c
      integer            ic    ( Max_Non_Zeros_B ),
     2                   ir    ( Max_Non_Zeros_B )
c
      double precision   b0    ( Max_Non_Zeros_B ),
     2                   bk    ( Max_Non_Zeros_B )
c
c-----------------------------------------------------------------------
c
c     Arrays for dense output
c
c-----------------------------------------------------------------------
c
      integer            ipt   ( Max_Row_Tab )
c
      double precision   Dense ( Max_Size, 
     2                   2 + ( Max_Row_Tab * ( Max_Row_Tab + 1 ) ) / 2 )
c
c-----------------------------------------------------------------------
c
c     Arrays for the Jacobian, the minimum of  memory requirements for
c     full or  banded matrix storage is determined  via some parameter
c     calculations.
c
c-----------------------------------------------------------------------
c
      integer            Pivot ( Max_Size )
c
      integer            LDJac, LDBJac, 
     2                   Max_Band_Size_1, Max_Band_Size_2,
     3                   mbs1, mbs2, mbs3, mbs4
c
      parameter        ( Max_Band_Size_1 =
     2                   1 + Max_Lower_Diagonals + Max_Upper_Diagonals )
c
      parameter        ( Max_Band_Size_2 =
     2                           Max_Band_Size_1 + Max_Lower_Diagonals )
c
      parameter        ( mbs1 = Max_Band_Size_2 / Max_Size, 
     2                   mbs2 = Max_Size / Max_Band_Size_2  )
c
      parameter        ( mbs3 = ( mbs1 + 2 ) / ( mbs1 + 1 ) - 1,
     2                   mbs4 = ( mbs2 + 2 ) / ( mbs2 + 1 ) - 1  )
c
      parameter        ( LDJac  = 
     2                        mbs3 * Max_Band_Size_1 + mbs4 * Max_Size,
     3                   LDBJac = 
     4                        mbs3 * Max_Band_Size_2 + mbs4 * Max_Size )
c
      double precision   Jac   ( LDJac,  Max_Size ),
     2                   B_Jac ( LDBJac, Max_Size )
c
c-----------------------------------------------------------------------
c
c     Internal parameters
c
c-----------------------------------------------------------------------
c
      integer            iCons1, iCons2, Max_Iter_Rich
c
      double precision   CStep, Epmin, ex1, ex2, fmin, h_Min, one, one1,
     2                   quart, RedFcn, RedPos, rmin, ro, safe, Small,
     3                   ten, ThMin, thresh, zero
c
      parameter ( CStep  = 1.0d-7,  Epmin  = 1.0d-10, ex1    = 0.6d0  , 
     2            ex2    = 1.5d0 ,  fmin   = 1.0d-2 , h_Min  = 1.0d-4 ,
     3            one    = 1.0d0 ,  one1   = 1.01d0 , quart  = 0.25d0 , 
     4            RedFcn = 0.5d0 ,  RedPos = 0.5d0  , rmin   = 0.9d0  ,
     5            ro     = 0.25d0,  safe   = 0.5d0  , Small  = 1.0d-30, 
     6            ten    = 10.0d0,  ThMin  = 2.0d-2 , thresh = 0.1d0  , 
     7            zero   = 0.0d0 ) 
c
      parameter ( iCons1 = 4, iCons2 = 4, Max_Iter_Rich = 10 )
c
c
c-----------------------------------------------------------------------
c
c     The following save statement is needed on some machines only for
c     the one step mode or for  continuation calls. Then all variables
c     should retain their values after a return from LIMEX.
c
c-----------------------------------------------------------------------
c
      save
c
c-----------------------------------------------------------------------
c
c     Define initial program status.
c
c-----------------------------------------------------------------------
c
      data Status / 0 /
c
c-----------------------------------------------------------------------
c
c     Check program status.
c
c-----------------------------------------------------------------------
c
      WarnCont = 0
c
      if ( Iopt(16) .eq. 0 ) then
c
         Status = 0
c
      else if ( Iopt(16) .eq. 1 ) then
c
         if ( Status .eq. 0 ) WarnCont = 1 
c
      else
c
         Status = 0
c
      end if
c
      if ( Status .eq. 1 ) go to 120 
c
      if ( Status .ge. 2 ) go to 310
c
      Status = 1
c
c-----------------------------------------------------------------------
c
c     Input data check and default settings.
c
c-----------------------------------------------------------------------
c
      IFail(1) = 0
      IFail(2) = 0
      IFail(3) = 0
c
      if ( Iopt(1) .lt. 0 .or. Iopt(1) .gt. 2 ) then
c
         if ( Iopt(2) .lt. 0 ) then
            write ( *, 9000 )
            write ( *, 9010 ) 1, Iopt(1)
            write ( *, 9020 ) 2
         else
            write ( Iopt(2), 9000 )
            write ( Iopt(2), 9010 ) 1, Iopt(1)
            write ( Iopt(2), 9020 ) 2
         end if
c
         IFail(1) = - 4
         go to 510
c
      end if
c
      if ( Iopt(1) .gt. 0 ) then
c
         if ( Iopt(2) .lt. 0 ) then
            write ( *, 9000 )
            write ( *, 9010 ) 2, Iopt(2)
            write ( *, 9030 )
            IFail(1) = - 5
            go to 510
         else if ( Iopt(2) .eq. 0 ) then
            MonOut = 6
         else
            MonOut = Iopt(2)
         end if
c
      end if
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9000 )
c
      if ( Max_Nr_of_Equations .lt. n .or. n .lt. 1 ) then
         if ( Iopt(1) .ne. 0 ) 
     2      write ( MonOut, 9040 ) n, Max_Nr_of_Equations
         IFail(1) = - 1
         go to 510
      end if
c
      if ( Iopt(3) .lt. 0 .or. Iopt(3) .gt. 2 ) then
         if ( Iopt(1) .ne. 0 ) then
            write ( MonOut, 9010 ) 3, Iopt(3)
            write ( MonOut, 9020 ) 2 
         end if
         IFail(1) = - 6
         go to 510
      end if
c
      if ( Iopt(3) .gt. 0 ) then
c
         if ( Iopt(4) .lt. 0 ) then
            if ( Iopt(1) .ne. 0 ) then
               write ( MonOut, 9010 ) 4, Iopt(4)
               write ( MonOut, 9030 ) 
            end if
            IFail(1) = - 7
            go to 510
         else if ( Iopt(4) .eq. 0 ) then
            SolOut = 6
         else
            SolOut = Iopt(4)
         end if
c
      end if
c
      do i = 5, 7
c
         if ( Iopt(i) .lt. 0 .or. Iopt(i) .gt. 1 ) then
            if ( Iopt(1) .ne. 0 ) then
               write ( MonOut, 9010 ) i, Iopt(i)
               write ( MonOut, 9020 ) 1
            end if
            IFail(1) = - ( i + 3 )
            go to 510
         end if
c
      end do
c
      if ( Iopt(8) .lt. 0 .or. Iopt(8) .gt. Max_Lower_Diagonals ) then
         if ( Iopt(1) .ne. 0 ) then
            write ( MonOut, 9010 ) 8, Iopt(8)
            write ( MonOut, 9020 ) Max_Lower_Diagonals 
         end if
         IFail(1) = - 11
         go to 510
      end if
c
      if ( Iopt(9) .lt. 0 .or. Iopt(9) .gt. Max_Upper_Diagonals ) then
         if ( Iopt(1) .ne. 0 ) then
            write ( MonOut, 9010 ) 9, Iopt(9)
            write ( MonOut, 9020 ) Max_Upper_Diagonals 
         end if
         IFail(1) = - 12
         go to 510
      end if
c
      do i = 10, 11
c
         if ( Iopt(i) .lt. 0 .or. Iopt(i) .gt. 1 ) then
            if ( Iopt(1) .ne. 0 ) then
               write ( MonOut, 9010 ) i, Iopt(i)
               write ( MonOut, 9020 ) 1
            end if
            IFail(1) = - ( i + 3 )
            go to 510
         end if
c
      end do
c
      if ( Iopt(12) .lt. 0 .or. Iopt(12) .gt. 2 ) then
         if ( Iopt(1) .ne. 0 ) then
            write ( MonOut, 9010 ) i, Iopt(i)
            write ( MonOut, 9020 ) 2
         end if
         IFail(1) = - 15
         go to 510
      end if
c
      if ( Iopt(13) .lt. 0 .or. Iopt(13) .gt. 3 ) then
         if ( Iopt(1) .ne. 0 ) then
            write ( MonOut, 9010 ) 13, Iopt(13)
            write ( MonOut, 9020 ) 3 
         end if
         IFail(1) = - 16
         go to 510
      end if
c
      if ( Iopt(13) .gt. 0 ) then
c
         if ( Iopt(14) .lt. 0 .and. Iopt(13) .ne. 3 ) then
            if ( Iopt(1) .ne. 0 ) then
               write ( MonOut, 9010 ) 14, Iopt(14)
               write ( MonOut, 9030 )
            end if
            IFail(1) = - 17
            go to 510
         end if
c
         if ( Iopt(15) .lt. 0 ) then
            if ( Iopt(1) .ne. 0 ) then
               write ( MonOut, 9010 ) 15, Iopt(15)
               write ( MonOut, 9030 )
            end if
            IFail(1) = - 18
            go to 510 
         else 
            DensOut = Iopt(15)
         end if
c
      end if
c
      if (        Iopt(13) .eq. 0
     2     .or. ( Iopt(13) .eq. 1 .and. Iopt(14) .le. 1 )
     3     .or. ( Iopt(13) .eq. 2 .and. Iopt(14) .eq. 0 ) ) then
         DensOut = 0
         NoDense = 0
      else
         NoDense = 1 
      end if
c
      do i = 16, 17
c
         if ( Iopt(i) .lt. 0 .or. Iopt(i) .gt. 1 ) then
            if ( Iopt(1) .ne. 0 ) then
               write ( MonOut, 9010 ) i, Iopt(i)
               write ( MonOut, 9020 ) 1 
            end if
            IFail(1) = - ( i + 3 )
            go to 510
         end if 
c
      end do
c
      if ( Iopt(18) .lt. - 1 ) then
         if (Iopt(1) .ne. 0 ) then
            write ( MonOut, 9010 ) 18, Iopt(18)
            write ( MonOut, 9050 )
         end if
         IFail(1) = - 21
         go to 510
      end if
c
      if ( Ropt(1) .lt. zero ) then
         if ( Iopt(1) .ne. 0 ) then
            write ( MonOut, 9060 ) 1, Ropt(1)
            write ( MonOut, 9030 )
         end if
         IFail(1) = - 27
         go to 510
      else if ( Ropt(1) .eq. zero ) then
         if ( Iopt(17) .eq. 1 ) then
            if ( Ropt(3) .ge. t_End ) then
               hMax = Ropt(3) - t_Begin
            else 
               hMax = zero
            end if
         else
            hMax = t_End - t_Begin
         end if
      else
         if ( Iopt(17) .eq. 1 .and. Ropt(3) .ge. t_End ) then
            hMax = min ( Ropt(1), Ropt(3) - t_Begin )
         else
            hMax = Ropt(1)
         end if
      end if
c
      if ( Ropt(2) .le. zero .and. Iopt(13) .eq. 3 ) then
         if ( Iopt(1) .ne. 0 ) then
            write ( MonOut, 9060 ) 2, Ropt(2)
            write ( MonOut, 9070 )
         end if
         IFail(1) = - 28
         go to 510
      end if
c
      if ( Iopt(11) .eq. 0 ) then
c
         if ( rTol(1) .le. zero ) then
            if ( Iopt(1) .ne. 0 ) write ( MonOut, 9080 ) 1, rTol(1)
            IFail(1) = - 2
            go to 510
         end if
c
         if ( aTol(1) .lt. zero ) then
            if ( Iopt(1) .ne. 0 ) then
               write ( MonOut, 9090 ) 1, aTol(1)
               write ( MonOut, 9030 )
            end if
            IFail(1) = - 3
            go to 510
         end if
c
      else
c
         do i = 1, n
c
            if ( rTol(i) .le. zero ) then
               if ( Iopt(1) .ne. 0 ) write ( MonOut, 9080 ) i, rTol(i)
               IFail(1) = - 2
               go to 510
            end if
c
         end do
c
         do i = 1, n
c
            if ( aTol(i) .lt. zero ) then
               if ( Iopt(1) .ne. 0 ) then
                  write ( MonOut, 9090 ) i, aTol(i)
                  write ( MonOut, 9030 )
               end if
               IFail(1) = - 3
               go to 510
            end if
c
         end do
c
      end if
c
      do i = 1, n
c
         if ( IPos(i) .eq. 1 .and. y(i) .lt. zero ) then
            if ( Iopt(1) .ne. 0 ) write ( MonOut, 9100 ) i, y(i)
            IFail(1) = - 32
            go to 510
         end if
c
      end do
c
c-----------------------------------------------------------------------
c
c     Print a protocol of the input parameters.
c
c-----------------------------------------------------------------------
c
      if ( Iopt(1) .gt. 0 ) then
c
         write ( MonOut, 9110 )
c
         if ( Iopt(1) .eq. 1 ) then
            write ( MonOut, 9120 ) MonOut
         else
            write ( MonOut, 9130 ) MonOut
         end if
c
         if ( Iopt(3) .eq. 0 ) then
            write ( MonOut, 9140 )
         else if ( Iopt(3) .eq. 1 ) then
            write ( MonOut, 9150 ) SolOut
         else
            write ( MonOut, 9160 ) SolOut
         end if
c
         write ( MonOut, 9170 ) n
c
         write ( MonOut, 9180 ) t_Begin, t_End
c
         if ( Iopt(5) .eq. 0 ) then
            write ( MonOut, 9190 )
         else
            write ( MonOut, 9200 )
         end if
c
         if ( Iopt(6) .eq. 0 ) then
            write ( MonOut, 9210 )
         else
            write ( MonOut, 9220 )
         end if
c
         if ( Iopt(7) .eq. 0 ) then
            write ( MonOut, 9230 )
         else
            write ( MonOut, 9240 )
         end if
c
         if ( Iopt(7) .eq. 0 )
     2      write ( MonOut, 9250 ) Iopt(8), Iopt(9)
c
         if ( Iopt(10) .eq. 0 ) then
            write ( MonOut, 9260 )
         else 
            write ( MonOut, 9270 ) ThMin  
         end if
c
         if ( Iopt(11) .eq. 0 ) then
            write ( MonOut, 9280 ) rTol(1), aTol(1)
         else
            write ( MonOut, 9290 )
         end if
c
         if ( Iopt(12) .eq. 0 ) then
            write ( MonOut, 9300 ) 
         else if ( Iopt(12) .eq. 1 ) then
            write ( MonOut, 9310 )
         else
            if ( NoDense .eq. 0 ) then
               write ( MonOut, 9300 ) 
            else if ( Iopt(13) .eq. 1 ) then
               write ( MonOut, 9320 ) Iopt(14)
            else if ( Iopt(13) .eq. 2 ) then
               write ( MonOut, 9330 ) Iopt(14)
            else
               write ( MonOut, 9340 ) Ropt(2)
            end if
         end if
c
         if ( Iopt(13) .eq. 0 .or. DensOut .eq. 0 ) then
            write ( MonOut, 9350 )
         else if ( Iopt(13) .eq. 1 ) then
            write ( MonOut, 9360 ) Iopt(14), DensOut
         else if ( Iopt(13) .eq. 2 ) then
            write ( MonOut, 9370 ) Iopt(14), DensOut
         else
            write ( MonOut, 9380 ) Ropt(2), DensOut
         end if
c
         if ( Iopt(17) .eq. 0 ) then
            write ( MonOut, 9390 ) 
         else 
            if ( Ropt(3) .ge. t_End ) then
               write ( MonOut, 9400 ) Ropt(3)
            else
               write ( MonOut, 9410 )
            end if
         end if
c
         if ( Iopt(18) .eq. 0 ) then
            write ( MonOut, 9420 )
            JacP = 0
         else 
            if ( Iopt(18) .gt. 0 ) then
               write ( MonOut, 9430 ) Iopt(18)
            else
               write ( MonOut, 9430 ) 0
            end if  
            JacP = Iopt(18)
         end if
c
         if ( hMax .ne. zero ) then
            write ( MonOut, 9440 ) hMax
         else
            write ( MonOut, 9450 )
         end if
c
         if ( WarnCont .eq. 1 ) write ( MonOut, 9460 )
c
      end if
c
c-----------------------------------------------------------------------
c
c     Some initial preparations.
c
c-----------------------------------------------------------------------
c
      nFcn  = 0
      nFcnE = 0
      nFcnJ = 0
      nDec  = 0
      nJac  = 0
      nSol  = 0
      nStep = 0
c
      HoldJac = 0
c
      t = t_Begin
c
      fn = dsqrt ( one / dble ( n ) )
c
c-----------------------------------------------------------------------
c
c     For h equal zero, the initial stepsize is set to h_Min.
c
c-----------------------------------------------------------------------
c
      if ( h .eq. zero ) h = h_Min
c
c-----------------------------------------------------------------------
c
c     Determine the maximum permitted order.
c
c-----------------------------------------------------------------------
c
      if ( Iopt(5) .eq. 1 ) then
         jmax = Max_Row_Tab
      else
         jmax = 5
      end if
c
c-----------------------------------------------------------------------
c
c     Determine the smallest positive real number (machine dependent).
c
c-----------------------------------------------------------------------
c
      call Get_eps_Mach(eps_Mach)
      eps_Mach10 = ten * eps_Mach
c
c-----------------------------------------------------------------------
c
c     Define delta_x for the numerical differentiation. 
c
c-----------------------------------------------------------------------
c
      EtaDif = sqrt ( eps_Mach10 )
c
c-----------------------------------------------------------------------
c
c     Define the stepsize sequence for LIMEX.
c
c-----------------------------------------------------------------------
c
      do i = 1, jmax
         nj(i) = i
      end do
c
c-----------------------------------------------------------------------
c
c     Preparations for consistent initial value (CIV) computation.
c
c-----------------------------------------------------------------------
c
      iCIV = Iopt(6)
      jCIV = 0
c
      if ( iCIV .eq. 1 ) then
c
         hSave = h
         h     = CStep 
c
         j_Cons = iCons1
         k_Cons = iCons2
c
      end if
c
c-----------------------------------------------------------------------
c
c     Initial definition of scaling parameters.
c
c-----------------------------------------------------------------------
c
      if ( Iopt(11) .eq. 0 ) then
c
         TolQuot = aTol(1) / rTol(1)
         rTol2   = 0.5d0 * rTol(1)
c
         do i = 1, n
            Scal(i) = abs ( y(i) ) + TolQuot
         end do
c
      else
c
         do i = 1, n
c
            TolQu(i) = aTol(i) / rTol(i)
            Scal(i)  = abs ( y(i) ) + TolQu(i)
c
         end do
c
      end if
c
c-----------------------------------------------------------------------
c
c     Determine  optimal handling of the matrices, i.e. full or banded
c     matrix representation.
c
c-----------------------------------------------------------------------
c
      ml = Iopt(8)
      mu = Iopt(9)
c
      if ( 2 * ml + mu + 1 .ge. n ) then
c
         Full_or_Band = 0
c
      else
c
         Full_or_Band = 1
c
         mb1 = mu + 1
         mb2 = ml + mu + 1
c
      end if
c
c-----------------------------------------------------------------------
c
c     Set the parameters for extrapolation and order control.
c
c-----------------------------------------------------------------------
c
      if ( Iopt(11) .eq. 0 ) then
         Mean_rTol = rTol(1)
      else
         Mean_rTol = fn * dnrm2 ( n, rTol, 1 )
      end if
c
      eph = ro * Mean_rTol
c
      aj(1) = dble ( nj(1) ) + sqrt ( dble ( n ) )
c
      do j = 2, jmax
c
         incr(j-1) = 0
c
         do i = 1, j - 1
            d(j,i) = dble ( nj(j) ) / dble ( nj(i) )
         end do
c
         aj(j) = aj(j-1) + dble ( nj(j) ) - one
c
         do i = 2, j-1
            al(j-1,i-1) = eph **
     2     ( ( aj(i) - aj(j) ) / ( ( aj(j) - aj(1) ) * dble ( i ) ) )
         end do
c
      end do
c
      do j = 2, jmax - 1
c
         koh = j
c
         if ( aj(j+1) * one1 .gt. aj(j) * al(j,j-1) ) go to 110
c
      end do
c
  110 continue
c
      joh       = koh + 1
      km        = koh
      jmh       = joh
      incr(jmh) = - 1
      omjo      = zero
      k         = 0
c
c-----------------------------------------------------------------------
c
c     Preparations for dense output.
c
c-----------------------------------------------------------------------
c
      if ( NoDense .eq. 1 .or. Iopt(17) .eq. 1 ) then
c
         ipt(1) = 3
c
         do j = 2, jmax
            ipt(j) = ipt(j-1) + nj(j)
         end do
c
      end if
c
c-----------------------------------------------------------------------
c
c     Print some control variables and/or initial values.
c
c-----------------------------------------------------------------------
c
      if ( Iopt(3) .ge. 1 ) write ( SolOut, 9470 ) t, ( y(i), i = 1, n )
c
      if (          DensOut .ne. 0 .and. iCIV .eq. 0 
     2      .and. ( Iopt(3) .lt. 2 .or. DensOut .ne. SolOut ) )
     3   write ( DensOut, 9480 ) t_Begin, ( y(i), i = 1, n )
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9490 )
c
      if (      abs ( t_End - t ) .lt. abs ( t ) * eps_Mach10 * ten
     2     .or. abs ( t_End - t ) .lt. abs ( h ) * eps_Mach10 ) 
     3   go to 360
c
c-----------------------------------------------------------------------
c
c     This is the begin of the basic integration step.
c
c-----------------------------------------------------------------------
c
  120 continue
c 
      if ( Iopt(13) .eq. 1 .and. Iopt(14) .gt. 1 ) then
c
         delta_Dense = ( t_End - t_Begin ) / dble ( Iopt(14) - 1 )
         tn_Dense    = t_Begin + delta_Dense
c
      end if
c
      if ( Iopt(17) .eq. 1 .and. Ropt(3) .ge. t_End ) then
         iCont = 1
      else
         iCont = 0
      end if
c
  130 continue
c
      t        = t_Begin
      tEnd_t   = t_End - t
      hMax_Int = abs ( hMax )
c
c-----------------------------------------------------------------------
c
c     Define order for next integration step resp. for the computation    
c     of consistent initial values.
c
c-----------------------------------------------------------------------
c
      if ( iCIV .eq. 0 ) then
c
         kc1 = k
         kc2 = koh
c
      else
c
         kc1 = j_Cons
         kc2 = k_Cons
c
         j_Cons = iCons1
         k_Cons = iCons2
c
      end if
c
c-----------------------------------------------------------------------
c
c     Monitor the current integration step.
c
c-----------------------------------------------------------------------
c
      if ( Iopt(16) .eq. 0 .or. jCIV .eq. 1 ) then
c
         Iopt(16) = 1
         jCIV     = 0
c
         if ( Iopt(1) .ge. 1 ) then
c
            if ( HoldJac .eq. 0 ) then
               write ( MonOut, 9500 ) nStep, '*', nFcn, t, h, kc1, kc2
            else
               write ( MonOut, 9500 ) nStep, ' ', nFcn, t, h, kc1, kc2
            end if
C
         end if
c
      end if
c
      if ( iCIV .eq. 0 .and. Iopt(3) .eq. 2 .and. nStep .gt. 0 ) 
     2   write ( SolOut, 9510 ) t, ( y(i), i = 1, n )
c
c-----------------------------------------------------------------------
c
c     Resize  h, if t_End will  be reached  in the current integration
c     step or can reached by a small enlargement of the stepsize. 
c
c-----------------------------------------------------------------------
c
      if ( Iopt(17) .ne. 1 ) then
c
         if ( abs ( tEnd_t ) .lt. one1 * abs ( h ) ) h = tEnd_t 
c
         hMax_Int = min ( abs ( hMax ), abs ( tEnd_t ) )
c
      else if ( iCont .eq. 1 ) then
c
         if ( abs ( Ropt(3) - t ) .lt. one1 * abs ( h ) ) 
     2      h = Ropt(3) - t
c
         hMax_Int = min ( abs ( hMax ), abs ( Ropt(3) - t ) )
c
      end if
c
c-----------------------------------------------------------------------
c
c     Reset some counters.
c
c-----------------------------------------------------------------------
c
      iRed  = 0
      iSing = 0
c
      do i = 1, km
         incr(i) = incr(i) + 1
      end do
c
c-----------------------------------------------------------------------
c
c     Compute f(y(0)) - B(y(0)) * y'(0) and check nz.
c
c-----------------------------------------------------------------------
c
      FcnInfo = 0
c
      call Fcn ( n, nz, t, y, dz, b0, ir, ic, FcnInfo )
c
      nFcn  = nFcn + 1
      nFcnE = nFcnE + 1
c
      if ( FcnInfo .lt. 0 ) go to 370
c
      if ( nz .gt. Max_Non_Zeros_B ) go to 380
c
      call dcopy ( n, ys, 1, Temp4, 1 )
c
c-----------------------------------------------------------------------
c
c     Label  for return after a stepsize reduction with a  non-current
c     Jacobian.
c
c-----------------------------------------------------------------------
c
  140 continue
c
      call dcopy ( n, dz, 1, yk, 1 )
c
      do i = 1, nz
         yk(ir(i)) = yk(ir(i)) - b0(i) * ys(ic(i)) 
      end do
c
c-----------------------------------------------------------------------
c
c     Compute new Jacobian matrix.
c
c-----------------------------------------------------------------------
c
      if ( HoldJac .eq. 0 ) then
c 
         JacAct = 1
         nJac   = nJac + 1
c
c-----------------------------------------------------------------------
c
c     Numerical difference approximation of the Jacobian  with control
c     of the discretization and scaling.
c
c-----------------------------------------------------------------------
c
         if ( Iopt(7) .eq. 0 ) then
c
c-----------------------------------------------------------------------
c
c     The Jacobian matrix is a full matrix.
c
c-----------------------------------------------------------------------
c
            if ( Full_or_Band .eq. 0 ) then
c
               do j = 1, n
c
                  iUpper = max ( 1, j - mu )
                  iLower = min ( n, j + ml )
c
                  y_orig  = y(j)
                  delta_y = sign ( Scal(j) * EtaDif, - yk(j) )
                  y(j)    = y_orig + delta_y
c
                  if ( IPos(j) .eq. 1 .and. y(j) .lt. zero ) then
c
                     delta_y = - delta_y
                     y(j)    = y_orig + delta_y
c
                  end if
c
                  FcnInfo = 2
c
                  call Fcn ( n, nz, t, y, del, bk, ir, ic, FcnInfo )
c
                  nFcn  = nFcn + 1
                  nFcnJ = nFcnJ + 1
c
                  y(j) = y_orig
c
                  if ( FcnInfo .lt. 0 ) go to 390
c
                  do i = 1, nz
                     del(ir(i)) = del(ir(i)) - bk(i) * ys(ic(i))
                  end do
c
                  delta_y_rec = Scal(j) / delta_y
c
                  do i = iUpper, iLower
c
                     if ( yk(i) .ne. del(i) ) then
                        Jac(i,j) = 
     2                    ( yk(i) - del(i) ) * ( delta_y_rec / Scal(i) )
                     else  
                        Jac(i,j) = zero
                     end if
c
                  end do
c
               end do
c
            else
c
c-----------------------------------------------------------------------
c
c     The Jacobian matrix is a banded matrix.
c
c-----------------------------------------------------------------------
c
               do j = 1, mb2
c
                  do k = j, n, mb2
c
                     Temp1(k) = y(k)
                     Temp2(k) = sign ( Scal(k) * EtaDif, - yk(k) )
                     y(k)     = y(k) + Temp2(k)
c
                     if ( IPos(k) .eq. 1 .and. y(k) .lt. zero ) then
c
                        Temp2(k) = - Temp2(k)
                        y(k)     = Temp1(k) + Temp2(k)
c
                     end if
c
                  end do
c
                  FcnInfo = 2
c
                  call Fcn ( n, nz, t, y, del, bk, ir, ic, FcnInfo )
c
                  nFcn  = nFcn + 1
                  nFcnJ = nFcnJ + 1
c
                  do k = j, n, mb2
                     y(k) = Temp1(k)
                  end do 
c
                  if ( FcnInfo .lt. 0 ) go to 390
c
                  do i = 1, nz
                     del(ir(i)) = del(ir(i)) - bk(i) * ys(ic(i))
                  end do
c
                  do k = j, n, mb2
c
                     iUpper = max ( 1, k - mu )
                     iLower = min ( n, k + ml )
                     k1     = mb1 - k
c
                     delta_y_rec = Scal(k) / Temp2(k)
c
                     do i = iUpper, iLower
c
                        if ( yk(i) .ne. del(i) ) then
                           Jac(k1+i,k) =
     2                  ( yk(i) - del(i) ) * ( delta_y_rec / Scal(i) )
                        else
                           Jac(k1+i,k) = zero
                        end if
c
                     end do
c
                  end do
c
               end do
c
            end if
c
c-----------------------------------------------------------------------
c
c     Analytical computation of the Jacobian.    
c
c-----------------------------------------------------------------------
c
         else
c
            call Jacobian ( n, t, y, ys, Jac, LDJac, ml, mu, 
     2                      Full_or_Band, JacInfo )
c
            if ( JacInfo .lt. 0 ) go to 400
c
c-----------------------------------------------------------------------
c
c     Scaling of the analytical Jacobian.
c
c-----------------------------------------------------------------------
c
            do j = 1, n
c
               iUpper = max ( 1, j - mu )
               iLower = min ( n, j + ml )
c
               if ( Full_or_Band .eq. 0 ) then
                  k = 0
               else
                  k = mb1 - j
               end if
c
               tmp = - Scal(j)
c
               do i = iUpper, iLower
                  Jac(i+k,j) = tmp * Jac(i+k,j) / Scal(i)
               end do
c
            end do
c
         end if
c
c-----------------------------------------------------------------------
c
c     Generate Postscript plot of the Jacobian matrix if required.
c
c-----------------------------------------------------------------------
c
         if (      ( JacP .eq. - 1   .and. nStep .eq. 0 )
     2        .or. ( JacP .eq. nStep .and. iCIV .eq. 0 
     3                               .and. nStep .gt. 0 ) ) then
c
            call JacPlot ( n, 0, Jac, LDJac, ml, mu, 0, 0, nStep, 
     2                     Full_or_Band )
c
            JacP = 0
c
         end if
c
c-----------------------------------------------------------------------
c
c     Otherwise the now used Jacobian is not the matrix of the current 
c     derivatives.
c
c-----------------------------------------------------------------------
c
      else
c
         JacAct = 0
c
      end if
c
c-----------------------------------------------------------------------
c
c     Scale f(y(0)).
c
c-----------------------------------------------------------------------
c
      if ( iCIV .eq. 1 ) call dcopy ( n, yk, 1, dz, 1 )
c
      do i = 1, n
         dz(i) = dz(i) / Scal(i)
      end do  
c
c-----------------------------------------------------------------------
c
c     Label  for return  after a  stepsize reduction  with the current
c     Jacobian.
c
c-----------------------------------------------------------------------
c
  150 continue
c
      Direct = .false.
c 
      if ( hMax_Int .gt. zero ) then
         fcm = max ( fmin, abs ( h ) / hMax_Int )
      else
         fcm = fmin
      end if
c
c-----------------------------------------------------------------------
c
c     Check existence  of dense output  points in the current step and
c     save y if search is successful.
c
c-----------------------------------------------------------------------
c
      n_Dense = 0
c
      if ( iCIV .eq. 0 .and. Iopt(13) .eq. 1 .and. Iopt(14) .gt. 2 )
     2   then
c
         do i = 0, Iopt(14) - 3
c
            t_Dense = tn_Dense + dble ( i ) * delta_Dense
c
            if ( t_Dense .ge. t .and. t_Dense .lt. t + h ) then
c
               n_Dense = 1
c
               go to 160
c 
            end if
c
         end do
c
      end if
c
      if (         iCIV .eq. 0 .and. n_Dense .eq. 0 
     2     .and. (      ( Iopt(13) .eq. 2 .and. Iopt(14) .gt. 0 )
     3             .or. ( Iopt(13) .eq. 3 .and. h .gt. Ropt(2) ) 
     4             .or. ( Iopt(17) .eq. 1 .and. t_End .le. t + h ) ) )
     5   n_Dense = 1
c
  160 continue
c
      if ( n_Dense .eq. 1 ) call dcopy ( n, y, 1, Dense, 1 )
c
c-----------------------------------------------------------------------
c
c     Basic discretization loop.
c
c-----------------------------------------------------------------------
c
      if ( iCIV .eq. 0 ) then
c
         jmm = jmh
c
      else
c
         jmm = j_Cons
c
      end if
c
      do j = 1, jmm
c
         g     = h / dble ( nj(j) )
         g_inv = one / g
c
c-----------------------------------------------------------------------
c
c     Compute B(0) / h - A.
c
c-----------------------------------------------------------------------
c
         if ( Full_or_Band .eq. 0 ) then
c
            do i = 1, n
               call dcopy ( n, Jac(1,i), 1, B_Jac(1,i), 1 )
            end do
c
         else
c
            do i = 1, n
               do k1 = max(1,i-mu) - i + mb1, min(n,i+ml) - i + mb1
                  B_Jac(k1+ml,i) = Jac(k1,i)
               end do
            end do
c
         end if
c
c-----------------------------------------------------------------------
c
c     Add entries of B(t,y). 
c
c-----------------------------------------------------------------------
c
         do i = 1, nz
c
            irr = ir(i)
            icc = ic(i)
c
            if ( Full_or_Band .eq. 0 ) then
               irh = irr
            else
               irh = irr - icc + mb2
            end if
c
            B_Jac(irh,icc) = B_Jac(irh,icc)
     2                    + Scal(icc) * b0(i) / ( g * Scal(irr) )
c
         end do
c
c-----------------------------------------------------------------------
c
c     Semi-implicit Euler starting step.
c
c-----------------------------------------------------------------------
c
c     Copy y(0) to y(k).
c
c-----------------------------------------------------------------------
c
         call dcopy ( n,  y, 1,  yk, 1 )
c
c-----------------------------------------------------------------------
c
c     Factorize B(0) / h - A. 
c   
c-----------------------------------------------------------------------
c
         if ( Full_or_Band .eq. 0 ) then
c
            call dgetrf ( n, n, B_Jac, LDBJac, Pivot, IFail(2) )
c
         else
c
            call dgbtrf ( n, n, ml, mu, B_Jac, LDBJac, Pivot, IFail(2) )
c
         end if
c
         nDec = nDec + 1
c
c-----------------------------------------------------------------------
c
c     A non-zero  error code from dgetrf or dgbtrf greater 0 signals a 
c     singular matrix. In  this case, an  empirical stepsize reduction 
c     is tried. All other values of IFail(2) are irrecoverable errors. 
c     If an error occurs during the computation of CIVs, the algorithm 
c     ends with an error exit.
c
c-----------------------------------------------------------------------
c
         if ( IFail(2) .lt. 0 ) go to 410
c
         if ( IFail(2) .gt. 0 ) then
c
            if ( iCIV .eq. 0 ) go to 260
c
            go to 460
c
         end if
c
c-----------------------------------------------------------------------
c
c     Solve  the linear  system applying the factorization from dgetrf 
c     or dgbtrf to h/nj(j) * f((y(0)).
c
c-----------------------------------------------------------------------
c
         call dcopy ( n, dz, 1, del, 1 )
c
         if ( Full_or_Band .eq. 0 ) then
c
            call dgetrs ( 'N', n, 1, B_Jac, LDBJac, Pivot, del, 
     2                    Max_Size, IFail(2) )
c
         else
c
            call dgbtrs ( 'N', n, ml, mu, 1, B_Jac, LDBJac, Pivot, del, 
     2                    Max_Size, IFail(2) )
c
         end if
c
         nSol = nSol + 1
c
c-----------------------------------------------------------------------
c
c     Check  for systems  with  index > 1 and/or inconsistent  initial
c     values.
c
c-----------------------------------------------------------------------
c
         if (       iCIV .eq. 0 .and. Iopt(5) .ne. 1 
     2        .and. j .le. 2    .and. nStep .le. 2 ) then
c
            if ( j .eq. 1 ) then
c
               delm = zero
c
               do i = 1, n
c
                  Temp3(i) = abs ( del(i) )
                  delm     = max ( delm, Temp3(i) )
c
               end do
c
            else    
c
               do i = 1, n
c
                  if ( Iopt(11) .eq. 1 ) rTol2 = 0.5d0 * rTol(i)
c
                  abs_del = abs ( del(i) ) - rTol2
c
                  if ( abs_del .gt. Small ) then
c
                     del_quot = ( Temp3(i) + rTol2 ) / abs_del 
c
                     if ( del_quot .lt. ex2 ) then
c
                        if (       Iopt(1) .eq. 2 .and. iRed .eq. 0 
     2                       .and. nStep .eq. 0  )
     3                     write ( MonOut, 9520 ) i
c
                        if (       Temp3(i) .gt. delm * thresh
     2                       .and. iRed. gt. 2 ) then
c
                           if ( del_quot .lt. ex1 ) go to 420 
c
                           go to 430 
c
                        end if
c
                     end if
c
                  end if
c
               end do  
c
            end if
c
         end if
c
c-----------------------------------------------------------------------
c
c     Semi-implicit Euler discretization.
c
c-----------------------------------------------------------------------
c
c     If Iopt(10) = 1 and j = 2 store  Delta(y(0)) for  the subsequent
c     computation of Theta for a possible Jacobian reuse.
c
c-----------------------------------------------------------------------
c
      if ( Iopt(10) .eq. 1 .and. j .eq. 2 .and. iCIV .eq. 0 ) then
c
         call dcopy ( n, del, 1, Temp2, 1 )
c
         Delta_y0 = dnrm2 ( n, Temp2, 1 )
c
      end if
c
c-----------------------------------------------------------------------
c
c     For CIV computation  store the  derivatives of  the differential 
c     variables.
c
c-----------------------------------------------------------------------
c
         if ( iCIV .eq. 1 ) then
c
            do i = 1, nz
               yk(ic(i)) = g * ys(ic(i))
            end do
c
         end if
c
c-----------------------------------------------------------------------
c
c     Save  intermediate y if dense output option is specified and the
c     current integration interval contains dense output points.
c
c-----------------------------------------------------------------------
c
         do i = 1, n
            Temp1(i) = Scal(i) * del(i)
         end do
c
         if ( n_Dense .eq. 1 ) then
c
            if ( j .eq. 1 ) then
               kl = 3
            else
               kl = ipt(j-1) + 1
            end if
c
            call dcopy ( n, Temp1, 1, Dense(1,kl), 1 )
c
         end if
c
c-----------------------------------------------------------------------
c
c     Compute the scaled approximation for y(x(0)+h).
c
c-----------------------------------------------------------------------
c
         do i = 1, n
            yk(i) = yk(i) + Temp1(i)
         end do
c
c-----------------------------------------------------------------------
c
c     Check positiveness of components, if required.
c
c-----------------------------------------------------------------------
c
         if ( iCIV .eq. 0 ) then
c
            do i = 1, n
c
               if ( IPos(i) .eq. 1 .and. yk(i) .lt. zero ) then
c
                  NonPosComp = i
c
                  go to 270 
c
               end if
c
            end do
c
         end if
c
c-----------------------------------------------------------------------
c
c     For CIV computation restore the derivatives of  the differential 
c     variables and estimate the convergence of the Newton iteration.
c
c-----------------------------------------------------------------------
c
         if ( iCIV .eq. 1 ) then
c
            do i = 1, nz 
               ys(ic(i)) = g_inv * yk(ic(i))
            end do
c
            do i = 1, nz
c
               yk(ic(i))  = y(ic(i))
               del(ic(i)) = zero
c
            end do
c
            do i = 1, n
c
               if ( IPos(i) .eq. 1 .and. yk(i) .lt. zero ) then
c
                  NonPosComp = i
c
                  go to 440
c
               end if
c
            end do

            CErr = dnrm2 ( n, del, 1 )
c
            CIV_Convergence = .false.
c
            if ( Iopt(11) .eq. 0 ) then
c
               do i = 1, n
                  if ( abs ( del(i) ) .gt. rTol(1) ) go to 170
               end do
c
            else
c
               do i = 1, n
                  if ( abs ( del(i) ) .gt. rTol(i) ) go to 170
               end do
c
            end if
c
            CIV_Convergence = .true.
c
            j_Cons = j
            k_Cons = 1
c
            go to 300
c
  170       continue
c
            if ( j .eq. 1 ) Cons_Err = CErr
c
         end if
c
c-----------------------------------------------------------------------
c
c     For k = 2,... compute higher order approximations.
c
c-----------------------------------------------------------------------
c
         if ( iCIV .eq. 0 ) then
c
            njj = nj(j)
c
         else
c
            njj = k_Cons
c
         end if
c
         do k = 1, njj - 1
c
            if ( iCIV .eq. 0 ) then
c
               th = t + dble ( k ) * g
c
            else
c
               th = t
c
            end if
c
c-----------------------------------------------------------------------
c
c     Label  for return, if the  Richardson iteration  method does not 
c     converge.
c
c-----------------------------------------------------------------------
c
  180       continue
c
c-----------------------------------------------------------------------
c
c    Compute the approximation f(y(k)) = f(y(0)+k*h).
c
c-----------------------------------------------------------------------
c 
            FcnInfo = 1
c
            call Fcn ( n, nz, th, yk, del, bk, ir, ic, FcnInfo )
c
            nFcn  = nFcn + 1
            nFcnE = nFcnE + 1
c
            if ( FcnInfo .lt. 0 ) go to 280 
c
            if ( iCIV .eq. 1 ) then
c
               do i = 1, nz
                  del(ir(i)) = del(ir(i)) - bk(i) * ys(ic(i))
               end do
c
            end if
c
c-----------------------------------------------------------------------
c
c     Scale f(y(k)).
c
c-----------------------------------------------------------------------
c
            do i = 1, n
               del(i) = del(i) / Scal(i)
            end do
c
c-----------------------------------------------------------------------
c
c     If j = 2 prepare Theta  computation for possible Jacobian matrix
c     reuse.
c
c-----------------------------------------------------------------------
c
            if ( Iopt(10) .eq. 1 .and. j .eq. 2 .and. iCIV .eq. 0 ) then
c
               do i = 1, n
                  Temp3(i) = zero
               end do
c
               do i = 1, nz
c
                  irr = ir(i)
                  icc = ic(i)
c
                  Temp3(irr) = Temp3(irr) 
     2              + Scal(icc) * bk(i) * Temp2(icc) / ( g * Scal(irr) )
c
               end do
c
            end if
c
c-----------------------------------------------------------------------
c
c     Direct solution of
c     y(k+1) - y(k) = ( B(k) / h - A )**(-1) * f(y(k)).
c
c-----------------------------------------------------------------------
c
            if ( Direct ) then
c
c-----------------------------------------------------------------------
c
c     Build matrix B(k) / h - A.
c
c-----------------------------------------------------------------------
c
               if ( Full_or_Band .eq. 0 ) then
c
                  do i = 1, n
                     call dcopy ( n, Jac(1,i), 1, B_Jac(1,i), 1 )
                  end do
c
               else
c
                  do i = 1, n
                     do k1 = max(1,i-mu) - i + mb1, 
     2                       min(n,i+ml) - i + mb1
                        B_Jac(k1+ml,i) = Jac(k1,i)
                     end do
                  end do
c
               end if
c
c-----------------------------------------------------------------------
c
c     Add entries of B(t,y).
c
c-----------------------------------------------------------------------
c
               do i = 1, nz
c
                  irr = ir(i)
                  icc = ic(i)
c
                  if ( Full_or_Band .eq. 0 ) then
                     irh = irr
                  else
                     irh = irr - icc + mb2
                  end if
c
                  B_Jac(irh,icc) =  B_Jac(irh,icc) 
     2                           + Scal(icc) * bk(i) / ( g * Scal(irr) )
c
               end do
c
c-----------------------------------------------------------------------
c
c     Factorize B(k) / h - A.
c
c-----------------------------------------------------------------------
c
               if ( Full_or_Band .eq. 0 ) then
c
                  call dgetrf ( n, n, B_Jac, LDBJac, Pivot, IFail(2) )
c
               else
c
                  call dgbtrf ( n, n, ml, mu, B_Jac, LDBJac, Pivot, 
     2                          IFail(2) )
c
               end if
c
               nDec = nDec + 1
c
c-----------------------------------------------------------------------
c
c     A non-zero  error code from dgetrf or dgbtrf greater 0 signals a
c     singular matrix. In  this case, an  empirical stepsize reduction
c     is tried. All other values of IFail(2) are irrecoverable errors.
c     If an error occurs during the computation of CIVs, the algorithm
c     ends with an error exit.
c
c-----------------------------------------------------------------------
c
               if ( IFail(2) .lt. 0 ) go to 410
c
               if ( IFail(2) .gt. 0 ) then 
c
                  if ( iCIV .eq. 0 ) go to 260
c
                  go to 460
c
               end if
c
c-----------------------------------------------------------------------
c
c     Solve the  linear system applying the factorization  from dgetrf 
c     or dgbtrf to h/nj(j) * f((y(k)).
c
c-----------------------------------------------------------------------
c
               if ( Full_or_Band .eq. 0 ) then
c
                  call dgetrs ( 'N', n, 1, B_Jac, LDBJac, Pivot, del,
     2                          Max_Size, IFail(2) )
c
                else
c
                  call dgbtrs ( 'N', n, ml, mu, 1, B_Jac, LDBJac, Pivot,
     2                          del, Max_Size, IFail(2) )
c
               end if
c
               nSol = nSol + 1
c
            end if
c
c-----------------------------------------------------------------------
c
c     Iterative solution of
c     y(k+1) - y(k) = ( B(k) / h - A )**(-1) * f(y(k)).
c
c-----------------------------------------------------------------------
c
            if ( .not. Direct ) then
c
c-----------------------------------------------------------------------
c
c     Compute  first ( B(0) / h - A )**(-1) * f(y(k)) using  the known
c     decomposition of B(0) / h - A.
c 
c-----------------------------------------------------------------------
c
               if ( Full_or_Band .eq. 0 ) then
c
                  call dgetrs ( 'N', n, 1, B_Jac, LDBJac, Pivot, del,
     2                          Max_Size, IFail(2) )
c
               else
c
                  call dgbtrs ( 'N', n, ml, mu, 1, B_Jac, LDBJac, Pivot,
     2                          del, Max_Size, IFail(2) )
c
               end if
c
               nSol = nSol + 1
c
c-----------------------------------------------------------------------
c
c     Start of the iterative Richardson method.
c
c-----------------------------------------------------------------------
c
               Iter_Rich = 0
c
  190          continue
c
               Iter_Rich = Iter_Rich + 1
c
c-----------------------------------------------------------------------
c
c     Computation of 
c     delta(B(k)) * delta(k+1) =
c                          ( B(y(k)) - B(0) ) * ( y(k+1) - y(k) ) / h.
c
c-----------------------------------------------------------------------
c
               if ( Iter_Rich .eq. 1 ) then
c
                  call dcopy ( n, del, 1, Temp2, 1 )
c
               else
c
                  call dcopy ( n, Temp1, 1, Temp2, 1 )
c
               end if
c
               do i = 1, n 
                  Temp1(i) = zero
               end do
c
               do i = 1, nz
c
                  icc = ic(i)
                  irr = ir(i)
c
                  Temp1(irr) = Temp1(irr) 
     2                   + g_inv * Temp2(icc) * ( b0(i) - bk(i) ) 
     3                           * Scal(icc) / Scal(irr)
c
               end do
c
c
c-----------------------------------------------------------------------
c
c     Check norm of scaled difference ( B(y(k)) - B(0) ) / h.
c
c-----------------------------------------------------------------------
c
               if ( dnrm2 ( nz, Temp1, 1 ) .lt. eps_Mach10 ) go to 210
c
c-----------------------------------------------------------------------
c
c     Solution of the system
c     ( B(0) / h - A )**(-1) * ( delta(B(k) / h ) * delta(k+1).
c
c-----------------------------------------------------------------------
c
               if ( Full_or_Band .eq. 0 ) then
c
                  call dgetrs ( 'N', n, 1, B_Jac, LDBJac, Pivot, Temp1,
     2                          Max_Size, IFail(2) )
c
               else
c
                  call dgbtrs ( 'N', n, ml, mu, 1, B_Jac, LDBJac, Pivot,
     2                          Temp1, Max_Size, IFail(2) )
c
               end if
c
               nSol = nSol + 1
c
c-----------------------------------------------------------------------
c
c     Error estimation of the Richardson iteration.
c
c-----------------------------------------------------------------------
c
               do i = 1, n
                  del(i) = del(i) + Temp1(i)
               end do
c
               if ( Iopt(11) .eq. 0 ) then
c
                  tmp = 1.0d-1 * rTol(1)
c
                  do i = 1, n
                     if ( abs ( Temp1(i) ) .gt. tmp ) go to 200
                  end do
c
               else
c
                  do i = 1, n
                     if ( abs ( Temp1(i) ) .gt. 1.0d-1 * rTol(i) ) 
     2                  go to 200
                  end do
c
               end if
c
c-----------------------------------------------------------------------
c
c     Error of  Richardson iteration is small enough to  guarantee the
c     required accuracy.
c
c-----------------------------------------------------------------------
c
               go to 210
c
c-----------------------------------------------------------------------
c
c     Otherwise check iteration counter, if less than 10, continue the
c     iteration.
c
c-----------------------------------------------------------------------
c
  200          continue
c
               if ( Iter_Rich .lt. Max_Iter_Rich ) go to 190
c
c-----------------------------------------------------------------------
c
c     Richardson iteration does not converge, try direct method.
c
c-----------------------------------------------------------------------
c
               Direct = .true.
c
               if ( Iopt(1) .eq. 2 ) write ( MonOut, 9530 )
c
               go to 180
c
            end if
c
c-----------------------------------------------------------------------
c
c     Estimate contractivity factor Theta and determine, whether a new
c     Jacobian is required in the next step.
c
c-----------------------------------------------------------------------
c
  210       continue
c
            if ( Iopt(10) .eq. 1. and. j .eq. 2 .and. iCIV .eq. 0 ) then
c
               if ( Full_or_Band .eq. 0 ) then
c
                  call dgetrs ( 'N', n, 1, B_Jac, LDBJac, Pivot, Temp3,
     2                          Max_Size, IFail(2) )
c
               else
c
                  call dgbtrs ( 'N', n, ml, mu, 1, B_Jac, LDBJac, Pivot,
     2                          Temp3, Max_Size, IFail(2) )
c
               end if
c
               nSol = nSol + 1
c
               tmp = zero
c
               do i = 1, n
                  tmp = tmp + ( del(i) - Temp3(i) ) ** 2
               end do
c
               if ( sqrt ( tmp ) .gt. Delta_y0 * ThMin ) then
c
                  NewJac = 1
c
               else
c
                  NewJac = 0
c
               end if
c
            end if
c
c-----------------------------------------------------------------------
c
c     For CIV  computation store  the derivatives  of the differential 
c     variables.
c
c-----------------------------------------------------------------------
c
            if ( iCIV .eq. 1 ) then
c
               if ( j .eq. 1 ) then
c
                  call dcopy ( n, ys, 1, Temp1, 1 )
                  call dcopy ( n, yk, 1, Temp2, 1 )
c
               end if
c
               do i = 1, nz
                  yk(ic(i)) = g * ys(ic(i))
               end do
c
            end if
c
c-----------------------------------------------------------------------
c
c     Save  intermediate y if dense output option is specified and the
c     current integration interval contains dense output points.
c
c-----------------------------------------------------------------------
c
            do i = 1, n
               Temp3(i) = Scal(i) * del(i)
            end do
c
            if ( n_Dense .eq. 1 )
     2         call dcopy ( n, Temp3, 1, Dense(1,ipt(j-1)+k+1), 1 )
c
c-----------------------------------------------------------------------
c
c     Compute the scaled approximation for y(x(0)+k*h).
c
c-----------------------------------------------------------------------
c
            do i = 1, n
               yk(i) = yk(i) + Temp3(i)
            end do
c
c-----------------------------------------------------------------------
c
c     Check positiveness of components, if required.
c
c-----------------------------------------------------------------------
c
            if ( iCIV .eq. 0 ) then
c
               do i = 1, n
c
                  if ( IPos(i) .eq. 1 .and. yk(i) .lt. zero ) then
c
                     NonPosComp = i
c
                     go to 270
c
                  end if
c
               end do
c
            end if
c
c-----------------------------------------------------------------------
c
c     For CIV computation  restore the derivatives of the differential 
c     variables and estimate the convergence of the Newton iteration.
c
c-----------------------------------------------------------------------
c
            if ( iCIV .eq. 1 ) then
c
               do i = 1, nz
                  ys(ic(i)) = g_inv * yk(ic(i))
               end do
c
               do i = 1, nz
c
                  yk(ic(i))  = y(ic(i))
                  del(ic(i)) = zero
c
               end do
c
               do i = 1, n
c
                  if ( IPos(i) .eq. 1 .and. yk(i) .lt. zero ) then
c
                     NonPosComp = i
c
                     go to 440
c
                  end if
c
               end do
c
               CErr = dnrm2 ( n, del, 1 )
c
               CIV_Convergence = .false.
c
               if ( Iopt(11) .eq. 0 ) then
c
                  do i = 1, n
                     if ( abs ( del(i) ) .gt. rTol(1) ) go to 220
                  end do
c
               else
c
                  do i = 1, n
                     if ( abs ( del(i) ) .gt. rTol(i) ) go to 220
                  end do
c
               end if
c
               CIV_Convergence = .true.
c
               j_Cons = j
               k_Cons = k
c
               go to 300
c
  220          continue
c
               if ( j .eq. 1 ) then
c
                  if ( CErr .lt. Cons_Err ) then
c
                     Cons_Err = CErr
c
                  else
c
                     k_Cons = k - 1
                     j_Cons = k
c
                     call dcopy ( n, Temp1, 1, ys, 1 )
                     call dcopy ( n, Temp2, 1, yk, 1 )
c
                     go to 230
c
                  end if
c
               end if
c
            end if
c
c-----------------------------------------------------------------------
c
c     End of the Euler discretization loop.
c
c-----------------------------------------------------------------------
c
         end do   
c
c-----------------------------------------------------------------------
c
c     Extrapolation of y'(0) as ( y(k) - y(k-1) ) / h with k = nj(j).
c
c-----------------------------------------------------------------------
c
  230    continue
c
c-----------------------------------------------------------------------
c
c     Within CIV computation consider only the differential variables.
c
c-----------------------------------------------------------------------
c
         if ( iCIV .eq. 1 ) then
c
            do i = 1, n
               del(i) = zero
            end do
c
            do i = 1, nz
               del(ic(i)) = g * ys(ic(i)) / Scal(ic(i))
            end do
c
         end if
c
c-----------------------------------------------------------------------
c
c     Recursive extrapolation process for the derivatives.
c
c-----------------------------------------------------------------------
c
         call dscal ( n, g_inv, del, 1 )
c
         if ( j .eq. 1 ) then
c
            call dcopy ( n, del, 1, dtp, 1 )
c
         else
c
            do i = 1, n
               Temp1(i) = del(i) - dtp(i,1)
            end do
c
            call dcopy ( n, del, 1, dtp, 1 )
c
            do k = 2, j - 1
c
               d1 = d(j,j-k+1)
               d2 = one / ( d1 - one )
c
               do i = 1, n
c
                  Temp2(i) = d2 * Temp1(i)     
                  del(i)   = del(i) + Temp2(i)
                  Temp1(i) = d1 * Temp2(i) - dtp(i,k)
                  dtp(i,k) = Temp2(i)
c
               end do   
c
            end do
c
            call dscal ( n, one / ( d(j,1) - one ), Temp1, 1 )
c
            call dcopy ( n, Temp1, 1, dtp(1,j), 1 )
c
            do i = 1, n
               ys(i) = Scal(i) * ( del(i) + Temp1(i) )
            end do  
c
         end if
c
c-----------------------------------------------------------------------
c
c     Recursive extrapolation process for the solution.
c
c-----------------------------------------------------------------------
c
         if ( j .eq. 1 ) then
c
            call dcopy ( n, yk, 1, dt, 1 )
c
         else
c
            do i = 1, n
c
               del(i)   = yk(i) - dt(i,1)
               dt(i,1)  = yk(i)
               Temp1(i) = zero
c
            end do
c
            do k = 2, j - 1
c
               d1 = d(j,j-k+1)
               d2 = one / ( d1 - one )
c
               do i = 1, n
c
                  Temp2(i) = d2 * del(i)
                  Temp1(i) = Temp1(i) + Temp2(i)
                  del(i)   = d1 * Temp2(i) - dt(i,k)
                  dt(i,k)  = Temp2(i)
c
               end do   
c
            end do
c
            call dscal ( n, one / ( d(j,1) - one ), del, 1 )
c
            call dcopy ( n, del, 1, dt(1,j), 1 )
c
            do i = 1, n
               yk(i) = yk(i) + ( Temp1(i) + del(i) )
            end do
c
c-----------------------------------------------------------------------
c
c     Check positiveness of components, if required.
c
c-----------------------------------------------------------------------
c
            do i = 1, n
c
               if ( IPos(i) .eq. 1 .and. yk(i) .lt. zero ) then  
c
                  if ( iCIV .eq. 0 ) then 
                     go to 270
                  else
                     go to 440
                  end if
c
                end if
c
            end do
c
c-----------------------------------------------------------------------
c
c     Error estimation by the sub-diagonal differences. 
c
c-----------------------------------------------------------------------
c
            if ( iCIV .eq. 0 ) then
c
               fc          = zero
               Convergence = .true.
c
               if ( Iopt(11) .eq. 0 ) then
c
                  do i = 1, n
c
                     fc = max ( fc, abs ( del(i) ) 
     2                              / ( abs ( yk(i) ) + TolQuot ) )
c
                  end do
c
                  if ( fc .gt. rTol(1) ) Convergence = .false.
c
               else
c
                  do i = 1, n
c
                     tmp =   abs ( del(i) ) 
     2                     / ( abs ( yk(i) ) + TolQu(i) ) 
c
                     fc = max ( fc, tmp )
c
                     if ( tmp .gt. rTol(i) ) Convergence = .false.
c
                  end do
c
               end if
c
            end if
c
         end if
c
         if ( iCIV .eq. 1 ) go to 240
c
c-----------------------------------------------------------------------
c
c     Determine the optimal order.
c
c-----------------------------------------------------------------------
c
         if ( j .gt. 1 ) then
c
            k   = j - 1
            fc  = max ( fcm, ( fc / eph ) ** ( one / dble ( j ) ) )
            omj = fc * aj(j)
c
            if (       ( j .le. 2 .or. omj * one1 .lt. omjo ) 
     2           .and. k .le. joh ) then
c
               ko   = k
               jo   = j
               omjo = omj
               fco  = fc
c
            end if
c
c-----------------------------------------------------------------------
c
c     Possible increase of order.
c
c-----------------------------------------------------------------------
c
            if ( j .ge. koh .or. nStep. eq. 0 ) then
c
               if ( Convergence ) then
c
                  if ( ko .ge. k .and. incr(j) .ge. 0 ) then
c
                     fc = max ( fcm, fco / al(j,k) )
c
                     if ( aj(j+1) * fc * one1 .lt. omjo ) then
c
                        fco = fc
                        ko  = jo
                        jo  = jo + 1
c
                     end if
c
                  end if
c
                  go to 300
c
               end if
c
c-----------------------------------------------------------------------
c
c     Convergence monitor.
c
c-----------------------------------------------------------------------
c
               red = one / fco
               jk  = min ( km, joh ) 
c
               if ( k .ge. jk)  go to 250 
c
               if ( ko .lt. koh ) red = al(koh,ko) * red
c
               if ( al(jk,ko) .lt. fco ) go to 250 
c
            end if
c
         end if
c
  240    continue
c
      end do   
c
c----------------------------------------------------------------------
c
c     This is the end of the basic discretization loop.
c
c-----------------------------------------------------------------------
c
c     Stepsize reduction due to extrapolation tableau.
c
c-----------------------------------------------------------------------
c
      if ( iCIV .eq. 1 ) go to 300
c
  250 continue
c
      red = min ( red * safe, rmin )
      h   = h * red
C
      if ( nStep .gt. 0 ) incr(koh) = - 2
c
      iRed = iRed + 1
c
      if ( Iopt(1) .eq. 2 ) write ( MonOut, 9540 ) iRed, red
c
      if ( iRed .gt. Max_Step_Red_Ex ) go to 450 
c
      go to 290
c
c-----------------------------------------------------------------------
c
c     Empirical stepsize reduction due to zero pivot.
c
c-----------------------------------------------------------------------
c
  260 continue
c
      if ( Iopt(5) .eq. 0 ) iSing = iSing + 1 
c
      if ( iSing .gt. Max_Step_Red_Dc ) go to 460 
c
      hMax_Int = g * quart * safe
      red      = hMax_Int / abs ( h )
      h        = hMax_Int 
c
      if ( iRed .eq. 0 .and. nStep .gt. 0 ) incr(koh) = - 2
c
      iRed = iRed + 1
c
      if ( Iopt(1) .eq. 2 ) write ( MonOut, 9550 ) iRed, red
c  
      if ( iRed .gt. Max_Step_Red_Ex ) go to 450 
c
      go to 290
c
c-----------------------------------------------------------------------
c
c     Empirical stepsize reduction due to non positive components.
c
c-----------------------------------------------------------------------
c
  270 continue
c
      h    = h * RedPos
      iRed = iRed + 1
c
      if ( Iopt(1) .eq. 2 ) 
     2   write ( MonOut, 9560 ) iRed, NonPosComp, RedPos
c
      if ( iRed .gt. Max_Step_Red_Ex ) go to 450
c
      go to 290
c
c-----------------------------------------------------------------------
c
c     Empirical stepsize reduction due to an internal error in Fcn.
c
c-----------------------------------------------------------------------
c
  280 continue
c
      h    = h * RedFcn
      iRed = iRed + 1
c
      if ( Iopt(1) .eq. 2 ) 
     2   write ( MonOut, 9570 ) iRed, FcnInfo, RedFcn
c
      if ( iRed .gt. Max_Step_Red_Ex ) go to 450
c
c-----------------------------------------------------------------------
c
c     If the  currently used Jacobian is not the matrix of the current
c     derivatives, rescale the solution and  restore y', then  perform
c     rescaling.
c
c-----------------------------------------------------------------------
c
  290 continue
c
      if ( JacAct .eq. 0 ) then
c
         do i = 1, n
            dz(i) = Scal(i) * dz(i)
         end do
c
         call dcopy ( n, Temp4, 1, ys, 1 )
c
         if ( Iopt(11) .eq. 0 ) then
c
            do i = 1, n
               Scal(i) = abs ( y(i) ) + TolQuot
            end do
c
         else
c
            do i = 1, n
               Scal(i) = abs ( y(i) ) + TolQu(i)
            end do
c
         end if
c
         HoldJac = 0
c
         go to 140
c
      else
c
         go to 150
c
      end if
c
c-----------------------------------------------------------------------
c
c     Compute the Hermite  interpolation polynomial and evaluate it if
c     the current integration interval contains dense output points.
c
c-----------------------------------------------------------------------
c
  300 continue
c
      if ( n_Dense .eq. 0 ) go to 350
c
      call dcopy ( n, yk, 1, Dense(1,2), 1 )
c
      call Comp_Herm ( n, Max_Size, Dense, k, ipt, nj, dt )
c
      k_Dense = 0
c
c-----------------------------------------------------------------------
c
c     Continuation label for one step mode on  dense output points and
c     continuation calls.
c
c-----------------------------------------------------------------------
c
  310 continue
c
c-----------------------------------------------------------------------
c
c     Dense  output on equidistant points within the whole integration
c     interval. Return to  the calling  program if desired in the  one
c     step mode or within a continuation call.
c
c-----------------------------------------------------------------------
c
      if ( Iopt(13) .eq. 1 .and. Iopt(14) .gt. 1 ) then
c
         if ( Status .eq. 3 ) then
c
            delta_Dense = ( t_End - t_Begin ) / dble ( Iopt(14) - 1 )
            tn_Dense    = t_Begin + delta_Dense
            k_Dense     = 0
c
         end if

  320    continue
c
         Status  = 2
         k_Dense = k_Dense + 1
c
         t_Dense = tn_Dense + dble ( k_Dense - 1 ) * delta_Dense
c
         if (       t_Dense .lt. t + h
     2        .and. t_End - t_Dense .gt. abs ( t ) * eps_Mach10 )
     3      then
c
            call Eval_Herm ( n, Max_Size, Dense, k, (t_Dense-t)/h, 
     2                       Temp1 )
c
            if ( DensOut .ne. 0 ) write ( DensOut, 9480 ) 
     2         t_Dense, ( Temp1(j), j = 1, n )
c
            if ( Iopt(12) .eq. 2 ) then
c
                call dcopy ( n, Temp1, 1, y, 1 )
c
                t_Begin = t_Dense
c
                Iopt(24) = nFcnE
                Iopt(25) = nFcnJ
                Iopt(26) = nDec
                Iopt(27) = nSol
                Iopt(28) = nStep
                Iopt(29) = nJac
c
                return
c
            end if
c
            go to 320
c
         else
c
            Status   = 1
            tn_Dense = t_Dense
c
         end if
c
      end if
c
c-----------------------------------------------------------------------
c
c     Dense output on equidistant points within every integration step
c     or at least at a minimum distance. Return to the calling program
c     if desired in the one step mode or within a continuation call.
c
c-----------------------------------------------------------------------
c
      if (      ( Iopt(13) .eq. 2 .and. Iopt(14) .gt. 2 )
     2     .or.   Iopt(13) .eq. 3 ) then
c
         if ( Status .eq. 1 ) then
c
            Status   = 2
            h_True   = min ( h, t_End - t )
            tn_Dense = t
c
            if ( Iopt(13) .eq. 2 ) then
               n_Dense = Iopt(14) - 2
            else
               n_Dense = int ( h_True / Ropt(2) )
            end if
c
            delta_Dense = h_True / dble ( n_Dense + 1 )
c
         end if
c
  330    continue
c
         k_Dense = k_Dense + 1
c
         if ( k_Dense .gt. n_Dense ) go to 340
c
         t_Dense = tn_Dense + dble ( k_Dense ) * delta_Dense 
c
         call Eval_Herm ( n, Max_Size, Dense, k, (t_Dense-t)/h,
     2                    Temp1 )
c
         if ( DensOut .ne. 0 ) 
     2      write ( DensOut, 9480 ) t_Dense, ( Temp1(j), j = 1, n )
c
         if ( Iopt(12) .eq. 2 ) then
c
             call dcopy ( n, Temp1, 1, y, 1 )
c
             t_Begin = t_Dense
c
             Iopt(24) = nFcnE
             Iopt(25) = nFcnJ
             Iopt(26) = nDec
             Iopt(27) = nSol
             Iopt(28) = nStep
             Iopt(29) = nJac
c
             return
c
         end if
c
         go to 330
c
  340    continue
c
         Status = 1
c
      end if
c
c-----------------------------------------------------------------------
c
c     Within continuation calls and t + h > t_End evaluate solution at
c     t_End and return to the calling program.
c
c-----------------------------------------------------------------------
c
      if ( t + h .ge. t_End .and. Iopt(17) .eq. 1 ) then 
c
         call Eval_Herm ( n, Max_Size, Dense, k, (t_End-t)/h, Temp1 )
c
         if (        DensOut .ne. 0 
     2        .and. ( Iopt(13) .eq. 1 .or. Iopt(13) .eq. 3 ) ) 
     3      write ( DensOut, 9480 ) t_End, ( Temp1(j), j = 1, n )
c
         if ( Ropt(3) .le. t_End ) then
c
            nStep = nStep + 1
c
            if ( Iopt(1) .ge. 1 ) then
c
               if ( HoldJac .eq. 0  ) then
                  write ( MonOut, 9500 ) 
     2               nStep, '*', nFcn, t_End, h, k, ko
               else
                  write ( MonOut, 9500 ) 
     2               nStep, ' ', nFcn, t_End, h, k, ko
               end if
c
            end if
c
         end if
c
         if ( Iopt(3) .gt. 0 )
     2      write ( SolOut, 9510 ) t_End, ( y(i), i = 1, n )
c
         if ( Iopt(13) .eq. 3 ) then
c
            k_Dense     = 0
            n_Dense     = int ( ( t + h - t_End ) / Ropt(2) )
            delta_Dense = ( t + h - t_End ) / dble ( n_Dense + 1 )
            tn_Dense    = t_End
c
         end if
c
         call dcopy ( n, Temp1, 1, y, 1 )
c
         Status  = 3
         t_Begin = t_End
c
         Iopt(24) = nFcnE
         Iopt(25) = nFcnJ
         Iopt(26) = nDec
         Iopt(27) = nSol
         Iopt(28) = nStep
         Iopt(29) = nJac
c
         return
c
      end if
c
      Status = 1
c
c-----------------------------------------------------------------------
c
c     Preparations for the next basic integration step.
c
c-----------------------------------------------------------------------
c
  350 continue
c
      if ( iCIV .eq. 0 ) t = t + h 
c
      t_Begin = t
      tEnd_t  = t_End - t
      nStep   = nStep + 1
c
      call dcopy ( n, yk, 1, y, 1 )
c
      if ( nStep .gt. Max_Int_Steps ) go to 470
c
      if ( iCIV .eq. 1 ) then
c
         if ( CIV_Convergence ) then
c
            if ( Iopt(1) .gt. 0 ) write ( MonOut, 9580 ) 
     2         nStep, CErr, j_Cons, k_Cons
c
            iCIV  = 0
            jCIV  = 1
            h     = hSave
            k     = 0
            nStep = 0
c
            if ( Iopt(3) .gt. 0 ) then
c
               write ( SolOut, 9590 ) ( yk(i), i = 1, n ) 
               write ( SolOut, 9600 ) ( ys(i), i = 1, n ) 
c
            end if
c
         else 
c
            if ( nStep .gt. Max_Step_CIV ) go to 480
c
         end if
c
      else 
c
c-----------------------------------------------------------------------
c
c     Stepsize prediction.
c
c-----------------------------------------------------------------------
c
         h   = h / fco
         koh = ko
         joh = koh + 1
c
         h_proposed = h
c
      end if
c
c-----------------------------------------------------------------------
c
c     Stop, it stepsize is too small.
c
c-----------------------------------------------------------------------
c
      if ( abs ( h ) .lt. abs ( t ) * eps_Mach10 * ten ) go to 490
c
c-----------------------------------------------------------------------
c
c     Perform rescaling if a new Jacobian is needed in the next step.
c
c-----------------------------------------------------------------------
c
      if ( Iopt(10) .eq. 0 .or. NewJac .eq. 1 ) then
c
         if ( Iopt(11) .eq. 0 ) then
c
            do i = 1, n 
               Scal(i) = abs ( y(i) ) + TolQuot
            end do
c
         else
c
            do i = 1, n 
               Scal(i) = abs ( y(i) ) + TolQu(i)
            end do
c
         end if
c
         HoldJac = 0
c
      else
c
         HoldJac = 1
c
      end if
c
c-----------------------------------------------------------------------
c
c     This is the end of the integration step. If t_End is reached, go
c     to the solution exit.
c
c-----------------------------------------------------------------------
c
      if ( Iopt(1) .ge. 1 .and. jCIV .eq. 0 ) then
c
         if ( iCIV .eq. 0 ) then
c
            kc1 = k
            kc2 = ko
c
         else
c
            kc1 = j_Cons
            kc2 = k_Cons
c
         end if
c
         if ( HoldJac .eq. 0  ) then
            write ( MonOut, 9500 ) nStep, '*', nFcn, t, h, kc1, kc2
         else
            write ( MonOut, 9500 ) nStep, ' ', nFcn, t, h, kc1, kc2
         end if
c
      end if
c
      if (       (        Iopt(13) .gt. 1 
     2             .or. ( Iopt(13) .eq. 1 .and. nStep .eq. 0 ) )
     3     .and. DensOut .ne. 0 
     4     .and. ( Iopt(3) .eq. 0 .or. SolOut .ne. DensOut ) )
     5   write ( DensOut, 9480 ) t, ( y(i), i = 1, n )
c
      if ( Iopt(17) .ne. 1 ) then
c
         if (      abs( tEnd_t ) .lt. abs( t ) * eps_Mach10 * ten
     2        .or. abs( tEnd_t ) .lt. abs( h ) * eps_Mach10 )
     3      go to 360
c
      else if ( iCont .eq. 1 ) then
c
         if (      abs( Ropt(3) - t ) .lt. abs( t ) * eps_Mach10 * ten
     2        .or. abs( Ropt(3) - t ) .lt. abs( h ) * eps_Mach10 )
     3      go to 360
c
      end if

      if (         jCIV .eq. 0 
     2     .and. (        Iopt(12) .eq. 1 
     3             .or. (       Iopt(12) .eq. 2 .and. NoDense .eq. 1 
     4                    .and. Iopt(13) .ge. 2 ) ) ) then
c
         Iopt(24) = nFcnE
         Iopt(25) = nFcnJ
         Iopt(26) = nDec
         Iopt(27) = nSol
         Iopt(28) = nStep
         Iopt(29) = nJac
c
         return
c
      end if
c
      go to 130 
c
c-----------------------------------------------------------------------
c
c     Solution exit.
c
c-----------------------------------------------------------------------
c
  360 continue
c
      h = h_proposed
c   
      Iopt(24) = nFcnE
      Iopt(25) = nFcnJ
      Iopt(26) = nDec
      Iopt(27) = nSol
      Iopt(28) = nStep
      Iopt(29) = nJac 
c
      t_Begin = t_End
c
      if (         Iopt(13) .eq. 1 .and. DensOut .ne. 0
     2     .and. ( Iopt(3) .eq. 0 .or. SolOut .ne. DensOut ) ) 
     3   write ( DensOut, 9480 ) t_End, ( y(i), i = 1, n )
c
      if ( Iopt(3) .gt. 0 ) 
     2   write ( SolOut, 9510 ) t_End, ( y(i), i = 1, n )
c
      return 
c
c-----------------------------------------------------------------------
c
c     Integration error exit.
c
c-----------------------------------------------------------------------
c
  370 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9620 ) FcnInfo
c
      IFail(1) = - 33
c
      go to 500
c
  380 continue
c
      if ( Iopt(1) .gt. 0 )
     2   write ( MonOut, 9610 ) nz, Max_Non_Zeros_B
c
      IFail(1) = - 34
c
      go to 500
c
  390 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9630 ) FcnInfo
c
      IFail(1) = - 35
c
      go to 500
c
  400 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9640 ) JacInfo
c
      IFail(1) = - 36
c
      go to 500
c
  410 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9650 ) IFail(2)
c
      IFail(1) = - 39
c
      go to 500
c
  420 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9660 )
c
      IFail(1) = - 43
c
      go to 500
c
  430 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9670 )
c
      IFail(1) = - 44
c
      go to 500
c
  440 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9680 ) NonPosComp
c
      IFail(1) = - 45
c
      go to 500
c
  450 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9690 ) Max_Step_Red_Ex
c
      IFail(1) = - 46
c
      go to 500
c
  460 continue
c
      if ( Iopt(1) .gt. 0 ) then
c
         if ( iCIV .eq. 0 ) then
            write ( MonOut, 9700 ) Max_Step_Red_Dc
         else 
            write ( MonOut, 9710 )
         end if
c
      end if
c
      IFail(1) = - 47
c
      go to 500
c
  470 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9720 ) Max_Int_Steps
c
      IFail(1) = - 48
c
      go to 500
c
  480 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9730 ) Max_Step_CIV
c
      IFail(1) = - 49
c
      go to 500
c
  490 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9740 ) h
c
      IFail(1) = - 50
c
  500 continue
c
      if ( Iopt(1) .gt. 0 ) write ( MonOut, 9750 ) IFail(1)
c
      if ( iCIV .eq. 1 ) then
c
         IFail(3) = 1
c
         if ( Iopt(1) .gt. 0 ) write ( MonOut, 9760 ) 
c
      end if
c
      h = zero
c
      Iopt(24) = nFcnE
      Iopt(25) = nFcnJ
      Iopt(26) = nDec
      Iopt(27) = nSol
      Iopt(28) = nStep
      Iopt(29) = nJac
c
      Status = 0
c
      return
c
c-----------------------------------------------------------------------
c
c     Wrong input data exit.
c
c-----------------------------------------------------------------------
c
  510 continue
c
      if ( Iopt(1) .ne. 0 ) write ( MonOut, 9750 ) IFail(1)
c
      Status = 0
c
      return
c
c-----------------------------------------------------------------------
c
c     Format-statements
c
c-----------------------------------------------------------------------
c
 9000 format ( /, 16x, 43('-'), /, 16x, '|', 41x, '|', /, 16x, '|', 11x,
     2         'L I M E X  V. 4.2A1', 11x, '|', /, 16x, '|', 41x, '|', 
     3         /, 16x, '| Copyright (C) 2000, Konrad-Zuse-Zentrum |', /,
     4         16x, '|  fuer Informationstechnik Berlin (ZIB)  |', /,
     5         16x, 43('-'), / )
 9010 format ( ' *** LIMEX error, wrong input for Iopt(', i2, ') : ',
     2         i5, ',', 16x, '***' )
 9020 format ( ' *** should be between 0 and ', i5, 33x, '***' )
 9030 format ( ' *** should be greater or equal 0', 34x, '***' )
 9040 format ( ' *** LIMEX error, wrong input for n : ', i5, ',', 23x,
     2         '***', /, ' *** should be between 1 and ', i5, 33x, 
     3         '***' )
 9050 format ( ' *** should be greater or equal - 1', 32x, '***' )
 9060 format ( ' *** LIMEX error, wrong input for Ropt(', i1, ') : ',
     2         d18.9, 5x, '***' )
 9070 format ( ' *** should be greater 0', 43x, '***' )
 9080 format ( ' *** LIMEX error, wrong input for rTol(', i5, ') :',
     2         d18.9, ', ***', /, ' *** should be greater 0.0d0', 39x,
     3         '***' )
 9090 format ( ' *** LIMEX error, wrong input for aTol(', i5, ') :',
     2         d18.9, ', ***' )
 9100 format ( ' *** LIMEX error, non positive comp. y(', i5, ') :',
     2         d18.9, ' ***' )
 9110 format ( 9x, 57('-'), /, 9x, '|', 55x, '|', /, 9x, 
     2         '| Integration options:', 34x, '|', /, 9x, '|', 55x,
     3         '|' )
 9120 format ( 9x, '| Standard integration monitor on unit', i3, 15x, 
     2         '|' ) 
 9130 format ( 9x, '| Enhanced integration monitor on unit', i3, 15x,
     2         '|' )
 9140 format ( 9x, '| No output of solution values', 26x, '|' )
 9150 format ( 9x, '| Output of initial and solution values on unit',
     2         i3, 6x, '|' )
 9160 format ( 9x, '| Output of intermediate solution values on unit',
     2         i3, 5x, '|' )
 9170 format ( 9x, '| Number of equations: ', i5, 28x, '|' )
 9180 format ( 9x, '| Integration between ', d13.6, ' and ', d13.6, 3x,
     2         '|' )
 9190 format ( 9x, '| Matrix B may be singular', 30x, '|' )
 9200 format ( 9x, '| Matrix B is non-singular', 30x, '|' )
 9210 format ( 9x, '| No determination of consistent initial values',
     2         9x, '|' )
 9220 format ( 9x, '| Determination of consistent initial values',
     2         12x, '|' )
 9230 format ( 9x, '| Numerical difference approximation for the ',
     2         'Jacobian', 3x, '|' )
 9240 format ( 9x, '| Analytical Jacobian', 35x, '|' )
 9250 format ( 9x, '| Lower bandwidth of the Jacobian:', i5, 17x, '|', 
     2         /, 9x, '| Upper bandwidth of the Jacobian:', i5, 17x, 
     3         '|' )
 9260 format ( 9x, '| No reuse of the Jacobian', 30x, '|' )
 9270 format ( 9x, '| Reuse of the Jacobian, ThMin = ', d7.1, 16x, '|' )
 9280 format ( 9x, '| rTol, aTol used as scalars:', d10.3, 2x, d10.3,
     2         5x, '|' ) 
 9290 format ( 9x, '| rTol, aTol used as vectors', 28x, '|' ) 
 9300 format ( 9x, '| One step mode is off', 34x, '|' )
 9310 format ( 9x, '| One step mode is on', 35x, '|' )
 9320 format ( 9x, '| One step mode on, return at ', i4, ' points',
     2         15x, '|' )
 9330 format ( 9x, '| One step mode on, return at ', i4, ' points',
     2         ' per interval  |' )
 9340 format ( 9x, '| One step mode is on, return in steps <= ', d10.3,
     2         4x, '|' )
 9350 format ( 9x, '| No dense output', 39x, '|' )
 9360 format ( 9x, '| Dense output at ', i4, ' points on unit ', i3,
     2         15x, '|' )
 9370 format ( 9x, '| Dense output at ', i4, ' points per interval ',
     2         'on unit ', i3, '  |' )
 9380 format ( 9x, '| Dense output in steps <= ', d10.3,
     2         ' on unit ', i3, 7x, '|' )
 9390 format ( 9x, '| Integration for t > t_End disabled', 20x, '|' ) 
 9400 format ( 9x, '| Integration for t > t_End enabled', 21x, '|', /,
     2         9x, '|', 4x, 'upper t-limit: ', d12.5, 24x, '|' ) 
 9410 format ( 9x, '| Integration for t > t_End enabled', 21x, '|', /,
     2         9x, '|', 4x, 'no upper t-limit', 35x, '|' ) 
 9420 format ( 9x, '| No PostScript plot of the Jacobian', 20x, '|' )
 9430 format ( 9x, '| PostScript plot of the Jacobian in step ', i6, 
     2         7x, ' |' )
 9440 format ( 9x, '| Maximum stepsize: ', d12.5, 24x, '|', /, 9x, '|',
     2         55x, '|', /, 9x, 57('-'), / )  
 9450 format ( 9x, '| No upper limit for stepsize', 27x, '|', /, 9x, 
     2         57('-'), / )
 9460 format ( 9x, '*** Warning: This LIMEX call was specified as a',
     2         7x, '***', /, 9x, '*** continuation call, but there ',
     3         'was no call before.  ***', / )
 9470 format ( '   Initial solution values at t = ', d18.9, /,
     2         ( 1x, 4d18.9 ) )
 9480 format ( ' t = ', d18.9, /, ( 1x, 4d18.9 ) )
 9490 format ( /, 3x, 'Step-Nr.', 5x, 'F-calls', 9x, 't', 12x, 
     2         'Stepsize', 5x, 'used/aimed stage', /, 1x, 73('-'), / )
 9500 format ( 4x, i5, 1x, a1, i10, 2x, d15.5, 2x, d15.5, i9, i8 )
 9510 format ( /, 11x, 'Solution values at t = ', d18.9, /,
     2         ( 1x, 4d18.9 ) )
 9520 format ( /, ' Warning: component', i5, ' does not have an ',
     2         'asymptotic', /, 9x, ' expansion in the first ',
     3         'integration step', / )
 9530 format ( /, ' Warning: switch to direct solution of higher ',
     2         'order approximations', / )
 9540 format ( /, 4x, i3, '. Stepsize reduction (no conv.in extrap.) ',
     2         'by ', d10.3, / )
 9550 format ( /, 4x, i3, '. Stepsize reduction ( singular matrix  ) ',
     2         'by ', d10.3, / )
 9560 format ( /, 4x, i3, '. Stepsize reduction (y(', i5, 
     2         ' ) not pos.) by ', d10.3, / )
 9570 format ( /, 4x, i3, '. Stepsize reduction ( FcnInfo : ', i5,
     2         ' ) by ', d10.3, / )
 9580 format ( 4x, 'CIVs computed in ', i2, ' steps, scaled res.: ',
     2         d12.5, ', step/stage:', 2i2 )
 9590 format ( /, '   Consistent values', /, ( 1x, 4d18.9 ) )
 9600 format ( ' ... and derivatives', /, ( 1x, 4d18.9 ) )
 9610 format ( /, ' *** LIMEX error: ', i5, ' non-zero entries in ',
     2         'matrix B,', 14x, '***', /, ' *** upper limit is ',
     3         i5, 42x, '***' )
 9620 format ( /, ' *** LIMEX error: FcnInfo = ', i4, ' in first call ',
     2         'of Fcn in this step ***' )
 9630 format ( /, ' *** LIMEX error: FcnInfo = ', i4, ' during ',
     2         'evaluation of the Jacobian ***' )
 9640 format ( /, ' *** LIMEX error: JacInfo = ', i4, ' in the ',
     2         'analytical Jacobian', 8x, '***' )
 9650 format ( /, ' *** LIMEX error in routine dgetrf or dgbtrf, error',
     2         ' code: ', i3, 6x, '***' )
 9660 format ( /, ' *** LIMEX error: problem not solvable with LIMEX, ',
     2         16x, '***', /, ' ***', 14x, 'probably index of system ',
     3         'is > 1', 18x, '***' )
 9670 format ( /, ' *** LIMEX error: problem not solvable with LIMEX, ',
     2         16x, '***', /, ' ***', 14x, 'initial data inconsistent ',
     3         'or', 21x, '***', /, ' ***', 14x, 'index of the system ',
     4         'is > 1', 23x, '***' )
 9680 format ( /, ' *** LIMEX error: non-positive component y(', i5,
     2         ' )', 17x, '***' )
 9690 format ( /, ' *** LIMEX error: more than ', i4, ' stepsize ',
     2         'reductions', 15x, '***' )
 9700 format ( /, ' *** LIMEX error: more than ', i4, ' stepsize ',
     2         'reductions due to', 8x, '***', /, ' *** singular ',
     3         'matrix, system is not solvable with LIMEX', 12x, '***' )
 9710 format ( /, ' *** LIMEX error: singular Jacobian matrix', 25x, 
     2         '***' )
 9720 format ( /, ' *** LIMEX error: more than ', i6, ' integration ',
     2         'steps', 15x, '***' )
 9730 format ( /, ' *** LIMEX error: more than ', i3, ' Newton steps',
     2         ' for CIV-computation', 3x, '***' )
 9740 format ( /, ' *** LIMEX error: stepsize h = ', d18.9, ' is too ',
     2         'small', 5x, '***' )
 9750 format ( ' *** LIMEX error code is: ', i3, 38x, '***' )
 9760 format ( ' *** This error occurred during CIV-computation', 20x,
     2         '***' )
c
c-----------------------------------------------------------------------
c
c     End of LIMEX
c
c-----------------------------------------------------------------------
c
      end
c




      subroutine JacPlot ( n, LenJac, Jac, LDJac, ml, mu, ir, ic, nStep,
     2                     mode )
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     This subroutine generates a Postscript plot of the Jacobian.
c
c-----------------------------------------------------------------------
c
c     Written by:
c
c     R. Ehrig, U. Nowak,
c     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
c     Takustrasse 7
c     D-14195 Berlin-Dahlem
c
c     phone : +49-30-84185-0
c     fax   : +49-30-84185-125
c     e-mail: ehrig@zib.de, nowak@zib.de
c     URL   : http://www.zib.de
c
c-----------------------------------------------------------------------
c
c     Arguments:
c
c-----------------------------------------------------------------------
c
c     n          Integer variable, must be set by caller on input, not
c                modified. Size of the system.
c
c     LenJac     Integer variable, must be set by caller on input, not
c                modified. Number of  non-zero entries in the Jacobian
c                if specified by parameter mode as a sparse matrix.  
c
c     Jac        Real array of size (LDJac,*), must  be set by  caller
c                on  input only  if  mode = 0 or 1, not  modified. The 
c                values of the Jacobian matrix in full or banded form. 
c
c     LDJac      Integer variable, must be set by caller on input only
c                if mode = 0 or 1, not  modified. Leading dimension of 
c                array Jac.
c
c     ml         Integer variable, must be set by caller on input only
c                if mode = 1, not  modified. For  banded Jacobians the
c                number of lower diagonals.
c
c     mu         Integer variable, must be set by caller on input only
c                if mode = 1, not  modified. For  banded Jacobians the
c                number of upper diagonals.
c
c     ir         Integer array, must be set by caller on input only if
c                mode = 2, not modified. Row  indices for the non-zero
c                entries of the Jacobian in sparse form.
c
c     ic         Integer array, must be set by caller on input only if
c                mode = 2, not  modified. Column indices for the  non-
c                zero entries of the Jacobian is sparse form.
c
c     nStep      Integer variable, must be set by caller on input, not
c                modified. Current step number, which  is displayed in
c                the title.
c
c     mode       Integer variable, must be set by caller on input, not
c                modified. Determines the Jacobian matrix format.     
c                = 0 : dense Jacobian
c                = 1 : dense banded Jacobian
c                = 2 : sparse Jacobian in coordinate format
c
c-----------------------------------------------------------------------
c
c     Subroutines and functions called:
c
c     none
c
c-----------------------------------------------------------------------
c
c     Declaration of the arguments
c
c-----------------------------------------------------------------------
c
      integer            n, LenJac, LDJac, ml, mu, ir ( * ), ic ( * ),
     2                   nStep, mode
c
      double precision   Jac ( LDJac, * )
c
c-----------------------------------------------------------------------
c
c     Local scalars
c
c-----------------------------------------------------------------------
c
      integer            i, j
c
      real               Bottom, cm_dot, Left, Margin, Right, Scaling, 
     2                   Top, xTitle, yTitle 
c
c-----------------------------------------------------------------------
c
c     Internal parameters, sizes may be changed.
c
c-----------------------------------------------------------------------
c
      real               PaperSize, PlotSize
c
      parameter        ( PaperSize = 21.0, PlotSize = 15.0 )
c
c-----------------------------------------------------------------------
c
c     Define conversion  factor from  cm to  dot, margin  and scaling 
c     factor.
c 
c-----------------------------------------------------------------------
c
      cm_dot  = 72.0 / 2.54
      Margin  = 0.5 * ( PaperSize - PlotSize )
      Scaling = 0.125 * PlotSize * cm_dot / float ( n + 1 ) + 10.0
c
c-----------------------------------------------------------------------
c
c     Define title coordinates and bounding box.
c
c-----------------------------------------------------------------------
c
      xTitle = 0.5 * PaperSize
      yTitle = 3.0 + PlotSize
c
      Left   = Margin * cm_dot - Scaling
      Right  = ( Margin + PlotSize ) * cm_dot + Scaling
      Bottom = 2.0 * cm_dot - Scaling
      Top    = ( 3.35 + PlotSize ) * cm_dot + Scaling
c
c-----------------------------------------------------------------------
c
c     Open file and write postscript header information.
c
c-----------------------------------------------------------------------
c
      open ( 99, file = 'Jacobian.ps' )
c
      write ( 99, 9000 ) '%!'
      write ( 99, 9000 ) '%%Creator: JacPlot in LIMEX4.2'
      write ( 99, 9010 ) '%%BoundingBox:', Left, Bottom, Right, Top
      write ( 99, 9000 ) '%%EndComments'
      write ( 99, 9000 ) '/cm {72 mul 2.54 div} def'
      write ( 99, 9000 ) '/mc {72 div 2.54 mul} def'
      write ( 99, 9000 ) '/pnum { 72 div 2.54 mul 20 string'
      write ( 99, 9000 ) 'cvs print ( ) print} def'
      write ( 99, 9000 )
     2   '/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def'
      write ( 99, 9000 ) 'gsave'
      write ( 99, 9000 ) 
     2   '/Helvetica findfont 0.5 cm scalefont setfont'
      write ( 99, 9020 ) xTitle, ' cm ', yTitle, ' cm moveto'
      write ( 99, 9030 ) '(Jacobian matrix from LIMEX4.2, step = ', 
     2   nStep, ') Cshow'
      write ( 99, 9040 ) Margin, 'cm 2.0 cm translate'
      write ( 99, 9050 ) PlotSize, 'cm', n + 1, 'div dup scale'
c
c-----------------------------------------------------------------------
c
c     Draw frame.
c
c-----------------------------------------------------------------------
c
      write ( 99, 9000 ) '0.25 setlinewidth'
      write ( 99, 9000 ) 'newpath'
      write ( 99, 9000 ) '0 0 moveto'
      write ( 99, 9060 ) n + 1, 0, 'lineto'
      write ( 99, 9060 ) n + 1, n + 1,'lineto'
      write ( 99, 9060 ) 0, n + 1, 'lineto'
      write ( 99, 9000 ) 'closepath stroke'
c
c-----------------------------------------------------------------------
c
c     Prepare matrix plot.
c
c-----------------------------------------------------------------------
c
      write ( 99, 9000 ) '1 1 translate'
      write ( 99, 9000 ) '0.8 setlinewidth'
      write ( 99, 9000 ) '/p {moveto 0 -.40 rmoveto'
      write ( 99, 9000 ) '           0  .80 rlineto stroke} def'
c     
c-----------------------------------------------------------------------
c
c     If mode = 0 the matrix is dense and full.
c
c-----------------------------------------------------------------------
c
      if ( mode .eq. 0 ) then
c
         do i = 1, n
            do j = 1, n
               if ( Jac(j,i) .ne. 0.0d0 ) 
     2            write ( 99, 9060 ) i - 1, n - j, 'p'
            end do
         end do
c
c-----------------------------------------------------------------------
c
c     If mode = 1 the matrix is dense and banded.
c
c-----------------------------------------------------------------------
c
      else if ( mode .eq. 1 ) then
c
         do i = 1, n
            do j = max ( 1, i - mu ), min ( n, i + ml )
               if ( Jac(mu+1+j-i,i) .ne. 0.0d0 )
     2            write ( 99, 9060 ) i - 1, n - j, 'p'
            end do
         end do
c
c-----------------------------------------------------------------------
c
c     If mode = 2 the matrix is sparse.
c
c-----------------------------------------------------------------------
c
      else if ( mode .eq. 2 ) then
c
         do i = 1, LenJac
            write ( 99, 9060 ) ic(i) - 1, n - ir(i), 'p'
         end do
c
      end if
c
c-----------------------------------------------------------------------
c
c     Finalize and close postscript file.
c
c-----------------------------------------------------------------------
c
      write ( 99, 9000 ) 'showpage'
c
      close ( 99 )
c
c-----------------------------------------------------------------------
c
c     Format-statements
c
c-----------------------------------------------------------------------
c
 9000 format ( a )
 9010 format ( a, 4 ( 1x, f9.2 ) )
 9020 format ( 2 ( 1x, f9.2, a ) )
 9030 format ( a, i6, a )
 9040 format ( f9.2, 1x, a )
 9050 format ( f9.2, 1x, a, 1x, i6, 1x, a )
 9060 format ( 2 ( i6, 1x ), a )
c
c-----------------------------------------------------------------------
c
      return
      end
      subroutine Get_eps_Mach(eps_Mach)
      implicit none
      double precision eps_Mach
      double precision one,sum
      parameter (one=1.d0)

c      eps_Mach=one
c      sum=one+eps_Mach
c      do while (sum .ne. one)
c         eps_Mach = 0.5d0*eps_Mach
c         sum=one+eps_Mach
c      end do
      eps_Mach=1.d-16

      end

