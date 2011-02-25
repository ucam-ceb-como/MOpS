

c-----------------------------------------------------------------------
c
c     This file will  be included by LIMEX version 4.2. It defines the
c     size of the problem to  be solved  by LIMEX. It should  have the 
c     name: 
c
c     xlimex_def.h
c
c     All dimension statements in LIMEX  refer to the size definitions
c     in the  include file 'xlimex_def.h'. The file has to
c     consist of four parameter statements of the following form:
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
c     parameter is not used in the LIMEX versions 4.2A1 and 4.2A2.
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
c     This  parameter  is  not used  in the LIMEX  versions 4.2B1  and 
c     4.2B2
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
c     Max_It_Vectors >= Iopt(20) + 4 (this  relation  is  checked). If  
c     the Krylov subspace is exhausted without  attaining the required
c     accuracy, a  restart is done. Reasonable values for Iopt(20) are
c     in many cases 10 or 20, then Max_It_Vectors should  be set to 15
c     resp. 25. The greater  the dimension  of the Krylov subspace is, 
c     the  faster GMRES  converges  (in theory), but  needs a  rapidly 
c     increasing number of floating point  operations and is also less
c     stable. If  BICGSTAB is  used, Max_It_Vectors must be 6 at least 
c     (this  condition is  checked). If  no iterative  method will  be
c     used, Max_It_Vectors may  be set equal 1. This parameter  is not
c     used in the LIMEX versions 4.2A1 and 4.2A2.
c
c     parameter ( Max_It_Vectors         =
c
c-----------------------------------------------------------------------
c
c     Current settings for LIMEX distribution examples
c
c-----------------------------------------------------------------------
c
c     
c		include 'constants.F'
      parameter ( Max_Nr_of_Equations    = 1000)
      parameter ( Max_Non_Zeros_Jacobian = 1)
      parameter ( Max_Non_Zeros_B        = 1000)
      parameter ( Max_Lower_Diagonals    = 1000)
      parameter ( Max_Upper_Diagonals    = 1000)
      parameter ( Max_It_Vectors         = 25)
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
