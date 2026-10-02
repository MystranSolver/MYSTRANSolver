! ##################################################################################################################################

      MODULE ARPACK_LANCZOS_EIG

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EIG_MSGLVL
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, SLU_INFO
      USE PARAMS, ONLY                :  SOLLIB
      USE ARPACK_UTIL
      USE SUPERLU_ADAPTERS, ONLY       :  FBS_SUPRLU,                                  &
     &                                   SYM_MAT_DECOMP_SUPRLU

      USE DATE_TIME_UTILS, ONLY     : OURTIM
      USE SPARSE_FULL_MULTIPLICATION, ONLY: MATMULT_SFF
      USE ARPACK_DIAGNOSTICS, ONLY : ARPACK_INFO_MSG
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY: LINK_MESSAGE

      character(1*byte), parameter   :: cr13_a = char(13)
      CHARACTER(44*BYTE)             :: MODNAM1            ! Name to write to screen to describe module being run.
      CHARACTER(44*BYTE)             :: MODNAM2            ! Name to write to screen to describe module being run.

c
c\SCCS Information: @(#)
c FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2
c
c     %---------------------------------%
c     | See debug.doc for documentation |
c     %---------------------------------%
      integer :: logfil = 2
      integer :: ndigit = 6
      integer  mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/
     &         logfil, ndigit, mgetv0,
     &         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd

c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
c\SCCS Information: @(#)
c FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2
c
      real      t0, t1, t2, t3, t4, t5
      save      t0, t1, t2, t3, t4, t5
c
      integer   nopx, nbx, nrorth, nitref, nrstrt
      real      tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &          tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &          tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &          tmvopx, tmvbx , tgetv0, titref, trvec
      common /timing/
     &          nopx, nbx, nrorth, nitref, nrstrt,
     &          tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &          tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &          tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &          tmvopx, tmvbx , tgetv0, titref, trvec

! This is the set of ARPACK routines that are used in the Lanczos algorithm. Below are listed the subroutines included in this
! module and the calls to other subroutines in this module

!     dsband: the driver for Lanczos. dsband calls:

!             dsaupd: reverse communication interface, which calls:

!                     dstats: Initialize statistic and timing information

!                     dsaup2: Intermediate level interface called by dsaupd, which calls

!                             dgetv0: Generate a random initial residual vector for the Arnoldi process. Force the residual vector
!                                     to be in the range of the operator OP.

!                             dsaitr: Reverse comm I/F for applying NP additional steps to a K step symmetric Arnoldi factorization.

!                                     dgetv0: see above

!                             dseigt:
!                                     dstqrb:

!                             dsgets: Given the eigenvalues of the symmetric tridiagonal matrix H, computes the NP shifts AMU that
!                                     are zeros of the polynomial of  degree NP which filters out components of the unwanted
!                                     eigenvectors  corresponding to the AMU's based on some given criteria.

!                                     dsortr: Sort an array

!                             dsconv: Convergence testing for the symmetric Arnoldi eigenvalue routine.

!                             dsortr: Sort an array

!                             dsapps: Given the Arnoldi factorization
!                                     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
!                                     apply NP shifts implicitly resulting in
!                                     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
!                                     where Q is an orthogonal matrix of order KEV+NP. Q is the product of rotations resulting from
!                                     the NP bulge chasing sweeps.  The updated Arnoldi factorization becomes:
!                                     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.


!             dseupd: returns the converged approximations to eigenvalues of A*z = lambda*B*z and (optionally):

!                     dsgets: see above

!                     dsteqr: see above

!                     dsesrt: Sort an array

!                     dsortr: see above

      CONTAINS

! ##################################################################################################################################
! 001 ARPACK_LANCZOS_EIG

c \BeginDoc
c
c \Name: dsband
c
c \Description:
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal (Lanczos) basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  Matrices A and B are stored in LAPACK-style band form.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  (Lanczos) basis is always computed.  There is an additional storage cost
c  of n*nev if both are requested (in this case a separate array Z must be
c  supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are called Ritz values and Ritz vectors respectively.  They are referred
c  to as such in the comments that follow.  The computed orthonormal basis
c  for the invariant subspace corresponding to these Ritz values is referred
c  to as a Lanczos basis.
c
c  dsband can be called with one of the following modes:
c
c  Mode 2:  AB*x = lambda*MB*x, AB symmetric, MB symmetric positive definite
c           ===>  OP = inv[MB]*AB
!                 NOTE: ARPACK's implementation of this was intended for:
!                       AB = K (stiffness) and MB = M (mass) which required
!                       M symm pos def which it generally is not

! **********************************************************************************************************************************
! **********************************************************************************************************************************
! **                                                              |                                                               **
! **      N O R M A L   E I G E N V A L U E   P R O B L E M       |                      B U C K L I N G                          **
! **      -------------------------------------------------       |                      ---------------                          **
! ** For mode 2, I changed the roles of AB and MB and added a     | For mode 2, I changed the roles of AB and MB and added a      **
! ** shift so that:                                               | shift so that (KD is differential stiffness matrix):          **
! **                                                              |                                                               **
! **            AB = M                                            |            AB = -KD                                           **
! **            MB = [K - sigma*M]                                |            MB = [K + sigma*KD]                                **
! ** then                                                         | then                                                          **
! **            OP = (inv[K - sigma*M])*M                         |            OP = -(inv[K + sigma*KD])*KD                       **
! **                                                              |                                                               **
! ** The generalized eigenvalue problem is:                       | The generalized eigenvalue problem is:                        **
! **                                                              |                                                               **
! **                [K - lam*M]*x = 0                      (1)    |                [K + lam*KD]*x = 0                      (1)    **
! ** define                                                       | define                                                        **
! **               lam = lam_bar + sigma                          |               lam = lam_bar + sigma                           **
! ** Then (1) is:                                                 | Then (1) is:                                                  **
! **            [K - (lam_bar + sigma)*M]x = 0                    |            [K + (lam_bar + sigma)*KD]x = 0                    **
! ** or                                                           | or                                                            **
! **             [K - sigma*M]x = lam_bar*Mx               (2)    |             [K + sigma*KD]x = -lam_bar*KDx             (2)    **
! ** define                                                       | define                                                        **
! **                 mu = 1/lam_bar                               |                 mu = 1/lam_bar                                **
! ** so that (2) is:                                              | so that (2) is:                                               **
! **             mu*x = (inv[K - sigma*M])*M*x                    |            -mu*x = (inv[K + sigma*KD])*KD*x                   **
! ** or                                                           | or                                                            **
! **                   A*x = mu*x                          (3)    |                   A*x = mu*x                           (3)    **
! ** where                                                        | where                                                         **
! **              A = (inv[K - sigma*M])*M                        |              A = -(inv[K + sigma*KD])*KD                      **
! **                                                              |                                                               **
! ** eqn (3) is the standard eigenvalue problem solved where      | eqn (3) is the standard eigenvalue problem solved where       **
! ** the eigenvalues, lam, of the gen eigen problem are:          | the eigenvalues, lam, of the gen eigen problem are:           **
! **                                                              |                                                               **
! **                lam = sigma + 1/mu                            |                lam = sigma + 1/mu                             **
! **                                                              |                                                               **
! ** The calc of lambda = sigma + 1/mu is done in my subr         | The calc of lambda = sigma + 1/mu is done in my subr          **
! ** INVERT_EIGENS after dsband returns to my subr EIG_LANCZOS    | INVERT_EIGENS after dsband returns to my subr EIG_LANCZOS     **
! **                                                              |                                                               **
! **********************************************************************************************************************************
! **********************************************************************************************************************************


c
c  Mode 3:  AB*x = lambda*MB*x, AB symmetric, MB symmetric positive semi-definite
c           ===>  OP = (inv[K - sigma*M])*M  and  B = M.
c           ===>  Shift-and-Invert mode
!                 NOTE: dseupd does the calc: lambda = sigma + 1/mu
c
c  The choice of mode must be specified in IPARAM(7) defined below.
c
c \Usage
c     CALL DSBAND ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA,
c    &            N, LDA, RFAC, KL, KU, WHICH, BMAT, NEV,
c    &            TOL, RESID, NCV, V, LDV, IPARAM, WORKD, WORKL,
c    &            LWORKL, IWORK, INFO, INFO_LAPACK, DTBSV_MSG, PITERS)
c
c \Arguments
c
c  RVEC    Logical (INPUT)
c          Specifies whether Ritz vectors corresponding to the Ritz value
c          approximations to the eigenproblem A*z = lambda*B*z are computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the associated Ritz vectors.
c
c  HOWMNY  Character*1  (INPUT)
c          Specifies how many Ritz vectors are wanted and the form of Z
c          the matrix of Ritz vectors. See remark 1 below.
c          = 'A': compute all Ritz vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE..
c          If HOWMNY = 'A' , SELECT is not referenced.
c
c  D       REAL(DOUBLE) array of dimension NEV.  (OUTPUT)
c          On exit, D contains the Ritz value approximations to the
c          eigenvalues of A*z = lambda*B*z. The values are returned
c          in ascending order. If IPARAM(7) = 3,4,5 then D represents
c          the Ritz values of OP computed by dsaupd transformed to
c          those of the original eigensystem A*z = lambda*B*z. If
c          IPARAM(7) = 1,2 then the Ritz values of OP are the same
c          as the those of A*z = lambda*B*z.
c
c  Z       REAL(DOUBLE) N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c          On exit, Z contains the B-orthonormal Ritz vectors of the
c          eigensystem A*z = lambda*B*z corresponding to the Ritz
c          value approximations.
c
c          If  RVEC = .FALSE. then Z is not referenced.
c          NOTE: The array Z may be set equal to first NEV columns of the
c          Lanczos basis array V computed by DSAUPD.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
c
c  SIGMA   REAL(DOUBLE)  (INPUT)
c          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c          IPARAM(7) = 1 or 2.
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  AB      REAL(DOUBLE) array of dimension LDA by N. (INPUT)
c          The matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c  MB      REAL(DOUBLE) array of dimension LDA by N. (INPUT)
c          The matrix M in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of M is stored in the j-th column of the
c          array AB as follows:
c          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c          Not referenced if IPARAM(7) = 1
c
c  LDA     Integer. (INPUT)
c          Leading dimension of AB, MB, RFAC.
c
c  RFAC    REAL(DOUBLE) array of LDA by N. (WORKSPACE/OUTPUT)
c          RFAC is used to store the LU factors of MB when IPARAM(7) = 2
c          is invoked.  It is used to store the LU factors of
c          (A-sigma*M) when IPARAM(7) = 3,4,5 is invoked.
c          It is not referenced when IPARAM(7) = 1.
c
c  KL      Integer. (INPUT)
c          Max(number of subdiagonals of A, number of subdiagonals of M)
c
c  KU      Integer. (OUTPUT)
c          Max(number of superdiagonals of A, number of superdiagonals of M)
c
c  WHICH   Character*2.  (INPUT)
c          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
c          the following.
c
c            'LM' -> want the NEV eigenvalues of largest magnitude.
c            'SM' -> want the NEV eigenvalues of smallest magnitude.
c            'LA' -> want the NEV eigenvalues of largest REAL part.
c            'SA' -> want the NEV eigenvalues of smallest REAL part.
c            'BE' -> Compute NEV eigenvalues, half from each end of the
c                    spectrum.  When NEV is odd, compute one more from
c                    the high end than from the low end.
c
c          When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

c  NEV     Integer. (INPUT)
c          Number of eigenvalues of OP to be computed.
c
c  TOL     REAL(DOUBLE) scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          Represents the dimension of the Lanczos basis constructed
c          by dsaupd for OP.
c
c  V       REAL(DOUBLE) array N by NCV.  (OUTPUT)
c          Upon INPUT: the NCV columns of V contain the Lanczos basis
c                      vectors as constructed by dsaupd for OP.
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       represent the Ritz vectors that span the desired
c                       invariant subspace.
c          NOTE: The array Z may be set equal to first NEV columns of the
c          Lanczos basis vector array V computed by dsaupd. In this case
c          if RVEC=.TRUE., the first NCONV=IPARAM(5) columns of V contain
c          the desired Ritz vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT:
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          It is set to 1 in this subroutine.  The user do not need
c          to set this parameter.
c          ------------------------------------------------------------
c          ISHIFT = 1: exact shifts with respect to the reduced
c                      tridiagonal matrix T.  This is equivalent to
c                      restarting the iteration with a starting vector
c                      that is a linear combination of Ritz vectors
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced.
c
c          IPARAM(3) = MXITER
c          On INPUT:  max number of Arnoldi update iterations allowed.
c          On OUTPUT: actual number of Arnoldi update iterations taken.
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" eigenvalues.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used.
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of dsband for the
c          five modes available.
c
c          IPARAM(8) = NP
c          Not referenced.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c WORKD    REAL(DOUBLE) work array of length at least 3*n. (WORKSPACE)
c
c WORKL    REAL(DOUBLE) work array of length LWORKL.  (WORKSPACE)
c
c LWORKL   Integer.  (INPUT)
c          LWORKL must be at least NCV**2 + 8*NCV.
c
c IWORK    Integer array of dimension at least N. (WORKSPACE)
c          Used when IPARAM(7)=2,3,4,5 to store the pivot information in the
c          factorization of M or (A-SIGMA*M).
c
c INFO     Integer.  (INPUT/OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)
c                returns the number of wanted converged Ritz values.
c          =  3: No shifts could be applied during a cycle of the
c                Implicitly restarted Arnoldi iteration. One possibility
c                is to increase the size of NCV relative to NEV.
c                See remark 4 in DSAUPD.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informational error from LAPACK routine dsteqr.
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: NEV and WHICH = 'BE' are incompatible.
c          = -13: HOWMNY must be one of 'A' or 'P'
c          = -14: DSAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current
c                   Arnoldi factorization.
c
!  info_lapack (output) INTEGER: Error from LAPACK

!  dtbsv_msg (input) CHARACTER
!           = 'Y', have subr DTBSV print Fwd, Back pass messages

!  piters (input) This is debug(47) to tell whether to write eigen
!                 info at each iteration
c \EndDoc
c
c------------------------------------------------------------------------
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Ph.D thesis, TR95-13, Rice Univ,
c     May 1995.
c
c\Routines called:
c     dsaupd  ARPACK reverse communication interface routine.
c     dseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dgbtrf  LAPACK band matrix factorization routine.
c     dgbtrs  LAPACK band linear system solve routine.
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the dot product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dgbmv   Level 2 BLAS that computes the band matrix vector product.
c
c\Remarks
c  1. The converged Ritz values are always returned in increasing
c     (algebraic) order.
c
c  2. Currently only HOWMNY = 'A' is implemented. It is included at this
c     stage for the user who wants to incorporate it.
c
c\Author
c     Danny Sorensen
c     Richard Lehoucq
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: sband.F   SID: 2.3   DATE OF SID: 10/17/00   RELEASE: 2
c
c\EndLib
c
c---------------------------------------------------------------------
c
      subroutine dsband( rvec, howmny, select, d, z, ldz, sigma,
     &           n, lda, rfac, kl, ku, which, bmat, nev,
     &           tol, resid, ncv, v, ldv, iparam, workd, workl,
     &           lworkl, iwork, info, info_lapack, dtbsv_msg, piters)

      USE SCONTR, ONLY                :  NTERM_KLLDn, NTERM_MLLn,
     &                                   NTERM_KMSMn, NTERM_ALL
      USE SPARSE_MATRICES, ONLY       :  I_KLLDn, J_KLLDn, KLLDn,
     &                                   I_MLLn , J_MLLn , MLLn,
     &                                   I_KMSMn, J_KMSMn, KMSMn
      USE MODEL_STUF, ONLY            :  EIG_LAP_MAT_TYPE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DSBAND'

c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character        which*2, bmat, howmny
      character*1      dtbsv_msg
      character*80     pagetitle
      integer          n, lda, kl, ku, nev, ncv, ldv,
     &                 ldz, lworkl, info
      integer          info_lapack, dsaupd_loop_count
      real(double)
     &                 tol, sigma
      logical          rvec
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer          iparam(*), iwork(*)
      logical          select(*)
      real(double)  :: workd1(n), workd2(n)
      real(double)
     &                 d(*), resid(*), v(ldv,*), z(ldz,*),
     &                 rfac(lda,*), workd(*), workl(*)
! B 05/24/04 //////////////////////////////////////////////////////////B
      integer nev_user
      real(double) eig_old(n)
! E ///////////////////////////////////////////////////////////////////E

! B 06/08/08 //////////////////////////////////////////////////////////B
      integer iter, iter_old
! E ///////////////////////////////////////////////////////////////////E

! B 11/30/09 //////////////////////////////////////////////////////////B
      integer piters
! E ///////////////////////////////////////////////////////////////////E

c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer          ipntr(14)
      REAL(DOUBLE)     KLLDn_DIAG(N), MLLn_DIAG(N), KMSMn_DIAG(N)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer          ido, i, j, type, imid, itop, ibot, ierr
c
c     %------------%
c     | Parameters |
c     %------------%
c
      REAL(DOUBLE)
     &                  one, zero
      parameter        (one = 1.0, zero = 0.0)
c
c
c     %-----------------------------%
c     | LAPACK & BLAS routines used |
c     %-----------------------------%
c
!:!   REAL(DOUBLE)
      REAL(DOUBLE)
     &                 ddot, dnrm2, dlapy2
      external         ddot, dnrm2, dlapy2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      ierr = 0
! **********************************************************************************************************************************
! B 05/24/04 //////////////////////////////////////////////////////////B

      do j=1,n
         eig_old(j) = 0.d0
      enddo

      nev_user = nev

! E ///////////////////////////////////////////////////////////////////E

c     %----------------------------------------------------------------%
c     | Set type of the problem to be solved.                          |
c     |                                                                |
c     | type = 3 --> Solving generalized problem in regular mode.      |
!     |              NOTE: see expalnation of my mode 2 above          |
c     | type = 4 --> Solving generalized problem in shift-invert mode. |
c     %----------------------------------------------------------------%
c
      if      ( iparam(7) .eq. 1 ) then
         type = 1                                                            ! ARPACK "mode" = 1
      else if ( iparam(7) .eq. 2 ) then
         type = 3                                                            ! ARPACK "mode" = 2
      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'G') then
         type = 4                                                            ! ARPACK "mode" = 3
      else                                                                   ! Err msg given in BD_EIGRL
         go to 9000
      end if
c
c     %------------------------%
c     | Initialize the reverse |
c     | communication flag.    |
c     %------------------------%
c
      ido   = 0
c
c     %----------------%
c     | Exact shift is |
c     | used.          |
c     %----------------%
c
      iparam(1) = 1
c
c     %-----------------------------------%
c     | Both matrices A and M are stored  |
c     | between rows itop and ibot.  Imid |
c     | is the index of the row that      |
c     | stores the diagonal elements.     |
c     %-----------------------------------%
c
      if      (eig_lap_mat_type(1:3) == 'DGB') then
         itop = kl + 1
         ibot = 2*kl + ku + 1
         imid = kl + ku + 1
      else if (eig_lap_mat_type(1:3) == 'DPB') then
         itop = 1
         ibot = ku + 1
         imid = ibot
      endif

! Factorize matrix [K - sigma*M]

!     if (sol_name(1:8) == 'BUCKLING') then
!        modnam1 = 'LAPACK DGB TRIANG FACTOR OF [KAA+sigma*KLLDn]'
!        modnam2 = 'LAPACK DPB TRIANG FACTOR OF [KAA+sigma*KLLDn]'
!     else
!        modnam1 = 'LAPACK DGB TRIANG FACTOR OF [KAA-sigma*MLLn]'
!        modnam2 = 'LAPACK DPB TRIANG FACTOR OF [KAA-sigma*MLLn]'
!     endif
      if ( type .eq. 3 ) then
c
c        %----------------------------------------------%
c        | Solving generalized eigenvalue problem in    |
c        | regular mode. Copy MB = [K - sigma*M] to rfac|
c        | and Call LAPACK routine dpbtrf to factor rfac|
c        %----------------------------------------------%
c

         IF(SOLLIB(1:6) == 'SPARSE') THEN

            !todo not sure about MATIN_SET = 'L '
            SLU_INFO = 0
            call SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME , 'KMSM',  'L ',
     &                                   n, NTERM_KMSMn, I_KMSMn,
     &                                   J_KMSMn,  KMSMn,  SLU_INFO )

         ELSE

            if      (eig_lap_mat_type(1:3) == 'DGB') then

!               CALL LINK_MESSAGE(modnam1)
               call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)

            else if (eig_lap_mat_type(1:3) == 'DPB') then

!               CALL LINK_MESSAGE(modnam2)
               call dpbtrf ( 'U', n, ku, rfac, ku+1, ierr )
               do i=1,n
                  iwork(i) = i                                                  ! Pivot indices (no pivoting in DPBTRF)
               enddo

            endif

            if (ierr .ne. 0) then
               info_lapack = ierr
               go to 9000
            end if

         ENDIF


c
      else if ( type .eq. 4 ) then
c
c        %-------------------------------------------%
c        | Solving generalized eigenvalue problem in |
c        | shift-invert, Buckling, or Cayley mode.   |
c        %-------------------------------------------%
c
c        %-------------------------------------%
c        | Construct and factor (A - sigma*M). |
c        %-------------------------------------%
c

         IF(SOLLIB(1:6) == 'SPARSE') THEN

            !todo not sure about MATIN_SET = 'L '
            SLU_INFO = 0
            call SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME , 'KMSM',  'L ',
     &                                   n, NTERM_KMSMn, I_KMSMn,
     &                                   J_KMSMn,  KMSMn,  SLU_INFO )

         ELSE

            if      (eig_lap_mat_type(1:3) == 'DGB') then

!               CALL LINK_MESSAGE(modnam1)
               call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)

            else if (eig_lap_mat_type(1:3) == 'DPB') then

!              CALL LINK_MESSAGE(modnam2)
               call dpbtrf ( 'U', n, ku, rfac, ku+1, ierr )
               do i=1,n
                  iwork(i) = i   ! Pivot indices (no pivoting in DPBTRF)
               enddo

            endif

            if ( ierr .ne. 0 )  then
                info_lapack = ierr
                go to 9000
            end if

         ENDIF
c
      end if

      IF (EIG_MSGLVL > 0) THEN
         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            DO I=1,N
               IF (I_KLLDn(I) == I_KLLDn(I+1)) THEN
                  KLLDn_DIAG(I) = ZERO
               ELSE
                  KLLDn_DIAG(I) = KLLDn(I_KLLDn(I))
               ENDIF
               IF (I_KMSMn(I) == I_KMSMn(I+1)) THEN
                  KMSMn_DIAG(I) = ZERO
               ELSE
                  KMSMn_DIAG(I) = KMSMn(I_KMSMn(I))
               ENDIF
            ENDDO
         ELSE
            DO I=1,N
               IF (I_MLLn(I) == I_MLLn(I+1)) THEN
                  MLLn_DIAG(I) = ZERO
               ELSE
                  MLLn_DIAG(I) = MLLn(I_MLLn(I))
               ENDIF
               IF (I_KMSMn(I) == I_KMSMn(I+1)) THEN
                  KMSMn_DIAG(I) = ZERO
               ELSE
                  KMSMn_DIAG(I) = KMSMn(I_KMSMn(I))
               ENDIF
            ENDDO
         ENDIF
      ENDIF
c
c     %--------------------------------------------%
c     |  M A I N   L O O P (reverse communication) |
c     %--------------------------------------------%
c
      iter_old          = 0
      dsaupd_loop_count = 0
      write(sc1, * )
  90  continue

      IF (EIG_MSGLVL > 0) THEN
         WRITE(F06,99990)
         WRITE(F06,98710) dsaupd_loop_count, ido, type,
     &                    eig_lap_mat_type(1:3)
      ENDIF
c
! B 04/16/10 //////////////////////////////////////////////////////////B
      iter = 0
! E ///////////////////////////////////////////////////////////////////E
      call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
     &              v, ldv, iparam, ipntr, workd, workl, lworkl,
! B 05/24/04 //////////////////////////////////////////////////////////B
!    &              info )
!    &              info, nev_user, sigma, eig_old )
! E ///////////////////////////////////////////////////////////////////E
! B 06/08/08 //////////////////////////////////////////////////////////B
     &              info, nev_user, sigma, eig_old, iter, piters )
! E ///////////////////////////////////////////////////////////////////E
c
      if (iter > iter_old) then
         dsaupd_loop_count = 1
         iter_old = iter
      else
         dsaupd_loop_count = dsaupd_loop_count + 1
      endif
      write(sc1,12345,advance='no') iter+1,dsaupd_loop_count,ido,cr13_a
c     Write(f06,*) 'In ARPACK_LANCZOS_EIG: type = ', type

! **********************************************************************
      if (ido .eq. -1) then
c
         if ( type .eq. 3 ) then
c
c           %-------------------------------------%
c           | Perform  y <--- OP*x = inv[MB]*AB*x |
!           | MB = [K - sigma*M] and AB = M       |
c           | to force the starting vector into   |
c           | the range of OP.                    |
c           %-------------------------------------%
c
            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn', n, n, nterm_klldn, 'N',
     &                             i_klldn, j_klldn, klldn,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn', n, n, nterm_mlln, 'N',
     &                             i_mlln, j_mlln, mlln,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo

            IF (EIG_MSGLVL > 0) CALL ARP_DEB(1,N,IDO,IPNTR)
            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)

            IF(SOLLIB(1:6) == 'SPARSE') THEN

               SLU_INFO = 0
               call FBS_SUPRLU ( SUBR_NAME, 'KMSMn', n,
     &                        NTERM_KMSMn, I_KMSMn, J_KMSMn, KMSMn,
     &                        0, workd(ipntr(2)), SLU_INFO )

            ELSE

               if      (eig_lap_mat_type(1:3) == 'DGB') then

                  call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                          iwork, workd(ipntr(2)), n, ierr )

               else if (eig_lap_mat_type(1:3) == 'DPB') then

                  call dpbtrs ( 'U', n, ku, 1, rfac, ku+1,
     &                          workd(ipntr(2)), n, ierr )

               endif
               IF (EIG_MSGLVL > 0) CALL ARP_DEB(2,N,IDO,IPNTR)

               if (ierr .ne. 0) then
                  info_lapack = ierr
                  go to 9000
               end if

            ENDIf

c
         else if ( type .eq. 4 ) then
c
c           %-----------------------------------------%
c           | Perform y <-- OP*x                      |
c           |           = inv[A-SIGMA*M]*M            |
c           | to force the starting vector into the   |
c           | range of OP.                            |
c           %-----------------------------------------%
c
            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn' , n, n, nterm_klldn, 'N',
     &                             i_klldn, j_klldn, klldn,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn'  , n, n, nterm_mlln , 'N',
     &                             i_mlln, j_mlln, mlln,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0) CALL ARP_DEB(1,N,IDO,IPNTR)

            IF(SOLLIB(1:6) == 'SPARSE') THEN

               SLU_INFO = 0
               call FBS_SUPRLU ( SUBR_NAME, 'KMSMn', n,
     &                        NTERM_KMSMn, I_KMSMn, J_KMSMn, KMSMn,
     &                        0, workd(ipntr(2)), SLU_INFO )

            ELSE

               if      (eig_lap_mat_type(1:3) == 'DGB') then

                  call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                         iwork, workd(ipntr(2)), n, ierr )

               else if (eig_lap_mat_type(1:3) == 'DPB') then

                  call dpbtrs ( 'U', n, ku, 1, rfac, ku+1,
     &                           workd(ipntr(2)), n, ierr )
               endif
               IF (EIG_MSGLVL > 0) CALL ARP_DEB(2,N,IDO,IPNTR)

               if (ierr .ne. 0) then
                  info_lapack = ierr
                  go to 9000
               end if

            ENDIF


         endif
c
! **********************************************************************
      else if (ido .eq. 1) then
c
         if ( type .eq. 3 ) then
c
c           %----------------------------------------------%
c           | Perform  y <--- OP*x = inv[MB]*AB*x          |
!           |                      = inv[K - sigma*M]*M*x  |
c           %----------------------------------------------%
c
            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn' , n, n, nterm_klldn, 'N',
     &                             i_klldn, j_klldn, klldn,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn'  , n, n, nterm_mlln , 'N',
     &                             i_mlln , j_mlln, mlln,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0) CALL ARP_DEB(1,N,IDO,IPNTR)
            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)

            IF(SOLLIB(1:6) == 'SPARSE') THEN

               SLU_INFO = 0
               call FBS_SUPRLU ( SUBR_NAME, 'KMSMn', n,
     &                        NTERM_KMSMn, I_KMSMn, J_KMSMn, KMSMn,
     &                        0, workd(ipntr(2)), SLU_INFO )

            ELSE

               if      (eig_lap_mat_type(1:3) == 'DGB') then

                  call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                          iwork, workd(ipntr(2)), n, ierr )

               else if (eig_lap_mat_type(1:3) == 'DPB') then

                  call dpbtrs ( 'U', n, ku, 1, rfac, ku+1,
     &                          workd(ipntr(2)), n, ierr )

               endif
               IF (EIG_MSGLVL > 0) CALL ARP_DEB(2,N,IDO,IPNTR)

               if (ierr .ne. 0) then
                  info_lapack = ierr
                  go to 9000
               end if

            ENDIF
c
         else if ( type .eq. 4 ) then
c
c           %----------------------------------------%
c           | Perform y <-- inv[AB-sigma*MB]*(MB*x). |
!           |               inv[K - sigma*M]*m*x     |
c           | (MB*x) has been computed and stored    |
c           | in workd(ipntr(3)).                    |
c           %----------------------------------------%
c

            call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)

            IF(SOLLIB(1:6) == 'SPARSE') THEN

               SLU_INFO = 0
               call FBS_SUPRLU ( SUBR_NAME, 'KMSMn', n,
     &                        NTERM_KMSMn, I_KMSMn, J_KMSMn, KMSMn,
     &                        0, workd(ipntr(2)), SLU_INFO )

            ELSE

               if      (eig_lap_mat_type(1:3) == 'DGB') then

                  call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                          iwork, workd(ipntr(2)), n, ierr )

               else if (eig_lap_mat_type(1:3) == 'DPB') then

                  call dpbtrs ( 'U', n, ku, 1, rfac, ku+1,
     &                          workd(ipntr(2)), n, ierr )

               endif
               IF (EIG_MSGLVL > 0) CALL ARP_DEB(2,N,IDO,IPNTR)

               if (ierr .ne. 0) then
                  info_lapack = ierr
                  go to 9000
               end if

            ENDIF

c
         end if
c
! **********************************************************************
      else if (ido .eq. 2) then
c
c        %----------------------------------%
c        |        Perform y <-- B*x         |
c        %----------------------------------%
c
         if (type == 3) then

            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo
            call matmult_sff ( 'KMSMn', n, n, nterm_kmsmn, 'N',
     &                          i_kmsmn, j_kmsmn, kmsmn,
     &                         'workd1', n, 1, workd1, 'N',
     &                         'workd2', one,  workd2 )

!           CALL DSBMV ( 'U', N, KU, 1.0D0, RFAC, KU+1, WORKD1, 1,  ! This does not work. RFAC is not K-sigma*M anymore
!    $                   0.D0, WORKD2, 1, 1 )                       ! It was overwritten in subr dpbtrf by its triangular factor
            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0)CALL ARP_DEB(3,N,IDO,IPNTR)

         else if (type == 4) then

            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn' , n, n, nterm_klldn, 'N',
     &                             i_klldn, j_klldn, klldn,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn'  , n, n, nterm_mlln , 'N',
     &                             i_mlln , j_mlln, mlln,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0) CALL ARP_DEB(1,N,IDO,IPNTR)

         endif

c
      else
c
c        %--------------------------------------%
c        | Either we have convergence, or error |
c        %--------------------------------------%
c
         if ( info .lt. 0) then
c
            call arpack_info_msg ('dsaupd',info,iparam,lworkl,nev,ncv)
            go to 9000
c
         else
c
            if ((info == 1) .or. (info == 3)) then
               call arpack_info_msg('dsaupd',info,iparam,lworkl,nev,ncv)
            endif
c
            if (iparam(5) .gt. 0) then
c
               call dseupd ( rvec, 'A', select, d, z, ldz, sigma,
     &                  bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &                  iparam, ipntr, workd, workl, lworkl, info )
c
               if ( info .ne. 0) then
c
              call arpack_info_msg ('dseupd',info,iparam,lworkl,nev,ncv)
                  go to 9000
c
               end if
c
            end if
c
         end if
c
         go to 9000
c
      end if
c
c     %----------------------------------------%
c     | L O O P  B A C K to call DSAUPD again. |
c     %----------------------------------------%
c
      go to 90
c
 9000 continue
c
! B 05/24/04 //////////////////////////////////////////////////////////B
      Write(f06,*) ! blank line after eigen iteration results
! E ///////////////////////////////////////////////////////////////////E
! **********************************************************************************************************************************
12345 format(5X,'Iteration',i4,' Rev comm loop',i4,' with IDO =',i3,a)

 4907 FORMAT(/,22X,A
     &      ,/,7X,'1',12X,'2',12X,'3',12X,'4',12X,'5',12X,'6',12X,
     &        '7',12X,'8',12X,'9',12X,'10')

 4908 FORMAT(10(1X,1ES12.5))

98710 FORMAT(' dsaupd loop count = ',I4,' ido = ',i4,', "type" = ',I3,
     &       ', using ',a,' LAPACK matrices',/)

99990 FORMAT('**********************************************************
     &***************************')


      RETURN

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE ARP_DEB ( WHICH, N, IDO, IPNTR )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  PROG_NAME, FATAL_ERR, WARN_ERR
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, ERR, F06

      IMPLICIT NONE

      INTEGER, INTENT(IN)            :: IDO                !
      INTEGER, INTENT(IN)            :: IPNTR(14)          !
      INTEGER, INTENT(IN)            :: N                  !
      INTEGER, INTENT(IN)            :: WHICH              !
      INTEGER                        :: I,J                ! DO loop

! **********************************************************************************************************************************
      IF (WHICH == 1) THEN

         WRITE(F06,98700) IDO, IPNTR(1), IPNTR(2)
         IF (EIG_MSGLVL > 1) THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,98711)
               DO J=1,N
                  WRITE(F06,98713) J, KLLDn_DIAG(J), WORKD1(J),WORKD2(J)
               ENDDO
            ELSE
               WRITE(F06,98712)
               DO J=1,N
                  WRITE(F06,98713) J, MLLn_DIAG(J), WORKD1(J), WORKD2(J)
               ENDDO
            ENDIF
            WRITE(F06,*)
         ENDIF

      ELSE IF (WHICH == 2) THEN

         WRITE(F06,98700) IDO, IPNTR(1), IPNTR(2)
         IF (EIG_MSGLVL > 1) THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,98721)
            ELSE
               WRITE(F06,98722)
            ENDIF
            DO J=1,N
               WRITE(F06,98723) J, RFAC(IMID,J), WORKD1(J), WORKD2(J)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ELSE IF (WHICH == 3) THEN

         WRITE(F06,98700) IDO, IPNTR(1), IPNTR(2)
         IF (EIG_MSGLVL > 1) THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,98731)
            ELSE
               WRITE(F06,98732)
            ENDIF
            DO J=1,N
               WRITE(F06,98733) J, KMSMn_DIAG(J), WORKD1(J), WORKD2(J)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ENDIF

! **********************************************************************************************************************************
98700 FORMAT(' ido, ipntr(1), ipntr(2) = ',I9,3x,2i3,/)

98711 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF                KLL                  WORKD1(J)      WOR
     &KD2(J)',/,' ---------        -------------------        ----------
     &--   ------------')

98712 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF                MLL                  WORKD1(J)      WOR
     &KD2(J)',/,' ---------        -------------------        ----------
     &--   ------------')

98713 FORMAT(I8,12X,1ES14.6,9X,1ES14.6,1X,1ES14.6)

98721 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF  tridiag factor of [KLL-sig*MLL]    WORKD1(J)      WOR
     &KD2(J)',/,' ---------  -------------------------------  ----------
     &--   ------------')

98722 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF  tridiag factor of [KLL-sig*KLL]    WORKD1(J)      WOR
     &KD2(J)',/,' ---------  -------------------------------  ----------
     &--   ------------')

98723 FORMAT(I8,12X,1ES14.6,9X,1ES14.6,1X,1ES14.6)

98731 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF          [KLL-sigma*MLL]            WORKD1(J)      WOR
     &KD2(J)',/,' ---------        -------------------        ----------
     &--   ------------')

98732 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF          [KLL-sigma*KLL]            WORKD1(J)      WOR
     &KD2(J)',/,' ---------        -------------------        ----------
     &--   ------------')

98733 FORMAT(I8,12X,1ES14.6,9X,1ES14.6,1X,1ES14.6)

! **********************************************************************************************************************************

      END SUBROUTINE ARP_DEB

      end subroutine dsband

! ##################################################################################################################################
! 002 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsaupd
c
c\Description:
c
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  Iteration.  For symmetric problems this reduces to a variant of the Lanczos
c  method.  This method has been designed to compute approximations to a
c  few eigenpairs of a linear operator OP that is real and symmetric
c  with respect to a real positive semi-definite symmetric matrix B,
c  i.e.
c
c       B*OP = (OP`)*B.
c
c  Another way to express this condition is
c
c       < x,OPy > = < OPx,y >  where < z,w > = z`Bw  .
c
c  In the standard eigenproblem B is the identity matrix.
c  ( A` denotes transpose of A)
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  dsaupd  is usually called iteratively to solve one of the
c  following problems:
c
c  Mode 1:  A*x = lambda*x, A symmetric
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
c           ===> OP = (inv[K - sigma*M])*M  and  B = M.
c           ===> Shift-and-Invert mode
c
c  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite,
c           KG symmetric indefinite
c           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
c           ===> Buckling mode
c
c  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
c           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
c           ===> Cayley transformed mode
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call dsaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO, PITERS )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first
c          call to dsaupd .  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          dsaupd  with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          (If Mode = 2 see remark 5 below)
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3,4 and 5, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute the IPARAM(8) shifts where
c                    IPNTR(11) is the pointer into WORKL for
c                    placing the shifts. See remark 6 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          Specify which of the Ritz values of OP to compute.
c
c          'LA' - compute the NEV largest (algebraic) eigenvalues.
c          'SA' - compute the NEV smallest (algebraic) eigenvalues.
c          'LM' - compute the NEV largest (in magnitude) eigenvalues.
c          'SM' - compute the NEV smallest (in magnitude) eigenvalues.
c          'BE' - compute NEV eigenvalues, half from each end of the
c                 spectrum.  When NEV is odd, compute one more from the
c                 high end than from the low end.
c           (see remark 1 below)
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N.
c
c  TOL     REAL(DOUBLE)  scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = DLAMCH ('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH ).
c
c  RESID   REAL(DOUBLE)  array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          This will indicate how many Lanczos vectors are generated
c          at each iteration.  After the startup phase in which NEV
c          Lanczos vectors are generated, the algorithm generates
c          NCV-NEV Lanczos vectors at each subsequent update iteration.
c          Most of the cost in generating each Lanczos vector is in the
c          matrix-vector product OP*x. (See remark 4 below).
c
c  V       REAL(DOUBLE)  N by NCV array.  (OUTPUT)
c          The NCV columns of V contain the Lanczos basis vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The NCV eigenvalues of
c                      the current tridiagonal matrix T are returned in
c                      the part of WORKL array corresponding to RITZ.
c                      See remark 6 below.
c          ISHIFT = 1: exact shifts with respect to the reduced
c                      tridiagonal matrix T.  This is equivalent to
c                      restarting the iteration with a starting vector
c                      that is a linear combination of Ritz vectors
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = LEVEC
c          No longer referenced. See remark 2 below.
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed.
c          On OUTPUT: actual number of Arnoldi update iterations taken.
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used.
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of dsaupd  for the
c          five modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), dsaupd  returns NP, the number
c          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
c          6 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c  IPNTR   Integer array of length 11.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Lanczos iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
c          IPNTR(6): pointer to the NCV RITZ values array in WORKL.
c          IPNTR(7): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZ in WORKL.
c          IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
c
c          Note: IPNTR(8:10) is only referenced by dseupd . See Remark 2.
c          IPNTR(8): pointer to the NCV RITZ values of the original system.
c          IPNTR(9): pointer to the NCV corresponding error bounds.
c          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
c                     of the tridiagonal matrix T. Only referenced by
c                     dseupd  if RVEC = .TRUE. See Remarks.
c          -------------------------------------------------------------
c
c  WORKD   REAL(DOUBLE)  work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If the Ritz vectors are desired
c          subroutine dseupd  uses this output.
c          See Data Distribution Note below.
c
c  WORKL   REAL(DOUBLE)  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least NCV**2 + 8*NCV .
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the
c                Implicitly restarted Arnoldi iteration. One possibility
c                is to increase the size of NCV relative to NEV.
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iterations allowed
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array WORKL is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informatinal error from LAPACK routine dsteqr .
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -13: NEV and WHICH = 'BE' are incompatable.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization. The user is advised to check that
c                   enough workspace and array storage has been allocated.
c
c
c\Remarks
c  1. The converged Ritz values are always returned in ascending
c     algebraic order.  The computed Ritz values are approximate
c     eigenvalues of OP.  The selection of WHICH should be made
c     with this in mind when Mode = 3,4,5.  After convergence,
c     approximate eigenvalues of the original problem may be obtained
c     with the ARPACK subroutine dseupd .
c
c  2. If the Ritz vectors corresponding to the converged Ritz values
c     are needed, the user must call dseupd  immediately following completion
c     of dsaupd . This is new starting with version 2.1 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL`
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular
c     linear systems should be solved with L and L` rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L`z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requrement is that NCV > NEV.
c     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.   The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically.
c
c  5. If IPARAM(7) = 2 then in the Reverse commuication interface the user
c     must do the following. When IDO = 1, Y = OP * X is to be computed.
c     When IPARAM(7) = 2 OP = inv(B)*A. After computing A*X the user
c     must overwrite X with A*X. Y is then the solution to the linear set
c     of equations B*Y = A*X.
c
c  6. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
c     NP = IPARAM(8) shifts in locations:
c     1   WORKL(IPNTR(11))
c     2   WORKL(IPNTR(11)+1)
c                        .
c                        .
c                        .
c     NP  WORKL(IPNTR(11)+NP-1).
c
c     The eigenvalues of the current tridiagonal matrix are located in
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are in the
c     order defined by WHICH. The associated Ritz estimates are located in
c     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note:
c
c  Fortran-D syntax:
c  ================
c  REAL       RESID(N), V(LDV,NCV), WORKD(3*N), WORKL(LWORKL)
c  DECOMPOSE  D1(N), D2(N,NCV)
c  ALIGN      RESID(I) with D1(I)
c  ALIGN      V(I,J)   with D2(I,J)
c  ALIGN      WORKD(I) with D1(I)     range (1:N)
c  ALIGN      WORKD(I) with D1(I-N)   range (N+1:2*N)
c  ALIGN      WORKD(I) with D1(I-2*N) range (2*N+1:3*N)
c  DISTRIBUTE D1(BLOCK), D2(BLOCK,:)
c  REPLICATED WORKL(LWORKL)
c
c  Cray MPP syntax:
c  ===============
c  REAL       RESID(N), V(LDV,NCV), WORKD(N,3), WORKL(LWORKL)
c  SHARED     RESID(BLOCK), V(BLOCK,:), WORKD(BLOCK,:)
c  REPLICATED WORKL(LWORKL)
c
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
c     1980.
c  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
c     Computer Physics Communications, 53 (1989), pp 169-179.
c  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
c     Implement the Spectral Transformation", Math. Comp., 48 (1987),
c     pp 663-673.
c  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos
c     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems",
c     SIAM J. Matr. Anal. Apps.,  January (1993).
c  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
c     for Updating the QR decomposition", ACM TOMS, December 1990,
c     Volume 16 Number 4, pp 369-377.
c  8. R.B. Lehoucq, D.C. Sorensen, "Implementation of Some Spectral
c     Transformations in a k-Step Arnoldi Method". In Preparation.
c
c\Routines called:
c     dsaup2   ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     dstats   ARPACK routine that initialize timing and other statistics
c             variables.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dvout    ARPACK utility routine that prints vectors.
c     dlamch   LAPACK routine that determines machine constants.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     12/15/93: Version ' 2.4'
c
c\SCCS Information: @(#)
c FILE: saupd.F   SID: 2.8   DATE OF SID: 04/10/01   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
! B 05/24/04 //////////////////////////////////////////////////////////B
!    &     ipntr, workd, workl, lworkl, info )
!    &     ipntr, workd, workl, lworkl, info, nev_user, sigma, eig_old )
! E ///////////////////////////////////////////////////////////////////E
! B 06/08/08 //////////////////////////////////////////////////////////B
     &     ipntr, workd, workl, lworkl, info, nev_user, sigma, eig_old,
     &     iter, piters )
! E ///////////////////////////////////////////////////////////////////E


! B 06/08/08 //////////////////////////////////////////////////////////B
      integer iter
! E ///////////////////////////////////////////////////////////////////E

! B 05/24/04 //////////////////////////////////////////////////////////B
      integer n, nev_user
      real(double) sigma, eig_old(n)
! E ///////////////////////////////////////////////////////////////////E

! B 11/30/09 //////////////////////////////////////////////////////////B
      integer piters
! E ///////////////////////////////////////////////////////////////////E

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSAUPD'
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cxxxx include   'debug.h'
cxxxx include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, ncv, nev
      REAL(DOUBLE)
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(11)
      REAL(DOUBLE)
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      REAL(DOUBLE)
     &           one, zero
      parameter (one = 1.0D+0 , zero = 0.0D+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, msglvl, mxiter, mode, nb,
     &           nev0, next, np, ritz, j
      INTEGER
     &           MXITER_ARRAY(1), NP_ARRAY(1)
      save       bounds, ierr, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, msglvl, mxiter, mode, nb,
     &           nev0, next, np, ritz
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
!     external   dsaup2 ,  dvout , ivout, second, dstats
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      REAL(DOUBLE)
     &           dlamch
      external   dlamch
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c


! **********************************************************************************************************************************

      if (ido .eq. 0) then
c
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call dstats
         call cpu_time (t0)
         msglvl = msaupd
c
         ierr   = 0
         ishift = iparam(1)
         mxiter = iparam(3)
c         nb     = iparam(4)
         nb     = 1
c
c        %--------------------------------------------%
c        | Revision 2 performs only implicit restart. |
c        %--------------------------------------------%
c
         iupd   = 1
         mode   = iparam(7)
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         if (n .le. 0) then
            ierr = -1
         else if (nev .le. 0) then
            ierr = -2
         else if (ncv .le. nev .or.  ncv .gt. n) then
            ierr = -3
         end if
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        %----------------------------------------------%
c
         np     = ncv - nev
c
         if (mxiter .le. 0)                     ierr = -4
         if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LA' .and.
     &       which .ne. 'SA' .and.
     &       which .ne. 'BE')                   ierr = -5
         if (bmat .ne. 'I' .and. bmat .ne. 'G') ierr = -6
c
         if (lworkl .lt. ncv**2 + 8*ncv)        ierr = -7
         if (mode .lt. 1 .or. mode .gt. 5) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
                                                ierr = -12
         else if (nev .eq. 1 .and. which .eq. 'BE') then
                                                ierr = -13
         end if
c
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)                         nb = 1
         if (tol .le. zero)                     tol = dlamch ('EpsMach')
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        | NEV0 is the local variable designating the   |
c        | size of the invariant subspace desired.      |
c        %----------------------------------------------%
c
         np     = ncv - nev
         nev0   = nev
c
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, ncv**2 + 8*ncv
            workl(j) = zero
 10      continue
c
c        %-------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
c        | etc... and the remaining workspace.                   |
c        | Also update pointer to be used on output.             |
c        | Memory is laid out as follows:                        |
c        | workl(1:2*ncv) := generated tridiagonal matrix        |
c        | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
c        | workl(3*ncv+1:3*ncv+ncv) := computed error bounds     |
c        | workl(4*ncv+1:4*ncv+ncv*ncv) := rotation matrix Q     |
c        | workl(4*ncv+ncv*ncv+1:7*ncv+ncv*ncv) := workspace     |
c        %-------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritz   = ih     + 2*ldh
         bounds = ritz   + ncv
         iq     = bounds + ncv
         iw     = iq     + ncv**2
         next   = iw     + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritz
         ipntr(7) = bounds
         ipntr(11) = iw
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Lanczos Iteration. |
c     %-------------------------------------------------------%
c
      call dsaup2
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritz),
     &     workl(bounds), workl(iq), ldq, workl(iw), ipntr, workd,
! B 05/24/04 //////////////////////////////////////////////////////////B
!    &     info )
!    &     info, nev_user, sigma, eig_old )
! B 06/08/08 //////////////////////////////////////////////////////////B
     &     info, nev_user, sigma, eig_old, iter, piters )

! E ///////////////////////////////////////////////////////////////////E
c
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP or shifts.    |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
c
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
c
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within dsaup2 .               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
      if (msglvl .gt. 0) then
         MXITER_ARRAY(1) = MXITER
         NP_ARRAY(1)     = NP
         call ivout (logfil, 1, MXITER_ARRAY, ndigit,
     &               '_saupd: number of update iterations taken')
         call ivout (logfil, 1, NP_ARRAY, ndigit,
     &               '_saupd: number of "converged" Ritz values')
         call dvout  (logfil, np, workl(Ritz), ndigit,
     &               '_saupd: final Ritz values')
         call dvout  (logfil, np, workl(Bounds), ndigit,
     &               '_saupd: corresponding error bounds')
      end if
c
      call cpu_time (t1)
      tsaupd = t1 - t0
c
      if (msglvl .gt. 0) then
c
c        %--------------------------------------------------------%
c        | Version Number & Version Date are defined in version.h |
c        %--------------------------------------------------------%
c
         write (logfil,1000)
         write (logfil,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt,
     &                  tmvopx, tmvbx, tsaupd, tsaup2, tsaitr, titref,
     &                  tgetv0, tseigt, tsgets, tsapps, tsconv
 1000    format (//,
     &      5x, '==========================================',/
     &      5x, '= Symmetric implicit Arnoldi update code =',/
     &      5x, '= Version Number:', ' 2.4' , 19x, ' =',/
     &      5x, '= Version Date:  ', ' 07/31/96' , 14x, ' =',/
     &      5x, '==========================================',/
     &      5x, '= Summary of timing statistics           =',/
     &      5x, '==========================================',//)
 1100    format (
     &      5x, 'Total number update iterations             = ', i5,/
     &      5x, 'Total number of OP*x operations            = ', i5,/
     &      5x, 'Total number of B*x operations             = ', i5,/
     &      5x, 'Total number of reorthogonalization steps  = ', i5,/
     &      5x, 'Total number of iterative refinement steps = ', i5,/
     &      5x, 'Total number of restart steps              = ', i5,/
     &      5x, 'Total time in user OP*x operation          = ', f12.6,/
     &      5x, 'Total time in user B*x operation           = ', f12.6,/
     &      5x, 'Total time in Arnoldi update routine       = ', f12.6,/
     &      5x, 'Total time in saup2 routine                = ', f12.6,/
     &      5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/
     &      5x, 'Total time in reorthogonalization phase    = ', f12.6,/
     &      5x, 'Total time in (re)start vector generation  = ', f12.6,/
     &      5x, 'Total time in trid eigenvalue subproblem   = ', f12.6,/
     &      5x, 'Total time in getting the shifts           = ', f12.6,/
     &      5x, 'Total time in applying the shifts          = ', f12.6,/
     &      5x, 'Total time in convergence testing          = ', f12.6)
      end if
c
 9000 continue
c

      RETURN

! **********************************************************************************************************************************

c
c     %---------------%
c     | End of dsaupd |
c     %---------------%
c
      end subroutine dsaupd

! ##################################################################################################################################
! 003 ARPACK_LANCZOS_EIG

c\BeginDoc
c
c\Name: dseupd
c
c\Description:
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) the corresponding approximate eigenvectors,
c
c      (2) an orthonormal (Lanczos) basis for the associated approximate
c          invariant subspace,
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  (Lanczos) basis is always computed.  There is an additional storage cost
c  of n*nev if both are requested (in this case a separate array Z must be
c  supplied).
c
c  These quantities are obtained from the Lanczos factorization computed
c  by DSAUPD  for the linear operator OP prescribed by the MODE selection
c  (see IPARAM(7) in DSAUPD  documentation.)  DSAUPD  must be called before
c  this routine is called. These approximate eigenvalues and vectors are
c  commonly called Ritz values and Ritz vectors respectively.  They are
c  referred to as such in the comments that follow.   The computed orthonormal
c  basis for the invariant subspace corresponding to these Ritz values is
c  referred to as a Lanczos basis.
c
c  See documentation in the header of the subroutine DSAUPD  for a definition
c  of OP as well as other terms and the relation of computed Ritz values
c  and vectors of OP with respect to the given problem  A*z = lambda*B*z.
c
c  The approximate eigenvalues of the original problem are returned in
c  ascending algebraic order.  The user may elect to call this routine
c  once for each desired Ritz vector and store it peripherally if desired.
c  There is also the option of computing a selected set of these vectors
c  with a single call.
c
c\Usage:
c  call dseupd
c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
c       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c  RVEC    LOGICAL  (INPUT)
c          Specifies whether Ritz vectors corresponding to the Ritz value
c          approximations to the eigenproblem A*z = lambda*B*z are computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute Ritz vectors.
c
c  HOWMNY  Character*1  (INPUT)
c          Specifies how many Ritz vectors are wanted and the form of Z
c          the matrix of Ritz vectors. See remark 1 below.
c          = 'A': compute NEV Ritz vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT/WORKSPACE)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE..
c          If HOWMNY = 'A' , SELECT is used as a workspace for
c          reordering the Ritz values.
c
c  D       REAL(DOUBLE)  array of dimension NEV.  (OUTPUT)
c          On exit, D contains the Ritz value approximations to the
c          eigenvalues of A*z = lambda*B*z. The values are returned
c          in ascending order. If IPARAM(7) = 3,4,5 then D represents
c          the Ritz values of OP computed by dsaupd  transformed to
c          those of the original eigensystem A*z = lambda*B*z. If
c          IPARAM(7) = 1,2 then the Ritz values of OP are the same
c          as the those of A*z = lambda*B*z.
c
c  Z       REAL(DOUBLE)  N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c          On exit, Z contains the B-orthonormal Ritz vectors of the
c          eigensystem A*z = lambda*B*z corresponding to the Ritz
c          value approximations.
c          If  RVEC = .FALSE. then Z is not referenced.
c          NOTE: The array Z may be set equal to first NEV columns of the
c          Arnoldi/Lanczos basis array V computed by DSAUPD .
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
c
c  SIGMA   REAL(DOUBLE)   (INPUT)
c          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c          IPARAM(7) = 1 or 2.
c
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to DSAUPD  that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to DSEUPD  following the last call
c         to DSAUPD .  These arguments MUST NOT BE MODIFIED between
c         the the last call to DSAUPD  and the call to DSEUPD .
c
c  Two of these parameters (WORKL, INFO) are also output parameters:
c
c  WORKL   REAL(DOUBLE)  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:4*ncv) contains information obtained in
c          dsaupd .  They are not changed by dseupd .
c          WORKL(4*ncv+1:ncv*ncv+8*ncv) holds the
c          untransformed Ritz values, the computed error estimates,
c          and the associated eigenvector matrix of H.
c
c          Note: IPNTR(8:10) contains the pointer into WORKL for addresses
c          of the above information computed by dseupd .
c          -------------------------------------------------------------
c          IPNTR(8): pointer to the NCV RITZ values of the original system.
c          IPNTR(9): pointer to the NCV corresponding error bounds.
c          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
c                     of the tridiagonal matrix T. Only referenced by
c                     dseupd  if RVEC = .TRUE. See Remarks.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Information error from LAPACK routine dsteqr .
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: NEV and WHICH = 'BE' are incompatible.
c          = -14: DSAUPD  did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: HOWMNY must be one of 'A' or 'S' if RVEC = .true.
c          = -16: HOWMNY = 'S' not yet implemented
c          = -17: DSEUPD  got a different count of the number of converged
c                 Ritz values than DSAUPD  got.  This indicates the user
c                 probably made an error in passing data from DSAUPD  to
c                 DSEUPD  or that the data was modified before entering
c                 DSEUPD .
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
c     1980.
c  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
c     Computer Physics Communications, 53 (1989), pp 169-179.
c  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
c     Implement the Spectral Transformation", Math. Comp., 48 (1987),
c     pp 663-673.
c  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos
c     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems",
c     SIAM J. Matr. Anal. Apps.,  January (1993).
c  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
c     for Updating the QR decomposition", ACM TOMS, December 1990,
c     Volume 16 Number 4, pp 369-377.
c
c\Remarks
c  1. The converged Ritz values are always returned in increasing
c     (algebraic) order.
c
c  2. Currently only HOWMNY = 'A' is implemented. It is included at this
c     stage for the user who wants to incorporate it.
c
c\Routines called:
c     dsesrt   ARPACK routine that sorts an array X, and applies the
c             corresponding permutation to a matrix A.
c     dsortr   dsortr   ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     dvout    ARPACK utility routine that prints vectors.
c     dgeqr2   LAPACK routine that computes the QR factorization of
c             a matrix.
c     dlacpy   LAPACK matrix copy routine.
c     dlamch   LAPACK routine that determines machine constants.
c     dorm2r   LAPACK routine that applies an orthogonal matrix in
c             factored form.
c     dsteqr   LAPACK routine that computes eigenvalues and eigenvectors
c             of a tridiagonal matrix.
c     dger     Level 2 BLAS rank one update to a matrix.
c     dcopy    Level 1 BLAS that copies one vector to another .
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c     dscal    Level 1 BLAS that scales a vector.
c     dswap    Level 1 BLAS that swaps the contents of two vectors.

c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Chao Yang                    Houston, Texas
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     12/15/93: Version ' 2.1'
c
c\SCCS Information: @(#)
c FILE: seupd.F   SID: 2.11   DATE OF SID: 04/10/01   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine dseupd (rvec  , howmny, select, d    ,
     &                   z     , ldz   , sigma , bmat ,
     &                   n     , which , nev   , tol  ,
     &                   resid , ncv   , v     , ldv  ,
     &                   iparam, ipntr , workd , workl,
     &                   lworkl, info )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSEUPD'
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cxxxx include   'debug.h'
cxxxx include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      REAL(DOUBLE)
     &           sigma, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      REAL(DOUBLE)
     &           d(nev)     , resid(n)  , v(ldv,ncv),
     &           z(ldz, nev), workd(2*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      REAL(DOUBLE)
     &           one, zero
      parameter (one = 1.0D+0 , zero = 0.0D+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      integer    bounds , ierr   , ih    , ihb   , ihd   ,
     &           iq     , iw     , j     , k     , ldh   ,
     &           ldq    , mode   , msglvl, nconv , next  ,
     &           ritz   , irz    , ibd   , np    , ishift,
     &           leftptr, rghtptr, numcnv, jj
      INTEGER
     &           NUMCNV_ARRAY(1), NCONV_ARRAY(1)
      REAL(DOUBLE)
     &           bnorm2 , rnorm, temp, temp1, eps23
      REAL(DOUBLE) TEMP_ARRAY(1)
      logical    reord
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
!:!   external   dcopy  , dger   , dgeqr2 , dlacpy , dorm2r , dscal ,
!:!  &           dsesrt , dsteqr , dswap  , dvout  , ivout , dsortr
      external   dswap  , dscal  , dsteqr
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
!:!   REAL(DOUBLE)
      REAL(DOUBLE)
     &           dnrm2
      external   dnrm2

      REAL(DOUBLE)
     &           dlamch
      external   dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%


! **********************************************************************************************************************************

      TEMP = 0.0D0
c
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      msglvl = mseupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
c
c     %--------------%
c     | Quick return |
c     %--------------%
c
      if (nconv .eq. 0) go to 9000
      ierr = 0
c
      if (nconv .le. 0)                        ierr = -14
      if (n .le. 0)                            ierr = -1
      if (nev .le. 0)                          ierr = -2
      if (ncv .le. nev .or.  ncv .gt. n)       ierr = -3
      if (which .ne. 'LM' .and.
     &    which .ne. 'SM' .and.
     &    which .ne. 'LA' .and.
     &    which .ne. 'SA' .and.
     &    which .ne. 'BE')                     ierr = -5
      if (bmat .ne. 'I' .and. bmat .ne. 'G')   ierr = -6
      if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec )
     &                                         ierr = -15
      if (rvec .and. howmny .eq. 'S')           ierr = -16
c
      if (rvec .and. lworkl .lt. ncv**2+8*ncv) ierr = -7
c
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 ) then
         type = 'SHIFTI'
      else if (mode .eq. 4 ) then
         type = 'BUCKLE'
      else if (mode .eq. 5 ) then
         type = 'CAYLEY'
      else
                                               ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')     ierr = -11
      if (nev .eq. 1 .and. which .eq. 'BE')    ierr = -12
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
c
c     %-------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
c     | etc... and the remaining workspace.                   |
c     | Also update pointer to be used on output.             |
c     | Memory is laid out as follows:                        |
c     | workl(1:2*ncv) := generated tridiagonal matrix H      |
c     |       The subdiagonal is stored in workl(2:ncv).      |
c     |       The dead spot is workl(1) but upon exiting      |
c     |       dsaupd  stores the B-norm of the last residual   |
c     | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
c     |       The wanted values are in the first NCONV spots. |
c     | workl(3*ncv+1:3*ncv+ncv) := computed Ritz estimates   |
c     |       The wanted values are in the first NCONV spots. |
c     | NOTE: workl(1:4*ncv) is set by dsaupd  and is not      |
c     |       modified by dseupd .                             |
c     %-------------------------------------------------------%
c
c     %-------------------------------------------------------%
c     | The following is used and set by dseupd .              |
c     | workl(4*ncv+1:4*ncv+ncv) := used as workspace during  |
c     |       computation of the eigenvectors of H. Stores    |
c     |       the diagonal of H. Upon EXIT contains the NCV   |
c     |       Ritz values of the original system. The first   |
c     |       NCONV spots have the wanted values. If MODE =   |
c     |       1 or 2 then will equal workl(2*ncv+1:3*ncv).    |
c     | workl(5*ncv+1:5*ncv+ncv) := used as workspace during  |
c     |       computation of the eigenvectors of H. Stores    |
c     |       the subdiagonal of H. Upon EXIT contains the    |
c     |       NCV corresponding Ritz estimates of the         |
c     |       original system. The first NCONV spots have the |
c     |       wanted values. If MODE = 1,2 then will equal    |
c     |       workl(3*ncv+1:4*ncv).                           |
c     | workl(6*ncv+1:6*ncv+ncv*ncv) := orthogonal Q that is  |
c     |       the eigenvector matrix for H as returned by     |
c     |       dsteqr . Not referenced if RVEC = .False.        |
c     |       Ordering follows that of workl(4*ncv+1:5*ncv)   |
c     | workl(6*ncv+ncv*ncv+1:6*ncv+ncv*ncv+2*ncv) :=         |
c     |       Workspace. Needed by dsteqr  and by dseupd .      |
c     | GRAND total of NCV*(NCV+8) locations.                 |
c     %-------------------------------------------------------%
c
c
      ih     = ipntr(5)
      ritz   = ipntr(6)
      bounds = ipntr(7)
      ldh    = ncv
      ldq    = ncv
      ihd    = bounds + ldh
      ihb    = ihd    + ldh
      iq     = ihb    + ldh
      iw     = iq     + ldh*ncv
      next   = iw     + 2*ncv
      ipntr(4)  = next
      ipntr(8)  = ihd
      ipntr(9)  = ihb
      ipntr(10) = iq
c
c     %----------------------------------------%
c     | irz points to the Ritz values computed |
c     |     by _seigt before exiting _saup2.   |
c     | ibd points to the Ritz estimates       |
c     |     computed by _seigt before exiting  |
c     |     _saup2.                            |
c     %----------------------------------------%
c
      irz = ipntr(11)+ncv
      ibd = irz+ncv
c
c
c     %---------------------------------%
c     | Set machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch ('Epsilon-Machine')
      eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
c     %---------------------------------------%
c     | RNORM is B-norm of the RESID(1:N).    |
c     | BNORM2 is the 2 norm of B*RESID(1:N). |
c     | Upon exit of dsaupd  WORKD(1:N) has    |
c     | B*RESID(1:N).                         |
c     %---------------------------------------%
c
      rnorm = workl(ih)
      if (bmat .eq. 'I') then
         bnorm2 = rnorm
      else if (bmat .eq. 'G') then
         bnorm2 = dnrm2 (n, workd, 1)
      end if
c
      if (msglvl .gt. 2) then
         call dvout (logfil, ncv, workl(irz), ndigit,
     &   '_seupd: Ritz values passed in from _SAUPD.')
         call dvout (logfil, ncv, workl(ibd), ndigit,
     &   '_seupd: Ritz estimates passed in from _SAUPD.')
      end if
c
      if (rvec) then
c
         reord = .false.
c
c        %---------------------------------------------------%
c        | Use the temporary bounds array to store indices   |
c        | These will be used to mark the select array later |
c        %---------------------------------------------------%
c
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
c
c        %-------------------------------------%
c        | Select the wanted Ritz values.      |
c        | Sort the Ritz values so that the    |
c        | wanted ones appear at the tailing   |
c        | NEV positions of workl(irr) and     |
c        | workl(iri).  Move the corresponding |
c        | error estimates in workl(bound)     |
c        | accordingly.                        |
c        %-------------------------------------%
c
         np     = ncv - nev
         ishift = 0
         call dsgets (ishift, which       , nev          ,
     &                np    , workl(irz)  , workl(bounds),
     &                workl)
c
         if (msglvl .gt. 2) then
            call dvout (logfil, ncv, workl(irz), ndigit,
     &      '_seupd: Ritz values after calling _SGETS.')
            call dvout (logfil, ncv, workl(bounds), ndigit,
     &      '_seupd: Ritz value indices after calling _SGETS.')
         end if
c
c        %-----------------------------------------------------%
c        | Record indices of the converged wanted Ritz values  |
c        | Mark the select array for possible reordering       |
c        %-----------------------------------------------------%
c
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23, abs(workl(irz+ncv-j)) )
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          workl(ibd+jj-1) .le. tol*temp1) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
c
c        %-----------------------------------------------------------%
c        | Check the count (numcnv) of converged Ritz values with    |
c        | the number (nconv) reported by _saupd.  If these two      |
c        | are different then there has probably been an error       |
c        | caused by incorrect passing of the _saupd data.           |
c        %-----------------------------------------------------------%
c
         NUMCNV_ARRAY(1) = NUMCNV
         NCONV_ARRAY(1)  = NCONV
         if (msglvl .gt. 2) then
             call ivout(logfil, 1, NUMCNV_ARRAY, ndigit,
     &            '_seupd: Number of specified eigenvalues')
             call ivout(logfil, 1, NCONV_ARRAY, ndigit,
     &            '_seupd: Number of "converged" eigenvalues')
         end if
c
         if (numcnv .ne. nconv) then
            info = -17
            go to 9000
         end if
c
c        %-----------------------------------------------------------%
c        | Call LAPACK routine _steqr to compute the eigenvalues and |
c        | eigenvectors of the final symmetric tridiagonal matrix H. |
c        | Initialize the eigenvector matrix Q to the identity.      |
c        %-----------------------------------------------------------%
c
         call dcopy (ncv-1, workl(ih+1), 1, workl(ihb), 1)
         call dcopy (ncv, workl(ih+ldh), 1, workl(ihd), 1)
c
         call dsteqr ('Identity', ncv, workl(ihd), workl(ihb),
     &                workl(iq) , ldq, workl(iw), ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c
         if (msglvl .gt. 1) then
            call dcopy (ncv, workl(iq+ncv-1), ldq, workl(iw), 1)
            call dvout (logfil, ncv, workl(ihd), ndigit,
     &          '_seupd: NCV Ritz values of the final H matrix')
            call dvout (logfil, ncv, workl(iw), ndigit,
     &           '_seupd: last row of the eigenvector matrix for H')
         end if
c
         if (reord) then
c
c           %---------------------------------------------%
c           | Reordered the eigenvalues and eigenvectors  |
c           | computed by _steqr so that the "converged"  |
c           | eigenvalues appear in the first NCONV       |
c           | positions of workl(ihd), and the associated |
c           | eigenvectors appear in the first NCONV      |
c           | columns.                                    |
c           %---------------------------------------------%
c
            leftptr = 1
            rghtptr = ncv
c
            if (ncv .eq. 1) go to 30
c
 20         if (select(leftptr)) then
c
c              %-------------------------------------------%
c              | Search, from the left, for the first Ritz |
c              | value that has not converged.             |
c              %-------------------------------------------%
c
               leftptr = leftptr + 1
c
            else if ( .not. select(rghtptr)) then
c
c              %----------------------------------------------%
c              | Search, from the right, the first Ritz value |
c              | that has converged.                          |
c              %----------------------------------------------%
c
               rghtptr = rghtptr - 1
c
            else
c
c              %----------------------------------------------%
c              | Swap the Ritz value on the left that has not |
c              | converged with the Ritz value on the right   |
c              | that has converged.  Swap the associated     |
c              | eigenvector of the tridiagonal matrix H as   |
c              | well.                                        |
c              %----------------------------------------------%
c
               temp = workl(ihd+leftptr-1)
               workl(ihd+leftptr-1) = workl(ihd+rghtptr-1)
               workl(ihd+rghtptr-1) = temp
               call dcopy (ncv, workl(iq+ncv*(leftptr-1)), 1,
     &                    workl(iw), 1)
               call dcopy (ncv, workl(iq+ncv*(rghtptr-1)), 1,
     &                    workl(iq+ncv*(leftptr-1)), 1)
               call dcopy (ncv, workl(iw), 1,
     &                    workl(iq+ncv*(rghtptr-1)), 1)
               leftptr = leftptr + 1
               rghtptr = rghtptr - 1
c
            end if
c
            if (leftptr .lt. rghtptr) go to 20
c
 30      end if
c
         if (msglvl .gt. 2) then
             call dvout  (logfil, ncv, workl(ihd), ndigit,
     &       '_seupd: The eigenvalues of H--reordered')
         end if
c
c        %----------------------------------------%
c        | Load the converged Ritz values into D. |
c        %----------------------------------------%
c
         call dcopy (nconv, workl(ihd), 1, d, 1)
c
      else
c
c        %-----------------------------------------------------%
c        | Ritz vectors not required. Load Ritz values into D. |
c        %-----------------------------------------------------%
c
         call dcopy (nconv, workl(ritz), 1, d, 1)
         call dcopy (ncv, workl(ritz), 1, workl(ihd), 1)
c
      end if
c
c     %------------------------------------------------------------------%
c     | Transform the Ritz values and possibly vectors and corresponding |
c     | Ritz estimates of OP to those of A*x=lambda*B*x. The Ritz values |
c     | (and corresponding data) are returned in ascending order.        |
c     %------------------------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
c        %---------------------------------------------------------%
c        | Ascending sort of wanted Ritz values, vectors and error |
c        | bounds. Not necessary if only Ritz values are desired.  |
c        %---------------------------------------------------------%
c
         if (rvec) then
            call dsesrt ('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call dcopy (ncv, workl(bounds), 1, workl(ihb), 1)
         end if
c
      else
c
c        %-------------------------------------------------------------%
c        | *  Make a copy of all the Ritz values.                      |
c        | *  Transform the Ritz values back to the original system.   |
c        |    For TYPE = 'SHIFTI' the transformation is                |
c        |             lambda = 1/theta + sigma                        |
c        |    For TYPE = 'BUCKLE' the transformation is                |
c        |             lambda = sigma * theta / ( theta - 1 )          |
c        |    For TYPE = 'CAYLEY' the transformation is                |
c        |             lambda = sigma * (theta + 1) / (theta - 1 )     |
c        |    where the theta are the Ritz values returned by dsaupd .  |
c        | NOTES:                                                      |
c        | *The Ritz vectors are not affected by the transformation.   |
c        |  They are only reordered.                                   |
c        %-------------------------------------------------------------%
c
         call dcopy  (ncv, workl(ihd), 1, workl(iw), 1)
         if (type .eq. 'SHIFTI') then
            do 40 k=1, ncv
               workl(ihd+k-1) = one / workl(ihd+k-1) + sigma
  40        continue
         else if (type .eq. 'BUCKLE') then
            do 50 k=1, ncv
               workl(ihd+k-1) = sigma * workl(ihd+k-1) /
     &                          (workl(ihd+k-1) - one)
  50        continue
         else if (type .eq. 'CAYLEY') then
            do 60 k=1, ncv
               workl(ihd+k-1) = sigma * (workl(ihd+k-1) + one) /
     &                          (workl(ihd+k-1) - one)
  60        continue
         end if
c
c        %-------------------------------------------------------------%
c        | *  Store the wanted NCONV lambda values into D.             |
c        | *  Sort the NCONV wanted lambda in WORKL(IHD:IHD+NCONV-1)   |
c        |    into ascending order and apply sort to the NCONV theta   |
c        |    values in the transformed system. We will need this to   |
c        |    compute Ritz estimates in the original system.           |
c        | *  Finally sort the lambda`s into ascending order and apply |
c        |    to Ritz vectors if wanted. Else just sort lambda`s into  |
c        |    ascending order.                                         |
c        | NOTES:                                                      |
c        | *workl(iw:iw+ncv-1) contain the theta ordered so that they  |
c        |  match the ordering of the lambda. We`ll use them again for |
c        |  Ritz vector purification.                                  |
c        %-------------------------------------------------------------%
c
         call dcopy (nconv, workl(ihd), 1, d, 1)
         call dsortr ('LA', .true., nconv, workl(ihd), workl(iw))
         if (rvec) then
            call dsesrt ('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call dcopy (ncv, workl(bounds), 1, workl(ihb), 1)
            call dscal (ncv, bnorm2/rnorm, workl(ihb), 1)
            call dsortr ('LA', .true., nconv, d, workl(ihb))
         end if
c
      end if
c
c     %------------------------------------------------%
c     | Compute the Ritz vectors. Transform the wanted |
c     | eigenvectors of the symmetric tridiagonal H by |
c     | the Lanczos basis matrix V.                    |
c     %------------------------------------------------%
c
      if (rvec .and. howmny .eq. 'A') then
c
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(iq,ldq).                                |
c        %----------------------------------------------------------%
c
         call dgeqr2 (ncv, nconv        , workl(iq) ,
     &                ldq, workl(iw+ncv), workl(ihb),
     &                ierr)
c
c        %--------------------------------------------------------%
c        | * Postmultiply V by Q.                                 |
c        | * Copy the first NCONV columns of VQ into Z.           |
c        | The N by NCONV matrix Z is now a matrix representation |
c        | of the approximate invariant subspace associated with  |
c        | the Ritz values in workl(ihd).                         |
c        %--------------------------------------------------------%
c
         call dorm2r ('Right', 'Notranspose', n        ,
     &                ncv    , nconv        , workl(iq),
     &                ldq    , workl(iw+ncv), v        ,
     &                ldv    , workd(n+1)   , ierr)
         call dlacpy ('All', n, nconv, v, ldv, z, ldz)
c
c        %-----------------------------------------------------%
c        | In order to compute the Ritz estimates for the Ritz |
c        | values in both systems, need the last row of the    |
c        | eigenvector matrix. Remember, it`s in factored form |
c        %-----------------------------------------------------%
c
         do 65 j = 1, ncv-1
            workl(ihb+j-1) = zero
  65     continue
         workl(ihb+ncv-1) = one
         temp_array(1) = temp
         call dorm2r ('Left', 'Transpose'  , ncv       ,
     &                1     , nconv        , workl(iq) ,
     &                ldq   , workl(iw+ncv), workl(ihb),
     &                ncv   , TEMP_ARRAY   , ierr)
c
      else if (rvec .and. howmny .eq. 'S') then
c
c     Not yet implemented. See remark 2 above.
c
      end if
c
      if (type .eq. 'REGULR' .and. rvec) then
c
            do 70 j=1, ncv
               workl(ihb+j-1) = rnorm * abs( workl(ihb+j-1) )
 70         continue
c
      else if (type .ne. 'REGULR' .and. rvec) then
c
c        %-------------------------------------------------%
c        | *  Determine Ritz estimates of the theta.       |
c        |    If RVEC = .true. then compute Ritz estimates |
c        |               of the theta.                     |
c        |    If RVEC = .false. then copy Ritz estimates   |
c        |              as computed by dsaupd .             |
c        | *  Determine Ritz estimates of the lambda.      |
c        %-------------------------------------------------%
c
         call dscal  (ncv, bnorm2, workl(ihb), 1)
         if (type .eq. 'SHIFTI') then
c
            do 80 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1) )
     &                        / workl(iw+k-1)**2
 80         continue
c
         else if (type .eq. 'BUCKLE') then
c
            do 90 k=1, ncv
               workl(ihb+k-1) = sigma * abs( workl(ihb+k-1) )
     &                        / (workl(iw+k-1)-one )**2
 90         continue
c
         else if (type .eq. 'CAYLEY') then
c
            do 100 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1)
     &                        / workl(iw+k-1)*(workl(iw+k-1)-one) )
 100        continue
c
         end if
c
      end if
c
      if (type .ne. 'REGULR' .and. msglvl .gt. 1) then
         call dvout (logfil, nconv, d, ndigit,
     &          '_seupd: Untransformed converged Ritz values')
         call dvout (logfil, nconv, workl(ihb), ndigit,
     &     '_seupd: Ritz estimates of the untransformed Ritz values')
      else if (msglvl .gt. 1) then
         call dvout (logfil, nconv, d, ndigit,
     &          '_seupd: Converged Ritz values')
         call dvout (logfil, nconv, workl(ihb), ndigit,
     &     '_seupd: Associated Ritz estimates')
      end if
c
c     %-------------------------------------------------%
c     | Ritz vector purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 3,4,5. See reference 7               |
c     %-------------------------------------------------%
c
      if (rvec .and. (type .eq. 'SHIFTI' .or. type .eq. 'CAYLEY')) then
c
         do 110 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1)
     &                  / workl(iw+k)
 110     continue
c
      else if (rvec .and. type .eq. 'BUCKLE') then
c
         do 120 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1)
     &                  / (workl(iw+k)-one)
 120     continue
c
      end if
c
      if (type .ne. 'REGULR')
     &   call dger  (n, nconv, one, resid, 1, workl(iw), 1, z, ldz)
c
 9000 continue
c

      RETURN

! **********************************************************************************************************************************

c
c     %---------------%
c     | End of dseupd |
c     %---------------%
c
      end subroutine dseupd

! ##################################################################################################################################
! 004 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsaup2
c
c\Description:
c  Intermediate level interface called by dsaupd.
c
c\Usage:
c  call dsaup2
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL,
c       IPNTR, WORKD, INFO, PITERS )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dsaupd.
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in dsaupd.
c
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during
c          each Arnoldi/Lanczos iteration.
c          If ISHIFT=1, NP is adjusted dynamically at each iteration
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user needs
c          to provide via reverse comunication. 0 < NP < NCV-NEV.
c          NP may be less than NCV-NEV since a leading block of the current
c          upper Tridiagonal matrix has split off and contains "unwanted"
c          Ritz values.
c          Upon termination of the IRA iteration, NP contains the number
c          of "converged" wanted Ritz values.
c
c  IUPD    Integer.  (INPUT)
c          IUPD .EQ. 0: use explicit restart instead implicit update.
c          IUPD .NE. 0: use implicit update.
c
c  V       REAL(DOUBLE) N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Lanczos basis vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       REAL(DOUBLE) (NEV+NP) by 2 array.  (OUTPUT)
c          H is used to store the generated symmetric tridiagonal matrix
c          The subdiagonal is stored in the first column of H starting
c          at H(2,1).  The main diagonal is stored in the second column
c          of H starting at H(1,2). If dsaup2 converges store the
c          B-norm of the final residual vector in H(1,1).
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RITZ    REAL(DOUBLE) array of length NEV+NP.  (OUTPUT)
c          RITZ(1:NEV) contains the computed Ritz values of OP.
c
c  BOUNDS  REAL(DOUBLE) array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to RITZ.
c
c  Q       REAL(DOUBLE) (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   REAL(DOUBLE) array of length at least 3*(NEV+NP).  (INPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in the computation of the
c          tridiagonal eigenvalue problem, the calculation and
c          application of the shifts and convergence checking.
c          If ISHIFT .EQ. O and IDO .EQ. 3, the first NP locations
c          of WORKL are used in reverse communication to hold the user
c          supplied shifts.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for
c          vectors used by the Lanczos iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in one of
c                    the spectral transformation modes.  X is the current
c                    operand.
c          -------------------------------------------------------------
c
c  WORKD   REAL(DOUBLE) work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Lanczos iteration
c          for reverse communication.  The user should not use WORKD
c          See Data Distribution Note in dsaupd.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: All possible eigenvalues of OP has been found.
c                   NP returns the size of the invariant subspace
c                   spanning the operator OP.
c          =     2: No shifts could be applied.
c          =    -8: Error return from trid. eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Lanczos factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
c     1980.
c  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
c     Computer Physics Communications, 53 (1989), pp 169-179.
c  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
c     Implement the Spectral Transformation", Math. Comp., 48 (1987),
c     pp 663-673.
c  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos
c     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems",
c     SIAM J. Matr. Anal. Apps.,  January (1993).
c  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
c     for Updating the QR decomposition", ACM TOMS, December 1990,
c     Volume 16 Number 4, pp 369-377.
c
c\Routines called:
c     dgetv0  ARPACK initial vector generation routine.
c     dsaitr  ARPACK Lanczos factorization routine.
c     dsapps  ARPACK application of implicit shifts routine.
c     dsconv  ARPACK convergence of Ritz values routine.
c     dseigt  ARPACK compute Ritz values and error bounds routine.
c     dsgets  ARPACK reorder Ritz values and error bounds routine.
c     dsortr  ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dscal   Level 1 BLAS that scales a vector.
c     dswap   Level 1 BLAS that swaps two vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     12/15/93: Version ' 2.4'
c     xx/xx/95: Version ' 2.4'.  (R.B. Lehoucq)
c
c\SCCS Information: @(#)
c FILE: saup2.F   SID: 2.7   DATE OF SID: 5/19/98   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsaup2
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, h, ldh, ritz, bounds,
! B 05/24/04 //////////////////////////////////////////////////////////B
!    &     q, ldq, workl, ipntr, workd, info )
!    &     q, ldq, workl, ipntr, workd, info, nev_user, sigma, eig_old )
! E ///////////////////////////////////////////////////////////////////E
! B 06/08/08 //////////////////////////////////////////////////////////B
     &     q, ldq, workl, ipntr, workd, info, nev_user, sigma, eig_old,
     &     iter, piters )   ! 11/30/09: add piters
! E ///////////////////////////////////////////////////////////////////E

! B 05/24/04 //////////////////////////////////////////////////////////B

      integer(long) n, numout, num_left, jj, jstart, nev_user
      real(double)  sigma
      real(double) eigout(n)
      real(double) eig_new(n)
      real(double) eig_old(n)
      real(double) eig_pc(n)

! E ///////////////////////////////////////////////////////////////////E

! B 06/08/08 //////////////////////////////////////////////////////////B
      integer iter_old
! E ///////////////////////////////////////////////////////////////////E

! B 11/30/09 //////////////////////////////////////////////////////////B
      integer piters
! E ///////////////////////////////////////////////////////////////////E

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSAUP2'
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cxxxx include   'debug.h'
cxxxx include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, ldh, ldq, ldv, mxiter,
     &           mode, nev, nev_array(1), np, np_array(1)
      REAL(DOUBLE)
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      REAL(DOUBLE)
     &           bounds(nev+np), h(ldh,2), q(ldq,nev+np), resid(n),
     &           ritz(nev+np), v(ldv,nev+np), workd(3*n),
     &           workl(3*(nev+np))
c
c     %------------%
c     | Parameters |
c     %------------%
c
      REAL(DOUBLE)
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  wprime*2
      logical    cnorm, getv0, initv, update, ushift
      integer    ierr, iter, j, kplusp, msglvl, nconv, nevbef, nev0,
     &           np0, nptemp, nevd2, nevm2, kp(3),
     &           iter_array(1), nconv_array(1)
      REAL(DOUBLE)
     &           rnorm, temp, eps23
      REAL(DOUBLE)
     &           RNORM_ARRAY(1)
      save       cnorm, getv0, initv, update, ushift,
! B 06/08/08 //////////////////////////////////////////////////////////B
!    &           iter, kplusp, msglvl, nconv, nev0, np0,
     &                 kplusp, msglvl, nconv, nev0, np0,
! E ///////////////////////////////////////////////////////////////////E
     &           rnorm, eps23
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
!:!   external   dcopy, dgetv0, dsaitr, dscal, dsconv, dseigt, dsgets,
!:!  &           dsapps, dsortr, dvout, ivout, second, dswap
      external   dswap, dscal

c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      REAL(DOUBLE)
     &           ddot, dnrm2
      external   ddot, dnrm2

      REAL(DOUBLE) dlamch
      external     dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c


! **********************************************************************************************************************************

      if (ido .eq. 0) then
c
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call cpu_time (t0)
         msglvl = msaup2
c
c        %---------------------------------%
c        | Set machine dependent constant. |
c        %---------------------------------%
c
         eps23 = dlamch('Epsilon-Machine')
         eps23 = eps23**(2.0D+0/3.0D+0)
c
c        %-------------------------------------%
c        | nev0 and np0 are integer variables  |
c        | hold the initial values of NEV & NP |
c        %-------------------------------------%
c
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvlues.      |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev0 + np0
         nconv  = 0
         iter   = 0
c
c        %--------------------------------------------%
c        | Set flags for computing the first NEV steps |
c        | of the Lanczos factorization.              |
c        %--------------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
c
c        %--------------------------------------------%
c        | User provides the initial residual vector. |
c        %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
c
   10 continue
c
      if (getv0) then
         call dgetv0 (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. zero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. |
c           %-----------------------------------------%
c
            info = -9
            go to 1200
         end if
         getv0 = .false.
         ido  = 0
      end if
c
c     %------------------------------------------------------------%
c     | Back from reverse communication: continue with update step |
c     %------------------------------------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Lanczos factorization |
c     %----------------------------------------------------------%
c
      call dsaitr (ido, bmat, n, 0, nev0, mode, resid, rnorm, v, ldv,
     &             h, ldh, ipntr, workd, info)
c
c     %---------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication  |
c     | to compute operations involving OP and possibly B |
c     %---------------------------------------------------%
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
c
c        %-----------------------------------------------------%
c        | dsaitr was unable to build an Lanczos factorization |
c        | of length NEV0. INFO is returned with the size of   |
c        | the factorization built. Exit main loop.            |
c        %-----------------------------------------------------%
c
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  LANCZOS  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Lanczos     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c
 1000 continue
c
         iter = iter + 1
c
         iter_array(1) = iter
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, iter_array, ndigit,
     &           '_saup2: **** Start of major iteration number ****')
         end if
         if (msglvl .gt. 1) then
            nev_array(1) = nev
            np_array(1)  = np
            call ivout (logfil, 1, nev_array, ndigit,
     &     '_saup2: The length of the current Lanczos factorization')
            call ivout (logfil, 1, np_array, ndigit,
     &           '_saup2: Extend the Lanczos factorization by')
         end if
c
c        %------------------------------------------------------------%
c        | Compute NP additional steps of the Lanczos factorization. |
c        %------------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call dsaitr (ido, bmat, n, nev, np, mode, resid, rnorm, v,
     &                ldv, h, ldh, ipntr, workd, info)
c
c        %---------------------------------------------------%
c        | ido .ne. 99 implies use of reverse communication  |
c        | to compute operations involving OP and possibly B |
c        %---------------------------------------------------%
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
c
c           %-----------------------------------------------------%
c           | dsaitr was unable to build an Lanczos factorization |
c           | of length NEV0+NP0. INFO is returned with the size  |
c           | of the factorization built. Exit main loop.         |
c           %-----------------------------------------------------%
c
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
         if (msglvl .gt. 1) then
            rnorm_array(1) = rnorm
            call dvout (logfil, 1, rnorm_array, ndigit,
     &           '_saup2: Current B-norm of residual for factorization')
         end if
c
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current symmetric tridiagonal matrix.           |
c        %--------------------------------------------------------%
c
         call dseigt (rnorm, kplusp, h, ldh, ritz, bounds, workl, ierr)

! B 05/24/04 //////////////////////////////////////////////////////////B
      if (piters > 0) then  ! 11/30/09: change debug(47) to piters

         do j=1,kplusp
            jj = kplusp - j + 1
            if (dabs(ritz(j)) > dlamch('s')) then
               if (mode == 2) then
                  eigout(jj) = sigma + one/ritz(j)
               else
                  eigout(jj) = ritz(j)
               endif
            else
               if (mode == 2) then
                  if (ritz(j) >= 0.d0) then
                     eigout(jj) = sigma + 1.0d0/dlamch('s')
                  else
                     eigout(jj) = sigma - 1.0d0/dlamch('s')
                  endif
               else
                  eigout(j) = ritz(j)
               endif
            endif
            eig_new(jj) = eigout(jj)
            eig_pc(jj)  = 100.d0*(eig_new(jj) - eig_old(jj))/eig_new(jj)
            eig_old(jj) = eig_new(jj)
         enddo

         numout = min(nev_user,kplusp)
         if (numout <= 10) then
            Write(f06,99101)iter,(eigout(j),j=1,numout)
         else
            Write(f06,99101)iter,(eigout(j),j=1,10)
            num_left = numout - 10
            jstart = 11
            do j=11,numout
               Write(f06,99102) (eigout(jj),jj=jstart,jstart+9)
               jstart   = jstart + 10
               num_left = num_left - 10
               if (num_left < 10) exit
            enddo
            if (num_left > 0) then
               Write(f06,99102) (eigout(j),j=jstart,numout)
            endif
            Write(f06,*)
         endif

      endif

99101 format(' Eigenvalue estimates in Lanczos iteration ',i5,': ',
     &         10(1es14.6))

99102 format(48x,10(1es14.6))

!9103 format(48x,10(a,1es13.5))

! E////////////////////////////////////////////////////////////////////E
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------%
c        | Make a copy of eigenvalues and corresponding error |
c        | bounds obtained from _seigt.                       |
c        %----------------------------------------------------%
c
         call dcopy(kplusp, ritz, 1, workl(kplusp+1), 1)
         call dcopy(kplusp, bounds, 1, workl(2*kplusp+1), 1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The selection is based on the requested number of |
c        | eigenvalues instead of the current NEV and NP to  |
c        | prevent possible misconvergence.                  |
c        | * Wanted Ritz values := RITZ(NP+1:NEV+NP)         |
c        | * Shifts := RITZ(1:NP) := WORKL(1:NP)             |
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
         call dsgets (ishift, which, nev, np, ritz, bounds, workl)
c
c        %-------------------%
c        | Convergence test. |
c        %-------------------%
c
         call dcopy (nev, bounds(np+1), 1, workl(np+1), 1)
         call dsconv (nev, ritz(np+1), workl(np+1), tol, nconv)
c
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = nconv
            call ivout (logfil, 3, kp, ndigit,
     &                  '_saup2: NEV, NP, NCONV are')
            call dvout (logfil, kplusp, ritz, ndigit,
     &           '_saup2: The eigenvalues of H')
            call dvout (logfil, kplusp, bounds, ndigit,
     &          '_saup2: Ritz estimates of the current NCV Ritz values')
         end if
c
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c
         if ( (nconv .ge. nev0) .or.
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP since we don't want to |
c           | swap overlapping locations.                    |
c           %------------------------------------------------%
c
            if (which .eq. 'BE') then
c
c              %-----------------------------------------------------%
c              | Both ends of the spectrum are requested.            |
c              | Sort the eigenvalues into algebraically decreasing  |
c              | order first then swap low end of the spectrum next  |
c              | to high end in appropriate locations.               |
c              | NOTE: when np < floor(nev/2) be careful not to swap |
c              | overlapping locations.                              |
c              %-----------------------------------------------------%
c
               wprime = 'SA'
               call dsortr (wprime, .true., kplusp, ritz, bounds)
               nevd2 = nev0 / 2
               nevm2 = nev0 - nevd2
               if ( nev .gt. 1 ) then
                  call dswap ( min(nevd2,np), ritz(nevm2+1), 1,
     &                 ritz( max(kplusp-nevd2+1,kplusp-np+1) ), 1)
                  call dswap ( min(nevd2,np), bounds(nevm2+1), 1,
     &                 bounds( max(kplusp-nevd2+1,kplusp-np+1)), 1)
               end if
c
            else
c
c              %--------------------------------------------------%
c              | LM, SM, LA, SA case.                             |
c              | Sort the eigenvalues of H into the an order that |
c              | is opposite to WHICH, and apply the resulting    |
c              | order to BOUNDS.  The eigenvalues are sorted so  |
c              | that the wanted part are always within the first |
c              | NEV locations.                                   |
c              %--------------------------------------------------%
c
               if (which .eq. 'LM') wprime = 'SM'
               if (which .eq. 'SM') wprime = 'LM'
               if (which .eq. 'LA') wprime = 'SA'
               if (which .eq. 'SA') wprime = 'LA'
c
               call dsortr (wprime, .true., kplusp, ritz, bounds)
c
            end if
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23,magnitude of the Ritz value).   |
c           %--------------------------------------------------%
c
            do 35 j = 1, nev0
               temp = max( eps23, abs(ritz(j)) )
               bounds(j) = bounds(j)/temp
 35         continue
c
c           %----------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz  |
c           | esitmates.  This will push all the converged ones  |
c           | towards the front of ritzr, ritzi, bounds          |
c           | (in the case when NCONV < NEV.)                    |
c           %----------------------------------------------------%
c
            wprime = 'LA'
            call dsortr(wprime, .true., nev0, bounds, ritz)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, nev0
                temp = max( eps23, abs(ritz(j)) )
                bounds(j) = bounds(j)*temp
 40         continue
c
c           %--------------------------------------------------%
c           | Sort the "converged" Ritz values again so that   |
c           | the "threshold" values and their associated Ritz |
c           | estimates appear at the appropriate position in  |
c           | ritz and bound.                                  |
c           %--------------------------------------------------%
c
            if (which .eq. 'BE') then
c
c              %------------------------------------------------%
c              | Sort the "converged" Ritz values in increasing |
c              | order.  The "threshold" values are in the      |
c              | middle.                                        |
c              %------------------------------------------------%
c
               wprime = 'LA'
               call dsortr(wprime, .true., nconv, ritz, bounds)
c
            else
c
c              %----------------------------------------------%
c              | In LM, SM, LA, SA case, sort the "converged" |
c              | Ritz values according to WHICH so that the   |
c              | "threshold" value appears at the front of    |
c              | ritz.                                        |
c              %----------------------------------------------%

               call dsortr(which, .true., nconv, ritz, bounds)
c
            end if
c
c           %------------------------------------------%
c           |  Use h( 1,1 ) as storage to communicate  |
c           |  rnorm to _seupd if needed               |
c           %------------------------------------------%
c
            h(1,1) = rnorm
c
            if (msglvl .gt. 1) then
               call dvout (logfil, kplusp, ritz, ndigit,
     &            '_saup2: Sorted Ritz values.')
               call dvout (logfil, kplusp, bounds, ndigit,
     &            '_saup2: Sorted ritz estimates.')
            end if
c
c           %------------------------------------%
c           | Max iterations have been exceeded. |
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. nev) info = 1
c
c           %---------------------%
c           | No shifts to apply. |
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. nev0) info = 2
c
            np = nconv
            go to 1100
c
         else if (nconv .lt. nev .and. ishift .eq. 1) then
c
c           %---------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.    |
c           | To prevent possible stagnation, adjust the number |
c           | of Ritz values and the shifts.                    |
c           %---------------------------------------------------%
c
            nevbef = nev
            nev = nev + min (nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 2) then
               nev = 2
            end if
            np  = kplusp - nev
c
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c
            if (nevbef .lt. nev)
     &         call dsgets (ishift, which, nev, np, ritz, bounds,
     &              workl)
c
         end if
c
         if (msglvl .gt. 0) then
            nconv_array(1) = nconv
            call ivout (logfil, 1, nconv_array, ndigit,
     &           '_saup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit,
     &              '_saup2: NEV and NP are')
               call dvout (logfil, nev, ritz(np+1), ndigit,
     &              '_saup2: "wanted" Ritz values.')
               call dvout (logfil, nev, bounds(np+1), ndigit,
     &              '_saup2: Ritz estimates of the "wanted" values ')
            end if
         end if

c
         if (ishift .eq. 0) then
c
c           %-----------------------------------------------------%
c           | User specified shifts: reverse communication to     |
c           | compute the shifts. They are returned in the first  |
c           | NP locations of WORKL.                              |
c           %-----------------------------------------------------%
c
            ushift = .true.
            ido = 3
            go to 9000
         end if
c
   50    continue
c
c        %------------------------------------%
c        | Back from reverse communication;   |
c        | User specified shifts are returned |
c        | in WORKL(1:*NP)                   |
c        %------------------------------------%
c
         ushift = .false.
c
c
c        %---------------------------------------------------------%
c        | Move the NP shifts to the first NP locations of RITZ to |
c        | free up WORKL.  This is for the non-exact shift case;   |
c        | in the exact shift case, dsgets already handles this.   |
c        %---------------------------------------------------------%
c
         if (ishift .eq. 0) call dcopy (np, workl, 1, ritz, 1)
c
         if (msglvl .gt. 2) then
            np_array(1) = np
            call ivout (logfil, 1, np_array, ndigit,
     &                  '_saup2: The number of shifts to apply ')
            call dvout (logfil, np, workl, ndigit,
     &                  '_saup2: shifts selected')
            if (ishift .eq. 1) then
               call dvout (logfil, np, bounds, ndigit,
     &                  '_saup2: corresponding Ritz estimates')
             end if
         end if
c
c        %---------------------------------------------------------%
c        | Apply the NP0 implicit shifts by QR bulge chasing.      |
c        | Each shift is applied to the entire tridiagonal matrix. |
c        | The first 2*N locations of WORKD are used as workspace. |
c        | After dsapps is done, we have a Lanczos                 |
c        | factorization of length NEV.                            |
c        %---------------------------------------------------------%
c
         call dsapps (n, nev, np, ritz, v, ldv, h, ldh, resid, q, ldq,
     &        workd)
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to dsaitr.  |
c        %---------------------------------------------%
c
         cnorm = .true.
         call cpu_time (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd, 1)
         end if
c
  100    continue
c
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c
         if (bmat .eq. 'G') then
            rnorm = ddot (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
         cnorm = .false.
  130    continue
c
         if (msglvl .gt. 2) then
            rnorm_array(1) = rnorm
            call dvout (logfil, 1, rnorm_array, ndigit,
     &      '_saup2: B-norm of residual for NEV factorization')
            call dvout (logfil, nev, h(1,2), ndigit,
     &           '_saup2: main diagonal of compressed H matrix')
            call dvout (logfil, nev-1, h(2,1), ndigit,
     &           '_saup2: subdiagonal of compressed H matrix')
         end if
c
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 1100 continue
c
      mxiter = iter
      nev = nconv
c
 1200 continue
      ido = 99
! B 05/24/04 //////////////////////////////////////////////////////////B
                                 ! 11/30/09: change debug(47) to piters
      if ((piters > 0) .and. (mxiter > 1)) then
         Write(f06,*)
         if (numout <= 10) then
            Write(f06,99201) (eig_pc(j),j=1,numout)
         else
            Write(f06,99201 ) (eig_pc(j),j=1,10)
            num_left = numout - 10
            jstart = 11
            do j=11,numout
               Write(f06,99202) (eig_pc(jj),jj=jstart,jstart+9)
               jstart   = jstart + 10
               num_left = num_left - 10
               if (num_left < 10) exit
            enddo
            if (num_left > 0) then
               Write(f06,99202) (eig_pc(j),j=jstart,numout)
            endif
         endif
      endif

99201 format(' % diffs in eigenvals from last 2 iterations    : ',
     &         10(1es14.6))

99202 format(48x,10(1es14.6))

! E////////////////////////////////////////////////////////////////////E
c
c     %------------%
c     | Error exit |
c     %------------%
c
      call cpu_time (t1)
      tsaup2 = t1 - t0
c
 9000 continue


      RETURN

! **********************************************************************************************************************************
c
c     %---------------%
c     | End of dsaup2 |
c     %---------------%
c
      end subroutine dsaup2

! ##################################################################################################################################
! 005 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsesrt
c
c\Description:
c  Sort the array X in the order specified by WHICH and optionally
c  apply the permutation to the columns of the matrix A.
c
c\Usage:
c  call dsesrt
c     ( WHICH, APPLY, N, X, NA, A, LDA)
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> X is sorted into increasing order of magnitude.
c          'SM' -> X is sorted into decreasing order of magnitude.
c          'LA' -> X is sorted into increasing order of algebraic.
c          'SA' -> X is sorted into decreasing order of algebraic.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to A.
c          APPLY = .FALSE. -> do not apply the sorted order to A.
c
c  N       Integer.  (INPUT)
c          Dimension of the array X.
c
c  X      REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
c          The array to be sorted.
c
c  NA      Integer.  (INPUT)
c          Number of rows of the matrix A.
c
c  A      REAL(DOUBLE) array of length NA by N.  (INPUT/OUTPUT)
c
c  LDA     Integer.  (INPUT)
c          Leading dimension of A.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Routines
c     dswap  Level 1 BLAS that swaps the contents of two vectors.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     12/15/93: Version ' 2.1'.
c               Adapted from the sort routine in LANSO and
c               the ARPACK code dsortr
c
c\SCCS Information: @(#)
c FILE: sesrt.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsesrt (which, apply, n, x, na, a, lda)

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSESRT'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    lda, n, na
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      REAL(DOUBLE)
     &           x(0:n-1), a(lda, 0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      REAL(DOUBLE)
     &           temp
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
!:!   external   dswap
      external   dswap
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c


! **********************************************************************************************************************************

      igap = n / 2
c
      if (which .eq. 'SA') then
c
c        X is sorted into decreasing order of algebraic.
c
   10    continue
         if (igap .eq. 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            if (x(j).lt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        X is sorted into decreasing order of magnitude.
c
   40    continue
         if (igap .eq. 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j.lt.0) go to 60
c
            if (abs(x(j)).lt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
c
      else if (which .eq. 'LA') then
c
c        X is sorted into increasing order of algebraic.
c
   70    continue
         if (igap .eq. 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c
            if (x(j).gt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
c
      else if (which .eq. 'LM') then
c
c        X is sorted into increasing order of magnitude.
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (abs(x(j)).gt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call dswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
      end if
c
 9000 continue


      RETURN

! **********************************************************************************************************************************
c
c     %---------------%
c     | End of dsesrt |
c     %---------------%
c
      end subroutine dsesrt

! ##################################################################################################################################
! 006 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsgets
c
c\Description:
c  Given the eigenvalues of the symmetric tridiagonal matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: This is called even in the case of user specified shifts in
c  order to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call dsgets
c     ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> KEV eigenvalues of largest magnitude are retained.
c          'SM' -> KEV eigenvalues of smallest magnitude are retained.
c          'LA' -> KEV eigenvalues of largest value are retained.
c          'SA' -> KEV eigenvalues of smallest value are retained.
c          'BE' -> KEV eigenvalues, half from each end of the spectrum.
c                  If KEV is odd, compute one more from the high end.
c
c  KEV      Integer.  (INPUT)
c          KEV+NP is the size of the matrix H.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be computed.
c
c  RITZ    REAL(DOUBLE) array of length KEV+NP.  (INPUT/OUTPUT)
c          On INPUT, RITZ contains the eigenvalues of H.
c          On OUTPUT, RITZ are sorted so that the unwanted eigenvalues
c          are in the first NP locations and the wanted part is in
c          the last KEV locations.  When exact shifts are selected, the
c          unwanted part corresponds to the shifts to be applied.
c
c  BOUNDS  REAL(DOUBLE) array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  SHIFTS  REAL(DOUBLE) array of length NP.  (INPUT/OUTPUT)
c          On INPUT:  contains the user specified shifts if ISHIFT = 0.
c          On OUTPUT: contains the shifts sorted into decreasing order
c          of magnitude with respect to the Ritz estimates contained in
c          BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dsortr  ARPACK utility sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dswap   Level 1 BLAS that swaps the contents of two vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/93: Version ' 2.1'
c
c\SCCS Information: @(#)
c FILE: sgets.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsgets ( ishift, which, kev, np, ritz, bounds, shifts )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSGETS'
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cxxxx include   'debug.h'
cxxxx include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      integer    ishift, kev, np
      INTEGER    KEV_ARRAY(1), NP_ARRAY(1)
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      REAL(DOUBLE)
     &           bounds(kev+np), ritz(kev+np), shifts(np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      REAL(DOUBLE)
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    kevd2, msglvl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
!:!   external   dswap, dcopy, dsortr, second
      external   dswap
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    max, min
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c


! **********************************************************************************************************************************

c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call cpu_time (t0)
      msglvl = msgets
c
      if (which .eq. 'BE') then
c
c        %-----------------------------------------------------%
c        | Both ends of the spectrum are requested.            |
c        | Sort the eigenvalues into algebraically increasing  |
c        | order first then swap high end of the spectrum next |
c        | to low end in appropriate locations.                |
c        | NOTE: when np < floor(kev/2) be careful not to swap |
c        | overlapping locations.                              |
c        %-----------------------------------------------------%
c
         call dsortr ('LA', .true., kev+np, ritz, bounds)
         kevd2 = kev / 2
         if ( kev .gt. 1 ) then
            call dswap ( min(kevd2,np), ritz, 1,
     &                   ritz( max(kevd2,np)+1 ), 1)
            call dswap ( min(kevd2,np), bounds, 1,
     &                   bounds( max(kevd2,np)+1 ), 1)
         end if
c
      else
c
c        %----------------------------------------------------%
c        | LM, SM, LA, SA case.                               |
c        | Sort the eigenvalues of H into the desired order   |
c        | and apply the resulting order to BOUNDS.           |
c        | The eigenvalues are sorted so that the wanted part |
c        | are always in the last KEV locations.               |
c        %----------------------------------------------------%
c
         call dsortr (which, .true., kev+np, ritz, bounds)
      end if
c
      if (ishift .eq. 1 .and. np .gt. 0) then
c
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first.       |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when the shifts  |
c        | are applied in subroutine dsapps.                     |
c        %-------------------------------------------------------%
c
         call dsortr ('SM', .true., np, bounds, ritz)
         call dcopy (np, ritz, 1, shifts, 1)
      end if
c
      call cpu_time (t1)
      tsgets = tsgets + (t1 - t0)
c
      if (msglvl .gt. 0) then
         kev_array(1) = kev
         np_array(1)  = np
         call ivout (logfil, 1, KEV_ARRAY, ndigit, '_sgets: KEV is')
         call ivout (logfil, 1, NP_ARRAY, ndigit, '_sgets: NP is')
         call dvout (logfil, kev+np, ritz, ndigit,
     &        '_sgets: Eigenvalues of current H matrix')
         call dvout (logfil, kev+np, bounds, ndigit,
     &        '_sgets: Associated Ritz estimates')
      end if
c

      RETURN

! **********************************************************************************************************************************
c
c     %---------------%
c     | End of dsgets |
c     %---------------%
c
      end subroutine dsgets

! ##################################################################################################################################
! 007 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsortr
c
c\Description:
c  Sort the array X1 in the order specified by WHICH and optionally
c  applies the permutation to the array X2.
c
c\Usage:
c  call dsortr
c     ( WHICH, APPLY, N, X1, X2 )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> X1 is sorted into increasing order of magnitude.
c          'SM' -> X1 is sorted into decreasing order of magnitude.
c          'LA' -> X1 is sorted into increasing order of algebraic.
c          'SA' -> X1 is sorted into decreasing order of algebraic.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to X2.
c          APPLY = .FALSE. -> do not apply the sorted order to X2.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  X1      REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
c          The array to be sorted.
c
c  X2      REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
c          Only referenced if APPLY = .TRUE.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     12/16/93: Version ' 2.1'.
c               Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#)
c FILE: sortr.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsortr (which, apply, n, x1, x2)

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSORTR'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      REAL(DOUBLE)
     &           x1(0:n-1), x2(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      REAL(DOUBLE)
     &           temp
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c


! **********************************************************************************************************************************

      igap = n / 2
c
      if (which .eq. 'SA') then
c
c        X1 is sorted into decreasing order of algebraic.
c
   10    continue
         if (igap .eq. 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            if (x1(j).lt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        X1 is sorted into decreasing order of magnitude.
c
   40    continue
         if (igap .eq. 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j.lt.0) go to 60
c
            if (abs(x1(j)).lt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
c
      else if (which .eq. 'LA') then
c
c        X1 is sorted into increasing order of algebraic.
c
   70    continue
         if (igap .eq. 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c
            if (x1(j).gt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
c
      else if (which .eq. 'LM') then
c
c        X1 is sorted into increasing order of magnitude.
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (abs(x1(j)).gt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
      end if
c
 9000 continue
c

      RETURN

! **********************************************************************************************************************************

c     %---------------%
c     | End of dsortr |
c     %---------------%
c
      end subroutine dsortr

! ##################################################################################################################################
! 008 ARPACK_LANCZOS_EIG

c
c\SCCS Information: @(#)
c FILE: stats.F   SID: 2.1   DATE OF SID: 4/19/96   RELEASE: 2
c     %---------------------------------------------%
c     | Initialize statistic and timing information |
c     | for symmetric Arnoldi code.                 |
c     %---------------------------------------------%

      subroutine dstats

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSTATS'

c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
cxxxx include   'stat.h'

c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%



! **********************************************************************************************************************************

      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0

      tsaupd = 0.0D+0
      tsaup2 = 0.0D+0
      tsaitr = 0.0D+0
      tseigt = 0.0D+0
      tsgets = 0.0D+0
      tsapps = 0.0D+0
      tsconv = 0.0D+0
      titref = 0.0D+0
      tgetv0 = 0.0D+0
      trvec  = 0.0D+0

c     %----------------------------------------------------%
c     | User time including reverse communication overhead |
c     %----------------------------------------------------%
      tmvopx = 0.0D+0
      tmvbx  = 0.0D+0


      RETURN

! **********************************************************************************************************************************
c
c     End of dstats
c
      end subroutine dstats

! ##################################################################################################################################
! 009 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dgetv0
c
c\Description:
c  Generate a random initial residual vector for the Arnoldi process.
c  Force the residual vector to be in the range of the operator OP.
c
c\Usage:
c  call dgetv0
c     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM,
c       IPNTR, WORKD, IERR )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first
c          call to dgetv0.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B in the (generalized)
c          eigenvalue problem A*x = lambda*B*x.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  ITRY    Integer.  (INPUT)
c          ITRY counts the number of times that dgetv0 is called.
c          It should be set to 1 on the initial call to dgetv0.
c
c  INITV   Logical variable.  (INPUT)
c          .TRUE.  => the initial residual vector is given in RESID.
c          .FALSE. => generate a random initial residual vector.
c
c  N       Integer.  (INPUT)
c          Dimension of the problem.
c
c  J       Integer.  (INPUT)
c          Index of the residual vector to be generated, with respect to
c          the Arnoldi process.  J > 1 in case of a "restart".
c
c  V       REAL(DOUBLE) N by J array.  (INPUT)
c          The first J-1 columns of V contain the current Arnoldi basis
c          if this is a "restart".
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  RESID   REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
c          Initial residual vector to be generated.  If RESID is
c          provided, force RESID into the range of the operator OP.
c
c  RNORM   REAL(DOUBLE) scalar.  (OUTPUT)
c          B-norm of the generated residual.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c
c  WORKD   REAL(DOUBLE) work array of length 2*N.  (REVERSE COMMUNICATION).
c          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
c
c  IERR    Integer.  (OUTPUT)
c          =  0: Normal exit.
c          = -1: Cannot generate a nontrivial restarted residual vector
c                in the range of the operator OP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine for vector output.
c     dlarnv  LAPACK routine for generating a random vector.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: getv0.F   SID: 2.7   DATE OF SID: 04/07/99   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dgetv0
     &   ( ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm,
     &     ipntr, workd, ierr )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DGETV0'

c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cxxxx include   'debug.h'
cxxxx include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      logical    initv
      integer    ido, ierr, itry, j, ldv, n
      REAL(DOUBLE)
     &           rnorm
      REAL(DOUBLE)  RNORM_ARRAY(1)
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      REAL(DOUBLE)
     &           resid(n), v(ldv,j), workd(2*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      REAL(DOUBLE)
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      logical    first, inits, orth
      integer    idist, iseed(4), iter, msglvl, jj
      REAL(DOUBLE)
     &           rnorm0
      REAL(DOUBLE)  RNORM0_ARRAY(1)
      save       first, iseed, inits, iter, msglvl, orth, rnorm0
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
!:!   external   dlarnv, dvout, dcopy, dgemv, second
      external   dgemv
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      REAL(DOUBLE)
     &           ddot, dnrm2
      external   ddot, dnrm2
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, sqrt
c
c     %-----------------%
c     | Data Statements |
c     %-----------------%
c
      data       inits /.true./
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c


! **********************************************************************************************************************************
c
c     %-----------------------------------%
c     | Initialize the seed of the LAPACK |
c     | random number generator           |
c     %-----------------------------------%
c
      if (inits) then
          iseed(1) = 1
          iseed(2) = 3
          iseed(3) = 5
          iseed(4) = 7
          inits = .false.
      end if
c
      if (ido .eq.  0) then
c
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call cpu_time (t0)
         msglvl = mgetv0
c
         ierr   = 0
         iter   = 0
         first  = .FALSE.
         orth   = .FALSE.
c
c        %-----------------------------------------------------%
c        | Possibly generate a random starting vector in RESID |
c        | Use a LAPACK random number generator used by the    |
c        | matrix generation routines.                         |
c        |    idist = 1: uniform (0,1)  distribution;          |
c        |    idist = 2: uniform (-1,1) distribution;          |
c        |    idist = 3: normal  (0,1)  distribution;          |
c        %-----------------------------------------------------%
c
         if (.not.initv) then
            idist = 2
            call dlarnv (idist, iseed, n, resid)
         end if
c
c        %----------------------------------------------------------%
c        | Force the starting vector into the range of OP to handle |
c        | the generalized problem when B is possibly (singular).   |
c        %----------------------------------------------------------%
c
         call cpu_time (t2)
         if (bmat .eq. 'G') then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call dcopy (n, resid, 1, workd, 1)
            ido = -1
            go to 9000
         end if
      end if
c
c     %-----------------------------------------%
c     | Back from computing OP*(initial-vector) |
c     %-----------------------------------------%
c
      if (first) go to 20
c
c     %-----------------------------------------------%
c     | Back from computing B*(orthogonalized-vector) |
c     %-----------------------------------------------%
c
      if (orth)  go to 40
c
      if (bmat .eq. 'G') then
         call cpu_time (t3)
         tmvopx = tmvopx + (t3 - t2)
      end if
c
c     %------------------------------------------------------%
c     | Starting vector is now in the range of OP; r = OP*r; |
c     | Compute B-norm of starting vector.                   |
c     %------------------------------------------------------%
c
      call cpu_time (t2)
      first = .TRUE.
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call dcopy (n, workd(n+1), 1, resid, 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call dcopy (n, resid, 1, workd, 1)
      end if
c
   20 continue
c
      if (bmat .eq. 'G') then
         call cpu_time (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
c
      first = .FALSE.
      if (bmat .eq. 'G') then
          rnorm0 = ddot (n, resid, 1, workd, 1)
          rnorm0 = sqrt(abs(rnorm0))
      else if (bmat .eq. 'I') then
           rnorm0 = dnrm2(n, resid, 1)
      end if
      rnorm  = rnorm0
c
c     %---------------------------------------------%
c     | Exit if this is the very first Arnoldi step |
c     %---------------------------------------------%
c
      if (j .eq. 1) go to 50
c
c     %----------------------------------------------------------------
c     | Otherwise need to B-orthogonalize the starting vector against |
c     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
c     | This is the case where an invariant subspace is encountered   |
c     | in the middle of the Arnoldi factorization.                   |
c     |                                                               |
c     |       s = V^{T}*B*r;   r = r - V*s;                           |
c     |                                                               |
c     | Stopping criteria used for iter. ref. is discussed in         |
c     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
c     %---------------------------------------------------------------%
c
      orth = .TRUE.
   30 continue
c
      call dgemv ('T', n, j-1, one, v, ldv, workd, 1,
     &            zero, workd(n+1), 1)
      call dgemv ('N', n, j-1, -one, v, ldv, workd(n+1), 1,
     &            one, resid, 1)
c
c     %----------------------------------------------------------%
c     | Compute the B-norm of the orthogonalized starting vector |
c     %----------------------------------------------------------%
c
      call cpu_time (t2)
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call dcopy (n, resid, 1, workd(n+1), 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call dcopy (n, resid, 1, workd, 1)
      end if
c
   40 continue
c
      if (bmat .eq. 'G') then
         call cpu_time (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
c
      if (bmat .eq. 'G') then
         rnorm = ddot (n, resid, 1, workd, 1)
         rnorm = sqrt(abs(rnorm))
      else if (bmat .eq. 'I') then
         rnorm = dnrm2(n, resid, 1)
      end if
c
c     %--------------------------------------%
c     | Check for further orthogonalization. |
c     %--------------------------------------%
c
      if (msglvl .gt. 2) then
          rnorm0_array(1) = rnorm0
          rnorm_array(1)  = rnorm
          call dvout (logfil, 1, RNORM0_ARRAY, ndigit,
     &                '_getv0: re-orthonalization ; rnorm0 is')
          call dvout (logfil, 1, RNORM_ARRAY, ndigit,
     &                '_getv0: re-orthonalization ; rnorm is')
      end if
c
      if (rnorm .gt. 0.717*rnorm0) go to 50
c
      iter = iter + 1
      if (iter .le. 5) then
c
c        %-----------------------------------%
c        | Perform iterative refinement step |
c        %-----------------------------------%
c
         rnorm0 = rnorm
         go to 30
      else
c
c        %------------------------------------%
c        | Iterative refinement step "failed" |
c        %------------------------------------%
c
         do 45 jj = 1, n
            resid(jj) = zero
   45    continue
         rnorm = zero
         ierr = -1
      end if
c
   50 continue
c
      if (msglvl .gt. 0) then
         rnorm_array(1) = rnorm
         call dvout (logfil, 1, RNORM_ARRAY, ndigit,
     &        '_getv0: B-norm of initial / restarted starting vector')
      end if
      if (msglvl .gt. 3) then
         call dvout (logfil, n, resid, ndigit,
     &        '_getv0: initial / restarted starting vector')
      end if
      ido = 99
c
      call cpu_time (t1)
      tgetv0 = tgetv0 + (t1 - t0)
c
 9000 continue
c

      RETURN

! **********************************************************************************************************************************

c     %---------------%
c     | End of dgetv0 |
c     %---------------%
c
      end subroutine dgetv0

! ##################################################################################################################################
! 010 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsaitr
c
c\Description:
c  Reverse communication interface for applying NP additional steps to
c  a K step symmetric Arnoldi factorization.
c
c  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
c
c          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
c
c  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
c
c          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
c
c  where OP and B are as in dsaupd.  The B-norm of r_{k+p} is also
c  computed and returned.
c
c\Usage:
c  call dsaitr
c     ( IDO, BMAT, N, K, NP, MODE, RESID, RNORM, V, LDV, H, LDH,
c       IPNTR, WORKD, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c                    This is for the restart phase to force the new
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y,
c                    IPNTR(3) is the pointer into WORK for B * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c          When the routine is used in the "shift-and-invert" mode, the
c          vector B * Q is already available and does not need to be
c          recomputed in forming OP * Q.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of matrix B that defines the
c          semi-inner product for the operator OP.  See dsaupd.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  K       Integer.  (INPUT)
c          Current order of H and the number of columns of V.
c
c  NP      Integer.  (INPUT)
c          Number of additional Arnoldi steps to take.
c
c  MODE    Integer.  (INPUT)
c          Signifies which form for "OP". If MODE=2 then
c          a reduction in the number of B matrix vector multiplies
c          is possible since the B-norm of OP*x is equivalent to
c          the inv(B)-norm of A*x.
c
c  RESID   REAL(DOUBLE) array of length N.  (INPUT/OUTPUT)
c          On INPUT:  RESID contains the residual vector r_{k}.
c          On OUTPUT: RESID contains the residual vector r_{k+p}.
c
c  RNORM   REAL(DOUBLE) scalar.  (INPUT/OUTPUT)
c          On INPUT the B-norm of r_{k}.
c          On OUTPUT the B-norm of the updated residual r_{k+p}.
c
c  V       REAL(DOUBLE) N by K+NP array.  (INPUT/OUTPUT)
c          On INPUT:  V contains the Arnoldi vectors in the first K
c          columns.
c          On OUTPUT: V contains the new NP Arnoldi vectors in the next
c          NP columns.  The first K columns are unchanged.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       REAL(DOUBLE) (K+NP) by 2 array.  (INPUT/OUTPUT)
c          H is used to store the generated symmetric tridiagonal matrix
c          with the subdiagonal in the first column starting at H(2,1)
c          and the main diagonal in the second column.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORK for
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c
c  WORKD   REAL(DOUBLE) work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The calling program should not
c          On INPUT, WORKD(1:N) = B*RESID where RESID is associated
c          with the K step Arnoldi factorization. Used to save some
c          computation at the first step.
c          On OUTPUT, WORKD(1:N) = B*RESID where RESID is associated
c          with the K+NP step Arnoldi factorization.
c
c  INFO    Integer.  (OUTPUT)
c          = 0: Normal exit.
c          > 0: Size of an invariant subspace of OP is found that is
c               less than K + NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dgetv0  ARPACK routine to generate the initial vector.
c     ivout   ARPACK utility routine that prints integers.
c     dmout   ARPACK utility routine that prints matrices.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dscal   Level 1 BLAS that scales a vector.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/93: Version ' 2.4'
c
c\SCCS Information: @(#)
c FILE: saitr.F   SID: 2.6   DATE OF SID: 8/28/96   RELEASE: 2
c
c\Remarks
c  The algorithm implemented is:
c
c  restart = .false.
c  Given V_{k} = [v_{1}, ..., v_{k}], r_{k};
c  r_{k} contains the initial residual vector even for k = 0;
c  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already
c  computed by the calling program.
c
c  betaj = rnorm ; p_{k+1} = B*r_{k} ;
c  For  j = k+1, ..., k+np  Do
c     1) if ( betaj < tol ) stop or restart depending on j.
c        if ( restart ) generate a new starting vector.
c     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];
c        p_{j} = p_{j}/betaj
c     3) r_{j} = OP*v_{j} where OP is defined as in dsaupd
c        For shift-invert mode p_{j} = B*v_{j} is already available.
c        wnorm = || OP*v_{j} ||
c     4) Compute the j-th step residual vector.
c        w_{j} =  V_{j}^T * B * OP * v_{j}
c        r_{j} =  OP*v_{j} - V_{j} * w_{j}
c        alphaj <- j-th component of w_{j}
c        rnorm = || r_{j} ||
c        betaj+1 = rnorm
c        If (rnorm > 0.717*wnorm) accept step and go back to 1)
c     5) Re-orthogonalization step:
c        s = V_{j}'*B*r_{j}
c        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
c        alphaj = alphaj + s_{j};
c     6) Iterative refinement step:
c        If (rnorm1 > 0.717*rnorm) then
c           rnorm = rnorm1
c           accept step and go back to 1)
c        Else
c           rnorm = rnorm1
c           If this is the first time in step 6), go to 5)
c           Else r_{j} lies in the span of V_{j} numerically.
c              Set r_{j} = 0 and rnorm = 0; go to 1)
c        EndIf
c  End Do
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsaitr
     &   (ido, bmat, n, k, np, mode, resid, rnorm, v, ldv, h, ldh,
     &    ipntr, workd, info)

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSAITR'
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cxxxx include   'debug.h'
cxxxx include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, mode, np
      REAL(DOUBLE)
     &           rnorm
      REAL(DOUBLE) RNORM_ARRAY(1)
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      REAL(DOUBLE)
     &           h(ldh,2), resid(n), v(ldv,k+np), workd(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      REAL(DOUBLE)
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    first, orth1, orth2, rstart, step3, step4
      integer    i, ierr, ipj, irj, ivj, iter, itry, j, msglvl,
     &           infol, jj
      INTEGER    J_ARRAY(1)
      REAL(DOUBLE)
     &           rnorm1, wnorm, safmin, temp1
      save       orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j, msglvl,
     &           rnorm1, safmin, wnorm
c
c     %-----------------------%
c     | Local Array Arguments |
c     %-----------------------%
c
      REAL(DOUBLE)
     &           xtemp(2)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
!:!   external   daxpy, dcopy, dscal, dgemv, dgetv0, dvout, dmout,
!:!  &           dlascl, ivout, second
      external   dgemv, dscal
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      REAL(DOUBLE) dlamch
      external     dlamch

      REAL(DOUBLE)
     &           ddot, dnrm2
      external   ddot, dnrm2
c
c     %-----------------%
c     | Data statements |
c     %-----------------%
c
      data      first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c


! **********************************************************************************************************************************

      if (first) then
         first = .false.
c
c        %--------------------------------%
c        | safmin = safe minimum is such  |
c        | that 1/sfmin does not overflow |
c        %--------------------------------%
c
         safmin = dlamch('safmin')
      end if
c
      if (ido .eq. 0) then
c
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call cpu_time (t0)
         msglvl = msaitr
c
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
c
c        %--------------------------------%
c        | Pointer to the current step of |
c        | the factorization to build     |
c        %--------------------------------%
c
         j      = k + 1
c
c        %------------------------------------------%
c        | Pointers used for reverse communication  |
c        | when using WORKD.                        |
c        %------------------------------------------%
c
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
c
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true.                                  |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         dgetv0.                                 |
c     %-------------------------------------------------%
c
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
c
c     %------------------------------%
c     | Else this is the first step. |
c     %------------------------------%
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%
c
 1000 continue
c
         if (msglvl .gt. 2) then
            j_array(1)     = j
            rnorm_array(1) = rnorm
            call ivout (logfil, 1, J_ARRAY, ndigit,
     &                  '_saitr: generating Arnoldi vector no.')
            call dvout (logfil, 1, RNORM_ARRAY, ndigit,
     &                  '_saitr: B-norm of the current residual =')
         end if
c
c        %---------------------------------------------------------%
c        | Check for exact zero. Equivalent to determing whether a |
c        | j-step Arnoldi factorization is present.                |
c        %---------------------------------------------------------%
c
         if (rnorm .gt. zero) go to 40
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c
            if (msglvl .gt. 0) then
               j_array(1) = j
               call ivout (logfil, 1, J_ARRAY, ndigit,
     &                     '_saitr: ****** restart at step ******')
            end if
c
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted. NRSTRT is used by stat.h         |
c           %---------------------------------------------%
c
            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c

            call dgetv0 (ido, bmat, itry, .false., n, j, v, ldv,
     &                   resid, rnorm, ipntr, workd, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
               info = j - 1
               call cpu_time (t1)
               tsaitr = tsaitr + (t1 - t0)
               ido = 99
               go to 9000
            end if
c
   40    continue
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
         call dcopy (n, resid, 1, v(1,j), 1)
         if (rnorm .ge. safmin) then
             temp1 = one / rnorm
             call dscal (n, temp1, v(1,j), 1)
             call dscal (n, temp1, workd(ipj), 1)
         else
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine SLASCL               |
c            %-----------------------------------------%
c
             call dlascl ('General', i, i, rnorm, one, n, 1,
     &                    v(1,j), n, infol)
             call dlascl ('General', i, i, rnorm, one, n, 1,
     &                    workd(ipj), n, infol)
         end if
c
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
         step3 = .true.
         nopx  = nopx + 1
         call cpu_time (t2)
         call dcopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
c
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c
         go to 9000
   50    continue
c
c        %-----------------------------------%
c        | Back from reverse communication;  |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}.   |
c        %-----------------------------------%
c
         call cpu_time (t3)
         tmvopx = tmvopx + (t3 - t2)
c
         step3 = .false.
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
         call dcopy (n, workd(irj), 1, resid, 1)
c
c        %-------------------------------------------%
c        | STEP 4:  Finish extending the symmetric   |
c        |          Arnoldi to length j. If MODE = 2 |
c        |          then B*OP = B*inv(B)*A = A and   |
c        |          we don't need to compute B*OP.   |
c        | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |
c        | assumed to have A*v_{j}.                  |
c        %-------------------------------------------%
c
         if (mode .eq. 2) go to 65
         call cpu_time (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
              call dcopy(n, resid, 1 , workd(ipj), 1)
         end if
   60    continue
c
c        %-----------------------------------%
c        | Back from reverse communication;  |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}. |
c        %-----------------------------------%
c
         if (bmat .eq. 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c
         step4 = .false.
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
   65    continue
         if (mode .eq. 2) then
c
c           %----------------------------------%
c           | Note that the B-norm of OP*v_{j} |
c           | is the inv(B)-norm of A*v_{j}.   |
c           %----------------------------------%
c
            wnorm = ddot (n, resid, 1, workd(ivj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'G') then
            wnorm = ddot (n, resid, 1, workd(ipj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'I') then
            wnorm = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c
         if (mode .ne. 2 ) then
            call dgemv('T', n, j, one, v, ldv, workd(ipj), 1, zero,
     &                  workd(irj), 1)
         else if (mode .eq. 2) then
            call dgemv('T', n, j, one, v, ldv, workd(ivj), 1, zero,
     &                  workd(irj), 1)
         end if
c
c        %--------------------------------------%
c        | Orthgonalize r_{j} against V_{j}.    |
c        | RESID contains OP*v_{j}. See STEP 3. |
c        %--------------------------------------%
c
         call dgemv('N', n, j, -one, v, ldv, workd(irj), 1, one,
     &               resid, 1)
c
c        %--------------------------------------%
c        | Extend H to have j rows and columns. |
c        %--------------------------------------%
c
         h(j,2) = workd(irj + j - 1)
         if (j .eq. 1  .or.  rstart) then
            h(j,1) = zero
         else
            h(j,1) = rnorm
         end if
         call cpu_time (t4)
c
         orth1 = .true.
         iter  = 0
c
         call cpu_time (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   70    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
         if (bmat .eq. 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c
         orth1 = .false.
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
         if (bmat .eq. 'G') then
            rnorm = ddot (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        %-----------------------------------------------------------%
c
         if (rnorm .gt. 0.717*wnorm) go to 100
         nrorth = nrorth + 1
c
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c
   80    continue
c
         if (msglvl .gt. 2) then
            xtemp(1) = wnorm
            xtemp(2) = rnorm
            call dvout (logfil, 2, xtemp, ndigit,
     &           '_saitr: re-orthonalization ; wnorm and rnorm are')
         end if
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
         call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, workd(irj), 1)
c
c        %----------------------------------------------%
c        | Compute the correction to the residual:      |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).  |
c        | The correction to H is v(:,1:J)*H(1:J,1:J) + |
c        | v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j, but only   |
c        | H(j,j) is updated.                           |
c        %----------------------------------------------%
c
         call dgemv ('N', n, j, -one, v, ldv, workd(irj), 1,
     &               one, resid, 1)
c
         if (j .eq. 1  .or.  rstart) h(j,1) = zero
         h(j,2) = h(j,2) + workd(irj + j - 1)
c
         orth2 = .true.
         call cpu_time (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   90    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
         if (bmat .eq. 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c
         if (bmat .eq. 'G') then
             rnorm1 = ddot (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat .eq. 'I') then
             rnorm1 = dnrm2(n, resid, 1)
         end if
c
         if (msglvl .gt. 0 .and. iter .gt. 0) then
            j_array(1) = j
            call ivout (logfil, 1, J_ARRAY, ndigit,
     &           '_saitr: Iterative refinement for Arnoldi residual')
            if (msglvl .gt. 2) then
                xtemp(1) = rnorm
                xtemp(2) = rnorm1
                call dvout (logfil, 2, xtemp, ndigit,
     &           '_saitr: iterative refinement ; rnorm and rnorm1 are')
            end if
         end if
c
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
         if (rnorm1 .gt. 0.717*rnorm) then
c
c           %--------------------------------%
c           | No need for further refinement |
c           %--------------------------------%
c
            rnorm = rnorm1
c
         else
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required. NITREF is used by stat.h     |
c           %-------------------------------------------%
c
            nitref = nitref + 1
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
c
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c
  100    continue
c
         rstart = .false.
         orth2  = .false.
c
         call cpu_time (t5)
         titref = titref + (t5 - t4)
c
c        %----------------------------------------------------------%
c        | Make sure the last off-diagonal element is non negative  |
c        | If not perform a similarity transformation on H(1:j,1:j) |
c        | and scale v(:,j) by -1.                                  |
c        %----------------------------------------------------------%
c
         if (h(j,1) .lt. zero) then
            h(j,1) = -h(j,1)
            if ( j .lt. k+np) then
               call dscal(n, -one, v(1,j+1), 1)
            else
               call dscal(n, -one, resid, 1)
            end if
         end if
c
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
         j = j + 1
         if (j .gt. k+np) then
            call cpu_time (t1)
            tsaitr = tsaitr + (t1 - t0)
            ido = 99
c
            if (msglvl .gt. 1) then
               call dvout (logfil, k+np, h(1,2), ndigit,
     &         '_saitr: main diagonal of matrix H of step K+NP.')
               if (k+np .gt. 1) then
               call dvout (logfil, k+np-1, h(2,1), ndigit,
     &         '_saitr: sub diagonal of matrix H of step K+NP.')
               end if
            end if
c
            go to 9000
         end if
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 9000 continue
c

      RETURN

! **********************************************************************************************************************************

c     %---------------%
c     | End of dsaitr |
c     %---------------%
c
      end subroutine dsaitr

! ##################################################################################################################################
! 011 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP shifts implicitly resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix of order KEV+NP. Q is the product of
c  rotations resulting from the NP bulge chasing sweeps.  The updated Arnoldi
c  factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call dsapps
c     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. dimension of matrix A.
c
c  KEV     Integer.  (INPUT)
c          INPUT: KEV+NP is the size of the input matrix H.
c          OUTPUT: KEV is the size of the updated matrix HNEW.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFT   REAL(DOUBLE) array of length NP.  (INPUT)
c          The shifts to be applied.
c
c  V       REAL(DOUBLE) N by (KEV+NP) array.  (INPUT/OUTPUT)
c          INPUT: V contains the current KEV+NP Arnoldi vectors.
c          OUTPUT: VNEW = V(1:n,1:KEV); the updated Arnoldi vectors
c          are in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       REAL(DOUBLE) (KEV+NP) by 2 array.  (INPUT/OUTPUT)
c          INPUT: H contains the symmetric tridiagonal matrix of the
c          Arnoldi factorization with the subdiagonal in the 1st column
c          starting at H(2,1) and the main diagonal in the 2nd column.
c          OUTPUT: H contains the updated tridiagonal matrix in the
c          KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   REAL(DOUBLE) array of length (N).  (INPUT/OUTPUT)
c          INPUT: RESID contains the the residual vector r_{k+p}.
c          OUTPUT: RESID is the updated residual vector rnew_{k}.
c
c  Q       REAL(DOUBLE) KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations during the bulge
c          chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKD   REAL(DOUBLE) work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlacpy  LAPACK matrix copy routine.
c     dlaset  LAPACK matrix initialization routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dscal   Level 1 BLAS that scales a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     12/16/93: Version ' 2.4'
c
c\SCCS Information: @(#)
c FILE: sapps.F   SID: 2.6   DATE OF SID: 3/28/97   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the subblocks of
c     the tridiagonal matrix H and not just to the submatrix that it
c     comes from. This routine assumes that the subdiagonal elements
c     of H that are stored in h(1:kev+np,1) are nonegative upon input
c     and enforce this condition upon output. This version incorporates
c     deflation. See code for documentation.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsapps
     &   ( n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, workd )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSAPPS'
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cxxxx include   'debug.h'
cxxxx include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    kev, ldh, ldq, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      REAL(DOUBLE)
     &           h(ldh,2), q(ldq,kev+np), resid(n), shift(np),
     &           v(ldv,kev+np), workd(2*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      REAL(DOUBLE)
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, iend, istart, itop, j, jj, kplusp, msglvl
      INTEGER    I_ARRAY(1), JJ_ARRAY(1)
      logical    first
      REAL(DOUBLE)
     &           a1, a2, a3, a4, big, c, epsmch, f, g, r, s
      save       epsmch, first
c
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
!:!   external   daxpy, dcopy, dscal, dlacpy, dlartg, dlaset, dvout,
!:!  &           ivout, second, dgemv
      external   dgemv, dscal
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      REAL(DOUBLE) dlamch
      external     dlamch
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs
c
c     %----------------%
c     | Data statments |
c     %----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c


! **********************************************************************************************************************************

      if (first) then
         epsmch = dlamch('Epsilon-Machine')
         first = .false.
      end if
      itop = 1
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call cpu_time (t0)
      msglvl = msapps
c
      kplusp = kev + np
c
c     %----------------------------------------------%
c     | Initialize Q to the identity matrix of order |
c     | kplusp used to accumulate the rotations.     |
c     %----------------------------------------------%
c
      call dlaset ('All', kplusp, kplusp, zero, one, q, ldq)
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
c     %----------------------------------------------------------%
c     | Apply the np shifts implicitly. Apply each shift to the  |
c     | whole matrix and not just to the submatrix from which it |
c     | comes.                                                   |
c     %----------------------------------------------------------%
c
      do 90 jj = 1, np
c
         istart = itop
c
c        %----------------------------------------------------------%
c        | Check for splitting and deflation. Currently we consider |
c        | an off-diagonal element h(i+1,1) negligible if           |
c        |         h(i+1,1) .le. epsmch*( |h(i,2)| + |h(i+1,2)| )   |
c        | for i=1:KEV+NP-1.                                        |
c        | If above condition tests true then we set h(i+1,1) = 0.  |
c        | Note that h(1:KEV+NP,1) are assumed to be non negative.  |
c        %----------------------------------------------------------%
c
   20    continue
c
c        %------------------------------------------------%
c        | The following loop exits early if we encounter |
c        | a negligible off diagonal element.             |
c        %------------------------------------------------%
c
         do 30 i = istart, kplusp-1
            big   = abs(h(i,2)) + abs(h(i+1,2))
            if (h(i+1,1) .le. epsmch*big) then
               if (msglvl .gt. 0) then
                  i_array(1)  = i
                  jj_array(1) = jj
                  call ivout (logfil, 1, I_ARRAY, ndigit,
     &                 '_sapps: deflation at row/column no.')
                  call ivout (logfil, 1, JJ_ARRAY, ndigit,
     &                 '_sapps: occured before shift number.')
                  call dvout (logfil, 1, h(i+1,1), ndigit,
     &                 '_sapps: the corresponding off diagonal element')
               end if
               h(i+1,1) = zero
               iend = i
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
c
         if (istart .lt. iend) then
c
c           %--------------------------------------------------------%
c           | Construct the plane rotation G'(istart,istart+1,theta) |
c           | that attempts to drive h(istart+1,1) to zero.          |
c           %--------------------------------------------------------%
c
             f = h(istart,2) - shift(jj)
             g = h(istart+1,1)
             call dlartg (f, g, c, s, r)
c
c            %-------------------------------------------------------%
c            | Apply rotation to the left and right of H;            |
c            | H <- G' * H * G,  where G = G(istart,istart+1,theta). |
c            | This will create a "bulge".                           |
c            %-------------------------------------------------------%
c
             a1 = c*h(istart,2)   + s*h(istart+1,1)
             a2 = c*h(istart+1,1) + s*h(istart+1,2)
             a4 = c*h(istart+1,2) - s*h(istart+1,1)
             a3 = c*h(istart+1,1) - s*h(istart,2)
             h(istart,2)   = c*a1 + s*a2
             h(istart+1,2) = c*a4 - s*a3
             h(istart+1,1) = c*a3 + s*a4
c
c            %----------------------------------------------------%
c            | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c            %----------------------------------------------------%
c
             do 60 j = 1, min(istart+jj,kplusp)
                a1            =   c*q(j,istart) + s*q(j,istart+1)
                q(j,istart+1) = - s*q(j,istart) + c*q(j,istart+1)
                q(j,istart)   = a1
   60        continue
c
c
c            %----------------------------------------------%
c            | The following loop chases the bulge created. |
c            | Note that the previous rotation may also be  |
c            | done within the following loop. But it is    |
c            | kept separate to make the distinction among  |
c            | the bulge chasing sweeps and the first plane |
c            | rotation designed to drive h(istart+1,1) to  |
c            | zero.                                        |
c            %----------------------------------------------%
c
             do 70 i = istart+1, iend-1
c
c               %----------------------------------------------%
c               | Construct the plane rotation G'(i,i+1,theta) |
c               | that zeros the i-th bulge that was created   |
c               | by G(i-1,i,theta). g represents the bulge.   |
c               %----------------------------------------------%
c
                f = h(i,1)
                g = s*h(i+1,1)
c
c               %----------------------------------%
c               | Final update with G(i-1,i,theta) |
c               %----------------------------------%
c
                h(i+1,1) = c*h(i+1,1)
                call dlartg (f, g, c, s, r)
c
c               %-------------------------------------------%
c               | The following ensures that h(1:iend-1,1), |
c               | the first iend-2 off diagonal of elements |
c               | H, remain non negative.                   |
c               %-------------------------------------------%
c
                if (r .lt. zero) then
                   r = -r
                   c = -c
                   s = -s
                end if
c
c               %--------------------------------------------%
c               | Apply rotation to the left and right of H; |
c               | H <- G * H * G',  where G = G(i,i+1,theta) |
c               %--------------------------------------------%
c
                h(i,1) = r
c
                a1 = c*h(i,2)   + s*h(i+1,1)
                a2 = c*h(i+1,1) + s*h(i+1,2)
                a3 = c*h(i+1,1) - s*h(i,2)
                a4 = c*h(i+1,2) - s*h(i+1,1)
c
                h(i,2)   = c*a1 + s*a2
                h(i+1,2) = c*a4 - s*a3
                h(i+1,1) = c*a3 + s*a4
c
c               %----------------------------------------------------%
c               | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c               %----------------------------------------------------%
c
                do 50 j = 1, min( i+jj, kplusp )
                   a1       =   c*q(j,i) + s*q(j,i+1)
                   q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                   q(j,i)   = a1
   50           continue
c
   70        continue
c
         end if
c
c        %--------------------------%
c        | Update the block pointer |
c        %--------------------------%
c
         istart = iend + 1
c
c        %------------------------------------------%
c        | Make sure that h(iend,1) is non-negative |
c        | If not then set h(iend,1) <-- -h(iend,1) |
c        | and negate the last column of Q.         |
c        | We have effectively carried out a        |
c        | similarity on transformation H           |
c        %------------------------------------------%
c
         if (h(iend,1) .lt. zero) then
             h(iend,1) = -h(iend,1)
             call dscal(kplusp, -one, q(1,iend), 1)
         end if
c
c        %--------------------------------------------------------%
c        | Apply the same shift to the next block if there is any |
c        %--------------------------------------------------------%
c
         if (iend .lt. kplusp) go to 20
c
c        %-----------------------------------------------------%
c        | Check if we can increase the the start of the block |
c        %-----------------------------------------------------%
c
         do 80 i = itop, kplusp-1
            if (h(i+1,1) .gt. zero) go to 90
            itop  = itop + 1
   80    continue
c
c        %-----------------------------------%
c        | Finished applying the jj-th shift |
c        %-----------------------------------%
c
   90 continue
c
c     %------------------------------------------%
c     | All shifts have been applied. Check for  |
c     | more possible deflation that might occur |
c     | after the last shift is applied.         |
c     %------------------------------------------%
c
      do 100 i = itop, kplusp-1
         big   = abs(h(i,2)) + abs(h(i+1,2))
         if (h(i+1,1) .le. epsmch*big) then
            if (msglvl .gt. 0) then
               i_array(1) = i
               call ivout (logfil, 1, I_ARRAY, ndigit,
     &              '_sapps: deflation at row/column no.')
               call dvout (logfil, 1, h(i+1,1), ndigit,
     &              '_sapps: the corresponding off diagonal element')
            end if
            h(i+1,1) = zero
         end if
 100  continue
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is not necessary if h(kev+1,1) = 0.         |
c     %-------------------------------------------------%
c
      if ( h(kev+1,1) .gt. zero )
     &   call dgemv ('N', n, kplusp, one, v, ldv,
     &                q(1,kev+1), 1, zero, workd(n+1), 1)
c
c     %-------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order    |
c     | taking advantage that Q is an upper triangular matrix |
c     | with lower bandwidth np.                              |
c     | Place results in v(:,kplusp-kev:kplusp) temporarily.  |
c     %-------------------------------------------------------%
c
      do 130 i = 1, kev
         call dgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call dcopy (n, workd, 1, v(1,kplusp-i+1), 1)
  130 continue
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
      call dlacpy ('All', n, kev, v(1,np+1), ldv, v, ldv)
c
c     %--------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the |
c     | appropriate place if h(kev+1,1) .ne. zero. |
c     %--------------------------------------------%
c
      if ( h(kev+1,1) .gt. zero )
     &     call dcopy (n, workd(n+1), 1, v(1,kev+1), 1)
c
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
      call dscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,1) .gt. zero)
     &   call daxpy (n, h(kev+1,1), v(1,kev+1), 1, resid, 1)
c
      if (msglvl .gt. 1) then
         call dvout (logfil, 1, q(kplusp,kev), ndigit,
     &      '_sapps: sigmak of the updated residual vector')
         call dvout (logfil, 1, h(kev+1,1), ndigit,
     &      '_sapps: betak of the updated residual vector')
         call dvout (logfil, kev, h(1,2), ndigit,
     &      '_sapps: updated main diagonal of H for next iteration')
         if (kev .gt. 1) then
         call dvout (logfil, kev-1, h(2,1), ndigit,
     &      '_sapps: updated sub diagonal of H for next iteration')
         end if
      end if
c
      call cpu_time (t1)
      tsapps = tsapps + (t1 - t0)
c
 9000 continue
c

      RETURN

! **********************************************************************************************************************************

c     %---------------%
c     | End of dsapps |
c     %---------------%
c
      end subroutine dsapps

! ##################################################################################################################################
! 012 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsconv
c
c\Description:
c  Convergence testing for the symmetric Arnoldi eigenvalue routine.
c
c\Usage:
c  call dsconv
c     ( N, RITZ, BOUNDS, TOL, NCONV )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Number of Ritz values to check for convergence.
c
c  RITZ    REAL(DOUBLE) array of length N.  (INPUT)
c          The Ritz values to be checked for convergence.
c
c  BOUNDS  REAL(DOUBLE) array of length N.  (INPUT)
c          Ritz estimates associated with the Ritz values in RITZ.
c
c  TOL     REAL(DOUBLE) scalar.  (INPUT)
c          Desired relative accuracy for a Ritz value to be considered
c          "converged".
c
c  NCONV   Integer scalar.  (OUTPUT)
c          Number of "converged" Ritz values.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Routines called:
c     second  ARPACK utility routine for timing.
c     dlamch  LAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
c
c\Remarks
c     1. Starting with version 2.4, this routine no longer uses the
c        Parlett strategy using the gap conditions.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsconv (n, ritz, bounds, tol, nconv)

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSCONV'
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cxxxx include   'debug.h'
cxxxx include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    n, nconv
      REAL(DOUBLE)
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      REAL(DOUBLE)
     &           ritz(n), bounds(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i
      REAL(DOUBLE)
     &           temp, eps23
c
c     %-------------------%
c     | External routines |
c     %-------------------%
c
      REAL(DOUBLE) dlamch
      external     dlamch

c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c


! **********************************************************************************************************************************

      call cpu_time (t0)
c
      eps23 = dlamch('Epsilon-Machine')
      eps23 = eps23**(2.0D+0 / 3.0D+0)
c
      nconv  = 0
      do 10 i = 1, n
c
c        %-----------------------------------------------------%
c        | The i-th Ritz value is considered "converged"       |
c        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |
c        %-----------------------------------------------------%
c
         temp = max( eps23, abs(ritz(i)) )
         if ( bounds(i) .le. tol*temp ) then
            nconv = nconv + 1
         end if
c
   10 continue
c
      call cpu_time (t1)
      tsconv = tsconv + (t1 - t0)
c

      RETURN

! **********************************************************************************************************************************

c     %---------------%
c     | End of dsconv |
c     %---------------%
c
      end subroutine dsconv

! ##################################################################################################################################
! 013 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dseigt
c
c\Description:
c  Compute the eigenvalues of the current symmetric tridiagonal matrix
c  and the corresponding error bounds given the current residual norm.
c
c\Usage:
c  call dseigt
c     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )
c
c\Arguments
c  RNORM   REAL(DOUBLE) scalar.  (INPUT)
c          RNORM contains the residual norm corresponding to the current
c          symmetric tridiagonal matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the symmetric tridiagonal matrix H.
c
c  H       REAL(DOUBLE) N by 2 array.  (INPUT)
c          H contains the symmetric tridiagonal matrix with the
c          subdiagonal in the first column starting at H(2,1) and the
c          main diagonal in second column.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  EIG     REAL(DOUBLE) array of length N.  (OUTPUT)
c          On output, EIG contains the N eigenvalues of H possibly
c          unsorted.  The BOUNDS arrays are returned in the
c          same sorted order as EIG.
c
c  BOUNDS  REAL(DOUBLE) array of length N.  (OUTPUT)
c          On output, BOUNDS contains the error estimates corresponding
c          to the eigenvalues EIG.  This is equal to RNORM times the
c          last components of the eigenvectors corresponding to the
c          eigenvalues in EIG.
c
c  WORKL   REAL(DOUBLE) work array of length 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from dstqrb.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dstqrb  ARPACK routine that computes the eigenvalues and the
c             last components of the eigenvectors of a symmetric
c             and tridiagonal matrix.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dcopy   Level 1 BLAS that copies one vector to another.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#)
c FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dseigt
     &   ( rnorm, n, h, ldh, eig, bounds, workl, ierr )


      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSEIGT'
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
cxxxx include   'debug.h'
cxxxx include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    ierr, ldh, n
      REAL(DOUBLE)
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      REAL(DOUBLE)
     &           eig(n), bounds(n), h(ldh,2), workl(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      REAL(DOUBLE)
     &           zero
      parameter (zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, k, msglvl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
!:!   external   dcopy, dstqrb, dvout, second
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c


! **********************************************************************************************************************************

c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call cpu_time (t0)
      msglvl = mseigt
c
      if (msglvl .gt. 0) then
         call dvout (logfil, n, h(1,2), ndigit,
     &              '_seigt: main diagonal of matrix H')
         if (n .gt. 1) then
         call dvout (logfil, n-1, h(2,1), ndigit,
     &              '_seigt: sub diagonal of matrix H')
         end if
      end if
c
      call dcopy  (n, h(1,2), 1, eig, 1)
      call dcopy  (n-1, h(2,1), 1, workl, 1)
      call dstqrb (n, eig, workl, bounds, workl(n+1), ierr)
      if (ierr .ne. 0) go to 9000
      if (msglvl .gt. 1) then
         call dvout (logfil, n, bounds, ndigit,
     &              '_seigt: last row of the eigenvector matrix for H')
      end if
c
c     %-----------------------------------------------%
c     | Finally determine the error bounds associated |
c     | with the n Ritz values of H.                  |
c     %-----------------------------------------------%
c
      do 30 k = 1, n
         bounds(k) = rnorm*abs(bounds(k))
   30 continue
c
      call cpu_time (t1)
      tseigt = tseigt + (t1 - t0)
c
 9000 continue
c

      RETURN

! **********************************************************************************************************************************

c     %---------------%
c     | End of dseigt |
c     %---------------%
c
      end subroutine dseigt

! ##################################################################################################################################
! 014 ARPACK_LANCZOS_EIG

c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dstqrb
c
c\Description:
c  Computes all eigenvalues and the last component of the eigenvectors
c  of a symmetric tridiagonal matrix using the implicit QL or QR method.
c
c  This is mostly a modification of the LAPACK routine dsteqr.
c  See Remarks.
c
c\Usage:
c  call dstqrb
c     ( N, D, E, Z, WORK, INFO )
c
c\Arguments
c  N       Integer.  (INPUT)
c          The number of rows and columns in the matrix.  N >= 0.
c
c  D       REAL(DOUBLE) array, dimension (N).  (INPUT/OUTPUT)
c          On entry, D contains the diagonal elements of the
c          tridiagonal matrix.
c          On exit, D contains the eigenvalues, in ascending order.
c          If an error exit is made, the eigenvalues are correct
c          for indices 1,2,...,INFO-1, but they are unordered and
c          may not be the smallest eigenvalues of the matrix.
c
c  E       REAL(DOUBLE) array, dimension (N-1).  (INPUT/OUTPUT)
c          On entry, E contains the subdiagonal elements of the
c          tridiagonal matrix in positions 1 through N-1.
c          On exit, E has been destroyed.
c
c  Z       REAL(DOUBLE) array, dimension (N).  (OUTPUT)
c          On exit, Z contains the last row of the orthonormal
c          eigenvector matrix of the symmetric tridiagonal matrix.
c          If an error exit is made, Z contains the last row of the
c          eigenvector matrix associated with the stored eigenvalues.
c
c  WORK    REAL(DOUBLE) array, dimension (max(1,2*N-2)).  (WORKSPACE)
c          Workspace used in accumulating the transformation for
c          computing the last components of the eigenvectors.
c
c  INFO    Integer.  (OUTPUT)
c          = 0:  normal return.
c          < 0:  if INFO = -i, the i-th argument had an illegal value.
c          > 0:  if INFO = +i, the i-th eigenvalue has not converged
c                              after a total of  30*N  iterations.
c
c\Remarks
c  1. None.
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dswap   Level 1 BLAS that swaps the contents of two vectors.
c     lsame   LAPACK character comparison routine.
c     dlae2   LAPACK routine that computes the eigenvalues of a 2-by-2
c             symmetric matrix.
c     dlaev2  LAPACK routine that eigendecomposition of a 2-by-2 symmetric
c             matrix.
c     dlamch  LAPACK routine that determines machine constants.
c     dlanst  LAPACK routine that computes the norm of a matrix.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dlaset  LAPACK matrix initialization routine.
c     dlasr   LAPACK routine that applies an orthogonal transformation to
c             a matrix.
c     dlasrt  LAPACK sorting routine.
c     dsteqr  LAPACK routine that computes eigenvalues and eigenvectors
c             of a symmetric tridiagonal matrix.
c     xerbla  LAPACK error handler routine.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: stqrb.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     1. Starting with version 2.5, this routine is a modified version
c        of LAPACK version 2.0 subroutine SSTEQR. No lines are deleted,
c        only commeted out and new lines inserted.
c        All lines commented out have "c$$$" at the beginning.
c        Note that the LAPACK version 1.0 subroutine SSTEQR contained
c        bugs.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dstqrb ( n, d, e, z, work, info )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))  :: SUBR_NAME = 'DSTQRB'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    info, n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      REAL(DOUBLE)
     &           d( n ), e( n-1 ), z( n ), work( 2*n-2 )
c
c     .. parameters ..
      REAL(DOUBLE)
     &                   zero, one, two, three
      parameter          ( zero = 0.0D+0, one = 1.0D+0,
     &                     two = 2.0D+0, three = 3.0D+0 )
      integer            maxit
      parameter          ( maxit = 30 )
c     ..
c     .. local scalars ..
      integer            i, icompz, ii, iscale, j, jtot, k, l, l1, lend,
     &                   lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1,
     &                   nm1, nmaxit
      REAL(DOUBLE)
     &                   anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2,
     &                   s, safmax, safmin, ssfmax, ssfmin, tst
c     ..
c     .. external functions ..
      logical            lsame
      REAL(DOUBLE)
     &                   dlapy2
      external           lsame, dlapy2

      REAL(DOUBLE)       dlamch, dlanst
      external           dlamch, dlanst
c     ..
c     .. external subroutines ..
!:!   external           dlae2, dlaev2, dlartg, dlascl, dlaset, dlasr,
!:!  &                   dlasrt, dswap, xerbla
      external   dswap
c     ..
c     .. intrinsic functions ..
      intrinsic          abs, max, sign, sqrt
c     ..
c     .. executable statements ..
c


! **********************************************************************************************************************************

c     test the input parameters.
c
      info = 0
c
c$$$      IF( LSAME( COMPZ, 'N' ) ) THEN
c$$$         ICOMPZ = 0
c$$$      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
c$$$         ICOMPZ = 1
c$$$      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
c$$$         ICOMPZ = 2
c$$$      ELSE
c$$$         ICOMPZ = -1
c$$$      END IF
c$$$      IF( ICOMPZ.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
c$$$     $         N ) ) ) THEN
c$$$         INFO = -6
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'SSTEQR', -INFO )
c$$$         RETURN
c$$$      END IF
c
c    *** New starting with version 2.5 ***
c
      icompz = 2
c    *************************************
c
c     quick return if possible
c
      if( n.eq.0 )
     $   return
c
      if( n.eq.1 ) then
         if( icompz.eq.2 )  z( 1 ) = one
         return
      end if
c
c     determine the unit roundoff and over/underflow thresholds.
c
      eps = dlamch( 'e' )
      eps2 = eps**2
      safmin = dlamch( 's' )
      safmax = one / safmin
      ssfmax = sqrt( safmax ) / three
      ssfmin = sqrt( safmin ) / eps2
c
c     compute the eigenvalues and eigenvectors of the tridiagonal
c     matrix.
c
c$$      if( icompz.eq.2 )
c$$$     $   call dlaset( 'full', n, n, zero, one, z, ldz )
c
c     *** New starting with version 2.5 ***
c
      if ( icompz .eq. 2 ) then
         do 5 j = 1, n-1
            z(j) = zero
  5      continue
         z( n ) = one
      end if
c     *************************************
c
      nmaxit = n*maxit
      jtot = 0
c
c     determine where the matrix splits and choose ql or qr iteration
c     for each block, according to whether top or bottom diagonal
c     element is smaller.
c
      l1 = 1
      nm1 = n - 1
c
   10 continue
      if( l1.gt.n )
     $   go to 160
      if( l1.gt.1 )
     $   e( l1-1 ) = zero
      if( l1.le.nm1 ) then
         do 20 m = l1, nm1
            tst = abs( e( m ) )
            if( tst.eq.zero )
     $         go to 30
            if( tst.le.( sqrt( abs( d( m ) ) )*sqrt( abs( d( m+
     $          1 ) ) ) )*eps ) then
               e( m ) = zero
               go to 30
            end if
   20    continue
      end if
      m = n
c
   30 continue
      l = l1
      lsv = l
      lend = m
      lendsv = lend
      l1 = m + 1
      if( lend.eq.l )
     $   go to 10
c
c     scale submatrix in rows and columns l to lend
c
      anorm = dlanst( 'i', lend-l+1, d( l ), e( l ) )
      iscale = 0
      if( anorm.eq.zero )
     $   go to 10
      if( anorm.gt.ssfmax ) then
         iscale = 1
         call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l ), n,
     $                info )
         call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l, 1, e( l ), n,
     $                info )
      else if( anorm.lt.ssfmin ) then
         iscale = 2
         call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l ), n,
     $                info )
         call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l, 1, e( l ), n,
     $                info )
      end if
c
c     choose between ql and qr iteration
c
      if( abs( d( lend ) ).lt.abs( d( l ) ) ) then
         lend = lsv
         l = lendsv
      end if
c
      if( lend.gt.l ) then
c
c        ql iteration
c
c        look for small subdiagonal element.
c
   40    continue
         if( l.ne.lend ) then
            lendm1 = lend - 1
            do 50 m = l, lendm1
               tst = abs( e( m ) )**2
               if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m+1 ) )+
     $             safmin )go to 60
   50       continue
         end if
c
         m = lend
c
   60    continue
         if( m.lt.lend )
     $      e( m ) = zero
         p = d( l )
         if( m.eq.l )
     $      go to 80
c
c        if remaining matrix is 2-by-2, use dlae2 or dlaev2
c        to compute its eigensystem.
c
         if( m.eq.l+1 ) then
            if( icompz.gt.0 ) then
               call dlaev2( d( l ), e( l ), d( l+1 ), rt1, rt2, c, s )
               work( l ) = c
               work( n-1+l ) = s
c$$$               call dlasr( 'r', 'v', 'b', n, 2, work( l ),
c$$$     $                     work( n-1+l ), z( 1, l ), ldz )
c
c              *** New starting with version 2.5 ***
c
               tst      = z(l+1)
               z(l+1) = c*tst - s*z(l)
               z(l)   = s*tst + c*z(l)
c              *************************************
            else
               call dlae2( d( l ), e( l ), d( l+1 ), rt1, rt2 )
            end if
            d( l ) = rt1
            d( l+1 ) = rt2
            e( l ) = zero
            l = l + 2
            if( l.le.lend )
     $         go to 40
            go to 140
         end if
c
         if( jtot.eq.nmaxit )
     $      go to 140
         jtot = jtot + 1
c
c        form shift.
c
         g = ( d( l+1 )-p ) / ( two*e( l ) )
         r = dlapy2( g, one )
         g = d( m ) - p + ( e( l ) / ( g+sign( r, g ) ) )
c
         s = one
         c = one
         p = zero
c
c        inner loop
c
         mm1 = m - 1
         do 70 i = mm1, l, -1
            f = s*e( i )
            b = c*e( i )
            call dlartg( g, f, c, s, r )
            if( i.ne.m-1 )
     $         e( i+1 ) = r
            g = d( i+1 ) - p
            r = ( d( i )-g )*s + two*c*b
            p = s*r
            d( i+1 ) = g + p
            g = c*r - b
c
c           if eigenvectors are desired, then save rotations.
c
            if( icompz.gt.0 ) then
               work( i ) = c
               work( n-1+i ) = -s
            end if
c
   70    continue
c
c        if eigenvectors are desired, then apply saved rotations.
c
         if( icompz.gt.0 ) then
            mm = m - l + 1
c$$$            call dlasr( 'r', 'v', 'b', n, mm, work( l ), work( n-1+l ),
c$$$     $                  z( 1, l ), ldz )
c
c             *** New starting with version 2.5 ***
c
              call dlasr( 'r', 'v', 'b', 1, mm, work( l ),
     &                    work( n-1+l ), z( l ), 1 )
c             *************************************
         end if
c
         d( l ) = d( l ) - p
         e( l ) = g
         go to 40
c
c        eigenvalue found.
c
   80    continue
         d( l ) = p
c
         l = l + 1
         if( l.le.lend )
     $      go to 40
         go to 140
c
      else
c
c        qr iteration
c
c        look for small superdiagonal element.
c
   90    continue
         if( l.ne.lend ) then
            lendp1 = lend + 1
            do 100 m = l, lendp1, -1
               tst = abs( e( m-1 ) )**2
               if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m-1 ) )+
     $             safmin )go to 110
  100       continue
         end if
c
         m = lend
c
  110    continue
         if( m.gt.lend )
     $      e( m-1 ) = zero
         p = d( l )
         if( m.eq.l )
     $      go to 130
c
c        if remaining matrix is 2-by-2, use dlae2 or dlaev2
c        to compute its eigensystem.
c
         if( m.eq.l-1 ) then
            if( icompz.gt.0 ) then
               call dlaev2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2, c, s )
c$$$               work( m ) = c
c$$$               work( n-1+m ) = s
c$$$               call dlasr( 'r', 'v', 'f', n, 2, work( m ),
c$$$     $                     work( n-1+m ), z( 1, l-1 ), ldz )
c
c               *** New starting with version 2.5 ***
c
                tst      = z(l)
                z(l)   = c*tst - s*z(l-1)
                z(l-1) = s*tst + c*z(l-1)
c               *************************************
            else
               call dlae2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2 )
            end if
            d( l-1 ) = rt1
            d( l ) = rt2
            e( l-1 ) = zero
            l = l - 2
            if( l.ge.lend )
     $         go to 90
            go to 140
         end if
c
         if( jtot.eq.nmaxit )
     $      go to 140
         jtot = jtot + 1
c
c        form shift.
c
         g = ( d( l-1 )-p ) / ( two*e( l-1 ) )
         r = dlapy2( g, one )
         g = d( m ) - p + ( e( l-1 ) / ( g+sign( r, g ) ) )
c
         s = one
         c = one
         p = zero
c
c        inner loop
c
         lm1 = l - 1
         do 120 i = m, lm1
            f = s*e( i )
            b = c*e( i )
            call dlartg( g, f, c, s, r )
            if( i.ne.m )
     $         e( i-1 ) = r
            g = d( i ) - p
            r = ( d( i+1 )-g )*s + two*c*b
            p = s*r
            d( i ) = g + p
            g = c*r - b
c
c           if eigenvectors are desired, then save rotations.
c
            if( icompz.gt.0 ) then
               work( i ) = c
               work( n-1+i ) = s
            end if
c
  120    continue
c
c        if eigenvectors are desired, then apply saved rotations.
c
         if( icompz.gt.0 ) then
            mm = l - m + 1
c$$$            call dlasr( 'r', 'v', 'f', n, mm, work( m ), work( n-1+m ),
c$$$     $                  z( 1, m ), ldz )
c
c           *** New starting with version 2.5 ***
c
            call dlasr( 'r', 'v', 'f', 1, mm, work( m ), work( n-1+m ),
     &                  z( m ), 1 )
c           *************************************
         end if
c
         d( l ) = d( l ) - p
         e( lm1 ) = g
         go to 90
c
c        eigenvalue found.
c
  130    continue
         d( l ) = p
c
         l = l - 1
         if( l.ge.lend )
     $      go to 90
         go to 140
c
      end if
c
c     undo scaling if necessary
c
  140 continue
      if( iscale.eq.1 ) then
         call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv+1, 1,
     $                d( lsv ), n, info )
         call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv, 1, e( lsv ),
     $                n, info )
      else if( iscale.eq.2 ) then
         call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv+1, 1,
     $                d( lsv ), n, info )
         call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv, 1, e( lsv ),
     $                n, info )
      end if
c
c     check for no convergence to an eigenvalue after a total
c     of n*maxit iterations.
c
      if( jtot.lt.nmaxit )
     $   go to 10
      do 150 i = 1, n - 1
         if( e( i ).ne.zero )
     $      info = info + 1
  150 continue
      go to 190
c
c     order eigenvalues and eigenvectors.
c
  160 continue
      if( icompz.eq.0 ) then
c
c        use quick sort
c
         call dlasrt( 'i', n, d, info )
c
      else
c
c        use selection sort to minimize swaps of eigenvectors
c
         do 180 ii = 2, n
            i = ii - 1
            k = i
            p = d( i )
            do 170 j = ii, n
               if( d( j ).lt.p ) then
                  k = j
                  p = d( j )
               end if
  170       continue
            if( k.ne.i ) then
               d( k ) = d( i )
               d( i ) = p
c$$$               call dswap( n, z( 1, i ), 1, z( 1, k ), 1 )
c           *** New starting with version 2.5 ***
c
               p    = z(k)
               z(k) = z(i)
               z(i) = p
c           *************************************
            end if
  180    continue
      end if
c
  190 continue
c

      RETURN

! **********************************************************************************************************************************

c     %---------------%
c     | End of dstqrb |
c     %---------------%
c
      end subroutine dstqrb

      END MODULE ARPACK_LANCZOS_EIG
