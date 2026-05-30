! ##################################################################################################################################

      MODULE LAPACK_LIN_EQN_DPB
! --- lapack_surgery begin --- !

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                   :  ERR, F06
      USE SCONTR, ONLY                   :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                   :  HOUR, MINUTE, SEC,
     &                                      SFRAC, TSEC
      USE LAPACK_BLAS_AUX
      USE LAPACK_DPBCON_HELPER, ONLY     :  DPBCON_HELPER => DPBCON
      USE LAPACK_DPBEQU_HELPER, ONLY     :  DPBEQU_HELPER => DPBEQU
      USE LAPACK_DPBTF2_HELPER, ONLY     :  DPBTF2_HELPER => DPBTF2
      USE LAPACK_DPBTRF_KERNEL, ONLY     :  DPBTRF_KERNEL => DPBTRF
      USE LAPACK_DPBTRS_HELPER, ONLY     :  DPBTRS_HELPER => DPBTRS
      USE LAPACK_DSYTF2_HELPER, ONLY     :  DSYTF2_HELPER => DSYTF2
      USE LAPACK_POTF2_HELPER, ONLY      :  DPOTF2_HELPER => DPOTF2
      USE PARAMS, ONLY                   :  NOCOUNTS

      character(1*byte), parameter      :: cr13_dpb = char(13)

! This is the set of LAPACK routines for solving equations

!                             Ax = B

! where matrix A is a dbl prec symmetric banded matrix.
! Matrix A is decomposed into an upper triangular matrix U such that:

!                         A = U(transp)*U

! This module contains:

!     DPBEQU to equilize A by scaling using diagonals of A

!     DPBTRF to do the factorization of A by calling

!         DPBTF2 if A is unblocked to do factorization of A, or

!         DPOTF2 if A is blocked   to do factorization of A

!     DPBCON to calculate the condition number of A, given it's triangular factors

!     DPBTRS to get the solution for x given the triangular factors of A

! In addition, files in module LAPACK_BLAS_AUX are also used

      CONTAINS

! ##################################################################################################################################
! 001 LAPACK_LINEAR_EQN_DPB

! --- lapack_peeloff begin --- !
      SUBROUTINE DPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL(DOUBLE)   AMAX, SCOND
      REAL(DOUBLE)   AB( LDAB, * ), S( * )
      CALL DPBEQU_HELPER( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )
      RETURN
      END SUBROUTINE DPBEQU
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 002 LAPACK_LINEAR_EQN_DPB

! --- lapack_peeloff begin --- !
      SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL(DOUBLE)   AB( LDAB, * )

      INTEGER            I, IB, NB, iblock, numblk
      INTEGER            ILAENV
      EXTERNAL           ILAENV
      INTRINSIC          MIN

      INFO = 0
      IF( N > 0 .AND. KD >= 0 .AND. LDAB >= KD+1 ) THEN
         NB = ILAENV( 1, 'DPBTRF', UPLO, N, KD, -1, -1 )
         NB = MIN( NB, 32 )
         IF( NB > 1 .AND. NB <= KD ) THEN
            numblk = int(n/nb) + 1
            iblock = 0
            WRITE (*,*)
            DO I = 1, N, NB
               iblock = iblock + 1
               IB = MIN( NB, N-I+1 )
               IF (NOCOUNTS .NE. 'Y') THEN
                  write(sc1,12345) iblock,numblk,ib
                ENDIF
            END DO
         END IF
      END IF

      CALL DPBTRF_KERNEL( UPLO, N, KD, AB, LDAB, INFO )

12345 format(5X,'Block ',i8,' of ',i8,'. Factoring rows 1 thru: ',i8)

      RETURN
      END SUBROUTINE DPBTRF
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 003 LAPACK_LINEAR_EQN_DPB

! --- lapack_peeloff begin --- !
      SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL(DOUBLE)   AB( LDAB, * )
      CALL DPBTF2_HELPER( UPLO, N, KD, AB, LDAB, INFO )
      RETURN
      END SUBROUTINE DPBTF2
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 004 LAPACK_LINEAR_EQN_DPB

! --- lapack_peeloff begin --- !
      SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      REAL(DOUBLE)   A( LDA, * )
      CALL DPOTF2_HELPER( UPLO, N, A, LDA, INFO )
      RETURN
      END SUBROUTINE DPOTF2
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 005 LAPACK_LINEAR_EQN_DPB

! --- lapack_peeloff begin --- !
      SUBROUTINE DPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,
     $                   IWORK, INFO, itmax, dtbsv_msg )
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      CHARACTER          UPLO
      CHARACTER*1        dtbsv_msg
      INTEGER            INFO, KD, LDAB, N
      REAL(DOUBLE)   ANORM, RCOND
      INTEGER            IWORK( * )
      REAL(DOUBLE)   AB( LDAB, * ), WORK( * )
      INTEGER            itmax
      CALL DPBCON_HELPER( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,
     $                    IWORK, INFO, itmax, dtbsv_msg )
      END SUBROUTINE DPBCON
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 006 LAPACK_LINEAR_EQN_DPB

! --- lapack_peeloff begin --- !
      SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO,
     $                   dtbsv_msg )     ! my addition

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER          UPLO
      character*1        dtbsv_msg
      INTEGER            INFO, KD, LDAB, LDB, N, NRHS
      REAL(DOUBLE)   AB( LDAB, * ), B( LDB, * )
      CALL DPBTRS_HELPER( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO,
     $                    dtbsv_msg )
      RETURN
      END SUBROUTINE DPBTRS
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 007 LAPACK_LINEAR_EQN_DPB

      SUBROUTINE DPTTRF_MYSTRAN( N, D, E, INFO )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DPTTRF'

!  This is my mod of LAPACK subr DPTTRF. I have changed the statements:

!      IF( D( I ).LE.ZERO ) THEN
!  to
!      IF( D( I ).EQ.ZERO ) THEN

!  in order to allow the subr to get the LDL decomp when diag elements
!  are neqative. I need this in LINK4 for subr EST_NUMBER_OF_EIGENS
!  which needs a count of the number of diag terms < 0 in order to
!  estimate the number of eigens less than a specified value.

*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      REAL(DOUBLE)   D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  DPTTRF computes the L*D*L' factorization of a real symmetric
*  positive definite tridiagonal matrix A.  The factorization may also
*  be regarded as having the form A = U'*D*U.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix
*          A.  On exit, the n diagonal elements of the diagonal matrix
*          D from the L*D*L' factorization of A.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix A.  On exit, the (n-1) subdiagonal elements of the
*          unit bidiagonal factor L from the L*D*L' factorization of A.
*          E can also be regarded as the superdiagonal of the unit
*          bidiagonal factor U from the U'*D*U factorization of A.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite; if k < N, the factorization could not
*               be completed, while if k = N, the factorization was
*               completed, but D(N) <= 0.
**
*     .. Parameters ..
      REAL(DOUBLE)   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I4
      REAL(DOUBLE)   EI
*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DPTTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Compute the L*D*L' (or U'*D*U) factorization of A.
*
      I4 = MOD( N-1, 4 )
      DO 10 I = 1, I4
         IF( D( I ).EQ.ZERO ) THEN
            INFO = I
            GO TO 30
         END IF
         EI = E( I )
         E( I ) = EI / D( I )
         D( I+1 ) = D( I+1 ) - E( I )*EI
   10 CONTINUE
*
      DO 20 I = I4 + 1, N - 4, 4
*
*        Drop out of the loop if d(i) <= 0: the matrix is not positive
*        definite.
*
         IF( D( I ).EQ.ZERO ) THEN
            INFO = I
            GO TO 30
         END IF
*
*        Solve for e(i) and d(i+1).
*
         EI = E( I )
         E( I ) = EI / D( I )
         D( I+1 ) = D( I+1 ) - E( I )*EI
*
         IF( D( I ).EQ.ZERO ) THEN
            INFO = I + 1
            GO TO 30
         END IF
*
*        Solve for e(i+1) and d(i+2).
*
         EI = E( I+1 )
         E( I+1 ) = EI / D( I+1 )
         D( I+2 ) = D( I+2 ) - E( I+1 )*EI
*
         IF( D( I ).EQ.ZERO ) THEN
            INFO = I + 2
            GO TO 30
         END IF
*
*        Solve for e(i+2) and d(i+3).
*
         EI = E( I+2 )
         E( I+2 ) = EI / D( I+2 )
         D( I+3 ) = D( I+3 ) - E( I+2 )*EI
*
         IF( D( I ).EQ.ZERO ) THEN
            INFO = I + 3
            GO TO 30
         END IF
*
*        Solve for e(i+3) and d(i+4).
*
         EI = E( I+3 )
         E( I+3 ) = EI / D( I+3 )
         D( I+4 ) = D( I+4 ) - E( I+3 )*EI
   20 CONTINUE
*
*     Check d(n) for positive definiteness.
*
      IF( D( N ).LE.ZERO )
     $   INFO = N
*
   30 CONTINUE
      RETURN
*
*     End of DPTTRF
*
! **********************************************************************************************************************************
 9000 continue            ! My lines

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE DPTTRF_MYSTRAN

! #################################################################################################################################
! 008 LAPACK_LINEAR_EQN_DPB

! --- lapack_peeloff begin --- !
      SUBROUTINE DSYTF2( UPLO, N, A, LDA, IPIV, INFO )

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DSYTF2'
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYTF2 computes the factorization of a real symmetric matrix A using
*  the Bunch-Kaufman diagonal pivoting method:
*
*     A = U*D*U'  or  A = L*D*L'
*
*  where U (or L) is a product of permutation and unit upper (lower)
*  triangular matrices, U' is the transpose of U, and D is symmetric and
*  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, the block diagonal matrix D and the multipliers used
*          to obtain the factor U or L (see below for further details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          Details of the interchanges and the block structure of D.
*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*          interchanged and D(k,k) is a 1-by-1 diagonal block.
*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
*               has been completed, but the block diagonal matrix D is
*               exactly singular, and division by zero will occur if it
*               is used to solve a system of equations.
*
*  Further Details
*  ===============
*
*  09-29-06 - patch from
*    Bobby Cheng, MathWorks
*
*    Replace l.204 and l.372
*         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*    by
*         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN
*
*  01-01-96 - Based on modifications by
*    J. Lewis, Boeing Computer Services Company
*    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*  1-96 - Based on modifications by J. Lewis, Boeing Computer Services
*         Company
*
*  If UPLO = 'U', then A = U*D*U', where
*     U = P(n)*U(n)* ... *P(k)U(k)* ...,
*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*
*             (   I    v    0   )   k-s
*     U(k) =  (   0    I    0   )   s
*             (   0    0    I   )   n-k
*                k-s   s   n-k
*
*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
*  and A(k,k), and v overwrites A(1:k-2,k-1:k).
*
*  If UPLO = 'L', then A = L*D*L', where
*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*
*             (   I    0     0   )  k-1
*     L(k) =  (   0    I     0   )  s
*             (   0    v     I   )  n-k-s+1
*                k-1   s  n-k-s+1
*
*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IMAX, J, JMAX, K, KK, KP, KSTEP
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, R1,
     $                   ROWMAX, T, WK, WKM1, WKP1
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME

*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      CALL DSYTF2_HELPER( UPLO, N, A, LDA, IPIV, INFO )
      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE DSYTF2
! --- lapack_peeloff end --- !

! --- lapack_surgery end --- !
      END MODULE LAPACK_LIN_EQN_DPB
