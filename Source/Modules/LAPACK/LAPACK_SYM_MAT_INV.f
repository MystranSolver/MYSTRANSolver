! ##################################################################################################################################

      MODULE LAPACK_SYM_MAT_INV
! --- lapack_surgery begin --- !

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE

      USE LAPACK_BLAS_AUX
      USE LAPACK_DLAUUM_HELPER, ONLY     :  DLAUUM_HELPER => DLAUUM
      USE LAPACK_DPOTRF_HELPER, ONLY     :  DPOTRF_HELPER => DPOTRF
      USE LAPACK_DPOTRI_HELPER, ONLY     :  DPOTRI_HELPER => DPOTRI
      USE LAPACK_DTRTI2_HELPER, ONLY     :  DTRTI2_HELPER => DTRTI2

! This is a set of LAPACK routines that are used in inverting symmetric matrices (not band matrices)

      CONTAINS

! ##################################################################################################################################
! 001 LAPACK_SYM_MAT_INV

! --- lapack_peeloff begin --- !
      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * )
      CALL DPOTRF_HELPER( UPLO, N, A, LDA, INFO )
      RETURN
      END SUBROUTINE DPOTRF
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 002 LAPACK_SYM_MAT_INV

! --- lapack_peeloff begin --- !
      SUBROUTINE DPOTRI( UPLO, N, A, LDA, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * )
      CALL DPOTRI_HELPER( UPLO, N, A, LDA, INFO )
      RETURN
      END SUBROUTINE DPOTRI
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 003 LAPACK_SYM_MAT_INV

! --- lapack_peeloff begin --- !
      SUBROUTINE DLAUUM( UPLO, N, A, LDA, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * )
      CALL DLAUUM_HELPER( UPLO, N, A, LDA, INFO )
      RETURN
      END SUBROUTINE DLAUUM
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 005 LAPACK_SYM_MAT_INV

! --- lapack_peeloff begin --- !
      SUBROUTINE DLAUU2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLAUU2 computes the product U * U' or L' * L, where the triangular
*  factor U or L is stored in the upper or lower triangular part of
*  the array A.
*
*  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
*  overwriting the factor U in A.
*  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
*  overwriting the factor L in A.
*
*  This is the unblocked form of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the triangular factor stored in the array A
*          is upper or lower triangular:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the triangular factor U or L.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular factor U or L.
*          On exit, if UPLO = 'U', the upper triangle of A is
*          overwritten with the upper triangle of the product U * U';
*          if UPLO = 'L', the lower triangle of A is overwritten with
*          the lower triangle of the product L' * L.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      DOUBLE PRECISION   AII
*     ..
*     .. External Functions ..

      LOGICAL            LSAME
      EXTERNAL           LSAME

*     ..
*     .. External Subroutines ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAUU2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Compute the product U * U'.
*
         DO 10 I = 1, N
            AII = A( I, I )
            IF( I.LT.N ) THEN
               A( I, I ) = DDOT( N-I+1, A( I, I ), LDA, A( I, I ), LDA )
               CALL DGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ),
     $                     LDA, A( I, I+1 ), LDA, AII, A( 1, I ), 1 )
            ELSE
               CALL DSCAL( I, AII, A( 1, I ), 1 )
            END IF
   10    CONTINUE
*
      ELSE
*
*        Compute the product L' * L.
*
         DO 20 I = 1, N
            AII = A( I, I )
            IF( I.LT.N ) THEN
               A( I, I ) = DDOT( N-I+1, A( I, I ), 1, A( I, I ), 1 )
               CALL DGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $                     A( I+1, I ), 1, AII, A( I, 1 ), LDA )
            ELSE
               CALL DSCAL( I, AII, A( I, 1 ), LDA )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DLAUU2
*
      END SUBROUTINE DLAUU2
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 006 LAPACK_SYM_MAT_INV

! --- lapack_peeloff begin --- !
      SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * )
      CALL DTRTI2_HELPER( UPLO, DIAG, N, A, LDA, INFO )
      RETURN
      END SUBROUTINE DTRTI2
! --- lapack_peeloff end --- !

! --- lapack_surgery end --- !
      END MODULE LAPACK_SYM_MAT_INV

