! ##################################################################################################################################

      MODULE LAPACK_DPOTRI_HELPER

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE

      USE LAPACK_DLAUUM_HELPER, ONLY     :  DLAUUM

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DPOTRI( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * )
      LOGICAL            LSAME
      EXTERNAL           LSAME
      EXTERNAL           DTRTRI
      INTRINSIC          MAX

      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRI', -INFO )
         RETURN
      END IF

      IF( N.EQ.0 )
     $   RETURN

      CALL DTRTRI( UPLO, 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN

      CALL DLAUUM( UPLO, N, A, LDA, INFO )

      RETURN
      END SUBROUTINE DPOTRI
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DPOTRI_HELPER
