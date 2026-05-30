! ##################################################################################################################################

      MODULE LAPACK_DTRTRS_HELPER

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
     $                   INFO )

      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      LOGICAL            NOUNIT
      LOGICAL            LSAME
      EXTERNAL           LSAME
      EXTERNAL           DTRSM, XERBLA
      INTRINSIC          MAX

      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.
     $         LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRS', -INFO )
         RETURN
      END IF

      IF( N.EQ.0 )
     $   RETURN

      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
      END IF
      INFO = 0

      CALL DTRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B,
     $            LDB )

      RETURN
      END SUBROUTINE DTRTRS
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DTRTRS_HELPER
