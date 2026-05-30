! ##################################################################################################################################

      MODULE LAPACK_DGETRS_HELPER

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE LAPACK_BLAS_AUX

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DGETRS_HELPER( TRANS, N, NRHS, A, LDA, IPIV, B, LDB,
     $                          INFO )

      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      REAL(DOUBLE)       A( LDA, * ), B( LDB, * )

      REAL(DOUBLE)       ONE
      PARAMETER          ( ONE = 1.0D+0 )

      LOGICAL            NOTRAN

      LOGICAL            LSAME
      EXTERNAL           LSAME

      INTRINSIC          MAX

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF

      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN

      IF( NOTRAN ) THEN

         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )

         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )

         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE

         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )

         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )

         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF

      RETURN

      END SUBROUTINE DGETRS_HELPER
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DGETRS_HELPER
