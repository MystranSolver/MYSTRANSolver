! ##################################################################################################################################

      MODULE LAPACK_DLACPY_HELPER

      USE PENTIUM_II_KIND, ONLY         :  DOUBLE

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DLACPY_HELPER( UPLO, M, N, A, LDA, B, LDB )

      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
      REAL(DOUBLE)   A( LDA, * ), B( LDB, * )
      INTEGER            I, J
      LOGICAL            LSAME
      EXTERNAL           LSAME

      INTRINSIC          MIN

      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN

      END SUBROUTINE DLACPY_HELPER

      END MODULE LAPACK_DLACPY_HELPER
