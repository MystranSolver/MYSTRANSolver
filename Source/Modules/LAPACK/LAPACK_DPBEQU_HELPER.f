! ##################################################################################################################################

      MODULE LAPACK_DPBEQU_HELPER

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE LAPACK_BLAS_AUX

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL(DOUBLE)       AMAX, SCOND
      REAL(DOUBLE)       AB( LDAB, * ), S( * )

      REAL(DOUBLE)       ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      LOGICAL            UPPER
      INTEGER            I, J
      REAL(DOUBLE)       SMIN

      LOGICAL            LSAME
      EXTERNAL           LSAME

      INTRINSIC          MAX, MIN, SQRT

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N < 0 ) THEN
         INFO = -2
      ELSE IF( KD < 0 ) THEN
         INFO = -3
      ELSE IF( LDAB < KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO /= 0 ) THEN
         CALL XERBLA( 'DPBEQU', -INFO )
         GO TO 9000
      END IF

      IF( N == 0 ) THEN
         SCOND = ONE
         AMAX = ZERO
         GO TO 9000
      END IF

      IF( UPPER ) THEN
         J = KD + 1
      ELSE
         J = 1
      END IF

      S( 1 ) = AB( J, 1 )
      SMIN = S( 1 )
      AMAX = S( 1 )

      DO 10 I = 2, N
         S( I ) = AB( J, I )
         SMIN = MIN( SMIN, S( I ) )
         AMAX = MAX( AMAX, S( I ) )
   10 CONTINUE

      IF( SMIN <= ZERO ) THEN
         DO 20 I = 1, N
            IF( S( I ) <= ZERO ) THEN
               INFO = I
               GO TO 9000
            END IF
   20    CONTINUE
      ELSE
         DO 30 I = 1, N
            S( I ) = ONE / SQRT( S( I ) )
   30    CONTINUE

         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      END IF

 9000 CONTINUE

      RETURN

      END SUBROUTINE DPBEQU
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DPBEQU_HELPER
