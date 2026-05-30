! ##################################################################################################################################

      MODULE LAPACK_DPBTF2_HELPER

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE LAPACK_BLAS_AUX

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL(DOUBLE)       AB( LDAB, * )

      REAL(DOUBLE)       ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )

      LOGICAL            UPPER
      INTEGER            J, KLD, KN
      REAL(DOUBLE)       AJJ

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
         CALL XERBLA( 'DPBTF2', -INFO )
         GO TO 9000
      END IF

      IF( N == 0 )
     &   GO TO 9000

      KLD = MAX( 1, LDAB-1 )

      IF( UPPER ) THEN
         DO 10 J = 1, N
            AJJ = AB( KD+1, J )
            IF( AJJ <= ZERO ) THEN
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ

            KN = MIN( KD, N-J )
            IF( KN > 0 ) THEN
               CALL DSCAL( KN, ONE / AJJ, AB( KD, J+1 ), KLD )
               CALL DSYR( 'Upper', KN, -ONE, AB( KD, J+1 ), KLD,
     $                    AB( KD+1, J+1 ), KLD )
            END IF
   10    CONTINUE
      ELSE
         DO 20 J = 1, N
            AJJ = AB( 1, J )
            IF( AJJ <= ZERO )
     $         GO TO 30
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ

            KN = MIN( KD, N-J )
            IF( KN > 0 ) THEN
               CALL DSCAL( KN, ONE / AJJ, AB( 2, J ), 1 )
               CALL DSYR( 'Lower', KN, -ONE, AB( 2, J ), 1,
     $                    AB( 1, J+1 ), KLD )
            END IF
   20    CONTINUE
      END IF
      GO TO 9000

   30 CONTINUE
      INFO = J

 9000 CONTINUE

      RETURN

      END SUBROUTINE DPBTF2
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DPBTF2_HELPER
