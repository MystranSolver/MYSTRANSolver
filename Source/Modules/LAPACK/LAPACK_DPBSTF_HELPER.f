! ##################################################################################################################################

      MODULE LAPACK_DPBSTF_HELPER

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DPBSTF_HELPER( UPLO, N, KD, AB, LDAB, INFO )

      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL(DOUBLE)       AB( LDAB, * )

      REAL(DOUBLE)       ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )

      LOGICAL            UPPER
      INTEGER            J, KLD, KM, M
      REAL(DOUBLE)       AJJ

      LOGICAL            LSAME
      EXTERNAL           LSAME

      INTRINSIC          MAX, MIN, SQRT

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPBSTF', -INFO )
         RETURN
      END IF

      IF( N.EQ.0 ) RETURN

      KLD = MAX( 1, LDAB-1 )
      M = ( N+KD ) / 2

      IF( UPPER ) THEN

         DO 10 J = N, M + 1, -1
            AJJ = AB( KD+1, J )
            IF( AJJ.LE.ZERO ) GO TO 50
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
            KM = MIN( J-1, KD )
            CALL DSCAL( KM, ONE / AJJ, AB( KD+1-KM, J ), 1 )
            CALL DSYR( 'Upper', KM, -ONE, AB( KD+1-KM, J ), 1,
     $                 AB( KD+1, J-KM ), KLD )
   10    CONTINUE

         DO 20 J = 1, M
            AJJ = AB( KD+1, J )
            IF( AJJ.LE.ZERO ) GO TO 50
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
            KM = MIN( KD, M-J )
            IF( KM.GT.0 ) THEN
               CALL DSCAL( KM, ONE / AJJ, AB( KD, J+1 ), KLD )
               CALL DSYR( 'Upper', KM, -ONE, AB( KD, J+1 ), KLD,
     $                    AB( KD+1, J+1 ), KLD )
            END IF
   20    CONTINUE
      ELSE

         DO 30 J = N, M + 1, -1
            AJJ = AB( 1, J )
            IF( AJJ.LE.ZERO ) GO TO 50
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
            KM = MIN( J-1, KD )
            CALL DSCAL( KM, ONE / AJJ, AB( KM+1, J-KM ), KLD )
            CALL DSYR( 'Lower', KM, -ONE, AB( KM+1, J-KM ), KLD,
     $                 AB( 1, J-KM ), KLD )
   30    CONTINUE

         DO 40 J = 1, M
            AJJ = AB( 1, J )
            IF( AJJ.LE.ZERO ) GO TO 50
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
            KM = MIN( KD, M-J )
            IF( KM.GT.0 ) THEN
               CALL DSCAL( KM, ONE / AJJ, AB( 2, J ), 1 )
               CALL DSYR( 'Lower', KM, -ONE, AB( 2, J ), 1,
     $                    AB( 1, J+1 ), KLD )
            END IF
   40    CONTINUE
      END IF
      RETURN

   50 CONTINUE
      INFO = J

      RETURN

      END SUBROUTINE DPBSTF_HELPER

      END MODULE LAPACK_DPBSTF_HELPER
