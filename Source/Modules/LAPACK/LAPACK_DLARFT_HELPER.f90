! ##################################################################################################################################

MODULE LAPACK_DLARFT_HELPER

USE PENTIUM_II_KIND, ONLY         :  DOUBLE

CONTAINS

! ##################################################################################################################################

SUBROUTINE DLARFT_HELPER( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )

CHARACTER          DIRECT, STOREV
INTEGER            K, LDT, LDV, N
REAL(DOUBLE)       T( LDT, * ), TAU( * ), V( LDV, * )
REAL(DOUBLE)       ONE, ZERO
PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
INTEGER            I, J
REAL(DOUBLE)       VII
LOGICAL            LSAME
EXTERNAL           LSAME
EXTERNAL           DGEMV, DTRMV

IF( N.EQ.0 ) RETURN

IF( LSAME( DIRECT, 'F' ) ) THEN
   DO 20 I = 1, K
      IF( TAU( I ).EQ.ZERO ) THEN
         DO 10 J = 1, I
            T( J, I ) = ZERO
10       CONTINUE
      ELSE
         VII = V( I, I )
         V( I, I ) = ONE
         IF( LSAME( STOREV, 'C' ) ) THEN
            CALL DGEMV( 'Transpose', N-I+1, I-1, -TAU( I ), &
                        V( I, 1 ), LDV, V( I, I ), 1, ZERO, &
                        T( 1, I ), 1 )
         ELSE
            CALL DGEMV( 'No transpose', I-1, N-I+1, -TAU( I ), &
                        V( 1, I ), LDV, V( I, I ), LDV, ZERO, &
                        T( 1, I ), 1 )
         END IF
         V( I, I ) = VII
         CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
                     LDT, T( 1, I ), 1 )
         T( I, I ) = TAU( I )
      END IF
20 CONTINUE
ELSE
   DO 40 I = K, 1, -1
      IF( TAU( I ).EQ.ZERO ) THEN
         DO 30 J = I, K
            T( J, I ) = ZERO
30       CONTINUE
      ELSE
         IF( I.LT.K ) THEN
            IF( LSAME( STOREV, 'C' ) ) THEN
               VII = V( N-K+I, I )
               V( N-K+I, I ) = ONE
               CALL DGEMV( 'Transpose', N-K+I, K-I, -TAU( I ), &
                           V( 1, I+1 ), LDV, V( 1, I ), 1, ZERO, &
                           T( I+1, I ), 1 )
               V( N-K+I, I ) = VII
            ELSE
               VII = V( I, N-K+I )
               V( I, N-K+I ) = ONE
               CALL DGEMV( 'No transpose', K-I, N-K+I, -TAU( I ), &
                           V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, &
                           T( I+1, I ), 1 )
               V( I, N-K+I ) = VII
            END IF
            CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
         END IF
         T( I, I ) = TAU( I )
      END IF
40 CONTINUE
END IF

END SUBROUTINE DLARFT_HELPER

END MODULE LAPACK_DLARFT_HELPER
