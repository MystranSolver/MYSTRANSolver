! ##################################################################################################################################

MODULE LAPACK_DLARFB_HELPER

USE PENTIUM_II_KIND, ONLY         :  DOUBLE

CONTAINS

! ##################################################################################################################################

SUBROUTINE DLARFB_HELPER( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                          T, LDT, C, LDC, WORK, LDWORK )

CHARACTER          SIDE, TRANS, DIRECT, STOREV
INTEGER            M, N, K, LDV, LDT, LDC, LDWORK
REAL(DOUBLE)       V( LDV, * ), T( LDT, * ), C( LDC, * ), WORK( LDWORK, * )
REAL(DOUBLE)       ONE
PARAMETER          ( ONE = 1.0D+0 )
CHARACTER          TRANST
INTEGER            I, J
LOGICAL            LSAME
EXTERNAL           LSAME
EXTERNAL           DCOPY, DTRMM, DGEMM

IF( M.LE.0 .OR. N.LE.0 ) RETURN

IF( LSAME( TRANS, 'N' ) ) THEN
   TRANST = 'T'
ELSE
   TRANST = 'N'
END IF

IF( LSAME( STOREV, 'C' ) ) THEN

   IF( LSAME( DIRECT, 'F' ) ) THEN

      IF( LSAME( SIDE, 'L' ) ) THEN

         DO 10 J = 1, K
            CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
10       CONTINUE

         CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                     K, ONE, V, LDV, WORK, LDWORK )
         IF( M.GT.K ) THEN
            CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                        ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                        ONE, WORK, LDWORK )
         END IF

         CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                     ONE, T, LDT, WORK, LDWORK )

         IF( M.GT.K ) THEN
            CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                        -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                        C( K+1, 1 ), LDC )
         END IF

         CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                     ONE, V, LDV, WORK, LDWORK )

         DO 30 J = 1, K
            DO 20 I = 1, N
               C( J, I ) = C( J, I ) - WORK( I, J )
20          CONTINUE
30       CONTINUE

      ELSE IF( LSAME( SIDE, 'R' ) ) THEN

         DO 40 J = 1, K
            CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
40       CONTINUE

         CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                     K, ONE, V, LDV, WORK, LDWORK )
         IF( N.GT.K ) THEN
            CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                        ONE, WORK, LDWORK )
         END IF

         CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                     ONE, T, LDT, WORK, LDWORK )

         IF( N.GT.K ) THEN
            CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                        -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                        C( 1, K+1 ), LDC )
         END IF

         CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                     ONE, V, LDV, WORK, LDWORK )

         DO 60 J = 1, K
            DO 50 I = 1, M
               C( I, J ) = C( I, J ) - WORK( I, J )
50          CONTINUE
60       CONTINUE
      END IF

   ELSE

      IF( LSAME( SIDE, 'L' ) ) THEN

         DO 70 J = 1, K
            CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
70       CONTINUE

         CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
         IF( M.GT.K ) THEN
            CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
         END IF

         CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                     ONE, T, LDT, WORK, LDWORK )

         IF( M.GT.K ) THEN
            CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
         END IF

         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                     ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )

         DO 90 J = 1, K
            DO 80 I = 1, N
               C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
80          CONTINUE
90       CONTINUE

      ELSE IF( LSAME( SIDE, 'R' ) ) THEN

         DO 100 J = 1, K
            CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
100      CONTINUE

         CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
         IF( N.GT.K ) THEN
            CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
         END IF

         CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                     ONE, T, LDT, WORK, LDWORK )

         IF( N.GT.K ) THEN
            CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
         END IF

         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                     ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )

         DO 120 J = 1, K
            DO 110 I = 1, M
               C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
110         CONTINUE
120      CONTINUE
      END IF
   END IF

ELSE IF( LSAME( STOREV, 'R' ) ) THEN

   IF( LSAME( DIRECT, 'F' ) ) THEN

      IF( LSAME( SIDE, 'L' ) ) THEN

         DO 130 J = 1, K
            CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
130      CONTINUE

         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                     ONE, V, LDV, WORK, LDWORK )
         IF( M.GT.K ) THEN
            CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
                        WORK, LDWORK )
         END IF

         CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                     ONE, T, LDT, WORK, LDWORK )

         IF( M.GT.K ) THEN
            CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                        C( K+1, 1 ), LDC )
         END IF

         CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, K, &
                     ONE, V, LDV, WORK, LDWORK )

         DO 150 J = 1, K
            DO 140 I = 1, N
               C( J, I ) = C( J, I ) - WORK( I, J )
140         CONTINUE
150      CONTINUE

      ELSE IF( LSAME( SIDE, 'R' ) ) THEN

         DO 160 J = 1, K
            CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
160      CONTINUE

         CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                     ONE, V, LDV, WORK, LDWORK )
         IF( N.GT.K ) THEN
            CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                        ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                        ONE, WORK, LDWORK )
         END IF

         CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                     ONE, T, LDT, WORK, LDWORK )

         IF( N.GT.K ) THEN
            CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                        C( 1, K+1 ), LDC )
         END IF

         CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, K, &
                     ONE, V, LDV, WORK, LDWORK )

         DO 180 J = 1, K
            DO 170 I = 1, M
               C( I, J ) = C( I, J ) - WORK( I, J )
170         CONTINUE
180      CONTINUE
      END IF

   ELSE

      IF( LSAME( SIDE, 'L' ) ) THEN

         DO 190 J = 1, K
            CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
190      CONTINUE

         CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
         IF( M.GT.K ) THEN
            CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                        C, LDC, V, LDV, ONE, WORK, LDWORK )
         END IF

         CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                     ONE, T, LDT, WORK, LDWORK )

         IF( M.GT.K ) THEN
            CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                        V, LDV, WORK, LDWORK, ONE, C, LDC )
         END IF

         CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, K, &
                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )

         DO 210 J = 1, K
            DO 200 I = 1, N
               C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
200         CONTINUE
210      CONTINUE

      ELSE IF( LSAME( SIDE, 'R' ) ) THEN

         DO 220 J = 1, K
            CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
220      CONTINUE

         CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
         IF( N.GT.K ) THEN
            CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
         END IF

         CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                     ONE, T, LDT, WORK, LDWORK )

         IF( N.GT.K ) THEN
            CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
         END IF

         CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, K, &
                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )

         DO 240 J = 1, K
            DO 230 I = 1, M
               C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
230         CONTINUE
240      CONTINUE
      END IF
   END IF
END IF

END SUBROUTINE DLARFB_HELPER

END MODULE LAPACK_DLARFB_HELPER
