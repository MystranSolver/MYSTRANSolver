! ##################################################################################################################################

      MODULE LAPACK_LANCZOS_EIG_HELPER

      USE PENTIUM_II_KIND, ONLY         :  BYTE, LONG, DOUBLE
      USE LAPACK_BLAS_AUX

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )

      INTEGER            INFO, LDA, M, N
      REAL(DOUBLE)       A( LDA, * ), TAU( * ), WORK( * )
      REAL(DOUBLE)       ONE
      PARAMETER          ( ONE = 1.0D+0 )
      INTEGER            I, K
      REAL(DOUBLE)       AII
      INTRINSIC          MAX, MIN

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQR2', -INFO )
         RETURN
      END IF

      K = MIN( M, N )

      DO 10 I = 1, K
         CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                TAU( I ) )
         IF( I.LT.N ) THEN
            AII = A( I, I )
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE

      RETURN
      END SUBROUTINE DGEQR2

! ##################################################################################################################################

      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )

      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
      REAL(DOUBLE)       A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      REAL(DOUBLE)       ONE
      PARAMETER          ( ONE = 1.0D+0 )
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      REAL(DOUBLE)       AII
      LOGICAL            LSAME
      EXTERNAL           LSAME
      INTRINSIC          MAX

      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )

      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF

      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN

      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF

      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF

      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
            MI = M - I + 1
            IC = I
         ELSE
            NI = N - I + 1
            JC = I
         END IF

         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ),
     $               LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE

      RETURN
      END SUBROUTINE DORM2R

      END MODULE LAPACK_LANCZOS_EIG_HELPER
