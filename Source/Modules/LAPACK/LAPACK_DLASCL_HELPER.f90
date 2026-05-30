! ##################################################################################################################################

MODULE LAPACK_DLASCL_HELPER

USE PENTIUM_II_KIND, ONLY         :  DOUBLE

CONTAINS

! ##################################################################################################################################

SUBROUTINE DLASCL_HELPER( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )

CHARACTER          TYPE
INTEGER            INFO, KL, KU, LDA, M, N
REAL(DOUBLE)       CFROM, CTO
REAL(DOUBLE)       A( LDA, * )
REAL(DOUBLE)       ZERO, ONE
PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
LOGICAL            DONE
INTEGER            I, ITYPE, J, K1, K2, K3, K4
REAL(DOUBLE)       BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
LOGICAL            LSAME
EXTERNAL           LSAME
REAL(DOUBLE)       DLAMCH
EXTERNAL           DLAMCH
INTRINSIC          ABS, MAX, MIN

INFO = 0

IF( LSAME( TYPE, 'G' ) ) THEN
   ITYPE = 0
ELSE IF( LSAME( TYPE, 'L' ) ) THEN
   ITYPE = 1
ELSE IF( LSAME( TYPE, 'U' ) ) THEN
   ITYPE = 2
ELSE IF( LSAME( TYPE, 'H' ) ) THEN
   ITYPE = 3
ELSE IF( LSAME( TYPE, 'B' ) ) THEN
   ITYPE = 4
ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
   ITYPE = 5
ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
   ITYPE = 6
ELSE
   ITYPE = -1
END IF

IF( ITYPE.EQ.-1 ) THEN
   INFO = -1
ELSE IF( CFROM.EQ.ZERO ) THEN
   INFO = -4
ELSE IF( M.LT.0 ) THEN
   INFO = -6
ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. &
         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
   INFO = -7
ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
   INFO = -9
ELSE IF( ITYPE.GE.4 ) THEN
   IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
      INFO = -2
   ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. &
            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) THEN
      INFO = -3
   ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. &
            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. &
            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
      INFO = -9
   END IF
END IF

IF( INFO.NE.0 ) THEN
   CALL XERBLA( 'DLASCL', -INFO )
   RETURN
END IF

IF( N.EQ.0 .OR. M.EQ.0 ) RETURN

SMLNUM = DLAMCH( 'S' )
BIGNUM = ONE / SMLNUM

CFROMC = CFROM
CTOC = CTO

10 CONTINUE
CFROM1 = CFROMC*SMLNUM
CTO1 = CTOC / BIGNUM
IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
   MUL = SMLNUM
   DONE = .FALSE.
   CFROMC = CFROM1
ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
   MUL = BIGNUM
   DONE = .FALSE.
   CTOC = CTO1
ELSE
   MUL = CTOC / CFROMC
   DONE = .TRUE.
END IF

IF( ITYPE.EQ.0 ) THEN
   DO 30 J = 1, N
      DO 20 I = 1, M
         A( I, J ) = A( I, J )*MUL
20    CONTINUE
30 CONTINUE
ELSE IF( ITYPE.EQ.1 ) THEN
   DO 50 J = 1, N
      DO 40 I = J, M
         A( I, J ) = A( I, J )*MUL
40    CONTINUE
50 CONTINUE
ELSE IF( ITYPE.EQ.2 ) THEN
   DO 70 J = 1, N
      DO 60 I = 1, MIN( J, M )
         A( I, J ) = A( I, J )*MUL
60    CONTINUE
70 CONTINUE
ELSE IF( ITYPE.EQ.3 ) THEN
   DO 90 J = 1, N
      DO 80 I = 1, MIN( J+1, M )
         A( I, J ) = A( I, J )*MUL
80    CONTINUE
90 CONTINUE
ELSE IF( ITYPE.EQ.4 ) THEN
   K3 = KL + 1
   K4 = N + 1
   DO 110 J = 1, N
      DO 100 I = 1, MIN( K3, K4-J )
         A( I, J ) = A( I, J )*MUL
100   CONTINUE
110 CONTINUE
ELSE IF( ITYPE.EQ.5 ) THEN
   K1 = KU + 2
   K3 = KU + 1
   DO 130 J = 1, N
      DO 120 I = MAX( K1-J, 1 ), K3
         A( I, J ) = A( I, J )*MUL
120   CONTINUE
130 CONTINUE
ELSE IF( ITYPE.EQ.6 ) THEN
   K1 = KL + KU + 2
   K2 = KL + 1
   K3 = 2*KL + KU + 1
   K4 = KL + KU + 1 + M
   DO 150 J = 1, N
      DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
         A( I, J ) = A( I, J )*MUL
140   CONTINUE
150 CONTINUE
END IF

IF( .NOT.DONE ) GO TO 10

RETURN

END SUBROUTINE DLASCL_HELPER

END MODULE LAPACK_DLASCL_HELPER
