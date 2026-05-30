! ##################################################################################################################################

MODULE LAPACK_DLACON_HELPER

USE PENTIUM_II_KIND, ONLY         :  DOUBLE

CONTAINS

! ##################################################################################################################################

SUBROUTINE DLACON_HELPER( N, V, X, ISGN, EST, KASE, itmax )

INTEGER            KASE, N
INTEGER            itmax
REAL(DOUBLE)       EST
INTEGER            ISGN( * )
REAL(DOUBLE)       V( * ), X( * )
REAL(DOUBLE)       ZERO, ONE, TWO
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
INTEGER            I, ITER, J, JLAST, JUMP
REAL(DOUBLE)       ALTSGN, ESTOLD, TEMP

INTRINSIC          ABS, DBLE, NINT, SIGN
SAVE

IF( KASE.EQ.0 ) THEN
   DO 10 I = 1, N
      X( I ) = ONE / DBLE( N )
10 CONTINUE
   KASE = 1
   JUMP = 1
   RETURN
END IF

GO TO ( 20, 40, 70, 110, 140 ) JUMP

20 CONTINUE
IF( N.EQ.1 ) THEN
   V( 1 ) = X( 1 )
   EST = ABS( V( 1 ) )
   GO TO 150
END IF
EST = DASUM( N, X, 1 )

DO 30 I = 1, N
   X( I ) = SIGN( ONE, X( I ) )
   ISGN( I ) = NINT( X( I ) )
30 CONTINUE
KASE = 2
JUMP = 2
RETURN

40 CONTINUE
J = IDAMAX( N, X, 1 )
ITER = 2

50 CONTINUE
DO 60 I = 1, N
   X( I ) = ZERO
60 CONTINUE
X( J ) = ONE
KASE = 1
JUMP = 3
RETURN

70 CONTINUE
CALL DCOPY( N, X, 1, V, 1 )
ESTOLD = EST
EST = DASUM( N, V, 1 )
DO 80 I = 1, N
   IF( NINT( SIGN( ONE, X( I ) ) ).NE.ISGN( I ) ) GO TO 90
80 CONTINUE
GO TO 120

90 CONTINUE
IF( EST.LE.ESTOLD ) GO TO 120

DO 100 I = 1, N
   X( I ) = SIGN( ONE, X( I ) )
   ISGN( I ) = NINT( X( I ) )
100 CONTINUE
KASE = 2
JUMP = 4
RETURN

110 CONTINUE
JLAST = J
J = IDAMAX( N, X, 1 )
IF( ( X( JLAST ).NE.ABS( X( J ) ) ) .AND. ( ITER.LT.ITMAX ) ) THEN
   ITER = ITER + 1
   GO TO 50
END IF

120 CONTINUE
ALTSGN = ONE
DO 130 I = 1, N
   X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) )
   ALTSGN = -ALTSGN
130 CONTINUE
KASE = 1
JUMP = 5
RETURN

140 CONTINUE
TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) )
IF( TEMP.GT.EST ) THEN
   CALL DCOPY( N, X, 1, V, 1 )
   EST = TEMP
END IF

150 CONTINUE
KASE = 0
RETURN

END SUBROUTINE DLACON_HELPER

END MODULE LAPACK_DLACON_HELPER
