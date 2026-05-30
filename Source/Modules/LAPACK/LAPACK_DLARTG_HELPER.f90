! ##################################################################################################################################

MODULE LAPACK_DLARTG_HELPER

USE PENTIUM_II_KIND, ONLY         :  DOUBLE

CONTAINS

! ##################################################################################################################################

SUBROUTINE DLARTG_HELPER( F, G, CS, SN, R )

REAL(DOUBLE)       CS, F, G, R, SN
REAL(DOUBLE)       ZERO
PARAMETER          ( ZERO = 0.0D0 )
REAL(DOUBLE)       ONE
PARAMETER          ( ONE = 1.0D0 )
REAL(DOUBLE)       TWO
PARAMETER          ( TWO = 2.0D0 )
LOGICAL            FIRST
INTEGER            COUNT, I
REAL(DOUBLE)       EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE

REAL(DOUBLE)       DLAMCH
EXTERNAL           DLAMCH
INTRINSIC          ABS, INT, LOG, MAX, SQRT

SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
DATA               FIRST / .TRUE. /

IF( FIRST ) THEN
   FIRST = .FALSE.
   SAFMIN = DLAMCH( 'S' )
   EPS = DLAMCH( 'E' )
   SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / &
            LOG( DLAMCH( 'B' ) ) / TWO )
   SAFMX2 = ONE / SAFMN2
END IF
IF( G.EQ.ZERO ) THEN
   CS = ONE
   SN = ZERO
   R = F
ELSE IF( F.EQ.ZERO ) THEN
   CS = ZERO
   SN = ONE
   R = G
ELSE
   F1 = F
   G1 = G
   SCALE = MAX( ABS( F1 ), ABS( G1 ) )
   IF( SCALE.GE.SAFMX2 ) THEN
      COUNT = 0
10    CONTINUE
      COUNT = COUNT + 1
      F1 = F1*SAFMN2
      G1 = G1*SAFMN2
      SCALE = MAX( ABS( F1 ), ABS( G1 ) )
      IF( SCALE.GE.SAFMX2 ) GO TO 10
      R = SQRT( F1**2+G1**2 )
      CS = F1 / R
      SN = G1 / R
      DO 20 I = 1, COUNT
         R = R*SAFMX2
20    CONTINUE
   ELSE IF( SCALE.LE.SAFMN2 ) THEN
      COUNT = 0
30    CONTINUE
      COUNT = COUNT + 1
      F1 = F1*SAFMX2
      G1 = G1*SAFMX2
      SCALE = MAX( ABS( F1 ), ABS( G1 ) )
      IF( SCALE.LE.SAFMN2 ) GO TO 30
      R = SQRT( F1**2+G1**2 )
      CS = F1 / R
      SN = G1 / R
      DO 40 I = 1, COUNT
         R = R*SAFMN2
40    CONTINUE
   ELSE
      R = SQRT( F1**2+G1**2 )
      CS = F1 / R
      SN = G1 / R
   END IF
   IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
      CS = -CS
      SN = -SN
      R = -R
   END IF
END IF

END SUBROUTINE DLARTG_HELPER

END MODULE LAPACK_DLARTG_HELPER
