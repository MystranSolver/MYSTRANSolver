! ##################################################################################################################################

MODULE LAPACK_DLARFG_HELPER

USE PENTIUM_II_KIND, ONLY         :  DOUBLE
USE LAPACK_DLAPY2_HELPER, ONLY    :  DLAPY2_HELPER

CONTAINS

! ##################################################################################################################################

SUBROUTINE DLARFG_HELPER( N, ALPHA, X, INCX, TAU )

INTEGER            INCX, N
REAL(DOUBLE)       ALPHA, TAU
REAL(DOUBLE)       X( * )
REAL(DOUBLE)       ONE, ZERO
PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
INTEGER            J, KNT
REAL(DOUBLE)       BETA, RSAFMN, SAFMIN, XNORM

INTRINSIC          ABS, SIGN
REAL(DOUBLE)       DLAMCH, DNRM2
EXTERNAL           DLAMCH, DNRM2
EXTERNAL           DSCAL

IF( N.LE.1 ) THEN
   TAU = ZERO
   RETURN
END IF

XNORM = DNRM2( N-1, X, INCX )

IF( XNORM.EQ.ZERO ) THEN
   TAU = ZERO
ELSE
   BETA = -SIGN( DLAPY2_HELPER( ALPHA, XNORM ), ALPHA )
   SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
   IF( ABS( BETA ).LT.SAFMIN ) THEN
      RSAFMN = ONE / SAFMIN
      KNT = 0
10    CONTINUE
      KNT = KNT + 1
      CALL DSCAL( N-1, RSAFMN, X, INCX )
      BETA = BETA*RSAFMN
      ALPHA = ALPHA*RSAFMN
      IF( ABS( BETA ).LT.SAFMIN ) GO TO 10

      XNORM = DNRM2( N-1, X, INCX )
      BETA = -SIGN( DLAPY2_HELPER( ALPHA, XNORM ), ALPHA )
      TAU = ( BETA-ALPHA ) / BETA
      CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )

      ALPHA = BETA
      DO 20 J = 1, KNT
         ALPHA = ALPHA*SAFMIN
20    CONTINUE
   ELSE
      TAU = ( BETA-ALPHA ) / BETA
      CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
      ALPHA = BETA
   END IF
END IF

END SUBROUTINE DLARFG_HELPER

END MODULE LAPACK_DLARFG_HELPER
