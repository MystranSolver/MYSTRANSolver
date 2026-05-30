! ##################################################################################################################################

MODULE LAPACK_DLAGTS_HELPER

USE PENTIUM_II_KIND, ONLY         :  DOUBLE

CONTAINS

! ##################################################################################################################################

SUBROUTINE DLAGTS_HELPER( JOB, N, A, B, C, D, IN, Y, TOL, INFO )

INTEGER            INFO, JOB, N
REAL(DOUBLE)       TOL
INTEGER            IN( * )
REAL(DOUBLE)       A( * ), B( * ), C( * ), D( * ), Y( * )
REAL(DOUBLE)       ZERO, ONE
PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
LOGICAL            NOTRAN
INTEGER            K
REAL(DOUBLE)       AK, PERT, SFMIN, TEMP, ABSAK, BIGNUM

INTRINSIC          ABS, SIGN
REAL(DOUBLE)       DLAMCH
EXTERNAL           DLAMCH

INFO = 0
IF( N.EQ.0 ) RETURN

IF( JOB.LT.-2 .OR. JOB.GT.2 ) THEN
   INFO = -1
   CALL XERBLA( 'DLAGTS', -INFO )
   RETURN
END IF

IF( JOB.EQ.0 ) THEN
   TOL = ZERO
END IF

NOTRAN = JOB.LT.0
IF( ABS(TOL).LE.ZERO ) THEN
   TOL = ABS( A( 1 ) )
   IF( N.GT.1 ) TOL = MAX( TOL, ABS( A( 2 ) ), ABS( B( 1 ) ) )
   DO 10 K = 3, N
      TOL = MAX( TOL, ABS( A( K ) ), ABS( B( K-1 ) ), ABS( D( K-2 ) ) )
10 CONTINUE
   TOL = TOL*DLAMCH( 'Epsilon' )
   IF( TOL.EQ.ZERO ) TOL = DLAMCH( 'Safe minimum' )
END IF

SFMIN = DLAMCH( 'Safe minimum' )
BIGNUM = ONE / SFMIN

IF( NOTRAN ) THEN
   DO 20 K = 2, N
      IF( IN( K-1 ).EQ.0 ) THEN
         Y( K ) = Y( K ) - C( K-1 )*Y( K-1 )
      ELSE
         TEMP = Y( K-1 )
         Y( K-1 ) = Y( K )
         Y( K ) = TEMP - C( K-1 )*Y( K )
      END IF
20 CONTINUE

   IF( JOB.EQ.-1 ) THEN
      DO 30 K = N, 1, -1
         IF( K.LE.N-2 ) THEN
            TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
         ELSE IF( K.EQ.N-1 ) THEN
            TEMP = Y( K ) - B( K )*Y( K+1 )
         ELSE
            TEMP = Y( K )
         END IF
         Y( K ) = TEMP / A( K )
30    CONTINUE
   ELSE
      DO 50 K = N, 1, -1
         IF( K.LE.N-2 ) THEN
            TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
         ELSE IF( K.EQ.N-1 ) THEN
            TEMP = Y( K ) - B( K )*Y( K+1 )
         ELSE
            TEMP = Y( K )
         END IF
         AK = A( K )
         PERT = SIGN( TOL, AK )
40       CONTINUE
         ABSAK = ABS( AK )
         IF( ABSAK.LT.ONE ) THEN
            IF( ABSAK.LT.SFMIN ) THEN
               IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) THEN
                  AK = AK + PERT
                  PERT = 2*PERT
                  GO TO 40
               ELSE
                  TEMP = TEMP*BIGNUM
                  AK = AK*BIGNUM
               END IF
            ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
               AK = AK + PERT
               PERT = 2*PERT
               GO TO 40
            END IF
         END IF
         Y( K ) = TEMP / AK
50    CONTINUE
   END IF
ELSE
   IF( JOB.EQ.2 ) THEN
      K = 1
      TEMP = Y( 1 )
      AK = A( 1 )
      ABSAK = ABS( AK )
      IF( ABSAK.LT.ONE ) THEN
         IF( ABSAK.LT.SFMIN ) THEN
            IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) THEN
               INFO = K
               RETURN
            ELSE
               TEMP = TEMP*BIGNUM
               AK = AK*BIGNUM
            END IF
         ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
            INFO = K
            RETURN
         END IF
      END IF
      Y( 1 ) = TEMP / AK
      IF( N.GE.2 ) THEN
         K = 2
         TEMP = Y( 2 ) - B( 1 )*Y( 1 )
         AK = A( 2 )
         ABSAK = ABS( AK )
         IF( ABSAK.LT.ONE ) THEN
            IF( ABSAK.LT.SFMIN ) THEN
               IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) THEN
                  INFO = K
                  RETURN
               ELSE
                  TEMP = TEMP*BIGNUM
                  AK = AK*BIGNUM
               END IF
            ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
               INFO = K
               RETURN
            END IF
         END IF
         Y( 2 ) = TEMP / AK
      END IF
      DO 60 K = 3, N
         TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
         AK = A( K )
         ABSAK = ABS( AK )
         IF( ABSAK.LT.ONE ) THEN
            IF( ABSAK.LT.SFMIN ) THEN
               IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) THEN
                  INFO = K
                  RETURN
               ELSE
                  TEMP = TEMP*BIGNUM
                  AK = AK*BIGNUM
               END IF
            ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
               INFO = K
               RETURN
            END IF
         END IF
         Y( K ) = TEMP / AK
60    CONTINUE
   ELSE
      K = 1
      TEMP = Y( 1 )
      AK = A( 1 )
      PERT = SIGN( TOL, AK )
70    CONTINUE
      ABSAK = ABS( AK )
      IF( ABSAK.LT.ONE ) THEN
         IF( ABSAK.LT.SFMIN ) THEN
            IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) THEN
               AK = AK + PERT
               PERT = 2*PERT
               GO TO 70
            ELSE
               TEMP = TEMP*BIGNUM
               AK = AK*BIGNUM
            END IF
         ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
            AK = AK + PERT
            PERT = 2*PERT
            GO TO 70
         END IF
      END IF
      Y( 1 ) = TEMP / AK
      IF( N.GE.2 ) THEN
         K = 2
         TEMP = Y( 2 ) - B( 1 )*Y( 1 )
         AK = A( 2 )
         PERT = SIGN( TOL, AK )
75       CONTINUE
         ABSAK = ABS( AK )
         IF( ABSAK.LT.ONE ) THEN
            IF( ABSAK.LT.SFMIN ) THEN
               IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) THEN
                  AK = AK + PERT
                  PERT = 2*PERT
                  GO TO 75
               ELSE
                  TEMP = TEMP*BIGNUM
                  AK = AK*BIGNUM
               END IF
            ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
               AK = AK + PERT
               PERT = 2*PERT
               GO TO 75
            END IF
         END IF
         Y( 2 ) = TEMP / AK
      END IF
      DO 80 K = 3, N
         TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
         AK = A( K )
         PERT = SIGN( TOL, AK )
77       CONTINUE
         ABSAK = ABS( AK )
         IF( ABSAK.LT.ONE ) THEN
            IF( ABSAK.LT.SFMIN ) THEN
               IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) THEN
                  AK = AK + PERT
                  PERT = 2*PERT
                  GO TO 77
               ELSE
                  TEMP = TEMP*BIGNUM
                  AK = AK*BIGNUM
               END IF
            ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
               AK = AK + PERT
               PERT = 2*PERT
               GO TO 77
            END IF
         END IF
         Y( K ) = TEMP / AK
80    CONTINUE
   END IF

   DO 90 K = N, 2, -1
      IF( IN( K-1 ).EQ.0 ) THEN
         Y( K-1 ) = Y( K-1 ) - C( K-1 )*Y( K )
      ELSE
         TEMP = Y( K-1 )
         Y( K-1 ) = Y( K )
         Y( K ) = TEMP - C( K-1 )*Y( K )
      END IF
90 CONTINUE
END IF

END SUBROUTINE DLAGTS_HELPER

END MODULE LAPACK_DLAGTS_HELPER
