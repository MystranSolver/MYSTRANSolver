! ##################################################################################################################################

MODULE LAPACK_DLAS_MISC_HELPER

USE PENTIUM_II_KIND, ONLY         :  DOUBLE

CONTAINS

! ##################################################################################################################################

DOUBLE PRECISION FUNCTION DLAMC3_HELPER( A, B )

REAL(DOUBLE)   A, B

DLAMC3_HELPER = A + B
RETURN

END FUNCTION DLAMC3_HELPER

! ##################################################################################################################################

SUBROUTINE DLASET_HELPER( UPLO, M, N, ALPHA, BETA, A, LDA )

CHARACTER          UPLO
INTEGER            LDA, M, N
REAL(DOUBLE)       ALPHA, BETA
REAL(DOUBLE)       A( LDA, * )
INTEGER            I, J
LOGICAL            LSAME
EXTERNAL           LSAME
INTRINSIC          MIN

IF( LSAME( UPLO, 'U' ) ) THEN
   DO 20 J = 2, N
      DO 10 I = 1, MIN( J-1, M )
         A( I, J ) = ALPHA
10    CONTINUE
20 CONTINUE
ELSE IF( LSAME( UPLO, 'L' ) ) THEN
   DO 40 J = 1, MIN( M, N )
      DO 30 I = J + 1, M
         A( I, J ) = ALPHA
30    CONTINUE
40 CONTINUE
ELSE
   DO 60 J = 1, N
      DO 50 I = 1, M
         A( I, J ) = ALPHA
50    CONTINUE
60 CONTINUE
END IF

DO 70 I = 1, MIN( M, N )
   A( I, I ) = BETA
70 CONTINUE

RETURN

END SUBROUTINE DLASET_HELPER

! ##################################################################################################################################

SUBROUTINE DLASWP_HELPER( N, A, LDA, K1, K2, IPIV, INCX )

INTEGER            INCX, K1, K2, LDA, N
INTEGER            IPIV( * )
REAL(DOUBLE)       A( LDA, * )
INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
REAL(DOUBLE)       TEMP

IF( INCX.GT.0 ) THEN
   IX0 = K1
   I1 = K1
   I2 = K2
   INC = 1
ELSE IF( INCX.LT.0 ) THEN
   IX0 = 1 + ( 1-K2 )*INCX
   I1 = K2
   I2 = K1
   INC = -1
ELSE
   RETURN
END IF

N32 = ( N / 32 )*32
IF( N32.NE.0 ) THEN
   DO 30 J = 1, N32, 32
      IX = IX0
      DO 20 I = I1, I2, INC
         IP = IPIV( IX )
         IF( IP.NE.I ) THEN
            DO 10 K = J, J + 31
               TEMP = A( I, K )
               A( I, K ) = A( IP, K )
               A( IP, K ) = TEMP
10          CONTINUE
         END IF
         IX = IX + INCX
20    CONTINUE
30 CONTINUE
END IF
IF( N32.NE.N ) THEN
   N32 = N32 + 1
   IX = IX0
   DO 50 I = I1, I2, INC
      IP = IPIV( IX )
      IF( IP.NE.I ) THEN
         DO 40 K = N32, N
            TEMP = A( I, K )
            A( I, K ) = A( IP, K )
            A( IP, K ) = TEMP
40       CONTINUE
      END IF
      IX = IX + INCX
50 CONTINUE
END IF

RETURN

END SUBROUTINE DLASWP_HELPER

END MODULE LAPACK_DLAS_MISC_HELPER
