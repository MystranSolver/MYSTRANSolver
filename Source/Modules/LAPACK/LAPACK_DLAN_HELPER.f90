! ##################################################################################################################################

MODULE LAPACK_DLAN_HELPER

USE PENTIUM_II_KIND, ONLY         :  BYTE, DOUBLE
USE IOUNT1, ONLY                  :  SC1
USE PARAMS, ONLY                  :  NOCOUNTS
USE LAPACK_DLASSQ_HELPER

CHARACTER(1*BYTE), PARAMETER      :: cr13_lba_dlan = CHAR(13)

CONTAINS

! ##################################################################################################################################

DOUBLE PRECISION FUNCTION DLANSB_HELPER( NORM, UPLO, N, K, AB, LDAB, WORK )

CHARACTER          NORM, UPLO
INTEGER            K, LDAB, N
REAL(DOUBLE)       AB( LDAB, * ), WORK( * )
REAL(DOUBLE)       ONE, ZERO
PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
INTEGER            I, J, L
REAL(DOUBLE)       ABSA, SCALE, SUM, VALUE
LOGICAL            LSAME
EXTERNAL           LSAME
INTRINSIC          ABS, MAX, MIN, SQRT

IF( N.EQ.0 ) THEN
   VALUE = ZERO
ELSE IF( LSAME( NORM, 'M' ) ) THEN
   VALUE = ZERO
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO 20 J = 1, N
         IF (NOCOUNTS .NE. 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') J, N, cr13_lba_dlan
         ENDIF
         DO 10 I = MAX( K+2-J, 1 ), K + 1
            VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
10       CONTINUE
20    CONTINUE
   ELSE
      DO 40 J = 1, N
         IF (NOCOUNTS .NE. 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') J, N, cr13_lba_dlan
         ENDIF
         DO 30 I = 1, MIN( N+1-J, K+1 )
            VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
30       CONTINUE
40    CONTINUE
   END IF
ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
   VALUE = ZERO
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO 60 J = 1, N
         IF (NOCOUNTS .NE. 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') J, N, cr13_lba_dlan
         ENDIF
         SUM = ZERO
         L = K + 1 - J
         DO 50 I = MAX( 1, J-K ), J - 1
            ABSA = ABS( AB( L+I, J ) )
            SUM = SUM + ABSA
            WORK( I ) = WORK( I ) + ABSA
50       CONTINUE
         WORK( J ) = SUM + ABS( AB( K+1, J ) )
60    CONTINUE
      DO 70 I = 1, N
         IF (NOCOUNTS .NE. 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') I, N, cr13_lba_dlan
         ENDIF
         VALUE = MAX( VALUE, WORK( I ) )
70    CONTINUE
   ELSE
      DO 80 I = 1, N
         WORK( I ) = ZERO
80    CONTINUE
      DO 100 J = 1, N
         IF (NOCOUNTS .NE. 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') J, N, cr13_lba_dlan
         ENDIF
         SUM = WORK( J ) + ABS( AB( 1, J ) )
         L = 1 - J
         DO 90 I = J + 1, MIN( N, J+K )
            ABSA = ABS( AB( L+I, J ) )
            SUM = SUM + ABSA
            WORK( I ) = WORK( I ) + ABSA
90       CONTINUE
         VALUE = MAX( VALUE, SUM )
100   CONTINUE
   END IF
ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
   SCALE = ZERO
   SUM = ONE
   IF( K.GT.0 ) THEN
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 110 J = 2, N
            IF (NOCOUNTS .NE. 'Y') THEN
               WRITE(SC1,12345,ADVANCE='NO') J, N, cr13_lba_dlan
            ENDIF
            CALL DLASSQ_HELPER( MIN( J-1, K ), AB( MAX( K+2-J, 1 ), J ), 1, SCALE, SUM )
110      CONTINUE
         L = K + 1
      ELSE
         DO 120 J = 1, N - 1
            IF (NOCOUNTS .NE. 'Y') THEN
               WRITE(SC1,12345,ADVANCE='NO') J, N, cr13_lba_dlan
            ENDIF
            CALL DLASSQ_HELPER( MIN( N-J, K ), AB( 2, J ), 1, SCALE, SUM )
120      CONTINUE
         L = 1
      END IF
      SUM = 2*SUM
   ELSE
      L = 1
   END IF
   CALL DLASSQ_HELPER( N, AB( L, 1 ), LDAB, SCALE, SUM )
   VALUE = SCALE*SQRT( SUM )
END IF

DLANSB_HELPER = VALUE
RETURN

12345 FORMAT(5X,'Row ',I8,' of ',I8, A)

END FUNCTION DLANSB_HELPER

! ##################################################################################################################################

DOUBLE PRECISION FUNCTION DLANSY_HELPER( NORM, UPLO, N, A, LDA, WORK )

CHARACTER          NORM, UPLO
INTEGER            LDA, N
REAL(DOUBLE)       A( LDA, * ), WORK( * )
REAL(DOUBLE)       ONE, ZERO
PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
INTEGER            I, J
REAL(DOUBLE)       ABSA, SCALE, SUM, VALUE
LOGICAL            LSAME
EXTERNAL           LSAME
INTRINSIC          ABS, MAX, SQRT

IF( N.EQ.0 ) THEN
   VALUE = ZERO
ELSE IF( LSAME( NORM, 'M' ) ) THEN
   VALUE = ZERO
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO 20 J = 1, N
         DO 10 I = 1, J
            VALUE = MAX( VALUE, ABS( A( I, J ) ) )
10       CONTINUE
20    CONTINUE
   ELSE
      DO 40 J = 1, N
         DO 30 I = J, N
            VALUE = MAX( VALUE, ABS( A( I, J ) ) )
30       CONTINUE
40    CONTINUE
   END IF
ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
   VALUE = ZERO
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO 60 J = 1, N
         SUM = ZERO
         DO 50 I = 1, J - 1
            ABSA = ABS( A( I, J ) )
            SUM = SUM + ABSA
            WORK( I ) = WORK( I ) + ABSA
50       CONTINUE
         WORK( J ) = SUM + ABS( A( J, J ) )
60    CONTINUE
      DO 70 I = 1, N
         VALUE = MAX( VALUE, WORK( I ) )
70    CONTINUE
   ELSE
      DO 80 I = 1, N
         WORK( I ) = ZERO
80    CONTINUE
      DO 100 J = 1, N
         SUM = WORK( J ) + ABS( A( J, J ) )
         DO 90 I = J + 1, N
            ABSA = ABS( A( I, J ) )
            SUM = SUM + ABSA
            WORK( I ) = WORK( I ) + ABSA
90       CONTINUE
         VALUE = MAX( VALUE, SUM )
100   CONTINUE
   END IF
ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
   SCALE = ZERO
   SUM = ONE
   IF( LSAME( UPLO, 'U' ) ) THEN
      DO 110 J = 2, N
         CALL DLASSQ_HELPER( J-1, A( 1, J ), 1, SCALE, SUM )
110   CONTINUE
   ELSE
      DO 120 J = 1, N - 1
         CALL DLASSQ_HELPER( N-J, A( J+1, J ), 1, SCALE, SUM )
120   CONTINUE
   END IF
   SUM = 2*SUM
   CALL DLASSQ_HELPER( N, A, LDA+1, SCALE, SUM )
   VALUE = SCALE*SQRT( SUM )
END IF

DLANSY_HELPER = VALUE
RETURN

END FUNCTION DLANSY_HELPER

END MODULE LAPACK_DLAN_HELPER
