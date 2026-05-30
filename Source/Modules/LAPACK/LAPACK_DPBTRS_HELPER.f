! ##################################################################################################################################

      MODULE LAPACK_DPBTRS_HELPER

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE LAPACK_BLAS_AUX

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO,
     &                   dtbsv_msg )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER          UPLO
      CHARACTER*1        dtbsv_msg
      INTEGER            INFO, KD, LDAB, LDB, N, NRHS
      REAL(DOUBLE)       AB( LDAB, * ), B( LDB, * )

      LOGICAL            UPPER
      INTEGER            J

      LOGICAL            LSAME
      EXTERNAL           LSAME

      INTRINSIC          MAX

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N < 0 ) THEN
         INFO = -2
      ELSE IF( KD < 0 ) THEN
         INFO = -3
      ELSE IF( NRHS < 0 ) THEN
         INFO = -4
      ELSE IF( LDAB < KD+1 ) THEN
         INFO = -6
      ELSE IF( LDB < MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO /= 0 ) THEN
         CALL XERBLA( 'DPBTRS', -INFO )
         GO TO 9000
      END IF

      IF( N == 0 .OR. NRHS == 0 )
     &   GO TO 9000

      IF( UPPER ) THEN
         DO 10 J = 1, NRHS
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1, dtbsv_msg )
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1, dtbsv_msg )
   10    CONTINUE
      ELSE
         DO 20 J = 1, NRHS
            CALL DTBSV( 'Lower', 'No transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1, dtbsv_msg )
            CALL DTBSV( 'Lower', 'Transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1, dtbsv_msg )
   20    CONTINUE
      END IF

 9000 CONTINUE

      RETURN

      END SUBROUTINE DPBTRS
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DPBTRS_HELPER
