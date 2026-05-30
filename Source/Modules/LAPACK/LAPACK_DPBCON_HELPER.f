! ##################################################################################################################################

      MODULE LAPACK_DPBCON_HELPER

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                   :  BLNK_SUB_NAM
      USE LAPACK_BLAS_AUX

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,
     $                   IWORK, INFO, itmax, dtbsv_msg )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: subr_name = 'DPBCON'

      CHARACTER          UPLO
      CHARACTER*1        dtbsv_msg
      INTEGER            INFO, KD, LDAB, N
      REAL(DOUBLE)       ANORM, RCOND
      INTEGER            IWORK( * )
      REAL(DOUBLE)       AB( LDAB, * ), WORK( * )

      REAL(DOUBLE)       ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )

      LOGICAL            UPPER
      CHARACTER          NORMIN
      INTEGER            IX, KASE
      INTEGER            iter_num, itmax
      REAL(DOUBLE)       AINVNM, SCALE, SCALEL, SCALEU, SMLNUM

      LOGICAL            LSAME
      EXTERNAL           LSAME

      REAL(DOUBLE)       DLAMCH
      EXTERNAL           DLAMCH

      INTRINSIC          ABS

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N < 0 ) THEN
         INFO = -2
      ELSE IF( KD < 0 ) THEN
         INFO = -3
      ELSE IF( LDAB < KD+1 ) THEN
         INFO = -5
      ELSE IF( ANORM < ZERO ) THEN
         INFO = -6
      END IF
      IF( INFO /= 0 ) THEN
         CALL XERBLA( 'DPBCON', -INFO )
         GO TO 9000
      END IF

      RCOND = ZERO
      IF( N == 0 ) THEN
         RCOND = ONE
         GO TO 9000
      ELSE IF( ANORM == ZERO ) THEN
         GO TO 9000
      END IF

      SMLNUM = DLAMCH( 'Safe minimum' )

      KASE = 0
      iter_num = 0
      NORMIN = 'N'
   10 CONTINUE
      iter_num = iter_num + 1
      CALL DLACON( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, itmax )
      IF( KASE /= 0 ) THEN
         IF( UPPER ) THEN
            CALL DLATBS( 'Upper', 'Transpose', 'Non-unit', NORMIN, N,
     $                   KD, AB, LDAB, WORK, SCALEL, WORK( 2*N+1 ),
     $                   INFO, iter_num, dtbsv_msg )
            NORMIN = 'Y'
            CALL DLATBS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N,
     $                   KD, AB, LDAB, WORK, SCALEU, WORK( 2*N+1 ),
     $                   INFO, iter_num, dtbsv_msg )
         ELSE
            CALL DLATBS( 'Lower', 'No transpose', 'Non-unit', NORMIN, N,
     $                   KD, AB, LDAB, WORK, SCALEL, WORK( 2*N+1 ),
     $                   INFO, iter_num, dtbsv_msg )
            NORMIN = 'Y'
            CALL DLATBS( 'Lower', 'Transpose', 'Non-unit', NORMIN, N,
     $                   KD, AB, LDAB, WORK, SCALEU, WORK( 2*N+1 ),
     $                   INFO, iter_num, dtbsv_msg )
         END IF

         SCALE = SCALEL*SCALEU
         IF( SCALE /= ONE ) THEN
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE < ABS( WORK( IX ) )*SMLNUM .OR. SCALE == ZERO )
     $         GO TO 20
            CALL DRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF

      IF( AINVNM /= ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM

   20 CONTINUE

 9000 CONTINUE

      RETURN

      END SUBROUTINE DPBCON
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DPBCON_HELPER
