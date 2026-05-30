! ##################################################################################################################################

      MODULE LAPACK_DPBTRF_KERNEL

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE LAPACK_BLAS_AUX
      USE LAPACK_DPBTF2_HELPER, ONLY     :  DPBTF2
      USE LAPACK_POTF2_HELPER, ONLY      :  DPOTF2

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL(DOUBLE)       AB( LDAB, * )

      REAL(DOUBLE)       ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 32, LDWORK = NBMAX+1 )

      INTEGER            I, I2, I3, IB, II, J, JJ, NB
      REAL(DOUBLE)       WORK( LDWORK, NBMAX )

      LOGICAL            LSAME
      EXTERNAL           LSAME

      INTEGER            ILAENV
      EXTERNAL           ILAENV

      INTRINSIC          MIN

      INFO = 0
      IF( ( .NOT.LSAME( UPLO, 'U' ) ) .AND.
     $    ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -1
      ELSE IF( N < 0 ) THEN
         INFO = -2
      ELSE IF( KD < 0 ) THEN
         INFO = -3
      ELSE IF( LDAB < KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO /= 0 ) THEN
         CALL XERBLA( 'DPBTRF', -INFO )
         GO TO 9000
      END IF

      IF( N == 0 )
     &   GO TO 9000

      NB = ILAENV( 1, 'DPBTRF', UPLO, N, KD, -1, -1 )
      NB = MIN( NB, NBMAX )

      IF( NB <= 1 .OR. NB > KD ) THEN
         CALL DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
      ELSE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, NB
               DO 10 I = 1, J - 1
                  WORK( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE

            DO 70 I = 1, N, NB
               IB = MIN( NB, N-I+1 )

               CALL DPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
               IF( II /= 0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB <= N ) THEN
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )

                  IF( I2 > 0 ) THEN
                     CALL DTRSM( 'Left', 'Upper', 'Transpose',
     $                           'Non-unit', IB, I2, ONE, AB( KD+1, I ),
     $                           LDAB-1, AB( KD+1-IB, I+IB ), LDAB-1 )

                     CALL DSYRK( 'Upper', 'Transpose', I2, IB, -ONE,
     $                           AB( KD+1-IB, I+IB ), LDAB-1, ONE,
     $                           AB( KD+1, I+IB ), LDAB-1 )
                  END IF

                  IF( I3 > 0 ) THEN
                     DO 40 JJ = 1, I3
                        DO 30 II = JJ, IB
                           WORK( II, JJ ) = AB( II-JJ+1, JJ+I+KD-1 )
   30                   CONTINUE
   40                CONTINUE

                     CALL DTRSM( 'Left', 'Upper', 'Transpose',
     $                           'Non-unit', IB, I3, ONE, AB( KD+1, I ),
     $                           LDAB-1, WORK, LDWORK )

                     IF( I2 > 0 )
     $                  CALL DGEMM( 'Transpose', 'No Transpose', I2, I3,
     $                              IB, -ONE, AB( KD+1-IB, I+IB ),
     $                              LDAB-1, WORK, LDWORK, ONE,
     $                              AB( 1+IB, I+KD ), LDAB-1 )

                     CALL DSYRK( 'Upper', 'Transpose', I3, IB, -ONE,
     $                           WORK, LDWORK, ONE, AB( KD+1, I+KD ),
     $                           LDAB-1 )

                     DO 60 JJ = 1, I3
                        DO 50 II = JJ, IB
                           AB( II-JJ+1, JJ+I+KD-1 ) = WORK( II, JJ )
   50                   CONTINUE
   60                CONTINUE
                  END IF
               END IF
   70       CONTINUE
         ELSE
            DO 90 J = 1, NB
               DO 80 I = J + 1, NB
                  WORK( I, J ) = ZERO
   80          CONTINUE
   90       CONTINUE

            DO 140 I = 1, N, NB
               IB = MIN( NB, N-I+1 )

               CALL DPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
               IF( II /= 0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF( I+IB <= N ) THEN
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )

                  IF( I2 > 0 ) THEN
                     CALL DTRSM( 'Right', 'Lower', 'Transpose',
     $                           'Non-unit', I2, IB, ONE, AB( 1, I ),
     $                           LDAB-1, AB( 1+IB, I ), LDAB-1 )

                     CALL DSYRK( 'Lower', 'No Transpose', I2, IB, -ONE,
     $                           AB( 1+IB, I ), LDAB-1, ONE,
     $                           AB( 1, I+IB ), LDAB-1 )
                  END IF

                  IF( I3 > 0 ) THEN
                     DO 110 JJ = 1, IB
                        DO 100 II = 1, MIN( JJ, I3 )
                           WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
  100                   CONTINUE
  110                CONTINUE

                     CALL DTRSM( 'Right', 'Lower', 'Transpose',
     $                           'Non-unit', I3, IB, ONE, AB( 1, I ),
     $                           LDAB-1, WORK, LDWORK )

                     IF( I2 > 0 )
     $                  CALL DGEMM( 'No transpose', 'Transpose', I3, I2,
     $                              IB, -ONE, WORK, LDWORK,
     $                              AB( 1+IB, I ), LDAB-1, ONE,
     $                              AB( 1+KD-IB, I+IB ), LDAB-1 )

                     CALL DSYRK( 'Lower', 'No Transpose', I3, IB, -ONE,
     $                           WORK, LDWORK, ONE, AB( 1, I+KD ),
     $                           LDAB-1 )

                     DO 130 JJ = 1, IB
                        DO 120 II = 1, MIN( JJ, I3 )
                           AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
  120                   CONTINUE
  130                CONTINUE
                  END IF
               END IF
  140       CONTINUE
         END IF
      END IF
      GO TO 9000

  150 CONTINUE

 9000 CONTINUE

      RETURN

      END SUBROUTINE DPBTRF
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DPBTRF_KERNEL
