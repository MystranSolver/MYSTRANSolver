! ##################################################################################################################################

      MODULE LAPACK_DGETRF_HELPER

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE LAPACK_BLAS_AUX
      USE LAPACK_DGETF2_HELPER

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DGETRF_HELPER( M, N, A, LDA, IPIV, INFO )

      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      REAL(DOUBLE)       A( LDA, * )

      REAL(DOUBLE)       ONE
      PARAMETER          ( ONE = 1.0D+0 )

      INTEGER            I, IINFO, J, JB, NB

      INTEGER            ILAENV
      EXTERNAL           ILAENV

      INTRINSIC          MAX, MIN

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF

      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN

      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN

         CALL DGETF2_HELPER( M, N, A, LDA, IPIV, INFO )
      ELSE

         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )

            CALL DGETF2_HELPER( M-J+1, JB, A( J, J ), LDA, IPIV( J ),
     $                          IINFO )

            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE

            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )

            IF( J+JB.LE.N ) THEN

               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )

               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN

                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF

      RETURN

      END SUBROUTINE DGETRF_HELPER
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DGETRF_HELPER
