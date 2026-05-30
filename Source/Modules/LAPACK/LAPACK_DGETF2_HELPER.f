! ##################################################################################################################################

      MODULE LAPACK_DGETF2_HELPER

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE PARAMS, ONLY                   :  EPSIL
      USE LAPACK_BLAS_AUX

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DGETF2_HELPER( M, N, A, LDA, IPIV, INFO )

      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      REAL(DOUBLE)       A( LDA, * )

      REAL(DOUBLE)       ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )

      INTEGER            J, JP

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
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF

      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN

      DO 10 J = 1, MIN( M, N )
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( DABS(A( JP, J )) > EPSIL(2)) THEN
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
         ELSE IF( INFO.EQ.0 ) THEN
            INFO = J
         END IF
         IF( J.LT.MIN( M, N ) ) THEN
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE

      RETURN

      END SUBROUTINE DGETF2_HELPER
! --- lapack_peeloff end --- !

      END MODULE LAPACK_DGETF2_HELPER
