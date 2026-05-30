! ##################################################################################################################################

      MODULE LAPACK_LANCZOS_EIG
! --- lapack_surgery begin --- !

      USE PENTIUM_II_KIND, ONLY              :  BYTE, LONG, DOUBLE
      USE LAPACK_LANCZOS_EIG_HELPER, ONLY    :  DGEQR2_HELPER => DGEQR2,                                       &
     &                                          DORM2R_HELPER => DORM2R

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
      INTEGER            INFO, LDA, M, N
      REAL(DOUBLE)       A( LDA, * ), TAU( * ), WORK( * )
      CALL DGEQR2_HELPER( M, N, A, LDA, TAU, WORK, INFO )
      RETURN
      END SUBROUTINE DGEQR2

! ##################################################################################################################################

      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
      REAL(DOUBLE)       A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      CALL DORM2R_HELPER( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                    WORK, INFO )
      RETURN
      END SUBROUTINE DORM2R

! --- lapack_surgery end --- !
      END MODULE LAPACK_LANCZOS_EIG
