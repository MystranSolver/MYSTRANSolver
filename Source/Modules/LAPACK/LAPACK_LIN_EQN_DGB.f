! ##################################################################################################################################

      MODULE LAPACK_LIN_EQN_DGB
! --- lapack_surgery begin --- !

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE LAPACK_LIN_EQN_DGB_KERNEL, ONLY: DGBTRF_KERNEL => DGBTRF,
     &                                     DGBTRS_KERNEL => DGBTRS,
     &                                     DGBTF2_KERNEL => DGBTF2

! This facade preserves the historical MYSTRAN module name while delegating
! all numerical work to LAPACK_LIN_EQN_DGB_KERNEL.

      CONTAINS

! ##################################################################################################################################
! 001 LAPACK_LINEAR_EQN_DGB

! --- lapack_peeloff begin --- !
      SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      INTEGER            INFO, KL, KU, LDAB, M, N
      INTEGER            IPIV( * )
      REAL(DOUBLE)       AB( LDAB, * )

      CALL DGBTRF_KERNEL( M, N, KL, KU, AB, LDAB, IPIV, INFO )

      RETURN

      END SUBROUTINE DGBTRF
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 002 LAPACK_LINEAR_EQN_DGB

! --- lapack_peeloff begin --- !
      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
     $                   INFO, dtbsv_msg )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      CHARACTER          TRANS, dtbsv_msg
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
      INTEGER            IPIV( * )
      REAL(DOUBLE)       AB( LDAB, * ), B( LDB, * )

      CALL DGBTRS_KERNEL( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B,
     $                    LDB, INFO, dtbsv_msg )

      RETURN

      END SUBROUTINE DGBTRS
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 003 LAPACK_LINEAR_EQN_DGB

! --- lapack_peeloff begin --- !
      SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      INTEGER            INFO, KL, KU, LDAB, M, N
      INTEGER            IPIV( * )
      REAL(DOUBLE)       AB( LDAB, * )

      CALL DGBTF2_KERNEL( M, N, KL, KU, AB, LDAB, IPIV, INFO )

      RETURN

      END SUBROUTINE DGBTF2
! --- lapack_peeloff end --- !

! --- lapack_surgery end --- !
      END MODULE LAPACK_LIN_EQN_DGB
