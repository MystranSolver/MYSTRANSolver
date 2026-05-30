! ##################################################################################################################################

      MODULE LAPACK_LIN_EQN_DGE
! --- lapack_surgery begin --- !

      USE PENTIUM_II_KIND, ONLY          :  DOUBLE
      USE LAPACK_DGETF2_HELPER
      USE LAPACK_DGETRF_HELPER
      USE LAPACK_DGETRI_HELPER
      USE LAPACK_DGETRS_HELPER

! This is the set of LAPACK routines for solving equations

!                             AX = B

! where matrix A is a dbl prec general matrix (i.e., not symmetric). Matrix A is decomposed into an upper triangular matrix U
! such that (P is a permutation matrix):

!                              A = P*L*U

! This module contains the following subroutines:

!        DGETRF: driver for the triangular decomp of A

!        DGETRS: does the forward-backward substition. A seperate call to DGETRS must be made for each right hand side

!        DGETF2: which is called by DGETRF to do the actual decomp

      CONTAINS

! ##################################################################################################################################
! 001 LAPACK_LINEAR_EQN_DGE

! --- lapack_peeloff begin --- !
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      REAL(DOUBLE)       A( LDA, * )

      CALL DGETRF_HELPER( M, N, A, LDA, IPIV, INFO )

      END SUBROUTINE DGETRF
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 002 LAPACK_LINEAR_EQN_DGE

! --- lapack_peeloff begin --- !
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      INTEGER            INFO, LDA, LWORK, N
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )

      CALL DGETRI_HELPER( N, A, LDA, IPIV, WORK, LWORK, INFO )

      END SUBROUTINE DGETRI
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 003 LAPACK_LINEAR_EQN_DGE

! --- lapack_peeloff begin --- !
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      REAL(DOUBLE)       A( LDA, * ), B( LDB, * )

      CALL DGETRS_HELPER( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

      END SUBROUTINE DGETRS
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 004 LAPACK_LINEAR_EQN_DGE

! --- lapack_peeloff begin --- !
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      REAL(DOUBLE)       A( LDA, * )

      CALL DGETF2_HELPER( M, N, A, LDA, IPIV, INFO )

      END SUBROUTINE DGETF2
! --- lapack_peeloff end --- !

! --- lapack_surgery end --- !
      END MODULE LAPACK_LIN_EQN_DGE
