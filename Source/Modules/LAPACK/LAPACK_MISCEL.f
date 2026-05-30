! ##################################################################################################################################

      MODULE LAPACK_MISCEL
! --- lapack_surgery begin --- !

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE LAPACK_BLAS_AUX
      USE LAPACK_DSTEV_HELPER, ONLY      :  DSTEV_HELPER  => DSTEV
      USE LAPACK_DTRTRS_HELPER, ONLY     :  DTRTRS_HELPER => DTRTRS

! This is a set of LAPACK routines that are used in several other modules but are not BLAS or auxiliary routines
! The routines included herein are:

!     DSTERF: computes all eigenvalues of a symmetric tridiagonal matrix using the Pal-Walker-Kahan variant of the QL or QR alg.

!     DSTEQR: computes all eigenvalues and, optionally, eigenvectors of a sym tridiag matrix using the implicit QL or QR method.

!     DSTEV : computes all eigenvalues and, optionally, eigenvectors of a real symmetric tridiagonal matrix A

!     DTRTRS: solves a triangular system of the form

!                A * X = B  or  A**T * X = B,

!             where A is a triangular matrix of order N, and B is an N-by-NRHS matrix.  A check is made to verify A is nonsingular.


      CONTAINS

! ##################################################################################################################################
! 003 LAPACK_MISCEL

! --- lapack_peeloff begin --- !
      SUBROUTINE DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
      CHARACTER          JOBZ
      INTEGER            INFO, LDZ, N
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
      CALL DSTEV_HELPER( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
      RETURN
      END SUBROUTINE DSTEV
! --- lapack_peeloff end --- !

! ##################################################################################################################################
! 004 LAPACK_MISCEL

! --- lapack_peeloff begin --- !
      SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
     $                   INFO )
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      CALL DTRTRS_HELPER( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,
     $                    INFO )
      RETURN
      END SUBROUTINE DTRTRS
! --- lapack_peeloff end --- !

! --- lapack_surgery end --- !
      END MODULE LAPACK_MISCEL

