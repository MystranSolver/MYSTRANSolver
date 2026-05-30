module LAPACK_MISCEL

  use PENTIUM_II_KIND, only : BYTE, LONG, DOUBLE

  implicit none

  interface
    subroutine DSTEQR(COMPZ, N, D, E, Z, LDZ, WORK, INFO)
      use PENTIUM_II_KIND, only : LONG, DOUBLE
      implicit none
      character(1), intent(in) :: COMPZ
      integer(LONG), intent(in) :: N, LDZ
      integer(LONG), intent(out) :: INFO
      real(DOUBLE), intent(inout) :: D(*), E(*), Z(LDZ,*), WORK(*)
    end subroutine DSTEQR

    subroutine DSTEV(JOBZ, N, D, E, Z, LDZ, WORK, INFO)
      use PENTIUM_II_KIND, only : LONG, DOUBLE
      implicit none
      character(1), intent(in) :: JOBZ
      integer(LONG), intent(in) :: N, LDZ
      integer(LONG), intent(out) :: INFO
      real(DOUBLE), intent(inout) :: D(*), E(*), Z(LDZ,*), WORK(*)
    end subroutine DSTEV

    subroutine DTRTRS(UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO)
      use PENTIUM_II_KIND, only : LONG, DOUBLE
      implicit none
      character(1), intent(in) :: UPLO, TRANS, DIAG
      integer(LONG), intent(in) :: N, NRHS, LDA, LDB
      integer(LONG), intent(out) :: INFO
      real(DOUBLE), intent(in) :: A(LDA,*)
      real(DOUBLE), intent(inout) :: B(LDB,*)
    end subroutine DTRTRS
  end interface

end module LAPACK_MISCEL
