module LAPACK_STD_EIG_1

  use PENTIUM_II_KIND, only : BYTE, LONG, DOUBLE

  implicit none

  interface
    subroutine DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
      use PENTIUM_II_KIND, only : LONG, DOUBLE
      implicit none
      character(1), intent(in) :: JOBZ, UPLO
      integer(LONG), intent(in) :: N, LDA, LWORK
      integer(LONG), intent(out) :: INFO
      real(DOUBLE), intent(inout) :: A(LDA,*)
      real(DOUBLE), intent(out) :: W(*)
      real(DOUBLE), intent(inout) :: WORK(*)
    end subroutine DSYEV
  end interface

end module LAPACK_STD_EIG_1
