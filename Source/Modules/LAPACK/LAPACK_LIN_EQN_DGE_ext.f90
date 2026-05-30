module LAPACK_LIN_EQN_DGE

  use PENTIUM_II_KIND, only : BYTE, LONG, DOUBLE

  implicit none

  interface
    subroutine DGETRF(M, N, A, LDA, IPIV, INFO)
      use PENTIUM_II_KIND, only : LONG, DOUBLE
      implicit none
      integer(LONG), intent(in) :: M, N, LDA
      integer(LONG), intent(out) :: IPIV(*), INFO
      real(DOUBLE), intent(inout) :: A(LDA,*)
    end subroutine DGETRF

    subroutine DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      use PENTIUM_II_KIND, only : LONG, DOUBLE
      implicit none
      character(1), intent(in) :: TRANS
      integer(LONG), intent(in) :: N, NRHS, LDA, LDB
      integer(LONG), intent(in) :: IPIV(*)
      integer(LONG), intent(out) :: INFO
      real(DOUBLE), intent(in) :: A(LDA,*)
      real(DOUBLE), intent(inout) :: B(LDB,*)
    end subroutine DGETRS

    subroutine DGETRI(N, A, LDA, IPIV, WORK, LWORK, INFO)
      use PENTIUM_II_KIND, only : LONG, DOUBLE
      implicit none
      integer(LONG), intent(in) :: N, LDA, LWORK
      integer(LONG), intent(in) :: IPIV(*)
      integer(LONG), intent(out) :: INFO
      real(DOUBLE), intent(inout) :: A(LDA,*), WORK(*)
    end subroutine DGETRI
  end interface

end module LAPACK_LIN_EQN_DGE
