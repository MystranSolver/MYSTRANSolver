! ##################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________

! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.

      MODULE EIGEN_MATRICES_1

! Matrices used in eigenvalue analyses

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      INTEGER(LONG), ALLOCATABLE      :: MODE_NUM(:)       ! Integer array of mode numbers
      INTEGER(LONG), ALLOCATABLE      :: MODES_OUT(:)      ! 0, 1 indicator of whether an eigenvec is requested for output

      REAL(DOUBLE) , ALLOCATABLE      :: EIGEN_VAL(:)      ! Matrix of eigenvalues
      REAL(DOUBLE) , ALLOCATABLE      :: EIGEN_VEC(:,:)    ! Matrix of eigenvectors
      REAL(DOUBLE) , ALLOCATABLE      :: GEN_MASS(:)       ! 1D array of gen masses; one for each eigenvalue
      REAL(DOUBLE) , ALLOCATABLE      :: MEFFMASS(:,:)     ! Modal effective masses
      REAL(DOUBLE) , ALLOCATABLE      :: MPFACTOR_N6(:,:)  ! Modal participation factors rel to 6 displ comps at MEFMGRID
      REAL(DOUBLE) , ALLOCATABLE      :: MPFACTOR_NR(:,:)  ! Modal participation factors for all NDOFR DOF's

      END MODULE EIGEN_MATRICES_1
