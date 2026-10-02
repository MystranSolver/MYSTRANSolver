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

      MODULE RIGID_BODY_DISP_MATS

! Rigid body displacement matrices

      USE PENTIUM_II_KIND, ONLY :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      REAL(DOUBLE), ALLOCATABLE       :: RBGLOBAL_GSET(:,:)! NDOFG x 6 matrix of displacements for the G-set
      REAL(DOUBLE), ALLOCATABLE       :: RBGLOBAL_NSET(:,:)! NDOFN x 6 matrix of displacements for the N-set
      REAL(DOUBLE), ALLOCATABLE       :: RBGLOBAL_FSET(:,:)! NDOFF x 6 matrix of displacements for the F-set
      REAL(DOUBLE), ALLOCATABLE       :: RBGLOBAL_ASET(:,:)! NDOFA x 6 matrix of displacements for the A-set
      REAL(DOUBLE), ALLOCATABLE       :: RBGLOBAL_LSET(:,:)! NDOFL x 6 matrix of displacements for the L-set

      REAL(DOUBLE), ALLOCATABLE       :: TR6_CG(:,:)       ! NDOFR x 6 rigid body displ matrix for R-set rel to model CG
      REAL(DOUBLE), ALLOCATABLE       :: TR6_MEFM(:,:)     ! NDOFR x 6 rigid body displ matrix for R-set rel to MEFMGRID
      REAL(DOUBLE), ALLOCATABLE       :: TR6_0(:,:)        ! NDOFR x 6 rigid body displ matrix for R-set rel to model basic coords

      END MODULE RIGID_BODY_DISP_MATS
