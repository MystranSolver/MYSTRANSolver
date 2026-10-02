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

      MODULE SPARSE_ALG_ARRAYS

! Arrays used in several of the sparse matrix operation routines (MATADD, MATMULT..., etc) written for MYSTRAN

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      LOGICAL      , ALLOCATABLE      :: LOGICAL_VEC(:)    ! LOGICAL_VEC will be used to show which cols of C have terms in the
!                                                            MATADD_SSS subr

      CHARACTER( 3*BYTE), ALLOCATABLE :: ALG(:)            ! Which algorithm is used (#2 for terms above diag when SYM_A='Y'

      INTEGER(LONG), ALLOCATABLE      :: J_AROW(:)         ! Col numbers of terms in real array AROW (used in sparse mult subrs)

      REAL(DOUBLE) , ALLOCATABLE      :: REAL_VEC(:)       ! REAL_VEC will be used to temp store terms which will be in C in the
!                                                            MATADD_SSS subr

      REAL(DOUBLE) , ALLOCATABLE      :: AROW(:)           ! Array containing the nonzero terms from one row of A
!                                                            (used in sparse mult subrs)

      END MODULE SPARSE_ALG_ARRAYS
