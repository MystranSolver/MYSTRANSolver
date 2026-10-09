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

      MODULE SCRATCH_MATRICES

! Arrays for scratch matrices. CRSi are matrices in sparse CRS format. CCSi are matrices in sparse CCS format

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      INTEGER(LONG), ALLOCATABLE      :: I_CCS1(:)         ! Row indicators for nonzero terms in scratch matrix CCS1
      INTEGER(LONG), ALLOCATABLE      :: I_CCS2(:)         ! Row indicators for nonzero terms in scratch matrix CCS2
      INTEGER(LONG), ALLOCATABLE      :: I_CCS3(:)         ! Row indicators for nonzero terms in scratch matrix CCS3

      INTEGER(LONG), ALLOCATABLE      :: I_CRS1(:)         ! Row indicators for nonzero terms in scratch matrix CRS1
      INTEGER(LONG), ALLOCATABLE      :: I_CRS2(:)         ! Row indicators for nonzero terms in scratch matrix CRS2
      INTEGER(LONG), ALLOCATABLE      :: I_CRS3(:)         ! Row indicators for nonzero terms in scratch matrix CRS3

      INTEGER(LONG), ALLOCATABLE      :: J_CCS1(:)         ! Col numbers for nonzero terms in scratch matrix CCS1
      INTEGER(LONG), ALLOCATABLE      :: J_CCS2(:)         ! Col numbers for nonzero terms in scratch matrix CCS2
      INTEGER(LONG), ALLOCATABLE      :: J_CCS3(:)         ! Col numbers for nonzero terms in scratch matrix CCS3

      INTEGER(LONG), ALLOCATABLE      :: J_CRS1(:)         ! Col numbers for nonzero terms in scratch matrix CRS1
      INTEGER(LONG), ALLOCATABLE      :: J_CRS2(:)         ! Col numbers for nonzero terms in scratch matrix CRS2
      INTEGER(LONG), ALLOCATABLE      :: J_CRS3(:)         ! Col numbers for nonzero terms in scratch matrix CRS3

      REAL(DOUBLE) , ALLOCATABLE      :: CCS1(:)           ! Compressed col storage scratch matrix
      REAL(DOUBLE) , ALLOCATABLE      :: CCS2(:)           ! Compressed col storage scratch matrix
      REAL(DOUBLE) , ALLOCATABLE      :: CCS3(:)           ! Compressed col storage scratch matrix

      REAL(DOUBLE) , ALLOCATABLE      :: CRS1(:)           ! Compressed row storage scratch matrix
      REAL(DOUBLE) , ALLOCATABLE      :: CRS2(:)           ! Compressed row storage scratch matrix
      REAL(DOUBLE) , ALLOCATABLE      :: CRS3(:)           ! Compressed row storage scratch matrix

      END MODULE SCRATCH_MATRICES
