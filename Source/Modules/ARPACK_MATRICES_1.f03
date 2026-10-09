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

      MODULE ARPACK_MATRICES_1

! Matrices used by ARPACK for Lanczos method.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      LOGICAL      , ALLOCATABLE      :: SELECT(:)         ! See ARPACK subroutine DSBAND for definition

      INTEGER(LONG), ALLOCATABLE      :: IWORK(:)          ! See ARPACK subroutine DSBAND for definition

      REAL(DOUBLE) , ALLOCATABLE      :: RESID(:)          ! See ARPACK subroutine DSBAND for definition
      REAL(DOUBLE) , ALLOCATABLE      :: RFAC(:,:)         ! See ARPACK subroutine DSBAND for definition
      REAL(DOUBLE) , ALLOCATABLE      :: VBAS(:,:)         ! See ARPACK subroutine DSBAND for definition
      REAL(DOUBLE) , ALLOCATABLE      :: WORKD(:)          ! See ARPACK subroutine DSBAND for definition
      REAL(DOUBLE) , ALLOCATABLE      :: WORKL(:)          ! See ARPACK subroutine DSBAND for definition

      END MODULE ARPACK_MATRICES_1
