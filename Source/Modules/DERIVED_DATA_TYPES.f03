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

MODULE DERIVED_DATA_TYPES

   USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

   IMPLICIT NONE

   TYPE                            :: INT2_REAL1
      INTEGER(LONG)                ::    Col_1
      INTEGER(LONG)                ::    Col_2
      REAL(DOUBLE)                 ::    Col_3
   END TYPE INT2_REAL1

   TYPE                            :: CHAR1_INT1
      CHARACTER(8*BYTE)            ::    Col_1
      INTEGER(LONG)                ::    Col_2
   END TYPE CHAR1_INT1

END MODULE DERIVED_DATA_TYPES
