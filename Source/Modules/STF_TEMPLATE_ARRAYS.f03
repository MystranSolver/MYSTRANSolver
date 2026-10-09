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

      MODULE STF_TEMPLATE_ARRAYS

! Arrays used when the user wants to get a template of the stiffness array written to the F06 file. Not useful for large problens as
! it would produce, possibly, an enormous amount of printed output.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      LOGICAL           , ALLOCATABLE :: TEMPLATE(:,:)     ! Template array that shows where there are nonzero stiff terms in KGG

      CHARACTER( 1*BYTE), ALLOCATABLE :: CROW(:)           ! Char array of 'K' or '_' to print a template of a row of KGG

      END MODULE STF_TEMPLATE_ARRAYS
