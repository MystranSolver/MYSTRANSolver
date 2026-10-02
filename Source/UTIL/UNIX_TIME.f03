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

   MODULE UNIX_TIME_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: UNIX_TIME

   CONTAINS

! Computes the current UNIX timestamp via C's time().

SUBROUTINE UNIX_TIME(T)

   USE PENTIUM_II_KIND, ONLY: LONG
   USE ISO_C_BINDING,   ONLY: C_PTR, C_NULL_PTR, C_LONG

   IMPLICIT NONE

   INTEGER(LONG), INTENT(OUT) :: T

   INTERFACE
      FUNCTION C_TIME(TLOC) BIND(C, NAME='time')
         IMPORT :: C_PTR, C_LONG
         TYPE(C_PTR), VALUE  :: TLOC
         INTEGER(C_LONG)     :: C_TIME
      END FUNCTION C_TIME
   END INTERFACE

   T = INT(C_TIME(C_NULL_PTR), LONG)

END SUBROUTINE UNIX_TIME

   END MODULE UNIX_TIME_MOD
