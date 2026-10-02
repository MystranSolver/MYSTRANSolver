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

      MODULE EMS_ARRAYS

! Initial, linked list, representation of the G-set mass matrix for elem mass. This is the form in which the individual element mass
! matrices are assembled into the system mass matrix in subr EMP.

! See module STF_ARRAYS for an example of this linked list format

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      INTEGER(LONG), ALLOCATABLE      :: EMSCOL(:)         ! See explanation below
      INTEGER(LONG), ALLOCATABLE      :: EMSPNT(:)         ! See explanation below
      INTEGER(LONG), ALLOCATABLE      :: EMSKEY(:)         ! See explanation below

      REAL(DOUBLE) , ALLOCATABLE      :: EMS(:)            ! See explanation below


!  EMS(I)    = 1D real array of the nonzero terms, above the diag, in the G-set elem mass matrix. Size is NTERM_EMG.

!  EMSKEY(I) = 1D integer array of pointers to where, in EMS, the 1st mass term for row I exists.
!              EMSKEY(I) = 0 if row I  is null. Size is NDOFA.

!  EMSPNT(I) = 1D integer array of pointers to where, in EMS, the 2nd and remaining mass terms for row I exist.
!              The 2nd term in row I is at EMSPNT(EMSKEY(I))          in EMS.
!              The 3rd term IN ROW i is at EMSPNT(EMSPNT(EMSKEY(I))), etc.
!              When a zero value of EMSPNT is reached in this scheme it means that there are no more terms in row I.
!              Size is NTERM_EMG.

!  EMSCOL(I) = 1D integer array of the column numbers of the terms in EMS(I). Size is NTERM_EMG.


      END MODULE EMS_ARRAYS
