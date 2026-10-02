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

      MODULE NONLINEAR_PARAMS

! Variables used in nonlinear analyses

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG

      IMPLICIT NONE

      SAVE

      CHARACTER(10*BYTE)              :: NL_NORM   = 'INFINITY  '! Method to use when getting the vector norm of UG_COL

      INTEGER(LONG)                   :: LOAD_ISTEP        = 0   ! In nonlinear statics, the load step number
!                                                                  In BUCKLING: 1 is the linear static sol'n,
!                                                                               2 is the BUCKLING eigen step

      INTEGER(LONG)                   :: NL_ITER_NUM       = 0   ! Iteration number within any 1 load step
      INTEGER(LONG)                   :: NL_MAXITER        = 25  ! Max number of iterations per load step
      INTEGER(LONG)                   :: NL_NUM_LOAD_STEPS = 1   ! Number of load steps
      INTEGER(LONG), ALLOCATABLE      :: NL_SID(:)               ! Set ID on the NLPARM entry

      END MODULE NONLINEAR_PARAMS
