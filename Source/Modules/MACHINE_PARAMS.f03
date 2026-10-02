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

      MODULE MACHINE_PARAMS

! These are parameters that will be set by a LAPACK routine that knows the capability of the processor of the users computer

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      SAVE

      REAL(DOUBLE)                    :: MACH_BASE         ! Base of the machine
      REAL(DOUBLE)                    :: MACH_EMAX         ! Largest exponent before overflow
      REAL(DOUBLE)                    :: MACH_EMIN         ! Minimum exponent before (gradual) underflow
      REAL(DOUBLE)                    :: MACH_EPS          ! Relative machine precision
      REAL(DOUBLE)                    :: MACH_LARGE_NUM    ! 1./MACH_SFMIN
      REAL(DOUBLE)                    :: MACH_PREC         ! eps*base
      REAL(DOUBLE)                    :: MACH_RMAX         ! Overflow threshold  - (base**emax)*(1-eps)
      REAL(DOUBLE)                    :: MACH_RMIN         ! Underflow threshold - base**(emin-1)
      REAL(DOUBLE)                    :: MACH_RND          ! 1.0 when rounding occurs in addition, 0.0 otherwise
      REAL(DOUBLE)                    :: MACH_SFMIN        ! Safe minimum, such that 1/sfmin does not overflow
      REAL(DOUBLE)                    :: MACH_T            ! Number of (base) digits in the mantissa

      END MODULE MACHINE_PARAMS
