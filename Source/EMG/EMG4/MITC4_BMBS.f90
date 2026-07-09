! #################################################################################################################################
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
SUBROUTINE MITC4_BMBS ( R, S, BM, BB, BS )

  USE PENTIUM_II_KIND, ONLY   : DOUBLE
  USE MODEL_STUF, ONLY    : ELGP, EPROP
  USE CONSTANTS_1, ONLY       : ZERO, ONE

  USE MITC4_B_Interface

  IMPLICIT NONE

  REAL(DOUBLE), INTENT(IN)    :: R, S
  REAL(DOUBLE), INTENT(OUT)   :: BM(3, 6*ELGP)
  REAL(DOUBLE), INTENT(OUT)   :: BB(3, 6*ELGP)
  REAL(DOUBLE), INTENT(OUT)   :: BS(2, 6*ELGP)

  REAL(DOUBLE)        :: BMEM(6, 6*ELGP)
  REAL(DOUBLE)        :: BBOT(6, 6*ELGP)
  REAL(DOUBLE)        :: BTOP(6, 6*ELGP)
  REAL(DOUBLE)        :: BSHR(6, 6*ELGP)

! **********************************************************************************************************************************
! Pure midsurface membrane operator.

  CALL MITC4_B( R, S, ZERO, .TRUE., .FALSE., .FALSE., BMEM )

  BM(1,:) = BMEM(1,:)  ! xx
  BM(2,:) = BMEM(2,:)  ! yy
  BM(3,:) = BMEM(4,:)  ! xy

! **********************************************************************************************************************************
! Pure curvature operator.
! Use bending-only calls so membrane terms cannot leak into curvature.
! The difference removes even-in-T bending terms.

  CALL MITC4_B( R, S, -ONE, .FALSE., .TRUE., .FALSE., BBOT )
  CALL MITC4_B( R, S, +ONE, .FALSE., .TRUE., .FALSE., BTOP )

  BB(1,:) = (BTOP(1,:) - BBOT(1,:)) / EPROP(1)  ! kxx-ish
  BB(2,:) = (BTOP(2,:) - BBOT(2,:)) / EPROP(1)  ! kyy-ish
  BB(3,:) = (BTOP(4,:) - BBOT(4,:)) / EPROP(1)  ! kxy-ish

! **********************************************************************************************************************************
! Pure transverse shear operator.
! MITC shear is already evaluated at tying points and is effectively midsurface shear here.
! Keep row order as zx, yz to match SHELL_T convention used by MITC4.f90.

  CALL MITC4_B( R, S, ZERO, .FALSE., .FALSE., .TRUE., BSHR )

  BS(1,:) = BSHR(6,:)  ! zx
  BS(2,:) = BSHR(5,:)  ! yz

  RETURN

END SUBROUTINE MITC4_BMBS
