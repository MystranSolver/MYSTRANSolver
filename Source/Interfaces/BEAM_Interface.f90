!--- cbeam_add begin ---!
! ###############################################################################################################################
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

   MODULE BEAM_Interface

   INTERFACE

      SUBROUTINE BEAM ( OPT, L, AREA, I1, I2, JTOR, CW, SCOEFF, K1, K2, I12, E, G, ALPHA, TREF )


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      IMPLICIT NONE

      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)

      REAL(DOUBLE), INTENT(IN)        :: ALPHA
      REAL(DOUBLE), INTENT(IN)        :: AREA
      REAL(DOUBLE), INTENT(IN)        :: CW
      REAL(DOUBLE), INTENT(IN)        :: E
      REAL(DOUBLE), INTENT(IN)        :: G
      REAL(DOUBLE), INTENT(IN)        :: I1
      REAL(DOUBLE), INTENT(IN)        :: I12
      REAL(DOUBLE), INTENT(IN)        :: I2
      REAL(DOUBLE), INTENT(IN)        :: JTOR
      REAL(DOUBLE), INTENT(IN)        :: K1
      REAL(DOUBLE), INTENT(IN)        :: K2
      REAL(DOUBLE), INTENT(IN)        :: L
      REAL(DOUBLE), INTENT(IN)        :: SCOEFF
      REAL(DOUBLE), INTENT(IN)        :: TREF

      END SUBROUTINE BEAM

   END INTERFACE

   END MODULE BEAM_Interface

!--- cbeam_add end ---!
