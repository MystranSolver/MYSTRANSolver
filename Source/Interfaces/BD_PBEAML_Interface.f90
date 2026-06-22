! ###############################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________
!
! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
! LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
! EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
! _______________________________________________________________________________________________________
!
! End MIT license text.

   MODULE BD_PBEAML_Interface

   INTERFACE

      SUBROUTINE BD_PBEAML ( CARD, LARGE_FLD_INP )


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, IN1
      USE PARAMS, ONLY                :  EPSIL, SUPINFO
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPBEAM, NPBEAM, MPBEAM_STATIONS, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PBEAM, PBEAM_NSTATIONS, PBEAM_XL, PBEAM_RPROPS, RPBEAM

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP



      END SUBROUTINE BD_PBEAML

   END INTERFACE

   END MODULE BD_PBEAML_Interface
