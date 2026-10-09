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

   MODULE SUSER1_RECOVERY

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: SUSER1

   CONTAINS

      SUBROUTINE SUSER1

!  Calc's stresses for user supplied subroutine elements

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  TYPE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SUSER1'





! **********************************************************************************************************************************
      WRITE(ERR,9204) SUBR_NAME,TYPE
      WRITE(F06,9204) SUBR_NAME,TYPE
      FATAL_ERR = FATAL_ERR + 1
      CALL OUTA_HERE ( 'Y' )



      RETURN

! **********************************************************************************************************************************
 9204 FORMAT(' *ERROR  9204: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CODE NOT WRITTEN FOR ELEMENT TYPE = ',A)

! **********************************************************************************************************************************

      END SUBROUTINE SUSER1

   END MODULE SUSER1_RECOVERY
