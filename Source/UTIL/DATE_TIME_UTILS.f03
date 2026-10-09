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

   MODULE DATE_TIME_UTILS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: OURDAT, OURTIM, TIME_INIT

   CONTAINS

      SUBROUTINE OURDAT

! Returns date info using Fortran DATE_AND_TIME intrinsic procedure

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE TIMDAT, ONLY                :  YEAR, MONTH, DAY

      IMPLICIT NONE

      CHARACTER( 8*BYTE)              :: DATE              ! Date returned from intrinsic function DATE_AND_TIME
      CHARACTER(10*BYTE)              :: TIME              ! Time returned from intrinsic function DATE_AND_TIME
      CHARACTER( 5*BYTE)              :: ZONE              ! Zone returned from intrinsic function DATE_AND_TIME

      INTEGER(LONG)                   :: VALUES(8)         ! Contains year, month, day from intrinsic function DATE_AND_TIME

      INTRINSIC                       :: DATE_AND_TIME

! **********************************************************************************************************************************
      CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)

      YEAR  = VALUES(1)
      MONTH = VALUES(2)
      DAY   = VALUES(3)

      RETURN

      END SUBROUTINE OURDAT


      SUBROUTINE OURTIM

! Returns time using Fortran DATE_AND_TIME intrinsic procedure

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC, DSEC

      IMPLICIT NONE

      CHARACTER( 8*BYTE)              :: DATE              ! Date returned from intrinsic function DATE_AND_TIME
      CHARACTER(10*BYTE)              :: TIME              ! Time returned from intrinsic function DATE_AND_TIME
      CHARACTER( 5*BYTE)              :: ZONE              ! Zone returned from intrinsic function DATE_AND_TIME

      INTEGER(LONG)                   :: OLD_HOUR          ! Value of HOUR   (module TIMDAT) when OURTIM was called
      INTEGER(LONG)                   :: OLD_MIN           ! Value of MINUTE (module TIMDAT) when OURTIM was called
      INTEGER(LONG)                   :: OLD_SEC           ! Value of SEC    (module TIMDAT) when OURTIM was called
      INTEGER(LONG)                   :: OLD_SFRAC         ! Value of SFRAC  (module TIMDAT) when OURTIM was called
      INTEGER(LONG)                   :: VALUES(8)         ! Contains year, month, day from intrinsic function DATE_AND_TIME

      INTRINSIC                       :: DATE_AND_TIME

! **********************************************************************************************************************************
      OLD_HOUR  = HOUR
      OLD_MIN   = MINUTE
      OLD_SEC   = SEC
      OLD_SFRAC = SFRAC

      CALL DATE_AND_TIME ( DATE, TIME, ZONE, VALUES )

      HOUR   = VALUES(5)
      MINUTE = VALUES(6)
      SEC    = VALUES(7)
      SFRAC  = VALUES(8)

! Reset OLD_HOUR for new day:

      IF (OLD_HOUR > HOUR) THEN
         OLD_HOUR = OLD_HOUR - 24
      ENDIF

! Check for first pass thru OURTIM (HOUR is initialized to -100 in the MYSTRAN.for) and calc DSEC, TSEC (module TIMDAT)

      IF (OLD_HOUR == -100) THEN
         DSEC = 0.0
         TSEC = 0.0
      ELSE
         DSEC = 3600*(HOUR-OLD_HOUR) + 60*(MINUTE-OLD_MIN)+(SEC-OLD_SEC) + .001*(SFRAC-OLD_SFRAC)
         TSEC = TSEC + DSEC
      ENDIF

      RETURN

      END SUBROUTINE OURTIM


      SUBROUTINE TIME_INIT

! Set time initializing parameters

      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC

      IMPLICIT NONE

! **********************************************************************************************************************************
! Set time initializing parameters

      HOUR   = -100                                        ! Used in subr OURTIM so subr begin times start, in a LINK, at 0.
      MINUTE =    0
      SEC    =    0
      SFRAC  =    0

! **********************************************************************************************************************************

      END SUBROUTINE TIME_INIT

   END MODULE DATE_TIME_UTILS
