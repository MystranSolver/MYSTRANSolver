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

   MODULE PROGRESS_COUNTERS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: COUNTER_INIT, COUNTER_PROGRESS

   CONTAINS

! This subroutine inits a counter so you can use COUNTER_PROGRESS.
! It's an alternative to way MYSTRAN used to do counters, which badly affect
! performance due to the blocking nature of terminal output.

! Arguments:
!   PREFIX: a string printed before the progress (e.g. "PROCESSING K-SET")
!    TOTAL: the maximum number this counter will reach.
SUBROUTINE COUNTER_INIT(PREFIX, TOTAL)

   USE PENTIUM_II_KIND, ONLY : LONG
   USE SCONTR, ONLY : COUNTER_VALUE, COUNTER_PERC, COUNTER_TOTAL, &
                      COUNTER_STARTED, COUNTER_PREFIX, FATAL_ERR, COUNTER_FMT
   USE IOUNT1, ONLY : ERR, F06
   USE PARAMS, ONLY : NOCOUNTS

      USE UNIX_TIME_MOD, ONLY:  UNIX_TIME
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

   IMPLICIT NONE

   CHARACTER(LEN=*), INTENT(IN) :: PREFIX
   INTEGER(LONG),    INTENT(IN) :: TOTAL

   ! Do nothing if NOCOUNTS is set
   IF (NOCOUNTS == 'Y') THEN
      RETURN
   END IF

   ! Check if the last counter was done
   IF (COUNTER_VALUE /= COUNTER_TOTAL) THEN
      FATAL_ERR = FATAL_ERR + 1
      WRITE(ERR,1414)
      WRITE(F06,1414)
      CALL OUTA_HERE ( 'Y' )
   END IF

   ! Reset everything
   COUNTER_VALUE = 0
   COUNTER_PERC = 0
   COUNTER_TOTAL = TOTAL

   ALLOCATE(CHARACTER(LEN(PREFIX)) :: COUNTER_PREFIX)

   ! Set the timestamp and the prefix
   CALL UNIX_TIME(COUNTER_STARTED)
   COUNTER_PREFIX = PREFIX
   WRITE(COUNTER_FMT, "(A, I0, A)") "(I", CEILING(LOG10(REAL(TOTAL))), ")"

   ! Call the progress subroutine and force printing
   CALL COUNTER_PROGRESS(0)

   1414 FORMAT(' *ERROR  1414: A COUNTER WAS STARTED BEFORE THE CURRENT ONE WAS DONE' &
                        ,/,14X,' THIS IS A DEVELOPER ERROR. PLEASE REPORT IT.')
END SUBROUTINE COUNTER_INIT


! This subroutine increases the progress of a counter by some amount.
! When the counter reaches the max value, it automatically resets itself.
! The counter is only printed when the progress changes, so there is no
! performance impact.

! Arguments:
!   NEW_VALUE_OR_NEGATIVE: the new value for the counter, or a negative for 100% (force-finish)
SUBROUTINE COUNTER_PROGRESS(NEW_VALUE_OR_NEGATIVE)

   USE PENTIUM_II_KIND, ONLY : LONG, DOUBLE
   USE SCONTR, ONLY : COUNTER_VALUE, COUNTER_PERC, COUNTER_TOTAL, &
                      COUNTER_STARTED, COUNTER_PREFIX, COUNTER_FMT, &
                      COUNTER_UPDATED, COUNTER_LIMITER
   USE IOUNT1, ONLY : SC1
   USE CONSTANTS_1, ONLY : ZERO
   USE PARAMS, ONLY : NOCOUNTS

      USE UNIX_TIME_MOD, ONLY:  UNIX_TIME

   IMPLICIT NONE

   INTEGER(LONG), INTENT(IN) :: NEW_VALUE_OR_NEGATIVE
   INTEGER(LONG)             :: NEW_PERC, ETA, NEW_TIME, ELAPSED, &
                                ETA_HOURS, ETA_MINS, ETA_SECS, &
                                ELAPSED_HOURS, ELAPSED_MINS, ELAPSED_SECS, &
                                NEW_VALUE
   REAL(DOUBLE)              :: SPEED
   LOGICAL                   :: FORCED_PRINT, PERC_CHANGED, TOO_MANY

   ! Do nothing if NOCOUNTS is set
   IF (NOCOUNTS == 'Y') THEN
      RETURN
   END IF

   IF (NEW_VALUE_OR_NEGATIVE < 0) THEN
      NEW_VALUE = COUNTER_TOTAL
   ELSE
      NEW_VALUE = NEW_VALUE_OR_NEGATIVE
   END IF

   ! Compute the new percentage
   NEW_PERC = FLOOR(100.0 * NEW_VALUE / COUNTER_TOTAL)

   FORCED_PRINT = NEW_VALUE == 0 .OR. NEW_VALUE == COUNTER_TOTAL
   PERC_CHANGED = NEW_PERC /= COUNTER_PERC
   TOO_MANY = .FALSE.
   ! Check if there has been a change, or if the amount is zero
   IF (FORCED_PRINT .OR. PERC_CHANGED) THEN
      ! Compute elapsed time, 0 means we're too fast
      CALL UNIX_TIME(NEW_TIME)
      ELAPSED = NEW_TIME - COUNTER_STARTED
      ! Zero the counter limiter if the current second has changed
      IF (COUNTER_UPDATED /= NEW_TIME) THEN
         COUNTER_LIMITER = 0
      END IF
      ! Limit to two non-forced updates
      TOO_MANY = COUNTER_LIMITER >= 1
   END IF

   IF (FORCED_PRINT .OR. (PERC_CHANGED .AND. .NOT. TOO_MANY)) THEN
      COUNTER_LIMITER = COUNTER_LIMITER + 1
      COUNTER_UPDATED = NEW_TIME
      ! Okay, percentage change. Let's do this.
      ! First, print the basic string.
      WRITE(SC1, 4000, ADVANCE="NO") CHAR(13)
      !WRITE(SC1, 1000, ADVANCE="NO") COUNTER_PREFIX, NEW_VALUE, COUNTER_TOTAL, NEW_PERC
      WRITE(SC1, 1001, ADVANCE="NO") COUNTER_PREFIX
      WRITE(SC1, COUNTER_FMT, ADVANCE="NO") NEW_VALUE
      WRITE(SC1, 1002, ADVANCE="NO")
      WRITE(SC1, COUNTER_FMT, ADVANCE="NO") COUNTER_TOTAL
      WRITE(SC1, 1003, ADVANCE="NO") NEW_PERC
      ! If the last amount wasn't zero, we can compute speed and ETA!
      IF (NEW_VALUE /= ZERO) THEN
         IF (NEW_VALUE < COUNTER_TOTAL) THEN
            IF (ELAPSED /= 0) THEN
               SPEED = NEW_VALUE / ELAPSED
               ! Write the speed
               WRITE(SC1, 2000, ADVANCE="NO") FLOOR(SPEED)
               ! Compute ETA components
               ETA = FLOOR((COUNTER_TOTAL - NEW_VALUE) / SPEED)
               ETA_HOURS = ETA / 3600
               ETA_MINS = MOD(ETA, 3600) / 60
               ETA_SECS = MOD(MOD(ETA, 3600), 60)
               IF (ETA_HOURS /= ZERO) THEN
                  WRITE(SC1, 3000, ADVANCE="NO") ETA_HOURS, "h"
               END IF
               IF (ETA_MINS /= ZERO) THEN
                  WRITE(SC1, 3000, ADVANCE="NO") ETA_MINS, "min"
               END IF
               WRITE(SC1, 3000, ADVANCE="NO") ETA_SECS, "s"
            END IF
         END IF
      END IF

      IF (NEW_VALUE < COUNTER_TOTAL) THEN
         WRITE(SC1, 4000, ADVANCE="NO") CHAR(13)
      ELSE
         WRITE(SC1, 2500, ADVANCE="NO")
         ELAPSED_HOURS = ELAPSED / 3600
         ELAPSED_MINS = MOD(ELAPSED, 3600) / 60
         ELAPSED_SECS = MOD(MOD(ELAPSED, 3600), 60)
         IF (ELAPSED_HOURS /= ZERO) THEN
            WRITE(SC1, 3000, ADVANCE="NO") ELAPSED_HOURS, "h"
         END IF
         IF (ELAPSED_MINS /= ZERO) THEN
            WRITE(SC1, 3000, ADVANCE="NO") ELAPSED_MINS, "min"
         END IF
         WRITE(SC1, 3000) ELAPSED_SECS, "s               "
         DEALLOCATE(COUNTER_PREFIX)
      END IF

      ! Update the values
      COUNTER_VALUE = NEW_VALUE
      COUNTER_PERC = NEW_PERC
   END IF

   !1000 FORMAT(A, " ", I8, " of ", I8, " (", I3, "%)")
   1001 FORMAT(A, " ")
   1002 FORMAT(" of ")
   1003 FORMAT(" (", I3, "%)")
   2000 FORMAT(" -- ~", I0, "/s, ETA = ")
   2500 FORMAT(" -- done in ~")
   3000 FORMAT(I2, A, " ")
   4000 FORMAT(A)
END SUBROUTINE COUNTER_PROGRESS

   END MODULE PROGRESS_COUNTERS
