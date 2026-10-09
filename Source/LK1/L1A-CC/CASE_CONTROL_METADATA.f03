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

   MODULE CASE_CONTROL_METADATA

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CC_ECHO, CC_LABE, CC_SUBT, CC_TITL

   CONTAINS

      SUBROUTINE CC_ECHO ( CARD )

! Processes Case Control ECHO cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN

      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_ECHO'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=LEN(CARD))        :: CARD1             ! CARD shifted to begin in col after "=" sign

      INTEGER(LONG)                   :: ECOL              ! Col, on CARD, where "=" sign is located
      INTEGER(LONG)                   :: IERR              ! Output from subr CSHIFT indicating an error




! **********************************************************************************************************************************
! Process ECHO card

      CALL CSHIFT ( CARD, '=', CARD1, ECOL, IERR )
      IF (IERR == 0) THEN
         IF      (CARD1(1:4) == 'NONE') THEN
            ECHO = 'NONE  '
         ELSE IF (CARD1(1:4) == 'BOTH') THEN
            ECHO = 'BOTH  '
         ELSE IF (CARD1(1:4) == 'SORT') THEN
            ECHO = 'SORT  '
         ELSE IF (CARD1(1:6) == 'UNSORT') THEN
            ECHO = 'UNSORT'
         ELSE
            ECHO = 'UNSORT'
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,8863) CARD1(1:8)
            IF (SUPWARN == 'N') THEN
               WRITE(F06,8863) CARD1(1:8)
            ENDIF
         ENDIF
      ELSE
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,8876)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,8876)
         ENDIF
         ECHO = 'UNSORT'
      ENDIF



      RETURN

! **********************************************************************************************************************************
 8863 FORMAT(' *WARNING    : ERROR ON CASE CONTROL DECK ECHO CARD.  ENTRY = ',A8,' INCORRECT. DEFAULT UNSORT WILL BE USED')

 8876 FORMAT(' *WARNING    : MISSING EQUAL (=) SIGN ON CASE CONTROL ECHO CARD. DEFAULT UNSORT WILL BE USED')

! **********************************************************************************************************************************

      END SUBROUTINE CC_ECHO


      SUBROUTINE CC_LABE ( CARD )

! Processes Case Control LABEL cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  WARN_ERR, LSUB, NSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  LABEL

      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_LABE'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=LEN(CARD))        :: CARD1             ! CARD shifted to begin in col after "=" sign

      INTEGER(LONG)                   :: ECOL              ! Col, on CARD, where "=" sign is located
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! Output from subr CSHIFT indicating an error




! **********************************************************************************************************************************
! Process LABEL card

      CALL CSHIFT ( CARD, '=', CARD1, ECOL, IERR )
      IF (IERR /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,8862)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,8862)
         ENDIF
      ENDIF
      IF (NSUB /= 0) THEN
         LABEL(NSUB) = CARD1
      ELSE
         DO I = 1,LSUB
            LABEL(I) = CARD1
        ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 8862 FORMAT(' *WARNING    : MISSING EQUAL (=) SIGN ON CASE CONTROL CARD: (TITLE, SUBTITLE OR LABEL).')

! **********************************************************************************************************************************

      END SUBROUTINE CC_LABE


      SUBROUTINE CC_SUBT ( CARD )

! Processes Case Control SUBTITLE cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  WARN_ERR, LSUB, NSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  STITLE

      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_SUBT'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=LEN(CARD))        :: CARD1             ! CARD shifted to begin in col after "=" sign

      INTEGER(LONG)                   :: ECOL              ! Col, on CARD, where "=" sign is located
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! Output from subr CSHIFT indicating an error




! **********************************************************************************************************************************
! Process SUBTITLE cards

      CALL CSHIFT ( CARD, '=', CARD1, ECOL, IERR )
      IF (IERR /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,8862)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,8862)
         ENDIF
      ENDIF
      IF (NSUB /= 0) THEN
         STITLE(NSUB) = CARD1
      ELSE
         DO I = 1,LSUB
            STITLE(I) = CARD1
        ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 8862 FORMAT(' *WARNING    : MISSING EQUAL (=) SIGN ON CASE CONTROL CARD: (TITLE, SUBTITLE OR LABEL).')

! **********************************************************************************************************************************

      END SUBROUTINE CC_SUBT


      SUBROUTINE CC_TITL ( CARD )

! Processes Case Control TITLE cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  WARN_ERR, LSUB, NSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  TITLE

      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_TITL'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=LEN(CARD))        :: CARD1             ! CARD shifted to begin in col after "=" sign

      INTEGER(LONG)                   :: ECOL              ! Col, on CARD, where "=" sign is located
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! Output from subr CSHIFT indicating an error




! **********************************************************************************************************************************
! Process TITLE card

      CALL CSHIFT ( CARD, '=', CARD1, ECOL, IERR )
      IF (IERR /= 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,8862)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,8862)
         ENDIF
      ENDIF
      IF (NSUB /= 0) THEN
         TITLE(NSUB) = CARD1
      ELSE
         DO I = 1,LSUB
            TITLE(I) = CARD1
        ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 8862 FORMAT(' *WARNING    : MISSING EQUAL (=) SIGN ON CASE CONTROL CARD: (TITLE, SUBTITLE OR LABEL).')

! **********************************************************************************************************************************

      END SUBROUTINE CC_TITL

   END MODULE CASE_CONTROL_METADATA
