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

      SUBROUTINE NEXTC ( CARD, ICONTINUE, IERR )

      ! Looks for a Bulk Data continuation card belonging to a parent card.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, IN1, INFILE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, ECHO, FATAL_ERR, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE NEXTC_USE_IFs

      IMPLICIT NONE
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'NEXTC'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A MYSTRAN data card
      CHARACTER(LEN=LEN(CARD))        :: CARD_IN           ! Version of CARD read here
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10), JCARD0(10) ! 10 fields of 8 characters of CARD
      CHARACTER(24*BYTE)              :: MESSAG            ! Message for output error purposes
      CHARACTER(1*BYTE)               :: NEWCHAR           ! first character of new line
      CHARACTER(3*BYTE)               :: RAW_PARENT        ! 'YES' if parent card came in as raw free-field text
      CHARACTER(LEN(JCARD))           :: NEWTAG            ! Field 1  of cont   card
      CHARACTER(LEN(JCARD))           :: OLDTAG            ! Field 10 of parent card
      CHARACTER(LEN=LEN(CARD))        :: TCARD             ! Temporary version of CARD

      INTEGER(LONG), INTENT(OUT)      :: ICONTINUE         ! =1 if next card is current card's continuation or =0 if not
      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator from subr FFIELD, called herein
      INTEGER(LONG)                   :: COMMENT_COL       ! Col on CARD where a comment begins (if one exists)
      INTEGER(LONG)                   :: NREADS            ! number of lines read
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from READ
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file. Input to subr READERR




! **********************************************************************************************************************************
      ! Initialize error indicator
      IERR = 0
      ICONTINUE = 0
      NEWCHAR = ' '
      RAW_PARENT = 'NO '
      IF ((INDEX(CARD,',') > 0) .OR. (INDEX(CARD,ACHAR(9)) > 0)) RAW_PARENT = 'YES'

      ! Make units for writing errors the error file and output file
      OUNT(1) = ERR
      OUNT(2) = F06

      ! Make JCARD for parent CARD
      ! split the line (CARD) into fields (JCARD)
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      ! copy jcard to jcard0 in case we have an error
      DO I=1,10
          JCARD0(I) = JCARD(I)
      ENDDO

      !------------------------
      ! Read next card.

      ! OLDTAG is field 10 of the current card coming into this subr
      OLDTAG = JCARD(10)

      ! Read next card
      CALL READ_BDF_LINE(IN1, IOCHK, TCARD)
      CARD_IN = TCARD
      NEWCHAR = TCARD(1:1)
      !
      ! Make JCARD for TCARD above and get FFIELD to left adjust and
      ! fix-field it (if necessary).
      !
      ! get a flag (ICONTINUE) that defines if we need to read the next line
      ! we do because we have OLDTAG

      ! we know know that the line doesn't start with a $, but it
      ! can be a different card (OLDTAG /= NEWTAG)
      !
      CALL FFIELD ( TCARD, IERR )

      ! split the continuation line (TCARD) into fields (JCARD)
      CALL MKJCARD ( SUBR_NAME, TCARD, JCARD )
      NEWTAG = JCARD(1)

      ! do i need to flag this?
      ! (NEWTAG(1:1) /= '*')
      IF ((NEWTAG(1:1) /= ' ') .AND. (NEWTAG(1:1) /= '+') .AND. (NEWTAG(1:1) /= '$')) THEN
         ! different card type (e.g., LOAD -> FORCE
         CARD = CARD_IN
         BACKSPACE(IN1)
         return

      ELSE IF ((RAW_PARENT == 'YES') .AND. (CARD_IN(1:1) == '+')) THEN
         ICONTINUE = 1
      ELSE IF (NEWTAG == OLDTAG) THEN
         ICONTINUE = 1
      ELSE IF ((OLDTAG(1:1) == '+') .AND. (NEWTAG(1:1) == ' ') .AND. (OLDTAG(2:8) == NEWTAG(2:8))) THEN
         ! small field
         ICONTINUE = 1
      ELSE IF ((OLDTAG(1:1) == ' ') .AND. (NEWTAG(1:1) == '+') .AND. (OLDTAG(2:8) == NEWTAG(2:8))) THEN
         ! small field
         ICONTINUE = 1
      ELSE IF ((NEWTAG(1:1) /= ' ') .AND. (NEWTAG(1:1) /= '+') .AND. (NEWTAG(1:1) /= '$')) THEN
         ! different card type (e.g., LOAD -> FORCE
         BACKSPACE(IN1)
         CARD = TCARD
         RETURN
      ELSE
         BACKSPACE(IN1)
         IF (RAW_PARENT == 'YES') THEN
            RETURN
         ENDIF
         IF (((JCARD0(1)(1:4) == 'CBAR') .OR. (JCARD0(1)(1:5) == 'CBEAM')) .AND.                                                    &
             (INDEX('0123456789',OLDTAG(1:1)) > 0)) THEN
            RETURN
         ENDIF
         ! can't find the continuation marker.  FATAL :)
         WRITE(F06,102) OLDTAG
         WRITE(ERR,102) OLDTAG
         WRITE(F06,103)
         WRITE(ERR,103)
         WRITE(F06,104) 'FIELDS1:', JCARD0
         WRITE(ERR,104) 'FIELDS1:', JCARD0
         WRITE(F06,104) 'FIELDS2:', JCARD
         WRITE(ERR,104) 'FIELDS2:', JCARD
         FLUSH(F06)
         FLUSH(ERR)
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE('Y')  ! FATAL error
         RETURN
      ENDIF
      IF (RAW_PARENT == 'YES') THEN
         CARD = CARD_IN
      ELSE
         CARD = TCARD
      ENDIF
      IF (ECHO(1:4) /= 'NONE') THEN
          WRITE(F06, 101) CARD_IN
      ENDIF
      FLUSH(ERR)



      RETURN

! **********************************************************************************************************************************
  101 FORMAT('ECHO:  nextc ', A)

  ! missing continuation fatal
  102 FORMAT(' *FATAL: COULD NOT FIND A CONTINUATION MARKER FOR ', A)
  103 FORMAT('         CHECK THAT THE CONTINUATION STARTS AT THE BEGINNING OF FIELD 10')

  ! missing continuation; name and 10 fields
  104 FORMAT('         ', A, ' 1: ', A, ' 2:', A, ' 3:', A, ' 4:', A, &
             ' 5:', A, ' 6:', A, ' 7:', A, ' 8:', A, ' 9:', A, ' 10:', A)
! **********************************************************************************************************************************

      END SUBROUTINE NEXTC
