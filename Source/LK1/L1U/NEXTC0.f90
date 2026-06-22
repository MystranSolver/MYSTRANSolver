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

      SUBROUTINE NEXTC0 ( CARD, ICONT, IERR )

      ! This version of NEXTC is used in BD_xxxx0 routines called by LOADB0
      ! and is the same as NEXTC except that it does not write CARD to F06
      ! under any circumstances
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, IN1, INFILE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE NEXTC0_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'NEXTC0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A MYSTRAN data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! 10 fields of characters of CARD
      CHARACTER(24*BYTE)              :: MESSAG            ! Message for output error purposes
      CHARACTER(3*BYTE)               :: RAW_PARENT        ! 'YES' if parent card came in as raw free-field text
      CHARACTER(LEN(JCARD))           :: NEWTAG            ! Field 1  of cont   card
      CHARACTER(LEN(JCARD))           :: OLDTAG            ! Field 10 of parent card
      CHARACTER(LEN=LEN(CARD))        :: CARD_IN           ! Raw continuation line before FFIELD alters it
      CHARACTER(LEN=LEN(CARD))        :: TCARD             ! Temporary version of CARD

      INTEGER(LONG), INTENT(OUT)      :: ICONT             ! =1 if next card is current card's continuation or =0 if not
      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator from subr FFIELD, called herein
      !INTEGER(LONG)                   :: COMMENT_COL       ! Col on CARD where a comment begins (if one exists)
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from READ
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file. Input to subr READERR




! **********************************************************************************************************************************
      ! Initialize error indicator and ICONT
      IERR  = 0
      ICONT = 0
      RAW_PARENT = 'NO '
      IF ((INDEX(CARD,',') > 0) .OR. (INDEX(CARD,ACHAR(9)) > 0)) RAW_PARENT = 'YES'

      ! Make units for writing errors the error file and output file
      OUNT(1) = ERR
      OUNT(2) = F06

      ! Make JCARD for parent CARD
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      ! Read next card. If it is the cont card for this parent, keep it.
      ! Otherwise backspace the input file.
      OLDTAG = JCARD(10)
      MESSAG = 'BULK DATA CARD          '
      CALL READ_BDF_LINE(IN1, IOCHK, TCARD)
      CARD_IN = TCARD
!      IF (IOCHK /= 0) THEN
!         REC_NO = -99
!         CALL READERR ( IOCHK, INFILE, MESSAG, REC_NO, OUNT, 'Y' )
!         FATAL_ERR = FATAL_ERR + 1
!      ENDIF
!
!      ! Remove any comments within the CARD by deleting everything fro $ on (after col 1)
!      COMMENT_COL = 1
!      DO I=2,BD_ENTRY_LEN
!         IF (TCARD(I:I) == '$') THEN
!            COMMENT_COL = I
!            EXIT
!         ENDIF
!      ENDDO
!
!      IF (COMMENT_COL > 1) THEN
!         TCARD(COMMENT_COL:) = ' '
!      ENDIF


! Make JCARD for TCARD above and get FFIELD to left adjust and fix-field it (if necessary).
!xx CODE COMMENTED OUT IS REPLACED WITH CODE BELOW IT
!xx   IF ((TCARD(1:1) /= '$')  .AND. (TCARD(1:) /= ' ')) THEN
!xx      CALL FFIELD ( TCARD, IERR )
!xx      CALL MKJCARD ( SUBR_NAME, TCARD, JCARD )
!xx      IF (OLDTAG == JCARD(1)) THEN
!xx         ICONT = 1
!xx         CARD = TCARD
!xx      ELSE
!xx         BACKSPACE(IN1)
!xx      ENDIF
!xx   ELSE
!xx      BACKSPACE(IN1)
!xx   ENDIF

      IF (TCARD(1:1) /= '$') THEN
         CALL FFIELD ( TCARD, IERR )
         CALL MKJCARD ( SUBR_NAME, TCARD, JCARD )
         NEWTAG = JCARD(1)
         IF ((RAW_PARENT == 'YES') .AND. (CARD_IN(1:1) == '+')) THEN
            ICONT = 1
         ELSE IF (NEWTAG == OLDTAG) THEN
            ICONT = 1
         ELSE IF ((OLDTAG(1:1) == '+') .AND. (NEWTAG(1:1) == ' ') .AND. (OLDTAG(2:8) == NEWTAG(2:8))) THEN
            ICONT = 1
         ELSE IF ((OLDTAG(1:1) == ' ') .AND. (NEWTAG(1:1) == '+') .AND. (OLDTAG(2:8) == NEWTAG(2:8))) THEN
            ICONT = 1
         ELSE
            BACKSPACE(IN1)
            RETURN
         ENDIF
         IF (RAW_PARENT == 'YES') THEN
            CARD = CARD_IN
         ELSE
            CARD = TCARD
         ENDIF
      ELSE
         BACKSPACE(IN1)
         REC_NO = -99
         CALL READERR (IOCHK, INFILE, MESSAG, REC_NO, OUNT)
         FATAL_ERR = FATAL_ERR + 1
      ENDIF




      RETURN

! **********************************************************************************************************************************
!  101 FORMAT(A)

! **********************************************************************************************************************************

      END SUBROUTINE NEXTC0
