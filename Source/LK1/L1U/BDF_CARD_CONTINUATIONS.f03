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

   MODULE BDF_CARD_CONTINUATIONS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: NEXTC, NEXTC0, NEXTC2, NEXTC20, MKCARD, MKJCARD, MKJCARD_08, TAGS_MATCH, LINE_TAG

   CONTAINS

      SUBROUTINE NEXTC ( CARD, ICONTINUE, IERR )

      ! Looks for a Bulk Data continuation card belonging to a parent card.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, IN1, INFILE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, ECHO, FATAL_ERR, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY        :  READERR

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      IMPLICIT NONE
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'NEXTC'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A MYSTRAN data card
      CHARACTER(LEN=LEN(CARD))        :: CARD_IN           ! Version of CARD read here
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10), JCARD0(10) ! 10 fields of 8 characters of CARD
      CHARACTER(24*BYTE)              :: MESSAG            ! Message for output error purposes
      CHARACTER(1*BYTE)               :: NEWCHAR           ! first character of new line
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

      ! A large field continuation ('*' and the rest of the marker in field 10 of CARD) may follow a small field entry: NEXTC2
      ! reads its two physical lines
      IF ((TCARD(1:1) == '*') .AND. TAGS_MATCH ( OLDTAG, TCARD(1:8) )) THEN
         BACKSPACE(IN1)
         CALL NEXTC2 ( CARD, ICONTINUE, IERR, TCARD )
         CARD = TCARD
         RETURN
      ENDIF
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

      ELSE IF (TAGS_MATCH ( OLDTAG, NEWTAG )) THEN       ! The first characters may differ: ' ', '+' or '*'
         ICONTINUE = 1
      ELSE IF ((NEWTAG(1:1) /= ' ') .AND. (NEWTAG(1:1) /= '+') .AND. (NEWTAG(1:1) /= '$')) THEN
         ! different card type (e.g., LOAD -> FORCE
         BACKSPACE(IN1)
         CARD = TCARD
         RETURN
      ELSE
         ! can't find the continuation marker.  FATAL :)
         BACKSPACE(IN1)
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
      CARD = TCARD
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


      SUBROUTINE NEXTC0 ( CARD, ICONT, IERR )

      ! This version of NEXTC is used in BD_xxxx0 routines called by LOADB0
      ! and is the same as NEXTC except that it does not write CARD to F06
      ! under any circumstances
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, IN1, INFILE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY        :  READERR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'NEXTC0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A MYSTRAN data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! 10 fields of characters of CARD
      CHARACTER(24*BYTE)              :: MESSAG            ! Message for output error purposes
      CHARACTER(LEN(JCARD))           :: NEWTAG            ! Field 1  of cont   card
      CHARACTER(LEN(JCARD))           :: OLDTAG            ! Field 10 of parent card
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

      ! The same continuations as NEXTC reads, so that this first pass counts what the second pass reads: a large field
      ! continuation ('*' line) of a small field entry is read by NEXTC20, and the markers match as in TAGS_MATCH
      IF ((TCARD(1:1) == '*') .AND. TAGS_MATCH ( OLDTAG, TCARD(1:8) )) THEN
         BACKSPACE(IN1)
         CALL NEXTC20 ( CARD, ICONT, IERR, TCARD )
         CARD = TCARD
         RETURN
      ENDIF

      IF (TCARD(1:1) /= '$') THEN
         CALL FFIELD ( TCARD, IERR )
         CALL MKJCARD ( SUBR_NAME, TCARD, JCARD )
         NEWTAG = JCARD(1)
         IF ((NEWTAG == OLDTAG) .OR. ((INDEX(' +', NEWTAG(1:1)) > 0) .AND. TAGS_MATCH ( OLDTAG, NEWTAG ))) THEN
            ICONT = 1
         ELSE
            BACKSPACE(IN1)
            RETURN
         ENDIF
         CARD = TCARD
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


      SUBROUTINE NEXTC2 ( PARENT, ICONTINUE, IERR, CHILD )

      ! Looks for 2 physical Bulk Data large field format continuation
      ! entries belonging to a large field parent.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, IN1, INFILE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, ECHO, FATAL_ERR, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY        :  READERR

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'NEXTC2'
      CHARACTER(LEN=*), INTENT(IN)    :: PARENT            !

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(OUT) :: CHILD    !

      CHARACTER(LEN=BD_ENTRY_LEN)     :: CHILD1            !
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CHILD2            !
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! 10 fields of 8 characters of PARENT
      CHARACTER(LEN(JCARD))           :: NEWTAG            ! Field 1  of cont   card
      CHARACTER(LEN(JCARD))           :: OLDTAG            ! Field 10 of parent card

      INTEGER(LONG), INTENT(OUT)      :: ICONTINUE         ! =1 if next card is current card's continuation or =0 if not
      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator from subr FFIELD, called herein
      INTEGER(LONG)                   :: COMMENT_COL       ! Col on PARENT where a comment begins (if one exists)
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from READ
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file. Input to subr READERR


! **********************************************************************************************************************************
      ! Initialize
      IERR = 0
      ICONTINUE = 0
      CHILD(1:) = 'z'
      OLDTAG(1:) = ' '

      IERR = 0

      ! Make units for writing errors the error file and output file
      OUNT(1) = ERR
      OUNT(2) = F06

      ! Make JCARD for PARENT and get the 1st 8 chars of field 10 (cont mnemonic)
      CALL MKJCARD ( SUBR_NAME, PARENT, JCARD )
      OLDTAG(1:8) = JCARD(10)(1:8)

      ! Read next card. If it is a continuation to the parent it will
      ! be the 1st half of the whole continuation
      CALL READ_BDF_LINE(IN1, IOCHK, CHILD1)

      NEWTAG = LINE_TAG ( CHILD1 )                          ! Field 1, also of a free field line
      ! A continuation: its field 1 is the parent's field 10 except for the first character, which is blank, '+' (small field)
      ! or '*' (large field) in either; the parent and its continuations may mix small and large field lines
      IF (TAGS_MATCH ( OLDTAG, NEWTAG )) THEN
         ICONTINUE = 1
      ELSE IF ((NEWTAG(1:1) == ' ') .OR. (NEWTAG(1:1) == '+')) THEN
         ! A continuation line whose marker is not the one in field 10 before it: it continues no entry (as in NEXTC)
         BACKSPACE(IN1)
         WRITE(F06,102) OLDTAG(1:8), NEWTAG(1:8)
         WRITE(ERR,102) OLDTAG(1:8), NEWTAG(1:8)
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE('Y')
         RETURN
      ELSE IF ((NEWTAG(1:1) /= '*') .AND. (NEWTAG(1:1) /= '$')) THEN
         ! different card type (e.g., LOAD -> FORCE
         BACKSPACE(IN1)
         CHILD = CHILD1
         RETURN
      ELSE
         ! can't find the continuation marker.  FATAL :)
         BACKSPACE(IN1)
         !WRITE(F06,102) OLDTAG
         !WRITE(ERR,102) OLDTAG
         !WRITE(F06,103)
         !WRITE(ERR,103)
         !WRITE(F06,104) 'CHILD1:', CHILD1
         !WRITE(ERR,104) 'CHILD1:', CHILD1
         FLUSH(F06)
         FLUSH(ERR)
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE('Y')  ! FATAL error
         RETURN
         RETURN
      ENDIF

      ! A small field continuation (field 1 blank or starting with '+') may follow a large field entry: it is one physical
      ! line with 8 character fields
      IF (CHILD1(1:1) /= '*') THEN
         CALL FFIELD ( CHILD1, IERR )
         CHILD = CHILD1
         IF (ECHO(1:4) /= 'NONE') THEN
            WRITE(F06,101) CHILD1
         ENDIF
         RETURN
      ENDIF

      ! Read 2nd half of continuation entry, if it exists
      CALL READ_BDF_LINE(IN1, IOCHK, CHILD2)

      OLDTAG = CHILD1(73:80)
      NEWTAG = CHILD2( 1: 8)
      ICONTINUE = 0


     ! The second half of a large field continuation starts with '*'; a line starting with a blank or '+' is the next
     ! (small field) continuation, read by the next call
     IF ((NEWTAG(1:1) == '*') .AND. TAGS_MATCH ( OLDTAG, NEWTAG )) THEN
        ICONTINUE = 1
     ELSE
        BACKSPACE(IN1)
        CHILD2(1:) = ' '
        ICONTINUE = 1
        CALL FFIELD2 ( CHILD1, CHILD2, CHILD, IERR )
        RETURN
     ENDIF

      ! Call FFIELD2 to put the 2 CHILDi's together and left justify
      CALL FFIELD2(CHILD1, CHILD2, CHILD, IERR)
      ICONTINUE = 1
      IF (ECHO(1:4) /= 'NONE') THEN
         WRITE(F06,101) CHILD1
         WRITE(F06,101) CHILD2
      ENDIF
! **********************************************************************************************************************************
  101 FORMAT('ECHO nextc2: ', A)
  102 FORMAT(' *FATAL: THE CONTINUATION MARKER "',A,'" IN FIELD 10 IS NOT CONTINUED BY THE NEXT LINE, WHOSE FIELD 1 IS "',A,'"')

! **********************************************************************************************************************************


      END SUBROUTINE NEXTC2


      SUBROUTINE NEXTC20 ( PARENT, ICONT, IERR, CHILD )

      ! Looks for 2 physical Bulk Data large field format continuation
      ! entries belonging to a large field parent.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, IN1, INFILE
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, ECHO, FATAL_ERR, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY        :  READERR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'NEXTC20'
      CHARACTER(LEN=*), INTENT(IN)    :: PARENT            !

      CHARACTER(LEN=BD_ENTRY_LEN), INTENT(OUT) :: CHILD    !

      CHARACTER(LEN=BD_ENTRY_LEN)     :: CHILD1            !
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CHILD2            !
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! 10 fields of 8 characters of PARENT
      CHARACTER(LEN(JCARD))           :: NEWTAG            ! Field 10 of cont   card
      CHARACTER(LEN(JCARD))           :: OLDTAG            ! Field 10 of parent card

      INTEGER(LONG), INTENT(OUT)      :: ICONT             ! =1 if next card is current card's continuation or =0 if not
      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator from subr FFIELD, called herein
      INTEGER(LONG)                   :: COMMENT_COL       ! Col on PARENT where a comment begins (if one exists)
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from READ
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file. Input to subr READERR


! **********************************************************************************************************************************
      ! Initialize
      CHILD(1:) = 'z'
      OLDTAG(1:) = ' '

      IERR  = 0
      ICONT = 0

      ! Make units for writing errors the error file and output file
      OUNT(1) = ERR
      OUNT(2) = F06

      ! Make JCARD for PARENT and get the 1st 8 chars of field 10 (cont mnemonic)
      CALL MKJCARD ( SUBR_NAME, PARENT, JCARD )
      OLDTAG(1:8) = JCARD(10)(1:8)

      ! Read next card. If it is a continuation to the parent it will be
      ! the 1st half of the whole continuation
      CALL READ_BDF_LINE(IN1, IOCHK, CHILD1)
      IF (IOCHK /= 0) THEN
         REC_NO = -99
         CALL READERR ( IOCHK, INFILE, 'BULK DATA CARD', REC_NO, OUNT )
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF
      NEWTAG = LINE_TAG ( CHILD1 )                          ! Field 1, also of a free field line

!xx CODE COMMENTED OUT IS REPLACED WITH CODE BELOW IT
!xx   IF (NEWTAG == OLDTAG) THEN
!xx      ICONT = 1
!xx   ELSE IF ((OLDTAG(1:8) == '*       ') .AND. (NEWTAG(1:8) == '       ')) THEN
!xx      ICONT = 1
!xx   ELSE IF ((OLDTAG(1:8) == '        ') .AND. (NEWTAG(1:8) == '*      ')) THEN
!xx      ICONT = 1
!xx   ELSE
!xx      BACKSPACE(IN1)
!xx      RETURN
!xx   ENDIF

      ! The same continuations as NEXTC2 reads (this first pass must count what the second pass reads)
      IF ((NEWTAG == OLDTAG) .OR. TAGS_MATCH ( OLDTAG, NEWTAG )) THEN
         ICONT = 1
      ELSE
         BACKSPACE(IN1)
         RETURN
      ENDIF

      ! A small field continuation (field 1 blank or starting with '+'): one line of 8 character fields
      IF (CHILD1(1:1) /= '*') THEN
         CALL FFIELD ( CHILD1, IERR )
         CHILD = CHILD1
         RETURN
      ENDIF

      ! Read 2nd half of continuation entry, if it exists
      CALL READ_BDF_LINE(IN1, IOCHK, CHILD2)
      IF (IOCHK /= 0) THEN
         REC_NO = -99
         CALL READERR ( IOCHK, INFILE, 'BULK DATA CARD', REC_NO, OUNT )
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF
      OLDTAG = CHILD1(73:80)
      NEWTAG = CHILD2( 1: 8)

      ! The second half starts with '*'; otherwise the first half is the whole continuation (as in NEXTC2)
      IF (.NOT. ((NEWTAG(1:1) == '*') .AND. TAGS_MATCH ( OLDTAG, NEWTAG ))) THEN
         BACKSPACE(IN1)
         CHILD2(1:) = ' '
         CALL FFIELD2 ( CHILD1, CHILD2, CHILD, IERR )
         ICONT = 1
         RETURN
      ENDIF


      ! Remove any comments within CHILD2 by deleting everything from $ on
      ! (after col 1).
      !
      ! NOTE: CHILD1 cannot have  comments since the last field is used for
      ! NEWTAG above
      COMMENT_COL = 1
      DO I=2,BD_ENTRY_LEN
         IF (CHILD2(I:I) == '$') THEN
            COMMENT_COL = I
            EXIT
         ENDIF
      ENDDO

      IF (COMMENT_COL > 1) THEN
         CHILD2(COMMENT_COL:) = ' '
      ENDIF

      ! Call FFIELD2 to put the 2 CHILDi's together and left justify
      CALL FFIELD2 ( CHILD1, CHILD2, CHILD, IERR )
      ICONT = 1

! **********************************************************************************************************************************
  101 FORMAT(A)

! **********************************************************************************************************************************

      END SUBROUTINE NEXTC20



      SUBROUTINE MKCARD ( JCARD, CARD )

! Routine to create CARD from the 10 CHAR input JCARD fields

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  JCARD_LEN

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)            :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN=10*JCARD_LEN), INTENT(OUT):: CARD              ! A MYSTRAN data card

      INTEGER(LONG)                           :: I                 ! DO loop index
      INTEGER(LONG)                           :: K1,K2             ! Range for setting CARD = JCARD

! **********************************************************************************************************************************
! Initialize outputs

      CARD(1:) = ' '

      DO I=1,10
         K1 = JCARD_LEN*(I-1) + 1
         K2 = K1 + JCARD_LEN - 1
         CARD(K1:K2) = JCARD(I)(1:JCARD_LEN)
      ENDDO

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MKCARD


      SUBROUTINE MKJCARD ( CALLING_SUBR, CARD, JCARD )

! Routine to create JCARD, a set of 10 CHAR fields (each of size JCARD_LEN), from input char CARD

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  JCARD_LEN

      IMPLICIT NONE

      CHARACTER(LEN=*)  , INTENT(IN)       :: CALLING_SUBR      ! Subr that called this one
      CHARACTER(LEN=*)  , INTENT(IN)       :: CARD              ! A MYSTRAN data card
      CHARACTER(LEN=JCARD_LEN), INTENT(OUT):: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                        :: I
      INTEGER(LONG)                        :: K1,K2             ! Range for setting CARD = JCARD

! **********************************************************************************************************************************
      DO I=1,10
         K1 = JCARD_LEN*(I-1) + 1
         K2 = K1 + JCARD_LEN - 1
         JCARD(I)(1:JCARD_LEN) = CARD(K1:K2)
      ENDDO
      RETURN

! **********************************************************************************************************************************

! **********************************************************************************************************************************

      END SUBROUTINE MKJCARD


      SUBROUTINE MKJCARD_08 ( CARD, JCARD_08 )

! Routine to create JCARD_08, a set of 10 CHAR fields (each of size 8 chars), from input char CARD

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG

!     USE MKJCARD_08_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=*)  , INTENT(IN)  :: CARD              ! A MYSTRAN data card
      CHARACTER( 8*BYTE), INTENT(OUT) :: JCARD_08(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: K1,K2             ! Range for setting CARD = JCARD_08

! **********************************************************************************************************************************
      DO I=1,10
         K1 = 8*(I-1) + 1
         K2 = K1 + 7
         JCARD_08(I)(1:8) = CARD(K1:K2)
      ENDDO

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MKJCARD_08

! ##################################################################################################################################

      LOGICAL FUNCTION TAGS_MATCH ( TAG_A, TAG_B )

! Continuation markers TAG_A (field 10 of a line) and TAG_B (field 1 of the next) match: equal after their first character,
! which is blank, '+' (small field) or '*' (large field) in each, so an entry may mix small and large field lines

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: TAG_A, TAG_B

      TAGS_MATCH = (INDEX(' +*', TAG_A(1:1)) > 0) .AND. (INDEX(' +*', TAG_B(1:1)) > 0) .AND. (TAG_A(2:8) == TAG_B(2:8))

      END FUNCTION TAGS_MATCH

! ##################################################################################################################################

      CHARACTER(LEN=8) FUNCTION LINE_TAG ( LINE )

! Field 1 (the continuation marker) of a Bulk Data line: its first 8 columns, or for a free field line (a comma before any
! comment, and not a large field '*' line) field 1 after FFIELD has made it fixed field ("+,50.,,23456" has field 1 "+")

      USE PENTIUM_II_KIND, ONLY       :  LONG

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: LINE
      CHARACTER(LEN=LEN(LINE))        :: T
      INTEGER(LONG)                   :: IERR, K

      LINE_TAG = LINE(1:8)
      IF (LINE(1:1) == '*') RETURN
      K = INDEX(LINE(2:), '$')                             ! A comment (after column 1, so at column K+1) ends the data:
      IF (K == 0) K = LEN(LINE)                            ! LINE(1:K)
      IF (INDEX(LINE(1:K), ',') == 0) RETURN
      T = LINE(1:K)
      CALL FFIELD ( T, IERR )
      LINE_TAG = T(1:8)

      END FUNCTION LINE_TAG

   END MODULE BDF_CARD_CONTINUATIONS
