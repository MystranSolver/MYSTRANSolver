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

   MODULE INPUT_FILE_MECHANICS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CSHIFT, FFIELD, FFIELD2, READ_INCLUDE_FILNAM, REPLACE_TABS_W_BLANKS, RW_INCLUDE_FILES, IS_THIS_A_RESTART

   CONTAINS

      SUBROUTINE CSHIFT ( CARD_IN, CHAR, CARD_SHIFTED, CHAR_COL, IERR )

! Shifts card string data on CARD_IN so that the data after character CHAR is shifted to start in col 1 (with blanks
! between CHAR and data on CARD_IN deleted). An error is indicated if CHAR is not found. The special case of CHAR = ' '
! input to this subr indicates we want to shift the card to begin in column 1

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))         :: SUBR_NAME = 'CSHIFT'
      CHARACTER(LEN=*) , INTENT(IN)            :: CARD_IN           ! Input Case Control card
      CHARACTER(LEN=LEN(CARD_IN)) , INTENT(OUT):: CARD_SHIFTED      ! C.C. card shifted to begin in 1st nonblank col after CHAR_COL
      CHARACTER(1*BYTE), INTENT(IN)            :: CHAR              ! Character to find in CARD

      INTEGER(LONG), INTENT(OUT)               :: IERR              ! Error indicator. If CHAR not found, IERR set to 1
      INTEGER(LONG), INTENT(OUT)               :: CHAR_COL          ! Column number on CARD where character CHAR is found
      INTEGER(LONG)                            :: CARD_IN_LEN       ! Length of CARD
      INTEGER(LONG)                            :: I                 ! DO loop index
      INTEGER(LONG)                            :: ISTART            ! The col on CARD where nonblank data begins after CHAR_COL


      INTRINSIC INDEX



! **********************************************************************************************************************************
      CARD_IN_LEN = LEN(CARD_IN)

! Initialize CARD_SHIFTED

      DO I=1,CARD_IN_LEN
         CARD_SHIFTED(I:I) = ' '
      ENDDO

      IERR = 0
      IF (CHAR == ' ') THEN                                ! Special case: shift card to begin in 1st nonblank col after col 1
         CHAR_COL = 0
      ELSE                                                 ! Regular case: shift card to begin in 1st nonblank col after CHAR_COL+1
         CHAR_COL = INDEX(CARD_IN(1:),CHAR)
         IF (CHAR_COL == 0) THEN
            IERR = 1
            RETURN
         ENDIF
      ENDIF

      ISTART = 1
      DO I=CHAR_COL+1,CARD_IN_LEN                          ! Skip over blanks to calc ISTART, the col in CARD_IN to shift to
         IF ((CARD_IN(I:I) == ' ') .OR. (CARD_IN(I:I) == ACHAR(9))) THEN
            CYCLE
         ELSE
            ISTART = I
            EXIT
         ENDIF
      ENDDO
      CARD_SHIFTED(1:) = CARD_IN(ISTART:CARD_IN_LEN)



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CSHIFT


      SUBROUTINE FFIELD ( CARD, IERR )

! Routine to handle only small field input CARD. The 8 col fields of CARD will be expanded to 16 col fields and the returned CARD
! will have 10 fields of 16 cols each:

!  1) Convert fields 2-9 of a free field Bulk Data card to left justified fixed field card. It is assumed that if input CARD has a
!     comma then it is a free-field card

!  2) Left justify fields 2 - 9 of cards that are fixed field

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, BD_ENTRY_LEN, FATAL_ERR, IMB_BLANK, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE BDF_CARD_CONTINUATIONS, ONLY:  MKCARD, MKJCARD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'FFIELD'
      CHARACTER(LEN=*),  INTENT(INOUT):: CARD              !
      CHARACTER( 1*BYTE)              :: FOUND_DATA        !
      CHARACTER( 3*BYTE)              :: FREEFLD           ! = 'Y' if CARD is free field form
      CHARACTER(LEN=JCARD_LEN)        :: LJCARD(10)        ! 10 fields of LCARD
      CHARACTER(LEN=BD_ENTRY_LEN)     :: LCARD             !
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! Fields of CARD
      CHARACTER(LEN=JCARD_LEN)        :: TJCARD(10)        ! Fields of TCARD
      CHARACTER(LEN=LEN(CARD))        :: TCARD             ! Temporary CARD

      INTEGER(LONG)                   :: CARD_LEN          ! Length of CARD
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG), INTENT(OUT)      :: IERR              ! = 1 if a field  is longer than 8 chars on a free field card
      INTEGER(LONG)                   :: IFD               ! Counter for the 10 fields of a Bulk Data CARD
      INTEGER(LONG)                   :: JCT               ! Column counter in free-field CARD
      INTEGER(LONG)                   :: K1S,K2S,K1L       ! Indices




! **********************************************************************************************************************************
      CARD_LEN = LEN(CARD)

! Initialize

      DO I = 1,CARD_LEN
         TCARD(I:I) = ' '
      ENDDO

      DO I=1,10
         TJCARD(I)(1:) = ' '
      ENDDO

      IERR = 0

! Look for ',' on input string 'CARD'; if found, we assume card is free-field

      IF ((INDEX(CARD,',') > 0) .OR. (INDEX(CARD,ACHAR(9)) > 0)) THEN
         FREEFLD = 'YES'
      ELSE
         FREEFLD = 'NO '
      ENDIF

! Process CARD

      IF (FREEFLD == 'NO ') THEN                           ! Expand small field card to 16 cols/field and left justify

         LCARD(1:) = ' '
         DO I=1,10
            K1S =  8*(I-1) + 1   ;   K2S = K1S + 7
            K1L = 16*(I-1) + 1
            LCARD(K1L:K1L+7) = CARD(K1S:K2S)
         ENDDO

         CALL MKJCARD ( SUBR_NAME, LCARD, LJCARD )         ! LJCARD are fields from LCARD (16 col fields)

         DO I = 2,9                                        ! Left justify fields
            IF (LJCARD(I)(1:) == ' ' .OR. LJCARD(I)(1:1) /= ' ') THEN
               TJCARD(I) = LJCARD(I)
            ELSE
               DO J=1,JCARD_LEN
                  IF (LJCARD(I)(J:J) /= ' ') THEN
                     TJCARD(I)(1:) = LJCARD(I)(J:)
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         DO I=2,9
            LJCARD(I) = TJCARD(I)
         ENDDO

         CALL MKCARD ( LJCARD, CARD )

      ELSE                                                 ! Convert free-field 'CARD' to fixed field 'TCARD'

         DO I=1,10
            TJCARD(I)(1:) = ' '
         ENDDO

         I    = 0
         IFD  = 1
         JCT  = 0
loop1:   DO

            I = I + 1
            IF (I > BD_ENTRY_LEN) EXIT loop1

            IF (CARD(I:I)  == ' ') CYCLE loop1

            IF ((CARD(I:I) /= ',') .AND. (CARD(I:I) /= ACHAR(9))) THEN

               JCT = JCT+1

               IF (JCT > 16) THEN                          ! NOTE: free field can only have <= 16 cols
                  WRITE(ERR,1002)
                  WRITE(F06,1002)
                  WRITE(ERR,129) CARD
                  WRITE(F06,129) CARD
                  IERR = 1
                  FATAL_ERR = FATAL_ERR + 1
                  EXIT loop1
               ELSE
                  TJCARD(IFD)(JCT:JCT) = CARD(I:I)
               ENDIF

            ELSE

               IFD = IFD+1
               JCT = 0
               IF (IFD > 10) EXIT loop1

            ENDIF

         ENDDO loop1

         CALL MKCARD ( TJCARD, CARD )

      ENDIF

! Check fields for any imbedded blanks and set error if any are found

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      DO I=2,9
         FOUND_DATA   = 'N'
         IMB_BLANK(I) = 'N'
         DO J=JCARD_LEN,1,-1
            IF (JCARD(I)(J:J) /= ' ') THEN
               FOUND_DATA = 'Y'
            ELSE
               IF(FOUND_DATA == 'Y') THEN
                  IMB_BLANK(I) = 'Y'
               ELSE
                  CYCLE
               ENDIF
            ENDIF
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************
 1002 FORMAT(' *ERROR  1002: TOO LONG AN ENTRY (MORE THAN 16 CHARS) ON THE FOLLOWING ENTRY (MAYBE A COMMA WAS FOUND WHERE ONE',    &
                           ' SHOULD NOT BE):')

  129 FORMAT(A)

! **********************************************************************************************************************************

      END SUBROUTINE FFIELD


      SUBROUTINE FFIELD2 ( CARD1, CARD2, CARD, IERR )

! Routine to process large field format BD entries (must be fixed field - free field not allowed for large field format:

! 1) Input 2 physical 80 col cards (read in LOADB) that form one logical entry

!    a) Card 1 has
!         i) BD entry name (CARD1_FLD1))
!        ii) 4 fields of data (CARD1_FLD2 - CARD1_FLD5))
!       iii) cont entry (CARD1_FLD6)

!    b) Card 2 has
!         i) Cont entry (CARD2_FLD1)
!        ii) 4 fields of data (CARD2_FLD2 - CARD2_FLD5)
!       iii) cont entry (CARD2_FLD6)

! 2) Make sure CARD1_FLD6 is same as CARD2_FLD1 (cont entries must match)

! 3) Expand fields in each entry from 8 to 16 cols

! 4) Create new entry with 10 fields of 16 cols each that has:

!    a) field  1: CARD1_FLD1 (BD card name)
!    b) field  2: CARD1_FLD2
!    c) field  3: CARD1_FLD3
!    d) field  4: CARD1_FLD4
!    e) field  5: CARD1_FLD5
!    f) field  6: CARD2_FLD2
!    g) field  7: CARD2_FLD3
!    h) field  8: CARD2_FLD4
!    i) field  9: CARD2_FLD5
!    j) field 10: CARD2_FLD6 (cont from field 6 of 2nd half of entry)

! 5) Left justify fields

! N O TE : each of the 2 physical entries making up 1 logical large field entry has 80 cols in 6 fields (1st and fields are 8 cols
!          and 2nd - 6th fields are large field format with 16 cols each)


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, BD_ENTRY_LEN, ECHO, FATAL_ERR, IMB_BLANK, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE BDF_CARD_CONTINUATIONS, ONLY:  MKCARD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'FFIELD2'
      CHARACTER(LEN=*),  INTENT(IN)   :: CARD1             ! 1st physical entry of the large field entry
      CHARACTER(LEN=*),  INTENT(IN)   :: CARD2             ! 2nd physical entry of the large field entry
      CHARACTER(LEN=*),  INTENT(OUT)  :: CARD              ! Card with 10 fields of 16 cols each with the data from CARD1, CARD2
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! 10 fields of CARD . All 10 fields are 16 cols wide
      CHARACTER(LEN=JCARD_LEN)        :: JCARD1(6)         !  6 fields of CARD1. Fields 1,6 are 8 cols. Fields 2,3,4,5 are 16 cols
      CHARACTER(LEN=JCARD_LEN)        :: JCARD2(6)         !  6 fields of CARD2. Fields 1,6 are 8 cols. Fields 2,3,4,5 are 16 cols
      CHARACTER(LEN=JCARD_LEN)        :: TJCARD(10)        ! Temporary JCARD's

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG), INTENT(OUT)      :: IERR              ! = 1 if a field  is longer than 8 chars on a free field card




! **********************************************************************************************************************************
! Initialize

      IERR = 0

      CARD(1:) = ' '

      DO I=1,6
         JCARD1(I)(1:) = ' '
         JCARD2(I)(1:) = ' '
      ENDDO

! Make sure the entries are not free field

      IF ((INDEX(CARD1,',') > 0) .OR. (INDEX(CARD1,ACHAR(9)) > 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         IF (ECHO == 'NONE  ') THEN
            IF (SUPWARN == 'N') THEN
               WRITE(F06,129) CARD1
            ENDIF
         ENDIF
         WRITE(ERR,1014)
         WRITE(F06,1014)
         RETURN
      ENDIF

! Put fields of CARD1 and CARD2 into 16 col fields of JCARD1(i=1,6), JCARD2(i=1,6)

      JCARD1(1)(1: 8) = CARD1( 1: 8)   ;   JCARD1(1)(9:16) = ' '
      JCARD1(2)(1:16) = CARD1( 9:24)
      JCARD1(3)(1:16) = CARD1(25:40)
      JCARD1(4)(1:16) = CARD1(41:56)
      JCARD1(5)(1:16) = CARD1(57:72)
      JCARD1(6)(1: 8) = CARD1(73:80)

      JCARD2(1)(1: 8) = CARD2( 1: 8)
      JCARD2(2)(1:16) = CARD2( 9:24)
      JCARD2(3)(1:16) = CARD2(25:40)
      JCARD2(4)(1:16) = CARD2(41:56)
      JCARD2(5)(1:16) = CARD2(57:72)
      JCARD2(6)(1: 8) = CARD2(73:80)   ;   JCARD2(6)(9:16) = ' '

! Make sure that CARD2 is the 2nd half of CARD1 (continuation entry from field 6 of CARD1 must match field 1 of CARD2)

      IERR = 0
      IF      (JCARD2(1)(1:) == JCARD1(6)(1:)) THEN
         CONTINUE
      ELSE IF ((JCARD2(1)(1:1) == '*') .AND. (JCARD1(6)(1:1) == ' ') .AND. (JCARD2(1)(2:8) == JCARD1(6)(2:8))) THEN
         CONTINUE
      ELSE IF ((JCARD2(1)(1:1) == ' ') .AND. (JCARD1(6)(1:1) == '*') .AND. (JCARD2(1)(2:8) == JCARD1(6)(2:8))) THEN
         CONTINUE
      ELSE
         IERR = IERR + 1
      ENDIF

      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         IF (ECHO == 'NONE  ') THEN
            IF (SUPWARN == 'N') THEN
               WRITE(F06,129) CARD1
               WRITE(F06,129) CARD2
            ENDIF
         ENDIF
         WRITE(ERR,1020) JCARD1(6), JCARD2(1)
         WRITE(F06,1020) JCARD1(6), JCARD2(1)
         RETURN
      ENDIF

! Assemble output CARD from fields of CARD1,2

      JCARD( 1) = JCARD1(1)
      JCARD( 2) = JCARD1(2)
      JCARD( 3) = JCARD1(3)
      JCARD( 4) = JCARD1(4)
      JCARD( 5) = JCARD1(5)

      JCARD( 6) = JCARD2(2)
      JCARD( 7) = JCARD2(3)
      JCARD( 8) = JCARD2(4)
      JCARD( 9) = JCARD2(5)
      JCARD(10) = JCARD2(6)

! Left justify fields of CARD

      DO I = 2,9                                        ! Left justify fields
         IF (JCARD(I)(1:) == ' ' .OR. JCARD(I)(1:1) /= ' ') THEN
            TJCARD(I) = JCARD(I)
         ELSE
            DO J=1,JCARD_LEN
               IF (JCARD(I)(J:J) /= ' ') THEN
                  TJCARD(I)(1:) = JCARD(I)(J:)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      DO I=2,9
         JCARD(I) = TJCARD(I)
      ENDDO

! Put left justified fields into CARD

      CALL MKCARD ( JCARD, CARD )



      RETURN

! **********************************************************************************************************************************
  129 FORMAT(A)

 1014 FORMAT(' *ERROR  1014: FREE FIELD FORMAT NOT ALLOWED FOR LARGE FIELD FORMAT BULK DATA ENTRIES')

 1020 FORMAT(' *ERROR  1020: 2ND PHYSICAL ENTRY OF A LARGE FIELD ENTRY IS NOT A CONTINUATION TO THE 1ST ENTRY.'                    &
                    ,/,14X,' THE LAST FIELD OF THE IST ENTRY = "',A,'" AND THE FIRST FIELD OF THE 2ND ENTRY = "',A,'"')

! ##################################################################################################################################

      END SUBROUTINE FFIELD2


      SUBROUTINE READ_INCLUDE_FILNAM ( CARD, IERR )

! If there is an INCLUDE entry this subr reads the file name from that entry. The entry will be of the form:
!     INCLUDE 'filename' with or without the ' marks

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, FILE_NAM_MAXLEN, INC, INCFIL
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, EC_ENTRY_LEN, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OPNERR

      IMPLICIT NONE

      LOGICAL                         :: LEXIST            ! 'T' if INCFIL exists

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'READ_INCLUDE_FILNAM'
      CHARACTER(LEN=EC_ENTRY_LEN), INTENT(IN)  :: CARD     ! An entry from an input data (DAT) file
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD1             ! CARD shifted to begin in col 1
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD2             ! CARD1 with the "INCLUDE" removed
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD3             ! CARD1 with the "INCLUDE" removed
      CHARACTER( 1*BYTE)              :: DONE              ! Indicator of having found start and end of file name

      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Local error count

      INTEGER(LONG)                   :: CHAR_COL          ! Column number on CARD where character CHAR is found
      INTEGER(LONG)                   :: DELTA_END_COL     ! Delta from START_ COL to END_COL
      INTEGER(LONG)                   :: END_COL           ! Col from CARD1 where the 2nd ' exists, if it does exist
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a Case Control card from unit IN1
      INTEGER(LONG)                   :: START_COL         ! Col from CARD1 where the 1st ' exists, if it does exist
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to




! **********************************************************************************************************************************
! Initialize

      INCFIL(1:) = ' '
      CALL CSHIFT ( CARD , ' ', CARD1, CHAR_COL, IERR )    ! This will remove leading blanks
      CARD2(1:)  = ' '
      CARD2(1:)  = CARD1(8:)                               ! CARD2 has only the file name plus, perhaps, leading/trailing blanks
      CALL CSHIFT ( CARD2, ' ', CARD3, CHAR_COL, IERR)     ! CARD3 is CARD2 without leading blanks
      DONE       = 'N'
      IERR       = 0

! Default units for writing errors the ERR, F06 files

      OUNT(1) = ERR
      OUNT(2) = F06

! Case 1: there are ' marks at the beginning and end of the file name

      START_COL = INDEX(CARD3,'''',.FALSE.) + 1
      IF (START_COL > 1) THEN                              ! 0 doesn't work, for some reason
         DELTA_END_COL = INDEX(CARD3(START_COL+1:),'''')
         IF (DELTA_END_COL > 0) THEN
            END_COL = START_COL + DELTA_END_COL - 1
            INCFIL(1:) = CARD3(START_COL:END_COL)
            DONE = 'Y'
         ELSE                                              ! There was a ' mark at the beginning but not at the end
            WRITE(ERR,1044) 'a', CARD
            WRITE(F06,1044) 'a', CARD
            IERR = 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

! Case 2: there are no quote marks

      IF (DONE == 'N') THEN

         DO I=EC_ENTRY_LEN,1,-1                            ! Find START_COL when there is no ' mark
            IF ((CARD3(I:I) == ' ') .OR. (CARD3(I:I) == ACHAR(9))) THEN
               CYCLE
            ELSE
               END_COL = I
               EXIT
            ENDIF
         ENDDO

         DO I=1,EC_ENTRY_LEN,1
            IF ((CARD3(I:I) == ' ') .OR. (CARD3(I:I) == ACHAR(9))) THEN
               CYCLE
            ELSE
               START_COL = I
               EXIT
            ENDIF
         ENDDO

         IF (END_COL > START_COL) THEN
            INCFIL(1:) = CARD3(START_COL:END_COL)
            DONE = 'Y'
         ELSE
            WRITE(ERR,1044) 'b', CARD
            WRITE(F06,1044) 'b', CARD
            IERR = 2
            FATAL_ERR = FATAL_ERR + 1
         ENDIF

      ENDIF

      IF (DONE == 'Y') THEN

         INCFIL(1:) = CARD3(START_COL:END_COL)

         INQUIRE(FILE=INCFIL,EXIST=LEXIST)
         IF (LEXIST) THEN
            OPEN(INC,FILE=INCFIL,STATUS='OLD',ACTION='READ',IOSTAT=IOCHK)
            IF (IOCHK /= 0) CALL OPNERR ( IOCHK, INCFIL, OUNT )
         ELSE
            WRITE(ERR,1042) INCFIL
            WRITE(F06,1042) INCFIL
            FATAL_ERR = FATAL_ERR + 1
            IERR = 3
         ENDIF

      ELSE

         WRITE(ERR,1044) 'c', CARD
         WRITE(F06,1044) 'c', CARD
         IERR = 4
         FATAL_ERR = FATAL_ERR + 1

      ENDIF

      IF (DEBUG(115) > 0) THEN
         CALL DEB_READ_INCL_FILNAM
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1042 FORMAT(' *ERROR  1042: THE FOLLOWING "INCLUDE" FILE  DOES NOT EXIST OR THE INCLUDE ENTRY WAS IN ERROR: "',A,'"')

 1044 FORMAT(' *ERROR  1044',A,':INCORRECT FORMAT FOR "INCLUDE" ENTRY: ',A)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEB_READ_INCL_FILNAM

      IMPLICIT NONE

! **********************************************************************************************************************************
      WRITE(F06,'(A,A,A)') ' INCLUDE statement   = "', CARD, '"'
      WRITE(F06,'(A,A,A)') ' INCLUDE 2nd word    = "', CARD3, '"'
      WRITE(F06,'(A,I4 )') ' START_COL           = ' , START_COL
      WRITE(F06,'(A,I4 )') ' END_COL             = ' , END_COL
      WRITE(F06,'(A,A,A)') ' INCLUDE filename is = "', INCFIL(1:END_COL-START_COL+1), '"'

! **********************************************************************************************************************************

      END SUBROUTINE DEB_READ_INCL_FILNAM

      END SUBROUTINE READ_INCLUDE_FILNAM



      SUBROUTINE REPLACE_TABS_W_BLANKS ( CARD )

! Searches input CARD for tab characters and replaces them with 1 white space character. Used primarily for Exec Control and
! Case Control entries (but not Bulk Data entries - which are handled differently).

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REPLACE_TABS_W_BLANKS'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! Input entry character line
      CHARACTER(LEN=LEN(CARD))        :: CARD0             ! Temporary CARD


      INTEGER(LONG)                   :: I                 ! DO loop index



! **********************************************************************************************************************************
! Strip all tab chars from input CARD

         CARD0(1:) = ' '
         DO I=1,LEN(CARD)
            IF (CARD(I:I) == ACHAR(9)) THEN
               CARD0(I:I) = ' '
            ELSE
               CARD0(I:I) = CARD(I:I)
            ENDIF
         ENDDO
         CARD(1:) = CARD0(1:)



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE REPLACE_TABS_W_BLANKS



      SUBROUTINE RW_INCLUDE_FILES ( UNIT_IN, UNIT_OUT )

! Reads card images from INCLUDE files and writes them out to the file that will have the complete input data (DAT file + INCLUDE
! files entries)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, FILE_NAM_MAXLEN, INCFIL
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, EC_ENTRY_LEN, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RW_INCLUDE_FILES'
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD              ! Entry from INCL_FILNAM


      INTEGER(LONG), INTENT(IN)       :: UNIT_IN           ! Unit number to read  INCLUDE entries from
      INTEGER(LONG), INTENT(IN)       :: UNIT_OUT          ! Unit number to write INCLUDE entries to
      INTEGER(LONG)                   :: ICNT        = 0   ! Counter
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading an entry from INCL_FILNAM



! **********************************************************************************************************************************
      WRITE(UNIT_OUT,201,IOSTAT=IOCHK) INCFIL
      IF (IOCHK > 0) THEN
         WRITE(ERR,1029) INCFIL
         WRITE(F06,1029) INCFIL
         WRITE(F06,'(A)') INCFIL
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

      ICNT = 0

main: DO

         READ(UNIT_IN,101,IOSTAT=IOCHK) CARD
         ICNT = ICNT + 1

         IF (IOCHK < 0) THEN                               ! When end of file is encountered, no more entries to read
            IF (ICNT == 1) THEN
               WRITE(ERR,1041) INCFIL
               WRITE(F06,1041) INCFIL
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
            EXIT main
         ENDIF

         IF (IOCHK > 0) THEN

            WRITE(ERR,1010) ICNT, INCFIL
            WRITE(F06,1010) ICNT, INCFIL
            FATAL_ERR = FATAL_ERR + 1
            CYCLE

         ELSE

            WRITE(UNIT_OUT,101,IOSTAT=IOCHK) CARD
            IF (IOCHK > 0) THEN
               WRITE(ERR,1029) ICNT, INCFIL
               WRITE(F06,1029) ICNT, INCFIL
               FATAL_ERR = FATAL_ERR + 1
               CYCLE
            ENDIF

         ENDIF

      ENDDO main

      WRITE(UNIT_OUT,202,IOSTAT=IOCHK) INCFIL

      IF (IOCHK > 0) THEN
         WRITE(ERR,1029) ICNT, INCFIL
         WRITE(F06,1029) ICNT, INCFIL
         FATAL_ERR = FATAL_ERR + 1
      ELSE
         WRITE(F06,301) ICNT-1, INCFIL
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

  201 FORMAT('$ -----------------------------------------------------------------------------',/,'$ Beg INCLUDE file ',A)

  202 FORMAT('$ End INCLUDE file ',A,/,'$ -----------------------------------------------------------------------------')

  301 FORMAT(' *INFORMATION: THERE WERE ',I8,' LINE(S) INCLUDED FROM FILE: ',A)

 1010 FORMAT(' *ERROR  1010: ERROR READING ENTRY NUMBER ',I8,' FROM INCLUDE FILE: ',A)

 1029 FORMAT(' *ERROR  1029: ERROR WRITING ENTRY NUMBER ',I8,' FROM INCLUDE FILE: ',A)

 1041 FORMAT(' *ERROR  1041: NO DATA IN INCLUDE FILE: ',A)

! **********************************************************************************************************************************

      END SUBROUTINE RW_INCLUDE_FILES



      SUBROUTINE IS_THIS_A_RESTART

! IS_THIS_A_RESTART reads in the EXEC CONTROL DECK to find if there is a RESTART entry

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  IN1, SC1
      USE SCONTR, ONLY                :  EC_ENTRY_LEN, FATAL_ERR, RESTART
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD              ! Exec Control deck card
      CHARACTER(LEN=EC_ENTRY_LEN)     :: CARD1             ! CARD shifted to begin in col 1
      CHARACTER(12*BYTE)              :: DECK_NAME   = 'EXEC CONTROL'
      CHARACTER( 4*BYTE), PARAMETER   :: END_CARD  = 'CEND'

      INTEGER(LONG)                   :: CHAR_COL          ! Column number on CARD where character CHAR is found
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator.
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a Case Control card from unit IN1

! **********************************************************************************************************************************
! Process EXECUTIVE CONTROL DECK

      RESTART        = 'N'

      DO

         READ(IN1,101,IOSTAT=IOCHK) CARD

         IF (IOCHK < 0) THEN                               ! Quit if EOF/EOR occurs during read
            WRITE(SC1,1011) END_CARD
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         IF (IOCHK > 0) THEN                               ! Check if error occurs during read.
            WRITE(SC1,1010) DECK_NAME
            WRITE(SC1,'(A)') CARD
            FATAL_ERR = FATAL_ERR + 1
            CYCLE
         ENDIF

         CALL CSHIFT ( CARD, ' ', CARD1, CHAR_COL, IERR )

         IF (CARD1(1:7) == 'RESTART'   ) THEN              ! No errors, so look for RESTART
            RESTART = 'Y'

         ELSE IF (CARD1(1:4) == 'CEND' ) THEN              ! Check for CEND card
            EXIT

         ENDIF

      ENDDO

! **********************************************************************************************************************************
  101 FORMAT(A)

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ',A,' ENTRY. ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

! **********************************************************************************************************************************

      END SUBROUTINE IS_THIS_A_RESTART

   END MODULE INPUT_FILE_MECHANICS

      SUBROUTINE CSHIFT ( CARD_IN, CHAR, CARD_SHIFTED, CHAR_COL, IERR )
      USE INPUT_FILE_MECHANICS, ONLY: CSHIFT_MOD => CSHIFT
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: CARD_IN
      CHARACTER(1), INTENT(IN) :: CHAR
      CHARACTER(LEN=*), INTENT(OUT) :: CARD_SHIFTED
      INTEGER, INTENT(OUT) :: CHAR_COL, IERR
      CALL CSHIFT_MOD ( CARD_IN, CHAR, CARD_SHIFTED, CHAR_COL, IERR )
      END SUBROUTINE CSHIFT

      SUBROUTINE FFIELD ( CARD, IERR )
      USE INPUT_FILE_MECHANICS, ONLY: FFIELD_MOD => FFIELD
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD
      INTEGER, INTENT(OUT) :: IERR
      CALL FFIELD_MOD ( CARD, IERR )
      END SUBROUTINE FFIELD

      SUBROUTINE FFIELD2 ( CARD1, CARD2, CARD, IERR )
      USE INPUT_FILE_MECHANICS, ONLY: FFIELD2_MOD => FFIELD2
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: CARD1, CARD2
      CHARACTER(LEN=*), INTENT(OUT) :: CARD
      INTEGER, INTENT(OUT) :: IERR
      CALL FFIELD2_MOD ( CARD1, CARD2, CARD, IERR )
      END SUBROUTINE FFIELD2
