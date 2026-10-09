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

   MODULE CASE_CONTROL_SETS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CC_SET0, CC_SET, CC_SUBC

   CONTAINS

      SUBROUTINE CC_SET0 ( CARD )

! Processes Case Control SET cards to determine LSETLN, the length of all SET characters that go into array
! ALL_STES_ARRAY. No error messages are written if, in trying to determine LSETLN, errors occur. This subr returns
! and, when CC_SET runs (when called by LOADC), it will discover the same errors and report them before attempting
! to write characters to ALL_SETS_ARRAY

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, IN1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_ENTRY_LEN, LSETLN
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_SET0'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=LEN(CARD))        :: CARD1             ! SET card read from C.C deck
      CHARACTER( 8*BYTE)              :: TOKEN             ! An 8 char string from the SET card where the set ID should be located
      CHARACTER( 8*BYTE)              :: TOKTYP            ! The type of the char string TOKEN

      INTEGER(LONG)                   :: ECOL      = 0     ! Col, on CARD, where "=" sign is located
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator if there is a cont card (yes if last entry is ",")
      INTEGER(LONG)                   :: MORE      = 0     ! Count of additional bytes to add to LSETLN as this card is read
      INTEGER(LONG)                   :: IOCHK     = 0     ! IOSTAT error number when reading a Case Control card from unit IN1
      INTEGER(LONG)                   :: K         = 0     ! Counter
      INTEGER(LONG)                   :: SETERR    = 0     ! Error indicator as set ID is read


      INTRINSIC INDEX



! **********************************************************************************************************************************
      CARD1 = CARD

! Process SET cards only to count LSETLN. Later processing in CC_SET will put data into ALL_SETS_ARRAY.

! Check for SET cards. These can have continuation if last entry is a comma. Need to make sure that equal sign is
! included for use in subr SETPRO which processes the character string ALL_SETS_ARRAY. A logical SET card consists
! of all physical cards in a set.

! Make sure equal sign is in. We use it later in subcase processor

      ECOL = INDEX(CARD1(1:),'=')
      IF (ECOL == 0) THEN
         RETURN
      ENDIF

! Now find SET ID and check for proper type

      TOKEN = '        '
      SETERR = 0
      K = 0
      DO I=5,ECOL-1
         IF (CARD1(I:I) == ' ') CYCLE
         K = K + 1
         IF (K > 8) THEN
            SETERR = 1
            EXIT
         ENDIF
         TOKEN(K:K) = CARD1(I:I)
      ENDDO
      IF (SETERR == 0) THEN
         CALL TOKCHK ( TOKEN, TOKTYP )
         IF (TOKTYP /= 'INTEGER ') THEN
            RETURN
         ENDIF
      ENDIF

! Find out if there are continuation cards (last non-blank entry on  this card will be a comma if there is more)

      DO

! Get rid of trailing blanks

         MORE = CC_ENTRY_LEN
         DO I=CC_ENTRY_LEN,1,-1
            IF (CARD1(I:I) == ' ') THEN
               MORE = MORE - 1
            ELSE
               EXIT
            ENDIF
         ENDDO

         LSETLN = LSETLN + MORE

! Find out if last entry is a ','

         ICONT = 0
         DO I=CC_ENTRY_LEN,1,-1
            IF (CARD1(I:I) == ' ') THEN
               CYCLE
            ELSE IF (CARD1(I:I) == ',') THEN
               ICONT = 1
               EXIT
            ELSE
               ICONT = 0
               EXIT
            ENDIF
         ENDDO

         IF (ICONT == 1) THEN
            CALL READ_BDF_LINE(IN1, IOCHK, CARD1)
            IF (IOCHK > 0) THEN
               RETURN
            ENDIF
            CYCLE
         ELSE
            EXIT
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

! **********************************************************************************************************************************

      END SUBROUTINE CC_SET0


      SUBROUTINE CC_SET ( CARD )

! Processes Case Control SET cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, IN1
      USE IOUNT1, ONLY                :  WRT_ERR
      USE SCONTR, ONLY                :  CC_ENTRY_LEN, FATAL_ERR, LSETS, LSETLN, MAX_TOKEN_LEN, NSETS, SETLEN, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ALL_SETS_ARRAY, SETS_IDS

      USE BDF_SET_SYNTAX, ONLY        :  STOKEN, TOKCHK
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_SET'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=LEN(CARD))        :: CARD1             ! SET card read from C.C deck
      CHARACTER( 3*BYTE), PARAMETER   :: CARD_NAME = 'SET' ! Name of the C.C. card
      CHARACTER( 1*BYTE)              :: COMMENT_FOUND     ! 'Y' if $ was found on an CARD
      CHARACTER(LSETLN*BYTE)          :: ERRTOK            ! Character string that holds part of an error message from subr STOKEN
      CHARACTER( 3*BYTE)              :: EXCEPT            ! An input/output for subr STOKEN
      CHARACTER( 3*BYTE)              :: THRU              ! An input/output for subr STOKEN
      CHARACTER( 8*BYTE)              :: TOKEN(3)          ! Char string output from subr STOKEN, called herein
      CHARACTER( 8*BYTE)              :: TOKTYP(3)         ! Type of the char TOKEN's output from subr STOKEN, called herein
      CHARACTER(LSETLN*BYTE)          :: TOKSTR_PART       ! TOKSTR_WHOLE without the "SET setid =" and following blanks
      CHARACTER(LSETLN*BYTE)          :: TOKSTR_WHOLE      ! All characters that are in the current SET, including "SET setid ="

      INTEGER(LONG)                   :: DATA_BEG          ! Posn in TOKSTR_WHOLE where data begins (after "SET setid =" and blanks)
      INTEGER(LONG)                   :: DATA_END          ! Posn in TOKSTR_WHOLE where data ends
      INTEGER(LONG)                   :: ECOL      = 0     ! Col, on CARD, where "=" sign is located
      INTEGER(LONG)                   :: ACT_ID            ! An actual grid or elem ID that appears in the SET definition
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: II                ! Temporary value for IEND
      INTEGER(LONG)                   :: I1                ! Lower  end of range in a "I1 THRU I2" portion of the SET definition
      INTEGER(LONG)                   :: I2                ! Higher end of range in a "I1 THRU I2" portion of the SET definition
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator if there is a cont card (yes if last entry is ",")
      INTEGER(LONG)                   :: IEND              ! The col no. on a physical card of the SET where data ends
      INTEGER(LONG)                   :: IERROR            ! An output from subr STOKEN, called herein
      INTEGER(LONG)                   :: IOCHK1    = 0     ! IOSTAT error number when reading a Case Control card from unit IN1
      INTEGER(LONG)                   :: IOCHK2    = 0     ! IOSTAT error number when reading a Case Control card from unit IN1
      INTEGER(LONG)                   :: ISTART    = 0     ! An input to subr STOKEN, called herein
      INTEGER(LONG)                   :: K         = 0     ! Counter
      INTEGER(LONG)                   :: NTOKEN            ! An output from subr STOKEN, called herein
      INTEGER(LONG)                   :: SETERR    = 0     ! Error indicator as set ID is read
      INTEGER(LONG)                   :: SETID     = 0     ! Set ID on this Case Control card
      INTEGER(LONG)                   :: TOKLEN    = 0     ! DATA_END - DATA_BEG + 1 (an input to subr STOKEN, called herein)


      INTRINSIC INDEX



! **********************************************************************************************************************************
! Process SET cards

      CARD1 = CARD

! Check for SET cards. These can have continuation if last entry is a comma. Need to make sure that equal sign is
! included for use in subroutine SETPRO which processes the character string ALL_SETS_ARRAY. A logical SET card
! consists of all physical cards in a set. This is written to ALL_SETS_ARRAY beginning with 'SET' and ending with
! the last entry on the last physical card for the set. Subsequent SET cards are written immediately after the
! previous entry in ALL_SETS_ARRAY.

! Make sure equal sign is in. We use it later in subcase processor

      ECOL = INDEX(CARD1(1:),'=')
      IF (ECOL == 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1267)
         WRITE(F06,1267)
         RETURN
      ENDIF

! Now find SET ID and check for proper type ('INTEGER ')

      SETID    = 0
      TOKEN(1) = '        '
      SETERR   = 0
      K = 0
      DO I=5,ECOL-1
         IF (CARD1(I:I) == ' ') CYCLE
         K = K + 1
         IF (K > MAX_TOKEN_LEN) THEN
            SETERR = 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1272) MAX_TOKEN_LEN
            WRITE(F06,1272) MAX_TOKEN_LEN
            EXIT
         ENDIF
         TOKEN(1)(K:K) = CARD1(I:I)
      ENDDO
      IF (SETERR == 0) THEN
         CALL TOKCHK ( TOKEN(1), TOKTYP(1) )
         IF (TOKTYP(1) /= 'INTEGER ') THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1273)
            WRITE(F06,1273)
            RETURN
         ELSE
            READ(TOKEN(1),'(I8)') SETID
         ENDIF
      ENDIF

      ! Update NSETS and check that other sets with same ID not already defined
      NSETS = NSETS + 1
      IF (NSETS > LSETS) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1163) SUBR_NAME,CARD_NAME,LSETS
         WRITE(F06,1163) SUBR_NAME,CARD_NAME,LSETS
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
      ENDIF
      DO I=1,NSETS-1
         IF (SETID == SETS_IDS(I)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1145) SETID
            WRITE(F06,1145) SETID
            EXIT
         ENDIF
      ENDDO
      SETS_IDS(NSETS) = SETID

! Write SET data to ALL_SETS_ARRAY and check for continuation card(s) in the following DO loop

      DATA_END = 0
main: DO

! Get rid of trailing blanks before storing string in ALL_SETS_ARRAY

         IEND = CC_ENTRY_LEN
i_do_1:  DO I=CC_ENTRY_LEN,1,-1
            IF (CARD1(I:I) == ' ') THEN
               IEND = IEND - 1
            ELSE
               EXIT i_do_1
            ENDIF
         ENDDO i_do_1

! Now check for comment (begins with $). If found, set IEND to col before

         COMMENT_FOUND = 'N'
         II = IEND
i_do_2:  DO I=1,II
            IF (CARD1(I:I) == '$') THEN
               COMMENT_FOUND = 'Y'
               IEND = I - 1
               EXIT i_do_2
            ENDIF
         ENDDO i_do_2

! Check length of character string against max allowable length

         IF ((SETLEN + IEND) > LSETLN) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1268) SUBR_NAME,LSETLN
            WRITE(F06,1268) SUBR_NAME,LSETLN
            CALL OUTA_HERE ( 'Y' )                                 ! Coding error, so quit
         ENDIF

! Put data in ALL_SETS_ARRAY (and also write the data to TOKSTR_WHOLE, which will contain only this set)

i_do_3:  DO I=1,IEND
            ALL_SETS_ARRAY(SETLEN+I)        = CARD1(I:I)
            DATA_END                        = DATA_END + 1
            TOKSTR_WHOLE(DATA_END:DATA_END) = CARD1(I:I)
         ENDDO i_do_3
         SETLEN = SETLEN + IEND

! If a comment was detected, exit main loop

         IF (COMMENT_FOUND == 'Y') THEN
            EXIT main
         ENDIF

! Is there a continuation card? If last entry is a comma, we assume that the next card is a continuation of the set
! definition. If it is not, then an error will be detected when this set's data is checked for syntax below.

         ICONT = 0
i_do_4:  DO I=CC_ENTRY_LEN,1,-1
            IF (CARD1(I:I) == ' ') THEN
               CYCLE
            ELSE IF (CARD1(I:I) == ',') THEN
               ICONT = 1
               EXIT i_do_4
            ELSE
               ICONT = 0
               EXIT i_do_4
            ENDIF
         ENDDO i_do_4

         ! If there is a continuation card, read it
         IF (ICONT == 1) THEN
            CALL READ_BDF_LINE(IN1, IOCHK1, CARD1)
            WRITE(F06,'(A)') CARD1
            IF (IOCHK1 > 0) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1269)
               WRITE(F06,1269)
               RETURN
            ENDIF
            CYCLE main
         ELSE
            EXIT main
         ENDIF

      ENDDO main

      ! Now check that we do not have a situation where a comment was
      ! entered on a card that was not the last card of the set.
      ! If the last character in ALL_SETS_ARRAY is a comma, then we
      ! give fatal error

i_do5:DO I=SETLEN,1,-1
         IF      (ALL_SETS_ARRAY(I) == ' ') THEN
            CYCLE
         ELSE IF (ALL_SETS_ARRAY(I) == ',') THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR, 1252)
            WRITE(ERR,12520) CARD1
            WRITE(ERR,12532)
            WRITE(F06, 1252)
            WRITE(F06,12520) CARD1
            WRITE(F06,12532)
         ELSE
            EXIT i_do5
         ENDIF
      ENDDO i_do5

! TOKSTR_WHOLE has everything in the current set, including "SET ID =". In TOKSTR_PART (below) we want only the set data,
! which begins after the equal sign (and blanks following the equal sign) in the set definition.
! Find where the data begins. The set data ends at DATA_END.

      DATA_BEG = ECOL + 1
      DO I=DATA_BEG,DATA_END
         IF (TOKSTR_WHOLE(I:I) == ' ') THEN
            DATA_BEG = DATA_BEG + 1
         ELSE
            EXIT
         ENDIF
      ENDDO

! Strip, from TOKSTR_WHOLE, everything up to the beginning of the set data (i.e., strip "SET ID = " ) and put into TOKSTR_PART

      TOKLEN                = DATA_END - DATA_BEG + 1
      TOKSTR_PART(1:TOKLEN) = TOKSTR_WHOLE(DATA_BEG:DATA_END)

! **********************************************************************************************************************************
      ! Check syntax of this set data
      ISTART = 1
      THRU   = 'OFF'
      EXCEPT = 'OFF'
      DO

         ! Call STOKEN in a DO loop until we run out of tokens in TOKSTR_PART.
         ! We run out of data in the set when ISTART > TOKLEN

         CALL STOKEN ( SUBR_NAME, TOKSTR_PART, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )

         ! Error: too long a token in the set definition
         IF (IERROR == 1) THEN
            FATAL_ERR = FATAL_ERR + 1

            WRITE(ERR,1252) SETID
            WRITE(F06,1252) SETID

            WRITE(ERR,12520) CARD1
            WRITE(F06,12520) CARD1

            WRITE(ERR,12521) MAX_TOKEN_LEN,ERRTOK(1:TOKLEN)
            WRITE(F06,12521) MAX_TOKEN_LEN,ERRTOK(1:TOKLEN)

            WRITE(F06,*)

            EXIT

! Error: "I1 THRU I2" is missing I2 in the set definition

         ELSE IF (IERROR == 2) THEN

            FATAL_ERR = FATAL_ERR + 1

            WRITE(ERR,1252) SETID
            WRITE(F06,1252) SETID

            WRITE(ERR,12520) CARD1
            WRITE(F06,12520) CARD1

            IF      (NTOKEN == 1) THEN
               WRITE(ERR,12530) NTOKEN,(TOKEN(I),I=1,NTOKEN)
               WRITE(F06,12530) NTOKEN,(TOKEN(I),I=1,NTOKEN)
            ELSE IF (NTOKEN == 2) THEN
               WRITE(ERR,12531) NTOKEN,(TOKEN(I),I=1,NTOKEN)
               WRITE(F06,12531) NTOKEN,(TOKEN(I),I=1,NTOKEN)
            ENDIF

            WRITE(F06,*)

            EXIT

         ! Error: found "EXCEPT" but EXCEPT has already been turned "ON "
         ELSE IF (IERROR == 3) THEN

            FATAL_ERR = FATAL_ERR + 1

            WRITE(ERR,1252) SETID
            WRITE(F06,1252) SETID

            WRITE(ERR,12520) CARD1
            WRITE(F06,12520) CARD1

            WRITE(ERR,12522)
            WRITE(F06,12522)

            WRITE(F06,*)

            EXIT

! Error: found "EXCEPT" but "THRU" had not been turned "ON" yet

         ELSE IF (IERROR == 4) THEN

            FATAL_ERR = FATAL_ERR + 1

            WRITE(ERR,1252) SETID
            WRITE(F06,1252) SETID

            WRITE(ERR,12520) CARD1
            WRITE(F06,12520) CARD1

            WRITE(ERR,12523)
            WRITE(F06,12523)

            WRITE(F06,*)

            EXIT

! Is the TOKEN a number not connected to a 'THRU'?

         ELSE IF (NTOKEN == 1) THEN
            IF (TOKTYP(1) == 'INTEGER ') THEN
               READ(TOKEN(1),'(I8)',IOSTAT=IOCHK1) ACT_ID
               IF (IOCHK1 > 0) THEN

                  FATAL_ERR = FATAL_ERR + 1

                  WRITE(ERR,1252) SETID
                  WRITE(F06,1252) SETID

                  WRITE(ERR,12520) CARD1
                  WRITE(F06,12520) CARD1

                  WRITE(ERR,12524)
                  WRITE(F06,12524)

                  WRITE(F06,*)

                  EXIT

               ELSE
                  IF (EXCEPT == 'ON ') THEN
                     IF ((ACT_ID >= I1) .AND. (ACT_ID <= I2)) THEN
                        CONTINUE
                     ELSE
                        EXCEPT = 'OFF'
                        THRU   = 'OFF'
                     ENDIF
                  ENDIF
               ENDIF

            ELSE IF (TOKTYP(1) == 'EXCEPT  ') THEN
               IF (THRU /= 'ON ') THEN

                  FATAL_ERR = FATAL_ERR + 1

                  WRITE(ERR,1252) SETID
                  WRITE(F06,1252) SETID

                  WRITE(ERR,12520) CARD1
                  WRITE(F06,12520) CARD1

                  WRITE(ERR,12525)
                  WRITE(F06,12525)

                  WRITE(F06,*)

                  EXIT

               ENDIF

            ELSE

               FATAL_ERR = FATAL_ERR + 1

               WRITE(ERR,1252) SETID
               WRITE(F06,1252) SETID

               WRITE(ERR,12520) CARD1
               WRITE(F06,12520) CARD1

               WRITE(ERR,12526) TOKEN(1)
               WRITE(F06,12526) TOKEN(1)

               WRITE(F06,*)

               EXIT

            ENDIF

! TOKEN's connected by 'THRU'

         ELSE IF (NTOKEN == 3) THEN
            IF ((TOKTYP(1) == 'INTEGER ') .AND. (TOKTYP(2) == 'THRU    ') .AND. (TOKTYP(3) == 'INTEGER ')) THEN
               READ(TOKEN(1),'(I8)',IOSTAT=IOCHK1) I1
               IF (IOCHK1 > 0) THEN

                  FATAL_ERR = FATAL_ERR + 1

                  WRITE(ERR,1252) SETID
                  WRITE(F06,1252) SETID

                  WRITE(ERR,12520) CARD1
                  WRITE(F06,12520) CARD1

                  WRITE(ERR,12524)
                  WRITE(F06,12524)

                  WRITE(F06,*)

                  EXIT

               ENDIF
               READ(TOKEN(3),'(I8)',IOSTAT=IOCHK2) I2
               IF (IOCHK2 > 0) THEN

                  FATAL_ERR = FATAL_ERR + 1

                  WRITE(ERR,1252) SETID
                  WRITE(F06,1252) SETID

                  WRITE(ERR,12520) CARD1
                  WRITE(F06,12520) CARD1

                  WRITE(ERR,12524)
                  WRITE(F06,12524)

                  WRITE(F06,*)

                  EXIT

               ENDIF

! Error: "I1 THRU I2" has I1 > I2

               IF ((IOCHK1 == 0) .AND. (IOCHK2 == 0)) THEN
                  IF (I1 > I2) THEN

                     FATAL_ERR = FATAL_ERR + 1

                     WRITE(ERR,1252) SETID
                     WRITE(F06,1252) SETID

                     WRITE(ERR,12520) CARD1
                     WRITE(F06,12520) CARD1

                     WRITE(ERR,12527) I1,I2
                     WRITE(F06,12527) I1,I2

                     WRITE(F06,*)

                     EXIT

                  ENDIF

               ENDIF

            ELSE

               FATAL_ERR = FATAL_ERR + 1

               WRITE(ERR,1252) SETID
               WRITE(F06,1252) SETID

               WRITE(ERR,12520) CARD1
               WRITE(F06,12520) CARD1

               WRITE(ERR,12528) TOKEN(1),TOKEN(2),TOKEN(3)
               WRITE(F06,12528) TOKEN(1),TOKEN(2),TOKEN(3)

               WRITE(F06,*)

               EXIT

            ENDIF

! Error: wrong tokens from set definition for some other reason

         ELSE

            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1251) SUBR_NAME, NTOKEN
            WRITE(F06,1251) SUBR_NAME, NTOKEN
            CALL OUTA_HERE ( 'Y' )                                 ! Coding error, so quit

            WRITE(F06,*)

            EXIT

         ENDIF

         IF (ISTART <= TOKLEN) THEN
            CYCLE
         ELSE
            EXIT
         ENDIF

      ENDDO



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1251 FORMAT(' *ERROR  1251: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VALUE FOR NOTKEN = ',I8,' NOT RECOGNIZED')

 1252 FORMAT(' *ERROR  1252: SYNTAX ERROR ON THE CASE CONTROL SET ENTRY WITH ID = ',I8,' LISTED BELOW:')

12520 FORMAT(15X,A)

12521 FORMAT(14X,' TOKENS ON SET ENTRIES MUST BE LESS THAN',I4,' CHARACTERS. BELOW IS ONE FOUND THAT EXCEEDS THIS REQUIREMENT:'    &
          ,/,15X,A)

12522 FORMAT(14X,' READ MORE THAN ONE "EXCEPT" BEFORE NEXT "THRU"')

12523 FORMAT(14X,' READ "EXCEPT" BEFORE AN "I1 THRU I2" WAS SPECIFIED')

12524 FORMAT(14X,' ERROR ATTEMPTING TO READ INTEGER ON SET ENTRY')

12525 FORMAT(14X,' THE QUALIFIER "EXCEPT" CAN ONLY MODIFY A PREVIOUS "THRU"')

12526 FORMAT(14X,' EXPECTING TO READ INTEGER TOKEN ON SET ENTRY BUT TOKEN IS: "',A8,'"')

12527 FORMAT(14X,' WHEN "I1 THRU I2" IS USED ON SET ENTRY, I1 MUST BE GREATER THEN I2. INPUT WAS I1 = ',I8,' AND I2 = ',I8)

12528 FORMAT(14X,' EXPECTING TO READ "I1 THRU I2" ON SET ENTRY BUT READ: ',3A8)

12530 FORMAT(14X,' EXPECTING TO READ 3 TOKENS: "I1 THRU I2", BUT ONLY ',I8,' TOKEN WAS FOUND:'                                     &
                   ,/,14X,'"',A8,'"')

12531 FORMAT(14X,' EXPECTING TO READ 3 TOKENS: "I1 THRU I2", BUT ONLY ',I8,' TOKENS WERE FOUND:'                                   &
                   ,/,14X,'"',A8,'" and "',A8,'"')

12532 FORMAT(14X,' THE LAST NONBLANK CHARACTER PRIOR TO ANY COMMENT (BEGINNING WITH $) CANNOT BE A COMMA')

 1267 FORMAT(' *ERROR  1267: FORMAT ERROR ON CASE CONTROL SET ENTRY. COULD NOT FIND EQUAL (=) SIGN')

 1268 FORMAT(' *ERROR  1268: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MUCH SET DATA. LIMIT IS ',I8,' CHARACTERS.')

 1269 FORMAT(' *ERROR  1269: ERROR READING CASE CONTROL SET CONTINUATION ENTRY. ENTRY IGNORED')

 1272 FORMAT(' *ERROR  1272: SET ID ON CASE CONTROL SET ENTRY CANNOT BE HAVE MORE THAN ',I4,' DIGITS')

 1273 FORMAT(' *ERROR  1273: THE SET ID ON ABOVE CASE CONTROL SET ENTRY IS NOT AN INTEGER')

 1145 FORMAT(' *ERROR  1145: DUPLICATE SET ENTRY WITH ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE CC_SET


      SUBROUTINE CC_SUBC ( CARD )

! Processes Case Control SUBCASE cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  CC_ENTRY_LEN, FATAL_ERR, LSUB, NSUB, BLNK_SUB_NAM, NUM_SUBC_CARDS
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  SCNUM

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE INPUT_FILE_MECHANICS, ONLY  :  CSHIFT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CC_SUBC'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=LEN(CARD))        :: CARD1             ! CARD shifted to begin in col after "=" sign
      CHARACTER(1*BYTE)               :: EQUAL_SIGN        ! 'Y' if CARD has an = sign between SUBCASE and the subcase number

      INTEGER(LONG)                   :: ECOL              ! Col, on CARD, where "E" is located (last letter of SUBCASE)
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! Output from subr CSHIFT indicating an error
      INTEGER(LONG)                   :: JERR              ! Error indicator if this subcase number is the same as a previous one
      INTEGER(LONG)                   :: SUBCASE_NUM       ! Subcase number from the SUBCASE card being read




! **********************************************************************************************************************************
! Process SUBCASE cards
      NUM_SUBC_CARDS = NUM_SUBC_CARDS + 1
      NSUB = NSUB + 1
      IF (NSUB > LSUB) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1279) SUBR_NAME,LSUB
         WRITE(F06,1279) SUBR_NAME,LSUB
         CALL OUTA_HERE ( 'Y' )                              ! Coding error, so quit
         RETURN
      ENDIF

! Is there an = sign?

      EQUAL_SIGN = 'N'
      DO I=1,CC_ENTRY_LEN
         IF (CARD(I:I) == '=') THEN
            EQUAL_SIGN = 'Y'
            EXIT
         ENDIF
      ENDDO

      JERR = 0

! There should be no problem finding the 'E' in SUBCASE, subr LOADC wouldn't have called this subr otherwise

      IF (EQUAL_SIGN == 'Y') THEN
         CALL CSHIFT ( CARD, '=', CARD1, ECOL, IERR )
      ELSE
         CALL CSHIFT ( CARD, 'E', CARD1, ECOL, IERR )
      ENDIF
      IF (IERR == 0) THEN
         READ(CARD1,'(I8)') SUBCASE_NUM
         DO I=1,NSUB-1
            IF (SUBCASE_NUM == SCNUM(I)) THEN
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
               WRITE(ERR,1277) SUBCASE_NUM
               WRITE(F06,1277) SUBCASE_NUM
               EXIT
            ENDIF
         ENDDO
         IF (JERR == 0) THEN
            SCNUM(NSUB) = SUBCASE_NUM
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1277 FORMAT(' *ERROR  1277: SUBCASE NUMBER ',I8,' IS A DUPLICATE SUBCASE NUMBER')

 1279 FORMAT(' *ERROR  1279: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY SUBCASES; LIMIT = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE CC_SUBC

   END MODULE CASE_CONTROL_SETS
