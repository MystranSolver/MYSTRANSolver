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

   MODULE BDF_SET_SYNTAX

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: STOKEN, TOKCHK, GET_SETID, GET_ANSID

   CONTAINS

      SUBROUTINE STOKEN ( CALLING_SUBR, TOKSTR, TOKEN_BEG, STRNG_END, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )

! Routine to tokenize character string data. The most common use is in extracting tokens from  the data contained in
! a SET Case Control command that can be of the form:

!         I1, I2, I3 THRU I4 EXCEPT I6, I7, ....

! On repeated calls to STOKEN the routine will find, either:

!  1) a single character token of max length MAX_TOKEN_LEN, or
!  2) a triad of char tokens of the form I1 THRU I2 where I1, I2 are integers

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  MAX_TOKEN_LEN, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)) :: SUBR_NAME = 'STOKEN'
      CHARACTER(LEN=*)          , INTENT(IN)   :: CALLING_SUBR! Character string to tokenize
      CHARACTER(LEN=*)          , INTENT(IN)   :: TOKSTR      ! Character string to tokenize
      CHARACTER( 3*BYTE)        , INTENT(INOUT):: EXCEPT      ! Flag indicating whether EXCEPT is "ON " or "OFF"
      CHARACTER( 3*BYTE)        , INTENT(INOUT):: THRU        ! Flag indicating whether THRU   is "ON " or "OFF"
      CHARACTER(LEN=LEN(TOKSTR)), INTENT(OUT)  :: ERRTOK      ! Char string with data for an error to be printed by calling subr
      CHARACTER( 8*BYTE), INTENT(OUT)          :: TOKEN(3)    ! Array of 3 char tokens (e.g. could contain I1, THRU, I2)
      CHARACTER( 8*BYTE), INTENT(OUT)          :: TOKTYP(3)   ! Array of 3 char indicators of what type of tokens are in TOKEN(1-3)

      INTEGER(LONG), INTENT(IN)                :: STRNG_END   ! Column of last character in TOKSTR
      INTEGER(LONG), INTENT(INOUT)             :: TOKEN_BEG   ! On entry, where to start to look for a token in TOKSTR
!                                                               During processing, it is where the current token starts in TOKSTR
!                                                               On return, it is the start of the next token in TOKSTR.

      INTEGER(LONG), INTENT(OUT)               :: IERROR      ! Integer error no. when an error occurs when processing tokens
      INTEGER(LONG), INTENT(OUT)               :: NTOKEN      ! The number of tokens found in this execution
      INTEGER(LONG)                            :: I           ! DO loop index
      INTEGER(LONG)                            :: TOKEN_END   ! Where, in TOKSTR, the end of the current token is located
      INTEGER(LONG)                            :: NUM_TOK_EXP ! No. of tokens we expect (if we find "THRU", we should find 3 tokens)
      INTEGER(LONG)                            :: PRINT_ITEM  ! An item number to print when DEBUG(19) is turned on




! **********************************************************************************************************************************
! Initialize outputs

      NTOKEN = 0
      IERROR = 0
      ERRTOK(1:) = ' '
      DO I=1,3
         TOKEN(I)(1:) = ' '
         TOKTYP(I)    = 'BLANK   '
      ENDDO
      IF (DEBUG(19) == 1) THEN                             ! Print debug data
         PRINT_ITEM = 0
         CALL DEB_STOKEN ( PRINT_ITEM )
      ENDIF

      NUM_TOK_EXP = 1                                      ! Initialize variables

      IF (DEBUG(19) == 1) THEN                             ! Print debug data
         PRINT_ITEM = 1
         CALL DEB_STOKEN ( PRINT_ITEM )
         WRITE(F06,*)
         WRITE(F06,*) '    Begin outer:DO'
         WRITE(F06,*) '    --------------'
      ENDIF

      DO I=TOKEN_BEG,STRNG_END                             ! Make sure we are positioned at beginning of a token
         IF ((TOKSTR(I:I) == ' ') .OR. (TOKSTR(I:I) == ',')) THEN
            TOKEN_BEG = TOKEN_BEG + 1
         ELSE
            EXIT
         ENDIF
      ENDDO

! Top of loop for processing tokens. If there is only 1 token, we will execute the loop once.
! If there is a triad ("I1 THRU I2") we will process the loop 3 times

outer:DO

         TOKEN_END = STRNG_END                             ! Find end of current token
i_loop1: DO I = TOKEN_BEG,STRNG_END
            IF ((TOKSTR(I:I) == ' ') .OR. (TOKSTR(I:I) == ',')) THEN
               TOKEN_END = I-1
               EXIT i_loop1
            ENDIF
         ENDDO i_loop1

         NTOKEN = NTOKEN + 1
                                                           ! Check for token too long and set error flag if so
         IF ((TOKEN_END - TOKEN_BEG + 1) > MAX_TOKEN_LEN) THEN
            IERROR                         = 1
            ERRTOK(1: )                    = TOKSTR(TOKEN_BEG:TOKEN_END)
            TOKEN(NTOKEN)(1:MAX_TOKEN_LEN) = TOKSTR(TOKEN_BEG:TOKEN_BEG+MAX_TOKEN_LEN-1)
         ELSE
            TOKEN(NTOKEN)(1:) = TOKSTR(TOKEN_BEG:TOKEN_END)
         ENDIF

         CALL TOKCHK ( TOKEN(NTOKEN), TOKTYP(NTOKEN) )     ! Process TOKEN

         IF (DEBUG(19) == 1) THEN                          ! Print debug data
            PRINT_ITEM = 2
            CALL DEB_STOKEN ( PRINT_ITEM )
         ENDIF

         IF (TOKEN(NTOKEN) == 'EXCEPT  ') THEN             ! Set EXCEPT flag if TOKEN = 'EXCEPT  ', but only if
!                                                            THRU is 'ON ' & EXCEPT is not already 'ON '.
            IF ((THRU == 'ON ') .AND. (EXCEPT == 'OFF')) THEN
               EXCEPT = 'ON '

            ELSE                                           ! Otherwise, set error number and exit outer

               IF      ((THRU == 'ON ') .AND. (EXCEPT == 'ON ')) THEN
                  IERROR = 3
               ELSE IF ((THRU == 'OFF') .AND. (EXCEPT == 'ON ')) THEN
                  IERROR = 4
               ELSE IF ((THRU == 'OFF') .AND. (EXCEPT == 'OFF')) THEN
                  IERROR = 4
               ENDIF

               EXIT outer                                  ! Exit outer due to errors

            ENDIF

         ENDIF

         IF (DEBUG(19) == 1) THEN                          ! Print debug data
            PRINT_ITEM = 3
            CALL DEB_STOKEN ( PRINT_ITEM )
         ENDIF

         IF (TOKEN_END <= STRNG_END) THEN

            TOKEN_BEG = STRNG_END + 1                      ! Find start of next token. Default it to STRNG_END+1 to start with
i_loop2:    DO I = TOKEN_END+1,STRNG_END                   ! just in case we are already at the end of TOKSTR.
               IF ((TOKSTR(I:I) == ' ') .OR. (TOKSTR(I:I) == ',')) THEN
                  CYCLE i_loop2
               ELSE
                  TOKEN_BEG = I                            ! This is where the next token starts
                  EXIT i_loop2
               ENDIF
            ENDDO i_loop2

            IF (TOKEN_BEG <= STRNG_END) THEN

               IF ((NTOKEN == 1) .AND. ((TOKEN_BEG+3) <= STRNG_END)) THEN
                                                           ! Peek ahead to see if 2nd token is "THRU"
                  IF (TOKSTR(TOKEN_BEG:TOKEN_BEG+3)=='THRU') THEN! We found "THRU"
                     NUM_TOK_EXP = 3
                     THRU   = 'ON '
                     EXCEPT = 'OFF'
                     TOKEN_END = TOKEN_BEG + 3

                     IF (DEBUG(19) == 1) THEN
                        PRINT_ITEM = 4
                        CALL DEB_STOKEN ( PRINT_ITEM )
                        WRITE(F06,*)
                        WRITE(F06,*) '    Cycling to top of outer:DO'
                        WRITE(F06,*) '    --------------------------'
                     ENDIF

                     CYCLE outer                           ! CYCLE to read/process 2nd token

                  ELSE

                     EXIT outer                            ! NTOKEN was 1 but we didn't find "THRU" next, so exit outer

                  ENDIF

               ELSE IF ((NTOKEN == 2) .AND. (NUM_TOK_EXP == 3)) THEN
                                                           ! Need to go back for 3rd token since we know we expect 3
                  IF (DEBUG(19) == 1) THEN
                     PRINT_ITEM = 5
                     CALL DEB_STOKEN ( PRINT_ITEM )
                     WRITE(F06,*)
                     WRITE(F06,*) '    Cycling to top of outer:DO'
                     WRITE(F06,*) '    --------------------------'
                  ENDIF

                  CYCLE outer                              ! CYCLE to read/process 3nd token

               ELSE

                  EXIT outer                               ! Need to exit outer (NTOKEN >= 3, etc)

               ENDIF

            ELSE

               EXIT outer                                  ! Need to exit outer (TOKEN_BEG > STRNG_END)

            ENDIF

         ELSE                                              ! TOKEN_END >= STRNG_END, so set TOKEN_BEG = STRNG_END + 1 and quit

            TOKEN_BEG = STRNG_END + 1

            IF ((NUM_TOK_EXP == 3) .AND. (NTOKEN < 3)) THEN
               IERROR = 2                                  ! Error: we found 'THRU' but didnt get 3 tokens
            ENDIF

            EXIT outer                                     ! Need to exit outer (TOKEN_END > STRNG_END)

         ENDIF

      ENDDO outer

      IF (DEBUG(19) == 1) THEN                             ! Print debug data following loop
         PRINT_ITEM = 6
         CALL DEB_STOKEN ( PRINT_ITEM )
      ENDIF



      RETURN

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEB_STOKEN ( PRINT_ITEM )

      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEB_STOKEN'

      INTEGER(LONG), INTENT(IN)       :: PRINT_ITEM        ! What item to print




! **********************************************************************************************************************************

      IF      (PRINT_ITEM == 0) THEN
         WRITE(F06,*)
         WRITE(F06,99000) CALLING_SUBR, '"', TOKSTR(1:STRNG_END), '"'
      ELSE IF (PRINT_ITEM == 1) THEN
         WRITE(F06,99001) TOKEN_BEG     ,THRU,EXCEPT,NUM_TOK_EXP,NTOKEN,(TOKTYP(I),I=1,3),(TOKEN(I),I=1,3),IERROR
      ELSE IF (PRINT_ITEM == 2) THEN
         WRITE(F06,99002) TOKEN_BEG,TOKEN_END,THRU,EXCEPT,NUM_TOK_EXP,NTOKEN,(TOKTYP(I),I=1,3),(TOKEN(I),I=1,3),IERROR
      ELSE IF (PRINT_ITEM == 3) THEN
         WRITE(F06,99003) TOKEN_BEG,TOKEN_END,THRU,EXCEPT,NUM_TOK_EXP,NTOKEN,(TOKTYP(I),I=1,3),(TOKEN(I),I=1,3),IERROR
      ELSE IF (PRINT_ITEM == 4) THEN
         WRITE(F06,99004) TOKEN_BEG,TOKEN_END,THRU,EXCEPT,NUM_TOK_EXP,NTOKEN,(TOKTYP(I),I=1,3),(TOKEN(I),I=1,3),IERROR
      ELSE IF (PRINT_ITEM == 5) THEN
         WRITE(F06,99005) TOKEN_BEG,TOKEN_END,THRU,EXCEPT,NUM_TOK_EXP,NTOKEN,(TOKTYP(I),I=1,3),(TOKEN(I),I=1,3),IERROR
      ELSE IF (PRINT_ITEM == 6) THEN
         WRITE(F06,99006) TOKEN_BEG     ,THRU,EXCEPT,NUM_TOK_EXP,NTOKEN,(TOKTYP(I),I=1,3),(TOKEN(I),I=1,3),IERROR
         WRITE(F06,99999)
      ENDIF

99000 FORMAT(' //////////////////////////////////////////////////////////////////////////////////////////////////////////////////',&
              '/////////////////',/,' Subr STOKEN (called by subr ',A,')  has array TOKSTR at the beg. of subr STOKEN =',/,1X,3A,//&
              ,45x,'Progress of subr STOKEN in parsing part of TOKSTR:',/                                                          &
              ,45x,'--------------------------------------------------',/                                                          &
              ,35x,'   TOKEN   THRU EXCEPT   Num Tokens  TOKTYP1  TOKTYP2  TOKTYP3  TOKEN1   TOKEN2   TOKEN3   IERROR',/           &
              ,35X,' ---------              ------------',/                                                                        &
              ,35X,' Beg   End              Expect Found')

99001 FORMAT(/,' (1) Beg STOKEN-before outer:DO:',I6,6X,2X,A3,2X,A3,7X,I2,3X,I2,3X,6(1X,A8),3X,I2)

99002 FORMAT(  ' (2) After calling TOKCHK      :',I6,I6,2X,A3,2X,A3,7X,I2,3X,I2,3X,6(1X,A8),3X,I2)

99003 FORMAT(  ' (3) After EXCEPT flag check   :',I6,I6,2X,A3,2X,A3,7X,I2,3X,I2,3X,6(1X,A8),3X,I2,/)

99004 FORMAT(  ' (4) CYCLE back to outer DO    :',I6,I6,2X,A3,2X,A3,7X,I2,3X,I2,3X,6(1X,A8),3X,I2,/)

99005 FORMAT(  ' (5) CYCLE back to outer DO    :',I6,I6,2X,A3,2X,A3,7X,I2,3X,I2,3X,6(1X,A8),3X,I2,/)

99006 FORMAT(  ' (6) End STOKEN-after outer DO :',I6,6X,2X,A3,2X,A3,7X,I2,3X,I2,3X,6(1X,A8),3X,I2)

99999 FORMAT(' \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\' &
            ,'\\\\\\\\\\\\\\\\\\\\\\\\', //)



      RETURN
! **********************************************************************************************************************************

      END SUBROUTINE DEB_STOKEN

      END SUBROUTINE STOKEN


      SUBROUTINE TOKCHK ( TOKEN, TOKTYPE )

! Determines the type of an 8 character input string, TOKEN. Output TOKTYPE is:

!   TOKTYPE = 'UNKNOWN ' if input TOKEN is none of the below
!   TOKTYPE = 'COMMENT ' if input TOKEN is $
!   TOKTYPE = 'BLANK   ' if input TOKEN is blank
!   TOKTYPE = 'INTEGER ' if input TOKEN is an integer
!   TOKTYPE = 'FL PT   ' if input TOKEN is a floating point number
!   TOKTYPE = 'ALL     ' if input TOKEN is "ALL"
!   TOKTYPE = 'THRU    ' if input TOKEN is "THRU"
!   TOKTYPE = 'NONE    ' if input TOKEN is "NONE"
!   TOKTYPE = 'PRINT   ' if input TOKEN is "PRINT"
!   TOKTYPE = 'PUNCH   ' if input TOKEN is "PUNCH"
!   TOKTYPE = 'BOTH    ' if input TOKEN is "BOTH"
!   TOKTYPE = 'EXCEPT  ' if input TOKEN is "EXCEPT"
!   TOKTYPE = 'FIJFIL  ' if input TOKEN is "FIJFIL"

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06

      IMPLICIT NONE

      CHARACTER(8*BYTE), INTENT(IN)   :: TOKEN             ! The input character string
      CHARACTER(8*BYTE), INTENT(OUT)  :: TOKTYPE           ! The type of TOKEN (see above)

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IEND              ! Col number (1-8) in TOKEN where TOKEN ends
      INTEGER(LONG)                   :: ISTART            ! Col number (1-8) in TOKEN where TOKEN begins

! **********************************************************************************************************************************
      TOKTYPE = 'UNKNOWN '

! Check for blank token

      IF (TOKEN == '        ') THEN
         TOKTYPE = 'BLANK   '
         RETURN
      ENDIF

! Check for comment token

      IF ((TOKEN(1:1) == '$') .OR. (TOKEN(2:2) == '$') .OR. (TOKEN(3:3) == '$') .OR. (TOKEN(4:4) == '$') .OR.                      &
          (TOKEN(5:5) == '$') .OR. (TOKEN(6:6) == '$') .OR. (TOKEN(7:7) == '$') .OR. (TOKEN(8:8) == '$')) THEN
         TOKTYPE = 'COMMENT '
         RETURN
      ENDIF

! Check for an integer. First, skip leading and trailing blanks then read digits until another blank (or field end)
! is encountered.

      IF (TOKTYPE == 'UNKNOWN') THEN

         ISTART = 1
         DO I=1,8
            IF(TOKEN(I:I) == ' ') THEN
               ISTART = ISTART + 1
               CYCLE
            ELSE
               EXIT
            ENDIF
         ENDDO

         IEND = 8
         DO I=8,1,-1
            IF(TOKEN(I:I) == ' ') THEN
               IEND = IEND - 1
               CYCLE
            ELSE
               EXIT
            ENDIF
         ENDDO

         IF ((TOKEN(ISTART:ISTART) == '+') .OR. (TOKEN(ISTART:ISTART) == '-')) THEN
            ISTART = ISTART+1
         ENDIF
         DO I=ISTART,IEND
            IF ((TOKEN(I:I)=='1') .OR. (TOKEN(I:I)=='2') .OR. (TOKEN(I:I)=='3') .OR. (TOKEN(I:I)=='4') .OR.  &
                (TOKEN(I:I)=='5') .OR. (TOKEN(I:I)=='6') .OR. (TOKEN(I:I)=='7') .OR. (TOKEN(I:I)=='8') .OR.  &
                (TOKEN(I:I)=='9') .OR. (TOKEN(I:I)=='0')) THEN
               TOKTYPE = 'INTEGER '
            ELSE
               TOKTYPE = 'UNKNOWN '
               EXIT
            ENDIF
         ENDDO
         IF (TOKTYPE == 'INTEGER ') THEN
            RETURN
         ENDIF

      ENDIF

! Check for floating point number

      IF (TOKTYPE == 'UNKNOWN') THEN                       ! First find a '.'
         DO I = 1,8
            IF (TOKEN(I:I) == '.') THEN
               TOKTYPE = '?FL PT? '
               EXIT
            ENDIF
         ENDDO

         IF (TOKTYPE == '?FL PT? ') THEN                   ! Check characters for floating point number
            CONTINUE
            DO I = 1,8
               IF ((TOKEN(I:I)=='1') .OR. (TOKEN(I:I)=='2') .OR. (TOKEN(I:I)=='3') .OR. (TOKEN(I:I)=='4') .OR.  &
                   (TOKEN(I:I)=='5') .OR. (TOKEN(I:I)=='6') .OR. (TOKEN(I:I)=='7') .OR. (TOKEN(I:I)=='8') .OR.  &
                   (TOKEN(I:I)=='9') .OR. (TOKEN(I:I)=='0') .OR. (TOKEN(I:I)=='+') .OR. (TOKEN(I:I)=='-') .OR.  &
                   (TOKEN(I:I)=='D') .OR. (TOKEN(I:I)=='E') .OR. (TOKEN(I:I)=='.') .OR. (TOKEN(I:I)==' ')) THEN
                   TOKTYPE = 'FL PT   '
                   CYCLE
               ELSE
                  TOKTYPE = 'UNKNOWN '
                  EXIT
               ENDIF
            ENDDO
            IF (TOKTYPE == 'FL PT   ') THEN
               RETURN
            ENDIF
         ENDIF

      ENDIF

! Check for 'ALL'

      IF (TOKTYPE == 'UNKNOWN') THEN
         IF               ((TOKEN(1:1)=='A') .OR. (TOKEN(1:1)=='a')) THEN
            IF            ((TOKEN(2:2)=='L') .OR. (TOKEN(2:2)=='l')) THEN
               IF         ((TOKEN(3:3)=='L') .OR. (TOKEN(3:3)=='l')) THEN
                  TOKTYPE = 'ALL     '
                  RETURN
               ENDIF
            ENDIF
         ENDIF
      ENDIF

! Check for 'THRU'

      IF (TOKTYPE == 'UNKNOWN') THEN
         IF               ((TOKEN(1:1)=='T') .OR. (TOKEN(1:1)=='t')) THEN
            IF            ((TOKEN(2:2)=='H') .OR. (TOKEN(2:2)=='h')) THEN
               IF         ((TOKEN(3:3)=='R') .OR. (TOKEN(3:3)=='r')) THEN
                  IF      ((TOKEN(4:4)=='U') .OR. (TOKEN(4:4)=='u')) THEN
                     TOKTYPE = 'THRU    '
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

! Check for 'NONE'

      IF (TOKTYPE == 'UNKNOWN') THEN
         IF               ((TOKEN(1:1)=='N') .OR. (TOKEN(1:1)=='n')) THEN
            IF            ((TOKEN(2:2)=='O') .OR. (TOKEN(2:2)=='o')) THEN
               IF         ((TOKEN(3:3)=='N') .OR. (TOKEN(3:3)=='n')) THEN
                  IF      ((TOKEN(4:4)=='E') .OR. (TOKEN(4:4)=='e')) THEN
                     TOKTYPE = 'NONE    '
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

! Check for 'PRINT'

      IF (TOKTYPE == 'UNKNOWN') THEN
         IF               ((TOKEN(1:1)=='P') .OR. (TOKEN(1:1)=='p')) THEN
            IF            ((TOKEN(2:2)=='R') .OR. (TOKEN(2:2)=='r')) THEN
               IF         ((TOKEN(3:3)=='I') .OR. (TOKEN(3:3)=='i')) THEN
                  IF      ((TOKEN(4:4)=='N') .OR. (TOKEN(4:4)=='n')) THEN
                     IF   ((TOKEN(5:5)=='T') .OR. (TOKEN(5:5)=='t')) THEN
                        TOKTYPE = 'PRINT   '
                        RETURN
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

! Check for 'PUNCH'

      IF (TOKTYPE == 'UNKNOWN') THEN
         IF               ((TOKEN(1:1)=='P') .OR. (TOKEN(1:1)=='p')) THEN
            IF            ((TOKEN(2:2)=='U') .OR. (TOKEN(2:2)=='u')) THEN
               IF         ((TOKEN(3:3)=='N') .OR. (TOKEN(3:3)=='n')) THEN
                  IF      ((TOKEN(4:4)=='C') .OR. (TOKEN(4:4)=='c')) THEN
                     IF   ((TOKEN(5:5)=='H') .OR. (TOKEN(5:5)=='h')) THEN
                        TOKTYPE = 'PUNCH   '
                        RETURN
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

! Check for 'BOTH'

      IF (TOKTYPE == 'UNKNOWN') THEN
         IF               ((TOKEN(1:1)=='B') .OR. (TOKEN(1:1)=='b')) THEN
            IF            ((TOKEN(2:2)=='O') .OR. (TOKEN(2:2)=='o')) THEN
               IF         ((TOKEN(3:3)=='T') .OR. (TOKEN(3:3)=='t')) THEN
                  IF      ((TOKEN(4:4)=='H') .OR. (TOKEN(4:4)=='h')) THEN
                     TOKTYPE = 'BOTH    '
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

! Check for 'EXCEPT'

      IF (TOKTYPE == 'UNKNOWN') THEN
         IF               ((TOKEN(1:1)=='E') .OR. (TOKEN(1:1)=='e')) THEN
            IF            ((TOKEN(2:2)=='X') .OR. (TOKEN(2:2)=='x')) THEN
               IF         ((TOKEN(3:3)=='C') .OR. (TOKEN(3:3)=='c')) THEN
                  IF      ((TOKEN(4:4)=='E') .OR. (TOKEN(4:4)=='e')) THEN
                     IF   ((TOKEN(5:5)=='P') .OR. (TOKEN(5:5)=='p')) THEN
                        IF((TOKEN(6:6)=='T') .OR. (TOKEN(6:6)=='t')) THEN
                           TOKTYPE = 'EXCEPT  '
                           RETURN
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

! Check for 'FIJFIL'

      IF (TOKTYPE == 'UNKNOWN') THEN
         IF               ((TOKEN(1:1)=='F') .OR. (TOKEN(1:1)=='f')) THEN
            IF            ((TOKEN(2:2)=='I') .OR. (TOKEN(2:2)=='i')) THEN
               IF         ((TOKEN(3:3)=='J') .OR. (TOKEN(3:3)=='j')) THEN
                  IF      ((TOKEN(4:4)=='F') .OR. (TOKEN(4:4)=='f')) THEN
                     IF   ((TOKEN(5:5)=='I') .OR. (TOKEN(5:5)=='i')) THEN
                        IF((TOKEN(6:6)=='L') .OR. (TOKEN(6:6)=='l')) THEN
                           TOKTYPE = 'FIJFIL  '
                           RETURN
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE TOKCHK


      SUBROUTINE GET_SETID ( CARD, SETID )

! Gets SET ID from CASE CONTROL cards:  LOAD, METHOD, MPC, NLPARM, SPC, TEMP

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  CC_ENTRY_LEN, FATAL_ERR, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC


      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_SETID'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Case Control card (can be modified by subr CSHIFT, called herein)
      CHARACTER(LEN=LEN(CARD)+1)      :: CARD1             ! CARD shifted to begin in col after "=" sign
      CHARACTER(LEN=LEN(CARD)+1)      :: ERRTOK            ! An output from subr STOKEN, called herein
      CHARACTER( 3*BYTE)              :: EXCEPT            ! An input/output to/from subr STOKEN, called herein
      CHARACTER( 3*BYTE)              :: THRU              ! An inputoutput to/from subr STOKEN, called herein
      CHARACTER( 8*BYTE)              :: TOKEN(3)          ! An output from subr STOKEN, called herein
      CHARACTER( 8*BYTE)              :: TOKTYP(3)         ! An output from subr STOKEN, called herein

      INTEGER(LONG), INTENT(OUT)      :: SETID             ! Set ID read from CARD after '=', if CARD contains an integer here.
!                                                            SETID is set to -1 if 'ALL' is found after '=' & to 0 if 'NONE' found
      INTEGER(LONG)                   :: ECOL              ! Column on CARD where '=' is located
      INTEGER(LONG)                   :: IERR              ! An output from subr CSHIFT, called herein
      INTEGER(LONG)                   :: IERROR            ! An output from subr STOKEN, called herein
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from READ
      INTEGER(LONG)                   :: ISTART            ! An input to subr STOKEN, called herein
      INTEGER(LONG)                   :: NTOKEN            ! An output from subr STOKEN, called herein
      INTEGER(LONG)                   :: TOKLEN            ! An input to subr STOKEN, called herein




! **********************************************************************************************************************************
! Get SETID

      SETID = 0
      CALL CSHIFT ( CARD, '=', CARD1, ECOL, IERR )
      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1703)
         WRITE(F06,1703)
      ELSE
         ISTART = 1
         THRU   = 'OFF'
         EXCEPT = 'OFF'
         TOKLEN = CC_ENTRY_LEN
         CALL STOKEN ( SUBR_NAME, CARD1, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )
         IF (NTOKEN > 1) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1704) NTOKEN
            WRITE(F06,1704) NTOKEN
         ELSE
            IF ((ISTART <= TOKLEN) .AND. (CARD1(ISTART:ISTART) /= '$')) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1704) NTOKEN
               WRITE(F06,1704) NTOKEN
            ELSE
               IF (TOKTYP(1) == 'INTEGER ') THEN
                  READ(TOKEN(1),'(I8)',IOSTAT=IOCHK) SETID
                  IF ((IOCHK < 0) .OR. (IOCHK > 0)) THEN
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1705)
                     WRITE(F06,1705)
                  ENDIF
                  IF (SETID <= 0) THEN
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1706)
                     WRITE(F06,1706)
                  ENDIF
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1707)
                  WRITE(F06,1707)
               ENDIF
            ENDIF
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1703 FORMAT(' *ERROR  1703: CANNOT FIND EQUAL SIGN (=) ON ABOVE CASE CONTROL ENTRY ')

 1704 FORMAT(' *ERROR  1704: ABOVE CASE CONTROL ENTRY MUST HAVE ONLY ONE INTEGER SET NUMBER FOLLOWING THE "=" SIGN.'               &
                    ,/,14X,' ABOVE ENTRY DOES NOT MEET THIS REQUIREMENT (NTOKEN = ', I8, ')')

 1705 FORMAT(' *ERROR  1705: ERROR READING SET ID ON PREVIOUS CASE CONTROL ENTRY')

 1706 FORMAT(' *ERROR  1706: ZERO OR NEGATIVE SET ID NOT ALLOWED ON PREVIOUS CASE CONTROL ENTRY')

 1707 FORMAT(' *ERROR  1707: SET ID MUST BE AN INTEGER OF <= 8 DIGITS')

! **********************************************************************************************************************************

      END SUBROUTINE GET_SETID


      SUBROUTINE GET_ANSID ( CARD, SETID )

! Gets 'ALL', 'NONE' or set ID from Case Control cards:  DISP, ELDATA, ELFORCE, GPFORCE, OLOAD, SPCFORCE, STRESS

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  CC_ENTRY_LEN, FATAL_ERR, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC


      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_ANSID'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Case Control card (can be modified by subr CSHIFT, called herein)
      CHARACTER(LEN=LEN(CARD)+1)      :: CARD1             ! CARD shifted to begin in col after "=" sign
      CHARACTER(LEN=LEN(CARD)+1)      :: ERRTOK            ! An output from subr STOKEN, called herein
      CHARACTER( 3*BYTE)              :: EXCEPT            ! An input/output to/from subr STOKEN, called herein
      CHARACTER( 3*BYTE)              :: THRU              ! An inputoutput to/from subr STOKEN, called herein
      CHARACTER( 8*BYTE)              :: TOKEN(3)          ! An output from subr STOKEN, called herein
      CHARACTER( 8*BYTE)              :: TOKTYP(3)         ! An output from subr STOKEN, called herein

      INTEGER(LONG), INTENT(OUT)      :: SETID             ! Set ID read from CARD after '=', if CARD contains an integer here.
!                                                            SETID is set to -1 if 'ALL' is found after '=' & to 0 if 'NONE' found
      INTEGER(LONG)                   :: ECOL              ! Column on CARD where '=' is located
      INTEGER(LONG)                   :: IERR              ! An output from subr CSHIFT, called herein
      INTEGER(LONG)                   :: IERROR            ! An output from subr STOKEN, called herein
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from READ
      INTEGER(LONG)                   :: ISTART            ! An input to subr STOKEN, called herein
      INTEGER(LONG)                   :: NTOKEN            ! An output from subr STOKEN, called herein
      INTEGER(LONG)                   :: TOKLEN            ! An input to subr STOKEN, called herein




! **********************************************************************************************************************************
! Get 'ALL', 'NONE' or SETID

      SETID = 0
      CALL CSHIFT ( CARD, '=', CARD1, ECOL, IERR )
      IF (IERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1703)
         WRITE(F06,1703)
      ELSE
         ISTART = 1
         THRU   = 'OFF'
         EXCEPT = 'OFF'
         TOKLEN = CC_ENTRY_LEN
         CALL STOKEN ( SUBR_NAME, CARD1, ISTART, TOKLEN, NTOKEN, IERROR, TOKTYP, TOKEN, ERRTOK, THRU, EXCEPT )
         IF (NTOKEN > 1) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1704) NTOKEN
            WRITE(F06,1704) NTOKEN
         ELSE
            IF ((ISTART <= TOKLEN) .AND. (CARD1(ISTART:ISTART) /= '$')) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1704) NTOKEN
               WRITE(F06,1704) NTOKEN
            ELSE
               IF      (TOKTYP(1) == 'INTEGER ') THEN
                  READ(TOKEN(1),'(I8)',IOSTAT=IOCHK) SETID
                  IF ((IOCHK < 0) .OR. (IOCHK > 0)) THEN
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1705)
                     WRITE(F06,1705)
                  ENDIF
                  IF (SETID <= 0) THEN
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1706)
                     WRITE(F06,1706)
                  ENDIF
               ELSE IF (TOKTYP(1) == 'ALL     ') THEN
                  SETID = -1
               ELSE IF (TOKTYP(1) == 'NONE    ') THEN
                  SETID = 0
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1708)
                  WRITE(F06,1708)
               ENDIF
            ENDIF
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1703 FORMAT(' *ERROR  1703: CANNOT FIND EQUAL SIGN (=) ON ABOVE CASE CONTROL ENTRY ')

 1704 FORMAT(' *ERROR  1704: ABOVE CASE CONTROL ENTRY MUST HAVE ONLY ONE INTEGER SET NUMBER FOLLOWING THE "=" SIGN.'               &
                    ,/,14X,' ABOVE ENTRY DOES NOT MEET THIS REQUIREMENT')

 1705 FORMAT(' *ERROR  1705: ERROR READING SET ID ON PREVIOUS CASE CONTROL ENTRY')

 1706 FORMAT(' *ERROR  1706: ZERO OR NEGATIVE SET ID NOT ALLOWED ON PREVIOUS CASE CONTROL ENTRY')

 1708 FORMAT(' *ERROR  1708: SET ID ON ABOVE CASE CONTROL ENTRY MUST BE: "ALL", "NONE", OR AN INTEGER SET ID OF <= 8 DIGITS')

! **********************************************************************************************************************************

      END SUBROUTINE GET_ANSID

   END MODULE BDF_SET_SYNTAX
