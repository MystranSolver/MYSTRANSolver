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

      SUBROUTINE BD_PLOAD2 ( CARD, CC_LOAD_FND )

! Processes PLOAD2 Bulk Data Cards. Reads and checks data and then writes CARD to file LINK1Q for later processing

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1Q
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPLOAD, LSUB, NPCARD, NPLOAD, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PRESS_SIDS, SUBLOD

      USE BD_PLOAD2_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PLOAD2'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD               ! A Bulk Data card
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_LOAD_FND(LSUB,2)! 'Y' if B.D load/temp card w/ same set ID (SID) as C.C. LOAD = SID
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)          ! The 10 fields of characters making up CARD
! --- cbeam_add begin --- !
      CHARACTER(LEN=JCARD_LEN)        :: PLOAD1_SCALE       ! PLOAD1 scale keyword
      CHARACTER(LEN=JCARD_LEN)        :: PLOAD1_TYPE        ! PLOAD1 load type keyword
! --- cbeam_add end --- !
      CHARACTER( 1*BYTE)              :: THRU               ! 'Y' if field 5 of parent card is "THRU"
      CHARACTER( 8*BYTE)              :: TOKEN              ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP             ! The type of token in a field of parent card. Output from subr TOKCHK

      INTEGER(LONG)                   :: PLOAD_ELID(6)      ! Elem ID's on parent card if "THRU" not used for input
      INTEGER(LONG)                   :: J                  ! DO loop index
      INTEGER(LONG)                   :: JERR               ! Error count
      INTEGER(LONG)                   :: SETID              ! Load set ID on PLOADi card


      REAL(DOUBLE)  :: RPRESS



! **********************************************************************************************************************************
! PLOAD2 Bulk Data card check

!   FIELD   ITEM
!   -----   ------------
!    2      SID
!    3      Pressure
!    4-9    Element ID's (PLOAD_ELID's)
! or:
!    4-6    ELID1 THRU ELID2


! Make JCARD from CARD

      JERR = 0
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow and increment NPLOAD

      NPLOAD = NPLOAD+1

! Check if load set ID on pressure card matches a Case Control request

      CALL I4FLD ( JCARD(2), JF(2), SETID )
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NSUB
            IF (SETID == SUBLOD(J,1)) THEN
               CC_LOAD_FND(J,1) = 'Y'
            ENDIF
         ENDDO
         PRESS_SIDS(NPLOAD) = SETID
      ELSE
         JERR = JERR + 1
      ENDIF

! Check for the 2 options on specifying data on the card. Either all data are PLOAD_ELID's or the THRU  option is used in
! which case field 5 will have "THRU".


! --- cbeam_add begin --- !
      IF      ((JCARD(1)(1:7) == 'PLOAD1 ') .OR. (JCARD(1)(1:7) == 'PLOAD1*')) THEN
         CALL I4FLD ( JCARD(3), JF(3), PLOAD_ELID(1) )
         IF (IERRFL(3) == 'N') THEN
            IF (PLOAD_ELID(1) <= 0) THEN
               JERR      = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1152) JCARD(1), JCARD(2)
               WRITE(F06,1152) JCARD(1), JCARD(2)
            ENDIF
         ELSE
            JERR = JERR + 1
         ENDIF

         CALL CHAR_FLD ( JCARD(4), JF(4), PLOAD1_TYPE )
         CALL LEFT_ADJ_BDFLD ( PLOAD1_TYPE )
         IF (PLOAD1_TYPE(1:2) == 'MX') PLOAD1_TYPE = 'MXE'
         IF (PLOAD1_TYPE(1:2) == 'MY') PLOAD1_TYPE = 'MYE'
         IF (PLOAD1_TYPE(1:2) == 'MZ') PLOAD1_TYPE = 'MZE'
         IF ((PLOAD1_TYPE /= 'FY' ) .AND. (PLOAD1_TYPE /= 'FYE') .AND. (PLOAD1_TYPE /= 'FZ' ) .AND. (PLOAD1_TYPE /= 'FZE') .AND.  &
             (PLOAD1_TYPE /= 'Y'  ) .AND. (PLOAD1_TYPE /= 'Z'  ) .AND. (PLOAD1_TYPE /= 'FXE') .AND. (PLOAD1_TYPE /= 'MXE') .AND.  &
             (PLOAD1_TYPE /= 'MYE') .AND. (PLOAD1_TYPE /= 'MZE')) THEN
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1199) JCARD(1), JCARD(2), PLOAD1_TYPE
            WRITE(F06,1199) JCARD(1), JCARD(2), PLOAD1_TYPE
         ENDIF

         CALL CHAR_FLD ( JCARD(5), JF(5), PLOAD1_SCALE )
         CALL LEFT_ADJ_BDFLD ( PLOAD1_SCALE )
         IF (PLOAD1_SCALE == ' ') PLOAD1_SCALE = 'FR'
         IF ((PLOAD1_SCALE /= 'FR'  ) .AND. (PLOAD1_SCALE /= 'LE'  ) .AND.                                             &
             (PLOAD1_SCALE /= 'FRPR') .AND. (PLOAD1_SCALE /= 'LEPR')) THEN
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1200) JCARD(1), JCARD(2), PLOAD1_SCALE
            WRITE(F06,1200) JCARD(1), JCARD(2), PLOAD1_SCALE
         ENDIF

         CALL R8FLD ( JCARD(6), JF(6), RPRESS )
         CALL R8FLD ( JCARD(7), JF(7), RPRESS )
         CALL R8FLD ( JCARD(8), JF(8), RPRESS )
         CALL R8FLD ( JCARD(9), JF(9), RPRESS )

         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
         CALL CRDERR ( CARD )
! --- cbeam_add end --- !
      ELSE IF ((JCARD(1)(1:7) == 'PLOAD2 ') .OR. (JCARD(1)(1:7) == 'PLOAD2*')) THEN

         CALL R8FLD ( JCARD(3), JF(3), RPRESS )

         THRU = 'N'
         TOKEN = JCARD(5)(1:8)                             ! Only send the 1st 8 chars of this JCARD. It has been left justified
         CALL TOKCHK ( TOKEN, TOKTYP )

         IF (TOKTYP == 'THRU    ') THEN
            THRU = 'Y'
         ENDIF

         IF (THRU == 'N') THEN
            DO J=4,9
               IF (JCARD(J)(1:) == ' ') EXIT
               CALL I4FLD ( JCARD(J), JF(J), PLOAD_ELID(J-3) )
               IF (IERRFL(J) == 'N') THEN
                  IF (PLOAD_ELID(J-3) <= 0) THEN
                     JERR      = JERR + 1
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1152) JCARD(1), JCARD(2)
                     WRITE(F06,1152) JCARD(1), JCARD(2)
                  ENDIF
               ELSE
                  JERR = JERR + 1
               ENDIF
            ENDDO
         ELSE
            CALL I4FLD ( JCARD(4), JF(4), PLOAD_ELID(1) )
            CALL I4FLD ( JCARD(6), JF(6), PLOAD_ELID(2) )
            IF ((IERRFL(4) == 'N') .AND. (IERRFL(4) == 'N')) THEN
               IF ((PLOAD_ELID(2) < PLOAD_ELID(1)) .OR. (PLOAD_ELID(1) <= 0) .OR. (PLOAD_ELID(2) <= 0)) THEN
                  JERR      = JERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1128) JCARD(1), JCARD(2)
                  WRITE(F06,1128) JCARD(1), JCARD(2)
               ENDIF
            ELSE
               JERR = JERR + 1
            ENDIF
         ENDIF

         IF (THRU == 'N') THEN
            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
         ELSE
            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-6
            CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,7,8,9 )! Issue warning if fieldS 7, 8, 9 not blank
         ENDIF
         CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      ENDIF


! Write data to file L1Q

      IF (JERR == 0) THEN
         WRITE(L1Q) CARD
      ENDIF

      NPCARD = NPCARD + 1



      RETURN

! **********************************************************************************************************************************
 1128 FORMAT(' *ERROR  1128: ON ',A,A,' THE IDs MUST BE IN INCREASING ORDER FOR THRU OPTION')

 1152 FORMAT(' *ERROR  1152: ON ',A,A,' ELEM IDs MUST BE > 0')

! --- cbeam_add begin --- !
 1199 FORMAT(' *ERROR      : ON ',A,A,' PLOAD1 TYPE "',A,'" IS NOT SUPPORTED. USE FYE/FZE/FXE/MXE/MYE/MZE',                 &
                    ' (legacy FY/FZ and Y/Z also accepted)')

 1200 FORMAT(' *ERROR      : ON ',A,A,' PLOAD1 SCALE "',A,'" IS NOT SUPPORTED. USE FR, LE, FRPR OR LEPR')
! --- cbeam_add end --- !

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PLOAD2
