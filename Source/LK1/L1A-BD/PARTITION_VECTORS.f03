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

   MODULE PARTITION_VECTORS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BD_PARVEC, BD_PARVEC1

   CONTAINS

      SUBROUTINE BD_PARVEC ( CARD )

! Processes PARVEC Bulk Data Cards. Reads and checks data and then writes a record to file LINK1V for later processing.
! Each record in file LINK1V has:

!          PARTVEC_NAME, COMPJ, GRIDJ, DOFSET

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1V
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, NUM_PARTVEC_RECORDS, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE OUTPUT4_MATRICES, ONLY      :  ACT_OU4_MYSTRAN_NAMES
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  SUPWARN
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, IP6CHK, LEFT_ADJ_BDFLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PARVEC'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: CHRFLD            ! A character field from the entry
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein

      CHARACTER(LEN(ACT_OU4_MYSTRAN_NAMES))                                                                                        &
                                      :: PARTVEC_NAME      ! Name in field 2 of the PARTVEC entry

      INTEGER(LONG)                   :: COMPJ     = 0     ! Displ components constrained at GRIDJ
      INTEGER(LONG)                   :: GRIDJ     = 0     ! Grid ID on PARTVEC card
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP6CHK not used herein
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count




! **********************************************************************************************************************************
!  PARTVEC Bulk Data Card routine

!    FIELD   ITEM           ARRAY ELEMENT
!    -----   ------------   -------------
!     2      PARTVEC name (i.e. "PHIG", "PHIGZ")
!     3      Grid ID
!     4      Comp. numbers
!     5      Grid ID
!     6      Comp.numbers
!     7      Grid ID
!     8      Comp.numbers


!  Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Read PARTVEC name

      CALL CHAR_FLD ( JCARD(2), JF(2), CHRFLD )            ! Read field 2: PARTVEC name
      IF (IERRFL(2) == 'N') THEN
         CALL LEFT_ADJ_BDFLD ( CHRFLD )
         PARTVEC_NAME = CHRFLD(1:)
      ENDIF

! If card has no data, write warning and return

      IF ((JCARD(3)(1:) == ' ') .AND. (JCARD(4)(1:) == ' ') .AND. (JCARD(5)(1:) == ' ') .AND.                                      &
          (JCARD(6)(1:) == ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,201) JCARD(1), JCARD(2), '3-8'
         WRITE(F06,201) JCARD(1), JCARD(2), '3-8'
      ENDIF

! Process data in fields 3-8

      JERR = 0
      DO I=3,7,2                                           ! There can be 3 sets of (Grid, Comp) on each card

         IF (JCARD(I)(1:) /= ' ') THEN

            CALL I4FLD ( JCARD(I), JF(I), GRIDJ )          ! Read Grid ID
                                                           ! Read displ components
            CALL IP6CHK ( JCARD(I+1), JCARDO, IP6TYP, IDUM )
            IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'ZERO    ') .OR. (IP6TYP == 'BLANK   ')) THEN
               CALL I4FLD ( JCARDO, JF(I+1), COMPJ )
            ELSE
               JERR      = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1124) JF(I+1), JCARD(1), JCARD(2), JF(I+1), JCARD(I+1)
               WRITE(F06,1124) JF(I+1), JCARD(1), JCARD(2), JF(I+1), JCARD(I+1)
            ENDIF

            IF ((JERR == 0 ) .AND. (IERRFL(I) == 'N') .AND. (IERRFL(I+1) == 'N')) THEN
               NUM_PARTVEC_RECORDS = NUM_PARTVEC_RECORDS + 1     ! Incr count of number of entries written to file LINK1V
               WRITE(L1V) PARTVEC_NAME, COMPJ, GRIDJ, GRIDJ
            ENDIF

         ELSE                                              ! Field 3, 5 or 7 is blank

            IF (JCARD(I+1)(1:) /= ' ') THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,202) JCARD(1), JCARD(2), JF(I+1), JF(I)
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,202) JCARD(1), JCARD(2), JF(I+1), JF(I)
               ENDIF
            ENDIF

         ENDIF

      ENDDO

      CALL BD_IMBEDDED_BLANK   ( JCARD,0,3,0,5,0,7,0,0 )   ! Make sure that there are no imbedded blanks in fields 2,3,5,6,8,9
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

  201 FORMAT(' *WARNING    : ON ',A,' SNAME = ',A,' FIELDS ',A,' ARE NULL. ENTRY IGNORED')

  202 FORMAT(' *WARNING    : ON ',A,' SNAME = ',A,' FIELD ',I3,' IS IGNORED SINCE GRID FIELD ',I3,' IS BLANK.')

 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE BD_PARVEC


      SUBROUTINE BD_PARVEC1 ( CARD, LARGE_FLD_INP )

! Processes PARTVEC1 Bulk Data Cards. Reads and checks data and then write a record to file LINK1V for later processing.
! Each record in file LINK1V has:

!          PARTVEC1_NAME, COMPJ, GRIDJ, DOFSET

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1V
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NUM_PARTVEC_RECORDS, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE OUTPUT4_MATRICES, ONLY      :  ACT_OU4_MYSTRAN_NAMES
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, IP6CHK, LEFT_ADJ_BDFLD
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PARVEC1'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN=JCARD_LEN)        :: CHRFLD            ! A character field from the entry
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! An output from subr TOKCHK called herein

      CHARACTER(LEN(ACT_OU4_MYSTRAN_NAMES))                                                                                        &
                                      :: PARTVEC1_NAME     ! Name in field 2 of the PARTVEC1 entry


      INTEGER(LONG)                   :: COMPJ     = 0     ! DOF's constrained at GRIDJ
      INTEGER(LONG)                   :: GRIDJ1    = 0     ! Grid ID on PARTVEC1 card
      INTEGER(LONG)                   :: GRIDJ2    = 0     ! Grid ID on PARTVEC1 card
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count




! **********************************************************************************************************************************
!  PARTVEC1 Bulk Data Card routine

!    FIELD   ITEM           ARRAY ELEMENT
!    -----   ------------   -------------
! Format #1:
!     2      PARTVEC1 name (e.g. "PHIG", "PHIGZ", etc)
!     3      Comp. numbers
!     4-9    Grid ID's
! Optional continuation cards
!     2-9    Grid ID's

! Format #2:
!     2      PARTVEC1 name (e.g. "U1", "U2", "R", etc)
!     3      Component numbers
!     4      Grid ID number 1
!     5      "THRU"
!     6      Grid ID number 2


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! If card has no data, write warning and return

      IF ((JCARD(3)(1:) == ' ') .AND. (JCARD(4)(1:) == ' ') .AND. (JCARD(5)(1:) == ' ') .AND. (JCARD(6)(1:) == ' ') .AND.          &
          (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ') .AND. (JCARD(9)(1:) == ' ')) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,201) JCARD(1), JCARD(2), '3-8'
         WRITE(F06,201) JCARD(1), JCARD(2), '3-8'
      ENDIF

! Read PARTVEC1 name

      CALL CHAR_FLD ( JCARD(2), JF(2), CHRFLD )            ! Read field 2: PARTVEC1 name
      IF (IERRFL(2) == 'N') THEN
         CALL LEFT_ADJ_BDFLD ( CHRFLD )
         PARTVEC1_NAME = CHRFLD(1:)
      ENDIF

! Get component numbers

      JERR = 0
      CALL IP6CHK ( JCARD(3), JCARDO, IP6TYP, IDUM )       ! Read field 3: components
      IF (IERRFL(3) == 'N') THEN
         IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'ZERO    ') .OR. (IP6TYP == 'BLANK   ')) THEN
            CALL I4FLD ( JCARDO, JF(3), COMPJ )
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1124) JF(3), JCARD(1), JCARD(2), JF(3), JCARD(3)
            WRITE(F06,1124) JF(3), JCARD(1), JCARD(2), JF(3), JCARD(3)
         ENDIF
      ELSE
         JERR = JERR + 1
      ENDIF

! Field 5 of PARTVEC1 must have "THRU" or a grid pt number or blank.

      TOKEN = JCARD(5)(1:8)                                ! Only send the 1st 8 chars of this JCARD. It has been left justified
      CALL TOKCHK ( TOKEN, TOKTYP )                        ! TOKTYP must be THRU', 'INTEGER', or 'BLANK'

! **********************************************************************************************************************************
! Format # 2

      IF (TOKTYP == 'THRU    ') THEN

         JERR = 0

         IF (JCARD(4)(1:) /= ' ') THEN                     ! Get 1st Grid ID, GRIDJ1
            CALL I4FLD ( JCARD(4), JF(4), GRIDJ1 )
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1125) 'GRID POINT', JF(4), JCARD(1)
            WRITE(F06,1125) 'GRID POINT', JF(4), JCARD(1)
         ENDIF

         IF (JCARD(6)(1:) /= ' ') THEN                     ! Get 2nd Grid ID, GRIDJ2
            CALL I4FLD ( JCARD(6), JF(6), GRIDJ2 )
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1125) 'GRID POINT', JF(6), JCARD(1)
            WRITE(F06,1125) 'GRID POINT', JF(6), JCARD(1)
         ENDIF

         IF ((IERRFL(4)=='N') .AND. (IERRFL(5)=='N')) THEN ! Check GRIDJ2 > GRIDJ1 if there were no errors reading them
            IF (GRIDJ2 < GRIDJ1) THEN
               JERR      = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1128) JCARD(1), JCARD(2)
               WRITE(F06,1128) JCARD(1), JCARD(2)
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK ( JCARD,2,0,4,5,6,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2,4,5,6
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,7,8,9 )! Issue warning if fields 7,8,9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

         IF (JERR == 0) THEN
            NUM_PARTVEC_RECORDS = NUM_PARTVEC_RECORDS + 1  ! Incr count of number of entries written to file LINK1V
            WRITE(L1V) PARTVEC1_NAME, COMPJ, GRIDJ1, GRIDJ2
         ENDIF

! **********************************************************************************************************************************
! Format # 1

      ELSE IF ((TOKTYP == 'INTEGER ') .OR. (TOKTYP == 'BLANK   ')) THEN

! Read and check data on parent card

         DO J=4,9                                          ! Read fields 4-9: Grid ID's.
            IF (JCARD(J)(1:) == ' ') THEN
               CYCLE
            ELSE
               CALL I4FLD ( JCARD(J), JF(J), GRIDJ1 )      ! Read another grid ID
               IF ((JERR == 0) .AND. (IERRFL(J) == 'N')) THEN
                  NUM_PARTVEC_RECORDS = NUM_PARTVEC_RECORDS + 1  ! Incr count of number of entries written to file LINK1V
                  WRITE(L1V) PARTVEC1_NAME, COMPJ, GRIDJ1, GRIDJ1
               ELSE
                  JERR = JERR + 1
               ENDIF
            ENDIF
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,2,0,4,5,6,7,8,9 )  ! Make sure there are no imbeded blanks fields 2-9, except 3 (components)
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

! Read and check data on optional continuation cards if JERR = 0 to this point

         IF (JERR == 0) THEN

            DO
               IF (LARGE_FLD_INP == 'N') THEN
                  CALL NEXTC  ( CARD, ICONT, IERR )
               ELSE
                  CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
                  CARD = CHILD
               ENDIF
               CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
               IF (ICONT == 1) THEN
                  DO J=2,9
                     IF (JCARD(J)(1:) == ' ') THEN
                        CYCLE
                     ELSE                                  ! Read another grid ID
                        CALL I4FLD ( JCARD(J), JF(J), GRIDJ1)
                        IF ((JERR == 0) .AND. (IERRFL(J) == 'N')) THEN
                           NUM_PARTVEC_RECORDS = NUM_PARTVEC_RECORDS + 1
                           WRITE(L1V) PARTVEC1_NAME, COMPJ, GRIDJ1, GRIDJ1
                        ENDIF
                     ENDIF
                  ENDDO
                                                           ! Make sure that there are no imbedded blanks in fields 2-9
                  CALL BD_IMBEDDED_BLANK (JCARD,2,3,4,5,6,7,8,9 )
                  CALL CRDERR ( CARD )                     ! CRDERR prints errors found when reading fields

                  CYCLE
               ELSE
                  EXIT
               ENDIF

            ENDDO

         ENDIF

         IF (NUM_PARTVEC_RECORDS == 0) THEN               ! No grids were specified on the PARTVEC1 entry, so error
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1100) JCARD(1)
            WRITE(F06,1100) JCARD(1)
         ENDIF

      ELSE                                                 ! Error - Field 5 did not have "THRU", or an integer, or was blank

         FATAL_ERR = FATAL_ERR+1
         WRITE(ERR,1127) JF(5), JCARD(1)
         WRITE(F06,1127) JF(5), JCARD(1)
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ENDIF



      RETURN

! **********************************************************************************************************************************
  201 FORMAT(' *WARNING    : ON ',A,' SNAME = ',A,' FIELDS ',A,' ARE NULL. ENTRY IGNORED')

 1100 FORMAT(' *ERROR  1100: ',A,' ENTRY MUST HAVE AT LEAST 1 GRID DEFINED (IN FIELDS 4-9 OF PARENT OR 2-9 OF CONTINUATION)')

 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1125 FORMAT(' *ERROR  1125: NO ',A,' SPECIFIED IN FIELD',I4,' ON ',A,' CARD')

 1127 FORMAT(' *ERROR  1127: INVALID DATA IN FIELD ',I2,' OF ',A,' CARD. FIELD MUST HAVE THRU OR A GRID NUMBER OR BE BLANK')

 1128 FORMAT(' *ERROR  1128: ON ',A,A,' THE IDs MUST BE IN INCREASING ORDER FOR THRU OPTION')

! **********************************************************************************************************************************

      END SUBROUTINE BD_PARVEC1

   END MODULE PARTITION_VECTORS
