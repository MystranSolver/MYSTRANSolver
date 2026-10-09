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

   MODULE DOF_SETS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BD_ASET, BD_ASET1, BD_USET, BD_USET1, BD_SUPORT

   CONTAINS

      SUBROUTINE BD_ASET ( CARD )

! Processes ASET, OMIT Bulk Data Cards
! Data is written to file LINK1N for later processing after checks on format of data.
! Each record contains:

!         COMPJ, GRIDJ, GRIDJ, SET

! GRIDJ is written twice to be compatible with the data written to file LINK1N for ASET1 data

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1N
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NAOCARD
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_ASET'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=LEN(TSET_CHR_LEN)):: SET               ! 'A' or 'O' depending on whether the B.D card is ASET or OMIT

      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0     ! Count of no. of errors when data fields are read from ASET/OMIT cards
      INTEGER(LONG)                   :: COMPJ     = 0     ! Displ component(s)  read from a B.D. ASET/OMIT card
      INTEGER(LONG)                   :: GRIDJ     = 0     ! A grid point number read from a B.D. ASET/OMIT card




! **********************************************************************************************************************************
! ASET, OMIT Bulk Data Card routine

!   FIELD   ITEM
!   -----   ------------
!   2,4,6,8 Grid ID's
!   3,5,7,9 Displ Components

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Find out if card is ASET or OMIT (must be one, otherwise this
! subroutine would not have been invoked):

      IF (JCARD(1)(1:4) == 'ASET') THEN
         SET = 'A'
      ELSE IF (JCARD(1)(1:4) == 'OMIT') THEN
         SET = 'O'
      ENDIF

! Process 8 fields of CARD (pairs of Grid ID's and DOF's)

      DO J=1,4
         IF ((JCARD(2*J)(1:) == ' ') .AND. (JCARD(2*J+1)(1:) == ' ')) THEN
            CYCLE
         ENDIF

! Get Grid ID. If Grid ID field blank, error

         IF (JCARD(2*J)(1:) == ' ') THEN
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1125) 'GRID POINT', JF(2*J), JCARD(1)
            WRITE(F06,1125) 'GRID POINT', JF(2*J), JCARD(1)
         ELSE
            CALL I4FLD ( JCARD(2*J), JF(2*J), GRIDJ )
         ENDIF

! Get DOF components (and put into integer COMPJ).

         CALL IP6CHK ( JCARD(2*J+1), JCARDO, IP6TYP, IDUM )
         IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'ZERO    ') .OR. (IP6TYP == 'BLANK   ')) THEN
            CALL I4FLD ( JCARDO, JF(2*J+1), COMPJ )
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1123) JF(2*J+1),JCARD(1),JF(2*J+1),JCARD(2*J+1)
            WRITE(F06,1123) JF(2*J+1),JCARD(1),JF(2*J+1),JCARD(2*J+1)
         ENDIF

! Write data to file LINK1N if there were no errors. Note, GRIDJ is written twice. This is done since
! ASET1 entries (processed in another subr), can have a format of GRID1 THRU GRID2.

         IF ((JERR == 0) .AND. (IERRFL(2*J) == 'N') .AND. (IERRFL(2*J+1) == 'N')) THEN
            WRITE(L1N) COMPJ,GRIDJ,GRIDJ,SET
            NAOCARD = NAOCARD + 1
         ENDIF

      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,0,4,0,6,0,8,0 )     ! Make sure that there are no imbedded blanks in fields 2, 4, 6, 8
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1123 FORMAT(' *ERROR  1123: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' CARD. MUST BE A COMBINATION OF DIGITS 1-6'                &
                    ,/,14X,' HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1125 FORMAT(' *ERROR  1125: NO ',A,' SPECIFIED IN FIELD',I4,' ON ',A,' CARD')

! **********************************************************************************************************************************

      END SUBROUTINE BD_ASET


      SUBROUTINE BD_ASET1 ( CARD, LARGE_FLD_INP )

! Processes ASET1, OMIT1 Bulk Data Cards
! Data is written to file LINK1N for later processing after checks on format of data.
! Each record contains:   COMPJ, GRIDJ1, GRIDJ2, SET

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1N
      USE SCONTR, ONLY                :  FATAL_ERR, IERRFL, JCARD_LEN, JF, NAOCARD, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_ASET1'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of 8 characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER(LEN(TSET_CHR_LEN))    :: SET               ! 'A' or 'O' depending on whether the B.D card is ASET or OMIT
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! An output from subr TOKCHK called herein

      INTEGER(LONG)                   :: ICONT             ! Indicator of whether a continuation card exists for this parent card
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: IERR              ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0     ! Error indicator for several types of error in format #2 of input
      INTEGER(LONG)                   :: COMPJ     = 0     ! Displ component(s)  read from a B.D. ASET/OMIT card
      INTEGER(LONG)                   :: GRIDJ     = 0     ! A grid point number read from a B.D. ASET/OMIT card in format #2
      INTEGER(LONG)                   :: GRIDJ1    = 0     ! 1st grid in format #1 of ASET/OMIT input
      INTEGER(LONG)                   :: GRIDJ2    = 0     ! 2nd grid in format #1 of ASET/OMIT input




! **********************************************************************************************************************************
! ASET1, OMIT1 Bulk Data Card routine

!   FIELD   ITEM
!   -----   ------------
! Format #1:
!    2      COMPJ, Displ component(s)
!   3-9     GRIDJ's, Grid ID's
! on optional continuation cards:
!   2-9     Grid ID's

! Format #2:
!    2      COMPJ , Displ component(s)
!    3      GRIDJ1, Grid ID number 1
!    4      "THRU"
!    5      GRIDJ2, Grid ID number 2

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Find out if card is ASET1 or OMIT1:

      IF      (JCARD(1)(1:5) == 'ASET1') THEN
         SET = 'A '
      ELSE IF (JCARD(1)(1:5) == 'OMIT1') THEN
         SET = 'O '
      ENDIF

! Field 4 of ASET1 or OMIT1 must have "THRU" or a grid pt number or blank.

      TOKEN = JCARD(4)(1:8)                                ! Only send the 1st 8 chars of this JCARD. It has been left justified
      CALL TOKCHK ( TOKEN, TOKTYP )                        ! TOKTYP must be THRU', 'INTEGR', or 'BLANK'

! **********************************************************************************************************************************
! Format # 2

      IF (TOKTYP == 'THRU    ') THEN

         JERR = 0

         CALL IP6CHK ( JCARD(2), JCARDO, IP6TYP, IDUM )    ! Get components (and put into integer COMPJ)
         IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'ZERO    ') .OR. (IP6TYP == 'BLANK   ')) THEN
            CALL I4FLD ( JCARDO, JF(2), COMPJ )
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1123) JF(2),JCARD(1),JF(2),JCARD(2)
            WRITE(F06,1123) JF(2),JCARD(1),JF(2),JCARD(2)
         ENDIF

         IF (JCARD(3)(1:) /= ' ') THEN                     ! Get 1st Grid ID, GRIDJ1
            CALL I4FLD ( JCARD(3), JF(3), GRIDJ1 )
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1125) 'GRID POINT', JF(3), JCARD(1)
            WRITE(F06,1125) 'GRID POINT', JF(3), JCARD(1)
         ENDIF

         IF (JCARD(5)(1:) /= ' ') THEN                     ! Get 2nd Grid ID, GRIDJ2
            CALL I4FLD ( JCARD(5), JF(5), GRIDJ2 )
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1125) 'GRID POINT', JF(5), JCARD(1)
            WRITE(F06,1125) 'GRID POINT', JF(5), JCARD(1)
         ENDIF

         IF ((IERRFL(3)=='N') .AND. (IERRFL(5)=='N')) THEN ! Check GRIDJ2 > GRIDJ1 if there were no errors reading them
            IF (GRIDJ2 <= GRIDJ1) THEN
               JERR      = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1128) JCARD(1)
               WRITE(F06,1128) JCARD(1)
            ENDIF
         ENDIF

         CALL BD_IMBEDDED_BLANK ( JCARD,0,3,0,5,0,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 3, 5
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )! Issue warning if fields 6, 7, 8, 9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

         IF ((JERR == 0) .AND. (IERRFL(2) == 'N') .AND. (IERRFL(3) == 'N') .AND. (IERRFL(5) == 'N')) THEN
            WRITE(L1N) COMPJ,GRIDJ1,GRIDJ2,SET             ! Write data to file LINK1N if no errors
            NAOCARD = NAOCARD + 1
         ENDIF


! **********************************************************************************************************************************
! Format #1

      ELSE IF ((TOKTYP == 'INTEGER ') .OR. (TOKTYP == 'BLANK   ')) THEN

         JERR = 0

         CALL IP6CHK ( JCARD(2), JCARDO, IP6TYP, IDUM )    ! Get components (and put into integer COMPJ)
         IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'ZERO    ') .OR. (IP6TYP == 'BLANK   ')) THEN
            CALL I4FLD ( JCARDO, JF(2), COMPJ )
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1123) JF(2),JCARD(1),JF(2),JCARD(2)
            WRITE(F06,1123) JF(2),JCARD(1),JF(2),JCARD(2)
         ENDIF

         DO J=3,9                                          ! Get Grid ID's in fields 3 - 9 and write data to LINK1N
            IF (JCARD(J)(1:) == ' ') THEN
               CYCLE
            ELSE
               CALL I4FLD ( JCARD(J), JF(J), GRIDJ )
               IF ((JERR == 0) .AND. (IERRFL(J) == 'N')) THEN
                  WRITE(L1N) COMPJ,GRIDJ,GRIDJ,SET         ! Note, GRIDJ is written twice to be compatible w/ Format #2
                  NAOCARD = NAOCARD + 1
               ENDIF
            ENDIF
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,0,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 3-9
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

         DO                                                ! Optional continuation cards w/ grid ID's, or blank, in fields 2-9
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
                     CYCLE                                 ! CYCLE to next field when a field is blank
                  ELSE
                     CALL I4FLD ( JCARD(J), JF(J), GRIDJ )
                     IF ((JERR == 0) .AND. (IERRFL(J) == 'N')) THEN
                        WRITE(L1N) COMPJ,GRIDJ,GRIDJ,SET   ! Note we don't write if error in either COMPJ or GRIDJ
                        NAOCARD = NAOCARD + 1
                     ENDIF
                  ENDIF
               ENDDO
               CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 ) ! Make sure that there are no imbedded blanks in fields 2-9
               CALL CRDERR ( CARD )                        ! CRDERR prints errors found when reading fields
               CYCLE
            ELSE
               EXIT
            ENDIF
         ENDDO

      ELSE                                                 ! Error - Field 4 did not have "THRU", or an integer, or was blank

         FATAL_ERR = FATAL_ERR+1
         WRITE(ERR,1127) JF(4), JCARD(1)
         WRITE(F06,1127) JF(4), JCARD(1)
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1123 FORMAT(' *ERROR  1123: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' CARD. MUST BE A COMBINATION OF DIGITS 1-6'                &
                    ,/,14X,' HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1125 FORMAT(' *ERROR  1125: NO ',A,' SPECIFIED IN FIELD',I4,' ON ',A,' CARD')

 1127 FORMAT(' *ERROR  1127: INVALID DATA IN FIELD ',I2,' OF ',A,' CARD. FIELD MUST HAVE THRU OR A GRID NUMBER OR BE BLANK')


 1128 FORMAT(' *ERROR  1128: ON ',A,' THE IDs MUST BE IN INCREASING ORDER FOR THRU OPTION')

! **********************************************************************************************************************************

      END SUBROUTINE BD_ASET1


      SUBROUTINE BD_USET ( CARD )

! Processes USET Bulk Data Cards. Reads and checks data and then write a record to file LINK1X for later processing.
! Each record in file LINK1X has:

!          USET_NAME, COMPJ, GRIDJ, DOFSET

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1X
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, NUM_USET_RECORDS, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  SUPWARN
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, IP6CHK, LEFT_ADJ_BDFLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_USET'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: CHRFLD            ! A character field from the entry
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER( 2*BYTE)              :: USET_NAME         ! Name in field 2 of the USET entry

      INTEGER(LONG)                   :: COMPJ     = 0     ! Displ components constrained at GRIDJ
      INTEGER(LONG)                   :: GRIDJ     = 0     ! Grid ID on USET card
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP6CHK not used herein
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count




! **********************************************************************************************************************************
!  USET Bulk Data Card routine

!    FIELD   ITEM           ARRAY ELEMENT
!    -----   ------------   -------------
!     2      USET name (i.e. "U1", "U2")
!     3      Grid ID
!     4      Comp. numbers
!     5      Grid ID
!     6      Comp.numbers
!     7      Grid ID
!     8      Comp.numbers


! Initialize

      JERR = 0

!  Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Read USET name

      CALL CHAR_FLD ( JCARD(2), JF(2), CHRFLD )            ! Read field 2: USET name
      IF (IERRFL(2) == 'N') THEN
         CALL LEFT_ADJ_BDFLD ( CHRFLD )
         USET_NAME = CHRFLD(1:2)
         IF ((USET_NAME /= 'U1      ') .AND. (USET_NAME /= 'U2      ')) THEN
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1199) SUBR_NAME, USET_NAME
            WRITE(F06,1199) SUBR_NAME, USET_NAME
         ENDIF
      ENDIF

! Process data in fields 3-8

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
               NUM_USET_RECORDS = NUM_USET_RECORDS + 1     ! Incr count of number of entries written to file LINK1X

               WRITE(L1X) USET_NAME(1:2), COMPJ, GRIDJ, GRIDJ
            ENDIF

         ELSE                                              ! Field 3, 5 or 7 is blank

            IF (JCARD(I+1)(1:) /= ' ') THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1144) JCARD(1), JCARD(2), JF(I+1), JF(I)
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1144) JCARD(1), JCARD(2), JF(I+1), JF(I)
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

 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1144 FORMAT(' *WARNING    : ON ',A,' SNAME = ',A,' FIELD ',I3,' IS IGNORED SINCE GRID FIELD ',I3,' IS BLANK.')

 1199 FORMAT(' *ERROR  1199: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' SNAME ON THE ABOVE USET ENTRY MUST BE "U1" OR "U2" BUT IS "',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE BD_USET


      SUBROUTINE BD_USET1 ( CARD, LARGE_FLD_INP )

! Processes USET1 Bulk Data Cards. Reads and checks data and then write a record to file LINK1X for later processing.
! Each record in file LINK1X has:

!          USET_NAME, COMPJ, GRIDJ, DOFSET

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1X
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NUM_USET_RECORDS
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, IP6CHK, LEFT_ADJ_BDFLD
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_USET1'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN=JCARD_LEN)        :: CHRFLD            ! A character field from the entry
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! An output from subr TOKCHK called herein
      CHARACTER( 2*BYTE)              :: USET_NAME         ! Name in field 2 of the USET entry

      INTEGER(LONG)                   :: COMPJ     = 0     ! DOF's constrained at GRIDJ
      INTEGER(LONG)                   :: GRIDJ1    = 0     ! Grid ID on USET1 card
      INTEGER(LONG)                   :: GRIDJ2    = 0     ! Grid ID on USET1 card
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count




! **********************************************************************************************************************************
!  USET1 Bulk Data Card routine

!    FIELD   ITEM           ARRAY ELEMENT
!    -----   ------------   -------------
! Format #1:
!     2      USET name (e.g. "U1", "U2", "R", etc)
!     3      Comp. numbers
!     4-9    Grid ID's
! Optional continuation cards
!     2-9    Grid ID's

! Format #2:
!     2      USET name (e.g. "U1", "U2", "R", etc)
!     3      Component numbers
!     4      Grid ID number 1
!     5      "THRU"
!     6      Grid ID number 2


! Initialize

      JERR = 0

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Read USET name

      CALL CHAR_FLD ( JCARD(2), JF(2), CHRFLD )            ! Read field 2: USET name
      IF (IERRFL(2) == 'N') THEN
         CALL LEFT_ADJ_BDFLD ( CHRFLD )
         USET_NAME = CHRFLD(1:2)
         IF ((USET_NAME /= 'U1      ') .AND. (USET_NAME /= 'U2      ')) THEN
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1199) SUBR_NAME, USET_NAME
            WRITE(F06,1199) SUBR_NAME, USET_NAME
         ENDIF
      ENDIF

! Get component numbers

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

! Field 5 of USET1 must have "THRU" or a grid pt number or blank.

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
            NUM_USET_RECORDS = NUM_USET_RECORDS + 1        ! Incr count of number of entries written to file LINK1X
            WRITE(L1X) USET_NAME, COMPJ, GRIDJ1, GRIDJ2
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
                  NUM_USET_RECORDS = NUM_USET_RECORDS + 1  ! Incr count of number of entries written to file LINK1X
                  WRITE(L1X) USET_NAME, COMPJ, GRIDJ1, GRIDJ1
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
                           NUM_USET_RECORDS = NUM_USET_RECORDS + 1
                           WRITE(L1X) USET_NAME, COMPJ, GRIDJ1, GRIDJ1
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

         IF (NUM_USET_RECORDS == 0) THEN                   ! No grids were specified on the USET1 entry, so error
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
 1100 FORMAT(' *ERROR  1100: ',A,' ENTRY MUST HAVE AT LEAST 1 GRID DEFINED (IN FIELDS 4-9 OF PARENT OR 2-9 OF CONTINUATION)')

 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1125 FORMAT(' *ERROR  1125: NO ',A,' SPECIFIED IN FIELD',I4,' ON ',A,' CARD')

 1127 FORMAT(' *ERROR  1127: INVALID DATA IN FIELD ',I2,' OF ',A,' CARD. FIELD MUST HAVE THRU OR A GRID NUMBER OR BE BLANK')

 1128 FORMAT(' *ERROR  1128: ON ',A,A,' THE IDs MUST BE IN INCREASING ORDER FOR THRU OPTION')

 1199 FORMAT(' *ERROR  1199: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' SNAME ON THE ABOVE USET ENTRY MUST BE "U1" OR "U2" BUT IS "',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE BD_USET1


      SUBROUTINE BD_SUPORT ( CARD )

! Processes SUPORT Bulk Data Cards
! Data is written to file LINK1T for later processing after checks on format of data.
! Each record contains:   GRIDJ1, COMPJ

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1T
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NUM_SUPT_CARDS
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_SUPORT'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein

      INTEGER(LONG)                   :: COMP(4)   = 0     ! Displ component(s)  read from a B.D. SUPORT card
      INTEGER(LONG)                   :: GRID(4)   = 0     ! A grid point number read from a B.D. SUPORT card
      INTEGER(LONG)                   :: IERR              ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: JERR      = 0     ! Error indicator for several types of error
      INTEGER(LONG)                   :: NUM_PAIRS         ! Number of pairs of grid/comp found on this SUPORT card




! **********************************************************************************************************************************
! SUPORT Bulk Data Card routine

!   FIELD   ITEM
!   -----   ------------
!    2      GRID, Grid ID number 1
!    3      COMP, Displ component(s) for GRID1
!    4      GRID, Grid ID number 2
!    5      COMP, Displ component(s) for GRID2
!    6      GRID, Grid ID number 3
!    7      COMP, Displ component(s) for GRID3
!    8      GRID, Grid ID number 4
!    9      COMP, Displ component(s) for GRID4

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! **********************************************************************************************************************************
      JERR = 0

      NUM_PAIRS = 0
      DO I=2,8,2

         IERR = 0
                                                           ! First make sure that if grid or comp is specified, other one is also
         IF ((JCARD(I)(1:) == ' ') .AND. (JCARD(I+1)(1:) /= ' ')) THEN
            WRITE(ERR,1125) 'GRID POINT', JF(I), JCARD(1)
            WRITE(F06,1125) 'GRID POINT', JF(I), JCARD(1)
            IERR = IERR + 1
         ENDIF
         IF ((JCARD(I)(1:) /= ' ') .AND. (JCARD(I+1)(1:) == ' ')) THEN
            WRITE(ERR,1125) 'DISPL COMPONENT', JF(I+1), JCARD(1)
            WRITE(F06,1125) 'DISPL COMPONENT', JF(I+1), JCARD(1)
            IERR = IERR + 1
         ENDIF

         IF (IERR == 0) THEN
            IF (JCARD(I)(1:) /= ' ') THEN                  ! Grid field not blank so read GRID, COMP
               NUM_PAIRS = NUM_PAIRS + 1
               CALL I4FLD ( JCARD(I), JF(I), GRID(NUM_PAIRS) )
               CALL IP6CHK ( JCARD(I+1), JCARDO, IP6TYP, IDUM )
               IF (IP6TYP == 'COMP NOS') THEN
                  CALL I4FLD ( JCARDO, JF(I+1), COMP(NUM_PAIRS) )
               ELSE
                  JERR      = JERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1123) JF(I+1),JCARD(1),JF(I+1),JCARD(I+1)
                  WRITE(F06,1123) JF(I+1),JCARD(1),JF(I+1),JCARD(I+1)
               ENDIF
            ENDIF
         ENDIF

      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,0,4,0,6,0,8,0 )  ! Make sure that there are no imbedded blanks in fields 2,4,6,8
      CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      IF (JERR == 0) THEN
         DO I=1,NUM_PAIRS
            WRITE(L1T) GRID(I), COMP(I)                 ! Write data to file LINK1T if no errors
            NUM_SUPT_CARDS = NUM_SUPT_CARDS + 1
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1123 FORMAT(' *ERROR  1123: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' CARD. MUST BE A COMBINATION OF DIGITS 1-6'                &
                    ,/,14X,' HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1125 FORMAT(' *ERROR  1125: NO ',A,' SPECIFIED IN FIELD',I4,' ON ',A,' CARD')

! **********************************************************************************************************************************

      END SUBROUTINE BD_SUPORT

   END MODULE DOF_SETS
