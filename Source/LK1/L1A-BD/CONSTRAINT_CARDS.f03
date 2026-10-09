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

   MODULE CONSTRAINT_CARDS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BD_SPC, BD_SPC1, BD_SPCADD, BD_SPCADD0, BD_MPC, BD_MPC0, BD_MPCADD, BD_MPCADD0

   CONTAINS

      SUBROUTINE BD_SPC ( CARD, CC_SPC_FND )

! Processes SPC Bulk Data Cards. Reads and checks data and then write a record to file LINK1O for later processing.
! Each record in file LINK1O has:

!          SETID, COMPJ, GRIDJ, RSPCJ, DOFSET

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1O
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, LSPC, NSPC, NUM_SPC_RECORDS, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL, SUPWARN
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN
      USE MODEL_STUF, ONLY            :  SPC_SIDS, SPCSET

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_SPC'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD              ! A Bulk Data card
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_SPC_FND        ! ='Y' if this SPC is a set requested in Case Control
      CHARACTER(LEN=LEN(TSET_CHR_LEN)):: DOFSET            ! 'SB' or 'SE' to denote the set that an SPC'd DOF belongs to
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein

      INTEGER(LONG)                   :: COMPJ     = 0     ! Displ components constrained at GRIDJ
      INTEGER(LONG)                   :: GRIDJ     = 0     ! Grid ID on SPC card
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      INTEGER(LONG)                   :: SETID     = 0     ! SPC set ID


      REAL(DOUBLE)                    :: DEPS1             ! A small positive number to compare real zero
      REAL(DOUBLE)                    :: RSPCJ     = ZERO  ! Enforced displ value



! **********************************************************************************************************************************
!  SPC Bulk Data Card routine

!    FIELD   ITEM           ARRAY ELEMENT
!    -----   ------------   -------------
!     2      Set ID
!     3      Grid ID
!     4      Comp. numbers
!     5      Displacement
!     6      Grid ID
!     7      Comp.numbers
!     8      Displacement


!xx   CC_SPC_FND = 'N'                                    ! ERROR. When this is in, then it is reset each time an SPC card is read

      DEPS1 = DABS(EPSIL(1))

!  Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NSPC = NSPC + 1
!xx   IF (NSPC > LSPC) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LSPC
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LSPC
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read set ID

      CALL I4FLD ( JCARD(2), JF(2), SETID )                ! Read field 2: Set ID
      IF (IERRFL(2) == 'N') THEN
         IF (SETID == SPCSET) THEN
            CC_SPC_FND = 'Y'
         ENDIF
         SPC_SIDS(NSPC) = SETID
      ELSE
         JERR = JERR + 1
      ENDIF

      DO I=1,2                                             ! There can be 2 sets of (Grid, Comp, Value) on each card

         IF (JCARD(3*I)(1:) /= ' ') THEN

            JERR = 0
            CALL I4FLD ( JCARD(3*I), JF(3*I), GRIDJ )      ! Read Grid ID
                                                           ! Read displ components
            CALL IP6CHK ( JCARD(3*I+1), JCARDO, IP6TYP, IDUM )
            IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'ZERO    ') .OR. (IP6TYP == 'BLANK   ')) THEN
               CALL I4FLD ( JCARDO, JF(3*I+1), COMPJ )
            ELSE
               JERR      = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1124) JF(3*I+1),JCARD(1),JCARD(2),JF(3*I+1),JCARD(3*I+1)
               WRITE(F06,1124) JF(3*I+1),JCARD(1),JCARD(2),JF(3*I+1),JCARD(3*I+1)
            ENDIF

            CALL R8FLD ( JCARD(3*I+2), JF(3*I+2), RSPCJ )  ! Read permanent SPC

            IF (DABS(RSPCJ) > DEPS1) THEN                  ! SPC are SE set if enforced displ > 0, SB set if = 0
               DOFSET = 'SE'
            ELSE
               DOFSET = 'SB'
            ENDIF

            IF ((JERR == 0 ) .AND. (IERRFL(3*I) == 'N') .AND. (IERRFL(3*I+1) == 'N') .AND. (IERRFL(3*I+2) == 'N')) THEN
               NUM_SPC_RECORDS = NUM_SPC_RECORDS + 1       ! Incr count of number of entries written to file LINK1O
               WRITE(L1O) SETID,COMPJ,GRIDJ,GRIDJ,RSPCJ,DOFSET
            ENDIF

         ELSE                                              ! Field 3 or 6 is blank

            IF ((JCARD(3*I+1)(1:) /= ' ') .OR. (JCARD(3*I+2)(1:) /= ' ')) THEN
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,101) CARD
               WRITE(ERR,1144) JCARD(1),JCARD(2),JF(3*I+1),JF(3*I+2),JF(3*I)
               IF (SUPWARN == 'N') THEN
                  IF (ECHO == 'NONE  ') THEN
                     WRITE(F06,101) CARD
                  ENDIF
                  WRITE(F06,1144) JCARD(1),JCARD(2),JF(3*I+1),JF(3*I+2),JF(3*I)
               ENDIF
            ENDIF

         ENDIF

      ENDDO

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,0,5,6,0,8,0 )   ! Make sure that there are no imbedded blanks in fields 2,3,5,6,8,9
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1144 FORMAT(' *WARNING    : ON ',A,' SET ID = ',A,' FIELDS ',I3,' AND ',I3,' ARE IGNORED SINCE GRID FIELD ',I3,' IS BLANK.')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_SPC


      SUBROUTINE BD_SPC1 ( CARD, LARGE_FLD_INP, CC_SPC_FND )

! Processes SPC Bulk Data Cards. Reads and checks data and then write a record to file LINK1O for later processing.
! Each record in file LINK1O has:

!          SETID, COMPJ, GRIDJ, RSPCJ = 0., DOFSET

! Note that RSPCJ is written to be compatible with the data written for Bulk Data card SPC

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1O
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LSPC1, NSPC1, NUM_SPC1_RECORDS
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN
      USE MODEL_STUF, ONLY            :  SPC1_SIDS, SPCSET

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_SPC1'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_SPC_FND        ! ='Y' if this SPC is a set requested in Case Control
      CHARACTER(LEN=LEN(TSET_CHR_LEN)):: DOFSET            ! 'SB' or 'SE' to denote the set that an SPC'd DOF belongs to
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! An output from subr TOKCHK called herein

      INTEGER(LONG)                   :: COMPJ     = 0     ! DOF's constrained at GRIDJ
      INTEGER(LONG)                   :: GRIDJ1    = 0     ! Grid ID on SPC card
      INTEGER(LONG)                   :: GRIDJ2    = 0     ! Grid ID on SPC card
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      INTEGER(LONG)                   :: SETID     = 0     ! SPC set ID


      REAL(DOUBLE) , PARAMETER        :: RSPCJ     = ZERO  ! Enforced displ value (always zero on SPC1). Included for file LINK1O
!                                                            with SPC format.



! **********************************************************************************************************************************
!  SPC1 Bulk Data Card routine

!    FIELD   ITEM           ARRAY ELEMENT
!    -----   ------------   -------------
! Format #1:
!     2      Set ID
!     3      Comp. numbers
!     4-9    Grid ID's
! Optional continuation cards
!     2-9    Grid ID's

! Format #2:
!     2      Set ID
!     3      Component numbers
!     4      Grid ID number 1
!     5      "THRU"
!     6      Grid ID number 2


! All SPC1 are SB set

      DOFSET = 'SB'

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NSPC1 = NSPC1 + 1
!xx   IF (NSPC1 > LSPC1) THEN                              ! Check for overflow
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LSPC1
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LSPC1
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Check if SPC set ID matches a Case Control request

      CALL I4FLD ( JCARD(2), JF(2), SETID )                ! Read field 2: Set ID
      IF (IERRFL(2) == 'N') THEN
         IF (SETID == SPCSET) THEN
            CC_SPC_FND = 'Y'
         ENDIF
         SPC1_SIDS(NSPC1) = SETID
      ELSE
         JERR = JERR + 1
      ENDIF

! Get component numbers that are SPC'd

      CALL IP6CHK ( JCARD(3), JCARDO, IP6TYP, IDUM )       ! Read field 3: components
      IF (IERRFL(3) == 'N') THEN
         IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'ZERO    ') .OR. (IP6TYP == 'BLANK   ')) THEN
            CALL I4FLD ( JCARDO, JF(3), COMPJ )
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1124) JF(3),JCARD(1),JCARD(2),JF(3),JCARD(3)
            WRITE(F06,1124) JF(3),JCARD(1),JCARD(2),JF(3),JCARD(3)
         ENDIF
      ELSE
         JERR = JERR + 1
      ENDIF

! Field 5 of SPC1 must have "THRU" or a grid pt number or blank.

      TOKEN = JCARD(5)(1:8)                                ! Only send the 1st 8 chars of this JCARD. It has been left justified
      CALL TOKCHK ( TOKEN, TOKTYP )                       ! TOKTYP must be THRU', 'INTEGER', or 'BLANK'

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

         IF ((IERRFL(4)=='N') .AND. (IERRFL(6)=='N')) THEN ! Check GRIDJ2 > GRIDJ1 if there were no errors reading them
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
            NUM_SPC1_RECORDS = NUM_SPC1_RECORDS + 1        ! Incr count of number of entries written to file LINK1O
            WRITE(L1O) SETID,COMPJ,GRIDJ1,GRIDJ2,RSPCJ,DOFSET
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
                  NUM_SPC1_RECORDS = NUM_SPC1_RECORDS + 1  ! Incr count of number of entries written to file LINK1O
                  WRITE(L1O) SETID,COMPJ,GRIDJ1,GRIDJ1,RSPCJ,DOFSET
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
                     ELSE
                        CALL I4FLD ( JCARD(J), JF(J), GRIDJ1)! Read another grid ID
                        IF ((JERR == 0) .AND. (IERRFL(J) == 'N')) THEN
                           NUM_SPC1_RECORDS = NUM_SPC1_RECORDS + 1
                           WRITE(L1O) SETID,COMPJ,GRIDJ1,GRIDJ1,RSPCJ,DOFSET
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

         IF (NUM_SPC1_RECORDS == 0) THEN                   ! No grids were specified on the SPC1 entry, so error
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

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_SPC1


      SUBROUTINE BD_SPCADD ( CARD, LARGE_FLD_INP, CC_SPC_FND )

! Processes SPCADD Bulk Data Cards. Reads and checks data and enters data into array SPCADD_SIDS

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LSPCADDR, LSUB, NSPCADD, LSPCADDC, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  SPCADD_SIDS, SPCSET, SUBLOD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_SPCADD'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_SPC_FND        ! 'Y' if B.D  card w/ same set ID as C.C. SPC = SID
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER( 3*BYTE)              :: NAME1   = 'SPC'   !
      CHARACTER( 6*BYTE)              :: NAME2   = 'SPCADD'!

      INTEGER(LONG)                   :: I,J,K             ! DO loop index
      INTEGER(LONG)                   :: NUM_SETIDS        ! Counter on number of set ID's on SPCADD card
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: SETID             ! Set ID for this SPCADD Bulk Data card




! **********************************************************************************************************************************
! SPCADD Bulk Data Card routine

!   FIELD   ITEM            ARRAY ELEMENT
!   -----   ------------    -------------
!    2      SPCADD set ID    SPCADD(nspcadd, 1)
!    3      SPC/SPC1 set ID  SPCADD(nspcadd, 2)
!    4      SPC/SPC1 set ID  SPCADD(nspcadd, 3)
!    5      SPC/SPC1 set ID  SPCADD(nspcadd, 4)
!    6      SPC/SPC1 set ID  SPCADD(nspcadd, 5)
!    7      SPC/SPC1 set ID  SPCADD(nspcadd, 6)
!    8      SPC/SPC1 set ID  SPCADD(nspcadd, 7)
!    9      SPC/SPC1 set ID  SPCADD(nspcadd, 8)

! 1st continuation card:
!
!    2      SPC/SPC1 set ID  SPCADD(nspcadd, 9)
!    3      SPC/SPC1 set ID  SPCADD(nspcadd,10)
!    4      SPC/SPC1 set ID  SPCADD(nspcadd,11)
!    5      SPC/SPC1 set ID  SPCADD(nspcadd,12)
!    6      SPC/SPC1 set ID  SPCADD(nspcadd,13)
!    7      SPC/SPC1 set ID  SPCADD(nspcadd,14)
!    8      SPC/SPC1 set ID  SPCADD(nspcadd,15)
!    9      SPC/SPC1 set ID  SPCADD(nspcadd,16)

! Subsequent con't cards follow the same patterm as the 1st


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      NSPCADD = NSPCADD+1

! Read and check data on parent card

      CALL I4FLD ( JCARD(2), JF(2), SETID )                ! Read set ID for this SPCADD Bulk Data card
      IF (IERRFL(2) == 'N') THEN
         SPCADD_SIDS(NSPCADD,1) = SETID
         DO I=1,NSUB
            IF (SETID == SPCSET) THEN
               CC_SPC_FND = 'Y'
            ENDIF
         ENDDO
      ENDIF

      NUM_SETIDS = 1                                       ! Read SPCADD ID's on parent card.
      DO J=3,9
         IF (JCARD(J)(1:) == ' ') THEN
            CYCLE                                          ! CYCLE if a set ID field is blank
         ELSE
            NUM_SETIDS = NUM_SETIDS + 1
            IF (NUM_SETIDS > LSPCADDC) THEN
               WRITE(ERR,1139) SETID, NAME1, NAME2, LSPCADDC
               WRITE(F06,1139) SETID, NAME1, NAME2, LSPCADDC
               CALL OUTA_HERE ( 'Y' )                              ! Coding error, so quit
            ENDIF
            CALL I4FLD ( JCARD(J), JF(J),  SPCADD_SIDS(NSPCADD,NUM_SETIDS) )
            DO K=1,NUM_SETIDS-1                            ! Check for duplicate set ID's
               IF (SPCADD_SIDS(NSPCADD,NUM_SETIDS) == SPCADD_SIDS(NSPCADD,K)) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1140) SPCADD_SIDS(NSPCADD,NUM_SETIDS),NAME2,SETID
                  WRITE(F06,1140) SPCADD_SIDS(NSPCADD,NUM_SETIDS),NAME2,SETID
               ENDIF
            ENDDO
            IF (IERRFL(J) == 'N') THEN
               IF (SPCADD_SIDS(NSPCADD,NUM_SETIDS) <= 0) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1175) NAME1,NAME2,SETID,J,SPCADD_SIDS(NSPCADD,NUM_SETIDS)
                  WRITE(F06,1175) NAME1,NAME2,SETID,J,SPCADD_SIDS(NSPCADD,NUM_SETIDS)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Read and check data on optional continuation cards

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
               ELSE
                  NUM_SETIDS = NUM_SETIDS + 1
                  IF (NUM_SETIDS > LSPCADDC) THEN
                     WRITE(ERR,1139) SETID,NAME1,NAME2,LSPCADDC
                     WRITE(F06,1139) SETID,NAME1,NAME2,LSPCADDC
                     CALL OUTA_HERE ( 'Y' )                        ! Coding error, so quit
                  ENDIF
                  CALL I4FLD ( JCARD(J), JF(J),  SPCADD_SIDS(NSPCADD,NUM_SETIDS) )
                  DO K=1,NUM_SETIDS-1                       ! Check for duplicate set ID's
                     IF (SPCADD_SIDS(NSPCADD,NUM_SETIDS) == SPCADD_SIDS(NSPCADD,K)) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1140) SPCADD_SIDS(NSPCADD,NUM_SETIDS),NAME2,SETID
                        WRITE(ERR,1140) SPCADD_SIDS(NSPCADD,NUM_SETIDS),NAME2,SETID
                     ENDIF
                  ENDDO
                  IF (IERRFL(J) == 'N') THEN
                     IF (SPCADD_SIDS(NSPCADD,NUM_SETIDS) <= 0) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1175) NAME1,NAME2,SETID,J,SPCADD_SIDS(NSPCADD,NUM_SETIDS)
                        WRITE(F06,1175) NAME1,NAME2,SETID,J,SPCADD_SIDS(NSPCADD,NUM_SETIDS)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9)! Make sure that there are no imbedded blanks in fields 2-9
            CALL CRDERR ( CARD )                           ! CRDERR prints errors found when reading fields

            CYCLE
         ELSE
            EXIT
         ENDIF

      ENDDO



      RETURN

! **********************************************************************************************************************************
 1139 FORMAT(' *ERROR  1139: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY PAIRS OF ',A,' SET IDs ON BULK DATA ',A,' ENTRY WITH SET ID = ',I8                           &
                    ,/,14X,' LIMIT IS = ',I12)

 1140 FORMAT(' *ERROR  1140: A DUPLICATE SET ID = ',I8,' WAS FOUND ON ',A,' BULK DATA ENTRY ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1175 FORMAT(' *ERROR  1175: ',A,' SET ID ON BULK DATA ',A,' ENTRY ID = ',I8,' IN FIELD ',I3,' IS = ',I8,'. MUST BE > 0')

! **********************************************************************************************************************************

      END SUBROUTINE BD_SPCADD


      SUBROUTINE BD_SPCADD0 ( CARD, LARGE_FLD_INP, ISPCADD )

! Processes SPCADD Bulk Data Cards to count the number of SPC (or SPC1) set ID's on this logical SPCADD card. The calling routine
! determines the max number od set ID's over all SPCADD cards in the data deck

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  SPCADD_SIDS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_SPCADD0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG), INTENT(OUT)      :: ISPCADD           ! Count of number of SPC or SPC1 set ID's defined on the SPCADD
      INTEGER(LONG)                   :: J                 ! DO loop index




! **********************************************************************************************************************************
! SPCADD Bulk Data Card routine

!   FIELD   ITEM            ARRAY ELEMENT
!   -----   ------------    -------------
!    2      SPCADD set ID    SPCADD(nspcadd, 1)
!    3      SPC/SPC1 set ID  SPCADD(nspcadd, 2)
!    4      SPC/SPC1 set ID  SPCADD(nspcadd, 3)
!    5      SPC/SPC1 set ID  SPCADD(nspcadd, 4)
!    6      SPC/SPC1 set ID  SPCADD(nspcadd, 5)
!    7      SPC/SPC1 set ID  SPCADD(nspcadd, 6)
!    8      SPC/SPC1 set ID  SPCADD(nspcadd, 7)
!    9      SPC/SPC1 set ID  SPCADD(nspcadd, 8)

! 1st continuation card:
!
!    2      SPC/SPC1 set ID  SPCADD(nspcadd, 9)
!    3      SPC/SPC1 set ID  SPCADD(nspcadd,10)
!    4      SPC/SPC1 set ID  SPCADD(nspcadd,11)
!    5      SPC/SPC1 set ID  SPCADD(nspcadd,12)
!    6      SPC/SPC1 set ID  SPCADD(nspcadd,13)
!    7      SPC/SPC1 set ID  SPCADD(nspcadd,14)
!    8      SPC/SPC1 set ID  SPCADD(nspcadd,15)
!    9      SPC/SPC1 set ID  SPCADD(nspcadd,16)

! Subsequent con't cards follow the same patterm as the 1st


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! ISPCADD will count the number of SPC or SPC1 set ID's on this SPCADD card. It starts with 1 since there is the SPCADD set ID.

! Count SPC/SPC1 set ID's on parent card

      ISPCADD = 1
      DO J=3,9
         IF (JCARD(J)(1:) == ' ') THEN
            CYCLE
         ELSE
            ISPCADD = ISPCADD + 1
         ENDIF
      ENDDO

! Count SPC/SPC1 set ID's on optional continuation cards

      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC0  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN
            DO J=2,9
               IF (JCARD(J)(1:) == ' ') THEN
                  CYCLE
               ELSE
                  ISPCADD = ISPCADD + 1
               ENDIF
            ENDDO
         ELSE
            EXIT
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_SPCADD0


      SUBROUTINE BD_MPC ( CARD, LARGE_FLD_INP, CC_MPC_FND )

! Processes MPC Bulk Data Cards. Writes MPC card data to file L1S

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1S
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LMPC, LSUB, MMPC, NMPC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  MPCSET, MPC_SIDS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_MPC'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_MPC_FND        ! ='Y' if this MPC is a set requested in Case Control
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: MPC_COMP(MMPC)    ! Array of GRID displ comp values (1-6) found on this logical MPC card
      INTEGER(LONG)                   :: MPC_GRID(MMPC)    ! Array of GRID ID's found on this logical MPC card
      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: NUM_TRIPLES       ! Counter on number of pairs of grid/comp/coeff triplets on this MPC
!                                                            card. Must be <= MMPC which was counted in subr BD_MPC0
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      INTEGER(LONG)                   :: SETID             ! Set ID for this LOAD Bulk Data card


      REAL(DOUBLE)                    :: MPC_COEFF(MMPC)   ! Array of MPC coeff values found on this MPC logical card



! **********************************************************************************************************************************
! MPC Bulk Data Card:

!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    2      SID             MPC set ID
!    3      DEP_GRD         Dependent grid
!    4      DCOMP           Component for dependent grid
!    5      DVAL            Coeff (value) for dep grid/comp
!    6      IND_GRD         1st independent grid
!    7      ICOMP           Component for 1st indep grid
!    8      IVAL            Coeff (value) for 1st indep grid/comp

! 1st continuation card:
!
!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    3      IND_GRD         2nd independent grid
!    4      ICOMP           Component for 2nd indep grid
!    5      IVAL            Coeff (value) for 2nd indep grid/comp
!    6      IND_GRD         3rd independent grid
!    7      ICOMP           Component for 3rd indep grid
!    8      IVAL            Coeff (value) for 3rd indep grid/comp

! Subsequent con't cards follow the same patterm as the 1st


! Initialize arrays

      DO J=1,MMPC
         MPC_GRID(J)  = 0
         MPC_COMP(J)  = 0
         MPC_COEFF(J) = ZERO
      ENDDO

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NMPC = NMPC+1
!xx   IF (NMPC > LMPC) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LMPC
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LMPC
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read MPC set ID in field 2

      CALL I4FLD ( JCARD(2), JF(2), SETID )
      IF (IERRFL(2) == 'N') THEN
         IF (SETID == MPCSET) THEN
            CC_MPC_FND = 'Y'
         ENDIF
         MPC_SIDS(NMPC) = SETID
      ELSE
         JERR = JERR + 1
      ENDIF

! Read and check dependent grid/comp/coeff in fields 3, 4, 5 on parent card. If fields 3, 4 or 5 are blank then dep
! grid, comp, or coeff was not input which is an error (must have dependent grid/comp/coeff)

      NUM_TRIPLES = 0

      IF ((JCARD(3)(1:) == ' ') .OR. (JCARD(4)(1:) == ' ') .OR. (JCARD(5)(1:) == ' ')) THEN
         WRITE(ERR,1126)
         WRITE(F06,1126)
         FATAL_ERR = FATAL_ERR + 1
      ELSE                                                 ! Increment NUM_TRIPLES and check for overflow
         NUM_TRIPLES = NUM_TRIPLES + 1
         IF (NUM_TRIPLES > MMPC) THEN
            WRITE(ERR,1120) SUBR_NAME,JCARD(1),LMPC
            WRITE(F06,1120) SUBR_NAME,JCARD(1),LMPC
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         CALL I4FLD (JCARD(3),JF(3),MPC_GRID(NUM_TRIPLES)) ! Read parent card field 3: dependent grid
         IF (IERRFL(3) == 'Y') THEN
            JERR = JERR + 1
         ENDIF

         CALL I4FLD ( JCARD(4), JF(4), I4INP )             ! Read parent card field 4: components at dependent grid
         IF ( IERRFL(4) == 'N' ) THEN
            IF ((I4INP >= 0) .AND. (I4INP <= 6)) THEN
               MPC_COMP(NUM_TRIPLES) = I4INP
            ELSE
               JERR      = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1124) JF(4),JCARD(1),JCARD(2),JF(4),JCARD(4)
               WRITE(F06,1124) JF(4),JCARD(1),JCARD(2),JF(4),JCARD(4)
            ENDIF
         ELSE
            JERR = JERR + 1
         ENDIF

         CALL R8FLD (JCARD(5),JF(5),MPC_COEFF(NUM_TRIPLES))! Read parent card field 5: dependent grid MPC coefficient
         IF (IERRFL(5) == 'Y') THEN
            JERR = JERR + 1
         ENDIF

      ENDIF

! Read and check independent grid/comp/coeff in fields 6, 7, 8 on parent card. We skip these 3 fields if all are blank

      IF ((JCARD(6)(1:) == ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN
         CONTINUE
      ELSE                                                 ! Increment NUM_TRIPLES and check for overflow
         NUM_TRIPLES = NUM_TRIPLES + 1
         IF (NUM_TRIPLES > MMPC) THEN
            WRITE(ERR,1120) SUBR_NAME,JCARD(1),LMPC
            WRITE(F06,1120) SUBR_NAME,JCARD(1),LMPC
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         CALL I4FLD (JCARD(6),JF(6),MPC_GRID (NUM_TRIPLES))! Read parent card field 6: independent grid
         IF (IERRFL(6) == 'Y') THEN
            JERR = JERR + 1
         ENDIF

         CALL I4FLD ( JCARD(7), JF(7), I4INP )             ! Read parent card field 7: components at independent grid
         IF ( IERRFL(7) == 'N' ) THEN
            IF ((I4INP >= 0) .AND. (I4INP <= 6)) THEN
               MPC_COMP(NUM_TRIPLES) = I4INP
            ELSE
               JERR      = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1124) JF(7),JCARD(1),JCARD(2),JF(7),JCARD(7)
               WRITE(F06,1124) JF(7),JCARD(1),JCARD(2),JF(7),JCARD(7)
            ENDIF
         ELSE
            JERR = JERR + 1
         ENDIF

         CALL R8FLD (JCARD(8),JF(8),MPC_COEFF(NUM_TRIPLES))! Read parent card field 8: independent grid MPC coefficient
         IF (IERRFL(8) == 'Y') THEN
            JERR = JERR + 1
         ENDIF

      ENDIF

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,0,5,6,0,8,0 )     ! Make sure there are no imbedded blanks (except 4,7: comps & 9: blank)
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )   ! Issue warning if field 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Read and check data on optional continuation cards

      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN

            DO J=1,2

               IF ((JCARD(3*J)(1:) == ' ') .AND. (JCARD(3*J+1)(1:) == ' ') .AND. (JCARD(3*J+2)(1:) == ' ')) THEN
                  CONTINUE
               ELSE                                        ! Increment NUM_TRIPLES and check for overflow
                  NUM_TRIPLES = NUM_TRIPLES + 1
                  IF (NUM_TRIPLES > MMPC) THEN
                     WRITE(ERR,1120) SUBR_NAME,JCARD(1),LMPC
                     WRITE(F06,1120) SUBR_NAME,JCARD(1),LMPC
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF

                  CALL I4FLD (JCARD(3*J),JF(3*J),MPC_GRID (NUM_TRIPLES))! Read independent grid
                  IF (IERRFL(3*J) == 'Y') THEN
                     JERR = JERR + 1
                  ENDIF
                                                           ! Read cont card: components at independent grid
                  CALL I4FLD ( JCARD(3*J+1), JF(3*J+1), I4INP )
                  IF ( IERRFL(3*J+1) == 'N' ) THEN
                     IF ((I4INP >= 0) .AND. (I4INP <= 6)) THEN
                        MPC_COMP(NUM_TRIPLES) = I4INP
                     ELSE
                        JERR      = JERR + 1
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1124) JF(3*J+1),JCARD(1),JCARD(2),JF(3*J+1),JCARD(3*J+1)
                        WRITE(F06,1124) JF(3*J+1),JCARD(1),JCARD(2),JF(3*J+1),JCARD(3*J+1)
                     ENDIF
                  ELSE
                     JERR = JERR + 1
                  ENDIF

                  CALL R8FLD (JCARD(3*J+2),JF(3*J+2),MPC_COEFF(NUM_TRIPLES))! Read independent grid MPC coefficient
                  IF (IERRFL(3*J+2) == 'Y') THEN
                     JERR = JERR + 1
                  ENDIF

               ENDIF

            ENDDO

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,0,5,6,0,8,0)! Make sure there are no imbedded blanks (except 4,7: comps & 9: blank)
            CALL CARD_FLDS_NOT_BLANK(JCARD,0,0,0,0,0,0,0,9)! Issue warning if field 9 not blank
            CALL CRDERR ( CARD )                           ! CRDERR prints errors found when reading fields

            CYCLE
         ELSE
            EXIT
         ENDIF

      ENDDO

! Write data to file L1S

      IF (JERR == 0) THEN
         WRITE(L1S) SETID
         WRITE(L1S) NUM_TRIPLES
         DO J=1,NUM_TRIPLES
            WRITE(L1S) MPC_GRID(J), MPC_COMP(J), MPC_COEFF(J)
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1120 FORMAT(' *ERROR  1120: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY TRIPLETS OF GRID/COMPONENT/MPC COEFF ON BULK DATA MPC CARD WITH SET ID = ',I8                &
                    ,/,14X,' LIMIT IS = ',I12)

 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')


 1126 FORMAT(' *ERROR  1126: FIELDS 3, 4 AND 5 ON MPC ENTRY (DEFINING THE DEPENDENT GRID, COMPONENT, MPC COEFF) MUST NOT BE BLANK')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_MPC


      SUBROUTINE BD_MPC0 ( CARD, LARGE_FLD_INP, IMPC )

! Processes MPC Bulk Data Cards to determine the number of triplets of grid/comp/coeff on logical MPC B.D.card. Blank entries
! for triplets are not counted.
! The number of triplets defined on this MPC card will be returned to the calling routine so that the max number of triplets
! over all MPC cards can be determined.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_MPC0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG), INTENT(OUT)      :: IMPC              ! Count of number of grid/comp/coeff triplets on this MPC logical card
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index




! **********************************************************************************************************************************
! MPC Bulk Data Card:

!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    2      SID             MPC set ID
!    3      DEP_GRD         Dependent grid
!    4      DCOMP           Component for dependent grid
!    5      DVAL            Coeff (value) for dep grid/comp
!    6      IND_GRD         1st independent grid
!    7      ICOMP           Component for 1st indep grid
!    8      IVAL            Coeff (value) for 1st indep grid/comp

!
!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    3      IND_GRD         2nd independent grid
!    4      ICOMP           Component for 2nd indep grid
!    5      IVAL            Coeff (value) for 2nd indep grid/comp
!    6      IND_GRD         3rd independent grid
!    7      ICOMP           Component for 3rd indep grid
!    8      IVAL            Coeff (value) for 3rd indep grid/comp

! Subsequent continuation cards follow the same pattern as the 1st continuation card


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! IMPC will count the number of grid/comp/coeff triplets. We will start with 1 since there must be a dependent grid/comp/coeff
! defined in fields 3, 4, 5 and count the independent triplets (the 1st of which is on the parent card in fields 6, 7, 8)

      IMPC = 1
      IF ((JCARD(6)(1:) /= ' ') .OR. (JCARD(7)(1:) /= ' ') .OR. (JCARD(8)(1:) /= ' ')) THEN
         IMPC = IMPC + 1
      ENDIF

! Count triplets of grid/comp/coeff on optional continuation cards

      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC0  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN
            DO J=3,8,3                                     ! Count up to 2 triplets per cont card (fields 3,4,5 and fields 6,7,8)
               IF ((JCARD(J)(1:) /= ' ') .OR. (JCARD(J+1)(1:) /= ' ') .OR. (JCARD(J+2)(1:) /= ' ')) THEN
                  IMPC = IMPC + 1
               ELSE
                  CYCLE
               ENDIF
            ENDDO
         ELSE
            EXIT
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_MPC0


      SUBROUTINE BD_MPCADD ( CARD, LARGE_FLD_INP, CC_MPC_FND )

! Processes MPCADD Bulk Data Cards. Reads and checks data and enters data into array MPCADD_SIDS

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LMPCADDR, LSUB, NMPCADD, LMPCADDC, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  MPCADD_SIDS, MPCSET, SUBLOD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_MPCADD'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_MPC_FND        ! 'Y' if B.D  card w/ same set ID as C.C. MPC = SID
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER( 3*BYTE)              :: NAME1   = 'MPC'   ! Name for output error message use
      CHARACTER( 6*BYTE)              :: NAME2   = 'MPCADD'! Name for output error message use

      INTEGER(LONG)                   :: I,J,K             ! DO loop index
      INTEGER(LONG)                   :: NUM_SETIDS        ! Counter on number of set ID's on MPCADD card. The SID defining this
!                                                            MPCADD card is counted  as the first one.
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: SETID             ! Set ID for this MPCADD Bulk Data card




! **********************************************************************************************************************************
! MPCADD Bulk Data Card routine

!   FIELD   ITEM            ARRAY ELEMENT
!   -----   ------------    -------------
!    2      MPCADD set ID    MPCADD(nspcadd, 1)
!    3      MPC/MPC1 set ID  MPCADD(nspcadd, 2)
!    4      MPC/MPC1 set ID  MPCADD(nspcadd, 3)
!    5      MPC/MPC1 set ID  MPCADD(nspcadd, 4)
!    6      MPC/MPC1 set ID  MPCADD(nspcadd, 5)
!    7      MPC/MPC1 set ID  MPCADD(nspcadd, 6)
!    8      MPC/MPC1 set ID  MPCADD(nspcadd, 7)
!    9      MPC/MPC1 set ID  MPCADD(nspcadd, 8)

! 1st continuation card:
!
!    2      MPC/MPC1 set ID  MPCADD(nspcadd, 9)
!    3      MPC/MPC1 set ID  MPCADD(nspcadd,10)
!    4      MPC/MPC1 set ID  MPCADD(nspcadd,11)
!    5      MPC/MPC1 set ID  MPCADD(nspcadd,12)
!    6      MPC/MPC1 set ID  MPCADD(nspcadd,13)
!    7      MPC/MPC1 set ID  MPCADD(nspcadd,14)
!    8      MPC/MPC1 set ID  MPCADD(nspcadd,15)
!    9      MPC/MPC1 set ID  MPCADD(nspcadd,16)

! Subsequent con't cards follow the same patterm as the 1st


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NMPCADD = NMPCADD+1
!xx   IF (NMPCADD > LMPCADDR) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LMPCADDR
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LMPCADDR
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data on parent card

      CALL I4FLD ( JCARD(2), JF(2), SETID )                ! Read set ID for this MPCADD Bulk Data card
      IF (IERRFL(2) == 'N') THEN
         MPCADD_SIDS(NMPCADD,1) = SETID
         DO I=1,NSUB
            IF (SETID == MPCSET) THEN
               CC_MPC_FND = 'Y'
            ENDIF
         ENDDO
      ENDIF

      NUM_SETIDS = 1                                       ! Read MPCADD ID's on parent card.
      DO J=3,9
         IF (JCARD(J)(1:) == ' ') THEN
            CYCLE                                          ! CYCLE if a set ID field is blank
         ELSE
            NUM_SETIDS = NUM_SETIDS + 1
            IF (NUM_SETIDS > LMPCADDC) THEN
               WRITE(ERR,1139) SETID,NAME1,NAME2,LMPCADDC
               WRITE(F06,1139) SETID,NAME1,NAME2,LMPCADDC
               CALL OUTA_HERE ( 'Y' )                      ! Coding error, so quit
            ENDIF
            CALL I4FLD ( JCARD(J), JF(J),  MPCADD_SIDS(NMPCADD,NUM_SETIDS) )
            DO K=1,NUM_SETIDS-1                            ! Check for duplicate set ID's
               IF (MPCADD_SIDS(NMPCADD,NUM_SETIDS) == MPCADD_SIDS(NMPCADD,K)) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1140) MPCADD_SIDS(NMPCADD,NUM_SETIDS),NAME2,SETID
                  WRITE(F06,1140) MPCADD_SIDS(NMPCADD,NUM_SETIDS),NAME2,SETID
               ENDIF
            ENDDO
            IF (IERRFL(J) == 'N') THEN
               IF (MPCADD_SIDS(NMPCADD,NUM_SETIDS) <= 0) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1175) NAME1,NAME2,SETID,J,MPCADD_SIDS(NMPCADD,NUM_SETIDS)
                  WRITE(F06,1175) NAME1,NAME2,SETID,J,MPCADD_SIDS(NMPCADD,NUM_SETIDS)
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Read and check data on optional continuation cards

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
               ELSE
                  NUM_SETIDS = NUM_SETIDS + 1
                  IF (NUM_SETIDS > LMPCADDC) THEN
                     WRITE(ERR,1139) SETID,NAME1,NAME2,LMPCADDC
                     WRITE(F06,1139) SETID,NAME1,NAME2,LMPCADDC
                     CALL OUTA_HERE ( 'Y' )                        ! Coding error, so quit
                  ENDIF
                  CALL I4FLD ( JCARD(J), JF(J),  MPCADD_SIDS(NMPCADD,NUM_SETIDS) )
                  DO K=1,NUM_SETIDS-1                       ! Check for duplicate set ID's
                     IF (MPCADD_SIDS(NMPCADD,NUM_SETIDS) == MPCADD_SIDS(NMPCADD,K)) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1140) MPCADD_SIDS(NMPCADD,NUM_SETIDS),NAME2,SETID
                        WRITE(ERR,1140) MPCADD_SIDS(NMPCADD,NUM_SETIDS),NAME2,SETID
                     ENDIF
                  ENDDO
                  IF (IERRFL(J) == 'N') THEN
                     IF (MPCADD_SIDS(NMPCADD,NUM_SETIDS) <= 0) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1175) NAME1,NAME2,SETID,J,MPCADD_SIDS(NMPCADD,NUM_SETIDS)
                        WRITE(F06,1175) NAME1,NAME2,SETID,J,MPCADD_SIDS(NMPCADD,NUM_SETIDS)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9)! Make sure that there are no imbedded blanks in fields 2-9
            CALL CRDERR ( CARD )                           ! CRDERR prints errors found when reading fields

            CYCLE
         ELSE
            EXIT
         ENDIF

      ENDDO


      RETURN

! **********************************************************************************************************************************
 1139 FORMAT(' *ERROR  1139: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY PAIRS OF ',A,' SET IDs ON BULK DATA ',A,' ENTRY WITH SET ID = ',I8                           &
                    ,/,14X,' LIMIT IS = ',I12)

 1140 FORMAT(' *ERROR  1140: A DUPLICATE SET ID = ',I8,' WAS FOUND ON ',A,' BULK DATA ENTRY ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1175 FORMAT(' *ERROR  1175: ',A,' SET ID ON BULK DATA ',A,' ENTRY ID = ',I8,' IN FIELD ',I3,' IS = ',I8,'. MUST BE > 0')

! **********************************************************************************************************************************

      END SUBROUTINE BD_MPCADD


      SUBROUTINE BD_MPCADD0 ( CARD, LARGE_FLD_INP, IMPCADD )

! Processes MPCADD Bulk Data Cards to count the number of MPC set ID's on this logical MPCADD card. The calling routine
! determines the max number of set ID's over all MPCADD cards in the data deck

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_MPCADD0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG), INTENT(OUT)      :: IMPCADD           ! Count of number of MPC set ID's defined on the MPCADD
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index




! **********************************************************************************************************************************
! MPCADD Bulk Data Card routine

!   FIELD   ITEM            ARRAY ELEMENT
!   -----   ------------    -------------
!    2      MPCADD set ID    MPCADD(nmpcadd, 1)
!    3      MPC MPC set ID   MPCADD(nmpcadd, 2)
!    4      MPC MPC set ID   MPCADD(nmpcadd, 3)
!    5      MPC MPC set ID   MPCADD(nmpcadd, 4)
!    6      MPC MPC set ID   MPCADD(nmpcadd, 5)
!    7      MPC MPC set ID   MPCADD(nmpcadd, 6)
!    8      MPC MPC set ID   MPCADD(nmpcadd, 7)
!    9      MPC MPC set ID   MPCADD(nmpcadd, 8)

! 1st continuation card:
!
!    2      MPC MPC set ID   MPCADD(nmpcadd, 9)
!    3      MPC MPC set ID   MPCADD(nmpcadd,10)
!    4      MPC MPC set ID   MPCADD(nmpcadd,11)
!    5      MPC MPC set ID   MPCADD(nmpcadd,12)
!    6      MPC MPC set ID   MPCADD(nmpcadd,13)
!    7      MPC MPC set ID   MPCADD(nmpcadd,14)
!    8      MPC MPC set ID   MPCADD(nmpcadd,15)
!    9      MPC MPC set ID   MPCADD(nmpcadd,16)

! Subsequent con't cards follow the same patterm as the 1st


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! IMPCADD will count the number of MPC set ID's on this MPCADD card. It starts with 1 since there is the MPCADD set ID.

! Count MPC set ID's on parent card

      IMPCADD = 1
      DO J=3,9
         IF (JCARD(J)(1:) == ' ') THEN
            CYCLE
         ELSE
            IMPCADD = IMPCADD + 1
         ENDIF
      ENDDO

! Count MPC set ID's on optional continuation cards

      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC0  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN
            DO J=2,9
               IF (JCARD(J)(1:) == ' ') THEN
                  CYCLE
               ELSE
                  IMPCADD = IMPCADD + 1
               ENDIF
            ENDDO
         ELSE
            EXIT
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_MPCADD0

   END MODULE CONSTRAINT_CARDS
