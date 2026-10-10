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

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, L1N
      USE SCONTR, ONLY                :  FATAL_ERR, JCARD_LEN, JF, NAOCARD, BLNK_SUB_NAM
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  I4FLD, IP6CHK

      USE BDF_ID_LISTS, ONLY          :  READ_ID_LIST

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_ASET1'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of 8 characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER(LEN(TSET_CHR_LEN))    :: SET               ! 'A' or 'O' depending on whether the B.D card is ASET or OMIT

      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0     ! Error indicator for several types of error in format #2 of input
      INTEGER(LONG)                   :: COMPJ     = 0     ! Displ component(s)  read from a B.D. ASET/OMIT card




      INTEGER(LONG), ALLOCATABLE      :: RANGES(:,:)       ! The grid IDs of the list: single (G, G) or ranges (G1, G2)
      INTEGER(LONG)                   :: NRANGE            ! Number of items in RANGES
      INTEGER(LONG)                   :: NBAD              ! Number of list fields in error

! **********************************************************************************************************************************
! ASET1, OMIT1 Bulk Data Card routine

!   FIELD   ITEM
!   -----   ------------
!    2      COMPJ, Displ component(s)
!   3-9     Grid IDs, and fields 2-9 of optional continuations: single IDs and ranges "G1 THRU G2" in any mix

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Find out if card is ASET1 or OMIT1:

      IF      (JCARD(1)(1:5) == 'ASET1') THEN
         SET = 'A '
      ELSE IF (JCARD(1)(1:5) == 'OMIT1') THEN
         SET = 'O '
      ENDIF

! Get the components (and put into integer COMPJ)

      JERR = 0
      CALL IP6CHK ( JCARD(2), JCARDO, IP6TYP, IDUM )
      IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'ZERO    ') .OR. (IP6TYP == 'BLANK   ')) THEN
         CALL I4FLD ( JCARDO, JF(2), COMPJ )
      ELSE
         JERR      = JERR + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1123) JF(2),JCARD(1),JF(2),JCARD(2)
         WRITE(F06,1123) JF(2),JCARD(1),JF(2),JCARD(2)
      ENDIF

! Fields 3-9 and the continuations: the grid IDs, single or as ranges "G1 THRU G2" in any mix (subr READ_ID_LIST). All records
! are written to file LINK1N only if the whole entry is correct.

      CALL READ_ID_LIST ( CARD, LARGE_FLD_INP, 3_LONG, .FALSE., NRANGE, RANGES, NBAD )

      IF ((NRANGE == 0) .AND. (NBAD == 0)) THEN            ! No grids on the entry
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1125) 'GRID POINT', JF(3), JCARD(1)
         WRITE(F06,1125) 'GRID POINT', JF(3), JCARD(1)
      ENDIF

      IF ((JERR == 0) .AND. (NBAD == 0)) THEN
         DO J=1,NRANGE
            WRITE(L1N) COMPJ,RANGES(1,J),RANGES(2,J),SET   ! A single grid G is written as the range G, G
            NAOCARD = NAOCARD + 1
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1123 FORMAT(' *ERROR  1123: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' CARD. MUST BE A COMBINATION OF DIGITS 1-6'                &
                    ,/,14X,' HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1125 FORMAT(' *ERROR  1125: NO ',A,' SPECIFIED IN FIELD',I4,' ON ',A,' CARD')


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

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, L1X
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NUM_USET_RECORDS
      USE CONSTANTS_1, ONLY           :  ZERO

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, I4FLD, IP6CHK, LEFT_ADJ_BDFLD

      USE BDF_ID_LISTS, ONLY          :  READ_ID_LIST

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_USET1'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN=JCARD_LEN)        :: CHRFLD            ! A character field from the entry
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER( 2*BYTE)              :: USET_NAME         ! Name in field 2 of the USET entry

      INTEGER(LONG)                   :: COMPJ     = 0     ! DOF's constrained at GRIDJ
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count




      INTEGER(LONG), ALLOCATABLE      :: RANGES(:,:)       ! The grid IDs of the list: single (G, G) or ranges (G1, G2)
      INTEGER(LONG)                   :: NRANGE            ! Number of items in RANGES
      INTEGER(LONG)                   :: NBAD              ! Number of list fields in error

! **********************************************************************************************************************************
!  USET1 Bulk Data Card routine

!    FIELD   ITEM           ARRAY ELEMENT
!    -----   ------------   -------------
!     2      USET name (e.g. "U1", "U2")
!     3      Component numbers
!     4-9    Grid IDs, and fields 2-9 of optional continuations: single IDs and ranges "G1 THRU G2" in any mix


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

! Fields 4-9 and the continuations: the grid IDs, single or as ranges "G1 THRU G2" in any mix (subr READ_ID_LIST). All records
! are written to file LINK1X only if the whole entry is correct.

      CALL BD_IMBEDDED_BLANK ( JCARD,2,0,0,0,0,0,0,0 )     ! Make sure that there are no imbedded blanks in field 2
      CALL READ_ID_LIST ( CARD, LARGE_FLD_INP, 4_LONG, .FALSE., NRANGE, RANGES, NBAD )

      IF ((NRANGE == 0) .AND. (NBAD == 0)) THEN            ! No grids were specified on the USET1 entry, so error
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1100) JCARD(1)
         WRITE(F06,1100) JCARD(1)
      ENDIF

      IF ((JERR == 0) .AND. (NBAD == 0)) THEN
         DO J=1,NRANGE
            NUM_USET_RECORDS = NUM_USET_RECORDS + 1        ! Incr count of number of entries written to file LINK1X
            WRITE(L1X) USET_NAME, COMPJ, RANGES(1,J), RANGES(2,J)
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1100 FORMAT(' *ERROR  1100: ',A,' ENTRY MUST HAVE AT LEAST 1 GRID DEFINED (IN FIELDS 4-9 OF PARENT OR 2-9 OF CONTINUATION)')

 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

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
