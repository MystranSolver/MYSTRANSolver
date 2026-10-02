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

   MODULE RIGID_ELEMENTS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BD_RBAR, BD_RBE1, BD_RBE2, BD_RBE3, BD_RBE30, BD_RSPLINE, BD_RSPLINE0

   CONTAINS

      SUBROUTINE BD_RBAR ( CARD )

! Processes RBAR Bulk Data Cards. Writes RBAR element records to file L1F for later processing.
! Two records are written for each dependent Grid/DOF pair:

!       1) Record 1 has the rigid element type: 'RBAR    '
!       2) Record 2 has:
!            RELID: Rigid element ID
!            GID1 : 1st Grid ID
!            IDOF1: Independent DOF's at GID1
!            DDOF1: Dependent   DOF's at GID1
!            GID2 : 2nd Grid ID
!            IDOF2: Independent DOF's at GID2
!            DDOF2: Dependent   DOF's at GID2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1F
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LRIGEL, NRBAR, NRIGEL, NRECARD
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  RIGID_ELEM_IDS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_RBAR'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER( 8*BYTE)              :: IP6TYP(5:8)       ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER( 8*BYTE), PARAMETER   :: RTYPE = 'RBAR    '! Rigid element type

      INTEGER(LONG)                   :: DOF_NOS(6)= 0     ! Fields 5 & 6 should fill this in to be: 1 2 3 4 5 6
      INTEGER(LONG)                   :: GID1      = 0     ! Grid ID for 1st RBAR grid
      INTEGER(LONG)                   :: GID2      = 0     ! Grid ID for 2nd RBAR grid
      INTEGER(LONG)                   :: INT1      = 0     ! An integer 1-6 read from internal file
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      INTEGER(LONG)                   :: IDOF_ERR  = 0     ! Count of the no. of DOF component errors in fields 5,6
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: RBDOF(4)          ! The DOF's in fields 5,6,7,8
      INTEGER(LONG)                   :: RELID     = 0     ! Rigid element ID




! **********************************************************************************************************************************
! RBAR Bulk Data Card routine

!   FIELD   ITEM
!   -----   ------------
!    2      RELID, Rigid Elem ID
!    3      GID1 , Grid ID for 1st RBAR grid
!    4      GID2 , Grid ID for 2nd RBAR grid
!    5      IDOF1, Indep. DOF's at GID1
!    6      IDOF2, Indep. DOF's at GID2
!    7      DDOF1, Depen. DOF's at GID1
!    8      DDOF2, Depen. DOF's at GID2


! Data is written to file LINK1F for later processing after checks on format of data.

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NRBAR  = NRBAR+1
      NRIGEL = NRIGEL + 1
!xx   IF (NRIGEL > LRIGEL) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LRIGEL
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LRIGEL
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data

      CALL I4FLD ( JCARD(2), JF(2), RELID )                 ! Read rigid element ID in field 2
      IF (IERRFL(2) /= 'N') THEN
         JERR = JERR + 1
      ELSE
         RIGID_ELEM_IDS(NRIGEL) = RELID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), GID1 )                 ! Read 1st grid in field 3
      IF (IERRFL(3) /= 'N') THEN
         JERR = JERR + 1
      ENDIF

      CALL I4FLD ( JCARD(4), JF(4), GID2 )                 ! Read 2nd grid in field 4
      IF (IERRFL(4) /= 'N') THEN
         JERR = JERR + 1
      ENDIF

      DO J=5,6                                             ! Read independent DOF's in fields 5, 6
         CALL I4FLD ( JCARD(J), JF(J), RBDOF(J-4) )
         IF (IERRFL(J) == 'N') THEN
            CALL IP6CHK ( JCARD(J), JCARDO, IP6TYP(J), IDUM )
            IF ((IP6TYP(J) == 'COMP NOS') .OR. (IP6TYP(J) == 'BLANK   ')) THEN
               CONTINUE
            ELSE
               FATAL_ERR  = FATAL_ERR + 1
               IDOF_ERR   = IDOF_ERR + 1
               WRITE(ERR,1124) J,JCARD(1),JCARD(2),J,JCARD(J)
               WRITE(F06,1124) J,JCARD(1),JCARD(2),J,JCARD(J)
            ENDIF
         ELSE
            IDOF_ERR = IDOF_ERR + 1
         ENDIF
      ENDDO

      DO J=7,8                                             ! Read dependent DOF's in fields 7, 8
         CALL I4FLD ( JCARD(J), JF(J), RBDOF(J-4) )
         IF (IERRFL(J) == 'N') THEN
            CALL IP6CHK ( JCARD(J), JCARDO, IP6TYP(J), IDUM )
            IF ((IP6TYP(J) == 'COMP NOS') .OR. (IP6TYP(J) == 'BLANK   ')) THEN
               CONTINUE
            ELSE
               FATAL_ERR  = FATAL_ERR + 1
               JERR       = JERR + 1
               WRITE(ERR,1124) J,JCARD(1),JCARD(2),J,JCARD(J)
               WRITE(F06,1124) J,JCARD(1),JCARD(2),J,JCARD(J)
            ENDIF
         ELSE
            JERR = JERR + 1
         ENDIF
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,0,0,0,0,0 )     ! Make sure that there are no imb blanks in fields 2-4. Flds 5-8 DOF's
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )   ! Issue warning if field 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! If IDOF fields 5,6 have valid DOF no's (or blank), check that the DOF's in those fields specify all 6 DOF no's

      IF (IDOF_ERR == 0) THEN
         DO J=1,JCARD_LEN
            IF (JCARD(5)(J:J) /= ' ') THEN
               READ (JCARD(5)(J:J),'(I8)') INT1
               IF (DOF_NOS(INT1) == 0) THEN
                  DOF_NOS(INT1) = INT1
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  IDOF_ERR  = IDOF_ERR + 1
               ENDIF
            ENDIF
         ENDDO
         DO J=1,JCARD_LEN
            IF (JCARD(6)(J:J) /= ' ') THEN
               READ (JCARD(6)(J:J),'(I8)') INT1
               IF (DOF_NOS(INT1) == 0) THEN
                  DOF_NOS(INT1) = INT1
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  IDOF_ERR  = IDOF_ERR + 1
               ENDIF
            ENDIF
         ENDDO
         DO J=1,6
            IF (DOF_NOS(J) == 0) THEN
               FATAL_ERR = FATAL_ERR + 1
               IDOF_ERR  = IDOF_ERR + 1
            ENDIF
         ENDDO
         IF (IDOF_ERR > 0) THEN
            WRITE(ERR,1142) RELID,JCARD(5),JCARD(6)
            WRITE(F06,1142) RELID,JCARD(5),JCARD(6)
         ENDIF
      ENDIF

! Write data to file L1F if JERR and IDOF_ERR = 0

      IF ((JERR == 0) .AND. (IDOF_ERR == 0)) THEN
         WRITE(L1F) RTYPE
         WRITE(L1F) RELID,GID1,RBDOF(1),RBDOF(3),GID2,RBDOF(2),RBDOF(4)
         NRECARD = NRECARD + 1
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1142 FORMAT(' *ERROR  1142: RBAR ELEM NUMBER ',I8,' HAS INCORRECT INDEP. DOFs IN FIELDS 5 AND/OR 6.'                              &
                    ,/,14X,' THESE 2 FIELDS MUST COMBINE TO SPECIFY DOFs 1,2,3,4,5,6 (EACH DOF ONCE, AND ONLY ONCE).'              &
                    ,/,14X,' HOWEVER, FIELD 5 WAS "',A,'", AND FIELD 6 WAS "',A,'"')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_RBAR


      SUBROUTINE BD_RBE1 ( CARD, LARGE_FLD_INP )

! Processes RBE1 Bulk Data Cards. Writes RBE1 element records to file L1F for later processing.
! Two records are written for each dependent Grid/DOF pair:

!       1) Record 1 has the rigid element type: 'RBE1    '
!       2) Record 2 has:
!             RELID           : Rigid element ID
!             IGID(i),IDOF(i) : Pairs of indep Grid/DOF (up to 6)
!             DGID, DDOF      : One pair of dependent Grid/DOF

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1F
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LRIGEL, NRBE1, NRIGEL, NRECARD
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  RIGID_ELEM_IDS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_RBE1'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: CHR_IDOF(6)
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN(JCARD))           :: JCARD1_PARENT     ! Field 1 of parent card
      CHARACTER(LEN(JCARD))           :: JCARD2_PARENT     ! Field 2 of parent card
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER( 1*BYTE)              :: MORE_IDOF = 'N'   ! = 'Y' if a cont card with indep DOF's is found
      CHARACTER( 1*BYTE)              :: UMFND             ! = 'Y' when we find "UM" in field 2 of cont card no. 1 or 2
      CHARACTER( 8*BYTE), PARAMETER   :: RTYPE = 'RBE1    '! Rigid element type

      INTEGER(LONG)                   :: CONT_NO   = 0     ! Count of the continuation cards for this parent
      INTEGER(LONG)                   :: DGID      = 0     ! A dependent grid
      INTEGER(LONG)                   :: DDOF      = 0     ! Dependent DOF's at DGID
      INTEGER(LONG)                   :: DOF_NOS(6)= 0     ! Fields 5 & 6 should fill this in to be: 1 2 3 4 5 6
      INTEGER(LONG)                   :: IDOF_ERR  = 0     ! Count of the no. of DOF component errors in fields 5,6
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IDOF(6)   = 0     ! Independent DOF's at IGID(i)
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: IGID(6)   = 0     ! Independent grid ID's
      INTEGER(LONG)                   :: INT1      = 0     ! An integer 1-6 read from internal file
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      INTEGER(LONG)                   :: JFLD1             ! A computed field number on the card
      INTEGER(LONG)                   :: JFLD2             ! A computed field number on the card
      INTEGER(LONG)                   :: NUM_IDOF_FLDS = 0 ! Number of fields that have independent DOF's specified
      INTEGER(LONG)                   :: RELID     = 0     ! This rigid elements' ID




! **********************************************************************************************************************************
! RBE1 Bulk Data Card routine

!   FIELD   ITEM
!   -----   ------------
!    2      RELID  , Rigid Elem ID
!    3      IGID(1), Grid ID for 1st independent grid
!    4      IDOF(1), DOF's at 1st independent grid
!    5      IGID(2), Grid ID for 2nd independent grid, if it exists
!    6      IDOF(2), DOF's at 2nd independent grid, if IGID(2) exists
!    7      IGID(3), Grid ID for 3rd independent grid, if it exists
!    8      IDOF(3), DOF's at 3rd independent grid, if IGID(3) exists

! Possible continuation card (if there are 4, 5 or 6 independent grids)
!   FIELD   ITEM
!   -----   ------------
!    3      IGID(4), Grid ID for 4th independent grid, if it exists
!    4      IDOF(4), DOF's at 4th independent grid, if IGID(4) exists
!    5      IGID(5), Grid ID for 5th independent grid, if it exists
!    6      IDOF(5), DOF's at 5th independent grid, if IGID(5) exists
!    7      IGID(6), Grid ID for 6th independent grid, if it exists
!    8      IDOF(6), DOF's at 6th independent grid, if IGID(6) exists
!    The collection of IDOF(i) must yield 6 DOF's that completely  describe a general rigid body motion of the
!    rigid element

! Mandatory continuation card
!   FIELD   ITEM
!   -----   ------------
!    2      "UM"
!    3      DGID, Grid ID for 1st dependent grid
!    4      DDOF, DOF's at 1st dependent grid
!    5-8    Up to 2 more pairs of DGID/DDOF, if they exist, followed by more continuation cards with up to 3 pairs
!           of DGID/DDOF in fields 3-8, if needed

! Data is written to file L1F for later processing after checks on
! format of data.

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      JCARD1_PARENT = JCARD(1)
      JCARD2_PARENT = JCARD(2)

! Check for overflow

      NRBE1  = NRBE1+1
      NRIGEL = NRIGEL+1
!xx   IF (NRIGEL > LRIGEL) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LRIGEL
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LRIGEL
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Initialize CHR_IDOF

      DO I=1,6
         CHR_IDOF(I)(1:) = ' '
      ENDDO

! Read Elem ID

      CALL I4FLD ( JCARD(2), JF(2), RELID )
      IF (IERRFL(2) /= 'N') THEN
         JERR = JERR + 1
      ELSE
         RIGID_ELEM_IDS(NRIGEL) = RELID
      ENDIF

! Read up to 3 pairs of independent grids and DOF's in fields 3-8 of the parent card.
! For any of the 3 pairs, if either of the 2 fields is blank, the other one must also be blank.

      NUM_IDOF_FLDS = 0
      DO J=1,3
         JFLD1 = 2*J+1
         JFLD2 = 2*J+2
         IF      ((JCARD(JFLD1)(1:) == ' ') .AND. (JCARD(JFLD2)(1:) == ' ')) THEN
            CYCLE
         ELSE IF ((JCARD(JFLD1)(1:) /= ' ') .AND. (JCARD(JFLD2)(1:) == ' ')) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1181) JCARD2_PARENT,JFLD1,JFLD2
            WRITE(F06,1181) JCARD2_PARENT,JFLD1,JFLD2
         ELSE IF ((JCARD(JFLD1)(1:) == ' ') .AND. (JCARD(JFLD2)(1:) /= ' ')) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1182) JCARD2_PARENT,JFLD1,JFLD2
            WRITE(F06,1182) JCARD2_PARENT,JFLD1,JFLD2
         ELSE IF ((JCARD(JFLD1)(1:) /= ' ') .AND. (JCARD(JFLD2)(1:) /= ' ')) THEN
            NUM_IDOF_FLDS = NUM_IDOF_FLDS + 1
            CHR_IDOF(NUM_IDOF_FLDS) = JCARD(JFLD2)
            CALL I4FLD ( JCARD(JFLD1), JF(JFLD1), IGID(NUM_IDOF_FLDS) )
            IF (IERRFL(JFLD1) /= 'N') THEN
               JERR = JERR + 1
            ENDIF
            CALL I4FLD ( JCARD(JFLD2), JF(JFLD2), IDOF(NUM_IDOF_FLDS) )
            IF (IERRFL(JFLD2) == 'N') THEN
               CALL IP6CHK ( JCARD(JFLD2), JCARDO, IP6TYP, IDUM )
               IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'BLANK   '))  THEN
                  CONTINUE
               ELSE
                  FATAL_ERR  = FATAL_ERR + 1
                  IDOF_ERR = IDOF_ERR + 1
                  WRITE(ERR,1124) JFLD2,JCARD1_PARENT,JCARD2_PARENT,JFLD2,JCARD(JFLD2)
                  WRITE(F06,1124) JFLD2,JCARD1_PARENT,JCARD2_PARENT,JFLD2,JCARD(JFLD2)
               ENDIF
            ELSE
               JERR = JERR + 1
            ENDIF
         ENDIF
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,0,5,0,7,0,0 )     ! Make sure no imbedded blanks in fields 2,3,5,7. Flds 4,6,8 ae DOF's
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )     ! Issue warning if field 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Now look for a continuation card. It may have up to 3 additional pairs of independent grid/DOF or it may begin the
! dependent grid/DOF list. If it is the later, then field 2 must have "UM"

      CONT_NO = 0
      UMFND   = 'N'
      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN
            CONT_NO = CONT_NO + 1
            IF (JCARD(2)(1:2) == 'UM') THEN
               UMFND = 'Y'
            ENDIF
            IF (UMFND =='N') THEN                          ! Read up to 3 additional pairs of indep grid/DOF's
               IF (CONT_NO == 1) THEN
                  MORE_IDOF = 'Y'
                  DO J=1,3
                     JFLD1 = 2*J+1
                     JFLD2 = 2*J+2
                     IF      ((JCARD(JFLD1)(1:) == ' ') .AND. (JCARD(JFLD2)(1:) == ' ')) THEN
                        CYCLE
                     ELSE IF ((JCARD(JFLD1)(1:) /= ' ') .AND. (JCARD(JFLD2)(1:) == ' ')) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1181) JCARD2_PARENT,JFLD1,JFLD2
                        WRITE(F06,1181) JCARD2_PARENT,JFLD1,JFLD2
                     ELSE IF ((JCARD(JFLD1)(1:) == ' ') .AND. (JCARD(JFLD2)(1:) /= ' ')) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1182) JCARD2_PARENT,JFLD1,JFLD2
                        WRITE(F06,1182) JCARD2_PARENT,JFLD1,JFLD2
                     ELSE IF ((JCARD(JFLD1)(1:) /= ' ') .AND. (JCARD(JFLD2)(1:) /= ' ')) THEN
                        NUM_IDOF_FLDS = NUM_IDOF_FLDS + 1
                        CHR_IDOF(NUM_IDOF_FLDS) = JCARD(JFLD2)
                        CALL I4FLD ( JCARD(JFLD1), JF(JFLD1), IGID(NUM_IDOF_FLDS) )
                        IF (IERRFL(JFLD1) /= 'N') THEN
                           JERR = JERR + 1
                        ENDIF
                        CALL I4FLD ( JCARD(JFLD2), JF(JFLD2), IDOF(NUM_IDOF_FLDS) )
                        IF (IERRFL(JFLD2) == 'N') THEN
                           CALL IP6CHK ( JCARD(JFLD2), JCARDO, IP6TYP, IDUM )
                           IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'BLANK   '))  THEN
                              CONTINUE
                           ELSE
                              FATAL_ERR  = FATAL_ERR + 1
                              IDOF_ERR = IDOF_ERR + 1
                              WRITE(ERR,1124) JFLD2,JCARD1_PARENT,JCARD2_PARENT,JFLD2,JCARD(JFLD2)
                              WRITE(F06,1124) JFLD2,JCARD1_PARENT,JCARD2_PARENT,JFLD2,JCARD(JFLD2)
                           ENDIF
                        ELSE
                           JERR = JERR + 1
                        ENDIF
                     ENDIF
                  ENDDO
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1135) JCARD1_PARENT, JCARD2_PARENT, CONT_NO
                  WRITE(F06,1135) JCARD1_PARENT, JCARD2_PARENT, CONT_NO
                  RETURN                                   ! Cannot continue, NUM_IDOF_FLDS cannot exceed 6 (for array CHR_IDOF)
               ENDIF

               CALL BD_IMBEDDED_BLANK ( JCARD,0,3,0,5,0,7,0,0 ) ! Make sure no imbedded blanks in fields 3,5,7 (4,6,8 can be IDOF's)
               CALL CARD_FLDS_NOT_BLANK ( JCARD,2,0,0,0,0,0,0,9 )
               CALL CRDERR ( CARD )                        ! CRDERR prints errors found when reading fields

            ELSE                                           ! Found "UM", so read up to 3 pairs of dep. grid/DOF on this con't card
               DO                                          ! There are any number of continuation cards w/dep. grid/DOF
                  DO J=1,3
                     JFLD1 = 2*J+1
                     JFLD2 = 2*J+2
                     IF      ((JCARD(JFLD1)(1:) == ' ') .AND. (JCARD(JFLD2)(1:) == ' ')) THEN
                        CYCLE
                     ELSE IF ((JCARD(JFLD1)(1:) /= ' ') .AND. (JCARD(JFLD2)(1:) == ' ')) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1181) JCARD2_PARENT,JFLD1,JFLD2
                        WRITE(F06,1181) JCARD2_PARENT,JFLD1,JFLD2
                     ELSE IF ((JCARD(JFLD1)(1:) == ' ') .AND. (JCARD(JFLD2)(1:) /= ' ')) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1182) JCARD2_PARENT,JFLD1,JFLD2
                        WRITE(F06,1182) JCARD2_PARENT,JFLD1,JFLD2
                     ELSE IF ((JCARD(JFLD1)(1:) /= ' ') .AND. (JCARD(JFLD2)(1:) /= ' ')) THEN
                        CALL I4FLD ( JCARD(JFLD1), JF(JFLD1), DGID )
                        IF (IERRFL(JFLD1) /= 'N') THEN
                           JERR = JERR + 1
                        ENDIF
                        CALL I4FLD ( JCARD(JFLD2), JF(JFLD2), DDOF )
                        IF (IERRFL(JFLD2) == 'N') THEN
                           CALL IP6CHK ( JCARD(JFLD2), JCARDO, IP6TYP, IDUM )
                           IF ((IP6TYP == 'COMP NOS') .OR. (IP6TYP == 'BLANK   '))  THEN
                              CONTINUE
                           ELSE
                              FATAL_ERR  = FATAL_ERR + 1
                              IDOF_ERR = IDOF_ERR + 1
                              WRITE(ERR,1124) JFLD2,JCARD1_PARENT,JCARD2_PARENT,JFLD2,JCARD(JFLD2)
                              WRITE(F06,1124) JFLD2,JCARD1_PARENT,JCARD2_PARENT,JFLD2,JCARD(JFLD2)
                           ENDIF
                        ELSE
                           JERR = JERR + 1
                        ENDIF
                        IF (JERR == 0) THEN
                           WRITE(L1F) RTYPE
                           WRITE(L1F) RELID,DGID,DDOF,NUM_IDOF_FLDS,(IGID(I),IDOF(I),I=1,NUM_IDOF_FLDS)
                           NRECARD = NRECARD + 1
                        ENDIF
                     ENDIF
                  ENDDO

                  CALL BD_IMBEDDED_BLANK ( JCARD,0,3,0,5,0,5,0,0 ) ! Make sure no imbed blanks in fields 3,5,7 (2,4,6 can be IDOF's)
                  CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )
                  CALL CRDERR ( CARD )

                  IF (LARGE_FLD_INP == 'N') THEN
                     CALL NEXTC  ( CARD, ICONT, IERR )
                  ELSE
                     CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
                     CARD = CHILD
                  ENDIF
                  CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
                  IF (ICONT == 1) THEN
                     CYCLE
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
         ELSE
            EXIT
         ENDIF
      ENDDO

      IF (UMFND == 'N') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1136) JCARD1_PARENT,JCARD2_PARENT
         WRITE(F06,1136) JCARD1_PARENT,JCARD2_PARENT
      ENDIF

! If IDOF fields all have valid DOF no's (or blank), check that the DOF's in those fields specify all 6 DOF no's

      IF (IDOF_ERR == 0) THEN
         DO I=1,NUM_IDOF_FLDS
            DO J=1,JCARD_LEN
               IF (CHR_IDOF(I)(J:J) /= ' ') THEN
                  READ (CHR_IDOF(I)(J:J),'(I8)') INT1
                  IF (DOF_NOS(INT1) == 0) THEN
                     DOF_NOS(INT1) = INT1
                  ELSE
                     FATAL_ERR = FATAL_ERR + 1
                     IDOF_ERR  = IDOF_ERR + 1
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         DO J=1,6
            IF (DOF_NOS(J) == 0) THEN
               FATAL_ERR = FATAL_ERR + 1
               IDOF_ERR  = IDOF_ERR + 1
            ENDIF
         ENDDO
         IF (IDOF_ERR > 0) THEN
            WRITE(ERR,1143) RELID,(CHR_IDOF(I),I=1,3)
            WRITE(F06,1143) RELID,(CHR_IDOF(I),I=1,3)
            IF (MORE_IDOF == 'N') THEN
               WRITE(ERR,11431)
               WRITE(F06,11431)
            ELSE
               WRITE(ERR,11432) (CHR_IDOF(I),I=4,6)
               WRITE(F06,11432) (CHR_IDOF(I),I=4,6)
            ENDIF
         ENDIF
      ENDIF

      CALL CRDERR ( CARD )



      RETURN

! **********************************************************************************************************************************
 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')

 1135 FORMAT(' *ERROR  1135: ',A,' ENTRY WITH ID = ',A,' SHOULD HAVE A MAX OF 1 CONTINUATION ENTRY WITHOUT "UM" IN FIELD 2.'       &
                    ,/,14X,' HOWEVER, CONTINUATION ENTRY NUMBER ',I8,' HAS BEEN FOUND')

 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

 1143 FORMAT(' *ERROR  1143: RBE1 ELEM NUMBER ',I8,' HAS INCORRECT INDEPENDENT DOFs IDENTIFIED.'                                   &
                    ,/,14X,' THESE (UP TO) 6 FIELDS MUST COMBINE TO SPECIFY DOFs 1,2,3,4,5,6 (EACH DOF ONCE, AND ONLY ONCE).'      &
                    ,/,14X,' HOWEVER, FIELDS 4,6,8 OF THE PARENT     ENTRY HAD: "',A,'", "',A,'", "',A,'"')

11431 FORMAT(          14X,' AND THERE WAS NO CONTINUATION ENTRY WITH ADDITIONAL INDEPENDENT DOFs IDENTIFIED')
11432 FORMAT(          14X,' AND      FIELDS 4,6,8 OF A CONTINUATION ENTRY HAD: "',A,'", "',A,'", "',A,'"')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1181 FORMAT(' *ERROR  1181: RBE1 ELEMENT NUMBER ',A,' HAS NON-BLANK GRID FIELD ',I2,' AND BLANK DOF FIELD ',I2                    &
                    ,/,14X,' THIS PAIR OF GRID/DOF MUST BE EITHER BOTH BLANK OR BOTH NON-BLANK')

 1182 FORMAT(' *ERROR  1182: RBE1 ELEMENT NUMBER ',A,' HAS BLANK GRID FIELD ',I2,' AND NON-BLANK DOF FIELD ',I2                    &
                    ,/,14X,' THIS PAIR OF GRID/DOF MUST BE EITHER BOTH BLANK OR BOTH NON-BLANK')

 1183 FORMAT(' *ERROR  1183: RBE1 ELEMENT NUMBER ',A,' HAS NO CONTINUATION ENTRY. AT LEAST 1 IS REQUIRED')

 1184 FORMAT(' *ERROR  1184: RBE1 ELEMENT NUMBER ',A,' HAS ',I8,' INDEPENDENT DOFs DEFINED. REQUIREMENT IS 6')

! **********************************************************************************************************************************

      END SUBROUTINE BD_RBE1


      SUBROUTINE BD_RBE2 ( CARD, LARGE_FLD_INP )

! Processes RBE2 Bulk Data Cards. Writes RBE2 element records to file L1F for later processing.
! Two records are written for each dependent Grid/DOF pair:

!       1) Record 1 has the rigid element type: 'RBAR    '
!       2) Record 2 has:
!            REID: Rigid element ID
!            DGID: Dependent   Grid ID
!            DDOF: Dependent   DOF's
!            IGID: Independent Grid ID

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1F
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LRIGEL, NRBE2, NRIGEL, NRECARD, NTERM_RMG
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  RIGID_ELEM_IDS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_RBE2'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER( 8*BYTE)              :: IP6TYP            ! An output from subr IP6CHK called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARDO            ! An output from subr IP6CHK called herein
      CHARACTER( 8*BYTE), PARAMETER   :: RTYPE = 'RBE2    '! Rigid element type

      INTEGER(LONG)                   :: DDOF      = 0     ! Dependent DOF's at DGID's
      INTEGER(LONG)                   :: DGID      = 0     ! Dependent grid ID's
      INTEGER(LONG)                   :: J                 ! DO loop indiex
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IDUM              ! Dummy arg in subr IP^CHK not used herein
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: IGID      = 0     ! Independent grid ID
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      INTEGER(LONG)                   :: NUM_COMP  = 0     ! Total number of components specified in DDOF
      INTEGER(LONG)                   :: RELID      = 0     ! This elements' ID




! **********************************************************************************************************************************
! RBE2 Bulk Data Card routine

!   FIELD   ITEM
!   -----   ------------
!    2      RELID, Rigid Elem ID
!    3      IGID, Independent Grid ID
!    4      DDOF, Dependent DOF's for the Grids following
!   5-9     DGID, Dependent Grid ID's

! on optional continuation cards:
!   2-9     DGID, Dependent grid ID's

! Data is written to file L1F for later processing after checks on format of data.

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NRBE2  = NRBE2+1
      NRIGEL = NRIGEL + 1
!xx   IF (NRIGEL > LRIGEL) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LRIGEL
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LRIGEL
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data

      CALL I4FLD ( JCARD(2), JF(2), RELID )                 ! Read Elem ID
      IF (IERRFL(2) /= 'N') THEN
         JERR = JERR + 1
      ELSE
         RIGID_ELEM_IDS(NRIGEL) = RELID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), IGID )                 ! Read independent Grid ID
      IF (IERRFL(3) /= 'N') THEN
         JERR = JERR + 1
      ENDIF

      NUM_COMP = 0                                         ! Get dependent DOF's in field 4 and count their number (NUM_COMP)
      CALL I4FLD ( JCARD(4), JF(4), DDOF )
      IF (IERRFL(4) == 'N') THEN
         CALL IP6CHK ( JCARD(4), JCARDO, IP6TYP, IDUM )
         IF (IP6TYP == 'COMP NOS') THEN
            DO J=1,JCARD_LEN
               IF(JCARD(4)(J:J) /= ' ') THEN
                 NUM_COMP = NUM_COMP + 1
               ENDIF
            ENDDO
         ELSE
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1124) JF(4),JCARD(1),JCARD(2),JF(4),JCARD(4)
            WRITE(F06,1124) JF(4),JCARD(1),JCARD(2),JF(4),JCARD(4)
         ENDIF
      ENDIF

      DO J=5,9
         IF (JCARD(J)(1:) == ' ') THEN
            CYCLE
         ELSE
            CALL I4FLD ( JCARD(J), JF(J), DGID )           ! Get dep. grid ID's in fields 5 - 9
            IF ((IERRFL(J) == 'N') .AND. (JERR == 0)) THEN
               WRITE(L1F) RTYPE                            ! Write element type to LINK1F
               WRITE(L1F) RELID,DGID,DDOF,IGID              ! Write data to LINK1F, one record per dependent Grid
               NRECARD = NRECARD + 1
               NTERM_RMG = NTERM_RMG + 7*NUM_COMP
            ENDIF
         ENDIF
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,0,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9,except 4(DOF)
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Read and check data on optional continuation cards with additional dep. grids

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
                  CALL I4FLD ( JCARD(J), JF(J), DGID )     ! Get dep. grid ID's in fields 5 - 9
                  IF ((IERRFL(J) == 'N') .AND. (JERR == 0)) THEN
                     WRITE(L1F) RTYPE                      ! Write element type to LINK1F
                     WRITE(L1F) RELID,DGID,DDOF,IGID        ! Write data to LINK1F, one record per dependent Grid
                     NRECARD = NRECARD + 1
                     NTERM_RMG = NTERM_RMG + 7*NUM_COMP
                  ENDIF
               ENDIF
            ENDDO

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 ) ! Make sure that there are no imbedded blanks in fields 2-9
            CALL CRDERR ( CARD )                           ! CRDERR prints errors found when reading fields

            CYCLE
         ELSE
            EXIT
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')


 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_RBE2


      SUBROUTINE BD_RBE3 ( CARD, LARGE_FLD_INP )

! Processes RBE3 Bulk Data Cards. Writes RBE3 card data to file L1F

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1F
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LRIGEL, MRBE3, NRECARD, NRIGEL
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  RIGID_ELEM_IDS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, IP6CHK, R8FLD
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK
      USE DOF_SET_CONSTRUCTION, ONLY  :  RDOF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_RBE3'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER( 1*BYTE)              :: CDOF(6)           ! An output from subr RDOF
      CHARACTER(LEN(JCARD))           :: CHRINP            ! An 8 character field returned from subr IP6CHK
      CHARACTER(LEN(JCARD))           :: CHAR_ELID         ! Char value for element ID in field 2 of parent entry
      CHARACTER( 1*BYTE)              :: FND_NEW_WGT       ! 'Y' if a weight entry has been found
      CHARACTER( 8*BYTE)              :: IP6TYP            ! Descriptor of what is in the 8 char field sent to subr IP6CHK
      CHARACTER( 1*BYTE)              :: NEXT_MUST_BE_GRID ! 'Y' if next field must be grid (i.e. need grid after wgt and comp)
      CHARACTER( 8*BYTE), PARAMETER   :: RTYPE = 'RBE3    '! Rigid element type
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! Type of TOKEN returned from subr TOKCHK

      INTEGER(LONG)                   :: CHK(2:9)          ! Array to tell which fields to check for imbedded blankS
      INTEGER(LONG)                   :: COMP(MRBE3)       ! Array of indep displ comp values (1-6) found on this logical RBE3 card
      INTEGER(LONG)                   :: GRID(MRBE3)       ! Array of indep GRID ID's found on this logical RBE3 card
      INTEGER(LONG)                   :: GRID1     = 0     ! One grid value
      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: INTDOF    = 0     ! Displ component numbers from a Ci field
      INTEGER(LONG)                   :: IRBE3             ! Count of triplets of indep grid/comp/weights read
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      INTEGER(LONG)                   :: NUM_Ci            ! Number of displ components in a Ci field
      INTEGER(LONG)                   :: REFC_NUM_Ci       ! Number of displ components in REFC field
      INTEGER(LONG)                   :: REFC      = 0     ! REFC value in field 5 of parent entry
      INTEGER(LONG)                   :: REFGRID   = 0     ! REFGRID value in field 4 of parent entry
      INTEGER(LONG)                   :: RELID     = 0     ! This elements' ID


      REAL(DOUBLE)                    :: R8INP             ! A real value read from a field on this RBE3 entry
      REAL(DOUBLE)                    :: WGT               ! A weight read from a field on this RBE3 entry
      REAL(DOUBLE)                    :: WTi(MRBE3)        ! Array of RBE3 weight values
      REAL(DOUBLE)                    :: WT_TOT            ! Total of all WTi(i)



! **********************************************************************************************************************************
! RBE3 Bulk Data Card:

!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    2      ELID            RBE3 elem ID
!    3      blank
!    4      REFGRID         Ref grid point ID (this grid will go into M-set)
!    5      REFC            Ref component number (1-6) (comp of REFGRID)
!    6      WT1             Weighting factor for 1st component of displ
!    7      C1              Comp number for associated with WT1
!    8      G1,1            Grid associated with WT1
!    9      G1,2 or next WTi or blank

! on optional continuation entries:
!    2-9    contain more of the same: weights, components, grida


! Subsequent entries have the same format where a WTi is specified followed by a displ component Ci followed by a list of grids that
! have the WTi for that component

      JERR = 0

! Initialize arrays

      DO J=1,MRBE3
         GRID(J) = 0
         COMP(J) = 0
         WTi(J)  = ZERO
      ENDDO
      WT_TOT = ZERO

      TOKTYP(1:) = ' '

      IRBE3 = 0

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Up the elem count

      NRIGEL = NRIGEL+1

! Read and check data on parent entry

      CHAR_ELID = JCARD(2)
      CALL I4FLD ( JCARD(2), JF(2), I4INP )                ! Field 2: Elem ID
      IF (IERRFL(2) == 'N') THEN
         RELID                   = I4INP
         RIGID_ELEM_IDS(NRIGEL) = RELID
      ELSE
         JERR = JERR + 1
      ENDIF

      CALL I4FLD ( JCARD(4), JF(4), I4INP )                ! Field 4: REFGRID
      IF (IERRFL(4) == 'N') THEN
         REFGRID = I4INP
      ELSE
         JERR = JERR + 1
      ENDIF

      CALL I4FLD ( JCARD(5), JF(5), I4INP )                ! Field 5: REFC
      IF (IERRFL(5) == 'N') THEN
         CALL IP6CHK ( JCARD(5), CHRINP, IP6TYP, REFC_NUM_Ci )
         IF (IP6TYP(1:8) == 'COMP NOS') THEN
            REFC = I4INP
         ELSE
            JERR = JERR + 1
            WRITE(ERR,1124) JF(5), 'RBE3', CHAR_ELID, JF(5), JCARD(5)
            WRITE(F06,1124) JF(5), 'RBE3', CHAR_ELID, JF(5), JCARD(5)
         ENDIF
      ELSE
         JERR = JERR + 1
      ENDIF

      WGT = ZERO
      CALL R8FLD ( JCARD(6), JF(6), R8INP )                ! Field 6: 1st weight
      IF (IERRFL(6) == 'N') THEN
         WGT = R8INP
      ELSE
         JERR = JERR + 1
      ENDIF

      NUM_Ci = 0                                           ! Field 7: 1st field with Ci components
      INTDOF = 0
      CALL I4FLD ( JCARD(7), JF(7), I4INP )
      IF (IERRFL(7) == 'N') THEN
         CALL IP6CHK ( JCARD(7), CHRINP, IP6TYP, NUM_Ci )
         IF (IP6TYP(1:8) == 'COMP NOS') THEN
            INTDOF = I4INP
         ELSE
            JERR = JERR + 1
            WRITE(ERR,1124) JF(7), 'RBE3', CHAR_ELID, JF(7), JCARD(7)
            WRITE(F06,1124) JF(7), 'RBE3', CHAR_ELID, JF(7), JCARD(7)
         ENDIF
      ELSE
         JERR = JERR + 1
      ENDIF

      NEXT_MUST_BE_GRID = 'Y'
      GRID1 = 0                                            ! Field 8: 1st grid
      CALL I4FLD ( JCARD(8), JF(8), I4INP )
      IF (IERRFL(8) == 'N') THEN
         GRID1 = I4INP
      ELSE
         JERR = JERR + 1
      ENDIF
      NEXT_MUST_BE_GRID = 'N'

      IRBE3 = IRBE3 + 1                                    ! Write 1st set of wgt/comp/indep grid to arrays
      IF (IRBE3 > MRBE3) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, IRBE3, 'WTi' )
      WTi(IRBE3)  = WGT
      WT_TOT = WT_TOT + WTi(IRBE3)
      COMP(IRBE3) = INTDOF
      GRID(IRBE3) = GRID1
                                                           ! Field 9: if not blank then it must be a new weight or another grid
      FND_NEW_WGT = 'N'                                    ! that uses the previously read WT1 and C1
      TOKEN = JCARD(9)(1:8)                                ! Only send the 1st 8 chars of this JCARD. It has been left justified
      CALL TOKCHK ( TOKEN, TOKTYP )
      IF (TOKTYP /= 'BLANK   ') THEN
         IF      (TOKTYP == 'FL PT   ') THEN
            CALL R8FLD ( JCARD(9), JF(9), R8INP )
            IF (IERRFL(9) == 'N') THEN
               WGT         = R8INP
               FND_NEW_WGT = 'Y'
            ELSE
               JERR = JERR + 1
            ENDIF
         ELSE IF (TOKTYP == 'INTEGER ') THEN               ! Since it is integer, assume a grid and incr IRBE3 and write to arrays
            CALL I4FLD ( JCARD(9), JF(9), I4INP )
            IRBE3 = IRBE3 + 1
            GRID(IRBE3) = 0
            IF (IRBE3 > MRBE3) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, IRBE3, 'WTi' )
            WTi(IRBE3)  = WGT
            WT_TOT = WT_TOT + WTi(IRBE3)
            COMP(IRBE3) = INTDOF
            IF (IERRFL(8) == 'N') THEN
               GRID(IRBE3) = I4INP
            ELSE
               JERR = JERR + 1
            ENDIF
         ELSE
            JERR = JERR + 1
            WRITE(ERR,1151) '9', 'RBE3', CHAR_ELID, JCARD(9)
            WRITE(F06,1151) '9', 'RBE3', CHAR_ELID, JCARD(9)
         ENDIF
      ENDIF

      CALL BD_IMBEDDED_BLANK ( JCARD,2,0,4,0,6,0,8,9 )     ! Make sure there are no imbedded blanks (except 3,5,7)
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,3,0,0,0,0,0,0 )   ! Issue warning if field 3 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Read and check data on optional continuation entries

      I = 0
      DO

         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN

            I = I + 1

            DO J=2,9
               CHK(J) = 0
               CALL CHAR_FLD ( JCARD(J), JF(J), CHRINP )
               IF ((CHRINP(1:2) == 'UM') .OR. (CHRINP(1:4) == '"UM"')) THEN
                  JERR = JERR + 1
                  WRITE(ERR,1155) JF(J), 'RBE3', CHAR_ELID, I
                  WRITE(F06,1155) JF(J), 'RBE3', CHAR_ELID, I
               ENDIF
            ENDDO

do_j:       DO J=2,9

               IF (FND_NEW_WGT == 'Y') THEN                ! Last field read was FL PT so next field must be a Ci (displ components)
                  TOKEN = JCARD(J)(1:8)                    ! Only send the 1st 8 chars of this JCARD. It has been left justified
                  CALL TOKCHK ( TOKEN, TOKTYP )
                  INTDOF = 0
                  CALL I4FLD ( JCARD(J), JF(J), I4INP )
                  IF (IERRFL(J) == 'N') THEN
                     CALL IP6CHK ( JCARD(J), CHRINP, IP6TYP, NUM_Ci )
                     IF (IP6TYP(1:8) == 'COMP NOS') THEN
                        INTDOF = I4INP
                     ELSE
                        JERR = JERR + 1
                        WRITE(ERR,1129) JF(J), 'RBE3', CHAR_ELID, I, JF(J), JCARD(J)
                        WRITE(F06,1129) JF(J), 'RBE3', CHAR_ELID, I, JF(J), JCARD(J)
                     ENDIF
                  ELSE
                     JERR = JERR + 1
                  ENDIF
                  CHK(J) = J
                  FND_NEW_WGT = 'N'
                  NEXT_MUST_BE_GRID = 'Y'
                  CYCLE do_j
               ENDIF

               TOKEN = JCARD(J)(1:8)                       ! Only send the 1st 8 chars of this JCARD. It has been left justified
               CALL TOKCHK ( TOKEN, TOKTYP )               ! Last field read was not a weight so see what it is
               IF      (TOKTYP == 'BLANK   ') THEN
                  CYCLE do_j
               ELSE IF (TOKTYP == 'FL PT   ') THEN         ! Found a new weight
                  IF (NEXT_MUST_BE_GRID == 'N') THEN
                     WGT = ZERO
                     FND_NEW_WGT = 'Y'
                     CALL R8FLD ( JCARD(J), JF(J), R8INP )
                     IF (IERRFL(6) == 'N') THEN
                        WGT = R8INP
                     ELSE
                        JERR = JERR + 1
                     ENDIF
                  ELSE
                     JERR = JERR + 1
                     WRITE(ERR,1154) J, 'RBE3', CHAR_ELID, I
                     WRITE(F06,1154) J, 'RBE3', CHAR_ELID, I
                     NEXT_MUST_BE_GRID = 'N'
                  ENDIF
               ELSE                                        ! Not a weight or Ci so this must be a grid
                  GRID1 = 0
                  NEXT_MUST_BE_GRID = 'N'
                  CALL I4FLD ( JCARD(J), JF(J), I4INP )
                  IF (IERRFL(J) == 'N') THEN
                     GRID1 = I4INP
                     IRBE3 = IRBE3 + 1
                     IF (IRBE3 > MRBE3) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, IRBE3, 'WTi' )
                     WTi(IRBE3)  = WGT
                     WT_TOT = WT_TOT + WTi(IRBE3)
                     COMP(IRBE3) = INTDOF
                     GRID(IRBE3) = GRID1
                  ELSE
                     JERR = JERR + 1
                  ENDIF
               ENDIF
            ENDDO do_j
            CALL BD_IMBEDDED_BLANK ( JCARD, CHK(2), CHK(3), CHK(4), CHK(5), CHK(6), CHK(7), CHK(8), CHK(9)  )
            CALL CRDERR ( CARD )

         ELSE

            EXIT

         ENDIF

      ENDDO

! Write data to file L1F

      IF (JERR == 0) THEN
         WRITE(L1F) RTYPE
         NRECARD = NRECARD + 1
         WRITE(L1F) RELID, REFGRID, REFC, IRBE3, WT_TOT
         DO I=1,IRBE3
            WRITE(L1F) GRID(I), COMP(I), WTi(I)
            CALL RDOF ( COMP(I), CDOF )
            DO J=1,6
               IF (CDOF(J) /= '0') THEN
               ENDIF
            ENDDO
         ENDDO
      ELSE
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

!xx   NTERM_RMG = REFC_NUM_Ci*(NTERM_RMG + 1)



      RETURN

! **********************************************************************************************************************************
 1120 FORMAT(' *ERROR  1120: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY TRIPLETS OF GRID/COMPONENT/RBE3 COEFF ON BULK DATA RBE3 CARD WITH SET ID = ',I8              &
                    ,/,14X,' LIMIT IS = ',I12)

 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')


 1129 FORMAT(' *ERROR  1129: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ID = ',A,' CONTINUATION ENTRY',I4                         &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')


 1151 FORMAT(' *ERROR  1151: FIELD ',A,' ON ',A,' ID = ',A                                                                         &
                    ,/,14X,' MUST BE EITHER BLANK, A WEIGHT, OR A GRID NUMBER BUT IS: "',A,'"')

 1154 FORMAT(' *ERROR  1154: FIELD ',I3,' ON ',A,' ID = ',A,' CONT ENTRY',I4,' MUST BE A GRID PT SINCE THE',                       &
                           ' PREVIOUS 2 ENTRIES WERE A WEIGHT AND COMPS')

 1155 FORMAT(' *ERROR  1155: FIELD ',I3,' ON ',A,' ID = ',A,' CONT ENTRY',I4,' HAS "UM". THIS OPTION NOT PROGRAMMED IN MYSTRAN'    &
                    ,/,14X,' NOTIFY AUTHOR IF YOU WOULD LIKE THIS OPTION ADDED')





! **********************************************************************************************************************************

      END SUBROUTINE BD_RBE3


      SUBROUTINE BD_RBE30 ( CARD, LARGE_FLD_INP, IRBE3 )

! Processes RBE3 Bulk Data Cards to estimate the number of triplets of grid/comp/coeff on logical RBE3 B.D.card.
! The number of triplets defined on this RBE3 card will be returned to the calling routine so that the max number of triplets
! over all RBE3 cards can be determined. This conservative estimate will be based on counting all non blank fields except ones with
! a floating point number

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_RBE30'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! Type of TOKEN returned from subr TOKCHK

      INTEGER(LONG), INTENT(OUT)      :: IRBE3             ! Count of number of grid/comp/coeff triplets on this RBE3 logical card
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index




! **********************************************************************************************************************************
! RBE3 Bulk Data Card:

!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    2      ELID            RBE3 elem ID
!    3      blank
!    4      REFGRID         Ref grid point ID (this DOF will go into M-set)
!    5      REFC            Ref component number (1-6) (comp of REFGRID)
!    6      WT1             Weighting factor for 1st component of displ
!    7      C1              Comp number for associated with WT1
!    8      G1,1            Grid associated with WT1
!    9      G1,2 or next WT or blank

! on optional continuation entries:
!    2-9    contain either a WT, Ci, or blank (to be conservative say that every field that is not floating pt is a grid and
!           therefore represents another triplet of (grid/comp/weight)


! Subsequent entries have the same format where a WTi is specified followed by a displ component Ci followed by a list of grids that
! have the WTi for that component

! Max number of grids on parent entry can be 2

      IRBE3 = 2

! Count continuation entries

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
               TOKEN = JCARD(J)(1:8)                       ! Only send the 1st 8 chars of this JCARD. It has been left justified
               CALL TOKCHK ( TOKEN, TOKTYP )
               IF (TOKTYP /= 'FL PT   ') THEN
                  IRBE3 = IRBE3 + 1
               ENDIF
            ENDDO
         ELSE
            EXIT
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_RBE30


      SUBROUTINE BD_RSPLINE ( CARD, LARGE_FLD_INP )

! Processes RSPLINE Bulk Data Cards. Writes RSPLINE card data to file L1F. Note that this code only recognizes 2 indep grids on a
! given RSPLINE unlike NASTRAN (thus user must enter a different RSPLINE for each pair of indep grids between which the cubic
! spline will be fitted to all of the dep grid/comps

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1F
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, MRSPLINE, NRSPLINE, NRECARD, NRIGEL
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  RIGID_ELEM_IDS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, IP6CHK, R8FLD
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_RSPLINE'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: CHRINP            ! A character field returned from subr IP6CHK
      CHARACTER( 1*BYTE)              :: ERR_1158_WRITTEN  ! Indicator of whether fatal error #1158 has been written

      CHARACTER(LEN(JCARD))           :: FLDS_C(MRSPLINE+1)! Char array of all fields on this entry that have grid or comp vals.
!                                                            It will have IGRID1, then pairs of DGRID(i)/COMP(i), then IGRID2.

      CHARACTER(LEN(JCARD))           :: ID                ! Char value for element ID in field 2 of parent entry
      CHARACTER( 8*BYTE)              :: IP6TYP            ! Descriptor of what is in the 8 char field sent to subr IP6CHK
      CHARACTER(LEN(JCARD))           :: NAME              ! JCARD(1) from parent entry
      CHARACTER( 8*BYTE), PARAMETER   :: RTYPE = 'RSPLINE '! Rigid element type

      INTEGER(LONG)                   :: DCOMP(MRSPLINE+1) ! Array of dep displ comp values (1-6) found on this logical RSPLINE card
      INTEGER(LONG)                   :: DGRID(MRSPLINE+1) ! Array of dep GRID ID's found on this logical RSPLINE card

      INTEGER(LONG)                   :: FLDS_I(MRSPLINE+1)! Integer array of all fields on this entry that have grid or comp vals.
!                                                            It will have IGRID1, then pairs of DGRID(i)/COMP(i), then IGRID2.

      INTEGER(LONG)                   :: FLDS_F(MRSPLINE+1)! Field number on the RSPLINE entry where a value exists (parent or cont)
      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: IGRID1,IGRID2     ! The 2 independent grids on the RSPLINE entry
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      INTEGER(LONG)                   :: JE                ! An intermediale variable
      INTEGER(LONG)                   :: NUM_DEPENDENTS    ! Count of number of pairs of dependent grids/components
      INTEGER(LONG)                   :: NUM_ENTRIES       ! Count of number of entries placed into array GC_FLDS
      INTEGER(LONG)                   :: NUM_Ci            ! Number of displ components in a DCOMP field
      INTEGER(LONG)                   :: ELID      = 0     ! This elements' ID


      REAL(DOUBLE)                    :: DL_RAT            ! Value in field 3 for D/L ratio
      REAL(DOUBLE)                    :: R8INP             ! A real value read from a field on this RSPLINE entry



! **********************************************************************************************************************************
! RSPLINE Bulk Data Card:

!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    2      ELID            RSPLINE elem ID
!    3      D/L             Diam/length ratio of the rod used in determining the beam deflection relationship
!    4      G1              Indep grid at one end of the RSPLINE line
!    5-9    Gi,Ci           Dep grid/comps (dep grid/comps)

! on optional continuation entries:
!    2-9    contain either a grid or comp number. To be conservative each cont entry will be assumed to have an 8 dep entries

! Subsequent entries have the same format.

      JERR = 0

! Initialize arrays

      IGRID1 = 0
      IGRID2 = 0
      DO J=1,MRSPLINE+1
         DGRID(J)      = 0
         DCOMP(J)      = 0
         FLDS_I(J)     = 0
         FLDS_F(J)     = 0
         FLDS_C(J)(1:) = ' '
      ENDDO

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      NAME = JCARD(1)
      ID   = JCARD(2)

! Up the elem count

      NRSPLINE = NRSPLINE + 1
      NRIGEL   = NRIGEL+1

! Read and check data on parent entry

      CALL I4FLD ( JCARD(2), JF(2), I4INP )                ! Field 2: Elem ID
      IF (IERRFL(2) == 'N') THEN
         ELID                   = I4INP
         RIGID_ELEM_IDS(NRIGEL) = ELID
      ELSE
         JERR = JERR + 1
      ENDIF

      CALL R8FLD ( JCARD(3), JF(3), R8INP )                ! Field 3: DL_RAT
      IF (IERRFL(3) == 'N') THEN
         DL_RAT = R8INP
      ELSE
         JERR = JERR + 1
      ENDIF

      NUM_ENTRIES = 0                                      ! Fields 4-9
      DO J=4,9
         CALL I4FLD ( JCARD(J), JF(J), I4INP )
         IF (IERRFL(J) == 'N') THEN
            NUM_ENTRIES = NUM_ENTRIES + 1
            IF (NUM_ENTRIES > MRSPLINE + 1) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, NUM_ENTRIES, 'FLDS ARRAYS' )
            FLDS_I(NUM_ENTRIES) = I4INP
            FLDS_C(NUM_ENTRIES) = JCARD(J)
            FLDS_F(NUM_ENTRIES) = JF(J)
         ELSE
            JERR = JERR + 1
         ENDIF
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure there are no imbedded blanks
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Read and check data on optional continuation entries

do_1: DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         IF (ICONT == 1) THEN
            CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
            DO J=2,9
               CALL I4FLD ( JCARD(J), JF(J), I4INP )
               IF (IERRFL(J) == 'N') THEN
                  NUM_ENTRIES = NUM_ENTRIES + 1
                  IF (NUM_ENTRIES > MRSPLINE + 1) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, NUM_ENTRIES, 'FLDS ARRAYS' )
                  FLDS_I(NUM_ENTRIES) = I4INP
                  FLDS_C(NUM_ENTRIES) = JCARD(J)
                  FLDS_F(NUM_ENTRIES) = JF(J)
               ELSE
                  JERR = JERR + 1
               ENDIF
            ENDDO
            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
            CALL CRDERR ( CARD )
         ELSE
            EXIT do_1
         ENDIF
      ENDDO do_1

! FLDS_I may have some number of 0 entries at the end which should not be counted in NUM_ENTRIES.
! Scan FLDS_I (from end to begining) to find where last nonzero entry exists. This value is the one we want for NUM_ENTRIES

      JE = NUM_ENTRIES
      DO J=JE,1,-1
         IF (FLDS_I(J) == 0) THEN
            NUM_ENTRIES = NUM_ENTRIES - 1
            CYCLE
         ELSE
            EXIT
         ENDIF
      ENDDO

! Now that we have an accurate count for NUM_ENTRIES we can make some sanity checks on FLDS_I. NUM_ENTRIES should be an even
! number (IGRID1, IGRID2 and pairs of DGRID(i), DCOMP(i)). Also, FLDS_I should not have any zero values

      ERR_1158_WRITTEN = 'N'

      JE = MODULO ( NUM_ENTRIES, 2 )
      IF (JE /= 0) THEN                                    ! NUM_ENTRIES is not an even number unless JE = 0
         JERR = JERR + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1158) NAME, ID, 'THIS RSPLINE ENTRY DOES NOT HAVE AN EVEN NUMBER OF: DEP GRIDS + PAIRS OF DEP GRID/COMP ENTRIES'
         WRITE(F06,1158) NAME, ID, 'THIS RSPLINE ENTRY DOES NOT HAVE AN EVEN NUMBER OF: DEP GRIDS + PAIRS OF DEP GRID/COMP ENTRIES'
         ERR_1158_WRITTEN = 'Y'
      ENDIF

      DO J=1,NUM_ENTRIES
         IF (FLDS_I(J) == 0) THEN                         ! Zero vals in FLDS_I
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            IF (ERR_1158_WRITTEN == 'N') THEN
               WRITE(ERR,1158) NAME, ID, 'THIS RSPLINE ENTRY HAS SOME FIELDS THAT ARE BLANK OR THAT CONTAIN 0 VALUES'
               WRITE(F06,1158) NAME, ID, 'THIS RSPLINE ENTRY HAS SOME FIELDS THAT ARE BLANK OR THAT CONTAIN 0 VALUES'
            ELSE
               WRITE(ERR,11582)
               WRITE(F06,11582)
            ENDIF
         ENDIF
      ENDDO

! Make sure all DCOMP(i) are proper component numbers

      DO J=3,NUM_ENTRIES,2
         CALL IP6CHK ( FLDS_C(J), CHRINP, IP6TYP, NUM_Ci )
         IF (IP6TYP(1:8) /= 'COMP NOS') THEN
            JERR = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1124) FLDS_F(J), NAME, ID, FLDS_F(J), FLDS_C(J)
            WRITE(F06,1124) FLDS_F(J), NAME, ID, FLDS_F(J), FLDS_C(J)
         ENDIF
      ENDDO

! Determine the fields that will be written to L1F

      IGRID1 = FLDS_I(1)
      IF (NUM_ENTRIES > 0) THEN
         IGRID2 = FLDS_I(NUM_ENTRIES)
      ENDIF
      NUM_DEPENDENTS = 0
      DO J=2,NUM_ENTRIES-1,2
         NUM_DEPENDENTS = NUM_DEPENDENTS + 1
         DGRID(NUM_DEPENDENTS) = FLDS_I(J)
         DCOMP(NUM_DEPENDENTS) = FLDS_I(J+1)
      ENDDO


! If JERR /= 0 we can write data to L1F

      IF (JERR == 0) THEN
         DO J=1,NUM_DEPENDENTS
            WRITE(L1F) RTYPE
            WRITE(L1F) ELID, J, NUM_DEPENDENTS, IGRID1, IGRID2, DGRID(J), DCOMP(J), DL_RAT
            NRECARD = NRECARD + 1
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1124 FORMAT(' *ERROR  1124: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ENTRY WITH ID = ',A                                       &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')


 1129 FORMAT(' *ERROR  1129: INVALID DOF NUMBER IN FIELD ',I3,' ON ',A,' ID = ',A,' CONTINUATION ENTRY',I4                         &
                    ,/,14X,' MUST BE A COMBINATION OF DIGITS 1-6. HOWEVER, FIELD ',I3, ' HAS: "',A,'"')


 1158 FORMAT(' *ERROR  1158: FORMAT ERROR ON ',A,A,'. FORMAT MUST HAVE:'                                                           &
                    ,/,14X,' (a) TWO AND ONLY TWO INDEP GRIDS (IN FIELD 4 OF THE PARENT ENTRY AND IN THE LAST FIELD ON THE',       &
                           ' LOGICAL RSPLINE ENTRY)'                                                                               &
                    ,/,14X,' (b) PAIRS OF DEP GRIDS/COMPS FOLLOWING INDEP GRID 1 AND BEFORE INDEP GRID 2 WITH NO BLANK FIELDS'     &
                    ,/,15x,A)


11582 FORMAT('               THIS RSPLINE ENTRY ALSO HAS SOME FIELDS THAT ARE BLANK OR THAT CONTAIN 0 VALUES')

! **********************************************************************************************************************************

      END SUBROUTINE BD_RSPLINE


      SUBROUTINE BD_RSPLINE0 ( CARD, LARGE_FLD_INP, IRSPLINE )

! Processes RSPLINE Bulk Data Cards to count the number of fields that can have a grid or comp number.
! The first and last entries will be the 2 indep grids. Everything in between should be pairs of dep grid/comp entries

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_RSPLINE0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein

      INTEGER(LONG), INTENT(OUT)      :: IRSPLINE          ! Count of number of grid/comp doublets on this RSPLINE logical card
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! RSPLINE Bulk Data Card:

!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    2      ELID            RSPLINE elem ID
!    3      D/L             Diam/length ratio of the rod used in determining the beam deflection relationship
!    4      G1              Indep grid at one end of the RSPLINE line
!    5-9    Gi,Ci           Dep grid/comps (dep grid/comps)

! on optional continuation entries:
!    2-9    contain either a grid or comp number. To be conservative each cont entry will be assumed to have an 8 dep entries

! Subsequent entries have the same format

! Max number of grids or comps on parent entry can be 6

      IRSPLINE = 6

! Count continuation entries. Each one can have up to 8 entries

      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC0  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         IF (ICONT == 1) THEN
            IRSPLINE = IRSPLINE + 8
         ELSE
            EXIT
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_RSPLINE0

   END MODULE RIGID_ELEMENTS
