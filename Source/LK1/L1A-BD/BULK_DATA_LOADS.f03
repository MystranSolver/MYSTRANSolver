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

   MODULE BULK_DATA_LOADS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BD_LOAD, BD_LOAD0, BD_FORMOM, BD_GRAV, BD_PLOAD2, BD_PLOAD4, BD_RFORCE, BD_SLOAD, BD_SLOAD0

   CONTAINS

      SUBROUTINE BD_LOAD ( CARD, LARGE_FLD_INP, CC_LOAD_FND )

! Processes LOAD Bulk Data Cards. Reads and checks data and enters data into arrays LOAD_SIDS and LOAD_FACS

!  1) Load set ID's from the LOAD Bulk Data card are entered into array LOAD_SIDS
!  2) Scale factors (overall for this LOAD card and individual for each set ID) are entered into array LOAD_FACS

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LLOADR, LSUB, NLOAD, LLOADC, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  LOAD_SIDS, LOAD_FACS, SUBLOD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_LOAD'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD                ! A Bulk Data card
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_LOAD_FND(LSUB,2) ! 'Y' if B.D load/temp card w/ same set ID (SID) as C.C. LOAD = SID
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)           ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: NAME                ! JCARD(1) from parent entry

      INTEGER(LONG)                   :: I,J,K               ! DO loop index
      INTEGER(LONG)                   :: NUM_PAIRS           ! Counter on number of pairs of set ID's and scale factors on LOAD card
      INTEGER(LONG)                   :: ICONT     = 0       ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0       ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: SETID               ! Set ID for this LOAD Bulk Data card




! **********************************************************************************************************************************
! LOAD Bulk Data Card routine

!   FIELD   ITEM            ARRAY ELEMENT
!   -----   ------------    -------------
!    2      SID              LOAD_SIDS(nload,1)
!    3      S0               LOAD_FACS(nload,1)
!    4      S1               LOAD_FACS(nload,2)
!    5      L1               LOAD_SIDS(nload,2)
!    6      S2               LOAD_FACS(nload,3)
!    7      L2               LOAD_SIDS(nload,3)
!    8      S3               LOAD_FACS(nload,4)
!    9      L3               LOAD_SIDS(nload,4)

! Optiona continuation cards:
!
!   FIELD   ITEM            ARRAY ELEMENT
!   -----   ------------    -------------
!    2      S4               LOAD_FACS(nload,5)
!    3      L4               LOAD_SIDS(nload,5)
!    4      S5               LOAD_FACS(nload,6)
!    5      L5               LOAD_SIDS(nload,6)
!    6      S6               LOAD_FACS(nload,7)
!    7      L6               LOAD_SIDS(nload,7)
!    8      S7               LOAD_FACS(nload,8)
!    9      L7               LOAD_SIDS(nload,8)

! Subsequent con't cards follow the same patterm as the 1st


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      NAME = JCARD(1)

! Check for overflow

      NLOAD = NLOAD+1
!xx   IF (NLOAD > LLOADR) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LLOADR
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LLOADR
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data on parent card

      CALL I4FLD ( JCARD(2), JF(2), LOAD_SIDS(NLOAD,1) )   ! Read set ID for this LOAD Bulk Data card
      IF (IERRFL(2) == 'N') THEN
         SETID = LOAD_SIDS(NLOAD,1)
         DO I=1,NSUB
            IF (SETID == SUBLOD(I,1)) THEN
               CC_LOAD_FND(I,1) = 'Y'
            ENDIF
         ENDDO
      ENDIF

      CALL R8FLD ( JCARD(3), JF(3), LOAD_FACS(NLOAD,1) )   ! Read overall scale factor on LOAD card

      NUM_PAIRS = 1                                        ! Read pairs of load mags and load ID's on parent card.
!                                                          ! 1st "pair" is for the load set ID on this LOAD B.D. entry
      DO J=4,8,2
         IF ((JCARD(J)(1:) == ' ') .AND. (JCARD(J+1)(1:) == ' ')) THEN
            CYCLE                                          ! CYCLE if a pair is blank
         ELSE
            NUM_PAIRS = NUM_PAIRS + 1
            IF (NUM_PAIRS > LLOADC) THEN
               WRITE(ERR,1139) SETID, NAME, NAME, LLOADC
               WRITE(F06,1139) SETID, NAME, NAME, LLOADC
               CALL OUTA_HERE ( 'Y' )                       ! Coding error, so quit
            ENDIF
            CALL R8FLD ( JCARD(J)  , JF(J)  , LOAD_FACS(NLOAD,NUM_PAIRS) )
            CALL I4FLD ( JCARD(J+1), JF(J+1), LOAD_SIDS(NLOAD,NUM_PAIRS) )
            DO K=1,NUM_PAIRS-1                             ! Check for duplicate set ID's
               IF (LOAD_SIDS(NLOAD,NUM_PAIRS) == LOAD_SIDS(NLOAD,K)) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1140) LOAD_SIDS(NLOAD,NUM_PAIRS),NAME,SETID
                  WRITE(F06,1140) LOAD_SIDS(NLOAD,NUM_PAIRS),NAME,SETID
               ENDIF
            ENDDO
            IF (IERRFL(J+1) == 'N') THEN
               IF (LOAD_SIDS(NLOAD,NUM_PAIRS) <= 0) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1176) SETID,J+1,LOAD_SIDS(NLOAD,NUM_PAIRS)
                  WRITE(F06,1176) SETID,J+1,LOAD_SIDS(NLOAD,NUM_PAIRS)
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
            DO J=2,8,2
               IF ((JCARD(J)(1:) == ' ') .AND. (JCARD(J+1)(1:) == ' ')) THEN
                  CYCLE
               ELSE
                  NUM_PAIRS = NUM_PAIRS + 1
                  IF (NUM_PAIRS > LLOADC) THEN
                     WRITE(ERR,1139) SETID, NAME, NAME, LLOADC
                     WRITE(F06,1139) SETID, NAME, NAME, LLOADC
                     CALL OUTA_HERE ( 'Y' )                 ! Coding error, so quit
                  ENDIF
                  CALL R8FLD ( JCARD(J)  , JF(J)  , LOAD_FACS(NLOAD,NUM_PAIRS) )
                  CALL I4FLD ( JCARD(J+1), JF(J+1), LOAD_SIDS(NLOAD,NUM_PAIRS) )
                  DO K=1,NUM_PAIRS-1                       ! Check for duplicate set ID's
                     IF (LOAD_SIDS(NLOAD,NUM_PAIRS) == LOAD_SIDS(NLOAD,K)) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1140) LOAD_SIDS(NLOAD,NUM_PAIRS),NAME,SETID
                        WRITE(F06,1140) LOAD_SIDS(NLOAD,NUM_PAIRS),NAME,SETID
                     ENDIF
                  ENDDO
                  IF (IERRFL(J+1) == 'N') THEN
                     IF (LOAD_SIDS(NLOAD,NUM_PAIRS) <= 0) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1176) SETID,J+1,LOAD_SIDS(NLOAD,NUM_PAIRS)
                        WRITE(F06,1176) SETID,J+1,LOAD_SIDS(NLOAD,NUM_PAIRS)
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

 1176 FORMAT(' *ERROR  1176: LOAD SET ID ON BULK DATA LOAD ENTRY ID = ',I8,' IN FIELD ',I3,' IS = ',I8,'. MUST BE > 0')

! **********************************************************************************************************************************

      END SUBROUTINE BD_LOAD


      SUBROUTINE BD_LOAD0 ( CARD, LARGE_FLD_INP, ILOAD )

! Processes LOAD Bulk Data Cards to determine the number of pairs of load SID's and magnitudes on the LOAD B.D. card. Blank entries
! for pairs are not counted (specifically, if the load SID is 0 the pair is not counted).
! The number of pairs defined on this LOAD card will be returned to the calling routine so that the max number of pairs
! over all LOAD cards can be determined.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_LOAD0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG), INTENT(OUT)      :: ILOAD             ! Count of no. real load factors on this card. Starts with 1
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index




! **********************************************************************************************************************************
! LOAD Bulk Data Card:

!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    2      SID             LOAD set ID
!    3      S0              Overall scale factor
!    4      S1              Scale factor for set whose ID is L1
!    5      L1              Set ID of a FORCE, MOMENT or GRAV load
!    6      S2              Scale factor for set whose ID is L2
!    7      L2              Set ID of a FORCE, MOMENT or GRAV load
!    8      S3              Scale factor for set whose ID is L3
!    9      L3              Set ID of a FORCE, MOMENT or GRAV load

! 1st continuation card:
!
!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    2      S4              Scale factor for set whose ID is L4
!    3      L4              Set ID of a FORCE, MOMENT or GRAV load
!    4      S5              Scale factor for set whose ID is L5
!    5      L5              Set ID of a FORCE, MOMENT or GRAV load
!    6      S6              Scale factor for set whose ID is L6
!    7      L6              Set ID of a FORCE, MOMENT or GRAV load
!    8      S7              Scale factor for set whose ID is L7
!    9      L7              Set ID of a FORCE, MOMENT or GRAV load

! Subsequent continuation cards follow the same pattern as the 1st continuation card


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! ILOAD will count the number of real load factors on this LOAD card. It starts with 1 since there is an overall factor on the LOAD
! card, and continues as factors for the individual loads defined on the card are read

! Count pairs of load ID's/factors on parent card

      ILOAD = 1
      DO J=5,9,2
         IF ((JCARD(J)(1:) == ' ') .AND. (JCARD(J+1)(1:) == ' ')) THEN
            CYCLE
         ELSE
            ILOAD = ILOAD + 1
         ENDIF
      ENDDO

! Count pairs of load ID's/factors on optional continuation cards

      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC0  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN
            DO J=2,8,2
               IF ((JCARD(J)(1:) == ' ') .AND. (JCARD(J+1)(1:) == ' ')) THEN
                  CYCLE
               ELSE
                  ILOAD = ILOAD + 1
               ENDIF
            ENDDO
         ELSE
            EXIT
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_LOAD0


      SUBROUTINE BD_FORMOM ( CARD, CC_LOAD_FND )

! Processes FORCE, MOMENT Bulk Data Cards. Also, the set ID is written to array FORMOM_SIDS which is checked
! in subroutine LOADB to make sure that all set ID's requested in Case Control were found in the Bulk Data.
! A record is written to file LINK1I for each FORCE or MOMENT Bulk Data card with the following data:

!   SETID, GRID_NO, CID, FORMON1, FORMON2, FORMON3, FOR_OR_MOM

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1I
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, LFORCE, LSUB, NFORCE, NSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL, SUPWARN
      USE MODEL_STUF, ONLY            :  FORMOM_SIDS, SUBLOD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_FORMOM'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD                ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)           ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: FOR_OR_MOM          ! = 'FORCE' or 'MOMENT' (JCARD(1) from parent entry)
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_LOAD_FND(LSUB,2) ! 'Y' if B.D load/temp card w/ same set ID (SID) as C.C. LOAD = SID

      INTEGER(LONG)                   :: CID       = 0       ! Coord ID on the FORCE/MOMENT card
      INTEGER(LONG)                   :: GRID_NO   = 0       ! Grid ID  on the FORCE/MOMENT card
      INTEGER(LONG)                   :: I                   ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0       ! A local error count
      INTEGER(LONG)                   :: SETID  = 0          ! Set ID on the FORCE/MOMENT card


      REAL(DOUBLE)                    :: EPS1                ! A small number to compare real zero
      REAL(DOUBLE)                    :: FORMON1   = ZERO    ! Force/moment magnitude for 1st dir in coord sys CID (= SCALEF*V1)
      REAL(DOUBLE)                    :: FORMON2   = ZERO    ! Force/moment magnitude for 2nd dir in coord sys CID (= SCALEF*V2)
      REAL(DOUBLE)                    :: FORMON3   = ZERO    ! Force/moment magnitude for 3rd dir in coord sys CID (= SCALEF*V3)
      REAL(DOUBLE)                    :: SCALEF    = ZERO    ! Scale factor on the FORCE/MOMENT card
      REAL(DOUBLE)                    :: V1        = ZERO    ! Amplitude of force/mom on FORCE/MOMENT card in 1st direction of CID
      REAL(DOUBLE)                    :: V2        = ZERO    ! Amplitude of force/mom on FORCE/MOMENT card in 2nd direction of CID
      REAL(DOUBLE)                    :: V3        = ZERO    ! Amplitude of force/mom on FORCE/MOMENT card in 3rd direction of CID
      REAL(DOUBLE)                    :: VMAG      = ZERO    ! V!**2 + V2**2 + V3**2

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! FORCE / MOMENT Bulk Data Card routine

!   FIELD   ITEM           VARIABLE
!   -----   ------------   -------------
!    2      Load set ID    FORMOM_SIDS(I), SETID
!    3      Grid ID        GRID_NO
!    4      Cord. ID       CID
!    5      Scale factor   SCALEF
!    6      Vect. comp. 1  V1
!    7      Vect. comp. 2  V2
!    8      Vect. comp. 3  V3

!    SCALEF, V1, V2, and V3 are processed to force component
!    form ( Fx,Fy,Fz or Mx,My,Mz ) and stored in FORMON1-3

      EPS1 = EPSIL(1)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      IF      (JCARD(1)(1:5) == 'FORCE' ) THEN
         FOR_OR_MOM = 'FORCE   '
      ELSE IF (JCARD(1)(1:6) == 'MOMENT') THEN
         FOR_OR_MOM = 'MOMENT  '
      ENDIF

! Check for overflow

      NFORCE = NFORCE+1
!xx   IF (NFORCE > LFORCE) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LFORCE
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LFORCE
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data

      CALL I4FLD ( JCARD(2), JF(2), SETID )
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NSUB
            IF (SETID == SUBLOD(I,1)) THEN
               CC_LOAD_FND(I,1) = 'Y'
            ENDIF
         ENDDO
         FORMOM_SIDS(NFORCE) = SETID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), GRID_NO )              ! Read grid that force is at
      CALL I4FLD ( JCARD(4), JF(4), CID )                  ! Read coord system force is described in
      CALL R8FLD ( JCARD(5), JF(5), SCALEF )               ! Read force scale factor
      CALL R8FLD ( JCARD(6), JF(6), V1 )                   ! Read magnitude in 1st direction of coord sys CID
      CALL R8FLD ( JCARD(7), JF(7), V2 )                   ! Read magnitude in 2nd direction of coord sys CID
      CALL R8FLD ( JCARD(8), JF(8), V3 )                   ! Read magnitude in 3rd direction of coord sys CID

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,0 )     ! Make sure that there are no imbedded blanks in fields 2-8
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )   ! Issue warning if field 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      DO I=2,8                                             ! Set JERR if any errors reading above data
         IF (IERRFL(I) == 'Y') THEN
            JERR = JERR + 1
         ENDIF
      ENDDO

! Write data to file LINK1I if there were no errors

      IF (JERR == 0) THEN
         VMAG = SQRT(V1*V1+V2*V2+V3*V3)
         IF (DABS(VMAG) < EPS1) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) CARD
            WRITE(ERR,1137) JCARD(1), JCARD(2)
            IF (SUPWARN == 'N') THEN
               IF (ECHO == 'NONE  ') THEN
                  WRITE(F06,101) CARD
               ENDIF
               WRITE(F06,1137) JCARD(1), JCARD(2)
            ENDIF
            FORMON1 = ZERO
            FORMON2 = ZERO
            FORMON3 = ZERO
         ELSE
            FORMON1 = SCALEF*V1
            FORMON2 = SCALEF*V2
            FORMON3 = SCALEF*V3
         ENDIF
         WRITE(L1I) SETID, GRID_NO, CID, FORMON1, FORMON2, FORMON3, FOR_OR_MOM
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1137 FORMAT(' *WARNING    : ',A,' ENTRY WITH SET ID = ',A,' HAS ZERO COMPONENTS')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_FORMOM


      SUBROUTINE BD_GRAV ( CARD, LARGE_FLD_INP, CC_LOAD_FND )

! Processes GRAV Bulk Data Cards.  Also, the set ID is written to array GRAV_SIDS which is checked in subroutine
! LOADB to make sure that all set ID's requested in Case Control were found in the Bulk Data.
! A record is written to file LINK1P for each GRAV Bulk Data card with the following data:

!   SETID, CID, ACCEL(1-6)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1P
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LGRAV, LSUB, NGRAV, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  GRAV_SIDS, SUBLOD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_GRAV'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD               ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)          ! The 10 fields of 8 characters making up CARD
      CHARACTER(LEN(JCARD))           :: NAME               ! JCARD(1) from parent entry
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_LOAD_FND(LSUB,2)! 'Y' if B.D load/temp card w/ same set ID (SID) as C.C. LOAD = SID

      INTEGER(LONG)                   :: CID        = 0     ! Coord ID on the GRAV card
      INTEGER(LONG)                   :: CONT_COUNT = 0     ! Count of number of continuation entries
      INTEGER(LONG)                   :: I                  ! DO loop index
      INTEGER(LONG)                   :: I4INP              ! An integer value read from GRAV entry
      INTEGER(LONG)                   :: ICONT      = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR       = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: GID        = 0     ! Grid ID (or 0) of the grid that the rotational grav accels refer to
      INTEGER(LONG)                   :: JERR       = 0     ! A local error count
      INTEGER(LONG)                   :: SETID      = 0     ! Set ID on the GRAV card


      REAL(DOUBLE)                    :: ACCEL(6)           ! Gravity magnitudes in the 3 translational and 3 rotational dirs
      REAL(DOUBLE)                    :: SCALEF     = ZERO  ! Scale factor on the GRAV card
      REAL(DOUBLE)                    :: VEC(6)             ! Vector components of gravity read from the GRAV entry

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! GRAV Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Load set ID    GRAV_SIDS(ngrav)
!    3      Cord. ID       CID
!    4      Grav - G       SCALEF
!    5-7    Vector comps   VEC(1-3)

! on optional second card:
!    2      GID
!    3-5    Vector comps   VEC(4-6)

! Initialize variables

      DO I=1,6
         ACCEL(I) = ZERO
         VEC(I)   = ZERO
      ENDDO

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      NAME = JCARD(1)

! Check for overflow

      NGRAV = NGRAV+1
!xx   IF (NGRAV > LGRAV) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LGRAV
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LGRAV
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data

      CALL I4FLD ( JCARD(2), JF(2), SETID )
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NSUB
            IF (SETID == SUBLOD(I,1)) THEN
               CC_LOAD_FND(I,1) = 'Y'
            ENDIF
         ENDDO
         GRAV_SIDS(NGRAV) = SETID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), CID )
      CALL R8FLD ( JCARD(4), JF(4), SCALEF )
      CALL R8FLD ( JCARD(5), JF(5), VEC(1) )
      CALL R8FLD ( JCARD(6), JF(6), VEC(2) )
      CALL R8FLD ( JCARD(7), JF(7), VEC(3) )

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,0,0 )     ! Make sure that there are no imbedded blanks in fields 2-7
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )   ! Issue warning if fields 8,9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      DO I=2,7
         IF (IERRFL(I) == 'Y') THEN
            JERR = JERR + 1
         ENDIF
      ENDDO

! Optional Second Card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         CONT_COUNT = 1
         CALL I4FLD ( JCARD(2), JF(2), I4INP )             ! Read grid ID
         IF (IERRFL(2) == 'N') THEN
            IF (I4INP >= 0) THEN
               GID = I4INP
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1194) NAME, JF(2), CONT_COUNT, JCARD(2)
               WRITE(F06,1194) NAME, JF(2), CONT_COUNT, JCARD(2)
            ENDIF
         ENDIF

         CALL R8FLD ( JCARD(3), JF(3), VEC(4) )
         CALL R8FLD ( JCARD(4), JF(4), VEC(5) )
         CALL R8FLD ( JCARD(5), JF(5), VEC(6) )
                                                           ! Get offsets, if present
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,0,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-5
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )! Issue warning if fields 6-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

         DO I=2,5
            IF (IERRFL(I) == 'Y') THEN
               JERR = JERR + 1
            ENDIF
         ENDDO

      ENDIF

! Write data to file LINK1P if there were no errors

      IF (JERR == 0) THEN
         DO I=1,6
            ACCEL(I) = SCALEF*VEC(I)
         ENDDO
         WRITE(L1P) SETID, CID, GID, (ACCEL(I),I=1,6)
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1194 FORMAT(' *ERROR  1194: NEGATIVE GRID ID NOT ALLOWED ON ',A,' ENTRY. VALUE IN FIELD ',I3,' OF CONTINUATION ENTRY NUMBER '     &
                            ,I3,' IS = ',A)

 ! *********************************************************************************************************************************

      END SUBROUTINE BD_GRAV


      SUBROUTINE BD_PLOAD2 ( CARD, CC_LOAD_FND )

! Processes PLOAD2 Bulk Data Cards. Reads and checks data and then writes CARD to file LINK1Q for later processing

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1Q
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPLOAD, LSUB, NPCARD, NPLOAD, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PRESS_SIDS, SUBLOD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PLOAD2'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD               ! A Bulk Data card
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_LOAD_FND(LSUB,2)! 'Y' if B.D load/temp card w/ same set ID (SID) as C.C. LOAD = SID
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)          ! The 10 fields of characters making up CARD
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

! Read pressure value

      CALL R8FLD ( JCARD(3), JF(3), RPRESS )

! Check for the 2 options on specifying data on the card. Either all data are PLOAD_ELID's or the THRU  option is used in
! which case field 5 will have "THRU".


      IF      ((JCARD(1)(1:7) == 'PLOAD1 ') .OR. (JCARD(1)(1:7) == 'PLOAD1*')) THEN
         WRITE(ERR,99)
         WRITE(F06,99)
         FATAL_ERR = FATAL_ERR + 1
      ELSE IF ((JCARD(1)(1:7) == 'PLOAD2 ') .OR. (JCARD(1)(1:7) == 'PLOAD2*')) THEN

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
   99 FORMAT(' *ERROR      : CODE NOT WRITTEN YET FOR PLOAD1')

 1128 FORMAT(' *ERROR  1128: ON ',A,A,' THE IDs MUST BE IN INCREASING ORDER FOR THRU OPTION')

 1152 FORMAT(' *ERROR  1152: ON ',A,A,' ELEM IDs MUST BE > 0')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PLOAD2


      SUBROUTINE BD_PLOAD4 ( CARD, CC_LOAD_FND )

! Processes PLOAD4 Bulk Data Cards. Reads and checks data and then writes CARD to file LINK1Q for later processing

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1Q
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPLOAD, LSUB, NPCARD, NPLOAD,             &
                                         NPLOAD4_3D, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PRESS_SIDS, SUBLOD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, LEFT_ADJ_BDFLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PLOAD4'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD               ! A Bulk Data card
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_LOAD_FND(LSUB,2)! 'Y' if B.D load/temp card w/ same set ID (SID) as C.C. LOAD = SID
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)          ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: ELID1,ELID2        ! Elem ID's on parent card. If "THRU" not in field 8, ELID2 is no present
      INTEGER(LONG)                   :: I4INP              ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: J                  ! DO loop index
      INTEGER(LONG)                   :: JERR               ! Error count
      INTEGER(LONG)                   :: SETID              ! Load set ID on PLOADi card


      REAL(DOUBLE)                    :: R8INP              ! A value read from input file that should be a real value



! **********************************************************************************************************************************
! PLOAD4 Bulk Data card check.
! Note that if both fields 8 and 9 are blank this is for a plate elem (and this does not need to be verified)

! Format 1:
! --------
!   FIELD   ITEM          Description
!   -----   ----       ------------------
!    2      SID        Load set ID
!    3      ELID       Element ID
!    4      P1         Pressure at grid 1
!    5      P2         Pressure at grid 2
!    6      P3         Pressure at grid 3
!    7      P4         Pressure at grid 4 (not used unless elem is a CHEXA, CQUAD)
!    8      G1         ID of grid connected to a corner of the face (3D elems only)
!    9      G3 or G4   G3 is the ID of a grid connected to a corner diagonally opposite to G1 on the same face of a CHEXA
!                      G4 is the ID of a CTETRA grid located at the corner (not on the face being loaded). CTETRA only

! Format 2: (QUAD4, TRIA3 elems only)
! --------
!   FIELD   ITEM          Description
!   -----   ----       ------------------
!    2      SID        Load set ID
!    3      ELID1      Element ID of 1st elem
!    4      P1         Pressure at grid 1
!    5      P2         Pressure at grid 2
!    6      P3         Pressure at grid 3
!    7      P4         Pressure at grid 4 (not used unless elem is a QUAD4)
!    8      "THRU"
!    9      ELID2      Element ID of 1st elem. Elems in the range ELID1 through ELID2 will be loaded.


! Make JCARD from CARD and up the count on NPLOAD

      JERR = 0
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

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


! Read data in fields 3-7 (same for either format). Only need to make sure data is correct format. We don't need values here

      ELID1 = 0                                             ! Element ID
      CALL I4FLD ( JCARD(3), JF(3), I4INP )
      IF (IERRFL(3) == 'N') THEN
         ELID1 = I4INP
      ENDIF

      IF (JCARD(4)(1:) /= ' ') THEN
         CALL R8FLD ( JCARD(4), JF(4), R8INP )
      ELSE
         JERR      = JERR + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1198) JF(4), JCARD(1), JCARD(2)
         WRITE(F06,1198) JF(4), JCARD(1), JCARD(2)
      ENDIF

      IF (JCARD(5)(1:) /= ' ') THEN                        ! Pressure at grid 2 and/or default to pressure at grid 1
         CALL R8FLD ( JCARD(5), JF(5), R8INP )
      ENDIF

      IF (JCARD(6)(1:) /= ' ') THEN                        ! Pressure at grid 3 and/or default to pressure at grid 1
         CALL R8FLD ( JCARD(6), JF(6), R8INP )
      ENDIF

      IF (JCARD(7)(1:) /= ' ') THEN                        ! Pressure at grid 4 and/or default to pressure at grid 1
         CALL R8FLD ( JCARD(7), JF(7), R8INP )
      ENDIF

! Format 2:
! --------
      CALL LEFT_ADJ_BDFLD ( JCARD(8) )
      IF (JCARD(8)(1:4) == 'THRU') THEN

         ELID2 = 0
         CALL I4FLD ( JCARD(9), JF(9), I4INP )
         IF (IERRFL(9) == 'N') THEN
            ELID2 = I4INP
            IF (ELID2 >= ELID1) THEN
               ELID2 = I4INP
            ELSE
               JERR      = JERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1128) JCARD(1), JCARD(2)
               WRITE(F06,1128) JCARD(1), JCARD(2)
            ENDIF
         ENDIF

! Format 1:
! --------                                                 ! If flds 8,9 are not blank and fld 8 /= THRU this must be for a 3D elem
      ELSE IF ((JCARD(8)(1:) /= ' ') .AND. (JCARD(9)(1:) /= ' ')) THEN

         jerr = jerr + 1
         fatal_err = fatal_err + 1
         Write(err,99)
         Write(f06,99)
         return

         CALL I4FLD ( JCARD(8), JF(8), I4INP )
         CALL I4FLD ( JCARD(9), JF(9), I4INP )
         NPLOAD4_3D = NPLOAD4_3D + 1                       ! Increment count of number of solid elems that have PLOAD4 pressure def

      ENDIF

! Check for card error and fields blank

      IF (JCARD(8)(1:4) == 'THRU') THEN
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-7
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )! Issue warning if fieldS 8, 9 not blank
      ELSE
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
      ENDIF
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Write data to file L1Q

      IF (JERR == 0) THEN
         WRITE(L1Q) CARD
      ENDIF

      NPCARD = NPCARD + 1



      RETURN

! **********************************************************************************************************************************
   99 FORMAT(' *ERROR      : CODE NOT WRITTEN YET FOR PLOAD4 OPTION USING G1 AND G3/G4')

 1128 FORMAT(' *ERROR  1128: ON ',A,A,' THE IDs MUST BE IN INCREASING ORDER FOR THRU OPTION')

 1152 FORMAT(' *ERROR  1152: ON ',A,A,' ELEM IDs MUST BE > 0')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1198 FORMAT(' *ERROR  1198: FIELD ',I2,' ON ',A,A,' CANNOT BE BLANK')

! **********************************************************************************************************************************

      END SUBROUTINE BD_PLOAD4


      SUBROUTINE BD_RFORCE ( CARD, LARGE_FLD_INP, CC_LOAD_FND )

! Processes RFORCE Bulk Data Cards.  Also, the set ID is written to array RFORCE_SIDS which is checked in subroutine
! LOADB to make sure that all set ID's requested in Case Control were found in the Bulk Data.
! A record is written to file LINK1U for each RFORCE Bulk Data card with the following data:

!   SETID, CID, ACCEL(1-6)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1U
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LRFORCE, LSUB, NRFORCE, NSUB
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  RFORCE_SIDS, SUBLOD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_RFORCE'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD               ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)          ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: NAME          ! The character name of this Bulk Data entry
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_LOAD_FND(LSUB,2)! 'Y' if B.D load/temp card w/ same set ID (SID) as C.C. LOAD = SID

      INTEGER(LONG)                   :: CID        = 0     ! Coord ID on the RFORCE card
      INTEGER(LONG)                   :: CONT_COUNT = 0     ! Count of number of continuation entries
      INTEGER(LONG)                   :: I                  ! DO loop index
      INTEGER(LONG)                   :: ICONT      = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR       = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: GID        = 0     ! Grid ID (or 0) of the grid that the rotational grav accels refer to
      INTEGER(LONG)                   :: JERR       = 0     ! A local error count
      INTEGER(LONG)                   :: SETID      = 0     ! Set ID on the RFORCE card


      REAL(DOUBLE)                    :: R8INP              ! A real value read from RFORCE entry
      REAL(DOUBLE)                    :: SCALEF_AA  = ZERO  ! Scale factor for angular accel    on the RFORCE card
      REAL(DOUBLE)                    :: SCALEF_AV  = ZERO  ! Scale factor for angular velocity on the RFORCE card
      REAL(DOUBLE)                    :: VEC(3)             ! Vector components of the angular vel and/or accel on the RFORCE entry

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! RFORCE Bulk Data Card routine

!   FIELD   ITEM                                        ARRAY ELEMENT
!   -----   ------------                                -------------
!    2      Load set ID                                 RFORCE_SIDS(nrforce)
!    3      Grid pt about which the rotation occurs     GID
!    4      Cord. ID that rotation vector is defined in CID
!    5      Scale factor for angular velocity           SCALEF_AV
!    6-8    Vector comps                                VEC(1-3)

! on optional second card:
!    2      Scale factor for angular acceleration       SCALEF_AA

! Initialize variables

      DO I=1,3
         VEC(I)   = ZERO
      ENDDO

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      NAME = JCARD(1)

! Check for overflow

      NRFORCE = NRFORCE+1
!xx   IF (NRFORCE > LRFORCE) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LRFORCE
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LRFORCE
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data

      CALL I4FLD ( JCARD(2), JF(2), SETID )
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NSUB
            IF (SETID == SUBLOD(I,1)) THEN
               CC_LOAD_FND(I,1) = 'Y'
            ENDIF
         ENDDO
         RFORCE_SIDS(NRFORCE) = SETID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3), GID )
      CALL I4FLD ( JCARD(4), JF(4), CID )
      CALL R8FLD ( JCARD(5), JF(5), SCALEF_AV )
      CALL R8FLD ( JCARD(6), JF(6), VEC(1) )
      CALL R8FLD ( JCARD(7), JF(7), VEC(2) )
      CALL R8FLD ( JCARD(8), JF(8), VEC(3) )

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,0 )     ! Make sure that there are no imbedded blanks in fields 2-8
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )   ! Issue warning if field 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      DO I=2,8
         IF (IERRFL(I) == 'Y') THEN
            JERR = JERR + 1
         ENDIF
      ENDDO

! Optional Second Card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         CONT_COUNT = 1
         CALL R8FLD ( JCARD(2), JF(2), R8INP )             ! Read grid ID
         IF (IERRFL(2) == 'N') THEN
            SCALEF_AA = R8INP
         ENDIF

         CALL BD_IMBEDDED_BLANK ( JCARD,2,0,0,0,0,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-5
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,3,4,5,6,7,8,9 )! Issue warning if fields 6-9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

         IF (IERRFL(2) == 'Y') THEN
            JERR = JERR + 1
         ENDIF

      ENDIF

! Write data to file LINK1U if there were no errors

      IF (JERR == 0) THEN
         WRITE(L1U) SETID, CID, GID, SCALEF_AV, SCALEF_AA, (VEC(I),I=1,3)
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 ! *********************************************************************************************************************************

      END SUBROUTINE BD_RFORCE


      SUBROUTINE BD_SLOAD ( CARD, CC_LOAD_FND )

! Processes SLOAD Bulk Data Cards. Also, the set ID is written to array SLOAD_SIDS which is checked
! in subroutine LOADB to make sure that all set ID's requested in Case Control were found in the Bulk Data.
! A record is written to file LINK1W for each SLOAD Bulk Data pairs of set ID/force mag with the following data:

!   SETID, scalar point, load mag

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1W
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, LFORCE, LSUB, NFORCE, NSLOAD, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL, SUPWARN
      USE MODEL_STUF, ONLY            :  SLOAD_SIDS, SUBLOD

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_SLOAD'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD                ! A Bulk Data card
      CHARACTER( 1*BYTE),INTENT(INOUT):: CC_LOAD_FND(LSUB,2) ! 'Y' if B.D SLOAD card w/ same set ID (SID) as C.C. LOAD = SID
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)           ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: GRID_NO(3)          ! Grid ID  on the FORCE/MOMENT card
      INTEGER(LONG)                   :: I                   ! DO loop index
      INTEGER(LONG)                   :: JERR      = 0       ! A local error count
      INTEGER(LONG)                   :: NUM_PAIRS = 0       ! Bumber of pairs of SPOINT/FMAG on a SLOAD entry (can be up to 3)
      INTEGER(LONG)                   :: SETID     = 0       ! Set ID on the FORCE/MOMENT card


      REAL(DOUBLE)                    :: FMAG(3)             ! Force magnitude



! **********************************************************************************************************************************
! SLOAD Bulk Data Card routine

!   FIELD   ITEM           VARIABLE
!   -----   ------------   -------------
!    2      Load set ID    SLOAD_SIDS(I), SETID
!    3      Scalar point   GRID_NO
!    4      Load magnitude FMAG

! Fields (3,4) can be repeated in fields (5,6) and (7,8)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Initialize

      DO I=1,3
         GRID_NO(I) = 0
         FMAG(I)    = ZERO
      ENDDO

! Read and check data

      CALL I4FLD ( JCARD(2), JF(2), SETID )
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NSUB
            IF (SETID == SUBLOD(I,1)) THEN
               CC_LOAD_FND(I,1) = 'Y'
            ENDIF
         ENDDO
      ENDIF

      NUM_PAIRS = 0
      DO I=1,3
         IF (JCARD(3+2*(I-1))(1:) /= ' ') THEN
            NSLOAD = NSLOAD + 1
            NUM_PAIRS = NUM_PAIRS + 1
            CALL I4FLD ( JCARD(3+2*(I-1)), JF(3+2*(I-1)), GRID_NO(I) )
            CALL R8FLD ( JCARD(4+2*(I-1)), JF(4+2*(I-1)), FMAG(I) )
            SLOAD_SIDS(NSLOAD) = SETID
         ENDIF
      ENDDO

      IF      (NUM_PAIRS == 1) THEN                        ! Check for imbedded blanks in fields where there is data
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,0,0,0,0,0 )
      ELSE IF (NUM_PAIRS == 2) THEN
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,0,0,0 )
      ELSE IF (NUM_PAIRS == 3) THEN
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,0 )
      ENDIF
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )   ! Issue warning if field 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      DO I=2,8                                             ! Set JERR if any errors reading above data
         IF (IERRFL(I) == 'Y') THEN
            JERR = JERR + 1
         ENDIF
      ENDDO

! Write data to file LINK1W if there were no errors

      IF (JERR == 0) THEN
         DO I=1,NUM_PAIRS
            WRITE(L1W) SETID, GRID_NO(I), FMAG(I)
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1137 FORMAT(' *WARNING    : ',A,' ENTRY WITH SET ID = ',A,' HAS ZERO COMPONENTS')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_SLOAD


      SUBROUTINE BD_SLOAD0 ( CARD, NUM_PAIRS )

! Processes SLOAD Bulk Data Cards to count the number of SLOAD points/force mags on one entry

!   SETID, scalar point, load mag

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, JCARD_LEN, JF
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_SLOAD0'
      CHARACTER(LEN=*),INTENT(IN)     :: CARD                ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)           ! The 10 fields of characters making up CARD

      INTEGER(LONG), INTENT(OUT)      :: NUM_PAIRS           ! Number of pairs of SPOINT/force MAG on a SLOAD entry (can be up to 3)
      INTEGER(LONG)                   :: I                   ! DO loop index




! **********************************************************************************************************************************
! SLOAD Bulk Data Card routine

!   FIELD   ITEM           VARIABLE
!   -----   ------------   -------------
!    2      Load set ID    SLOAD_SIDS(I), SETID
!    3      Scalar point   GRID_NO
!    4      Load magnitude FMAG

! Fields (3,4) can be repeated in fields (5,6) and (7,8)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      NUM_PAIRS = 0
      DO I=1,3
         IF (JCARD(3+2*(I-1))(1:) /= ' ') THEN
            NUM_PAIRS = NUM_PAIRS + 1
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_SLOAD0

   END MODULE BULK_DATA_LOADS
