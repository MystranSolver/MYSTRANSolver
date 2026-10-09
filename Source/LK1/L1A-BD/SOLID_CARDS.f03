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

   MODULE SOLID_CARDS

   IMPLICIT NONE

   PRIVATE

   EXTERNAL :: ELEPRO

   PUBLIC :: BD_CHEXA, BD_CHEXA0, BD_CPENTA, BD_CPENTA0, BD_CTETRA, BD_CTETRA0, BD_PSOLID

   CONTAINS

      SUBROUTINE BD_CHEXA ( CARD, LARGE_FLD_INP, NUM_GRD )

! Processes CHEXA Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, JCARD_LEN, NCHEXA8, NCHEXA20, NEDAT, NELE
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ETYPE

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CHEXA'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: ID                ! Character value of element ID (field 2 of parent card)
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN=JCARD_LEN)        :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO
      CHARACTER(LEN=JCARD_LEN)        :: NAME              ! Field 1 of CARD

      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD           ! Number of GRID's + SPOINT's for the elem
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! CHEXA element Bulk Data Card routine

!   FIELD   ITEM                   ARRAY ELEMENT
!   -----   ------------   ---------------------------------

! Mandatory parent card

!    1      Element type   ETYPE(nele) = 'HEXA8' or 'HEXA20'
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+2)
!    4-9    Grids 1-6      EDAT(nedat+3) thru EDAT(nedat+8)

! Mandatory 2nd card (for 8 or 20 node HEXA):

!    2,3    Grid 7, 8      EDAT(nedat+9, 10)
!    4-9    Grids 9-14     EDAT(nedat+11) thru EDAT(nedat+16) if this is a HEXA 20 node element

! Possible 3rd card if this is a 20 node HEXA:

!    2-7    Grids 15-20    EDAT(nedat+17) thru EDAT(nedat+22) if this is a HEXA 20 node element


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      NAME = JCARD(1)
      ID   = JCARD(2)

! Set JCARD_EDAT to JCARD

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Read and check data

      CALL ELEPRO ( 'Y', JCARD_EDAT, 8, 8, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y' )

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,8,9 )   ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Required 2nd card:

      ETYPE(NELE) = '        '
      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN

         IF (JCARD(4)(1:) == ' ') THEN
            NCHEXA8  = NCHEXA8 + 1
            ETYPE(NELE)   = 'HEXA8   '
            NUM_GRD = 8
         ELSE
            NCHEXA20 = NCHEXA20 + 1
            ETYPE(NELE)   = 'HEXA20  '
            NUM_GRD = 20
         ENDIF

         DO I=1,10
            JCARD_EDAT(I) = JCARD(I)
         ENDDO

         IF      (ETYPE(NELE) == 'HEXA8   ') THEN
            CALL ELEPRO ( 'N', JCARD_EDAT, 2, 2, 'Y', 'Y', 'N', 'N', 'N', 'N', 'N', 'N' )
            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,0,0,0,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-3
            CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,8,9 )! Issue warning if fields 4-6 not blank
         ELSE IF (ETYPE(NELE) == 'HEXA20  ') THEN
            CALL ELEPRO ( 'N', JCARD_EDAT, 8, 8, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y' )
            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
         ENDIF
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields


      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1136) NAME, ID
         WRITE(F06,1136) NAME, ID

      ENDIF

! Required 3rd card (if TYPE = HEXA20):

      IF (ETYPE(NELE) == 'HEXA20') THEN

         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN

            DO I=1,10
               JCARD_EDAT(I) = JCARD(I)
            ENDDO

            CALL ELEPRO ( 'N', JCARD_EDAT, 6, 6, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'N', 'N' )

            CALL BD_IMBEDDED_BLANK( JCARD,2,3,4,5,6,7,0,0 )! Make sure that there are no imbedded blanks in fields 2-7
            CALL CARD_FLDS_NOT_BLANK(JCARD,0,0,0,0,0,0,8,9)! Issue warning if fields 8, 9 not blank
            CALL CRDERR ( CARD )                           ! CRDERR prints errors found when reading fields

         ELSE

            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1136) NAME, ID
            WRITE(F06,1136) NAME, ID

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

! **********************************************************************************************************************************

      END SUBROUTINE BD_CHEXA


      SUBROUTINE BD_CHEXA0 ( CARD, LARGE_FLD_INP, DELTA_LEDAT )

! Processes CHEXA Bulk Data Cards to determine how many words to allocate to array EDAT for this element

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, MEDAT_CHEXA8, MEDAT_CHEXA20
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CHEXA'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein

      INTEGER(LONG), INTENT(OUT)      :: DELTA_LEDAT       ! Delta number of words to add to LEDAT for this element




! **********************************************************************************************************************************
! This element must have at least 1 continuation card since there must be at leas 8 grids defined (all corner nodes). The parent
! card has the element ID, property ID and the first 6 grids. The first continuation card must contain nodes 7 and 8 (in fields
! 2 and 3) but can also define any or all of nodes 9 thru 14 (in fields 4 - 9). If there is a 2nd continuation card it can define
! any or all of nodes 15 thru 20 (in fields 2 - 7). If we find any more than the mandatory 8 nodes we will assume a full complement
! of 20 nodes for safety

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC0  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      IF (ICONT == 1) THEN                                 ! There was a 1st continuation card (mandatory)
         DELTA_LEDAT = MEDAT_CHEXA8
         DO I=4,9
            IF (JCARD(I)(1:) /= ' ') THEN                  ! Some of fields 4 thru 9 have data so we will assume a 20 node HEXA
               DELTA_LEDAT = MEDAT_CHEXA20
               EXIT
            ENDIF
         ENDDO
      ELSE                                                 ! There was no 1st contin card. This is error that will be caught later
         DELTA_LEDAT = MEDAT_CHEXA20                       ! For this error, set DELTA_LEDAT to largest until error is caught
      ENDIF


      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CHEXA0


      SUBROUTINE BD_CPENTA ( CARD, LARGE_FLD_INP, NUM_GRD )

! Processes CPENTA Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, JCARD_LEN, NCPENTA6, NCPENTA15, NEDAT, NELE, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ETYPE

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CPENTA'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: ID                ! Character value of element ID (field 2 of parent card)
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO
      CHARACTER(LEN(JCARD))           :: NAME              ! JCARD(1) from parent entry

      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD           ! Number of GRID's + SPOINT's for the elem
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! CPENTA element Bulk Data Card routine

!   FIELD   ITEM                   ARRAY ELEMENT
!   -----   ------------   ---------------------------------

! Mandatory parent card

!    1      Element type   ETYPE(nele) = 'PENTA6' or 'PENTA15'
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+1)
!    4-9    Grids 1-6      EDAT(nedat+3) thru EDAT(nedat+8)

! Possible 2nd card if this is a 15 node PENTA:

!    2-9    Grids 7-14     EDAT(nedat+9) thru EDAT(nedat+16) if this is a PENTA 15 node element

! Possible 3rd card if this is a 15 node PENTA:

!    2      Grid  15       EDAT(nedat+17) if this is a PENTA 15 node element


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      NAME = JCARD(1)
      ID   = JCARD(2)

! Set JCARD_EDAT to JCARD

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Read and check data

      CALL ELEPRO ( 'Y', JCARD_EDAT, 8, 8, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y' )

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,8,9 )   ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Check for continuation (if so must be a PENTA15 which then requires 3rd card also)

      ETYPE(NELE)(1:) = ' '
      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 0) THEN                                 ! Assume PENTA6 since there was no 2nd card

         NCPENTA6    = NCPENTA6 + 1                        ! Default values unless there is a non-blank continuation
         ETYPE(NELE) = 'PENTA6   '
         NUM_GRD     = 6

      ELSE

         IF (CARD(1:) /= ' ') THEN                         ! Only process continuation entries if 1st one is not totally blank
                                                           ! Assume PENTA15 since there was a 2nd card
            NCPENTA15   = NCPENTA15 + 1
            ETYPE(NELE) = 'PENTA15  '
            NUM_GRD     = 15

            DO I=1,10
               JCARD_EDAT(I) = JCARD(I)
            ENDDO

            CALL ELEPRO ( 'N', JCARD_EDAT, 8, 8, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y' )
            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
            CALL CRDERR ( CARD )

            IF (LARGE_FLD_INP == 'N') THEN
               CALL NEXTC  ( CARD, ICONT, IERR )           ! Required 3rd card for PENTA15
            ELSE
               CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
               CARD = CHILD
            ENDIF
            CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
            IF (ICONT == 1) THEN

               DO I=1,10
                  JCARD_EDAT(I) = JCARD(I)
               ENDDO

               CALL ELEPRO ( 'N', JCARD_EDAT, 1, 1, 'Y', 'N', 'N', 'N', 'N', 'N', 'N', 'N' )

               CALL BD_IMBEDDED_BLANK( JCARD,2,0,0,0,0,0,0,0 )
               CALL CARD_FLDS_NOT_BLANK(JCARD,0,3,4,5,6,7,8,9)
               CALL CRDERR ( CARD )

            ELSE

               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1136) NAME, ID
               WRITE(F06,1136) NAME, ID

            ENDIF

         ELSE                                              ! Since no cont data was read, default to PENTA6

            NCPENTA6    = NCPENTA6 + 1
            ETYPE(NELE) = 'PENTA6   '
            NUM_GRD     = 6
            ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

! **********************************************************************************************************************************

      END SUBROUTINE BD_CPENTA


      SUBROUTINE BD_CPENTA0 ( CARD, LARGE_FLD_INP, DELTA_LEDAT )

! Processes CPENTA Bulk Data Cards to determine how many words to allocate to array EDAT for this element

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, MEDAT_CPENTA6, MEDAT_CPENTA15
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CPENTA'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein

      INTEGER(LONG), INTENT(OUT)      :: DELTA_LEDAT       ! Delta number of words to add to LEDAT for this element




! **********************************************************************************************************************************
! This element does not need any continuation cards. The parent must define 6 nodes of the PENTA. Continuation cards can define
! optional nodes 7 thru 15 (in fields 2-9 of continuation card 1 abd node 15 on continuation card 2). If there are no continuation
! cards we assume a 6 node PENTA. If there is a continuation card we will assume it is a 15 node PENTA for safety

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC0  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      IF (ICONT == 1) THEN                                 ! There was a 1st continuation card
         DELTA_LEDAT = MEDAT_CPENTA15
      ELSE                                                 ! There was no 1st contin card so assume a 6 node PENTA
         DELTA_LEDAT = MEDAT_CPENTA6
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CPENTA0


      SUBROUTINE BD_CTETRA ( CARD, LARGE_FLD_INP, NUM_GRD )

! Processes CTETRA Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, FATAL_ERR, NCTETRA4, NCTETRA10, NEDAT, NELE
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ETYPE

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CTETRA'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: ID                ! Character value of element ID (field 2 of parent card)
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO
      CHARACTER(LEN(JCARD))           :: NAME              ! JCARD(1) from parent entry

      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD           ! Number of GRID's + SPOINT's for the elem
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
! CTETRA element Bulk Data Card routine

!   FIELD   ITEM                   ARRAY ELEMENT
!   -----   ------------   ---------------------------------

! Mandatory parent card

!    1      Element type   ETYPE(nele) = 'TETRA4' or 'TETRA10'
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+1)
!    4-9    Grids 1-6      EDAT(nedat+3) thru EDAT(nedat+8)

! Possible 2nd card if this is a 10 node TETRA:

!    2-5    Grids 7-10     EDAT(nedat+9) thru EDAT(nedat+16) if this is a TETRA 10 node element


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      NAME = JCARD(1)
      ID   = JCARD(2)

! Set JCARD_EDAT to JCARD

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Read and check data

      IF (JCARD(8)(1:) == ' ') THEN                        ! If field 8 is blank then this must be a TETRA 4
         CALL ELEPRO ( 'Y', JCARD_EDAT, 6, 6, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'N', 'N' )
         ETYPE(NELE) = 'TETRA4  '
         NCTETRA4 = NCTETRA4 + 1
         NUM_GRD = 4
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,0,0 )
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )
      ELSE                                                 ! This must be a TETRA 10 if field 8 is not blank
         CALL ELEPRO ( 'Y', JCARD_EDAT, 8, 8, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y' )
         ETYPE(NELE) = 'TETRA10 '
         NCTETRA10 = NCTETRA10 + 1
         NUM_GRD = 10
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )
      ENDIF
      CALL CRDERR ( CARD )

! 2nd card - required if TETRA 10:

      IF (ETYPE(NELE) == 'TETRA10 ') THEN

         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN

            DO I=1,10
               JCARD_EDAT(I) = JCARD(I)
            ENDDO

            IF (ETYPE(NELE) == 'TETRA10  ') THEN
               CALL ELEPRO ( 'N', JCARD_EDAT, 4, 4, 'Y', 'Y', 'Y', 'Y', 'N', 'N', 'N', 'N' )
               CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
            ENDIF


         ELSE

            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1136) NAME, ID
            WRITE(F06,1136) NAME, ID

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

! **********************************************************************************************************************************

      END SUBROUTINE BD_CTETRA


      SUBROUTINE BD_CTETRA0 ( CARD, LARGE_FLD_INP, DELTA_LEDAT )

! Processes CTETRA Bulk Data Cards to determine how many words to allocate to array EDAT for this element

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, MEDAT_CTETRA4, MEDAT_CTETRA10
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CTETRA'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG), INTENT(OUT)      :: DELTA_LEDAT       ! Delta number of words to add to LEDAT for this element
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! This element does not need any continuation cards. The parent must define the 4 mandatory corner nodes of the TETRA and can also
! define nodes 5 and 6 (in fields 8-9) of a higher order TETRA.  A continuation card can define optional nodes 7 thru 10
! of a higher order TETRA (in fields 2-5. If there is a continuation card we assume a 10 node TETRA for safety. If not, we check
! fields 8 and 9 of parent card and if they are blank, assume a 4 node TETRA; otherwise assume a 10 node TETRA (again for safety)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC0  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      IF (ICONT == 1) THEN
         DELTA_LEDAT = MEDAT_CTETRA10
      ELSE
         IF ((JCARD(8)(1:) == ' ') .AND. (JCARD(9)(1:) == ' ')) THEN
            DELTA_LEDAT = MEDAT_CTETRA4
         ELSE
            DELTA_LEDAT = MEDAT_CTETRA10
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CTETRA0


      SUBROUTINE BD_PSOLID ( CARD, IOR3D )

! Processes PSOLID Bulk Data Cards. Reads and checks:

!  1) Property ID and material ID and enter them into array PSOLID

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPSOLID, NPSOLID, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  PSOLID

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, LEFT_ADJ_BDFLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PSOLID'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: CHR_FLD           ! Character data from 1 card field that has been left adjusted

      INTEGER(LONG), INTENT(OUT)      :: IOR3D             ! Integration order for this PSOLID entry
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: ID        = 0     ! An integer ID read from a field of this card




! **********************************************************************************************************************************
! PSOLID Bulk Data Card routine

!   FIELD   ITEM            ARRAY ELEMENT
!   -----   ------------    -------------
!    2      Prop ID          PSOLID(npsolid,1)
!    3      Mat ID           PSOLID(npsolid,2)
!    4      Matl coord sys   PSOLID(npsolid,3)
!    5      Int order        PSOLID(npsolid,4)
!    6      Stress locations PSOLID(npsolid,5) ***** Not currently used *****
!    7      Int scheme       PSOLID(npsolid,6)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NPSOLID = NPSOLID+1
!xx   IF (NPSOLID > LPSOLID) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LPSOLID
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LPSOLID
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data on parent card

      PSOLID(NPSOLID,1) = 0
      CALL I4FLD ( JCARD(2), JF(2), ID )                   ! Read property ID and enter into array PSOLID
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NPSOLID-1
            IF (ID == PSOLID(J,1)) THEN                    ! Make sure that this is a unique PSOLID entry
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),ID
               WRITE(F06,1145) JCARD(1),ID
               EXIT
            ENDIF
         ENDDO
         PSOLID(NPSOLID,1) = ID
      ENDIF

      PSOLID(NPSOLID,2) = 0
      CALL I4FLD ( JCARD(3), JF(3), ID )                   ! Read material ID and enter into array PSOLID
      IF (IERRFL(3) == 'N') THEN
         IF (ID <= 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', ID
            WRITE(F06,1192) JF(3), JCARD(1), JCARD(2), ' > 0 ', ID
         ELSE
            PSOLID(NPSOLID,2) = ID
         ENDIF
      ENDIF

      PSOLID(NPSOLID,3) = -1                               ! Mat'l coord sys ID default value = -1 to distinguish from CID > or = 0
      IF (JCARD(4)(1:) /= ' ') THEN                        ! Read material coord system ID and enter into array PSOLID
         CALL I4FLD ( JCARD(4), JF(4), ID )
         IF (IERRFL(4) == 'N') THEN
            IF (ID < -1) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1192) JF(4), JCARD(1), JCARD(2), ' >= -1 ', ID
               WRITE(F06,1192) JF(4), JCARD(1), JCARD(2), ' >= -1 ', ID
            ELSE
               PSOLID(NPSOLID,3) = ID
            ENDIF
         ENDIF
      ENDIF

      PSOLID(NPSOLID,4) = 0                                ! Read integration order and enter into array PSOLID
      CHR_FLD = JCARD(5)
      CALL LEFT_ADJ_BDFLD ( CHR_FLD )
      IF       (CHR_FLD(1:) == ' ') THEN
         PSOLID(NPSOLID,4) = -1                            ! Temporary value that will be changed in ELMDAT1
      ELSE IF ((CHR_FLD(1:8) == 'TWO     ') .OR. (CHR_FLD(1:8) == '2       ')) THEN
         PSOLID(NPSOLID,4) = 2
      ELSE IF ((CHR_FLD(1:8) == 'THREE   ') .OR. (CHR_FLD(1:8) == '3       ')) THEN
         PSOLID(NPSOLID,4) = 3
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1104) JCARD(5),JF(5),JCARD(1),JCARD(2)
         WRITE(F06,1104) JCARD(5),JF(5),JCARD(1),JCARD(2)
      ENDIF
      IOR3D = PSOLID(NPSOLID,4)                            ! Return this    integration order to calling subr

      PSOLID(NPSOLID,5) = 0                                ! This feature - stress location definition, not currently used)

      PSOLID(NPSOLID,6) = 0                                ! Read integration scheme and enter into array PSOLID
      CHR_FLD = JCARD(7)
      CALL LEFT_ADJ_BDFLD ( CHR_FLD )
      IF      (CHR_FLD(1:8) == '1       ') THEN
         ID = -3                                           ! Temporary value that will be changed in ELMDAT1
      ELSE IF (CHR_FLD(1:8) == '0       ') THEN
         ID = -2                                           ! Temporary value that will be changed in ELMDAT1
      ELSE IF (CHR_FLD(1:8) == '        ') THEN
         ID = -1                                           ! Temporary value that will be changed in ELMDAT1
      ELSE IF (CHR_FLD(1:8) == 'REDUCED ') THEN
         ID = 0
      ELSE IF (CHR_FLD(1:8) == 'FULL    ') THEN
         ID = 1
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1104) JCARD(7),JF(7),JCARD(1),JCARD(2)
         WRITE(F06,1104) JCARD(7),JF(7),JCARD(1),JCARD(2)
      ENDIF
      PSOLID(NPSOLID,6) = ID

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,0,7,0,0 )     ! Make sure that there are no imbedded blanks in fields 2-5,7
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,0,8,9 )   ! Issue warning if fields 6, 8 and 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1104 FORMAT(' *ERROR  1104: INVALID ENTRY = "',A,'" IN FIELD ',I3,' OF ',A,' ENTRY WITH ID = ',A)

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1192 FORMAT(' *ERROR  1192: ID IN FIELD ',I3,' OF ',A,A,' MUST BE ',A,' BUT IS = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PSOLID

   END MODULE SOLID_CARDS
