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

   MODULE SPRING_BUSH_MASS

   IMPLICIT NONE

   PRIVATE

   EXTERNAL :: ELEPRO

   PUBLIC :: BD_CELAS1, BD_CELAS2, BD_CELAS3, BD_CELAS4, BD_PELAS, BD_CBUSH, BD_CBUSH0, BD_PBUSH, BD_CMASS1, BD_CMASS2, BD_CMASS3, BD_CMASS4, BD_PMASS, BD_CONM2

   CONTAINS

      SUBROUTINE BD_CELAS1 ( CARD )

! Processes CELAS1 Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, MEDAT_CELAS1, NCELAS1, NELE, NEDAT
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CELAS1'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: CELAS_ELID        ! Field 2 of CELAS1 card (this CELAS1's elem ID)
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD but with fields 5 and 6 switched to get G.P.'s together in EDAT

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IDOF              ! Displ component (1,2,3,4,5 or 6) that one end of CELSA conn. to




! **********************************************************************************************************************************
! CELAS1 scalar spring element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele) =E1 for CELAS1
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+2)
!    4      Grid-A         EDAT(nedat+3)
!    5      Comp-A         EDAT(nedat+5)
!    6      Grid-B         EDAT(nedat+4)
!    7      Comp-B         EDAT(nedat+6)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      CELAS_ELID = JCARD(2)

! Make JCARD_EDAT, which is the version that will have JCARD fields 5, 6 switched when subr ELEPRO called

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Check property ID field. Set to element ID if blank

      IF (JCARD(3)(1:) == ' ') THEN
         JCARD_EDAT(3) = JCARD(2)
      ENDIF

! Flip Comp-A and Grid-B in JCARD_EDAT so when ELEPRO runs it will have Grid-A and Grid-B back-to-back

      JCARD_EDAT(5) = JCARD(6)
      JCARD_EDAT(6) = JCARD(5)

      CALL ELEPRO ( 'Y', JCARD_EDAT, 6, MEDAT_CELAS1, 'Y', 'Y', 'Y', 'Y', 'N', 'N', 'N', 'N' )
      NCELAS1 = NCELAS1+1
      ETYPE(NELE) = 'ELAS1   '

! Check to make sure that numbers in fields 5 and 7 are valid component numbers

      IF (IERRFL(JF(5)) == 'N') THEN
         CALL I4FLD ( JCARD(5), JF(5), IDOF )
         IF ((IDOF <= 0) .OR. (IDOF > 6)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1133) IDOF, JF(5), CELAS_ELID
            WRITE(F06,1133) IDOF, JF(5), CELAS_ELID
         ENDIF
      ENDIF

      IF (IERRFL(JF(7)) == 'N') THEN
         CALL I4FLD ( JCARD(7), JF(7), IDOF )
         IF ((IDOF <= 0) .OR. (IDOF > 6)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1133) IDOF, JF(7), CELAS_ELID
            WRITE(F06,1133) IDOF, JF(7), CELAS_ELID
         ENDIF
      ENDIF

! Issue warning if fields 8, 9 are not blank

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,0,0 )   ! Make sure that there are no imbedded blanks in fields 2-7
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )   ! Issue warning if fields 8, 9 are not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1133 FORMAT(' *ERROR  1133: INVALID COMPONEMT NUMBER = ',I8,' IN FIELD ',I2,' ON CELAS1 ID = ',A8,' .MUST BE SINGLE DIGIT 1-6')

! **********************************************************************************************************************************

      END SUBROUTINE BD_CELAS1


      SUBROUTINE BD_CELAS2 ( CARD )

! Processes CELAS2 Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, MEDAT_CELAS2, NCELAS2, NELE, NEDAT, NPELAS
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE, PELAS, RPELAS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CELAS2'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: CELAS_ELID        ! Field 2 of CELAS2 card (this CELAS2's elem ID)
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD but with fields 5 and 6 switched to get G.P.'s together in EDAT

      INTEGER(LONG)                   :: ELEM_ID           ! Elem ID from field 2
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: I4INP             ! An integer read
      INTEGER(LONG)                   :: IDOF              ! Displ component (1,2,3,4,5 or 6) that one end of CELSA conn. to
      INTEGER(LONG)                   :: IERR              ! Error count


      REAL(DOUBLE)                    :: R8INP             ! A real value read



! **********************************************************************************************************************************
! CELAS2 scalar spring element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele) =E1 for CELAS2
!    2      Element ID     EDAT(nedat+1)
!    3      Stiffness, K   RPELAS(npelas,1)
!    4      Grid-A         EDAT(nedat+3)
!    5      Comp-A         EDAT(nedat+5)
!    6      Grid-B         EDAT(nedat+4)
!    7      Comp-B         EDAT(nedat+6)
!    8      Damping, GE    RPELAS(npelas,2)
!    9      Str rec, S     RPELAS(npelas,3)
!  none     Prop ID, PID   EDAT(nedat,2) (created: PID = -EID)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! First, check that fields 2-9 have the proper data type (we are going to have to rearrange the fields prior to calling ELEPRO).
! If any erors, return

      IERR = 0

      CALL I4FLD ( JCARD(2), JF(2), I4INP )
      CALL R8FLD ( JCARD(3), JF(3), R8INP )
      CALL I4FLD ( JCARD(4), JF(4), I4INP )
      CALL I4FLD ( JCARD(5), JF(5), I4INP )
      CALL I4FLD ( JCARD(6), JF(6), I4INP )
      CALL I4FLD ( JCARD(7), JF(7), I4INP )
      CALL R8FLD ( JCARD(8), JF(8), R8INP )
      CALL R8FLD ( JCARD(9), JF(9), R8INP )

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )  ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )
      IF (IERR > 0) THEN
         RETURN
      ENDIF

! Get elem ID from field 2 so we can use negative of it as property ID

      CELAS_ELID = JCARD(2)
      NPELAS = NPELAS + 1
      CALL I4FLD ( JCARD(2), JF(2), ELEM_ID )
      PELAS(NPELAS,1) = -ELEM_ID

! Write real data to array RPELAS (so we can rewrite JCARD(3) to be a prop ID)

      CALL R8FLD ( JCARD(3), JF(3), RPELAS(NPELAS,1) )
      CALL R8FLD ( JCARD(8), JF(8), RPELAS(NPELAS,2) )
      CALL R8FLD ( JCARD(9), JF(9), RPELAS(NPELAS,3) )

! Make JCARD_EDAT, which is the version that will have JCARD sent to subr ELEPRO

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Now change JCARD(3) to be a property ID so that subr ELEPRO will handle EDAT data correctly. We want PID = -EID but we send
! JCARD(3) = JCARD(2) (which has PID = EID) to ELEPRO. When ELEPRO returns change term in EDAT for PID to be -PID (i.e. PID = -EID)

      JCARD_EDAT(3) = JCARD_EDAT(2)

! Flip Comp-A and Grid-B in JCARD_EDAT so when ELEPRO runs it will have Grid-A and Grid-B back-to-back

      JCARD_EDAT(5) = JCARD(6)
      JCARD_EDAT(6) = JCARD(5)
                                                           ! Do not check fields. That was already done above
      CALL ELEPRO ( 'Y', JCARD_EDAT, 6, MEDAT_CELAS2, 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N' )
      NCELAS2 = NCELAS2+1
      ETYPE(NELE) = 'ELAS2   '

! Now change PID in EDAT to -PID

      EDAT(NEDAT-4) = -EDAT(NEDAT-4)

! Check to make sure that numbers in fields 5 and 7 are valid component numbers

      IF (IERRFL(JF(5)) == 'N') THEN
         CALL I4FLD ( JCARD(5), JF(5), IDOF )
         IF ((IDOF <= 0) .OR. (IDOF > 6)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1133) IDOF,JF(5),CELAS_ELID
            WRITE(F06,1133) IDOF,JF(5),CELAS_ELID
         ENDIF
      ENDIF

      IF (IERRFL(JF(7)) == 'N') THEN
         CALL I4FLD ( JCARD(7), JF(7), IDOF )
         IF ((IDOF <= 0) .OR. (IDOF > 6)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1133) IDOF,JF(7),CELAS_ELID
            WRITE(F06,1133) IDOF,JF(7),CELAS_ELID
         ENDIF
      ENDIF

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,8,9 )   ! Make sure that there are no imbedded blanks in fields 2-7
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1133 FORMAT(' *ERROR  1133: INVALID COMPONEMT NUMBER = ',I8,' IN FIELD ',I2,' ON CELAS2 ID = ',A8,' .MUST BE SINGLE DIGIT 1-6')

! **********************************************************************************************************************************

      END SUBROUTINE BD_CELAS2


      SUBROUTINE BD_CELAS3 ( CARD )

! Processes CELAS3 Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, MEDAT_CELAS3, NCELAS3, NELE, NEDAT
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CELAS3'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD but with fields 5 and 6 switched to get G.P.'s together in EDAT

      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
! CELAS3 scalar spring element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele) =E1 for CELAS3
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+2)
!    4      Scalar point A EDAT(nedat+3)
!    5      Scalar point B EDAT(nedat+5)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Make JCARD_EDAT, which is the version that will have JCARD fields 5, 6 switched when subr ELEPRO called

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Check property ID field. Set to element ID if blank

      IF (JCARD(3)(1:) == ' ') THEN
         JCARD_EDAT(3) = JCARD(2)
      ENDIF

      CALL ELEPRO ( 'Y', JCARD_EDAT, 4, MEDAT_CELAS3, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'N', 'N' )
      NCELAS3 = NCELAS3+1
      ETYPE(NELE) = 'ELAS3   '

! Issue warning if fields 6-9 are not blank

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,0,0,0,0 )   ! Make sure that there are no imbedded blanks in fields 2-7
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )   ! Issue warning if fields 6-9 are not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CELAS3


      SUBROUTINE BD_CELAS4 ( CARD )

! Processes CELAS4 Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, MEDAT_CELAS4, NCELAS4, NELE, NEDAT, NPELAS
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE, PELAS, RPELAS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CELAS4'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD but with fields 5 and 6 switched to get G.P.'s together in EDAT

      INTEGER(LONG)                   :: ELEM_ID           ! Elem ID from field 2
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: I4INP             ! Integer value read from a field of the CELAS4 entry
      INTEGER(LONG)                   :: IERR              ! Error count


      REAL(DOUBLE)                    :: R8INP             ! Real value read from a field on the PSHEAR entry



! **********************************************************************************************************************************
! CELAS4 scalar spring element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele) =E1 for CELAS4
!    2      Element ID     EDAT(nedat+1)
!           Property ID    EDAT(nedat+2) (fictitous property PID = -EID)
!    3      Stiffness      PELAS(npelas,1)
!    4      Scalar point A EDAT(nedat+3)
!    5      Scalar point B EDAT(nedat+5)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! First, check that fields 2-5 have the proper data type (we are going to have to rearrange the fields prior to calling ELEPRO).
! If any erors, return

      IERR = 0

      CALL I4FLD ( JCARD(2), JF(2), I4INP )
      CALL R8FLD ( JCARD(3), JF(3), R8INP )
      CALL I4FLD ( JCARD(4), JF(4), I4INP )
      CALL I4FLD ( JCARD(5), JF(5), I4INP )

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,0,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-5
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )! Issue warning if fields 6, 7, 8, 9 not blank
      CALL CRDERR ( CARD )
      IF (IERR > 0) THEN
         RETURN
      ENDIF

! Get elem ID from field 2 so we can use negative of it as property ID

      NPELAS = NPELAS + 1
      CALL I4FLD ( JCARD(2), JF(2), ELEM_ID )
      IF (IERRFL(2) == 'N') THEN
         PELAS(NPELAS,1) = -ELEM_ID
      ENDIF

! Write real data to array RPELAS (so we can rewrite JCARD(3) to be a prop ID)

      CALL R8FLD ( JCARD(3), JF(3), RPELAS(NPELAS,1) )

! Make JCARD_EDAT, which is the version that will have JCARD fields 5, 6 switched when subr ELEPRO called

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Now change JCARD(3) to be a property ID so that subr ELEPRO will handle EDAT data correctly. We want PID = -EID but we send
! JCARD(3) = JCARD(2) (which has PID = EID) to ELEPRO. When ELEPRO returns we change term in EDAT for PID to be -PID (i.e. -EID)

      JCARD_EDAT(3) = JCARD_EDAT(2)
                                                           ! Do not check fields. That was already done above
      CALL ELEPRO ( 'Y', JCARD_EDAT, 4, MEDAT_CELAS4, 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N' )
      NCELAS4 = NCELAS4+1
      ETYPE(NELE) = 'ELAS4   '

! Now change PID in EDAT to -PID

      EDAT(NEDAT-2) = -EDAT(NEDAT-2)



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CELAS4


      SUBROUTINE BD_PELAS ( CARD )

! Processes PELAS Bulk Data Cards. Read and check

!  1) Property ID and enter into array PELAS
!  2) Stiffness, damping, stress recovery coeff. and enter into array RPELAS

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPELAS, NPELAS
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PELAS, RPELAS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PELAS'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: PROP_ID   = 0     ! Property ID (field 2 of this property card)




! **********************************************************************************************************************************
! PELAS Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Prop ID         PELAS(npelas,1)
!    3      Sp. Rate - K   RPELAS(npelas,1)
!    4      Damping - GE   RPELAS(npelas,2)
!    5      Stress recov.  RPELAS(npelas,3)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NPELAS = NPELAS+1

! Read and check data

      CALL I4FLD ( JCARD(2), JF(2), PROP_ID )              ! Read property ID and enter into array PELAS
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NPELAS-1
            IF (PROP_ID == PELAS(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1), PROP_ID
               WRITE(F06,1145) JCARD(1), PROP_ID
               EXIT
            ENDIF
         ENDDO
         PELAS(NPELAS,1) = PROP_ID
      ENDIF

      DO J = 1,3                                           ! Read real property values in fields 3-5
         CALL R8FLD ( JCARD(J+2), JF(J+2), RPELAS(NPELAS,J) )
      ENDDO

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,0,0,0,0 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )   ! Issue warning if fields 6, 7, 8, 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PELAS


      SUBROUTINE BD_CBUSH ( CARD, LARGE_FLD_INP )

! Processes CBUSH Bulk Data card:

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LBUSHOFF, LVVEC, MEDAT_CBUSH,&
                                         NBUSHOFF, NCBUSH, NEDAT, NELE, NVVEC, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, HALF
      USE PARAMS, ONLY                :  EPSIL, SUPWARN
      USE MODEL_STUF, ONLY            :  BUSHOFF, EDAT, ETYPE, VVEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, LEFT_ADJ_BDFLD, R8FLD
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CBUSH'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: EID               ! Field 2 of CBUSH card
      CHARACTER( 1*BYTE)              :: FOUND     = 'N'   ! 'Y' if the V vec is one that is already stored in array VVEC
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of 8 characters making up CARD
      CHARACTER(LEN=JCARD_LEN)        :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO
      CHARACTER(LEN=JCARD_LEN)        :: NAME              ! Name of this entry (field 1)
      CHARACTER( 9*BYTE)              :: VVEC_TYPE         ! Type of V vector on this CBUSH

      INTEGER(LONG)                   :: CID               ! Coord sys ID for v vector (required under some circumstances)
      INTEGER(LONG)                   :: G0   = 0          ! Grid specifying V vector for this CBUSH, if input
      INTEGER(LONG)                   :: I4INP     = 0     ! A value read from input file that should be an integer value
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: JERR      = 0     ! A local error count
      INTEGER(LONG)                   :: NEDAT_START       ! Value of NEDAT at start of this subr
      INTEGER(LONG)                   :: OCID              ! Coord sys ID for offsets
      INTEGER(LONG)                   :: VVEC_NUM  = 0     ! V vector number


      REAL(DOUBLE)                    :: VV(3)             ! The 3 components of the V vector for this CBUSH elem
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: R8INP     = ZERO  ! A value read from input file that should be a real value

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! CBUSH element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele)
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+2)
!    -      Num of grids   EDAT)nedat+3
!    4      Grid A         EDAT(nedat+4)
!    5      Grid B         EDAT(nedat+5)
!    6-8    V-Vector       EDAT(nedat+6) (see VVEC explanation below)
!    9      CID            EDAT(nedat+7) Elem coord sys identification. 0 is basic. blank means to use G0 or X1,2,3

! on optional second card:
!    2      S              Location of spring/damper (def = 0.5). This is a relative distance = offset/BUSH length
!    3      OCID           EDAT(nedat+8) Offset coord sys ID
!   none    NBUSHOFF       EDAT(nedat+9)
!    4-6    S1,2,3         Offset vector. The Si are actual offsets

! NOTES:

! If V-vector is specified via a grid point then EDAT(nedat+5) is set to that grid number.
! If V-vector is specified via an actual vector, the vector is loaded into array VVEC(NVVEC,J) (J=1,2,3) unless
! a vector equal to it has been put in VVEC. EDAT(nedat+5) is set equal to -NVVEC, where NVVEC is the row number in array VVEC.

! Offsets are in fields 4 - 9 of the first continuation card. If there are no offsets for this element, a zero is entered
! in array EDAT(nedat7).

      NEDAT_START = NEDAT
      EPS1 = EPSIL(1)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      NAME = JCARD(1)
      EID  = JCARD(2)

! Set JCARD_EDAT to JCARD

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Check property ID field. Set to EID if blank

      IF (JCARD(3)(1:) == ' ') THEN
         JCARD_EDAT(3) = JCARD(2)
      ENDIF

! Make sure that grids in fields 4 and 5 are different

      CALL LEFT_ADJ_BDFLD ( JCARD(4) )
      CALL LEFT_ADJ_BDFLD ( JCARD(5) )
      IF (JCARD(4) == JCARD(5)) THEN
         WRITE(F06,*) ' * ERROR : GRIDS ON CBUSH CANNOT BE SAME'
      ENDIF

! Call ELEPRO to increment NELE and load some of the connection data into array EDAT

      CALL ELEPRO ( 'Y', JCARD_EDAT, 4, MEDAT_CBUSH , 'Y', 'Y', 'Y', 'Y', 'N', 'N', 'N', 'N' )
      NCBUSH = NCBUSH+1
      ETYPE(NELE)(1:4) = 'BUSH'

! Get the V vector for this CBUSH. If field 9 is not blank, ignore fields 6-8 since v vector will be defined using CID in field 9
! The following sets slots 5 and 6 in EDAT since VVEC can be defined by:
!   (1) A vector defined in field 6,7,8 (field 6 can be a grid or 6-8 can be a vector), or
!   (2) CID in field 9. In this case the VVEC will be along the directions defined by CID and will be determined in a later subr
! If CID is blank then V must be defined by G0 in field 6 or Xi in fields 6-8

      VVEC_TYPE = 'UNDEFINED'
vec:  IF (JCARD(9)(1:) == ' ') THEN                        ! CID field is blank so VVEC should be defined in fields 6,7,8

         DO J=1,3
            VV(J) = ZERO
         ENDDO

         DO J=1,JCARD_LEN                                  ! See if there is an actual V vector.
            IF ((JCARD(6)(J:J) == '.') .OR. (JCARD(7)(J:J) == '.') .OR. (JCARD(8)(J:J) == '.')) THEN
               VVEC_TYPE = 'VECTOR   '
               EXIT
            ENDIF
         ENDDO
         IF (VVEC_TYPE == 'VECTOR   ') THEN                ! If there is an actual V vector, get components
            LVVEC = LVVEC + 1
            JERR = 0
            DO J=1,3
               CALL R8FLD ( JCARD(J+5), JF(J+5), R8INP )
               IF (IERRFL(J+5) == 'N') THEN
                  VV(J) = R8INP
               ELSE
                  JERR = JERR + 1
               ENDIF
            ENDDO
            IF (JERR /= 0) THEN
               VVEC_TYPE = 'ERROR    '                     ! Found error in V vector components, so reset VVEC_TYPE
            ENDIF
         ELSE                                              ! Check to see if there is a grid no. for specifying VVEC
            IF ((JCARD(6)(1:) /= ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN
               VVEC_TYPE = 'GRID     '
               CALL I4FLD ( JCARD(6), JF(6), I4INP )
               IF (IERRFL(6) == 'N') THEN
                  G0 = I4INP
                  IF (G0 < 0) THEN
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1187) NAME, EID, G0
                     WRITE(F06,1187) NAME, EID, G0
                  ENDIF
               ELSE
                  VVEC_TYPE = 'ERROR    '                  ! Found error in field 6, so reset VVEC_TYPE
               ENDIF
!xx         ELSE
!xx            IF ((JCARD(6)(1:) == ' ') .AND. (JCARD(7)(1:) == ' ') .AND. (JCARD(8)(1:) == ' ')) THEN
!xx               VVEC_TYPE = 'UNDEFINED'
!xx               FATAL_ERR = FATAL_ERR + 1
!xx               WRITE(ERR,1188) NAME, EID
!xx               WRITE(F06,1188) NAME, EID
!xx            ELSE
!xx               VVEC_TYPE = 'ERROR    '
!xx               FATAL_ERR = FATAL_ERR + 1
!xx               WRITE(ERR,1186) NAME, EID
!xx               WRITE(F06,1186) NAME, EID
!xx            ENDIF
            ENDIF
         ENDIF
                                                           ! Load V vector data into EDAT and into VVEC, if not already there
         IF (VVEC_TYPE == 'GRID     ') THEN

            EDAT(NEDAT_START+6) = G0                       ! --- Slot 6 in EDAT is VVEC defined by grid G0

         ELSE IF (VVEC_TYPE == 'VECTOR   ') THEN

            FOUND = 'N'                                    ! --- See if there is already a VVEC with these components
            DO J=1,NVVEC
               IF ((DABS(VV(1) - VVEC(J,1)) < EPS1) .AND. (DABS(VV(2) - VVEC(J,2)) < EPS1) .AND.                                   &
                   (DABS(VV(3) - VVEC(J,3)) < EPS1)) THEN
                  VVEC_NUM = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO

            IF (FOUND == 'N') THEN
               NVVEC = NVVEC + 1
               IF (NVVEC > LVVEC) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1132) SUBR_NAME,LVVEC
                  WRITE(F06,1132) SUBR_NAME,LVVEC
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
               VVEC_NUM = NVVEC
               DO J=1,3
                  VVEC(NVVEC,J) = VV(J)
               ENDDO
            ENDIF

            EDAT(NEDAT_START+6) = -VVEC_NUM                ! --- Slot 6 in EDAT is VVEC defined by neg of a vector number

         ELSE

            EDAT(NEDAT_START+6) = 0                        ! --- Slot 6 in EDAT is undefined VVEC

         ENDIF

         EDAT(NEDAT_START+7) = -99                         ! --- Slot 7 in EDAT is for CID. Use CID = -99 when CID field is blank

      ELSE                                                 ! --- CID field is not blank so VVEC is defined by CID

         EDAT(NEDAT_START+6) = 0                           ! --- Load 0 into EDAT slot 5 for VVEC since field 9 not blank

         EDAT(NEDAT_START+7) = 0                           ! --- Initialize. 0 (does not mean basic)
         CALL I4FLD ( JCARD(9), JF(9), CID )
         IF (IERRFL(9) == 'N') THEN
            IF (CID >= 0) THEN
               EDAT(NEDAT_START+7) = CID                   ! --- Slot 6 in EDAT is for CID. Use actual value if no read error
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1189) NAME, EID, JF(9), CID
               WRITE(F06,1189) NAME, EID, JF(9), CID
            ENDIF
         ENDIF
                                                           ! Issue warning since both 6-8 and 9 are non blank
         IF ((JCARD(6)(1:) /= ' ') .OR. (JCARD(7)(1:) /= ' ') .OR. (JCARD(8)(1:) /= ' ')) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,1002) NAME, EID
            IF (SUPWARN == 'N') THEN
               WRITE(F06,1002) NAME, EID
            ENDIF
         ENDIF

      ENDIF vec

! Write warnings and errors for parent entry, if any

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Optional Second Card: NOTE: The offset on this entry is relative to grid A and has only 3 components.
! The offset from grid B will be to the same location in space as where the offset from grid A is (since the CBUSH is a zero
! length element). Thus, we really do not need terms in array OFFDIS for grid B. However, put the grid A offset values in the
! grid B offset position anyway.

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      OCID = -99                                           ! Initially set OCID = -99 since we will test for -1 and > 0 below
      IF (ICONT == 1) THEN

         IF (CARD(1:) /= ' ') THEN                         ! Only process continuation entries if 1st one is not totally blank

            CALL LEFT_ADJ_BDFLD ( JCARD(3) )               ! Get OCID. OCID=-1 is default and means use elem sys and S, not S1,2,3
            IF ((JCARD(3)(1:) == ' ') .OR. (JCARD(3)(1:2) == '-1')) THEN
               OCID = -1                                   ! OCID = -1 is default and means use elem system and S, not S1,2,3
            ELSE                                           ! Field 3 not blank so get value to use for OCID
               CALL I4FLD ( JCARD(3), JF(3), I4INP )
               IF (IERRFL(2) == 'N') THEN
                  OCID = I4INP
               ELSE
                  OCID = -99                               ! Use OCID = -99 to indicate an error
               ENDIF
               IF (OCID < -1) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1180) JCARD(3), NAME, EID
                  WRITE(F06,1180) JCARD(3), NAME, EID
               ENDIF
            ENDIF
            EDAT(NEDAT_START+8) = OCID                     ! Slot 8 in EDAT is the OCID value. Slot 7 will be NBUSHOFF

            IF      (OCID == -1) THEN                      ! Get components of BUSHOFF
               CALL R8FLD ( JCARD(2), JF(2), R8INP )
               IF (IERRFL(1) == 'N') THEN
                  NBUSHOFF = NBUSHOFF + 1
                  EDAT(NEDAT_START+9) = NBUSHOFF           ! Slot 9 in EDAT is for the offset key
                  BUSHOFF(NBUSHOFF,1) = R8INP
                  BUSHOFF(NBUSHOFF,2) = ZERO
                  BUSHOFF(NBUSHOFF,3) = ZERO
               ENDIF
            ELSE IF (OCID >= 0) THEN
               NBUSHOFF = NBUSHOFF + 1
               EDAT(NEDAT_START+9) = NBUSHOFF              ! Slot 9 in EDAT is for the offset key
               DO J=1,3
                  CALL R8FLD ( JCARD(J+3), JF(J+3), R8INP )
                  IF (IERRFL(J+3) == 'N') THEN
                     BUSHOFF(NBUSHOFF,J  ) = R8INP
                  ENDIF
               ENDDO
            ENDIF
            CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,0,0,0 )
            CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,7,8,9 )
            CALL CRDERR ( CARD )

         ELSE

            EDAT(NEDAT_START+8) = -1                       ! Con't entry was blank so default OCID and null offset flag
            EDAT(NEDAT_START+7) =  0

            NBUSHOFF = NBUSHOFF + 1
            EDAT(NEDAT_START+9) = NBUSHOFF                 ! Slot 9 in EDAT is for the offset key
            BUSHOFF(NBUSHOFF,1) = HALF
            BUSHOFF(NBUSHOFF,2) = ZERO
            BUSHOFF(NBUSHOFF,3) = ZERO

         ENDIF

      ELSE

         EDAT(NEDAT_START+8) = -1                          ! There was no con't entry so default OCID
         NBUSHOFF = NBUSHOFF + 1
         EDAT(NEDAT_START+9) = NBUSHOFF                    ! Slot 9 in EDAT is for the offset key
         BUSHOFF(NBUSHOFF,1) = HALF
         BUSHOFF(NBUSHOFF,2) = ZERO
         BUSHOFF(NBUSHOFF,3) = ZERO

      ENDIF

! Write warnings and errors if any

! Update NEDAT

      NEDAT = NEDAT_START + MEDAT_CBUSH



      RETURN

! **********************************************************************************************************************************
 1002 FORMAT(' *WARNING    : ',A,A,' HAS V VEC DEFINED TWO WAYS: FIELDS 6,7,AND/OR 8 AND ALSO FIELD 9. CID IN FIELD 9 WILL BE USED')

 1132 FORMAT(' *ERROR  1132: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY V VECTORS. LIMIT IS ',I8)

 1186 FORMAT(' *ERROR  1186: ERROR IN SPECIFYING V VECTOR ON ',A,A,'. EITHER FIELD 6 MUST BE A POSITIVE INTEGER GRID POINT'        &
                    ,/,14X,' OR FIELDS 6, 7, 8 MUST CONTAIN REAL VECTOR COMPONENTS (WITH DECIMAL POINTS)')

 1180 FORMAT(' *ERROR  1180: INVALID OCID = ',A,' ON ',A,A)

 1187 FORMAT(' *ERROR  1187: GRID SPECIFYING V VECTOR ON ',A,A,' MUST BE > 0. VALUE IS = ',I8)

 1188 FORMAT(' *ERROR  1188: NO V VECTOR SPECIFIED FOR ',A,' ELEMENT ID = ',A)

 1189 FORMAT(' *ERROR  1189: ',A,A,' MUST HAVE NON-NEGATIVE VALUE FOR CID IF FIELD ',I2,' IS NOT BLANK. VALUE READ WAS ',I8)




! **********************************************************************************************************************************

      END SUBROUTINE BD_CBUSH


      SUBROUTINE BD_CBUSH0 ( CARD, LARGE_FLD_INP )

! Processes CBUSH Bulk Data Cards to increment LVVEC and LBUSHOFF if the CBUSH entry has a V vector or offsets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, LBUSHOFF, LVVEC
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CBUSH0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER( 1*BYTE)              :: FND_VVEC          ! Indicator of whether there is an actual V vec on this CBUSH card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: J
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! See if there is an actual V vector (not a grid point). If so, increment LVVEC

      FND_VVEC = 'N'
      DO J=1,JCARD_LEN
         IF ((JCARD(6)(J:J) == '.') .OR. (JCARD(7)(J:J) == '.') .OR. (JCARD(8)(J:J) == '.')) THEN
            FND_VVEC = 'Y'
            EXIT
         ENDIF
      ENDDO
      IF (FND_VVEC == 'Y') THEN
         LVVEC = LVVEC + 1
      ENDIF

! Optional Second Card - see if there are any offsets. If so, increment LBUSHOFF

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC0  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      LBUSHOFF = LBUSHOFF + 1                              ! Even if no cont entry need one BUSHOFF per CBUSH since S = 0.5 default
      IF (ICONT == 1) THEN
         IF (JCARD(2)(1:) /= ' ') THEN
            LBUSHOFF = LBUSHOFF + 1
         ENDIF
         IF ((JCARD(4)(1:) /= ' ') .OR. (JCARD(5)(1:) == ' ') .OR. (JCARD(6)(1:) == ' ')) THEN
            LBUSHOFF = LBUSHOFF + 1                        ! This adds one more that may not be needed
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CBUSH0


      SUBROUTINE BD_PBUSH ( CARD, LARGE_FLD_INP )

! Processes PBUSH Bulk Data Cards. Reads and checks:

!  1) Prop ID and Material ID and enter into array PBUSH
!  2) Area, moments of inertia, torsional constant ans nonstructural mass and enter into array RPBUSH
!  3) From 1st continuation card (if present): coords of 4 points for stress recovery and enter into array RPBUSH
!  4) From 2nd continuation card (if present): area factors for transverse shear and I12 and enter into array RPBUSH

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE PARAMS, ONLY                :  EPSIL
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPBUSH, NPBUSH, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PBUSH, RPBUSH

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, LEFT_ADJ_BDFLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   =   'BD_PBUSH'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(1*BYTE)               :: FOUND_RCV = 'N'   ! 'Y' if RCV continuation entry is present
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: I4INP             ! An integer value read from a field on this BD entry
      INTEGER(LONG)                   :: ICONT       = 0   ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR        = 0   ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: NUM_ENTRIES       ! Num of quantities to read depending on field 3 of parent or cont entry
      INTEGER(LONG)                   :: OFFSET            ! Array index offset
      INTEGER(LONG)                   :: PROPERTY_ID = 0   ! Property ID (field 2 of this property card)


      REAL(DOUBLE)                    :: R8INP             ! A real value read from a field on this BD entry



! **********************************************************************************************************************************
! PBUSH Bulk Data Card routine

!   FIELD   ITEM                                            ARRAY ELEMENT
!   -----   ----                                            -------------
!    2      Property ID                                     PBUSH(npbush,1)
!    3      "K", "B", 'GE" or "RCV"
!   If field 3 = "K"  :
!    4-9    Ki                                             RPBUSH(npbush,1- 6)
!   If field 3 = "B"  :
!    4-9    Bi                                             RPBUSH(npbush,7-12)
!   If field 3 = "GE" :
!    4-9    GEi                                            RPBUSH(npbush,13)
!   If field 3 = "RCV":
!    4-9    RCVi                                           RPBUSH(npbush,14-17)

! Cont entries have the vals for "K", "B", "GE", "RCV" that are not on the 1st entry

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      NPBUSH = NPBUSH+1

! Read and check data on parent card

      CALL I4FLD ( JCARD(2), JF(2), I4INP )                ! Read property ID and enter into array PBUSH
      IF (IERRFL(2) == 'N') THEN
         PROPERTY_ID = I4INP
         DO J=1,NPBUSH-1
            IF (PROPERTY_ID == PBUSH(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),PROPERTY_ID
               WRITE(F06,1145) JCARD(1),PROPERTY_ID
               EXIT
             ENDIF
         ENDDO
         PBUSH(NPBUSH,1) = PROPERTY_ID
      ENDIF

      CALL LEFT_ADJ_BDFLD ( JCARD(3) )                     ! Determine which of the inputs are on the parent entry
      OFFSET      = 0
      NUM_ENTRIES = 6
      IF      (JCARD(3)(1:3) == 'K  ') THEN
         OFFSET      = 0
         NUM_ENTRIES = 6
         CALL BD_IMBEDDED_BLANK   ( JCARD,2,0,4,5,6,7,8,9 )
      ELSE IF (JCARD(3)(1:3) == 'B  ') THEN
         OFFSET      = 6
         NUM_ENTRIES = 6
         CALL BD_IMBEDDED_BLANK   ( JCARD,2,0,4,5,6,7,8,9 )
      ELSE IF (JCARD(3)(1:3) == 'GE ') THEN
         OFFSET      = 12
         NUM_ENTRIES = 6
         CALL BD_IMBEDDED_BLANK   ( JCARD,2,0,4,0,0,0,0,0 )
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,5,6,7,8,9 )
      ELSE IF (JCARD(3)(1:3) == 'RCV') THEN
         OFFSET      = 18
         NUM_ENTRIES = 4
         CALL BD_IMBEDDED_BLANK   ( JCARD,2,0,4,5,6,7,0,0 )
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1178) JCARD(1), JCARD(2), JCARD(3)
         WRITE(F06,1178) JCARD(1), JCARD(2), JCARD(3)
      ENDIF

      DO J=1,NUM_ENTRIES                                   ! Read real property values in fields 4-9
         CALL R8FLD ( JCARD(J+3), JF(J+3), R8INP )
         IF (IERRFL(J+3) == 'N') THEN
            RPBUSH(NPBUSH,J+OFFSET) = R8INP
         ENDIF
      ENDDO

      CALL CRDERR ( CARD )

! Read and check data on 3 optional con't entries

      FOUND_RCV = 'N'
pcont:DO I=1,4

         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         IF (ICONT == 1) THEN

            CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
            CALL LEFT_ADJ_BDFLD ( JCARD(3) )               ! Determine whaich of the inputs are on the parent entry

            OFFSET      = 0
            NUM_ENTRIES = 0
            IF      (JCARD(3)(1:3) == 'K  ') THEN
               OFFSET      = 0
               NUM_ENTRIES = 6
               CALL BD_IMBEDDED_BLANK   ( JCARD,2,0,4,5,6,7,8,9 )
            ELSE IF (JCARD(3)(1:3) == 'B  ') THEN
               OFFSET      = 6
               NUM_ENTRIES = 6
               CALL BD_IMBEDDED_BLANK   ( JCARD,2,0,4,5,6,7,8,9 )
            ELSE IF (JCARD(3)(1:3) == 'GE ') THEN
               OFFSET      = 12
               NUM_ENTRIES = 6
               CALL BD_IMBEDDED_BLANK   ( JCARD,2,0,4,0,0,0,0,0 )
               CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,5,6,7,8,9 )
            ELSE IF (JCARD(3)(1:3) == 'RCV') THEN
               FOUND_RCV = 'Y'
               OFFSET      = 18
               NUM_ENTRIES = 4
               CALL BD_IMBEDDED_BLANK   ( JCARD,2,0,4,5,6,7,0,0 )
               CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )
            ELSE IF (JCARD(3)(1:3) == 'M  ') THEN
               OFFSET      = 22
               NUM_ENTRIES = 1
               CALL BD_IMBEDDED_BLANK   ( JCARD,2,0,4,5,6,7,0,0 )
               CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1148) PROPERTY_ID, 'K, B, GE, or RCV', JCARD(3)
               WRITE(F06,1148) PROPERTY_ID, 'K, B, GE, or RCV', JCARD(3)
               EXIT
            ENDIF
            DO J=1,NUM_ENTRIES                             ! Read real property values in fields 4-9
               CALL R8FLD ( JCARD(J+3), JF(J+3), R8INP )
               IF (IERRFL(J) == 'N') THEN
                  RPBUSH(NPBUSH,J+OFFSET) = R8INP
               ENDIF
            ENDDO

            CALL CRDERR ( CARD )

         ELSE

            EXIT pcont

         ENDIF

      ENDDO pcont

      IF (FOUND_RCV == 'N') THEN
         DO J=19,22                                     ! If no RCV con't entry, default stress/strain RCV's to 1.0
            RPBUSH(NPBUSH,J) = ONE
         ENDDO
      ENDIF




      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1148 FORMAT(' *ERROR  1148: FIELD 3 OF PBUSH ', I8, ' CONTINUATION ENTRY SHOULD BE ', A, ' BUT IS: ', A)

 1178 FORMAT(' *ERROR  1178: ',A,A,' HAS INCORRECT VALUE "',A,'"  IN FIELD 3')

! **********************************************************************************************************************************

      END SUBROUTINE BD_PBUSH


      SUBROUTINE BD_CMASS1 ( CARD )

! Processes CMASS1 Bulk Data Cards. NOTE: MYSTRAN scalar masses must be attached to only 1 point

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NCMASS
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CMASS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CMASS1'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: CMASS_ELID        ! Element ID
      INTEGER(LONG)                   :: GPOINT1,GPOINT2   ! 2 grid points (1 must be blank or zero)
      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
! CMASS1 scalar spring element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Element ID     CMASS(ncmass,1)
!  none     Type (1,2,3,4) CMASS(ncmass,2)
!    3      Prop ID, PID   CMASS(ncmass,3)
!    4      Grid-A         CMASS(ncmass,4)
!    5      Comp-A         CMASS(ncmass,5)
!    6      Grid-B         CMASS(ncmass,6)
!    7      Comp-B         CMASS(ncmass,7)


      NCMASS = NCMASS + 1

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Get element ID and check for duplicate

      CALL I4FLD ( JCARD(2), JF(2), CMASS_ELID )
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NCMASS-1
            IF (CMASS_ELID == CMASS(I,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),CMASS_ELID
               WRITE(F06,1145) JCARD(1),CMASS_ELID
               EXIT
            ENDIF
         ENDDO
         CMASS(NCMASS,1) = CMASS_ELID
      ENDIF

! Get type and put into CMASS col 2

      IF      (JCARD(1)(1:6) == 'CMASS1') THEN
         CMASS(NCMASS,2) = 1
      ELSE IF (JCARD(1)(1:6) == 'CMASS2') THEN
         CMASS(NCMASS,2) = 2
      ELSE IF (JCARD(1)(1:6) == 'CMASS3') THEN
         CMASS(NCMASS,2) = 3
      ELSE IF (JCARD(1)(1:6) == 'CMASS4') THEN
         CMASS(NCMASS,2) = 4
      ENDIF

! Get prop ID and the 2 grids and components

      DO I=3,7
         CALL I4FLD ( JCARD(I), JF(I), CMASS(NCMASS,I) )
      ENDDO

      GPOINT1 = 0
      IF (IERRFL(4) == 'N') THEN
         GPOINT1 = CMASS(NCMASS,4)
      ENDIF

      GPOINT2 = 0
      IF (IERRFL(6) == 'N') THEN
         GPOINT2 = CMASS(NCMASS,6)
      ENDIF

! Special case for MYSTRAN: 1 of the points must be 0 (so that this scalar mass defines only a 1x1 matrix, NOT 2x2)

      IF ((GPOINT1 /= 0) .AND. (GPOINT2 /= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1138) JCARD(1), JCARD(2), GPOINT1, GPOINT2
         WRITE(F06,1138) JCARD(1), JCARD(2), GPOINT1, GPOINT2
      ENDIF

      IF ((GPOINT1 <= 0) .AND. (GPOINT2 <= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1138) JCARD(1), JCARD(2), GPOINT1, GPOINT2
         WRITE(F06,1138) JCARD(1), JCARD(2), GPOINT1, GPOINT2
      ENDIF

! Issue warning if fields 8, 9 are not blank

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,0,0 )   ! Make sure that there are no imbedded blanks in fields 2-7
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )   ! Issue warning if fields 8, 9 are not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1138 FORMAT(' *ERROR  1138: ',2A,' MUST HAVE 1 AND ONLY 1 GRID POINT DEFINED BUT ENTRY HAS: ',I8,' AND ',I8,                      &
                                  ' (SPECIAL CASE FOR MYSTRAN)')

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_CMASS1


      SUBROUTINE BD_CMASS2 ( CARD )

! Processes CMASS2 Bulk Data Cards. NOTE: MYSTRAN scalar masses must be attached to only 1 point

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NCMASS, NPMASS
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CMASS, PMASS, RPMASS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CMASS2'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: CMASS_ELID        ! Element ID
      INTEGER(LONG)                   :: GPOINT1,GPOINT2   ! 2 grid points (1 must be blank or zero)
      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
! CMASS2 scalar spring element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Element ID     CMASS(ncmass,1)
!  none     Type (1,2,3,4) CMASS(ncmass,2)
!    3      Mass, M        RPMASS(npmass,1)
!    4      Grid-A         CMASS(ncmass,4)
!    5      Comp-A         CMASS(ncmass,5)
!    6      Grid-B         CMASS(ncmass,6)
!    7      Comp-B         CMASS(ncmass,7)
!  none     Prop ID, PID   PMASS(npmass,1) (created: PID = -EID)

      NCMASS = NCMASS + 1
      NPMASS = NPMASS + 1

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Get element ID and check for duplicate. Set prop ID = -elem ID and put mass in field 3 into RPMASS

      CALL I4FLD ( JCARD(2), JF(2), CMASS_ELID )
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NCMASS-1
            IF (CMASS_ELID == CMASS(I,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),CMASS_ELID
               WRITE(F06,1145) JCARD(1),CMASS_ELID
               EXIT
            ENDIF
         ENDDO
         CMASS(NCMASS,1) =  CMASS_ELID                      ! CMASS2 elem ID
         CMASS(NCMASS,3) = -CMASS_ELID                      ! Prop ID of CMASS2 set = -elem ID
         PMASS(NPMASS,1) = -CMASS_ELID                      ! Mass value
         CALL R8FLD ( JCARD(3), JF(3), RPMASS(NPMASS,1) )
      ENDIF

! Get type and put into CMASS col 2

      IF      (JCARD(1)(1:6) == 'CMASS1') THEN
         CMASS(NCMASS,2) = 1
      ELSE IF (JCARD(1)(1:6) == 'CMASS2') THEN
         CMASS(NCMASS,2) = 2
      ELSE IF (JCARD(1)(1:6) == 'CMASS3') THEN
         CMASS(NCMASS,2) = 3
      ELSE IF (JCARD(1)(1:6) == 'CMASS4') THEN
         CMASS(NCMASS,2) = 4
      ENDIF

! Get the 2 grids and components

      DO I=4,7
         CALL I4FLD ( JCARD(I), JF(I), CMASS(NCMASS,I) )
      ENDDO

      GPOINT1 = 0
      IF (IERRFL(4) == 'N') THEN
         GPOINT1 = CMASS(NCMASS,4)
      ENDIF

      GPOINT2 = 0
      IF (IERRFL(6) == 'N') THEN
         GPOINT2 = CMASS(NCMASS,6)
      ENDIF

! Special case for MYSTRAN: 1 of the points must be 0 (so that this scalar mass defines only a 1x1 matrix, NOT 2x2)

      IF ((GPOINT1 /= 0) .AND. (GPOINT2 /= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1138) JCARD(1), JCARD(2), GPOINT1, GPOINT2
         WRITE(F06,1138) JCARD(1), JCARD(2), GPOINT1, GPOINT2
      ENDIF

      IF ((GPOINT1 <= 0) .AND. (GPOINT2 <= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1138) JCARD(1), JCARD(2), GPOINT1, GPOINT2
         WRITE(F06,1138) JCARD(1), JCARD(2), GPOINT1, GPOINT2
      ENDIF

! Issue warning if fields 8, 9 are not blank

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,0,0 )   ! Make sure that there are no imbedded blanks in fields 2-7
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )   ! Issue warning if fields 8, 9 are not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1138 FORMAT(' *ERROR  1138: ',2A,' MUST HAVE 1 AND ONLY 1 GRID POINT DEFINED BUT ENTRY HAS: ',I8,' AND ',I8,                      &
                                  ' (SPECIAL CASE FOR MYSTRAN)')

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_CMASS2


      SUBROUTINE BD_CMASS3 ( CARD )

! Processes CMASS3 Bulk Data Cards. NOTE: MYSTRAN scalar masses must be attached to only 1 scalar point

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NCMASS
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CMASS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CMASS3'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: CMASS_ELID        ! Element ID
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SPOINT1,SPOINT2   ! 2 scalar points (1 must be blank or zero)




! **********************************************************************************************************************************
! CMASS3 scalar spring element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Element ID     CMASS(ncmass,1)
!  none     Type (1,2,3,4) CMASS(ncmass,2)
!    3      Prop ID, PID   CMASS(ncmass,3)
!    4      Grid-A         CMASS(ncmass,4)
!    5      Grid-B         CMASS(ncmass,6)


      NCMASS = NCMASS + 1

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Get element ID and check for duplicate

      CALL I4FLD ( JCARD(2), JF(2), CMASS_ELID )
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NCMASS-1
            IF (CMASS_ELID == CMASS(I,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),CMASS_ELID
               WRITE(F06,1145) JCARD(1),CMASS_ELID
               EXIT
            ENDIF
         ENDDO
         CMASS(NCMASS,1) = CMASS_ELID
      ENDIF

! Get type and put into CNASS col 2

      IF      (JCARD(1)(1:6) == 'CMASS1') THEN
         CMASS(NCMASS,2) = 1
      ELSE IF (JCARD(1)(1:6) == 'CMASS2') THEN
         CMASS(NCMASS,2) = 2
      ELSE IF (JCARD(1)(1:6) == 'CMASS3') THEN
         CMASS(NCMASS,2) = 3
      ELSE IF (JCARD(1)(1:6) == 'CMASS4') THEN
         CMASS(NCMASS,2) = 4
      ENDIF

! Get prop ID and the 2 grids. Set components to 1

      CALL I4FLD ( JCARD(3), JF(3), CMASS(NCMASS,3) )

      SPOINT1 = 0
      CALL I4FLD ( JCARD(4), JF(4), CMASS(NCMASS,4) )
      IF (IERRFL(4) == 'N') THEN
         SPOINT1 = CMASS(NCMASS,4)
      ENDIF

      SPOINT2 = 0
      CALL I4FLD ( JCARD(5), JF(5), CMASS(NCMASS,6) )
      IF (IERRFL(5) == 'N') THEN
         SPOINT2 = CMASS(NCMASS,6)
      ENDIF

      CMASS(NCMASS,5) = 1
      CMASS(NCMASS,7) = 1

! Special case for MYSTRAN: 1 of the scalar points must be 0 (so that this scalar mass defines only a 1x1 matrix, NOT 2x2)

      IF ((SPOINT1 /= 0) .AND. (SPOINT2 /= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1138) JCARD(1), JCARD(2), SPOINT1, SPOINT2
         WRITE(F06,1138) JCARD(1), JCARD(2), SPOINT1, SPOINT2
      ENDIF

! Make sure that at least 1 SPOINT was defined

      IF ((SPOINT1 <= 0) .AND. (SPOINT2 <= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1138) JCARD(1), JCARD(2), SPOINT1, SPOINT2
         WRITE(F06,1138) JCARD(1), JCARD(2), SPOINT1, SPOINT2
      ENDIF

! Issue warning if fields 6-9 are not blank

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,0,0,0,0 )   ! Make sure that there are no imbedded blanks in fields 2-5
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )   ! Issue warning if fields 6-9 are not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1138 FORMAT(' *ERROR  1138: ',2A,' MUST HAVE 1 AND ONLY 1 SCALAR POINT DEFINED BUT ENTRY HAS: ',I8,' AND ',I8,                    &
                                  ' (SPECIAL CASE FOR MYSTRAN)')

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_CMASS3


      SUBROUTINE BD_CMASS4 ( CARD )

! Processes CMASS4 Bulk Data Cards. NOTE: MYSTRAN scalar masses must be attached to only 1 scalar point

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NCMASS, NPMASS
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CMASS, PMASS, RPMASS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CMASS4'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: CMASS_ELID        ! Element ID
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: SPOINT1,SPOINT2   ! 2 scalar points (1 must be blank or zero)




! **********************************************************************************************************************************
! CMASS4 scalar spring element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Element ID     CMASS(ncmass,1)
!  none     Type (1,2,3,4) CMASS(ncmass,2)
!    3      Mass, M        RPMASS(npmass,1)
!    4      Grid-A         CMASS(ncmass,4)
!    6      Grid-B         CMASS(ncmass,6)
!  none     Prop ID, PID   PMASS(npmass,1) (created: PID = -EID)

      NCMASS = NCMASS + 1
      NPMASS = NPMASS + 1

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Get element ID and check for duplicate. Set prop ID = -elem ID and put mass in field 3 into RPMASS

      CALL I4FLD ( JCARD(2), JF(2), CMASS_ELID )
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NCMASS-1
            IF (CMASS_ELID == CMASS(I,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),CMASS_ELID
               WRITE(F06,1145) JCARD(1),CMASS_ELID
               EXIT
            ENDIF
         ENDDO
         CMASS(NCMASS,1) =  CMASS_ELID                      ! CMASS4 elem ID
         CMASS(NCMASS,3) = -CMASS_ELID                      ! Prop ID of CMASS4 set = -elem ID
         PMASS(NPMASS,1) = -CMASS_ELID                      ! Mass value
         CALL R8FLD ( JCARD(3), JF(3), RPMASS(NPMASS,1) )
      ENDIF

! Get type and put into CMASS col 2

      IF      (JCARD(1)(1:6) == 'CMASS1') THEN
         CMASS(NCMASS,2) = 1
      ELSE IF (JCARD(1)(1:6) == 'CMASS2') THEN
         CMASS(NCMASS,2) = 2
      ELSE IF (JCARD(1)(1:6) == 'CMASS3') THEN
         CMASS(NCMASS,2) = 3
      ELSE IF (JCARD(1)(1:6) == 'CMASS4') THEN
         CMASS(NCMASS,2) = 4
      ENDIF

! Get 2 grids. Set components to 1

      SPOINT1 = 0
      CALL I4FLD ( JCARD(4), JF(4), CMASS(NCMASS,4) )
      IF (IERRFL(4) == 'N') THEN
         SPOINT1 = CMASS(NCMASS,4)
      ENDIF

      SPOINT2 = 0
      CALL I4FLD ( JCARD(5), JF(5), CMASS(NCMASS,6) )
      IF (IERRFL(5) == 'N') THEN
         SPOINT2 = CMASS(NCMASS,6)
      ENDIF

      CMASS(NCMASS,5) = 1
      CMASS(NCMASS,7) = 1

! Special case for MYSTRAN: 1 of the scalar points must be 0 (so that this scalar mass defines only a 1x1 matrix, NOT 2x2)

      IF ((SPOINT1 /= 0) .AND. (SPOINT2 /= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1138) JCARD(1), JCARD(2), SPOINT1, SPOINT2
         WRITE(F06,1138) JCARD(1), JCARD(2), SPOINT1, SPOINT2
      ENDIF

! Make sure that at least 1 SPOINT was defined

      IF ((SPOINT1 <= 0) .AND. (SPOINT2 <= 0)) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1138) JCARD(1), JCARD(2), SPOINT1, SPOINT2
         WRITE(F06,1138) JCARD(1), JCARD(2), SPOINT1, SPOINT2
      ENDIF

! Issue warning if fields 8, 9 are not blank

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,0,0 )   ! Make sure that there are no imbedded blanks in fields 2-7
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )   ! Issue warning if fields 8, 9 are not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1138 FORMAT(' *ERROR  1138: ',2A,' MUST HAVE 1 AND ONLY 1 SCALAR POINT DEFINED BUT ENTRY HAS: ',I8,' AND ',I8,                    &
                                  ' (SPECIAL CASE FOR MYSTRAN)')

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_CMASS4


      SUBROUTINE BD_PMASS ( CARD )

! Processes PMASS Bulk Data Cards. Read and check

!  1) Property ID and enter into array PMASS
!  2) Mass value and enter into array RPMASS

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, NPMASS
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PMASS, RPMASS

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PMASS'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: PMASS_PID         ! Prop number from field 2,4,6,8




! **********************************************************************************************************************************
! PMASS Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Prop ID         PMASS(npelas,1)
!    3      Mass, MI       RPMASS(npelas,1)
!    4      Prop ID         PMASS(npelas,2)
!    5      Mass, MI       RPMASS(npelas,2)
!    6      Prop ID         PMASS(npelas,3)
!    7      Mass, MI       RPMASS(npelas,3)
!    8      Prop ID         PMASS(npelas,4)
!    9      Mass, MI       RPMASS(npelas,4)


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Read and check data

      DO I=1,4

         IF (JCARD(2*I)(1:) /= ' ') THEN
            NPMASS = NPMASS + 1

            CALL I4FLD ( JCARD(2*I), JF(2*I), PMASS_PID )
            IF (IERRFL(2*I) == 'N') THEN
               DO J=1,NPMASS-1
                  IF (PMASS_PID == PMASS(J,1)) THEN
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1145) JCARD(1),PMASS_PID
                     WRITE(F06,1145) JCARD(1),PMASS_PID
                     EXIT
                  ENDIF
               ENDDO
               PMASS(NPMASS,1) = PMASS_PID
            ENDIF
            CALL R8FLD ( JCARD(2*I+1), JF(2*I+1), RPMASS(NPMASS,1) )
         ENDIF

      ENDDO

      IF (JCARD(2)(1:) /= ' ') THEN
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,0,0,0,0,0,0 )
      ENDIF

      IF (JCARD(4)(1:) /= ' ') THEN
         CALL BD_IMBEDDED_BLANK ( JCARD,0,0,4,5,0,0,0,0 )
      ENDIF

      IF (JCARD(6)(1:) /= ' ') THEN
         CALL BD_IMBEDDED_BLANK ( JCARD,0,0,0,0,6,7,0,0 )
      ENDIF

      IF (JCARD(8)(1:) /= ' ') THEN
         CALL BD_IMBEDDED_BLANK ( JCARD,0,0,0,0,0,0,8,9 )
      ENDIF

      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PMASS


      SUBROUTINE BD_CONM2 ( CARD, LARGE_FLD_INP )

! Processes CONM2 Bulk Data Cards
!  1) Reads CONM2 ID, grid ID and coord system ID and puts them into array CONM2
!  2) Reads mass, offsets, and moments of inertia and puts them into array RCONM2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, LCONM2, NCONM2, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  CONM2, RCONM2

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CONM2'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN=LEN(CARD))        :: CARDP             ! Parent card
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_ID          ! The COMN2 ID
      CHARACTER(16*BYTE)              :: NAME              ! Name for output error purposes

      INTEGER(LONG)                   :: CONM2_ID  = 0     ! The ID for this CONM2 (field 2)
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! CONM2 Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    2      Element ID     CONM2(nconm2,1)
!    3      Grid ID        CONM2(nconm2,2)
!    4      Cord. ID       CONM2(nconm2,3)
!    5      Mass           RCONM2(nconm2,1)
!    6      D1 offset      RCONM2(nconm2,2)
!    7      D2 offset      RCONM2(nconm2,3)
!    8      D3 offset      RCONM2(nconm2,4)
! on optional second card:
!    2      I11 moi-11     RCONM2(nconm2,5)
!    3      I21 moi-21     RCONM2(nconm2,6)
!    4      I22 moi-22     RCONM2(nconm2,7)
!    5      I31 moi-31     RCONM2(nconm2,8)
!    6      I32 moi-32     RCONM2(nconm2,9)
!    7      I33 moi-33     RCONM2(nconm2,10)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      CARDP = CARD
      JCARD_ID = JCARD(2)

! Check for overflow

      NCONM2 = NCONM2+1

! Read data

      CALL I4FLD ( JCARD(2), JF(2), CONM2_ID )             ! Check for duplicate ID number
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NCONM2-1
            IF (CONM2_ID == CONM2(I,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),CONM2_ID
               WRITE(F06,1145) JCARD(1),CONM2_ID
               EXIT
            ENDIF
         ENDDO
         CONM2(NCONM2,1) = CONM2_ID
      ENDIF

      CALL I4FLD ( JCARD(3), JF(3),  CONM2(NCONM2,2) )
      CALL I4FLD ( JCARD(4), JF(4),  CONM2(NCONM2,3) )
      CALL R8FLD ( JCARD(5), JF(5), RCONM2(NCONM2,1) )
      CALL R8FLD ( JCARD(6), JF(6), RCONM2(NCONM2,2) )
      CALL R8FLD ( JCARD(7), JF(7), RCONM2(NCONM2,3) )
      CALL R8FLD ( JCARD(8), JF(8), RCONM2(NCONM2,4) )

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,0 )     ! Make sure that there are no imbedded blanks in fields 2-8
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )   ! Issue warning if field 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Optional Second Card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      RCONM2(NCONM2, 5) = 0.0
      RCONM2(NCONM2, 6) = 0.0
      RCONM2(NCONM2, 7) = 0.0
      RCONM2(NCONM2, 8) = 0.0
      RCONM2(NCONM2, 9) = 0.0
      RCONM2(NCONM2,10) = 0.0

      IF (ICONT == 1) THEN

         IF (CARD(1:) /= ' ') THEN                        ! Only process continuation entries if 1st one is not totally blank

            CALL R8FLD ( JCARD(2), JF(2), RCONM2(NCONM2, 5) )
            CALL R8FLD ( JCARD(3), JF(3), RCONM2(NCONM2, 6) )
            CALL R8FLD ( JCARD(4), JF(4), RCONM2(NCONM2, 7) )
            CALL R8FLD ( JCARD(5), JF(5), RCONM2(NCONM2, 8) )
            CALL R8FLD ( JCARD(6), JF(6), RCONM2(NCONM2, 9) )
            CALL R8FLD ( JCARD(7), JF(7), RCONM2(NCONM2,10) )

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,0,0 )
            CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9)
            CALL CRDERR ( CARD )

         ENDIF

      ENDIF

! Check for negative mass/MOI's

      IF((RCONM2(NCONM2, 1) < ZERO) .OR. (RCONM2(NCONM2, 1) < ZERO) .OR.                                                           &
         (RCONM2(NCONM2, 1) < ZERO) .OR. (RCONM2(NCONM2, 1) < ZERO)) THEN
         WRITE(ERR,101) CARDP
         IF (ECHO == 'NONE  ') THEN
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) CARDP
               WRITE(F06,101) CARD
            ENDIF
         ENDIF
      ENDIF

      IF (RCONM2(NCONM2, 1) < ZERO) THEN
         WARN_ERR = WARN_ERR + 1
         NAME  = 'MASS            '
         WRITE(ERR,1134) JCARD_ID,NAME,RCONM2(NCONM2, 1)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1134) JCARD_ID,NAME,RCONM2(NCONM2, 1)
         ENDIF
      ENDIF

      IF (RCONM2(NCONM2, 5) < ZERO) THEN
         WARN_ERR = WARN_ERR + 1
         NAME = 'PRINCIPAL MOI-11'
         WRITE(ERR,1134) JCARD_ID,NAME,RCONM2(NCONM2, 5)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1134) JCARD_ID,NAME,RCONM2(NCONM2, 5)
         ENDIF
      ENDIF

      IF (RCONM2(NCONM2, 7) < ZERO) THEN
         WARN_ERR = WARN_ERR + 1
         NAME = 'PRINCIPAL MOI-22'
         WRITE(ERR,1134) JCARD_ID,NAME,RCONM2(NCONM2, 7)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1134) JCARD_ID,NAME,RCONM2(NCONM2, 7)
         ENDIF
      ENDIF

      IF (RCONM2(NCONM2,10) < ZERO) THEN
         WARN_ERR = WARN_ERR + 1
         NAME = 'PRINCIPAL MOI-33'
         WRITE(ERR,1134) JCARD_ID, NAME, RCONM2(NCONM2,10)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1134) JCARD_ID, NAME, RCONM2(NCONM2,10)
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1134 FORMAT(' *WARNING    : NEGATIVE MASS OR PRINCIPAL MOI VALUE ON CONM2 ',A,': ',A16,' = ',ES14.6)

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_CONM2

   END MODULE SPRING_BUSH_MASS
