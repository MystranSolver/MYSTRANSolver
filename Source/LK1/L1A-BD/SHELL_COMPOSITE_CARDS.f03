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

   MODULE SHELL_COMPOSITE_CARDS

   IMPLICIT NONE

   PRIVATE

   EXTERNAL :: ELEPRO

   PUBLIC :: BD_CQUAD, BD_CQUAD0, BD_CQUAD8, BD_CQUAD80, BD_CTRIA, BD_CTRIA0, BD_CSHEAR, BD_PSHEAR, BD_PSHEL, BD_PCOMP, BD_PCOMP0, BD_PCOMP1, BD_PCOMP10

   CONTAINS

      SUBROUTINE BD_CQUAD ( CARD, LARGE_FLD_INP, NUM_GRD )

! Processes CQUADi, CQDPLTi, CQDMEMi Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, FATAL_ERR, JCARD_LEN, JF, LMATANGLE, LPLATEOFF, LPLATETHICK,        &
                                         MEDAT_CQUAD, NCQUAD4K, NCQUAD4, NEDAT, NELE, NMATANGLE, NPLATEOFF, NPLATETHICK
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE, MATANGLE, PLATEOFF, PLATETHICK

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CQUAD'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! Indicator of the type of val found in a B.D. field (e.g. int, real,...)

      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD           ! Number of GRID's + SPOINT's for the elem
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: I4INP             ! A value read from input file that should be an integer
      INTEGER(LONG)                   :: INT41,INT42        ! An integer used in getting MATANGLE
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein


      REAL(DOUBLE)                    :: R8INP     = ZERO  ! A value read from input file that should be a real value



! **********************************************************************************************************************************
! CQUADi element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele) = Q1, Q2, Q3, QA, QB
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+2)
!    4      Grid A         EDAT(nedat+3)
!    5      Grid B         EDAT(nedat+4)
!    6      Grid C         EDAT(nedat+5)
!    7      Grid D         EDAT(nedat+6)
!    8      Matl angle     NMATANGLE (MATANGLE key) goes in EDAT(nedat+7)
!    9      Offset         NPLATEOFF (PLATEOFF key) goes in EDAT(nedat+8)
!                          EDAT(nedat+9)
! on optional second card:
!   4-7     Ti             Membrane thicknes at grids 1-4. These will go into array PLATETHICK
!                          EDAT(nedat+10) will hold NPLATETHICK the row in PLATETHICK where Ti is located

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Set JCARD_EDAT to JCARD

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Check property ID field. Set to element ID if blank

      IF (JCARD(3)(1:) == ' ') THEN
         JCARD_EDAT(3) = JCARD(2)
      ENDIF

! Read and check data
                                                           ! Load 6 items into EDAT
      CALL ELEPRO ( 'Y', JCARD_EDAT, 6, MEDAT_CQUAD, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'N', 'N' )

      NUM_GRD = 4
      IF       (JCARD(1)(1:7) == 'CQUAD4K') THEN
         NCQUAD4K = NCQUAD4K + 1
         ETYPE(NELE) = 'QUAD4K  '
      ELSE IF ((JCARD(1)(1:7) == 'CQUAD4 ') .OR. (JCARD(1)(1:7) == 'CQUAD4*')) THEN
         NCQUAD4 = NCQUAD4 + 1
         ETYPE(NELE) = 'QUAD4   '
      ENDIF

! Read material property orientation angle. It takes 2 values put into EDAT to cover all of the possibilities of field 8:
!  (a) If field 8 is a real value it is the angle of the material axis relative to the element x axis.
!         (1) the 2 values to put into EDAT are: the value in field 8 (this will be the row in MATANGLE to get the angle), and a 0
!  (b) If field 8 is an integer it means that the angle is identified by a coord system.
!         (2) if field 8 is a positive number the 2 values to put into EDAT are: the neg of field 8 integer and a 2
!         (3) if field 8 is an actual 0 (basic coodr sys) the 2 values to put into EDAT are: the neg of field 8 integer and a 1
!  (c) If field 8 is blank:
!         (4) put into EDAT the 2 values: 0 and 0

      IF (JCARD(8)(1:) /= ' ') THEN                        ! Matl angle field has something in it so there is a angle defined

         TOKEN = JCARD(8)(1:8)                             ! Only send the 1st 8 chars of this JCARD. It has been left justified
         CALL TOKCHK ( TOKEN, TOKTYP )
         IF      ((TOKTYP /= 'INTEGER ') .AND. (TOKTYP /= 'FL PT   ')) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE (ERR,1196) JCARD(1), JCARD(2), JCARD(8)
            WRITE (F06,1196) JCARD(1), JCARD(2), JCARD(8)

         ELSE IF (TOKTYP == 'FL PT   ') THEN               ! -- Set keys based on "FL PT" number in material angle field
            NMATANGLE = NMATANGLE + 1
            IF (NMATANGLE > LMATANGLE) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1141) SUBR_NAME,LMATANGLE
               WRITE(F06,1141) SUBR_NAME,LMATANGLE
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            INT41 = NMATANGLE                              ! ---- Mtrl angle key is row in MATANGLE
            INT42 = 0                                      ! ---- Basic or not key (0 is undefined)
            CALL R8FLD ( JCARD(8), JF(8), R8INP )
            IF (IERRFL(8) == 'N') THEN
               MATANGLE(NMATANGLE) = R8INP
            ENDIF

         ELSE IF (TOKTYP == 'INTEGER ') THEN               ! -- Set keys for "INTEGER" number in material angle field
            CALL I4FLD ( JCARD(8), JF(8), I4INP )
            INT41 = -I4INP                                 ! ---- Mtrl angle key is set to neg of the coord sys ID
            IF      (I4INP > 0) THEN
               INT42 = 2                                   ! -------- If field was > 0 set "Basic or not" key to 2
            ELSE IF (I4INP == 0) THEN
               INT42 = 1                                   ! -------- If field was = 0 set "Basic or not" key to 1
            ELSE IF (I4INP < 0) THEN
               INT42 = 0                                   ! -------- If field was < 0 then error
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1197) JCARD(1), JCARD(2), JF(8), JCARD(8)
               WRITE(F06,1197) JCARD(1), JCARD(2), JF(8), JCARD(8)
            ENDIF

         ENDIF

      ELSE                                                 ! Matl angle field is blank so no angle specified

         INT41 = 0
         INT42 = 0

      ENDIF

      NEDAT = NEDAT + 1
      EDAT(NEDAT) = INT41
      NEDAT = NEDAT + 1
      EDAT(NEDAT) = INT42

! Read plate offset

      IF (JCARD(9)(1:) /= ' ') THEN                        ! 8: Load plate offset key (or 0 if none) into EDAT
         NPLATEOFF = NPLATEOFF + 1
         IF (NPLATEOFF > LPLATEOFF) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1144) SUBR_NAME,' TOO MANY PLATE OFFSETS. LIMIT IS NPLATEOFF =  ',LPLATEOFF
            WRITE(F06,1144) SUBR_NAME,' TOO MANY PLATE OFFSETS. LIMIT IS NPLATEOFF =  ',LPLATEOFF
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         NEDAT = NEDAT + 1
         EDAT(NEDAT) = NPLATEOFF
         CALL R8FLD ( JCARD(9), JF(9), R8INP )
         IF (IERRFL(9) == 'N') THEN
            PLATEOFF(NPLATEOFF) = R8INP
         ENDIF
      ELSE

         NEDAT = NEDAT + 1
         EDAT(NEDAT) = 0

      ENDIF

! Load a 0 into EDAT as a flag for whether this element references a PSHELL or a PCOMP.
! This will be decided in subr ELEM_PROP_MATL_IIDS when the actual PID in EDAT is converted to an internal PID.

      NEDAT = NEDAT + 1                                    ! 8: PSHELL/PCOMP flag (to be set in subr ELEM_PROP_MATL_IIDS)
      EDAT(NEDAT) = 0

! Load a 0 into EDAT as a flag for whether this element has thicknesses defined on a continuation entry

      NEDAT = NEDAT + 1                                    ! 8: PSHELL/PCOMP flag (to be set in subr ELEM_PROP_MATL_IIDS)
      EDAT(NEDAT) = 0

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,8,9 )   ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Optional second card (to define membrane thicknesses as grid values):

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      IF (ICONT == 1) THEN                                 ! Since there is a cont entry we assume Ti are here

         IF (CARD(1:) /= ' ') THEN                         ! Only process continuation entries if 1st one is not totally blank

            EDAT(NEDAT) = NPLATETHICK + 1

            DO J=4,7                                       ! Read 4 thicknesses
               NPLATETHICK = NPLATETHICK + 1
               IF (NPLATETHICK > LPLATETHICK) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1144) SUBR_NAME,' TOO MANY PLATE THICKNESSES. LIMIT IS NPLATETHICK = ',LPLATETHICK
                  WRITE(F06,1144) SUBR_NAME,' TOO MANY PLATE THICKNESSES. LIMIT IS NPLATETHICK = ',LPLATETHICK
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
               CALL R8FLD ( JCARD(J), JF(J), R8INP )
               IF (IERRFL(J) == 'N') THEN
                  PLATETHICK(NPLATETHICK) = R8INP
               ENDIF
            ENDDO

            CALL BD_IMBEDDED_BLANK ( JCARD,0,0,4,5,6,7,0,0 )
            CALL CARD_FLDS_NOT_BLANK ( JCARD,2,3,0,0,0,0,8,9 )
            CALL CRDERR ( CARD )

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1141 FORMAT(' *ERROR  1141: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY PLATE ELEMENT MATERIAL PROPERTY ANGLES. LIMIT IS NMATANGLE =  ',I8)

 1144 FORMAT(' *ERROR  1144: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A,I8)

 1196 FORMAT(' *ERROR  1196: VALUE FOR MATERIAL ANGLE ON ',A,A,' MUST BE AN INTEGER OR REAL NUMBER BUT VALUE READ WAS ',A)

 1197 FORMAT(' *ERROR  1197: FOR ',A,A,' THE COORD SYS ID IN FIELD ',I2,' MUST BE >= 0. HOWEVER, THE VALUE INPUT WAS ',A)

! **********************************************************************************************************************************

      END SUBROUTINE BD_CQUAD


      SUBROUTINE BD_CQUAD0 ( CARD, LARGE_FLD_INP )

! Processes CQUAD Bulk Data Cards to increment:
!   (1) LMATANGLE   if the elem has a material prop angle
!   (2) LPLATEOFF   if the elem has an offset
!   (3) LPLATETHICK if the elem has thicknesses defined on a continuation entry

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, LMATANGLE, LPLATEOFF, LPLATETHICK
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CQUAD0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! See if there is a material orientation angle. If so, increment LMATANGLE

      IF (JCARD(8)(1:) /= ' ') THEN
         LMATANGLE = LMATANGLE + 1
      ENDIF

! See if there is a plate offset. If so, increment LPLATEOFF

      IF (JCARD(9)(1:) /= ' ') THEN
         LPLATEOFF = LPLATEOFF + 1
      ENDIF

! Optional second card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC0  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      IF (ICONT == 1) THEN                                ! Since there is a cont entry we assume 4 Ti are here
         IF (CARD(1:) /= ' ') THEN                        ! Only process continuation entries if 1st one is not totally blank
            LPLATETHICK = LPLATETHICK + 4
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CQUAD0


      SUBROUTINE BD_CQUAD8 ( CARD, LARGE_FLD_INP, NUM_GRD )

! Processes CQUAD8 Bulk Data Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, FATAL_ERR, JCARD_LEN, JF, LMATANGLE, LPLATEOFF, LPLATETHICK,        &
                                         MEDAT_CQUAD8, NCQUAD8, NEDAT, NELE, NMATANGLE, NPLATEOFF, NPLATETHICK
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE, MATANGLE, PLATEOFF, PLATETHICK

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CQUAD8'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! Indicator of the type of val found in a B.D. field (e.g. int, real,...)
      CHARACTER(LEN=JCARD_LEN)        :: ID                ! Character value of element ID (field 2 of parent card)
      CHARACTER(LEN=JCARD_LEN)        :: NAME              ! Field 1 of CARD

      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD           ! Number of GRID's + SPOINT's for the elem
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: I4INP             ! A value read from input file that should be an integer
      INTEGER(LONG)                   :: INT41,INT42        ! An integer used in getting MATANGLE
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein

      REAL(DOUBLE)                    :: R8INP     = ZERO  ! A value read from input file that should be a real value


! **********************************************************************************************************************************
! CQUAD8 element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele)
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+2)
!    4      Grid 1         EDAT(nedat+3)
!    5      Grid 2         EDAT(nedat+4)
!    6      Grid 3         EDAT(nedat+5)
!    7      Grid 4         EDAT(nedat+6)
!    8      Grid 5         EDAT(nedat+7)
!    9      Grid 6         EDAT(nedat+8)
! on required second card:
!    2      Grid 7         EDAT(nedat+9)
!    3      Grid 8         EDAT(nedat+10)
!   4-7     Ti             Membrane thicknes at grids 1-4. These will go into array PLATETHICK
!                          EDAT(nedat+14) will hold NPLATETHICK the row in PLATETHICK where Ti is located
!    8      Matl angle     NMATANGLE (MATANGLE key) goes in EDAT(nedat+11)
!    9      Offset         NPLATEOFF (PLATEOFF key) goes in EDAT(nedat+12)
!                          EDAT(nedat+13)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      NAME = JCARD(1)
      ID = JCARD(2)

! Set JCARD_EDAT to JCARD

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Read and check data
                                                           ! Load 8 items into EDAT
      CALL ELEPRO ( 'Y', JCARD_EDAT, 8, MEDAT_CQUAD8, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y' )

      NUM_GRD = 8
      NCQUAD8 = NCQUAD8 + 1
      ETYPE(NELE) = 'QUAD8   '

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,8,9 )   ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Required 2nd card:

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

        CALL ELEPRO ( 'N', JCARD_EDAT, 2, 2, 'Y', 'Y', 'N', 'N', 'N', 'N', 'N', 'N' )

! Read material property orientation angle. It takes 2 values put into EDAT to cover all of the possibilities of field 8:
!  (a) If field 8 is a real value it is the angle of the material axis relative to the element x axis.
!         (1) the 2 values to put into EDAT are: the value in field 8 (this will be the row in MATANGLE to get the angle), and a 0
!  (b) If field 8 is an integer it means that the angle is identified by a coord system.
!         (2) if field 8 is a positive number the 2 values to put into EDAT are: the neg of field 8 integer and a 2
!         (3) if field 8 is an actual 0 (basic coodr sys) the 2 values to put into EDAT are: the neg of field 8 integer and a 1
!  (c) If field 8 is blank:
!         (4) put into EDAT the 2 values: 0 and 0

        IF (JCARD(8)(1:) /= ' ') THEN                      ! Matl angle field has something in it so there is a angle defined

          TOKEN = JCARD(8)(1:8)                            ! Only send the 1st 8 chars of this JCARD. It has been left justified
          CALL TOKCHK ( TOKEN, TOKTYP )
          IF      ((TOKTYP /= 'INTEGER ') .AND. (TOKTYP /= 'FL PT   ')) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE (ERR,1196) JCARD(1), JCARD(2), JCARD(8)
            WRITE (F06,1196) JCARD(1), JCARD(2), JCARD(8)

          ELSE IF (TOKTYP == 'FL PT   ') THEN              ! -- Set keys based on "FL PT" number in material angle field
            NMATANGLE = NMATANGLE + 1
            IF (NMATANGLE > LMATANGLE) THEN
              FATAL_ERR = FATAL_ERR + 1
              WRITE(ERR,1141) SUBR_NAME,LMATANGLE
              WRITE(F06,1141) SUBR_NAME,LMATANGLE
              CALL OUTA_HERE ( 'Y' )
            ENDIF
            INT41 = NMATANGLE                              ! ---- Mtrl angle key is row in MATANGLE
            INT42 = 0                                      ! ---- Basic or not key (0 is undefined)
            CALL R8FLD ( JCARD(8), JF(8), R8INP )
            IF (IERRFL(8) == 'N') THEN
              MATANGLE(NMATANGLE) = R8INP
            ENDIF

          ELSE IF (TOKTYP == 'INTEGER ') THEN              ! -- Set keys for "INTEGER" number in material angle field
            CALL I4FLD ( JCARD(8), JF(8), I4INP )
            INT41 = -I4INP                                 ! ---- Mtrl angle key is set to neg of the coord sys ID
            IF      (I4INP > 0) THEN
              INT42 = 2                                    ! -------- If field was > 0 set "Basic or not" key to 2
            ELSE IF (I4INP == 0) THEN
              INT42 = 1                                    ! -------- If field was = 0 set "Basic or not" key to 1
            ELSE IF (I4INP < 0) THEN
              INT42 = 0                                    ! -------- If field was < 0 then error
              FATAL_ERR = FATAL_ERR + 1
              WRITE(ERR,1197) JCARD(1), JCARD(2), JF(8), JCARD(8)
              WRITE(F06,1197) JCARD(1), JCARD(2), JF(8), JCARD(8)
            ENDIF

          ENDIF

        ELSE                                               ! Matl angle field is blank so no angle specified

          INT41 = 0
          INT42 = 0

        ENDIF

        NEDAT = NEDAT + 1
        EDAT(NEDAT) = INT41
        NEDAT = NEDAT + 1
        EDAT(NEDAT) = INT42

! Read plate offset

        IF (JCARD(9)(1:) /= ' ') THEN                      ! 12: Load plate offset key (or 0 if none) into EDAT
          NPLATEOFF = NPLATEOFF + 1
          IF (NPLATEOFF > LPLATEOFF) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1144) SUBR_NAME,' TOO MANY PLATE OFFSETS. LIMIT IS LPLATEOFF =  ',LPLATEOFF
            WRITE(F06,1144) SUBR_NAME,' TOO MANY PLATE OFFSETS. LIMIT IS LPLATEOFF =  ',LPLATEOFF
            CALL OUTA_HERE ( 'Y' )
          ENDIF
          NEDAT = NEDAT + 1
          EDAT(NEDAT) = NPLATEOFF
          CALL R8FLD ( JCARD(9), JF(9), R8INP )
          IF (IERRFL(9) == 'N') THEN
            IF(R8INP /= 0.0) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,*) ' *ERROR : ZOFFS NOT ALLOWED ON CQUAD8'
               WRITE(F06,*) ' *ERROR : ZOFFS NOT ALLOWED ON CQUAD8'
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            PLATEOFF(NPLATEOFF) = R8INP
          ENDIF
        ELSE

          NEDAT = NEDAT + 1
          EDAT(NEDAT) = 0

        ENDIF


! Load a 0 into EDAT as a flag for whether this element references a PSHELL or a PCOMP.
! This will be decided in subr ELEM_PROP_MATL_IIDS when the actual PID in EDAT is converted to an internal PID.

        NEDAT = NEDAT + 1                                  ! 13: PSHELL/PCOMP flag (to be set in subr ELEM_PROP_MATL_IIDS)
        EDAT(NEDAT) = 0

! Load a 0 into EDAT as a flag for whether this element has thicknesses defined on a continuation entry

        NEDAT = NEDAT + 1                                  ! 14: PSHELL/PCOMP flag (to be set in subr ELEM_PROP_MATL_IIDS)
        EDAT(NEDAT) = 0


! Optional fields 4-7 (to define membrane thicknesses as grid values):
         IF ((JCARD(4)(1:) /= ' ') .OR. (JCARD(5)(1:) /= ' ') .OR. (JCARD(6)(1:) /= ' ') .OR. (JCARD(7)(1:) /= ' ')) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,*) ' *ERROR : GRID POINT THICKNESSES T1, T2, T3, T4 ARE NOT ALLOWED ON CQUAD8'
            WRITE(F06,*) ' *ERROR : GRID POINT THICKNESSES T1, T2, T3, T4 ARE NOT ALLOWED ON CQUAD8'
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         IF (JCARD(4)(1:) /= ' ') THEN

            EDAT(NEDAT) = NPLATETHICK + 1

            DO J=4,7                                         ! Read 4 thicknesses
               NPLATETHICK = NPLATETHICK + 1
               IF (NPLATETHICK > LPLATETHICK) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1144) SUBR_NAME,' TOO MANY PLATE THICKNESSES. LIMIT IS LPLATETHICK = ',LPLATETHICK
                  WRITE(F06,1144) SUBR_NAME,' TOO MANY PLATE THICKNESSES. LIMIT IS LPLATETHICK = ',LPLATETHICK
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
               CALL R8FLD ( JCARD(J), JF(J), R8INP )
               IF (IERRFL(J) == 'N') THEN
                  PLATETHICK(NPLATETHICK) = R8INP
               ENDIF
            ENDDO

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )

         ELSE

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,0,0,0,0,8,9 )  ! Make sure that there are no imbedded blanks in fields 2,3,8,9
            CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,4,5,6,7,0,0 )! Issue warning if fields 4-7 not blank

         ENDIF

        CALL CRDERR ( CARD )                               ! CRDERR prints errors found when reading fields


      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1136) NAME, ID
         WRITE(F06,1136) NAME, ID

      ENDIF


! **********************************************************************************************************************************

      RETURN

! **********************************************************************************************************************************
 1136 FORMAT(' *ERROR  1136: REQUIRED CONTINUATION FOR ',A,' ID = ',A,' MISSING')

 1141 FORMAT(' *ERROR  1141: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY PLATE ELEMENT MATERIAL PROPERTY ANGLES. LIMIT IS NMATANGLE =  ',I8)

 1144 FORMAT(' *ERROR  1144: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A,I8)

 1196 FORMAT(' *ERROR  1196: VALUE FOR MATERIAL ANGLE ON ',A,A,' MUST BE AN INTEGER OR REAL NUMBER BUT VALUE READ WAS ',A)

 1197 FORMAT(' *ERROR  1197: FOR ',A,A,' THE COORD SYS ID IN FIELD ',I2,' MUST BE >= 0. HOWEVER, THE VALUE INPUT WAS ',A)

! **********************************************************************************************************************************

      END SUBROUTINE BD_CQUAD8


      SUBROUTINE BD_CQUAD80 ( CARD, LARGE_FLD_INP )

! Processes CQUAD8 Bulk Data Cards to increment:
!   (1) LMATANGLE   if the elem has a material prop angle
!   (2) LPLATEOFF   if the elem has an offset
!   (3) LPLATETHICK if the elem has thicknesses defined

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, LMATANGLE, LPLATEOFF, LPLATETHICK

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CQUAD80'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein

! **********************************************************************************************************************************

! Skip the first card
      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC0  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF

      IF (ICONT == 1) THEN                                 ! There was a 1st continuation card (mandatory)

! See if there is a material orientation angle. If so, increment LMATANGLE
        CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
        IF (JCARD(8)(1:) /= ' ') THEN
           LMATANGLE = LMATANGLE + 1
        ENDIF

! See if there is a plate offset. If so, increment LPLATEOFF

        IF (JCARD(9)(1:) /= ' ') THEN
           LPLATEOFF = LPLATEOFF + 1
        ENDIF

        LPLATETHICK = LPLATETHICK + 4

      ENDIF


! **********************************************************************************************************************************

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CQUAD80


      SUBROUTINE BD_CTRIA ( CARD, LARGE_FLD_INP, NUM_GRD )

! Processes CTRIAi, CTRPLTi, CTRMEMi BULK DATA Cards
!  1) Sets ETYPE for this element type
!  2) Calls subr ELEPRO to read element ID, property ID and connection data into array EDAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, FATAL_ERR, JCARD_LEN, JF, LMATANGLE, LPLATEOFF, LPLATETHICK,        &
                                         MEDAT_CTRIA, NCTRIA3K, NCTRIA3, NEDAT, NELE, NMATANGLE, NPLATEOFF, NPLATETHICK
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE, MATANGLE, PLATEOFF, PLATETHICK

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CTRIA'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO
      CHARACTER( 8*BYTE)              :: TOKEN             ! The 1st 8 characters from a JCARD
      CHARACTER( 8*BYTE)              :: TOKTYP            ! Indicator of the type of val found in a B.D. field (e.g. int, real,...)

      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD           ! Number of GRID's + SPOINT's for the elem
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: I4INP             ! A value read from input file that should be an integer
      INTEGER(LONG)                   :: INT41,INT42        ! An integer used in getting MATANGLE
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein


      REAL(DOUBLE)                    :: R8INP     = ZERO  ! A value read from input file that should be a real value



! **********************************************************************************************************************************
! CTRIAi, CTRPLTi, CTRMEM element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele) = T1, T2, T3, TA, TB
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+2)
!    4      Grid A         EDAT(nedat+3)
!    5      Grid B         EDAT(nedat+4)
!    6      Grid C         EDAT(nedat+5)
!    7      Matl angle     NMATANGLE (MATANGLE key) goes in EDAT(nedat+6)
!    8      Offset         NPLATEOFF (PLATEOFF key) goes in EDAT(nedat+7)
!                          EDAT(nedat+8)
! on optional second card:
!   4-7     Ti             Membrane thicknes at grids 1-4. These will go into array PLATETHICK
!                          EDAT(nedat+9) will hold NPLATETHICK the row in PLATETHICK where Ti is located


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Set JCARD_EDAT to JCARD

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Check property ID field. Set to element ID if blank

      IF (JCARD(3)(1:) == ' ') THEN
         JCARD_EDAT(3) = JCARD(2)
      ENDIF

! Read and check data
                                                           ! Load 5 items into EDAT
      CALL ELEPRO ( 'Y', JCARD_EDAT, 5, MEDAT_CTRIA, 'Y', 'Y', 'Y', 'Y', 'Y', 'N', 'N', 'N' )

      NUM_GRD = 3
      IF       (JCARD(1)(1:7) == 'CTRIA3K') THEN
         NCTRIA3K = NCTRIA3K + 1
         ETYPE(NELE) = 'TRIA3K  '
      ELSE IF ((JCARD(1)(1:7) == 'CTRIA3 ') .OR. (JCARD(1)(1:7) == 'CTRIA3*')) THEN
         NCTRIA3 = NCTRIA3 + 1
         ETYPE(NELE) = 'TRIA3   '
      ENDIF

! Read material property orientation angle. It takes 2 values put into EDAT to cover all of the possibilities of field 7:
!  (a) If field 7 is a real value it is the angle of the material axis relative to the element x axis.
!         (1) the 2 values to put into EDAT are: the value in field 7 (this will be the row in MATANGLE to get the angle), and a 0
!  (b) If field 7 is an integer it means that the angle is identified by a coord system.
!         (2) if field 7 is a positive number the 2 values to put into EDAT are: the neg of field 7 integer and a 2
!         (3) if field 7 is an actual 0 (basic coodr sys) the 2 values to put into EDAT are: the neg of field 7 integer and a 1
!  (c) If field 7 is blank:
!         (4) put into EDAT the 2 values: 0 and 0

      IF (JCARD(7)(1:) /= ' ') THEN                        ! Matl angle field has something in it so there is a angle defined

         TOKEN = JCARD(7)(1:8)                             ! Only send the 1st 8 chars of this JCARD. It has been left justified
         CALL TOKCHK ( TOKEN, TOKTYP )
         IF      ((TOKTYP /= 'INTEGER ') .AND. (TOKTYP /= 'FL PT   ')) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE (ERR,1196) JCARD(1), JCARD(2), JCARD(7)
            WRITE (F06,1196) JCARD(1), JCARD(2), JCARD(7)

         ELSE IF (TOKTYP == 'FL PT   ') THEN               ! -- Set keys based on "FL PT" number in material angle field
            NMATANGLE = NMATANGLE + 1
            IF (NMATANGLE > LMATANGLE) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1141) SUBR_NAME,LMATANGLE
               WRITE(F06,1141) SUBR_NAME,LMATANGLE
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            INT41 = NMATANGLE                              ! ---- Mtrl angle key is row in MATANGLE
            INT42 = 0                                      ! ---- Basic or not key (0 is undefined)
            CALL R8FLD ( JCARD(7), JF(7), R8INP )
            IF (IERRFL(7) == 'N') THEN
               MATANGLE(NMATANGLE) = R8INP
            ENDIF

         ELSE IF (TOKTYP == 'INTEGER ') THEN               ! -- Set keys for "INTEGER" number in material angle field
            CALL I4FLD ( JCARD(7), JF(7), I4INP )
            INT41 = -I4INP                                 ! ---- Mtrl angle key is set to neg of the coord sys ID
            IF      (I4INP > 0) THEN
               INT42 = 2                                   ! -------- If field was > 0 set "Basic or not" key to 2
            ELSE IF (I4INP == 0) THEN
               INT42 = 1                                   ! -------- If field was = 0 set "Basic or not" key to 1
            ELSE IF (I4INP < 0) THEN
               INT42 = 0                                   ! -------- If field was < 0 then error
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1197) JCARD(1), JCARD(2), JF(7), JCARD(7)
               WRITE(F06,1197) JCARD(1), JCARD(2), JF(7), JCARD(7)
            ENDIF

         ENDIF

      ELSE                                                 ! Matl angle field is blank so no angle specified

         INT41 = 0
         INT42 = 0

      ENDIF

      NEDAT = NEDAT + 1
      EDAT(NEDAT) = INT41
      NEDAT = NEDAT + 1
      EDAT(NEDAT) = INT42

! Read plate offset

      IF (JCARD(8)(1:) /= ' ') THEN                        ! 7: Load plate offset key (or 0 if none) into EDAT
         NPLATEOFF = NPLATEOFF + 1
         IF (NPLATEOFF > LPLATEOFF) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1144) SUBR_NAME,' TOO MANY PLATE OFFSETS. LIMIT IS NPLATEOFF =  ',LPLATEOFF
            WRITE(F06,1144) SUBR_NAME,' TOO MANY PLATE OFFSETS. LIMIT IS NPLATEOFF =  ',LPLATEOFF
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         NEDAT = NEDAT + 1
         EDAT(NEDAT) = NPLATEOFF
         CALL R8FLD ( JCARD(8), JF(8), R8INP )
         IF (IERRFL(8) == 'N') THEN
            PLATEOFF(NPLATEOFF) = R8INP
         ENDIF
      ELSE

         NEDAT = NEDAT + 1
         EDAT(NEDAT) = 0

      ENDIF

! Load a 0 into EDAT as a flag for whether this element references a PSHELL or a PCOMP.
! This will be decided in subr ELEM_PROP_MATL_IIDS when the actual PID in EDAT is converted to an internal PID.

      NEDAT = NEDAT + 1                                    ! 8: PSHELL/PCOMP flag (to be set in subr ELEM_PROP_MATL_IIDS)
      EDAT(NEDAT) = 0

! Load a 0 into EDAT as a flag for whether this element has thicknesses defined on a continuation entry

      NEDAT = NEDAT + 1                                    ! 8: PSHELL/PCOMP flag (to be set in subr ELEM_PROP_MATL_IIDS)
      EDAT(NEDAT) = 0

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,8,0 )   ! Make sure that there are no imbedded blanks in fields 2-6
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,0,9 )   ! Issue warning if field 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Optional second card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      ! differential thickness is not supported
      IF (.FALSE. .AND. (ICONT == 1)) THEN
         ! Since there is a cont entry we assume Ti are here
         EDAT(NEDAT) = NPLATETHICK + 1

         DO J=4,6                                         ! Read 3 thicknesses
            NPLATETHICK = NPLATETHICK + 1
            IF (NPLATETHICK > LPLATETHICK) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1144) SUBR_NAME,' TOO MANY PLATE THICKNESSES. LIMIT IS NPLATETHICK = ',LPLATETHICK
               WRITE(F06,1144) SUBR_NAME,' TOO MANY PLATE THICKNESSES. LIMIT IS NPLATETHICK = ',LPLATETHICK
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            CALL R8FLD ( JCARD(J), JF(J), R8INP )
            IF (IERRFL(J) == 'N') THEN
               PLATETHICK(NPLATETHICK) = R8INP
            ENDIF
         ENDDO

         CALL BD_IMBEDDED_BLANK ( JCARD,0,0,4,5,6,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-4
         CALL CARD_FLDS_NOT_BLANK ( JCARD,2,3,0,0,0,7,8,9 )! Issue warning if fields 2, 3. 7, 8, 9 not blank
         CALL CRDERR ( CARD )                              ! CRDERR prints errors found when reading fields

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1141 FORMAT(' *ERROR  1141: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY PLATE ELEMENT MATERIAL PROPERTY ANGLES. LIMIT IS NMATANGLE =  ',I8)

 1144 FORMAT(' *ERROR  1144: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A,I8)

 1196 FORMAT(' *ERROR  1196: VALUE FOR MATERIAL ANGLE ON ',A,A,' MUST BE AN INTEGER OR REAL NUMBER BUT VALUE READ WAS ',A)

 1197 FORMAT(' *ERROR  1197: FOR ',A,A,' THE COORD SYS ID IN FIELD ',I2,' MUST BE >= 0. HOWEVER, THE VALUE INPUT WAS ',A)

! **********************************************************************************************************************************

      END SUBROUTINE BD_CTRIA


      SUBROUTINE BD_CTRIA0 ( CARD, LARGE_FLD_INP )

! Processes CTRIA Bulk Data Cards to increment LMATANGLE if the elem has a material property angle

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, LMATANGLE, LPLATEOFF, LPLATETHICK
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CTRIA0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! See if there is a material orientation angle. If so, increment LMATANGLE

      IF (JCARD(7)(1:) /= ' ') THEN
         LMATANGLE = LMATANGLE + 1
      ENDIF

! See if there is a plate offset. If so, increment LPLATEOFF

      IF (JCARD(8)(1:) /= ' ') THEN
         LPLATEOFF = LPLATEOFF + 1
      ENDIF

! Optional second card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC0  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      IF (ICONT == 1) THEN                                ! Since there is a cont entry we assume 3 Ti are here
         LPLATETHICK = LPLATETHICK + 3
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CTRIA0


      SUBROUTINE BD_CSHEAR ( CARD, NUM_GRD )

! Processes CSHEAR Bulk Data Cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IERRFL, JCARD_LEN, JF, MEDAT_CSHEAR, NCSHEAR, NELE
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, ETYPE

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_CSHEAR'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: JCARD_EDAT(10)    ! JCARD values sent to subr ELEPRO

      INTEGER(LONG), INTENT(OUT)      :: NUM_GRD           ! Number of GRID's + SPOINT's for the elem
      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
! CSHEAR element Bulk Data Card routine

!   FIELD   ITEM           ARRAY ELEMENT
!   -----   ------------   -------------
!    1      Element type   ETYPE(nele) = Q1, Q2, Q3, QA, QB
!    2      Element ID     EDAT(nedat+1)
!    3      Property ID    EDAT(nedat+2)
!    4      Grid A         EDAT(nedat+3)
!    5      Grid B         EDAT(nedat+4)
!    6      Grid C         EDAT(nedat+5)
!    7      Grid D         EDAT(nedat+6)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Set JCARD_EDAT to JCARD

      DO I=1,10
         JCARD_EDAT(I) = JCARD(I)
      ENDDO

! Check property ID field. Set to EID if blank

      IF (JCARD(3)(1:) == ' ') THEN
         JCARD_EDAT(3) = JCARD(2)
      ENDIF

! Read and check data
                                                           ! Load 6 items into EDAT
      CALL ELEPRO ( 'Y', JCARD_EDAT, 6, MEDAT_CSHEAR, 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'N', 'N' )

      ETYPE(NELE) = 'SHEAR   '
      NCSHEAR = NCSHEAR + 1
      NUM_GRD = 4

      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,6,7,0,0 )   ! Make sure that there are no imbedded blanks in fields 2-7
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,0,0,8,9 )! Issue warning if fields 8, 9 not blank
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_CSHEAR


      SUBROUTINE BD_PSHEAR ( CARD )

! Processes PSHEAR Bulk Data Cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, MPSHEAR, MRPSHEAR, NPSHEAR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  PSHEAR, RPSHEAR

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PSHEAR'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: MATERIAL_ID   = 0 ! Material ID
      INTEGER(LONG)                   :: PROPERTY_ID   = 0 ! Property ID (field 2 of this parent property card)


      REAL(DOUBLE)                    :: R8INP             ! Real value read from a field on the PSHEAR entry



! **********************************************************************************************************************************
! PSHEAR element Bulk Data Card routine

!    FIELD   ITEM           ARRAY ELEMENT
!    -----   ------------   -------------
!     2      Prop ID         PSHEAR(npshear,1)
!     3      MID1            PSHEAR(npshear,2)
!     4      Thickness      RPSHEAR(npshear,1)
!     5      NSM            RPSHEAR(npshear,2)
!     6      F1             RPSHEAR(npshear,3)
!     7      F2             RPSHEAR(npshear,4)

!  Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NPSHEAR = NPSHEAR+1

! Read and check data on parent card

      CALL I4FLD ( JCARD(2), JF(2), PROPERTY_ID )          ! Read PSHEAR ID
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NPSHEAR-1
            IF (PROPERTY_ID == PSHEAR(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),PROPERTY_ID
               WRITE(F06,1145) JCARD(1),PROPERTY_ID
               EXIT
            ENDIF
         ENDDO
         PSHEAR(NPSHEAR,1) = PROPERTY_ID
      ENDIF

! Read material ID from field 3

      CALL I4FLD ( JCARD(3), JF(3), MATERIAL_ID )
      IF (IERRFL(3) == 'N') THEN
         IF (MATERIAL_ID > 0) THEN
            PSHEAR(NPSHEAR,2) = MATERIAL_ID
         ELSE
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1193) JF(J),JCARD(1),MATERIAL_ID
            WRITE(F06,1193) JF(J),JCARD(1),MATERIAL_ID
         ENDIF
      ENDIF

! Read real values in fields 4-7

      DO J=4,7
         CALL R8FLD ( JCARD(J), JF(J), R8INP )
         IF (IERRFL(J) == 'N') THEN
            RPSHEAR(NPSHEAR,J-3) = R8INP
         ENDIF
      ENDDO

!cc NOTE: CHANGE THIS SO THAT FIELDS 6,7 CAN BE READ ONCE I ADD F1, F2 CAPABILITY IN SUBR QSHEAR
      CALL BD_IMBEDDED_BLANK   ( JCARD,2,3,4,5,0,0,0,0 )
      CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,0,6,7,8,9 )
      CALL CRDERR ( CARD )



      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1193 FORMAT(' *ERROR  1193: MATERIAL ID IN FIELD ',I3,' OF ',A,' MUST BE > 0 OR BLANK BUT IS = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PSHEAR


      SUBROUTINE BD_PSHEL ( CARD, LARGE_FLD_INP )

!  Processes PSHELL Bulk Data Cards. Reads and checks data:

!  1) Property ID and 3 material ID's (membrane, bending, transverse shear) and enter into array PSHEL
!  2) Membrane thick., bending moment of inertia ratio, shear thick. ratio, nonstr mass and enter into array RPSHEL
!  3) From 1st cont card (if present): material ID for bending/membrane coupling
!  4) From 1st cont card (if present): locations for stress recovery and offset and enter into array RPSHEL

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  IERRFL, FATAL_ERR, JCARD_LEN, JF, LPSHEL, NPSHEL, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  PSHEL, RPSHEL

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CRDERR, I4FLD, R8FLD
      USE TEXT_FIELD_UTILS, ONLY      :  CARD_FLDS_NOT_BLANK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME     = 'BD_PSHEL'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD

      INTEGER(LONG)                   :: ICONT         = 0 ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR          = 0 ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: MATERIAL_ID   = 0 ! Material ID
      INTEGER(LONG)                   :: N             = 1 ! Counter
      INTEGER(LONG)                   :: PROPERTY_ID   = 0 ! Property ID (field 2 of this parent property card)


      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
!  PSHELL Bulk Data Card routine

!    FIELD   ITEM           ARRAY ELEMENT
!    -----   ------------   -------------
!  on first card:
!     2      Prop ID         PSHEL(npshel,1)
!     3      MID1            PSHEL(npshel,2)
!     4      Thickness - TM RPSHEL(npshel,1)
!     5      MID2            PSHEL(npshel,3)
!     6      Bend. Inertia  RPSHEL(npshel,2)
!     7      MID3            PSHEL(npshel,4)
!     8      Shear - TS/TM  RPSHEL(npshel,3)
!     9      NSM            RPSHEL(npshel,4)
!  on optional second card:
!     2      Stress - Z1    RPSHEL(npshel,5)
!     3      Stress - Z2    RPSHEL(npshel,6)
!     4      MID4            PSHEL(npshel,5)
!            TS/TM indicator PSHEL(npshel,6) (set to 1 to use default value if field 8 of parent card is blank, otherwise 0 to
!                                             indicate that there is a TS/TM value in field 8)

! Blank material ID fields will be interpreted as:
!     MID1 : no membrane stiffness, no membrane/coupling stiffness
!     MID2 : no bending stiffness, no transverse shear stiffness, no membrane/bending coupling stiffness
!     MID3 : no transverse shear flexibility
!     MID4 : no membrane/coupling

!  Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Check for overflow

      NPSHEL = NPSHEL+1

! Read and check data on parent card

      CALL I4FLD ( JCARD(2), JF(2), PROPERTY_ID )          ! Read PSHEL ID
      IF (IERRFL(2) == 'N') THEN
         DO J=1,NPSHEL-1
            IF (PROPERTY_ID == PSHEL(J,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),PROPERTY_ID
               WRITE(F06,1145) JCARD(1),PROPERTY_ID
               EXIT
            ENDIF
         ENDDO
         PSHEL(NPSHEL,1) = PROPERTY_ID
      ENDIF

! Read data from fields 3,5,7 (material ID's) and 2,4,6 (TM, 12I/(TM^3), TS/TM)

      N = 1                                                ! Read 3 material ID's into array PSHEL
      DO J = 3,7,2
         N = N + 1
         CALL I4FLD ( JCARD(J), JF(J), MATERIAL_ID )
         IF (IERRFL(J) == 'N') THEN
            IF (JCARD(J)(1:) == ' ') THEN                  ! Set MID1-3 to zero if field is blank (see below for reset of MID3),
               PSHEL(NPSHEL,N) = 0
            ELSE                                           ! otherwise, read numerical value of MID1-3 and make sure it is > 0
               IF (MATERIAL_ID <= 0) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1193) JF(J),JCARD(1),MATERIAL_ID
                  WRITE(F06,1193) JF(J),JCARD(1),MATERIAL_ID
               ELSE
                  PSHEL(NPSHEL,N) = MATERIAL_ID
               ENDIF
            ENDIF
         ENDIF
         CALL R8FLD (JCARD(J+1),JF(J+1),RPSHEL(NPSHEL,N-1))! Read TM, 12I/(TM^3), TS/TM
      ENDDO

!xx   IF (JCARD(7)(1:) == ' ') THEN                        ! Reset MID3 to MID2 if MID3 field is blank
!xx      PSHEL(NPSHEL,4) = PSHEL(NPSHEL,3)
!xx   ENDIF

      IF ((IERRFL(5)=='N') .AND. (IERRFL(7)=='N')) THEN    ! Check that MID3 is not > 0 when MID2 is blank (error)
         IF ((PSHEL(NPSHEL,4) > 0) .AND. (JCARD(5)(1:) == ' ')) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1160) PSHEL(NPSHEL,1),PSHEL(NPSHEL,3),PSHEL(NPSHEL,4)
            WRITE(F06,1160) PSHEL(NPSHEL,1),PSHEL(NPSHEL,3),PSHEL(NPSHEL,4)
         ENDIF
      ENDIF

      IF (JCARD(6)(1:) == ' ') THEN                        ! Set default for bending inertia when field 6 is blank
         RPSHEL(NPSHEL,2) = ONE
      ENDIF

      IF (JCARD(8)(1:) == ' ') THEN                        ! Set PSHEL(npshel,6) if default value for TS/TM is to be used
         PSHEL(NPSHEL,6) = 1                               ! Use default value
      ELSE
         PSHEL(NPSHEL,6) = 0                               ! Do not use default value, a value was on the PSHEL card for TS/TM
      ENDIF

      CALL R8FLD ( JCARD(9), JF(9), RPSHEL(NPSHEL,4) )     ! Read NSM, nonstructural mass

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

! Read and check data on optional continuation card:

      IF (LARGE_FLD_INP == 'N') THEN
         CALL NEXTC  ( CARD, ICONT, IERR )
      ELSE
         CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
         CARD = CHILD
      ENDIF
      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

      PSHEL (NPSHEL,5) =  0                                ! Default mat'l ID for bending/membrane coupling = 0
      RPSHEL(NPSHEL,5) = -RPSHEL(NPSHEL,1)/TWO             ! Default Z1 = -TM/2
      RPSHEL(NPSHEL,6) =  RPSHEL(NPSHEL,1)/TWO             ! Default Z2 =  TM/2
      RPSHEL(NPSHEL,7) =  ZERO                             ! Default offset = 0

      IF (ICONT == 1) THEN                                 ! Read values from cont card:
         IF (JCARD(2)(1:) /= ' ') THEN                     ! Read Z1
            CALL R8FLD ( JCARD(2), JF(2), RPSHEL(NPSHEL,5) )
         ENDIF
         IF (JCARD(3)(1:) /= ' ') THEN                     ! Read Z2
            CALL R8FLD ( JCARD(3), JF(3), RPSHEL(NPSHEL,6) )
         ENDIF
         IF (JCARD(4)(1:) /= ' ') THEN                     ! Read mat'l ID for bending/membrane coupling
            CALL I4FLD ( JCARD(4), JF(4),  PSHEL(NPSHEL,5) )
            IF( PSHEL(NPSHEL,5) /= ZERO ) THEN             ! MID4 is not working yet so don't allow it.
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1194) JF(4),JCARD(1),PSHEL(NPSHEL,5)
               WRITE(F06,1194) JF(4),JCARD(1),PSHEL(NPSHEL,5)
            ENDIF
         ENDIF
         CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,0,0,0,0,0 )  ! Make sure that there are no imbedded blanks in fields 2-4
         CALL CARD_FLDS_NOT_BLANK ( JCARD,0,0,0,5,6,7,8,9 )! Issue warning if fields 5, 6, 7, 8, 9 not blank
         CALL CRDERR ( CARD )
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1160 FORMAT(' *ERROR  1160: PSHEL ID ',I8,' HAS BENDING MATERIAL ID = ',I8,' AND TRANSVERSE SHEAR MATERIAL ID = ',I8              &
                    ,/,14X,' WHEN BENDING MATERIAL ID IS 0, THE TRANSVERSE SHEAR MATERIAL ID CANNOT BE > 0')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

 1193 FORMAT(' *ERROR  1193: MATERIAL ID IN FIELD ',I3,' OF ',A,' MUST BE > 0 OR BLANK BUT IS = ',I8)

 1194 FORMAT(' *ERROR  1194: MATERIAL ID IN FIELD ',I3,' (MID4) OF PSHELL ',A,' MUST BE 0 OR BLANK BUT IS = ',I8)


! **********************************************************************************************************************************

      END SUBROUTINE BD_PSHEL


      SUBROUTINE BD_PCOMP ( CARD, LARGE_FLD_INP )

! Processes PCOMP Bulk Data Cards

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPCOMP, MPCOMP0, MRPCOMP0, MPCOMP_PLIES,  &
                                         MRPCOMP_PLIES, NPCOMP
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, TWO
      USE MODEL_STUF, ONLY            :  PCOMP, RPCOMP
      USE PARAMS, ONLY                :  EPSIL

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, R8FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PCOMP'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD               ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)          ! The 10 fields of characters making up CARD
      CHARACTER( 1*BYTE)              :: CALC_Z0     = 'N'  ! If 'Y' then calculate ZO from ply thicknesses
      CHARACTER( 1*BYTE)              :: CRD_ERR     = 'N'  ! If 'Y' then this B.D. entry has error
      CHARACTER(LEN(JCARD))           :: FT                 ! Field 6 of parent entry
      CHARACTER(LEN(JCARD))           :: LAM                ! Field 9 of parent entry
      CHARACTER(LEN(JCARD))           :: SOUT               ! Field 5 or 9 of cont entry

      INTEGER(LONG)                   :: CONT_NUM    = 0    ! Count of continuation entries for this PCOMP B.D. entry
      INTEGER(LONG)                   :: ICONT       = 0    ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR        = 0    ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: I4INP              ! Integer read from field on PCOMP entry
      INTEGER(LONG)                   :: I,J                ! DO loop indices
      INTEGER(LONG)                   :: I1,I2              ! Counters
      INTEGER(LONG)                   :: J1,J2              ! Counters
      INTEGER(LONG)                   :: PCOMP_PLIES0       ! Count of number of plies on PCOMP entry
      INTEGER(LONG)                   :: PCOMP_PLIES        ! No. of plies in 1 PCOMP entry incl sym plies not explicitly defined
      INTEGER(LONG)                   :: PROPERTY_ID = 0    ! Property ID (field 2 of this parent property card)


      REAL(DOUBLE)                    :: EPS1               ! A small number
      REAL(DOUBLE)                    :: GE                 ! Damping coeff
      REAL(DOUBLE)                    :: NSM                ! Non structural mass
      REAL(DOUBLE)                    :: R8INP              ! Real value resd from a PCOMP1 field
      REAL(DOUBLE)                    :: SB          = ZERO ! Allowable shear stress
      REAL(DOUBLE)                    :: TT          = ZERO ! Total plate thickness
      REAL(DOUBLE)                    :: TREF        = ZERO ! Ref temp
      REAL(DOUBLE)                    :: Z0          = ZERO ! Dist (+/-) from ref plane to bottom surface
      REAL(DOUBLE)                    :: ZI          = ZERO ! Dist (+/-) from ref plane to middle of ply i



! **********************************************************************************************************************************
! PCOMP Bulk Data Card:

!   FIELD   ITEM     ARRAY ELEMENT          EXPLANATION
!   -----   ----    ---------------         -------------
!    2      PID      PCOMP(npcomp,1)        Prop ID
!    3      Z0      RPCOMP(npcomp,1)        Dist from ref plane to bottom surface (default = -0.5 times layer thickness)
!    4      NSM     RPCOMP(npcomp,2)        Non structural mass per unit area
!    5      SB      RPCOMP(npcomp,3)        Allowable interlaminar shear stress. Required if FT is specified
!    6      S FT       PCOMP(npcomp,3)        Failure theory (1=HILL, 2=HOFF, 3=TSAI, 4=STRN)
!    7      TREF    RPCOMP(npcomp,4)        Ref temperature
!    8      GE      RPCOMP(npcomp,5)        Damping coeff
!    9      LAM      PCOMP(npcomp,4)        Symm lamination opt (for 1=SYM only plies on one side of elem centerline are specified)
!                                           (plies are numbered starting with 1 at the bottom layer. If an odd number of plies is
!                                           desired with SYM option the center ply thickness should be 1/2 the actual thickness
!  none              PCOMP(npcomp,2)        Type (0 for PCOMP)
!  none PCOMP_PLIES  PCOMP(npcomp,5)        Number of plies for this PCOMP (not an input - determined herein)
!  none              PCOMP(npcomp,6)        Indicator whether equiv PSHELL, MAT2 entries have been written to F06 for this PCOMP
!  none      TT     RPCOMP(npcomp,6)        Total plate thickness

! continuation cards (2 plies specified per continuation entry:
!
!   FIELD   ITEM     ARRAY ELEMENT          EXPLANATION
!   -----   -----   -----------------       -------------
!    2      MIDi     PCOMP(npcomp,6+i)      Material ID of ply i. See Remark (a)
!    3      Ti      RPCOMP(npcomp,6+i)      Thickness of ply i (real or blank but T1 -1st ply - must be specified). See Remark (b)
!    4      THETAi  RPCOMP(npcomp,7+i)      Orientation angle of longitudinal dir of ply i wrt material axis for the composite elem
!    5      SOUTi    PCOMP(npcomp,7+i)      Stress or strain output request (1=YES or 0=NO) for ply i

!    6      MIDj     repeat w/ j=i+1        Same as above for ply j=i+1
!    7      Tj
!    8      THETAj
!    9      SOUTj

!  none     Zi      RPCOMP(npcomp,8+i)      z coord (+/-) from ref plane to center of ply

! Subsequent continuation cards follow the same pattern as the 1st continuation card

! Remarks:
! -------
! (a) Blank entries for Ti   (after ply 1 will indicate thickness for ply i, Ti  ,  equals value from previous ply (i-1)
! (b) Blank entries for MIDi (after ply 1 will indicate mat'l ID  for ply i, MIDi,  equals value from previous ply (i-1)


      EPS1 = EPSIL(1)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Initialize

      LAM(1:)  = ' '
      SOUT(1:) = ' '

! Check for overflow

      NPCOMP = NPCOMP+1
!xx   IF (NPCOMP > LPCOMP) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LPCOMP
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LPCOMP
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data on parent card

      CALL I4FLD ( JCARD(2), JF(2), PROPERTY_ID )          ! Read PID  into  PCOMP(npcomp,1) and type into PCOMP(npcomp,2)
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NPCOMP-1
            IF (PROPERTY_ID == PCOMP(I,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),PROPERTY_ID
               WRITE(F06,1145) JCARD(1),PROPERTY_ID
               EXIT
            ENDIF
         ENDDO
         PCOMP(NPCOMP,1) = PROPERTY_ID
         PCOMP(NPCOMP,2) = 0                               ! 0 for PCOMP B.D entry and 1 for PCOMP1
      ENDIF

      Z0      = ZERO
      CALC_Z0 = 'N'
      IF (JCARD(3)(1:) /= ' ') THEN                        ! Read Z0   into RPCOMP(npcomp,1)
         CALL R8FLD ( JCARD(3), JF(3), Z0 )
         IF (IERRFL(3) == 'N') THEN
            RPCOMP(NPCOMP,1) = Z0
         ENDIF
      ELSE
         IF (IERRFL(3) == 'N') THEN
            CALC_Z0 = 'Y'
         ENDIF
      ENDIF

      IF (JCARD(4)(1:) /= ' ') THEN                        ! Read NSM  into RPCOMP(npcomp,2)
         CALL R8FLD ( JCARD(4), JF(4), NSM )
         IF (IERRFL(4) == 'N') THEN
            RPCOMP(NPCOMP,2) = NSM
         ENDIF
      ENDIF

      FT(1:) = ' '
      IF (JCARD(6)(1:) /= ' ') THEN                        ! Read FT   into  PCOMP(npcomp,3)
         CALL CHAR_FLD ( JCARD(6), JF(6), FT )
         IF (IERRFL(6) == 'N') THEN
            IF      (FT(1:5) == 'HILL ') THEN
               PCOMP(NPCOMP,3) = 1
            ELSE IF (FT(1:5) == 'HOFF ') THEN
               PCOMP(NPCOMP,3) = 2
            ELSE IF (FT(1:5) == 'TSAI ') THEN
               PCOMP(NPCOMP,3) = 3
            ELSE IF (FT(1:4) == 'STRE ') THEN
               PCOMP(NPCOMP,3) = 4
            ELSE IF (FT(1:5) == 'STRN ') THEN
               PCOMP(NPCOMP,3) = 5
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1146) JF(6), JCARD(6)
               WRITE(F06,1146) JF(6), JCARD(6)
            ENDIF
         ENDIF
      ENDIF

      SB = ZERO
      IF (JCARD(5)(1:) /= ' ') THEN                        ! Read SB   into RPCOMP(npcomp,3)
         CALL R8FLD ( JCARD(5), JF(5), SB )
         IF (IERRFL(6) == 'N') THEN
            RPCOMP(NPCOMP,3) = SB
         ENDIF
      ELSE                                                 ! SB field is blank so make sure FT not specified
         IF ((FT(1:) /= ' ') .AND. (DABS(SB) <= EPS1)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1112) JCARD(5)
            WRITE(F06,1112) JCARD(5)
         ENDIF
      ENDIF

      IF (JCARD(7)(1:) /= ' ') THEN                        ! Read TREF into RPCOMP(npcomp,4)
         CALL R8FLD ( JCARD(7), JF(7), TREF )
         IF (IERRFL(7) == 'N') THEN
            RPCOMP(NPCOMP,4) = TREF
         ENDIF
      ENDIF

      IF (JCARD(8)(1:) /= ' ') THEN                        ! Read GE   into RPCOMP(npcomp,5)
         CALL R8FLD ( JCARD(8), JF(8), GE )
         IF (IERRFL(8) == 'N') THEN
            RPCOMP(NPCOMP,5) = GE
         ENDIF
      ENDIF

      LAM(1:) = ' '
      CALL CHAR_FLD ( JCARD(9), JF(9), LAM )                  ! Read LAM  into  PCOMP(npcomp,4)
      IF (LAM(1:4) == 'SYM ') THEN
         PCOMP(NPCOMP,4) = 1
      ELSE
         PCOMP(NPCOMP,4) = 0
      ENDIF

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      IF ((IERRFL(2) == 'Y') .OR. (IERRFL(3) == 'Y') .OR. (IERRFL(4) == 'Y') .OR. (IERRFL(5) == 'Y') .OR. (IERRFL(6) == 'Y') .OR.  &
          (IERRFL(7) == 'Y') .OR. (IERRFL(8) == 'Y') .OR. (IERRFL(9) == 'Y')) THEN
         CRD_ERR = 'Y'
      ENDIF

! Get ply data from continuation entries. There must be at least 1 continuation entry

      I1       = 0
      I2       = 0
      CONT_NUM = 0
      PCOMP_PLIES0 = 0
      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN

            CONT_NUM = CONT_NUM + 1
                                                           ! Read 1st set of  ply data: MID, T, THRTA, SOUT
            IF((JCARD(2)(1:) /= ' ') .OR. (JCARD(3)(1:) /= ' ') .OR. (JCARD(4)(1:) /= ' ') .OR. (JCARD(5)(1:) /= ' ')) THEN

               PCOMP_PLIES0 = PCOMP_PLIES0 + 1

               I1 = I1 + 1
               IF (JCARD(2)(1:) /= ' ') THEN               ! Read MID   into  PCOMP
                  CALL I4FLD ( JCARD(2), JF(2), I4INP )
                  IF (IERRFL(2) == 'N') THEN
                     PCOMP(NPCOMP,MPCOMP0+I1) = I4INP
                  ENDIF
               ENDIF

               I2 = I2 + 1
               IF (JCARD(3)(1:) /= ' ') THEN               ! Read T     into RPCOMP
                  CALL R8FLD ( JCARD(3), JF(3), R8INP )
                  IF (IERRFL(3) == 'N') THEN
                     RPCOMP(NPCOMP,MRPCOMP0+I2) = R8INP
                  ENDIF
               ENDIF

               I2 = I2 + 1
               IF (JCARD(4)(1:) /= ' ') THEN               ! Read THETA into RPCOMP
                  CALL R8FLD ( JCARD(4), JF(4), R8INP )
                  IF (IERRFL(4) == 'N') THEN
                     RPCOMP(NPCOMP,MRPCOMP0+I2) = R8INP
                  ENDIF
               ENDIF

               I2 = I2 + 1                                 ! Make space for ZI to be calc'd below

               I1 = I1 + 1
               SOUT = 'NO      '
               IF (JCARD(5)(1:) /= ' ') THEN               ! Read SOUT  into  PCOMP
                  CALL CHAR_FLD ( JCARD(5), JF(5),  SOUT )
                  IF (IERRFL(5) == 'N') THEN
                     IF      (SOUT(1:4) == 'NO  ') THEN
                        PCOMP(NPCOMP,MPCOMP0+I1) = 0
                     ELSE IF (SOUT(1:4) == 'YES ') THEN
                        PCOMP(NPCOMP,MPCOMP0+I1) = 1
                     ELSE
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1134) SOUT
                        WRITE(F06,1134) SOUT
                     ENDIF
                  ENDIF
               ENDIF

            ENDIF
                                                           ! Read 2nd set of ply data: MID, T, THRTA, SOUT
            IF((JCARD(6)(1:) /= ' ') .OR. (JCARD(7)(1:) /= ' ') .OR. (JCARD(8)(1:) /= ' ') .OR. (JCARD(9)(1:) /= ' ')) THEN

               PCOMP_PLIES0 = PCOMP_PLIES0 + 1

               I1 = I1 + 1
               IF (JCARD(6)(1:) /= ' ') THEN               ! Read MID   into  PCOMP
                  CALL I4FLD ( JCARD(6), JF(6), I4INP )
                  IF (IERRFL(6) == 'N') THEN
                     PCOMP(NPCOMP,MPCOMP0+I1) = I4INP
                  ENDIF
               ENDIF

               I2 = I2 + 1
               IF (JCARD(7)(1:) /= ' ') THEN               ! Read T     into RPCOMP
                  CALL R8FLD ( JCARD(7), JF(7), R8INP )
                  IF (IERRFL(7) == 'N') THEN
                     RPCOMP(NPCOMP,MRPCOMP0+I2) = R8INP
                  ENDIF
               ENDIF

               I2 = I2 + 1
               IF (JCARD(8)(1:) /= ' ') THEN               ! Read THETA into RPCOMP
                  CALL R8FLD ( JCARD(8), JF(8), R8INP )
                  IF (IERRFL(8) == 'N') THEN
                     RPCOMP(NPCOMP,MRPCOMP0+I2) = R8INP
                  ENDIF
               ENDIF

               I2 = I2 + 1                                 ! Make space for ZI to be calc'd below

               I1 = I1 + 1
               SOUT(1:4) = 'NO  '
               IF (JCARD(9)(1:) /= ' ') THEN               ! Read SOUT  into  PCOMP
                  CALL CHAR_FLD ( JCARD(9), JF(9),  SOUT )
                  IF (IERRFL(9) == 'N') THEN
                     IF      (SOUT(1:4) == 'NO  ') THEN
                        PCOMP(NPCOMP,MPCOMP0+I1) = 0
                     ELSE IF (SOUT(1:4) == 'YES ') THEN
                        PCOMP(NPCOMP,MPCOMP0+I1) = 1
                     ELSE
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1134) SOUT
                        WRITE(F06,1134) SOUT
                     ENDIF
                  ENDIF
               ENDIF

            ENDIF

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9)! Make sure that there are no imbedded blanks in fields 2-9
            CALL CRDERR ( CARD )                           ! CRDERR prints errors found when reading fields

         ELSE

            EXIT

         ENDIF

      ENDDO

      IF ((IERRFL(2) == 'Y') .OR. (IERRFL(3) == 'Y') .OR. (IERRFL(4) == 'Y') .OR. (IERRFL(5) == 'Y') .OR. (IERRFL(6) == 'Y') .OR.  &
          (IERRFL(7) == 'Y') .OR. (IERRFL(8) == 'Y') .OR. (IERRFL(9) == 'Y')) THEN
         CRD_ERR = 'Y'
      ENDIF

! Make sure there is at least 1 continuation entry and, if so, that MID1 = PCOMP(NPCOMP,MPCOMP0+1)
! and T1 = RPCOMP(NPCOMP,MRPCOMP0+1) are specified

      IF (CONT_NUM >= 1) THEN

         IF (PCOMP(NPCOMP,MPCOMP0+1) == 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1114) PROPERTY_ID, PCOMP(NPCOMP,MPCOMP0+1)
            WRITE(F06,1114) PROPERTY_ID, PCOMP(NPCOMP,MPCOMP0+1)
         ENDIF

         IF (RPCOMP(NPCOMP,MRPCOMP0+1) == ZERO) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1131) PROPERTY_ID, RPCOMP(NPCOMP,MRPCOMP0+1)
            WRITE(F06,1131) PROPERTY_ID, RPCOMP(NPCOMP,MRPCOMP0+1)
         ENDIF

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1147) 'PCOMP'
         WRITE(F06,1147) 'PCOMP'

      ENDIF

! Reset MID and T that are 0 to the values from the previous ply

      IF ((CRD_ERR == 'N') .AND. (CONT_NUM >=1)) THEN

         DO I=2,PCOMP_PLIES0

            J1 = MPCOMP0 + MPCOMP_PLIES*(I-1) + 1
            IF (PCOMP(NPCOMP,J1) == 0) THEN
               PCOMP(NPCOMP,J1) = PCOMP(NPCOMP,J1-MPCOMP_PLIES)
            ENDIF

            J2 = MRPCOMP0 + MRPCOMP_PLIES*(I-1)+ 1
            IF (RPCOMP(NPCOMP,J2) == ZERO) THEN
               RPCOMP(NPCOMP,J2) = RPCOMP(NPCOMP,J2-MRPCOMP_PLIES)
            ENDIF

         ENDDO

! Total plate thickness is TT or 2*TT depending on LAM

         TT = ZERO
         DO I=1,PCOMP_PLIES0
            J2 = MRPCOMP0 + MRPCOMP_PLIES*(I-1) + 1
            TT = TT + RPCOMP(NPCOMP,J2)
         ENDDO

         IF (LAM(1:4) == 'SYM ') THEN
            TT = TWO*TT
         ENDIF
         RPCOMP(NPCOMP,6) = TT

! Rewrite Z0 into RPCOMP if it was calculated

         IF (CALC_Z0 == 'Y') THEN
            Z0 = -HALF*TT
            RPCOMP(NPCOMP,1) = Z0
         ENDIF

! Calc Zi = coord from ref plane to center of ply and add to RPCOMP for eah ply data (after THETAi)

         DO I=1,PCOMP_PLIES0
            I1 = MRPCOMP0 + MRPCOMP_PLIES*I                   ! Index where Di goes
            I2 = MRPCOMP0 + MRPCOMP_PLIES*(I-1) + 1           ! Index where ply i thickness is located
            ZI = Z0
            DO J=1,I-1
               J2 = MRPCOMP0 + MRPCOMP_PLIES*(J-1) + 1        ! Index where ply j thickness is located
               ZI = ZI + RPCOMP(NPCOMP,J2)                    ! At end of j loop, ZI = coord of bottom of ply i
            ENDDO
            ZI = ZI + HALF*RPCOMP(NPCOMP,I2)                       ! Now ZI is coordof middle of ply i
            RPCOMP(NPCOMP,I1) = ZI
         ENDDO

! If layup was symmetric (LAM = 'SYM') then duplicate ply data for symmetric plies

         IF (LAM(1:4) == 'SYM ') THEN

            DO I=1,PCOMP_PLIES0
               I1 = MPCOMP0 + MPCOMP_PLIES*(I-1) + 1
               I2 = MPCOMP0 + 2*PCOMP_PLIES0*MPCOMP_PLIES - MPCOMP_PLIES*(I-1) - (MPCOMP_PLIES - 1)
               PCOMP(NPCOMP,I2)   = PCOMP(NPCOMP,I1)
               PCOMP(NPCOMP,I2+1) = PCOMP(NPCOMP,I1+1)
            ENDDO

            DO I=1,PCOMP_PLIES0
               I1 = MRPCOMP0 + MRPCOMP_PLIES*(I-1) + 1
               I2 = MRPCOMP0 + 2*PCOMP_PLIES0*MRPCOMP_PLIES - MRPCOMP_PLIES*(I-1) - (MRPCOMP_PLIES - 1)
               RPCOMP(NPCOMP,I2)   =  RPCOMP(NPCOMP,I1)
               RPCOMP(NPCOMP,I2+1) =  RPCOMP(NPCOMP,I1+1)
               RPCOMP(NPCOMP,I2+2) = -RPCOMP(NPCOMP,I1+2)
            ENDDO

         ENDIF

! Now that we know the number of PLIES, enter this into array PCOMP

         IF (LAM(1:4) == 'SYM ') THEN
            PCOMP_PLIES = 2*PCOMP_PLIES0
         ELSE
            PCOMP_PLIES = PCOMP_PLIES0
         ENDIF
         PCOMP(NPCOMP,5) = PCOMP_PLIES

! Set flag in PCOMP col 5 that will indicate if this PCOMP's equiv PSHELL and MAT2 entries have been written to F06

         PCOMP(NPCOMP,6) = 0                               ! Will be reset to 1 when equiv PSHELL, MAT2 entries have been written

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1112 FORMAT(' *ERROR  1112: FT IN FIELD 6 IS > 0 BUT SB IN FIELD 5 = "',A,'". MUST BE > 0')

 1114 FORMAT(' *ERROR  1114: MID1 IN FIELD 2 OF 1ST CONTINUATION OF PCOMP ',I8,' MUST BE > 0 BUT IS = ',I8)

 1131 FORMAT(' *ERROR  1131: T1 IN FIELD 3 OF 1ST CONTINUATION OF PCOMP ',I8,' MUST BE > 0 BUT IS = ',1ES9.1)

 1134 FORMAT(' *ERROR  1134: SOUT IN FIELD 5 OR 9 MUST BE "YES" OF "NO" BUT IS ',A)

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1146 FORMAT(' *ERROR  1146: FAILURE THEORY IN FIELD ',I2,' MUST BE "HILL", "HOFF", "TSAI", "STRN" OR MUST BE BLANK. HOWEVER'      &
                            ,' "',A,'" WAS ENTERED')

 1147 FORMAT(' *ERROR  1147: THERE MUST BE AT LEAST 1 CONTINUATION ENTRY FOR BULK DATA ',A,' ENTRY')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PCOMP


      SUBROUTINE BD_PCOMP0 ( CARD, LARGE_FLD_INP, IPLIES )

! Processes PCOMP Bulk Data Cards to determine the number of plies there are defined for this PCOMP entry

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  f06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PCOMP0'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: LAM               ! Field 9 of parent entry. Symmetry option

      INTEGER(LONG), INTENT(OUT)      :: IPLIES            ! Count of number of plies defined by this PCOMP
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! PCOMP Bulk Data Card:

!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    2      PID             Prop ID
!    3      Z0              Dist from ref plane to bottom surface (default = -0.5 times layer thickness)
!    4      NSM             Non structural mass per unit area
!    5      SB              Allowable interlaminar shear stress. Required if FT is specified
!    6      FT              Failure theory ("HILL", "HOFF", "TSAI", "STRN")
!    7      TREF            Ref temperature
!    8      GE              Damping coeff
!    9      LAM             Symm lamination option (if "SYM" only plies on one side of elem centerline are specified)
!                           (plies are numbered starting with 1 at the bottom layer. If an odd number of plies is
!                           desired with SYM option the center ply thickness should be 1/2 the actual thickness

! continuation cards (2 plies specified per continuation entry:
!
!   FIELD   ITEM            EXPLANATION
!   -----   ------------    -------------
!    2      MID1            Material ID of ply 1
!    3      T1              Thickness of ply 1
!    4      THETA1          Orientation angle of longitudinal direction of ply 1 wrt material axis for the composite element
!    5      SOUT1           Stress or strain output request ("YES" or "NO")
!    6      MID2            Same as above for 2nd ply
!    7      T2
!    8      THETA2
!    9      SOUT2

! Subsequent continuation cards follow the same pattern as the 1st continuation card


! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      LAM = JCARD(9)

! Count number of plies on continuation cards. There can be 2 plies/cont card and a ply is assumed to exist if any of the 4
! fields for the ply have any data

      IPLIES = 0

      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC0  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN
                                                           ! See if any data is in fields 2-5. If so, up ply count
            IF((JCARD(2)(1:) /= ' ') .OR. (JCARD(3)(1:) /= ' ') .OR. (JCARD(4)(1:) /= ' ') .OR. (JCARD(5)(1:) /= ' ')) THEN
               IPLIES = IPLIES + 1
            ENDIF
                                                           ! See if any data is in fields 6-9. If so, up ply count
            IF((JCARD(6)(1:) /= ' ') .OR. (JCARD(7)(1:) /= ' ') .OR. (JCARD(8)(1:) /= ' ') .OR. (JCARD(9)(1:) /= ' ')) THEN
               IPLIES = IPLIES + 1
            ENDIF

         ELSE
            EXIT
         ENDIF
      ENDDO

      IF (LAM(1:4) == 'SYM ') THEN
         IPLIES = 2*IPLIES
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_PCOMP0


      SUBROUTINE BD_PCOMP1 ( CARD, LARGE_FLD_INP )

! Processes PCOMP1 Bulk Data Cards. Data is somewhat different than on the more general PCOMP B.D. entry but the data from this
! PCOMP1 will be put into the same PCOMP, RPCOMP arrays.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IERRFL, JCARD_LEN, JF, LPCOMP_PLIES, LPCOMP, MPCOMP0, MRPCOMP0,  &
                                         MPCOMP_PLIES, MRPCOMP_PLIES, NPCOMP
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, TWO
      USE MODEL_STUF, ONLY            :  PCOMP, RPCOMP
      USE PARAMS, ONLY                :  EPSIL

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC2
      USE BDF_FIELD_VALIDATION, ONLY  :  BD_IMBEDDED_BLANK, CHAR_FLD, CRDERR, I4FLD, R8FLD

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PCOMP1'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD               ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)          ! The 10 fields of characters making up CARD
      CHARACTER( 1*BYTE)              :: CALC_Z0     = 'N'  ! If 'Y' then calculate ZO from ply thicknesses
      CHARACTER( 1*BYTE)              :: CRD_ERR     = 'N'  ! If 'Y' then this B.D. entry has error
      CHARACTER(LEN(JCARD))           :: FT                 ! Field 6 of parent entry
      CHARACTER(LEN(JCARD))           :: LAM                ! Field 9 of parent entry

      INTEGER(LONG)                   :: CONT_NUM    = 0    ! Count of continuation entries for this PCOMP1 B.D. entry
      INTEGER(LONG)                   :: ICONT       = 0    ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR        = 0    ! Error indicator returned from subr NEXTC called herein
      INTEGER(LONG)                   :: I,J                ! DO loop indices
      INTEGER(LONG)                   :: I1,I2              ! Counters
      INTEGER(LONG)                   :: J2                 ! Counters
      INTEGER(LONG)                   :: MID                ! Material ID for all plies
      INTEGER(LONG)                   :: PCOMP_PLIES0       ! Count of number of plies on PCOMP1 entry
      INTEGER(LONG)                   :: PCOMP_PLIES        ! No. of plies in 1 PCOMP1 entry incl sym plies not explicitly defined
      INTEGER(LONG)                   :: PROPERTY_ID = 0    ! Property ID (field 2 of this parent property card)
      INTEGER(LONG)                   :: SOUT_INT    = 0    ! Entry in array PCOMP (not defined on PCOMP1)


      REAL(DOUBLE)                    :: EPS1               ! A small number
      REAL(DOUBLE)                    :: NSM                ! Non structural mass
      REAL(DOUBLE)                    :: R8INP              ! Real value resd from a PCOMP1 field
      REAL(DOUBLE)                    :: SB          = ZERO ! Allowable shear stress
      REAL(DOUBLE)                    :: THETA(LPCOMP_PLIES)! Angle of longitudinal axis of ply
      REAL(DOUBLE)                    :: THICK       = ZERO ! Thickness of each ply
      REAL(DOUBLE)                    :: TT          = ZERO ! Total plate thickness
      REAL(DOUBLE)                    :: Z0          = ZERO ! Dist (+/-) from ref plane to bottom surface
      REAL(DOUBLE)                    :: ZI          = ZERO ! Dist (+/-) from ref plane to middle of ply i



! **********************************************************************************************************************************
! PCOMP1 Bulk Data Card:

!   FIELD   ITEM     ARRAY ELEMENT          EXPLANATION
!   -----   ----    ---------------         -------------
!    2      PID      PCOMP(npcomp,1)        Prop ID
!    3      Z0      RPCOMP(npcomp,1)        Dist from ref plane to bottom surface (default = -0.5 times layer thickness)
!    4      NSM     RPCOMP(npcomp,2)        Non structural mass per unit area
!    5      SB      RPCOMP(npcomp,3)        Allowable interlaminar shear stress. Required if FT is specified
!    6      FT       PCOMP(npcomp,3)        Failure theory (1=HILL, 2=HOFF, 3=TSAI, 4=STRN)
!    7      MID      PCOMP(npcomp,6+i)      Material ID of all plies
!    8      T       RPCOMP(npcomp,7+i)      Thickness of all plies
!    9      LAM      PCOMP(npcomp,4)        Symm lamination opt (for 1=SYM only plies on one side of elem centerline are specified)
!                                           (plies are numbered starting with 1 at the bottom layer. If an odd number of plies is
!                                           desired with SYM option the center ply thickness should be 1/2 the actual thickness
!  none              PCOMP(npcomp,2)        Type (1 for PCOMP1)
!  none PCOMP_PLIES  PCOMP(npcomp,5)        Number of plies for this PCOMP1 (not an input - determined herein)
!  none              PCOMP(npcomp,6)        Indicator whether equiv PSHELL, MAT2 entries have been written to F06 for this PCOMP1
!  none             RPCOMP(npcomp,4)        TREF not defined on PCOMP1
!  none             RPCOMP(npcomp,5)        GE not defined on PCOMP1
!  none     TT      RPCOMP(npcomp,6)        Total plate thickness

! continuation cards:
!
!   FIELD   ITEM     ARRAY ELEMENT          EXPLANATION
!   -----   -----   -----------------       -------------
!    2-9    THETAi  RPCOMP(npcomp,7+i)      Orientation angle of longitudinal dir of ply i wrt material axis for the composite elem

!  none     Zi      RPCOMP(npcomp,8+i)      z coord (+/-) from ref plane to center of ply

! Subsequent continuation cards follow the same pattern as the 1st continuation card

      EPS1 = EPSIL(1)

! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

! Initialize

      LAM(1:) = ' '

! Check for overflow

      NPCOMP = NPCOMP+1
!xx   IF (NPCOMP > LPCOMP) THEN
!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,1163) SUBR_NAME,JCARD(1),LPCOMP
!xx      WRITE(F06,1163) SUBR_NAME,JCARD(1),LPCOMP
!xx      CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
!xx   ENDIF

! Read and check data on parent card

      CALL I4FLD ( JCARD(2), JF(2), PROPERTY_ID )          ! Read PID  into  PCOMP(npcomp,1) and type into PCOMP(npcomp,2)
      IF (IERRFL(2) == 'N') THEN
         DO I=1,NPCOMP-1
            IF (PROPERTY_ID == PCOMP(I,1)) THEN
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1145) JCARD(1),PROPERTY_ID
               WRITE(F06,1145) JCARD(1),PROPERTY_ID
               EXIT
            ENDIF
         ENDDO
         PCOMP(NPCOMP,1) = PROPERTY_ID
         PCOMP(NPCOMP,2) = 1                               ! 1 for PCOMP1 B.D entry and 0 for PCOMP
      ENDIF

      FT(1:) = ' '
      IF (JCARD(6)(1:) /= ' ') THEN                        ! Read FT   into  PCOMP(npcomp,3)
         CALL CHAR_FLD ( JCARD(6), JF(6), FT )
         IF (IERRFL(6) == 'N') THEN
            IF      (FT(1:5) == 'HILL ') THEN
               PCOMP(NPCOMP,3) = 1
            ELSE IF (FT(1:5) == 'HOFF ') THEN
               PCOMP(NPCOMP,3) = 2
            ELSE IF (FT(1:5) == 'TSAI ') THEN
               PCOMP(NPCOMP,3) = 3
            ELSE IF (FT(1:5) == 'STRE ') THEN
               PCOMP(NPCOMP,3) = 4
            ELSE IF (FT(1:5) == 'STRN ') THEN
               PCOMP(NPCOMP,3) = 5
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1146) JF(6), JCARD(6)
               WRITE(F06,1146) JF(6), JCARD(6)
            ENDIF
         ENDIF
      ENDIF

      IF (JCARD(7)(1:) /= ' ') THEN                        ! Read MID. Save for later entry into array PCOMP when we know # plies
         CALL I4FLD ( JCARD(7), JF(7), MID )
      ENDIF

      IF (JCARD(8)(1:) /= ' ') THEN                        ! Read THICK. Save for later entry into array RPCOMP when we know # plies
         CALL R8FLD ( JCARD(8), JF(8), THICK )
      ENDIF

      LAM(1:) = ' '
      CALL CHAR_FLD ( JCARD(9), JF(9), LAM )               ! Read LAM  into  PCOMP(npcomp,4)
      IF (LAM(1:4) == 'SYM ') THEN
         PCOMP(NPCOMP,4) = 1
      ELSE
         PCOMP(NPCOMP,4) = 0
      ENDIF

      Z0      = ZERO
      CALC_Z0 = 'N'
      IF (JCARD(3)(1:) /= ' ') THEN                        ! Read Z0   into RPCOMP(npcomp,1)
         CALL R8FLD ( JCARD(3), JF(3), Z0 )
         IF (IERRFL(3) == 'N') THEN
            RPCOMP(NPCOMP,1) = Z0
         ENDIF
      ELSE
         IF (IERRFL(3) == 'N') THEN
            CALC_Z0 = 'Y'
         ENDIF
      ENDIF

      IF (JCARD(4)(1:) /= ' ') THEN                        ! Read NSM  into RPCOMP(npcomp,2)
         CALL R8FLD ( JCARD(4), JF(4), NSM )
         IF (IERRFL(4) == 'N') THEN
            RPCOMP(NPCOMP,2) = NSM
         ENDIF
      ENDIF

      SB = ZERO
      IF (JCARD(5)(1:) /= ' ') THEN                        ! Read SB   into RPCOMP(npcomp,3)
         CALL R8FLD ( JCARD(5), JF(5), SB )
         IF (IERRFL(6) == 'N') THEN
            RPCOMP(NPCOMP,3) = SB
         ENDIF
      ELSE                                                 ! SB field is blank so make sure FT not specified
         IF ((FT(1:) /= ' ') .AND. (DABS(SB) <= EPS1)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1112) JCARD(5)
            WRITE(F06,1112) JCARD(5)
         ENDIF
      ENDIF

      RPCOMP(NPCOMP,4) = -999.D0                           ! Fictitous TREF and GE entries to make RPCOMP structure the same as
      RPCOMP(NPCOMP,5) =  ZERO                             ! for the B.D. PCOMP entry since PCOMP1 doesn't define these

      CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9 )     ! Make sure that there are no imbedded blanks in fields 2-9
      CALL CRDERR ( CARD )                                 ! CRDERR prints errors found when reading fields

      IF ((IERRFL(2) == 'Y') .OR. (IERRFL(3) == 'Y') .OR. (IERRFL(4) == 'Y') .OR. (IERRFL(5) == 'Y') .OR. (IERRFL(6) == 'Y') .OR.  &
          (IERRFL(7) == 'Y') .OR. (IERRFL(8) == 'Y') .OR. (IERRFL(9) == 'Y')) THEN
         CRD_ERR = 'Y'
      ENDIF

! Get ply data from continuation entries. There must be at least 1 continuation entry

      CONT_NUM     = 0
      PCOMP_PLIES0 = 0
      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
         IF (ICONT == 1) THEN

            CONT_NUM = CONT_NUM + 1

            DO I=2,9                                       ! Read THETA's on cont entries
               IF (JCARD(I)(1:) /= ' ') THEN
                  PCOMP_PLIES0 = PCOMP_PLIES0 + 1          ! Count of number of plies defined by this PCOMP1 entry (not incl SYM)
                  CALL R8FLD ( JCARD(I), JF(I), R8INP )
                  IF (IERRFL(I) == 'N') THEN
                     THETA(PCOMP_PLIES0) = R8INP
                  ENDIF
               ENDIF
            ENDDO

            CALL BD_IMBEDDED_BLANK ( JCARD,2,3,4,5,6,7,8,9)! Make sure that there are no imbedded blanks in fields 2-9
            CALL CRDERR ( CARD )                           ! CRDERR prints errors found when reading fields

         ELSE

            EXIT

         ENDIF

      ENDDO

      IF ((IERRFL(2) == 'Y') .OR. (IERRFL(3) == 'Y') .OR. (IERRFL(4) == 'Y') .OR. (IERRFL(5) == 'Y') .OR. (IERRFL(6) == 'Y') .OR.  &
          (IERRFL(7) == 'Y') .OR. (IERRFL(8) == 'Y') .OR. (IERRFL(9) == 'Y')) THEN
         CRD_ERR = 'Y'
      ENDIF

! Make sure there was at least 1 continuation entry

      IF (CONT_NUM < 1) THEN

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1147) 'PCOMP'
         WRITE(F06,1147) 'PCOMP'

      ENDIF

! Now that we know the number of PLIES, enter this into array PCOMP

      IF (LAM(1:4) == 'SYM ') THEN
         PCOMP_PLIES = 2*PCOMP_PLIES0
      ELSE
         PCOMP_PLIES = PCOMP_PLIES0
      ENDIF
      PCOMP(NPCOMP,5) = PCOMP_PLIES
                                                           ! Flag to indicate if this PCOMP's equiv PSHELL/MAT2 been written to F06
      PCOMP(NPCOMP,6) = 0                                  ! It will be reset to 1 when equiv PSHELL, MAT2 entries have been written

! If there were no errors on 1st entry and at least 1 cont was read, process other info for arrays PCOMP, RPCOMP

      IF ((CRD_ERR == 'N') .AND. (CONT_NUM >= 1)) THEN     ! Put MID and SOUT into PCOMP (repeated for each ply so that array PCOMP
                                                           ! will be same data structure as for BD_PCOMP
         DO I=1,PCOMP_PLIES0
            I1 = MPCOMP0 + MPCOMP_PLIES*(I-1) + 1          ! I1 = 7+2*(I-1) for params MPCOMP0 = 6, MPCOMP_PLIES = 2 (module SCONTR)
            PCOMP(NPCOMP,I1)   = MID                          ! row 7,  9, 11 etc in array PCOMP for above SCONTR params
            PCOMP(NPCOMP,I1+1) = SOUT_INT                     ! row 8, 10, 12 etc in array PCOMP for above SCONTR params
         ENDDO

         TT = PCOMP_PLIES0*THICK                           ! Tot plate thick is PCOMP_PLIES0*THICK, or twice that, depending on LAM
         IF (LAM(1:4) == 'SYM ') THEN
            TT = TWO*TT
         ENDIF
         RPCOMP(NPCOMP,6) = TT

         IF (CALC_Z0 == 'Y') THEN                          ! Rewrite Z0 into RPCOMP if it was calculated
            Z0 = -HALF*TT
            RPCOMP(NPCOMP,1) = Z0
         ENDIF

         DO I=1,PCOMP_PLIES0                               ! Put THICK, THETA into RPCOMP for each of the PCOMP_PLIES0 plies
            I1 = MRPCOMP0 + MRPCOMP_PLIES*(I-1) + 1           ! I1 = 6+3*(I-1)+1 for params MRPCOMP0 = 6, MRPCOMP_PLIES = 3 (SCONTR)
            RPCOMP(NPCOMP,I1)   = THICK
            RPCOMP(NPCOMP,I1+1) = THETA(I)
         ENDDO
                                                           ! Calc Zi = coord from ref plane to center of ply. Put into array RPCOMP
         DO I=1,PCOMP_PLIES0                                  ! I1 = 6+3*I for params MRPCOMP0 = 6, MRPCOMP_PLIES = 3 (mod SCONTR)
            I1 = MRPCOMP0 + MRPCOMP_PLIES*I                   ! Index where ZI goes
            I2 = MRPCOMP0 + MRPCOMP_PLIES*(I-1) + 1           ! Index where ply i thickness is located
            ZI = Z0
            DO J=1,I-1
               J2 = MRPCOMP0 + MRPCOMP_PLIES*(J-1) + 1        ! Index where ply j thickness is located
               ZI = ZI + RPCOMP(NPCOMP,J2)                    ! At end of j loop, ZI = coord of bottom of ply i
            ENDDO
            ZI = ZI + HALF*RPCOMP(NPCOMP,I2)                  ! Now ZI is coord of middle of ply i
            RPCOMP(NPCOMP,I1) = ZI
         ENDDO

         IF (LAM(1:4) == 'SYM ') THEN                      ! If layup was sym (LAM = 'SYM') then duplicate ply data for sym plies

            DO I=1,PCOMP_PLIES0
               I1 = MPCOMP0 + MPCOMP_PLIES*(I-1) + 1
               I2 = MPCOMP0 + 2*PCOMP_PLIES0*MPCOMP_PLIES - MPCOMP_PLIES*(I-1) - (MPCOMP_PLIES - 1)
               PCOMP(NPCOMP,I2)   = PCOMP(NPCOMP,I1)
               PCOMP(NPCOMP,I2+1) = PCOMP(NPCOMP,I1+1)
            ENDDO

            DO I=1,PCOMP_PLIES0
               I1 = MRPCOMP0 + MRPCOMP_PLIES*(I-1) + 1
               I2 = MRPCOMP0 + 2*PCOMP_PLIES0*MRPCOMP_PLIES - MRPCOMP_PLIES*(I-1) - (MRPCOMP_PLIES - 1)
               RPCOMP(NPCOMP,I2)   =  RPCOMP(NPCOMP,I1)
               RPCOMP(NPCOMP,I2+1) =  RPCOMP(NPCOMP,I1+1)
               RPCOMP(NPCOMP,I2+2) = -RPCOMP(NPCOMP,I1+2)
            ENDDO

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1112 FORMAT(' *ERROR  1112: FT IN FIELD 6 IS > 0 BUT SB IN FIELD 5 = "',A,'". MUST BE > 0')

 1114 FORMAT(' *ERROR  1114: MID1 IN FIELD 2 OF 1ST CONTINUATION OF PCOMP ',I8,' MUST BE > 0 BUT IS = ',I8)

 1131 FORMAT(' *ERROR  1131: T1 IN FIELD 3 OF 1ST CONTINUATION OF PCOMP ',I8,' MUST BE > 0 BUT IS = ',1ES9.1)

 1134 FORMAT(' *ERROR  1134: SOUT IN FIELD 5 OR 9 MUST BE "YES" OF "NO" BUT IS ',A)

 1145 FORMAT(' *ERROR  1145: DUPLICATE ',A,' ENTRY WITH ID = ',I8)

 1146 FORMAT(' *ERROR  1146: FAILURE THEORY IN FIELD ',I2,' MUST BE "HILL", "HOFF", "TSAI", "STRN" OR MUST BE BLANK. HOWEVER'      &
                            ,' "',A,'" WAS ENTERED')

 1147 FORMAT(' *ERROR  1147: THERE MUST BE AT LEAST 1 CONTINUATION ENTRY FOR BULK DATA ',A,' ENTRY')

 1163 FORMAT(' *ERROR  1163: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY ',A,' ENTRIES; LIMIT = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE BD_PCOMP1


      SUBROUTINE BD_PCOMP10 ( CARD, LARGE_FLD_INP, IPLIES )

! Processes PCOMP1 Bulk Data Cards to determine the number of plies for a B.D. PCOMP1 entry

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  f06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC0, NEXTC20

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_PCOMP10'
      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! A Bulk Data card
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      CHARACTER(LEN(CARD))            :: CHILD             ! "Child" card read in subr NEXTC, called herein
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER(LEN(JCARD))           :: LAM               ! Field 9 of parent entry. Symmetry option

      INTEGER(LONG), INTENT(OUT)      :: IPLIES            ! Count of number of plies defined by this PCOMP1
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICONT     = 0     ! Indicator of whether a cont card exists. Output from subr NEXTC
      INTEGER(LONG)                   :: IERR      = 0     ! Error indicator returned from subr NEXTC called herein




! **********************************************************************************************************************************
! Make JCARD from CARD

      CALL MKJCARD ( SUBR_NAME, CARD, JCARD )
      LAM = JCARD(9)

! Count number of plies on continuation cards. There can be 2 plies/cont card and a ply is assumed to exist if any of the 4
! fields for the ply have any data

      IPLIES = 0

      DO
         IF (LARGE_FLD_INP == 'N') THEN
            CALL NEXTC0  ( CARD, ICONT, IERR )
         ELSE
            CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
            CARD = CHILD
         ENDIF
         CALL MKJCARD ( SUBR_NAME, CARD, JCARD )

         IF (ICONT == 1) THEN
                                                           ! See if any data is in fields 2-5. If so, up ply count
            DO I=2,9
               IF (JCARD(I)(1:) /= ' ') THEN
                  IPLIES = IPLIES + 1
               ENDIF
            ENDDO
                                                           ! See if any data is in fields 6-9. If so, up ply count
         ELSE

            EXIT

         ENDIF

      ENDDO

      IF (LAM(1:4) == 'SYM ') THEN
         IPLIES = 2*IPLIES
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BD_PCOMP10

   END MODULE SHELL_COMPOSITE_CARDS
