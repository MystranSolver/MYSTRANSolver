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

   MODULE EXECUTIVE_CONTROL_IO

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: EC_DEBUG, EC_IN4FIL, EC_OUTPUT4, EC_PARTN

   CONTAINS

      SUBROUTINE EC_DEBUG ( CARD )

! Processes DEBUG in Exec Control. This is done to allow DEBUG to be read before some of the Bulk Data is input. Some of the DEBUG
! capability needs to be processed prior to the complete input data file being read

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  ECHO, IERRFL, JCARD_LEN, JF, WARN_ERR
      USE PARAMS, ONLY                :  SUPWARN
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG, NDEBUG

      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD_08
      USE BDF_FIELD_VALIDATION, ONLY  :  I4FLD

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=JCARD_LEN)        :: CHARFLD2          !
      CHARACTER(LEN=JCARD_LEN)        :: CHARFLD3          !
      CHARACTER( 8*BYTE)              :: JCARD_08(10)      ! The 10 fields of 8 characters making up CARD

      INTEGER(LONG)                   :: INDEX             ! An index into array DEBUG read on B.D. DEBUG card
      INTEGER(LONG), PARAMETER        :: LOWER    = 1      ! Lower allowable value for an integer parameter
      INTEGER(LONG)                   :: UPPER    = NDEBUG ! Upper allowable value for an integer parameter
      INTEGER(LONG)                   :: VALUE             ! Value for DEBUG(INDEX) read on B.D. DEBUG card

! **********************************************************************************************************************************
      CHARFLD2(1:) = ' '
      CHARFLD3(1:) = ' '

! Make JCARD_08 from CARD

      CALL MKJCARD_08 ( CARD, JCARD_08 )
      CHARFLD2(1:) = JCARD_08(2)(1:)
      CHARFLD3(1:) = JCARD_08(3)(1:)

! Read DEBUG index and DEBUG(INDEX) value

      CALL I4FLD ( CHARFLD2, JF(2), INDEX )
      IF (IERRFL(2) == 'N') THEN
         IF ((INDEX >= LOWER) .AND. (INDEX <= UPPER)) THEN
            CALL I4FLD ( CHARFLD3, JF(3), VALUE )
            IF (IERRFL(3) == 'N') THEN
               DEBUG(INDEX) = VALUE
            ENDIF
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) CARD
            WRITE(ERR,1120) LOWER,UPPER,INDEX
            IF (SUPWARN == 'N') THEN
               IF (ECHO == 'NONE  ') THEN
                  WRITE(F06,101) CARD
               ENDIF
               WRITE(F06,1120) LOWER,UPPER,INDEX
            ENDIF
         ENDIF
      ENDIF

! **********************************************************************************************************************************
  101 FORMAT(A)

 1120 FORMAT(' *WARNING    : DEBUG INDEX MUST BE >= ',I4,' AND <= ',I4,' BUT INPUT VALUE IS: ',I8,'. ENTRY IGNORED')

! **********************************************************************************************************************************

      END SUBROUTINE EC_DEBUG


      SUBROUTINE EC_IN4FIL ( CARD )

! Processes Executive Control IN4 entries that define IN4 files to be read

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, FILE_NAM_MAXLEN, IN4FIL, IN4FIL_NUM, NUM_IN4_FILES
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_TOKEN_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EC_IN4FIL'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(LEN=LEN(CARD))        :: CARD1             ! Part of CARD
      CHARACTER(LEN=LEN(CARD))        :: CARD2             ! Part of CARD
      CHARACTER( 8*BYTE)              :: TOKEN(3)          ! Char string output from subr STOKEN, called herein
      CHARACTER( 8*BYTE)              :: TOKTYP(3)         ! Type of the char TOKEN's output from subr STOKEN, called herein

      INTEGER(LONG)                   :: DOLLAR_COL        ! Col, on CARD, where "$" sign is located
      INTEGER(LONG)                   :: ECOL      = 0     ! Col, on CARD, where "=" sign is located
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IEND              !
      INTEGER(LONG)                   :: ISTART            !
      INTEGER(LONG)                   :: JEND      = 0     !
      INTEGER(LONG)                   :: K         = 0     ! Counter
      INTEGER(LONG)                   :: SETERR    = 0     ! Error indicator as set ID is read


      INTRINSIC INDEX



! **********************************************************************************************************************************
! Process SET cards

      CARD1 = CARD

! Find equal sign

      ECOL = INDEX(CARD1(1:),'=')
      IF (ECOL == 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1267)
         WRITE(F06,1267)
         RETURN
      ENDIF

! Now find IN4FIL_NUM and check for proper type ('INTEGER ') and value within limits

      NUM_IN4_FILES = NUM_IN4_FILES + 1
      TOKEN(1)  = '        '
      SETERR    = 0
      K = 0
      DO I=5,ECOL-1
         IF (CARD1(I:I) == ' ') CYCLE
         K = K + 1
         IF (K > MAX_TOKEN_LEN) THEN
            SETERR = 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1272) MAX_TOKEN_LEN
            WRITE(F06,1272) MAX_TOKEN_LEN
            EXIT
         ENDIF
         TOKEN(1)(K:K) = CARD1(I:I)
      ENDDO
      IF (SETERR == 0) THEN
         CALL TOKCHK ( TOKEN(1), TOKTYP(1) )
         IF (TOKTYP(1) == 'INTEGER ') THEN
            READ(TOKEN(1),'(I8)') IN4FIL_NUM(NUM_IN4_FILES)
         ELSE
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1273)
            WRITE(F06,1273)
            RETURN
         ENDIF
      ENDIF

! Remainder of card has IN4 file name (and also possibly leading/trailing blanks and trailing comment)

      CARD2 = CARD1(ECOL+1:)
      DOLLAR_COL = INDEX(CARD2(1:),'$')

      DO I=1,FILE_NAM_MAXLEN
         IF (CARD2(I:I) /= ' ') THEN
            ISTART = I
            EXIT
         ENDIF
      ENDDO

      JEND = FILE_NAM_MAXLEN
      IF (DOLLAR_COL > 0) THEN
         JEND = DOLLAR_COL-1
      ENDIF

      DO I=JEND,ISTART,-1
         IF (CARD2(I:I) == ' ') THEN
            CYCLE
         ELSE
            IEND = I
            EXIT
         ENDIF
      ENDDO

      IN4FIL(NUM_IN4_FILES)(1:) = CARD2(ISTART:IEND)



      RETURN

! **********************************************************************************************************************************
 1267 FORMAT(' *ERROR  1267: FORMAT ERROR ON CASE CONTROL SET ENTRY. COULD NOT FIND EQUAL (=) SIGN')

 1272 FORMAT(' *ERROR  1272: SET ID ON CASE CONTROL SET ENTRY CANNOT BE HAVE MORE THAN ',I4,' DIGITS')

 1273 FORMAT(' *ERROR  1273: THE SET ID ON ABOVE CASE CONTROL SET ENTRY IS NOT AN INTEGER')

97865 format(' I, DOLLAR_COL, JEND, ISTART, IEND, IN4FIL(I) = ',5i8,'"',a,'"')

! **********************************************************************************************************************************

      END SUBROUTINE EC_IN4FIL


      SUBROUTINE EC_OUTPUT4 ( CARD1, IERR, ANY_OU4_NAME_BAD )

! EC_OUTPUT4 reads in the Exec Control entry OUTPUT4. The form of the entry is:

!     OUTPUT4 M1,M2,M3,M4,M5//ITAPE/IUNIT       (ITAPE not used)

! From 1 to 5 matrices can be requested. All 4 commas must be present. An example requesting MLL and KLL is:

!     OUTPUT4 MLL,KLL,,,//-1/21
!     OUTPUT4 MLL,KLL,,,//

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, EC_ENTRY_LEN
      USE IOUNT1, ONLY                :  ERR, F06, MOU4, OU4, OU4_ELM_OTM, OU4_GRD_OTM, SC1
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE OUTPUT4_MATRICES, ONLY      :  NUM_OU4_VALID_NAMES, TAPE_ACTION_MAX_VAL, TAPE_ACTION_MIN_VAL, NUM_OU4_REQUESTS,          &
                                         OU4_FILE_UNITS, OU4_TAPE_ACTION, ACT_OU4_MYSTRAN_NAMES, ACT_OU4_OUTPUT_NAMES,             &
                                         ALLOW_OU4_MYSTRAN_NAMES, ALLOW_OU4_OUTPUT_NAMES

      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE TEXT_FIELD_UTILS, ONLY      :  TO_UPPER

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EC_OUTPUT4'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD1             ! Card read in LOADE and shifted to begin in col 1
      CHARACTER(LEN=*), INTENT(OUT)   :: ANY_OU4_NAME_BAD  ! 'Y'/'N' if requested OUTPUT4 matrix name is valid
      CHARACTER(16*BYTE)              :: MYSTRAN_NAME(5)   !
      CHARACTER(16*BYTE)              :: OUTPUT_NAME(5)    !
      CHARACTER(LEN=LEN(CARD1))       :: CARD2             ! CARD1 truncated at $ (trailing comment) if there is one
      CHARACTER(LEN=EC_ENTRY_LEN)     :: DATA_80(7)        ! Temp slot for holding data until lead/trail blanks stripped
      CHARACTER(16*BYTE)              :: DATA_16(7)        ! Matrix name read from OUTPUT4 entry
      CHARACTER( 1*BYTE)              :: DUPLICATE         ! 'Y'/'N' if requested matrix is a duplicate of any previous one
      CHARACTER( 1*BYTE)              :: ITAPE_OK          ! 'Y'/'N' if ITAPE value is OK
      CHARACTER( 1*BYTE)              :: IUNIT_OK          ! 'Y'/'N' if IUNIT value is OK
      CHARACTER( 1*BYTE)              :: VALID_OU4_NAME    ! 'Y'/'N' if requested OUTPUT4 matrix name is valid

      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator. If CHAR not found, IERR set to 1
      INTEGER(LONG)                   :: DATA_BEG          ! Column where data begins (after OUTPUT4)
      INTEGER(LONG)                   :: DELTA_COL         ! Column where data begins (after OUTPUT4)
      INTEGER(LONG)                   :: DOLLAR_COL        ! Column where $ is located (for comments)
      INTEGER(LONG)                   :: COMMA_COL(5)      ! Column where comma is found in CARD2
      INTEGER(LONG)                   :: I,J               ! DO loop index
      INTEGER(LONG)                   :: ITAPE             ! Tape action (rewind before write, etc)
      INTEGER(LONG)                   :: IUNIT             ! OUTPUT4 unit number to write the matrices to
      INTEGER(LONG)                   :: JBEG              ! Beg col in data
      INTEGER(LONG)                   :: ROW_NUM           !
      INTEGER(LONG)                   :: SLASH1_COL        ! Col in matrix name where character  "/"  is  found
      INTEGER(LONG)                   :: SLASH2_COL        ! Col in matrix name where characters "//" are found




! **********************************************************************************************************************************
! Initialize

      IERR             =  0
      VALID_OU4_NAME   = 'Y'
      ANY_OU4_NAME_BAD = 'N'

      DO I=1,5
         DATA_16(I)(1:) = ' '
         DATA_80(I)(1:) = ' '
      ENDDO

! Find $ in CARD1 (if it is there) and truncate to get CARD2 (i.e. get rid of all data from $ sign on)

      CARD2 = CARD1
      DOLLAR_COL = INDEX ( CARD1, "$" )
      IF (DOLLAR_COL > 0) THEN
         CARD2(1:) = CARD1(1:DOLLAR_COL-1)
      ENDIF

! Find where data begins after OUTPUT4. NOTE: input CARD1 was shifted in subr LOADE so that the "OUTPUT4" began in col 1. Thus the
! data will begin after col 8 (after "OUTPUT4")

      DO I=8,EC_ENTRY_LEN
         IF (CARD2(I:I) == ' ') THEN
            CYCLE
         ELSE
            DATA_BEG = I
            EXIT
         ENDIF
      ENDDO

! Find commas (there should be no 5th comma - this is checked below)

      COMMA_COL(1) =                INDEX ( CARD2             , "," )
      COMMA_COL(2) = COMMA_COL(1) + INDEX ( CARD2(COMMA_COL(1)+1:), "," )
      COMMA_COL(3) = COMMA_COL(2) + INDEX ( CARD2(COMMA_COL(2)+1:), "," )
      COMMA_COL(4) = COMMA_COL(3) + INDEX ( CARD2(COMMA_COL(3)+1:), "," )
      COMMA_COL(5) = COMMA_COL(4) + INDEX ( CARD2(COMMA_COL(4)+1:), "," )

! Make sure there are 4 commas

      IF ((COMMA_COL(1) == 0) .OR. (COMMA_COL(2) == COMMA_COL(1)) .OR. (COMMA_COL(3) == COMMA_COL(2)) .OR.                         &
                                   (COMMA_COL(4) == COMMA_COL(3))) THEN
         IERR = IERR + 1
         WRITE(ERR,1031) '4'
         WRITE(F06,1031) '4'
      ENDIF

! Make sure there are not more than 4 commas

      IF (COMMA_COL(5) /= COMMA_COL(4)) THEN
         IERR = IERR + 1
         WRITE(ERR,1032) '4'
         WRITE(F06,1032) '4'
      ENDIF

! Make sure there is // somewhat after the 4th comma followed later by /

      SLASH1_COL = 0
      SLASH2_COL = 0
      DELTA_COL = INDEX ( CARD2(COMMA_COL(4)+1:), "//" )
      IF (DELTA_COL /= 0) THEN                             ! There is // so now see if there is a later /
         SLASH2_COL = COMMA_COL(4) + DELTA_COL
         SLASH1_COL = 0
         DELTA_COL = INDEX ( CARD2(SLASH2_COL+2:), "/" )
         IF (DELTA_COL /= 0) THEN
            SLASH1_COL = (SLASH2_COL + 2) + DELTA_COL - 1
         ELSE
            IERR = IERR + 1
            WRITE(ERR,1035)
            WRITE(F06,1035)
         ENDIF
      ELSE                                                 ! No // so give error
         IERR = IERR + 1
         WRITE(ERR,1034)
         WRITE(F06,1034)
      ENDIF

! Process data from CARD2 if there were no errors

nerr: IF (IERR == 0) THEN
                                                           ! DATA_80(1-5) is the data from CARD2 bet commas from data beg to //:
         DATA_80(1)(1:) = CARD2 ( DATA_BEG       : COMMA_COL(1)-1 )
         DATA_80(2)(1:) = CARD2 ( COMMA_COL(1)+1 : COMMA_COL(2)-1 )
         DATA_80(3)(1:) = CARD2 ( COMMA_COL(2)+1 : COMMA_COL(3)-1 )
         DATA_80(4)(1:) = CARD2 ( COMMA_COL(3)+1 : COMMA_COL(4)-1 )
         DATA_80(5)(1:) = CARD2 ( COMMA_COL(4)+1 : SLASH2_COL-1)

         DATA_80(6) = CARD2(SLASH2_COL+2:SLASH1_COL-1)     ! ITAPE entry
         DATA_80(7) = CARD2(SLASH1_COL+1:)                 ! IUNIT entry

         DO I=1,7                                          ! Get 8 char matrix names by stripping leading/trailing blanks
            IF (DATA_80(I)(1:) /= ' ') THEN
               DO J=1,EC_ENTRY_LEN
                  IF (DATA_80(I)(J:J) == ' ') THEN
                     CYCLE
                  ELSE
                     JBEG = J
                     EXIT
                  ENDIF
               ENDDO
               DO J=EC_ENTRY_LEN,JBEG,-1
                  IF (DATA_80(I)(J:J) == ' ') THEN
                     CYCLE
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
               DATA_16(I)(1:) = DATA_80(I)(JBEG:JBEG+15)
            ELSE
               DATA_16(I) = '                '
            ENDIF
         ENDDO

         READ(DATA_16(6),'(I2)') ITAPE                     ! Check ITAPE against valid values for OUTPUT4 matrices
         ITAPE_OK = 'N'
         IF ((ITAPE <= TAPE_ACTION_MAX_VAL) .AND. (ITAPE >= TAPE_ACTION_MIN_VAL)) THEN
            ITAPE_OK = 'Y'
         ENDIF
         IF (ITAPE_OK == 'N') THEN
            IERR = IERR + 1
            WRITE(ERR,1036 ) ITAPE, TAPE_ACTION_MAX_VAL,TAPE_ACTION_MIN_VAL
            WRITE(F06,1036 ) ITAPE, TAPE_ACTION_MAX_VAL,TAPE_ACTION_MIN_VAL
         ENDIF

         READ(DATA_16(7),'(I2)') IUNIT                     ! Check IUNIT against valid unit no's for OUTPUT4 matrices
         IUNIT_OK = 'N'
         DO I=1,MOU4
            IF (IUNIT == OU4(I)) THEN
               IUNIT_OK = 'Y'
               EXIT
            ELSE
               CYCLE
            ENDIF
         ENDDO
         IF ((IUNIT == OU4_ELM_OTM) .OR. (IUNIT == OU4_GRD_OTM)) IUNIT_OK = 'N'
         IF (IUNIT_OK == 'N') THEN
            IERR = IERR + 1
            WRITE(ERR,1037 ) IUNIT
            WRITE(ERR,10371) (OU4(I),I=1,MOU4-1)
            WRITE(ERR,10372)  OU4(MOU4)
            WRITE(F06,1037 ) IUNIT
            WRITE(F06,10371) (OU4(I),I=1,MOU4-1)
            WRITE(F06,10372)  OU4(MOU4)
         ENDIF

         DO I=1,5                                          ! Check to see if this is an alias name
            IF (DATA_16(I)(1:) /= ' ') THEN
               CALL OU4_NAME_ALIAS (DATA_16(I), MYSTRAN_NAME(I), OUTPUT_NAME(I) )
            ENDIF
         ENDDO

         IF (IERR == 0) THEN
            DO I=1,5                                       ! Update list of requested OUTPUT4 matrices. Don't enter duplicates
               IF (DATA_16(I)(1:16) /= ' ') THEN
                  IF (NUM_OU4_REQUESTS+1 <= NUM_OU4_VALID_NAMES) THEN
                     CALL CHECK_MATRIX_NAME ( MYSTRAN_NAME(I), ROW_NUM )
                     IF ((VALID_OU4_NAME == 'Y') .AND. (DUPLICATE == 'N')) THEN
                        NUM_OU4_REQUESTS = NUM_OU4_REQUESTS + 1
                        ACT_OU4_OUTPUT_NAMES(NUM_OU4_REQUESTS)  = OUTPUT_NAME(I)
                        ACT_OU4_MYSTRAN_NAMES(NUM_OU4_REQUESTS) = ALLOW_OU4_MYSTRAN_NAMES(ROW_NUM)
                        OU4_FILE_UNITS (NUM_OU4_REQUESTS) = IUNIT
                        OU4_TAPE_ACTION(NUM_OU4_REQUESTS) = ITAPE
                     ENDIF
                  ELSE
                     IERR = IERR + 1
                     WRITE(ERR,1033) SUBR_NAME, NUM_OU4_VALID_NAMES
                     WRITE(F06,1033) SUBR_NAME, NUM_OU4_VALID_NAMES
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

      ENDIF nerr


      IF ((DEBUG(197) == 1) .OR. (DEBUG(197) == 3)) THEN

         WRITE(F06,'(A)') '********************************************************************************************************'
         WRITE(F06,'(A)') 'Debug output from subr EC_OUTPUT4'
         WRITE(F06,'(A)') '---------------------------------'
         WRITE(F06,*)
         WRITE(F06,*) CARD1
         WRITE(F06,*) 'Names for the ',num_ou4_valid_names,' files that can be processed for OUTPUT4 are:'
         DO I=1, NUM_OU4_VALID_NAMES
            WRITE(F06,88987) I, ALLOW_OU4_OUTPUT_NAMES(I), ALLOW_OU4_MYSTRAN_NAMES(I)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,*) 'So far there are ',NUM_OU4_REQUESTS,'files requested for OUTPUT4:'
         WRITE(F06,88988)
         DO I=1, NUM_OU4_REQUESTS
            WRITE(F06,88989) I, OU4_FILE_UNITS(I), ACT_OU4_OUTPUT_NAMES(I), ACT_OU4_MYSTRAN_NAMES(I)
         ENDDO
         WRITE(F06,*)
                                                           ! F06 messages
         WRITE(F06,*)
         WRITE(F06,99880) CARD1
         WRITE(F06,*)
         WRITE(F06,99881) CARD2
         WRITE(F06,*)

         WRITE(F06,99882) DATA_BEG
         WRITE(F06,*)

         WRITE(F06,99883) (COMMA_COL(I),I=1,4)
         WRITE(F06,*)

         IF (COMMA_COL(5) /= COMMA_COL(4)) THEN
            WRITE(F06,99884) COMMA_COL(5)
            WRITE(F06,*)
         ENDIF

         WRITE(F06,99885) SLASH2_COL
         WRITE(F06,99886) SLASH1_COL
         WRITE(F06,*)

         IF (IERR == 0) THEN

            WRITE(F06,99890) (I, DATA_80(I),I=1,7)
            WRITE(F06,*)

            DO I=1,7
               WRITE(F06,99891) I, DATA_16(I)
            ENDDO
            WRITE(F06,*)

            DO I=1,NUM_OU4_REQUESTS
               WRITE(F06,99892) I, ACT_OU4_MYSTRAN_NAMES(I), OU4_FILE_UNITS(I), OU4_TAPE_ACTION(I)
            ENDDO
            WRITE(F06,*)

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1031 FORMAT(' *ERROR  1031: EXEC CONTROL ENTRY OUTPUT4 MUST HAVE ',A,' COMMAS.')

 1032 FORMAT(' *ERROR  1032: EXEC CONTROL ENTRY OUTPUT4 CANNOT HAVE MORE THAN ',A,' COMMAS.')

 1033 FORMAT(' *ERROR  1033: PROGRAMMING ERROR IN SUBR ',A                                                                         &
                    ,/,14X,' THE TOTAL NUMBER OF MATRICES REQUESTED FOR OUTPUT4 EXCEEDS THE NUMBER OF AVAILABLE MATRICES = ',I4)

 1034 FORMAT(' *ERROR  1034: THE OUTPUT4 ENTRY FORMAT REQUIRES A DOUBLE SLASH (//) AFTER THE MATRICES BUT NONE WAS FOUND')

 1035 FORMAT(' *ERROR  1035: THE OUTPUT4 ENTRY FORMAT REQUIRES A SINGLE SLASH (/) BEFORE THE FILE UNIT NUMBER BUT NONE WAS FOUND')

 1036 FORMAT(' *ERROR  1036: REQUESTED OUTPUT4 ITAPE ',I2,' IS NOT VALID. IT MUST BE <= ',I3,' AND >= ',I3)

 1037 FORMAT(' *ERROR  1037: REQUESTED OUTPUT4 UNIT ',I2,' IS NOT VALID. THE VALID UNITS ARE LISTED BELOW:')

10371 FORMAT('               UNIT ',I2)

10372 FORMAT('               UNIT ',I2,' reserved for OTM outputs requested via CASE CONTROL entries')

88987 FORMAT(' I, ALLOW_OU4_OUTPUT_NAMES(I), ALLOW_OU4_MYSTRAN_NAMES(I) =',i3,', "',a,'"',', "',a,'"')

88988 FORMAT('         File unit    Output name        MYSTRAN name')

88989 FORMAT(3X,'(',I3,')',3X,I3,10X,A,3X,A)

99880 FORMAT('CARD1:',/,A)

99881 FORMAT('CARD2:',/,A)

99882 FORMAT('Data begins in OUTPUT4 entry in col ',I3)

99883 FORMAT('COMMA(1:4)                      = ',3(I4,','),I4,' (cols where commas exist)')

99884 FORMAT('There are more than 4 commas. A 5th one was found in col ',I4)

99885 FORMAT('Col where // exists             = ',I4,' (0 indicates no //)')

99886 FORMAT('Col where /  exists             = ',I4,' (0 indicates no / following //)')

99890 FORMAT('I, DATA_80(I)                   = ',I4,',',3X,'"',A,'"')

99891 FORMAT('I, DATA_16(I)                   = ',I4,',',3X,'"',A,'"')

99892 FORMAT('I, OU4_MYSTRAN_NAMES(I), Unt No = ',I4,',',3X,'"',A,'"',',',2(3X,I2))

99893 FORMAT('   ITAPE                        = ',i4)

99894 FORMAT('   IUNIT                        = ',I4)

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE OU4_NAME_ALIAS ( INPUT_NAME, MYSTRAN_NAME, OUTPUT_NAME )

      USE PENTIUM_II_KIND, ONLY       :  BYTE

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(INOUT) :: INPUT_NAME
      CHARACTER(LEN=*), INTENT(OUT)   :: MYSTRAN_NAME
      CHARACTER(LEN=*), INTENT(OUT)   :: OUTPUT_NAME

      INTEGER(LONG)                   :: JJ                ! DO loop index

! **********************************************************************************************************************************
      OUTPUT_NAME  = INPUT_NAME
      MYSTRAN_NAME = INPUT_NAME
      DO JJ=1,NUM_OU4_VALID_NAMES
         IF (INPUT_NAME == ALLOW_OU4_OUTPUT_NAMES(JJ)) THEN
            OUTPUT_NAME  = ALLOW_OU4_OUTPUT_NAMES(JJ)
            MYSTRAN_NAME = ALLOW_OU4_MYSTRAN_NAMES(JJ)
         ENDIF
      ENDDO

! **********************************************************************************************************************************

      END SUBROUTINE OU4_NAME_ALIAS

! ##################################################################################################################################

      SUBROUTINE CHECK_MATRIX_NAME ( MYSTRAN_NAME, INDEX )

      USE PENTIUM_II_KIND

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: MYSTRAN_NAME
      CHARACTER( 1*BYTE)              :: FOUND             ! 'Y'/'N' if something was found

      INTEGER(LONG), INTENT(OUT)      :: INDEX             ! Row in array ALLOW_OU4_MYSTRAN_NAMES where name was found
      INTEGER(LONG)                   :: JJ                ! DO loop index

      CHARACTER(16*BYTE)              :: ALLOW_OU4_MYSTRAN_NAME_UPPER
      CHARACTER(16*BYTE)              :: ALLOW_OU4_OUTPUT_NAME_UPPER

! **********************************************************************************************************************************
! Check requested OUTPUT4 names to make sure they are in the list of valid names

      VALID_OU4_NAME = 'Y'
      DO JJ=1,NUM_OU4_VALID_NAMES
         FOUND = 'N'
         ! MYSTRAN_NAME was upper-cased when reading the card.
         ! Here, upper-case the allowed names so that comparison is case-insensitive
         ! and thus allows KRRcb, MRRcb, KLR(t).
         ALLOW_OU4_MYSTRAN_NAME_UPPER = TO_UPPER(ALLOW_OU4_MYSTRAN_NAMES(JJ))
         ALLOW_OU4_OUTPUT_NAME_UPPER = TO_UPPER(ALLOW_OU4_OUTPUT_NAMES(JJ))
         IF ((MYSTRAN_NAME == ALLOW_OU4_MYSTRAN_NAME_UPPER) .OR. (MYSTRAN_NAME == ALLOW_OU4_OUTPUT_NAME_UPPER)) THEN
            FOUND = 'Y'
            INDEX = JJ
            EXIT
         ENDIF
      ENDDO
      IF (FOUND == 'N') THEN
         IERR = IERR + 1
         VALID_OU4_NAME   = 'N'
         ANY_OU4_NAME_BAD = 'Y'
         WRITE(ERR,1026) MYSTRAN_NAME
         WRITE(F06,1026) MYSTRAN_NAME
      ENDIF

! Check if name is a duplicate of one alreay requested

      DUPLICATE = 'N'
      DO JJ=1,NUM_OU4_REQUESTS
         IF ((MYSTRAN_NAME == ACT_OU4_MYSTRAN_NAMES(JJ))  .AND. (IUNIT == OU4_FILE_UNITS(JJ)))THEN
            DUPLICATE = 'Y'
            EXIT
         ENDIF
      ENDDO

! **********************************************************************************************************************************
 1026 FORMAT(' *ERROR  1026: OUTPUT4 REQUESTED MATRIX NAMED ',A,' IS NOT IN THE LIST OF VALID OUTPUT4 NAMES.')

! **********************************************************************************************************************************

      END SUBROUTINE CHECK_MATRIX_NAME
     
      
      END SUBROUTINE EC_OUTPUT4


      SUBROUTINE EC_PARTN ( CARD1, IERR )

! EC_PARTN reads in the Exec Control entry PARTN (partition a matrix). The NASTRAN form of the entry is:

!     PARTN A,CP,RP/A11,A21,A12,A22/SYM/TYPE/F11,F21,F12,F22

! The MYSTRAN form will be the same except only the data up to the 1st slash (/), if it exists, is processed:

!     PARTN A,CP,RP

! Matrix (A) will be partitioned using the column (CP) and row (RP) partitioning vectors.  The partition will be what is output as
! the OUTPUT4 matrix. PARTN must follow an OUTPUT4 requesting output of the matrix A

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  EC_ENTRY_LEN
      USE IOUNT1, ONLY                :  ERR, F06, MOU4, OU4, OU4_ELM_OTM, OU4_GRD_OTM
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE OUTPUT4_MATRICES, ONLY      :  NUM_OU4_REQUESTS, NUM_PARTN_REQUESTS, OU4_PART_VEC_NAMES, OU4_PART_MAT_NAMES,             &
                                         ACT_OU4_MYSTRAN_NAMES, ACT_OU4_OUTPUT_NAMES

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE TEXT_FIELD_UTILS, ONLY      :  TO_UPPER

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: CARD1             ! Card read in LOADE and shifted to begin in col 1
      CHARACTER(LEN=LEN(CARD1))       :: CARD2             ! CARD1 truncated at $ (trailing comment) if there is one
      CHARACTER(LEN=EC_ENTRY_LEN)     :: DATA_80(3)        ! Temp slot for holding data until lead/trail blanks stripped
      CHARACTER(16*BYTE)              :: DATA_16(3)        ! Matrix name read from OUTPUT4 entry
      CHARACTER( 1*BYTE)              :: FOUND             ! 'Y' if we found something we were looking for

      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator. If CHAR not found, IERR set to 1
      INTEGER(LONG)                   :: DATA_BEG          ! Column where data begins (after OUTPUT4)
      INTEGER(LONG)                   :: DATA_END          ! Column where data ends (after ist slash or at $)
      INTEGER(LONG)                   :: COMMA_COL(3)      ! Column where comma is found in CARD2
      INTEGER(LONG)                   :: I,J               ! DO loop index
      INTEGER(LONG)                   :: JBEG              ! Beg col in data




! **********************************************************************************************************************************
      IF ((DEBUG(197) == 2) .OR. (DEBUG(197) == 3)) THEN
         WRITE(F06,'(A)') '******************************************************************************************************'
      ENDIF

! Initialize

      IERR =  0

      DO I=1,3
         DATA_16(I)(1:) = ' '
         DATA_80(I)(1:) = ' '
      ENDDO

      CARD2(1:) = ' '

! Find where data begins after PARTN. NOTE: input CARD1 was shifted in subr LOADE so that the "PARTN" began in col 1. Thus the
! data will begin after col 6 (after "PARTN")

      DATA_BEG = 1
      DO I=6,EC_ENTRY_LEN
         IF (CARD1(I:I) == ' ') THEN
            CYCLE
         ELSE
            DATA_BEG = I
            EXIT
         ENDIF
      ENDDO

! Find where data we want ends (no data will be kept past a "/" or a "$" character)

      DATA_END = EC_ENTRY_LEN
      DO I=6,EC_ENTRY_LEN
         IF ((CARD1(I:I) == "/") .OR. (CARD2(I:I) == "$")) THEN
            DATA_END = I - 1
            EXIT
         ENDIF
      ENDDO

! The data we want is in CARD1 from DATA_BEG thru DATA_END

      CARD2(1:) = CARD1(DATA_BEG:DATA_END)

! Find commas (there should be no 3rd comma - this is checked below)

      COMMA_COL(1) =                INDEX ( CARD2             , "," )
      COMMA_COL(2) = COMMA_COL(1) + INDEX ( CARD2(COMMA_COL(1)+1:), "," )
      COMMA_COL(3) = COMMA_COL(2) + INDEX ( CARD2(COMMA_COL(2)+1:), "," )

! Make sure there are 2 commas

      IF ((COMMA_COL(1) == 0) .OR. (COMMA_COL(2) == COMMA_COL(1))) THEN
         IERR = IERR + 1
         WRITE(ERR,1031) '2'
         WRITE(F06,1031) '2'
      ENDIF

! Make sure there are not more than 2 commas

      IF (COMMA_COL(3) /= COMMA_COL(2)) THEN
         IERR = IERR + 1
         WRITE(ERR,1032) '2'
         WRITE(F06,1032) '2'
      ENDIF

! Process data from CARD2 if there were no errors

nerr: IF (IERR == 0) THEN

         NUM_PARTN_REQUESTS = NUM_PARTN_REQUESTS + 1
                                                           ! DATA_80(1-3) is the data from CARD2 bet commas from data beg to //:
         DATA_80(1)(1:) = CARD2 ( 1              : COMMA_COL(1)-1 ) ! This is the name of the matrix to be partitioned
         DATA_80(2)(1:) = CARD2 ( COMMA_COL(1)+1 : COMMA_COL(2)-1 ) ! This is the name of the col partitioning vector
         DATA_80(3)(1:) = CARD2 ( COMMA_COL(2)+1 : DATA_END       ) ! This is the name of the row partitioning vector

         DO I=1,3                                          ! Get 16 char matrix names by stripping leading/trailing blanks
            IF (DATA_80(I)(1:) /= ' ') THEN
               DO J=1,EC_ENTRY_LEN
                  IF (DATA_80(I)(J:J) == ' ') THEN
                     CYCLE
                  ELSE
                     JBEG = J
                     EXIT
                  ENDIF
               ENDDO
               DO J=EC_ENTRY_LEN,JBEG,-1
                  IF (DATA_80(I)(J:J) == ' ') THEN
                     CYCLE
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
               DATA_16(I)(1:) = DATA_80(I)(JBEG:JBEG+15)
            ELSE
               DATA_16(I)(1:) = ' '
            ENDIF
         ENDDO

         FOUND = 'N'
         DO I=1,NUM_OU4_REQUESTS                           ! Set names of the matrix to be partitioned, the partitions and the vecs
            ! DATA_16(1) was upper-cased when reading the card.
            ! Here, upper-case the allowed names so that comparison is case-insensitive
            ! and thus allows KRRcb, MRRcb, KLR(t).
            IF (DATA_16(1) == TO_UPPER(ACT_OU4_MYSTRAN_NAMES(I))) THEN
               FOUND = 'Y'
               OU4_PART_MAT_NAMES(I,1) = DATA_16(1)
               OU4_PART_VEC_NAMES(I,1) = DATA_16(2)
               OU4_PART_VEC_NAMES(I,2) = DATA_16(3)
            ENDIF
         ENDDO

         IF (FOUND == 'N') THEN                            ! Matrix to be partitioned is not an OU4 matrix or an OU4 request
            IERR = IERR + 1
            WRITE(ERR,1040) DATA_16(1)
            WRITE(F06,1040) DATA_16(1)
         ENDIF

      ENDIF nerr

      IF ((DEBUG(197) == 2) .OR. (DEBUG(197) == 3)) THEN

         WRITE(F06,'(A)') 'Debug output from subr EC_PARTN'
         WRITE(F06,'(A)') '---------------------------------'
         WRITE(F06,*)
         WRITE(F06,99881) 'CARD1', CARD1
         WRITE(F06,*)
         WRITE(F06,99881) 'CARD2', CARD2
         WRITE(F06,*)

         WRITE(F06,99882) DATA_BEG
         WRITE(F06,*)

         WRITE(F06,99883) (COMMA_COL(I),I=1,2)
         WRITE(F06,*)

         IF (COMMA_COL(3) /= COMMA_COL(2)) THEN
            WRITE(F06,99884) COMMA_COL(3)
            WRITE(F06,*)
         ENDIF

         DO I=1,NUM_OU4_REQUESTS

            IF (ACT_OU4_MYSTRAN_NAMES(I)(1:16) == DATA_80(1)(1:16)) THEN
               WRITE(F06,99901) I, 'Matrix to be partitioned: OU4_PART_MAT_NAMES(I,1)', OU4_PART_MAT_NAMES(I,1)
               WRITE(F06,*)
               WRITE(F06,99901) I, 'Partitioning vec 1 name : OU4_PART_VEC_NAMES(I,1)', OU4_PART_VEC_NAMES(I,1)
               WRITE(F06,99901) I, 'Partitioning vec 2 name : OU4_PART_VEC_NAMES(I,2)', OU4_PART_VEC_NAMES(I,2)
            ENDIF

         ENDDO

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1025 FORMAT(' *ERROR  1025: THE PARTN ENTRY FORMAT REQUIRES ONE SLASH (/) AFTER THE SPECIFICATION OF THE PARTITIONING VECTORS')

 1031 FORMAT(' *ERROR  1031: EXEC CONTROL ENTRY OUTPUT4 MUST HAVE ',A,' COMMAS (AS IN PARTN A,CP,RP).')

 1032 FORMAT(' *ERROR  1032: EXEC CONTROL ENTRY OUTPUT4 CANNOT HAVE MORE THAN ',A,' COMMAS.')

 1039 FORMAT(' *ERROR  1039: A SLASH (/) MUST FOLLOW THE 3 PIECES OF DATA (MATRIX NAME, 2 PARTITION VECTORS) SEPARATED BY 2 COMMAS')

 1040 FORMAT(' *ERROR  1040: PARTN REQUEST HAS THE "',A,'" MATRIX REQUESTED TO BE PARTITIONED. HOWEVER, EITHER THIS IS NOT A'      &
                    ,/,14X,' VALID OUTPUT4 MATRIX NAME OR THE OUTPUT4 REQUEST FOR IT WAS NOT FOUND BEFORE THE PARTN REQUEST')

99881 FORMAT(A,':'/,A)

99882 FORMAT('Data begins in OUTPUT4 entry in col ',I3)

99883 FORMAT('COMMA_COL(1:2) = ',3(I4,','),I4,' (cols where commas exist)')

99884 FORMAT('There are more than 4 commas. A 5th one was found in col ',I4)

99901 FORMAT(' I =',I3,2X,A,2X,'"',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE EC_PARTN

   END MODULE EXECUTIVE_CONTROL_IO
