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

   MODULE BDF_FIELD_VALIDATION

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CHAR_FLD, I4FLD, R8FLD, IP6CHK, LEFT_ADJ_BDFLD, CRDERR, BD_IMBEDDED_BLANK

   CONTAINS

      SUBROUTINE CHAR_FLD ( JCARDI, IFLD, CHAR_INP )

! Reads a field of CHARACTER data that can be 1 to LEN(JCARDI) chars in length

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  IERRFL, FATAL_ERR

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)       :: JCARDI            ! The field of characters to read
      CHARACTER(LEN(JCARDI)), INTENT(OUT):: CHAR_INP          ! The character variable to read

      INTEGER(LONG), INTENT(IN)          :: IFLD              ! Field (2 - 9) of a Bulk Data card to read
      INTEGER(LONG)                      :: IOCHK             ! IOSTAT error value from READ

! **********************************************************************************************************************************
      CHAR_INP(1:) = ' '

      READ(JCARDI,'(A)',IOSTAT=IOCHK) CHAR_INP

      IF      (IOCHK < 0) THEN                             ! EOF/EOR during read

         IERRFL(IFLD) = 'Y'
         FATAL_ERR    = FATAL_ERR + 1

      ELSE IF (IOCHK == 0) THEN                            ! READ was OK

         IERRFL(IFLD) = 'N'

      ELSE IF (IOCHK > 0) THEN                             ! Error during READ

         IERRFL(IFLD) = 'Y'
         FATAL_ERR    = FATAL_ERR + 1

      ENDIF

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CHAR_FLD


      SUBROUTINE I4FLD ( JCARDI, IFLD, I4INP )

! Reads 8 column field of INTEGER*4 data

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  IERRFL, FATAL_ERR, JCARD_LEN, MAX_INTEGER_LEN

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: JCARDI            ! The field of 8 characters to read
      CHARACTER( 1*BYTE)              :: DEC_PT            ! 'Y'/'N' indicator of whether a decimal point was founr in JCARDI

      INTEGER(LONG), INTENT(IN)       :: IFLD              ! Field (2 - 9) of a Bulk Data card to read
      INTEGER(LONG), INTENT(OUT)      :: I4INP             ! The 4 byte integer value read
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: INTEGER_LEN       ! Length in digits of the integer in JCARDI
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from READ

! **********************************************************************************************************************************
      I4INP = 0
      IERRFL(IFLD) = 'N'

! Make sure integer number is not larger than (2^32)/2 = 2,147,483,648. For conservatism do not allow integers with more than 10
! digits and if it has 10 digits, make sure that the leading digit is <= 2

      INTEGER_LEN = JCARD_LEN
      DO I=JCARD_LEN,1,-1
         IF (JCARDI(I:I) == ' ') THEN
            INTEGER_LEN = INTEGER_LEN - 1
         ELSE
            EXIT
         ENDIF
      ENDDO

      IF (INTEGER_LEN > MAX_INTEGER_LEN) THEN              ! First, make sure integer has no more than 10 digits
         IERRFL(IFLD) = 'Y'
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1709) IFLD, INTEGER_LEN, MAX_INTEGER_LEN
         WRITE(F06,1709) IFLD, INTEGER_LEN, MAX_INTEGER_LEN

      ELSE IF (INTEGER_LEN == MAX_INTEGER_LEN) THEN

         IF ((JCARDI(1:1) == '1') .OR. (JCARDI(1:1) == '2')) THEN
            CONTINUE
         ELSE
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1710) IFLD, JCARDI
            WRITE(F06,1710) IFLD, JCARDI
         ENDIF

      ENDIF

! Read integer data in JCARDI

      READ(JCARDI,'(I8)',IOSTAT=IOCHK) I4INP

      IF (IOCHK /= 0) THEN
         IERRFL(IFLD) = 'Y'
         FATAL_ERR    = FATAL_ERR + 1
      ENDIF

! Scan to make sure there was not a decimal point. Don't set IERRFL, since an error message is written here.

      IF (JCARDI /= '        ') THEN
         DEC_PT = 'N'
         DO I=1,JCARD_LEN
            IF (JCARDI(I:I) == '.') THEN
               DEC_PT = 'Y'
               EXIT
            ENDIF
         ENDDO
         IF (DEC_PT == 'Y') THEN
            IERRFL(IFLD) = 'Y'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1700) IFLD
            WRITE(F06,1700) IFLD
         ENDIF
      ENDIF

      RETURN

! **********************************************************************************************************************************
 1700 FORMAT(' *ERROR  1700: A  DECIMAL POINT WAS FOUND IN WHAT IS SUPPOSED TO BE A INTEGER NUMBER IN FIELD ',I3,' OF THE',        &
                           ' PREVIOUS BULK DATA CARD')

 1709 FORMAT(' *ERROR  1709: FIELD ',I3,' OF THE PREVIOUS ENTRY HAS AN INTEGER WITH ',I3,' DIGITS. ALL INTEGERS MUST HAVE NO MORE',&
                           ' THAN ',I3,' DIGITS')

 1710 FORMAT(' *ERROR  1710: FIELD ',I3,' OF THE PREVIOUS ENTRY HAS AN INTEGER = ',A,' GREATER THAN 2,000,000,000. NOT ALLOWED')

! **********************************************************************************************************************************

      END SUBROUTINE I4FLD


      SUBROUTINE R8FLD ( JCARDI, IFLD, R8INP )

! Reads 8 column field of REAL DOUBLE data

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  IERRFL, FATAL_ERR, JCARD_LEN
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: JCARDI            ! The field of 8 characters to read
      CHARACTER( 1*BYTE)              :: DEC_PT            ! 'Y'/'N' indicator of whether a decimal point was founr in JCARDI

      INTEGER(LONG), INTENT(IN)       :: IFLD              ! Field (2 - 9) of a Bulk Data card to read
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error value from READ

      REAL(DOUBLE) , INTENT(OUT)      :: R8INP             ! The 8 byte real value read

! **********************************************************************************************************************************
      R8INP = ZERO
      IERRFL(IFLD) = 'N'

      READ(JCARDI,'(F16.0)',IOSTAT=IOCHK) R8INP

      IF (IOCHK /= 0) THEN
         IERRFL(IFLD) = 'Y'
         FATAL_ERR    = FATAL_ERR + 1
      ENDIF

! Scan to make sure there was a decimal point. Don't set IERRFL, since an error message is written here.

      IF (JCARDI /= '        ') THEN
         DEC_PT = 'N'
         DO I=1,JCARD_LEN
            IF (JCARDI(I:I) == '.') THEN
               DEC_PT = 'Y'
               EXIT
            ENDIF
         ENDDO
         IF (DEC_PT == 'N') THEN
            IERRFL(IFLD) = 'Y'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1701) IFLD
            WRITE(F06,1701) IFLD
         ENDIF
      ENDIF


      RETURN

! **********************************************************************************************************************************
 1701 FORMAT(' *ERROR  1701: NO DECIMAL POINT WAS FOUND IN WHAT IS SUPPOSED TO BE A REAL    NUMBER IN FIELD ',I3,' OF THE',        &
                           ' PREVIOUS BULK DATA CARD')

! **********************************************************************************************************************************

      END SUBROUTINE R8FLD


      SUBROUTINE IP6CHK ( JCARDI, JCARDO, IP6TYP, TOTAL_NUM_DIGITS )

! Routine to check DOF component and PINFLG fields on Bulk Data cards to make sure that they contain valid entries:
! Output JCARDO is:

!   (1) If input JCARDI has valid 1,2,3,4,5, and/or 6 then output JCARDO has these integers in ascending order
!                        If   JCARDI = '  5 2  4'
!                        Then JCARDO = '245     '
!       NOTE: if a DOF number is repeated, that is OK (e.g. 1355 is the same as 135); however, a warning is printed

!   (2) If input JCARDI does not have valid 1-6, output JCARDO is left blank and an error is indicated in IP6TYP

! Output IP6TYP is a descriptor of JCARDI and is:

!   (1) 'BLANK   ' if JCARDI is blank
!   (2) 'COMP NOS' if JCARDI has 1,2,3,4,5, and/or 6's
!   (3) 'ZERO    ' if JCARDI has one zero
!   (4) 'ERROR   ' if JCARDI has anything else


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))    :: SUBR_NAME = 'IP6CHK'
      CHARACTER(LEN=JCARD_LEN), INTENT(IN):: JCARDI            ! Input 8 character field
      CHARACTER(8*BYTE), INTENT(OUT)      :: IP6TYP            ! Descriptor of JCARDI, see above
      CHARACTER(LEN(JCARDI)), INTENT(OUT) :: JCARDO            ! Output 8 character field, described above
      CHARACTER(LEN(JCARDI))              :: JCARDO_TMP        ! Output 8 character field, described above

      INTEGER(LONG), INTENT(OUT)          :: TOTAL_NUM_DIGITS  ! Total of NUM_DIGITS(I)
      INTEGER(LONG)                       :: I                 ! DO loop index
      INTEGER(LONG)                       :: NUM_DIGITS(6)     ! NUM_DIGITS(I) is a count of the num of digits found in JCARDI
      INTEGER(LONG)                       :: POSN              ! An position in JCARDO




! **********************************************************************************************************************************
! Initialize

      JCARDO(1:)     = ' '
      JCARDO_TMP(1:) = ' '
      IP6TYP(1:)     = ' '
      DO I=1,6
         NUM_DIGITS(I) = 0
      ENDDO

      TOTAL_NUM_DIGITS = 0

! Check for all blank input:

      IF (JCARDI == '        ') THEN
         IP6TYP = 'BLANK   '

! Check for only 1 zero in the field

      ELSE IF ((JCARDI == '0       ') .OR. (JCARDI == ' 0      ') .OR. (JCARDI == '  0     ') .OR. (JCARDI == '   0    ') .OR.  &
               (JCARDI == '    0   ') .OR. (JCARDI == '     0  ') .OR. (JCARDI == '      0 ') .OR. (JCARDI == '       0')) THEN
         IP6TYP = 'ZERO    '

! Check for digits 1, 2, 3, 4, 5 and/or 6 (blanks among them are OK). If anything else found, then error.

      ELSE

i_loop:  DO I = 1,JCARD_LEN

            IF (JCARDI(I:I) == ' ') CYCLE i_loop

            IF      (JCARDI(I:I) == '1') THEN

               IP6TYP = 'COMP NOS'
               JCARDO_TMP(1:1) = '1'
               NUM_DIGITS(1) = NUM_DIGITS(1) + 1

            ELSE IF (JCARDI(I:I) == '2') THEN

               IP6TYP = 'COMP NOS'
               JCARDO_TMP(2:2) = '2'
               NUM_DIGITS(2) = NUM_DIGITS(2) + 1

            ELSE IF (JCARDI(I:I) == '3') THEN

               IP6TYP = 'COMP NOS'
               JCARDO_TMP(3:3) = '3'
               NUM_DIGITS(3) = NUM_DIGITS(3) + 1

            ELSE IF (JCARDI(I:I) == '4') THEN

               IP6TYP = 'COMP NOS'
               JCARDO_TMP(4:4) = '4'
               NUM_DIGITS(4) = NUM_DIGITS(4) + 1

            ELSE IF (JCARDI(I:I) == '5') THEN

               IP6TYP = 'COMP NOS'
               JCARDO_TMP(5:5) = '5'
               NUM_DIGITS(5) = NUM_DIGITS(5) + 1

            ELSE IF (JCARDI(I:I) == '6') THEN

               IP6TYP = 'COMP NOS'
               JCARDO_TMP(6:6) = '6'
               NUM_DIGITS(6) = NUM_DIGITS(6) + 1

            ELSE

               IP6TYP = 'ERROR   '
               JCARDO = '        '

               RETURN

            ENDIF

         ENDDO i_loop

         POSN = 0
         DO I=1,JCARD_LEN
            IF (JCARDO_TMP(I:I) == ' ') THEN
               CYCLE
            ELSE
               POSN = POSN + 1
               JCARDO(POSN:POSN) = JCARDO_TMP(I:I)
            ENDIF
         ENDDO

      ENDIF

      IF ((NUM_DIGITS(1) > 1) .OR. (NUM_DIGITS(2) > 1) .OR. (NUM_DIGITS(3) > 1) .OR.                                               &
          (NUM_DIGITS(4) > 1) .OR. (NUM_DIGITS(5) > 1) .OR. (NUM_DIGITS(6) > 1)) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,1738)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,1738)
         ENDIF
      ELSE
         DO I=1,6
            TOTAL_NUM_DIGITS = TOTAL_NUM_DIGITS + NUM_DIGITS(I)
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1738 FORMAT(' *WARNING    : SOME DOF COMPONENT NUMBERS 1 THROUGH 6 ARE REPEATED ON ABOVE CARD')

! **********************************************************************************************************************************

      END SUBROUTINE IP6CHK


      SUBROUTINE LEFT_ADJ_BDFLD ( CHR_FLD )

! Shifts a character string so that it is left adjusted

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))       :: SUBR_NAME = 'LEFT_ADJ_BDFLD'
      CHARACTER(LEN=JCARD_LEN), INTENT(INOUT):: CHR_FLD           ! Char field to left adjust and return
      CHARACTER(LEN=JCARD_LEN)               :: TCHR_FLD          ! Temporary char field

      INTEGER(LONG)                          :: I                 ! DO loop index




! **********************************************************************************************************************************
      IF (CHR_FLD(1:1) == ' ') THEN                        ! We need to shift:

         TCHR_FLD(1:) = CHR_FLD(1:)                        ! Set temporary field to CHR8_FLD

         DO I = 2,JCARD_LEN                                ! Perform shift
            IF (CHR_FLD(I:I) /= ' ') THEN
               TCHR_FLD(1:) = CHR_FLD(I:)
               EXIT
            ENDIF
         ENDDO

         CHR_FLD(1:) = TCHR_FLD(1:)                        ! Reset CHR_FLD and return

      ENDIF



      RETURN

! **********************************************************************************************************************************


      END SUBROUTINE LEFT_ADJ_BDFLD


      SUBROUTINE CRDERR ( CARD )

! Prints Bulk Data card errors and warnings

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, IERRFL
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CRDERR'
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER( 1*BYTE)              :: CARD_ERR          ! = 'Y' if IERRFL is 'Y' for any Bulk Data card field

      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
      CARD_ERR = 'N'
      DO I=1,10
         IF (IERRFL(I) == 'Y') THEN
            CARD_ERR = 'Y'
            EXIT
         ENDIF
      ENDDO

      IF (CARD_ERR == 'Y') THEN
         IF (ECHO == 'NONE  ') THEN
            WRITE(ERR,101) CARD
            WRITE(F06,101) CARD
         ENDIF
         DO I=1,10
            IF (IERRFL(I) == 'Y') THEN
               WRITE(ERR,1702) I
               WRITE(F06,1702) I
            ENDIF
         ENDDO
      ENDIF

      DO I=1,10
         IERRFL(I) = 'N'
      ENDDO



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1702 FORMAT(' *ERROR  1702: FORMAT ERROR IN FIELD',I3,' OF PREVIOUS CARD')

! **********************************************************************************************************************************

      END SUBROUTINE CRDERR


      SUBROUTINE BD_IMBEDDED_BLANK ( JCARD, CF2, CF3, CF4, CF5, CF6, CF7, CF8, CF9 )

! Prepares message when some fields of a B.D card have imbedded blanks when they should not (but field can be completely blank)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, BLNK_SUB_NAM, JCARD_LEN
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BD_IMBEDDED_BLANK'
      CHARACTER(LEN=*), INTENT(IN)    :: JCARD(10)         ! The 10 fields of characters making up CARD
      CHARACTER( 1*BYTE)              :: IMB_BLANK(2:9)    ! 'Y'/'N' indicator of whether fields 2-9 have imbedded blanks

      INTEGER(LONG), INTENT(IN)       :: CF2               ! = 2 if field 2 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)       :: CF3               ! = 3 if field 3 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)       :: CF4               ! = 4 if field 4 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)       :: CF5               ! = 5 if field 5 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)       :: CF6               ! = 6 if field 6 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)       :: CF7               ! = 7 if field 7 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)       :: CF8               ! = 8 if field 8 is to be checked, or 0 otherwise
      INTEGER(LONG), INTENT(IN)       :: CF9               ! = 9 if field 9 is to be checked, or 0 otherwise
      INTEGER(LONG)                   :: CHK_FLD(2:9)      ! Array containing CF2 through CF9
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: JCARDI_BEG        ! Position where data begins in one JCARD
      INTEGER(LONG)                   :: JCARDI_END        ! Position where data ends   in one JCARD
      INTEGER(LONG)                   :: NUMBER(2:9)       ! Number of imbedded blanks found in a Bulk Data card field




! **********************************************************************************************************************************
! Load CF2 through CF9 into array CHK_FLD

      CHK_FLD(2) = CF2
      CHK_FLD(3) = CF3
      CHK_FLD(4) = CF4
      CHK_FLD(5) = CF5
      CHK_FLD(6) = CF6
      CHK_FLD(7) = CF7
      CHK_FLD(8) = CF8
      CHK_FLD(9) = CF9

! Initialize IMB_BLANK array

      DO I=2,9
         IMB_BLANK(I) = 'N'
      ENDDO

! Find beginning and end of data in

! Check fields for any imbedded blanks and set error if so

      DO I=2,9

         JCARDI_BEG = 1
         JCARDI_END = JCARD_LEN

         NUMBER(I) = 0

         IF      (CHK_FLD(I) == 0) THEN                    ! We don't want to check field I

            CYCLE

         ELSE IF (CHK_FLD(I) == I) THEN                    ! We do want to check field I

            IF (JCARD(I)(1:) /= ' ') THEN                  ! Don't need to check blank fields

j_do1:         DO J=1,JCARD_LEN                            ! Find where data begins in this field
                  IF (JCARD(I)(J:J) == ' ') THEN
                     CYCLE j_do1
                  ELSE
                     JCARDI_BEG = J
                     EXIT j_do1
                  ENDIF
               ENDDO j_do1

j_do2:         DO J=JCARD_LEN,1,-1                         ! Find where data ends in this field
                  IF (JCARD(I)(J:J) == ' ') THEN
                     CYCLE j_do2
                  ELSE
                     JCARDI_END = J
                     EXIT j_do2
                  ENDIF
               ENDDO j_do2

               IMB_BLANK(I) = 'N'                          ! Check between data begin/end for blanks
               DO J=JCARDI_BEG,JCARDI_END
                  IF(JCARD(I)(J:J) == ' ') THEN
                     IMB_BLANK(I) = 'Y'
                     NUMBER(I)    = NUMBER(I) + 1
                  ENDIF
               ENDDO

            ENDIF

         ELSE                                              ! Coding error, CHK_FLD(I) must be 0 or I

            WRITE(ERR,1121) SUBR_NAME, CHK_FLD(I)
            WRITE(F06,1121) SUBR_NAME, CHK_FLD(I)
            CALL OUTA_HERE ( 'Y' )

         ENDIF

      ENDDO

! Write error if fields checked have imbedded blanks

      DO I=2,9
         IF(IMB_BLANK(I) == 'Y') THEN
            WRITE(ERR,1122) NUMBER(I),I,JCARD(1),JCARD(2)
            WRITE(F06,1122) NUMBER(I),I,JCARD(1),JCARD(2)
            FATAL_ERR = FATAL_ERR + 1
            CYCLE
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
 1121 FORMAT(' *ERROR  1121: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CHK_FLD(I) MUST BE EITHER 0 OR I BUT IS ',I8)

 1122 FORMAT(' *ERROR  1122: THERE WERE ',I2,' IMBEDDED BLANKS FOUND IN FIELD ',I2,' OF BULK DATA ENTRY ',A,' ',A,                 &
                                            '. IMBEDDED BLANKS NOT ALLOWED HERE')


! **********************************************************************************************************************************

      END SUBROUTINE BD_IMBEDDED_BLANK

   END MODULE BDF_FIELD_VALIDATION
