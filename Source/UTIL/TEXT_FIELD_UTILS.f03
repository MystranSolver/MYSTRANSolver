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

   MODULE TEXT_FIELD_UTILS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CARD_FLDS_NOT_BLANK, CONVERT_INT_TO_CHAR, FMT_ES14_6, FMT_I8_RJ, GET_CHAR_STRING_END, PARSE_CHAR_STRING, REAL_DATA_TO_C8FLD, TO_UPPER

   CONTAINS

      SUBROUTINE CARD_FLDS_NOT_BLANK ( JCARD, FLD2, FLD3, FLD4, FLD5, FLD6, FLD7, FLD8, FLD9 )

! Prepares message when some fields of a Bulk data card that should be blank, aren't

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, JCARD_LEN, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SUPWARN

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))    :: SUBR_NAME = 'CARD_FLDS_NOT_BLANK'
      CHARACTER(LEN=JCARD_LEN), INTENT(IN):: JCARD(10)         ! The 10 fields of 8 characters making up CARD
      CHARACTER( 1*BYTE)                  :: COMMENT           ! 'Y' or 'N' depending on whether non-blank fields are a comment
      CHARACTER( 8*BYTE)                  :: MSSG8             ! Message with all fields that are not blank that should be blank
      CHARACTER( 1*BYTE)                  :: MSSG1             ! Message that has the field number in it

      INTEGER(LONG), INTENT(IN)           :: FLD2              ! Refers to field 2 of a B.D. card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)           :: FLD3              ! Refers to field 3 of a B.D. card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)           :: FLD4              ! Refers to field 4 of a B.D. card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)           :: FLD5              ! Refers to field 5 of a B.D. card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)           :: FLD6              ! Refers to field 6 of a B.D. card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)           :: FLD7              ! Refers to field 7 of a B.D. card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)           :: FLD8              ! Refers to field 8 of a B.D. card. If /= 0, then check this field
      INTEGER(LONG), INTENT(IN)           :: FLD9              ! Refers to field 9 of a B.D. card. If /= 0, then check this field
      INTEGER(LONG)                       :: ALL_FLDS(2:9)     ! Array of the FLDi (2 through 9)
      INTEGER(LONG)                       :: I,J               ! Do loop indices




! **********************************************************************************************************************************
! Set ALL_FLDS

      ALL_FLDS(2) = FLD2
      ALL_FLDS(3) = FLD3
      ALL_FLDS(4) = FLD4
      ALL_FLDS(5) = FLD5
      ALL_FLDS(6) = FLD6
      ALL_FLDS(7) = FLD7
      ALL_FLDS(8) = FLD8
      ALL_FLDS(9) = FLD9

! Scan for 1st non-blank field that should not have data and see if the 1st non-blank character is '$'.
! If it is, then the non-blank fields have a comment, so return

      COMMENT = 'N'
i_do: DO I=2,9
         IF (ALL_FLDS(I) /= I) THEN                        ! CYCLE until we find 1st field to check
            CYCLE
         ELSE                                              ! Found a field that should be blank (unless it begins a comment)
j_do:       DO J=1,JCARD_LEN                               ! CYCLE through chars of this field to see if 1st non-blank char is '$'
               IF (JCARD(I)(J:J) == ' ') THEN
                  CYCLE j_do
               ELSE
                  IF (JCARD(I)(J:J) == '$') THEN           ! Found a '$' (with blanks preceeding it) so we found a comment
                     COMMENT = 'Y'
                     EXIT i_do                             ! Set COMMENT = 'Y' and EXIT outer loop
                  ENDIF
               ENDIF
            ENDDO j_do
         ENDIF
      ENDDO i_do

! Issue warning if fields not blank that should be.

      IF (COMMENT == 'N') THEN

         MSSG8 = '        '
         IF ((JCARD(2)(1:) /= ' ') .AND. (FLD2 == 2)) THEN
            MSSG1 = '2'
            MSSG8 =                  MSSG1(1:1) // MSSG8(2:8)
         ENDIF
         IF ((JCARD(3)(1:) /= ' ') .AND. (FLD3 == 3)) THEN
            MSSG1 = '3'
            MSSG8 = MSSG8(1:1) // MSSG1(1:1) // MSSG8(3:8)
         ENDIF
         IF ((JCARD(4)(1:) /= ' ') .AND. (FLD4 == 4)) THEN
            MSSG1 = '4'
            MSSG8 = MSSG8(1:2) // MSSG1(1:1) // MSSG8(4:8)
         ENDIF
         IF ((JCARD(5)(1:) /= ' ') .AND. (FLD5 == 5)) THEN
            MSSG1 = '5'
            MSSG8 = MSSG8(1:3) // MSSG1(1:1) // MSSG8(5:8)
         ENDIF
         IF ((JCARD(6)(1:) /= ' ') .AND. (FLD6 == 6)) THEN
            MSSG1 = '6'
            MSSG8 = MSSG8(1:4) // MSSG1(1:1) // MSSG8(6:8)
         ENDIF
         IF ((JCARD(7)(1:) /= ' ') .AND. (FLD7 == 7)) THEN
            MSSG1 = '7'
            MSSG8 = MSSG8(1:5) // MSSG1(1:1) // MSSG8(7:8)
         ENDIF
         IF ((JCARD(8)(1:) /= ' ') .AND. (FLD8 == 8)) THEN
            MSSG1 = '8'
            MSSG8 = MSSG8(1:6) // MSSG1(1:1) // MSSG8(8:8)
         ENDIF
         IF ((JCARD(9)(1:) /= ' ') .AND. (FLD9 == 9)) THEN
            MSSG1 = '9'
            MSSG8 = MSSG8(1:7) // MSSG1(1:1)
         ENDIF
         IF (MSSG8 /= '        ') THEN
            WARN_ERR = WARN_ERR+1
            WRITE(ERR,1726) MSSG8
            IF (SUPWARN == 'N') THEN
               WRITE(F06,1726) MSSG8
            ENDIF
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1726 FORMAT(' *WARNING    : FIELD(s) ',A8,' ON PREVIOUS ENTRY SHOULD BE BLANK AND ARE IGNORED')

! **********************************************************************************************************************************

      END SUBROUTINE CARD_FLDS_NOT_BLANK


      SUBROUTINE CONVERT_INT_TO_CHAR ( INT_NUM, CHAR_VALUE )

! Convert an integer 1, 2, 3, 4, 5 or 6 to character '1', '2' ... '6'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CONVERT_INT_TO_CHAR'
      CHARACTER(1*BYTE), INTENT(OUT)  :: CHAR_VALUE        ! If INT_NUM = 1, then CHAR_VALUE = '1', etc

      INTEGER(LONG), INTENT(IN)       :: INT_NUM           ! Integer 1, 2, 3, 4, 5 O5 6




! **********************************************************************************************************************************
! Initialize outputs

      CHAR_VALUE = ' '

! Make sure that INT_NUM is in the range 1-6

      IF ((INT_NUM < 1) .OR. (INT_NUM > 6)) THEN
         WRITE(ERR,935) INT_NUM
         WRITE(F06,935) INT_NUM
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Convert the integer to character

      IF      (INT_NUM == 1) THEN
         CHAR_VALUE = '1'
      ELSE IF (INT_NUM == 2) THEN
         CHAR_VALUE = '2'
      ELSE IF (INT_NUM == 3) THEN
         CHAR_VALUE = '3'
      ELSE IF (INT_NUM == 4) THEN
         CHAR_VALUE = '4'
      ELSE IF (INT_NUM == 5) THEN
         CHAR_VALUE = '5'
      ELSE IF (INT_NUM == 6) THEN
         CHAR_VALUE = '6'
      ENDIF



      RETURN

! **********************************************************************************************************************************
  935 FORMAT(' *ERROR   935: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ILLEGAL INTEGER VALUE = ',I8,' INPUT. VALUE MUST BE AN INTEGER IN THE RANGE 1-6')

! **********************************************************************************************************************************

      END SUBROUTINE CONVERT_INT_TO_CHAR


      SUBROUTINE FMT_ES14_6 ( V, OUT )

! Hand-rolled formatter that produces the same 14-character output as Fortran's 1ES14.6 edit descriptor,
! but bypasses the (relatively expensive) internal WRITE machinery. Used in the LK9 output pipeline to
! accelerate writing of the F06 file. Exact zeros are emitted as ' 0.000000E+00' to match Fortran's
! native 1ES14.6 output byte-for-byte. (Callers that prefer the WRT_REAL_TO_CHAR_VAR style '  0.0         '
! substitution must perform that replacement themselves.)
!
! Layout of OUT (positions 1..14):
!   1 : leading space
!   2 : sign ('-' or ' ')
!   3 : single mantissa digit
!   4 : '.'
!   5-10 : 6 fractional digits
!   11 : 'E'
!   12 : exponent sign ('+' or '-')
!   13-14 : 2-digit exponent
!
! Assumes |decimal exponent| < 100. Values exceeding that range are rare in FEA stress/strain output;
! if encountered the routine falls back to a Fortran internal WRITE so the field is still well formed.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: V                    ! Real value to format
      CHARACTER(14*BYTE), INTENT(OUT) :: OUT                  ! 14-char formatted result

      REAL(DOUBLE)                    :: AV, M
      INTEGER(LONG)                   :: E10, IM, K
      INTEGER(LONG)                   :: DIG(0:6)
      CHARACTER(1*BYTE)               :: ESIGN
      LOGICAL                         :: NEG

      INTEGER(LONG), PARAMETER        :: ZERO_CHAR = IACHAR('0')

! **********************************************************************************************************************************
      IF (V == ZERO) THEN
         OUT = '  0.000000E+00'
         RETURN
      ENDIF

      NEG = V < ZERO
      AV  = ABS(V)

      E10 = FLOOR(LOG10(AV))
      M   = AV * (10.0D0 ** (-E10))

      ! Guard against floating-point rounding in LOG10/FLOOR that could put M outside [1,10).
      IF (M < 1.0D0) THEN
         M   = M * 10.0D0
         E10 = E10 - 1
      ELSE IF (M >= 10.0D0) THEN
         M   = M * 0.1D0
         E10 = E10 + 1
      ENDIF

      IM = NINT(M * 1.0D6, KIND=LONG)
      IF (IM >= 10000000) THEN          ! carry from rounding 9.9999995 -> 10.000000
         IM  = IM / 10
         E10 = E10 + 1
      ENDIF

      ! Fall back to Fortran formatting for the rare |E10| >= 100 case so the field still has 14 chars.
      IF ((E10 >= 100) .OR. (E10 <= -100)) THEN
         WRITE(OUT,'(1ES14.6)') V
         RETURN
      ENDIF

      DO K = 6, 0, -1
         DIG(K) = MOD(IM, 10_LONG)
         IM     = IM / 10
      ENDDO

      IF (E10 < 0) THEN
         ESIGN = '-'
         E10   = -E10
      ELSE
         ESIGN = '+'
      ENDIF

      OUT(1:1)   = ' '
      IF (NEG) THEN
         OUT(2:2) = '-'
      ELSE
         OUT(2:2) = ' '
      ENDIF
      OUT(3:3)   = ACHAR(ZERO_CHAR + DIG(0))
      OUT(4:4)   = '.'
      OUT(5:5)   = ACHAR(ZERO_CHAR + DIG(1))
      OUT(6:6)   = ACHAR(ZERO_CHAR + DIG(2))
      OUT(7:7)   = ACHAR(ZERO_CHAR + DIG(3))
      OUT(8:8)   = ACHAR(ZERO_CHAR + DIG(4))
      OUT(9:9)   = ACHAR(ZERO_CHAR + DIG(5))
      OUT(10:10) = ACHAR(ZERO_CHAR + DIG(6))
      OUT(11:11) = 'E'
      OUT(12:12) = ESIGN
      OUT(13:13) = ACHAR(ZERO_CHAR + E10 / 10)
      OUT(14:14) = ACHAR(ZERO_CHAR + MOD(E10, 10_LONG))

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE FMT_ES14_6


      SUBROUTINE FMT_I8_RJ ( V, OUT )

! Right-justify a signed integer into an 8-character field, padded with spaces. Produces the same output
! as Fortran's I8 edit descriptor but avoids the internal WRITE machinery. Used in the LK9 output
! pipeline alongside FMT_ES14_6 to accelerate F06 line assembly. Values with more than 8 significant
! digits (counting the sign) fall back to a Fortran internal WRITE so the field stays well formed.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: V                    ! Integer value to format
      CHARACTER(8*BYTE), INTENT(OUT)  :: OUT                  ! 8-char right-justified result

      INTEGER(LONG)                   :: X, K
      INTEGER(LONG), PARAMETER        :: ZERO_CHAR = IACHAR('0')

! **********************************************************************************************************************************
      OUT = '        '

      IF (V == 0) THEN
         OUT(8:8) = '0'
         RETURN
      ENDIF

      ! Overflow guard: an I8 field holds at most 8 chars (sign + 7 digits for negatives, 8 digits for non-negatives).
      IF ((V > 99999999_LONG) .OR. (V < -9999999_LONG)) THEN
         WRITE(OUT,'(I8)') V
         RETURN
      ENDIF

      X = ABS(V)
      K = 8
      DO WHILE ((X > 0) .AND. (K >= 1))
         OUT(K:K) = ACHAR(ZERO_CHAR + MOD(X, 10_LONG))
         X        = X / 10
         K        = K - 1
      ENDDO
      IF ((V < 0) .AND. (K >= 1)) OUT(K:K) = '-'

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE FMT_I8_RJ


      SUBROUTINE GET_CHAR_STRING_END ( CHAR_STRING, IEND )

! Searches integer array ARRAY to find column where data ends (IEND)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_CHAR_STRING_END'
      CHARACTER(LEN=*) , INTENT(IN)   :: CHAR_STRING       ! String to get ending of

      INTEGER(LONG)    , INTENT(OUT)  :: IEND              ! Col where CHAR_STRING stops having non blanks
      INTEGER(LONG)                   :: I                 ! DO loop index




! **********************************************************************************************************************************
      IEND = LEN(CHAR_STRING)
      DO I=LEN(CHAR_STRING),1,-1
         IF (CHAR_STRING(I:I) /= ' ') THEN
            IEND = I
            EXIT
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE GET_CHAR_STRING_END


      SUBROUTINE GET_FORMATTED_INTEGER ( INT, CHAR_INT, NUM_CHARS, NUM_DIGITS )

! Converts an integer to a character value with comma format (e.g. 12345 becomes char value 12,345) and writes result to unit UNT

      USE PENTIUM_II_KIND, ONLY             :  BYTE, LONG
      USE SCONTR, ONLY                      :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                      :  TSEC

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER              :: WORD_LEN  = 13    ! Length of character string that INT will be entered into
!                                                                  This allows for a number up to 9,999,999,999

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))      :: SUBR_NAME = 'GET_FORMATTED_INTEGER'
      CHARACTER(WORD_LEN*BYTE), INTENT(OUT) :: CHAR_INT          ! Integer formatted to have comma's (36879 becomes 36,879)
      CHARACTER(WORD_LEN*BYTE)              :: TEMP_CHAR_INT     ! Temporary value of CHAR_INT


      INTEGER(LONG), INTENT(IN)             :: INT               ! Integer to be converted to formated value in CHAR_INT
      INTEGER(LONG), INTENT(OUT)            :: NUM_CHARS         ! Num of non blank chars in CHAR_INT after formatting w/ commas
      INTEGER(LONG), INTENT(OUT)            :: NUM_DIGITS        ! Number of digits in INT
      INTEGER(LONG)                         :: I,J,K             ! DO loop indices or counters



! **********************************************************************************************************************************
! Initialize

      NUM_CHARS = WORD_LEN

      TEMP_CHAR_INT(1:) = ' '
      CHAR_INT(1:)      = ' '
      IF (WORD_LEN == 13) THEN                             ! This code to make sure the format of the WRITE is same length as
         WRITE(TEMP_CHAR_INT,'(I13)') INT                  ! above declaratioin for INTEGER PARAMETER WORD_LEN
      ELSE
         CHAR_INT(1:) = '*'
         RETURN
      ENDIF

! Find out haw many digits are in INT

      DO I=WORD_LEN,1,-1
         IF (TEMP_CHAR_INT(I:I) == ' ') THEN
            NUM_DIGITS = WORD_LEN - I
            EXIT
         ENDIF
      ENDDO

! Move digits from TEMP_CHAR_INT to CHAR_INT inserting commas

      IF (NUM_DIGITS > 3) THEN

         K = WORD_LEN
         DO I=WORD_LEN,WORD_LEN-NUM_DIGITS,-3

            DO J=1,3
               CHAR_INT(K-J+1:K-J+1) = TEMP_CHAR_INT(I-J+1:I-J+1)
            ENDDO
            K = K - 3

            IF (TEMP_CHAR_INT(I-3:I-3) /= ' ') THEN
               CHAR_INT(K:K) = ','
               K = K - 1
            ENDIF

         ENDDO

      ELSE

         CHAR_INT(1:) = TEMP_CHAR_INT(1:)

      ENDIF
!xxLeft adjust CHAR_INT
!xx
!xx   TEMP_CHAR_INT(1:) = ' '
!xx   IF (CHAR_INT(1:1) == ' ') THEN                       ! We need to shift:
!xx
!xx      TEMP_CHAR_INT(1:) = CHAR_INT(1:)                  ! Set temporary field to CHR8_FLD
!xx
!xx      DO I = 2,WORD_LEN                                       ! Perform shift
!xx         IF (CHAR_INT(I:I) /= ' ') THEN
!xx            TEMP_CHAR_INT(1:) = CHAR_INT(I:)
!xx            EXIT
!xx         ENDIF
!xx      ENDDO
!xx
!xx      CHAR_INT(1:) = TEMP_CHAR_INT(1:)                  ! Reset CHR_FLD and return
!xx
!xx   ENDIF
!xx
!xxCount nonblank characters in CHAR_INT
!xx
!xx   NUM_CHARS = WORD_LEN
!xx   DO I=WORD_LEN,1,-1
!xx      IF ((CHAR_INT(I:I) == ' ') .OR. (CHAR_INT(I:I) == ',')) THEN
!xx         CYCLE
!xx      ELSE
!xx         NUM_CHARS = WORD_LEN - I
!xx      write(f06,'(a,i3,a,a,a,i3)') ' I, CHAR_INT(I:I)= ', i, '  "', char_int(i:i), '"'
!xx      ENDIF
!xx   ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE GET_FORMATTED_INTEGER



      SUBROUTINE PARSE_CHAR_STRING ( CHAR_STRING, STRING_LEN, MAX_WORDS, MWLEN, NUM_WORDS, WORDS, IERR )

! Parses a character string whose words are separated by blanks and/or commas into a 1-D array, WORDS, of the words in the string.
! For example, if CHAR_STRING is the string (without quotes): "  SORT1 ,   REAL,PRINT  ,,  ,  QAZ  VONMISES", then this subr
! parses CHAR_STRING into 5 words in the array WORDS:

!                            WORDS(1) = SORT1
!                            WORDS(2) = REAL
!                            WORDS(3) = PRINT
!                            WORDS(4) = QAZ
!                            WORDS(5) = VONMISES

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER         :: MAX_LEN_BAD_WRD=150!

      INTEGER(LONG), INTENT(IN)        :: MAX_WORDS          ! Dim of WORDS (number of words in CHAR_STRING cannot exceed this)
      INTEGER(LONG), INTENT(IN)        :: MWLEN              ! Maximum length, in characters, of the entries in array WORDS

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)) :: SUBR_NAME = 'PARSE_CHAR_STRING'
      CHARACTER(LEN=*)    , INTENT(IN) :: CHAR_STRING        ! Character string to be parsed
      CHARACTER(LEN=MWLEN), INTENT(OUT):: WORDS(MAX_WORDS)   ! Array of the words parsed from CHAR_STRING.
      CHARACTER(LEN=MWLEN)             :: WORD               ! One word that will go into array WORDS

      INTEGER(LONG), INTENT(IN)        :: STRING_LEN         ! Length, in characters, of CHAR_STRING
      INTEGER(LONG), INTENT(OUT)       :: IERR               ! Error designator
      INTEGER(LONG), INTENT(OUT)       :: NUM_WORDS          ! Number of distinct words in CHAR_STRING
      INTEGER(LONG)                    :: NUM         = 0    ! Lesser of NUM_WORDS and MAX_WORDS
      INTEGER(LONG)                    :: CHAR_COUNT         ! Index into CHAR_STRING to a character in that string (not ' ' or ',')
      INTEGER(LONG)                    :: I,J                ! DO loop indices
      INTEGER(LONG)                    :: WORD_LEN           ! Length of one of the words in CHAR_STRING (must be <= MWLEN)



! **********************************************************************************************************************************
! Initialize

      DO I=1,MAX_WORDS
         WORDS(I)(1:) = ' '
      ENDDO

      IERR = 0

! Call debug code if requested

      IF (DEBUG(173) >= 1) CALL DEBUG_PARSE_CHR_STRNG ( 1, '173' )

! Parse words from CHAR_STRING

      CHAR_COUNT = 0
      NUM_WORDS  = 0
nwrds:DO

         CHAR_COUNT = CHAR_COUNT + 1

         IF (CHAR_COUNT > STRING_LEN) THEN

            EXIT nwrds

         ELSE

            IF ((CHAR_STRING(CHAR_COUNT:CHAR_COUNT) == ',') .OR. (CHAR_STRING(CHAR_COUNT:CHAR_COUNT) == ' ')) THEN

               CYCLE nwrds

            ELSE

               WORD_LEN  = 0
               NUM_WORDS = NUM_WORDS + 1

               WORD(1:) = ' '
one_wrd:       DO
                  WORD_LEN = WORD_LEN + 1
                  IF (DEBUG(173) > 1) CALL DEBUG_PARSE_CHR_STRNG ( 2, '173' )
                  IF (WORD_LEN > MWLEN) THEN
                     IERR = 1
                     CALL PARSE_CHAR_STRING_MSG ( 1, WORD )
                     CYCLE one_wrd
                  ENDIF
                  WORD(WORD_LEN:WORD_LEN) = CHAR_STRING(CHAR_COUNT:CHAR_COUNT)
                  CHAR_COUNT = CHAR_COUNT + 1
                  IF (CHAR_COUNT > STRING_LEN+1) THEN      ! Need +1 in order to get last WORD put into array WORDS
                     EXIT nwrds
                  ENDIF
                  IF ((CHAR_STRING(CHAR_COUNT:CHAR_COUNT) == ',') .OR. (CHAR_STRING(CHAR_COUNT:CHAR_COUNT) == ' ')) THEN
                     EXIT one_wrd
                  ELSE
                     IF (WORD_LEN < MWLEN) THEN
                        CYCLE one_wrd
                     ELSE IF (WORD_LEN == MWLEN) then      ! The word is already as long as is allowed. See if it continues
j_do:                   DO J=1,MAX_LEN_BAD_WRD             ! Loop trying to get to end of this long word to see if there are others

                           CHAR_COUNT = CHAR_COUNT + 1
                           IF ((CHAR_STRING(CHAR_COUNT:CHAR_COUNT) == ',') .OR. (CHAR_STRING(CHAR_COUNT:CHAR_COUNT) == ' ')) THEN
                              EXIT one_wrd
                           ENDIF

                           IF (J < MAX_LEN_BAD_WRD) THEN
                              IERR = 1
                              CALL PARSE_CHAR_STRING_MSG ( 1, WORD )
                              CYCLE j_do
                           ENDIF

                           IF (J == MAX_LEN_BAD_WRD) THEN  ! Give up and get out of this loop. This word is absurdly too long
                              EXIT one_wrd
                           ENDIF

                        ENDDO j_do
                     ELSE
                        EXIT one_wrd
                     ENDIF
                  ENDIF
               ENDDO one_wrd

               IF (NUM_WORDS > MAX_WORDS) THEN
                  IERR = 2
                  CALL PARSE_CHAR_STRING_MSG ( 2, WORD )
                  EXIT nwrds
               ENDIF

               WORDS(NUM_WORDS) = WORD

               IF (DEBUG(173) > 1) CALL DEBUG_PARSE_CHR_STRNG ( 3, '173' )

            ENDIF

         ENDIF

      ENDDO nwrds

      IF (DEBUG(173) >= 1) CALL DEBUG_PARSE_CHR_STRNG ( 4, '173' )

      IF (IERR > 0) THEN
         WRITE(F06,9998)
         NUM = MAX_WORDS
         IF (NUM_WORDS <= MAX_WORDS) NUM = NUM_WORDS
         DO I=1,NUM
            WRITE(F06,9999) I, WORDS(I)
         ENDDO
         WRITE(F06,*)
      ENDIF



! **********************************************************************************************************************************
 9998 FORMAT('               THE WORDS FROM THE STRING ARE PRINTED BELOW:',/)

 9999 FORMAT(I16,2X,'"',A,'"')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE PARSE_CHAR_STRING_MSG ( OPT, WORD )

      IMPLICIT NONE

      CHARACTER(LEN=MWLEN)            :: WORD              ! One word that will go into array WORDS

      INTEGER(LONG)                   :: OPT               ! Tells which option to use here

      WARN_ERR = WARN_ERR + 1

      IF (OPT == 1) THEN
         WRITE(ERR,961) MWLEN, WORD
         WRITE(F06,961) MWLEN, WORD
      ELSE IF (OPT == 2) THEN
         WRITE(ERR,962) MAX_WORDS
         WRITE(F06,962) MAX_WORDS
      ENDIF

! **********************************************************************************************************************************
  961 FORMAT(' *WARNING    : THE LENGTH OF A WORD IN ABOVE STRING EXCEEDED ',I5,' CHARACTERS. THAT WORD, ','"',A,'...", IS IGNORED')

  962 FORMAT(' *WARNING    : THE NUMBER OF WORDS IN ABOVE CHARACTER STRING EXCEEDED',I5,'. REMAINING WORDS IGNORED')

! **********************************************************************************************************************************

      END SUBROUTINE PARSE_CHAR_STRING_MSG

! ##################################################################################################################################

      SUBROUTINE DEBUG_PARSE_CHR_STRNG ( WHAT, DEB_NUM )

      CHARACTER( 3*BYTE)              :: DEB_NUM           !

      INTEGER(LONG)                   :: WHAT

! **********************************************************************************************************************************
      IF (WHAT == 1) THEN
         WRITE(F06,10)
         WRITE(F06,11) DEB_NUM
         WRITE(F06,12) CHAR_STRING(1:STRING_LEN)
         WRITE(F06,*)
      ENDIF

      IF (WHAT == 2) THEN
         WRITE(F06,21) NUM_WORDS, CHAR_COUNT, CHAR_STRING(CHAR_COUNT:CHAR_COUNT), word_len
      ENDIF

      IF (WHAT == 3) THEN
         WRITE(F06,*)
      ENDIF

      IF (WHAT == 4) THEN
         WRITE(F06,41) NUM_WORDS, STRING_LEN
         NUM = NUM_WORDS
         IF (NUM_WORDS > MAX_WORDS) NUM = MAX_WORDS
         DO I=1,NUM
            WRITE(F06,42) I, WORDS(I)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,10)
      ENDIF

! **********************************************************************************************************************************
   10 FORMAT('******************************************************************************')

   11 FORMAT('In subr PARSE_CHAR_STRING for DEBUG(', A,')' /)

   12 FORMAT('   CHAR_STRING = "',A,'"')

   21 FORMAT('   Word num, I, CHAR_STRING(I:I) = ',I3,I4,2X,A, i8)

   41 FORMAT('   There are ',I2,' word(s) in the ',I3,' character variable CHAR_STRING:',/,                                        &
             '   --------------------------------------------------------------')

   42 FORMAT('     Word ',I2,' = "',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE DEBUG_PARSE_CHR_STRNG

      END SUBROUTINE PARSE_CHAR_STRING



      SUBROUTINE REAL_DATA_TO_C8FLD ( REAL_INP, CHAR8_OUT )

! Converts a real number to 8 character field:

!     -1.45367+05 converts to -1.453+5
!     -1.65743+12 converts to -1.65+12
!     +1.43678+05 converts to 1.4367+5
!     +1.43567+12 converts to 1.435+12

! This is used to create 8 char Bulk Data fields for real numbers with as many significant digits as possible. This subr is called
! when subr WRITE_PCOMP_EQUIV writes the PSHELL equivalent of a PCOMP to the F06 file (based on user request via Bulk Data PARAM
! PCOMPEQ

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      IMPLICIT NONE

      CHARACTER( 8*BYTE), INTENT(OUT) :: CHAR8_OUT         ! 8 character representation of REAL_INP
      CHARACTER(11*BYTE)              :: TEMP_CHAR         ! Temporary char field to store REAL_INP in 1ES11.4 format

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IBEG,IEND         ! Locations in char strings

      REAL(DOUBLE), INTENT(IN)        :: REAL_INP          ! Double precision real input number to be converted

! **********************************************************************************************************************************
      WRITE(TEMP_CHAR,'(1ES11.4)') REAL_INP

      CHAR8_OUT(1:) = ' '

      IBEG = 1
      IF (TEMP_CHAR(1:1) == '-') THEN
         CHAR8_OUT(1:1) = '-'
         IBEG = IBEG + 1
      ENDIF
      CHAR8_OUT(IBEG:IBEG) = TEMP_CHAR(2:2)
      IBEG = IBEG + 1
      CHAR8_OUT(IBEG:IBEG) = '.'

      CHAR8_OUT(8:8) = TEMP_CHAR(11:11)
      IEND = 7
      IF (TEMP_CHAR(10:10) /= '0') THEN
         CHAR8_OUT(IEND:IEND) = TEMP_CHAR(10:10)
         IEND = IEND - 1
      ENDIF
      CHAR8_OUT(IEND:IEND) = TEMP_CHAR(9:9)
      IEND=IEND-1

      DO I=4,11
         IBEG = IBEG+1
         IF ((CHAR8_OUT(IBEG:IBEG) == '-') .OR. (CHAR8_OUT(IBEG:IBEG) == '+')) exit
         CHAR8_OUT(IBEG:IBEG) = TEMP_CHAR(I:I)
      ENDDO

! **********************************************************************************************************************************

      END SUBROUTINE REAL_DATA_TO_C8FLD



      FUNCTION TO_UPPER(IN_TEXT) RESULT(rslt)

      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT (IN)   :: IN_TEXT
      CHARACTER(len=len(IN_TEXT))     :: rslt
      INTEGER                         :: I, J
      CHARACTER(26), parameter        :: UPP = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      CHARACTER(26), parameter        :: LOW = 'abcdefghijklmnopqrstuvwxyz'

      DO I = 1,len(IN_TEXT)
         J = index(LOW, IN_TEXT(I:I))
         IF (J>0) THEN
            rslt(I:I) = UPP(J:J)
         ELSE
            rslt(I:I) = IN_TEXT(I:I)
         ENDIF
      ENDDO

      END FUNCTION TO_UPPER

   END MODULE TEXT_FIELD_UTILS
