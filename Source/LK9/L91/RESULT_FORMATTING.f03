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

   MODULE RESULT_FORMATTING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: WRT_REAL_TO_CHAR_VAR

   CONTAINS

      SUBROUTINE WRT_REAL_TO_CHAR_VAR ( REAL_VAR, NROWS, NCOLS, ROW_NUM, CHAR_VAR )

! Writes real values in one row of array REAL_VAR to a character variable, CHAR_VAR. This is done so that we can change real
! zero's to character values. Then the calling routine that will write REAL_VAR values to the F06 file can have values such as
! '  0.0         ' written rather than the 1ES14.6 format: 0.000000E+00. This subr is written for word lengths of 14 bytes only

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR
      USE CONSTANTS_1, ONLY           :  ZERO

      USE TEXT_FIELD_UTILS, ONLY      :  FMT_ES14_6

      IMPLICIT NONE

      CHARACTER(14*BYTE), INTENT(OUT) :: CHAR_VAR(NCOLS)      ! Character representation of the real data in one row of REAL_VAR

      INTEGER(LONG), INTENT(IN)       :: NCOLS                ! Number of cols in array REAL_VAR
      INTEGER(LONG), INTENT(IN)       :: NROWS                ! Number of rows in array REAL_VAR
      INTEGER(LONG), INTENT(IN)       :: ROW_NUM              ! Row number in array REAL_VAR to write
      INTEGER(LONG)                   :: J                    ! DO loop index

      REAL(DOUBLE) , INTENT(IN)       :: REAL_VAR(NROWS,NCOLS)!

! **********************************************************************************************************************************
! FMT_ES14_6 produces the same 14-char output as Fortran's 1ES14.6 edit descriptor but avoids the
! per-value internal WRITE overhead. Exact zeros are then overwritten with this routine's historical
! '  0.0         ' substitution so the F06 output stays bit-identical.
      DO J=1,NCOLS
         IF (ABS(REAL_VAR(ROW_NUM,J)) == ZERO) THEN
            CHAR_VAR(J) = '  0.0         '
         ELSE
            CALL FMT_ES14_6 ( REAL_VAR(ROW_NUM,J), CHAR_VAR(J) )
         ENDIF
      ENDDO

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE WRT_REAL_TO_CHAR_VAR

   END MODULE RESULT_FORMATTING
