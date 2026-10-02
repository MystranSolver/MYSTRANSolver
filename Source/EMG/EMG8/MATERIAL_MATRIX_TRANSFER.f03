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

   MODULE MATERIAL_MATRIX_TRANSFER

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MATGET, MATPUT

   CONTAINS

      SUBROUTINE MATGET ( A, NROWA, NCOLA, BEG_ROW, BEG_COL, NROW, NCOL, B )

! Gets a NROW x NCOL partition of a matrix starting at row BEG_ROW and column BEG_COL

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATGET'

      INTEGER(LONG), INTENT(IN)       :: BEG_COL           ! Beginning row of input matrix to get partition from
      INTEGER(LONG), INTENT(IN)       :: BEG_ROW           ! Beginning col of input matrix to get partition from
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG), INTENT(IN)       :: NCOLA             ! Number of cols in input matrix
      INTEGER(LONG), INTENT(IN)       :: NROWA             ! Number of rows in input matrix
      INTEGER(LONG), INTENT(IN)       :: NCOL              ! No. of cols to get from input matrix
      INTEGER(LONG), INTENT(IN)       :: NROW              ! No. of rows to get from input matrix
      INTEGER(LONG)                   :: ICNT              ! A computed index into array A
      INTEGER(LONG)                   :: ICNT0             ! Part of ICNT
      INTEGER(LONG)                   :: II                ! Counter


      REAL(DOUBLE) , INTENT(IN)       :: A(NROWA*NCOLA)    ! Input matrix from which a partition will be extracted
      REAL(DOUBLE) , INTENT(OUT)      :: B(NROW*NCOL)      ! Output matrix, which is the partition extracted from A



! **********************************************************************************************************************************
! Initialize outputs

      B = ZERO

      ICNT0 = NROWA*(BEG_COL - 2) + BEG_ROW - 1
      II = 0
      DO J=1,NCOL
         ICNT0 = ICNT0 + NROWA
         DO I=1,NROW
            ICNT = ICNT0 + I
            II = II + 1
            B(II) = A(ICNT)
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MATGET


      SUBROUTINE MATPUT ( B, NROWA, NCOLA, BEG_ROW, BEG_COL, NROW, NCOL, A )

! Puts a NROW x NCOL partition of a matrix into another matrix starting at row BEG_ROW and column BEG_COL

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATPUT'

      INTEGER(LONG), INTENT(IN)       :: BEG_COL           ! Beginning row of input matrix to get partition from
      INTEGER(LONG), INTENT(IN)       :: BEG_ROW           ! Beginning col of input matrix to get partition from
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG), INTENT(IN)       :: NCOLA             ! Number of cols in input matrix
      INTEGER(LONG), INTENT(IN)       :: NROWA             ! Number of rows in input matrix
      INTEGER(LONG), INTENT(IN)       :: NCOL              ! No. of cols to get from input matrix
      INTEGER(LONG), INTENT(IN)       :: NROW              ! No. of rows to get from input matrix
      INTEGER(LONG)                   :: ICNT              ! A computed index into array A
      INTEGER(LONG)                   :: ICNT0             ! Part of ICNT
      INTEGER(LONG)                   :: II                ! Counter


      REAL(DOUBLE) , INTENT(IN)       :: B(NROW*NCOL)      ! Input matrix that will be put into A
      REAL(DOUBLE) , INTENT(INOUT)    :: A(NROWA*NCOLA)    ! Output matrix, containing inserted terms from B



! **********************************************************************************************************************************
      ICNT0 = NROWA*(BEG_COL-2) + BEG_ROW - 1
      II = 0
      DO J=1,NCOL
         ICNT0 = ICNT0 + NROWA
         DO I=1,NROW
            ICNT = ICNT0 + I
            II = II + 1
            A(ICNT) = B(II)
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MATPUT

   END MODULE MATERIAL_MATRIX_TRANSFER
