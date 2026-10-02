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

   MODULE MATRIX_TEXT_OUTPUT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: WRITE_MATRIX_BY_ROWS, WRITE_GRID_COORDS, WRITE_INTEGER_VEC, WRITE_VECTOR

   CONTAINS

      SUBROUTINE WRITE_MATRIX_BY_COLS ( MAT_DESCR, MATOUT, NROWS, NCOLS, OUT_UNT )

! Writes a matrix one column at a time in a format that has 10 terms across the page (repeated un til col is completely written)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_MATRIX_BY_COLS'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_DESCR          ! Character descriptor of the matrix to be printed
      CHARACTER(131*BYTE)             :: HEADER             ! MAT_DESCRIPTOR centered in a line of output

      INTEGER(LONG), INTENT(IN)       :: NROWS              ! Number of rows in matrix MATOUT
      INTEGER(LONG), INTENT(IN)       :: NCOLS              ! Number of cols in matrix MATOUT
      INTEGER(LONG), INTENT(IN)       :: OUT_UNT            ! Output unit number
      INTEGER(LONG)                   :: I,J,K              ! DO loop indices
      INTEGER(LONG)                   :: NUM_LEFT           !
      INTEGER(LONG)                   :: PAD                ! Number of spaces to pad in HEADER to center MAT_DESCRIPTOR


      REAL(DOUBLE) , INTENT(IN)       :: MATOUT(NROWS,NCOLS)! Matrix to write out
      REAL(DOUBLE)                    :: MAT_LINE(10)       ! Up to 10 terms from one col of MATOUT

      INTRINSIC                       :: LEN



! **********************************************************************************************************************************
      PAD = (132 - LEN(MAT_DESCR))/2
      HEADER(1:) = ' '
      HEADER(PAD+1:) = MAT_DESCR
      WRITE(OUT_UNT,101) HEADER

      DO J=1,NCOLS
         WRITE(OUT_UNT,102) J
         NUM_LEFT = NROWS
         DO I=1,NROWS,10
            IF (NUM_LEFT >= 10) THEN
               DO K=1,10
                  MAT_LINE(K) = MATOUT(I+K-1,J)
               ENDDO
               WRITE(OUT_UNT,103) (MAT_LINE(K),K=1,10)
            ELSE
               DO K=1,NUM_LEFT
                  MAT_LINE(K) = MATOUT(I+K-1,J)
               ENDDO
               WRITE(OUT_UNT,103) (MAT_LINE(K),K=1,NUM_LEFT)
            ENDIF
            NUM_LEFT = NUM_LEFT - 10
         ENDDO
         WRITE(OUT_UNT,*)
      ENDDO



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(1X,/,1X,131A,/)

  102 FORMAT(1X,'COLUMN ',I8,/,1X,'---------------')

  103 FORMAT(1X,10(1ES14.6))

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_MATRIX_BY_COLS


      SUBROUTINE WRITE_MATRIX_BY_ROWS ( MAT_DESCR, MATOUT, NROWS, NCOLS, OUT_UNT )

! Writes a matrix one row at a time in a format that has 10 terms across the page (repeated until row is completely written)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_MATRIX_BY_ROWS'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_DESCR          ! Character descriptor of the matrix to be printed
      CHARACTER(131*BYTE)             :: HEADER             ! MAT_DESCRIPTOR centered in a line of output

      INTEGER(LONG), INTENT(IN)       :: NROWS              ! Number of rows in matrix MATOUT
      INTEGER(LONG), INTENT(IN)       :: NCOLS              ! Number of cols in matrix MATOUT
      INTEGER(LONG), INTENT(IN)       :: OUT_UNT            ! Output unit number
      INTEGER(LONG)                   :: I,J,K              ! DO loop indices
      INTEGER(LONG)                   :: NUM_LEFT           !
      INTEGER(LONG)                   :: PAD                ! Number of spaces to pad in HEADER to center MAT_DESCRIPTOR


      REAL(DOUBLE) , INTENT(IN)       :: MATOUT(NROWS,NCOLS)! Matrix to write out
      REAL(DOUBLE)                    :: MAT_LINE(10)       ! Up to 10 terms from one col of MATOUT

      INTRINSIC                       :: LEN



! **********************************************************************************************************************************
      PAD = (132 - LEN(MAT_DESCR))/2
      HEADER(1:) = ' '
      HEADER(PAD+1:) = MAT_DESCR
      WRITE(OUT_UNT,101) HEADER

      DO I=1,NROWS
         WRITE(OUT_UNT,102) I
         NUM_LEFT = NCOLS
         DO J=1,NCOLS,10
            IF (NUM_LEFT >= 10) THEN
               DO K=1,10
                  MAT_LINE(K) = MATOUT(I,J+K-1)
               ENDDO
               WRITE(OUT_UNT,103) (MAT_LINE(K),K=1,10)
            ELSE
               DO K=1,NUM_LEFT
                  MAT_LINE(K) = MATOUT(I,J+K-1)
               ENDDO
               WRITE(OUT_UNT,103) (MAT_LINE(K),K=1,NUM_LEFT)
            ENDIF
            NUM_LEFT = NUM_LEFT - 10
         ENDDO
         WRITE(OUT_UNT,*)
      ENDDO



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(1X,/,1X,131A,/)

  102 FORMAT(1X,'ROW ',I8,/,1X,'------------')

  103 FORMAT(1X,10(1ES14.6))

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_MATRIX_BY_ROWS


      SUBROUTINE WRITE_GRID_COORDS

! Writes grid coordinates in basic coords to the F06 file if user has Bulk Data PARAM PRTBASIC defined

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  NGRID
      USE MODEL_STUF, ONLY            :  GRID, RGRID

      IMPLICIT NONE

      INTEGER(LONG)                   :: I, J              ! DO loop indices

! **********************************************************************************************************************************
      WRITE(F06,1002)
      WRITE(F06,1003)
      DO I = 1,NGRID
         WRITE(F06,1004) GRID(I,1),(RGRID(I,J),J = 1,3)
      ENDDO
      WRITE(F06,*)
      WRITE(F06,*)

! **********************************************************************************************************************************
 1002 FORMAT(//,48X,'GRID POINT COORDINATES IN BASIC COORDINATE SYSTEM')

 1003 FORMAT(/,36X,'GRID ID           X                 Y                 Z')

 1004 FORMAT(35X,I8,1X,1ES17.6,1X,1ES17.6,1X,1ES17.6)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_GRID_COORDS



      SUBROUTINE WRITE_INTEGER_VEC ( ARRAY_DESCR, INT_VEC, NROWS )

! Writes an integer vector to F06 in a format that has 10 terms across the page (repeated until vector is completely written)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  WRT_ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_INTEGER_VEC'
      CHARACTER(LEN=*), INTENT(IN)    :: ARRAY_DESCR       ! Character descriptor of the integer array to be printed
      CHARACTER(131*BYTE)             :: HEADER            ! MAT_DESCRIPTOR centered in a line of output

      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in matrix MATOUT
      INTEGER(LONG), INTENT(IN)       :: INT_VEC(NROWS)    ! Integer vector to write out
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: INT_VEC_LINE(10)  ! 10 members of INT_VEC
      INTEGER(LONG)                   :: NUM_LEFT          ! Count of the number of rows of INT_VEC left to write out
      INTEGER(LONG)                   :: PAD               ! Number of spaces to pad in HEADER to center MAT_DESCRIPTOR


      INTRINSIC                       :: LEN



! **********************************************************************************************************************************
      PAD = (132 - LEN(ARRAY_DESCR))/2
      HEADER(1:) = ' '
      HEADER(PAD+1:) = ARRAY_DESCR
      WRITE(F06,101) HEADER

      NUM_LEFT = NROWS
      DO I=1,NROWS,10
         IF (NUM_LEFT >= 10) THEN
            DO J=1,10
               INT_VEC_LINE(J) = INT_VEC(I+J-1)
            ENDDO
            WRITE(F06,103) (INT_VEC_LINE(J),J=1,10)
         ELSE
            DO J=1,NUM_LEFT
               INT_VEC_LINE(J) = INT_VEC(I+J-1)
            ENDDO
            WRITE(F06,103) (INT_VEC_LINE(J),J=1,NUM_LEFT)
         ENDIF
         NUM_LEFT = NUM_LEFT - 10
      ENDDO
      WRITE(F06,*)



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(1X,/,1X,A,/)

  103 FORMAT(16X,10(I10))

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_INTEGER_VEC


      SUBROUTINE WRITE_VECTOR ( VEC_NAME, WHAT, NUM, UX )

! Writes a vector in full format to the F06 file

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_VECTOR'
      CHARACTER(LEN=*), INTENT(IN)    :: VEC_NAME          ! Name of vector being output
      CHARACTER(LEN=*), INTENT(IN)    :: WHAT              ! Title over output vector (e.g. DISPL, FORCE, etc.)
      CHARACTER(132*BYTE)             :: LINE_OUT          ! Line to print out (to describe matrix) that is centered

      INTEGER(LONG), INTENT(IN)       :: NUM               ! Size of vector UX to write out
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: INDEX             ! Index into character array LINE_OUT
      INTEGER(LONG)                   :: VEC_NAME_LEN      ! Length of char array VEC_NAME. On input, it is the length as defined
!                                                            in the calling subr. In this subr, VEC_NAME is striped of trailing
!                                                            blanks to get only the actual message. On exit VEC_NAME_LEN is the
!                                                            length of the finite message in VEC_NAME (i.e. without trailing blanks)


      REAL(DOUBLE) , INTENT(IN)       :: UX(NUM)           ! Vector to write out



! **********************************************************************************************************************************
! Strip out trailing blanks from VEC_NAME and put remainder centered in array LINE_OUT

      VEC_NAME_LEN = LEN(VEC_NAME)                         ! This is the length of VEC_NAME as input (includes trailing blanks)

      DO  I=VEC_NAME_LEN,1,-1                              ! Calc length of description in VEC_NAME, excluding trailing blanks
         IF (VEC_NAME(I:I) == ' ') THEN
            CYCLE
         ELSE
            VEC_NAME_LEN = I
            EXIT
         ENDIF
      ENDDO

      LINE_OUT(1:) = ' '                                   ! Center VEC_NAME (w/0 trailing blanks) in LINE_OUT
      INDEX = (LEN(LINE_OUT) - VEC_NAME_LEN)/2
      LINE_OUT(INDEX:) = VEC_NAME(1:VEC_NAME_LEN)

      WRITE(F06,2101) LINE_OUT, WHAT

      DO I=1,NUM
         WRITE(F06,2102) I, UX(I)
      ENDDO
      WRITE(F06,*)



      RETURN

! **********************************************************************************************************************************
 2101 FORMAT(A,//,54X,'I            ',A,/,51X,'(DOF)')


 2102 FORMAT(43X,I12,8X,1ES13.6)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_VECTOR

   END MODULE MATRIX_TEXT_OUTPUT
