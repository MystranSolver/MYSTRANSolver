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

   MODULE DOF_LOOKUP_UTILS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: GET_GRID_AND_COMP, TDOF_COL_NUM

   CONTAINS

      SUBROUTINE TDOF_COL_NUM ( CHAR_SET, COL_NUM )

! Converts character representation of displ set (G, N, F, etc) to a column number in array TDOF

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MTDOF, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TDOF_COL_NUM'
      CHARACTER(LEN=*), INTENT(IN)    :: CHAR_SET          ! The char description of the displ set that was input (e.g. 'G ', 'SB')

      INTEGER(LONG), INTENT(OUT)      :: COL_NUM           ! Col number in array TDOF where displ set CHAR_SET exists
      INTEGER(LONG), PARAMETER        :: OFFSET    = 4     ! Columns of TDOF prior to where the G-set begins



! **********************************************************************************************************************************
! Initialize outputs

      COL_NUM = 1

! Get output column number (COL_NUM) in TDOF for the requested set

      IF      (CHAR_SET == 'G ') THEN   ;   COL_NUM = OFFSET +  1
      ELSE IF (CHAR_SET == 'M ') THEN   ;   COL_NUM = OFFSET +  2
      ELSE IF (CHAR_SET == 'N ') THEN   ;   COL_NUM = OFFSET +  3
      ELSE IF (CHAR_SET == 'SA') THEN   ;   COL_NUM = OFFSET +  4
      ELSE IF (CHAR_SET == 'SB') THEN   ;   COL_NUM = OFFSET +  5
      ELSE IF (CHAR_SET == 'SG') THEN   ;   COL_NUM = OFFSET +  6
      ELSE IF (CHAR_SET == 'SZ') THEN   ;   COL_NUM = OFFSET +  7
      ELSE IF (CHAR_SET == 'SE') THEN   ;   COL_NUM = OFFSET +  8
      ELSE IF (CHAR_SET == 'S ') THEN   ;   COL_NUM = OFFSET +  9
      ELSE IF (CHAR_SET == 'F ') THEN   ;   COL_NUM = OFFSET + 10
      ELSE IF (CHAR_SET == 'O ') THEN   ;   COL_NUM = OFFSET + 11
      ELSE IF (CHAR_SET == 'A ') THEN   ;   COL_NUM = OFFSET + 12
      ELSE IF (CHAR_SET == 'R ') THEN   ;   COL_NUM = OFFSET + 13
      ELSE IF (CHAR_SET == 'L ') THEN   ;   COL_NUM = OFFSET + 14
      ELSE IF (CHAR_SET == 'U1') THEN   ;   COL_NUM = OFFSET + 15
      ELSE IF (CHAR_SET == 'U2') THEN   ;   COL_NUM = OFFSET + 16
      ELSE                                                 ! Incorrect set designation
         WRITE(ERR,1327) SUBR_NAME,CHAR_SET
         WRITE(F06,1327) SUBR_NAME,CHAR_SET
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
      ENDIF

! Make sure that we have not coded for a column number greater than MTDOF (which was the number of columns
! allocated to array TDOF)

      IF (COL_NUM > MTDOF) THEN
         WRITE(ERR,1328) SUBR_NAME,COL_NUM,MTDOF
         WRITE(F06,1328) SUBR_NAME,COL_NUM,MTDOF
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1327 FORMAT(' *ERROR  1327: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INPUT VARIABLE CHAR_SET = "',A2,'" IS NOT ONE OF THE CORRECT DESIGNATIONS FOR A DISPL SET')

 1328 FORMAT(' *ERROR  1328: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' OUTPUT VARIABLE COL_NUM = ',I3,' CANNOT BE GREATER THAN MTDOF = ',I3)

! **********************************************************************************************************************************

      END SUBROUTINE TDOF_COL_NUM


      SUBROUTINE GET_GRID_AND_COMP ( X_SET, DOF_NUM, GRIDV, COMPV )

! Gets the grid and displ component (1-6) numbers for a DOF listed in array TDOFI.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TDOFI

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_GRID_AND_COMP'
      CHARACTER(LEN=*), INTENT(IN)    :: X_SET             ! Displ set designator (if one exists) for the col in TDOFI

      INTEGER(LONG), INTENT(IN)       :: DOF_NUM           ! DOF number in TDOF
      INTEGER(LONG), INTENT(OUT)      :: COMPV             ! Comp. num corresponding to DOF_NUM in array TDOFI, col X_SET_COL_NUM
      INTEGER(LONG), INTENT(OUT)      :: GRIDV             ! Grid num corresponding to DOF_NUM in array TDOFI, col X_SET_COL_NUM
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: X_SET_COL_NUM     ! Col number, in TDOFI array, of the X-set DOF list



! **********************************************************************************************************************************
! Initialize outputs

      GRIDV = 0
      COMPV = 0

! Calc outputs

      CALL TDOF_COL_NUM ( X_SET, X_SET_COL_NUM )
      DO I=1,NDOFG
         IF (TDOFI(I, X_SET_COL_NUM) == DOF_NUM) THEN
            GRIDV = TDOFI(I,1)
            COMPV = TDOFI(I,2)
            EXIT
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

! **********************************************************************************************************************************

      END SUBROUTINE GET_GRID_AND_COMP

   END MODULE DOF_LOOKUP_UTILS
