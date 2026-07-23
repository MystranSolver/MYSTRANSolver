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

      SUBROUTINE CALC_ELEM_STRESSES ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

! Calls routines that process the STRESS array calculated in subr ELEM_STRE_STRN_ARRAYS to obtain element specific stress values
! that will be written to the F06 file

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  TYPE
      USE CONSTANTS_1, ONLY           :  ZERO

      USE CALC_ELEM_STRESSES_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_ELEM_STRESSES'
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_OGEL        ! If 'Y' then write data to array OGEL
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_FEMAP       ! If 'Y' then write data to array FEMAP arrays

      INTEGER(LONG), INTENT(IN)       :: SIZE_ALLOCATED    ! No. of rows allocated to array that will be written to in a subr
!                                                            called here so we can check that we don't try to write more rows
      INTEGER(LONG), INTENT(IN)       :: NUM_FEMAP_ROWS    ! Number of rows that will be written to FEMAP arrays
      INTEGER(LONG), INTENT(INOUT)    :: NUM1              ! Cum rows written to OGEL prior to running this subr




! **********************************************************************************************************************************
! Calculate STRESS(4-6) for elements that have nonzero STRESS(4-6)

      IF      ((TYPE(1:3) == 'BAR'  ) .OR. (TYPE(1:4) == 'BEAM') .OR. (TYPE(1:4) == 'BUSH') .OR. (TYPE(1:4) == 'ELAS') .OR.        &
               (TYPE(1:3) == 'ROD'  )) THEN
         CALL ONE_D_STRESS_OUTPUTS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

      ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'QUAD8') .OR.                                 &
               (TYPE(1:5) == 'SHEAR') .OR. (TYPE(1:6) == 'USERIN')) THEN
         CALL SHELL_STRESS_OUTPUTS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

      ELSE IF ((TYPE(1:4) == 'HEXA' ) .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
         CALL SOLID_STRESS_OUTPUTS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9203) SUBR_NAME, TYPE
         WRITE(F06,9203) SUBR_NAME, TYPE
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
! not used
 9200 FORMAT(' *ERROR  9200: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ARRAY OGEL WAS ALLOCATED TO HAVE ',I12,' ROWS. ATTEMPT TO WRITE TO OGEL BEYOND THIS')

 9203 FORMAT(' *ERROR  9203: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INCORRECT ELEMENT TYPE = "',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE CALC_ELEM_STRESSES
