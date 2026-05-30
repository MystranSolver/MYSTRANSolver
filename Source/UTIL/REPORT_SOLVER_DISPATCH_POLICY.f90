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

      SUBROUTINE REPORT_SOLVER_DISPATCH_POLICY ( MATRIX_NAME, CALLING_SUBR )

! --- BANDED_optimizisation -begin-- !
! Reports the effective solver dispatch policy. This routine is diagnostic only; it does not change solver selection.
! --- BANDED_optimizisation -end-- !

      USE PENTIUM_II_KIND, ONLY       :  BYTE
      USE IOUNT1, ONLY                :  ERR, F06
      USE PARAMS, ONLY                :  SOLLIB, SPARSE_FLAVOR, SUPINFO

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR
      CHARACTER(LEN=*), INTENT(IN)    :: MATRIX_NAME

! **********************************************************************************************************************************

      IF (SOLLIB == 'BANDED  ') THEN
         WRITE(ERR,9001) MATRIX_NAME, CALLING_SUBR
         IF (SUPINFO == 'N') WRITE(F06,9001) MATRIX_NAME, CALLING_SUBR
      ELSE IF (SOLLIB == 'SPARSE  ') THEN
         WRITE(ERR,9002) MATRIX_NAME, CALLING_SUBR, SPARSE_FLAVOR
         IF (SUPINFO == 'N') WRITE(F06,9002) MATRIX_NAME, CALLING_SUBR, SPARSE_FLAVOR
      ELSE
         WRITE(ERR,9003) MATRIX_NAME, CALLING_SUBR, SOLLIB
         IF (SUPINFO == 'N') WRITE(F06,9003) MATRIX_NAME, CALLING_SUBR, SOLLIB
      ENDIF

      RETURN

! **********************************************************************************************************************************
 9001 FORMAT(' *INFORMATION: SOLVER DISPATCH POLICY FOR MATRIX ',A,' IN ',A,': BANDED_ONLY',                                      &
                    /,14X,' SOLLIB=BANDED WILL USE THE BANDED PATH; NO AUTOMATIC SUPERLU FALLBACK IS ENABLED.')

 9002 FORMAT(' *INFORMATION: SOLVER DISPATCH POLICY FOR MATRIX ',A,' IN ',A,': SPARSE',                                           &
                    /,14X,' SOLLIB=SPARSE WILL USE SPARSE_FLAVOR=',A8,'.')

 9003 FORMAT(' *WARNING    : SOLVER DISPATCH POLICY FOR MATRIX ',A,' IN ',A,' HAS UNKNOWN SOLLIB=',A8)

! **********************************************************************************************************************************

      END SUBROUTINE REPORT_SOLVER_DISPATCH_POLICY
