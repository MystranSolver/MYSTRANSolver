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

   MODULE BAR_PROPERTY_VALIDATION

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CHECK_BAR_MOIS

   CONTAINS

      SUBROUTINE CHECK_BAR_MOIs ( NAME, ID, I1, I2, I12, IERR )

! Checks sensibility of the 3 MOI's of a BAR or BEAM element and replaces zero values with small finite ones

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  EPSIL, SUPINFO
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CHECK_BAR_MOIs'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Either PBAR, PBARL or PBEAM
      CHARACTER(LEN=*), INTENT(IN)    :: ID                ! Character value of the bar's ID

      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator


      REAL(DOUBLE), INTENT(INOUT)     :: I1                ! MOI of the bar or beam
      REAL(DOUBLE), INTENT(INOUT)     :: I2                ! MOI of the bar or beam
      REAL(DOUBLE), INTENT(INOUT)     :: I12               ! MOI of the bar or beam
      REAL(DOUBLE)                    :: EPS1              ! A small number


! **********************************************************************************************************************************
! Initialize

      IERR = 0

      EPS1 = EPSIL(1)


! If I12 is zero replace a zero value of I1 and/or I2 with a small positive value. If I12 is not zero, make sure I1*I2 > I12^2

      IF (DABS(I12) <= EPS1) THEN
         IF (DABS(I1) < EPS1) THEN
            I1 = 10*EPS1
            WRITE(ERR,1001) NAME, ID, 'I1', I1
            IF (SUPINFO == 'N') THEN
               WRITE(F06,1001) NAME, ID, 'I1', I1
            ENDIF
         ENDIF

         IF (DABS(I2) <= EPS1) THEN
            I2 = 10*EPS1
            WRITE(ERR,1001) NAME, ID, 'I2', I2
            IF (SUPINFO == 'N') THEN
               WRITE(F06,1001) NAME, ID, 'I2', I2
            ENDIF
         ENDIF

      ELSE                                                 ! I12^2 must be <= I1*I2 is proven by the Cauchy-Schwarz inequality

         IF (I12*I12 > I1*I2) THEN
            IERR = IERR + 1
            WRITE(ERR,1195) NAME, ID, I1, I2, I12
            WRITE(F06,1195) NAME, ID, I1, I2, I12
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1001 FORMAT(' *INFORMATION: FOR ',A,'   ',A8,' MOMENT OF INERTIA ',A,' HAS BEEN CHANGED FROM 0 TO A SMALL NUMBER = ',1ES10.3)

 1195 FORMAT(' *ERROR  1195: THE MOMENTS AND PRODUCTS OF INERTIA ON ',A,'   ',A8,' DO NOT SATISIFY THE REQUIREMENT THAT:'          &
                    ,/,14X,' I12^2  <=  I1*12  WHERE: I1 = ',1ES10.3,', I2 = ',1ES10.3,', I12 = ',1ES10.3)

! **********************************************************************************************************************************

      END SUBROUTINE CHECK_BAR_MOIs


   END MODULE BAR_PROPERTY_VALIDATION
