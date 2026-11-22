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

! Checks sensibility of the 3 MOI's of a BAR or BEAM element. Zero I1 and/or I2 are allowed (no bending stiffness in that plane)
! and are left as zero. The BAR element handles them when it divides by I1*I2 - I12^2.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CHECK_BAR_MOIs'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME              ! Either PBAR, PBARL or PBEAM
      CHARACTER(LEN=*), INTENT(IN)    :: ID                ! Character value of the bar's ID

      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator


      REAL(DOUBLE), INTENT(IN)        :: I1                ! MOI of the bar or beam
      REAL(DOUBLE), INTENT(IN)        :: I2                ! MOI of the bar or beam
      REAL(DOUBLE), INTENT(IN)        :: I12               ! MOI of the bar or beam


! **********************************************************************************************************************************
! Initialize

      IERR = 0


! If I12 is not zero, make sure I1*I2 > I12^2. I12^2 <= I1*I2 is proven by the Cauchy-Schwarz inequality, and equality is not
! allowed because the BAR element divides by I1*I2 - I12^2 when I12 is not zero. So I1*I2 - I12^2 can only be zero if I12 = 0.

      IF (I12 /= ZERO) THEN

         IF (I12*I12 >= I1*I2) THEN
            IERR = IERR + 1
            WRITE(ERR,1195) NAME, ID, I1, I2, I12
            WRITE(F06,1195) NAME, ID, I1, I2, I12
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1195 FORMAT(' *ERROR  1195: THE MOMENTS AND PRODUCTS OF INERTIA ON ',A,'   ',A8,' DO NOT SATISIFY THE REQUIREMENT THAT:'          &
                    ,/,14X,' I12^2  <  I1*I2  WHERE: I1 = ',1ES10.3,', I2 = ',1ES10.3,', I12 = ',1ES10.3)

! **********************************************************************************************************************************

      END SUBROUTINE CHECK_BAR_MOIs


   END MODULE BAR_PROPERTY_VALIDATION
