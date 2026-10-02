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

   MODULE QUADRATURE

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ORDER_GAUSS, ORDER_TRIA, ORDER_TETRA

   CONTAINS

      SUBROUTINE ORDER_GAUSS ( KORDER, SSS, HHH )

! Calculates abscissa and weight coefficients for Gaussian integration of order KORDER = 1 to 10.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, MEFE
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO
      USE CONSTANTS_GAUSS, ONLY       :  HHV, SSV
      USE MODEL_STUF, ONLY            :  EMG_IFE, ERR_SUB_NAM, NUM_EMG_FATAL_ERRS

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ORDER_GAUSS'

      INTEGER(LONG), INTENT(IN)       :: KORDER               ! Gaussian integration order to use
      INTEGER(LONG)                   :: I                    ! DO loop index
      INTEGER(LONG)                   :: II                   ! A term in a computed index into SSS, HHH arrays
      INTEGER(LONG)                   :: JJ                   ! A computed index into SSS, HHH arrays
      INTEGER(LONG)                   :: KK                   ! A computed index into SSS, HHH arrays
      INTEGER(LONG)                   :: LL                   ! A computed index into SSS, HHH arrays
      INTEGER(LONG)                   :: MM                   ! A computed index into SSS, HHH arrays
      INTEGER(LONG)                   :: NN                   ! A term in a computed index into SSS, HHH arrays
      INTEGER(LONG)                   :: IBEGIN(11) = (/0, 1, 2, 4, 6, 9,12,16,20,25,30/)


      REAL(DOUBLE) ,INTENT(OUT)       :: SSS(MAX_ORDER_GAUSS) ! Gauss abscissa's
      REAL(DOUBLE) ,INTENT(OUT)       :: HHH(MAX_ORDER_GAUSS) ! Gauss weight coeffs

      INTRINSIC MOD



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,MAX_ORDER_GAUSS
         SSS(I) = ZERO
         HHH(I) = ZERO
      ENDDO

! Check KORDER to make sure it no less than 1 nor greater than MAX_ORDER.  Error if not.

      IF ((KORDER >= 1) .AND. (KORDER <= MAX_ORDER_GAUSS)) THEN

! Abscissa and weight coefficients for Gaussian integ. of order KORDER

         IF (KORDER == 1) THEN
            SSS(1) = ZERO
            HHH(1) = TWO
         ELSE
            II = IBEGIN(KORDER)
            LL = IBEGIN(KORDER+1)
            JJ = LL - II
            NN = 2*JJ
            IF (MOD(KORDER,2) /= 0) THEN
               SSS(JJ) = SSV(LL-1)
               HHH(JJ) = HHV(LL-1)
               NN = NN - 1
               JJ = JJ - 1
            ENDIF
            DO KK=1,JJ
               LL = II + KK - 1
               MM = NN - KK + 1
               SSS(MM) =  SSV(LL)
               HHH(MM) =  HHV(LL)
               SSS(KK) = -SSV(LL)
               HHH(KK) =  HHV(LL)
             ENDDO
         ENDIF

      ELSE

         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1931) SUBR_NAME, MAX_ORDER_GAUSS, KORDER
            WRITE(F06,1931) SUBR_NAME, MAX_ORDER_GAUSS, KORDER
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1931
               EMG_IFE(NUM_EMG_FATAL_ERRS,2)   = KORDER
            ENDIF
         ENDIF
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit

      ENDIF


      RETURN

! **********************************************************************************************************************************
 1931 FORMAT(' *ERROR  1931: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' GAUSSIAN INTEGRATION ORDER CANNOT BE < 1 OR >',I3,' BUT VALUE IS ',I8)


! **********************************************************************************************************************************

      END SUBROUTINE ORDER_GAUSS


      SUBROUTINE ORDER_TRIA ( KORDER, SS_I, SS_J, HH_IJ )

! Calculates abscissa and weight coefficients for triangular integration for the PENTA element

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_TRIA
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, SIXTH, THIRD, HALF, TWO

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ORDER_TRIA'

      INTEGER(LONG), INTENT(IN)       :: KORDER                ! Triangular integration order to use
      INTEGER(LONG)                   :: I                     ! DO loop index


      REAL(DOUBLE) ,INTENT(OUT)       :: SS_I(MAX_ORDER_TRIA)  ! Triangular integration abscissa's
      REAL(DOUBLE) ,INTENT(OUT)       :: SS_J(MAX_ORDER_TRIA)  ! Triangular integration abscissa's
      REAL(DOUBLE) ,INTENT(OUT)       :: HH_IJ(MAX_ORDER_TRIA) ! Triangular integration weight coeffs
      REAL(DOUBLE) , PARAMETER        :: A1 = .0597158717D0    ! Intermediate constant
      REAL(DOUBLE) , PARAMETER        :: A2 = .7974269853D0    ! Intermediate constant
      REAL(DOUBLE) , PARAMETER        :: B1 = .4701420641D0    ! Intermediate constant
      REAL(DOUBLE) , PARAMETER        :: B2 = .1012865073D0    ! Intermediate constant
      REAL(DOUBLE) , PARAMETER        :: W1 = .1125D0          ! Intermediate constant
      REAL(DOUBLE) , PARAMETER        :: W2 = .0661970763D0    ! Intermediate constant
      REAL(DOUBLE) , PARAMETER        :: W3 = .0629695902D0    ! Intermediate constant



! **********************************************************************************************************************************
      DO I=1,MAX_ORDER_TRIA
         SS_I(I)  = ZERO
         SS_J(I)  = ZERO
         HH_IJ(I) = ZERO
      ENDDO

      IF      (KORDER == 1) THEN

         SS_I(1) = THIRD    ;     SS_J(1) = THIRD    ;     HH_IJ(1) = HALF

      ELSE IF (KORDER == 3) THEN

         SS_I(1) = SIXTH    ;     SS_J(1) = TWO*THIRD;     HH_IJ(1) = SIXTH
         SS_I(2) = SIXTH    ;     SS_J(2) = SIXTH    ;     HH_IJ(2) = SIXTH
         SS_I(3) = TWO*THIRD;     SS_J(3) = SIXTH    ;     HH_IJ(3) = SIXTH

      ELSE IF (KORDER == 7) THEN

         SS_I(1) = THIRD;     SS_J(1) = THIRD;     HH_IJ(1) = W1
         SS_I(2) = A1   ;     SS_J(2) = B1   ;     HH_IJ(2) = W2
         SS_I(3) = B1   ;     SS_J(3) = A1   ;     HH_IJ(3) = W2
         SS_I(4) = B1   ;     SS_J(4) = B1   ;     HH_IJ(4) = W2
         SS_I(5) = A2   ;     SS_J(5) = B2   ;     HH_IJ(5) = W3
         SS_I(6) = B2   ;     SS_J(6) = A2   ;     HH_IJ(6) = W3
         SS_I(7) = B2   ;     SS_J(7) = B2   ;     HH_IJ(7) = W3

      ELSE

         WRITE(ERR,1931) SUBR_NAME, KORDER
         WRITE(F06,1931) SUBR_NAME, KORDER
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1931 FORMAT(' *ERROR  1931: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TRIANGULAR INTEGRATION ORDER MUST BE EITHER 1, 3 OR 7 BUT VALUE IS ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE ORDER_TRIA


      SUBROUTINE ORDER_TETRA ( KORDER, SSS_I, SSS_J, SSS_K, HHH_IJK )

! Calculates abscissa and weight coefficients for triangular integration for the TETRA element

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_TETRA
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, SIXTH, QUARTER, HALF, ONE, TWO, TWELVE

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ORDER_TETRA'

      INTEGER(LONG), INTENT(IN)       :: KORDER                   ! Triangular integration order to use
      INTEGER(LONG)                   :: I                        ! DO loop index


      REAL(DOUBLE) , INTENT(OUT)      :: SSS_I (MAX_ORDER_TETRA)  ! Gauss abscissa's
      REAL(DOUBLE) , INTENT(OUT)      :: SSS_J (MAX_ORDER_TETRA)  ! Gauss abscissa's
      REAL(DOUBLE) , INTENT(OUT)      :: SSS_K (MAX_ORDER_TETRA)  ! Gauss abscissa's
      REAL(DOUBLE) , INTENT(OUT)      :: HHH_IJK(MAX_ORDER_TETRA) ! Gauss weight coeffs
      REAL(DOUBLE) , PARAMETER        :: ALPHA = .58541020D0      ! Intermediate constant
      REAL(DOUBLE) , PARAMETER        :: BETA  = .13819660D0      ! Intermediate constant



! **********************************************************************************************************************************
      DO I=1,MAX_ORDER_TETRA
         SSS_I(I)   = ZERO
         SSS_J(I)   = ZERO
         SSS_K(I)   = ZERO
         HHH_IJK(I) = ZERO
      ENDDO

      IF      (KORDER == 1) THEN

         SSS_I(1) = QUARTER  ;     SSS_J(1) = QUARTER  ;     SSS_K(1) = QUARTER  ;     HHH_IJK(1) = SIXTH

      ELSE IF (KORDER == 4) THEN

         SSS_I(1) = ALPHA    ;     SSS_J(1) = BETA     ;     SSS_K(1) = BETA     ;     HHH_IJK(1) = ONE/(TWO*TWELVE)
         SSS_I(2) = BETA     ;     SSS_J(2) = ALPHA    ;     SSS_K(2) = BETA     ;     HHH_IJK(2) = ONE/(TWO*TWELVE)
         SSS_I(3) = BETA     ;     SSS_J(3) = BETA     ;     SSS_K(3) = ALPHA    ;     HHH_IJK(3) = ONE/(TWO*TWELVE)
         SSS_I(4) = BETA     ;     SSS_J(4) = BETA     ;     SSS_K(4) = BETA     ;     HHH_IJK(4) = ONE/(TWO*TWELVE)

      ELSE

         WRITE(ERR,1931) SUBR_NAME, KORDER
         WRITE(F06,1931) SUBR_NAME, KORDER
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1931 FORMAT(' *ERROR  1931: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TETRAHEDRAL INTEGRATION ORDER MUST BE EITHER 1 OR 4 BUT VALUE IS ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE ORDER_TETRA

   END MODULE QUADRATURE
