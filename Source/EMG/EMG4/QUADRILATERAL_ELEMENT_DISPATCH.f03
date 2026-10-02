! #################################################################################################################################
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

   MODULE QUADRILATERAL_ELEMENT_DISPATCH

   USE QUAD_MEMBRANE_SHEAR, ONLY :  QMEM1, QSHEAR
   USE DKQ, ONLY :  QPLT1
   USE MIN4, ONLY :  QPLT2
   USE MIN4T, ONLY :  QPLT3

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: QDEL1

   CONTAINS

      SUBROUTINE QDEL1 ( OPT, INT_ELEM_ID, WRITE_WARN )

! Calculates, or calls subr's to calculate, quadrilateral element matrices:

!  1) ME        = element mass matrix                  , if OPT(1) = 'Y'
!  2) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  3) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  4) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  5) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  6) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, ERR, F06, WRT_BUG, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, MEFE
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, FOUR, TWELVE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL, IORQ1B, IORQ1M, IORQ1S, IORQ2B, QUAD4TYP
      USE MODEL_STUF, ONLY            :  EID, ELDOF, EMG_IFE, EMG_RFE, EMAT, ERR_SUB_NAM, EB, INTL_MID, KE,                        &
                                         MASS_PER_UNIT_AREA, NUM_EMG_FATAL_ERRS, ME, PCOMP_LAM, PCOMP_PROPS, SHELL_B, TYPE, XEL
      USE MODEL_STUF, ONLY            :  BENSUM, SHRSUM, PHI_SQ, PSI_HAT

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE QUADRATURE, ONLY            :  ORDER_GAUSS
      USE JACOBIAN, ONLY               :  JAC2D
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE


      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'QDEL1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)               ! 'Y'/'N' flags for whether to calc certain elem matrices
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID          ! Internal element ID
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN           ! If 'Y" write warning messages, otherwise do not
      CHARACTER( 1*BYTE)              :: RED_INT_SHEAR        ! If 'Y', use Gaussian weighted average of B matrices for shear terms

      INTEGER(LONG)                   :: GAUSS_PT             ! Gauss point number (used for DEBUG output in subr SHP2DQ
      INTEGER(LONG)                   :: I,J,K,L              ! DO loop indices
      INTEGER(LONG)                   :: IORD                 ! Gaussian integration order for QMEM1 portion of element
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Do not change IORD_PCOMP. It must be such that it squared = number of nodes for the QUAD4 (MIN4T)

      INTEGER(LONG), PARAMETER        :: IORD_PCOMP = 2       ! Int order for nonsym layup PCOMP must be 2 (checked in subr
!                                                               SHELL_ABD_MATRICES)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


      REAL(DOUBLE)                    :: AREA                 ! Elem area
      REAL(DOUBLE)                    :: AR                   ! Elem aspect ratio

      REAL(DOUBLE)                    :: BIG_BB(3,ELDOF,IORQ2B*IORQ2B)
                                                              ! Strain-displ matrix for bending for all Gauss points and all DOF's

      REAL(DOUBLE)                    :: BIG_BBI(3,ELDOF)     ! BIG_BB for 1 Gauss point

      REAL(DOUBLE)                    :: BIG_BM(3,ELDOF,IORQ1M*IORQ1M)
                                                              ! Strain-displ matrix for this elem for all Gauss points/all DOF's

      REAL(DOUBLE)                    :: BIG_BMI(3,ELDOF)     ! BIG_BM for 1 Gauss point

      REAL(DOUBLE)                    :: D1(3)                ! Vector from G.P. 1 to G.P. 3 (a diagonal)
      REAL(DOUBLE)                    :: D2(3)                ! Vector from G.P. 2 to G.P. 4 (a diagonal)
      REAL(DOUBLE)                    :: D1M                  ! Mag of D1
      REAL(DOUBLE)                    :: D2M                  ! Mag of D2
      REAL(DOUBLE)                    :: DETJ                 ! An output from subr JAC2D4, called herein. Determinant of JAC
      REAL(DOUBLE)                    :: DUM1(3,ELDOF)        ! Intermediate result in calc SHELL_B effect on KE
      REAL(DOUBLE)                    :: DUM2(ELDOF,ELDOF)    ! Intermediate result in calc SHELL_B effect on KE
      REAL(DOUBLE)                    :: EPS1                 ! A small number to compare to real zero
      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS) ! An output from subr ORDER, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: INTFAC               ! An integration factor (constant multiplier for the Gauss integr)
      REAL(DOUBLE)                    :: JAC(2,2)             ! An output from subr JAC2D4, called herein. 2 x 2 Jacobian matrix.
      REAL(DOUBLE)                    :: JACI(2,2)            ! An output from subr JAC2D4, called herein. 2 x 2 Jacobian inverse.
      REAL(DOUBLE)                    :: M0                   ! An intermediate variable used in calc elem mass, ME
      REAL(DOUBLE)                    :: SSS(MAX_ORDER_GAUSS) ! An output from subr ORDER, called herein. Gauss abscissa's.
      REAL(DOUBLE)                    :: XSD(4)               ! Diffs in x coords of quad sides in local coords
      REAL(DOUBLE)                    :: YSD(4)               ! Diffs in y coords of quad sides in local coords

      INTRINSIC DSQRT



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Initialize

      BENSUM  = ZERO
      SHRSUM  = ZERO
      PSI_HAT = ZERO
      PHI_SQ  = ZERO

! Calculate side diffs

      XSD(1) = XEL(1,1) - XEL(2,1)                         ! x coord diffs (in local elem coords)
      XSD(2) = XEL(2,1) - XEL(3,1)
      XSD(3) = XEL(3,1) - XEL(4,1)
      XSD(4) = XEL(4,1) - XEL(1,1)

      YSD(1) = XEL(1,2) - XEL(2,2)                         ! y coord diffs (in local elem coords)
      YSD(2) = XEL(2,2) - XEL(3,2)
      YSD(3) = XEL(3,2) - XEL(4,2)
      YSD(4) = XEL(4,2) - XEL(1,2)

      IF ((DEBUG(6) > 0) .AND. (WRT_BUG(0) > 0)) THEN
         WRITE(BUG,*) ' Element side differences in x, y coords:'
         WRITE(BUG,*) ' ---------------------------------------'
         WRITE(BUG,98761) XSD(1), YSD(1)
         WRITE(BUG,98762) XSD(2), YSD(2)
         WRITE(BUG,98763) XSD(3), YSD(3)
         WRITE(BUG,98764) XSD(4), YSD(4)
         WRITE(BUG,*)
      ENDIF

! Calculate area by Gaussian integration

      AREA = ZERO
      CALL ORDER_GAUSS ( 2, SSS, HHH )
      DO I=1,2
         DO J=1,2
            CALL JAC2D ( SSS(I), SSS(J), XSD, YSD, 'N', JAC, JACI, DETJ )
            AREA = AREA + HHH(I)*HHH(J)*DETJ
         ENDDO
      ENDDO

! If AREA <= 0, set error and return

      IF (AREA < EPS1) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1925) EID, TYPE, 'AREA', AREA
            WRITE(F06,1925) EID, TYPE, 'AREA', AREA
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1925
               EMG_RFE(NUM_EMG_FATAL_ERRS,1) = AREA
            ENDIF
         ENDIF
         RETURN
      ENDIF

! Calculate aspect ratio. Note that D1M and D2M have been checked to be > 0 in subr ELMGM2

      DO I=1,2
         D1(I) = XEL(3,I) - XEL(1,I)
         D2(I) = XEL(4,I) - XEL(2,I)
      ENDDO
      D1M = DSQRT(D1(1)*D1(1) + D1(2)*D1(2))
      D2M = DSQRT(D2(1)*D2(1) + D2(2)*D2(2))
      AR = D1M/D2M
      IF (AR < ONE) THEN
         AR = ONE/AR
      ENDIF

! **********************************************************************************************************************************
! Generate the mass matrix for this element. For the pure bending element the mass is based only on the non-structural mass.

      IF (OPT(1) == 'Y') THEN

         M0 = MASS_PER_UNIT_AREA*AREA/FOUR
         ME( 1 ,1) = M0
         ME( 2 ,2) = M0
         ME( 3 ,3) = M0
         ME( 7 ,7) = M0
         ME( 8 ,8) = M0
         ME( 9 ,9) = M0
         ME(13,13) = M0
         ME(14,14) = M0
         ME(15,15) = M0
         ME(19,19) = M0
         ME(20,20) = M0
         ME(21,21) = M0

      ENDIF

! **********************************************************************************************************************************
      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(5) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         IF (TYPE == 'SHEAR   ') THEN
            IORD = IORQ1S
            RED_INT_SHEAR = 'N'
            CALL QSHEAR ( OPT, IORD, RED_INT_SHEAR, XSD, YSD )
         ENDIF

         IF (TYPE(1:5) == 'QUAD4') THEN
            IF (INTL_MID(1) /= 0) THEN
               IORD = IORQ1M
               IF (IORQ1S < IORQ1M) THEN
                  RED_INT_SHEAR = 'Y'
               ELSE
                  RED_INT_SHEAR = 'N'
               ENDIF
               CALL QMEM1 ( OPT, INT_ELEM_ID, IORD, RED_INT_SHEAR, AREA, XSD, YSD, BIG_BM )
            ENDIF
         ENDIF

         IF (TYPE(1:6) == 'QUAD4K') THEN

            IF (INTL_MID(2) /= 0) THEN

               if ((pcomp_props == 'Y') .and. (pcomp_lam == 'NON')) then

                  WRITE(ERR,*) ' *ERROR: Code not written for SHELL_B effect on KE yet for QUAD4K elements'
                  WRITE(ERR,*) '         Or, if QUAD4, make sure that the element has nonzero transverse shear moduli, G1Z, G2Z'
                  WRITE(F06,*) ' *ERROR :Code not written for SHELL_B effect on KE yet for QUAD4K elements'
                  WRITE(F06,*) '         Or, if QUAD4, make sure that the element has nonzero transverse shear moduli, G1Z, G2Z'
                  call outa_here ( 'Y' )

               endif

               CALL QPLT1 ( OPT, AREA, XSD, YSD )

            ENDIF

         ELSE IF (TYPE(1:6) == 'QUAD4 ') THEN
            IF (INTL_MID(2) /= 0) THEN                     ! If MID2 = 0, do not calculate bending or transverse shear stiffness

               IF      (QUAD4TYP == 'MIN4  ') THEN
                  CALL QPLT2 ( OPT, AREA, XSD, YSD, BIG_BB )
               ELSE IF (QUAD4TYP == 'MIN4T ') THEN
                  CALL QPLT3 ( OPT, AREA, XSD, YSD, BIG_BB )
               ELSE
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR = FATAL_ERR + 1
                  IF (WRT_ERR > 0) THEN
                     WRITE(ERR,1927) SUBR_NAME, QUAD4TYP
                     WRITE(F06,1927) SUBR_NAME, QUAD4TYP
                  ELSE
                     IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
                        ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
                        EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1927
                     ENDIF
                  ENDIF
                  CALL OUTA_HERE ( 'Y' )
               ENDIF

            ENDIF

         ENDIF

      ENDIF


! **********************************************************************************************************************************
! Calc BM'*SHELL_B*BB (and its transpose) and add to KE. Only do this if this is a composite element with nonsym layup

      IF (OPT(4) == 'Y') THEN

         IF ((PCOMP_PROPS == 'Y') .AND. (PCOMP_LAM == 'NON')) THEN

            IF      (TYPE(1:5) == 'QUAD4') THEN

               if (type == 'QUAD4K  ') then
                  WRITE(ERR,*) ' *ERROR: Code not written for SHELL_B effect on KE yet for QUAD4K elements'
                  WRITE(F06,*) ' *ERROR :Code not written for SHELL_B effect on KE yet for QUAD4K elements'
                  call outa_here ( 'Y' )
               endif

               CALL ORDER_GAUSS ( IORD_PCOMP, SSS, HHH )

               GAUSS_PT = 0
               DO I=1,IORD_PCOMP

                  DO J=1,IORD_PCOMP

                     GAUSS_PT = GAUSS_PT + 1

                     DO K=1,3
                        DO L=1,ELDOF
                           BIG_BMI(K,L) = BIG_BM(K,L,GAUSS_PT)
                           BIG_BMI(K,L) = BIG_BM(K,L,GAUSS_PT)
                        ENDDO
                     ENDDO

                     DO K=1,3
                        DO L=1,ELDOF
                           BIG_BBI(K,L) = BIG_BB(K,L,GAUSS_PT)
                           BIG_BBI(K,L) = BIG_BB(K,L,GAUSS_PT)
                        ENDDO
                     ENDDO

                     CALL JAC2D ( SSS(I), SSS(J), XSD, YSD, 'N', JAC, JACI, DETJ )
                     CALL MATMULT_FFF ( SHELL_B, BIG_BBI, 3, 3, ELDOF, DUM1 )
                     CALL MATMULT_FFF_T ( BIG_BMI, DUM1, 3, ELDOF, ELDOF, DUM2 )
                     INTFAC = DETJ*HHH(I)*HHH(J)

                     DO K=1,ELDOF
                        DO L=1,ELDOF
                           KE(K,L) = KE(K,L) + INTFAC*(DUM2(K,L) + DUM2(L,K))
                        ENDDO
                     ENDDO

                  ENDDO

               ENDDO

            ELSE IF (TYPE(1:5) == 'TRIA3') THEN

               if (type == 'TRIA3K  ') then
                  WRITE(ERR,*) ' *ERROR: Code not written for SHELL_B effect on KE yet for TRIA3K elements'
                  WRITE(F06,*) ' *ERROR :Code not written for SHELL_B effect on KE yet for TRIA3K elements'
                  call outa_here ( 'Y' )
               endif

               DO K=1,3
                  DO L=1,ELDOF
                     BIG_BMI(K,L) = BIG_BM(K,L,1)
                     BIG_BMI(K,L) = BIG_BM(K,L,1)
                  ENDDO
               ENDDO

               DO K=1,3
                  DO L=1,ELDOF
                     BIG_BBI(K,L) = BIG_BB(K,L,1)
                     BIG_BBI(K,L) = BIG_BB(K,L,1)
                  ENDDO
               ENDDO

               CALL MATMULT_FFF ( SHELL_B, BIG_BBI, 3, 3, ELDOF, DUM1 )
               CALL MATMULT_FFF_T ( BIG_BMI, DUM1, 3, ELDOF, ELDOF, DUM2 )

               DO K=1,ELDOF
                  DO L=1,ELDOF
                     KE(K,L) = KE(K,L) + (DUM2(K,L) + DUM2(L,K))
                  ENDDO
               ENDDO

            ENDIF

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1925 FORMAT(' *ERROR  1925: ELEMENT ',I8,', TYPE ',A,', HAS ZERO OR NEGATIVE ',A,' = ',1ES9.1)

 1927 FORMAT(' *ERROR  1927: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CHAR PARAMETER QUAD4TYP MUST BE EITHER "MIN4T " OR "MIN4  " BUT IS "',A,'"')

 1948 FORMAT(' *ERROR  1948: ',A,I8,' MUST HAVE INTEGRATION ORDERS FOR PARAMS ',A,' = ',I3,' IF THE ELEMENT IS A PCOMP'            &
                             ,/,14X,' WITH SYM LAYUP. HOWEVER, THE TWO INTEGRATION ORDERS WERE: ',A,' = ',I3,' AND ',A,' = ',I3)

98761 FORMAT('   X1-X2 = ',1ES14.6,'   Y1-Y2 = ',1ES14.6)

98762 FORMAT('   X2-X3 = ',1ES14.6,'   Y2-Y3 = ',1ES14.6)

98763 FORMAT('   X3-X4 = ',1ES14.6,'   Y3-Y4 = ',1ES14.6)

98764 FORMAT('   X4-X1 = ',1ES14.6,'   Y4-Y1 = ',1ES14.6)

! **********************************************************************************************************************************

      END SUBROUTINE QDEL1

   END MODULE QUADRILATERAL_ELEMENT_DISPATCH
