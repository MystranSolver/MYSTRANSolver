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

   MODULE DKQ

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: QPLT1

   CONTAINS

      SUBROUTINE QPLT1 ( OPT, AREA, XSD, YSD )

! DKQ quadrilateral thin (Kirchoff) plate bending element. This element is based on the following work:

! "Evaluation Of A New Quadrilateral Thin Plate Bending Element", by Jean_louis Batoz and Marbrouk Ben Tahar,
! International Journal For Numerical Methods In Engineering, Vol 18 (1982) pp 1655-1677

! Element matrices calculated are:

!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  4) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  5) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MAX_ORDER_GAUSS, NSUB, NTSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, FOUR
      USE PARAMS, ONLY                :  IORQ2B
      USE MODEL_STUF, ONLY            :  ALPVEC, BE2, DT, EB, EID, KE, PRESS, PPE, PTE, SE2, STE2, SHELL_D, SHELL_DALP

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE QUADRATURE, ONLY            :  ORDER_GAUSS
      USE QUADRILATERAL_SHAPE_FUNCTIONS, ONLY:  SHP2DQ
      USE JACOBIAN, ONLY               :  JAC2D
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE QUADRILATERAL_B_MATRICES, ONLY:  BBDKQ

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'QPLT1'
      CHARACTER(46*BYTE)              :: IORD_MSG          ! Character name of the integration order (used for debug output)
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices
      INTEGER(LONG)                   :: GAUSS_PT          ! Gauss point number (used for DEBUG output in subr SHP2DQ
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG), PARAMETER        :: ID(12) =  (/ 3, & ! ID(1) =  3 means virgin 12x12 elem DOF  1 is MYSTRAN 24X24 elem DOF  3
                                                      4, & ! ID(2) =  4 means virgin 12x12 elem DOF  2 is MYSTRAN 24X24 elem DOF  4
                                                      5, & ! ID(2) =  5 means virgin 12x12 elem DOF  3 is MYSTRAN 24X24 elem DOF  5
                                                      9, & ! ID(3) =  9 means virgin 12x12 elem DOF  4 is MYSTRAN 24X24 elem DOF  9
                                                     10, & ! ID(4) = 10 means virgin 12x12 elem DOF  5 is MYSTRAN 24X24 elem DOF 10
                                                     11, & ! ID(4) = 11 means virgin 12x12 elem DOF  6 is MYSTRAN 24X24 elem DOF 11
                                                     15, & ! ID(5) = 15 means virgin 12x12 elem DOF  7 is MYSTRAN 24X24 elem DOF 15
                                                     16, & ! ID(6) = 16 means virgin 12x12 elem DOF  8 is MYSTRAN 24X24 elem DOF 16
                                                     17, & ! ID(6) = 17 means virgin 12x12 elem DOF  9 is MYSTRAN 24X24 elem DOF 17
                                                     21, & ! ID(7) = 21 means virgin 12x12 elem DOF 10 is MYSTRAN 24X24 elem DOF 21
                                                     22, & ! ID(8) = 22 means virgin 12x12 elem DOF 11 is MYSTRAN 24X24 elem DOF 22
                                                     23 /) ! ID(8) = 23 means virgin 12x12 elem DOF 12 is MYSTRAN 24X24 elem DOF 23

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Do not change IORD_STRESS_Q4. The algorithm to find Gauss point coords in elem x,y space requires there to be the same number of
! shape functions and Gauss points as elem nodes.

      INTEGER(LONG), PARAMETER        :: IORD_STRESS_Q4 = 2! Gauss integration order for stress/strain recovery matrices
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INTEGER(LONG), PARAMETER        :: NUM_NODES = 8     ! DKQ element has 8 nodes (4 are internal)
                                                           ! Indicator of no output of elem data to BUG file


      REAL(DOUBLE) , INTENT(IN)       :: AREA              ! Element area
      REAL(DOUBLE) , INTENT(IN)       :: XSD(4)            ! Diffs in x coords of quad sides in local coords
      REAL(DOUBLE) , INTENT(IN)       :: YSD(4)            ! Diffs in y coords of quad sides in local coords
      REAL(DOUBLE)                    :: ALP(3)            ! Col of ALPVEC
      REAL(DOUBLE)                    :: BB(3,12)          ! Output from subr BBDKQ, called herein. Strain-displ matrix
      REAL(DOUBLE)                    :: DETJ              ! An output from subr JAC2D, called herein. Determinant of JAC
      REAL(DOUBLE)                    :: DPSHG(2,8)        ! Output from subr SHP2DQ. Derivatives of PSH wrt elem isopar coords.
      REAL(DOUBLE)                    :: DPSHX(2,8)        ! Derivatives of PSH wrt elem x, y coords.
      REAL(DOUBLE)                    :: DUM1(3,12)        ! Intermediate matrix used in solving for PTE therm lds & KE stiffness
      REAL(DOUBLE)                    :: DUM2(12,12)       ! Intermediate matrix used in solving for KE element stiffness
      REAL(DOUBLE)                    :: DUM3(12)          ! Intermediate matrix used in solving for PTE thermal loads
      REAL(DOUBLE)                    :: DUM4(12)          ! Intermediate matrix used in solving for PTE thermal loads
      REAL(DOUBLE)                    :: EALP(3)           ! Intermed var used in calc STEi therm stress coeffs
      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS)
                                                           ! An output from subr ORDER, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: INTFAC            ! An integration factor (constant multiplier for the Gauss integration)
      REAL(DOUBLE)                    :: JAC(2,2)          ! An output from subr JAC2D, called herein. 2 x 2 Jacobian matrix.
      REAL(DOUBLE)                    :: JACI(2,2)         ! An output from subr JAC2D, called herein. 2 x 2 Jacobian inverse.
      REAL(DOUBLE)                    :: PSH(8)            ! Output from subr SHP2DQ. Shape fcn for 8 node serendipity element.
      REAL(DOUBLE)                    :: SLN(4)            ! Quad side lengths
      REAL(DOUBLE)                    :: SSI               ! A particular value of SSS
      REAL(DOUBLE)                    :: SSJ               ! A particular value of SSS
      REAL(DOUBLE)                    :: SSS(MAX_ORDER_GAUSS)
                                                           ! An output from subr ORDER, called herein. Gauss abscissa's.



! **********************************************************************************************************************************
! Quad side lengths

      DO I=1,4
         SLN(I) = DSQRT(XSD(I)*XSD(I) + YSD(I)*YSD(I))
      ENDDO

! **********************************************************************************************************************************
! Determine element thermal loads.

      IF (OPT(2) == 'Y') THEN

         CALL ORDER_GAUSS ( IORQ2B, SSS, HHH )

         DO I=1,12
            DUM3(I) = ZERO
         ENDDO

         GAUSS_PT = 0
         IORD_MSG = 'for plate bending strains,           IORQ2B = '
         DO I=1,IORQ2B
            DO J=1,IORQ2B
               GAUSS_PT = GAUSS_PT + 1
               SSI = SSS(I)
               SSJ = SSS(J)
               CALL SHP2DQ ( I, J, NUM_NODES, SUBR_NAME, IORD_MSG, IORQ2B, SSI, SSJ, 'N', PSH, DPSHG )
               CALL JAC2D ( SSI, SSJ, XSD, YSD, 'N', JAC, JACI, DETJ )
               CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 8, DPSHX )
               CALL BBDKQ ( DPSHX, XSD, YSD, SLN, I, J, 'bending strains', 'N', BB )
               CALL MATMULT_FFF_T ( BB, SHELL_DALP, 3, 12, 1, DUM4 )
               INTFAC = DETJ*HHH(I)*HHH(J)
               DO K=1,12
                  DUM3(K) = DUM3(K) + INTFAC*DUM4(K)
               ENDDO
            ENDDO
         ENDDO

         DO J=1,NTSUB
            DO K=1,12
               PTE(ID(K),J) = DUM3(K)*DT(5,J)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Determine element pressure loads.

      IF (OPT(5) == 'Y') THEN
         DO J=1,NSUB
            PPE( 3,J) = AREA*PRESS(3,J)/FOUR
            PPE( 9,J) = AREA*PRESS(3,J)/FOUR
            PPE(15,J) = AREA*PRESS(3,J)/FOUR
            PPE(21,J) = AREA*PRESS(3,J)/FOUR
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Calculate SE matrix (3 x 24) for stress data recovery.
! Note: stress recovery matrices only make sense for individual plies (or whole elem if only 1 "ply")

      IF (OPT(3) == 'Y') THEN

! Bending moment terms (calculated at center of element). Note that engineering forces (moments) are determined from
! the SE2 matrices using FCONV, which was set in calling routine and is negative.

         SSI = ZERO                                        ! BE2 at element center
         SSJ = ZERO
         GAUSS_PT = 1
         IORD_MSG = 'for plate bending strains,                  = '
         CALL SHP2DQ ( 1, 1, NUM_NODES, SUBR_NAME, IORD_MSG, 1, SSI, SSJ, 'N', PSH, DPSHG )
         CALL JAC2D ( SSI, SSJ, XSD, YSD, 'N', JAC, JACI, DETJ )
         CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 8, DPSHX )
         CALL BBDKQ ( DPSHX, XSD, YSD, SLN, 1, 1, 'bending strains', 'N', BB )

         DO I=1,3
            DO J=1,12
               BE2(I,ID(J),1) = BB(I,J)
            ENDDO
         ENDDO

! SE2, STE2 generated in elem coords. Then, in LINK9 the stresses, calc'd in elem coords, will be transformed to ply coords

         CALL MATMULT_FFF ( EB, BB, 3, 3, 12, DUM1 )       ! SE2 at element center
         DO I=1,3
            DO J=1,12
               SE2(I,ID(J),1) = DUM1(I,J)
            ENDDO
         ENDDO

         ALP(1) = ALPVEC(1,1)
         ALP(2) = ALPVEC(2,1)
         ALP(3) = ALPVEC(3,1)

         CALL MATMULT_FFF ( EB, ALP, 3, 3, 1, EALP )
         DO J=1,NTSUB
            DO I=1,3
               STE2(I,J,1) = EALP(I)*DT(5,J)
            ENDDO
         ENDDO

! Generate BE2, SE2 at Gauss points

         CALL ORDER_GAUSS ( IORQ2B, SSS, HHH )

         GAUSS_PT = 0
         DO I=1,IORD_STRESS_Q4

            DO J=1,IORD_STRESS_Q4

               GAUSS_PT = GAUSS_PT + 1
               IORD_MSG = 'for plate bending strains,   IORD_STRESS_Q4 = '
               CALL SHP2DQ ( I, J, NUM_NODES, SUBR_NAME, IORD_MSG, IORD_STRESS_Q4, SSS(I), SSS(J), 'N', PSH, DPSHG )
               CALL JAC2D ( SSS(I), SSS(J), XSD, YSD, 'N', JAC, JACI, DETJ )
               CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 8, DPSHX )
               CALL BBDKQ ( DPSHX, XSD, YSD, SLN, I, J, 'bending strains', 'N', BB )

               DO L=1,3
                  DO M=1,12
                     BE2(L,ID(M),GAUSS_PT+1) = BB(L,M)
                  ENDDO
               ENDDO

               CALL MATMULT_FFF ( EB, BB, 3, 3, 12, DUM1 )        ! Generate SE2 in element coords (at this point EM is elem coords)
               DO L=1,3
                  DO M=1,12
                     SE2(L,ID(M),GAUSS_PT+1) = DUM1(L,M)
                  ENDDO
               ENDDO

            ENDDO

         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE.

      IF(OPT(4) == 'Y') THEN

         CALL ORDER_GAUSS ( IORQ2B, SSS, HHH )

         GAUSS_PT = 0
         IORD_MSG = 'for plate bending strains,           IORQ2B = '
         DO I=1,IORQ2B
            DO J=1,IORQ2B
               SSI = SSS(I)
               SSJ = SSS(J)
               GAUSS_PT = GAUSS_PT + 1
               CALL SHP2DQ ( I, J, NUM_NODES, SUBR_NAME, IORD_MSG, IORQ2B, SSI, SSJ, 'Y', PSH, DPSHG )
               CALL JAC2D ( SSI, SSJ, XSD, YSD, 'Y', JAC, JACI, DETJ )
               CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 8, DPSHX )
               CALL BBDKQ ( DPSHX, XSD, YSD, SLN, I, J, 'bending strains', 'Y', BB )
               CALL MATMULT_FFF ( SHELL_D, BB, 3, 3, 12, DUM1 )
               CALL MATMULT_FFF_T ( BB, DUM1, 3, 12, 12, DUM2 )
               INTFAC = DETJ*HHH(I)*HHH(J)
               DO K=1,12
                  DO L=K,12
                     KE(ID(K),ID(L)) = KE(ID(K),ID(L)) + INTFAC*DUM2(K,L)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

! Set lower triangular partition equal to upper partition

         DO I=2,24
            DO J=1,I-1
               KE(I,J) = KE(J,I)
            ENDDO
         ENDDO

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE QPLT1

   END MODULE DKQ
