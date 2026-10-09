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

   MODULE MIN4

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: QPLT2

   CONTAINS

      SUBROUTINE QPLT2 ( OPT, AREA, XSD, YSD, BIG_BB )

! MIN4 quadrilateral thick (Mindlin) plate bending plate element. This element is based on the following work:

! "An Improved Treatment Of Transverse Shear In The Mindlin-Type Four-Node Quadrilateral Element", by Alexander Tessler and
!  Thomas J.R. Hughes, Computer Methods In Applied Mechanics And Engineering 39 (1983) pp 311-335

! Subroutine calculates:

!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  4) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  5) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MAX_ORDER_GAUSS, NSUB, NTSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, ONE, FOUR
      USE PARAMS, ONLY                :  EPSIL, IORQ2B, IORQ2T
      USE MODEL_STUF, ONLY            :  ALPVEC, BE2, BE3, BENSUM, DT, EID, ELDOF, EB, ET,                                         &
                                         ERR_SUB_NAM, FCONV, KE, INTL_MID, PCOMP_LAM, PCOMP_PROPS, PHI_SQ, PPE,                    &
                                         PRESS, PTE, SE2, SE3, SHELL_D, SHELL_DALP, SHELL_T, SHRSUM, STE2, TYPE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE QUADRATURE, ONLY            :  ORDER_GAUSS
      USE QUADRILATERAL_SHAPE_FUNCTIONS, ONLY:  SHP2DQ
      USE JACOBIAN, ONLY               :  JAC2D
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE MIN4_B_MATRICES, ONLY       :  BBMIN4, BSMIN4
      USE SHELL_ROTATION_SUPPORT, ONLY:  CALC_PHI_SQ

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'QPLT2'
      CHARACTER(46*BYTE)              :: IORD_MSG          ! Character name of the integration order (used for debug output)
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices

      INTEGER(LONG)                   :: GAUSS_PT          ! Gauss point number (used for DEBUG output in subr SHP2DQ
      INTEGER(LONG), PARAMETER        :: IDB( 8) = (/ 4, & ! IDB(1) =  4 means bending 8X8  elem DOF  1 is MYSTRAN 24X24 elem DOF  4
                                                      5, & ! IDB(2) =  5 means bending 8X8  elem DOF  2 is MYSTRAN 24X24 elem DOF  5
                                                     10, & ! IDB(3) = 10 means bending 8X8  elem DOF  3 is MYSTRAN 24X24 elem DOF 10
                                                     11, & ! IDB(4) = 11 means bending 8X8  elem DOF  4 is MYSTRAN 24X24 elem DOF 11
                                                     16, & ! IDB(5) = 16 means bending 8X8  elem DOF  5 is MYSTRAN 24X24 elem DOF 16
                                                     17, & ! IDB(6) = 17 means bending 8X8  elem DOF  6 is MYSTRAN 24X24 elem DOF 17
                                                     22, & ! IDB(7) = 22 means bending 8X8  elem DOF  7 is MYSTRAN 24X24 elem DOF 22
                                                     23 /) ! IDB(8) = 23 means bending 8x8  elem DOF  8 is MYSTRAN 24x24 elem DOF 23

      INTEGER(LONG), PARAMETER        :: IDS(12) = (/ 3, & ! IDS(1) =  3 means shear  12x12 elem DOF  1 is MYSTRAN 24X24 elem DOF  3
                                                      4, & ! IDS(2) =  4 means shear  12x12 elem DOF  2 is MYSTRAN 24X24 elem DOF  4
                                                      5, & ! IDS(3) =  5 means shear  12x12 elem DOF  3 is MYSTRAN 24X24 elem DOF  5
                                                      9, & ! IDS(4) =  9 means shear  12x12 elem DOF  4 is MYSTRAN 24X24 elem DOF  9
                                                     10, & ! IDS(4) = 10 means shear  12x12 elem DOF  5 is MYSTRAN 24X24 elem DOF 10
                                                     11, & ! IDS(4) = 11 means shear  12x12 elem DOF  6 is MYSTRAN 24X24 elem DOF 11
                                                     15, & ! IDS(5) = 15 means shear  12x12 elem DOF  7 is MYSTRAN 24X24 elem DOF 15
                                                     16, & ! IDS(6) = 16 means shear  12x12 elem DOF  8 is MYSTRAN 24X24 elem DOF 16
                                                     17, & ! IDS(6) = 17 means shear  12x12 elem DOF  9 is MYSTRAN 24X24 elem DOF 17
                                                     21, & ! IDS(7) = 21 means shear  12x12 elem DOF 10 is MYSTRAN 24X24 elem DOF 21
                                                     22, & ! IDS(8) = 22 means shear  12x12 elem DOF 11 is MYSTRAN 24X24 elem DOF 22
                                                     23 /) ! IDS(8) = 23 means shear  12x12 elem DOF 12 is MYSTRAN 24X24 elem DOF 23

      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Do not change IORD_STRESS_Q4. The algorithm to find Gauss point coords in elem x,y space requires there to be the same number of
! shape functions and Gauss points as elem nodes.

      INTEGER(LONG), PARAMETER        :: IORD_STRESS_Q4 = 2! Gauss integration order for stress/strain recovery matrices
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INTEGER(LONG)                   :: IORDXX            ! Gaussian integration order to use when subr ORDER is called
      INTEGER(LONG), PARAMETER        :: NUM_NODES = 4     ! Quad has 4 nodes
                                                           ! Indicator of no output of elem data to BUG file


      REAL(DOUBLE) , INTENT(IN)       :: AREA              ! Element area
      REAL(DOUBLE) , INTENT(IN)       :: XSD(4)            ! Diffs in x coords of quad sides in local coords
      REAL(DOUBLE) , INTENT(IN)       :: YSD(4)            ! Diffs in y coords of quad sides in local coords
      REAL(DOUBLE) , INTENT(OUT)      :: BIG_BB(3,ELDOF,IORQ2B*IORQ2B)
                                                           ! Strain-displ matrix for bending for all Gauss points and all DOF's

      REAL(DOUBLE)                    :: ALP(3)            ! Part of a col of SHELL_DALP
      REAL(DOUBLE)                    :: BB(3,8)           ! Strain-displ matrix for bending. Output from subr BBMIN4
      REAL(DOUBLE)                    :: BS(2,12)          ! Strain-displ matrix for bending. Output from subr BBMIN4
      REAL(DOUBLE)                    :: DETJ              ! Determinant of JAC. An output from subr JAC2D.
      REAL(DOUBLE)                    :: DNXSHG(2,4)       ! Deriv's of NXSH wrt xi, eta elem isopar coords. Output from subr MIN4SH
      REAL(DOUBLE)                    :: DNYSHG(2,4)       ! Deriv's of NYSH wrt xi, eta elem isopar coords. Output from subr MIN4SH
      REAL(DOUBLE)                    :: DNXSHX(2,4)       ! Deriv's of NXSH wrt x, y elem local coords. Output from subr JAC2D.
      REAL(DOUBLE)                    :: DNYSHX(2,4)       ! Deriv's of NYSH wrt  x, y elem local coords. Output from subr JAC2D.
      REAL(DOUBLE)                    :: DPSHG(2,4)        ! Derivatives of PSH wrt elem isopar coords. Output from subr SHP2DQ.
      REAL(DOUBLE)                    :: DPSHX(2,4)        ! Derivatives of PSH wrt elem x, y elem local coords.
      REAL(DOUBLE)                    :: DUM1(3,8)         ! Intermediate result in calc KB bend stiff & SEi stress recov matrices
      REAL(DOUBLE)                    :: DUM2(8,8)         ! Intermediate result in calc KB elem bending stiffness
      REAL(DOUBLE)                    :: DUM3(2,12)        ! Intermediate result in calc KS shear stiff & SEi stress recov matrices
      REAL(DOUBLE)                    :: DUM4(12,12)       ! Intermediate result in calc KS elem bending stiffness
      REAL(DOUBLE)                    :: DUM6(8)           ! Intermediate result in calc PTE elem thermal loads
      REAL(DOUBLE)                    :: DUM8(12)          ! Intermediate result in calc PPE elem pressure loads
      REAL(DOUBLE)                    :: DUM7(8)           ! Intermediate result in calc PTE elem thermal loads
      REAL(DOUBLE)                    :: DUM9(12)          ! Intermediate result in calc PPE elem pressure loads
      REAL(DOUBLE)                    :: EALP(3)           ! Intermed var used in calc STEi therm stress coeffs
      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS)
                                                           ! Gauss weights. An output from subr ORDER.
      REAL(DOUBLE)                    :: INTFAC            ! An integration factor (constant multiplier for the Gauss integration)
      REAL(DOUBLE)                    :: JAC(2,2)          ! 2 x 2 Jacobian matrix. An output from subr JAC2D.
      REAL(DOUBLE)                    :: JACI(2,2)         ! 2 x 2 Jacobian inverse. An output from subr JAC2D.
      REAL(DOUBLE)                    :: KB(8,8)           ! Bending stiffness contribution to KE for the MIN4 elem
      REAL(DOUBLE)                    :: KS(12,12)         ! PHISQ*KS is the shear stiff contribution to KE for the MIN4 elem
      REAL(DOUBLE)                    :: NXSH(4)           ! Constrained Nx shape functions. Output from subr MIN4SH.
      REAL(DOUBLE)                    :: NYSH(4)           ! Constrained Ny shape functions. Output from subr MIN4SH.
      REAL(DOUBLE)                    :: PSH(4)            ! Shape fcn at Gauss pts SSI, SSJ. Output from subr SHP2DQ.
      REAL(DOUBLE)                    :: SDETJ             ! Sum of DETJ's when calc BS (to get Gauss weighted BS matrix)
      REAL(DOUBLE)                    :: SSI               ! A particular value of SSS
      REAL(DOUBLE)                    :: SSJ               ! A particular value of SSS
      REAL(DOUBLE)                    :: SSS(MAX_ORDER_GAUSS)
                                                           ! Gauss abscissa's. An output from subr ORDER.

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Determine element thermal loads.

      IF (OPT(2) == 'Y') THEN

         CALL ORDER_GAUSS ( IORQ2B, SSS, HHH )

         DO I=1,8
            DUM7(I) = ZERO
         ENDDO

         GAUSS_PT = 0
         IORD_MSG = 'for plate bending strains,           IORQ2B = '
         DO I=1,IORQ2B
            DO J=1,IORQ2B
               GAUSS_PT = GAUSS_PT + 1
               SSI  = SSS(I)
               SSJ  = SSS(J)
               CALL SHP2DQ ( I, J, NUM_NODES, SUBR_NAME, IORD_MSG, IORQ2B, SSI, SSJ, 'N', PSH, DPSHG )
               CALL JAC2D ( SSI, SSJ, XSD, YSD, 'N', JAC, JACI, DETJ )
               CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 4, DPSHX )
               CALL BBMIN4 ( DPSHX, I, J, 'bending strains', 'N', BB )
               CALL MATMULT_FFF_T ( BB, SHELL_DALP, 3, 8, 1, DUM6 )
               INTFAC = DETJ*HHH(I)*HHH(J)
               DO K=1,8
                  DUM7(K) = DUM7(K) + INTFAC*DUM6(K)
               ENDDO
            ENDDO
         ENDDO

         DO I=1,8
            DO J=1,NTSUB
               PTE(IDB(I),J) = DUM7(I)*DT(5,J)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Determine element pressure loads

      IF (OPT(5) == 'Y') THEN

         IF (DEBUG(16) == 0) THEN                          ! Generate PPE as work equilavent loads

            CALL ORDER_GAUSS ( IORQ2B, SSS, HHH )

            DO I=1,12
               DUM9(I) = ZERO
            ENDDO

            GAUSS_PT = 0
            IORD_MSG = 'for plate bending strains,           IORQ2B = '
            DO I=1,IORQ2B
               DO J=1,IORQ2B
                  SSI = SSS(I)
                  SSJ = SSS(J)
                  GAUSS_PT = GAUSS_PT + 1
                  CALL SHP2DQ ( I, J, NUM_NODES, SUBR_NAME, IORD_MSG, IORQ2B, SSI, SSJ, 'N', PSH, DPSHG )
                  CALL MIN4SH ( SSI, SSJ, XSD, YSD, 'N', NXSH, NYSH, DNXSHG, DNYSHG )
                  CALL JAC2D ( SSI, SSJ, XSD, YSD, 'N', JAC,  JACI, DETJ )
                  DUM8( 1) =  PSH(1)
                  DUM8( 2) = -NXSH(1)
                  DUM8( 3) =  NYSH(1)

                  DUM8( 4) =  PSH(2)
                  DUM8( 5) = -NXSH(2)
                  DUM8( 6) =  NYSH(2)

                  DUM8( 7) =  PSH(3)
                  DUM8( 8) = -NXSH(3)
                  DUM8( 9) =  NYSH(3)

                  DUM8(10) =  PSH(4)
                  DUM8(11) = -NXSH(4)
                  DUM8(12) =  NYSH(4)
                  INTFAC   =  DETJ*HHH(I)*HHH(J)
                  DO K=1,12
                     DUM9(K) = DUM9(K) + INTFAC*DUM8(K)
                  ENDDO
               ENDDO
            ENDDO

            DO I=1,12
               DO J=1,NSUB
                  PPE(IDS(I),J) = DUM9(I)*PRESS(3,J)
               ENDDO
            ENDDO

         ELSE                                              ! Generate PPE as static equilavent loads

            DO J=1,NSUB
               PPE( 3,J) = AREA*PRESS(3,J)/FOUR
               PPE( 9,J) = AREA*PRESS(3,J)/FOUR
               PPE(15,J) = AREA*PRESS(3,J)/FOUR
               PPE(21,J) = AREA*PRESS(3,J)/FOUR
            ENDDO

         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE. Note that we need to calc KE if the stress recovery matrices (for OPT(3)) are to be cal'd
! since the BE strain recovery matrices use PHI_SQ

      IF ((OPT(4) == 'Y') .OR. (OPT(3) == 'Y')) THEN

! Bending stiffness terms:

         CALL ORDER_GAUSS ( IORQ2B, SSS, HHH )

         DO I=1,8
            DO J=1,8
               KB(I,J) = ZERO
            ENDDO
         ENDDO

         GAUSS_PT = 0
         IORD_MSG = 'for plate bending strains,           IORQ2B = '

         DO I=1,IORQ2B
            DO J=1,IORQ2B
               SSI = SSS(I)
               SSJ = SSS(J)
               GAUSS_PT = GAUSS_PT + 1
               CALL SHP2DQ ( I, J, NUM_NODES, SUBR_NAME, IORD_MSG, IORQ2B, SSI, SSJ, 'Y', PSH, DPSHG )
               CALL JAC2D ( SSI, SSJ, XSD, YSD, 'Y', JAC, JACI, DETJ )
               CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 4, DPSHX )
               CALL BBMIN4 ( DPSHX, I, J, 'bending strains', 'Y', BB )
               CALL MATMULT_FFF ( SHELL_D, BB, 3, 3, 8, DUM1 )
               CALL MATMULT_FFF_T ( BB, DUM1, 3, 8, 8, DUM2 )
               INTFAC = DETJ*HHH(I)*HHH(J)
               DO K=1,8
                  DO L=K,8
                     KB(K,L) = KB(K,L) + INTFAC*DUM2(K,L)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         BENSUM = ZERO                                     ! Add all diag terms from KB
         DO I=1,8
            BENSUM = BENSUM + KB(I,I)
         ENDDO

! Shear stiffness terms. If IORQ2T is positive, then normal Gaussian integration of order IORQ2T is performed.
! If IORQ2T is neqative, a Jacobian weighted value for BS is used for all of the Gauss pts.

         IF (IORQ2T > 0) THEN
            IORDXX =  IORQ2T
         ELSE
            IORDXX = -IORQ2T
         ENDIF
         CALL ORDER_GAUSS ( IORDXX, SSS, HHH )

         DO I=1,12
            DO J=1,12
               KS(I,J) = ZERO
            ENDDO
         ENDDO

         IF (IORQ2T > 0) THEN
            GAUSS_PT = 0
            IORD_MSG = 'for transverse shear strains,        IORQ2T = '
            DO I=1,IORDXX
               DO J=1,IORDXX
                  SSI = SSS(I)
                  SSJ = SSS(J)
                  GAUSS_PT = GAUSS_PT + 1
                  CALL SHP2DQ ( I, J, NUM_NODES, SUBR_NAME, IORD_MSG, IORQ2T, SSI, SSJ, 'Y', PSH, DPSHG )
                  CALL MIN4SH ( SSI, SSJ, XSD, YSD, 'Y', NXSH, NYSH, DNXSHG, DNYSHG )
                  CALL JAC2D ( SSI, SSJ, XSD, YSD, 'Y', JAC, JACI, DETJ )
                  CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 4, DPSHX )
                  CALL MATMULT_FFF ( JACI, DNXSHG, 2, 2, 4, DNXSHX )
                  CALL MATMULT_FFF ( JACI, DNYSHG, 2, 2, 4, DNYSHX )
                  CALL BSMIN4 ( PSH, DPSHX, DNXSHX, DNYSHX, I, J, 'transverse shear strains', 'Y', BS )
                  CALL MATMULT_FFF ( SHELL_T, BS, 2, 2, 12, DUM3 )
                  CALL MATMULT_FFF_T ( BS, DUM3, 2, 12, 12, DUM4 )
                  INTFAC = DETJ*HHH(I)*HHH(J)
                  DO K=1,12
                     DO L=K,12
                        KS(K,L) = KS(K,L) + INTFAC*DUM4(K,L)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO

         ELSE

            DO I=1,2
               DO J=1,12
                  BS(I,J) = ZERO
               ENDDO
            ENDDO
            SDETJ = ZERO
            GAUSS_PT = 0
            IORD_MSG = 'for transverse shear strains,        IORQ2T = '
            DO I=1,IORDXX
               DO J=1,IORDXX
                  SSI = SSS(I)
                  SSJ = SSS(J)
                  GAUSS_PT = GAUSS_PT + 1
                  CALL SHP2DQ ( I, J, NUM_NODES, SUBR_NAME, IORD_MSG, IORQ2T, SSI, SSJ, 'Y', PSH, DPSHG )
                  CALL MIN4SH ( SSI, SSJ, XSD, YSD, 'Y', NXSH, NYSH, DNXSHG, DNYSHG )
                  CALL JAC2D ( SSI, SSJ, XSD, YSD, 'Y', JAC, JACI, DETJ )
                  SDETJ = SDETJ + DETJ
                  CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 4, DPSHX )
                  CALL MATMULT_FFF ( JACI, DNXSHG, 2, 2, 4, DNXSHX )
                  CALL MATMULT_FFF ( JACI, DNYSHG, 2, 2, 4, DNYSHX )
                  CALL BSMIN4 ( PSH, DPSHX, DNXSHX, DNYSHX, I, J, 'transverse shear strains', 'Y', DUM3 )
                  DO K=1,2
                     DO L=1,12
                        BS(K,L) = BS(K,L) + DETJ*DUM3(K,L)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO

            DO I=1,2
               DO J=1,12
                  BS(I,J) = BS(I,J)/SDETJ
               ENDDO
            ENDDO

            DO I=1,IORDXX
               DO J=1,IORDXX
                  SSI = SSS(I)
                  SSJ = SSS(J)
                  CALL JAC2D ( SSI, SSJ, XSD, YSD, 'N', JAC, JACI, DETJ )
                  CALL MATMULT_FFF ( SHELL_T, BS, 2, 2, 12, DUM3 )
                  CALL MATMULT_FFF_T ( BS, DUM3, 2, 12, 12, DUM4 )
                  INTFAC = DETJ*HHH(I)*HHH(J)
                  DO K=1,12
                     DO L=K,12
                        KS(K,L) = KS(K,L) + INTFAC*DUM4(K,L)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO

         ENDIF
                                                        ! Add all diagonal terms from KS for rotational DOF's to get SHRSUM
         SHRSUM = KS(2,2) + KS(3,3) + KS(5,5) + KS(6,6) + KS(8,8) + KS(9,9) + KS(11,11) + KS(12,12)

! Now calculate the finite elem shear factor, PHI_SQ

         CALL CALC_PHI_SQ ( IERROR )

! Return if IERROR > 0

         IF (IERROR > 0) RETURN

         DO I=1,8
            DO J=I,8
               KE(IDB(I),IDB(J)) = KE(IDB(I),IDB(J)) + KB(I,J)
            ENDDO
         ENDDO

         DO I=1,12
            DO J=I,12
               KE(IDS(I),IDS(J)) = KE(IDS(I),IDS(J)) + PHI_SQ*KS(I,J)
            ENDDO
         ENDDO

! Set lower triangular partition equal to upper partition

         DO I=2,24
            DO J=1,I-1
               KE(I,J) = KE(J,I)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate BE, SE matrices (3 x 24) for strain/stress data recovery.
! Note: stress recovery matrices only make sense for individual plies (or whole elem if only 1 "ply")

      IF (OPT(3) == 'Y') THEN

! Bending moment terms (calculated at center of element). Note that engineering forces (moments) are determined from
! the SE2 matrices using FCONV, which was set in calling routine and is negative.

         SSI = ZERO
         SSJ = ZERO
         GAUSS_PT = 1
         IORD_MSG = 'for plate bending strains,                  = '
         CALL SHP2DQ ( 1, 1, NUM_NODES, SUBR_NAME, IORD_MSG, 1, SSI, SSJ, 'N', PSH, DPSHG )
         CALL JAC2D ( SSI, SSJ, XSD, YSD, 'N', JAC, JACI, DETJ )
         CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 4, DPSHX )
         CALL BBMIN4 ( DPSHX, 1, 1, 'bending strains', 'N', BB )

         DO I=1,3                                          ! Bending strain-displ matrix
            DO J=1,8
               BE2(I,IDB(J),1) = BB(I,J)
            ENDDO
         ENDDO

! SE2, STE2 generated in elem coords. Then, in LINK9 the stresses, calc'd in elem coords, will be transformed to ply coords

         CALL MATMULT_FFF ( EB, BB, 3, 3, 8, DUM1 )        ! Generate SE2 in element coords (at this point EB is elem coords)
         DO I=1,3
            DO J=1,8
               SE2(I,IDB(J),1) = DUM1(I,J)
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
               CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 4, DPSHX )
               CALL BBMIN4 ( DPSHX, I, J, 'bending strains', 'N', BB )

               DO L=1,3
                  DO M=1,8
                     BE2(L,IDB(M),GAUSS_PT+1) = BB(L,M)
                  ENDDO
               ENDDO

               CALL MATMULT_FFF ( EB, BB, 3, 3, 8, DUM1 )  ! Generate SE2 in element coords (at this point EM is elem coords)
               DO L=1,3
                  DO M=1,8
                     SE2(L,IDB(M),GAUSS_PT+1) = DUM1(L,M)
                  ENDDO
               ENDDO

            ENDDO

         ENDDO

! Transverse shear terms. Calculate one set based on average at four Gauss points and 4 others at the Gauss points.

         IORDXX = 2
         CALL ORDER_GAUSS ( IORDXX, SSS, HHH )

         GAUSS_PT = 0
         DO I=1,2
            DO J=1,2

               SSI = SSS(I)
               SSJ = SSS(J)
               GAUSS_PT = GAUSS_PT + 1
               IORD_MSG = 'for transverse shear strains,               = '
               CALL SHP2DQ ( I, J, NUM_NODES, SUBR_NAME, IORD_MSG, 2, SSI, SSJ, 'N', PSH, DPSHG )
               CALL MIN4SH ( SSI, SSJ, XSD, YSD, 'N', NXSH, NYSH, DNXSHG, DNYSHG )
               CALL JAC2D ( SSI, SSJ, XSD, YSD, 'N', JAC, JACI, DETJ )
               CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 4, DPSHX )
               CALL MATMULT_FFF ( JACI, DNXSHG, 2, 2, 4, DNXSHX )
               CALL MATMULT_FFF ( JACI, DNYSHG, 2, 2, 4, DNYSHX )
               CALL BSMIN4 ( PSH, DPSHX, DNXSHX, DNYSHX, I, J, 'transverse shear strains', 'N', BS )

               DO K=1,2                                    ! Transverse shear strain-displ terms - avg of 4 Gauss points
                  DO L=1,12
                     BE3(K,IDS(L),1) = BE3(K,IDS(L),1) + BS(K,L)/FOUR
                  ENDDO
               ENDDO

               DO K=1,2                                    ! Transverse shear strain-displ terms - one for each Gauss point
                  DO L=1,12
                     BE3(K,IDS(L),GAUSS_PT+1) = BS(K,L)
                  ENDDO
               ENDDO

! SE3 generated in elem coords. Then, in LINK9, the stresses calc'd in elem coords, will be transformed to ply coords

               CALL MATMULT_FFF ( ET, BS, 2, 2, 12, DUM3 ) ! Transverse shear stress-displ terms

               DO K=1,2
                  DO L=1,12
                     SE3(K,IDS(L),1) = SE3(K,IDS(L),1) + PHI_SQ*DUM3(K,L)/FOUR
                  ENDDO
               ENDDO

               DO K=1,2
                  DO L=1,12
                     SE3(K,IDS(L),GAUSS_PT+1) = PHI_SQ*DUM3(K,L)
                  ENDDO
               ENDDO

            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! If element is a composite and if it is a nonsym layup we need to calc BIG_BB for later use

      DO K=1,3
         DO L=1,ELDOF
            DO M=1,IORQ2B*IORQ2B
               BIG_BB(K,L,M) = ZERO
            ENDDO
         ENDDO
      ENDDO

      IF ((PCOMP_PROPS == 'Y') .AND. (PCOMP_LAM == 'NON')) THEN
         CALL ORDER_GAUSS ( IORQ2B, SSS, HHH )
         GAUSS_PT = 0
         IORD_MSG = 'for plate bending strains,           IORQ2B = '
         DO I=1,IORQ2B
            DO J=1,IORQ2B
               SSI = SSS(I)
               SSJ = SSS(J)
               GAUSS_PT = GAUSS_PT + 1
               CALL SHP2DQ ( I, J, NUM_NODES, SUBR_NAME, IORD_MSG, IORQ2B, SSI, SSJ, 'N', PSH, DPSHG )
               CALL JAC2D ( SSI, SSJ, XSD, YSD, 'Y', JAC, JACI, DETJ )
               CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 4, DPSHX )
               CALL BBMIN4 ( DPSHX, I, J, 'bending strains', 'N', BB )
               DO L=1,3
                  DO M=1,8
                     BIG_BB(L,IDB(M),GAUSS_PT) = BB(L,M)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE QPLT2


      SUBROUTINE MIN4SH ( SSI, SSJ, XSD, YSD, WRT_BUG_THIS_TIME, NXSH, NYSH, DNXSHG, DNYSHG )

! Generates constrained shapes for use with MIN4 plate element for transverse shear. This is based on the paper referenced in subr
! QPLT2.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, F06, WRT_BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, EIGHT

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MIN4SH'
      CHARACTER(17*BYTE)              :: NAME(2)           ! Used for BUG output annotation
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG)                   :: I,J               ! DO loop indices


      REAL(DOUBLE) , INTENT(IN)       :: SSI               ! Gauss point coordinate
      REAL(DOUBLE) , INTENT(IN)       :: SSJ               ! Gauss point coordinate
      REAL(DOUBLE) , INTENT(IN)       :: XSD(4)            ! 1-D arrays of differences in x side dimensions (local)
      REAL(DOUBLE) , INTENT(IN)       :: YSD(4)            ! 1-D arrays of differences in y side dimensions (local)
      REAL(DOUBLE) , INTENT(OUT)      :: NXSH(4)           ! Constrained Nx shape functions for the MIN4 quad element
      REAL(DOUBLE) , INTENT(OUT)      :: NYSH(4)           ! Constrained Ny shape functions for the MIN4 quad element
      REAL(DOUBLE) , INTENT(OUT)      :: DNXSHG(2,4)       ! Derivatives of NXSH wrt xi, eta.
      REAL(DOUBLE) , INTENT(OUT)      :: DNYSHG(2,4)       ! Derivatives of NYSH wrt xi, eta.
      REAL(DOUBLE)                    :: N5,N6,N7,N8       ! N5 thru N8 are the virgin shape functions used in finding the
                                                             ! constrained shape functions (NXSH, NYSH)
      REAL(DOUBLE)                    :: N5X,N6X,N7X,N8X   ! N5X thru N8X are derivatives of N5 thru N8 wrt xi
      REAL(DOUBLE)                    :: N5Y,N6Y,N7Y,N8Y   ! N5Y thru N8Y are derivatives of N5 thru N8 wrt eta
      REAL(DOUBLE)                    :: X12               ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: X23               ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: X34               ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: X41               ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: Y12               ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: Y23               ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: Y34               ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: Y41               ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: XM                ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: XP                ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: X2M               ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: YM                ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: YP                ! Intermediate variable used in calculating outputs
      REAL(DOUBLE)                    :: Y2M               ! Intermediate variable used in calculating outputs


! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,4
         NXSH(I) = ZERO
         NYSH(I) = ZERO
      ENDDO

      DO I=1,2
         DO J=1,4
            DNXSHG(I,J) = ZERO
            DNYSHG(I,J) = ZERO
         ENDDO
      ENDDO

! Get constants needed from local G.P. geometry. Xij = Xi - Xj  and  Yij = Yi - Yj where Xi, Yi, etc., are element
! node coords in local element coord system

      X12 = XSD(1)
      X23 = XSD(2)
      X34 = XSD(3)
      X41 = XSD(4)

      Y12 = YSD(1)
      Y23 = YSD(2)
      Y34 = YSD(3)
      Y41 = YSD(4)

      XM  = ONE - SSI
      XP  = ONE + SSI
      YM  = ONE - SSJ
      YP  = ONE + SSJ
      X2M = ONE - SSI*SSI
      Y2M = ONE - SSJ*SSJ

! N5 thru N8 are the virgin shape functions used in finding the  constrained shape functions (NXSH, NYSH)

      N5 = X2M*YM/TWO
      N6 = Y2M*XP/TWO
      N7 = X2M*YP/TWO
      N8 = Y2M*XM/TWO

! N5X thru N8Y are derivatives of N5 thru N8 wrt xi, eta

      N5X = -SSI*YM
      N6X =  Y2M/TWO
      N7X = -SSI*YP
      N8X = -Y2M/TWO

      N5Y = -X2M/TWO
      N6Y = -SSJ*XP
      N7Y =  X2M/TWO
      N8Y = -SSJ*XM

! Constrained shapes:

      NXSH(1) = (-Y41*N8 + Y12*N5)/EIGHT
      NXSH(2) = (-Y12*N5 + Y23*N6)/EIGHT
      NXSH(3) = (-Y23*N6 + Y34*N7)/EIGHT
      NXSH(4) = (-Y34*N7 + Y41*N8)/EIGHT

      NYSH(1) = (-X41*N8 + X12*N5)/EIGHT
      NYSH(2) = (-X12*N5 + X23*N6)/EIGHT
      NYSH(3) = (-X23*N6 + X34*N7)/EIGHT
      NYSH(4) = (-X34*N7 + X41*N8)/EIGHT

! Derivatives of NXSH wrt xi, eta:

      DNXSHG(1,1) = (-Y41*N8X + Y12*N5X)/EIGHT
      DNXSHG(1,2) = (-Y12*N5X + Y23*N6X)/EIGHT
      DNXSHG(1,3) = (-Y23*N6X + Y34*N7X)/EIGHT
      DNXSHG(1,4) = (-Y34*N7X + Y41*N8X)/EIGHT

      DNXSHG(2,1) = (-Y41*N8Y + Y12*N5Y)/EIGHT
      DNXSHG(2,2) = (-Y12*N5Y + Y23*N6Y)/EIGHT
      DNXSHG(2,3) = (-Y23*N6Y + Y34*N7Y)/EIGHT
      DNXSHG(2,4) = (-Y34*N7Y + Y41*N8Y)/EIGHT

! Derivatives of NYSH wrt xi, eta:

      DNYSHG(1,1) = (-X41*N8X + X12*N5X)/EIGHT
      DNYSHG(1,2) = (-X12*N5X + X23*N6X)/EIGHT
      DNYSHG(1,3) = (-X23*N6X + X34*N7X)/EIGHT
      DNYSHG(1,4) = (-X34*N7X + X41*N8X)/EIGHT

      DNYSHG(2,1) = (-X41*N8Y + X12*N5Y)/EIGHT
      DNYSHG(2,2) = (-X12*N5Y + X23*N6Y)/EIGHT
      DNYSHG(2,3) = (-X23*N6Y + X34*N7Y)/EIGHT
      DNYSHG(2,4) = (-X34*N7Y + X41*N8Y)/EIGHT

! **********************************************************************************************************************************
! Debug output:

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(7) > 0)) THEN

         NAME(1) = 'Nodes  1 thru  4:'
         NAME(2) = 'Nodes  5 thru  8:'

         WRITE(BUG,1112)

         WRITE(BUG,1104)
         WRITE(BUG,1105) NAME(2), N5, N6, N7, N8
         WRITE(BUG,*)
         WRITE(BUG,*)

         WRITE(BUG,1115)
         WRITE(BUG,1105) NAME(2), N5X, N6X, N7X, N8X
         WRITE(BUG,*)
         WRITE(BUG,*)

         WRITE(BUG,1125)
         WRITE(BUG,1105) NAME(2), N5Y, N6Y, N7Y, N8Y
         WRITE(BUG,*)
         WRITE(BUG,*)

         WRITE(BUG,1204)
         WRITE(BUG,1105) NAME(2), (NXSH(I), I=1,4)
         WRITE(BUG,*)
         WRITE(BUG,*)

         WRITE(BUG,1205)
         WRITE(BUG,1105) NAME(2), (NYSH(I),I=1,4)
         WRITE(BUG,*)
         WRITE(BUG,*)

         WRITE(BUG,1304)
         WRITE(BUG,1105) NAME(1), (DNXSHG(1,J), J=1,4)
         WRITE(BUG,1105) NAME(2), (DNXSHG(2,J), J=1,4)
         WRITE(BUG,*)
         WRITE(BUG,*)

         WRITE(BUG,1305)
         WRITE(BUG,1105) NAME(1), (DNYSHG(1,J), J=1,4)
         WRITE(BUG,1105) NAME(2), (DNYSHG(2,J), J=1,4)
         WRITE(BUG,*)
         WRITE(BUG,*)

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1102 FORMAT(45X,A,' Element ',I8,', Gauss point ',I3,//)

 1112 format(51x,'Outputs from subroutine MIN4SH',/,51x,'------------------------------',/)

 1104 FORMAT(48X,'Unconstrained shape functions N5 - N8')

 1105 FORMAT(23X,A,1X,4(1ES16.6))

 1115 FORMAT(35X,'Derivatives of unconstrained shape functions with respect to XI')

 1125 FORMAT(35X,'Derivatives of unconstrained shape functions with respect to ET')

 1204 FORMAT(49X,'Constrained shape functions NXSH(i)')

 1205 FORMAT(51X,'Constrained shape functions NYSH(i)')

 1304 FORMAT(36X,'Derivatives of constrained shape functions with respect to XI')

 1305 FORMAT(36X,'Derivatives of constrained shape functions with respect to ET')

! **********************************************************************************************************************************

      END SUBROUTINE MIN4SH

   END MODULE MIN4
