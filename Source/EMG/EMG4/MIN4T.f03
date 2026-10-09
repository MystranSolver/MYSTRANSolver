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

   MODULE MIN4T

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: QPLT3

   CONTAINS

      SUBROUTINE QPLT3 ( OPT, AREA_QUAD, XSD, YSD, BIG_BB )

! MIN4T quadrilateral thick (Mindlin) plate bending plate element. This is made of 4 non-overlapping MIN3 trianglau elements with
! the central point reduced out using either Kirchoff restraints or static condensation (user controlled via PARAM MIN4TRED).
! This element is based on the following work:

! "A four-node, shear deformable shell element developed via explicit Kirchoff constraints", Int J For Numerical Methods In Engr, 49
!  by Alex Tessler, et al

! The quad is made up of 4 non-overlapping MIN3 triangular elements as shown below

!                         node 4 -------------------- node 3
!                               |\                  /|
!                               | \                / |
!                               |  \    tria      /  |  MIN3 triangle 1 is connected to nodes 1-2-5
!                               |   \     3      /   |
!                               |    \          /    |  MIN3 triangle 2 is connected to nodes 2-3-5
!                               |     \        /     |
!                               |      \      /      |  MIN3 triangle 3 is connected to nodes 3-4-5
!                               |       \node/       |
!                               | tria   \5 / tria   |  MIN3 triangle 4 is connected to nodes 4-1-5
!                               |   4     \/    2    |
!                               |         /\         |  Node 5 is reduced out using:
!                               |        /  \        |      1) Kirchoff constraints if MIN4TRED = 'B54' (see module PARAMS)
!                               |       /    \       |      2) static condensation  if MIN4TRED = 'STC' (see module PARAMS)
!                               |      /      \      |
!                               |     /        \     |
!                               |    /  tria    \    |
!                               |   /     1      \   |
!                               |  /              \  |
!                               | /                \ |
!                               |/                  \|
!                        node 1  -------------------- node 2

! Subroutine calculates (by calling TPLT2):

!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  4) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  5) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MEFE, MIN4T_QUAD4_TRIA_NO, NSUB, NTSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, QUARTER, HALF, ONE, TWO, FOUR, CONV_RAD_DEG, PI
      USE PARAMS, ONLY                :  MIN4TRED
      USE MACHINE_PARAMS, ONLY        :  MACH_SFMIN
      USE MODEL_STUF, ONLY            :  BE2, BE3, DT, EID, ELDOF, EMG_IFE, ERR_SUB_NAM,NUM_EMG_FATAL_ERRS,                        &
                                         FCONV, KE, PHI_SQ, PPE, PTE, SE2, SE3, TE, XEB, XEL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FULL_MATRIX_ALGEBRA, ONLY   :  CHECK_MAT_INVERSE, MATADD_FFF, MATMULT_FFF, MATMULT_FFF_T
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE TRIANGULAR_SHELL_ELEMENTS, ONLY:  TPLT2

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'QPLT3'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices.
      CHARACTER(1*BYTE)               :: OPT_MIN4T(6)      ! Values of OPT to use in this subr. We need OPT(4) = 'Y' if the
!                                                            MIN4T reduction technique is static condensation (MIN4TRED = 'STC')

      CHARACTER(LEN=1), PARAMETER     :: MN4T_QD   = 'Y'   ! Arg used in call to TPLT2 to say the triangular elem is part of a QUAD4

      INTEGER(LONG), PARAMETER        :: IDV(9)  = (/ 1, & ! IDV( 1) =  1 means tria elem virgin DOF 1 is MYSTRAN elem DOF  1
                                                      4, & ! IDV( 2) =  4 means tria elem virgin DOF 2 is MYSTRAN elem DOF  4
                                                      7, & ! IDV( 3) =  7 means tria elem virgin DOF 3 is MYSTRAN elem DOF  7
                                                      2, & ! IDV( 4) =  2 means tria elem virgin DOF 4 is MYSTRAN elem DOF  2
                                                      5, & ! IDV( 5) =  5 means tria elem virgin DOF 5 is MYSTRAN elem DOF  5
                                                      8, & ! IDV( 6) =  8 means tria elem virgin DOF 6 is MYSTRAN elem DOF  8
                                                      3, & ! IDV( 7) =  3 means tria elem virgin DOF 7 is MYSTRAN elem DOF  3
                                                      6, & ! IDV( 8) =  6 means tria elem virgin DOF 8 is MYSTRAN elem DOF  6
                                                      9 /) ! IDV( 9) =  9 means tria elem virgin DOF 9 is MYSTRAN elem DOF  9

      INTEGER(LONG)                   :: IDI(4,9)

      INTEGER(LONG), PARAMETER        :: IDM(12) = (/ 3, & ! IDM( 1) =  3 means quad elem DOF  1 is MYSTRAN elem DOF  3
                                                      4, & ! IDM( 2) =  4 means quad elem DOF  2 is MYSTRAN elem DOF  4
                                                      5, & ! IDM( 3) =  5 means quad elem DOF  3 is MYSTRAN elem DOF  5
                                                      9, & ! IDM( 4) =  8 means quad elem DOF  4 is MYSTRAN elem DOF  9
                                                     10, & ! IDM( 5) = 10 means quad elem DOF  5 is MYSTRAN elem DOF 10
                                                     11, & ! IDM( 6) = 11 means quad elem DOF  6 is MYSTRAN elem DOF 11
                                                     15, & ! IDM( 7) = 15 means quad elem DOF  7 is MYSTRAN elem DOF 15
                                                     16, & ! IDM( 8) = 16 means quad elem DOF  8 is MYSTRAN elem DOF 16
                                                     17, & ! IDM( 9) = 17 means quad elem DOF  9 is MYSTRAN elem DOF 17
                                                     21, & ! IDM(10) = 21 means quad elem DOF 10 is MYSTRAN elem DOF 21
                                                     22, & ! IDM(11) = 22 means quad elem DOF 11 is MYSTRAN elem DOF 22
                                                     23 /) ! IDM(12) = 23 means quad elem DOF 12 is MYSTRAN elem DOF 23

      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: II,JJ,KK          ! DO loop indices or counters
      INTEGER(LONG)                   :: IERROR    = 0     ! Local error indicator

      INTEGER(LONG), PARAMETER        :: NUM_TRIAS = 4     ! DO NOT CHANGE THIS. Num of triangles that subdivide the QUAD4
!                                                            ------------------

      INTEGER(LONG)                   :: PROG_ERR  = 0     ! Local error indicator
      INTEGER(LONG)                   :: TRIA_NUM          ! 1, 2, 3, or 4 designator of a subtriangle of the quad


      REAL(DOUBLE) , INTENT(IN)       :: AREA_QUAD         ! Element area
      REAL(DOUBLE) , INTENT(IN)       :: XSD(4)            ! Diffs in x coords of quad sides in local coords
      REAL(DOUBLE) , INTENT(IN)       :: YSD(4)            ! Diffs in y coords of quad sides in local coords

                                                           ! Strain-displ matrix for bending for all DOF's
      REAL(DOUBLE) , INTENT(OUT)      :: BIG_BB(3,ELDOF,4)

      REAL(DOUBLE)                    :: ALPHA(4)          ! Angle (rad) between side 1-2 of a triangle and side 1-3 of the triangle
      REAL(DOUBLE)                    :: BETA(4)           ! Angle (rad) between side 1-2 of a triangle and side 2-3 of the triangle
      REAL(DOUBLE)                    :: GAMMA(4)          ! Angle (rad) between side 1-3 of a triangle and side 2-3 of the triangle

      REAL(DOUBLE)                    :: PSI               ! Angle that one of the triangles needs to be rotated to align with
!                                                            the QUAD axis. This will be used in subr TPLT2 to rotate orthotropic
!                                                            mat'l matrices so that ortho mat'l props for all triangles are aligned

      REAL(DOUBLE)                    :: ALPHA_K           ! ALPHA(K)
      REAL(DOUBLE)                    :: AREA              ! Area of all 4 triangles
      REAL(DOUBLE)                    :: AREA_TRIA_(4)     ! Area of a triangle in the quad
      REAL(DOUBLE)                    :: B54(15,12)        ! Matrix to transform 5 node quad to 4 node quad

      REAL(DOUBLE)                    :: B2V_TT(3,9)       ! 3x9 virgin strain-displ matrix for one MIN3 in tria local coords
!                                                            Virgin  DOF order: w1, w2, w3, th-x1, th-x2, th-x3, th-y1, th-y2, th-y3
      REAL(DOUBLE)                    :: B2M_TT(3,9)       ! 3x9 strain-displ matrix for one MIN3 in tria local coords
!                                                            MYSTRAN DOF order: w1, th-x1, th-y1, w2, th-x2, th-y2, w3, th-x3, th-y3
      REAL(DOUBLE)                    :: B2M_TQ(3,9)       ! B2M_TT for a triangle transformed to the local elem system for the quad
      REAL(DOUBLE)                    :: B2M_QQk(3,9)      ! B2M_TQ strain terms transformed to quad stress axes
      REAL(DOUBLE)                    :: B2M_QQ(3,15,4)    ! B2M_TQ for a tria expanded to 12 DOF for the quad
      REAL(DOUBLE)                    :: B2M_QQ_5(3,15,5)  ! 1st one is avg of B2M_QQ for all 4 triangles. 2nd-5th are B2M_QQ above
      REAL(DOUBLE)                    :: B2M_QQ_4(3,12,5)  ! B2M_QQ_5 reduced from 5 to 4 nodes of the quad for the 5 B2M_QQ_5's

      REAL(DOUBLE)                    :: B3V_TT(3,9)       ! 3x9 virgin strain-displ matrix for one MIN3 in tria local coords
!                                                            Virgin  DOF order: w1, w2, w3, th-x1, th-x2, th-x3, th-y1, th-y2, th-y3
      REAL(DOUBLE)                    :: B3M_TT(3,9)       ! 3x9 strain-displ matrix for one MIN3 in tria local coords
!                                                            MYSTRAN DOF order: w1, th-x1, th-y1, w2, th-x2, th-y2, w3, th-x3, th-y3
      REAL(DOUBLE)                    :: B3M_TQ(3,9)       ! B3M_TT for a triangle transformed to the local elem system for the quad
      REAL(DOUBLE)                    :: B3M_QQk(3,9)      ! B3M_TQ strain terms transformed to quad stress axes
      REAL(DOUBLE)                    :: B3M_QQ(3,15,4)    ! B3M_TQ for a tria expanded to 12 DOF for the quad
      REAL(DOUBLE)                    :: B3M_QQ_5(3,15,5)  ! 1st one is avg of B3M_QQ for all 4 triangles. 2nd-5th are B3M_QQ above
      REAL(DOUBLE)                    :: B3M_QQ_4(3,12,5)  ! B3M_QQ_5 reduced from 5 to 4 nodes of the quad for the 5 B3M_QQ_5's

      REAL(DOUBLE)                    :: CTH2              ! COS(2*THETA)
      REAL(DOUBLE)                    :: D(9,9)            ! Matrix to transform a tria stiff matrix from tria to quad local axes
      REAL(DOUBLE)                    :: DEN               ! An intermediate variable
      REAL(DOUBLE)                    :: DT_QUAD(5,NTSUB)  ! Array DT (elem temps/gradient) for the quad so we can send appropriate
!                                                            DT array to subr TPLT2 for each of the 4 triangles
      REAL(DOUBLE)                    :: DUM               ! Intermediate constant
      REAL(DOUBLE)                    :: DUM1(9,9)         ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM2(15,12)       ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM3(3,18,1)      ! Matrix returned from subr TPLT2 not used here (BIG_BB for a tria but we
!                                                            generate BIG_BB for the quad here from BV_TT (same as all other mats)
      REAL(DOUBLE)                    :: FCONV3_AVG        ! Average of the FCONV(3) values from the 4 triangles

      REAL(DOUBLE)                    :: KV_TT(9,9)        ! 9x9 virgin stiffness matrix for one MIN3 in tria local coords
!                                                            Virgin  DOF order: w1, w2, w3, th-x1, th-x2, th-x3, th-y1, th-y2, th-y3
      REAL(DOUBLE)                    :: KM_TT(9,9)        ! 9x9  stiffness matrix for one MIN3 in tria local coords
!                                                            MYSTRAN DOF order: w1, th-x1, th-y1, w2, th-x2, th-y2, w3, th-x3, th-y3
      REAL(DOUBLE)                    :: KM_TQ(9,9)        ! KM_TT for a triangle transformed to the local elem system for the quad
      REAL(DOUBLE)                    :: KM_QQ(4,15,15)    ! KM_TQ for a tria expanded to 15x15 size for the 5 nodes of the quad
      REAL(DOUBLE)                    :: KM_QQ_5(15,15)    ! Sum all KM_QQi. Stiff matrix for the 5 node quad in local quad coords
      REAL(DOUBLE)                    :: KM_QQ_4(12,12)    ! KM_QQ_5 reduced to 4 nodes

      REAL(DOUBLE)                    :: PHI_SQ_TRIA(4)    ! PHI_SQ for each TRIA3 from subr TRPLT2 called herein
      REAL(DOUBLE)                    :: L12               ! Magnitude of vector V12 (length of side 1-2 of sub-triangle)

      REAL(DOUBLE)                    :: L13               ! Magnitude of vector V13 (length of side 1-3 of sub-triangle)
      REAL(DOUBLE)                    :: NUM               ! An intermediate variable

      REAL(DOUBLE)                    :: PPV_TT(9,NSUB)    ! 9xNSUB virgin thermal load matrix for one MIN3 in tria local coords
!                                                            Virgin  DOF order: w1, w2, w3, th-x1, th-x2, th-x3, th-y1, th-y2, th-y3
      REAL(DOUBLE)                    :: PPM_TT(9,NSUB)    ! 9xNSUB pressure load matrix for one MIN3 in tria local coords
!                                                            MYSTRAN DOF order: w1, th-x1, th-y1, w2, th-x2, th-y2, w3, th-x3, th-y3
      REAL(DOUBLE)                    :: PPM_TQ(9,NSUB)    ! PPM_TT for a triangle transformed to the local elem system for the quad
      REAL(DOUBLE)                    :: PPM_QQ(4,15,NSUB) ! PPM_TQ for a tria expanded to 15xNSUB size for the 5 nodes of the quad
      REAL(DOUBLE)                    :: PPM_QQ_5(15,NSUB) ! Sum all PPM_QQi. Press loads for the 5 node quad in local quad coords
      REAL(DOUBLE)                    :: PPM_QQ_4(12,NSUB) ! PPM_QQ_5 reduced to 4 nodes

      REAL(DOUBLE)                    :: PTV_TT(9,NSUB)    ! 9xNTSUB virgin thermal load matrix for one MIN3 in tria local coords
!                                                            Virgin  DOF order: w1, w2, w3, th-x1, th-x2, th-x3, th-y1, th-y2, th-y3
      REAL(DOUBLE)                    :: PTM_TT(9,NSUB)    ! 9xNTSUB thermal load matrix for one MIN3 in tria local coords
!                                                            MYSTRAN DOF order: w1, th-x1, th-y1, w2, th-x2, th-y2, w3, th-x3, th-y3
      REAL(DOUBLE)                    :: PTM_TQ(9,NTSUB)   ! PTM_TT for a triangle transformed to the local elem system for the quad
      REAL(DOUBLE)                    :: PTM_QQ(4,15,NTSUB)! PTM_TQ for a tria expanded to 15xNTSUB size for the 5 nodes of the quad
      REAL(DOUBLE)                    :: PTM_QQ_5(15,NTSUB)! Sum all PTM_QQi. Therm loads for the 5 node quad in local quad coords
      REAL(DOUBLE)                    :: PTM_QQ_4(12,NTSUB)! PTM_QQ_5 reduced to 4 nodes

      REAL(DOUBLE)                    :: RAT               ! An intermediate variable

      REAL(DOUBLE)                    :: S2V_TT(3,9)       ! 3x9 virgin stress recovery matrix for one MIN3 in tria local coords
!                                                            Virgin  DOF order: w1, w2, w3, th-x1, th-x2, th-x3, th-y1, th-y2, th-y3
      REAL(DOUBLE)                    :: S2M_TT(3,9)       ! 3x9 stress recovery matrix for one MIN3 in tria local coords
!                                                            MYSTRAN DOF order: w1, th-x1, th-y1, w2, th-x2, th-y2, w3, th-x3, th-y3
      REAL(DOUBLE)                    :: S2M_TQ(3,9)       ! S2M_TT for a triangle transformed to the local elem system for the quad
      REAL(DOUBLE)                    :: S2M_QQk(3,9)      ! S2M_TQ stress terms transformed to quad stress axes
      REAL(DOUBLE)                    :: S2M_QQ(3,15,4)    ! S2M_TQ for each of the tria's expanded to the 15 DOF for the quad
      REAL(DOUBLE)                    :: S2M_QQ_5(3,15,5)  ! 1st one is avg of S2M_QQ for all 4 triangles. 2nd-5th are S2M_QQ above
      REAL(DOUBLE)                    :: S2M_QQ_4(3,12,5)  ! S2M_QQ_5 reduced from 5 to 4 nodes of the quad for the 5 S2M_QQ_5's

      REAL(DOUBLE)                    :: S3V_TT(3,9)       ! 3x9 virgin stress recovery matrix for one MIN3 in tria local coords
!                                                            Virgin  DOF order: w1, w2, w3, th-x1, th-x2, th-x3, th-y1, th-y2, th-y3
      REAL(DOUBLE)                    :: S3M_TT(3,9)       ! 3x9 stress recovery matrix for one MIN3 in tria local coords
!                                                            MYSTRAN DOF order: w1, th-x1, th-y1, w2, th-x2, th-y2, w3, th-x3, th-y3
      REAL(DOUBLE)                    :: S3M_TQ(3,9)       ! S3M_TT for a triangle transformed to the local elem system for the quad
      REAL(DOUBLE)                    :: S3M_QQk(3,9)      ! S3M_TQ stress terms transformed to quad stress axes
      REAL(DOUBLE)                    :: S3M_QQ(3,15,4)    ! S3M_TQ for a tria expanded to 12 DOF for the quad
      REAL(DOUBLE)                    :: S3M_QQ_5(3,15,5)  ! 1st one is avg of S3M_QQ for all 4 triangles. 2nd-5th are S3M_QQ above
      REAL(DOUBLE)                    :: S3M_QQ_4(3,12,5)  ! S3M_QQ_5 reduced from 5 to 4 nodes of the quad for the 5 S3M_QQ_5's

      REAL(DOUBLE)                    :: STH2              ! SIN(2*THETA)
      REAL(DOUBLE)                    :: SUM_DIAGS         ! Sum of diag terms from KOO*KOOI for static cond to elim node 5
      REAL(DOUBLE)                    :: SUM_OFF_DIAGS     ! Sum of off-diag terms from KOO*KOOI for static cond to elim node 5
      REAL(DOUBLE)                    :: TEMP_NODE_5(NTSUB)! Calculated temperature of the 5th node of the quad
      REAL(DOUBLE)                    :: TE_SHEAR(3,3)     ! Coord transf for transverse stresses from tria axes to quad axes
      REAL(DOUBLE)                    :: TE_STRESS(3,3)    ! Coord transf for bend &in-plane shear stresses from tria to quad axes
      REAL(DOUBLE)                    :: TE_TRIA(3,3)      ! Coord transf between tria basic and local coords: UEL = TE_TRIA*UEB
      REAL(DOUBLE)                    :: TET(3,3)          ! Transpose of coord transformation matrix TE
      REAL(DOUBLE)                    :: THETA             ! Angle (rad) between side 1-2 of sub-triangle and quad x axis

                                                           ! 1st 4 rows of XQB are equal to XEB
      REAL(DOUBLE)                    :: XQB(5,3)          ! X, Y coords of the 5 quad points in basic coords
      REAL(DOUBLE)                    :: XQB5(3)           ! Row 5 of XQB (equal to basic coords of 5th node)

                                                           ! 1st 4 rows of XQL are equal to XEL
      REAL(DOUBLE)                    :: XQL(5,3)          ! X, Y coords of the 5 quad points in quad local coords
      REAL(DOUBLE)                    :: XQL5(3)           ! Row 5 of XQL (equal to local coords of 5th node)

      REAL(DOUBLE)                    :: X2QL, Y2QL        ! X, Y coords of element node 2 in quad local elem coord system
      REAL(DOUBLE)                    :: X3QL, Y3QL        ! X, Y coords of element node 3 in quad local elem coord system
      REAL(DOUBLE)                    :: X4QL, Y4QL        ! X, Y coords of element node 4 in quad local elem coord system
      REAL(DOUBLE)                    :: X5QL, Y5QL        ! X, Y coords of element node 5 in quad local elem coord system
      REAL(DOUBLE)                    :: X2TL              ! x coord of elem node 2 of one of the 4 triangles (in tria local coords)
      REAL(DOUBLE)                    :: X3TL              ! x coord of elem node 3 of one of the 4 triangles (in tria local coords)
      REAL(DOUBLE)                    :: Y3TL              ! y coord of elem node 3 of one of the 4 triangles (in tria local coords)

      REAL(DOUBLE)                    :: V12(2)            ! Components of a vector along side 1-2 of a sub-triangle
      REAL(DOUBLE)                    :: V13(2)            ! Components of a vector along side 1-3 of a sub-triangle



! **********************************************************************************************************************************
! Initialize

      CALL QPLT3_INIT

! Set OPT_MIN4T to OPT. Then set OPT_MIN4T(4) depending on the MIN4T reduction technIque. This is done since the stiffness
! partitions are needed to reduce the load matrices if MIN4TRED = 'STC' (static condensation for reducing 5 node quad to 4 nodes)

      DO I=1,5
         OPT_MIN4T(I) = OPT(I)
      ENDDO

      IF (MIN4TRED == 'STC') THEN
         OPT_MIN4T(4) = 'Y'
      ENDIF

! Set the ID arrays that correlate triangle DOF with quad DOF

      IDI(1,1) =  1                                        ! means tria 1 elem DOF 1 is quad elem DOF  1
      IDI(1,2) =  2                                        ! means tria 1 elem DOF 2 is quad elem DOF  2
      IDI(1,3) =  3                                        ! means tria 1 elem DOF 3 is quad elem DOF  3
      IDI(1,4) =  4                                        ! means tria 1 elem DOF 4 is quad elem DOF  4
      IDI(1,5) =  5                                        ! means tria 1 elem DOF 5 is quad elem DOF  5
      IDI(1,6) =  6                                        ! means tria 1 elem DOF 6 is quad elem DOF  6
      IDI(1,7) = 13                                        ! means tria 1 elem DOF 7 is quad elem DOF 13
      IDI(1,8) = 14                                        ! means tria 1 elem DOF 8 is quad elem DOF 14
      IDI(1,9) = 15                                        ! means tria 1 elem DOF 9 is quad elem DOF 15

      IDI(2,1) =  4                                        ! means tria 2 elem DOF 1 is quad elem DOF  4
      IDI(2,2) =  5                                        ! means tria 2 elem DOF 2 is quad elem DOF  5
      IDI(2,3) =  6                                        ! means tria 2 elem DOF 3 is quad elem DOF  6
      IDI(2,4) =  7                                        ! means tria 2 elem DOF 4 is quad elem DOF  7
      IDI(2,5) =  8                                        ! means tria 2 elem DOF 5 is quad elem DOF  8
      IDI(2,6) =  9                                        ! means tria 2 elem DOF 6 is quad elem DOF  9
      IDI(2,7) = 13                                        ! means tria 2 elem DOF 7 is quad elem DOF 13
      IDI(2,8) = 14                                        ! means tria 2 elem DOF 8 is quad elem DOF 14
      IDI(2,9) = 15                                        ! means tria 2 elem DOF 9 is quad elem DOF 15

      IDI(3,1) =  7                                        ! means tria 3 elem DOF 1 is quad elem DOF  7
      IDI(3,2) =  8                                        ! means tria 3 elem DOF 2 is quad elem DOF  8
      IDI(3,3) =  9                                        ! means tria 3 elem DOF 3 is quad elem DOF  9
      IDI(3,4) = 10                                        ! means tria 3 elem DOF 4 is quad elem DOF 10
      IDI(3,5) = 11                                        ! means tria 3 elem DOF 5 is quad elem DOF 11
      IDI(3,6) = 12                                        ! means tria 3 elem DOF 6 is quad elem DOF 12
      IDI(3,7) = 13                                        ! means tria 3 elem DOF 7 is quad elem DOF 13
      IDI(3,8) = 14                                        ! means tria 3 elem DOF 8 is quad elem DOF 14
      IDI(3,9) = 15                                        ! means tria 3 elem DOF 9 is quad elem DOF 15

      IDI(4,1) = 10                                        ! means tria 4 elem DOF 1 is quad elem DOF 10
      IDI(4,2) = 11                                        ! means tria 4 elem DOF 2 is quad elem DOF 11
      IDI(4,3) = 12                                        ! means tria 4 elem DOF 3 is quad elem DOF 12
      IDI(4,4) =  1                                        ! means tria 4 elem DOF 4 is quad elem DOF  1
      IDI(4,5) =  2                                        ! means tria 4 elem DOF 5 is quad elem DOF  2
      IDI(4,6) =  3                                        ! means tria 4 elem DOF 6 is quad elem DOF  3
      IDI(4,7) = 13                                        ! means tria 4 elem DOF 7 is quad elem DOF 13
      IDI(4,8) = 14                                        ! means tria 4 elem DOF 8 is quad elem DOF 14
      IDI(4,9) = 15                                        ! means tria 4 elem DOF 9 is quad elem DOF 15

! Set DT_QUAD array here so we can mod the DT array before calling subr TPLT2. Later we will reset DT to DT_QUAD so it is the same
! when we leave this subr as it was when it began

      DO I=1,5
         DO J=1,NTSUB
            DT_QUAD(I,J) = DT(I,J)
         ENDDO
      ENDDO

! Calc temperature of node 5. Use avg of the diagonal averages

      DO J=1,NTSUB
         TEMP_NODE_5(J) = HALF*(HALF*(DT(1,J) + DT(3,J)) + HALF*(DT(2,J) + DT(4,J)))
      ENDDO

! We need the coords of the intersection of the 2 diagonals. This will be the central node to which the 4 non-overlapping
! triangles will all connect. In local coords call these X5QL, Y5QL. First get coords of points 2-4 (pt 1 is origin)

      X2QL = XEL(2,1)   ;   Y2QL = XEL(2,2)
      X3QL = XEL(3,1)   ;   Y3QL = XEL(3,2)
      X4QL = XEL(4,1)   ;   Y4QL = XEL(4,2)

      NUM  = Y4QL*(X2QL - X4QL) - X4QL*(Y2QL - Y4QL)
      DEN  = Y3QL*(X2QL - X4QL) - X3QL*(Y2QL - Y4QL)
      RAT  = NUM/DEN
      X5QL = X3QL*RAT
      Y5QL = Y3QL*RAT
      XQL(1,1) = ZERO   ;   XQL(1,2) = ZERO   ;   XQL(1,3) = ZERO
      XQL(2,1) = X2QL   ;   XQL(2,2) = Y2QL   ;   XQL(2,3) = ZERO
      XQL(3,1) = X3QL   ;   XQL(3,2) = Y3QL   ;   XQL(3,3) = ZERO
      XQL(4,1) = X4QL   ;   XQL(4,2) = Y4QL   ;   XQL(4,3) = ZERO
      XQL(5,1) = X5QL   ;   XQL(5,2) = Y5QL   ;   XQL(5,3) = ZERO

! Calculate basic coords of the 5th node from the local coords of that node using the transformation:

! Since            XEL     = (XEB - XEB1)*TEt    (where XEB1 are the coords of the 1st grid point of the quad),
! then:            XEL*TE  = (XEB - XEB1)*TEt*TE (where TEt is the transpose of TE.
! however,         TEt*TE  = I                   (I is the identity matrix).
! Thus             XEB     = XEL*TE + XEB1  ,    however, we will use the 1st 4 of 5 rows of XQB to hold the 1st 4 rows of XEB and
! calculate        XQB     = XQL*TE + XEB1

! This can be used to solve for basic coords of the 5th quad point (the point near the center of the quad where the 4 subtri's meet.
! The 1st 4 rows of the above XEB will be the same as XEB and the 5th row are the coords of node 5

      DO I=1,4                                             ! XQB: basic coords of all 5 quad nodes. XEB:basic coords of corner nodes
         DO J=1,3
            XQB(I,J) = XEB(I,J)
         ENDDO
      ENDDO

      DO I=1,3                                             ! Transpose TE to TEt
         DO J=1,3
            TET(I,J) = TE(J,I)
         ENDDO
      ENDDO

      DO I=1,3                                             ! XQL5 are the local coords of node 5 (5th row of XQL)
         XQL5(I) = XQL(5,I)
      ENDDO

      CALL MATMULT_FFF ( XQL5, TE, 1, 3, 3, XQB5 )         ! Mult XQL*TEt
      DO J=1,3                                             ! and add basic coords of node 1
         XQB5(J) = XQB5(J) + XEB(1,J)                      ! to get basic coords of node 5
      ENDDO

      DO J=1,3                                             ! Put 5th node basic coords into array XQB
         XQB(5,J) = XQB5(J)
      ENDDO


      IF (DEBUG(198) > 0) CALL DEB_QPLT3_1 ( OPT_MIN4T, 0 )

! **********************************************************************************************************************************
      FCONV3_AVG = ZERO

! Get matrices for the triangles (I = 1, 2, 3, 4)

      DO K=1,NUM_TRIAS
         ALPHA(K) = ZERO
         BETA(K)  = ZERO
         GAMMA(K) = ZERO
      ENDDO

      PSI = ZERO

trias:DO K=1,NUM_TRIAS

         TRIA_NUM = K
         CALL ELMGM_TRIA ( TRIA_NUM, XQB, XQL, V12, V13, L12, L13, ALPHA, BETA, GAMMA, PSI, THETA, X2TL, X3TL, Y3TL, TE_TRIA )
         IF (PROG_ERR > 0) THEN
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         AREA_TRIA_(K) = X2TL*Y3TL/TWO
         ALPHA_K = ALPHA(K)

         MIN4T_QUAD4_TRIA_NO = K
         CALL TPLT2 ( OPT_MIN4T, AREA_TRIA_(K), X2TL, X3TL, Y3TL, 'N', IERROR, KV_TT, PTV_TT, PPV_TT, B2V_TT, B3V_TT, S2V_TT,      &
                      S3V_TT, DUM3, MN4T_QD, TRIA_NUM, PSI )

         PHI_SQ_TRIA(K) = PHI_SQ

         FCONV3_AVG = FCONV3_AVG + FCONV(3)/NUM_TRIAS
         IF (IERROR > 0) THEN
            WRITE(ERR,99999) TRIA_NUM
            WRITE(F06,99999) TRIA_NUM
         ENDIF

         DO I=1,9
            DO J=1,9
               D(I,J) = ZERO
            ENDDO
         ENDDO
                                                           ! Matrix to transform a tria stiff matrix from tria to quad local axes
         D(1,1) = TE_TRIA(1,1)   ;   D(1,2) = TE_TRIA(1,2)   ;   D(1,3) = TE_TRIA(1,3)
         D(2,1) = TE_TRIA(2,1)   ;   D(2,2) = TE_TRIA(2,2)   ;   D(2,3) = TE_TRIA(2,3)
         D(3,1) = TE_TRIA(3,1)   ;   D(3,2) = TE_TRIA(3,2)   ;   D(3,3) = TE_TRIA(3,3)

         D(4,4) = TE_TRIA(1,1)   ;   D(4,5) = TE_TRIA(1,2)   ;   D(4,6) = TE_TRIA(1,3)
         D(5,4) = TE_TRIA(2,1)   ;   D(5,5) = TE_TRIA(2,2)   ;   D(5,6) = TE_TRIA(2,3)
         D(6,4) = TE_TRIA(3,1)   ;   D(6,5) = TE_TRIA(3,2)   ;   D(6,6) = TE_TRIA(3,3)

         D(7,7) = TE_TRIA(1,1)   ;   D(7,8) = TE_TRIA(1,2)   ;   D(7,9) = TE_TRIA(1,3)
         D(8,7) = TE_TRIA(2,1)   ;   D(8,8) = TE_TRIA(2,2)   ;   D(8,9) = TE_TRIA(2,3)
         D(9,7) = TE_TRIA(3,1)   ;   D(9,8) = TE_TRIA(3,2)   ;   D(9,9) = TE_TRIA(3,3)

         DO I=1,3
            DO J=1,3
               TE_STRESS(I,J) = ZERO
            ENDDO
         ENDDO

         CTH2 = DCOS(TWO*THETA)                            ! Calc matrix which will transform stresses thru angle THETA
         STH2 = DSIN(TWO*THETA)
         TE_STRESS(1,1) = HALF*(ONE + CTH2)   ;   TE_STRESS(1,2) =  HALF*(ONE - CTH2)   ;   TE_STRESS(1,3) =  STH2
         TE_STRESS(2,1) = HALF*(ONE - CTH2)   ;   TE_STRESS(2,2) =  HALF*(ONE + CTH2)   ;   TE_STRESS(2,3) = -STH2
         TE_STRESS(3,1) = HALF*STH2           ;   TE_STRESS(3,2) = -HALF*STH2           ;   TE_STRESS(3,3) =  CTH2

! **********************************************************************************************************************************
! Expand/transform strain-displ matrices from virgin tria 9 DOF's in tria coords to quad 15 DOF's in quad coords
! Do this indep of OPT since some of the the strain displ matrices are needed for BIG_BB

         DO I=1,3                                       ! Reorder B2M from tria virgin DOF to local tria DOF order
            DO J=1,9
               B2M_TT(I,IDV(J)) = B2V_TT(I,J)
            ENDDO
         ENDDO
                                                        ! Transform B2M to quad local axes
         CALL MATMULT_FFF ( B2M_TT, D, 3, 9, 9, B2M_TQ )
         CALL MATMULT_FFF ( TE_STRESS, B2M_TQ, 3, 3, 9, B2M_QQk )

         DO I=1,3                                       ! Expand B2M up to 15 DOF for the 5 node quad
            DO J=1,9
               B2M_QQ(I,IDI(K,J),K) = B2M_QQk(I,J)
            ENDDO
         ENDDO

         DO I=1,3                                       ! Reorder B3M from tria virgin DOF to local tria DOF order
            DO J=1,9
               B3M_TT(I,IDV(J)) = B3V_TT(I,J)
            ENDDO
         ENDDO
                                                        ! Transform B3M to quad local axes
         CALL MATMULT_FFF ( B3M_TT, D, 3, 9, 9, B3M_TQ )

         DO I=1,3                                       ! Expand B3M up to 15 DOF for the 5 node quad
            DO J=1,3
               TE_SHEAR(I,J) = ZERO
            ENDDO
         ENDDO
         TE_SHEAR(1,1) = TE_TRIA(2,2)   ;   TE_SHEAR(1,2) = TE_TRIA(3,2)
         TE_SHEAR(2,1) = TE_TRIA(2,3)   ;   TE_SHEAR(2,2) = TE_TRIA(3,3)
         TE_SHEAR(3,3) = ONE
         CALL MATMULT_FFF ( TE_SHEAR, B3M_TQ, 3, 3, 9, B3M_QQk )

         DO I=1,3                                       ! Expand to 15 DOF quad (5 nodes/3 DOF each)
            DO J=1,9
               B3M_QQ(I,IDI(K,J),K) = B3M_QQk(I,J)
            ENDDO
         ENDDO







! Expand/transform thermal loads from virgin tria 9 DOF's in tria coords to quad 15 DOF's in quad coords

         IF (OPT_MIN4T(2) == 'Y') THEN

            DO I=1,9                                       ! Reorder from tria virgin DOF to local tria DOF order
               DO J=1,NTSUB
                  PTM_TT(IDV(I),J) = PTV_TT(I,J)
               ENDDO
            ENDDO
                                                           ! Transform to quad local axes
            CALL MATMULT_FFF_T ( D, PTM_TT, 9, 9, NTSUB, PTM_TQ )

            DO I=1,9                                       ! Expand up to 15 DOF for the 5 node quad
               DO J=1,NTSUB
                  PTM_QQ(K,IDI(K,I),J) = PTM_TQ(I,J)
               ENDDO
            ENDDO

         ENDIF

! Expand/transform stress recovery matrices (S2, S3) from virgin tria 9 DOF's in tria coords to quad 15 DOF's in quad coords

         IF (OPT_MIN4T(3) == 'Y') THEN

            DO I=1,3                                       ! Reorder S2M from tria virgin DOF to local tria DOF order
               DO J=1,9
                  S2M_TT(I,IDV(J)) = S2V_TT(I,J)
               ENDDO
            ENDDO
                                                           ! Transform S2M to quad local axes
            CALL MATMULT_FFF ( S2M_TT, D, 3, 9, 9, S2M_TQ )
            CALL MATMULT_FFF ( TE_STRESS, S2M_TQ, 3, 3, 9, S2M_QQk )

            DO I=1,3                                       ! Expand S2M up to 15 DOF for the 5 node quad
               DO J=1,9
                  S2M_QQ(I,IDI(K,J),K) = S2M_QQk(I,J)
               ENDDO
            ENDDO

            DO I=1,3                                       ! Reorder S3M from tria virgin DOF to local tria DOF order
               DO J=1,9
                  S3M_TT(I,IDV(J)) = S3V_TT(I,J)
               ENDDO
            ENDDO
                                                           ! Transform S3M to quad local axes
            CALL MATMULT_FFF ( S3M_TT, D, 3, 9, 9, S3M_TQ )

            DO I=1,3                                       ! Expand S3M up to 15 DOF for the 5 node quad
               DO J=1,3
                  TE_SHEAR(I,J) = ZERO
               ENDDO
            ENDDO
            TE_SHEAR(1,1) = TE_TRIA(2,2)   ;   TE_SHEAR(1,2) = TE_TRIA(3,2)
            TE_SHEAR(2,1) = TE_TRIA(2,3)   ;   TE_SHEAR(2,2) = TE_TRIA(3,3)
            TE_SHEAR(3,3) = ONE
            CALL MATMULT_FFF ( TE_SHEAR, S3M_TQ, 3, 3, 9, S3M_QQk )

            DO I=1,3                                       ! Expand to 15 DOF quad (5 nodes/3 DOF each)
               DO J=1,9
                  S3M_QQ(I,IDI(K,J),K) = S3M_QQk(I,J)
               ENDDO
            ENDDO

         ENDIF

! Expand/transform stiffness matrices from virgin tria 9 DOF's in tria coords to quad 15 DOF's in quad coords

         IF (OPT_MIN4T(4) == 'Y') THEN

            DO I=1,9                                       ! Reorder from tria virgin DOF to local tria DOF order
               DO J=1,9
                  KM_TT(IDV(I),IDV(J)) = KV_TT(I,J)
               ENDDO
            ENDDO
                                                           ! Transform to quad local axes
            CALL MATMULT_FFF   ( KM_TT, D , 9, 9, 9, DUM1 )
            CALL MATMULT_FFF_T ( D, DUM1, 9, 9, 9, KM_TQ )

            DO I=1,9                                       ! Expand to 15 DOF quad (5 nodes/3 DOF each)
               DO J=1,9
                  KM_QQ(K,IDI(K,I),IDI(K,J)) = KM_TQ(I,J)
               ENDDO
            ENDDO

         ENDIF

! Expand/transform pressure loads from virgin tria 9 DOF's in tria coords to quad 15 DOF's in quad coords

         IF (OPT_MIN4T(5) == 'Y') THEN

            DO I=1,9                                       ! Reorder from tria virgin DOF to local tria DOF order
               DO J=1,NSUB
                  PPM_TT(IDV(I),J) = PPV_TT(I,J)
               ENDDO
            ENDDO
                                                           ! Transform to quad local axes
            CALL MATMULT_FFF_T ( D, PPM_TT, 9, 9, NSUB, PPM_TQ )

            DO I=1,9                                       ! Expand to 15 DOF quad (5 nodes/3 DOF each)
               DO J=1,NSUB
                  PPM_QQ(K,IDI(K,I),J) = PPM_TQ(I,J)
               ENDDO
            ENDDO

         ENDIF

         IF (DEBUG(198) > 0) CALL DEB_QPLT3_1 ( OPT_MIN4T, K )

      ENDDO trias

! Now calc PHI_SQ for the QUAD based on the TRIA3 values

      PHI_SQ = ZERO
      DUM    = ZERO

      IF (DEBUG(172) == 0) THEN
         DO I=1,NUM_TRIAS                                  ! Use simple average of the 4 TRIA3 PHI_SQ values
            DUM = DUM + PHI_SQ_TRIA(I)
         ENDDO
         PHI_SQ = DUM/NUM_TRIAS
      ELSE
         DO I=1,NUM_TRIAS                                  ! Use area weighted average of the 4 TRIA3 PHI_SQ values
            DUM = DUM + PHI_SQ_TRIA(I)*AREA_TRIA_(I)
         ENDDO
         PHI_SQ = DUM/AREA_QUAD
      ENDIF


! **********************************************************************************************************************************
! Add matrices for all triangles. They are all in quad coords for 3 DOF/5 nodes per triangle

      AREA = ZERO
      DO I=1,NUM_TRIAS
         AREA = AREA + AREA_TRIA_(I)
      ENDDO
      IF (AREA /= AREA_QUAD) THEN
      ENDIF

      DO I=1,3                                             ! B2M_QQ_5(i,j,1) is average of all 4 triangles
         DO J=1,15
            B2M_QQ_5(I,J,1) = ZERO
            DO K=1,NUM_TRIAS
               B2M_QQ_5(I,J,1) = B2M_QQ_5(I,J,1) + QUARTER* B2M_QQ(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      DO I=1,3                                             ! B2M_QQ_5(i,j,k) for k=2,3,4,5 are indiv triangle S2M_QQ(i,j,k)
         DO J=1,15
            DO K=2,5
               B2M_QQ_5(I,J,K) = B2M_QQ(I,J,K-1)
            ENDDO
         ENDDO
      ENDDO

      DO I=1,3                                             ! B3M_QQ_5(i,j,1) is average of all 4 triangles
         DO J=1,15
            B3M_QQ_5(I,J,1) = ZERO
            DO K=1,NUM_TRIAS
               B3M_QQ_5(I,J,1) = B3M_QQ_5(I,J,1) + QUARTER*B3M_QQ(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      DO I=1,3                                             ! B3M_QQ_5(i,j,k) for k=2,3,4,5 are indiv triangle S3M_QQ(i,j,k)
         DO J=1,15
            DO K=2,5
               B3M_QQ_5(I,J,K) = B3M_QQ(I,J,K-1)
            ENDDO
         ENDDO
      ENDDO

      IF (OPT_MIN4T(2) == 'Y') THEN                        ! Add thermal load matrices for all 4 triangles

         DO I=1,15
            DO J=1,NTSUB
               PTM_QQ_5(I,J) = ZERO
               DO K=1,NUM_TRIAS
                  PTM_QQ_5(I,J) = PTM_QQ_5(I,J) + PTM_QQ(K,I,J)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      IF (OPT_MIN4T(3) == 'Y') THEN                        ! Stress recovery matrices. 5 per quad with 1st being avg of the tria's
!                                                            to get avg stresses and the other 4 are for each indiv triangle

         DO I=1,3                                          ! S2M_QQ_5(i,j,1) is average of all 4 triangles
            DO J=1,15
               S2M_QQ_5(I,J,1) = ZERO
               DO K=1,NUM_TRIAS
                  S2M_QQ_5(I,J,1) = S2M_QQ_5(I,J,1) + QUARTER*S2M_QQ(I,J,K)
               ENDDO
            ENDDO
         ENDDO
         DO I=1,3                                          ! S2M_QQ_5(i,j,k) for k=2,3,4,5 are indiv triangle S2M_QQ(i,j,k)
            DO J=1,15
               DO K=2,5
                  S2M_QQ_5(I,J,K) = S2M_QQ(I,J,K-1)
               ENDDO
            ENDDO
         ENDDO

         DO I=1,3                                          ! S3M_QQ_5(i,j,1) is average of all 4 triangles
            DO J=1,15
               S3M_QQ_5(I,J,1) = ZERO
               DO K=1,NUM_TRIAS
                  S3M_QQ_5(I,J,1) = S3M_QQ_5(I,J,1) + QUARTER*S3M_QQ(I,J,K)
               ENDDO
            ENDDO
         ENDDO
         DO I=1,3                                          ! S3M_QQ_5(i,j,k) for k=2,3,4,5 are indiv triangle S3M_QQ(i,j,k)
            DO J=1,15
               DO K=2,5
                  S3M_QQ_5(I,J,K) = S3M_QQ(I,J,K-1)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      IF (OPT_MIN4T(4) == 'Y') THEN                        ! Add stiffness matrices for all 4 triangles

         DO I=1,15
            DO J=1,15
               KM_QQ_5(I,J) = ZERO
               DO K=1,NUM_TRIAS
                  KM_QQ_5(I,J) = KM_QQ_5(I,J) + KM_QQ(K,I,J)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      IF (OPT_MIN4T(5) == 'Y') THEN                        ! Add pressure load matrices for all 4 triangles

         DO I=1,15
            DO J=1,NSUB
               PPM_QQ_5(I,J) = ZERO
               DO K=1,NUM_TRIAS
                  PPM_QQ_5(I,J) = PPM_QQ_5(I,J) + PPM_QQ(K,I,J)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      IF (DEBUG(198) > 0) CALL DEB_QPLT3_1 ( OPT_MIN4T, 5 )

! **********************************************************************************************************************************
! Reduce 5 node quad to 4 node quad

      IF      (MIN4TRED == 'B54') THEN

         CALL B54_REDUCTION ( XQL, IERROR )
         IF (IERROR > 0) THEN
            WRITE(ERR,99999) TRIA_NUM
            WRITE(F06,99999) TRIA_NUM
            RETURN
         ENDIF

      ELSE IF (MIN4TRED == 'STC') THEN

         CALL STATIC_CONDENSATION

      ELSE

         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1929) MIN4TRED
            WRITE(F06,1929) MIN4TRED
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1929
            ENDIF
         ENDIF

      ENDIF

      IF (DEBUG(198) > 0) CALL DEB_QPLT3_1 ( OPT_MIN4T, 7 )

! **********************************************************************************************************************************
! Expand 12 DOF (4 grids with w, theta-x, theta-y DOF's) quad to MYSTRAN 4 grids @ 6 DOF each.

      DO I=1,3
         DO J=1,12
            DO K=1,4
               BIG_BB(I,IDM(J),K) = B2M_QQ_4(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO I=1,3
         DO J=1,12
            DO K=1,5
               BE2(I,IDM(J),K) = B2M_QQ_4(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      DO I=1,3
         DO J=1,12
            DO K=1,5
               BE3(I,IDM(J),K) = B3M_QQ_4(I,J,K)
            ENDDO
         ENDDO
      ENDDO

      IF (OPT_MIN4T(2) == 'Y') THEN

         DO I=1,12
            DO J=1,NTSUB
               PTE(IDM(I),J) = PTM_QQ_4(I,J)
            ENDDO
         ENDDO

      ENDIF

      IF (OPT_MIN4T(3) == 'Y') THEN

         DO I=1,3
            DO J=1,12
               DO K=1,5
                  SE2(I,IDM(J),K) = S2M_QQ_4(I,J,K)
               ENDDO
            ENDDO
         ENDDO

         DO I=1,3
            DO J=1,12
               DO K=1,5
                  SE3(I,IDM(J),K) = S3M_QQ_4(I,J,K)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      IF (OPT_MIN4T(4) == 'Y') THEN

         DO I=1,12
            DO J=1,12
               KE(IDM(I),IDM(J)) = KE(IDM(I),IDM(J)) + KM_QQ_4(I,J)
            ENDDO
         ENDDO

         FCONV(3) = FCONV3_AVG

      ENDIF

      IF (OPT_MIN4T(5) == 'Y') THEN

         DO I=1,12
            DO J=1,NSUB
               PPE(IDM(I),J) = PPM_QQ_4(I,J)
            ENDDO
         ENDDO

      ENDIF

      IF (DEBUG(198) > 0) CALL DEB_QPLT3_1 ( OPT_MIN4T, 8 )

! **********************************************************************************************************************************
! Reset DT to DT_QUAD

      DO I=1,5
         DO J=1,NTSUB
            DT(I,J) = DT_QUAD(I,J)
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************
 1929 FORMAT(' *ERROR  1929: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CHAR PARAMETER MIN4TRED MUST BE EITHER "B54" OR "STC" BUT IS "',A,'"')

99999 FORMAT(14X,' THIS ERROR OCCURRED IN SUB-TRIANGLE ',I2,' OF THE QUAD ELEMENT')

! **********************************************************************************************************************************

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE QPLT3_INIT

      IMPLICIT NONE

      DO I=1,3
         DO J=1,ELDOF
            DO K=1,4
               BIG_BB(I,J,K) = ZERO
            ENDDO
         ENDDO
      ENDDO

      DO I=1,12
         DO J=1,NTSUB
            PTM_QQ_4(I,J) = ZERO
         ENDDO
         DO J=1,12
            KM_QQ_4(I,J) = ZERO
         ENDDO
         DO J=1,NSUB
            PPM_QQ_4(I,J) = ZERO
         ENDDO
      ENDDO

      DO I=1,3
         DO J=1,12
            DO K=1,5
               B2M_QQ_4(I,J,K) = ZERO
               S2M_QQ_4(I,J,K) = ZERO
            ENDDO
         ENDDO
      ENDDO

      DO I=1,3
         DO J=1,12
            DO K=1,5
               B3M_QQ_4(I,J,K) = ZERO
               S3M_QQ_4(I,J,K) = ZERO
            ENDDO
         ENDDO
      ENDDO

      DO I=1,3
         DO J=1,15
            DO K=1,5
               B2M_QQ_5(I,J,K) = ZERO
               S2M_QQ_5(I,J,K) = ZERO
            ENDDO
         ENDDO
      ENDDO

      DO I=1,3
         DO J=1,15
            DO K=1,5
               B3M_QQ_5(I,J,K) = ZERO
               S3M_QQ_5(I,J,K) = ZERO
            ENDDO
         ENDDO
      ENDDO

      DO I=1,15
         DO J=1,NTSUB
            PTM_QQ_5(I,J) = ZERO
         ENDDO
         DO J=1,15
            KM_QQ_5(I,J) = ZERO
         ENDDO
         DO J=1,NSUB
            PPM_QQ_5(I,J) = ZERO
         ENDDO
      ENDDO

      DO K=1,4
         DO I=1,15
            DO J=1,NTSUB
               PTM_QQ(K,I,J) = ZERO
            ENDDO
            DO J=1,15
               KM_QQ(K,I,J) = ZERO
            ENDDO
            DO J=1,NSUB
               PPM_QQ(K,I,J) = ZERO
            ENDDO
         ENDDO
      ENDDO

      DO K=1,4
         DO I=1,3
            DO J=1,15
               S2M_QQ(I,J,K) = ZERO
               B2M_QQ(I,J,K) = ZERO
            ENDDO
         ENDDO
      ENDDO

      DO K=1,4
         DO I=1,3
            DO J=1,15
               B3M_QQ(I,J,K) = ZERO
               S3M_QQ(I,J,K) = ZERO
            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE QPLT3_INIT

! ##################################################################################################################################

      SUBROUTINE ELMGM_TRIA ( TRIA, XQB, XQL, V12, V13, L12, L13, ALPHA, BETA, GAMMA, PSI, THETA, X2TL, X3TL, Y3TL, TE_TRIA )

! Generates transformation matrix from one of the sub-triangles of the quad to the quad local coord system

!                             | T3-tria |   | 1      0           0      | | T3-quad |
!                             | R1-tria | = | 0  cos(THETA)  sin(THETA) | | R1-quad |
!                             | R2-tria |   | 0 -sin(THETA   cos(THETA) | | R2-quad |


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE IOUNT1, ONLY                :  ERR, F06
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MODEL_STUF, ONLY            :  NUM_EMG_FATAL_ERRS, XTB, XTL

      IMPLICIT NONE

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG), INTENT(IN)       :: TRIA              ! 1, 2, 3, or 4 designator of a subtriangle of the quad

      REAL(DOUBLE), INTENT(IN)        :: XQB(5,3)          ! coords of 5 quad nodes in quad basic coords
      REAL(DOUBLE), INTENT(IN)        :: XQL(5,3)          ! coords of 5 quad nodes in quad local coords
      REAL(DOUBLE), INTENT(OUT)       :: X2TL              ! X coord of this triangle node 2 in this triangle local coords
      REAL(DOUBLE), INTENT(OUT)       :: X3TL              ! X coord of this triangle node 3 in this triangle local coords
      REAL(DOUBLE), INTENT(OUT)       :: Y3TL              ! Y coord of this triangle node 3 in this triangle local coords
      REAL(DOUBLE), INTENT(OUT)       :: TE_TRIA(3,3)      ! Coord transf between tria local and quad coords
      REAL(DOUBLE), INTENT(INOUT)     :: ALPHA(4)          ! Angle (rad) between side 1-2 of a triangle and side 1-3 of the triangle
      REAL(DOUBLE), INTENT(INOUT)     :: BETA(4)           ! Angle (rad) between side 1-2 of a triangle and side 2-3 of the triangle
      REAL(DOUBLE), INTENT(INOUT)     :: GAMMA(4)          ! Angle (rad) between side 1-3 of a triangle and side 2-3 of the triangle
      REAL(DOUBLE), INTENT(OUT)       :: PSI               ! Angle to rotate triangle ortho mat'l props (see declarations in main)
      REAL(DOUBLE)                    :: CALPHA            ! cos(ALPHA)
      REAL(DOUBLE)                    :: CBETA             ! cos(BETA)
      REAL(DOUBLE), INTENT(OUT)       :: THETA             ! Angle (rad) between side 1-2 of triangle and quad x axis
      REAL(DOUBLE)                    :: CTHETA            ! cos(THETA)
      REAL(DOUBLE)                    :: STHETA            ! sin(THETA)
      REAL(DOUBLE), INTENT(OUT)       :: L12               ! Magnitude of vector V12 (length of side 1-2 of sub-triangle)
      REAL(DOUBLE), INTENT(OUT)       :: L13               ! Magnitude of vector V13 (length of side 1-3 of sub-triangle)
      REAL(DOUBLE)                    :: L23               ! Magnitude of vector V23 (length of side 2-3 of sub-triangle)
      REAL(DOUBLE), INTENT(OUT)       :: V12(2)            ! Components of a vector along side 1-2 of a sub-triangle
      REAL(DOUBLE), INTENT(OUT)       :: V13(2)            ! Components of a vector along side 1-3 of a sub-triangle
      REAL(DOUBLE)                    :: V23(2)            ! Components of a vector along side 2-3 of a sub-triangle

! **********************************************************************************************************************************
! Initialize

      DO I=1,3
         DO J=1,3
            TE_TRIA = ZERO
         ENDDO
      ENDDO
      TE_TRIA(1,1) = ONE

! Calculate a vector from end of tria element at node 1 to end of element at node 2

      IF      (TRIA == 1) THEN                             ! This triangle connects to nodes 1-2-5

         V12(1) = XQL(2,1) - XQL(1,1)
         V12(2) = XQL(2,2) - XQL(1,2)

         V13(1) = XQL(5,1) - XQL(1,1)
         V13(2) = XQL(5,2) - XQL(1,2)

         V23(1) = XQL(5,1) - XQL(2,1)
         V23(2) = XQL(5,2) - XQL(2,2)

         DO J=1,NTSUB
            DT(1,J) = DT_QUAD(1,J)                         ! Node 1 of tria 1 is node 1 of the quad
            DT(2,J) = DT_QUAD(2,J)                         ! Node 2 of tria 1 is node 2 of the quad
            DT(3,J) = TEMP_NODE_5(J)                       ! Node 3 of tria 1 is node 5 of the quad
            DT(4,J) = DT_QUAD(5,J)                         ! Gradient is DT(3) for tria and DT(4) from quad
         ENDDO

         DO J=1,3
            XTB(1,J) = XQB(1,J)
            XTB(2,J) = XQB(2,J)
            XTB(3,J) = XQB(5,J)
         ENDDO

      ELSE IF (TRIA == 2) THEN                             ! This triangle connects to nodes 2-3-5

         V12(1) = XQL(3,1) - XQL(2,1)
         V12(2) = XQL(3,2) - XQL(2,2)

         V13(1) = XQL(5,1) - XQL(2,1)
         V13(2) = XQL(5,2) - XQL(2,2)

         V23(1) = XQL(5,1) - XQL(3,1)
         V23(2) = XQL(5,2) - XQL(3,2)

         DO J=1,NTSUB
            DT(1,J) = DT_QUAD(2,J)                         ! Node 1 of tria 2 is node 2 of the quad
            DT(2,J) = DT_QUAD(3,J)                         ! Node 2 of tria 2 is node 3 of the quad
            DT(3,J) = TEMP_NODE_5(J)                       ! Node 3 of tria 2 is node 5 of the quad
            DT(4,J) = DT_QUAD(5,J)                         ! Gradient is DT(3) for tria and DT(4) from quad
         ENDDO

         DO J=1,3
            XTB(1,J) = XQB(2,J)
            XTB(2,J) = XQB(3,J)
            XTB(3,J) = XQB(5,J)
         ENDDO

      ELSE IF (TRIA == 3) THEN                             ! This triangle connects to nodes 3-4-5

         V12(1) = XQL(4,1) - XQL(3,1)
         V12(2) = XQL(4,2) - XQL(3,2)

         V13(1) = XQL(5,1) - XQL(3,1)
         V13(2) = XQL(5,2) - XQL(3,2)

         V23(1) = XQL(5,1) - XQL(4,1)
         V23(2) = XQL(5,2) - XQL(4,2)

         DO J=1,NTSUB
            DT(1,J) = DT_QUAD(3,J)                         ! Node 1 of tria 1 is node 3 of the quad
            DT(2,J) = DT_QUAD(4,J)                         ! Node 2 of tria 1 is node 4 of the quad
            DT(3,J) = TEMP_NODE_5(J)                       ! Node 3 of tria 1 is node 5 of the quad
            DT(4,J) = DT_QUAD(5,J)                         ! Gradient is DT(3) for tria and DT(4) from quad
         ENDDO

         DO J=1,3
            XTB(1,J) = XQB(3,J)
            XTB(2,J) = XQB(4,J)
            XTB(3,J) = XQB(5,J)
         ENDDO

      ELSE IF (TRIA == 4) THEN                             ! This triangle connects to nodes 4-1-5

         V12(1) = XQL(1,1) - XQL(4,1)
         V12(2) = XQL(1,2) - XQL(4,2)

         V13(1) = XQL(5,1) - XQL(4,1)
         V13(2) = XQL(5,2) - XQL(4,2)

         V23(1) = XQL(5,1) - XQL(1,1)
         V23(2) = XQL(5,2) - XQL(1,2)

         DO J=1,NTSUB
            DT(1,J) = DT_QUAD(4,J)                         ! Node 1 of tria 1 is node 4 of the quad
            DT(2,J) = DT_QUAD(1,J)                         ! Node 2 of tria 1 is node 1 of the quad
            DT(3,J) = TEMP_NODE_5(J)                       ! Node 3 of tria 1 is node 5 of the quad
            DT(4,J) = DT_QUAD(5,J)                         ! Gradient is DT(3) for tria and DT(4) from quad
         ENDDO

         DO J=1,3
            XTB(1,J) = XQB(4,J)
            XTB(2,J) = XQB(1,J)
            XTB(3,J) = XQB(5,J)
         ENDDO

      ELSE

         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         PROG_ERR = PROG_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1907) SUBR_NAME, TRIA
            WRITE(F06,1907) SUBR_NAME, TRIA
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1907
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = TRIA
            ENDIF
         ENDIF

      ENDIF

      L23          = DSQRT(V23(1)*V23(1) + V23(2)*V23(2))

      L12          = DSQRT(V12(1)*V12(1) + V12(2)*V12(2))
      CTHETA       = V12(1)/L12
      STHETA       = V12(2)/L12
      THETA        = ATAN2(V12(2),V12(1))
      TE_TRIA(2,2) = CTHETA   ;   TE_TRIA(2,3) = STHETA
      TE_TRIA(3,2) =-STHETA   ;   TE_TRIA(3,3) = CTHETA

      L13          = DSQRT(V13(1)*V13(1) + V13(2)*V13(2))
!!    CALPHA       = (V12(1)*V13(1) + V12(2)*V13(2))/(L12*L13)   ! Replace this with Law of Cosines even though it apparently works
      CALPHA       = (L13*L13 + L12*L12 - L23*L23)/(TWO*L13*L12) ! Calculated using the Law of Cosines
      ALPHA(TRIA)  = ACOS(CALPHA)
      X3TL         = L13*CALPHA
      Y3TL         = L13*SIN(ALPHA(TRIA))
      X2TL         = L12

      CBETA        = (L23*L23 + L12*L12 - L13*L13)/(TWO*L23*L12)
      BETA(TRIA)   = ACOS(CBETA)

      GAMMA(TRIA)  = PI - ALPHA(TRIA) - BETA(TRIA)

      DO I=1,3
         DO J=1,3
            XTL(I,J) = ZERO
         ENDDO
      ENDDO

      XTL(2,1) = X2TL
      XTL(3,1) = X3TL
      XTL(3,2) = Y3TL

! Calculate PSI, angle to rotate ortho mat'l props for this triangle (in subr TPLT2)

      IF      (TRIA == 1) THEN
         PSI = ZERO
      ELSE IF (TRIA == 2) THEN
         PSI = BETA(1) + ALPHA(2)
      ELSE IF (TRIA == 3) THEN
         PSI = BETA(1) + ALPHA(2) + BETA(2) + ALPHA(3)
      ELSE IF (TRIA == 4) THEN
         PSI = BETA(1) + ALPHA(2) + BETA(2) + ALPHA(3) + BETA(3) + ALPHA(4)
      ENDIF





! **********************************************************************************************************************************
 1907 FORMAT(' *ERROR  1907: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VARIABLE TRIANGLE MUST BE EITHER 1, 2, 3, OR 4 BUT IS "',I8,'"')


! **********************************************************************************************************************************

      END SUBROUTINE ELMGM_TRIA

! ##################################################################################################################################

      SUBROUTINE B54_REDUCTION ( XQL, IERROR )

      USE CONSTANTS_1, ONLY           :  ZERO, QUARTER, HALF, ONE
      USE SCONTR, ONLY                :  FATAL_ERR
      USE IOUNT1, ONLY                :  ERR, F06
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF,ONLY             :  NUM_EMG_FATAL_ERRS

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: IERROR            ! Error indicator

      REAL(DOUBLE) , INTENT(IN)       :: XQL(5,2)          ! X, Y coords of the 5 quad points in quad local coords
      REAL(DOUBLE)                    :: A12,A23,A34,A41   ! Coefficients in matrix A_5
      REAL(DOUBLE)                    :: A15,A25,A35,A45   ! Coefficients in matrix A_TH
      REAL(DOUBLE)                    :: B12,B23,B34,B41   ! Coefficients in matrix A_5
      REAL(DOUBLE)                    :: B15,B25,B35,B45   ! Coefficients in matrix A_TH
      REAL(DOUBLE)                    :: A_5(4,2)          ! Matrix A5 in Alex Tessler paper (see ref above)
      REAL(DOUBLE)                    :: A_TH(4,8)         ! Matrix Atheta in Alex Tessler paper (see ref above)
      REAL(DOUBLE)                    :: A_5T(2,4)         ! Transpose of A_5
      REAL(DOUBLE)                    :: ATA(2,2)          ! A5'xA5
      REAL(DOUBLE)                    :: ATAI(2,2)         ! Inverse of ATA
      REAL(DOUBLE)                    :: C(2,8)            ! Part of B54
      REAL(DOUBLE)                    :: DET_ATA           ! determinant of matrix ATA
      REAL(DOUBLE)                    :: DUM28(2,8)        ! Intermediate matrix result
      REAL(DOUBLE)                    :: DUM312(3,12)      ! Intermediate matrix result
      REAL(DOUBLE)                    :: DUM315(3,15)      ! Intermediate matrix result
      REAL(DOUBLE)                    :: SA5,SB5   ! Sum of Ai5, Bi5 terms respectively
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

! **********************************************************************************************************************************
      IERROR = 0
      EPS1 = EPSIL(1)

      A12 = HALF*(XQL(1,1) - XQL(2,1))   ;   A15 = HALF*(XQL(1,1) - XQL(5,1))
      A23 = HALF*(XQL(2,1) - XQL(3,1))   ;   A25 = HALF*(XQL(2,1) - XQL(5,1))
      A34 = HALF*(XQL(3,1) - XQL(4,1))   ;   A35 = HALF*(XQL(3,1) - XQL(5,1))
      A41 = HALF*(XQL(4,1) - XQL(1,1))   ;   A45 = HALF*(XQL(4,1) - XQL(5,1))

      B12 = HALF*(XQL(1,2) - XQL(2,2))   ;   B15 = HALF*(XQL(1,2) - XQL(5,2))
      B23 = HALF*(XQL(2,2) - XQL(3,2))   ;   B25 = HALF*(XQL(2,2) - XQL(5,2))
      B34 = HALF*(XQL(3,2) - XQL(4,2))   ;   B35 = HALF*(XQL(3,2) - XQL(5,2))
      B41 = HALF*(XQL(4,2) - XQL(1,2))   ;   B45 = HALF*(XQL(4,2) - XQL(5,2))

      A_5(1,1) = B12   ;   A_5(1,2) = A12
      A_5(2,1) = B23   ;   A_5(2,2) = A23
      A_5(3,1) = B34   ;   A_5(3,2) = A34
      A_5(4,1) = B41   ;   A_5(4,2) = A41

      SA5 = +(A15 + A25 + A35 + A45)
      SB5 = +(B15 + B25 + B35 + B45)

      DO I=1,4
         DO J=1,2
            A_5T(J,I) = A_5(I,J)
         ENDDO
      ENDDO
      CALL MATMULT_FFF ( A_5T, A_5, 2, 4, 2, ATA )

      DET_ATA = ATA(1,1)*ATA(2,2) - ATA(1,2)*ATA(2,1)
      IF (DABS(DET_ATA) < EPS1) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IERROR = IERROR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1922) DET_ATA, SUBR_NAME
            WRITE(F06,1922) DET_ATA, SUBR_NAME
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1922
            ENDIF
         ENDIF
      ENDIF

      ATAI(1,1) =  ATA(2,2)/DET_ATA
      ATAI(1,2) = -ATA(1,2)/DET_ATA
      ATAI(2,1) = -ATA(2,1)/DET_ATA
      ATAI(2,2) =  ATA(1,1)/DET_ATA

      DO I=1,4
         DO J=1,8
            A_TH(I,J) = ZERO
         ENDDO
      ENDDO
      A_TH(1,1) = -B25   ;   A_TH(1,2) =  B15   ;   A_TH(1,5) = -A25   ;   A_TH(1,6) =  A15

      A_TH(2,2) = -B35   ;   A_TH(2,3) =  B25   ;   A_TH(2,6) = -A35   ;   A_TH(2,7) =  A25

      A_TH(3,3) = -B45   ;   A_TH(3,4) =  B35   ;   A_TH(3,7) = -A45   ;   A_TH(3,8) =  A35

      A_TH(4,1) =  B45   ;   A_TH(4,4) = -B15   ;   A_TH(4,5) =  A45   ;   A_TH(4,8) = -A15

      CALL MATMULT_FFF ( A_5T, A_TH , 2, 4, 8, DUM28 )
      CALL MATMULT_FFF ( ATAI, DUM28, 2, 2, 8, C )

! Formulate B54

      DO I=1,15
         DO J=1,12
            B54(I,J) = ZERO
         ENDDO
      ENDDO

      DO I=1,12
         B54(I,I) = ONE
      ENDDO

      B54(13, 1) = QUARTER
      B54(13, 2) =-QUARTER*(B15 + SB5*C(1,1) + SA5*C(2,1))
      B54(13, 3) =+QUARTER*(A15 + SB5*C(1,5) + SA5*C(2,5))

      B54(13, 4) = QUARTER
      B54(13, 5) =-QUARTER*(B25 + SB5*C(1,2) + SA5*C(2,2))
      B54(13, 6) =+QUARTER*(A25 + SB5*C(1,6) + SA5*C(2,6))

      B54(13, 7) = QUARTER
      B54(13, 8) =-QUARTER*(B35 + SB5*C(1,3) + SA5*C(2,3))
      B54(13, 9) =+QUARTER*(A35 + SB5*C(1,7) + SA5*C(2,7))

      B54(13,10) = QUARTER
      B54(13,11) =-QUARTER*(B45 + SB5*C(1,4) + SA5*C(2,4))
      B54(13,12) =+QUARTER*(A45 + SB5*C(1,8) + SA5*C(2,8))

      B54(14, 2) = C(1,1)
      B54(14, 3) = C(1,5)
      B54(14, 5) = C(1,2)
      B54(14, 6) = C(1,6)
      B54(14, 8) = C(1,3)
      B54(14, 9) = C(1,7)
      B54(14,11) = C(1,4)
      B54(14,12) = C(1,8)

      B54(15, 2) = C(2,1)
      B54(15, 3) = C(2,5)
      B54(15, 5) = C(2,2)
      B54(15, 6) = C(2,6)
      B54(15, 8) = C(2,3)
      B54(15, 9) = C(2,7)
      B54(15,11) = C(2,4)
      B54(15,12) = C(2,8)

! Reduce strain-displ matrices from 5 to 4 nodes

      DO K=1,5

         DO I=1,3
            DO J=1,15
               DUM315(I,J) = B2M_QQ_5(I,J,K)
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( DUM315, B54, 3, 15, 12, DUM312 )
         DO I=1,3
            DO J=1,12
               B2M_QQ_4(I,J,K) = DUM312(I,J)
            ENDDO
         ENDDO

         DO I=1,3
            DO J=1,15
               DUM315(I,J) = B3M_QQ_5(I,J,K)
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( DUM315, B54, 3, 15, 12, DUM312 )
         DO I=1,3
            DO J=1,12
               B3M_QQ_4(I,J,K) = DUM312(I,J)
            ENDDO
         ENDDO

      ENDDO

! Reduce thermal load matrix from 5 to 4 nodes

      IF (OPT_MIN4T(2) == 'Y') THEN

         CALL MATMULT_FFF_T ( B54, PTM_QQ_5, 15, 12, NTSUB, PTM_QQ_4 )

      ENDIF

! Reduce stress recovery matrices from 5 to 4 nodes

      IF (OPT_MIN4T(3) == 'Y') THEN

         DO K=1,5

            DO I=1,3
               DO J=1,15
                  DUM315(I,J) = S2M_QQ_5(I,J,K)
               ENDDO
            ENDDO

            CALL MATMULT_FFF ( DUM315, B54, 3, 15, 12, DUM312 )
            DO I=1,3
               DO J=1,12
                  S2M_QQ_4(I,J,K) = DUM312(I,J)
               ENDDO
            ENDDO

            DO I=1,3
               DO J=1,15
                  DUM315(I,J) = S3M_QQ_5(I,J,K)
               ENDDO
            ENDDO

            CALL MATMULT_FFF ( DUM315, B54, 3, 15, 12, DUM312 )
            DO I=1,3
               DO J=1,12
                  S3M_QQ_4(I,J,K) = DUM312(I,J)
               ENDDO
            ENDDO

         ENDDO

      ENDIF

! Reduce stiffness matrix from 5 to 4 nodes

      IF (OPT_MIN4T(4) == 'Y') THEN

         CALL MATMULT_FFF   ( KM_QQ_5, B54, 15, 15, 12, DUM2 )
         CALL MATMULT_FFF_T ( B54, DUM2, 15, 12, 12, KM_QQ_4 )

      ENDIF

! Reduce pressure load matrix from 5 to 4 nodes

      IF (OPT_MIN4T(5) == 'Y') THEN

         CALL MATMULT_FFF_T ( B54, PPM_QQ_5, 15, 12, NSUB, PPM_QQ_4 )

      ENDIF

      RETURN

! **********************************************************************************************************************************
 1922 FORMAT(' *ERROR  1922: THE DETERMINANT OF MATRIX ATA = ',1ES14.6,' IS TOO CLOSE TO ZERO IN SUBR ',A                          &
                    ,/,14X,' CANNOT INVERT MATRIX ATA IN SOLVING FOR MATRIX B54')

! **********************************************************************************************************************************

      END SUBROUTINE B54_REDUCTION

! ##################################################################################################################################

      SUBROUTINE STATIC_CONDENSATION

      USE CONSTANTS_1, ONLY           :  ZERO, ONE

      IMPLICIT NONE

      INTEGER(LONG)                   :: IERR              ! Error indicator from subr INVERT_33_MAT, called herein

      REAL(DOUBLE)                    :: B2AB( 3,12)       ! First 12 cols of B2M_QQ_5
      REAL(DOUBLE)                    :: B3AB( 3,12)       ! First 12 cols of B3M_QQ_5
      REAL(DOUBLE)                    :: B2O ( 3, 3)       ! Last   3 cols of B2M_QQ_5
      REAL(DOUBLE)                    :: B3O ( 3, 3)       ! Last   3 cols of B3M_QQ_5
      REAL(DOUBLE)                    :: DET_KOO           ! Determinant of A
      REAL(DOUBLE)                    :: DUM2(12,12)       ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM3(12,NTSUB)    ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM4(12,NSUB)     ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM51( 3,12)      ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM52( 3,12)      ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM61( 3,12)      ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM62( 3,12)      ! Intermediate matrix
      REAL(DOUBLE)                    :: GOA ( 3,12)       ! KOOI*KOA
      REAL(DOUBLE)                    :: GOAT(12, 3)       ! GOA(t)
      REAL(DOUBLE)                    :: KAAB(12,12)       ! The 12 x12 upper left  partition from KM_QQ_5 for node 5 of the elem
      REAL(DOUBLE)                    :: KOO ( 3, 3)       ! The  3 x 3 lower right partition from KM_QQ_5 for node 5 of the elem
      REAL(DOUBLE)                    :: KOA ( 3,12)       ! The  3 x12 lower left  partition from KM_QQ_5
      REAL(DOUBLE)                    :: KOOI( 3, 3)       ! Inverse of K00
      REAL(DOUBLE)                    :: PPAB(12,NSUB)     ! Upper 12 rows of PPM_QQ_5
      REAL(DOUBLE)                    :: PPO ( 3,NSUB)     ! Lower 3 rows of PPM_QQ_5
      REAL(DOUBLE)                    :: PTAB(12,NTSUB)    ! Upper 12 rows of PTM_QQ_5
      REAL(DOUBLE)                    :: PTO ( 3,NTSUB)    ! Lower 3 rows of PTM_QQ_5
      REAL(DOUBLE)                    :: S2AB( 3,12)       ! First 12 cols of S2M_QQ_5
      REAL(DOUBLE)                    :: S2O ( 3, 3)       ! Last   3 cols of S2M_QQ_5
      REAL(DOUBLE)                    :: S3AB( 3,12)       ! First 12 cols of S3M_QQ_5
      REAL(DOUBLE)                    :: S3O ( 3, 3)       ! Last   3 cols of S3M_QQ_5

! **********************************************************************************************************************************
! Calc KOA(t)*KOOI = GOAT

      DO I=1,12                                            ! Initialize
         DO J=1,12
            KAAB(I,J) = KM_QQ_5(I,J)
         ENDDO
      ENDDO

      DO I=13,15
         DO J=13,15
            KOO (I-12,J-12) = KM_QQ_5(I,J)
         ENDDO
      ENDDO

      DO I=13,15
         DO J=1,12
            KOA(I-12,J) = KM_QQ_5(I,J)
         ENDDO
      ENDDO

! Invert KOO

      CALL INVERT_33_MAT ( 'KOO', KOO, KOOI, 1.0D50*MACH_SFMIN, DET_KOO, IERR, SUM_DIAGS, SUM_OFF_DIAGS )
      IF (IERR /= 0) THEN
         WRITE(ERR,1947) EID, SUBR_NAME
         WRITE(F06,1947) EID, SUBR_NAME
         IF (DEBUG(188) == 0) THEN
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDIF

! Solve for GOA and transpose to get GOAT

      CALL MATMULT_FFF ( KOOI, KOA, 3, 3, 12, GOA )

      DO I=1,12
         DO J=1,3
            GOAT(I,J) = GOA(J,I)
         ENDDO
      ENDDO

! Reduce strain-displ matrices from 5 to 4 nodes by static condensation

      DO K=1,5

         DO I=1,3                                       ! Reduce B2M_QQ_5's
            DO J=1,12
               B2AB(I,J) = B2M_QQ_5(I,J,K)
            ENDDO
         ENDDO

         DO I=1,3
            DO J=13,15
               B2O(I,J-12) = B2M_QQ_5(I,J,K)
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( B2O , GOA , 3,  3, 12, DUM51 )
         CALL MATADD_FFF  ( B2AB, DUM51, 3, 12, ONE, -ONE, 0, DUM52 )
         DO I=1,3
            DO J=1,12
               B2M_QQ_4(I,J,K) = DUM52(I,J)
            ENDDO
         ENDDO

         DO I=1,3                                       ! Reduce B3M_QQ_5's
            DO J=1,12
               B3AB(I,J) = B3M_QQ_5(I,J,K)
            ENDDO
         ENDDO

         DO I=1,3
            DO J=13,15
               B3O(I,J-12) = B3M_QQ_5(I,J,K)
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( B3O , GOA , 3,  3, 12, DUM61 )
         CALL MATADD_FFF  ( B3AB, DUM61, 3, 12, ONE, -ONE, 0, DUM62 )
         DO I=1,3
            DO J=1,12
               B3M_QQ_4(I,J,K) = DUM62(I,J)
            ENDDO
         ENDDO

      ENDDO

! Reduce thermal loads from 5 to 4 nodes by static condensation

      IF (OPT_MIN4T(2) == 'Y') THEN

         DO I=1,12
            DO J=1,NTSUB
               PTAB(I,J) = PTM_QQ_5(I,J)
            ENDDO
         ENDDO

         DO I=13,15
            DO J=1,NTSUB
               PTO(I-12,J) = PTM_QQ_5(I,J)
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( GOAT, PTO, 12, 3, NTSUB, DUM3 )
         CALL MATADD_FFF  ( PTAB, DUM3, 12, NTSUB, ONE, -ONE, 0, PTM_QQ_4 )

      ENDIF

! Reduce stress recovery matrices from 5 to 4 nodes by static condensation

      IF (OPT_MIN4T(3) == 'Y') THEN

         DO K=1,5

            DO I=1,3                                       ! Reduce S2M_QQ_5's
               DO J=1,12
                  S2AB(I,J) = S2M_QQ_5(I,J,K)
               ENDDO
            ENDDO

            DO I=1,3
               DO J=13,15
                  S2O(I,J-12) = S2M_QQ_5(I,J,K)
               ENDDO
            ENDDO

            CALL MATMULT_FFF ( S2O , GOA , 3,  3, 12, DUM51 )
            CALL MATADD_FFF  ( S2AB, DUM51, 3, 12, ONE, -ONE, 0, DUM52 )
            DO I=1,3
               DO J=1,12
                  S2M_QQ_4(I,J,K) = DUM52(I,J)
               ENDDO
            ENDDO

            DO I=1,3                                       ! Reduce S3M_QQ_5's
               DO J=1,12
                  S3AB(I,J) = S3M_QQ_5(I,J,K)
               ENDDO
            ENDDO

            DO I=1,3
               DO J=13,15
                  S3O(I,J-12) = S3M_QQ_5(I,J,K)
               ENDDO
            ENDDO

            CALL MATMULT_FFF ( S3O , GOA , 3,  3, 12, DUM61 )
            CALL MATADD_FFF  ( S3AB, DUM61, 3, 12, ONE, -ONE, 0, DUM62 )
            DO I=1,3
               DO J=1,12
                  S3M_QQ_4(I,J,K) = DUM62(I,J)
               ENDDO
            ENDDO

         ENDDO

      ENDIF

! Reduce stiffness matrix from 5 to 4 nodes by static condensation

      IF (OPT_MIN4T(4) == 'Y') THEN

         CALL MATMULT_FFF ( GOAT, KOA , 12,  3,  12, DUM2 )
         CALL MATADD_FFF  ( KAAB, DUM2, 12, 12, ONE, -ONE, 0, KM_QQ_4 )

      ENDIF

! Reduce pressure loads from 5 to 4 nodes by static condensation

      IF (OPT_MIN4T(5) == 'Y') THEN

         DO I=1,12
            DO J=1,NSUB
               PPAB(I,J) = PPM_QQ_5(I,J)
            ENDDO
         ENDDO

         DO I=13,15
            DO J=1,NSUB
               PPO(I-12,J) = PPM_QQ_5(I,J)
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( GOAT, PPO, 12, 3, NSUB, DUM4 )
         CALL MATADD_FFF  ( PPAB, DUM4, 12, NSUB, ONE, -ONE, 0, PPM_QQ_4 )


      ENDIF

      IF (DEBUG(198) > 0) CALL DEB_QPLT3_2 ( OPT_MIN4T, 6, 'STC', GOAT, DUM2 , KAAB, KOO, KOA, KOOI,                               &
                                                            PTAB, PTO , DUM3 , PPAB, PPO, DUM4,                                    &
                                                            S2AB, S2O , DUM51, S3AB, S3O, DUM61)

      RETURN

! **********************************************************************************************************************************
 1947 FORMAT(' *ERROR  1947: CANNOT INVERT KOO MATRIX FOR CENTRAL NODE 5 IN THE MIN4T FORMULATION OF THE QUAD4 ELEMENT ',I8,       &
                           ' IN SUBR ',A                                                                                           &
                           ,/,14X,' USER CAN OVERRIDE THE PROGRAM ABORTING DUE TO THIS ERROR BY USING BULK DATA DEBUG 188 > 0')

! **********************************************************************************************************************************

      END SUBROUTINE STATIC_CONDENSATION

! ##################################################################################################################################

      SUBROUTINE INVERT_33_MAT ( MAT_A_NAME, A, AI, SMALL_NUM, DETA, IERR, SUM_DIAGS, SUM_OFF_DIAGS )

! Inverts an input 3x3 matrix A (and calls the inverse AI) and, optionally, checks that A*AI is the identity matrix

      USE PENTIUM_II_KIND, ONLY       :  BYTE, DOUBLE, LONG
      USE IOUNT1, ONLY                :  F06
      USE CONSTANTS_1, ONLY           :  ZERO, THREE

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME

      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator:
!                                                            IERR = 1 indicates DETA <= SMALL_NUM
!                                                            IERR = 2 indicates SUM_DIAGS     > SMALL_NUM
!                                                            IERR = 3 indicates SUM_OFF_DIAGS > SMALL_NUM
!                                                            IERR = 5 indicates SUM_DIAGS and SUM_OFF_DIAGS > SMALL_NUM
      INTEGER(LONG)                   :: I,J               ! DO loop indices

      REAL(DOUBLE) , INTENT(IN)       :: A(3,3)            ! Matrix to invert
      REAL(DOUBLE) , INTENT(OUT)      :: DETA              ! Determinant of A
      REAL(DOUBLE) , INTENT(IN)       :: SMALL_NUM         ! Small number for comparing to near zero
      REAL(DOUBLE) , INTENT(OUT)      :: AI(3,3)           ! Inverse of A
      REAL(DOUBLE) , INTENT(OUT)      :: SUM_DIAGS         ! Sum of diag terms from IDENT
      REAL(DOUBLE) , INTENT(OUT)      :: SUM_OFF_DIAGS     ! Sum of off-diag terms from IDENT
      REAL(DOUBLE)                    :: Z(3,3)            ! Adjoint of A

! **********************************************************************************************************************************
! Initialize

      IERR = 0
      DO I=1,3
         DO J=1,3
            AI(I,J) = ZERO
         ENDDO
      ENDDO

      SUM_DIAGS     = ZERO
      SUM_OFF_DIAGS = ZERO

      DETA = ( A(1,1)*A(2,2)*A(3,3) + A(2,1)*A(3,2)*A(1,3) + A(3,1)*A(2,3)*A(1,2) )                                                &
            -( A(1,3)*A(2,2)*A(3,1) + A(1,2)*A(2,1)*A(3,3) + A(1,1)*A(3,2)*A(2,3) )

      IF (DABS(DETA) <= SMALL_NUM) THEN

         IERR = 1
         WRITE(ERR,9000) EID
         WRITE(ERR,9001) DETA, SMALL_NUM
         WRITE(F06,9000) EID
         WRITE(F06,9001) DETA, SMALL_NUM

      ELSE

         Z(1,1)  =  A(2,2)*A(3,3) - A(2,3)*A(3,2)
         Z(1,2)  =  A(2,1)*A(3,3) - A(2,3)*A(3,1)
         Z(1,3)  =  A(2,1)*A(3,2) - A(2,2)*A(3,1)

         Z(2,1)  =  A(1,2)*A(3,3) - A(1,3)*A(3,2)
         Z(2,2)  =  A(1,1)*A(3,3) - A(1,3)*A(3,1)
         Z(2,3)  =  A(1,1)*A(3,2) - A(1,2)*A(3,1)

         Z(3,1)  =  A(1,2)*A(2,3) - A(1,3)*A(2,2)
         Z(3,2)  =  A(1,1)*A(2,3) - A(1,3)*A(2,1)
         Z(3,3)  =  A(1,1)*A(2,2) - A(1,2)*A(2,1)

         AI(1,1) =  Z(1,1)/DETA
         AI(1,2) = -Z(2,1)/DETA
         AI(1,3) =  Z(3,1)/DETA

         AI(2,1) = -Z(1,2)/DETA
         AI(2,2) =  Z(2,2)/DETA
         AI(2,3) = -Z(3,2)/DETA

         AI(3,1) =  Z(1,3)/DETA
         AI(3,2) = -Z(2,3)/DETA
         AI(3,3) =  Z(3,3)/DETA

         IF (DEBUG(199) > 0) THEN
            CALL CHECK_MAT_INVERSE ( MAT_A_NAME, A, AI, 3 )
         ENDIF

      ENDIF

      RETURN

!***********************************************************************************************************************************
!1001 FORMAT(' Input matrix is A with name =',A,//,' Input matrix A:',/,' --------------')

!1002 FORMAT(' AI = inverse of the input matrix :',/,' ---------------------------------')

!1003 FORMAT(' IDENT = A*AI :',/,' ------------')

!2001 FORMAT(5X,3(A,I3))

!2002 FORMAT(' Row',I3,':',3(1ES14.6))

!3001 FORMAT(' The determinant of A is DETA = ',1ES14.6,/)

!4001 FORMAT(' Sum of diagonal terms in IDENT = ',1ES14.6)

!4002 FORMAT(' Avg    diagonal term  in IDENT = ',1ES14.6)

!4003 FORMAT(' Sum of off diag terms in IDENT = ',1ES14.2,'% of a unity term in an identity matrix (which A*AI should be)')

 9000 FORMAT(' *INFORMATION: KOO MATRIX FOR CENTRAL NODE 5 ON MIN4T QUAD4 ELEMENT ',I8,' CANNOT BE INVERTED, OR THE INVERSE FAILED'&
                        ,/,14X,' THE CHECK THAT KOO TIMES ITS INVERSE IS THE IDENTITY MATRIX, DUE TO THE FOLLOWING:')

 9001 FORMAT('               IT CANNOT BE INVERTED SINCE ITS DETERMINANT = ',1ES9.2,' AND THE THRESHOLD FOR ZERO IS = ',1ES9.2)

!***********************************************************************************************************************************

      END SUBROUTINE INVERT_33_MAT

! ##################################################################################################################################

      SUBROUTINE DEB_QPLT3_1 ( OPT_MIN4T, WHICH )

      USE CONSTANTS_1, ONLY            : PI, ONE80
      USE PARAMS, ONLY                 : QUADAXIS
      USE MODEL_STUF, ONLY             : QUAD_DELTA

      CHARACTER(LEN=*) , INTENT(IN)   :: OPT_MIN4T(5)      ! 'Y'/'N' flags for whether certain elem matrices were calc'd

      INTEGER(LONG)    , INTENT(IN)   :: WHICH             ! Indicator of which debug output to write in the current call
      INTEGER(LONG)                   :: I,J,K             ! Local loop indices

! **********************************************************************************************************************************
      IF (WHICH == 0) THEN

         Write(f06,1101)

         Write(f06,*) 'Data for the quad element:'
         Write(f06,*) '-------------------------'
         Write(f06,*) '-------------------------'
         Write(f06,*)

         if      (quadaxis == 'SPLITD') then
            Write(f06,2101) (180/pi)*quad_delta
         else if (quadaxis == 'SIDE12') then
            Write(f06,2102)
         endif
         Write(f06,*)

         Write(f06,*) 'Matrix TE for the quad (from subr ELMGM2):'
         Write(f06,*) '-----------------------------------------'
         do i=1,3
            Write(f06,99887) (te(i,j),j=1,3)
         enddo
         Write(f06,*)

         Write(f06,*) 'Array XQL - local coords of the 5 nodes of the quad in the quad local coord system:'
         Write(f06,*) '----------------------------------------------------------------------------------'
         Write(f06,87651) xql(1,1), xql(1,2)
         Write(f06,87652) xql(2,1), xql(2,2)
         Write(f06,87653) xql(3,1), xql(3,2)
         Write(f06,87654) xql(4,1), xql(4,2)
         Write(f06,87655) xql(5,1), xql(5,2)
         Write(f06,*)

      ENDIF

! **********************************************************************************************************************************
! Output for triangle K

      IF ((WHICH > 0) .AND. (WHICH <= NUM_TRIAS)) THEN

         Write(f06,1102)
         Write(f06,*) 'Data for sub-triangle ', which,':'
         Write(f06,*) '------------------------'
         Write(f06,*) '------------------------'
         Write(f06,*)

         Write(f06,98712) tria_num, (v12(i),i=1,2), l12
         Write(f06,98713) tria_num, (v13(i),i=1,2), l13
         Write(f06,*)
         Write(f06,98714) tria_num, conv_rad_deg*theta
         Write(f06,*)
         Write(f06,98715) tria_num, conv_rad_deg*alpha_k
         Write(f06,*)
         Write(f06,98716) tria_num, x2tl, x3tl, y3tl
         Write(f06,*)

         Write(f06,*) 'Matrix TE_TRIA for triangle ',tria_num
            Write(f06,*) '------------------------------'
         do i=1,3
            Write(f06,98717) (te_tria(i,j),j=1,3)
         enddo
         Write(f06,*)

         Write(f06,*) 'Matrix TE_STRESS for triangle ',tria_num
            Write(f06,*) '--------------------------------'
         do i=1,3
            Write(f06,98717) (te_stress(i,j),j=1,3)
         enddo
         Write(f06,*)

         Write(f06,*) 'Matrix D  - triangle ', tria_num
            Write(f06,*) '-----------------------'
         do i=1,9,3
            Write(f06,98761) (d(i  ,j),j=1,9)
            Write(f06,98761) (d(i+1,j),j=1,9)
            Write(f06,98761) (d(i+2,j),j=1,9)
            Write(f06,*)
         enddo

         IF (OPT_MIN4T(2) == 'Y') THEN                                 ! Calculate thermal loads

            do j=1,ntsub
               Write(f06,76859) j,tria_num,(dt(i,j),i=1,4)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix PTV_TT  - triangle ',which,' load matrix in local tria coords and virgin DOF order:'
            Write(f06,*) '-----------------------------------------------------------------------------------------'
            do i=1,9,3
               Write(f06,98761) (ptv_tt(i  ,j),j=1,ntsub)
               Write(f06,98761) (ptv_tt(i+1,j),j=1,ntsub)
               Write(f06,98761) (ptv_tt(i+2,j),j=1,ntsub)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix PTM_TT  - triangle ',which,' thermal load matrix in local tria coords and MYSTRAN DOF order:'
            Write(f06,*) '--------------------------------------------------------------------------------------------------'
            do i=1,9,3
               Write(f06,98761) (ptm_tt(i  ,j),j=1,ntsub)
               Write(f06,98761) (ptm_tt(i+1,j),j=1,ntsub)
               Write(f06,98761) (ptm_tt(i+2,j),j=1,ntsub)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix PTM_TQ - triangle ',which,'  thermal load matrix in local quad coords and MYSTRAN DOF order:'
            Write(f06,*) '--------------------------------------------------------------------------------------------------'
            do i=1,9,3
               Write(f06,98761) (ptm_tq(i  ,j),j=1,ntsub)
               Write(f06,98761) (ptm_tq(i+1,j),j=1,ntsub)
               Write(f06,98761) (ptm_tq(i+2,j),j=1,ntsub)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix PTM_QQ - quad thermal load matrix for tria ',which,'  in quad 15 DOF coords:'
            Write(f06,*) '----------------------------------------------------------------------------------'
            do i=1,15,3
               Write(f06,98763) (ptm_qq(k,i  ,j),j=1,ntsub)
               Write(f06,98763) (ptm_qq(k,i+1,j),j=1,ntsub)
               Write(f06,98763) (ptm_qq(k,i+2,j),j=1,ntsub)
               Write(f06,*)
            enddo

         ENDIF

         IF (OPT_MIN4T(3) == 'Y') THEN                                 ! Stress recovery matrices

            Write(f06,*) 'Matrix S2V_TT - tria ',which,' bend stress recovery matrix in local tria coords and virgin DOF order:'
            Write(f06,*) '----------------------------------------------------------------------------------------------------'
            do i=1,3
               Write(f06,98761) (s2v_tt(i,j),j=1,9)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix S2M_TT - tria ',which,' bend stress recovery matrix in local tria coords and MYSTRAN DOF order:'
            Write(f06,*) '-----------------------------------------------------------------------------------------------------'
            do i=1,3
               Write(f06,98761) (s2m_tt(i,j),j=1,9)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix S2M_TQ - tria ',which,' bend stress recovery matrix in local quad coords and MYSTRAN DOF order:'
            Write(f06,*) '-----------------------------------------------------------------------------------------------------'
            do i=1,3
               Write(f06,98761) (s2m_tq(i,j),j=1,9)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix S2M_QQ - quad bend stress recovery matrix for tria ',which,'  in quad 15 DOF coords:'
            Write(f06,*) '------------------------------------------------------------------------------------------'
            do i=1,3
               Write(f06,98763) (s2m_qq(i,j,k),j=1,15)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix S3V_TT - tria ',which,' transv shear stress rec matrix in local tria coords and virgin DOF order:'
            Write(f06,*) '-------------------------------------------------------------------------------------------------------'
            do i=1,2
               Write(f06,98761) (s3v_tt(i,j),j=1,9)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix S3M_TT - tria ',which,' transv shear stress rec matrix in local tria coords and MYSTRAN DOF order:'
            Write(f06,*) '--------------------------------------------------------------------------------------------------------'
            do i=1,2
               Write(f06,98761) (s3m_tt(i,j),j=1,9)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix S3M_TQ - tria ',which,' transv shear stress rec matrix in local quad coords and MYSTRAN DOF order:'
            Write(f06,*) '--------------------------------------------------------------------------------------------------------'
            do i=1,2
               Write(f06,98761) (s3m_tq(i,j),j=1,9)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix S3M_QQ - quad transv shear stress recovery matrix for tria ',which,'  in quad 15 DOF coords:'
            Write(f06,*) '--------------------------------------------------------------------------------------------------'
            do i=1,2
               Write(f06,98763) (s3m_qq(i,j,k),j=1,15)
            enddo
            Write(f06,*)

         ENDIF

         IF (OPT_MIN4T(4) == 'Y') THEN                              ! Calculate stiffness matrix

            Write(f06,*) 'Matrix KV_TT  - triangle ',which,' stiffness matrix in local tria coords and virgin DOF order:'
            Write(f06,*) '---------------------------------------------------------------------------------------------'
            do i=1,9,3
               Write(f06,98761) (kv_tt(i  ,j),j=1,9)
               Write(f06,98761) (kv_tt(i+1,j),j=1,9)
               Write(f06,98761) (kv_tt(i+2,j),j=1,9)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix KM_TT  - triangle ',which,' stiffness matrix in local tria coords and MYSTRAN DOF order:'
            Write(f06,*) '----------------------------------------------------------------------------------------------'
            do i=1,9,3
               Write(f06,98761) (km_tt(i  ,j),j=1,9)
               Write(f06,98761) (km_tt(i+1,j),j=1,9)
               Write(f06,98761) (km_tt(i+2,j),j=1,9)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix KM_TQ - triangle ',which,' stiffness matrix in local quad coords and MYSTRAN DOF order:'
            Write(f06,*) '---------------------------------------------------------------------------------------------'
            do i=1,9,3
               Write(f06,98761) (km_tq(i  ,j),j=1,9)
               Write(f06,98761) (km_tq(i+1,j),j=1,9)
               Write(f06,98761) (km_tq(i+2,j),j=1,9)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix KM_QQ - quad stiffness matrix for tria ',which,' in quad 15 DOF coords:'
            Write(f06,*) '-----------------------------------------------------------------------------'
            do i=1,15,3
               Write(f06,98763) (km_qq(k,i  ,j),j=1,15)
               Write(f06,98763) (km_qq(k,i+1,j),j=1,15)
               Write(f06,98763) (km_qq(k,i+2,j),j=1,15)
               Write(f06,*)
            enddo

         ENDIF

         IF (OPT_MIN4T(5) == 'Y') THEN                              ! Calculate pressure loads

            Write(f06,*) 'Matrix PPV_TT  - triangle ',which,' pressure load matrix in local tria coords and virgin DOF order:'
            Write(f06,*) '--------------------------------------------------------------------------------------------------'
            do i=1,9,3
               Write(f06,98761) (ppv_tt(i  ,j),j=1,nsub)
               Write(f06,98761) (ppv_tt(i+1,j),j=1,nsub)
               Write(f06,98761) (ppv_tt(i+2,j),j=1,nsub)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix PPM_TT  - triangle ',which,' pressure load matrix in local tria coords and MYSTRAN DOF order:'
            Write(f06,*) '---------------------------------------------------------------------------------------------------'
            do i=1,9,3
               Write(f06,98761) (ppm_tt(i  ,j),j=1,nsub)
               Write(f06,98761) (ppm_tt(i+1,j),j=1,nsub)
               Write(f06,98761) (ppm_tt(i+2,j),j=1,nsub)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix PPM_TQ - triangle ',which,' pressure load matrix in local quad coords and MYSTRAN DOF order:'
            Write(f06,*) '--------------------------------------------------------------------------------------------------'
            do i=1,9,3
               Write(f06,98761) (ppm_tq(i  ,j),j=1,nsub)
               Write(f06,98761) (ppm_tq(i+1,j),j=1,nsub)
               Write(f06,98761) (ppm_tq(i+2,j),j=1,nsub)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix PPM_QQ - quad pressure load matrix for tria ',which,' in quad 15 DOF coords:'
            Write(f06,*) '----------------------------------------------------------------------------------'
            do i=1,15,3
               Write(f06,98763) (ppm_qq(k,i  ,j),j=1,nsub)
               Write(f06,98763) (ppm_qq(k,i+1,j),j=1,nsub)
               Write(f06,98763) (ppm_qq(k,i+2,j),j=1,nsub)
               Write(f06,*)
            enddo
         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Output after adding matrices for all triangles

      IF (WHICH == 5) THEN

         Write(f06,1102)
         Write(f06,*) 'Matrices after adding all 4 triangle matrices (for the 5 nodes):'
         Write(f06,*) '---------------------------------------------------------------'
         Write(f06,*) '---------------------------------------------------------------'
         Write(f06,*)

         IF (OPT_MIN4T(2) == 'Y') THEN                              ! Add thermal load matrices for all 4 triangles

            Write(f06,*) 'Matrix PTM_QQ_5 - quad thermal load for the 5 nodes:'
            Write(f06,*) '---------------------------------------------------'
            do i=1,15,3
               Write(f06,98763) (ptm_qq_5(i  ,j),j=1,ntsub)
               Write(f06,98763) (ptm_qq_5(i+1,j),j=1,ntsub)
               Write(f06,98763) (ptm_qq_5(i+2,j),j=1,ntsub)
               Write(f06,*)
            enddo
            Write(f06,*) 'Matrix PTM_QQ_5 - quad thermal load for the 5 nodes:'
            Write(f06,*) '---------------------------------------------------'
            k = 0
            do i=1,15
               do j=I,ntsub
                  if (dabs(ptm_qq_5(i,j)) > 0.d0) then
                     k = k + 1
                     Write(f06,98764) i, j, k, ptm_qq_5(i,j)
                  endif
               enddo
               Write(f06,*)
            enddo
         ENDIF

         IF (OPT_MIN4T(3) == 'Y') THEN

            do k=1,5

               Write(f06,*) 'Matrix S2M_QQ_5 - 5 node quad bending stress recovery matrix for stress recovery point ',k,':'
               Write(f06,*) '-----------------------------------------------------------------------------------------'
               do i=1,3
                  Write(f06,98763) (s2m_qq_5(i,j,k),j=1,15)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix S3M_QQ_5 - 5 node quad transverse shear stress recovery matrix for stress recovery point ',k,':'
               Write(f06,*) '--------------------------------------------------------------------------------------------------'
               do i=1,2
                  Write(f06,98763) (s3m_qq_5(i,j,k),j=1,15)
               enddo
               Write(f06,*)

            enddo

         ENDIF

         IF (OPT_MIN4T(4) == 'Y') THEN                              ! Add stiffness matrices for all 4 triangles

            Write(f06,*) 'Matrix KM_QQ_5 - quad stiffness for the 5 nodes:'
            Write(f06,*) '-----------------------------------------------'
            do i=1,15,3
               Write(f06,98763) (km_qq_5(i  ,j),j=1,15)
               Write(f06,98763) (km_qq_5(i+1,j),j=1,15)
               Write(f06,98763) (km_qq_5(i+2,j),j=1,15)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix KM_QQ_5 - quad stiffness for the 5 nodes:'
            Write(f06,*) '-----------------------------------------------'
            k = 0
            do i=1,15
               do j=I,15
                  if (dabs(km_qq_5(i,j)) > 0.d0) then
                     k = k + 1
                     Write(f06,98765) i, j, k, km_qq_5(i,j)
                  endif
               enddo
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix KM_QQ_5 - checking symmetry:'
            Write(f06,*) '----------------------------------'
            k = 0
            do i=1,15
               do j=i+1,15
                  if (i /= j) then
                     Write(f06,97765) i, j, km_qq_5(i,j),km_qq_5(j,i),(km_qq_5(i,j)-km_qq_5(j,i))
                  endif
               enddo
               Write(f06,*)
            enddo

         ENDIF

         IF (OPT_MIN4T(5) == 'Y') THEN                              ! Add pressure load matrices for all 4 triangles

            Write(f06,*) 'Matrix PPM_QQ_5 - quad pressure load for the 5 nodes:'
            Write(f06,*) '----------------------------------------------------'
            do i=1,15,3
               Write(f06,98763) (ppm_qq_5(i  ,j),j=1,nsub)
               Write(f06,98763) (ppm_qq_5(i+1,j),j=1,nsub)
               Write(f06,98763) (ppm_qq_5(i+2,j),j=1,nsub)
               Write(f06,*)
            enddo

            Write(f06,*) 'Matrix PPM_QQ_5 - quad pressure loads for the 5 nodes:'
            Write(f06,*) '-----------------------------------------------------'
            k = 0
            do i=1,15
               do j=I,nsub
                  if (dabs(ppm_qq_5(i,j)) > 0.d0) then
                     k = k + 1
                     Write(f06,98766) i, j, k, ppm_qq_5(i,j)
                  endif
               enddo
               Write(f06,*)
            enddo

         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Output after reducing 5 node quad to 4 node quad

      IF (WHICH == 7) THEN

         Write(f06,1102)
         Write(f06,*) 'Matrices after reducing 5 node quad to 4 nodes:'
         Write(f06,*) '-----------------------------------------------'
         Write(f06,*) '-----------------------------------------------'
         Write(f06,*)

         IF (OPT_MIN4T(2) == 'Y') THEN

            Write(f06,*) 'Matrix PTM_QQ_4 - quad thermal loads for the 4 nodes:'
            Write(f06,*) '----------------------------------------------------'
            do i=1,12,3
               Write(f06,98763) (ptm_qq_4(i  ,j),j=1,ntsub)
               Write(f06,98763) (ptm_qq_4(i+1,j),j=1,ntsub)
               Write(f06,98763) (ptm_qq_4(i+2,j),j=1,ntsub)
               Write(f06,*)
            enddo

         ENDIF

         IF (OPT_MIN4T(3) == 'Y') THEN

            do k=1,5

               Write(f06,*) 'Matrix S2M_QQ_4 - 4 node quad bending stress recovery matrix for stress recovery point ',k,':'
               Write(f06,*) '-----------------------------------------------------------------------------------------'
               do i=1,3
                  Write(f06,98763) (s2m_qq_4(i,j,k),j=1,12)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix S3M_QQ_4 - 4 node quad transverse shear stress recovery matrix for stress recovery point ',k,':'
               Write(f06,*) '--------------------------------------------------------------------------------------------------'
               do i=1,2
                  Write(f06,98763) (s3m_qq_4(i,j,k),j=1,12)
               enddo
               Write(f06,*)

            enddo

         ENDIF

         IF (OPT_MIN4T(4) == 'Y') THEN

            Write(f06,*) 'Matrix KM_QQ_4 - quad stiffness for the 4 nodes:'
            Write(f06,*) '-----------------------------------------------'
            do i=1,12,3
               Write(f06,98763) (km_qq_4(i  ,j),j=1,12)
               Write(f06,98763) (km_qq_4(i+1,j),j=1,12)
               Write(f06,98763) (km_qq_4(i+2,j),j=1,12)
               Write(f06,*)
            enddo

         ENDIF

         IF (OPT_MIN4T(5) == 'Y') THEN

            Write(f06,*) 'Matrix PPM_QQ_4 - quad pressure loads for the 4 nodes:'
            Write(f06,*) '-----------------------------------------------------'
            do i=1,12,3
               Write(f06,98763) (ppm_qq_4(i  ,j),j=1,nsub)
               Write(f06,98763) (ppm_qq_4(i+1,j),j=1,nsub)
               Write(f06,98763) (ppm_qq_4(i+2,j),j=1,nsub)
               Write(f06,*)
            enddo

         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Final matrices for the element in MYSTRAN DOF order and 6 DOF/grid

      IF (WHICH == 8) THEN

         Write(f06,1102)
         Write(f06,*) 'Final matrices for the 4 grid and 6 DOF/grid quad :'
         Write(f06,*) '--------------------------------------------------'
         Write(f06,*) '--------------------------------------------------'
         Write(f06,*)

         IF (OPT_MIN4T(2) == 'Y') THEN

            Write(f06,*) 'Matrix PTE - quad thermal loads in MYSTRAN 6 DOF/grid format:'
            Write(f06,*) '------------------------------------------------------------'
            do i=1,24,6
               Write(f06,98773) (pte(i  ,j),j=1,ntsub)
               Write(f06,98773) (pte(i+1,j),j=1,ntsub)
               Write(f06,98773) (pte(i+2,j),j=1,ntsub)
               Write(f06,98773) (pte(i+3,j),j=1,ntsub)
               Write(f06,98773) (pte(i+4,j),j=1,ntsub)
               Write(f06,98773) (pte(i+5,j),j=1,ntsub)
               Write(f06,*)
            enddo

         ENDIF

         IF (OPT_MIN4T(3) == 'Y') THEN

            Write(f06,*) 'Matrix BE2 - quad bending strain recovery matrix in MYSTRAN 6 DOF/grid format:'
            Write(f06,*) '------------------------------------------------------------------------------'
            do i=1,3
               Write(f06,98773) (be2(i,j,1),j=1,24)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix BE3 - quad transverse shear strain recovery matrix in MYSTRAN 6 DOF/grid format:'
            Write(f06,*) '---------------------------------------------------------------------------------------'
            do i=1,2
               Write(f06,98773) (be3(i,j,1),j=1,24)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix SE2 - quad bending stress recovery matrix in MYSTRAN 6 DOF/grid format:'
            Write(f06,*) '-----------------------------------------------------------------------------'
            do i=1,3
               Write(f06,98773) (se2(i,j,1),j=1,24)
            enddo
            Write(f06,*)

            Write(f06,*) 'Matrix SE3 - quad transverse shear stress recovery matrix in MYSTRAN 6 DOF/grid format:'
            Write(f06,*) '--------------------------------------------------------------------------------------'
            do i=1,2
               Write(f06,98773) (se3(i,j,1),j=1,24)
            enddo
            Write(f06,*)

         ENDIF

         IF (OPT_MIN4T(4) == 'Y') THEN

            Write(f06,*) 'Matrix KE - quad stiffness in MYSTRAN 6 DOF/grid format:'
            Write(f06,*) '-------------------------------------------------------'
            do i=1,24,6
               Write(f06,98773) (ke(i  ,j),j=1,24)
               Write(f06,98773) (ke(i+1,j),j=1,24)
               Write(f06,98773) (ke(i+2,j),j=1,24)
               Write(f06,98773) (ke(i+3,j),j=1,24)
               Write(f06,98773) (ke(i+4,j),j=1,24)
               Write(f06,98773) (ke(i+5,j),j=1,24)
               Write(f06,*)
            enddo

         ENDIF

         IF (OPT_MIN4T(5) == 'Y') THEN

            Write(f06,*) 'Matrix PE - quad pressure loads in MYSTRAN 6 DOF/grid format:'
            Write(f06,*) '------------------------------------------------------------'
            do i=1,24,6
               Write(f06,98773) (ppe(i  ,j),j=1,nsub)
               Write(f06,98773) (ppe(i+1,j),j=1,nsub)
               Write(f06,98773) (ppe(i+2,j),j=1,nsub)
               Write(f06,98773) (ppe(i+3,j),j=1,nsub)
               Write(f06,98773) (ppe(i+4,j),j=1,nsub)
               Write(f06,98773) (ppe(i+5,j),j=1,nsub)
               Write(f06,*)
            enddo

         ENDIF

         Write(f06,3102)

      ENDIF

! **********************************************************************************************************************************
 1101 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' :::::::::::::::::::::::::::::::::::::::::::START DEBUG(198) OUTPUT FROM SUBROUTINE QPLT3::::::::::::::::::::::::::',&
             ':::::::::::::::::',/)

 1102 FORMAT(' ******************************************************************************************************************',&
             '*****************')

 2101 FORMAT(' The local x axis for the quad splits the angle between the 2 diagonals.',                                           &
             ' The angle between side 1-2 and the x axis = ',f7.3)

 2102 FORMAT(' The local x axis is along side 1-2 of the quad')

 3102 FORMAT(' ::::::::::::::::::::::::::::::::::::::::::::END DEBUG(198) OUTPUT FROM SUBROUTINE QPLT3:::::::::::::::::::::::::::',&
             ':::::::::::::::::'                                                                                                ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

76859 format(' For internal subcase ',i3,' DT(1-4) for triangle ',i2,' = ',4(1es14.6))

87651 FORMAT('  node 1 = ',2(1ES14.6))

87652 FORMAT('  node 2 = ',2(1ES14.6))

87653 FORMAT('  node 3 = ',2(1ES14.6))

87654 FORMAT('  node 4 = ',2(1ES14.6))

87655 FORMAT('  node 5 = ',2(1ES14.6))

97765 format(2i4,3(1es14.6))

98712 format(' Vector along side 1-2 of sub-triangle ',i2,' is: V12(1-2)     = ',2(1es13.6),'. Length of side 1-2 = ',1es14.6)

98713 format(' Vector along side 1-3 of sub-triangle ',i2,' is: V13(1-2)     = ',2(1es13.6),'. Length of side 1-3 = ',1es14.6)

98714 format(' Angle THETA between side 1-2 of sub-triangle ',i2,' and quad x axis = ',f9.3)

98715 format(' Angle ALPHA between sides 1-2 and 1-3 of sub-triangle ',i2,' = ',f9.3)

98716 format(' Local coords of sub-triangle ',i2,' nodes are: X2TL, X3TL, Y3TL = ',3(1es14.6))

98717 format(3(1es14.6))

98761 format(3(3(1es14.6),3X))

98763 format(5(3(1es14.6),3X))

98764 format(' I, J, count, PTM_QQ_5(I,J) = ',3i4,3x,1es21.14)

98765 format(' I, J, count, KM_QQ_5(I,J)  = ',3i4,3x,1es21.14)

98766 format(' I, J, count, PPM_QQ_5(I,J) = ',3i4,3x,1es21.14)

98773 format(4(6(1es14.6),3X))

99887 format(3(1es14.6))

! **********************************************************************************************************************************

      END SUBROUTINE DEB_QPLT3_1

! ##################################################################################################################################

      SUBROUTINE DEB_QPLT3_2 ( OPT_MIN4T, WHICH, REDUCTION_METH,                                                                   &
                               GOAT, DUM2, KAAB, KOO, KOA, KOOI,                                                                   &
                               PTAB, PTO , DUM3, PPAB, PPO, DUM4,                                                                  &
                               S2AB, S2O , DUM5, S3AB, S3O, DUM6)

      CHARACTER(LEN=*), INTENT(IN)    :: OPT_MIN4T(5)            ! 'Y'/'N' flags for whether certain elem matrices were calc'd
      CHARACTER(LEN=*), INTENT(IN)    :: REDUCTION_METH    ! B54 or STC for the method used to reduce the 5 node matrices to 4 nodes

      INTEGER(LONG)   , INTENT(IN)    :: WHICH             ! Indicator of which debug output to Write in the current call

      REAL(DOUBLE)    , INTENT(IN)    :: DUM2(12,12)       ! Intermediate matrix
      REAL(DOUBLE)    , INTENT(IN)    :: DUM3(12,NTSUB)    ! Intermediate matrix
      REAL(DOUBLE)    , INTENT(IN)    :: DUM4(12,NSUB)     ! Intermediate matrix
      REAL(DOUBLE)    , INTENT(IN)    :: DUM5( 3,12)       ! Intermediate matrix
      REAL(DOUBLE)    , INTENT(IN)    :: DUM6( 3,12)       ! Intermediate matrix
      REAL(DOUBLE)    , INTENT(IN)    :: GOAT(12, 3)       ! Intermediate matrix
      REAL(DOUBLE)    , INTENT(IN)    :: KAAB(12,12)       ! The 12 x12 upper left  partition from KM_QQ_5 for node 5 of the elem
      REAL(DOUBLE)    , INTENT(IN)    :: KOO ( 3, 3)       ! The  3 x 3 lower right partition from KM_QQ_5 for node 5 of the elem
      REAL(DOUBLE)    , INTENT(IN)    :: KOA ( 3,12)       ! The  3 x12 lower left  partition from KM_QQ_5
      REAL(DOUBLE)    , INTENT(IN)    :: KOOI( 3, 3)       ! Inverse of K00
      REAL(DOUBLE)    , INTENT(IN)    :: PPAB(12,NSUB)     ! Upper 12 rows of PPM_QQ_5
      REAL(DOUBLE)    , INTENT(IN)    :: PPO ( 3,NSUB)     ! Lower 3 rows of PPM_QQ_5
      REAL(DOUBLE)    , INTENT(IN)    :: PTAB(12,NTSUB)    ! Upper 12 rows of PTM_QQ_5
      REAL(DOUBLE)    , INTENT(IN)    :: PTO ( 3,NTSUB)    ! Lower 3 rows of PTM_QQ_5
      REAL(DOUBLE)    , INTENT(IN)    :: S2AB( 3,12)       ! First 12 cols of S2M_QQ_5
      REAL(DOUBLE)    , INTENT(IN)    :: S2O ( 3, 3)       ! Last   3 cols of S2M_QQ_5
      REAL(DOUBLE)    , INTENT(IN)    :: S3AB( 2,12)       ! First 12 cols of S2M_QQ_5
      REAL(DOUBLE)    , INTENT(IN)    :: S3O ( 2, 3)       ! Last   3 cols of S2M_QQ_5

! **********************************************************************************************************************************
! Output in process of reducing 5 node quad to 4 node quad

      IF (WHICH == 6) THEN

         IF (REDUCTION_METH == 'B54') THEN


         ENDIF

         IF (REDUCTION_METH == 'STC') THEN

            Write(f06,1102)
            Write(f06,*) 'Matrices in the reduction of 5 to 4 nodes via the static condensation method:'
            Write(f06,*) '----------------------------------------------------------------------------'
            Write(f06,*) '----------------------------------------------------------------------------'
            Write(f06,*)

            Write(f06,*) 'Matrix GOAT = KOA(t)*KOOI:'
            Write(f06,*) '-------------------------'
            do i=1,12,3
               Write(f06,97763) (goat(i  ,j),j=1,3)
               Write(f06,97763) (goat(i+1,j),j=1,3)
               Write(f06,97763) (goat(i+2,j),j=1,3)
               Write(f06,*)
            enddo
            Write(f06,*)

            IF (OPT_MIN4T(2) == 'Y') THEN

               Write(f06,*) 'Matrix PTAB (columns are for each thermal case):'
               Write(f06,*) '-----------------------------------------------'
               do i=1,12
                  Write(f06,95001) (ptab(i,j),j=1,ntsub)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix PTO (columns are for each thermal case):'
               Write(f06,*) '----------------------------------------------'
               do i=1,3
                  Write(f06,95001) (pto(i,j),j=1,ntsub)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix DUM3 = GOAT*PTO (columns are for each thermal case):'
               Write(f06,*) '----------------------------------------------------------'
               do i=1,12
                  Write(f06,95001) (dum3(i,j),j=1,ntsub)
               enddo
               Write(f06,*)

            ENDIF

            IF (OPT_MIN4T(3) == 'Y') THEN

               Write(f06,*) 'Matrix S2AB:'
               Write(f06,*) '------------'
               do i=1,3
                  Write(f06,96001) (s2ab(i,j),j=1,12)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix S2O:'
               Write(f06,*) '-----------'
               do i=1,3
                  Write(f06,96001) (s2o(i,j),j=1,3)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix DUM5 = S2O*GOA:'
               Write(f06,*) '----------------------'
               do i=1,3
                  Write(f06,96001) (dum5(i,j),j=1,12)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix S3AB:'
               Write(f06,*) '------------'
               do i=1,2
                  Write(f06,96001) (s3ab(i,j),j=1,12)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix S3O:'
               Write(f06,*) '-----------'
               do i=1,2
                  Write(f06,96001) (s3o(i,j),j=1,3)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix DUM5 = S2O*GOA:'
               Write(f06,*) '----------------------'
               do i=1,2
                  Write(f06,96001) (dum6(i,j),j=1,12)
               enddo
               Write(f06,*)

            ENDIF

            IF (OPT_MIN4T(4) == 'Y') THEN

               Write(f06,*) 'Matrix KAAB:'
               Write(f06,*) '-----------'
               do i=1,12,3
                  Write(f06,97763) (kaab(i  ,j),j=1,12)
                  Write(f06,97763) (kaab(i+1,j),j=1,12)
                  Write(f06,97763) (kaab(i+2,j),j=1,12)
                  Write(f06,*)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix KOO:'
               Write(f06,*) '----------'
               do i=1,3,3
                  Write(f06,97763) (koo(i  ,j),j=1,3)
                  Write(f06,97763) (koo(i+1,j),j=1,3)
                  Write(f06,97763) (koo(i+2,j),j=1,3)
                  Write(f06,*)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix KOA:'
               Write(f06,*) '----------'
               do i=1,3,3
                  Write(f06,97763) (koa(i  ,j),j=1,12)
                  Write(f06,97763) (koa(i+1,j),j=1,12)
                  Write(f06,97763) (koa(i+2,j),j=1,12)
                  Write(f06,*)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix KOOI:'
               Write(f06,*) '-----------'
               do i=1,3,3
                  Write(f06,97763) (kooi(i  ,j),j=1,3)
                  Write(f06,97763) (kooi(i+1,j),j=1,3)
                  Write(f06,97763) (kooi(i+2,j),j=1,3)
                  Write(f06,*)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix DUM2 = KOA(t)*KOOI*KOA:'
               Write(f06,*) '-----------------------------'
               do i=1,12,3
                  Write(f06,97763) (dum2(i  ,j),j=1,12)
                  Write(f06,97763) (dum2(i+1,j),j=1,12)
                  Write(f06,97763) (dum2(i+2,j),j=1,12)
                  Write(f06,*)
               enddo
               Write(f06,*)

            ENDIF

            IF (OPT_MIN4T(5) == 'Y') THEN

               Write(f06,*) 'Matrix PPAB (columns are for each subcase):'
               Write(f06,*) '------------------------------------------'
               do i=1,12
                  Write(f06,98001) (ppab(i,j),j=1,nsub)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix PPO (columns are for each subcase):'
               Write(f06,*) '-----------------------------------------'
               do i=1,3
                  Write(f06,98001) (ppo(i,j),j=1,nsub)
               enddo
               Write(f06,*)

               Write(f06,*) 'Matrix DUM4 = GOAT*PPO (columns are for each subcase):'
               Write(f06,*) '-----------------------------------------------------'
               do i=1,12
                  Write(f06,98001) (dum4(i,j),j=1,nsub)
               enddo
               Write(f06,*)

            ENDIF

         ENDIF

      ENDIF

! **********************************************************************************************************************************
 1102 FORMAT(' ******************************************************************************************************************',&
             '*****************')

95001 format(32767(1es14.6))

96001 format(12(1es14.6))

97763 format(4(6(1es14.6),3X))

98001 format(32767(1es14.6))

! **********************************************************************************************************************************

      END SUBROUTINE DEB_QPLT3_2

      END SUBROUTINE QPLT3

   END MODULE MIN4T
