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

   MODULE SOLID_ELEMENTS

   USE RESULT_COORDINATES, ONLY: ELMDIS

   USE QUADRATURE, ONLY          :  ORDER_GAUSS, ORDER_TRIA, ORDER_TETRA
   USE SHAPE_FUNCTIONS_3D, ONLY  :  SHP3DH, SHP3DP, SHP3DT
   USE JACOBIAN, ONLY             :  JAC2D, JAC3D
   USE BMATRIX_3D, ONLY           :  B3D_ISOPARAMETRIC

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: HEXA, PENTA, TETRA

   CONTAINS

      SUBROUTINE HEXA ( OPT, INT_ELEM_ID,IORD, RED_INT_SHEAR, BUBBLE, WRITE_WARN )

! Isoparametric hexahedron solid element (8 or 20 nodes)

! Formulation (selected by PSOLID IN and ISOP as for the Nastran CHEXA; the defaults are BUBBLE for HEXA8 and THREE for HEXA20,
! with ISOP blank):

!  HEXA8 bubble element     : IN = BUBBLE, or blank with ISOP blank/REDUCED (BUBBLE = 'Y').
!                             Enhanced assumed strain element Q1/E21 (Andelfinger & Ramm 1993). 21 internal strain modes,
!                             defined in the natural (xi,eta,zeta) frame and symmetric in the three directions, are pushed
!                             forward with the centroid Jacobian and scaled by det(J0)/det(J) (Simo & Rifai 1990) so the
!                             element passes the patch test on distorted meshes. The modes are statically condensed at the
!                             element level. Relieves shear (bending) and volumetric locking without any preferred direction,
!                             is invariant to rigid rotation and to node relabeling, and has no spurious zero-energy modes.
!  ISOP blank/REDUCED       : Reduced shear (RED_INT_SHEAR = 'Y').
!                             IN = TWO: Jacobian-weighted averages of the shear rows of B over the Gauss points of opposite faces.
!                             IN = THREE: shear strain energy integrated with 2x2x2 points, the rest with 3x3x3.
!  ISOP = FULL              : Standard displacement-based isoparametric element with full Gaussian integration.

! For the HEXA8 bubble element the same enhanced strain-displacement operator is used consistently for the stiffness, the thermal
! loads and the stress/strain recovery matrices.

! Subroutine calculates:

!  1) ME        = element mass matrix                  , if OPT(1) = 'Y'
!  2) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  3) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  4) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  5) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  6) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, MELDOF, MPLOAD4_3D_DATA, NPLOAD4_3D, NSUB, NTSUB
      USE CONSTANTS_1, ONLY           :  QUARTER, HALF, ZERO, ONE, TWO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  AGRID, ALPVEC, BE1, BE2, DT, EID, ELGP, NUM_EMG_FATAL_ERRS, ES, KE, KED, ME,              &
                                         NUM_EMG_FATAL_ERRS, PLOAD4_3D_DATA, PPE, PRESS, PTE, RHO, SE1, SE2, STE1, STE2, STE3,     &
                                         STRESS, TREF, TYPE, XEL
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE RESULT_COORDINATES, ONLY    :  ELMDIS
      USE ELEMENT_RECOVERY_SUPPORT, ONLY:  ELEM_STRE_STRN_ARRAYS
      USE MASS_DOF_EXPANSION, ONLY    :  EXPAND_MASS_DOFS


      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'HEXA'
      CHARACTER( 1*BYTE), INTENT(IN)  :: RED_INT_SHEAR          ! 'Y' for ISOP blank/REDUCED: reduced shear formulation (not used
!                                                                 by the HEXA8 bubble element). 'N' for ISOP = FULL: standard
!                                                                 displacement element
      CHARACTER( 1*BYTE), INTENT(IN)  :: BUBBLE                 ! 'Y' for the HEXA8 bubble element: Q1/E21 enhanced assumed strain
      CHARACTER( 1*BYTE), INTENT(IN)  :: OPT(6)                 ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN             ! If 'Y" write warning messages, otherwise do not

      CHARACTER( 1*BYTE)              :: GENERATE_PRESS_LOAD    ! If 'Y' this element will have a pressure load on at least 1 face
!                                                                 for at least 1 subcase

      CHARACTER(46*BYTE)              :: IORD_MSG               ! Character name of an integration order (used for debug output)

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID             ! Internal element ID
      INTEGER(LONG), INTENT(IN)       :: IORD                   ! Gaussian integration order for element
      INTEGER(LONG), PARAMETER        :: EAS_NUM_MODES = 21_LONG ! Q1/E21 enhanced strain modes for HEXA8
      INTEGER(LONG)                   :: ELIDA,ELIDI            ! Actual and internal elem ID's read from array PLOAD4_3D_DATA
      INTEGER(LONG)                   :: FACE_COUNT             ! Count of the num of HEXA faces processed for PLOAD4 pressure loads
      INTEGER(LONG)                   :: FACE_NODES(6,4)        ! Array of actual   grid numbers on the 6 faces of the HEXA
      INTEGER(LONG)                   :: FACE_AGRDS(6,4)        ! Array of internal grid numbers on the 6 faces of the HEXA
      INTEGER(LONG)                   :: FACE_NUM               ! Number of the face that the pressure acts on (may be more than 1)
      INTEGER(LONG)                   :: GAUSS_PT               ! Gauss point number (used for DEBUG output in subr SHP3DH
      INTEGER(LONG)                   :: GA,GB                  ! Actual grid numbers read from array PLOAD4_3D_DATA
      INTEGER(LONG)                   :: I,J,K,L,M,N            ! DO loop indices
      INTEGER(LONG)                   :: IDELT                  !
      INTEGER(LONG)                   :: IDOF                   !
      INTEGER(LONG)                   :: IERR                   ! Local error count
      INTEGER(LONG)                   :: II,JJ                  ! Counters
      INTEGER(LONG)                   :: ID(3*ELGP)             ! Array which shows equivalence of DOF's in virgin element with the
!                                                                 6 DOF/grid of the final element stiffness matrix
      INTEGER(LONG)                   :: IORD_SH                ! Gaussian integration order for shear for HEXA20 element
      INTEGER(LONG)                   :: ISCNUM                 ! Internal subcase number read from array PLOAD4_3D_DATA
      INTEGER(LONG)                   :: K1,K2,K3,K4            ! Array indices
      INTEGER(LONG)                   :: K5,K6,K7,K8            ! Array indices

      INTEGER(LONG)                   :: STR_PT_NUM             ! Stress point number. 1 is center, 2+ are element nodes 1+.

      REAL(DOUBLE)                    :: ALP(6)                 ! First col of ALPVEC

      REAL(DOUBLE)                    :: B(6,3*ELGP,IORD*IORD*IORD)
                                                                ! Strain-displ matrix for this element for all Gauss points

      REAL(DOUBLE)                    :: BI(6,3*ELGP)           ! Strain-displ matrix for this element for one Gauss point
      REAL(DOUBLE)                    :: B_SH(6,3*ELGP,8)       ! Shear rows of B at the 2x2x2 points (reduced shear integration, IORD = 3)
      REAL(DOUBLE)                    :: DETJ_SH(8)             ! Jacobian determinants at the 2x2x2 points
      LOGICAL                         :: RED3                   ! True for IORD = 3 with reduced integration of the shear terms

      REAL(DOUBLE)                    :: CBAR(3,3*ELGP)         ! Derivatives of shape fcns wrt x,y,z used in diff stiff matrix
!                                                                 (contains terms from DPSHX matrices for each grid of the HEXA)

      REAL(DOUBLE)                    :: DETJ(IORD*IORD*IORD)   ! Determinant of JAC for all Gauss points
      REAL(DOUBLE)                    :: DPSHG(3,ELGP)          ! Output from subr SHP3DH. Derivatives of PSH wrt elem isopar coords
      REAL(DOUBLE)                    :: DPSHX(3,ELGP)          ! Derivatives of PSH wrt elem x, y coords.
      REAL(DOUBLE)                    :: DUM0(3*ELGP)           ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM1(3*ELGP)           ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM2(6,3*ELGP)         ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM3(3*ELGP,3*ELGP)    ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM4(6,3*ELGP)         ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM5(3*ELGP,3*ELGP)    ! Intermediate matrix used in solving for KE elem matrices
      REAL(DOUBLE)                    :: DUM6(3,3*ELGP)         ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM7(ELGP,ELGP)        ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM8(ELGP,ELGP)        ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM9(3,ELGP)           ! Intermediate matrix used in solving for elem matrices

      REAL(DOUBLE)                    :: EAS_G(EAS_NUM_MODES,3*ELGP)
                                                                ! HEXA8 EAS coupling matrix  G = integral(M' * ES * B dV)
      REAL(DOUBLE)                    :: EAS_HFAC(EAS_NUM_MODES,EAS_NUM_MODES)
                                                                ! Cholesky factor of the equilibrated internal-mode stiffness
!                                                                 D*H*D, H = integral(M' * ES * M dV), D = diag(H)^(-1/2)
      REAL(DOUBLE)                    :: EAS_HSCL(EAS_NUM_MODES)! Diagonal equilibration D of H
      REAL(DOUBLE)                    :: EAS_HINV_G(EAS_NUM_MODES,3*ELGP)
                                                                ! inv(H)*G (internal modes in terms of the element displacements)
      REAL(DOUBLE)                    :: EAS_HINV_Q(EAS_NUM_MODES,NTSUB)
                                                                ! inv(H)*q for the internal-mode thermal load of each thermal subcase,
!                                                                 q = integral(M' * ES * eps_th dV)
      REAL(DOUBLE)                    :: EAS_MODE(6,EAS_NUM_MODES,IORD*IORD*IORD)
                                                                ! Enhanced strain modes M (Cartesian, det-scaled) at each Gauss pt
      REAL(DOUBLE)                    :: EAS_MODE_PT(6,EAS_NUM_MODES)
                                                                ! Enhanced strain modes at a recovery point
      REAL(DOUBLE)                    :: EAS_T0(6,6)            ! Natural -> Cartesian strain transformation built from J0

      REAL(DOUBLE)                    :: EALP(6)                ! Variable used in calc PTE therm loads & STEi therm stress coeffs
      REAL(DOUBLE)                    :: EPS1                   ! A small number to compare to real zero
      REAL(DOUBLE)                    :: EPS_Q_PT(6,NTSUB)      ! HEXA8 EAS: internal-mode thermal strains M*inv(H)*q at a recovery pt
      REAL(DOUBLE)                    :: FACE_AREA              ! Area of a face of the HEXA where a PLOAD4 pressure acts

                                                                ! Array of all DT values at the grids GRID_DT_ARRAY(i,j) = DT(i,j)
      REAL(DOUBLE)                    :: GRID_DT_ARRAY(ELGP,NTSUB)

      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS)   ! An output from subr ORDER_GAUSS, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: HHH_SH(MAX_ORDER_GAUSS)! An output from subr ORDER_GAUSS, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: INTFAC                 ! An integration factor (constant multiplier for the Gauss integ)
      REAL(DOUBLE)                    :: JAC(3,3)               ! An output from subr JAC3D, called herein. 3 x 3 Jacobian matrix.
      REAL(DOUBLE)                    :: JACI(3,3)              ! An output from subr JAC3D, called herein. 3 x 3 Jacobian inverse.
      REAL(DOUBLE)                    :: KWW(3,3)               ! Portion of differential stiffness matrix

      LOGICAL                         :: HEXA8_EAS              ! True for the HEXA8 bubble element (EAS formulation)

      REAL(DOUBLE)                    :: JAC0_INV(3,3)          ! Jacobian inverse at element centroid (Simo & Rifai J0)
      REAL(DOUBLE)                    :: DETJ0                  ! Determinant of centroid Jacobian
      REAL(DOUBLE)                    :: PSH(ELGP)              ! Output from subr SHP3DH. Shape fcn at Gauss pts SSI, SSJ
      REAL(DOUBLE)                    :: PSIGN                  !
      REAL(DOUBLE)                    :: SIGxx                  ! Normal stress in the elem x  direction
      REAL(DOUBLE)                    :: SIGyy                  ! Normal stress in the elem y  direction
      REAL(DOUBLE)                    :: SIGzz                  ! Normal stress in the elem z  direction
      REAL(DOUBLE)                    :: SIGxy                  ! Shear  stress in the elem xy direction
      REAL(DOUBLE)                    :: SIGyz                  ! Shear  stress in the elem yz direction
      REAL(DOUBLE)                    :: SIGzx                  ! Shear  stress in the elem zx direction
      REAL(DOUBLE)                    :: SSS(MAX_ORDER_GAUSS)   ! An output from subr ORDER_GAUSS, called herein. Gauss abscissa's.
      REAL(DOUBLE)                    :: SSS_SH(MAX_ORDER_GAUSS)! An output from subr ORDER_GAUSS, called herein. Gauss abscissa's.
      REAL(DOUBLE)                    :: SUMB                   ! An intermediate variable used in calc B matrix for reduced integ
      REAL(DOUBLE)                    :: SUMD                   ! An intermediate variable used in calc B matrix for reduced integ
      REAL(DOUBLE)                    :: TBAR(NTSUB)            ! Average elem temperature for each subcase
      REAL(DOUBLE)                    :: TEMP                   ! Temperature to use in PTE calc
      REAL(DOUBLE)                    :: TREF1                  ! TREF(1)
      REAL(DOUBLE)                    :: TGAUSS(1,NTSUB)        ! Temp at a Gauss point for a theral subcase
      REAL(DOUBLE)                    :: VOLUME                 ! 3D element volume
      REAL(DOUBLE)                    :: SSI,SSJ,SSK            ! Isoparametric coordinates of a point.
      REAL(DOUBLE)                    :: M_1DOF(ELGP,ELGP)      ! Consistent mass matrix with 1 DOF per node.



! **********************************************************************************************************************************

      EPS1      = EPSIL(1)
      HEXA8_EAS = (TYPE == 'HEXA8   ') .AND. (BUBBLE == 'Y')

! Calculate ID array

      JJ = 0
      DO I=1,ELGP
         II = 6*(I - 1)
         DO J=1,3
            II = II + 1
            JJ = JJ + 1
            ID(JJ) = II
         ENDDO
      ENDDO

! ALP is first col of ALPVEC
      DO I=1,6
        ALP(I) = ALPVEC(I,1)
      ENDDO

      TREF1 = TREF(1)

! EALP is needed to calculate both PTE and STE2
      EALP = MATMUL(ES,ALP)

! Calc TBAR (used for PTE, STEi)

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         DO J=1,NTSUB
            TBAR(J) = ZERO
            DO K=1,ELGP
               TBAR(J) = TBAR(J) + DT(K,J)
            ENDDO
            TBAR(J) = TBAR(J)/ELGP - TREF1
         ENDDO
      ENDIF

! Calculate volume by Gaussian integration

      IORD_MSG = 'for 3-D solid strains,           input IORD = '
      VOLUME = ZERO
      GAUSS_PT = 0
      CALL ORDER_GAUSS ( IORD, SSS, HHH )
      DO K=1,IORD
         DO J=1,IORD
            DO I=1,IORD
               GAUSS_PT = GAUSS_PT + 1
               CALL SHP3DH ( I, J, K, ELGP, SUBR_NAME, IORD_MSG, IORD, SSS(I), SSS(J), SSS(K), 'Y', PSH, DPSHG )
               CALL JAC3D ( SSS(I), SSS(J), SSS(K), DPSHG, 'Y', JAC, JACI, DETJ(GAUSS_PT) )
               VOLUME = VOLUME + HHH(I)*HHH(J)*HHH(K)*DETJ(GAUSS_PT)
            ENDDO
         ENDDO
      ENDDO

! If VOLUME <= 0, write error and return. det(J) at every Gauss point has already been checked in subr JAC3D relative to the
! element's own dimensions, so no absolute tolerance (which would depend on the element size and model units) is used here.

      IF (VOLUME <= ZERO) THEN
         WRITE(ERR,1925) EID, TYPE, 'VOLUME', VOLUME
         WRITE(F06,1925) EID, TYPE, 'VOLUME', VOLUME
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

! **********************************************************************************************************************************
! Generate the mass matrix for this element.

      IF (OPT(1) == 'Y') THEN

         ! Consistent mass matrix
         ! ME = ∫ N' ρ N det(J) dv

         M_1DOF(:,:) = ZERO
         GAUSS_PT = 0
         CALL ORDER_GAUSS ( IORD, SSS, HHH )
         DO K=1,IORD
            DO J=1,IORD
               DO I=1,IORD
                  GAUSS_PT = GAUSS_PT + 1

                  CALL SHP3DH ( I, J, K, ELGP, SUBR_NAME, IORD_MSG, IORD, SSS(I), SSS(J), SSS(K), 'N', PSH, DPSHG )
                  INTFAC = DETJ(GAUSS_PT)*HHH(I)*HHH(J)*HHH(K)

                  DO L=1,ELGP
                     DO M=1,ELGP
                        M_1DOF(L,M) = M_1DOF(L,M) + PSH(L) * PSH(M) * RHO(1) * INTFAC
                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO


         CALL EXPAND_MASS_DOFS( M_1DOF )

      ENDIF

! **********************************************************************************************************************************
! Generate B matrices. Dimensions 1 and 2 of B store a element B matrix for 1 Gauss point. The 3rd dimension has B for all other
! Gauss points (denoted by IGAUSS index)

       IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         B = ZERO

         CALL ORDER_GAUSS ( IORD, SSS, HHH )
         IORD_MSG = 'for 3-D solid strains,           input IORD = '
         GAUSS_PT = 0
         DO K=1,IORD
            DO J=1,IORD
               DO I=1,IORD
                  GAUSS_PT = GAUSS_PT + 1
                  CALL SHP3DH ( I, J, K, ELGP, SUBR_NAME, IORD_MSG, IORD, SSS(I), SSS(J), SSS(K), 'N', PSH,DPSHG )
                  CALL JAC3D ( SSS(I), SSS(J), SSS(K), DPSHG, 'N', JAC, JACI, DETJ(GAUSS_PT) )
                  DPSHX = MATMUL(JACI, DPSHG)
                  CALL B3D_ISOPARAMETRIC ( DPSHX, GAUSS_PT, I, J, K, 'direct strains', 'Y', BI )
                  DO L=1,6
                     DO M=1,3*ELGP
                        B(L,M,GAUSS_PT) = BI(L,M)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         IF ((RED_INT_SHEAR == 'Y') .AND. (.NOT. HEXA8_EAS)) THEN      ! Reduced shear (the bubble element uses EAS instead)

            IF (IORD == 2) THEN                               ! Use selective substitution for the 2x2x2 integration

               DO L=4,6                                       ! NOTE: the following numbering for the Ki is not the CW or CCW num
!                                                               of element nodes. It is the numbering of Gauss points consistent w/
!                                                               the DO loops (K,J,I) which use Gauss coords SSS(I), SSS(J), SSS(K).

                  IF      (L == 4) THEN                       ! Numbering for selecive  substitution for Gamma-xy
                     K1 = 1
                     K2 = 2
                     K3 = 3
                     K4 = 4
                     K5 = 5
                     K6 = 6
                     K7 = 7
                     K8 = 8
                  ELSE IF (L == 5) THEN                       ! Numbering for selecive  substitution for Gamma-yz
                     K1 = 1
                     K2 = 3
                     K3 = 5
                     K4 = 7
                     K5 = 2
                     K6 = 4
                     K7 = 6
                     K8 = 8
                  ELSE IF (L == 6) THEN                       ! Numbering for selecive  substitution for Gamma-zx
                     K1 = 1
                     K2 = 2
                     K3 = 5
                     K4 = 6
                     K5 = 3
                     K6 = 4
                     K7 = 7
                     K8 = 8
                  ENDIF

                  SUMD = DETJ(K1) + DETJ(K2) + DETJ(K3) + DETJ(K4)
                  DO J=1,3*ELGP                                  ! Jacobian weighted average for shear terms on one face
                     SUMB = DETJ(K1)*B(L,J,K1) + DETJ(K2)*B(L,J,K2) + DETJ(K3)*B(L,J,K3) + DETJ(K4)*B(L,J,K4)
                     B(L,J,K1) = SUMB/SUMD
                     B(L,J,K2) = B(L,J,K1)
                     B(L,J,K3) = B(L,J,K1)
                     B(L,J,K4) = B(L,J,K1)
                  ENDDO

                  SUMD = DETJ(K5) + DETJ(K6) + DETJ(K7) + DETJ(K8)
                  DO J=1,3*ELGP                                  ! Jacobian weighted average for shear terms on opposite face
                     SUMB = DETJ(K5)*B(L,J,K5) + DETJ(K6)*B(L,J,K6) + DETJ(K7)*B(L,J,K7) + DETJ(K8)*B(L,J,K8)
                     B(L,J,K5) = SUMB/SUMD
                     B(L,J,K6) = B(L,J,K5)
                     B(L,J,K7) = B(L,J,K5)
                     B(L,J,K8) = B(L,J,K5)
                  ENDDO
               ENDDO

            ELSE                                              ! IORD = 3 so for reduced integration use 2x2x2 pattern for shear

               IORD_SH = 2
               CALL ORDER_GAUSS ( IORD_SH, SSS_SH, HHH_SH )
               IORD_MSG = 'for reduced integration,            IORD_SH = '
               GAUSS_PT = 0
               DO K=1,IORD_SH
                  DO J=1,IORD_SH
                     DO I=1,IORD_SH
                        GAUSS_PT = GAUSS_PT + 1
                        CALL SHP3DH ( I, J, K, ELGP, SUBR_NAME, IORD_MSG, IORD_SH, SSS_SH(I), SSS_SH(J), SSS_SH(K), 'N', PSH,DPSHG )
                        CALL JAC3D ( SSS_SH(I), SSS_SH(J), SSS_SH(K), DPSHG, 'N', JAC, JACI, DETJ_SH(GAUSS_PT) )
                        DPSHX = MATMUL(JACI, DPSHG)
                        CALL B3D_ISOPARAMETRIC ( DPSHX, GAUSS_PT, I, J, K, 'direct strains', 'Y', BI )
                        DO L=4,6
                           DO M=1,3*ELGP
                              B_SH(L,M,GAUSS_PT) = BI(L,M)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF

         ENDIF

         IF (HEXA8_EAS) THEN
            CALL BUILD_HEXA8_ASSUMED_STRAIN_DATA ( IERR )
            IF (IERR > 0) RETURN
            IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) CALL BUILD_HEXA8_EAS_THERMAL_DATA
         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Calculate element thermal loads.

      IF (OPT(2) == 'Y') THEN

         GRID_DT_ARRAY = DT(:ELGP, :NTSUB)

         DO N=1,NTSUB

            DUM0 = ZERO
            DUM1 = ZERO

            IORD_MSG = 'for 3-D solid strains,           input IORD = '
            GAUSS_PT = 0
            DO K=1,IORD
               DO J=1,IORD
                  DO I=1,IORD
                     GAUSS_PT = GAUSS_PT + 1

                     BI = B(:,:,GAUSS_PT)
                     DUM0 = MATMUL(TRANSPOSE(BI), EALP)
                     INTFAC = DETJ(GAUSS_PT)*HHH(I)*HHH(J)*HHH(K)
                     IF (DEBUG(191) == 0) THEN             ! Use temperatures at Gauss points for PTE
                        CALL SHP3DH ( I, J, K, ELGP, SUBR_NAME, IORD_MSG, IORD, SSS(I), SSS(J), SSS(K), 'N', PSH, DPSHG )
                        TGAUSS(1,:) = MATMUL(PSH, GRID_DT_ARRAY)
                        TEMP = TGAUSS(1,N) - TREF1
                     ELSE                                  ! Use avg element temperature for PTE
                        TEMP = TBAR(N)
                     ENDIF
                     DUM1 = DUM1 + DUM0*TEMP*INTFAC
                  ENDDO
               ENDDO
            ENDDO

            IF (HEXA8_EAS) THEN                         ! Condense the internal-mode thermal load: PTE = PTE - G'*inv(H)*q
               DO L=1,3*ELGP
                  DO M=1,EAS_NUM_MODES
                     DUM1(L) = DUM1(L) - EAS_G(M,L)*EAS_HINV_Q(M,N)
                  ENDDO
               ENDDO
            ENDIF

            PTE(ID(:),N) = DUM1

         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate SEi, STEi, matrices for stress data recovery and BE matrices for strain recovery.

      IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN

        DO STR_PT_NUM=1,9

                                                           ! Isoparametric coordinates of the point
! These coordinates are also defined in SHP3DH. Could try to put them somewhere common.
          SELECT CASE (STR_PT_NUM)
            CASE (1)                                       ! Center
              SSI = 0; SSJ = 0; SSK = 0
            CASE (2)                                       ! Node 1
              SSI = -1; SSJ = -1; SSK = -1
            CASE (3)                                       ! Node 2
              SSI = +1; SSJ = -1; SSK = -1
            CASE (4)                                       ! Node 3
              SSI = +1; SSJ = +1; SSK = -1
            CASE (5)                                       ! Node 4
              SSI = -1; SSJ = +1; SSK = -1
            CASE (6)                                       ! Node 5
              SSI = -1; SSJ = -1; SSK = +1
            CASE (7)                                       ! Node 6
              SSI = +1; SSJ = -1; SSK = +1
            CASE (8)                                       ! Node 7
              SSI = +1; SSJ = +1; SSK = +1
            CASE (9)                                       ! Node 8
              SSI = -1; SSJ = +1; SSK = +1
          END SELECT

          DUM2 = ZERO
                                                           ! Calc SE1,2,3
          IORD_MSG = 'for 3-D solid strains,                      = '
          CALL SHP3DH ( 0, 0, 0, ELGP, SUBR_NAME, IORD_MSG, 1, SSI, SSJ, SSK, 'N', PSH, DPSHG )
          CALL HEXA_EXTRAPOLATED_STRAIN_DISPL ( SSI, SSJ, SSK, BI, EPS_Q_PT )
          DUM2 = MATMUL(ES,BI)

          DO I=1,3                                         ! Stress-displ matrices
            DO J=1,3*ELGP
               SE1(I,ID(J),STR_PT_NUM) = DUM2(I  ,J)
               SE2(I,ID(J),STR_PT_NUM) = DUM2(I+3,J)
            ENDDO
          ENDDO

          DO J=1,NTSUB                                     ! Thermal recovery data (see ELEM_STRE_STRN_ARRAYS): internal-mode
            STE1(1:3,J,STR_PT_NUM) = EPS_Q_PT(1:3,J)       ! thermal strains and the temperature at the point
            STE2(1:3,J,STR_PT_NUM) = EPS_Q_PT(4:6,J)
            STE3(1  ,J,STR_PT_NUM) = SOLID_POINT_TEMPERATURE ( PSH, J )
          ENDDO

          DO I=1,3                                         ! Strain-displ matrices
            DO J=1,3*ELGP
               BE1(I,ID(J),STR_PT_NUM) = BI(I  ,J)
               BE2(I,ID(J),STR_PT_NUM) = BI(I+3,J)
            ENDDO
          ENDDO

        ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE.

      IF (OPT(4) == 'Y') THEN

         DUM3 = ZERO

         RED3 = ((RED_INT_SHEAR == 'Y') .AND. (IORD == 3) .AND. (.NOT. HEXA8_EAS))
         IORD_MSG = ' '
         GAUSS_PT = 0
         DO K=1,IORD
            DO J=1,IORD
               DO I=1,IORD
                  GAUSS_PT = GAUSS_PT + 1
                  BI = B(:,:,GAUSS_PT)
                  IF (RED3) BI(4:6,:) = ZERO                      ! shear terms are integrated separately with 2x2x2 points
                  DUM4 = MATMUL(ES,BI)
                  DUM5 = MATMUL(TRANSPOSE(BI),DUM4)
                  INTFAC = DETJ(GAUSS_PT)*HHH(I)*HHH(J)*HHH(K)
                  DUM3 = DUM3 + DUM5*INTFAC
               ENDDO
            ENDDO
         ENDDO

         IF (RED3) THEN                                       ! Shear strain energy with 2x2x2 Gauss points
            CALL ORDER_GAUSS ( 2, SSS_SH, HHH_SH )
            GAUSS_PT = 0
            DO K=1,2
               DO J=1,2
                  DO I=1,2
                     GAUSS_PT = GAUSS_PT + 1
                     BI = ZERO
                     BI(4:6,:) = B_SH(4:6,:,GAUSS_PT)
                     DUM4 = MATMUL(ES,BI)
                     DUM5 = MATMUL(TRANSPOSE(BI),DUM4)
                     INTFAC = DETJ_SH(GAUSS_PT)*HHH_SH(I)*HHH_SH(J)*HHH_SH(K)
                     DO L=1,3*ELGP
                        DO M=1,3*ELGP
                           DUM3(L,M) = DUM3(L,M) + DUM5(L,M)*INTFAC
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

         IF (HEXA8_EAS) THEN                                  ! Static condensation of the EAS internal modes: K - G'*inv(H)*G
            DUM5 = MATMUL(TRANSPOSE(EAS_G), EAS_HINV_G)
            DUM3 = DUM3 - DUM5
         ENDIF

         DO I=1,3*ELGP
            DO J=1,3*ELGP
               KE(ID(I),ID(J)) = DUM3(I,J)
            ENDDO
         ENDDO

! Set lower triangular portion of KE equal to upper portion

         DO I=2,6*ELGP
            DO J=1,I-1
               KE(I,J) = KE(J,I)
            ENDDO
         ENDDO


      ENDIF

! **********************************************************************************************************************************
! Calculate element PLOAD4 pressure loads.

      IF (OPT(5) == 'Y') THEN                              ! PPE was init to all 0's when EMG was called prior to calling this subr

         GENERATE_PRESS_LOAD = 'N'
         DO J=1,NSUB
            DO I=1,NPLOAD4_3D
               ELIDI  = PLOAD4_3D_DATA(I,2)
               ISCNUM = PLOAD4_3D_DATA(I,3)
               IF ((ELIDI == INT_ELEM_ID) .AND. (ISCNUM == J)) THEN
                  GENERATE_PRESS_LOAD = 'Y'
               ENDIF
            ENDDO
         ENDDO

         IF (GENERATE_PRESS_LOAD == 'Y') THEN
            IF (DEBUG(186) > 0) CALL PRESS_LOAD_DEBUG ( '1' )
            FACE_COUNT = 0
            DO J=1,NSUB

               DO I=1,NPLOAD4_3D

                  ELIDI  = PLOAD4_3D_DATA(I,2)
                  ISCNUM = PLOAD4_3D_DATA(I,3)
                  IF ((ELIDI == INT_ELEM_ID) .AND. (ISCNUM == J)) THEN

                     FACE_COUNT = FACE_COUNT + 1
                     ELIDA = PLOAD4_3D_DATA(I,1)
                     IF (ELIDA /= EID) THEN
                        FATAL_ERR = FATAL_ERR + 1
                        WRITE(ERR,1921) SUBR_NAME, ELIDA, EID
                        WRITE(F06,1921) SUBR_NAME, ELIDA, EID
                        CALL OUTA_HERE ( 'Y' )
                     ENDIF
                     GA = PLOAD4_3D_DATA(I,4)
                     GB = PLOAD4_3D_DATA(I,5)
                     CALL PRESS_FACE_GRIDS ( IERR )
                     IF (IERR > 0) RETURN

                     IF (DEBUG(186) > 0) CALL PRESS_LOAD_DEBUG ( '2' )

                     CALL CALC_FACE_AREA

                     DO K=1,4
                        IF (FACE_NUM == 1) PSIGN =  ONE   ;   IF (FACE_NUM == 1) IDELT = 1
                        IF (FACE_NUM == 2) PSIGN = -ONE   ;   IF (FACE_NUM == 2) IDELT = 1
                        IF (FACE_NUM == 3) PSIGN =  ONE   ;   IF (FACE_NUM == 3) IDELT = 2
                        IF (FACE_NUM == 4) PSIGN = -ONE   ;   IF (FACE_NUM == 4) IDELT = 2
                        IF (FACE_NUM == 5) PSIGN =  ONE   ;   IF (FACE_NUM == 5) IDELT = 3
                        IF (FACE_NUM == 6) PSIGN = -ONE   ;   IF (FACE_NUM == 6) IDELT = 3
                        IDOF = 6*(FACE_NODES(FACE_NUM,K) - 1) + IDELT
                        PPE(IDOF,J) = PPE(IDOF,J) + PSIGN*QUARTER*FACE_AREA*PRESS(3,J)
                     ENDDO

                     IF (DEBUG(186) > 0) CALL PRESS_LOAD_DEBUG ( '3' )

                  ENDIF

               ENDDO

            ENDDO

            IF ((GENERATE_PRESS_LOAD == 'Y') .AND. (DEBUG(186) > 0)) WRITE(F06,*)

         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Calculate linear differential stiffness matrix

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

        CALL ELMDIS
        CALL ELEM_STRE_STRN_ARRAYS ( 1, 'Y' )              ! Stress at the center, in local elem coords like B (not STR_CID coords)

        SIGxx = STRESS(1)
        SIGyy = STRESS(2)
        SIGzz = STRESS(3)
        SIGxy = STRESS(4)
        SIGyz = STRESS(5)
        SIGzx = STRESS(6)


        IF (.FALSE.) THEN                                  ! Use the full G' S G form (ELSE branch)

          KWW(1,1) = SIGyy + SIGzz   ;   KWW(1,2) = -SIGxy   ;   KWW(1,3) = -SIGzx
          KWW(2,2) = SIGxx + SIGzz   ;   KWW(2,3) = -SIGyz
          KWW(3,3) = SIGxx + SIGyy
          KWW(2,1) = KWW(1,2)
          KWW(3,1) = KWW(1,3)
          KWW(3,2) = KWW(2,3)

         DUM3 = ZERO

                                                           ! KED = int( CBAR^T * KWW * CBAR , dV)
          IORD_MSG = ' '
          GAUSS_PT = 0
          DO K=1,IORD
            DO J=1,IORD
              DO I=1,IORD
                GAUSS_PT = GAUSS_PT + 1
                DO L=1,ELGP
                  DPSHX(1,L) = B(1,3*(L-1)+1,GAUSS_PT)     ! Collect shape function derivatives at this Gauss point.
                  DPSHX(2,L) = B(2,3*(L-1)+2,GAUSS_PT)
                  DPSHX(3,L) = B(3,3*(L-1)+3,GAUSS_PT)
                  CBAR(1,3*(L-1)+1) =  ZERO            ; CBAR(1,3*(L-1)+2) = -HALF*DPSHX(3,L) ; CBAR(1,3*(L-1)+3)=  HALF*DPSHX(2,L)
                  CBAR(2,3*(L-1)+1) =  HALF*DPSHX(3,L) ; CBAR(2,3*(L-1)+2) =  ZERO            ; CBAR(2,3*(L-1)+3)= -HALF*DPSHX(1,L)
                  CBAR(3,3*(L-1)+1) = -HALF*DPSHX(2,L) ; CBAR(3,3*(L-1)+2) =  HALF*DPSHX(1,L) ; CBAR(3,3*(L-1)+3)=  ZERO
                ENDDO
                DUM6 = MATMUL(KWW, CBAR)
                DUM5 = MATMUL(TRANSPOSE(CBAR), DUM6)
                INTFAC = DETJ(GAUSS_PT)*HHH(I)*HHH(J)*HHH(K)
                DO L=1,3*ELGP
                  DO M=1,3*ELGP
                    DUM3(L,M) = DUM3(L,M) + DUM5(L,M)*INTFAC
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          DO I=1,3*ELGP
            DO J=1,3*ELGP
              KED(ID(I),ID(J)) = DUM3(I,J)
            ENDDO
          ENDDO

          DO I=2,6*ELGP                                    ! Set lower triangular portion of KE equal to upper portion
            DO J=1,I-1
              KED(I,J) = KED(J,I)
            ENDDO
          ENDDO

        ELSE

          ! Alternative method described in ANSYS's manual
          ! "Theory Reference for the Mechanical APDL and Mechanical Applications"
          ! section 3.4.3 Stress Stiffening Implementation.
          ! The results are about the same as the other method so this probably isn't useful.

          KWW(1,1) = SIGxx   ;   KWW(1,2) = SIGxy   ;   KWW(1,3) = SIGzx
          KWW(2,2) = SIGyy   ;   KWW(2,3) = SIGyz
          KWW(3,3) = SIGzz
          KWW(2,1) = KWW(1,2)
          KWW(3,1) = KWW(1,3)
          KWW(3,2) = KWW(2,3)

         DUM7 = ZERO

                                                           ! KED = int( DPSHX^T * KWW * DPSHX , dV)
          IORD_MSG = ' '
          GAUSS_PT = 0
          DO K=1,IORD
            DO J=1,IORD
               DO I=1,IORD
                  GAUSS_PT = GAUSS_PT + 1
                  DO L=1,ELGP                              ! Collect shape function derivatives at this Gauss point.
                    DPSHX(1,L) = B(1,3*(L-1)+1,GAUSS_PT)
                    DPSHX(2,L) = B(2,3*(L-1)+2,GAUSS_PT)
                    DPSHX(3,L) = B(3,3*(L-1)+3,GAUSS_PT)
                  ENDDO
                  DUM9 = MATMUL(KWW, DPSHX)
                  DUM8 = MATMUL(TRANSPOSE(DPSHX), DUM9)
                  INTFAC = DETJ(GAUSS_PT)*HHH(I)*HHH(J)*HHH(K)
                  DUM7 = DUM7 + DUM8*INTFAC
               ENDDO
            ENDDO
          ENDDO

                                                           ! Initialize KED to zero
         KED = ZERO
                                                           ! Copy submatrix into KED for each DOF.
          DO I=1,ELGP
            DO J=1,ELGP
               KED(ID(3*(I-1)+1),ID(3*(J-1)+1)) = DUM7(I,J)
               KED(ID(3*(I-1)+2),ID(3*(J-1)+2)) = DUM7(I,J)
               KED(ID(3*(I-1)+3),ID(3*(J-1)+3)) = DUM7(I,J)
            ENDDO
          ENDDO


        ENDIF



        IF (DEBUG(178) >= 1) THEN
          WRITE(F06,2001) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
          K = 0
          DO I=1,ELGP
            K = 6*(I-1) + 1
            WRITE(F06,2101) (KED(K  ,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+1,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+2,J),J=1,3*ELGP)
            WRITE(F06,*)
            WRITE(F06,2101) (KED(K+3,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+4,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+5,J),J=1,3*ELGP)
            WRITE(F06,*)
          ENDDO
          WRITE(F06,*)
        ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1921 FORMAT(' *ERROR  1921: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ELEMENT ID READ FROM DATA SET PLOAD4_3D_DATA WAS ',I8,' BUT SUBR IS WORKING ON ELEMENT ',I8)

 1925 FORMAT(' *ERROR  1925: ELEMENT ',I8,', TYPE ',A,', HAS ZERO OR NEGATIVE ',A,' = ',1ES9.1)

 2001 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': KED Differential Stiffness Matrix for ', A,            &
             ' element number ', I8)

 2101 FORMAT(1000(1ES14.6))

89215 format('  L  J      SUMB          SUMD        SUMB/SUMD',/,'  ----------------------------------------------')

89216 format(2i3,7(1es14.6))

91826 format('  I  J  K  Gauss Pt     SSS(I)        SSS(J)        SSS(K)')

91827 format(3i3,i10,3(1es14.6))

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE HEXA_STRAIN_DISPL ( XI, ETA, ZETA, BMAT, EPS_Q )

! Strain-displacement matrix used for stress/strain recovery at isoparametric point (XI,ETA,ZETA): the compatible B including the
! condensed enhanced strains for HEXA8 EAS. The condensed internal modes alpha = inv(H)*(q - G*u) also carry the thermal part
! inv(H)*q, so the element strain at the point is BMAT*u + EPS_Q

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: XI, ETA, ZETA       ! Isoparametric coordinates of the point
      REAL(DOUBLE), INTENT(OUT)       :: BMAT(6,3*ELGP)      ! Strain-displacement matrix
      REAL(DOUBLE), INTENT(OUT)       :: EPS_Q(6,NTSUB)      ! Internal-mode thermal strains M*inv(H)*q (zero unless HEXA8 EAS)

      REAL(DOUBLE)                    :: DETJ_PT             ! det(J) at the point
      REAL(DOUBLE)                    :: DPSHG_PT(3,ELGP)    ! Shape function derivatives w.r.t. xi, eta, zeta
      REAL(DOUBLE)                    :: DPSHX_PT(3,ELGP)    ! Shape function derivatives w.r.t. x, y, z
      REAL(DOUBLE)                    :: EAS_CORR(6,3*ELGP)  ! Condensed enhanced strain-displacement terms
      REAL(DOUBLE)                    :: JAC_PT(3,3)         ! Jacobian at the point
      REAL(DOUBLE)                    :: JACI_PT(3,3)        ! Inverse Jacobian at the point
      REAL(DOUBLE)                    :: PSH_PT(ELGP)        ! Shape functions at the point

! **********************************************************************************************************************************

      CALL SHP3DH ( 0, 0, 0, ELGP, SUBR_NAME, IORD_MSG, 1, XI, ETA, ZETA, 'N', PSH_PT, DPSHG_PT )
      CALL JAC3D ( XI, ETA, ZETA, DPSHG_PT, 'N', JAC_PT, JACI_PT, DETJ_PT )
      DPSHX_PT = MATMUL(JACI_PT, DPSHG_PT)
      CALL B3D_ISOPARAMETRIC ( DPSHX_PT, 0, 1, 1, 1, 'all strains', 'N', BMAT )

      EPS_Q = ZERO
      IF (HEXA8_EAS) THEN
         CALL HEXA8_EAS_MODES ( XI, ETA, ZETA, DETJ_PT, EAS_MODE_PT )
         EAS_CORR = MATMUL(EAS_MODE_PT, EAS_HINV_G)
         BMAT     = BMAT - EAS_CORR
         IF (NTSUB > 0) EPS_Q = MATMUL(EAS_MODE_PT, EAS_HINV_Q)
      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE HEXA_STRAIN_DISPL

! ##################################################################################################################################

      SUBROUTINE HEXA_EXTRAPOLATED_STRAIN_DISPL ( XI, ETA, ZETA, BMAT, EPS_Q )

! Strain-displacement matrix at (XI,ETA,ZETA) obtained by trilinear extrapolation of the recovery B from the 2x2x2 Gauss points
! (xi, eta, zeta = +/- 1/sqrt(3)). Used at points where the Jacobian is not positive. The extrapolation reproduces any trilinear
! strain field, so constant strain (patch test) and the strain of an affine (parallelepiped) HEXA8 are recovered exactly. The
! internal-mode thermal strains are extrapolated the same way.

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: XI, ETA, ZETA       ! Isoparametric coordinates of the point
      REAL(DOUBLE), INTENT(OUT)       :: BMAT(6,3*ELGP)      ! Strain-displacement matrix
      REAL(DOUBLE), INTENT(OUT)       :: EPS_Q(6,NTSUB)      ! Internal-mode thermal strains (see HEXA_STRAIN_DISPL)

      INTEGER(LONG)                   :: IC                  ! Corner (sampling point) index

      REAL(DOUBLE)                    :: BG(6,3*ELGP)        ! Recovery B at one sampling point
      REAL(DOUBLE)                    :: EPS_QG(6,NTSUB)     ! Internal-mode thermal strains at one sampling point
      REAL(DOUBLE)                    :: GSIGN(3)            ! Signs of the sampling point coordinates
      REAL(DOUBLE)                    :: LG                  ! Trilinear interpolation weight of the sampling point
      REAL(DOUBLE)                    :: S3                  ! sqrt(3)

! **********************************************************************************************************************************

      S3    = SQRT(3.0_DOUBLE)
      BMAT  = ZERO
      EPS_Q = ZERO
      DO IC=1,8
         GSIGN(1) = MERGE( ONE, -ONE, (IC == 2) .OR. (IC == 3) .OR. (IC == 6) .OR. (IC == 7))
         GSIGN(2) = MERGE( ONE, -ONE, (IC == 3) .OR. (IC == 4) .OR. (IC == 7) .OR. (IC == 8))
         GSIGN(3) = MERGE( ONE, -ONE, IC >= 5)
         CALL HEXA_STRAIN_DISPL ( GSIGN(1)/S3, GSIGN(2)/S3, GSIGN(3)/S3, BG, EPS_QG )
         LG    = (ONE + S3*GSIGN(1)*XI)*(ONE + S3*GSIGN(2)*ETA)*(ONE + S3*GSIGN(3)*ZETA)/8.0_DOUBLE
         BMAT  = BMAT  + LG*BG
         EPS_Q = EPS_Q + LG*EPS_QG
      ENDDO

! **********************************************************************************************************************************

      END SUBROUTINE HEXA_EXTRAPOLATED_STRAIN_DISPL

! ##################################################################################################################################

      SUBROUTINE BUILD_HEXA8_ASSUMED_STRAIN_DATA ( IERR_OUT )

! Build the HEXA8 Q1/E21 enhanced assumed strain data used by the stiffness, thermal load and recovery calculations:

!   M(xi)  = det(J0)/det(J(xi)) * T0 * Mnat(xi)        enhanced strain modes (Simo & Rifai 1990; Andelfinger & Ramm 1993)
!   G      = integral( M' * ES * B dV )
!   H      = integral( M' * ES * M dV )
!   KE     = integral( B' * ES * B dV ) - G' * inv(H) * G

! Every mode of Mnat has a zero integral over the parent cube, so with the det(J0)/det(J) scaling the enhanced strains are
! orthogonal to constant stress on any element shape (patch test). H is solved with a Cholesky factorization of the diagonally
! equilibrated matrix D*H*D, which removes the dependence of its condition number on element size and aspect ratio.

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: IERR_OUT            ! 0 on success, 1 if H is not positive definite (degenerate element)

      INTEGER(LONG)                   :: IG, JG, KG          ! Gauss loop indices
      INTEGER(LONG)                   :: IGP                 ! Gauss point number
      INTEGER(LONG)                   :: II, JJ_             ! Loop indices
      INTEGER(LONG)                   :: INFO                ! LAPACK status

      REAL(DOUBLE)                    :: DPSHG0(3,ELGP)      ! Shape function derivatives at the centroid
      REAL(DOUBLE)                    :: EAS_H(EAS_NUM_MODES,EAS_NUM_MODES)
                                                             ! Internal-mode stiffness H
      REAL(DOUBLE)                    :: ESM(6,EAS_NUM_MODES)! ES*M at one Gauss point
      REAL(DOUBLE)                    :: G_GP(EAS_NUM_MODES,3*ELGP)
                                                             ! M'*ES*B at one Gauss point
      REAL(DOUBLE)                    :: INTFAC2             ! Integration factor at one Gauss point
      REAL(DOUBLE)                    :: JAC0(3,3)           ! Jacobian at the element centroid
      REAL(DOUBLE)                    :: PSH0(ELGP)          ! Shape functions at the centroid (required by SHP3DH, not used)

! **********************************************************************************************************************************

      IERR_OUT = 0

      CALL SHP3DH ( 0, 0, 0, ELGP, SUBR_NAME, IORD_MSG, 1, ZERO, ZERO, ZERO, 'N', PSH0, DPSHG0 )
      CALL JAC3D ( ZERO, ZERO, ZERO, DPSHG0, 'N', JAC0, JAC0_INV, DETJ0 )
      CALL SOLID_STRAIN_TRANSFORM ( JAC0_INV, EAS_T0 )

      EAS_G = ZERO
      EAS_H = ZERO

      IGP = 0
      DO KG=1,IORD
         DO JG=1,IORD
            DO IG=1,IORD
               IGP = IGP + 1
               INTFAC2 = DETJ(IGP)*HHH(IG)*HHH(JG)*HHH(KG)
               CALL HEXA8_EAS_MODES ( SSS(IG), SSS(JG), SSS(KG), DETJ(IGP), EAS_MODE(:,:,IGP) )
               ESM   = MATMUL(ES, EAS_MODE(:,:,IGP))
               G_GP  = MATMUL(TRANSPOSE(ESM), B(:,:,IGP))
               EAS_G = EAS_G + G_GP*INTFAC2
               EAS_H = EAS_H + MATMUL(TRANSPOSE(EAS_MODE(:,:,IGP)), ESM)*INTFAC2
            ENDDO
         ENDDO
      ENDDO

      INFO = 0
      DO II=1,EAS_NUM_MODES
         IF (EAS_H(II,II) > ZERO) THEN
            EAS_HSCL(II) = ONE/SQRT(EAS_H(II,II))
         ELSE
            INFO = II
            EXIT
         ENDIF
      ENDDO

      IF (INFO == 0) THEN
         DO JJ_=1,EAS_NUM_MODES
            DO II=1,EAS_NUM_MODES
               EAS_HFAC(II,JJ_) = EAS_HSCL(II)*EAS_H(II,JJ_)*EAS_HSCL(JJ_)
            ENDDO
         ENDDO
         CALL DPOTRF ( 'U', EAS_NUM_MODES, EAS_HFAC, EAS_NUM_MODES, INFO )
      ENDIF

      IF (INFO /= 0) THEN
         WRITE(ERR,1966) EID, TYPE, 'EAS INTERNAL MODE STIFFNESS H', INFO
         WRITE(F06,1966) EID, TYPE, 'EAS INTERNAL MODE STIFFNESS H', INFO
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR          = FATAL_ERR + 1
         IERR_OUT           = 1
         RETURN
      ENDIF

      EAS_HINV_G = EAS_G
      CALL HEXA8_EAS_SOLVE ( 3*ELGP, EAS_HINV_G )

! **********************************************************************************************************************************
 1966 FORMAT(' *ERROR  1966: ELEMENT ',I8,', TYPE ',A,': ',A,' IS NOT POSITIVE DEFINITE (CHOLESKY PIVOT ',I3,').',               &
                    /,14X,' THE ELEMENT IS DEGENERATE (COINCIDENT NODES, OR ZERO/NEGATIVE JACOBIAN AT AN INTEGRATION POINT)',     &
                    /,14X,' OR ITS MATERIAL MATRIX IS NOT POSITIVE DEFINITE. CORRECT THE MODEL BEFORE RERUNNING.')

      END SUBROUTINE BUILD_HEXA8_ASSUMED_STRAIN_DATA

! ##################################################################################################################################

      SUBROUTINE BUILD_HEXA8_EAS_THERMAL_DATA

! EAS_HINV_Q = inv(H)*q for the internal-mode thermal load q = integral( M' * ES * alpha*(T - TREF) dV ) of every thermal subcase,
! with T at the Gauss points (the element average with DEBUG(191) > 0, as for PTE). The condensed thermal load PTE and the
! recovered strains B*u + M*inv(H)*q both use it, so the recovered strain is that of the element's own internal modes
! alpha = inv(H)*(q - G*u). It is built for every call that needs either, including recovery calls without OPT(2).

      IMPLICIT NONE

      INTEGER(LONG)                   :: IG, JG, KG          ! Gauss loop indices
      INTEGER(LONG)                   :: IGP                 ! Gauss point number
      INTEGER(LONG)                   :: NT                  ! Thermal subcase index

      REAL(DOUBLE)                    :: DPSHG_GP(3,ELGP)    ! Shape function derivatives at the Gauss point (not used)
      REAL(DOUBLE)                    :: INTFAC2             ! Integration factor at one Gauss point
      REAL(DOUBLE)                    :: MEALP(EAS_NUM_MODES)! M' * ES * alpha at one Gauss point
      REAL(DOUBLE)                    :: PSH_GP(ELGP)        ! Shape functions at the Gauss point
      REAL(DOUBLE)                    :: TEMP_GP             ! Temperature minus TREF at the Gauss point

! **********************************************************************************************************************************

      EAS_HINV_Q = ZERO
      IF (NTSUB <= 0) RETURN

      IGP = 0
      DO KG=1,IORD
         DO JG=1,IORD
            DO IG=1,IORD
               IGP = IGP + 1
               INTFAC2 = DETJ(IGP)*HHH(IG)*HHH(JG)*HHH(KG)
               MEALP   = MATMUL(TRANSPOSE(EAS_MODE(:,:,IGP)), EALP)
               CALL SHP3DH ( IG, JG, KG, ELGP, SUBR_NAME, IORD_MSG, IORD, SSS(IG), SSS(JG), SSS(KG), 'N', PSH_GP, DPSHG_GP )
               DO NT=1,NTSUB
                  TEMP_GP = SOLID_POINT_TEMPERATURE ( PSH_GP, NT )
                  EAS_HINV_Q(:,NT) = EAS_HINV_Q(:,NT) + MEALP*TEMP_GP*INTFAC2
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL HEXA8_EAS_SOLVE ( NTSUB, EAS_HINV_Q )

! **********************************************************************************************************************************

      END SUBROUTINE BUILD_HEXA8_EAS_THERMAL_DATA

! ##################################################################################################################################

      FUNCTION SOLID_POINT_TEMPERATURE ( PSH_PT, NT )

! Temperature minus TREF at a point with shape functions PSH_PT, thermal subcase NT. With DEBUG(191) > 0 the element average
! TBAR is used instead, as for the thermal loads PTE.

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: NT                  ! Thermal subcase index
      REAL(DOUBLE) , INTENT(IN)       :: PSH_PT(ELGP)        ! Shape functions at the point
      REAL(DOUBLE)                    :: SOLID_POINT_TEMPERATURE

! **********************************************************************************************************************************

      IF (DEBUG(191) == 0) THEN
         SOLID_POINT_TEMPERATURE = DOT_PRODUCT(PSH_PT, DT(1:ELGP,NT)) - TREF1
      ELSE
         SOLID_POINT_TEMPERATURE = TBAR(NT)
      ENDIF

! **********************************************************************************************************************************

      END FUNCTION SOLID_POINT_TEMPERATURE

! ##################################################################################################################################

      SUBROUTINE HEXA8_EAS_SOLVE ( NRHS, RHS )

! Overwrite RHS with inv(H)*RHS using the equilibrated Cholesky factor from BUILD_HEXA8_ASSUMED_STRAIN_DATA:
! inv(H) = D * inv(D*H*D) * D

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: NRHS                ! Number of right hand side columns
      REAL(DOUBLE) , INTENT(INOUT)    :: RHS(EAS_NUM_MODES,NRHS)

      INTEGER(LONG)                   :: II                  ! Loop index
      INTEGER(LONG)                   :: INFO                ! LAPACK status

! **********************************************************************************************************************************

      DO II=1,EAS_NUM_MODES
         RHS(II,:) = EAS_HSCL(II)*RHS(II,:)
      ENDDO

      CALL DPOTRS ( 'U', EAS_NUM_MODES, NRHS, EAS_HFAC, EAS_NUM_MODES, RHS, EAS_NUM_MODES, INFO )
      IF (INFO /= 0) THEN                                  ! DPOTRS only fails for invalid arguments
         WRITE(ERR,1921) SUBR_NAME, EID, TYPE, INFO, EAS_NUM_MODES, NRHS
         WRITE(F06,1921) SUBR_NAME, EID, TYPE, INFO, EAS_NUM_MODES, NRHS
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      DO II=1,EAS_NUM_MODES
         RHS(II,:) = EAS_HSCL(II)*RHS(II,:)
      ENDDO

! **********************************************************************************************************************************
 1921 FORMAT(' *ERROR  1921: PROGRAMMING ERROR IN SUBROUTINE ',A,' (HEXA8_EAS_SOLVE) FOR ELEMENT ',I8,', TYPE ',A,              &
                    /,14X,' LAPACK DPOTRS RETURNED INFO = ',I8,' WITH N = ',I8,', NRHS = ',I8)

      END SUBROUTINE HEXA8_EAS_SOLVE

! ##################################################################################################################################

      SUBROUTINE HEXA8_EAS_MODES ( XI, ETA, ZETA, DETJ_PT, MODES )

! Q1/E21 enhanced strain modes at isoparametric point (XI,ETA,ZETA), in the Cartesian element frame and scaled by
! det(J0)/det(J). The natural-frame modes Mnat (engineering strains, ordered xx,yy,zz,xy,yz,zx in xi,eta,zeta) are

!   e_xixi    : xi,  xi*eta,  xi*zeta        e_xieta   : xi,  eta,  xi*zeta,  eta*zeta
!   e_etaeta  : eta, eta*xi,  eta*zeta       e_etazeta : eta, zeta, eta*xi,   zeta*xi
!   e_zetazeta: zeta, zeta*xi, zeta*eta      e_zetaxi  : zeta, xi,  zeta*eta, xi*eta

! Modes 1-9 span the Wilson incompatible-mode (Q1/E9) strains. Each normal mode contains its own coordinate and no mode
! contains a pure function of the other two coordinates, so no compatible strain field of the trilinear element can be
! represented by the modes; the condensed stiffness therefore has only the 6 rigid body zero-energy modes.

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: XI, ETA, ZETA       ! Isoparametric coordinates of the point
      REAL(DOUBLE), INTENT(IN)        :: DETJ_PT             ! det(J) at the point
      REAL(DOUBLE), INTENT(OUT)       :: MODES(6,EAS_NUM_MODES)

      REAL(DOUBLE)                    :: MNAT(6,EAS_NUM_MODES)! Natural-frame modes

! **********************************************************************************************************************************

      MNAT = ZERO

      MNAT(1, 1) = XI         ;   MNAT(2, 2) = ETA        ;   MNAT(3, 3) = ZETA
      MNAT(4, 4) = XI         ;   MNAT(4, 5) = ETA
      MNAT(5, 6) = ETA        ;   MNAT(5, 7) = ZETA
      MNAT(6, 8) = ZETA       ;   MNAT(6, 9) = XI

      MNAT(1,10) = XI*ETA     ;   MNAT(1,11) = XI*ZETA
      MNAT(2,12) = ETA*XI     ;   MNAT(2,13) = ETA*ZETA
      MNAT(3,14) = ZETA*XI    ;   MNAT(3,15) = ZETA*ETA

      MNAT(4,16) = XI*ZETA    ;   MNAT(4,17) = ETA*ZETA
      MNAT(5,18) = ETA*XI     ;   MNAT(5,19) = ZETA*XI
      MNAT(6,20) = ZETA*ETA   ;   MNAT(6,21) = XI*ETA

      MODES = (DETJ0/DETJ_PT)*MATMUL(EAS_T0, MNAT)

! **********************************************************************************************************************************

      END SUBROUTINE HEXA8_EAS_MODES

! ##################################################################################################################################

      SUBROUTINE SOLID_STRAIN_TRANSFORM ( JAC_INV, T )

! Transformation T of engineering strains (xx,yy,zz,xy,yz,zx) from the isoparametric (xi,eta,zeta) frame to the Cartesian frame,
! eps_kl = sum_pq (dxi_p/dx_k)*(dxi_q/dx_l)*eps_pq, with JAC_INV(k,p) = dxi_p/dx_k (inverse of JAC from subr JAC3D).

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: JAC_INV(3,3)        ! Inverse Jacobian
      REAL(DOUBLE), INTENT(OUT)       :: T(6,6)              ! Strain transformation

      INTEGER(LONG), PARAMETER        :: IP(6) = (/ 1, 2, 3, 1, 2, 3 /)
      INTEGER(LONG), PARAMETER        :: JP(6) = (/ 1, 2, 3, 2, 3, 1 /)
      INTEGER(LONG)                   :: R, C, K_, L_, P_, Q_ ! Loop indices and strain component index pairs

! **********************************************************************************************************************************

      DO R=1,6
         K_ = IP(R)   ;   L_ = JP(R)
         DO C=1,6
            P_ = IP(C)   ;   Q_ = JP(C)
            IF (R <= 3) THEN                               ! Normal Cartesian strain
               IF (C <= 3) THEN
                  T(R,C) = JAC_INV(K_,P_)*JAC_INV(K_,P_)
               ELSE
                  T(R,C) = JAC_INV(K_,P_)*JAC_INV(K_,Q_)
               ENDIF
            ELSE                                           ! Engineering shear Cartesian strain
               IF (C <= 3) THEN
                  T(R,C) = TWO*JAC_INV(K_,P_)*JAC_INV(L_,P_)
               ELSE
                  T(R,C) = JAC_INV(K_,P_)*JAC_INV(L_,Q_) + JAC_INV(K_,Q_)*JAC_INV(L_,P_)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

! **********************************************************************************************************************************

      END SUBROUTINE SOLID_STRAIN_TRANSFORM

! ##################################################################################################################################

      SUBROUTINE PRESS_FACE_GRIDS ( IERR )

! (1) Generates 6x4 arrays of the 4 corner grids on the 6 faces of a HEXA
! (2) Determines FACE_NUM, the HEXA face on which both GA and GB (from the PLOAD4 entry) are located

! The corner grid numbering is shown below:

!                   In plane ZI = -1

!                         ETA
!                          |
!                          |
!                          |
!              4 . . . . . . . . . . . 3
!              .           |           .
!              .           |           .
!              .           |           .
!              .           -----------------> XI     NOTE: for HEXA8, only nodes 1-4 are used in plane ZI = -1
!              .                       .             ----
!              .                       .
!              .                       .
!              1 . . . . . . . . . . . 2



!                   In plane ZI = 1

!                         ETA
!                          |
!                          |
!                          |
!              8. . . . .  .. . . . .  7
!              .           |           .
!              .           |           .
!              .           |           .
!              .           -----------------> XI     NOTE: for HEXA8, only nodes 5-8 are used in plane ZI = 1
!              .                       .             ----
!              .                       .
!              .                       .
!              5 . . . . . . . . . . . 6

! The 6 faces of the HEXA are defined as follows:

!               1) Face 1 contains the 4 grids in the plane at  XI = -1: nodes 1-5-8-4
!               2) Face 2 contains the 4 gtids in the plane at  XI =  1: nodes 2-3-7-6
!               3) Face 3 contains the 4 gtids in the plane at ETA = -1: nodes 1-2-6-5
!               4) Face 4 contains the 4 gtids in the plane at ETA =  1: nodes 3-4-8-7
!               5) Face 5 contains the 4 gtids in the plane at ZET = -1: nodes 1-4-3-2
!               6) Face 6 contains the 4 gtids in the plane at ZET =  1: nodes 5-6-7-8


      USE PENTIUM_II_KIND
      USE IOUNT1, ONLY                :  ERR, F06

      IMPLICIT NONE

      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Local error count

! **********************************************************************************************************************************
! Row i of FACE_NODES has the corner node numbers for face i.
! Row i of FACE_AGRDS has the corner grid numbers (actual) for face i.

      FACE_NODES(1,1) =       4    ;   FACE_NODES(1,2) =       1    ;   FACE_NODES(1,3) =       5    ;   FACE_NODES(1,4) =       8
      FACE_NODES(2,1) =       2    ;   FACE_NODES(2,2) =       3    ;   FACE_NODES(2,3) =       7    ;   FACE_NODES(2,4) =       6
      FACE_NODES(3,1) =       1    ;   FACE_NODES(3,2) =       2    ;   FACE_NODES(3,3) =       6    ;   FACE_NODES(3,4) =       5
      FACE_NODES(4,1) =       3    ;   FACE_NODES(4,2) =       4    ;   FACE_NODES(4,3) =       8    ;   FACE_NODES(4,4) =       7
      FACE_NODES(5,1) =       2    ;   FACE_NODES(5,2) =       1    ;   FACE_NODES(5,3) =       4    ;   FACE_NODES(5,4) =       3
      FACE_NODES(6,1) =       5    ;   FACE_NODES(6,2) =       6    ;   FACE_NODES(6,3) =       7    ;   FACE_NODES(6,4) =       8

      FACE_AGRDS(1,1) = AGRID(4)   ;   FACE_AGRDS(1,2) = AGRID(1)   ;   FACE_AGRDS(1,3) = AGRID(5)   ;   FACE_AGRDS(1,4) = AGRID(8)
      FACE_AGRDS(2,1) = AGRID(2)   ;   FACE_AGRDS(2,2) = AGRID(3)   ;   FACE_AGRDS(2,3) = AGRID(7)   ;   FACE_AGRDS(2,4) = AGRID(6)
      FACE_AGRDS(3,1) = AGRID(1)   ;   FACE_AGRDS(3,2) = AGRID(2)   ;   FACE_AGRDS(3,3) = AGRID(6)   ;   FACE_AGRDS(3,4) = AGRID(5)
      FACE_AGRDS(4,1) = AGRID(3)   ;   FACE_AGRDS(4,2) = AGRID(4)   ;   FACE_AGRDS(4,3) = AGRID(8)   ;   FACE_AGRDS(4,4) = AGRID(7)
      FACE_AGRDS(5,1) = AGRID(2)   ;   FACE_AGRDS(5,2) = AGRID(1)   ;   FACE_AGRDS(5,3) = AGRID(4)   ;   FACE_AGRDS(5,4) = AGRID(3)
      FACE_AGRDS(6,1) = AGRID(5)   ;   FACE_AGRDS(6,2) = AGRID(6)   ;   FACE_AGRDS(6,3) = AGRID(7)   ;   FACE_AGRDS(6,4) = AGRID(8)


      IERR     = 0
      FACE_NUM = 0

! If GA and GB are grids that are both on face i then set FACE_NUM to i

      IF ((GA == FACE_AGRDS(1,1)) .OR. (GA == FACE_AGRDS(1,2)) .OR. (GA == FACE_AGRDS(1,3)) .OR. (GA == FACE_AGRDS(1,4))) THEN
         IF   ((GB == FACE_AGRDS(1,1)) .OR. (GB == FACE_AGRDS(1,2)) .OR. (GB == FACE_AGRDS(1,3)) .OR. (GB == FACE_AGRDS(1,4))) THEN
            FACE_NUM = 1
         ENDIF
      ENDIF

      IF ((GA == FACE_AGRDS(2,1)) .OR. (GA == FACE_AGRDS(2,2)) .OR. (GA == FACE_AGRDS(2,3)) .OR. (GA == FACE_AGRDS(2,4))) THEN
         IF   ((GB == FACE_AGRDS(2,1)) .OR. (GB == FACE_AGRDS(2,2)) .OR. (GB == FACE_AGRDS(2,3)) .OR. (GB == FACE_AGRDS(2,4))) THEN
            FACE_NUM = 2
         ENDIF
      ENDIF

      IF ((GA == FACE_AGRDS(3,1)) .OR. (GA == FACE_AGRDS(3,2)) .OR. (GA == FACE_AGRDS(3,3)) .OR. (GA == FACE_AGRDS(3,4))) THEN
         IF   ((GB == FACE_AGRDS(3,1)) .OR. (GB == FACE_AGRDS(3,2)) .OR. (GB == FACE_AGRDS(3,3)) .OR. (GB == FACE_AGRDS(3,4))) THEN
            FACE_NUM = 3
         ENDIF
      ENDIF

      IF ((GA == FACE_AGRDS(4,1)) .OR. (GA == FACE_AGRDS(4,2)) .OR. (GA == FACE_AGRDS(4,3)) .OR. (GA == FACE_AGRDS(4,4))) THEN
         IF   ((GB == FACE_AGRDS(4,1)) .OR. (GB == FACE_AGRDS(4,2)) .OR. (GB == FACE_AGRDS(4,3)) .OR. (GB == FACE_AGRDS(4,4))) THEN
            FACE_NUM = 4
         ENDIF
      ENDIF

      IF ((GA == FACE_AGRDS(5,1)) .OR. (GA == FACE_AGRDS(5,2)) .OR. (GA == FACE_AGRDS(5,3)) .OR. (GA == FACE_AGRDS(5,4))) THEN
         IF   ((GB == FACE_AGRDS(5,1)) .OR. (GB == FACE_AGRDS(5,2)) .OR. (GB == FACE_AGRDS(5,3)) .OR. (GB == FACE_AGRDS(5,4))) THEN
            FACE_NUM = 5
         ENDIF
      ENDIF

      IF ((GA == FACE_AGRDS(6,1)) .OR. (GA == FACE_AGRDS(6,2)) .OR. (GA == FACE_AGRDS(6,3)) .OR. (GA == FACE_AGRDS(6,4))) THEN
         IF   ((GB == FACE_AGRDS(6,1)) .OR. (GB == FACE_AGRDS(6,2)) .OR. (GB == FACE_AGRDS(6,3)) .OR. (GB == FACE_AGRDS(6,4))) THEN
            FACE_NUM = 6
         ENDIF
      ENDIF

      IF (FACE_NUM == 0) THEN
         IERR = IERR + 1
         FATAL_ERR = FATAL_ERR + 1
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         WRITE(ERR,1961) TYPE, EID, GA, GB
         WRITE(F06,1961) TYPE, EID, GA, GB
         RETURN
      ENDIF

! **********************************************************************************************************************************
 1961 FORMAT(' *ERROR  1961: A PLOAD4 BULK DATA ENTRY FOR ',A,' ELEMENT ',I8,' REFERENCES GRIDS ',I8,' AND ',I8                    &
                    ,/,14X,' AT LEAST ONE OF THESE GRIDS IS NOT CONNECTED TO THIS ELEMENT')

! **********************************************************************************************************************************

      END SUBROUTINE PRESS_FACE_GRIDS

! ##################################################################################################################################

      SUBROUTINE CALC_FACE_AREA

      IMPLICIT NONE

      INTEGER(LONG)                   :: II,JJ                ! DO loop indices

      REAL(DOUBLE)                    :: DETJ                 ! An output from subr JAC2D4, called herein. Determinant of JAC
      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS) ! An output from subr ORDER, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: JAC(2,2)             ! An output from subr JAC2D4, called herein. 2 x 2 Jacobian matrix.
      REAL(DOUBLE)                    :: SSS(MAX_ORDER_GAUSS) ! An output from subr ORDER, called herein. Gauss abscissa's.
      REAL(DOUBLE)                    :: XSD(4)               ! Diffs in x coords of HEXA sides in local coords
      REAL(DOUBLE)                    :: YSD(4)               ! Diffs in y coords of HEXA sides in local coords
      REAL(DOUBLE)                    :: ZSD(4)               ! Diffs in z coords of HEXA sides in local coords

! **********************************************************************************************************************************
! Calculate side diffs

      XSD(1) = XEL(FACE_NODES(FACE_NUM,1),1) - XEL(FACE_NODES(FACE_NUM,2),1)
      XSD(2) = XEL(FACE_NODES(FACE_NUM,2),1) - XEL(FACE_NODES(FACE_NUM,3),1)
      XSD(3) = XEL(FACE_NODES(FACE_NUM,3),1) - XEL(FACE_NODES(FACE_NUM,4),1)
      XSD(4) = XEL(FACE_NODES(FACE_NUM,4),1) - XEL(FACE_NODES(FACE_NUM,1),1)

      YSD(1) = XEL(FACE_NODES(FACE_NUM,1),2) - XEL(FACE_NODES(FACE_NUM,2),2)
      YSD(2) = XEL(FACE_NODES(FACE_NUM,2),2) - XEL(FACE_NODES(FACE_NUM,3),2)
      YSD(3) = XEL(FACE_NODES(FACE_NUM,3),2) - XEL(FACE_NODES(FACE_NUM,4),2)
      YSD(4) = XEL(FACE_NODES(FACE_NUM,4),2) - XEL(FACE_NODES(FACE_NUM,1),2)

      ZSD(1) = XEL(FACE_NODES(FACE_NUM,1),3) - XEL(FACE_NODES(FACE_NUM,2),3)
      ZSD(2) = XEL(FACE_NODES(FACE_NUM,2),3) - XEL(FACE_NODES(FACE_NUM,3),3)
      ZSD(3) = XEL(FACE_NODES(FACE_NUM,3),3) - XEL(FACE_NODES(FACE_NUM,4),3)
      ZSD(4) = XEL(FACE_NODES(FACE_NUM,4),3) - XEL(FACE_NODES(FACE_NUM,1),3)


! Calculate area by Gaussian integration

      FACE_AREA = ZERO
      CALL ORDER_GAUSS ( 2, SSS, HHH )
      DO II=1,2
         DO JJ=1,2
            IF      (FACE_NUM == 1) THEN
               CALL JAC2D ( SSS(II), SSS(JJ),-YSD, ZSD, 'N', JAC, JACI, DETJ )
            ELSE IF (FACE_NUM == 2) THEN
               CALL JAC2D ( SSS(II), SSS(JJ), YSD, ZSD, 'N', JAC, JACI, DETJ )
            ELSE IF (FACE_NUM == 3) THEN
               CALL JAC2D ( SSS(II), SSS(JJ), XSD, ZSD, 'N', JAC, JACI, DETJ )
            ELSE IF (FACE_NUM == 4) THEN
               CALL JAC2D ( SSS(II), SSS(JJ),-XSD, ZSD, 'N', JAC, JACI, DETJ )
            ELSE IF (FACE_NUM == 5) THEN
               CALL JAC2D ( SSS(II), SSS(JJ),-XSD, YSD, 'N', JAC, JACI, DETJ )
            ELSE IF (FACE_NUM == 6) THEN
               CALL JAC2D ( SSS(II), SSS(JJ), XSD, YSD, 'N', JAC, JACI, DETJ )
            ENDIF
            FACE_AREA = FACE_AREA + HHH(II)*HHH(JJ)*DETJ
         ENDDO
      ENDDO

! **********************************************************************************************************************************




! **********************************************************************************************************************************

      END SUBROUTINE CALC_FACE_AREA

! ##################################################################################################################################

      SUBROUTINE PRESS_LOAD_DEBUG ( WHAT )

      USE PENTIUM_II_KIND
      USE IOUNT1, ONLY                :  ERR, F06

      IMPLICIT NONE

      CHARACTER( 1*BYTE)              :: WHAT              ! What debug info to print out

      INTEGER(LONG)                   :: KK                ! Counter

! **********************************************************************************************************************************
      IF      (WHAT == '1') THEN
         WRITE(F06,101) TYPE, EID

      ELSE IF (WHAT == '2') THEN
         WRITE(F06,201) FACE_COUNT, ISCNUM, FACE_NUM
         WRITE(F06,202) I, EID, INT_ELEM_ID, ISCNUM, GA, GB
         WRITE(F06,203)
         DO KK=1,4
            WRITE(F06,204) KK, KK, FACE_AGRDS(FACE_NUM,KK), FACE_NODES(FACE_NUM,KK)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == '3') THEN
         WRITE(F06,301) FACE_AREA, PRESS(3,J)
         WRITE(F06,302)
         DO KK=1,6*ELGP
            IF (PPE(KK,J) /= 0.D0) THEN
               WRITE(F06,303) KK, J, PPE(KK,J)
            ENDIF
         ENDDO
         WRITE(F06,*)

      ENDIF

! **********************************************************************************************************************************
  101 FORMAT(1X,A,' element ',I8,' PLOAD4 pressure data generation (subr HEXA)',/,                                                 &
             ' =====================================================================',/)

  201 FORMAT( 4X,'(',I1,') Load data and generation for:',//,                                                                      &
             14X,'+ + + + + + + + + + + + + + + + + + + + + +'                ,/,                                                  &
             14X,'+  Internal subcase (S/C) number',I8,  '  +'                ,/,                                                  &
             14X,'+       element face number     ',I8,  '  +'                ,/,                                                  &
             14X,'+ + + + + + + + + + + + + + + + + + + + + +'                ,/)

  202 FORMAT(11X,'PLOAD4_3D_DATA array, row ',I1,':'             ,/,                                                               &
             11X,'   Element ID (EID)................ =',I9,/,                                                                     &
             11X,'   Internal elem ID (INT_ELEM_ID).. =',I9,/,                                                                     &
             11X,'   Internal S/C number (ISCNUM).... =',I9,/,                                                                     &
             11X,'   PLOAD4 field 8 grid (GA)........ =',I9,/,                                                                     &
             11X,'   PLOAD4 field 8 grid (GA)........ =',I9,/)

  203 FORMAT(11X,'Element actual grid and internal node numbers:')

  204 FORMAT(11X,'   Actual grid(',I1,'), internal node(',I1,') = ',I8,I3)

  301 FORMAT(11X,'Face area and average pressure:',/,                                                                              &
             11X,'   Face area....................... = ',1ES13.6,/,                                                               &
             11X,'   Avg face pressure............... = ',1ES13.6,/)

  302 FORMAT(11X,'Nonzero PLOAD4 pressure loads:')

  303 FORMAT(11X,'   PPE(',I2,',',I2,')...................... = ',1ES13.6)

! **********************************************************************************************************************************

      END SUBROUTINE PRESS_LOAD_DEBUG

      END SUBROUTINE HEXA


      SUBROUTINE PENTA ( OPT, INT_ELEM_ID, IORD_IJ, IORD_K, RED_INT_SHEAR, WRITE_WARN )

! Isoparametric pentahedron solid element (6 or 15 nodes and with full or reduced gaussian integration)

! Subroutine calculates:

!  1) ME        = element mass matrix                  , if OPT(1) = 'Y'
!  2) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  3) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  4) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  5) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_TRIA, MAX_ORDER_GAUSS, NTSUB
      USE CONSTANTS_1, ONLY           :  HALF, THIRD, ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, BE1, BE2, DT, EID, ELGP, NUM_EMG_FATAL_ERRS, ES, KE, KED, ME, PTE, RHO,           &
                                         SE1, SE2, STE3, STRESS, TREF, TYPE

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE ELEMENT_RECOVERY_SUPPORT, ONLY:  ELEM_STRE_STRN_ARRAYS
      USE MASS_DOF_EXPANSION, ONLY    :  EXPAND_MASS_DOFS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'PENTA'
      CHARACTER( 1*BYTE), INTENT(IN)  :: RED_INT_SHEAR           ! If 'Y', use Gaussian weighted avg of B matrices for shear terms
      CHARACTER( 1*BYTE), INTENT(IN)  :: OPT(6)                  ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN              ! If 'Y" write warning messages, otherwise do not
      CHARACTER(46*BYTE)              :: IORD_MSG                ! Character name of an integration order (used for debug output)

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID             ! Internal element ID
      INTEGER(LONG), INTENT(IN)       :: IORD_IJ                 ! Integration order in the triangular plane
      INTEGER(LONG), INTENT(IN)       :: IORD_K                  ! Integration order in Z direction
      INTEGER(LONG)                   :: IORD_IJ_SH              ! Integration order in the triang plane for red integr for shear
      INTEGER(LONG)                   :: IORD_K_SH               ! Integration order in Z direction for red integr for shear
      INTEGER(LONG)                   :: GAUSS_PT                ! Gauss point number (used for DEBUG output in subr SHP3DP
      INTEGER(LONG)                   :: I,J,IJ,K,L,M,N          ! DO loop indices
      INTEGER(LONG)                   :: II,JJ                   ! Counters
      INTEGER(LONG)                   :: ID(3*ELGP)              ! Array which shows equivalence of DOF's in virgin element with the
   !                                                             6 DOF/grid of the final element stiffness matrix

      INTEGER(LONG)                   :: STR_PT_NUM              ! Stress point number. 1 is center, 2+ are element nodes 1+.

      REAL(DOUBLE)                    :: ALP(6)                  ! First col of ALPVEC
                                                                 ! Strain-displ matrix for this element for all Gauss points
      REAL(DOUBLE)                    :: B(6,3*ELGP,IORD_IJ*IORD_K)

      REAL(DOUBLE)                    :: BI(6,3*ELGP)            ! Strain-displ matrix for this element for one Gauss point
      REAL(DOUBLE)                    :: CBAR(3,3*ELGP)          ! Derivatives of shape fcns wrt x,y,z used in diff stiff matrix
!                                                                  (contains terms from DPSHX matrices for each grid of the PENTA)

      REAL(DOUBLE)                    :: DETJ(IORD_IJ*IORD_K)    ! Determinant of JAC for all Gauss points
      REAL(DOUBLE)                    :: DUM_DETJ
      REAL(DOUBLE)                    :: DPSHG(3,ELGP)           ! Output from subr SHP3DP. Derivs of PSH wrt elem isopar coords.
      REAL(DOUBLE)                    :: DPSHX(3,ELGP)           ! Derivatives of PSH wrt elem x, y coords.
      REAL(DOUBLE)                    :: DUM0(3*ELGP)            ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM1(3*ELGP)            ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM2(6,3*ELGP)          ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM3(3*ELGP,3*ELGP)     ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM4(6,3*ELGP)          ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM5(3*ELGP,3*ELGP)     ! Intermediate matrix used in solving for KE elem matrices
      REAL(DOUBLE)                    :: DUM6(3,3*ELGP)          ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: EALP(6)                 ! Interm var used in calc PTE therm loads & STEi stress coeffs

                                                                 ! Array of all DT values at the grids GRID_DT_ARRAY(i,j) = DT(i,j)
      REAL(DOUBLE)                    :: GRID_DT_ARRAY(ELGP,NTSUB)

      REAL(DOUBLE)                    :: HH_IJ(MAX_ORDER_TRIA)   ! An output from subr ORDER_TRIA, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: HH_K(MAX_ORDER_GAUSS)   ! An output from subr ORDER_TRIA, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: HH_IJ_SH(MAX_ORDER_TRIA)! An output from subr ORDER_TRIA, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: HH_K_SH(MAX_ORDER_GAUSS)! An output from subr ORDER_TRIA, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: INTFAC                  ! An integration factor (constant multiplier for the Gauss integ)
      REAL(DOUBLE)                    :: JAC(3,3)                ! An output from subr JAC3D, called herein. 3 x 3 Jacobian matrix.
      REAL(DOUBLE)                    :: JACI(3,3)               ! An output from subr JAC3D, called herein. 3 x 3 Jacobian inverse.
      REAL(DOUBLE)                    :: KWW(3,3)                ! Portion of differential stiffness matrix
      REAL(DOUBLE)                    :: DUM8(ELGP,ELGP)         ! G' S G for one Gauss point (differential stiffness)
      REAL(DOUBLE)                    :: DUM9(3,ELGP)            ! S G for one Gauss point (differential stiffness)
      REAL(DOUBLE)                    :: PSH(ELGP)               ! Output from subr SHP3DP. Shape fcn at Gauss pts SSI, SSJ
      REAL(DOUBLE)                    :: SIGxx                   ! Normal stress in the elem x  direction
      REAL(DOUBLE)                    :: SIGyy                   ! Normal stress in the elem y  direction
      REAL(DOUBLE)                    :: SIGzz                   ! Normal stress in the elem z  direction
      REAL(DOUBLE)                    :: SIGxy                   ! Shear  stress in the elem xy direction
      REAL(DOUBLE)                    :: SIGyz                   ! Shear  stress in the elem yz direction
      REAL(DOUBLE)                    :: SIGzx                   ! Shear  stress in the elem zx direction
      REAL(DOUBLE)                    :: SS_I(MAX_ORDER_TRIA)    ! An output from subr ORDER_TRIA, called herein. Gauss abscissa's
      REAL(DOUBLE)                    :: SS_J(MAX_ORDER_TRIA)    ! An output from subr ORDER_TRIA, called herein. Gauss abscissa's
      REAL(DOUBLE)                    :: SS_K(MAX_ORDER_GAUSS)   ! An output from subr ORDER_GAUSS, called herein. Gauss abscissa's
      REAL(DOUBLE)                    :: SS_I_SH(MAX_ORDER_TRIA) ! An output from subr ORDER_TRIA, called herein. Gauss abscissa's
      REAL(DOUBLE)                    :: SS_J_SH(MAX_ORDER_TRIA) ! An output from subr ORDER_TRIA, called herein. Gauss abscissa's
      REAL(DOUBLE)                    :: SS_K_SH(MAX_ORDER_GAUSS)! An output from subr ORDER_GAUSS, called herein. Gauss abscissa's
      REAL(DOUBLE)                    :: SUMB                    ! An intermediate variable used in calc B matrix for reduced integ
      REAL(DOUBLE)                    :: SUMD                    ! An intermediate variable used in calc B matrix for reduced integ
      REAL(DOUBLE)                    :: TBAR(NTSUB)             ! Average elem temperature for each subcase
      REAL(DOUBLE)                    :: TEMP                    ! Temperature to use in PTE calc
      REAL(DOUBLE)                    :: TREF1                   ! TREF(1)
      REAL(DOUBLE)                    :: TGAUSS(1,NTSUB)         ! Temp at a Gauss point for a theral subcase
      REAL(DOUBLE)                    :: VOLUME                  ! 3D element volume
      REAL(DOUBLE)                    :: SSI,SSJ,SSK             ! Isoparametric coordinates of a point.
      REAL(DOUBLE)                    :: M_1DOF(ELGP,ELGP)      ! Consistent mass matrix with 1 DOF per node.



! **********************************************************************************************************************************

! Calculate volume by Gaussian integration

      CALL ORDER_TRIA  ( IORD_IJ, SS_I, SS_J, HH_IJ )
      CALL ORDER_GAUSS ( IORD_K , SS_K , HH_K  )

! Calculate volume by Gaussian integration

      IORD_MSG = 'for triangular-cross section x z direction  = '
      VOLUME = ZERO
      GAUSS_PT = 0
      DO K=1,IORD_K
         DO IJ=1,IORD_IJ
            GAUSS_PT = GAUSS_PT + 1
            CALL SHP3DP ( I, J, K, ELGP, SUBR_NAME, IORD_MSG, IORD_IJ, IORD_K, SS_I(IJ), SS_J(IJ), SS_K(K), 'Y', PSH, DPSHG )
            CALL JAC3D ( SS_I(IJ), SS_J(IJ), SS_K(K), DPSHG, 'Y', JAC, JACI, DETJ(GAUSS_PT) )
            VOLUME = VOLUME + HH_IJ(IJ)*HH_K(K)*DETJ(GAUSS_PT)
         ENDDO
      ENDDO

      IF (VOLUME <= ZERO) THEN                             ! If VOLUME <= 0, write error and return (det(J) is checked in JAC3D)
         WRITE(ERR,1925) EID, TYPE, 'VOLUME', VOLUME
         WRITE(F06,1925) EID, TYPE, 'VOLUME', VOLUME
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

! Calculate ID array

      JJ = 0
      DO I=1,ELGP
         II = 6*(I - 1)
         DO J=1,3
            II = II + 1
            JJ = JJ + 1
            ID(JJ) = II
         ENDDO
      ENDDO

! ALP is first col of ALPVEC
      DO I=1,6
        ALP(I) = ALPVEC(I,1)
      ENDDO

      TREF1 = TREF(1)

! EALP is needed to calculate both PTE and STE2

      CALL MATMULT_FFF ( ES, ALP, 6, 6, 1, EALP )

! Calc TBAR (used for PTE, STEi)

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         DO J=1,NTSUB
            TBAR(J) = ZERO
            DO K=1,ELGP
               TBAR(J) = TBAR(J) + DT(K,J)
            ENDDO
            TBAR(J) = TBAR(J)/ELGP - TREF1
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Generate the mass matrix for this element.

      IF (OPT(1) == 'Y') THEN

         ! Consistent mass matrix
         ! ME = ∫ N' ρ N det(J) dv

         M_1DOF(:,:) = ZERO
         GAUSS_PT = 0
         DO K=1,IORD_K
            DO IJ=1,IORD_IJ
               GAUSS_PT = GAUSS_PT + 1

               CALL SHP3DP ( I, J, K, ELGP, SUBR_NAME, IORD_MSG, IORD_IJ, IORD_K, SS_I(IJ), SS_J(IJ), SS_K(K), 'N', PSH, DPSHG )
               INTFAC = DETJ(GAUSS_PT)*HH_IJ(IJ)*HH_K(K)

               DO L=1,ELGP
                  DO M=1,ELGP
                     M_1DOF(L,M) = M_1DOF(L,M) + PSH(L) * PSH(M) * RHO(1) * INTFAC
                  ENDDO
               ENDDO
            ENDDO
         ENDDO


         CALL EXPAND_MASS_DOFS( M_1DOF )

      ENDIF

! **********************************************************************************************************************************
! Generate B matrices. Dimensions 1 and 2 of B store a element B matrix for 1 Gauss point. The 3rd dimension has B for all other
! Gauss points (denoted by IGAUSS index)

opt234:IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         B = ZERO

         IORD_MSG = 'for triangular-cross section x z direction  = '
         GAUSS_PT = 0
         DO K=1,IORD_K
            DO IJ=1,IORD_IJ
               GAUSS_PT = GAUSS_PT + 1
               CALL SHP3DP ( IJ, IJ, K, ELGP, SUBR_NAME, IORD_MSG, IORD_IJ, IORD_K, SS_I(IJ), SS_J(IJ), SS_K(K), 'N', PSH, DPSHG )
               CALL JAC3D ( SS_I(IJ), SS_J(IJ), SS_K(K), DPSHG, 'N', JAC, JACI, DETJ(GAUSS_PT) )
               CALL MATMULT_FFF ( JACI, DPSHG, 3, 3, ELGP, DPSHX )
               CALL B3D_ISOPARAMETRIC ( DPSHX, GAUSS_PT, IJ, IJ, K, 'direct strains', 'Y', BI )
               DO L=1,6
                  DO M=1,3*ELGP
                     B(L,M,GAUSS_PT) = BI(L,M)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         IF (RED_INT_SHEAR == 'Y') THEN

            IF (IORD_IJ == 3) THEN                            ! Full integ is 2x3 so calc B using selective substit for shear terms
!                                                               -------------------------------------------------------------------

               DO L=4,6                                       ! NOTE: the following numbering for the Ki is not the CW or CCW num
!                                                               of element nodes. It is the numbering of Gauss points consistent w/
!                                                               the DO loops (K,J,I) which use Gauss coords SSS(I), SSS(J), SSS(K).
                  IF      (L == 4) THEN                       ! Numbering for selecive  substitution for Gamma-xy

                     SUMD = DETJ(1) + DETJ(2) + DETJ(3)
                     DO J=1,3*ELGP                            !      Jacobian weighted average for xy shear terms on one z face
                        SUMB = DETJ(1)*B(L,J,1) + DETJ(2)*B(L,J,2) + DETJ(3)*B(L,J,3)
                        B(L,J,1) = SUMB/SUMD
                        B(L,J,2) = B(L,J,1)
                        B(L,J,3) = B(L,J,1)
                     ENDDO

                     SUMD = DETJ(4) + DETJ(5) + DETJ(6)
                     DO J=1,3*ELGP                            !      Jacobian weighted average for xy shear terms on other z face
                        SUMB = DETJ(4)*B(L,J,4) + DETJ(5)*B(L,J,5) + DETJ(6)*B(L,J,6)
                        B(L,J,4) = SUMB/SUMD
                        B(L,J,5) = B(L,J,4)
                        B(L,J,6) = B(L,J,4)
                     ENDDO

                  ELSE IF (L == 5) THEN                       ! Numbering for selecive  substitution for Gamma-yz

                     SUMD = DETJ(1) + DETJ(2) + DETJ(4) + DETJ(5)
                     DO J=1,3*ELGP                            !      Jacobian weighted average for yz shear terms
                        SUMB = DETJ(1)*B(L,J,1) + DETJ(2)*B(L,J,2) + DETJ(4)*B(L,J,4) + DETJ(5)*B(L,J,5)
                        B(L,J,1) = SUMB/SUMD
                        B(L,J,2) = B(L,J,1)
                        B(L,J,4) = B(L,J,1)
                        B(L,J,5) = B(L,J,1)
                     ENDDO

                     SUMD = DETJ(3) + DETJ(6)
                     DO J=1,3*ELGP                            !      Jacobian weighted average for yz shear terms
                        SUMB = DETJ(3)*B(L,J,3) + DETJ(6)*B(L,J,6)
                        B(L,J,3) = SUMB/SUMD
                        B(L,J,6) = B(L,J,3)
                     ENDDO

                  ELSE IF (L == 6) THEN                       ! Numbering for selecive  substitution for Gamma-zx

                     SUMD = DETJ(2) + DETJ(3) + DETJ(5) + DETJ(6)
                     DO J=1,3*ELGP                            !      Jacobian weighted average for xz shear terms
                        SUMB = DETJ(2)*B(L,J,2) + DETJ(3)*B(L,J,3) + DETJ(5)*B(L,J,5) + DETJ(6)*B(L,J,6)
                        B(L,J,2) = SUMB/SUMD
                        B(L,J,3) = B(L,J,2)
                        B(L,J,5) = B(L,J,2)
                        B(L,J,6) = B(L,J,2)
                     ENDDO

                     SUMD = DETJ(1) + DETJ(4)
                     DO J=1,3*ELGP                            !      Jacobian weighted average for yz shear terms
                        SUMB = DETJ(1)*B(L,J,1) + DETJ(4)*B(L,J,4)
                        B(L,J,1) = SUMB/SUMD
                        B(L,J,4) = B(L,J,1)
                     ENDDO

                  ENDIF

               ENDDO

            ELSE                                              ! IORD_IJ = 7, so calc B using 2x3 integration for shear terms
!                                                               ------------------------------------------------------------
               IORD_IJ_SH = 3
               IORD_K_SH  = 2
               CALL ORDER_TRIA  ( IORD_IJ_SH, SS_I_SH, SS_J_SH, HH_IJ_SH )
               CALL ORDER_GAUSS ( IORD_K_SH , SS_K_SH , HH_K_SH  )

               GAUSS_PT = 0
               IORD_MSG = 'reduced: tria-cross-section x z direction   = '
               DO K=1,IORD_K_SH
                  DO IJ=1,IORD_IJ_SH
                     GAUSS_PT = GAUSS_PT + 1
                     CALL SHP3DP ( IJ, IJ, K, ELGP, SUBR_NAME, IORD_MSG, IORD_IJ_SH, IORD_K_SH,                                    &
                                   SS_I_SH(IJ), SS_J_SH(IJ), SS_K_SH(K), 'N', PSH, DPSHG )
                     CALL JAC3D ( SS_I_SH(IJ), SS_J_SH(IJ), SS_K_SH(K), DPSHG, 'N', JAC, JACI, DETJ(GAUSS_PT) )
                     CALL MATMULT_FFF ( JACI, DPSHG, 3, 3, ELGP, DPSHX )
                     CALL B3D_ISOPARAMETRIC ( DPSHX, GAUSS_PT, IJ, IJ, K, 'direct strains', 'Y', BI )
                     DO L=4,6
                        DO M=1,3*ELGP
                           B(L,M,GAUSS_PT) = BI(L,M)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF

         ENDIF

      ENDIF opt234

! **********************************************************************************************************************************
! Calculate element thermal loads.

      IF (OPT(2) == 'Y') THEN

         DO N=1,NTSUB
            DO L=1,ELGP
               GRID_DT_ARRAY(L,N) = DT(L,N)
            ENDDO
         ENDDO

         DO N=1,NTSUB

            DUM0 = ZERO
            DUM1 = ZERO

            IORD_MSG = ' '
            GAUSS_PT = 0
            IORD_MSG = 'for triangular-cross section x z direction  = '
            DO K=1,IORD_K
               DO IJ=1,IORD_IJ
                  GAUSS_PT = GAUSS_PT + 1
                  DO L=1,6
                     DO M=1,3*ELGP
                        BI(L,M) = B(L,M,GAUSS_PT)
                     ENDDO
                  ENDDO
                  CALL MATMULT_FFF_T ( BI, EALP, 6, 3*ELGP, 1, DUM0 )
                  INTFAC = DETJ(GAUSS_PT)*HH_IJ(IJ)*HH_K(K)
                  IF (DEBUG(191) == 0) THEN                ! Use temperatures at Gauss points for PTE
                     CALL SHP3DP ( IJ, IJ, K, ELGP, SUBR_NAME, IORD_MSG, IORD_IJ, IORD_K, SS_I(IJ), SS_J(IJ), SS_K(K),             &
                                   'N', PSH, DPSHG )
                     CALL MATMULT_FFF ( PSH, GRID_DT_ARRAY, 1, ELGP, NTSUB, TGAUSS )
                     TEMP = TGAUSS(1,N) - TREF1
                  ELSE                                     ! Use avg element temperature for PTE
                     TEMP = TBAR(N)
                  ENDIF
                  DO L=1,3*ELGP
                     DUM1(L) = DUM1(L) + DUM0(L)*TEMP*INTFAC
                  ENDDO
               ENDDO
            ENDDO

            DO L=1,3*ELGP
               PTE(ID(L),N) = DUM1(L)
            ENDDO

         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate SE1, SE2 matrices for stress data recovery. All stresses calculated at center of element

      IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN

        DO STR_PT_NUM=1,7

                                                            ! Isoparametric coordinates of the point
          SELECT CASE (STR_PT_NUM)
            CASE (1)                                       ! Center
              SSI = THIRD; SSJ = THIRD; SSK = ZERO
            CASE (2)                                       ! Node 1
              SSI =  0; SSJ =  0; SSK = -1
            CASE (3)                                       ! Node 2
              SSI = +1; SSJ =  0; SSK = -1
            CASE (4)                                       ! Node 3
              SSI =  0; SSJ = +1; SSK = -1
            CASE (5)                                       ! Node 4
              SSI =  0; SSJ =  0; SSK = +1
            CASE (6)                                       ! Node 5
              SSI = +1; SSJ =  0; SSK = +1
            CASE (7)                                       ! Node 6
              SSI =  0; SSJ = +1; SSK = +1
          END SELECT

         DUM2 = ZERO

          IORD_MSG = 'for triangular-cross section x z direction  = '
          CALL SHP3DP ( 1, 1, 1, ELGP, SUBR_NAME, IORD_MSG, 1, 1, SSI, SSJ, SSK, 'N', PSH, DPSHG )
          CALL JAC3D ( SSI, SSJ, SSK, DPSHG, 'N', JAC, JACI, DUM_DETJ )
          CALL MATMULT_FFF ( JACI, DPSHG, 3, 3, ELGP, DPSHX )
          CALL B3D_ISOPARAMETRIC ( DPSHX, 0, 1, 1, 1, 'all strains', 'N', BI )
          CALL MATMULT_FFF ( ES, BI, 6, 6, 3*ELGP, DUM2 )

          DO I=1,3
            DO J=1,3*ELGP
               SE1(I,ID(J),STR_PT_NUM) = DUM2(I  ,J)
               SE2(I,ID(J),STR_PT_NUM) = DUM2(I+3,J)
            ENDDO
          ENDDO

          DO J=1,NTSUB                                     ! Thermal recovery data (see ELEM_STRE_STRN_ARRAYS): temperature
            IF (DEBUG(191) == 0) THEN                      ! minus TREF at the point (element average, as for PTE, with
               STE3(1,J,STR_PT_NUM) = DOT_PRODUCT(PSH, DT(1:ELGP,J)) - TREF1   ! DEBUG(191) > 0)
            ELSE
               STE3(1,J,STR_PT_NUM) = TBAR(J)
            ENDIF
          ENDDO

          DO I=1,3                                         ! Strain-displ matrix
            DO J=1,3*ELGP
               BE1(I,ID(J),STR_PT_NUM) = BI(I  ,J)
               BE2(I,ID(J),STR_PT_NUM) = BI(I+3,J)
            ENDDO
          ENDDO

        ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE.

      IF (OPT(4) == 'Y') THEN

         DUM3 = ZERO

         IORD_MSG = ' '
         GAUSS_PT = 0
         DO K=1,IORD_K
            DO IJ=1,IORD_IJ
               GAUSS_PT = GAUSS_PT + 1
               DO L=1,6
                  DO M=1,3*ELGP
                     BI(L,M) = B(L,M,GAUSS_PT)
                  ENDDO
               ENDDO
               CALL MATMULT_FFF ( ES, BI, 6, 6, 3*ELGP, DUM4 )
               CALL MATMULT_FFF_T ( BI, DUM4, 6, 3*ELGP, 3*ELGP, DUM5 )
               INTFAC = DETJ(GAUSS_PT)*HH_IJ(IJ)*HH_K(K)
               DO L=1,3*ELGP
                  DO M=1,3*ELGP
                     DUM3(L,M) = DUM3(L,M) + DUM5(L,M)*INTFAC
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         DO I=1,3*ELGP
            DO J=1,3*ELGP
               KE(ID(I),ID(J)) = DUM3(I,J)
            ENDDO
         ENDDO

! Set lower triangular portion of KE equal to upper portion

         DO I=2,6*ELGP
            DO J=1,I-1
               KE(I,J) = KE(J,I)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate linear differential stiffness matrix

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

        CALL ELMDIS
        CALL ELEM_STRE_STRN_ARRAYS ( 1, 'Y' )              ! Stress at the center, in local elem coords like B (not STR_CID coords)

        SIGxx = STRESS(1)
        SIGyy = STRESS(2)
        SIGzz = STRESS(3)
        SIGxy = STRESS(4)
        SIGyz = STRESS(5)
        SIGzx = STRESS(6)

        KWW(1,1) = SIGxx   ;   KWW(1,2) = SIGxy   ;   KWW(1,3) = SIGzx    ! Full initial stress tensor: KED = int(G' S G dV)
        KWW(2,2) = SIGyy   ;   KWW(2,3) = SIGyz                      ! for each displacement component (rotation-only
        KWW(3,3) = SIGzz                                           ! form omega' K omega dropped the strain terms)
        KWW(2,1) = KWW(1,2)
        KWW(3,1) = KWW(1,3)
        KWW(3,2) = KWW(2,3)

        DUM3 = ZERO

                                                           ! KED = int( CBAR^T * KWW * CBAR , dV)
        IORD_MSG = ' '
        GAUSS_PT = 0
        DO K=1,IORD_K
          DO IJ=1,IORD_IJ
            GAUSS_PT = GAUSS_PT + 1
            DO L=1,ELGP
              DPSHX(1,L) = B(1,3*(L-1)+1,GAUSS_PT)         ! Collect shape function derivatives at this Gauss point.
              DPSHX(2,L) = B(2,3*(L-1)+2,GAUSS_PT)
              DPSHX(3,L) = B(3,3*(L-1)+3,GAUSS_PT)
              CBAR(1,3*(L-1)+1) =  ZERO            ; CBAR(1,3*(L-1)+2) = -HALF*DPSHX(3,L) ; CBAR(1,3*(L-1)+3)=  HALF*DPSHX(2,L)
              CBAR(2,3*(L-1)+1) =  HALF*DPSHX(3,L) ; CBAR(2,3*(L-1)+2) =  ZERO            ; CBAR(2,3*(L-1)+3)= -HALF*DPSHX(1,L)
              CBAR(3,3*(L-1)+1) = -HALF*DPSHX(2,L) ; CBAR(3,3*(L-1)+2) =  HALF*DPSHX(1,L) ; CBAR(3,3*(L-1)+3)=  ZERO
            ENDDO
            CALL MATMULT_FFF ( KWW, DPSHX, 3, 3, ELGP, DUM9 )
            CALL MATMULT_FFF_T ( DPSHX, DUM9, 3, ELGP, ELGP, DUM8 )
            DUM5 = ZERO
            DO L=1,ELGP
              DO M=1,ELGP
                DUM5(3*(L-1)+1,3*(M-1)+1) = DUM8(L,M)
                DUM5(3*(L-1)+2,3*(M-1)+2) = DUM8(L,M)
                DUM5(3*(L-1)+3,3*(M-1)+3) = DUM8(L,M)
              ENDDO
            ENDDO
            INTFAC = DETJ(GAUSS_PT)*HH_IJ(IJ)*HH_K(K)
            DO L=1,3*ELGP
              DO M=1,3*ELGP
                DUM3(L,M) = DUM3(L,M) + DUM5(L,M)*INTFAC
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        DO I=1,3*ELGP
          DO J=1,3*ELGP
            KED(ID(I),ID(J)) = DUM3(I,J)
          ENDDO
        ENDDO

        DO I=2,6*ELGP                                      ! Set lower triangular portion of KE equal to upper portion
          DO J=1,I-1
            KED(I,J) = KED(J,I)
          ENDDO
        ENDDO

        IF (DEBUG(178) >= 1) THEN
          WRITE(F06,2001) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
          K = 0
          DO I=1,ELGP
            K = 6*(I-1) + 1
            WRITE(F06,2101) (KED(K  ,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+1,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+2,J),J=1,3*ELGP)
            WRITE(F06,*)
            WRITE(F06,2101) (KED(K+3,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+4,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+5,J),J=1,3*ELGP)
            WRITE(F06,*)
          ENDDO
          WRITE(F06,*)
        ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1925 FORMAT(' *ERROR  1925: ELEMENT ',I8,', TYPE ',A,', HAS ZERO OR NEGATIVE ',A,' = ',1ES9.1)

 1935 FORMAT(' *ERROR  1935: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ELEMENT ',I8,', TYPE ',A,' HAS MATERIAL COUPLING BETWEEN DIRECT AND SHEAR STRESSES. THIS SUBR IS NOT',&
                           ' PROGRAMMED FOR THIS')

 1941 FORMAT(' *ERROR  1941: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ELEMENT ',I8,', TYPE ',A,' HAS AN INVALID NUMBER OF NODES = ',I8)

 2001 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': KED Differential Stiffness Matrix for ', A,            &
             ' element number ', I8)

 2101 FORMAT(1000(1ES14.6))

! **********************************************************************************************************************************

      END SUBROUTINE PENTA


      SUBROUTINE TETRA ( OPT, INT_ELEM_ID, IORD, RED_INT_SHEAR, WRITE_WARN )

! Isoparametric tetrahedron solid element (4 or 10 nodes). Full gaussian integration is the only option (no reduced integration)

! Subroutine calculates:

!  1) ME        = element mass matrix                  , if OPT(1) = 'Y'
!  2) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  3) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  4) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  5) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_TETRA, NTSUB
      USE CONSTANTS_1, ONLY           :  HALF, QUARTER, ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, BE1, BE2, DT, EID, ELGP, NUM_EMG_FATAL_ERRS, ES, KE, KED, ME, PTE, RHO,           &
                                         SE1, SE2, STE3, STRESS, TREF, TYPE

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE ELEMENT_RECOVERY_SUPPORT, ONLY:  ELEM_STRE_STRN_ARRAYS
      USE MASS_DOF_EXPANSION, ONLY    :  EXPAND_MASS_DOFS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TETRA'
      CHARACTER( 1*BYTE), INTENT(IN)  :: RED_INT_SHEAR           ! If 'Y', use Gaussian weighted avg of B matrices for shear terms
      CHARACTER( 1*BYTE), INTENT(IN)  :: OPT(6)                  ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN              ! If 'Y" write warning messages, otherwise do not
      INTEGER(LONG)                   :: GAUSS_PT                ! Gauss point number (used for DEBUG output in subr SHP3DP
      CHARACTER(46*BYTE)              :: IORD_MSG                ! Character name of an integ order (used for debug output)

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID             ! Internal element ID
      INTEGER(LONG), INTENT(IN)       :: IORD                    ! Gaussian integ order for element
      INTEGER(LONG)                   :: I,J,K,L,M,N             ! DO loop indices
      INTEGER(LONG)                   :: II,JJ                   ! Counters
      INTEGER(LONG)                   :: ID(3*ELGP)              ! Array which shows equivalence of DOF's in virgin element with the
     !                                                             6 DOF/grid of the final element stiffness matrix
                                                                 ! Indicator of no output of elem data to BUG file

      INTEGER(LONG)                   :: STR_PT_NUM              ! Stress point number. 1 is center, 2+ are element nodes 1+.

      REAL(DOUBLE)                    :: ALP(6)                  ! First col of ALPVEC

      REAL(DOUBLE)                    :: B(6,3*ELGP,IORD*IORD*IORD)
                                                                 ! Strain-displ matrix for this element for all Gauss points

      REAL(DOUBLE)                    :: BI(6,3*ELGP)            ! Strain-displ matrix for this element for one Gauss point

      REAL(DOUBLE)                    :: CBAR(3,3*ELGP)          ! Derivatives of shape fcns wrt x,y,z used in diff stiff matrix
!                                                                  (contains terms from DPSHX matrices for each grid of the HEXA)

      REAL(DOUBLE)                    :: DETJ(IORD*IORD*IORD)    ! Determinant of JAC for all Gauss points
      REAL(DOUBLE)                    :: DUM_DETJ
      REAL(DOUBLE)                    :: DPSHG(3,ELGP)           ! Output from subr SHP3DT. Derivatives of PSH wrt elem isopar coord
      REAL(DOUBLE)                    :: DPSHX(3,ELGP)           ! Derivatives of PSH wrt elem x, y coords.
      REAL(DOUBLE)                    :: DUM0(3*ELGP)            ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM1(3*ELGP)            ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM2(6,3*ELGP)          ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM3(3*ELGP,3*ELGP)     ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM4(6,3*ELGP)          ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: DUM5(3*ELGP,3*ELGP)     ! Intermediate matrix used in solving for KE elem matrices
      REAL(DOUBLE)                    :: DUM6(3,3*ELGP)          ! Intermediate matrix used in solving for elem matrices
      REAL(DOUBLE)                    :: EALP(6)                 ! Intermed var used in calc PTE therm lds & STEi therm stress coeff
                                                                 ! Array of all DT values at the grids GRID_DT_ARRAY(i,j) = DT(i,j)
      REAL(DOUBLE)                    :: GRID_DT_ARRAY(ELGP,NTSUB)

      REAL(DOUBLE)                    :: HHH_IJK(MAX_ORDER_TETRA)! An output from subr ORDER_TRIA, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: INTFAC                  ! An integ factor (constant multiplier for the Gauss integ)
      REAL(DOUBLE)                    :: JAC(3,3)                ! An output from subr JAC3D, called herein. 3 x 3 Jacobian matrix.
      REAL(DOUBLE)                    :: JACI(3,3)               ! An output from subr JAC3D, called herein. 3 x 3 Jacobian inverse.
      REAL(DOUBLE)                    :: KWW(3,3)                ! Portion of differential stiffness matrix
      REAL(DOUBLE)                    :: DUM8(ELGP,ELGP)         ! G' S G for one Gauss point (differential stiffness)
      REAL(DOUBLE)                    :: DUM9(3,ELGP)            ! S G for one Gauss point (differential stiffness)
      REAL(DOUBLE)                    :: PSH(ELGP)               ! Output from subr SHP3DT. Shape fcn at Gauss pts SSI, SSJ
      REAL(DOUBLE)                    :: SIGxx                   ! Normal stress in the elem x  direction
      REAL(DOUBLE)                    :: SIGyy                   ! Normal stress in the elem y  direction
      REAL(DOUBLE)                    :: SIGzz                   ! Normal stress in the elem z  direction
      REAL(DOUBLE)                    :: SIGxy                   ! Shear  stress in the elem xy direction
      REAL(DOUBLE)                    :: SIGyz                   ! Shear  stress in the elem yz direction
      REAL(DOUBLE)                    :: SIGzx                   ! Shear  stress in the elem zx direction
      REAL(DOUBLE)                    :: SSS_I(MAX_ORDER_TETRA)  ! An output from subr ORDER_TRIA, called herein. Gauss abscissa's.
      REAL(DOUBLE)                    :: SSS_J(MAX_ORDER_TETRA)  ! An output from subr ORDER_TRIA, called herein. Gauss abscissa's.
      REAL(DOUBLE)                    :: SSS_K(MAX_ORDER_TETRA)  ! An output from subr ORDER_TRIA, called herein. Gauss abscissa's.
      REAL(DOUBLE)                    :: TBAR(NTSUB)             ! Average elem temperature for each subcase
      REAL(DOUBLE)                    :: TEMP                    ! Temperature to use in PTE calc
      REAL(DOUBLE)                    :: TGAUSS(1,NTSUB)         ! Temp at a Gauss point for a theral subcase
      REAL(DOUBLE)                    :: TREF1                   ! TREF(1)
      REAL(DOUBLE)                    :: VOLUME                  ! 3D element volume
      REAL(DOUBLE)                    :: SSI,SSJ,SSK             ! Isoparametric coordinates of a point.
      REAL(DOUBLE)                    :: M_1DOF(ELGP,ELGP)      ! Consistent mass matrix with 1 DOF per node.



! **********************************************************************************************************************************

! Calculate volume by Gaussian integration

      CALL ORDER_TETRA ( IORD, SSS_I, SSS_J, SSS_K, HHH_IJK )
      IORD_MSG = ' '
      IORD_MSG = 'for 3-D solid strains,           input IORD = '
      VOLUME = ZERO
      DO I=1,IORD
         CALL SHP3DT ( I, ELGP, SUBR_NAME, IORD_MSG, IORD, SSS_I(I), SSS_J(I), SSS_K(I), 'Y', PSH, DPSHG )
         CALL JAC3D ( SSS_I(I), SSS_J(I), SSS_K(I), DPSHG, 'Y', JAC, JACI, DETJ(I) )
         VOLUME = VOLUME + HHH_IJK(I)*DETJ(I)
      ENDDO

! If VOLUME <= 0, write error and return. det(J) at every Gauss point has already been checked in subr JAC3D relative to the
! element's own dimensions, so no absolute tolerance (which would depend on the element size and model units) is used here.

      IF (VOLUME <= ZERO) THEN
         WRITE(ERR,1925) EID, TYPE, 'VOLUME', VOLUME
         WRITE(F06,1925) EID, TYPE, 'VOLUME', VOLUME
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

! Calculate ID array

      JJ = 0
      DO I=1,ELGP
         II = 6*(I - 1)
         DO J=1,3
            II = II + 1
            JJ = JJ + 1
            ID(JJ) = II
         ENDDO
      ENDDO

! ALP is first col of ALPVEC
      DO I=1,6
        ALP(I) = ALPVEC(I,1)
      ENDDO

      TREF1 = TREF(1)

! EALP is needed to calculate both PTE and STE2

      CALL MATMULT_FFF ( ES, ALP, 6, 6, 1, EALP )

! Calc TBAR (used for PTE, STEi)

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         DO J=1,NTSUB
            TBAR(J) = ZERO
            DO K=1,ELGP
               TBAR(J) = TBAR(J) + DT(K,J)
            ENDDO
            TBAR(J) = TBAR(J)/ELGP - TREF1
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Generate the mass matrix for this element.

      IF (OPT(1) == 'Y') THEN

         ! Consistent mass matrix
         ! ME = ∫ N' ρ N det(J) dv

         M_1DOF(:,:) = ZERO
         GAUSS_PT = 0
         DO I=1,IORD
            GAUSS_PT = GAUSS_PT + 1

            CALL SHP3DT ( I, ELGP, SUBR_NAME, IORD_MSG, IORD, SSS_I(I), SSS_J(I), SSS_K(I), 'N', PSH, DPSHG )
            INTFAC = DETJ(I)*HHH_IJK(I)

            DO L=1,ELGP
               DO M=1,ELGP
                  M_1DOF(L,M) = M_1DOF(L,M) + PSH(L) * PSH(M) * RHO(1) * INTFAC
               ENDDO
            ENDDO
         ENDDO


         CALL EXPAND_MASS_DOFS( M_1DOF )

      ENDIF

! **********************************************************************************************************************************
! Generate B matrices. The B matrix terms for Gauss point k are B(i,j,k)

      B = ZERO

      IORD_MSG = 'for 3-D solid strains,           input IORD = '
      DO I=1,IORD
         CALL SHP3DT ( I, ELGP, SUBR_NAME, IORD_MSG, IORD, SSS_I(I), SSS_J(I), SSS_K(I), 'N', PSH, DPSHG )
         CALL JAC3D ( SSS_I(I), SSS_J(I), SSS_K(I), DPSHG, 'N', JAC, JACI, DETJ(I) )
         CALL MATMULT_FFF ( JACI, DPSHG, 3, 3, ELGP, DPSHX )
         CALL B3D_ISOPARAMETRIC ( DPSHX, I, I, I, I, 'direct strains', 'Y', BI )
         DO L=1,6
            DO M=1,3*ELGP
               B(L,M,I) = BI(L,M)
            ENDDO
         ENDDO
      ENDDO

! **********************************************************************************************************************************
! Calculate element thermal loads.

      IF (OPT(2) == 'Y') THEN

         DO N=1,NTSUB
            DO L=1,ELGP
               GRID_DT_ARRAY(L,N) = DT(L,N)
            ENDDO
         ENDDO

         DO N=1,NTSUB

            DUM0 = ZERO
            DUM1 = ZERO

            IORD_MSG = 'for 3-D solid strains,           input IORD = '
            GAUSS_PT = 0
            DO I=1,IORD
               GAUSS_PT = GAUSS_PT + 1
               DO L=1,6
                  DO M=1,3*ELGP
                     BI(L,M) = B(L,M,GAUSS_PT)
                  ENDDO
               ENDDO
               CALL MATMULT_FFF_T ( BI, EALP, 6, 3*ELGP, 1, DUM0 )
               INTFAC = DETJ(GAUSS_PT)*HHH_IJK(I)
               IF (DEBUG(191) == 0) THEN                ! Use temperatures at Gauss points for PTE
                  CALL SHP3DT ( I, ELGP, SUBR_NAME, IORD_MSG, IORD, SSS_I(I), SSS_J(I), SSS_K(I), 'N', PSH, DPSHG )
                  CALL MATMULT_FFF ( PSH, GRID_DT_ARRAY, 1, ELGP, NTSUB, TGAUSS )
                  TEMP = TGAUSS(1,N) - TREF1
               ELSE                                     ! Use avg element temperature for PTE
                  TEMP = TBAR(N)
               ENDIF
               DO L=1,3*ELGP
                  DUM1(L) = DUM1(L) + DUM0(L)*TEMP*INTFAC
               ENDDO
            ENDDO

            DO L=1,3*ELGP
               PTE(ID(L),N) = DUM1(L)
            ENDDO

         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate BEi, SEi matrices for strain, stress data recovery.

      IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN

        DO STR_PT_NUM=1,9
                                                           ! Isoparametric coordinates of the point
          SELECT CASE (STR_PT_NUM)
            CASE (1)                                       ! Center
              SSI = QUARTER; SSJ = QUARTER; SSK = QUARTER
            CASE (2)                                       ! Node 1
              SSI =  0; SSJ =  0; SSK =  0
            CASE (3)                                       ! Node 2
              SSI = +1; SSJ =  0; SSK =  0
            CASE (4)                                       ! Node 3
              SSI =  0; SSJ = +1; SSK =  0
            CASE (5)                                       ! Node 4
              SSI =  0; SSJ =  0; SSK = +1
          END SELECT

         DUM2 = ZERO

          GAUSS_PT = 1                                     ! Calc SE1,2,3
          IORD_MSG = 'for 3-D solid strains,                      = '
          CALL SHP3DT ( 1, ELGP, SUBR_NAME, IORD_MSG, 1, SSI, SSJ, SSK, 'N', PSH, DPSHG )
          CALL JAC3D ( SSI, SSJ, SSK, DPSHG, 'N', JAC, JACI, DUM_DETJ )
          CALL MATMULT_FFF ( JACI, DPSHG, 3, 3, ELGP, DPSHX )
          CALL B3D_ISOPARAMETRIC ( DPSHX, 0, 1, 1, 1, 'all strains', 'N', BI )
          CALL MATMULT_FFF ( ES, BI, 6, 6, 3*ELGP, DUM2 )

          DO I=1,3
            DO J=1,3*ELGP
               SE1(I,ID(J),STR_PT_NUM) = DUM2(I  ,J)
               SE2(I,ID(J),STR_PT_NUM) = DUM2(I+3,J)
            ENDDO
          ENDDO

          DO J=1,NTSUB                                     ! Thermal recovery data (see ELEM_STRE_STRN_ARRAYS): temperature
            IF ((DEBUG(191) == 0) .AND. (ELGP /= 4)) THEN  ! minus TREF at the point. TETRA4 has a constant strain, so it
               STE3(1,J,STR_PT_NUM) = DOT_PRODUCT(PSH, DT(1:ELGP,J)) - TREF1   ! always uses the element average, which is
            ELSE                                           ! also what DEBUG(191) > 0 selects for every solid (as for PTE)
               STE3(1,J,STR_PT_NUM) = TBAR(J)
            ENDIF
          ENDDO

          DO I=1,3                                         ! Strain-displ matrix
            DO J=1,3*ELGP
               BE1(I,ID(J),STR_PT_NUM) = BI(I  ,J)
               BE2(I,ID(J),STR_PT_NUM) = BI(I+3,J)
            ENDDO
          ENDDO

        ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE.

      IF (OPT(4) == 'Y') THEN

         DUM3 = ZERO

         IORD_MSG = ' '
         DO I=1,IORD


            BI = B(:,:,I)


            CALL MATMULT_FFF ( ES, BI, 6, 6, 3*ELGP, DUM4 )

            CALL MATMULT_FFF_T ( BI, DUM4, 6, 3*ELGP, 3*ELGP, DUM5 )

            INTFAC = DETJ(I)*HHH_IJK(I)

            DO L=1,3*ELGP
               DO M=1,3*ELGP
                  DUM3(L,M) = DUM3(L,M) + DUM5(L,M)*INTFAC
               ENDDO
            ENDDO

         ENDDO

         DO I=1,3*ELGP
            DO J=1,3*ELGP
               KE(ID(I),ID(J)) = DUM3(I,J)
            ENDDO
         ENDDO

! Set lower triangular portion of KE equal to upper portion

         DO I=2,6*ELGP
            DO J=1,I-1
               KE(I,J) = KE(J,I)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate linear differential stiffness matrix
      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

        CALL ELMDIS
        CALL ELEM_STRE_STRN_ARRAYS ( 1, 'Y' )              ! Stress at the center, in local elem coords like B (not STR_CID coords)

        SIGxx = STRESS(1)
        SIGyy = STRESS(2)
        SIGzz = STRESS(3)
        SIGxy = STRESS(4)
        SIGyz = STRESS(5)
        SIGzx = STRESS(6)

        KWW(1,1) = SIGxx   ;   KWW(1,2) = SIGxy   ;   KWW(1,3) = SIGzx    ! Full initial stress tensor: KED = int(G' S G dV)
        KWW(2,2) = SIGyy   ;   KWW(2,3) = SIGyz                      ! for each displacement component (rotation-only
        KWW(3,3) = SIGzz                                           ! form omega' K omega dropped the strain terms)
        KWW(2,1) = KWW(1,2)
        KWW(3,1) = KWW(1,3)
        KWW(3,2) = KWW(2,3)

        DO I=1,3*ELGP
          DO J=1,3*ELGP
            DUM3(I,J) = ZERO
          ENDDO
        ENDDO

                                                           ! KED = int( CBAR^T * KWW * CBAR , dV)
        IORD_MSG = ' '
        GAUSS_PT = 0
        DO I=1,IORD
          GAUSS_PT = GAUSS_PT + 1
          DO L=1,ELGP
            DPSHX(1,L) = B(1,3*(L-1)+1,GAUSS_PT)           ! Collect shape function derivatives at this Gauss point.
            DPSHX(2,L) = B(2,3*(L-1)+2,GAUSS_PT)
            DPSHX(3,L) = B(3,3*(L-1)+3,GAUSS_PT)
            CBAR(1,3*(L-1)+1) =  ZERO             ;  CBAR(1,3*(L-1)+2) = -HALF*DPSHX(3,L)  ;  CBAR(1,3*(L-1)+3)=  HALF*DPSHX(2,L)
            CBAR(2,3*(L-1)+1) =  HALF*DPSHX(3,L)  ;  CBAR(2,3*(L-1)+2) =  ZERO             ;  CBAR(2,3*(L-1)+3)= -HALF*DPSHX(1,L)
            CBAR(3,3*(L-1)+1) = -HALF*DPSHX(2,L)  ;  CBAR(3,3*(L-1)+2) =  HALF*DPSHX(1,L)  ;  CBAR(3,3*(L-1)+3)=  ZERO
          ENDDO
          CALL MATMULT_FFF ( KWW, DPSHX, 3, 3, ELGP, DUM9 )
          CALL MATMULT_FFF_T ( DPSHX, DUM9, 3, ELGP, ELGP, DUM8 )
          DUM5 = ZERO
          DO L=1,ELGP
            DO M=1,ELGP
              DUM5(3*(L-1)+1,3*(M-1)+1) = DUM8(L,M)
              DUM5(3*(L-1)+2,3*(M-1)+2) = DUM8(L,M)
              DUM5(3*(L-1)+3,3*(M-1)+3) = DUM8(L,M)
            ENDDO
          ENDDO
          INTFAC = DETJ(I)*HHH_IJK(I)
          DO L=1,3*ELGP
            DO M=1,3*ELGP
              DUM3(L,M) = DUM3(L,M) + DUM5(L,M)*INTFAC
            ENDDO
          ENDDO
        ENDDO

        DO I=1,3*ELGP
          DO J=1,3*ELGP
            KED(ID(I),ID(J)) = DUM3(I,J)
          ENDDO
        ENDDO

        DO I=2,6*ELGP                                      ! Set lower triangular portion of KE equal to upper portion
          DO J=1,I-1
            KED(I,J) = KED(J,I)
          ENDDO
        ENDDO

        IF (DEBUG(178) >= 1) THEN
          WRITE(F06,2001) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
          K = 0
          DO I=1,ELGP
            K = 6*(I-1) + 1
            WRITE(F06,2101) (KED(K  ,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+1,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+2,J),J=1,3*ELGP)
            WRITE(F06,*)
            WRITE(F06,2101) (KED(K+3,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+4,J),J=1,3*ELGP)
            WRITE(F06,2101) (KED(K+5,J),J=1,3*ELGP)
            WRITE(F06,*)
          ENDDO
          WRITE(F06,*)
        ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1925 FORMAT(' *ERROR  1925: ELEMENT ',I8,', TYPE ',A,', HAS ZERO OR NEGATIVE ',A,' = ',1ES9.1)

 2001 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': KED Differential Stiffness Matrix for ', A,            &
             ' element number ', I8)

 2101 FORMAT(1000(1ES14.6))


! **********************************************************************************************************************************

      END SUBROUTINE TETRA

   END MODULE SOLID_ELEMENTS
