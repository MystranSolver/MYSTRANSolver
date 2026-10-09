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

      SUBROUTINE HEXA ( OPT, INT_ELEM_ID,IORD, RED_INT_SHEAR, WRITE_WARN )

! Isoparametric hexahedron solid element (8 or 20 nodes and with full or reduced gaussian integration)

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
      USE CONSTANTS_1, ONLY           :  QUARTER, HALF, ZERO, ONE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  AGRID, ALPVEC, BE1, BE2, DT, EID, ELGP, NUM_EMG_FATAL_ERRS, ES, KE, KED, ME,              &
                                         NUM_EMG_FATAL_ERRS, PLOAD4_3D_DATA, PPE, PRESS, PTE, RHO, SE1, SE2, STE1, STRESS, TREF,   &
                                         TYPE, XEL
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE RESULT_COORDINATES, ONLY    :  ELMDIS
      USE ELEMENT_RECOVERY_SUPPORT, ONLY:  ELEM_STRE_STRN_ARRAYS
      USE MASS_DOF_EXPANSION, ONLY    :  EXPAND_MASS_DOFS


      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'HEXA'
      CHARACTER( 1*BYTE), INTENT(IN)  :: RED_INT_SHEAR          ! If 'Y', use Gaussian weighted avg of B matrices for shear terms
      CHARACTER( 1*BYTE), INTENT(IN)  :: OPT(6)                 ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN             ! If 'Y" write warning messages, otherwise do not

      CHARACTER( 1*BYTE)              :: GENERATE_PRESS_LOAD    ! If 'Y' this element will have a pressure load on at least 1 face
!                                                                 for at least 1 subcase

      CHARACTER(46*BYTE)              :: IORD_MSG               ! Character name of an integration order (used for debug output)

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID             ! Internal element ID
      INTEGER(LONG), INTENT(IN)       :: IORD                   ! Gaussian integration order for element
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
                                                                ! Indicator of no output of elem data to BUG file

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
      REAL(DOUBLE)                    :: DUM_DETJ
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

      REAL(DOUBLE)                    :: EALP(6)                ! Variable used in calc PTE therm loads & STEi therm stress coeffs
      REAL(DOUBLE)                    :: EPS1                   ! A small number to compare to real zero
      REAL(DOUBLE)                    :: FACE_AREA              ! Area of a face of the HEXA where a PLOAD4 pressure acts

                                                                ! Array of all DT values at the grids GRID_DT_ARRAY(i,j) = DT(i,j)
      REAL(DOUBLE)                    :: GRID_DT_ARRAY(ELGP,NTSUB)

      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS)   ! An output from subr ORDER_GAUSS, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: HHH_SH(MAX_ORDER_GAUSS)! An output from subr ORDER_GAUSS, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: INTFAC                 ! An integration factor (constant multiplier for the Gauss integ)
      REAL(DOUBLE)                    :: JAC(3,3)               ! An output from subr JAC3D, called herein. 3 x 3 Jacobian matrix.
      REAL(DOUBLE)                    :: JACI(3,3)              ! An output from subr JAC3D, called herein. 3 x 3 Jacobian inverse.
      REAL(DOUBLE)                    :: KWW(3,3)               ! Portion of differential stiffness matrix
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

      EPS1 = EPSIL(1)

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

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y')) THEN
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

! If VOLUME <= 0, write error and return

      IF (VOLUME < EPS1) THEN
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

         IF (RED_INT_SHEAR == 'Y') THEN

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
                     DO L=1,6
                        DO M=1,3*ELGP
                           BI(L,M) = B(L,M,GAUSS_PT)
                        ENDDO
                     ENDDO
                     DUM0 = MATMUL(TRANSPOSE(BI), EALP)
                     INTFAC = DETJ(GAUSS_PT)*HHH(I)*HHH(J)*HHH(K)
                     IF (DEBUG(191) == 0) THEN             ! Use temperatures at Gauss points for PTE
                        CALL SHP3DH ( I, J, K, ELGP, SUBR_NAME, IORD_MSG, IORD, SSS(I), SSS(J), SSS(K), 'N', PSH, DPSHG )
                        TGAUSS(1,:) = MATMUL(PSH, GRID_DT_ARRAY)
                        TEMP = TGAUSS(1,N) - TREF1
                     ELSE                                  ! Use avg element temperature for PTE
                        TEMP = TBAR(N)
                     ENDIF
                     DO L=1,3*ELGP
                        DUM1(L) = DUM1(L) + DUM0(L)*TEMP*INTFAC
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO

            DO L=1,3*ELGP
               PTE(ID(L),N) = DUM1(L)
            ENDDO

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
          CALL JAC3D ( SSI, SSJ, SSK, DPSHG, 'N', JAC, JACI, DUM_DETJ )
          DPSHX = MATMUL(JACI, DPSHG)
          CALL B3D_ISOPARAMETRIC ( DPSHX, 0, 1, 1, 1, 'all strains', 'N', BI )
          DUM2 = MATMUL(ES,BI)

          DO I=1,3                                         ! Stress-displ matrices
            DO J=1,3*ELGP
               SE1(I,ID(J),STR_PT_NUM) = DUM2(I  ,J)
               SE2(I,ID(J),STR_PT_NUM) = DUM2(I+3,J)
            ENDDO
          ENDDO

          DO J=1,NTSUB                                     ! STE thermal stress terms
            STE1(1,J,STR_PT_NUM) = EALP(1)*TBAR(J)
            STE1(2,J,STR_PT_NUM) = EALP(2)*TBAR(J)
            STE1(3,J,STR_PT_NUM) = EALP(3)*TBAR(J)
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

         RED3 = ((RED_INT_SHEAR == 'Y') .AND. (IORD == 3))
         IORD_MSG = ' '
         GAUSS_PT = 0
         DO K=1,IORD
            DO J=1,IORD
               DO I=1,IORD
                  GAUSS_PT = GAUSS_PT + 1
                  DO L=1,6
                     DO M=1,3*ELGP
                        BI(L,M) = B(L,M,GAUSS_PT)
                     ENDDO
                  ENDDO
                  IF (RED3) BI(4:6,:) = ZERO                      ! shear terms are integrated separately with 2x2x2 points                  
                  DUM4 = MATMUL(ES,BI)
                  DUM5 = MATMUL(TRANSPOSE(BI),DUM4)
                  INTFAC = DETJ(GAUSS_PT)*HHH(I)*HHH(J)*HHH(K)
                  DO L=1,3*ELGP
                     DO M=1,3*ELGP
                        DUM3(L,M) = DUM3(L,M) + DUM5(L,M)*INTFAC
                     ENDDO
                  ENDDO
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
        CALL ELEM_STRE_STRN_ARRAYS ( 1 )

        SIGxx = STRESS(1)
        SIGyy = STRESS(2)
        SIGzz = STRESS(3)
        SIGxy = STRESS(4)
        SIGyz = STRESS(5)
        SIGzx = STRESS(6)


        IF (.TRUE.) THEN

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
      USE PARAMS, ONLY                :  EPSIL
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, BE1, BE2, DT, EID, ELGP, NUM_EMG_FATAL_ERRS, ES, KE, KED, ME, PTE, RHO,           &
                                         SE1, SE2, STE1, STRESS, TREF, TYPE

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
      REAL(DOUBLE)                    :: EPS1                    ! A small number to compare to real zero

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
      EPS1 = EPSIL(1)

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

      IF (VOLUME < EPS1) THEN                              ! If VOLUME <= 0, write error and return
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

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y')) THEN
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

          DO J=1,NTSUB                                     ! STE thermal stress terms
            STE1(1,J,STR_PT_NUM) = EALP(1)*TBAR(J)
            STE1(2,J,STR_PT_NUM) = EALP(2)*TBAR(J)
            STE1(3,J,STR_PT_NUM) = EALP(3)*TBAR(J)
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
        CALL ELEM_STRE_STRN_ARRAYS ( 1 )

        SIGxx = STRESS(1)
        SIGyy = STRESS(2)
        SIGzz = STRESS(3)
        SIGxy = STRESS(4)
        SIGyz = STRESS(5)
        SIGzx = STRESS(6)

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
            CALL MATMULT_FFF ( KWW, CBAR, 3, 3, 3*ELGP, DUM6 )
            CALL MATMULT_FFF_T ( CBAR, DUM6, 3, 3*ELGP, 3*ELGP, DUM5 )
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
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  ALPVEC, BE1, BE2, DT, EID, ELGP, NUM_EMG_FATAL_ERRS, ES, KE, KED, ME, PTE, RHO,           &
                                         SE1, SE2, STE1, STRESS, TREF, TYPE

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
      REAL(DOUBLE)                    :: EPS1                    ! A small number to compare to real zero
                                                                 ! Array of all DT values at the grids GRID_DT_ARRAY(i,j) = DT(i,j)
      REAL(DOUBLE)                    :: GRID_DT_ARRAY(ELGP,NTSUB)

      REAL(DOUBLE)                    :: HHH_IJK(MAX_ORDER_TETRA)! An output from subr ORDER_TRIA, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: INTFAC                  ! An integ factor (constant multiplier for the Gauss integ)
      REAL(DOUBLE)                    :: JAC(3,3)                ! An output from subr JAC3D, called herein. 3 x 3 Jacobian matrix.
      REAL(DOUBLE)                    :: JACI(3,3)               ! An output from subr JAC3D, called herein. 3 x 3 Jacobian inverse.
      REAL(DOUBLE)                    :: KWW(3,3)                ! Portion of differential stiffness matrix
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
      EPS1 = EPSIL(1)

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

! If VOLUME <= 0, write error and return

      IF (VOLUME < EPS1) THEN
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

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y')) THEN
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

          DO J=1,NTSUB                                     ! STE thermal stress terms
            STE1(1,J,STR_PT_NUM) = EALP(1)*TBAR(J)
            STE1(2,J,STR_PT_NUM) = EALP(2)*TBAR(J)
            STE1(3,J,STR_PT_NUM) = EALP(3)*TBAR(J)
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
        CALL ELEM_STRE_STRN_ARRAYS ( 1 )

        SIGxx = STRESS(1)
        SIGyy = STRESS(2)
        SIGzz = STRESS(3)
        SIGxy = STRESS(4)
        SIGyz = STRESS(5)
        SIGzx = STRESS(6)

        KWW(1,1) = SIGyy + SIGzz   ;   KWW(1,2) = -SIGxy   ;   KWW(1,3) = -SIGzx
        KWW(2,2) = SIGxx + SIGzz   ;   KWW(2,3) = -SIGyz
        KWW(3,3) = SIGxx + SIGyy
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
          CALL MATMULT_FFF ( KWW, CBAR, 3, 3, 3*ELGP, DUM6 )
          CALL MATMULT_FFF_T ( CBAR, DUM6, 3, 3*ELGP, 3*ELGP, DUM5 )
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
