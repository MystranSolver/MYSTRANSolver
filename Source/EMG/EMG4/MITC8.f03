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

   MODULE MITC8_MOD

   USE MITC_KERNELS, ONLY :  MITC_ADD_TO_B, MITC_CONTRAVARIANT_BASIS, MITC_COVARIANT_BASIS, MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION, MITC_DETJ, MITC_ELASTICITY, MITC_INITIALIZE, MITC_SHAPE_FUNCTIONS, MITC_TRANSFORM_B

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MITC8

   CONTAINS
      SUBROUTINE MITC8 ( OPT, INT_ELEM_ID )

! Calculates, or calls subr's to calculate, quadrilateral element matrices:

!  1) ME        = element mass matrix                  , if OPT(1) = 'Y'
!  2) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  3) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  4) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  5) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  6) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, MAX_STRESS_POINTS
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  NUM_EMG_FATAL_ERRS, PCOMP_PROPS, ELGP, ES, KE, EM, ET, BE1, BE2, BE3, PHI_SQ, FCONV,      &
                                         EPROP, SHELL_STR_ANGLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
      USE PARAMS, ONLY                :  TSTM_DEF
      USE MITC_STUF, ONLY             :  GP_RS

      USE QUADRATURE, ONLY            :  ORDER_GAUSS
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      USE VECTOR_GEOMETRY, ONLY       :  CROSS
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MITC8'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID
      INTEGER(LONG), PARAMETER        :: IORD_IJ = 3       ! Integration order for stiffness matrix
      INTEGER(LONG), PARAMETER        :: IORD_K = 2        ! Integration order for stiffness matrix in thickness direction
      INTEGER(LONG), PARAMETER        :: IORD_STRESS_Q8 = 2! Gauss integration order for stress/strain recovery matrices
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG)                   :: STR_PT_NUM        ! Stress recovery point number

      REAL(DOUBLE)                    :: HH_IJ(MAX_ORDER_GAUSS) ! Gauss weights for integration in in-layer directions
      REAL(DOUBLE)                    :: SS_IJ(MAX_ORDER_GAUSS) ! Gauss abscissa's for integration in in-layer directions
      REAL(DOUBLE)                    :: HH_K(MAX_ORDER_GAUSS)  ! Gauss weights for integration in thickness direction
      REAL(DOUBLE)                    :: SS_K(MAX_ORDER_GAUSS)  ! Gauss abscissa's for integration in thickness direction
      REAL(DOUBLE)                    :: R, S, T                ! Isoparametric coordinates of a point
      REAL(DOUBLE)                    :: BI(6,6*ELGP)      ! Strain-displ matrix for this element for one Gauss point
      REAL(DOUBLE)                    :: BI1(6,6*ELGP)     ! Strain-displ matrix for this element for one Gauss point bottom
      REAL(DOUBLE)                    :: BI2(6,6*ELGP)     ! Strain-displ matrix for this element for one Gauss point top
      REAL(DOUBLE)                    :: DUM1(6,6*ELGP)    ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM2(6*ELGP,6*ELGP)    ! Intermediate matrix
      REAL(DOUBLE)                    :: INTFAC            ! An integration factor (constant multiplier for the Gauss integration)
      REAL(DOUBLE)                    :: DETJ              ! Jacobian determinant
      REAL(DOUBLE)                    :: E(6,6)            ! Elasticity matrix in the material coordinate system.
      REAL(DOUBLE)                    :: EE(6,6)           ! Elasticity matrix in the cartesian local coordinate system.
      REAL(DOUBLE)                    :: LOCAL_BASIS(3,3)  ! Cartesian local basis
      REAL(DOUBLE)                    :: ELEMENT_BASIS(3,3)! Element coordinate system basis
      REAL(DOUBLE)                    :: XL(3)
      REAL(DOUBLE)                    :: ZL(3)
      REAL(DOUBLE)                    :: XE(3)
      REAL(DOUBLE)                    :: CROSS_XLE(3)

! **********************************************************************************************************************************

! COORDINATE SYSTEMS
! ==================
!
! Cartesian local
!  e^_1, e^_2, e^_3 in Bathe but they are oriented differently there.
!  e^_3 is the midsurface normal which is the same as the director vector for MITC8 because it doesn't support SNORM.
!  Used for strain in the strain-displacement matrix and the material elasticity matrix is transformed to this to integrate KE.
!  Defined the same way as the zero THETA material coordinate system in Siemens SimCenter. This definition is used because it
!  has uniform orientation on distorted flat elements and is only non-uniform as needed to accomodate out-of-plane curvature.
!  The uniformity allows stress to be interpolated and extrapolated to different locations conveniently.
!  Orthogonal
!
! Element local (Nastran definition)
!  x_l, y_l, z_l
!  Used for element stress, strain, and force outputs
!  Defined the same way as MSC (element coordinate system) and SimCenter (local coordinate system).
!  Defined by x_l being the bisection of the R, S isoparametric basis vectors rotated about the normal by -45 degrees.
!  Orthogonal
!
! XEL element (internal use)
!  Used for the grid point DOFs of the strain-displacement and the element stiffness matrices.
!  Used for extrapolating stress and strain from Gauss points to corners.
!  Grid point coordinates stored in XEL are in a coordinate system which is flat, with the normal being the cross product
!  of vectors from grid points 1-3 and 2-4. The x axis is an arbitrary direction in this plane. The flat coordinate system
!  is used for grid point coordinates for extrapolating stress because the polynomial curve fit code to extrapolate stress/strain
!  from Gauss points to corners is only 2D.
!  Orthogonal
!
! Material
!  Used for material elasticity read from the input file.
!  Currently, this is the same as the cartesian local coordinate system. To allow non-isotropic materials, it
!  should find the angle between the two systems's x axes at each integration point and rotate the material
!  elasticity matrix about that when building the stiffness matrix KE.
!  Orthogonal
!
! Isoparametric (natural)
!  R, S, T in code. r_1, r_2, r_3 in Bathe.
!  Each coordinate has range [-1,1].
!  T is parallel to the director vector.
!  Not orthogonal
!
! Covariant
!  g_r, g_s, g_t. g_1, g_2, g_3 in Bathe.
!  Parallel to the isoparametric coordinates but scaled by the element size. Eg. |g_t| = half thickness.
!  Not orthogonal
!
! Contravariant
!  g^r, g^s, g^t. g^1, g^2, g^3 in Bathe.
!  Contravariant to the covariant.
!  Not orthogonal
!

! **********************************************************************************************************************************

! Initialize
      PHI_SQ  = ONE                                        ! Not used for this element
      CALL MITC_INITIALIZE ()


      IF (PCOMP_PROPS == 'Y') THEN
        WRITE(ERR,*) ' *ERROR: Code not written for composite material with QUAD8'
        WRITE(F06,*) ' *ERROR: Code not written for composite material with QUAD8'
        NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
        FATAL_ERR = FATAL_ERR + 1
        CALL OUTA_HERE ( 'Y' )
      ENDIF



! **********************************************************************************************************************************
! Generate the mass matrix for this element.

      IF (OPT(1) == 'Y') THEN
        !Not implememented yet but we can't make it a fatal error because this gets called even when it doesn't need it.

      ENDIF



! **********************************************************************************************************************************
! Calculate element thermal loads.

      IF (OPT(2) == 'Y') THEN

        WRITE(ERR,*) ' *ERROR: Code not written for QUAD8 thermal loads'
        WRITE(F06,*) ' *ERROR: Code not written for QUAD8 thermal loads'
        NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
        FATAL_ERR = FATAL_ERR + 1
        CALL OUTA_HERE ( 'Y' )

      ENDIF

! **********************************************************************************************************************************
! BE1 matrix (3 x 48) for membrane strain/stress/force data recovery.
! BE2 matrix (3 x 48) for bending strain/stress/force data recovery.
! BE3 matrix (2 x 48) for transverse shear strain/stress/force data recovery.
! All calculated at Gauss points and not center.
! The displacements are in basic coordinates and the strains are in element coordinates.


! There's a possible bug where the numbering of the element grid points (eg. 1-2-3-4 vs 2-3-4-1) affects the stress/strain/elforce,
! when the element is distorted. This includes von Mises which should be invariant to rotation.
!
! - It doesn't occur if we skip the extrapolation from Gauss points to corners so Gauss point strains are probably OK.
! - It makes no differences if the cartesian local coordiante system is changed to be G1G2 or the bisected diagonals projected
!   onto the surface everywhere or like Siemens material coordinates with an intermediate reference plane.
! - It works OK if the cartesian local coordinate system is defined using the same vector for each element, eg. (1,0,0), instead
!   of G1G2 with either the Siemens or direct projection. However, this won't generalize to elements in any orientation and might
!   just be hiding the problem.
! - The problem is probably the extrapolation from Gauss points to grid points. Maybe it should be done in covariant coordinates
!   the way strains are interpolated by MITC. Somehow.

      IF (OPT(3) == 'Y') THEN

         STR_PT_NUM = 1

         CALL ORDER_GAUSS ( IORD_STRESS_Q8, SS_IJ, HH_IJ )

         DO I=1,IORD_STRESS_Q8
            DO J=1,IORD_STRESS_Q8

               STR_PT_NUM = STR_PT_NUM + 1

               R = SS_IJ(I)
               S = SS_IJ(J)

               CALL MITC8_B( R, S, -ONE, .TRUE., .TRUE., BI1)
               CALL MITC8_B( R, S, +ONE, .TRUE., .TRUE., BI2)

                                                  ! Membrane strain is the average of the strains at the two t points.
               BE1(1,:,STR_PT_NUM) = (BI2(1,:) + BI1(1,:)) / TWO           ! xx
               BE1(2,:,STR_PT_NUM) = (BI2(2,:) + BI1(2,:)) / TWO           ! yy
               BE1(3,:,STR_PT_NUM) = (BI2(4,:) + BI1(4,:)) / TWO           ! xy

                                                  ! Curvature is (strain_top - strain_bottom) / thickness
                                                  ! To allow grid point thicknesses, this should be the thickness
                                                  ! interpolated at the Gauss point.
               BE2(1,:,STR_PT_NUM) = (BI2(1,:) - BI1(1,:)) / EPROP(1)      ! xx
               BE2(2,:,STR_PT_NUM) = (BI2(2,:) - BI1(2,:)) / EPROP(1)      ! yy
               BE2(3,:,STR_PT_NUM) = (BI2(4,:) - BI1(4,:)) / EPROP(1)      ! xy

                                                  ! Transverse shear strain. Note reversed order of rows.
               BE3(1,:,STR_PT_NUM) = (BI2(6,:) + BI1(6,:)) / TWO           ! zx
               BE3(2,:,STR_PT_NUM) = (BI2(5,:) + BI1(5,:)) / TWO           ! yz

            ENDDO
         ENDDO

                                                           ! Find angle of the element coordinate system's x axis from
                                                           ! the cartesian local coordinate system's x axis at each
                                                           ! corner.
                                                           ! This will be used to transform stress and strain to the
                                                           ! element coordinate system after extrapolating to corners.
         DO STR_PT_NUM=2,5

            R = GP_RS(1, STR_PT_NUM - 1)
            S = GP_RS(2, STR_PT_NUM - 1)

            LOCAL_BASIS = MITC8_CARTESIAN_LOCAL_BASIS( R, S )
            XL = LOCAL_BASIS(:,1)                          ! X axis of cartesian local basis
            ZL = LOCAL_BASIS(:,3)                          ! Normal
            ELEMENT_BASIS = MITC8_ELEMENT_CS_BASIS( R, S )
            XE = ELEMENT_BASIS(:,1)                        ! X axis of element coordinate system

            CALL CROSS( XL, XE, CROSS_XLE )
            SHELL_STR_ANGLE( STR_PT_NUM ) = ATAN2(DOT_PRODUCT( ZL, CROSS_XLE ), DOT_PRODUCT( XL, XE ))

         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE.

      IF(OPT(4) == 'Y') THEN

! Based on
! MITC4 paper "A continuum mechanics based four-node shell element for general nonlinear analysis"
!   by Dvorkin and Bathe
! MITC8 paper "A FORMULATION OF GENERAL SHELL ELEMENTS-THE USE OF MIXED INTERPOLATION OF TENSORIAL COMPONENTS"
!   by Dvorkin and Bathe, 1986


         ! K = int( [B]^T [EE] [B] dV )
         !    dV = |det(J)|dr ds dt
         ! K = int( [B]^T [EE] [B] |det(J)| dr ds dt )
         !
         ! [EE] is material elasticity matrix in cartesian local coordinates
         ! stress = [EE] * strain
         ! strain = [B] * displacement
         ! K is in the basic coordinate system

         E = MITC_ELASTICITY()

         KE(1:6*ELGP,1:6*ELGP) = ZERO

         CALL ORDER_GAUSS ( IORD_IJ, SS_IJ, HH_IJ )
         CALL ORDER_GAUSS ( IORD_K, SS_K, HH_K )

         DO I=1,IORD_IJ
            DO J=1,IORD_IJ
               DO K=1,IORD_K
                  R = SS_IJ(I)
                  S = SS_IJ(J)
                  T = SS_K(K)
                  CALL MITC8_B( R, S, T, .TRUE., .TRUE., BI)

                  ! For non-isotropic materials, this should be rotated from the material coordinate system to the cartesian local
                  ! coordinate system here. The rotation angle may be different at each Gauss point.
                  EE(:,:) = E(:,:)

                  CALL MATMULT_FFF ( EE, BI, 6, 6, 6*ELGP, DUM1 )
                  CALL MATMULT_FFF_T ( BI, DUM1, 6, 6*ELGP, 6*ELGP, DUM2 )
                  DETJ = MITC_DETJ ( R, S, T )
                  INTFAC = DETJ*HH_IJ(I)*HH_IJ(J)*HH_K(K)
                  KE(1:6*ELGP,1:6*ELGP) = KE(1:6*ELGP,1:6*ELGP) + DUM2(:,:)*INTFAC
               ENDDO
            ENDDO
         ENDDO


      ENDIF


! **********************************************************************************************************************************
! Determine element pressure loads

      IF (OPT(5) == 'Y') THEN

        WRITE(ERR,*) ' *ERROR: Code not written for QUAD8 pressure loads'
        WRITE(F06,*) ' *ERROR: Code not written for QUAD8 pressure loads'
        NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
        FATAL_ERR = FATAL_ERR + 1
        CALL OUTA_HERE ( 'Y' )

      ENDIF

! **********************************************************************************************************************************
! Calculate linear differential stiffness matrix

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

        WRITE(ERR,*) ' *ERROR: Code not written for QUAD8 differential stiffness matrix'
        WRITE(F06,*) ' *ERROR: Code not written for QUAD8 differential stiffness matrix'
        NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
        FATAL_ERR = FATAL_ERR + 1
        CALL OUTA_HERE ( 'Y' )

      ENDIF






      RETURN

! **********************************************************************************************************************************


! **********************************************************************************************************************************

      END SUBROUTINE MITC8

      SUBROUTINE MITC8_B ( R, S, T, INLAYER, SHEAR, B )

! Calculates the strain-displacement matrix in the cartesian local coordinate system
! for MITC8 shell at one point in isoparametric coordinates.
! Based on
! MITC8 paper "A FORMULATION OF GENERAL SHELL ELEMENTS-THE USE OF MIXED INTERPOLATION OF TENSORIAL COMPONENTS"
!   by Dvorkin and Bathe, 1986


      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP
      USE CONSTANTS_1, ONLY           :  ZERO, QUARTER, ONE, TWO, THREE

      USE FULL_MATRIX_ALGEBRA, ONLY   :  OUTER_PRODUCT

      IMPLICIT NONE

      LOGICAL, INTENT(IN)             :: INLAYER           ! TRUE for in-layer rows (1-4)
      LOGICAL, INTENT(IN)             :: SHEAR             ! TRUE for transverse shear rows (5-6)

      INTEGER(LONG)                   :: I, J              ! Basis vector indices. 1=R, 2=S
      INTEGER(LONG)                   :: K, L              ! DO loop indices
      INTEGER(LONG)                   :: ROW               ! A row of B. 1-6.
      INTEGER(LONG)                   :: COL               ! A column (element DOF) of B. 1-48.
      INTEGER(LONG)                   :: POINT             ! Sampling point number
      INTEGER(LONG)                   :: POINT_A           ! Corner sampling point number adjacent to midside one.
      INTEGER(LONG)                   :: POINT_B           ! Corner sampling point number adjacent to midside one.

      REAL(DOUBLE) , INTENT(IN)       :: R, S, T           ! Isoparametric coordinates
      REAL(DOUBLE) , INTENT(OUT)      :: B(6, 6*ELGP)      ! Strain-displacement matrix
      REAL(DOUBLE)                    :: E(6, 6*ELGP)      ! Strain-displacement matrix directly interpolated
      REAL(DOUBLE)                    :: A
      REAL(DOUBLE)                    :: POINT_R(8)        ! Sampling point isoparametric R coordinates
      REAL(DOUBLE)                    :: POINT_S(8)        ! Sampling point isoparametric S coordinates
      REAL(DOUBLE)                    :: POINT_COORDS_1(8) ! Shear sampling point isoparametric R or S coordinates
      REAL(DOUBLE)                    :: POINT_COORDS_2(8) ! Shear sampling point isoparametric S or R coordinates
      REAL(DOUBLE)                    :: PSH(ELGP)         ! Shape functions (interpolation functions)
      REAL(DOUBLE)                    :: DPSHG(2,ELGP)     ! Derivatives of shape functions with respect to R and S.
      REAL(DOUBLE)                    :: H_IS(8)           ! Interpolation functions indexed by in-layer sampling point number.
      REAL(DOUBLE)                    :: H_IT(5)           ! Interpolation functions indexed by shear sampling point number.
      REAL(DOUBLE)                    :: G(3,3)            ! Array of 3 covariant basis vectors in basic coordinates
      REAL(DOUBLE)                    :: G_CONTRA(3,3)     ! Array of 3 contravariant basis vectors in basic coordinates
      REAL(DOUBLE)                    :: G_J_NORMALIZED(3)
      REAL(DOUBLE)                    :: GG(3,3)           ! Outer product of two contravariant basis vectors
      REAL(DOUBLE)                    :: B_1(6,6*ELGP,4)   ! Part of strain-displacement matrix for each of 4 sampling points.
      REAL(DOUBLE)                    :: B_2(6,6*ELGP,1)   ! Part of strain-displacement matrix for a sampling point.
      REAL(DOUBLE)                    :: E_AVERAGE(6,6*ELGP)
      REAL(DOUBLE)                    :: EJJ(6,6*ELGP)     ! Part of strain-displacement matrix with only row J used.
      REAL(DOUBLE)                    :: EIT(6,6*ELGP)     ! Part of strain-displacement matrix with only row RT or ST.
      REAL(DOUBLE)                    :: EIT_RA(6,6*ELGP)  ! Part of strain-displacement matrix with only row RT or ST at RA.
      REAL(DOUBLE)                    :: EIT_RB(6,6*ELGP)  ! Part of strain-displacement matrix with only row RT or ST at RB.
      REAL(DOUBLE)                    :: G_JJ
      REAL(DOUBLE)                    :: G_RS
      REAL(DOUBLE)                    :: SCALAR
      REAL(DOUBLE)                    :: DUM1(3)
      REAL(DOUBLE)                    :: COORD_1           ! R or S coordinate of a sampling point.
      REAL(DOUBLE)                    :: COORD_2           ! S or R coordinate of a sampling point.
      REAL(DOUBLE)                    :: TRANSFORM(3,3)    ! Transformation matrix

      INTRINSIC                       :: DSQRT

! **********************************************************************************************************************************
! Initialize empty matrix

      B(:,:) = ZERO


! **********************************************************************************************************************************
! Add in-layer strain-displacement terms


      IF ( INLAYER ) THEN

                                                           ! 8 in-layer sampling points
        !         ^ s
        !         |
        !  4      7      3          o sampling point
        !   +-----+-----+           + node
        !   | o   o   o |
        !   |  2  5  1  |
        ! 8 + o6     8o + 6 --> r
        !   |  3  7  4  |
        !   | o   o   o |
        !   +-----+-----+
        !  1      5      2
        !

        A = ONE / DSQRT(THREE)
        POINT_R = (/ A, -A, -A,  A, ZERO, -A,     ZERO, A    /)
        POINT_S = (/ A,  A, -A, -A, A,     ZERO, -A,    ZERO /)

                                                           ! Interpolation functions for the 8 sampling points.
        CALL MITC_SHAPE_FUNCTIONS(R/A, S/A, PSH, DPSHG)


                                                            ! Convert interpolation functions from native node numbering
                                                            ! to Bathe interpolation point/node numbering.
        H_IS(1) = PSH(3)
        H_IS(2) = PSH(4)
        H_IS(3) = PSH(1)
        H_IS(4) = PSH(2)
        H_IS(5) = PSH(7)
        H_IS(6) = PSH(8)
        H_IS(7) = PSH(5)
        H_IS(8) = PSH(6)


        DO POINT=1,4

          B_1(:,:,POINT) = ZERO

          ! The contribution to ε in the global (basic) basis is
          ! = hIS(i) * ( ε~_rr g^r g^r|_i  +  ε~_ss g^s g^s|_i  +  ε~_rs [g^r g^s|_i  + g^r g^s|_i^T] )

          CALL MITC_COVARIANT_BASIS( POINT_R(POINT), POINT_S(POINT), T, G )
          CALL MITC_CONTRAVARIANT_BASIS( G, G_CONTRA )
          CALL MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION( POINT_R(POINT), POINT_S(POINT), T, 1, 4, E )

                                                           ! ε~_rr g^r g^r|_i
          CALL OUTER_PRODUCT( G_CONTRA(:,1), G_CONTRA(:,1), 3, 3, GG )
          DO COL=1,6*ELGP
            CALL MITC_ADD_TO_B( B_1, POINT, COL, E(1, COL), GG )
          ENDDO

                                                           ! ε~_ss g^s g^s|_i
          CALL OUTER_PRODUCT( G_CONTRA(:,2), G_CONTRA(:,2), 3, 3, GG )
          DO COL=1,6*ELGP
            CALL MITC_ADD_TO_B( B_1, POINT, COL, E(2, COL), GG )
          ENDDO

                                                           ! ε~_rs [g^r g^s|_i  + g^r g^s|_i^T]
          CALL OUTER_PRODUCT( G_CONTRA(:,1), G_CONTRA(:,2), 3, 3, GG )
          GG = GG + TRANSPOSE(GG)
          DO COL=1,6*ELGP
            CALL MITC_ADD_TO_B( B_1, POINT, COL, E(4, COL), GG )
          ENDDO

          B = B + B_1(:,:,POINT) * H_IS(POINT)

        ENDDO



        DO POINT=5,8

          SELECT CASE (POINT)
            CASE (5); I=1; J=2; POINT_A=1; POINT_B=2
            CASE (6); I=2; J=1; POINT_A=2; POINT_B=3
            CASE (7); I=1; J=2; POINT_A=3; POINT_B=4
            CASE (8); I=2; J=1; POINT_A=4; POINT_B=1
          END SELECT

                                                           ! ε_ss g^s g^s |_SamplingPoint if sampling point is 5 or 7
                                                           ! ε_rr g^r g^r |_SamplingPoint if sampling point is 6 or 8
          CALL MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION( POINT_R(POINT), POINT_S(POINT), T, J, J, EJJ )
          CALL MITC_COVARIANT_BASIS(POINT_R(POINT), POINT_S(POINT), T, G)

          !              _
          ! Convert g to g
          !
          ! This only has an effect for non-rectangular elements where r and s aren't perpendicular.
          ! I assume g_rs means the component of g_r in the s direction.
          ! G_JJ = |G_J| where J is R or S.
          G_JJ = DSQRT( G(1,J)*G(1,J) + G(2,J)*G(2,J) + G(3,J)*G(3,J) )

          ! I assume g_rs in both formulas is a mistake and it should be g_sr in the other one.
          ! This is because g_rs in both wouldn't be symmetric between the two sets of sampling points.
          ! I decided which should be g_rs and which g_sr from slightly better convergence in a skewed element test
          ! and in a skewed curved thin strip test.
          ! This choice was also slightly better than making no adjustments to g at all.
          G_J_NORMALIZED(:) = G(:,J) / G_JJ
          G_RS = G(1,I)*G_J_NORMALIZED(1) + G(2,I)*G_J_NORMALIZED(2) + G(3,I)*G_J_NORMALIZED(3)
          G(:,I) = G(:,I) - G(:,J) * G_RS / G_JJ

          CALL MITC_CONTRAVARIANT_BASIS( G, G_CONTRA )

          B_2(:,:,1) = ZERO

          CALL OUTER_PRODUCT( G_CONTRA(:,J), G_CONTRA(:,J), 3, 3, GG )
          DO COL=1,6*ELGP
            CALL MITC_ADD_TO_B( B_2, 1, COL, EJJ(J, COL), GG )
          ENDDO

                                                           ! [1/2 (ε|_A + ε|_B)]
          E_AVERAGE = (B_1(:,:,POINT_A) + B_1(:,:,POINT_B)) / TWO

                      ! {g_r · [1/2 (ε|_A + ε|_B)] · g_r} g^r g^r |_SamplingPoint if sampling point is 5 or 7
                      ! {g_s · [1/2 (ε|_A + ε|_B)] · g_s} g^s g^s |_SamplingPoint if sampling point is 6 or 8
          CALL OUTER_PRODUCT( G_CONTRA(:,I), G_CONTRA(:,I), 3, 3, GG )
          DO COL=1,6*ELGP
                                                           ! g_i · [...]
            DUM1(1) = G(1,I) * E_AVERAGE(1,COL) + G(2,I) * E_AVERAGE(4,COL) + G(3,I) * E_AVERAGE(6,COL)
            DUM1(2) = G(1,I) * E_AVERAGE(4,COL) + G(2,I) * E_AVERAGE(2,COL) + G(3,I) * E_AVERAGE(5,COL)
            DUM1(3) = G(1,I) * E_AVERAGE(6,COL) + G(2,I) * E_AVERAGE(5,COL) + G(3,I) * E_AVERAGE(3,COL)
                                                           ! [...] · g_i
            SCALAR = DUM1(1) * G(1,I) + DUM1(2) * G(2,I) + DUM1(3) * G(3,I)
            CALL MITC_ADD_TO_B( B_2, 1, COL, SCALAR, GG )
          ENDDO

                      ! {g_r · [1/2 (ε|_A + ε|_B)] · g_s} (g^r g^s + [g^r g^s]^T |_SamplingPoint
          CALL OUTER_PRODUCT( G_CONTRA(:,1), G_CONTRA(:,2), 3, 3, GG )
          GG = GG + TRANSPOSE(GG)
          DO COL=1,6*ELGP
                                                           ! g_r · [...]
            DUM1(1) = G(1,1) * E_AVERAGE(1,COL) + G(2,1) * E_AVERAGE(4,COL) + G(3,1) * E_AVERAGE(6,COL)
            DUM1(2) = G(1,1) * E_AVERAGE(4,COL) + G(2,1) * E_AVERAGE(2,COL) + G(3,1) * E_AVERAGE(5,COL)
            DUM1(3) = G(1,1) * E_AVERAGE(6,COL) + G(2,1) * E_AVERAGE(5,COL) + G(3,1) * E_AVERAGE(3,COL)
                                                           ! [...] · g_s
            SCALAR = DUM1(1) * G(1,2) + DUM1(2) * G(2,2) + DUM1(3) * G(3,2)
            CALL MITC_ADD_TO_B( B_2, 1, COL, SCALAR, GG )
          ENDDO

          B = B + B_2(:,:,1) * H_IS(POINT)

        ENDDO


      ENDIF




! **********************************************************************************************************************************
! Add transverse shear strain-displacement terms

      IF ( SHEAR ) THEN

        DO I=1,2

          SELECT CASE (I)
            CASE (1) ; COORD_1 = R; COORD_2 = S; ROW = 6   ! ε_rt
            CASE (2) ; COORD_1 = S; COORD_2 = R; ROW = 5   ! ε_st
          END SELECT


          ! 7 sampling points for ε_rt (I=1)
          !         ^ s
          !         |
          !  4      7      3          o sampling point
          !   +-o---+---o-+           + node
          !   | 2       1 |
          !   |           |
          ! 8 + o   o   o + 6 --> r
          !   | RB  5  RA |
          !   | 3       4 |
          !   +-o---+---o-+
          !  1      5      2
          !

          ! 7 sampling points for ε_st (I=2)
          ! Points 2 and 4 are interchanged compared to Bathe because this makes it symmetric with ε_rt to reuse code.
          !         ^ s
          !         |
          !  7      6      5          o sampling point
          !   +-----+-----+           + node
          !   o4  SAo    1o
          !   |           |
          ! 8 +   5 o     + 4 --> r
          !   |           |
          !   o3  SBo    2o
          !   +-----+-----+
          !  1      2      3
          !

          ! Points 6 and 7 are called RA and RB respectively in Bathe.
          A = ONE / DSQRT(THREE)
          POINT_COORDS_1 = (/ A,   -A,   -A,    A,   ZERO, A,   -A,    ZERO /)
          POINT_COORDS_2 = (/ ONE,  ONE, -ONE, -ONE, ZERO, ZERO, ZERO, ZERO /)
          SELECT CASE (I)
            CASE (1) ; POINT_R = POINT_COORDS_1; POINT_S = POINT_COORDS_2
            CASE (2) ; POINT_R = POINT_COORDS_2; POINT_S = POINT_COORDS_1
          END SELECT

          ! Interpolation functions for the 5 numbered sampling points.
          H_IT(5) = (ONE - (COORD_1 / A) ** TWO) * (ONE - COORD_2 ** TWO)
          H_IT(1) = QUARTER * (ONE + COORD_1 / A) * (ONE + COORD_2) - QUARTER * H_IT(5)
          H_IT(2) = QUARTER * (ONE - COORD_1 / A) * (ONE + COORD_2) - QUARTER * H_IT(5)
          H_IT(3) = QUARTER * (ONE - COORD_1 / A) * (ONE - COORD_2) - QUARTER * H_IT(5)
          H_IT(4) = QUARTER * (ONE + COORD_1 / A) * (ONE - COORD_2) - QUARTER * H_IT(5)


          ! The contribution to ε in the global basis is
          !   ε~_rt (g^r g^t + g^t g^r) + ε~_st (g^s g^t + g^t g^s)
          ! = ε~_rt (g^r g^t + [g^r g^t]^T) + ε~_st (g^s g^t + [g^s g^t]^T)
          ! ε~_rt g^r g^t = sum over sampling points 1-4 of (h ε~_rt g^r g^t) evaluated at the sampling point
          ! plus a special term for sampling point 5.

          DO POINT=1,5

            ! I guess that the transverse shear sampling points are at the nodal surface (t=0, T=0). Same as MITC4.
            CALL MITC_COVARIANT_BASIS( POINT_R(POINT), POINT_S(POINT), ZERO, G )
            CALL MITC_CONTRAVARIANT_BASIS( G, G_CONTRA )

            IF(POINT == 5) THEN
              CALL MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION( POINT_R(6), POINT_S(6), ZERO, ROW, ROW, EIT_RA )
              CALL MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION( POINT_R(7), POINT_S(7), ZERO, ROW, ROW, EIT_RB )
              EIT = (EIT_RA + EIT_RB) / TWO
            ELSE
              CALL MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION( POINT_R(POINT), POINT_S(POINT), ZERO, ROW, ROW, EIT )
            ENDIF

            B_2(:,:,1) = ZERO

            CALL OUTER_PRODUCT( G_CONTRA(:,I), G_CONTRA(:,3), 3, 3, GG )
            ! Evaluate both the g^r g^t term and the [g^r g^t]^T term together so their sum is symmetric since only
            ! the unique elements are stored in the B matrix.
            GG = GG + TRANSPOSE(GG)
            DO COL=1,6*ELGP
              CALL MITC_ADD_TO_B( B_2, 1, COL, EIT(ROW, COL), GG )
            ENDDO

            B = B + B_2(:,:,1) * H_IT(POINT)

          ENDDO




        ENDDO


      ENDIF


! **********************************************************************************************************************************
! Transform strain from the global cartesian basis (basic) to the cartesian local basis.

      TRANSFORM = TRANSPOSE(MITC8_CARTESIAN_LOCAL_BASIS(R, S))
      CALL MITC_TRANSFORM_B( TRANSFORM, B)

      ! Double shear terms because it's now treated as vectors instead of tensors.
      B(4:6,:) = B(4:6,:) * 2

      RETURN

! **********************************************************************************************************************************


! **********************************************************************************************************************************

      END SUBROUTINE MITC8_B

      FUNCTION MITC8_CARTESIAN_LOCAL_BASIS ( R, S )

! Finds the basis vectors of the cartesian local coordinate system expressed in the basic coordinate system.
! This is defined the same way as the material coordinate system in Simcenter Nastran with THETA=0.
! This definition is chosen because it is uniform over the surface of the element so stress and strain outputs
! can be interpolated/extrapolated in this coordinate system then transformed to the local element coordinate
! system at the corner grid points and center for output.
!
! First index of the result (row) is a vector component in basic coordinates (x,y,z)
! Second index of the result (column) is basis vector (x_l, y_l, normal)

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP, XEL, TYPE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO

      USE VECTOR_GEOMETRY, ONLY       :  CROSS

      IMPLICIT NONE

      INTEGER(LONG)                   :: I                 ! DO loop indices

      REAL(DOUBLE)                    :: MITC8_CARTESIAN_LOCAL_BASIS(3,3)
      REAL(DOUBLE) , INTENT(IN)       :: R
      REAL(DOUBLE) , INTENT(IN)       :: S
      REAL(DOUBLE)                    :: PSH(ELGP)
      REAL(DOUBLE)                    :: DPSHG(2,ELGP)     ! Derivatives of shape functions with respect to R and S.
      REAL(DOUBLE)                    :: E_XI(3)
      REAL(DOUBLE)                    :: E_ETA(3)
      REAL(DOUBLE)                    :: Z_REF(3)
      REAL(DOUBLE)                    :: R_G1G2(3)
      REAL(DOUBLE)                    :: X(3)
      REAL(DOUBLE)                    :: Y(3)
      REAL(DOUBLE)                    :: Z(3)


! **********************************************************************************************************************************


      CALL MITC_SHAPE_FUNCTIONS(R, S, PSH, DPSHG)

                                                           ! Unit normal to the reference plane
      CALL CROSS(XEL(3,:) - XEL(1,:), XEL(4,:) - XEL(2,:), Z_REF)
      Z_REF = Z_REF / DSQRT(DOT_PRODUCT(Z_REF, Z_REF))

                                                           ! Project R_G1G2 onto the reference plane
      R_G1G2 = XEL(2,:) - XEL(1,:)
      R_G1G2 = R_G1G2 - Z_REF * DOT_PRODUCT(R_G1G2, Z_REF) / DOT_PRODUCT(Z_REF, Z_REF)

                                                           ! Unit normal to shell surface (Z)
      E_XI(:)=ZERO
      E_ETA(:)=ZERO
      DO I=1,ELGP
        E_XI(:) = E_XI(:) + XEL(I,:) * DPSHG(1,I)
        E_ETA(:) = E_ETA(:) + XEL(I,:) * DPSHG(2,I)
      ENDDO
      CALL CROSS(E_XI, E_ETA, Z)
      Z = Z / DSQRT(DOT_PRODUCT(Z, Z))

                                                           ! Y tangent to the surface
      CALL CROSS(Z, R_G1G2, Y)
      Y = Y / DSQRT(DOT_PRODUCT(Y, Y))

                                                           ! Rotate the projected R_G1G2 about Y to be tangent to the surface
      CALL CROSS(Y, Z, X)

      MITC8_CARTESIAN_LOCAL_BASIS(:,1) = X
      MITC8_CARTESIAN_LOCAL_BASIS(:,2) = Y
      MITC8_CARTESIAN_LOCAL_BASIS(:,3) = Z


      RETURN


! **********************************************************************************************************************************

      END FUNCTION MITC8_CARTESIAN_LOCAL_BASIS

      FUNCTION MITC8_ELEMENT_CS_BASIS ( R, S )

! Finds the basis vectors of the local element coordinate system expressed in the basic coordinate system.
! This is defined the same way as in MSC Nastran and is used for element stress, strain and force output.
! First index of the result is a vector component in basic coordinates (x,y,z)
! Second index of the result is basis vector (x_l, y_l, normal)

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP, XEL, TYPE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO

      USE VECTOR_GEOMETRY, ONLY       :  CROSS

      IMPLICIT NONE

      INTEGER(LONG)                   :: I                 ! DO loop indices

      REAL(DOUBLE)                    :: MITC8_ELEMENT_CS_BASIS(3,3)
      REAL(DOUBLE) , INTENT(IN)       :: R
      REAL(DOUBLE) , INTENT(IN)       :: S
      REAL(DOUBLE)                    :: PSH(ELGP)
      REAL(DOUBLE)                    :: DPSHG(2,ELGP)     ! Derivatives of shape functions with respect to R and S.
      REAL(DOUBLE)                    :: E_XI(3)
      REAL(DOUBLE)                    :: E_ETA(3)
      REAL(DOUBLE)                    :: A(3)
      REAL(DOUBLE)                    :: B(3)
      REAL(DOUBLE)                    :: X_L_ACB(3)
      REAL(DOUBLE)                    :: Y_L_ACB(3)
      REAL(DOUBLE)                    :: T(3,3)

      INTRINSIC                       :: DSQRT


! **********************************************************************************************************************************


      CALL MITC_SHAPE_FUNCTIONS(R, S, PSH, DPSHG)

                                                           ! e_ξ(r, s) = d/dR X = sum over nodes[ dN/dR X ]
                                                           ! e_η(r, s) = d/dS X = sum over nodes[ dN/dS X ]
      E_XI(:)=ZERO
      E_ETA(:)=ZERO
      DO I=1,ELGP
        E_XI(:) = E_XI(:) + XEL(I,:) * DPSHG(1,I)
        E_ETA(:) = E_ETA(:) + XEL(I,:) * DPSHG(2,I)
      ENDDO

                                                           !Normalize e_ξ and e_η
      E_XI = E_XI / DSQRT(DOT_PRODUCT(E_XI, E_XI))
      E_ETA = E_ETA / DSQRT(DOT_PRODUCT(E_ETA, E_ETA))

                                                           ! A = bisection of e_ξ and e_η
      A = E_XI + E_ETA
      A = A / DSQRT(DOT_PRODUCT(A, A))

                                                           ! B = common normal of e_ξ and e_η
      CALL CROSS(E_XI, E_ETA, B)
      B = B / DSQRT(DOT_PRODUCT(B, B))

                                                           ! x_l and y_l in the A C B coordinate system.
      X_L_ACB = (/ ONE/DSQRT(TWO), -ONE/DSQRT(TWO), ZERO /)
      Y_L_ACB = (/ ONE/DSQRT(TWO),  ONE/DSQRT(TWO), ZERO /)

                                                           ! Rotation matrix from A C B to basic coordinates
      T(:,1) = A
      CALL CROSS(B, A, T(:,2))
      T(:,3) = B

                                                           ! Transform x_l and y_l to the basic coordinate system
      MITC8_ELEMENT_CS_BASIS(:,1) = MATMUL(T, X_L_ACB)
      MITC8_ELEMENT_CS_BASIS(:,2) = MATMUL(T, Y_L_ACB)
      MITC8_ELEMENT_CS_BASIS(:,3) = B


      RETURN


! **********************************************************************************************************************************

      END FUNCTION MITC8_ELEMENT_CS_BASIS

   END MODULE MITC8_MOD
