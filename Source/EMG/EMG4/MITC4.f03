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

   MODULE MITC4_MOD

   USE MITC_KERNELS, ONLY :  MITC4_CARTESIAN_LOCAL_BASIS, MITC_CONTRAVARIANT_BASIS, MITC_COVARIANT_BASIS, MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION, MITC_DETJ, MITC_ELASTICITY, MITC_INITIALIZE, MITC_SHAPE_FUNCTIONS, MITC_TRANSFORM_B, MITC_TRANSFORM_CONTRAVARIANT_TO_LOCAL

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MITC4

   CONTAINS
      SUBROUTINE MITC4 ( OPT, INT_ELEM_ID )

! Calculates, or calls subr's to calculate, quadrilateral element matrices:

!  1) ME        = element mass matrix                  , if OPT(1) = 'Y'
!  2) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  3) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  4) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  5) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  6) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, MAX_STRESS_POINTS, NTSUB, NSUB
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  NUM_EMG_FATAL_ERRS, PCOMP_PROPS, ELGP, ES, KE, EM, EB, ET, BE1, BE2, BE3, PHI_SQ,         &
                                         FCONV, EPROP, PTE, ALPVEC, TREF, DT, PPE, PRESS, MASS_PER_UNIT_AREA,                      &
                                         NUM_PLIES, PCOMP_LAM, PLY_NUM, TPLY, STRESS, KED,                                         &
                                         SHELL_A, SHELL_B, SHELL_D, SHELL_T, SHELL_AALP, SHELL_BALP, SHELL_DALP, SHELL_TALP
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, ONE, TWO, FOUR, QUARTER
      USE MITC_STUF, ONLY             :  THERM_KAPPA, MSTRPT
      USE RESULT_COORDINATES, ONLY     :  ELMDIS, ELMDIS_PLY
      USE ELEMENT_RECOVERY_SUPPORT, ONLY:  ELEM_STRE_STRN_ARRAYS

      USE QUADRATURE, ONLY            :  ORDER_GAUSS
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE VECTOR_GEOMETRY, ONLY       :  CROSS, PLANE_COORD_TRANS_21
      USE MATERIAL_TRANSFORMATIONS, ONLY:  MATL_TRANSFORM_MATRIX
      USE MASS_DOF_EXPANSION, ONLY    :  EXPAND_MASS_DOFS

      USE COMPOSITE_SHELL_PREPARATION, ONLY:  GET_ELEM_NUM_PLIES, SHELL_ABD_MATRICES
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MITC4'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID
      INTEGER(LONG), PARAMETER        :: IORD_IJ = 2       ! Integration order for stiffness matrix
      INTEGER(LONG), PARAMETER        :: IORD_K = 2        ! Integration order for stiffness matrix in thickness direction
      INTEGER(LONG), PARAMETER        :: IORD_STRESS_Q4 = 2! Gauss integration order for stress/strain recovery matrices
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG)                   :: STR_PT_NUM        ! Stress recovery point number
      INTEGER(LONG)                   :: GP                ! Element grid point number

      REAL(DOUBLE)                    :: HH_IJ(MAX_ORDER_GAUSS) ! Gauss weights for integration in in-layer directions
      REAL(DOUBLE)                    :: SS_IJ(MAX_ORDER_GAUSS) ! Gauss abscissa's for integration in in-layer directions
      REAL(DOUBLE)                    :: HH_K(MAX_ORDER_GAUSS)  ! Gauss weights for integration in thickness direction
      REAL(DOUBLE)                    :: SS_K(MAX_ORDER_GAUSS)  ! Gauss abscissa's for integration in thickness direction
      REAL(DOUBLE)                    :: R, S, T                ! Isoparametric coordinates of a point
      REAL(DOUBLE)                    :: BI(6,6*ELGP)      ! Strain-displ matrix for this element for one Gauss point
      REAL(DOUBLE)                    :: DUM1(6,6*ELGP)    ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM2(6*ELGP,6*ELGP)    ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM3(6*ELGP,6)    ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM4(3)           ! Intermediate matrix
      REAL(DOUBLE)                    :: INTFAC            ! An integration factor (constant multiplier for the Gauss integration)
      REAL(DOUBLE)                    :: DETJ              ! Jacobian determinant
      REAL(DOUBLE)                    :: E2(6,6)           ! Membrane and shear elasticity matrix in the element coordinate system.
      REAL(DOUBLE)                    :: E3(6,6)           ! Membrane and shear elasticity matrix in the cartesian local coordinate system.
      REAL(DOUBLE)                    :: EM2(6,6)          ! Membrane elasticity matrix in the element coordinate system.
      REAL(DOUBLE)                    :: EM3(6,6)          ! Membrane elasticity matrix in the cartesian local coordinate system.
      REAL(DOUBLE)                    :: EB2(6,6)          ! Bending elasticity matrix in the element coordinate system.
      REAL(DOUBLE)                    :: EB3(6,6)          ! Bending elasticity matrix in the cartesian local coordinate system.
      REAL(DOUBLE)                    :: CLB(3,3)          ! Cartesian local basis basis vectors
      REAL(DOUBLE)                    :: MATL_AXES_ROTATE
      REAL(DOUBLE)                    :: TRANSFORM(3,3)
      REAL(DOUBLE)                    :: DUM66(6,6)        ! Intermediate matrix in calculating outputs
      REAL(DOUBLE)                    :: T66(6,6)          ! 6x6 transformation matrix for elasticity
      REAL(DOUBLE)                    :: CTE(6)            ! Coefficient of thermal expansion vector
      REAL(DOUBLE)                    :: THERMAL_STRAIN(6) ! Thermal strain vector
      REAL(DOUBLE)                    :: TBAR              ! Average elem temperature
      REAL(DOUBLE)                    :: UNIT_PTE(6*ELGP)  ! Thermal load vector for unit temperature change.
      REAL(DOUBLE)                    :: UNIT_PPE(6*ELGP)  ! Pressure load vector for unit pressure.
      REAL(DOUBLE)                    :: PSH(ELGP)
      REAL(DOUBLE)                    :: DPSHG(2,ELGP)     ! Derivatives of shape functions with respect to R and S.
      REAL(DOUBLE)                    :: G(3,3)
      REAL(DOUBLE)                    :: NORMAL(3)

      REAL(DOUBLE)                    :: BMI(6,6*ELGP)     ! Strain-displ matrix for membrane for one Gauss point
      REAL(DOUBLE)                    :: BBI(6,6*ELGP)     ! Strain-displ matrix for bending for one Gauss point
      REAL(DOUBLE)                    :: BSI(6,6*ELGP)     ! Strain-displ matrix for shear for one Gauss point
      REAL(DOUBLE)                    :: BMI3(3,6*ELGP)    ! Mid-surface membrane strain-disp operator (MITC4_BMBS)
      REAL(DOUBLE)                    :: BBI3(3,6*ELGP)    ! Mid-surface curvature strain-disp operator (MITC4_BMBS)
      REAL(DOUBLE)                    :: BSI2(2,6*ELGP)    ! Mid-surface trans shear strain-disp operator (MITC4_BMBS)
      REAL(DOUBLE)                    :: DUM1_3(3,6*ELGP)  ! Intermediate matrix for 3-row (membrane/bending) matmults
      REAL(DOUBLE)                    :: DUM1_2(2,6*ELGP)  ! Intermediate matrix for 2-row (transverse shear) matmults
      REAL(DOUBLE)                    :: M_1DOF(ELGP,ELGP) ! Consistent mass matrix with 1 DOF per node.
      REAL(DOUBLE)                    :: DENSITY
      REAL(DOUBLE)                    :: FORCEx(IORD_STRESS_Q4*IORD_STRESS_Q4) ! Engineering force in the elem x direction at Gauss points
      REAL(DOUBLE)                    :: FORCEy(IORD_STRESS_Q4*IORD_STRESS_Q4) ! Engineering force in the elem x direction at Gauss points
      REAL(DOUBLE)                    :: FORCExy(IORD_STRESS_Q4*IORD_STRESS_Q4)! Engineering force in the elem xy direction at Gauss points
      REAL(DOUBLE)                    :: KS(ELGP,ELGP)     ! KED matrix for one DOF
      REAL(DOUBLE)                    :: DUM12(ELGP,2)     ! Intermediate matrix used in solving for KED matrices
      INTEGER(LONG)                   :: GAUSS_PT          ! Gauss point number (used for output in subr SHP2DQ
                                                           ! An output from subr ORDER, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS)
      INTEGER(LONG)                   :: JPLY              ! PLY_NUM in a DO loop
      REAL(DOUBLE)                    :: SSS(MAX_ORDER_GAUSS)
      REAL(DOUBLE)                    :: DUM11(2,2)        ! Intermediate matrix used in solving for KED matrices
      INTEGER(LONG)                   :: KI, KJ            ! For converting grid point number to element DOF number
      REAL(DOUBLE)                    :: DUM13(ELGP,ELGP)  ! Intermediate matrix used in solving for KED matrices
      REAL(DOUBLE)                    :: DPSHX(2,4)        ! Derivatives of PSH wrt elem x, y coords.
      REAL(DOUBLE)                    :: JACI(2,2)         ! An output from subr JAC2D, called herein. 2 x 2 Jacobian inverse.
      REAL(DOUBLE)                    :: DUM14(3,3)
      REAL(DOUBLE)                    :: DUM33(3,3)
      REAL(DOUBLE)                    :: JAC2x2(2,2)

      ! Thermal curvature change caused by uniform mid-surface scaling. Built by the contained subroutine
      ! THERM_CURV_AT_RS and used by both the thermal load vector (OPT(2)) and the stress recovery (OPT(3)).
      REAL(DOUBLE)                    :: KAP_TH(3)         ! Thermal curvature per unit temperature (kxx, kyy, kxy eng.)

! **********************************************************************************************************************************

! COORDINATE SYSTEMS
! ==================
!
! Basic
!  SNORM vector components are stored in this and transformed to element coordinates before use.
!
! Cartesian local
!  e^_1, e^_2, e^_3 in Bathe.
!  e^_3 is parallel to the director vector.
!  Used for strain in the strain-displacement matrix and the material elasticity matrix is transformed to this to integrate KE.
!  Orthogonal
!
! Element
!  x_element, y_element, z_element
!  Used for the grid point DOFs of the strain-displacement and the element stiffness matrices.
!  Used for extrapolating stress and strain from Gauss points to corners.
!  Used for element stress, strain, and force outputs
!  Defined the same way as MSC (element coordinate system).
!  Defined by x_element being the bisection of the diagonals and z_element being normal to both diagonals.
!  It's flat even when the element is warped.
!  Orthogonal
!
! Material
!  Used for material elasticity read from the input file.
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
!  Parallel to the isoparametric coordinates but scaled by the element size. Eg. |g_t| = half thickness in the direction of
!  the director vector (SNORM).
!  Not orthogonal
!
! Contravariant
!  g^r, g^s, g^t. g^1, g^2, g^3 in Bathe.
!  Contravariant to the covariant.
!  Not orthogonal
!
! V1, V2, Vn
! Used to express node rotations when building the strain-displacement matrix before being transformed to the element
! coordinate system.
! Vn is the director vector. V1 and V2 are in arbitrary orthogonal directions.
! Orthogonal

! **********************************************************************************************************************************

! Initialize
      PHI_SQ  = ONE                                        ! Not used for this element
      CALL MITC_INITIALIZE ()



! **********************************************************************************************************************************
! Generate the mass matrix for this element.

      IF (OPT(1) == 'Y') THEN

         ! Consistent mass matrix
         ! ME = ∫ N' ρ N det(J) dv

         M_1DOF(:,:) = ZERO

         DENSITY = MASS_PER_UNIT_AREA / EPROP(1)

         CALL ORDER_GAUSS ( IORD_IJ, SS_IJ, HH_IJ )
         CALL ORDER_GAUSS ( IORD_K, SS_K, HH_K )

                                                           ! Make mass matrix with 1 DOF per node
         DO I=1,IORD_IJ
            DO J=1,IORD_IJ
               DO K=1,IORD_K
                  R = SS_IJ(I)
                  S = SS_IJ(J)
                  T = SS_K(K)

                  CALL MITC_SHAPE_FUNCTIONS(R, S, PSH, DPSHG)

                  DETJ = MITC_DETJ ( R, S, T )
                  INTFAC = DETJ*HH_IJ(I)*HH_IJ(J)*HH_K(K)  ! det(J) * Gauss point weight

                  DO L=1,ELGP
                     DO M=1,ELGP
                        M_1DOF(L,M) = M_1DOF(L,M) + PSH(L) * PSH(M) * DENSITY * INTFAC
                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO

         CALL EXPAND_MASS_DOFS( M_1DOF )


      ENDIF



! **********************************************************************************************************************************
! Calculate element thermal loads.

      IF (OPT(2) == 'Y') THEN

! Thermal load, assembled using the ply-summed thermal force/moment resultants SHELL_AALP
! (membrane) and SHELL_BALP (membrane-bending coupling, nonzero for unsymmetric composite
! layups) from SHELL_ABD_MATRICES, instead of a homogeneous-material (R,S,T) volume Gauss
! integration -- same conversion, and same rationale, as the OPT(4) stiffness assembly above.
!
! Only a spatially-uniform temperature change (TBAR, averaged over the 4 corner grid points) is
! supported here, matching what this routine supported before this conversion; SHELL_DALP (the
! weight-z^2 thermal resultant that would respond to a through-thickness temperature *gradient*)
! is not used, since there is no gradient input to pair it with.
!
!    PTE = ( [Bm]' [A_alpha] + [Bb]' [B_alpha] ) * (TBAR - TREF)   integrated over mid-surface area
!
! On a curved reference surface the free thermal strain also contains the change of the generalized curvature produced by
! uniform mid-surface scaling. If the mid-surface is scaled by the free strain e = alpha*(TBAR - TREF) while the director field
! is held fixed, the mid-surface tangents become (1+e) x_,a, so the curvature x_,a . d_,b changes by e x_,a . d_,b. Without this
! term a uniformly heated cylinder can avoid artificial bending energy by opening its seam instead of expanding radially.
! The free transverse shear strain remains zero. Including the coupling to the complete section tangent gives
!
!    PTE = ( [Bm]' ( [A_alpha] + [B] {k_th} ) + [Bb]' ( [B_alpha] + [D] {k_th} ) ) * (TBAR - TREF)
!
! where {k_th} is the thermal curvature per unit temperature, {k_th} = (k_xx, k_yy, k_xy(engineering)), built here for a general
! in-plane free strain tensor {e_th} = [A]^-1 {A_alpha} (for an isotropic material {e_th} = alpha*(1,1,0), i.e. the scalar
! uniform mid-surface scaling above). With the free strain tensor e_th acting on the mid-surface tangents g_i (i = r,s), the
! T-linear covariant strain coefficient is
!
!    q_ij = 1/4 ( (e_th g_i) . w_j + (e_th g_j) . w_i ),     w_j = d(director*thickness)/d(xi_j) = g_j(T=+1) - g_j(T=-1)
!
! and it is carried to element coordinates and normalised by the thickness in exactly the same way that MITC4_BMBS builds the
! curvature operator (top minus bottom surface, each with the contravariant basis at its own T, divided by EPROP(1)), so that
! {k_th} is in the same basis and sign convention as [Bb]. For a flat element with the director normal to it, w_j = 0 and the
! term vanishes identically, so flat-plate thermal loads are unchanged.

         UNIT_PTE(:) = ZERO

         CALL ORDER_GAUSS ( IORD_IJ, SS_IJ, HH_IJ )

         DO I=1,IORD_IJ
            DO J=1,IORD_IJ
               R = SS_IJ(I)
               S = SS_IJ(J)

               CALL MITC4_BMBS( R, S, BMI3, BBI3, BSI2 )

                                                           ! Mid-surface area Jacobian (same construction used
                                                           ! for OPT(4) and the pressure load below).
               CALL MITC_COVARIANT_BASIS( R, S, ZERO, G )
               CALL CROSS( G(:,1), G(:,2), DUM4 )
               DETJ = SQRT(DOT_PRODUCT(DUM4, DUM4))
               INTFAC = DETJ*HH_IJ(I)*HH_IJ(J)

               UNIT_PTE(1:6*ELGP) = UNIT_PTE(1:6*ELGP)                                                                             &
                                   + MATMUL( TRANSPOSE(BMI3), SHELL_AALP ) * INTFAC                                                &
                                   + MATMUL( TRANSPOSE(BBI3), SHELL_BALP ) * INTFAC

                                                           ! Change of curvature from uniform mid-surface scaling.
                                                           ! KAP_TH is zero for a flat element whose director is normal to it.
               CALL THERM_CURV_AT_RS ( R, S, KAP_TH )

               UNIT_PTE(1:6*ELGP) = UNIT_PTE(1:6*ELGP)                                                                             &
                                   + MATMUL( TRANSPOSE(BMI3), MATMUL( SHELL_B, KAP_TH ) ) * INTFAC                                 &
                                   + MATMUL( TRANSPOSE(BBI3), MATMUL( SHELL_D, KAP_TH ) ) * INTFAC

            ENDDO
         ENDDO

                                                  ! Scale the thermal load vector by the temperature differences in each
                                                  ! subcase with thermal load.
         DO J=1,NTSUB
                                                  ! Use constant temperature to match the strain field of
                                                  ! linear elements.
            TBAR = (DT(1,J) + DT(2,J) + DT(3,J) + DT(4,J))/FOUR

            PTE(1:6*ELGP,J) = UNIT_PTE(1:6*ELGP) * (TBAR - TREF(1))
         ENDDO


      ENDIF

! **********************************************************************************************************************************


      IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         STR_PT_NUM = 0

         CALL ORDER_GAUSS ( IORD_STRESS_Q4, SS_IJ, HH_IJ )

                                                           ! Free thermal curvature per unit temperature at each stress point.
                                                           ! ELEM_STRE_STRN_ARRAYS subtracts it from the recovered curvature so
                                                           ! that the bending stress it reports is the mechanical part only, in
                                                           ! the same way ALPTM removes the free membrane thermal strain. Without
                                                           ! it a free, uniformly heated curved shell - which now deforms
                                                           ! correctly because OPT(2) puts the same term in the load vector -
                                                           ! would report a spurious bending stress equal to [D]{k_th}.
         THERM_KAPPA(:,:) = ZERO

         DO STR_PT_NUM = 1,5

                                                           ! Account for Bathe's R,S coordinates vs node numbering being different
                                                           ! from Mystran's
               SELECT CASE (STR_PT_NUM)
                  CASE (1); R=ZERO    ; S=ZERO             ! Center
                  CASE (2); R=SS_IJ(2); S=SS_IJ(2)         ! Gauss point 1
                  CASE (3); R=SS_IJ(2); S=SS_IJ(1)         ! Gauss point 2
                  CASE (4); R=SS_IJ(1); S=SS_IJ(2)         ! Gauss point 3
                  CASE (5); R=SS_IJ(1); S=SS_IJ(1)         ! Gauss point 4
               END SELECT

                                                           ! Get the mid-surface membrane, curvature, and transverse shear
                                                           ! strain-displacement operators at this (R,S), via the shared
                                                           ! top/bottom extraction in MITC4_BMBS (see that routine for the
                                                           ! rationale and the caveat about warped-element local bases).
               CALL MITC4_BMBS( R, S, BE1(1:3,1:6*ELGP,STR_PT_NUM), BE2(1:3,1:6*ELGP,STR_PT_NUM), BE3(1:2,1:6*ELGP,STR_PT_NUM) )

               IF (STR_PT_NUM <= MSTRPT) THEN
                  CALL THERM_CURV_AT_RS ( R, S, THERM_KAPPA(1:3,STR_PT_NUM) )
               ENDIF

         ENDDO


      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE.

      IF(OPT(4) == 'Y') THEN

! Element stiffness matrix, assembled using the A/B/D/T (classical-lamination-theory idealized)
! approach shared with the rest of MYSTRAN's shell elements (MIN4T's QDEL1/QPLT3 for QUAD4,
! TREL1/TPLT2 for TRIA3), instead of MITC4's original full (R,S,T) volume-integrated,
! single-homogeneous-material formulation.
!
! SHELL_A, SHELL_B, SHELL_D, SHELL_T are already populated for this element by
! SHELL_ABD_MATRICES (called once from EMG.f03 before this routine is reached), summed over
! PCOMP plies via classical lamination theory where applicable. Using them here (rather than
! re-deriving a homogeneous elasticity tensor and volume-integrating through T) is what removes
! the PCOMP restriction on MITC4/MITC4+: see the MITC4/MITC4+ ABD conversion plan.
!
!    KE = int over mid-surface area of:
!       [Bm]' [A] [Bm]  +  [Bm]' [B] [Bb] + [Bb]' [B]' [Bm]  +  [Bb]' [D] [Bb]  +  [Bs]' [T] [Bs]  dA
!
! Bm, Bb, Bs (membrane, curvature, and transverse-shear strain-displacement operators,
! evaluated at the mid-surface) come from MITC4_BMBS, which extracts them from the full 3-D
! MITC4_B by evaluating at the top/bottom surfaces (T = +-1) and combining -- exact, because
! MITC4_B's in-plane strain rows are linear in T and its shear rows are already independent of
! T. This preserves the MITC tying-point interpolation (the mechanism that avoids shear
! locking); only the through-thickness material integration changes, from numerical (a T Gauss
! loop applied to one homogeneous material) to analytic (SHELL_ABD_MATRICES, already correct
! for a ply stack).

         KE(1:6*ELGP,1:6*ELGP) = ZERO

         CALL ORDER_GAUSS ( IORD_IJ, SS_IJ, HH_IJ )

         DO I=1,IORD_IJ
            DO J=1,IORD_IJ
               R = SS_IJ(I)
               S = SS_IJ(J)

               CALL MITC4_BMBS( R, S, BMI3, BBI3, BSI2 )

                                                           ! Mid-surface area Jacobian, same construction as used
                                                           ! for the pressure load below: dA = |g_r x g_s| dR dS.
               CALL MITC_COVARIANT_BASIS( R, S, ZERO, G )
               CALL CROSS( G(:,1), G(:,2), DUM4 )
               DETJ = SQRT(DOT_PRODUCT(DUM4, DUM4))
               INTFAC = DETJ*HH_IJ(I)*HH_IJ(J)

                                                           ! Membrane
                                                           ! ∫ Bm' A Bm dA
               CALL MATMULT_FFF   ( SHELL_A, BMI3, 3, 3, 6*ELGP, DUM1_3 )
               CALL MATMULT_FFF_T ( BMI3, DUM1_3, 3, 6*ELGP, 6*ELGP, DUM2 )
               KE(1:6*ELGP,1:6*ELGP) = KE(1:6*ELGP,1:6*ELGP) + DUM2(:,:)*INTFAC

                                                           ! Bending
                                                           ! ∫ Bb' D Bb dA
               CALL MATMULT_FFF   ( SHELL_D, BBI3, 3, 3, 6*ELGP, DUM1_3 )
               CALL MATMULT_FFF_T ( BBI3, DUM1_3, 3, 6*ELGP, 6*ELGP, DUM2 )
               KE(1:6*ELGP,1:6*ELGP) = KE(1:6*ELGP,1:6*ELGP) + DUM2(:,:)*INTFAC

                                                           ! Transverse shear
                                                           ! ∫ Bs' T Bs dA
               CALL MATMULT_FFF   ( SHELL_T, BSI2, 2, 2, 6*ELGP, DUM1_2 )
               CALL MATMULT_FFF_T ( BSI2, DUM1_2, 2, 6*ELGP, 6*ELGP, DUM2 )
               KE(1:6*ELGP,1:6*ELGP) = KE(1:6*ELGP,1:6*ELGP) + DUM2(:,:)*INTFAC

                                                           ! Membrane-bending coupling
                                                           ! Zero for a symmetric PSHELL/PCOMP layup with no offset;
                                                           ! nonzero for offset shells, unsymmetric composite layups,
                                                           ! or (if enabled) PSHELL MID4 coupling.
                                                           ! ∫ Bm' B Bb dA + ∫ Bb' B' Bm dA
               CALL MATMULT_FFF   ( SHELL_B, BBI3, 3, 3, 6*ELGP, DUM1_3 )
               CALL MATMULT_FFF_T ( BMI3, DUM1_3, 3, 6*ELGP, 6*ELGP, DUM2 )
               KE(1:6*ELGP,1:6*ELGP) = KE(1:6*ELGP,1:6*ELGP) + DUM2(:,:)*INTFAC
               CALL MATMULT_FFF_T ( SHELL_B, BMI3, 3, 3, 6*ELGP, DUM1_3 )
               CALL MATMULT_FFF_T ( BBI3, DUM1_3, 3, 6*ELGP, 6*ELGP, DUM2 )
               KE(1:6*ELGP,1:6*ELGP) = KE(1:6*ELGP,1:6*ELGP) + DUM2(:,:)*INTFAC

            ENDDO
         ENDDO


      ENDIF


! **********************************************************************************************************************************
! Determine element pressure loads

      IF (OPT(5) == 'Y') THEN

         UNIT_PPE(:) = ZERO

         CALL ORDER_GAUSS ( IORD_IJ, SS_IJ, HH_IJ )

         DO I=1,IORD_IJ
            DO J=1,IORD_IJ
               R = SS_IJ(I)
               S = SS_IJ(J)

               CALL MITC_SHAPE_FUNCTIONS(R, S, PSH, DPSHG)
                                                           ! Normalized normal vector at R,S
                                                           ! This is the interpolated director vector
                                                           ! and follows SNORM if specified.
               CALL MITC_COVARIANT_BASIS( R, S, ZERO, G )
               NORMAL(:) = G(:,3) / SQRT(DOT_PRODUCT(G(:,3), G(:,3)))

               CALL CROSS(G(:,1),G(:,2),DUM4)
               DETJ = SQRT(DOT_PRODUCT(DUM4, DUM4))
               INTFAC = DETJ * HH_IJ(I) * HH_IJ(J)         ! Contribution to area for this Gauss point

               DO GP=1,ELGP
                  K = (GP-1) * 6
                  UNIT_PPE(K+1:K+3) = UNIT_PPE(K+1:K+3) + NORMAL * PSH(GP) * INTFAC
               ENDDO

            ENDDO
         ENDDO

                                                           ! Scale the unit pressure load vector by the pressure in each subcase.
         DO J=1,NSUB
            PPE(1:6*ELGP,J) = UNIT_PPE(1:6*ELGP) * PRESS(3,J)
         ENDDO



      ENDIF

! **********************************************************************************************************************************
! Calculate linear differential stiffness matrix

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

                                                           ! Find membrane engineering forces at each Gauss point
                                                           ! by summing the forces in each ply.
         FORCEx(:)  = 0
         FORCEy(:)  = 0
         FORCExy(:) = 0

         CALL GET_ELEM_NUM_PLIES ( INT_ELEM_ID )           ! Get NUM_PLIES

         DO JPLY=1,NUM_PLIES
                                                           ! Get UEL, EM, TPLY for this ply.
            IF (PCOMP_PROPS == 'N') THEN
                                                           ! EM is already set above.

               TPLY = EPROP(1)                             ! Element thickness
               CALL ELMDIS

            ELSE


               IF (PCOMP_LAM == 'NON') THEN                ! Delete this IF block to allow the LAM field set to nonsymmetric.
                  FATAL_ERR          = FATAL_ERR + 1
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  WRITE(ERR,*) ' *ERROR: Code not written for non-symmetric composite buckling or differential stiffness.'
                  WRITE(F06,*) ' *ERROR: Code not written for non-symmetric composite buckling or differential stiffness.'
                  CALL OUTA_HERE ( 'Y' )
               ENDIF

               PLY_NUM = JPLY                              ! Used by SHELL_ABD_MATRICES
               CALL SHELL_ABD_MATRICES ( INT_ELEM_ID, 'N' )! Get EM, ZPLY, TPLY, ALPVEC for this ply
               CALL ELMDIS
               CALL ELMDIS_PLY                             ! Adjust UEL using ZPLY

            ENDIF

            GAUSS_PT = 0
            DO I=1,IORD_STRESS_Q4
               DO J=1,IORD_STRESS_Q4
                  GAUSS_PT = GAUSS_PT + 1

                                                           ! Stress at this Gauss point using UEL, BE1, EM, ALPVEC, DT
                  CALL ELEM_STRE_STRN_ARRAYS ( GAUSS_PT+1 )

                  FORCEx(GAUSS_PT)  = FORCEx(GAUSS_PT)  + TPLY*STRESS(1)
                  FORCEy(GAUSS_PT)  = FORCEy(GAUSS_PT)  + TPLY*STRESS(2)
                  FORCExy(GAUSS_PT) = FORCExy(GAUSS_PT) + TPLY*STRESS(3)
               ENDDO
            ENDDO

         ENDDO


                                                           ! Transform force from element coordinates to cartesian local
                                                           ! This isn't quite right for non-flat elements becuase all the z
                                                           ! components are omitted in both coordinate systems.
                                                           ! It would be better to obtain force directly in cartesian local
                                                           ! coordinates instead.
         CALL ORDER_GAUSS ( IORD_STRESS_Q4, SSS, HHH )

         DO GAUSS_PT = 1,4
                                                           ! 3x3 force tensor
            DUM14(:,:) = ZERO
            DUM14(1,1) = FORCEx(GAUSS_PT)
            DUM14(2,2) = FORCEy(GAUSS_PT)
            DUM14(1,2) = FORCExy(GAUSS_PT)
            DUM14(2,1) = FORCExy(GAUSS_PT)

            SELECT CASE (GAUSS_PT)
               CASE (1); R=SSS(2); S=SSS(2)             ! Gauss point 1
               CASE (2); R=SSS(2); S=SSS(1)             ! Gauss point 2
               CASE (3); R=SSS(1); S=SSS(2)             ! Gauss point 3
               CASE (4); R=SSS(1); S=SSS(1)             ! Gauss point 4
            END SELECT

            CLB = MITC4_CARTESIAN_LOCAL_BASIS( R, S, ZERO )

            CALL MATMULT_FFF (DUM14, CLB, 3, 3, 3, DUM33 )
            CALL MATMULT_FFF_T (CLB, DUM33, 3, 3, 3, DUM14 )

            FORCEx(GAUSS_PT) = DUM14(1,1)
            FORCEy(GAUSS_PT) = DUM14(2,2)
            FORCExy(GAUSS_PT) = DUM14(1,2)
         ENDDO


! Accoring to:
!   Robert D. Cook, David S. Malkus, Michael E. Plesha Concepts and Applications of Finite Element Analysis, 3rd Edition  1989
!   Section 14.3 Stress Stiffness Matrix Of A Plate Element

!         +1  +1
! [  ]   ⌠   ⌠  [   ]T [   ]-T [ Nx  Nxy ] [   ]-1 [   ]
! [k ] = |   |  [ G ]  [ J ]   [ Nxy Ny  ] [ J ]   [ G ]  |J|  dξ dη
! [ σ]   ⌡   ⌡  [  I]  [   ]               [   ]   [  I]
!        -1  -1

! k_σ is the stress stiffness (differential stiffness) matrix
! Nx, Ny, Nxy are membrane engineering forces
! |J| is the Jacobian determinant
! G_I is a 2xELGP matrix of shape function derivatives with respect to isoparametric coordinates ξ  and η.

! DPSHX = J^-1 G_I is the 2 x ELGP matrix of shape function derivatives with respect to element coordinates x and y.

!        +1  +1
! [k ]   ⌠   ⌠        T [ Nx  Nxy ]
! [ σ] = ⌡   ⌡   DPSHX  [ Nxy Ny  ]  DPSHX  |J|  dξ dη
!        -1  -1

         KS(:,:) = ZERO



         GAUSS_PT = 0
         DO I=1,IORD_STRESS_Q4
            DO J=1,IORD_STRESS_Q4
               GAUSS_PT = GAUSS_PT + 1

                                                           ! Account for Bathe's R,S coordinates vs node numbering being different
                                                           ! from Mystran's
               SELECT CASE (GAUSS_PT)
                  CASE (1); R=SSS(2); S=SSS(2)             ! Gauss point 1
                  CASE (2); R=SSS(2); S=SSS(1)             ! Gauss point 2
                  CASE (3); R=SSS(1); S=SSS(2)             ! Gauss point 3
                  CASE (4); R=SSS(1); S=SSS(1)             ! Gauss point 4
               END SELECT

               DUM11(1,1) = FORCEx(GAUSS_PT)  ; DUM11(1,2) = FORCExy(GAUSS_PT)
               DUM11(2,1) = FORCExy(GAUSS_PT) ; DUM11(2,2) = FORCEy(GAUSS_PT)


                                                           ! 2D inverse Jacobian
               CALL MITC_COVARIANT_BASIS ( R, S, ZERO, G )
               CLB = MITC4_CARTESIAN_LOCAL_BASIS( R, S, ZERO )
               DUM14 = MATMUL(TRANSPOSE(G), CLB)
               JAC2x2(1:2,1:2) = DUM14(1:2,1:2)
               DETJ = JAC2x2(1,1)*JAC2x2(2,2) - JAC2x2(1,2)*JAC2x2(2,1)
               JACI(1,1) =  JAC2x2(2,2)/DETJ
               JACI(1,2) = -JAC2x2(1,2)/DETJ
               JACI(2,1) = -JAC2x2(2,1)/DETJ
               JACI(2,2) =  JAC2x2(1,1)/DETJ
               CALL MITC_SHAPE_FUNCTIONS ( R, S, PSH, DPSHG )
                                                           ! Shape function derivatives at this Gauss point.
               CALL MATMULT_FFF ( JACI, DPSHG, 2, 2, 4, DPSHX )

               CALL MATMULT_FFF_T ( DPSHX, DUM11, 2, ELGP, 2, DUM12 )
               CALL MATMULT_FFF ( DUM12, DPSHX, ELGP, 2, ELGP, DUM13 )

               INTFAC = DETJ*HHH(I)*HHH(J)

                                                           ! Accumulate integrand into the result
               DO K=1,ELGP
                  DO L=1,ELGP
                     KS(K,L) = KS(K,L) + DUM13(K,L) * INTFAC
                  ENDDO
               ENDDO

            ENDDO
         ENDDO


                                                           ! Copy KS into KED for each translational DOF.
         KED(1:6*ELGP,1:6*ELGP) = 0
         DO I=1,ELGP
            DO J=1,ELGP
               KI = (I-1) * 6
               KJ = (J-1) * 6
               KED(KI + 1, KJ + 1) = KS(I,J)
               KED(KI + 2, KJ + 2) = KS(I,J)
               KED(KI + 3, KJ + 3) = KS(I,J)
            ENDDO
         ENDDO


      ENDIF




      RETURN

! **********************************************************************************************************************************

      CONTAINS

! **********************************************************************************************************************************

      SUBROUTINE THERM_CURV_AT_RS ( RR, SS, KAP )

! Free thermal curvature per unit temperature at one point of the mid-surface, in element coordinates, in the same basis and sign
! convention as the curvature operator [Bb] returned by MITC4_BMBS.
!
! On a curved reference surface the free thermal strain contains a change of the generalized curvature produced by uniform
! mid-surface scaling. If the mid-surface is scaled by the free strain e = alpha*(TBAR - TREF) while the director field is held
! fixed, the mid-surface tangents become (1+e) x_,a, so the curvature x_,a . d_,b changes by e x_,a . d_,b.
!
! With the free strain tensor e_th = [A]^-1 {A_alpha} acting on the mid-surface tangents g_i (i = r,s), the T-linear covariant
! strain coefficient is
!
!    q_ij = 1/4 ( (e_th g_i) . w_j + (e_th g_j) . w_i ),     w_j = d(director*thickness)/d(xi_j) = g_j(T=+1) - g_j(T=-1)
!
! carried to element coordinates and normalised by the thickness exactly as MITC4_BMBS builds the curvature operator (top and
! bottom surface, each with the contravariant basis at its own T, divided by EPROP(1)).
!
! For a flat element whose director is normal to it, w_j = 0 and KAP is identically zero, so flat plates are unaffected. KAP is
! also returned zero if [A] is singular.

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: RR, SS            ! Isoparametric coordinates of the point
      REAL(DOUBLE), INTENT(OUT)       :: KAP(3)            ! Thermal curvature per unit temperature (kxx, kyy, kxy engineering)

      INTEGER(LONG)                   :: IA, IB, ID        ! Tensor / DO loop indices
      INTEGER(LONG)                   :: KK, LL            ! Covariant in-plane indices

      REAL(DOUBLE)                    :: SIDE              ! T = +1 or -1 surface used to build the curvature operator
      REAL(DOUBLE)                    :: DETA              ! Determinant of SHELL_A
      REAL(DOUBLE)                    :: ADJA(3,3)         ! Adjugate of SHELL_A
      REAL(DOUBLE)                    :: EPS_TH(3)         ! Free mid-surface thermal strain per unit temperature (xx,yy,xy eng.)
      REAL(DOUBLE)                    :: EPS_TEN(3,3)      ! EPS_TH as a 3x3 tensor in element coordinates
      REAL(DOUBLE)                    :: GMID(3,3)         ! Covariant basis at T = 0
      REAL(DOUBLE)                    :: GTOP(3,3)         ! Covariant basis at T = +1
      REAL(DOUBLE)                    :: GBOT(3,3)         ! Covariant basis at T = -1
      REAL(DOUBLE)                    :: GSIDE(3,3)        ! Covariant basis at T = SIDE
      REAL(DOUBLE)                    :: GCON(3,3)         ! Contravariant basis at T = SIDE
      REAL(DOUBLE)                    :: DIRW(3,2)         ! d(director*thickness)/dR, /dS
      REAL(DOUBLE)                    :: EPSX(3,2)         ! EPS_TEN * g_r, EPS_TEN * g_s
      REAL(DOUBLE)                    :: QTH(2,2)          ! Covariant T-linear thermal strain coefficient
      REAL(DOUBLE)                    :: KTEN(3,3)         ! Thermal curvature tensor in element coordinates

! **********************************************************************************************************************************

      KAP(:) = ZERO

                                                           ! Free mid-surface thermal strain per unit temperature:
                                                           ! {e_th} = [A]^-1 {A_alpha}, via the adjugate of the 3x3 SHELL_A.
      DETA = SHELL_A(1,1)*(SHELL_A(2,2)*SHELL_A(3,3) - SHELL_A(2,3)*SHELL_A(3,2))                                                  &
           - SHELL_A(1,2)*(SHELL_A(2,1)*SHELL_A(3,3) - SHELL_A(2,3)*SHELL_A(3,1))                                                  &
           + SHELL_A(1,3)*(SHELL_A(2,1)*SHELL_A(3,2) - SHELL_A(2,2)*SHELL_A(3,1))

      IF (DABS(DETA) <= 1.0D-12*DABS(SHELL_A(1,1)*SHELL_A(2,2)*SHELL_A(3,3))) THEN
         RETURN
      ENDIF

      ADJA(1,1) = SHELL_A(2,2)*SHELL_A(3,3) - SHELL_A(2,3)*SHELL_A(3,2)
      ADJA(1,2) = SHELL_A(1,3)*SHELL_A(3,2) - SHELL_A(1,2)*SHELL_A(3,3)
      ADJA(1,3) = SHELL_A(1,2)*SHELL_A(2,3) - SHELL_A(1,3)*SHELL_A(2,2)
      ADJA(2,1) = SHELL_A(2,3)*SHELL_A(3,1) - SHELL_A(2,1)*SHELL_A(3,3)
      ADJA(2,2) = SHELL_A(1,1)*SHELL_A(3,3) - SHELL_A(1,3)*SHELL_A(3,1)
      ADJA(2,3) = SHELL_A(1,3)*SHELL_A(2,1) - SHELL_A(1,1)*SHELL_A(2,3)
      ADJA(3,1) = SHELL_A(2,1)*SHELL_A(3,2) - SHELL_A(2,2)*SHELL_A(3,1)
      ADJA(3,2) = SHELL_A(1,2)*SHELL_A(3,1) - SHELL_A(1,1)*SHELL_A(3,2)
      ADJA(3,3) = SHELL_A(1,1)*SHELL_A(2,2) - SHELL_A(1,2)*SHELL_A(2,1)

      EPS_TH(1:3) = MATMUL( ADJA, SHELL_AALP ) / DETA

      EPS_TEN(:,:) = ZERO                                  ! Tensor form of {e_th} in element coordinates (engineering xy -> /2)
      EPS_TEN(1,1) = EPS_TH(1)
      EPS_TEN(2,2) = EPS_TH(2)
      EPS_TEN(1,2) = HALF*EPS_TH(3)
      EPS_TEN(2,1) = HALF*EPS_TH(3)

      CALL MITC_COVARIANT_BASIS( RR, SS,  ZERO, GMID )
      CALL MITC_COVARIANT_BASIS( RR, SS,  +ONE, GTOP )
      CALL MITC_COVARIANT_BASIS( RR, SS,  -ONE, GBOT )

      DO IA=1,2
         DIRW(:,IA) = GTOP(:,IA) - GBOT(:,IA)
         EPSX(:,IA) = MATMUL( EPS_TEN, GMID(:,IA) )
      ENDDO

      DO IA=1,2
         DO IB=1,2
            QTH(IA,IB) = QUARTER*( DOT_PRODUCT( EPSX(:,IA), DIRW(:,IB) ) + DOT_PRODUCT( EPSX(:,IB), DIRW(:,IA) ) )
         ENDDO
      ENDDO

                                                           ! Top minus bottom surface, as in MITC4_BMBS. The covariant
                                                           ! coefficient is +Q at T=+1 and -Q at T=-1, so the two surface
                                                           ! contributions add.
      KTEN(:,:) = ZERO
      DO ID=1,2
         SIDE = ONE
         IF (ID == 2) SIDE = -ONE
         CALL MITC_COVARIANT_BASIS( RR, SS, SIDE, GSIDE )
         CALL MITC_CONTRAVARIANT_BASIS( GSIDE, GCON )
         DO IA=1,3
            DO IB=1,3
               DO KK=1,2
                  DO LL=1,2
                     KTEN(IA,IB) = KTEN(IA,IB) + GCON(IA,KK)*GCON(IB,LL)*QTH(KK,LL)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      KAP(1) =     KTEN(1,1)/EPROP(1)
      KAP(2) =     KTEN(2,2)/EPROP(1)
      KAP(3) = TWO*KTEN(1,2)/EPROP(1)                      ! Engineering xy component, same as row 3 of BBI3

      RETURN

      END SUBROUTINE THERM_CURV_AT_RS

! **********************************************************************************************************************************

      END SUBROUTINE MITC4

      SUBROUTINE MITC4_B ( R, S, T, MEMBRANE, BENDING, SHEAR, B )

! Calculates the strain-displacement matrix in the cartesian local coordinate system
! for MITC4 shell at one point in isoparametric coordinates.
!
! Reference [1]:
!  "A new MITC4+ shell element" by Ko, Lee, Bathe, 2016
!
! Reference [2]:
!  MITC4 paper "A continuum mechanics based four-node shell element for general nonlinear analysis"
!     by Dvorkin and Bathe


      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP, XEL
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, ONE, QUARTER, TWO, FOUR
      USE PARAMS, ONLY                :  QUAD4TYP
      USE MITC_STUF, ONLY             :  GP_RS

      USE VECTOR_GEOMETRY, ONLY       :  CROSS

      IMPLICIT NONE

      INTEGER(LONG)                   :: COL               ! A column (element DOF) of B. 1-24.
      INTEGER(LONG)                   :: GP                ! Grid point number. 1-4.

      REAL(DOUBLE) , INTENT(IN)       :: R, S, T           ! Isoparametric coordinates
      REAL(DOUBLE) , INTENT(OUT)      :: B(6, 6*ELGP)      ! Strain-displacement matrix
      REAL(DOUBLE)                    :: X_R(3)            ! Characteristic geometry vector x_r
      REAL(DOUBLE)                    :: X_S(3)            ! Characteristic geometry vector x_s
      REAL(DOUBLE)                    :: X_D(3)            ! Characteristic geometry vector x_d (distortion vector)
      REAL(DOUBLE)                    :: BM(6, 6*ELGP)     ! Strain-displacement matrix for membrane
      REAL(DOUBLE)                    :: BB(6, 6*ELGP)     ! Strain-displacement matrix for bending
      REAL(DOUBLE)                    :: BS(6, 6*ELGP)     ! Strain-displacement matrix for shear
      REAL(DOUBLE)                    :: BM_A(6, 6*ELGP)
      REAL(DOUBLE)                    :: BM_B(6, 6*ELGP)
      REAL(DOUBLE)                    :: BM_C(6, 6*ELGP)
      REAL(DOUBLE)                    :: BM_D(6, 6*ELGP)
      REAL(DOUBLE)                    :: BM_E(6, 6*ELGP)
      REAL(DOUBLE)                    :: E(6, 6*ELGP)      ! Strain-displacement matrix directly interpolated
      REAL(DOUBLE)                    :: BS_A(6, 6*ELGP)
      REAL(DOUBLE)                    :: BS_B(6, 6*ELGP)
      REAL(DOUBLE)                    :: BS_C(6, 6*ELGP)
      REAL(DOUBLE)                    :: BS_D(6, 6*ELGP)
      REAL(DOUBLE)                    :: XRxXS(3)
      REAL(DOUBLE)                    :: MR(3)
      REAL(DOUBLE)                    :: MS(3)
      REAL(DOUBLE)                    :: DUM1(3)
      REAL(DOUBLE)                    :: c_r, c_s, d       ! Distortion variables used in ref [1]
                                                           ! Intermediate variables used in ref [1]
      REAL(DOUBLE)                    :: a_A, a_B, a_C, a_D, a_E

      LOGICAL      , INTENT(IN)       :: MEMBRANE          ! If true, generate membrane parts of B (rows 1,2,4)
      LOGICAL      , INTENT(IN)       :: BENDING           ! If true, generate bending parts of B (rows 1,2,4)
      LOGICAL      , INTENT(IN)       :: SHEAR             ! If true, generate shear parts of B (rows 5,6)


! **********************************************************************************************************************************
! Initialize empty matrix

      B(:,:) = ZERO


! **********************************************************************************************************************************
! Add in-layer strain-displacement terms

                                                           ! Characteristic geometry vectors
      X_R(:) = ZERO
      X_S(:) = ZERO
      X_D(:) = ZERO
      DO GP=1,ELGP
         X_R(:) = X_R(:) + QUARTER * GP_RS(1, GP)                * XEL(GP, :)
         X_S(:) = X_S(:) + QUARTER *                GP_RS(2, GP) * XEL(GP, :)
         X_D(:) = X_D(:) + QUARTER * GP_RS(1, GP) * GP_RS(2, GP) * XEL(GP, :)
      ENDDO


      IF(QUAD4TYP == 'MITC4+') THEN
                                                           ! MITC4+ according to ref [1]

         IF(MEMBRANE) THEN
                                                           ! BM at each membrane strain tying point.
            !
            !Membrane strain tying points A,B,C,D,E
            ! 2     A     1
            !  +----o----+
            !  |    ^s   |
            !  |    |    |
            !D o    +->r o C
            !  |   E     |
            !  |         |
            !  +----o----+
            ! 3     B     4
            !
            CALL MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION( ZERO,  ONE , T, X_R, X_S, X_D, .TRUE., .FALSE., 1, 1, BM_A )
            CALL MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION( ZERO, -ONE , T, X_R, X_S, X_D, .TRUE., .FALSE., 1, 1, BM_B )
            CALL MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION( ONE ,  ZERO, T, X_R, X_S, X_D, .TRUE., .FALSE., 2, 2, BM_C )
            CALL MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION(-ONE ,  ZERO, T, X_R, X_S, X_D, .TRUE., .FALSE., 2, 2, BM_D )
            CALL MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION( ZERO,  ZERO, T, X_R, X_S, X_D, .TRUE., .FALSE., 4, 4, BM_E )

                                                           ! Dual basis vectors m^r and m^s to the characteristic
                                                           ! geometry vectors x_r and x_s. From eqn (11) in ref [1].
            CALL CROSS(X_R, X_S, XRxXS)
            CALL CROSS(X_S, XRxXS, DUM1)
            MR = DUM1 / DOT_PRODUCT(X_R, DUM1)
            CALL CROSS(XRxXS, X_R, DUM1)
            MS = DUM1 / DOT_PRODUCT(X_S, DUM1)

                                                           ! c_r, c_s, d from eqn (24) in ref [1].
            c_r = DOT_PRODUCT(X_D, MR)
            c_s = DOT_PRODUCT(X_D, MS)
            d = c_r * c_r + c_s * c_s - ONE

            a_A = c_r * (c_r - 1) / (TWO * d)
            a_B = c_r * (c_r + 1) / (TWO * d)
            a_C = c_s * (c_s - 1) / (TWO * d)
            a_D = c_s * (c_s + 1) / (TWO * d)
            a_E = 2 * c_r * c_s / d

                                                           ! Eqn (27a) in ref [1]
            BM(1,:) = HALF * (ONE - TWO * a_A + S + 2 * a_A * S*S ) * BM_A(1,:)                                                    &
                    + HALF * (ONE - TWO * a_B - S + 2 * a_B * S*S ) * BM_B(1,:)                                                    &
                    + a_C * (-ONE + S*S) * BM_C(2,:)                                                                               &
                    + a_D * (-ONE + S*S) * BM_D(2,:)                                                                               &
                    + a_E * (-ONE + S*S) * BM_E(4,:)
                                                           ! Eqn (27b) in ref [1]
            BM(2,:) = a_A * (-ONE + R*R) * BM_A(1,:)                                                                               &
                    + a_B * (-ONE + R*R) * BM_B(1,:)                                                                               &
                    + HALF * (ONE - TWO * a_C + R + 2 * a_C * R*R ) * BM_C(2,:)                                                    &
                    + HALF * (ONE - TWO * a_D - R + 2 * a_D * R*R ) * BM_D(2,:)                                                    &
                    + a_E * (-ONE + R*R) * BM_E(4,:)

            BM(3,:) = ZERO
                                                           ! Eqn (27c) in ref [1]
            BM(4,:) = QUARTER * ( R + FOUR * a_A * R * S) * BM_A(1,:)                                                              &
                    + QUARTER * (-R + FOUR * a_B * R * S) * BM_B(1,:)                                                              &
                    + QUARTER * ( S + FOUR * a_C * R * S) * BM_C(2,:)                                                              &
                    + QUARTER * (-S + FOUR * a_D * R * S) * BM_D(2,:)                                                              &
                    + (1 + a_E * R * S) * BM_E(4,:)

            B(1:4,:) = B(1:4,:) + BM(1:4,:)

         ENDIF

         IF(BENDING) THEN
                                                          ! Bending is the same as the MITC4+ form of MITC4
            CALL MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION( R, S, T, X_R, X_S, X_D, .FALSE., .TRUE., 1, 4, BB )
            B(1:4,:) = B(1:4,:) + BB(1:4,:)
         ENDIF

      ELSEIF(QUAD4TYP == 'MITC4 ') THEN

         IF(.TRUE.) THEN

            IF(MEMBRANE) THEN
                                                           ! MITC4+ form of MITC4 according to ref [1]
               CALL MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION( R, S, T, X_R, X_S, X_D, .TRUE., .FALSE., 1, 4, BM )
               B(1:4,:) = B(1:4,:) + BM(1:4,:)
            ENDIF

            IF(BENDING) THEN
               CALL MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION( R, S, T, X_R, X_S, X_D, .FALSE., .TRUE., 1, 4, BB )
               B(1:4,:) = B(1:4,:) + BB(1:4,:)
            ENDIF

         ELSE

                                                           ! MITC4 according to ref [2]
                                                           ! Equivalent to the MITC4+ form of MITC4
                                                           ! but can't separate membrane and bending.
                                                           ! Could be removed and this branch is never reached.
            IF(MEMBRANE .AND. BENDING) THEN

               CALL MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION( R, S, T, 1, 4, E )
               B(1:4,:) = B(1:4,:) + E(1:4,:)

            ENDIF

         ENDIF

      ENDIF



! **********************************************************************************************************************************
! Add transverse shear strain-displacement terms

      IF(SHEAR) THEN

         ! According to ref [2]. Tying point labels are different from ref [1] but it's otherwise equivalent.
         ! The same in MITC4 and MITC4+.

         !
         !Tying points A,B,C,D are the same as in Bathe wrt R and S (Bathe's r_1 and r_2) and same node numbering:
         ! 2     A     1
         !  +----o----+
         !  |    ^s   |
         !  |    |    |
         !B o    +->r o D
         !  |         |
         !  |         |
         !  +----o----+
         ! 3     C     4
         !

         CALL MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION( ZERO, ONE,  ZERO, 5, 6, BS_A )
         CALL MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION(-ONE,  ZERO, ZERO, 5, 6, BS_B )
         CALL MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION( ZERO,-ONE,  ZERO, 5, 6, BS_C )
         CALL MITC_COVARIANT_STRAIN_DIRECT_INTERPOLATION( ONE,  ZERO, ZERO, 5, 6, BS_D )

         DO COL=1,6*ELGP
           !e_st
           B(5, COL) = HALF * (ONE + R) * BS_D(5, COL) + HALF * (ONE - R) * BS_B(5, COL)
           !e_rt
           B(6, COL) = HALF * (ONE + S) * BS_A(6, COL) + HALF * (ONE - S) * BS_C(6, COL)
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Transform covariant strain components from the contravariant basis to the cartesian local basis.

      CALL MITC_TRANSFORM_CONTRAVARIANT_TO_LOCAL( R, S, T, B )


      ! Double shear terms because it's now treated as vectors instead of tensors.
      B(4:6,:) = B(4:6,:) * 2

      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE MITC4_B

SUBROUTINE MITC4_BMBS ( R, S, BM, BB, BS )

  USE PENTIUM_II_KIND, ONLY   : DOUBLE
  USE MODEL_STUF, ONLY    : ELGP, EPROP
  USE CONSTANTS_1, ONLY       : ZERO, HALF, ONE


  IMPLICIT NONE

  REAL(DOUBLE), INTENT(IN)    :: R, S
  REAL(DOUBLE), INTENT(OUT)   :: BM(3, 6*ELGP)
  REAL(DOUBLE), INTENT(OUT)   :: BB(3, 6*ELGP)
  REAL(DOUBLE), INTENT(OUT)   :: BS(2, 6*ELGP)

  REAL(DOUBLE)        :: BMEM(6, 6*ELGP)
  REAL(DOUBLE)        :: BBOT(6, 6*ELGP)
  REAL(DOUBLE)        :: BTOP(6, 6*ELGP)
  REAL(DOUBLE)        :: BSHR(6, 6*ELGP)
  REAL(DOUBLE)        :: G(3,3)            ! Covariant basis at the mid-surface, in element coordinates
  REAL(DOUBLE)        :: M1, M2            ! In-plane slope of the director relative to the element facet

! **********************************************************************************************************************************
! Subr MITC4_B returns the strain-displacement matrix in the CARTESIAN LOCAL basis, whose x axis lies along the covariant g_r
! direction. Everything that uses the operators returned here - the elasticity matrices built from EM, EB and ET, the stress,
! strain and engineering force output, and the extrapolation of Gauss point values to the corners - works in the ELEMENT
! coordinate system. So the rows have to be rotated out of the cartesian local basis before they are handed back.
!
! For a rectangular element the two systems differ by 180 degrees about z, because the element x axis starts along side 1-2
! while the cartesian local x axis lies along g_r, and side 1-2 runs in the negative r direction in Bathe's node ordering.
!
! MITC4_B doubles rows 4 to 6 on the way out to make them engineering shear strains. MITC_TRANSFORM_B rotates tensor
! components, so those rows are halved before the rotation and doubled again afterwards. This matches what subr MITC4 on the
! dev branch does at each of these points.

! **********************************************************************************************************************************
! Pure midsurface membrane operator.

  CALL MITC4_B( R, S, ZERO, .TRUE., .FALSE., .FALSE., BMEM )
  CALL TO_ELEMENT_BASIS( ZERO, BMEM )

  BM(1,:) = BMEM(1,:)  ! xx
  BM(2,:) = BMEM(2,:)  ! yy
  BM(3,:) = BMEM(4,:)  ! xy

! **********************************************************************************************************************************
! Pure curvature operator.
! Use bending-only calls so membrane terms cannot leak into curvature.
! The difference removes even-in-T bending terms.

  CALL MITC4_B( R, S, -ONE, .FALSE., .TRUE., .FALSE., BBOT )
  CALL TO_ELEMENT_BASIS( -ONE, BBOT )

  CALL MITC4_B( R, S, +ONE, .FALSE., .TRUE., .FALSE., BTOP )
  CALL TO_ELEMENT_BASIS( +ONE, BTOP )

  BB(1,:) = (BTOP(1,:) - BBOT(1,:)) / EPROP(1)  ! kxx-ish
  BB(2,:) = (BTOP(2,:) - BBOT(2,:)) / EPROP(1)  ! kyy-ish
  BB(3,:) = (BTOP(4,:) - BBOT(4,:)) / EPROP(1)  ! kxy-ish

! **********************************************************************************************************************************
! Pure transverse shear operator.
! MITC shear is already evaluated at tying points and is effectively midsurface shear here.
! Keep row order as zx, yz to match SHELL_T convention used by MITC4.f03.

  CALL MITC4_B( R, S, ZERO, .FALSE., .FALSE., .TRUE., BSHR )
  CALL TO_ELEMENT_BASIS( ZERO, BSHR )

  BS(1,:) = BSHR(6,:)  ! zx
  BS(2,:) = BSHR(5,:)  ! yz



! Shear coupling correction removed because it worsens curved shell thermal loads.
! Parallelogram elements wouldn't happen naturally with averaged normals except at
! the transition between +ve and -ve curvature which is only a 1D line so the effect
! should vanish with mesh refinement.

! **********************************************************************************************************************************
! Remove the spurious membrane to transverse shear coupling that a director which is not normal to the element facet produces.
!
! The reference geometry of the degenerated shell is X(r,s,t) = Xbar(r,s) + (t*h/2) * V, so when the director V is not normal to
! the facet the through-thickness fibre is slanted and a point at height z sits in-plane offset by z*m, where m is the in-plane
! slope of V in element coordinates. Straining the mid-surface then drags the top of the fibre relative to the bottom by
! (du/dx) * m * z, which the covariant strain registers as transverse shear even though the fibre has not rotated relative to the
! facet at all. In element coordinates that spurious shear is exactly the in-plane strain contracted with the slope,
!
!    gamma_xz = eps_xx * m1 + eps_xy * m2 ,   gamma_yz = eps_xy * m1 + eps_yy * m2
!
! written in the sign convention of BS as it is returned above, and it is removed below.
!
! For a shell whose geometry is modelled exactly the director is the surface normal, m is zero and the term does not exist. It is
! produced purely by the mismatch between the flat facet and the nodal normals, so it is faceting error, not physics, and it grows
! linearly with the tilt. Subtracting it leaves the element free of membrane-shear coupling. Only the symmetric (strain) part of
! the in-plane displacement gradient is removed. The antisymmetric part is the in-plane rigid rotation, whose contribution is
! cancelled by the corresponding fibre rotation, so removing it as well would destroy rigid body invariance.

  ! CALL MITC_COVARIANT_BASIS( R, S, ZERO, G )

  ! IF (DABS(G(3,3)) > 1.0D-12 * DSQRT(DOT_PRODUCT(G(:,3), G(:,3)))) THEN

     ! M1 = G(1,3) / G(3,3)
     ! M2 = G(2,3) / G(3,3)

     ! BS(1,:) = BS(1,:) - ( BM(1,:) * M1 + HALF * BM(3,:) * M2 )
     ! BS(2,:) = BS(2,:) - ( HALF * BM(3,:) * M1 + BM(2,:) * M2 )

  ! ENDIF






  RETURN

! **********************************************************************************************************************************

CONTAINS

! **********************************************************************************************************************************

  SUBROUTINE TO_ELEMENT_BASIS ( T, B )

  ! Rotate one strain-displacement matrix from the cartesian local basis to the element coordinate system.

  REAL(DOUBLE), INTENT(IN)    :: T                 ! Isoparametric thickness coordinate the matrix was evaluated at
  REAL(DOUBLE), INTENT(INOUT) :: B(6, 6*ELGP)
  REAL(DOUBLE)                :: TRANSFORM(3,3)    ! Cartesian local basis vectors, in element coordinates

  TRANSFORM = MITC4_CARTESIAN_LOCAL_BASIS( R, S, T )

  B(4:6,:) = B(4:6,:) / 2                          ! Remove the engineering shear factor of 2 to rotate as a tensor
  CALL MITC_TRANSFORM_B( TRANSFORM, B )
  B(4:6,:) = B(4:6,:) * 2                          ! Reinstate it

  RETURN

  END SUBROUTINE TO_ELEMENT_BASIS

END SUBROUTINE MITC4_BMBS

      SUBROUTINE MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION ( R, S, T, X_R, X_S, X_D, MEMBRANE, BENDING, ROW_FROM, ROW_TO, B )

! Reference [1]:
! "A new MITC4+ shell element" by Ko, Lee, Bathe, 2016
!
! Reference [2]:
!  MITC4 paper "A continuum mechanics based four-node shell element for general nonlinear analysis"
!     by Dvorkin and Bathe


! Covariant strain-displacement components at point (R,S,T) directly evaluated from the displacement and rotation interpolations.
! Only for in-layer strains.
!
!        Grid point 1        Grid point 2      ...
!      ux uy uz rx ry rz   ux uy uz rx ry rz
! 11 [ #  #  #  #  #  #  | #  #  #  #  #  #  |     ]
! 22 [ #  #  #  #  #  #  | #  #  #  #  #  #  |     ]
! 33 [ 0  0  0  0  0  0  | 0  0  0  0  0  0  |     ]
! 12 [ #  #  #  #  #  #  | #  #  #  #  #  #  | ... ]
! 23 [ 0  0  0  0  0  0  | 0  0  0  0  0  0  |     ]
! 13 [ 0  0  0  0  0  0  | 0  0  0  0  0  0  |     ]


      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP, TYPE
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, ONE, TWO, FOUR, QUARTER
      USE SCONTR, ONLY                :  FATAL_ERR
      USE MITC_STUF, Only             :  DIRECTOR, DIR_THICKNESS, GP_RS

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE VECTOR_GEOMETRY, ONLY       :  CROSS

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: ROW_FROM          ! First row of B to generate. Strain component index 1-4.
      INTEGER(LONG), INTENT(IN)       :: ROW_TO            ! Last row of B to generate. Strain component index 1-4.
      INTEGER(LONG)                   :: GP                ! Grid point number. 1-4.
      INTEGER(LONG)                   :: I, J              ! Tensor indices.
      INTEGER(LONG)                   :: ROW               ! Row number of B
      INTEGER(LONG)                   :: K                 ! Column of B before the column for DOF 1 of the current node.

      REAL(DOUBLE) , INTENT(IN)       :: R,S,T             ! Isparametric coordinates
      REAL(DOUBLE) , INTENT(IN)       :: X_R(3)            ! Characteristic geometry vector x_r
      REAL(DOUBLE) , INTENT(IN)       :: X_S(3)            ! Characteristic geometry vector x_s
      REAL(DOUBLE) , INTENT(IN)       :: X_D(3)            ! Characteristic geometry vector x_d (distortion vector)
      REAL(DOUBLE) , INTENT(OUT)      :: B(6, 6*ELGP)      ! Strain-displacement matrix.
      REAL(DOUBLE)                    :: PSH(ELGP)         ! Shape functions
      REAL(DOUBLE)                    :: DPSHG(2,ELGP)     ! Derivatives of shape functions with respect to R and S.
      REAL(DOUBLE)                    :: DXMDRS(3,2)       ! Partial derivatives of x_m with respect to r and s
      REAL(DOUBLE)                    :: DXBDRS(3,2)       ! Partial derivatives of x_b with respect to r and s
      REAL(DOUBLE)                    :: V1(3,ELGP)        ! Basis vector orthogonal to the director vector.
      REAL(DOUBLE)                    :: V2(3,ELGP)        ! Basis vector orthogonal to the director vector and V1.
      REAL(DOUBLE)                    :: TRANSFORM(3,3)    ! Transformation matrix.

      LOGICAL      , INTENT(IN)       :: MEMBRANE          ! If true, generate membrane parts of B
      LOGICAL      , INTENT(IN)       :: BENDING           ! If true, generate bending parts of B

! **********************************************************************************************************************************

                                                           ! Shape function derivatives at R,S
      CALL MITC_SHAPE_FUNCTIONS(R, S, PSH, DPSHG)

                                                           ! Initialize B
      B(ROW_FROM:ROW_TO,:) = ZERO

                                                           ! Eqn (9) of ref [1].
      DXMDRS(:,1) = X_R + S * X_D                          ! ∂x_m/∂r
      DXMDRS(:,2) = X_S + R * X_D                          ! ∂x_m/∂s


      IF(BENDING) THEN

                                                           ! ∂x_b/∂r
                                                           ! ∂x_b/∂s
                                                           ! From eqns (8a) and (2) of ref [1].
         DXBDRS(:,:) = ZERO
         DO GP=1,ELGP
            DXBDRS(:,1) = DXBDRS(:,1) + HALF * DIR_THICKNESS(GP) * DIRECTOR(:,GP) * DPSHG(1,GP)
            DXBDRS(:,2) = DXBDRS(:,2) + HALF * DIR_THICKNESS(GP) * DIRECTOR(:,GP) * DPSHG(2,GP)
         ENDDO

                                                           ! Find a V1 and V2 for each node which form an
                                                           ! orthogonal right-handed coordinate system V1, V2, Vn
                                                           ! where Vn is the director vector.
         DO GP=1,ELGP
                                                           ! X_R is a convenient vector that's never parallel to Vn.
                                                           ! Project X_R onto the plane normal to the director vector.
            V1(:,GP) = X_R - DIRECTOR(:,GP) * DOT_PRODUCT(X_R, DIRECTOR(:,GP)) / DOT_PRODUCT(DIRECTOR(:,GP), DIRECTOR(:,GP))
                                                           ! Normalize V1
            V1(:,GP) = V1(:,GP) / DSQRT(DOT_PRODUCT(V1(:,GP), V1(:,GP)))
                                                           ! Calculate V2
            CALL CROSS(DIRECTOR(:,GP), V1(:,GP), V2(:,GP))
         ENDDO

      ENDIF

      DO ROW=ROW_FROM,ROW_TO

                                                           ! Tensor indices for the row
         SELECT CASE (ROW)
            CASE (1); I=1; J=1                             ! In-layer normal strain
            CASE (2); I=2; J=2                             ! In-layer normal strain
            CASE (3); I=0; J=0; CYCLE                      ! No zz strain
            CASE (4); I=1; J=2                             ! In-layer shear strain
            CASE DEFAULT
               I=0; J=0
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
         END SELECT

         IF(MEMBRANE) THEN
                                                           ! Membrane e^m_xx, e^m_yy, e^m_xy terms of eqn (7a)
                                                           ! described in eqn (7b) in ref [1]

                                                           !              1  / ∂x_m     ∂u_m \
                                                           ! B(ROW,:) +=  - (  ---- dot ----  )
                                                           !              2  \ ∂r_i     ∂r_j /
            CALL ADD_TERM_M(ROW, DXMDRS(:,I), J, ONE)
                                                           !              1  / ∂x_m     ∂u_m \
                                                           ! B(ROW,:) +=  - (  ---- dot ----  )
                                                           !              2  \ ∂r_j     ∂r_i /
            CALL ADD_TERM_M(ROW, DXMDRS(:,J), I, ONE)


         ENDIF

         IF(BENDING) THEN
                                                           ! Bending e^b1_xx, e^b1_yy, e^b1_xy terms of eqn (7a)
                                                           ! described in eqn (7c) in ref [1]

                                                           !              t  / ∂x_m     ∂u_b \
                                                           ! B(ROW,:) +=  - (  ---- dot ----  )
                                                           !              2  \ ∂r_i     ∂r_j /
            CALL ADD_TERM_B(ROW, DXMDRS(:,I), J, T)
                                                           !              t  / ∂x_m     ∂u_b \
                                                           ! B(ROW,:) +=  - (  ---- dot ----  )
                                                           !              2  \ ∂r_j     ∂r_i /
            CALL ADD_TERM_B(ROW, DXMDRS(:,J), I, T)


                                                           !              t  / ∂x_b     ∂u_m \
                                                           ! B(ROW,:) +=  - (  ---- dot ----  )
                                                           !              2  \ ∂r_i     ∂r_j /
            CALL ADD_TERM_M(ROW, DXBDRS(:,I), J, T)
                                                           !              t  / ∂x_b     ∂u_m \
                                                           ! B(ROW,:) +=  - (  ---- dot ----  )
                                                           !              2  \ ∂r_j     ∂r_i /
            CALL ADD_TERM_M(ROW, DXBDRS(:,J), I, T)

                                                           ! Bending e^b2_xx, e^b2_yy, e^b2_xy terms of eqn (7a)
                                                           ! described in eqn (7d) in ref [1]

                                                           !              t^2  / ∂x_b     ∂u_b \
                                                           ! B(ROW,:) +=  --- (  ---- dot ----  )
                                                           !               2   \ ∂r_i     ∂r_j /
            CALL ADD_TERM_B(ROW, DXBDRS(:,I), J, T*T)
                                                           !              t^2  / ∂x_b     ∂u_b \
                                                           ! B(ROW,:) +=  --- (  ---- dot ----  )
                                                           !               2   \ ∂r_j     ∂r_i /
            CALL ADD_TERM_B(ROW, DXBDRS(:,J), I, T*T)

         ENDIF

      ENDDO


      IF(BENDING) THEN
                                                           ! Transform the rotational dof terms of B from V1,V2,Vn
                                                           ! coordinates to basic x,y,z.
         DO GP=1,ELGP
            K = (GP-1) * 6
            TRANSFORM(:,1) = V1(:,GP)
            TRANSFORM(:,2) = V2(:,GP)
            TRANSFORM(:,3) = DIRECTOR(:,GP)
            DO ROW=ROW_FROM,ROW_TO
               IF(ROW /= 3) THEN
                  B(ROW,K+4:K+6) = MATMUL(TRANSFORM, B(ROW,K+4:K+6))
               ENDIF
            ENDDO
         ENDDO

      ENDIF


      RETURN

! **********************************************************************************************************************************

      CONTAINS

! **********************************************************************************************************************************

      SUBROUTINE ADD_TERM_M(ROW, LEFT, IU, COEFFICIENT)

      !              COEFFICIENT  /           ∂u_m  \
      ! B(ROW,:) +=  ----------- (  LEFT  dot -----  )
      !                   2       \           ∂r_IU /

      REAL(DOUBLE) , INTENT(IN)       :: LEFT(3)           ! The vector on the left of the dot product
      REAL(DOUBLE) , INTENT(IN)       :: COEFFICIENT       ! Scalar to multiply each term by. Coefficient in eqn (7a) of ref [1]
      REAL(DOUBLE)                    :: DUMDRS(3)         ! One grid point's term in the sum for the coefficients of the partial
                                                           ! derivatives of u_m with respect to R or S.

      INTEGER(LONG), INTENT(IN)       :: ROW               ! Row of B to add the result to. 1, 2, or 4.
      INTEGER(LONG), INTENT(IN)       :: IU                ! Index of dr in the u derivative. 1 or 2.

      DO GP=1,ELGP
         K = (GP-1) * 6
                                                           ! Eqn (9) of ref [1]. This node's term of:
         IF (IU == 1) THEN                                 ! ∂u_m/∂r = u_r + s * u_d
            DUMDRS = ( GP_RS(1,GP) + S * GP_RS(1,GP) * GP_RS(2,GP) ) / FOUR
         ELSE                                              ! IU == 2
            DUMDRS = ( GP_RS(2,GP) + R * GP_RS(1,GP) * GP_RS(2,GP) ) / FOUR
         ENDIF

         B(ROW, K+1) = B(ROW, K+1) + COEFFICIENT / TWO * LEFT(1) * DUMDRS(1)
         B(ROW, K+2) = B(ROW, K+2) + COEFFICIENT / TWO * LEFT(2) * DUMDRS(2)
         B(ROW, K+3) = B(ROW, K+3) + COEFFICIENT / TWO * LEFT(3) * DUMDRS(3)
      ENDDO

      END SUBROUTINE ADD_TERM_M

! **********************************************************************************************************************************

      SUBROUTINE ADD_TERM_B(ROW, LEFT, IU, COEFFICIENT)

      !              COEFFICIENT  /           ∂u_b  \
      ! B(ROW,:) +=  ----------- (  LEFT  dot -----  )
      !                   2       \           ∂r_IU /

      REAL(DOUBLE) , INTENT(IN)       :: LEFT(3)           ! The vector on the left of the dot product
      REAL(DOUBLE) , INTENT(IN)       :: COEFFICIENT       ! Scalar to multiply each term by. Coefficient in eqn (7a) of ref [1]
      REAL(DOUBLE)                    :: DUMa(3)           ! Coefficient of alpha in ∂u_b/∂r_IU for one node.
      REAL(DOUBLE)                    :: DUMb(3)           ! Coefficient of beta  in ∂u_b/∂r_IU for one node.

      INTEGER(LONG), INTENT(IN)       :: ROW               ! Row of B to add the result to. 1, 2, or 4.
      INTEGER(LONG), INTENT(IN)       :: IU                ! Index of dr in the u derivative. 1 or 2.

      DO GP=1,ELGP
         K = (GP-1) * 6

                                                           ! Put coefficients of alpha in DOF 4.
         DUMa = DIR_THICKNESS(GP) / TWO * DPSHG(IU,GP) * (-V2(:,GP))
         B(ROW, K+4) = B(ROW, K+4) + COEFFICIENT / TWO * DOT_PRODUCT(LEFT, DUMa)

                                                           ! Put coefficients of beta in DOF 5.
         DUMb = DIR_THICKNESS(GP) / TWO * DPSHG(IU,GP) * ( V1(:,GP))
         B(ROW, K+5) = B(ROW, K+5) + COEFFICIENT / TWO * DOT_PRODUCT(LEFT, DUMb)

      ENDDO

      END SUBROUTINE ADD_TERM_B

! **********************************************************************************************************************************


      END SUBROUTINE MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION

   END MODULE MITC4_MOD
