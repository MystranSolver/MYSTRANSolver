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

   MODULE MITC3_MOD

! MITC3+ 3-node triangular shell element (PARAM TRIA3TYP = MITC3+).
!
! Reference [1]:
!  "The MITC3+ shell element and its performance", Y. Lee, P.S. Lee, K.J. Bathe, Computers & Structures 138 (2014) 12-23.
! Reference [2]:
!  "A continuum mechanics based four-node shell element for general nonlinear analysis", E.N. Dvorkin, K.J. Bathe (1984).
!
! The element is the continuum-based (degenerated) shell of ref [2] on a 3-node triangle, with the rotation field enriched by a
! cubic bubble at an internal node (ref [1] eqns 6-8) and the transverse shear strains tied at the interior points A-F of ref [1]
! (eqn 17). The two bubble rotations are condensed out at element level, so the element presents the usual 18 DOF to the rest of
! MYSTRAN. The in-plane strains are displacement based (constant strain membrane), exactly as in ref [1].
!
! The through-thickness material integration is done the same way as for the MITC4/MITC4+ element (subr MITC4): the mid-surface
! membrane, curvature and transverse shear strain-displacement operators are extracted from the 3-D covariant strain and combined
! with the laminate matrices SHELL_A, SHELL_B, SHELL_D, SHELL_T (and SHELL_AALP, SHELL_BALP, SHELL_DALP for thermal loads) from
! subr SHELL_ABD_MATRICES. That gives PSHELL (including MID4 membrane-bending coupling), PCOMP (symmetric or not, with Z0) and
! thermal loading on all of them with one code path. Differential stiffness is built from per-ply membrane forces as in MITC4.
!
! DOF layout used internally while the matrices are built (NDOFB = 20 columns):
!    columns  1-18: the 3 grids, 6 DOF each, in element coordinates (as everywhere else in MYSTRAN)
!    columns 19-20: alpha4, beta4 - the bubble node rotations about its V1 and V2 vectors (ref [1] eqn 8)
! Everything handed back in KE, PTE, BEi is condensed to 18 columns.

   USE MITC_KERNELS, ONLY :  MITC_CONTRAVARIANT_BASIS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MITC3

   INTEGER, PARAMETER, PRIVATE     :: NDOFB = 20            ! 18 grid DOF + 2 bubble rotations

   CONTAINS
! ##################################################################################################################################

      SUBROUTINE MITC3 ( OPT, INT_ELEM_ID )

! Calculates, or calls subr's to calculate, MITC3+ triangular element matrices:

!  1) ME        = element mass matrix                  , if OPT(1) = 'Y'
!  2) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  3) BEi       = element strain data recovery matrices, if OPT(3) = 'Y'
!  4) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  5) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  6) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_TRIA, NTSUB, NSUB
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  NUM_EMG_FATAL_ERRS, PCOMP_PROPS, ELGP, KE, EM, BE1, BE2, BE3, PHI_SQ, EID, TYPE,          &
                                         EPROP, PTE, TREF, DT, PPE, PRESS, MASS_PER_UNIT_AREA, NUM_PLIES, PLY_NUM, TPLY, STRESS,   &
                                         KED, SHELL_A, SHELL_B, SHELL_D, SHELL_T, SHELL_AALP, SHELL_BALP, SHELL_DALP
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, ONE, TWO, THREE, THIRD
      USE MITC_STUF, ONLY             :  THERM_KAPPA, MSTRPT
      USE RESULT_COORDINATES, ONLY    :  ELMDIS, ELMDIS_PLY
      USE ELEMENT_RECOVERY_SUPPORT, ONLY:  ELEM_STRE_STRN_ARRAYS
      USE QUADRATURE, ONLY            :  ORDER_TRIA
      USE FILE_LIFECYCLE, ONLY        :  OUTA_HERE
      USE MASS_DOF_EXPANSION, ONLY    :  EXPAND_MASS_DOFS
      USE COMPOSITE_SHELL_PREPARATION, ONLY:  GET_ELEM_NUM_PLIES, SHELL_ABD_MATRICES

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MITC3'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID
      INTEGER(LONG), PARAMETER        :: IORD_IJ = 3       ! In-plane triangular integration order (ref [1] uses 3 points)
      INTEGER(LONG)                   :: I,J,K,L           ! DO loop indices
      INTEGER(LONG)                   :: GP                ! Element grid point number
      INTEGER(LONG)                   :: JPLY              ! PLY_NUM in a DO loop
      INTEGER(LONG)                   :: KI, KJ            ! For converting grid point number to element DOF number

      REAL(DOUBLE)                    :: SS_I(MAX_ORDER_TRIA) ! Triangular integration abscissa's (r)
      REAL(DOUBLE)                    :: SS_J(MAX_ORDER_TRIA) ! Triangular integration abscissa's (s)
      REAL(DOUBLE)                    :: HH_IJ(MAX_ORDER_TRIA)! Triangular integration weights (sum to 1/2)
      REAL(DOUBLE)                    :: R, S              ! Isoparametric coordinates of a point
      REAL(DOUBLE)                    :: DETJ              ! Mid-surface area Jacobian |g_r x g_s|
      REAL(DOUBLE)                    :: INTFAC            ! Integration factor DETJ*weight
      REAL(DOUBLE)                    :: AREA              ! Element area
      REAL(DOUBLE)                    :: G(3,3)            ! Covariant basis vectors
      REAL(DOUBLE)                    :: NORMAL(3)         ! Unit director at a point
      REAL(DOUBLE)                    :: H(3), DH(2,3)     ! Corner shape functions and derivatives
      REAL(DOUBLE)                    :: F(4), DF(2,4)     ! Enriched rotation shape functions and derivatives

      REAL(DOUBLE)                    :: BMI3(3,NDOFB)     ! Mid-surface membrane strain-disp operator
      REAL(DOUBLE)                    :: BBI3(3,NDOFB)     ! Mid-surface curvature strain-disp operator
      REAL(DOUBLE)                    :: BSI2(2,NDOFB)     ! Mid-surface trans shear strain-disp operator
      REAL(DOUBLE)                    :: KB(NDOFB,NDOFB)   ! Element stiffness incl bubble DOF, before condensation
      REAL(DOUBLE)                    :: KBB(2,2)          ! Bubble-bubble partition of KB
      REAL(DOUBLE)                    :: KBBI(2,2)         ! Inverse of KBB
      REAL(DOUBLE)                    :: CMAT(2,18)        ! KBB^-1 * KBA: bubble DOF per unit grid DOF
      REAL(DOUBLE)                    :: DETKBB            ! Determinant of KBB
      REAL(DOUBLE)                    :: TRKAA             ! Trace of grid DOF partition of KB (scale for the KBB check)
      LOGICAL                         :: CONDENSE          ! .TRUE. if the bubble DOF carry stiffness and are condensed

      REAL(DOUBLE)                    :: UNIT_PTE(NDOFB)   ! Thermal load vector for unit temperature change.
      REAL(DOUBLE)                    :: UNIT_PTE_G(NDOFB) ! Thermal load vector for unit through-thickness temperature gradient.
      REAL(DOUBLE)                    :: KAP_TH(3)         ! Free thermal curvature per unit temperature (curved shells w/ SNORM)
      REAL(DOUBLE)                    :: TBAR              ! Average elem temperature
      REAL(DOUBLE)                    :: UNIT_PPE(18)      ! Pressure load vector for unit pressure.
      REAL(DOUBLE)                    :: M_1DOF(ELGP,ELGP) ! Consistent mass matrix with 1 DOF per node.

      REAL(DOUBLE)                    :: FORCEx, FORCEy, FORCExy ! Membrane engineering forces at the centroid (summed over plies)
      REAL(DOUBLE)                    :: DUM11(2,2)        ! Membrane force tensor
      REAL(DOUBLE)                    :: DPSHX(2,3)        ! Derivatives of corner shape functions wrt elem x, y
      REAL(DOUBLE)                    :: JAC2x2(2,2)       ! Mid-surface 2x2 Jacobian
      REAL(DOUBLE)                    :: JACI(2,2)         ! Inverse of JAC2x2
      REAL(DOUBLE)                    :: KS(3,3)           ! KED for one DOF

! **********************************************************************************************************************************
! COORDINATE SYSTEMS
! ==================
!
! Element: x_element along side 1-2, z_element normal to the facet (subr ELMGM1 puts XEL in this system). Used for the grid point
!          DOFs of all element matrices and for stress, strain and force output. This is also MSC's CTRIA3 element system.
!
! Isoparametric: R, S in [0,1] on the parent triangle (h1 = 1-R-S, h2 = R, h3 = S), T in [-1,1] along the director.
!
! Covariant / contravariant: g_r, g_s, g_t and g^r, g^s, g^t of the degenerated shell, as in subr MITC4.
!
! Because the triangle is flat and the director is normally the facet normal, the covariant strain components are carried
! straight to the element Cartesian system with the contravariant basis (eps_kl = e_ij g^i_k g^j_l); there is no separate
! "cartesian local" system as in MITC4. With SNORM directors that are not normal to the facet the same transformation applies.
! **********************************************************************************************************************************

! Initialize

      PHI_SQ  = ONE                                        ! No shear correction factor for this element (MIN3 only)

      CALL MITC3_INITIALIZE

      CALL ORDER_TRIA ( IORD_IJ, SS_I, SS_J, HH_IJ )

      CALL MITC3_COVARIANT_BASIS ( THIRD, THIRD, ZERO, G )
      DETJ = DSQRT( DOT_PRODUCT( CROSS3(G(:,1),G(:,2)), CROSS3(G(:,1),G(:,2)) ) )
      AREA = HALF*DETJ

! **********************************************************************************************************************************
! Generate the mass matrix for this element (translational, at the grid plane, no rotary inertia - same as MITC4).

      IF (OPT(1) == 'Y') THEN

         M_1DOF(:,:) = ZERO

         DO I=1,IORD_IJ
            R = SS_I(I)
            S = SS_J(I)
            CALL MITC3_SHAPE_FUNCTIONS ( R, S, H, DH, F, DF )
            CALL MITC3_COVARIANT_BASIS ( R, S, ZERO, G )
            DETJ = DSQRT( DOT_PRODUCT( CROSS3(G(:,1),G(:,2)), CROSS3(G(:,1),G(:,2)) ) )
            INTFAC = DETJ*HH_IJ(I)
            DO K=1,ELGP
               DO L=1,ELGP
                  M_1DOF(K,L) = M_1DOF(K,L) + H(K)*H(L)*MASS_PER_UNIT_AREA*INTFAC
               ENDDO
            ENDDO
         ENDDO

         CALL EXPAND_MASS_DOFS ( M_1DOF )

      ENDIF

! **********************************************************************************************************************************
! Element stiffness matrix with the bubble DOF, and its static condensation. Needed for KE (OPT(4)), and for the condensation of the
! thermal loads (OPT(2)) and of the strain recovery matrices (OPT(3), OPT(6)).
!
!    KE = int over mid-surface area of:
!       [Bm]' [A] [Bm]  +  [Bm]' [B] [Bb] + [Bb]' [B]' [Bm]  +  [Bb]' [D] [Bb]  +  [Bs]' [T] [Bs]  dA
!
! with Bm, Bb, Bs from MITC3_BMBS. SHELL_B is used whenever it is nonzero (unsymmetric PCOMP, PCOMP LAM=SYM with Z0, PSHELL MID4),
! there is no PCOMP_LAM gate as in TREL1.

      CONDENSE  = .FALSE.
      CMAT(:,:) = ZERO
      KB(:,:)   = ZERO

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         DO I=1,IORD_IJ
            R = SS_I(I)
            S = SS_J(I)

            CALL MITC3_BMBS ( R, S, BMI3, BBI3, BSI2 )

            CALL MITC3_COVARIANT_BASIS ( R, S, ZERO, G )
            DETJ = DSQRT( DOT_PRODUCT( CROSS3(G(:,1),G(:,2)), CROSS3(G(:,1),G(:,2)) ) )
            INTFAC = DETJ*HH_IJ(I)

            KB = KB + ( MATMUL( TRANSPOSE(BMI3), MATMUL( SHELL_A, BMI3 ) )                                                           &
                      + MATMUL( TRANSPOSE(BBI3), MATMUL( SHELL_D, BBI3 ) )                                                           &
                      + MATMUL( TRANSPOSE(BSI2), MATMUL( SHELL_T, BSI2 ) )                                                           &
                      + MATMUL( TRANSPOSE(BMI3), MATMUL( SHELL_B, BBI3 ) )                                                           &
                      + MATMUL( TRANSPOSE(BBI3), MATMUL( TRANSPOSE(SHELL_B), BMI3 ) ) ) * INTFAC
         ENDDO

                                                           ! Static condensation of the 2 bubble rotations:
                                                           ! KE = Kaa - Kab Kbb^-1 Kba. The bubble only carries bending and transverse
                                                           ! shear stiffness, so Kbb is singular for a membrane-only element (no MID2):
                                                           ! then Kab is zero as well and there is nothing to condense.
         KBB(1:2,1:2) = KB(19:20,19:20)
         DETKBB = KBB(1,1)*KBB(2,2) - KBB(1,2)*KBB(2,1)
         TRKAA  = ZERO
         DO I=1,18
            TRKAA = TRKAA + DABS(KB(I,I))
         ENDDO

         IF (DABS(DETKBB) > 1.0D-24*TRKAA*TRKAA) THEN
            CONDENSE  = .TRUE.
            KBBI(1,1) =  KBB(2,2)/DETKBB
            KBBI(1,2) = -KBB(1,2)/DETKBB
            KBBI(2,1) = -KBB(2,1)/DETKBB
            KBBI(2,2) =  KBB(1,1)/DETKBB
            CMAT = MATMUL( KBBI, KB(19:20,1:18) )
         ENDIF

         IF (OPT(4) == 'Y') THEN
            KE(1:18,1:18) = KB(1:18,1:18)
            IF (CONDENSE) THEN
               KE(1:18,1:18) = KE(1:18,1:18) - MATMUL( KB(1:18,19:20), CMAT )
            ENDIF
                                                           ! Symmetrize (round-off only)
            DO I=2,18
               DO J=1,I-1
                  KE(I,J) = HALF*(KE(I,J) + KE(J,I))
                  KE(J,I) = KE(I,J)
               ENDDO
            ENDDO
         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Calculate element thermal loads. Same assembly as subr MITC4:
!
!    PTE = ( [Bm]' ( [A_alpha] + [B] {k_th} ) + [Bb]' ( [B_alpha] + [D] {k_th} ) ) * (TBAR - TREF) + [Bb]' [D_alpha] * TPRIME
!
! integrated over the mid-surface area. SHELL_BALP is the thermal moment of a coupled laminate (nonsymmetric PCOMP, PCOMP with Z0,
! PSHELL MID4); {k_th} is the free thermal curvature of a curved reference surface (nonzero only with SNORM directors that are not
! normal to the facet, see MITC3_THERM_CURV). TPRIME is the TEMPP1 through-thickness gradient, stored by ELMDAT2 in DT(ELGP+1,J).
! The bubble entries of the load vector are condensed with the stiffness (they vanish for a uniform load: the bubble curvature
! operator integrates to zero over the element).

      IF (OPT(2) == 'Y') THEN

         UNIT_PTE(:)   = ZERO
         UNIT_PTE_G(:) = ZERO

         DO I=1,IORD_IJ
            R = SS_I(I)
            S = SS_J(I)

            CALL MITC3_BMBS ( R, S, BMI3, BBI3, BSI2 )

            CALL MITC3_COVARIANT_BASIS ( R, S, ZERO, G )
            DETJ = DSQRT( DOT_PRODUCT( CROSS3(G(:,1),G(:,2)), CROSS3(G(:,1),G(:,2)) ) )
            INTFAC = DETJ*HH_IJ(I)

            CALL MITC3_THERM_CURV ( R, S, KAP_TH )

            UNIT_PTE   = UNIT_PTE   + ( MATMUL( TRANSPOSE(BMI3), SHELL_AALP + MATMUL( SHELL_B, KAP_TH ) )                          &
                                      + MATMUL( TRANSPOSE(BBI3), SHELL_BALP + MATMUL( SHELL_D, KAP_TH ) ) ) * INTFAC
            UNIT_PTE_G = UNIT_PTE_G +   MATMUL( TRANSPOSE(BBI3), SHELL_DALP ) * INTFAC
         ENDDO

         IF (CONDENSE) THEN
            UNIT_PTE(1:18)   = UNIT_PTE(1:18)   - MATMUL( KB(1:18,19:20), MATMUL( KBBI, UNIT_PTE(19:20)   ) )
            UNIT_PTE_G(1:18) = UNIT_PTE_G(1:18) - MATMUL( KB(1:18,19:20), MATMUL( KBBI, UNIT_PTE_G(19:20) ) )
         ENDIF

         DO J=1,NTSUB
            TBAR = (DT(1,J) + DT(2,J) + DT(3,J))/THREE
            PTE(1:18,J) = UNIT_PTE(1:18) * (TBAR - TREF(1)) + UNIT_PTE_G(1:18) * DT(ELGP+1,J)
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Strain recovery matrices at the element centroid (the TRIA3 has 1 stress point). The bubble DOF are condensed into the grid DOF:
! B_recovered = B_a - B_b Kbb^-1 Kba (the bubble carries no load in linear statics, see the thermal load above). The bubble has no
! curvature at the centroid, so only the transverse shear operator is actually changed by the condensation.

      IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         CALL MITC3_BMBS ( THIRD, THIRD, BMI3, BBI3, BSI2 )

         BE1(1:3,1:18,1) = BMI3(1:3,1:18)
         BE2(1:3,1:18,1) = BBI3(1:3,1:18)
         BE3(1:2,1:18,1) = BSI2(1:2,1:18)
         IF (CONDENSE) THEN
            BE1(1:3,1:18,1) = BE1(1:3,1:18,1) - MATMUL( BMI3(1:3,19:20), CMAT )
            BE2(1:3,1:18,1) = BE2(1:3,1:18,1) - MATMUL( BBI3(1:3,19:20), CMAT )
            BE3(1:2,1:18,1) = BE3(1:2,1:18,1) - MATMUL( BSI2(1:2,19:20), CMAT )
         ENDIF

         THERM_KAPPA(:,:) = ZERO
         CALL MITC3_THERM_CURV ( THIRD, THIRD, THERM_KAPPA(1:3,1) )

      ENDIF

! **********************************************************************************************************************************
! Determine element pressure loads: work equivalent nodal forces of the pressure along the (interpolated) director, plus the in-plane
! components PRESS(1), PRESS(2) along element x, y for parity with TMEM1.

      IF (OPT(5) == 'Y') THEN

         UNIT_PPE(:) = ZERO
         PPE(1:18,1:NSUB) = ZERO

         DO I=1,IORD_IJ
            R = SS_I(I)
            S = SS_J(I)

            CALL MITC3_SHAPE_FUNCTIONS ( R, S, H, DH, F, DF )
            CALL MITC3_COVARIANT_BASIS ( R, S, ZERO, G )
            NORMAL(:) = G(:,3) / DSQRT(DOT_PRODUCT(G(:,3), G(:,3)))
            DETJ = DSQRT( DOT_PRODUCT( CROSS3(G(:,1),G(:,2)), CROSS3(G(:,1),G(:,2)) ) )
            INTFAC = DETJ*HH_IJ(I)

            DO GP=1,ELGP
               K = (GP-1)*6
               UNIT_PPE(K+1:K+3) = UNIT_PPE(K+1:K+3) + NORMAL*H(GP)*INTFAC
            ENDDO
         ENDDO

         DO J=1,NSUB
            PPE(1:18,J) = UNIT_PPE(1:18)*PRESS(3,J)
            DO GP=1,ELGP
               K = (GP-1)*6
               PPE(K+1,J) = PPE(K+1,J) + AREA*PRESS(1,J)/THREE
               PPE(K+2,J) = PPE(K+2,J) + AREA*PRESS(2,J)/THREE
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate linear differential stiffness matrix from the membrane forces at the centroid, summed over the plies from the ply
! strains e0 + z*k (subr ELMDIS_PLY), so that N = [A]{e0} + [B]{k} exactly for coupled laminates (same as MITC4). The membrane strain
! is constant over the element and the curvature varies linearly, so the centroid value is the exact element average.

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

         FORCEx  = ZERO
         FORCEy  = ZERO
         FORCExy = ZERO

         CALL GET_ELEM_NUM_PLIES ( INT_ELEM_ID )           ! Get NUM_PLIES

         DO JPLY=1,NUM_PLIES

            IF (PCOMP_PROPS == 'N') THEN                   ! EM is already set
               TPLY = EPROP(1)
               CALL ELMDIS
            ELSE
               PLY_NUM = JPLY                              ! Used by SHELL_ABD_MATRICES
               CALL SHELL_ABD_MATRICES ( INT_ELEM_ID, 'N' )! Get EM, ZPLY, TPLY, ALPVEC for this ply
               CALL ELMDIS
               CALL ELMDIS_PLY                             ! Adjust UEL using ZPLY
            ENDIF

            CALL ELEM_STRE_STRN_ARRAYS ( 1 )               ! Stress at the centroid using UEL, BE1, EM, ALPVEC, DT

            FORCEx  = FORCEx  + TPLY*STRESS(1)
            FORCEy  = FORCEy  + TPLY*STRESS(2)
            FORCExy = FORCExy + TPLY*STRESS(3)

         ENDDO

! According to:
!   Robert D. Cook, David S. Malkus, Michael E. Plesha Concepts and Applications of Finite Element Analysis, 3rd Edition  1989
!   Section 14.3 Stress Stiffness Matrix Of A Plate Element
!
! [k_sigma] = int  DPSHX' [ Nx Nxy ; Nxy Ny ] DPSHX dA,  DPSHX = J^-1 dN/d(r,s): derivatives of the corner shape functions wrt x, y.
! Everything is constant over the flat triangle, so the integral is AREA times the integrand.

         DUM11(1,1) = FORCEx  ; DUM11(1,2) = FORCExy
         DUM11(2,1) = FORCExy ; DUM11(2,2) = FORCEy

         CALL MITC3_SHAPE_FUNCTIONS ( THIRD, THIRD, H, DH, F, DF )
         CALL MITC3_COVARIANT_BASIS ( THIRD, THIRD, ZERO, G )
         JAC2x2(1,1) = G(1,1) ; JAC2x2(1,2) = G(2,1)       ! Row i = d(x,y)/d(xi_i) in element coords (the facet is the x-y plane)
         JAC2x2(2,1) = G(1,2) ; JAC2x2(2,2) = G(2,2)
         DETJ = JAC2x2(1,1)*JAC2x2(2,2) - JAC2x2(1,2)*JAC2x2(2,1)
         JACI(1,1) =  JAC2x2(2,2)/DETJ
         JACI(1,2) = -JAC2x2(1,2)/DETJ
         JACI(2,1) = -JAC2x2(2,1)/DETJ
         JACI(2,2) =  JAC2x2(1,1)/DETJ
         DPSHX = MATMUL( JACI, DH )

         KS = MATMUL( TRANSPOSE(DPSHX), MATMUL( DUM11, DPSHX ) ) * AREA

         KED(1:18,1:18) = ZERO
         DO I=1,ELGP
            DO J=1,ELGP
               KI = (I-1)*6
               KJ = (J-1)*6
               KED(KI+1,KJ+1) = KS(I,J)
               KED(KI+2,KJ+2) = KS(I,J)
               KED(KI+3,KJ+3) = KS(I,J)
            ENDDO
         ENDDO

      ENDIF

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MITC3

! ##################################################################################################################################

      FUNCTION CROSS3 ( A, B )

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: A(3), B(3)
      REAL(DOUBLE)                    :: CROSS3(3)

      CROSS3(1) = A(2)*B(3) - A(3)*B(2)
      CROSS3(2) = A(3)*B(1) - A(1)*B(3)
      CROSS3(3) = A(1)*B(2) - A(2)*B(1)

      END FUNCTION CROSS3

! ##################################################################################################################################

      SUBROUTINE MITC3_INITIALIZE

! Director vector and thickness at the 3 grids and at the internal bubble node, stored in MITC_STUF (slots 1-3 grids, slot 4 bubble):
!
!   DIRECTOR(:,i), DIR_THICKNESS(i) = V_n^i, a_i          of ref [1]
!   a_4 V_n^4 = (a_1 V_n^1 + a_2 V_n^2 + a_3 V_n^3)/3      ref [1] eqn 6
!
! The director is the facet normal (0,0,1 in element coordinates) unless the grid has an SNORM, in which case the SNORM is used in
! the same way as for the MITC4 QUAD4 (subr MITC_INITIALIZE) and the thickness along the director is scaled accordingly.
! V1 and V2 at the bubble node (BUB_V1, BUB_V2), needed for its 2 rotation DOF, are built from the side 1-2 direction.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR
      USE MODEL_STUF, ONLY            :  ELGP, EPROP, XEL, BGRID, GRID_SNORM, TE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, THIRD
      USE MITC_STUF, ONLY             :  DIRECTOR, DIR_THICKNESS, GP_RS, BUB_V1, BUB_V2
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF
      USE FILE_LIFECYCLE, ONLY        :  OUTA_HERE

      IMPLICIT NONE

      INTEGER(LONG)                   :: GP                ! Element grid point number

      REAL(DOUBLE)                    :: NORMAL(3)         ! Facet normal in element coordinates
      REAL(DOUBLE)                    :: SNORM(3)          ! SNORM of a grid in basic coordinates
      REAL(DOUBLE)                    :: THICKNESS_FACTOR
      REAL(DOUBLE)                    :: AV4(3)            ! a_4 V_n^4
      REAL(DOUBLE)                    :: XR(3)             ! Side 1-2 vector

! **********************************************************************************************************************************

      GP_RS(1,1) = ZERO  ;  GP_RS(2,1) = ZERO
      GP_RS(1,2) = ONE   ;  GP_RS(2,2) = ZERO
      GP_RS(1,3) = ZERO  ;  GP_RS(2,3) = ONE

      NORMAL = CROSS3( XEL(2,:) - XEL(1,:), XEL(3,:) - XEL(1,:) )
      NORMAL = NORMAL / DSQRT(DOT_PRODUCT(NORMAL, NORMAL))

      DO GP=1,ELGP

         SNORM = GRID_SNORM(BGRID(GP),:)

         IF (ANY(SNORM /= ZERO)) THEN                      ! Transform SNORM from basic to XEL element coordinates.
            CALL MATMULT_FFF ( TE, SNORM, 3, 3, 1, DIRECTOR(:,GP) )
            DIRECTOR(:,GP) = DIRECTOR(:,GP) / DSQRT(DOT_PRODUCT(DIRECTOR(:,GP), DIRECTOR(:,GP)))
         ELSE
            DIRECTOR(:,GP) = NORMAL
         ENDIF

         THICKNESS_FACTOR = DOT_PRODUCT(DIRECTOR(:,GP), NORMAL)

         IF (THICKNESS_FACTOR < 0.01D0) THEN
            WRITE(ERR,*) ' *ERROR: SNORM IS TOO FAR FROM MIDSURFACE NORMAL'
            WRITE(F06,*) ' *ERROR: SNORM IS TOO FAR FROM MIDSURFACE NORMAL'
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         DIR_THICKNESS(GP) = EPROP(1) / THICKNESS_FACTOR

      ENDDO

      AV4 = THIRD*( DIR_THICKNESS(1)*DIRECTOR(:,1) + DIR_THICKNESS(2)*DIRECTOR(:,2) + DIR_THICKNESS(3)*DIRECTOR(:,3) )
      DIR_THICKNESS(4) = DSQRT(DOT_PRODUCT(AV4, AV4))
      DIRECTOR(:,4)    = AV4 / DIR_THICKNESS(4)

      XR = XEL(2,:) - XEL(1,:)
      BUB_V1 = XR - DIRECTOR(:,4)*DOT_PRODUCT(XR, DIRECTOR(:,4))
      BUB_V1 = BUB_V1 / DSQRT(DOT_PRODUCT(BUB_V1, BUB_V1))
      BUB_V2 = CROSS3( DIRECTOR(:,4), BUB_V1 )

      RETURN

      END SUBROUTINE MITC3_INITIALIZE

! ##################################################################################################################################

      SUBROUTINE MITC3_SHAPE_FUNCTIONS ( R, S, H, DH, F, DF )

! Corner shape functions h_i and the enriched rotation shape functions f_i of ref [1] eqn 7:
!    h1 = 1-r-s, h2 = r, h3 = s
!    f_i = h_i - f4/3 (i=1,3),  f4 = 27 r s (1-r-s)   (cubic bubble, zero on all edges, 1 at the centroid)

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THIRD

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: R, S
      REAL(DOUBLE), INTENT(OUT)       :: H(3)              ! h_i
      REAL(DOUBLE), INTENT(OUT)       :: DH(2,3)           ! d h_i / d(r,s)
      REAL(DOUBLE), INTENT(OUT)       :: F(4)              ! f_i
      REAL(DOUBLE), INTENT(OUT)       :: DF(2,4)           ! d f_i / d(r,s)

      REAL(DOUBLE), PARAMETER         :: C27 = 27.0D0
      REAL(DOUBLE)                    :: F4, DF4(2)

! **********************************************************************************************************************************

      H(1) = ONE - R - S  ;  DH(1,1) = -ONE  ;  DH(2,1) = -ONE
      H(2) = R            ;  DH(1,2) =  ONE  ;  DH(2,2) = ZERO
      H(3) = S            ;  DH(1,3) = ZERO  ;  DH(2,3) =  ONE

      F4     = C27*R*S*(ONE - R - S)
      DF4(1) = C27*S*(ONE - TWO*R - S)
      DF4(2) = C27*R*(ONE - R - TWO*S)

      F(1:3)    = H(1:3)    - THIRD*F4
      DF(1,1:3) = DH(1,1:3) - THIRD*DF4(1)
      DF(2,1:3) = DH(2,1:3) - THIRD*DF4(2)
      F(4)      = F4
      DF(:,4)   = DF4

      RETURN

      END SUBROUTINE MITC3_SHAPE_FUNCTIONS

! ##################################################################################################################################

      SUBROUTINE MITC3_COVARIANT_BASIS ( R, S, T, G )

! Covariant basis g_r, g_s, g_t (columns of G) in element coordinates, ref [1] eqn 6:
!    x(r,s,t) = sum_i h_i x_i + (t/2) sum_{i=1}^{4} f_i a_i V_n^i

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  XEL
      USE CONSTANTS_1, ONLY           :  ZERO, HALF
      USE MITC_STUF, ONLY             :  DIRECTOR, DIR_THICKNESS

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: R, S, T
      REAL(DOUBLE), INTENT(OUT)       :: G(3,3)

      INTEGER(LONG)                   :: I
      REAL(DOUBLE)                    :: H(3), DH(2,3), F(4), DF(2,4)
      REAL(DOUBLE)                    :: AV(3)

! **********************************************************************************************************************************

      CALL MITC3_SHAPE_FUNCTIONS ( R, S, H, DH, F, DF )

      G(:,:) = ZERO
      DO I=1,3
         G(:,1) = G(:,1) + DH(1,I)*XEL(I,:)
         G(:,2) = G(:,2) + DH(2,I)*XEL(I,:)
      ENDDO
      DO I=1,4
         AV = DIR_THICKNESS(I)*DIRECTOR(:,I)
         G(:,1) = G(:,1) + HALF*T*DF(1,I)*AV
         G(:,2) = G(:,2) + HALF*T*DF(2,I)*AV
         G(:,3) = G(:,3) + HALF*F(I)*AV
      ENDDO

      RETURN

      END SUBROUTINE MITC3_COVARIANT_BASIS

! ##################################################################################################################################

      SUBROUTINE MITC3_DISP_GRAD ( R, S, T, UM, UB )

! Displacement gradient operators, ref [1] eqn 8:
!    u(r,s,t) = sum_i h_i u_i + (t/2) sum_{i=1}^{4} f_i a_i (theta_i x V_n^i)
! so that  d u / d xi_J = ( UM(:,:,J) + UB(:,:,J) ) * q,  J = 1 (r), 2 (s), 3 (t), with q the NDOFB element DOF.
! UM holds the grid translation terms (the "membrane" displacement u_m of MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION), UB the
! rotation terms (u_b). For the grids the rotation is the rotation vector theta in element coordinates (DOF 4-6), and
! theta x (a V_n) = -[a V_n]_x theta. For the bubble node theta_4 = alpha_4 V1 + beta_4 V2, so theta_4 x V_n^4 = -alpha_4 V2 + beta_4 V1.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, HALF
      USE MITC_STUF, ONLY             :  DIRECTOR, DIR_THICKNESS, BUB_V1, BUB_V2

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: R, S, T
      REAL(DOUBLE), INTENT(OUT)       :: UM(3,NDOFB,3)
      REAL(DOUBLE), INTENT(OUT)       :: UB(3,NDOFB,3)

      INTEGER(LONG)                   :: I, J, K
      REAL(DOUBLE)                    :: H(3), DH(2,3), F(4), DF(2,4)
      REAL(DOUBLE)                    :: AV(3)
      REAL(DOUBLE)                    :: CX(3,3)           ! theta x AV = CX * theta
      REAL(DOUBLE)                    :: CB(3,2)           ! Bubble: (alpha,beta) -> theta_4 x (a_4 V_n^4)

! **********************************************************************************************************************************

      CALL MITC3_SHAPE_FUNCTIONS ( R, S, H, DH, F, DF )

      UM(:,:,:) = ZERO
      UB(:,:,:) = ZERO

      DO I=1,3
         K = 6*(I-1)
         DO J=1,2
            UM(1,K+1,J) = DH(J,I)
            UM(2,K+2,J) = DH(J,I)
            UM(3,K+3,J) = DH(J,I)
         ENDDO
      ENDDO

      DO I=1,3
         K  = 6*(I-1)
         AV = DIR_THICKNESS(I)*DIRECTOR(:,I)
         CX(1,1) =  ZERO   ;  CX(1,2) =  AV(3)  ;  CX(1,3) = -AV(2)
         CX(2,1) = -AV(3)  ;  CX(2,2) =  ZERO   ;  CX(2,3) =  AV(1)
         CX(3,1) =  AV(2)  ;  CX(3,2) = -AV(1)  ;  CX(3,3) =  ZERO
         DO J=1,2
            UB(:,K+4:K+6,J) = UB(:,K+4:K+6,J) + HALF*T*DF(J,I)*CX
         ENDDO
         UB(:,K+4:K+6,3) = UB(:,K+4:K+6,3) + HALF*F(I)*CX
      ENDDO

      CB(:,1) = -DIR_THICKNESS(4)*BUB_V2
      CB(:,2) =  DIR_THICKNESS(4)*BUB_V1
      DO J=1,2
         UB(:,19:20,J) = UB(:,19:20,J) + HALF*T*DF(J,4)*CB
      ENDDO
      UB(:,19:20,3) = UB(:,19:20,3) + HALF*F(4)*CB

      RETURN

      END SUBROUTINE MITC3_DISP_GRAD

! ##################################################################################################################################

      SUBROUTINE MITC3_COVARIANT_STRAIN_ROW ( R, S, T, I, J, MEMBRANE, BENDING, ROW )

! One covariant strain component e_ij = 1/2 ( g_i . u_,j + g_j . u_,i ) at (R,S,T), directly from the interpolations, as a row of the
! strain-displacement matrix (NDOFB columns). I, J = 1 (r), 2 (s), 3 (t).
!
! The geometry is split as x = x_m + x_b into the mid-surface part (g_r, g_s at T = 0) and the director part (everything that
! carries T, plus g_t), and the displacement as u = u_m + u_b (translations / rotations), so that the membrane strain
! (x_m . u_m) and the bending strain (x_m . u_b, x_b . u_m, x_b . u_b) can be returned separately, exactly as in subr
! MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, HALF

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: R, S, T
      INTEGER(LONG), INTENT(IN)       :: I, J
      LOGICAL, INTENT(IN)             :: MEMBRANE, BENDING
      REAL(DOUBLE), INTENT(OUT)       :: ROW(NDOFB)

      REAL(DOUBLE)                    :: G(3,3), G0(3,3), GB(3,3)
      REAL(DOUBLE)                    :: UM(3,NDOFB,3), UB(3,NDOFB,3)

! **********************************************************************************************************************************

      CALL MITC3_COVARIANT_BASIS ( R, S, ZERO, G0 )        ! Mid-surface tangents (g_t set to zero below)
      CALL MITC3_COVARIANT_BASIS ( R, S, T   , G  )
      GB(:,1:2) = G(:,1:2) - G0(:,1:2)                     ! T-dependent part of g_r, g_s
      GB(:,3)   = G(:,3)                                   ! g_t is director only
      G0(:,3)   = ZERO

      CALL MITC3_DISP_GRAD ( R, S, T, UM, UB )

      ROW(:) = ZERO

      IF (MEMBRANE) THEN
         ROW = ROW + HALF*( MATMUL( G0(:,I), UM(:,:,J) ) + MATMUL( G0(:,J), UM(:,:,I) ) )
      ENDIF

      IF (BENDING) THEN
         ROW = ROW + HALF*( MATMUL( G0(:,I), UB(:,:,J) ) + MATMUL( G0(:,J), UB(:,:,I) ) )                                           &
                   + HALF*( MATMUL( GB(:,I), UM(:,:,J) ) + MATMUL( GB(:,J), UM(:,:,I) ) )                                           &
                   + HALF*( MATMUL( GB(:,I), UB(:,:,J) ) + MATMUL( GB(:,J), UB(:,:,I) ) )
      ENDIF

      RETURN

      END SUBROUTINE MITC3_COVARIANT_STRAIN_ROW

! ##################################################################################################################################

      SUBROUTINE MITC3_B ( R, S, T, MEMBRANE, BENDING, SHEAR, B )

! Strain-displacement matrix (6 x NDOFB) in the element Cartesian system at one point (R,S,T). Rows: xx, yy, zz, xy, yz, zx with the
! shear rows as engineering strains (same row convention as MITC4_B).
!
! In-plane covariant strains e_rr, e_ss, e_rs are displacement based. The transverse shear strains are the MITC3+ assumed strains of
! ref [1] eqn 17, tied at the points of ref [1] Table 1 with the tying distance d = 1/10000:
!
!    e_rt = 2/3 ( e_rt^B - 1/2 e_st^B ) + 1/3 ( e_rt^C + e_st^C ) + 1/3 c^ (3s - 1)
!    e_st = 2/3 ( e_st^A - 1/2 e_rt^A ) + 1/3 ( e_rt^C + e_st^C ) + 1/3 c^ (1 - 3r)
!    c^   = e_rt^F - e_rt^D - e_st^F + e_st^E
!
!      A = (1/6, 2/3)   B = (2/3, 1/6)   C = (1/6, 1/6)
!      D = (1/3+d, 1/3-2d)   E = (1/3-2d, 1/3+d)   F = (1/3+d, 1/3+d)
!
! The c^ term is O(d) by design: it gives the in-plane twisting mode a small but nonzero stiffness (ref [1] section 3).
! The covariant components are then carried to the element system with the contravariant basis, eps_kl = e_ij (g^i)_k (g^j)_l.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, HALF, THIRD, SIXTH

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: R, S, T
      LOGICAL, INTENT(IN)             :: MEMBRANE, BENDING, SHEAR
      REAL(DOUBLE), INTENT(OUT)       :: B(6,NDOFB)

      REAL(DOUBLE), PARAMETER         :: DTIE = 1.0D-4     ! Tying distance d of ref [1]
      REAL(DOUBLE), PARAMETER         :: TWO3 = TWO/THREE
      REAL(DOUBLE)                    :: BCOV(6,NDOFB)     ! Covariant (tensor) strain components: rr, ss, tt, rs, st, rt
      REAL(DOUBLE)                    :: ERT_A(NDOFB), EST_A(NDOFB), ERT_B(NDOFB), EST_B(NDOFB), ERT_C(NDOFB), EST_C(NDOFB)
      REAL(DOUBLE)                    :: ERT_D(NDOFB), ERT_F(NDOFB), EST_E(NDOFB), EST_F(NDOFB), CHAT(NDOFB)
      REAL(DOUBLE)                    :: G(3,3), GCON(3,3)

      INTEGER(LONG)                   :: ROW, K, L, I, J
      INTEGER(LONG), PARAMETER        :: INDEX1(6) = (/ 1, 2, 3, 1, 2, 1 /)
      INTEGER(LONG), PARAMETER        :: INDEX2(6) = (/ 1, 2, 3, 2, 3, 3 /)
      INTEGER(LONG), PARAMETER        :: ROWS(3,3) = RESHAPE( (/ 1, 4, 6, 4, 2, 5, 6, 5, 3 /), (/3,3/) )

! **********************************************************************************************************************************

      BCOV(:,:) = ZERO

      IF (MEMBRANE .OR. BENDING) THEN
         CALL MITC3_COVARIANT_STRAIN_ROW ( R, S, T, 1, 1, MEMBRANE, BENDING, BCOV(1,:) )
         CALL MITC3_COVARIANT_STRAIN_ROW ( R, S, T, 2, 2, MEMBRANE, BENDING, BCOV(2,:) )
         CALL MITC3_COVARIANT_STRAIN_ROW ( R, S, T, 1, 2, MEMBRANE, BENDING, BCOV(4,:) )
      ENDIF

      IF (SHEAR) THEN

         CALL MITC3_COVARIANT_STRAIN_ROW ( SIXTH         , TWO3          , T, 1, 3, .TRUE., .TRUE., ERT_A )
         CALL MITC3_COVARIANT_STRAIN_ROW ( SIXTH         , TWO3          , T, 2, 3, .TRUE., .TRUE., EST_A )
         CALL MITC3_COVARIANT_STRAIN_ROW ( TWO3          , SIXTH         , T, 1, 3, .TRUE., .TRUE., ERT_B )
         CALL MITC3_COVARIANT_STRAIN_ROW ( TWO3          , SIXTH         , T, 2, 3, .TRUE., .TRUE., EST_B )
         CALL MITC3_COVARIANT_STRAIN_ROW ( SIXTH         , SIXTH         , T, 1, 3, .TRUE., .TRUE., ERT_C )
         CALL MITC3_COVARIANT_STRAIN_ROW ( SIXTH         , SIXTH         , T, 2, 3, .TRUE., .TRUE., EST_C )
         CALL MITC3_COVARIANT_STRAIN_ROW ( THIRD+DTIE    , THIRD-TWO*DTIE, T, 1, 3, .TRUE., .TRUE., ERT_D )
         CALL MITC3_COVARIANT_STRAIN_ROW ( THIRD-TWO*DTIE, THIRD+DTIE    , T, 2, 3, .TRUE., .TRUE., EST_E )
         CALL MITC3_COVARIANT_STRAIN_ROW ( THIRD+DTIE    , THIRD+DTIE    , T, 1, 3, .TRUE., .TRUE., ERT_F )
         CALL MITC3_COVARIANT_STRAIN_ROW ( THIRD+DTIE    , THIRD+DTIE    , T, 2, 3, .TRUE., .TRUE., EST_F )

         CHAT = ERT_F - ERT_D - EST_F + EST_E

         BCOV(6,:) = TWO3*( ERT_B - HALF*EST_B ) + THIRD*( ERT_C + EST_C ) + THIRD*CHAT*(THREE*S - ONE)   ! e_rt
         BCOV(5,:) = TWO3*( EST_A - HALF*ERT_A ) + THIRD*( ERT_C + EST_C ) + THIRD*CHAT*(ONE - THREE*R)   ! e_st

      ENDIF

! Transform covariant components to the element Cartesian system: eps_kl = e_ij (g^i)_k (g^j)_l

      CALL MITC3_COVARIANT_BASIS ( R, S, T, G )
      CALL MITC_CONTRAVARIANT_BASIS ( G, GCON )

      B(:,:) = ZERO
      DO ROW=1,6
         K = INDEX1(ROW)
         L = INDEX2(ROW)
         DO I=1,3
            DO J=1,3
               B(ROW,:) = B(ROW,:) + GCON(K,I)*GCON(L,J)*BCOV(ROWS(I,J),:)
            ENDDO
         ENDDO
      ENDDO

      B(4:6,:) = TWO*B(4:6,:)                              ! Engineering shear strains

      RETURN

      END SUBROUTINE MITC3_B

! ##################################################################################################################################

      SUBROUTINE MITC3_BMBS ( R, S, BM, BB, BS )

! Mid-surface membrane (BM: xx, yy, xy), curvature (BB: kxx, kyy, kxy) and transverse shear (BS: zx, yz) strain-displacement
! operators in element coordinates, extracted from the 3-D strain-displacement matrix in the same way as MITC4_BMBS: the membrane
! part at T = 0, the curvature as (top - bottom)/thickness of the bending-only part, and the shear at T = 0. Rows and sign
! conventions are those of SHELL_A/B/D/T and of ELMDIS_PLY (ply strain = e0 + z*k).

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE
      USE MODEL_STUF, ONLY            :  EPROP
      USE CONSTANTS_1, ONLY           :  ZERO, ONE

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: R, S
      REAL(DOUBLE), INTENT(OUT)       :: BM(3,NDOFB)
      REAL(DOUBLE), INTENT(OUT)       :: BB(3,NDOFB)
      REAL(DOUBLE), INTENT(OUT)       :: BS(2,NDOFB)

      REAL(DOUBLE)                    :: BMEM(6,NDOFB), BTOP(6,NDOFB), BBOT(6,NDOFB), BSHR(6,NDOFB)

! **********************************************************************************************************************************

      CALL MITC3_B ( R, S, ZERO, .TRUE., .FALSE., .FALSE., BMEM )
      BM(1,:) = BMEM(1,:)
      BM(2,:) = BMEM(2,:)
      BM(3,:) = BMEM(4,:)

      CALL MITC3_B ( R, S, -ONE, .FALSE., .TRUE., .FALSE., BBOT )
      CALL MITC3_B ( R, S, +ONE, .FALSE., .TRUE., .FALSE., BTOP )
      BB(1,:) = (BTOP(1,:) - BBOT(1,:)) / EPROP(1)
      BB(2,:) = (BTOP(2,:) - BBOT(2,:)) / EPROP(1)
      BB(3,:) = (BTOP(4,:) - BBOT(4,:)) / EPROP(1)

      CALL MITC3_B ( R, S, ZERO, .FALSE., .FALSE., .TRUE., BSHR )
      BS(1,:) = BSHR(6,:)                                  ! zx
      BS(2,:) = BSHR(5,:)                                  ! yz

      RETURN

      END SUBROUTINE MITC3_BMBS

! ##################################################################################################################################

      SUBROUTINE MITC3_THERM_CURV ( R, S, KAP )

! Free thermal curvature per unit temperature at one point of the mid-surface, in element coordinates, in the same basis and sign
! convention as the curvature operator BB of MITC3_BMBS. This is the MITC3+ counterpart of THERM_CURV_AT_RS in subr MITC4 (see
! there for the derivation): with the free strain tensor e_th = [A]^-1 {A_alpha} acting on the mid-surface tangents g_i,
!
!    q_ij = 1/4 ( (e_th g_i) . w_j + (e_th g_j) . w_i ),     w_j = g_j(T=+1) - g_j(T=-1)
!
! carried to element coordinates with the contravariant basis at T = +1 and T = -1 (top minus bottom) and divided by the thickness.
! For a flat triangle whose directors are the facet normal, w_j = 0 and KAP is identically zero; it is nonzero only with SNORM
! directors that are not normal to the facet (faceted curved shells). KAP is also returned zero if [A] is singular.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  EPROP, SHELL_A, SHELL_AALP
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, ONE, TWO, QUARTER

      IMPLICIT NONE

      REAL(DOUBLE), INTENT(IN)        :: R, S
      REAL(DOUBLE), INTENT(OUT)       :: KAP(3)

      INTEGER(LONG)                   :: IA, IB, ID, KK, LL
      REAL(DOUBLE)                    :: SIDE, DETA
      REAL(DOUBLE)                    :: ADJA(3,3), EPS_TH(3), EPS_TEN(3,3)
      REAL(DOUBLE)                    :: GMID(3,3), GTOP(3,3), GBOT(3,3), GSIDE(3,3), GCON(3,3)
      REAL(DOUBLE)                    :: DIRW(3,2), EPSX(3,2), QTH(2,2), KTEN(3,3)

! **********************************************************************************************************************************

      KAP(:) = ZERO

      CALL MITC3_COVARIANT_BASIS ( R, S, +ONE, GTOP )
      CALL MITC3_COVARIANT_BASIS ( R, S, -ONE, GBOT )
      DIRW(:,1) = GTOP(:,1) - GBOT(:,1)
      DIRW(:,2) = GTOP(:,2) - GBOT(:,2)
      IF (MAXVAL(DABS(DIRW)) <= 1.0D-14*EPROP(1)) RETURN   ! Flat element with directors normal to it: no thermal curvature

      DETA = SHELL_A(1,1)*(SHELL_A(2,2)*SHELL_A(3,3) - SHELL_A(2,3)*SHELL_A(3,2))                                                  &
           - SHELL_A(1,2)*(SHELL_A(2,1)*SHELL_A(3,3) - SHELL_A(2,3)*SHELL_A(3,1))                                                  &
           + SHELL_A(1,3)*(SHELL_A(2,1)*SHELL_A(3,2) - SHELL_A(2,2)*SHELL_A(3,1))

      IF (DABS(DETA) <= 1.0D-12*DABS(SHELL_A(1,1)*SHELL_A(2,2)*SHELL_A(3,3))) RETURN

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

      EPS_TEN(:,:) = ZERO
      EPS_TEN(1,1) = EPS_TH(1)
      EPS_TEN(2,2) = EPS_TH(2)
      EPS_TEN(1,2) = HALF*EPS_TH(3)
      EPS_TEN(2,1) = HALF*EPS_TH(3)

      CALL MITC3_COVARIANT_BASIS ( R, S, ZERO, GMID )
      DO IA=1,2
         EPSX(:,IA) = MATMUL( EPS_TEN, GMID(:,IA) )
      ENDDO

      DO IA=1,2
         DO IB=1,2
            QTH(IA,IB) = QUARTER*( DOT_PRODUCT( EPSX(:,IA), DIRW(:,IB) ) + DOT_PRODUCT( EPSX(:,IB), DIRW(:,IA) ) )
         ENDDO
      ENDDO

      KTEN(:,:) = ZERO
      DO ID=1,2
         SIDE = ONE
         IF (ID == 2) SIDE = -ONE
         CALL MITC3_COVARIANT_BASIS ( R, S, SIDE, GSIDE )
         CALL MITC_CONTRAVARIANT_BASIS ( GSIDE, GCON )
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
      KAP(3) = TWO*KTEN(1,2)/EPROP(1)

      RETURN

      END SUBROUTINE MITC3_THERM_CURV

   END MODULE MITC3_MOD
