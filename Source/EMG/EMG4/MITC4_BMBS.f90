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
      SUBROUTINE MITC4_BMBS ( R, S, BM, BB, BS )

! Extracts the mid-surface membrane, curvature, and transverse shear strain-displacement
! operators for the MITC4/MITC4+ element at in-plane location (R,S), for use by an
! A/B/D/T (classical-lamination-theory style) stiffness/load assembly.
!
! MITC4's strain-displacement matrix B(R,S,T) (see MITC4_B) is linear in the through-thickness
! isoparametric coordinate T for the in-plane strain components (rows 1,2,4), and the transverse
! shear components (rows 5,6) are already evaluated independent of T inside MITC4_B (constant
! through thickness, per the MITC tying-point construction). This means the membrane, curvature,
! and shear operators can be recovered exactly by evaluating B at the top (T=+1) and bottom
! (T=-1) surfaces and combining, with no loss of information relative to full through-thickness
! (R,S,T) Gauss integration.
!
! This is the same extraction already used (independently, inline) by MITC4's stress/strain
! recovery and differential-stiffness (KED) logic; it is factored out here so that the main
! stiffness (KE) and thermal-load (PTE) assembly can reuse it too, as part of converting MITC4
! from a volume-integrated formulation to the A/B/D/T idealization used elsewhere in MYSTRAN
! (MIN4T, TREL1/TPLT2, QDEL1/QPLT3). See the MITC4/MITC4+ PCOMP conversion plan.
!
! NOTE (carried over from the original inline code): using EPROP(1) as the through-thickness
! divisor for curvature assumes a uniform thickness across the element. For a PSHELL or PCOMP
! with grid-point-varying thickness this should strictly be the thickness interpolated at
! (R,S), not the element-constant EPROP(1). Not addressed by this extraction; flagged for
! follow-up.
!
! NOTE: the local Cartesian basis used to transform contravariant strains (via
! MITC4_CARTESIAN_LOCAL_BASIS) is evaluated separately at T=+1 and T=-1 and is, in general,
! not the same at both surfaces for a warped element. This is an existing, documented
! approximation (see original stress-recovery comment) and is preserved as-is here.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  ELGP, EPROP
      USE CONSTANTS_1, ONLY           :  ONE, TWO

      USE MITC4_B_Interface
      USE MITC4_CARTESIAN_LOCAL_BASIS_Interface
      USE MITC_TRANSFORM_B_Interface

      IMPLICIT NONE

      REAL(DOUBLE) , INTENT(IN)       :: R, S              ! Isoparametric in-plane coordinates
      REAL(DOUBLE) , INTENT(OUT)      :: BM(3, 6*ELGP)     ! Membrane strain-disp matrix   (rows: xx, yy, xy)
      REAL(DOUBLE) , INTENT(OUT)      :: BB(3, 6*ELGP)     ! Curvature strain-disp matrix  (rows: xx, yy, xy)
      REAL(DOUBLE) , INTENT(OUT)      :: BS(2, 6*ELGP)     ! Trans shear strain-disp matrix(rows: zx, yz)

      REAL(DOUBLE)                    :: BI1(6, 6*ELGP)    ! Full strain-displ matrix at T = -1 (bottom)
      REAL(DOUBLE)                    :: BI2(6, 6*ELGP)    ! Full strain-displ matrix at T = +1 (top)
      REAL(DOUBLE)                    :: TRANSFORM(3,3)    ! Contravariant-to-local-cartesian transform


! **********************************************************************************************************************************
! Get the strain-displacement matrix at top and bottom and transform strain terms to the
! element (local Cartesian) coordinate system. Identical procedure to what MITC4's stress
! recovery previously did inline.

      CALL MITC4_B( R, S, -ONE, .TRUE., .TRUE., .TRUE., BI1 )
      TRANSFORM = MITC4_CARTESIAN_LOCAL_BASIS( R, S, -ONE )
      BI1(4:6,:) = BI1(4:6,:) / 2
      CALL MITC_TRANSFORM_B( TRANSFORM, BI1 )
      BI1(4:6,:) = BI1(4:6,:) * 2

      CALL MITC4_B( R, S, +ONE, .TRUE., .TRUE., .TRUE., BI2 )
      TRANSFORM = MITC4_CARTESIAN_LOCAL_BASIS( R, S, +ONE )
      BI2(4:6,:) = BI2(4:6,:) / 2
      CALL MITC_TRANSFORM_B( TRANSFORM, BI2 )
      BI2(4:6,:) = BI2(4:6,:) * 2


! **********************************************************************************************************************************
! Membrane strain is the average of the strains at the two T points.

      BM(1,1:6*ELGP) = (BI2(1,:) + BI1(1,:)) / TWO           ! xx
      BM(2,1:6*ELGP) = (BI2(2,:) + BI1(2,:)) / TWO           ! yy
      BM(3,1:6*ELGP) = (BI2(4,:) + BI1(4,:)) / TWO           ! xy


! **********************************************************************************************************************************
! Curvature is (strain_top - strain_bottom) / thickness.

      BB(1,1:6*ELGP) = (BI2(1,:) - BI1(1,:)) / EPROP(1)      ! xx
      BB(2,1:6*ELGP) = (BI2(2,:) - BI1(2,:)) / EPROP(1)      ! yy
      BB(3,1:6*ELGP) = (BI2(4,:) - BI1(4,:)) / EPROP(1)      ! xy


! **********************************************************************************************************************************
! Transverse shear strain. Note reversed order of rows (zx then yz) to match SHELL_T convention.

      BS(1,1:6*ELGP) = (BI2(6,:) + BI1(6,:)) / TWO           ! zx
      BS(2,1:6*ELGP) = (BI2(5,:) + BI1(5,:)) / TWO           ! yz

      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE MITC4_BMBS
