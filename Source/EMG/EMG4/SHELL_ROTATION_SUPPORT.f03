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

   MODULE SHELL_ROTATION_SUPPORT

   USE COMPOSITE_SHELL_PREPARATION, ONLY :  GET_PCOMP_SECT_PROPS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CALC_K6ROT, CALC_PHI_SQ

   CONTAINS
      SUBROUTINE CALC_K6ROT ()

! Builds the stiffness matrix for a spring connecting the drilling DOF to the translational DOFs of the adjacent nodes
! of the element. Adds that to the existing stiffness matrix KE in the element coordinate system.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE MODEL_STUF, ONLY            :  TYPE, ELGP, INTL_MID, XEL, SHELL_A, KE, HBAR
      USE PARAMS, ONLY                :  K6ROT, QUAD4TYP
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE SCONTR, ONLY                :  MAX_ORDER_GAUSS
      USE JACOBIAN, ONLY               :  JAC2D
      USE QUADRATURE, ONLY            :  ORDER_GAUSS
      USE VECTOR_GEOMETRY, ONLY       :  CROSS
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF_T

      IMPLICIT NONE

      REAL(DOUBLE)                    :: KROT(6*ELGP,6*ELGP)  ! Stifness matrix of K6ROT.
      REAL(DOUBLE)                    :: N(3)                 ! Normal at one grid point
      REAL(DOUBLE)                    :: X_PREV(3)            ! Coordinates of the previous grid point (n-1)
      REAL(DOUBLE)                    :: X_NEXT(3)            ! Coordinates of the next grid point (n+1)
      REAL(DOUBLE)                    :: TERM_PREV(3)
      REAL(DOUBLE)                    :: TERM_NEXT(3)
      REAL(DOUBLE)                    :: ZG(ELGP)             ! Height of each grid point above the element x-y plane
      REAL(DOUBLE)                    :: B(6*ELGP)            ! Strain-displacement matrix for one grid point's K6ROT spring "element"
      REAL(DOUBLE)                    :: STIFFNESS            ! Spring stiffness of one grid point's K6ROT "element"
      REAL(DOUBLE)                    :: AREA                 ! Elem area
      REAL(DOUBLE)                    :: DETJ                 ! An output from subr JAC2D4, called herein. Determinant of JAC
      INTEGER(LONG)                   :: I,J
      INTEGER(LONG)                   :: GP                   ! Element grid point number (1 to ELGP).
      INTEGER(LONG)                   :: GP_PREV
      INTEGER(LONG)                   :: GP_NEXT
      REAL(DOUBLE)                    :: JAC(2,2)             ! An output from subr JAC2D4, called herein. 2 x 2 Jacobian matrix.
      REAL(DOUBLE)                    :: JACI(2,2)            ! An output from subr JAC2D4, called herein. 2 x 2 Jacobian inverse.
      REAL(DOUBLE)                    :: X2E                  ! x coord of elem node 2
      REAL(DOUBLE)                    :: Y3E                  ! y coord of elem node 3
      REAL(DOUBLE)                    :: XSD(4)               ! Diffs in x coords of quad sides in local coords
      REAL(DOUBLE)                    :: YSD(4)               ! Diffs in y coords of quad sides in local coords
      REAL(DOUBLE)                    :: HHH(MAX_ORDER_GAUSS) ! An output from subr ORDER, called herein.  Gauss weights.
      REAL(DOUBLE)                    :: SSS(MAX_ORDER_GAUSS) ! An output from subr ORDER, called herein. Gauss abscissa's.


! **********************************************************************************************************************************


                                                     ! No K6ROT for shells that only use MID1.
      IF (INTL_MID(2) > 0) THEN

         AREA = ZERO

         IF ((TYPE(1:5) == "QUAD4")) THEN

            XSD(1) = XEL(1,1) - XEL(2,1)             ! x coord diffs (in local elem coords)
            XSD(2) = XEL(2,1) - XEL(3,1)
            XSD(3) = XEL(3,1) - XEL(4,1)
            XSD(4) = XEL(4,1) - XEL(1,1)

            YSD(1) = XEL(1,2) - XEL(2,2)             ! y coord diffs (in local elem coords)
            YSD(2) = XEL(2,2) - XEL(3,2)
            YSD(3) = XEL(3,2) - XEL(4,2)
            YSD(4) = XEL(4,2) - XEL(1,2)

            CALL ORDER_GAUSS ( 2, SSS, HHH )
            DO I=1,2
               DO J=1,2
                  CALL JAC2D ( SSS(I), SSS(J), XSD, YSD, 'N', JAC, JACI, DETJ )
                  AREA = AREA + HHH(I)*HHH(J)*DETJ
               ENDDO
            ENDDO

         ELSEIF (TYPE(1:5) == "TRIA3") THEN

            X2E  = XEL(2,1)
            Y3E  = XEL(3,2)
                                                     ! Actual area is half this but using this value
                                                     ! gives the same stiffness as MSC.
            AREA = X2E*Y3E

         ENDIF


         !According to:
         !https://www.dynalook.com/conferences/9th-european-ls-dyna-conference/drilling-rota
         !tion-constraint-for-shell-elements-in-implicit-and-explicit-analyses
         !but scaled to match MSC Nastran.

         ! SHELL_A(3,3) is membrane shear modulus times thickness
         STIFFNESS = 10.0**(-6.0) * K6ROT * SHELL_A(3,3) * ABS(AREA)

         ! Use the uniform element normal at all grid points.
         ! This might be supposed to be the shell normal but it hardly seems to make a difference and this way is simpler.
         N = [ZERO, ZERO, ONE]

         ! Grid heights above the element x-y plane. XEL holds them for MITC4; for MIN4/MIN4T XEL is the mean plane and the
         ! grids of a warped quad are at -HBAR, +HBAR, -HBAR, +HBAR from it.
         ZG = XEL(1:ELGP,3)
         IF ((TYPE(1:5) == 'QUAD4') .AND. (QUAD4TYP(1:4) == 'MIN4')) THEN
            ZG = [-HBAR, HBAR, -HBAR, HBAR]
         ENDIF

         DO GP=1,ELGP

            B = ZERO

            ! The spring at one grid point is effectively a 3-node element with nodes and DOFs of:
            ! Node n:          tx, ty, tz, rx, ry, rz
            ! Node n+1 (next): tx, ty, tz
            ! Node n-1 (prev): tx, ty, tz

            GP_PREV = GP - 1
            IF (GP_PREV < 1) THEN
               GP_PREV = ELGP
            ENDIF

            GP_NEXT = GP + 1
            IF (GP_NEXT > ELGP) THEN
               GP_NEXT = 1
            ENDIF

            X_PREV = [XEL(GP_PREV,1:2) - XEL(GP,1:2), ZG(GP_PREV) - ZG(GP)]
            X_NEXT = [XEL(GP_NEXT,1:2) - XEL(GP,1:2), ZG(GP_NEXT) - ZG(GP)]

            !Contribution of previous node's displacement
            !        - n × (x_n-1 - x_n)
            !ε_n  =  ------------------- * u_n-1
            !        2 * |x_n-1 - x_n|^2
            CALL CROSS(N, X_PREV, TERM_PREV)
            TERM_PREV = TERM_PREV * 1 / (2 * (X_PREV(1)**2 + X_PREV(2)**2 + X_PREV(3)**2) )

            B((GP_PREV - 1) * 6 + 1) = -TERM_PREV(1)
            B((GP_PREV - 1) * 6 + 2) = -TERM_PREV(2)
            B((GP_PREV - 1) * 6 + 3) = -TERM_PREV(3)

            !Contribution of next node's displacement
            !        - n × (x_n+1 - x_n)
            !ε_n += ------------------- * u_n+1
            !        2 * |x_n+1 - x_n|^2
            CALL CROSS(N, X_NEXT, TERM_NEXT)
            TERM_NEXT = TERM_NEXT * 1 / (2 * (X_NEXT(1)**2 + X_NEXT(2)**2 + X_NEXT(3)**2) )

            B((GP_NEXT - 1) * 6 + 1) = -TERM_NEXT(1)
            B((GP_NEXT - 1) * 6 + 2) = -TERM_NEXT(2)
            B((GP_NEXT - 1) * 6 + 3) = -TERM_NEXT(3)

            !Contribution of current node's displacement and rotation
            !          n × (x_n-1 - x_n)     n × (x_n+1 - x_n)
            !ε_n += ( ------------------- + ------------------- ) * u_n  +  n · r_n
            !         2 * |x_n-1 - x_n|^2   2 * |x_n+1 - x_n|^2
            B((GP - 1) * 6 + 1) = TERM_PREV(1) + TERM_NEXT(1)
            B((GP - 1) * 6 + 2) = TERM_PREV(2) + TERM_NEXT(2)
            B((GP - 1) * 6 + 3) = TERM_PREV(3) + TERM_NEXT(3)
            !Rotation of node n. On a warped element the edges are not normal to n, and a rigid rotation w about an in-plane axis
            !moves the next and previous nodes by w x (x - x_n), which the terms above take as a rotation about n of
            !  sum over the 2 edges of  (n . (x - x_n)) (w . (x - x_n)) / (2 |x - x_n|^2).
            !Subtracting that through the rotation of node n makes the spring strain zero for every rigid motion (it is n . r_n
            !alone when the edges lie in the plane).
            B((GP - 1) * 6 + 4:(GP - 1) * 6 + 6) = N - DOT_PRODUCT(N, X_PREV) * X_PREV / (2 * DOT_PRODUCT(X_PREV, X_PREV))   &
                                                     - DOT_PRODUCT(N, X_NEXT) * X_NEXT / (2 * DOT_PRODUCT(X_NEXT, X_NEXT))

            ! stiffness * B' * B
            CALL MATMULT_FFF_T(B, B, 1, 6*ELGP, 6*ELGP, KROT)
            KROT = KROT * STIFFNESS

            KE(1:6*ELGP, 1:6*ELGP) = KE(1:6*ELGP, 1:6*ELGP) + KROT(1:6*ELGP, 1:6*ELGP)

         ENDDO

      ENDIF


! **********************************************************************************************************************************

      END SUBROUTINE CALC_K6ROT


      SUBROUTINE CALC_PHI_SQ ( IERROR )

! Calculates the PHI_SQ shear correction factor used in the MIN3 and MIN4 elements (TRIA3, QUAD4 Mindlin elements) in subrs
! TRPLT2, QPLT2. PHI_SQ is used as a multiplier of the shear stiffness and the element stresses

!   PHI_SQ  = CBMIN*PSI_HAT/(1 + CBMIN*PSI_HAT) where

!   PSI_HAT = BENSUM/SHRSUM

!   (1) CBMIN is a constant from the published papers by Alexander TEssler (see refs in the MYSTRAN User's Manual) for the MIN3
!       triangular element (TRIA3) and MIN4 quadrilateral element (QUAD4)

!   (2) BENSUM is the sum of the diagonal terms for rotation DOF's from the bending portion of the element stiffness matrix

!   (3) SHRSUM is the sum of the diagonal terms from the transverse shear portion of the element stiffness matrix

! BENSUM and SHRSUM are determined in subrs TPLT2 (for the MIN3 TRIA3 elem) and QPLT2 (for the MIN4 QUAD4 elem) prior to calling
! this subr.

! For composite elements (PCOMP properties), a correction to PHI_SQ has to be made when PHI_SQ is calculated during stress recovery
! of individual plies. In these cases, the PHI_SQ calculated when the stiffness and stress recovery matrices is based on the
! individual ply thicknesses (bending and transverse shear) but should really be based on the whole element - i.e. the integrated
! effect of all plies. It is known that BENSUM is proportional to the bending inertia of the element and SHRSUM is proportional to
! the transverse shear thickness of the element. Thus, for composite elements, the BENSUM and SHRSUM calculated for 1 ply are
! adjusted by the ratio of these quantities between 1 ply and the composite layup


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MEFE
      USE IOUNT1, ONLY                :  ERR, F06, WRT_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWELVE
      USE PARAMS, ONLY                :  CBMIN3, CBMIN4, CBMIN4T, EPSIL, PCMPTSTM, QUAD4TYP
      USE MODEL_STUF, ONLY            :  BENSUM, EID, EMG_IFE, EMG_RFE, ERR_SUB_NAM, NUM_EMG_FATAL_ERRS, INTL_MID, PHI_SQ,  &
                                         PCOMP_PROPS, PLY_NUM, PSI_HAT, SHRSUM, TPLY, TYPE
      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_PHI_SQ'

      INTEGER(LONG), INTENT(OUT)      :: IERROR            ! Local error indicator


      REAL(DOUBLE)                    :: CBMIN  = ZERO     ! Either CBMIN3 or CBMIN4
      REAL(DOUBLE)                    :: DEN               ! Denominator term in calculating PHI_SQ
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare to real zero
      REAL(DOUBLE)                    :: PCOMP_TM          ! Membrane thick of PCOMP for equiv PSHELL (see subr SHELL_ABD_MATRICES)
      REAL(DOUBLE)                    :: PCOMP_IB          ! Bending MOI of PCOMP for equiv PSHELL (see subr SHELL_ABD_MATRICES)
      REAL(DOUBLE)                    :: PCOMP_TS          ! Transv shear thick of PCOMP for equiv PSHELL (subr SHELL_ABD_MATRICES)
      REAL(DOUBLE)                    :: PLY_IB            ! Bending MOI of a ply
      REAL(DOUBLE)                    :: PLY_TS            ! Transv shear thick of a ply



! **********************************************************************************************************************************
      IERROR = 0
      EPS1   = EPSIL(1)

! For composite elements, adjust BENSUM and SHRSUM to be based on the whole composite element properties since they were calculated
! based on only 1 ply

      IF (PCOMP_PROPS == 'Y') THEN
         IF (PLY_NUM /= 0) THEN
            CALL GET_PCOMP_SECT_PROPS ( PCOMP_TM, PCOMP_IB, PCOMP_TS )
            PLY_IB = TPLY*TPLY*TPLY/TWELVE
            PLY_TS = TPLY*PCMPTSTM
            BENSUM = (PCOMP_IB/PLY_IB)*BENSUM
            SHRSUM = (PCOMP_TS/PLY_TS)*SHRSUM
         ENDIF
      ENDIF

! Now calculate PHI_SQ (but only if the element has transverse shear flexibility)

      IF      (TYPE(1:5) == 'QUAD4') THEN                  ! Regardless of who calls this subr, TYPE will det which CBMIN to use
         IF (QUAD4TYP == 'MIN4T ') THEN                    ! quad TYPE's will use either CBMIN4 or CBMIN4T
            CBMIN = CBMIN4T
         ELSE
            CBMIN = CBMIN4
         ENDIF
      ELSE IF (TYPE(1:5) == 'TRIA3') THEN                  ! tria TYPE's will use CBMIN3 (NOTE: MIN4T quads will use CBMIN4T, above)
         CBMIN = CBMIN3
      ELSE
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IERROR  = IERROR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1943) TYPE, EID
            WRITE(F06,1943) TYPE, EID
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1943
            ENDIF
         ENDIF
      ENDIF

      IERROR  = 0
      PSI_HAT = ZERO
      IF (DABS(SHRSUM) > EPS1) THEN
         PSI_HAT = BENSUM/SHRSUM
      ELSE
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IERROR  = IERROR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1923) SHRSUM, TYPE, EID
            WRITE(F06,1923) SHRSUM, TYPE, EID
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1923
               EMG_RFE(NUM_EMG_FATAL_ERRS,1) = SHRSUM
            ENDIF
         ENDIF
      ENDIF

      DEN = ONE + CBMIN*PSI_HAT
      IF (DABS(DEN) > EPS1) THEN
         PHI_SQ = CBMIN*PSI_HAT/DEN
      ELSE
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IERROR  = IERROR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1942) DEN, TYPE, EID
            WRITE(F06,1942) DEN, TYPE, EID
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1942
               EMG_RFE(NUM_EMG_FATAL_ERRS,1) = DEN
            ENDIF
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1923 FORMAT(' *ERROR  1923: SHRSUM PARAMETER = ',1ES9.2,' IS TOO CLOSE TO ZERO FOR ',A,' ELEMENT ',I8                             &
                    ,/,14X,' CANNOT CALCULATE PSI_HAT FACTOR NEEDED FOR TRANSVERSE SHEAR STIFFNESS CALCULATION')

 1942 FORMAT(' *ERROR  1942: 1 + CBMIN*PSI_HAT = ',1ES9.2,' IS TOO CLOSE TO ZERO FOR ',A,' ELEMENT ',I8                            &
                    ,/,14X,' CANNOT CALCULATE PHI_SQ FACTOR NEEDED FOR TRANSVERSE SHEAR STIFFNESS CALCULATION')

 1943 FORMAT(' *ERROR  1943: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,15X,  A, 'ELEMENT ',I8,' IS NOT A VALID ELEMENT FOR THIS SUBR')

! **********************************************************************************************************************************

      END SUBROUTINE CALC_PHI_SQ

   END MODULE SHELL_ROTATION_SUPPORT
