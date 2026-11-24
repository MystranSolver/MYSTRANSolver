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

   MODULE LINE_ELEMENTS

   USE RESULT_COORDINATES, ONLY: ELMDIS
   USE ELEMENT_RECOVERY_SUPPORT, ONLY: CALC_ELEM_NODE_FORCES

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BREL1

   CONTAINS

      SUBROUTINE BREL1 ( OPT, WRITE_WARN )

! Calculates, or calls subr's to calculate, quadrilateral element matrices:

!  1) ME        = element mass matrix                  , if OPT(1) = 'Y'
!  2) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  3) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  4) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  5) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO
      USE DEBUG_PARAMETERS
      USE MODEL_STUF, ONLY            :  DOFPIN, EID, ELEM_LEN_AB, EMAT, NUM_EMG_FATAL_ERRS, EPROP, FCONV, ME, ULT_STRE,      &
                                         ULT_STRN, TYPE, ZS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BREL1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not
      CHARACTER(1*BYTE)               :: CALC_BAR          ! 'Y' if the BAR stiffness routine needs to be called

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: NUM_PFLAG_DOFS    ! The number of pin flagged DOF's for this element (from BAR1/BART)



      REAL(DOUBLE)                    :: ALPHA             ! Coefficient of thermal expansion
      REAL(DOUBLE)                    :: AREA              ! Cross-sectional area
      REAL(DOUBLE)                    :: E                 ! Youngs modulus
      REAL(DOUBLE)                    :: G                 ! Shear modulus
      REAL(DOUBLE)                    :: GE                ! Material damping coeff
      REAL(DOUBLE)                    :: I1                ! Bending inertia in plane 1
      REAL(DOUBLE)                    :: I12               ! Product of inertia
      REAL(DOUBLE)                    :: I2                ! Bending inertia in plane 2
      REAL(DOUBLE)                    :: JTOR              ! Torsional constant
      REAL(DOUBLE)                    :: K1                ! Shear constant for plane 1 (used in K1*AREA*G)
      REAL(DOUBLE)                    :: K2                ! Shear constant for plane 2 (used in K1*AREA*G)
      REAL(DOUBLE)                    :: M0                ! Intermediate variable in calculating element mass matrix, ME
      REAL(DOUBLE)                    :: NSM               ! Nonstructural mass
      REAL(DOUBLE)                    :: RHO               ! Material density
      REAL(DOUBLE)                    :: TREF              ! Element reference temperature
      REAL(DOUBLE)                    :: TPIN(12,12)       ! Pin flag transformation from BAR1/BART



! **********************************************************************************************************************************

! Set element property and material constants

      IF (TYPE == 'ROD     ') THEN

         AREA     = EPROP(1)                               ! Cross-sectional area
         JTOR     = EPROP(2)                               ! Torsional constant
         ZS(1)    = EPROP(3)                               ! C (Tors. stress recovery coeff) on PROD
         NSM      = EPROP(4)                               ! Non-structural mass
         FCONV(1) = AREA

      ELSE IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BART    ')) THEN

         AREA     = EPROP( 1)                              ! Cross-sectional area
         I1       = EPROP( 2)                              ! Plane 1 moment of inertia
         I2       = EPROP( 3)                              ! Plane 2 moment of inertia
         JTOR     = EPROP( 4)                              ! Torsional constant
         NSM      = EPROP( 5)                              ! Non-structural mass
         ZS(1)    = EPROP( 6)                              ! y coord of 1st point for stress recovery
         ZS(2)    = EPROP( 7)                              ! z coord of 1st point for stress recovery
         ZS(3)    = EPROP( 8)                              ! y coord of 2nd point for stress recovery
         ZS(4)    = EPROP( 9)                              ! z coord of 2nd point for stress recovery
         ZS(5)    = EPROP(10)                              ! y coord of 3rd point for stress recovery
         ZS(6)    = EPROP(11)                              ! z coord of 3rd point for stress recovery
         ZS(7)    = EPROP(12)                              ! y coord of 4th point for stress recovery
         ZS(8)    = EPROP(13)                              ! z coord of 4th point for stress recovery
         K1       = EPROP(14)                              ! Plane 1 shear factor
         K2       = EPROP(15)                              ! Plane 2 shear factor
         I12      = EPROP(16)                              ! Product of inertia
         ZS(9)    = EPROP(17)                              ! Torsional stress recovery coefficient
         FCONV(1) = AREA

      ELSE IF (TYPE == 'BEAM    ') THEN

      ENDIF

! Need to set some values for materials here since subr for material properties not called for these 1D elements

      E             = EMAT( 1,1)                           ! Young's modulus
      G             = EMAT( 2,1)                           ! Shear modulus
      RHO           = EMAT( 4,1)                           ! Mass density
      ALPHA         = EMAT( 5,1)                           ! Coefficient of thermal expansion
      TREF          = EMAT( 6,1)                           ! Reference temperature for thermal expaqnsion
      GE            = EMAT( 7,1)                           ! Structural damping coefficient
      ULT_STRE(1,1) = EMAT( 8,1)                           ! Max allowable stress in tension
      ULT_STRE(2,1) = EMAT( 9,1)                           ! Max allowable stress in compression
      ULT_STRE(3,1) = EMAT(10,1)                           ! Max allowable stress in shear

! **********************************************************************************************************************************
! Generate the mass matrix for this element (array was initialized in subr EMG).

      IF (OPT(1) == 'Y') THEN
         M0 = (RHO*AREA + NSM)*(ELEM_LEN_AB)/TWO
         ME(1,1) = M0
         ME(2,2) = M0
         ME(3,3) = M0
         ME(7,7) = M0
         ME(8,8) = M0
         ME(9,9) = M0
      ENDIF

! **********************************************************************************************************************************
! Call routines to calc element matrices (stiffness, etc.).
! For BAR's with pin flags the stiffness routine is also called when only ME is requested, since the pin flag transformation, TPIN,
! is generated from the linear stiffness matrix and must also be applied to ME

      NUM_PFLAG_DOFS = 0

      CALC_BAR = 'N'
      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(5) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         CALC_BAR = 'Y'
      ELSE IF ((OPT(1) == 'Y') .AND. (TYPE == 'BAR     ')) THEN
         DO I=1,12
            IF (DOFPIN(I) > 0) THEN
               CALC_BAR = 'Y'
               EXIT
            ENDIF
         ENDDO
      ENDIF

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(5) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         IF      (TYPE == 'ROD     ') THEN

            CALL ROD1 ( OPT, ELEM_LEN_AB, AREA, JTOR, ZS(1), E, G, ALPHA, TREF )

         ELSE IF (TYPE == 'BEAM    ') THEN                 ! General beam

            CALL BEAM

         ENDIF

      ENDIF

      IF ((CALC_BAR == 'Y') .AND. (TYPE == 'BAR     ')) THEN ! Bernoulli-Euler prismatic beam

         IF (DEBUG(249) == 0) THEN

            CALL BAR1 ( OPT, ELEM_LEN_AB, AREA, I1, I2, JTOR, ZS(9), K1, K2, I12, E, G, ALPHA, TREF, NUM_PFLAG_DOFS, TPIN )

         ELSE
            IF (I12 == ZERO) THEN
               CALL BART ( OPT, ELEM_LEN_AB, AREA, I1, I2, JTOR, ZS(9), K1, K2, I12, E, G, ALPHA, TREF, NUM_PFLAG_DOFS, TPIN )
            ELSE
               WRITE(ERR,1963) EID
               WRITE(F06,1963) EID
               RETURN
            ENDIF

         ENDIF

! Apply the pin flag condensation to the mass matrix

         IF ((OPT(1) == 'Y') .AND. (NUM_PFLAG_DOFS > 0)) THEN
            CALL APPLY_PINFLG ( TPIN, ME(1:12,1:12), 12 )
         ENDIF

      ENDIF

! **********************************************************************************************************************************
 1963 FORMAT(' *ERROR  1962: TIMOSHENKO BAR ELEMENT ',A,' CANNOT HAVE NONZERO I12. IT WILL BE SET TO I12 = 0.')



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BREL1


      SUBROUTINE ROD1 ( OPT, L, AREA, JTOR, SCOEFF, E, G, ALPHA, TREF )

! Calculates, for 1-D ROD element

!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  4) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  NTSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  TWO
      USE MODEL_STUF, ONLY            :  DT, KE, KED, PEL, PTE, SE1, STE1, UEL

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ROD1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices
      INTEGER(LONG)                   :: J                 ! Do loop index


      REAL(DOUBLE) , INTENT(IN)       :: ALPHA             ! Coefficient of thermal expansion
      REAL(DOUBLE) , INTENT(IN)       :: AREA              ! Cross-sectional area
      REAL(DOUBLE) , INTENT(IN)       :: E                 ! Youngs modulus
      REAL(DOUBLE) , INTENT(IN)       :: G                 ! Shear modulus
      REAL(DOUBLE) , INTENT(IN)       :: JTOR              ! Torsional constant
      REAL(DOUBLE) , INTENT(IN)       :: L                 ! Elem length
      REAL(DOUBLE) , INTENT(IN)       :: SCOEFF            ! Stress recovery coeff for torsion
      REAL(DOUBLE) , INTENT(IN)       :: TREF              ! Element reference temperature
      REAL(DOUBLE)                    :: C01               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: C02               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: CT0               ! Intermediate variable used in calc thermal loads, PTE
      REAL(DOUBLE)                    :: Fx                ! Axizl load in the ROD element as a result of the applied loads
      REAL(DOUBLE)                    :: KD0               ! Intermediate variable used in calc KED
      REAL(DOUBLE)                    :: TBAR              ! Average elem temperature



! **********************************************************************************************************************************
! Determine element thermal loads.

      IF (OPT(2) == 'Y') THEN
         DO J=1,NTSUB
            TBAR = (DT(1,J) + DT(2,J))/TWO
            CT0 = AREA*E*ALPHA*(TBAR  - TREF)
            PTE( 1,J) = -CT0
            PTE( 7,J) =  CT0
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Calculate SE matrix (2 x 12) for stress data recovery.

      IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         SE1(1, 1,1) = -E/L
         SE1(1, 7,1) =  E/L
         SE1(2, 4,1) = -SCOEFF*G/L                         ! SCOEFF*G/L = (SCOEFF/J)*KE(7,7)
         SE1(2,10,1) =  SCOEFF*G/L
         DO J=1,NTSUB
            TBAR = (DT(1,J) + DT(2,J))/TWO
            STE1(1,J,1) = E*ALPHA*(TBAR - TREF)
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE(12,12).

      IF ((OPT(4) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         C01 = AREA*E/L
         C02 = JTOR*G/L
         KE( 1, 1) = C01
         KE( 1, 7) =-C01
         KE( 7, 1) = KE( 1, 7)
         KE( 4, 4) = C02
         KE( 4,10) =-C02
         KE(10, 4) = KE( 4,10)
         KE( 7, 7) = C01
         KE(10,10) = C02
      ENDIF

! **********************************************************************************************************************************
! Calculate element differential stiffness matrix KED(12,12).

      IF (OPT(6) == 'Y') THEN

         CALL ELMDIS
         CALL CALC_ELEM_NODE_FORCES
         Fx = -PEL(1)                                      ! Fx  (axial force for BAR or ROD) - NASTRAN Fx
         KD0 = Fx/L

         KED( 2, 2) =  KD0
         KED( 2, 8) = -KD0
         KED( 3, 3) =  KD0
         KED( 3, 9) = -KD0

         KED( 8, 2) =  KED( 2, 8)
         KED( 8, 8) =  KED( 2, 2)
         KED( 9, 3) =  KED( 3, 9)
         KED( 9, 9) =  KED( 3, 3)
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE ROD1


      SUBROUTINE BAR1 ( OPT, L, AREA, I1, I2, JTOR, SCOEFF, K1, K2, I12, E, G, ALPHA, TREF, NUM_PFLAG_DOFS, TPIN )

! Calculates, for 1-D BAR element

!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linear stiffness matrix      , if OPT(6) = 'N' (i.e. always calc KE linear unless OPT(6) = 'Y')
!  4) KED       = element differen stiff matrix        , if OPT(6) = 'Y'
!  5) NUM_PFLAG_DOFS, TPIN = number of pin flagged DOF's and the pin flag transformation (TPIN is only defined if NUM_PFLAG_DOFS > 0)


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, NTSUB, BLNK_SUB_NAM, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX, TEN, TWELVE
      USE DEBUG_PARAMETERS
      USE PARAMS, ONLY                :  EPSIL, ART_KED, ART_ROT_KED, ART_TRAN_KED
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ELDOF, DOFPIN, DT, EID, NUM_EMG_FATAL_ERRS, KE, KED, PEL, PTE, SE1, SE2, STE1, STE2, TYPE,&
                                         UEL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF
      USE RESULT_COORDINATES, ONLY    :  ELMDIS
      USE ELEMENT_RECOVERY_SUPPORT, ONLY:  CALC_ELEM_NODE_FORCES

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BAR1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices

      INTEGER(LONG)                   :: I,J               ! DO loop induces
      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG), INTENT(OUT)      :: NUM_PFLAG_DOFS    ! The number of pin flagged DOF's for this element


      REAL(DOUBLE) , INTENT(IN)       :: ALPHA             ! Coefficient of thermal expansion
      REAL(DOUBLE) , INTENT(IN)       :: AREA              ! Cross-sectional area
      REAL(DOUBLE) , INTENT(IN)       :: E                 ! Youngs modulus
      REAL(DOUBLE) , INTENT(IN)       :: G                 ! Shear modulus
      REAL(DOUBLE) , INTENT(IN)       :: I1                ! Bending inertia in plane 1
      REAL(DOUBLE) , INTENT(IN)       :: I12               ! Product of inertia
      REAL(DOUBLE) , INTENT(IN)       :: I2                ! Bending inertia in plane 2
      REAL(DOUBLE) , INTENT(IN)       :: JTOR              ! Torsional constant
      REAL(DOUBLE) , INTENT(IN)       :: K1                ! Shear constant for plane 1 (used in K1*AREA*G)
      REAL(DOUBLE) , INTENT(IN)       :: K2                ! Shear constant for plane 2 (used in K1*AREA*G)
      REAL(DOUBLE) , INTENT(IN)       :: L                 ! Elem length
      REAL(DOUBLE) , INTENT(IN)       :: SCOEFF            ! Stress recovery coeff for torsion
      REAL(DOUBLE) , INTENT(IN)       :: TREF              ! Element reference temperature
      REAL(DOUBLE)                    :: ABAR(6,5)         ! Thermal displacements per unit TPRIME
      REAL(DOUBLE)                    :: BETA              ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: B1(3,6)           ! Intermediate array used in calculating SE1 stress matrix
      REAL(DOUBLE)                    :: B2(3,6)           ! Intermediate array used in calculating SE2 stress matrix
      REAL(DOUBLE)                    :: BT1(3,5)          ! Intermediate array used in calculating STE1 thermal stress effects
      REAL(DOUBLE)                    :: BT2(3,5)          ! Intermediate array used in calculating STE2 thermal stress effects
      REAL(DOUBLE)                    :: BTA(6,5)          ! = KAA*ABAR. Each col gives thermal loads at end a per unit TPRIME value
      REAL(DOUBLE)                    :: BTB(6,5)          ! = KBA*ABAR. Each col gives thermal loads at end a per unit TPRIME value
      REAL(DOUBLE)                    :: C01               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: C02               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: C03               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: DELTA1            ! I2/DEN (1/I1, or 0 if I1 = 0, when DEN = 0)
      REAL(DOUBLE)                    :: DELTA2            ! I1/DEN (1/I2, or 0 if I2 = 0, when DEN = 0)
      REAL(DOUBLE)                    :: DELTA12           ! I12/DEN (0 when DEN = 0)
      REAL(DOUBLE)                    :: D1I1, D2I2        ! DELTA1*I1, DELTA2*I2 (1 when DEN = 0)
      REAL(DOUBLE)                    :: D12I1, D12I2      ! DELTA12*I1, DELTA12*I2 (0 when DEN = 0)
      REAL(DOUBLE)                    :: DEN               ! I1*I2 - I12**2. Only zero if I12 = 0 and I1 and/or I2 = 0
      REAL(DOUBLE)                    :: DUM1(3,NTSUB)     ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM2(3,NTSUB)     ! Intermediate matrix
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare for real zero
      REAL(DOUBLE)                    :: G1                ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: G2                ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: KAA(6,6)          ! Upper left  6 x 6 partition of KE
      REAL(DOUBLE) , INTENT(OUT)      :: TPIN(12,12)       ! Pin flag transformation, u = TPIN*u_r (from PINFLG_TRANSFORM).
                                                           ! Returned so the caller can apply it to other elem matrices (e.g. ME)
      REAL(DOUBLE)                    :: KBA(6,6)          ! Lower left  6 x 6 partition of KE
      REAL(DOUBLE)                    :: KAB(6,6)          ! Upper right 6 x 6 partition of KE
      REAL(DOUBLE)                    :: K11               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: K12               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: K1T               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: K21               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: K22               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: K2T               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: K3T               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: K4T               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: PTA(6,NTSUB)      ! Thermal loads at end a of BAR (rows 1- 6 of PTE)
      REAL(DOUBLE)                    :: PTb(6,NTSUB)      ! Thermal loads at end b of BAR (rows 7-12 of PTE)
      REAL(DOUBLE)                    :: R0                ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: R1                ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: R1D               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: R2                ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: R2D               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: RA                ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: RG                ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: S11(3,6)          ! Intermediate matrix used in calculating SE1 stress matrix
      REAL(DOUBLE)                    :: S12(3,6)          ! Intermediate matrix used in calculating SE1 stress matrix
      REAL(DOUBLE)                    :: S21(3,6)          ! Intermediate matrix used in calculating SE2 stress matrix
      REAL(DOUBLE)                    :: S22(3,6)          ! Intermediate matrix used in calculating SE2 stress matrix
      REAL(DOUBLE)                    :: TBAR              ! Average elem temperature

! The following are used for the differential stiffness matrix calc. See NASTRAN Prog's Manual (COSMIC 1972) page 4.87.30

      REAL(DOUBLE)                    :: M1a               ! Bend mom, plane 1, end a (also Maz, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: M2a               ! Bend mom, plane 2, end a (also May, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: M1b               ! Bend mom, plane 1, end b (also Mbz, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: M2b               ! Bend mom, plane 2, end b (also Mby, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: V1                ! Plane 1 shear (also Vy, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: V2                ! Plane 2 shear (also Vz, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: Fx                ! Axial force (also Fx, in NASTRAN Prog Man differ stiff)

      INTRINSIC DABS

       REAL(DOUBLE)                    :: TPRIME(5,NTSUB)   ! Matrix where each col has the 5 temperature/gradients for the BAR elem



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Calculate element stiffness matrix KE(12,12), regardless of OPT(4) since the calculation of PTE and/or SEi, STEi
! uses values from KE. Initialize KE here to make sure it starts as zero for the BAR element

      DO I=1,12
         DO J=1,12
            KE(I,J) = ZERO
         ENDDO
      ENDDO

      G1   = K1*AREA*G
      G2   = K2*AREA*G
      R0   = TWELVE*E/(L*L*L)
      R1   = R0*I1
      R2   = R0*I2
      BETA = R0*I12
      DEN  = I1*I2 - I12*I12                               !    I1*I2 > I12^2 was checked when PBAR's were read in subr LOADB
      IF (DEBUG(203) > 0) CALL DEBUG_BAR1 ( 1 )

! If I12 > 0 then R1, R2 remain as is. If I12 = 0 then R1, R2 depend on G1 and G2. Quit if R1D or R2D = 0.
! If DEBUG(12) = 1, then use R1 and R2 from below even when I12 = 0

      IERROR = 0
      IF ((I12 == ZERO) .OR. (DEBUG(12) == 1)) THEN
         IF (DABS(G1) >= EPS1) THEN
            R1D = ONE + R1*L/G1
            IF (R1D <= EPS1) THEN
               WRITE(ERR,1918) TYPE, EID
               WRITE(F06,1918) TYPE, EID
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               FATAL_ERR = FATAL_ERR + 1
               IERROR = IERROR + 1
            ELSE
               R1 = R1/R1D
            ENDIF
         ENDIF
         IF (DABS(G2) >= EPS1) THEN
            R2D = ONE + R2*L/G2
            IF (R2D <= EPS1) THEN
               WRITE(ERR,1918) TYPE, EID
               WRITE(F06,1918) TYPE, EID
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               FATAL_ERR = FATAL_ERR + 1
               IERROR = IERROR + 1
            ELSE
               R2 = R2/R2D
            ENDIF
         ENDIF
      ENDIF

! Return if IERROR > 0

      IF (IERROR > 0) RETURN

! Calc intermediate variables

      CALL BAR_INERTIA_FACTORS ( I1, I2, I12, DELTA1, DELTA2, DELTA12, D1I1, D2I2, D12I1, D12I2 )
      RA     = AREA*E/L
      RG     = JTOR*G/L
      K11    = L*L*R1/FOUR
      K12    = E*I1/L
      K21    = L*L*R2/FOUR
      K22    = E*I2/L
      K1T    = K11 + K12
      K2T    = K21 + K22
      K3T    = K11 - K12
      K4T    = K21 - K22
      C01    = L/TWO
      C02    = L*L/THREE
      C03    = C02/TWO
      IF (DEBUG(203) > 0) CALL DEBUG_BAR1 ( 2 )

! Generate KE

      KE( 1, 1) = RA                                       ! A-A partition of KE
      KE( 2, 2) = R1
      KE( 2, 3) = BETA
      KE( 2, 5) =-C01*BETA
      KE( 2, 6) = C01*R1
      KE( 3, 3) = R2
      KE( 3, 5) =-C01*R2
      KE( 3, 6) = C01*BETA
      KE( 4, 4) = RG
      KE( 5, 5) = K2T
      KE( 5, 6) =-C02*BETA
      KE( 6, 6) = K1T

      KE( 1, 7) =-RA                                       ! A-B partition of KE
      KE( 2, 8) =-R1
      KE( 2, 9) =-BETA
      KE( 2,11) =-C01*BETA
      KE( 2,12) = C01*R1
      KE( 3, 8) =-BETA
      KE( 3, 9) =-R2
      KE( 3,11) =-C01*R2
      KE( 3,12) = C01*BETA
      KE( 4,10) =-RG
      KE( 5, 8) = C01*BETA
      KE( 5, 9) = C01*R2
      KE( 5,11) = K4T
      KE( 5,12) =-C03*BETA
      KE( 6, 8) =-C01*R1
      KE( 6, 9) =-C01*BETA
      KE( 6,11) =-C03*BETA
      KE( 6,12) = K3T

      KE( 7, 7) = RA                                       ! B-B partition of KE
      KE( 8, 8) = R1
      KE( 8, 9) = BETA
      KE( 8,11) = C01*BETA
      KE( 8,12) =-C01*R1
      KE( 9, 9) = R2
      KE( 9,11) = C01*R2
      KE( 9,12) =-C01*BETA
      KE(10,10) = RG
      KE(11,11) = K2T
      KE(11,12) =-C02*BETA
      KE(12,12) = K1T

      DO I=2,12                                            ! Set lower triangular partition of KE using symmetry
         DO J=1,I-1
            KE(I,J) = KE(J,I)
         ENDDO
      ENDDO

! Process Pin Flags. NUM_PFLAG_DOFS is a count of the total number of DOF pin flagged. DOFPIN(i) generated in ELMDAT is an
! integer array of the DOF numbers of the pin flagged DOF. TPIN is generated from the linear KE and then applied to KE here, and
! to any other element matrix (e.g. KED below) that must have the same pin flag condensation.

      NUM_PFLAG_DOFS = 0
      DO I=1,12
         IF (DOFPIN(I) > 0) THEN
            NUM_PFLAG_DOFS = NUM_PFLAG_DOFS + 1
         ENDIF
      ENDDO
      IF (NUM_PFLAG_DOFS /= 0) THEN
         CALL PINFLG_TRANSFORM ( NUM_PFLAG_DOFS, TPIN )
         CALL APPLY_PINFLG ( TPIN, KE(1:12,1:12), 12 )
      ENDIF

      DO I=1,6                                             ! Upper left partition of KE is KAA. Need this for calc'ing PTE, SEi
         DO J=1,6
            KAA(I,J) = KE(I,J)
         ENDDO
      ENDDO

      DO I=1,6                                             ! Upper right partition of KE is KAB. Need this for calculating SEi
         DO J=7,12
            KAB(I,J-6) = KE(I,J)
         ENDDO
      ENDDO

      DO I=7,12                                            ! Lower left partition of KE is KBA. Need this for calculating PTE
         DO J=1,6
            KBA(I-6,J) = KE(I,J)
         ENDDO
      ENDDO

! The following matrices are needed if either OPT(2) or OPT(3) or OPT(6) is called for

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         IF (NTSUB > 0) THEN

            DO J=1,NTSUB
               TBAR        = (DT(1,J) + DT(2,J))/TWO
               TPRIME(1,J) = TBAR - TREF
               TPRIME(2,J) = DT(3,J)
               TPRIME(3,J) = DT(4,J)
               TPRIME(4,J) = DT(5,J)
               TPRIME(5,J) = DT(6,J)
            ENDDO

         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Determine element thermal loads.

!     IF ((OPT(2) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         IF (NTSUB > 0) THEN

            DO I=1,6
               DO J=1,5
                  ABAR(I,J) = ZERO
               ENDDO
            ENDDO

            ABAR(1,1) =  ONE
            ABAR(2,2) =  D1I1*L/SIX
            ABAR(2,3) =  D1I1*L/THREE
            ABAR(2,4) = -D12I2*L/SIX
            ABAR(2,5) = -D12I2*L/THREE
            ABAR(3,2) = -D12I1*L/SIX
            ABAR(3,3) = -D12I1*L/THREE
            ABAR(3,4) =  D2I2*L/SIX
            ABAR(3,5) =  D2I2*L/THREE
            ABAR(5,2) = -D12I1/TWO
            ABAR(5,3) = -D12I1/TWO
            ABAR(5,4) =  D2I2/TWO
            ABAR(5,5) =  D2I2/TWO
            ABAR(6,2) = -D1I1/TWO
            ABAR(6,3) = -D1I1/TWO
            ABAR(6,4) =  D12I2/TWO
            ABAR(6,5) =  D12I2/TWO

            DO I=1,6
               DO J=1,5
                  ABAR(I,J) = -ALPHA*L*ABAR(I,J)
               ENDDO
            ENDDO

            CALL MATMULT_FFF ( KAA, ABAR, 6, 6, 5, BTA )
            CALL MATMULT_FFF ( KBA, ABAR, 6, 6, 5, BTB )

            CALL MATMULT_FFF ( BTA, TPRIME, 6, 5, NTSUB, PTA )
            CALL MATMULT_FFF ( BTB, TPRIME, 6, 5, NTSUB, PTB )
            DO I=1,6
               DO J=1,NTSUB
                  PTE(I,J)   = PTA(I,J)
                  PTE(I+6,J) = PTB(I,J)
               ENDDO
            ENDDO

         ENDIF

!     ENDIF

! **********************************************************************************************************************************
! Calculate SE matrices for stress data recovery.

!     IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         DO I=1,3
            DO J=1,6
               B1(I,J) = ZERO
               B2(I,J) = ZERO
            ENDDO
         ENDDO

         IF (DABS(AREA) > EPS1) THEN                       ! Can't calculate axial stress if the area = 0
            B1(1,1) = -ONE/AREA
         ENDIF

         B1(2,5) = -DELTA12
         B1(2,6) = -DELTA1

         B1(3,5) =  DELTA2
         B1(3,6) =  DELTA12

         B2(1,2) =  DELTA1*L
         B2(1,3) = -DELTA12*L
         B2(1,5) = -DELTA12
         B2(1,6) = -DELTA1

         B2(2,2) = -DELTA12*L
         B2(2,3) =  DELTA2*L
         B2(2,5) =  DELTA2
         B2(2,6) =  DELTA12

         IF (JTOR /= ZERO) THEN                            ! Can't calc torsional stress if torsional constant = 0
            B2(3,4) = -SCOEFF/JTOR
         ENDIF

         CALL MATMULT_FFF ( B1, KAA, 3, 6, 6, S11 )
         CALL MATMULT_FFF ( B1, KAB, 3, 6, 6, S12 )
         CALL MATMULT_FFF ( B2, KAA, 3, 6, 6, S21 )
         CALL MATMULT_FFF ( B2, KAB, 3, 6, 6, S22 )

         DO I=1,3
            DO J=1,6
               SE1(I,J,1) = S11(I,J)
               SE2(I,J,1) = S21(I,J)
            ENDDO
            DO J=7,12
               SE1(I,J,1) = S12(I,J-6)
               SE2(I,J,1) = S22(I,J-6)
            ENDDO
         ENDDO

         IF (NTSUB > 0) THEN                               ! Calculate STEi for thermal stress data recovery

            DO I=1,3
               DO J=1,5
                  BT1(I,J) = ZERO
                  BT2(I,J) = ZERO
               ENDDO
            ENDDO

            CALL MATMULT_FFF ( S11, ABAR, 3, 6, 5, BT1 )
            CALL MATMULT_FFF ( S21, ABAR, 3, 6, 5, BT2 )

            CALL MATMULT_FFF ( BT1, TPRIME, 3, 5, NTSUB, DUM1 )
            CALL MATMULT_FFF ( BT2, TPRIME, 3, 5, NTSUB, DUM2 )
            DO I=1,3
               DO J=1,NTSUB
                  STE1(I,J,1) = DUM1(I,J)
                  STE2(I,J,1) = DUM2(I,J)
               ENDDO
            ENDDO

         ENDIF

!     ENDIF

! **********************************************************************************************************************************
! Calculate linear differential stiffness matrix

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

         CALL ELMDIS

                                                           ! Calc BAR forces
         CALL CALC_ELEM_NODE_FORCES
         M1a = -PEL(6)                                     ! M1a (bending moment, plane 1, end a for BAR) - NASTRAN Maz
         M2a =  PEL(5)                                     ! M2a (bending moment, plane 2, end a for BAR) - NASTRAN May
         M1b = -PEL(6) + PEL(2)*L                          ! M1b (bending moment, plane 1, end b for BAR) - NASTRAN Mbz
         M2b =  PEL(5) + PEL(3)*L                          ! M2b (bending moment, plane 2, end b for BAR) - NASTRAN Mby
         V1  = -PEL(2)                                     ! V1  (plane 1 shear for BAR)                  - NASTRAN Vy
         V2  = -PEL(3)                                     ! V2  (plane 2 shear for BAR)                  - NASTRAN Vz
         Fx  = -PEL(1)                                     ! Fx  (axial force for BAR or ROD)             - NASTRAN Fx

         IF (ART_KED == 'Y') THEN
            KED( 1, 1) = ART_TRAN_KED
            KED( 4, 4) = ART_ROT_KED
            KED( 7, 7) = ART_TRAN_KED
            KED(10,10) = ART_ROT_KED
         ENDIF

         C01 = Fx/L

! Bending moment terms. M1 and M2 are internal moments (EI1*v" and EI2*w") with the same sign at both ends for constant moment. The
! coupling of twist to transverse displ must vanish for rigid body twist with constant moment so the end a and b terms have opposite
! signs. The M/2 terms coupling twist to bending rotation make the end moments semitangential so that they are consistent between
! elements at joints where not all rotations are continuous. See McGuire, Gallagher and Ziemian, Matrix Structural Analysis.

         KED( 2, 2) =  (SIX/FIVE)*Fx/L
         KED( 2, 4) =  M2a/L
         KED( 2, 6) =  Fx/TEN
         KED( 2, 8) = -KED( 2, 2)
         KED( 2,10) = -M2b/L
         KED( 2,12) =  KED( 2, 6)

         KED( 3, 3) =  KED( 2, 2)
         KED( 3, 4) = -M1a/L
         KED( 3, 5) = -KED( 2, 6)
         KED( 3, 9) = -KED( 2, 2)
         KED( 3,10) =  M1b/L
         KED( 3,11) = -KED( 2, 6)

         IF (AREA > ZERO) THEN                             ! Fx = 0 if AREA = 0, so the term is zero
            KED( 4, 4) =  (I1 + I2)*Fx/(AREA*L)            ! Polar moment of area
         ENDIF
         KED( 4, 5) = -V1*L/SIX + M1a/TWO
         KED( 4, 6) = -V2*L/SIX + M2a/TWO
         KED( 4, 8) = -KED( 2, 4)
         KED( 4, 9) = -KED( 3, 4)
         KED( 4,10) = -KED( 4, 4)
         KED( 4,11) =  V1*L/SIX
         KED( 4,12) =  V2*L/SIX

         KED( 5, 5) =  FOUR*Fx*L/(THREE*TEN)
         KED( 5, 9) =  KED( 2, 6)
         KED( 5,10) =  V1*L/SIX
         KED( 5,11) = -Fx*L/(THREE*TEN)

         KED( 6, 6) =  KED( 5, 5)
         KED( 6, 8) = -KED( 2, 6)
         KED( 6,10) =  V2*L/SIX
         KED( 6,12) =  KED( 5,11)

         KED( 8, 8) =  KED( 2, 2)
         KED( 8,10) = -KED( 2,10)
         KED( 8,12) = -KED( 2, 6)

         KED( 9, 9) =  KED( 2, 2)
         KED( 9,10) = -KED( 3,10)
         KED( 9,11) =  KED( 2, 6)

         KED(10,10) =  KED( 4, 4)
         KED(10,11) = -V1*L/SIX - M1b/TWO
         KED(10,12) = -V2*L/SIX - M2b/TWO

         KED(11,11) =  KED( 5, 5)

         KED(12,12) =  KED( 5, 5)

         DO I=2,12
            DO J=1,I-1
               KED(I,J) = KED(J,I)
            ENDDO
         ENDDO

! Apply the same pin flag condensation to KED that was applied to KE

         IF (NUM_PFLAG_DOFS /= 0) THEN
            CALL APPLY_PINFLG ( TPIN, KED(1:12,1:12), 12 )
         ENDIF

         IF (DEBUG(178) >= 1) THEN
            WRITE(F06,2001) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
            DO I=1,12
               WRITE(F06,2101) (KED(I,J),J=1,12)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ENDIF

! If this is DIFFEREN or NLSTATIC, add linear and differential stiffness matrices

      IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN

         DO I=1,ELDOF
            DO J=1,ELDOF
               KE(I,J) = KE(I,J) + KED(I,J)
            ENDDO
         ENDDO

         IF (DEBUG(178) >= 1) THEN
            IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN
               WRITE(F06,2002) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
            ELSE
               WRITE(F06,2003) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
            ENDIF
            DO I=1,12
               WRITE(F06,2101) (KE(I,J),J=1,12)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1918 FORMAT(' *ERROR  1918: NEGATIVE SHEAR FACTOR, K1 OR K2, ON ',A,' ELEMENT ',I8,' IS CAUSING A SINGULAR STIFFNESS')

 2001 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': Differential Stiffness Matrix for ', A,                &
             ' element number ', I8)

 2002 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': Linear + Differential Stiffness Matrix for ', A,       &
             ' element number ', I8)

 2003 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': Linear Stiffness Matrix for ', A,' element number ', I8)

 2101 FORMAT(12(1ES14.6))


! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEBUG_BAR1 (WHAT)

      USE PENTIUM_II_KIND
      USE IOUNT1, ONLY                :  ERR, F06

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: WHAT

! **********************************************************************************************************************************
      IF (WHAT == 1) THEN

         WRITE(F06,*)
         WRITE(F06,1997)
         WRITE(F06,'(A,I8)') 'In subr BAR1 with BAR element ',EID
         WRITE(F06,*) '-------------------------------------'
         WRITE(F06,1998) 'L       = ',L
         WRITE(F06,1998) 'AREA    = ',AREA
         WRITE(F06,1998) 'I1      = ',I1
         WRITE(F06,1998) 'I2      = ',I2
         WRITE(F06,1998) 'I12     = ',I12
         WRITE(F06,1998) 'JTOR    = ',JTOR
         WRITE(F06,1998) 'K1      = ',K1
         WRITE(F06,1998) 'K2      = ',K2
         WRITE(F06,1998) 'E       = ',E
         WRITE(F06,1998) 'G       = ',G
         WRITE(F06,1998) 'SCOEFF  = ',SCOEFF
         WRITE(F06,1998) 'ALPHA   = ',ALPHA
         WRITE(F06,1998) 'TREF    = ',TREF
         WRITE(F06,1998) 'G1      = ',G1
         WRITE(F06,1998) 'G2      = ',G2
         WRITE(F06,1998) 'R0      = ',R0
         WRITE(F06,1999) 'R1      = ',R1, ' initially (before checking I12, G1)'
         WRITE(F06,1999) 'R2      = ',R2, ' initially (before checking I12, G2)'
         WRITE(F06,1998) 'BETA    = ',BETA
         WRITE(F06,1998) 'DEN     = ',DEN

      ELSE IF (WHAT == 2) THEN

         WRITE(F06,1999) 'R1      = ',R1, ' finally'
         WRITE(F06,1999) 'R2      = ',R2, ' finally'
         WRITE(F06,1998) 'DELTA1  = ',DELTA1
         WRITE(F06,1998) 'DELTA2  = ',DELTA2
         WRITE(F06,1998) 'DELTA12 = ',DELTA12
         WRITE(F06,1998) 'RA      = ',RA
         WRITE(F06,1998) 'RG      = ',RG
         WRITE(F06,1998) 'K11     = ',K11
         WRITE(F06,1998) 'K12     = ',K12
         WRITE(F06,1998) 'K21     = ',K21
         WRITE(F06,1998) 'K22     = ',K22
         WRITE(F06,1998) 'K1T     = ',K1T
         WRITE(F06,1998) 'K2T     = ',K2T
         WRITE(F06,1998) 'K3T     = ',K3T
         WRITE(F06,1998) 'K4T     = ',K4T
         WRITE(F06,1998) 'C01     = ',C01
         WRITE(F06,1998) 'C02     = ',C02
         WRITE(F06,1998) 'C03     = ',C03
         WRITE(F06,*)
         WRITE(F06,1997)
         WRITE(F06,*)

      ENDIF

! **********************************************************************************************************************************
 1997 FORMAT('************************************************************')

 1998 FORMAT(A, 1ES14.6)

 1999 FORMAT(A, 1ES14.6, A)

! **********************************************************************************************************************************

      END SUBROUTINE DEBUG_BAR1

      END SUBROUTINE BAR1


      SUBROUTINE BART ( OPT, L, AREA, I1, I2, JTOR, SCOEFF, K1, K2, I12, E, G, ALPHA, TREF, NUM_PFLAG_DOFS, TPIN )

! Calculates, for 1-D Timoshenko beam element

!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linear stiffness matrix      , if OPT(6) = 'N' (i.e. always calc KE linear unless OPT(6) = 'Y')
!  4) KED       = element differen stiff matrix        , if OPT(6) = 'Y'
!  5) NUM_PFLAG_DOFS, TPIN = number of pin flagged DOF's and the pin flag transformation (TPIN is only defined if NUM_PFLAG_DOFS > 0)


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, NTSUB, BLNK_SUB_NAM, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX, TEN, TWELVE
      USE DEBUG_PARAMETERS
      USE PARAMS, ONLY                :  EPSIL, ART_KED, ART_ROT_KED, ART_TRAN_KED
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ELDOF, DOFPIN, DT, EID, NUM_EMG_FATAL_ERRS, KE, KED, PEL, PTE, SE1, SE2, STE1, STE2, TYPE,&
                                         UEL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF
      USE ELEMENT_RECOVERY_SUPPORT, ONLY:  CALC_ELEM_NODE_FORCES

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BART'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices

      INTEGER(LONG)                   :: I,J               ! DO loop induces
      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG), INTENT(OUT)      :: NUM_PFLAG_DOFS    ! The number of pin flagged DOF's for this element


      REAL(DOUBLE) , INTENT(IN)       :: ALPHA             ! Coefficient of thermal expansion
      REAL(DOUBLE) , INTENT(IN)       :: AREA              ! Cross-sectional area
      REAL(DOUBLE) , INTENT(IN)       :: E                 ! Youngs modulus
      REAL(DOUBLE) , INTENT(IN)       :: G                 ! Shear modulus
      REAL(DOUBLE) , INTENT(IN)       :: I1                ! Bending inertia in plane 1
      REAL(DOUBLE) , INTENT(IN)       :: I12               ! Product of inertia
      REAL(DOUBLE) , INTENT(IN)       :: I2                ! Bending inertia in plane 2
      REAL(DOUBLE) , INTENT(IN)       :: JTOR              ! Torsional constant
      REAL(DOUBLE) , INTENT(IN)       :: K1                ! Shear constant for plane 1 (used in K1*AREA*G)
      REAL(DOUBLE) , INTENT(IN)       :: K2                ! Shear constant for plane 2 (used in K1*AREA*G)
      REAL(DOUBLE) , INTENT(IN)       :: L                 ! Elem length
      REAL(DOUBLE) , INTENT(IN)       :: SCOEFF            ! Stress recovery coeff for torsion
      REAL(DOUBLE) , INTENT(IN)       :: TREF              ! Element reference temperature
      REAL(DOUBLE)                    :: ABAR(6,5)         ! Thermal displacements per unit TPRIME
      REAL(DOUBLE)                    :: B1(3,6)           ! Intermediate array used in calculating SE1 stress matrix
      REAL(DOUBLE)                    :: B2(3,6)           ! Intermediate array used in calculating SE2 stress matrix
      REAL(DOUBLE)                    :: BT1(3,5)          ! Intermediate array used in calculating STE1 thermal stress effects
      REAL(DOUBLE)                    :: BT2(3,5)          ! Intermediate array used in calculating STE2 thermal stress effects
      REAL(DOUBLE)                    :: BTA(6,5)          ! = KAA*ABAR. Each col gives thermal loads at end a per unit TPRIME value
      REAL(DOUBLE)                    :: BTB(6,5)          ! = KBA*ABAR. Each col gives thermal loads at end a per unit TPRIME value
      REAL(DOUBLE)                    :: C01               ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: DELTA1            ! I2/DEN (1/I1, or 0 if I1 = 0, when DEN = 0)
      REAL(DOUBLE)                    :: DELTA2            ! I1/DEN (1/I2, or 0 if I2 = 0, when DEN = 0)
      REAL(DOUBLE)                    :: DELTA12           ! I12/DEN (0 when DEN = 0)
      REAL(DOUBLE)                    :: D1I1, D2I2        ! DELTA1*I1, DELTA2*I2 (1 when DEN = 0)
      REAL(DOUBLE)                    :: D12I1, D12I2      ! DELTA12*I1, DELTA12*I2 (0 when DEN = 0)
      REAL(DOUBLE)                    :: DEN               ! I1*I2 - I12**2. Only zero if I12 = 0 and I1 and/or I2 = 0
      REAL(DOUBLE)                    :: DUM1(3,NTSUB)     ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM2(3,NTSUB)     ! Intermediate matrix
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare for real zero
      REAL(DOUBLE)                    :: G1                ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: G2                ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: C1,C2,C3,C4,C5,C6,C7,C8,C9
      REAL(DOUBLE)                    :: KAA(6,6)          ! Upper left  6 x 6 partition of KE
      REAL(DOUBLE) , INTENT(OUT)      :: TPIN(12,12)       ! Pin flag transformation, u = TPIN*u_r (from PINFLG_TRANSFORM).
                                                           ! Returned so the caller can apply it to other elem matrices (e.g. ME)
      REAL(DOUBLE)                    :: KBA(6,6)          ! Lower left  6 x 6 partition of KE
      REAL(DOUBLE)                    :: KAB(6,6)          ! Upper right 6 x 6 partition of KE
      REAL(DOUBLE)                    :: PTA(6,NTSUB)      ! Thermal loads at end a of BAR (rows 1- 6 of PTE)
      REAL(DOUBLE)                    :: PTb(6,NTSUB)      ! Thermal loads at end b of BAR (rows 7-12 of PTE)
      REAL(DOUBLE)                    :: RG                ! Intermediate variable used in calc terms for the stiff matrix, KE
      REAL(DOUBLE)                    :: S11(3,6)          ! Intermediate matrix used in calculating SE1 stress matrix
      REAL(DOUBLE)                    :: S12(3,6)          ! Intermediate matrix used in calculating SE1 stress matrix
      REAL(DOUBLE)                    :: S21(3,6)          ! Intermediate matrix used in calculating SE2 stress matrix
      REAL(DOUBLE)                    :: S22(3,6)          ! Intermediate matrix used in calculating SE2 stress matrix
      REAL(DOUBLE)                    :: TBAR              ! Average elem temperature

! The following are used for the differential stiffness matrix calc. See NASTRAN Prog's Manual (COSMIC 1972) page 4.87.30

      REAL(DOUBLE)                    :: M1a               ! Bend mom, plane 1, end a (also Maz, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: M2a               ! Bend mom, plane 2, end a (also May, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: M1b               ! Bend mom, plane 1, end b (also Mbz, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: M2b               ! Bend mom, plane 2, end b (also Mby, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: V1                ! Plane 1 shear (also Vy, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: V2                ! Plane 2 shear (also Vz, in NASTRAN Prog Man differ stiff)
      REAL(DOUBLE)                    :: Fx                ! Axial force (also Fx, in NASTRAN Prog Man differ stiff)

      INTRINSIC DABS

       REAL(DOUBLE)                    :: TPRIME(5,NTSUB)   ! Matrix where each col has the 5 temperature/gradients for the BAR elem



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Calculate element stiffness matrix KE(12,12), regardless of OPT(4) since the calculation of PTE and/or SEi, STEi
! uses values from KE. Initialize KE here to make sure it starts as zero for the BAR element

      DO I=1,12
         DO J=1,12
            KE(I,J) = ZERO
         ENDDO
      ENDDO

      G1 = K1*AREA*G
      G2 = K2*AREA*G

      C1 =  AREA*E/L

      C2 =  G1/L
      C3 =  G1/2                                           ! NOTE: TESSLER HAS OPPOSITE SIGN BUT NEGATIVE, NOT POSITIVE, WORKS
      C4 =  E*I2/L + G1*L/4
      C5 = -E*I2/L + G1*L/4

      C6 =  G2/L
      C7 =  G2/2
      C8 =  E*I1/L + G2*L/4
      C9 = -E*I1/L + G2*L/4

      RG =  G*JTOR/L

      DEN    = I1*I2 - I12*I12                             ! I1*I2 > I12^2 was checked when PBAR's were read in subr LOADB
      CALL BAR_INERTIA_FACTORS ( I1, I2, ZERO, DELTA1, DELTA2, DELTA12, D1I1, D2I2, D12I1, D12I2 )   ! I12 = 0 for BART

      IF (DEBUG(203) > 0) CALL DEBUG_BART ( 1 )

! Generate KE

      KE( 1, 1) = C1
      KE( 1, 7) =-C1

      KE( 2, 2) = C6
      KE( 2, 6) = C7
      KE( 2, 8) =-C6
      KE( 2,12) = C7

      KE( 3, 3) = C2
      KE( 3, 5) =-C3
      KE( 3, 9) =-C2
      KE( 3,11) =-C3

      KE( 4, 4) = RG
      KE( 4,10) =-RG

      KE( 5, 5) = C4
      KE( 5, 9) = C3
      KE( 5,11) = C5

      KE( 6, 6) = C8
      KE( 6, 8) =-C7
      KE( 6,12) = C9

      KE( 7, 7) = C1

      KE( 8, 8) = C6
      KE( 8,12) =-C7

      KE( 9, 9) = C2
      KE( 9,11) = C3

      KE(10,10) = RG

      KE(11,11) = C4

      KE(12,12) = C8

      DO I=2,12                                            ! Set lower triangular partition of KE using symmetry
         DO J=1,I-1
            KE(I,J) = KE(J,I)
         ENDDO
      ENDDO

! Process Pin Flags. NUM_PFLAG_DOFS is a count of the total number of DOF pin flagged. DOFPIN(i) generated in ELMDAT is an
! integer array of the DOF numbers of the pin flagged DOF. TPIN is generated from the linear KE and then applied to KE here, and
! to any other element matrix (e.g. KED below) that must have the same pin flag condensation.

      NUM_PFLAG_DOFS = 0
      DO I=1,12
         IF (DOFPIN(I) > 0) THEN
            NUM_PFLAG_DOFS = NUM_PFLAG_DOFS + 1
         ENDIF
      ENDDO
      IF (NUM_PFLAG_DOFS /= 0) THEN
         CALL PINFLG_TRANSFORM ( NUM_PFLAG_DOFS, TPIN )
         CALL APPLY_PINFLG ( TPIN, KE(1:12,1:12), 12 )
      ENDIF

      DO I=1,6                                             ! Upper left partition of KE is KAA. Need this for calc'ing PTE, SEi
         DO J=1,6
            KAA(I,J) = KE(I,J)
         ENDDO
      ENDDO

      DO I=1,6                                             ! Upper right partition of KE is KAB. Need this for calculating SEi
         DO J=7,12
            KAB(I,J-6) = KE(I,J)
         ENDDO
      ENDDO

      DO I=7,12                                            ! Lower left partition of KE is KBA. Need this for calculating PTE
         DO J=1,6
            KBA(I-6,J) = KE(I,J)
         ENDDO
      ENDDO

! The following matrices are needed if either OPT(2) or OPT(3) or OPT(6) is called for

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN
         IF (NTSUB > 0) THEN

            DO J=1,NTSUB
               TBAR        = (DT(1,J) + DT(2,J))/TWO
               TPRIME(1,J) = TBAR - TREF
               TPRIME(2,J) = DT(3,J)
               TPRIME(3,J) = DT(4,J)
               TPRIME(4,J) = DT(5,J)
               TPRIME(5,J) = DT(6,J)
            ENDDO

         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Determine element thermal loads.

!     IF ((OPT(2) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         IF (NTSUB > 0) THEN

            DO I=1,6
               DO J=1,5
                  ABAR(I,J) = ZERO
               ENDDO
            ENDDO

            ABAR(1,1) =  ONE
            ABAR(2,2) =  D1I1*L/SIX
            ABAR(2,3) =  D1I1*L/THREE
            ABAR(2,4) = -D12I2*L/SIX
            ABAR(2,5) = -D12I2*L/THREE
            ABAR(3,2) = -D12I1*L/SIX
            ABAR(3,3) = -D12I1*L/THREE
            ABAR(3,4) =  D2I2*L/SIX
            ABAR(3,5) =  D2I2*L/THREE
            ABAR(5,2) = -D12I1/TWO
            ABAR(5,3) = -D12I1/TWO
            ABAR(5,4) =  D2I2/TWO
            ABAR(5,5) =  D2I2/TWO
            ABAR(6,2) = -D1I1/TWO
            ABAR(6,3) = -D1I1/TWO
            ABAR(6,4) =  D12I2/TWO
            ABAR(6,5) =  D12I2/TWO

            DO I=1,6
               DO J=1,5
                  ABAR(I,J) = -ALPHA*L*ABAR(I,J)
               ENDDO
            ENDDO

            CALL MATMULT_FFF ( KAA, ABAR, 6, 6, 5, BTA )
            CALL MATMULT_FFF ( KBA, ABAR, 6, 6, 5, BTB )

            CALL MATMULT_FFF ( BTA, TPRIME, 6, 5, NTSUB, PTA )
            CALL MATMULT_FFF ( BTB, TPRIME, 6, 5, NTSUB, PTB )
            DO I=1,6
               DO J=1,NTSUB
                  PTE(I,J)   = PTA(I,J)
                  PTE(I+6,J) = PTB(I,J)
               ENDDO
            ENDDO

         ENDIF

!     ENDIF

! **********************************************************************************************************************************
! Calculate SE matrices for stress data recovery.

!     IF ((OPT(3) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         DO I=1,3
            DO J=1,6
               B1(I,J) = ZERO
               B2(I,J) = ZERO
            ENDDO
         ENDDO

         IF (DABS(AREA) > EPS1) THEN                       ! Can't calculate axial stress if the area = 0
            B1(1,1) = -ONE/AREA
         ENDIF

         B1(2,5) = -DELTA12
         B1(2,6) = -DELTA1

         B1(3,5) =  DELTA2
         B1(3,6) =  DELTA12

         B2(1,2) =  DELTA1*L
         B2(1,3) = -DELTA12*L
         B2(1,5) = -DELTA12
         B2(1,6) = -DELTA1

         B2(2,2) = -DELTA12*L
         B2(2,3) =  DELTA2*L
         B2(2,5) =  DELTA2
         B2(2,6) =  DELTA12

         IF (JTOR /= ZERO) THEN                            ! Can't calc torsional stress if torsional constant = 0
            B2(3,4) = -SCOEFF/JTOR
         ENDIF

         CALL MATMULT_FFF ( B1, KAA, 3, 6, 6, S11 )
         CALL MATMULT_FFF ( B1, KAB, 3, 6, 6, S12 )
         CALL MATMULT_FFF ( B2, KAA, 3, 6, 6, S21 )
         CALL MATMULT_FFF ( B2, KAB, 3, 6, 6, S22 )

         DO I=1,3
            DO J=1,6
               SE1(I,J,1) = S11(I,J)
               SE2(I,J,1) = S21(I,J)
            ENDDO
            DO J=7,12
               SE1(I,J,1) = S12(I,J-6)
               SE2(I,J,1) = S22(I,J-6)
            ENDDO
         ENDDO

         IF (NTSUB > 0) THEN                               ! Calculate STEi for thermal stress data recovery

            DO I=1,3
               DO J=1,5
                  BT1(I,J) = ZERO
                  BT2(I,J) = ZERO
               ENDDO
            ENDDO

            CALL MATMULT_FFF ( S11, ABAR, 3, 6, 5, BT1 )
            CALL MATMULT_FFF ( S21, ABAR, 3, 6, 5, BT2 )

            CALL MATMULT_FFF ( BT1, TPRIME, 3, 5, NTSUB, DUM1 )
            CALL MATMULT_FFF ( BT2, TPRIME, 3, 5, NTSUB, DUM2 )
            DO I=1,3
               DO J=1,NTSUB
                  STE1(I,J,1) = DUM1(I,J)
                  STE2(I,J,1) = DUM2(I,J)
               ENDDO
            ENDDO

         ENDIF

!     ENDIF

! **********************************************************************************************************************************
! Calculate linear differential stiffness matrix

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

         CALL ELMDIS

                                                           ! Calc BAR forces
         CALL CALC_ELEM_NODE_FORCES
         M1a = -PEL(6)                                     ! M1a (bending moment, plane 1, end a for BAR) - NASTRAN Maz
         M2a =  PEL(5)                                     ! M2a (bending moment, plane 2, end a for BAR) - NASTRAN May
         M1b = -PEL(6) + PEL(2)*L                          ! M1b (bending moment, plane 1, end b for BAR) - NASTRAN Mbz
         M2b =  PEL(5) + PEL(3)*L                          ! M2b (bending moment, plane 2, end b for BAR) - NASTRAN Mby
         V1  = -PEL(2)                                     ! V1  (plane 1 shear for BAR)                  - NASTRAN Vy
         V2  = -PEL(3)                                     ! V2  (plane 2 shear for BAR)                  - NASTRAN Vz
         Fx  = -PEL(1)                                     ! Fx  (axial force for BAR or ROD)             - NASTRAN Fx

         IF (ART_KED == 'Y') THEN
            KED( 1, 1) = ART_TRAN_KED
            KED( 4, 4) = ART_ROT_KED
            KED( 7, 7) = ART_TRAN_KED
            KED(10,10) = ART_ROT_KED
         ENDIF

         C01 = Fx/L

! Bending moment terms. M1 and M2 are internal moments (EI1*v" and EI2*w") with the same sign at both ends for constant moment. The
! coupling of twist to transverse displ must vanish for rigid body twist with constant moment so the end a and b terms have opposite
! signs. The M/2 terms coupling twist to bending rotation make the end moments semitangential so that they are consistent between
! elements at joints where not all rotations are continuous. See McGuire, Gallagher and Ziemian, Matrix Structural Analysis.

         KED( 2, 2) =  (SIX/FIVE)*Fx/L
         KED( 2, 4) =  M2a/L
         KED( 2, 6) =  Fx/TEN
         KED( 2, 8) = -KED( 2, 2)
         KED( 2,10) = -M2b/L
         KED( 2,12) =  KED( 2, 6)

         KED( 3, 3) =  KED( 2, 2)
         KED( 3, 4) = -M1a/L
         KED( 3, 5) = -KED( 2, 6)
         KED( 3, 9) = -KED( 2, 2)
         KED( 3,10) =  M1b/L
         KED( 3,11) = -KED( 2, 6)

         IF (AREA > ZERO) THEN                             ! Fx = 0 if AREA = 0, so the term is zero
            KED( 4, 4) =  (I1 + I2)*Fx/(AREA*L)            ! Polar moment of area
         ENDIF
         KED( 4, 5) = -V1*L/SIX + M1a/TWO
         KED( 4, 6) = -V2*L/SIX + M2a/TWO
         KED( 4, 8) = -KED( 2, 4)
         KED( 4, 9) = -KED( 3, 4)
         KED( 4,10) = -KED( 4, 4)
         KED( 4,11) =  V1*L/SIX
         KED( 4,12) =  V2*L/SIX

         KED( 5, 5) =  FOUR*Fx*L/(THREE*TEN)
         KED( 5, 9) =  KED( 2, 6)
         KED( 5,10) =  V1*L/SIX
         KED( 5,11) = -Fx*L/(THREE*TEN)

         KED( 6, 6) =  KED( 5, 5)
         KED( 6, 8) = -KED( 2, 6)
         KED( 6,10) =  V2*L/SIX
         KED( 6,12) =  KED( 5,11)

         KED( 8, 8) =  KED( 2, 2)
         KED( 8,10) = -KED( 2,10)
         KED( 8,12) = -KED( 2, 6)

         KED( 9, 9) =  KED( 2, 2)
         KED( 9,10) = -KED( 3,10)
         KED( 9,11) =  KED( 2, 6)

         KED(10,10) =  KED( 4, 4)
         KED(10,11) = -V1*L/SIX - M1b/TWO
         KED(10,12) = -V2*L/SIX - M2b/TWO

         KED(11,11) =  KED( 5, 5)

         KED(12,12) =  KED( 5, 5)

         DO I=2,12
            DO J=1,I-1
               KED(I,J) = KED(J,I)
            ENDDO
         ENDDO

! Apply the same pin flag condensation to KED that was applied to KE

         IF (NUM_PFLAG_DOFS /= 0) THEN
            CALL APPLY_PINFLG ( TPIN, KED(1:12,1:12), 12 )
         ENDIF

         IF (DEBUG(178) >= 1) THEN
            WRITE(F06,2001) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
            DO I=1,12
               WRITE(F06,2101) (KED(I,J),J=1,12)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ENDIF

! If this is DIFFEREN or NLSTATIC, add linear and differential stiffness matrices

      IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN

         DO I=1,ELDOF
            DO J=1,ELDOF
               KE(I,J) = KE(I,J) + KED(I,J)
            ENDDO
         ENDDO

         IF (DEBUG(178) >= 1) THEN
            IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN
               WRITE(F06,2002) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
            ELSE
               WRITE(F06,2003) TRIM(SUBR_NAME), OPT(6), LOAD_ISTEP, TYPE, EID
            ENDIF
            DO I=1,12
               WRITE(F06,2101) (KE(I,J),J=1,12)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1918 FORMAT(' *ERROR  1918: NEGATIVE SHEAR FACTOR, K1 OR K2, ON ',A,' ELEMENT ',I8,' IS CAUSING A SINGULAR STIFFNESS')

 2001 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': Differential Stiffness Matrix for ', A,                &
             ' element number ', I8)

 2002 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': Linear + Differential Stiffness Matrix for ', A,       &
             ' element number ', I8)

 2003 FORMAT(' In ', A, ' with OPT(6) = ', A, ', and LOAD_ISTEP = ', I8, ': Linear Stiffness Matrix for ', A,' element number ', I8)

 2101 FORMAT(12(1ES14.6))


! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEBUG_BART (WHAT)

      USE PENTIUM_II_KIND
      USE IOUNT1, ONLY                :  ERR, F06

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: WHAT

! **********************************************************************************************************************************
      IF (WHAT == 1) THEN

         WRITE(F06,*)
         WRITE(F06,1997)
         WRITE(F06,'(A,I8)') 'In subr BART with BAR element ',EID
         WRITE(F06,*) '-------------------------------------'
         WRITE(F06,1998) 'L       = ',L
         WRITE(F06,1998) 'AREA    = ',AREA
         WRITE(F06,1998) 'I1      = ',I1
         WRITE(F06,1998) 'I2      = ',I2
         WRITE(F06,1998) 'I12     = ',I12
         WRITE(F06,1998) 'JTOR    = ',JTOR
         WRITE(F06,1998) 'K1      = ',K1
         WRITE(F06,1998) 'K2      = ',K2
         WRITE(F06,1998) 'E       = ',E
         WRITE(F06,1998) 'G       = ',G
         WRITE(F06,1998) 'SCOEFF  = ',SCOEFF
         WRITE(F06,1998) 'ALPHA   = ',ALPHA
         WRITE(F06,1998) 'TREF    = ',TREF
         WRITE(F06,1998) 'G1      = ',G1
         WRITE(F06,1998) 'G2      = ',G2
         WRITE(F06,1998) 'C1      = ',C1
         WRITE(F06,1998) 'C2      = ',C2
         WRITE(F06,1998) 'C3      = ',C3
         WRITE(F06,1998) 'C4      = ',C4
         WRITE(F06,1998) 'C5      = ',C5
         WRITE(F06,1998) 'C6      = ',C6
         WRITE(F06,1998) 'C7      = ',C7
         WRITE(F06,1998) 'C8      = ',C8
         WRITE(F06,1998) 'C9      = ',C9
         WRITE(F06,1998) 'RG      = ',RG
         WRITE(F06,1998) 'DEN     = ',DEN
         WRITE(F06,1998) 'DELTA1  = ',DELTA1
         WRITE(F06,1998) 'DELTA2  = ',DELTA2
         WRITE(F06,1998) 'DELTA12 = ',DELTA12

      ELSE IF (WHAT == 2) THEN


      ENDIF

! **********************************************************************************************************************************
 1997 FORMAT('************************************************************')

 1998 FORMAT(A, 1ES14.6)

! **********************************************************************************************************************************

      END SUBROUTINE DEBUG_BART

      END SUBROUTINE BART


      SUBROUTINE BEAM

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  NUM_EMG_FATAL_ERRS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BEAM'




! **********************************************************************************************************************************
      FATAL_ERR = FATAL_ERR + 1
      NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
      WRITE(ERR,1)
      WRITE(F06,1)



      RETURN

! **********************************************************************************************************************************
    1 FORMAT(' *ERROR     : CODE NOT WRITTEN FOR BEAM ELEMENT MATRICES AS YET')

! **********************************************************************************************************************************

      END SUBROUTINE BEAM


! ##################################################################################################################################

      SUBROUTINE BAR_INERTIA_FACTORS ( I1, I2, I12, DELTA1, DELTA2, DELTA12, D1I1, D2I2, D12I1, D12I2 )

! Calculates the terms of the inverse of the BAR section inertia matrix used for stress recovery, and their products with I1 and I2
! used for thermal loads. With DEN = I1*I2 - I12^2:
!     DELTA1 = I2/DEN, DELTA2 = I1/DEN, DELTA12 = I12/DEN,
!     D1I1 = DELTA1*I1, D2I2 = DELTA2*I2, D12I1 = DELTA12*I1, D12I2 = DELTA12*I2
! DEN can only be zero when I12 = 0 and I1 and/or I2 = 0, since subr CHECK_BAR_MOIs requires I12^2 < I1*I2 when I12 is not zero.
! A zero I1 or I2 means the BAR has no bending stiffness in that plane, so it carries no bending moment and has no bending stress
! there, and the DELTA for that plane is set to zero. The products are given their I12 = 0 values (D1I1 = D2I2 = 1,
! D12I1 = D12I2 = 0) rather than 0*(1/0). The thermal load and thermal stress terms they multiply are zero for the plane with no
! bending stiffness.

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE

      IMPLICIT NONE

      REAL(DOUBLE) , INTENT(IN)       :: I1                ! Bending inertia in plane 1
      REAL(DOUBLE) , INTENT(IN)       :: I2                ! Bending inertia in plane 2
      REAL(DOUBLE) , INTENT(IN)       :: I12               ! Product of inertia
      REAL(DOUBLE) , INTENT(OUT)      :: DELTA1            ! I2/DEN
      REAL(DOUBLE) , INTENT(OUT)      :: DELTA2            ! I1/DEN
      REAL(DOUBLE) , INTENT(OUT)      :: DELTA12           ! I12/DEN
      REAL(DOUBLE) , INTENT(OUT)      :: D1I1              ! DELTA1*I1
      REAL(DOUBLE) , INTENT(OUT)      :: D2I2              ! DELTA2*I2
      REAL(DOUBLE) , INTENT(OUT)      :: D12I1             ! DELTA12*I1
      REAL(DOUBLE) , INTENT(OUT)      :: D12I2             ! DELTA12*I2
      REAL(DOUBLE)                    :: DEN               ! I1*I2 - I12^2

! **********************************************************************************************************************************
      DEN = I1*I2 - I12*I12

      IF (DEN /= ZERO) THEN

         DELTA1  = I2/DEN
         DELTA2  = I1/DEN
         DELTA12 = I12/DEN
         D1I1    = DELTA1*I1
         D2I2    = DELTA2*I2
         D12I1   = DELTA12*I1
         D12I2   = DELTA12*I2

      ELSE                                                 ! I12 = 0 and I1 and/or I2 = 0

         DELTA1 = ZERO
         IF (I1 /= ZERO) DELTA1 = ONE/I1
         DELTA2 = ZERO
         IF (I2 /= ZERO) DELTA2 = ONE/I2
         DELTA12 = ZERO
         D1I1    = ONE
         D2I2    = ONE
         D12I1   = ZERO
         D12I2   = ZERO

      ENDIF

      RETURN

      END SUBROUTINE BAR_INERTIA_FACTORS

! ##################################################################################################################################

      SUBROUTINE PINFLG_TRANSFORM ( NUM_PFLAG_DOFS, TPIN )

! Generates the transformation, TPIN, that statically condenses the pin flagged DOF's out of the element linear stiffness matrix KE.
! With the elem DOF's partitioned into retained (r) and pin flagged (p) sets:

!     u = TPIN*u_r     where     u_p = -Kpp^-1 * Kpr * u_r

! TPIN is ELDOF x ELDOF: identity for the retained DOF's, rows p hold -Kpp^-1*Kpr and columns p are zero. An element matrix X is
! condensed with X = TPIN'*X*TPIN (subr APPLY_PINFLG), which gives zero rows/cols for the pin flagged DOF's. For KE this is the
! usual static condensation. KE is not modified here. If the pin flags cannot be applied TPIN is returned as identity.

! Kpp is factored with Cholesky. The squared diagonal terms of the factor are the same pivots that would be obtained by condensing
! the pin flagged DOF's out of KE one at a time, so they are checked against EPSIL(1) to detect a pin flagged DOF that has lost all
! of its stiffness due to the other pin flags (e.g. torsion released at both ends).

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  EID, ELDOF, NUM_EMG_FATAL_ERRS, KE, DOFPIN, TYPE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'PINFLG_TRANSFORM'

      INTEGER(LONG), INTENT(IN)       :: NUM_PFLAG_DOFS    ! The number of pin flagged DOF's for this element
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IERROR            ! Count of errors. Error occurs if a diag KE for a pinflaged DOF is zero
      INTEGER(LONG)                   :: INFO              ! Output from LAPACK DPOTRF/DPOTRS
      INTEGER(LONG)                   :: PDOF              ! A DOF component number (1 digit) from array DOFPIN

      LOGICAL                         :: PINNED(ELDOF)     ! .TRUE. for elem DOF's that are pin flagged

      REAL(DOUBLE) , INTENT(OUT)      :: TPIN(ELDOF,ELDOF) ! Transformation from retained to all elem DOF's due to pin flags
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare for real zero
                                                           ! Partition of KE for the pin flagged DOF's (Cholesky factor after DPOTRF)
      REAL(DOUBLE)                    :: KPP(NUM_PFLAG_DOFS,NUM_PFLAG_DOFS)
                                                           ! Rows of KE for pin flagged DOF's with cols of pin flagged DOF's zero.
      REAL(DOUBLE)                    :: KPR(NUM_PFLAG_DOFS,ELDOF)          ! Kpp^-1*Kpr after DPOTRS

      INTRINSIC DABS

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Initialize TPIN to identity (no condensation)

      DO J=1,ELDOF
         DO K=1,ELDOF
            TPIN(J,K) = ZERO
         ENDDO
         TPIN(J,J) = ONE
      ENDDO

! Check to make sure that the diagonal stiffness for the pin flagged DOF's are not zero

      IERROR = 0
      DO I=1,NUM_PFLAG_DOFS
         PDOF = DOFPIN(I)
         IF (DABS(KE(PDOF,PDOF)) <= EPS1) THEN
            IERROR = IERROR + 1
            WRITE(ERR,1921) PDOF, TYPE, EID
            WRITE(F06,1921) PDOF, TYPE, EID
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO
      IF (IERROR > 0) THEN
         RETURN
      ENDIF

! Partition KE into Kpp and Kpr

      DO K=1,ELDOF
         PINNED(K) = .FALSE.
      ENDDO
      DO I=1,NUM_PFLAG_DOFS
         PINNED(DOFPIN(I)) = .TRUE.
      ENDDO

      DO I=1,NUM_PFLAG_DOFS
         DO J=1,NUM_PFLAG_DOFS
            KPP(I,J) = KE(DOFPIN(I),DOFPIN(J))
         ENDDO
         DO K=1,ELDOF
            IF (PINNED(K)) THEN
               KPR(I,K) = ZERO
            ELSE
               KPR(I,K) = KE(DOFPIN(I),K)
            ENDIF
         ENDDO
      ENDDO

! Factor Kpp and check pivots

      CALL DPOTRF ( 'U', NUM_PFLAG_DOFS, KPP, NUM_PFLAG_DOFS, INFO )
      IF (INFO == 0) THEN
         DO I=1,NUM_PFLAG_DOFS
            IF (KPP(I,I)*KPP(I,I) <= EPS1) THEN
               INFO = I
               EXIT
            ENDIF
         ENDDO
      ENDIF
      IF (INFO /= 0) THEN
         IF (INFO > 0) THEN
            PDOF = DOFPIN(INFO)
         ELSE
            PDOF = 0
         ENDIF
         WRITE(ERR,1937) TYPE, EID, PDOF
         WRITE(F06,1937) TYPE, EID, PDOF
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

! Solve Kpp*X = Kpr (DPOTRS only fails for invalid arguments, which cannot occur here) and put -X in rows p of TPIN. The cols of
! X for the pin flagged DOF's are zero, so this also zeros the diagonal terms of TPIN for the pin flagged DOF's

      CALL DPOTRS ( 'U', NUM_PFLAG_DOFS, ELDOF, KPP, NUM_PFLAG_DOFS, KPR, NUM_PFLAG_DOFS, INFO )

      DO I=1,NUM_PFLAG_DOFS
         DO K=1,ELDOF
            TPIN(DOFPIN(I),K) = -KPR(I,K)
         ENDDO
      ENDDO

      RETURN

! **********************************************************************************************************************************
 1921 FORMAT(' *ERROR  1921: DISPL COMPONENT ',I3,' ON ',A,' ELEMENT ',I8,' HAS NO STIFFNESS AND CANNOT BE PIN FLAGGED')

 1937 FORMAT(' *ERROR  1937: DURING THE PROCESS OF APPLYING PIN FLAGS ON ',A,' ELEMENT ',I8,' THE DIAGONAL STIFFNESS FOR DISPL '   &
                    ,/,14X,' COMPONENT ',I2,' HAS BECOME ZERO. CANNOT PIN FLAG THIS COMPONENT')

! **********************************************************************************************************************************

      END SUBROUTINE PINFLG_TRANSFORM

! ##################################################################################################################################

      SUBROUTINE APPLY_PINFLG ( TPIN, X, N )

! Applies the pin flag transformation generated in subr PINFLG_TRANSFORM to element matrix X: X = TPIN'*X*TPIN.
! X is symmetrized afterwards to remove round-off asymmetry from the triple product.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE CONSTANTS_1, ONLY           :  HALF
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'APPLY_PINFLG'

      INTEGER(LONG), INTENT(IN)       :: N                 ! Row/col size of TPIN and X
      INTEGER(LONG)                   :: I,J               ! DO loop indices

      REAL(DOUBLE) , INTENT(IN)       :: TPIN(N,N)         ! Pin flag transformation from subr PINFLG_TRANSFORM
      REAL(DOUBLE) , INTENT(INOUT)    :: X(N,N)            ! Element matrix to transform
      REAL(DOUBLE)                    :: XT(N,N)           ! Intermediate matrix X*TPIN

! **********************************************************************************************************************************
      CALL MATMULT_FFF   ( X, TPIN, N, N, N, XT )
      CALL MATMULT_FFF_T ( TPIN, XT, N, N, N, X )

      DO I=2,N
         DO J=1,I-1
            X(I,J) = HALF*(X(I,J) + X(J,I))
            X(J,I) = X(I,J)
         ENDDO
      ENDDO

      RETURN

      END SUBROUTINE APPLY_PINFLG

   END MODULE LINE_ELEMENTS
