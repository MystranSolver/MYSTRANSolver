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
      USE CONSTANTS_1, ONLY           :  HALF, ONE, TWO, ZERO
      USE PARAMS, ONLY                :  EPSIL
      USE DEBUG_PARAMETERS
      USE MODEL_STUF, ONLY            :  EID, ELEM_LEN_AB, EMAT, NUM_EMG_FATAL_ERRS, EPROP, FCONV, ME, ULT_STRE, ULT_STRN, &
                                         TYPE, ZS
! --- cbeam_add begin --- !
      USE PARAMS, ONLY                :  BEAMAMO, BEAMAMO_PID, BEAMAMO_VAL, BEAMM1MO, BEAMM1MO_PID, BEAMM1MO_VAL,          &
                                         BEAMM2MO, BEAMM2MO_PID, BEAMM2MO_VAL, BEAMTMO, BEAMTMO_PID, BEAMTMO_VAL,          &
                                         BEAMV1MO, BEAMV1MO_PID, BEAMV1MO_VAL, BEAMV2MO, BEAMV2MO_PID, BEAMV2MO_VAL,       &
                                         CBEAMAREA, CBEAMAREA_PID, CBEAMAREA_VAL, CBEAMSHR, CBEAMSHR_PID,                  &
                                         CBEAMSHR_VAL, EPSIL, MBEAMAMO_PID, MBEAMM1MO_PID, MBEAMM2MO_PID, MBEAMTMO_PID,    &
                                         MBEAMV1MO_PID, MBEAMV2MO_PID, NCBEAMAREA_PID, NCBEAMSHR_PID, NBEAMAMO_PID,        &
                                         NBEAMM1MO_PID, NBEAMM2MO_PID, NBEAMTMO_PID, NBEAMV1MO_PID, NBEAMV2MO_PID
      USE MODEL_STUF, ONLY            :  CBEAM_ACTIVE_AREA_SCALE, CBEAM_ACTIVE_NSTATIONS, CBEAM_ACTIVE_RPROPS,             &
                                         CBEAM_ACTIVE_XL, INTL_PID, PBAR, PBEAM
! --- cbeam_add end --- !

      USE BREL1_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BREL1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not

      REAL(DOUBLE)                    :: ALPHA             ! Coefficient of thermal expansion
      REAL(DOUBLE)                    :: AREA              ! Cross-sectional area
      REAL(DOUBLE)                    :: E                 ! Youngs modulus
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare for real zero
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
! --- cbeam_add begin --- !
      REAL(DOUBLE)                    :: AREA_EFF          ! Area x Stiffness Modifier
      REAL(DOUBLE)                    :: CW                ! Warping coefficient
      REAL(DOUBLE)                    :: DXI
      REAL(DOUBLE)                    :: K1_EFF
      REAL(DOUBLE)                    :: K2_EFF
      REAL(DOUBLE)                    :: AMOD
      REAL(DOUBLE)                    :: V1MOD
      REAL(DOUBLE)                    :: V2MOD
      REAL(DOUBLE)                    :: M1MOD
      REAL(DOUBLE)                    :: M2MOD
      REAL(DOUBLE)                    :: TMOD
      REAL(DOUBLE)                    :: XI1
      REAL(DOUBLE)                    :: XI2
      INTEGER(LONG)                   :: ISTA
      INTEGER(LONG)                   :: NSTA
      INTEGER(LONG)                   :: PID_EXT
! --- cbeam_add end --- !


! **********************************************************************************************************************************
      EPS1 = EPSIL(1)
! --- cbeam_add begin --- !
      NSM  = ZERO
      CW   = ZERO
      AREA_EFF = ZERO
      K1_EFF = ZERO
      K2_EFF = ZERO
      AMOD  = ONE
      V1MOD = ONE
      V2MOD = ONE
      M1MOD = ONE
      M2MOD = ONE
      TMOD  = ONE
! --- cbeam_add end --- !

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
! --- cbeam_add begin --- !
         PID_EXT  = PBAR(INTL_PID,1)
         AMOD     = GET_BEAMAMO_FOR_PID(PID_EXT)
         V1MOD    = GET_BEAMV1MO_FOR_PID(PID_EXT)
         V2MOD    = GET_BEAMV2MO_FOR_PID(PID_EXT)
         M1MOD    = GET_BEAMM1MO_FOR_PID(PID_EXT)
         M2MOD    = GET_BEAMM2MO_FOR_PID(PID_EXT)
         TMOD     = GET_BEAMTMO_FOR_PID(PID_EXT)
! --- cbeam_add end --- !
         K1       = EPROP(14)                              ! Plane 1 shear factor
         K2       = EPROP(15)                              ! Plane 2 shear factor
! --- cbeam_add begin --- !
         AREA_EFF = AMOD*AREA
         K1_EFF   = V1MOD*K1
         K2_EFF   = V2MOD*K2
         I1       = M1MOD*I1
         I2       = M2MOD*I2
         JTOR     = TMOD*JTOR
         I12      = EPROP(16)                              ! Product of inertia
         ZS(9)    = EPROP(17)                              ! Torsional stress recovery coefficient
         FCONV(1) = AREA
! --- cbeam_add end --- !
      ELSE IF (TYPE == 'BEAM    ') THEN
!! --- cbeam_add begin --- !
         AREA     = EPROP( 1)
         I1       = EPROP( 2)
         I2       = EPROP( 3)
         I12      = EPROP( 4)
         JTOR     = EPROP( 5)
         NSM      = EPROP( 6)
         ZS(1)    = EPROP( 7)
         ZS(2)    = EPROP( 8)
         ZS(3)    = EPROP( 9)
         ZS(4)    = EPROP(10)
         ZS(5)    = EPROP(11)
         ZS(6)    = EPROP(12)
         ZS(7)    = EPROP(13)
         ZS(8)    = EPROP(14)
         PID_EXT  = PBEAM(INTL_PID,1)
         AMOD     = GET_BEAMAMO_FOR_PID(PID_EXT)
         V1MOD    = GET_BEAMV1MO_FOR_PID(PID_EXT)
         V2MOD    = GET_BEAMV2MO_FOR_PID(PID_EXT)
         M1MOD    = GET_BEAMM1MO_FOR_PID(PID_EXT)
         M2MOD    = GET_BEAMM2MO_FOR_PID(PID_EXT)
         TMOD     = GET_BEAMTMO_FOR_PID(PID_EXT)
         K1       = EPROP(30)
         K2       = EPROP(31)
         AREA_EFF = AMOD*AREA
         K1_EFF   = V1MOD*K1
         K2_EFF   = V2MOD*K2
         I1       = M1MOD*I1
         I2       = M2MOD*I2
         JTOR     = TMOD*JTOR
         CBEAM_ACTIVE_AREA_SCALE = GET_CBEAMAREA_FOR_PID(PID_EXT)
         CW       = (EPROP(36) + EPROP(37))/TWO
         ZS(9)    = ZERO
         FCONV(1) = AREA
         NSTA = CBEAM_ACTIVE_NSTATIONS
         IF (NSTA > 1) THEN
            AREA = ZERO
            I1   = ZERO
            I2   = ZERO
            I12  = ZERO
            JTOR = ZERO
            NSM  = ZERO
            DO ISTA=1,NSTA-1
               XI1 = CBEAM_ACTIVE_XL(ISTA)
               XI2 = CBEAM_ACTIVE_XL(ISTA+1)
               DXI = XI2 - XI1
               IF (DXI > EPS1) THEN
                  AREA = AREA + HALF*DXI*(CBEAM_ACTIVE_RPROPS(ISTA,1) + CBEAM_ACTIVE_RPROPS(ISTA+1,1))
                  I1   = I1   + HALF*DXI*(CBEAM_ACTIVE_RPROPS(ISTA,2) + CBEAM_ACTIVE_RPROPS(ISTA+1,2))
                  I2   = I2   + HALF*DXI*(CBEAM_ACTIVE_RPROPS(ISTA,3) + CBEAM_ACTIVE_RPROPS(ISTA+1,3))
                  I12  = I12  + HALF*DXI*(CBEAM_ACTIVE_RPROPS(ISTA,4) + CBEAM_ACTIVE_RPROPS(ISTA+1,4))
                  JTOR = JTOR + HALF*DXI*(CBEAM_ACTIVE_RPROPS(ISTA,5) + CBEAM_ACTIVE_RPROPS(ISTA+1,5))
                  NSM  = NSM  + HALF*DXI*(CBEAM_ACTIVE_RPROPS(ISTA,6) + CBEAM_ACTIVE_RPROPS(ISTA+1,6))
               ENDIF
            ENDDO
            IF (AREA <= EPS1) AREA = EPROP(1)
            IF (I1   <= EPS1) I1   = EPROP(2)
            IF (I2   <= EPS1) I2   = EPROP(3)
            IF (DABS(I12) <= EPS1) I12 = EPROP(4)
            IF (JTOR <= EPS1) JTOR = EPROP(5)
            AREA_EFF = CBEAM_ACTIVE_AREA_SCALE*AREA
            FCONV(1) = AREA
         ENDIF
! --- cbeam_add end --- !
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
! Call routines to calc element matrices (stiffness, etc.)

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(5) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         IF      (TYPE == 'ROD     ') THEN

            CALL ROD1 ( OPT, ELEM_LEN_AB, AREA, JTOR, ZS(1), E, G, ALPHA, TREF )

         ELSE IF (TYPE == 'BAR     ') THEN                 ! Bernoulli-Euler prismatic beam

            IF (DEBUG(249) == 0) THEN
! added area_eff, k1 eff, k2 eff
               CALL BAR1 ( OPT, ELEM_LEN_AB, AREA_EFF, I1, I2, JTOR, ZS(9), K1_EFF, K2_EFF, I12, E, G, ALPHA, TREF )

            ELSE
               IF (DABS(I12) < EPS1) THEN
! added area_eff, k1 eff, k2 eff
               CALL BART ( OPT, ELEM_LEN_AB, AREA_EFF, I1, I2, JTOR, ZS(9), K1_EFF, K2_EFF, I12, E, G, ALPHA, TREF )
               ELSE
                  WRITE(ERR,'(A,I8,A)') ' *ERROR  1962: TIMOSHENKO BAR ELEMENT ',EID,' CANNOT HAVE NONZERO I12. IT WILL BE SET TO I12 = 0.'
                  WRITE(F06,'(A,I8,A)') ' *ERROR  1962: TIMOSHENKO BAR ELEMENT ',EID,' CANNOT HAVE NONZERO I12. IT WILL BE SET TO I12 = 0.'
                  RETURN
               ENDIF

            ENDIF

         ELSE IF (TYPE == 'BEAM    ') THEN                 ! General beam

! --- cbeam_add begin --- !
            CALL BEAM ( OPT, ELEM_LEN_AB, AREA_EFF, I1, I2, JTOR, CW, ZS(9), K1_EFF, K2_EFF, I12, E, G, ALPHA, TREF )
! --- cbeam_end begin --- !

         ENDIF

      ENDIF

      RETURN

! **********************************************************************************************************************************

      CONTAINS
! --- cbeam_add begin --- !
      REAL(DOUBLE) FUNCTION GET_CBEAMSHR_FOR_PID ( PID_IN )

      INTEGER(LONG), INTENT(IN)       :: PID_IN
      INTEGER(LONG)                   :: ILOC

      GET_CBEAMSHR_FOR_PID = CBEAMSHR
      DO ILOC=1,NCBEAMSHR_PID
         IF (CBEAMSHR_PID(ILOC) == PID_IN) THEN
            GET_CBEAMSHR_FOR_PID = CBEAMSHR_VAL(ILOC)
            EXIT
         ENDIF
      ENDDO

      END FUNCTION GET_CBEAMSHR_FOR_PID

! **********************************************************************************************************************************
! cbeam_add shear area
      REAL(DOUBLE) FUNCTION GET_CBEAMAREA_FOR_PID ( PID_IN )

      INTEGER(LONG), INTENT(IN)       :: PID_IN
      INTEGER(LONG)                   :: ILOC

      GET_CBEAMAREA_FOR_PID = CBEAMAREA
      DO ILOC=1,NCBEAMAREA_PID
         IF (CBEAMAREA_PID(ILOC) == PID_IN) THEN
            GET_CBEAMAREA_FOR_PID = CBEAMAREA_VAL(ILOC)
            EXIT
         ENDIF
      ENDDO

      END FUNCTION GET_CBEAMAREA_FOR_PID

! **********************************************************************************************************************************
! cbeam_add axial stiffness modifier
      REAL(DOUBLE) FUNCTION GET_BEAMAMO_FOR_PID ( PID_IN )

      INTEGER(LONG), INTENT(IN)       :: PID_IN
      INTEGER(LONG)                   :: ILOC

      GET_BEAMAMO_FOR_PID = BEAMAMO
      DO ILOC=1,NBEAMAMO_PID
         IF (BEAMAMO_PID(ILOC) == PID_IN) THEN
            GET_BEAMAMO_FOR_PID = BEAMAMO_VAL(ILOC)
            EXIT
         ENDIF
      ENDDO

      END FUNCTION GET_BEAMAMO_FOR_PID

! **********************************************************************************************************************************
! cbeam_add shear minor stiffness modifier
      REAL(DOUBLE) FUNCTION GET_BEAMV1MO_FOR_PID ( PID_IN )

      INTEGER(LONG), INTENT(IN)       :: PID_IN
      INTEGER(LONG)                   :: ILOC

      GET_BEAMV1MO_FOR_PID = BEAMV1MO
      DO ILOC=1,NBEAMV1MO_PID
         IF (BEAMV1MO_PID(ILOC) == PID_IN) THEN
            GET_BEAMV1MO_FOR_PID = BEAMV1MO_VAL(ILOC)
            EXIT
         ENDIF
      ENDDO

      END FUNCTION GET_BEAMV1MO_FOR_PID

! **********************************************************************************************************************************
! cbeam_add shear minor tiffness modifier
      REAL(DOUBLE) FUNCTION GET_BEAMV2MO_FOR_PID ( PID_IN )

      INTEGER(LONG), INTENT(IN)       :: PID_IN
      INTEGER(LONG)                   :: ILOC

      GET_BEAMV2MO_FOR_PID = BEAMV2MO
      DO ILOC=1,NBEAMV2MO_PID
         IF (BEAMV2MO_PID(ILOC) == PID_IN) THEN
            GET_BEAMV2MO_FOR_PID = BEAMV2MO_VAL(ILOC)
            EXIT
         ENDIF
      ENDDO

      END FUNCTION GET_BEAMV2MO_FOR_PID

! **********************************************************************************************************************************
! cbeam_add minor moment stiffness modifier
      REAL(DOUBLE) FUNCTION GET_BEAMM1MO_FOR_PID ( PID_IN )

      INTEGER(LONG), INTENT(IN)       :: PID_IN
      INTEGER(LONG)                   :: ILOC

      GET_BEAMM1MO_FOR_PID = BEAMM1MO
      DO ILOC=1,NBEAMM1MO_PID
         IF (BEAMM1MO_PID(ILOC) == PID_IN) THEN
            GET_BEAMM1MO_FOR_PID = BEAMM1MO_VAL(ILOC)
            EXIT
         ENDIF
      ENDDO

      END FUNCTION GET_BEAMM1MO_FOR_PID

! **********************************************************************************************************************************
! cbeam_add minor moment stiffness modifier
      REAL(DOUBLE) FUNCTION GET_BEAMM2MO_FOR_PID ( PID_IN )

      INTEGER(LONG), INTENT(IN)       :: PID_IN
      INTEGER(LONG)                   :: ILOC

      GET_BEAMM2MO_FOR_PID = BEAMM2MO
      DO ILOC=1,NBEAMM2MO_PID
         IF (BEAMM2MO_PID(ILOC) == PID_IN) THEN
            GET_BEAMM2MO_FOR_PID = BEAMM2MO_VAL(ILOC)
            EXIT
         ENDIF
      ENDDO

      END FUNCTION GET_BEAMM2MO_FOR_PID

! **********************************************************************************************************************************
! cbeam_add torque stiffness modifier
      REAL(DOUBLE) FUNCTION GET_BEAMTMO_FOR_PID ( PID_IN )

      INTEGER(LONG), INTENT(IN)       :: PID_IN
      INTEGER(LONG)                   :: ILOC

      GET_BEAMTMO_FOR_PID = BEAMTMO
      DO ILOC=1,NBEAMTMO_PID
         IF (BEAMTMO_PID(ILOC) == PID_IN) THEN
            GET_BEAMTMO_FOR_PID = BEAMTMO_VAL(ILOC)
            EXIT
         ENDIF
      ENDDO

      END FUNCTION GET_BEAMTMO_FOR_PID
! --- cbeam_add end --- !
      END SUBROUTINE BREL1
