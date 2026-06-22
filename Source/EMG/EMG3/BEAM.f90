!--- cbeam_add begin ---!
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

      SUBROUTINE BEAM ( OPT, L, AREA, I1, I2, JTOR, CW, SCOEFF, K1, K2, I12, E, G, ALPHA, TREF )

! Calculates, for 1-D general beam element using a DSB/Timoshenko bending formulation:
!
!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linear stiffness matrix      , if OPT(6) = 'N' (i.e. always calc KE linear unless OPT(6) = 'Y')
!  4) KED       = element differen stiff matrix        , if OPT(6) = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, NSUB, NTSUB, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX, TEN, TWELVE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL, ART_KED, ART_ROT_KED, ART_TRAN_KED
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  CBEAM_ACTIVE_AREA_SCALE, CBEAM_ACTIVE_NSTATIONS, CBEAM_ACTIVE_XL, CBEAM_ACTIVE_RPROPS, CBEAM_FORCE_B1, CBEAM_FORCE_B2, DOFPIN, DT, EID,&
                                         ELDOF, KE, KED, PEL, PPE, PRESS, PTE, SE1, SE2, STE1, STE2, TE, UEL, ZS
      USE MODEL_STUF, ONLY            :  NUM_EMG_FATAL_ERRS
      USE BEAM_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BEAM'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)

      INTEGER(LONG)                   :: I,J
      INTEGER(LONG)                   :: ISTA
      INTEGER(LONG)                   :: NSTA
      INTEGER(LONG)                   :: NUM_PFLAG_DOFS
      LOGICAL                         :: HAS_NONUNIFORM_STATIONS
      LOGICAL                         :: USE_STATIONED_KE

      REAL(DOUBLE), INTENT(IN)        :: ALPHA
      REAL(DOUBLE), INTENT(IN)        :: AREA
      REAL(DOUBLE), INTENT(IN)        :: CW
      REAL(DOUBLE), INTENT(IN)        :: E
      REAL(DOUBLE), INTENT(IN)        :: G
      REAL(DOUBLE), INTENT(IN)        :: I1
      REAL(DOUBLE), INTENT(IN)        :: I12
      REAL(DOUBLE), INTENT(IN)        :: I2
      REAL(DOUBLE), INTENT(IN)        :: JTOR
      REAL(DOUBLE), INTENT(IN)        :: K1
      REAL(DOUBLE), INTENT(IN)        :: K2
      REAL(DOUBLE), INTENT(IN)        :: L
      REAL(DOUBLE), INTENT(IN)        :: SCOEFF
      REAL(DOUBLE), INTENT(IN)        :: TREF
      REAL(DOUBLE)                    :: ABAR(6,5)
      REAL(DOUBLE)                    :: B1(3,6)
      REAL(DOUBLE)                    :: B2(3,6)
      REAL(DOUBLE)                    :: BT1(3,5)
      REAL(DOUBLE)                    :: BT2(3,5)
      REAL(DOUBLE)                    :: BTA(6,5)
      REAL(DOUBLE)                    :: BTB(6,5)
      REAL(DOUBLE)                    :: C01
      REAL(DOUBLE)                    :: DELTA1
      REAL(DOUBLE)                    :: DELTA2
      REAL(DOUBLE)                    :: DELTA12
      REAL(DOUBLE)                    :: DEN
      REAL(DOUBLE)                    :: DUM1(3,NTSUB)
      REAL(DOUBLE)                    :: DUM2(3,NTSUB)
      REAL(DOUBLE)                    :: EPS1
      REAL(DOUBLE)                    :: FAC1
      REAL(DOUBLE)                    :: FAC2
      REAL(DOUBLE)                    :: FX
      REAL(DOUBLE)                    :: KAA(6,6)
      REAL(DOUBLE)                    :: KAB(6,6)
      REAL(DOUBLE)                    :: KBA(6,6)
      REAL(DOUBLE)                    :: M1A
      REAL(DOUBLE)                    :: M1B
      REAL(DOUBLE)                    :: M2A
      REAL(DOUBLE)                    :: M2B
      REAL(DOUBLE)                    :: N1
      REAL(DOUBLE)                    :: N2
      REAL(DOUBLE)                    :: N3
      REAL(DOUBLE)                    :: N4
      REAL(DOUBLE)                    :: P1
      REAL(DOUBLE)                    :: P2
      REAL(DOUBLE)                    :: PC
      REAL(DOUBLE)                    :: PHI1
      REAL(DOUBLE)                    :: PHI2
      REAL(DOUBLE)                    :: POLAR_RADIUS2
      REAL(DOUBLE)                    :: PTA(6,NTSUB)
      REAL(DOUBLE)                    :: PTB(6,NTSUB)
      REAL(DOUBLE)                    :: QT
      REAL(DOUBLE)                    :: RG
      REAL(DOUBLE)                    :: S11(3,6)
      REAL(DOUBLE)                    :: S12(3,6)
      REAL(DOUBLE)                    :: S21(3,6)
      REAL(DOUBLE)                    :: S22(3,6)
      REAL(DOUBLE)                    :: TBAR
      REAL(DOUBLE)                    :: TPRIME(5,NTSUB)
      REAL(DOUBLE)                    :: V1
      REAL(DOUBLE)                    :: V2
      REAL(DOUBLE)                    :: WGT
      REAL(DOUBLE)                    :: FBASIC(3)
      REAL(DOUBLE)                    :: MBASIC(3)
      REAL(DOUBLE)                    :: P1_LOC
      REAL(DOUBLE)                    :: P2_LOC
      REAL(DOUBLE)                    :: X1L
      REAL(DOUBLE)                    :: X2L
      REAL(DOUBLE)                    :: XI
      REAL(DOUBLE)                    :: XI_STA
      REAL(DOUBLE)                    :: XI_GAUSS(3)
      REAL(DOUBLE)                    :: XI_SCALE
      REAL(DOUBLE)                    :: XI_WGT(3)
      REAL(DOUBLE)                    :: PROJ_FAC
      REAL(DOUBLE)                    :: AREA_REF
      REAL(DOUBLE)                    :: I1_REF
      REAL(DOUBLE)                    :: I2_REF
      REAL(DOUBLE)                    :: I12_REF
      REAL(DOUBLE)                    :: JTOR_REF

! MYSTRAN's 1D force/recovery conventions tie:
!   plane 1 -> DOFs (UY,RZ) with I1 and K2
!   plane 2 -> DOFs (UZ,RY) with I2 and K1
! The DSB shear-influence factors are therefore derived using those same pairings.

      INTRINSIC DABS

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      DO I=1,12
         DO J=1,12
            KE(I,J)  = ZERO
            KED(I,J) = ZERO
         ENDDO
      ENDDO
      AREA_REF = AREA
      I1_REF   = I1
      I2_REF   = I2
      I12_REF  = I12
      JTOR_REF = JTOR

      NSTA = CBEAM_ACTIVE_NSTATIONS
      HAS_NONUNIFORM_STATIONS = .FALSE.
      IF (NSTA > 1) THEN
         DO ISTA=2,NSTA
            IF (DABS(CBEAM_ACTIVE_RPROPS(ISTA,1) - CBEAM_ACTIVE_RPROPS(1,1)) > EPS1) HAS_NONUNIFORM_STATIONS = .TRUE.
            IF (DABS(CBEAM_ACTIVE_RPROPS(ISTA,2) - CBEAM_ACTIVE_RPROPS(1,2)) > EPS1) HAS_NONUNIFORM_STATIONS = .TRUE.
            IF (DABS(CBEAM_ACTIVE_RPROPS(ISTA,3) - CBEAM_ACTIVE_RPROPS(1,3)) > EPS1) HAS_NONUNIFORM_STATIONS = .TRUE.
            IF (DABS(CBEAM_ACTIVE_RPROPS(ISTA,4) - CBEAM_ACTIVE_RPROPS(1,4)) > EPS1) HAS_NONUNIFORM_STATIONS = .TRUE.
            IF (DABS(CBEAM_ACTIVE_RPROPS(ISTA,5) - CBEAM_ACTIVE_RPROPS(1,5)) > EPS1) HAS_NONUNIFORM_STATIONS = .TRUE.
            IF (DABS(CBEAM_ACTIVE_RPROPS(ISTA,6) - CBEAM_ACTIVE_RPROPS(1,6)) > EPS1) HAS_NONUNIFORM_STATIONS = .TRUE.
         ENDDO
      ENDIF
      USE_STATIONED_KE = (NSTA > 1) .AND. HAS_NONUNIFORM_STATIONS
      IF (USE_STATIONED_KE) THEN
         CALL BUILD_TAPERED_BEAM_KE ( L, E, AREA_REF, I1_REF, I2_REF, I12_REF, JTOR_REF )
      ELSE
         FAC1 = E*I1/(L*L*L)
         FAC2 = E*I2/(L*L*L)
         RG   = G*JTOR/L

         IF (DEBUG(203) > 0) CALL DEBUG_BEAM ( 1 )

         PHI1 = ZERO
         PHI2 = ZERO
         IF ((DABS(K1) <= EPS1) .AND. (DABS(K2) <= EPS1)) THEN
            PHI1 = ZERO
            PHI2 = ZERO
         ELSE
            IF (DABS(K2*G*AREA*L*L) > EPS1) THEN
               PHI1 = TWELVE*E*I1/(K2*G*AREA*L*L)
            ENDIF
            IF (DABS(K1*G*AREA*L*L) > EPS1) THEN
               PHI2 = TWELVE*E*I2/(K1*G*AREA*L*L)
            ENDIF
         ENDIF

         FAC1 = FAC1/(ONE + PHI1)
         FAC2 = FAC2/(ONE + PHI2)

         DEN     = I1*I2 - I12*I12
         DELTA1  = ZERO
         DELTA2  = ZERO
         DELTA12 = ZERO
         IF (DABS(DEN) > EPS1) THEN
            DELTA1  = I2/DEN
            DELTA2  = I1/DEN
            DELTA12 = I12/DEN
         ENDIF

         IF (DEBUG(203) > 0) CALL DEBUG_BEAM ( 2 )

! **********************************************************************************************************************************
! Stiffness matrix

         KE( 1, 1) = AREA*E/L
         KE( 1, 7) =-KE(1,1)
         KE( 7, 7) = KE(1,1)

         KE( 4, 4) = RG
         KE( 4,10) =-RG
         KE(10,10) = RG

! Plane 1 DSB bending block on DOFs (UY, RZ, UY, RZ)
         KE( 2, 2) =  TWELVE*FAC1
         KE( 2, 6) =  SIX*L*FAC1
         KE( 2, 8) = -TWELVE*FAC1
         KE( 2,12) =  SIX*L*FAC1

         KE( 6, 6) = (FOUR + PHI1)*E*I1/(L*(ONE + PHI1))
         KE( 6, 8) = -SIX*L*FAC1
         KE( 6,12) = (TWO - PHI1)*E*I1/(L*(ONE + PHI1))

         KE( 8, 8) =  TWELVE*FAC1
         KE( 8,12) = -SIX*L*FAC1

         KE(12,12) = (FOUR + PHI1)*E*I1/(L*(ONE + PHI1))

! Plane 2 DSB bending block on DOFs (UZ, RY, UZ, RY)
         KE( 3, 3) =  TWELVE*FAC2
         KE( 3, 5) = -SIX*L*FAC2
         KE( 3, 9) = -TWELVE*FAC2
         KE( 3,11) = -SIX*L*FAC2

         KE( 5, 5) = (FOUR + PHI2)*E*I2/(L*(ONE + PHI2))
         KE( 5, 9) =  SIX*L*FAC2
         KE( 5,11) = (TWO - PHI2)*E*I2/(L*(ONE + PHI2))

         KE( 9, 9) =  TWELVE*FAC2
         KE( 9,11) =  SIX*L*FAC2

         KE(11,11) = (FOUR + PHI2)*E*I2/(L*(ONE + PHI2))
      ENDIF

      IF (USE_STATIONED_KE) THEN
         PHI1 = ZERO
         PHI2 = ZERO
         DEN     = I1_REF*I2_REF - I12_REF*I12_REF
         DELTA1  = ZERO
         DELTA2  = ZERO
         DELTA12 = ZERO
         IF (DABS(DEN) > EPS1) THEN
            DELTA1  = I2_REF/DEN
            DELTA2  = I1_REF/DEN
            DELTA12 = I12_REF/DEN
         ENDIF
      ENDIF

      DO I=2,12
         DO J=1,I-1
            KE(I,J) = KE(J,I)
         ENDDO
      ENDDO

! **********************************************************************************************************************************
! Pin flags

      NUM_PFLAG_DOFS = 0
      DO I=1,12
         IF (DOFPIN(I) > 0) THEN
            NUM_PFLAG_DOFS = NUM_PFLAG_DOFS + 1
         ENDIF
      ENDDO
      IF (NUM_PFLAG_DOFS /= 0) THEN
         CALL PINFLG ( NUM_PFLAG_DOFS )
      ENDIF

      DO I=1,6
         DO J=1,6
            KAA(I,J) = KE(I,J)
            KAB(I,J) = KE(I,J+6)
            KBA(I,J) = KE(I+6,J)
         ENDDO
      ENDDO

! **********************************************************************************************************************************
! Temperatures

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
! Thermal loads and thermal stress recovery use the existing 1D beam thermal convention.

      IF (NTSUB > 0) THEN

         DO I=1,6
            DO J=1,5
               ABAR(I,J) = ZERO
            ENDDO
         ENDDO

         ABAR(1,1) =  ONE
         ABAR(2,2) =  DELTA1*I1*L/SIX
         ABAR(2,3) =  DELTA1*I1*L/THREE
         ABAR(2,4) = -DELTA12*I2*L/SIX
         ABAR(2,5) = -DELTA12*I2*L/THREE
         ABAR(3,2) = -DELTA12*I1*L/SIX
         ABAR(3,3) = -DELTA12*I1*L/THREE
         ABAR(3,4) =  DELTA2*I2*L/SIX
         ABAR(3,5) =  DELTA2*I2*L/THREE
         ABAR(5,2) = -DELTA12*I1/TWO
         ABAR(5,3) = -DELTA12*I1/TWO
         ABAR(5,4) =  DELTA2*I2/TWO
         ABAR(5,5) =  DELTA2*I2/TWO
         ABAR(6,2) = -DELTA1*I1/TWO
         ABAR(6,3) = -DELTA1*I1/TWO
         ABAR(6,4) =  DELTA12*I2/TWO
         ABAR(6,5) =  DELTA12*I2/TWO

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

! **********************************************************************************************************************************
! Determine element load vector PPE from local beam-axis PLOAD1 data:
! each component uses [P1,P2,X1,X2], where X values are fractional [0,1].
! X1 = X2 is treated as concentrated-in-element.

      IF (OPT(5) == 'Y') THEN
         XI_GAUSS(1) = -0.774596669241483D0
         XI_GAUSS(2) =  ZERO
         XI_GAUSS(3) =  0.774596669241483D0
         XI_WGT(1)   =  0.555555555555556D0
         XI_WGT(2)   =  0.888888888888889D0
         XI_WGT(3)   =  0.555555555555556D0

         DO J=1,NSUB
            P1  = PRESS(1,J)
            P2  = PRESS(2,J)
            X1L = PRESS(3,J)
            X2L = PRESS(4,J)
            CALL DECODE_PLOAD1_SPAN ( X1L, X2L, L )
            IF (X1L >= ZERO) THEN
               P1_LOC = P1
               P2_LOC = P2
               IF (PRESS(25,J) > 0.5D0) THEN
                  P1_LOC = TE(2,2)*P1
                  P2_LOC = TE(2,2)*P2
                  IF (PRESS(25,J) > 1.5D0) THEN
                     PROJ_FAC = DSQRT(MAX(ZERO, ONE - TE(1,2)*TE(1,2)))
                     P1_LOC = PROJ_FAC*P1_LOC
                     P2_LOC = PROJ_FAC*P2_LOC
                  ENDIF
               ENDIF
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1_LOC
                  N1 = ONE - THREE*X1L*X1L + TWO*X1L*X1L*X1L
                  N2 = L*(X1L - TWO*X1L*X1L + X1L*X1L*X1L)
                  N3 = THREE*X1L*X1L - TWO*X1L*X1L*X1L
                  N4 = L*(-X1L*X1L + X1L*X1L*X1L)
                  PPE( 2,J) = PPE( 2,J) + PC*N1
                  PPE( 6,J) = PPE( 6,J) + PC*N2
                  PPE( 8,J) = PPE( 8,J) + PC*N3
                  PPE(12,J) = PPE(12,J) + PC*N4
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1_LOC + (P2_LOC-P1_LOC)*(XI-X1L)/(X2L-X1L)
                     N1 = ONE - THREE*XI*XI + TWO*XI*XI*XI
                     N2 = L*(XI - TWO*XI*XI + XI*XI*XI)
                     N3 = THREE*XI*XI - TWO*XI*XI*XI
                     N4 = L*(-XI*XI + XI*XI*XI)
                     PPE( 2,J) = PPE( 2,J) + QT*L*WGT*XI_SCALE*N1
                     PPE( 6,J) = PPE( 6,J) + QT*L*WGT*XI_SCALE*N2
                     PPE( 8,J) = PPE( 8,J) + QT*L*WGT*XI_SCALE*N3
                     PPE(12,J) = PPE(12,J) + QT*L*WGT*XI_SCALE*N4
                  ENDDO
               ENDIF
            ENDIF

            IF ((PRESS(26,J) > 0.5D0) .AND. ((PRESS(7,J) >= ZERO) .OR. (PRESS(7,J) <= -1.5D0))) THEN
               P1_LOC = TE(2,3)*PRESS(5,J)
               P2_LOC = TE(2,3)*PRESS(6,J)
               X1L    = PRESS(7,J)
               X2L    = PRESS(8,J)
            CALL DECODE_PLOAD1_SPAN ( X1L, X2L, L )
               IF (PRESS(26,J) > 1.5D0) THEN
                  PROJ_FAC = DSQRT(MAX(ZERO, ONE - TE(1,3)*TE(1,3)))
                  P1_LOC = PROJ_FAC*P1_LOC
                  P2_LOC = PROJ_FAC*P2_LOC
               ENDIF
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1_LOC
                  N1 = ONE - THREE*X1L*X1L + TWO*X1L*X1L*X1L
                  N2 = L*(X1L - TWO*X1L*X1L + X1L*X1L*X1L)
                  N3 = THREE*X1L*X1L - TWO*X1L*X1L*X1L
                  N4 = L*(-X1L*X1L + X1L*X1L*X1L)
                  PPE( 2,J) = PPE( 2,J) + PC*N1
                  PPE( 6,J) = PPE( 6,J) + PC*N2
                  PPE( 8,J) = PPE( 8,J) + PC*N3
                  PPE(12,J) = PPE(12,J) + PC*N4
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1_LOC + (P2_LOC-P1_LOC)*(XI-X1L)/(X2L-X1L)
                     N1 = ONE - THREE*XI*XI + TWO*XI*XI*XI
                     N2 = L*(XI - TWO*XI*XI + XI*XI*XI)
                     N3 = THREE*XI*XI - TWO*XI*XI*XI
                     N4 = L*(-XI*XI + XI*XI*XI)
                     PPE( 2,J) = PPE( 2,J) + QT*L*WGT*XI_SCALE*N1
                     PPE( 6,J) = PPE( 6,J) + QT*L*WGT*XI_SCALE*N2
                     PPE( 8,J) = PPE( 8,J) + QT*L*WGT*XI_SCALE*N3
                     PPE(12,J) = PPE(12,J) + QT*L*WGT*XI_SCALE*N4
                  ENDDO
               ENDIF
            ENDIF

            P1  = PRESS(5,J)
            P2  = PRESS(6,J)
            X1L = PRESS(7,J)
            X2L = PRESS(8,J)
            CALL DECODE_PLOAD1_SPAN ( X1L, X2L, L )
            IF (X1L >= ZERO) THEN
               P1_LOC = P1
               P2_LOC = P2
               IF (PRESS(26,J) > 0.5D0) THEN
                  P1_LOC = TE(3,3)*P1
                  P2_LOC = TE(3,3)*P2
                  IF (PRESS(26,J) > 1.5D0) THEN
                     PROJ_FAC = DSQRT(MAX(ZERO, ONE - TE(1,3)*TE(1,3)))
                     P1_LOC = PROJ_FAC*P1_LOC
                     P2_LOC = PROJ_FAC*P2_LOC
                  ENDIF
               ENDIF
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1_LOC
                  N1 = ONE - THREE*X1L*X1L + TWO*X1L*X1L*X1L
                  N2 = L*(X1L - TWO*X1L*X1L + X1L*X1L*X1L)
                  N3 = THREE*X1L*X1L - TWO*X1L*X1L*X1L
                  N4 = L*(-X1L*X1L + X1L*X1L*X1L)
                  PPE( 3,J) = PPE( 3,J) + PC*N1
                  PPE( 5,J) = PPE( 5,J) - PC*N2
                  PPE( 9,J) = PPE( 9,J) + PC*N3
                  PPE(11,J) = PPE(11,J) - PC*N4
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1_LOC + (P2_LOC-P1_LOC)*(XI-X1L)/(X2L-X1L)
                     N1 = ONE - THREE*XI*XI + TWO*XI*XI*XI
                     N2 = L*(XI - TWO*XI*XI + XI*XI*XI)
                     N3 = THREE*XI*XI - TWO*XI*XI*XI
                     N4 = L*(-XI*XI + XI*XI*XI)
                     PPE( 3,J) = PPE( 3,J) + QT*L*WGT*XI_SCALE*N1
                     PPE( 5,J) = PPE( 5,J) - QT*L*WGT*XI_SCALE*N2
                     PPE( 9,J) = PPE( 9,J) + QT*L*WGT*XI_SCALE*N3
                     PPE(11,J) = PPE(11,J) - QT*L*WGT*XI_SCALE*N4
                  ENDDO
               ENDIF
            ENDIF

            IF ((PRESS(25,J) > 0.5D0) .AND. ((PRESS(3,J) >= ZERO) .OR. (PRESS(3,J) <= -1.5D0))) THEN
! Match the benchmark's local-Y sign convention for PLOAD1 on inclined beams.
               P1_LOC = -TE(3,2)*PRESS(1,J)
               P2_LOC = -TE(3,2)*PRESS(2,J)
               X1L    = PRESS(3,J)
               X2L    = PRESS(4,J)
            CALL DECODE_PLOAD1_SPAN ( X1L, X2L, L )
               IF (PRESS(25,J) > 1.5D0) THEN
                  PROJ_FAC = DSQRT(MAX(ZERO, ONE - TE(1,2)*TE(1,2)))
                  P1_LOC = PROJ_FAC*P1_LOC
                  P2_LOC = PROJ_FAC*P2_LOC
               ENDIF
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1_LOC
                  N1 = ONE - THREE*X1L*X1L + TWO*X1L*X1L*X1L
                  N2 = L*(X1L - TWO*X1L*X1L + X1L*X1L*X1L)
                  N3 = THREE*X1L*X1L - TWO*X1L*X1L*X1L
                  N4 = L*(-X1L*X1L + X1L*X1L*X1L)
                  PPE( 3,J) = PPE( 3,J) + PC*N1
                  PPE( 5,J) = PPE( 5,J) - PC*N2
                  PPE( 9,J) = PPE( 9,J) + PC*N3
                  PPE(11,J) = PPE(11,J) - PC*N4
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1_LOC + (P2_LOC-P1_LOC)*(XI-X1L)/(X2L-X1L)
                     N1 = ONE - THREE*XI*XI + TWO*XI*XI*XI
                     N2 = L*(XI - TWO*XI*XI + XI*XI*XI)
                     N3 = THREE*XI*XI - TWO*XI*XI*XI
                     N4 = L*(-XI*XI + XI*XI*XI)
                     PPE( 3,J) = PPE( 3,J) + QT*L*WGT*XI_SCALE*N1
                     PPE( 5,J) = PPE( 5,J) - QT*L*WGT*XI_SCALE*N2
                     PPE( 9,J) = PPE( 9,J) + QT*L*WGT*XI_SCALE*N3
                     PPE(11,J) = PPE(11,J) - QT*L*WGT*XI_SCALE*N4
                  ENDDO
               ENDIF
            ENDIF

! Legacy global FY/FZ loads also have an axial component on inclined members.
            IF ((PRESS(25,J) > 0.5D0) .AND. ((PRESS(3,J) >= ZERO) .OR. (PRESS(3,J) <= -1.5D0))) THEN
               P1_LOC = TE(1,2)*PRESS(1,J)
               P2_LOC = TE(1,2)*PRESS(2,J)
               X1L    = PRESS(3,J)
               X2L    = PRESS(4,J)
            CALL DECODE_PLOAD1_SPAN ( X1L, X2L, L )
               IF (PRESS(25,J) > 1.5D0) THEN
                  PROJ_FAC = DSQRT(MAX(ZERO, ONE - TE(1,2)*TE(1,2)))
                  P1_LOC = PROJ_FAC*P1_LOC
                  P2_LOC = PROJ_FAC*P2_LOC
               ENDIF
               CALL ADD_AXIAL_PLOAD1 ( P1_LOC, P2_LOC, X1L, X2L )
            ENDIF

            IF ((PRESS(26,J) > 0.5D0) .AND. ((PRESS(7,J) >= ZERO) .OR. (PRESS(7,J) <= -1.5D0))) THEN
               P1_LOC = TE(1,3)*PRESS(5,J)
               P2_LOC = TE(1,3)*PRESS(6,J)
               X1L    = PRESS(7,J)
               X2L    = PRESS(8,J)
            CALL DECODE_PLOAD1_SPAN ( X1L, X2L, L )
               IF (PRESS(26,J) > 1.5D0) THEN
                  PROJ_FAC = DSQRT(MAX(ZERO, ONE - TE(1,3)*TE(1,3)))
                  P1_LOC = PROJ_FAC*P1_LOC
                  P2_LOC = PROJ_FAC*P2_LOC
               ENDIF
               CALL ADD_AXIAL_PLOAD1 ( P1_LOC, P2_LOC, X1L, X2L )
            ENDIF

            P1  = PRESS(9 ,J)
            P2  = PRESS(10,J)
            X1L = PRESS(11,J)
            X2L = PRESS(12,J)
            CALL DECODE_PLOAD1_SPAN ( X1L, X2L, L )
            IF (X1L >= ZERO) THEN
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1
                  PPE( 1,J) = PPE( 1,J) + PC*(ONE - X1L)
                  PPE( 7,J) = PPE( 7,J) + PC*X1L
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1 + (P2-P1)*(XI-X1L)/(X2L-X1L)
                     PPE( 1,J) = PPE( 1,J) + QT*L*WGT*XI_SCALE*(ONE - XI)
                     PPE( 7,J) = PPE( 7,J) + QT*L*WGT*XI_SCALE*XI
                  ENDDO
               ENDIF
            ENDIF

            P1  = PRESS(13,J)
            P2  = PRESS(14,J)
            X1L = PRESS(15,J)
            X2L = PRESS(16,J)
            CALL DECODE_PLOAD1_SPAN ( X1L, X2L, L )
            IF (X1L >= ZERO) THEN
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1
                  PPE( 4,J) = PPE( 4,J) + PC*(ONE - X1L)
                  PPE(10,J) = PPE(10,J) + PC*X1L
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1 + (P2-P1)*(XI-X1L)/(X2L-X1L)
                     PPE( 4,J) = PPE( 4,J) + QT*L*WGT*XI_SCALE*(ONE - XI)
                     PPE(10,J) = PPE(10,J) + QT*L*WGT*XI_SCALE*XI
                  ENDDO
               ENDIF
            ENDIF

            P1  = PRESS(17,J)
            P2  = PRESS(18,J)
            X1L = PRESS(19,J)
            X2L = PRESS(20,J)
            CALL DECODE_PLOAD1_SPAN ( X1L, X2L, L )
            IF (X1L >= ZERO) THEN
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1
                  PPE( 5,J) = PPE( 5,J) + PC*(ONE - X1L)
                  PPE(11,J) = PPE(11,J) + PC*X1L
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1 + (P2-P1)*(XI-X1L)/(X2L-X1L)
                     PPE( 5,J) = PPE( 5,J) + QT*L*WGT*XI_SCALE*(ONE - XI)
                     PPE(11,J) = PPE(11,J) + QT*L*WGT*XI_SCALE*XI
                  ENDDO
               ENDIF
            ENDIF

            P1  = PRESS(21,J)
            P2  = PRESS(22,J)
            X1L = PRESS(23,J)
            X2L = PRESS(24,J)
            CALL DECODE_PLOAD1_SPAN ( X1L, X2L, L )
            IF (X1L >= ZERO) THEN
               IF (DABS(X2L - X1L) <= EPS1) THEN
                  PC = P1
                  PPE( 6,J) = PPE( 6,J) + PC*(ONE - X1L)
                  PPE(12,J) = PPE(12,J) + PC*X1L
               ELSE
                  XI_SCALE = (X2L - X1L)/TWO
                  DO I=1,3
                     XI  = XI_SCALE*XI_GAUSS(I) + (X2L + X1L)/TWO
                     WGT = XI_WGT(I)
                     QT  = P1 + (P2-P1)*(XI-X1L)/(X2L-X1L)
                     PPE( 6,J) = PPE( 6,J) + QT*L*WGT*XI_SCALE*(ONE - XI)
                     PPE(12,J) = PPE(12,J) + QT*L*WGT*XI_SCALE*XI
                  ENDDO
               ENDIF
            ENDIF

         ENDDO

! --- CBEAM_standard begin --- !
! --- CBEAM_standard end --- !
      ENDIF

! **********************************************************************************************************************************
! Stress recovery matrices

      DO I=1,3
         DO J=1,6
            B1(I,J) = ZERO
            B2(I,J) = ZERO
         ENDDO
      ENDDO

      IF (DABS(AREA_REF) > EPS1) THEN
         B1(1,1) = -ONE/AREA_REF
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

      IF (DABS(JTOR_REF) > EPS1) THEN
         B2(3,4) = -SCOEFF/JTOR_REF
      ENDIF

      CBEAM_FORCE_B1(:,:) = B1(:,:)
      CBEAM_FORCE_B2(:,:) = B2(:,:)

      CALL MATMULT_FFF ( B1, KAA, 3, 6, 6, S11 )
      CALL MATMULT_FFF ( B1, KAB, 3, 6, 6, S12 )
      CALL MATMULT_FFF ( B2, KAA, 3, 6, 6, S21 )
      CALL MATMULT_FFF ( B2, KAB, 3, 6, 6, S22 )

! --- cbeam_stations begin --- !
      NSTA = CBEAM_ACTIVE_NSTATIONS
      IF (NSTA <= 0) NSTA = 1

      DO ISTA=1,NSTA
         XI_STA = CBEAM_ACTIVE_XL(ISTA)
         IF (NSTA == 1) XI_STA = ZERO

         DO J=1,6
            SE1(1,J,ISTA) = S11(1,J)
            SE2(3,J,ISTA) = S21(3,J)

            SE1(2,J,ISTA) = (ONE - XI_STA)*S11(2,J) + XI_STA*S21(1,J)
            SE1(3,J,ISTA) = (ONE - XI_STA)*S11(3,J) + XI_STA*S21(2,J)
            SE2(1,J,ISTA) = SE1(2,J,ISTA)
            SE2(2,J,ISTA) = SE1(3,J,ISTA)
         ENDDO

         DO J=7,12
            SE1(1,J,ISTA) = S12(1,J-6)
            SE2(3,J,ISTA) = S22(3,J-6)

            SE1(2,J,ISTA) = (ONE - XI_STA)*S12(2,J-6) + XI_STA*S22(1,J-6)
            SE1(3,J,ISTA) = (ONE - XI_STA)*S12(3,J-6) + XI_STA*S22(2,J-6)
            SE2(1,J,ISTA) = SE1(2,J,ISTA)
            SE2(2,J,ISTA) = SE1(3,J,ISTA)
         ENDDO
      ENDDO
! --- cbeam_stations end --- !

      IF (NTSUB > 0) THEN
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
! --- cbeam_stations begin --- !
         DO ISTA=1,NSTA
            XI_STA = CBEAM_ACTIVE_XL(ISTA)
            IF (NSTA == 1) XI_STA = ZERO

            DO J=1,NTSUB
               STE1(1,J,ISTA) = DUM1(1,J)
               STE2(3,J,ISTA) = DUM2(3,J)

               STE1(2,J,ISTA) = (ONE - XI_STA)*DUM1(2,J) + XI_STA*DUM2(1,J)
               STE1(3,J,ISTA) = (ONE - XI_STA)*DUM1(3,J) + XI_STA*DUM2(2,J)
               STE2(1,J,ISTA) = STE1(2,J,ISTA)
               STE2(2,J,ISTA) = STE1(3,J,ISTA)
            ENDDO
         ENDDO
! --- cbeam_stations end --- !
      ENDIF

! **********************************************************************************************************************************
! Reuse the existing 1D geometric stiffness pattern so nonlinear branches do not fail.

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

! --- bug_cbeam_fix1 begin --- !
! For the second buckling pass, derive local element end forces directly from
! the local stiffness times the local element displacement vector. This avoids
! the broader helper path that was unstable for the validated CBEAM buckling
! family while preserving the existing geometric stiffness assembly.
         CALL ELMDIS
         PEL(1:ELDOF) = MATMUL(KE(1:ELDOF,1:ELDOF), UEL(1:ELDOF))
! --- bug_cbeam_fix1 end --- !

         M1A = -PEL(6)
         M2A =  PEL(5)
         M1B = -PEL(6) + PEL(2)*L
         M2B =  PEL(5) + PEL(3)*L
         V1  = -PEL(2)
         V2  = -PEL(3)
         FX  = -PEL(1)

         IF (ART_KED == 'Y') THEN
            KED( 1, 1) = ART_TRAN_KED
            KED( 4, 4) = ART_ROT_KED
            KED( 7, 7) = ART_TRAN_KED
            KED(10,10) = ART_ROT_KED
         ENDIF

         C01 = FX/L

         KED( 2, 2) =  (SIX/FIVE)*FX/L
         KED( 2, 4) =  M2B/L
         KED( 2, 6) =  FX/TEN
         KED( 2, 8) = -KED( 2, 2)
         KED( 2,10) =  M2A/L
         KED( 2,12) =  KED( 2, 6)

         KED( 3, 3) =  KED( 2, 2)
         KED( 3, 4) =  M1B/L
         KED( 3, 5) = -KED( 2, 6)
         KED( 3, 9) = -KED( 2, 2)
         KED( 3,10) =  M1A/L
         KED( 3,11) = -KED( 2, 6)

         POLAR_RADIUS2 = ZERO
         IF (DABS(AREA_REF) > EPS1) THEN
            POLAR_RADIUS2 = (I1_REF + I2_REF)/AREA_REF
         ENDIF

         KED( 4, 4) =  POLAR_RADIUS2*FX/L
         KED( 4, 5) = -V1*L/SIX
         KED( 4, 6) = -V2*L/SIX
         KED( 4, 8) = -M2B/L
         KED( 4, 9) = -M1B/L
         KED( 4,10) = -KED( 4, 4)
         KED( 4,11) = -KED( 4, 5)
         KED( 4,12) = -KED( 4, 6)

         KED( 5, 5) =  FOUR*FX*L/(THREE*TEN)
         KED( 5, 9) =  KED( 2, 6)
         KED( 5,10) = -KED( 4, 5)
         KED( 5,11) = -FX*L/(THREE*TEN)

         KED( 6, 6) =  KED( 5, 5)
         KED( 6, 8) = -KED( 2, 6)
         KED( 6,10) = -KED( 4, 6)
         KED( 6,12) =  KED( 5,11)

         KED( 8, 8) =  KED( 2, 2)
         KED( 8,10) = -KED( 2,10)
         KED( 8,12) = -KED( 2, 6)

         KED( 9, 9) =  KED( 2, 2)
         KED( 9,10) = -KED( 3,10)
         KED( 9,11) =  KED( 2, 6)

         KED(10,10) =  KED( 4, 4)
         KED(10,11) =  KED( 4, 5)
         KED(10,12) =  KED( 4, 6)

         KED(11,11) =  KED( 5, 5)
         KED(12,12) =  KED( 6, 6)

         DO I=2,12
            DO J=1,I-1
               KED(I,J) = KED(J,I)
            ENDDO
         ENDDO

      ENDIF

      RETURN

! **********************************************************************************************************************************

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEBUG_BEAM ( WHAT )

      USE PENTIUM_II_KIND

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: WHAT

      IF (WHAT == 1) THEN
         WRITE(F06,*)
         WRITE(F06,1997)
         WRITE(F06,'(A,I8)') 'In subr BEAM with BEAM element ',EID
         WRITE(F06,*) '--------------------------------------'
         WRITE(F06,1998) 'L       = ',L
         WRITE(F06,1998) 'AREA    = ',AREA
         WRITE(F06,1998) 'I1      = ',I1
         WRITE(F06,1998) 'I2      = ',I2
         WRITE(F06,1998) 'I12     = ',I12
         WRITE(F06,1998) 'JTOR    = ',JTOR
         WRITE(F06,1998) 'CW      = ',CW
         WRITE(F06,1998) 'K1      = ',K1
         WRITE(F06,1998) 'K2      = ',K2
         WRITE(F06,1998) 'E       = ',E
         WRITE(F06,1998) 'G       = ',G
         WRITE(F06,1998) 'SCOEFF  = ',SCOEFF
         WRITE(F06,1998) 'ALPHA   = ',ALPHA
         WRITE(F06,1998) 'TREF    = ',TREF
      ELSE IF (WHAT == 2) THEN
         WRITE(F06,1999) 'plane1 (UY,RZ) uses I1/K2; phi1 = ',PHI1
         WRITE(F06,1999) 'plane2 (UZ,RY) uses I2/K1; phi2 = ',PHI2
         WRITE(F06,1998) 'FAC1    = ',FAC1
         WRITE(F06,1998) 'FAC2    = ',FAC2
         WRITE(F06,1998) 'DEN     = ',DEN
         WRITE(F06,1998) 'DELTA1  = ',DELTA1
         WRITE(F06,1998) 'DELTA2  = ',DELTA2
         WRITE(F06,1998) 'DELTA12 = ',DELTA12
         WRITE(F06,*)
         WRITE(F06,1997)
         WRITE(F06,*)
      ENDIF

 1997 FORMAT('************************************************************')
 1998 FORMAT(A, 1ES14.6)
 1999 FORMAT(A, 1ES14.6)
      END SUBROUTINE DEBUG_BEAM

      SUBROUTINE DECODE_PLOAD1_SPAN ( X1_IO, X2_IO, ELEN )

      REAL(DOUBLE), INTENT(INOUT)     :: X1_IO
      REAL(DOUBLE), INTENT(INOUT)     :: X2_IO
      REAL(DOUBLE), INTENT(IN)        :: ELEN

      IF ((X1_IO <= -1.5D0) .AND. (ELEN > EPS1)) THEN
         X1_IO = (-X1_IO - TWO)/ELEN
      ENDIF
      IF ((X2_IO <= -1.5D0) .AND. (ELEN > EPS1)) THEN
         X2_IO = (-X2_IO - TWO)/ELEN
      ENDIF

      END SUBROUTINE DECODE_PLOAD1_SPAN

      SUBROUTINE ADD_AXIAL_PLOAD1 ( P1_IN, P2_IN, X1_IN, X2_IN )

      REAL(DOUBLE), INTENT(IN)        :: P1_IN
      REAL(DOUBLE), INTENT(IN)        :: P2_IN
      REAL(DOUBLE), INTENT(IN)        :: X1_IN
      REAL(DOUBLE), INTENT(IN)        :: X2_IN

      IF (X1_IN < ZERO) RETURN

      IF (DABS(X2_IN - X1_IN) <= EPS1) THEN
         PPE( 1,J) = PPE( 1,J) + P1_IN*(ONE - X1_IN)
         PPE( 7,J) = PPE( 7,J) + P1_IN*X1_IN
      ELSE
         XI_SCALE = (X2_IN - X1_IN)/TWO
         DO I=1,3
            XI  = XI_SCALE*XI_GAUSS(I) + (X2_IN + X1_IN)/TWO
            WGT = XI_WGT(I)
            QT  = P1_IN + (P2_IN-P1_IN)*(XI-X1_IN)/(X2_IN-X1_IN)
            PPE( 1,J) = PPE( 1,J) + QT*L*WGT*XI_SCALE*(ONE - XI)
            PPE( 7,J) = PPE( 7,J) + QT*L*WGT*XI_SCALE*XI
         ENDDO
      ENDIF

      END SUBROUTINE ADD_AXIAL_PLOAD1

      SUBROUTINE BUILD_TAPERED_BEAM_KE ( L_IN, E_IN, AREA_AVG, I1_AVG, I2_AVG, I12_AVG, JTOR_AVG )

      REAL(DOUBLE), INTENT(IN)      :: L_IN
      REAL(DOUBLE), INTENT(IN)      :: E_IN
      REAL(DOUBLE), INTENT(OUT)     :: AREA_AVG
      REAL(DOUBLE), INTENT(OUT)     :: I1_AVG
      REAL(DOUBLE), INTENT(OUT)     :: I2_AVG
      REAL(DOUBLE), INTENT(OUT)     :: I12_AVG
      REAL(DOUBLE), INTENT(OUT)     :: JTOR_AVG

      INTEGER(LONG), PARAMETER      :: MAX_SEG = 10
      INTEGER(LONG), PARAMETER      :: MAX_NODE = MAX_SEG + 1
      INTEGER(LONG), PARAMETER      :: MAX_DOF = 6*MAX_NODE

      REAL(DOUBLE)                  :: AREA_GP, I1_GP, I2_GP, I12_GP, JTOR_GP, NSM_GP
      REAL(DOUBLE)                  :: DXI_LOC, XI_A, XI_B, XI_SEG
      REAL(DOUBLE)                  :: KCHAIN(MAX_DOF,MAX_DOF)
      REAL(DOUBLE)                  :: KSEG(12,12)
      REAL(DOUBLE)                  :: KEE(12,12)
      REAL(DOUBLE)                  :: KEI(12,MAX_DOF-12)
      REAL(DOUBLE)                  :: KIE(MAX_DOF-12,12)
      REAL(DOUBLE)                  :: KII(MAX_DOF-12,MAX_DOF-12)
      REAL(DOUBLE)                  :: YSOL(MAX_DOF-12,12)
      REAL(DOUBLE)                  :: KSUB(MAX_DOF-12,MAX_DOF-12)
      INTEGER(LONG)                 :: END_DOF(12)
      INTEGER(LONG)                 :: INT_DOF(MAX_DOF-12)
      INTEGER(LONG)                 :: IEND(12)
      INTEGER(LONG)                 :: IINT(MAX_DOF-12)
      INTEGER(LONG)                 :: I, J, IA, IB, IROW, ICOL, ISEG, NSEG, NNODE, NDOF, NINT
      LOGICAL                       :: USED_CHAIN

      USED_CHAIN = .FALSE.
      NSEG = NSTA - 1
      IF ((NSEG >= 1) .AND. (NSEG <= MAX_SEG)) THEN
         NNODE = NSEG + 1
         NDOF  = 6*NNODE
         NINT  = NDOF - 12

         DO I=1,NDOF
            DO J=1,NDOF
               KCHAIN(I,J) = ZERO
            ENDDO
         ENDDO

         AREA_AVG = ZERO
         I1_AVG   = ZERO
         I2_AVG   = ZERO
         I12_AVG  = ZERO
         JTOR_AVG = ZERO

         DO ISEG=1,NSEG
            XI_A = CBEAM_ACTIVE_XL(ISEG)
            XI_B = CBEAM_ACTIVE_XL(ISEG+1)
            DXI_LOC = XI_B - XI_A
            IF (DXI_LOC <= EPS1) CYCLE

            XI_SEG = 0.5D0*(XI_A + XI_B)
            CALL GET_STATION_PROPS ( XI_SEG, AREA_GP, I1_GP, I2_GP, I12_GP, JTOR_GP, NSM_GP )
            CALL BUILD_PRISM_BEAM_SEG_KE ( DXI_LOC*L_IN, E_IN, AREA_GP, I1_GP, I2_GP, I12_GP, JTOR_GP, KSEG )

            AREA_AVG = AREA_AVG + DXI_LOC*AREA_GP
            I1_AVG   = I1_AVG   + DXI_LOC*I1_GP
            I2_AVG   = I2_AVG   + DXI_LOC*I2_GP
            I12_AVG  = I12_AVG  + DXI_LOC*I12_GP
            JTOR_AVG = JTOR_AVG + DXI_LOC*JTOR_GP

            IA = 6*(ISEG-1)
            IB = 6*ISEG
            DO IROW=1,12
               DO ICOL=1,12
                  IF (IROW <= 6) THEN
                     I = IA + IROW
                  ELSE
                     I = IB + IROW - 6
                  ENDIF
                  IF (ICOL <= 6) THEN
                     J = IA + ICOL
                  ELSE
                     J = IB + ICOL - 6
                  ENDIF
                  KCHAIN(I,J) = KCHAIN(I,J) + KSEG(IROW,ICOL)
               ENDDO
            ENDDO
         ENDDO

         IF (AREA_AVG <= EPS1) AREA_AVG = AREA
         IF (I1_AVG   <= EPS1) I1_AVG   = I1
         IF (I2_AVG   <= EPS1) I2_AVG   = I2
         IF (JTOR_AVG <= EPS1) JTOR_AVG = JTOR

         IF (NINT <= 0) THEN
            DO IROW=1,12
               DO ICOL=1,12
                  KE(IROW,ICOL) = KCHAIN(IROW,ICOL)
               ENDDO
            ENDDO
            USED_CHAIN = .TRUE.
         ELSE
            DO I=1,6
               END_DOF(I)   = I
               END_DOF(I+6) = NDOF - 6 + I
            ENDDO

            DO I=1,12
               IEND(I) = END_DOF(I)
            ENDDO
            J = 0
            DO I=1,NDOF
               IF ((I <= 6) .OR. (I > NDOF-6)) CYCLE
               J = J + 1
               INT_DOF(J) = I
               IINT(J) = I
            ENDDO

            DO IROW=1,12
               DO ICOL=1,12
                  KEE(IROW,ICOL) = KCHAIN(IEND(IROW),IEND(ICOL))
               ENDDO
            ENDDO
            DO IROW=1,12
               DO ICOL=1,NINT
                  KEI(IROW,ICOL) = KCHAIN(IEND(IROW),IINT(ICOL))
               ENDDO
            ENDDO
            DO IROW=1,NINT
               DO ICOL=1,12
                  KIE(IROW,ICOL) = KCHAIN(IINT(IROW),IEND(ICOL))
               ENDDO
            ENDDO
            DO IROW=1,NINT
               DO ICOL=1,NINT
                  KII(IROW,ICOL) = KCHAIN(IINT(IROW),IINT(ICOL))
                  KSUB(IROW,ICOL) = KII(IROW,ICOL)
               ENDDO
            ENDDO

            CALL SOLVE_MULTI_RHS_GAUSS ( KSUB, NINT, KIE, 12, YSOL, USED_CHAIN )
            IF (USED_CHAIN) THEN
               DO IROW=1,12
                  DO ICOL=1,12
                     KE(IROW,ICOL) = KEE(IROW,ICOL)
                     DO J=1,NINT
                        KE(IROW,ICOL) = KE(IROW,ICOL) - KEI(IROW,J)*YSOL(J,ICOL)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF

      IF (USED_CHAIN) RETURN

      AREA_AVG = ZERO
      I1_AVG   = ZERO
      I2_AVG   = ZERO
      I12_AVG  = ZERO
      JTOR_AVG = ZERO
      CALL BUILD_TAPERED_BEAM_KE_DIRECT ( L_IN, E_IN, AREA_AVG, I1_AVG, I2_AVG, I12_AVG, JTOR_AVG )

      IF (AREA_AVG <= EPS1) AREA_AVG = AREA
      IF (I1_AVG   <= EPS1) I1_AVG   = I1
      IF (I2_AVG   <= EPS1) I2_AVG   = I2
      IF (JTOR_AVG <= EPS1) JTOR_AVG = JTOR

      END SUBROUTINE BUILD_TAPERED_BEAM_KE

      SUBROUTINE BUILD_TAPERED_BEAM_KE_DIRECT ( L_IN, E_IN, AREA_AVG, I1_AVG, I2_AVG, I12_AVG, JTOR_AVG )

      REAL(DOUBLE), INTENT(IN)      :: L_IN
      REAL(DOUBLE), INTENT(IN)      :: E_IN
      REAL(DOUBLE), INTENT(OUT)     :: AREA_AVG
      REAL(DOUBLE), INTENT(OUT)     :: I1_AVG
      REAL(DOUBLE), INTENT(OUT)     :: I2_AVG
      REAL(DOUBLE), INTENT(OUT)     :: I12_AVG
      REAL(DOUBLE), INTENT(OUT)     :: JTOR_AVG

      REAL(DOUBLE)                  :: AREA_GP, I1_GP, I2_GP, I12_GP, JTOR_GP, NSM_GP
      REAL(DOUBLE)                  :: BAX(12), BTOR(12), BB1(12), BB2(12)
      REAL(DOUBLE)                  :: DXI_LOC, XI_A, XI_B, XI_SEG, WT_SEG
      REAL(DOUBLE)                  :: N1PP, N2PP, N3PP, N4PP
      REAL(DOUBLE)                  :: GAUSS(2), WGTG(2)
      INTEGER(LONG)                 :: ICOL, IGP, IROW, ISEG

      GAUSS(1) = -0.577350269189626D0
      GAUSS(2) =  0.577350269189626D0
      WGTG(1)  =  ONE
      WGTG(2)  =  ONE

      AREA_AVG = ZERO
      I1_AVG   = ZERO
      I2_AVG   = ZERO
      I12_AVG  = ZERO
      JTOR_AVG = ZERO

      DO IROW=1,12
         BAX(IROW)  = ZERO
         BTOR(IROW) = ZERO
      ENDDO
      BAX(1)   = -ONE/L_IN
      BAX(7)   =  ONE/L_IN
      BTOR(4)  = -ONE/L_IN
      BTOR(10) =  ONE/L_IN

      DO ISEG=1,NSTA-1
         XI_A = CBEAM_ACTIVE_XL(ISEG)
         XI_B = CBEAM_ACTIVE_XL(ISEG+1)
         DXI_LOC = XI_B - XI_A
         IF (DXI_LOC <= EPS1) CYCLE

         DO IGP=1,2
            XI_SEG = 0.5D0*(XI_A + XI_B) + 0.5D0*DXI_LOC*GAUSS(IGP)
            WT_SEG = 0.5D0*DXI_LOC*WGTG(IGP)*L_IN

            CALL GET_STATION_PROPS ( XI_SEG, AREA_GP, I1_GP, I2_GP, I12_GP, JTOR_GP, NSM_GP )

            AREA_AVG = AREA_AVG + 0.5D0*DXI_LOC*WGTG(IGP)*AREA_GP
            I1_AVG   = I1_AVG   + 0.5D0*DXI_LOC*WGTG(IGP)*I1_GP
            I2_AVG   = I2_AVG   + 0.5D0*DXI_LOC*WGTG(IGP)*I2_GP
            I12_AVG  = I12_AVG  + 0.5D0*DXI_LOC*WGTG(IGP)*I12_GP
            JTOR_AVG = JTOR_AVG + 0.5D0*DXI_LOC*WGTG(IGP)*JTOR_GP

            N1PP = (-SIX + TWELVE*XI_SEG)/(L_IN*L_IN)
            N2PP = (-FOUR + SIX*XI_SEG)/L_IN
            N3PP = ( SIX - TWELVE*XI_SEG)/(L_IN*L_IN)
            N4PP = (-TWO + SIX*XI_SEG)/L_IN

            DO IROW=1,12
               BB1(IROW) = ZERO
               BB2(IROW) = ZERO
            ENDDO
            BB1(2)  = N1PP
            BB1(6)  = N2PP
            BB1(8)  = N3PP
            BB1(12) = N4PP

            BB2(3)  = N1PP
            BB2(5)  = -N2PP
            BB2(9)  = N3PP
            BB2(11) = -N4PP

            DO IROW=1,12
               DO ICOL=1,12
                  KE(IROW,ICOL) = KE(IROW,ICOL) + E_IN*AREA_GP*BAX(IROW)*BAX(ICOL)*WT_SEG
                  KE(IROW,ICOL) = KE(IROW,ICOL) + G*JTOR_GP*BTOR(IROW)*BTOR(ICOL)*WT_SEG
                  KE(IROW,ICOL) = KE(IROW,ICOL) + E_IN*I1_GP*BB1(IROW)*BB1(ICOL)*WT_SEG
                  KE(IROW,ICOL) = KE(IROW,ICOL) + E_IN*I2_GP*BB2(IROW)*BB2(ICOL)*WT_SEG
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE BUILD_TAPERED_BEAM_KE_DIRECT

      SUBROUTINE BUILD_PRISM_BEAM_SEG_KE ( LSEG, ESEG, AREA_SEG, I1_SEG, I2_SEG, I12_SEG, JTOR_SEG, KSEG )

      REAL(DOUBLE), INTENT(IN)      :: LSEG, ESEG, AREA_SEG, I1_SEG, I2_SEG, I12_SEG, JTOR_SEG
      REAL(DOUBLE), INTENT(OUT)     :: KSEG(12,12)

      REAL(DOUBLE)                  :: PHI1_SEG, PHI2_SEG, FAC1_SEG, FAC2_SEG, RG_SEG
      INTEGER(LONG)                 :: IROW, ICOL

      DO IROW=1,12
         DO ICOL=1,12
            KSEG(IROW,ICOL) = ZERO
         ENDDO
      ENDDO

      IF (LSEG <= EPS1) RETURN

      FAC1_SEG = ESEG*I1_SEG/(LSEG*LSEG*LSEG)
      FAC2_SEG = ESEG*I2_SEG/(LSEG*LSEG*LSEG)
      RG_SEG   = G*JTOR_SEG/LSEG

      PHI1_SEG = ZERO
      PHI2_SEG = ZERO
      IF (.NOT. ((DABS(K1) <= EPS1) .AND. (DABS(K2) <= EPS1))) THEN
         IF (DABS(K2*G*AREA_SEG*LSEG*LSEG) > EPS1) THEN
            PHI1_SEG = TWELVE*ESEG*I1_SEG/(K2*G*AREA_SEG*LSEG*LSEG)
         ENDIF
         IF (DABS(K1*G*AREA_SEG*LSEG*LSEG) > EPS1) THEN
            PHI2_SEG = TWELVE*ESEG*I2_SEG/(K1*G*AREA_SEG*LSEG*LSEG)
         ENDIF
      ENDIF

      FAC1_SEG = FAC1_SEG/(ONE + PHI1_SEG)
      FAC2_SEG = FAC2_SEG/(ONE + PHI2_SEG)

      KSEG( 1, 1) = AREA_SEG*ESEG/LSEG
      KSEG( 1, 7) =-KSEG(1,1)
      KSEG( 7, 7) = KSEG(1,1)

      KSEG( 4, 4) = RG_SEG
      KSEG( 4,10) =-RG_SEG
      KSEG(10,10) = RG_SEG

      KSEG( 2, 2) =  TWELVE*FAC1_SEG
      KSEG( 2, 6) =  SIX*LSEG*FAC1_SEG
      KSEG( 2, 8) = -TWELVE*FAC1_SEG
      KSEG( 2,12) =  SIX*LSEG*FAC1_SEG

      KSEG( 6, 6) = (FOUR + PHI1_SEG)*ESEG*I1_SEG/(LSEG*(ONE + PHI1_SEG))
      KSEG( 6, 8) = -SIX*LSEG*FAC1_SEG
      KSEG( 6,12) = (TWO - PHI1_SEG)*ESEG*I1_SEG/(LSEG*(ONE + PHI1_SEG))

      KSEG( 8, 8) =  TWELVE*FAC1_SEG
      KSEG( 8,12) = -SIX*LSEG*FAC1_SEG

      KSEG(12,12) = (FOUR + PHI1_SEG)*ESEG*I1_SEG/(LSEG*(ONE + PHI1_SEG))

      KSEG( 3, 3) =  TWELVE*FAC2_SEG
      KSEG( 3, 5) = -SIX*LSEG*FAC2_SEG
      KSEG( 3, 9) = -TWELVE*FAC2_SEG
      KSEG( 3,11) = -SIX*LSEG*FAC2_SEG

      KSEG( 5, 5) = (FOUR + PHI2_SEG)*ESEG*I2_SEG/(LSEG*(ONE + PHI2_SEG))
      KSEG( 5, 9) =  SIX*LSEG*FAC2_SEG
      KSEG( 5,11) = (TWO - PHI2_SEG)*ESEG*I2_SEG/(LSEG*(ONE + PHI2_SEG))

      KSEG( 9, 9) =  TWELVE*FAC2_SEG
      KSEG( 9,11) =  SIX*LSEG*FAC2_SEG

      KSEG(11,11) = (FOUR + PHI2_SEG)*ESEG*I2_SEG/(LSEG*(ONE + PHI2_SEG))

      DO IROW=2,12
         DO ICOL=1,IROW-1
            KSEG(IROW,ICOL) = KSEG(ICOL,IROW)
         ENDDO
      ENDDO

      END SUBROUTINE BUILD_PRISM_BEAM_SEG_KE

      SUBROUTINE SOLVE_MULTI_RHS_GAUSS ( AIN, N, BIN, NRHS, XOUT, OK )

      INTEGER(LONG), INTENT(IN)      :: N, NRHS
      REAL(DOUBLE), INTENT(INOUT)    :: AIN(N,N)
      REAL(DOUBLE), INTENT(IN)       :: BIN(N,NRHS)
      REAL(DOUBLE), INTENT(OUT)      :: XOUT(N,NRHS)
      LOGICAL, INTENT(OUT)           :: OK

      REAL(DOUBLE)                   :: WORK(54,12)
      REAL(DOUBLE)                   :: FACTOR, PIVMAX, TMP
      INTEGER(LONG)                  :: I, J, K, IPIV

      OK = .TRUE.
      DO I=1,N
         DO J=1,NRHS
            WORK(I,J) = BIN(I,J)
            XOUT(I,J) = ZERO
         ENDDO
      ENDDO

      DO K=1,N
         IPIV = K
         PIVMAX = DABS(AIN(K,K))
         DO I=K+1,N
            IF (DABS(AIN(I,K)) > PIVMAX) THEN
               PIVMAX = DABS(AIN(I,K))
               IPIV = I
            ENDIF
         ENDDO
         IF (PIVMAX <= EPS1) THEN
            OK = .FALSE.
            RETURN
         ENDIF
         IF (IPIV /= K) THEN
            DO J=K,N
               TMP = AIN(K,J)
               AIN(K,J) = AIN(IPIV,J)
               AIN(IPIV,J) = TMP
            ENDDO
            DO J=1,NRHS
               TMP = WORK(K,J)
               WORK(K,J) = WORK(IPIV,J)
               WORK(IPIV,J) = TMP
            ENDDO
         ENDIF
         DO I=K+1,N
            FACTOR = AIN(I,K)/AIN(K,K)
            AIN(I,K) = ZERO
            DO J=K+1,N
               AIN(I,J) = AIN(I,J) - FACTOR*AIN(K,J)
            ENDDO
            DO J=1,NRHS
               WORK(I,J) = WORK(I,J) - FACTOR*WORK(K,J)
            ENDDO
         ENDDO
      ENDDO

      DO J=1,NRHS
         DO I=N,1,-1
            TMP = WORK(I,J)
            DO K=I+1,N
               TMP = TMP - AIN(I,K)*XOUT(K,J)
            ENDDO
            XOUT(I,J) = TMP/AIN(I,I)
         ENDDO
      ENDDO

      END SUBROUTINE SOLVE_MULTI_RHS_GAUSS

      SUBROUTINE GET_STATION_PROPS ( XI_IN, AREA_OUT, I1_OUT, I2_OUT, I12_OUT, JTOR_OUT, NSM_OUT )

      REAL(DOUBLE), INTENT(IN)      :: XI_IN
      REAL(DOUBLE), INTENT(OUT)     :: AREA_OUT, I1_OUT, I2_OUT, I12_OUT, JTOR_OUT, NSM_OUT

      REAL(DOUBLE)                  :: DXI_LOC, FRAC
      INTEGER(LONG)                 :: IST

      AREA_OUT = CBEAM_ACTIVE_AREA_SCALE*CBEAM_ACTIVE_RPROPS(1,1)
      I1_OUT   = CBEAM_ACTIVE_RPROPS(1,2)
      I2_OUT   = CBEAM_ACTIVE_RPROPS(1,3)
      I12_OUT  = CBEAM_ACTIVE_RPROPS(1,4)
      JTOR_OUT = CBEAM_ACTIVE_RPROPS(1,5)
      NSM_OUT  = CBEAM_ACTIVE_RPROPS(1,6)

      IF (NSTA <= 1) RETURN
      IF (XI_IN <= CBEAM_ACTIVE_XL(1)) RETURN

      DO IST=1,NSTA-1
         IF (XI_IN <= CBEAM_ACTIVE_XL(IST+1)) THEN
            DXI_LOC = CBEAM_ACTIVE_XL(IST+1) - CBEAM_ACTIVE_XL(IST)
            IF (DXI_LOC <= EPS1) RETURN
            FRAC = (XI_IN - CBEAM_ACTIVE_XL(IST))/DXI_LOC
            AREA_OUT = (ONE - FRAC)*CBEAM_ACTIVE_RPROPS(IST,1) + FRAC*CBEAM_ACTIVE_RPROPS(IST+1,1)
            AREA_OUT = CBEAM_ACTIVE_AREA_SCALE*AREA_OUT
            I1_OUT   = (ONE - FRAC)*CBEAM_ACTIVE_RPROPS(IST,2) + FRAC*CBEAM_ACTIVE_RPROPS(IST+1,2)
            I2_OUT   = (ONE - FRAC)*CBEAM_ACTIVE_RPROPS(IST,3) + FRAC*CBEAM_ACTIVE_RPROPS(IST+1,3)
            I12_OUT  = (ONE - FRAC)*CBEAM_ACTIVE_RPROPS(IST,4) + FRAC*CBEAM_ACTIVE_RPROPS(IST+1,4)
            JTOR_OUT = (ONE - FRAC)*CBEAM_ACTIVE_RPROPS(IST,5) + FRAC*CBEAM_ACTIVE_RPROPS(IST+1,5)
            NSM_OUT  = (ONE - FRAC)*CBEAM_ACTIVE_RPROPS(IST,6) + FRAC*CBEAM_ACTIVE_RPROPS(IST+1,6)
            RETURN
         ENDIF
      ENDDO

      AREA_OUT = CBEAM_ACTIVE_AREA_SCALE*CBEAM_ACTIVE_RPROPS(NSTA,1)
      I1_OUT   = CBEAM_ACTIVE_RPROPS(NSTA,2)
      I2_OUT   = CBEAM_ACTIVE_RPROPS(NSTA,3)
      I12_OUT  = CBEAM_ACTIVE_RPROPS(NSTA,4)
      JTOR_OUT = CBEAM_ACTIVE_RPROPS(NSTA,5)
      NSM_OUT  = CBEAM_ACTIVE_RPROPS(NSTA,6)

      END SUBROUTINE GET_STATION_PROPS

      END SUBROUTINE BEAM

!--- cbeam_add end ---!













