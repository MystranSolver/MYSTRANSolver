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

   MODULE MATERIAL_PROPERTIES

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MATERIAL_PROPS_2D, MATERIAL_PROPS_3D

   CONTAINS

      SUBROUTINE MATERIAL_PROPS_2D ( WRITE_WARN )

! Calculates stress/strain material matrices for isotropic and orthotropic plane stress shell elements with 3 different material
! properties (membrane, bending, transverse shear, bending/membrane coupling)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MEFE, MEMATC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  EPSIL, QUAD4TYP
      USE MODEL_STUF, ONLY            :  ALPVEC, EID, EMG_IFE, EMG_RFE, ERR_SUB_NAM, EB, EBM, EM, ET, NUM_EMG_FATAL_ERRS, EMAT,    &
                                         MTRL_TYPE, QUAD_DELTA, RHO, ULT_STRE, ULT_STRN, THETAM, TREF, TYPE
      use debug_parameters

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATERIAL_PROPS_2D'
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not
      CHARACTER(16*BYTE)              :: MESSAG            ! Part of message printed if material type is invalid

      INTEGER(LONG)                   :: IERROR        = 0 ! Local error indicator meaning some calcs cannot be done
      INTEGER(LONG)                   :: PROG_ERR      = 0 ! Coding error indicator for invalid material type
      INTEGER(LONG)                   :: I,j               ! DO loop index


      REAL(DOUBLE)                    :: ALPHA             ! Isotropic coefficient of thermal expansion
      REAL(DOUBLE)                    :: ALPHA1            ! Orthotropic/Anisotropic coeff of thermal expansion in direction 1
      REAL(DOUBLE)                    :: ALPHA2            ! Orthotropic/Anisotropic coeff of thermal expansion in direction 2
      REAL(DOUBLE)                    :: ALPHA3            ! Anisotropic coeff of thermal expansion in direction 3
      REAL(DOUBLE)                    :: DEN1              ! An intermaediate variable in calculating outputs
      REAL(DOUBLE)                    :: E                 ! Isotropic Young's modulus
      REAL(DOUBLE)                    :: E1                ! Orthotropic Young's modulus in the 1 direction
      REAL(DOUBLE)                    :: E2                ! Orthotropic Young's modulus in the 2 direction
      REAL(DOUBLE)                    :: E0                ! An intermaediate variable in calculating outputs
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare for real zero
      REAL(DOUBLE)                    :: G                 ! Shear modulus
      REAL(DOUBLE)                    :: G12               ! Orthotropic in-plane shear modulus
      REAL(DOUBLE)                    :: G1Z               ! Orthotropic transverse shear modulus for shear in the 1Z plane
      REAL(DOUBLE)                    :: G2Z               ! Orthotropic transverse shear modulus for shear in the 2Z plane
      REAL(DOUBLE)                    :: NU                ! Poisson's ratio
      REAL(DOUBLE)                    :: NU12              ! Orthotropic Poisson's ratio 12
      REAL(DOUBLE)                    :: NU21              ! Orthotropic Poisson's ratio 21
      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1   = EPSIL(1)

! **********************************************************************************************************************************
! Initialize arrays
      ! Material membrane and bending matrices
      EM = ZERO
      EB = ZERO
      ! Material transverse shear matrix
      ET = ZERO
      ! Vector of CTE's
      ALPVEC = ZERO

! **********************************************************************************************************************************
! In-plane stress material props.

mem:  IF (MTRL_TYPE(1) /= 0) THEN

         IF (MTRL_TYPE(1) == 1) THEN                       ! Isotropic properties
!                                                            --------------------
            E             = EMAT( 1,1)
            G             = EMAT( 2,1)
            NU            = EMAT( 3,1)
            RHO(1)        = EMAT( 4,1)
            ALPHA         = EMAT( 5,1)
            TREF(1)       = EMAT( 6,1)

            ULT_STRE(1,1) = EMAT( 8,1)                     ! Dir 1 tension  allowable
            ULT_STRE(2,1) = EMAT( 9,1)                     ! Dir 1 compr    allowable
            ULT_STRE(3,1) = ULT_STRE(1,1)                  ! Dir 2 tension  allowable
            ULT_STRE(4,1) = ULT_STRE(2,1)                  ! Dir 2 compr    allowable
            ULT_STRE(7,1) = EMAT(10,1)                     ! In-plane shear allowable

            ULT_STRN(1,1) = ZERO
            ULT_STRN(2,1) = ZERO
            ULT_STRN(3,1) = ZERO
            ULT_STRN(4,1) = ZERO
            ULT_STRN(7,1) = ZERO

            IF (DABS(E) > EPS1) THEN
               ULT_STRN(1,1) = ULT_STRE(1,1)/E
               ULT_STRN(2,1) = ULT_STRE(2,1)/E
               ULT_STRN(3,1) = ULT_STRE(3,1)/E
               ULT_STRN(4,1) = ULT_STRE(4,1)/E
            ENDIF

            IF (DABS(G) > EPS1) THEN
               ULT_STRN(7,1) = ULT_STRE(7,1)/G
            ENDIF

            DEN1 = ONE - NU*NU
            IF (DABS(DEN1) > EPS1) THEN                    ! If 1-NU^2 = 0, error

               E0 = E/DEN1
               EM(1,1) = E0
               EM(2,2) = EM(1,1)
               EM(3,3) = G
               EM(1,2) = E0*NU
               EM(2,1) = EM(1,2)

               ALPVEC(1,1) = ALPHA
               ALPVEC(2,1) = ALPHA

            ELSE

               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               IERROR = IERROR + 1
               IF (WRT_ERR > 0) THEN
                  WRITE(ERR,1924) DEN1, EID, TYPE
                  WRITE(F06,1924) DEN1, EID, TYPE
               ELSE
                  ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
                  EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1924
                  EMG_RFE(NUM_EMG_FATAL_ERRS,1)   = DEN1
               ENDIF

            ENDIF

         ELSE IF (MTRL_TYPE(1) == 2) THEN                  ! Anisotropic properties
!                                                            ----------------------
            EM(1,1)       = EMAT( 1,1)
            EM(1,2)       = EMAT( 2,1)
            EM(1,3)       = EMAT( 3,1)

            EM(2,1)       = EM(1,2)
            EM(2,2)       = EMAT( 4,1)
            EM(2,3)       = EMAT( 5,1)

            EM(3,1)       = EM(1,3)
            EM(3,2)       = EM(2,3)
            EM(3,3)       = EMAT( 6,1)

            RHO(1)        = EMAT( 7,1)
            ALPHA1        = EMAT( 8,1)
            ALPHA2        = EMAT( 9,1)
            ALPHA3        = EMAT(10,1)
            TREF(1)       = EMAT(11,1)

            ULT_STRE(1,1) = EMAT(13,1)                     ! Dir 1 tension  allowable
            ULT_STRE(2,1) = EMAT(14,1)                     ! Dir 1 compr    allowable
            ULT_STRE(3,1) = ULT_STRE(1,1)                  ! Dir 2 tension  allowable
            ULT_STRE(4,1) = ULT_STRE(2,1)                  ! Dir 2 compr    allowable
            ULT_STRE(7,1) = EMAT(15,1)                     ! In-plane shear allowable

            ULT_STRN(1,1) = ZERO
            ULT_STRN(2,1) = ZERO
            ULT_STRN(3,1) = ZERO
            ULT_STRN(4,1) = ZERO
            ULT_STRN(7,1) = ZERO

            IF (DABS(EM(1,1)) > EPS1) THEN
               ULT_STRN(1,1) = ULT_STRE(1,1)/EM(1,1)
               ULT_STRN(2,1) = ULT_STRE(2,1)/EM(1,1)
            ENDIF

            IF (DABS(EM(2,2)) > EPS1) THEN
               ULT_STRN(3,1) = ULT_STRE(3,1)/EM(2,2)
               ULT_STRN(4,1) = ULT_STRE(4,1)/EM(2,2)
            ENDIF

            IF (DABS(EM(3,3)) > EPS1) THEN
               ULT_STRN(7,1) = ULT_STRE(7,1)/EM(3,3)
            ENDIF

            ALPVEC(1,1)   = ALPHA1
            ALPVEC(2,1)   = ALPHA2
            ALPVEC(3,1)   = ALPHA3

         ELSE IF (MTRL_TYPE(1) == 8) THEN                  ! Orthotropic properties
!                                                            ----------------------                                                           ! Remove the following l.c. code when MIN4T prob w/ ortho mat'l is fixed
            E1            = EMAT( 1,1)
            E2            = EMAT( 2,1)
            NU12          = EMAT( 3,1)
            NU21          = NU12*(E2/E1)                   ! E1, E2 /= 0 checked in BD_MAT8 subr
            G12           = EMAT( 4,1)
            G1Z           = EMAT( 5,1)
            G2Z           = EMAT( 6,1)
            RHO(1)        = EMAT( 7,1)
            ALPHA1        = EMAT( 8,1)
            ALPHA2        = EMAT( 9,1)
            TREF(1)       = EMAT(10,1)

            IF (EMAT(18,1) == 0) THEN                         ! This indicates allowables on the MAT8 entry were stress allowables

               ULT_STRE(1,1) = EMAT(11,1)                     ! Dir 1 tension  allowable
               ULT_STRE(2,1) = EMAT(12,1)                     ! Dir 1 compr    allowable
               ULT_STRE(3,1) = EMAT(13,1)                     ! Dir 2 tension  allowable
               ULT_STRE(4,1) = EMAT(14,1)                     ! Dir 2 compr    allowable
               ULT_STRE(7,1) = EMAT(15,1)                     ! In-plane shear allowable

               ULT_STRN(1,1) = ZERO
               ULT_STRN(2,1) = ZERO
               ULT_STRN(3,1) = ZERO
               ULT_STRN(4,1) = ZERO
               ULT_STRN(7,1) = ZERO

               IF (DABS(E1) > EPS1) THEN
                  ULT_STRN(1,1) = ULT_STRE(1,1)/E1
               ENDIF

               IF (DABS(E2) > EPS1) THEN
                  ULT_STRN(2,1) = ULT_STRE(2,1)/E2
               ENDIF

               IF (DABS(G12) > EPS1) THEN
                  ULT_STRN(3,1) = ULT_STRE(3,1)/G12
               ENDIF

               IF (DABS(G1Z) > EPS1) THEN
                  ULT_STRN(4,1) = ULT_STRE(3,1)/G1Z
               ENDIF

               IF (DABS(G2Z) > EPS1) THEN
                  ULT_STRN(7,1) = ULT_STRE(7,1)/G2Z
               ENDIF

            ELSE                                           ! Allowables on MAT8 entry were strain allowables

               ULT_STRN(1,1) = EMAT(11,1)                     ! Dir 1 tension  allowable
               ULT_STRN(2,1) = EMAT(12,1)                     ! Dir 1 compr    allowable
               ULT_STRN(3,1) = EMAT(13,1)                     ! Dir 2 tension  allowable
               ULT_STRN(4,1) = EMAT(14,1)                     ! Dir 2 compr    allowable
               ULT_STRN(7,1) = EMAT(15,1)                     ! In-plane shear allowable

               ULT_STRE(1,1) = ZERO
               ULT_STRE(2,1) = ZERO
               ULT_STRE(3,1) = ZERO
               ULT_STRE(4,1) = ZERO
               ULT_STRE(7,1) = ZERO

               IF (DABS(E1) > EPS1) THEN
                  ULT_STRE(1,1) = ULT_STRN(1,1)*E1
               ENDIF

               IF (DABS(E2) > EPS1) THEN
                  ULT_STRE(2,1) = ULT_STRN(2,1)*E2
               ENDIF

               IF (DABS(G12) > EPS1) THEN
                  ULT_STRE(3,1) = ULT_STRN(3,1)*G12
               ENDIF

               IF (DABS(G1Z) > EPS1) THEN
                  ULT_STRE(4,1) = ULT_STRN(4,1)*G1Z
               ENDIF

               IF (DABS(G2Z) > EPS1) THEN
                  ULT_STRE(7,1) = ULT_STRN(7,1)*G2Z
               ENDIF

            ENDIF

            DEN1 = ONE - NU12*NU21
            IF (DABS(DEN1) > EPS1) THEN                    ! If 1-NU12*NU21 = 0, error

               EM(1,1) = E1/DEN1
               EM(2,2) = E2/DEN1
               EM(3,3) = G12
               EM(1,2) = NU21*E1/DEN1
               EM(2,1) = EM(1,2)


               ALPVEC(1,1) = ALPHA1
               ALPVEC(2,1) = ALPHA2


            ELSE

               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               IERROR = IERROR + 1
               IF (WRT_ERR > 0) THEN
                  WRITE(ERR,1930) DEN1, EID, TYPE
                  WRITE(F06,1930) DEN1, EID, TYPE
               ELSE
                  ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
                  EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1930
                  EMG_RFE(NUM_EMG_FATAL_ERRS,1)   = DEN1
               ENDIF

            ENDIF

         ELSE

            PROG_ERR = PROG_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            MESSAG = 'MEMBRANE'
            IERROR = IERROR + 1
            IF (WRT_ERR > 0) THEN
               WRITE(ERR,1903) SUBR_NAME, MESSAG, MTRL_TYPE(1)
               WRITE(F06,1903) SUBR_NAME, MESSAG, MTRL_TYPE(1)
            ELSE
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1903
               EMG_IFE(NUM_EMG_FATAL_ERRS,2)   = 1
               EMG_IFE(NUM_EMG_FATAL_ERRS,3)   = MTRL_TYPE(1)
            ENDIF

         ENDIF

      ENDIF mem

! **********************************************************************************************************************************
! Bending stress material props.

bend: IF (MTRL_TYPE(2) /= 0) THEN

         IF (MTRL_TYPE(2) == 1) THEN                       ! Isotropic properties
!                                                            --------------------
            E             = EMAT( 1,2)
            G             = EMAT( 2,2)
            NU            = EMAT( 3,2)
            RHO(2)        = EMAT( 4,2)
            ALPHA         = EMAT( 5,2)
            TREF(2)       = EMAT( 6,2)

            ULT_STRE(1,2) = EMAT( 8,2)                     ! Dir 1 tension  allowable
            ULT_STRE(2,2) = EMAT( 9,2)                     ! Dir 1 compr    allowable
            ULT_STRE(3,2) = ULT_STRE(1,2)                  ! Dir 2 tension  allowable
            ULT_STRE(4,2) = ULT_STRE(2,2)                  ! Dir 2 compr    allowable
            ULT_STRE(7,2) = EMAT(10,2)                     ! In-plane shear allowable

            ULT_STRN(1,2) = ZERO
            ULT_STRN(2,2) = ZERO
            ULT_STRN(3,2) = ZERO
            ULT_STRN(4,2) = ZERO
            ULT_STRN(7,2) = ZERO

            IF (DABS(E) > EPS1) THEN
               ULT_STRN(1,2) = ULT_STRE(1,2)/E
               ULT_STRN(2,2) = ULT_STRE(2,2)/E
               ULT_STRN(3,2) = ULT_STRE(3,2)/E
               ULT_STRN(4,2) = ULT_STRE(4,2)/E
            ENDIF

            IF (DABS(G) > EPS1) THEN
               ULT_STRN(7,2) = ULT_STRE(7,2)/G
            ENDIF

            DEN1 = ONE - NU*NU
            IF (DABS(DEN1) > EPS1) THEN                    ! If 1-NU^2 = 0, error and stop.

               E0 = E/DEN1
               EB(1,1) = E0
               EB(2,2) = EB(1,1)
               EB(3,3) = G
               EB(1,2) = E0*NU
               EB(2,1) = EB(1,2)

               ALPVEC(1,2) = ALPHA
               ALPVEC(2,2) = ALPHA

            ELSE

               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               IERROR = IERROR + 1
               IF (WRT_ERR > 0) THEN
                  WRITE(ERR,1924) DEN1, EID, TYPE
                  WRITE(F06,1924) DEN1, EID, TYPE
               ELSE
                  ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
                  EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1924
                  EMG_RFE(NUM_EMG_FATAL_ERRS,1)   = DEN1
               ENDIF

            ENDIF

         ELSE IF (MTRL_TYPE(2) == 2) THEN                  ! Anisotropic properties
!                                                            ----------------------
            EB(1,1) = EMAT( 1,2)
            EB(1,2) = EMAT( 2,2)
            EB(1,3) = EMAT( 3,2)

            EB(2,1) = EB(1,2)
            EB(2,2) = EMAT( 4,2)
            EB(2,3) = EMAT( 5,2)

            EB(3,1) = EB(1,3)
            EB(3,2) = EB(2,3)
            EB(3,3) = EMAT( 6,2)

            RHO(2)        = EMAT( 7,2)
            ALPHA1        = EMAT( 8,2)
            ALPHA2        = EMAT( 9,2)
            ALPHA3        = EMAT(10,2)
            TREF(2)       = EMAT(11,2)

            ULT_STRE(1,2) = EMAT(13,2)                     ! Dir 1 tension  allowable
            ULT_STRE(2,2) = EMAT(14,2)                     ! Dir 1 compr    allowable
            ULT_STRE(3,2) = ULT_STRE(1,2)                  ! Dir 2 tension  allowable
            ULT_STRE(4,2) = ULT_STRE(2,2)                  ! Dir 2 compr    allowable
            ULT_STRE(7,2) = EMAT(15,2)                     ! In-plane shear allowable

            ULT_STRN(1,2) = ZERO
            ULT_STRN(2,2) = ZERO
            ULT_STRN(3,2) = ZERO
            ULT_STRN(4,2) = ZERO
            ULT_STRN(7,2) = ZERO

            IF (DABS(EM(1,1)) > EPS1) THEN
               ULT_STRN(1,2) = ULT_STRE(1,2)/EM(1,1)
               ULT_STRN(2,2) = ULT_STRE(2,2)/EM(1,1)
            ENDIF

            IF (DABS(EM(2,2)) > EPS1) THEN
               ULT_STRN(3,2) = ULT_STRE(3,2)/EM(2,2)
               ULT_STRN(4,2) = ULT_STRE(4,2)/EM(2,2)
            ENDIF

            IF (DABS(EM(3,3)) > EPS1) THEN
               ULT_STRN(7,2) = ULT_STRE(7,2)/EM(3,3)
            ENDIF

            ALPVEC(1,2)   = ALPHA1
            ALPVEC(2,2)   = ALPHA2
            ALPVEC(3,2)   = ALPHA3

         ELSE IF (MTRL_TYPE(2) == 8) THEN                  ! Orthotropic properties
!                                                            ----------------------
            E1            = EMAT( 1,2)
            E2            = EMAT( 2,2)
            NU12          = EMAT( 3,2)
            NU21          = NU12*(E2/E1)                   ! E1, E2 /= 0 checked in BD_MAT8 subr
            G12           = EMAT( 4,2)
            G1Z           = EMAT( 5,2)
            G2Z           = EMAT( 6,2)
            RHO(2)        = EMAT( 7,2)
            ALPHA1        = EMAT( 8,2)
            ALPHA2        = EMAT( 9,2)
            TREF(2)       = EMAT(10,2)

            IF (EMAT(18,2) == 0) THEN                         ! This indicates allowables on the MAT8 entry were stress allowables

               ULT_STRE(1,2) = EMAT(11,2)                     ! Dir 1 tension  allowable
               ULT_STRE(2,2) = EMAT(12,2)                     ! Dir 1 compr    allowable
               ULT_STRE(3,2) = EMAT(13,2)                     ! Dir 2 tension  allowable
               ULT_STRE(4,2) = EMAT(14,2)                     ! Dir 2 compr    allowable
               ULT_STRE(7,2) = EMAT(15,2)                     ! In-plane shear allowable

               ULT_STRN(1,2) = ZERO
               ULT_STRN(2,2) = ZERO
               ULT_STRN(3,2) = ZERO
               ULT_STRN(4,2) = ZERO
               ULT_STRN(7,2) = ZERO

               IF (DABS(E1) > EPS1) THEN
                  ULT_STRN(1,2) = ULT_STRE(1,2)/E1
               ENDIF

               IF (DABS(E2) > EPS1) THEN
                  ULT_STRN(2,2) = ULT_STRE(2,2)/E2
               ENDIF

               IF (DABS(G12) > EPS1) THEN
                  ULT_STRN(3,2) = ULT_STRE(3,2)/G12
               ENDIF

               IF (DABS(G1Z) > EPS1) THEN
                  ULT_STRN(4,2) = ULT_STRE(3,2)/G1Z
               ENDIF

               IF (DABS(G2Z) > EPS1) THEN
                  ULT_STRN(7,2) = ULT_STRE(7,2)/G2Z
               ENDIF

            ELSE                                           ! Allowables on MAT8 entry were strain allowables

               ULT_STRN(1,2) = EMAT(11,2)                     ! Dir 1 tension  allowable
               ULT_STRN(2,2) = EMAT(12,2)                     ! Dir 1 compr    allowable
               ULT_STRN(3,2) = EMAT(13,2)                     ! Dir 2 tension  allowable
               ULT_STRN(4,2) = EMAT(14,2)                     ! Dir 2 compr    allowable
               ULT_STRN(7,2) = EMAT(15,2)                     ! In-plane shear allowable

               ULT_STRE(1,2) = ZERO
               ULT_STRE(2,2) = ZERO
               ULT_STRE(3,2) = ZERO
               ULT_STRE(4,2) = ZERO
               ULT_STRE(7,2) = ZERO

               IF (DABS(E1) > EPS1) THEN
                  ULT_STRE(1,2) = ULT_STRN(1,2)*E1
               ENDIF

               IF (DABS(E2) > EPS1) THEN
                  ULT_STRE(2,2) = ULT_STRN(2,2)*E2
               ENDIF

               IF (DABS(G12) > EPS1) THEN
                  ULT_STRE(3,2) = ULT_STRN(3,2)*G12
               ENDIF

               IF (DABS(G1Z) > EPS1) THEN
                  ULT_STRE(4,2) = ULT_STRN(4,2)*G1Z
               ENDIF

               IF (DABS(G2Z) > EPS1) THEN
                  ULT_STRE(7,2) = ULT_STRN(7,2)*G2Z
               ENDIF

            ENDIF

            DEN1 = ONE - NU12*NU21
            IF (DABS(DEN1) > EPS1) THEN                    ! If 1-NU12*NU21 = 0, error

               EB(1,1) = E1/DEN1
               EB(2,2) = E2/DEN1
               EB(3,3) = G12
               EB(1,2) = NU12*E2/DEN1
               EB(2,1) = EB(1,2)

               ALPVEC(1,2) = ALPHA1
               ALPVEC(2,2) = ALPHA2

            ELSE

               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               IERROR = IERROR + 1
               IF (WRT_ERR > 0) THEN
                  WRITE(ERR,1930) DEN1, EID, TYPE
                  WRITE(F06,1930) DEN1, EID, TYPE
               ELSE
                  ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
                  EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1930
                  EMG_RFE(NUM_EMG_FATAL_ERRS,1)   = DEN1
               ENDIF

            ENDIF

         ELSE

            PROG_ERR = PROG_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            MESSAG = 'BENDING'
            IERROR = IERROR + 1
            IF (WRT_ERR > 0) THEN
               WRITE(ERR,1903) SUBR_NAME, MESSAG, MTRL_TYPE(2)
               WRITE(F06,1903) SUBR_NAME, MESSAG, MTRL_TYPE(2)
            ELSE
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1903
               EMG_IFE(NUM_EMG_FATAL_ERRS,2)   = 2
               EMG_IFE(NUM_EMG_FATAL_ERRS,3)   = MTRL_TYPE(2)
            ENDIF

         ENDIF

      ENDIF bend

! **********************************************************************************************************************************
! Transverse shear stress material properties

tshr: IF (MTRL_TYPE(3) /= 0) THEN

         IF (MTRL_TYPE(3) == 1) THEN                       ! Isotropic properties
!                                                            --------------------
            G       = EMAT(2,3)
            RHO(3)  = EMAT(4,3)
            ALPHA   = EMAT(5,3)
            TREF(3) = EMAT(6,3)

            ET(1,1) = G
            ET(2,2) = ET(1,1)

            ALPVEC(1,3) = ALPHA
            ALPVEC(2,3) = ALPHA

            ULT_STRE(8,3) = EMAT(19,3)                     ! EMAT(19,3) is either same as EMAT(10,3) or SB (if PCOMP)
            ULT_STRE(9,3) = ULT_STRE(8,3)

            IF (DABS(G) > EPS1) THEN
               ULT_STRN(8,3) = ULT_STRE(8,3)/G
               ULT_STRN(9,3) = ULT_STRE(9,3)/G
            ENDIF

         ELSE IF (MTRL_TYPE(3) == 2) THEN                  ! Anisotropic properties

            ET(1,1) = EMAT( 1,3)                           ! Use same values as EM - see "Linear Static Analysis" guide, p 137
            ET(1,2) = EMAT( 2,3)
            ET(2,2) = EMAT( 4,3)
            ET(2,1) = ET(1,2)

            ULT_STRE(8,3) = EMAT(19,3)                     ! EMAT(19,3) is either same as EMAT(10,3) or SB (if PCOMP)
            ULT_STRE(9,3) = ULT_STRE(8,3)

            IF (DABS(ET(1,1)) > EPS1) THEN
               ULT_STRN(8,3) = ULT_STRE(8,3)/ET(1,1)
            ENDIF
            IF (DABS(ET(2,2)) > EPS1) THEN
               ULT_STRN(9,3) = ULT_STRE(9,3)/ET(2,2)
            ENDIF

         ELSE IF (MTRL_TYPE(3) == 8) THEN                  ! Orthotropic properties
!                                                            ----------------------
            G1Z = EMAT(5,3)
            G2Z = EMAT(6,3)

            ET(1,1) = G1Z
            ET(2,2) = G2Z

            IF (EMAT(18,3) == 0) THEN

               ULT_STRE(8,3) = EMAT(19,3)                  ! EMAT(19,3) is either same as EMAT(10,3) or SB (if PCOMP)
               ULT_STRE(9,3) = ULT_STRE(8,3)

               ULT_STRN(8,3) = ZERO
               ULT_STRN(9,3) = ZERO

               IF (DABS(G1Z) > EPS1) THEN
                  ULT_STRN(8,3) = ULT_STRE(8,3)/G1Z
               ENDIF

               IF (DABS(G2Z) > EPS1) THEN
                  ULT_STRN(9,3) = ULT_STRE(9,3)/G2Z
               ENDIF

            ELSE

            ENDIF

         ELSE

            PROG_ERR = PROG_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            MESSAG = 'TRANSVERSE SHEAR'
            IERROR = IERROR + 1
            IF (WRT_ERR > 0) THEN
               WRITE(ERR,1903) SUBR_NAME, MESSAG, MTRL_TYPE(3)
               WRITE(F06,1903) SUBR_NAME, MESSAG, MTRL_TYPE(3)
            ELSE
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1903
               EMG_IFE(NUM_EMG_FATAL_ERRS,2)   = 3
               EMG_IFE(NUM_EMG_FATAL_ERRS,3)   = MTRL_TYPE(3)
            ENDIF

         ENDIF

      ENDIF tshr

! **********************************************************************************************************************************
! Bending/membrane coupling stress material props.

coup: IF (MTRL_TYPE(4) /= 0) THEN

         IF (MTRL_TYPE(4) == 1) THEN                       ! Isotropic properties
!                                                            --------------------
            E             = EMAT( 1,4)
            G             = EMAT( 2,4)
            NU            = EMAT( 3,4)
            RHO(4)        = EMAT( 4,4)
            ALPHA         = EMAT( 5,4)
            TREF(4)       = EMAT( 6,4)

            ULT_STRE(1,4) = EMAT( 8,4)                     ! Dir 1 tension  allowable
            ULT_STRE(2,4) = EMAT( 9,4)                     ! Dir 1 compr    allowable
            ULT_STRE(3,4) = ULT_STRE(1,2)                  ! Dir 2 tension  allowable
            ULT_STRE(4,4) = ULT_STRE(2,4)                  ! Dir 2 compr    allowable
            ULT_STRE(7,4) = EMAT(10,4)                     ! In-plane shear allowable

            ULT_STRN(1,2) = ZERO
            ULT_STRN(2,2) = ZERO
            ULT_STRN(3,2) = ZERO
            ULT_STRN(4,2) = ZERO
            ULT_STRN(7,2) = ZERO

            IF (DABS(E) > EPS1) THEN
               ULT_STRN(1,4) = ULT_STRE(1,4)/E
               ULT_STRN(2,4) = ULT_STRE(2,4)/E
               ULT_STRN(3,4) = ULT_STRE(3,4)/E
               ULT_STRN(4,4) = ULT_STRE(4,4)/E
            ENDIF

            IF (DABS(G) > EPS1) THEN
               ULT_STRN(7,4) = ULT_STRE(7,4)/G
            ENDIF

            DEN1 = ONE - NU*NU
            IF (DABS(DEN1) > EPS1) THEN                    ! If 1-NU^2 = 0, error and stop.

               E0 = E/DEN1
               EBM(1,1) = E0
               EBM(2,2) = EBM(1,1)
               EBM(3,3) = G
               EBM(1,2) = E0*NU
               EBM(2,1) = EBM(1,2)

               ALPVEC(1,4) = ALPHA
               ALPVEC(2,4) = ALPHA

            ELSE

               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               IERROR = IERROR + 1
               IF (WRT_ERR > 0) THEN
                  WRITE(ERR,1924) DEN1, EID, TYPE
                  WRITE(F06,1924) DEN1, EID, TYPE
               ELSE
                  ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
                  EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1924
                  EMG_RFE(NUM_EMG_FATAL_ERRS,1)   = DEN1
               ENDIF

            ENDIF

         ELSE IF (MTRL_TYPE(4) == 2) THEN                  ! Anisotropic properties
!                                                            ----------------------
            EBM(1,1) = EMAT( 1,4)
            EBM(1,2) = EMAT( 2,4)
            EBM(1,3) = EMAT( 3,4)

            EBM(2,1) = EBM(1,2)
            EBM(2,2) = EMAT( 4,4)
            EBM(2,3) = EMAT( 5,4)

            EBM(3,1) = EBM(1,3)
            EBM(3,2) = EBM(2,3)
            EBM(3,3) = EMAT( 6,4)

            RHO(4)        = EMAT( 8,4)
            ALPHA1        = EMAT( 8,4)
            ALPHA2        = EMAT( 9,4)
            ALPHA3        = EMAT(10,4)
            TREF(4)       = EMAT(11,4)

            ULT_STRE(1,4) = EMAT(13,4)                     ! Dir 1 tension  allowable
            ULT_STRE(2,4) = EMAT(14,4)                     ! Dir 1 compr    allowable
            ULT_STRE(3,4) = ULT_STRE(1,4)                  ! Dir 2 tension  allowable
            ULT_STRE(4,4) = ULT_STRE(2,4)                  ! Dir 2 compr    allowable
            ULT_STRE(7,4) = EMAT(15,4)                     ! In-plane shear allowable

            ULT_STRN(1,4) = ZERO
            ULT_STRN(2,4) = ZERO
            ULT_STRN(3,4) = ZERO
            ULT_STRN(4,4) = ZERO
            ULT_STRN(7,4) = ZERO

            IF (DABS(EM(1,1)) > EPS1) THEN
               ULT_STRN(1,4) = ULT_STRE(1,4)/EM(1,1)
               ULT_STRN(2,4) = ULT_STRE(2,4)/EM(1,1)
            ENDIF

            IF (DABS(EM(2,2)) > EPS1) THEN
               ULT_STRN(3,4) = ULT_STRE(3,4)/EM(2,2)
               ULT_STRN(4,4) = ULT_STRE(4,4)/EM(2,2)
            ENDIF

            IF (DABS(EM(3,3)) > EPS1) THEN
               ULT_STRN(7,4) = ULT_STRE(7,4)/EM(3,3)
            ENDIF

            ALPVEC(1,4) = ALPHA1
            ALPVEC(2,4) = ALPHA2
            ALPVEC(3,4) = ALPHA3

         ELSE IF (MTRL_TYPE(4) == 8) THEN                  ! Orthotropic properties
!                                                            ----------------------
            E1            = EMAT( 1,4)
            E2            = EMAT( 2,4)
            NU12          = EMAT( 3,4)
            NU21          = NU12*(E2/E1)                   ! E1, E2 /= 0 checked in BD_MAT8 subr
            G12           = EMAT( 4,4)
            G1Z           = EMAT( 5,4)
            G2Z           = EMAT( 6,4)
            RHO(4)        = EMAT( 7,4)
            ALPHA1        = EMAT( 8,4)
            ALPHA2        = EMAT( 9,4)
            TREF(4)       = EMAT(10,4)

            IF (EMAT(18,4) == 0) THEN                         ! This indicates allowables on the MAT8 entry were stress allowables

               ULT_STRE(1,4) = EMAT(11,4)                     ! Dir 1 tension  allowable
               ULT_STRE(2,4) = EMAT(12,4)                     ! Dir 1 compr    allowable
               ULT_STRE(3,4) = EMAT(13,4)                     ! Dir 2 tension  allowable
               ULT_STRE(4,4) = EMAT(14,4)                     ! Dir 2 compr    allowable
               ULT_STRE(7,4) = EMAT(15,4)                     ! In-plane shear allowable

               ULT_STRN(1,4) = ZERO
               ULT_STRN(2,4) = ZERO
               ULT_STRN(3,4) = ZERO
               ULT_STRN(4,4) = ZERO
               ULT_STRN(7,4) = ZERO

               IF (DABS(E1) > EPS1) THEN
                  ULT_STRN(1,4) = ULT_STRE(1,4)/E1
               ENDIF

               IF (DABS(E2) > EPS1) THEN
                  ULT_STRN(2,4) = ULT_STRE(2,4)/E2
               ENDIF

               IF (DABS(G12) > EPS1) THEN
                  ULT_STRN(3,4) = ULT_STRE(3,4)/G12
               ENDIF

               IF (DABS(G1Z) > EPS1) THEN
                  ULT_STRN(4,4) = ULT_STRE(4,4)/G1Z
               ENDIF

               IF (DABS(G2Z) > EPS1) THEN
                  ULT_STRN(7,4) = ULT_STRE(7,4)/G2Z
               ENDIF

            ELSE                                           ! Allowables on MAT8 entry were strain allowables

               ULT_STRN(1,4) = EMAT(11,4)                     ! Dir 1 tension  allowable
               ULT_STRN(2,4) = EMAT(12,4)                     ! Dir 1 compr    allowable
               ULT_STRN(3,4) = EMAT(13,4)                     ! Dir 2 tension  allowable
               ULT_STRN(4,4) = EMAT(14,4)                     ! Dir 2 compr    allowable
               ULT_STRN(7,4) = EMAT(15,4)                     ! In-plane shear allowable

               ULT_STRE(1,4) = ZERO
               ULT_STRE(2,4) = ZERO
               ULT_STRE(3,4) = ZERO
               ULT_STRE(4,4) = ZERO
               ULT_STRE(7,4) = ZERO

               IF (DABS(E1) > EPS1) THEN
                  ULT_STRE(1,4) = ULT_STRN(1,4)*E1
               ENDIF

               IF (DABS(E2) > EPS1) THEN
                  ULT_STRE(2,4) = ULT_STRN(2,4)*E2
               ENDIF

               IF (DABS(G12) > EPS1) THEN
                  ULT_STRE(3,4) = ULT_STRN(3,4)*G12
               ENDIF

               IF (DABS(G1Z) > EPS1) THEN
                  ULT_STRE(4,4) = ULT_STRN(3,4)*G1Z
               ENDIF

               IF (DABS(G2Z) > EPS1) THEN
                  ULT_STRE(7,4) = ULT_STRN(7,4)*G2Z
               ENDIF

            ENDIF

            DEN1 = ONE - NU12*NU21
            IF (DABS(DEN1) > EPS1) THEN                    ! If 1-NU12*NU21 = 0, error

               EBM(1,1) = E1/DEN1
               EBM(2,2) = E2/DEN1
               EBM(3,3) = G12
               EBM(1,2) = NU12*E2/DEN1
               EBM(2,1) = EBM(1,2)

               ALPVEC(1,4) = ALPHA1
               ALPVEC(2,4) = ALPHA2

            ELSE

               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               IERROR = IERROR + 1
               IF (WRT_ERR > 0) THEN
                  WRITE(ERR,1930) DEN1, EID, TYPE
                  WRITE(F06,1930) DEN1, EID, TYPE
               ELSE
                  ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
                  EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1930
                  EMG_RFE(NUM_EMG_FATAL_ERRS,1)   = DEN1
               ENDIF

            ENDIF

         ELSE

            PROG_ERR = PROG_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            MESSAG = 'BENDING'
            IERROR = IERROR + 1
            IF (WRT_ERR > 0) THEN
               WRITE(ERR,1903) SUBR_NAME, MESSAG, MTRL_TYPE(4)
               WRITE(F06,1903) SUBR_NAME, MESSAG, MTRL_TYPE(4)
            ELSE
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1)   = 1903
               EMG_IFE(NUM_EMG_FATAL_ERRS,2)   = 4
               EMG_IFE(NUM_EMG_FATAL_ERRS,3)   = MTRL_TYPE(4)
            ENDIF

         ENDIF

      ENDIF coup

! **********************************************************************************************************************************
! Quit if there was an invalid material type

      IF (PROG_ERR > 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
! Return if IERROR > 0

      IF (IERROR > 0) RETURN



      RETURN

! **********************************************************************************************************************************
 1903 FORMAT(' *ERROR  1903: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATERIAL TYPE FOR ',A,' MUST BE 1 OR 8 FOR 2D ELEMENTS BUT TYPE WAS ',I8)

 1924 FORMAT(' *ERROR  1924: MATERIAL ERROR. 1 - NU^2 = ',1ES9.2,' IS TOO CLOSE TO ZERO FOR ELEMENT ',I8,', TYPE ',A)

 1930 FORMAT(' *ERROR  1930: MATERIAL ERROR. 1 - NU12*NU21 = ',1ES9.2,' IS TOO CLOSE TO ZERO FOR ELEMENT ',I8,', TYPE ',A)

 1998 format(' *ERROR  1999: Current code for ', a, ' element ', i8, ' is incomplete for Orthotropic material properties.'         &
                    ,/,14x,' Try a PARAM, QUAD4TYP, MIN4 for an alternate QUAD4 if your elements are nearly rectangular.'          &
                    ,/,14x,' Otherwise you will have to use CTRIA3 elements'                                                       &
                    ,/,14x,' You can override this and continue with a bdf entry: DEBUG, 248, 1')

 1999 format(          14x,' Fatal error stop overridden with DEBUG(248) > 0')

! **********************************************************************************************************************************

      END SUBROUTINE MATERIAL_PROPS_2D


      SUBROUTINE MATERIAL_PROPS_3D ( WRITE_WARN )

! Calculates material stress/strain matrices for isotropic or anisotropic 3-D solid elements

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  ALPVEC, EID, ES, EMAT, NUM_EMG_FATAL_ERRS, MTRL_TYPE, RHO, ULT_STRE, ULT_STRN, TREF, TYPE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATERIAL_PROPS_3D'
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not

      INTEGER(LONG)                   :: IERROR       = 0  ! Local error indicator
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: K                 ! Counter


      REAL(DOUBLE)                    :: ALPHA             ! Isotropic coefficient of thermal expansion
      REAL(DOUBLE)                    :: DEN1              ! An intermaediate variable in calculating outputs
      REAL(DOUBLE)                    :: DEN2              ! An intermaediate variable in calculating outputs
      REAL(DOUBLE)                    :: E                 ! Isotropic Young's modulus
      REAL(DOUBLE)                    :: E0                ! An intermaediate variable in calculating outputs
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare for real zero
      REAL(DOUBLE)                    :: G                 ! Shear modulus
      REAL(DOUBLE)                    :: NU                ! Poisson's ratio
      REAL(DOUBLE)                    :: NU0               ! An intermaediate variable in calculating outputs

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1   = EPSIL(1)

      ! 3D material matrix
      ES = ZERO

      IF (MTRL_TYPE(1) == 1) THEN                          ! Isotropic properties

         E         = EMAT( 1,1)
         G         = EMAT( 2,1)
         NU        = EMAT( 3,1)
         RHO(1)    = EMAT( 4,1)
         ALPHA     = EMAT( 5,1)
         TREF(1)   = EMAT( 6,1)

         ULT_STRE(1,1) = EMAT( 8,1)                        ! Dir 1 tension  allowable
         ULT_STRE(2,1) = EMAT( 9,1)                        ! Dir 1 compr    allowable
         ULT_STRE(3,1) = ULT_STRE(1,1)                     ! Dir 2 tension  allowable
         ULT_STRE(4,1) = ULT_STRE(2,1)                     ! Dir 2 compr    allowable
         ULT_STRE(5,1) = ULT_STRE(1,1)                     ! Dir 3 tension  allowable
         ULT_STRE(6,1) = ULT_STRE(2,1)                     ! Dir 3 compr    allowable
         ULT_STRE(7,1) = EMAT(10,1)                        ! Plane 12 shear allowable
         ULT_STRE(8,1) = EMAT(10,1)                        ! Plane 23 shear allowable
         ULT_STRE(9,1) = EMAT(10,1)                        ! Plane 31 shear allowable

         ULT_STRN(1,1) = ZERO
         ULT_STRN(2,1) = ZERO
         ULT_STRN(3,1) = ZERO
         ULT_STRN(4,1) = ZERO
         ULT_STRN(5,1) = ZERO
         ULT_STRN(6,1) = ZERO
         ULT_STRN(7,1) = ZERO
         ULT_STRN(8,1) = ZERO
         ULT_STRN(9,1) = ZERO

         IF (DABS(E) > EPS1) THEN
            ULT_STRN(1,1) = ULT_STRE(1,1)/E
            ULT_STRN(2,1) = ULT_STRE(2,1)/E
            ULT_STRN(3,1) = ULT_STRE(3,1)/E
            ULT_STRN(4,1) = ULT_STRE(4,1)/E
            ULT_STRN(5,1) = ULT_STRE(5,1)/E
            ULT_STRN(6,1) = ULT_STRE(6,1)/E
         ENDIF

         IF (DABS(G) > EPS1) THEN
            ULT_STRN(7,1) = ULT_STRE(7,1)/G
            ULT_STRN(8,1) = ULT_STRE(8,1)/G
            ULT_STRN(9,1) = ULT_STRE(9,1)/G
         ENDIF

         DEN1 = ONE + NU
         IF (DABS(DEN1) <= EPS1) THEN
            WRITE(ERR,1919) EID,TYPE
            WRITE(F06,1919) EID,TYPE
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
            IERROR = IERROR + 1
         ENDIF

         DEN2 = ONE - TWO*NU
         IF (DABS(DEN2) <= EPS1) THEN
            WRITE(ERR,1936) EID,TYPE
            WRITE(F06,1936) EID,TYPE
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
            IERROR = IERROR + 1
         ENDIF

         E0  = E/DEN1
         NU0 = (ONE - NU)/DEN2

         ES(1,1) = E0*NU0
         ES(1,2) = E0*NU/DEN2
         ES(1,3) = ES(1,2)

         ES(2,1) = ES(1,2)
         ES(2,2) = ES(1,1)
         ES(2,3) = ES(1,2)

         ES(3,1) = ES(1,2)
         ES(3,2) = ES(1,2)
         ES(3,3) = ES(1,1)

         ES(4,4) = G
         ES(5,5) = ES(4,4)
         ES(6,6) = ES(4,4)

         ALPVEC(1,1) = ALPHA
         ALPVEC(2,1) = ALPVEC(1,1)
         ALPVEC(3,1) = ALPVEC(1,1)
         ALPVEC(4,1) = ZERO
         ALPVEC(5,1) = ZERO
         ALPVEC(6,1) = ZERO

      ELSE IF (MTRL_TYPE(1) == 9) THEN                          ! Anisotropic properties

         K = 0
         DO I=1,6
            DO J=I,6
               K = K + 1
               ES(I,J) = EMAT(K,1)
               ES(J,I) = ES(I,J)
            ENDDO
         ENDDO

         RHO(1)      = EMAT(22,1)

         ALPVEC(1,1) = EMAT(23,1)
         ALPVEC(2,1) = EMAT(24,1)
         ALPVEC(3,1) = EMAT(25,1)
         ALPVEC(4,1) = EMAT(26,1)
         ALPVEC(5,1) = EMAT(27,1)
         ALPVEC(6,1) = EMAT(28,1)

         TREF(1)     = EMAT(29,1)

         DO I=1,9
            ULT_STRE(I,1) = ZERO
            ULT_STRN(I,1) = ZERO
         ENDDO

      ELSE

         WRITE(ERR,1926) MTRL_TYPE(1)
         WRITE(F06,1926) MTRL_TYPE(1)
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
         IERROR = IERROR + 1

      ENDIF

! Return if IERROR > 0

      IF (IERROR > 0) RETURN

! **********************************************************************************************************************************
! Return if IERROR > 0

      IF (IERROR > 0) RETURN



      RETURN

! **********************************************************************************************************************************
 1919 FORMAT(' *ERROR  1919: MATERIAL ERROR. 1 + NU IS TOO CLOSE TO ZERO FOR ELEMENT ',I8,', TYPE ',A)

 1926 FORMAT(' *ERROR  1926: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATERIAL TYPE MUST BE 1 FOR 3D ELEMENTS BUT TYPE WAS ',I8)

 1936 FORMAT(' *ERROR  1936: MATERIAL ERROR. 1 - 2*NU IS TOO CLOSE TO ZERO FOR ELEMENT ',I8,', TYPE ',A)


! **********************************************************************************************************************************
      END SUBROUTINE MATERIAL_PROPS_3D

   END MODULE MATERIAL_PROPERTIES
