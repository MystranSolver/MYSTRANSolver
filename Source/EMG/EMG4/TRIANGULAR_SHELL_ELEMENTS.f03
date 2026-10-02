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

   MODULE TRIANGULAR_SHELL_ELEMENTS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: TREL1, TPLT2

   CONTAINS

      SUBROUTINE TREL1 ( OPT, WRITE_WARN )

! Calculates, or calls subr's to calculate, triangular element matrices:

!  1) ME        = element mass matrix                  , if OPT(1) = 'Y'
!  2) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  3) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  4) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  5) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  6) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MEWE, NSUB, NTSUB, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TENTH, ONE, TWO, THREE, TWELVE
      USE PARAMS, ONLY                :  SUPWARN
      USE MODEL_STUF, ONLY            :  EID, ELDOF, EMG_IWE, EMG_RWE, INTL_MID, KE, MASS_PER_UNIT_AREA, ME,                       &
                                         NUM_EMG_FATAL_ERRS, PCOMP_LAM, PCOMP_PROPS, SHELL_B, TYPE, XEB, XEL
      USE MODEL_STUF, ONLY            :  BENSUM, SHRSUM, PHI_SQ, PSI_HAT, XTB, XTL

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TREL1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not

      INTEGER(LONG)                   :: IERROR            ! Local error indicator from one of the subrs called
      INTEGER(LONG)                   :: K,L               ! DO loop indices


      REAL(DOUBLE)                    :: AR                ! Elem aspect ratio
      REAL(DOUBLE)                    :: AREA              ! Elem area

      REAL(DOUBLE)                    :: B2V(3,9)          ! The 3x9     virgin strain   recovery matrix for MIN3 for bending
      REAL(DOUBLE)                    :: B3V(3,9)          ! The 3x9     virgin strain   recovery matrix for MIN3 for transv shear

      REAL(DOUBLE)                    :: BIG_BB(3,ELDOF,1) ! Strain-displ matrix for bending for all Gauss points and all DOF's

      REAL(DOUBLE)                    :: BIG_BBI(3,ELDOF)  ! BIG_BB for 3rd subscript

      REAL(DOUBLE)                    :: BIG_BM(3,ELDOF,1) ! Strain-displ matrix for this elem for all Gauss points (for all DOF's)

      REAL(DOUBLE)                    :: BIG_BMI(3,ELDOF)  ! BIG_BM for 3rd subscript

      REAL(DOUBLE)                    :: DUM1(3,ELDOF)     ! Intermediate result in calc SHELL_B effect on KE
      REAL(DOUBLE)                    :: DUM2(ELDOF,ELDOF) ! Intermediate result in calc SHELL_B effect on KE
      REAL(DOUBLE)                    :: KV(9,9)           ! KB + PHISQ*KS (the 9x9 virgin stiffness matrix for MIN3)
      REAL(DOUBLE)                    :: M0                ! An intermediate variable used in calc elem mass, ME
      REAL(DOUBLE)                    :: PPV(9,NSUB)       ! The 9xNSUB  virgin thermal  load     matrix for MIN3
      REAL(DOUBLE)                    :: PTV(9,NTSUB)      ! The 9xNTSUB virgin pressure load     matrix for MIN3
      REAL(DOUBLE)                    :: S2V(3,9)          ! The 3x9     virgin stress   recovery matrix for MIN3 for bending
      REAL(DOUBLE)                    :: S3V(3,9)          ! The 3x9     virgin stress   recovery matrix for MIN3 for transv shear
      REAL(DOUBLE)                    :: X2E               ! x coord of elem node 2
      REAL(DOUBLE)                    :: X3E               ! x coord of elem node 3
      REAL(DOUBLE)                    :: Y3E               ! y coord of elem node 3

! The following 3 args are needed when subr TPLT2 is called when that triangular shell element is used in a MIN4T QUAD4 which is
! made up of 4 non-overlapping TPLT2 elements. Since TPLT2 can also be a stand-alone element, we need these args when that occurs.
! The only time TPLT2 is called from this subr is when it is a stand-alone element.

! We give them the PARAMETER attributes below to make sure this subr won't call TPLT2 with any other values. For the MIN4T QUAD4,
! TPLT2 is called in subr QPLT3 with appropriate values of TRIA_NUM, PSI.

      CHARACTER(LEN=1), PARAMETER     :: MN4T_QD   = 'N'
      INTEGER(LONG)   , PARAMETER     :: TRIA_NUM  = 1
      REAL(DOUBLE)    , PARAMETER     :: PSI       = 0.0D0



! **********************************************************************************************************************************
! Initialize

      BENSUM  = ZERO
      SHRSUM  = ZERO
      PSI_HAT = ZERO
      PHI_SQ  = ZERO

! Calculate element geometry parameters from data block XEL

      X2E  = XEL(2,1)
      X3E  = XEL(3,1)
      Y3E  = XEL(3,2)
      AREA = X2E*Y3E/TWO

! XTB, XTL may be needed when TPLT2 calls BBMIN3, BSMIN3. Since TPLT2 is also called from QPLT3, which is made up of 4 TPLT2's,
! we cannot use XEB and XEL since, in that case, they are the values for the MIN4T QUAD4 element geometry, not for the 4 triangles
! making up that quad eleent

      DO K=1,3
         DO L=1,3
            XTB(K,L) = XEB(K,L)
            XTL(K,L) = XEL(K,L)
         ENDDO
      ENDDO

! Calculate and check element aspect ratio, AR. Print warning if AR > 2.0 for TMEM1

      AR = X2E/Y3E
      IF(AR < ONE) THEN
         AR = ONE/AR
      ENDIF
      IF (AR > TWO + TENTH) THEN
         WARN_ERR = WARN_ERR + 1
         IF ((WRT_ERR > 0) .AND. (WRITE_WARN == 'Y')) THEN
            WRITE(ERR,1924) TYPE, EID, AR, TWO
            IF (SUPWARN == 'N') THEN
               WRITE(F06,1924) TYPE, EID, AR, TWO
            ENDIF
         ELSE
            IF (WARN_ERR <= MEWE) THEN
               EMG_IWE(WARN_ERR,1) = 1924
               EMG_RWE(WARN_ERR,1) = AR
               EMG_RWE(WARN_ERR,2) = TWO
            ENDIF
         ENDIF
      ENDIF

! **********************************************************************************************************************************
! Generate the mass matrix for this element. For the pure bending element the mass is based only on the non-structural mass.
! The mass matrix was initialized in subr EMG

      IF (OPT(1) == 'Y') THEN
         M0 = MASS_PER_UNIT_AREA*AREA/THREE
         ME( 1 ,1) = M0
         ME( 2 ,2) = M0
         ME( 3 ,3) = M0
         ME( 7 ,7) = M0
         ME( 8 ,8) = M0
         ME( 9 ,9) = M0
         ME(13,13) = M0
         ME(14,14) = M0
         ME(15,15) = M0
      ENDIF

! **********************************************************************************************************************************
! If TYPE is 'TRMEM' or 'TRIA3K' or 'TRIA3' generate the membrane stiffness
! If TYPE is 'TRPLT1', 'TRPLT2', 'TRIA3K', or 'TRIA3' generate the bending stiffness

      IF ((OPT(2) == 'Y') .OR. (OPT(3) == 'Y') .OR. (OPT(4) == 'Y') .OR. (OPT(5) == 'Y') .OR. (OPT(6) == 'Y')) THEN

         IF (TYPE(1:5) == 'TRIA3') THEN
            IF (INTL_MID(1) /= 0) THEN
               CALL TMEM1 ( OPT, AREA, X2E, X3E, Y3E, 'Y', BIG_BM )
            ENDIF
         ENDIF

         IF (TYPE == 'TRIA3K  ') THEN
            IF (INTL_MID(2) /= 0) THEN
               CALL TPLT1 ( OPT, AREA, X2E, X3E, Y3E )
            ENDIF
         ENDIF

         IF (TYPE == 'TRIA3   ') THEN
            IF (INTL_MID(2) /= 0) THEN
               CALL TPLT2 (OPT, AREA, X2E, X3E, Y3E, 'Y', IERROR, KV, PTV, PPV, B2V, B3V, S2V, S3V, BIG_BB, MN4T_QD, TRIA_NUM, PSI)
            ENDIF
         ENDIF

      ENDIF


! **********************************************************************************************************************************
! Calc BM'*SHELL_B*BB (and its transpose) and add to KE. Only do this if this is a composite element with nonsym layup

      IF (OPT(4) == 'Y') THEN

         IF (TYPE(1:5) == 'TRIA3') THEN

            IF ((PCOMP_PROPS == 'Y') .AND. (PCOMP_LAM == 'NON')) THEN

               if (type == 'TRIA3K  ') then
                  WRITE(ERR,*) ' *ERROR: Code not written for SHELL_B effect on KE yet for TRIA3K elements'
                  WRITE(F06,*) ' *ERROR :Code not written for SHELL_B effect on KE yet for TRIA3K elements'
                  call outa_here ( 'Y' )
               endif

               DO K=1,3
                  DO L=1,ELDOF
                     BIG_BMI(K,L) = BIG_BM(K,L,1)
                     BIG_BMI(K,L) = BIG_BM(K,L,1)
                  ENDDO
               ENDDO

               DO K=1,3
                  DO L=1,ELDOF
                     BIG_BBI(K,L) = BIG_BB(K,L,1)
                     BIG_BBI(K,L) = BIG_BB(K,L,1)
                  ENDDO
               ENDDO

               CALL MATMULT_FFF ( SHELL_B, BIG_BBI, 3, 3, ELDOF, DUM1 )
               CALL MATMULT_FFF_T ( BIG_BMI, DUM1, 3, ELDOF, ELDOF, DUM2 )

               DO K=1,ELDOF
                  DO L=1,ELDOF
                     KE(K,L) = KE(K,L) + AREA*(DUM2(K,L) + DUM2(L,K))
                  ENDDO
               ENDDO

            ENDIF

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1924 FORMAT(' *WARNING    : ASPECT RATIO OF ',A,' ELEMENT ',I8,' IS:',F7.1,'. IT SHOULD BE < ',F3.0)

! **********************************************************************************************************************************

      END SUBROUTINE TREL1


      SUBROUTINE TMEM1 ( OPT, AREA, X2E, X3E, Y3E, WRT_BUG_THIS_TIME, BIG_BM )

! Constant strain membrane triangle

! Subroutine calculates:

!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  4) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  5) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, BUG, WRT_BUG, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_BCHK_BIT, ELDT_BUG_BMAT_BIT, NSUB, NTSUB, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, THREE
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, BE1, EID, DT, EM, ELDOF, KE, PCOMP_LAM, PCOMP_PROPS, PRESS, PPE, PTE, SE1, STE1,  &
                                         SHELL_AALP, SHELL_A, SHELL_PROP_ALP, TREF, TYPE, XEB, XEL, ELGP, FCONV, STRESS, KED,      &
                                         NUM_EMG_FATAL_ERRS
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE RESULT_COORDINATES, ONLY     :  ELMDIS
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: ELEM_STRE_STRN_ARRAYS

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE QUADRILATERAL_B_MATRICES, ONLY:  BCHECK_2D
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TMEM1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: ID(18)


      REAL(DOUBLE) , INTENT(IN)       :: AREA              ! Element area
      REAL(DOUBLE) , INTENT(IN)       :: X2E               ! x coord of elem node 2
      REAL(DOUBLE) , INTENT(IN)       :: X3E               ! x coord of elem node 3
      REAL(DOUBLE) , INTENT(IN)       :: Y3E               ! y coord of elem node 3

      REAL(DOUBLE) , INTENT(OUT)      :: BIG_BM(3,ELDOF,1) ! Strain-displ matrix for this elem for all Gauss points (for all DOF's)

      REAL(DOUBLE)                    :: BW(3,14)          ! Output from subr BCHECK (matrix of 3 elem strains for 14 various elem
!                                                            rigid body motions/constant strain distortions)

      REAL(DOUBLE)                    :: ALP(3)            ! Col of ALPVEC
      REAL(DOUBLE)                    :: KS(ELGP,ELGP)     ! KED matrix for one DOF
      REAL(DOUBLE)                    :: BM(3,ELDOF)       ! Strain-displ matrix for this elem
      REAL(DOUBLE)                    :: AMB(3,ELDOF)      ! SHELL_A matrix times strain-displ matrix for this elem
      REAL(DOUBLE)                    :: DPSHX(2,ELGP)     ! Derivatives of PSH wrt elem x, y coords.
      REAL(DOUBLE)                    :: DUM(ELDOF,ELDOF)  ! Needed for calc 18 x 18 KE  using MATMULT, since KE is MELDOF x MELDOF
      REAL(DOUBLE)                    :: DUM1(ELDOF,1)     ! Intermediate matrix used in determining PTE thermal loads
      REAL(DOUBLE)                    :: DUM11(2,2)        ! Intermediate matrix used in solving for KED matrices
      REAL(DOUBLE)                    :: DUM12(ELGP,2)     ! Intermediate matrix used in solving for KED matrices
      REAL(DOUBLE)                    :: DUM13(ELGP,ELGP)  ! Intermediate matrix used in solving for KED matrices
      REAL(DOUBLE)                    :: EALP(3)           ! Intermed var used in calc STEi therm stress coeffs
      REAL(DOUBLE)                    :: EMB(3,ELDOF)      ! Mat'l matrix times strain-displ matrix for this elem
      REAL(DOUBLE)                    :: C01               ! Intermediate variable used in calc PTE, SEi, STEi, KE
      REAL(DOUBLE)                    :: C02               ! Intermediate variable used in calc PTE, SEi, STEi, KE
      REAL(DOUBLE)                    :: C03               ! Intermediate variable used in calc PTE, SEi, STEi, KE
      REAL(DOUBLE)                    :: C04               ! Intermediate variable used in calc PTE, SEi, STEi, KE
      REAL(DOUBLE)                    :: CT0               ! Intermediate variable used in calc PTE thermal loads
      REAL(DOUBLE)                    :: TBAR              ! Average elem temperature
      REAL(DOUBLE)                    :: FORCEx            ! Engineering force in the elem x direction
      REAL(DOUBLE)                    :: FORCEy            ! Engineering force in the elem x direction
      REAL(DOUBLE)                    :: FORCExy           ! Engineering force in the elem xy direction




! **********************************************************************************************************************************
! Determine element strain-displacement matrix.

      DO I=1,3
         DO J=1,ELDOF
            BM(I,J) = ZERO
         ENDDO
      ENDDO

      C01 = ONE/X2E
      C02 = ONE/Y3E
      C03 = (X3E - X2E)*C01*C02
      C04 = X3E*C01*C02

      BM(1, 1) = -C01
      BM(1, 7) =  C01

      BM(2, 2) =  C03
      BM(2, 8) = -C04
      BM(2,14) =  C02

      BM(3, 1) =  C03
      BM(3, 2) = -C01
      BM(3, 7) = -C04
      BM(3, 8) =  C01
      BM(3,13) =  C02

      DO I=1,18
         ID(I) = I
      ENDDO

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(8) > 0)) THEN

         WRITE(BUG,1101) ELDT_BUG_BMAT_BIT, TYPE, EID
         WRITE(BUG,8901) SUBR_NAME
         DO I=1,3
            WRITE(BUG,8902) I,(BM(I,J),J=1,ELDOF)
            WRITE(BUG,*)
         ENDDO
         WRITE(BUG,*)

      ENDIF

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(9) > 0)) THEN
        IF (DEBUG(202) > 0) THEN
           WRITE(BUG,1101) ELDT_BUG_BCHK_BIT, TYPE, EID
           WRITE(BUG,9100)
           WRITE(BUG,9101)
           WRITE(BUG,9102)
           WRITE(BUG,9103)
           WRITE(BUG,9104)
           CALL BCHECK_2D ( BM, 'M', ID, 3, 18, 3, XEL, XEB, BW )
        ENDIF
      ENDIF

! **********************************************************************************************************************************
! If element is a composite and if it is a nonsym layup we need to calc BIG_BB for later use

      DO I=1,3
         DO J=1,ELDOF
            BIG_BM(I,J,1) = ZERO
         ENDDO
      ENDDO

      IF ((PCOMP_PROPS == 'Y') .AND. (PCOMP_LAM == 'NON')) THEN

         DO I=1,3
            DO J=1,18
               BIG_BM(I,J,1) = BM(I,J)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Determine element thermal loads.

      IF (OPT(2) == 'Y') THEN

         CALL MATMULT_FFF_T ( BM, SHELL_AALP, 3, ELDOF, 1, DUM1 )

         DO J=1,NTSUB
            TBAR = (DT(1,J) + DT(2,J) + DT(3,J))/THREE
            CT0 = AREA*(TBAR - TREF(1))
            DO I=1,ELDOF
               PTE(I,J) = CT0*DUM1(I,1)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate BE1, SE1 matrices (3 x ELDOF) for strain/stress data recovery.
! Note: strain/stress recovery matrices only make sense for individual plies (or whole elem if only 1 "ply")

      IF (OPT(3) == 'Y' .OR. OPT(6) == "Y") THEN

         DO I=1,3
            DO J=1,ELDOF
               BE1(I,J,1) = BM(I,J)
            ENDDO
         ENDDO

! SE1, STE1 generated in elem coords. Then, in LINK9 the stresses, calc'd in elem coords, will be transformed to ply coords

         CALL MATMULT_FFF ( EM, BM, 3, 3, ELDOF, EMB )     ! Generate SE1 in element coords (at this point EM is elem coords)
         DO I=1,3
            DO J=1,ELDOF
               SE1(I,J,1) = EMB(I,J)
            ENDDO
         ENDDO

         ALP(1) = ALPVEC(1,1)
         ALP(2) = ALPVEC(2,1)
         ALP(3) = ALPVEC(4,1)

         CALL MATMULT_FFF ( EM, ALP, 3, 3, 1, EALP )
         DO J=1,NTSUB
            TBAR = (DT(1,J) + DT(2,J) + DT(3,J))/THREE
            DO I=1,3
               STE1(I,J,1) = EALP(I)*(TBAR - TREF(1))
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE.

      IF (OPT(4) == 'Y') THEN

         CALL MATMULT_FFF ( SHELL_A, BM, 3, 3, ELDOF, AMB )
         CALL MATMULT_FFF_T ( BM, AMB, 3, ELDOF, ELDOF, DUM )
         DO I=1,ELDOF
            DO J=1,ELDOF
               KE(I,J) = KE(I,J) + AREA*DUM(I,J)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate element pressure load matrix PPE.
! NOTE: for this element work equivalent and static equivalent loads are the same

      IF (OPT(5) == 'Y') THEN

         DO J=1,NSUB
            PPE( 1,J) = AREA*PRESS(1,J)/THREE
            PPE( 2,J) = AREA*PRESS(2,J)/THREE
            PPE( 7,J) = AREA*PRESS(1,J)/THREE
            PPE( 8,J) = AREA*PRESS(2,J)/THREE
            PPE(13,J) = AREA*PRESS(1,J)/THREE
            PPE(14,J) = AREA*PRESS(2,J)/THREE
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate linear differential stiffness matrix

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

        IF (PCOMP_PROPS == 'Y') THEN
          FATAL_ERR          = FATAL_ERR + 1
          NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
          WRITE(ERR,*) ' *ERROR: Code not written for CTRIA3 composite buckling or differential stiffness.'
          WRITE(F06,*) ' *ERROR: Code not written for CTRIA3 composite buckling or differential stiffness.'
          CALL OUTA_HERE ( 'Y' )
        ENDIF

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

        CALL ELMDIS

        DO I=1,ELGP
          DO J=1,ELGP
            KS(I,J) = ZERO
          ENDDO
        ENDDO

        CALL ELEM_STRE_STRN_ARRAYS (1)                     ! Stress at the Gauss point

        FORCEx  = FCONV(1)*STRESS(1)                       ! Engineering forces at the Gauss point
        FORCEy  = FCONV(1)*STRESS(2)
        FORCExy = FCONV(1)*STRESS(3)
        DUM11(1,1) = FORCEx  ; DUM11(1,2) = FORCExy
        DUM11(2,1) = FORCExy ; DUM11(2,2) = FORCEy

        DO I=1,ELGP                                        ! Shape function derivatives at the Gauss point.
          DPSHX(1,I) = BM(1,6*(I-1)+1)
          DPSHX(2,I) = BM(2,6*(I-1)+2)
        ENDDO

        CALL MATMULT_FFF_T ( DPSHX, DUM11, 2, ELGP, 2, DUM12 )
        CALL MATMULT_FFF ( DUM12, DPSHX, ELGP, 2, ELGP, DUM13 )

                                                           ! Accumulate integrand into the result
        DO I=1,ELGP
          DO J=1,ELGP
            KS(I,J) = KS(I,J) + DUM13(I,J) * AREA
          ENDDO
        ENDDO
                                                           ! Copy KS into KED for each translational DOF.
        DO I=1,6*ELGP
          DO J=1,6*ELGP
            KED(I,J) = 0
          ENDDO
        ENDDO
        DO I=1,ELGP
          DO J=1,ELGP
            KED(6*(I-1) + 1,6*(J-1) + 1) = KS(I,J)
            KED(6*(I-1) + 2,6*(J-1) + 2) = KS(I,J)
            KED(6*(I-1) + 3,6*(J-1) + 3) = KS(I,J)
          ENDDO
        ENDDO


      ENDIF




      RETURN

! **********************************************************************************************************************************

 1101 FORMAT(' ------------------------------------------------------------------------------------------------------------------',&
             '-----------------',/,                                                                                                &
             ' ELDATA(',I2,',PRINT) requests for ',A,' element number ',I8,/,                                                      &
             ' ==============================================================',/)

 8901 FORMAT(' Strain-displacement matrix BM for membrane portion of element in subr ',A,/)

 8902 FORMAT(' Row ',I2,/,9(1ES14.6))

 9100 FORMAT(14X,'Check on strain-displacement matrix BM for membrane portion of the element in subr BCHECK'/)

 9101 FORMAT(63X,'S T R A I N S'/,62X,'(direct strains)')

 9102 FORMAT('                                                     Exx            Eyy            Exy')

 9103 FORMAT(7X,'Element displacements consistent with:')

 9104 FORMAT(7X,'---------------------------------------')

! **********************************************************************************************************************************

      END SUBROUTINE TMEM1


      SUBROUTINE TPLT1 ( OPT, AREA, X2E, X3E, Y3E )

! DKT triangular thin (Kirchoff) plate bending element. This element is based on the following work:

! "An Explicit Formulation For An Efficient Triangular Plate-Bending Element", by Jean-Louis Batoz,
! International Journal For Numerical Methods In Engineering, Vol 18 (1982) pp 1655-1677

! Element matrices calculated are:

!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  4) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  5) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y' = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  f06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NSUB, NTSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, FOUR, SIX, TWELVE
      USE MODEL_STUF, ONLY            :  ALPVEC, BE2, DT, EB, KE, PRESS, PPE, PTE, SHELL_DALP, SHELL_D, SHELL_PROP_ALP, SE2, STE2

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TPLT1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices

      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: I1                ! A computed index into array S
      INTEGER(LONG)                   :: I2                ! Part of a computed index into array S
      INTEGER(LONG)                   :: J1                ! A computed index into array KE
      INTEGER(LONG)                   :: K1                ! A computed index into array KE


      REAL(DOUBLE) , INTENT(IN)       :: AREA              ! Element area
      REAL(DOUBLE) , INTENT(IN)       :: X2E               ! x coord of elem node 2
      REAL(DOUBLE) , INTENT(IN)       :: X3E               ! x coord of elem node 3
      REAL(DOUBLE) , INTENT(IN)       :: Y3E               ! y coord of elem node 3
      REAL(DOUBLE)                    :: ALP(3)            ! Col of ALPVEC
      REAL(DOUBLE)                    :: ALPHA(3,3,3,3)    ! 9 sets of 3x3 submatrices needed for calc of elem stiff matrix, KE
      REAL(DOUBLE)                    :: AREAF             ! 48*AREA
      REAL(DOUBLE)                    :: C11               ! Intermediate variable used in calc SEi stress recovery matrices
      REAL(DOUBLE)                    :: C12               ! Intermediate variable used in calc SEi stress recovery matrices
      REAL(DOUBLE)                    :: C21               ! Intermediate variable used in calc SEi stress recovery matrices
      REAL(DOUBLE)                    :: C22               ! Intermediate variable used in calc SEi stress recovery matrices
      REAL(DOUBLE)                    :: C33               ! Intermediate variable used in calc SEi stress recovery matrices
      REAL(DOUBLE)                    :: CT0               ! Intermediate variable used in calc PTE stress recovery matrices
      REAL(DOUBLE)                    :: D1(3,3)           ! Output from subr ATRA, called herein (used in calc KE)
      REAL(DOUBLE)                    :: D2(3,3)           ! Output from subr ATRA, called herein (used in calc KE)
      REAL(DOUBLE)                    :: D3(3,3)           ! Output from subr ATRA, called herein (used in calc KE)
      REAL(DOUBLE)                    :: D4(3,3)           ! Output from subr ATRA, called herein (used in calc KE)
      REAL(DOUBLE)                    :: D5(3,3)           ! Output from subr ATRA, called herein (used in calc KE)
      REAL(DOUBLE)                    :: E1                ! Intermediate variable used in calc KE elem stiffness
      REAL(DOUBLE)                    :: E2                ! Intermediate variable used in calc KE elem stiffness
      REAL(DOUBLE)                    :: E4                ! Intermediate variable used in calc KE elem stiffness
      REAL(DOUBLE)                    :: EALP(3)           ! Intermed var used in calc STEi therm stress coeffs
      REAL(DOUBLE)                    :: L12S              ! Intermediate variable used in calc Pi, Ti, Qi, Ri variables
      REAL(DOUBLE)                    :: L23S              ! Intermediate variable used in calc Pi, Ti, Qi, Ri variables
      REAL(DOUBLE)                    :: L31S              ! Intermediate variable used in calc Pi, Ti, Qi, Ri variables
      REAL(DOUBLE)                    :: P4                ! Intermediate variable used in calc array S
      REAL(DOUBLE)                    :: P5                ! Intermediate variable used in calc array S
      REAL(DOUBLE)                    :: P6                ! Intermediate variable used in calc array S
      REAL(DOUBLE)                    :: Q4                ! Intermediate variable used in calc array S
      REAL(DOUBLE)                    :: Q5                ! Intermediate variable used in calc array S
      REAL(DOUBLE)                    :: R4                ! Intermediate variable used in calc array S
      REAL(DOUBLE)                    :: R5                ! Intermediate variable used in calc array S
      REAL(DOUBLE)                    :: S(3,9)            ! Column sums of the ALPHA-ij submatrices of the 4 dimensional array A
      REAL(DOUBLE)                    :: T4                ! Intermediate variable used in calc array S
      REAL(DOUBLE)                    :: T5                ! Intermediate variable used in calc array S
      REAL(DOUBLE)                    :: X12               ! Diff in x coords of elem nodes 1 and 2
      REAL(DOUBLE)                    :: X23               ! Diff in x coords of elem nodes 2 and 3
      REAL(DOUBLE)                    :: X31               ! Diff in x coords of elem nodes 3 and 1
      REAL(DOUBLE)                    :: Y23               ! Diff in y coords of elem nodes 2 and 3
      REAL(DOUBLE)                    :: Y31               ! Diff in y coords of elem nodes 3 and 1



! **********************************************************************************************************************************
! Generate element parameters

      AREAF = FOUR*TWELVE*AREA
      E1 = SHELL_D(1,1)/AREAF
      E2 = SHELL_D(1,2)/AREAF
      E4 = SHELL_D(3,3)/AREAF

      X12  = -X2E
      X31  =  X3E
      X23  =  X2E - X3E
      Y31  =  Y3E
      Y23  = -Y3E
      L12S =  X12*X12
      L31S =  X31*X31 + Y23*Y23
      L23S =  X23*X23 + Y23*Y23
      P4   = -SIX*X23/L23S
      P5   = -SIX*X3E/L31S
      P6   = -SIX*X12/L12S
      T4   = -SIX*Y23/L23S
      T5   = -SIX*Y3E/L31S
      Q4   =  THREE*X23*Y23/L23S
      Q5   =  THREE*X3E*Y3E/L31S
      R4   =  THREE*Y23*Y23/L23S
      R5   =  THREE*Y31*Y31/L31S

! **********************************************************************************************************************************
! Determine element thermal loads.

      IF (OPT(2) == 'Y') THEN

         DO J=1,NTSUB
            CT0 = SHELL_DALP(1)*DT(4,J)/TWO
            PTE(4,J)  = CT0*X23
            PTE(5,J)  = CT0*Y23
            PTE(10,J) = CT0*X31
            PTE(11,J) = CT0*Y31
            PTE(16,J) = CT0*X12
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Determine element pressure loads.

      IF (OPT(5) == 'Y') THEN
         DO J=1,NSUB
            PPE( 3,J) = AREA*PRESS(3,J)/THREE
            PPE( 9,J) = AREA*PRESS(3,J)/THREE
            PPE(15,J) = AREA*PRESS(3,J)/THREE
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Calculate column sums of the ALPHA-ij submatrices of the 4 dimensional array A. These are needed for SE2 ,STE2 and
! KE calculation. The sums are calculated explicitly instead of summing the ALPHA(i,j,k,l) terms calculated later since
! several terms cancel in the summation.

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y')) THEN

         S(1,1) = Y3E*P5                                   ! ALPHA-11 Column sums
         S(1,2) =-Y3E*Q5
         S(1,3) =-Y3E*R5

         S(2,1) =-X3E*T5                                   ! ALPHA-21 Column sums
         S(2,2) = X3E*R5 + THREE*X23
         S(2,3) =-X3E*Q5

         S(3,1) =-X3E*P5 - X2E*P6 + Y3E*T5                 ! ALPHA-31 Column sums
         S(3,2) = X3E*Q5 + Y3E*(THREE - R5)
         S(3,3) = X3E*R5 + Y3E*Q5

         S(1,4) = Y3E*P4                                   ! ALPHA-12 Column sums
         S(1,5) = Y3E*Q4
         S(1,6) = Y3E*R4

         S(2,4) = X23*T4                                   ! ALPHA-22 Column sums
         S(2,5) = X23*R4 + THREE*X3E
         S(2,6) =-X23*Q4

         S(3,4) = X23*P4 + X2E*P6 + Y3E*T4                 ! ALPHA-32 Column sums
         S(3,5) = X23*Q4 + Y3E*(R4 - THREE)
         S(3,6) = X23*R4 - Y3E*Q4

         S(1,7) =-Y3E*(P4 + P5)                            ! ALPHA-13 Column sums
         S(1,8) = Y3E*(Q4 - Q5)
         S(1,9) = Y3E*(R4 - R5)

         S(2,7) =-X23*T4 + X3E*T5                          ! ALPHA-23 Column sums
         S(2,8) = X23*R4 + X3E*R5 - THREE*X2E
         S(2,9) =-X23*Q4 - X3E*Q5

         S(3,7) =-X23*P4 + X3E*P5 - Y3E*(T4 + T5)          ! ALPHA-33 Column sums
         S(3,8) = X23*Q4 + X3E*Q5 + Y3E*(R4 - R5)
         S(3,9) = X23*R4 + X3E*R5 + Y3E*(Q5 - Q4)

      ENDIF

! **********************************************************************************************************************************
! Calculate SE2, STE2 matrices for stress data recovery.
! Note: stress recovery matrices only make sense for individual plies (or whole elem if only 1 "ply")

      IF (OPT(3) == 'Y') THEN
                                                           ! Strain recovery matrix BE2
         C11 = ONE/(SIX*AREA)
         C12 = ONE/(SIX*AREA)
         C21 = ONE/(SIX*AREA)
         C22 = ONE/(SIX*AREA)
         C33 = ONE/(SIX*AREA)

         BE2(1, 3,1) = (C11*S(1,1) + C12*S(2,1))
         BE2(1, 4,1) = (C11*S(1,2) + C12*S(2,2))
         BE2(1, 5,1) = (C11*S(1,3) + C12*S(2,3))
         BE2(1, 9,1) = (C11*S(1,4) + C12*S(2,4))
         BE2(1,10,1) = (C11*S(1,5) + C12*S(2,5))
         BE2(1,11,1) = (C11*S(1,6) + C12*S(2,6))
         BE2(1,15,1) = (C11*S(1,7) + C12*S(2,7))
         BE2(1,16,1) = (C11*S(1,8) + C12*S(2,8))
         BE2(1,17,1) = (C11*S(1,9) + C12*S(2,9))

         BE2(2, 3,1) = (C21*S(1,1) + C22*S(2,1))
         BE2(2, 4,1) = (C21*S(1,2) + C22*S(2,2))
         BE2(2, 5,1) = (C21*S(1,3) + C22*S(2,3))
         BE2(2, 9,1) = (C21*S(1,4) + C22*S(2,4))
         BE2(2,10,1) = (C21*S(1,5) + C22*S(2,5))
         BE2(2,11,1) = (C21*S(1,6) + C22*S(2,6))
         BE2(2,15,1) = (C21*S(1,7) + C22*S(2,7))
         BE2(2,16,1) = (C21*S(1,8) + C22*S(2,8))
         BE2(2,17,1) = (C21*S(1,9) + C22*S(2,9))

         BE2(3, 3,1) = C33*S(3,1)
         BE2(3, 4,1) = C33*S(3,2)
         BE2(3, 5,1) = C33*S(3,3)
         BE2(3, 9,1) = C33*S(3,4)
         BE2(3,10,1) = C33*S(3,5)
         BE2(3,11,1) = C33*S(3,6)
         BE2(3,15,1) = C33*S(3,7)
         BE2(3,16,1) = C33*S(3,8)
         BE2(3,17,1) = C33*S(3,9)

! SE2, STE2 generated in elem coords. Then, in LINK9 the stresses, calc'd in elem coords, will be transformed to ply coords

         C11 = EB(1,1)/(SIX*AREA)
         C12 = EB(1,2)/(SIX*AREA)
         C21 = EB(2,1)/(SIX*AREA)
         C22 = EB(2,2)/(SIX*AREA)
         C33 = EB(3,3)/(SIX*AREA)

         SE2(1, 3,1) = (C11*S(1,1) + C12*S(2,1))
         SE2(1, 4,1) = (C11*S(1,2) + C12*S(2,2))
         SE2(1, 5,1) = (C11*S(1,3) + C12*S(2,3))
         SE2(1, 9,1) = (C11*S(1,4) + C12*S(2,4))
         SE2(1,10,1) = (C11*S(1,5) + C12*S(2,5))
         SE2(1,11,1) = (C11*S(1,6) + C12*S(2,6))
         SE2(1,15,1) = (C11*S(1,7) + C12*S(2,7))
         SE2(1,16,1) = (C11*S(1,8) + C12*S(2,8))
         SE2(1,17,1) = (C11*S(1,9) + C12*S(2,9))

         SE2(2, 3,1) = (C21*S(1,1) + C22*S(2,1))
         SE2(2, 4,1) = (C21*S(1,2) + C22*S(2,2))
         SE2(2, 5,1) = (C21*S(1,3) + C22*S(2,3))
         SE2(2, 9,1) = (C21*S(1,4) + C22*S(2,4))
         SE2(2,10,1) = (C21*S(1,5) + C22*S(2,5))
         SE2(2,11,1) = (C21*S(1,6) + C22*S(2,6))
         SE2(2,15,1) = (C21*S(1,7) + C22*S(2,7))
         SE2(2,16,1) = (C21*S(1,8) + C22*S(2,8))
         SE2(2,17,1) = (C21*S(1,9) + C22*S(2,9))

         SE2(3, 3,1) = C33*S(3,1)
         SE2(3, 4,1) = C33*S(3,2)
         SE2(3, 5,1) = C33*S(3,3)
         SE2(3, 9,1) = C33*S(3,4)
         SE2(3,10,1) = C33*S(3,5)
         SE2(3,11,1) = C33*S(3,6)
         SE2(3,15,1) = C33*S(3,7)
         SE2(3,16,1) = C33*S(3,8)
         SE2(3,17,1) = C33*S(3,9)

         ALP(1) = ALPVEC(1,1)
         ALP(2) = ALPVEC(2,1)
         ALP(3) = ALPVEC(3,1)

         CALL MATMULT_FFF ( EB, ALP, 3, 3, 1, EALP )
         DO J=1,NTSUB
            DO I=1,3
               STE2(I,J,1) = EALP(I)*DT(4,J)
            ENDDO
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Determine ALPHA-ij matrices needed for KE

      IF(OPT(4) == 'Y') THEN

         ALPHA(1,1,1,1) = Y3E*P6                           ! ALPHA-11, Col 1
         ALPHA(2,1,1,1) =-ALPHA(1,1,1,1)
         ALPHA(3,1,1,1) = Y3E*P5

         ALPHA(1,2,1,1) = ZERO                             !    "      Col 2
         ALPHA(2,2,1,1) = ZERO
         ALPHA(3,2,1,1) =-Y3E*Q5

         ALPHA(1,3,1,1) =-FOUR*Y3E                         !    "      Col 3
         ALPHA(2,3,1,1) = Y3E + Y3E
         ALPHA(3,3,1,1) = Y3E*(TWO - R5)

         ALPHA(1,1,2,1) =-X2E*T5                           ! ALPHA-21, Col 1
         ALPHA(2,1,2,1) = ZERO
         ALPHA(3,1,2,1) = X23*T5

         ALPHA(1,2,2,1) = X23 + X2E*R5                     !    "      Col 2
         ALPHA(2,2,2,1) = X23
         ALPHA(3,2,2,1) = X23*(ONE - R5)

         ALPHA(1,3,2,1) =-X2E*Q5                           !    "      Col 3
         ALPHA(2,3,2,1) = ZERO
         ALPHA(3,3,2,1) = X23*Q5

         ALPHA(1,1,3,1) =-X3E*P6 - X2E*P5                  ! ALPHA-31, Col 1
         ALPHA(2,1,3,1) =-X23*P6
         ALPHA(3,1,3,1) = X23*P5 + Y3E*T5

         ALPHA(1,2,3,1) = X2E*Q5 + Y3E                     !    "      Col 2
         ALPHA(2,2,3,1) = Y3E
         ALPHA(3,2,3,1) = -X23*Q5 + (ONE - R5)*Y3E

         ALPHA(1,3,3,1) = -FOUR*X23 +X2E*R5                !    "      Col 2
         ALPHA(2,3,3,1) = X23 + X23
         ALPHA(3,3,3,1) = (TWO - R5)*X23 + Y3E*Q5

         ALPHA(1,1,1,2) =-Y3E*P6                           ! ALPHA-12, Col 1
         ALPHA(2,1,1,2) = Y3E*P6
         ALPHA(3,1,1,2) = Y3E*P4

         ALPHA(1,2,1,2) = ZERO                             !    "      Col 2
         ALPHA(2,2,1,2) = ZERO
         ALPHA(3,2,1,2) = Y3E*Q4

         ALPHA(1,3,1,2) =-Y3E - Y3E                        !    "      Col 3
         ALPHA(2,3,1,2) = FOUR*Y3E
         ALPHA(3,3,1,2) = Y3E*(R4 - TWO)

         ALPHA(1,1,2,2) = ZERO                             ! ALPHA-22, Col 1
         ALPHA(2,1,2,2) = X2E*T4
         ALPHA(3,1,2,2) =-X3E*T4

         ALPHA(1,2,2,2) = X3E                              !    "      Col 2
         ALPHA(2,2,2,2) = X3E + X2E*R4
         ALPHA(3,2,2,2) = X3E*(ONE - R4)

         ALPHA(1,3,2,2) = ZERO                             !    "      Col 3
         ALPHA(2,3,2,2) =-X2E*Q4
         ALPHA(3,3,2,2) = X3E*Q4

         ALPHA(1,1,3,2) = X3E*P6                           ! ALPHA-32, Col 1
         ALPHA(2,1,3,2) = X23*P6 + X2E*P4
         ALPHA(3,1,3,2) =-X3E*P4 + Y3E*T4

         ALPHA(1,2,3,2) =-Y3E                              !    "      Col 2
         ALPHA(2,2,3,2) = -Y3E + X2E*Q4
         ALPHA(3,2,3,2) = (R4 - ONE)*Y3E - X3E*Q4

         ALPHA(1,3,3,2) = X3E + X3E                        !    "      Col 3
         ALPHA(2,3,3,2) = -FOUR*X3E + X2E*R4
         ALPHA(3,3,3,2) = (TWO - R4)*X3E - Y3E*Q4

         ALPHA(1,1,1,3) = ZERO                             ! ALPHA-13, Col 1
         ALPHA(2,1,1,3) = ZERO
         ALPHA(3,1,1,3) = -Y3E*(P4 + P5)

         ALPHA(1,2,1,3) = ZERO                             !    "      Col 2
         ALPHA(2,2,1,3) = ZERO
         ALPHA(3,2,1,3) = Y3E*(Q4 - Q5)

         ALPHA(1,3,1,3) = ZERO                             !    "      Col 3
         ALPHA(2,3,1,3) = ZERO
         ALPHA(3,3,1,3) = Y3E*(R4 - R5)

         ALPHA(1,1,2,3) = X2E*T5                           ! ALPHA-23, Col 1
         ALPHA(2,1,2,3) =-X2E*T4
         ALPHA(3,1,2,3) =-X23*T5 + X3E*T4

         ALPHA(1,2,2,3) = X2E*(R5 - ONE)                   !    "      Col 2
         ALPHA(2,2,2,3) = X2E*(R4 - ONE)
         ALPHA(3,2,2,3) =-X23*R5 - X3E*R4 - X2E

         ALPHA(1,3,2,3) =-X2E*Q5                           !    "      Col 3
         ALPHA(2,3,2,3) =-X2E*Q4
         ALPHA(3,3,2,3) = X3E*Q4 + X23*Q5

         ALPHA(1,1,3,3) = X2E*P5                           ! ALPHA-33, Col 1
         ALPHA(2,1,3,3) =-X2E*P4
         ALPHA(3,1,3,3) =-X23*P5 + X3E*P4 - Y3E*(T4 + T5)

         ALPHA(1,2,3,3) = X2E*Q5                           !    "      Col 2
         ALPHA(2,2,3,3) = X2E*Q4
         ALPHA(3,2,3,3) =-X23*Q5 - X3E*Q4 + Y3E*(R4 - R5)

         ALPHA(1,3,3,3) = X2E*(R5 - TWO)                   !    "      Col 3
         ALPHA(2,3,3,3) = X2E*(R4 - TWO)
         ALPHA(3,3,3,3) =-X23*R5 - X3E*R4 + FOUR*X2E + Y3E*(Q5 - Q4)

! Calculate the 9 - 3x3 partitions of the element stiffness matrix. Since it is symmetric, only 6 of the 3x3's need to
! be calculated. These 3x3's are put into a global size stiffness matrix for this element which has 18 global DOF.
! The only nonzero's are for DOF's 3,4,5. The resulting 18x18 matrix is in elem coords for the 6 DOF's per grid point.

! Each of the 6 - 3x3's has 5 terms in it. Each of these 5 terms has a triple matrix product consisting of:
!                     T
!           (ALPHA-mi) R (ALPHA-kj)

! where i,j range over the 3 grid points to which the elem connects and m,k are 1,1  2,1  1,2  2,2  3,3 for the 5 terms
! for each i,j pair.

         DO I=1,3

            I1 = 3*I - 2
            I2 = 6*I - 4

            DO J=I,3

               J1 = 6*J - 4
                                                           ! Calculate the 5 triple matrix products using subroutine ATRA (explicit)
               CALL ATRA ( ALPHA(1,1,1,I), ALPHA(1,1,1,J), S(1,I1), S(1,I1+1), S(1,I1+2), D1 )

               CALL ATRA ( ALPHA(1,1,2,I), ALPHA(1,1,1,J), S(2,I1), S(2,I1+1), S(2,I1+2), D2 )

               CALL ATRA ( ALPHA(1,1,2,I), ALPHA(1,1,2,J), S(2,I1), S(2,I1+1), S(2,I1+2), D4 )

               CALL ATRA ( ALPHA(1,1,3,I), ALPHA(1,1,3,J), S(3,I1), S(3,I1+1), S(3,I1+2), D5 )

               IF(I == J) THEN                             ! Calculate the 3x3 element stiffness matrix partition for I,J

                  DO K=1,3
                     K1 = K + I2
                     KE(K1,J1+1) = KE(K1,J1+1) + E1*(D1(K,1)+D4(K,1)) + E2*(D2(K,1)+D2(1,K)) + E4*D5(K,1)
                     KE(K1,J1+2) = KE(K1,J1+2) + E1*(D1(K,2)+D4(K,2)) + E2*(D2(K,2)+D2(2,K)) + E4*D5(K,2)
                     KE(K1,J1+3) = KE(K1,J1+3) + E1*(D1(K,3)+D4(K,3)) + E2*(D2(K,3)+D2(3,K)) + E4*D5(K,3)
                  ENDDO

               ELSE

                  CALL ATRA ( ALPHA(1,1,1,I), ALPHA(1,1,2,J), S(1,I1), S(1,I1+1), S(1,I1+2), D3 )
                  DO K=1,3
                     K1 = K + I2
                     KE(K1,J1+1) = KE(K1,J1+1) + E1*(D1(K,1)+D4(K,1)) + E2*(D2(K,1)+D3(K,1)) + E4*D5(K,1)
                     KE(K1,J1+2) = KE(K1,J1+2) + E1*(D1(K,2)+D4(K,2)) + E2*(D2(K,2)+D3(K,2)) + E4*D5(K,2)
                     KE(K1,J1+3) = KE(K1,J1+3) + E1*(D1(K,3)+D4(K,3)) + E2*(D2(K,3)+D3(K,3)) + E4*D5(K,3)
                  ENDDO

               ENDIF

            ENDDO

         ENDDO

         DO I=2,18                                         ! Calculate sub-diagonal portion by enforcing symmetry
            DO J=1,I-1
               KE(I,J) = KE(J,I)
            ENDDO
         ENDDO

      ENDIF



      RETURN

! **********************************************************************************************************************************

! ######################################################################

      CONTAINS

! ######################################################################

      SUBROUTINE ATRA ( A1, A2, SL1, SL2, SL3, D )

! Subroutine to calculate the triple matrix product, below, needed for the DKT elem stiff matrix:

!                       T
!             (ALPHA-mi) R (ALPHA-kj)

! The product is evaluated explicitly since R is a simple form. R is a 3x3 matrix whose diagonals are all 2.0 and all
! other terms are 1.0

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ATRA'



      REAL(DOUBLE) , INTENT(IN)       :: A1(3,3)           ! ALPHA-mi matrix
      REAL(DOUBLE) , INTENT(IN)       :: A2(3,3)           ! ALPHA-kj matrix
      REAL(DOUBLE) , INTENT(IN)       :: SL1               ! The 3 column sums of ALPHA-mi (k=1 column)
      REAL(DOUBLE) , INTENT(IN)       :: SL2               ! The 3 column sums of ALPHA-mi (k=2 column)
      REAL(DOUBLE) , INTENT(IN)       :: SL3               ! The 3 column sums of ALPHA-mi (k=3 column)
      REAL(DOUBLE) , INTENT(OUT)      :: D(3,3)            ! Triple matrix product discussed above
      REAL(DOUBLE)                    :: W11               ! Intermediate variable used in calculating array D
      REAL(DOUBLE)                    :: W12               ! Intermediate variable used in calculating array D
      REAL(DOUBLE)                    :: W13               ! Intermediate variable used in calculating array D
      REAL(DOUBLE)                    :: W21               ! Intermediate variable used in calculating array D
      REAL(DOUBLE)                    :: W22               ! Intermediate variable used in calculating array D
      REAL(DOUBLE)                    :: W23               ! Intermediate variable used in calculating array D
      REAL(DOUBLE)                    :: W31               ! Intermediate variable used in calculating array D
      REAL(DOUBLE)                    :: W32               ! Intermediate variable used in calculating array D
      REAL(DOUBLE)                    :: W33               ! Intermediate variable used in calculating array D



! **********************************************************************************************************************************
! Wij are the values in ALPHA-mi (transpose) times R

      W11 = A1(1,1) + SL1
      W12 = A1(2,1) + SL1
      W13 = A1(3,1) + SL1

      W21 = A1(1,2) + SL2
      W22 = A1(2,2) + SL2
      W23 = A1(3,2) + SL2

      W31 = A1(1,3) + SL3
      W32 = A1(2,3) + SL3
      W33 = A1(3,3) + SL3

! D is the triple matrix product ALPHA-mi (transp) R ALPHA-kj

      D(1,1) = W11*A2(1,1) + W12*A2(2,1) + W13*A2(3,1)
      D(1,2) = W11*A2(1,2) + W12*A2(2,2) + W13*A2(3,2)
      D(1,3) = W11*A2(1,3) + W12*A2(2,3) + W13*A2(3,3)

      D(2,1) = W21*A2(1,1) + W22*A2(2,1) + W23*A2(3,1)
      D(2,2) = W21*A2(1,2) + W22*A2(2,2) + W23*A2(3,2)
      D(2,3) = W21*A2(1,3) + W22*A2(2,3) + W23*A2(3,3)

      D(3,1) = W31*A2(1,1) + W32*A2(2,1) + W33*A2(3,1)
      D(3,2) = W31*A2(1,2) + W32*A2(2,2) + W33*A2(3,2)
      D(3,3) = W31*A2(1,3) + W32*A2(2,3) + W33*A2(3,3)



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE ATRA

      END SUBROUTINE TPLT1


      SUBROUTINE TPLT2(OPT, AREA, X2E, X3E, Y3E, CALC_EMATS, IERROR, KV, PTV, PPV, B2V, B3V, S2V, S3V, BIG_BB, MN4T_QD,TRIA_NUM,PSI)

! MIN3 triangular thick (Mindlin) plate bending element. This element is based on the following work:

! "A Three-Node Mindlin Plate Element With Improved Transverse Shear", by Alexander Tessler and Thomas J.R. Hughes,
! Computer Methods In Applied Mechanics And Engineering 50 (1985) pp 71-101

! Subroutine calculates:

!  1) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  2) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  3) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  4) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  5) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MEMATC, NSUB, NTSUB
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, THREE, FOUR, SIX, EIGHT, TWELVE, CONV_RAD_DEG
      USE MODEL_STUF, ONLY            :  ALPVEC, BE2, BE3, BENSUM, DT, EB, EBM, EID, ET, ELDOF, FCONV, KE,                         &
                                         MTRL_TYPE, PCOMP_LAM, PCOMP_PROPS, PHI_SQ, PPE, PRESS, PTE, SE2, SE3, SHELL_B, SHELL_DALP,&
                                         SHELL_D, SHELL_T, SHRSUM, STE2, TYPE
      USE PARAMS, ONLY                :  EPSIL, CBMIN3, CBMIN4T
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE VECTOR_GEOMETRY, ONLY       :  PLANE_COORD_TRANS_21
      USE MATERIAL_TRANSFORMATIONS, ONLY:  MATL_TRANSFORM_MATRIX
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T
      USE MIN3_B_MATRICES, ONLY       :  BBMIN3, BSMIN3
      USE SHELL_ROTATION_SUPPORT, ONLY:  CALC_PHI_SQ

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'TPLT2'
      CHARACTER(1*BYTE), INTENT(IN)   :: CALC_EMATS        ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*) , INTENT(IN)   :: MN4T_QD           ! Arg used to say whether the triangular elem is part of a QUAD4

      INTEGER(LONG), INTENT(OUT)      :: IERROR            ! Local error indicator
      INTEGER(LONG), INTENT(IN)       :: TRIA_NUM          ! Tria number (1, 2, 3 or 4) for the subtriangles of a MIN4T QUAD4
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG), PARAMETER        :: ID(9) =   (/ 3, & ! ID(1) =  3 means virgin 9x9 elem DOF 1 is MYSTRAN 18x18 elem DOF  3
                                                      9, & ! ID(2) =  9 means virgin 9x9 elem DOF 2 is MYSTRAN 18x18 elem DOF  9
                                                     15, & ! ID(3) = 15 means virgin 9x9 elem DOF 3 is MYSTRAN 18x18 elem DOF 15
                                                      4, & ! ID(4) =  4 means virgin 9x9 elem DOF 4 is MYSTRAN 18x18 elem DOF  4
                                                     10, & ! ID(5) = 10 means virgin 9x9 elem DOF 5 is MYSTRAN 18x18 elem DOF 10
                                                     16, & ! ID(6) = 16 means virgin 9x9 elem DOF 6 is MYSTRAN 18x18 elem DOF 16
                                                      5, & ! ID(7) =  5 means virgin 9x9 elem DOF 7 is MYSTRAN 18x18 elem DOF  5
                                                     11, & ! ID(8) = 11 means virgin 9x9 elem DOF 8 is MYSTRAN 18x18 elem DOF 11
                                                     17 /) ! ID(9) = 17 means virgin 9x9 elem DOF 9 is MYSTRAN 18x18 elem DOF 17


      REAL(DOUBLE) , INTENT(IN)       :: AREA              ! Element area
      REAL(DOUBLE) , INTENT(IN)       :: PSI               ! Angle to rotate orthotropic mat'l matrix of a sub-tria to align w QUAD
      REAL(DOUBLE) , INTENT(IN)       :: X2E               ! x coord of elem node 2
      REAL(DOUBLE) , INTENT(IN)       :: X3E               ! x coord of elem node 3
      REAL(DOUBLE) , INTENT(IN)       :: Y3E               ! y coord of elem node 3
      REAL(DOUBLE) , INTENT(OUT)      :: BIG_BB(3,18,1)    ! Strain-displ matrix for bending for all DOF's
      REAL(DOUBLE) , INTENT(OUT)      :: B2V(3,9)          ! Strain recovery matrix for virgin DOF's for bending
      REAL(DOUBLE) , INTENT(OUT)      :: B3V(3,9)          ! Strain recovery matrix for virgin DOF's for transverse shear
      REAL(DOUBLE) , INTENT(OUT)      :: KV(9,9)           ! KB + PHISQ*KS (the 9x9 virgin stiffness matrix for MIN3)
      REAL(DOUBLE) , INTENT(OUT)      :: PPV(9,NSUB)       ! The 9xNSUB  virgin thermal  load matrix for MIN3
      REAL(DOUBLE) , INTENT(OUT)      :: PTV(9,NTSUB)      ! The 9xNTSUB virgin pressure load matrix for MIN3
      REAL(DOUBLE) , INTENT(OUT)      :: S2V(3,9)          ! Stress recovery matrix for virgin DOF's for bending
      REAL(DOUBLE) , INTENT(OUT)      :: S3V(3,9)          ! Stress recovery matrix for virgin DOF's for transverse shear
      REAL(DOUBLE)                    :: ALP_TRIA(3)       ! Col of ALPVEC_TRIA

      REAL(DOUBLE)                    :: A(3)              ! 3 diffs in two x coords of the triangle (some x coords are 0)
      REAL(DOUBLE)                    :: A4                ! 4*AREA
      REAL(DOUBLE)                    :: A42               ! 4*AREA**2
      REAL(DOUBLE)                    :: A1(3,3)           ! Intermediate variables used in calc KS (Alex Tessler's matrix Alpha*a)
      REAL(DOUBLE)                    :: A2(3,3)           ! Intermediate variables used in calc KS (Alex Tessler's matrix Alpha*b)
      REAL(DOUBLE)                    :: B(3)              ! 3 diffs in two y coords of the triangle (some y coords are 0)
      REAL(DOUBLE)                    :: B1(3,3)           ! Intermediate variables used in calc KS (Alex Tessler's matrix Beta*a)
      REAL(DOUBLE)                    :: B2(3,3)           ! Intermediate variables used in calc KS (Alex Tessler's matrix Beta*b)
      REAL(DOUBLE)                    :: BB(3,9)           ! Bending strain-displ matrix for the MIN3 elem (from subr BBMIN3)
      REAL(DOUBLE)                    :: BS(2,9)           ! Shear strain-displ matrix for the MIN3 elem (from subr BSMIN3)
      REAL(DOUBLE)                    :: DUM0(9)           ! Intermediate variables used in calc PTE, PPE (thermal, pressure loads)
      REAL(DOUBLE)                    :: DUM31(3)          ! Intermadiate result in calc KE elem stiffness
      REAL(DOUBLE)                    :: DUM1(3,3)         ! Intermadiate result in calc KE elem stiffness
      REAL(DOUBLE)                    :: DUM2(3,3)         ! Intermadiate result in calc KE elem stiffness
      REAL(DOUBLE)                    :: DUM3(3,3)         ! Intermadiate result in calc KE elem stiffness
      REAL(DOUBLE)                    :: DUM4(3,3)         ! Intermadiate result in calc KE elem stiffness
      REAL(DOUBLE)                    :: DUM5(3,9)         ! Intermadiate result in calc SEi stress recovery matrices
      REAL(DOUBLE)                    :: DUM6(2,9)         ! Intermadiate result in calc SEi stress recovery matrices
      REAL(DOUBLE)                    :: DUM22(2,2)        ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM33(3,3)        ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM64(6,4)        ! Intermediate matrix
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare to real zero
      REAL(DOUBLE)                    :: FXX(3,3)          ! Intermadiate result in calc BB (Alex Tessler matrix fxx)
      REAL(DOUBLE)                    :: FXY(3,3)          ! Intermadiate result in calc BB (Alex Tessler matrix fxy)
      REAL(DOUBLE)                    :: FYY(3,3)          ! Intermadiate result in calc BB (Alex Tessler matrix fyy)
      REAL(DOUBLE)                    :: I00(3,3)          ! Intermadiate result in calc KS shear stiffness
      REAL(DOUBLE)                    :: IX0(3,3)          ! Intermadiate result in calc KS shear stiffness
      REAL(DOUBLE)                    :: IY0(3,3)          ! Intermadiate result in calc KS shear stiffness
      REAL(DOUBLE)                    :: KB(9,9)           ! Bending stiffness contribution to KE for the MIN3 elem
      REAL(DOUBLE)                    :: KS(9,9)           ! PHISQ*KS is the shear stiff contribution to KE for the MIN3 elem
      REAL(DOUBLE)                    :: QCONS             ! = AREA/24, used in calc PPE pressure loads
      REAL(DOUBLE)                    :: S1(3,3)           ! Intermadiate result in calc shear stiffness, KS
      REAL(DOUBLE)                    :: S2(3,3)           ! Intermadiate result in calc shear stiffness, KS
      REAL(DOUBLE)                    :: T1(3,3)           ! Intermadiate result in calc shear stiffness, KS
      REAL(DOUBLE)                    :: T2(3,3)           ! Intermadiate result in calc shear stiffness, KS
      REAL(DOUBLE)                    :: T3(3,3)           ! Intermadiate result in calc shear stiffness, KS
      REAL(DOUBLE)                    :: T4(3,3)           ! Intermadiate result in calc shear stiffness, KS
      REAL(DOUBLE)                    :: TF(6,6)           ! Transforms 6 stress and 6x6 material mats from material to element axes
!                                                            (TF' transforms strains)
      REAL(DOUBLE)                    :: TME(3,3)          ! Coord transf matrix which will rotate a vector in local element coord
!                                                            system to a vector in the MIN4T QUAD4 matl coord system (Um = TME*Ue)
      REAL(DOUBLE)                    :: TF_MB(3,3)        ! Portion of TF: transforms 3x3 EM, EB, EBM from material to elem axes
      REAL(DOUBLE)                    :: TF_TS(2,2)        ! Portion of TF: transforms 3x3 ET from material to elem axes
      REAL(DOUBLE)                    :: XI(3)

! Following  are matl matrices used when sub-trias of a MIN4T QUAD4 need to be transformed to align with orthotropic matl angles

      REAL(DOUBLE)                    :: ALPVEC_TRIA(6,MEMATC)
      REAL(DOUBLE)                    :: EALP_TRIA(3)      ! Intermed var used in calc STEi therm stress coeffs
      REAL(DOUBLE)                    :: EB0(3,3)          ! Plane stress matl matrix for bending before coord transformation
      REAL(DOUBLE)                    :: EBM0(3,3)         ! Bend/membr coupling matl matrix before coord transformation
      REAL(DOUBLE)                    :: ET0(2,2)          ! 2D transverse shear matl matrix before coord transformation
      REAL(DOUBLE)                    :: EB_TRIA(3,3)      ! Plane stress matl matrix for bending after coord transformation
      REAL(DOUBLE)                    :: EBM_TRIA(3,3)     ! Bend/membr coupling matl matrix before after transformation
      REAL(DOUBLE)                    :: ET_TRIA(2,2)      ! 2D transverse shear matl matrix before after transformation

      REAL(DOUBLE)                    :: SHELL_B0_TRIA(3,3)! SHELL_B_TRIA before coord transformation
      REAL(DOUBLE)                    :: SHELL_D0_TRIA(3,3)! SHELL_D_TRIA before coord transformation
      REAL(DOUBLE)                    :: SHELL_T0_TRIA(2,2)! SHELL_T_TRIA before coord transformation
      REAL(DOUBLE)                    :: SHELL_B_TRIA(3,3) ! SHELL_B_TRIA after  coord transformation
      REAL(DOUBLE)                    :: SHELL_D_TRIA(3,3) ! SHELL_D_TRIA after  coord transformation
      REAL(DOUBLE)                    :: SHELL_T_TRIA(2,2) ! SHELL_T_TRIA after  coord transformation

                                                           ! SHELL_D_TRIA before coord transformation
      REAL(DOUBLE)                    :: SHELL_DALP0_TRIA(3)

                                                           ! SHELL_D_TRIA after  coord transformation
      REAL(DOUBLE)                    :: SHELL_DALP_TRIA(3)

      INTRINSIC DABS



! **********************************************************************************************************************************
! Initialize

      EPS1   = EPSIL(1)

      IERROR = 0

      DO I=1,9
         DO J=1,9
            KV(I,J) = ZERO
         ENDDO
         DO J=1,NSUB
            PPV(I,J) = ZERO
         ENDDO
         DO J=1,NTSUB
            PTV(I,J) = ZERO
         ENDDO
      ENDDO

! **********************************************************************************************************************************
! If this subr is being called because the triangular shell element is part of a MIN4T QUAD4 with orthotropic mat'l properties then
! we need to re-orient those ortho material properties to be aligned with the quad.

! Set initial values before coord transform (NOTE: Only needed for DEBUG(53)

      DO I=1,3
         SHELL_DALP0_TRIA(I) = SHELL_DALP(I)
         DO J=1,3
            EB0(I,J)           = EB(I,J)
            SHELL_D0_TRIA(I,J) = SHELL_D(I,J)
            EBM0(I,J)          = EBM(I,J)
            SHELL_B0_TRIA(I,J) = SHELL_B(I,J)
         ENDDO
      ENDDO

      DO I=1,2
         DO J=1,2
            ET0(I,J)           = ET(I,J)
            SHELL_T0_TRIA(I,J) = SHELL_T(I,J)
         ENDDO
      ENDDO

      DO I=1,3
         SHELL_DALP_TRIA(I) = SHELL_DALP(I)
         DO J=1,3
            EB_TRIA(I,J)      = EB(I,J)
            SHELL_D_TRIA(I,J) = SHELL_D(I,J)
            EBM_TRIA(I,J)     = EBM(I,J)
            SHELL_B_TRIA(I,J) = SHELL_B(I,J)
         ENDDO
      ENDDO

      DO I=1,2
         DO J=1,2
            ET_TRIA(I,J)      = ET(I,J)
            SHELL_T_TRIA(I,J) = SHELL_T(I,J)
         ENDDO
      ENDDO

      DO I=1,6
         DO J=1,MEMATC
            ALPVEC_TRIA(I,J) = ALPVEC(I,J)
         ENDDO
      ENDDO

                                                           ! The following 3 conditions have to be met before we look at mat'l props
      IF ((MN4T_QD == 'Y') .AND. (TRIA_NUM >= 1) .AND. (TRIA_NUM <= 4)) THEN
                                                           ! If either bending or transverse shear props are ortho we need TF matrix
         CBMIN4T = CBMIN3
         IF ((MTRL_TYPE(2) == 8) .OR. (MTRL_TYPE(3) == 8)) THEN

            CALL PLANE_COORD_TRANS_21 ( PSI, TME, SUBR_NAME )
            CALL MATL_TRANSFORM_MATRIX ( TME, TF )
                                                        ! TF_MB is for Sxx, Syy, Sxy which are rows and cols 1,2,4 from TF
            TF_MB(1,1) = TF(1,1)     ;     TF_MB(1,2) = TF(1,2)     ;     TF_MB(1,3) = TF(1,4)
            TF_MB(2,1) = TF(2,1)     ;     TF_MB(2,2) = TF(2,2)     ;     TF_MB(2,3) = TF(2,4)
            TF_MB(3,1) = TF(4,1)     ;     TF_MB(3,2) = TF(4,2)     ;     TF_MB(3,3) = TF(4,4)

                                                        ! TF_ET is for Syz, Szx which are rows and cols 5,6 from TF
            TF_TS(1,1) = TF(5,5)     ;     TF_TS(1,2) = TF(5,6)
            TF_TS(2,1) = TF(6,5)     ;     TF_TS(2,2) = TF(6,6)

!           ------------------------------------------------------------------------------------------------------------------------
            IF ((MTRL_TYPE(2) == 2) .OR. (MTRL_TYPE(2) == 8)) THEN   ! Transform bending material matrix
               CALL MATMULT_FFF   ( EB_TRIA , TF_MB  , 3, 3, 3, DUM33)
               CALL MATMULT_FFF_T ( TF_MB   , DUM33  , 3, 3, 3, EB_TRIA)
            ENDIF

            IF ((MTRL_TYPE(2) == 2) .OR. (MTRL_TYPE(2) == 8)) THEN   ! Transform SHELL_D matrix
               CALL MATMULT_FFF   ( SHELL_D_TRIA , TF_MB  , 3, 3, 3, DUM33)
               CALL MATMULT_FFF_T ( TF_MB        , DUM33  , 3, 3, 3, SHELL_D_TRIA)
            ENDIF
!           ------------------------------------------------------------------------------------------------------------------------
            IF ((MTRL_TYPE(3) == 2) .OR. (MTRL_TYPE(3) == 8)) THEN   ! Transform transverse shear material matrix
               CALL MATMULT_FFF   ( ET_TRIA , TF_TS  , 2, 2, 2, DUM22)
               CALL MATMULT_FFF_T ( TF_TS   , DUM22  , 2, 2, 2, ET_TRIA)
            ENDIF

            IF ((MTRL_TYPE(3) == 2) .OR. (MTRL_TYPE(3) == 8)) THEN   ! Transform SHELL_T matrix
               CALL MATMULT_FFF   ( SHELL_T_TRIA , TF_TS  , 2, 2, 2, DUM22)
               CALL MATMULT_FFF_T ( TF_TS        , DUM22  , 2, 2, 2, SHELL_T_TRIA)
            ENDIF
!           ------------------------------------------------------------------------------------------------------------------------
            IF ((MTRL_TYPE(4) == 2) .OR. (MTRL_TYPE(4) == 8)) THEN   ! Transform coupled bending/membrane material matrix
               CALL MATMULT_FFF   ( EBM_TRIA, TF_MB  , 3, 3, 3, DUM33)
               CALL MATMULT_FFF_T ( TF_MB   , DUM33  , 3, 3, 3, EBM_TRIA)
            ENDIF

            IF ((MTRL_TYPE(4) == 2) .OR. (MTRL_TYPE(4) == 8)) THEN   ! Transform SHELL_B matrix
               CALL MATMULT_FFF   ( SHELL_B_TRIA , TF_MB  , 3, 3, 3, DUM33)
               CALL MATMULT_FFF_T ( TF_MB        , DUM33  , 3, 3, 3, SHELL_B_TRIA)
            ENDIF
!           ------------------------------------------------------------------------------------------------------------------------

            CALL MATMULT_FFF_T (TF, ALPVEC_TRIA, 6,6, MEMATC, DUM64) ! Transform CTE matrix
            DO I=1,6
               DO J=1,MEMATC
                  ALPVEC_TRIA(I,J) = DUM64(I,J)
               ENDDO
            ENDDO

            CALL MATMULT_FFF_T (TF_MB, SHELL_DALP_TRIA, 3, 3, 1, DUM31) ! Transform SHELL_DALP matrix
            DO I=1,3
               SHELL_DALP_TRIA(I) = DUM31(I)
            ENDDO
!           ------------------------------------------------------------------------------------------------------------------------

            IF (DEBUG(53) > 0) THEN
               CALL DEBUG_ROT_AXES_2
            ENDIF

         ENDIF

      ENDIF

! **********************************************************************************************************************************
      A(1) =  X3E - X2E
      A(2) = -X3E
      A(3) =  X2E
      B(1) = -Y3E
      B(2) =  Y3E
      B(3) =  ZERO

      A4  = FOUR*AREA
      A42 = A4*AREA
! BB is used in several places below:

      CALL BBMIN3 ( A, B, AREA, '(bending strains)', 'Y', BB )

! **********************************************************************************************************************************
! Determine element thermal loads. Code is only valid for materials with ALPHA-12 = 0.

      IF (OPT(2) == 'Y') THEN

         CALL MATMULT_FFF_T ( BB, SHELL_DALP, 3, 9, 1, DUM0 )
         DO J=1,NTSUB
            DO I=1,9
               PTV(I,J) = AREA*DUM0(I)*DT(4,J)
               IF (CALC_EMATS == 'Y') THEN
                  PTE(ID(I),J) = PTV(I,J)
               ENDIF
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE. Note that we need to calc KE if the stress recovery matrices (for OPT(3)) are to be cal'd
! since the BE strain recovery matrices use PHI_SQ

      IF ((OPT(4) == 'Y') .OR. (OPT(3) == 'Y')) THEN

         DO I=1,9
            DO J=1,9
               KB(I,J) = ZERO
            ENDDO
         ENDDO

! Bending stiffness terms (KB)

         DO I=1,3
            DO J=1,3
               FXX(I,J) = B(I)*B(J)/A42
               FYY(I,J) = A(I)*A(J)/A42
               FXY(I,J) = B(I)*A(J)/A42
            ENDDO
         ENDDO
         DO I=1,3
            DO J=1,3
               KB(I+3,J+3)= AREA*(SHELL_D_TRIA(2,2)*FYY(I,J) + SHELL_D_TRIA(2,3)*(FXY(I,J) + FXY(J,I)) + SHELL_D_TRIA(3,3)*FXX(I,J))
            ENDDO
         ENDDO

         DO I=1,3
            DO J=1,3
               KB(I+3,J+6) = -AREA*(SHELL_D_TRIA(1,2)*FXY(J,I) + SHELL_D_TRIA(2,3)*FYY(I,J) + SHELL_D_TRIA(1,3)*FXX(I,J) +         &
                                    SHELL_D_TRIA(3,3)*FXY(I,J))
               KB(J+6,I+3) = KB(I+3,J+6)
            ENDDO
         ENDDO

         DO I=1,3
            DO J=1,3
               KB(I+6,J+6)= AREA*(SHELL_D_TRIA(1,1)*FXX(I,J) + SHELL_D_TRIA(1,3)*(FXY(I,J) + FXY(J,I)) + SHELL_D_TRIA(3,3)*FYY(I,J))
            ENDDO
         ENDDO

         BENSUM = ZERO
         DO I=4,9
            BENSUM = BENSUM + KB(I,I)
         ENDDO
! Shear stress terms (KS)

         DO I=1,9
            DO J=1,9
               KS(I,J) = ZERO
            ENDDO
         ENDDO

         B2(1,1) =   ZERO
         B2(2,2) =   ZERO
         B2(3,3) =   ZERO
         B2(1,2) =  -B(2)*B(3)/A4
         B2(1,3) =  -B2(1,2)
         B2(2,1) =   B(1)*B(3)/A4
         B2(2,3) =  -B2(2,1)
         B2(3,1) =  -B(1)*B(2)/A4
         B2(3,2) =  -B2(3,1)

         A2(1,1) = ( A(2)*B(3) - A(3)*B(2))/A4
         A2(2,2) = (-A(1)*B(3) + A(3)*B(1))/A4
         A2(3,3) = (-A(2)*B(1) + A(1)*B(2))/A4
         A2(1,2) =  -A(2)*B(3)/A4
         A2(1,3) =   A(3)*B(2)/A4
         A2(2,1) =   A(1)*B(3)/A4
         A2(2,3) =  -A(3)*B(1)/A4
         A2(3,1) =  -A(1)*B(2)/A4
         A2(3,2) =   A(2)*B(1)/A4

         B1(1,1) = (-B(2)*A(3) + B(3)*A(2))/A4
         B1(2,2) = ( B(1)*A(3) - B(3)*A(1))/A4
         B1(3,3) = ( B(2)*A(1) - B(1)*A(2))/A4
         B1(1,2) =   B(2)*A(3)/A4
         B1(1,3) =  -B(3)*A(2)/A4
         B1(2,1) =  -B(1)*A(3)/A4
         B1(2,3) =   B(3)*A(1)/A4
         B1(3,1) =   B(1)*A(2)/A4
         B1(3,2) =  -B(2)*A(1)/A4

         A1(1,1) =   ZERO
         A1(2,2) =   ZERO
         A1(3,3) =   ZERO
         A1(1,2) =   A(2)*A(3)/A4
         A1(1,3) =  -A1(1,2)
         A1(2,1) =  -A(1)*A(3)/A4
         A1(2,3) =  -A1(2,1)
         A1(3,1) =   A(1)*A(2)/A4
         A1(3,2) =  -A1(3,1)

         I00(1,2) = AREA/TWELVE
         I00(1,3) = I00(1,2)
         I00(2,1) = I00(1,2)
         I00(2,3) = I00(1,2)
         I00(3,1) = I00(1,2)
         I00(3,2) = I00(1,2)
         I00(1,1) = TWO*I00(1,2)
         I00(2,2) = I00(1,1)
         I00(3,3) = I00(1,1)

         IX0(1,1) = B(1)/SIX
         IX0(1,2) = IX0(1,1)
         IX0(1,3) = IX0(1,1)
         IX0(2,1) = B(2)/SIX
         IX0(2,2) = IX0(2,1)
         IX0(2,3) = IX0(2,1)
         IX0(3,1) = B(3)/SIX
         IX0(3,2) = IX0(3,1)
         IX0(3,3) = IX0(3,1)

         IY0(1,1) = A(1)/SIX
         IY0(1,2) = IY0(1,1)
         IY0(1,3) = IY0(1,1)
         IY0(2,1) = A(2)/SIX
         IY0(2,2) = IY0(2,1)
         IY0(2,3) = IY0(2,1)
         IY0(3,1) = A(3)/SIX
         IY0(3,2) = IY0(3,1)
         IY0(3,3) = IY0(3,1)

         DO I=1,3
            DO J=1,3
               IF (I == J) THEN
                  S1(I,J) = A2(I,J) + ONE
                  S2(I,J) = B1(I,J) + ONE
               ELSE
                  S1(I,J) = A2(I,J)
                  S2(I,J) = B1(I,J)
               ENDIF
            ENDDO
         ENDDO

         DO I=1,3
            DO J=1,3
               T1(I,J) = (SHELL_T_TRIA(1,1)*B2(I,J) + SHELL_T_TRIA(1,2)*S1(I,J))
               T2(I,J) = (SHELL_T_TRIA(1,2)*B2(I,J) + SHELL_T_TRIA(2,2)*S1(I,J))
               T3(I,J) = (SHELL_T_TRIA(1,1)*S2(I,J) + SHELL_T_TRIA(1,2)*A1(I,J))
               T4(I,J) = (SHELL_T_TRIA(1,2)*S2(I,J) + SHELL_T_TRIA(2,2)*A1(I,J))
            ENDDO
         ENDDO
                                                        ! 3x3 KS-11 Partition
         DO I=1,3
            DO J=1,3
               KS(I,J) = AREA*(SHELL_T_TRIA(1,1)*FXX(I,J) + SHELL_T_TRIA(1,2)*(FXY(I,J) + FXY(J,I)) + SHELL_T_TRIA(2,2)*FYY(I,J))
            ENDDO
         ENDDO
                                                        ! 3x3 KS-12, 21 Partition
         CALL MATMULT_FFF ( IX0, T1, 3, 3, 3, DUM1 )
         CALL MATMULT_FFF ( IY0, T2, 3, 3, 3, DUM2 )
         DO I=1,3
            DO J=1,3
               KS(I  ,J+3) = -DUM1(I,J) - DUM2(I,J)
               KS(J+3,I  ) = KS(I,J+3)
            ENDDO
         ENDDO
                                                        ! 3x3 KS-13, 31 Partition
         CALL MATMULT_FFF ( IX0, T3, 3, 3, 3, DUM1 )
         CALL MATMULT_FFF ( IY0, T4, 3, 3, 3, DUM2 )
         DO I=1,3
            DO J=1,3
               KS(I  ,J+6) = DUM1(I,J) + DUM2(I,J)
               KS(J+6,I  ) = KS(I,J+6)
            ENDDO
         ENDDO
                                                        ! 3x3 KS-22 Partition
         CALL MATMULT_FFF   ( I00, T1  , 3, 3, 3, DUM3 )
         CALL MATMULT_FFF_T ( B2 , DUM3, 3, 3, 3, DUM1 )
         CALL MATMULT_FFF   ( I00, T2  , 3, 3, 3, DUM4 )
         CALL MATMULT_FFF_T ( S1 , DUM4, 3, 3, 3, DUM2 )
         DO I=1,3
            DO J=1,3
               KS(I+3,J+3) = DUM1(I,J) + DUM2(I,J)
            ENDDO
         ENDDO
                                                        ! 3x3 KS-23, 32 Partition
         CALL MATMULT_FFF   ( I00, T3  , 3, 3, 3, DUM3 )
         CALL MATMULT_FFF_T ( B2 , DUM3, 3, 3, 3, DUM1 )
         CALL MATMULT_FFF   ( I00, T4  , 3, 3, 3, DUM4 )
         CALL MATMULT_FFF_T ( S1 , DUM4, 3, 3, 3, DUM2 )
         DO I=1,3
            DO J=1,3
               KS(I+3,J+6) = -DUM1(I,J) - DUM2(I,J)
               KS(J+6,I+3) = KS(I+3,J+6)
            ENDDO
         ENDDO
                                                        ! 3x3 KS-33 Partition
         CALL MATMULT_FFF   ( I00, T3  , 3, 3, 3, DUM3 )
         CALL MATMULT_FFF_T ( S2 , DUM3, 3, 3, 3, DUM1 )
         CALL MATMULT_FFF   ( I00, T4  , 3, 3, 3, DUM4 )
         CALL MATMULT_FFF_T ( A1 , DUM4, 3, 3, 3, DUM2 )
         DO I=1,3
            DO J=1,3
               KS(I+6,J+6) = DUM1(I,J) + DUM2(I,J)
            ENDDO
         ENDDO

         SHRSUM = ZERO
         DO I=4,9
            SHRSUM = SHRSUM + KS(I,I)
         ENDDO
  ! Now calculate the finite elem shear factor, PHI_SQ

         IF (SHRSUM > EPS1) THEN
            CALL CALC_PHI_SQ ( IERROR )
         ELSE
            PHI_SQ = ZERO
         ENDIF
! Return if IERROR > 0

         IF (IERROR > 0) RETURN

         DO I=1,9
            DO J=1,9
               KV(I,J) = KB(I,J) + PHI_SQ*KS(I,J)
               IF (CALC_EMATS == 'Y') THEN
                  KE(ID(I),ID(J)) = KE(ID(I),ID(J)) + KV(I,J)
               ENDIF
            ENDDO
         ENDDO
         IF (DEBUG(54) == 1) THEN
            IF (MN4T_QD == 'Y') THEN
               WRITE(F06,'(A,I2,A,A,I8))') ' KV, in TPLT2 for MIN4T QUAD4: tria ', TRIA_NUM, ' of 4 for ',TRIM(TYPE), EID
            ELSE
               WRITE(F06,'(A,I8)') ' KV in TPLT2 for TRIA3 ', EID
            ENDIF
            WRITE(F06,*) '    Column  3      Column  4      Column  5      Column  9      Column 10      Column 11      Column 15',&
                       '      Column 16      Column 17'
            DO I=1,9
               WRITE(F06,'(9ES15.6)') (KV(I,J),J=1,9)
            ENDDO
            WRITE(F06,*)
         ENDIF

! Set lower triangular partition equal to upper partition
         IF (CALC_EMATS == 'Y') THEN
            DO I=2,ELDOF
               DO J=1,I-1
                  KE(I,J) = KE(J,I)
               ENDDO
            ENDDO
         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Calculate BE2, SE2 matrix (3 x 18) for strain/stress data recovery.
! Note: stress recovery matrices only make sense for individual plies (or whole elem if only 1 "ply")

! BS (transverse shear strain-displ matrices) is used in several places below:

      CALL BBMIN3 ( A, B, AREA, '(bending strains)', 'N', BB )

      XI(1) = ONE/THREE
      XI(2) = ONE/THREE
      XI(3) = ONE/THREE
      CALL BSMIN3 ( XI, A, B, AREA, '(transverse shear strains)', 'Y', BS )

! Calc BEi matrices regardless of OPT

      DO I=1,3
         DO J=1,9
            B2V(I,J) = ZERO
            B3V(I,J) = ZERO
         ENDDO
      ENDDO

      DO I=1,3
         DO J=1,9
            B2V(I,J)       = BB(I,J)
            BE2(I,ID(J),1) = BB(I,J)
         ENDDO
      ENDDO

      DO I=1,2
         DO J=1,9
            B3V(I,J)       = BS(I,J)
            BE3(I,ID(J),1) = BS(I,J)
         ENDDO
      ENDDO

! SE2, STE2 generated in elem coords. Then, in LINK9 the stresses, calc'd in elem coords, will be transformed to ply coords

      IF (OPT(3) == 'Y') THEN

         DO I=1,3
            DO J=1,9
               S2V(I,J) = ZERO
               S3V(I,J) = ZERO
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( EB_TRIA, BB, 3, 3, 9, DUM5 )   ! Generate SE2 in element coords (at this point EB is elem coords)
         DO I=1,3
            DO J=1,9
               S2V(I,J)       = DUM5(I,J)
               SE2(I,ID(J),1) = DUM5(I,J)
            ENDDO
         ENDDO

         ALP_TRIA(1) = ALPVEC_TRIA(1,1)
         ALP_TRIA(2) = ALPVEC_TRIA(2,1)
         ALP_TRIA(3) = ALPVEC_TRIA(3,1)

         CALL MATMULT_FFF ( EB_TRIA, ALP_TRIA, 3, 3, 1, EALP_TRIA )  ! Mult EB (ply coords) times ALP (ply coords)

         DO J=1,NTSUB                                      ! Thermal stress terms
            DO I=1,3
               STE2(I,J,1) = EALP_TRIA(I)*DT(4,J)
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( ET_TRIA, BS, 2, 2, 9, DUM6 )
         DO I=1,2
            DO J=1,9
               S3V(I,J)       = PHI_SQ*DUM6(I,J)
               SE3(I,ID(J),1) = PHI_SQ*DUM6(I,J)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! If element is a composite and if it is a nonsym layup we need to calc BIG_BB for later use

      DO I=1,3
         DO J=1,18
            BIG_BB(I,J,1) = ZERO
         ENDDO
      ENDDO

      IF ((PCOMP_PROPS == 'Y') .AND. (PCOMP_LAM == 'NON')) THEN

         DO I=1,3
            DO J=1,9
               BIG_BB(I,ID(J),1) = BB(I,J)
            ENDDO
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Determine element pressure loads.

      IF (OPT(5) == 'Y') THEN

         IF (DEBUG(16) == 0) THEN                          ! Generate PPE as work equilavent loads

            QCONS   = AREA/(TWO*TWELVE)
            DUM0(1) = QCONS*EIGHT
            DUM0(2) = QCONS*EIGHT
            DUM0(3) = QCONS*EIGHT
            DUM0(4) = QCONS*(-B(3) + B(2))
            DUM0(5) = QCONS*( B(3) - B(1))
            DUM0(6) = QCONS*( B(1) - B(2))
            DUM0(7) = QCONS*(-A(3) + A(2))
            DUM0(8) = QCONS*( A(3) - A(1))
            DUM0(9) = QCONS*( A(1) - A(2))

            DO J=1,NSUB
               DO I=1,9
                  PPV(I,J)     = PRESS(3,J)*DUM0(I)
                  IF (CALC_EMATS == 'Y') THEN
                     PPE(ID(I),J) = PPV(I,J)
                  ENDIF
               ENDDO
            ENDDO

         ELSE                                              ! Generate PPE as static equilavent loads

            DO J=1,NSUB
               PPV( 1,J) = AREA*PRESS(3,J)/THREE
               PPV( 2,J) = AREA*PRESS(3,J)/THREE
               PPV( 3,J) = AREA*PRESS(3,J)/THREE
               IF (CALC_EMATS == 'Y') THEN
                  PPE( 3,J) = AREA*PRESS(3,J)/THREE
                  PPE( 9,J) = AREA*PRESS(3,J)/THREE
                  PPE(15,J) = AREA*PRESS(3,J)/THREE
               ENDIF
            ENDDO

         ENDIF

      ENDIF



      RETURN

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEBUG_ROT_AXES_2

      USE IOUNT1, ONLY                :  F06
      USE MODEL_STUF, ONLY            :  EID

      IMPLICIT NONE

      CHARACTER(LEN=100)              :: MI(9)             ! Messages to be written out

! **********************************************************************************************************************************
      DO I=1,9
         MI(I)(1:) = ' '
      ENDDO

! Write outputs

      WRITE(F06,98720)

      WRITE(F06,'(A,I2,A,I8,A,1ES14.6,A)') ' Data for TRIA_NUM ', tria_num,' of MIN4T QUAD4 elem ', eid, ' which is to be rotated',&
                                          conv_rad_deg*psi,' deg'
      WRITE(F06,55566)

      MI(1) = ' Transforms 6 stress and 6x6 matl matrices from tria to quad axes (TF transpose transforms strains)'
      WRITE(F06,99664) 'TF', MI(1)
      DO I=1,3
         WRITE(F06,99667) (TF(I,J),J=1,6)
      ENDDO
      WRITE(F06,*)
      DO I=4,6
         WRITE(F06,99667) (TF(I,J),J=1,6)
      ENDDO
      WRITE(F06,*)

      MI(2) = ' Portion of TF: transforms 3x3 EB_TRIA, EBM_TRIA from tria to quad axes'
      WRITE(F06,99664) 'TF_MB', MI(2)
      DO I=1,3
         WRITE(F06,99668) (TF_MB(I,J),J=1,3)
      ENDDO
      WRITE(F06,*)

      MI(3) = '  Portion of TF: transforms 2x2 ET_TRIA from tria to quad axes'
      WRITE(F06,99664) 'TF_TS', MI(3)
      DO I=1,2
         WRITE(F06,99669) (TF_TS(I,J),J=1,2)
      ENDDO
      WRITE(F06,*)

      MI(4) = ' Coord transf matrix which will rotate a vector in local elem'
      MI(5) = ' coord system to a vector in the material coord system (Um = TME*Ue)'
      WRITE(F06,99664) 'TME', MI(4)   ;   WRITE(F06,99663) MI(5)
      DO I=1,3
         WRITE(F06,99668) (TME(I,J),J=1,3)
      ENDDO
      WRITE(F06,*)


bend: IF ((MTRL_TYPE(2) == 2) .OR. (MTRL_TYPE(2) == 8)) THEN

         WRITE(F06,99665) 'EB_TRIA '
         DO I=1,3
            WRITE(F06,99668) (EB0(I,J),J=1,3), (EB_TRIA(I,J),J=1,3)
         ENDDO
         WRITE(F06,*)   ;   WRITE(F06,*)

         WRITE(F06,99665) 'SHELL_D_TRIA'
         DO I=1,3
            WRITE(F06,99668) (SHELL_D0_TRIA(I,J),J=1,3), (SHELL_D_TRIA(I,J),J=1,3)
         ENDDO
         WRITE(F06,*)

      ENDIF bend

tshr: IF ((MTRL_TYPE(3) == 2) .OR. (MTRL_TYPE(3) == 8)) THEN

         write(f06,99665) 'ET_TRIA '
         do i=1,2
            write(f06,99669) (et0(i,j),j=1,2), (et_tria(i,j),j=1,2)
         enddo
         write(f06,*)   ;   write(f06,*)

         write(f06,99665) 'SHELL_T_TRIA'
         DO I=1,2
            WRITE(F06,99669) (SHELL_T0_TRIA(I,J),J=1,2), (SHELL_T_TRIA(I,J),J=1,2)
         ENDDO
         WRITE(F06,*)

      ENDIF tshr

mbend:IF ((MTRL_TYPE(4) == 2) .OR. (MTRL_TYPE(4) == 8)) THEN

         WRITE(F06,99665) 'EBM_TRIA'
         DO I=1,3
            WRITE(F06,99668) (EBM0(I,J),J=1,3), (EBM_TRIA(I,J),J=1,3)
         ENDDO
         WRITE(F06,*)   ;   WRITE(F06,*)

         WRITE(F06,99665) 'SHELL_B_TRIA'
         DO I=1,3
            WRITE(F06,99668) (SHELL_B0_TRIA(I,J),J=1,3), (SHELL_B_TRIA(I,J),J=1,3)
         ENDDO
         WRITE(F06,*)

      ENDIF mbend

      WRITE(F06,*)

      WRITE(F06,98799)

      WRITE(F06,*)

! **********************************************************************************************************************************
55566 FORMAT(1X,'------------------------------------------------------------------------------------------',/)

67549 FORMAT(6(1ES14.6))

99663 FORMAT(1X,A)

99664 FORMAT('  Transformation matrix ',2a)

99665 FORMAT(16x,'Material matrix ',a12,' before/after TF transformation',/,                                                       &
             19x,'before',41X,'after',/,'  ----------------------------------------      ----------------------------------------')

99667 FORMAT(3(1ES14.6), 4X,3(1ES14.6))

99668 FORMAT(3(1ES14.6), 4X,3(1ES14.6))

99669 FORMAT(7X,2(1ES14.6),18X,2(1ES14.6))

98720 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::START DEBUG(53) OUTPUT FROM SUBROUTINE ROT_AXES_MATL_TO_LOC:::::::::::::::::::',&
             ':::::::::::::::::',/)

98799 FORMAT(' :::::::::::::::::::::::::::::::::::::END DEBUG(53) OUTPUT FROM SUBROUTINE ROT_AXES_MATL_TO_LOC::::::::::::::::::::',&
             ':::::::::::::::::'                                                                                                ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE DEBUG_ROT_AXES_2

      END SUBROUTINE TPLT2

   END MODULE TRIANGULAR_SHELL_ELEMENTS
