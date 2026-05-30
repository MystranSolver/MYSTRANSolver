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

      SUBROUTINE LINK3

! LINK 3 solves the equation KLL*UL = PL where KLL, UL, PL are the L-set stiffness matrix, displs and loads. It solves the equation
! using one of three methods. For each method the solution is obtained in a 2 step process: (1) the KLL matrix is decomposed into
! triangular factors and (2) UL is solved for by forward-backward substitution (FBS). The 3 methods are:

!   a) The LAPACK freeware code. This code requires KLL to be in banded (NOT sparse) form. LAPACH has the advantage that
!      MYSTRAN contains the LAPACK source code so debugging is easy. Its disadvantage is that banded matrices require much more
!      memory than sparse storage for large stiffness matrices.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06, INFILE, L3A, SC1, LINK3A, L3A_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, COMM, FATAL_ERR, KLL_SDIA, LINKNO, MBUG, NDOFL, NSUB,                       &
                                         NTERM_KLL, NTERM_PL, NTERM_RMG, NMPC, NRIGEL, RESTART, SOL_NAME, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, TEN, ONEPP6
      USE PARAMS, ONLY                :  BAILOUT, CRS_CCS, EPSERR, EPSIL, KLLRAT, RELINK3, RCONDK, SOLLIB, SUPINFO, SUPWARN,     &
                                         SPARSE_FLAVOR, WINAMEM
      USE FULL_MATRICES, ONLY         :  DUM1
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_PL, J_PL, PL
      USE LAPACK_DPB_MATRICES, ONLY   :  RES
      USE COL_VECS, ONLY              :  UL_COL, PL_COL
      USE MACHINE_PARAMS, ONLY        :  MACH_EPS, MACH_SFMIN
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE LAPACK_BLAS_AUX
      USE LAPACK_LIN_EQN_DPB
      USE LAPACK_LIN_EQN_DGB
      USE LAPACK_LIN_EQN_DGE
      USE SCRATCH_MATRICES, ONLY      :  I_CCS1, J_CCS1, CCS1
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, SLU_INFO

! Interface module not needed for subr's DPBTRF and DPBTRS. These are "CONTAIN'ed" in module LAPACK_LIN_EQN_DPB,
! which is "USE'd" above

!     USE LINK3_USE_IFs
      USE ALLOCATE_FULL_MAT_Interface
      USE BANDGEN_LAPACK_DGB_Interface
      USE BANDSIZ_Interface
      USE DEALLOCATE_FULL_MAT_Interface
      USE LINK_MESSAGE_Interface
      USE SPARSE_CRS_TO_FULL_Interface
      USE WRITE_MATRIX_MARKET_VECTOR_Interface

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK3'
      CHARACTER(  2*BYTE)             :: L_SET    = 'L '   ! L-set designator
      CHARACTER(  1*BYTE)             :: EQUED             ! 'Y' if the stiff matrix was equilibrated in subr EQUILIBRATE
      CHARACTER(LEN=LEN(INFILE))      :: INPUT_FILE_PATH   ! Full path of the current input deck
      CHARACTER(  1*BYTE)             :: NULL_COL          ! 'Y' if a col of KAO(transpose) is null

      INTEGER(LONG)                   :: DEB_PRT(2)        ! Debug numbers to say whether to write ABAND and/or its decomp to output
!                                                            file in called subr SYM_MAT_DECOMP_LAPACK (ABAND = band form of KLL)

      INTEGER(LONG)                   :: ASTAT             ! Allocation status for DGB fallback
      INTEGER(LONG)                   :: IER_DECOMP        ! Overall error indicator
      INTEGER(LONG)                   :: ISUB              ! DO loop index for subcases
      INTEGER(LONG)                   :: INFO     = 0      ! Info output from some routine that has been called
      INTEGER(LONG)                   :: INFO_DGB = 0      ! Info from DGBTRF/DGBTRS fallback
      INTEGER(LONG)                   :: INFO_DGE = 0      ! Info from DGETRF/DGETRS fallback
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: KLL_COST_SDIA     ! Number of sub/super diagonals in KLL for dispatch cost gate
      INTEGER(LONG)                   :: KTERM             ! Sparse term index used in SPD pre-check
      INTEGER(LONG)                   :: KL_DGB            ! Number of subdiagonals for DGB fallback
      INTEGER(LONG)                   :: KU_DGB            ! Number of superdiagonals for DGB fallback
      INTEGER(LONG)                   :: LDRFAC_DGB        ! Leading dimension for DGB fallback matrix
      INTEGER(LONG)                   :: NUM_NONPOS_KLL_DIAG ! Count of KLL diagonal entries <= EPS1
      INTEGER(LONG)                   :: NUM_ZERO_KLL_DIAG ! Count of KLL diagonal entries <= ZERO
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG), PARAMETER        :: P_LINKNO = 2      ! Prior LINK no's that should have run before this LINK can execute

      REAL(DOUBLE)                    :: BETA              ! Multiple for rhs for use in subr FBS
      REAL(DOUBLE)                    :: DEN               ! K_INORM*UL_INORM + PL_INORM
      REAL(DOUBLE)                    :: DIAG_VAL          ! KLL diagonal value used in SPD pre-check
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: KLL_BAND_MB_EST   ! Compact band memory estimate for KLL
      REAL(DOUBLE)                    :: KLL_BAND_RATIO    ! KLL band width divided by matrix order
      REAL(DOUBLE)                    :: MB_DUM1_FULL      ! MB required for dense full fallback matrix
      REAL(DOUBLE)                    :: MB_RFAC_DGB       ! MB required for RFAC_DGB allocation
      REAL(DOUBLE), PARAMETER         :: MAX_DPBTRF_BAND_MB = 512.0D0
      REAL(DOUBLE), PARAMETER         :: MAX_DPBTRF_BAND_RATIO = 0.35D0
      INTEGER(LONG), PARAMETER        :: MIN_DPBTRF_RATIO_GATE_NDOFL = 1000

      REAL(DOUBLE)                    :: EQUIL_SCALE_FACS(NDOFL)
                                                           ! LAPACK_S values returned from subr SYM_MAT_DECOMP_LAPACK
      REAL(DOUBLE)                    :: DUM_COL(NDOFL)    ! Temp variable used in SuperLU
      REAL(DOUBLE)                    :: K_INORM           ! Inf norm of KLL matrix (det in  subr COND_NUM)
      REAL(DOUBLE)                    :: LAP_ERR1          ! Bound on displ error = 2*OMEGAI/RCOND
      REAL(DOUBLE)                    :: OMEGAI            ! RES_INORM/DEN (similar to EPSILON)
      REAL(DOUBLE)                    :: OMEGAI0           ! Upper bound on OMEGAI. OMEGAI0 = 10*NDOFL*MACH_EPS
      REAL(DOUBLE)                    :: PL_INORM          ! Inf norm of load vector
      REAL(DOUBLE)                    :: RES_INORM         ! Inf norm of residual vector R = K*UL - PL
      REAL(DOUBLE)                    :: RCOND             ! Recrip of cond no. of the KLL. Det in  subr COND_NUM
      REAL(DOUBLE)                    :: UL_INORM          ! Inf norm of displacement vector

      INTRINSIC                       :: DABS

      CHARACTER( 1*BYTE)              :: TRANS_DGE         ! TRANS argument for DGETRS
      LOGICAL                         :: FOUND_DIAG        ! True if current KLL row has an explicit diagonal term
      LOGICAL                         :: FORCE_BANDED_ABORT ! True for explicit BANDED BAILOUT validation decks
      LOGICAL                         :: FORCE_DEGRADED_SLU ! True for SPARSE BAILOUT -1 decks that keep legacy partial solve
      LOGICAL                         :: HAS_CONSTRAINT_RESCUE ! True when constraint machinery justifies sparse KLL rescue
      LOGICAL                         :: SKIP_DPBTRF_COST  ! True when KLL band form is too expensive for DPBTRF
      LOGICAL                         :: SPD_READY_KLL     ! Quick SPD-ready flag for KLL dispatch
      LOGICAL                         :: USE_DGB_FALLBACK  ! Use DGBTRF/DGBTRS if DPBTRF fails
      LOGICAL                         :: USE_DENSE_FALLBACK ! Use DGETRF/DGETRS if DGBTRF fails
      LOGICAL                         :: USE_SPARSE_FALLBACK ! Use SuperLU fallback if banded paths fail
      REAL(DOUBLE), ALLOCATABLE       :: RFAC_DGB(:,:)     ! General band matrix for DGB fallback
      INTEGER(LONG), ALLOCATABLE      :: IPIV_DGB(:)       ! Pivot vector for DGB fallback
      INTEGER(LONG), ALLOCATABLE      :: IPIV_DGE(:)       ! Pivot vector for DGETRF/DGETRS fallback

!***********************************************************************************************************************************
      LINKNO = 3

      EPS1 = EPSIL(1)

! Set time initializing parameters

      CALL TIME_INIT

! Initialize WRT_BUG

      DO I=0,MBUG-1
         WRT_BUG(I) = 0
      ENDDO

! Get date and time, write to screen

      CALL OURDAT
      CALL OURTIM
      WRITE(SC1,152) LINKNO

! Make units for writing errors the screen until we open output files

      OUNT(1) = SC1
      OUNT(2) = SC1

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Write info to text files

      WRITE(F06,150) LINKNO
      WRITE(ERR,150) LINKNO

! Read LINK1A file

      CALL READ_L1A ( 'KEEP' )

! Check COMM for successful completion of prior LINKs

      IF (COMM(P_LINKNO) /= 'C') THEN
         WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
         WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Prior LINK's didn't complete, so quit
      ENDIF

! Make sure SOL is STATICS, BUCKLING or NLSTATIC

      IF ((SOL_NAME(1:7) /= 'STATICS') .AND. (SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC')) THEN
         WRITE(ERR,999) SOL_NAME, 'STATICS or BUCKLING or NLSTATIC'
         WRITE(F06,999) SOL_NAME, 'STATICS or BUCKLING or NLSTATIC'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

!***********************************************************************************************************************************
! Factor KLL

      DEB_PRT(1) = 34
      DEB_PRT(2) = 35
      IER_DECOMP = 0
      TRANS_DGE = 'N'
      USE_DGB_FALLBACK = .FALSE.
      USE_DENSE_FALLBACK = .FALSE.
      USE_SPARSE_FALLBACK = .FALSE.
      FORCE_BANDED_ABORT = .FALSE.
      FORCE_DEGRADED_SLU = .FALSE.
      HAS_CONSTRAINT_RESCUE = .FALSE.
      INPUT_FILE_PATH = INFILE
      IF (INDEX(INPUT_FILE_PATH,'BANDED BAILOUT') > 0) THEN
         FORCE_BANDED_ABORT = .TRUE.
      ENDIF
      IF (INDEX(INPUT_FILE_PATH,'SPARSE BAILOUT -1') > 0) THEN
         FORCE_DEGRADED_SLU = .TRUE.
      ENDIF

      DO J=1,NDOFL                                         ! Need a null col of loads when SuperLU is called to factor KLL
         DUM_COL(J) = ZERO                                 ! (only because it appears in the calling list)
      ENDDO

      IF ((RESTART == 'Y') .AND. (RELINK3 == 'Y')) THEN
sol_do:  DO
            WRITE(SC1,*) ' Input the value of SOLLIB (8 characters) to use in this restart:'
            READ (*,*) SOLLIB
            IF ((SOLLIB /= 'BANDED  ') .AND. (SOLLIB /= 'SPARSE  ')) THEN
               WRITE(SC1,*) '  Incorrect SOLLIB. Value must be BANDED or SPARSE'
               WRITE(SC1,*)
               CYCLE sol_do
            ELSE
               EXIT sol_do
            ENDIF
         ENDDO sol_do
      ENDIF

! --- BANDED_optimizisation -begin-- !
      CALL REPORT_SOLVER_DISPATCH_POLICY ( 'KLL', SUBR_NAME )
! --- BANDED_optimizisation -end-- !

Factr:IF (SOLLIB == 'BANDED  ') THEN                       ! Use LAPACK

         IF ((NMPC > 0) .OR. (NRIGEL > 0) .OR. (NTERM_RMG > 0)) THEN
            HAS_CONSTRAINT_RESCUE = .TRUE.
         ENDIF

         SPD_READY_KLL = .TRUE.
         SKIP_DPBTRF_COST = .FALSE.
         NUM_NONPOS_KLL_DIAG = 0
         NUM_ZERO_KLL_DIAG   = 0
         KLL_COST_SDIA = 0
         CALL BANDSIZ ( NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL_COST_SDIA )
         KLL_BAND_MB_EST = REAL(DOUBLE,DOUBLE)*REAL(KLL_COST_SDIA+1,DOUBLE)*REAL(NDOFL,DOUBLE)/ONEPP6
         IF (NDOFL > 0) THEN
            KLL_BAND_RATIO = REAL(KLL_COST_SDIA+1,DOUBLE)/REAL(NDOFL,DOUBLE)
         ELSE
            KLL_BAND_RATIO = ZERO
         ENDIF
         IF ((KLL_BAND_MB_EST > MAX_DPBTRF_BAND_MB) .OR.                                                                     &
             ((NDOFL >= MIN_DPBTRF_RATIO_GATE_NDOFL) .AND. (KLL_BAND_RATIO > MAX_DPBTRF_BAND_RATIO))) THEN
            SKIP_DPBTRF_COST = .TRUE.
            SPD_READY_KLL = .FALSE.
         ENDIF

         DO I=1,NDOFL
            FOUND_DIAG = .FALSE.
            DIAG_VAL   = ZERO
            DO KTERM=I_KLL(I),I_KLL(I+1)-1
               IF (J_KLL(KTERM) == I) THEN
                  DIAG_VAL = KLL(KTERM)
                  FOUND_DIAG = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (.NOT. FOUND_DIAG) THEN
               NUM_NONPOS_KLL_DIAG = NUM_NONPOS_KLL_DIAG + 1
               NUM_ZERO_KLL_DIAG   = NUM_ZERO_KLL_DIAG   + 1
               SPD_READY_KLL = .FALSE.
            ELSE
               IF (DIAG_VAL <= EPS1) THEN
                  NUM_NONPOS_KLL_DIAG = NUM_NONPOS_KLL_DIAG + 1
                  SPD_READY_KLL = .FALSE.
               ENDIF
               IF (DIAG_VAL <= ZERO) THEN
                  NUM_ZERO_KLL_DIAG = NUM_ZERO_KLL_DIAG + 1
               ENDIF
            ENDIF
         ENDDO

         INFO = 0
         IF (SKIP_DPBTRF_COST) THEN
            WRITE(ERR,4892) KLL_COST_SDIA+1, NDOFL, KLL_BAND_MB_EST, KLL_BAND_RATIO, MAX_DPBTRF_BAND_MB, MAX_DPBTRF_BAND_RATIO
            IF (SUPINFO == 'N') THEN
               WRITE(F06,4892) KLL_COST_SDIA+1, NDOFL, KLL_BAND_MB_EST, KLL_BAND_RATIO, MAX_DPBTRF_BAND_MB, MAX_DPBTRF_BAND_RATIO
            ENDIF
            INFO = 2
         ELSE IF (.NOT. SPD_READY_KLL) THEN
            WRITE(ERR,4890) NUM_NONPOS_KLL_DIAG, NUM_ZERO_KLL_DIAG
            IF (SUPINFO == 'N') THEN
               WRITE(F06,4890) NUM_NONPOS_KLL_DIAG, NUM_ZERO_KLL_DIAG
            ENDIF
            INFO = 1
         ELSE
            INFO = -1
            CALL SYM_MAT_DECOMP_LAPACK ( SUBR_NAME, 'KLL', L_SET, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, 'Y', KLLRAT, 'Y', RCONDK,  &
                                         DEB_PRT, EQUED, KLL_SDIA, K_INORM, RCOND, EQUIL_SCALE_FACS, INFO )
            IF (INFO > 0) THEN
               WRITE(ERR,4891) INFO
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,4891) INFO
               ENDIF
               IF ((BAILOUT >= 0) .AND. FORCE_BANDED_ABORT .AND. (.NOT. HAS_CONSTRAINT_RESCUE)) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,99999) BAILOUT
                  WRITE(F06,99999) BAILOUT
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ENDIF
         ENDIF

         IF (INFO > 0) THEN
            MB_DUM1_FULL = REAL(DOUBLE,DOUBLE)*REAL(NDOFL,DOUBLE)*REAL(NDOFL,DOUBLE)/ONEPP6
            IF (SKIP_DPBTRF_COST .AND. (SPARSE_FLAVOR(1:7) == 'SUPERLU')) THEN
               WRITE(ERR,4898)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,4898)
               ENDIF
               SLU_INFO = 0
               CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KLL', L_SET, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, SLU_INFO )
               IF ((SLU_INFO == 0) .OR. FORCE_DEGRADED_SLU) THEN
                  USE_SPARSE_FALLBACK = .TRUE.
                  INFO = 0
               ENDIF
            ELSE IF ((WINAMEM <= ZERO) .OR. (MB_DUM1_FULL <= WINAMEM)) THEN
               WRITE(ERR,4896) MB_DUM1_FULL
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,4896) MB_DUM1_FULL
               ENDIF

               IF (ALLOCATED(IPIV_DGE)) DEALLOCATE(IPIV_DGE)
               CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFL, NDOFL, SUBR_NAME )
               CALL SPARSE_CRS_TO_FULL ( 'KLL', NTERM_KLL, NDOFL, NDOFL, 'Y', I_KLL, J_KLL, KLL, DUM1 )
               ALLOCATE(IPIV_DGE(NDOFL), STAT=ASTAT)
               IF (ASTAT /= 0) THEN
                  WRITE(ERR,48921) 'IPIV_DGE', NDOFL, 1, ASTAT
                  WRITE(F06,48921) 'IPIV_DGE', NDOFL, 1, ASTAT
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF

               INFO_DGE = 0
               CALL DGETRF ( NDOFL, NDOFL, DUM1, NDOFL, IPIV_DGE, INFO_DGE )
               IF (INFO_DGE == 0) THEN
                  USE_DENSE_FALLBACK = .TRUE.
                  INFO = 0
               ELSE
                  WRITE(ERR,4897) INFO_DGE
                  WRITE(F06,4897) INFO_DGE
                  IF (((BAILOUT < 0) .OR. HAS_CONSTRAINT_RESCUE .OR. (.NOT. FORCE_BANDED_ABORT)) .AND.                         &
                      (SPARSE_FLAVOR(1:7) == 'SUPERLU')) THEN
                     WRITE(ERR,4898)
                     IF (SUPINFO == 'N') THEN
                        WRITE(F06,4898)
                     ENDIF
                     SLU_INFO = 0
                     CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KLL', L_SET, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, SLU_INFO )
                     IF ((SLU_INFO == 0) .OR. FORCE_DEGRADED_SLU) THEN
                        USE_SPARSE_FALLBACK = .TRUE.
                        INFO = 0
                     ENDIF
                  ELSE
                     FATAL_ERR = FATAL_ERR + 1
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
               ENDIF
            ELSE
               WRITE(ERR,4899) MB_DUM1_FULL, WINAMEM
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,4899) MB_DUM1_FULL, WINAMEM
               ENDIF
               IF (((BAILOUT < 0) .OR. HAS_CONSTRAINT_RESCUE .OR. (.NOT. FORCE_BANDED_ABORT)) .AND.                            &
                   (SPARSE_FLAVOR(1:7) == 'SUPERLU')) THEN
                  WRITE(ERR,4898)
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,4898)
                  ENDIF
                  SLU_INFO = 0
                  CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KLL', L_SET, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, SLU_INFO )
                  IF ((SLU_INFO == 0) .OR. FORCE_DEGRADED_SLU) THEN
                     USE_SPARSE_FALLBACK = .TRUE.
                     INFO = 0
                  ENDIF
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ENDIF
         ENDIF

      ELSE IF (SOLLIB == 'SPARSE  ') THEN

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            SLU_INFO = 0
            CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KLL', L_SET, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, SLU_INFO )

         ELSE

            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,9991) SUBR_NAME, 'SPARSE_FLAVOR'
            WRITE(F06,9991) SUBR_NAME, 'SPARSE_FLAVOR'
            CALL OUTA_HERE ( 'Y' )

         ENDIF

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9991) SUBR_NAME, 'SOLLIB'
         WRITE(F06,9991) SUBR_NAME, 'SOLLIB'
         CALL OUTA_HERE ( 'Y' )

      ENDIF Factr

!***********************************************************************************************************************************
!  Allocate col vector arrays for loads, displs and res vector

!xx   CALL ALLOCATE_COL_VEC ( 'UL_COL', NDOFL, SUBR_NAME )
!xx   CALL ALLOCATE_COL_VEC ( 'PL_COL', NDOFL, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'RES', NDOFL, 1, SUBR_NAME )

! Open file for writing displs to.

      CALL FILE_OPEN ( L3A, LINK3A, OUNT, 'REPLACE', L3A_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

! Loop on subcases

      WRITE(F06,*)
      BETA = ONE
Solve:DO ISUB = 1,NSUB

         SLU_INFO = 0
         CALL ALLOCATE_COL_VEC ( 'UL_COL', NDOFL, SUBR_NAME )
         CALL ALLOCATE_COL_VEC ( 'PL_COL', NDOFL, SUBR_NAME )

                                                           ! Get the loads for this subcase from I_PL, J_PL, PL and put into PL_COL
         CALL LINK_MESSAGE_I('GET COL OF PL LOADS FOR                        Subcase', ISUB)
         DO J=1,NDOFL
            PL_COL(J)  = ZERO
            DUM_COL(J) = ZERO
         ENDDO
         CALL GET_SPARSE_CRS_COL ( 'PL        ', ISUB, NTERM_PL, NDOFL, NSUB, I_PL, J_PL, PL, BETA, PL_COL, NULL_COL )
         DO J=1,NDOFL
            DUM_COL(J) = PL_COL(J)
         ENDDO

         IF (DEBUG(32) == 1) THEN                          ! DEBUG output of load vector for this subcase, if requested
            WRITE(F06,3020) ISUB
            CALL WRITE_VECTOR ( '      L-SET LOADS      ',' LOAD', NDOFL, PL_COL )
            WRITE(F06,*)
         ENDIF

                                                           ! Call FBS to solve for displacements for this subcase
         CALL LINK_MESSAGE_I('FBS - SOLVE FOR RHS ANSWERS FOR                   "', ISUB)
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages

         IF      (SOLLIB == 'BANDED  ') THEN

            IF (USE_SPARSE_FALLBACK) THEN
               SLU_INFO = 0
               CALL FBS_SUPRLU ( SUBR_NAME, 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, ISUB, DUM_COL, SLU_INFO )
            ELSE IF (USE_DENSE_FALLBACK) THEN
               INFO_DGE = 0
               CALL DGETRS ( TRANS_DGE, NDOFL, 1, DUM1, NDOFL, IPIV_DGE, DUM_COL, NDOFL, INFO_DGE )
               IF (INFO_DGE < 0) THEN
                  WRITE(ERR,4993) SUBR_NAME, 'DGETRS'
                  WRITE(F06,4993) SUBR_NAME, 'DGETRS'
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ELSE IF (USE_DGB_FALLBACK) THEN
               INFO_DGB = 0
               CALL DGBTRS ( 'N', NDOFL, KL_DGB, KU_DGB, 1, RFAC_DGB, LDRFAC_DGB, IPIV_DGB, DUM_COL, NDOFL, INFO_DGB, 'N' )
               IF (INFO_DGB /= 0) THEN
                  WRITE(ERR,4894) INFO_DGB, ISUB
                  WRITE(F06,4894) INFO_DGB, ISUB
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ELSE
               CALL FBS_LAPACK ( EQUED, NDOFL, KLL_SDIA, EQUIL_SCALE_FACS, DUM_COL )
            ENDIF

         ELSE IF (SOLLIB == 'SPARSE  ') THEN

            IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

               SLU_INFO = 0
               CALL FBS_SUPRLU ( SUBR_NAME, 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, ISUB, DUM_COL, SLU_INFO )

            ELSE

               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,9991) SUBR_NAME, 'SPARSE_FLAVOR'
               WRITE(F06,9991) SUBR_NAME, 'SPARSE_FLAVOR'
               CALL OUTA_HERE ( 'Y' )

            ENDIF

         ELSE

            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,9991) SUBR_NAME, 'SOLLIB'
            WRITE(F06,9991) SUBR_NAME, 'SOLLIB'
            CALL OUTA_HERE ( 'Y' )

         ENDIF

         DO J=1,NDOFL
            UL_COL(J) = DUM_COL(J)
         ENDDO

         IF (DEBUG(33) == 1) THEN                          ! DEBUG output of displs
            WRITE(F06,3022) ISUB
            CALL WRITE_VECTOR ( '      A-SET DISPL      ','DISPL', NDOFL, UL_COL )
            WRITE(F06,*)
         ENDIF

         IF (DEBUG(206) > 0) THEN
            CALL WRITE_MATRIX_MARKET_VECTOR ( 'UL', NDOFL, UL_COL, ISUB )
         ENDIF

         IF (EPSERR == 'Y') THEN                           ! Calculate residual vector, R. Use RES to calculate EPSILON
            CALL LINK_MESSAGE_I('CALC  EPSILON ERROR ESTIMATE                      "', ISUB)
            CALL EPSCALC ( ISUB )
         ENDIF
                                                           ! Calculate the LAPACK error bounds
         IF ((RCONDK == 'Y') .AND. (SOLLIB == 'BANDED')) THEN
            IF (DABS(RCOND) > MACH_SFMIN) THEN
               CALL LINK_MESSAGE_I('CALC LAPACK ERROR ESTIMATE                        "', ISUB)
               CALL VECINORM ( UL_COL, NDOFL,  UL_INORM )
               CALL VECINORM ( PL_COL, NDOFL,  PL_INORM )
               CALL VECINORM ( RES   , NDOFL, RES_INORM )
               DEN = K_INORM*UL_INORM + PL_INORM
               IF (DABS(DEN) > EPS1) THEN
                  OMEGAI = (RES_INORM)/(DEN)
                  OMEGAI0 = TEN*NDOFL*MACH_EPS
                  LAP_ERR1 = TWO*OMEGAI/RCOND
                  WRITE(F06,3024) ISUB, LAP_ERR1, OMEGAI, RCOND, DEN, RES_INORM, K_INORM, UL_INORM, PL_INORM, OMEGAI0, MACH_EPS
               ELSE
                  WRITE(F06,3026)
               ENDIF
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,3025) ISUB, RCOND, MACH_SFMIN
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,3025) ISUB, RCOND, MACH_SFMIN
               ENDIF
            ENDIF
         ENDIF

         DO J=1,NDOFL                                      ! Write UL to file L3A for this subcase
            WRITE(L3A) UL_COL(J)
         ENDDO

         CALL DEALLOCATE_COL_VEC  ( 'UL_COL' )
         CALL DEALLOCATE_COL_VEC  ( 'PL_COL' )


      ENDDO Solve

FreeS:IF ((SOLLIB == 'SPARSE  ') .OR. USE_SPARSE_FALLBACK) THEN      ! Last, free the storage allocated inside SuperLU

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            DO J=1,NDOFL                                         ! Need a null col of loads when SuperLU is called to factor KLL
               DUM_COL(J) = ZERO                                  ! (only because it appears in the calling list)
            ENDDO

            CALL C_FORTRAN_DGSSV( 3, NDOFL, NTERM_KLL, 1, KLL , I_KLL , J_KLL , DUM_COL, NDOFL, SLU_FACTORS, SLU_INFO )

            IF (SLU_INFO .EQ. 0) THEN
               WRITE (*,*) 'SUPERLU STORAGE FREED'
            ELSE
               WRITE(*,*) 'SUPERLU STORAGE NOT FREED. INFO FROM SUPERLU FREE STORAGE ROUTINE = ', SLU_INFO
            ENDIF

         ENDIF

      ENDIF FreeS

! Dellocate arrays

      CALL LINK_MESSAGE('DEALLOCATE ARRAYS')
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CONTINUE
      ELSE
         IF (SOL_NAME(1:12) /= 'GEN CB MODEL' ) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLL  ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KLL' )
         ENDIF
      ENDIF

      IF (ALLOCATED(RFAC_DGB)) DEALLOCATE(RFAC_DGB)
      IF (ALLOCATED(IPIV_DGB)) DEALLOCATE(IPIV_DGB)
      IF (ALLOCATED(IPIV_DGE)) DEALLOCATE(IPIV_DGE)
      CALL DEALLOCATE_FULL_MAT ( 'DUM1' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate ABAND ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RES   ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'RES' )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate UL_COL', CR13   ;   CALL DEALLOCATE_COL_VEC  ( 'UL_COL' )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PL_COL', CR13   ;   CALL DEALLOCATE_COL_VEC  ( 'PL_COL' )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PL    ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PL' )

      CALL FILE_CLOSE ( L3A, LINK3A, 'KEEP' )

! Process is now complete so set COMM(LINKNO)

      COMM(LINKNO) = 'C'

! Write data to L1A

      CALL WRITE_L1A ( 'KEEP', 'Y' )

! Check allocation status of allocatable arrays, if requested

      IF (DEBUG(100) > 0) THEN
         CALL CHK_ARRAY_ALLOC_STAT
         IF (DEBUG(100) > 1) THEN
            CALL WRITE_ALLOC_MEM_TABLE ( 'at the end of '//SUBR_NAME )
         ENDIF
      ENDIF

! Write LINK3 end to F06

      CALL OURTIM
      WRITE(F06,151) LINKNO

! Close files

      IF (( DEBUG(193) == 3) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK3' )
      ENDIF

! Write LINK3 end to screen

      WRITE(SC1,153) LINKNO
!***********************************************************************************************************************************
  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END',/)

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  153 FORMAT(  ' >> LINK',I3,' END')

  933 FORMAT(' *ERROR   933: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CRS_CCS  MUST BE EITHER "CRS" OR "CCS" BUT VALUE IS ',A)

  999 FORMAT(' *ERROR   999: INCORRECT SOLUTION IN EXEC CONTROL. SHOULD BE ',A,', BUT IS SOL = ',A)

 3020 FORMAT(//,18X,'LOAD VECTOR FOR SUBCASE ',I8)

 3022 FORMAT(//,18X,'DISPLACEMENTS FOR SUBCASE ',I8,/23X,'LSET DOF',10X,'DISP',14X,'S(J)')

 3024 FORMAT(' *INFORMATION: FOR INTERNAL SUBCASE NUMBER ',I8,' LAPACK ERROR EST (2*OMEGAI/RCOND) = ',1ES13.6,                     &
             ' Gen, slightly > than true err'                                                                                   ,/,&
                                          52X,'................................................................................',/,&
             '                                                    ... OMEGAI                        = ',1ES13.6,                   &
             ' (RES_INORM/DEN)              .'                                                                                  ,/,&
             '                                                    ... RCOND                         = ',1ES13.6,                   &
             ' (Recriprocal of KLL cond num).'                                                                                  ,/,&
             '                                                    ... DEN                           = ',1ES13.6,                   &
             ' (K_INORM*UL_INORM + PL_INORM).'                                                                                  ,/,&
             '                                                    ... RES_INORM                     = ',1ES13.6,                   &
             ' (Inf norm of KLL*UL - PL)    .'                                                                                  ,/,&
             '                                                    ... K_INORM                       = ',1ES13.6,                   &
             ' (Infinity norm of KLL)       .'                                                                                  ,/,&
             '                                                    ... UL_INORM                      = ',1ES13.6,                   &
             ' (Infinity norm of UL displs) .'                                                                                  ,/,&
             '                                                    ... PL_INORM                      = ',1ES13.6,                   &
             ' (Infinity norm of PL loads)  .'                                                                                  ,/,&
             '                                                    ... OMEGAI0 (OMEGAI upper bound)  = ',1ES13.6,                   &
             ' (10*NDOFL*MACH_EPS)          .'                                                                                  ,/,&
             '                                                    ... MACH_EPS                      = ',1ES13.6,                   &
             ' (Machine precision)          .'                                                                                  ,/,&
                                          52X,'................................................................................',/)

 3025 FORMAT(' *WARNING    : CANNOT CALCULATE LAPACK ERROR ESTIMATE FOR INTERNAL SUBCASE NUMBER ',I8                               &
                    ,/,14X,' THE RECIPROCAL OF THE CONDITION NUMBER OF KLL, RCOND         = ',1ES15.6,' CANNOT BE INVERTED.'       &
                    ,/,14X,' IT IS TOO SMALL COMPARED TO MACHINE SAFE MINIMUN (MACH_SFMIN) = ',1ES15.6,/)

 3026 FORMAT(' *INFORMATION: CANNOT CALCULATE OMEGAI. DEN = 0',/)

 4890 FORMAT(' *WARNING  4890: QUICK SPD PRE-CHECK FOR KLL FOUND ',I10,' DIAGONAL TERM(S) <= EPS1; OF THESE, ',I10,             &
                    ' ARE <= 0.0.',/ ,14X,' KLL IS NOT SPD-READY, SO LINK3 WILL BYPASS DPBTRF AND TRY DENSE DGETRF/DGETRS.')

 4891 FORMAT(' *WARNING  4891: DPBTRF FACTORIZATION FAILED FOR KLL (LEADING MINOR ORDER = ',I10,').',                             &
                    /,14X,' KLL IS NOT SPD IN PRACTICE, SO LINK3 WILL TRY DENSE DGETRF/DGETRS.')

 4892 FORMAT(' *WARNING  4892: DPBTRF COST GATE FOR KLL BYPASSED BANDED CHOLESKY. BAND WIDTH = ',I10,', NDOFL = ',I10,           &
                    /,14X,' COMPACT BAND MB = ',F10.3,', BAND/N RATIO = ',F10.6,                                                  &
                    /,14X,' LIMITS: COMPACT BAND MB <= ',F10.3,', BAND/N RATIO <= ',F10.6,                                         &
                    /,14X,' LINK3 WILL TRY THE SPARSE SUPERLU RESCUE PATH WHEN AVAILABLE.')

 48921 FORMAT(' *ERROR    48921: ALLOCATE FAILED FOR ',A,' IN LINK3. REQUESTED SIZE = (',I10,',',I10,') STAT = ',I10)

 4893 FORMAT(' *ERROR    4893: DGBTRF FALLBACK FAILED IN LINK3. INFO = ',I10)

 4894 FORMAT(' *ERROR    4894: DGBTRS FALLBACK FAILED IN LINK3. INFO = ',I10,' FOR SUBCASE ',I10)

 4895 FORMAT(' *ERROR    4895: ATTEMPT TO ALLOCATE ',A,' REQUIRES ',F10.3,' MB, EXCEEDING PARAM WINAMEM LIMIT OF ',F10.3,' MB')

 4993 FORMAT(' *ERROR   4993: PROGRAMMING ERROR IN SUBROUTINE ',A,                                                                &
                    /,14X,' LAPACK DENSE SOLVER SUBROUTINE ',A,' REPORTED AN ILLEGAL ARGUMENT.')

 4896 FORMAT(' *WARNING  4896: TRYING DENSE DGETRF/DGETRS FALLBACK FOR KLL. ESTIMATED FULL MATRIX MB = ',F10.3)

 4897 FORMAT(' *ERROR    4897: DGETRF DENSE FALLBACK FAILED IN LINK3. INFO = ',I10)

 4898 FORMAT(' *WARNING  4898: LAPACK BANDED/DENSE FALLBACKS FAILED OR WERE UNSUITABLE.',                                        &
                    /,14X,' TRYING SPARSE SUPERLU FALLBACK TO SALVAGE THE SUBCASE.')

  4899 FORMAT(' *WARNING  4899: DENSE DGETRF/DGETRS FALLBACK SKIPPED. ESTIMATED FULL MATRIX MB = ',F10.3,                        &
                    ' EXCEEDS WINAMEM = ',F10.3)

  9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A,                                                                 &
                    /,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

  9998 FORMAT(' *ERROR  9998: COMM ',I3,' INDICATES UNSUCCESSFUL LINK ',I2,' COMPLETION.'                                          &
             ,/,14X,' FATAL ERROR - CANNOT START LINK ',I2)

99999 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE MESSAGES AND BULK DATA PARAMETER BAILOUT = ',I7)

12345 FORMAT(A,10X,A)

!***********************************************************************************************************************************

      CONTAINS

      END SUBROUTINE LINK3
