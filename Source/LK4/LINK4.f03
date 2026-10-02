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

   MODULE LINK4_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: LINK4

   CONTAINS

   SUBROUTINE LINK4

! Calculates system eigenvalues, eigenvectors. There are 4 eigenvalue extraction methods in MYSTRAN, none of which seem suited to
! very large eigenvalue problems for one reason or another:

!   (1) LANCZOS method:

!          calculates some eigenvalues and some eigenvectors of KLL, MLL (or KLL, KLLD for BUCKLING). This is the most widely used
!          eigenvalue extraction method for large problems. The LANCZOS method in MYSTRAN uses 1 of 2 algorithms:
!             (a) LAPACK: requires KLL, MLL (or KLL, KLLD for buckling) to be in band storage - sparse storage can NOT be used

!          Tis Lanczos algorithms is not practical for very large eigenvalue problems since LAPACK/ARPACK will require
!          large amounts of memory to store the banded KLL, MLL matrices.

!   (2) GIV (Givens) method:
!          calculates all eigenvalues and some eigenvactors of KLL, MLL (or KLL, KLLD for buckling). This method is only practical
!          for relatively small problems. It requires MLL (or KLLD for buckling) to be a positive definite matrix. The algorithm
!          performs a Cholesky decomp of matrix MLL (or KLLD for buckling) which can be time consuming for large problems

!   (3) MGIV (modified Givens) method:
!          calculates all eigenvalues and some eigenvactors of KLL, MLL (or KLL, KLLD for buckling). This method is only practical
!          for relatively small problems. It requires KLL to be a positive definite matrix. The algorithm performs a Cholesky
!          decomp of matrix KLL which can be time consuming for large problems

!   (4) INV (Inverse Power) method:
!          calculates only the lowest eigenvalue and its eigenvector of KLL, MLL (or KLL, KLLD for buckling).


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, ERR, ERRSTAT, F06, L1M, L3A, SC1
      USE IOUNT1, ONLY                :  LINK1M,  LINK2I,  LINK3A, L1M_MSG, L3A_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, COMM, FATAL_ERR, LINKNO, MBUG, NDOFG, NDOFL, NSUB,                           &
                                         NTERM_KLL, NTERM_KLLD, NTERM_KLLDn,                                                       &
                                         NTERM_MLL, NTERM_MLLn,                                                                    &
                                         NVEC, NUM_EIGENS, NUM_KLLD_DIAG_ZEROS, NUM_MLL_DIAG_ZEROS, SOL_NAME, WARN_ERR,            &
                                         MODE_SUBCASE, NUM_MODES_SUBS, NUM_BUCKLING_SUBS, TOTAL_MODES, NUM_SUBC_CARDS
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  EPSIL, SOLLIB, SPARSTOR, SUPINFO
      USE MODEL_STUF, ONLY            :  EIG_COMP, EIG_CRIT, EIG_FRQ1, EIG_FRQ2, EIG_GRID, EIG_METH, EIG_MSGLVL, EIG_LAP_MAT_TYPE, &
                                         EIG_MODE, EIG_N1, EIG_N2, EIG_NCVFACL, EIG_NORM, EIG_SID, EIG_SIGMA, EIG_VECS, MAXMIJ,    &
                                         MIJ_COL, MIJ_ROW, NUM_FAIL_CRIT, EIG_PARAMS, IS_MODES_SUBCASE, IS_BUCKLING_SUBCASE,       &
                                         NUM_EIGENS_SUB, CC_EIGR_SID, SCNUM

      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_KLLD, J_KLLD, KLLD, I_KLLDn, J_KLLDn, KLLDn,                         &
                                         I_MLL, J_MLL, MLL, I_MLLn, J_MLLn, MLLn
      USE OUTPUT4_MATRICES, ONLY      :  NUM_OU4_REQUESTS
      USE EIGEN_MATRICES_1, ONLY      :  GEN_MASS, MODE_NUM, EIGEN_VAL, EIGEN_VEC
      USE LAPACK_DPB_MATRICES, ONLY   :  ABAND, BBAND
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE GIVENS, ONLY                :  EIG_GIV_MGIV
      USE INVERSE_POWER, ONLY         :  EIG_INV_PWR
      USE LANCZOS, ONLY               :  EIG_LANCZOS_ARPACK, EIG_LANCZOS_ARPACK_ADAPTIVE
      USE EIGEN_SUPPORT, ONLY         :  CALC_GEN_MASS, RENORM_ON_MASS, EIG_SUMMARY

      USE DATE_TIME_UTILS, ONLY       :  OURDAT, OURTIM, TIME_INIT
      USE TEMP_FILE_READERS, ONLY     :  READ_L1A, READ_L1M
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE, WRITE_L1A
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_MAT_DIAG_ZEROS
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_SYM_TO_CRS_NONSYM
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE EIGEN_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_EIGEN1_MAT, DEALLOCATE_EIGEN1_MAT
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_L1M
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_INQUIRE, FILE_OPEN
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_VECTOR
      USE OUTPUT4_FILE_IO, ONLY       :  OUTPUT4_PROC
      USE LAPACK_MATRIX_LIFECYCLE, ONLY:  DEALLOCATE_LAPACK_MAT
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_ARRAY_ALLOC_STAT, LINK_MESSAGE, WRITE_ALLOC_MEM_TABLE
      USE RESTART_FILE_IO, ONLY       :  READ_L5A_UG_FOR_SUBCASE
      USE LINK2_MOD, ONLY             :  REBUILD_KLLD_FROM_KGGD
      USE COL_VEC_LIFECYCLE, ONLY     :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
      USE MODEL_STORAGE_DEALLOCATION, ONLY:  DEALLOCATE_MODEL_STUF

   USE MATRIX_FILE_IO, ONLY           :  WRITE_SPARSE_CRS
      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK4'

      INTEGER(LONG)                   :: I,J                 ! DO loop indices or counters.
      INTEGER(LONG)                   :: IERROR              ! Error count when reading records from a file.
      INTEGER(LONG)                   :: OUNT(2)             ! File units to write messages to. Input to subr UNFORMATTED_OPEN.
      INTEGER(LONG), PARAMETER        :: P_LINKNO = 2        ! Prior LINK no's that should have run before this LINK can execute.

      ! Multi-subcase modes-loop locals (per-subcase eigensolve dispatch for SOL 103 with multiple METHODs)
      INTEGER(LONG)                   :: CUR_ISUB            ! Subcase index currently being solved
      INTEGER(LONG)                   :: CANONICAL_ISUB      ! Subcase whose EIGR/EIGRL params were written to L1M
      INTEGER(LONG)                   :: N_MODES_ITER        ! Number of solver iterations (= max(1, NUM_MODES_SUBS))
      INTEGER(LONG)                   :: ITER, KCNT, KMODE, IDX, TOTAL_MODES_LOCAL
      INTEGER(LONG)                   :: NTERM_KLL_BAK       ! Backup of NTERM_KLL for multi-iter restore
      INTEGER(LONG), ALLOCATABLE      :: I_KLL_BAK(:)        ! Shadow of I_KLL across solver iterations (eigensolver deallocates KLL)
      INTEGER(LONG), ALLOCATABLE      :: J_KLL_BAK(:)        ! Shadow of J_KLL across solver iterations
      REAL(DOUBLE),  ALLOCATABLE      :: KLL_BAK(:)          ! Shadow of KLL across solver iterations

      ! BUCKLING multi-subcase per-iter preload swap state (SOL 105 with multiple buckling subcases / distinct STATSUBs)
      LOGICAL                         :: IS_BUCK_MULTI       ! True when SOL 105 with NUM_BUCKLING_SUBS > 1
      INTEGER(LONG)                   :: CURRENT_PRELOAD_ISUB! Static-subcase ISUB whose UG_COL is currently loaded / whose KLLD is built
      INTEGER(LONG)                   :: TARGET_PRELOAD      ! EIG_PARAMS(CUR_ISUB)%STATSUB_REF for the iter being solved
      INTEGER(LONG)                   :: IERR_RELOAD         ! Return status from READ_L5A_UG_FOR_SUBCASE

      REAL(DOUBLE)                    :: EPS1                ! Small number to compare variables against zero.
      REAL(DOUBLE)                    :: EIGEN_VEC_COL(NDOFL)! One eigenvector put into a 1-D array.
      LOGICAL                         :: WRITE_MLL           ! write the MLL matrix
! **********************************************************************************************************************************
      LINKNO = 4

      EPS1   = EPSIL(1)
      WRITE_MLL = (DEBUG(42) == 2)


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

      ! Make sure we have correct SOL
      IF ((SOL_NAME(1:5) /= 'MODES') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL') .AND. (SOL_NAME(1:8) /= 'BUCKLING')) THEN
         WRITE(ERR,999) 'MODES or BUCKLING or GEN CB MODEL', SOL_NAME
         WRITE(F06,999) 'MODES or BUCKLING or GEN CB MODEL', SOL_NAME
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      ! Read data from file LINK1M
      CALL READ_L1M ( IERROR )

      IF (DEBUG(184) > 0) THEN
         WRITE(F06,*   ) ' Data written to file L1M'
         WRITE(F06,9102) '   EIG_SID         ', EIG_SID
         WRITE(F06,9101) '   EIG_METH        ', EIG_METH
         WRITE(F06,9103) '   EIG_FRQ1        ', EIG_FRQ1
         WRITE(F06,9103) '   EIG_FRQ2        ', EIG_FRQ2
         WRITE(F06,9102) '   EIG_N1          ', EIG_N1
         WRITE(F06,9102) '   EIG_N2          ', EIG_N2
         WRITE(F06,9101) '   EIG_VECS        ', EIG_VECS
         WRITE(F06,9103) '   EIG_CRIT        ', EIG_CRIT
         WRITE(F06,9101) '   EIG_NORM        ', EIG_NORM
         WRITE(F06,9102) '   EIG_GRID        ', EIG_GRID
         WRITE(F06,9102) '   EIG_COMP        ', EIG_COMP
         WRITE(F06,9102) '   EIG_MODE        ', EIG_MODE
         WRITE(F06,9103) '   EIG_SIGMA       ', EIG_SIGMA
         WRITE(F06,9101) '   EIG_LAP_MAT_TYPE', EIG_LAP_MAT_TYPE
         WRITE(F06,9102) '   EIG_MSGLVL      ', EIG_MSGLVL
         WRITE(F06,9102) '   EIG_NCVFACL     ', EIG_NCVFACL
         WRITE(F06,9102) '   NUM_FAIL_CRIT   ', NUM_FAIL_CRIT
         WRITE(F06,9103) '   MAXMIJ          ', MAXMIJ
         WRITE(F06,9102) '   MIJ_ROW         ', MIJ_ROW
         WRITE(F06,9102) '   MIJ_COL         ', MIJ_COL
         WRITE(F06,*)
      ENDIF

      IF (IERROR > 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      ! NUM_MLL_DIAG_ZEROS will be used for a message written
      ! when the eigen summary is printed in subr EIG_SUMMARY
      ! (if more than this number of eigens are requested)
      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CONTINUE
      ELSE
         CALL SPARSE_MAT_DIAG_ZEROS ( 'MLL', NDOFL, NTERM_MLL, I_MLL, J_MLL, NUM_MLL_DIAG_ZEROS )
      ENDIF

! Generate nonsymmetric storage form for KLLD (if BUCKLING soln) or MLL. This is done since subr MATMULT_SFF, used when subr DSBAND
! is called herein, will run faster. MATMULT_SFF is called in each "Reverse commumication loop" in DSBAND for the LANCZOS method.

      IF      (SPARSTOR == 'SYM   ') THEN

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN

            CALL SPARSE_MAT_DIAG_ZEROS ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, NUM_KLLD_DIAG_ZEROS )
            NTERM_KLLDn = 2*NTERM_KLLD  - (NDOFL - NUM_KLLD_DIAG_ZEROS)

            CALL LINK_MESSAGE('ALLOCATE SPARSE KLLDn ARRAYS')
            CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )

            CALL LINK_MESSAGE('CONVERT SYM CRS KLLD TO NONSYM CRS KLLDn')
            CALL CRS_SYM_TO_CRS_NONSYM ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, KLLD, 'KLLDn', NTERM_KLLDn,                    &
                                         I_KLLDn, J_KLLDn, KLLDn, 'Y' )

         ELSE

            CALL SPARSE_MAT_DIAG_ZEROS ( 'MLL', NDOFL, NTERM_MLL, I_MLL, J_MLL, NUM_MLL_DIAG_ZEROS )
            NTERM_MLLn = 2*NTERM_MLL  - (NDOFL - NUM_MLL_DIAG_ZEROS)

            CALL LINK_MESSAGE('ALLOCATE SPARSE MLLn ARRAYS')
            CALL ALLOCATE_SPARSE_MAT ( 'MLLn', NDOFL, NTERM_MLLn, SUBR_NAME )

            CALL LINK_MESSAGE('CONVERT SYM CRS MLL TO NONSYM CRS MLLn')
            CALL CRS_SYM_TO_CRS_NONSYM ( 'MLL', NDOFL, NTERM_MLL, I_MLL, J_MLL, MLL, 'MLLn', NTERM_MLLn, I_MLLn, J_MLLn, MLLn, 'Y' )

         ENDIF

      ELSE IF (SPARSTOR == 'NONSYM') THEN

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN

            CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR NONSYM STORAGE OF KLLD')
            NTERM_KLLDn = NTERM_KLLD
            CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )

            CALL LINK_MESSAGE('GET VALUES FOR NONSYM FORM OF KLLD')
            DO I=1,NDOFL+1
               I_KLLDn(I) = I_KLLD(I)
            ENDDO
            DO J=1,NTERM_KLLDn
               J_KLLDn(J) = J_KLLD(J)
                 KLLDn(J) =   KLLD(J)
            ENDDO

         ELSE

            CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR NONSYM STORAGE OF MLL')
            NTERM_MLLn = NTERM_MLL
            CALL ALLOCATE_SPARSE_MAT ( 'MLLn', NDOFL, NTERM_MLLn, SUBR_NAME )

            CALL LINK_MESSAGE('GET VALUES FOR NONSYM FORM OF MLL')
            DO I=1,NDOFL+1
               I_MLLn(I) = I_MLL(I)
            ENDDO
            DO J=1,NTERM_MLLn
               J_MLLn(J) = J_MLL(J)
                 MLLn(J) =   MLL(J)
            ENDDO

         ENDIF

      ELSE
         !      Error - incorrect SPARSTOR
         WRITE(ERR,932) SUBR_NAME, SPARSTOR
         WRITE(F06,932) SUBR_NAME, SPARSTOR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF

      IF (WRITE_MLL) THEN
         CALL WRITE_SPARSE_CRS ( ' MLLn', 'A ', 'A ', NTERM_MLLn, NDOFL, I_MLLn, J_MLLn, MLLn )
      ENDIF


! **********************************************************************************************************************************
      ! Identify modes-subcases for SOL 103. For BUCKLING / GEN CB MODEL the loop below runs exactly once with whatever scalars are
      ! already loaded from L1M (legacy behavior preserved). For SOL 103, we iterate per modes-subcase and capture each subcase's
      ! eigenresults into EIG_PARAMS(:); after the loop the scratch EIGEN_VAL/VEC/MODE_NUM/GEN_MASS are concatenated across subcases
      ! so the downstream L3A write and LINK5/LINK9 can see the full mode set with MODE_SUBCASE giving per-mode subcase attribution.
      NUM_MODES_SUBS = 0
      CANONICAL_ISUB = 0
      IS_BUCK_MULTI  = .FALSE.
      ! For SOL 105 LOADC populates IS_MODES_SUBCASE in lockstep with IS_BUCKLING_SUBCASE, so the same iteration logic
      ! drives both the modes-multi (SOL 103) and the buckling-multi (SOL 105) paths.
      IF (((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:8) == 'BUCKLING')) .AND. ALLOCATED(IS_MODES_SUBCASE)) THEN
         DO I=1,NSUB
            IF (IS_MODES_SUBCASE(I) == 'Y') THEN
               NUM_MODES_SUBS = NUM_MODES_SUBS + 1
               ! Canonical = first modes-subcase whose EIG_PARAMS%SID matches the L1M scalar CC_EIGR_SID
               IF ((CANONICAL_ISUB == 0) .AND. (EIG_PARAMS(I)%SID == CC_EIGR_SID)) THEN
                  CANONICAL_ISUB = I
               ENDIF
            ENDIF
         ENDDO
         ! Fallback: if no SID match found (shouldn't happen given LOADC sync), use first modes-subcase
         IF (CANONICAL_ISUB == 0) THEN
            DO I=1,NSUB
               IF (IS_MODES_SUBCASE(I) == 'Y') THEN
                  CANONICAL_ISUB = I
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      IF (CANONICAL_ISUB == 0) CANONICAL_ISUB = 1

      N_MODES_ITER = MAX(1, NUM_MODES_SUBS)

      ! Detect SOL 105 multi-buckling-subcase mode. The canonical preload was loaded by LINK5 step 1 (= first buckling subcase's
      ! STATSUB_REF) so iter 1 will not need to rebuild KLLD; iter 2+ may swap UG_COL and re-run REBUILD_KLLD_FROM_KGGD.
      IS_BUCK_MULTI = ((SOL_NAME(1:8) == 'BUCKLING') .AND. (NUM_BUCKLING_SUBS > 1))
      ! Initialize to 0 (unknown) so ITER=1 always reloads UG_COL from L5A. We cannot assume UG_COL holds the canonical
      ! preload here because LINK9 (LOAD_ISTEP=1) may have advanced past multiple preload vectors, leaving UG_COL pointing
      ! at the last one. LINK2 (LOAD_ISTEP=2) then built KLLD from that wrong UG_COL; ITER=1 must correct it.
      CURRENT_PRELOAD_ISUB = 0

      ! For multi-iter MODES solves we need to preserve KLL across iterations because EIG_LANCZOS_ARPACK destructively
      ! deallocates KLL mid-solve. Snapshot the CSR triple once here; restore at the head of iterations 2+.
      IF (N_MODES_ITER > 1) THEN
         NTERM_KLL_BAK = NTERM_KLL
         ALLOCATE(I_KLL_BAK(NDOFL+1))
         ALLOCATE(J_KLL_BAK(NTERM_KLL))
         ALLOCATE(KLL_BAK(NTERM_KLL))
         I_KLL_BAK = I_KLL
         J_KLL_BAK = J_KLL
         KLL_BAK   = KLL
      ENDIF

      ! Modes-subcase solver loop
m_lp: DO ITER = 1, N_MODES_ITER

         ! Restore KLL from shadow at the start of iterations 2+ (EIG_LANCZOS_ARPACK deallocated it during iter 1)
         IF (ITER > 1) THEN
            NTERM_KLL = NTERM_KLL_BAK
            IF (.NOT. ALLOCATED(I_KLL)) ALLOCATE(I_KLL(NDOFL+1))
            IF (.NOT. ALLOCATED(J_KLL)) ALLOCATE(J_KLL(NTERM_KLL))
            IF (.NOT. ALLOCATED(KLL))   ALLOCATE(KLL(NTERM_KLL))
            I_KLL = I_KLL_BAK
            J_KLL = J_KLL_BAK
            KLL   = KLL_BAK
         ENDIF

         ! For SOL 105 multi-buckling-subcase, iter>1 needs KLLD/KLLDn rebuilt (the inline KLLD/KLLDn deallocs at the end of
         ! the previous iter freed them). If the target subcase's STATSUB preload differs from the one currently loaded, also
         ! reload UG_COL from L5A before rebuilding. CUR_ISUB is determined just below from IS_MODES_SUBCASE / ITER; for the
         ! SOL 105 path that mapping is identical because LOADC sets IS_MODES_SUBCASE == IS_BUCKLING_SUBCASE.
         ! NOTE: ITER==1 also does an explicit reload so that LINK9 advancing past multiple static preload vectors in L5A
         !       cannot leave UG_COL pointing at the wrong preload for the first buckling subcase.
         IF (IS_BUCK_MULTI) THEN
            ! Map ITER -> CUR_ISUB locally so we can resolve TARGET_PRELOAD before the canonical scalar reload below
            KCNT = 0
            CUR_ISUB = CANONICAL_ISUB
            DO I=1,NSUB
               IF (IS_MODES_SUBCASE(I) == 'Y') THEN
                  KCNT = KCNT + 1
                  IF (KCNT == ITER) THEN
                     CUR_ISUB = I
                     EXIT
                  ENDIF
               ENDIF
            ENDDO
            TARGET_PRELOAD = 0
            IF (ALLOCATED(EIG_PARAMS)) TARGET_PRELOAD = EIG_PARAMS(CUR_ISUB)%STATSUB_REF
            IF ((TARGET_PRELOAD > 0) .AND. (TARGET_PRELOAD /= CURRENT_PRELOAD_ISUB)) THEN
               CALL LINK_MESSAGE('RELOAD UG_COL FROM L5A FOR NEXT PRELOAD')
               CALL DEALLOCATE_COL_VEC ( 'UG_COL' )
               CALL ALLOCATE_COL_VEC ( 'UG_COL', NDOFG, SUBR_NAME )
               IERR_RELOAD = 0
               CALL READ_L5A_UG_FOR_SUBCASE ( TARGET_PRELOAD, IERR_RELOAD )
               IF (IERR_RELOAD /= 0) THEN
                  WRITE(ERR,9994) SUBR_NAME, TARGET_PRELOAD, IERR_RELOAD
                  WRITE(F06,9994) SUBR_NAME, TARGET_PRELOAD, IERR_RELOAD
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
               CURRENT_PRELOAD_ISUB = TARGET_PRELOAD
               ! For ITER=1 the preload changed: need to rebuild KLLD (LINK2 built it from the wrong UG_COL)
               IF (ITER == 1) THEN
                  CALL LINK_MESSAGE('REBUILD KLLD FROM CORRECTED UG_COL (STATSUB ITER=1)')
                  CALL DEALLOCATE_SPARSE_MAT ( 'KLLDn' )  ! LINK2 left KLLDn allocated; free before REBUILD+realloc
                  CALL REBUILD_KLLD_FROM_KGGD
                  IF      (SPARSTOR == 'SYM   ') THEN
                     CALL SPARSE_MAT_DIAG_ZEROS ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, NUM_KLLD_DIAG_ZEROS )
                     NTERM_KLLDn = 2*NTERM_KLLD - (NDOFL - NUM_KLLD_DIAG_ZEROS)
                     CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )
                     CALL CRS_SYM_TO_CRS_NONSYM ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, KLLD, 'KLLDn', NTERM_KLLDn,          &
                                                  I_KLLDn, J_KLLDn, KLLDn, 'Y' )
                  ELSE IF (SPARSTOR == 'NONSYM') THEN
                     NTERM_KLLDn = NTERM_KLLD
                     CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )
                     DO I=1,NDOFL+1
                        I_KLLDn(I) = I_KLLD(I)
                     ENDDO
                     DO J=1,NTERM_KLLDn
                        J_KLLDn(J) = J_KLLD(J)
                          KLLDn(J) =   KLLD(J)
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF
            IF (ITER > 1) THEN
               CALL LINK_MESSAGE('REBUILD KLLD FROM CURRENT UG_COL (STATSUB)')
               CALL REBUILD_KLLD_FROM_KGGD
               ! Redo KLLDn conversion / copy (mirrors the pre-loop SPARSTOR block for BUCKLING)
               IF      (SPARSTOR == 'SYM   ') THEN
                  CALL SPARSE_MAT_DIAG_ZEROS ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, NUM_KLLD_DIAG_ZEROS )
                  NTERM_KLLDn = 2*NTERM_KLLD - (NDOFL - NUM_KLLD_DIAG_ZEROS)
                  CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )
                  CALL CRS_SYM_TO_CRS_NONSYM ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, KLLD, 'KLLDn', NTERM_KLLDn,             &
                                               I_KLLDn, J_KLLDn, KLLDn, 'Y' )
               ELSE IF (SPARSTOR == 'NONSYM') THEN
                  NTERM_KLLDn = NTERM_KLLD
                  CALL ALLOCATE_SPARSE_MAT ( 'KLLDn', NDOFL, NTERM_KLLDn, SUBR_NAME )
                  DO I=1,NDOFL+1
                     I_KLLDn(I) = I_KLLD(I)
                  ENDDO
                  DO J=1,NTERM_KLLDn
                     J_KLLDn(J) = J_KLLD(J)
                       KLLDn(J) =   KLLD(J)
                  ENDDO
               ENDIF
            ENDIF
         ENDIF

         ! Determine which subcase this iteration solves
         IF (NUM_MODES_SUBS == 0) THEN
            CUR_ISUB = CANONICAL_ISUB                      ! BUCKLING / GEN CB MODEL / single-shot legacy
         ELSE
            KCNT = 0
            CUR_ISUB = CANONICAL_ISUB
            DO I=1,NSUB
               IF (IS_MODES_SUBCASE(I) == 'Y') THEN
                  KCNT = KCNT + 1
                  IF (KCNT == ITER) THEN
                     CUR_ISUB = I
                     EXIT
                  ENDIF
               ENDIF
            ENDDO

            ! Load EIG_* scalars from EIG_PARAMS(CUR_ISUB) so the eigensolvers see this subcase's params
            EIG_SID          = EIG_PARAMS(CUR_ISUB)%SID
            EIG_METH         = EIG_PARAMS(CUR_ISUB)%METHOD
            EIG_NORM         = EIG_PARAMS(CUR_ISUB)%NORM
            EIG_GRID         = EIG_PARAMS(CUR_ISUB)%GRID
            EIG_COMP         = EIG_PARAMS(CUR_ISUB)%COMP
            EIG_FRQ1         = EIG_PARAMS(CUR_ISUB)%FRQ1
            EIG_FRQ2         = EIG_PARAMS(CUR_ISUB)%FRQ2
            EIG_N1           = EIG_PARAMS(CUR_ISUB)%N1
            EIG_N2           = EIG_PARAMS(CUR_ISUB)%N2
            EIG_NCVFACL      = EIG_PARAMS(CUR_ISUB)%NCVFACL
            EIG_MSGLVL       = EIG_PARAMS(CUR_ISUB)%MSGLVL
            EIG_MODE         = EIG_PARAMS(CUR_ISUB)%MODE
            EIG_VECS         = EIG_PARAMS(CUR_ISUB)%VECS
            EIG_CRIT         = EIG_PARAMS(CUR_ISUB)%CRIT
            EIG_SIGMA        = EIG_PARAMS(CUR_ISUB)%SIGMA
            EIG_LAP_MAT_TYPE = EIG_PARAMS(CUR_ISUB)%LAP_MAT_TYPE
         ENDIF

         ! Solve eigenvalue problem
         IF ((EIG_METH(1:3) == 'GIV') .OR. (EIG_METH(1:4) == 'MGIV')) THEN
            CALL EIG_GIV_MGIV

         ELSE IF (EIG_METH(1:3) == 'INV') THEN
            CALL EIG_INV_PWR

         ELSE IF (EIG_METH(1:7) == 'LANCZOS') THEN
            ! Use adaptive version if frequency range specified and not BUCKLING/GEN CB MODEL
            IF ((EIG_FRQ2 > EPS1) .AND. (SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')) THEN
               CALL EIG_LANCZOS_ARPACK_ADAPTIVE
            ELSE
               CALL EIG_LANCZOS_ARPACK
            ENDIF

         ELSE

            WRITE(ERR,4005) SUBR_NAME, EIG_METH
            WRITE(F06,4005) SUBR_NAME, EIG_METH
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit

         ENDIF

         ! For BUCKLING (always single-iter) we can deallocate KLLD inline. For MODES we keep MLL alive across iterations
         ! and dealloc it after the loop, since each iteration's eigensolver needs MLL/MLLn for the eigenproblem.
         IF (SOL_NAME(1:12) /= 'GEN CB MODEL') THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLLD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KLLD' )
            ENDIF
         ENDIF

         ! Calc generalized masses and renorm eigenvectors to mass (users renorm is done in LINK5)
         NUM_FAIL_CRIT = 0
         MAXMIJ        = 0
         MIJ_ROW       = 0
         MIJ_COL       = 0

         CALL ALLOCATE_EIGEN1_MAT ( 'GEN_MASS', NUM_EIGENS, 1, SUBR_NAME )

         IF (NVEC > 0) THEN
                                                              ! Calc gen mass
            CALL LINK_MESSAGE('CALCULATE GENERALIZED MASS')
            CALL CALC_GEN_MASS

            IF (EIG_NORM == 'MASS') THEN
                                                              ! Renorm vecs to mass if user asked for 'MASS'.
               CALL LINK_MESSAGE('RENORMALIZE EIGENVECTORS TO UNIT GEN MASS')
               CALL RENORM_ON_MASS ( NVEC, EPS1 )
            ENDIF

         ELSE

            DO I=1,NUM_EIGENS
               GEN_MASS(I) = ZERO
            ENDDO

         ENDIF

         ! BUCKLING is single-iter, so KLLDn dealloc is safe here. MLLn is deferred to after the loop for MODES.
         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLLDn', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KLLDn' )
         ENDIF

         ! Capture this subcase's eigenresults into EIG_PARAMS(CUR_ISUB)
         IF (NUM_MODES_SUBS > 0) THEN
            IF (ALLOCATED(EIG_PARAMS(CUR_ISUB)%EIGEN_VAL))  DEALLOCATE(EIG_PARAMS(CUR_ISUB)%EIGEN_VAL)
            IF (ALLOCATED(EIG_PARAMS(CUR_ISUB)%MODE_NUM))   DEALLOCATE(EIG_PARAMS(CUR_ISUB)%MODE_NUM)
            IF (ALLOCATED(EIG_PARAMS(CUR_ISUB)%GEN_MASS))   DEALLOCATE(EIG_PARAMS(CUR_ISUB)%GEN_MASS)
            IF (ALLOCATED(EIG_PARAMS(CUR_ISUB)%EIGEN_VEC))  DEALLOCATE(EIG_PARAMS(CUR_ISUB)%EIGEN_VEC)
            ALLOCATE(EIG_PARAMS(CUR_ISUB)%EIGEN_VAL(NUM_EIGENS))
            ALLOCATE(EIG_PARAMS(CUR_ISUB)%MODE_NUM(NUM_EIGENS))
            ALLOCATE(EIG_PARAMS(CUR_ISUB)%GEN_MASS(NUM_EIGENS))
            ALLOCATE(EIG_PARAMS(CUR_ISUB)%EIGEN_VEC(NDOFL, MAX(1,NVEC)))
            EIG_PARAMS(CUR_ISUB)%EIGEN_VAL(1:NUM_EIGENS) = EIGEN_VAL(1:NUM_EIGENS)
            EIG_PARAMS(CUR_ISUB)%MODE_NUM(1:NUM_EIGENS)  = MODE_NUM(1:NUM_EIGENS)
            EIG_PARAMS(CUR_ISUB)%GEN_MASS(1:NUM_EIGENS)  = GEN_MASS(1:NUM_EIGENS)
            IF (NVEC > 0) THEN
               EIG_PARAMS(CUR_ISUB)%EIGEN_VEC(1:NDOFL,1:NVEC) = EIGEN_VEC(1:NDOFL,1:NVEC)
            ENDIF
            EIG_PARAMS(CUR_ISUB)%NUM_EIGENS    = NUM_EIGENS
            EIG_PARAMS(CUR_ISUB)%NVEC          = NVEC
            EIG_PARAMS(CUR_ISUB)%NUM_FAIL_CRIT = NUM_FAIL_CRIT
            EIG_PARAMS(CUR_ISUB)%MAXMIJ        = MAXMIJ
            EIG_PARAMS(CUR_ISUB)%MIJ_ROW       = MIJ_ROW
            EIG_PARAMS(CUR_ISUB)%MIJ_COL       = MIJ_COL
            NUM_EIGENS_SUB(CUR_ISUB)           = NUM_EIGENS
         ENDIF

         ! Write eigenvalue analysis summary to output file (per-subcase summary in the multi-METHOD case)
         IF ((EIG_NORM == 'MASS    ') .OR. (EIG_NORM == 'NONE')) THEN
            CALL LINK_MESSAGE('WRITE EIGENVALUE SUMMARY TO OUTFIL')
            CALL EIG_SUMMARY(CUR_ISUB)
         ENDIF

         ! If more iterations remain, free scratch eigen arrays and Lanczos workspaces so the next iteration can re-allocate cleanly
         IF (ITER < N_MODES_ITER) THEN
            CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )
            CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC' )
            CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM' )
            CALL DEALLOCATE_EIGEN1_MAT ( 'GEN_MASS' )
            CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )
            CALL DEALLOCATE_LAPACK_MAT ( 'BBAND' )
            CALL DEALLOCATE_LAPACK_MAT ( 'RFAC' )
         ENDIF

      ENDDO m_lp

      ! Deferred MLL/MLLn deallocations for the MODES path (BUCKLING already did KLLD/KLLDn inline above; GEN CB MODEL skips entirely)
      IF ((SOL_NAME(1:5) == 'MODES')) THEN
         IF (ALLOCATED(MLL))  THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MLL ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MLL' )
         ENDIF
         IF (ALLOCATED(MLLn)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MLLn ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MLLn' )
         ENDIF
      ENDIF

      ! Free the multi-iter KLL shadow
      IF (ALLOCATED(I_KLL_BAK)) DEALLOCATE(I_KLL_BAK)
      IF (ALLOCATED(J_KLL_BAK)) DEALLOCATE(J_KLL_BAK)
      IF (ALLOCATED(KLL_BAK))   DEALLOCATE(KLL_BAK)

      ! Concatenate per-subcase eigenresults into the scratch arrays so downstream L3A / LINK5 / LINK9 see the full mode set.
      ! Also populate MODE_SUBCASE so LINK9 can attribute each mode to its owning subcase.
      IF (NUM_MODES_SUBS > 1) THEN
         TOTAL_MODES_LOCAL = 0
         DO I=1,NSUB
            IF (IS_MODES_SUBCASE(I) == 'Y') TOTAL_MODES_LOCAL = TOTAL_MODES_LOCAL + NUM_EIGENS_SUB(I)
         ENDDO
         TOTAL_MODES = TOTAL_MODES_LOCAL

         CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )
         CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC' )
         CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM' )
         CALL DEALLOCATE_EIGEN1_MAT ( 'GEN_MASS' )
         CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', TOTAL_MODES_LOCAL, 1, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC', NDOFL, TOTAL_MODES_LOCAL, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM', TOTAL_MODES_LOCAL, 1, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'GEN_MASS', TOTAL_MODES_LOCAL, 1, SUBR_NAME )

         IF (ALLOCATED(MODE_SUBCASE)) DEALLOCATE(MODE_SUBCASE)
         ALLOCATE(MODE_SUBCASE(TOTAL_MODES_LOCAL))

         IDX = 0
         DO I=1,NSUB
            IF (IS_MODES_SUBCASE(I) /= 'Y') CYCLE
            DO KMODE = 1, NUM_EIGENS_SUB(I)
               IDX = IDX + 1
               EIGEN_VAL(IDX)        = EIG_PARAMS(I)%EIGEN_VAL(KMODE)
               MODE_NUM(IDX)         = IDX
               GEN_MASS(IDX)         = EIG_PARAMS(I)%GEN_MASS(KMODE)
               EIGEN_VEC(1:NDOFL,IDX) = EIG_PARAMS(I)%EIGEN_VEC(1:NDOFL,KMODE)
               MODE_SUBCASE(IDX)     = I
            ENDDO
         ENDDO
         NUM_EIGENS = TOTAL_MODES_LOCAL
         NVEC       = TOTAL_MODES_LOCAL
      ELSE
         ! Single-iter path: trivial MODE_SUBCASE mapping
         IF (ALLOCATED(MODE_SUBCASE)) DEALLOCATE(MODE_SUBCASE)
         ALLOCATE(MODE_SUBCASE(MAX(1, NUM_EIGENS)))
         MODE_SUBCASE = CANONICAL_ISUB
         TOTAL_MODES  = NUM_EIGENS
      ENDIF

      ! End-of-loop cleanup for the SOL 105 multi-buckling-subcase path. LINK1 step 2 deferred the MPC_IND_GRIDS dealloc so that
      ! REBUILD_KLLD_FROM_KGGD / BUILD_KGGD_FROM_UG could be re-invoked per iter. Free it here, plus the residual UG_COL.
      IF (IS_BUCK_MULTI) THEN
         CALL DEALLOCATE_MODEL_STUF ( 'MPC_IND_GRIDS' )
         CALL DEALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS' )
         CALL DEALLOCATE_MODEL_STUF ( 'SUBLOD' )
         CALL DEALLOCATE_COL_VEC ( 'UG_COL' )
      ENDIF

      ! Restore canonical-subcase EIG_* scalars and write L1M (only once, after the loop, so the file holds the global mode count
      ! and canonical params). LINK5 / LINK9 read L1M for the scalars and the global eigen list, then use EIG_PARAMS(MODE_SUBCASE(J))
      ! for per-mode subcase attribution where needed.
      IF (NUM_MODES_SUBS > 0) THEN
         EIG_SID          = EIG_PARAMS(CANONICAL_ISUB)%SID
         EIG_METH         = EIG_PARAMS(CANONICAL_ISUB)%METHOD
         EIG_NORM         = EIG_PARAMS(CANONICAL_ISUB)%NORM
         EIG_GRID         = EIG_PARAMS(CANONICAL_ISUB)%GRID
         EIG_COMP         = EIG_PARAMS(CANONICAL_ISUB)%COMP
         EIG_FRQ1         = EIG_PARAMS(CANONICAL_ISUB)%FRQ1
         EIG_FRQ2         = EIG_PARAMS(CANONICAL_ISUB)%FRQ2
         EIG_N1           = EIG_PARAMS(CANONICAL_ISUB)%N1
         EIG_N2           = EIG_PARAMS(CANONICAL_ISUB)%N2
         EIG_NCVFACL      = EIG_PARAMS(CANONICAL_ISUB)%NCVFACL
         EIG_MSGLVL       = EIG_PARAMS(CANONICAL_ISUB)%MSGLVL
         EIG_MODE         = EIG_PARAMS(CANONICAL_ISUB)%MODE
         EIG_VECS         = EIG_PARAMS(CANONICAL_ISUB)%VECS
         EIG_CRIT         = EIG_PARAMS(CANONICAL_ISUB)%CRIT
         EIG_SIGMA        = EIG_PARAMS(CANONICAL_ISUB)%SIGMA
         EIG_LAP_MAT_TYPE = EIG_PARAMS(CANONICAL_ISUB)%LAP_MAT_TYPE
         NUM_FAIL_CRIT    = EIG_PARAMS(CANONICAL_ISUB)%NUM_FAIL_CRIT
         MAXMIJ           = EIG_PARAMS(CANONICAL_ISUB)%MAXMIJ
         MIJ_ROW          = EIG_PARAMS(CANONICAL_ISUB)%MIJ_ROW
         MIJ_COL          = EIG_PARAMS(CANONICAL_ISUB)%MIJ_COL
      ENDIF
      CALL WRITE_L1M

      ! Open and set up file L3A (used to hold eigenvectors)
      CALL FILE_OPEN ( L3A, LINK3A, OUNT, 'REPLACE', L3A_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

      ! Write out computed eigenvectors to L3A
      CALL LINK_MESSAGE('WRITE EIGENVECTORS TO DISK FILE')
      DO J=1,NVEC
         DO I=1,NDOFL
           WRITE(L3A) EIGEN_VEC(I,J)
         ENDDO
      ENDDO
      CALL FILE_CLOSE ( L3A, LINK3A, 'KEEP' )

      ! Optional eigenvector debug output
      IF (DEBUG(43) == 1) THEN
         DO J=1,NVEC
            DO I=1,NDOFL
               EIGEN_VEC_COL(I) = EIGEN_VEC(I,J)
            ENDDO
            WRITE(F06,'(//,1X,''EIGENVECTOR'',I8/)') J
            CALL WRITE_VECTOR ('    A-SET EIGENVECTOR   ','DISPL',NDOFL, EIGEN_VEC_COL )
         ENDDO
      ENDIF

      ! Call OUTPUT4 processor to process output requests for OUTPUT4 matrices generated in this link
      IF (NUM_OU4_REQUESTS > 0) THEN
         CALL LINK_MESSAGE('WRITE OUTPUT4 NATRICES      ')
         WRITE(F06,*)
         CALL OUTPUT4_PROC ( SUBR_NAME )
      ENDIF

      ! Deallocate arrays
      CALL DEALLOCATE_LAPACK_MAT ( 'RFAC' )

      ! leave EIGEN_VAL until LINK9 since it may be needed there
!xx   CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )
      CALL DEALLOCATE_EIGEN1_MAT ( 'GEN_MASS' )
      CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC' )
      CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM' )

      CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )
      CALL DEALLOCATE_LAPACK_MAT ( 'BBAND' )

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

      ! Write LINK4 end to F06
      CALL OURTIM
      WRITE(F06,151) LINKNO

      ! Close files
      IF (( DEBUG(193) == 4) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK4' )
      ENDIF

      ! Write LINK4 end to screen
      WRITE(SC1,153) LINKNO

      RETURN

! **********************************************************************************************************************************
  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END',/)

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  153 FORMAT(  ' >> LINK',I3,' END')

  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

  999 FORMAT(' *ERROR   999: INCORRECT SOLUTION IN EXEC CONTROL. SHOULD BE "',A,'", BUT IS "',A,'"')

 4005 FORMAT(' *ERROR  4005: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CODE ONLY WRITTEN FOR METHOD = GIV, MGIV, OR LANCZOS BUT METHOD IS = ',A8)

 9101 FORMAT(1X,A,' =  ','"',A,'"')

 9102 FORMAT(1X,A,' =  ',I13)

 9103 FORMAT(1X,A,' =  ',1ES13.6)

 9994 FORMAT(' *ERROR  9994: SUBROUTINE ',A,' FAILED TO RELOAD UG_COL FROM FILE L5A FOR PRELOAD SUBCASE ',I8,                     &
                    ' (IOSTAT = ',I8,').')

 9998 FORMAT(' *ERROR  9998: COMM ',I3,' INDICATES UNSUCCESSFUL LINK ',I2,' COMPLETION.'                                           &
                    ,/,14X,' FATAL ERROR - CANNOT START LINK ',I2)

12345 FORMAT(A,10X,A)

99001 FORMAT(1X,6(1ES14.6))

! **********************************************************************************************************************************

   END SUBROUTINE LINK4

   END MODULE LINK4_MOD
