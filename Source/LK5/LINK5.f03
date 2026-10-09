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

   MODULE LINK5_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: LINK5

   CONTAINS

      SUBROUTINE LINK5

! LINK5 takes the L-set displacements solved for in LINK3 (statics) or LINK4 (eigenvalues) and builds it back up to the G-set.
! See Appendix B to the MYSTRAN User's Reference Guide for an explanation of how this is done.

! In addition, for Craig-Bampton model generation (SOL = GEN CB MODEL or 31), array PHIXA is expanded to G-set size

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, ERR, F06, L1H, L2A, L2E, L2F, L3A, L5A, L5B, SC1
      USE IOUNT1, ONLY                :  LINK1H, LINK2A, LINK2E, LINK2F, LINK3A, LINK5A, LINK5B
      USE IOUNT1, ONLY                :  L1H_MSG, L2A_MSG, L2E_MSG, L2F_MSG, L3A_MSG, L5A_MSG, L5B_MSG
      USE IOUNT1, ONLY                :  ERRSTAT, L1HSTAT, L2ESTAT, L2FSTAT, L3ASTAT
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, COMM, FATAL_ERR, LINKNO, MBUG, NDOFA, NDOFF, NDOFG, NDOFL, NDOFM,           &
                                         NDOFN, NDOFO, NDOFR, NDOFS, NDOFSE, NGRID, NSUB, NTERM_GMN, NTERM_GOA, NTERM_PO,          &
                                         NUM_CB_DOFS, NUM_EIGENS, NVEC, SOL_NAME, WARN_ERR, MODE_SUBCASE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  EIGNORM2, SUPINFO, SUPWARN
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE EIGEN_MATRICES_1 , ONLY     :  EIGEN_VAL, EIGEN_VEC, GEN_MASS, MODE_NUM
      USE FULL_MATRICES, ONLY         :  PHIZG_FULL
      USE SPARSE_MATRICES, ONLY       :  I_GMN, J_GMN, GMN, I_GOA, J_GOA, GOA
      USE OUTPUT4_MATRICES, ONLY      :  NUM_OU4_REQUESTS
      USE MISC_MATRICES, ONLY         :  UG_T123_MAT
      USE COL_VECS, ONLY              :  UG_COL, YSe, UO0_COL, UL_COL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE DOF_TABLES, ONLY            :  TDOF, TDOFI
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID, INV_GRID_SEQ, EIG_COMP, EIG_GRID, EIG_NORM, MAXMIJ, MIJ_COL, MIJ_ROW,      &
                                         EIG_PARAMS, IS_BUCKLING_SUBCASE, IS_MODES_SUBCASE

      USE DATE_TIME_UTILS, ONLY       :  OURDAT, OURTIM, TIME_INIT
      USE TEMP_FILE_READERS, ONLY     :  READ_L1A, READ_L1M
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE, WRITE_L1A
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE COL_VEC_LIFECYCLE, ONLY     :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
      USE MATRIX_FILE_IO, ONLY        :  READ_MATRIX_1
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_INQUIRE, FILE_OPEN, READERR, WRITE_FILNAM
      USE EIGEN_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_EIGEN1_MAT, DEALLOCATE_EIGEN1_MAT
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE MISC_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_MISC_MAT, DEALLOCATE_MISC_MAT
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_GRID_NUM_COMPS, GET_UG_123_IN_GRD_ORD
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_L1M
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE EIGEN_SUPPORT, ONLY         :  EIG_SUMMARY
      USE OUTPUT4_FILE_IO, ONLY       :  OUTPUT4_PROC
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_ARRAY_ALLOC_STAT, LINK_MESSAGE, LINK_MESSAGE_I, WRITE_ALLOC_MEM_TABLE
      USE BDF_FIELD_VALIDATION, ONLY  :  I4FLD
      USE RESTART_FILE_IO, ONLY       :  READ_L5A_UG_FOR_SUBCASE

      IMPLICIT NONE

      LOGICAL                         :: VEC_SIGN_CHG(NDOFL) ! Indicators of whether user wants to change sign of an eigenvector

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK5'
      CHARACTER( 1*BYTE)              :: CLOSE_IT          ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER( 8*BYTE)              :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER( 1*BYTE)              :: DO_IT             ! If 'Y' execute some code
      CHARACTER( 1*BYTE)              :: MIJ_COL_FOUND='N' ! 'Y' if MIJ_ROW is processed as a solution vector in this LINK
      CHARACTER( 1*BYTE)              :: MIJ_ROW_FOUND='N' ! 'Y' if MIJ_ROW is processed as a solution vector in this LINK
      CHARACTER( 1*BYTE)              :: READ_NTERM        ! 'Y' or 'N' Input to subr READ_MATRIX_1
      CHARACTER( 1*BYTE)              :: OPND              ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not
      CHARACTER( 1*BYTE)              :: READ_UO0          ! If 'Y' then read UO0 data from file L2F

      INTEGER(LONG)                   :: COL_NUM           ! Arg passed to subr BUILD_A_LR
      INTEGER(LONG)                   :: EIG_NORM_GSET_DOF ! A-set DOF no. for EIG_GRID/EIG_COMP
      INTEGER(LONG)                   :: EIGNORM2_ERR        ! Error indicator for reading param EIGNORM2 data
      INTEGER(LONG)                   :: G_SET_COL         ! Col number in TDOF, TDOFI where G-set DOF's exist
      INTEGER(LONG)                   :: I,J,K,L           ! DO loop indices
      INTEGER(LONG)                   :: I_ESUB            ! Subcase loop index for EIG_SUMMARY header writes
      INTEGER(LONG)                   :: IERROR            ! Error count
      INTEGER(LONG)                   :: IGRID             ! Internal grid numbER for EIG_GRID
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID_NUM is an physical grid, 1 if an SPOINT
      INTEGER(LONG)                   :: NUM_SOLNS    = 0  ! No. of solutions to process (e.g. NSUB for STATICS)
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: P_LINKNO          ! Prior LINK no's that should have run before this LINK can execute
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file

      REAL(DOUBLE)                    :: MIJ_COL_SCALE=ZERO! Scale fac for a col of gen mass matrix to renorm MAXMIJ from LINK4
      REAL(DOUBLE)                    :: MIJ_ROW_SCALE=ZERO! Scale fac for a col of gen mass matrix to renorm MAXMIJ from LINK4
      REAL(DOUBLE)                    :: PHI_SCALE_FAC     ! Scale factor that the eigenvector was renormalized to in subr RENORM

! **********************************************************************************************************************************
      LINKNO = 5

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

      IF (SOL_NAME(1:7) == 'STATICS') THEN

         P_LINKNO = 3
         IF (COMM(P_LINKNO) /= 'C') THEN
            WRITE(SC1,9998) P_LINKNO,P_LINKNO,LINKNO
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ELSE IF (SOL_NAME(1:5) == 'MODES') THEN

         P_LINKNO = 4
         IF (COMM(P_LINKNO) /= 'C') THEN
            WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
            WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN

         P_LINKNO = 6
         IF (COMM(P_LINKNO) /= 'C') THEN
            WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
            WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ELSE IF (SOL_NAME(1:8) == 'BUCKLING') THEN

         IF      (LOAD_ISTEP == 1) THEN
            P_LINKNO = 3
            IF (COMM(P_LINKNO) /= 'C') THEN
               WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
               WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
         ELSE IF (LOAD_ISTEP == 2) THEN
            P_LINKNO = 4
            IF (COMM(P_LINKNO) /= 'C') THEN
               WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
               WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
         ELSE
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,5002) SUBR_NAME, LOAD_ISTEP
            WRITE(F06,5002) SUBR_NAME, LOAD_ISTEP
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ELSE IF (SOL_NAME(1:8) == 'NLSTATIC') THEN

         P_LINKNO = 3
         IF (COMM(P_LINKNO) /= 'C') THEN
            WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
            WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ELSE

         WRITE(ERR,5001) SOL_NAME
         WRITE(F06,5001) SOL_NAME
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF

! **********************************************************************************************************************************
! Allocate arrays for YSe, GMN, GOA

      CALL LINK_MESSAGE('ALLOCATE SEVERAL ARRAYS')
      CALL ALLOCATE_SPARSE_MAT ( 'GMN', NDOFM, NTERM_GMN, SUBR_NAME )
      CALL ALLOCATE_SPARSE_MAT ( 'GOA', NDOFO, NTERM_GOA, SUBR_NAME )
      CALL ALLOCATE_COL_VEC ( 'YSe' , NDOFS, SUBR_NAME )

! Read GMN matrix if there are MPC's

      IF (NTERM_GMN > 0) THEN

         CALL LINK_MESSAGE('READ GMN MATRIX')
         READ_NTERM = 'Y'
         OPND       = 'N'
         CLOSE_IT   = 'Y'
         CALL READ_MATRIX_1 ( LINK2A, L2A, OPND, CLOSE_IT, 'KEEP', L2A_MSG, 'GMN', NTERM_GMN, READ_NTERM, NDOFM                    &
                            , I_GMN, J_GMN, GMN )
      ENDIF

! Read GOA matrix if there are omitted DOFs

      IF (NTERM_GOA > 0) THEN
         CALL LINK_MESSAGE('READ GOA MATRIX')
         READ_NTERM = 'Y'
         OPND       = 'N'
         CLOSE_IT   = 'Y'
!        CLOSE_STAT = L2ESTAT
         CLOSE_STAT = 'KEEP'
         CALL READ_MATRIX_1 ( LINK2E, L2E, OPND, CLOSE_IT, CLOSE_STAT, L2E_MSG, 'GOA', NTERM_GOA, READ_NTERM, NDOFO                &
                            , I_GOA, J_GOA, GOA )
      ENDIF

! Read enforced displ's

      IF (NDOFSE > 0) THEN

         IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                                  &
            ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN

            CALL FILE_OPEN ( L1H, LINK1H, OUNT, 'OLD', L1H_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

            CALL LINK_MESSAGE('READ YSe ENFORCED DISPLACEMENTS')

            IERROR = 0
            DO I=1,NDOFSE
               READ(L1H,IOSTAT=IOCHK) YSe(I)
               IF (IOCHK /= 0) THEN
                  IERROR = IERROR + 1
                  REC_NO = I
                  CALL READERR ( IOCHK, LINK1H, L1H_MSG, REC_NO, OUNT )
               ENDIF
            ENDDO
            IF (IERROR /= 0) THEN
               WRITE(ERR,9995) LINKNO,IERROR
               WRITE(F06,9995) LINKNO,IERROR
               CALL OUTA_HERE ( 'Y' )
            ENDIF

            IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN
               ! ensure L1H survives for the second round
               CALL FILE_CLOSE ( L1H, LINK1H, 'KEEP' )
            ELSE
               CALL FILE_CLOSE ( L1H, LINK1H, L1HSTAT )
            END IF

         ENDIF

      ENDIF

! Read eigenvalue data from file L1M if this is an eigenvalue problem. Check to see if user wants to change sign of any vector.
! EIGEN_VAL was not deallocated in LINK4 (see LINK4 comment 01/11/19) so we do not allocate it here anymore

      DO_IT = 'N'
      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
         DO_IT = 'Y'
      ENDIF
      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         DO_IT = 'Y'
      ENDIF

      IF (DO_IT == 'Y') THEN
         IERROR = 0
         CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM' , NUM_EIGENS, 1, SUBR_NAME )
!xx      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', NUM_EIGENS, 1, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'GEN_MASS' , NUM_EIGENS, 1, SUBR_NAME )
         CALL READ_L1M ( IERROR )

         IF (IERROR /= 0) THEN
            WRITE(ERR,9995) LINKNO,IERROR
            WRITE(F06,9995) LINKNO,IERROR
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ENDIF

! Open file that has L-set displs, or eigenvectors ('MODES') or PHIZL ('GEN CB MODEL')

      IF (NDOFL > 0) THEN
         CALL FILE_OPEN ( L3A, LINK3A, OUNT, 'OLD', L3A_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
      ENDIF

! Open file for writing displs to.

      CALL FILE_CLOSE ( L5A, LINK5A, 'KEEP' )
      CALL FILE_OPEN  ( L5A, LINK5A, OUNT, 'REPLACE', L5A_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

! Open file that has UO0

      IF (NTERM_PO > 0) THEN
         CALL FILE_OPEN ( L2F, LINK2F, OUNT, 'OLD', L2F_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
      ENDIF

! Set NUM_SOLNS for use in loop (below) to get outputs for each subcase/solution vector

      IF      ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
         NUM_SOLNS = NSUB
      ELSE IF (SOL_NAME(1:5) == 'MODES') THEN
         NUM_SOLNS = NVEC
      ELSE IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         IF (LOAD_ISTEP == 1) THEN
            ! Process every subcase's static solution so file LINK5A holds one full UG_COL per subcase.
            ! The step-2 KGGD assembly (and any per-buckling-subcase preload selection done via STATSUB) seeks into
            ! L5A by subcase index and reads back the appropriate preload UG, so we must have all NSUB columns on disk.
            NUM_SOLNS = NSUB
         ELSE IF (LOAD_ISTEP == 2) THEN
            NUM_SOLNS = NVEC
         ENDIF
      ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
         NUM_SOLNS = NUM_CB_DOFS
      ENDIF

! Allocate memory to EIGEN_VEC so that we can write the G-set eigenvectors to it for possible OUTPUT4 later

      IF (SOL_NAME(1:5) == 'MODES') THEN
         CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC', NDOFG, NUM_SOLNS, SUBR_NAME )
      ENDIF

! Allocate memory for CB matrix PHIZG_FULL

      IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
         CALL ALLOCATE_FULL_MAT ( 'PHIZG_FULL', NDOFG, NUM_CB_DOFS, SUBR_NAME )
      ENDIF

! Check for renorm = MAX or POINT.

      DO_IT = 'N'
      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
         DO_IT = 'Y'
      ENDIF
      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         DO_IT = 'Y'
      ENDIF

      IF (DO_IT == 'Y') THEN
         EIG_NORM_GSET_DOF = 0
         IGRID = 0
         IF (EIG_NORM == 'POINT   ') THEN                  ! User requested to renormalize eigenvectors on POINT

            CALL TDOF_COL_NUM ( 'G ',  G_SET_COL )
i_do:       DO I=1,NDOFG
               IF (TDOF(I,1) == EIG_GRID) THEN
                  IGRID = TDOF(I,3)
                  EIG_NORM_GSET_DOF = TDOF(I,G_SET_COL) + EIG_COMP - 1
                  EXIT i_do
               ENDIF
            ENDDO i_do

            IF (IGRID > 0) THEN
               WRITE(ERR,5102) EIG_GRID, EIG_COMP
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,5102) EIG_GRID, EIG_COMP
               ENDIF
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,5111) EIG_GRID
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,5111) EIG_GRID
               ENDIF
               PHI_SCALE_FAC = ONE
            ENDIF

         ELSE IF (EIG_NORM == 'MAX     ') THEN             ! User requested to renormalize eigenvectors on MAX

            WRITE(ERR,5103)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,5103)
            ENDIF

         ELSE                                              ! No renorm needed here (was done in LINK4 if not POINT or MASS)

            WRITE(ERR,5104)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,5104)
            ENDIF

         ENDIF

      ENDIF

      IF ((EIGNORM2 == 'Y') .AND. (NVEC > 0)) THEN
         CALL READ_EIGNORM2
      ENDIF

! **********************************************************************************************************************************
      READ_UO0 = 'Y'
      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         READ_UO0 = 'N'
      ENDIF

! Begin loop for reading L-set displs and building up to G-set displs one subcase/solution vector at a time

j_do: DO J = 1,NUM_SOLNS

         CALL ALLOCATE_COL_VEC ('UL_COL', NDOFL, SUBR_NAME)! Allocate array UL_COL

                                                           ! Read UL displs for the current subcase/vector from LINK3A
         IF     ((SOL_NAME(1: 7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
            CALL LINK_MESSAGE_I('READ  L-SET DISPLACEMENTS                      Subcase', J)
         ELSE IF (SOL_NAME(1: 5) == 'MODES') THEN
            CALL LINK_MESSAGE_I('READ  L-SET EIGENVECTORS                       Vector', J)
         ELSE IF (SOL_NAME(1: 8) == 'BUCKLING') THEN
            IF (LOAD_ISTEP == 1) THEN
               CALL LINK_MESSAGE_I('READ  L-SET DISPLACEMENTS                      Subcase', J)
            ELSE IF (LOAD_ISTEP == 2) THEN
               CALL LINK_MESSAGE_I('READ  L-SET EIGENVECTORS                       Vector', J)
            ENDIF
         ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
            CALL LINK_MESSAGE_I('READ  L-SET CB VECTORS (PHIZL)                 CB vec', J)
         ENDIF

         REC_NO = 0
         IERROR = 0
         DO I=1,NDOFL
            REC_NO = REC_NO + 1
            READ(L3A,IOSTAT=IOCHK) UL_COL(I)               ! For CB, a col of PHIZL. So UG_COL, calc'd in this subr, is a PHIZG col
            IF (IOCHK /= 0) THEN
               CALL READERR ( IOCHK, LINK3A, L3A_MSG, REC_NO, OUNT )
               IERROR = IERROR + 1
            ENDIF
         ENDDO

         IF (IERROR /= 0) THEN
            WRITE(ERR,9995) LINKNO,IERROR
            WRITE(F06,9995) LINKNO,IERROR
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         ! Multi-METHOD MODES: switch EIG_NORM / EIG_GRID / EIG_COMP and recompute EIG_NORM_GSET_DOF to this mode's
         ! owning subcase before the per-mode renorm call below. For legacy single-METHOD this is a no-op (EIG_PARAMS
         ! entries all reference the canonical subcase's params). Nested IFs avoid evaluating SIZE/index on
         ! MODE_SUBCASE when it is unallocated (gfortran does not guarantee short-circuit evaluation of .AND.).
         IF (SOL_NAME(1:5) == 'MODES') THEN
            IF (ALLOCATED(MODE_SUBCASE)) THEN
               IF (J <= SIZE(MODE_SUBCASE)) THEN
                  IF ((MODE_SUBCASE(J) >= 1) .AND. ALLOCATED(EIG_PARAMS)) THEN
                     IF (EIG_PARAMS(MODE_SUBCASE(J))%SID /= 0) THEN
                        EIG_NORM = EIG_PARAMS(MODE_SUBCASE(J))%NORM
                        EIG_GRID = EIG_PARAMS(MODE_SUBCASE(J))%GRID
                        EIG_COMP = EIG_PARAMS(MODE_SUBCASE(J))%COMP
                        IF (EIG_NORM == 'POINT   ') THEN
                           EIG_NORM_GSET_DOF = 0
                           CALL TDOF_COL_NUM ( 'G ',  G_SET_COL )
                           DO I=1,NDOFG
                              IF (TDOF(I,1) == EIG_GRID) THEN
                                 EIG_NORM_GSET_DOF = TDOF(I,G_SET_COL) + EIG_COMP - 1
                                 EXIT
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
                                                           ! Build UA from UL and UR
         CALL ALLOCATE_COL_VEC ( 'UA_COL', NDOFA, SUBR_NAME )
         CALL ALLOCATE_COL_VEC ( 'UR_COL', NDOFR, SUBR_NAME )
         CALL LINK_MESSAGE_I('BUILD UA DISPLS FROM UL, UR:                      "', J)
         COL_NUM = 0
         IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
            IF ((J > NDOFR+nvec) .AND. (J <= NUM_CB_DOFS)) THEN
               COL_NUM = J
            ENDIF
         ENDIF
         CALL BUILD_A_LR ( COL_NUM )
         CALL DEALLOCATE_COL_VEC ( 'UL_COL' )
         CALL DEALLOCATE_COL_VEC ( 'UR_COL' )

                                                           ! Solve for UO and Build UF from UA and UO
         CALL ALLOCATE_COL_VEC ( 'UF_COL' , NDOFF, SUBR_NAME )
         CALL ALLOCATE_COL_VEC ( 'UO_COL' , NDOFO, SUBR_NAME )
         CALL ALLOCATE_COL_VEC ( 'UO0_COL', NDOFO, SUBR_NAME )
         CALL LINK_MESSAGE_I('BUILD UF DISPLS FROM UA, UO:                      "', J)
         IF (READ_UO0 == 'Y') THEN
            IF (NDOFO > 0) THEN
               IF (NTERM_PO > 0) THEN
                  CALL LINK_MESSAGE_I('  READ UO0 DISPLS,                                "', J)

                  IERROR = 0
                  DO I=1,NDOFO
                     READ(L2F,IOSTAT=IOCHK) UO0_COL(I)
                     IF (IOCHK /= 0) THEN
                        REC_NO = I+1
                        CALL READERR ( IOCHK, LINK2F, L2F_MSG, REC_NO, OUNT )
                        IERROR = IERROR + 1
                     ENDIF
                  ENDDO
                  IF (IERROR /= 0) THEN
                     WRITE(ERR,9995) LINKNO,IERROR
                     WRITE(F06,9995) LINKNO,IERROR
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
               ELSE
                  DO I=1,NDOFO
                     UO0_COL(I) = ZERO
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
         CALL BUILD_F_AO
         CALL DEALLOCATE_COL_VEC ( 'UA_COL' )
         CALL DEALLOCATE_COL_VEC ( 'UO_COL' )
         CALL DEALLOCATE_COL_VEC ( 'UO0_COL' )
                                                           ! Build UN from UF and US
         CALL ALLOCATE_COL_VEC ( 'UN_COL', NDOFN , SUBR_NAME)
         CALL ALLOCATE_COL_VEC ( 'US_COL', NDOFS, SUBR_NAME )
         CALL LINK_MESSAGE_I('BUILD UN DISPLS FROM UF, US:                      "', J)
         CALL BUILD_N_FS
         CALL DEALLOCATE_COL_VEC ( 'UF_COL' )
         CALL DEALLOCATE_COL_VEC ( 'US_COL' )
                                                           ! Build UG from UN and UM
         CALL DEALLOCATE_COL_VEC ( 'UG_COL' )
         CALL ALLOCATE_COL_VEC ( 'UG_COL', NDOFG, SUBR_NAME )
         CALL ALLOCATE_COL_VEC ( 'UM_COL', NDOFM, SUBR_NAME )
         CALL LINK_MESSAGE_I('BUILD UG DISPLS FROM UN, UM:                      "', J)
         CALL BUILD_G_NM
         CALL DEALLOCATE_COL_VEC ( 'UN_COL' )
         CALL DEALLOCATE_COL_VEC ( 'UM_COL' )

         IERROR = 0
!xx      IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
            IF (ALLOCATED(UG_T123_MAT)) THEN
               CALL DEALLOCATE_MISC_MAT ( 'UG_T123_MAT' )
            ENDIF
            CALL ALLOCATE_MISC_MAT ( 'UG_T123_MAT', NGRID, 3, SUBR_NAME )
            CALL GET_UG_123_IN_GRD_ORD ( IERROR )
!xx      ENDIF
         IF (IERROR /= 0) THEN
            WRITE(ERR,9995) LINKNO,IERROR
            WRITE(F06,9995) LINKNO,IERROR
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         IF      (SOL_NAME(1: 5) == 'MODES'       ) THEN   ! For modes, write UG_COL to jth col of EIGEN_VEC (for later OUTPUT4)
            DO I=1,NDOFG
               EIGEN_VEC(I,J)  = UG_COL(I)
            ENDDO
         ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN   ! For CB   , write UG_COL to jth col of PHiZG     (for later OUTPUT4)
            DO I=1,NDOFG
               PHIZG_FULL(I,J) = UG_COL(I)
            ENDDO
         ENDIF
                                                           ! Renorm eigenvecs, if requested
         IF (J <= NVEC) THEN

            DO_IT = 'N'
            IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
               DO_IT = 'Y'
            ENDIF
            IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
               DO_IT = 'Y'
            ENDIF

            IF (DO_IT == 'Y') THEN
               IF ((EIG_NORM == 'POINT   ') .OR. (EIG_NORM == 'MAX     ')) THEN
                  CALL RENORM (J, EIG_GRID, EIG_COMP, EIG_NORM, EIG_NORM_GSET_DOF, GEN_MASS(J), PHI_SCALE_FAC )
                  IF      (J == MIJ_ROW) THEN
                     MIJ_ROW_FOUND = 'Y'
                     MIJ_ROW_SCALE = PHI_SCALE_FAC
                  ELSE IF (J == MIJ_COL) THEN
                     MIJ_COL_FOUND = 'Y'
                     MIJ_COL_SCALE = PHI_SCALE_FAC
                  ENDIF
                  CALL WRITE_L1M                           ! Need to update GEN_MASS in L1M if vecs are renorm'd here
               ENDIF

            ENDIF

         ENDIF
                                                           ! See if user requested eigenvector sign changes
         IF ((EIGNORM2 == 'Y') .AND. (EIGNORM2_ERR == 0)) THEN
            IF (VEC_SIGN_CHG(J)) THEN
               DO I=1,NDOFG
                  UG_COL(I) = -UG_COL(I)
               ENDDO
            ENDIF
         ENDIF

                                                           ! Write UG displs for this subcase to file LINK5A
         CALL LINK_MESSAGE_I('WRITE UG DISPLS TO FILE,                          "', J)
         WRITE(SC1, * )                                    ! Separator between UG_COL calcs
         DO I=1,NDOFG
            WRITE(L5A) UG_COL(I)                           ! For CB this is a col of PHIZG (which is never processed as an array)
         ENDDO

      ENDDO j_do                                           ! End of loop on NUM_SOLNS

! For SOL 105 step 1, the j_do loop above iterates over every subcase, including the buckling subcases (which carry no load,
! so their UG_COL is zero). Without intervention the UG_COL left in memory after the loop is whichever subcase happened to be
! processed last. The step-2 LINK1 ESP path uses that residual UG_COL to assemble the differential stiffness KGGD. To preserve
! legacy single-preload behaviour (and to give multi-buckling decks a sensible default until the explicit per-buckling-subcase
! preload selection lands in a later phase), we reload UG_COL with the canonical preload subcase's column from L5A.
! The canonical choice is the STATSUB_REF resolved for the first buckling subcase. If for some reason none of that information
! is available (defensive fallback only -- LOADC always resolves it for valid decks) we leave UG_COL untouched.

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN
         BUCKLING_PRELOAD_RELOAD : BLOCK
            INTEGER(LONG) :: I_BUCK, ISUB_PRELOAD, IERR_RELOAD
            ISUB_PRELOAD = 0
            IF (ALLOCATED(IS_BUCKLING_SUBCASE) .AND. ALLOCATED(EIG_PARAMS)) THEN
               DO I_BUCK = 1, NSUB
                  IF (IS_BUCKLING_SUBCASE(I_BUCK) == 'Y') THEN
                     IF (EIG_PARAMS(I_BUCK)%STATSUB_REF > 0) THEN
                        ISUB_PRELOAD = EIG_PARAMS(I_BUCK)%STATSUB_REF
                        EXIT
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            IF (ISUB_PRELOAD > 0) THEN
               CALL DEALLOCATE_COL_VEC ( 'UG_COL' )
               CALL ALLOCATE_COL_VEC ( 'UG_COL', NDOFG, SUBR_NAME )
               IERR_RELOAD = 0
               CALL READ_L5A_UG_FOR_SUBCASE ( ISUB_PRELOAD, IERR_RELOAD )
               IF (IERR_RELOAD /= 0) THEN
                  WRITE(ERR,9995) LINKNO, IERR_RELOAD
                  WRITE(F06,9995) LINKNO, IERR_RELOAD
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ENDIF
         END BLOCK BUCKLING_PRELOAD_RELOAD
      ENDIF

! If CB soln, expand PHIXA to G-set size and write to file unit L5B

     IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                                                           ! Open file for writing cols of PHIXG
         CALL FILE_OPEN ( L5B, LINK5B, OUNT, 'REPLACE', L5B_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

         CALL DEALLOCATE_COL_VEC ( 'UG_COL' )
         CALL EXPAND_PHIXA_TO_PHIXG                        ! Expand PHIXA to PHIXG and write cols to file L5B
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PHIXA', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PHIXA' )

      ENDIF

! Code to write PHIZG in full format

      IF ((DEBUG(55) == 2) .OR. (DEBUG(55) == 3)) THEN
         IF (SOL_NAME == 'GEN CB MODEL') THEN
            WRITE(F06,*)
            WRITE(F06,99885)
            WRITE(F06,99886) (J,J=1,NUM_CB_DOFS)
            L = 0
            DO I=1,NGRID
               CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
               DO K=1,NUM_COMPS
                  L = L + 1
                  IF (K == 1) THEN
                     WRITE(F06,99887) TDOFI(L,1), K, (PHIZG_FULL(L,J),J=1,NUM_CB_DOFS)
                  ELSE
                     WRITE(F06,99888)             K, (PHIZG_FULL(L,J),J=1,NUM_CB_DOFS)
                  ENDIF
               ENDDO
               WRITE(F06,*)
            ENDDO
            WRITE(F06,*)
         ENDIF
      ENDIF

! **********************************************************************************************************************************
! Write eigen analysis summary if we have renormed vectors here in LINK5

      DO_IT = 'N'
      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
         DO_IT = 'Y'
      ENDIF
      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         DO_IT = 'Y'
      ENDIF

      IF (DO_IT == 'Y') THEN
         IF ((EIG_NORM == 'POINT   ') .OR. (EIG_NORM == 'MAX     ')) THEN

            IF ((MIJ_ROW_FOUND == 'Y') .AND. (MIJ_COL_FOUND == 'Y')) THEN
               MAXMIJ = MAXMIJ/(MIJ_ROW_SCALE * MIJ_COL_SCALE)
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,5101) LINKNO
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,5101) LINKNO
               ENDIF
            ENDIF

            IF (SOL_NAME(1:5) == 'MODES') THEN
               IF (ALLOCATED(EIG_PARAMS) .AND. ALLOCATED(IS_MODES_SUBCASE)) THEN
                  DO I_ESUB = 1, NSUB
                     IF (IS_MODES_SUBCASE(I_ESUB) == 'Y') THEN
                        IF ((EIG_PARAMS(I_ESUB)%NORM == 'POINT   ') .OR.                           &
                            (EIG_PARAMS(I_ESUB)%NORM == 'MAX     ')) THEN
                           CALL EIG_SUMMARY(I_ESUB)
                        ENDIF
                     ENDIF
                  ENDDO
               ELSE
                  CALL EIG_SUMMARY(1)
               ENDIF
            ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
               IF (ALLOCATED(EIG_PARAMS) .AND. ALLOCATED(IS_BUCKLING_SUBCASE)) THEN
                  DO I_ESUB = 1, NSUB
                     IF (IS_BUCKLING_SUBCASE(I_ESUB) == 'Y') THEN
                        IF ((EIG_PARAMS(I_ESUB)%NORM == 'POINT   ') .OR.                           &
                            (EIG_PARAMS(I_ESUB)%NORM == 'MAX     ')) THEN
                           CALL EIG_SUMMARY(I_ESUB)
                        ENDIF
                     ENDIF
                  ENDDO
               ELSE
                  CALL EIG_SUMMARY(1)
               ENDIF
            ELSE
               CALL EIG_SUMMARY(1)                       ! GEN CB MODEL: no subcase cards
            ENDIF
         ENDIF
      ENDIF

! Close files

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN
         CLOSE_STAT = 'KEEP'
      ELSE
!xx      CLOSE_STAT = L2FSTAT
         CLOSE_STAT = 'KEEP'
      ENDIF
      CALL FILE_CLOSE ( L2F, LINK2F, CLOSE_STAT )

      IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
         CALL FILE_CLOSE ( L3A, LINK3A, 'KEEP' )
      ELSE
         CALL FILE_CLOSE ( L3A, LINK3A, L3ASTAT )
      ENDIF

      CALL FILE_CLOSE ( L5A, LINK5A, 'KEEP' )
      CALL FILE_CLOSE ( L5B, LINK5B, 'KEEP' )

! Call OUTPUT4 processor to process output requests for OUTPUT4 matrices generated in this link

      IF (NUM_OU4_REQUESTS > 0) THEN
         CALL LINK_MESSAGE('WRITE OUTPUT4 NATRICES      ')
         WRITE(F06,*)
         CALL OUTPUT4_PROC ( SUBR_NAME )
      ENDIF

! Deallocate arrays (except EIGEN_VAL, may be needed later)

!xx   WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GMN      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'GMN' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GOA      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'GOA' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MODE_NUM ', CR13  ;   CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM' )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate EIGEN_VAL', CR13  ;   CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GEN_MASS ', CR13  ;   CALL DEALLOCATE_EIGEN1_MAT ( 'GEN_MASS' )
      CALL DEALLOCATE_COL_VEC    ( 'YSe' )
      CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC' )

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

! Write LINK5 end to F06

      CALL OURTIM
      WRITE(F06,151) LINKNO

! Close files

      IF (( DEBUG(193) == 5) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK5' )
      ENDIF

! Write LINK5 end to screen
      WRITE(SC1,153) LINKNO

! **********************************************************************************************************************************
  101 format(32767(1es22.14))

  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END',/)

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  153 FORMAT(  ' >> LINK',I3,' END')

 5001 FORMAT(' *ERROR  5001: INVALID SOLUTION NAME: SOL = ',A)

 5002 FORMAT(' *ERROR  5002: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                     ,/,14X,'VARIABLE LOAD_ISTEP MUST BE 1 OR 2 BUT VALUE IS = ',I8)

 5094 FORMAT(/,' >> LINK',I2,' END',19X,I2,':',I2,':',I2,'.',I3,/)

 5101 FORMAT(' *WARNING    : THE LARGEST OFF-DIAGONAL GENERALIZED MASS TERM REPORTED IN THE EIGENVALUE ANALYSIS SUMMARY CANNOT BE' &
                    ,/,14X,' RESCALED IN LINK ',I2,' BASED ON THE EIGENVALUE RENORMALIZATION REQUESTED.'                           &
                    ,/,14X,' IT IS RENORMALIZED BASED ON UNIT GENERALIZED MASS',/)

 5102 FORMAT(' *INFORMATION: ATTEMPTING TO RENORMALIZE EIGENVECTORS BASED ON GRID POINT-COMPONENT ',2I8                            &
                    ,/,14X,' THE RENORMALIZATION WILL BE DONE FOR ALL EIGENVECTORS WHOSE AMPLITUDE IS NONZERO AT THIS DOF',/)

 5103 FORMAT(' *INFORMATION: EIGENVECTORS WILL BE RENORMALIZED BASED ON MAX VALUE = 1.0',/)

 5104 FORMAT(' *INFORMATION: EIGENVECTORS WERE NORMALIZED TO MASS IN LINK4',/)

 5111 FORMAT(' *WARNING    : REQUEST TO NORMALIZE EIGENVECTORS BASED ON GRID POINT ',I8,', HOWEVER, THIS NOT A VALID GRID POINT.'  &
                    ,/,14X,' EIGENVECTORS WILL NOT BE RENORMALIZED',/)

 9995 FORMAT(/,' PROCESSING ENDED IN LINK ',I3,' DUE TO ABOVE ',I8,' ERRORS')

 9998 FORMAT(' *ERROR  9998: COMM ',I3,' INDICATES UNSUCCESSFUL LINK ',I2,' COMPLETION.'                                           &
                    ,/,14X,' FATAL ERROR - CANNOT START LINK ',I2)

98001 FORMAT(41X,'R E A L   E I G E N V A L U E S')

98002 FORMAT(82X,'BEFORE RENORMALIZATION OF EIGENVECTORS')

98003 FORMAT(3X,' MODE  EXTRACTION     EIGENVALUE            RADIANS              CYCLES            GENERALIZED         GENERALIZED&
&        ',/,3X,'NUMBER   ORDER                                                                        MASS              STIFFNESS'&
          ,/)

98004 FORMAT(1X,2I8,5(1ES20.6))

99885 FORMAT(82X,'MATRIX PHIZG',/,82X,'------------')

99886 FORMAT(5X,32676(I14))

99887 FORMAT(I8,'-',I1,32767(1ES14.6))

99888 FORMAT(8X,'-',I1,32767(1ES14.6))

12345 FORMAT(A,10X,A)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE READ_EIGNORM2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  EIN, EINFIL, ERR, F06
      USE SCONTR, ONLY                :  IERRFL, JCARD_LEN

      IMPLICIT NONE

      LOGICAL                         :: FILE_EXIST

      CHARACTER(JCARD_LEN*BYTE)       :: DATA_FIELD
      CHARACTER(80*BYTE)              :: TITLE             ! First record in EINFIL

      INTEGER(LONG)                   :: IOCHK             ! Vaue of IOSTAT in file open
      INTEGER(LONG)                   :: NUM_CHANGED       ! Number of eigenvectors to have their sign changed
      INTEGER(LONG)                   :: REC_NUM           ! Number of the record read from file
      INTEGER(LONG)                   :: VEC_NUM           ! Number of a eigenvector read from EINFIL
      INTEGER(LONG)                   :: VECS_CHANGED(NVEC)! Numbers of eigenvectors to have their sign changed

! **********************************************************************************************************************************
      DO I=1,NDOFL
         VEC_SIGN_CHG(I) = .FALSE.
      ENDDO

      EIGNORM2_ERR = 0
      INQUIRE ( FILE=EINFIL, EXIST=FILE_EXIST )

      IF (FILE_EXIST) THEN

         EIGNORM2_ERR = 0

         OPEN (EIN, FILE=EINFIL, STATUS='OLD', IOSTAT=IOCHK)

         IF (IOCHK == 0) THEN                              ! File was opened successfully

            WRITE(ERR,*)
            WRITE(ERR,5996) EIGNORM2
            IF (SUPINFO == 'N') THEN
               WRITE(F06,*)
               WRITE(F06,5996) EIGNORM2
            ENDIF
            CALL WRITE_FILNAM ( EINFIL, F06, 15 )
            READ(EIN,'(A)') TITLE

            NUM_CHANGED = 0

            REC_NUM = 0
            DO

               READ(EIN,'(A)',IOSTAT=IOCHK) DATA_FIELD   ;   REC_NUM = REC_NUM + 1

               IF (IOCHK == 0) THEN

                  IF (DATA_FIELD(1:3) == 'END') EXIT

                  CALL I4FLD (DATA_FIELD(1:JCARD_LEN), 1, VEC_NUM )
                  IF (IERRFL(1) == 'N') THEN
                     IF ((VEC_NUM >= 1) .AND. (VEC_NUM <= NVEC)) THEN
                        VEC_SIGN_CHG(VEC_NUM) = .TRUE.
                        NUM_CHANGED = NUM_CHANGED + 1
                        VECS_CHANGED(NUM_CHANGED) = VEC_NUM
                     ELSE
                        WRITE(ERR,5992) REC_NUM, DATA_FIELD, NVEC
                        WRITE(F06,5992) REC_NUM, DATA_FIELD, NVEC
                     ENDIF
                  ELSE
                     WRITE(ERR,5995) REC_NUM, DATA_FIELD
                     WRITE(F06,5995) REC_NUM, DATA_FIELD
                     EXIT
                  ENDIF

                  CYCLE

               ELSE

                  WRITE(ERR,5991)
                  CALL WRITE_FILNAM ( EINFIL, ERR, 15 )
                  WRITE(ERR,59911)

                  WRITE(F06,5991)
                  CALL WRITE_FILNAM ( EINFIL, F06, 15 )
                  WRITE(F06,59911)

                  EIGNORM2_ERR = EIGNORM2_ERR + 1
                  EXIT

               ENDIF

               IF (REC_NUM > NDOFL) THEN
                  EIGNORM2_ERR = EIGNORM2_ERR + 1
                  WRITE(ERR,5998) EINFIL, REC_NUM, NVEC
                  WRITE(F06,5998) EINFIL, REC_NUM, NVEC
               ENDIF

            ENDDO

         ELSE                                              ! File not opened successfully, so no vecs will have signs changed

            WRITE(ERR,5997)
            CALL WRITE_FILNAM ( EINFIL, ERR, 15 )
            WRITE(ERR,59971)

            WRITE(F06,5997)
            CALL WRITE_FILNAM ( EINFIL, F06, 15 )
            WRITE(F06,59971)

            EIGNORM2_ERR = EIGNORM2_ERR + 1

         ENDIF

      ELSE

         WRITE(ERR,5990)
         CALL WRITE_FILNAM ( EINFIL, ERR, 15 )
         WRITE(ERR,59901)

         WRITE(F06,5990)
         CALL WRITE_FILNAM ( EINFIL, F06, 15 )
         WRITE(F06,59901)

         EIGNORM2_ERR = EIGNORM2_ERR + 1

      ENDIF

      IF (EIGNORM2_ERR == 0) THEN
         WRITE(ERR,*)
         WRITE(ERR,5993)
         WRITE(ERR,5994) (VECS_CHANGED(I),I=1,NUM_CHANGED)
         WRITE(ERR,*)
         IF (SUPINFO == 'N') THEN
            WRITE(F06,*)
            WRITE(F06,5993)
            WRITE(F06,5994) (VECS_CHANGED(I),I=1,NUM_CHANGED)
            WRITE(F06,*)
         ENDIF
      ELSE
         WRITE(F06,5999)
         CALL WRITE_FILNAM ( EINFIL, F06, 15 )
      ENDIF

      RETURN

! **********************************************************************************************************************************
 5990 FORMAT(' *WARNING    : PARAMETER EIGNORM2 REQUESTS THAT SOME EIGENVECTORS ARE TO HAVE THEIR SIGN CHANGED, BUT FILE:')

59901 FORMAT('               WHICH SHOULD HAVE THE NUMBERS OF THE VECTORS TO HAVE THEIR SIGN CHANGED CANNOT BE FOUND')

 5991 FORMAT(' *WARNING    : EOF/ERR READING FILE:')

59911 FORMAT('               ON RECORD NUMBER ',I8,'. FIELD THAT HAS ERROR = "',A,'"')

 5992 FORMAT(' *WARNING    : ON REC NUM ',I8,' EIGENVEC ',A,' FOR SIGN CHANGE WAS OUT OF RANGE.',                                  &
                           ' MUST BE 1 TO NVEC = ',I8,'. VALUE IGNORED')

 5993 FORMAT(' *INFORMATION: THE FOLLOWING EIGENVECTORS HAVE HAD THEIR SIGN CHANGED:')

 5994 FORMAT(10I8)

 5995 FORMAT(' *WARNING    : ERROR READING RECORD NUMBER ',I8,' = ',A,'. REMAINDER OF EIGENVECTOR SIGN CHANGE FILE IGNORED')

 5996 FORMAT(' *INFORMATION: BASED ON PARAMETER EIGNORM2 = ',A,' SOME EIGENVECS WILL HAVE THEIR SIGN CHANGED AS DEFINED IN FILE:')

 5997 FORMAT(' *WARNING    : PARAMETER EIGNORM2 REQUESTS THAT SOME EIGENVECTORS ARE TO HAVE THEIR SIGN CHANGED, BUT FILE:')

59971 FORMAT('               WHICH SHOULD HAVE THE NUMBERS OF THE VECTORS TO HAVE THEIR SIGN CHANGED COULD NOT BE OPENED')

 5998 FORMAT(' *WARNING    : ERROR READING FILE: ',A                                                                               &
                    ,/,14X,' ATTEMPTING TO READ RECORD NUMBER ',I8,' BUT THERE SHOULD ONLY BE NVEC = ',I8,' RECORDS')

 5999 FORMAT(' *WARNING    : NO EIGENVECTORS WILL HAVE THEIR SIGN CHANGED DUE TO ABOVE ERRORS READING FILE:')

! **********************************************************************************************************************************

      END SUBROUTINE READ_EIGNORM2

      END SUBROUTINE LINK5


      SUBROUTINE BUILD_A_LR ( COL_NUM )

! For one subcase:

!   1) Merge UL and UR to get UA where UL was read into subr LINK5

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NDOFL, NDOFA, NDOFR, NVEC, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  PRTDISP
      USE COL_VECS, ONLY              :  UL_COL, UA_COL, UR_COL

      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE MATRIX_MERGING, ONLY        :  MERGE_COL_VECS
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_VECTOR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'BUILD_A_LR'

      INTEGER(LONG), INTENT(IN)       :: COL_NUM
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: A_SET_COL         ! Col no. in TDOF for A  displ set definition
      INTEGER(LONG)                   :: L_SET_COL         ! Col no. in TDOF for L  displ set definition
      INTEGER(LONG)                   :: R_SET_COL         ! Col no. in TDOF for R  displ set definition




! **********************************************************************************************************************************
! Get column numbers for various DOF sets

      IF (NDOFR > 0) THEN
                                                           ! Merge UL and UR to get UA
         CALL TDOF_COL_NUM('A ', A_SET_COL)
         CALL TDOF_COL_NUM('L ', L_SET_COL)
         CALL TDOF_COL_NUM('R ', R_SET_COL)

         DO I=1,NDOFR
            UR_COL(I) = ZERO
         ENDDO

         IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (COL_NUM > 0)) THEN
            UR_COL(COL_NUM-(NDOFR+NVEC)) = ONE
         ENDIF

         CALL MERGE_COL_VECS ( L_SET_COL, NDOFL, UL_COL, R_SET_COL, NDOFR, UR_COL, A_SET_COL, NDOFA, UA_COL )

      ELSE
                                                           ! Set UA = UL if no R set DOF's
         DO I=1,NDOFA
            UA_COL(I) = UL_COL(I)
         ENDDO

      ENDIF

! Print out displ matrices if PRTDISP says to

      IF ((PRTDISP(3) == 1) .OR. (PRTDISP(3) == 3)) THEN
         IF (NDOFL  > 0) THEN
            CALL WRITE_VECTOR ( '   F-SET DISPL VECTOR   ', 'DISPL', NDOFL, UL_COL)
         ENDIF
      ENDIF

      IF ((PRTDISP(3) == 2) .OR. (PRTDISP(3) == 3)) THEN
         IF (NDOFR  > 0) THEN
            CALL WRITE_VECTOR ( '   S-SET DISPL VECTOR   ', 'DISPL', NDOFR, UR_COL)
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BUILD_A_LR


      SUBROUTINE BUILD_F_AO

! For one subcase:

!   1) Calcs UO displacements: UO = GOA*UA + UO0 where:

!          UA  = from calculation in subr BUILD_A_LR
!          UO0 = KOO(-1)*PO from LINK2

!   2) Merge UO and UA to get UF

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NDOFA, NDOFF, NDOFO, NTERM_GOA, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE PARAMS, ONLY                :  PRTDISP
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE SPARSE_MATRICES, ONLY       :  I_GOA, J_GOA, GOA, SYM_GOA
      USE COL_VECS, ONLY              :  UA_COL, UF_COL, UO_COL, UO0_COL

      USE SPARSE_FULL_MULTIPLICATION, ONLY:  MATMULT_SFF
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE MATRIX_MERGING, ONLY        :  MERGE_COL_VECS
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_VECTOR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'BUILD_F_AO'

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: F_SET_COL         ! Col no. in TDOF for F  displ set definition
      INTEGER(LONG)                   :: A_SET_COL         ! Col no. in TDOF for A  displ set definition
      INTEGER(LONG)                   :: O_SET_COL         ! Col no. in TDOF for O  displ set definition
      INTEGER(LONG), PARAMETER        :: NUMCOLS     = 1   ! Variable for number of cols of an array




! **********************************************************************************************************************************
! Multiply GOA x UA to recover part of UO

      IF (NDOFO > 0) THEN

         CALL MATMULT_SFF ( 'GOA', NDOFO, NDOFA, NTERM_GOA, SYM_GOA, I_GOA, J_GOA, GOA, 'UA', NDOFA, NUMCOLS, UA_COL, 'Y',         &
                            'UO', ONE, UO_COL )

! Add UO0 (= KOO(-1) x PO) to get final UO but only if statics solution

         IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                                  &
            ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
            DO I=1,NDOFO
               UO_COL(I) = UO_COL(I) + UO0_COL(I)
            ENDDO
         ENDIF

! Merge UA and UO to get UF

         CALL TDOF_COL_NUM ( 'F ', F_SET_COL )
         CALL TDOF_COL_NUM ( 'A ', A_SET_COL )
         CALL TDOF_COL_NUM ( 'O ', O_SET_COL )

         CALL MERGE_COL_VECS ( A_SET_COL, NDOFA, UA_COL, O_SET_COL, NDOFO, UO_COL, F_SET_COL, NDOFF, UF_COL )

      ELSE

         DO I=1,NDOFF
            UF_COL(I) = UA_COL(I)
         ENDDO

      ENDIF

! Print out displ matrices if PRTDISP says to

      IF ((PRTDISP(4) == 1) .OR. (PRTDISP(4) == 3)) THEN
         IF (NDOFA  > 0) THEN
            CALL WRITE_VECTOR ( '   A-SET DISPL VECTOR   ', 'DISPL', NDOFA, UA_COL)
         ENDIF
      ENDIF

      IF ((PRTDISP(4) == 2) .OR. (PRTDISP(4) == 3)) THEN
         IF (NDOFO  > 0) THEN
            CALL WRITE_VECTOR ( '   O-SET DISPL VECTOR   ', 'DISPL', NDOFO, UO_COL)
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BUILD_F_AO


      SUBROUTINE BUILD_N_FS

! For one subcase:

!   1) Merge UF and US to get UN where UF is calc'd in subr BUILD_F_AO

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  NDOFF, NDOFN, NDOFS, NDOFSE, NDOFSZ, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  PRTDISP
      USE COL_VECS, ONLY              :  UF_COL, UN_COL, US_COL, YSe

      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE MATRIX_MERGING, ONLY        :  MERGE_COL_VECS
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_VECTOR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'BUILD_N_FS'

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: N_SET_COL         ! Col no. in TDOF for N  displ set definition
      INTEGER(LONG)                   :: F_SET_COL         ! Col no. in TDOF for F  displ set definition
      INTEGER(LONG)                   :: S_SET_COL         ! Col no. in TDOF for S  displ set definition
      INTEGER(LONG)                   :: SZ_SET_COL        ! Col no. in TDOF for SZ displ set definition
      INTEGER(LONG)                   :: SE_SET_COL        ! Col no. in TDOF for SE displ set definition


      REAL(DOUBLE)                    :: USZ_COL(NDOFSZ)   ! Array of zero displs for the SZ set



! **********************************************************************************************************************************
! Get column numbers for various DOF sets

      IF (NDOFS > 0) THEN

         CALL TDOF_COL_NUM('N ', N_SET_COL)
         CALL TDOF_COL_NUM('F ', F_SET_COL)
         CALL TDOF_COL_NUM('S ', S_SET_COL)
         CALL TDOF_COL_NUM('SZ',SZ_SET_COL)
         CALL TDOF_COL_NUM('SE',SE_SET_COL)

! Merge zeros for USZ ( the SZ-set) with YSe to get US ( the S-set)

         DO I=1,NDOFSZ
            USZ_COL(I) = ZERO
         ENDDO

         IF (NDOFSE > 0) THEN
            CALL MERGE_COL_VECS ( SZ_SET_COL, NDOFSZ, USZ_COL, SE_SET_COL, NDOFSE, YSe, S_SET_COL, NDOFS, US_COL )
         ELSE
            DO I=1,NDOFS
               US_COL(I) = ZERO
            ENDDO
         ENDIF

! Merge UF and US to get UN

         CALL MERGE_COL_VECS ( F_SET_COL, NDOFF, UF_COL, S_SET_COL, NDOFS, US_COL, N_SET_COL, NDOFN, UN_COL )

      ELSE
         DO I=1,NDOFN
            UN_COL(I) = UF_COL(I)
         ENDDO

      ENDIF


! Print out displ matrices if PRTDISP says to

      IF ((PRTDISP(3) == 1) .OR. (PRTDISP(3) == 3)) THEN
         IF (NDOFF  > 0) THEN
            CALL WRITE_VECTOR ( '   F-SET DISPL VECTOR   ', 'DISPL', NDOFF, UF_COL)
         ENDIF
      ENDIF

      IF ((PRTDISP(3) == 2) .OR. (PRTDISP(3) == 3)) THEN
         IF (NDOFS  > 0) THEN
            CALL WRITE_VECTOR ( '   S-SET DISPL VECTOR   ', 'DISPL', NDOFS, US_COL)
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BUILD_N_FS


      SUBROUTINE BUILD_G_NM

! For one subcase:

!   1) Calcs UM displacements: UM = GMN*UN where:

!          UN  = Displs calc'd in subr BUILD_N_FS

!   2) Merge UM and UN to get UG

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  NDOFG, NDOFM, NDOFN, NTERM_GMN, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE PARAMS, ONLY                :  PRTDISP
      USE SPARSE_MATRICES, ONLY       :  I_GMN, J_GMN, GMN, SYM_GMN
      USE COL_VECS, ONLY              :  UG_COL, UM_COL, UN_COL

      USE SPARSE_FULL_MULTIPLICATION, ONLY:  MATMULT_SFF
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE MATRIX_MERGING, ONLY        :  MERGE_COL_VECS
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_VECTOR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'BUILD_G_NM'

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: G_SET_COL         ! Col no. in TDOF for G  displ set definition
      INTEGER(LONG)                   :: N_SET_COL         ! Col no. in TDOF for N  displ set definition
      INTEGER(LONG)                   :: M_SET_COL         ! Col no. in TDOF for M  displ set definition
      INTEGER(LONG), PARAMETER        :: NUMCOLS     = 1   ! Variable for number of cols of an array




! **********************************************************************************************************************************
! Recover UM from GMN x UN

      IF (NDOFM > 0) THEN

         CALL MATMULT_SFF ( 'GMN', NDOFM, NDOFN, NTERM_GMN, SYM_GMN, I_GMN, J_GMN, GMN, 'UN', NDOFN, NUMCOLS, UN_COL, 'Y',         &
                            'UM', ONE, UM_COL )
! Merge UN and UM to get UG

         CALL TDOF_COL_NUM ( 'G ', G_SET_COL )
         CALL TDOF_COL_NUM ( 'N ', N_SET_COL )
         CALL TDOF_COL_NUM ( 'M ', M_SET_COL )

         CALL MERGE_COL_VECS ( N_SET_COL, NDOFN, UN_COL, M_SET_COL, NDOFM, UM_COL,G_SET_COL, NDOFG, UG_COL )

      ELSE

         DO I=1,NDOFG
            UG_COL(I) = UN_COL(I)
         ENDDO

      ENDIF


! Print out displ matrices if PRTDISP says to

      IF (PRTDISP(1) == 1) THEN
         IF (NDOFG  > 0) THEN
            CALL WRITE_VECTOR ( '   G-SET DISPL VECTOR   ', 'DISPL', NDOFG, UG_COL)
         ENDIF
      ENDIF

      IF ((PRTDISP(2) == 1) .OR. (PRTDISP(2) == 3)) THEN
         IF (NDOFN  > 0) THEN
            CALL WRITE_VECTOR ( '   N-SET DISPL VECTOR   ', 'DISPL', NDOFN, UN_COL)
         ENDIF
      ENDIF

      IF ((PRTDISP(2) == 2) .OR. (PRTDISP(2) == 3)) THEN
         IF (NDOFM  > 0) THEN
            CALL WRITE_VECTOR ( '   M-SET DISPL VECTOR   ', 'DISPL', NDOFM, UM_COL)
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BUILD_G_NM


      SUBROUTINE EXPAND_PHIXA_TO_PHIXG

! Expand array PHIXA (whose cols are stored in array UA_COL) to G-set columns (UG_COL). Each UG_COL is a column of matrix PHIXG.
! The Craig-Bampton mode array PHIXA is described in the MYSTRAN User's Reference Manual, Appendix D

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ONE
      USE IOUNT1, ONLY                :  ERR, F06, L5B, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, LINKNO, NDOFA, NDOFF, NDOFG, NDOFM, NDOFN, NDOFO, NDOFR, NDOFS, NTERM_PHIXA,&
                                         NTERM_PHIXG, NVEC, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE COL_VECS, ONLY              :  UA_COL, UG_COL
      USE PARAMS, ONLY                :  EPSIL, TINY
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE SPARSE_MATRICES, ONLY       :  I_PHIXA, J_PHIXA, PHIXA, I_PHIXG, J_PHIXG, PHIXG
      USE COL_VEC_LIFECYCLE, ONLY     :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
      USE SPARSE_CRS_ACCESS, ONLY     :  GET_SPARSE_CRS_COL
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CNT_NONZ_IN_FULL_MAT, LINK_MESSAGE, LINK_MESSAGE_I
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  FULL_TO_SPARSE_CRS


      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EXPAND_PHIXA_TO_PHIXG'
      CHARACTER( 1*BYTE)              :: NULL_COL          ! = 'Y' if col of PHIXA is null

      INTEGER(LONG)                   :: I,J               ! DO loop indices


      REAL(DOUBLE)                    :: PHIXG_FULL(NDOFG,NDOFR+NVEC)
!                                                          ! Full representation of matrix PHIXG before converting to sparse matrix

      REAL(DOUBLE)                    :: SMALL             ! A number used in filtering out small numbers from a full matrix



! **********************************************************************************************************************************
! Expand PHIXA (cols stored in UA_COL) to G-set columns (UG_COL). Each UG_COL is a column of matrix PHIXG

      DO J=1,NDOFR+NVEC
                                                           ! Put a col of PHIXA (A-set) into UA_COL
         CALL ALLOCATE_COL_VEC ( 'UA_COL', NDOFA, SUBR_NAME )
         CALL GET_SPARSE_CRS_COL ( 'PHIXA', J, NTERM_PHIXA, NDOFA, NDOFR+NVEC, I_PHIXA, J_PHIXA, PHIXA, ONE, UA_COL, NULL_COL )
                                                           ! Build F-set from A and O-set
         CALL ALLOCATE_COL_VEC ( 'UF_COL' , NDOFF, SUBR_NAME )
         CALL ALLOCATE_COL_VEC ( 'UO_COL' , NDOFO, SUBR_NAME )
         CALL LINK_MESSAGE_I('BUILD UF DISPLS FROM UA, UO:                      "', J)
         CALL BUILD_F_AO
         CALL DEALLOCATE_COL_VEC ( 'UA_COL' )
         CALL DEALLOCATE_COL_VEC ( 'UO_COL' )
                                                           ! Build N-set from F and S-set
         CALL ALLOCATE_COL_VEC ( 'UN_COL', NDOFN, SUBR_NAME)
         CALL ALLOCATE_COL_VEC ( 'US_COL', NDOFS, SUBR_NAME )
         CALL LINK_MESSAGE_I('BUILD UN DISPLS FROM UF, US:                      "',J)
         CALL BUILD_N_FS
         CALL DEALLOCATE_COL_VEC ( 'UF_COL' )
         CALL DEALLOCATE_COL_VEC ( 'US_COL' )
                                                           ! Build G-set from N and M-set
         CALL ALLOCATE_COL_VEC ( 'UG_COL', NDOFG, SUBR_NAME )
         CALL ALLOCATE_COL_VEC ( 'UM_COL', NDOFM, SUBR_NAME )
         CALL LINK_MESSAGE_I('BUILD UG DISPLS FROM UN, UM:                      "', J)
         CALL BUILD_G_NM
         CALL DEALLOCATE_COL_VEC ( 'UN_COL' )
         CALL DEALLOCATE_COL_VEC ( 'UM_COL' )

                                                           ! Write UG displs for this subcase to file LINK5A
         CALL LINK_MESSAGE_I('WRITE PHIXG DISPLS TO FILE,                       "', J)
   !xx   WRITE(SC1, * )                                    ! Separator between UG_COL calcs
         DO I=1,NDOFG
            WRITE(L5B) UG_COL(I)                           ! For CB this is a col of PHIXG (which is never processed as an array)
         ENDDO

         DO I=1,NDOFG
            PHIXG_FULL(I,J) = UG_COL(I)
         ENDDO

         CALL DEALLOCATE_COL_VEC ( 'UG_COL' )

      ENDDO

! Convert full PHIXG to sparse

      CALL CNT_NONZ_IN_FULL_MAT ( 'PHIZG_FULL', PHIXG_FULL, NDOFG, NDOFR+NVEC, 'N', NTERM_PHIXG, SMALL )
      CALL ALLOCATE_SPARSE_MAT  ( 'PHIXG', NDOFG, NTERM_PHIXG, SUBR_NAME )
      CALL FULL_TO_SPARSE_CRS   ( 'PHIXG_FULL', NDOFG, NDOFR+NVEC, PHIXG_FULL, NTERM_PHIXG, SMALL, SUBR_NAME, 'N',                 &
                                   I_PHIXG, J_PHIXG, PHIXG )



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *INFORMATION: TERMS WHOSE ABS VALUE ARE < MACH_PREC =',1ES10.3,' ARE NOT INCLUDED IN MATRIX ',A,' IN SUBR ',A       &
                    ,/,14X,' AS THIS FULL MATRIX IS BEING CONVERTED TO A SPARSE MATRIX')

  102 FORMAT(' *INFORMATION: TERMS WHOSE ABS VALUE ARE < PARAM TINY =',1ES10.3,' ARE NOT INCLUDED IN MATRIX ',A,' IN SUBR ',A      &
                    ,/,14X,' AS THIS FULL MATRIX IS BEING CONVERTED TO A SPARSE MATRIX')

99885 FORMAT(82X,'MATRIX PHIXG',/,82X,'------------')

99886 FORMAT(5X,32676(I14))

99887 FORMAT(I8,'-',I1,32767(1ES14.6))

99888 FORMAT(8X,'-',I1,32767(1ES14.6))

! **********************************************************************************************************************************

      END SUBROUTINE EXPAND_PHIXA_TO_PHIXG




      SUBROUTINE RENORM ( VEC_NUM, NORM_GRD, NORM_COMP, NORM, NORM_GSET_DOF, GEN_MASS1, PHI_SCALE_FAC )

! Renormalizes eigenves based on NORM = POINT or MAX if requested on Bulk Data entry EIGR or EIGRL

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NDOFG, NDOFG, NGRID, WARN_ERR
      USE PARAMS, ONLY                :  EPSIL, SUPWARN
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE COL_VECS, ONLY              :  UG_COL

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RENORM'
      CHARACTER( 8*BYTE), INTENT(IN)  :: NORM              ! Eigenvector renormalization methof from EIGR card (e.g. 'MAX     ')

      INTEGER(LONG), INTENT(IN)       :: NORM_COMP         ! Comp. (1-6) for renormalizing eigenvectors (from EIGR card)
      INTEGER(LONG), INTENT(IN)       :: NORM_GRD          ! Grid Point  for renormalizing eigenvectors (from EIGR card)
      INTEGER(LONG), INTENT(IN)       :: NORM_GSET_DOF     ! G-set DOF no. for NORM_GRD/NORM_COMP
      INTEGER(LONG), INTENT(IN)       :: VEC_NUM           ! Number used to control an output message (only want this information
!                                                            message written if tyhis is the first call to this subr).
      INTEGER(LONG)                   :: I                 ! DO loop index


      REAL(DOUBLE) , INTENT(INOUT)    :: GEN_MASS1         ! Generalized mass for 1 eigenvector
      REAL(DOUBLE) , INTENT(OUT)      :: PHI_SCALE_FAC     ! Scale factor for the eigenvector to renormalize it
      REAL(DOUBLE)                    :: DPHI_MAX          ! Absolute value of PHI_MAX
      REAL(DOUBLE)                    :: EPS1              ! Small number to compare variables against zero
      REAL(DOUBLE)                    :: PHI_MAX           ! Largest DZIJ for all DOF'si for one eigenvector
      REAL(DOUBLE)                    :: PHI_POINT         ! Variable used when normalizing gen. mass and eigenvectors

      INTRINSIC DSQRT,DABS



! **********************************************************************************************************************************
! Initialize outputs

      PHI_SCALE_FAC = ONE

! Check for renorm = MAX or POINT and renormalize.

      EPS1 = EPSIL(1)

      IF (NORM == 'POINT   ') THEN                         ! Renormalize eigenvector on POINT

         IF (NORM_GSET_DOF > 0) THEN                       ! If not, msg was written in LINK5 for no renorm, so return

            PHI_POINT = UG_COL(NORM_GSET_DOF)
            IF (DABS(PHI_POINT) > EPS1) THEN
               DO I=1,NDOFG
                  UG_COL(I) = UG_COL(I)/PHI_POINT
               ENDDO
               PHI_SCALE_FAC = PHI_POINT
               GEN_MASS1 = GEN_MASS1/(PHI_SCALE_FAC * PHI_SCALE_FAC)
            ELSE
               PHI_SCALE_FAC = ONE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,4113) VEC_NUM,NORM_GRD,NORM_COMP
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,4113) VEC_NUM,NORM_GRD,NORM_COMP
               ENDIF
            ENDIF

         ELSE

            RETURN                                         ! No renorm if NORM_GSET_DOF undefined (msg written in LINK5)

         ENDIF

      ELSE                                                 ! NORM is MAX

         PHI_MAX  = ZERO
         DPHI_MAX = ZERO
         DO I=1,NDOFG                                      ! Scan eigenvector to find largest value (+ or -)
            IF (DABS(UG_COL(I)) > DPHI_MAX) THEN
               PHI_MAX  = UG_COL(I)
               DPHI_MAX = DABS(PHI_MAX)
            ENDIF
         ENDDO

         IF (DPHI_MAX > EPS1) THEN                         ! Renormalize the eigenvector if PHI_MAX is not 0
            DO I=1,NDOFG
               UG_COL(I) = UG_COL(I)/PHI_MAX
            ENDDO
            PHI_SCALE_FAC = PHI_MAX
            GEN_MASS1 = GEN_MASS1/(PHI_SCALE_FAC * PHI_SCALE_FAC)
         ELSE
            PHI_SCALE_FAC = ONE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,4114) VEC_NUM
            IF (SUPWARN == 'N') THEN
               WRITE(F06,4114) VEC_NUM
            ENDIF
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 4113 FORMAT(' *WARNING    : EIGENVECTOR ',I8,' IS ZERO FOR GRID POINT-COMPONENT ',2I8,'. THIS VECTOR WILL NOT BE RENORMALIZED')

 4114 FORMAT(' *WARNING    : EIGENVECTOR ',I8,' MAX VALUE IS ZERO. THIS VECTOR CANNOT BE RENORMALIZED')

! **********************************************************************************************************************************

      END SUBROUTINE RENORM

   END MODULE LINK5_MOD
