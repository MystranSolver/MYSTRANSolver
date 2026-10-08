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

   MODULE REDUCTION_G_TO_N

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: REDUCE_G_NM

   CONTAINS

      SUBROUTINE REDUCE_G_NM

! Call routines to reduce stiffness, mass, loads and constraint matrices from G-set to N, M-sets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1C, LINK1C, L1C_MSG, SC1, WRT_ERR

      USE SCONTR, ONLY                :  LINKNO    , NDOFG, NDOFN, NDOFM, NGRID, NSUB,                                             &
                                         NTERM_KGG , NTERM_KNN , NTERM_KNM , NTERM_KMM ,                                           &
                                         NTERM_KGGD, NTERM_KNND, NTERM_KNMD, NTERM_KMMD,                                           &
                                         NTERM_MGG , NTERM_MNN , NTERM_MNM , NTERM_MMM , NTERM_PG, NTERM_PN, NTERM_PM ,            &
                                         NTERM_GMN , NTERM_RMN , NTERM_RMM ,                                                       &
                                         PROG_NAME , SOL_NAME, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE PARAMS, ONLY                :  AUTOSPC, AUTOSPC_NSET, EQCHK_OUTPUT, MATSPARS, PRTSTIFD, PRTSTIFF, PRTMASS, PRTFOR,       &
                                         SUPINFO
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE DOF_TABLES, ONLY            :  TDOF, TDOFI
      USE MODEL_STUF, ONLY            :  GRID_ID
      USE RIGID_BODY_DISP_MATS, ONLY  :  RBGLOBAL_GSET, RBGLOBAL_NSET
      USE SPARSE_MATRICES, ONLY       :  I_KGG , J_KGG , KGG , I_KGGD, J_KGGD, KGGD,                                               &
                                         I_KNN , J_KNN , KNN , I_KNM , J_KNM , KNM , I_KMM , J_KMM , KMM ,                         &
                                         I_KNND, J_KNND, KNND, I_KNMD, J_KNMD, KNMD, I_KMMD, J_KMMD, KMMD,                         &
                                         I_MGG , J_MGG , MGG , I_MNN , J_MNN , MNN , I_MNM , J_MNM , MNM , I_MMM , J_MMM , MMM ,   &
                                         I_PG  , J_PG  , PG  , I_PN  , J_PN  , PN  , I_PM  , J_PM  , PM  ,                         &
                                         I_RMG , J_RMG , RMG

      USE SPARSE_MATRICES, ONLY       :  SYM_KNN
      USE OUTPUT4_MATRICES, ONLY      :  ACT_OU4_MYSTRAN_NAMES, NUM_OU4_REQUESTS
      USE SCRATCH_MATRICES

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_VEC
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  AUTOSPC_SUMMARY_MSGS, GET_MATRIX_DIAG_STATS
      USE RIGID_BODY_STORAGE_LIFECYCLE, ONLY:  ALLOCATE_RBGLOBAL
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM, TDOF_PROC
      USE REDUCTION_CHECKS, ONLY      :  STIFF_MAT_EQUIL_CHK
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  FILERR, FILE_CLOSE, FILE_OPEN
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_DOF_TABLES
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_G_NM'
      CHARACTER(  8*BYTE)             :: ASPC_SUM_MSG1       ! Message to be printed out in the AUTOSPC summary table
      CHARACTER(100*BYTE)             :: ASPC_SUM_MSG2       ! Message to be printed out in the AUTOSPC summary table
      CHARACTER( 13*BYTE)             :: ASPC_SUM_MSG3       ! Message to be printed out in the AUTOSPC summary table
      CHARACTER(  1*BYTE)             :: DEALLOCATE_KGG = 'Y'! Indicator of whether we need to keep KGG  allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_KGGD= 'Y'! Indicator of whether we need to keep KGGD allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_MGG = 'Y'! Indicator of whether we need to keep MGG allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_PG  = 'Y'! Indicator of whether we need to keep PG  allocated for OU4 output
      CHARACTER(132*BYTE)             :: MATRIX_NAME         ! Name of matrix for printout
      CHARACTER(44*BYTE)              :: MODNAM              ! Name to write to screen to describe module being run

      INTEGER(LONG)                   :: DO_WHICH_CODE_FRAG    ! 1 or 2 depending on which seg of code to run (depends on BUCKLING)
      INTEGER(LONG)                   :: I,J,K               ! DO loop indices
      INTEGER(LONG)                   :: N_SET_COL           ! Col no. in array TDOFI where the N-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: N_SET_DOF           ! N-set DOF number
      INTEGER(LONG)                   :: NUM_ASPC_BY_COMP(6) ! Number of AUTOSPC's by component number
      INTEGER(LONG)                   :: NUM_COMPS           ! 6 if GRID_NUM is an physical grid, 1 if an SPOINT
      INTEGER(LONG)                   :: PART_VEC_G_NM(NDOFG)! Partitioning vector (G set into N and M sets)
      INTEGER(LONG)                   :: PART_VEC_M(NDOFM)   ! Partitioning vector (1's for all M set DOF's)
      INTEGER(LONG)                   :: PART_VEC_SUB(NSUB)  ! Partitioning vector (1's for all subcases)
      INTEGER(LONG)                   :: SA_SET_COL          ! Col no. in array TDOF where the SA-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: TOT_NUM_ASPC        ! Sum of NUM_ASPC_BY_COMP(6)

      REAL(DOUBLE)                    :: KNN_DIAG(NDOFN)     ! Diagonal terms from KNN
      REAL(DOUBLE)                    :: KNN_MAX_DIAG        ! Max diag term from  KNN
      REAL(DOUBLE)                    :: KNND_DIAG(NDOFN)    ! Diagonal terms from KNND
      REAL(DOUBLE)                    :: KNND_MAX_DIAG       ! Max diag term from  KNND

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Determine if we need to keep any OUTPUT4 matrices allocated until after they are processed in LINK2

      IF (NUM_OU4_REQUESTS > 0) THEN
         DO I=1,NUM_OU4_REQUESTS
            IF      (ACT_OU4_MYSTRAN_NAMES(I)(1:3) == 'KGG') THEN
               DEALLOCATE_KGG = 'N'
            ELSE IF (ACT_OU4_MYSTRAN_NAMES(I)(1:3) == 'MGG') THEN
               DEALLOCATE_MGG = 'N'
            ELSE IF (ACT_OU4_MYSTRAN_NAMES(I)(1:3) == 'PG' ) THEN
               DEALLOCATE_PG  = 'N'
            ENDIF
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Depending on whether this is a BUCKLING soln (and LOAD_ISTEP value) or not, one or another segment of code will be run

      IF ((SOL_NAME(1:8) == 'BUCKLING')) THEN
         IF      (LOAD_ISTEP == 1) THEN
            DO_WHICH_CODE_FRAG = 1
         ELSE IF (LOAD_ISTEP == 2) THEN
            DO_WHICH_CODE_FRAG = 2
         ENDIF
      ELSE
         DO_WHICH_CODE_FRAG = 1
      ENDIF

! **********************************************************************************************************************************
      IF (DO_WHICH_CODE_FRAG == 1) THEN                    ! This is for all except BUCKLING w LOAD_ISTEP=2 (eigen part of BUCKLING)

! If there is an M-set, reduce KGG to KNN, MGG to MNN, PG to PN using UM = GMN*UN, where GMN = -RMM(-1)*RMN (partitions of RMG)
! If there is no M-set, then equate KNN to KGG, MNN to MGG, PN to PG

! First, need to create partitioning vectors used in the reduction (if NDOFM > 0)

         IF (NDOFM > 0) THEN

            CALL PARTITION_VEC (NDOFG,'G ','N ','M ',PART_VEC_G_NM)

            DO I=1,NDOFM
               PART_VEC_M(I) = 1
            ENDDO

            DO I=1,NSUB
               PART_VEC_SUB = 1
            ENDDO

            CALL OURTIM
            MODNAM = '  SOLVE FOR GMN CONSTRAINT MATRIX'
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            CALL SOLVE_GMN ( PART_VEC_G_NM, PART_VEC_M )   ! First, solve for GMN
      !xx   WRITE(SC1,  * )                                ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RMG', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'RMG' )

            IF (NTERM_KGG > 0) THEN                        ! Reduce KGG to KNN

               CALL OURTIM
               IF (MATSPARS == 'Y') THEN
                  MODNAM = '  REDUCE KGG TO KNN (SPARSE MATRIX ROUTINES)'
               ELSE
                  MODNAM = '  REDUCE KGG TO KNN (FULL MATRIX ROUTINES)'
               ENDIF
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_KGG_TO_KNN ( PART_VEC_G_NM )

            ELSE

               NTERM_KNN = 0
               NTERM_KNM = 0
               NTERM_KMM = 0
               CALL ALLOCATE_SPARSE_MAT ( 'KNN', NDOFN, NTERM_KNN, SUBR_NAME )

            ENDIF

            IF (NTERM_MGG > 0) THEN                        ! Reduce MGG to MNN

               CALL OURTIM
               IF (MATSPARS == 'Y') THEN
                  MODNAM = '  REDUCE MGG TO MNN (SPARSE MATRIX ROUTINES)'
               ELSE
                  MODNAM = '  REDUCE MGG TO MNN (FULL MATRIX ROUTINES)'
               ENDIF
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_MGG_TO_MNN ( PART_VEC_G_NM )

            ELSE

               NTERM_MNN = 0
               NTERM_MNM = 0
               NTERM_MMM = 0
               CALL ALLOCATE_SPARSE_MAT ( 'MNN', NDOFN, NTERM_MNN, SUBR_NAME )

            ENDIF

            IF ((SOL_NAME(1:5) /= 'MODES') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')) THEN

               IF (NTERM_PG > 0) THEN                      ! Reduce PG to PN

                  CALL OURTIM
                  IF (MATSPARS == 'Y') THEN
                     MODNAM = '  REDUCE PG  TO PN  (SPARSE MATRIX ROUTINES)'
                  ELSE
                     MODNAM = '  REDUCE PG  TO PN  (FULL MATRIX ROUTINES)'
                  ENDIF
                  WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

                  CALL REDUCE_PG_TO_PN ( PART_VEC_G_NM, PART_VEC_SUB )

               ELSE

                  NTERM_PN = 0
                  NTERM_PM = 0
                  CALL ALLOCATE_SPARSE_MAT ( 'PN', NDOFN, NTERM_PN, SUBR_NAME )

               ENDIF

            ENDIF

         ELSE                                              ! There is no M-set, so equate N and G sets

            CALL OURTIM
            MODNAM = '  EQUATING N-SET TO G-SET'
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            NDOFN     = NDOFG

            NTERM_KNN = NTERM_KGG
            NTERM_KNM = 0
            NTERM_KMM = 0

            NTERM_MNN = NTERM_MGG
            NTERM_MNM = 0
            NTERM_MMM = 0

            NTERM_PN  = NTERM_PG
            NTERM_PM  = 0

            NTERM_GMN = 0

            CALL ALLOCATE_SPARSE_MAT ( 'KNN', NDOFN, NTERM_KNN, SUBR_NAME )

!xx         DO I=1,NDOFN+1
!xx            I_KNN(I) = I_KGG(I)
!xx         ENDDO
!xx
!xx         DO I=1,NTERM_KNN
!xx            J_KNN(I) = J_KGG(I)
!xx              KNN(I) =   KGG(I)
!xx         ENDDO
!xx
            CALL ALLOCATE_SPARSE_MAT ( 'MNN', NDOFN, NTERM_MNN, SUBR_NAME )

!xx         DO I=1,NDOFN+1
!xx            I_MNN(I) = I_MGG(I)
!xx         ENDDO
!xx
!xx         DO I=1,NTERM_MNN
!xx            J_MNN(I) = J_MGG(I)
!xx              MNN(I) =   MGG(I)
!xx         ENDDO
!xx
            IF ((SOL_NAME(1:5) /= 'MODES') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')) THEN

               CALL ALLOCATE_SPARSE_MAT ( 'PN', NDOFN, NTERM_PN, SUBR_NAME )

!xx            DO I=1,NDOFN+1
!xx               I_PN(I) = I_PG(I)
!xx            ENDDO
!xx
!xx            DO I=1,NTERM_PN
!xx            J_PN(I) = J_PG(I)
!xx              PN(I) =   PG(I)
!xx            ENDDO
!xx
            ENDIF

         ENDIF

! Deallocate G-set arrays

         MODNAM = '  DEALLOCATE G-SET ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages

         IF (DEALLOCATE_KGG == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KGG', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KGG' )
         ENDIF

         IF (DEALLOCATE_MGG == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MGG', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MGG' )
         ENDIF

         IF (DEALLOCATE_PG == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PG ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PG' )
         ENDIF
         WRITE(SC1,12345,ADVANCE='NO')    '       Deallocate GMN', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'GMN' )

! Print out stiffness matrix partitions, if requested

         IF (( PRTSTIFF(2) == 1) .OR. ( PRTSTIFF(2) == 3)) THEN
            IF (NTERM_KNN > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KNN'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'N ', 'N ', NTERM_KNN, NDOFN, I_KNN, J_KNN, KNN )
            ENDIF
         ENDIF

         IF (( PRTSTIFF(2) == 2) .OR. ( PRTSTIFF(2) == 3)) THEN
            IF (NTERM_KNM > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KNM'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'N ', 'M ', NTERM_KNM, NDOFN, I_KNM, J_KNM, KNM )
            ENDIF
            IF (NTERM_KMM > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KMM'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'M ', 'M ', NTERM_KMM, NDOFM, I_KMM, J_KMM, KMM )
            ENDIF
         ENDIF

         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KNM', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KNM' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KMM', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KMM' )

! Write matrix diagonal and stats, if requested.
! NOTE: call this subr even if PRTSTIFFD(2) = 0 since we need KNN_DIAG, KNN_MAX_DIAG for the equilibrium check

         CALL GET_MATRIX_DIAG_STATS ( 'KNN', 'N ', NDOFN, NTERM_KNN, I_KNN, J_KNN, KNN, PRTSTIFD(2), KNN_DIAG, KNN_MAX_DIAG )

! Print out mass matrix partitions, if requested

         IF (( PRTMASS(2) == 1) .OR. ( PRTMASS(2) == 3)) THEN
            IF (NTERM_MNN > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MNN'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'N ', 'N ', NTERM_MNN, NDOFN, I_MNN, J_MNN, MNN )
            ENDIF
         ENDIF

         IF (( PRTMASS(2) == 2) .OR. ( PRTMASS(2) == 3)) THEN
            IF (NTERM_MNM > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MNM'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'N ', 'M ', NTERM_MNM, NDOFN, I_MNM, J_MNM, MNM )
            ENDIF
            IF (NTERM_MMM > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MMM'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'M ', 'M ', NTERM_MMM, NDOFM, I_MMM, J_MMM, MMM )
            ENDIF
         ENDIF

         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MMN', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MMN' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MNM', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MNM' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MMM', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MMM' )

! Print out load matrix partitions, if requested

         IF (( PRTFOR(2) == 1) .OR. ( PRTFOR(2) == 3)) THEN
            IF (NTERM_PN  > 0) THEN
               MATRIX_NAME = 'LOAD MATRIX PN'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'N ', 'SUBCASE', NTERM_PN, NDOFN, I_PN, J_PN, PN )
            ENDIF
         ENDIF

         IF (( PRTFOR(2) == 2) .OR. ( PRTFOR(2) == 3)) THEN
            IF (NTERM_PM  > 0) THEN
               MATRIX_NAME = 'LOAD MATRIX PM'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'M ', 'SUBCASE', NTERM_PM, NDOFM, I_PM, J_PM, PM )
            ENDIF
         ENDIF

         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PM ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PM' )
         WRITE(SC1,*) CR13

! Do equilibrium check on the N-set stiffness matrix, if requested

         IF ((EQCHK_OUTPUT(2) > 0) .OR. (EQCHK_OUTPUT(3) > 0) .OR. (EQCHK_OUTPUT(4) > 0) .OR. (EQCHK_OUTPUT(5) > 0)) THEN
            CALL ALLOCATE_RBGLOBAL ( 'N ', SUBR_NAME )
            IF (NDOFM > 0) THEN
               CALL TDOF_COL_NUM ( 'N ', N_SET_COL )
               DO I=1,NDOFG
                  N_SET_DOF = TDOFI(I,N_SET_COL)
                  IF (N_SET_DOF > 0) THEN
                     DO J=1,6
                        RBGLOBAL_NSET(N_SET_DOF,J) = RBGLOBAL_GSET(I,J)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               DO I=1,NDOFN
                  DO J=1,6
                     RBGLOBAL_NSET(I,J) = RBGLOBAL_GSET(I,J)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF

         IF (EQCHK_OUTPUT(2) > 0) THEN
            CALL OURTIM
            MODNAM = '  EQUILIBRIUM CHECK ON KNN                '
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
            CALL STIFF_MAT_EQUIL_CHK ( EQCHK_OUTPUT(2),'N ', SYM_KNN, NDOFN, NTERM_KNN, I_KNN, J_KNN, KNN, KNN_DIAG, KNN_MAX_DIAG,&
                                       RBGLOBAL_NSET)
         ENDIF

! If AUTOSPC = 'Y' check to see if any rows of KNN are null for DOF's that are not already in the S or O sets
! (in KGG_SINGULARIYT_PROC we passed by grids that were indep on MPC and some of these may have zero stiffness)
! Also, if user requested, check for small terms on the diag of KNN and AUTOSPC if their ratio to max diag term is small enough

         IF (AUTOSPC == 'Y') THEN
            IF ((AUTOSPC_NSET == 1) .OR. (AUTOSPC_NSET == 3)) THEN ! Check for null rows
            CALL OURTIM
            MODNAM = '  N SET AUTOSPC PROC #1                   '
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
            CALL N_SET_AUTOSPC_PROC_1
            ENDIF
            IF ((AUTOSPC_NSET == 2) .OR. (AUTOSPC_NSET == 3)) THEN ! Check for small diag terms
            CALL OURTIM
            MODNAM = '  N SET AUTOSPC PROC #2                   '
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
               CALL N_SET_AUTOSPC_PROC_2
            ENDIF
         ENDIF

! Now print final AUTOSPC summary table. Need to calc NUM_ASPC_BY_COMP from TDOF table

         CALL TDOF_COL_NUM ( 'SA', SA_SET_COL )

         DO J=1,6
            NUM_ASPC_BY_COMP(J) = 0
         ENDDO

         CALL OURTIM
         MODNAM = '  AUTOSPC SUMMARY TABLE                   '
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
         K = 0
         DO I=1,NGRID
            CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               K = K + 1
               IF (TDOF(K,SA_SET_COL) /= 0) THEN
                  NUM_ASPC_BY_COMP(J) = NUM_ASPC_BY_COMP(J) + 1
               ENDIF
            ENDDO
         ENDDO

         TOT_NUM_ASPC = 0
         DO J=1,6
            TOT_NUM_ASPC = TOT_NUM_ASPC + NUM_ASPC_BY_COMP(J)
         ENDDO

         ASPC_SUM_MSG1(1:) = 'Overall:'
         ASPC_SUM_MSG2(1:) = 'after identification of all AUTOSPC''s'
         ASPC_SUM_MSG3(1:) = 'overall      '
         CALL AUTOSPC_SUMMARY_MSGS ( ASPC_SUM_MSG1, ASPC_SUM_MSG2, ASPC_SUM_MSG3, 'Y', NUM_ASPC_BY_COMP )

! **********************************************************************************************************************************
      ELSE                                                 ! This is BUCKLING with LOAD_ISTEP = 2 (eigen part of BUCKLING)

         IF (NDOFM > 0) THEN

            CALL PARTITION_VEC (NDOFG,'G ','N ','M ',PART_VEC_G_NM)

            DO I=1,NDOFM
               PART_VEC_M(I) = 1
            ENDDO

            DO I=1,NSUB
               PART_VEC_SUB = 1
            ENDDO

! Reduce KGG to KNN

            IF (NTERM_KGGD > 0) THEN                          ! Reduce KGGD to KNND

               CALL OURTIM
               MODNAM = '  REDUCE KGGD TO KNND (SPARSE MATRIX ROUTINES)'
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_KGGD_TO_KNND ( PART_VEC_G_NM )

            ELSE

               NTERM_KNND = 0
               NTERM_KNMD = 0
               NTERM_KMMD = 0
               CALL ALLOCATE_SPARSE_MAT ( 'KNND', NDOFN, NTERM_KNND, SUBR_NAME )

            ENDIF

            CALL DEALLOCATE_SPARSE_MAT ( 'RMG' )

! There is no M-set, so equate N and G sets

         ELSE

            CALL OURTIM
            MODNAM = '  EQUATING N-SET TO G-SET'
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            NDOFN      = NDOFG

            NTERM_KNND = NTERM_KGGD
            NTERM_KNMD = 0
            NTERM_KMMD = 0
            CALL ALLOCATE_SPARSE_MAT ( 'KNND', NDOFN, NTERM_KNND, SUBR_NAME )

         ENDIF

! Deallocate G-set arrays

         MODNAM = '  DEALLOCATE G-SET ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages

         IF (DEALLOCATE_KGGD == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KGGD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KGGD' )
         ENDIF

! Print out stiffness matrix partitions, if requested

         IF (( PRTSTIFF(2) == 1) .OR. ( PRTSTIFF(2) == 3)) THEN
            IF (NTERM_KNND > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KNND'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'N ', 'N ', NTERM_KNND, NDOFN, I_KNND, J_KNND, KNND )
            ENDIF
         ENDIF

         IF (( PRTSTIFF(2) == 2) .OR. ( PRTSTIFF(2) == 3)) THEN
            IF (NTERM_KNMD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KNMD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'N ', 'M ', NTERM_KNMD, NDOFN, I_KNMD, J_KNMD, KNMD )
            ENDIF
            IF (NTERM_KMMD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KMMD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'M ', 'M ', NTERM_KMMD, NDOFM, I_KMMD, J_KMMD, KMMD )
            ENDIF
         ENDIF

         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KNMD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KNMD' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KMMD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KMMD' )

! Write matrix diagonal and stats, if requested.
! NOTE: call this subr even if PRTSTIFFD(2) = 0 since we need KNN_DIAG, KNN_MAX_DIAG for the equilibrium check

         CALL GET_MATRIX_DIAG_STATS ( 'KNND', 'N ', NDOFN, NTERM_KNND, I_KNND, J_KNND, KNND, PRTSTIFD(2), KNND_DIAG, KNND_MAX_DIAG )

      ENDIF



      RETURN

! **********************************************************************************************************************************
 2092 FORMAT(4X,A44,20X,I2,':',I2,':',I2,'.',I3)

12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE N_SET_AUTOSPC_PROC_1

! Checks KNN to see if any rows are null for DOF's not already in the S or O-sets, and, if so, puts these in the SA set and
! reruns subr TDOF_PROC and writes the new TSET, TDOF, TDOFI tables to file L1C

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  DATA_NAM_LEN, FATAL_ERR, NDOFG, NDOFSA, NGRID, NUM_PCHD_SPC1
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1C, L1C_MSG, LINK1C, SPC, SPCFIL
      USE PARAMS, ONLY                :  AUTOSPC, AUTOSPC_INFO, AUTOSPC_NSET, PCHSPC1, PRTTSET, SPC1SID
      USE DOF_TABLES, ONLY            :  TDOF, TDOFI, TSET
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID, GRID_SEQ

      IMPLICIT NONE

      CHARACTER(  7*BYTE)             :: ASPC_SUM_MSG1      ! Message to be printed out in the AUTOSPC summary table
      CHARACTER(100*BYTE)             :: ASPC_SUM_MSG2      ! Message to be printed out in the AUTOSPC summary table
      CHARACTER( 13*BYTE)             :: ASPC_SUM_MSG3      ! Message to be printed out in the AUTOSPC summary table
      CHARACTER(132*BYTE)             :: TDOF_MSG           ! Msg to be printed out regarding at what point in the run the TDOF,I
!                                                             tables are printed out

      INTEGER(LONG)                   :: AGRID              ! Actual grid ID for IGRID
      INTEGER(LONG)                   :: COMP               ! DOF component number (1-6)
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM    ! Row number in array GRID_ID where AGRID is found
      INTEGER(LONG)                   :: I,J                ! DO loop indices
      INTEGER(LONG)                   :: IOCHK              ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: NUM_ASPC_BY_COMP(6)! Number of AUTOSPC's by component number
      INTEGER(LONG)                   :: NUM_N_SET_ROWS_NULL! Number of rows in KNN that are null and are not S or O-set members
      INTEGER(LONG)                   :: N_SET_COL          ! Col no. in array TDOF where the  N-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG), ALLOCATABLE      :: N_SET_TDOFI_ROW(:) ! Row in TDOFI for each N-set DOF number
      INTEGER(LONG)                   :: R_SET_COL          ! Col no. in array TDOF where the  R-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: S_SET_COL          ! Col no. in array TDOF where the  S-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: OUNT(2)            ! File units to write messages to. Input to subr UNFORMATTED_OPEN
! **********************************************************************************************************************************
      OUNT(1) = ERR
      OUNT(2) = F06

! Open file SPC to write SPC1 records if this subr finds singularities and user wants SPCFIL written

      IF (NUM_PCHD_SPC1 > 0) THEN                          ! Subr KGG_SINGULARITY_PROC already opened and wrote to this file
         OPEN (SPC, FILE=SPCFIL, STATUS='OLD', POSITION='APPEND', IOSTAT=IOCHK)
      ELSE                                                 ! File has not been written to, so open as replace
         OPEN (SPC, FILE=SPCFIL, STATUS='REPLACE', IOSTAT=IOCHK)
      ENDIF
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SPCFIL, OUNT )
         CALL FILERR ( OUNT )
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL TDOF_COL_NUM ( 'N ',  N_SET_COL )
      CALL TDOF_COL_NUM ( 'R ',  R_SET_COL )
      CALL TDOF_COL_NUM ( 'S ',  S_SET_COL )

      IF (NDOFN > 0) THEN
         ALLOCATE ( N_SET_TDOFI_ROW(NDOFN), STAT=IOCHK )
         IF (IOCHK /= 0) THEN
            WRITE(ERR,*) ' *ERROR: ALLOCATING N_SET_TDOFI_ROW IN ', SUBR_NAME
            WRITE(F06,*) ' *ERROR: ALLOCATING N_SET_TDOFI_ROW IN ', SUBR_NAME
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         N_SET_TDOFI_ROW = 0
         DO J=1,NDOFG
            N_SET_DOF = TDOFI(J,N_SET_COL)
            IF (N_SET_DOF > 0) N_SET_TDOFI_ROW(N_SET_DOF) = J
         ENDDO
      ENDIF

      WRITE(ERR,101) AUTOSPC_NSET, PROG_NAME
      IF (SUPINFO == 'N') THEN
         WRITE(F06,101) AUTOSPC_NSET, PROG_NAME
      ENDIF

! Look for null rows in KNN. If found, move that DOF to the SA set

      DO I=1,6                                             ! Initialize NUM_ASPC_BY_COMP
         NUM_ASPC_BY_COMP(I) = 0
      ENDDO

      NUM_N_SET_ROWS_NULL = 0
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      CALL COUNTER_INIT('       Proc N-set DOF ', NDOFN)
i_do: DO I=1,NDOFN
         IF (I_KNN(I+1) == I_KNN(I)) THEN                  ! If true, row i is null
            J = N_SET_TDOFI_ROW(I)
            IF (J > 0) THEN
               IF (TDOFI(J,N_SET_COL) == I) THEN
                  IF ((TDOFI(J,S_SET_COL) == 0) .AND. (TDOFI(J,R_SET_COL) == 0)) THEN
                     NUM_N_SET_ROWS_NULL = NUM_N_SET_ROWS_NULL + 1
                     AGRID = TDOFI(J,1)
                     COMP  = TDOFI(J,2)
                     NUM_ASPC_BY_COMP(COMP) = NUM_ASPC_BY_COMP(COMP) + 1
                     CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, GRID_ID_ROW_NUM )
                     TSET ( GRID_ID_ROW_NUM, COMP ) = 'SA'
                     NDOFSA = NDOFSA + 1
                     IF (AUTOSPC_INFO == 'Y') THEN
                        WRITE(ERR,102) AGRID, COMP, AUTOSPC
                        IF (SUPINFO == 'N') THEN
                           WRITE(F06,102) AGRID, COMP, AUTOSPC
                        ENDIF
                     ENDIF
                     IF (PCHSPC1 == 'Y') THEN
                        WRITE(SPC,109) SPC1SID, COMP, AGRID
                        NUM_PCHD_SPC1 = NUM_PCHD_SPC1 + 1
                     ENDIF
                  ENDIF
               ELSE
                  WRITE(ERR,*) ' *ERROR: N_SET_AUTOSPC_PROC_1 LOOKUP MISMATCH FOR N-SET DOF ', I
                  WRITE(F06,*) ' *ERROR: N_SET_AUTOSPC_PROC_1 LOOKUP MISMATCH FOR N-SET DOF ', I
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ELSE
               WRITE(ERR,*) ' *ERROR: N_SET_AUTOSPC_PROC_1 LOOKUP FAILED FOR N-SET DOF ', I
               WRITE(F06,*) ' *ERROR: N_SET_AUTOSPC_PROC_1 LOOKUP FAILED FOR N-SET DOF ', I
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
         ENDIF
         CALL COUNTER_PROGRESS(I)
      ENDDO i_do
      WRITE(SC1,*) CR13

! Close SPC file

      IF (NUM_PCHD_SPC1 > 0) THEN
         CALL FILE_CLOSE ( SPC, SPCFIL,  'KEEP' )
      ELSE
         CALL FILE_CLOSE ( SPC, SPCFIL,  'DELETE' )
      ENDIF

! IF we changed some DOF's from the N-set to the SA-set regenerate TDOF, TDOFI tables and write them to L1C

      WRITE(F06,*)
      IF  (NUM_N_SET_ROWS_NULL > 0) THEN

         WRITE(ERR,103) PROG_NAME, NUM_N_SET_ROWS_NULL
         IF (SUPINFO == 'N') THEN
            WRITE(F06,103) PROG_NAME, NUM_N_SET_ROWS_NULL
         ENDIF

         IF (PRTTSET > 0) THEN
            WRITE(F06,56)
            WRITE(F06,57)
            DO J = 1,NGRID
               WRITE(F06,58) GRID(J,1), GRID_SEQ(J), (TSET(J,K),K = 1,6)
            ENDDO
            WRITE(F06,'(//)')
         ENDIF

         ASPC_SUM_MSG1(1:) = 'Stage 2:'
         ASPC_SUM_MSG2(1:) = 'after identification of AUTOSPC''s to eliminate null rows in the N-set stiffness matrix'
         ASPC_SUM_MSG3(1:) = 'in this stage'
         CALL AUTOSPC_SUMMARY_MSGS ( ASPC_SUM_MSG1, ASPC_SUM_MSG2, ASPC_SUM_MSG3, 'N', NUM_ASPC_BY_COMP )

         TDOF_MSG(1:)  = ' '
         TDOF_MSG(22:) = ASPC_SUM_MSG2(1:)
         CALL TDOF_PROC ( TDOF_MSG )

         OUNT(1) = ERR
         OUNT(2) = F06
         CALL FILE_OPEN ( L1C, LINK1C, OUNT, 'REPLACE', L1C_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL WRITE_DOF_TABLES
         CALL FILE_CLOSE ( L1C, LINK1C, 'KEEP' )

      ELSE

         WRITE(ERR,104) PROG_NAME
         IF (SUPINFO == 'N') THEN
            WRITE(F06,104) PROG_NAME
         ENDIF

      ENDIF
! **********************************************************************************************************************************
   56 FORMAT(64X,'DEGREE OF FREEDOM SET TABLE (TSET)')

   57 FORMAT(33x,'     GRID SEQUENCE       T1       T2       T3       R1       R2       R3',/)

   58 FORMAT(33x,2(1X,I8),6(7X,A2))

  101 FORMAT(' *INFORMATION: BASED ON PARAMETER AUTOSPC_NSET = ',I2,1X,A,' IS CHECKING KNN TO SEE IF THERE ARE NULL ROWS',         &
                           ' THAT SHOULD BE AUTOSPC''d',/)

  102 FORMAT(' *INFORMATION: GRID POINT ',I8,' HAS SINGULARITY FOR DISPL COMPONENT    ',5X,I1,'. SINCE PARAM AUTOSPC = ',A,        &
                          ', THIS  WILL BE AUTOSPC''d')

  103 FORMAT(' *INFORMATION: ',A,' WILL AUTOSPC ',I8,' DOF''s FROM THE N-SET THAT WERE PREVIOUSLY NOT MEMBERS OF THE S-SET',/)

  104 FORMAT(' *INFORMATION: ',A,' FOUND NO N-SET DOF''s THAT WERE SINGULAR AND THAT WERE NOT ALREADY MEMBERS OF THE S-SET',/)

  109 FORMAT('SPC1    ',3I8)

! **********************************************************************************************************************************

      END SUBROUTINE N_SET_AUTOSPC_PROC_1

! ##################################################################################################################################

      SUBROUTINE N_SET_AUTOSPC_PROC_2

! Checks KNN to see if any diag terms are small (compared to AUTOSPC_RAT) for DOF's not already in the S or O-sets, and, if so, puts
! these in the SA set and reruns subr TDOF_PROC and writes the new TSET, TDOF, TDOFI tables to file L1C

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  DATA_NAM_LEN, NDOFN, NDOFG, NDOFSA, NGRID, NUM_PCHD_SPC1, PROG_NAME
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L1C, L1C_MSG, LINK1C, SPC, SPCFIL
      USE PARAMS, ONLY                :  AUTOSPC, AUTOSPC_INFO, AUTOSPC_NSET, AUTOSPC_RAT, PCHSPC1, PRTTSET, SPC1SID
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DOF_TABLES, ONLY            :  TDOF, TDOFI, TSET
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID, GRID_SEQ
      USE SPARSE_MATRICES, ONLY       :  I_KFF, KFF

      IMPLICIT NONE

      CHARACTER(  7*BYTE)             :: ASPC_SUM_MSG1      ! Message to be printed out in the AUTOSPC summary table
      CHARACTER(100*BYTE)             :: ASPC_SUM_MSG2      ! Message to be printed out in the AUTOSPC summary table
      CHARACTER( 13*BYTE)             :: ASPC_SUM_MSG3      ! Message to be printed out in the AUTOSPC summary table
      CHARACTER(132*BYTE)             :: TDOF_MSG           ! Msg to be printed out regarding at what point in the run the TDOF,I
!                                                             tables are printed out

      INTEGER(LONG)                   :: AGRID              ! Actual grid ID for IGRID
      INTEGER(LONG)                   :: COMP               ! DOF component number (1-6)
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM    ! Row number in array GRID_ID where AGRID is found
      INTEGER(LONG)                   :: I,J                ! DO loop indices
      INTEGER(LONG)                   :: IOCHK              ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: JSTART             ! DO loop start point
      INTEGER(LONG)                   :: NUM_ASPC_BY_COMP(6)! Number of AUTOSPC's by component number
      INTEGER(LONG)                   :: NUM_NSET_DOFS_SPCD ! Number of rows in KNN that are null and are not S or O-set members
      INTEGER(LONG)                   :: N_SET_COL          ! Col no. in array TDOF where the  N-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: R_SET_COL          ! Col no. in array TDOF where the  R-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: S_SET_COL          ! Col no. in array TDOF where the  S-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: OUNT(2)            ! File units to write messages to. Input to subr UNFORMATTED_OPEN
! **********************************************************************************************************************************
      OUNT(1) = ERR
      OUNT(2) = F06

! Open file SPC to write SPC1 records if this subr finds singularities

      IF (NUM_PCHD_SPC1 > 0) THEN                          ! Subr KGG_SINGULARITY_PROC already opened and wrote to this file
         OPEN (SPC, FILE=SPCFIL, STATUS='OLD', POSITION='APPEND', IOSTAT=IOCHK)
      ELSE                                                 ! File has not been written to, so open as replace
         OPEN (SPC, FILE=SPCFIL, STATUS='REPLACE', IOSTAT=IOCHK)
      ENDIF
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SPCFIL, OUNT )
         CALL FILERR ( OUNT )
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL TDOF_COL_NUM ( 'N ',  N_SET_COL )
      CALL TDOF_COL_NUM ( 'R ',  R_SET_COL )
      CALL TDOF_COL_NUM ( 'S ',  S_SET_COL )

      WRITE(ERR,101) AUTOSPC_NSET, PROG_NAME, AUTOSPC_RAT
      IF (SUPINFO == 'N') THEN
         WRITE(F06,101) AUTOSPC_NSET, PROG_NAME, AUTOSPC_RAT
      ENDIF

! Check ratios of diag to max diag term. If smaller than requirement, move that DOF to the SA set

      DO I=1,6                                             ! Initialize NUM_ASPC_BY_COMP
         NUM_ASPC_BY_COMP(I) = 0
      ENDDO

      NUM_NSET_DOFS_SPCD = 0
      JSTART = 1
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      CALL COUNTER_INIT('       Proc N-set DOF ', NDOFN)
i_do: DO I=1,NDOFN
         IF ((DABS(KNN_DIAG(I)/KNN_MAX_DIAG) < AUTOSPC_RAT) .OR. (KNN_DIAG(I) < ZERO)) THEN
j_do:       DO J=JSTART,NDOFG                               ! Loop over rows of TDOFI to find where this N-set row is null
               IF (TDOFI(J,N_SET_COL) == I) THEN
                  IF ((TDOFI(J,S_SET_COL) == 0) .AND. (TDOFI(J,R_SET_COL) == 0)) THEN
                     NUM_NSET_DOFS_SPCD = NUM_NSET_DOFS_SPCD + 1
                     AGRID = TDOFI(J,1)
                     COMP  = TDOFI(J,2)
                     NUM_ASPC_BY_COMP(COMP) = NUM_ASPC_BY_COMP(COMP) + 1
                     CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, GRID_ID_ROW_NUM )
                     TSET ( GRID_ID_ROW_NUM, COMP ) = 'SA'
                     NDOFSA = NDOFSA + 1
                     IF (AUTOSPC_INFO == 'Y') THEN
                        WRITE(ERR,102) AGRID, COMP, AUTOSPC
                        IF (SUPINFO == 'N') THEN
                           WRITE(F06,102) AGRID, COMP, AUTOSPC
                        ENDIF
                     ENDIF
                     IF (PCHSPC1 == 'Y') THEN
                        WRITE(SPC,109) SPC1SID, COMP, AGRID
                        NUM_PCHD_SPC1 = NUM_PCHD_SPC1 + 1
                     ENDIF
                     JSTART = J
                     EXIT j_do
                  ENDIF
               ENDIF
            ENDDO j_do
         ENDIF
         CALL COUNTER_PROGRESS(I)
      ENDDO i_do
      WRITE(SC1,*) CR13

! Close SPC file

      IF (NUM_PCHD_SPC1 > 0) THEN
         CALL FILE_CLOSE ( SPC, SPCFIL,  'KEEP' )
      ELSE
         CALL FILE_CLOSE ( SPC, SPCFIL,  'DELETE' )
      ENDIF

! IF we changed some DOF's from the N-set to the SA-set regenerate TDOF, TDOFI tables and write them to L1C

      WRITE(F06,*)
      IF  (NUM_NSET_DOFS_SPCD > 0) THEN

         WRITE(ERR,103) PROG_NAME, NUM_NSET_DOFS_SPCD
         IF (SUPINFO == 'N') THEN
            WRITE(F06,103) PROG_NAME, NUM_NSET_DOFS_SPCD
         ENDIF

         IF (PRTTSET > 0) THEN
            WRITE(F06,56)
            WRITE(F06,57)
            DO J = 1,NGRID
               WRITE(F06,58) GRID(J,1), GRID_SEQ(J), (TSET(J,K),K = 1,6)
            ENDDO
            WRITE(F06,'(//)')
         ENDIF

         ASPC_SUM_MSG1(1:) = 'Stage 3:'
         ASPC_SUM_MSG2(1:) ='after identification of AUTOSPC''s to eliminate N-set DOFs with stiffness ratios < PARAM AUTOSPC_RAT'
         ASPC_SUM_MSG3(1:) = 'in this stage'
         CALL AUTOSPC_SUMMARY_MSGS ( ASPC_SUM_MSG1, ASPC_SUM_MSG2, ASPC_SUM_MSG3, 'Y', NUM_ASPC_BY_COMP )

         TDOF_MSG(1:)  = ' '
         TDOF_MSG(16:) = ASPC_SUM_MSG2(1:)
         CALL TDOF_PROC ( TDOF_MSG )

         OUNT(1) = ERR
         OUNT(2) = F06
         CALL FILE_OPEN ( L1C, LINK1C, OUNT, 'REPLACE', L1C_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

         CALL WRITE_DOF_TABLES

         CALL FILE_CLOSE ( L1C, LINK1C, 'KEEP' )

      ELSE

         WRITE(ERR,104) PROG_NAME, AUTOSPC_RAT
         IF (SUPINFO == 'N') THEN
            WRITE(F06,104) PROG_NAME, AUTOSPC_RAT
         ENDIF

      ENDIF
! **********************************************************************************************************************************
  56 FORMAT(64X,'DEGREE OF FREEDOM SET TABLE (TSET)')

   57 FORMAT(33x,'     GRID SEQUENCE       T1       T2       T3       R1       R2       R3',/)

   58 FORMAT(33x,2(1X,I8),6(7X,A2))

  101 FORMAT(' *INFORMATION: BASED ON PARAMETER AUTOSPC_NSET = ',I2,1X,A,' IS CHECKING KNN TO SEE IF THERE ARE DOF''s THAT ARE'    &
                          ,' NOT ALREADY IN THE'                                                                                   &
                    ,/,14X,' S-SET BUT SHOULD BE AUTOSPC''d BASED ON SMALL DIAGONAL TERMS WHOSE RATIO WITH MAX DIAGONAL TERM IS < '&
                    ,1ES13.6,/)

  102 FORMAT(' *INFORMATION: GRID POINT ',I8,' HAS SINGULARITY FOR DISPL COMPONENT    ',5X,I1,'. SINCE PARAM AUTOSPC = ',A,        &
                          ', THIS  WILL BE AUTOSPC''d')

  103 FORMAT(' *INFORMATION: ',A,' HAS AUTOSPC''d ',I8,' DOF''s FROM THE N-SET THAT WERE PREVIOUSLY NOT MEMBERS OF THE S-SET',/)

  104 FORMAT(' *INFORMATION: ',A,' FOUND NO N-SET DOF''s THAT HAD SMALL DIAG TERMS (RATIO TO MAX DIAG TERM < ',1ES15.6,')'         &
                    ,/,14X,' AND THAT WERE NOT ALREADY MEMBERS OF THE S-SET',/)

  105 FORMAT('               AUTOSPC_RAT = ',1ES13.6)

  109 FORMAT('SPC1    ',3I8)

! **********************************************************************************************************************************

      END SUBROUTINE N_SET_AUTOSPC_PROC_2

      END SUBROUTINE REDUCE_G_NM


      SUBROUTINE REDUCE_KGG_TO_KNN ( PART_VEC_G_NM )

! Call routines to reduce the KGG linear stiffness matrix from the G-set to the N, M-sets. See Appendix B to the MYSTRAN User's
! Reference Manual for the derivation of the reduction equations.

! NOTE: This subr has code for sparse matrices as well as full matrices, (i.e. Bulk Data PARAM MATSPARS = 'Y' for sparse and 'N'
! for full). The code for full matrices was put in originally in order that the sparse code could be thoroughly checked. That task
! is complete and the remaining full matrix code has not been maintained since around 2005. In addition, new capability added to
! MYSTRAN since that approx time does not have full matrix code.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L2J, LINK2J, L2J_MSG, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NDOFN, NDOFM, NTERM_HMN, NTERM_KGG, NTERM_KNN,            &
                                         NTERM_KNM, NTERM_KMM, NTERM_GMN
      USE PARAMS, ONLY                :  EPSIL, MATSPARS, SPARSTOR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE SPARSE_MATRICES, ONLY       :  I_HMN, J_HMN, HMN, I_KGG, J_KGG, KGG, I_KNN, J_KNN, KNN, I_KNM, J_KNM, KNM,               &
                                         I_KMM, J_KMM, KMM,I_KMN, J_KMN, KMN, I_GMN, J_GMN, GMN,  I_GMNt, J_GMNt, GMNt
      USE SPARSE_MATRICES, ONLY       :  SYM_GMN, SYM_HMN, SYM_KGG, SYM_KNN, SYM_KNM, SYM_KMM, SYM_KMN
      USE FULL_MATRICES, ONLY         :  KNN_FULL, KNM_FULL, KMM_FULL, GMN_FULL, DUM1, DUM2, DUM3
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM, MATTRNSP_SS
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_NONSYM_TO_CRS_SYM, SPARSE_CRS_SPARSE_CCS, SPARSE_CRS_TO_FULL
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_CRS_TERM_COUNT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_MATRIX_1
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATADD_FFF, MATMULT_FFF, MATMULT_FFF_T
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CNT_NONZ_IN_FULL_MAT

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_KGG_TO_KNN'
      CHARACTER(  1*BYTE)             :: SYM_CRS1            ! Storage format for matrix CRS1 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)
      CHARACTER(  1*BYTE)             :: SYM_CRS3            ! Storage format for matrix CRS3 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_G_NM(NDOFG)! Partitioning vector (G set into N and M sets)
      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG)                   :: ITRNSPB             ! Transpose indicator for matrix multiply routine
      INTEGER(LONG)                   :: KNN_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KNM_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
!xx   INTEGER(LONG)                   :: KMN_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KMM_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: NTERM_CCS1          ! Number of terms in matrix CCS1
      INTEGER(LONG)                   :: NTERM_CRS1          ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2          ! Number of terms in matrix CRS2
      INTEGER(LONG)                   :: NTERM_CRS3          ! Number of terms in matrix CRS3
      INTEGER(LONG)                   :: NTERM_KMN           ! Number of nonzeros in sparse matrix KMN (should = NTERM_KNM)
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices


      REAL(DOUBLE)                    :: ALPHA = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: BETA  = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: SMALL             ! A number used in filtering out small numbers from a full matrix

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Partition KNN from KGG (This is KNN before reduction, or KNN(bar) )

      IF (NDOFN > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KGG', NTERM_KGG, NDOFG, NDOFG, SYM_KGG, I_KGG, J_KGG,      PART_VEC_G_NM, PART_VEC_G_NM,       &
                                    NUM1, NUM1, KNN_ROW_MAX_TERMS, 'KNN', NTERM_KNN, SYM_KNN )

         CALL ALLOCATE_SPARSE_MAT ( 'KNN', NDOFN, NTERM_KNN, SUBR_NAME )

         IF (NTERM_KNN > 0) THEN
            CALL PARTITION_SS ( 'KGG', NTERM_KGG, NDOFG, NDOFG, SYM_KGG, I_KGG, J_KGG, KGG, PART_VEC_G_NM, PART_VEC_G_NM,          &
                                 NUM1, NUM1, KNN_ROW_MAX_TERMS, 'KNN', NTERM_KNN, NDOFN, SYM_KNN, I_KNN, J_KNN, KNN )
         ENDIF

      ENDIF

! Partition KNM from KGG

      IF ((NDOFN > 0) .AND. (NDOFM > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'KGG', NTERM_KGG, NDOFG, NDOFG, SYM_KGG, I_KGG, J_KGG,      PART_VEC_G_NM, PART_VEC_G_NM,       &
                                    NUM1, NUM2, KNM_ROW_MAX_TERMS, 'KNM', NTERM_KNM, SYM_KNM )

         CALL ALLOCATE_SPARSE_MAT ( 'KNM', NDOFN, NTERM_KNM, SUBR_NAME )

         IF (NTERM_KNM > 0) THEN
            CALL PARTITION_SS ( 'KGG', NTERM_KGG, NDOFG, NDOFG, SYM_KGG, I_KGG, J_KGG, KGG, PART_VEC_G_NM, PART_VEC_G_NM,          &
                                 NUM1, NUM2, KNM_ROW_MAX_TERMS, 'KNM', NTERM_KNM, NDOFN, SYM_KNM, I_KNM, J_KNM, KNM )
         ENDIF

      ENDIF

! Partition KMN from KGG

      IF ((NDOFN > 0) .AND. (NDOFM > 0)) THEN

!xx      CALL PARTITION_SS_NTERM ( 'KGG', NTERM_KGG, NDOFG, NDOFG, SYM_KGG, I_KGG, J_KGG,      PART_VEC_G_NM, PART_VEC_G_NM,       &
!xx                                 NUM2, NUM1, KMN_ROW_MAX_TERMS, 'KMN', NTERM_KMN, SYM_KMN )
!xx      IF (NTERM_KMN /= NTERM_KNM) THEN
!xx         FATAL_ERR = FATAL_ERR + 1
!xx         WRITE(ERR,936) SUBR_NAME, NTERM_KMN, NTERM_KNM
!xx         WRITE(F06,936) SUBR_NAME, NTERM_KMN, NTERM_KNM
!xx         CALL OUTA_HERE ( 'Y' )
!xx      ENDIF
!xx
!xx      CALL ALLOCATE_SPARSE_MAT ( 'KMN', NDOFM, NTERM_KMN, SUBR_NAME )
!xx
!xx      IF (NTERM_KMN > 0) THEN
!xx         CALL PARTITION_SS ( 'KGG', NTERM_KGG, NDOFG, NDOFG, SYM_KGG, I_KGG, J_KGG, KGG, PART_VEC_G_NM, PART_VEC_G_NM,          &
!xx                              NUM2, NUM1, KMN_ROW_MAX_TERMS, 'KMN', NTERM_KMN, NDOFM, SYM_KMN, I_KMN, J_KMN, KMN )
!xx      ENDIF
!xx
         NTERM_KMN = NTERM_KNM
         CALL ALLOCATE_SPARSE_MAT ( 'KMN', NDOFM, NTERM_KNM, SUBR_NAME )

         IF (NTERM_KMN > 0) THEN
            CALL MATTRNSP_SS ( NDOFN, NDOFM, NTERM_KNM, 'KNM', I_KNM, J_KNM, KNM, 'KMN', I_KMN, J_KMN, KMN )

         ENDIF

      ENDIF

! Partition KMM from KGG

      IF (NDOFM > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KGG', NTERM_KGG, NDOFG, NDOFG, SYM_KGG, I_KGG, J_KGG,      PART_VEC_G_NM, PART_VEC_G_NM,       &
                                    NUM2, NUM2, KMM_ROW_MAX_TERMS, 'KMM', NTERM_KMM, SYM_KMM )

         CALL ALLOCATE_SPARSE_MAT ( 'KMM', NDOFM, NTERM_KMM, SUBR_NAME )

         IF (NTERM_KMM > 0) THEN
            CALL PARTITION_SS ( 'KGG', NTERM_KGG, NDOFG, NDOFG, SYM_KGG, I_KGG, J_KGG, KGG, PART_VEC_G_NM, PART_VEC_G_NM,          &
                                 NUM2, NUM2, KMM_ROW_MAX_TERMS, 'KMM', NTERM_KMM, NDOFM, SYM_KMM, I_KMM, J_KMM, KMM )
         ENDIF

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Reduce KGG to KNN = KNN(bar) + KNM*GMN + (KNM*GMN)' + GMN'*KMM*GMN.
! If PARAM MATSPARS = 'Y', then we use sparse matrix operations (multiply/add/transpose). If not, use full matrix operations

      IF (MATSPARS == 'Y') THEN                              ! Reduce KGG to KNN using sparse matrix operations

         CALL ALLOCATE_SPARSE_MAT ( 'GMNt', NDOFN, NTERM_GMN, SUBR_NAME )
         CALL MATTRNSP_SS ( NDOFM, NDOFN, NTERM_GMN, 'GMN', I_GMN, J_GMN, GMN, 'GMNt', I_GMNt, J_GMNt, GMNt )
                                                           ! CCS1 will be sparse CCS format version of sparse CRS matrix GMN
         CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFN, NTERM_GMN, SUBR_NAME )
         CALL SPARSE_CRS_SPARSE_CCS ( NDOFM, NDOFN, NTERM_GMN, 'GMN', I_GMN, J_GMN, GMN, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

         IF (NTERM_KNM > 0) THEN                           ! Part I of reduced KNN: calc KNM*GMN & add it & transpose to orig KNN

                                                           ! I-1, sparse multiply to get CRS1 = KNM*GMN. Use CCS1 for GMN CCS
            CALL MATMULT_SSS_NTERM ( 'KNM' , NDOFN, NTERM_KNM , SYM_KNM, I_KNM , J_KNM ,                                           &
                                     'GMN' , NDOFN, NTERM_GMN , SYM_GMN, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1',        NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'KNM' , NDOFN, NTERM_KNM , SYM_KNM, I_KNM , J_KNM , KNM ,                                           &
                               'GMN' , NDOFN, NTERM_GMN , SYM_GMN, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )
            NTERM_CRS2 = NTERM_CRS1                        ! I-2, allocate memory to array CRS2 which will hold transpose of CRS1
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFN, NTERM_CRS2, SUBR_NAME )

                                                           ! I-3, transpose CRS1 to get CRS2 = (KNM*GMN)t
            CALL MATTRNSP_SS ( NDOFN, NDOFN, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CRS2', I_CRS2, J_CRS2, CRS2 )
                                                           ! I-4, sparse add to get CRS3 = CRS1 + CRS2 = (KNM*GMN) + (KNM*GMN)t
            CALL MATADD_SSS_NTERM (NDOFN,'KNM*GMN', NTERM_CRS1, I_CRS1, J_CRS1, 'N', '(KNM*GMN)t', NTERM_CRS2, I_CRS2, J_CRS2, 'N',&
                                         'CRS3', NTERM_CRS3)
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
            CALL MATADD_SSS ( NDOFN, 'KNM*GMN', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE, '(KNM*GMN)t', NTERM_CRS2,                   &
                              I_CRS2, J_CRS2, CRS2, ONE, 'CRS3', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! I-5, deallocate CRS1 which was KNM*GMN
            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! I-6, deallocate CRS2 which was (KNM*GMN)t

                                                           ! I-7, CRS3 = (KNM*GMN) + (KNM*GMN)t has all nonzero terms in it.
            IF      (SPARSTOR == 'SYM   ') THEN            !      If SPARSTOR == 'SYM   ', rewrite CRS3 as sym in CRS1

               CALL SPARSE_CRS_TERM_COUNT ( NDOFN, NTERM_CRS3, '(KNM*GMN) + (KNM*GMN)t', I_CRS3, J_CRS3, NTERM_CRS1 )
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )
               CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS3 = (KNM*GMN) + (KNM*GMN)t all nonzeros', NDOFN, NTERM_CRS3, I_CRS3, J_CRS3, CRS3, &
                                            'CRS1 = (KNM*GMN) + (KNM*GMN)t stored sym'  ,        NTERM_CRS1, I_CRS1, J_CRS1, CRS1 )
               SYM_CRS1 = 'Y'

            ELSE IF (SPARSTOR == 'NONSYM') THEN            !      If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

               NTERM_CRS1 = NTERM_CRS3
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )
               DO I=1,NDOFN+1
                  I_CRS1(I) = I_CRS3(I)
               ENDDO
               DO I=1,NTERM_CRS1
                  J_CRS1(I) = J_CRS3(I)
                    CRS1(I) =   CRS3(I)
               ENDDO
               SYM_CRS1 = 'N'

            ELSE                                           !      Error - incorrect SPARSTOR

               WRITE(ERR,932) SUBR_NAME, SPARSTOR
               WRITE(F06,932) SUBR_NAME, SPARSTOR
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF

            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! I-8, deallocate CRS3 which was KNM*GMN + (KNM*GMN)t and now is in CRS1

                                                           ! I-9, add to get CRS3 = KNN-bar + CRS1 = KNN-bar + KNM*GMN + (KNM*GMN)t
            CALL MATADD_SSS_NTERM ( NDOFN, 'KNN-bar', NTERM_KNN , I_KNN , J_KNN , SYM_KNN , 'KNM*GMN + (KNM*GMN)t',                &
                                                      NTERM_CRS1, I_CRS1, J_CRS1, SYM_CRS1, 'CRS3', NTERM_CRS3 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
            CALL MATADD_SSS ( NDOFN, 'KNN-bar', NTERM_KNN, I_KNN, J_KNN, KNN, ONE, 'KNM*GMN + (KNM*GMN)t', NTERM_CRS1,             &
                                      I_CRS1, J_CRS1, CRS1, ONE, 'CRS1', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )
            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! I-10, deallocate CRS1 = KNM*GMN + (KNM*GMN)t

            NTERM_KNN = NTERM_CRS3                         ! I-11, reallocate KNN to be size of CRS3
            WRITE(SC1, * ) '    Reallocate KNN'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KNN', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KNN' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KNN', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KNN', NDOFN, NTERM_KNN, SUBR_NAME )

            DO I=1,NDOFN+1                                 ! I-12, set KNN = CRS3
               I_KNN(I) = I_CRS3(I)
            ENDDO
            DO J=1,NTERM_KNN
               J_KNN(J) = J_CRS3(J)
                 KNN(J) =   CRS3(J)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! I-13, deallocate CRS3
                                                           ! At this point, CRS1, CRS2, CRS3 are deallocated, CCS1 is being used
         ENDIF

         IF (NTERM_KMM > 0) THEN                           ! Part II of reduced KNN: calc GMN(t)*KMM*GMN and add to KNN

                                                           ! II-1, sparse multiply to get CRS1 = KMM*GMN using CCS1 for GMN CCS
            CALL MATMULT_SSS_NTERM ( 'KMM' , NDOFM, NTERM_KMM , SYM_KMM, I_KMM , J_KMM ,                                           &
                                     'GMN' , NDOFN, NTERM_GMN , SYM_GMN, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1',        NTERM_CRS1 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFM, NTERM_CRS1, SUBR_NAME )
            CALL MATMULT_SSS ( 'KMM' , NDOFM, NTERM_KMM , SYM_KMM, I_KMM , J_KMM , KMM ,                                           &
                               'GMN' , NDOFN, NTERM_GMN , SYM_GMN, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )
            IF (NTERM_KMN > 0) THEN
               CALL MATADD_SSS_NTERM (NDOFM, 'KMM*GMN', NTERM_CRS1, I_CRS1, J_CRS1, 'N', 'KMN', NTERM_KMN, I_KMN, J_KMN, SYM_KMN,  &
                                             'HMN', NTERM_HMN)

               CALL ALLOCATE_SPARSE_MAT ( 'HMN', NDOFM, NTERM_HMN, SUBR_NAME )

               CALL MATADD_SSS ( NDOFM, 'KMM*GM', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE, 'KMN', NTERM_KMN, I_KMN, J_KMN, KMN,      &
                                 ONE, 'HMN', NTERM_HMN, I_HMN, J_HMN, HMN )
            ELSE
               NTERM_HMN = NTERM_CRS1
               CALL ALLOCATE_SPARSE_MAT ( 'HMN', NDOFM, NTERM_HMN, SUBR_NAME )
               DO I=1,NDOFM+1
                  I_HMN(I) = I_CRS1(I)
               ENDDO
               DO I=1,NTERM_HMN
                  J_HMN(I) = J_CRS1(I)
                  HMN(I) =   CRS1(I)
               ENDDO
            ENDIF

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )             ! II-2, deallocate CCS1 which was CCS version of GMN

            NTERM_CCS1 = NTERM_CRS1                        ! II-3, allocate CCS1 to be same as CRS1 = KMM*GMN but in CCS format
            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFN, NTERM_CCS1, SUBR_NAME )
            CALL SPARSE_CRS_SPARSE_CCS ( NDOFM, NDOFN, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! II-4, deallocate CRS1 which was KMM*GMN
                                                           ! II-5, sparse multiply to get CRS1 = GMNt*CCS1  with CCS1 = KMM*GMN
!                                                           (note: use SYM_GMN for sym indicator of KMM*GMN)
            CALL MATMULT_SSS_NTERM ( 'GMNt', NDOFN, NTERM_GMN , SYM_GMN, I_GMNt, J_GMNt,                                           &
                                     'CCS1', NDOFN, NTERM_CCS1, SYM_GMN, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1',        NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'GMNt', NDOFN, NTERM_GMN , SYM_GMN, I_GMNt, J_GMNt, GMNt,                                           &
                               'CCS1', NDOFN, NTERM_CCS1, SYM_GMN, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )
            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )             ! II-6, deallocate CCS1

                                                           ! II-7, CRS1 = GMNt*KMM*GMN has all nonzero terms in it.
            IF      (SPARSTOR == 'SYM   ') THEN            !      If SPARSTOR == 'SYM   ', rewrite CRS1 as sym in CRS3

               CALL SPARSE_CRS_TERM_COUNT ( NDOFN, NTERM_CRS1, 'GMNt*KMM*GMN all nonzeros', I_CRS1, J_CRS1, NTERM_CRS3 )
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
               CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS1 = GMNt*KMM*GMN all nonzeros', NDOFN, NTERM_CRS1, I_CRS1, J_CRS1, CRS1,           &
                                            'CRS3 = GMNt*KMM*GMN stored sym'  ,        NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )
               SYM_CRS3 = 'Y'

            ELSE IF (SPARSTOR == 'NONSYM') THEN            !      If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

               NTERM_CRS3 = NTERM_CRS1
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
               DO I=1,NDOFN+1
                  I_CRS3(I) = I_CRS1(I)
               ENDDO
               DO I=1,NTERM_CRS3
                  J_CRS3(I) = J_CRS1(I)
                    CRS3(I) =   CRS1(I)
               ENDDO
               SYM_CRS3 = 'N'

            ELSE                                           !      Error - incorrect SPARSTOR

               WRITE(ERR,932) SUBR_NAME, SPARSTOR
               WRITE(F06,932) SUBR_NAME, SPARSTOR
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF
                                                           ! II-8, sparse add to get CRS2 = KNN + CRS3 = KNN + GMNt*KMM*GMN
            CALL MATADD_SSS_NTERM ( NDOFN, 'KNN-bar + KNM*GMN + (KNM*GMN)t', NTERM_KNN, I_KNN, J_KNN, SYM_KNN, 'GMNt*KMM*GMN',     &
                                    NTERM_CRS3, I_CRS3, J_CRS3, SYM_CRS3, 'CRS2', NTERM_CRS2 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFN, NTERM_CRS2, SUBR_NAME )
            CALL MATADD_SSS ( NDOFN, 'KNN-bar + KNM*GMN + (KNM*GMN)t' , NTERM_KNN ,I_KNN, J_KNN, KNN, ONE, 'GMNt*KMM*GMN',         &
                              NTERM_CRS3, I_CRS3, J_CRS3, CRS3, ONE, 'CRS2', NTERM_CRS2, I_CRS2, J_CRS2, CRS2 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! II-9, deallocate CRS1
            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! II-10, deallocate CRS3

            NTERM_KNN = NTERM_CRS2                         ! II-11, reallocate KNN to be size of CRS2
            WRITE(SC1, * ) '    Reallocate KNN'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KNN', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KNN' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KNN', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KNN', NDOFN, NTERM_KNN, SUBR_NAME )

            DO I=1,NDOFN+1                                 ! II-12, set reduced KNN = CRS2 until we see if KNM is null, below
               I_KNN(I) = I_CRS2(I)
            ENDDO
            DO I=1,NTERM_KNN
               J_KNN(I) = J_CRS2(I)
                 KNN(I) =   CRS2(I)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! II-13, Deallocate CRS2
            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

         ELSE

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

            IF (NTERM_KMN > 0) THEN                        ! Set HMN - KMN if KMN nonzero, else HMN is null
               NTERM_HMN = NTERM_KMN
               CALL ALLOCATE_SPARSE_MAT ( 'HMN', NDOFM, NTERM_HMN, SUBR_NAME )
               DO I=1,NDOFM+1
                  I_HMN(I) = I_KMN(I)
               ENDDO
               DO I=1,NTERM_HMN
                  J_HMN(I) = J_KMN(I)
                  HMN(I) =   KMN(I)
               ENDDO
            ELSE
               NTERM_HMN = 0
            ENDIF

         ENDIF

         IF (NTERM_HMN > 0) THEN
            CALL WRITE_MATRIX_1 ( LINK2J, L2J, 'Y', 'KEEP', L2J_MSG, 'HMN', NTERM_HMN, NDOFM, I_HMN, J_HMN, HMN )
         ENDIF

         WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GMNt', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'GMNt' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate HMN ', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'HMN' )
         CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (MATSPARS == 'N') THEN                       ! Reduce KGG to KNN using full matrix operations

         CALL ALLOCATE_FULL_MAT ( 'KNN_FULL', NDOFN, NDOFN, SUBR_NAME )

         IF (NTERM_KNN > 0) THEN                           ! Put KNN-bar into KNN_FULL
            CALL SPARSE_CRS_TO_FULL ( 'KNN', NTERM_KNN, NDOFN, NDOFN, SYM_KNN, I_KNN, J_KNN, KNN, KNN_FULL )
         ENDIF

         IF (NTERM_KNM > 0) THEN                           ! Part 1: calc KNM*GMN and add it & it's transpose to KNN_FULL

            CALL ALLOCATE_FULL_MAT ( 'KNM_FULL', NDOFN, NDOFM, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'KNM', NTERM_KNM, NDOFN, NDOFM, SYM_KNM, I_KNM, J_KNM, KNM, KNM_FULL )

            CALL ALLOCATE_FULL_MAT ( 'GMN_FULL', NDOFM, NDOFN, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'GMN', NTERM_GMN, NDOFM, NDOFN, SYM_GMN, I_GMN, J_GMN, GMN, GMN_FULL )

            CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFN, NDOFN, SUBR_NAME )

            CALL MATMULT_FFF (  KNM_FULL, GMN_FULL, NDOFN, NDOFM, NDOFN, DUM1 )

            CALL DEALLOCATE_FULL_MAT ( 'GMN_FULL' )
            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFN, NDOFN, SUBR_NAME )

            ITRNSPB = 0                                    ! Calc KNN-bar + KNM*GMN = KNN + DUM1 and put into DUM1
            CALL MATADD_FFF ( KNN_FULL, DUM1, NDOFN, NDOFN, ALPHA, BETA, ITRNSPB, DUM2 )
            ITRNSPB = 1                                    ! Calc KNN-bar + KNM*GMN + (KNM*GMN)' = DUM2+DUM1' and put into KNN_FULL
            CALL MATADD_FFF ( DUM2    , DUM1, NDOFN, NDOFN, ALPHA, BETA, ITRNSPB, KNN_FULL )


            CALL DEALLOCATE_FULL_MAT ( 'DUM1' )
            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

         ENDIF

         IF (NTERM_KMM > 0) THEN                           ! Part 2: calc GMN(t)*KMM*GMN and add to KNN_FULL

            CALL ALLOCATE_FULL_MAT ( 'KMM_FULL', NDOFM, NDOFM, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'KMM', NTERM_KMM, NDOFM, NDOFM, SYM_KMM, I_KMM, J_KMM, KMM, KMM_FULL )

            CALL ALLOCATE_FULL_MAT ( 'GMN_FULL', NDOFM, NDOFN, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'GMN', NTERM_GMN, NDOFM, NDOFN, SYM_GMN, I_GMN, J_GMN, GMN, GMN_FULL )

            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFM, NDOFN, SUBR_NAME )

            CALL MATMULT_FFF ( KMM_FULL, GMN_FULL, NDOFM, NDOFM, NDOFN, DUM2 )

            CALL DEALLOCATE_FULL_MAT ( 'KMM_FULL' )

            CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFN, NDOFN, SUBR_NAME )

            CALL MATMULT_FFF_T (GMN_FULL, DUM2, NDOFM, NDOFN, NDOFN, DUM1 )

            CALL DEALLOCATE_FULL_MAT ( 'GMN_FULL' )
            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

            ITRNSPB = 0                                    ! Add GMN'*KMM*GMN to what is already in KNN_FULL (from above)
            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFN, NDOFN, SUBR_NAME )
            CALL MATADD_FFF ( KNN_FULL, DUM1, NDOFN, NDOFN, ALPHA, BETA, ITRNSPB, DUM2 )
            DO I=1,NDOFN
               DO J=1,NDOFN
                  KNN_FULL(I,J) = DUM2(I,J)
               ENDDO
            ENDDO

            CALL DEALLOCATE_FULL_MAT ( 'DUM1' )
            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

         ENDIF

         CALL CNT_NONZ_IN_FULL_MAT ( 'KNN_FULL  ', KNN_FULL, NDOFN, NDOFN, SYM_KNN, NTERM_KNN, SMALL )

         WRITE(SC1, * ) '    Reallocate KNN'
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KNN', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'KNN' )
         WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KNN', CR13
         CALL ALLOCATE_SPARSE_MAT ( 'KNN', NDOFN, NTERM_KNN, SUBR_NAME )

         IF (NTERM_KNN > 0) THEN                           ! Create new sparse arrays from KNN_FULL
!xx         CALL FULL_TO_SPARSE_CRS ( 'KNN_FULL  ', NDOFN, NDOFN, KNN_FULL, NTERM_KNN, SYM_KNN, I_KNN, J_KNN, KNN )
            CALL DEALLOCATE_FULL_MAT ( 'KNN_FULL' )
         ENDIF

      ELSE

         WRITE(ERR,911) SUBR_NAME,MATSPARS
         WRITE(F06,911) SUBR_NAME,MATSPARS
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
  911 FORMAT(' *ERROR   911: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER MATSPARS MUST BE EITHER "Y" OR "N" BUT VALUE IS ',A)

  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

  936 FORMAT(' *ERROR   936: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF NONZERO TERMS IN SPARSE MATRIX MMN = ',I12,' IS NOT EQUAL TO THOSE IN MATRIX MNM = ',I12)

12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_KGG_TO_KNN


      SUBROUTINE REDUCE_KGGD_TO_KNND ( PART_VEC_G_NM )

! Call routines to reduce the KGGD differential stiffness matrix from the G-set to the N, M-sets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, LINK2A, L2A, L2ASTAT, L2A_MSG, L2J, LINK2J, L2J_MSG, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NDOFN, NDOFM, NTERM_HMN, NTERM_KGGD, NTERM_KNND,          &
                                         NTERM_KNMD, NTERM_KMMD, NTERM_GMN
      USE PARAMS, ONLY                :  EPSIL, SPARSTOR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE SPARSE_MATRICES, ONLY       :  I_HMN, J_HMN, HMN, I_KGGD, J_KGGD, KGGD, I_KNND, J_KNND, KNND, I_KNMD, J_KNMD, KNMD,      &
                                         I_KMMD, J_KMMD, KMMD, I_KMND, J_KMND, KMND, I_GMN, J_GMN, GMN,  I_GMNt, J_GMNt, GMNt
      USE SPARSE_MATRICES, ONLY       :  SYM_GMN, SYM_HMN, SYM_KGGD, SYM_KNND, SYM_KNMD, SYM_KMMD, SYM_KMND
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM, MATTRNSP_SS
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_NONSYM_TO_CRS_SYM, SPARSE_CRS_SPARSE_CCS
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_CRS_TERM_COUNT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  READ_MATRIX_1, WRITE_MATRIX_1

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_KGGD_TO_KNND'
      CHARACTER(  1*BYTE)             :: SYM_CRS1            ! Storage format for matrix CRS1 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)
      CHARACTER(  1*BYTE)             :: SYM_CRS3            ! Storage format for matrix CRS3 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_G_NM(NDOFG)! Partitioning vector (G set into N and M sets)
      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
!                                                              the ones on and above the diagonal (controlled by param SPARSTOR)
      INTEGER(LONG)                   :: KNND_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KNMD_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
!xx   INTEGER(LONG)                   :: KMND_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KMMD_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: NTERM_CCS1          ! Number of terms in matrix CCS1
      INTEGER(LONG)                   :: NTERM_CRS1          ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2          ! Number of terms in matrix CRS2
      INTEGER(LONG)                   :: NTERM_CRS3          ! Number of terms in matrix CRS3
      INTEGER(LONG)                   :: NTERM_KMND           ! Number of nonzeros in sparse matrix KMND (should = NTERM_KNMD)
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices


      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Partition KNND from KGGD (This is KNND before reduction, or KNND(bar) )

      IF (NDOFN > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KGGD', NTERM_KGGD, NDOFG, NDOFG, SYM_KGGD, I_KGGD, J_KGGD,      PART_VEC_G_NM, PART_VEC_G_NM,  &
                                    NUM1, NUM1, KNND_ROW_MAX_TERMS, 'KNND', NTERM_KNND, SYM_KNND )

         CALL ALLOCATE_SPARSE_MAT ( 'KNND', NDOFN, NTERM_KNND, SUBR_NAME )

         IF (NTERM_KNND > 0) THEN
            CALL PARTITION_SS ( 'KGGD', NTERM_KGGD, NDOFG, NDOFG, SYM_KGGD, I_KGGD, J_KGGD, KGGD, PART_VEC_G_NM, PART_VEC_G_NM,    &
                                 NUM1, NUM1, KNND_ROW_MAX_TERMS, 'KNND', NTERM_KNND, NDOFN, SYM_KNND, I_KNND, J_KNND, KNND )
         ENDIF

      ENDIF

! Partition KNMD from KGGD

      IF ((NDOFN > 0) .AND. (NDOFM > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'KGGD', NTERM_KGGD, NDOFG, NDOFG, SYM_KGGD, I_KGGD, J_KGGD,      PART_VEC_G_NM, PART_VEC_G_NM,  &
                                    NUM1, NUM2, KNMD_ROW_MAX_TERMS, 'KNMD', NTERM_KNMD, SYM_KNMD )

         CALL ALLOCATE_SPARSE_MAT ( 'KNMD', NDOFN, NTERM_KNMD, SUBR_NAME )

         IF (NTERM_KNMD > 0) THEN
            CALL PARTITION_SS ( 'KGGD', NTERM_KGGD, NDOFG, NDOFG, SYM_KGGD, I_KGGD, J_KGGD, KGGD, PART_VEC_G_NM, PART_VEC_G_NM,    &
                                 NUM1, NUM2, KNMD_ROW_MAX_TERMS, 'KNMD', NTERM_KNMD, NDOFN, SYM_KNMD, I_KNMD, J_KNMD, KNMD )
         ENDIF

      ENDIF

! Partition KMND from KGGD

      IF ((NDOFN > 0) .AND. (NDOFM > 0)) THEN

!xx      CALL PARTITION_SS_NTERM ( 'KGGD', NTERM_KGGD, NDOFG, NDOFG, SYM_KGGD, I_KGGD, J_KGGD,      PART_VEC_G_NM, PART_VEC_G_NM,  &
!xx                                 NUM2, NUM1, KMND_ROW_MAX_TERMS, 'KMND', NTERM_KMND, SYM_KMND )

!xx      IF (NTERM_KMND /= NTERM_KNMD) THEN
!xx         FATAL_ERR = FATAL_ERR + 1
!xx         WRITE(ERR,936) SUBR_NAME, NTERM_KMND, NTERM_KNMD
!xx         WRITE(F06,936) SUBR_NAME, NTERM_KMND, NTERM_KNMD
!xx         CALL OUTA_HERE ( 'Y' )
!xx      ENDIF

!xx      CALL ALLOCATE_SPARSE_MAT ( 'KMND', NDOFM, NTERM_KMND, SUBR_NAME )

!xx      IF (NTERM_KMND > 0) THEN
!xx         CALL PARTITION_SS ( 'KGGD', NTERM_KGGD, NDOFG, NDOFG, SYM_KGGD, I_KGGD, J_KGGD, KGGD, PART_VEC_G_NM, PART_VEC_G_NM,    &
!xx                              NUM2, NUM1, KMND_ROW_MAX_TERMS, 'KMND', NTERM_KMND, NDOFM, SYM_KMND, I_KMND, J_KMND, KMND )
!xx      ENDIF

         NTERM_KMND = NTERM_KNMD
         CALL ALLOCATE_SPARSE_MAT ( 'KMND', NDOFM, NTERM_KNMD, SUBR_NAME )

         IF (NTERM_KMND > 0) THEN
            CALL MATTRNSP_SS ( NDOFN, NDOFM, NTERM_KNMD, 'KNMD', I_KNMD, J_KNMD, KNMD, 'KMND', I_KMND, J_KMND, KMND )

         ENDIF

      ENDIF

! Partition KMMD from KGGD

      IF (NDOFM > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KGGD', NTERM_KGGD, NDOFG, NDOFG, SYM_KGGD, I_KGGD, J_KGGD,      PART_VEC_G_NM, PART_VEC_G_NM,  &
                                    NUM2, NUM2, KMMD_ROW_MAX_TERMS, 'KMMD', NTERM_KMMD, SYM_KMMD )

         CALL ALLOCATE_SPARSE_MAT ( 'KMMD', NDOFM, NTERM_KMMD, SUBR_NAME )

         IF (NTERM_KMMD > 0) THEN
            CALL PARTITION_SS ( 'KGGD', NTERM_KGGD, NDOFG, NDOFG, SYM_KGGD, I_KGGD, J_KGGD, KGGD, PART_VEC_G_NM, PART_VEC_G_NM,    &
                                 NUM2, NUM2, KMMD_ROW_MAX_TERMS, 'KMMD', NTERM_KMMD, NDOFM, SYM_KMMD, I_KMMD, J_KMMD, KMMD )
         ENDIF

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Reduce KGGD to KNND = KNND(bar) + KNMD*GMN + (KNMD*GMN)' + GMN'*KMMD*GMN.
! If PARAM MATSPARS = 'Y', then we use sparse matrix operations (multiply/add/transpose). If not, use full matrix operations


         IF (.NOT. ALLOCATED(GMN)) THEN
            CALL ALLOCATE_SPARSE_MAT ( 'GMN' , NDOFN, NTERM_GMN, SUBR_NAME )

            CALL READ_MATRIX_1 ( LINK2A, L2A, 'N', 'N', L2ASTAT, L2A_MSG, &
                                 'GMN', NTERM_GMN, 'Y', NDOFM,           &
                                 I_GMN, J_GMN, GMN )
         ENDIF
         CALL ALLOCATE_SPARSE_MAT ( 'GMNt', NDOFN, NTERM_GMN, SUBR_NAME )
         CALL MATTRNSP_SS ( NDOFM, NDOFN, NTERM_GMN, 'GMN', I_GMN, J_GMN, GMN, 'GMNt', I_GMNt, J_GMNt, GMNt )

                                                           ! CCS1 will be sparse CCS format version of sparse CRS matrix GMN
         CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFN, NTERM_GMN, SUBR_NAME )
         CALL SPARSE_CRS_SPARSE_CCS ( NDOFM, NDOFN, NTERM_GMN, 'GMN', I_GMN, J_GMN, GMN, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

         IF (NTERM_KNMD > 0) THEN                           ! Part I of reduced KNND: calc KNMD*GMN & add & transpose to orig KNND

                                                           ! I-1, sparse multiply to get CRS1 = KNMD*GMN. Use CCS1 for GMN CCS
            CALL MATMULT_SSS_NTERM ( 'KNMD', NDOFN, NTERM_KNMD, SYM_KNMD, I_KNMD, J_KNMD,                                          &
                                     'GMN' , NDOFN, NTERM_GMN , SYM_GMN , J_CCS1, I_CCS1, AROW_MAX_TERMS,                          &
                                     'CRS1',        NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'KNMD' , NDOFN, NTERM_KNMD , SYM_KNMD, I_KNMD , J_KNMD , KNMD ,                                     &
                               'GMN' , NDOFN, NTERM_GMN , SYM_GMN, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            NTERM_CRS2 = NTERM_CRS1                        ! I-2, allocate memory to array CRS2 which will hold transpose of CRS1
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFN, NTERM_CRS2, SUBR_NAME )

                                                           ! I-3, transpose CRS1 to get CRS2 = (KNMD*GMN)t
            CALL MATTRNSP_SS ( NDOFN, NDOFN, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CRS2', I_CRS2, J_CRS2, CRS2 )

                                                           ! I-4, sparse add to get CRS3 = CRS1 + CRS2 = (KNMD*GMN) + (KNMD*GMN)t
            CALL MATADD_SSS_NTERM (NDOFN,'KNMD*GMN', NTERM_CRS1, I_CRS1, J_CRS1,'N','(KNMD*GMN)t', NTERM_CRS2, I_CRS2, J_CRS2, 'N',&
                                         'CRS3', NTERM_CRS3)
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
            CALL MATADD_SSS ( NDOFN, 'KNMD*GMN', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE, '(KNMD*GMN)t', NTERM_CRS2,                 &
                              I_CRS2, J_CRS2, CRS2, ONE, 'CRS3', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! I-5, deallocate CRS1 which was KNMD*GMN
            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! I-6, deallocate CRS2 which was (KNMD*GMN)t

                                                           ! I-7, CRS3 = (KNMD*GMN) + (KNMD*GMN)t has all nonzero terms in it.
            IF      (SPARSTOR == 'SYM   ') THEN            !      If SPARSTOR == 'SYM   ', rewrite CRS3 as sym in CRS1

               CALL SPARSE_CRS_TERM_COUNT ( NDOFN, NTERM_CRS3, '(KNMD*GMN) + (KNMD*GMN)t', I_CRS3, J_CRS3, NTERM_CRS1 )
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )
               CALL CRS_NONSYM_TO_CRS_SYM ('CRS3 = (KNMD*GMN) + (KNMD*GMN)t all nonzeros', NDOFN, NTERM_CRS3, I_CRS3, J_CRS3, CRS3,&
                                           'CRS1 = (KNMD*GMN) + (KNMD*GMN)t stored sym'  ,        NTERM_CRS1, I_CRS1, J_CRS1, CRS1 )
               SYM_CRS1 = 'Y'

            ELSE IF (SPARSTOR == 'NONSYM') THEN            !      If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

               NTERM_CRS1 = NTERM_CRS3
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )
               DO I=1,NDOFN+1
                  I_CRS1(I) = I_CRS3(I)
               ENDDO
               DO I=1,NTERM_CRS1
                  J_CRS1(I) = J_CRS3(I)
                    CRS1(I) =   CRS3(I)
               ENDDO
               SYM_CRS1 = 'N'

            ELSE                                           !      Error - incorrect SPARSTOR

               WRITE(ERR,932) SUBR_NAME, SPARSTOR
               WRITE(F06,932) SUBR_NAME, SPARSTOR
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF

            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! I-8, deallocate CRS3 which was KNMD*GMN + (KNMD*GMN)t. Now is in CRS1

                                                           ! I-9, add: CRS3 = KNND-bar + CRS1 = KNND-bar + KNMD*GMN + (KNMD*GMN)t
            CALL MATADD_SSS_NTERM ( NDOFN, 'KNND-bar', NTERM_KNND , I_KNND , J_KNND , SYM_KNND , 'KNMD*GMN + (KNMD*GMN)t',         &
                                                      NTERM_CRS1, I_CRS1, J_CRS1, SYM_CRS1, 'CRS3', NTERM_CRS3 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
            CALL MATADD_SSS ( NDOFN, 'KNND-bar', NTERM_KNND, I_KNND, J_KNND, KNND, ONE, 'KNMD*GMN + (KNMD*GMN)t', NTERM_CRS1,      &
                                      I_CRS1, J_CRS1, CRS1, ONE, 'CRS1', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! I-10, deallocate CRS1 = KNMD*GMN + (KNMD*GMN)t

            NTERM_KNND = NTERM_CRS3                        ! I-11, reallocate KNND to be size of CRS3
            WRITE(SC1, * ) '    Reallocate KNND'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KNND', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KNND' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KNND', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KNND', NDOFN, NTERM_KNND, SUBR_NAME )

            DO I=1,NDOFN+1                                 ! I-12, set KNND = CRS3
               I_KNND(I) = I_CRS3(I)
            ENDDO
            DO J=1,NTERM_KNND
               J_KNND(J) = J_CRS3(J)
                 KNND(J) =   CRS3(J)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! I-13, deallocate CRS3
                                                           ! At this point, CRS1, CRS2, CRS3 are deallocated, CCS1 is being used
         ENDIF

         IF (NTERM_KMMD > 0) THEN                           ! Part II of reduced KNND: calc GMN(t)*KMMD*GMN and add to KNND

                                                           ! II-1, sparse multiply to get CRS1 = KMMD*GMN using CCS1 for GMN CCS

            CALL MATMULT_SSS_NTERM ( 'KMMD', NDOFM, NTERM_KMMD, SYM_KMMD, I_KMMD, J_KMMD ,                                         &
                                     'GMN' , NDOFN, NTERM_GMN , SYM_GMN , J_CCS1, I_CCS1, AROW_MAX_TERMS,                          &
                                     'CRS1',        NTERM_CRS1 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFM, NTERM_CRS1, SUBR_NAME )
            CALL MATMULT_SSS ( 'KMMD' , NDOFM, NTERM_KMMD , SYM_KMMD, I_KMMD , J_KMMD , KMMD ,                                     &
                               'GMN' , NDOFN, NTERM_GMN , SYM_GMN, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            IF (NTERM_KMND > 0) THEN
               CALL MATADD_SSS_NTERM (NDOFM, 'KMMD*GMN', NTERM_CRS1, I_CRS1, J_CRS1, 'N', 'KMND', NTERM_KMND, I_KMND, J_KMND,      &
                                              SYM_KMND, 'HMN', NTERM_HMN)

               CALL ALLOCATE_SPARSE_MAT ( 'HMN', NDOFM, NTERM_HMN, SUBR_NAME )

               CALL MATADD_SSS ( NDOFM, 'KMMD*GM', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE, 'KMND', NTERM_KMND, I_KMND, J_KMND, KMND,&
                                 ONE, 'HMN', NTERM_HMN, I_HMN, J_HMN, HMN )
            ELSE
               NTERM_HMN = NTERM_CRS1
               CALL ALLOCATE_SPARSE_MAT ( 'HMN', NDOFM, NTERM_HMN, SUBR_NAME )
               DO I=1,NDOFM+1
                  I_HMN(I) = I_CRS1(I)
               ENDDO
               DO I=1,NTERM_HMN
                  J_HMN(I) = J_CRS1(I)
                  HMN(I) =   CRS1(I)
               ENDDO
            ENDIF


            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )             ! II-2, deallocate CCS1 which was CCS version of GMN

            NTERM_CCS1 = NTERM_CRS1                        ! II-3, allocate CCS1 to be same as CRS1 = KMMD*GMN but in CCS format
            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFN, NTERM_CCS1, SUBR_NAME )
            CALL SPARSE_CRS_SPARSE_CCS ( NDOFM, NDOFN, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! II-4, deallocate CRS1 which was KMMD*GMN
                                                           ! II-5, sparse multiply to get CRS1 = GMNt*CCS1  with CCS1 = KMMD*GMN
!                                                           (note: use SYM_GMN for sym indicator of KMMD*GMN)
            CALL MATMULT_SSS_NTERM ( 'GMNt', NDOFN, NTERM_GMN , SYM_GMN, I_GMNt, J_GMNt,                                           &
                                     'CCS1', NDOFN, NTERM_CCS1, SYM_GMN, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1',        NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'GMNt', NDOFN, NTERM_GMN , SYM_GMN, I_GMNt, J_GMNt, GMNt,                                           &
                               'CCS1', NDOFN, NTERM_CCS1, SYM_GMN, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )             ! II-6, deallocate CCS1

                                                           ! II-7, CRS1 = GMNt*KMMD*GMN has all nonzero terms in it.
            IF      (SPARSTOR == 'SYM   ') THEN            !      If SPARSTOR == 'SYM   ', rewrite CRS1 as sym in CRS3

               CALL SPARSE_CRS_TERM_COUNT ( NDOFN, NTERM_CRS1, 'GMNt*KMMD*GMN all nonzeros', I_CRS1, J_CRS1, NTERM_CRS3 )
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
               CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS1 = GMNt*KMMD*GMN all nonzeros', NDOFN, NTERM_CRS1, I_CRS1, J_CRS1, CRS1,          &
                                            'CRS3 = GMNt*KMMD*GMN stored sym'  ,        NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )
               SYM_CRS3 = 'Y'

            ELSE IF (SPARSTOR == 'NONSYM') THEN            !      If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

               NTERM_CRS3 = NTERM_CRS1
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
               DO I=1,NDOFN+1
                  I_CRS3(I) = I_CRS1(I)
               ENDDO
               DO I=1,NTERM_CRS3
                  J_CRS3(I) = J_CRS1(I)
                    CRS3(I) =   CRS1(I)
               ENDDO
               SYM_CRS3 = 'N'

            ELSE                                           !      Error - incorrect SPARSTOR

               WRITE(ERR,932) SUBR_NAME, SPARSTOR
               WRITE(F06,932) SUBR_NAME, SPARSTOR
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF
                                                           ! II-8, sparse add to get CRS2 = KNND + CRS3 = KNND + GMNt*KMMD*GMN
            CALL MATADD_SSS_NTERM ( NDOFN, 'KNND-bar + KNMD*GMN + (KNMD*GMN)t', NTERM_KNND, I_KNND, J_KNND, SYM_KNND,              &
                                   'GMNt*KMMD*GMN',NTERM_CRS3, I_CRS3, J_CRS3, SYM_CRS3, 'CRS2', NTERM_CRS2 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFN, NTERM_CRS2, SUBR_NAME )
            CALL MATADD_SSS ( NDOFN, 'KNND-bar + KNMD*GMN + (KNMD*GMN)t' , NTERM_KNND ,I_KNND, J_KNND, KNND, ONE, 'GMNt*KMMD*GMN', &
                              NTERM_CRS3, I_CRS3, J_CRS3, CRS3, ONE, 'CRS2', NTERM_CRS2, I_CRS2, J_CRS2, CRS2 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! II-9, deallocate CRS1
            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! II-10, deallocate CRS3

            NTERM_KNND = NTERM_CRS2                        ! II-11, reallocate KNND to be size of CRS2
            WRITE(SC1, * ) '    Reallocate KNND'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KNND', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KNND' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KNND', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KNND', NDOFN, NTERM_KNND, SUBR_NAME )

            DO I=1,NDOFN+1                                 ! II-12, set reduced KNND = CRS2 until we see if KNMD is null, below
               I_KNND(I) = I_CRS2(I)
            ENDDO
            DO I=1,NTERM_KNND
               J_KNND(I) = J_CRS2(I)
                 KNND(I) =   CRS2(I)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! II-13, Deallocate CRS2
            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

         ELSE

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

            IF (NTERM_KMND > 0) THEN                        ! Set HMN - KMND if KMND nonzero, else HMN is null
               NTERM_HMN = NTERM_KMND
               CALL ALLOCATE_SPARSE_MAT ( 'HMN', NDOFM, NTERM_HMN, SUBR_NAME )
               DO I=1,NDOFM+1
                  I_HMN(I) = I_KMND(I)
               ENDDO
               DO I=1,NTERM_HMN
                  J_HMN(I) = J_KMND(I)
                  HMN(I) =   KMND(I)
               ENDDO
            ELSE
               NTERM_HMN = 0
            ENDIF

         ENDIF

         IF (NTERM_HMN > 0) THEN
            CALL WRITE_MATRIX_1 ( LINK2J, L2J, 'Y', 'KEEP', L2J_MSG, 'HMN', NTERM_HMN, NDOFM, I_HMN, J_HMN, HMN )
         ENDIF

         WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GMNt', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'GMNt' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate HMN ', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'HMN' )
         CALL DEALLOCATE_SCR_MAT ( 'CCS1' )




      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

  936 FORMAT(' *ERROR   936: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF NONZERO TERMS IN SPARSE MATRIX MMN = ',I12,' IS NOT EQUAL TO THOSE IN MATRIX MNM = ',I12)

12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_KGGD_TO_KNND


      SUBROUTINE REDUCE_MGG_TO_MNN ( PART_VEC_G_NM )

! Call routines to reduce the MGG mass matrix from the G-set to the N, M-sets. See Appendix B to the MYSTRAN User's Reference Manual
! for the derivation of the reduction equations.

! NOTE: This subr has code for sparse matrices as well as full matrices, (i.e. Bulk Data PARAM MATSPARS = 'Y' for sparse and 'N'
! for full). The code for full matrices was put in originally in order that the sparse code could be thoroughly checked. That task
! is complete and the remaining full matrix code has not been maintained since around 2005. In addition, new capability added to
! MYSTRAN since that approx time does not have full matrix code.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L2R, LINK2R, L2R_MSG, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NDOFN, NDOFM, NTERM_MGG, NTERM_MNN, NTERM_MNM, NTERM_MMM, &
                                         NTERM_GMN, NTERM_LMN
      USE PARAMS, ONLY                :  EPSIL, MATSPARS, SPARSTOR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE SPARSE_MATRICES, ONLY       :  I_LMN, J_LMN, LMN, I_MGG, J_MGG, MGG, I_MNN, J_MNN, MNN, I_MNM , J_MNM , MNM ,            &
                                         I_MMN, J_MMN, MMN, I_MMM, J_MMM, MMM, I_GMN, J_GMN, GMN, I_GMNt, J_GMNt, GMNt
      USE SPARSE_MATRICES, ONLY       :  SYM_GMN, SYM_LMN, SYM_MGG, SYM_MNN, SYM_MNM, SYM_MMN, SYM_MMM
      USE FULL_MATRICES, ONLY         :  MNN_FULL, MNM_FULL, MMM_FULL, GMN_FULL, DUM1, DUM2, DUM3
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM, MATTRNSP_SS
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_NONSYM_TO_CRS_SYM, SPARSE_CRS_SPARSE_CCS, SPARSE_CRS_TO_FULL
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_CRS_TERM_COUNT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_MATRIX_1
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATADD_FFF, MATMULT_FFF, MATMULT_FFF_T
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CNT_NONZ_IN_FULL_MAT

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_MGG_TO_MNN'
      CHARACTER(  1*BYTE)             :: SYM_CRS1            ! Storage format for matrix CRS1 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)
      CHARACTER(  1*BYTE)             :: SYM_CRS3            ! Storage format for matrix CRS3 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_G_NM(NDOFG)! Partitioning vector (G set into N and M sets)
      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG)                   :: ITRNSPB             ! Transpose indicator for matrix multiply routine
      INTEGER(LONG)                   :: MNN_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: MNM_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
!xx   INTEGER(LONG)                   :: MMN_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: MMM_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: NTERM_CCS1          ! Number of terms in matrix CCS1
      INTEGER(LONG)                   :: NTERM_CRS1          ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2          ! Number of terms in matrix CRS2
      INTEGER(LONG)                   :: NTERM_CRS3          ! Number of terms in matrix CRS3
      INTEGER(LONG)                   :: NTERM_MMN           ! Number of nonzeros in sparse matrix MMN (should = NTERM_MNM)
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices


      REAL(DOUBLE)                    :: ALPHA = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: BETA  = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: SMALL             ! A number used in filtering out small numbers from a full matrix

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Partition MNN from MGG (This is MNN before reduction, or MNN(bar) )

      IF (NDOFN > 0) THEN

         CALL PARTITION_SS_NTERM ( 'MGG', NTERM_MGG, NDOFG, NDOFG, SYM_MGG, I_MGG, J_MGG,      PART_VEC_G_NM, PART_VEC_G_NM,       &
                                    NUM1, NUM1, MNN_ROW_MAX_TERMS, 'MNN', NTERM_MNN, SYM_MNN )

         CALL ALLOCATE_SPARSE_MAT ( 'MNN', NDOFN, NTERM_MNN, SUBR_NAME )

         IF (NTERM_MNN > 0) THEN
            CALL PARTITION_SS ( 'MGG', NTERM_MGG, NDOFG, NDOFG, SYM_MGG, I_MGG, J_MGG, MGG, PART_VEC_G_NM, PART_VEC_G_NM,          &
                                 NUM1, NUM1, MNN_ROW_MAX_TERMS, 'MNN', NTERM_MNN, NDOFN, SYM_MNN, I_MNN, J_MNN, MNN )
         ENDIF

      ENDIF

! Partition MNM from MGG

      IF ((NDOFN > 0) .AND. (NDOFM > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'MGG', NTERM_MGG, NDOFG, NDOFG, SYM_MGG, I_MGG, J_MGG,      PART_VEC_G_NM, PART_VEC_G_NM,       &
                                    NUM1, NUM2, MNM_ROW_MAX_TERMS, 'MNM', NTERM_MNM, SYM_MNM )

         CALL ALLOCATE_SPARSE_MAT ( 'MNM', NDOFN, NTERM_MNM, SUBR_NAME )

         IF (NTERM_MNM > 0) THEN
            CALL PARTITION_SS ( 'MGG', NTERM_MGG, NDOFG, NDOFG, SYM_MGG, I_MGG, J_MGG, MGG, PART_VEC_G_NM, PART_VEC_G_NM,          &
                                 NUM1, NUM2, MNM_ROW_MAX_TERMS, 'MNM', NTERM_MNM, NDOFN, SYM_MNM, I_MNM, J_MNM, MNM )
         ENDIF

      ENDIF

! Partition MMN from MGG

      IF ((NDOFN > 0) .AND. (NDOFM > 0)) THEN

!xx      CALL PARTITION_SS_NTERM ( 'MGG', NTERM_MGG, NDOFG, NDOFG, SYM_MGG, I_MGG, J_MGG,      PART_VEC_G_NM, PART_VEC_G_NM,       &
!xx                                 NUM2, NUM1, MMN_ROW_MAX_TERMS, 'MMN', NTERM_MMN, SYM_MMN )

!xx      IF (NTERM_MMN /= NTERM_MNM) THEN
!xx         FATAL_ERR = FATAL_ERR + 1
!xx         WRITE(ERR,936) SUBR_NAME, NTERM_MMN, NTERM_MNM
!xx         WRITE(F06,936) SUBR_NAME, NTERM_MMN, NTERM_MNM
!xx         CALL OUTA_HERE ( 'Y' )
!xx      ENDIF

!xx      CALL ALLOCATE_SPARSE_MAT ( 'MMN', NDOFN, NTERM_MMN, SUBR_NAME )

!xx      IF (NTERM_MMN > 0) THEN
!xx         CALL PARTITION_SS ( 'MGG', NTERM_MGG, NDOFG, NDOFG, SYM_MGG, I_MGG, J_MGG, MGG, PART_VEC_G_NM, PART_VEC_G_NM,          &
!xx                              NUM2, NUM1, MMN_ROW_MAX_TERMS, 'MMN', NTERM_MMN, NDOFN, SYM_MMN, I_MMN, J_MMN, MMN )
!xx      ENDIF

         CALL ALLOCATE_SPARSE_MAT ( 'MMN', NDOFM, NTERM_MNM, SUBR_NAME )

         NTERM_MMN = NTERM_MNM
         IF (NTERM_MMN > 0) THEN
            CALL MATTRNSP_SS ( NDOFN, NDOFM, NTERM_MNM, 'MNM', I_MNM, J_MNM, MNM, 'MMN', I_MMN, J_MMN, MMN )

         ENDIF

      ENDIF

! Partition MMM from MGG

      IF (NDOFM > 0) THEN

         CALL PARTITION_SS_NTERM ( 'MGG', NTERM_MGG, NDOFG, NDOFG, SYM_MGG, I_MGG, J_MGG,      PART_VEC_G_NM, PART_VEC_G_NM,       &
                                    NUM2, NUM2, MMM_ROW_MAX_TERMS, 'MMM', NTERM_MMM, SYM_MMM )

         CALL ALLOCATE_SPARSE_MAT ( 'MMM', NDOFM, NTERM_MMM, SUBR_NAME )

         IF (NTERM_MMM > 0) THEN
            CALL PARTITION_SS ( 'MGG', NTERM_MGG, NDOFG, NDOFG, SYM_MGG, I_MGG, J_MGG, MGG, PART_VEC_G_NM, PART_VEC_G_NM,          &
                                 NUM2, NUM2, MMM_ROW_MAX_TERMS, 'MMM', NTERM_MMM, NDOFM, SYM_MMM, I_MMM, J_MMM, MMM )
         ENDIF

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Reduce MGG to MNN = MNN(bar) + MNM*GMN + (MNM*GMN)' + GMN'*MMM*GMN.
! If PARAM MATSPARS = 'Y', then we use sparse matrix operations (multiply/add/transpose). If not, use full matrix operations

      IF (MATSPARS == 'Y') THEN

         CALL ALLOCATE_SPARSE_MAT ( 'GMNt', NDOFN, NTERM_GMN, SUBR_NAME )
         CALL MATTRNSP_SS ( NDOFM, NDOFN, NTERM_GMN, 'GMN', I_GMN, J_GMN, GMN, 'GMNt', I_GMNt, J_GMNt, GMNt )

                                                           ! CCS1 will be sparse CCS format version of sparse CRS matrix GMN
         CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFN, NTERM_GMN, SUBR_NAME )
         CALL SPARSE_CRS_SPARSE_CCS ( NDOFM, NDOFN, NTERM_GMN, 'GMN', I_GMN, J_GMN, GMN, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

         IF (NTERM_MNM > 0) THEN                           ! Part I of reduced MNN: calc MNM*GMN & add it & transpose to orig MNN

                                                           ! I-1, sparse multiply to get CRS1 = MNM*GMN. Use CCS1 for GMN CCS
            CALL MATMULT_SSS_NTERM ( 'MNM' , NDOFN, NTERM_MNM , SYM_MNM, I_MNM , J_MNM ,                                           &
                                     'GMN' , NDOFN, NTERM_GMN , SYM_GMN, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1' ,       NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'MNM' , NDOFN, NTERM_MNM , SYM_MNM, I_MNM , J_MNM , MNM ,                                           &
                               'GMN' , NDOFN, NTERM_GMN , SYM_GMN, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            NTERM_CRS2 = NTERM_CRS1                        ! I-2, allocate memory to array CRS2 which will hold transpose of CRS1
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFN, NTERM_CRS2, SUBR_NAME )

                                                           ! I-3, transpose CRS1 to get CRS2 = (MNM*GMN)t
            CALL MATTRNSP_SS ( NDOFN, NDOFN, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CRS2', I_CRS2, J_CRS2, CRS2 )

                                                           ! I-4, sparse add to get CRS3 = CRS1 + CRS2 = (MNM*GMN) + (MNM*GMN)t
            CALL MATADD_SSS_NTERM (NDOFN,'MNM*GMN', NTERM_CRS1, I_CRS1, J_CRS1, 'N', '(MNM*GMN)t', NTERM_CRS2, I_CRS2, J_CRS2, 'N',&
                                         'CRS3', NTERM_CRS3)
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
            CALL MATADD_SSS ( NDOFN, 'MNM*GMN', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE, '(MNM*GMN)t', NTERM_CRS2, I_CRS2, J_CRS2,   &
                              CRS2, ONE, 'CRS3', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! I-5, deallocate CRS1, CRS2
            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! I-6, deallocate CRS2 which was (MNM*GMN)t

                                                           ! I-7, CRS3 = (MNM*GMN) + (MNM*GMN)t has all nonzero terms in it.
            IF      (SPARSTOR == 'SYM   ') THEN            !      If SPARSTOR == 'SYM   ', rewrite CRS3 as sym in CRS1

               CALL SPARSE_CRS_TERM_COUNT ( NDOFN, NTERM_CRS3, '(MNM*GMN) + (MNM*GMN)t', I_CRS3, J_CRS3, NTERM_CRS1 )
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )
               CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS3 = (MNM*GMN) + (MNM*GMN)t all nonzeros', NDOFN, NTERM_CRS3, I_CRS3, J_CRS3, CRS3,&
                                            'CRS1 = (MNM*GMN) + (MNM*GMN)t stored sym'  ,        NTERM_CRS1, I_CRS1, J_CRS1, CRS1 )
               SYM_CRS1 = 'Y'

            ELSE IF (SPARSTOR == 'NONSYM') THEN            !      If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

               NTERM_CRS1 = NTERM_CRS3
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )
               DO I=1,NDOFN+1
                  I_CRS1(I) = I_CRS3(I)
               ENDDO
               DO I=1,NTERM_CRS1
                  J_CRS1(I) = J_CRS3(I)
                    CRS1(I) =   CRS3(I)
               ENDDO
               SYM_CRS1 = 'N'

            ELSE                                           !      Error - incorrect SPARSTOR

               WRITE(ERR,932) SUBR_NAME, SPARSTOR
               WRITE(F06,932) SUBR_NAME, SPARSTOR
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF

            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! I-8, deallocate CRS3

                                                           ! I-9, sparse add to get CRS3 = MNN-bar + CRS1 = (MNM*GMN) + (MNM*GMN)t
            CALL MATADD_SSS_NTERM ( NDOFN, 'MNN-bar', NTERM_MNN, I_MNN, J_MNN, SYM_MNN, 'MNM*GMN + (MNM*GMN)t',                    &
                                    NTERM_CRS1, I_CRS1, J_CRS1, SYM_CRS1, 'CRS3', NTERM_CRS3 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
            CALL MATADD_SSS ( NDOFN, 'MNN-bar', NTERM_MNN, I_MNN, J_MNN, MNN, ONE, 'MNM*GMN + (MNM*GMN)t', NTERM_CRS1,             &
                              I_CRS1, J_CRS1, CRS1, ONE, 'CRS1', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! I-10, deallocate CRS1 and CRS3

            NTERM_MNN = NTERM_CRS3                         ! I-11, reallocate MNN to be size of CRS1
            WRITE(SC1, * ) '    Reallocate MNN'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MNN', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'MNN' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   MNN', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'MNN', NDOFN, NTERM_MNN, SUBR_NAME )
                                                           ! I-12, set MNN = CRS1
            DO I=1,NDOFN+1
               I_MNN(I) = I_CRS3(I)
            ENDDO
            DO J=1,NTERM_MNN
               J_MNN(J) = J_CRS3(J)
                 MNN(J) =   CRS3(J)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! I-13, deallocate CRS3
                                                           ! At this point, CRS1, CRS2, CRS3 are deallocated, CCS1 is being used
         ENDIF

         IF (NTERM_MMM > 0) THEN                           ! Part II of reduced MNN: calc GMN(t)*MMM*GMN and add to MNN

                                                           ! II-1 sparse multiply to get CRS1 = MMM*GMN using CCS1 for GMN CCS
            CALL MATMULT_SSS_NTERM ( 'MMM' , NDOFM, NTERM_MMM , SYM_MMM, I_MMM , J_MMM ,                                           &
                                     'GMN' , NDOFN, NTERM_GMN , SYM_GMN, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1' ,       NTERM_CRS1 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFM, NTERM_CRS1, SUBR_NAME )
            CALL MATMULT_SSS ( 'MMM' , NDOFM, NTERM_MMM , SYM_MMM, I_MMM , J_MMM , MMM ,                                           &
                               'GMN' , NDOFN, NTERM_GMN , SYM_GMN, J_CCS1, I_CCS1, CCS1,  AROW_MAX_TERMS,                          &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )


            IF (NTERM_MMN > 0) THEN
               CALL MATADD_SSS_NTERM (NDOFM, 'MMM*GMN', NTERM_CRS1, I_CRS1, J_CRS1, 'N', 'MMN', NTERM_MMN, I_MMN, J_MMN, SYM_MMN,  &
                                             'LMN', NTERM_LMN)

               CALL ALLOCATE_SPARSE_MAT ( 'LMN', NDOFM, NTERM_LMN, SUBR_NAME )

               CALL MATADD_SSS ( NDOFM, 'MMM*GM', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE, 'MMN', NTERM_MMN, I_MMN, J_MMN, MMN,      &
                                 ONE, 'LMN', NTERM_LMN, I_LMN, J_LMN, LMN )
            ELSE
               NTERM_LMN = NTERM_CRS1
               CALL ALLOCATE_SPARSE_MAT ( 'LMN', NDOFM, NTERM_LMN, SUBR_NAME )
               DO I=1,NDOFM+1
                  I_LMN(I) = I_CRS1(I)
               ENDDO
               DO I=1,NTERM_LMN
                  J_LMN(I) = J_CRS1(I)
                  LMN(I) =   CRS1(I)
               ENDDO
            ENDIF


            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )             ! II-2 , deallocate CCS1 which was CCS version of GMN

            NTERM_CCS1 = NTERM_CRS1                        ! II-3 , allocate CCS1 to be same as CRS1 but in CCS format
            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFN, NTERM_CCS1, SUBR_NAME )
            CALL SPARSE_CRS_SPARSE_CCS ( NDOFM, NDOFN, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! II-4 , deallocate CRS1
                                                           ! II-5 , sparse multiply to get CRS1 = GMNt*CCS1 with CCS1 = MMM*GMN
!                                                            (note: use SYM_GMN for sym indicator for CCS1)
            CALL MATMULT_SSS_NTERM ( 'GMNt', NDOFN, NTERM_GMN , SYM_GMN, I_GMNt, J_GMNt,                                           &
                                     'CCS1', NDOFN, NTERM_CCS1, SYM_GMN, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1 ',       NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'GMNt', NDOFN, NTERM_GMN , SYM_GMN, I_GMNt, J_GMNt, GMNt,                                           &
                               'CCS1', NDOFN, NTERM_CCS1, SYM_GMN, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )             ! II-6 , deallocate CCS1

                                                           ! II-7 , CRS1 = GMNt*MMM*GMN has all nonzero terms in it.
            IF      (SPARSTOR == 'SYM   ') THEN            !      If SPARSTOR == 'SYM   ', rewrite CRS1 as sym in CRS3

               CALL SPARSE_CRS_TERM_COUNT ( NDOFN, NTERM_CRS1, 'GMNt*MMM*GMN all nonzeros', I_CRS1, J_CRS1, NTERM_CRS3 )
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
               CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS1 = GMNt*MMM*GMN all nonzeros', NDOFN, NTERM_CRS1, I_CRS1, J_CRS1, CRS1,           &
                                            'CRS3 = GMNt*MMM*GMN stored sym'  ,        NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )
               SYM_CRS3 = 'Y'

            ELSE IF (SPARSTOR == 'NONSYM') THEN            !      If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

               NTERM_CRS3 = NTERM_CRS1
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFN, NTERM_CRS3, SUBR_NAME )
               DO I=1,NDOFN+1
                  I_CRS3(I) = I_CRS1(I)
               ENDDO
               DO I=1,NTERM_CRS3
                  J_CRS3(I) = J_CRS1(I)
                    CRS3(I) =   CRS1(I)
               ENDDO
               SYM_CRS3 = 'N'

            ELSE                                           !      Error - incorrect SPARSTOR

               WRITE(ERR,932) SUBR_NAME, SPARSTOR
               WRITE(F06,932) SUBR_NAME, SPARSTOR
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF
                                                           ! II-8 , sparse add to get CRS2 = MNN + CRS3 = MNN + GMNt*MMM*GMN
            CALL MATADD_SSS_NTERM ( NDOFN, 'MNN-bar + MNM*GMN + (MNM*GMN)t', NTERM_MNN, I_MNN, J_MNN, SYM_MNN, 'GMNt*MMM*GMN',     &
                                    NTERM_CRS3, I_CRS3, J_CRS3, SYM_CRS3, 'CRS2', NTERM_CRS2 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFN, NTERM_CRS2, SUBR_NAME )
            CALL MATADD_SSS ( NDOFN, 'MNN-bar + MNM*GMN + (MNM*GMN)t' , NTERM_MNN , I_MNN , J_MNN , MNN, ONE, 'GMNt*MMM*GMN',      &
                              NTERM_CRS3, I_CRS3, J_CRS3, CRS3, ONE, 'CRS2', NTERM_CRS2, I_CRS2, J_CRS2, CRS2 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! II-9 , deallocate CRS1 and CRS3
            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! II-10, deallocate CRS3

            NTERM_MNN = NTERM_CRS2                         ! II-11, reallocate MNN to be size of CRS2
            WRITE(SC1, * ) '    Reallocate MNN'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MNN', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'MNN' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   MNN', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'MNN', NDOFN, NTERM_MNN, SUBR_NAME )

            DO I=1,NDOFN+1                                 ! II-12, set reduced MNN = CRS2 until we see if MNM is null, below
               I_MNN(I) = I_CRS2(I)
            ENDDO
            DO I=1,NTERM_MNN
               J_MNN(I) = J_CRS2(I)
                 MNN(I) =   CRS2(I)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! II-13, Deallocate CRS2

         ELSE                                              ! MMM is null

            CALL DEALLOCATE_SCR_MAT ( 'CCS1')

            IF (NTERM_MMN > 0) THEN                        ! Set LMN - MMN if MMN nonzero, else LMN is null
               NTERM_LMN = NTERM_MMN
               CALL ALLOCATE_SPARSE_MAT ( 'LMN', NDOFM, NTERM_LMN, SUBR_NAME )
               DO I=1,NDOFM+1
                  I_LMN(I) = I_MMN(I)
               ENDDO
               DO I=1,NTERM_LMN
                  J_LMN(I) = J_MMN(I)
                  LMN(I) =   MMN(I)
               ENDDO
            ELSE
               NTERM_LMN = 0
            ENDIF

         ENDIF

         IF (NTERM_LMN > 0) THEN
            CALL WRITE_MATRIX_1 ( LINK2R, L2R, 'Y', 'KEEP', L2R_MSG, 'LMN', NTERM_LMN, NDOFM, I_LMN, J_LMN, LMN )
         ENDIF

         WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GMNt', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'GMNt' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate LMN ', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'LMN' )

! ----------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (MATSPARS == 'N') THEN

         CALL ALLOCATE_FULL_MAT ( 'MNN_FULL', NDOFN, NDOFN, SUBR_NAME )

         IF (NTERM_MNN > 0) THEN
            CALL SPARSE_CRS_TO_FULL ( 'MNN', NTERM_MNN, NDOFN, NDOFN, SYM_MNN, I_MNN, J_MNN,MNN, MNN_FULL )
         ENDIF

         IF (NTERM_MNM > 0) THEN                           ! Part 1: calc MNM*GMN and add it & it's transpose to MNN_FULL

            CALL ALLOCATE_FULL_MAT ( 'MNM_FULL', NDOFN, NDOFM, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'MNM', NTERM_MNM, NDOFN, NDOFM, SYM_MNM, I_MNM, J_MNM, MNM, MNM_FULL )

            CALL ALLOCATE_FULL_MAT ( 'GMN_FULL', NDOFM, NDOFN, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'GMN', NTERM_GMN, NDOFM, NDOFN, SYM_GMN, I_GMN, J_GMN, GMN, GMN_FULL )

            CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFN, NDOFN, SUBR_NAME )

            CALL MATMULT_FFF (  MNM_FULL, GMN_FULL, NDOFN, NDOFM, NDOFN, DUM1 )

            CALL DEALLOCATE_FULL_MAT ( 'MNM_FULL' )
            CALL DEALLOCATE_FULL_MAT ( 'GMN_FULL' )
            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFN, NDOFN, SUBR_NAME )

            ITRNSPB = 0                                    ! Calc MNN-bar + MNM*GMN = MNN + DUM1 and put into DUM1
            CALL MATADD_FFF ( MNN_FULL, DUM1, NDOFN, NDOFN, ALPHA, BETA, ITRNSPB, DUM2 )
            ITRNSPB = 1                                    ! Calc MNN-bar + MNM*GMN + (MNM*GMN)' = DUM2+DUM1' and put into MNN_FULL
            CALL MATADD_FFF ( DUM2    , DUM1, NDOFN, NDOFN, ALPHA, BETA, ITRNSPB, MNN_FULL )


            CALL DEALLOCATE_FULL_MAT ( 'DUM1' )
            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

         ENDIF

         IF (NTERM_MMM > 0) THEN                           ! Part 2: calc GMN(t)*MMM*GMN and add to MNN_FULL

            CALL ALLOCATE_FULL_MAT ( 'MMM_FULL', NDOFM, NDOFM, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'MMM', NTERM_MMM, NDOFM, NDOFM, SYM_MMM, I_MMM, J_MMM, MMM, MMM_FULL )

            CALL ALLOCATE_FULL_MAT ( 'GMN_FULL', NDOFM, NDOFN, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'GMN', NTERM_GMN, NDOFM, NDOFN, SYM_GMN, I_GMN, J_GMN, GMN, GMN_FULL )

            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFM, NDOFN, SUBR_NAME )

            CALL MATMULT_FFF ( MMM_FULL, GMN_FULL, NDOFM, NDOFM, NDOFN, DUM2 )

            CALL DEALLOCATE_FULL_MAT ( 'MMM_FULL' )

            CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFN, NDOFN, SUBR_NAME )

            CALL MATMULT_FFF_T (GMN_FULL, DUM2, NDOFM, NDOFN, NDOFN, DUM1 )

            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

            CALL DEALLOCATE_FULL_MAT ( 'GMN_FULL' )

            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFN, NDOFN, SUBR_NAME )

            ITRNSPB = 0
            CALL MATADD_FFF ( MNN_FULL, DUM1, NDOFN, NDOFN, ALPHA, BETA, ITRNSPB, DUM2 )
            DO I=1,NDOFN
               DO J=1,NDOFN
                  MNN_FULL(I,J) = DUM2(I,J)
               ENDDO
            ENDDO

            CALL DEALLOCATE_FULL_MAT ( 'DUM1' )
            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

         ENDIF

         CALL CNT_NONZ_IN_FULL_MAT ( 'MNN_FULL  ', MNN_FULL, NDOFN, NDOFN, SYM_MNN, NTERM_MNN, SMALL )

         WRITE(SC1, * ) '    Reallocate MNN'
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MNN', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'MNN' )
         WRITE(SC1,12345,ADVANCE='NO') '       Allocate   MNN', CR13
         CALL ALLOCATE_SPARSE_MAT ( 'MNN', NDOFN, NTERM_MNN, SUBR_NAME )

         IF (NTERM_MNN > 0) THEN                            ! Create new sparse arrays from MNN_FULL
!xx         CALL FULL_TO_SPARSE_CRS ( 'MNN_FULL  ', NDOFN, NDOFN, MNN_FULL, NTERM_MNN, SYM_MNN, I_MNN, J_MNN, MNN )
            CALL DEALLOCATE_FULL_MAT ( 'MNN_FULL' )
         ENDIF

      ELSE

         WRITE(ERR,911) SUBR_NAME,MATSPARS
         WRITE(F06,911) SUBR_NAME,MATSPARS
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
  911 FORMAT(' *ERROR   911: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER MATSPARS MUST BE EITHER ','Y',' OR ','N',' BUT VALUE IS ',A)

  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

  936 FORMAT(' *ERROR   936: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF NONZERO TERMS IN SPARSE MATRIX MMN = ',I12,' IS NOT EQUAL TO THOSE IN MATRIX MNM = ',I12)


12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_MGG_TO_MNN


      SUBROUTINE REDUCE_PG_TO_PN ( PART_VEC_G_NM, PART_VEC_SUB )

! Call routines to reduce the PG grid point load matrix from the G-set to the N, M-sets. See Appendix B to the MYSTRAN User's
! Reference Manual for the derivation of the reduction equations.

! NOTE: This subr has code for sparse matrices as well as full matrices, (i.e. Bulk Data PARAM MATSPARS = 'Y' for sparse and 'N'
! for full). The code for full matrices was put in originally in order that the sparse code could be thoroughly checked. That task
! is complete and the remaining full matrix code has not been maintained since around 2005. In addition, new capability added to
! MYSTRAN since that approx time does not have full matrix code.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NDOFN, NDOFM, NSUB, NTERM_GMN, NTERM_PG, NTERM_PN, NTERM_PM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE PARAMS, ONLY                :  EPSIL, MATSPARS
      USE SPARSE_MATRICES, ONLY       :  I_PG, J_PG, PG, I_PN, J_PN, PN, I_PM, J_PM, PM, I_GMN, J_GMN, GMN, I_GMNt, J_GMNt, GMNt
      USE SPARSE_MATRICES, ONLY       :  SYM_GMN, SYM_PG, SYM_PN, SYM_PM
      USE FULL_MATRICES, ONLY         :  PN_FULL, PM_FULL, GMN_FULL, DUM1, DUM2
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM, MATTRNSP_SS
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  SPARSE_CRS_SPARSE_CCS, SPARSE_CRS_TO_FULL
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATADD_FFF, MATMULT_FFF_T
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CNT_NONZ_IN_FULL_MAT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_PG_TO_PN'

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_G_NM(NDOFG)! Partitioning vector (G set into N and M sets)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_SUB(NSUB)  ! Partitioning vector (1's for all subcases)
      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG), PARAMETER        :: ITRNSPB     = 0     ! Transpose indicator for matrix multiply routine
      INTEGER(LONG)                   :: NTERM_CRS1          ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2          ! Number of terms in matrix CRS2
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices
      INTEGER(LONG)                   :: PN_ROW_MAX_TERMS    ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: PM_ROW_MAX_TERMS    ! Output from subr PARTITION_SIZE (max terms in any row of matrix)


      REAL(DOUBLE)                    :: ALPHA = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: BETA  = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: SMALL               ! A number used in filtering out small numbers from a full matrix

      INTRINSIC                       :: DABS




! **********************************************************************************************************************************
! Partition PN from PG (This is PN before reduction, or PN(bar) )

      IF (NDOFN > 0) THEN

         CALL PARTITION_SS_NTERM ( 'PG' , NTERM_PG , NDOFG, NSUB , SYM_PG , I_PG , J_PG ,      PART_VEC_G_NM, PART_VEC_SUB,        &
                                    NUM1, NUM1, PN_ROW_MAX_TERMS, 'PN', NTERM_PN, SYM_PN )

         CALL ALLOCATE_SPARSE_MAT ( 'PN', NDOFN, NTERM_PN, SUBR_NAME )

         IF (NTERM_PN  > 0) THEN
            CALL PARTITION_SS ( 'PG' , NTERM_PG , NDOFG, NSUB , SYM_PG , I_PG , J_PG , PG , PART_VEC_G_NM, PART_VEC_SUB,           &
                                 NUM1, NUM1, PN_ROW_MAX_TERMS, 'PN', NTERM_PN , NDOFN, SYM_PN, I_PN , J_PN , PN  )

         ENDIF

      ENDIF

! Partition PM from PG

      IF (NDOFM > 0) THEN

         CALL PARTITION_SS_NTERM ( 'PG' , NTERM_PG , NDOFG, NSUB , SYM_PG , I_PG , J_PG ,      PART_VEC_G_NM, PART_VEC_SUB,        &
                                    NUM2, NUM1, PM_ROW_MAX_TERMS, 'PM', NTERM_PM, SYM_PM )

         CALL ALLOCATE_SPARSE_MAT ( 'PM', NDOFM, NTERM_PM, SUBR_NAME )

         IF (NTERM_PM  > 0) THEN
            CALL PARTITION_SS ( 'PG' , NTERM_PG , NDOFG, NSUB , SYM_PG , I_PG , J_PG , PG , PART_VEC_G_NM, PART_VEC_SUB,           &
                                 NUM2, NUM1, PM_ROW_MAX_TERMS, 'PM', NTERM_PM , NDOFM, SYM_PM, I_PM , J_PM , PM  )
         ENDIF

      ENDIF

! Reduce PG to PN = PN(bar) + GMN'*PM
! If PARAM MATSPARS = 'Y', then we use sparse matrix operations (multiply/add/transpose). If not, use full matrix operations

      IF (MATSPARS == 'Y') THEN


         IF (NTERM_PM > 0) THEN                            ! Calc GMNt*PM & add it orig PN

            CALL ALLOCATE_SPARSE_MAT ( 'GMNt', NDOFN, NTERM_GMN, SUBR_NAME )
            CALL MATTRNSP_SS ( NDOFM, NDOFN, NTERM_GMN, 'GMN', I_GMN, J_GMN, GMN, 'GMNt', I_GMNt, J_GMNt, GMNt )

                                                           ! CCS1 will be sparse CCS format version of sparse CRS matrix PM
! *********************************************************************************************************************************
! ERROR: NDOFM should be NSUB for col size of PM
!!          CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFM, NTERM_PM, SUBR_NAME )
            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NSUB, NTERM_PM, SUBR_NAME )
! *********************************************************************************************************************************

            CALL SPARSE_CRS_SPARSE_CCS ( NDOFM, NSUB, NTERM_PM, 'PM', I_PM, J_PM, PM, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

                                                           ! Sparse multiply to get CRS1 = GMNt*PM = GMNt*CCS1.
! *********************************************************************************************************************************
! ERROR: NDOFM should be NSUB for col size of PM
!!          CALL MATMULT_SSS_NTERM ( 'GMNt', NDOFN, NTERM_GMN , SYM_GMN, I_GMNt, J_GMNt,                                           &
!!                                   'PM'  , NDOFM, NTERM_PM  , SYM_PM , J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
!!                                   'CRS1' ,       NTERM_CRS1 )

            CALL MATMULT_SSS_NTERM ( 'GMNt', NDOFN, NTERM_GMN , SYM_GMN, I_GMNt, J_GMNt,                                           &
                                     'PM'  , NSUB , NTERM_PM  , SYM_PM , J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1' ,       NTERM_CRS1 )
! *********************************************************************************************************************************

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFN, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'GMNt', NDOFN, NTERM_GMN , SYM_GMN, I_GMNt, J_GMNt, GMNt,                                           &
                               'PM'  , NSUB , NTERM_PM  , SYM_PM , J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

                                                           ! Sparse add to get CRS2 = PN-bar + CRS1
            CALL MATADD_SSS_NTERM ( NDOFN, 'PN-bar', NTERM_PN, I_PN, J_PN, SYM_PN, 'GMNt*PM', NTERM_CRS1, I_CRS1, J_CRS1, 'N',     &
                                         'CRS2', NTERM_CRS2 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFN, NTERM_CRS2, SUBR_NAME )
            CALL MATADD_SSS ( NDOFN, 'PN-bar', NTERM_PN, I_PN, J_PN, PN, ONE, 'GMNt*PM', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE,    &
                                     'CRS2', NTERM_CRS2, I_CRS2, J_CRS2, CRS2 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! Deallocate CRS1

            NTERM_PN = NTERM_CRS2                          ! Reallocate KAA to be size of CRS2
            WRITE(SC1, * ) '    Reallocate PN '
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PN ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'PN' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PN ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'PN', NDOFN, NTERM_PN, SUBR_NAME )

            DO I=1,NDOFN+1
               I_PN(I) = I_CRS2(I)
            ENDDO
            DO J=1,NTERM_PN
               J_PN(J) = J_CRS2(J)
                 PN(J) =   CRS2(J)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! Deallocate CRS2 and GMNt
            WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GMNt', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'GMNt' )

         ENDIF

      ELSE IF (MATSPARS == 'N') THEN

         CALL ALLOCATE_FULL_MAT ( 'PN_FULL', NDOFN, NSUB, SUBR_NAME )  ! Calc reduced PN

         IF (NTERM_PN > 0) THEN
            CALL SPARSE_CRS_TO_FULL ( 'PN', NTERM_PN, NDOFN, NSUB, SYM_PN, I_PN, J_PN, PN, PN_FULL )
         ENDIF

         IF (NTERM_PM > 0) THEN                            ! Calc GMN(t)*PM and add to partitioned PN to get reduced PN

            CALL ALLOCATE_FULL_MAT ( 'PM_FULL', NDOFM, NSUB, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'PM', NTERM_PM, NDOFM, NSUB, SYM_PM, I_PM, J_PM, PM, PM_FULL )

            CALL ALLOCATE_FULL_MAT ( 'GMN_FULL', NDOFM, NDOFN, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'GMN', NTERM_GMN, NDOFM, NDOFN, SYM_GMN, I_GMN, J_GMN, GMN, GMN_FULL )

            CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFN, NSUB, SUBR_NAME )

            CALL MATMULT_FFF_T (GMN_FULL, PM_FULL, NDOFM, NDOFN, NSUB, DUM1 )

            CALL DEALLOCATE_FULL_MAT ( 'PM_FULL' )
            CALL DEALLOCATE_FULL_MAT ( 'GMN_FULL' )

            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFN, NSUB, SUBR_NAME )

            CALL MATADD_FFF ( PN_FULL, DUM1, NDOFN, NSUB, ALPHA, BETA, ITRNSPB, DUM2 )
            DO I=1,NDOFN
               DO J=1,NSUB
                  PN_FULL(I,J) = DUM2(I,J)
               ENDDO
            ENDDO

            CALL DEALLOCATE_FULL_MAT ( 'DUM1' )
            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

         ENDIF


         CALL CNT_NONZ_IN_FULL_MAT ( 'PN_FULL  ', PN_FULL, NDOFN, NSUB, SYM_PN, NTERM_PN, SMALL )

         WRITE(SC1, * ) '    Reallocate PN'
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PN', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'PN' )
         WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PN', CR13
         CALL ALLOCATE_SPARSE_MAT ( 'PN', NDOFN, NTERM_PN, SUBR_NAME )

         IF (NTERM_PN > 0) THEN                            ! Create new sparse arrays from PN_FULL
!xx         CALL FULL_TO_SPARSE_CRS ( 'PN_FULL   ', NDOFN, NSUB , PN_FULL,  NTERM_PN,  SYM_PN,  I_PN,  J_PN,  PN  )
            CALL DEALLOCATE_FULL_MAT ( 'PN_FULL' )
         ENDIF

      ELSE

         WRITE(ERR,911) SUBR_NAME,MATSPARS
         WRITE(F06,911) SUBR_NAME,MATSPARS
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
  911 FORMAT(' *ERROR   911: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER MATSPARS MUST BE EITHER ','Y',' OR ','N',' BUT VALUE IS ',A)


12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_PG_TO_PN



      SUBROUTINE SOLVE_GMN ( PART_VEC_G_NM, PART_VEC_M )

! Solves the sustem of equations: RMM*GMN = -RMN for matrix GMN which is used in the reduction of the G set stiffness, mass and
! load matrices from the G-set to the N, M_sets. If RMM is diagonal, a simple algorithm is used. If it is not, routines
! are called to do the decomp of RMM and the forward-backward substitution (FBS) to obtain GMN

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SCR, L2A, LINK2A, L2A_MSG, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NDOFM, NTERM_RMG, NTERM_RMN, NTERM_RMM, NTERM_GMN
      USE PARAMS, ONLY                :  EPSIL, PRTRMG, PRTGMN, SOLLIB, SPARSE_FLAVOR, SUPINFO
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE SPARSE_MATRICES, ONLY       :  I_RMG, J_RMG, RMG, I_RMN, J_RMN, RMN, I_RMM, J_RMM, RMM, I_GMN, J_GMN, GMN
      USE SPARSE_MATRICES, ONLY       :  SYM_RMG, SYM_RMN, SYM_RMM
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE MATRIX_FILE_IO, ONLY        :  READ_MATRIX_1, READ_MATRIX_2, WRITE_MATRIX_1, WRITE_SPARSE_CRS
      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  SPARSE_CRS_SPARSE_CCS, SPARSE_CRS_TO_FULL
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_GRID_AND_COMP
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, DEALLOCATE_SCR_MAT
      USE SUPERLU_ADAPTERS, ONLY      :  FBS_SUPRLU, SYM_MAT_DECOMP_SUPRLU
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE
      USE SPARSE_CRS_ACCESS, ONLY     :  GET_SPARSE_CRS_COL
      USE SORTING, ONLY               :  SORT_INT2_REAL1
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_GMN'
      CHARACTER(  1*BYTE)             :: CLOSE_IT            ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER(  8*BYTE)             :: CLOSE_STAT          ! Char constant for the CLOSE status of a file
      CHARACTER(  1*BYTE)             :: RMM_DIAG            ! 'Y' if matrix RMM is diagonal.
      CHARACTER(  1*BYTE)             :: RMM_IDENTITY        ! 'Y' if matrix RMM is an identity matrix.

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_G_NM(NDOFG)! Partitioning vector (G set into N and M sets)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_M(NDOFM)   ! Partitioning vector (1's for all M set DOF's)
      INTEGER(LONG)                   :: I,J,K               ! DO loop indices
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices
      INTEGER(LONG)                   :: RMN_ROW_I_NTERMS    ! No. terms in row I of matrix RMN
      INTEGER(LONG)                   :: RMN_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: RMM_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)


      REAL(DOUBLE)                    :: EPS1                ! A small number to compare real zero

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Print out constraint matrix RMG if requested

      IF (( PRTRMG == 1) .OR. ( PRTRMG == 3)) THEN
         IF (NTERM_RMG > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CONSTRAINT MATRIX RMG', 'M ', 'G ', NTERM_RMG, NDOFM, I_RMG, J_RMG, RMG )
         ENDIF
      ENDIF

! Partition RMN from RMG. If no terms in RMN, write error and quit

      CALL PARTITION_SS_NTERM ( 'RMG', NTERM_RMG, NDOFM, NDOFG, SYM_RMG, I_RMG, J_RMG,      PART_VEC_M, PART_VEC_G_NM,             &
                                 NUM1, NUM1, RMN_ROW_MAX_TERMS,'RMN', NTERM_RMN, SYM_RMN )

      IF (NTERM_RMN > 0) THEN
         CALL ALLOCATE_SPARSE_MAT ( 'RMN', NDOFM, NTERM_RMN, SUBR_NAME )
         CALL PARTITION_SS ( 'RMG', NTERM_RMG, NDOFM, NDOFG, SYM_RMG, I_RMG, J_RMG, RMG, PART_VEC_M, PART_VEC_G_NM,                &
                              NUM1, NUM1, RMN_ROW_MAX_TERMS, 'RMN', NTERM_RMN, NDOFM, SYM_RMN, I_RMN, J_RMN, RMN )
      ELSE
         WRITE(ERR,2201) NTERM_RMN
         WRITE(F06,2201) NTERM_RMN
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Partition RMM from RMG. If no terms in RMM, write error and quit

      CALL PARTITION_SS_NTERM ( 'RMG', NTERM_RMG, NDOFM, NDOFG, SYM_RMG, I_RMG, J_RMG ,     PART_VEC_M, PART_VEC_G_NM,             &
                                 NUM1, NUM2, RMM_ROW_MAX_TERMS, 'RMM', NTERM_RMM, SYM_RMM )

      IF (NTERM_RMM > 0) THEN
         CALL ALLOCATE_SPARSE_MAT ( 'RMM', NDOFM, NTERM_RMM, SUBR_NAME )
         CALL PARTITION_SS ( 'RMG', NTERM_RMG, NDOFM, NDOFG, SYM_RMG, I_RMG, J_RMG, RMG, PART_VEC_M, PART_VEC_G_NM,                &
                              NUM1, NUM2, RMM_ROW_MAX_TERMS, 'RMM', NTERM_RMM, NDOFM, SYM_RMM, I_RMM, J_RMM, RMM )
      ELSE
         WRITE(ERR,2202) SUBR_NAME,NTERM_RMM
         WRITE(F06,2202) SUBR_NAME,NTERM_RMM
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Print out constraint matrix RMG partitions, if requested

      IF (( PRTRMG == 2) .OR. ( PRTRMG == 3)) THEN
         IF (NTERM_RMN > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CONSTRAINT PARTITION RMN', 'M ', 'N ', NTERM_RMN, NDOFM, I_RMN, J_RMN, RMN )
         ENDIF
         IF (NTERM_RMM > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CONSTRAINT PARTITION RMM', 'M ', 'M ', NTERM_RMM, NDOFM, I_RMM, J_RMM, RMM )
         ENDIF
      ENDIF

! Find out if RMM is a diagonal matrix. Getting sol'n for GMN will then be trivial

      RMM_DIAG     = 'Y'                                   ! Find out if RMM is a diagonal or identity matrix
      RMM_IDENTITY = 'Y'
      IF (NTERM_RMM == NDOFM) THEN                         ! There are as many terms in RMM as rows so maybe diag or identity
         DO I=1,NDOFM
            IF (J_RMM(I) /= I) THEN                        ! The i-th term in RMM is not a diagonal term
               RMM_DIAG     = 'N'
               RMM_IDENTITY = 'N'
               EXIT
            ENDIF
         ENDDO
         IF (RMM_DIAG == 'Y') THEN                         ! If RMM is diagonal, check for 1.0 on diagonal (identity matrix)
            DO I=1,NDOFM
               IF (DABS(RMM(I) - ONE) > EPS1) THEN
                  RMM_IDENTITY = 'N'
               ENDIF
            ENDDO
         ENDIF
      ELSE                                                 ! RMM is not identity or diagonal
         RMM_DIAG     = 'N'
         RMM_IDENTITY = 'N'
      ENDIF
! Now solve for GMN using either simple algorithm (RMM diagonal) or using BANDED or SPARSE equation solver

      IF ((RMM_DIAG == 'Y') .AND. (DEBUG(20) /= 1)) THEN  ! We can do simple inverse of diagonal matrix RMM

         WRITE(ERR,2293)
         IF (SUPINFO == 'N') THEN
            WRITE(F06,2293)
         ENDIF
         NTERM_GMN = NTERM_RMN

         CALL ALLOCATE_L2_GMN_2 ( SUBR_NAME )
         CALL ALLOCATE_SPARSE_MAT ( 'GMN', NDOFM, NTERM_GMN, SUBR_NAME )

         DO I=1,NDOFM+1
            I_GMN(I) = I_RMN(I)
         ENDDO
         K = 0
         DO I=1,NDOFM
            RMN_ROW_I_NTERMS = I_RMN(I+1) - I_RMN(I)
            DO J=1,RMN_ROW_I_NTERMS
               K = K + 1
               J_GMN(K) = J_RMN(K)
               IF (RMM_IDENTITY == 'Y') THEN
                  GMN(K) = -RMN(K)
               ELSE
                  IF (DABS(RMM(I)) > EPS1) THEN
                     J_GMN(K)  = J_RMN(K)
                       GMN(K)  = -RMN(K)/RMM(I)
                  ELSE
                     WRITE(ERR,2203) K
                     WRITE(F06,2203) K
                     FATAL_ERR = FATAL_ERR + 1
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ELSE                                                 ! Either RMM is not diagonal or DEBUG(20) =1 so we will do full sol'n
!                                                            for GMN from eqn RMM*GMN = -RMN
         IF (RMM_DIAG == 'Y') THEN
            WRITE(ERR,2294)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2294)
            ENDIF
         ENDIF
         CALL SOLVE_GMN_SOLVER

      ENDIF

! Check to maks sure GMN has nonzero terma

      IF (NTERM_GMN <= 0 )THEN
         WRITE(ERR,2204) SUBR_NAME,NTERM_GMN
         WRITE(F06,2204) SUBR_NAME,NTERM_GMN
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error (NTERM_GMN <= 0 ), so quit
      ENDIF

! Close GMN

      IF (NTERM_GMN > 0) THEN
         CLOSE_IT   = 'Y'
         CLOSE_STAT = 'KEEP'
         CALL WRITE_MATRIX_1 ( LINK2A, L2A, CLOSE_IT, CLOSE_STAT, L2A_MSG, 'GMN', NTERM_GMN, NDOFM, I_GMN, J_GMN, GMN )
      ENDIF

! Print out constraint matrix GMN, if requested

      IF ( PRTGMN == 1) THEN
         IF (NTERM_GMN > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CONSTRAINT MATRIX GMN', 'M ', 'N ', NTERM_GMN, NDOFM, I_GMN, J_GMN, GMN )
         ENDIF
      ENDIF

! Deallocate partitions of RMG: RMN, RMM. Keep GMN, it is needed in the reduction of KGG, MGG and PG

      WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RMN', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'RMN' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RMM', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'RMM' )

      CALL DEALLOCATE_L2_GMN_2



      RETURN

! **********************************************************************************************************************************
 2201 FORMAT(' *ERROR  2201: THE RMN PARTITION OF CONSTRAINT MATRIX RMG HAS ',I12,' TERMS IN IT.'                                  &
                    ,/,14X,' THE NUMBER MUST BE > 0.'                                                                              &
                    ,/,14X,' THIS PROBABLY MEANS THAT THERE WERE NO N-SET DOFs THAT THE RIGID ELEMS & MPCs WERE DEPENDENT ON')

 2202 FORMAT(' *ERROR  2202: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE RMM PARTITION OF CONSTRAINT MATRIX RMG HAS ',I12,' TERMS IN IT. THE NUMBER MUST BE > 0.')

 2203 FORMAT(' *ERROR  2203: ZERO DIAGONAL IN MATRIX RMM FOR M-SET DOF ',I8)

 2204 FORMAT(' *ERROR  2204: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THERE IS AN M-SET BUT CONSTRAINT MATRIX GMN HAS ',I12,' TERMS IN IT. MUST BE > 0')

 2293 FORMAT(' *INFORMATION: THE RMM CONSTRAINT MATRIX IS DIAGONAL.'                                                               &
                    ,/,14X,' A SIMPLE SOLUTION FOR THE GMN CONSTRAINT MATRIX WILL BE USED AVOIDING CALLING SUBR SOLVE_GMN',/)

 2294 FORMAT(' *INFORMATION: THE RMM CONSTRAINT MATRIX IS DIAGONAL. HOWEVER, SINCE DEBUG(20) = 1'                                  &
                    ,/,14X,' SUBR SOLVE_GMN_SOLVER WILL BE CALLED TO SOLVE FOR THE GMN CONSTRAINT MATRIX',/)

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(A,10X)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE SOLVE_GMN_SOLVER

! Solves RMM x GMN = -RMN for matrix GMN using unsymmetric decomp from LAPACK (SOLLIB BANDED), or (SOLLIB SPARSE) block by block:
! RMM falls apart into independent blocks (each rigid element or MPC couples only its own dependent DOFs), which module
! BLOCK_DIAGONAL_SOLVE factors separately, so that each column of RMN is solved only in the blocks it touches.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  NDOFG, NDOFM, NDOFN, NTERM_GMN, NTERM_RMM, NTERM_RMN, BLNK_SUB_NAM
      USE PARAMS, ONLY                :  EPSIL, SOLLIB, SPARSE_FLAVOR
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE SPARSE_MATRICES, ONLY       :  I_RMN, J_RMN, RMN, I_RMM, J_RMM, RMM, I2_GMN, I_GMN, J_GMN, GMN
      USE SCRATCH_MATRICES, ONLY      :  I_CCS1, J_CCS1, CCS1
      USE FULL_MATRICES, ONLY         :  RMM_FULL
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, SLU_INFO, SLU_SYMMETRIC
      USE BLOCK_DIAGONAL_SOLVE, ONLY  :  BDS_FACTOR, BDS_FREE, BDS_SOLVE

! Interface module not needed for subr's DGETRF and DGETRS. These are "CONTAIN'ed" in module LAPACK_LIN_EQN_DPB, which
! is "USE'd" above

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_GMN_SOLVER'
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: CALLED_SUBR = ' ' ! Name of a called subr (for output error purposes)
      CHARACTER( 1*BYTE)              :: CLOSE_IT          ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER( 8*BYTE)              :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER(24*BYTE)              :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN
      CHARACTER(44*BYTE)              :: MODNAM            ! Name to write to screen to describe module being run
      CHARACTER(24*BYTE)              :: MODNAM1           ! Name to write to screen to describe module being run
      CHARACTER( 1*BYTE)              :: READ_NTERM        ! 'Y' or 'N' Input to subr READ_MATRIX_1
      CHARACTER( 1*BYTE)              :: NULL_COL          ! 'Y' if a col of RMN is null
      CHARACTER( 1*BYTE)              :: OPND              ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SCRFIL            ! File name
      CHARACTER( 1*BYTE)              :: TRANS             ! 'Y' if

      INTEGER(LONG)                   :: COMPV             ! Component number (1-6) of a grid DOF
      INTEGER(LONG)                   :: GRIDV             ! Grid number
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: INFO      = 0     ! Output from factorization routines
      INTEGER(LONG)                   :: IPIV(NDOFM)       ! Pivot indices from factorization of RMM
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: NRHS              ! No. of RHS's in solving (RMM)*(GMN) = -RMN
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN


      REAL(DOUBLE)                    :: BETA              ! Multiple for rhs for use in subr FBS
      REAL(DOUBLE)                    :: DUM_COL(NDOFM)    ! Temp variable used in SuperLU
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: GMN_COL(NDOFM)    ! A column of GMN solved for herein
      REAL(DOUBLE)                    :: RMN_COL(NDOFM)    ! A column of RMN. The solution for GMN_COL is from RMM*GMN_COL = RMN_COL
      REAL(DOUBLE) , ALLOCATABLE      :: RMN_CCS(:)        ! RMN by cols (CCS), so that getting a col of RMN is not a search of
      INTEGER(LONG), ALLOCATABLE      :: J_RMN_CCS(:)      !   all of RMN (NDOFN cols times NTERM_RMN terms)
      INTEGER(LONG), ALLOCATABLE      :: I_RMN_CCS(:)
      INTEGER(LONG), ALLOCATABLE      :: X_ROW(:)          ! Rows of a column of GMN from BDS_SOLVE
      REAL(DOUBLE) , ALLOCATABLE      :: X_VAL(:)          ! and its values
      INTEGER(LONG)                   :: NB, NX            ! Terms in a column of RMN, and in the solution
      INTEGER(LONG)                   :: NUM_BLOCKS        ! Number of independent blocks of RMM
      INTEGER(LONG)                   :: LARGEST           ! Size of the largest block

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

      IF (SOLLIB == 'BANDED  ') THEN
                                                           ! Create full matrix RMM_FULL from sparse RMM
         CALL ALLOCATE_FULL_MAT  ( 'RMM_FULL', NDOFM, NDOFM, SUBR_NAME )
         CALL SPARSE_CRS_TO_FULL ( 'RMM       ', NTERM_RMM, NDOFM, NDOFM, SYM_RMM, I_RMM, J_RMM, RMM, RMM_FULL )

                                                           ! Perform factorization of RMM_FULL matrix.

         CALL OURTIM
         MODNAM = '    Lapack factorization of RMM'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
         CALL DGETRF (NDOFM, NDOFM, RMM_FULL, NDOFM, IPIV, INFO )

         CALLED_SUBR = 'DGETRF'
         IF      (INFO < 0) THEN                           ! LAPACK subr XERBLA should have reported error on an illegal argument
!                                                            in a called LAPACK subr, so we should not have gotten here
            WRITE(ERR,993) SUBR_NAME, CALLED_SUBR
            WRITE(F06,993) SUBR_NAME, CALLED_SUBR
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                         ! Coding error, so quit

         ELSE IF (INFO > 0) THEN                           ! 0 diag in RMM

            CALL GET_GRID_AND_COMP ( 'M ', INFO, GRIDV, COMPV  )

            WRITE(ERR,2501) CALLED_SUBR, SUBR_NAME, INFO
            WRITE(F06,2501) CALLED_SUBR, SUBR_NAME, INFO
            FATAL_ERR = FATAL_ERR + 1
            IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
               WRITE(ERR,25012) GRIDV, COMPV
               WRITE(F06,25012) GRIDV, COMPV
            ENDIF
            CALL OUTA_HERE ( 'Y' )

         ENDIF

      ELSE IF (SOLLIB == 'SPARSE  ') THEN

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            CALL BDS_FACTOR ( NDOFM, NTERM_RMM, I_RMM, J_RMM, RMM, INFO, NUM_BLOCKS, LARGEST )
            IF (INFO > 0) THEN                             ! Zero pivot: RMM is singular
               CALL GET_GRID_AND_COMP ( 'M ', INFO, GRIDV, COMPV  )
               CALLED_SUBR = 'DGETRF'
               WRITE(ERR,2501) CALLED_SUBR, SUBR_NAME, INFO
               WRITE(F06,2501) CALLED_SUBR, SUBR_NAME, INFO
               FATAL_ERR = FATAL_ERR + 1
               IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
                  WRITE(ERR,25012) GRIDV, COMPV
                  WRITE(F06,25012) GRIDV, COMPV
               ENDIF
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            WRITE(F06,2502) NUM_BLOCKS, LARGEST
            ALLOCATE ( X_ROW(MAX(NDOFM,1)), X_VAL(MAX(NDOFM,1)) )

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

! **********************************************************************************************************************************
! Open a scratch file that will be used to write GMN nonzero terms to as we solve for columns of GMN. After all col's
! of GMN have been solved for, and we have a count on NTERM_GMN, we will allocate memory to the GMN arrays and read
! the scratch file values into those arrays. Then, in the calling subroutine, we will write NTERM_GMN, followed by
! GMN row/col/value to a permanent file

      SCRFIL(1:)  = ' '
      SCRFIL(1:9) = 'SCRATCH-991'
      OPEN (SCR(1),STATUS='SCRATCH',FORM='UNFORMATTED',ACTION='READWRITE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SCRFIL, OUNT )
         CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )
         CALL OUTA_HERE ( 'Y' )                            ! Can't open scratch file, so quit
      ENDIF
      REWIND (SCR(1))

! Loop on columns of RMN

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      ALLOCATE ( J_RMN_CCS(NDOFN+1), I_RMN_CCS(MAX(NTERM_RMN,1)), RMN_CCS(MAX(NTERM_RMN,1)) )
      CALL SPARSE_CRS_SPARSE_CCS ( NDOFM, NDOFN, NTERM_RMN, 'RMN', I_RMN, J_RMN, RMN, 'RMN_CCS', J_RMN_CCS, I_RMN_CCS, RMN_CCS,  &
                                   'N' )

      NTERM_GMN = 0
      CALL COUNTER_INIT('      Solve for GMN col ', NDOFN)
      DO J = 1,NDOFN

         IF (SOLLIB == 'SPARSE  ') THEN                    ! Block by block: only the blocks that col J of RMN touches
            NB = J_RMN_CCS(J+1) - J_RMN_CCS(J)
            IF (NB > 0) THEN
               CALL BDS_SOLVE ( NB, I_RMN_CCS(J_RMN_CCS(J):J_RMN_CCS(J+1)-1), -RMN_CCS(J_RMN_CCS(J):J_RMN_CCS(J+1)-1), NX,     &
                                X_ROW, X_VAL )
               DO K=1,NX                                   ! Count NTERM_GMN and write nonzero GMN to scratch file (rows ascending)
                  IF (DABS(X_VAL(K)) > EPS1) THEN
                     NTERM_GMN = NTERM_GMN + 1
                     WRITE(SCR(1)) X_ROW(K),J,X_VAL(K)
                  ENDIF
               ENDDO
            ENDIF
            CALL COUNTER_PROGRESS(J)
            CYCLE
         ENDIF

         !CALL OURTIM
         !MODNAM1 = '      Solve for GMN col '



! Set RMN_COL to the negative of i-th col of array RMN. First, initialize RMN_COL to zero
! Keep track of whether this col is null, so we can avoid FBS if it is.

         NULL_COL = 'Y'
         BETA = -ONE
         IF (J_RMN_CCS(J+1) > J_RMN_CCS(J)) THEN           ! Col J of RMN has terms (same col as GET_SPARSE_CRS_COL would get)
            NULL_COL = 'N'
            DO I=1,NDOFM
               RMN_COL(I) = ZERO
               GMN_COL(I) = ZERO
            ENDDO
            DO K=J_RMN_CCS(J),J_RMN_CCS(J+1)-1
               RMN_COL(I_RMN_CCS(K)) = BETA*RMN_CCS(K)
            ENDDO
         ENDIF

! Calculate GMN_COL via forward/backward substitution. Remember that rhs is -RMN.

         IF (NULL_COL == 'N') THEN                         ! DGETRS will solve for GMN_COL & load it into GMN array

            IF      (SOLLIB == 'BANDED  ') THEN
               TRANS = 'N'
               NRHS = 1
               CALL DGETRS ( TRANS, NDOFM, NRHS ,RMM_FULL, NDOFM, IPIV, RMN_COL, NDOFM, INFO )

               CALLED_SUBR = 'DGETRS'
               IF      (INFO < 0) THEN                     ! LAPACK subr XERBLA should have reported error on an illegal argument
!                                                            in calling a LAPACK subr, so we should not have gotten here
                  WRITE(ERR,993) SUBR_NAME, CALLED_SUBR
                  WRITE(F06,993) SUBR_NAME, CALLED_SUBR
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )

               ENDIF

            ELSE IF (SOLLIB == 'SPARSE  ') THEN
               IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

                  SLU_INFO = 0
                  CALL FBS_SUPRLU ( SUBR_NAME, 'RMM', NDOFM, NTERM_RMM, J_CCS1, I_CCS1, CCS1, J, RMN_COL, SLU_INFO )

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

            DO I=1,NDOFM
               GMN_COL(I) = RMN_COL(I)
            ENDDO

            DO I=1,NDOFM                                   ! Count NTERM_GMN and write nonzero GMN to scratch file
               IF (DABS(GMN_COL(I)) > EPS1) THEN
                  NTERM_GMN = NTERM_GMN + 1
                  WRITE(SCR(1)) I,J,GMN_COL(I)
               ENDIF
            ENDDO
         ENDIF
         CALL COUNTER_PROGRESS(J)
      ENDDO

      WRITE(SC1,*) CR13
      DEALLOCATE ( J_RMN_CCS, I_RMN_CCS, RMN_CCS )

      CALL DEALLOCATE_SCR_MAT ( 'CCS1' )
      CALL DEALLOCATE_FULL_MAT ( 'RMM_FULL' )

      IF (SOLLIB == 'SPARSE  ') THEN                       ! Last, free the factors of the blocks of RMM
         CALL BDS_FREE
         DEALLOCATE ( X_ROW, X_VAL )
      ENDIF

! The GMN data in SCRATCH-991 is written 1 col at a time. We need it to be written for 1 row at a time with rows in numerical order

      CALL ALLOCATE_L2_GMN_2 ( SUBR_NAME )
      CALL ALLOCATE_SPARSE_MAT ( 'GMN', NDOFM, NTERM_GMN, SUBR_NAME )
      REWIND (SCR(1))
      MESSAG = 'SCRATCH: GMN ROW/COL/VAL'
      READ_NTERM = 'N'
      OPND       = 'Y'
      CLOSE_IT   = 'N'
      CLOSE_STAT = 'KEEP    '
      CALL READ_MATRIX_2 (SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG,'GMN',NDOFM, NTERM_GMN, READ_NTERM, I2_GMN, J_GMN, GMN)
      CALL SORT_INT2_REAL1 ( SUBR_NAME, 'I2_GMN, J_GMN, GMN', NTERM_GMN, I2_GMN, J_GMN, GMN )
      REWIND (SCR(1))
      WRITE(SCR(1)) NTERM_GMN
      DO K=1,NTERM_GMN
         WRITE(SCR(1)) I2_GMN(K),J_GMN(K),GMN(K)
      ENDDO

! Reallocate memory to GMN based on the NTERM_GMN counted above and read values from scratch file into GMN arrays

      CALL DEALLOCATE_L2_GMN_2

      WRITE(SC1, * ) '    Reallocate GMN'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GMN', CR13   ;  CALL DEALLOCATE_SPARSE_MAT ( 'GMN' )
      WRITE(SC1,12345,ADVANCE='NO') '       Allocate   GMN', CR13   ;  CALL ALLOCATE_SPARSE_MAT ('GMN', NDOFM, NTERM_GMN, SUBR_NAME)

      REWIND (SCR(1))
      MESSAG = 'SCRATCH: GMN ROW/COL/VAL'
      READ_NTERM = 'Y'
      OPND       = 'Y'
      CLOSE_IT   = 'Y'
      CLOSE_STAT = 'DELETE  '
      CALL READ_MATRIX_1 ( SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG, 'GMN', NTERM_GMN, READ_NTERM, NDOFM,                &
                           I_GMN, J_GMN, GMN)

      CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )



      RETURN

! **********************************************************************************************************************************
 2502 FORMAT(' *INFORMATION: RMM HAS ',I8,' INDEPENDENT BLOCKS (THE LARGEST HAS ',I8,' DOFS), SOLVED BLOCK BY BLOCK FOR GMN')

 2501 FORMAT(' *ERROR  2501: LAPACK SUBROUTINE, ',A8,' CALLED BY SUBROUTINE ',A                                                    &
                    ,/,14X,' HAS DETECTED A ZERO ON THE DIAG IN ROW ',I12,' OF THE TRIANG FACTOR IN THE DECOMP OF MATRIX RMM.')

25012 FORMAT('               THIS CORRESPONDS TO THE ROW & COL IN RMM FOR GRID POINT ',I8,' COMPONENT ',I3,'.'                     &
                    ,/,14X,' TO CORRECT THIS SITUATION, REMOVE THAT COMPONENT FROM REFC IN FIELD 5 OF THE OFFENDING RBE3(s)')

 2092 FORMAT(4X,A44,20X,I2,':',I2,':',I2,'.',I3)

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(A,10X,A)






  993 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' LAPACK SUBR XERBLA SHOULD HAVE REPORTED AN ERROR ON AN ILLEGAL ARGUMENT IN A CALL TO LAPACK SUBR '    &
                    ,/,15X,A,' (OR A SUBR CALLED BY IT) AND THEN ABORTED')

! **********************************************************************************************************************************

      END SUBROUTINE SOLVE_GMN_SOLVER

      END SUBROUTINE SOLVE_GMN


      SUBROUTINE ALLOCATE_L2_GMN_2 ( CALLING_SUBR )

! Allocate some arrays for use in LINK2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPP6
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFM, NTERM_GMN, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE SPARSE_MATRICES, ONLY       :  I2_GMN

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_L2_GMN_2'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format
      CHARACTER(24*BYTE)              :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: NROWS             ! Number of rows in array
      INTEGER(LONG), PARAMETER        :: NCOLS     = 1     ! Number of cols in array


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
! Allocate array for I2_GMN
! NOTE: I2_GMN has row no's for ALL terms in GMN - not the sparse format that is compressed row storage.
! This is a temporary allocate for subr SOLVE_GMN

      MB_ALLOCATED = ZERO
      NROWS = NTERM_GMN
      JERR = 0

! Allocate array I2_GMN

      NAME = 'I2_GMN                  '
      IF (ALLOCATED(I2_GMN)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (I2_GMN(NTERM_GMN),STAT=IERR)
         IF (IERR == 0) THEN
            DO I=1,NTERM_GMN
               I2_GMN(I) = 0
            ENDDO
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1699) TRIM(SUBR_NAME), CALLING_SUBR
         WRITE(F06,1699) SUBR_NAME, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS)/ONEPP6
      CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )


      RETURN

! **********************************************************************************************************************************
  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_L2_GMN_2


      SUBROUTINE DEALLOCATE_L2_GMN_2

! Deallocate some arrays used in LINK

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE SPARSE_MATRICES, ONLY       :  I2_GMN

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_L2_GMN_2'
      CHARACTER(24*BYTE)              :: NAME              ! Array name (used for output error message)
      CHARACTER(14*BYTE)              :: NAMEL             ! First 14 bytes of NAMEO

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

! Deallocate I2_GMN

      IF (ALLOCATED(I2_GMN)) THEN
         DEALLOCATE (I2_GMN,STAT=IERR)
         NAME = 'I2_GMN                  '
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
         ELSE
            CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         ENDIF
      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)

! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_L2_GMN_2

   END MODULE REDUCTION_G_TO_N
