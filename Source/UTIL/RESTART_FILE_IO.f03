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

   MODULE RESTART_FILE_IO

   USE COL_VEC_LIFECYCLE, ONLY :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
   USE SPARSE_MATRIX_ALLOCATION, ONLY :  ALLOCATE_SPARSE_MAT
   USE SPARSE_MATRIX_DEALLOCATION, ONLY :  DEALLOCATE_SPARSE_MAT
   USE FILE_LIFECYCLE, ONLY :  FILE_CLOSE, FILE_OPEN, READERR
   USE MATRIX_FILE_IO, ONLY :  READ_MATRIX_1, WRITE_SPARSE_CRS
   USE MATRIX_TEXT_OUTPUT, ONLY :  WRITE_VECTOR
   USE DIAGNOSTICS_MEMORY_REPORTING, ONLY :  GET_MATRIX_DIAG_STATS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: PRT_MATS_ON_RESTART, READ_L5A_UG_FOR_SUBCASE, RESTART_DATA_FOR_L3

   CONTAINS

      SUBROUTINE PRT_MATS_ON_RESTART

! Upon user request (via Bulk Data PARAM PRTijk entries) writes matrices to F06 on a restart

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      USE IOUNT1, ONLY                :  ERR, F06, SC1

      USE IOUNT1, ONLY                :  L1E    , L1H    , L1J    , L1L    , L1R    , L2A    , L2B    , L2C    , L2D    , L2E    , &
                                         L2F    , L2G    , L2H    , L2I    , L2J    , L2K    , L2L    , L2M    , L2N    , L2O    , &
                                         L2P    , L2Q    , L3A    , L5A

      USE IOUNT1, ONLY                :  LINK1E , LINK1H , LINK1J , LINK1L , LINK1R , LINK2A , LINK2B , LINK2C , LINK2D , LINK2E , &
                                         LINK2F , LINK2G , LINK2H , LINK2I , LINK2J , LINK2K , LINK2L , LINK2M , LINK2N , LINK2O , &
                                         LINK2P , LINK2Q , LINK3A , LINK5A

      USE IOUNT1, ONLY                :  L1ESTAT, L1HSTAT, L1JSTAT, L1LSTAT, L1RSTAT, L2ASTAT, L2BSTAT, L2CSTAT, L2DSTAT, L2ESTAT, &
                                         L2FSTAT, L2GSTAT, L2HSTAT, L2ISTAT, L2JSTAT, L2KSTAT, L2LSTAT, L2MSTAT, L2NSTAT, L2OSTAT, &
                                         L2PSTAT, L2QSTAT, L3ASTAT, L5ASTAT

      USE IOUNT1, ONLY                :  L1E_MSG, L1H_MSG, L1J_MSG, L1L_MSG, L1R_MSG, L2A_MSG, L2B_MSG, L2C_MSG, L2D_MSG, L2E_MSG, &
                                         L2F_MSG, L2G_MSG, L2H_MSG, L2I_MSG, L2J_MSG, L2K_MSG, L2L_MSG, L2M_MSG, L2N_MSG, L2O_MSG, &
                                         L2P_MSG, L2Q_MSG, L3A_MSG, L5A_MSG

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NSUB, NVEC, SOL_NAME, WARN_ERR

      USE SCONTR, ONLY                :  NDOFA    , NDOFF    , NDOFG    , NDOFL    , NDOFM    , NDOFO    , NDOFR    , NDOFS

      USE SCONTR, ONLY                :  NTERM_KFS , NTERM_GMN , NTERM_GOA , NTERM_HMN , NTERM_HMN , NTERM_KAA , NTERM_KGG ,       &
                                         NTERM_KLL,  NTERM_KRL , NTERM_KRR , NTERM_MAA , NTERM_MGG , NTERM_MLL , NTERM_MRL ,       &
                                         NTERM_MRR , NTERM_PA  , NTERM_PG  , NTERM_PL  , NTERM_PS  , NTERM_QSYS, NTERM_RMG

      USE TIMDAT, ONLY                :  STIME, TSEC

      USE PARAMS, ONLY                :  PRTDISP, PRTFOR, PRTGMN, PRTGOA, PRTHMN, PRTMASS, PRTQSYS, PRTRMG,                        &
                                         PRTSTIFD, PRTSTIFF, PRTUO0, PRTYS, SUPWARN

      USE COL_VECS, ONLY              :  UG_COL, UL_COL, UO0_COL, YSe

      USE SPARSE_MATRICES, ONLY       :  I_GMN , J_GMN , GMN ,I_GOA , J_GOA , GOA ,I_HMN , J_HMN , HMN ,                           &
                                         I_KAA , J_KAA , KAA ,I_KGG , J_KGG , KGG ,I_KLL , J_KLL , KLL ,                           &
                                         I_KRL , J_KRL , KRL ,I_KRR , J_KRR , KRR ,I_KSF , J_KSF , KSF ,                           &
                                         I_MAA , J_MAA , MAA ,I_MGG , J_MGG , MGG ,I_MLL , J_MLL , MLL ,                           &
                                         I_MRL , J_MRL , MRL ,I_MRR , J_MRR , MRR ,                                                &
                                         I_PA  , J_PA  , PA  ,I_PG  , J_PG  , PG  ,I_PL  , J_PL  , PL  ,I_PS  , J_PS  , PS  ,      &
                                         I_QSYS, J_QSYS, QSYS,I_RMG , J_RMG , RMG



      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY: READ_MATRIX_1
      USE MATRIX_FILE_IO, ONLY: WRITE_SPARSE_CRS
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE FILE_LIFECYCLE, ONLY: FILE_OPEN
      USE COL_VEC_LIFECYCLE, ONLY     :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
      USE FILE_LIFECYCLE, ONLY: READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE MATRIX_TEXT_OUTPUT, ONLY: WRITE_VECTOR
      USE FILE_LIFECYCLE, ONLY: FILE_CLOSE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_MATRIX_DIAG_STATS

      IMPLICIT NONE

      LOGICAL                         :: FILE_EXIST        ! Result from INQUIRE is true if a file exists

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME  = 'PRT_MATS_ON_RESTART'
      CHARACTER(  1*BYTE)             :: CLOSE_IT   = 'Y'  ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER(  1*BYTE)             :: NTERM_RD   = 'Y'  ! 'Y' or 'N' Input to subr READ_MATRIX_1
      CHARACTER(  1*BYTE)             :: OPND       = 'N'  ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERROR            ! Error count
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: NTERM_KSF         ! Number of nonzeros in sparse matrix KSF (= NTERM_KFS)
      INTEGER(LONG)                   :: NUM_SOLNS         ! NSUB for statics, NVEC for eigenvalues, etc
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file


      REAL(DOUBLE)                    :: KAA_DIAG(NDOFA)   ! Diagonal of KAA
      REAL(DOUBLE)                    :: KGG_DIAG(NDOFG)   ! Diagonal of KGG
      REAL(DOUBLE)                    :: KLL_DIAG(NDOFL)   ! Diagonal of KLL
      REAL(DOUBLE)                    :: KRR_DIAG(NDOFR)   ! Diagonal of KRR

      REAL(DOUBLE)                    :: KAA_MAX_DIAG      ! Max diag term from KAA
      REAL(DOUBLE)                    :: KGG_MAX_DIAG      ! Max diag term from KGG
      REAL(DOUBLE)                    :: KLL_MAX_DIAG      ! Max diag term from KLL
      REAL(DOUBLE)                    :: KRR_MAX_DIAG      ! Max diag term from KRR



! **********************************************************************************************************************************
      IF      ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
         NUM_SOLNS = NSUB
      ELSE IF (SOL_NAME == 'MODES           ') THEN
         NUM_SOLNS = NVEC
      ELSE IF (SOL_NAME == 'GEN CB MODEL    ') THEN
         NUM_SOLNS = 2*NDOFR + NVEC
      ELSE
      ENDIF

      OUNT(1) = ERR
      OUNT(2) = F06

!xx   WRITE(SC1, * ) '    ALLOCATE/DEALLOCATE SOME ARRAYS'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

! ----------------------------------------------------------------------------------------------------------------------------------
! Write PG load matrix

      IF (PRTFOR(1) >= 1) THEN
         INQUIRE ( FILE=LINK1E, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_PG > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PG  ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'PG', NDOFG, NTERM_PG, SUBR_NAME)
            CALL READ_MATRIX_1 ( LINK1E, L1E, OPND, CLOSE_IT, L1ESTAT, L1E_MSG,'PG', NTERM_PG, NTERM_RD, NDOFG, I_PG, J_PG, PG)
            CALL WRITE_SPARSE_CRS ( L1E_MSG, 'G ', 'SUBCASE', NTERM_PG, NDOFG, I_PG, J_PG, PG )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PG  ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'PG' )
            WRITE(SC1,*) CR13
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L1E_MSG, LINK1E
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L1E_MSG, LINK1E
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write YSe displ vector

      IF (PRTYS == 1) THEN

         INQUIRE ( FILE=LINK1H, EXIST=FILE_EXIST )
         IF (FILE_EXIST) THEN

            CALL FILE_OPEN ( L1H, LINK1H, OUNT, 'OLD', L1H_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
            CALL ALLOCATE_COL_VEC ( 'YSe', NDOFS, SUBR_NAME )

            REC_NO = 0
            IERROR = 0
            DO I=1,NDOFS
               REC_NO = REC_NO + 1
               READ(L1H,IOSTAT=IOCHK) YSe(I)
               IF (IOCHK /= 0) THEN
                  CALL READERR ( IOCHK, LINK1H, L1H_MSG, REC_NO, OUNT )
                  IERROR = IERROR + 1
               ENDIF
            ENDDO

            IF (IERROR /= 0) THEN
               WRITE(ERR,9995) IERROR
               WRITE(F06,9995) IERROR
               CALL OUTA_HERE ( 'Y' )
            ENDIF

            CALL WRITE_VECTOR ( 'S-SET ENFORCED DISPL VECTOR', 'DISPL', NDOFS, YSe)
            CALL FILE_CLOSE ( L1H, LINK1H, L1HSTAT )
            CALL DEALLOCATE_COL_VEC ( 'YSe' )

         ELSE

            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L1H_MSG, LINK1H
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L1H_MSG, LINK1H
            ENDIF

         ENDIF

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write RMG

      IF (PRTRMG > 0) THEN
         INQUIRE ( FILE=LINK1J, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_RMG > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   RMG ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'RMG', NDOFM, NTERM_RMG, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK1J, L1J, OPND, CLOSE_IT, L1JSTAT, L1J_MSG,'RMG', NTERM_RMG, NTERM_RD, NDOFM, I_RMG, J_RMG, RMG)
            CALL WRITE_SPARSE_CRS ( L1J_MSG, 'M ', 'G ', NTERM_RMG, NDOFM, I_RMG, J_RMG, RMG )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RMG ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'RMG' )
            WRITE(SC1,*) CR13
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L1J_MSG, LINK1J
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L1J_MSG, LINK1J
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write KGG stiffness matrix and matrix diagonal and stats, if requested
! NOTE: KGG was allocated and read in LINK0, if needed, prior to this (if RESTART = 'Y')
! ----

      IF ((PRTSTIFF(1) >= 1) .OR. (PRTSTIFD(1) >= 1)) THEN
         INQUIRE ( FILE=LINK1L, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_KGG > 0)) THEN
!xx         CALL ALLOCATE_SPARSE_MAT ( 'KGG', NDOFG, NTERM_KGG, SUBR_NAME )
!xx         CALL READ_MATRIX_1 ( LINK1L, L1L, OPND, CLOSE_IT, L1LSTAT, L1L_MSG,'KGG', NTERM_KGG, NTERM_RD, NDOFG, I_KGG, J_KGG, KGG)
            IF (PRTSTIFF(1) >= 1) THEN
               CALL WRITE_SPARSE_CRS ( L1L_MSG, 'G ', 'G ', NTERM_KGG, NDOFG, I_KGG, J_KGG, KGG )
            ENDIF
            IF (PRTSTIFD(1) > 0) THEN
               CALL GET_MATRIX_DIAG_STATS ( 'KGG', 'G ', NDOFG, NTERM_KGG, I_KGG, J_KGG, KGG, PRTSTIFD(1), KGG_DIAG, KGG_MAX_DIAG  )
            ENDIF
!xx         CALL DEALLOCATE_SPARSE_MAT ( 'KGG' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L1L_MSG, LINK1L
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L1L_MSG, LINK1L
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write MGG mass matrix

      IF (PRTMASS(1) >= 1) THEN
         INQUIRE ( FILE=LINK1R, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_MGG > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   MGG ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'MGG', NDOFG, NTERM_MGG, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK1R, L1R, OPND, CLOSE_IT, L1RSTAT, L1R_MSG,'MGG', NTERM_MGG, NTERM_RD, NDOFG, I_MGG, J_MGG, MGG)
            CALL WRITE_SPARSE_CRS ( L1R_MSG, 'G ', 'G ', NTERM_MGG, NDOFG, I_MGG, J_MGG, MGG )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MGG ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'MGG' )
            WRITE(SC1,*) CR13
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L1R_MSG, LINK1R
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L1R_MSG, LINK1R
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write GMN

      IF (PRTGMN > 0) THEN
         INQUIRE ( FILE=LINK2A, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_GMN > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   GMN ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'GMN', NDOFM, NTERM_GMN, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2A, L2A, OPND, CLOSE_IT, L2ASTAT, L2A_MSG,'GMN', NTERM_GMN, NTERM_RD, NDOFM, I_GMN, J_GMN, GMN)
            CALL WRITE_SPARSE_CRS ( L2A_MSG, 'M ', 'N ', NTERM_GMN, NDOFM, I_GMN, J_GMN, GMN )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GMN ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'GMN' )
            WRITE(SC1,*) CR13
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2A_MSG, LINK2A
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2A_MSG, LINK2A
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write KSF stiffness matrix

      IF ((PRTSTIFF(3) == 2) .OR. (PRTSTIFF(3) == 3)) THEN
         INQUIRE ( FILE=LINK2B, EXIST=FILE_EXIST )
         NTERM_KSF = NTERM_KFS
         IF ((FILE_EXIST) .AND. (NTERM_KSF > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KSF ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KSF', NDOFS, NTERM_KSF, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2B, L2B, OPND, CLOSE_IT, L2BSTAT, L2B_MSG,'KSF', NTERM_KSF, NTERM_RD, NDOFS, I_KSF, J_KSF, KSF)
            CALL WRITE_SPARSE_CRS ( L2B_MSG, 'S ', 'F ', NTERM_KSF, NDOFS, I_KSF, J_KSF, KSF )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KSF ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KSF' )
            WRITE(SC1,*) CR13
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2B_MSG, LINK2B
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2B_MSG, LINK2B
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write QSYS matrix

      IF (PRTQSYS > 0) THEN
         INQUIRE ( FILE=LINK2C, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_QSYS > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   QSYS', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'QSYS', NDOFS, NTERM_QSYS, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2C, L2C, OPND, CLOSE_IT, L2CSTAT, L2C_MSG,'QSYS', NTERM_QSYS, NTERM_RD, NDOFS,                &
                                 I_QSYS, J_QSYS, QSYS)
            CALL WRITE_SPARSE_CRS ( L2C_MSG, 'S ', '  ', NTERM_QSYS, NDOFS, I_QSYS, J_QSYS, QSYS )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate QSYS', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'QSYS' )
            WRITE(SC1,*) CR13
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2C_MSG, LINK2C
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2C_MSG, LINK2C
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write PS load matrix

      IF ((PRTFOR(3) == 2) .OR. (PRTFOR(3) == 3)) THEN
         INQUIRE ( FILE=LINK2D, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_PS > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PS  ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'PS', NDOFS, NTERM_PS, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2D, L2D, OPND, CLOSE_IT, L2DSTAT, L2D_MSG,'PS', NTERM_PS, NTERM_RD, NDOFS, I_PS, J_PS, PS)
            CALL WRITE_SPARSE_CRS ( L2D_MSG, 'S ', 'SUBCASE', NTERM_PS, NDOFS, I_PS, J_PS, PS )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PS  ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'PS' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2D_MSG, LINK2D
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2D_MSG, LINK2D
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write GOA

      IF (PRTGOA > 0) THEN
         INQUIRE ( FILE=LINK2E, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_GOA > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   GOA ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'GOA', NDOFO, NTERM_GOA, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2E, L2E, OPND, CLOSE_IT, L2ESTAT, L2E_MSG,'GOA', NTERM_GOA, NTERM_RD, NDOFO, I_GOA, J_GOA, GOA)
            CALL WRITE_SPARSE_CRS ( L2E_MSG, 'O ', 'A ', NTERM_GOA, NDOFO, I_GOA, J_GOA, GOA )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GOA ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'GOA' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2E_MSG, LINK2E
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2E_MSG, LINK2E
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write UO0 displ vector

      IF (PRTUO0 == 1) THEN

         INQUIRE ( FILE=LINK2F, EXIST=FILE_EXIST )
         IF (FILE_EXIST) THEN
            CALL FILE_OPEN ( L2F, LINK2F, OUNT, 'OLD', L2F_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
            DO J=1,NUM_SOLNS

               CALL ALLOCATE_COL_VEC ( 'UO0_COL', NDOFO, SUBR_NAME )

               REC_NO = 0
               IERROR = 0
               DO I=1,NDOFO
                  REC_NO = REC_NO + 1
                  READ(L2F,IOSTAT=IOCHK) UO0_COL(I)
                  IF (IOCHK /= 0) THEN
                     CALL READERR ( IOCHK, LINK2F, L2F_MSG, REC_NO, OUNT )
                     IERROR = IERROR + 1
                  ENDIF
               ENDDO

               IF (IERROR /= 0) THEN
                  WRITE(ERR,9995) IERROR
                  WRITE(F06,9995) IERROR
                  CALL OUTA_HERE ( 'Y' )
               ENDIF

               CALL WRITE_VECTOR ( 'UO0 = KOO(-1)*PO', 'DISPL', NDOFO, UO0_COL)
               CALL DEALLOCATE_COL_VEC ( 'UO0_COL' )

            ENDDO
            CALL FILE_CLOSE ( L2F, LINK2F, L2FSTAT )

         ELSE

            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2F_MSG, LINK2F
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2F_MSG, LINK2F
            ENDIF

         ENDIF

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write KLL stiffness matrix and matrix diagonal and stats, if requested

      IF ((PRTSTIFF(5) == 1) .OR. (PRTSTIFF(5) == 3) .OR.(PRTSTIFD(5) > 0)) THEN
         INQUIRE ( FILE=LINK2G, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_KLL > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KLL ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KLL', NDOFL, NTERM_KLL, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2G, L2G, OPND, CLOSE_IT, L2GSTAT, L2G_MSG,'KLL', NTERM_KLL, NTERM_RD, NDOFL, I_KLL, J_KLL, KLL)
            IF ((PRTSTIFF(5) == 1) .OR. (PRTSTIFF(5) == 3)) THEN
               CALL WRITE_SPARSE_CRS ( L2G_MSG, 'L ', 'L ', NTERM_KLL, NDOFL, I_KLL, J_KLL, KLL )
            ENDIF
            IF (PRTSTIFD(5) > 0) THEN
               CALL GET_MATRIX_DIAG_STATS ( 'KLL', 'L ', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, PRTSTIFD(5), KLL_DIAG, KLL_MAX_DIAG  )
            ENDIF
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLL ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KLL' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2G_MSG, LINK2G
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2G_MSG, LINK2G
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write PL load matrix

      IF ((PRTFOR(5) == 1) .OR. (PRTFOR(5) == 3)) THEN
         INQUIRE ( FILE=LINK2H, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_PL > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PL  ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'PL', NDOFL, NTERM_PL, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2H, L2H, OPND, CLOSE_IT, L2HSTAT, L2H_MSG,'PL', NTERM_PL, NTERM_RD, NDOFL, I_PL, J_PL, PL)
            CALL WRITE_SPARSE_CRS ( L2H_MSG, 'L ', 'SUBCASE', NTERM_PL, NDOFL, I_PL, J_PL, PL )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PL  ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'PL' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2H_MSG, LINK2H
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2H_MSG, LINK2H
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write MLL mass matrix

      IF ((PRTMASS(5) == 1) .OR. (PRTMASS(5) == 3)) THEN
         INQUIRE ( FILE=LINK2I, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_MLL > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   MLL ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'MLL', NDOFL, NTERM_MLL, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2I, L2I, OPND, CLOSE_IT, L2ISTAT, L2I_MSG,'MLL', NTERM_MLL, NTERM_RD, NDOFL, I_MLL, J_MLL, MLL)
            CALL WRITE_SPARSE_CRS ( L2I_MSG, 'L ', 'L ', NTERM_MLL, NDOFL, I_MLL, J_MLL, MLL )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MLL ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'MLL' )

            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2I_MSG, LINK2I
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2I_MSG, LINK2I
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write HMN

      IF (PRTHMN > 0) THEN
         INQUIRE ( FILE=LINK2J, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_HMN > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   HMN ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'HMN', NDOFM, NTERM_HMN, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2J, L2J, OPND, CLOSE_IT, L2JSTAT, L2J_MSG,'HMN', NTERM_HMN, NTERM_RD, NDOFM, I_HMN, J_HMN, HMN)
            CALL WRITE_SPARSE_CRS ( L2J_MSG, 'M ', 'N ', NTERM_HMN, NDOFM, I_HMN, J_HMN, HMN )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate HMN ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'HMN' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2J_MSG, LINK2J
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2J_MSG, LINK2J
            ENDIF
         ENDIF
      ENDIF

! -------------------------------------------------------------------------------
! Write KRL stiffness matrix

      IF ((PRTSTIFF(5) == 2) .OR. (PRTSTIFF(5) == 3)) THEN
         INQUIRE ( FILE=LINK2K, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_KRL > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KRL ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KRL', NDOFS, NTERM_KRL, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2K, L2K, OPND, CLOSE_IT, L2KSTAT, L2K_MSG,'KRL', NTERM_KRL, NTERM_RD, NDOFS, I_KRL, J_KRL, KRL)
            CALL WRITE_SPARSE_CRS ( L2K_MSG, 'R ', 'L ', NTERM_KRL, NDOFS, I_KRL, J_KRL, KRL )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KRL ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KRL' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2K_MSG, LINK2K
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2K_MSG, LINK2K
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write KRR stiffness matrix and matrix diagonal and stats, if requested

      IF ((PRTSTIFF(5) == 2) .OR. (PRTSTIFF(5) == 3) .OR. (PRTSTIFD(5) == 1)) THEN
         INQUIRE ( FILE=LINK2L, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_KRR > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KRR ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KRR', NDOFR, NTERM_KRR, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2L, L2L, OPND, CLOSE_IT, L2LSTAT, L2L_MSG,'KRR', NTERM_KRR, NTERM_RD, NDOFR, I_KRR, J_KRR, KRR)
            IF ((PRTSTIFF(5) == 2) .OR. (PRTSTIFF(5) == 3)) THEN
               CALL WRITE_SPARSE_CRS ( L2L_MSG, 'R ', 'R ', NTERM_KRR, NDOFR, I_KRR, J_KRR, KRR )
            ENDIF
            IF (PRTSTIFD(5) > 0) THEN
               CALL GET_MATRIX_DIAG_STATS ( 'KRR', 'R ', NDOFR, NTERM_KRR, I_KRR, J_KRR, KRR, PRTSTIFD(5), KRR_DIAG, KRR_MAX_DIAG  )
            ENDIF
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KRR ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KRR' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2L_MSG, LINK2L
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2L_MSG, LINK2L
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write MRL mass matrix

      IF ((PRTMASS(5) == 2) .OR. (PRTMASS(5) == 3)) THEN
         INQUIRE ( FILE=LINK2M, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_MRL > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   MRL ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'MRL', NDOFS, NTERM_MRL, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2M, L2M, OPND, CLOSE_IT, L2MSTAT, L2M_MSG,'MRL', NTERM_MRL, NTERM_RD, NDOFS, I_MRL, J_MRL, MRL)
            CALL WRITE_SPARSE_CRS ( L2M_MSG, 'R ', 'L ', NTERM_MRL, NDOFR, I_MRL, J_MRL, MRL )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MRL ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'MRL' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2M_MSG, LINK2M
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2M_MSG, LINK2M
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write MRR mass matrix

      IF ((PRTMASS(5) == 2) .OR. (PRTMASS(5) == 3)) THEN
         INQUIRE ( FILE=LINK2N, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_MRR > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   MRR ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'MRR', NDOFR, NTERM_MRR, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2N, L2N, OPND, CLOSE_IT, L2NSTAT, L2N_MSG,'MRR', NTERM_MRR, NTERM_RD, NDOFR, I_MRR, J_MRR, MRR)
            CALL WRITE_SPARSE_CRS ( L2N_MSG, 'R ', 'R ', NTERM_MRR, NDOFR, I_MRR, J_MRR, MRR )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MRR ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'MRR' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2N_MSG, LINK2N
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2N_MSG, LINK2N
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write KAA stiffness matrix and matrix diagonal and stats, if requested

      IF ((PRTSTIFF(4) == 1) .OR. (PRTSTIFF(4) == 3) .OR. (PRTSTIFD(4) == 1)) THEN
         INQUIRE ( FILE=LINK2O, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_KAA > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KAA ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KAA', NDOFA, NTERM_KAA, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2O, L2O, OPND, CLOSE_IT, L2OSTAT, L2O_MSG,'KAA', NTERM_KAA, NTERM_RD, NDOFA, I_KAA, J_KAA, KAA)
            IF ((PRTSTIFF(4) == 1) .OR. (PRTSTIFF(4) == 3)) THEN
               CALL WRITE_SPARSE_CRS ( L2O_MSG, 'A ', 'A ', NTERM_KAA, NDOFA, I_KAA, J_KAA, KAA )
            ENDIF
            IF (PRTSTIFD(4) > 0) THEN
               CALL GET_MATRIX_DIAG_STATS ( 'KAA', 'A ', NDOFA, NTERM_KAA, I_KAA, J_KAA, KAA, PRTSTIFD(4), KAA_DIAG, KAA_MAX_DIAG  )
            ENDIF
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KAA ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KAA' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2O_MSG, LINK2O
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2O_MSG, LINK2O
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write MAA mass matrix

      IF ((PRTMASS(4) >= 1) .OR. (PRTMASS(4) >= 3)) THEN
         INQUIRE ( FILE=LINK2P, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_MAA > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   MAA ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'MAA', NDOFA, NTERM_MAA, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2P, L2P, OPND, CLOSE_IT, L2PSTAT, L2P_MSG,'MAA', NTERM_MAA, NTERM_RD, NDOFA, I_MAA, J_MAA, MAA)
            CALL WRITE_SPARSE_CRS ( L2P_MSG, 'A ', 'A ', NTERM_MAA, NDOFA, I_MAA, J_MAA, MAA )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MAA ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'MAA' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2P_MSG, LINK2P
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2P_MSG, LINK2P
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write PA load matrix

      IF ((PRTFOR(4) == 1) .OR. (PRTFOR(4) == 3)) THEN
         INQUIRE ( FILE=LINK2Q, EXIST=FILE_EXIST )
         IF ((FILE_EXIST) .AND. (NTERM_PA > 0)) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PA  ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'PA', NDOFA, NTERM_PA, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK2Q, L2Q, OPND, CLOSE_IT, L2QSTAT, L2Q_MSG,'PA', NTERM_PA, NTERM_RD, NDOFA, I_PA, J_PA, PA)
            CALL WRITE_SPARSE_CRS ( L2Q_MSG, 'A ', 'SUBCASE', NTERM_PA, NDOFA, I_PA, J_PA, PA )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PA  ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'PA' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L2Q_MSG, LINK2Q
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L2Q_MSG, LINK2Q
            ENDIF
         ENDIF
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write UL displ vector

      IF (PRTDISP(5) == 1) THEN

         INQUIRE ( FILE=LINK3A, EXIST=FILE_EXIST )
         IF (FILE_EXIST) THEN
            CALL FILE_OPEN ( L3A, LINK3A, OUNT, 'OLD', L3A_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
            DO J=1,NUM_SOLNS

               CALL ALLOCATE_COL_VEC ( 'UL_COL', NDOFL, SUBR_NAME )

               REC_NO = 0
               IERROR = 0
               DO I=1,NDOFL
                  REC_NO = REC_NO + 1
                  READ(L3A,IOSTAT=IOCHK) UL_COL(I)
                  IF (IOCHK /= 0) THEN
                     CALL READERR ( IOCHK, LINK3A, L3A_MSG, REC_NO, OUNT )
                     IERROR = IERROR + 1
                  ENDIF
               ENDDO

               IF (IERROR /= 0) THEN
                  WRITE(ERR,9995) IERROR
                  WRITE(F06,9995) IERROR
                  CALL OUTA_HERE ( 'Y' )
               ENDIF

               CALL WRITE_VECTOR ( 'L-SET DISPL VECTOR', 'DISPL', NDOFL, UL_COL)
               CALL DEALLOCATE_COL_VEC ( 'UL_COL' )

            ENDDO
            CALL FILE_CLOSE ( L3A, LINK3A, L3ASTAT )

         ELSE

            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L3A_MSG, LINK3A
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L3A_MSG, LINK3A
            ENDIF

         ENDIF

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Write UG displ vector

      IF (PRTDISP(1) == 1) THEN

         INQUIRE ( FILE=LINK5A, EXIST=FILE_EXIST )
         IF (FILE_EXIST) THEN
            CALL FILE_OPEN ( L5A, LINK5A, OUNT, 'OLD', L5A_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
            DO J=1,NUM_SOLNS

               CALL ALLOCATE_COL_VEC ( 'UG_COL', NDOFG, SUBR_NAME )

               REC_NO = 0
               IERROR = 0
               DO I=1,NDOFG
                  REC_NO = REC_NO + 1
                  READ(L5A,IOSTAT=IOCHK) UG_COL(I)
                  IF (IOCHK /= 0) THEN
                     CALL READERR ( IOCHK, LINK5A, L5A_MSG, REC_NO, OUNT )
                     IERROR = IERROR + 1
                  ENDIF
               ENDDO

               IF (IERROR /= 0) THEN
                  WRITE(ERR,9995) IERROR
                  WRITE(F06,9995) IERROR
                  CALL OUTA_HERE ( 'Y' )
               ENDIF

               CALL WRITE_VECTOR ( 'G-SET DISPL VECTOR', 'DISPL', NDOFG, UG_COL)
               CALL DEALLOCATE_COL_VEC ( 'UG_COL' )

            ENDDO
            CALL FILE_CLOSE ( L5A, LINK5A, L5ASTAT )

         ELSE

            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,101) L5A_MSG, LINK5A
            IF (SUPWARN == 'N') THEN
               WRITE(F06,101) L5A_MSG, LINK5A
            ENDIF

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *WARNING    : THE FOLLOWING FILE FOR THE ',A                                                                        &
                    ,/,14X,' EITHER DOES NOT EXIST OR IS NULL: '                                                                   &
                    ,/,15X,A)

 9995 FORMAT(/,' PROCESSING ENDED DUE TO ABOVE ',I8,' ERRORS')

12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE PRT_MATS_ON_RESTART


      SUBROUTINE READ_L5A_UG_FOR_SUBCASE ( ISUB, IERROR )

! Reads the G-set displacement vector (UG_COL) for a specific internal subcase index from file LINK5A.
! File LINK5A is written by LINK5 during the linear-static portion of an analysis: for each of NUM_SOLNS solution columns,
! LINK5 writes NDOFG unformatted records, one per G-set DOF, holding the column's UG values. For SOL 105 (BUCKLING) the
! step-1 LINK5 invocation writes one full UG vector per subcase (NUM_SOLNS == NSUB), so seeking past (ISUB-1)*NDOFG records
! and reading the next NDOFG records yields subcase ISUB's UG displacements.
!
! Caller responsibilities:
!   * UG_COL must already be allocated to size NDOFG (use ALLOCATE_COL_VEC('UG_COL',NDOFG,SUBR_NAME) if not).
!   * The L5A file must contain at least ISUB columns of NDOFG records each (i.e. LINK5 step 1 must have run with
!     NUM_SOLNS >= ISUB).
! On error, IERROR is incremented; on success it is left unchanged.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L5A, LINK5A, L5A_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG
      USE COL_VECS, ONLY              :  UG_COL

      USE FILE_LIFECYCLE, ONLY: FILE_OPEN
      USE FILE_LIFECYCLE, ONLY: FILE_CLOSE
      USE FILE_LIFECYCLE, ONLY: READERR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'READ_L5A_UG_FOR_SUBCASE'

      INTEGER(LONG), INTENT(IN)       :: ISUB              ! Internal subcase index (1..NSUB) to fetch UG_COL for
      INTEGER(LONG), INTENT(INOUT)    :: IERROR            ! Cumulative error count, incremented on any IO failure

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT from READ
      INTEGER(LONG)                   :: REC_NO            ! Record number for READERR diagnostics
      INTEGER(LONG)                   :: SKIP              ! Number of records to skip before reaching subcase ISUB
      INTEGER(LONG)                   :: OUNT(2)           ! File units for READERR



! **********************************************************************************************************************************
      OUNT(1) = ERR
      OUNT(2) = F06

      IF (ISUB < 1) THEN
         WRITE(ERR,9101) SUBR_NAME, ISUB
         WRITE(F06,9101) SUBR_NAME, ISUB
         FATAL_ERR = FATAL_ERR + 1
         IERROR = IERROR + 1
         RETURN
      ENDIF

      IF (.NOT. ALLOCATED(UG_COL)) THEN
         WRITE(ERR,9102) SUBR_NAME
         WRITE(F06,9102) SUBR_NAME
         FATAL_ERR = FATAL_ERR + 1
         IERROR = IERROR + 1
         RETURN
      ENDIF

      IF (SIZE(UG_COL) < NDOFG) THEN
         WRITE(ERR,9103) SUBR_NAME, SIZE(UG_COL), NDOFG
         WRITE(F06,9103) SUBR_NAME, SIZE(UG_COL), NDOFG
         FATAL_ERR = FATAL_ERR + 1
         IERROR = IERROR + 1
         RETURN
      ENDIF

      ! Open L5A fresh in READ/REWIND mode. If it is already open elsewhere we close it first so we can rewind cleanly.
      CALL FILE_CLOSE ( L5A, LINK5A, 'KEEP' )
      CALL FILE_OPEN  ( L5A, LINK5A, OUNT, 'OLD', L5A_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

      ! Skip (ISUB-1)*NDOFG records to position at the start of subcase ISUB's UG vector.
      SKIP = (ISUB - 1) * NDOFG
      REC_NO = 0
      DO I = 1, SKIP
         REC_NO = REC_NO + 1
         READ(L5A,IOSTAT=IOCHK)
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK5A, L5A_MSG, REC_NO, OUNT )
            IERROR = IERROR + 1
            CALL FILE_CLOSE ( L5A, LINK5A, 'KEEP' )
            RETURN
         ENDIF
      ENDDO

      ! Read NDOFG records into UG_COL.
      DO I = 1, NDOFG
         REC_NO = REC_NO + 1
         READ(L5A,IOSTAT=IOCHK) UG_COL(I)
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK5A, L5A_MSG, REC_NO, OUNT )
            IERROR = IERROR + 1
            CALL FILE_CLOSE ( L5A, LINK5A, 'KEEP' )
            RETURN
         ENDIF
      ENDDO

      ! Leave L5A closed so the next writer (e.g. step-2 LINK5) can REPLACE it cleanly.
      CALL FILE_CLOSE ( L5A, LINK5A, 'KEEP' )



      RETURN

! **********************************************************************************************************************************
 9101 FORMAT(' *ERROR  9101: PROGRAMMING ERROR IN ',A,': ISUB MUST BE >= 1 BUT IS ',I8)
 9102 FORMAT(' *ERROR  9102: PROGRAMMING ERROR IN ',A,': UG_COL IS NOT ALLOCATED')
 9103 FORMAT(' *ERROR  9103: PROGRAMMING ERROR IN ',A,': SIZE(UG_COL)=',I8,' IS LESS THAN NDOFG=',I8)

! **********************************************************************************************************************************

      END SUBROUTINE READ_L5A_UG_FOR_SUBCASE


      SUBROUTINE READ_XTIME ( UNT, FILNAM, MESSAG, OUNT )

! Reads STIME from an unformatted file and calls READERR if there is an error condition

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY: READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'READ_XTIME'
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! File description. Used for error messaging
      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM            ! File name

      INTEGER(LONG), INTENT(IN)       :: UNT               ! File unit number
      INTEGER(LONG), INTENT(IN)       :: OUNT(2)           ! File units to write messages to
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: XTIME             ! Time stamp read from an unformatted file




! *********************************************************************************************************************************
      REWIND (UNT)
      READ(UNT,IOSTAT=IOCHK) XTIME
      IF (IOCHK /= 0) THEN
         CALL READERR ( IOCHK, FILNAM, MESSAG, 1, OUNT )
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE READ_XTIME



      SUBROUTINE RESTART_DATA_FOR_L3

! Reads matrices needed when a restart is made in LINK3

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE IOUNT1, ONLY                :  L2G, LINK2G, L2G_MSG, L2GSTAT
      USE IOUNT1, ONLY                :  L2H, LINK2H, L2H_MSG, L2HSTAT
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFL, NTERM_KLL, NTERM_PL
      USE TIMDAT, ONLY                :  TSEC
      USE SPARSE_MATRICES, ONLY       :  I_KLL , J_KLL , KLL ,I_PL , J_PL , PL

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY: READ_MATRIX_1
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      LOGICAL                         :: FILE_EXIST        ! Result from INQUIRE is true if a file exists

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RESTART_DATA_FOR_L3'
      CHARACTER(  1*BYTE)             :: CLOSE_IT   = 'Y'  ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER(  1*BYTE)             :: NTERM_RD   = 'Y'  ! 'Y' or 'N' Input to subr READ_MATRIX_1
      CHARACTER(  1*BYTE)             :: OPND       = 'N'  ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not

      INTEGER(LONG)                   :: IERR              ! Local error count





! **********************************************************************************************************************************
      IERR = 0

! Read KLL

      INQUIRE ( FILE=LINK2G, EXIST=FILE_EXIST )
      IF ((FILE_EXIST) .AND. (NTERM_KLL > 0)) THEN
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KLL ', CR13
         WRITE(SC1,*) CR13
         CALL ALLOCATE_SPARSE_MAT ( 'KLL', NDOFL, NTERM_KLL, SUBR_NAME )
         CALL READ_MATRIX_1 ( LINK2G, L2G, OPND, CLOSE_IT, L2GSTAT, L2G_MSG,'KLL', NTERM_KLL, NTERM_RD, NDOFL, I_KLL, J_KLL, KLL)
      ELSE
         IERR = IERR + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,101) L2G_MSG, LINK2G
         WRITE(F06,101) L2G_MSG, LINK2G
      ENDIF

! Read PL

      INQUIRE ( FILE=LINK2H, EXIST=FILE_EXIST )
      IF ((FILE_EXIST) .AND. (NTERM_PL > 0)) THEN
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PL  ', CR13
         WRITE(SC1,*) CR13
         CALL ALLOCATE_SPARSE_MAT ( 'PL', NDOFL, NTERM_PL, SUBR_NAME )
         CALL READ_MATRIX_1 ( LINK2H, L2H, OPND, CLOSE_IT, L2HSTAT, L2H_MSG,'PL', NTERM_PL, NTERM_RD, NDOFL, I_PL, J_PL, PL)
      ELSE
         IERR = IERR + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,101) L2G_MSG, LINK2G
         WRITE(F06,101) L2G_MSG, LINK2G
      ENDIF

      IF (IERR > 0) THEN
         WRITE(SC1, * ) ' PROCESSING STOPPED. SEE F06 FILE FOR DETAILS'
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *ERROR      : THE FOLLOWING FILE FOR THE ',A                                                                        &
                    ,/,14X,' EITHER DOES NOT EXIST OR IS NULL: '                                                                   &
                    ,/,15X,A)

12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE RESTART_DATA_FOR_L3

   END MODULE RESTART_FILE_IO
