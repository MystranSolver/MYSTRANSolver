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

   MODULE REDUCTION_N_TO_F

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: REDUCE_N_FS

   CONTAINS

      SUBROUTINE REDUCE_N_FS

! Call routines to reduce stiffness, mass, loads from N-set to F, S-sets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1H, LINK1H, L1H_MSG, L2C, LINK2C, L2C_MSG, SC1, WRT_ERR
      USE SCONTR, ONLY                :  LINKNO    , NDOFF, NDOFG, NDOFN, NDOFS, NDOFSE, NSUB,                                     &
                                         NTERM_KNN , NTERM_KFF , NTERM_KFS , NTERM_KSS , NTERM_KSSe ,                              &
                                         NTERM_KNND, NTERM_KFFD, NTERM_KFSD, NTERM_KSSD, NTERM_KSSDe,                              &
                                         NTERM_QSYS, NTERM_PN  , NTERM_PF  , NTERM_PS  , NTERM_MNN  ,                              &
                                         NTERM_MFF , NTERM_MFS , NTERM_MSS , SOL_NAME, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC, YEAR, MONTH, DAY, HOUR, MINUTE, SEC, SFRAC, STIME
      USE CONSTANTS_1, ONLY           :  ONE
      USE DOF_TABLES, ONLY            :  TDOFI
      USE RIGID_BODY_DISP_MATS, ONLY  :  RBGLOBAL_GSET, RBGLOBAL_NSET, RBGLOBAL_FSET
      USE PARAMS, ONLY                :  EQCHK_OUTPUT, MATSPARS, PRTSTIFD, PRTSTIFF, PRTMASS, PRTFOR
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE SPARSE_MATRICES, ONLY       :  I_KNN  , J_KNN  , KNN  , I_KFF  , J_KFF  , KFF  , I_KFS  , J_KFS  , KFS  ,                &
                                         I_KSS  , J_KSS  , KSS  , I_KSSe , J_KSSe , KSSe ,                                         &
                                         I_KNND , J_KNND , KNND , I_KFFD , J_KFFD , KFFD , I_KFSD , J_KFSD , KFSD ,                &
                                         I_KSSD , J_KSSD , KSSD , I_KSSDe, J_KSSDe, KSSDe,                                         &
                                         I_MNN  , J_MNN  , MNN  , I_MFF  , J_MFF  , MFF  , I_MFS  , J_MFS  , MFS  ,                &
                                         I_MSS  , J_MSS  , MSS  ,                                                                  &
                                         I_PN   , J_PN   , PN   , I_PF   , J_PF   , PF   , I_PS   , J_PS   , PS   ,                &
                                         I_MSF  , J_MSF  , MSF  ,                                                                  &
                                         I_QSYS , J_QSYS , QSYS
      USE SPARSE_MATRICES, ONLY       :  SYM_KFF, SYM_KSSe
      USE COL_VECS, ONLY              :  YSe
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE SCRATCH_MATRICES

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_VEC
      USE COL_VEC_LIFECYCLE, ONLY     :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN, READERR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_VECTOR
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE SPARSE_FULL_MULTIPLICATION, ONLY:  MATMULT_SFS, MATMULT_SFS_NTERM
      USE MATRIX_FILE_IO, ONLY        :  WRITE_MATRIX_1, WRITE_SPARSE_CRS
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_MATRIX_DIAG_STATS
      USE RIGID_BODY_STORAGE_LIFECYCLE, ONLY:  ALLOCATE_RBGLOBAL, DEALLOCATE_RBGLOBAL
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE REDUCTION_CHECKS, ONLY      :  STIFF_MAT_EQUIL_CHK

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_N_FS'
                                                               ! If 'Y' then matrix is stored as all nonzeros on & above diag.
      CHARACTER(  1*BYTE)             :: CLOSE_IT              ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close file or not
      CHARACTER(  8*BYTE)             :: CLOSE_STAT            ! Char constant for the CLOSE status of a file
      CHARACTER( 44*BYTE)             :: MODNAM                ! Name to write to screen to describe module being run
      CHARACTER(132*BYTE)             :: MATRIX_NAME           ! Name of matrix for printout

      INTEGER(LONG)                   :: AROW_MAX_TERMS        ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: DO_WHICH_CODE_FRAG    ! 1 or 2 depending on which seg of code to run (depends on BUCKLING)
      INTEGER(LONG)                   :: F_SET_COL             ! Col no. in array TDOFI where the F-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: F_SET_DOF             ! F-set DOF number
      INTEGER(LONG)                   :: IERROR                ! Error count
      INTEGER(LONG)                   :: IOCHK                 ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: I,J                   ! DO loop indices
      INTEGER(LONG)     , PARAMETER   :: NUM_YS_COLS = 1       ! Variable for number of cols in array YSe
      INTEGER(LONG)                   :: OUNT(2)               ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: PART_VEC_F(NDOFF)     ! Partitioning vector (1's for all of F set)
      INTEGER(LONG)                   :: PART_VEC_N_FS(NDOFN)  ! Partitioning vector (N set into F and S sets)
      INTEGER(LONG)                   :: PART_VEC_S(NDOFS)     ! Partitioning vector (1's for all of S set)
      INTEGER(LONG)                   :: PART_VEC_S_SzSe(NDOFS)! Partitioning vector (S set into SZ and SE sets)
      INTEGER(LONG)                   :: PART_VEC_SUB(NSUB)    ! Partitioning vector (1's for all subcases)
      INTEGER(LONG)                   :: REC_NO                ! Record number when reading a file


      REAL(DOUBLE)                    :: KFF_DIAG(NDOFF)       ! Diagonal terms from KFF
      REAL(DOUBLE)                    :: KFF_MAX_DIAG          ! Max diag term from KFF


      ! ensure output units are set
      OUNT(1) = ERR
      OUNT(2) = F06



! **********************************************************************************************************************************
! Set partitioning vectors

      CALL PARTITION_VEC ( NDOFN, 'N ', 'F ', 'S ', PART_VEC_N_FS )

      CALL PARTITION_VEC ( NDOFS, 'S ', 'SZ', 'SE', PART_VEC_S_SzSe)

      DO I=1,NDOFF
         PART_VEC_F(I) = 1
      ENDDO

      DO I=1,NDOFS
         PART_VEC_S(I) = 1
      ENDDO

      DO I=1,NSUB
         PART_VEC_SUB = 1
      ENDDO

! Read enforced displ's if there are any

      CALL ALLOCATE_COL_VEC ( 'YSe', NDOFSE, SUBR_NAME )
      IF (NDOFSE > 0) THEN

         CALL FILE_OPEN ( L1H, LINK1H, OUNT, 'OLD', L1H_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

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
            CALL OUTA_HERE ( 'Y' )                                 ! Quit due to read errors in YSe array file
         ENDIF

         CALL FILE_CLOSE ( L1H, LINK1H, 'KEEP' )

         IF (DEBUG(26) == 1) THEN
            CALL WRITE_VECTOR ( 'SE-SET YS ENFORCED DISPLS','DISPL', NDOFSE, YSe )
         ENDIF

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

         IF (NDOFS > 0) THEN                                  ! If NDOFS > 0 reduce KNN to KFF

! Reduce KNN to KFF

            IF (NTERM_KNN > 0) THEN

               CALL OURTIM
               MODNAM = '  REDUCE KNN TO KFF (PARTITION, ONLY)'
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_KNN_TO_KFF ( PART_VEC_N_FS, PART_VEC_S_SzSe, PART_VEC_F, PART_VEC_S )
            ELSE

               NTERM_KFF = 0
               NTERM_KFS = 0
               NTERM_KSS = 0
               CALL ALLOCATE_SPARSE_MAT ( 'KFF', NDOFF, NTERM_KFF, SUBR_NAME )
               !WRITE(SC1,*) CR13

            ENDIF

! Reduce MNN to MFF

            IF (NTERM_MNN > 0) THEN

               CALL OURTIM
               MODNAM = '  REDUCE MNN TO MFF (PARTITION, ONLY)'
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_MNN_TO_MFF ( PART_VEC_N_FS )

            ELSE

               NTERM_MFF = 0
               NTERM_MFS = 0
               NTERM_MSS = 0
               CALL ALLOCATE_SPARSE_MAT ( 'MFF', NDOFF, NTERM_MFF, SUBR_NAME )

            ENDIF

! Reduce PN to PF.

            IF ((SOL_NAME(1:5) /= 'MODES') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')) THEN

               CALL OURTIM
               IF (MATSPARS == 'Y') THEN
                  MODNAM = '  REDUCE PN TO PF (SPARSE MATRIX ROUTINES)'
               ELSE
                  MODNAM = '  REDUCE PN TO PF (FULL MATRIX ROUTINES)'
               ENDIF
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_PN_TO_PF ( PART_VEC_N_FS, PART_VEC_SUB )

            ENDIF

         ELSE

! There is no S-set, so equate N, F sets

            CALL OURTIM
            MODNAM = '  EQUATING F-SET TO N-SET'
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            NDOFF     = NDOFN

            NTERM_KFF = NTERM_KNN
            NTERM_KFS = 0
            NTERM_KSS = 0

            NTERM_MFF = NTERM_MNN
            NTERM_MFS = 0
            NTERM_MSS = 0

            NTERM_PF  = NTERM_PN
            NTERM_PS  = 0

            CALL ALLOCATE_SPARSE_MAT ( 'KFF', NDOFF, NTERM_KFF, SUBR_NAME )

!xx         DO I=1,NDOFF+1
!xx            I_KFF(I) = I_KNN(I)
!xx         ENDDO
!xx
!xx         DO I=1,NTERM_KFF
!xx            J_KFF(I) = J_KNN(I)
!xx              KFF(I) =   KNN(I)
!xx         ENDDO
!xx
            CALL ALLOCATE_SPARSE_MAT ( 'MFF', NDOFF, NTERM_MFF, SUBR_NAME )

!xx         DO I=1,NDOFF+1
!xx            I_MFF(I) = I_MNN(I)
!xx         ENDDO
!xx
!xx         DO I=1,NTERM_MFF
!xx            J_MFF(I) = J_MNN(I)
!xx              MFF(I) =   MNN(I)
!xx         ENDDO

            IF ((SOL_NAME(1:5) /= 'MODES') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')) THEN

               CALL ALLOCATE_SPARSE_MAT ( 'PF', NDOFF, NTERM_PF, SUBR_NAME )

!xx            DO I=1,NDOFF+1
!xx               I_PF(I)  = I_PN(I)
!xx            ENDDO
!xx
!xx            DO I=1,NTERM_PF
!xx               J_PF(I) = J_PN(I)
!xx                 PF(I) =   PN(I)
!xx            ENDDO

            ENDIF

         ENDIF

! Deallocate N-set arrays

         MODNAM = '   DEALLOCATE SOME ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KNN', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KNN' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MNN', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MNN' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PN ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PN' )
         WRITE(SC1,*) CR13

! Calc QSYS = KSSe*YSe and write to L2C for constraint force recovery in LINK9.

            IF (NTERM_KSSe > 0) THEN                          ! Calc QSYS = KSSe * YSe

               CALL MATMULT_SFS_NTERM ( 'KSSe', NDOFS, NTERM_KSSe, SYM_KSSe, I_KSSe, J_KSSe                                        &
                                        ,'YSe', NDOFSE, NUM_YS_COLS, YSe, AROW_MAX_TERMS, 'QSYS', NTERM_QSYS )

               CALL ALLOCATE_SPARSE_MAT ( 'QSYS', NDOFS, NTERM_QSYS, SUBR_NAME )

               IF ((DEBUG(24) == 2) .OR.(DEBUG(24) == 3)) THEN
                  MATRIX_NAME = 'STIFFNESS PARTITION KSSe (columns of KSS for enforced displs)'
                  CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'S ', 'SE', NTERM_KSSe, NDOFS, I_KSSe, J_KSSe, KSSe )
               ENDIF

               IF (NTERM_QSYS > 0) THEN

                  CALL MATMULT_SFS ( 'KSSe', NDOFS, NTERM_KSSe, SYM_KSSe, I_KSSe, J_KSSe, KSSe                                     &
                                    ,'YSe', NDOFSE, NUM_YS_COLS, YSe, AROW_MAX_TERMS, 'QSYS', ONE, NTERM_QSYS, I_QSYS, J_QSYS, QSYS)

                  CLOSE_IT   = 'Y'
                  CLOSE_STAT = 'KEEP'
                  CALL WRITE_MATRIX_1 ( LINK2C, L2C, CLOSE_IT, CLOSE_STAT, L2C_MSG, 'QSYS', NTERM_QSYS, NDOFS, I_QSYS, J_QSYS, QSYS)

               ENDIF

            ENDIF

! Print out QSYS

         IF ((DEBUG(25) == 2) .OR.(DEBUG(25) == 3)) THEN      ! DEBUG(25) controls output of QSYS
            IF (NTERM_QSYS > 0) THEN
               MATRIX_NAME = 'MATRIX QSYS = KSSe*YSe (portion of SPC forces due to enforced displs)'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'S ', '  ', NTERM_QSYS, NDOFS, I_QSYS, J_QSYS, QSYS )
            ENDIF
         ENDIF

         MODNAM = '   DEALLOCATE SOME ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KSSe', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KSSe' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate QSYS', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'QSYS' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate YSe ', CR13   ;   CALL DEALLOCATE_COL_VEC ( 'YSe' )
         WRITE(SC1,*) CR13

! Print out stiffness matrix partitions, if requested

         IF (( PRTSTIFF(3) == 1) .OR. ( PRTSTIFF(3) == 3)) THEN
            IF (NTERM_KFF > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KFF'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'F ', 'F ', NTERM_KFF, NDOFF, I_KFF, J_KFF, KFF )
            ENDIF
         ENDIF

         IF (( PRTSTIFF(3) == 2) .OR. ( PRTSTIFF(3) == 3)) THEN
            IF (NTERM_KFS > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KFS'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'F ', 'S ', NTERM_KFS, NDOFF, I_KFS, J_KFS, KFS )
            ENDIF
            IF (NTERM_KSS > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KSS'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'S ', 'S ', NTERM_KSS, NDOFS, I_KSS, J_KSS, KSS )
            ENDIF
         ENDIF

         WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KFS', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KFS' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KSF', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KSF' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KSS', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KSS' )
         WRITE(SC1,*) CR13

! Write matrix diagonal and stats, if requested.
! NOTE: call this subr even if PRTSTIFFD(3) = 0 since we need KFF_DIAG, KFF_MAX_DIAG for the equilibrium check

         IF (NDOFF > 0) THEN
            CALL GET_MATRIX_DIAG_STATS ( 'KFF', 'F ', NDOFF, NTERM_KFF, I_KFF, J_KFF, KFF, PRTSTIFD(3), KFF_DIAG, KFF_MAX_DIAG )
         ENDIF

! Print out mass matrix partitions, if requested

         IF (( PRTMASS(3) == 1) .OR. ( PRTMASS(3) == 3)) THEN
            IF (NTERM_MFF > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MFF'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'F ', 'F ', NTERM_MFF, NDOFF, I_MFF, J_MFF, MFF )
            ENDIF
         ENDIF

         IF (( PRTMASS(3) == 2) .OR. ( PRTMASS(3) == 3)) THEN
            IF (NTERM_MFS > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MFS'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'F ', 'S ', NTERM_MFS, NDOFF, I_MFS, J_MFS, MFS )
            ENDIF
            IF (NTERM_MSS > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MSS'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'S ', 'S ', NTERM_MSS, NDOFS, I_MSS, J_MSS, MSS )
            ENDIF
         ENDIF

         WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MFS', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MFS' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MSF', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MSF' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MSS', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MSS' )
         WRITE(SC1,*) CR13

! Print out load matrix partitions, if requested

         IF (( PRTFOR(3) == 1) .OR. ( PRTFOR(3) == 3)) THEN
            IF (NTERM_PF  > 0) THEN
               MATRIX_NAME = 'LOAD MATRIX PF'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'F ', 'SUBCASE', NTERM_PF, NDOFF, I_PF, J_PF, PF )
            ENDIF
         ENDIF

         IF (( PRTFOR(3) == 2) .OR. ( PRTFOR(3) == 3)) THEN
            IF (NTERM_PS  > 0) THEN
               MATRIX_NAME = 'LOAD MATRIX PS'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'S ', 'SUBCASE', NTERM_PS, NDOFS, I_PS, J_PS, PS )
            ENDIF
         ENDIF

         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PS ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PS' )
         WRITE(SC1,*) CR13

! Do equilibrium check on the F-set stiffness matrix, if requested

         IF ((EQCHK_OUTPUT(3) > 0) .OR. (EQCHK_OUTPUT(4) > 0) .OR. (EQCHK_OUTPUT(5) > 0)) THEN
            CALL ALLOCATE_RBGLOBAL ( 'F ', SUBR_NAME )
            IF (NDOFS > 0) THEN
               CALL TDOF_COL_NUM ( 'F ', F_SET_COL )
               DO I=1,NDOFG
                  F_SET_DOF = TDOFI(I,F_SET_COL)
                  IF (F_SET_DOF > 0) THEN
                     DO J=1,6
                        RBGLOBAL_FSET(F_SET_DOF,J) = RBGLOBAL_GSET(I,J)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               DO I=1,NDOFF
                  DO J=1,6
                     RBGLOBAL_FSET(I,J) = RBGLOBAL_NSET(I,J)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF

         IF (EQCHK_OUTPUT(3) > 0) THEN
            CALL OURTIM
            MODNAM = '   EQUILIBRIUM CHECK ON KFF                '
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
            CALL STIFF_MAT_EQUIL_CHK ( EQCHK_OUTPUT(3),'F ', SYM_KFF, NDOFF, NTERM_KFF, I_KFF, J_KFF, KFF, KFF_DIAG, KFF_MAX_DIAG, &
                                       RBGLOBAL_FSET)
         ENDIF
         CALL DEALLOCATE_RBGLOBAL ( 'N ' )

! **********************************************************************************************************************************
      ELSE                                                 ! This is BUCKLING with LOAD_ISTEP = 2 (eigen part of BUCKLING)

         IF (NDOFS > 0) THEN                                  ! If NDOFS > 0 reduce KNN to KFF

! Reduce KNND to KFFD

            IF (NTERM_KNND > 0) THEN

               CALL OURTIM
               MODNAM = '  REDUCE KNND TO KFFD (PARTITION, ONLY)'
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_KNND_TO_KFFD ( PART_VEC_N_FS, PART_VEC_S_SzSe, PART_VEC_F, PART_VEC_S )

            ELSE

               NTERM_KFFD = 0
               NTERM_KFSD = 0
               NTERM_KSSD = 0
               CALL ALLOCATE_SPARSE_MAT ( 'KFFD', NDOFF, NTERM_KFFD, SUBR_NAME )

            ENDIF

         ELSE

! There is no S-set, so equate N, F sets

            CALL OURTIM
            MODNAM = '  EQUATING F-SET TO N-SET'
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            NDOFF     = NDOFN

            NTERM_KFFD = NTERM_KNND
            NTERM_KFSD = 0
            NTERM_KSSD = 0

            CALL ALLOCATE_SPARSE_MAT ( 'KFFD', NDOFF, NTERM_KFFD, SUBR_NAME )

         ENDIF

! Deallocate N-set arrays

         MODNAM = '   DEALLOCATE SOME ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KNND', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KNND' )
         WRITE(SC1,*) CR13

! Print out stiffness matrix partitions, if requested

         IF (( PRTSTIFF(3) == 1) .OR. ( PRTSTIFF(3) == 3)) THEN
            IF (NTERM_KFFD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KFFD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'F ', 'F ', NTERM_KFFD, NDOFF, I_KFFD, J_KFFD, KFFD )
            ENDIF
         ENDIF

         IF (( PRTSTIFF(3) == 2) .OR. ( PRTSTIFF(3) == 3)) THEN
            IF (NTERM_KFSD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KFSD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'F ', 'S ', NTERM_KFSD, NDOFF, I_KFSD, J_KFSD, KFSD )
            ENDIF
            IF (NTERM_KSSD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KSSD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'S ', 'S ', NTERM_KSSD, NDOFS, I_KSSD, J_KSSD, KSSD )
            ENDIF
         ENDIF

         WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate YSe ', CR13   ;   CALL DEALLOCATE_COL_VEC    ( 'YSe' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KFSD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KFSD' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KSFD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KSFD' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KSSD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KSSD' )
         WRITE(SC1,*) CR13

      ENDIF



      RETURN

! **********************************************************************************************************************************
 2092 FORMAT(4X,A44,20X,I2,':',I2,':',I2,'.',I3)

 9995 FORMAT(/,' PROCESSING ENDED IN LINK ',I3,' DUE TO ABOVE ',I8,' ERRORS')

12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_N_FS


      SUBROUTINE REDUCE_KNND_TO_KFFD ( PART_VEC_N_FS, PART_VEC_S_SzSe, PART_VEC_F, PART_VEC_S )

! Call routines to reduce the KNND differential stiffness matrix from the N-set to the F, S-sets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L2B, LINK2B, L2B_MSG
      USE SCONTR, ONLY                :  FATAL_ERR, NDOFN, NDOFF, NDOFS, NDOFSE, NTERM_KNND, NTERM_KFFD, NTERM_KFSD, NTERM_KSSD,   &
                                         NTERM_KFSDe, NTERM_KSSDe, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE SPARSE_MATRICES, ONLY       :  I_KNND, J_KNND, KNND, I_KFFD, J_KFFD, KFFD, I_KFSD, J_KFSD, KFSD, I_KFSDe, J_KFSDe, KFSDe,&
                                         I_KSFD, J_KSFD, KSFD, I_KSSD, J_KSSD, KSSD, I_KSSDe, J_KSSDe, KSSDe
      USE SPARSE_MATRICES, ONLY       :  SYM_KNND, SYM_KFFD, SYM_KFSD, SYM_KFSDe, SYM_KSSD, SYM_KSSD, SYM_KSSDe
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATTRNSP_SS
      USE MATRIX_FILE_IO, ONLY        :  WRITE_MATRIX_1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_KNND_TO_KFFD'
      CHARACTER(  1*BYTE)             :: CLOSE_IT               ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close file or not
      CHARACTER(  8*BYTE)             :: CLOSE_STAT             ! Char constant for the CLOSE status of a file

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_F(NDOFF)      ! Partitioning vector (1's for all of F set)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_N_FS(NDOFN)   ! Partitioning vector (N set into F and S sets)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_S(NDOFS)      ! Partitioning vector (1's for all of S set)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_S_SzSe(NDOFS) ! Partitioning vector (S set into SZ and SE sets)
      INTEGER(LONG)                   :: KFFD_ROW_MAX_TERMS     ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KFSD_ROW_MAX_TERMS     ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
!xx   INTEGER(LONG)                   :: KSFD_ROW_MAX_TERMS     ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KFSDe_ROW_MAX_TERMS    ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KSSD_ROW_MAX_TERMS     ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KSSDe_ROW_MAX_TERMS    ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: NTERM_KSFD             ! Number of nonzeros in sparse matrix KSFD (should = NTERM_KFSD)
      INTEGER(LONG), PARAMETER        :: NUM1        = 1        ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2        ! Used in subr's that partition matrices


      INTRINSIC                       :: DABS



! **********************************************************************************************************************************

! Partition KFFD from KNND. This is final KFFD.

      IF (NDOFF > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KNND', NTERM_KNND, NDOFN, NDOFN, SYM_KNND, I_KNND, J_KNND,      PART_VEC_N_FS, PART_VEC_N_FS,  &
                                    NUM1, NUM1, KFFD_ROW_MAX_TERMS, 'KFFD', NTERM_KFFD, SYM_KFFD )

         CALL ALLOCATE_SPARSE_MAT ( 'KFFD', NDOFF, NTERM_KFFD, SUBR_NAME )

         IF (NTERM_KFFD > 0) THEN
            CALL PARTITION_SS ( 'KNND', NTERM_KNND, NDOFN, NDOFN, SYM_KNND, I_KNND, J_KNND, KNND, PART_VEC_N_FS, PART_VEC_N_FS,    &
                                 NUM1, NUM1, KFFD_ROW_MAX_TERMS, 'KFFD', NTERM_KFFD, NDOFF, SYM_KFFD, I_KFFD, J_KFFD, KFFD )
         ENDIF

      ENDIF

! Partition KFSD from KNND. Then partition KFSDe from KFSD

      IF ((NDOFF > 0) .AND. (NDOFS > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'KNND', NTERM_KNND, NDOFN, NDOFN, SYM_KNND, I_KNND, J_KNND,      PART_VEC_N_FS, PART_VEC_N_FS,  &
                                    NUM1, NUM2, KFSD_ROW_MAX_TERMS, 'KFSD', NTERM_KFSD, SYM_KFSD )

         CALL ALLOCATE_SPARSE_MAT ( 'KFSD', NDOFF, NTERM_KFSD, SUBR_NAME )

         IF (NTERM_KFSD > 0) THEN
            CALL PARTITION_SS ( 'KNND', NTERM_KNND, NDOFN, NDOFN, SYM_KNND, I_KNND, J_KNND, KNND, PART_VEC_N_FS, PART_VEC_N_FS,    &
                                 NUM1, NUM2, KFSD_ROW_MAX_TERMS, 'KFSD', NTERM_KFSD, NDOFF, SYM_KFSD, I_KFSD, J_KFSD, KFSD )

            IF (NDOFSE > 0) THEN

               CALL PARTITION_SS_NTERM ( 'KFSD', NTERM_KFSD, NDOFF, NDOFS, SYM_KFSD, I_KFSD, J_KFSD, PART_VEC_F, PART_VEC_S_SzSe   &
                                         ,NUM1, NUM2, KFSDe_ROW_MAX_TERMS, 'KFSDe', NTERM_KFSDe, SYM_KFSDe )

               CALL ALLOCATE_SPARSE_MAT ( 'KFSDe', NDOFF, NTERM_KFSDe, SUBR_NAME )

               IF (NTERM_KFSDe > 0) THEN
                  CALL PARTITION_SS ('KFSD', NTERM_KFSD, NDOFF, NDOFS, SYM_KFSD, I_KFSD, J_KFSD, KFSD, PART_VEC_F, PART_VEC_S_SzSe,&
                                  NUM1, NUM2, KFSDe_ROW_MAX_TERMS, 'KFSDe', NTERM_KFSDe, NDOFF, SYM_KFSDe, I_KFSDe, J_KFSDe, KFSDe )
               ENDIF

            ENDIF

         ENDIF

      ENDIF

! Partition KSFD from KNND and write KSFD to L2B for constraint force recovery in LINK9.

      IF ((NDOFF > 0) .AND. (NDOFS > 0)) THEN

!xx      CALL PARTITION_SS_NTERM ( 'KNND', NTERM_KNND, NDOFN, NDOFN, SYM_KNND, I_KNND, J_KNND,      PART_VEC_N_FS, PART_VEC_N_FS,  &
!xx                                 NUM2, NUM1, KSFD_ROW_MAX_TERMS, 'KSFD', NTERM_KSFD, SYM_KFSD )

!xx      IF (NTERM_KSFD /= NTERM_KFSD) THEN
!xx         FATAL_ERR = FATAL_ERR + 1
!xx         WRITE(ERR,936) SUBR_NAME, NTERM_KSFD, NTERM_KFSD
!xx         WRITE(F06,936) SUBR_NAME, NTERM_KSFD, NTERM_KFSD
!xx         CALL OUTA_HERE ( 'Y' )
!xx      ENDIF

!xx      CALL ALLOCATE_SPARSE_MAT ( 'KSFD', NDOFS, NTERM_KSFD, SUBR_NAME )

!xx      IF (NTERM_KFSD > 0) THEN
!xx         CALL PARTITION_SS ( 'KNND', NTERM_KNND, NDOFN, NDOFN, SYM_KNND, I_KNND, J_KNND, KNND, PART_VEC_N_FS, PART_VEC_N_FS,    &
!xx                              NUM2, NUM1, KSFD_ROW_MAX_TERMS, 'KSFD', NTERM_KSFD, NDOFS, SYM_KFSD, I_KSFD, J_KSFD, KSFD )

!xx         CLOSE_IT   = 'Y'
!xx         CLOSE_STAT = 'KEEP'
!xx         CALL WRITE_MATRIX_1 ( LINK2B, L2B, CLOSE_IT, CLOSE_STAT, L2B_MSG, 'KSFD', NTERM_KSFD, NDOFS, I_KSFD, J_KSFD, KSFD )

!xx      ENDIF

         NTERM_KSFD = NTERM_KFSD
         CALL ALLOCATE_SPARSE_MAT ( 'KSFD', NDOFS, NTERM_KSFD, SUBR_NAME )

         IF (NTERM_KSFD > 0) THEN
            CALL MATTRNSP_SS ( NDOFF, NDOFS, NTERM_KFSD, 'KFSD', I_KFSD, J_KFSD, KFSD, 'KSFD', I_KSFD, J_KSFD, KSFD )

            CLOSE_IT   = 'Y'
            CLOSE_STAT = 'KEEP'
            CALL WRITE_MATRIX_1 ( LINK2B, L2B, CLOSE_IT, CLOSE_STAT, L2B_MSG, 'KSFD', NTERM_KSFD, NDOFS, I_KSFD, J_KSFD, KSFD )

         ENDIF

      ENDIF

! Partition KSSD from KNND. Then partition KSSDe from KSSD

      IF (NDOFS > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KNND', NTERM_KNND, NDOFN, NDOFN, SYM_KNND, I_KNND, J_KNND,      PART_VEC_N_FS, PART_VEC_N_FS,  &
                                    NUM2, NUM2, KSSD_ROW_MAX_TERMS, 'KSSD', NTERM_KSSD, SYM_KSSD )

         CALL ALLOCATE_SPARSE_MAT ( 'KSSD', NDOFS, NTERM_KSSD, SUBR_NAME )

         IF (NTERM_KSSD > 0) THEN

            CALL PARTITION_SS ( 'KNND', NTERM_KNND, NDOFN, NDOFN, SYM_KNND, I_KNND, J_KNND, KNND, PART_VEC_N_FS, PART_VEC_N_FS,    &
                                 NUM2, NUM2, KSSD_ROW_MAX_TERMS, 'KSSD', NTERM_KSSD, NDOFS, SYM_KSSD, I_KSSD, J_KSSD, KSSD )

            IF (NDOFSE > 0) THEN

               CALL PARTITION_SS_NTERM ( 'KSSD', NTERM_KSSD, NDOFS, NDOFS, SYM_KSSD, I_KSSD, J_KSSD, PART_VEC_S, PART_VEC_S_SzSe   &
                                         ,NUM1, NUM2, KSSDe_ROW_MAX_TERMS, 'KSSDe', NTERM_KSSDe, SYM_KSSDe )

               CALL ALLOCATE_SPARSE_MAT ( 'KSSDe', NDOFS, NTERM_KSSDe, SUBR_NAME )

               IF (NTERM_KSSDe > 0) THEN
                  CALL PARTITION_SS ('KSSD', NTERM_KSSD, NDOFS, NDOFS, SYM_KSSD, I_KSSD, J_KSSD, KSSD, PART_VEC_S, PART_VEC_S_SzSe,&
                                  NUM1, NUM2, KSSDe_ROW_MAX_TERMS, 'KSSDe', NTERM_KSSDe, NDOFS, SYM_KSSDe, I_KSSDe, J_KSSDe, KSSDe )
               ENDIF

            ENDIF

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
  936 FORMAT(' *ERROR   936: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF NONZERO TERMS IN SPARSE MATRIX KSFD = ',I12,' IS NOT EQUAL TO THOSE IN MATRIX KFSD = ', &
                       I12)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_KNND_TO_KFFD


      SUBROUTINE REDUCE_KNN_TO_KFF ( PART_VEC_N_FS, PART_VEC_S_SzSe, PART_VEC_F, PART_VEC_S )

! Call routines to reduce the KNN linear stiffness matrix from the N-set to the F, S-sets. See Appendix B to the MYSTRAN User's
! Reference Manual for the derivation of the reduction equations.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L2B, LINK2B, L2B_MSG, SC1
      USE SCONTR, ONLY                :  FATAL_ERR, NDOFN, NDOFF, NDOFS, NDOFSE, NTERM_KNN, NTERM_KFF, NTERM_KFS, NTERM_KSS,       &
                                         NTERM_KFSe, NTERM_KSSe, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE SPARSE_MATRICES, ONLY       :  I_KNN, J_KNN, KNN, I_KFF, J_KFF, KFF, I_KFS, J_KFS, KFS, I_KFSe, J_KFSe, KFSe,            &
                                         I_KSF, J_KSF, KSF, I_KSS, J_KSS, KSS, I_KSSe, J_KSSe, KSSe
      USE SPARSE_MATRICES, ONLY       :  SYM_KNN, SYM_KFF, SYM_KFS, SYM_KFSe, SYM_KSS, SYM_KSS, SYM_KSSe
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATTRNSP_SS
      USE MATRIX_FILE_IO, ONLY        :  WRITE_MATRIX_1

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_KNN_TO_KFF'
      CHARACTER(  1*BYTE)             :: CLOSE_IT               ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close file or not
      CHARACTER(  8*BYTE)             :: CLOSE_STAT             ! Char constant for the CLOSE status of a file

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_F(NDOFF)      ! Partitioning vector (1's for all of F set)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_N_FS(NDOFN)   ! Partitioning vector (N set into F and S sets)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_S(NDOFS)      ! Partitioning vector (1's for all of S set)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_S_SzSe(NDOFS) ! Partitioning vector (S set into SZ and SE sets)
      INTEGER(LONG)                   :: KFF_ROW_MAX_TERMS      ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KFS_ROW_MAX_TERMS      ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KSF_ROW_MAX_TERMS      ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KFSe_ROW_MAX_TERMS     ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KSS_ROW_MAX_TERMS      ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KSSe_ROW_MAX_TERMS     ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: NTERM_KSF              ! Number of nonzeros in sparse matrix KSF (should = NTERM_KFS)
      INTEGER(LONG), PARAMETER        :: NUM1        = 1        ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2        ! Used in subr's that partition matrices


      INTRINSIC                       :: DABS



! **********************************************************************************************************************************

! Partition KFF from KNN. This is final KFF.

      IF (NDOFF > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KNN', NTERM_KNN, NDOFN, NDOFN, SYM_KNN, I_KNN, J_KNN,      PART_VEC_N_FS, PART_VEC_N_FS,       &
                                    NUM1, NUM1, KFF_ROW_MAX_TERMS, 'KFF', NTERM_KFF, SYM_KFF )

         CALL ALLOCATE_SPARSE_MAT ( 'KFF', NDOFF, NTERM_KFF, SUBR_NAME )

         IF (NTERM_KFF > 0) THEN
            CALL PARTITION_SS ( 'KNN', NTERM_KNN, NDOFN, NDOFN, SYM_KNN, I_KNN, J_KNN, KNN, PART_VEC_N_FS, PART_VEC_N_FS,          &
                                 NUM1, NUM1, KFF_ROW_MAX_TERMS, 'KFF', NTERM_KFF, NDOFF, SYM_KFF, I_KFF, J_KFF, KFF )
         ENDIF

      ENDIF

! Partition KFS from KNN. Then partition KFSe from KFS

      IF ((NDOFF > 0) .AND. (NDOFS > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'KNN', NTERM_KNN, NDOFN, NDOFN, SYM_KNN, I_KNN, J_KNN,      PART_VEC_N_FS, PART_VEC_N_FS,       &
                                    NUM1, NUM2, KFS_ROW_MAX_TERMS, 'KFS', NTERM_KFS, SYM_KFS )
         CALL ALLOCATE_SPARSE_MAT ( 'KFS', NDOFF, NTERM_KFS, SUBR_NAME )

         IF (NTERM_KFS > 0) THEN
            CALL PARTITION_SS ( 'KNN', NTERM_KNN, NDOFN, NDOFN, SYM_KNN, I_KNN, J_KNN, KNN, PART_VEC_N_FS, PART_VEC_N_FS,          &
                                 NUM1, NUM2, KFS_ROW_MAX_TERMS, 'KFS', NTERM_KFS, NDOFF, SYM_KFS, I_KFS, J_KFS, KFS )

            IF (NDOFSE > 0) THEN

               CALL PARTITION_SS_NTERM ( 'KFS', NTERM_KFS, NDOFF, NDOFS, SYM_KFS, I_KFS, J_KFS,      PART_VEC_F,    PART_VEC_S_SzSe&
                                         ,NUM1, NUM2, KFSe_ROW_MAX_TERMS, 'KFSe', NTERM_KFSe, SYM_KFSe )

               CALL ALLOCATE_SPARSE_MAT ( 'KFSe', NDOFF, NTERM_KFSe, SUBR_NAME )

               IF (NTERM_KFSe > 0) THEN
                  CALL PARTITION_SS ( 'KFS', NTERM_KFS, NDOFF, NDOFS, SYM_KFS, I_KFS, J_KFS, KFS, PART_VEC_F   , PART_VEC_S_SzSe,  &
                                       NUM1, NUM2, KFSe_ROW_MAX_TERMS, 'KFSe', NTERM_KFSe, NDOFF, SYM_KFSe, I_KFSe, J_KFSe, KFSe )
               ENDIF

            ENDIF

         ENDIF

      ENDIF

! Partition KSF from KNN and write KSF to L2B for constraint force recovery in LINK9.

      IF ((NDOFF > 0) .AND. (NDOFS > 0)) THEN

!xx      CALL PARTITION_SS_NTERM ( 'KNN', NTERM_KNN, NDOFN, NDOFN, SYM_KNN, I_KNN, J_KNN,      PART_VEC_N_FS, PART_VEC_N_FS,       &
!xx                                 NUM2, NUM1, KSF_ROW_MAX_TERMS, 'KSF', NTERM_KSF, SYM_KFS )

!xx      IF (NTERM_KSF /= NTERM_KFS) THEN
!xx         FATAL_ERR = FATAL_ERR + 1
!xx         WRITE(ERR,936) SUBR_NAME, NTERM_KSF, NTERM_KFS
!xx         WRITE(F06,936) SUBR_NAME, NTERM_KSF, NTERM_KFS
!xx         CALL OUTA_HERE ( 'Y' )
!xx      ENDIF

!xx      CALL ALLOCATE_SPARSE_MAT ( 'KSF', NDOFS, NTERM_KSF, SUBR_NAME )

!xx      IF (NTERM_KFS > 0) THEN
!xx         CALL PARTITION_SS ( 'KNN', NTERM_KNN, NDOFN, NDOFN, SYM_KNN, I_KNN, J_KNN, KNN, PART_VEC_N_FS, PART_VEC_N_FS,          &
!xx                              NUM2, NUM1, KSF_ROW_MAX_TERMS, 'KSF', NTERM_KSF, NDOFS, SYM_KFS, I_KSF, J_KSF, KSF )

!xx         CLOSE_IT   = 'Y'
!xx         CLOSE_STAT = 'KEEP'
!xx         CALL WRITE_MATRIX_1 ( LINK2B, L2B, CLOSE_IT, CLOSE_STAT, L2B_MSG, 'KSF', NTERM_KSF, NDOFS, I_KSF, J_KSF, KSF )

!xx      ENDIF

         NTERM_KSF = NTERM_KFS
         CALL ALLOCATE_SPARSE_MAT ( 'KSF', NDOFS, NTERM_KSF, SUBR_NAME )

         IF (NTERM_KSF > 0) THEN
            CALL MATTRNSP_SS ( NDOFF, NDOFS, NTERM_KFS, 'KFS', I_KFS, J_KFS, KFS, 'KSF', I_KSF, J_KSF, KSF )

            CLOSE_IT   = 'Y'
            CLOSE_STAT = 'KEEP'
            CALL WRITE_MATRIX_1 ( LINK2B, L2B, CLOSE_IT, CLOSE_STAT, L2B_MSG, 'KSF', NTERM_KSF, NDOFS, I_KSF, J_KSF, KSF )

         ENDIF

      ENDIF

! Partition KSS from KNN. Then partition KSSe from KSS

      IF (NDOFS > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KNN', NTERM_KNN, NDOFN, NDOFN, SYM_KNN, I_KNN, J_KNN,      PART_VEC_N_FS, PART_VEC_N_FS,       &
                                    NUM2, NUM2, KSS_ROW_MAX_TERMS, 'KSS', NTERM_KSS, SYM_KSS )

         CALL ALLOCATE_SPARSE_MAT ( 'KSS', NDOFS, NTERM_KSS, SUBR_NAME )

         IF (NTERM_KSS > 0) THEN

            CALL PARTITION_SS ( 'KNN', NTERM_KNN, NDOFN, NDOFN, SYM_KNN, I_KNN, J_KNN, KNN, PART_VEC_N_FS, PART_VEC_N_FS,          &
                                 NUM2, NUM2, KSS_ROW_MAX_TERMS, 'KSS', NTERM_KSS, NDOFS, SYM_KSS, I_KSS, J_KSS, KSS )

            IF (NDOFSE > 0) THEN

               CALL PARTITION_SS_NTERM ( 'KSS', NTERM_KSS, NDOFS, NDOFS, SYM_KSS, I_KSS, J_KSS,      PART_VEC_S,    PART_VEC_S_SzSe&
                                         ,NUM1, NUM2, KSSe_ROW_MAX_TERMS, 'KSSe', NTERM_KSSe, SYM_KSSe )

               CALL ALLOCATE_SPARSE_MAT ( 'KSSe', NDOFS, NTERM_KSSe, SUBR_NAME )

               IF (NTERM_KSSe > 0) THEN
                  CALL PARTITION_SS ( 'KSS', NTERM_KSS, NDOFS, NDOFS, SYM_KSS, I_KSS, J_KSS, KSS, PART_VEC_S   , PART_VEC_S_SzSe,  &
                                       NUM1, NUM2, KSSe_ROW_MAX_TERMS, 'KSSe', NTERM_KSSe, NDOFS, SYM_KSSe, I_KSSe, J_KSSe, KSSe )
               ENDIF

            ENDIF

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
  936 FORMAT(' *ERROR   936: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF NONZERO TERMS IN SPARSE MATRIX KSF = ',I12,' IS NOT EQUAL TO THOSE IN MATRIX KFS = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_KNN_TO_KFF


      SUBROUTINE REDUCE_MNN_TO_MFF ( PART_VEC_N_FS )

! Call routines to reduce the MNN mass matrix from the N-set to the F, S-sets. See Appendix B to the MYSTRAN User's Reference Manual
! for the derivation of the reduction equations.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L2S, LINK2S, L2S_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFN, NDOFF, NDOFS, NTERM_MNN, NTERM_MFF, NTERM_MFS, NTERM_MSS
      USE TIMDAT, ONLY                :  TSEC
      USE SPARSE_MATRICES, ONLY       :  I_MNN, J_MNN, MNN, I_MFF, J_MFF, MFF, I_MFS, J_MFS, MFS, I_MSF, J_MSF, MSF,               &
                                         I_MSS, J_MSS, MSS
      USE SPARSE_MATRICES, ONLY       :  SYM_MNN, SYM_MFF, SYM_MFS, SYM_MSS
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATTRNSP_SS
      USE MATRIX_FILE_IO, ONLY        :  WRITE_MATRIX_1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_MNN_TO_MFF'
      CHARACTER(  1*BYTE)             :: CLOSE_IT               ! Input to subr READ_MATRIX_1. 'Y'/'N' whether to close a file
      CHARACTER(  8*BYTE)             :: CLOSE_STAT             ! Char constant for the CLOSE status of a file

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_N_FS(NDOFN)   ! Partitioning vector (N set into F and S sets)
      INTEGER(LONG)                   :: MFF_ROW_MAX_TERMS      ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: MFS_ROW_MAX_TERMS      ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
!xx   INTEGER(LONG)                   :: MSF_ROW_MAX_TERMS      ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: MSS_ROW_MAX_TERMS      ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG), PARAMETER        :: NUM1        = 1        ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2        ! Used in subr's that partition matrices
      INTEGER(LONG)                   :: NTERM_MSF              ! Number of nonzeros in sparse matrix MSF (should = NTERM_MFS)


      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Partition MFF from MNN. This is final MFF.

      IF (NDOFF > 0) THEN

         CALL PARTITION_SS_NTERM ( 'MNN', NTERM_MNN, NDOFN, NDOFN, SYM_MNN, I_MNN, J_MNN,      PART_VEC_N_FS, PART_VEC_N_FS,       &
                                    NUM1, NUM1, MFF_ROW_MAX_TERMS, 'MFF', NTERM_MFF, SYM_MFF )

         CALL ALLOCATE_SPARSE_MAT ( 'MFF', NDOFF, NTERM_MFF, SUBR_NAME )

         IF (NTERM_MFF > 0) THEN
            CALL PARTITION_SS ( 'MNN', NTERM_MNN, NDOFN, NDOFN, SYM_MNN, I_MNN, J_MNN, MNN, PART_VEC_N_FS, PART_VEC_N_FS,          &
                                 NUM1, NUM1, MFF_ROW_MAX_TERMS, 'MFF', NTERM_MFF, NDOFF, SYM_MFF, I_MFF, J_MFF, MFF )
         ENDIF

      ENDIF

! Partition MFS from MNN

      IF ((NDOFF > 0) .AND. (NDOFS > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'MNN', NTERM_MNN, NDOFN, NDOFN, SYM_MNN, I_MNN, J_MNN,      PART_VEC_N_FS, PART_VEC_N_FS,       &
                                    NUM1, NUM2, MFS_ROW_MAX_TERMS, 'MFS', NTERM_MFS, SYM_MFS )

         CALL ALLOCATE_SPARSE_MAT ( 'MFS', NDOFF, NTERM_MFS, SUBR_NAME )

         IF (NTERM_MFS > 0) THEN
            CALL PARTITION_SS ( 'MNN', NTERM_MNN, NDOFN, NDOFN, SYM_MNN, I_MNN, J_MNN, MNN, PART_VEC_N_FS, PART_VEC_N_FS,          &
                                 NUM1, NUM2, MFS_ROW_MAX_TERMS, 'MFS', NTERM_MFS, NDOFF, SYM_MFS, I_MFS, J_MFS, MFS )
         ENDIF

      ENDIF

! Partition MSS from MNN.

      IF (NDOFS > 0) THEN

         CALL PARTITION_SS_NTERM ( 'MNN', NTERM_MNN, NDOFN, NDOFN, SYM_MNN, I_MNN, J_MNN,      PART_VEC_N_FS, PART_VEC_N_FS,       &
                                    NUM2, NUM2, MSS_ROW_MAX_TERMS, 'MSS', NTERM_MSS, SYM_MSS )

         CALL ALLOCATE_SPARSE_MAT ( 'MSS', NDOFS, NTERM_MSS, SUBR_NAME )

         IF (NTERM_MFF > 0) THEN
            CALL PARTITION_SS ( 'MNN', NTERM_MNN, NDOFN, NDOFN, SYM_MNN, I_MNN, J_MNN, MNN, PART_VEC_N_FS, PART_VEC_N_FS,          &
                                 NUM2, NUM2, MSS_ROW_MAX_TERMS, 'MSS', NTERM_MSS, NDOFS, SYM_MSS, I_MSS, J_MSS, MSS )
         ENDIF

      ENDIF

! Partition MSF from MNN and write MSF to L2S for use in LINK9.

      IF ((NDOFF > 0) .AND. (NDOFS > 0)) THEN

!xx      CALL PARTITION_SS_NTERM ( 'MNN', NTERM_MNN, NDOFN, NDOFN, SYM_MNN, I_MNN, J_MNN,      PART_VEC_N_FS, PART_VEC_N_FS,       &
!xx                                 NUM2, NUM1, MSF_ROW_MAX_TERMS, 'MSF', NTERM_MSF, SYM_MFS )

!xx      IF (NTERM_MSF /= NTERM_MFS) THEN
!xx         FATAL_ERR = FATAL_ERR + 1
!xx         WRITE(ERR,936) SUBR_NAME, NTERM_MSF, NTERM_MFS
!xx         WRITE(F06,936) SUBR_NAME, NTERM_MSF, NTERM_MFS
!xx         CALL OUTA_HERE ( 'Y' )
!xx      ENDIF

!xx      CALL ALLOCATE_SPARSE_MAT ( 'MSF', NDOFS, NTERM_MSF, SUBR_NAME )

!xx      IF (NTERM_MSF > 0) THEN
!xx         CALL PARTITION_SS ( 'MNN', NTERM_MNN, NDOFN, NDOFN, SYM_MNN, I_MNN, J_MNN, MNN, PART_VEC_N_FS, PART_VEC_N_FS,          &
!xx                              NUM2, NUM1, MSF_ROW_MAX_TERMS, 'MSF', NTERM_MSF, NDOFS, SYM_MFS, I_MSF, J_MSF, MSF )

!xx         CLOSE_IT   = 'Y'
!xx         CLOSE_STAT = 'KEEP'
!xx         CALL WRITE_MATRIX_1 ( LINK2S, L2S, CLOSE_IT, CLOSE_STAT, L2S_MSG, 'MSF', NTERM_MSF, NDOFS, I_MSF, J_MSF, MSF )

!xx      ENDIF

         NTERM_MSF = NTERM_MFS
         CALL ALLOCATE_SPARSE_MAT ( 'MSF', NDOFS, NTERM_MSF, SUBR_NAME )

         IF (NTERM_MSF > 0) THEN
            CALL MATTRNSP_SS ( NDOFF, NDOFS, NTERM_MFS, 'MFS', I_MFS, J_MFS, MFS, 'MSF', I_MSF, J_MSF, MSF )

            CLOSE_IT   = 'Y'
            CLOSE_STAT = 'KEEP'
            CALL WRITE_MATRIX_1 ( LINK2S, L2S, CLOSE_IT, CLOSE_STAT, L2S_MSG, 'MSF', NTERM_MSF, NDOFS, I_MSF, J_MSF, MSF )

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
  936 FORMAT(' *ERROR   936: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF NONZERO TERMS IN SPARSE MATRIX MSF = ',I12,' IS NOT EQUAL TO THOSE IN MATRIX MFS = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_MNN_TO_MFF


      SUBROUTINE REDUCE_PN_TO_PF ( PART_VEC_N_FS, PART_VEC_SUB )

! Call routines to reduce the PN grid point load matrix from the N-set to the F, S-sets. See Appendix B to the MYSTRAN User's
! Reference Manual for the derivation of the reduction equations.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L2D, LINK2D, L2D_MSG, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFN, NDOFF, NDOFS, NDOFSE, NSUB, NTERM_KFSe, NTERM_PN,         &
                                         NTERM_PF, NTERM_PFYS, NTERM_PS
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE PARAMS, ONLY                :  MATSPARS
      USE SPARSE_MATRICES, ONLY       :  I_KFSe, J_KFSe, KFSe, I_PN, J_PN, PN, I_PF, J_PF, PF, I_PS, J_PS, PS, I_PF_TMP, J_PF_TMP, &
                                         PF_TMP, I_PFYS, J_PFYS, PFYS, I_PFYS1, J_PFYS1, PFYS1, I_QSYS, J_QSYS, QSYS
      USE SPARSE_MATRICES, ONLY       :  SYM_KFSe, SYM_PN, SYM_PF, SYM_PFYS, SYM_PF_TMP, SYM_PS
      USE COL_VECS, ONLY              :  YSe
      USE FULL_MATRICES, ONLY         :  KFSe_FULL, PF_FULL, PFYS_FULL, DUM1
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE SPARSE_FULL_MULTIPLICATION, ONLY:  MATMULT_SFS, MATMULT_SFS_NTERM
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM
      USE MATRIX_FILE_IO, ONLY        :  WRITE_MATRIX_1, WRITE_SPARSE_CRS
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  SPARSE_CRS_TO_FULL
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATADD_FFF, MATMULT_FFF
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_PN_TO_PF'
      CHARACTER(  1*BYTE)             :: CLOSE_IT            ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER(  8*BYTE)             :: CLOSE_STAT          ! Char constant for the CLOSE status of a file
      CHARACTER(132*BYTE)             :: MATRIX_NAME         ! Name of matrix for printout

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_N_FS(NDOFN)! Partitioning vector (F set into A and O sets)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_SUB(NSUB)  ! Partitioning vector (1's for all subcases)
      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG), PARAMETER        :: ITRNSPB     = 0     ! Transpose indicator for matrix multiply routine
      INTEGER(LONG)                   :: K                   ! Counter or DO loop index
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM_YS_COLS = 1     ! Variable for number of cols in array YSe
      INTEGER(LONG)                   :: NTERM_PF_TMP   = 0  ! No. of terms in matrix PF_TMP
      INTEGER(LONG)                   :: NTERM_PFYS1 = 0     ! No. of terms in matrix PFYS1
      INTEGER(LONG)                   :: PF_ROW_MAX_TERMS    ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: PS_ROW_MAX_TERMS    ! Output from subr PARTITION_SIZE (max terms in any row of matrix)


      REAL(DOUBLE)                    :: ALPHA               ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: BETA                ! Scalar multiplier for matrix

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Partition PF from PN (This is PF before reduction, or PF(bar) )

      IF (NDOFF > 0) THEN

         CALL PARTITION_SS_NTERM ( 'PN' , NTERM_PN , NDOFN, NSUB , SYM_PN , I_PN , J_PN ,      PART_VEC_N_FS, PART_VEC_SUB,        &
                                    NUM1, NUM1, PF_ROW_MAX_TERMS, 'PF', NTERM_PF, SYM_PF )

         CALL ALLOCATE_SPARSE_MAT ( 'PF', NDOFF, NTERM_PF, SUBR_NAME )

         IF (NTERM_PF  > 0) THEN
            CALL PARTITION_SS ( 'PN' , NTERM_PN , NDOFN, NSUB , SYM_PN , I_PN , J_PN , PN , PART_VEC_N_FS, PART_VEC_SUB,           &
                                 NUM1, NUM1, PF_ROW_MAX_TERMS, 'PF', NTERM_PF , NDOFF, SYM_PF, I_PF , J_PF , PF  )
         ENDIF

      ENDIF

! Partition PS from PN

      IF (NDOFS > 0) THEN

         CALL PARTITION_SS_NTERM ( 'PN' , NTERM_PN , NDOFN, NSUB , SYM_PN , I_PN , J_PN ,      PART_VEC_N_FS, PART_VEC_SUB,        &
                                    NUM2, NUM1, PS_ROW_MAX_TERMS, 'PS', NTERM_PS, SYM_PS )

         CALL ALLOCATE_SPARSE_MAT ( 'PS', NDOFS, NTERM_PS, SUBR_NAME )

         IF (NTERM_PS  > 0) THEN
            CALL PARTITION_SS ( 'PN' , NTERM_PN , NDOFN, NSUB , SYM_PN , I_PN , J_PN , PN , PART_VEC_N_FS, PART_VEC_SUB,           &
                                 NUM2, NUM1, PS_ROW_MAX_TERMS, 'PS', NTERM_PS , NDOFS, SYM_PS, I_PS , J_PS , PS  )
         ENDIF

      ENDIF

! Deallocate PN

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PN ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PN' )

! Reduce PN to PF = PF(bar) - PFYS where PF(bar) was original F-set partition from PN. If there are no enforced displs then
! PF = PF(bar) and no reduction is required. If there are enforced displs (NDOFSE is used to indicate if there are enforced
! displs) then reduction is required. So, when NDOFSE > 0 the reduction is done in 3 steps:

!     (1) calc PFYS1 = KFSe*YSe. This is 1 col since YSe is for all subcases.
!     (2) expand PFYS1 from 1 col to NSUB cols in PFYS with all cols = KFSe*YSe
!     (3) Calc reduced PF from PF(bar) and PFYS. This step is done 1 of 2 ways:
!         (a) If there are applied loads on the F-set DOF's (i.e. if NTERM_PF > 0), then we calc PF = PF(bar) - PFYS
!         (b) If there are no applied loads on the F-set DOF's (i.e. if NTERM_PF = 0), then we calc PF = -PFYS

! If PARAM MATSPARS = 'Y', then we use sparse matrix operations (multiply/add/transpose). If not, use full matrix operations

      IF (MATSPARS == 'Y') THEN

         IF ((NDOFF > 0) .AND. (NDOFSE > 0)) THEN          ! Do reduction only if there is an SE-set and an F-set
                                                           ! Step (1), calc PFYS1 = KFSe * YSe
            CALL MATMULT_SFS_NTERM ( 'KFSe', NDOFF, NTERM_KFSe, SYM_KFSe, I_KFSe, J_KFSe                                           &
                                     ,'YSe', NDOFSE, NUM_YS_COLS, YSe, AROW_MAX_TERMS, 'PFYS1', NTERM_PFYS1 )

            CALL ALLOCATE_SPARSE_MAT ( 'PFYS1', NDOFF, NTERM_PFYS1, SUBR_NAME )

            CALL MATMULT_SFS ( 'KFSe', NDOFF, NTERM_KFSe, SYM_KFSe, I_KFSe, J_KFSe, KFSe                                           &
                              ,'YSe', NDOFSE, NUM_YS_COLS, YSe, AROW_MAX_TERMS, 'PFYS1', ONE, NTERM_PFYS1, I_PFYS1, J_PFYS1, PFYS1 )


            NTERM_PFYS = NSUB*NTERM_PFYS1                  ! Step (2), expand PFS1 (1 col) into PFYS (NSUB cols, all = KFSe*YSe)

            CALL ALLOCATE_SPARSE_MAT ( 'PFYS', NDOFF, NTERM_PFYS, SUBR_NAME )

            I_PFYS(1) = 1
            DO I=2,NDOFF+1
               I_PFYS(I) = I_PFYS(I-1) + NSUB*(I_PFYS1(I) - I_PFYS1(I-1))
            ENDDO

            K = 0
            DO I=1,NTERM_PFYS1
               DO J=1,NSUB
                  K         = K + 1
                  J_PFYS(K) = J_PFYS1(I) + (J-1)
                    PFYS(K) = PFYS1(I)
               ENDDO
            ENDDO

      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PFYS1', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PFYS1' )


            IF (NTERM_PF > 0) THEN                         ! Step (3a), NTERM_PF > 0 so reduced PF = PF(bar) - PFYS

               IF (NTERM_PFYS > 0) THEN                    ! Only do step (3a) if there are nonzero's in PFYS, otherwise PF=PF(bar)
                  ALPHA =  ONE
                  BETA  = -ONE

                  CALL ALLOCATE_SPARSE_MAT ( 'PF_TMP', NDOFF, NTERM_PF, SUBR_NAME )

                  DO I=1,NDOFF+1                           ! PF_TMP = PF(bar) so we can use MATADD to get PF = PF_TMP - PFYS
                     I_PF_TMP(I) = I_PF(I)
                  ENDDO
                  NTERM_PF_TMP = NTERM_PF
                  DO K=1,NTERM_PF_TMP
                     J_PF_TMP(K) = J_PF(K)
                       PF_TMP(K) =   PF(K)
                  ENDDO
                                                           ! Recalc NTERM_PF for new PF = PF_TMP + PFYS
                  CALL MATADD_SSS_NTERM ( NDOFF, 'PF-bar', NTERM_PF_TMP, I_PF_TMP, J_PF_TMP, SYM_PF_TMP,                           &
                                                 'KFse*YSe', NTERM_PFYS  , I_PFYS  , J_PFYS  , SYM_PFYS  , 'PFYS', NTERM_PF )
            !xx   WRITE(SC1, * )                           ! Advance 1 line for screen messages
                  WRITE(SC1, * ) '    Reallocate PF '
            !xx   WRITE(SC1, * )                           ! Advance 1 line for screen messages
                  WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PF ', CR13
                  CALL DEALLOCATE_SPARSE_MAT ( 'PF' )
                  WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PF ', CR13
                  CALL ALLOCATE_SPARSE_MAT ( 'PF', NDOFF, NTERM_PF, SUBR_NAME )

                  CALL MATADD_SSS ( NDOFF, 'PF-bar', NTERM_PF_TMP, I_PF_TMP, J_PF_TMP, PF_TMP, ALPHA,                              &
                                   'KFse*YSe'  , NTERM_PFYS  , I_PFYS  , J_PFYS  , PFYS  , BETA ,'PF', NTERM_PF, I_PF, J_PF, PF)
            !xx   WRITE(SC1, * )                           ! Advance 1 line for screen messages
                  WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PF_TMP', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PF_TMP' )

               ENDIF

            ELSE                                           ! Step (3b), NTERM_PF = 0 so reduced PF = PF(bar)

               IF (NTERM_PFYS > 0) THEN                    ! Only do step (3b) if there are nonzero's in PFYS, otherwise PF = 0

                  NTERM_PF = NTERM_PFYS*NSUB
                  WRITE(SC1, * ) '    Reallocate PF '
            !xx   WRITE(SC1, * )                           ! Advance 1 line for screen messages
                  WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PF ', CR13
                  CALL DEALLOCATE_SPARSE_MAT ( 'PF' )
                  WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PF ', CR13
                  CALL ALLOCATE_SPARSE_MAT ( 'PF', NDOFF, NTERM_PF, SUBR_NAME )

                  DO I=1,NDOFF+1
                     I_PF(I) = I_PFYS(I)
                  ENDDO
                  K = 0
                  DO I=1,NTERM_PFYS
                     DO J=1,NSUB
                        K = K + 1
                        J_PF(K) = J_PFYS(I)
                          PF(K) =  -PFYS(I)
                     ENDDO
                  ENDDO


                  IF ((DEBUG(24) == 1) .OR.(DEBUG(24) == 3)) THEN
                     MATRIX_NAME = 'STIFFNESS PARTITION KFSe (partition of KFS for cols of enforced displs)   '
                     CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'F ', 'SE', NTERM_KFSe, NDOFF, I_KFSe, J_KFSe, KFSe )
                  ENDIF

                  IF ((DEBUG(25) == 1) .OR.(DEBUG(25) == 3)) THEN
                     MATRIX_NAME = 'LOAD MATRIX PFYS = KFSe*YSe (portion of F-set loads due to enforced displs)'
                     CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'F ', '  ', NTERM_PFYS, NDOFF, I_PFYS, J_PFYS, PFYS )
                  ENDIF

               ENDIF

            ENDIF

            WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PFYS', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PFYS' )
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KFSe', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KFSe' )

         ENDIF

      ELSE IF (MATSPARS == 'N') THEN
                                                           ! Allocate and set PF_FULL to sparse PF
         CALL ALLOCATE_FULL_MAT ( 'PF_FULL', NDOFF, NSUB, SUBR_NAME )

         IF (NTERM_PN > 0) THEN
            CALL SPARSE_CRS_TO_FULL ( 'PF', NTERM_Pf, NDOFf, NSUB, SYM_Pf, I_Pf, J_Pf, Pf, Pf_FULL )
         ENDIF

         IF ((NDOFF > 0) .AND. (NDOFSE > 0)) THEN          ! Do reduction only if there is an SE-set and an F-set

            CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFSE, NSUB, SUBR_NAME )
            DO I=1,NDOFSE
               DO J=1,NSUB
                  DUM1(I,J) = YSe(I)
               ENDDO
            ENDDO
            CALL ALLOCATE_FULL_MAT ( 'KFSe_FULL', NDOFF, NDOFSE, SUBR_NAME )
            CALL ALLOCATE_FULL_MAT ( 'PFYS_FULL', NDOFF, NSUB, SUBR_NAME )
            CALL MATMULT_FFF (KFSe_FULL, DUM1, NDOFF, NDOFSE, NSUB, PFYS_FULL )
            CALL DEALLOCATE_FULL_MAT ( 'DUM1' )

            IF (NTERM_PF > 0) THEN                         ! Step (3a), NTERM_PF > 0 so reduced PF = PF(bar) - PFYS

               IF (NTERM_PFYS > 0) THEN                    ! Only do step (3a) if there are nonzero's in PFYS, otherwise PF=PF(bar)

                  CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFF, NSUB, SUBR_NAME )

                  DO I=1,NDOFF
                     DO J=1,NSUB
                        DUM1(I,J) = PF_FULL(I,J)
                     ENDDO
                  ENDDO

                  ALPHA =  ONE
                  BETA  = -ONE
                  CALL MATADD_FFF ( DUM1, PFYS_FULL, NDOFF, NSUB, ALPHA, BETA, ITRNSPB, PF_FULL)

                  CALL DEALLOCATE_FULL_MAT ( 'DUM1' )

               ENDIF

            ELSE                                           ! Step (3b), NTERM_PF = 0 so reduced PF = PF(bar)

               IF (NTERM_PFYS > 0) THEN                    ! Only do step (3b) if there are nonzero's in PFYS, otherwise PF = 0

                  DO I=1,NDOFF
                     DO J=1,NSUB
                        PF_FULL(I,J) = -PFYS_FULL(I,J)
                     ENDDO
                  ENDDO

               ENDIF

            ENDIF

            CALL DEALLOCATE_FULL_MAT ( 'PFYS_FULL' )
            CALL DEALLOCATE_FULL_MAT ( 'KFSe_FULL' )

         ENDIF

         IF (NTERM_PF > 0) THEN                            ! Create new sparse arrays from PF_FULL
!xx         CALL FULL_TO_SPARSE_CRS ( 'PF_FULL   ', NDOFF, NSUB , PF_FULL,  NTERM_PF,  SYM_PF,  I_PF,  J_PF,  PF  )
            CALL DEALLOCATE_FULL_MAT ( 'PF_FULL' )
         ENDIF

      ELSE

         WRITE(ERR,911) SUBR_NAME,MATSPARS
         WRITE(F06,911) SUBR_NAME,MATSPARS
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF

! Write PS to L2D for constraint force recovery in LINK9.

      IF (NTERM_PS > 0) THEN
         CLOSE_IT   = 'Y'
         CLOSE_STAT = 'KEEP'
         CALL WRITE_MATRIX_1 ( LINK2D, L2D, CLOSE_IT, CLOSE_STAT, L2D_MSG, 'PS ', NTERM_PS , NDOFS, I_PS , J_PS , PS  )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  911 FORMAT(' *ERROR   911: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER MATSPARS MUST BE EITHER ','Y',' OR ','N',' BUT VALUE IS ',A)


12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_PN_TO_PF


   END MODULE REDUCTION_N_TO_F
