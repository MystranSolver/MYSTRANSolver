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

   MODULE OUTPUT4_FILE_IO

   USE SPARSE_CRS_ACCESS, ONLY :  SPARSE_MAT_DIAG_ZEROS
   USE SPARSE_FORMAT_CONVERSION, ONLY :  CRS_SYM_TO_CRS_NONSYM, SPARSE_CRS_SPARSE_CCS
   USE MATRIX_PARTITIONING, ONLY :  PARTITION_FF, PARTITION_SS, PARTITION_SS_NTERM
   USE SCRATCH_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
   USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
   USE FILE_LIFECYCLE, ONLY :  FILE_CLOSE, FILE_OPEN
   USE TEMP_FILE_WRITERS, ONLY :  WRITE_PARTND_MAT_HDRS
   USE DIAGNOSTICS_MEMORY_REPORTING, ONLY :  GET_OU4_MAT_STATS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: OUTPUT4_MATRIX_MSGS, OUTPUT4_PROC, READ_IN4_FULL_MAT, WRITE_OU4_FULL_MAT

   CONTAINS

      SUBROUTINE OUTPUT4_MATRIX_MSGS ( OUNT )

! Writes messages to F06 file regarding matrices that will be output in OUTOUT4 format

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, LEN_INPUT_FNAME, MOU4, OU4, OU4FIL, OU4_MSG, OU4STAT
      USE OUTPUT4_MATRICES, ONLY      :  ACT_OU4_MYSTRAN_NAMES, ACT_OU4_MYSTRAN_NAMES, ACT_OU4_OUTPUT_NAMES,                       &
                                         ALLOW_OU4_MYSTRAN_NAMES, ALLOW_OU4_OUTPUT_NAMES, NUM_OU4_REQUESTS, NUM_OU4_VALID_NAMES,   &
                                         OU4_FILE_UNITS
      USE PARAMS, ONLY                :  SUPINFO

      USE FILE_LIFECYCLE, ONLY: FILE_OPEN
      USE FILE_LIFECYCLE, ONLY: FILE_CLOSE

      IMPLICIT NONE

      CHARACTER(LEN(ALLOW_OU4_MYSTRAN_NAMES))                                                                                      &
                                      :: MYSTRAN_NAMES_I(NUM_OU4_VALID_NAMES)
                                                                     ! OUTPUT4 matrix names written to one OU4 unit

      CHARACTER(LEN(ALLOW_OU4_MYSTRAN_NAMES))                                                                                      &
                                      :: OUTPUT_NAMES_I(NUM_OU4_VALID_NAMES)
                                                                     ! OUTPUT4 matrix names written to one OU4 unit

      INTEGER(LONG), INTENT(IN)       :: OUNT(2)             ! File units to write messages to. Input to subr UNFORMATTED_OPEN.
      INTEGER(LONG)                   :: I,J                         ! DO loop indices
      INTEGER(LONG)                   :: NOU4_UNITS                  ! Number of OU4 units for files requested
      INTEGER(LONG)                   :: NOU4_FILES                  ! Number of matrices to write to a specific OPi file

! **********************************************************************************************************************************
      DO I=1,MOU4                                          ! Initially set all to close status of 'DELETE'
         OU4STAT(I) = 'DELETE'                             ! (should have been done in module IOUNT1 but do it here to make sure)
      ENDDO

      NOU4_UNITS = 0                                       ! Reset to "KEEP' close status of units where OUTPUT4 files written
      DO I=1,NUM_OU4_REQUESTS
         DO J=1,MOU4
            IF (OU4_FILE_UNITS(I) == OU4(J)) THEN
               OU4STAT(J) = 'KEEP'
               EXIT
            ENDIF
         ENDDO
      ENDDO

      NOU4_UNITS = 0                                       ! Count the number of OUTPUT4 files that need to be kept
      DO I=1,MOU4
         IF (OU4STAT(I) == 'KEEP') THEN
            NOU4_UNITS = NOU4_UNITS + 1
         ENDIF
      ENDDO

      IF (NOU4_UNITS > 0) THEN
         WRITE(F06,*)
         WRITE(ERR,288) NUM_OU4_REQUESTS, NOU4_UNITS
!xx      IF (SUPINFO == 'N') THEN
            WRITE(F06,288) NUM_OU4_REQUESTS, NOU4_UNITS
!xx      ENDIF
         DO I=1,MOU4                                       ! Open and then close as "KEEP' the requested OUTPUT4 files
            NOU4_FILES = 0
            IF (OU4STAT(I) == 'KEEP') THEN
               WRITE(F06,289) OU4(I), OU4FIL(I)(1:LEN_INPUT_FNAME+3)
               DO J=1,NUM_OU4_REQUESTS
                  IF (OU4_FILE_UNITS(J) == OU4(I)) THEN
                     NOU4_FILES = NOU4_FILES + 1
                     MYSTRAN_NAMES_I(NOU4_FILES) = ACT_OU4_MYSTRAN_NAMES(J)
                     OUTPUT_NAMES_I(NOU4_FILES)  = ACT_OU4_OUTPUT_NAMES(J)
                  ENDIF
               ENDDO
               DO J=1,NOU4_FILES
                  IF (OUTPUT_NAMES_I(J) == MYSTRAN_NAMES_I(J)) THEN
                     WRITE(F06,290) J, OUTPUT_NAMES_I(J)
                  ELSE
                     WRITE(F06,291) J, OUTPUT_NAMES_I(J), MYSTRAN_NAMES_I(J)
                  ENDIF
               ENDDO
               WRITE(F06,*)
               CALL FILE_OPEN  ( OU4(I), OU4FIL(I), OUNT, 'REPLACE', OU4_MSG(I), 'NEITHER', 'UNFORMATTED', 'WRITE', 'REWIND',      &
                                 'Y', 'N')
               CALL FILE_CLOSE ( OU4(I), OU4FIL(I), 'KEEP' )
            ENDIF
         ENDDO
         WRITE(F06,*)
      ENDIF

! **********************************************************************************************************************************
  288 FORMAT(' *INFORMATION: THE FOLLOWING ',I3,' MATRICES HAVE BEEN REQUESTED TO BE WRITTEN TO ',I2,' OUTPUT4 FILES IN THE ORDER',&
                           ' LISTED BELOW:',/)

  289 FORMAT(14X,' OUTPUT4 file on unit ',I3,' has been created as: ',A,' and will contain the matrices:')

  290 FORMAT(23X,'(',I2,') ',A)

  291 FORMAT(23X,'(',I2,') ',A,3X,': this is MYSTRAN matrix ',A)

! **********************************************************************************************************************************

      END SUBROUTINE OUTPUT4_MATRIX_MSGS



      SUBROUTINE OUTPUT4_PROC ( CALLING_SUBR )

! Checks whether a matrix is requested for OUTPUT4 and, if so:

!   - calls OU4_PARTVEC_PROC to calculate the row/col partitioning vectors (OU4_PARTVEC_ROW, OU4_PARTVEC_COL)
!   - calls PARTITION_SS to do the actual partitioning
!   - calls WRITE_OU4_SPARSE_MAT or WRITE_OU4_FULL_MAT to write the matrix to an unformatted disk file

! This subr does not process the grid and/or element related Output Transformation Matrices (OTM's). That is done in LINK9

      USE PENTIUM_II_KIND, ONLY        :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                 :  ERR, F06, MOU4, OU4, OU4_MSG, OU4FIL

      USE SCONTR, ONLY                 :  BLNK_SUB_NAM, FATAL_ERR   ,                                                              &
                                          NTERM_CG_LTM, NTERM_DLR   , NTERM_IF_LTM, NTERM_KLL   ,                                  &
                                          NTERM_KRL   , NTERM_KRR   , NTERM_KRRcb , NTERM_KXX   ,                                  &
                                          NTERM_LTM   ,                                                                            &
                                          NTERM_MLL   , NTERM_MRL   , NTERM_MRN   , NTERM_MRR   ,                                  &
                                          NTERM_MRRcb , NTERM_MXX   , NTERM_PHIXG , NTERM_PHIZL ,                                  &
                                          NTERM_KAA, NTERM_MAA, NTERM_KGG, NTERM_MGG, NTERM_PA, NTERM_PG, NTERM_PL, SOL_NAME

      USE SCONTR, ONLY                 :  NTERM_KAA

      USE MODEL_STUF, ONLY             :  MCG
      USE PARAMS, ONLY                 :  MPFOUT

      USE EIGEN_MATRICES_1, ONLY       :  EIGEN_VAL, EIGEN_VEC, GEN_MASS,  MEFFMASS,  MPFACTOR_N6,  MPFACTOR_NR

      USE FULL_MATRICES, ONLY          :  DUM1, PHIZG_FULL

      USE MODEL_STUF, ONLY             :  MCG

      USE RIGID_BODY_DISP_MATS, ONLY   :  TR6_0, TR6_CG

      USE SPARSE_MATRICES, ONLY        :  I_CG_LTM  ,I_DLR     ,I_IF_LTM  ,I_KAA     ,I_KGG     ,I_KLL     ,I_KRL     ,I_KRR     , &
                                          I_KRRcb   ,I_KXX     ,I_LTM     ,I_MAA     ,I_MGG     ,I_MLL     ,I_MRL     ,I_MRN     , &
                                          I_MRR     ,I_MRRcb   ,I_MXX     ,I_PA      ,I_PG      ,I_PL      ,I_PHIXG

      USE SPARSE_MATRICES, ONLY        :  J_CG_LTM  ,J_DLR     ,J_IF_LTM  ,J_KAA     ,J_KGG     ,J_KLL     ,J_KRL     ,J_KRR     , &
                                          J_KRRcb   ,J_KXX     ,J_LTM     ,J_MAA     ,J_MGG     ,J_MLL     ,J_MRL     ,J_MRN     , &
                                          J_MRR     ,J_MRRcb   ,J_MXX     ,J_PA      ,J_PG      ,J_PL      ,J_PHIXG

      USE SPARSE_MATRICES, ONLY        :  CG_LTM    ,DLR       ,IF_LTM    ,KAA       ,KGG       ,KLL       ,KRL       ,KRR       , &
                                          KRRcb     ,KXX       ,LTM       ,MAA       ,MGG       ,MLL       ,MRL       ,MRN       , &
                                          MRR       ,MRRcb     ,MXX       ,PA        ,PG        ,PL        ,PHIXG

      USE SPARSE_MATRICES, ONLY        :  SYM_CG_LTM,SYM_DLR   ,SYM_IF_LTM,SYM_KAA   ,SYM_KGG   ,SYM_KLL   ,SYM_KRL   ,SYM_KRR   , &
                                          SYM_KRRcb ,SYM_KXX   ,SYM_LTM   ,SYM_MAA   ,SYM_MGG   ,SYM_MLL   ,SYM_MRL   ,SYM_MRN   , &
                                          SYM_MRR   ,SYM_MRRcb ,SYM_MXX   ,SYM_PA    ,SYM_PG    ,SYM_PL    ,SYM_PHIXG

      USE SCRATCH_MATRICES, ONLY       :  I_CRS1, J_CRS1, CRS1

      USE OUTPUT4_MATRICES, ONLY       :  ACT_OU4_MYSTRAN_NAMES, HAS_OU4_MAT_BEEN_PROCESSED, NUM_OU4_REQUESTS, OU4_FILE_UNITS,     &
                                          OU4_PARTVEC_COL, OU4_PARTVEC_ROW,                                                        &
                                          OU4_PART_MAT_NAMES, OU4_PART_VEC_NAMES, RBM0, SUBR_WHEN_TO_WRITE_OU4_MATS

      USE TIMDAT, ONLY                 :  TSEC
      USE RIGID_BODY_DISP_MATS, ONLY   :  TR6_CG, TR6_0
      USE EIGEN_MATRICES_1, ONLY       :  GEN_MASS, EIGEN_VAL, EIGEN_VEC, MEFFMASS, MPFACTOR_N6, MPFACTOR_NR

      USE SPARSE_MATRICES, ONLY        :  I_CG_LTM, J_CG_LTM, CG_LTM,      I_DLR   , J_DLR   , DLR   ,                             &
                                          I_IF_LTM, J_IF_LTM, IF_LTM,      I_KLL   , J_KLL   , KLL   ,                             &
                                          I_KRL   , J_KRL   , KRL   ,      I_KRR   , J_KRR   , KRR   ,                             &
                                          I_KRRcb , J_KRRcb , KRRcb ,      I_KXX   , J_KXX   , KXX   ,                             &
                                          I_LTM   , J_LTM   , LTM   ,                                                              &
                                          I_MLL   , J_MLL   , MLL   ,      I_MRL   , J_MRL   , MRL   ,                             &
                                          I_MRN   , J_MRN   , MRN   ,      I_MRR   , J_MRR   , MRR   ,                             &
                                          I_MRRcb , J_MRRcb , MRRcb ,      I_MXX   , J_MXX   , MXX   ,                             &
                                          I_PHIXG , J_PHIXG , PHIXG

      USE SPARSE_MATRICES, ONLY        :  I_KAA, J_KAA, KAA, I_KGG, J_KGG, KGG, I_MAA, J_MAA, MAA, I_MGG, J_MGG, MGG,              &
                                          I_PA , J_PA , PA , I_PG , J_PG , PG , I_PL , J_PL , PL

      USE FULL_MATRICES, ONLY          :  PHIZG_FULL

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_OU4_MAT_STATS
      USE OUTPUT4_PARTITIONING, ONLY  :  OU4_PARTVEC_PROC
      USE TEMP_FILE_WRITERS, ONLY: WRITE_PARTND_MAT_HDRS
      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_FF, PARTITION_SS, PARTITION_SS_NTERM
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE


      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OUTPUT4_PROC'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this one
      CHARACTER( 1*BYTE)              :: CAN_PARTN         ! 'Y' if matrix can be partitioned (rows, cols or both)

      CHARACTER(LEN(ACT_OU4_MYSTRAN_NAMES))                                                                                        &
                                      :: MAT_NAME          ! Name of matrix requested for OUTPUT4 for this call to OUTPUT4_PROC

      CHARACTER(1*BYTE)               :: SM                ! 'Y' if OUTPUT4 matrix is stored symmetrically

      INTEGER(LONG)                   :: AROW_MAX_TERMS    ! Max number of terms in any row of partitioned matrix
      INTEGER(LONG)                   :: FORM              ! Format of OUTPUT4 matrix
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: NCOLS_F           ! Number of cols in the complete OUTPUT4 matrix
      INTEGER(LONG)                   :: NROWS_F           ! Number of cols in the complete OUTPUT4 matrix
      INTEGER(LONG)                   :: NCOLS_P           ! Number of cols in that will be in the OUTPUT4 matrix when partitioned
      INTEGER(LONG)                   :: NROWS_P           ! Number of cols in that will be in the OUTPUT4 matrix when partitioned
      INTEGER(LONG)                   :: NTERM_CRS1        ! Number of terms in partitioned sparse matrix CRS1
      INTEGER(LONG)                   :: VAL_C             ! Non-zero vals in OU4_PARTVEC_ROWS
      INTEGER(LONG)                   :: VAL_R             ! Non-zero vals in OU4_PARTVEC_COLS
      INTEGER(LONG)                   :: UNT               !




! **********************************************************************************************************************************
      DO I=1,NUM_OU4_REQUESTS

         MAT_NAME = ACT_OU4_MYSTRAN_NAMES(I)
         UNT      = OU4_FILE_UNITS(I)

         IF      (MAT_NAME == 'CG_LTM          ') THEN     ! ( 1)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS( 1) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(CG_LTM)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, '--', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, '--', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'CG_LTM', NTERM_CG_LTM, NROWS_F, NCOLS_F, 'N', I_CG_LTM, J_CG_LTM,                &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'CG_LTM', NTERM_CG_LTM, NROWS_F, NCOLS_F, 'N', I_CG_LTM, J_CG_LTM, CG_LTM,        &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_CG_LTM, I_CG_LTM, J_CG_LTM, CG_LTM, UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'DLR             ') THEN     ! ( 2)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS( 2) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(DLR)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'L ', 'R ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'L ', 'R ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'DLR', NTERM_DLR, NROWS_F, NCOLS_F, 'N', I_DLR, J_DLR,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'DLR', NTERM_DLR, NROWS_F, NCOLS_F, 'N', I_DLR, J_DLR, DLR,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_DLR   , I_DLR   , J_DLR   , DLR   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'EIGEN_VAL       ') THEN     ! ( 3)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS( 3) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(EIGEN_VAL)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, '--', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, '--', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL ALLOCATE_FULL_MAT ( 'DUM1', NROWS_P, NCOLS_P, SUBR_NAME )
                     CALL PARTITION_FF ( 'EIGEN_VAL', NROWS_F, NCOLS_F, EIGEN_VAL, OU4_PARTVEC_ROW, OU4_PARTVEC_COL,             &
                                          VAL_R, VAL_C, 'DUM1', NROWS_P, NCOLS_P, DUM1 )
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, EIGEN_VAL , UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, EIGEN_VAL, UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'EIGEN_VEC       ') THEN     ! ( 4)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS( 4) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(EIGEN_VEC)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'G ', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'G ', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL ALLOCATE_FULL_MAT ( 'DUM1', NROWS_P, NCOLS_P, SUBR_NAME )
                     CALL PARTITION_FF ( 'EIGEN_VEC', NROWS_F, NCOLS_F, EIGEN_VEC, OU4_PARTVEC_ROW, OU4_PARTVEC_COL,             &
                                          VAL_R, VAL_C, 'DUM1', NROWS_P, NCOLS_P, DUM1 )
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, EIGEN_VEC , UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, EIGEN_VEC, UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'GEN_MASS        ') THEN     ! ( 5)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS( 5) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(GEN_MASS)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, '--', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, '--', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL ALLOCATE_FULL_MAT ( 'DUM1', NROWS_P, NCOLS_P, SUBR_NAME )
                     CALL PARTITION_FF ( 'GEN_MASS', NROWS_F, NCOLS_F, GEN_MASS, OU4_PARTVEC_ROW, OU4_PARTVEC_COL,                 &
                                          VAL_R, VAL_C, 'DUM1', NROWS_P, NCOLS_P, DUM1 )
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, GEN_MASS , UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, GEN_MASS , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'IF_LTM          ') THEN     ! ( 6)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS( 6) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(IF_LTM)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'R ', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'R ', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'IF_LTM', NTERM_IF_LTM, NROWS_F, NCOLS_F, 'N', I_IF_LTM, J_IF_LTM,                &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'IF_LTM', NTERM_IF_LTM, NROWS_F, NCOLS_F, 'N', I_IF_LTM, J_IF_LTM, IF_LTM,        &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_IF_LTM, I_IF_LTM, J_IF_LTM, IF_LTM, UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'KAA             ') THEN     ! ( 7)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS( 7) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(KAA)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'A ', 'A ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'A ', 'A ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'KAA', NTERM_KAA, NROWS_F, NCOLS_F, 'N', I_KAA, J_KAA,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'KAA', NTERM_KAA, NROWS_F, NCOLS_F, 'N', I_KAA, J_KAA, KAA,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_KAA, I_KAA, J_KAA, KAA, UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'KGG             ') THEN     ! ( 8)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS( 8) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(KGG)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'G ', 'G ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'G ', 'G ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'KGG', NTERM_KGG, NROWS_F, NCOLS_F, 'N', I_KGG, J_KGG,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'KGG', NTERM_KGG, NROWS_F, NCOLS_F, 'N', I_KGG, J_KGG, KGG,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_KGG   , I_KGG   , J_KGG   , KGG   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'KLL             ') THEN     ! ( 9)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS( 9) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(KLL)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'L ', 'L ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'L ', 'L ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'KLL', NTERM_KLL, NROWS_F, NCOLS_F, 'N', I_KLL, J_KLL,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'KLL', NTERM_KLL, NROWS_F, NCOLS_F, 'N', I_KLL, J_KLL, KLL,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_KLL   , I_KLL   , J_KLL   , KLL   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'KRL             ') THEN     ! (10)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(10) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(KRL)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'R ', 'L ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'R ', 'L ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'KRL', NTERM_KRL, NROWS_F, NCOLS_F, 'N', I_KRL, J_KRL,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'KRL', NTERM_KRL, NROWS_F, NCOLS_F, 'N', I_KRL, J_KRL, KRL,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_KRL   , I_KRL   , J_KRL   , KRL   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'KRR             ') THEN     ! (11)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(11) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(KRR)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'R ', 'R ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'R ', 'R ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'KRR', NTERM_KRR, NROWS_F, NCOLS_F, 'N', I_KRR, J_KRR,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'KRR', NTERM_KRR, NROWS_F, NCOLS_F, 'N', I_KRR, J_KRR, KRR,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_KRR   , I_KRR   , J_KRR   , KRR   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'KRRcb           ') THEN     ! (12)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(12) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(KRRcb)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'R ', 'R ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'R ', 'R ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'KRRcb', NTERM_KRRcb, NROWS_F, NCOLS_F, 'N', I_KRRcb, J_KRRcb,                    &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'KRRcb', NTERM_KRRcb, NROWS_F, NCOLS_F, 'N', I_KRRcb, J_KRRcb, KRRcb,             &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_KRRcb , I_KRRcb , J_KRRcb , KRRcb , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'KXX             ') THEN     ! (13)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(13) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(KXX)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, '--', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, '--', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'KXX', NTERM_KXX, NROWS_F, NCOLS_F, 'N', I_KXX, J_KXX,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'KXX', NTERM_KXX, NROWS_F, NCOLS_F, 'N', I_KXX, J_KXX, KXX,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_KXX   , I_KXX   , J_KXX   , KXX   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'LTM             ') THEN     ! (14)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(14) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(LTM)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, '--', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, '--', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'LTM', NTERM_LTM, NROWS_F, NCOLS_F, 'N', I_LTM, J_LTM,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'LTM', NTERM_LTM, NROWS_F, NCOLS_F, 'N', I_LTM, J_LTM, LTM,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_LTM   , I_LTM   , J_LTM   , LTM   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'MCG             ') THEN     ! (15)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(15) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                  WRITE(ERR,101) MAT_NAME
                  WRITE(F06,101) MAT_NAME
               ENDIF
               NROWS_P = NROWS_F
               NCOLS_P = NCOLS_F
               CALL WRITE_OU4_FULL_MAT   ( MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, MCG      , UNT)
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
               HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
            ENDIF

         ELSE IF (MAT_NAME == 'MEFFMASS        ') THEN     ! (16)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(16) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(MEFFMASS)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, '--', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, '--', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL ALLOCATE_FULL_MAT ( 'DUM1', NROWS_P, NCOLS_P, SUBR_NAME )
                     CALL PARTITION_FF ( 'MEFFMASS', NROWS_F, NCOLS_F, MEFFMASS, OU4_PARTVEC_ROW, OU4_PARTVEC_COL,                 &
                                          VAL_R, VAL_C, 'DUM1', NROWS_P, NCOLS_P, DUM1 )
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, MEFFMASS , UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, MEFFMASS, UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'MPFACTOR        ') THEN     ! (17)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(17) == CALLING_SUBR) THEN
               IF (MPFOUT == '6') THEN
                  CALL GET_OU4_MAT_STATS    ('MPFACTOR_N6', NROWS_F, NCOLS_F, FORM, SM )
                  IF (ALLOCATED(MPFACTOR_N6)) THEN
                     CAN_PARTN = 'N'
                     IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                        CALL OU4_PARTVEC_PROC (I, MAT_NAME, NROWS_F, NCOLS_F, '--', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C)
                        CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, '--', '--', NROWS_P, NCOLS_P )
                     ELSE
                        NROWS_P = NROWS_F
                        NCOLS_P = NCOLS_F
                     ENDIF
                     IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                        CALL ALLOCATE_FULL_MAT ( 'DUM1', NROWS_P, NCOLS_P, SUBR_NAME )
                        CALL PARTITION_FF ( 'MPFACTOR_N6', NROWS_F, NCOLS_F, MPFACTOR_N6, OU4_PARTVEC_ROW, OU4_PARTVEC_COL,        &
                                             VAL_R, VAL_C, 'DUM1', NROWS_P, NCOLS_P, DUM1 )
                        CALL WRITE_OU4_FULL_MAT ( 'MPFACTOR_N6', NROWS_P, NCOLS_P, FORM, SM, MPFACTOR_N6 , UNT)
                     ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                        CALL WRITE_OU4_FULL_MAT ( 'MPFACTOR_N6', NROWS_F, NCOLS_F, FORM, SM, MPFACTOR_N6, UNT)
                     ENDIF
                     HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
                  ELSE
                     IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                        WRITE(ERR,201) MAT_NAME, SOL_NAME
                        WRITE(F06,201) MAT_NAME, SOL_NAME
                     ENDIF
                     HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
                  ENDIF
               ELSE
                  CALL GET_OU4_MAT_STATS    ('MPFACTOR_NR', NROWS_F, NCOLS_F, FORM, SM )
                  IF (ALLOCATED(MPFACTOR_NR)) THEN
                     CAN_PARTN = 'N'
                     IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                        CALL OU4_PARTVEC_PROC (I, MAT_NAME, NROWS_F, NCOLS_F, '--', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C)
                        CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, '--', '--', NROWS_P, NCOLS_P )
                     ELSE
                        NROWS_P = NROWS_F
                        NCOLS_P = NCOLS_F
                     ENDIF
                     IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                        CALL ALLOCATE_FULL_MAT ( 'DUM1', NROWS_P, NCOLS_P, SUBR_NAME )
                        CALL PARTITION_FF ( 'MPFACTOR_NR', NROWS_F, NCOLS_F, MPFACTOR_NR, OU4_PARTVEC_ROW, OU4_PARTVEC_COL,        &
                                             VAL_R, VAL_C, 'DUM1', NROWS_P, NCOLS_P, DUM1 )
                        CALL WRITE_OU4_FULL_MAT ( 'MPFACTOR_NR', NROWS_P, NCOLS_P, FORM, SM, MPFACTOR_NR , UNT)
                     ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                        CALL WRITE_OU4_FULL_MAT ( 'MPFACTOR_NR', NROWS_F, NCOLS_F, FORM, SM, MPFACTOR_NR, UNT)
                     ENDIF
                     HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
                  ELSE
                     IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                        WRITE(ERR,201) MAT_NAME, SOL_NAME
                        WRITE(F06,201) MAT_NAME, SOL_NAME
                     ENDIF
                     HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
                  ENDIF
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'MAA             ') THEN     ! (18)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(18) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(MAA)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'A ', 'A ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'A ', 'A ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'MAA', NTERM_MAA, NROWS_F, NCOLS_F, 'N', I_MAA, J_MAA,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'MAA', NTERM_MAA, NROWS_F, NCOLS_F, 'N', I_MAA, J_MAA, MAA,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_MAA   , I_MAA   , J_MAA   , MAA   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'MGG             ') THEN     ! (19)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(19) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(MGG)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'G ', 'G ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'G ', 'G ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'MGG', NTERM_MGG, NROWS_F, NCOLS_F, 'N', I_MGG, J_MGG,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'MGG', NTERM_MGG, NROWS_F, NCOLS_F, 'N', I_MGG, J_MGG, MGG,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_MGG   , I_MGG   , J_MGG   , MGG   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'MLL             ') THEN     ! (20)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(20) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(MLL)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'L ', 'L ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'L ', 'L ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'MLL', NTERM_MLL, NROWS_F, NCOLS_F, 'N', I_MLL, J_MLL,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'MLL', NTERM_MLL, NROWS_F, NCOLS_F, 'N', I_MLL, J_MLL, MLL,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_MLL   , I_MLL   , J_MLL   , MLL   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'MRL             ') THEN     ! (21)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(21) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(MRL)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'R ', 'L ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'R ', 'L ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'MRL', NTERM_MRL, NROWS_F, NCOLS_F, 'N', I_MRL, J_MRL,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'MRL', NTERM_MRL, NROWS_F, NCOLS_F, 'N', I_MRL, J_MRL, MRL,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_MRL   , I_MRL   , J_MRL   , MRL   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'MRN             ') THEN     ! (22)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(22) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(MRN)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'R ', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'R ', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'MRN', NTERM_MRN, NROWS_F, NCOLS_F, 'N', I_MRN, J_MRN,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'MRN', NTERM_MRN, NROWS_F, NCOLS_F, 'N', I_MRN, J_MRN, MRN,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_MRN   , I_MRN   , J_MRN   , MRN   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'MRR             ') THEN     ! (23)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(23) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(MRR)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'R ', 'R ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'R ', 'R ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'MRR', NTERM_MRR, NROWS_F, NCOLS_F, 'N', I_MRR, J_MRR,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'MRR', NTERM_MRR, NROWS_F, NCOLS_F, 'N', I_MRR, J_MRR, MRR,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_MRR   , I_MRR   , J_MRR   , MRR   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'MRRcb           ') THEN     ! (24)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(24) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(MRRcb)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'R ', 'R ', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'R ', 'R ', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'MRRcb', NTERM_MRRcb, NROWS_F, NCOLS_F, 'N', I_MRRcb, J_MRRcb,                    &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'MRRcb', NTERM_MRRcb, NROWS_F, NCOLS_F, 'N', I_MRRcb, J_MRRcb, MRRcb,             &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_MRRcb , I_MRRcb , J_MRRcb , MRRcb , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'MXX             ') THEN     ! (25)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(25) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(MXX)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, '--', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, '--', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'MXX', NTERM_MXX, NROWS_F, NCOLS_F, 'N', I_MXX, J_MXX,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'MXX', NTERM_MXX, NROWS_F, NCOLS_F, 'N', I_MXX, J_MXX, MXX,                       &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_MXX   , I_MXX   , J_MXX   , MXX   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'PA              ') THEN     ! (26)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(26) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(PA)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'A ', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'A ', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'PA', NTERM_PA, NROWS_F, NCOLS_F, 'N', I_PA, J_PA,                                &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'PA', NTERM_PA, NROWS_F, NCOLS_F, 'N', I_PA, J_PA, PA,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_PA    , I_PA    , J_PA    , PA    , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'PG              ') THEN     ! (27)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(27) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(PG)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'G ', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'G ', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'PG', NTERM_PG, NROWS_F, NCOLS_F, 'N', I_PG, J_PG,                                &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'PG', NTERM_PG, NROWS_F, NCOLS_F, 'N', I_PG, J_PG, PG,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_PG    , I_PG    , J_PG    , PG    , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'PL              ') THEN     ! (28)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(28) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(PL)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'L ', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'L ', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'PL', NTERM_PL, NROWS_F, NCOLS_F, 'N', I_PL, J_PL,                                &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'PL', NTERM_PL, NROWS_F, NCOLS_F, 'N', I_PL, J_PL, PL,                            &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_PL    , I_PL    , J_PL    , PL    , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'PHIXG           ') THEN     ! (29)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(29) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(PHIXG)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'G ', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'G ', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL PARTITION_SS_NTERM   ( 'PHIXG', NTERM_PHIXG, NROWS_F, NCOLS_F, 'N', I_PHIXG, J_PHIXG,                    &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, 'N' )

                     CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS_P, NTERM_CRS1, SUBR_NAME )

                     CALL PARTITION_SS         ( 'PHIXG', NTERM_PHIXG, NROWS_F, NCOLS_F, 'N', I_PHIXG, J_PHIXG, PHIXG,             &
                                                  OU4_PARTVEC_ROW, OU4_PARTVEC_COL, VAL_R, VAL_C, AROW_MAX_TERMS,                  &
                                                 'CRS1', NTERM_CRS1, NROWS_P, 'N', &
                                                  I_CRS1, J_CRS1, CRS1 )

                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, NTERM_CRS1, I_CRS1, J_CRS1, CRS1, UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_SPARSE_MAT (MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, NTERM_PHIXG , I_PHIXG , J_PHIXG , PHIXG , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'PHIZG           ') THEN     ! (30)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(30) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(PHIZG_FULL)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'G ', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'G ', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL ALLOCATE_FULL_MAT ( 'DUM1', NROWS_P, NCOLS_P, SUBR_NAME )
                     CALL PARTITION_FF ( 'PHIZG_FULL', NROWS_F, NCOLS_F, PHIZG_FULL, OU4_PARTVEC_ROW, OU4_PARTVEC_COL,             &
                                          VAL_R, VAL_C, 'DUM1', NROWS_P, NCOLS_P, DUM1 )
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, PHIZG_FULL , UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, PHIZG_FULL , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'RBM0            ') THEN     ! (31)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(31) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                  WRITE(ERR,101) MAT_NAME
                  WRITE(F06,101) MAT_NAME
               ENDIF
               NROWS_P = NROWS_F
               NCOLS_P = NCOLS_F
               CALL WRITE_OU4_FULL_MAT   ( MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, RBM0     , UNT)
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
               HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
            ENDIF

         ELSE IF (MAT_NAME == 'TR6_0           ') THEN     ! (32)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(32) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(TR6_0)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'R ', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'R ', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL ALLOCATE_FULL_MAT ( 'DUM1', NROWS_P, NCOLS_P, SUBR_NAME )
                     CALL PARTITION_FF ( 'TR6_0', NROWS_F, NCOLS_F, TR6_0, OU4_PARTVEC_ROW, OU4_PARTVEC_COL,                       &
                                          VAL_R, VAL_C, 'DUM1', NROWS_P, NCOLS_P, DUM1 )
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, TR6_0 , UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, TR6_0    , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE IF (MAT_NAME == 'TR6_CG          ') THEN     ! (33)
            IF (SUBR_WHEN_TO_WRITE_OU4_MATS(33) == CALLING_SUBR) THEN
               CALL GET_OU4_MAT_STATS    ( MAT_NAME    , NROWS_F, NCOLS_F, FORM, SM )
               IF (ALLOCATED(TR6_CG)) THEN
                  CAN_PARTN = 'N'
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     CALL OU4_PARTVEC_PROC ( I, MAT_NAME, NROWS_F, NCOLS_F, 'R ', '--', CAN_PARTN, NROWS_P, NCOLS_P, VAL_R, VAL_C )
                     CALL WRITE_PARTNd_MAT_HDRS ( MAT_NAME, 'R ', '--', NROWS_P, NCOLS_P )
                  ELSE
                     NROWS_P = NROWS_F
                     NCOLS_P = NCOLS_F
                  ENDIF
                  IF (CAN_PARTN == 'Y') THEN               ! Partition the matrix and write it out to the OU4 file
                     CALL ALLOCATE_FULL_MAT ( 'DUM1', NROWS_P, NCOLS_P, SUBR_NAME )
                     CALL PARTITION_FF ( 'TR6_CG', NROWS_F, NCOLS_F, TR6_CG, OU4_PARTVEC_ROW, OU4_PARTVEC_COL,                     &
                                          VAL_R, VAL_C, 'DUM1', NROWS_P, NCOLS_P, DUM1 )
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_P, NCOLS_P, FORM, SM, TR6_CG , UNT)
                  ELSE                                     ! Matrix cannot be partitioned so just write it out to the OU4 file
                     CALL WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS_F, NCOLS_F, FORM, SM, TR6_CG   , UNT)
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'Y'
               ELSE
                  IF (OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') THEN
                     WRITE(ERR,201) MAT_NAME, SOL_NAME
                     WRITE(F06,201) MAT_NAME, SOL_NAME
                  ENDIF
                  HAS_OU4_MAT_BEEN_PROCESSED( I,1) = 'N'
               ENDIF
            ENDIF

         ELSE

            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,946) SUBR_NAME, MAT_NAME, CALLING_SUBR
            WRITE(F06,946) SUBR_NAME, MAT_NAME, CALLING_SUBR
            CALL OUTA_HERE ( 'Y' )

         ENDIF

         CALL DEALLOCATE_SCR_MAT  ( 'CRS1' )
         CALL DEALLOCATE_FULL_MAT ( 'DUM1' )

      ENDDO



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *WARNING    : OUTPUT4 MATRIX "',A,'" IS NOT AVAILABLE FOR PARTITIONING. IT DOES NOT HAVE EITHER COLS OR ROWS THAT', &
                           ' ARE GRID/COMP SET ORIENTED. IT WILL NOT BE PARTITIONED')

  201 FORMAT(' *WARNING    : OUTPUT4 MATRIX "',A,'" IS NOT CURRENTLY ALLOCATED AND IS UNAVAILABLE FOR PARTITIONING IN SOL "',A,'"')

  945 FORMAT(' *ERROR   945: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' FILE UNIT NUMBER ',I8,' IS NOT ONE OF THE NUMBERS ASSIGNED TO OUTPUT4 FILES')

  946 FORMAT(' *ERROR   946: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INVALID INPUT: MAT_NAME = ',A,' CALLING SUBR = ',A)


! **********************************************************************************************************************************

      END SUBROUTINE OUTPUT4_PROC


      SUBROUTINE READ_IN4_FULL_MAT ( ELEM_TYP, ELEM_ID, MAT_NAME_IN, NRI, NCI, UNT, FILNAM, MAT_FULL, IERRT, CALLING_SUBR )

! Reads a matrix that is in full NASTRAN OUTPUT4 format from file FILNAM attached to unit UNT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  NUM_EMG_FATAL_ERRS

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'READ_IN4_FULL_MAT'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Name of subr that called this one
      CHARACTER(LEN=*), INTENT(IN)    :: ELEM_TYP          ! Elem type for which this subr was called
      CHARACTER(LEN=*), INTENT(IN)    :: FILNAM            ! Name of file data is to be read from
      CHARACTER( 1*BYTE)              :: MAT_FND           ! 'Y' if matrix was found on FILNAM
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME_IN       ! Name of matrix to read from UNT
      CHARACTER( 8*BYTE)              :: MAT_NAME          ! Name of matrix in file that is read

      INTEGER(LONG), INTENT(IN)       :: UNT               ! I/O unit number from which to read MAT
      INTEGER(LONG), INTENT(IN)       :: ELEM_ID           ! ID of element for which this subr was called
      INTEGER(LONG), INTENT(IN)       :: NRI               ! Number of rows expected in MAT_FULL
      INTEGER(LONG), INTENT(IN)       :: NCI               ! Number of cols expected in MAT
      INTEGER(LONG), INTENT(OUT)      :: IERRT             ! IERR1+IERR2
      INTEGER(LONG)                   :: FORM              !
      INTEGER(LONG)                   :: ICOL              ! The column number being read in a record from FILNAM
      INTEGER(LONG)                   :: IERR1             ! Local error count
      INTEGER(LONG)                   :: IERR2             ! Local error count
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: IROW              ! Starting row number for terms read in a record from FILNAM
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: NCOLS             ! Number of cols for matrix read from IN4
      INTEGER(LONG)                   :: NROWS             ! Number of rows for matrix read from IN4
      INTEGER(LONG)                   :: NWRDS             ! Number of 4 byte words in a column read from FILNAM
      INTEGER(LONG)                   :: MAT_NUM           ! Number of matrix read from file
      INTEGER(LONG)                   :: NC                ! From matrix trailer. Should be NCOLS+1
      INTEGER(LONG)                   :: PREC              ! Matrix precision (2 indicates double precision)
      INTEGER(LONG)                   :: REC_NUM           !


      REAL(DOUBLE), ALLOCATABLE       :: CCS1_COL(:)       ! One column of MAT
      REAL(DOUBLE), INTENT(OUT)       :: MAT_FULL(NRI,NCI) ! Array of terms in matrix MAT
      REAL(DOUBLE)                    :: RJUNK             ! Values read from file for matrix other than the one we want
!xx   REAL(DOUBLE)                    :: Z0                ! Zero values read from matrix trailer



! **********************************************************************************************************************************
      MAT_NUM = 0
      MAT_FND = 'N'

      DO I=1,NRI
         DO J=1,NCI
            MAT_FULL(I,J) = ZERO
         ENDDO
      ENDDO

do_1: DO                                                   ! Loop over unknown number of matrices in this IN4 file

         IERR1 = 0
         IERR2 = 0
         IERRT = 0
         REC_NUM = 0
         MAT_NUM = MAT_NUM + 1                             ! Read header for a matrix
         READ(UNT,IOSTAT=IOCHK) NCOLS, NROWS, FORM, PREC, MAT_NAME        ;  REC_NUM = REC_NUM + 1
         IF (PREC /= 2) THEN
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR          = FATAL_ERR + 1
            WRITE(ERR,956) PREC, FILNAM
            WRITE(F06,956) PREC, FILNAM
            RETURN
         ENDIF

         IF      (IOCHK == 0) THEN                         ! No problems reading header, so proceed to read matrix values


            IF (MAT_NAME == MAT_NAME_IN) THEN

               IF ((NROWS /= NRI) .OR. (NCOLS /= NCI)) THEN
                  IERR1 = IERR1 + 1
               ENDIF

               IF (IERR1 == 0) THEN

                  MAT_FND = 'Y'

                  ALLOCATE ( CCS1_COL(NROWS) )

                  DO J=1,NCOLS+1                           ! Read cols, 1 at a time. Extra one is not used
                     READ(UNT,IOSTAT=IOCHK) ICOL, IROW, NWRDS, (CCS1_COL(K),K=IROW,IROW+(NWRDS/PREC)-1)     ;  REC_NUM = REC_NUM + 1
                      IF      (IOCHK == 0) THEN            ! No problem reading matrix row, so enter data into MAT_FULL
                        IF (J <= NCOLS) THEN
                           DO K=IROW,IROW+(NWRDS/2)-1
                              MAT_FULL(K,ICOL) = CCS1_COL(K)
                           ENDDO
                        ENDIF
                     ELSE IF (IOCHK > 0) THEN              ! Error reading matrix row so set error condition
                        IERR2 = IERR2 + 1
                        IERRT = IERRT + IERR2
                        CALL IN4_READ_ERR ( 21 )
                     ELSE IF (IOCHK < 0) THEN              ! EOR or EOF before finished reading all rows, so set error condition
                        IERR2 = IERR2 + 1
                        IERRT = IERRT + IERR2
                        CALL IN4_READ_ERR ( 22 )
                     ENDIF
                  ENDDO
                  IF (IERR2 > 0) THEN
                     EXIT do_1
                  ENDIF
                                                           ! Read trailer
!xx               READ(UNT,IOSTAT=IOCHK) NC, IROW, PREC, (Z0,J=1,PREC)                     ;  REC_NUM = REC_NUM + 1
!xx               IF      (IOCHK > 0) THEN                 ! Error reading trailer so set error condition
!xx                  IERR2 = IERR2 + 1
!xx                  IERRT = IERRT + IERR2
!xx                  CALL IN4_READ_ERR ( 31 )
!xx                  EXIT do_1
!xx               ELSE IF (IOCHK < 0) THEN                 ! EOR or EOF  reading trailer so set error condition
!xx                  IERR2 = IERR2 + 1
!xx                  IERRT = IERRT + IERR2
!xx                  CALL IN4_READ_ERR ( 32 )
!xx                  EXIT do_1
!xx               ENDIF


                  DEALLOCATE ( CCS1_COL )

                  EXIT do_1                                ! Matrix was found and read, so exit loop and return to calling subr

               ELSE

                  WRITE(ERR,953) ELEM_TYP, ELEM_ID, MAT_NAME, FILNAM, NRI, NCI, NROWS, NCOLS
                  WRITE(F06,953) ELEM_TYP, ELEM_ID, MAT_NAME, FILNAM, NRI, NCI, NROWS, NCOLS
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR          = FATAL_ERR + 1
                  IERRT              = IERRT + IERR1
                  RETURN

               ENDIF

            ELSE

               DO J=1,NCOLS
                  READ(UNT,IOSTAT=IOCHK) ICOL, IROW, NROWS, (RJUNK,K=1,NROWS)
               ENDDO
               READ(UNT,IOSTAT=IOCHK) NC, IROW, PREC, (RJUNK,J=1,PREC)  ! Read trailer

            ENDIF

         ELSE IF (IOCHK > 0) THEN                          ! Error reading header so set error condition

            IERR2 = IERR2 + 1
            IERRT = IERRT + IERR2
            CALL IN4_READ_ERR ( 11 )
            EXIT do_1

         ELSE IF (IOCHK < 0) THEN                          ! EOR or EOF, so assume no more matrices and exit

            EXIT do_1

         ENDIF

      ENDDO do_1

      IF (ALLOCATED(CCS1_COL)) THEN
         DEALLOCATE ( CCS1_COL )
      ENDIF

      IF (IERR2 > 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9999) IERR2, SUBR_NAME, CALLING_SUBR
         WRITE(F06,9999) IERR2, SUBR_NAME, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF

       IF (MAT_FND == 'N') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,950) MAT_NAME_IN, FILNAM, SUBR_NAME, CALLING_SUBR
         WRITE(F06,950) MAT_NAME_IN, FILNAM, SUBR_NAME, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
       ENDIF



      RETURN

! **********************************************************************************************************************************



  950 FORMAT(' *ERROR   950: INPUT MATRIX "',A,'" WAS NOT FOUND ON FILE:'                                                          &
                    ,/,15X,  A                                                                                                     &
                    ,/,14X,' IN SUBR: ',A,' CALLED BY SUBR: ',A)

  953 FORMAT(' *ERROR   953: FOR ',A,I8,' THE NUMBER OF ROWS AND COLS FOR INPUT MATRIX ',A                                         &
                    ,/,14X,' FROM FILE: ',A                                                                                        &
                    ,/,14X,' SHOULD BE ROWS           = ',I8,', COLS = ',I8,' (BASED ON THE ELEMENT''s B.D CUSERIN ENTRY)'&
                    ,/,14X,' BUT WAS FOUND TO BE ROWS = ',I8,', COLS = ',I8,' (BASED ON THE HEADER RECORD OF THE ABOVE FILE)')

  956 FORMAT(' *ERROR   956: THE PRECISION OF THE MATRIX READ FROM THE FILE BELOW WAS ',I8,'. ONLY PREC -= 2 IS ALLOWED',/,15X,A)

 9999 FORMAT(' PROCESSING STOPPED DUE TO ',I8,' ERROR(S) READING IN4 MATRICES IN SUBR ',A,/,' THIS SUBR WAS CALLED BY ',A)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE IN4_READ_ERR ( II )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06

      IMPLICIT NONE

      CHARACTER(30*BYTE)              :: MESSAGE           ! Error message to write

      INTEGER(LONG), INTENT(IN)       :: II                ! Error message indicator

! **********************************************************************************************************************************
      IF      (II == 11) THEN
         MESSAGE = 'ERROR READING MATRIX HEADER'
      ELSE IF (II == 21) THEN
         MESSAGE = 'ERROR READING MATRIX '
      ELSE IF (II == 22) THEN
         MESSAGE = 'EOR/EOF READING MATRIX '
      ELSE IF (II == 31) THEN
         MESSAGE = 'ERROR READING MATRIX TRAILER'
      ELSE IF (II == 32) THEN
         MESSAGE = 'EOR/EOF READING MATRIX TRAILER'
      ENDIF

      WRITE(ERR,947) MESSAGE, REC_NUM, MAT_NUM, FILNAM, SUBR_NAME
      WRITE(F06,947) MESSAGE, REC_NUM, MAT_NUM, FILNAM, SUBR_NAME

! **********************************************************************************************************************************
  947 FORMAT(' *ERROR   947: ',A,' RECORD NUMBER ',I8,' IN MATRIX NUMBER ',I8,'. ERROR OCCURRED WHILE READING FROM FILE:'          &
               ,/,15X,A,/,14X,' IN SUBR ',A)

! **********************************************************************************************************************************

      END SUBROUTINE IN4_READ_ERR

      END SUBROUTINE READ_IN4_FULL_MAT


      SUBROUTINE WRITE_OU4_FULL_MAT ( MAT_NAME, NROWS, NCOLS, FORM, SYM, MAT, UNT )

! Writes a matrix that is in full format to unformatted file attached to unit UNT in NASTRAN OUTPUT4 format.
! Used for OUTPUT4 matrices

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06, LEN_INPUT_FNAME, OU4, OU4FIL, MOU4
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  PRTOU4

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE


      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_OU4_FULL_MAT'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME          ! Matrix name (only 1st 8 characters will be written)
      CHARACTER(LEN=*), INTENT(IN)    :: SYM               ! 'Y' if input matrix is symmetric
      CHARACTER(LEN_INPUT_FNAME+3)    :: FILNAM            ! The filename for an OUTPUT4 file corresponding to input unit number UNT
      CHARACTER(16*BYTE)              :: MAT_OUT_NAME      ! 16 chars of MAT_NAME (or padded w/ blanks)

      INTEGER(LONG), INTENT(IN)       :: FORM              ! NASTRAN matrix FORM (not really used in MYSTRAN but needed for OUTPUT4)
      INTEGER(LONG), INTENT(IN)       :: NCOLS             ! Number of cols in MAT
      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in MAT
      INTEGER(LONG), INTENT(IN)       :: UNT               ! Unit number where to write matrix
      INTEGER(LONG)                   :: I,J               ! DO loop indices or counters
      INTEGER(LONG), PARAMETER        :: IROW        = 1   ! A term written to UNT for the trailer record (just to be like NASTRAN)
      INTEGER(LONG), PARAMETER        :: PREC        = 2   ! Matrix precision (2 indicates double precision)
      INTEGER(LONG), PARAMETER        :: ROW_BEG     = 1   ! 1st row of matrix output to UNT is row 1


      REAL(DOUBLE) , INTENT(IN)       :: MAT(NROWS,NCOLS)  ! Array of terms in matrix MAT



! **********************************************************************************************************************************
! Get file name for unit UNT

      FILNAM(1:) = ' '
      DO I=1,MOU4
         IF (OU4(I) == UNT) THEN
            FILNAM = OU4FIL(I)
         ENDIF
      ENDDO

      IF (LEN(MAT_NAME) > 16) THEN
         MAT_OUT_NAME(1:) = MAT_NAME(1:16)
      ELSE
         MAT_OUT_NAME(1:) = ' '
         MAT_OUT_NAME(1:) = MAT_NAME(1:)
      ENDIF

! Write matrix header

      WRITE(UNT) NCOLS, NROWS, FORM, PREC, MAT_OUT_NAME(1:4), MAT_OUT_NAME(5:8)

! Write matrix data

      DO J=1,NCOLS
         WRITE(UNT) J, ROW_BEG, 2*NROWS, (MAT(I,J),I=1,NROWS)
      ENDDO

! Write matrix trailer

      WRITE(UNT) NCOLS+1, IROW, PREC, ZERO                    ! Write matrix trailer (PREC words = one DOUBLE)

! Write matrix to f06 file, if requested

      IF (PRTOU4 > 0) THEN

         WRITE(F06,101) MAT_NAME, NCOLS, NROWS, FORM, PREC

         WRITE(F06,104) (J,J=1,NCOLS)
         DO I=1,NROWS
            WRITE(F06,105) I,(MAT(I,J),J=1,NCOLS)
         ENDDO
         WRITE(F06,*)

      ENDIF

      WRITE(F06,1001) MAT_NAME, NROWS, NCOLS, FILNAM



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(14X,A,8X,'NCOLS = ',I8,8X,'NROWS = ',I8,8X,'FORM  = ',I8,8X,'PREC  = ',I8)

  104 FORMAT(32767(21X,I3))

  105 FORMAT(I10,32767(2X,1ES22.14))

 1001 FORMAT(' *INFORMATION: MATRIX ',A,' WITH ',I8,' ROWS AND ',I8,' COLS HAS BEEN WRITTEN IN OUTPUT4 FORMAT TO FILE ',A)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_OU4_FULL_MAT


      SUBROUTINE WRITE_OU4_SPARSE_MAT ( MAT_NAME, NROWS, NCOLS, FORM, SYM, NTERM_MAT, I_MAT, J_MAT, MAT, UNT )

! Writes a matrix that is in sparse CRS format to unformatted file attached to unit UNT in NASTRAN OUTPUT4 format.
! Used for OUTPUT4 matrices

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, LEN_INPUT_FNAME, OU4, OU4FIL, mou4
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  PRTOU4, SPARSTOR
      USE SCRATCH_MATRICES, ONLY      :  I_CRS1, J_CRS1, CRS1, I_CCS1, J_CCS1, CCS1

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_MAT_DIAG_ZEROS
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_SYM_TO_CRS_NONSYM, SPARSE_CRS_SPARSE_CCS

      IMPLICIT NONE


      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_OU4_SPARSE_MAT'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME          ! Matrix name (only 1st 8 characters will be written)
      CHARACTER(LEN=*), INTENT(IN)    :: SYM               ! 'Y' if input matrix is symmetric
      CHARACTER(LEN(MAT_NAME))        :: CCS_MAT_NAME      ! Name for CCS form of MAT
      CHARACTER(LEN(MAT_NAME))        :: CRS_MAT_NAME      ! Name for CRS form of MAT
      CHARACTER(LEN_INPUT_FNAME+3)    :: FILNAM
      CHARACTER(16*BYTE)              :: MAT_OUT_NAME      ! 16 chars of MAT_NAME (or padded w/ blanks)

      INTEGER(LONG), INTENT(IN)       :: FORM              ! NASTRAN matrix FORM (not really used in MYSTRAN but needed for OUTPUT4)
      INTEGER(LONG), INTENT(IN)       :: NCOLS             ! Number of cols in MAT
      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in MAT
      INTEGER(LONG), INTENT(IN)       :: NTERM_MAT         ! Number of nonzero terms in MAT
      INTEGER(LONG), INTENT(IN)       :: I_MAT(NROWS+1)    ! Row indicators for MAT: I_MAT(I+1)-I_MAT(i) = no. terms in row I
      INTEGER(LONG), INTENT(IN)       :: J_MAT(NTERM_MAT)  ! Col numbers in MAT
      INTEGER(LONG), INTENT(IN)       :: UNT               ! Unit number where to write matrix
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG), PARAMETER        :: IROW        = 1   !
!xx   INTEGER(LONG)                   :: NTERM_ROW_CCS1    ! Number of nonzero terms in row I of CRS1
      INTEGER(LONG)                   :: NUM_MAT_DIAG_0    !
      INTEGER(LONG)                   :: NTERM_CCS1        !
      INTEGER(LONG)                   :: NTERM_CRS1        !
      INTEGER(LONG)                   :: NTERM_COL_J       !
      INTEGER(LONG), PARAMETER        :: PREC        = 2   ! Matrix precision (2 indicates double precision)
      INTEGER(LONG), PARAMETER        :: ROW_BEG     = 1   ! 1st row of matrix output to UNT is row 1


      REAL(DOUBLE) , INTENT(IN)       :: MAT(NTERM_MAT)    ! Array of terms in matrix MAT
      REAL(DOUBLE)                    :: CCS1_COL(NROWS)   ! One column of CCS1 in full format



! **********************************************************************************************************************************
! Get file name for unit UNT

      FILNAM(1:) = ' '
      DO I=1,MOU4
         IF (OU4(I) == UNT) THEN
            FILNAM = OU4FIL(I)
         ENDIF
      ENDDO

      IF (LEN(MAT_NAME) > 16) THEN
         MAT_OUT_NAME(1:) = MAT_NAME(1:16)
      ELSE
         MAT_OUT_NAME(1:) = ' '
         MAT_OUT_NAME(1:) = MAT_NAME(1:)
      ENDIF

      CCS_MAT_NAME(1:)   = MAT_NAME(1:)                    ! // ' in CCS format'
      CRS_MAT_NAME(1:)   = MAT_NAME(1:)                    ! // ' in nonsym format'


! Convert to CCS format

      NTERM_CRS1 = NTERM_MAT
                                                           ! If matrix stored as sym (therefore is square), first convert to nonsym
      IF ((SYM == 'Y') .AND. (SPARSTOR == 'SYM')) THEN

         CALL SPARSE_MAT_DIAG_ZEROS ( MAT_NAME, NROWS, NTERM_MAT, I_MAT, J_MAT, NUM_MAT_DIAG_0 )
         NTERM_CRS1 = 2*NTERM_MAT  - (NROWS - NUM_MAT_DIAG_0)
         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NROWS, NTERM_CRS1, SUBR_NAME )
         CALL CRS_SYM_TO_CRS_NONSYM ( MAT_NAME, NROWS, NTERM_MAT, I_MAT, J_MAT, MAT, 'CRS1', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, 'N' )

         NTERM_CCS1 = NTERM_CRS1
         CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NCOLS, NTERM_CCS1, SUBR_NAME )
         CALL SPARSE_CRS_SPARSE_CCS (NROWS,NCOLS,NTERM_CRS1, CRS_MAT_NAME, I_CRS1, J_CRS1, CRS1, CCS_MAT_NAME, J_CCS1, I_CCS1,     &
                                     CCS1, 'N' )

      ELSE                                                 ! Matrix was stored nonsym, so convert to CCS directly

         NTERM_CCS1 = NTERM_MAT
         CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NCOLS, NTERM_CCS1, SUBR_NAME )
         CALL SPARSE_CRS_SPARSE_CCS (NROWS,NCOLS,NTERM_MAT , MAT_NAME    , I_MAT , J_MAT , MAT , CCS_MAT_NAME, J_CCS1, I_CCS1,     &
                                     CCS1, 'N' )

      ENDIF

! Write matrix header

!xx   WRITE(UNT) NCOLS, NROWS, FORM, PREC, MAT_NAME(1:16)
      WRITE(UNT) NCOLS, NROWS, FORM, PREC, MAT_OUT_NAME(1:4), MAT_OUT_NAME(5:8)

! Write matrix data


      IF (PRTOU4 > 0) THEN                                 ! Write matrix to f06 file, if requested
         WRITE(F06,101) MAT_NAME, NCOLS, NROWS, FORM, PREC
         WRITE(F06,104) (J,J=1,NCOLS)
      ENDIF

      K = 0
      DO J=1,NCOLS
         NTERM_COL_J = J_CCS1(J+1) - J_CCS1(J)
         DO I=1,NROWS
            CCS1_COL(I) = ZERO
         ENDDO
         IF (NTERM_COL_J > 0) THEN
            DO I=1,NTERM_COL_J
               K = K + 1
               CCS1_COL(I_CCS1(K)) = CCS1(K)
            ENDDO
         ENDIF                                             ! Always write whole matrix
         WRITE(UNT) J, ROW_BEG, 2*NROWS, (CCS1_COL(I),I=1,NROWS)
         IF (PRTOU4 > 0) THEN
            WRITE(F06,105) J,(CCS1_COL(I),I=1,NROWS)
         ENDIF
      ENDDO

      IF (PRTOU4 > 0) THEN
         WRITE(F06,*)
      ENDIF

      WRITE(UNT) NCOLS+1, IROW, PREC, ZERO                    ! Write matrix trailer (PREC words = one DOUBLE)

      CALL DEALLOCATE_SCR_MAT ( 'CCS1' )
      CALL DEALLOCATE_SCR_MAT ( 'CRS1' )

      WRITE(F06,1001) MAT_NAME, NROWS, NCOLS, FILNAM



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(14X,A,8X,'NCOLS = ',I8,8X,'NROWS = ',I8,8X,'FORM  = ',I8,8X,'PREC  = ',I8)

  104 FORMAT(32767(21X,I3))

  105 FORMAT(I10,32767(2X,1ES22.14))

 1001 FORMAT(' *INFORMATION: MATRIX ',A,' WITH ',I8,' ROWS AND ',I8,' COLS HAS BEEN WRITTEN IN OUTPUT4 FORMAT TO FILE ',A)





! **********************************************************************************************************************************

      END SUBROUTINE WRITE_OU4_SPARSE_MAT

   END MODULE OUTPUT4_FILE_IO
