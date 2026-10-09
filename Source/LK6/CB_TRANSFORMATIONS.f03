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

   MODULE CB_TRANSFORMATIONS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CALC_PHIZL, MERGE_PHIXA

   CONTAINS

      SUBROUTINE CALC_PHIZL

! Merges matrices to get L-set displ transfer matrix:

!            PHIZL    = [ PHIZL1  |  PHIZL2  |  DLR ]
! where
!            PHIZL1   = -KLL(-1)*(MLR + MLL*DLR)
! and
!            PHIZL2   = -EIGEN_VEC*diag(EIGEN_VAL)

! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFL, NDOFR,                                                    &
                                         NTERM_DLR, NTERM_PHIZL, NTERM_PHIZL1, NTERM_PHIZL2 , NTERM_MLL, NTERM_MLR, NTERM_MRL,     &
                                         NUM_CB_DOFS, NVEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE PARAMS, ONLY                :  PRTPHIZL, SPARSTOR
      USE TIMDAT, ONLY                :  TSEC
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VEC, EIGEN_VAL
      USE SPARSE_MATRICES, ONLY       :  SYM_MLL, SYM_MLR, SYM_MRL, SYM_KLL, SYM_DLR, SYM_PHIZL  , SYM_PHIZL1  , SYM_PHIZL2

      USE SPARSE_MATRICES, ONLY       :  I_MLL   , J_MLL   , MLL   , I_MLR   , J_MLR   , MLR   , I_MRL   , J_MRL   , MRL,          &
                                         I_KLL   , J_KLL   , KLL   , I_DLR   , J_DLR   , DLR   ,                                   &
                                         I_PHIZL , J_PHIZL , PHIZL , I_PHIZL1, J_PHIZL1, PHIZL1, I_PHIZL2, J_PHIZL2, PHIZL2

      USE SCRATCH_MATRICES, ONLY      :  I_CRS1, J_CRS1, CRS1, I_CRS2, J_CRS2, CRS2, I_CRS3, J_CRS3, CRS3, I_CCS1, J_CCS1, CCS1


      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM, MATTRNSP_SS
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  FULL_TO_SPARSE_CRS, SPARSE_CRS_SPARSE_CCS
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CNT_NONZ_IN_FULL_MAT
      USE MATRIX_MERGING, ONLY        :  MERGE_MAT_COLS_SSS
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_PHIZL'

      INTEGER(LONG)                   :: AROW_MAX_TERMS    ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: NTERM_CRS1        ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2        ! Number of terms in matrix CRS2
      INTEGER(LONG)                   :: NTERM_CRS3        ! Number of terms in matrix CRS3


      REAL(DOUBLE)                    :: DUM1(NDOFL,NVEC)  ! Intermediate matrix
      REAL(DOUBLE)                    :: SMALL             ! A number used in filtering out small numbers from a full matrix



! **********************************************************************************************************************************
! Part 1: Calculate PHIZL1 = -KLL(-1)*(MLR + MLL*DLR). Use CRS3 to hold (MLR + MLL*DLR)
! -------------------------------------------------------------------------------------

      NTERM_MLR = NTERM_MRL
      CALL ALLOCATE_SPARSE_MAT ( 'MLR', NDOFL, NTERM_MLR, SUBR_NAME )
      IF (NTERM_MLR > 0) THEN                              ! Transpose MRL to get MLR

         CALL MATTRNSP_SS ( NDOFR, NDOFL, NTERM_MRL, 'MRL', I_MRL, J_MRL, MRL, 'MLR', I_MLR, J_MLR, MLR )

      ENDIF

      IF (NTERM_MLL > 0) THEN                              ! There will be term MLL*DLR in CRS3
                                                           ! Put DLR into CCS format in array CCS1
         CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFR, NTERM_DLR, SUBR_NAME )
         CALL SPARSE_CRS_SPARSE_CCS ( NDOFL, NDOFR, NTERM_DLR, 'DLR', I_DLR, J_DLR, DLR, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

                                                           ! Sparse multiply to get CRS1 = MLL*DLR = MLL*CCS1
!                                                            (note: use SYM_DLR for sym indicator for CCS1)
         CALL MATMULT_SSS_NTERM ( 'MLL' , NDOFL, NTERM_MLL , SYM_MLL, I_MLL, J_MLL,                                                &
                                  'DLR' , NDOFR, NTERM_DLR , SYM_DLR, J_CCS1, I_CCS1, AROW_MAX_TERMS,                              &
                                  'CRS1',        NTERM_CRS1 )

         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFL, NTERM_CRS1, SUBR_NAME )

         CALL MATMULT_SSS ( 'MLL' , NDOFL, NTERM_MLL , SYM_MLL, I_MLL , J_MLL , MLL,                                               &
                            'DLR' , NDOFR, NTERM_DLR , SYM_DLR, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                              &
                            'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

         CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

         IF (NTERM_MLR > 0) THEN                           ! Sparse add to get CRS3 = MLR + MLL*DLR = MLR + CRS1

            CALL MATADD_SSS_NTERM ( NDOFL, 'MLR', NTERM_MLR, I_MLR, J_MLR, SYM_MLR, 'MLL*DLR', NTERM_CRS1, I_CRS1, J_CRS1, 'N',    &
                                   'CRS3', NTERM_CRS3 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFL, NTERM_CRS3, SUBR_NAME )
            CALL MATADD_SSS ( NDOFL, 'MLR', NTERM_MLR, I_MLR, J_MLR, MLR, ONE, 'MLL*DLR', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE,   &
                             'CRS3', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

         ELSE                                              ! MLL > 0 but MLR = 0 so set CRS3 to CRS1 (which is MLL*DLR)

            NTERM_CRS3 = NTERM_CRS1
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFL, NTERM_CRS3, SUBR_NAME )
            DO I=1,NDOFL+1
               I_CRS3(I) = I_CRS1(I)
            ENDDO
            DO I=1,NTERM_CRS3
               J_CRS3(I) = J_CRS1(I)
                 CRS3(I) =   CRS1(I)
            ENDDO

         ENDIF

         CALL DEALLOCATE_SCR_MAT ( 'CRS1' )

      ELSE                                                 ! MLL = 0 so check MLR

         IF (NTERM_MLR > 0) THEN                           ! CRS3 = MLR

            NTERM_CRS3 = NTERM_MLR
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFL, NTERM_CRS3, SUBR_NAME )
            DO I=1,NDOFL+1
               I_CRS3(I) = I_MLR(I)
            ENDDO
            DO I=1,NTERM_CRS3
               J_CRS3(I) = J_MLR(I)
                 CRS3(I) =   MLR(I)
            ENDDO

         ELSE                                              ! CRS3 = 0 so PHIZL1 = 0

            NTERM_CRS3 = 0

         ENDIF

      ENDIF

! Now solve for PHIZL1

      IF (NTERM_CRS3 > 0) THEN                             ! CRS3 = MLR + MLL*DLR > 0 so solve KLL*PHIZL1 = CRS3 for PHIZL1
         CALL SOLVE_PHIZL1 ( NTERM_CRS3 )
      ELSE
         NTERM_PHIZL1 = 0
      ENDIF
!xx   CALL ALLOCATE_SPARSE_MAT ( 'PHIZL1', NDOFL, NTERM_PHIZL1, SUBR_NAME )

      CALL DEALLOCATE_SCR_MAT ( 'CRS3' )

! Part 2: Calculate PHIZL2 = -EIGEN_VEC*diag(EIGEN_VAL)
! -----------------------------------------------------

      DO I=1,NDOFL                                         ! 1st calc full matrix with PHIZL2 terms
         DO J=1,NVEC
            DUM1(I,J) = -EIGEN_VEC(I,J)/EIGEN_VAL(J)
         ENDDO
      ENDDO
                                                           ! Now convert to sparse format
      CALL CNT_NONZ_IN_FULL_MAT ( 'PHIZL2', DUM1, NDOFL, NVEC, 'N', NTERM_PHIZL2, SMALL )
      CALL ALLOCATE_SPARSE_MAT ( 'PHIZL2', NDOFL, NTERM_PHIZL2, SUBR_NAME )
      CALL FULL_TO_SPARSE_CRS ( '-EIGEN_VEC*diag(EIGEN_VAL)', NDOFL, NVEC, DUM1, NTERM_PHIZL2, SMALL, SUBR_NAME, 'N',              &
                                 I_PHIZL2, J_PHIZL2, PHIZL2 )

! Part 3: Merge PHIZL1 and PHIZL2 into CRS2
! -------------------------------------

      NTERM_CRS2 = NTERM_PHIZL1 + NTERM_PHIZL2

      CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFL, NTERM_CRS2, SUBR_NAME )

      CALL MERGE_MAT_COLS_SSS   ( 'PHIZL1', NTERM_PHIZL1, I_PHIZL1, J_PHIZL1, PHIZL1, SYM_PHIZL1, NDOFR,                           &
                                  'PHIZL2', NTERM_PHIZL2, I_PHIZL2, J_PHIZL2, PHIZL2, SYM_PHIZL2, NDOFL,                           &
                                  'CRS2',             I_CRS2, J_CRS2, CRS2, 'N'        )

! Part 4: Merge CRS2 (= PHIZL1, PHIZL2)  with DLR to get final PHIZL
! ------------------------------------------------------------

      NTERM_PHIZL = NTERM_CRS2 + NTERM_DLR

      CALL ALLOCATE_SPARSE_MAT ( 'PHIZL', NDOFL, NTERM_PHIZL, SUBR_NAME )

      CALL MERGE_MAT_COLS_SSS  ( 'CRS2', NTERM_CRS2, I_CRS2, J_CRS2, CRS2, 'N'     , NDOFR+NVEC,                                   &
                                 'DLR' , NTERM_DLR , I_DLR , J_DLR , DLR , SYM_DLR , NDOFL,                                        &
                                 'PHIZL' ,             I_PHIZL , J_PHIZL , PHIZL , 'N'        )

      CALL DEALLOCATE_SCR_MAT ( 'CRS2' )

      IF (PRTPHIZL > 0) THEN
         CALL WRITE_SPARSE_CRS ( 'PHIZL','  ','  ', NTERM_PHIZL, NDOFL, I_PHIZL, J_PHIZL, PHIZL )
      ENDIF

! Part 5: Deallocate PHIZL1 and PHIZL2

      CALL DEALLOCATE_SPARSE_MAT ( 'PHIZL1' )
      CALL DEALLOCATE_SPARSE_MAT ( 'PHIZL2' )



      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)


! **********************************************************************************************************************************

      END SUBROUTINE CALC_PHIZL


      SUBROUTINE SOLVE_PHIZL1 ( NTERM_CRS3 )


! Solves KLL*PHIZL1 = CRS3 for PHIZL1   where CRS3 = (MLR + MLL*DLR).

! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_ERR, ERR, F06, SCR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FACTORED_MATRIX, FATAL_ERR, KLL_SDIA, NDOFR, NDOFL, NTERM_DLR,              &
                                         NTERM_PHIZL1, NTERM_KLL, NTERM_KLLs
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE PARAMS, ONLY                :  EPSIL, SOLLIB, SPARSE_FLAVOR, SPARSTOR
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE SCRATCH_MATRICES, ONLY      :  I_CRS3, J_CRS3, CRS3
      USE SPARSE_MATRICES, ONLY       :  I2_PHIZL1, I_PHIZL1, J_PHIZL1, PHIZL1, I2_PHIZL1t, I_PHIZL1t, J_PHIZL1t, PHIZL1t,         &
                                         I_KLL, I2_KLL, J_KLL, KLL, I_KLLs, I2_KLLs, J_KLLs, KLLs

      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE
      USE SPARSE_CRS_ACCESS, ONLY     :  GET_SPARSE_CRS_COL
      USE LAPACK_ADAPTERS, ONLY       :  FBS_LAPACK
      USE SUPERLU_ADAPTERS, ONLY      :  FBS_SUPRLU
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  READ_MATRIX_2
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_I_MAT_FROM_I2_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATTRNSP_SS
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS
      USE L6_WORKSPACE, ONLY          :  ALLOCATE_L6_2

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_PHIZL1  '
      CHARACTER(  1*BYTE)             :: CLOSE_IT          ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER(  8*BYTE)             :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER(  1*BYTE)             :: EQUED             ! 'Y' if KLL stiff matrix was equilibrated in subr EQUILIBRATE
      CHARACTER( 24*BYTE)             :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN
      CHARACTER( 24*BYTE)             :: MODNAM1           ! Name to write to screen to describe module being run
      CHARACTER(  1*BYTE)             :: READ_NTERM        ! 'Y' or 'N' Input to subr READ_MATRIX_1
      CHARACTER(  1*BYTE)             :: NULL_COL          ! 'Y' if a col of CRS3 is null
      CHARACTER(  1*BYTE)             :: OPND              ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SCRFIL            ! File name

      INTEGER(LONG), INTENT(IN)       :: NTERM_CRS3        ! Number of terms in matrix CRS3
      INTEGER(LONG)                   :: DEB_PRT(2)        ! Debug numbers to say whether to write ABAND and/or its decomp to output
!                                                            file in called subr SYM_MAT_DECOMP_LAPACK (ABAND = band form of KOO)

      INTEGER(LONG)                   :: I,J               ! DO loop indices or counters
      INTEGER(LONG)                   :: INFO        = 0   ! Input value for subr SYM_MAT_DECOMP_LAPACK (quit on sing KRRCB)
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      REAL(DOUBLE)                    :: EQUIL_SCALE_FACS(NDOFL)
                                                           ! LAPACK_S values returned from subr SYM_MAT_DECOMP_LAPACK

      REAL(DOUBLE)                    :: INOUT_COL(NDOFL)  ! Temp variable for subr FBS
      REAL(DOUBLE)                    :: K_INORM           ! Inf norm of KLL matrix (det in  subr COND_NUM)
      REAL(DOUBLE)                    :: PHIZL1_COL(NDOFL)   ! A column of PHIZL1   solved for herein
      REAL(DOUBLE)                    :: RCOND             ! Recrip of cond no. of the KLL. Det in  subr COND_NUM

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Solve for PHIZL1

! Open a scratch file that will be used to write PHIZL1   nonzero terms to as we solve for columns of PHIZL1  . After all col's
! of PHIZL1   have been solved for, and we have a count on NTERM_PHIZL1  , we will allocate memory to the PHIZL1   arrays and read
! the scratch file values into those arrays. Then, in the calling subroutine, we will write NTERM_PHIZL1  , followed by
! PHIZL1    row/col/value to a permanent file

      OUNT(1) = ERR
      OUNT(2) = F06
      SCRFIL(1:)  = ' '
      SCRFIL(1:9) = 'SCRATCH-991'
      OPEN (SCR(1),STATUS='SCRATCH',FORM='UNFORMATTED',ACTION='READWRITE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SCRFIL, OUNT )
         CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )
         CALL OUTA_HERE ( 'Y' )                            ! Can't open scratch file, so quit
      ENDIF

! Loop on columns of CRS3

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      NTERM_PHIZL1   = 0
      CALL COUNTER_INIT('   Solve for PHIZL1 col ', NDOFR)
      DO J = 1,NDOFR

! To solve for the j-th col of PHIZL1, use the j-th col of CRS3 (= MLR + MLL*DLR) as a rhs vector. Get the j-th col of CRS3 and put
! the negative of it into array INOUT_COL:

         NULL_COL = 'Y'
         DO I=1,NDOFL
            INOUT_COL(I) = ZERO
            PHIZL1_COL(I)  = ZERO
         ENDDO
         CALL GET_SPARSE_CRS_COL ( 'MLR + MLL*DLR', J,  NTERM_CRS3, NDOFL, NDOFR, I_CRS3, J_CRS3, CRS3, -ONE, INOUT_COL, NULL_COL )

! Calculate PHIZL1_COL via forward/backward substitution.

         IF (NULL_COL == 'N') THEN                         ! FBS will solve for PHIZL1_COL & load it into PHIZL1   array
                                                           ! DPBTRS will return PHIZL1_COL = -KLL(-1)*RHS_col
!                                                            Note 1st arg = 'N' assures that EQUIL_SCAL_FACS will not be used
            IF      (SOLLIB == 'BANDED  ') THEN

               CALL FBS_LAPACK ( 'N', NDOFL, KLL_SDIA, EQUIL_SCALE_FACS, INOUT_COL )

            ELSE IF (SOLLIB == 'SPARSE  ') THEN

               IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

                  INFO = 0
                  CALL FBS_SUPRLU ( SUBR_NAME, 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, J, INOUT_COL, INFO )
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

            DO I=1,NDOFL
               PHIZL1_COL(I) = INOUT_COL(I)
            ENDDO
            DO I=1,NDOFL                                   ! Count NTERM_PHIZL1   and write nonzero PHIZL1   to scratch file
               IF (DABS(PHIZL1_COL(I)) > EPS1) THEN
                  NTERM_PHIZL1   = NTERM_PHIZL1   + 1
                  WRITE(SCR(1)) I, J, PHIZL1_COL(I)
               ENDIF
            ENDDO
         ENDIF
         CALL COUNTER_PROGRESS(J)
      ENDDO

      call deallocate_sparse_mat ( 'KLLs' )

! The PHIZL1 data in SCRATCH-991 is written one col at a time for PHIZL1. Therefore it is rows of PHIZL1t

      REWIND (SCR(1))
      MESSAG = 'SCRATCH: PHIZL1 ROW/COL/VAL'
      READ_NTERM = 'N'
      OPND       = 'Y'
      CLOSE_IT   = 'N'
      CLOSE_STAT = 'KEEP    '

      CALL ALLOCATE_SPARSE_MAT ( 'PHIZL1t', NDOFR, NTERM_PHIZL1, SUBR_NAME )
      CALL ALLOCATE_L6_2 ( 'PHIZL1t', SUBR_NAME )

      CALL ALLOCATE_SPARSE_MAT ( 'PHIZL1', NDOFL, NTERM_PHIZL1, SUBR_NAME )
      CALL ALLOCATE_L6_2 ( 'PHIZL1', SUBR_NAME )
                                                           ! J_PHIZL1t is same as I2_PHIZL1 and I2_PHIZL1t is same as J_PHIZL1
      CALL READ_MATRIX_2 ( SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG, 'PHIZL1t', NDOFL, NTERM_PHIZL1, READ_NTERM,         &
                           J_PHIZL1t, I2_PHIZL1t, PHIZL1t)

! Now get PHIZL1 from PHIZL1t

      CALL GET_I_MAT_FROM_I2_MAT ( 'PHIZL1t', NDOFR, NTERM_PHIZL1, I2_PHIZL1t, I_PHIZL1t )

      CALL MATTRNSP_SS ( NDOFR, NDOFL, NTERM_PHIZL1, 'PHIZL1t', I_PHIZL1t, J_PHIZL1t, PHIZL1t, 'PHIZL1', I_PHIZL1, J_PHIZL1, PHIZL1)

      CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )



      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

 2092 FORMAT(4X,A44,20X,I2,':',I2,':',I2,'.',I3)

 6001 FORMAT(' *ERROR  6001: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX KLL WAS EQUILIBRATED: EQUED = ',A,'. CODE NOT WRITTEN TO ALLOW THIS AS YET')

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(3X,A,I8,' of ',I8,A)








! **********************************************************************************************************************************

      END SUBROUTINE SOLVE_PHIZL1


      SUBROUTINE MERGE_PHIXA ( PART_VEC_A_LR )

! Merges matrices to get CB matrix PHIXA:

!                        |IRR  0RN |     IRR = R-set identity matrix, 0RN  = RxN null matrix (N = num modes)
!                PHIXA = |         |
!                        |DLR  PHIL|     DLR = C.B. boundary modes  , PHIL = EIGEN_VECS() LxN eigenvectors for the L-set x N modes


! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFA, NDOFR, NVEC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VEC
      USE SPARSE_MATRICES, ONLY       :  I_DLR , J_DLR , DLR , I_IRR , J_IRR , IRR , I_PHIXA, J_PHIXA, PHIXA

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MERGE_PHIXA'

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_A_LR(NDOFA)! Partitioning vector (N set into F and S sets)
      INTEGER(LONG)                   :: I,J                 ! DO loop indices or counters
      INTEGER(LONG)                   :: KTERM_IRR           ! Count of number terms in arrays J_IRR and A
      INTEGER(LONG)                   :: KTERM_DLR           ! Count of number terms in arrays J_DLR and B
      INTEGER(LONG)                   :: KTERM_PHIXA         ! Count of number terms in arrays J_C and C
!xx   INTEGER(LONG)                   :: NTERM_PHIXA         ! Number of nonzero terms in output matrix C
      INTEGER(LONG)                   :: NUM_IN_ROW_OF_IRR   ! Num terms in a row of IRR matrix
      INTEGER(LONG)                   :: NUM_IN_ROW_OF_DLR   ! Num terms in a row of DLR matrix
!xx   INTEGER(LONG)                   :: NUM_IN_ROW_OF_PHIXA ! Num terms in a row of IRR
      INTEGER(LONG)                   :: ROW_NUM_DLR         ! Row number in matrix DLR
      INTEGER(LONG)                   :: ROW_NUM_EV          ! Row number in matrix EIGEN_VEC (L-set eigenvectors)
      INTEGER(LONG)                   :: ROW_NUM_IRR         ! Row number in matrix IRR (R-set identity matrix)




! **********************************************************************************************************************************
      ROW_NUM_DLR        = 0
      ROW_NUM_EV         = 0
      ROW_NUM_IRR        = 0
      NUM_IN_ROW_OF_DLR  = 0
      NUM_IN_ROW_OF_IRR  = 0
!xx   NUM_IN_ROW_OF_PHIXA = 0
      KTERM_IRR   = 0
      KTERM_DLR   = 0
      KTERM_PHIXA = 0
      I_PHIXA(1)  = 1
      DO I=1,NDOFA
         I_PHIXA(I+1) = I_PHIXA(I)
         IF      (PART_VEC_A_LR(I) == 2) THEN              ! Get a row of matrix IRR and put it into PHIXA
            ROW_NUM_IRR         = ROW_NUM_IRR + 1
            NUM_IN_ROW_OF_IRR   = I_IRR(ROW_NUM_IRR+1) - I_IRR(ROW_NUM_IRR)
!xx         NUM_IN_ROW_OF_PHIXA = NUM_IN_ROW_OF_IRR
            DO J=1,NUM_IN_ROW_OF_IRR
               KTERM_IRR            = KTERM_IRR + 1
               I_PHIXA(I+1)         = I_PHIXA(I+1) + 1
               KTERM_PHIXA          = KTERM_PHIXA + 1
               J_PHIXA(KTERM_PHIXA) = J_IRR(KTERM_IRR)
                 PHIXA(KTERM_PHIXA) =   IRR(KTERM_IRR)
            ENDDO

         ELSE IF (PART_VEC_A_LR(I) == 1) THEN              ! Get a row of matrix DLR and EIGEN_VEC and put it into PHIXA
            ROW_NUM_DLR             = ROW_NUM_DLR + 1
            ROW_NUM_EV              = ROW_NUM_EV  + 1
            NUM_IN_ROW_OF_DLR       = I_DLR(ROW_NUM_DLR+1) - I_DLR(ROW_NUM_DLR)
!xx         NUM_IN_ROW_OF_PHIXA     = NUM_IN_ROW_OF_DLR
            DO J=1,NUM_IN_ROW_OF_DLR
               KTERM_DLR            = KTERM_DLR + 1
               I_PHIXA(I+1)         = I_PHIXA(I+1) + 1
               KTERM_PHIXA          = KTERM_PHIXA + 1
               J_PHIXA(KTERM_PHIXA) = J_DLR(KTERM_DLR)
                 PHIXA(KTERM_PHIXA) =   DLR(KTERM_DLR)
            ENDDO
            DO J=1,NVEC
               I_PHIXA(I+1)         = I_PHIXA(I+1) +1
               KTERM_PHIXA          = KTERM_PHIXA + 1
               J_PHIXA(KTERM_PHIXA) = NDOFR + J
                 PHIXA(KTERM_PHIXA) = EIGEN_VEC(ROW_NUM_EV,J)
            ENDDO
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MERGE_PHIXA


   END MODULE CB_TRANSFORMATIONS
