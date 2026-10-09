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

   MODULE CB_REDUCED_MATRICES

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CALC_KRRCB, CALC_MRRCB, CALC_MRN, MERGE_KXX, MERGE_MXX

   CONTAINS

      SUBROUTINE CALC_KRRcb

! Calculates the R-set row and col matrix KRRcb in the CB transformation matrix:

!                                  KRRcb = KRR + KRL*DLR

! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FACTORED_MATRIX, FATAL_ERR, KRRcb_SDIA,                                     &
                                         NDOFL, NDOFR, NTERM_DLR, NTERM_KRL, NTERM_KRR, NTERM_KRRcb, NTERM_KRRcbs
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  SPARSTOR, SUPINFO
      USE SPARSE_MATRICES , ONLY      :  SYM_DLR, SYM_KRL, SYM_KRR
      USE SPARSE_MATRICES , ONLY      :  I_KRL, J_KRL, KRL, I_KRR, J_KRR, KRR, I_DLR, J_DLR, DLR,                                  &
                                         I_KRRcb, J_KRRcb, KRRcb, I_KRRcbs, J_KRRcbs, KRRcbs
      USE SCRATCH_MATRICES
      USE LAPACK_DPB_MATRICES, ONLY   :  ABAND

      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_NONSYM_TO_CRS_SYM, SPARSE_CRS_SPARSE_CCS
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_CRS_TERM_COUNT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE LAPACK_MATRIX_LIFECYCLE, ONLY:  DEALLOCATE_LAPACK_MAT
      USE LAPACK_ADAPTERS, ONLY       :  SYM_MAT_DECOMP_LAPACK

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_KRRcb'
      CHARACTER(  1*BYTE)             :: EQUED             ! 'Y' if KRRcb stiff matrix was equilibrated in subr EQUILIBRATE
      CHARACTER(  1*BYTE)             :: SYM_CRS3          ! Storage format for matrix CRS3 (either 'Y' for sym storage or
!                                                            'N' for nonsymmetric storage)

      INTEGER(LONG)                   :: AROW_MAX_TERMS    ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: DEB_PRT(2)        ! Debug numbers to say whether to write ABAND and/or its decomp to output
!                                                            file in called subr SYM_MAT_DECOMP_LAPACK
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: INFO        = -1  ! Input value for subr SYM_MAT_DECOMP_LAPACK (don't quit on sing KRRCB)

      INTEGER(LONG)                   :: NTERM_CRS1        ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS3        ! Number of terms in matrix CRS3


      REAL(DOUBLE)                    :: EQUIL_SCALE_FACS(NDOFR)
                                                           ! LAPACK_S values returned from subr SYM_MAT_DECOMP_LAPACK
      REAL(DOUBLE)                    :: K_INORM           ! Inf norm of KRRcb matrix (det in  subr COND_NUM)
      REAL(DOUBLE)                    :: RCOND             ! Recrip of cond no. of the KLL. Det in  subr COND_NUM



! **********************************************************************************************************************************
! Calc KRRcb = KRR + KRL*DLR
                                                           ! CCS1 will be sparse CCS format version of sparse CRS matrix DLR
      IF (NTERM_KRL > 0) THEN                              ! Part I of KRRcb: calc KRL*DLR add it to KRR

         CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFR, NTERM_DLR, SUBR_NAME )
         CALL SPARSE_CRS_SPARSE_CCS ( NDOFL, NDOFR, NTERM_DLR, 'DLR', I_DLR, J_DLR, DLR, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )


                                                           ! I-1 , sparse multiply to get CRS1 = KRL*DLR. Use CCS1 for DLR CCS
         CALL MATMULT_SSS_NTERM ( 'KRL' , NDOFR, NTERM_KRL , SYM_KRL, I_KRL , J_KRL ,                                              &
                                  'DLR' , NDOFR, NTERM_DLR , SYM_DLR, J_CCS1, I_CCS1, AROW_MAX_TERMS,                              &
                                  'CRS1',        NTERM_CRS1 )

         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFR, NTERM_CRS1, SUBR_NAME )

         CALL MATMULT_SSS ( 'KRL' , NDOFR, NTERM_KRL, SYM_KRL, I_KRL , J_KRL , KRL ,                                               &
                            'DLR' , NDOFR, NTERM_DLR, SYM_DLR, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                               &
                            'CRS1', ONE  , NTERM_CRS1,         I_CRS1, J_CRS1, CRS1 )

         CALL DEALLOCATE_SCR_MAT ( 'CCS1')

         IF      (SPARSTOR == 'SYM   ') THEN               !      If SPARSTOR == 'SYM   ', rewrite CRS1 (= KRL*DLR) as sym in CRS3

            CALL SPARSE_CRS_TERM_COUNT ( NDOFR, NTERM_CRS1, 'CRS1 = KRL*DLR', I_CRS1, J_CRS1, NTERM_CRS3 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFR, NTERM_CRS3, SUBR_NAME )
            CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS1 = KRL*DLR all nonzeros', NDOFR, NTERM_CRS1, I_CRS1, J_CRS1, CRS1,                   &
                                         'CRS3 = KRL*DLR stored sym'  ,        NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )
            SYM_CRS3 = 'Y'

         ELSE IF (SPARSTOR == 'NONSYM') THEN               ! If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

            NTERM_CRS3 = NTERM_CRS1
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFR, NTERM_CRS3, SUBR_NAME )
            DO I=1,NDOFR+1
               I_CRS3(I) = I_CRS1(I)
            ENDDO
            DO I=1,NTERM_CRS3
               J_CRS3(I) = J_CRS1(I)
                 CRS3(I) =   CRS1(I)
            ENDDO
            SYM_CRS3 = 'N'

         ELSE                                              ! Error - incorrect SPARSTOR

            WRITE(ERR,932) SUBR_NAME, SPARSTOR
            WRITE(F06,932) SUBR_NAME, SPARSTOR
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )

         ENDIF
                                                           ! I-8 , sparse add to get KRRcb = KRR + CRS3 = KRR + KRL*DLR
         CALL MATADD_SSS_NTERM ( NDOFR, 'KRR', NTERM_KRR, I_KRR, J_KRR, SYM_KRR, 'KRL*DLR', NTERM_CRS3, I_CRS3, J_CRS3, SYM_CRS3,  &
                                      'KRRcb', NTERM_KRRcb )
         CALL ALLOCATE_SPARSE_MAT ( 'KRRcb', NDOFR, NTERM_KRRcb, SUBR_NAME )
         CALL MATADD_SSS ( NDOFR, 'KRR', NTERM_KRR, I_KRR, J_KRR, KRR, ONE, 'KRL*DLR', NTERM_CRS3, I_CRS3, J_CRS3, CRS3, ONE,      &
                                  'KRRcb', NTERM_KRRcb, I_KRRcb, J_KRRcb, KRRcb )

         CALL DEALLOCATE_SCR_MAT ( 'CRS1' )                ! I-9 , deallocate CRS1


         CALL DEALLOCATE_SCR_MAT ( 'CRS3' )                ! I-12, deallocate CRS3

      ELSE

      NTERM_KRRcb = NTERM_KRR                              ! Allocate KRRcb and equate to KRR since there is no KRL term
      CALL ALLOCATE_SPARSE_MAT ( 'KRRcb', NDOFR, NTERM_KRRcb, SUBR_NAME )

      DO I=1,NDOFR+1
         I_KRRcb(I) = I_KRR(I)
      ENDDO
      DO J=1,NTERM_KRRcb
         J_KRRcb(J) = J_KRR(J)
           KRRcb(J) =   KRR(J)
      ENDDO

      ENDIF

! If DEBUG(104) > 0, check if KRRcb is singular. It should be singular regardless of the number of boundary DOF's.
! (KRR is singular if NDOFR = 6 and is a determinant set of supports. KRRcb should be singular always)
! (CODE ONLY IMPLEMENTED FOR SOLLIB == 'BANDED  ')

      IF (DEBUG(104) > 0) THEN

         CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )
         DEB_PRT(1) = 66
         DEB_PRT(2) = 67

         INFO = -1                                      ! Set INFO, on input, -1 so that SYM_MAT_DECOMP_LAPACK will not abort
         CALL SYM_MAT_DECOMP_LAPACK ( SUBR_NAME, 'KRRcb', 'R ', NDOFR, NTERM_KRRcb, I_KRRcb, J_KRRcb, KRRcb, 'N', 'N', 'N', 'N',&
                                      DEB_PRT, EQUED, KRRcb_SDIA, K_INORM, RCOND, EQUIL_SCALE_FACS, INFO )
         IF (INFO > 0) THEN                             ! KRRcb was singular as it should be
            WRITE(ERR,9971)
            WRITE(ERR,*)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,9971)
               WRITE(F06,*)
            ENDIF
         ELSE                                           ! KRRcb was not singular. Model must be constrained from RB motion
            WRITE(ERR,9972)
            WRITE(ERR,*)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,9972)
               WRITE(F06,*)
            ENDIF
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

 9971 FORMAT(' *INFORMATION: MATRIX KRRcb WAS FOUND TO BE SINGULAR AS IT SHOULD BE. ANY ERRORS REGARDING SINGULARITY OF KRRcb',    &
                           ' SHOULD BE IGNORED')

 9972 FORMAT(' *INFORMATION: MATRIX KRRcb WAS CHECKED FOR SINGULARITY BUT WAS FOUND TO BE NONSINGULAR.',                           &
                           ' USER SHOULD CHECK MODEL TO MAKE SURE IT IS NOT RESTRAINED FROM RIGID BODY MOTION')

! **********************************************************************************************************************************

      END SUBROUTINE CALC_KRRcb


      SUBROUTINE CALC_MRRcb

! Calculates the R-set row and col matrix MRRcb in the CB transformation matrix:

!            MRRcb = MRR + MRL*DLR + (MRL*DLR)' + DLR'*MLL*DLR

! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFL, NDOFR, NTERM_DLR, NTERM_MLL, NTERM_MRL, NTERM_MRR,        &
                                         NTERM_MRRcb, NTERM_MRRcbn
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE PARAMS, ONLY                :  SPARSTOR, WTMASS
      USE RIGID_BODY_DISP_MATS, ONLY  :  TR6_MEFM
      USE MODEL_STUF, ONLY            :  MEFM_RB_MASS
      USE SPARSE_MATRICES , ONLY      :  SYM_DLR, SYM_MLL, SYM_MRL, SYM_MRR, SYM_MRRcb

      USE SPARSE_MATRICES , ONLY      :  I_MLL  , J_MLL  , MLL  , I_MRL  , J_MRL  , MRL  , I_MRR  , J_MRR  , MRR  ,                &
                                         I_DLR  , J_DLR  , DLR  , I_DLRt , J_DLRt , DLRt , I_MRRcb, J_MRRcb, MRRcb,                &
                                         SYM_MRRcb

      USE SCRATCH_MATRICES

      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_NONSYM_TO_CRS_SYM, SPARSE_CRS_SPARSE_CCS, SPARSE_CRS_TO_FULL
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM, MATTRNSP_SS
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_CRS_TERM_COUNT, SPARSE_MAT_DIAG_ZEROS
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF, MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_MRRcb'
      CHARACTER(  1*BYTE)             :: SYM_CRS1            ! Storage format for matrix CRS1 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)
      CHARACTER(  1*BYTE)             :: SYM_CRS3            ! Storage format for matrix CRS3 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)

      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG)                   :: NTERM_CCS1          ! Number of terms in matrix CCS1
      INTEGER(LONG)                   :: NTERM_CRS1          ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2          ! Number of terms in matrix CRS2
      INTEGER(LONG)                   :: NTERM_CRS3          ! Number of terms in matrix CRS3
      INTEGER(LONG)                   :: NUM_MRRcb_DIAG_0    ! Number of zero diagonal terms in MRRcb


      REAL(DOUBLE)                    :: DUMR6(NDOFR,6)      ! Intermediate matrix
                                                             ! Full representation of MRRcb
      REAL(DOUBLE)                    :: MRRcb_FULL(NDOFR,NDOFR)



! **********************************************************************************************************************************
! Calc MRRcb = MRR + MRL*DLR + (MRL*DLR)' + DLR'*MLL*DLR

      NTERM_MRRcb = NTERM_MRR                              ! First, allocate MRRcb and equate to MRR until we get other terms later
      CALL ALLOCATE_SPARSE_MAT ( 'MRRcb', NDOFR, NTERM_MRRcb, SUBR_NAME )

      DO I=1,NDOFR+1
         I_MRRcb(I) = I_MRR(I)
      ENDDO
      DO J=1,NTERM_MRRcb
         J_MRRcb(J) = J_MRR(J)
           MRRcb(J) =   MRR(J)
      ENDDO
                                                           ! CCS1 will be sparse CCS format version of sparse CRS matrix DLR
      CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFR, NTERM_DLR, SUBR_NAME )
      CALL SPARSE_CRS_SPARSE_CCS ( NDOFL, NDOFR, NTERM_DLR, 'DLR', I_DLR, J_DLR, DLR, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

      IF (NTERM_MRL > 0) THEN                              ! Part I of MRRcb: calc MRL*DLR & add it & transpose to MRR
                                                           ! I-1 , sparse multiply to get CRS1 = MRL*DLR. Use CCS1 for DLR CCS
         CALL MATMULT_SSS_NTERM ( 'MRL' , NDOFR, NTERM_MRL , SYM_MRL, I_MRL, J_MRL,                                                &
                                  'DLR' , NDOFR, NTERM_DLR , SYM_DLR, J_CCS1, I_CCS1, AROW_MAX_TERMS,                              &
                                  'CRS1',        NTERM_CRS1 )

         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFR, NTERM_CRS1, SUBR_NAME )

         CALL MATMULT_SSS ( 'MRL' , NDOFR, NTERM_MRL , SYM_MRL, I_MRL , J_MRL , MRL ,                                              &
                            'DLR' , NDOFR, NTERM_DLR , SYM_DLR, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                              &
                            'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

         NTERM_CRS2 = NTERM_CRS1                           ! I-2 , allocate memory to array CRS2 which will hold transpose of CRS1
         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFR, NTERM_CRS2, SUBR_NAME )

                                                           ! I-3 , transpose CRS1 to get CRS2 = (MRL*DLR)t
         CALL MATTRNSP_SS ( NDOFR, NDOFR, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CRS2', I_CRS2, J_CRS2, CRS2 )

                                                           ! I-4 , sparse add to get CRS3 = CRS1 + CRS2 = (MRL*DLR) + (MRL*DLR)t
         CALL MATADD_SSS_NTERM (NDOFR,'MRL*DLR', NTERM_CRS1, I_CRS1, J_CRS1, 'N', '(MRL*DLR)t', NTERM_CRS2, I_CRS2, J_CRS2, 'N',   &
                                      'CRS3', NTERM_CRS3 )
         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFR, NTERM_CRS3, SUBR_NAME )
         CALL MATADD_SSS ( NDOFR, 'MRL*DLR', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE, '(MRL*DLR)t', NTERM_CRS2, I_CRS2, J_CRS2, CRS2,&
                                  ONE, 'CRS3', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

         CALL DEALLOCATE_SCR_MAT ( 'CRS1' )                ! I-5 , deallocate CRS1, CRS2
         CALL DEALLOCATE_SCR_MAT ( 'CRS2' )

                                                           ! I-6 , CRS3 = (MRL*DLR) + (MRL*DLR)t has all nonzero terms in it.
         IF      (SPARSTOR == 'SYM   ') THEN               !       If SPARSTOR == 'SYM   ', rewrite CRS3 as sym in CRS1

            CALL SPARSE_CRS_TERM_COUNT ( NDOFR, NTERM_CRS3, '(MRL*DLR) + (MRL*DLR)t', I_CRS3, J_CRS3, NTERM_CRS1 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFR, NTERM_CRS1, SUBR_NAME )
            CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS3 = (MRL*DLR) + (MRL*DLR)t all nonzeros', NDOFR, NTERM_CRS3, I_CRS3, J_CRS3, CRS3,    &
                                         'CRS1 = (MRL*DLR) + (MRL*DLR)t stored sym'  ,        NTERM_CRS1, I_CRS1, J_CRS1, CRS1 )
            SYM_CRS1 = 'Y'

         ELSE IF (SPARSTOR == 'NONSYM') THEN               !      If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

            NTERM_CRS1 = NTERM_CRS3
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFR, NTERM_CRS1, SUBR_NAME )
            DO I=1,NDOFR+1
               I_CRS1(I) = I_CRS3(I)
            ENDDO
            DO I=1,NTERM_CRS1
               J_CRS1(I) = J_CRS3(I)
                 CRS1(I) =   CRS3(I)
            ENDDO
            SYM_CRS1 = 'N'

         ELSE

            WRITE(ERR,932) SUBR_NAME, SPARSTOR
            WRITE(F06,932) SUBR_NAME, SPARSTOR
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )

         ENDIF

         CALL DEALLOCATE_SCR_MAT ( 'CRS3' )                ! 1-7 , Now CRS1 = (MRL*DLR) + (MRL*DLR)t so deallocate CRS3
                                                           ! I-8 , sparse add: CRS3 = MRR + CRS1 = MRR + (MRL*DLR) + (MRL*DLR)t
         CALL MATADD_SSS_NTERM ( NDOFR, 'MRR', NTERM_MRR, I_MRR, J_MRR, SYM_MRR, 'MRL*DLR + (MRL*DLR)t', NTERM_CRS1,               &
                                 I_CRS1, J_CRS1, SYM_CRS1, 'CRS3', NTERM_CRS3 )
         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFR, NTERM_CRS3, SUBR_NAME )
         CALL MATADD_SSS ( NDOFR, 'MRR', NTERM_MRR, I_MRR, J_MRR, MRR, ONE, 'MRL*DLR + (MRL*DLR)t', NTERM_CRS1,                    &
                           I_CRS1, J_CRS1, CRS1, ONE, 'CRS1', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

         CALL DEALLOCATE_SCR_MAT ( 'CRS1' )                ! I-9 , Now CRS3 = MRR + (MRL*DLR) + (MRL*DLR)t so deallocate CRS1

         NTERM_MRRcb = NTERM_CRS3                          ! I-10, allocate MRRcb to be size of CRS1

         CALL DEALLOCATE_SPARSE_MAT ( 'MRRcb' )            ! Reset MRRcb to CRS3 =
         CALL ALLOCATE_SPARSE_MAT ( 'MRRcb', NDOFR, NTERM_MRRcb, SUBR_NAME )
                                                           ! I-11, set MRRcb = CRS3 = MRR + (MRL*DLR) + (MRL*DLR)t
         DO I=1,NDOFR+1
            I_MRRcb(I) = I_CRS3(I)
         ENDDO
         DO J=1,NTERM_MRRcb
            J_MRRcb(J) = J_CRS3(J)
              MRRcb(J) =   CRS3(J)
         ENDDO

         CALL DEALLOCATE_SCR_MAT ( 'CRS3' )                ! I-12, deallocate CRS3

      ENDIF

      IF (NTERM_MLL > 0) THEN                              ! Part II of MRRcb: calc DLR(t)*MLL*DLR and add to MRRcb

                                                           ! II-1 , sparse multiply to get CRS1 = MLL*DLR using CCS1 for DLR CCS
         CALL MATMULT_SSS_NTERM ( 'MLL' , NDOFL, NTERM_MLL , SYM_MLL, I_MLL , J_MLL,                                               &
                                  'DLR' , NDOFR, NTERM_DLR , SYM_DLR, J_CCS1, I_CCS1, AROW_MAX_TERMS,                              &
                                  'CRS1',        NTERM_CRS1 )

         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFL, NTERM_CRS1, SUBR_NAME )

         CALL MATMULT_SSS ( 'MLL ', NDOFL, NTERM_MLL , SYM_MLL, I_MLL, J_MLL, MLL,                                                 &
                            'DLR' , NDOFR, NTERM_DLR , SYM_DLR, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                              &
                            'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

         CALL DEALLOCATE_SCR_MAT ( 'CCS1' )                ! II-2 , deallocate CCS1

         NTERM_CCS1 = NTERM_CRS1                           ! II-3 , allocate CCS1 to be same as CRS1 but in CCS format
         CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFR, NTERM_CCS1, SUBR_NAME )
         CALL SPARSE_CRS_SPARSE_CCS ( NDOFL, NDOFR, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

         CALL DEALLOCATE_SCR_MAT ( 'CRS1' )                ! II-4 , deallocate CRS1

                                                           ! II-5 , sparse multiply to get CRS1 = DLRt*CCS1 with CCS1 = MLL*DLR
!                                                                   (note: SYM_DLR used for sym indicator for CCS1)
         CALL MATMULT_SSS_NTERM ( 'DLRt', NDOFR, NTERM_DLR , SYM_DLR, I_DLRt, J_DLRt,                                              &
                                  'CCS1', NDOFR, NTERM_CCS1, SYM_DLR, J_CCS1, I_CCS1, AROW_MAX_TERMS,                              &
                                  'CRS1 ',       NTERM_CRS1 )

         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFR, NTERM_CRS1, SUBR_NAME )

         CALL MATMULT_SSS ( 'DLRt', NDOFR, NTERM_DLR , SYM_DLR, I_DLRt, J_DLRt, DLRt,                                              &
                            'CCS1', NDOFR, NTERM_CCS1, SYM_DLR, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                              &
                            'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

         CALL DEALLOCATE_SCR_MAT ( 'CCS1' )                ! II-6 , deallocate CCS1

                                                           ! II-7 , CRS1 = DLRt*MLL*DLR has all nonzero terms in it.
         IF      (SPARSTOR == 'SYM   ') THEN               !        If SPARSTOR == 'SYM   ', rewrite CRS1 as sym in CRS3

            CALL SPARSE_CRS_TERM_COUNT ( NDOFR, NTERM_CRS1, 'DLRt*MLL*DLR all nonzeros', I_CRS1, J_CRS1, NTERM_CRS3 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFR, NTERM_CRS3, SUBR_NAME )
            CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS1 = DLRt*MLL*DLR all nonzeros', NDOFR, NTERM_CRS1, I_CRS1, J_CRS1, CRS1,              &
                                         'CRS3 = DLRt*MLL*DLR stored sym'  ,        NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )
            SYM_CRS3 = 'Y'

         ELSE IF (SPARSTOR == 'NONSYM') THEN               ! If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

            NTERM_CRS3 = NTERM_CRS1
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFR, NTERM_CRS3, SUBR_NAME )
            DO I=1,NDOFR+1
               I_CRS3(I) = I_CRS1(I)
            ENDDO
            DO I=1,NTERM_CRS3
               J_CRS3(I) = J_CRS1(I)
                 CRS3(I) =   CRS1(I)
            ENDDO
            SYM_CRS3 = 'N'

         ELSE

            WRITE(ERR,932) SUBR_NAME, SPARSTOR
            WRITE(F06,932) SUBR_NAME, SPARSTOR
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )

         ENDIF
                                                           ! II-8 , sparse add to get CRS2 = prior MRRcb + CRS3 or
!                                                                   CRS2= MRR + (MRL*DLR) + (MRL*DLR)t + DLRt*MLL*DLR which is MRRcb
         CALL MATADD_SSS_NTERM ( NDOFR, 'MRR', NTERM_MRRcb, I_MRRcb, J_MRRcb, SYM_MRRcb,                                           &
                                        '(MRL*DLR) + (MRL*DLR)t + DLRt*MLL*DLR' , NTERM_CRS3 , I_CRS3 , J_CRS3 , SYM_CRS3,         &
                                        'CRS2' , NTERM_CRS2 )
         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFR, NTERM_CRS2, SUBR_NAME )
         CALL MATADD_SSS ( NDOFR, 'MRR', NTERM_MRRcb, I_MRRcb, J_MRRcb, MRRcb, ONE, 'CRS3', NTERM_CRS3, I_CRS3, J_CRS3, CRS3,      &
                               ONE,'(MRL*DLR) + (MRL*DLR)t + DLRt*MLL*DLR', NTERM_CRS2,  I_CRS2,  J_CRS2,  CRS2 )

         CALL DEALLOCATE_SCR_MAT ( 'CRS1' )                ! II-9 , deallocate CRS1 and CRS3
         CALL DEALLOCATE_SCR_MAT ( 'CRS3' )                ! II-10, deallocate CRS3

         NTERM_MRRcb = NTERM_CRS2                          ! II-11, reallocate MRRcb to be size of CRS2
         CALL DEALLOCATE_SPARSE_MAT ( 'MRRcb' )
         CALL ALLOCATE_SPARSE_MAT ( 'MRRcb', NDOFR, NTERM_MRRcb, SUBR_NAME )

         DO I=1,NDOFR+1                                    ! II-12, set MRRcb = CRS2 until we see if MRL is null, below
            I_MRRcb(I) = I_CRS2(I)
         ENDDO
         DO I=1,NTERM_MRRcb
            J_MRRcb(I) = J_CRS2(I)
              MRRcb(I) =   CRS2(I)
         ENDDO

         CALL DEALLOCATE_SCR_MAT ( 'CRS2' )                ! II-13, deallocate CRS2

      ELSE

         CALL DEALLOCATE_SCR_MAT ( 'CCS1')

      ENDIF

! Calc the 6x6 mass matrix relative to the MEFMGRID and convert it to same units as input mass

      CALL SPARSE_CRS_TO_FULL ( 'MRRcb', NTERM_MRRcb, NDOFR, NDOFR, SYM_MRRcb, I_MRRcb, J_MRRcb, MRRcb, MRRcb_FULL )

      CALL MATMULT_FFF   ( MRRcb_FULL, TR6_MEFM, NDOFR, NDOFR, 6, DUMR6   )
      CALL MATMULT_FFF_T ( TR6_MEFM  , DUMR6   , NDOFR, 6    , 6, MEFM_RB_MASS )
      DO I=1,6
         DO J=1,6
            MEFM_RB_MASS(I,J) = MEFM_RB_MASS(I,J)/WTMASS
         ENDDO
      ENDDO

! The following code sets NTERM_MRRcbn, needed in subr MERGE_MXX. NTERM_MRRcbn cannot wait to be calculated there since several
! arrays have to be dimensioned using it in that subr

      IF      (SPARSTOR == 'SYM   ') THEN                  ! Convert MRRcb (stored symmetric) to MRRcbn (stored nonsymmetric)

         CALL SPARSE_MAT_DIAG_ZEROS ( 'MRRcb', NDOFR, NTERM_MRRcb, I_MRRcb, J_MRRcb, NUM_MRRcb_DIAG_0 )
         NTERM_MRRcbn = 2*NTERM_MRRcb  - (NDOFR - NUM_MRRcb_DIAG_0)

      ELSE IF (SPARSTOR == 'NONSYM') THEN

         NTERM_MRRcbn = NTERM_MRRcb                        ! If SPARSTOR is nonsym, then MRRcb is also stored nonsym

      ELSE

         WRITE(ERR,932) SUBR_NAME, SPARSTOR
         WRITE(F06,932) SUBR_NAME, SPARSTOR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

! **********************************************************************************************************************************

      END SUBROUTINE CALC_MRRcb


      SUBROUTINE CALC_MRN

! Calculates the R-set row by NVEC (number of eigenvectors) col matrix MRN   in the CB transformation matrix:

!                        MRN   = (MRL + DLR'*MLL)*EIGEN_VEC ....... (NOTE: EIGEN_VEC is L-set rows by N=NVEC cols)

! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR,                                                                  &
                                         NDOFL, NDOFR, NTERM_DLR, NTERM_MLL, NTERM_MLLn, NTERM_MPF0, NTERM_MRL, NTERM_MRN,         &
                                         NUM_MLL_DIAG_ZEROS, NVEC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  SPARSTOR
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VEC
      USE SPARSE_MATRICES , ONLY      :  SYM_DLR, SYM_MLL, SYM_MLLn, SYM_MRL, SYM_MRL, SYM_MRN

      USE SPARSE_MATRICES , ONLY      :  I_MLL , J_MLL , MLL , I_MLLn, J_MLLn, MLLn, I_MRL , J_MRL , MRL ,                         &
                                         I_DLR , J_DLR , DLR , I_DLRt, J_DLRt, DLRt, I_MRN, J_MRN, MRN,                            &
                                         I_MPF0, J_MPF0, MPF0

      USE SCRATCH_MATRICES, ONLY      :  I_CCS1, J_CCS1, CCS1, I_CRS1, J_CRS1, CRS1, I_CRS2, J_CRS2, CRS2, I_CRS3, J_CRS3, CRS3

      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_MAT_DIAG_ZEROS
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_SYM_TO_CRS_NONSYM, SPARSE_CRS_SPARSE_CCS
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_FULL_MULTIPLICATION, ONLY:  MATMULT_SFS, MATMULT_SFS_NTERM
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_MRN  '

      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG)                   :: NTERM_CRS1          ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2          ! Number of terms in matrix CRS2
      INTEGER(LONG)                   :: NTERM_CRS3          ! Number of terms in matrix CRS3




! **********************************************************************************************************************************
! Calc MRN = (MRL + DLR'*MLL)*EIGEN_VEC

      NTERM_MRN   = 0                                      ! First, allocate MRN until we get other terms later
      CALL ALLOCATE_SPARSE_MAT ( 'MRN', NDOFR, NTERM_MRN, SUBR_NAME )

      IF (NTERM_MLL > 0) THEN                              ! Part I of MRRcb: calc DLR'*MLL

         IF      (SPARSTOR == 'SYM   ') THEN               ! I-a: Convert MLL to nonsym MLLn, if required. Need for MATMULT_SSS

            CALL SPARSE_MAT_DIAG_ZEROS ( 'MLL', NDOFL, NTERM_MLL, I_MLL, J_MLL, NUM_MLL_DIAG_ZEROS )
            NTERM_MLLn = 2*NTERM_MLL - (NDOFL - NUM_MLL_DIAG_ZEROS)

            CALL ALLOCATE_SPARSE_MAT ( 'MLLn', NDOFL, NTERM_MLLn, SUBR_NAME )

            CALL CRS_SYM_TO_CRS_NONSYM ( 'MLL', NDOFL, NTERM_MLL, I_MLL, J_MLL, MLL, 'MLLn', NTERM_MLLn, I_MLLn, J_MLLn, MLLn, 'Y' )
                                                           ! I-a-1: CCS1 will be sparse CCS format version of sparse CRS matrix MLLn
            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFL, NTERM_MLLn, SUBR_NAME )
            CALL SPARSE_CRS_SPARSE_CCS ( NDOFL, NDOFL, NTERM_MLLn, 'MLLn', I_MLLn, J_MLLn, MLLn, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )
                                                           ! I-a-2: Mult DLRt*CCS1 where CCs1 is MLLn
            CALL MATMULT_SSS_NTERM ( 'DLRt', NDOFR, NTERM_DLR , SYM_DLR , I_DLRt, J_DLRt,                                          &
                                     'MLLn', NDOFL, NTERM_MLLn, SYM_MLLn, J_CCS1, I_CCS1, AROW_MAX_TERMS,                          &
                                     'CRS1' ,       NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFR, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'DLRt', NDOFR, NTERM_DLR , SYM_DLR , I_DLRt, J_DLRt, DLRt,                                          &
                               'MLLn', NDOFL, NTERM_MLLn, SYM_MLLn, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                          &
                               'CRS1', ONE  , NTERM_CRS1,           I_CRS1, J_CRS1, CRS1 )

         ELSE IF (SPARSTOR == 'NONSYM') THEN               ! I-b: Use MLL as is since it is nonsym (needed for MATMULT_SSS)

            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFL, NTERM_MLL, SUBR_NAME )
                                                           ! I-b-1: CCS1 will be sparse CCS format version of sparse CRS matrix MLL
            CALL SPARSE_CRS_SPARSE_CCS ( NDOFL, NDOFL, NTERM_MLL, 'MLL', I_MLL, J_MLL, MLL, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )
                                                           ! I-b-2: Mult DLRt*CCS1 where CCs1 is MLL
            CALL MATMULT_SSS_NTERM ( 'DLRt', NDOFR, NTERM_DLR , SYM_DLR, I_DLRt, J_DLRt,                                           &
                                     'MLL' , NDOFL, NTERM_MLL , SYM_MLL, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1' ,       NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFR, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'DLRt', NDOFR, NTERM_DLR, SYM_DLR, I_DLRt, J_DLRt, DLRt,                                            &
                               'MLL' , NDOFL, NTERM_MLL, SYM_MLL, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                            &
                               'CRS1', ONE  , NTERM_CRS1        , I_CRS1, J_CRS1, CRS1 )

         ELSE

            WRITE(ERR,932) SUBR_NAME, SPARSTOR
            WRITE(F06,932) SUBR_NAME, SPARSTOR
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )

         ENDIF

         IF (DEBUG(103) > 0) THEN                          ! Algorithm for calculating MPF's will not use MRL (or MLR)

            NTERM_CRS2 = NTERM_CRS1                        ! Store CRS1 = DLRt*MLL in CRS2 for use below in calculating MPF0
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFR, NTERM_CRS2, SUBR_NAME )

            DO I=1,NDOFR+1
               I_CRS2 (I) = I_CRS1(I)
            ENDDO
            DO J=1,NTERM_CRS2
               J_CRS2 (J) = J_CRS1(J)
                 CRS2 (J) =   CRS1(J)
            ENDDO


         ENDIF

      ELSE

         NTERM_CRS1 = 0
         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFR, NTERM_CRS1, SUBR_NAME )

      ENDIF

      CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

      IF (NTERM_MRL > 0) THEN                              ! Part II of MRN  : add MRL to CRS1 to get MRL + DLR'*MLL

         CALL MATADD_SSS_NTERM ( NDOFR, 'MRL', NTERM_MRL, I_MRL, J_MRL, SYM_MRL, 'DLRt*MLL', NTERM_CRS1, I_CRS1, J_CRS1, 'N',      &
                                       'CRS3', NTERM_CRS3 )
         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFR, NTERM_CRS3, SUBR_NAME )
         CALL MATADD_SSS ( NDOFR, 'MRL', NTERM_MRL, I_MRL, J_MRL, MRL, ONE, 'DLRt*MLL', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE,     &
                                  'CRS3', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

         CALL DEALLOCATE_SCR_MAT ( 'CRS1' )                ! II-3 , deallocate CRS1 and then reallocate it
         NTERM_CRS1 = NTERM_CRS3
         CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFR, NTERM_CRS1, SUBR_NAME )

         DO I=1,NDOFR+1                                    ! Now CRS1 = MRL + DLR'*MLL
            I_CRS1(I) = I_CRS3(I)
         ENDDO
         DO J=1,NTERM_CRS1
            J_CRS1(J) = J_CRS3(J)
              CRS1(J) =   CRS3(J)
         ENDDO

         CALL DEALLOCATE_SCR_MAT ( 'CRS3' )

      ENDIF

      IF (( NTERM_MLL > 0) .OR. (NTERM_MRL > 0)) THEN      ! Multiply CRS1 = MRL + DLR'*MLL  times EIGEN_VEC to get MRN

         CALL MATMULT_SFS_NTERM ('MRL + DLRt*MLL', NDOFR, NTERM_CRS1, 'N', I_CRS1, J_CRS1, 'EIGEN_VEC', NDOFL, NVEC, EIGEN_VEC     &
                                , AROW_MAX_TERMS, 'MRN  ', NTERM_MRN   )
         CALL DEALLOCATE_SPARSE_MAT ( 'MRN' )
         CALL ALLOCATE_SPARSE_MAT ( 'MRN', NDOFR, NTERM_MRN, SUBR_NAME )
         CALL MATMULT_SFS ('MRL + DLRt*MLL', NDOFR, NTERM_CRS1, 'N', I_CRS1, J_CRS1, CRS1, 'EIGEN_VEC', NDOFL, NVEC, EIGEN_VEC,    &
                            AROW_MAX_TERMS, 'MRN', ONE, NTERM_MRN, I_MRN, J_MRN, MRN  )
         NTERM_MPF0 = NTERM_MRN

         IF (DEBUG(103) == 0) THEN                         ! Include MRL (or MLR) in MPF (modal participation factor)  calculation

            NTERM_MPF0 = NTERM_MRN                         ! MPF0 = MRN = (MRL + DLR'*MLL)*EIGEN_VEC
            CALL ALLOCATE_SPARSE_MAT ( 'MPF0', NDOFR, NTERM_MPF0, SUBR_NAME )
            DO I=1,NDOFR+1
               I_MPF0(I) = I_MRN(I)
            ENDDO
            DO J=1,NTERM_MPF0
               J_MPF0(J) = J_MRN(J)
                 MPF0(J) =   MRN(J)
            ENDDO

         ELSE                                              ! Do not include MRL (or MLR) in MPF calculation

            NTERM_MPF0 = NTERM_CRS2                        ! MPF0 = DLR'*MLL*EIGEN_VEC if MRL not included
            CALL ALLOCATE_SPARSE_MAT ( 'MPF0', NDOFR, NTERM_MPF0, SUBR_NAME )
            CALL MATMULT_SFS ('DLRt*MLL', NDOFR, NTERM_CRS2, 'N', I_CRS2, J_CRS2, CRS2, 'EIGEN_VEC', NDOFL, NVEC, EIGEN_VEC,       &
                               AROW_MAX_TERMS, 'MPF0', ONE, NTERM_MPF0, I_MPF0, J_MPF0, MPF0  )
!xx         DO I=1,NDOFR+1
!xx            I_MPF0(I) = I_CRS2(I)
!xx         ENDDO
!xx         DO J=1,NTERM_CRS2
!xx            J_MPF0(J) = J_CRS2(J)
!xx              MPF0(J) =   CRS2(J)
!xx         ENDDO
            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )

         ENDIF

      ELSE                                                 ! MLL, MRL are null so MRN   is null also

         NTERM_MRN   = 0
         CALL ALLOCATE_SPARSE_MAT ( 'MRN  ', NDOFR, NTERM_MRN  , SUBR_NAME )

      ENDIF

      CALL DEALLOCATE_SCR_MAT ( 'CRS1' )




      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

 97531 format(' J, J_CCS1(J)          = ',2i8)

 97532 format(' I, I_CCS1(I), CCS1(I) = ', 2i8,1es14.6)

! **********************************************************************************************************************************

      END SUBROUTINE CALC_MRN


      SUBROUTINE MERGE_KXX

! Merges matrices to get CB stifness matrix:

!                    | KRRcb    0   |        KRRcb = KRR + KLR'*DLR
!            KXX   = |              |
!                    |  0      Kee  |         Kee  = gen stiffnesses for the eigenvecs of the L-set = EIGEN_VAL(i)*GEN_MASS(i)
!                                                    (this is a diagonal matrix)


! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NDOFR, NTERM_KRRcb, NTERM_KXX  , NVEC
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  PRTKXX
      USE EIGEN_MATRICES_1, ONLY      :  GEN_MASS, EIGEN_VAL
      USE SPARSE_MATRICES , ONLY      :  I_KRRcb, J_KRRcb, KRRcb, I_KXX  , J_KXX  , KXX

      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MERGE_KXX  '

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: K                 ! Counter




! **********************************************************************************************************************************
      NTERM_KXX   = NTERM_KRRcb + NVEC
      CALL ALLOCATE_SPARSE_MAT ( 'KXX', NDOFR+NVEC, NTERM_KXX  , SUBR_NAME )
      DO I=1,NDOFR+1
         I_KXX  (I) = I_KRRcb(I)
      ENDDO
      DO I=NDOFR+2,NDOFR+NVEC+1                            ! This is for the diagonal terms in Kee
         I_KXX  (I) = I_KXX  (I-1) + 1
      ENDDO

      DO J=1,NTERM_KRRcb
         J_KXX  (J) = J_KRRcb(J)
           KXX  (J) = KRRcb(J)
      ENDDO
      K = 0
      DO J=NTERM_KRRcb+1,NTERM_KXX
         K          = K + 1
         J_KXX  (J) = NDOFR + K
           KXX  (J) = EIGEN_VAL(K)*GEN_MASS(K)
      ENDDO


      IF (PRTKXX > 0) THEN
         CALL WRITE_SPARSE_CRS ( 'KXX  ','  ','  ', NTERM_KXX  , NDOFR+NVEC, I_KXX  , J_KXX  , KXX   )
      ENDIF



      RETURN

! **********************************************************************************************************************************


! **********************************************************************************************************************************

      END SUBROUTINE MERGE_KXX


      SUBROUTINE MERGE_MXX

! Merges matrices to get CB stifness matrix:

!                    | MRRcb  MRN   |        MRRcb = MRR + MRL*DLR + (MRL*DLR)' + DLR'*MLL*DLR
!            MXX   = |              |        MRN   = (MRL + DLR'*MLL)*EIGEN_VEC
!                    |  sym    Mee  |        Mee   = generalized MASSES for the eigenvectors of the L-set


! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFR, NVEC, NTERM_MRRcb, NTERM_MRRcbn, NTERM_MRN, NTERM_MXX,    &
                                         NTERM_MXXn
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  PRTMXX, SPARSTOR
      USE EIGEN_MATRICES_1, ONLY      :  GEN_MASS
      USE SPARSE_MATRICES, ONLY       :  SYM_MRRcbn, SYM_MRN  , SYM_MXX  , SYM_MXXn
      USE SPARSE_MATRICES, ONLY       :  I_MRRcb, J_MRRcb, MRRcb, I_MRRcbn, J_MRRcbn, MRRcbn, I_MRN  , J_MRN  , MRN  ,             &
                                         I_MXX  , J_MXX  , MXX  , I_MXXn  , J_MXXn  , MXXn

      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_NONSYM_TO_CRS_SYM, CRS_SYM_TO_CRS_NONSYM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE MATRIX_MERGING, ONLY        :  MERGE_MAT_COLS_SSS, MERGE_MAT_ROWS_SSS
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATTRNSP_SS
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MERGE_MXX  '

      INTEGER(LONG)                   :: I,J                ! DO loop indices
      INTEGER(LONG)                   :: I_GEN_MASS2(NVEC+1)!
      INTEGER(LONG)                   :: I_MNR(NVEC+1)      !
      INTEGER(LONG)                   :: I_MXXa(NDOFR+1)    ! Upper NDOFR rows of MXX
      INTEGER(LONG)                   :: I_MXXb(NVEC+1)     ! Lower NVEC  rows of MXX
      INTEGER(LONG)                   :: J_GEN_MASS2(NVEC)  !
      INTEGER(LONG)                   :: J_MNR(NTERM_MRN)   !
      INTEGER(LONG)                   :: J_MXXa(NTERM_MRRcbn+NTERM_MRN)
      INTEGER(LONG)                   :: J_MXXb(NTERM_MRN+NVEC)
!xx   INTEGER(LONG)                   :: K1,K2              ! Counters
      INTEGER(LONG)                   :: MXXn_MERGE_VEC(NDOFR+NVEC)
      INTEGER(LONG)                   :: NTERM_MNR          !
      INTEGER(LONG)                   :: NTERM_MXXa         !
      INTEGER(LONG)                   :: NTERM_MXXb         !
!xx   INTEGER(LONG)                   :: NUM_MNR_IN_ROW_I   ! Number of terms in row i of MNR


      REAL(DOUBLE)                    :: GEN_MASS2(NVEC)
      REAL(DOUBLE)                    :: MNR(NTERM_MRN)
      REAL(DOUBLE)                    :: MXXa(NTERM_MRRcbn+NTERM_MRN)
      REAL(DOUBLE)                    :: MXXb(NTERM_MRN+NVEC)



! **********************************************************************************************************************************
      NTERM_MNR  = NTERM_MRN
      NTERM_MXXa = NTERM_MRRcbn + NTERM_MRN
      NTERM_MXXb = NTERM_MNR + NVEC

! NTERM_MRRcbn was set in subr CALC_MRRcb

      CALL ALLOCATE_SPARSE_MAT ( 'MRRcbn', NDOFR, NTERM_MRRcbn, SUBR_NAME )

      IF      (SPARSTOR == 'SYM   ') THEN                  ! Convert MRRcb (stored symmetric) to MRRcbn (stored nonsymmetric)

         CALL CRS_SYM_TO_CRS_NONSYM ( 'MRRcb' , NDOFR, NTERM_MRRcb , I_MRRcb , J_MRRcb , MRRcb,                                    &
                                      'MRRcbn',        NTERM_MRRcbn, I_MRRcbn, J_MRRcbn, MRRcbn, 'Y' )
      ELSE IF (SPARSTOR == 'NONSYM') THEN

         DO I=1,NDOFR+1
            I_MRRcbn(I) = I_MRRcb(I)
         ENDDO
         DO J=1,NTERM_MRRcbn
            J_MRRcbn(J) = J_MRRcb(J)
              MRRcbn(J) =   MRRcb(J)
         ENDDO

      ELSE

         WRITE(ERR,932) SUBR_NAME, SPARSTOR
         WRITE(F06,932) SUBR_NAME, SPARSTOR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF

! Allocate enough memory for complete MXXn matrix (nonsym format of MXX  )

      NTERM_MXXn = NTERM_MRRcbn + 2*NTERM_MRN + NVEC
      CALL ALLOCATE_SPARSE_MAT ( 'MXXn', NDOFR+NVEC, NTERM_MXXn, SUBR_NAME )

! Merge MXXa <-- MRRcbn and MRN (both in nonsym format) into temporary MXXa

      CALL MERGE_MAT_COLS_SSS ( 'MRRcbn' , NTERM_MRRcbn, I_MRRcbn, J_MRRcbn, MRRcbn, SYM_MRRcbn, NDOFR,                            &
                                'MRN  '  , NTERM_MRN   , I_MRN   , J_MRN   , MRN   , SYM_MRN   , NDOFR,                            &
                                'MXXa'                 , I_MXXa  , J_MXXa  , MXXa  , 'N' )

! Transpose MRN to MNR

      CALL MATTRNSP_SS ( NDOFR, NVEC, NTERM_MRN, 'MRN', I_MRN, J_MRN, MRN, 'MNR', I_MNR, J_MNR, MNR )

! Put gen mass data into sparse array GEN_MASS2:

      DO I=1,NVEC
         I_GEN_MASS2(I) = I
         J_GEN_MASS2(I) = I
           GEN_MASS2(I) = GEN_MASS(I)
      ENDDO
      I_GEN_MASS2(NVEC+1) = NVEC+1

! Merge MXXb <-- MNR and GEN_MASS

!xx    I_MXXb(1) = I_MNR(1)
!xx    DO I=1,NVEC
!xx       I_MXXb(I+1) = I_MNR(I+1) + I
!xx    ENDDO

!xx    K1 = 0
!xx    K2 = 0
!xx    DO I=1,NVEC
!xx       NUM_MNR_IN_ROW_I = I_MNR(I+1) - I_MNR(I)
!xx       DO J=1,NUM_MNR_IN_ROW_I
!xx          K1 = K1 + 1
!xx          K2 = K2 + 1
!xx          J_MXXb(K2) = J_MNR(K1)
!xx            MXXb(K2) =   MNR(K1)
!xx       ENDDO
!xx       K2 = K2 + 1
!xx       J_MXXb(K2) = J_MNR(K1) + I
!xx         MXXb(K2) = GEN_MASS(I)
!xx    ENDDO

      CALL MERGE_MAT_COLS_SSS ( 'MNR'     , NTERM_MNR, I_MNR      , J_MNR      , MNR      , 'N', NDOFR,                           &
                                'GEN_MASS', NVEC     , I_GEN_MASS2, J_GEN_MASS2, GEN_MASS2, 'N', NVEC,                            &
                                'MXXb'               , I_MXXb     , J_MXXb     , MXXb     , 'N' )

! Merge MXXa rows with rows of MXXb into MXXn

      DO I=1,NDOFR
         MXXn_MERGE_VEC(I) = 1
      ENDDO
      DO I=NDOFR+1,NDOFR+NVEC
         MXXn_MERGE_VEC(I) = 2
      ENDDO

      CALL MERGE_MAT_ROWS_SSS ( 'MXXa', NDOFR, NTERM_MXXa, I_MXXa, J_MXXa, MXXa, 1,                                                &
                                'MXXb', NVEC , NTERM_MXXb, I_MXXb, J_MXXb, MXXb, 2, MXXn_MERGE_VEC,                                &
                                'MXXn'                   , I_MXXn, J_MXXn, MXXn )

! Convert MXXn to symmetric format MXX

      IF (SPARSTOR == 'SYM   ') THEN

         NTERM_MXX = NTERM_MRRcb + NTERM_MRN + NVEC
         CALL ALLOCATE_SPARSE_MAT   ( 'MXX' , NDOFR+NVEC, NTERM_MXX , SUBR_NAME )
         CALL CRS_NONSYM_TO_CRS_SYM ( 'MXXn', NDOFR+NVEC, NTERM_MXXn, I_MXXn, J_MXXn, MXXn,                                        &
                                      'MXX' ,             NTERM_MXX , I_MXX , J_MXX , MXX   )

      ELSE

         NTERM_MXX   = NTERM_MXXn                          ! If SPARSTOR is nonsym, then MXX is also stored nonsym
         CALL ALLOCATE_SPARSE_MAT ( 'MXX', NDOFR+NVEC, NTERM_MXX  , SUBR_NAME )
         DO I=1,NDOFR+NVEC+1
            I_MXX(I) = I_MXXn(I)
         ENDDO
         DO J=1,NTERM_MXX
            J_MXX(J) = J_MXXn(J)
              MXX(J) =   MXXn(J)
         ENDDO

      ENDIF

      CALL DEALLOCATE_SPARSE_MAT ( 'MXXn' )

      IF (PRTMXX > 0) THEN
         CALL WRITE_SPARSE_CRS ( 'MXX  ','  ','  ', NTERM_MXX  , NDOFR+NVEC, I_MXX  , J_MXX  , MXX   )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

! **********************************************************************************************************************************

      END SUBROUTINE MERGE_MXX

   END MODULE CB_REDUCED_MATRICES
