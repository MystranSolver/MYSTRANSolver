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

   MODULE CB_LOAD_TRANSFORMATIONS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: INTERFACE_FORCE_LTM, NET_CG_LOADS_LTM, MERGE_LTM, SOLVE_DLR

   CONTAINS

      SUBROUTINE INTERFACE_FORCE_LTM

! Merges matrices to get the interface force Loads Transformation Matrix (LTM):

!            IF_LTM   = | MRRcb  MRN    KRRcb |

! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFR, NTERM_KRRcb, NTERM_KRRcbn, NTERM_MRRcbn, NTERM_MRN  ,     &
                                         NTERM_IF_LTM  , NVEC
      USE PARAMS, ONLY                :  PRTIFLTM, SPARSTOR
      USE TIMDAT, ONLY                :  TSEC

      USE SPARSE_MATRICES, ONLY       :  SYM_KRRcb, SYM_KRRcbn, SYM_MRN  , SYM_MRRcbn, SYM_IF_LTM

      USE SPARSE_MATRICES, ONLY       :  I_MRRcbn   , J_MRRcbn   , MRRcbn   , I_MRN      , J_MRN      , MRN      ,                 &
                                         I_KRRcb    , J_KRRcb    , KRRcb    , I_KRRcbn   , J_KRRcbn   , KRRcbn   ,                 &
                                         I_IF_LTM   , J_IF_LTM   , IF_LTM

      USE SCRATCH_MATRICES, ONLY      :  I_CRS1, J_CRS1, CRS1


      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_MAT_DIAG_ZEROS
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_SYM_TO_CRS_NONSYM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE MATRIX_MERGING, ONLY        :  MERGE_MAT_COLS_SSS
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'INTERFACE_FORCE_LTM'

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: NCOL_CRS1         ! Number of cols in scratch matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS1        ! Number of nonzero terms in scratch matrix CRS1
      INTEGER(LONG)                   :: NUM_KRRcb_DIAG_0  ! Number of zeros on the diagonal of matrix KRRcb





! **********************************************************************************************************************************
! Set KRRcbn based on SPARSTOR (we need KRRcb in nonsym format since IF_LTM   will be nonsym

      IF      (SPARSTOR == 'SYM   ') THEN                  ! Convert KRRcb (stored symmetric) to KRRcbn (stored nonsymmetric)

         CALL SPARSE_MAT_DIAG_ZEROS ( 'KRRcb', NDOFR, NTERM_KRRcb, I_KRRcb, J_KRRcb, NUM_KRRcb_DIAG_0 )
         NTERM_KRRcbn = 2*NTERM_KRRcb  - (NDOFR - NUM_KRRcb_DIAG_0)
         CALL ALLOCATE_SPARSE_MAT ( 'KRRcbn', NDOFR, NTERM_KRRcbn, SUBR_NAME )
         CALL CRS_SYM_TO_CRS_NONSYM ( 'KRRcb' , NDOFR, NTERM_KRRcb , I_KRRcb , J_KRRcb , KRRcb,                                    &
                                      'KRRcbn',        NTERM_KRRcbn, I_KRRcbn, J_KRRcbn, KRRcbn, 'Y' )
      ELSE IF (SPARSTOR == 'NONSYM') THEN

         NTERM_KRRcbn = NTERM_KRRcb                        ! If SPARSTOR is nonsym, then KRRcb is also stored nonsym
         CALL ALLOCATE_SPARSE_MAT ( 'KRRcbn', NDOFR, NTERM_KRRcbn, SUBR_NAME )
         DO I=1,NDOFR+1
            I_KRRcbn(I) = I_KRRcb(I)
         ENDDO
         DO J=1,NTERM_KRRcbn
            J_KRRcbn(J) = J_KRRcb(J)
              KRRcbn(J) =   KRRcb(J)
         ENDDO

      ELSE

         WRITE(ERR,932) SUBR_NAME, SPARSTOR
         WRITE(F06,932) SUBR_NAME, SPARSTOR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF


! Allocate enough memory for merge of cols of MRRcbn with MRN

      NTERM_CRS1 = NTERM_MRRcbn+NTERM_MRN
      CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFR, NTERM_CRS1, SUBR_NAME )

! Merge MRRcbn and MRN   (both in nonsym format) into nonsym format temporary scratch matrix CRS1.

      CALL MERGE_MAT_COLS_SSS ( 'MRRcbn' , NTERM_MRRcbn   , I_MRRcbn, J_MRRcbn, MRRcbn, SYM_MRRcbn, NDOFR,                         &
                                'MRN  '  , NTERM_MRN      , I_MRN   , J_MRN   , MRN   , SYM_MRN   , NDOFR,                         &
                                'MRRcbn merged with MRN  ', I_CRS1  , J_CRS1  , CRS1  , 'N'        )

! Merge CRS1 with KRRcb to get IF_LTM

      NCOL_CRS1      = NDOFR + NVEC
      NTERM_IF_LTM   = NTERM_CRS1 + NTERM_KRRcbn
      CALL ALLOCATE_SPARSE_MAT ( 'IF_LTM  ', NDOFR, NTERM_IF_LTM  , SUBR_NAME )
      CALL MERGE_MAT_COLS_SSS ( 'CRS1'    , NTERM_CRS1    , I_CRS1    , J_CRS1    , CRS1    , 'N'         , NCOL_CRS1,             &
                                'KRRcbn'  , NTERM_KRRcbn  , I_KRRcbn  , J_KRRcbn  , KRRcbn  , SYM_KRRcbn  , NDOFR,                 &
                                'IF_LTM',                   I_IF_LTM  , J_IF_LTM  , IF_LTM  , SYM_IF_LTM   )
      CALL DEALLOCATE_SCR_MAT ( 'CRS1' )
      CALL DEALLOCATE_SPARSE_MAT ( 'KRRcbn' )

      IF (PRTIFLTM > 0) THEN
         CALL WRITE_SPARSE_CRS ( 'IF_LTM  ','  ','  ', NTERM_IF_LTM  , NDOFR, I_IF_LTM  , J_IF_LTM  , IF_LTM   )
      ENDIF




      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)



! **********************************************************************************************************************************

      END SUBROUTINE INTERFACE_FORCE_LTM


      SUBROUTINE NET_CG_LOADS_LTM

! Merges matrices to get the net CG Loads Transformation Matrix (LTM):

!            CG_LTM   = MCG(-1)*TR6_CG'*| MRRcb  MRN  0RR |     6x(2R+N)

! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFR, NTERM_MRRcbn, NTERM_MRN, NTERM_CG_LTM, NUM_CB_DOFS
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE PARAMS, ONLY                :  PRTCGLTM, WTMASS
      USE RIGID_BODY_DISP_MATS, ONLY  :  TR6_CG, TR6_0
      USE MODEL_STUF, ONLY            :  MCG
      USE OUTPUT4_MATRICES, ONLY      :  RBM0
      USE SPARSE_MATRICES, ONLY       :  SYM_MRN   , SYM_MRRcbn, SYM_CG_LTM
      USE SPARSE_MATRICES, ONLY       :  I_MRRcbn  , J_MRRcbn  , MRRcbn   ,  I_MRN      , J_MRN      , MRN      ,                  &
                                         I_CG_LTM  , J_CG_LTM  , CG_LTM

      USE SCRATCH_MATRICES, ONLY      :  I_CRS1, J_CRS1, CRS1, I_CRS2, J_CRS2, CRS2, I_CCS1, J_CCS1, CCS1


      USE SPARSE_FULL_MULTIPLICATION, ONLY:  MATMULT_SFF
      USE FULL_MATRIX_ALGEBRA, ONLY   :  INVERT_FF_MAT, MATMULT_FFF, MATMULT_FFF_T
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CNT_NONZ_IN_FULL_MAT
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  FULL_TO_SPARSE_CRS, SPARSE_CRS_SPARSE_CCS
      USE MATRIX_MERGING, ONLY        :  MERGE_MAT_COLS_SSS
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATMULT_SSS, MATMULT_SSS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'NET_CG_LOADS_LTM'

      INTEGER(LONG)                   :: AROW_MAX_TERMS    ! Max number of terms in any row of matrix A sent to subr MATMULT_SSS
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: INFO              ! Info designator from subrs DPOTRF and DPOTRI when INVERT_FF_MAT called
      INTEGER(LONG)                   :: NTERM_CCS1        ! Number of nonzero terms in scratch matrix CCS1
      INTEGER(LONG)                   :: NTERM_CRS1        ! Number of nonzero terms in scratch matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2        ! Number of nonzero terms in scratch matrix CRS2


      REAL(DOUBLE)                    :: DUM1(NDOFR,6)     ! MRRcbn*TR6_CG
      REAL(DOUBLE)                    :: DUM2(6,NDOFR)     !
      REAL(DOUBLE)                    :: WCG(6,6)          ! MCG/WTMASS
      REAL(DOUBLE)                    :: WBASIC(6,6)       ! RBM0/WTMASS
      REAL(DOUBLE)                    :: MCGI(6,6)         ! MCG inverse
      REAL(DOUBLE)                    :: SMALL             ! A number used in filtering out small numbers from a full matrix
      REAL(DOUBLE)                    :: TR6_CGt(6,NDOFR)  ! TR6_CG'



! **********************************************************************************************************************************

! Calc MCG = TR6_CG'*MRRcbn*TR6_CG
                                                           ! 1st, multiply MRRcbn*TR6_CG to get DUM1
      CALL MATMULT_SFF('MRRcbn', NDOFR, NDOFR, NTERM_MRRcbn, SYM_MRRcbn, I_MRRcbn, J_MRRcbn, MRRcbn, 'TR6_CG', NDOFR, 6, TR6_CG,   &
                        'N','DUM1', ONE, DUM1)
      CALL MATMULT_FFF_T ( TR6_CG, DUM1, NDOFR, 6, 6, MCG )

! Calc and print WCG

      DO I=1,6
         DO J=1,6
            WCG(I,J)  = MCG(I,J)/WTMASS
         ENDDO
      ENDDO

      WRITE(F06,101)
      WRITE(F06,103)
      DO I=1,3
         WRITE(F06,109) (WCG(I,J),J=1,6)
      ENDDO
      WRITE(F06,103)
      DO I=4,6
         WRITE(F06,109) (WCG(I,J),J=1,6)
      ENDDO
      WRITE(F06,103)
      WRITE(F06,*)



! Calc RBM0 = TR6_0'*MRRcbn*TR6_0

      CALL MATMULT_SFF('MRRcbn', NDOFR, NDOFR, NTERM_MRRcbn, SYM_MRRcbn, I_MRRcbn, J_MRRcbn, MRRcbn, 'TR6_0', NDOFR, 6, TR6_0,    &
                        'N','DUM1', ONE, DUM1)
      CALL MATMULT_FFF_T ( TR6_0, DUM1, NDOFR, 6, 6, RBM0 )

! Calc and print WBASIC

      DO I=1,6
         DO J=1,6
            WBASIC(I,J)  = RBM0(I,J)/WTMASS
         ENDDO
      ENDDO

      WRITE(F06,102)
      WRITE(F06,103)
      DO I=1,3
         WRITE(F06,109) (WBASIC(I,J),J=1,6)
      ENDDO
      WRITE(F06,103)
      DO I=4,6
         WRITE(F06,109) (WBASIC(I,J),J=1,6)
      ENDDO
      WRITE(F06,103)
      WRITE(F06,*)

! Invert MCG (to MCGI). Note INFO is returned from INVERT_FF_MAT but any necessary action was taken there

      DO I=1,6                                             ! 1st set MCGI = MCG since INVERT will write over the input matrix
         DO J=1,6
            MCGI(I,J) = MCG(I,J)
         ENDDO
      ENDDO

      CALL INVERT_FF_MAT ( SUBR_NAME, 'MCG', MCGI, 6, INFO )

! Calc DUM2 =  MCG(-1)*TR6_CG'. First, transpose TR6_CG to TR6_CGt. Then rewrite DUM2 as a sparse matrix CRS1 so we can use
! MATMULT_SSS to multiply it times the sparse matrix CRS2 below (to get final CG loads LTM)
! If user wants cg LTM calculated such that translational terms are in G's, scale the upper 3 rows of the LTM (but do it on DUM2
! since it is easier here)

      DO I=1,6
         DO J=1,NDOFR
            TR6_CGt(I,J) = TR6_CG(J,I)
         ENDDO
      ENDDO

      CALL MATMULT_FFF ( MCGI, TR6_CGt, 6, 6, NDOFR, DUM2 )
!     IF (SC_CGLTM == 'Y') THEN
         DO I=1,3
            DO J=1,NDOFR
               DUM2(I,J) = WTMASS*DUM2(I,J)
            ENDDO
         ENDDO
!     ENDIF

      CALL CNT_NONZ_IN_FULL_MAT ( 'MCGI*TR6_CGt', DUM2, 6, NDOFR, 'N', NTERM_CRS1, SMALL )
      CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', 6, NTERM_CRS1, SUBR_NAME )
      CALL FULL_TO_SPARSE_CRS ( 'MCGI*TR6_CGt', 6, NDOFR, DUM2, NTERM_CRS1, SMALL, SUBR_NAME, 'N', I_CRS1, J_CRS1, CRS1 )


! Merge MRRcbn and MRN (both in nonsym format) into nonsym format temporary scratch matrix CRS2.
! First allocate enough memory for merge of cols of MRRcbn with MRN  .

      NTERM_CRS2 = NTERM_MRRcbn + NTERM_MRN
      CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFR, NTERM_CRS2, SUBR_NAME )
      CALL MERGE_MAT_COLS_SSS ( 'MRRcbn', NTERM_MRRcbn, I_MRRcbn, J_MRRcbn, MRRcbn, SYM_MRRcbn, NDOFR,                             &
                                'MRN  ' , NTERM_MRN   , I_MRN   , J_MRN   , MRN   , SYM_MRN   , NDOFR,                             &
                                'CRS2'  ,               I_CRS2  , J_CRS2  , CRS2  , 'N'        )

! Calc CG loads LTM: Mult DUM2 = MCG(-1)*TR6_CG' times CRS2 = | MRRcbn  MRN 0 |. First, convert CRS2 to sparse col storage in CCS1

      NTERM_CCS1 = NTERM_CRS2                              ! CCS1 will be CCS strage of | MRRcbn  MRN 0 |
      CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NUM_CB_DOFS, NTERM_CCS1, SUBR_NAME )
      CALL SPARSE_CRS_SPARSE_CCS ( NDOFR, NUM_CB_DOFS, NTERM_CRS2, 'CRS2', I_CRS2, J_CRS2, CRS2, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y')
      CALL DEALLOCATE_SCR_MAT ( 'CRS2' )

      CALL MATMULT_SSS_NTERM ( 'MCGI*TR6_CGt'     , 6          , NTERM_CRS1  , 'N', I_CRS1, J_CRS1,                                &
                               '| MRRcbn  MRN 0 |', NUM_CB_DOFS, NTERM_CCS1  , 'N', J_CCS1, I_CCS1, AROW_MAX_TERMS,                &
                               'CG_LTM'           ,              NTERM_CG_LTM   )

      CALL ALLOCATE_SPARSE_MAT ( 'CG_LTM', NDOFR, NTERM_CG_LTM  , SUBR_NAME )

      CALL MATMULT_SSS ( 'MCGI*TR6_CGt'     , 6          , NTERM_CRS1  , 'N', I_CRS1  , J_CRS1  , CRS1  ,                          &
                         '| MRRcbn  MRN 0 |', NUM_CB_DOFS, NTERM_CCS1  , 'N', J_CCS1  , I_CCS1  , CCS1  ,  AROW_MAX_TERMS,         &
                         'CG_LTM'          , ONE         , NTERM_CG_LTM,      I_CG_LTM, J_CG_LTM, CG_LTM   )

      CALL DEALLOCATE_SCR_MAT ( 'CRS1' )
      CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

      IF (PRTCGLTM > 0) THEN
         CALL WRITE_SPARSE_CRS ( 'CG_LTM','  ','  ', NTERM_CG_LTM, NDOFR, I_CG_LTM, J_CG_LTM, CG_LTM   )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(40X,'RIGID BODY WEIGHT MATRIX RELATIVE TO THE CG OF THE MODEL',/,                                                     &
             40X,'                       MATRIX WCG',/)

  102 FORMAT(35X,'RIGID BODY WEIGHT MATRIX RELATIVE TO THE MODEL BASIC SYSTEM ORIGIN',/,                                           &
             40X,'                      MATRIX WBASIC',/)

  103 FORMAT(19X,'-----------------------------------------------------------------------------------------------------')

  109 FORMAT(19X,'|',2(1ES13.6,4X),1ES13.6,'  | ',2(1ES13.6,4X),1ES13.6,' |')

  142 FORMAT(21X,6(1ES15.6))


! **********************************************************************************************************************************

      END SUBROUTINE NET_CG_LOADS_LTM


      SUBROUTINE MERGE_LTM

! Merges CG_LTM and IF_LTM into LTM matrix:

!                    | CG_LTM |        6 x NUM_CB_DOFS      1st 6 rows of LTM
!            LTM   = |        |
!                    | IF_LTM |        NDOFR x NUM_CB_DOFS  last NDOFR rows of LTM


! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFR, NTERM_CG_LTM, NTERM_IF_LTM, NTERM_LTM, NUM_CB_DOFS
      USE TIMDAT, ONLY                :  TSEC
      USE SPARSE_MATRICES, ONLY       :  I_CG_LTM, J_CG_LTM, CG_LTM, I_IF_LTM, J_IF_LTM, IF_LTM, I_LTM, J_LTM, LTM

      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_MERGING, ONLY        :  MERGE_MAT_ROWS_SSS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MERGE_LTM  '

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: LTM_MERGE_VEC(6+NDOFR)




! **********************************************************************************************************************************
! Merge CG_LTM rows with rows ofIF_LTM into LTM

      NTERM_LTM = NTERM_CG_LTM + NTERM_IF_LTM

      CALL ALLOCATE_SPARSE_MAT ( 'LTM', 6+NDOFR, NTERM_LTM, SUBR_NAME )

      DO I=1,6
         LTM_MERGE_VEC(I) = 1
      ENDDO
      DO I=7,6+NDOFR
         LTM_MERGE_VEC(I) = 2
      ENDDO

      CALL MERGE_MAT_ROWS_SSS ( 'CG_LTM', 6    , NTERM_CG_LTM, I_CG_LTM, J_CG_LTM, CG_LTM, 1,                                      &
                                'IF_LTM', NDOFR, NTERM_IF_LTM, I_IF_LTM, J_IF_LTM, IF_LTM, 2, LTM_MERGE_VEC,                       &
                                'LTM'                        , I_LTM   , J_LTM   , LTM )


      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MERGE_LTM


      SUBROUTINE SOLVE_DLR

! Solves KLL*DLR = -KLR for matrix DLR. However, we will use rows of KRL instead of cols of KLR in the solution.

! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_ERR, ERR, F06, SCR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FACTORED_MATRIX, FATAL_ERR, KLL_SDIA, NDOFR, NDOFL, NTERM_DLR, NTERM_KLL,   &
                                         NTERM_KRL
      USE PARAMS, ONLY                :  EPSIL, PRTDLR, SOLLIB, SPARSE_FLAVOR, SPARSTOR
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE SPARSE_MATRICES, ONLY       :  I2_DLR, I_DLR, J_DLR, DLR, I_DLRt, I2_DLRt, J_DLRt, DLRt, I_KRL, J_KRL, KRL,              &
                                         I_KLL, I2_KLL, J_KLL, KLL


      USE LAPACK_ADAPTERS, ONLY       :  FBS_LAPACK, SYM_MAT_DECOMP_LAPACK
      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE SUPERLU_ADAPTERS, ONLY      :  FBS_SUPRLU, SYM_MAT_DECOMP_SUPRLU
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE
      USE SPARSE_CRS_ACCESS, ONLY     :  GET_SPARSE_CRS_ROW
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  READ_MATRIX_2, WRITE_SPARSE_CRS
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_I_MAT_FROM_I2_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATTRNSP_SS
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS
      USE L6_WORKSPACE, ONLY          :  ALLOCATE_L6_2

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_DLR'
      CHARACTER(  1*BYTE)             :: CLOSE_IT          ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER(  8*BYTE)             :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER(  1*BYTE)             :: EQUED             ! 'Y' if KLL stiff matrix was equilibrated in subr EQUILIBRATE
      CHARACTER( 24*BYTE)             :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN
      CHARACTER( 22*BYTE)             :: MODNAM1           ! Name to write to screen to describe module being run
      CHARACTER(  1*BYTE)             :: READ_NTERM        ! 'Y' or 'N' Input to subr READ_MATRIX_1
      CHARACTER(  1*BYTE)             :: NULL_COL          ! 'Y' if a col of KLR(transpose) is null
      CHARACTER(  1*BYTE)             :: OPND              ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SCRFIL            ! File name

      INTEGER(LONG)                   :: DEB_PRT(2)        ! Debug numbers to say whether to write ABAND and/or its decomp to output
!                                                            file in called subr SYM_MAT_DECOMP_LAPACK
      INTEGER(LONG)                   :: I,J               ! DO loop indices or counters
      INTEGER(LONG)                   :: INFO        = 0   ! Info on success of factorization or solve
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      REAL(DOUBLE)                    :: EQUIL_SCALE_FACS(NDOFL)
                                                           ! LAPACK_S values returned from subr SYM_MAT_DECOMP_LAPACK

      REAL(DOUBLE)                    :: DLR_COL(NDOFL)    ! A column of DLR solved for herein
      REAL(DOUBLE)                    :: INOUT_COL(NDOFL)  ! Temp variable for subr FBS
      REAL(DOUBLE)                    :: K_INORM           ! Inf norm of KLL matrix (det in  subr COND_NUM)
      REAL(DOUBLE)                    :: RCOND             ! Recrip of cond no. of the KLL. Det in  subr COND_NUM

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Decomp KLL

      DEB_PRT(1) = 64
      DEB_PRT(2) = 65

      IF (SOLLIB == 'BANDED') THEN

         INFO  = 0
         EQUED = 'N'
         CALL SYM_MAT_DECOMP_LAPACK ( SUBR_NAME, 'KLL', 'L ', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, 'Y', 'N', 'N', 'N', DEB_PRT, &
                                      EQUED, KLL_SDIA, K_INORM, RCOND, EQUIL_SCALE_FACS, INFO )
         IF (EQUED == 'Y') THEN                         ! If EQUED == 'Y' then error. We don't want KLL equilibrated
            WRITE(ERR,6001) SUBR_NAME, EQUED
            WRITE(F06,6001) SUBR_NAME, EQUED
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ELSE IF (SOLLIB == 'SPARSE  ') THEN

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            INFO = 0
            CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KLL', 'L ', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, INFO )

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

      DO I=1,NDOFL                                         ! Make sure that scale factors are one. We don't want any equil scaling
         EQUIL_SCALE_FACS(I) = ONE                         ! of KLL in this subr. FBS below has EQUIL_SCALE_FACS as an input but
      ENDDO                                                ! they shouldn't be used as EQUED = 'N' is also input there (1st arg)

!***********************************************************************************************************************************
! Solve for DLR

! Open a scratch file that will be used to write DLR nonzero terms to as we solve for columns of DLR. After all col's
! of DLR have been solved for, and we have a count on NTERM_DLR, we will allocate memory to the DLR arrays and read
! the scratch file values into those arrays. Then, in the calling subroutine, we will write NTERM_DLR, followed by
! DLR  row/col/value to a permanent file

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

! Loop on columns of KLR using rows of KRL

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      NTERM_DLR = 0
      CALL COUNTER_INIT('    Solve for DLR col ', NDOFR)
      DO J = 1,NDOFR

! To solve for the j-th col of DLR, use the j-th row of KRL (j-th col of KLR) as a rhs vector. Get the j-th row of KRL and put
! the negative of it into array INOUT_COL:

         NULL_COL = 'Y'
         DO I=1,NDOFL
            INOUT_COL(I) = ZERO
            DLR_COL(I)   = ZERO
         ENDDO
         CALL GET_SPARSE_CRS_ROW ( 'KRL', J,  NTERM_KRL, NDOFR, NDOFL, I_KRL, J_KRL, KRL, -ONE, INOUT_COL, NULL_COL )

! Calculate DLR_COL via forward/backward substitution.

         IF (NULL_COL == 'N') THEN                         ! FBS will solve for DLR_COL & load it into DLR array
                                                           ! DPBTRS will return DLR_COL = -KLL(-1)*RHS_col
!                                                            Note 1st arg = 'N' assures that EQUIL_SCAL_FACS will not be used
            IF      (SOLLIB == 'BANDED') THEN

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
               DLR_COL(I) = INOUT_COL(I)
            ENDDO
            DO I=1,NDOFL                                   ! Count NTERM_DLR and write nonzero DLR to scratch file
               IF (DABS(DLR_COL(I)) > EPS1) THEN
                  NTERM_DLR = NTERM_DLR + 1
                  WRITE(SCR(1)) I, J, DLR_COL(I)
               ENDIF
            ENDDO
         ENDIF
         CALL COUNTER_PROGRESS(J)
      ENDDO

! The DLR data in SCRATCH-991 is written one col at a time for DLR. Therefore it is rows of DLRt

      REWIND (SCR(1))
      MESSAG = 'SCRATCH: DLR ROW/COL/VAL'
      READ_NTERM = 'N'
      OPND       = 'Y'
      CLOSE_IT   = 'N'
      CLOSE_STAT = 'KEEP    '

      CALL ALLOCATE_SPARSE_MAT ( 'DLRt', NDOFR, NTERM_DLR, SUBR_NAME )
      CALL ALLOCATE_L6_2 ( 'DLRt', SUBR_NAME )

      CALL ALLOCATE_SPARSE_MAT ( 'DLR', NDOFL, NTERM_DLR, SUBR_NAME )
      CALL ALLOCATE_L6_2 ( 'DLR', SUBR_NAME )
                                                           ! J_DLRt is same as I2_DLR and I2_DLRt is same as J_DLR
      CALL READ_MATRIX_2 ( SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG, 'DLRt', NDOFR, NTERM_DLR, READ_NTERM,               &
                           J_DLRt, I2_DLRt, DLRt)

! Now get DLR from DLRt

      CALL GET_I_MAT_FROM_I2_MAT ( 'DLRt', NDOFR, NTERM_DLR, I2_DLRt, I_DLRt )

      CALL MATTRNSP_SS ( NDOFR, NDOFL, NTERM_DLR, 'DLRt', I_DLRt, J_DLRt, DLRt, 'DLR', I_DLR, J_DLR, DLR )

      CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )

! Print out constraint matrix DLR, if requested

      IF ( PRTDLR == 1) THEN
         IF (NTERM_DLR > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CB BOUNDARY MODE MATRIX DLR', 'L ', 'R ', NTERM_DLR, NDOFL, I_DLR, J_DLR, DLR )
         ENDIF
      ENDIF



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

      END SUBROUTINE SOLVE_DLR

   END MODULE CB_LOAD_TRANSFORMATIONS
