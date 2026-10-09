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

   MODULE SPARSE_MATRIX_ALGEBRA

   USE SPARSE_CRS_ACCESS, ONLY :  ROW_AT_COLJ_BEGEND

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM, MATTRNSP_SS

   CONTAINS

      SUBROUTINE MATADD_SSS ( NROWS, MAT_A_NAME, NTERM_A, I_A, J_A, A, ALPHA,                          &
                                     MAT_B_NAME, NTERM_B, I_B, J_B, B, BETA,                           &
                                     MAT_C_NAME, NTERM_C, I_C, J_C, C )

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! Subroutine MATADD_SSS_NTERM must be run before this subroutine to calculate NTERM_C, an input to this subroutine, that is the
! number of nonzero terms in C. Then memory can be allocated to C before this subroutine is called
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

! Performs the matrix operation C = A + B. Input matrices A and B are in sparse (compressed row storage) format.
! Output matrix C is in sparse format.

! Matrices A and B must have the same number of rows and cols. Ensuring this is the responsibility of the user.

! NOTE: Both of the 2 input matrices, as well as the matrix resulting from the addition, must be stored the same as far as
! symmetry is concerned. That is, if A (or B) is a symmetric matrix and is stored with only the terms on and above the diagonal,
! then both input matrices must be stored with only the terms on and above the diagonal. Matrix C will, by implication, be
! symmetric and have only terms on and above its diagonal in array C. Thus, this subr cannot add 2 matrices where one is stored
! symmetric and the other is not. The user is required to ensure that this is the case.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATADD_SSS'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME        ! Name of matrix A
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_B_NAME        ! Name of matrix B
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_C_NAME        ! Name of matrix C

      INTEGER(LONG), INTENT(IN )      :: NROWS             ! Number of rows in input matrices A and B
      INTEGER(LONG), INTENT(IN )      :: NTERM_A           ! Number of nonzero terms in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_B           ! Number of nonzero terms in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NTERM_C           ! Number of nonzero terms in output matrix C
      INTEGER(LONG), INTENT(IN )      :: I_A(NROWS+1)      ! I_A(I+1) - I_A(I) = no. terms in row I of matrix A
      INTEGER(LONG), INTENT(IN )      :: I_B(NROWS+1)      ! I_B(I+1) - I_B(I) = no. terms in row I of matrix B
      INTEGER(LONG), INTENT(IN )      :: J_A(NTERM_A)      ! Col no's for nonzero terms in matrix A
      INTEGER(LONG), INTENT(IN )      :: J_B(NTERM_B)      ! Col no's for nonzero terms in matrix B
      INTEGER(LONG), INTENT(OUT)      :: I_C(NROWS+1)      ! I_C(I+1) - I_C(I) = no. terms in row I of matrix C
      INTEGER(LONG), INTENT(OUT)      :: J_C(NTERM_C)      ! Col no's for nonzero terms in matrix C


      REAL(DOUBLE) , INTENT(IN )      :: A(NTERM_A)        ! Nonzero terms in matrix A
      REAL(DOUBLE) , INTENT(IN )      :: B(NTERM_B)        ! Nonzero terms in matrix B
      REAL(DOUBLE) , INTENT(IN )      :: ALPHA             ! Scalar multiplier for matrix A
      REAL(DOUBLE) , INTENT(IN )      :: BETA              ! Scalar multiplier for matrix B
      REAL(DOUBLE) , INTENT(OUT)      :: C(NTERM_C)        ! Nonzero terms in matrix C

      INTEGER(LONG)                   :: ROW
      INTEGER(LONG)                   :: P_A
      INTEGER(LONG)                   :: P_B
      INTEGER(LONG)                   :: COL_A
      INTEGER(LONG)                   :: COL_B
      INTEGER(LONG)                   :: CNT
      REAL(DOUBLE)                    :: V




! **********************************************************************************************************************************

      CNT = 0
      I_C(1) = 1

      DO ROW=1,NROWS
         P_A = I_A(ROW)
         P_B = I_B(ROW)

         DO WHILE(P_A < I_A(ROW+1) .OR. P_B < I_B(ROW+1))

                                                           ! Sentinel when A's row is exhausted
            IF (P_A < I_A(ROW+1)) then
               COL_A = J_A(P_A)
            ELSE
               COL_A = HUGE(0)
            ENDIF
                                                           ! Sentinel when B's row is exhausted
            IF (P_B < I_B(ROW+1)) then
               COL_B = J_B(P_B)
            ELSE
               COL_B = HUGE(0)
            ENDIF


            IF (COL_A < COL_B) THEN                        ! Only A has an entry in this column
               CNT = CNT + 1
               C(CNT) = ALPHA * A(P_A)
               J_C(CNT) = COL_A
               P_A = P_A + 1
            ELSE IF (COL_B < COL_A) THEN                   ! Only B has an entry in this column
               CNT = CNT + 1
               C(CNT) = BETA * B(P_B)
               J_C(CNT) = COL_B
               P_B = P_B + 1
            ELSE                                           ! Both have an entry -- add
               V = ALPHA * A(P_A) + BETA * B(P_B)
               CNT = CNT + 1
               C(CNT) = V
               J_C(CNT) = COL_A
               P_A = P_A + 1
               P_B = P_B + 1
            ENDIF

         ENDDO

         I_C(ROW+1) = CNT + 1

      ENDDO




      RETURN

! **********************************************************************************************************************************


      END SUBROUTINE MATADD_SSS


      SUBROUTINE MATADD_SSS_NTERM ( NROWS, MAT_A_NAME, NTERM_A, I_A, J_A, SYM_A, MAT_B_NAME, NTERM_B, I_B, J_B, SYM_B,             &
                                           MAT_C_NAME, NTERM_C )

! Setup routine for performing the sparse matrix add operation C = A + B where A, B and C are in stored in sparse CRS format.
! This subr must be run prior to the subr that actually does the add (MATADD_SSS) in order to calc NTERM_C, the number of terms that
! will be in C (so that memory could be allocated, prior to this MATADD_SSS, for arrays J_C and C)

! Matrices A and B must have the same number of rows and cols. Ensuring this is the responsibility of the user.

! NOTE: Both of the 2 input matrices, as well as the matrix resulting from the addition, must be stored the same as far as
! symmetry is concerned. That is, if A (or B) is a symmetric matrix and is stored with only the terms on and above the diagonal,
! then both input matrices must be stored with only the terms on and above the diagonal. Matrix C will, by implication, be
! symmetric and have only terms on and above its diagonal in array C. Thus, this subr cannot add 2 matrices where one is stored
! symmetric and the other is not. The user is required to ensure that this is the case.

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATADD_SSS_NTERM'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME        ! Name of matrix A
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_B_NAME        ! Name of matrix B
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_C_NAME        ! Name of matrix C
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_A             ! Flag for whether matrix A is stored sym (terms on and above diag)
!                                                            or nonsym (all terms)
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_B             ! Flag for whether matrix B is stored sym (terms on and above diag)
!                                                            or nonsym (all terms)

      INTEGER(LONG), INTENT(IN )      :: NROWS             ! Number of rows in input matrices A and B
      INTEGER(LONG), INTENT(IN )      :: NTERM_A           ! Number of nonzero terms in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_B           ! Number of nonzero terms in input matrix B
      INTEGER(LONG), INTENT(IN )      :: I_A(NROWS+1)      ! I_A(I+1) - I_A(I) = no. terms in row I of matrix A
      INTEGER(LONG), INTENT(IN )      :: I_B(NROWS+1)      ! I_B(I+1) - I_B(I) = no. terms in row I of matrix B
      INTEGER(LONG), INTENT(IN )      :: J_A(NTERM_A)      ! Col no's for nonzero terms in matrix A
      INTEGER(LONG), INTENT(IN )      :: J_B(NTERM_B)      ! Col no's for nonzero terms in matrix B
      INTEGER(LONG), INTENT(OUT)      :: NTERM_C           ! Number of nonzero terms in output matrix C



      INTEGER(LONG)                   :: ROW
      INTEGER(LONG)                   :: P_A
      INTEGER(LONG)                   :: P_B
      INTEGER(LONG)                   :: COL_A
      INTEGER(LONG)                   :: COL_B
      INTEGER(LONG)                   :: CNT




! **********************************************************************************************************************************
! Make sure that input matrices A and B are stored in same format

      IF (SYM_A /= SYM_B) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,941) SUBR_NAME, MAT_A_NAME, MAT_B_NAME, MAT_A_NAME, SYM_A, MAT_B_NAME, SYM_B
         WRITE(F06,941) SUBR_NAME, MAT_A_NAME, MAT_B_NAME, MAT_A_NAME, SYM_A, MAT_B_NAME, SYM_B
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      CNT = 0

      DO ROW=1,NROWS
         P_A = I_A(ROW)
         P_B = I_B(ROW)

         DO WHILE(P_A < I_A(ROW+1) .OR. P_B < I_B(ROW+1))

                                                           ! Sentinel when A's row is exhausted
            IF (P_A < I_A(ROW+1)) then
               COL_A = J_A(P_A)
            ELSE
               COL_A = HUGE(0)
            ENDIF
                                                           ! Sentinel when B's row is exhausted
            IF (P_B < I_B(ROW+1)) then
               COL_B = J_B(P_B)
            ELSE
               COL_B = HUGE(0)
            ENDIF


            IF (COL_A < COL_B) THEN                        ! Only A has an entry in this column
               CNT = CNT + 1
               P_A = P_A + 1
            ELSE IF (COL_B < COL_A) THEN                   ! Only B has an entry in this column
               CNT = CNT + 1
               P_B = P_B + 1
            ELSE                                           ! Both have an entry
               CNT = CNT + 1
               P_A = P_A + 1
               P_B = P_B + 1
            ENDIF

         ENDDO

      ENDDO

      NTERM_C = CNT




      RETURN

! **********************************************************************************************************************************
  941 FORMAT(' *ERROR   941: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INPUT MATRICES ',A,' AND ',A,' MUST BOTH BE STORED IN THE SAME FORMAT (SYM MUST BE BOTH "Y" OR "N").' &
                    ,/,14X,' HOWEVER, MATRIX ',A,' HAS SYM = ',A,' AND MATRIX ',A,' HAS SYM = ',A)


! **********************************************************************************************************************************

      END SUBROUTINE MATADD_SSS_NTERM


      SUBROUTINE MATMULT_SSS ( MAT_A_NAME, NROW_A, NTERM_A, SYM_A, I_A, J_A, A,                                                    &
                               MAT_B_NAME, NCOL_B, NTERM_B, SYM_B, J_B, I_B, B, AROW_MAX_TERMS, MAT_C_NAME, CONS,                  &
                                                   NTERM_C,        I_C, J_C, C )

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! Subroutine MATMULT_SSS_NTERM must be run before this subroutine to calculate NTERM_C, an input to this subroutine, that is the
! number of nonzero terms in C. Then memory can be allocated to C before this subroutine is called
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

! Sparse matrix multiply to obtain C = cons*A*B with A, B and C in sparse format. A and C are stored in compressed row storage
! (CRS) format. In addition, if A is symmetric it can be stored  with only the terms on, and above, the diagonal.
! Matrix B must be stored in compressed col storage (CCS) and cannot be stored symmetric (i.e. all terms in B must be stored)

! Input  matrix A is stored in compressed row storage (CRS) format using arrays I_A(NROW_A+1), J_A(NTERM_A) A(NTERM_A) where NROW_A
! is the number of rows in matrix A and NTERM_A are the number of nonzero terms in matrix A:

!      I_A is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          I_A(I+1) - I_A(I) are the number of nonzero terms in row I of matrix A

!      J_A is an integer array giving the col numbers of the NTERM_A nonzero terms in matrix A

!        A is a real array of the nonzero terms in matrix A. If SYM_A='Y' then only the terms on, and above, the diag are stored.

! Input  matrix B is stored in compressed col storage (CCS) format using arrays J_B(NCOL_B+1), I_B(NTERM_B), B(NTERM_B) where NCOL_B
! is the number of columns in matrix B and NTERM_B are the number of nonzero terms in matrix B:

!      J_B is an array of NCOL_B+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          J_B(I+1) - J_B(I) are the number of nonzero terms in col I of matrix B

!      I_B is an integer array giving the row numbers of the NTERM_B nonzero terms in matrix B

!        B is a real array of the nonzero terms in matrix B. All of the nonzero terms of B must be included (i.e., no SYM option)

! Output matrix C, which must have the same number of rows as matrix A and the same number of columns as matrix B
! is stored in compressed row storage (CRS) format using arrays I_C(NROW_A+1), J_C(NTERM_C) C(NTERM_C) where NROW_A is
! the number of rows in matrix A and NTERM_C are the number of nonzero terms in matrix C:

!      I_C is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix C. That is:
!          I_C(I+1) - I_C(I) are the number of nonzero terms in row I of matrix C

!      J_C is an integer array giving the col numbers of the NTERM_C nonzero terms in matrix C

!        C is a real array of the nonzero terms in matrix C. All of the nonzero terms of C are stored (i.e., no SYM option for C)

! This subr determines integer arrays I_C and J_C and real array C.

! In order to handle symmetric A matrices, which do not have all terms in a row (only those from the diagonal out), array AROW is
! used. AROW is a 1D array that contains all nonzero terms in one row of A (including those that are not explicitly in array A due
! to symmetry). The multiplication of matrix A times matrix B is then accomplished (row by row of result matrix C) by multiplying
! AROW times matrix B. AROW is a compact array containing only the nonzero terms from one row of matrix A. Thus, integer array
! J_AROW is needed to give the column numbers, from matrix A (for one row), that the terms in array AROW are for.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_CRS_ACCESS, ONLY     :  ROW_AT_COLJ_BEGEND
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULT_SSS'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME            ! Name of matrix A
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_B_NAME            ! Name of matrix B
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_C_NAME            ! Name of matrix C
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_A                 ! ='Y' if matrix A is input symmetric (terms on and above diag only)
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_B                  ! ='Y' if matrix B is input sym (terms on and above diag only)

      INTEGER(LONG), INTENT(IN )      :: AROW_MAX_TERMS        ! Max number of terms in any row of A
      INTEGER(LONG), INTENT(IN )      :: NCOL_B                ! Number of cols in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NROW_A                ! Num rows in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_A               ! Num non 0's in input  matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_B               ! Num non 0's in input  matrix B
      INTEGER(LONG), INTENT(IN )      :: NTERM_C               ! Size of arrays J_C and C (MUST be determined by subr MATMULT_SSS)
      INTEGER(LONG), INTENT(IN )      :: I_A(NROW_A+1)         ! I_A(I+1) - I_A(I) = num nonzeros in row I of matrix A (CRS format)
      INTEGER(LONG), INTENT(IN )      :: J_A(NTERM_A)          ! Col no's for nonzero terms in matrix A
      INTEGER(LONG), INTENT(IN )      :: J_B(NCOL_B+1)         ! J_B(I+1) - J_B(I) = num nonzeros in col I of matrix B (CCS format)
      INTEGER(LONG), INTENT(IN )      :: I_B(NTERM_B)          ! Row no's for nonzero terms in matrix B
      INTEGER(LONG), INTENT(OUT)      :: I_C(NROW_A+1)         ! I_C(I+1) - I_C(I) = num nonzeros in row I of matrix C (CRS format)
      INTEGER(LONG), INTENT(OUT)      :: J_C(NTERM_C)          ! Col no's for nonzero terms in matrix C
      INTEGER(LONG)                   :: A_ROW_BEG             ! Index into array A where a row of matrix A begins
      INTEGER(LONG)                   :: A_ROW_END             ! Index into array A where a row of matrix A ends
      INTEGER(LONG)                   :: A_COL_NUM             ! A col number in matrix A
      INTEGER(LONG)                   :: B_COL_BEG             ! Index into array B where a col of matrix B begins
      INTEGER(LONG)                   :: B_COL_END             ! Index into array B where a col of matrix B ends
      INTEGER(LONG)                   :: B_ROW_NUM             ! A row number in matrix B
      INTEGER(LONG)                   :: A_NTERM_ROW_I         ! Number of terms in row I of matrix A
      INTEGER(LONG)                   :: B_NTERM_COL_J         ! Number of terms in col J of matrix B
      INTEGER(LONG)                   :: DELTA_KTERM_C         ! Incr in KTERM_C (0 or 1) when mult row I of A times col J of B
      INTEGER(LONG)                   :: I,J,K,L,II            ! DO loop indices
      INTEGER(LONG)                   :: I1,I2                 ! DO loop range
      INTEGER(LONG)                   :: J_AROW(AROW_MAX_TERMS)! Col numbers of terms in real array AROW (see below)
      INTEGER(LONG)                   :: KTERM_C               ! Count of number of nonzero terms put into output matrix C
      INTEGER(LONG)                   :: NHITS                 ! Number of "hits" of terms in a row of A existing where terms in a
!                                                                col of B exist when a row of A is multiplied by a col of B
      INTEGER(LONG)                   :: NTERM_AROW            ! Number of nonzero terms in AROW (one row of A)
      INTEGER(LONG)                   :: A_ROW_COLJ_BEG(NROW_A)! jth term is row number in array A where col j nonzeros begin
      INTEGER(LONG)                   :: A_ROW_COLJ_END(NROW_A)! jth term is row number in MATIN where col j nonzeros end


      REAL(DOUBLE) , INTENT(IN )      :: CONS                  ! Constant multiplier in cons*A*B to get C
      REAL(DOUBLE) , INTENT(IN )      :: A(NTERM_A)            ! Nonzero values in matrix A
      REAL(DOUBLE) , INTENT(IN )      :: B(NTERM_B)            ! Nonzero values in matrix B
      REAL(DOUBLE) , INTENT(OUT)      :: C(NTERM_C)            ! Nonzero values in matrix C
      REAL(DOUBLE)                    :: CTEMP                 ! A value accumulated as the nonzero terms from one row of A are
!                                                                multiplied by the corresponding nonzero terms from one col of B
      REAL(DOUBLE)                    :: AROW(AROW_MAX_TERMS)  ! Array containing the nonzero terms from one row of A

      INTRINSIC                       :: MAX



! **********************************************************************************************************************************
! Make sure B is stored as nonsym (all terms)

      IF (SYM_B /= 'N') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,948) SUBR_NAME, MAT_B_NAME, SYM_B
         WRITE(F06,948) SUBR_NAME, MAT_B_NAME, SYM_B
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Initialize outputs

      DO I=1,NROW_A+1
         I_C(I) = 0
      ENDDO

      DO I=1,NTERM_C
         J_C(I) = 0
           C(I) = ZERO
      ENDDO

      IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '1', '   ' )

! Create arrays that give the row numbers at which col j begins and ends. When SYM_A = 'Y'', we have to get
! terms for the A matrix that are not explicitly in A. This is done by getting A terms in the column (above the
! diagonal of a row) as well as the explicit terms from A that are there from the diagonal out to the end of the row. The two
! arrays: A_ROW_COLJ_BEG and A_ROW_COLJ_END are used to aid in getting the terms in the column above the diagonal.
! A_ROW_COLJ_BEG is an array that gives, for each col of A, the starting row number of nonzero terms in that column.
! A_ROW_COLJ_END is an array that gives, for each col of A, the ending   row number of nonzero terms in that column.
! The span: A_ROW_COLJ_BEG to A_ROW_COLJ_END is used when we search for terms in the columns.
! We only need A_ROW_COLJ_BEG and A_ROW_COLJ_END when A is input symmetric and MATOUT is not to be output as symmetric

      IF (SYM_A == 'Y') THEN                              ! The number of cols in A is NROW_A (due to SYM_A = 'Y')

         CALL ROW_AT_COLJ_BEGEND ( MAT_A_NAME, NROW_A, NROW_A, NTERM_A, I_A, J_A, A_ROW_COLJ_BEG, A_ROW_COLJ_END )

      ENDIF

! Do the multiply, using values put into AROW and J_AROW for each row of A. This is done to facilitate the SYM option for matrix A

      DELTA_KTERM_C = 0                                    ! Initialize variables used in the matrix multiplication in the DO I loop
      KTERM_C       = 0
      NHITS         = 0
      A_ROW_BEG        = 1
      CTEMP         = ZERO
      I_C(1)        = 1
i_do: DO I=1,NROW_A                                        ! Matrix multiply loop. Range over the rows in A

         I_C(I+1) = I_C(I)                                 ! End value in I_C for next this row is initially set at beginning value
         A_NTERM_ROW_I = I_A(I+1) - I_A(I)                 ! Number of terms in matrix A in row I
         IF (A_NTERM_ROW_I == 0) CYCLE i_do
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1         ! A_ROW_BEG to A_ROW_END is range of indices of terms in A for row I of A
         IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '2', '   ' )

         DO K=1,AROW_MAX_TERMS                             ! Null J_AROW and AROW each time we begin a new row of A
            AROW(K)   = ZERO
            J_AROW(K) = 0
         ENDDO
                                                           ! Build AROW for one row of matrix A
         NTERM_AROW = 0                                    ! 1st, look for terms that would be in this row, but are not, due to SYM
         IF (SYM_A == 'Y') THEN
            DO K=1,A_ROW_BEG-1
               IF (J_A(K) == I) THEN
                  NTERM_AROW         = NTERM_AROW + 1
                  AROW(NTERM_AROW)   = A(K)
                  I1 = A_ROW_COLJ_BEG(I)
                  I2 = A_ROW_COLJ_END(I)
                  DO II=I1,I2
                     IF ((K >= I_A(II)) .AND. (K < I_A(II+1))) THEN
                        J_AROW(NTERM_AROW) = II
                        IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '5', ' #2' )
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

         DO K=A_ROW_BEG,A_ROW_END                          ! 2nd, get terms from this row of A from the diagonal out
            NTERM_AROW = NTERM_AROW + 1
            AROW(NTERM_AROW)   = A(K)
            J_AROW(NTERM_AROW) = J_A(K)
            IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '6', ' #1' )
         ENDDO

         B_COL_BEG = 1
j_do:    DO J=1,NCOL_B                                     ! J loops over the number of columns in B
            B_NTERM_COL_J = J_B(J+1) - J_B(J)
            IF (B_NTERM_COL_J == 0) CYCLE j_do
            B_COL_END = B_COL_BEG + B_NTERM_COL_J - 1      ! B_COL_BEG to B_COL_END is range of indices of terms in B for col J of B

k_do:       DO K=1,NTERM_AROW                              ! The following 2 loops produce the ij-th term of C
               A_COL_NUM = J_AROW(K)
l_do:          DO L=B_COL_BEG,B_COL_END
                  B_ROW_NUM = I_B(L)
                  IF (A_COL_NUM == B_ROW_NUM) THEN
                     NHITS = NHITS + 1
                     DELTA_KTERM_C = 1
                     CTEMP = CTEMP + CONS*AROW(K)*B(L)     ! This is the nonzero ij-th term in matrix C.
                  ENDIF
               ENDDO l_do
            ENDDO k_do
            B_COL_BEG   = B_COL_END + 1

            KTERM_C  = KTERM_C + DELTA_KTERM_C             ! Now update sparse CRS representation of C
            IF (KTERM_C > NTERM_C) CALL ARRAY_SIZE_ERROR_1( SUBR_NAME, NTERM_C, MAT_C_NAME )
             IF (NHITS > 0) THEN
               I_C(I+1) = I_C(I+1) + 1
               J_C(KTERM_C) = J
                 C(KTERM_C) = CTEMP
               IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '7', '   ' )
            ENDIF
            DELTA_KTERM_C = 0
            CTEMP = ZERO
            NHITS = 0

         ENDDO j_do
         A_ROW_BEG = A_ROW_END + 1
         IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) THEN
            WRITE(F06,*)
         ENDIF

      ENDDO i_do

      IF ((DEBUG(84) == 2) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_DEB ( '9', '   ' )



      RETURN

! **********************************************************************************************************************************
  948 FORMAT(' *ERROR   948: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX B = ',A,' MUST BE STORED AS SYM_B = ','''N''',' FOR THIS SUBR TO WORK BUT IS STORED AS',   &
                           ' SYM_B = ',''',A,''')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MATMULT_SSS_DEB ( WHICH, ALG )

      CHARACTER(LEN=*), INTENT(IN)    :: ALG                    ! Which algorithm is used (#1 for terms above diag when SYM_A='Y'
!                                                                 or #2 for terms in row from diag out)
      CHARACTER( 1*BYTE)              :: WHICH                  ! Decides what to print out for this call to this subr
      INTEGER(LONG)                   :: I,J,K                  ! Local loop indices

! **********************************************************************************************************************************
      IF      (WHICH == '1') THEN

         WRITE(F06,*)
         WRITE(F06,1011)
         WRITE(F06,1012)
         WRITE(F06,1013)
         WRITE(F06,1014) MAT_A_NAME, MAT_B_NAME, MAT_C_NAME
         WRITE(F06,1015) CONS
         WRITE(F06,1016) NROW_A, NTERM_A, NCOL_B, NTERM_B, NTERM_C
         IF (SYM_A == 'Y') THEN
            WRITE(F06,1017)
         ELSE
            WRITE(F06,1018)
         ENDIF
         WRITE(F06,1019)
         WRITE(F06,*)

      ELSE IF (WHICH == '2') THEN

         WRITE(F06,1021)
         WRITE(F06,1022) I
         WRITE(F06,1023) I, I, A_ROW_BEG, A_ROW_END
         WRITE(F06,1024)
         WRITE(F06,*)

      ELSE IF (WHICH == '3') THEN

      ELSE IF (WHICH == '4') THEN

      ELSE IF (WHICH == '5') THEN

         WRITE(F06,1051) ALG,                     K , II, J_A(K) , A(K), NTERM_AROW, J_AROW(NTERM_AROW), AROW(NTERM_AROW)

      ELSE IF (WHICH == '6') THEN

         WRITE(F06,1061) ALG, K, I, J_A(K), A(K),                         NTERM_AROW, J_AROW(NTERM_AROW), AROW(NTERM_AROW)

      ELSE IF (WHICH == '7') THEN

         IF (NHITS > 0) THEN
            IF (J == 1) THEN
               WRITE(F06,*)
               WRITE(F06,1071) I
            ENDIF
            WRITE(F06,1072) I, J, NHITS, KTERM_C, I, J_C(KTERM_C), C(KTERM_C)
         ENDIF

      ELSE IF (WHICH == '8') THEN

      ELSE IF (WHICH == '9') THEN

         WRITE(F06,*)
         WRITE(F06,1091) MAT_C_NAME
         DO I=1,NROW_A+1                                   ! The number of rows in C is the same as that in A
            WRITE(F06,9192) I,I_C(I)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,1093)
         DO K=1,NTERM_C
            WRITE(F06,1094) K, J_C(K), C(K)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,1095)

      ENDIF

! **********************************************************************************************************************************
 1011 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' :::::::::::::::::::::::::::::::::::::::START DEBUG(84) OUTPUT FROM SUBROUTINE MATMULT_SSS:::::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1012 FORMAT(' SSS SPARSE MATRIX MULTIPLY ROUTINE: Multiply matrix A, stored in sparse Compressed Row Storage (CRS) format, times',&
' matrix B, stored in',/,' -----------------------------------',/,' sparse Compressed Column Storage (CCS) format, to obtain'     ,&
' sparse CRS matrix C',/)

 1013 FORMAT(' A may be stored as symmetric (only terms on and above the diagonal) or with all nonzero terms included.'         ,/,&
' Matrix B must be stored with all nonzero terms. Result matrix C will also be stored with all nonzero terms',/)

 1014 FORMAT(40X,' The name of CRS formatted matrix A is: ',A                                                                   ,/,&
             40x,' The name of CRS formatted matrix B is: ',A                                                                   ,/,&
             40x,' The name of CRS formatted matrix C is: ',A,/)

 1015 FORMAT(' Multiply ',1ES14.6,' times the product of matrix A and matrix B to obtain matrix C',/)

 1016 FORMAT(36X,' Matrix A has ',I8,' rows and '  ,I12,' nonzero terms'                                                        ,/,&
             36X,' Matrix B has ',I8,' cols and '  ,I12,' nonzero terms'                                                        ,/,&
             36X,' Matrix C will have             ',I12,' nonzero terms*'                                                       ,/,&
             22X,'*(as detrmined by subr MATMULT_SSS_NTERM which had to have been run prior to this subr)'/)

 1017 FORMAT(' Matrix A was input as a symmetric CRS array (only those terms on and above the diagonal are stored in array A)',/)

 1018 FORMAT(' Matrix A was input as a CRS array with all nonzero terms stored in array A',/)

 1019 FORMAT(                                                                                                                      &
' In order to handle symmetric A matrices, which do not have all terms in a row (only those from the diagonal out), arrays AROW',/,&
' and J_AROW are used. AROW is a 1D array that contains all nonzero terms from one row of A (including those that are not'      ,/,&
' explicitly in array A due to symmetry). The multiplication of matrix A times matrix B is then accomplished (row by row of'    ,/,&
' result matrix C) by multiplying AROW times matrix B. Since AROW is a compact array containing only the nonzero terms from one',/,&
' row of matrix A, integer array J_AROW is also needed to give the col numbers, from matrix A (for one row), for the terms in'  ,/,&
' array AROW.'                                                                                                                  ,//&
' Alg #1 (below) gets data for arrays J_AROW and AROW directly from the Compressed Row Storage (CRS) format of array A'         ,/,&
' Alg #2 (below) is only needed if matrix A is input as symmetric (only terms on and above the diagonal) and gets terms for'    ,/,&
'         J_AROW and AROW from column I of matrix A while working on row I of matrix A. These are the terms that would be below',/,&
'         the diag in matrix A but are not explicitly in the array due to symmetry storage'                                    ,//,&
' For each row of matrix A, the following shows the development of arrays J_AROW and AROW and the result of mult AROW times'   ,   &
' matrix B to get one row of result matrix C. Output is only given for non null rows of matrix A and non null cols of matrix B',/)

 1021 FORMAT(' ******************************************************************************************************************',&
              '*****************')

 1022 FORMAT(30X,' W O R K I N G   O N   R O W ',I8,'   O F   O U T P U T   M A T R I X   C',/)

 1023 FORMAT(' Multiply row ',I8,' of matrix A times all columns of matrix B to get row ',I8,' of matrix C'                     ,/,&
             ' This row of A begins in array A(K) at index K  = ',I8,' and ends at index K = ',I8,//)


 1024 FORMAT(16X,'Data from input array A             Data below diag of matrix A not in array A             Data for array AROW',/&
,9X,'----------------------------------------   ------------------------------------------      --------------------------------',/&
,' Alg     Index       Row     Col        Value         Index       Row     Col        Value          Index    Col No        Value'&
,/,13X,'K         I    J_A(K)       A(K)             K        II    J_A(K)       A(K)              M  J_AROW(M)      AROW(M)')

 1051 FORMAT(1X,A3,10X,10X,10X,15X    ,I10,I10,I10,1ES15.6,I11,I11,1ES16.6)

 1061 FORMAT(1X,A3,I10,I10,I10,1ES15.6,10X,10X,10X,15X    ,I11,I11,1ES16.6)

 1071 FORMAT('                                                                          Data for row',I8,' of output matrix C'  ,/,&
             '                                                                       --------------------------------------------' &
          ,/,'                                                                          Index       Row     Col        Value'   ,/,&
             '                                                                              K         I    J_C(K)       C(K)',/)

 1072 FORMAT(' Row ',I8,' of A times col ',I8,' of B gets ',I8,' hits and:   ',I10,I10,i10,1ES16.6)


 1091 FORMAT(' ******************************************************************************************************************',&
              '*****************'                                                                                               ,/,&
             ' SUMMARY: Compressed Row Storage (CRS) format of matrix C = ',A,':',/,' -------'                                 ,//,&
             ' 1) Index, L, and array I_C(L) for matrix C, where I_C(L+1) - I_C(L) is the number of nonzero terms in row L of',    &
             ' matrix C.',/,'    (also, I_C(L) is the index, K, in array C(K) where row L begins - up to, but not including, the', &
             ' last entry in I_C(L)).',/)

 9192 FORMAT('    L, I_C(L)       = ',2I12)

 1093 FORMAT(' 2) Index, K, and arrays J_C(K) and C(K). C(K) are the nonzeros in matrix C and J_C(K) is the col number in matrix',&
             ' C for term C(K).',/)

 1094 FORMAT('    K, J_C(K), C(K) = ',2I12,1ES15.6)

 1095 FORMAT(' ::::::::::::::::::::::::::::::::::::::END DEBUG(84) OUTPUT FROM SUBROUTINE MATMULT_SSS::::::::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE MATMULT_SSS_DEB

      END SUBROUTINE MATMULT_SSS


      SUBROUTINE MATMULT_SSS_NTERM ( MAT_A_NAME, NROW_A, NTERM_A, SYM_A, I_A, J_A,                                                 &
                                     MAT_B_NAME, NCOL_B, NTERM_B, SYM_B, J_B, I_B, AROW_MAX_TERMS, MAT_C_NAME, NTERM_C )

! Setup routine for subr MATMULT_SSS (which performs the matrix operation C = cons*A*B With A and B in sparse format
! Matrices A and C are stored in compressed row storage (CRS) format. In addition, if A is symmetric it can be stored  with only the
! terms on, and above, the diagonal. Matrix B must be stored in compressed col storage (CCS) format.

! Input  matrix A is stored in compressed row storage (CRS) format using arrays I_A(NROW_A+1), J_A(NTERM_A) A(NTERM_A) where NROW_A
! is the number of rows in matrix A and NTERM_A are the number of nonzero terms in matrix A:

!      I_A is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          I_A(I+1) - I_A(I) are the number of nonzero terms in row I of matrix A

!      J_A is an integer array giving the col numbers of the NTERM_A nonzero terms in matrix A

!        A is a real array of the nonzero terms in matrix A. If SYM_A='Y' then only the terms on, and above, the diag are stored.
!          array A is not used in this subr - it is used in the actual matrix multiplication routine.

! Input  matrix B is stored in compressed col storage (CCS) format using arrays J_B(NCOL_B+1), I_B(NTERM_B), B(NTERM_B) where NCOL_B
! is the number of columns in matrix B and NTERM_B are the number of nonzero terms in matrix B:

!      J_B is an array of NCOL_B+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          J_B(I+1) - J_B(I) are the number of nonzero terms in col I of matrix B

!      I_B is an integer array giving the row numbers of the NTERM_B nonzero terms in matrix B

!        B is a real array of the nonzero terms in matrix B. All of the nonzero terms of B must be included (i.e., no SYM option)
!          array B is not used in this subr - it is used in the actual matrix multiplication routine.

! This subr determines the storage required (NTERM_C) for array J_C (and C) from the data in arrays I_A, J_A, I_B and J_B. This
! information is used to allocate memory for arrays J_C and C prior to calling subr MATMULT_SSS so that it can it can do the sparse
! matrix multiply

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE SPARSE_ALG_ARRAYS, ONLY     :  J_AROW
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_ALGORITHM_LIFECYCLE, ONLY:  ALLOCATE_SPARSE_ALG, DEALLOCATE_SPARSE_ALG
      USE SPARSE_CRS_ACCESS, ONLY     :  ROW_AT_COLJ_BEGEND

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULT_SSS_NTERM'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME             ! Name of matrix A
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_B_NAME             ! Name of matrix B
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_C_NAME             ! Name of matrix C
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_A                  ! ='Y' if matrix A is input sym (terms on and above diag only)
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_B                  ! ='Y' if matrix B is input sym (terms on and above diag only)

      INTEGER(LONG), INTENT(IN )      :: NCOL_B                 ! Number of cols in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NROW_A                 ! Number of rows in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_A                ! Number of nonzero terms in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_B                ! Number of nonzero terms in input matrix B
      INTEGER(LONG), INTENT(IN )      :: I_A(NROW_A+1)          ! I_A(I+1) - I_A(I) = no. terms in row I of matrix A
      INTEGER(LONG), INTENT(IN )      :: J_A(NTERM_A)           ! Col no's for nonzero terms in matrix A
      INTEGER(LONG), INTENT(IN )      :: J_B(NCOL_B+1)          ! J_B(I+1) - J_B(I) = no. terms in row I of matrix B
      INTEGER(LONG), INTENT(IN )      :: I_B(NTERM_B)           ! Row no's for nonzero terms in matrix B
      INTEGER(LONG), INTENT(OUT)      :: AROW_MAX_TERMS         ! Max number of terms in any row of A
      INTEGER(LONG), INTENT(OUT)      :: NTERM_C                ! Number of nonzero terms in output matrix C
      INTEGER(LONG)                   :: A_ROW_COLJ_BEG(NROW_A) ! jth term is row number in A where col j nonzeros begin
      INTEGER(LONG)                   :: A_ROW_COLJ_END(NROW_A) ! jth term is row number in A where col j nonzeros end
      INTEGER(LONG)                   :: A_NTERM_ROW_I          ! Number of terms in row I of matrix A
      INTEGER(LONG)                   :: A_COL_NUM              ! A col number in matrix A
      INTEGER(LONG)                   :: A_ROW_BEG              ! Index into array I_A where a row of matrix A begins
      INTEGER(LONG)                   :: A_ROW_END              ! Index into array I_A where a row of matrix A ends
      INTEGER(LONG)                   :: B_NTERM_COL_J          ! Number of terms in col J of matrix B
      INTEGER(LONG)                   :: B_ROW_NUM              ! A row number in matrix B
      INTEGER(LONG)                   :: B_COL_BEG              ! Index into array J_B where a col of matrix B begins
      INTEGER(LONG)                   :: B_COL_END              ! Index into array J_B where a col of matrix B ends
      INTEGER(LONG)                   :: DELTA_NTERM_C = 0      ! Incr in NTERM_C (0,1) resulting from mult row of A times col of B
      INTEGER(LONG)                   :: I,J,K,L,II             ! DO loop indices
      INTEGER(LONG)                   :: I1,I2                  ! DO loop range
      INTEGER(LONG)                   :: NHITS                  ! Num of "hits" of terms in a col of A existing where terms in a
!                                                                 row of B exist when a row of A is multiplied by a col of B
      INTEGER(LONG)                   :: NHITS_TOT_FOR_ROW_OF_A ! Num of "hits" of terms in a col of A existing where terms in any
!                                                                 row of B exist when a row of A is multiplied by al cols of B
      INTEGER(LONG)                   :: NTERM_AROW             ! Max number of nonzero terms in one row of A


      INTRINSIC                       :: MAX



! **********************************************************************************************************************************
! Make sure B is stored as nonsym (all terms)

      IF (SYM_B /= 'N') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,948) SUBR_NAME, MAT_B_NAME, SYM_B
         WRITE(F06,948) SUBR_NAME, MAT_B_NAME, SYM_B
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Initialize outputs

      AROW_MAX_TERMS = 0
      NTERM_C        = 0

! First, set up array J_AROW so it can handle the row from A that has the most terms in it

      A_ROW_BEG         = 1
      DO I=1,NROW_A
         A_NTERM_ROW_I = I_A(I+1) - I_A(I)
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1
         NTERM_AROW = 0
         IF (SYM_A == 'Y') THEN
            DO K=1,A_ROW_BEG-1
               IF (J_A(K) == I) THEN
                  NTERM_AROW = NTERM_AROW + 1
               ENDIF
            ENDDO
         ENDIF
         DO K=A_ROW_BEG,A_ROW_END
            NTERM_AROW = NTERM_AROW + 1
         ENDDO
         IF (NTERM_AROW > AROW_MAX_TERMS) THEN
            AROW_MAX_TERMS = NTERM_AROW
         ENDIF
         A_ROW_BEG = A_ROW_END + 1
      ENDDO
      NTERM_AROW = 0
                                                           ! Allocate integer vector to hold the column no's of the terms in the row
      CALL ALLOCATE_SPARSE_ALG ( 'J_AROW', AROW_MAX_TERMS, 0, SUBR_NAME)

      IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '1', '   ' )

! Create arrays that give the row numbers at which col j begins and ends. When SYM_A = 'Y'', we have to get
! terms for the A matrix that are not explicitly in A. This is done by getting A terms in the column (above the
! diagonal of a row) as well as the explicit terms from A that are there from the diagonal out to the end of the row. The two
! arrays: A_ROW_COLJ_BEG and A_ROW_COLJ_END are used to aid in getting the terms in the column above the diagonal.
! A_ROW_COLJ_BEG is an array that gives, for each col of A, the starting row number of nonzero terms in that column.
! A_ROW_COLJ_END is an array that gives, for each col of A, the ending   row number of nonzero terms in that column.
! The span: A_ROW_COLJ_BEG to A_ROW_COLJ_END is used when we search for terms in the columns.
! We only need A_ROW_COLJ_BEG and A_ROW_COLJ_END when A is input symmetric and MATOUT is not to be output as symmetric

      IF (SYM_A == 'Y') THEN

         CALL ROW_AT_COLJ_BEGEND ( MAT_A_NAME, NROW_A, NROW_A, NTERM_A, I_A, J_A, A_ROW_COLJ_BEG, A_ROW_COLJ_END )

      ENDIF

! Now count the terms that go into C, using values put into AROW and J_AROW for each row of A. This is done to facilitate the
! SYM option for matrix A.

      NHITS   = 0                                          ! Initialize variables used in the matrix multiplication in the DO I loop
      NHITS_TOT_FOR_ROW_OF_A = 0
      A_ROW_BEG  = 1
i_do: DO I=1,NROW_A                                        ! Matrix multiply loop. Range over the rows in A

         A_NTERM_ROW_I = I_A(I+1) - I_A(I)                 ! Number of terms in matrix A in row I
         IF (A_NTERM_ROW_I == 0) CYCLE i_do
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1         ! A_ROW_BEG to A_ROW_END is the range of indices of terms in A for row I
         IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '2', '   ' )

         DO K=1,AROW_MAX_TERMS                             ! Null J_AROW and AROW each time we begin a new row of A
            J_AROW(K) = 0
         ENDDO
                                                           ! Formulate J_AROW
         NTERM_AROW = 0                                    ! 1st, look for terms that would be in this row, but are not, due to SYM
         IF (SYM_A == 'Y') THEN
            DO K=1,A_ROW_BEG-1
               IF (J_A(K) == I) THEN
                  NTERM_AROW = NTERM_AROW + 1
                  I1 = A_ROW_COLJ_BEG(I)
                  I2 = A_ROW_COLJ_END(I)
                  DO II=I1,I2
                     IF ((K >= I_A(II)) .AND. (K < I_A(II+1))) THEN
                        J_AROW(NTERM_AROW) = II
                        IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '5', ' #2' )
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
         DO K=A_ROW_BEG,A_ROW_END                          ! 2nd, get terms from this row of A from the diagonal out
            NTERM_AROW = NTERM_AROW + 1
            J_AROW(NTERM_AROW) = J_A(K)
            IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '6', ' #1' )
         ENDDO

         B_COL_BEG = 1
j_do:    DO J=1,NCOL_B                                     ! J loops over the number of columns in B
            B_NTERM_COL_J = J_B(J+1) - J_B(J)
            IF (B_NTERM_COL_J == 0) CYCLE j_do
            B_COL_END = B_COL_BEG + B_NTERM_COL_J - 1

k_do:       DO K=1,NTERM_AROW                              ! The following 2 loops det  the ij-th term of C
               A_COL_NUM = J_AROW(K)
l_do:          DO L=B_COL_BEG,B_COL_END
                  B_ROW_NUM = I_B(L)
                  IF (A_COL_NUM == B_ROW_NUM) THEN
                     NHITS = NHITS + 1
                     NHITS_TOT_FOR_ROW_OF_A = NHITS_TOT_FOR_ROW_OF_A + 1
                     DELTA_NTERM_C = 1
                  ENDIF
               ENDDO l_do
            ENDDO k_do
            B_COL_BEG = B_COL_END + 1

            NTERM_C = NTERM_C + DELTA_NTERM_C
            IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '7', '   ' )
            DELTA_NTERM_C = 0
            NHITS = 0

         ENDDO j_do

         IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) THEN
            WRITE(F06,*)
         ENDIF

         IF (NHITS_TOT_FOR_ROW_OF_A == 0) THEN
             IF ((DEBUG(83) == 1) .OR. (DEBUG(83) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '8', '   ' )
         ENDIF
         NHITS_TOT_FOR_ROW_OF_A = 0
         A_ROW_BEG = A_ROW_END + 1

      ENDDO i_do

      IF ((DEBUG(84) == 1) .OR. (DEBUG(84) == 3)) CALL MATMULT_SSS_NTERM_DEB ( '9', '   ' )

      CALL DEALLOCATE_SPARSE_ALG ( 'J_AROW' )



      RETURN

! **********************************************************************************************************************************
  948 FORMAT(' *ERROR   948: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX B WHICH IS "',A,'" MUST BE STORED AS SYM_B = "N" FOR THIS SUBR TO WORK BUT IS STORED AS',      &
                           ' SYM_B = "',A,'"')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MATMULT_SSS_NTERM_DEB ( WHICH, ALG )

      CHARACTER( 3*BYTE), INTENT(IN)  :: ALG                    ! Which algorithm is used (#1 for terms above diag when SYM_A='Y'
!                                                                 or #2 for terms in row from diag out)
      CHARACTER( 1*BYTE)              :: WHICH                  ! Decides what to print out for this call to this subr

! **********************************************************************************************************************************
      IF      (WHICH == '1') THEN

         WRITE(F06,*)
         WRITE(F06,2011)
         WRITE(F06,2012)
         WRITE(F06,2013)
         WRITE(F06,2014) MAT_A_NAME, MAT_B_NAME, MAT_C_NAME
         WRITE(F06,2016) NROW_A, NTERM_A, NCOL_B, NTERM_B
         IF (SYM_A == 'Y') THEN
            WRITE(F06,2017)
         ELSE
            WRITE(F06,2018)
         ENDIF
         WRITE(F06,2019)
         WRITE(F06,*)

      ELSE IF (WHICH == '2') THEN

         WRITE(F06,2021)
         WRITE(F06,2022) I
         WRITE(F06,2023) I, I, A_ROW_BEG, A_ROW_END
         WRITE(F06,2024)
         WRITE(F06,*)

      ELSE IF (WHICH == '3') THEN

      ELSE IF (WHICH == '4') THEN

      ELSE IF (WHICH == '5') THEN

         WRITE(F06,2051) ALG,               K , II, J_A(K) , NTERM_AROW, J_AROW(NTERM_AROW)

      ELSE IF (WHICH == '6') THEN

         WRITE(F06,2061) ALG, K, I, J_A(K),                  NTERM_AROW, J_AROW(NTERM_AROW)

      ELSE IF (WHICH == '7') THEN

         IF (NHITS > 0) THEN
            IF (J == 1) THEN
               WRITE(F06,*)
            ENDIF
            WRITE(F06,2071) I, J, NHITS, NTERM_C
         ENDIF

      ELSE IF (WHICH == '8') THEN

         WRITE(F06,2081) I
         WRITE(F06,*)

      ELSE IF (WHICH == '9') THEN

         WRITE(F06,*)
         WRITE(F06,2021)
         WRITE(F06,2091) MAT_C_NAME, NTERM_C
         WRITE(F06,9100)
         WRITE(F06,*)

      ENDIF

! **********************************************************************************************************************************
 2011 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::START DEBUG(84) OUTPUT FROM SUBROUTINE MATMULT_SSS_NTERM::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 2012 FORMAT(' SSS SETUP FOR SPARSE MATRIX MULTIPLY ROUTINE: Determine memory required for sparse Compressed Row Storage (CRS)'   ,&
' formatted matrix C',/,' --------------------------------------------',/,' resulting from the multiplication of sparse CRS'      ,&
' matrix A and sparse Compressed Column Storage (CCS) formatted matrix B',/)

 2013 FORMAT(' A may be stored as symmetric (only terms on and above the diagonal) or with all nonzero terms included.'         ,/,&
' Matrix B must be stored with all nonzero terms.',/)

 2014 FORMAT(40X,' The name of CRS formatted matrix A is: ',A                                                                   ,/,&
             40x,' The name of CRS formatted matrix B is: ',A                                                                   ,/,&
             40x,' The name of CRS formatted matrix C is: ',A,/)

 2016 FORMAT(36X,' Matrix A has ',I8,' rows and '  ,I12,' nonzero terms'                                                        ,/,&
             36X,' Matrix B has ',I8,' cols and '  ,I12,' nonzero terms'                                                        ,/,&
             36X,' Matrix C will have                       20 nonzero terms*'                                                  ,/,&
             22X,' *(as detrmined by subr MATMULT_SSS_NTERM which had to have been run prior to this subr',/)

 2017 FORMAT(' Matrix A was flagged as being a symmetric CRS array (only nonzero terms on and above the diagonal)',/)

 2018 FORMAT(' Matrix A was flagged as a CRS array that contains all nonzero terms',/)

 2019 FORMAT(                                                                                                                      &
' In order to handle symmetric A matrices, which do not have all terms in a row (only those from the diag out), array J_AROW is',/,&
' used. J_AROW is a 1D array that contains the column numbers for all nonzero terms in one row of matrix A including those that',/,&
' may not actually be in array A due to symmetry storage. This array is used in a simulated multiply of one row of matrix A',/    ,&
' times columns of matrix B to determine whether there are any "hits" of terms in the row of A with terms in a column of B. In' ,/,&
' this manner, the number of nonzero terms (NTERM_C) that will be in array C can be estimated. The estimate is exact except for',/,&
' the case where a row of A multiplies a column of B and has "hits" that turn out to sum to zero when subr MATMULT_SSS is run' ,//,&
' Alg #1 (below) gets data for array J_AROW directly from the Compressed Row Storage (CRS) format of array A'                   ,/,&
' Alg #2 (below) is only needed if matrix A is input as symmetric (only terms on and above the diagonal) and gets terms for'    ,/,&
'         J_AROW from column I of matrix A while working on row I of matrix A. These are the terms that would be below the diag',/,&
'         in matrix A but are not explicitly in the array due to symmetry storage',//,                                             &
' For each row of matrix A, the following shows the development of array AROW and the result of multiplying AROW times matrix B',/,&
' to get one row of result matrix C. Output is only given for non null rows of matrix A and non null cols of matrix B',/)

 2021 FORMAT(' ******************************************************************************************************************',&
              '*****************')

 2022 FORMAT(30X,' W O R K I N G   O N   R O W ',I8,'   O F   O U T P U T   M A T R I X   C',/)

 2023 FORMAT(' Simulate the multiplication of row ',I8,' of matrix A times all cols of matrix B to get row ',I8,' of matrix C'  ,/,&
             ' This row of A begins in array A(K) at index K  = ',I8,' and ends at index K = ',I8,//)


 2024 FORMAT(11X,'Data from input array A           Data below diag of matrix A not in array A             Data for array AROW' ,/,&
 9X,'--------------------------          ------------------------------------------             -------------------'            ,/,&
 ' Alg     Index       Row     Col                      Index       Row     Col                         Index    Col No'           &
,/,13X,'K         I    J_A(K)                        K        II    J_A(K)                           M  J_AROW(M)')

 2051 FORMAT(1X,A3,10X,10X,10X,I25,I10,I10,I28,I11)

 2061 FORMAT(1X,A3,I10,I10,I10,25X,10X,10X,I28,I11)

 2071 FORMAT(' Row ',I8,' of A times col ',I8,' of B results in ',I8,' hits and accum number of nonzero terms in C = ',I9)

 2081 FORMAT(' There were no "hits" of any terms in row ',I8,' of matrix A multiplying terms from any col in matrix B')

 2091 FORMAT(' Total number of nonzero terms needed for allocating arrays J_C and C for matrix C = ',A,' is NTERM_C = ',I8,//)

 9100 FORMAT(' :::::::::::::::::::::::::::::::::::::END DEBUG(84) OUTPUT FROM SUBROUTINE MATMULT_SSS_NTERM:::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE MATMULT_SSS_NTERM_DEB

      END SUBROUTINE MATMULT_SSS_NTERM


      SUBROUTINE MATTRNSP_SS ( NROWA, NCOLA, NTERM, MAT_A_NAME, I_A, J_A, A, MAT_AT_NAME, I_AT, J_AT, AT )

! Transposes input matrix defined by arrays I_A, J_A, A (CRS sparse format) and puts result into output arrays
! I_AT, J_AT, AT (also in CRS sparse format)

! Input matrix A is stored in compressed row storage (CRS) format using arrays I_A(NROWA+1), J_A(NTERM_A) A(NTERM_A) where NROWA
! is the number of rows in matrix A and NTERM_A are the number of nonzero terms in matrix A:

!      I_A is an array of NROWA+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          I_A(I+1) - I_A(I) are the number of nonzero terms in row I of matrix A

!      J_A is an integer array giving the col numbers of the NTERM_A nonzero terms in matrix A

!        A is a real array of the nonzero terms in matrix A.

! Output matrix AT is stored in compressed row storage (CRS) format using arrays I_AT(NCOLA+1), J_AT(NTERM_A) AT(NTERM_A).
! The number of nonzero terms in matrix AT is the same as in matrix A:

!      I_AT is an array of NCOLA+1 integers that is used to specify the number of nonzero terms in rows of matrix AT. That is:
!          I_AT(I+1) - I_AT(I) are the number of nonzero terms in row I of matrix AT

!      J_AT is an integer array giving the col numbers of the nonzero terms in matrix AT

!        AT is a real array of the nonzero terms in matrix AT (same values as in A but arranged differently).

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATTRNSP_SS'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME        ! Name of matrix to be transposed
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_AT_NAME       ! Name of matrix that is transposed

      INTEGER(LONG), INTENT(IN)       :: NCOLA             ! Number of cols in input matrix, A
      INTEGER(LONG), INTENT(IN)       :: NROWA             ! Number of rows in input matrix, A
      INTEGER(LONG), INTENT(IN)       :: NTERM             ! Number of nonzero terms in input matrix, A
      INTEGER(LONG), INTENT(IN)       :: I_A(NROWA+1)      ! I_A(I+1) - I_A(I) are the number of nonzeros in A row I
      INTEGER(LONG), INTENT(IN)       :: J_A(NTERM)        ! Col numbers for nonzero terms in A
      INTEGER(LONG), INTENT(OUT)      :: I_AT(NCOLA+1)     ! I_AT(I+1) - I_AT(I) are the num of nonzeros in AT row I
      INTEGER(LONG), INTENT(OUT)      :: J_AT(NTERM)       ! Col numbers for nonzero terms in AT
      INTEGER(LONG)                   :: I,J               ! DO loop indices or counters
      INTEGER(LONG)                   :: ISTART            ! Starting value of I when looking for row number of a term in MATIN


      REAL(DOUBLE) , INTENT(IN)       :: A(NTERM)          ! Real nonzero values in input  matrix A
      REAL(DOUBLE) , INTENT(OUT)      :: AT(NTERM)         ! Real nonzero values in output matrix AT

      INTEGER(LONG)                   :: TMP, COL          ! temp variables for storage in loops
      INTEGER(LONG)                   :: CUMSUM            ! cumulative sum


! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NCOLA+1
         I_AT(I) = 0
      ENDDO

      DO I=1,NTERM
         J_AT(I) = 0
           AT(I) = ZERO
      ENDDO

      IF ((DEBUG(85) == 1) .OR. (DEBUG(85) == 3)) CALL MATTRNSP_SS_DEB ( '1', '   ' )

      ! count the entries per column
      TMP=0
      DO I=1,NTERM
         TMP = J_A(I)
         I_AT(TMP) = I_AT(TMP) + 1
      ENDDO

      ! cumulative sum along columns
      TMP = 0
      CUMSUM = 0
      DO I=1,NCOLA
         TMP = I_AT(I)
         I_AT(I) = CUMSUM
         CUMSUM = CUMSUM + TMP
      ENDDO
      I_AT(NCOLA+1) = NTERM

      ! do the transpose
      TMP=0
      DO I=1,NROWA
         DO J=I_A(I),I_A(I+1)-1
            COL = J_A(J)
            TMP = I_AT(COL)+1
            J_AT(TMP) = I
            AT(TMP) = A(J)
            I_AT(COL) = I_AT(COL) + 1
         ENDDO
      ENDDO

      ! right shift the array
      TMP = 0
      CUMSUM = 0
      DO I=1,NCOLA
         TMP = I_AT(I)
         I_AT(I) = CUMSUM+1
         CUMSUM = TMP
      ENDDO
      I_AT(NCOLA+1)=CUMSUM+1

      IF ((DEBUG(85) == 1) .OR. (DEBUG(85) == 3)) CALL MATTRNSP_SS_DEB ( '2', '   ' )



      RETURN

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MATTRNSP_SS_DEB ( WHICH, ALG )

      CHARACTER(LEN=*), INTENT(IN)    :: ALG                    ! Which algorithm is used (#1 for terms above diag when SYM_A='Y'
!                                                                 or #2 for terms in row from diag out)
      CHARACTER( 1*BYTE)              :: WHICH                  ! Decides what to print out for this call to this subr

      INTEGER(LONG)                   :: II                     ! Local DO loop index

! **********************************************************************************************************************************
      IF      (WHICH == '1') THEN

         WRITE(F06,*)
         WRITE(F06,1011)
         WRITE(F06,1012)
         WRITE(F06,1014) MAT_A_NAME, MAT_AT_NAME
         WRITE(F06,1016) NROWA, NCOLA, NTERM
         WRITE(F06,*)

         WRITE(F06,2800)
         WRITE(F06,*)
         WRITE(F06,3001)
         DO II=1,NROWA+1
            WRITE(F06,3002) II,I_A(II)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,3003)
         DO II=1,NTERM
            WRITE(F06,3004) II, J_A(II), A(II)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHICH == '2') THEN

         WRITE(F06,2800)
         WRITE(F06,3020) MAT_AT_NAME, NCOLA, NROWA, NTERM
         WRITE(F06,*)
         WRITE(F06,3021)
         DO II=1,NCOLA+1
            WRITE(F06,3022) II,I_AT(II)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,3023)
         DO II=1,NTERM
            WRITE(F06,3024) II, J_AT(II), AT(II)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,1095)

      ENDIF

! **********************************************************************************************************************************
 1011 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' :::::::::::::::::::::::::::::::::::::::START DEBUG(85) OUTPUT FROM SUBROUTINE MATTRNSP_SS:::::::::::::::::::::::::' &
             ,':::::::::::::::::',/)

 1012 FORMAT(' SSS SPARSE MATRIX TRANSPOSITION ROUTINE: Transpose matrix A, stored in sparse Compressed Row Storage (CRS) format,',&
' to matrix AT,',/,' ----------------------------------------',/,' also stored in sparse CRS format.',/)

 1014 FORMAT(40X,' The name of CRS formatted matrix A  is: ',A                                                                  ,/,&
             40x,' The name of CRS formatted matrix AT is: ',A,/)

 1016 FORMAT(26X,' Matrix A  has ',I8,' rows and ',I8,' cols and '  ,I12,' nonzero terms',/)

 2800 FORMAT(1X,'****************************************************************************************************************',&
                '*******************')

 3001 FORMAT(' Compressed Row Storage (CRS) format of input matrix A:'                                                          ,/,&
             ' ------------------------------------------------------'                                                         ,//,&
             ' 1) Index, L, and array I_A(L) for matrix A, where I_A(L+1) - I_A(L) is the number of nonzero terms in row L of',    &
             ' matrix A.',/,'    (also, I_A(L) is the index, K, in array A(K) where row L begins - up to, but not including, the', &
             ' last entry in I_A(L)).',/)

 3002 FORMAT('    L, I_A(L)       = ',2I12)

 3003 FORMAT(' 2) Index, K, and arrays J_A(K) and A(K). A(K) are the nonzeros in matrix A and J_A(K) is the col number in matrix A'&
            ,' for term A(K).',/)

 3004 FORMAT('    K, J_A(K), A(K) = ',2I12,1ES15.6)

 3020 FORMAT(' Output matrix AT = ',A,' has ',I8,' rows and ',I8,' cols with ',I8,' nonzero terms')

 3021 FORMAT(' Compressed Row Storage (CRS) format of output matrix AT:'                                                        ,/,&
             ' --------------------------------------------------------'                                                       ,//,&
             ' 1) Index, L, and array I_AT(L) for matrix AT, where I_AT(L+1) - I_AT(L) is the number of nonzero terms in row L of',&
             ' matrix AT.',/,'    (also, I_AT(L) is the index, K, in array AT(K) where row L begins - up to, but not including,',  &
             ' the last entry in I_AT(L)).',/)

 3022 FORMAT('    L, I_C(L)       = ',2I12)

 3023 FORMAT(' 2) Index, K, and arrays J_AT(K) and AT(K). AT(K) are the nonzeros in matrix AT and J_AT(K) is the col number in',   &
             '  matrix ATfor term AT(K).',/)

 3024 FORMAT('    K, J_AT(K), AT(K) = ',2I12,1ES15.6)

 1095 FORMAT(' ::::::::::::::::::::::::::::::::::::::END DEBUG(85) OUTPUT FROM SUBROUTINE MATTRNSP_SS::::::::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' ___________________________________________________________________________________________________________________'&
            ,'________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE MATTRNSP_SS_DEB

      END SUBROUTINE MATTRNSP_SS

   END MODULE SPARSE_MATRIX_ALGEBRA
