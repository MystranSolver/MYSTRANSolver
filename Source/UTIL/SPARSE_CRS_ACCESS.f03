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

   MODULE SPARSE_CRS_ACCESS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: GET_SPARSE_CRS_COL, GET_SPARSE_CRS_ROW, ROW_AT_COLJ_BEGEND, SPARSE_CRS_TERM_COUNT, SPARSE_MAT_DIAG_ZEROS

   CONTAINS

      SUBROUTINE GET_SPARSE_CRS_COL ( MATIN_NAME, COL_NUM, NTERM, NROWS, NCOLS, I_MATIN, J_MATIN, MATIN, BETA, OUT_VEC, NULL_COL )

! Gets col number COL_NUM from a matrix in sparse (compressed row storage) format described by I_MATIN, J_MATIN, MATIN
! arrays, multiplies it by BETA, and puts result into array OUT_VEC. Sets NULL_COL to 'Y' if result is null.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_SPARSE_CRS_COL'
      CHARACTER(LEN=*), INTENT(IN )   :: MATIN_NAME        ! Name of input matrix to be partitioned
      CHARACTER(1*BYTE),INTENT(OUT)   :: NULL_COL          ! = 'Y' if OUT_VEC is null

      INTEGER(LONG), INTENT(IN )      :: NROWS             ! No. rows in MATIN
      INTEGER(LONG), INTENT(IN )      :: NTERM             ! No. terms in MATIN
      INTEGER(LONG), INTENT(IN )      :: I_MATIN(NROWS+1)  ! Starting locations in MATIN for each row
      INTEGER(LONG), INTENT(IN )      :: J_MATIN(NTERM)    ! Col numbers for terms in MATIN
      INTEGER(LONG), INTENT(IN )      :: NCOLS             ! No. cols in MATIN
      INTEGER(LONG), INTENT(IN )      :: COL_NUM           ! Col number for the col to get in MATIN
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: NUM_TERMS_IN_ROW  ! No. terms in a row of MATIN. Each term will be checked to see if it
!                                                            belongs to col number COL_NUM


      REAL(DOUBLE) , INTENT(IN)       :: MATIN(NTERM)      ! Nonzero terms in matrix MATIN
      REAL(DOUBLE) , INTENT(IN)       :: BETA              ! Scalar multiplier for row from MATIN
      REAL(DOUBLE) , INTENT(OUT)      :: OUT_VEC(NROWS)    ! Output vector containing the terms from col COL_NUM of MATIN



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NROWS
         OUT_VEC(I) = ZERO
      ENDDO

! Make sure COL_NUM is a col of MATIN

      IF ((COL_NUM < 1) .OR. (COL_NUM > NCOLS)) THEN
         WRITE(ERR,930) SUBR_NAME,COL_NUM,MATIN_NAME,NROWS
         WRITE(F06,930) SUBR_NAME,COL_NUM,MATIN_NAME,NROWS
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error (invalid COL_NUM), so quit
      ENDIF

! Null OUT_VEC

      DO J=1,NROWS
         OUT_VEC(J) = ZERO
      ENDDO

! Load nonzero terms from col COL_NUM of MATIN into OUT_VEC

      NULL_COL = 'Y'
      K = 0
      DO I=1,NROWS
         NUM_TERMS_IN_ROW = I_MATIN(I+1) - I_MATIN(I)
         DO J=1,NUM_TERMS_IN_ROW                           ! Check each term to see if it is in col number COL_NUM
            K = K + 1
            IF (J_MATIN(K) == COL_NUM) THEN
               NULL_COL = 'N'
               OUT_VEC(I) = BETA*MATIN(K)
            ENDIF
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************
  930 FORMAT(' *ERROR   930: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ATTEMPT TO GET COL = ',I12,' FROM MATRIX ',A,' WHEN IT ONLY HAS COL NUMBERS 1 THRU ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE GET_SPARSE_CRS_COL


      SUBROUTINE GET_SPARSE_CRS_ROW ( MATIN_NAME, ROW_NUM, NTERM, NROWS, NCOLS, I_MATIN, J_MATIN, MATIN, BETA, OUT_VEC, NULL_ROW )

! Gets a row (ROW_NUM) from a matrix in sparse (compressed row storage) format described by I_MATIN, J_MATIN, MATIN
! arrays, multiplies it by BETA, and puts result into array OUT_VEC. Sets NULL_ROW to 'Y' if result is null.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_SPARSE_CRS_ROW'
      CHARACTER(LEN=*), INTENT(IN )   :: MATIN_NAME        ! Name of input matrix to be partitioned
      CHARACTER(LEN=*), INTENT(OUT)   :: NULL_ROW          ! = 'Y' if OUT_VEC is null

      INTEGER(LONG), INTENT(IN )      :: NCOLS             ! No. cols in MATIN
      INTEGER(LONG), INTENT(IN )      :: NROWS             ! No. rows in MATIN
      INTEGER(LONG), INTENT(IN )      :: NTERM             ! No. terms in MATIN
      INTEGER(LONG), INTENT(IN )      :: I_MATIN(NROWS+1)  ! Starting locations in MATIN for each row
      INTEGER(LONG), INTENT(IN )      :: J_MATIN(NTERM)    ! Col numbers for terms in MATIN
      INTEGER(LONG), INTENT(IN )      :: ROW_NUM           ! Row number for the row to get in MATIN
      INTEGER(LONG)                   :: J,K               ! DO loop indices or counters
      INTEGER(LONG)                   :: NUM_TERMS_IN_ROW  ! No. terms in row ROW_NUM of MATIN


      REAL(DOUBLE) , INTENT(IN)       :: MATIN(NTERM)      ! Nonzero terms in matrix MATIN
      REAL(DOUBLE) , INTENT(IN)       :: BETA              ! Scalar multiplier for row from MATIN
      REAL(DOUBLE) , INTENT(OUT)      :: OUT_VEC(NCOLS)    ! Output vector containing the terms from row ROW_NUM of MATIN



! **********************************************************************************************************************************
! Initialize outputs

      DO J=1,NCOLS
         OUT_VEC(J) = ZERO
      ENDDO

! Make sure ROW NUM is a row of MATIN

      IF ((ROW_NUM < 1) .OR. (ROW_NUM > NROWS)) THEN
         WRITE(ERR,929) SUBR_NAME,ROW_NUM,MATIN_NAME,NROWS
         WRITE(F06,929) SUBR_NAME,ROW_NUM,MATIN_NAME,NROWS
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                              ! Coding error (invalid ROW_NUM), so quit
      ENDIF

! Load nonzero terms from row ROW_NUM of MATIN into OUT_VEC

      NULL_ROW = 'N'
      NUM_TERMS_IN_ROW = I_MATIN(ROW_NUM+1) - I_MATIN(ROW_NUM)
      IF (NUM_TERMS_IN_ROW == 0) THEN
         NULL_ROW = 'Y'
      ELSE
         K = I_MATIN(ROW_NUM)
         DO J=1,NUM_TERMS_IN_ROW
            OUT_VEC(J_MATIN(K)) = BETA*MATIN(K)
            K = K + 1
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
  929 FORMAT(' *ERROR   929: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ATTEMPT TO GET ROW = ',I12,' FROM MATRIX ',A,' WHEN IT ONLY HAS ROW NUMBERS 1 THRU ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE GET_SPARSE_CRS_ROW


      SUBROUTINE GET_SPARSE_MAT_TERM ( MATIN_NAME, I_MATIN, J_MATIN, MATIN, IROW, JCOL, N, NTERMS, MATIN_VAL )

! Given a row/col index, gets the real value from a sparse CRS matrix

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_SPARSE_MAT_TERM'
      CHARACTER(LEN=*), INTENT(IN )   :: MATIN_NAME        ! Name of input matrix MATIN

      INTEGER(LONG), INTENT(IN)       :: N                 ! Row/col size of MATIN
      INTEGER(LONG), INTENT(IN)       :: NTERMS            ! Number of nonzero terms in sparse matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: IROW              ! Row index of the term to retrieve from sparse MATIN
      INTEGER(LONG), INTENT(IN)       :: JCOL              ! Col index of the term to retrieve from sparse MATIN
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(N+1)      ! Indices of the beginning terms in each row for MATIN values
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERMS)   ! Col numbers of nonzero term in MATIN
      INTEGER(LONG)                   :: NUM_TERMS_IN_ROW  ! No. terms in row IROW of MATIN
      INTEGER(LONG)                   :: J                 ! DO loop index
      INTEGER(LONG)                   :: K                 ! Counter


      REAL(DOUBLE) , INTENT(IN)       :: MATIN(NTERMS)     ! Real vals in sparse matrix MATIN
      REAL(DOUBLE) , INTENT(OUT)      :: MATIN_VAL



! **********************************************************************************************************************************
! Initialize output value

      MATIN_VAL = ZERO

! Make sure IROW, JCOL are valid row/col numbers of MATIN

      IF ((IROW < 1) .OR. (IROW > N)) THEN
         WRITE(ERR,929) SUBR_NAME, 'ROW', IROW, MATIN_NAME, 'ROW', N
         WRITE(F06,929) SUBR_NAME, IROW, MATIN_NAME, N
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF ((JCOL < 1) .OR. (JCOL > N)) THEN
         WRITE(ERR,929) SUBR_NAME, 'COL', JCOL, MATIN_NAME, 'COL', N
         WRITE(F06,929) SUBR_NAME, JCOL, MATIN_NAME, N
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      NUM_TERMS_IN_ROW = I_MATIN(IROW+1) - I_MATIN(IROW)
      K = I_MATIN(IROW)
      DO J=1,NUM_TERMS_IN_ROW
         IF (K > NTERMS) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, NTERMS, MATIN_NAME )
         IF (J_MATIN(K) == JCOL) THEN
            MATIN_VAL = MATIN(K)
         ENDIF
         K = K + 1
      ENDDO



      RETURN

! **********************************************************************************************************************************
  929 FORMAT(' *ERROR   929: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ATTEMPT TO GET ',A,' = ',I12,' FROM MATRIX ',A,' WHEN IT ONLY HAS ',A,'  NUMBERS 1 THRU ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE GET_SPARSE_MAT_TERM



      SUBROUTINE ROW_AT_COLJ_BEGEND ( NAME, NROWS, NCOLS, NTERM, I_A, J_A, ROW_AT_COLJ_BEG, ROW_AT_COLJ_END )

! Creates arrays ROW_AT_COLJ_BEG and ROW_AT_COLJ_END which are:

! ROW_AT_COLJ_BEG is an array that gives, for each col of sparse CRS input matrix A, the start row no. of nonzero terms in that col.
! ROW_AT_COLJ_END is an array that gives, for each col of sparse CRS input matrix A, the last  row no. of nonzero terms in that col.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ROW_AT_COLJ_BEGEND'
      CHARACTER(LEN=*), INTENT(IN )   :: NAME                  ! Name of input matrix A

      INTEGER(LONG), INTENT(IN )      :: NTERM                 ! No. terms in MATIN
      INTEGER(LONG), INTENT(IN )      :: NROWS                 ! No. rows in MATIN
      INTEGER(LONG), INTENT(IN )      :: NCOLS                 ! No. cols in MATIN
      INTEGER(LONG), INTENT(IN )      :: I_A(NROWS+1)          ! I_A(i+1) - I_A(i) is no. terms in row i of matrix A
      INTEGER(LONG), INTENT(IN )      :: J_A(NTERM)            ! Array of column numbers for matrix A
      INTEGER(LONG), INTENT(OUT)      :: ROW_AT_COLJ_BEG(NCOLS)! jth term is row number in MATIN where col j nonzeros begin
      INTEGER(LONG), INTENT(OUT)      :: ROW_AT_COLJ_END(NCOLS)! jth term is row number in MATIN where col j nonzeros end
      INTEGER(LONG)                   :: COL_NUM               ! A column number from J_MATIN
      INTEGER(LONG)                   :: I,J,K                 ! DO loop indices or counters
      INTEGER(LONG)                   :: NTERM_ROW_I           ! Number of terms in matrix A row I




! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NCOLS
         ROW_AT_COLJ_END(I) = 0
      ENDDO


      DO J=1,NCOLS                                         ! Initialize arrays
         ROW_AT_COLJ_BEG(J) = 0                            ! This is the array that has the row nos at which each col of NAME begins
         ROW_AT_COLJ_END(J) = 0                            ! This is the array that has the row nos at which each col of NAME ends
      ENDDO

      K = 0
      DO I=1,NROWS
         NTERM_ROW_I = I_A(I+1) - I_A(I)
         DO J=1,NTERM_ROW_I
            K = K + 1
            IF (K > NTERM) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, NTERM, NAME)
            COL_NUM = J_A(K)
            IF (ROW_AT_COLJ_BEG(COL_NUM) == 0) THEN
               ROW_AT_COLJ_BEG(COL_NUM) = I
            ENDIF
         ENDDO
      ENDDO

      K = NTERM + 1
      DO I=NROWS,1,-1
         NTERM_ROW_I = I_A(I+1) - I_A(I)
         DO J=NTERM_ROW_I,1,-1
            K = K - 1
            IF (K < 1) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, NTERM, NAME)
            COL_NUM = J_A(K)
            IF (ROW_AT_COLJ_END(COL_NUM) == 0) THEN
               ROW_AT_COLJ_END(COL_NUM) = I
            ENDIF
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE ROW_AT_COLJ_BEGEND



      SUBROUTINE SPARSE_CRS_TERM_COUNT ( NROWS, NTERM_IN, MATIN_NAME, I_MATIN, J_MATIN, NTERM_OUT )

! Counts terms in sparse compressed row (CRS) storage format matrices that are on, or above, the diagonal. Used to get the number
! of terms from a matrix stored as sparse nonsym that will be in the same matrix stored as sparse sym.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_CRS_TERM_COUNT'
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME          ! Name of input matrix

      INTEGER(LONG), INTENT(IN)       :: NROWS               ! Number of rows in input matrix, MATIN
      INTEGER(LONG), INTENT(IN)       :: NTERM_IN            ! Number of nonzero terms in input matrix, MATIN
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(NROWS+1)    ! I_MATIN(I+1) - I_MATIN(I) are the number of nonzeros in MATIN row I
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERM_IN)   ! Col numbers for nonzero terms in MATIN
      INTEGER(LONG), INTENT(OUT)      :: NTERM_OUT           ! Number of nonzero terms in output matrix, MATOUT
      INTEGER(LONG)                   :: I,K                 ! DO loop indices or counters




! **********************************************************************************************************************************
! Initialize outputs

      NTERM_OUT = 0

! Calc outputs

      NTERM_OUT   = 0
      DO I=1,NROWS
         DO K=I_MATIN(I),I_MATIN(I+1)-1
            IF (J_MATIN(K) >= I) THEN                      ! This is a term on or above the diag that will go into the output matrix
               NTERM_OUT = NTERM_OUT + 1
            ELSE
               CYCLE
            ENDIF
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE SPARSE_CRS_TERM_COUNT


       SUBROUTINE SPARSE_MAT_DIAG_ZEROS ( NAME, NROWS_A, NTERM_A, I_A, J_A, NUM_A_DIAG_ZEROS )

! Determines the number of zero diagonal terms in an input matrix that is stored in compressed row storage format (CRS format)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  WRT_ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_MAT_DIAG_ZEROS'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME               ! Name of input matrix
      CHARACTER( 1*BYTE)              :: FND_DIAG_TERM='x'  ! If 'Y' we found a diag term in a row of matrix NAME

      INTEGER(LONG), INTENT(IN)       :: NROWS_A            ! Number of rows in input matrix A
      INTEGER(LONG), INTENT(IN)       :: NTERM_A            ! Number of nonzero terms in input matrix A
      INTEGER(LONG), INTENT(IN)       :: I_A(NROWS_A+1)     ! Array of row no's for terms in input matrix A
      INTEGER(LONG), INTENT(IN)       :: J_A(NTERM_A)       ! Array of col no's for terms in input matrix A
      INTEGER(LONG)                   :: A_NTERM_ROW_I      ! Number of terms in row I of matrix A
      INTEGER(LONG)                   :: A_ROW_BEG          ! Index into array I_A where a row of matrix A begins
      INTEGER(LONG)                   :: A_ROW_END          ! Index into array I_A where a row of matrix A ends
      INTEGER(LONG), INTENT(OUT)      :: NUM_A_DIAG_ZEROS   ! Number of zero diagonal terms in input matrix A
      INTEGER(LONG)                   :: I,K                ! DO loop indices
      INTEGER(LONG)                   :: ZERO_DIAGS(NROWS_A)! Row numbers where there are zero diag terms




! **********************************************************************************************************************************
! Initialize outputs

      NUM_A_DIAG_ZEROS = 0

! Calc outputs

      DO I=1,NROWS_A
         ZERO_DIAGS(I) = 0
      ENDDO

      A_ROW_BEG        = 1
i_do: DO I=1,NROWS_A
         A_NTERM_ROW_I = I_A(I+1) - I_A(I)
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1         ! A_ROW_BEG to A_ROW_END is range of indices of terms in A for row I of A
         IF (A_NTERM_ROW_I == 0) THEN                      ! If there are no terms in row I then there is a zero diag term
            NUM_A_DIAG_ZEROS = NUM_A_DIAG_ZEROS + 1
            ZERO_DIAGS(NUM_A_DIAG_ZEROS) = I
            FND_DIAG_TERM = 'Y'
         ELSE
k_do:       DO K=A_ROW_BEG,A_ROW_END
               FND_DIAG_TERM = 'Y'
               IF (J_A(K) == I) THEN                       ! There is a diag term in this row, so cycle on rows
                  A_ROW_BEG = A_ROW_END + 1
                  CYCLE i_do
               ELSE
                  FND_DIAG_TERM = 'N'
               ENDIF
            ENDDO k_do
            IF (FND_DIAG_TERM == 'N') THEN
               NUM_A_DIAG_ZEROS = NUM_A_DIAG_ZEROS + 1     ! In a row with some nonzero terms we have found no diag term
               ZERO_DIAGS(NUM_A_DIAG_ZEROS) = I
            ENDIF
         ENDIF
         A_ROW_BEG = A_ROW_END + 1
      ENDDO i_do

      IF (DEBUG(89) > 0) THEN
         WRITE(F06,100) NAME, NUM_A_DIAG_ZEROS
         WRITE(F06,101) (ZERO_DIAGS(I),I=1,NUM_A_DIAG_ZEROS)
      ENDIF



      RETURN

! **********************************************************************************************************************************
  100 FORMAT(' *INFORMATION: MATRIX ',A,' HAS ',I8,' ZERO DIAGONAL TERMS IN ROWS:')

  101 FORMAT(16I8)

! **********************************************************************************************************************************

      END SUBROUTINE SPARSE_MAT_DIAG_ZEROS

   END MODULE SPARSE_CRS_ACCESS
