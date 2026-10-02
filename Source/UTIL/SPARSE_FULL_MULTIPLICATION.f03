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

   MODULE SPARSE_FULL_MULTIPLICATION

   USE SPARSE_CRS_ACCESS, ONLY :  ROW_AT_COLJ_BEGEND

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: MATMULT_SFF, MATMULT_SFS, MATMULT_SFS_NTERM

   CONTAINS

      SUBROUTINE MATMULT_SFF ( MAT_A_NAME, NROWS_A, NCOLS_A, NTERM_A, SYM_A, I_A, J_A, A, MAT_B_NAME, NROWS_B, NCOLS_B, B,         &
                               WRITE_SC1, MAT_C_NAME, CONS, C )

! Sparse matrix multiply to obtain C = cons*A*B with A in sparse format and B in full format. Input A is stored in compressed
! row storage (CRS) format. In addition, if A is symmetric it can be stored  with only the terms on, and above, the diagonal

! Input matrix A is stored in compressed row storage (CRS) format using arrays I_A(NROWS_A+1), J_A(NTERM_A) A(NTERM_A) where NROWS_A
! IS the number of rows in matrix A and NTERM_A are the number of nonzero terms in matrix A:

!      I_A is an array of NROWS_A+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          I_A(I+1) - I_A(I) are the number of nonzero terms in row I of matrix A

!      J_A is an integer array giving the col numbers of the NTERM_A nonzero terms in matrix A

!        A is a real array of the nonzero terms in matrix A. If SYM_A='Y' then only the terms on, and above, the diag are stored.

! Input  matrix B is stored in B(NROWS_B,NCOLS_B) where NROWS_B, NCOLS_B are the number of rows and columns of matrix B

! Output matrix C, which must have the same number of rows as matrix A and the same number of columns as matrix B,
! is stored in full format.

! This subr determines real array C

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE SPARSE_ALG_ARRAYS, ONLY     :  AROW, J_AROW

      USE SPARSE_ALGORITHM_LIFECYCLE, ONLY:  ALLOCATE_SPARSE_ALG, DEALLOCATE_SPARSE_ALG

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULT_SFF'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME        ! Name of matrix A
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_B_NAME        ! Name of matrix B
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_C_NAME        ! Name of matrix C
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_A             ! ='Y' if matrix A is input symmetric (terms on and above diag only)
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_SC1         ! ='Y' if we need to write to screen to advance 1 line before write msg

      INTEGER(LONG), INTENT(IN )      :: NROWS_B           ! Number of rows in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NCOLS_A           ! NUMBER OF COLS IN INPUT MATRIX A
      INTEGER(LONG), INTENT(IN )      :: NCOLS_B           ! Number of cols in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NROWS_A           ! Number of rows in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_A           ! Number of nonzero terms in input  matrix A
      INTEGER(LONG), INTENT(IN )      :: I_A(NROWS_A+1)    ! I_A(I+1) - I_A(I) = num nonzeros in row I of matrix A (CRS)
      INTEGER(LONG), INTENT(IN )      :: J_A(NTERM_A)      ! Col no's for nonzero terms in matrix A
      INTEGER(LONG)                   :: A_COL_NUM         ! A col number in matrix A
      INTEGER(LONG)                   :: A_ROW_BEG         ! Index into array I_A where a row of matrix A begins
      INTEGER(LONG)                   :: A_ROW_END         ! Index into array I_A where a row of matrix A ends
      INTEGER(LONG)                   :: A_NTERM_ROW_I     ! Number of terms in row I of matrix A
      INTEGER(LONG)                   :: AROW_MAX_TERMS    ! Max number of terms in any row of A
      INTEGER(LONG)                   :: DELTA_KTERM_C     ! Incr in KTERM_C (0,1) resulting from mult row of A x col of B
      INTEGER(LONG)                   :: I,J,K,II          ! DO loop indices
      INTEGER(LONG)                   :: KTERM_C           ! Count of number of terms put into output matrix C
      INTEGER(LONG)                   :: NCOLS_C           ! Number of cols in output matrix A (= NCOLS_B)
      INTEGER(LONG)                   :: NROWS_C           ! Number of rows in output matrix A (= NROWS_A)
      INTEGER(LONG)                   :: NHITS             ! Number of "hits" of terms in a row of A existing where terms in
!                                                            a col of B exist when a row of A is multiplied by a col of B
      INTEGER(LONG)                   :: NTERM_AROW        ! Number of nonzero terms in AROW (one row of A)


      REAL(DOUBLE) , INTENT(IN )      :: A(NTERM_A)        ! Nonzero values in matrix A
      REAL(DOUBLE) , INTENT(IN )      :: B(NROWS_B,NCOLS_B)! Real values in matrix B
      REAL(DOUBLE) , INTENT(OUT)      :: C(NROWS_A,NCOLS_B)! Real values in matrix c
      REAL(DOUBLE) , INTENT(IN )      :: CONS              ! Constant multiplier in cons*A*B to get C
      REAL(DOUBLE)                    :: EPS1              ! A small value to compare to zero



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NROWS_A
         DO J=1,NCOLS_B
            C(I,J) = ZERO
         ENDDO
      ENDDO

! Calc outputs

      IF (WRITE_SC1 /= 'Y') THEN
         !WRITE(SC1, * )                                    ! Advance 1 line for screen messages
      ENDIF

      EPS1 = EPSIL(1)

      NROWS_C = NROWS_A
      NCOLS_C = NCOLS_B

! Calc AROW_MAX_TERMS (used to allocate AROW and J_AROW)

      A_ROW_BEG      = 1
      AROW_MAX_TERMS = 0                                   ! Find the max number of nonzero terms in any row of input matrix A
      !CALL COUNTER_INIT("     setting up row ", NROWS_A)
      DO I=1,NROWS_A
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
         !CALL COUNTER_PROGRESS(I)
      ENDDO
!xx   WRITE(SC1,*) CR13
                                                           ! Allocate  vectors long enough to hold the max terms in a row of A
      CALL ALLOCATE_SPARSE_ALG ( 'J_AROW', AROW_MAX_TERMS, 0, SUBR_NAME )
      CALL ALLOCATE_SPARSE_ALG (   'AROW', AROW_MAX_TERMS, 0, SUBR_NAME )

      IF (DEBUG(82) == 1) CALL MATMULT_SFF_DEB ( '1', '   ' )

! Now do the multiply, using values put into AROW and J_AROW for each row of A. This is done to facilitate the SYM option for A

      DELTA_KTERM_C = 0                                    ! Initialize variables used in the matrix multiplication in the DO I loop
      KTERM_C       = 0
      NHITS         = 0
      A_ROW_BEG = 1
i_do: DO I=1,NROWS_A                                       ! Matrix multiply loop. Range over the rows in A

         A_NTERM_ROW_I = I_A(I+1) - I_A(I)                 ! Number of terms in matrix A in row I
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1         ! A_ROW_BEG to A_ROW_END is range of indices of terms in A for row I of A
         IF (DEBUG(82) == 1) CALL MATMULT_SFF_DEB ( '2', '   ' )

         DO K=1,AROW_MAX_TERMS                             ! Null J_AROW and AROW each time we begin a new row of A
            J_AROW(K) = 0
              AROW(K) = ZERO
         ENDDO

         NTERM_AROW = 0                                    ! Formulate J_AROW, AROW - a CRS representation of one row of A
         IF (SYM_A == 'Y') THEN                            ! 1st, look for terms that would be in this row, but are not, due to SYM
            DO K=1,A_ROW_BEG-1
              IF (J_A(K) == I) THEN
                  NTERM_AROW         = NTERM_AROW + 1
                  AROW(NTERM_AROW)   = A(K)
                  DO II=1,I-1
                     IF ((K >= I_A(II)) .AND. (K < I_A(II+1))) THEN
                        J_AROW(NTERM_AROW) = II
                        IF (DEBUG(82) == 1) CALL MATMULT_SFF_DEB ( '5', ' #2' )
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF

         DO K=A_ROW_BEG,A_ROW_END                          ! 2nd, get terms from this row of A from the diagonal out
            NTERM_AROW = NTERM_AROW + 1
            AROW(NTERM_AROW)   = A(K)
            J_AROW(NTERM_AROW) = J_A(K)
            IF (DEBUG(82) == 1) CALL MATMULT_SFF_DEB ( '6', ' #1' )
         ENDDO

         !CALL COUNTER_INIT("     calc col ", NCOLS_B)
j_do:    DO J=1,NCOLS_B                                    ! J loops over the number of columns in B

            C(I,J) = ZERO

k_do:       DO K=1,NTERM_AROW                              ! The following 2 loops produce the ij-th term of C
               A_COL_NUM = J_AROW(K)
               C(I,J) = C(I,J) + CONS*AROW(K)*B(A_COL_NUM,J)
!_do:          DO L=1,NROWS_B
!                 IF (A_COL_NUM == L) THEN
!***********         IF (DABS(B(L,J)) > EPS1) THEN !! REMOVE EPS1 TEST - IGNORING SMALL TERMS CAUSED PROBS W/ MATMULT_SFF IN LANCZOS
!                       NHITS = NHITS + 1
!                       DELTA_KTERM_C = 1
!                       C(I,J) = C(I,J) + CONS*AROW(K)*B(L,J)
!***********         ENDIF
!                 ENDIF
!              ENDDO l_do
            ENDDO k_do
            KTERM_C  = KTERM_C + DELTA_KTERM_C
            IF (NHITS > 0) THEN
               IF (DEBUG(82) == 1) CALL MATMULT_SFF_DEB ( '7', '   ' )
            ENDIF

            NHITS = 0

            !CALL COUNTER_PROGRESS(J)

         ENDDO j_do

         A_ROW_BEG = A_ROW_END + 1
         IF (DEBUG(82) == 1) THEN
            WRITE(F06,*)
         ENDIF

      ENDDO i_do
      !WRITE(SC1,*) CR13

      CALL DEALLOCATE_SPARSE_ALG ( 'J_AROW' )
      CALL DEALLOCATE_SPARSE_ALG ( 'AROW' )

      IF (DEBUG(82) == 1) CALL MATMULT_SFF_DEB ( '9', '   ' )


      RETURN

! **********************************************************************************************************************************
12345 FORMAT(5X,'calc row ',I8,' of ',I8,' col ',I8,' of ',I8,A)









! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MATMULT_SFF_DEB ( WHICH, ALG )

      CHARACTER(LEN=*), INTENT(IN)    :: ALG                    ! Which algorithm is used (#1 for terms above diag when SYM_A='Y'
!                                                                 or #2 for terms in row from diag out)
      CHARACTER( 1*BYTE)              :: WHICH                  ! Decides what to print out for this call to this subr

      INTEGER(LONG)                   :: II,JJ,KK               ! Local DO loop indices

! **********************************************************************************************************************************
      IF      (WHICH == '1') THEN

         WRITE(F06,*)
         WRITE(F06,1011)
         WRITE(F06,1012)
         WRITE(F06,1013)
         WRITE(F06,1014) MAT_A_NAME, MAT_B_NAME, MAT_C_NAME
         WRITE(F06,1015) CONS
         WRITE(F06,1016) NROWS_A, NTERM_A, NCOLS_B, NROWS_B, NCOLS_C, NROWS_C
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
            WRITE(F06,1072) I, J, NHITS, KTERM_C, I, J, C(I,J)
         ENDIF

      ELSE IF (WHICH == '8') THEN

      ELSE IF (WHICH == '9') THEN

         WRITE(F06,*)
         WRITE(F06,1091) MAT_C_NAME
         WRITE(F06,*)
         DO II=1,NROWS_C
            WRITE(F06,9192) II
            WRITE(F06,1093) (C(II,JJ),JJ=1,NCOLS_C)
            WRITE(F06,*)
         ENDDO
         WRITE(F06,1094)

      ENDIF

! **********************************************************************************************************************************
 1011 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' :::::::::::::::::::::::::::::::::::::::START DEBUG(82) OUTPUT FROM SUBROUTINE MATMULT_SFF:::::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1012 FORMAT(' SFF SPARSE MATRIX MULTIPLY ROUTINE: Multiply matrix A, stored in sparse Compressed Row Storage (CRS) format, times',&
' matrix B, stored',/,' -----------------------------------',/,' in full format, to obtain matrix C stored in full format',/)

 1013 FORMAT(' A may be stored as symmetric (only terms on and above the diagonal) or with all nonzero terms included.'         ,/,&
' Matrices B and C must be stored with all terms.',/)

 1014 FORMAT(40X,' The name of CRS  formatted matrix A is: ',A                                                                  ,/,&
             40x,' The name of full formatted matrix B is: ',A                                                                  ,/,&
             40x,' The name of full formatted matrix C is: ',A,/)

 1015 FORMAT(' Multiply ',1ES14.6,' times the product of matrix A and matrix B to obtain matrix C',/)

 1016 FORMAT(36X,' Matrix A has       ',I8,' rows and '  ,I12,' nonzero terms'                                                  ,/,&
             36X,' Matrix B has       ',I8,' cols and '  ,I12,' rows'                                                           ,/,&
             36X,' Matrix C will have ',I8,' cols and '  ,I12,' rows'/)

 1017 FORMAT(' Matrix A was input as a symmetric CRS array (only those terms on and above the diagonal are stored in array A)',/)

 1018 FORMAT(' Matrix A was input as a CRS array with all nonzero terms stored in array A',/)

 1019 FORMAT(                                                                                                                      &
' In order to handle symmetric A matrices, which do not have all terms in a row (only those from the diagonal out), arrays AROW' ,/&
' and J_AROW are used. AROW is a 1D array that contains all nonzero terms from one row of A (including those that are not'      ,/,&
' explicitly in array A due to symmetry). The multiplication of matrix A times matrix B is then accomplished (row by row of'    ,/,&
' result matrix C) by multiplying AROW times matrix B. Since AROW is a compact array containing only the nonzero terms from one',/,&
' row of matrix A, integer array J_AROW is also needed to give the col numbers, from matrix A (for one row), for the terms in'  ,/,&
' array AROW.'                                                                                                                  ,//&
' Alg #1 (below) gets data for arrays J_AROW and AROW directly from the Compressed Row Storage (CRS) format of array A'         ,/,&
' Alg #2 (below) is only needed if matrix A is input as symmetric (only terms on and above the diagonal) and gets terms for'    ,/,&
'         J_AROW and AROW from column I of matrix A while working on row I of matrix A. These are the terms that would be below',/,&
'         the diag in matrix A but are not explicitly in the array due to symmetry storage'                                    ,//,&
' For each row of matrix A, the following shows the development of arrays J_AROW and AROW and the result of mult AROW times'      ,&
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
             ' SUMMARY: Full formatted matrix C = ',A,' is:',/)

 9192 FORMAT(' Row ',I10,/,' --------------')

 1093 FORMAT(10(1ES15.6))

 1094 FORMAT(' ::::::::::::::::::::::::::::::::::::::END DEBUG(82) OUTPUT FROM SUBROUTINE MATMULT_SFF::::::::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE MATMULT_SFF_DEB

      END SUBROUTINE MATMULT_SFF


      SUBROUTINE MATMULT_SFS ( MAT_A_NAME, NROW_A, NTERM_A, SYM_A, I_A, J_A, A, MAT_B_NAME, NROW_B, NCOL_B, B,  &
                               AROW_MAX_TERMS, MAT_C_NAME, CONS,   NTERM_C, I_C, J_C, C )

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! Subroutine MATMULT_SFS_NTERM must be run before this subroutine to calculate NTERM_C, an input to this subroutine, that is the
! number of nonzero terms in C. Then memory can be allocated to C before this subroutine is called
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

! Sparse matrix multiply to obtain C = cons*A*B with A and C in sparse format and B in full format. Input A and output C are stored
! in compressed row storage (CRS) format. In addition, if A is symmetric it can be stored  with only the terms on, and above,
! the diagonal

! Input matrix A is stored in compressed row storage (CRS) format using arrays I_A(NROW_A+1), J_A(NTERM_A) A(NTERM_A) where NROW_A
! is the number of rows in matrix A and NTERM_A are the number of nonzero terms in matrix A:

!      I_A is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          I_A(I+1) - I_A(I) are the number of nonzero terms in row I of matrix A

!      J_A is an integer array giving the col numbers of the NTERM_A nonzero terms in matrix A

!        A is a real array of the nonzero terms in matrix A. If SYM_A='Y' then only the terms on, and above, the diag are stored.

! Input  matrix B is stored in B(NROW_B,NCOL_B) where NROW_B, NCOL_B are the number of rows and columns of matrix B

! Output matrix C, which must have the same number of rows as matrix A and the same number of columns as matrix B,
! is stored in compressed row storage (CRS) format using arrays I_C(NROW_A+1), J_C(NTERM_C) C(NTERM_C) where NROW_A is
! the number of rows in matrix A and NTERM_C are the number of nonzero terms in matrix C:

!      I_C is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix C. That is:
!          I_C(I+1) - I_C(I) are the number of nonzero terms in row I of matrix C

!      J_C is an integer array giving the col numbers of the NTERM_C nonzero terms in matrix C

!        C is a real array of the nonzero terms in matrix C.

! This subr determines integer arrays I_C and J_C and  real array C

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE SPARSE_CRS_ACCESS, ONLY     :  ROW_AT_COLJ_BEGEND
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULT_SFS'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME            ! Name of matrix A
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_B_NAME            ! Name of matrix B
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_C_NAME            ! Name of matrix C
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_A                 ! ='Y' if matrix A is input symmetric (terms on and above diag only)

      INTEGER(LONG), INTENT(IN )      :: NROW_B                ! Number of rows in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NCOL_B                ! Number of cols in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NROW_A                ! Number of rows in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_A               ! Number of nonzero terms in input  matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_C               ! Size of arrays J_C and C (MUST be determined by subr MATMULT_SFS)
      INTEGER(LONG), INTENT(IN )      :: AROW_MAX_TERMS        ! Max number of terms in any row of A
      INTEGER(LONG), INTENT(IN )      :: I_A(NROW_A+1)         ! I_A(I+1) - I_A(I) = num nonzeros in row I of matrix A (CRS)
      INTEGER(LONG), INTENT(IN )      :: J_A(NTERM_A)          ! Col no's for nonzero terms in matrix A
      INTEGER(LONG), INTENT(OUT)      :: I_C(NROW_A+1)         ! I_C(I+1) - I_C(I) = num nonzeros in row I of matrix C (CRS)
      INTEGER(LONG), INTENT(OUT)      :: J_C(NTERM_C)          ! Col no's for nonzero terms in matrix C
      INTEGER(LONG)                   :: A_COL_NUM             ! A col number in matrix A
      INTEGER(LONG)                   :: A_NTERM_ROW_I         ! Number of terms in row I of matrix A
      INTEGER(LONG)                   :: DELTA_KTERM_C         ! Incr in KTERM_C (0,1) resulting from mult row of A x col of B
      INTEGER(LONG)                   :: I,J,K,L,II            ! DO loop indices
      INTEGER(LONG)                   :: I1,I2                 ! DO loop range
      INTEGER(LONG)                   :: J_AROW(AROW_MAX_TERMS)! Col numbers of terms in real array AROW (see below)
      INTEGER(LONG)                   :: A_ROW_BEG             ! Index into array I_A where a row of matrix A begins
      INTEGER(LONG)                   :: A_ROW_END             ! Index into array I_A where a row of matrix A ends
      INTEGER(LONG)                   :: KTERM_C               ! Count of number of nonzero terms put into output matrix C
      INTEGER(LONG)                   :: NHITS                 ! Number of "hits" of terms in a row of A existing where terms in
!                                                                a col of B exist when a row of A is multiplied by a col of B
      INTEGER(LONG)                   :: NTERM_AROW            ! Number of nonzero terms in AROW (one row of A)
      INTEGER(LONG)                   :: ROW_AT_COLJ_BEG(NROW_A)! jth term is row number in MATIN where col j nonzeros begin
      INTEGER(LONG)                   :: ROW_AT_COLJ_END(NROW_A)! jth term is row number in MATIN where col j nonzeros end


      REAL(DOUBLE) , INTENT(IN )      :: CONS                  ! Constant multiplier in cons*A*B to get C
      REAL(DOUBLE) , INTENT(IN )      :: A(NTERM_A)            ! Nonzero values in matrix A
      REAL(DOUBLE) , INTENT(IN )      :: B(NROW_B,NCOL_B)      ! Real values in matrix B
      REAL(DOUBLE) , INTENT(OUT)      :: C(NTERM_C)            ! Nonzero values in matrix C
      REAL(DOUBLE)                    :: CTEMP                 ! A value accumulated as the nonzero terms from one row of A are
!                                                                multiplied by the corresponding nonzero terms from one col of B
      REAL(DOUBLE)                    :: AROW(AROW_MAX_TERMS)  ! Array containing the nonzero terms from one row of A
      REAL(DOUBLE)                    :: EPS1                  ! A small value to compare to zero

      INTRINSIC                       :: MAX



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NROW_A+1
         I_C(I) = 0
      ENDDO

      DO I=1,NTERM_C
         J_C(I) = 0
           C(I) = ZERO
      ENDDO

      EPS1 = EPSIL(1)

      IF ((DEBUG(83) == 2) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_DEB ( '1', '   ' )

! Create arrays that give the row numbers at which col j begins and ends. When SYM_A = 'Y'', we have to get
! terms for the A matrix that are not explicitly in a. This is done by getting A terms in the column (above the
! diagonal of a row) as well as the explicit terms from A that are there from the diagonal out to the end of the row. The two
! arrays: ROW_AT_COLJ_BEG and ROW_AT_COLJ_END are used to aid in getting the terms in the column above the diagonal.
! ROW_AT_COLJ_BEG is an array that gives, for each col of A, the starting row number of nonzero terms in that column.
! ROW_AT_COLJ_END is an array that gives, for each col of A, the ending   row number of nonzero terms in that column.
! The span: ROW_AT_COLJ_BEG to ROW_AT_COLJ_END is used when we search for terms in the columns.
! We only need ROW_AT_COLJ_BEG and ROW_AT_COLJ_END when A is input symmetric and MATOUT is not to be output as symmetric

      IF (SYM_A == 'Y') THEN

         CALL ROW_AT_COLJ_BEGEND ( MAT_A_NAME, NROW_A, NROW_A, NTERM_A, I_A, J_A, ROW_AT_COLJ_BEG, ROW_AT_COLJ_END )

      ENDIF

! Do the multiply using values put into AROW and J_AROW for each row of A. This is done to facilitate the SYM option for matrix A

      DELTA_KTERM_C = 0                                    ! Initialize variables used in the matrix multiplication in the DO I loop
      KTERM_C       = 0
      NHITS         = 0
      A_ROW_BEG     = 1
      CTEMP         = ZERO
      I_C(1)        = 1
i_do: DO I=1,NROW_A                                        ! Matrix multiply loop. Range over the rows in A

         I_C(I+1) = I_C(I)                                 ! End value in I_C for next this row is initially set at beginning value
         A_NTERM_ROW_I = I_A(I+1) - I_A(I)                 ! Number of terms in matrix A in row I
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1         ! A_ROW_BEG to _END is the range of indices of terms in A for row I of A
         IF ((DEBUG(83) == 2) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_DEB ( '2', '   ' )

         DO K=1,AROW_MAX_TERMS                             ! Null J_AROW and AROW each time we begin a new row of A
            AROW(K)   = ZERO
            J_AROW(K) = 0
         ENDDO

         NTERM_AROW = 0                                    ! Formulate J_AROW, AROW - a CRS representation of one row of A
         IF (SYM_A == 'Y') THEN                            ! 1st, look for terms that would be in this row, but are not, due to SYM
            DO K=1,A_ROW_BEG-1
               IF (J_A(K) == I) THEN
                  NTERM_AROW         = NTERM_AROW + 1
                  AROW(NTERM_AROW)   = A(K)
                  I1 = ROW_AT_COLJ_BEG(I)
                  I2 = ROW_AT_COLJ_END(I)
                  DO II=I1,I2
                     IF ((K >= I_A(II)) .AND. (K < I_A(II+1))) THEN
                        J_AROW(NTERM_AROW) = II
                        IF ((DEBUG(83) == 2) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_DEB ( '5', ' #2' )
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
         DO K=A_ROW_BEG,A_ROW_END                                ! 2nd, get terms from this row of A from the diagonal out
            NTERM_AROW = NTERM_AROW + 1
            AROW(NTERM_AROW)   = A(K)
            J_AROW(NTERM_AROW) = J_A(K)
            IF ((DEBUG(83) == 2) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_DEB ( '6', ' #1' )
         ENDDO

j_do:    DO J=1,NCOL_B                                     ! J loops over the number of columns in B

k_do:       DO K=1,NTERM_AROW                              ! The following 2 loops produce the ij-th term of C
               A_COL_NUM = J_AROW(K)
l_do:          DO L=1,NROW_B
                  IF (A_COL_NUM == L) THEN
!***********         IF (DABS(B(L,J)) > EPS1) THEN !! REMOVE EPS1 TEST - IGNORING SMALL TERMS CAUSED PROBS W/ MATMULT_SFF IN LANCZOS
                        NHITS = NHITS + 1
                        DELTA_KTERM_C = 1
                        CTEMP = CTEMP + CONS*AROW(K)*B(L,J)! This is the nonzero ij-th term in matrix C.
!***********         ENDIF
                  ENDIF
               ENDDO l_do
            ENDDO k_do

            KTERM_C  = KTERM_C + DELTA_KTERM_C             ! Now update sparse CRS representation of C
            IF (KTERM_C > NTERM_C) CALL ARRAY_SIZE_ERROR_1( SUBR_NAME, NTERM_C, MAT_C_NAME )
            IF (NHITS > 0) THEN
               I_C(I+1) = I_C(I+1) + 1
               J_C(KTERM_C) = J
                 C(KTERM_C) = CTEMP
               IF ((DEBUG(83) == 2) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_DEB ( '7', '   ' )
            ENDIF
            DELTA_KTERM_C = 0
            CTEMP = ZERO
            NHITS = 0

         ENDDO j_do
         A_ROW_BEG = A_ROW_END + 1
         IF ((DEBUG(83) == 2) .OR. (DEBUG(83) == 3)) THEN
            WRITE(F06,*)
         ENDIF

      ENDDO i_do

      IF ((DEBUG(83) == 2) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_DEB ( '9', '   ' )



      RETURN

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MATMULT_SFS_DEB ( WHICH, ALG )

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
         WRITE(F06,1016) NROW_A, NTERM_A, NCOL_B, NROW_B, NTERM_C
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
             ' :::::::::::::::::::::::::::::::::::::::START DEBUG(83) OUTPUT FROM SUBROUTINE MATMULT_SFS:::::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1012 FORMAT(' SFS SPARSE MATRIX MULTIPLY ROUTINE: Multiply matrix A, stored in sparse Compressed Row Storage (CRS) format, times',&
' matrix B, stored',/,' -----------------------------------',/,' in full format, to obtain sparse CRS matrix C',/)

 1013 FORMAT(' A may be stored as symmetric (only terms on and above the diagonal) or with all nonzero terms included.'         ,/,&
' Matrix B must be stored with all terms. Result matrix C will be stored with all nonzero terms',/)

 1014 FORMAT(40X,' The name of CRS formatted matrix A is: ',A                                                                   ,/,&
             40x,' The name of CRS formatted matrix B is: ',A                                                                   ,/,&
             40x,' The name of CRS formatted matrix C is: ',A,/)

 1015 FORMAT(' Multiply ',1ES14.6,' times the product of matrix A and matrix B to obtain matrix C',/)

 1016 FORMAT(36X,' Matrix A has ',I8,' rows and '  ,I12,' nonzero terms'                                                        ,/,&
             36X,' Matrix B has ',I8,' cols and '  ,I12,' rows'                                                        ,/,&
             36X,' Matrix C will have             ',I12,' nonzero terms*'                                                       ,/,&
             22X,'*(as detrmined by subr MATMULT_SFS_NTERM which had to have been run prior to this subr)'/)

 1017 FORMAT(' Matrix A was input as a symmetric CRS array (only those terms on and above the diagonal are stored in array A)',/)

 1018 FORMAT(' Matrix A was input as a CRS array with all nonzero terms stored in array A',/)

 1019 FORMAT(                                                                                                                      &
' In order to handle symmetric A matrices, which do not have all terms in a row (only those from the diagonal out), arrays AROW' ,/&
' and J_AROW are used. AROW is a 1D array that contains all nonzero terms from one row of A (including those that are not'      ,/,&
' explicitly in array A due to symmetry). The multiplication of matrix A times matrix B is then accomplished (row by row of'    ,/,&
' result matrix C) by multiplying AROW times matrix B. Since AROW is a compact array containing only the nonzero terms from one',/,&
' row of matrix A, integer array J_AROW is also needed to give the col numbers, from matrix A (for one row), for the terms in'  ,/,&
' array AROW.'                                                                                                                  ,//&
' Alg #1 (below) gets data for arrays J_AROW and AROW directly from the Compressed Row Storage (CRS) format of array A'         ,/,&
' Alg #2 (below) is only needed if matrix A is input as symmetric (only terms on and above the diagonal) and gets terms for'    ,/,&
'         J_AROW and AROW from column I of matrix A while working on row I of matrix A. These are the terms that would be below',/,&
'         the diag in matrix A but are not explicitly in the array due to symmetry storage'                                    ,//,&
' For each row of matrix A, the following shows the development of arrays J_AROW and AROW and the result of mult AROW times'      ,&
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
             ' 1) Index, L, and array I_C(L) for matrix C, where I_C(L+1) - I_C(L) is the number of nonzero terms in row L of'    ,&
             ' matrix C.',/,'    (also, I_C(L) is the index, K, in array C(K) where row L begins - up to, but not including, the' ,&
             ' last entry in I_C(L)).',/)

 9192 FORMAT('    L, I_C(L)       = ',2I12)

 1093 FORMAT(' 2) Index, K, and arrays J_C(K) and C(K). C(K) are the nonzeros in matrix C and J_C(K) is the col number in matrix', &
             ' C for term C(K).',/)

 1094 FORMAT('    K, J_C(K), C(K) = ',2I12,1ES15.6)

 1095 FORMAT(' ::::::::::::::::::::::::::::::::::::::END DEBUG(83) OUTPUT FROM SUBROUTINE MATMULT_SFS::::::::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE MATMULT_SFS_DEB

      END SUBROUTINE MATMULT_SFS


      SUBROUTINE MATMULT_SFS_NTERM ( MAT_A_NAME, NROW_A, NTERM_A, SYM_A, I_A, J_A, MAT_B_NAME, NROW_B, NCOL_B, B,  &
                                     AROW_MAX_TERMS, MAT_C_NAME, NTERM_C )

! Setup routine for subr MATMULT_SFS (which performs the matrix operation C = cons*a*b with A and C in sparse format and B full
! Matrices A and C are stored in compressed row storage (CRS) format. In addition, if A is symmetric it can be stored  with only the
! terms on, and above, the diagonal.

! Matrix A, for example, is stored using arrays I_A(NROW_A+1),
! J_A(NTERM_A) and A(NTERM_A) where NROW_A is the number of rows in matrix A and NTERM_A are the number of nonzero terms in matrix A
! (similar definition for B and C):

!      I_A is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          I_A(I+1) - I_A(I) are the number of nonzero terms in row I of matrix A

!      J_A is an integer array giving the col numbers of the NTERM_A nonzero terms in matrix A

!        A is a real array of the nonzero terms in matrix A. If SYM_A='Y' then only the terms on, and above, the diag are stored.

! Input  matrix B is stored in B(NROW_B,NCOL_B) where NROW_B, NCOL_B are the number of rows and columns of matrix B

! This subr determines the storage required (NTERM_C) for array J_C from the data in arrays I_A, J_A, and B. This information is
! used to allocate memory for arrays J_C and C prior to calling subr MATMULT_SFS so that it can it can do the sparse matrix multiply

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  EPSIL
      USE SPARSE_ALG_ARRAYS, ONLY     :  J_AROW
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE SPARSE_ALGORITHM_LIFECYCLE, ONLY:  ALLOCATE_SPARSE_ALG, DEALLOCATE_SPARSE_ALG
      USE SPARSE_CRS_ACCESS, ONLY     :  ROW_AT_COLJ_BEGEND

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULT_SFS_NTERM'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME             ! Name of matrix A
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_B_NAME             ! Name of matrix B
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_C_NAME             ! Name of matrix C
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_A                  ! ='Y' if matrix A is input symmetric (terms on and above diag only)

      INTEGER(LONG), INTENT(IN )      :: NROW_B                 ! Number of rows in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NCOL_B                 ! Number of cols in input matrix B
      INTEGER(LONG), INTENT(IN )      :: NROW_A                 ! Number of rows in input matrix A
      INTEGER(LONG), INTENT(IN )      :: NTERM_A                ! Number of nonzero terms in input matrix A
      INTEGER(LONG), INTENT(IN )      :: I_A(NROW_A+1)          ! I_A(I+1) - I_A(I) = no. terms in row I of matrix A
      INTEGER(LONG), INTENT(IN )      :: J_A(NTERM_A)           ! Col no's for nonzero terms in matrix A
      INTEGER(LONG), INTENT(OUT)      :: AROW_MAX_TERMS         ! Max number of terms in any row of A
      INTEGER(LONG), INTENT(OUT)      :: NTERM_C                ! Number of nonzero terms in output matrix C
      INTEGER(LONG)                   :: A_COL_NUM              ! A col number in matrix A
      INTEGER(LONG)                   :: A_NTERM_ROW_I          ! Number of terms in row I of matrix A
      INTEGER(LONG)                   :: A_ROW_BEG              ! Index into array I_A where a row of matrix A begins
      INTEGER(LONG)                   :: A_ROW_END              ! Index into array I_A where a row of matrix A ends
      INTEGER(LONG)                   :: DELTA_NTERM_C = 0      ! Incr in NTERM_C (0,1) resulting from mult row of A x col of B
      INTEGER(LONG)                   :: I,J,K,L,II             ! DO loop indices
      INTEGER(LONG)                   :: I1,I2                  ! DO loop range
      INTEGER(LONG)                   :: NHITS                  ! Number of "hits" of terms in a row of A existing where terms in
!                                                            a col of B exist when a row of A is multiplied by a col of B
      INTEGER(LONG)                   :: NHITS_TOT_FOR_ROW_OF_A ! Num of "hits" of terms in a col of A existing where terms in any
!                                                                 row of B exist when a row of A is multiplied by al cols of B
      INTEGER(LONG)                   :: NTERM_AROW             ! Number of nonzero terms in AROW (one row of A)
      INTEGER(LONG)                   :: ROW_AT_COLJ_BEG(NROW_A)! jth term is row number in MATIN where col j nonzeros begin
      INTEGER(LONG)                   :: ROW_AT_COLJ_END(NROW_A)! jth term is row number in MATIN where col j nonzeros end


      REAL(DOUBLE) , INTENT(IN )      :: B(NROW_B,NCOL_B)       ! Real values in matrix B


      INTRINSIC                       :: MAX



! **********************************************************************************************************************************
! Initialize outputs

      AROW_MAX_TERMS = 0
      NTERM_C        = 0

! First, set up array J_AROW so it can handle the row from A that has the most terms in it

      A_ROW_BEG      = 1
      AROW_MAX_TERMS = 0                                   ! Find the max number of nonzero terms in any row of input matrix A
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
      CALL ALLOCATE_SPARSE_ALG ( 'J_AROW', AROW_MAX_TERMS, 0, SUBR_NAME )

      IF ((DEBUG(83) == 1) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_NTERM_DEB ( '1', '   ' )

! Create arrays that give the row numbers at which col j begins and ends. When SYM_A = 'Y'', we have to get
! terms for the A matrix that are not explicitly in a. This is done by getting A terms in the column (above the
! diagonal of a row) as well as the explicit terms from A that are there from the diagonal out to the end of the row. The two
! arrays: ROW_AT_COLJ_BEG and ROW_AT_COLJ_END are used to aid in getting the terms in the column above the diagonal.
! ROW_AT_COLJ_BEG is an array that gives, for each col of A, the starting row number of nonzero terms in that column.
! ROW_AT_COLJ_END is an array that gives, for each col of A, the ending   row number of nonzero terms in that column.
! The span: ROW_AT_COLJ_BEG to ROW_AT_COLJ_END is used when we search for terms in the columns.
! We only need ROW_AT_COLJ_BEG and ROW_AT_COLJ_END when A is input symmetric and MATOUT is not to be output as symmetric

      IF (SYM_A == 'Y') THEN

         CALL ROW_AT_COLJ_BEGEND ( MAT_A_NAME, NROW_A, NROW_A, NTERM_A, I_A, J_A, ROW_AT_COLJ_BEG, ROW_AT_COLJ_END )

      ENDIF

! Now count the terms that go into C, using values put into AROW and J_AROW for each row of A. This is done to facilitate the
! SYM option for matrix A.

      NHITS   = 0
      NHITS_TOT_FOR_ROW_OF_A = 0
      NTERM_C = 0
      A_ROW_BEG  = 1
i_do: DO I=1,NROW_A

         A_NTERM_ROW_I = I_A(I+1) - I_A(I)
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1
         IF ((DEBUG(83) == 1) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_NTERM_DEB ( '2', '   ' )

         DO K=1,AROW_MAX_TERMS                             ! Null J_AROW and AROW each time we begin a new row of A
            J_AROW(K) = 0
         ENDDO

         NTERM_AROW = 0                                    ! Formulate J_AROW, AROW - a CRS representation of one row of A
         IF (SYM_A == 'Y') THEN                            ! 1st, look for terms that would be in this row, but are not, due to SYM
            DO K=1,A_ROW_BEG-1
               IF (J_A(K) == I) THEN
                  NTERM_AROW         = NTERM_AROW + 1
                  I1 = ROW_AT_COLJ_BEG(I)
                  I2 = ROW_AT_COLJ_END(I)
                  DO II=I1,I2
                     IF ((K >= I_A(II)) .AND. (K < I_A(II+1))) THEN
                        J_AROW(NTERM_AROW) = II
                        IF ((DEBUG(83) == 1) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_NTERM_DEB ( '5', '   ' )
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
         DO K=A_ROW_BEG,A_ROW_END                          ! 2nd, get terms from this row of A from the diagonal out
            NTERM_AROW = NTERM_AROW + 1
            J_AROW(NTERM_AROW) = J_A(K)
            IF ((DEBUG(83) == 1) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_NTERM_DEB ( '6', '   ' )

         ENDDO

j_do:    DO J=1,NCOL_B                                     ! J loops over the number of columns in B

k_do:       DO K=1,NTERM_AROW                              ! The following 2 loops produce the ij-th term of C
               A_COL_NUM = J_AROW(K)
l_do:          DO L=1,NROW_B
                  IF (A_COL_NUM == L) THEN
!******************* IF (DABS(B(L,J)) > 0.D0) THEN!! REMOVE EPS1 TEST - IGNORING SMALL TERMS CAUSED PROBS W/ MATMULT_SFF IN LANCZOS
                        NHITS = NHITS + 1
                        NHITS_TOT_FOR_ROW_OF_A = NHITS_TOT_FOR_ROW_OF_A + 1
                        DELTA_NTERM_C = 1
!******************* ENDIF
                  ENDIF
               ENDDO l_do
            ENDDO k_do

            NTERM_C = NTERM_C + DELTA_NTERM_C
            IF ((DEBUG(83) == 1) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_NTERM_DEB ( '7', '   ' )
            DELTA_NTERM_C = 0
            NHITS = 0

         ENDDO j_do

         IF ((DEBUG(83) == 1) .OR. (DEBUG(83) == 3)) THEN
            WRITE(F06,*)
         ENDIF

         IF (NHITS_TOT_FOR_ROW_OF_A == 0) THEN
             IF ((DEBUG(83) == 1) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_NTERM_DEB ( '8', '   ' )
         ENDIF
         NHITS_TOT_FOR_ROW_OF_A = 0
         A_ROW_BEG = A_ROW_END + 1

      ENDDO i_do

      IF ((DEBUG(83) == 1) .OR. (DEBUG(83) == 3)) CALL MATMULT_SFS_NTERM_DEB ( '9', '   ' )

      CALL DEALLOCATE_SPARSE_ALG ( 'J_AROW' )



      RETURN

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MATMULT_SFS_NTERM_DEB ( WHICH, ALG )

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
         WRITE(F06,2016) NROW_A, NTERM_A, NCOL_B, NROW_B
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
             ' ::::::::::::::::::::::::::::::::::::START DEBUG(83) OUTPUT FROM SUBROUTINE MATMULT_SFS_NTERM::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 2012 FORMAT(' SFS SETUP FOR SPARSE MATRIX MULTIPLY ROUTINE: Determine memory required for sparse Compressed Row Storage (CRS)'   ,&
' formatted matrix C',/,' ---------------------------------------------',/,' resulting from the multiplication of sparse CRS'     ,&
' formatted matrix A and full formatted matrix B')

 2013 FORMAT(' A may be stored as symmetric (only terms on and above the diagonal) or with all nonzero terms included.'         ,/,&
' Matrix B must be stored with all terms.',/)

 2014 FORMAT(40X,' The name of CRS  formatted matrix A is: ',A                                                                  ,/,&
             40x,' The name of full formatted matrix B is: ',A                                                                  ,/,&
             40x,' The name of CRS  formatted matrix C is: ',A,/)

 2016 FORMAT(36X,' Matrix A has ',I8,' rows and '  ,I12,' nonzero terms'                                                        ,/,&
             36X,' Matrix B has ',I8,' cols and '  ,I12,' rows'                                                                 ,/,&
             36X,' Matrix C will have                       20 nonzero terms*'                                                  ,/,&
             22X,' *(as detrmined by subr MATMULT_SFS_NTERM which had to have been run prior to this subr',/)

 2017 FORMAT(' Matrix A was flagged as being a symmetric CRS array (only nonzero terms on and above the diagonal)',/)

 2018 FORMAT(' Matrix A was flagged as a CRS array that contains all nonzero terms',/)

 2019 FORMAT(                                                                                                                      &
' In order to handle symmetric A matrices, which do not have all terms in a row (only those from the diag out), array J_AROW is',/,&
' used. J_AROW is a 1D array that contains the column numbers for all nonzero terms in one row of matrix A including those that',/,&
' may not actually be in array A due to symmetry storage. This array is used in a simulated multiply of one row of matrix A',/    ,&
' times columns of matrix B to determine whether there are any "hits" of terms in the row of A with terms in a column of B. In' ,/,&
' this manner, the number of nonzero terms (NTERM_C) that will be in array C can be estimated. The estimate is exact except for',/,&
' the case where a row of A multiplies a column of B and has "hits" that turn out to sum to zero when subr MATMULT_SFS is run' ,//,&
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

 9100 FORMAT(' :::::::::::::::::::::::::::::::::::::END DEBUG(83) OUTPUT FROM SUBROUTINE MATMULT_SFS_NTERM:::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE MATMULT_SFS_NTERM_DEB

      END SUBROUTINE MATMULT_SFS_NTERM

   END MODULE SPARSE_FULL_MULTIPLICATION
