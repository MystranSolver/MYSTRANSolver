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

   MODULE SPARSE_FORMAT_CONVERSION

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CRS_NONSYM_TO_CRS_SYM, CRS_SYM_TO_CRS_NONSYM, FULL_TO_SPARSE_CRS, SPARSE_CRS_SPARSE_CCS, SPARSE_CRS_TO_FULL, SET_SPARSE_MAT_SYM

   CONTAINS

      SUBROUTINE CRS_NONSYM_TO_CRS_SYM ( NAME_A, NROW_A, NTERM_A, I_A, J_A, A, NAME_B, NTERM_B, I_B, J_B, B )

! Transforms a square symmetric input matrix, A, that is stored in non-symmetric sparse CRS form (i.e. all nonzero terms stored)
! into output matrix B that is stored as symmetric sparse CRS form (i.e. only terms on and above the diagonal stored)

! Input matrix A is stored as:

!      I_A is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          I_A(I+1) - I_A(I) are the number of nonzero terms in row I of matrix A

!      J_A is an integer array giving the col numbers of the NTERM_A nonzero terms in matrix A

!        A is a real array of all of the nonzero terms in matrix A.

! Output matrix B is stored as

!      I_B is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix B. That is:
!          I_B(I+1) - I_B(I) are the number of nonzero terms in row I of matrix B

!      J_B is an integer array giving the col numbers of the nonzero terms in matrix B

!        B is a real array of the nonzero terms on and above the diagonal of input matrix A.

! The number of terms in B (i.e. NTERM_B) is related to NTERM_A and NROW_A as: NTERM_B = (NTERM_A - NROW_A)/2 + NROW_A. This assumes
! a square input matrix. The relationship is not checked herein. In addition, symmetry of the input matrix is assumed and only the
! terms on, and above, the diagonal of A are stored in B.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CRS_NONSYM_TO_CRS_SYM'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME_A            ! Name of input matrix
      CHARACTER(LEN=*), INTENT(IN)    :: NAME_B            ! Name of output matrix

      INTEGER(LONG), INTENT(IN)       :: NROW_A            ! Number of rows in input matrix, A
      INTEGER(LONG), INTENT(IN)       :: NTERM_A           ! Number of nonzero terms in input  matrix, A
      INTEGER(LONG), INTENT(IN)       :: NTERM_B           ! Number of nonzero terms in output matrix, B
      INTEGER(LONG), INTENT(IN)       :: I_A(NROW_A+1)     ! I_A(I+1) - I_A(I) are the number of nonzeros in A row I
      INTEGER(LONG), INTENT(IN)       :: J_A(NTERM_A)      ! Col numbers for nonzero terms in A
      INTEGER(LONG), INTENT(OUT)      :: I_B(NROW_A+1)     ! I_B(I+1) - I_B(I) are the num of nonzeros in B row I
      INTEGER(LONG), INTENT(OUT)      :: J_B(NTERM_B)      ! Col numbers for nonzero terms in B
      INTEGER(LONG)                   :: I,K               ! DO loop indices or counters
      INTEGER(LONG)                   :: KBEG_A            ! Index into array I_A where a row of matrix A begins
      INTEGER(LONG)                   :: KEND_A            ! Index into array I_A where a row of matrix A ends
      INTEGER(LONG)                   :: KTERM_B           ! Count of number of nonzero terms put into output matrix B
      INTEGER(LONG)                   :: A_NTERM_ROW_I     ! Number of terms in a row of input matrix A


      REAL(DOUBLE) , INTENT(IN)       :: A(NTERM_A)        ! Real nonzero values in input  matrix A
      REAL(DOUBLE) , INTENT(OUT)      :: B(NTERM_B)        ! Real nonzero values in output matrix B



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NROW_A+1
         I_B(I) = 0
      ENDDO

      DO I=1,NTERM_B
         J_B(I) = 0
           B(I) = ZERO
      ENDDO

      I_B(1)  = 1
      KTERM_B = 0
      KBEG_A  = 1
      DO I=1,NROW_A
         I_B(I+1) = I_B(I)
         A_NTERM_ROW_I = I_A(I+1) - I_A(I)
         KEND_A = KBEG_A + A_NTERM_ROW_I - 1               ! KBEG_A to KEND_A is the range of indices of terms in A for row I of A
         DO K=KBEG_A,KEND_A
            IF (J_A(K) >= I) THEN                          ! This is a term from A that is on, or above, the diagonal
               KTERM_B = KTERM_B + 1                       ! Increment the counter for the total number of terms in output matrix B
               IF (KTERM_B > NTERM_B) CALL ARRAY_SIZE_ERROR_1( SUBR_NAME, NTERM_B, NAME_B )
               I_B(I+1) = I_B(I+1) + 1                       ! Increment I_B(I+1) for the number of terms in row I of B
               J_B(KTERM_B) = J_A(K)                       ! Set column number for the term to go into B
                 B(KTERM_B) =   A(K)                       ! Set value for the term to go into B
            ENDIF
         ENDDO
         KBEG_A = KEND_A + 1
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CRS_NONSYM_TO_CRS_SYM


      SUBROUTINE CRS_SYM_TO_CRS_NONSYM ( NAME_A, NROW_A, NTERM_A, I_A, J_A, A, NAME_B, NTERM_B, I_B, J_B, B, WRT_SCREEN )

! Transforms a square symmetric input matrix, A, that is stored in symmetric sparse CRS form (i.e. only nonzero terms on and above
! diagonal stored) into output matrix B that is stored as nonsymmetric sparse CRS form (i.e. all nonzero terms stored)

! Input matrix A is stored as:

!      I_A is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix A. That is:
!          I_A(I+1) - I_A(I) are the number of nonzero terms in row I of matrix A

!      J_A is an integer array giving the col numbers of the NTERM_A nonzero terms in matrix A

!        A is a real array of the nonzero terms on and above the diagonal in matrix A.

! Output matrix B is stored as

!      I_B is an array of NROW_A+1 integers that is used to specify the number of nonzero terms in rows of matrix B. That is:
!          I_B(I+1) - I_B(I) are the number of nonzero terms in row I of matrix B

!      J_B is an integer array giving the col numbers of the nonzero terms in matrix B

!        B is a real array of all of the nonzero terms that would be in input matrix A.

! The number of terms in B (i.e. NTERM_B) is related to NTERM_A and NROW_A as: NTERM_B = 2*NTERM_A - NDIAG_A_NZ where
! NDIAG_A_NZ are the number of nonzero diagonal terms in input matriz A. This assumes
! a square input matrix. The relationship is not checked herein. In addition, symmetry of the input matrix is assumed.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CRS_SYM_TO_CRS_NONSYM'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME_A            ! Name of input matrix
      CHARACTER(LEN=*), INTENT(IN)    :: NAME_B            ! Name of output matrix
      CHARACTER(LEN=*), INTENT(IN)    :: WRT_SCREEN        ! If 'Y' then write msgs to screen

      INTEGER(LONG), INTENT(IN)       :: NROW_A            ! Number of rows in input matrix, A
      INTEGER(LONG), INTENT(IN)       :: NTERM_A           ! Number of nonzero terms in input  matrix, A
      INTEGER(LONG), INTENT(IN)       :: NTERM_B           ! Number of nonzero terms in output matrix, B
      INTEGER(LONG), INTENT(IN)       :: I_A(NROW_A+1)     ! I_A(I+1) - I_A(I) are the number of nonzeros in A row I
      INTEGER(LONG), INTENT(IN)       :: J_A(NTERM_A)      ! Col numbers for nonzero terms in A
      INTEGER(LONG), INTENT(OUT)      :: I_B(NROW_A+1)     ! I_B(I+1) - I_B(I) are the num of nonzeros in B row I
      INTEGER(LONG), INTENT(OUT)      :: J_B(NTERM_B)      ! Col numbers for nonzero terms in B
      INTEGER(LONG)                   :: A_NTERM_ROW_I     ! Number of terms in a row of input matrix A
      INTEGER(LONG)                   :: A_ROW_BEG         ! Index into array I_A where a row of matrix A begins
      INTEGER(LONG)                   :: A_ROW_END         ! Index into array I_A where a row of matrix A ends
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: I2_A(NTERM_A)     ! Row numbers of the terms in A
      INTEGER(LONG)                   :: K                 ! Counter


      REAL(DOUBLE) , INTENT(IN)       :: A(NTERM_A)        ! Real nonzero values in input  matrix A
      REAL(DOUBLE) , INTENT(OUT)      :: B(NTERM_B)        ! Real nonzero values in output matrix B

      CHARACTER(LEN=LEN(NAME_A)+7+LEN("Calculating : row")) :: COUNTER_TEMPLATE



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NROW_A+1
         I_B(I) = 0
      ENDDO

      DO I=1,NTERM_B
         J_B(I) = 0
           B(I) = ZERO
      ENDDO

      IF (WRT_SCREEN == 'Y') THEN
   !xx   WRITE(SC1, * )                                    ! Write blank line to screen
      ENDIF

! Generate I2

      K = 0
      DO I=1,NROW_A
         A_NTERM_ROW_I = I_A(I+1) - I_A(I)
         DO J = 1,A_NTERM_ROW_I
            K = K + 1
            I2_A(K) = I
         ENDDO
      ENDDO

      I_B(1) = 1
      K      = 0
      A_ROW_BEG = 1

      WRITE(COUNTER_TEMPLATE, 12345) NAME_B
      CALL COUNTER_INIT(COUNTER_TEMPLATE, NROW_A)
i_do: DO I=1,NROW_A                                        ! Matrix multiply loop. Range over the rows in A

         I_B(I+1) = I_B(I)

         A_NTERM_ROW_I = I_A(I+1) - I_A(I)                 ! Number of terms in matrix A in row I
         A_ROW_END = A_ROW_BEG + A_NTERM_ROW_I - 1         ! A_ROW_BEG to A_ROW_END is range of indices of terms in A for row I of A

         DO J=1,A_ROW_BEG-1                                ! 1st, look for terms that would be in this row, but are not, due to sym
            IF (J_A(J) == I) THEN
               I_B(I+1) = I_B(I+1) + 1
               K        = K + 1
               J_B(K)   = I2_A(J)
                 B(K)   =    A(J)
            ENDIF
         ENDDO

         DO J=A_ROW_BEG,A_ROW_END                          ! 2nd, get terms from this row of A from the diagonal out
            I_B(I+1) = I_B(I+1) + 1
            K        = K + 1
            J_B(K)   = J_A(J)
              B(K)   =   A(J)
         ENDDO

         A_ROW_BEG = A_ROW_END + 1

         CALL COUNTER_PROGRESS(I)
      ENDDO i_do
      WRITE(SC1,*) CR13



      RETURN

! **********************************************************************************************************************************
12345 FORMAT("       Calculating ", A, ": row")

! **********************************************************************************************************************************

      END SUBROUTINE CRS_SYM_TO_CRS_NONSYM


      SUBROUTINE FULL_TO_SPARSE_CRS ( MATIN_NAME, N, M, MATIN_FULL, NTERM_ALLOC, SMALL, CALLING_SUBR, SYM_OUT,                     &
                                      I_MATOUT, J_MATOUT, MATOUT )

! Converts matrices in full format to sparse (compressed row storage) format

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'FULL_TO_SPARSE_CRS'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR         ! Name of subr that called this one
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME           ! Name of matrix
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_OUT              ! 'Y' or 'N' symmetry indicator for output matrix.

      INTEGER(LONG), INTENT(IN)       :: N                    ! Number of rows in input matrix, MATIN_FULL
      INTEGER(LONG), INTENT(IN)       :: M                    ! Number of cols in input matrix, MATIN_FULL
      INTEGER(LONG), INTENT(IN)       :: NTERM_ALLOC          ! Number of nonzero terms allocated to MATOUT in calling subr
      INTEGER(LONG), INTENT(OUT)      :: I_MATOUT(N+1)        ! I_MATOUT(I+1) - I_MATOUT(I) = number of nonzeros in MATOUT row I
      INTEGER(LONG), INTENT(OUT)      :: J_MATOUT(NTERM_ALLOC)! Col numbers for nonzero terms in MATOUT
      INTEGER(LONG)                   :: I,J                  ! DO loop indices
      INTEGER(LONG)                   :: JSTART               ! Starting value for a DO loop
      INTEGER(LONG)                   :: KTERM                ! Counter
      INTEGER(LONG)                   :: ROW_I_NTERMS         ! No. terms in row I of output matrix MATOUT


      REAL(DOUBLE) , INTENT(IN)       :: MATIN_FULL(N,M)      ! Real nonzero values in input matrix MATIN
      REAL(DOUBLE) , INTENT(IN)       :: SMALL                ! Terms < SMALL are filtered out (both here and in calling subr)
      REAL(DOUBLE) , INTENT(OUT)      :: MATOUT(NTERM_ALLOC)  ! Real nonzero values in output matrix MATOUT

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
34568 FORMAT(' I, J, MATIN_FULL(I,J) = ', 2i8, 1es14.6)
! Initialize outputs

      DO I=1,N+1
         I_MATOUT(I) = 0
      ENDDO

      DO I=1,NTERM_ALLOC
         J_MATOUT(I) = 0
           MATOUT(I) = ZERO
      ENDDO

      KTERM = 0
      I_MATOUT(1) = 1
      DO I=1,N
         ROW_I_NTERMS = 0
         I_MATOUT(I+1) = I_MATOUT(I)
         IF      (SYM_OUT == 'Y') THEN
            JSTART = I
         ELSE IF (SYM_OUT == 'N') THEN
            JSTART = 1
         ENDIF
         DO J=JSTART,M
            IF (DABS(MATIN_FULL(I,J)) > SMALL) THEN
               KTERM = KTERM + 1
               ROW_I_NTERMS = ROW_I_NTERMS + 1
               IF (KTERM > NTERM_ALLOC) CALL ARRAY_SIZE_ERROR_1( SUBR_NAME, NTERM_ALLOC, MATIN_NAME )
               I_MATOUT(I+1)   = I_MATOUT(I+1) + 1
               J_MATOUT(KTERM) = J
                 MATOUT(KTERM) = MATIN_FULL(I,J)
            ENDIF
         ENDDO
      ENDDO

      IF (KTERM /= NTERM_ALLOC) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,963) SUBR_NAME, MATIN_NAME, KTERM, NTERM_ALLOC, CALLING_SUBR
         WRITE(F06,963) SUBR_NAME, MATIN_NAME, KTERM, NTERM_ALLOC, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  963 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE CALCD NUMBER OF TERMS PUT INTO MATRIX "',A,'" WAS ',I12,' BUT THERE WERE ',I12,' ALLOCATED TO IT' &
                    ,/,14X,' BY THE CALLING SUBR: ',A)

! **********************************************************************************************************************************

      END SUBROUTINE FULL_TO_SPARSE_CRS


      SUBROUTINE SPARSE_CRS_SPARSE_CCS ( NROWS_A, NCOLS_A, NTERMS_A, MAT_A_NAME, I_A, J_A, A, MAT_B_NAME, J_B, I_B, B, WRT_SCREEN )

! Converts matrices in sparse compressed row storage (CRS) format to sparse compressed column storage (CCS) format

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_CRS_SPARSE_CCS'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_A_NAME        ! Name of input  matrix in CRS format
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_B_NAME        ! Name of output matrix in CCS format
      CHARACTER(LEN=*), INTENT(IN)    :: WRT_SCREEN        ! If 'Y' then write msgs to screen

      INTEGER(LONG), INTENT(IN)       :: NCOLS_A           ! Number of cols in input matrix, A (and output matrix B)
      INTEGER(LONG), INTENT(IN)       :: NROWS_A           ! Number of rows in input matrix, A (and output matrix B)
      INTEGER(LONG), INTENT(IN)       :: NTERMS_A          ! Number of nonzero terms in input matrix, A (and output matrix B)
      INTEGER(LONG), INTENT(IN)       :: I_A(NROWS_A+1)    ! I_A(I+1) - I_A(I) are the number of nonzeros in A row I
      INTEGER(LONG), INTENT(IN)       :: J_A(NTERMS_A)     ! Col numbers for nonzero terms in A
      INTEGER(LONG), INTENT(OUT)      :: I_B(NTERMS_A)     ! Row numbers for nonzero terms in B
      INTEGER(LONG), INTENT(OUT)      :: J_B(NCOLS_A+1)    ! J_B(I+1) - J_B(I) are the number of nonzeros in B col I
      INTEGER(LONG)                   :: I,J,K,L           ! DO loop indices or counters
      INTEGER(LONG)                   :: I2_A(NTERMS_A)    ! Array of row numbers for each term in A
      INTEGER(LONG)                   :: COL_J_NUM_TERMS   ! Number of terms in col J of output matrix B
      INTEGER(LONG)                   :: ROW_I_NUM_TERMS   ! Number of terms in row I of input  matrix A


      REAL(DOUBLE) , INTENT(IN)       :: A(NTERMS_A)       ! Real nonzero values in input  matrix A
      REAL(DOUBLE) , INTENT(OUT)      :: B(NTERMS_A)       ! Real nonzero values in output matrix B

      CHARACTER(LEN=LEN(MAT_A_NAME)+LEN(MAT_B_NAME)+7+LEN("Extracting -> col")) :: COUNTER_TEMPLATE



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NTERMS_A
         I_B(I) = 0
           B(I) = ZERO
      ENDDO

      DO I=1,NCOLS_A+1
         J_B(I) = 0
      ENDDO

      IF (WRT_SCREEN == 'Y') THEN
   !xx   WRITE(SC1, * )                                    ! Write blank line to screen
      ENDIF

      IF ((DEBUG(87) == 1) .OR. (DEBUG(87) == 3)) CALL CRS_CCS_DEB ( '1' )

! Create I2_A array of row numbers for terms in A. These can be deduced from I_A but we can eliminate an internal
! DO loop if we create I2_A from I_A

      K = 0
      DO I=1,NROWS_A
         ROW_I_NUM_TERMS = I_A(I+1) - I_A(I)
         DO J = 1,ROW_I_NUM_TERMS
            K = K + 1
            I2_A(K) = I
         ENDDO
      ENDDO

      L = 0                                                ! Counter for terms going into B
      J_B(1) = 1
      IF (WRT_SCREEN == 'Y') THEN
         WRITE(COUNTER_TEMPLATE, 12345) MAT_A_NAME, MAT_B_NAME
         CALL COUNTER_INIT(COUNTER_TEMPLATE, NCOLS_A)
      END IF
      DO J=1,NCOLS_A
         COL_J_NUM_TERMS = 0
         DO K=1,NTERMS_A
            IF (J_A(K) == J) THEN                          ! We found a term that belongs in col J
               COL_J_NUM_TERMS = COL_J_NUM_TERMS + 1       ! Update the number of terms counted that belong to this column
               L = L + 1
               IF (L > NTERMS_A) CALL ARRAY_SIZE_ERROR_1( SUBR_NAME, NTERMS_A, MAT_B_NAME )
               I_B(L) = I2_A(K)                            ! Array I_B has row numbers of the NTERMS_A terms going into B
                 B(L) = A(K)
            ENDIF
         ENDDO
         IF (WRT_SCREEN == 'Y') THEN
            CALL COUNTER_PROGRESS(J)
         ENDIF
         J_B(J+1) = J_B(J) + COL_J_NUM_TERMS               ! J_B used to tell how many terms there are in each col of B
      ENDDO
      WRITE(SC1,*) CR13

      IF ((DEBUG(87) == 1) .OR. (DEBUG(87) == 3)) CALL CRS_CCS_DEB ( '2' )



      RETURN

! **********************************************************************************************************************************
!12345 FORMAT(7X,'Extracting col  ',I8,' of ',I8,' from matrix ',A,' for matrix ',A,A)
 12345 FORMAT("       Extracting ", A, "->", A, " col")

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE CRS_CCS_DEB ( WHICH )

      CHARACTER( 1*BYTE)              :: WHICH                  ! Decides what to print out for this call to this subr

      INTEGER(LONG)                   :: II                     ! Local DO loop index
      INTEGER(LONG)                   :: NCOLS_B                ! Number of cols in matrix B
      INTEGER(LONG)                   :: NTERMS_B               ! Number of nonzero terms in matrix B

! **********************************************************************************************************************************
      IF      (WHICH == '1') THEN

         WRITE(F06,*)
         WRITE(F06,1011)
         WRITE(F06,1012)
         WRITE(F06,1014) MAT_A_NAME, MAT_B_NAME
         WRITE(F06,1016) NROWS_A, NCOLS_A, NTERMS_A
         WRITE(F06,*)

         WRITE(F06,2800)
         WRITE(F06,*)
         WRITE(F06,3001)
         DO II=1,NROWS_A+1
            WRITE(F06,3002) II,I_A(II)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,3003)
         DO II=1,NTERMS_A
            WRITE(F06,3004) II, J_A(II), A(II)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHICH == '2') THEN

         NCOLS_B  = NCOLS_A
         NTERMS_B = NTERMS_A

         WRITE(F06,2800)
         WRITE(F06,*)
         WRITE(F06,3021)
         DO II=1,NCOLS_B+1
            WRITE(F06,3022) II,J_B(II)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,3023)
         DO II=1,NTERMS_B
            WRITE(F06,3024) II, I_B(II), B(II)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,1095)

      ENDIF

! **********************************************************************************************************************************
 1011 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' :::::::::::::::::::::::::::::::::::::::::START DEBUG(87) OUTPUT FROM SUBROUTINE CRS_CCS:::::::::::::::::::::::'    ,&
              ':::::::::::::::::::::',/)

 1012 FORMAT(' SPARSE MATRIX CRS TO CCS CONVERSION ROUTINE: Convert input matrix A, stored in sparse Compressed Row Storage'      ,&
' (CRS) format, to',/,' -------------------------------------------',/,&
' output matrix B, stored in sparse Compressed Column (CCS) format.',/)

 1014 FORMAT(40X,' The name of CRS formatted matrix A is: ',A                                                                   ,/,&
             40x,' The name of CCS formatted matrix B is: ',A,/)

 1016 FORMAT(26X,' Matrix A  has ',I8,' rows and ',I8,' cols and '  ,I12,' nonzero terms',/)

 2800 FORMAT(1X,'****************************************************************************************************************',&
                '*******************')

 3001 FORMAT(' Compressed Row Storage (CRS) format of input matrix A:'                                                          ,/,&
             ' ------------------------------------------------------'                                                         ,//,&
             ' 1) Index, L, and array I_A(L) for matrix A, where I_A(L+1) - I_A(L) is the number of nonzero terms in row L of'    ,&
             ' matrix A.',/,'    (also, I_A(L) is the index, K, in array A(K) where row L begins - up to, but not including, the' ,&
             ' last entry in I_A(L)).',/)

 3002 FORMAT('    L, I_A(L)       = ',2I12)

 3003 FORMAT(' 2) Index, K, and arrays J_A(K) and A(K). A(K) are the nonzeros in matrix A and J_A(K) is the col number in matrix', &
             ' A for term A(K).',/)

 3004 FORMAT('    K, J_A(K), A(K) = ',2I12,1ES15.6)

 3021 FORMAT(' Compressed Col Storage (CCS) format of output matrix B:'                                                         ,/,&
             ' -------------------------------------------------------'                                                        ,//,&
             ' 1) Index, L, and array J_B(L) for matrix B, where J_B(L+1) - J_B(L) is the number of nonzero terms in row L of'    ,&
             ' matrix B.',/,'    (also, J_B(L) is the index, K in array B(K), where row L begins - up to, but not including, the' ,&
             ' last entry in J_B(L)).',/)

 3022 FORMAT('    L, J_B(L)       = ',2I12)

 3023 FORMAT(' 2) Index, K, and arrays I_B(K) and B(K). B(K) are the nonzeros in matrix B and I_B(K) is the col number in matrix', &
             ' B for term B(K).',/)

 3024 FORMAT('    K, I_B(K), B(K) = ',2I12,1ES15.6)

 1095 FORMAT(' ::::::::::::::::::::::::::::::::::::::::END DEBUG(87) OUTPUT FROM SUBROUTINE CRS_CCS::::::::::::::::::::::::::',    &
              ':::::::::::::::::::::'                                                                                           ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE CRS_CCS_DEB

      END SUBROUTINE SPARSE_CRS_SPARSE_CCS


      SUBROUTINE SPARSE_CRS_TO_FULL ( MATIN_NAME, NTERM_IN, NROWS, NCOLS, SYM_IN, I_MATIN, J_MATIN, MATIN, MATOUT )

! Converts matrices in sparse compressed row storage format to full format

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE CONSTANTS_1, ONLY           :  ZERO

      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1
      USE RESULT_FORMATTING, ONLY     :  WRT_REAL_TO_CHAR_VAR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_CRS_TO_FULL'
      CHARACTER(LEN=*), INTENT(IN)    :: SYM_IN              ! 'Y' or 'N' symmetry indicator for input matrix.
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME          ! Name of matrix

      INTEGER(LONG), INTENT(IN)       :: NCOLS               ! Number of cols in input matrix, MATIN
      INTEGER(LONG), INTENT(IN)       :: NROWS               ! Number of rows in input matrix, MATIN
      INTEGER(LONG), INTENT(IN)       :: NTERM_IN            ! Number of nonzero terms in input matrix, MATIN
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(NROWS+1)    ! I_MATIN(I+1) - I_MATIN(I) are the number of nonzeros in MATIN row I
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERM_IN)   ! Col numbers for nonzero terms in MATIN
      INTEGER(LONG)                   :: I,J,K               ! DO loop indices or counters
      INTEGER(LONG)                   :: ROW_I_NTERMS        ! No. terms in row I of input matrix MATIN


      REAL(DOUBLE) , INTENT(IN)       :: MATIN(NTERM_IN)     ! Real nonzero values in input  matrix MATIN
      REAL(DOUBLE) , INTENT(OUT)      :: MATOUT(NROWS,NCOLS) ! Real nonzero values in output matrix MATOUT



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NROWS
         DO J=1,NCOLS
            MATOUT(I,J) = ZERO
         ENDDO
      ENDDO

! Calc outputs

      DO I=1,NROWS
         DO J=1,NCOLS
            MATOUT(I,J) = ZERO
         ENDDO
      ENDDO

! Create full matrix MATOUT from sparse input matrix MATIN

      K = 0
      DO I=1,NROWS
         ROW_I_NTERMS = I_MATIN(I+1) - I_MATIN(I)
         DO J=1,ROW_I_NTERMS
            K = K + 1
            IF (K > NTERM_IN) CALL ARRAY_SIZE_ERROR_1( SUBR_NAME, NTERM_IN, MATIN_NAME )
            MATOUT(I,J_MATIN(K)) = MATIN(K)
         ENDDO
      ENDDO

! If input matrix was tagged as symmetric, then lower triang portion was not in MATIN, so set lower triang portion:

      IF (SYM_IN == 'Y') THEN
         DO I=1,NROWS
            DO J=1,I-1
               MATOUT(I,J) = MATOUT(J,I)
            ENDDO
         ENDDO
      ENDIF


      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE SPARSE_CRS_TO_FULL


      SUBROUTINE SET_SPARSE_MAT_SYM

! Sets symmetry indicators for sparse matrices depending on Bulk Data PARAM SPARSTOR

      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06

      USE PARAMS, ONLY                :  SPARSTOR, SUPINFO

      USE SPARSE_MATRICES, ONLY       :  SYM_KGG    , SYM_MGG    , SYM_MGGC   , SYM_MGGE   , SYM_MGGS   , SYM_PG     , SYM_RMG    ,&
                                         SYM_KGGD

      USE SPARSE_MATRICES, ONLY       :  SYM_KNN    , SYM_KNM    , SYM_KMM    , SYM_KMN    ,                                       &
                                         SYM_KNND   , SYM_KNMD   , SYM_KMMD   , SYM_KMND   ,                                       &
                                         SYM_MNN    , SYM_MNM    , SYM_MMN    , SYM_MMM    ,                                       &
                                         SYM_PN     , SYM_PM     ,                                                                 &
                                         SYM_RMN    , SYM_RMM    , SYM_GMN    , SYM_GMNt   , SYM_HMN    , SYM_LMN

      USE SPARSE_MATRICES, ONLY       :  SYM_KFF    , SYM_KFS    , SYM_KSF    , SYM_KSS    , SYM_KFSe   , SYM_KSSe   ,             &
                                         SYM_KFFD   , SYM_KFSD   , SYM_KSFD   , SYM_KSSD   , SYM_KFSDe  , SYM_KSSDe  ,             &
                                         SYM_MFF    , SYM_MSF    , SYM_MFS    , SYM_MSS    ,                                       &
                                         SYM_PF     , SYM_PF_TMP , SYM_PFYS   , SYM_PFYS1  , SYM_PS     , SYM_QSYS

      USE SPARSE_MATRICES, ONLY       :  SYM_KAA    , SYM_KAO    , SYM_KOO    , SYM_KMSM   , SYM_KMSMn  ,                          &
                                         SYM_KAAD   , SYM_KAOD   , SYM_KOOD   ,                                                    &
                                         SYM_MAA    , SYM_MAO    , SYM_MOO    ,                                                    &
                                         SYM_PA     , SYM_PO     ,                                                                 &
                                         SYM_GOA    , SYM_GOAt

      USE SPARSE_MATRICES, ONLY       :  SYM_KLL    , SYM_KRL    , SYM_KRR    , SYM_KLLs   , SYM_KRRcb  , SYM_KRRcbn ,&
                                         SYM_KRRcbs , SYM_KXX    , SYM_KLLD   , SYM_KRLD   , SYM_KRRD   , SYM_KLLDs  ,             &
                                         SYM_MPF0   , SYM_MLL    , SYM_MLLn   , SYM_MRL    , SYM_MLR    , SYM_MRR    , SYM_MLLs   ,&
                                         SYM_MRN    , SYM_MRRcb  , SYM_MRRcbn , SYM_MXX    , SYM_MXXn   ,                          &
                                         SYM_PL     , SYM_PR     ,                                                                 &
                                         SYM_DLR    , SYM_DLRt   , SYM_IRR    , SYM_PHIXA  , SYM_IF_LTM , SYM_CG_LTM , SYM_PHIZL  ,&
                                         SYM_PHIZL1 , SYM_PHIZL2 , SYM_LTM    , SYM_PHIZL1t

      IMPLICIT NONE

! **********************************************************************************************************************************
      WRITE(ERR,101) SPARSTOR
      IF (SUPINFO == 'N') THEN
         WRITE(F06,101) SPARSTOR
      ENDIF

      IF (SPARSTOR == 'SYM   ') THEN

         SYM_KGG      = 'Y'
         SYM_KGGD     = 'Y'
         SYM_MGG      = 'Y'
         SYM_MGGC     = 'Y'
         SYM_MGGE     = 'Y'
         SYM_MGGS     = 'Y'
         SYM_PG       = 'N'
         SYM_RMG      = 'N'

         SYM_KNN      = 'Y'
         SYM_KNM      = 'N'
         SYM_KMM      = 'Y'
         SYM_KMN      = 'N'
         SYM_KNND     = 'Y'
         SYM_KNMD     = 'N'
         SYM_KMMD     = 'Y'
         SYM_KMND     = 'N'
         SYM_MNN      = 'Y'
         SYM_MNM      = 'N'
         SYM_MMN      = 'N'
         SYM_MMM      = 'Y'
         SYM_PN       = 'N'
         SYM_PM       = 'N'
         SYM_RMN      = 'N'
         SYM_RMM      = 'N'
         SYM_GMN      = 'N'
         SYM_GMNt     = 'N'
         SYM_HMN      = 'N'
         SYM_LMN      = 'N'

         SYM_KFF      = 'Y'
         SYM_KFS      = 'N'
         SYM_KSS      = 'Y'
         SYM_KFSe     = 'N'
         SYM_KSSe     = 'N'
         SYM_KFFD     = 'Y'
         SYM_KFSD     = 'N'
         SYM_KSSD     = 'Y'
         SYM_KFSDe    = 'N'
         SYM_KSSDe    = 'N'
         SYM_MFF      = 'Y'
         SYM_MFS      = 'N'
         SYM_MSF      = 'N'
         SYM_MSS      = 'Y'
         SYM_PF       = 'N'
         SYM_PF_TMP   = 'N'
         SYM_PFYS     = 'N'
         SYM_PFYS1    = 'N'
         SYM_PS       = 'N'
         SYM_QSYS     = 'N'

         SYM_KAA      = 'Y'
         SYM_KAO      = 'N'
         SYM_KOO      = 'Y'
         SYM_KAAD     = 'Y'
         SYM_KAOD     = 'N'
         SYM_KOOD     = 'Y'
         SYM_MAA      = 'Y'
         SYM_MAO      = 'N'
         SYM_MOO      = 'Y'
         SYM_PA       = 'N'
         SYM_PO       = 'N'
         SYM_GOA      = 'N'
         SYM_GOAt     = 'N'

         SYM_KLL      = 'Y'
         SYM_KLLs     = 'Y'
         SYM_KRL      = 'N'
         SYM_KRR      = 'Y'
         SYM_KLLD     = 'Y'
         SYM_KLLDs    = 'Y'
         SYM_KRLD     = 'N'
         SYM_KRRD     = 'Y'
         SYM_MPF0     = 'N'
         SYM_MLL      = 'Y'
         SYM_MLLn     = 'N'
         SYM_MLLs     = 'Y'
         SYM_MLR      = 'N'
         SYM_MRL      = 'N'
         SYM_MRR      = 'Y'
         SYM_DLR      = 'N'
         SYM_DLRt     = 'N'
         SYM_PHIZL    = 'N'
         SYM_PHIZL1   = 'N'
         SYM_PHIZL1t  = 'N'
         SYM_PHIZL2   = 'N'
         SYM_CG_LTM   = 'N'
         SYM_IF_LTM   = 'N'
         SYM_LTM      = 'N'
         SYM_IRR      = 'Y'
         SYM_PHIXA    = 'N'
         SYM_KRRcb    = 'Y'
         SYM_KRRcbn   = 'N'
         SYM_KRRcbs   = 'Y'
         SYM_KXX      = 'Y'
         SYM_MRN      = 'N'
         SYM_MRRcb    = 'Y'
         SYM_MRRcbn   = 'N'
         SYM_MXX      = 'Y'
         SYM_MXXn     = 'N'
         SYM_PL       = 'N'
         SYM_PR       = 'N'

         SYM_KMSM     = 'Y'
         SYM_KMSMn    = 'N'

      ELSE IF (SPARSTOR == 'NONSYM') THEN

         SYM_KGG      = 'N'
         SYM_KGGD     = 'N'
         SYM_MGG      = 'N'
         SYM_MGGC     = 'N'
         SYM_MGGe     = 'N'
         SYM_MGGS     = 'N'
         SYM_PG       = 'N'
         SYM_RMG      = 'N'

         SYM_KNN      = 'N'
         SYM_KNM      = 'N'
         SYM_KMM      = 'N'
         SYM_KMN      = 'N'
         SYM_KNND     = 'N'
         SYM_KNMD     = 'N'
         SYM_KMMD     = 'N'
         SYM_KMND     = 'N'
         SYM_MNN      = 'N'
         SYM_MNM      = 'N'
         SYM_MMN      = 'N'
         SYM_MMM      = 'N'
         SYM_PN       = 'N'
         SYM_PM       = 'N'
         SYM_RMN      = 'N'
         SYM_RMM      = 'N'
         SYM_GMN      = 'N'
         SYM_GMNt     = 'N'
         SYM_HMN      = 'N'
         SYM_LMN      = 'N'

         SYM_KFF      = 'N'
         SYM_KFS      = 'N'
         SYM_KSS      = 'N'
         SYM_KFSe     = 'N'
         SYM_KSSe     = 'N'
         SYM_KFFD     = 'N'
         SYM_KFSD     = 'N'
         SYM_KSSD     = 'N'
         SYM_KFSDe    = 'N'
         SYM_KSSDe    = 'N'
         SYM_MFF      = 'N'
         SYM_MFS      = 'N'
         SYM_MSF      = 'N'
         SYM_MSS      = 'N'
         SYM_PF       = 'N'
         SYM_PF_TMP   = 'N'
         SYM_PFYS     = 'N'
         SYM_PFYS1    = 'N'
         SYM_PS       = 'N'
         SYM_QSYS     = 'N'

         SYM_KAA      = 'N'
         SYM_KAO      = 'N'
         SYM_KOO      = 'N'
         SYM_KAAD     = 'N'
         SYM_KAOD     = 'N'
         SYM_KOOD     = 'N'
         SYM_MAA      = 'N'
         SYM_MAO      = 'N'
         SYM_MOO      = 'N'
         SYM_PA       = 'N'
         SYM_PO       = 'N'
         SYM_GOA      = 'N'
         SYM_GOAt     = 'N'

         SYM_KLL      = 'N'
         SYM_KLLs     = 'Y'                                ! KLLs is always symmetric
         SYM_KRL      = 'N'
         SYM_KRR      = 'N'
         SYM_KLLD     = 'N'
         SYM_KLLDs    = 'Y'                                ! KLLs is always symmetric
         SYM_KRLD     = 'N'
         SYM_KRRD     = 'N'
         SYM_MPF0     = 'N'
         SYM_MLL      = 'N'
         SYM_MLLn     = 'N'
         SYM_MLLs     = 'Y'                                ! MLLs is always symmetric
         SYM_MLR      = 'N'
         SYM_MRL      = 'N'
         SYM_MRR      = 'N'
         SYM_DLR      = 'N'
         SYM_DLRt     = 'N'
         SYM_PHIZL    = 'N'
         SYM_PHIZL1   = 'N'
         SYM_PHIZL1t  = 'N'
         SYM_PHIZL2   = 'N'
         SYM_CG_LTM   = 'N'
         SYM_IF_LTM   = 'N'
         SYM_LTM      = 'N'
         SYM_IRR      = 'N'
         SYM_PHIXA    = 'N'
         SYM_KRRcb    = 'N'
         SYM_KRRcbn   = 'N'
         SYM_KRRcbs   = 'Y'
         SYM_KXX      = 'N'
         SYM_MRN      = 'N'
         SYM_MRRcb    = 'N'
         SYM_MRRcbn   = 'N'
         SYM_MXX      = 'N'
         SYM_MXXn     = 'N'
         SYM_PL       = 'N'
         SYM_PR       = 'N'

         SYM_KMSM     = 'N'
         SYM_KMSMn    = 'N'

      ENDIF

! **********************************************************************************************************************************
  101 FORMAT(' *INFORMATION: SPARSE MATRICES ARE STORED IN ',A,' FORMAT',/)

! **********************************************************************************************************************************

      END SUBROUTINE SET_SPARSE_MAT_SYM

   END MODULE SPARSE_FORMAT_CONVERSION
