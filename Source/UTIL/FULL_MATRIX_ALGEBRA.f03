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

   MODULE FULL_MATRIX_ALGEBRA

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CHECK_MAT_INVERSE, INVERT_FF_MAT, MATADD_FFF, MATMULT_FFF, MATMULT_FFF_T, OUTER_PRODUCT

   CONTAINS

      SUBROUTINE CHECK_MAT_INVERSE ( MAT_NAME, A, AI, NSIZE )

! Checks whether a matrix (A) times its calculated inverse (AI) is the identity matrix

      USE PENTIUM_II_KIND, ONLY        : BYTE, DOUBLE, LONG
      USE IOUNT1, ONLY                 : ERR, F06
      USE CONSTANTS_1, ONLY            : ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                 :  EPSIL

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME          ! Name of input matrix
      CHARACTER(11*BYTE)              :: COL(3)

      INTEGER(LONG), INTENT(IN)       :: NSIZE             ! Row/col size of input matrices
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! Error indicator:
!                                                            IERR = 2 indicates SUM_DIAGS     > SMALL_NUM
!                                                            IERR = 3 indicates SUM_OFF_DIAGS > SMALL_NUM
!                                                            IERR = 5 indicates SUM_DIAGS and SUM_OFF_DIAGS > SMALL_NUM

      REAL(DOUBLE) , INTENT(IN)       :: A(NSIZE,NSIZE)    ! Matrix to invert
      REAL(DOUBLE) , INTENT(IN)       :: AI(NSIZE,NSIZE)   ! Inverse of A
      REAL(DOUBLE)                    :: IDENT(NSIZE,NSIZE)! A*AI should be the identity matrix
      REAL(DOUBLE)                    :: SMALL_NUM         ! Small number for comparing to near zero
      REAL(DOUBLE)                    :: SUM_DIAGS         ! Sum of diag terms from IDENT
      REAL(DOUBLE)                    :: SUM_OFF_DIAGS     ! Sum of off-diag terms from IDENT

! **********************************************************************************************************************************
! Initialize

      IERR = 0
      SMALL_NUM = EPSIL(1)

      CALL MATMULT_FFF ( A, AI, NSIZE, NSIZE, NSIZE, IDENT )

      IF (DEBUG(199) > 1) THEN

         WRITE(F06,1001) MAT_NAME                             ! Write input matrix
         DO I=1,NSIZE
            COL(I) = '        Col'
         ENDDO
         WRITE(F06,2001) (COL(I),I, I=1,NSIZE)
         DO I=1,NSIZE
            WRITE(F06,2002) I, (A(I,J),J=1,NSIZE)
         ENDDO
         WRITE (F06,*)

         WRITE(F06,1002)                                      ! Write inverse of input matrix
         DO I=1,NSIZE
            COL(I) = '        Col'
         ENDDO
         WRITE(F06,2001) (COL(I),I, I=1,NSIZE)
         DO I=1,NSIZE
            WRITE(F06,2002) I, (AI(I,J),J=1,NSIZE)
         ENDDO
         WRITE (F06,*)

         WRITE(F06,1003)                                      ! Write IDENT matrix
         DO I=1,NSIZE
            COL(I) = '        Col'
         ENDDO
         WRITE(F06,2001) (COL(I),I, I=1,NSIZE)
         DO I=1,NSIZE
            WRITE(F06,2002) I, (IDENT(I,J),J=1,NSIZE)
         ENDDO
         WRITE (F06,*)

      ENDIF

      SUM_DIAGS = ZERO
      DO I=1,NSIZE
         SUM_DIAGS = SUM_DIAGS + IDENT(I,I)
      ENDDO

      SUM_OFF_DIAGS = ZERO
      DO I=1,NSIZE
         DO J=1,NSIZE
            IF (J /= I) THEN
               SUM_OFF_DIAGS = SUM_OFF_DIAGS + IDENT(I,J)
            ENDIF
         ENDDO
      ENDDO

      WRITE(F06,4001) SUM_DIAGS
      WRITE(F06,4002) SUM_DIAGS/NSIZE
      WRITE(F06,4003) 100.D0*SUM_OFF_DIAGS
      WRITE(F06,*)

      IF ( DABS(SUM_DIAGS - NSIZE) > SMALL_NUM ) THEN
         IERR = 2
      ENDIF

      IF ( DABS(SUM_OFF_DIAGS) > SMALL_NUM ) THEN
         IERR = IERR + 3
      ENDIF

      IF (IERR > 0) THEN
         WRITE(ERR,9000) MAT_NAME
         WRITE(F06,9000) MAT_NAME
      ENDIF

      IF      (IERR == 2) THEN
         WRITE(ERR,9002)
         WRITE(F06,9002)
      ELSE IF (IERR == 3) THEN
         WRITE(ERR,9003)
         WRITE(F06,9003)
      ELSE IF (IERR == 5) THEN
         WRITE(ERR,9002)
         WRITE(F06,9002)
         WRITE(ERR,9003)
         WRITE(F06,9003)
      ENDIF

      RETURN

!**********************************************************************************************************************************
 1001 FORMAT(' ***************************************************************************************************************',/,&
             ' Input matrix is A with name =',A,//,' Input matrix A:',/,' --------------')

 1002 FORMAT(' AI = inverse of the input matrix:',/,' ---------------------------------')

 1003 FORMAT(' IDENT = A*AI:',/,' ------------')

 2001 FORMAT(5X,3(A,I3))

 2002 FORMAT(' Row',I3,':',3(1ES14.6))

 4001 FORMAT(' Sum of diagonal terms in IDENT = ',1ES14.6)

 4002 FORMAT(' Avg    diagonal term  in IDENT = ',1ES14.6)

 4003 FORMAT(' Sum of off diag terms in IDENT = ',1ES14.2,'% of a unity term in an identity matrix (which A*AI should be)')

 9000 FORMAT(' *INFORMATION: MATRIX ',A,' CANNOT BE INVERTED, OR THE INVERSE FAILED THE CHECK THAT IT TIMES ITS INVERSE IS THE',   &
                           ' IDENTITY MATRIX, DUE TO THE FOLLOWING:')

 9002 FORMAT('               IT FAILED THE INVERT CHECK IN THAT THE SUM OF THE DIAG TERMS IN THE 3x3 CHECK MATRIX WAS NOT = 3')

 9003 FORMAT('               IT FAILED THE INVERT CHECK IN THAT THE SUM OF THE OFF-DIAG TERMS IN THE 3x3 CHECK MATRIX WAS NOT = 0')

!**********************************************************************************************************************************

      END SUBROUTINE CHECK_MAT_INVERSE



      SUBROUTINE INVERT_FF_MAT ( CALLING_SUBR, MAT_A_NAME, A, NROWS, INFO )

! Invert symmetric matrix A which is stored in full format. The return has the inverse of the matrix in array A

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE, LONG
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'INVERT_FF_MAT'
      CHARACTER(LEN=*), INTENT(IN )   :: CALLING_SUBR      ! Name of subr that called this subr
      CHARACTER(LEN=*), INTENT(IN )   :: MAT_A_NAME        ! Name of input matrix to be inverted

      INTEGER(LONG)   , INTENT(IN)    :: NROWS             ! Row/col size of input matrix A
      INTEGER(LONG)   , INTENT(OUT)   :: INFO              ! Output from LAPACK routines to do factorization of Lapack band matrix
                                                           !   0:  successful exit
                                                           ! < 0:  if INFO = -i, the i-th argument had an illegal value
                                                           ! > 0:  if INFO =  i, the leading minor of order i is not pos definite
      INTEGER(LONG)                   :: I,J               ! DO loop indices

      REAL(DOUBLE)    , INTENT(INOUT) :: A(NROWS,NROWS)    ! Matrix to invert. Inverted matrix returned in A



! **********************************************************************************************************************************
! In DPOTRF A has the matrix to invert as input and the triangular factor of A coming out. In DPOTRI A has the tria factor of the
! matrix to invert going in and the inverse of the matrix coming out

      CALL DPOTRF ( 'U', NROWS, A, NROWS, INFO )           ! Factor A

      IF (INFO > 0) THEN                                ! Error factoring A

         WRITE(ERR,994) MAT_A_NAME, INFO, CALLING_SUBR
         WRITE(F06,994) MAT_A_NAME, INFO, CALLING_SUBR

      ELSE

         CALL DPOTRI ( 'U', NROWS, A, NROWS, INFO )     ! No error factoring A so invert

         IF (INFO > 0) THEN                             ! Error inverting A

            WRITE(ERR,995) MAT_A_NAME, INFO, CALLING_SUBR
            WRITE(F06,995) MAT_A_NAME, INFO, CALLING_SUBR

         ELSE                                           ! Matrix inverted OK so set lower triangle of A based on symmetry

            DO I=1,NROWS
               DO J=1,I-1
                  A(I,J) = A(J,I)
               ENDDO
            ENDDO

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
  994 FORMAT(' *ERROR   994: THE FACTORIZATION OF MATRIX NAMED "',A,'" COULD NOT BE COMPLETED BY LAPACK SUBR DPOTRF '              &
                    ,/,14X,' A SINGULARITY WAS FOUND IN ROW ',I8,'. THIS SUBR WAS CALLED BY SUBR ',A)

  995 FORMAT(' *ERROR   995: INVERSION OF MATRIX ',A,' COULD NOT BE COMPLETED BY LAPACK SUBR DPOTRI '                              &
                    ,/,14X,' THE DIAG TERM IN ROW ',I8,' OF THE TRIANG FACTOR OF THE MATRIX IS 0. THIS SUBR WAS CALLED BY SUBR ',A)

98764 format(32767(1es14.6))

98765 format(32767a14)

! **********************************************************************************************************************************

      END SUBROUTINE INVERT_FF_MAT



      SUBROUTINE MATADD_FFF ( A, B, NROW, NCOL, ALPHA, BETA, ITRNSPB, C)

! Adds two matrices: A + B (if ITRNSPB = 0), or A + B' (if ITRNSPB = 1 and A and B are square).
! Returns result, matrix C. All matrices are in full format
! User must make certain that matrices A and B have the same number of rows and cols

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATADD_FFF'

      INTEGER(LONG), INTENT(IN)       :: NROW              ! Number of rows in matrces A, B, C
      INTEGER(LONG), INTENT(IN)       :: NCOL              ! Number of cols in matrces A, B, C
      INTEGER(LONG), INTENT(IN)       :: ITRNSPB           ! Transpose indicator for matrix B
      INTEGER(LONG)                   :: I,J               ! DO loop indices or counters


      REAL(DOUBLE) , INTENT(IN)       :: A(NROW,NCOL)      ! Input  matrix A
      REAL(DOUBLE) , INTENT(IN)       :: B(NROW,NCOL)      ! Input  matrix B
      REAL(DOUBLE) , INTENT(IN)       :: ALPHA             ! Scalar multiplier for matrix A
      REAL(DOUBLE) , INTENT(IN)       :: BETA              ! Scalar multiplier for matrix B

      REAL(DOUBLE) , INTENT(OUT)      :: C(NROW,NCOL)      ! Output matrix C



! **********************************************************************************************************************************


! Initialize outputs

      DO I=1,NROW
         DO J=1,NCOL
            C(I,J) = ZERO
         ENDDO
      ENDDO

! Check for coding error

      IF ((ITRNSPB /=0) .AND. (ITRNSPB /=1)) THEN
         WRITE(ERR,933) SUBR_NAME,ITRNSPB
         WRITE(F06,933) SUBR_NAME,ITRNSPB
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF ((ITRNSPB == 1) .AND. (NROW /= NCOL)) THEN
         WRITE(ERR,934) SUBR_NAME,ITRNSPB, NROW, NCOL
         WRITE(F06,934) SUBR_NAME,ITRNSPB, NROW, NCOL
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Add matrices A and B (or A and B')

      IF (ITRNSPB == 0) THEN                               ! Add A, B

         DO I =1,NROW
            DO J = 1,NCOL
               C(I,J) = ALPHA*A(I,J) + BETA*B(I,J)
            ENDDO
         ENDDO

      ELSE                                                 ! Add A, B'

         DO I =1,NROW
            DO J = 1,NCOL
               C(I,J) = ALPHA*A(I,J) + BETA*B(J,I)
            ENDDO
         ENDDO

      ENDIF




      RETURN

! **********************************************************************************************************************************
  933 FORMAT(' *ERROR   933: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INPUT PARAMETER ITRNSPB MUST BE 0 OR 1 BUT VALUE IS ',I8)

  934 FORMAT(' *ERROR   934: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' WHEN INPUT PARAMETER ITRNSPB =',I3,' THEN INPUT PARAMETERS NROW AND NCOL MUST BE EQUAL. '             &
                    ,/,14X,' HOWEVER, AS INPUT TO THIS SUBR, NROW = ',I8,' AND NCOL = ',I8)


! **********************************************************************************************************************************

      END SUBROUTINE MATADD_FFF


      SUBROUTINE MATMULT_FFF ( A, B, NROWA, NCOLA, NCOLB, C )

! Multiplies two matrices: A x B. Returns result, matrix C. All matrices are in full format
! NOTE: User is responsible for making sure that A and B are conformable

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULT_FFF'

      INTEGER(LONG), INTENT(IN)       :: NROWA             ! No. rows in input matrix A
      INTEGER(LONG), INTENT(IN)       :: NCOLA             ! No. cols in input matrix A
      INTEGER(LONG), INTENT(IN)       :: NCOLB             ! No. cols in input matrix B
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: NROWB             !


      REAL(DOUBLE) , INTENT(IN)       :: A(NROWA,NCOLA)    ! Input  matrix A
      REAL(DOUBLE) , INTENT(IN)       :: B(NCOLA,NCOLB)    ! Input  matrix B
      REAL(DOUBLE) , INTENT(OUT)      :: C(NROWA,NCOLB)    ! Output matrix C



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NROWA
         DO J=1,NCOLB
            C(I,J) = ZERO
         ENDDO
      ENDDO

      NROWB = NCOLA

! Multiply A x B

      DO I =1,NROWA
         DO J = 1,NCOLB
            C(I,J) = ZERO
            DO K = 1,NROWB
               C(I,J) = C(I,J) + A(I,K)*B(K,J)
            ENDDO
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MATMULT_FFF


      SUBROUTINE MATMULT_FFF_T ( A, B, NROWA, NCOLA, NCOLB, C )

! Multiplies two matrices: A' x B (A is transposed). Returns result, matrix C. All matrices are in full format
! NOTE: User is responsible for making sure that A(t) and B are conformable

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MATMULT_FFF_T'

      INTEGER(LONG), INTENT(IN)       :: NROWA             ! No. rows in input matrix A (NOT A')
      INTEGER(LONG), INTENT(IN)       :: NCOLA             ! No. cols in input matrix A (NOT A')
      INTEGER(LONG), INTENT(IN)       :: NCOLB             ! No. cols in input matrix B
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: NROWB             ! No. rows in input matrix B
      INTEGER(LONG)                   :: NROWA_T           ! No. rows in A' (ranspose)
      INTEGER(LONG)                   :: NCOLA_T           ! No. cols in A' (ranspose)


      REAL(DOUBLE) , INTENT(IN)       :: A(NROWA,NCOLA)    ! Input  matrix A
      REAL(DOUBLE) , INTENT(IN)       :: B(NROWA,NCOLB)    ! Input  matrix B
      REAL(DOUBLE) , INTENT(OUT)      :: C(NCOLA,NCOLB)    ! Output matrix C



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NCOLA
         DO J=1,NCOLB
            C(I,J) = ZERO
         ENDDO
      ENDDO

      NROWA_T = NCOLA
      NCOLA_T = NROWA
      NROWB   = NCOLA_T

! Multiply A' x B

      DO I =1,NROWA_T
         DO J = 1,NCOLB
            C(I,J) = ZERO
            DO K = 1,NROWB
               C(I,J) = C(I,J) + A(K,I)*B(K,J)
            ENDDO
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE MATMULT_FFF_T


      SUBROUTINE OUTER_PRODUCT ( A, B, NA, NB, C )

! Computes the outer product of two vectors: A and B. Returns result, matrix C.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: NA                ! No. elements in input vector A
      INTEGER(LONG), INTENT(IN)       :: NB                ! No. elements in input vector B
      INTEGER(LONG)                   :: I,J               ! DO loop indices

      REAL(DOUBLE) , INTENT(IN)       :: A(NA)             ! Input  vector A
      REAL(DOUBLE) , INTENT(IN)       :: B(NB)             ! Input  vector B
      REAL(DOUBLE) , INTENT(OUT)      :: C(NA,NB)          ! Output matrix C

! *********************************************************************************************************************************

      DO I =1,NA
        DO J = 1,NB
          C(I,J) = A(I)*B(J)
        ENDDO
      ENDDO

! **********************************************************************************************************************************

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE OUTER_PRODUCT

   END MODULE FULL_MATRIX_ALGEBRA
