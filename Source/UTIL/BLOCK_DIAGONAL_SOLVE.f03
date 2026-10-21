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

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
! LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
! NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.

      MODULE BLOCK_DIAGONAL_SOLVE

! Solves A*x = b for a square, unsymmetric, sparse A (CRS) that falls apart into many independent blocks, with sparse right hand
! sides. Made for the constraint matrix RMM of SOLVE_GMN: each rigid element or MPC couples only its own dependent DOFs, so RMM is
! block diagonal, and each column of RMN touches a few blocks.
!   BDS_FACTOR   finds the blocks (connected components of the structure of A + A'), factors each block of up to NDENSE_MAX
!                DOFs with LAPACK DGETRF, and all larger blocks together with SuperLU;
!   BDS_SOLVE    solves for one sparse b, working only in the blocks that b touches, and returns x on those blocks;
!   BDS_FREE     frees the factors.
! The work for one b is the size of the blocks it touches, not the size of A.

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE, DBL_LONG

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: BDS_FACTOR, BDS_SOLVE, BDS_FREE

      INTEGER(LONG), PARAMETER        :: NDENSE_MAX = 200  ! Largest block factored as a dense matrix

      INTEGER(LONG)                   :: N = 0             ! Order of A
      INTEGER(LONG)                   :: NCOMP = 0         ! Number of blocks
      INTEGER(LONG), ALLOCATABLE      :: COMP(:)           ! Block of each row
      INTEGER(LONG), ALLOCATABLE      :: LOC(:)            ! Index of each row in its block
      INTEGER(LONG), ALLOCATABLE      :: CPTR(:), CDOF(:)  ! Rows of each block, ascending: CDOF(CPTR(c):CPTR(c+1)-1)
      INTEGER(DBL_LONG), ALLOCATABLE  :: LUOFF(:)          ! Offset of the dense LU of each block in LU (0: in the SuperLU part)
      REAL(DOUBLE) , ALLOCATABLE      :: LU(:)             ! Dense LU factors of the small blocks, column by column
      INTEGER      , ALLOCATABLE      :: PIV(:)            ! Their pivots, by row (at CPTR(c) ...)

      INTEGER(LONG)                   :: NBIG = 0          ! Rows in blocks larger than NDENSE_MAX (factored together)
      INTEGER(LONG)                   :: NBIGC = 0         ! Number of such blocks
      INTEGER(LONG), ALLOCATABLE      :: BIGROW(:)         ! Row of A of each of them; BIGLOC(i): index of row i among them, or 0
      INTEGER(LONG), ALLOCATABLE      :: BIGLOC(:)
      INTEGER(LONG)                   :: NTERM_BIG = 0
      INTEGER(LONG), ALLOCATABLE      :: J_BIG(:), I_BIG(:)! CCS of the large blocks: column starts, row numbers
      REAL(DOUBLE) , ALLOCATABLE      :: A_BIG(:)
      INTEGER(DBL_LONG)               :: BIG_FACTORS = 0   ! SuperLU handle

      REAL(DOUBLE) , ALLOCATABLE      :: W(:)              ! Work vector of length N, zero between solves
      INTEGER(LONG), ALLOCATABLE      :: MARK(:)           ! Last solve that touched each block
      INTEGER(LONG)                   :: NSOLVE = 0

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE BDS_FACTOR ( NA, NTERM, I_A, J_A, A, INFO, NUM_BLOCKS, LARGEST )

! Finds the blocks of A and factors them. INFO = 0, or the row of A of a zero pivot (A singular).

      USE SPARSE_FORMAT_CONVERSION, ONLY:  SPARSE_CRS_SPARSE_CCS

      INTEGER(LONG), INTENT(IN)       :: NA, NTERM, I_A(NA+1), J_A(NTERM)
      REAL(DOUBLE) , INTENT(IN)       :: A(NTERM)
      INTEGER(LONG), INTENT(OUT)      :: INFO, NUM_BLOCKS, LARGEST
      INTEGER(LONG), ALLOCATABLE      :: PARENT(:), CNT(:), ROOT_COMP(:), IBC(:), JBC(:)
      REAL(DOUBLE) , ALLOCATABLE      :: ABC(:)
      INTEGER(LONG)                   :: I, J, K, C, R1, R2, NC, NT, M
      INTEGER(DBL_LONG)               :: OFF
      INTEGER                         :: NI, INFO4, SYM0
      REAL(DOUBLE)                    :: DUM(1)

      CALL BDS_FREE
      N = NA
      INFO = 0
      ALLOCATE ( PARENT(N), COMP(N), LOC(N), W(N) )
      W = 0.0D0
      DO I=1,N                                             ! Union-find over the structure of A
         PARENT(I) = I
      ENDDO
      DO I=1,N
         DO K=I_A(I),I_A(I+1)-1
            R1 = FIND_ROOT ( PARENT, I )
            R2 = FIND_ROOT ( PARENT, J_A(K) )
            IF (R1 /= R2) PARENT(MAX(R1,R2)) = MIN(R1,R2)
         ENDDO
      ENDDO

      ALLOCATE ( ROOT_COMP(N) )                            ! Number the blocks in the order of their first row
      ROOT_COMP = 0
      NCOMP = 0
      DO I=1,N
         R1 = FIND_ROOT ( PARENT, I )
         IF (ROOT_COMP(R1) == 0) THEN
            NCOMP = NCOMP + 1
            ROOT_COMP(R1) = NCOMP
         ENDIF
         COMP(I) = ROOT_COMP(R1)
      ENDDO
      ALLOCATE ( CPTR(NCOMP+1), CDOF(N), CNT(NCOMP), MARK(NCOMP), LUOFF(NCOMP) )
      MARK = 0
      CNT  = 0
      DO I=1,N
         CNT(COMP(I)) = CNT(COMP(I)) + 1
      ENDDO
      CPTR(1) = 1
      DO C=1,NCOMP
         CPTR(C+1) = CPTR(C) + CNT(C)
      ENDDO
      CNT = 0
      DO I=1,N                                             ! Rows of each block in ascending order
         C = COMP(I)
         LOC(I) = CNT(C) + 1
         CDOF(CPTR(C) + CNT(C)) = I
         CNT(C) = CNT(C) + 1
      ENDDO

      LARGEST = 0                                          ! Dense storage for the small blocks
      OFF  = 0
      NBIG = 0
      DO C=1,NCOMP
         NC = CPTR(C+1) - CPTR(C)
         LARGEST = MAX(LARGEST, NC)
         IF (NC <= NDENSE_MAX) THEN
            LUOFF(C) = OFF
            OFF = OFF + INT(NC,DBL_LONG)*NC
         ELSE
            LUOFF(C) = -1
            NBIG  = NBIG + NC
            NBIGC = NBIGC + 1
         ENDIF
      ENDDO
      NUM_BLOCKS = NCOMP
      ALLOCATE ( LU(MAX(OFF,1_DBL_LONG)), PIV(N), BIGLOC(N), BIGROW(MAX(NBIG,1)) )
      LU = 0.0D0
      BIGLOC = 0
      NT = 0
      M  = 0
      DO C=1,NCOMP                                         ! Rows of the large blocks, numbered among themselves
         IF (LUOFF(C) >= 0) CYCLE
         DO K=CPTR(C),CPTR(C+1)-1
            M = M + 1
            BIGROW(M) = CDOF(K)
            BIGLOC(CDOF(K)) = M
         ENDDO
      ENDDO
      DO I=1,N                                             ! Terms of A into the dense blocks, or counted for the large ones
         C = COMP(I)
         NC = CPTR(C+1) - CPTR(C)
         DO K=I_A(I),I_A(I+1)-1
            IF (LUOFF(C) >= 0) THEN
               OFF = LUOFF(C) + INT(LOC(J_A(K))-1,DBL_LONG)*NC + LOC(I)
               LU(OFF) = LU(OFF) + A(K)
            ELSE
               NT = NT + 1
            ENDIF
         ENDDO
      ENDDO

      DO C=1,NCOMP                                         ! Dense LU of each small block
         IF (LUOFF(C) < 0) CYCLE
         NC = CPTR(C+1) - CPTR(C)
         NI = INT(NC)
         CALL DGETRF ( NI, NI, LU(LUOFF(C)+1), NI, PIV(CPTR(C)), INFO4 )
         IF (INFO4 /= 0) THEN
            INFO = CDOF(CPTR(C) + MAX(INFO4,1) - 1)
            RETURN
         ENDIF
      ENDDO

      IF (NBIG > 0) THEN                                   ! The large blocks together, with SuperLU (general matrix)
         NTERM_BIG = NT
         ALLOCATE ( IBC(NBIG+1), JBC(MAX(NT,1)), ABC(MAX(NT,1)), J_BIG(NBIG+1), I_BIG(MAX(NT,1)), A_BIG(MAX(NT,1)) )
         IBC(1) = 1
         NT = 0
         DO M=1,NBIG                                       ! CRS of the large blocks in their own numbering
            I = BIGROW(M)
            DO K=I_A(I),I_A(I+1)-1
               NT = NT + 1
               JBC(NT) = BIGLOC(J_A(K))
               ABC(NT) = A(K)
            ENDDO
            IBC(M+1) = NT + 1
         ENDDO
         CALL SPARSE_CRS_SPARSE_CCS ( NBIG, NBIG, NTERM_BIG, 'RMM BLOCKS', IBC, JBC, ABC, 'CCS', J_BIG, I_BIG, A_BIG, 'N' )
         DEALLOCATE ( IBC, JBC, ABC )
         SYM0 = 0
         CALL C_FORTRAN_DGSSV_SYMMETRIC ( SYM0 )
         DUM = 0.0D0
         CALL C_FORTRAN_DGSSV ( 1, NBIG, NTERM_BIG, 1, A_BIG, I_BIG, J_BIG, DUM, NBIG, BIG_FACTORS, INFO )
         SYM0 = 1
         CALL C_FORTRAN_DGSSV_SYMMETRIC ( SYM0 )
         IF (INFO /= 0) THEN
            IF ((INFO > 0) .AND. (INFO <= NBIG)) INFO = BIGROW(INFO)
            RETURN
         ENDIF
      ENDIF

      DEALLOCATE ( PARENT, CNT, ROOT_COMP )

      END SUBROUTINE BDS_FACTOR

! ##################################################################################################################################

      SUBROUTINE BDS_SOLVE ( NB, B_ROW, B_VAL, NX, X_ROW, X_VAL )

! x = A^-1 * b for b with NB nonzeros (rows B_ROW, values B_VAL). Returns x on all rows of the blocks b touches (NX values),
! in ascending row order. X_ROW and X_VAL must have room for N values.

      INTEGER(LONG), INTENT(IN)       :: NB, B_ROW(NB)
      REAL(DOUBLE) , INTENT(IN)       :: B_VAL(NB)
      INTEGER(LONG), INTENT(OUT)      :: NX, X_ROW(:)
      REAL(DOUBLE) , INTENT(OUT)      :: X_VAL(:)
      INTEGER(LONG)                   :: I, K, C, NC, NTOUCH, INFO
      INTEGER(LONG), ALLOCATABLE      :: TOUCH(:)
      INTEGER                         :: NI, INFO4
      LOGICAL                         :: ANY_BIG
      REAL(DOUBLE) , ALLOCATABLE      :: XB(:)

      NSOLVE = NSOLVE + 1
      ALLOCATE ( TOUCH(MAX(NB+NBIGC,1)) )
      NTOUCH  = 0
      ANY_BIG = .FALSE.
      DO K=1,NB                                            ! b into W; the blocks touched
         W(B_ROW(K)) = W(B_ROW(K)) + B_VAL(K)
         C = COMP(B_ROW(K))
         IF (MARK(C) /= NSOLVE) THEN
            MARK(C) = NSOLVE
            NTOUCH = NTOUCH + 1
            TOUCH(NTOUCH) = C
            IF (LUOFF(C) < 0) ANY_BIG = .TRUE.
         ENDIF
      ENDDO

      DO K=1,NTOUCH                                        ! Small blocks: dense solve in place on their rows
         C = TOUCH(K)
         IF (LUOFF(C) < 0) CYCLE
         NC = CPTR(C+1) - CPTR(C)
         NI = INT(NC)
         X_VAL(1:NC) = W(CDOF(CPTR(C):CPTR(C+1)-1))
         CALL DGETRS ( 'N', NI, 1, LU(LUOFF(C)+1), NI, PIV(CPTR(C)), X_VAL, NI, INFO4 )
         W(CDOF(CPTR(C):CPTR(C+1)-1)) = X_VAL(1:NC)
      ENDDO
      IF (ANY_BIG) THEN                                    ! Large blocks: one SuperLU solve over all of their rows
         ALLOCATE ( XB(NBIG) )
         XB = W(BIGROW(1:NBIG))
         INFO = 0
         CALL C_FORTRAN_DGSSV ( 2, NBIG, NTERM_BIG, 1, A_BIG, I_BIG, J_BIG, XB, NBIG, BIG_FACTORS, INFO )
         W(BIGROW(1:NBIG)) = XB
         DEALLOCATE ( XB )
      ENDIF

      IF (ANY_BIG) THEN                                    ! All large blocks were solved together: all of them are in x
         DO C=1,NCOMP
            IF ((LUOFF(C) < 0) .AND. (MARK(C) /= NSOLVE)) THEN
               MARK(C) = NSOLVE
               NTOUCH = NTOUCH + 1
               TOUCH(NTOUCH) = C
            ENDIF
         ENDDO
      ENDIF
      CALL SORT_BLOCKS ( NTOUCH, TOUCH )                   ! x on the touched blocks, rows ascending; W back to zero
      CALL MERGE_BLOCK_ROWS ( NTOUCH, TOUCH, NX, X_ROW, X_VAL )

      DEALLOCATE ( TOUCH )

      END SUBROUTINE BDS_SOLVE

! ##################################################################################################################################

      SUBROUTINE MERGE_BLOCK_ROWS ( NTOUCH, TOUCH, NX, X_ROW, X_VAL )

! Rows of the touched blocks in ascending order (each block's rows are ascending: a k-way merge by repeated selection of the
! smallest head, NTOUCH is small), with their values from W; W is set back to zero

      INTEGER(LONG), INTENT(IN)       :: NTOUCH, TOUCH(:)
      INTEGER(LONG), INTENT(OUT)      :: NX, X_ROW(:)
      REAL(DOUBLE) , INTENT(OUT)      :: X_VAL(:)
      INTEGER(LONG), ALLOCATABLE      :: HEAD(:)
      INTEGER(LONG)                   :: K, KMIN, RMIN, R

      ALLOCATE ( HEAD(MAX(NTOUCH,1)) )
      DO K=1,NTOUCH
         HEAD(K) = CPTR(TOUCH(K))
      ENDDO
      NX = 0
      DO
         KMIN = 0
         RMIN = HUGE(RMIN)
         DO K=1,NTOUCH
            IF (HEAD(K) < CPTR(TOUCH(K)+1)) THEN
               R = CDOF(HEAD(K))
               IF (R < RMIN) THEN
                  RMIN = R
                  KMIN = K
               ENDIF
            ENDIF
         ENDDO
         IF (KMIN == 0) EXIT
         NX = NX + 1
         X_ROW(NX) = RMIN
         X_VAL(NX) = W(RMIN)
         W(RMIN) = 0.0D0
         HEAD(KMIN) = HEAD(KMIN) + 1
      ENDDO
      DEALLOCATE ( HEAD )

      END SUBROUTINE MERGE_BLOCK_ROWS

! ##################################################################################################################################

      SUBROUTINE SORT_BLOCKS ( NT, T )

! Insertion sort of block numbers (blocks are numbered by their first row, so this orders them by first row)

      INTEGER(LONG), INTENT(IN)       :: NT
      INTEGER(LONG), INTENT(INOUT)    :: T(:)
      INTEGER(LONG)                   :: I, J, V

      DO I=2,NT
         V = T(I)
         J = I - 1
         DO WHILE (J >= 1)
            IF (T(J) <= V) EXIT
            T(J+1) = T(J)
            J = J - 1
         ENDDO
         T(J+1) = V
      ENDDO

      END SUBROUTINE SORT_BLOCKS

! ##################################################################################################################################

      SUBROUTINE BDS_FREE

      INTEGER(LONG)                   :: INFO
      REAL(DOUBLE)                    :: DUM(1)

      IF ((NBIG > 0) .AND. (BIG_FACTORS /= 0)) THEN
         DUM  = 0.0D0
         INFO = 0
         CALL C_FORTRAN_DGSSV ( 3, NBIG, NTERM_BIG, 1, A_BIG, I_BIG, J_BIG, DUM, NBIG, BIG_FACTORS, INFO )
      ENDIF
      BIG_FACTORS = 0
      NBIG  = 0
      NBIGC = 0
      NCOMP = 0
      N     = 0
      IF (ALLOCATED(COMP))   DEALLOCATE ( COMP )
      IF (ALLOCATED(LOC))    DEALLOCATE ( LOC )
      IF (ALLOCATED(CPTR))   DEALLOCATE ( CPTR )
      IF (ALLOCATED(CDOF))   DEALLOCATE ( CDOF )
      IF (ALLOCATED(LUOFF))  DEALLOCATE ( LUOFF )
      IF (ALLOCATED(LU))     DEALLOCATE ( LU )
      IF (ALLOCATED(PIV))    DEALLOCATE ( PIV )
      IF (ALLOCATED(BIGROW)) DEALLOCATE ( BIGROW )
      IF (ALLOCATED(BIGLOC)) DEALLOCATE ( BIGLOC )
      IF (ALLOCATED(J_BIG))  DEALLOCATE ( J_BIG )
      IF (ALLOCATED(I_BIG))  DEALLOCATE ( I_BIG )
      IF (ALLOCATED(A_BIG))  DEALLOCATE ( A_BIG )
      IF (ALLOCATED(W))      DEALLOCATE ( W )
      IF (ALLOCATED(MARK))   DEALLOCATE ( MARK )

      END SUBROUTINE BDS_FREE

! ##################################################################################################################################

      INTEGER(LONG) FUNCTION FIND_ROOT ( PARENT, I )

! Root of I in the union-find forest, with path halving

      INTEGER(LONG), INTENT(INOUT)    :: PARENT(:)
      INTEGER(LONG), INTENT(IN)       :: I
      INTEGER(LONG)                   :: R

      R = I
      DO WHILE (PARENT(R) /= R)
         PARENT(R) = PARENT(PARENT(R))
         R = PARENT(R)
      ENDDO
      FIND_ROOT = R

      END FUNCTION FIND_ROOT

      END MODULE BLOCK_DIAGONAL_SOLVE
