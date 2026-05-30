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

      SUBROUTINE SOLVE_GMN ( PART_VEC_G_NM, PART_VEC_M )

! Solves the sustem of equations: RMM*GMN = -RMN for matrix GMN which is used in the reduction of the G set stiffness, mass and
! load matrices from the G-set to the N, M_sets. If RMM is diagonal, a simple algorithm is used. If it is not, routines
! are called to do the decomp of RMM and the forward-backward substitution (FBS) to obtain GMN

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SCR, L2A, LINK2A, L2A_MSG, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NDOFM, NTERM_RMG, NTERM_RMN, NTERM_RMM, NTERM_GMN
      USE PARAMS, ONLY                :  EPSIL, PRTRMG, PRTGMN, SOLLIB, SPARSE_FLAVOR, SUPINFO
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE SPARSE_MATRICES, ONLY       :  I_RMG, J_RMG, RMG, I_RMN, J_RMN, RMN, I_RMM, J_RMM, RMM, I_GMN, J_GMN, GMN
      USE SPARSE_MATRICES, ONLY       :  SYM_RMG, SYM_RMN, SYM_RMM
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE SOLVE_GMN_USE_IFs

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_GMN'
      CHARACTER(  1*BYTE)             :: CLOSE_IT            ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER(  8*BYTE)             :: CLOSE_STAT          ! Char constant for the CLOSE status of a file
      CHARACTER(  1*BYTE)             :: RMM_DIAG            ! 'Y' if matrix RMM is diagonal.
      CHARACTER(  1*BYTE)             :: RMM_IDENTITY        ! 'Y' if matrix RMM is an identity matrix.

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_G_NM(NDOFG)! Partitioning vector (G set into N and M sets)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_M(NDOFM)   ! Partitioning vector (1's for all M set DOF's)
      INTEGER(LONG)                   :: I,J,K               ! DO loop indices
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices
      INTEGER(LONG)                   :: RMN_ROW_I_NTERMS    ! No. terms in row I of matrix RMN
      INTEGER(LONG)                   :: RMN_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: RMM_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)


      REAL(DOUBLE)                    :: EPS1                ! A small number to compare real zero

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Print out constraint matrix RMG if requested

      IF (( PRTRMG == 1) .OR. ( PRTRMG == 3)) THEN
         IF (NTERM_RMG > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CONSTRAINT MATRIX RMG', 'M ', 'G ', NTERM_RMG, NDOFM, I_RMG, J_RMG, RMG )
         ENDIF
      ENDIF

! Partition RMN from RMG. If no terms in RMN, write error and quit

      CALL PARTITION_SS_NTERM ( 'RMG', NTERM_RMG, NDOFM, NDOFG, SYM_RMG, I_RMG, J_RMG,      PART_VEC_M, PART_VEC_G_NM,             &
                                 NUM1, NUM1, RMN_ROW_MAX_TERMS,'RMN', NTERM_RMN, SYM_RMN )

      IF (NTERM_RMN > 0) THEN
         CALL ALLOCATE_SPARSE_MAT ( 'RMN', NDOFM, NTERM_RMN, SUBR_NAME )
         CALL PARTITION_SS ( 'RMG', NTERM_RMG, NDOFM, NDOFG, SYM_RMG, I_RMG, J_RMG, RMG, PART_VEC_M, PART_VEC_G_NM,                &
                              NUM1, NUM1, RMN_ROW_MAX_TERMS, 'RMN', NTERM_RMN, NDOFM, SYM_RMN, I_RMN, J_RMN, RMN )
      ELSE
         WRITE(ERR,2201) NTERM_RMN
         WRITE(F06,2201) NTERM_RMN
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Partition RMM from RMG. If no terms in RMM, write error and quit

      CALL PARTITION_SS_NTERM ( 'RMG', NTERM_RMG, NDOFM, NDOFG, SYM_RMG, I_RMG, J_RMG ,     PART_VEC_M, PART_VEC_G_NM,             &
                                 NUM1, NUM2, RMM_ROW_MAX_TERMS, 'RMM', NTERM_RMM, SYM_RMM )

      IF (NTERM_RMM > 0) THEN
         CALL ALLOCATE_SPARSE_MAT ( 'RMM', NDOFM, NTERM_RMM, SUBR_NAME )
         CALL PARTITION_SS ( 'RMG', NTERM_RMG, NDOFM, NDOFG, SYM_RMG, I_RMG, J_RMG, RMG, PART_VEC_M, PART_VEC_G_NM,                &
                              NUM1, NUM2, RMM_ROW_MAX_TERMS, 'RMM', NTERM_RMM, NDOFM, SYM_RMM, I_RMM, J_RMM, RMM )
      ELSE
         WRITE(ERR,2202) SUBR_NAME,NTERM_RMM
         WRITE(F06,2202) SUBR_NAME,NTERM_RMM
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Print out constraint matrix RMG partitions, if requested

      IF (( PRTRMG == 2) .OR. ( PRTRMG == 3)) THEN
         IF (NTERM_RMN > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CONSTRAINT PARTITION RMN', 'M ', 'N ', NTERM_RMN, NDOFM, I_RMN, J_RMN, RMN )
         ENDIF
         IF (NTERM_RMM > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CONSTRAINT PARTITION RMM', 'M ', 'M ', NTERM_RMM, NDOFM, I_RMM, J_RMM, RMM )
         ENDIF
      ENDIF

! Find out if RMM is a diagonal matrix. Getting sol'n for GMN will then be trivial

      RMM_DIAG     = 'Y'                                   ! Find out if RMM is a diagonal or identity matrix
      RMM_IDENTITY = 'Y'
      IF (NTERM_RMM == NDOFM) THEN                         ! There are as many terms in RMM as rows so maybe diag or identity
         DO I=1,NDOFM
            IF (J_RMM(I) /= I) THEN                        ! The i-th term in RMM is not a diagonal term
               RMM_DIAG     = 'N'
               RMM_IDENTITY = 'N'
               EXIT
            ENDIF
         ENDDO
         IF (RMM_DIAG == 'Y') THEN                         ! If RMM is diagonal, check for 1.0 on diagonal (identity matrix)
            DO I=1,NDOFM
               IF (DABS(RMM(I) - ONE) > EPS1) THEN
                  RMM_IDENTITY = 'N'
               ENDIF
            ENDDO
         ENDIF
      ELSE                                                 ! RMM is not identity or diagonal
         RMM_DIAG     = 'N'
         RMM_IDENTITY = 'N'
      ENDIF
! Now solve for GMN using either simple algorithm (RMM diagonal) or using BANDED or SPARSE equation solver

      IF ((RMM_DIAG == 'Y') .AND. (DEBUG(20) /= 1)) THEN  ! We can do simple inverse of diagonal matrix RMM

         WRITE(ERR,2293)
         IF (SUPINFO == 'N') THEN
            WRITE(F06,2293)
         ENDIF
         NTERM_GMN = NTERM_RMN

         CALL ALLOCATE_L2_GMN_2 ( SUBR_NAME )
         CALL ALLOCATE_SPARSE_MAT ( 'GMN', NDOFM, NTERM_GMN, SUBR_NAME )

         DO I=1,NDOFM+1
            I_GMN(I) = I_RMN(I)
         ENDDO
         K = 0
         DO I=1,NDOFM
            RMN_ROW_I_NTERMS = I_RMN(I+1) - I_RMN(I)
            DO J=1,RMN_ROW_I_NTERMS
               K = K + 1
               J_GMN(K) = J_RMN(K)
               IF (RMM_IDENTITY == 'Y') THEN
                  GMN(K) = -RMN(K)
               ELSE
                  IF (DABS(RMM(I)) > EPS1) THEN
                     J_GMN(K)  = J_RMN(K)
                       GMN(K)  = -RMN(K)/RMM(I)
                  ELSE
                     WRITE(ERR,2203) K
                     WRITE(F06,2203) K
                     FATAL_ERR = FATAL_ERR + 1
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ELSE                                                 ! Either RMM is not diagonal or DEBUG(20) =1 so we will do full sol'n
!                                                            for GMN from eqn RMM*GMN = -RMN
         IF (RMM_DIAG == 'Y') THEN
            WRITE(ERR,2294)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2294)
            ENDIF
         ENDIF
         CALL SOLVE_GMN_SOLVER

      ENDIF

! Check to maks sure GMN has nonzero terma

      IF (NTERM_GMN <= 0 )THEN
         WRITE(ERR,2204) SUBR_NAME,NTERM_GMN
         WRITE(F06,2204) SUBR_NAME,NTERM_GMN
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error (NTERM_GMN <= 0 ), so quit
      ENDIF

! Close GMN

      IF (NTERM_GMN > 0) THEN
         CLOSE_IT   = 'Y'
         CLOSE_STAT = 'KEEP'
         CALL WRITE_MATRIX_1 ( LINK2A, L2A, CLOSE_IT, CLOSE_STAT, L2A_MSG, 'GMN', NTERM_GMN, NDOFM, I_GMN, J_GMN, GMN )
      ENDIF

! Print out constraint matrix GMN, if requested

      IF ( PRTGMN == 1) THEN
         IF (NTERM_GMN > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CONSTRAINT MATRIX GMN', 'M ', 'N ', NTERM_GMN, NDOFM, I_GMN, J_GMN, GMN )
         ENDIF
      ENDIF

! Deallocate partitions of RMG: RMN, RMM. Keep GMN, it is needed in the reduction of KGG, MGG and PG

      WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RMN', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'RMN' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RMM', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'RMM' )

      CALL DEALLOCATE_L2_GMN_2



      RETURN

! **********************************************************************************************************************************
 2201 FORMAT(' *ERROR  2201: THE RMN PARTITION OF CONSTRAINT MATRIX RMG HAS ',I12,' TERMS IN IT.'                                  &
                    ,/,14X,' THE NUMBER MUST BE > 0.'                                                                              &
                    ,/,14X,' THIS PROBABLY MEANS THAT THERE WERE NO N-SET DOFs THAT THE RIGID ELEMS & MPCs WERE DEPENDENT ON')

 2202 FORMAT(' *ERROR  2202: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE RMM PARTITION OF CONSTRAINT MATRIX RMG HAS ',I12,' TERMS IN IT. THE NUMBER MUST BE > 0.')

 2203 FORMAT(' *ERROR  2203: ZERO DIAGONAL IN MATRIX RMM FOR M-SET DOF ',I8)

 2204 FORMAT(' *ERROR  2204: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THERE IS AN M-SET BUT CONSTRAINT MATRIX GMN HAS ',I12,' TERMS IN IT. MUST BE > 0')

 2293 FORMAT(' *INFORMATION: THE RMM CONSTRAINT MATRIX IS DIAGONAL.'                                                               &
                    ,/,14X,' A SIMPLE SOLUTION FOR THE GMN CONSTRAINT MATRIX WILL BE USED AVOIDING CALLING SUBR SOLVE_GMN',/)

 2294 FORMAT(' *INFORMATION: THE RMM CONSTRAINT MATRIX IS DIAGONAL. HOWEVER, SINCE DEBUG(20) = 1'                                  &
                    ,/,14X,' SUBR SOLVE_GMN_SOLVER WILL BE CALLED TO SOLVE FOR THE GMN CONSTRAINT MATRIX',/)

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(A,10X)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE SOLVE_GMN_SOLVER

! Solves RMM x GMN = -RMN for matrix GMN using unsymmetric decomp from LAPACK

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  NDOFG, NDOFM, NDOFN, NTERM_GMN, NTERM_RMM, NTERM_RMN, BLNK_SUB_NAM
      USE PARAMS, ONLY                :  EPSIL, SOLLIB, SPARSE_FLAVOR
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE SPARSE_MATRICES, ONLY       :  I_RMN, J_RMN, RMN, I_RMM, J_RMM, RMM, I2_GMN, I_GMN, J_GMN, GMN
      USE SCRATCH_MATRICES, ONLY      :  I_CCS1, J_CCS1, CCS1
      USE FULL_MATRICES, ONLY         :  RMM_FULL
      USE LAPACK_LIN_EQN_DGE
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, SLU_INFO

! Interface module not needed for subr's DGETRF and DGETRS. These are "CONTAIN'ed" in module LAPACK_LIN_EQN_DPB, which
! is "USE'd" above

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_GMN_SOLVER'
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: CALLED_SUBR = ' ' ! Name of a called subr (for output error purposes)
      CHARACTER( 1*BYTE)              :: CLOSE_IT          ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER( 8*BYTE)              :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER(24*BYTE)              :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN
      CHARACTER(44*BYTE)              :: MODNAM            ! Name to write to screen to describe module being run
      CHARACTER(24*BYTE)              :: MODNAM1           ! Name to write to screen to describe module being run
      CHARACTER( 1*BYTE)              :: READ_NTERM        ! 'Y' or 'N' Input to subr READ_MATRIX_1
      CHARACTER( 1*BYTE)              :: NULL_COL          ! 'Y' if a col of RMN is null
      CHARACTER( 1*BYTE)              :: OPND              ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SCRFIL            ! File name
      CHARACTER( 1*BYTE)              :: TRANS             ! 'Y' if
      LOGICAL                         :: RMM_SOLVED_WITH_SUPERLU
      LOGICAL                         :: RMM_FULL_ALLOCATED
      LOGICAL                         :: CCS1_ALLOCATED

      INTEGER(LONG)                   :: COMPV             ! Component number (1-6) of a grid DOF
      INTEGER(LONG)                   :: GRIDV             ! Grid number
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: INFO      = 0     ! Output from factorization routines
      INTEGER(LONG)                   :: IPIV(NDOFM)       ! Pivot indices from factorization of RMM
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: NRHS              ! No. of RHS's in solving (RMM)*(GMN) = -RMN
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN


      REAL(DOUBLE)                    :: BETA              ! Multiple for rhs for use in subr FBS
      REAL(DOUBLE)                    :: DUM_COL(NDOFM)    ! Temp variable used in SuperLU
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: GMN_COL(NDOFM)    ! A column of GMN solved for herein
      REAL(DOUBLE)                    :: RMN_COL(NDOFM)    ! A column of RMN. The solution for GMN_COL is from RMM*GMN_COL = RMN_COL

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06
      RMM_SOLVED_WITH_SUPERLU = .FALSE.
      RMM_FULL_ALLOCATED      = .FALSE.
      CCS1_ALLOCATED          = .FALSE.

      IF (SOLLIB == 'BANDED  ') THEN
                                                            ! Create full matrix RMM_FULL from sparse RMM
         CALL ALLOCATE_FULL_MAT  ( 'RMM_FULL', NDOFM, NDOFM, SUBR_NAME )
         RMM_FULL_ALLOCATED = .TRUE.
         CALL SPARSE_CRS_TO_FULL ( 'RMM       ', NTERM_RMM, NDOFM, NDOFM, SYM_RMM, I_RMM, J_RMM, RMM, RMM_FULL )

                                                           ! Perform factorization of RMM_FULL matrix.

         CALL OURTIM
         MODNAM = '    Lapack factorization of RMM'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
         CALL DGETRF (NDOFM, NDOFM, RMM_FULL, NDOFM, IPIV, INFO )

         CALLED_SUBR = 'DGETRF'
         IF      (INFO < 0) THEN                           ! LAPACK subr XERBLA should have reported error on an illegal argument
!                                                            in a called LAPACK subr, so we should not have gotten here
            WRITE(ERR,993) SUBR_NAME, CALLED_SUBR
            WRITE(F06,993) SUBR_NAME, CALLED_SUBR
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                         ! Coding error, so quit

         ELSE IF (INFO > 0) THEN                           ! 0 diag in RMM

             CALL GET_GRID_AND_COMP ( 'M ', INFO, GRIDV, COMPV  )

            IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN
               WRITE(ERR,2501) CALLED_SUBR, SUBR_NAME, INFO
               WRITE(F06,2501) CALLED_SUBR, SUBR_NAME, INFO
               IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
                  WRITE(ERR,25012) GRIDV, COMPV
                  WRITE(F06,25012) GRIDV, COMPV
               ENDIF
               WRITE(ERR,25013)
               WRITE(F06,25013)
               SLU_INFO = 0
               CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFM, NTERM_RMM, SUBR_NAME )
               CCS1_ALLOCATED = .TRUE.
               CALL SPARSE_CRS_SPARSE_CCS ( NDOFM, NDOFM, NTERM_RMM, 'RMM', I_RMM, J_RMM, RMM, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y')
               CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'RMM', 'M ', NDOFM, NTERM_RMM, J_CCS1, I_CCS1, CCS1, SLU_INFO )
               RMM_SOLVED_WITH_SUPERLU = .TRUE.
            ELSE
               WRITE(ERR,2501) CALLED_SUBR, SUBR_NAME, INFO
               WRITE(F06,2501) CALLED_SUBR, SUBR_NAME, INFO
               FATAL_ERR = FATAL_ERR + 1
               IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
                  WRITE(ERR,25012) GRIDV, COMPV
                  WRITE(F06,25012) GRIDV, COMPV
               ENDIF
               CALL OUTA_HERE ( 'Y' )
            ENDIF

         ENDIF

      ELSE IF (SOLLIB == 'SPARSE  ') THEN

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            SLU_INFO = 0
            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFM, NTERM_RMM, SUBR_NAME )
            CCS1_ALLOCATED = .TRUE.
            CALL SPARSE_CRS_SPARSE_CCS ( NDOFM, NDOFM, NTERM_RMM, 'RMM', I_RMM, J_RMM, RMM, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y')
            CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'RMM', 'M ', NDOFM, NTERM_RMM, J_CCS1, I_CCS1, CCS1, SLU_INFO )
            RMM_SOLVED_WITH_SUPERLU = .TRUE.

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

! **********************************************************************************************************************************
! Open a scratch file that will be used to write GMN nonzero terms to as we solve for columns of GMN. After all col's
! of GMN have been solved for, and we have a count on NTERM_GMN, we will allocate memory to the GMN arrays and read
! the scratch file values into those arrays. Then, in the calling subroutine, we will write NTERM_GMN, followed by
! GMN row/col/value to a permanent file

      SCRFIL(1:)  = ' '
      SCRFIL(1:9) = 'SCRATCH-991'
      OPEN (SCR(1),STATUS='SCRATCH',FORM='UNFORMATTED',ACTION='READWRITE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SCRFIL, OUNT )
         CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )
         CALL OUTA_HERE ( 'Y' )                            ! Can't open scratch file, so quit
      ENDIF
      REWIND (SCR(1))

! Loop on columns of RMN

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      NTERM_GMN = 0
      CALL COUNTER_INIT('      Solve for GMN col ', NDOFN)
      DO J = 1,NDOFN

         !CALL OURTIM
         !MODNAM1 = '      Solve for GMN col '



! Set RMN_COL to the negative of i-th col of array RMN. First, initialize RMN_COL to zero
! Keep track of whether this col is null, so we can avoid FBS if it is.

         NULL_COL = 'Y'
         DO I=1,NDOFM
            RMN_COL(I) = ZERO
            gmn_col(i) = zero
         ENDDO

         BETA = -ONE
         CALL GET_SPARSE_CRS_COL ( 'RMN       ',J, NTERM_RMN, NDOFM, NDOFN, I_RMN, J_RMN, RMN, BETA, RMN_COL, NULL_COL )

! Calculate GMN_COL via forward/backward substitution. Remember that rhs is -RMN.

         IF (NULL_COL == 'N') THEN                         ! DGETRS will solve for GMN_COL & load it into GMN array

            IF (RMM_SOLVED_WITH_SUPERLU) THEN
               SLU_INFO = 0
               CALL FBS_SUPRLU ( SUBR_NAME, 'RMM', NDOFM, NTERM_RMM, J_CCS1, I_CCS1, CCS1, J, RMN_COL, SLU_INFO )

            ELSE IF (SOLLIB == 'BANDED  ') THEN
               TRANS = 'N'
               NRHS = 1
               CALL DGETRS ( TRANS, NDOFM, NRHS ,RMM_FULL, NDOFM, IPIV, RMN_COL, NDOFM, INFO )

               CALLED_SUBR = 'DGETRS'
               IF      (INFO < 0) THEN                     ! LAPACK subr XERBLA should have reported error on an illegal argument
!                                                            in calling a LAPACK subr, so we should not have gotten here
                  WRITE(ERR,993) SUBR_NAME, CALLED_SUBR
                  WRITE(F06,993) SUBR_NAME, CALLED_SUBR
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )

               ENDIF

            ELSE IF (SOLLIB == 'SPARSE  ') THEN
               IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

                  SLU_INFO = 0
                  CALL FBS_SUPRLU ( SUBR_NAME, 'RMM', NDOFM, NTERM_RMM, J_CCS1, I_CCS1, CCS1, J, RMN_COL, SLU_INFO )

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

            DO I=1,NDOFM
               GMN_COL(I) = RMN_COL(I)
            ENDDO

            DO I=1,NDOFM                                   ! Count NTERM_GMN and write nonzero GMN to scratch file
               IF (DABS(GMN_COL(I)) > EPS1) THEN
                  NTERM_GMN = NTERM_GMN + 1
                  WRITE(SCR(1)) I,J,GMN_COL(I)
               ENDIF
            ENDDO
         ENDIF
         CALL COUNTER_PROGRESS(J)
      ENDDO

      WRITE(SC1,*) CR13

      IF (CCS1_ALLOCATED) THEN
         CALL DEALLOCATE_SCR_MAT ( 'CCS1' )
      ENDIF
      IF (RMM_FULL_ALLOCATED) THEN
         CALL DEALLOCATE_FULL_MAT ( 'RMM_FULL' )
      ENDIF

 FreeS:IF (RMM_SOLVED_WITH_SUPERLU) THEN                    ! Last, free the storage allocated inside SuperLU

          IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            DO J=1,NDOFM                                         ! Need a null col of loads when SuperLU is called to factor KLL
               DUM_COL(J) = ZERO                                  ! (only because it appears in the calling list)
            ENDDO

            CALL C_FORTRAN_DGSSV( 3, NDOFM, NTERM_RMM, 1, RMM, I_RMM, J_RMM, DUM_COL, NDOFM, SLU_FACTORS, SLU_INFO )

            IF (SLU_INFO .EQ. 0) THEN
               WRITE (*,*) 'SUPERLU STORAGE FREED'
            ELSE
               WRITE(*,*) 'SUPERLU STORAGE NOT FREED. INFO FROM SUPERLU FREE STORAGE ROUTINE = ', SLU_INFO
            ENDIF

         ENDIF

      ENDIF FreeS

! The GMN data in SCRATCH-991 is written 1 col at a time. We need it to be written for 1 row at a time with rows in numerical order

      CALL ALLOCATE_L2_GMN_2 ( SUBR_NAME )
      CALL ALLOCATE_SPARSE_MAT ( 'GMN', NDOFM, NTERM_GMN, SUBR_NAME )
      REWIND (SCR(1))
      MESSAG = 'SCRATCH: GMN ROW/COL/VAL'
      READ_NTERM = 'N'
      OPND       = 'Y'
      CLOSE_IT   = 'N'
      CLOSE_STAT = 'KEEP    '
      CALL READ_MATRIX_2 (SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG,'GMN',NDOFM, NTERM_GMN, READ_NTERM, I2_GMN, J_GMN, GMN)
      CALL SORT_INT2_REAL1 ( SUBR_NAME, 'I2_GMN, J_GMN, GMN', NTERM_GMN, I2_GMN, J_GMN, GMN )
      REWIND (SCR(1))
      WRITE(SCR(1)) NTERM_GMN
      DO K=1,NTERM_GMN
         WRITE(SCR(1)) I2_GMN(K),J_GMN(K),GMN(K)
      ENDDO

! Reallocate memory to GMN based on the NTERM_GMN counted above and read values from scratch file into GMN arrays

      CALL DEALLOCATE_L2_GMN_2

      WRITE(SC1, * ) '    Reallocate GMN'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GMN', CR13   ;  CALL DEALLOCATE_SPARSE_MAT ( 'GMN' )
      WRITE(SC1,12345,ADVANCE='NO') '       Allocate   GMN', CR13   ;  CALL ALLOCATE_SPARSE_MAT ('GMN', NDOFM, NTERM_GMN, SUBR_NAME)

      REWIND (SCR(1))
      MESSAG = 'SCRATCH: GMN ROW/COL/VAL'
      READ_NTERM = 'Y'
      OPND       = 'Y'
      CLOSE_IT   = 'Y'
      CLOSE_STAT = 'DELETE  '
      CALL READ_MATRIX_1 ( SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG, 'GMN', NTERM_GMN, READ_NTERM, NDOFM,                &
                           I_GMN, J_GMN, GMN)

      CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )



      RETURN

! **********************************************************************************************************************************
 2501 FORMAT(' *ERROR  2501: LAPACK SUBROUTINE, ',A8,' CALLED BY SUBROUTINE ',A                                                    &
                    ,/,14X,' HAS DETECTED A ZERO ON THE DIAG IN ROW ',I12,' OF THE TRIANG FACTOR IN THE DECOMP OF MATRIX RMM.')

25012 FORMAT('               THIS CORRESPONDS TO THE ROW & COL IN RMM FOR GRID POINT ',I8,' COMPONENT ',I3,'.'                     &
                    ,/,14X,' TO CORRECT THIS SITUATION, REMOVE THAT COMPONENT FROM REFC IN FIELD 5 OF THE OFFENDING RBE3(s)')

25013 FORMAT(' *WARNING    : LAPACK RMM FACTORIZATION FAILED. MYSTRAN WILL FALL BACK TO SUPERLU FOR THE GMN CONSTRAINT SOLVE.')

 2092 FORMAT(4X,A44,20X,I2,':',I2,':',I2,'.',I3)

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(A,10X,A)






  993 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' LAPACK SUBR XERBLA SHOULD HAVE REPORTED AN ERROR ON AN ILLEGAL ARGUMENT IN A CALL TO LAPACK SUBR '    &
                    ,/,15X,A,' (OR A SUBR CALLED BY IT) AND THEN ABORTED')

! **********************************************************************************************************************************

      END SUBROUTINE SOLVE_GMN_SOLVER

      END SUBROUTINE SOLVE_GMN
