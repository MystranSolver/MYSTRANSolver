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

   MODULE LAPACK_ADAPTERS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BANDGEN_LAPACK_DGB, BANDGEN_LAPACK_DPB, BANDSIZ, FBS_LAPACK, SYM_MAT_DECOMP_LAPACK

   CONTAINS

       SUBROUTINE BANDGEN_LAPACK_DGB ( MATIN_NAME, N, KD, NTERM_MATIN, I_MATIN, J_MATIN, MATIN, MATOUT, CALLING_SUBR )

! Puts sparse matrix MATIN into ARPACK banded form in matrix MATOUT. MATIN is symmetric and has KD super diagonals (determined in
! subr BANDSIZ). The terms from MATIN (all terms, not just the upper triangle) are stored in rows KD+1 through 3KD+1. The first KD
! rows of MATOUT are not used in storing MATIN terms (they must be needed in the ARPACK algorithm for other purposes?)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  SPARSTOR

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BANDGEN_LAPACK_DGB'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR         ! Name of subr calling this one
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME           ! Name of matrix input

      INTEGER(LONG)                   :: ROW_NO,COL_NO        ! Row and col no's of terms in banded matrix MATOUT

      INTEGER(LONG), INTENT(IN)       :: N                    ! Number of cols (or rows) of symmetric matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: NTERM_MATIN          ! No. of terms in sparse matrix
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(N+1)         ! Array of row no's for terms in matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERM_MATIN) ! Array of col no's for terms in matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: KD                   ! Number of sub (or super) diagonals in matrix MATIN.
                                                              ! The total band width of MATIN is 2*KD + 1. However, the upper
                                                              ! or lower triangle can be stored in an array of size KD+1 x N
      INTEGER(LONG)                   :: I,J                  ! DO loop indices
      INTEGER(LONG)                   :: K                    ! Counter
      INTEGER(LONG)                   :: MATOUT_DIAG_ROW_NUM  ! Number of the row, in mATOUT, where the diagonal of MATIN goes
      INTEGER(LONG)                   :: NUM_TERMS_ROW_I      ! Number of terms in MATIN matrix in row I


      REAL(DOUBLE) , INTENT(IN)       :: MATIN(NTERM_MATIN)   ! Array of terms in sparse matrix MATIN
      REAL(DOUBLE) , INTENT(INOUT)    :: MATOUT(3*KD+1,N)     ! Array of terms in band matrix MATOUT



! **********************************************************************************************************************************
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

! First, put terms from upper triangle of MATIN (which has only the upper triangle stored) into MATOUT. Note that this storage
! scheme (required by ARPACK Lanczos subr's) has the first KD rows unused in MATOUT.

      K  = 0
      CALL COUNTER_INIT("       Orig matrix row", N)
      DO I=1,N
         NUM_TERMS_ROW_I = I_MATIN(I+1) - I_MATIN(I)       ! Number of terms in row I
         DO J=1,NUM_TERMS_ROW_I
            K = K + 1
            IF (K > NTERM_MATIN) THEN
               WRITE(ERR,923) SUBR_NAME,K,NTERM_MATIN
               WRITE(F06,923) SUBR_NAME,K,NTERM_MATIN
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Coding error (attempt to exceed allocated array size), so quit
            ENDIF
            COL_NO = J_MATIN(K)
            ROW_NO = 2*KD + 1 + I - J_MATIN(K)
            IF ((ROW_NO > 0) .AND. (ROW_NO <= 3*KD+1) .AND. (COL_NO > 0) .AND. (COL_NO <= N)) THEN
               MATOUT(ROW_NO,COL_NO) = MATIN(K)
            ELSE
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,955) SUBR_NAME, MATIN_NAME, CALLING_SUBR, KD+1, ROW_NO, N, COL_NO
               WRITE(F06,955) SUBR_NAME, MATIN_NAME, CALLING_SUBR, N   , ROW_NO, N, COL_NO
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            MATOUT(ROW_NO,COL_NO) = MATIN(K)
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO
      WRITE(SC1,*) CR13

! Second, symmetric terms from MATIN into MATOUT

      IF (SPARSTOR == 'SYM') THEN
   !xx   WRITE(SC1, * )
         MATOUT_DIAG_ROW_NUM = 2*KD + 1
         DO I=1,KD
            CALL COUNTER_INIT("       Band matrix row", N-I)
            DO J=1,N-I
               MATOUT(MATOUT_DIAG_ROW_NUM+I,J) = MATOUT(MATOUT_DIAG_ROW_NUM-I,J+I)
               CALL COUNTER_PROGRESS(J)
            ENDDO
         ENDDO
         WRITE(SC1,*) CR13
      ENDIF



      RETURN

! **********************************************************************************************************************************

  923 FORMAT(' *ERROR   923: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INDEX K = ',I12,' IS GREATER THAN NTERM_MATIN = ',I12)

  955 FORMAT(' *ERROR   955: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INVALID ROW OR COL NUMBER CALCULATED FOR BAND FORM OF MATRIX ',A,' FROM CALLING SUBR ',A              &
                    ,/,14X,' ROW NUMBERS MUST BE > 0 AND <= ',I8,'. ROW NUMBER CALCULATED WAS ',I8                                 &
                    ,/,14X,' COL NUMBERS MUST BE > 0 AND <= ',I8,'. COL NUMBER CALCULATED WAS ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE BANDGEN_LAPACK_DGB


       SUBROUTINE BANDGEN_LAPACK_DPB ( MATIN_NAME, N, KD, NTERM_MATIN, I_MATIN, J_MATIN, MATIN, MATOUT, CALLING_SUBR )

! Puts sparse matrix MATIN into LAPACK DPB banded form in matrix MATOUT. MATIN is symmetric and has KD super diagonals
! (determined in subr BANDSIZ) and can be stored (upper triangle) in array MATOUT with KD+1 rows and N cols

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  SPARSTOR

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BANDGEN_LAPACK_DPB'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR         ! Name of subr calling this one
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME           ! Name of matrix input

      INTEGER(LONG)                   :: ROW_NO,COL_NO        ! Row and col no's of terms in banded matrix MATOUT

      INTEGER(LONG), INTENT(IN)       :: N                    ! Number of cols (or rows) of symmetric matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: NTERM_MATIN          ! No. of terms in sparse matrix
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(N+1)         ! Array of row no's for terms in matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERM_MATIN) ! Array of col no's for terms in matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: KD                   ! Number of sub (or super) diagonals in matrix MATIN.
                                                              ! The total band width of MATIN is 2*KD + 1. However, the upper
                                                              ! or lower triangle can be stored in an array of size KD+1 x N
      INTEGER(LONG)                   :: I,J                  ! DO loop indices
      INTEGER(LONG)                   :: K                    ! Counter
      INTEGER(LONG)                   :: NUM_TERMS_ROW_I      ! Number of terms in MATIN matrix in row I


      REAL(DOUBLE) , INTENT(IN)       :: MATIN(NTERM_MATIN)   ! Array of terms in sparse matrix MATIN
      REAL(DOUBLE) , INTENT(INOUT)    :: MATOUT(KD+1,N)       ! Array of terms in band matrix MATOUT



! **********************************************************************************************************************************
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      ROW_NO = 0
      K  = 0
      CALL COUNTER_INIT("       Orig matrix row", N)
      DO I=1,N
         NUM_TERMS_ROW_I = I_MATIN(I+1) - I_MATIN(I)       ! Number of terms in row I
         DO J=1,NUM_TERMS_ROW_I
            K = K + 1
            IF (K > NTERM_MATIN) THEN
               WRITE(ERR,923) SUBR_NAME, K, NTERM_MATIN
               WRITE(F06,923) SUBR_NAME, K, NTERM_MATIN
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Coding error (attempt to exceed allocated array size), so quit
            ENDIF
            COL_NO = J_MATIN(K)
            IF (J_MATIN(K) >= I) THEN                      ! New test that should make this subr work whether MATIN is stored sym
               ROW_NO = (KD + 1) - (J_MATIN(K) - I)        ! (nonzeros on/above diag) or nonsym (all nonzero terms in each row)
               IF ((ROW_NO > 0) .AND. (ROW_NO <= KD+1) .AND. (COL_NO > 0) .AND. (COL_NO <= N)) THEN
                  MATOUT(ROW_NO,COL_NO) = MATIN(K)
               ELSE
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,955) SUBR_NAME, MATIN_NAME, CALLING_SUBR, KD+1, ROW_NO, N, COL_NO
                  WRITE(F06,955) SUBR_NAME, MATIN_NAME, CALLING_SUBR, N   , ROW_NO, N, COL_NO
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ENDIF
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO
      WRITE(SC1,*) CR13



      RETURN

! **********************************************************************************************************************************
  923 FORMAT(' *ERROR   923: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INDEX K = ',I12,' IS GREATER THAN NTERM_MATIN = ',I12)

  955 FORMAT(' *ERROR   955: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INVALID ROW OR COL NUMBER CALCULATED FOR BAND FORM OF MATRIX ',A,' FROM CALLING SUBR ',A              &
                    ,/,14X,' ROW NUMBERS MUST BE > 0 AND <= ',I8,'. ROW NUMBER CALCULATED WAS ',I8                                 &
                    ,/,14X,' COL NUMBERS MUST BE > 0 AND <= ',I8,'. COL NUMBER CALCULATED WAS ',I8)

12345 format(7X,'Row ',I8,' of ',I8,' of orig matrix with ',I8,' term(s)',10X,A)

98788 format('       I       J       K  J_MATIN(K)      KD  ROW_NO  COL_NO    MATIN(K)'                                            &
          ,/,'       -       -       -  ----------      --  ------  ------   --------')

98789 format(3i8,i12,3i8,1es14.6)

! **********************************************************************************************************************************

      END SUBROUTINE BANDGEN_LAPACK_DPB


       SUBROUTINE BANDSIZ ( N, NTERM_MATIN, I_MATIN, J_MATIN, KD )

! Determines the band size of a matrix so that it can be put into the banded form required by the
! LAPACK routines that will be used to decompose it.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BANDSIZ'

      INTEGER(LONG), INTENT(IN)       :: N                    ! Col or row size of matrix MATIN (no. of A-set DOF's)
      INTEGER(LONG), INTENT(IN)       :: NTERM_MATIN          ! No. of terms in sparse matrix
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(N+1)         ! Array of row no's for terms in matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERM_MATIN) ! Array of col no's for terms in matrix MATIN
      INTEGER(LONG), INTENT(OUT)      :: KD                   ! Number of sub (or super) diagonals in matrix MATIN.
                                                              ! The total band width of MATIN is 2*KD + 1. However, the upper
                                                              ! or lower triangle can be stored in an array of size KD+1 x N
      INTEGER(LONG)                   :: I,J                  ! DO loop index
      INTEGER(LONG)                   :: K                    ! Counter
      INTEGER(LONG)                   :: KD_TEMP              ! Temporary value of in calculation of KD
      INTEGER(LONG)                   :: NUM_TERMS_ROW_I      ! Number of terms in MATIN matrix in row I




! **********************************************************************************************************************************
! Initialize outputs

      KD = 0

      K  = 0
      DO I=1,N
         NUM_TERMS_ROW_I = I_MATIN(I+1) - I_MATIN(I)  ! Number of terms in row I
         DO J=1,NUM_TERMS_ROW_I
            K = K + 1
            IF (K > NTERM_MATIN) THEN
               WRITE(ERR,927) SUBR_NAME,K,NTERM_MATIN
               WRITE(F06,927) SUBR_NAME,K,NTERM_MATIN
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            KD_TEMP = J_MATIN(K) - I
            IF (KD_TEMP >  KD) THEN
               KD = KD_TEMP
            ENDIF
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************
  927 FORMAT(' *ERROR   927: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INDEX K = ',I12,' IS GREATER THAN NTERM_MATIN = ',I12)




! **********************************************************************************************************************************

      END SUBROUTINE BANDSIZ


      SUBROUTINE FBS_LAPACK ( EQUED, NROWS, MATIN_SDIA, EQUIL_SCALE_FACS, INOUT_COL )

! FBS_LAPACK performs the forward-backward substitution to get displacements after the stiffness matrix has been decomposed using
! subr SYM_MAT_DECOMP_LAPACK

! (1) Scales the INOUT (e.g. load) vector if requested
! (2) Does the forward-backward solution to solve for a left-hand side vector given 1 right-hand side input vector (INOUT_COL)
!     The reult is returned in INOUT_COL
! (3) Scales the solution vector if requested

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, LINKNO
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, STIME, TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSERR, RCONDK
      USE LAPACK_DPB_MATRICES, ONLY   :  ABAND, LAPACK_S, RES
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG, NDEBUG
      USE MACHINE_PARAMS, ONLY        :  MACH_EPS, MACH_SFMIN

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE


      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'FBS_LAPACK'
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: CALLED_SUBR = ' ' ! Name of a called subr (for output error purposes)
      CHARACTER(  1*BYTE), INTENT(IN) :: EQUED             ! 'Y' if MATIN was equilibrated in subr EQUILIBRATE (called herein)
                                                           !       and the factorization (in DPBTRS) could not be completed.
      CHARACTER(  1*BYTE), PARAMETER  :: UPLO        = 'U' ! Indicates upper triang part of matrix is stored


      INTEGER(LONG), INTENT(IN)       :: MATIN_SDIA        ! No. of superdiags in the MATIN upper triangle
      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in sparse matrix MATIN
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: INFO        = 0   ! Output from LAPACK routine to do factorization of ABAND
                                                           !   0:  successful exit
                                                           ! < 0:  if INFO = -i, the i-th argument had an illegal value
                                                           ! > 0:  if INFO = i, the leading minor of order i is not pos def
                                                           !       and the factorization (in DPBTRS) could not be completed.
      INTEGER(LONG), PARAMETER        :: NUM_COLS    = 1   ! Number of vectors to solve in this call


      REAL(DOUBLE) , INTENT(IN)       :: EQUIL_SCALE_FACS(NROWS)
                                                           ! LAPACK_S values to return to calling subr

      REAL(DOUBLE) , INTENT(INOUT)    :: INOUT_COL(NROWS)    ! INOUT input  vector



! **********************************************************************************************************************************
! Scale the INOUT vector if requested

      IF (EQUED == 'Y') THEN
         DO I=1,NROWS
            INOUT_COL(I) = EQUIL_SCALE_FACS(I)*INOUT_COL(I)
         ENDDO
      ENDIF

! Calculate the answer via forward/backward substitution. Subr DPBTRS returns the answer in the workspace (INOUT_COL)
! for the right-hand side (which, at entry, was INOUT_COL)

      CALL DPBTRS ( UPLO, NROWS, MATIN_SDIA, NUM_COLS, ABAND, MATIN_SDIA+1, INOUT_COL, NROWS, INFO )

      CALLED_SUBR = 'DPBTRS'
      IF (INFO < 0) THEN                                   ! LAPACK subr XERBLA should have reported error on an illegal argument
         WRITE(ERR,993) SUBR_NAME, CALLED_SUBR             ! in calling a LAPACK subr, so we should not have gotten here
         WRITE(F06,993) SUBR_NAME, CALLED_SUBR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Scale the solution vector if requested (now overwritten in INOUT_COL)

      IF (EQUED == 'Y') THEN
         DO I=1,NROWS
            INOUT_COL(I) = EQUIL_SCALE_FACS(I)*INOUT_COL(I)
         ENDDO
      ENDIF

!***********************************************************************************************************************************
  993 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' LAPACK SUBR XERBLA SHOULD HAVE REPORTED AN ERROR ON AN ILLEGAL ARGUMENT IN A CALL TO LAPACK SUBR '    &
                    ,/,15X,A,' (OR A SUBR CALLED BY IT) AND THEN ABORTED')



      RETURN

!***********************************************************************************************************************************

      END SUBROUTINE FBS_LAPACK


      SUBROUTINE SYM_MAT_DECOMP_LAPACK ( CALLING_SUBR, MATIN_NAME, MATIN_SET, NROWS, NTERMS, I_MATIN, J_MATIN, MATIN, PRT_ERRS,    &
                                         MATIN_DIAG_RAT, EQUIL_MATIN, CALC_COND_NUM, DEB_PRT, EQUED, MATIN_SDIA, K_INORM, RCOND,   &
                                         EQUIL_SCALE_FACS, INFO )

! Decomposes a symmetric band matrix into triangular factors. The input matrix, MATIN, is stored in CRS sparse format and is
! converted, in this subr, to band matrix ABAND (stored in module LAPACK_DPB_MATRICES) needed for the LAPACK routines that do the
! actual work

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FACTORED_MATRIX, FATAL_ERR, LINKNO
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, ONEPP6
      USE PARAMS, ONLY                :  BAILOUT, EPSIL, SUPINFO
      USE LAPACK_DPB_MATRICES, ONLY   :  ABAND, LAPACK_S
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG, NDEBUG

      USE LAPACK_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_LAPACK_MAT, DEALLOCATE_LAPACK_MAT
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_MATRIX_BY_ROWS
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  BAILOUT_CHECK, COND_NUM, GET_GRID_AND_COMP, LINK_MESSAGE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SYM_MAT_DECOMP_LAPACK'

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: CALLED_SUBR = ' ' ! Name of a called subr (for output error purposes)
      CHARACTER(LEN=*) , INTENT(IN)   :: CALC_COND_NUM     ! If "Y" calc RCOND (reciprocal of condition number of MATIN)
      CHARACTER(LEN=*) , INTENT(IN)   :: CALLING_SUBR      ! The subr that called this subr (used for output error purposes)
      CHARACTER(LEN=*) , INTENT(IN)   :: EQUIL_MATIN       ! If "Y" attempt to equilibrate MATIN (if it needs it)
      CHARACTER(1*BYTE), INTENT(OUT)  :: EQUED             ! 'Y' if MATIN was equilibrated in subr EQUILIBRATE (called herein)
      CHARACTER(LEN=*) , INTENT(IN)   :: MATIN_DIAG_RAT    ! If "Y" calculate max ratio of matrix diagonal to factor diagonal
      CHARACTER(LEN=*) , INTENT(IN)   :: MATIN_NAME        ! Name of matrix to be decomposed
      CHARACTER(LEN=*) , INTENT(IN)   :: MATIN_SET         ! Set designator for the input matrix. If it corresponds to a MYSTRAN
!                                                            displ set (e.g. 'L ' set) then error messages about singulatities
!                                                            can reference the grid/comp that is singular (otherwise the row/col
!                                                            where the singularity occurs is referenced). If it is not a MYSTRAN
!                                                            set designator it should be blank
      CHARACTER(LEN=*) , INTENT(IN)   :: PRT_ERRS          ! If not 'N', print singularity errors

      CHARACTER( 1*BYTE), PARAMETER   :: INORM    = 'I'    ! Indicates to calculate the infinity norm via LAPACK function DLANSB
      CHARACTER( 1*BYTE)              :: QUIT_ON_POS_INFO  ! Indicator of whether to quit if output value of INFO is found to be > 0
      CHARACTER( 1*BYTE), PARAMETER   :: UPLO     = 'U'    ! Indicates upper triang part of matrix is stored

      INTEGER(LONG), INTENT(IN)       :: DEB_PRT(2)        ! Debug numbers to say whether to write ABAND and/or its decomp to file
      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in sparse matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: NTERMS            ! Number of nonzeros in sparse matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: I_MATIN(NROWS+1)  ! Indicators of number of nonzero terms in rows of matrix MATIN
      INTEGER(LONG), INTENT(IN)       :: J_MATIN(NTERMS)   ! Col numberts of nonzero terms in matrix MATIN

      INTEGER(LONG), INTENT(INOUT)    :: INFO              ! Output from LAPACK routine to do factorization of ABAND
                                                           !   0:  successful exit
                                                           ! < 0:  if INFO = -i, the i-th argument had an illegal value
                                                           ! > 0:  if INFO = i, the leading minor of order i is not pos def
                                                           !       and the factorization (in DPBTRS) could not be completed.

      INTEGER(LONG), INTENT(OUT)      :: MATIN_SDIA        ! No. of superdiags in the MATIN upper triangle
      INTEGER(LONG)                   :: COMPV             ! Component number (1-6) of a grid DOF
      INTEGER(LONG)                   :: GRIDV             ! Grid number
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IIMAX             ! Row/Col in MATIN where max diagonal term occurs
      INTEGER(LONG)                   :: IIMIN             ! Row/Col in MATIN where min diagonal term occurs


      REAL(DOUBLE) , INTENT(IN)       :: MATIN(NTERMS)     ! A small number to compare real zero
      REAL(DOUBLE) , INTENT(OUT)      :: RCOND             ! Recrip of cond no. of MATIN. Determined in  subr COND_NUM
      REAL(DOUBLE) , INTENT(OUT)      :: K_INORM           ! Inf norm of MATIN matrix (det in  subr COND_NUM)
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      REAL(DOUBLE) , INTENT(OUT)      :: EQUIL_SCALE_FACS(NROWS)
                                                           ! LAPACK_S values to return to calling subr

      REAL(DOUBLE)                    :: MATIN_DIAG(NROWS) ! Diagonal terms from MATIN matrix
      REAL(DOUBLE)                    :: FACTOR_DIAG(NROWS)! Diagonal terms from factor
      REAL(DOUBLE)                    :: KRATIO            ! Ratio: MAXKII/MINKII
      REAL(DOUBLE)                    :: MAXKII            ! Maximum diagonal term in MATIN
      REAL(DOUBLE)                    :: MAXIMAX_RATIO     ! Largest of the ratios of matrix diagonal to factor diagonal
      REAL(DOUBLE)                    :: MB_TO_ALLOCATE    ! MB of memory to allocate
      REAL(DOUBLE)                    :: MINKII            ! Minimum diagonal term in MATIN
!xx   REAL(DOUBLE)                    :: SCOND             ! Ratio of min to max scaling factors, LAPACK_S(i), if MATIN is equil'ed.

      LOGICAL                         :: FACTORIZATION_PROBLEM

      REAL(DOUBLE), EXTERNAL          :: DLANSB

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Deallocate ABAND in case it is already allocated

      CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )

! We will abort if errors occur and if INFO, on input, is not -1

      QUIT_ON_POS_INFO = 'Y'
      IF (INFO == -1) THEN
         QUIT_ON_POS_INFO = 'N'
      ENDIF

      EPS1 = EPSIL(1)

! Determine bandwidth of matrix

      CALL LINK_MESSAGE('CALC BANDWIDTH OF MATRIX ' // MATIN_NAME(1:))
      CALL BANDSIZ ( NROWS, NTERMS, I_MATIN, J_MATIN, MATIN_SDIA )
      MB_TO_ALLOCATE = (REAL(DOUBLE))*(REAL(MATIN_SDIA+1))*(REAL(NROWS))/ONEPP6
      WRITE(SC1,3094) MATIN_NAME, MATIN_SDIA+1, MB_TO_ALLOCATE
      WRITE(ERR,3002) MATIN_NAME, MATIN_SDIA+1
      IF (SUPINFO == 'N') THEN
         WRITE(F06,3002) MATIN_NAME, MATIN_SDIA+1
      ENDIF
      IF (MB_TO_ALLOCATE <= ONE) THEN
         WRITE(ERR,3003) MATIN_NAME, MB_TO_ALLOCATE
         IF (SUPINFO == 'N') THEN
            WRITE(F06,3003) MATIN_NAME, MB_TO_ALLOCATE
         ENDIF
      ELSE
         WRITE(ERR,3004) MATIN_NAME, MB_TO_ALLOCATE
         IF (SUPINFO == 'N') THEN
            WRITE(F06,3004) MATIN_NAME, MB_TO_ALLOCATE
         ENDIF
      ENDIF

! Allocate array ABAND (matrix in band form for LAPACK)

      CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR LAPACK BAND FORM OF ' // MATIN_NAME(1:))
      CALL ALLOCATE_LAPACK_MAT ( 'ABAND', MATIN_SDIA+1, NROWS, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'LAPACK_S', NROWS, 1, SUBR_NAME )

! Put MATIN matrix into ABAND form required by LAPACK band matrix.

      CALL LINK_MESSAGE('PUT INTO LAPACK BAND FORM: MATRIX ' // MATIN_NAME(1:))
      CALL BANDGEN_LAPACK_DPB ( MATIN_NAME, NROWS, MATIN_SDIA, NTERMS, I_MATIN, J_MATIN, MATIN, ABAND, SUBR_NAME )

! Output ABAND, if requested

      IF ((DEB_PRT(1) > 0) .AND. (DEB_PRT(1) <= NDEBUG)) THEN
         IF ((DEBUG(DEB_PRT(1)) == 1) .OR. (DEBUG(DEB_PRT(1)) == 3)) THEN
            CALL WRITE_MATRIX_BY_ROWS ( 'LAPACK BAND FORM FOR MATRIX ' // MATIN_NAME(1:), ABAND, MATIN_SDIA+1, NROWS, F06 )
         ENDIF
      ENDIF

! Calc the infinity norm of the matrix using the LAPACK function DLANSB. K_INORM is needed later for estimates
! of errors in the solution. Use array S for workspace in the calculation.

      IF (CALC_COND_NUM == 'Y') THEN
         CALL LINK_MESSAGE('CALC INFINITY NORM OF MATRIX ' // MATIN_NAME(1:))
   !xx   WRITE(SC1, * )
         K_INORM = DLANSB ( INORM, UPLO, NROWS, MATIN_SDIA, ABAND, MATIN_SDIA+1, LAPACK_S )
         WRITE(F06,3005) MATIN_NAME, K_INORM
      ENDIF

! Get max & min diagonals from the original matrix. Code assumes all diag terms positive

      IF (MATIN_DIAG_RAT == 'Y') THEN
         CALL LINK_MESSAGE('GET MAX/MIN DIAGONALS OF MATRIX ' // MATIN_NAME(1:))

         MAXKII = ZERO
         IIMAX  = 0
         DO I=1,NROWS
            IF (ABAND(MATIN_SDIA+1,I) > MAXKII) THEN
               MAXKII = ABAND(MATIN_SDIA+1,I)
               IIMAX  = I
            ENDIF
         ENDDO
         WRITE(F06,3006) MATIN_NAME, MAXKII,IIMAX

         MINKII = MAXKII
         IIMIN  = IIMAX
         DO I=1,NROWS
            IF (ABAND(MATIN_SDIA+1,I) < MINKII) THEN
               MINKII = ABAND(MATIN_SDIA+1,I)
               IIMIN  = I
            ENDIF
         ENDDO
         WRITE(F06,3007) MATIN_NAME, MINKII,IIMIN

         IF (DABS(MINKII) > EPS1) THEN
            KRATIO = MAXKII/MINKII
            WRITE(F06,3008) MATIN_NAME, KRATIO
         ENDIF
      ENDIF

      EQUED = 'N'
! Equilibrate matrix, if user requested it via input arg EQUIL_MATIN
! TEMPORARILY REMOVE THIS CODE. IT WAS CAUSING ERRORS - FAILURES DUE TO RATIO OF MATRIX DIAG TO FACTOR DIAG WHEN EQUILIBRATED
!     IF (EQUIL_MATIN == 'Y') THEN
!        CALL LINK_MESSAGE('EQUILIBRATING (IF NEEDED) MATRIX ' // MATIN_NAME(1:))
!        CALL EQUILIBRATE( MATIN_NAME, MATIN_SET, NROWS, MATIN_SDIA, ABAND, LAPACK_S, EQUED, SCOND )
!     ENDIF

! Set equilibrate scale factors to return to calling subr. We do this since we need to deallocate LAPACK_S here since this
! subr (SYM_MAT_DECOMP_LAPACK) is called other times in MYSTRAN.

      DO I=1,NROWS
         EQUIL_SCALE_FACS(I) = LAPACK_S(I)
      ENDDO

! Deallocate LAPACK_S

      CALL DEALLOCATE_LAPACK_MAT ( 'LAPACK_S' )

! Get max & min diagonals from the equilibrated matrix. Code assumes all diag terms positive

      IF ((EQUED == 'Y') .AND. (MATIN_DIAG_RAT == 'Y')) THEN
         CALL LINK_MESSAGE('GET MAX/MIN DIAGONALS OF EQUILIBRATED MATRIX' // MATIN_NAME(1:))

         MAXKII = ZERO
         IIMAX  = 0
         DO I=1,NROWS
            IF (ABAND(MATIN_SDIA+1,I) > MAXKII) THEN
               MAXKII = ABAND(MATIN_SDIA+1,I)
               IIMAX  = I
            ENDIF
         ENDDO
         WRITE(F06,3009) MATIN_NAME, MAXKII,IIMAX

         MINKII = MAXKII
         IIMIN  = 0
         DO I=1,NROWS
            IF (ABAND(MATIN_SDIA+1,I) < MINKII) THEN
               MINKII = ABAND(MATIN_SDIA+1,I)
               IIMIN  = I
            ENDIF
         ENDDO
         WRITE(F06,3010) MATIN_NAME, MINKII,IIMIN

         IF (DABS(MINKII) > EPS1) THEN
            KRATIO = MAXKII/MINKII
            WRITE(F06,3011) MATIN_NAME, KRATIO
         ENDIF

         IF ((DEB_PRT(1) > 0) .AND. (DEB_PRT(1) <= NDEBUG)) THEN
            IF ((DEBUG(DEB_PRT(1)) == 2) .OR. (DEBUG(DEB_PRT(1)) == 3)) THEN
               CALL WRITE_MATRIX_BY_ROWS ('LAPACK BAND FORM FOR EQUILIBRATED MATRIX' // MATIN_NAME(1:), ABAND, MATIN_SDIA+1,       &
                                           NROWS, F06)
            ENDIF
         ENDIF

      ENDIF

! Perform factorization of matrix. ABAND is the original matrix going into the decomp routine and is the upper triangular factor on
! exit.

      CALL LINK_MESSAGE('LAPACK TRIANGULAR FACTORIZATION OF MATRIX ' // MATIN_NAME(1:))
      CALL DPBTRF ( UPLO, NROWS, MATIN_SDIA, ABAND, MATIN_SDIA+1, INFO )

      CALLED_SUBR = 'DPBTRF'
      IF (INFO == 0) THEN

         FACTORED_MATRIX(1:) = ' '
         FACTORED_MATRIX     = MATIN_NAME

      ELSE IF (INFO < 0) THEN                              ! LAPACK subr XERBLA should have reported error on an illegal argument
!                                                            in calling a LAPACK subr, so we should not have gotten here
         WRITE(ERR,993) SUBR_NAME, CALLED_SUBR, CALLING_SUBR
         WRITE(F06,993) SUBR_NAME, CALLED_SUBR, CALLING_SUBR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit

      ELSE IF (INFO > 0) THEN                              ! Leading minor of MATIN not positive definite

         CALL GET_GRID_AND_COMP ( MATIN_SET, INFO, GRIDV, COMPV  )

         IF (PRT_ERRS /= 'N') THEN
            WRITE(ERR,981) MATIN_NAME, CALLED_SUBR, INFO
            WRITE(F06,981) MATIN_NAME, CALLED_SUBR, INFO
            IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
               WRITE(ERR,9811) GRIDV, COMPV, CALLING_SUBR
               WRITE(F06,9811) GRIDV, COMPV, CALLING_SUBR
            ELSE
               WRITE(ERR,9812) INFO, CALLING_SUBR
               WRITE(F06,9812) INFO, CALLING_SUBR
            ENDIF
         ENDIF

      ENDIF

! Output ABAND (factor of MATIN now), if requested

      IF ((DEB_PRT(2) > 0) .AND. (DEB_PRT(2) <= NDEBUG)) THEN
         IF (DEBUG(DEB_PRT(2)) == 1) THEN
            CALL WRITE_MATRIX_BY_ROWS('TRIANGULAR FACTOR IN LAPACK BAND FORM FOR MATRIX ' // MATIN_NAME(1:), ABAND, MATIN_SDIA+1,  &
                                       NROWS, F06)
         ENDIF
      ENDIF


      DO I=1,NROWS
         FACTOR_DIAG(I) = ABAND(MATIN_SDIA+1,I)
      ENDDO

      FACTORIZATION_PROBLEM = BAILOUT_CHECK( CALLING_SUBR, MATIN_NAME, MATIN_SET, NROWS, NTERMS, I_MATIN, MATIN, PRT_ERRS,         &
                                             FACTOR_DIAG )


      IF (FACTORIZATION_PROBLEM .OR. (INFO > 0)) THEN
                                                           ! If BAILOUT >= 0 then quit. Otherwise, continue processing.
         IF ((BAILOUT >= 0) .AND. (QUIT_ON_POS_INFO == 'Y')) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,99999) BAILOUT
            WRITE(F06,99999) BAILOUT
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ENDIF

! If CALC_COND_NUM = 'Y', then calc the reciprocal of the condition number of the matrix. This is done in
! subr COND_NUM using LAPACK subroutine DPBCON and is used for a better estimate of solution errors later.

      RCOND = ZERO
      IF (CALC_COND_NUM == 'Y') THEN
         CALL LINK_MESSAGE('CALC RECIP OF COND NUM OF MATRIX ' // MATIN_NAME(1:))
         CALL COND_NUM ( MATIN_NAME, NROWS, MATIN_SDIA, K_INORM, ABAND, RCOND )
      ENDIF

!***********************************************************************************************************************************
  981 FORMAT(' *ERROR   981: THE FACTORIZATION OF THE MATRIX ',A,' COULD NOT BE COMPLETED BY LAPACK SUBR ',A                       &
                    ,/,14X,' THE LEADING MINOR OF ORDER ',I12,' IS NOT POSITIVE DEFINITE')

 9811 FORMAT('               THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT ',I8,' COMP ',I3,'. THE CALLING SUBR WAS: ',A,/)

 9812 FORMAT('               THIS IS FOR ROW AND COL ',I8,' IN THE MATRIX. THE CALLING SUBR WAS: ',A,/)

  993 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' LAPACK SUBR XERBLA SHOULD HAVE REPORTED AN ERROR ON AN ILLEGAL ARGUMENT IN A CALL TO LAPACK SUBR: '   &
                    ,/,15X,A,/,' (OR A SUBR CALLED BY IT) AND THEN ABORTED. THE CALLING SUBR WAS: ',A)

 3002 FORMAT(' *INFORMATION: BANDWIDTH OF MATRIX ',A11,'                                        = ',I13,' (BW)',/)

 3003 FORMAT(' *INFORMATION: MEMORY REQUIRED FOR TRIANGULAR DECOMPOSITION OF MATRIX   ',A11,'   = ',F13.6,' MB ',                  &
                           '(8*NDOF*BW)',/,14X,' (using LAPACK band matrix algorithm)',/)

 3004 FORMAT(' *INFORMATION: MEMORY REQUIRED FOR TRIANGULAR DECOMPOSITION OF MATRIX   ',A11,'   = ',F13.3,' MB ',                  &
                           '(8*NDOF*BW)',/,14X,' (using LAPACK band matrix algorithm)',/)

 3005 FORMAT(' *INFORMATION: INFINITY NORM OF MATRIX ',A11,'                                    = ',1ES13.6,/)

 3006 FORMAT(' *INFORMATION: MAXIMUM DIAGONAL TERM IN MATRIX ',A11,'                            = ',1ES13.6,                       &
                           ' Occurs in row/col no. ',I8)

 3007 FORMAT(' *INFORMATION: MINIMUM DIAGONAL TERM IN MATRIX ',A11,'                            = ',1ES13.6,                       &
                           ' Occurs in row/col no. ',I8)

 3008 FORMAT(' *INFORMATION: RATIO OF MAX TO MIN DIAGONALS IN MATRIX ',A11,'                    = ',1ES13.6,/)

 3009 FORMAT(' *INFORMATION: MAXIMUM DIAGONAL TERM IN THE EQUILIBRATED MATRIX ',A11,'           = ',1ES13.6,                       &
                           ' Occurs in row/col no. ',I8)

 3010 FORMAT(' *INFORMATION: MINIMUM DIAGONAL TERM IN THE EQUILIBRATED MATRIX ',A11,'           = ',1ES13.6,                       &
                           ' Occurs in row/col no. ',I8)

 3011 FORMAT(' *INFORMATION: RATIO OF MAX TO MIN DIAGONALS IN THE EQUILIBRATED MATRIX ',A11,'   = ',1ES13.6,/)

 3094 FORMAT(5X,' Bandwidth of ',A,'  = ',I8,' and requires ',F10.3,' MB of memory')


99999 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE MESSAGES AND BULK DATA PARAMETER BAILOUT = ',I7)



      RETURN

!***********************************************************************************************************************************

      END SUBROUTINE SYM_MAT_DECOMP_LAPACK

   END MODULE LAPACK_ADAPTERS
