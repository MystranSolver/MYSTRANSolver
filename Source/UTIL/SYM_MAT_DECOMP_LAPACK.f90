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
      USE PARAMS, ONLY                :  BAILOUT, EPSIL, SPARSTOR, SUPINFO
      USE LAPACK_DPB_MATRICES, ONLY   :  ABAND, LAPACK_S
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG, NDEBUG
      USE LAPACK_LIN_EQN_DPB

      USE SYM_MAT_DECOMP_LAPACK_USE_IFs
      USE LINK_MESSAGE_Interface

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
      INTEGER(LONG)                   :: NUM_NONPOS_DIAG   ! Number of diagonal terms <= EPS1 (not SPD-ready)
      INTEGER(LONG)                   :: NUM_ZERO_DIAG     ! Number of diagonal terms <= ZERO


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
! --- BANDED_optimizisation -begin-- !
      REAL(DOUBLE)                    :: BAND_TERMS_EST    ! Number of terms in compact symmetric band storage
! --- BANDED_optimizisation -end-- !
      REAL(DOUBLE)                    :: MB_TO_ALLOCATE    ! MB of memory to allocate
      REAL(DOUBLE)                    :: MINKII            ! Minimum diagonal term in MATIN
!xx   REAL(DOUBLE)                    :: SCOND             ! Ratio of min to max scaling factors, LAPACK_S(i), if MATIN is equil'ed.

      LOGICAL                         :: FACTORIZATION_PROBLEM

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
! --- BANDED_optimizisation -begin-- !
      BAND_TERMS_EST = REAL(MATIN_SDIA+1,DOUBLE)*REAL(NROWS,DOUBLE)
      MB_TO_ALLOCATE = REAL(DOUBLE,DOUBLE)*BAND_TERMS_EST/ONEPP6
! --- BANDED_optimizisation -end-- !
      WRITE(SC1,3094) MATIN_NAME, MATIN_SDIA+1, MB_TO_ALLOCATE
      WRITE(ERR,3002) MATIN_NAME, MATIN_SDIA+1
      IF (SUPINFO == 'N') THEN
         WRITE(F06,3002) MATIN_NAME, MATIN_SDIA+1
      ENDIF
! --- BANDED_optimizisation -begin-- !
      CALL REPORT_BANDED_STORAGE_ESTIMATE ( MATIN_NAME, NROWS, NTERMS, I_MATIN, J_MATIN, MATIN_SDIA+1, SUPINFO )
! --- BANDED_optimizisation -end-- !
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

! Quick SPD readiness diagnostic for LAPACK band Cholesky (DPBTRF/DPBTRS)

      NUM_NONPOS_DIAG = 0
      NUM_ZERO_DIAG   = 0
      DO I=1,NROWS
         IF (ABAND(MATIN_SDIA+1,I) <= EPS1) NUM_NONPOS_DIAG = NUM_NONPOS_DIAG + 1
         IF (ABAND(MATIN_SDIA+1,I) <= ZERO) NUM_ZERO_DIAG   = NUM_ZERO_DIAG + 1
      ENDDO
      IF (NUM_NONPOS_DIAG > 0) THEN
         WRITE(ERR,3090) MATIN_NAME, NUM_NONPOS_DIAG, NUM_ZERO_DIAG, SPARSTOR
         IF (SUPINFO == 'N') THEN
            WRITE(F06,3090) MATIN_NAME, NUM_NONPOS_DIAG, NUM_ZERO_DIAG, SPARSTOR
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

 3090 FORMAT(' *WARNING    3090: MATRIX ',A11,' HAS ',I10,' DIAGONAL TERM(S) <= EPS1; OF THESE, ',I10,' ARE <= 0.0',            &
                    /,14X,' THIS MATRIX MAY NOT BE SPD, SO LAPACK BAND CHOLESKY (DPBTRF/DPBTRS) CAN FAIL.',                      &
                    /,14X,' CURRENT PARAM SPARSTOR = ',A6,'. IF THIS IS EXPECTED, USE A NON-SPD PATH (E.G. SPARSE SOLVER).',/)

 3094 FORMAT(5X,' Bandwidth of ',A,'  = ',I8,' and requires ',F10.3,' MB of memory')

99999 FORMAT(/,' PROCESSING TERMINATED DUE TO ABOVE MESSAGES AND BULK DATA PARAMETER BAILOUT = ',I7)



      RETURN

!***********************************************************************************************************************************

      END SUBROUTINE SYM_MAT_DECOMP_LAPACK
