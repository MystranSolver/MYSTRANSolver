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

       SUBROUTINE EIG_GIV_MGIV

! Solves for eigenvalues and eigenvectors when method is GIV (Givens) or MGIV (modified Givens)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, KLL_SDIA, KLLD_SDIA, MLL_SDIA, NDOFL, NTERM_KLL, NTERM_KLLD,     &
                                         NTERM_MLL, NUM_EIGENS, NUM_KLLD_DIAG_ZEROS, NUM_MLL_DIAG_ZEROS, NVEC, SOL_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  BAILOUT, EPSIL, SUPINFO, SUPWARN
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, PI
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, EIGEN_VEC, MODE_NUM
      USE MODEL_STUF, ONLY            :  EIG_FRQ1, EIG_FRQ2, EIG_METH, EIG_N1, EIG_N2, EIG_VECS
      USE LAPACK_DPB_MATRICES, ONLY   :  ABAND, BBAND
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_KLLD, J_KLLD, KLLD, I_MLL, J_MLL, MLL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE LAPACK_GIV_MGIV_EIG

      USE EIG_GIV_MGIV_USE_IFs
      USE LINK_MESSAGE_Interface

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_GIV_MGIV'
      CHARACTER( 1*BYTE)              :: JOBZ                ! 'V' or 'N' input to subr DSBGVX ( if 'V', calc vecs, 'N' do not).
      CHARACTER( 1*BYTE)              :: IFAIL_NULL          ! 'Y'/'N' indicator if array IFAIL has any nonzero values.
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: CALLED_SUBR = ' '   ! Name of a called subr (for output error purposes)
      CHARACTER( 8*BYTE)              :: NAME2       = ' '   ! Name for output purposes.
      CHARACTER( 1*BYTE)              :: RANGE               ! 'V' or 'I' indicator for LAPACK of whether eigen range is based on
!                                                                eigenvalues or eigenvector numbers.

      INTEGER(LONG)                   :: A_SDIA              ! No. of superdiags in the ABAND matrix upper triangle.
      INTEGER(LONG)                   :: B_SDIA              ! No. of superdiags in the BBAND (MLL or KLLD) matrix upper triangle.
      INTEGER(LONG)                   :: COMPV               ! Component number (1-6) of a grid DOF
      INTEGER(LONG)                   :: GRIDV               ! Grid number
      INTEGER(LONG)                   :: I                   ! DO loop index
      INTEGER(LONG)                   :: IFAIL(NDOFL)        ! For LAPACK - integer numbers of eigenvectors that failed to converge
      INTEGER(LONG)                   :: IFAIL_IND(NDOFL)    ! The integer numbers from array IFAIL that are not 0.
      INTEGER(LONG)                   :: IL, IU              ! For LAPACK - the lower/upper eigenvector numbers in over
!                                                               which LAPACK will calc eigenvectors.
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! NOTE: LAPACK says to dimension         IWORK(5*NDOFL),       This failed on some small DOF problems. Seems like 8*NDOFL works.
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      INTEGER(LONG)                   :: IWORK(8*NDOFL)      ! Integer workspace used by LAPACK.
      INTEGER(LONG)                   :: INFO        = 0     ! Output from LAPACK routine to do eigenvalue, vector calculation:
!                                                                 0:  successful exit
!                                                               < 0:  if INFO = -i, the i-th argument had an illegal value
!                                                               > 0:  if INFO = i, the leading minor of order i is not pos def
!                                                                     and the factorization (in DPBTRS) could not be completed.
      INTEGER(LONG)                   :: KLLD_NULL_ROWS      ! Number of null rows in the KLLD diff stiff matrix.
      INTEGER(LONG)                   :: LDAB                ! Leading dimension of ABAND (banded KLL).
      INTEGER(LONG)                   :: LDBB                ! Leading dimension of BBAND (banded MLL, or KLLD).
      INTEGER(LONG)                   :: LDQ                 ! For LAPACK - leading dimension of matrix Q.
      INTEGER(LONG)                   :: LDZ                 ! For LAPACK - leading dimension of eigenvector matrix.
      INTEGER(LONG)                   :: LMINOR              ! Leading minor number (for error returned from LAPACK)
      INTEGER(LONG)                   :: MLL_NULL_ROWS       ! Number of null rows in the MLL mass matrix.
      INTEGER(LONG)                   :: NUM_FAIL            ! Number of eigenvalues/vectors that failed to converge.
      INTEGER(LONG)                   :: NUM1                ! Number to use for max no. of eigens to find. Must be NUM1 <= NDOFL


      REAL(DOUBLE)                    :: ABSTOL              ! Tolerance number for LAPACK routines.
      REAL(DOUBLE)                    :: EPS1                ! Small number to compare variables against zero.
      REAL(DOUBLE)                    :: VL                  ! EIG_FRQ1 written as an eigenvalue.
      REAL(DOUBLE)                    :: VU                  ! EIG_FRQ2 written as an eigenvalue.
      REAL(DOUBLE)                    :: LAMBDA1             ! Eigenvalue corresponding to EIG_FRQ1.
      REAL(DOUBLE)                    :: LAMBDA2             ! Eigenvalue corresponding to EIG_FRQ2.
      REAL(DOUBLE)                    :: Q(NDOFL,NDOFL)      ! Matrix used in LAPACK reduction of eigen problem to standard form.
      REAL(DOUBLE)                    :: WORK(7*NDOFL)       ! Real workspace for METH = GIV.

      REAL(DOUBLE)                    :: DLAMCH
      EXTERNAL                        :: DLAMCH

      INTRINSIC                       :: DSQRT, MIN



! **********************************************************************************************************************************
      EPS1   = EPSIL(1)

! Determine bandwidth of stiffness and mass matrices so BANDGEN can put them in LAPACK band form

      CALL LINK_MESSAGE('CALCULATE BANDWIDTH OF KLL MATRIX')
      CALL BANDSIZ ( NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL_SDIA )
      WRITE(ERR,4904) KLL_SDIA
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4904) KLL_SDIA
      ENDIF
! --- BANDED_optimizisation -begin-- !
      CALL REPORT_BANDED_STORAGE_ESTIMATE ( 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL_SDIA+1, SUPINFO )
! --- BANDED_optimizisation -end-- !

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL LINK_MESSAGE('CALCULATE BANDWIDTH OF KLLD MATRIX')
         CALL BANDSIZ ( NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, KLLD_SDIA )
         WRITE(ERR,4905) 'KLLD', KLLD_SDIA
         IF (SUPINFO == 'N') THEN
            WRITE(F06,4905) 'KLLD', KLLD_SDIA
         ENDIF
! --- BANDED_optimizisation -begin-- !
         CALL REPORT_BANDED_STORAGE_ESTIMATE ( 'KLLD', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, KLLD_SDIA+1, SUPINFO )
! --- BANDED_optimizisation -end-- !
      ELSE
         CALL LINK_MESSAGE('CALCULATE BANDWIDTH OF MLL MATRIX')
         CALL BANDSIZ ( NDOFL, NTERM_MLL, I_MLL, J_MLL, MLL_SDIA )
         WRITE(ERR,4905) 'MLL', MLL_SDIA
         IF (SUPINFO == 'N') THEN
            WRITE(F06,4905) 'MLL', MLL_SDIA
         ENDIF
! --- BANDED_optimizisation -begin-- !
         CALL REPORT_BANDED_STORAGE_ESTIMATE ( 'MLL', NDOFL, NTERM_MLL, I_MLL, J_MLL, MLL_SDIA+1, SUPINFO )
! --- BANDED_optimizisation -end-- !
      ENDIF

! A_SDIA and B_SDIA are the number of super-diags in the band form of the ABAND stiffness and BBAND mass matrices.
! LAPACK (for GIV or MGIV) requires that A_SDIA >= B_SDIA. ARPACK (for LANCZOS) requires A_SDIA and B_SDIA to be equal.

      IF      (EIG_METH(1:3) == 'GIV ') THEN

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            IF (KLL_SDIA >= KLLD_SDIA) THEN
               A_SDIA = KLL_SDIA
               B_SDIA = KLL_SDIA
            ELSE
               A_SDIA = KLLD_SDIA
               B_SDIA = KLLD_SDIA
            ENDIF
            LDAB = A_SDIA + 1
            LDBB = B_SDIA + 1
         ELSE
            IF (KLL_SDIA >= MLL_SDIA) THEN
               A_SDIA = KLL_SDIA
               B_SDIA = KLL_SDIA
            ELSE
               A_SDIA = MLL_SDIA
               B_SDIA = MLL_SDIA
            ENDIF
            LDAB = A_SDIA + 1
            LDBB = B_SDIA + 1
         ENDIF

      ELSE IF (EIG_METH(1:4) == 'MGIV') THEN               ! ABAND/BBAND roles will be reversed when DSBGVX_GIV_MGIV is called below

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            IF (KLLD_SDIA >= KLL_SDIA) THEN
               A_SDIA = KLLD_SDIA
               B_SDIA = KLLD_SDIA
            ELSE
               A_SDIA = KLL_SDIA
               B_SDIA = KLL_SDIA
            ENDIF
            LDAB = A_SDIA + 1
            LDBB = B_SDIA + 1
         ELSE
            IF (MLL_SDIA >= KLL_SDIA) THEN
               A_SDIA = MLL_SDIA
               B_SDIA = MLL_SDIA
            ELSE
               A_SDIA = KLL_SDIA
               B_SDIA = KLL_SDIA
            ENDIF
            LDAB = A_SDIA + 1
            LDBB = B_SDIA + 1
         ENDIF

      ENDIF

! Allocate arrays ABAND and BBAND (stiffness, mass matrices in band form for LAPACK)

! --- BANDED_optimizisation -begin-- !
      CALL REPORT_SOLVER_DISPATCH_POLICY ( 'KLL/MLL', SUBR_NAME )
! --- BANDED_optimizisation -end-- !
      CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR LAPACK BAND FORM OF KLL')
      CALL ALLOCATE_LAPACK_MAT ( 'ABAND', LDAB, NDOFL, SUBR_NAME )

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR LAPACK BAND FORM OF KLLD')
      ELSE
         CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR LAPACK BAND FORM OF MLL')
      ENDIF
      CALL ALLOCATE_LAPACK_MAT ( 'BBAND', LDBB, NDOFL, SUBR_NAME )

! Put stiffness and mass matrices in form required by LAPACK band matrix and write them out, if requested.

      CALL LINK_MESSAGE('PUT KLL MATRIX IN LAPACK BAND FORM')
      CALL BANDGEN_LAPACK_DPB ( 'KLL', NDOFL, A_SDIA, NTERM_KLL, I_KLL, J_KLL, KLL, ABAND, SUBR_NAME )
      IF ((DEBUG(40) == 1) .OR. (DEBUG(40) == 3)) THEN
         CALL WRITE_MATRIX_BY_ROWS ( 'STIFFNESS MATRIX KLL IN LAPACK BAND FORM', ABAND, LDAB, NDOFL, F06 )
      ENDIF

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL LINK_MESSAGE('PUT KLLD MATRIX IN LAPACK BAND FORM')
         CALL BANDGEN_LAPACK_DPB ( 'KLLD', NDOFL, B_SDIA, NTERM_KLLD, I_KLLD, J_KLLD, KLLD, BBAND, SUBR_NAME )
         IF ((DEBUG(40) == 2) .OR. (DEBUG(40) == 3)) THEN
            CALL WRITE_MATRIX_BY_ROWS ( 'DIFF STIFF MATRIX KLLD IN LAPACK BAND FORM'     , BBAND, LDBB, NDOFL, F06 )
         ENDIF
      ELSE
         CALL LINK_MESSAGE('PUT MLL MATRIX IN LAPACK BAND FORM')
         CALL BANDGEN_LAPACK_DPB ( 'MLL', NDOFL, B_SDIA, NTERM_MLL, I_MLL, J_MLL, MLL, BBAND, SUBR_NAME )
         IF ((DEBUG(40) == 2) .OR. (DEBUG(40) == 3)) THEN
            CALL WRITE_MATRIX_BY_ROWS ( 'MASS MATRIX MLL IN LAPACK BAND FORM'     , BBAND, LDBB, NDOFL, F06 )
         ENDIF
      ENDIF

! If this is not a CB or BUCKLING soln, dellocate arrays for KLL.      ! Keep arrays MLL, KLLD. Need them later to calc gen mass

      IF ((SOL_NAME(1:12) /= 'GEN CB MODEL' ) .AND. (SOL_NAME(1:8) /= 'BUCKLING')) THEN
         CALL LINK_MESSAGE('DEALLOCATE SPARSE KLL ARRAYS')
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLL', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'KLL' )
      ENDIF

      ABSTOL = TWO*DLAMCH( 'S' )

      CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM' , NDOFL, 1, SUBR_NAME )
      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', NDOFL, 1, SUBR_NAME )

      IF (EIG_VECS == 'Y') THEN                            ! Eigenvectors are requested, so set some variables
         JOBZ = 'V'
         LDZ  = NDOFL
         LDQ  = NDOFL
      ELSE
         JOBZ = 'N'
         LDZ  = 1
         LDQ  = 1
      ENDIF

      IF (DEBUG(185) == 0) THEN                            ! If 0, only find finite eigens within the range requested
         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            NUM1 = NDOFL - NUM_KLLD_DIAG_ZEROS
         ELSE
            NUM1 = NDOFL - NUM_MLL_DIAG_ZEROS
         ENDIF
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,48006) NUM1
         IF (SUPWARN == 'N') THEN
            WRITE(F06,48006) NUM1
         ENDIF
      ELSE
         NUM1 = NDOFL
      ENDIF

      IF (EIG_N2 > 0) THEN                                 ! The IL-th through IU-th eigenvalues will be found
         RANGE = 'I'
         IF      (EIG_METH(1:3) == 'GIV' ) THEN
            IL = MIN( EIG_N1, NUM1 )
            IU = MIN( NDOFL, EIG_N2 )
            VL = ZERO
            VU = ZERO
         ELSE IF (EIG_METH(1:4) == 'MGIV') THEN            ! Get IL, IU in reverse order for MGIV (NOTE: here IL > IU)
            IL = NDOFL - (MIN( NUM1 , EIG_N2 ) - 1)
            IU = NDOFL - (MIN( NDOFL, EIG_N1 ) - 1)
            VL = ZERO
            VU = ZERO
         ENDIF
      ELSE                                                 ! All eigenvalues in the half-open interval (VL,VU) will be found
         RANGE = 'V'
         LAMBDA1 = (TWO*PI*EIG_FRQ1)*(TWO*PI*EIG_FRQ1)
         LAMBDA2 = (TWO*PI*EIG_FRQ2)*(TWO*PI*EIG_FRQ2)
         VL = LAMBDA1
         VU = LAMBDA2
         IF (EIG_METH(1:4) == 'MGIV') THEN
            IF (DABS(LAMBDA2) > EPS1) THEN
               VL = ONE/(LAMBDA2)                          ! Flip VL, VU for MGIV
            ELSE
               VL = ONE/DLAMCH('S')
            ENDIF
            IF (DABS(LAMBDA1) > EPS1) THEN
               VU = ONE/(LAMBDA1)                          ! Flip VL, VU for MGIV
            ELSE
               VU = ONE/DLAMCH('S')
            ENDIF
         ENDIF
         IL = 0
         IU = 0
      ENDIF

      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC', NDOFL, LDZ, SUBR_NAME )

      IF (EIG_METH(1:3) == 'GIV') THEN

         CALL LINK_MESSAGE('SOLVE FOR EIGENVALUES/VECTORS - GIV METHOD')
         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            CALL DSBGVX_GIV_MGIV ( JOBZ, RANGE, 'U', NDOFL, A_SDIA, B_SDIA, ABAND, LDAB, -BBAND, LDBB, Q, LDQ, VL, VU, IL, IU,     &
                                   ABSTOL, NUM_EIGENS, EIGEN_VAL, EIGEN_VEC, LDZ, WORK, IWORK, IFAIL, INFO, EIG_METH, MODE_NUM,    &
                                   NVEC )
         ELSE
            CALL DSBGVX_GIV_MGIV ( JOBZ, RANGE, 'U', NDOFL, A_SDIA, B_SDIA, ABAND, LDAB,  BBAND, LDBB, Q, LDQ, VL, VU, IL, IU,     &
                                   ABSTOL, NUM_EIGENS, EIGEN_VAL, EIGEN_VEC, LDZ, WORK, IWORK, IFAIL, INFO, EIG_METH, MODE_NUM,    &
                                   NVEC )
         ENDIF

      ELSE IF (EIG_METH(1:4) == 'MGIV') THEN

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            KLLD_NULL_ROWS = 0
            DO I=2,NDOFL+1
               IF ( I_KLLD(I) == I_KLLD(I-1) ) THEN
                  KLLD_NULL_ROWS = KLLD_NULL_ROWS + 1
               ENDIF
            ENDDO
         ELSE
            MLL_NULL_ROWS = 0
            DO I=2,NDOFL+1
               IF ( I_MLL(I) == I_MLL(I-1) ) THEN
                  MLL_NULL_ROWS = MLL_NULL_ROWS + 1
               ENDIF
            ENDDO
         ENDIF

         CALL LINK_MESSAGE('SOLVE FOR EIGENVALUES/VECTORS - MGIV METHOD')
         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            CALL DSBGVX_GIV_MGIV ( JOBZ, RANGE, 'U', NDOFL, B_SDIA, A_SDIA, -BBAND, LDBB, ABAND, LDAB, Q, LDQ, VL, VU, IL, IU,     &
                                   ABSTOL, NUM_EIGENS, EIGEN_VAL, EIGEN_VEC, LDZ, WORK, IWORK, IFAIL, INFO, EIG_METH, MODE_NUM,    &
                                   NVEC )
         ELSE
            CALL DSBGVX_GIV_MGIV ( JOBZ, RANGE, 'U', NDOFL, B_SDIA, A_SDIA,  BBAND, LDBB, ABAND, LDAB, Q, LDQ, VL, VU, IL, IU,     &
                                   ABSTOL, NUM_EIGENS, EIGEN_VAL, EIGEN_VEC, LDZ, WORK, IWORK, IFAIL, INFO, EIG_METH, MODE_NUM,    &
                                   NVEC )
         ENDIF

         CALL INVERT_EIGENS ( NUM_EIGENS, NDOFL, EIGEN_VAL, EIGEN_VEC, MODE_NUM )

      ENDIF


! If INFO /= 0 there was some error so write message and quit

      IF      (INFO < 0) THEN                              ! LAPACK subr XERBLA should have reported error on an illegal argument
!                                                            in a call to a LAPACK subr, so we should not have gotten here
         CALLED_SUBR = 'DSBGVX_GIV_MGIV'
         WRITE(ERR,993) SUBR_NAME, CALLED_SUBR
         WRITE(F06,993) SUBR_NAME, CALLED_SUBR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit

      ELSE IF (INFO > 0) THEN

         IF (INFO <= NDOFL) THEN                           ! Some eigenvecs failed to converge (eigenval converge fail is reported
!                                                            in subr DSBGVX_GIV_MGIV to avoid confusion here of meaning of INFO)
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,4005) EIG_METH
            IF (SUPWARN == 'N') THEN
               WRITE(f06,4005) EIG_METH
            ENDIF

            NUM_FAIL   = 0
            IFAIL_NULL = 'Y'
            DO I=1,NDOFL
               IF (IFAIL(I) /= 0) THEN
                  IFAIL_NULL = 'N'
                  NUM_FAIL = NUM_FAIL + 1
                  IFAIL_IND(NUM_FAIL) = IFAIL(I)
               ENDIF
            ENDDO
            IF (IFAIL_NULL == 'N') THEN
               WRITE(ERR,40051)
               WRITE(ERR,40052) (IFAIL_IND(I),I=1,NUM_FAIL)
               WRITE(F06,40051)
               WRITE(F06,40052) (IFAIL_IND(I),I=1,NUM_FAIL)
            ENDIF

         ELSE                                              ! INFO > NDOFL indicates lead minor of MLL (or KLLD) or KLL not pos def

            NAME2 = 'No name'
            IF      (EIG_METH(1:3) == 'GIV' ) THEN
               IF (SOL_NAME(1:8) == 'BUCKLING') THEN
                  NAME2 = 'KLLD'
               ELSE
                  NAME2 = 'MLL'
               ENDIF
            ELSE IF (EIG_METH(1:4) == 'MGIV') THEN
               NAME2 = 'KLL'
            ENDIF

            CALLED_SUBR = 'DPBSTF'
            LMINOR = INFO - NDOFL
            CALL GET_GRID_AND_COMP ( 'A ', LMINOR, GRIDV, COMPV  )

            WRITE(ERR,989) NAME2, CALLED_SUBR, LMINOR
            WRITE(F06,989) NAME2, CALLED_SUBR, LMINOR
            FATAL_ERR = FATAL_ERR + 1
            IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
               WRITE(ERR,9892) GRIDV, COMPV
               WRITE(F06,9892) GRIDV, COMPV
            ENDIF

            IF (BAILOUT >= 0) THEN                         ! If BAILOUT >= 0 then quit. Otherwise, continue processing.
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
  989 FORMAT(' *ERROR   989: THE FACTORIZATION OF THE MATRIX ',A,' COULD NOT BE COMPLETED BY LAPACK SUBR ',A                       &
                    ,/,14X,' THE LEADING MINOR OF ORDER ',I12,' IS NOT POSITIVE DEFINITE')

 9892 FORMAT('               THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT ',I8,' COMPONENT ',I3)

  993 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' LAPACK SUBR XERBLA SHOULD HAVE REPORTED AN ERROR ON AN ILLEGAL ARGUMENT IN A CALL TO LAPACK SUBR '    &
                    ,/,15X,A,' (OR A SUBR CALLED BY IT) AND THEN ABORTED')

 4005 FORMAT(' *WARNING    : THE ',A,' METHOD OF EIGENVALUE EXTRACTION FAILED TO CONVERGE ON SOME EIGENVECTORS')

40051 FORMAT('               THEIR INDICES ARE:')

40052 FORMAT(15X,10I8)

 4904 FORMAT(' *INFORMATION: NUMBER OF SUPERDIAGONALS IN THE KLL STIFFNESS MATRIX UPPER TRIANGLE IS = ',I12,/)

 4905 FORMAT(' *INFORMATION: NUMBER OF SUPERDIAGONALS IN THE ',A,' MATRIX UPPER TRIANGLE IS      = ',I12,/)

48006 FORMAT(' *WARNING    : THE L-SET MASS MATRIX HAS ONLY ',I8,' NONZEROS ON ITS DIAGONAL. THERE ARE NO MORE FINITE EIGENVALUES',&
                           ' BEYOND THIS NUMBER'                                                                                   &
                    ,/,14x,' (NOTE: USE OF BULK DATA PARAM ART_MASS WITH SMALL VALUE MAY ALLOW MGIV TO FIND MORE EIGENVALUES)')


12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE EIG_GIV_MGIV
