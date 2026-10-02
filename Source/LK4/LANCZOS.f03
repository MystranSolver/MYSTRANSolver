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

   MODULE LANCZOS

   USE GIVENS, ONLY: INVERT_EIGENS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: EIG_LANCZOS_ARPACK, EIG_LANCZOS_ARPACK_ADAPTIVE

   CONTAINS

      SUBROUTINE EIG_LANCZOS_ARPACK

! Solves for eigenvalues and eigenvectors when the Lanczos method is requested (on Bulk Data EIGRL entry)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, KMSM_SDIA, LINKNO, NDOFL, NTERM_KLL, NTERM_KLLD, NTERM_KMSM,     &
                                         NTERM_KMSMn, NTERM_KMSMs, NTERM_MLL, NTERM_ULL, NUM_EIGENS, NUM_KLLD_DIAG_ZEROS,          &
                                         NUM_MLL_DIAG_ZEROS, NVEC, SOL_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, PI
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  ARP_TOL, BAILOUT, DARPACK, EIGESTL, EPSIL, MXITERL, SOLLIB, SPARSTOR, SUPINFO,            &
                                         SUPWARN
      USE DOF_TABLES, ONLY            :  TDOFI
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, EIGEN_VEC, MODE_NUM
      USE MODEL_STUF, ONLY            :  EIG_FRQ1, EIG_FRQ2, EIG_LANCZOS_NEV_DELT, EIG_LAP_MAT_TYPE, EIG_MODE, EIG_N1, EIG_N2,     &
                                         EIG_NCVFACL, EIG_SIGMA
      USE ARPACK_MATRICES_1, ONLY     :  IWORK, RESID, RFAC, SELECT, VBAS, WORKD, WORKL
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_KLLD, J_KLLD, KLLD, I_MLL, J_MLL, MLL, SYM_KLL, SYM_KLLD, SYM_MLL,   &
                                         I_KMSM, J_KMSM, KMSM, I_KMSMn, J_KMSMn, KMSMn, I_KMSMs, J_KMSMs, KMSMs

      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, SLU_INFO

      USE ARPACK_LANCZOS_EIG
      USE MYSTRAN_LAPACK_EXT

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE LAPACK_ADAPTERS, ONLY       :  BANDGEN_LAPACK_DGB, BANDGEN_LAPACK_DPB, BANDSIZ
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE LAPACK_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_LAPACK_MAT, DEALLOCATE_LAPACK_MAT
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_MATRIX_BY_ROWS
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_MAT_DIAG_ZEROS
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_SYM_TO_CRS_NONSYM
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS
      USE EIGEN_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_EIGEN1_MAT
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_GRID_AND_COMP, LINK_MESSAGE

      IMPLICIT NONE

      LOGICAL                         :: RVEC              ! = .TRUE. or .FALSE. Specifies whether eigenvectors are to be calculated

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_LANCZOS_ARPACK'
      CHARACTER( 1*BYTE)              :: BMAT              ! BMAT specifies the type of the matrix B that defines the semi-inner
!                                                            product for the ARPACK operator OP (see ARPACK subr DSBAND):
!                                                            = 'I' -> standard eigenvalue problem:    A*x = lambda*x
!                                                            = 'G' -> generalized eigenvalue problem: A*x = lambda*M*x
      CHARACTER( 1*BYTE)              :: HOWMNY            ! How many eigenvecs are wanted and the form of the eigenvec matrix Z
!                                                            = 'A': compute all eigenvecs
!                                                            = 'S': compute some of the eigenvecs, specified by logical array SELECT
!                                                                   (see documentation for ARPACK subr DSBAND)
      CHARACTER( 2*BYTE)              :: WHICH             ! When IPARAM(7)= 1 or 2,  WHICH can be set to any one of the following:
!                                                            = 'LM' -> want the NEV eigenvalues of largest magnitude.
!                                                            = 'SM' -> want the NEV eigenvalues of smallest magnitude.
!                                                            = 'LA' -> want the NEV eigenvalues of largest REAL part.
!                                                            = 'SA' -> want the NEV eigenvalues of smallest REAL part.
!                                                            = 'BE' -> Compute  NEV eigenvalues, half from each end of the spectrum
!                                                                      If NEV is odd, compute 1 more from high end than from low end
!                                                            When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only.

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: CALLED_SUBR = ' ' ! Name of a called subr (for output error purposes)

      INTEGER(LONG)                   :: COMPV             ! Component number (1-6) of a grid DOF
      INTEGER(LONG)                   :: GRIDV             ! Grid number
      INTEGER(LONG)                   :: I,J               ! DO loop index.
      INTEGER(LONG)                   :: INFO_ARPACK = 0   ! Error msg number from ARPACK routine to do eigenvalue/vec calc
      INTEGER(LONG)                   :: INFO_LAPACK = 0   ! Error msg number from ARPACK routine to do eigenvalue/vec calc
      INTEGER(LONG)                   :: IPARAM(11)        ! Integer array input to subr DSBAND (in module ARPACK_LANCZOS_1)
      INTEGER(LONG)                   :: KL                ! Max(number of subdiags of KLL, number of subdiags of MLL or KLLD)
      INTEGER(LONG)                   :: KU                ! Max(number of superdiags of KLL, number of superdiags of MLL or KLLD)
      INTEGER(LONG)                   :: LDRFAC            ! Leading dimension of RFAC (band form of KLL - EIG_SIGMA*MLL, or KLLD)
      INTEGER(LONG)                   :: LWORKL            ! Used to dimension a work array
      INTEGER(LONG)                   :: NEV               ! Number of eigenvalues to find
      INTEGER(LONG)                   :: NCV               ! Number of columns of the matrix V. Represents the dimension of the
!                                                            Lanczos basis constructed by ARPACK subr dsaupd
      INTEGER(LONG)                   :: NUM1              ! Number to use for max no. of eigens to find. Must be NUM1 <= NDOFL
      INTEGER(LONG)                   :: NUM_NEG_TERMS1    ! Number of negative terms on the diagonal of RFAC for EIG_FRQ1
      INTEGER(LONG)                   :: NUM_NEG_TERMS2    ! Number of negative terms on the diagonal of RFAC for EIG_FRQ2
      INTEGER(LONG)                   :: NUM_EST_EIGENS    ! Number of estimated eigens in the freq interval (EIG_FRQ2 - EIG_FRQ1)
      INTEGER(LONG)                   :: NUM_KMSM_DIAG_0   ! Number of zero diagonal terms in KMSM

      INTEGER(LONG)                   :: MIN_NCV, MAX_NCV, LNONZEROS

      REAL(DOUBLE)                    :: EPS1              ! A small number to compare zero to

      INTRINSIC                       :: MIN



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      NUM_EST_EIGENS = 0
      NUM_NEG_TERMS1 = 0
      NUM_NEG_TERMS2 = 0

!xx   IF      (SOLLIB == 'BANDED  ') THEN

         IF (NDOFL < EIGESTL) THEN
            IF (EIG_FRQ2 > EPS1) THEN
               IF (EIG_FRQ1 > EPS1) THEN
                  CALL EST_NUM_EIGENS_BANDED ( EIG_FRQ1, NUM_NEG_TERMS1 )
               ENDIF
               CALL EST_NUM_EIGENS_BANDED ( EIG_FRQ2, NUM_NEG_TERMS2 )
               IF (NUM_NEG_TERMS2 > 0) THEN
                  NUM_EST_EIGENS = NUM_NEG_TERMS2 - NUM_NEG_TERMS1
               ELSE
                  ! ?
               ENDIF
            ENDIF
         ENDIF

!xx   ELSE

!xx      FATAL_ERR = FATAL_ERR + 1
!xx      WRITE(ERR,9991) SUBR_NAME, 'SOLLIB'
!xx      WRITE(F06,9991) SUBR_NAME, 'SOLLIB'
!xx      CALL OUTA_HERE ( 'Y' )

!xx   ENDIF

      IF (EIG_FRQ2 > EPS1) THEN
         WRITE(SC1,101) EIG_FRQ2, NUM_EST_EIGENS
         WRITE(F06,101) EIG_FRQ2, NUM_EST_EIGENS
      ENDIF

! Calc KMSM = KLL - EIG_SIGMA*MLL (or + EIG_SIGMA*KLLD for BUCKLING) where EIG_SIGMA = shift freq

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL MATADD_SSS_NTERM ( NDOFL, 'KLL' , NTERM_KLL , I_KLL , J_KLL , SYM_KLL ,  'eig_sigma*KLLD',                           &
                                                NTERM_KLLD, I_KLLD, J_KLLD, SYM_KLLD, 'KMSM', NTERM_KMSM )
         CALL ALLOCATE_SPARSE_MAT ( 'KMSM', NDOFL, NTERM_KMSM, SUBR_NAME )
         CALL MATADD_SSS       ( NDOFL, 'KLL' , NTERM_KLL , I_KLL , J_KLL , KLL , ONE, 'eig_sigma*KLLD',                           &
                                                NTERM_KLLD, I_KLLD, J_KLLD, KLLD, EIG_SIGMA,                                       &
                                        'KMSM', NTERM_KMSM, I_KMSM, J_KMSM, KMSM )
      ELSE
         CALL MATADD_SSS_NTERM ( NDOFL, 'KLL',  NTERM_KLL , I_KLL , J_KLL , SYM_KLL ,  '-eig_sigma*MLL',                           &
                                                NTERM_MLL , I_MLL , J_MLL , SYM_MLL , 'KMSM', NTERM_KMSM )
         CALL ALLOCATE_SPARSE_MAT ( 'KMSM', NDOFL, NTERM_KMSM, SUBR_NAME )
         CALL MATADD_SSS       ( NDOFL, 'KLL' , NTERM_KLL , I_KLL , J_KLL , KLL , ONE, '-eig_sigma*MLL',                           &
                                                NTERM_MLL , I_MLL , J_MLL , MLL, -EIG_SIGMA,                                       &
                                        'KMSM', NTERM_KMSM, I_KMSM, J_KMSM, KMSM )
      ENDIF


! Det bandwidth of KMSM so BANDGEN can put it in LAPACK band form. KMSM_SDIA is the number of super-diags in the band form of KMSM

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL LINK_MESSAGE('CALCULATE BANDWIDTH OF [KLL + sigma*KLLD]')
      ELSE
         CALL LINK_MESSAGE('CALCULATE BANDWIDTH OF [KLL - sigma*MLL]')
      ENDIF
      CALL BANDSIZ ( NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM_SDIA )
      WRITE(ERR,4905) KMSM_SDIA
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4905) KMSM_SDIA
      ENDIF

! EIG_LAP_MAT_TYPE was checked in BD_EIGRL for correctness, but make sure, here, that it is correct

      IF(SOLLIB(1:6) == 'SPARSE') THEN
         LDRFAC = 1 ! make it small because it's not used but parameter has to be allocated.
      ELSEIF      (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
         LDRFAC = KMSM_SDIA + 1
      ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DGB') THEN
         LDRFAC = 3*KMSM_SDIA + 1
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,4003) SUBR_NAME, EIG_LAP_MAT_TYPE
         WRITE(F06,4003) SUBR_NAME, EIG_LAP_MAT_TYPE
         CALL OUTA_HERE ( 'Y' )
      ENDIF


! Allocate array RFAC = (KLL - EIG_SIGMA*MLL, or KLL + EIG_SIGMA*KLLD) for ARACK

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CALL LINK_MESSAGE('ALLOCATE ARPACK BAND MAT: RFAC = KLL + sigma*KLLD')
      ELSE
         CALL LINK_MESSAGE('ALLOCATE ARPACK BAND MAT: RFAC = KLL - sigma*MLL')
      ENDIF
      CALL ALLOCATE_LAPACK_MAT ( 'RFAC', LDRFAC, NDOFL, SUBR_NAME )


      IF(SOLLIB(1:6) == 'BANDED') THEN

! Put KMSM in form required by LAPACK band matrix. Call result array RFAC

         CALL LINK_MESSAGE('PUT RFAC MATRIX IN ARPACK BAND FORM')
         IF      (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
            CALL BANDGEN_LAPACK_DPB ( 'KMSM', NDOFL, KMSM_SDIA, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, RFAC, SUBR_NAME )
         ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DGB') THEN
            CALL BANDGEN_LAPACK_DGB ( 'KMSM', NDOFL, KMSM_SDIA, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, RFAC, SUBR_NAME )
         ENDIF

! Write RFAC, if requested

         IF (DEBUG(40) == 1) THEN
            IF      (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
               CALL WRITE_MATRIX_BY_ROWS ( 'MATRIX [KLL - sigma*KLLD] IN LAPACK BAND FORM', RFAC, LDRFAC, NDOFL, F06 )
            ELSE
               CALL WRITE_MATRIX_BY_ROWS ( 'MATRIX [KLL - sigma*MLL ] IN LAPACK BAND FORM', RFAC, LDRFAC, NDOFL, F06 )
            ENDIF
         ENDIF

      ENDIF

! If this is not a CB or BUCKLING soln, dellocate arrays for KLL.      ! Keep arrays MLL, KLLD. Need them later to calc gen mass

      IF ((SOL_NAME(1:12) /= 'GEN CB MODEL' ) .AND. (SOL_NAME(1:8) /= 'BUCKLING')) THEN
         CALL LINK_MESSAGE('DEALLOCATE SPARSE KLL ARRAYS')
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLL', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'KLL' )
      ENDIF

! Generate nonsymmetric storage form for KMSM. This is done since subr MATMULT_SFF, used when subr DSBAND is called herein,
! will run faster. MATMULT_SFF is called in each "Reverse commumication loop" in DSBAND.

      IF      (SPARSTOR == 'SYM   ') THEN

         CALL SPARSE_MAT_DIAG_ZEROS ( 'KMSM', NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, NUM_KMSM_DIAG_0 )
         NTERM_KMSMn = 2*NTERM_KMSM  - (NDOFL - NUM_KMSM_DIAG_0)
         CALL LINK_MESSAGE('ALLOCATE SPARSE KMSMn ARRAYS')
         CALL ALLOCATE_SPARSE_MAT ( 'KMSMn', NDOFL, NTERM_KMSMn, SUBR_NAME )

         CALL LINK_MESSAGE('CONVERT SYM CRS KMSM TO NONSYM CRS KMSMn')
         CALL CRS_SYM_TO_CRS_NONSYM ( 'KMSM', NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, 'KMSMn', NTERM_KMSMn, I_KMSMn, J_KMSMn,     &
                                       KMSMn, 'Y' )

      ELSE IF (SPARSTOR == 'NONSYM') THEN

         NTERM_KMSMn = NTERM_KMSM
         CALL ALLOCATE_SPARSE_MAT ( 'KMSMn', NDOFL, NTERM_KMSMn, SUBR_NAME )
         DO I=1,NDOFL+1
            I_KMSMn(I) = I_KMSM(I)
         ENDDO
         DO J=1,NTERM_KMSMn
            J_KMSMn(J) = J_KMSM(J)
              KMSMn(J) =   KMSM(J)
         ENDDO

      ELSE                                                 !      Error - incorrect SPARSTOR

         WRITE(ERR,932) SUBR_NAME, SPARSTOR
         WRITE(F06,932) SUBR_NAME, SPARSTOR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF

      IF (DEBUG(49) > 0) THEN
         CALL WRITE_SPARSE_CRS ( ' KMSMn', 'A ', 'A ', NTERM_KMSMn, NDOFL, I_KMSMn, J_KMSMn, KMSMn )
      ENDIF

! Now we can deallocate KMSM (since KMSMn will be used in subr DSBAND)

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KMSM', CR13
      CALL DEALLOCATE_SPARSE_MAT ( 'KMSM' )

! Check that user did not ask for more than NDOFL-1 eigens (Lanczos can't be used to find all). If request is > NDOFL-1, decrease
! request to NDOFL-1 and give warning
      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         LNONZEROS = NDOFL - NUM_KLLD_DIAG_ZEROS
      ELSE
         LNONZEROS = NDOFL - NUM_MLL_DIAG_ZEROS
      ENDIF


      IF (NUM_EST_EIGENS > 0) THEN
         NEV = NUM_EST_EIGENS + EIG_LANCZOS_NEV_DELT
      ELSE
         ! prevent the addition of DARPACK from crashing ARPACK
         IF ((EIG_N2 + DARPACK) > (LNONZEROS - 4)) THEN
            WRITE(ERR,9775) NDOFL, DARPACK
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9775) NDOFL, DARPACK
            ENDIF
            DARPACK = 0
         END IF
         NEV = EIG_N2 + DARPACK
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
      IF (NEV > (NUM1-1)) THEN
!xx      WARN_ERR = WARN_ERR + 1
!xx      WRITE(ERR,4901) EIG_N2, NDOFL-1-DARPACK, NDOFL-1-DARPACK
!xx      IF (SUPWARN == 'N') THEN
!xx         WRITE(F06,4901) EIG_N2, NDOFL-1-DARPACK, NDOFL-1-DARPACK
!xx      ENDIF
         NEV = NUM1 - 1
      ENDIF

      MIN_NCV = NEV + 2
      MAX_NCV = LNONZEROS - 2

      ! sanity check on feasible NCV range
      IF (MIN_NCV > MAX_NCV) THEN
         WRITE(ERR,9776) NDOFL, NEV
         WRITE(F06,9776) NDOFL, NEV
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      END IF

      ! NCV determination step (try user factor first)
      DO I=EIG_NCVFACL,2,-1
         NCV = I*NEV
         IF (NCV >= MIN_NCV .AND. NCV <= MAX_NCV) THEN
            EXIT
         END IF
      END DO

      ! still invalid? (too large OR too small)
      IF (NCV > MAX_NCV .OR. NCV < MIN_NCV) THEN
         DO I=5,2,-1
            NCV = NEV+I
            IF (NCV >= MIN_NCV .AND. NCV <= MAX_NCV) THEN
               EXIT
            END IF
         END DO
      END IF

      ! no valid NCV.
      IF (NCV > MAX_NCV .OR. NCV < MIN_NCV) THEN
         WRITE(ERR,9777) NDOFL, NEV
         WRITE(F06,9777) NDOFL, NEV
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      END IF

      ! inform the user if we didn't use EIG_NCFACL*NEV
      IF (NCV < (EIG_NCVFACL*NEV)) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,9778) NEV, NDOFL, NCV, (EIG_NCVFACL*NEV)
         IF (SUPWARN == 'N') THEN
            WRITE(F06,9778) NEV, NDOFL, NCV, (EIG_NCVFACL*NEV)
         ENDIF
      END IF


      LWORKL = NCV*(NCV + 8)                               ! Size of WORKL workspace used by ARPACK

      KL = KMSM_SDIA
      KU = KL

! Set IPARAM, used in the ARPACK code

      DO I=1,11
         IPARAM(I) = 0                                     ! IPARAM(1) = 1 is set in subr DSBAND
      ENDDO
      IPARAM(3) = MXITERL                                  ! Max number of Arnoldi update iterations allowed
      IPARAM(4) = 1                                        ! NB, the blocksize used in the recurrence. ARPACK only works for NB = 1
      IPARAM(7) = EIG_MODE                                 ! If EIG_MODE = 2 solve (1/lambda)*K*x = M*x, M symmetric, K sym pos def
!                                                            If EIG_MODE = 3 solve (K - EIG_SIGMA*M)x = lbar*Mx (shift invert mode)
!                                                            where, for EIG_MODE = 3 the eigenvalues are lambda = EIG_SIGMA + lbar
!                                                            NOTE: EIG_MODE was checked in subr BD_EIGRL for correct value
! Set other inputs to DSBAND (see above for definitions)

      RVEC   = .TRUE.
      HOWMNY = 'A'
      IF (EIG_MODE == 1) THEN
         BMAT   = 'I'
      ELSE
         BMAT   = 'G'
      ENDIF

      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC', NDOFL, NEV, SUBR_NAME )
      CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM' , NDOFL, 1, SUBR_NAME )
      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', NDOFL, 1, SUBR_NAME )

      CALL ALLOCATE_LAPACK_MAT ( 'IWORK' , NDOFL  , 1, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'RESID' , NDOFL  , 1, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'SELECT', NCV    , 1, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'VBAS'  , NDOFL  , NCV, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'WORKD' , 3*NDOFL, 1, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'WORKL' , LWORKL , 1, SUBR_NAME )

      DO I=1,NCV                                           ! With HOWMNY = 'A' we are calc'ing eigenvecs for all eigenvalues found
         SELECT(I) = .FALSE.                               ! so all members of SELECT are .FALSE.
      ENDDO

      CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - LANCZOS METH')
!xx   WRITE(SC1, * )                                       ! Make new line for DTBSV pass messages that overwrite each other

      IF (DEBUG(50) == 1) CALL DEBUG_EIG_LANCZOS
      WHICH = 'LM'
      CALL DSBAND ( RVEC, HOWMNY, SELECT, EIGEN_VAL, EIGEN_VEC, NDOFL, EIG_SIGMA, NDOFL, LDRFAC, RFAC, KL, KU, WHICH, BMAT,        &
                    NEV, ARP_TOL, RESID, NCV, VBAS, NDOFL, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO_ARPACK,                       &
                    INFO_LAPACK, 'Y', DEBUG(47) )

      NVEC       = IPARAM(5)                               ! With HOWMNY = 'A' we are calc'ing eigenvecs for all eigenvalues found
      NUM_EIGENS = IPARAM(5)

      ! Free SuperLU factorization stored in SLU_FACTORS by SYM_MAT_DECOMP_SUPRLU inside DSBAND (SPARSE solver).
      ! Required to avoid leaking the factor across multi-METHOD MODES iterations.
      IF (SOLLIB(1:6) == 'SPARSE') THEN
         BLOCK
            REAL(DOUBLE) :: DUM_COL(1)
            DUM_COL(1) = ZERO
            CALL C_FORTRAN_DGSSV ( 3, NDOFL, NTERM_KMSMn, 1, KMSMn, J_KMSMn, I_KMSMn, DUM_COL, NDOFL, SLU_FACTORS, SLU_INFO )
         END BLOCK
      ENDIF

!xx   WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KMSMn ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KMSMn' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate IWORK ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'IWORK' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RFAC  ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'RFAC' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RESID ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'RESID' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate SELECT', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'SELECT' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate VBAS  ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'VBAS' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate WORKD ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'WORKD' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate WORKL ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'WORKL' )

      WRITE(ERR,4002) IPARAM(3)                            ! DSBAND returns the actual number of Arnoldi update iterations taken
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4002) IPARAM(3)
      ENDIF

      IF ((EIG_MODE == 1) .OR. (EIG_MODE == 2)) THEN       ! Modes came from subr DSBAND in reverse order since 1, 2 means we solved
         DO I=1,NUM_EIGENS                                 ! for 1/lambda and when we invert the eigens will go from high to low
            MODE_NUM(I) = NUM_EIGENS + 1 - I
         ENDDO
      ELSE IF (EIG_MODE == 3) THEN                         ! Modes came from subr DSBAND in correct order since we solved for lambda
         DO I=1,NUM_EIGENS                                 ! not 1/lambda
            MODE_NUM(I) = I
         ENDDO
      ENDIF

      IF (INFO_ARPACK < 0) THEN
         WRITE(ERR,9996)
         WRITE(F06,9996)
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALLED_SUBR = 'DPBTRF'
      IF      (INFO_LAPACK < 0) THEN                       ! LAPACK subr XERBLA should have reported error on an illegal argument
!                                                            in calling a LAPACK subr, so we should not have gotten here
         WRITE(ERR,993) SUBR_NAME, CALLED_SUBR
         WRITE(F06,993) SUBR_NAME, CALLED_SUBR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ELSE IF (INFO_LAPACK > 0) THEN                       ! NOTE: this will only occur from call to LAPACK dgbtrs in ARPACK

         CALL GET_GRID_AND_COMP ( 'A ', INFO_LAPACK, GRIDV, COMPV  )

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            WRITE(ERR,989) '[KLL + SIGMA*KLLD]', CALLED_SUBR, INFO_LAPACK
            WRITE(F06,989) '[KLL + SIGMA*KLLD]', CALLED_SUBR, INFO_LAPACK
         ELSE
            WRITE(ERR,989) '[KLL - SIGMA*MLL]' , CALLED_SUBR, INFO_LAPACK
            WRITE(F06,989) '[KLL - SIGMA*MLL]' , CALLED_SUBR, INFO_LAPACK
         ENDIF
         FATAL_ERR = FATAL_ERR + 1
         IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
            WRITE(ERR,9892) GRIDV, COMPV
            WRITE(F06,9892) GRIDV, COMPV
         ENDIF
         IF (BAILOUT >= 0) THEN                            ! If BAILOUT >= 0 then quit (leading minor <= 0)
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ENDIF

      IF ((EIG_MODE == 1) .OR. (EIG_MODE == 2)) THEN       ! We found inverses of eigens, so invert them
         CALL INVERT_EIGENS ( NUM_EIGENS, NDOFL, EIGEN_VAL, EIGEN_VEC, MODE_NUM )
      ENDIF

! With HOWMNY = 'A' we are calculating eigenvecs for all eigenvalues found

      NVEC = NUM_EIGENS - DARPACK                          ! Get rid of the higher DARPACK modes in LANCZOS
      NUM_EIGENS = NVEC




      RETURN

! **********************************************************************************************************************************
  101 FORMAT(6X,'Estimated number of eigens below ',F8.0,' Hz is ',I7)

  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

  989 FORMAT(' *ERROR   989: THE FACTORIZATION OF THE MATRIX ',A,' COULD NOT BE COMPLETED BY LAPACK SUBR ',A                       &
                    ,/,14X,' THE LEADING MINOR OF ORDER ',I12,' IS NOT POSITIVE DEFINITE')

  993 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' LAPACK SUBR XERBLA SHOULD HAVE REPORTED AN ERROR ON AN ILLEGAL ARGUMENT IN A CALL TO LAPACK SUBR '    &
                    ,/,15X,A,' (OR A SUBR CALLED BY IT) AND THEN ABORTED',/)

 4002 FORMAT(' *INFORMATION: THE ACTUAL NUMBER OF ARNOLDI UPDATE ITERATIONS TAKEN IN THE SOLUTION FOR THE EIGENVALUES WAS = ',I8)

 4003 FORMAT(' *ERROR  4003: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' EIG_LAP_MAT_TYPE MUST BE EITHER "DGB" OR "DPB" BUT IS = ',A)

 4901 FORMAT(' *WARNING    : REQUEST FOR ',I8,' EIGENVALUES CANNOT BE HONORED. LANCZOS CAN BE USED TO FIND NO MORE THAN',          &
                           ' NEV = '                                                                                               &
                 ,I8,/,14X,' ATTEMPT WILL BE MADE TO FIND ',I8,' EIGENVALUES (SEE MANUAL FOR DEFINITION OF PARAMETER DARPACK)')

 4905 FORMAT(' *INFORMATION: NUMBER OF SUPERDIAGONALS IN THE KMSM = [KLL - sigma*MLL] MATRIX UPPER TRIANGLE IS      = ',I12,/)

 9892 FORMAT('               THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT ',I8,' COMPONENT ',I3)

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

 9996 FORMAT('  PROCESSING STOPPED DUE TO ABOVE ARPACK ERRORS')

12345 FORMAT(A,10X,A)

48006 FORMAT(' *WARNING    : THE L-SET MASS MATRIX HAS ONLY ',I8,' NONZEROS ON ITS DIAGONAL. THERE ARE NO MORE FINITE EIGENVALUES',&
                           ' BEYOND THIS NUMBER'                                                                                   &
                    ,/,14x,' (NOTE: USE OF BULK DATA PARAM ART_MASS WITH SMALL VALUE MAY ALLOW MGIV TO FIND MORE EIGENVALUES)')

 9775 FORMAT(' *WARNING    : DUE TO THE SIZE OF THE PROBLEM (NDOFL=',I0,'), PARAMETER DARPACK',/,15X,&
             'WAS SET TO ZERO TO PREVENT ARPACK FROM CRASHING. ORIGINAL VALUE: ',I0)

 9776 FORMAT(' *ERROR  9776: TOO MANY EIGENVALUES REQUESTED FOR THIS PROBLEM SIZE: NDOFL=',I0,', NEV=',I0,'.')

 9777 FORMAT(' *ERROR  9777: UNABLE TO DETERMINE KRYLOV SUBSPACE SIZE NCV. NDOFL=',I0,', NEV=',I0,'.',/,15X,&
             'USER HAS LIKELY REQUESTED TOO MANY EIGENVALUES FOR THIS MODEL.',/,15X,&
             'NCV MUST BE AT LEAST NEV+2, BUT NO MORE THAN NDOFL-2.')

 9778 FORMAT(' *WARNING    : DUE TO THE REQUESTED EIGENVALUES (NEV=',I0,') RELATIVE TO MODEL SIZE (NDOFL=',I0,'),',/,15X,&
             'NCVFACL WAS NOT USED FOR DETERMINATION OF THE KRYLOV SUBSPACE SIZE NCV.',/,15X,&
             'NCV WAS SET TO ',I0,' INSTEAD OF NEV*NCVFACL=',I0,'.')


! **********************************************************************************************************************************

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE EST_NUM_EIGENS_BANDED ( FREQ, NUM_NEG_TERMS )


      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS
      IMPLICIT NONE

      INTEGER(LONG)                   :: INFO              !
      INTEGER(LONG)                   :: NUM_NEG_TERMS     ! Number of negative terms on the diagonal of RFAC

      REAL(DOUBLE)                    :: DIAG(NDOFL)       !
      REAL(DOUBLE)                    :: FREQ              !
      REAL(DOUBLE)                    :: OFF_DIAG(NDOFL-1) !
      REAL(DOUBLE)                    :: QMAT(1,NDOFL)     !
      REAL(DOUBLE)                    :: SIGMA             !
      REAL(DOUBLE)                    :: WORK(NDOFL)       !

! **********************************************************************************************************************************
      NUM_NEG_TERMS = 0
      SIGMA = (TWO*PI*FREQ)**2

! Add KLL minus SIGMA*MLL used in Lanczos

      CALL MATADD_SSS_NTERM ( NDOFL, 'KLL', NTERM_KLL,  I_KLL,  J_KLL, SYM_KLL, '-SIGMA*MLL', NTERM_MLL, I_MLL, J_MLL, SYM_MLL,&
                                    'KMSM', NTERM_KMSM )
      CALL ALLOCATE_SPARSE_MAT ( 'KMSM', NDOFL, NTERM_KMSM, SUBR_NAME )
      CALL MATADD_SSS      ( NDOFL, 'KLL',  NTERM_KLL,  I_KLL,  J_KLL,  KLL,  ONE, '-SIGMA*MLL', NTERM_MLL, I_MLL, J_MLL, MLL, &
                            -SIGMA, 'KMSM', NTERM_KMSM, I_KMSM, J_KMSM, KMSM )

! Det bandwidth of KMSM so BANDGEN can put it in LAPACK band form. KMSM_SDIA is the number of super-diags in the band form of KMSM

      CALL BANDSIZ ( NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM_SDIA )

      IF      (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
         LDRFAC = KMSM_SDIA + 1
      ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DGB') THEN
         LDRFAC = 3*KMSM_SDIA + 1
      ENDIF

! Allocate array RFAC = (KLL - SIGMA*MLL) for ARACK

      CALL ALLOCATE_LAPACK_MAT ( 'RFAC', LDRFAC, NDOFL, SUBR_NAME )

! Put KMSM in form required by LAPACK band matrix. Call result array RFAC

      IF      (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
         CALL BANDGEN_LAPACK_DPB ( 'KMSM', NDOFL, KMSM_SDIA, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, RFAC, SUBR_NAME )
      ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DGB') THEN
         CALL BANDGEN_LAPACK_DGB ( 'KMSM', NDOFL, KMSM_SDIA, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, RFAC, SUBR_NAME )
      ENDIF

! Reduce KMSM to tridiag form and do Cholesky L*D*L' decomp on it

      CALL DSBTRD ( 'N', 'U', NDOFL, KMSM_SDIA, RFAC, KMSM_SDIA+1, DIAG, OFF_DIAG, QMAT, 1, WORK, INFO )
      CALL DPTTRF_MYSTRAN ( NDOFL, DIAG, OFF_DIAG, INFO )
      NUM_NEG_TERMS = 0
      DO I=1,NDOFL
         IF (DIAG(I) < ZERO) THEN
            NUM_NEG_TERMS = NUM_NEG_TERMS + 1
         ENDIF
      ENDDO

      CALL DEALLOCATE_SPARSE_MAT ( 'KMSM' )
      CALL DEALLOCATE_LAPACK_MAT ( 'RFAC' )

      RETURN

! **********************************************************************************************************************************
12345 FORMAT(A,10X)

! **********************************************************************************************************************************

      END SUBROUTINE EST_NUM_EIGENS_BANDED

! ##################################################################################################################################

      SUBROUTINE DEBUG_EIG_LANCZOS

      WRITE(F06,'(A,I8)')    ' IN EIG_LANCZOS_ARPACK: NUM_EST_EIGENS       = ', NUM_EST_EIGENS
      WRITE(F06,'(A,I8)')    '                      : EIG_LANCZOS_NEV_DELT = ', EIG_LANCZOS_NEV_DELT
      WRITE(F06,'(A,I8)')    '                      : EIG_N2               = ', EIG_N2
      WRITE(F06,'(A,I8)')    '                      : DARPACK              = ', DARPACK
      WRITE(F06,'(A,I8)')    '                      : NDOFL                = ', NDOFL
      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         WRITE(F06,'(A,I8)') '                      : NUM_KLLD_DIAG_ZEROS  = ', NUM_KLLD_DIAG_ZEROS
      ELSE
         WRITE(F06,'(A,I8)') '                      : NUM_MLL_DIAG_ZEROS   = ', NUM_MLL_DIAG_ZEROS
      ENDIF
      WRITE(F06,'(A,I8)')    '                      : NUM1                 = ', NUM1
      WRITE(F06,'(A,I8)')    '                      : EIG_NCVFACL          = ', EIG_NCVFACL
      WRITE(F06,'(A,I8)')    '                      : NCV                  = ', NCV
      WRITE(F06,'(A,I8)')    '                      : NEV                  = ', NEV
      WRITE(F06,*)

! **********************************************************************************************************************************

      END SUBROUTINE DEBUG_EIG_LANCZOS

      END SUBROUTINE EIG_LANCZOS_ARPACK

! End MIT license text.

      SUBROUTINE EIG_LANCZOS_ARPACK_ADAPTIVE

! Solves for eigenvalues and eigenvectors using Lanczos method with ADAPTIVE NEV selection.
! Instead of using inertia-based estimation, this routine starts with a small NEV and doubles it
! until eigenvalues are found outside the specified frequency range [EIG_FRQ1, EIG_FRQ2] on both sides.
!
! This routine:
!   - Sets sigma at the CENTER of the frequency range for optimal shift-invert convergence
!   - Iteratively increases NEV (starting at 10, doubling each iteration, max 10 doublings)
!   - Stops when eigenvalues are found below EIG_FRQ1 AND above EIG_FRQ2, or when EIG_N2 is reached
!   - Filters final results to return only eigenvalues within [EIG_FRQ1, EIG_FRQ2]
!
! NOTE: This routine only supports normal modes (not BUCKLING or GEN CB MODEL)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, KMSM_SDIA, LINKNO, NDOFL, NTERM_KLL, NTERM_KMSM,                 &
                                         NTERM_KMSMn, NTERM_MLL, NUM_EIGENS, NUM_MLL_DIAG_ZEROS, NVEC, SOL_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, PI
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  ARP_TOL, BAILOUT, EPSIL, MXITERL, SOLLIB, SPARSTOR, SUPINFO, SUPWARN
      USE DOF_TABLES, ONLY            :  TDOFI
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, EIGEN_VEC, MODE_NUM
      USE MODEL_STUF, ONLY            :  EIG_FRQ1, EIG_FRQ2, EIG_LAP_MAT_TYPE, EIG_N2, EIG_NCVFACL
      USE ARPACK_MATRICES_1, ONLY     :  IWORK, RESID, RFAC, SELECT, VBAS, WORKD, WORKL
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_MLL, J_MLL, MLL, SYM_KLL, SYM_MLL,                                   &
                                         I_KMSM, J_KMSM, KMSM, I_KMSMn, J_KMSMn, KMSMn
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, SLU_INFO

      USE ARPACK_LANCZOS_EIG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE LAPACK_ADAPTERS, ONLY       :  BANDGEN_LAPACK_DGB, BANDGEN_LAPACK_DPB, BANDSIZ
      USE LAPACK_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_LAPACK_MAT, DEALLOCATE_LAPACK_MAT
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_MAT_DIAG_ZEROS
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_SYM_TO_CRS_NONSYM
      USE EIGEN_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_EIGEN1_MAT, DEALLOCATE_EIGEN1_MAT
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_GRID_AND_COMP, LINK_MESSAGE
      USE DSBAND_PREFAC_MOD, ONLY :  DSBAND_PREFAC
      USE SUPERLU_ADAPTERS, ONLY      :  SYM_MAT_DECOMP_SUPRLU

      IMPLICIT NONE

      LOGICAL                         :: RVEC              ! Specifies whether eigenvectors are to be calculated
      LOGICAL                         :: FOUND_BELOW_FRQ1  ! True if we found an eigenvalue below EIG_FRQ1
      LOGICAL                         :: FOUND_ABOVE_FRQ2  ! True if we found an eigenvalue above EIG_FRQ2
      LOGICAL                         :: SEARCH_COMPLETE   ! True when adaptive search should stop
      LOGICAL                         :: ZERO_LOWER_BOUND  ! True when EIG_FRQ1 ~ 0 (one-sided search)

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! Carriage return for screen output
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_LANCZOS_ARPACK_ADAPTIVE'
      CHARACTER( 1*BYTE)              :: BMAT              ! 'G' for generalized eigenvalue problem
      CHARACTER( 1*BYTE)              :: HOWMNY            ! 'A' to compute all eigenvectors
      CHARACTER( 2*BYTE)              :: WHICH             ! 'LM' for largest magnitude (closest to sigma in shift-invert)
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: CALLED_SUBR = ' ' ! Name of called subr for error messages

      INTEGER(LONG)                   :: COMPV             ! Component number (1-6) of a grid DOF
      INTEGER(LONG)                   :: GRIDV             ! Grid number
      INTEGER(LONG)                   :: I, J, K           ! Loop indices
      INTEGER(LONG)                   :: INFO_ARPACK = 0   ! Error from ARPACK
      INTEGER(LONG)                   :: INFO_LAPACK = 0   ! Error from LAPACK
      INTEGER(LONG)                   :: IPARAM(11)        ! ARPACK parameters
      INTEGER(LONG)                   :: KL                ! Number of subdiagonals
      INTEGER(LONG)                   :: KU                ! Number of superdiagonals
      INTEGER(LONG)                   :: LDRFAC            ! Leading dimension of RFAC
      INTEGER(LONG)                   :: LWORKL            ! Size of WORKL workspace
      INTEGER(LONG)                   :: NEV               ! Number of eigenvalues to find (current trial)
      INTEGER(LONG)                   :: NCV               ! Lanczos basis dimension
      INTEGER(LONG)                   :: NUM_KMSM_DIAG_0   ! Number of zero diagonal terms in KMSM
      INTEGER(LONG)                   :: NUM_CONVERGED     ! Number of converged eigenvalues from DSBAND
      INTEGER(LONG)                   :: NUM_IN_RANGE      ! Number of eigenvalues within [FRQ1, FRQ2]
      INTEGER(LONG)                   :: PREV_IN_RANGE     ! In-range count from previous iteration
      INTEGER(LONG)                   :: NUM_ABOVE_FRQ2    ! Number of eigenvalues above FRQ2
      INTEGER(LONG)                   :: MIN_ABOVE_FOR_STOP ! Minimum modes above FRQ2 to confirm coverage
      INTEGER(LONG)                   :: NUM_DOUBLINGS     ! Count of NEV doublings
      INTEGER(LONG)                   :: MAX_DOUBLINGS     ! Maximum number of doublings allowed
      INTEGER(LONG)                   :: INITIAL_NEV       ! Starting value for NEV
      INTEGER(LONG)                   :: MAX_NEV           ! Maximum NEV based on problem size and EIG_N2


      REAL(DOUBLE)                    :: EPS1              ! Small number for comparisons
      REAL(DOUBLE)                    :: SIGMA             ! Shift value (center of frequency range in omega^2)
      REAL(DOUBLE)                    :: OMEGA_FRQ1_SQ     ! (2*PI*EIG_FRQ1)^2
      REAL(DOUBLE)                    :: OMEGA_FRQ2_SQ     ! (2*PI*EIG_FRQ2)^2
      REAL(DOUBLE)                    :: FREQ              ! Temporary for frequency conversion
      REAL(DOUBLE)                    :: MIN_FREQ_FOUND    ! Minimum frequency found in current iteration
      REAL(DOUBLE)                    :: MAX_FREQ_FOUND    ! Maximum frequency found in current iteration
      REAL(DOUBLE)                    :: MIN_EIGEN_FOUND   ! Minimum eigenvalue found in current iteration
      REAL(DOUBLE)                    :: MAX_EIGEN_FOUND   ! Maximum eigenvalue found in current iteration
      REAL(DOUBLE)                    :: DELTA_SHIFT       ! Small shift for zero-lower-bound case
      INTEGER(LONG)                   :: NUM_DISCARDED     ! Number of eigenvalues discarded (outside range)
      INTEGER(LONG)                   :: IERR              ! Error return from LAPACK factorization
      REAL(DOUBLE)                    :: DUM_COL(1)        ! Dummy column for SuperLU deallocation
      INTEGER(LONG)                   :: MIN_NCV, MAX_NCV

      ! Temporary arrays for filtering
      REAL(DOUBLE), ALLOCATABLE       :: TEMP_EIGEN_VAL(:) ! Temporary eigenvalues for filtering
      REAL(DOUBLE), ALLOCATABLE       :: TEMP_EIGEN_VEC(:,:) ! Temporary eigenvectors for filtering
      INTEGER(LONG), ALLOCATABLE      :: TEMP_MODE_NUM(:)  ! Temporary mode numbers for filtering
      INTEGER(LONG), ALLOCATABLE      :: IN_RANGE_IDX(:)   ! Indices of eigenvalues in range

      INTRINSIC                       :: MIN, MAX, SQRT, ABS



! **********************************************************************************************************************************
! NOTE: This routine should only be called when:
!   - A frequency range is specified (EIG_FRQ2 > 0)
!   - SOL_NAME is not BUCKLING or GEN CB MODEL
! These conditions are checked in LINK4 before calling this routine.
! The checks below are defensive programming to catch any programming errors.

      EPS1 = EPSIL(1)

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9101) SUBR_NAME
         WRITE(F06,9101) SUBR_NAME
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9102) SUBR_NAME
         WRITE(F06,9102) SUBR_NAME
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (EIG_FRQ2 <= EPS1) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9103) SUBR_NAME
         WRITE(F06,9103) SUBR_NAME
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
! Initialize adaptive search parameters

      INITIAL_NEV = 10
      MAX_DOUBLINGS = 10
      NUM_DOUBLINGS = 0

      ! Calculate maximum NEV: MAX[?] of (NDOFL-1, EIG_N2 if specified, or NDOFL-1)
      MAX_NEV = NDOFL - NUM_MLL_DIAG_ZEROS - 1
      IF (EIG_N2 > 0) THEN
         MAX_NEV = MAX(MAX_NEV, EIG_N2)
      ENDIF
      MAX_NEV = MIN(MAX_NEV, INITIAL_NEV*(2**MAX_DOUBLINGS))

      ! Compute sigma (shift point in omega^2 = eigenvalue space)
      OMEGA_FRQ1_SQ = (TWO * PI * EIG_FRQ1)**2
      OMEGA_FRQ2_SQ = (TWO * PI * EIG_FRQ2)**2

      ! Determine if this is a zero-lower-bound (one-sided) search
      ZERO_LOWER_BOUND = (EIG_FRQ1 <= EPS1)

      IF (ZERO_LOWER_BOUND) THEN
         ! FRQ1 ~ 0: ONE-SIDED SEARCH targeting lowest eigenvalues
         ! Use small NEGATIVE sigma to avoid singularity if K has rigid-body modes
         ! and to make ARPACK's "largest magnitude" find smallest eigenvalues cleanly.
         ! Delta = (2*PI*1.0)^2 corresponds to ~1 Hz shift below zero
         DELTA_SHIFT = (TWO * PI * ONE)**2
         SIGMA = -DELTA_SHIFT
         PREV_IN_RANGE = 0
         MIN_ABOVE_FOR_STOP = 3  ! Require at least 3 modes above FRQ2 for confidence
         WRITE(F06,1033) EIG_FRQ1, EIG_FRQ2, DELTA_SHIFT
         IF (SUPINFO == 'N') THEN
            WRITE(SC1,1033) EIG_FRQ1, EIG_FRQ2, DELTA_SHIFT
         ENDIF
      ELSE
         ! Normal case: center of frequency range (two-sided search)
         SIGMA = (OMEGA_FRQ1_SQ + OMEGA_FRQ2_SQ) / TWO
         PREV_IN_RANGE = -1  ! Not used in two-sided mode
         MIN_ABOVE_FOR_STOP = 1
         WRITE(F06,1001) EIG_FRQ1, EIG_FRQ2, SQRT(SIGMA)/(TWO*PI)
         IF (SUPINFO == 'N') THEN
            WRITE(SC1,1001) EIG_FRQ1, EIG_FRQ2, SQRT(SIGMA)/(TWO*PI)
         ENDIF
      ENDIF

! **********************************************************************************************************************************
! Build KMSM = KLL - SIGMA*MLL (this is done ONCE since sigma is fixed)

      CALL LINK_MESSAGE('BUILD SHIFTED MATRIX [KLL - sigma*MLL]')

      CALL MATADD_SSS_NTERM ( NDOFL, 'KLL', NTERM_KLL, I_KLL, J_KLL, SYM_KLL, '-sigma*MLL',                                        &
                                            NTERM_MLL, I_MLL, J_MLL, SYM_MLL, 'KMSM', NTERM_KMSM )
      CALL ALLOCATE_SPARSE_MAT ( 'KMSM', NDOFL, NTERM_KMSM, SUBR_NAME )
      CALL MATADD_SSS ( NDOFL, 'KLL', NTERM_KLL, I_KLL, J_KLL, KLL, ONE, '-sigma*MLL',                                             &
                                      NTERM_MLL, I_MLL, J_MLL, MLL, -SIGMA,                                                        &
                               'KMSM', NTERM_KMSM, I_KMSM, J_KMSM, KMSM )

! Calculate bandwidth of KMSM
      CALL LINK_MESSAGE('CALCULATE BANDWIDTH OF [KLL - sigma*MLL]')
      CALL BANDSIZ ( NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM_SDIA )
      WRITE(ERR,4905) KMSM_SDIA
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4905) KMSM_SDIA
      ENDIF

! Determine LDRFAC based on matrix type
      IF (SOLLIB(1:6) == 'SPARSE') THEN
         LDRFAC = 1
      ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
         LDRFAC = KMSM_SDIA + 1
      ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DGB') THEN
         LDRFAC = 3*KMSM_SDIA + 1
      ELSE
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,4003) SUBR_NAME, EIG_LAP_MAT_TYPE
         WRITE(F06,4003) SUBR_NAME, EIG_LAP_MAT_TYPE
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Build KMSMn (nonsymmetric form) - needed for matrix-vector products in DSBAND
      IF (SPARSTOR == 'SYM   ') THEN
         CALL SPARSE_MAT_DIAG_ZEROS ( 'KMSM', NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, NUM_KMSM_DIAG_0 )
         NTERM_KMSMn = 2*NTERM_KMSM - (NDOFL - NUM_KMSM_DIAG_0)
         CALL ALLOCATE_SPARSE_MAT ( 'KMSMn', NDOFL, NTERM_KMSMn, SUBR_NAME )
         CALL CRS_SYM_TO_CRS_NONSYM ( 'KMSM', NDOFL, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, 'KMSMn', NTERM_KMSMn,                       &
                                       I_KMSMn, J_KMSMn, KMSMn, 'Y' )
      ELSE IF (SPARSTOR == 'NONSYM') THEN
         NTERM_KMSMn = NTERM_KMSM
         CALL ALLOCATE_SPARSE_MAT ( 'KMSMn', NDOFL, NTERM_KMSMn, SUBR_NAME )
         DO I=1,NDOFL+1
            I_KMSMn(I) = I_KMSM(I)
         ENDDO
         DO J=1,NTERM_KMSMn
            J_KMSMn(J) = J_KMSM(J)
            KMSMn(J) = KMSM(J)
         ENDDO
      ELSE
         WRITE(ERR,932) SUBR_NAME, SPARSTOR
         WRITE(F06,932) SUBR_NAME, SPARSTOR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Deallocate KLL (not needed anymore, and saves memory)
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLL', CR13
      CALL DEALLOCATE_SPARSE_MAT ( 'KLL' )

! **********************************************************************************************************************************
! FACTORIZATION (done ONCE before the adaptive loop)
! For BANDED: allocate RFAC, copy banded matrix, factor with dgbtrf/dpbtrf
! For SPARSE: call SYM_MAT_DECOMP_SUPRLU to factor using SuperLU

      KL = KMSM_SDIA
      KU = KL

      CALL LINK_MESSAGE('FACTOR SHIFTED MATRIX [KLL - sigma*MLL]')

      ! Allocate RFAC and IWORK (kept across all iterations)
      CALL ALLOCATE_LAPACK_MAT ( 'RFAC', LDRFAC, NDOFL, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'IWORK', NDOFL, 1, SUBR_NAME )

      IF (SOLLIB(1:6) == 'SPARSE') THEN
         ! Factor using SuperLU - factorization stored in SLU_FACTORS
         SLU_INFO = 0
         CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KMSM', 'L ',                                                                     &
                                      NDOFL, NTERM_KMSMn, I_KMSMn, J_KMSMn, KMSMn, SLU_INFO )
         IF (SLU_INFO /= 0) THEN
            WRITE(ERR,9903) SLU_INFO, SUBR_NAME
            WRITE(F06,9903) SLU_INFO, SUBR_NAME
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ELSE
         ! BANDED solver: copy KMSM to RFAC in band format, then factor
         IF (EIG_LAP_MAT_TYPE(1:3) == 'DPB') THEN
            CALL BANDGEN_LAPACK_DPB ( 'KMSM', NDOFL, KMSM_SDIA, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, RFAC, SUBR_NAME )
            IERR = 0
            CALL DPBTRF ( 'U', NDOFL, KU, RFAC, KU+1, IERR )
            DO I = 1, NDOFL
               IWORK(I) = I   ! No pivoting in Cholesky
            ENDDO
         ELSE IF (EIG_LAP_MAT_TYPE(1:3) == 'DGB') THEN
            CALL BANDGEN_LAPACK_DGB ( 'KMSM', NDOFL, KMSM_SDIA, NTERM_KMSM, I_KMSM, J_KMSM, KMSM, RFAC, SUBR_NAME )
            IERR = 0
            CALL DGBTRF ( NDOFL, NDOFL, KL, KU, RFAC, LDRFAC, IWORK, IERR )
         ENDIF
         IF (IERR /= 0) THEN
            CALL GET_GRID_AND_COMP ( 'A ', IERR, GRIDV, COMPV )
            WRITE(ERR,989) '[KLL - SIGMA*MLL]', 'DPBTRF/DGBTRF', IERR
            WRITE(F06,989) '[KLL - SIGMA*MLL]', 'DPBTRF/DGBTRF', IERR
            FATAL_ERR = FATAL_ERR + 1
            IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
               WRITE(ERR,9892) GRIDV, COMPV
               WRITE(F06,9892) GRIDV, COMPV
            ENDIF
            IF (BAILOUT >= 0) THEN
               CALL OUTA_HERE ( 'Y' )
            ENDIF
         ENDIF
      ENDIF

      WRITE(F06,1040)
      IF (SUPINFO == 'N') THEN
         WRITE(SC1,1040)
      ENDIF

! **********************************************************************************************************************************
! ADAPTIVE NEV LOOP
! Start with small NEV and double until we find eigenvalues outside [FRQ1, FRQ2] on both sides

      NEV = INITIAL_NEV
      SEARCH_COMPLETE = .FALSE.

      WRITE(SC1,1002)
      WRITE(SC1,1003) INITIAL_NEV, MAX_DOUBLINGS, MAX_NEV

      ADAPTIVE_LOOP: DO WHILE (.NOT. SEARCH_COMPLETE)

         ! Cap NEV at maximum
         IF (NEV > MAX_NEV) THEN
            NEV = MAX_NEV
         ENDIF

         WRITE(SC1,1004) NUM_DOUBLINGS + 1, NEV, MAX_NEV
         WRITE(SC1,1015) EIG_FRQ1, EIG_FRQ2
         WRITE(SC1,1016) OMEGA_FRQ1_SQ, OMEGA_FRQ2_SQ, SIGMA

         MIN_NCV = NEV + 2
         MAX_NCV = NDOFL - NUM_MLL_DIAG_ZEROS - 2

         ! sanity check on feasible NCV range
         IF (MIN_NCV > MAX_NCV) THEN
            WRITE(ERR,9776) NDOFL, NEV
            WRITE(F06,9776) NDOFL, NEV
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         END IF

         ! NCV determination step (try user factor first)
         DO I=EIG_NCVFACL,2,-1
            NCV = I*NEV
            IF (NCV >= MIN_NCV .AND. NCV <= MAX_NCV) THEN
               EXIT
            END IF
         END DO

         ! still invalid? (too large OR too small)
         IF (NCV > MAX_NCV .OR. NCV < MIN_NCV) THEN
            DO I=5,2,-1
               NCV = NEV+I
               IF (NCV >= MIN_NCV .AND. NCV <= MAX_NCV) THEN
                  EXIT
               END IF
            END DO
         END IF

         ! no valid NCV.
         IF (NCV > MAX_NCV .OR. NCV < MIN_NCV) THEN
            WRITE(ERR,9777) NDOFL, NEV
            WRITE(F06,9777) NDOFL, NEV
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         END IF

         ! Calculate NCV and workspace size
         !NCV = MIN(EIG_NCVFACL * NEV, NDOFL)
         !IF (NCV < NEV + 2) NCV = MIN(NEV + 2, NDOFL)
         LWORKL = NCV * (NCV + 8)

         KL = KMSM_SDIA
         KU = KL

         ! Set IPARAM for ARPACK
         DO I = 1, 11
            IPARAM(I) = 0
         ENDDO
         IPARAM(3) = MXITERL     ! Max iterations
         IPARAM(4) = 1           ! Block size (must be 1)
         IPARAM(7) = 3           ! Mode 3: shift-invert for generalized problem

         ! Set other DSBAND_PREFAC parameters
         RVEC = .TRUE.
         HOWMNY = 'A'
         BMAT = 'G'
         WHICH = 'LM'

         ! NOTE: RFAC and IWORK are already allocated and factored before the loop

         ! Allocate eigenvalue/eigenvector arrays
         CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC', NDOFL, NEV, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM', NDOFL, 1, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', NDOFL, 1, SUBR_NAME )

         ! Allocate ARPACK work arrays (these change size with NCV, so reallocate each iteration)
         CALL ALLOCATE_LAPACK_MAT ( 'RESID', NDOFL, 1, SUBR_NAME )
         CALL ALLOCATE_LAPACK_MAT ( 'SELECT', NCV, 1, SUBR_NAME )
         CALL ALLOCATE_LAPACK_MAT ( 'VBAS', NDOFL, NCV, SUBR_NAME )
         CALL ALLOCATE_LAPACK_MAT ( 'WORKD', 3*NDOFL, 1, SUBR_NAME )
         CALL ALLOCATE_LAPACK_MAT ( 'WORKL', LWORKL, 1, SUBR_NAME )

         DO I = 1, NCV
            SELECT(I) = .FALSE.
         ENDDO

         ! Call DSBAND_PREFAC to compute eigenvalues (uses pre-factored RFAC/SLU_FACTORS)
         CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - LANCZOS METH')

         INFO_ARPACK = 0
         INFO_LAPACK = 0

         CALL DSBAND_PREFAC ( RVEC, HOWMNY, SELECT, EIGEN_VAL, EIGEN_VEC, NDOFL, SIGMA, NDOFL, LDRFAC, RFAC, KL, KU, WHICH, BMAT,  &
                              NEV, ARP_TOL, RESID, NCV, VBAS, NDOFL, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO_ARPACK,             &
                              INFO_LAPACK, 'Y', DEBUG(47) )

         NUM_CONVERGED = IPARAM(5)

         WRITE(SC1,1005) NUM_CONVERGED

         ! Check for errors
         IF (INFO_ARPACK < 0) THEN
            WRITE(ERR,9996)
            WRITE(F06,9996)
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         IF (INFO_LAPACK < 0) THEN
            WRITE(ERR,993) SUBR_NAME, 'DPBTRF/DGBTRF'
            WRITE(F06,993) SUBR_NAME, 'DPBTRF/DGBTRF'
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ELSE IF (INFO_LAPACK > 0) THEN
            CALL GET_GRID_AND_COMP ( 'A ', INFO_LAPACK, GRIDV, COMPV )
            WRITE(ERR,989) '[KLL - SIGMA*MLL]', 'DPBTRF/DGBTRF', INFO_LAPACK
            WRITE(F06,989) '[KLL - SIGMA*MLL]', 'DPBTRF/DGBTRF', INFO_LAPACK
            FATAL_ERR = FATAL_ERR + 1
            IF ((GRIDV > 0) .AND. (COMPV > 0)) THEN
               WRITE(ERR,9892) GRIDV, COMPV
               WRITE(F06,9892) GRIDV, COMPV
            ENDIF
            IF (BAILOUT >= 0) THEN
               CALL OUTA_HERE ( 'Y' )
            ENDIF
         ENDIF

         ! Analyze results: find min and max frequencies and eigenvalues
         ! Also count modes in range and above FRQ2 for stopping logic
         IF (ZERO_LOWER_BOUND) THEN
            FOUND_BELOW_FRQ1 = .TRUE.  ! Automatically satisfied for one-sided search
            WRITE(SC1,1032)
         ELSE
            FOUND_BELOW_FRQ1 = .FALSE.
         ENDIF
         FOUND_ABOVE_FRQ2 = .FALSE.
         MIN_FREQ_FOUND = 1.0D30
         MAX_FREQ_FOUND = -1.0D30
         MIN_EIGEN_FOUND = 1.0D30
         MAX_EIGEN_FOUND = -1.0D30
         NUM_IN_RANGE = 0
         NUM_ABOVE_FRQ2 = 0

         WRITE(SC1,1017) NUM_CONVERGED

         DO I = 1, NUM_CONVERGED
            IF (EIGEN_VAL(I) >= ZERO) THEN
               FREQ = SQRT(EIGEN_VAL(I)) / (TWO * PI)
            ELSE
               ! Negative eigenvalue - treat as zero frequency (rigid body mode)
               FREQ = ZERO
            ENDIF

            IF (EIGEN_VAL(I) < MIN_EIGEN_FOUND) MIN_EIGEN_FOUND = EIGEN_VAL(I)
            IF (EIGEN_VAL(I) > MAX_EIGEN_FOUND) MAX_EIGEN_FOUND = EIGEN_VAL(I)
            IF (FREQ < MIN_FREQ_FOUND) MIN_FREQ_FOUND = FREQ
            IF (FREQ > MAX_FREQ_FOUND) MAX_FREQ_FOUND = FREQ

            IF (FREQ < EIG_FRQ1 - EPS1) FOUND_BELOW_FRQ1 = .TRUE.
            IF (FREQ > EIG_FRQ2 + EPS1) THEN
               FOUND_ABOVE_FRQ2 = .TRUE.
               NUM_ABOVE_FRQ2 = NUM_ABOVE_FRQ2 + 1
            ENDIF
            IF (FREQ >= EIG_FRQ1 - EPS1 .AND. FREQ <= EIG_FRQ2 + EPS1) THEN
               NUM_IN_RANGE = NUM_IN_RANGE + 1
            ENDIF
         ENDDO

         ! Detailed debug output
         WRITE(SC1,1018) MIN_EIGEN_FOUND, MIN_FREQ_FOUND
         WRITE(SC1,1019) MAX_EIGEN_FOUND, MAX_FREQ_FOUND
         WRITE(SC1,1020) OMEGA_FRQ1_SQ, EIG_FRQ1
         WRITE(SC1,1021) OMEGA_FRQ2_SQ, EIG_FRQ2
         WRITE(SC1,1022) FOUND_BELOW_FRQ1, FOUND_ABOVE_FRQ2
         WRITE(SC1,1034) NUM_IN_RANGE, NUM_ABOVE_FRQ2

         ! Check stopping conditions
         WRITE(SC1,1023)

         IF (ZERO_LOWER_BOUND) THEN
            ! ONE-SIDED SEARCH: Stop when we have modes above FRQ2 AND count stabilized
            IF (NUM_CONVERGED < NEV) THEN
               ! ARPACK exhausted - all modes found
               SEARCH_COMPLETE = .TRUE.
               WRITE(SC1,1024)
               WRITE(SC1,1009) NUM_CONVERGED, NEV
            ELSE IF (NEV >= MAX_NEV) THEN
               SEARCH_COMPLETE = .TRUE.
               WRITE(SC1,1025)
               WRITE(SC1,1010) MAX_NEV
            ELSE IF (NUM_DOUBLINGS >= MAX_DOUBLINGS) THEN
               SEARCH_COMPLETE = .TRUE.
               WRITE(SC1,1026)
               WRITE(SC1,1011) MAX_DOUBLINGS
            ELSE IF (NUM_ABOVE_FRQ2 >= MIN_ABOVE_FOR_STOP .AND. NUM_IN_RANGE == PREV_IN_RANGE) THEN
               ! Have enough modes above FRQ2 AND in-range count stabilized
               SEARCH_COMPLETE = .TRUE.
               WRITE(SC1,1035)
               WRITE(SC1,1036) NUM_ABOVE_FRQ2, MIN_ABOVE_FOR_STOP, NUM_IN_RANGE
            ELSE IF (NUM_ABOVE_FRQ2 >= MIN_ABOVE_FOR_STOP .AND. NUM_IN_RANGE > 0) THEN
               ! Have modes above FRQ2 but count still changing - continue one more time
               WRITE(SC1,1037)
               WRITE(SC1,1038) NUM_ABOVE_FRQ2, PREV_IN_RANGE, NUM_IN_RANGE
            ELSE
               WRITE(SC1,1027)
               IF (NUM_ABOVE_FRQ2 < MIN_ABOVE_FOR_STOP) WRITE(SC1,1039) NUM_ABOVE_FRQ2, MIN_ABOVE_FOR_STOP
            ENDIF
            PREV_IN_RANGE = NUM_IN_RANGE
         ELSE
            ! TWO-SIDED SEARCH: Original logic
            IF (FOUND_BELOW_FRQ1 .AND. FOUND_ABOVE_FRQ2) THEN
               ! Found eigenvalues on both sides of the range - we're done
               SEARCH_COMPLETE = .TRUE.
               WRITE(SC1,1008)
            ELSE IF (NUM_CONVERGED < NEV) THEN
               ! ARPACK found fewer eigenvalues than requested - all eigenvalues found
               SEARCH_COMPLETE = .TRUE.
               WRITE(SC1,1024)
               WRITE(SC1,1009) NUM_CONVERGED, NEV
            ELSE IF (NEV >= MAX_NEV) THEN
               ! Reached maximum NEV
               SEARCH_COMPLETE = .TRUE.
               WRITE(SC1,1025)
               WRITE(SC1,1010) MAX_NEV
            ELSE IF (NUM_DOUBLINGS >= MAX_DOUBLINGS) THEN
               ! Reached maximum number of doublings
               SEARCH_COMPLETE = .TRUE.
               WRITE(SC1,1026)
               WRITE(SC1,1011) MAX_DOUBLINGS
            ELSE
               ! Not stopping - explain why
               WRITE(SC1,1027)
               IF (.NOT. FOUND_BELOW_FRQ1) WRITE(SC1,1028) EIG_FRQ1
               IF (.NOT. FOUND_ABOVE_FRQ2) WRITE(SC1,1029) EIG_FRQ2
            ENDIF
         ENDIF

         IF (.NOT. SEARCH_COMPLETE) THEN
            ! Need to continue - deallocate work arrays (but NOT RFAC/IWORK) and double NEV
            ! RFAC and IWORK contain the factorization, kept across all iterations
            CALL DEALLOCATE_LAPACK_MAT ( 'RESID' )
            CALL DEALLOCATE_LAPACK_MAT ( 'SELECT' )
            CALL DEALLOCATE_LAPACK_MAT ( 'VBAS' )
            CALL DEALLOCATE_LAPACK_MAT ( 'WORKD' )
            CALL DEALLOCATE_LAPACK_MAT ( 'WORKL' )
            CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC' )
            CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM' )
            CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )

            NEV = NEV * 2
            NUM_DOUBLINGS = NUM_DOUBLINGS + 1
         ENDIF

      ENDDO ADAPTIVE_LOOP

! **********************************************************************************************************************************
! Filter results to keep only eigenvalues within [EIG_FRQ1, EIG_FRQ2]

      WRITE(SC1,1012)

      ! Count eigenvalues in range
      NUM_IN_RANGE = 0
      DO I = 1, NUM_CONVERGED
         IF (EIGEN_VAL(I) >= ZERO) THEN
            FREQ = SQRT(EIGEN_VAL(I)) / (TWO * PI)
         ELSE
            FREQ = ZERO
         ENDIF
         IF (FREQ >= EIG_FRQ1 - EPS1 .AND. FREQ <= EIG_FRQ2 + EPS1) THEN
            NUM_IN_RANGE = NUM_IN_RANGE + 1
         ENDIF
      ENDDO

      NUM_DISCARDED = NUM_CONVERGED - NUM_IN_RANGE
      WRITE(SC1,1030) NUM_CONVERGED, NUM_IN_RANGE, NUM_DISCARDED
      WRITE(SC1,1031) EIG_FRQ1, EIG_FRQ2

      IF (NUM_IN_RANGE == 0) THEN
         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,9104) EIG_FRQ1, EIG_FRQ2
         WRITE(F06,9104) EIG_FRQ1, EIG_FRQ2
         NUM_EIGENS = 0
         NVEC = 0
      ELSE
         ! Allocate temporary arrays for filtering
         ALLOCATE(TEMP_EIGEN_VAL(NUM_IN_RANGE))
         ALLOCATE(TEMP_EIGEN_VEC(NDOFL, NUM_IN_RANGE))
         ALLOCATE(TEMP_MODE_NUM(NUM_IN_RANGE))
         ALLOCATE(IN_RANGE_IDX(NUM_IN_RANGE))

         ! Identify indices of in-range eigenvalues
         K = 0
         DO I = 1, NUM_CONVERGED
            IF (EIGEN_VAL(I) >= ZERO) THEN
               FREQ = SQRT(EIGEN_VAL(I)) / (TWO * PI)
            ELSE
               FREQ = ZERO
            ENDIF
            IF (FREQ >= EIG_FRQ1 - EPS1 .AND. FREQ <= EIG_FRQ2 + EPS1) THEN
               K = K + 1
               IN_RANGE_IDX(K) = I
            ENDIF
         ENDDO

         ! Copy filtered results to temporary arrays
         DO K = 1, NUM_IN_RANGE
            I = IN_RANGE_IDX(K)
            TEMP_EIGEN_VAL(K) = EIGEN_VAL(I)
            TEMP_MODE_NUM(K) = K  ! Renumber modes sequentially
            DO J = 1, NDOFL
               TEMP_EIGEN_VEC(J, K) = EIGEN_VEC(J, I)
            ENDDO
         ENDDO

         ! Deallocate original arrays
         CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC' )
         CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM' )
         CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )

         ! Reallocate with correct size
         CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC', NDOFL, NUM_IN_RANGE, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM', NDOFL, 1, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', NDOFL, 1, SUBR_NAME )

         ! Copy filtered results back
         DO K = 1, NUM_IN_RANGE
            EIGEN_VAL(K) = TEMP_EIGEN_VAL(K)
            MODE_NUM(K) = TEMP_MODE_NUM(K)
            DO J = 1, NDOFL
               EIGEN_VEC(J, K) = TEMP_EIGEN_VEC(J, K)
            ENDDO
         ENDDO

         ! Free temporary arrays
         DEALLOCATE(TEMP_EIGEN_VAL)
         DEALLOCATE(TEMP_EIGEN_VEC)
         DEALLOCATE(TEMP_MODE_NUM)
         DEALLOCATE(IN_RANGE_IDX)

         NUM_EIGENS = NUM_IN_RANGE
         NVEC = NUM_IN_RANGE
      ENDIF

! **********************************************************************************************************************************
! Cleanup

! Free SuperLU factorization (for SPARSE solver)
      IF (SOLLIB(1:6) == 'SPARSE') THEN
         DUM_COL(1) = ZERO
         CALL C_FORTRAN_DGSSV ( 3, NDOFL, NTERM_KMSMn, 1, KMSMn, J_KMSMn, I_KMSMn, DUM_COL, NDOFL, SLU_FACTORS, SLU_INFO )
      ENDIF

      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KMSMn ', CR13
      CALL DEALLOCATE_SPARSE_MAT ( 'KMSMn' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KMSM  ', CR13
      CALL DEALLOCATE_SPARSE_MAT ( 'KMSM' )
      CALL DEALLOCATE_LAPACK_MAT ( 'IWORK' )
      CALL DEALLOCATE_LAPACK_MAT ( 'RFAC' )
      CALL DEALLOCATE_LAPACK_MAT ( 'RESID' )
      CALL DEALLOCATE_LAPACK_MAT ( 'SELECT' )
      CALL DEALLOCATE_LAPACK_MAT ( 'VBAS' )
      CALL DEALLOCATE_LAPACK_MAT ( 'WORKD' )
      CALL DEALLOCATE_LAPACK_MAT ( 'WORKL' )

      WRITE(SC1,1014) NUM_EIGENS


      RETURN

! **********************************************************************************************************************************
! Format statements

 1001 FORMAT(/,' *INFORMATION: ADAPTIVE LANCZOS EIGENVALUE SEARCH',                                                                &
             /,'               Frequency range: ',F12.4,' Hz to ',F12.4,' Hz',                                                     &
             /,'               Shift (sigma) at center frequency: ',F12.4,' Hz')

 1002 FORMAT(/,' *INFORMATION: BEGINNING ADAPTIVE NEV SEARCH')

 1003 FORMAT('               Initial NEV = ',I8,', Max doublings = ',I4,', Max NEV = ',I8)

 1004 FORMAT(/,'  ======================================================',                                         &
             /,'  ADAPTIVE ITERATION ',I3,': NEV = ',I8,' (max = ',I8,')',                                                         &
             /,'  ======================================================')

 1005 FORMAT('               ARPACK converged ',I8,' eigenvalues')

 1006 FORMAT('               Frequency range found: ',F12.4,' Hz to ',F12.4,' Hz')

 1007 FORMAT('               Found below FRQ1: ',L1,', Found above FRQ2: ',L1)

 1008 FORMAT('    --> STOP REASON: Both bounds satisfied (eigenvalues found outside [FRQ1,FRQ2] on both sides)')

 1009 FORMAT('                     ARPACK returned ',I8,' eigenvalues, requested ',I8)

 1010 FORMAT('                     Maximum NEV reached: ',I8)

 1011 FORMAT('                     Maximum doublings reached: ',I4)

 1012 FORMAT(/,'  ======================================================',                                         &
             /,'  FILTERING RESULTS',                                                                                              &
             /,'  ======================================================')

 1013 FORMAT('               Total converged: ',I8,', In range [',F10.3,', ',F10.3,' Hz]: ',I8)

 1014 FORMAT(/,'  ======================================================',                                         &
             /,'  ADAPTIVE LANCZOS COMPLETE',                                                                                      &
             /,'  Final eigenvalue count (in range): ',I8,                                                                         &
             /,'  ======================================================',/)

 1015 FORMAT('    Search frequency band:   [',ES14.6,' Hz, ',ES14.6,' Hz]')

 1016 FORMAT('    Search eigenvalue band:  [',ES14.6,', ',ES14.6,'] (omega^2)',                                                    &
             /,'    Shift sigma (center):    ',ES14.6,' (omega^2)')

 1017 FORMAT(/,'    ARPACK returned ',I8,' converged eigenvalues')

 1018 FORMAT('    Lowest  eigenvalue found: ',ES14.6,' (omega^2) = ',ES14.6,' Hz')

 1019 FORMAT('    Highest eigenvalue found: ',ES14.6,' (omega^2) = ',ES14.6,' Hz')

 1020 FORMAT('    Lower  bound (FRQ1):      ',ES14.6,' (omega^2) = ',ES14.6,' Hz')

 1021 FORMAT('    Upper  bound (FRQ2):      ',ES14.6,' (omega^2) = ',ES14.6,' Hz')

 1022 FORMAT('    Found eigenvalue below FRQ1? ',L1,'    Found eigenvalue above FRQ2? ',L1)

 1023 FORMAT(/,'    Checking stopping conditions...')

 1024 FORMAT('    --> STOP REASON: ARPACK exhausted (returned fewer eigenvalues than requested)')

 1025 FORMAT('    --> STOP REASON: Maximum NEV reached (early stop)')

 1026 FORMAT('    --> STOP REASON: Maximum doublings reached (early stop)')

 1027 FORMAT('    --> CONTINUING: Bounds not yet satisfied')

 1028 FORMAT('        - Still need eigenvalue BELOW FRQ1 = ',ES14.6,' Hz')

 1029 FORMAT('        - Still need eigenvalue ABOVE FRQ2 = ',ES14.6,' Hz')

 1030 FORMAT('    Total converged: ',I8,', In range: ',I8,', Discarded: ',I8)

 1031 FORMAT('    Keeping eigenvalues in [',ES14.6,' Hz, ',ES14.6,' Hz]')

 1032 FORMAT('    Lower bound FRQ1 ~ 0: ONE-SIDED SEARCH (lower bound auto-satisfied)')

 1033 FORMAT(/,' *INFORMATION: ADAPTIVE LANCZOS - ONE-SIDED SEARCH (FRQ1 ~ 0)',                                                    &
             /,'               Frequency range: ',F12.4,' Hz to ',F12.4,' Hz',                                                     &
             /,'               Using NEGATIVE sigma = -delta, delta = ',ES12.4,' (omega^2)',                                       &
             /,'               This targets lowest eigenvalues first to ensure complete coverage')

 1034 FORMAT('    Modes in range: ',I8,', Modes above FRQ2: ',I8)

 1035 FORMAT('    --> STOP REASON: One-sided search complete (count stabilized with modes above FRQ2)')

 1036 FORMAT('                     Modes above FRQ2: ',I4,' >= ',I4,' required, In-range count stable at ',I8)

 1037 FORMAT('    --> CONTINUING: Have modes above FRQ2 but in-range count still changing')

 1038 FORMAT('                     Modes above FRQ2: ',I4,', In-range: ',I8,' -> ',I8)

 1039 FORMAT('        - Need more modes above FRQ2: have ',I4,', need >= ',I4)

 1040 FORMAT(' *INFORMATION: MATRIX [KLL - sigma*MLL] FACTORIZATION COMPLETE (REUSED FOR ALL ITERATIONS)')

  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A,                                                                  &
                    /,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

  989 FORMAT(' *ERROR   989: THE FACTORIZATION OF THE MATRIX ',A,' COULD NOT BE COMPLETED BY LAPACK SUBR ',A,                      &
                    /,14X,' THE LEADING MINOR OF ORDER ',I12,' IS NOT POSITIVE DEFINITE')

  993 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A,                                                                  &
                    /,14X,' LAPACK SUBR XERBLA REPORTED AN ERROR IN A CALL TO ',A)

 4003 FORMAT(' *ERROR  4003: PROGRAMMING ERROR IN SUBROUTINE ',A,                                                                  &
                    /,14X,' EIG_LAP_MAT_TYPE MUST BE EITHER "DGB" OR "DPB" BUT IS = ',A)

 4905 FORMAT(' *INFORMATION: NUMBER OF SUPERDIAGONALS IN THE KMSM = [KLL - sigma*MLL] MATRIX UPPER TRIANGLE IS      = ',I12,/)

 9101 FORMAT(' *ERROR  9101: SUBROUTINE ',A,' DOES NOT SUPPORT BUCKLING SOLUTIONS.',                                               &
                    /,14X,' USE EIG_LANCZOS_ARPACK INSTEAD.')

 9102 FORMAT(' *ERROR  9102: SUBROUTINE ',A,' DOES NOT SUPPORT GEN CB MODEL SOLUTIONS.',                                           &
                    /,14X,' USE EIG_LANCZOS_ARPACK INSTEAD.')

 9103 FORMAT(' *ERROR  9103: SUBROUTINE ',A,' REQUIRES A FREQUENCY RANGE (EIG_FRQ2 > 0).',                                         &
                    /,14X,' USE EIG_LANCZOS_ARPACK FOR NON-FREQUENCY-RANGE SEARCHES.')

 9104 FORMAT(' *WARNING 9104: NO EIGENVALUES FOUND IN THE FREQUENCY RANGE ',F12.4,' Hz TO ',F12.4,' Hz')

 9892 FORMAT('               THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT ',I8,' COMPONENT ',I3)

 9903 FORMAT(' *ERROR  9903: SUPERLU SPARSE SOLVER HAS FAILED WITH INFO = ', I12,' IN SUBR ', A)

 9996 FORMAT('  PROCESSING STOPPED DUE TO ARPACK ERRORS')

 9776 FORMAT(' *ERROR  9776: TOO MANY EIGENVALUES REQUESTED FOR THIS PROBLEM SIZE: NDOFL=',I0,', NEV=',I0,'.')

 9777 FORMAT(' *ERROR  9777: UNABLE TO DETERMINE KRYLOV SUBSPACE SIZE NCV. NDOFL=',I0,', NEV=',I0,'.',/,15X,&
             'USER HAS LIKELY REQUESTED TOO MANY EIGENVALUES FOR THIS MODEL.',/,15X,&
             'NCV MUST BE AT LEAST NEV+2, BUT NO MORE THAN NDOFL-2.')


12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE EIG_LANCZOS_ARPACK_ADAPTIVE

   END MODULE LANCZOS
