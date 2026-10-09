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

   MODULE EIGEN_SUPPORT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CALC_GEN_MASS, RENORM_ON_MASS, EIG_SUMMARY

   CONTAINS

      SUBROUTINE CALC_GEN_MASS

! Generates generalized mass from mass matrix and eigenvectors:

!   The generalized mass matrix is a square matrix of terms:

!         MIJ = EIGEN_VEC(i)'*MLL*EIGEN_VEC(j)   where EIGEN_VEC(i) is the ith eigenvector and MLL is the L-set mass matrix
!                                                The ' indicates a transpose of EIGEN_VEC(i)

!   Array GEN_MASS is a 1-D array of the diagonal terms, MIJ (i = j) from the square generalized mass matrix outlined above.
!   This subr calculates all NDOFL diagonal terms of the generalized mass matrix plus all off diagonal terms below the diagonal.
!   The diagonal terms go into array GEN_MASS. The off diagonal terms are not stored, only the largest one, MAXMIJ is kept, and
!   output later, so that the user will know to what accuracy the eigenvectors were calculated


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NDOFL, NTERM_KLLDn, NTERM_MLLn, NVEC, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL
      USE EIGEN_MATRICES_1, ONLY      :  GEN_MASS, EIGEN_VEC
      USE MODEL_STUF, ONLY            :  EIG_CRIT, MAXMIJ, MIJ_COL, MIJ_ROW, NUM_FAIL_CRIT
      USE SPARSE_MATRICES, ONLY       :  I_KLLDn, J_KLLDn, KLLDn, I_MLLn, J_MLLn, MLLn
      USE SPARSE_MATRICES, ONLY       :  SYM_MLLn

      USE SPARSE_FULL_MULTIPLICATION, ONLY:  MATMULT_SFF
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_GEN_MASS'


      INTEGER(LONG)                   :: I,J,K             ! DO loop indices

      REAL(DOUBLE)                    :: DMIJ              ! DABS of MIJ
      REAL(DOUBLE)                    :: MAX               ! Temporary variable used in finding MAXMIJ
      REAL(DOUBLE)                    :: MIJ               ! The i,j-th value from gen. mass matrix. Used to find MAXMIJ
      REAL(DOUBLE)                    :: OUTVECI(NDOFL,1)  ! One eigenvector
      REAL(DOUBLE)                    :: OUTVECJ(NDOFL,1)  ! One eigenvector
      REAL(DOUBLE)                    :: ZVEC(NDOFL,1)     ! Intermediate matrix in the calculation of GEN_MASS

      REAL(DOUBLE), EXTERNAL          :: DDOT              ! BLAS dot-product function

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      NUM_FAIL_CRIT = 0
      MIJ_ROW       = 1
      MIJ_COL       = 1
      MAX           = ZERO
      MAXMIJ        = ZERO
      CALL COUNTER_INIT('     Diag term for eigenvector ', NVEC)
      DO I=1,NVEC

         DO K=1,NDOFL                                      ! Calc diag terms
            OUTVECI(K,1) = EIGEN_VEC(K,I)
         ENDDO

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            CALL MATMULT_SFF ( 'KLLDn', NDOFL, NDOFL, NTERM_KLLDn, 'N'     , I_KLLDn, J_KLLDn, KLLDn, 'OUTVECI', NDOFL, 1,         &
                                OUTVECI, 'N', 'ZVEC', ONE, ZVEC )
         ELSE
            CALL MATMULT_SFF ( 'MLLn' , NDOFL, NDOFL, NTERM_MLLn , SYM_MLLn, I_MLLn , J_MLLn , MLLn , 'OUTVECI', NDOFL, 1,         &
                                OUTVECI, 'N', 'ZVEC', ONE, ZVEC )
         ENDIF

         GEN_MASS(I) = DDOT ( NDOFL, OUTVECI, 1, ZVEC, 1 )
         GEN_MASS(I) = ABS(GEN_MASS(I))
         IF (DEBUG(48) == 0) THEN                          ! Calc off-diag terms

            !CALL COUNTER_INIT('Off-diag term ', I-1)
            DO J=1,I-1

               DO K=1,NDOFL
                  OUTVECJ(K,1) = EIGEN_VEC(K,J)
               ENDDO

               IF (SOL_NAME(1:8) == 'BUCKLING') THEN
                  CALL MATMULT_SFF ( 'KLLDn', NDOFL, NDOFL, NTERM_KLLDn, 'N'     , I_KLLDn, J_KLLDn, KLLDn, 'OUTVECJ', NDOFL, 1,   &
                                      OUTVECJ, 'N', 'ZVEC', ONE, ZVEC )
               ELSE
                  CALL MATMULT_SFF ( 'MLLn' , NDOFL, NDOFL, NTERM_MLLn , SYM_MLLn, I_MLLn , J_MLLn , MLLn , 'OUTVECJ', NDOFL, 1,   &
                                      OUTVECJ,'N', 'ZVEC', ONE, ZVEC )
               ENDIF

               MIJ = DDOT ( NDOFL, OUTVECI, 1, ZVEC, 1 )

               DMIJ = DABS(MIJ)
               IF (DMIJ > MAX) THEN
                  MAXMIJ  = MIJ
                  MAX     = DMIJ
                  MIJ_ROW = I
                  MIJ_COL = J
               ENDIF
               IF (DMIJ > EIG_CRIT) THEN
                  NUM_FAIL_CRIT = NUM_FAIL_CRIT + 1
               ENDIF
               !CALL COUNTER_PROGRESS(J)
            ENDDO

         ENDIF
         CALL COUNTER_PROGRESS(I)
      ENDDO



      RETURN

! **********************************************************************************************************************************


! **********************************************************************************************************************************

      END SUBROUTINE CALC_GEN_MASS


      SUBROUTINE RENORM_ON_MASS ( NVC, EPS1 )

! Renormalizes eigenvectors to unit generalized mass

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  NDOFL, BLNK_SUB_NAM, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  EPSIL, SUPINFO, SUPWARN
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE EIGEN_MATRICES_1 , ONLY     :  GEN_MASS, EIGEN_VEC
      USE MODEL_STUF, ONLY            :  EIG_NORM, MAXMIJ, MIJ_COL, MIJ_ROW
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RENORM_ON_MASS'

      INTEGER(LONG), INTENT(IN)       :: NVC               ! Number of eigenvectors to be renormalized.
      INTEGER(LONG)                   :: I,J               ! DO loop index


      REAL(DOUBLE) , INTENT(IN)       :: EPS1              ! Small number to compare variables against zero
      REAL(DOUBLE)                    :: DEN               ! Normalizing factor in gen mass matrix normalization

      INTRINSIC DSQRT,DABS



! **********************************************************************************************************************************
      IF (EIG_NORM /= 'MASS    ') THEN
         WRITE(ERR,1001) EIG_NORM
         IF (SUPINFO == 'N') THEN
            WRITE(F06,1001) EIG_NORM
         ENDIF
      ENDIF

      DO I=1,NVC
         IF (DABS(GEN_MASS(I)) < EPS1) THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,4301) I, GEN_MASS(I)
            IF (SUPWARN == 'N') THEN
               WRITE(F06,4301) I, GEN_MASS(I)
            ENDIF
            RETURN
         ENDIF
      ENDDO

! Adjust MAXMIJ, the largest off-diag gen mass term. It was originally calculated in subr CALC_GEN_MASS and will change if the
! gen masses have changed as a result of this renormalization

      MAXMIJ = MAXMIJ/(GEN_MASS(MIJ_ROW)*GEN_MASS(MIJ_COL))! NOTE: all gen mass terms checked above for > 0.

! Normalize the eigenvectors so that they produce unit generalized mass and reset the gen masses to unity

      DO J=1,NVC
         DEN = DSQRT(GEN_MASS(J))
         DO I=1,NDOFL
            EIGEN_VEC(I,J) = EIGEN_VEC(I,J)/DEN
         ENDDO
         GEN_MASS(J) = ONE                                 ! Now reset generalized masses to unity
      ENDDO



      RETURN

! **********************************************************************************************************************************
 1001 FORMAT(' *INFORMATION: EIGENVECTORS WILL BE RENORMALIZED BASED ON GEN MASS IN LINK4. THEY WILL BE RENORMALIZED TO ',         &
                             A,' LATER IN LINK5',/)

 4301 FORMAT(' *WARNING    : THE GENERALIZED MASS MATRIX HAS A DIAGONAL TERM THAT IS TOO SMALL TO ALLOW RENORMALIZATION OF THE',   &
                           ' EIGENVECTORS.'                                                                                        &
                    ,/,14X,' THE SMALL TERM IS FOR EIGENVECTOR ',I8,' AND ITS VALUE IS ',1ES9.2                                    &
                    ,/,14X,' EIGENVECTORS WILL NOT BE RENORMALIZED TO UNIT MASS IN SUBR RENORM_ON_MASS')

99001 FORMAT(1X,10(1ES13.6))

! **********************************************************************************************************************************

      END SUBROUTINE RENORM_ON_MASS


      SUBROUTINE EIG_SUMMARY ( ISUB )

! Prints eigenvalue analysis summary table

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NDOFL, NUM_EIGENS, NVEC, NUM_KLLD_DIAG_ZEROS, NUM_MLL_DIAG_ZEROS, NUM_SUBC_CARDS, SOL_NAME, &
                                         WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  ART_MASS, ART_ROT_MASS, ART_TRAN_MASS, DARPACK, SOLLIB, SUPINFO, SUPWARN
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, PI
      USE EIGEN_MATRICES_1, ONLY      :  GEN_MASS, MODE_NUM, EIGEN_VAL
      USE MODEL_STUF, ONLY            :  EIG_COMP, EIG_CRIT, EIG_GRID, EIG_LAP_MAT_TYPE, EIG_METH, EIG_MODE, EIG_N2, EIG_NORM,     &
                                         EIG_SIGMA, LABEL, MAXMIJ, MIJ_COL, MIJ_ROW, NUM_FAIL_CRIT, SCNUM, STITLE, TITLE

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: ISUB              ! Internal subcase index (for header and buckling load factor label)

      LOGICAL                         :: FILE_OPND         ! .TRUE. if a file is opened

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EIG_SUMMARY'
      CHARACTER( 1*BYTE)              :: ASTERISK = '*'    ! Used for denoting negative eigenvalues

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: NUM_FINITE_EIGENS ! Number of eigenvalues that are finite (excluding zero mass modes)
      INTEGER(LONG)                   :: MAX_LANCZOS_EIGENS! Max number of eigenvalues that can be found by Lanczos method
      INTEGER(LONG)                   :: NUM_NEG_EIGENS    ! Number of eigenvalues that are negative
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN


      REAL(DOUBLE)                    :: CYCLES1           ! Circular frequency of a mode
      REAL(DOUBLE)                    :: GEN_STIFF1        ! Generalized stiffness for a mode
      REAL(DOUBLE)                    :: RADS1             ! Radian frequency of a mode



! **********************************************************************************************************************************
      OUNT(1) = ERR
      OUNT(2) = F06

      IF (NUM_SUBC_CARDS > 0)           WRITE(F06,9001) SCNUM(ISUB)
      IF (TITLE(ISUB)(1:)  /= ' ')      WRITE(F06,9011) TITLE(ISUB)
      IF (STITLE(ISUB)(1:) /= ' ')      WRITE(F06,9011) STITLE(ISUB)
      IF (LABEL(ISUB)(1:)  /= ' ')      WRITE(F06,9011) LABEL(ISUB)
      WRITE(F06,*)

      IF (EIG_METH == 'LANCZOS') THEN

         IF      (SOLLIB == 'BANDED  ') THEN
            WRITE(F06,90001) EIG_METH, EIG_MODE, TRIM(EIG_LAP_MAT_TYPE), EIG_SIGMA, '(BANDED solution)'
         ELSE IF (SOLLIB == 'SPARSE  ') THEN
            WRITE(F06,90001) EIG_METH, EIG_MODE, TRIM(EIG_LAP_MAT_TYPE), EIG_SIGMA, '(SPARSE solution)'
         ENDIF

      ELSE

         IF      (SOLLIB == 'BANDED  ') THEN
            WRITE(F06,90003) EIG_METH, '(BANDED solution)'
         ELSE IF (SOLLIB == 'SPARSE  ') THEN
            WRITE(F06,90003) EIG_METH, '(SPARSE solution)'
         ENDIF

      ENDIF
      WRITE(F06,90004) NUM_EIGENS

      IF (NVEC > 0) THEN

         IF (DEBUG(48) == 0) THEN                          ! Off diag terms were calculated, so write results

            IF (EIG_NORM == 'MASS') THEN
               WRITE(F06,91001) MAXMIJ
            ELSE IF (EIG_NORM == 'MAX') THEN
               WRITE(F06,91002) MAXMIJ
            ELSE IF (EIG_NORM == 'POINT') THEN
               WRITE(F06,91003) MAXMIJ,EIG_GRID,EIG_COMP
            ELSE IF (EIG_NORM == 'NONE') THEN
               WRITE(F06,91004) MAXMIJ
            ENDIF

            WRITE(F06,92004) MIJ_ROW
            WRITE(F06,92005)
            WRITE(F06,92006) MIJ_COL
            WRITE(F06,92007)
            WRITE(F06,92008) EIG_CRIT, NUM_FAIL_CRIT

         ELSE

            WRITE(F06,91005) DEBUG(48)

         ENDIF

         WRITE(F06,*)
         WRITE(F06,*)

      ELSE

         WRITE(F06,92010)

      ENDIF

      IF (NUM_SUBC_CARDS > 0)           WRITE(F06,9001) SCNUM(ISUB)
      IF (TITLE(ISUB)(1:)  /= ' ')      WRITE(F06,9011) TITLE(ISUB)
      IF (STITLE(ISUB)(1:) /= ' ')      WRITE(F06,9011) STITLE(ISUB)
      IF (LABEL(ISUB)(1:)  /= ' ')      WRITE(F06,9011) LABEL(ISUB)
      WRITE(F06,*)
      WRITE(F06,*)

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         WRITE(F06,94101) SCNUM(ISUB)
         WRITE(F06,94102)
      ELSE
         WRITE(F06,94201)
         WRITE(F06,94202)
      ENDIF

      NUM_NEG_EIGENS = 0
      DO I=1,NUM_EIGENS
         RADS1      = DSQRT(DABS(EIGEN_VAL(I)))
         CYCLES1    = RADS1/(TWO*PI)
         GEN_STIFF1 = EIGEN_VAL(I)*GEN_MASS(I)
         IF (EIGEN_VAL(I) < ZERO) THEN
            NUM_NEG_EIGENS = NUM_NEG_EIGENS + 1
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,95301) MODE_NUM(I),I,EIGEN_VAL(I),ASTERISK
            ELSE
               WRITE(F06,95302) MODE_NUM(I),I,EIGEN_VAL(I),ASTERISK,RADS1,CYCLES1,GEN_MASS(I),GEN_STIFF1
            ENDIF
         ELSE
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,95401) MODE_NUM(I),I,EIGEN_VAL(I)
            ELSE
               WRITE(F06,95402) MODE_NUM(I),I,EIGEN_VAL(I),         RADS1,CYCLES1,GEN_MASS(I),GEN_STIFF1
            ENDIF
         ENDIF
      ENDDO

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         NUM_FINITE_EIGENS = NDOFL - NUM_KLLD_DIAG_ZEROS
      ELSE
         NUM_FINITE_EIGENS = NDOFL - NUM_MLL_DIAG_ZEROS
      ENDIF

      IF  (EIG_N2 > NUM_FINITE_EIGENS) THEN

         WRITE(F06,*)
         WARN_ERR = WARN_ERR + 1

         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            WRITE(ERR,98006) EIG_N2, NUM_FINITE_EIGENS, NUM_KLLD_DIAG_ZEROS
            IF (SUPWARN == 'N') THEN
               WRITE(F06,98006) EIG_N2, NUM_FINITE_EIGENS, NUM_MLL_DIAG_ZEROS
            ENDIF
         ELSE
            WRITE(ERR,98006) EIG_N2, NUM_FINITE_EIGENS, NUM_MLL_DIAG_ZEROS
            IF (SUPWARN == 'N') THEN
               WRITE(F06,98006) EIG_N2, NUM_FINITE_EIGENS, NUM_MLL_DIAG_ZEROS
            ENDIF
         ENDIF

      ELSE

         IF (ART_MASS == 'Y') THEN
            WRITE(F06,*)
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,98007) ART_TRAN_MASS, ART_ROT_MASS
            IF (SUPWARN == 'N') THEN
               WRITE(F06,98007) ART_TRAN_MASS, ART_ROT_MASS
            ENDIF
         ENDIF

      ENDIF

      IF ((EIG_METH == 'LANCZOS') .AND. (DEBUG(185) == 0)) THEN
         MAX_LANCZOS_EIGENS = NUM_FINITE_EIGENS - 1 - DARPACK
         IF (EIG_N2 > MAX_LANCZOS_EIGENS) THEN
            WRITE(F06,*)
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,98008) NUM_EIGENS, NUM_FINITE_EIGENS, DARPACK
            IF (SUPWARN == 'N') THEN
               WRITE(F06,98008) NUM_EIGENS, NUM_FINITE_EIGENS, DARPACK
            ENDIF
         ENDIF
      ENDIF

      IF (NUM_NEG_EIGENS > 0) THEN
         WRITE(F06,*)
         WRITE(F06,99000) NUM_NEG_EIGENS
      ENDIF
      WRITE(F06,*)
      WRITE(F06,*)




      RETURN

! **********************************************************************************************************************************
99001 FORMAT(A1)
 9001 FORMAT(' OUTPUT FOR SUBCASE ',I8)
 9011 FORMAT(1X,A)

90001 FORMAT(/,27X,'E I G E N V A L U E   A N A L Y S I S   S U M M A R Y',3X,'(',A8,' Mode',I2,1X,A,', Shift eigen = ',1ES9.2,')',&
             /,70X,A,/)

90003 FORMAT(/,27X,'E I G E N V A L U E   A N A L Y S I S   S U M M A R Y',3X,'(',A8,')',/,A,/)

90004 FORMAT(32X,'NUMBER OF EIGENVALUES EXTRACTED  . . . . . .',2X,I8,/)

91001 FORMAT(32X,'LARGEST OFF-DIAGONAL GENERALIZED MASS TERM  ',1ES10.1,' (Vecs renormed to 1.0 for gen masses)',/)

91002 FORMAT(32X,'LARGEST OFF-DIAGONAL GENERALIZED MASS TERM  ',1ES10.1,' (Vecs renormed to 1.0 for max value)',/)

91003 FORMAT(32X,'LARGEST OFF-DIAGONAL GENERALIZED MASS TERM  ',1ES10.1,' (Vecs renormed to 1.0 at grid-comp',I8,'-',I1,')',/)

91004 FORMAT(32X,'LARGEST OFF-DIAGONAL GENERALIZED MASS TERM  ',1ES10.1,' (Vecs not renormalized)',/)

91005 FORMAT(32X,'OFF-DIAGONAL GENERALIZED MASS TERMS NOT CALCULATED BASED ON DEBUG(48) = ',I8)

92004 FORMAT(32X,'                                       . . .',2X,I8)

92005 FORMAT(32X,'          MODE PAIR . . . . . . . . . .')

92006 FORMAT(32X,'                                       . . .',2X,I8,/)

92007 FORMAT(32X,'NUMBER OF OFF DIAGONAL GENERALIZED MASS')

92008 FORMAT(32X,'TERMS FAILING CRITERION OF ',1ES8.1,'. . . . .',2X,I8,/)

92010 FORMAT(4X,'NO EIGENVECTORS WERE REQUESTED TO BE OUTPUT, SO NO GENERALIZED MASS OR GENERALIZED STIFFNESS HAS BEEN CALCULATED',&
             /)

94101 FORMAT(44X,'R E A L   E I G E N V A L U E S',/,38X,'(subcase',I8,' buckling load factors)',/)

94102 FORMAT(40X,' MODE  EXTRACTION      EIGENVALUE',/,40X,'NUMBER   ORDER',/)

94201 FORMAT(44X,'R E A L   E I G E N V A L U E S')

94202 FORMAT(3X,' MODE  EXTRACTION      EIGENVALUE           RADIANS              CYCLES            GENERALIZED         GENERALIZED&
&        ',/,3X,'NUMBER   ORDER                                                                        MASS              STIFFNESS'&
          ,/)

95301 FORMAT(38X,2I8,1ES20.6,A)

95302 FORMAT(1X,2I8,1ES20.6,A,1ES19.6,3(1ES20.6))

95401 FORMAT(38X,2I8,1ES20.6)

95402 FORMAT(1X,2I8,5(1ES20.6))

98006 FORMAT(' *WARNING    : THE BULK DATA EIGR/EIGRL ENTRY ASKED FOR MODES UP TO NUMBER',I8,'. HOWEVER, THIS MODEL HAS ONLY',I8   &
                    ,/,14X,' FINITE EIGENVALUES DUE TO THE FACT THAT THE L-SET MASS MATRIX HAS',I8,' ZERO MASS DEGREES OF FREEDOM.'&
                    ,/,14x,' (USE OF BULK DATA PARAM ART_MASS WITH SMALL VALUE MAY HELP TO AVOID EIGENVALUES THAT ARE',      &
                           ' THEORETICALLY INFINITE)')

98007 FORMAT(' *WARNING    : THE USE OF BULK DATA PARAM ART_MASS MAY HAVE MASKED SOME OTHERWISE INFINITE EIGENVALUES. IF A SMALL', &
                           ' VALUE FOR'                                                                                            &
                    ,/,14X,' ART_MASS WAS USED THERE MAY BE SOME LARGE EIGENVALUES AND, IF SO, THESE SHOULD BE IGNORED.',          &
                           ' IF A LARGE VALUE FOR'                                                                                 &
                    ,/,14X,' ART_MASS WAS USED THEN ALL EIGENVALUES MAY BE SUSPECT. THIS RUN USED THE FOLLOWING ART_MASS VALUES:',/&
                    ,/,14X,'                   Artificial mass for translational G-set DOFs = ',1ES14.6                            &
                    ,/,14X,'                   Artificial mass for rotational    G-set DOFs = ',1ES14.6)

98008 FORMAT(' *WARNING    : LANCZOS ONLY FOUND',I8,' EIGENVALUES DUE TO:'                                                         &
                    ,/,14X,' (1) IT CAN NEVER FIND ALL',I8,' FINITE EIGENVALUES (IT IS ONLY POSSIBLE TO FIND 1 LESS THAN ALL)'     &
                    ,/,14X,' (2) IT ALSO DOES NOT PRINT THE ',I2,' HIGHEST MODES (DEFAULT VALUE FOR BULK DATA PARAM DARPACK) DUE', &
                               ' TO ROUNDOFF ERROR.'                                                                               &
                    ,/,18X,    ' THE USER CAN CHANGE THIS NUMBER VIA BULK DATA PARAM DARPACK.')


99000 FORMAT(' *INFORMATION: ASTERISK INDICATES ' ,I8, ' NEGATIVE EIGENVALUES. MYSTRAN TOOK THE ABSOLUTE VALUES OF THESE TO CALC', &
                                  ' RADIAN FREQUENCY.',/,                                                                          &
                    '               IF NEGATIVE EIGENVALUES ARE NOT SMALL THEY MAY BE IN ERROR',/,                                 &
                    '               (NOTE THAT LARGE EIGENVALUE MAGNITUDES IN MGIV WILL RESULT IF THERE ARE MASSLESS DOF''s IN THE'&
                                  ,' A-SET)')

! **********************************************************************************************************************************

      END SUBROUTINE EIG_SUMMARY

   END MODULE EIGEN_SUPPORT
