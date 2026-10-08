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

   MODULE LINK3_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: LINK3

   CONTAINS

      SUBROUTINE LINK3

! LINK 3 solves the equation KLL*UL = PL where KLL, UL, PL are the L-set stiffness matrix, displs and loads. It solves the equation
! using one of three methods. For each method the solution is obtained in a 2 step process: (1) the KLL matrix is decomposed into
! triangular factors and (2) UL is solved for by forward-backward substitution (FBS). The 3 methods are:

!   a) The LAPACK freeware code. This code requires KLL to be in banded (NOT sparse) form. LAPACH has the advantage that
!      MYSTRAN contains the LAPACK source code so debugging is easy. Its disadvantage is that banded matrices require much more
!      memory than sparse storage for large stiffness matrices.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06, L3A, SC1, LINK3A, L3A_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, COMM, FATAL_ERR, KLL_SDIA, LINKNO, MBUG, NDOFL, NSUB,                       &
                                         NTERM_KLL, NTERM_PL, RESTART,  SOL_NAME, WARN_ERR
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, TEN
      USE PARAMS, ONLY                :  CRS_CCS, EPSERR, EPSIL, KLLRAT, RELINK3, RCONDK, SOLLIB, SUPWARN, SPARSE_FLAVOR
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL, I_PL, J_PL, PL
      USE LAPACK_DPB_MATRICES, ONLY   :  RES
      USE COL_VECS, ONLY              :  UL_COL, PL_COL
      USE MACHINE_PARAMS, ONLY        :  MACH_EPS, MACH_SFMIN
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE SCRATCH_MATRICES, ONLY      :  I_CCS1, J_CCS1, CCS1
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, SLU_INFO, SLU_DIAG_RATIO

! Interface module not needed for subr's DPBTRF and DPBTRS. These are "CONTAIN'ed" in module LAPACK_LIN_EQN_DPB,
! which is "USE'd" above

      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_ARRAY_ALLOC_STAT, LINK_MESSAGE, LINK_MESSAGE_I, WRITE_ALLOC_MEM_TABLE

      USE SPARSE_CRS_ACCESS, ONLY     :  GET_SPARSE_CRS_COL
      USE LAPACK_ADAPTERS, ONLY       :  FBS_LAPACK, SYM_MAT_DECOMP_LAPACK
      USE SUPERLU_ADAPTERS, ONLY      :  FBS_SUPRLU, SYM_MAT_DECOMP_SUPRLU
      USE DATE_TIME_UTILS, ONLY       :  OURDAT, OURTIM, TIME_INIT
      USE COL_VEC_LIFECYCLE, ONLY     :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
      USE LAPACK_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_LAPACK_MAT, DEALLOCATE_LAPACK_MAT
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_INQUIRE, FILE_OPEN
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_VECTOR
      USE TEMP_FILE_READERS, ONLY     :  READ_L1A
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE, WRITE_L1A
      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK3'
      CHARACTER(  2*BYTE)             :: L_SET    = 'L '   ! L-set designator
      CHARACTER(  1*BYTE)             :: EQUED             ! 'Y' if the stiff matrix was equilibrated in subr EQUILIBRATE
      CHARACTER(  1*BYTE)             :: NULL_COL          ! 'Y' if a col of KAO(transpose) is null

      INTEGER(LONG)                   :: DEB_PRT(2)        ! Debug numbers to say whether to write ABAND and/or its decomp to output
!                                                            file in called subr SYM_MAT_DECOMP_LAPACK (ABAND = band form of KLL)

      INTEGER(LONG)                   :: IER_DECOMP        ! Overall error indicator
      INTEGER(LONG)                   :: ISUB              ! DO loop index for subcases
      INTEGER(LONG)                   :: INFO     = 0      ! Info output from some routine that has been called
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG), PARAMETER        :: P_LINKNO = 2      ! Prior LINK no's that should have run before this LINK can execute

      REAL(DOUBLE)                    :: BETA              ! Multiple for rhs for use in subr FBS
      REAL(DOUBLE)                    :: DEN               ! K_INORM*UL_INORM + PL_INORM
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      REAL(DOUBLE)                    :: EQUIL_SCALE_FACS(NDOFL)
                                                           ! LAPACK_S values returned from subr SYM_MAT_DECOMP_LAPACK
      REAL(DOUBLE)                    :: DUM_COL(NDOFL)    ! Temp variable used in SuperLU
      REAL(DOUBLE)                    :: K_INORM           ! Inf norm of KLL matrix (det in  subr COND_NUM)
      REAL(DOUBLE)                    :: LAP_ERR1          ! Bound on displ error = 2*OMEGAI/RCOND
      REAL(DOUBLE)                    :: OMEGAI            ! RES_INORM/DEN (similar to EPSILON)
      REAL(DOUBLE)                    :: OMEGAI0           ! Upper bound on OMEGAI. OMEGAI0 = 10*NDOFL*MACH_EPS
      REAL(DOUBLE)                    :: PL_INORM          ! Inf norm of load vector
      REAL(DOUBLE)                    :: RES_INORM         ! Inf norm of residual vector R = K*UL - PL
      REAL(DOUBLE)                    :: RCOND             ! Recrip of cond no. of the KLL. Det in  subr COND_NUM
      REAL(DOUBLE)                    :: UL_INORM          ! Inf norm of displacement vector

      INTRINSIC                       :: DABS

!***********************************************************************************************************************************
      LINKNO = 3

      EPS1 = EPSIL(1)

! Set time initializing parameters

      CALL TIME_INIT

! Initialize WRT_BUG

      DO I=0,MBUG-1
         WRT_BUG(I) = 0
      ENDDO

! Get date and time, write to screen

      CALL OURDAT
      CALL OURTIM
      WRITE(SC1,152) LINKNO

! Make units for writing errors the screen until we open output files

      OUNT(1) = SC1
      OUNT(2) = SC1

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Write info to text files

      WRITE(F06,150) LINKNO
      WRITE(ERR,150) LINKNO

! Read LINK1A file

      CALL READ_L1A ( 'KEEP' )

! Check COMM for successful completion of prior LINKs

      IF (COMM(P_LINKNO) /= 'C') THEN
         WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
         WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Prior LINK's didn't complete, so quit
      ENDIF

! Make sure SOL is STATICS, BUCKLING or NLSTATIC

      IF ((SOL_NAME(1:7) /= 'STATICS') .AND. (SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC')) THEN
         WRITE(ERR,999) SOL_NAME, 'STATICS or BUCKLING or NLSTATIC'
         WRITE(F06,999) SOL_NAME, 'STATICS or BUCKLING or NLSTATIC'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

!***********************************************************************************************************************************
! Factor KLL

      DEB_PRT(1) = 34
      DEB_PRT(2) = 35
      IER_DECOMP = 0

      DO J=1,NDOFL                                         ! Need a null col of loads when SuperLU is called to factor KLL
         DUM_COL(J) = ZERO                                 ! (only because it appears in the calling list)
      ENDDO

      IF ((RESTART == 'Y') .AND. (RELINK3 == 'Y')) THEN
sol_do:  DO
            WRITE(SC1,*) ' Input the value of SOLLIB (8 characters) to use in this restart:'
            READ (*,*) SOLLIB
            IF ((SOLLIB /= 'BANDED  ') .AND. (SOLLIB /= 'SPARSE  ')) THEN
               WRITE(SC1,*) '  Incorrect SOLLIB. Value must be BANDED or SPARSE'
               WRITE(SC1,*)
               CYCLE sol_do
            ELSE
               EXIT sol_do
            ENDIF
         ENDDO sol_do
      ENDIF

Factr:IF (SOLLIB == 'BANDED  ') THEN                       ! Use LAPACK

         INFO = 0
         CALL SYM_MAT_DECOMP_LAPACK ( SUBR_NAME, 'KLL', L_SET, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, 'Y', KLLRAT, 'Y', RCONDK,      &
                                      DEB_PRT, EQUED, KLL_SDIA, K_INORM, RCOND, EQUIL_SCALE_FACS, INFO )

      ELSE IF (SOLLIB == 'SPARSE  ') THEN

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            SLU_INFO = 0
            SLU_DIAG_RATIO = 'Y'                           ! Check matrix diag / factor diag (MAXRATIO), as on the LAPACK path
            CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KLL', L_SET, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, SLU_INFO )
            SLU_DIAG_RATIO = 'N'

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

      ENDIF Factr

!***********************************************************************************************************************************
!  Allocate col vector arrays for loads, displs and res vector

!xx   CALL ALLOCATE_COL_VEC ( 'UL_COL', NDOFL, SUBR_NAME )
!xx   CALL ALLOCATE_COL_VEC ( 'PL_COL', NDOFL, SUBR_NAME )
      CALL ALLOCATE_LAPACK_MAT ( 'RES', NDOFL, 1, SUBR_NAME )

! Open file for writing displs to.

      CALL FILE_OPEN ( L3A, LINK3A, OUNT, 'REPLACE', L3A_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

! Loop on subcases

      WRITE(F06,*)
      BETA = ONE
Solve:DO ISUB = 1,NSUB

         SLU_INFO = 0
         CALL ALLOCATE_COL_VEC ( 'UL_COL', NDOFL, SUBR_NAME )
         CALL ALLOCATE_COL_VEC ( 'PL_COL', NDOFL, SUBR_NAME )

                                                           ! Get the loads for this subcase from I_PL, J_PL, PL and put into PL_COL
         CALL LINK_MESSAGE_I('GET COL OF PL LOADS FOR                        Subcase', ISUB)
         DO J=1,NDOFL
            PL_COL(J)  = ZERO
            DUM_COL(J) = ZERO
         ENDDO
         CALL GET_SPARSE_CRS_COL ( 'PL        ', ISUB, NTERM_PL, NDOFL, NSUB, I_PL, J_PL, PL, BETA, PL_COL, NULL_COL )
         DO J=1,NDOFL
            DUM_COL(J) = PL_COL(J)
         ENDDO

         IF (DEBUG(32) == 1) THEN                          ! DEBUG output of load vector for this subcase, if requested
            WRITE(F06,3020) ISUB
            CALL WRITE_VECTOR ( '      L-SET LOADS      ',' LOAD', NDOFL, PL_COL )
            WRITE(F06,*)
         ENDIF

                                                           ! Call FBS to solve for displacements for this subcase
         CALL LINK_MESSAGE_I('FBS - SOLVE FOR RHS ANSWERS FOR                   "', ISUB)
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages

         IF      (SOLLIB == 'BANDED  ') THEN

            CALL FBS_LAPACK ( EQUED, NDOFL, KLL_SDIA, EQUIL_SCALE_FACS, DUM_COL )

         ELSE IF (SOLLIB == 'SPARSE  ') THEN

            IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

               SLU_INFO = 0
               CALL FBS_SUPRLU ( SUBR_NAME, 'KLL', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, ISUB, DUM_COL, SLU_INFO )

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

         DO J=1,NDOFL
            UL_COL(J) = DUM_COL(J)
         ENDDO

         IF (DEBUG(33) == 1) THEN                          ! DEBUG output of displs
            WRITE(F06,3022) ISUB
            CALL WRITE_VECTOR ( '      A-SET DISPL      ','DISPL', NDOFL, UL_COL )
            WRITE(F06,*)
         ENDIF

         IF (EPSERR == 'Y') THEN                           ! Calculate residual vector, R. Use RES to calculate EPSILON
            CALL LINK_MESSAGE_I('CALC  EPSILON ERROR ESTIMATE                      "', ISUB)
            CALL EPSCALC ( ISUB )
         ENDIF
                                                           ! Calculate the LAPACK error bounds
         IF ((RCONDK == 'Y') .AND. (SOLLIB == 'BANDED')) THEN
            IF (DABS(RCOND) > MACH_SFMIN) THEN
               CALL LINK_MESSAGE_I('CALC LAPACK ERROR ESTIMATE                        "', ISUB)
               CALL VECINORM ( UL_COL, NDOFL,  UL_INORM )
               CALL VECINORM ( PL_COL, NDOFL,  PL_INORM )
               CALL VECINORM ( RES   , NDOFL, RES_INORM )
               DEN = K_INORM*UL_INORM + PL_INORM
               IF (DABS(DEN) > EPS1) THEN
                  OMEGAI = (RES_INORM)/(DEN)
                  OMEGAI0 = TEN*NDOFL*MACH_EPS
                  LAP_ERR1 = TWO*OMEGAI/RCOND
                  WRITE(F06,3024) ISUB, LAP_ERR1, OMEGAI, RCOND, DEN, RES_INORM, K_INORM, UL_INORM, PL_INORM, OMEGAI0, MACH_EPS
               ELSE
                  WRITE(F06,3026)
               ENDIF
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,3025) ISUB, RCOND, MACH_SFMIN
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,3025) ISUB, RCOND, MACH_SFMIN
               ENDIF
            ENDIF
         ENDIF

         DO J=1,NDOFL                                      ! Write UL to file L3A for this subcase
            WRITE(L3A) UL_COL(J)
         ENDDO

         CALL DEALLOCATE_COL_VEC  ( 'UL_COL' )
         CALL DEALLOCATE_COL_VEC  ( 'PL_COL' )


      ENDDO Solve

FreeS:IF (SOLLIB == 'SPARSE  ') THEN                       ! Last, free the storage allocated inside SuperLU

         IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

            DO J=1,NDOFL                                         ! Need a null col of loads when SuperLU is called to factor KLL
               DUM_COL(J) = ZERO                                  ! (only because it appears in the calling list)
            ENDDO

            CALL C_FORTRAN_DGSSV( 3, NDOFL, NTERM_KLL, 1, KLL , I_KLL , J_KLL , DUM_COL, NDOFL, SLU_FACTORS, SLU_INFO )

            IF (SLU_INFO .EQ. 0) THEN
               WRITE (*,*) 'SUPERLU STORAGE FREED'
            ELSE
               WRITE(*,*) 'SUPERLU STORAGE NOT FREED. INFO FROM SUPERLU FREE STORAGE ROUTINE = ', SLU_INFO
            ENDIF

         ENDIF

      ENDIF FreeS

! Dellocate arrays

      CALL LINK_MESSAGE('DEALLOCATE ARRAYS')
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         CONTINUE
      ELSE
         IF (SOL_NAME(1:12) /= 'GEN CB MODEL' ) THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLL  ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KLL' )
         ENDIF
      ENDIF

      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate ABAND ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate RES   ', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'RES' )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate UL_COL', CR13   ;   CALL DEALLOCATE_COL_VEC  ( 'UL_COL' )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PL_COL', CR13   ;   CALL DEALLOCATE_COL_VEC  ( 'PL_COL' )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PL    ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PL' )

      CALL FILE_CLOSE ( L3A, LINK3A, 'KEEP' )

! Process is now complete so set COMM(LINKNO)

      COMM(LINKNO) = 'C'

! Write data to L1A

      CALL WRITE_L1A ( 'KEEP', 'Y' )

! Check allocation status of allocatable arrays, if requested

      IF (DEBUG(100) > 0) THEN
         CALL CHK_ARRAY_ALLOC_STAT
         IF (DEBUG(100) > 1) THEN
            CALL WRITE_ALLOC_MEM_TABLE ( 'at the end of '//SUBR_NAME )
         ENDIF
      ENDIF

! Write LINK3 end to F06

      CALL OURTIM
      WRITE(F06,151) LINKNO

! Close files

      IF (( DEBUG(193) == 3) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK3' )
      ENDIF

! Write LINK3 end to screen

      WRITE(SC1,153) LINKNO
!***********************************************************************************************************************************
  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END',/)

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  153 FORMAT(  ' >> LINK',I3,' END')

  933 FORMAT(' *ERROR   933: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CRS_CCS  MUST BE EITHER "CRS" OR "CCS" BUT VALUE IS ',A)

  999 FORMAT(' *ERROR   999: INCORRECT SOLUTION IN EXEC CONTROL. SHOULD BE ',A,', BUT IS SOL = ',A)

 3020 FORMAT(//,18X,'LOAD VECTOR FOR SUBCASE ',I8)

 3022 FORMAT(//,18X,'DISPLACEMENTS FOR SUBCASE ',I8,/23X,'LSET DOF',10X,'DISP',14X,'S(J)')

 3024 FORMAT(' *INFORMATION: FOR INTERNAL SUBCASE NUMBER ',I8,' LAPACK ERROR EST (2*OMEGAI/RCOND) = ',1ES13.6,                     &
             ' Gen, slightly > than true err'                                                                                   ,/,&
                                          52X,'................................................................................',/,&
             '                                                    ... OMEGAI                        = ',1ES13.6,                   &
             ' (RES_INORM/DEN)              .'                                                                                  ,/,&
             '                                                    ... RCOND                         = ',1ES13.6,                   &
             ' (Recriprocal of KLL cond num).'                                                                                  ,/,&
             '                                                    ... DEN                           = ',1ES13.6,                   &
             ' (K_INORM*UL_INORM + PL_INORM).'                                                                                  ,/,&
             '                                                    ... RES_INORM                     = ',1ES13.6,                   &
             ' (Inf norm of KLL*UL - PL)    .'                                                                                  ,/,&
             '                                                    ... K_INORM                       = ',1ES13.6,                   &
             ' (Infinity norm of KLL)       .'                                                                                  ,/,&
             '                                                    ... UL_INORM                      = ',1ES13.6,                   &
             ' (Infinity norm of UL displs) .'                                                                                  ,/,&
             '                                                    ... PL_INORM                      = ',1ES13.6,                   &
             ' (Infinity norm of PL loads)  .'                                                                                  ,/,&
             '                                                    ... OMEGAI0 (OMEGAI upper bound)  = ',1ES13.6,                   &
             ' (10*NDOFL*MACH_EPS)          .'                                                                                  ,/,&
             '                                                    ... MACH_EPS                      = ',1ES13.6,                   &
             ' (Machine precision)          .'                                                                                  ,/,&
                                          52X,'................................................................................',/)

 3025 FORMAT(' *WARNING    : CANNOT CALCULATE LAPACK ERROR ESTIMATE FOR INTERNAL SUBCASE NUMBER ',I8                               &
                    ,/,14X,' THE RECIPROCAL OF THE CONDITION NUMBER OF KLL, RCOND         = ',1ES15.6,' CANNOT BE INVERTED.'       &
                    ,/,14X,' IT IS TOO SMALL COMPARED TO MACHINE SAFE MINIMUN (MACH_SFMIN) = ',1ES15.6,/)

 3026 FORMAT(' *INFORMATION: CANNOT CALCULATE OMEGAI. DEN = 0',/)

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

 9998 FORMAT(' *ERROR  9998: COMM ',I3,' INDICATES UNSUCCESSFUL LINK ',I2,' COMPLETION.'                                           &
                    ,/,14X,' FATAL ERROR - CANNOT START LINK ',I2)

12345 FORMAT(A,10X,A)

!***********************************************************************************************************************************

      END SUBROUTINE LINK3





      SUBROUTINE EPSCALC ( ISUB )

   ! Calculates, and prints, EPSILON, an indicator of numerical accuracy in the soln for UA:

   !   EPSILON = UL(t)*[ PL - KLL*UL ]/[ UL(t)*PL ],  UL: displ's, PL: loads, KLL: stiff matrix for the L-set

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NDOFL, NTERM_KLl, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE PARAMS, ONLY                :  EPSIL, SUPINFO, SUPWARN
      USE MACHINE_PARAMS, ONLY        :  MACH_SFMIN
      USE LAPACK_DPB_MATRICES, ONLY   :  RES
      USE SPARSE_MATRICES, ONLY       :  I_KLL, J_KLL, KLL
      USE SPARSE_MATRICES, ONLY       :  SYM_KLL
      USE COL_VECS, ONLY              :  UL_COL, PL_COL

      USE SPARSE_FULL_MULTIPLICATION, ONLY:  MATMULT_SFF
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATADD_FFF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EPSCALC'

      INTEGER(LONG), INTENT(IN)       :: ISUB              ! Internal subcase no. (1 to NSUB)


      REAL(DOUBLE) , PARAMETER        :: ALPHA     =  ONE  ! Scalar multiplier for KLL in calc'ing residual vector, RES
      REAL(DOUBLE) , PARAMETER        :: BETA      = -ONE  ! Scalar multiplier for PL in calc'ing residual vector, RES
      REAL(DOUBLE)                    :: DEN               ! Denominator in EPSILON calculation
      REAL(DOUBLE)                    :: EPSILON           ! The indicator of numerical accuracy of the displ soln, UL
      REAL(DOUBLE)                    :: KU(NDOFL)         ! Result of multiplying KLL and UL_COL
      REAL(DOUBLE)                    :: NUM               ! Numerator in EPSILON calculation

      REAL(DOUBLE), EXTERNAL          :: DDOT              ! BLAS dot-product function

      INTRINSIC                       :: DABS



   ! **********************************************************************************************************************************
   ! Calculate residual vector. First, multiply KLL x UL_COL and then add -PL_COL:

      CALL MATMULT_SFF ( 'KLL', NDOFL, NDOFL, NTERM_KLL, SYM_KLL, I_KLL, J_KLL, KLL, 'UL_COL', NDOFL, 1, UL_COL, 'Y',             &
                'KLL*UL_COL',  ONE, KU )
      CALL MATADD_FFF  ( KU, PL_COL, NDOFL, 1, ALPHA, BETA, 0, RES )


   ! Calculate dot product of UL displ vector and RES residual vector

      NUM = DDOT (NDOFL, UL_COL, 1, RES, 1)

   ! Calculate dot product of UL displ vector and PL load vector

      DEN = DDOT (NDOFL, UL_COL, 1, PL_COL, 1)

   ! Calculate EPSILON and print

      IF (DABS(DEN) > MACH_SFMIN) THEN
         EPSILON = NUM/DEN
         WRITE(F06,3701) ISUB,EPSILON
      ELSE
         WARN_ERR = WARN_ERR + 1
         WRITE(F06,3702) ISUB, DEN, MACH_SFMIN
      ENDIF



      RETURN

   ! **********************************************************************************************************************************
    3701 FORMAT(' *INFORMATION: FOR INTERNAL SUBCASE NUMBER ',I8,' EPSILON ERROR ESTIMATE            = ',1ES13.6,                     &
                  ' Based on U''*(K*U - P)/(U''*P)',/)

    3702 FORMAT(' *WARNING    : CANNOT CALCULATE EPSILON ERROR ESTIMATE FOR INTERNAL SUBCASE NUMBER ',I8                              &
              ,/,14X,' THE DOT PRODUCT OF DISPL AND LOAD VECTORS                     = ',1ES15.6,' CANNOT BE INVERTED.'      &
              ,/,14X,' IT IS TOO SMALL COMPARED TO MACHINE SAFE MINIMUM (MACH_SFMIN) = ',1ES15.6,/)

   ! **********************************************************************************************************************************

      END SUBROUTINE EPSCALC





      SUBROUTINE VECINORM ( X, N, X_INORM )

   ! Calculates the infinity norm, X_INORM, of an input vector, X, of dimension N. The infinity norm of a vector is the absolute value
   ! of the numerically largest term in the vector.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'VECINORM'

      INTEGER(LONG), INTENT(IN)       :: N                 ! Dimension of the input vector X
      INTEGER(LONG)                   :: I                 ! DO loop index


      REAL(DOUBLE),  INTENT(IN)       :: X(N)              ! The input vector for which the infinity norm is calc'd
      REAL(DOUBLE),  INTENT(OUT)      :: X_INORM           ! The calc'd infinity norm of X

      INTRINSIC                       :: DABS



   ! **********************************************************************************************************************************
   ! Calculate infinity norm of a vector

      X_INORM = ZERO
      DO I=1,N
         IF (DABS(X(I)) > X_INORM) THEN
         X_INORM = DABS(X(I))
         ENDIF
      ENDDO



      RETURN

   ! **********************************************************************************************************************************

      END SUBROUTINE VECINORM





      END MODULE LINK3_MOD
