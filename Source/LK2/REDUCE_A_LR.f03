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

   MODULE REDUCTION_A_TO_L

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: REDUCE_A_LR

   CONTAINS

      SUBROUTINE REDUCE_A_LR

! Call routines to reduce stiffness, mass, loads from A-set to L, R-sets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, LINKNO,   NDOFA, NDOFG, NDOFL, NDOFR, NSUB, SOL_NAME,                       &
                                         NTERM_KAA , NTERM_KLL , NTERM_KRL , NTERM_KRR ,                                           &
                                         NTERM_KAAD, NTERM_KLLD, NTERM_KRLD, NTERM_KRRD,                                           &
                                         NTERM_MAA , NTERM_MLL , NTERM_MRL , NTERM_MRR ,                                           &
                                         NTERM_PA  , NTERM_PL  , NTERM_PR
      USE TIMDAT, ONLY                :  TSEC, YEAR, MONTH, DAY, HOUR, MINUTE, SEC, SFRAC, STIME
      USE CONSTANTS_1, ONLY           :  ONE
      USE DOF_TABLES, ONLY            :  TDOFI
      USE RIGID_BODY_DISP_MATS, ONLY  :  RBGLOBAL_ASET, RBGLOBAL_GSET, RBGLOBAL_LSET
      USE PARAMS, ONLY                :  EQCHK_OUTPUT, MATSPARS, PRTSTIFD, PRTSTIFF, PRTMASS, PRTFOR
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE SPARSE_MATRICES, ONLY       :  I_KAA , J_KAA , KAA , I_KLL , J_KLL , KLL , I_KRL , J_KRL , KRL , I_KRR , J_KRR , KRR ,   &
                                         I_KAAD, J_KAAD, KAAD, I_KLLD, J_KLLD, KLLD, I_KRLD, J_KRLD, KRLD, I_KRRD, J_KRRD, KRRD,   &
                                         I_MAA , J_MAA , MAA , I_MLL , J_MLL , MLL , I_MRL , J_MRL , MRL , I_MRR , J_MRR , MRR ,   &
                                         I_PA  , J_PA  , PA  , I_PL  , J_PL  , PL  , I_PR  , J_PR  , PR
      USE SPARSE_MATRICES, ONLY       :  SYM_KLL
      USE FULL_MATRICES, ONLY         :  KAA_FULL
      USE OUTPUT4_MATRICES, ONLY      :  ACT_OU4_MYSTRAN_NAMES, NUM_OU4_REQUESTS
      USE DEBUG_PARAMETERS

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_VEC
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_MATRIX_DIAG_STATS
      USE RIGID_BODY_STORAGE_LIFECYCLE, ONLY:  ALLOCATE_RBGLOBAL, DEALLOCATE_RBGLOBAL
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE REDUCTION_CHECKS, ONLY      :  STIFF_MAT_EQUIL_CHK

      USE SPARSE_FORMAT_CONVERSION, ONLY:  SPARSE_CRS_TO_FULL
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT
      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_A_LR'
      CHARACTER(  1*BYTE)             :: DEALLOCATE_KAA = 'Y'  ! Indicator of whether we need to keep KAA allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_KRL = 'Y'  ! Indicator of whether we need to keep KRL allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_KRR = 'Y'  ! Indicator of whether we need to keep KRR allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_KAAD= 'Y'  ! Indicator of whether we need to keep KAA allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_KRLD= 'Y'  ! Indicator of whether we need to keep KRL allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_MAA = 'Y'  ! Indicator of whether we need to keep MAA allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_MRL = 'Y'  ! Indicator of whether we need to keep MRL allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_MRR = 'Y'  ! Indicator of whether we need to keep MRR allocated for OU4 output
      CHARACTER(  1*BYTE)             :: DEALLOCATE_PA  = 'Y'  ! Indicator of whether we need to keep PA  allocated for OU4 output
      CHARACTER(132*BYTE)             :: MATRIX_NAME           ! Name of matrix for printout
      CHARACTER( 44*BYTE)             :: MODNAM                ! Name to write to screen to describe module being run

      INTEGER(LONG)                   :: DO_WHICH_CODE_FRAG    ! 1 or 2 depending on which seg of code to run (depends on BUCKLING)
      INTEGER(LONG)                   :: L_SET_COL             ! Col no. in array TDOFI where the F-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: L_SET_DOF             ! F-set DOF number
      INTEGER(LONG)                   :: I,J                   ! DO loop indices
      INTEGER(LONG)                   :: PART_VEC_A_LR(NDOFA)  ! Partitioning vector (N set into F and S sets)
      INTEGER(LONG)                   :: PART_VEC_SUB(NSUB)    ! Partitioning vector (1's for all subcases)


      REAL(DOUBLE)                    :: KLL_DIAG(NDOFL)       ! Diagonal terms from KLL
      REAL(DOUBLE)                    :: KLL_MAX_DIAG          ! Max diag term from  KLL
      REAL(DOUBLE)                    :: KLLD_DIAG(NDOFL)      ! Diagonal terms from KLLD
      REAL(DOUBLE)                    :: KLLD_MAX_DIAG         ! Max diag term from  KLLD



! **********************************************************************************************************************************
! Determine if we need to keep any OUTPUT4 matrices allocated until after they are processed in LINK2

      IF (NUM_OU4_REQUESTS > 0) THEN
         DO I=1,NUM_OU4_REQUESTS
            IF      (ACT_OU4_MYSTRAN_NAMES(I)(1:3) == 'KAA') THEN
               DEALLOCATE_KAA = 'N'
            ELSE IF (ACT_OU4_MYSTRAN_NAMES(I)(1:3) == 'KRL') THEN
               DEALLOCATE_KRL = 'N'
            ELSE IF (ACT_OU4_MYSTRAN_NAMES(I)(1:3) == 'KRR') THEN
               DEALLOCATE_KRR = 'N'
            ELSE IF (ACT_OU4_MYSTRAN_NAMES(I)(1:3) == 'MAA') THEN
               DEALLOCATE_MAA = 'N'
            ELSE IF (ACT_OU4_MYSTRAN_NAMES(I)(1:3) == 'MRL') THEN
               DEALLOCATE_MRL = 'N'
            ELSE IF (ACT_OU4_MYSTRAN_NAMES(I)(1:3) == 'MRR') THEN
               DEALLOCATE_MRR = 'N'
            ELSE IF (ACT_OU4_MYSTRAN_NAMES(I)(1:3) == 'PA' ) THEN
               DEALLOCATE_PA  = 'N'
            ENDIF
         ENDDO
      ENDIF

! **********************************************************************************************************************************
! Depending on whether this is a BUCKLING soln (and LOAD_ISTEP value) or not, one or another segment of code will be run

      IF ((SOL_NAME(1:8) == 'BUCKLING')) THEN
         IF      (LOAD_ISTEP == 1) THEN
            DO_WHICH_CODE_FRAG = 1
         ELSE IF (LOAD_ISTEP == 2) THEN
            DO_WHICH_CODE_FRAG = 2
         ENDIF
      ELSE
         DO_WHICH_CODE_FRAG = 1
      ENDIF

! **********************************************************************************************************************************
      IF (DO_WHICH_CODE_FRAG == 1) THEN                    ! This is for all except BUCKLING w LOAD_ISTEP=2 (eigen part of BUCKLING)

! Do reduction

         IF (NDOFR > 0) THEN                               ! If NDOFR > 0 reduce KAA to KLL

! Reduce KAA to KLL

            CALL PARTITION_VEC ( NDOFA, 'A ', 'L ', 'R ', PART_VEC_A_LR )

            DO I=1,NSUB
               PART_VEC_SUB = 1
            ENDDO

            IF (NTERM_KAA > 0) THEN

               CALL OURTIM
               MODNAM = '  REDUCE KAA TO KLL (PARTITION, ONLY)'
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_KAA_TO_KLL ( PART_VEC_A_LR )

            ELSE

               NTERM_KLL = 0
               NTERM_KRL = 0
               NTERM_KRR = 0
               CALL ALLOCATE_SPARSE_MAT ( 'KLL', NDOFL, NTERM_KLL, SUBR_NAME )

            ENDIF

! Reduce MAA to MLL

            IF (NTERM_MAA > 0) THEN

               CALL OURTIM
               MODNAM = '  REDUCE MAA TO MLL(PARTITION, ONLY)'
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_MAA_TO_MLL ( PART_VEC_A_LR )

            ELSE

               NTERM_MLL = 0
               NTERM_MRL = 0
               NTERM_MRR = 0
               CALL ALLOCATE_SPARSE_MAT ( 'MLL', NDOFL, NTERM_MLL, SUBR_NAME )

            ENDIF

! Reduce PA to PL.

            IF ((SOL_NAME(1:5) /= 'MODES') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')) THEN

               CALL OURTIM
               IF (MATSPARS == 'Y') THEN
                  MODNAM = '  REDUCE PA TO PL (SPARSE MATRIX ROUTINES)'
               ELSE
                  MODNAM = '  REDUCE PA TO PL (FULL MATRIX ROUTINES)'
               ENDIF
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_PA_TO_PL ( PART_VEC_A_LR, PART_VEC_SUB )

            ENDIF

         ELSE

! There is no R-set, so equate A, L sets

            CALL OURTIM
            MODNAM = '  EQUATING L-SET TO A-SET'
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            NDOFL     = NDOFA

            NTERM_KLL = NTERM_KAA
            NTERM_KRL = 0
            NTERM_KRR = 0

            NTERM_MLL = NTERM_MAA
            NTERM_MRL = 0
            NTERM_MRR = 0

            NTERM_PL  = NTERM_PA
            NTERM_PR  = 0

            CALL ALLOCATE_SPARSE_MAT ( 'KLL', NDOFL, NTERM_KLL, SUBR_NAME )

            CALL ALLOCATE_SPARSE_MAT ( 'MLL', NDOFL, NTERM_MLL, SUBR_NAME )

            IF ((SOL_NAME(1:5) /= 'MODES') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')) THEN

               CALL ALLOCATE_SPARSE_MAT ( 'PL', NDOFL, NTERM_PL, SUBR_NAME )

            ENDIF

         ENDIF

! Deallocate A-set arrays

         IF (DEBUG(57) > 0) THEN
            CALL ALLOCATE_FULL_MAT ( 'KAA_FULL', NDOFA, NDOFA, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'KAA', NTERM_KAA, NDOFA, NDOFA, 'Y', I_KAA, J_KAA, KAA, KAA_FULL )
            WRITE(F06,*) ' A-Set Stiffness Matrix, KAA'
            DO I=1,NDOFA
               IF (DEBUG(57) == 2) THEN
                  WRITE(F06,*) '  ROW',I
               ENDIF
               WRITE(F06,3006) (KAA_FULL(I,J),J=1,NDOFA)
               WRITE(F06,*)
            ENDDO
 3006       FORMAT(6(1ES15.6))
         ENDIF

         MODNAM = '  DEALLOCATE A-SET ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages

         IF (DEALLOCATE_KAA == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KAA', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KAA' )
            WRITE(SC1,*) CR13
         ENDIF

         IF (DEALLOCATE_MAA == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MAA', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MAA' )
            WRITE(SC1,*) CR13
         ENDIF

         IF (DEALLOCATE_PA == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PA ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PA' )
            WRITE(SC1,*) CR13
         ENDIF

! Print out stiffness matrix partitions, if requested

         IF (( PRTSTIFF(5) == 1) .OR. ( PRTSTIFF(5) == 3)) THEN
            IF (NTERM_KLL > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KLL'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'L ', 'L ', NTERM_KLL, NDOFL, I_KLL, J_KLL, KLL )
            ENDIF
         ENDIF

         IF (( PRTSTIFF(5) == 2) .OR. ( PRTSTIFF(5) == 3)) THEN
            IF (NTERM_KRL > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KRL'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'R ', 'L ', NTERM_KRL, NDOFR, I_KRL, J_KRL, KRL )
            ENDIF
            IF (NTERM_KRR > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KRR'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'R ', 'R ', NTERM_KRR, NDOFR, I_KRR, J_KRR, KRR )
            ENDIF
         ENDIF

         MODNAM = '  DEALLOCATE R-SET ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages

         IF (DEALLOCATE_KRL == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KRL', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KRL' )
            WRITE(SC1,*) CR13
         ENDIF

         IF (DEALLOCATE_KRR == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KRR', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KRR' )
            WRITE(SC1,*) CR13
         ENDIF

! Write matrix diagonal and stats, if requested.
! NOTE: call this subr even if PRTSTIFFD(2) = 0 since we need KNN_DIAG, KNN_MAX_DIAG for the equilibrium check

         CALL GET_MATRIX_DIAG_STATS ( 'KLL', 'L ', NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, PRTSTIFD(5), KLL_DIAG, KLL_MAX_DIAG )

! Print out mass matrix partitions, if requested

         IF (( PRTMASS(5) == 1) .OR. ( PRTMASS(5) == 3)) THEN
            IF (NTERM_MLL > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MLL'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'L ', 'L ', NTERM_MLL, NDOFL, I_MLL, J_MLL, MLL )
            ENDIF
         ENDIF

         IF (( PRTMASS(5) == 2) .OR. ( PRTMASS(5) == 3)) THEN
            IF (NTERM_MRL > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MRL'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'R ', 'L ', NTERM_MRL, NDOFR, I_MRL, J_MRL, MRL )
            ENDIF
            IF (NTERM_MRR > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MRR'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'R ', 'R ', NTERM_MRR, NDOFR, I_MRR, J_MRR, MRR )
            ENDIF
         ENDIF

         IF (DEALLOCATE_MRL == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MRL', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MRL' )
            WRITE(SC1,*) CR13
         ENDIF

         IF (DEALLOCATE_MRR == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MRR', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MRR' )
            WRITE(SC1,*) CR13
         ENDIF

! Print out load matrix partitions, if requested

         IF (( PRTFOR(5) == 1) .OR. ( PRTFOR(5) == 3)) THEN
            IF (NTERM_PL  > 0) THEN
               MATRIX_NAME = 'LOAD MATRIX PL'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'L ', 'SUBCASE', NTERM_PL, NDOFL, I_PL, J_PL, PL )
            ENDIF
         ENDIF

         IF (( PRTFOR(5) == 2) .OR. ( PRTFOR(5) == 3)) THEN
            IF (NTERM_PR  > 0) THEN
               MATRIX_NAME = 'LOAD MATRIX PR'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'R ', 'SUBCASE', NTERM_PR, NDOFR, I_PR, J_PR, PR )
            ENDIF
         ENDIF

         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PR ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PR' )
         WRITE(SC1,*) CR13

! Do equilibrium check on the L-set stiffness matrix, if requested

         IF (EQCHK_OUTPUT(5) > 0) THEN
            CALL ALLOCATE_RBGLOBAL ( 'L ', SUBR_NAME )
            IF (NDOFR > 0) THEN
               CALL TDOF_COL_NUM ( 'L ', L_SET_COL )
               DO I=1,NDOFG
                  L_SET_DOF = TDOFI(I,L_SET_COL)
                  IF (L_SET_DOF > 0) THEN
                     DO J=1,6
                        RBGLOBAL_LSET(L_SET_DOF,J) = RBGLOBAL_GSET(I,J)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               DO I=1,NDOFL
                  DO J=1,6
                     RBGLOBAL_LSET(I,J) = RBGLOBAL_ASET(I,J)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF

         IF (EQCHK_OUTPUT(5) > 0) THEN
            CALL OURTIM
            MODNAM = '  EQUILIBRIUM CHECK ON KLL                '
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
            CALL STIFF_MAT_EQUIL_CHK ( EQCHK_OUTPUT(3),'L ', SYM_KLL, NDOFL, NTERM_KLL, I_KLL, J_KLL, KLL, KLL_DIAG, KLL_MAX_DIAG, &
                                       RBGLOBAL_LSET)
         ENDIF
         CALL DEALLOCATE_RBGLOBAL ( 'A ' )
         CALL DEALLOCATE_RBGLOBAL ( 'L ' )



! **********************************************************************************************************************************
      ELSE                                                 ! This is BUCKLING with LOAD_ISTEP = 2 (eigen part of BUCKLING)

! Do reduction

         IF (NDOFR > 0) THEN                               ! If NDOFR > 0 reduce KAA to KLL

! Reduce KAAD to KLLD

            CALL PARTITION_VEC ( NDOFA, 'A ', 'L ', 'R ', PART_VEC_A_LR )

            DO I=1,NSUB
               PART_VEC_SUB = 1
            ENDDO

            IF (NTERM_KAAD > 0) THEN

               CALL OURTIM
               MODNAM = '  REDUCE KAAD TO KLLD (PARTITION, ONLY)'
               WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

               CALL REDUCE_KAAD_TO_KLLD ( PART_VEC_A_LR )

            ELSE

               NTERM_KLLD = 0
               NTERM_KRLD = 0
               NTERM_KRRD = 0
               CALL ALLOCATE_SPARSE_MAT ( 'KLLD', NDOFL, NTERM_KLLD, SUBR_NAME )

            ENDIF

         ELSE

! There is no R-set, so equate A, L sets

            CALL OURTIM
            MODNAM = '  EQUATING L-SET TO A-SET'
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            NDOFL     = NDOFA

            NTERM_KLLD = NTERM_KAAD
            NTERM_KRLD = 0
            NTERM_KRRD = 0

            CALL ALLOCATE_SPARSE_MAT ( 'KLLD', NDOFL, NTERM_KLLD, SUBR_NAME )

         ENDIF

! Deallocate A-set arrays

         MODNAM = '  DEALLOCATE A-SET ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages

         IF (DEALLOCATE_KAAD == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KAAD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KAAD' )
            WRITE(SC1,*) CR13
         ENDIF

! Print out stiffness matrix partitions, if requested

         IF (( PRTSTIFF(5) == 1) .OR. ( PRTSTIFF(5) == 3)) THEN
            IF (NTERM_KLLD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KLLD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'L ', 'L ', NTERM_KLLD, NDOFL, I_KLLD, J_KLLD, KLLD )
            ENDIF
         ENDIF

         IF (( PRTSTIFF(5) == 2) .OR. ( PRTSTIFF(5) == 3)) THEN
            IF (NTERM_KRLD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KRLD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'R ', 'L ', NTERM_KRLD, NDOFR, I_KRLD, J_KRLD, KRLD )
            ENDIF
            IF (NTERM_KRRD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KRRD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'R ', 'R ', NTERM_KRRD, NDOFR, I_KRRD, J_KRRD, KRRD )
            ENDIF
         ENDIF

         MODNAM = '  DEALLOCATE R-SET ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages

         IF (DEALLOCATE_KRLD == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KRLD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KRLD' )
            WRITE(SC1,*) CR13
         ENDIF

         IF (DEALLOCATE_KRR == 'Y') THEN
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KRRD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KRRD' )
            WRITE(SC1,*) CR13
         ENDIF

! Write matrix diagonal and stats, if requested.
! NOTE: call this subr even if PRTSTIFFD(2) = 0 since we need KNN_DIAG, KNN_MAX_DIAG for the equilibrium check

         CALL GET_MATRIX_DIAG_STATS ( 'KLLD', 'L ', NDOFL, NTERM_KLLD, I_KLLD, J_KLLD, KLLD, PRTSTIFD(5), KLLD_DIAG, KLLD_MAX_DIAG )

      ENDIF



      RETURN

! **********************************************************************************************************************************
 2092 FORMAT(4X,A44,20X,I2,':',I2,':',I2,'.',I3)

 9995 FORMAT(/,' PROCESSING ENDED IN LINK ',I3,' DUE TO ABOVE ',I8,' ERRORS')



12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_A_LR


      SUBROUTINE REDUCE_KAAD_TO_KLLD ( PART_VEC_A_LR )

! Call routines to reduce the KAAD differential stiffness matrix from the A-set to the L, R-sets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L2K, L2L, LINK2K, LINK2L, L2K_MSG, L2L_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFA, NDOFL, NDOFR, NTERM_KAAD, NTERM_KLLD, NTERM_KRLD,         &
                                         NTERM_KRRD,  SOL_NAME
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE SPARSE_MATRICES, ONLY       :  I_KAAD, J_KAAD, KAAD, I_KLLD, J_KLLD, KLLD, I_KRLD, J_KRLD, KRLD, I_KRRD, J_KRRD, KRRD,   &
                                         SYM_KAAD, SYM_KLLD, SYM_KRLD, SYM_KRRD
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_KAAD_TO_KLLD'

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_A_LR(NDOFA)! Partitioning vector (F set into A and O sets)
      INTEGER(LONG)                   :: KLLD_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KRLD_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KRRD_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices




! **********************************************************************************************************************************
! Partition KLLD from KAAD (This is KLLD before reduction, or KLLD(bar) )

      IF (NDOFL > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KAAD', NTERM_KAAD, NDOFA, NDOFA, SYM_KAAD, I_KAAD, J_KAAD,      PART_VEC_A_LR, PART_VEC_A_LR,  &
                                    NUM1, NUM1, KLLD_ROW_MAX_TERMS, 'KLLD', NTERM_KLLD, SYM_KLLD )

         CALL ALLOCATE_SPARSE_MAT ( 'KLLD', NDOFL, NTERM_KLLD, SUBR_NAME )

         IF (NTERM_KLLD > 0) THEN
            CALL PARTITION_SS ( 'KAAD', NTERM_KAAD, NDOFA, NDOFA, SYM_KAAD, I_KAAD, J_KAAD, KAAD, PART_VEC_A_LR, PART_VEC_A_LR,    &
                                 NUM1, NUM1, KLLD_ROW_MAX_TERMS, 'KLLD', NTERM_KLLD, NDOFL, SYM_KLLD, I_KLLD, J_KLLD, KLLD )
         ENDIF

      ENDIF

! Partition KRLD from KAAD

      IF ((NDOFL > 0) .AND. (NDOFR > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'KAAD', NTERM_KAAD, NDOFA, NDOFA, SYM_KAAD, I_KAAD, J_KAAD,      PART_VEC_A_LR, PART_VEC_A_LR,  &
                                    NUM2, NUM1, KRLD_ROW_MAX_TERMS, 'KRLD', NTERM_KRLD, SYM_KRLD )

         CALL ALLOCATE_SPARSE_MAT ( 'KRLD', NDOFR, NTERM_KRLD, SUBR_NAME )

         IF (NTERM_KRLD > 0) THEN
            CALL PARTITION_SS ( 'KAAD', NTERM_KAAD, NDOFA, NDOFA, SYM_KAAD, I_KAAD, J_KAAD, KAAD, PART_VEC_A_LR, PART_VEC_A_LR,    &
                                 NUM2, NUM1, KRLD_ROW_MAX_TERMS, 'KRLD', NTERM_KRLD, NDOFR, SYM_KRLD, I_KRLD, J_KRLD, KRLD )
         ENDIF

      ENDIF

! Partition KRRD from KAAD

      IF (NDOFR > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KAAD', NTERM_KAAD, NDOFA, NDOFA, SYM_KAAD, I_KAAD, J_KAAD,      PART_VEC_A_LR, PART_VEC_A_LR,  &
                                    NUM2, NUM2, KRRD_ROW_MAX_TERMS, 'KRRD', NTERM_KRRD, SYM_KRRD )

         CALL ALLOCATE_SPARSE_MAT ( 'KRRD', NDOFR, NTERM_KRRD, SUBR_NAME )

         IF (NTERM_KRRD > 0) THEN
            CALL PARTITION_SS ( 'KAAD', NTERM_KAAD, NDOFA, NDOFA, SYM_KAAD, I_KAAD, J_KAAD, KAAD, PART_VEC_A_LR, PART_VEC_A_LR,    &
                                 NUM2, NUM2, KRRD_ROW_MAX_TERMS, 'KRRD', NTERM_KRRD, NDOFR, SYM_KRRD, I_KRRD, J_KRRD, KRRD )
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_KAAD_TO_KLLD


      SUBROUTINE REDUCE_KAA_TO_KLL ( PART_VEC_A_LR )

! Call routines to reduce the KAA linear stiffness matrix from the A-set to the L, R-sets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L2K, L2L, LINK2K, LINK2L, L2K_MSG, L2L_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFA, NDOFL, NDOFR, NTERM_KAA, NTERM_KLL, NTERM_KRL, NTERM_KRR, &
                                         SOL_NAME
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE SPARSE_MATRICES, ONLY       :  I_KAA, J_KAA, KAA, I_KLL, J_KLL, KLL, I_KRL, J_KRL, KRL, I_KRR, J_KRR, KRR,               &
                                         SYM_KAA, SYM_KLL, SYM_KRL, SYM_KRR
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_MATRIX_1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_KAA_TO_KLL'

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_A_LR(NDOFA)! Partitioning vector (F set into A and O sets)
      INTEGER(LONG)                   :: KLL_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KRL_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KRR_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices




! **********************************************************************************************************************************
! Partition KLL from KAA (This is KLL before reduction, or KLL(bar) )

      IF (NDOFL > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KAA', NTERM_KAA, NDOFA, NDOFA, SYM_KAA, I_KAA, J_KAA,      PART_VEC_A_LR, PART_VEC_A_LR,       &
                                    NUM1, NUM1, KLL_ROW_MAX_TERMS, 'KLL', NTERM_KLL, SYM_KLL )

         CALL ALLOCATE_SPARSE_MAT ( 'KLL', NDOFL, NTERM_KLL, SUBR_NAME )

         IF (NTERM_KLL > 0) THEN
            CALL PARTITION_SS ( 'KAA', NTERM_KAA, NDOFA, NDOFA, SYM_KAA, I_KAA, J_KAA, KAA, PART_VEC_A_LR, PART_VEC_A_LR,          &
                                 NUM1, NUM1, KLL_ROW_MAX_TERMS, 'KLL', NTERM_KLL, NDOFL, SYM_KLL, I_KLL, J_KLL, KLL )
         ENDIF

      ENDIF

! Partition KRL from KAA

      IF ((NDOFL > 0) .AND. (NDOFR > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'KAA', NTERM_KAA, NDOFA, NDOFA, SYM_KAA, I_KAA, J_KAA,      PART_VEC_A_LR, PART_VEC_A_LR,       &
                                    NUM2, NUM1, KRL_ROW_MAX_TERMS, 'KRL', NTERM_KRL, SYM_KRL )

         CALL ALLOCATE_SPARSE_MAT ( 'KRL', NDOFR, NTERM_KRL, SUBR_NAME )

         IF (NTERM_KRL > 0) THEN
            CALL PARTITION_SS ( 'KAA', NTERM_KAA, NDOFA, NDOFA, SYM_KAA, I_KAA, J_KAA, KAA, PART_VEC_A_LR, PART_VEC_A_LR,          &
                                 NUM2, NUM1, KRL_ROW_MAX_TERMS, 'KRL', NTERM_KRL, NDOFR, SYM_KRL, I_KRL, J_KRL, KRL )
         ENDIF

      ENDIF

! Partition KRR from KAA

      IF (NDOFR > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KAA', NTERM_KAA, NDOFA, NDOFA, SYM_KAA, I_KAA, J_KAA,      PART_VEC_A_LR, PART_VEC_A_LR,       &
                                    NUM2, NUM2, KRR_ROW_MAX_TERMS, 'KRR', NTERM_KRR, SYM_KRR )

         CALL ALLOCATE_SPARSE_MAT ( 'KRR', NDOFR, NTERM_KRR, SUBR_NAME )

         IF (NTERM_KRR > 0) THEN
            CALL PARTITION_SS ( 'KAA', NTERM_KAA, NDOFA, NDOFA, SYM_KAA, I_KAA, J_KAA, KAA, PART_VEC_A_LR, PART_VEC_A_LR,          &
                                 NUM2, NUM2, KRR_ROW_MAX_TERMS, 'KRR', NTERM_KRR, NDOFR, SYM_KRR, I_KRR, J_KRR, KRR )
         ENDIF

      ENDIF

! Write matrices needed for Craig-Bampton, if this is a CB soln

      IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
         CALL WRITE_MATRIX_1 ( LINK2K, L2K, 'Y', 'KEEP', L2K_MSG, 'KRL', NTERM_KRL, NDOFR, I_KRL, J_KRL, KRL )
         CALL WRITE_MATRIX_1 ( LINK2L, L2L, 'Y', 'KEEP', L2L_MSG, 'KRR', NTERM_KRR, NDOFR, I_KRR, J_KRR, KRR )
      ENDIF



      RETURN

! **********************************************************************************************************************************
 2092 FORMAT(4X,A44,20X,I2,':',I2,':',I2,'.',I3)

 2400 FORMAT(' *ERROR  2400: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THERE IS AN O-SET BUT GUYAN REDUCTION MATRIX GOA HAS ',I12,' TERMS IN IT. MUST BE > 0')

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_KAA_TO_KLL


      SUBROUTINE REDUCE_MAA_TO_MLL ( PART_VEC_A_LR )

! Call routines to reduce the MAA mass matrix from the A-set to the L, R-sets. See Appendix B to the MYSTRAN User's Reference Manual
! Reference Manual for the derivation of the reduction equations.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, L2M, L2N, LINK2M, LINK2N, L2M_MSG, L2N_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFA, NDOFL, NDOFR, NTERM_MAA, NTERM_MLL, NTERM_MRL, NTERM_MRR, &
                                         SOL_NAME
      USE PARAMS, ONLY                :  EPSIL
      USE TIMDAT, ONLY                :  TSEC
      USE SPARSE_MATRICES, ONLY       :  I_MAA, J_MAA, MAA, I_MLL, J_MLL, MLL, I_MRL, J_MRL, MRL, I_MRR, J_MRR, MRR,               &
                                         SYM_MAA, SYM_MLL, SYM_MRL, SYM_MRR
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_MATRIX_1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_MAA_TO_MLL'

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_A_LR(NDOFA)! Partitioning vector (F set into A and O sets)
      INTEGER(LONG)                   :: MLL_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: MRL_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: MRR_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices




! **********************************************************************************************************************************
! Partition MLL from MAA

      IF (NDOFL > 0) THEN

         CALL PARTITION_SS_NTERM ( 'MAA', NTERM_MAA, NDOFA, NDOFA, SYM_MAA, I_MAA, J_MAA,      PART_VEC_A_LR, PART_VEC_A_LR,       &
                                    NUM1, NUM1, MLL_ROW_MAX_TERMS, 'MLL', NTERM_MLL, SYM_MLL )

         CALL ALLOCATE_SPARSE_MAT ( 'MLL', NDOFL, NTERM_MLL, SUBR_NAME )

         IF (NTERM_MLL > 0) THEN
            CALL PARTITION_SS ( 'MAA', NTERM_MAA, NDOFA, NDOFA, SYM_MAA, I_MAA, J_MAA, MAA, PART_VEC_A_LR, PART_VEC_A_LR,          &
                                 NUM1, NUM1, MLL_ROW_MAX_TERMS, 'MLL', NTERM_MLL, NDOFL, SYM_MLL, I_MLL, J_MLL, MLL )
         ENDIF

      ENDIF

! Partition MRL from MAA

      IF ((NDOFL > 0) .AND. (NDOFR > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'MAA', NTERM_MAA, NDOFA, NDOFA, SYM_MAA, I_MAA, J_MAA,      PART_VEC_A_LR, PART_VEC_A_LR,       &
                                    NUM2, NUM1, MRL_ROW_MAX_TERMS, 'MRL', NTERM_MRL, SYM_MRL )

         CALL ALLOCATE_SPARSE_MAT ( 'MRL', NDOFR, NTERM_MRL, SUBR_NAME )

         IF (NTERM_MRL > 0) THEN
            CALL PARTITION_SS ( 'MAA', NTERM_MAA, NDOFA, NDOFA, SYM_MAA, I_MAA, J_MAA, MAA, PART_VEC_A_LR, PART_VEC_A_LR,          &
                                 NUM2, NUM1, MRL_ROW_MAX_TERMS, 'MRL', NTERM_MRL, NDOFR, SYM_MRL, I_MRL, J_MRL, MRL )
         ENDIF

      ENDIF

! Partition MRR from MAA

      IF (NDOFR > 0) THEN

         CALL PARTITION_SS_NTERM ( 'MAA', NTERM_MAA, NDOFA, NDOFA, SYM_MAA, I_MAA, J_MAA,      PART_VEC_A_LR, PART_VEC_A_LR,       &
                                    NUM2, NUM2, MRR_ROW_MAX_TERMS, 'MRR', NTERM_MRR, SYM_MRR )

         CALL ALLOCATE_SPARSE_MAT ( 'MRR', NDOFR, NTERM_MRR, SUBR_NAME )

         IF (NTERM_MRR > 0) THEN
            CALL PARTITION_SS ( 'MAA', NTERM_MAA, NDOFA, NDOFA, SYM_MAA, I_MAA, J_MAA, MAA, PART_VEC_A_LR, PART_VEC_A_LR,          &
                                 NUM2, NUM2, MRR_ROW_MAX_TERMS, 'MRR', NTERM_MRR, NDOFR, SYM_MRR, I_MRR, J_MRR, MRR )
         ENDIF

      ENDIF

! Write matrices needed for Craig-Bampton, if this is a CB soln

      IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
         CALL WRITE_MATRIX_1 ( LINK2M, L2M, 'Y', 'KEEP', L2M_MSG, 'MRL', NTERM_MRL, NDOFR, I_MRL, J_MRL, MRL )
         CALL WRITE_MATRIX_1 ( LINK2N, L2N, 'Y', 'KEEP', L2N_MSG, 'MRR', NTERM_MRR, NDOFR, I_MRR, J_MRR, MRR )
      ENDIF



      RETURN

! **********************************************************************************************************************************

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_MAA_TO_MLL


      SUBROUTINE REDUCE_PA_TO_PL ( PART_VEC_A_LR, PART_VEC_SUB )

! Call routines to reduce the PA grid point load matrix from the A-set to the L, R-sets. See Appendix B to the MYSTRAN User's
! Reference Manual for the derivation of the reduction equations.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFA, NDOFL, NDOFR, NSUB, NTERM_GOA, NTERM_PA, NTERM_PL, NTERM_PR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE SPARSE_MATRICES, ONLY       :  I_PA, J_PA, PA, I_PL, J_PL, PL, I_PR, J_PR, PR, I_GOA, J_GOA, GOA, I_GOAt, J_GOAt, GOAt
      USE SPARSE_MATRICES, ONLY       :  SYM_PA, SYM_PL, SYM_PR

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)) :: SUBR_NAME = 'REDUCE_PA_TO_PL'

      INTEGER(LONG), INTENT(IN)        :: PART_VEC_A_LR(NDOFA)! Partitioning vector (F set into A and O sets)
      INTEGER(LONG), INTENT(IN)        :: PART_VEC_SUB(NSUB)  ! Partitioning vector (1's for all subcases)
      INTEGER(LONG), PARAMETER         :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER         :: NUM2        = 2     ! Used in subr's that partition matrices
      INTEGER(LONG)                    :: PL_ROW_MAX_TERMS    ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                    :: PR_ROW_MAX_TERMS    ! Output from subr PARTITION_SIZE (max terms in any row of matrix)




! **********************************************************************************************************************************
! Partition PL from PA

      IF (NDOFL > 0) THEN

         CALL PARTITION_SS_NTERM ( 'PA' , NTERM_PA , NDOFA, NSUB , SYM_PA , I_PA , J_PA      , PART_VEC_A_LR, PART_VEC_SUB,        &
                                    NUM1, NUM1, PL_ROW_MAX_TERMS, 'PL', NTERM_PL, SYM_PL )

         CALL ALLOCATE_SPARSE_MAT ( 'PL', NDOFL, NTERM_PL, SUBR_NAME )

         IF (NTERM_PL  > 0) THEN
            CALL PARTITION_SS ( 'PA' , NTERM_PA , NDOFA, NSUB , SYM_PA , I_PA , J_PA , PA , PART_VEC_A_LR, PART_VEC_SUB,           &
                                 NUM1, NUM1, PL_ROW_MAX_TERMS, 'PL', NTERM_PL , NDOFL, SYM_PL, I_PL , J_PL , PL  )
         ENDIF

      ENDIF

! Partition PR from PA

      IF (NDOFR > 0) THEN

         CALL PARTITION_SS_NTERM ( 'PA' , NTERM_PA , NDOFA, NSUB , SYM_PA , I_PA , J_PA ,      PART_VEC_A_LR, PART_VEC_SUB,        &
                                    NUM2, NUM1, PR_ROW_MAX_TERMS, 'PR', NTERM_PR, SYM_PR )

         CALL ALLOCATE_SPARSE_MAT ( 'PR', NDOFR, NTERM_PR, SUBR_NAME )

         IF (NTERM_PR  > 0) THEN
            CALL PARTITION_SS ( 'PA' , NTERM_PA , NDOFA, NSUB , SYM_PA , I_PA , J_PA , PA , PART_VEC_A_LR, PART_VEC_SUB,           &
                                 NUM2, NUM1, PR_ROW_MAX_TERMS, 'PR', NTERM_PR , NDOFR, SYM_PR, I_PR , J_PR , PR  )
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_PA_TO_PL


   END MODULE REDUCTION_A_TO_L
