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

   MODULE REDUCTION_F_TO_A

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: REDUCE_F_AO

   CONTAINS

      SUBROUTINE REDUCE_F_AO

! Call routines to reduce stiffness, mass, loads from F-set to A, O-sets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE CONSTANTS_1, ONLY           :  ZERO
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, LINKNO, KOO_SDIA, NDOFF, NDOFG, NDOFA, NDOFO, NSUB, SOL_NAME,               &
                                         NTERM_KFF , NTERM_KAA , NTERM_KAO , NTERM_KOO ,                                           &
                                         NTERM_KFFD, NTERM_KAAD, NTERM_KAOD, NTERM_KOOD,                                           &
                                         NTERM_MFF , NTERM_MAA , NTERM_MAO , NTERM_MOO ,                                           &
                                         NTERM_PF  , NTERM_PA  , NTERM_PO  , NTERM_GOA
      USE PARAMS, ONLY                :  EQCHK_OUTPUT, MATSPARS, PRTSTIFD, PRTSTIFF, PRTMASS, PRTFOR, SOLLIB, SPARSE_FLAVOR
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE DOF_TABLES, ONLY            :  TDOFI
      USE RIGID_BODY_DISP_MATS, ONLY  :  RBGLOBAL_GSET, RBGLOBAL_FSET, RBGLOBAL_ASET
      USE SPARSE_MATRICES, ONLY       :  I_KFF , J_KFF , KFF , I_KAA , J_KAA , KAA , I_KAO , J_KAO , KAO , I_KOO , J_KOO , KOO ,   &
                                         I_KFFD, J_KFFD, KFFD, I_KAAD, J_KAAD, KAAD, I_KAOD, J_KAOD, KAOD, I_KOOD, J_KOOD, KOOD,   &
                                         I_MFF , J_MFF , MFF , I_MAA , J_MAA , MAA , I_MAO , J_MAO , MAO , I_MOO , J_MOO , MOO ,   &
                                         I_PF  , J_PF  , PF  , I_PA  , J_PA  , PA  , I_PO  , J_PO  , PO
      USE SPARSE_MATRICES, ONLY       :  SYM_KAA
      USE SCRATCH_MATRICES
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, SLU_INFO

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_VEC
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE LAPACK_MATRIX_LIFECYCLE, ONLY:  DEALLOCATE_LAPACK_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  GET_MATRIX_DIAG_STATS
      USE RIGID_BODY_STORAGE_LIFECYCLE, ONLY:  ALLOCATE_RBGLOBAL, DEALLOCATE_RBGLOBAL
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE REDUCTION_CHECKS, ONLY      :  STIFF_MAT_EQUIL_CHK

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_F_AO'
      CHARACTER(132*BYTE)             :: MATRIX_NAME         ! Name of matrix for printout
      CHARACTER(44*BYTE)              :: MODNAM              ! Name to write to screen to describe module being run

      INTEGER(LONG)                   :: A_SET_COL           ! Col no. in array TDOFI where the A-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: A_SET_DOF           ! A-set DOF number
      INTEGER(LONG)                   :: DO_WHICH_CODE_FRAG    ! 1 or 2 depending on which seg of code to run (depends on BUCKLING)
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG)                   :: PART_VEC_F_AO(NDOFF)! Partitioning vector (G set into N and M sets)
      INTEGER(LONG)                   :: PART_VEC_SUB(NSUB)  ! Partitioning vector (1's for all subcases)


      REAL(DOUBLE)                    :: DUM_COL(NDOFO)      ! Temp variable used in SuperLU
      REAL(DOUBLE)                    :: KAA_DIAG(NDOFA)     ! Diagonal terms from KAA
      REAL(DOUBLE)                    :: KAA_MAX_DIAG        ! Max diag term  from KAA
      REAL(DOUBLE)                    :: KAAD_DIAG(NDOFA)    ! Diagonal terms from KAAD
      REAL(DOUBLE)                    :: KAAD_MAX_DIAG       ! Max diag term  from KAAD

      INTRINSIC                       :: DABS



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

! If there is an O-set, reduce KFF to KAA, MFF to MAA, PF to PA using UO = GOA*UA + UO0, where GOA = -KOO(-1)*KAO', UO0 = KOO(-1)*PO
! If there is no O-set, then equate KAA to KFF, MAA to MFF, PA to PF

         IF (NDOFO > 0) THEN
                                                           ! First, need to create part vectors used in the reduction (if NDOFO > 0)
            CALL PARTITION_VEC (NDOFF,'F ','A ','O ',PART_VEC_F_AO)

            DO I=1,NSUB
               PART_VEC_SUB = 1
            ENDDO

            CALL OURTIM                                    ! Reduce KFF to KAA
            IF (MATSPARS == 'Y') THEN
               MODNAM = 'REDUCE KFF TO KAA (SPARSE MATRIX ROUTINES)'
            ELSE
               MODNAM = 'REDUCE KFF TO KAA (FULL MATRIX ROUTINES)'
            ENDIF
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            CALL REDUCE_KFF_TO_KAA ( PART_VEC_F_AO )

            CALL OURTIM                                    ! Reduce MFF to MAA
            IF (MATSPARS == 'Y') THEN
               MODNAM = 'REDUCE MFF TO MAA (SPARSE MATRIX ROUTINES)'
            ELSE
               MODNAM = 'REDUCE MFF TO MAA (FULL MATRIX ROUTINES)'
            ENDIF
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            CALL REDUCE_MFF_TO_MAA ( PART_VEC_F_AO )


            IF ((SOL_NAME(1:5) /= 'MODES') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')) THEN

               IF (NTERM_PF > 0) THEN                      ! Reduce PF to PA

                  CALL OURTIM
                  IF (MATSPARS == 'Y') THEN
                     MODNAM = 'REDUCE PF  TO PA  (SPARSE MATRIX ROUTINES)'
                  ELSE
                     MODNAM = 'REDUCE PF  TO PA  (FULL MATRIX ROUTINES)'
                  ENDIF
                  WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

                  CALL REDUCE_PF_TO_PA ( PART_VEC_F_AO, PART_VEC_SUB )

               ELSE

                  NTERM_PA = 0
                  NTERM_PO = 0
                  CALL ALLOCATE_SPARSE_MAT ( 'PA', NDOFA, NTERM_PA, SUBR_NAME )

               ENDIF

            ENDIF

FreeS:      IF (SOLLIB == 'SPARSE  ') THEN                       ! Last, free the storage allocated inside SuperLU

               IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

                  DO J=1,NDOFO
                     DUM_COL(J) = ZERO
                  ENDDO

                  CALL C_FORTRAN_DGSSV( 3, NDOFO, NTERM_KOO, 1, KOO , I_KOO , J_KOO , DUM_COL, NDOFO, SLU_FACTORS, SLU_INFO )

                  IF (SLU_INFO .EQ. 0) THEN
                     WRITE (*,*) 'SUPERLU STORAGE FREED'
                  ELSE
                     WRITE(*,*) 'SUPERLU STORAGE NOT FREED. INFO FROM SUPERLU FREE STORAGE ROUTINE = ', SLU_INFO
                  ENDIF

               ENDIF

            ENDIF FreeS

         ELSE                                              ! There is no O-set, so equate F, A sets

            CALL OURTIM
            MODNAM = 'EQUATING A-SET TO F-SET'
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            NDOFA     = NDOFF

            NTERM_KAA = NTERM_KFF
            NTERM_KAO = 0
            NTERM_KOO = 0

            NTERM_MAA = NTERM_MFF
            NTERM_MAO = 0
            NTERM_MOO = 0

            NTERM_PA  = NTERM_PF
            NTERM_PO  = 0

            CALL ALLOCATE_SPARSE_MAT ( 'KAA', NDOFA, NTERM_KAA, SUBR_NAME )

!xx      DO I=1,NDOFA+1
!xx         I_KAA(I) = I_KFF(I)
!xx      ENDDO
!xx
!xx      DO I=1,NTERM_KAA
!xx         J_KAA(I) = J_KFF(I)
!xx           KAA(I) =   KFF(I)
!xx      ENDDO

            CALL ALLOCATE_SPARSE_MAT ( 'MAA', NDOFA, NTERM_MAA, SUBR_NAME )

!xx      DO I=1,NDOFA+1
!xx         I_MAA(I) = I_MFF(I)
!xx      ENDDO
!xx
!xx      DO I=1,NTERM_MAA
!xx         J_MAA(I) = J_MFF(I)
!xx           MAA(I) =   MFF(I)
!xx      ENDDO

            IF ((SOL_NAME(1:5) /= 'MODES') .AND. (SOL_NAME(1:12) /= 'GEN CB MODEL')) THEN

               CALL ALLOCATE_SPARSE_MAT ( 'PA', NDOFA, NTERM_PA, SUBR_NAME )

!xx         DO I=1,NDOFA+1
!xx            I_PA(I)  = I_PF(I)
!xx         ENDDO
!xx
!xx         DO I=1,NTERM_PA
!xx            J_PA(I) = J_PF(I)
!xx              PA(I) =   PF(I)
!xx         ENDDO

            ENDIF

         ENDIF

! Deallocate F set arrays

         MODNAM = 'DEALLOCATE F-SET ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KFF  ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KFF' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MFF  ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MFF' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PF   ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PF' )
         WRITE(SC1,*) CR13

! Deallocate ABAND (which was KOO before decomp and TOO later) and GOA

         MODNAM = 'DEALLOCATE GOA, ABAND ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GOA  ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'GOA' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate ABAND', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )
         WRITE(SC1,*) CR13

         IF (NDOFO == 0) THEN
            NTERM_GOA = 0
         ENDIF

! Print out stiffness matrix partitions, if requested

         IF (( PRTSTIFF(4) == 1) .OR. ( PRTSTIFF(4) == 3)) THEN
            IF (NTERM_KAA > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KAA'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'A ', 'A ', NTERM_KAA, NDOFA, I_KAA, J_KAA, KAA )
            ENDIF
         ENDIF

         IF (( PRTSTIFF(4) == 2) .OR. ( PRTSTIFF(4) == 3)) THEN
            IF (NTERM_KAO > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KAO'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'A ', 'O ', NTERM_KAO, NDOFA, I_KAO, J_KAO, KAO )
            ENDIF
            IF (NTERM_KOO > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KOO'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'O ', 'O ', NTERM_KOO, NDOFO, I_KOO, J_KOO, KOO )
            ENDIF
         ENDIF

         MODNAM = 'DEALLOCATE O SET ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KAO', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KAO' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KOO', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KOO' )
         WRITE(SC1,*) CR13

! Write matrix diagonal and stats, if requested.
! NOTE: call this subr even if PRTSTIFFD(4) = 0 since we need KAA_DIAG, KAA_MAX_DIAG for the equilibrium check

         CALL GET_MATRIX_DIAG_STATS ( 'KAA', 'A ', NDOFA, NTERM_KAA, I_KAA, J_KAA, KAA, PRTSTIFD(4), KAA_DIAG, KAA_MAX_DIAG )

! Print out mass matrix partitions, if requested

         IF (( PRTMASS(4) == 1) .OR. ( PRTMASS(4) == 3)) THEN
            IF (NTERM_MAA > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MAA'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'A ', 'A ', NTERM_MAA, NDOFA, I_MAA, J_MAA, MAA )
            ENDIF
         ENDIF

         IF (( PRTMASS(4) == 2) .OR. ( PRTMASS(4) == 3)) THEN
            IF (NTERM_MAO > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MAO'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'A ', 'O ', NTERM_MAO, NDOFA, I_MAO, J_MAO, MAO )
            ENDIF
            IF (NTERM_MOO > 0) THEN
               MATRIX_NAME = 'MASS MATRIX MOO'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'O ', 'O ', NTERM_MOO, NDOFO, I_MOO, J_MOO, MOO )
            ENDIF
         ENDIF

         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MAO', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MAO' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MOO', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'MOO' )
         WRITE(SC1,*) CR13

! Print out load matrix partitions, if requested

         IF (( PRTFOR(4) == 1) .OR. ( PRTFOR(4) == 3)) THEN
            IF (NTERM_PA  > 0) THEN
               MATRIX_NAME = 'LOAD MATRIX PA'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'A ', 'SUBCASE', NTERM_PA, NDOFA, I_PA, J_PA, PA )
            ENDIF
         ENDIF

         IF (( PRTFOR(4) == 2) .OR. ( PRTFOR(4) == 3)) THEN
            IF (NTERM_PO  > 0) THEN
               MATRIX_NAME = 'LOAD MATRIX PO'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'O ', 'SUBCASE', NTERM_PO, NDOFO, I_PO, J_PO, PO )
            ENDIF
         ENDIF

         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PO ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'PO' )
         WRITE(SC1,*) CR13

! Do equilibrium check on the A-set stiffness matrix, if requested

         IF ((EQCHK_OUTPUT(4) > 0) .OR. (EQCHK_OUTPUT(5) > 0)) THEN
            CALL ALLOCATE_RBGLOBAL ( 'A ', SUBR_NAME )
            IF (NDOFO > 0) THEN
               CALL TDOF_COL_NUM ( 'A ', A_SET_COL )
               DO I=1,NDOFG
                  A_SET_DOF = TDOFI(I,A_SET_COL)
                  IF (A_SET_DOF > 0) THEN
                     DO J=1,6
                        RBGLOBAL_ASET(A_SET_DOF,J) = RBGLOBAL_GSET(I,J)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               DO I=1,NDOFA
                  DO J=1,6
                     RBGLOBAL_ASET(I,J) = RBGLOBAL_FSET(I,J)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF

         IF (EQCHK_OUTPUT(4) > 0) THEN
            CALL OURTIM
            MODNAM = 'EQUILIBRIUM CHECK ON KAA                '
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
            CALL STIFF_MAT_EQUIL_CHK ( EQCHK_OUTPUT(4),'A ', SYM_KAA, NDOFA, NTERM_KAA, I_KAA, J_KAA, KAA, KAA_DIAG, KAA_MAX_DIAG, &
                                       RBGLOBAL_ASET)
         ENDIF
!xx      CALL DEALLOCATE_RBGLOBAL ( 'A ' )
         CALL DEALLOCATE_RBGLOBAL ( 'F ' )

! **********************************************************************************************************************************
      ELSE                                                 ! This is BUCKLING with LOAD_ISTEP = 2 (eigen part of BUCKLING)

         IF (NDOFO > 0) THEN
                                                           ! First, need to create part vectors used in the reduction (if NDOFO > 0)
            CALL PARTITION_VEC (NDOFF,'F ','A ','O ',PART_VEC_F_AO)

            DO I=1,NSUB
               PART_VEC_SUB = 1
            ENDDO

            CALL OURTIM                                       ! Reduce KFF to KAA
            IF (MATSPARS == 'Y') THEN
               MODNAM = 'REDUCE KFFD TO KAAD (SPARSE MATRIX ROUTINES)'
            ELSE
               MODNAM = 'REDUCE KFFD TO KAAD (FULL MATRIX ROUTINES)'
            ENDIF
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            CALL REDUCE_KFFD_TO_KAAD ( PART_VEC_F_AO )

         ELSE                                                 ! There is no O-set, so equate F, A sets

            CALL OURTIM
            MODNAM = 'EQUATING A-SET TO F-SET'
            WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC

            NDOFA     = NDOFF

            NTERM_KAAD = NTERM_KFFD
            NTERM_KAOD = 0
            NTERM_KOOD = 0

            CALL ALLOCATE_SPARSE_MAT ( 'KAAD', NDOFA, NTERM_KAAD, SUBR_NAME )

         ENDIF

! Deallocate F set arrays

         MODNAM = 'DEALLOCATE F-SET ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KFFD  ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KFFD' )
         WRITE(SC1,*) CR13

! Deallocate ABAND (which was KOO before decomp and TOO later) and GOA

         MODNAM = 'DEALLOCATE GOA, ABAND ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GOA  ', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'GOA' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate ABAND', CR13   ;   CALL DEALLOCATE_LAPACK_MAT ( 'ABAND' )
         WRITE(SC1,*) CR13

         IF (NDOFO == 0) THEN
            NTERM_GOA = 0
         ENDIF

! Print out stiffness matrix partitions, if requested

         IF (( PRTSTIFF(4) == 1) .OR. ( PRTSTIFF(4) == 3)) THEN
            IF (NTERM_KAAD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KAAD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'A ', 'A ', NTERM_KAAD, NDOFA, I_KAAD, J_KAAD, KAAD )
            ENDIF
         ENDIF

         IF (( PRTSTIFF(4) == 2) .OR. ( PRTSTIFF(4) == 3)) THEN
            IF (NTERM_KAOD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KAOD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'A ', 'O ', NTERM_KAOD, NDOFA, I_KAOD, J_KAOD, KAOD )
            ENDIF
            IF (NTERM_KOOD > 0) THEN
               MATRIX_NAME = 'STIFFNESS MATRIX KOOD'
               CALL WRITE_SPARSE_CRS ( MATRIX_NAME, 'O ', 'O ', NTERM_KOOD, NDOFO, I_KOOD, J_KOOD, KOOD )
            ENDIF
         ENDIF

         MODNAM = 'DEALLOCATE O SET ARRAYS'
         WRITE(SC1,2092) MODNAM,HOUR,MINUTE,SEC,SFRAC
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KAOD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KAOD' )
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KOOD', CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'KOOD' )
         WRITE(SC1,*) CR13

! Write matrix diagonal and stats, if requested.
! NOTE: call this subr even if PRTSTIFFD(4) = 0 since we need KAA_DIAG, KAA_MAX_DIAG for the equilibrium check

         CALL GET_MATRIX_DIAG_STATS ( 'KAAD', 'A ', NDOFA, NTERM_KAAD, I_KAAD, J_KAAD, KAAD, PRTSTIFD(4), KAAD_DIAG, KAAD_MAX_DIAG )

      ENDIF



      RETURN

! **********************************************************************************************************************************
 2092 FORMAT(6X,A44,18X,I2,':',I2,':',I2,'.',I3)

12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_F_AO


      SUBROUTINE REDUCE_KFFD_TO_KAAD ( PART_VEC_F_AO )

! Call routines to reduce the KFFD differential stiffness matrix from the F-set to the A, O-sets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L2E, L2ESTAT, LINK2E, L2E_MSG, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FACTORED_MATRIX, FATAL_ERR, NDOFF, NDOFA, NDOFO, NTERM_KFFD, NTERM_KAAD,    &
                                         NTERM_KAOD, NTERM_KOOD, NTERM_KOODs, NTERM_GOA
      USE PARAMS, ONLY                :  EPSIL, KOORAT, SPARSTOR, RCONDK
      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE SPARSE_MATRICES, ONLY       :  I_KFFD, J_KFFD, KFFD, I_KAAD, J_KAAD, KAAD, I_KAOD, J_KAOD, KAOD, I_GOA, J_GOA, GOA,      &
                                         I_KOOD, I2_KOOD, J_KOOD, KOOD, I_KOODs, I2_KOODs, J_KOODs, KOODs

      USE SPARSE_MATRICES, ONLY       :  SYM_GOA, SYM_KFFD, SYM_KAAD, SYM_KAOD, SYM_KOOD
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  READ_MATRIX_1
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_NONSYM_TO_CRS_SYM, SPARSE_CRS_SPARSE_CCS
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_CRS_TERM_COUNT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_KFFD_TO_KAAD'
      CHARACTER(  1*BYTE)             :: SYM_CRS2            ! Storage format for matrix CRS2 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_F_AO(NDOFF)! Partitioning vector (F set into A and O sets)
      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS

      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG)                   :: KAAD_ROW_MAX_TERMS  ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KAOD_ROW_MAX_TERMS  ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KOOD_ROW_MAX_TERMS  ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: NTERM_CRS1          ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2          ! Number of terms in matrix CRS2
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices


      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Partition KAAD from KFFD (This is KAAD before reduction, or KAAD(bar) )

      IF (NDOFA > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KFFD', NTERM_KFFD, NDOFF, NDOFF, SYM_KFFD, I_KFFD, J_KFFD,      PART_VEC_F_AO, PART_VEC_F_AO,  &
                                    NUM1, NUM1, KAAD_ROW_MAX_TERMS, 'KAAD', NTERM_KAAD, SYM_KAAD )

         CALL ALLOCATE_SPARSE_MAT ( 'KAAD', NDOFA, NTERM_KAAD, SUBR_NAME )

         IF (NTERM_KAAD > 0) THEN
            CALL PARTITION_SS ( 'KFFD', NTERM_KFFD, NDOFF, NDOFF, SYM_KFFD, I_KFFD, J_KFFD, KFFD, PART_VEC_F_AO, PART_VEC_F_AO,    &
                                 NUM1, NUM1, KAAD_ROW_MAX_TERMS, 'KAAD', NTERM_KAAD, NDOFA, SYM_KAAD, I_KAAD, J_KAAD, KAAD )
         ENDIF

      ENDIF

! Partition KAOD from KFFD

      IF ((NDOFA > 0) .AND. (NDOFO > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'KFFD', NTERM_KFFD, NDOFF, NDOFF, SYM_KFFD, I_KFFD, J_KFFD,      PART_VEC_F_AO, PART_VEC_F_AO,  &
                                    NUM1, NUM2, KAOD_ROW_MAX_TERMS, 'KAOD', NTERM_KAOD, SYM_KAOD )

         CALL ALLOCATE_SPARSE_MAT ( 'KAOD', NDOFA, NTERM_KAOD, SUBR_NAME )

         IF (NTERM_KAOD > 0) THEN
            CALL PARTITION_SS ( 'KFFD', NTERM_KFFD, NDOFF, NDOFF, SYM_KFFD, I_KFFD, J_KFFD, KFFD, PART_VEC_F_AO, PART_VEC_F_AO,    &
                                 NUM1, NUM2, KAOD_ROW_MAX_TERMS, 'KAOD', NTERM_KAOD, NDOFA, SYM_KAOD, I_KAOD, J_KAOD, KAOD )
         ENDIF

      ENDIF

! Partition KOOD from KFFD

      IF (NDOFO > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KFFD', NTERM_KFFD, NDOFF, NDOFF, SYM_KFFD, I_KFFD, J_KFFD,      PART_VEC_F_AO, PART_VEC_F_AO,  &
                                    NUM2, NUM2, KOOD_ROW_MAX_TERMS, 'KOOD', NTERM_KOOD, SYM_KOOD )

         CALL ALLOCATE_SPARSE_MAT ( 'KOOD', NDOFO, NTERM_KOOD, SUBR_NAME )

         IF (NTERM_KOOD > 0) THEN
            CALL PARTITION_SS ( 'KFFD', NTERM_KFFD, NDOFF, NDOFF, SYM_KFFD, I_KFFD, J_KFFD, KFFD, PART_VEC_F_AO, PART_VEC_F_AO,    &
                                 NUM2, NUM2, KOOD_ROW_MAX_TERMS, 'KOOD', NTERM_KOOD, NDOFO, SYM_KOOD, I_KOOD, J_KOOD, KOOD )
         ENDIF

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Reduce KFFD to KAAD = KAAD(bar) + KAOD*GOA. Note: (KAOD*GOA)' + GOA'*KOOD*GOA = 0 by definition of GOA = -KOOD(-1)*KAOD'
! If PARAM MATSPARS = 'Y', then we use sparse matrix operations (multiply/add/transpose). If not, use full matrix operations

         IF (NTERM_KAOD > 0) THEN                          ! Calc KAOD*GOA & add it orig KAAD

            IF (.NOT. ALLOCATED(GOA)) THEN                 ! Self-load GOA from L2E if not in memory (needed for re-entrant call)
               CALL ALLOCATE_SPARSE_MAT ( 'GOA', NDOFO, NTERM_GOA, SUBR_NAME )
               CALL READ_MATRIX_1 ( LINK2E, L2E, 'N', 'Y', 'KEEP', L2E_MSG,                                                        &
                                    'GOA', NTERM_GOA, 'Y', NDOFO, I_GOA, J_GOA, GOA )
            ENDIF

                                                           ! CCS1 will be sparse CCS format version of sparse CRS matrix GOA
            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFA, NTERM_GOA, SUBR_NAME )
            CALL SPARSE_CRS_SPARSE_CCS ( NDOFO, NDOFA, NTERM_GOA, 'GOA', I_GOA, J_GOA, GOA, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

                                                           ! Sparse multiply to get CRS1 = KAOD*GOA. Use CCS1 for GOA CCS
            CALL MATMULT_SSS_NTERM ( 'KAOD' , NDOFA, NTERM_KAOD, SYM_KAOD, I_KAOD , J_KAOD ,                                       &
                                     'GOA' , NDOFA, NTERM_GOA, SYM_GOA, J_CCS1, I_CCS1, AROW_MAX_TERMS,                            &
                                     'CRS1' ,       NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFA, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'KAOD' , NDOFA, NTERM_KAOD , SYM_KAOD, I_KAOD , J_KAOD , KAOD ,                                     &
                               'GOA' , NDOFA, NTERM_GOA , SYM_GOA, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

                                                           ! CRS1 = KAOD*GOA has all nonzero terms in it.
            IF      (SPARSTOR == 'SYM   ') THEN            !      If SPARSTOR == 'SYM   ', rewrite CRS1 as sym in CRS2

               CALL SPARSE_CRS_TERM_COUNT ( NDOFA, NTERM_CRS1, 'CRS1 = KAOD*GOA all nonzeros', I_CRS1, J_CRS1, NTERM_CRS2 )
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFA, NTERM_CRS2, SUBR_NAME )
               CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS1 = KAOD*GOA all nonzeros', NDOFA, NTERM_CRS1, I_CRS1, J_CRS1, CRS1,               &
                                            'CRS2 = KAOD*GOA stored sym'  ,        NTERM_CRS2, I_CRS2, J_CRS2, CRS2 )
               SYM_CRS2 = 'Y'

            ELSE IF (SPARSTOR == 'NONSYM') THEN            !      If SPARSTOR == 'NONSYM', rewrite CRS1 in CRS2 with NTERM_CRS2

               NTERM_CRS2 = NTERM_CRS1
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFA, NTERM_CRS2, SUBR_NAME )
               DO I=1,NDOFA+1
                  I_CRS2(I) = I_CRS1(I)
               ENDDO
               DO I=1,NTERM_CRS2
                  J_CRS2(I) = J_CRS1(I)
                    CRS2(I) =   CRS1(I)
               ENDDO
!xx            SYM_CRS2 = 'Y'      ! NO: SHOULD BE 'N' HERE SINCE SPARSTOR = 'NONSYM' IN THIS BRANCH OF THE IF BLOCK
               SYM_CRS2 = 'N'

            ELSE                                           ! Error - incorrect SPARSTOR

               WRITE(ERR,932) SUBR_NAME, SPARSTOR
               WRITE(F06,932) SUBR_NAME, SPARSTOR
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! Deallocate CRS1
                                                           ! Sparse add to get CRS1 = KAAD-bar + CRS2
            CALL MATADD_SSS_NTERM ( NDOFA, 'KAAD-bar', NTERM_KAAD, I_KAAD , J_KAAD , SYM_KAAD , 'KOOD*GOA', NTERM_CRS2,            &
                                                                 I_CRS2, J_CRS2, SYM_CRS2, 'CRS1', NTERM_CRS1 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFA, NTERM_CRS1, SUBR_NAME )
            CALL MATADD_SSS (NDOFA, 'KAAD-bar' , NTERM_KAAD , I_KAAD , J_KAAD , KAAD , ONE, 'KOOD*GOA', NTERM_CRS2,                &
                             I_CRS2, J_CRS2, CRS2, ONE, 'CRS1', NTERM_CRS1, I_CRS1, J_CRS1, CRS1 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! Deallocate CRS2

            NTERM_KAAD = NTERM_CRS1                        ! Reallocate KAAD to be size of CRS1
            WRITE(SC1, * ) '    Reallocate KAAD'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KAAD', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KAAD' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KAAD', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KAAD', NDOFA, NTERM_KAAD, SUBR_NAME )
                                                           ! Set KAAD = CRS1
            DO I=1,NDOFA+1
               I_KAAD(I) = I_CRS1(I)
            ENDDO
            DO J=1,NTERM_KAAD
               J_KAAD(J) = J_CRS1(J)
                 KAAD(J) =   CRS1(J)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! Deallocate CRS1

         ENDIF



      RETURN

! **********************************************************************************************************************************
  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

 2400 FORMAT(' *ERROR  2400: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THERE IS AN O-SET BUT GUYAN REDUCTION MATRIX GOA HAS ',I12,' TERMS IN IT. MUST BE > 0')

12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_KFFD_TO_KAAD


      SUBROUTINE REDUCE_KFF_TO_KAA ( PART_VEC_F_AO )

! Call routines to reduce the KFF linear stiffness matrix from the F-set to the A, O-sets

! NOTE: This subr has code for sparse matrices as well as full matrices, (i.e. Bulk Data PARAM MATSPARS = 'Y' for sparse and 'N'
! for full). The code for full matrices was put in originally in order that the sparse code could be thoroughly checked. That task
! is complete and the remaining full matrix code has not been maintained since around 2005. In addition, new capability added to
! MYSTRAN since that approx time does not have full matrix code.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L2E, LINK2E, L2E_MSG, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, KOO_SDIA, NDOFF, NDOFA, NDOFO, NTERM_KFF,       &
                                         NTERM_KAA, NTERM_KAO, NTERM_KOO, NTERM_GOA
      USE PARAMS, ONLY                :  KOORAT, MATSPARS, SOLLIB, SPARSTOR, SPARSE_FLAVOR, RCONDK
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE FULL_MATRICES, ONLY         :  KAA_FULL, KAO_FULL, GOA_FULL, DUM1, DUM2
      USE SPARSE_MATRICES, ONLY       :  I_KFF, J_KFF, KFF, I_KAA, J_KAA, KAA, I_KAO, J_KAO, KAO, I_GOA, J_GOA, GOA,               &
                                         I_KOO, J_KOO, KOO

      USE SPARSE_MATRICES, ONLY       :  SYM_GOA, SYM_KFF, SYM_KAA, SYM_KAO, SYM_KOO
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE LAPACK_ADAPTERS, ONLY       :  SYM_MAT_DECOMP_LAPACK
      USE SUPERLU_ADAPTERS, ONLY      :  SYM_MAT_DECOMP_SUPRLU
      USE SuperLU_STUF, ONLY          :  SLU_DIAG_RATIO
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE MATRIX_FILE_IO, ONLY        :  WRITE_MATRIX_1
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_NONSYM_TO_CRS_SYM, SPARSE_CRS_SPARSE_CCS, SPARSE_CRS_TO_FULL
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_CRS_TERM_COUNT
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATADD_FFF, MATMULT_FFF
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CNT_NONZ_IN_FULL_MAT

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_KFF_TO_KAA'
      CHARACTER(  1*BYTE)             :: CLOSE_IT            ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER(  8*BYTE)             :: CLOSE_STAT          ! Char constant for the CLOSE status of a file
      CHARACTER(  1*BYTE)             :: EQUED               ! 'Y' if the stiff matrix was equilibrated in subr EQUILIBRATE
      CHARACTER(  1*BYTE)             :: EQUIL_KOO           ! 'Y'/'N' for whether to equilibrate KOO in subr SYM_MAT_DECOMP_LAPACK
      CHARACTER(  1*BYTE)             :: SYM_CRS2            ! Storage format for matrix CRS2 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_F_AO(NDOFF)! Partitioning vector (F set into A and O sets)
      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS

      INTEGER(LONG)                   :: DEB_PRT(2)          ! Debug numbers to say whether to write ABAND and/or its decomp
!                                                              in called subr SYM_MAT_DECOMP_LAPACK

      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG)                   :: INFO        = 0     ! Input value for subr SYM_MAT_DECOMP_LAPACK (quit on sing KRRCB)
      INTEGER(LONG), PARAMETER        :: ITRNSPB     = 0     ! Transpose indicator for matrix multiply routine
      INTEGER(LONG)                   :: KAA_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KAO_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: KOO_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: NTERM_CRS1          ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2          ! Number of terms in matrix CRS2
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices


      REAL(DOUBLE)                    :: ALPHA = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: BETA  = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: K_INORM             ! Inf norm of KOO matrix

      REAL(DOUBLE)                    :: KOO_SCALE_FACS(NDOFO)
                                                             ! KOO equilibration scale factors

      REAL(DOUBLE)                    :: RCOND               ! Recrip of cond no. of the KLL. Det in  subr COND_NUM
      REAL(DOUBLE)                    :: SMALL               ! A number used in filtering out small numbers from a full matrix

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Initialize outputs

      KOO_SDIA   = 0

! Partition KAA from KFF (This is KAA before reduction, or KAA(bar) )

      IF (NDOFA > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KFF', NTERM_KFF, NDOFF, NDOFF, SYM_KFF, I_KFF, J_KFF,      PART_VEC_F_AO, PART_VEC_F_AO,       &
                                    NUM1, NUM1, KAA_ROW_MAX_TERMS, 'KAA', NTERM_KAA, SYM_KAA )

         CALL ALLOCATE_SPARSE_MAT ( 'KAA', NDOFA, NTERM_KAA, SUBR_NAME )

         IF (NTERM_KAA > 0) THEN
            CALL PARTITION_SS ( 'KFF', NTERM_KFF, NDOFF, NDOFF, SYM_KFF, I_KFF, J_KFF, KFF, PART_VEC_F_AO, PART_VEC_F_AO,          &
                                 NUM1, NUM1, KAA_ROW_MAX_TERMS, 'KAA', NTERM_KAA, NDOFA, SYM_KAA, I_KAA, J_KAA, KAA )
         ENDIF

      ENDIF

! Partition KAO from KFF

      IF ((NDOFA > 0) .AND. (NDOFO > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'KFF', NTERM_KFF, NDOFF, NDOFF, SYM_KFF, I_KFF, J_KFF,      PART_VEC_F_AO, PART_VEC_F_AO,       &
                                    NUM1, NUM2, KAO_ROW_MAX_TERMS, 'KAO', NTERM_KAO, SYM_KAO )

         CALL ALLOCATE_SPARSE_MAT ( 'KAO', NDOFA, NTERM_KAO, SUBR_NAME )

         IF (NTERM_KAO > 0) THEN
            CALL PARTITION_SS ( 'KFF', NTERM_KFF, NDOFF, NDOFF, SYM_KFF, I_KFF, J_KFF, KFF, PART_VEC_F_AO, PART_VEC_F_AO,          &
                                 NUM1, NUM2, KAO_ROW_MAX_TERMS, 'KAO', NTERM_KAO, NDOFA, SYM_KAO, I_KAO, J_KAO, KAO )
         ENDIF

      ENDIF

! Partition KOO from KFF

      IF (NDOFO > 0) THEN

         CALL PARTITION_SS_NTERM ( 'KFF', NTERM_KFF, NDOFF, NDOFF, SYM_KFF, I_KFF, J_KFF,      PART_VEC_F_AO, PART_VEC_F_AO,       &
                                    NUM2, NUM2, KOO_ROW_MAX_TERMS, 'KOO', NTERM_KOO, SYM_KOO )

         CALL ALLOCATE_SPARSE_MAT ( 'KOO', NDOFO, NTERM_KOO, SUBR_NAME )

         IF (NTERM_KOO > 0) THEN
            CALL PARTITION_SS ( 'KFF', NTERM_KFF, NDOFF, NDOFF, SYM_KFF, I_KFF, J_KFF, KFF, PART_VEC_F_AO, PART_VEC_F_AO,          &
                                 NUM2, NUM2, KOO_ROW_MAX_TERMS, 'KOO', NTERM_KOO, NDOFO, SYM_KOO, I_KOO, J_KOO, KOO )
         ENDIF

      ENDIF

! Decompose KOO, if there are any O-set DOF's

      IF (NDOFO > 0) THEN

         DEB_PRT(1) = 28
         DEB_PRT(2) = 29

         IF (SOLLIB == 'BANDED  ') THEN                    ! Use LAPACK

            KOO_SDIA   = 0
            EQUIL_KOO  = 'N'
            INFO = 0
            CALL SYM_MAT_DECOMP_LAPACK ( SUBR_NAME, 'KOO', 'O ', NDOFO, NTERM_KOO, I_KOO, J_KOO, KOO, 'Y', KOORAT, EQUIL_KOO,      &
                                         RCONDK, DEB_PRT, EQUED, KOO_SDIA, K_INORM, RCOND, KOO_SCALE_FACS, INFO )
         ELSE IF (SOLLIB == 'SPARSE  ') THEN

            IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

               INFO = 0
               SLU_DIAG_RATIO = 'Y'                        ! Check matrix diag / factor diag (MAXRATIO), as on the LAPACK path
               CALL SYM_MAT_DECOMP_SUPRLU ( SUBR_NAME, 'KOO', 'O ', NDOFO, NTERM_KOO, I_KOO, J_KOO, KOO, INFO )
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

         ENDIF

      ENDIF

! Solve for GOA = -KOO(-1)*KAO'

      IF (NTERM_KAO > 0) THEN

         CALL SOLVE_GOA                                      ! Solve for GOA matrix
         CLOSE_IT   = 'Y'
         CLOSE_STAT = 'KEEP'
         CALL WRITE_MATRIX_1 ( LINK2E, L2E, CLOSE_IT, CLOSE_STAT, L2E_MSG, 'GOA', NTERM_GOA, NDOFO, I_GOA, J_GOA, GOA )

      ELSE

         NTERM_GOA = 0                                     ! GOA is null
         CALL ALLOCATE_SPARSE_MAT ( 'GOA', NDOFO, NTERM_GOA, SUBR_NAME )

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Reduce KFF to KAA = KAA(bar) + KAO*GOA. Note: (KAO*GOA)' + GOA'*KOO*GOA = 0 by definition of GOA = -KOO(-1)*KAO'
! If PARAM MATSPARS = 'Y', then we use sparse matrix operations (multiply/add/transpose). If not, use full matrix operations

      IF (MATSPARS == 'Y') THEN                            ! Reduce KFF to KAA using sparse matrix operations

         IF (NTERM_KAO > 0) THEN                           ! Calc KAO*GOA & add it orig KAA


                                                           ! CCS1 will be sparse CCS format version of sparse CRS matrix GOA
            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFA, NTERM_GOA, SUBR_NAME )
            CALL SPARSE_CRS_SPARSE_CCS ( NDOFO, NDOFA, NTERM_GOA, 'GOA', I_GOA, J_GOA, GOA, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

                                                           ! Sparse multiply to get CRS1 = KAO*GOA. Use CCS1 for GOA CCS
            CALL MATMULT_SSS_NTERM ( 'KAO' , NDOFA, NTERM_KAO, SYM_KAO, I_KAO , J_KAO ,                                            &
                                     'GOA' , NDOFA, NTERM_GOA, SYM_GOA, J_CCS1, I_CCS1, AROW_MAX_TERMS,                            &
                                     'CRS1' ,       NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFA, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'KAO' , NDOFA, NTERM_KAO , SYM_KAO, I_KAO , J_KAO , KAO ,                                           &
                               'GOA' , NDOFA, NTERM_GOA , SYM_GOA, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

                                                           ! CRS1 = KAO*GOA has all nonzero terms in it.
            IF      (SPARSTOR == 'SYM   ') THEN            !      If SPARSTOR == 'SYM   ', rewrite CRS1 as sym in CRS2

               CALL SPARSE_CRS_TERM_COUNT ( NDOFA, NTERM_CRS1, 'CRS1 = KAO*GOA all nonzeros', I_CRS1, J_CRS1, NTERM_CRS2 )
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFA, NTERM_CRS2, SUBR_NAME )
               CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS1 = KAO*GOA all nonzeros', NDOFA, NTERM_CRS1, I_CRS1, J_CRS1, CRS1,                &
                                            'CRS2 = KAO*GOA stored sym'  ,        NTERM_CRS2, I_CRS2, J_CRS2, CRS2 )
               SYM_CRS2 = 'Y'

            ELSE IF (SPARSTOR == 'NONSYM') THEN            !      If SPARSTOR == 'NONSYM', rewrite CRS1 in CRS2 with NTERM_CRS2

               NTERM_CRS2 = NTERM_CRS1
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFA, NTERM_CRS2, SUBR_NAME )
               DO I=1,NDOFA+1
                  I_CRS2(I) = I_CRS1(I)
               ENDDO
               DO I=1,NTERM_CRS2
                  J_CRS2(I) = J_CRS1(I)
                    CRS2(I) =   CRS1(I)
               ENDDO
               SYM_CRS2 = 'N'

            ELSE                                           ! Error - incorrect SPARSTOR

               WRITE(ERR,932) SUBR_NAME, SPARSTOR
               WRITE(F06,932) SUBR_NAME, SPARSTOR
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! Deallocate CRS1
                                                           ! Sparse add to get CRS1 = KAA-bar + CRS2
            CALL MATADD_SSS_NTERM ( NDOFA, 'KAA-bar', NTERM_KAA, I_KAA , J_KAA , SYM_KAA , 'KOO*GOA', NTERM_CRS2,                  &
                                                                 I_CRS2, J_CRS2, SYM_CRS2, 'CRS1', NTERM_CRS1 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFA, NTERM_CRS1, SUBR_NAME )
            CALL MATADD_SSS (NDOFA, 'KAA-bar' , NTERM_KAA , I_KAA , J_KAA , KAA , ONE, 'KOO*GOA', NTERM_CRS2, I_CRS2, J_CRS2, CRS2,&
                            ONE, 'CRS1', NTERM_CRS1, I_CRS1, J_CRS1, CRS1 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! Deallocate CRS2

            NTERM_KAA = NTERM_CRS1                         ! Reallocate KAA to be size of CRS1
            WRITE(SC1, * ) '    Reallocate KAA'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KAA', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'KAA' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KAA', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'KAA', NDOFA, NTERM_KAA, SUBR_NAME )
                                                           ! Set KAA = CRS1
            DO I=1,NDOFA+1
               I_KAA(I) = I_CRS1(I)
            ENDDO
            DO J=1,NTERM_KAA
               J_KAA(J) = J_CRS1(J)
                 KAA(J) =   CRS1(J)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! Deallocate CRS1

         ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (MATSPARS == 'N') THEN
                                                           ! Calc reduced KAA <= KAA + KAO*GOA
         CALL ALLOCATE_FULL_MAT ('KAA_FULL', NDOFA, NDOFA, SUBR_NAME )

         IF (NTERM_KAA > 0) THEN
            CALL SPARSE_CRS_TO_FULL ( 'KAA', NTERM_KAA, NDOFA, NDOFA, SYM_KAA, I_KAA, J_KAA, KAA, KAA_FULL )
         ENDIF

         IF (NTERM_KAO > 0) THEN                           ! Part 1: calc KAO*GOA and add it & it's transpose to KNN_FULL

            CALL ALLOCATE_FULL_MAT ( 'KAO_FULL', NDOFA, NDOFO, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'KAO', NTERM_KAO, NDOFA, NDOFO, SYM_KAO, I_KAO, J_KAO, KAO, KAO_FULL )

            CALL ALLOCATE_FULL_MAT ( 'GOA_FULL', NDOFO, NDOFA, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'GOA', NTERM_GOA, NDOFO, NDOFA, SYM_GOA, I_GOA, J_GOA, GOA, GOA_FULL )

            CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFA, NDOFA, SUBR_NAME )

            CALL MATMULT_FFF (  KAO_FULL, GOA_FULL, NDOFA, NDOFO, NDOFA, DUM1 )
            CALL DEALLOCATE_FULL_MAT ( 'KAO_FULL' )
            CALL DEALLOCATE_FULL_MAT ( 'GOA_FULL' )

            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFA, NDOFA, SUBR_NAME )
            CALL MATADD_FFF ( KAA_FULL, DUM1, NDOFA, NDOFA, ALPHA, BETA, ITRNSPB, DUM2 )
            DO I=1,NDOFA
               DO J=1,NDOFA
                  KAA_FULL(I,J) = DUM2(I,J)
               ENDDO
            ENDDO

            CALL DEALLOCATE_FULL_MAT ( 'DUM1' )
            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

         ENDIF

         CALL CNT_NONZ_IN_FULL_MAT ( 'KAA_FULL  ', KAA_FULL, NDOFA, NDOFA, SYM_KAA, NTERM_KAA, SMALL )

         WRITE(SC1, * ) '    Reallocate KAA'
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KAA', CR13
         CALL DEALLOCATE_SPARSE_MAT ( 'KAA' )
         WRITE(SC1,12345,ADVANCE='NO') '       Allocate   KAA', CR13
         CALL ALLOCATE_SPARSE_MAT ( 'KAA', NDOFA, NTERM_KAA, SUBR_NAME )

         IF (NTERM_KAA > 0) THEN                           ! Create new sparse arrays from KAA_FULL
!xx         CALL FULL_TO_SPARSE_CRS ( 'KAA_FULL  ', NDOFA, NDOFA, KAA_FULL, NTERM_KAA, SYM_KAA, I_KAA, J_KAA, KAA )
            CALL DEALLOCATE_FULL_MAT ( 'KAA_FULL' )
         ENDIF

      ELSE

         WRITE(ERR,911) SUBR_NAME,MATSPARS
         WRITE(F06,911) SUBR_NAME,MATSPARS
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
  911 FORMAT(' *ERROR   911: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER MATSPARS MUST BE EITHER ','Y',' OR ','N',' BUT VALUE IS ',A)

  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_KFF_TO_KAA


      SUBROUTINE REDUCE_MFF_TO_MAA ( PART_VEC_F_AO )

! Call routines to reduce the MFF mass matrix from the F-set to the A, O-sets. See Appendix B to the MYSTRAN User's Reference Manual
! for the derivation of the reduction equations.

! NOTE: This subr has code for sparse matrices as well as full matrices, (i.e. Bulk Data PARAM MATSPARS = 'Y' for sparse and 'N'
! for full). The code for full matrices was put in originally in order that the sparse code could be thoroughly checked. That task
! is complete and the remaining full matrix code has not been maintained since around 2005. In addition, new capability added to
! MYSTRAN since that approx time does not have full matrix code.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFF, NDOFA, NDOFO, NTERM_MFF, NTERM_MAA, NTERM_MAO, NTERM_MOO, &
                                         NTERM_GOA
      USE PARAMS, ONLY                :  EPSIL, MATSPARS, SPARSTOR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE SPARSE_MATRICES, ONLY       :  I_MFF, J_MFF, MFF, I_MAA, J_MAA, MAA, I_MAO, J_MAO, MAO, I_MOO, J_MOO, MOO,               &
                                         I_GOA, J_GOA, GOA,  I_GOAt, J_GOAt, GOAt
      USE SPARSE_MATRICES, ONLY       :  SYM_GOA, SYM_MFF, SYM_MAA, SYM_MAO, SYM_MOO
      USE FULL_MATRICES, ONLY         :  MAA_FULL, MAO_FULL, MOO_FULL, GOA_FULL, DUM1, DUM2, DUM3
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM, MATTRNSP_SS
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  CRS_NONSYM_TO_CRS_SYM, SPARSE_CRS_SPARSE_CCS, SPARSE_CRS_TO_FULL
      USE SPARSE_CRS_ACCESS, ONLY     :  SPARSE_CRS_TERM_COUNT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATADD_FFF, MATMULT_FFF, MATMULT_FFF_T
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CNT_NONZ_IN_FULL_MAT

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_MFF_TO_MAA'
      CHARACTER(  1*BYTE)             :: SYM_CRS1            ! Storage format for matrix CRS1 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)
      CHARACTER(  1*BYTE)             :: SYM_CRS3            ! Storage format for matrix CRS3 (either 'Y' for sym storage or
!                                                              'N' for nonsymmetric storage)

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_F_AO(NDOFF)! Partitioning vector (F set into A and O sets)
      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
!                                                              the ones on and above the diagonal (controlled by param SPARSTOR)
      INTEGER(LONG)                   :: ITRNSPB             ! Transpose indicator for matrix multiply routine
      INTEGER(LONG)                   :: MAA_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: MAO_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: MOO_ROW_MAX_TERMS   ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: NTERM_CCS1          ! Number of terms in matrix CCS1
      INTEGER(LONG)                   :: NTERM_CRS1          ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2          ! Number of terms in matrix CRS2
      INTEGER(LONG)                   :: NTERM_CRS3          ! Number of terms in matrix CRS3
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices


      REAL(DOUBLE)                    :: ALPHA = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: BETA  = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: SMALL             ! A number used in filtering out small numbers from a full matrix

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Partition MAA from MFF (This is MAA before reduction, or MAA(bar) )

      IF (NDOFA > 0) THEN

         CALL PARTITION_SS_NTERM ( 'MFF', NTERM_MFF, NDOFF, NDOFF, SYM_MFF, I_MFF, J_MFF,      PART_VEC_F_AO, PART_VEC_F_AO,       &
                                    NUM1, NUM1, MAA_ROW_MAX_TERMS, 'MAA', NTERM_MAA, SYM_MAA )

         CALL ALLOCATE_SPARSE_MAT ( 'MAA', NDOFA, NTERM_MAA, SUBR_NAME )

         IF (NTERM_MAA > 0) THEN
            CALL PARTITION_SS ( 'MFF', NTERM_MFF, NDOFF, NDOFF, SYM_MFF, I_MFF, J_MFF, MFF, PART_VEC_F_AO, PART_VEC_F_AO,          &
                                 NUM1, NUM1, MAA_ROW_MAX_TERMS, 'MAA', NTERM_MAA, NDOFA, SYM_MAA, I_MAA, J_MAA, MAA )
         ENDIF

      ENDIF

! Partition MAO from MFF

      IF ((NDOFA > 0) .AND. (NDOFO > 0)) THEN

         CALL PARTITION_SS_NTERM ( 'MFF', NTERM_MFF, NDOFF, NDOFF, SYM_MFF, I_MFF, J_MFF,      PART_VEC_F_AO, PART_VEC_F_AO,       &
                                    NUM1, NUM2, MAO_ROW_MAX_TERMS, 'MAO', NTERM_MAO, SYM_MAO )

         CALL ALLOCATE_SPARSE_MAT ( 'MAO', NDOFA, NTERM_MAO, SUBR_NAME )

         IF (NTERM_MAO > 0) THEN
            CALL PARTITION_SS ( 'MFF', NTERM_MFF, NDOFF, NDOFF, SYM_MFF, I_MFF, J_MFF, MFF, PART_VEC_F_AO, PART_VEC_F_AO,          &
                                 NUM1, NUM2, MAO_ROW_MAX_TERMS, 'MAO', NTERM_MAO, NDOFA, SYM_MAO, I_MAO, J_MAO, MAO )
         ENDIF

      ENDIF

! Partition MOO from MFF

      IF (NDOFO > 0) THEN

         CALL PARTITION_SS_NTERM ( 'MFF', NTERM_MFF, NDOFF, NDOFF, SYM_MFF, I_MFF, J_MFF,      PART_VEC_F_AO, PART_VEC_F_AO,       &
                                    NUM2, NUM2, MOO_ROW_MAX_TERMS, 'MOO', NTERM_MOO, SYM_MOO )

         CALL ALLOCATE_SPARSE_MAT ( 'MOO', NDOFO, NTERM_MOO, SUBR_NAME )

         IF (NTERM_MOO > 0) THEN
            CALL PARTITION_SS ( 'MFF', NTERM_MFF, NDOFF, NDOFF, SYM_MFF, I_MFF, J_MFF, MFF, PART_VEC_F_AO, PART_VEC_F_AO,          &
                                 NUM2, NUM2, MOO_ROW_MAX_TERMS, 'MOO', NTERM_MOO, NDOFO, SYM_MOO, I_MOO, J_MOO, MOO )
         ENDIF

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Reduce MFF to MAA = MAA(bar) + MAO*GOA + (MAO*GOA)' + GOA'*MOO*GOA.
! If PARAM MATSPARS = 'Y', then we use sparse matrix operations (multiply/add/transpose). If not, use full matrix operations

      IF (MATSPARS == 'Y') THEN

         CALL ALLOCATE_SPARSE_MAT ( 'GOAt', NDOFA, NTERM_GOA, SUBR_NAME )
         CALL MATTRNSP_SS ( NDOFO, NDOFA, NTERM_GOA, 'GOA', I_GOA, J_GOA, GOA, 'GOAt', I_GOAt, J_GOAt, GOAt )

                                                           ! CCS1 will be sparse CCS format version of sparse CRS matrix GOA
         CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFA, NTERM_GOA, SUBR_NAME )
         CALL SPARSE_CRS_SPARSE_CCS ( NDOFO, NDOFA, NTERM_GOA, 'GOA', I_GOA, J_GOA, GOA, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

         IF (NTERM_MAO > 0) THEN                           ! Part I of reduced MAA: calc MAO*GOA & add it & transpose to orig MAA

                                                           ! I-1 , sparse multiply to get CRS1 = MAO*GOA. Use CCS1 for GOA CCS
            CALL MATMULT_SSS_NTERM ( 'MAO' , NDOFA, NTERM_MAO , SYM_MAO, I_MAO , J_MAO ,                                           &
                                     'GOA' , NDOFA, NTERM_GOA , SYM_GOA, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1',        NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFA, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'MAO' , NDOFA, NTERM_MAO , SYM_MAO, I_MAO , J_MAO , MAO ,                                           &
                               'GOA' , NDOFA, NTERM_GOA , SYM_GOA, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            NTERM_CRS2 = NTERM_CRS1                        ! I-2 , allocate memory to array CRS2 which will hold transpose of CRS1
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFA, NTERM_CRS2, SUBR_NAME )

                                                           ! I-3 , transpose CRS1 to get CRS2 = (MAO*GOA)t
            CALL MATTRNSP_SS ( NDOFA, NDOFA, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CRS2', I_CRS2, J_CRS2, CRS2 )

                                                           ! I-4 , sparse add to get CRS3 = CRS1 + CRS2 = (MAO*GOA) + (MAO*GOA)t
            CALL MATADD_SSS_NTERM (NDOFA,'MAO*GOA', NTERM_CRS1, I_CRS1, J_CRS1, 'N', '(MAO*GOA)t', NTERM_CRS2, I_CRS2, J_CRS2, 'N',&
                                         'CRS3', NTERM_CRS3)
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFA, NTERM_CRS3, SUBR_NAME )
            CALL MATADD_SSS ( NDOFA, 'MAO*GOA', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE, '(MAO*GOA)t', NTERM_CRS2,                   &
                              I_CRS2, J_CRS2, CRS2, ONE, 'CRS3', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! I-5 , deallocate CRS1, CRS2
            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )

                                                           ! I-6 , CRS3 = (MAO*GOA) + (MAO*GOA)t has all nonzero terms in it.
            IF      (SPARSTOR == 'SYM   ') THEN            !      If SPARSTOR == 'SYM   ', rewrite CRS3 as sym in CRS1

               CALL SPARSE_CRS_TERM_COUNT ( NDOFA, NTERM_CRS3, '(MAO*GOA) + (MAO*GOA)t all nonzeros', I_CRS3, J_CRS3, NTERM_CRS1 )
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFA, NTERM_CRS1, SUBR_NAME )
               CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS3 = (MAO*GOA) + (MAO*GOA)t all nonzeros', NDOFA, NTERM_CRS3, I_CRS3, J_CRS3, CRS3, &
                                            'CRS1 = (MAO*GOA) + (MAO*GOA)t stored sym'  ,        NTERM_CRS1, I_CRS1, J_CRS1, CRS1 )
               SYM_CRS1 = 'Y'

            ELSE IF (SPARSTOR == 'NONSYM') THEN            !      If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

               NTERM_CRS1 = NTERM_CRS3
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFA, NTERM_CRS1, SUBR_NAME )
               DO I=1,NDOFA+1
                  I_CRS1(I) = I_CRS3(I)
               ENDDO
               DO I=1,NTERM_CRS1
                  J_CRS1(I) = J_CRS3(I)
                    CRS1(I) =   CRS3(I)
               ENDDO
               SYM_CRS1 = 'N'

            ELSE                                           !      Error - incorrect SPARSTOR

               WRITE(ERR,932) SUBR_NAME, SPARSTOR
               WRITE(F06,932) SUBR_NAME, SPARSTOR
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF

            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! I-7 , deallocate CRS3

                                                           ! I-8 , sparse add to get CRS3 = MAA-bar + CRS1
            CALL MATADD_SSS_NTERM ( NDOFA, 'MAA-bar', NTERM_MAA , I_MAA , J_MAA , SYM_MAA , '(MAO*GOA) + (MAO*GOA)t',              &
                                                      NTERM_CRS1, I_CRS1, J_CRS1, SYM_CRS1, 'CRS3', NTERM_CRS3 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFA, NTERM_CRS3, SUBR_NAME )
            CALL MATADD_SSS ( NDOFA, 'MAA-bar', NTERM_MAA, I_MAA, J_MAA, MAA, ONE, '(MAO*GOA) + (MAO*GOA)t', NTERM_CRS1,           &
                              I_CRS1, J_CRS1, CRS1, ONE, 'CRS1', NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! I-9 , deallocate CRS1 and CRS3

            NTERM_MAA = NTERM_CRS3                         ! I-10, reallocate MAA to be size of CRS1
            WRITE(SC1, * ) '    Reallocate MAA'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MAA', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'MAA' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   MAA', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'MAA', NDOFA, NTERM_MAA, SUBR_NAME )
                                                           ! I-11, set MAA = CRS1
            DO I=1,NDOFA+1
               I_MAA(I) = I_CRS3(I)
            ENDDO
            DO J=1,NTERM_MAA
               J_MAA(J) = J_CRS3(J)
                 MAA(J) =   CRS3(J)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! I-12, deallocate CRS3

         ENDIF

         IF (NTERM_MOO > 0) THEN                           ! Part II of reduced MAA: calc GOA(t)*MOO*GOA and add to MAA

                                                           ! II-1 , sparse multiply to get CRS1 = MOO*GOA using CCS1 for GOA CCS
            CALL MATMULT_SSS_NTERM ( 'MOO' , NDOFO, NTERM_MOO , SYM_MOO, I_MOO , J_MOO ,                                           &
                                     'GOA' , NDOFA, NTERM_GOA , SYM_GOA, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1',        NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFO, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'MOO' , NDOFO, NTERM_MOO , SYM_MOO, I_MOO , J_MOO , MOO ,                                           &
                               'GOA' , NDOFA, NTERM_GOA , SYM_GOA, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )             ! II-2 , deallocate CCS1

            NTERM_CCS1 = NTERM_CRS1                        ! II-3 , allocate CCS1 to be same as CRS1 but in CCS format
            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFA, NTERM_CCS1, SUBR_NAME )
            CALL SPARSE_CRS_SPARSE_CCS ( NDOFO, NDOFA, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! II-4 , deallocate CRS1
                                                           ! II-5 , sparse multiply to get CRS1 = GOAt*CCS1 with CCS1= MOO*GOA
!                                                            (note: use SYM_GOA for sym indicator of CCS1)
            CALL MATMULT_SSS_NTERM ( 'GOAt', NDOFA, NTERM_GOA , SYM_GOA, I_GOAt, J_GOAt,                                           &
                                     'CCS1', NDOFA, NTERM_CCS1, SYM_GOA, J_CCS1, I_CCS1, AROW_MAX_TERMS,                           &
                                     'CRS1 ',       NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFA, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'GOAt', NDOFA, NTERM_GOA , SYM_GOA, I_GOAt, J_GOAt, GOAt,                                           &
                               'CCS1', NDOFA, NTERM_CCS1, SYM_GOA, J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )             ! II-6 , deallocate CCS1

                                                           ! II-7 , CRS1 = GOAt*MOO*GOA has all nonzero terms in it.
            IF      (SPARSTOR == 'SYM   ') THEN            !      If SPARSTOR == 'SYM   ', rewrite CRS1 as sym in CRS3

               CALL SPARSE_CRS_TERM_COUNT ( NDOFA, NTERM_CRS1, 'GOAt*MOO*GOA all nonzeros', I_CRS1, J_CRS1, NTERM_CRS3 )
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFA, NTERM_CRS3, SUBR_NAME )
!xx            CALL SPARSE_NONSYM_TO_SYM ( NDOFA, NTERM_CRS1, 'CRS1', I_CRS1, J_CRS1, CRS1,                                        &
!xx                                               NTERM_CRS3, 'CRS3', I_CRS3, J_CRS3, CRS3 )
               CALL CRS_NONSYM_TO_CRS_SYM ( 'CRS1 = GOAt*MOO*GOA all nonzeros', NDOFA, NTERM_CRS1, I_CRS1, J_CRS1, CRS1,           &
                                            'CRS3 = GOAt*MOO*GOA stored sym'  ,        NTERM_CRS3, I_CRS3, J_CRS3, CRS3 )
               SYM_CRS3 = 'Y'

            ELSE IF (SPARSTOR == 'NONSYM') THEN            !      If SPARSTOR == 'NONSYM', rewrite CRS3 in CRS1 with NTERM_CRS3

               NTERM_CRS3 = NTERM_CRS1
               CALL ALLOCATE_SCR_CRS_MAT ( 'CRS3', NDOFA, NTERM_CRS3, SUBR_NAME )
               DO I=1,NDOFA+1
                  I_CRS3(I) = I_CRS1(I)
               ENDDO
               DO I=1,NTERM_CRS3
                  J_CRS3(I) = J_CRS1(I)
                    CRS3(I) =   CRS1(I)
               ENDDO
               SYM_CRS3 = 'N'

            ELSE                                           !      Error - incorrect SPARSTOR

               WRITE(ERR,932) SUBR_NAME, SPARSTOR
               WRITE(F06,932) SUBR_NAME, SPARSTOR
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )

            ENDIF
                                                           ! II-8 , sparse add to get CRS2 = MAA + CRS3 = MAA + GOAt*MOO*GOA
            CALL MATADD_SSS_NTERM ( NDOFA, 'MAA-bar + (MAO*GOA) + (MAO*GOA)t', NTERM_MAA, I_MAA, J_MAA, SYM_MAA, 'GOAt*MOO*GOA',   &
                                    NTERM_CRS3, I_CRS3, J_CRS3, SYM_CRS3, 'CRS2', NTERM_CRS2 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFA, NTERM_CRS2, SUBR_NAME )
            CALL MATADD_SSS ( NDOFA, 'MAA-bar + (MAO*GOA) + (MAO*GOA)t' , NTERM_MAA , I_MAA , J_MAA , MAA, ONE, 'GOAt*MOO*GOA',    &
                              NTERM_CRS3, I_CRS3, J_CRS3, CRS3, ONE, 'CRS2', NTERM_CRS2, I_CRS2, J_CRS2, CRS2 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! II-9 , deallocate CRS1 and CRS3
            CALL DEALLOCATE_SCR_MAT ( 'CRS3' )             ! II-10, deallocate CRS3

            NTERM_MAA = NTERM_CRS2                         ! II-11, reallocate MAA to be size of CRS2
            WRITE(SC1, * ) '    Reallocate MAA'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MAA', CR13
             CALL DEALLOCATE_SPARSE_MAT ( 'MAA' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   MAA', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'MAA', NDOFA, NTERM_MAA, SUBR_NAME )

            DO I=1,NDOFA+1                                 ! II-12, set reduced MAA = CRS2 until we see if MAO is null, below
               I_MAA(I) = I_CRS2(I)
            ENDDO
            DO I=1,NTERM_MAA
               J_MAA(I) = J_CRS2(I)
                 MAA(I) =   CRS2(I)
            ENDDO

            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! II-13, deallocate CRS2

         ELSE

            CALL DEALLOCATE_SCR_MAT ( 'CCS1')

         ENDIF

         CALL DEALLOCATE_SPARSE_MAT ( 'GOAt' )

! ----------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (MATSPARS == 'N') THEN

         CALL ALLOCATE_FULL_MAT ( 'MAA_FULL', NDOFA, NDOFA, SUBR_NAME )

         IF (NTERM_MAA > 0) THEN
            CALL SPARSE_CRS_TO_FULL ( 'MAA', NTERM_MAA, NDOFA, NDOFA, SYM_MAA, I_MAA, J_MAA,MAA, MAA_FULL )
         ENDIF

         IF (NTERM_MAO > 0) THEN                           ! Part 1: calc MAO*GOA and add it & it's transpose to MAA_FULL

            CALL ALLOCATE_FULL_MAT ( 'MAO_FULL', NDOFA, NDOFO, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'MAO', NTERM_MAO, NDOFA, NDOFO, SYM_MAO, I_MAO, J_MAO, MAO, MAO_FULL )

            CALL ALLOCATE_FULL_MAT ( 'GOA_FULL', NDOFO, NDOFA, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'GOA', NTERM_GOA, NDOFO, NDOFA, SYM_GOA, I_GOA, J_GOA, GOA, GOA_FULL )

            CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFA, NDOFA, SUBR_NAME )

            CALL MATMULT_FFF (  MAO_FULL, GOA_FULL, NDOFA, NDOFO, NDOFA, DUM1 )

            CALL DEALLOCATE_FULL_MAT ( 'MAO_FULL' )
            CALL DEALLOCATE_FULL_MAT ( 'GOA_FULL' )
            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFA, NDOFA, SUBR_NAME )

            ITRNSPB = 0                                    ! Calc MAA-bar + MAO*GOA = MAA + DUM1 and put into DUM1
            CALL MATADD_FFF ( MAA_FULL, DUM1, NDOFA, NDOFA, ALPHA, BETA, ITRNSPB, DUM2 )
            ITRNSPB = 1                                    ! Calc MAA-bar + MAO*GOA + (OAM*GOA)' = DUM2+DUM1' and put into MAA_FULL
            CALL MATADD_FFF ( DUM2    , DUM1, NDOFA, NDOFA, ALPHA, BETA, ITRNSPB, MAA_FULL )


            CALL DEALLOCATE_FULL_MAT ( 'DUM1' )
            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

         ENDIF

         IF (NTERM_MOO > 0) THEN                           ! Part 2: calc GOA(t)*MOO*GOA and add to MAA_FULL

            CALL ALLOCATE_FULL_MAT ( 'MOO_FULL', NDOFO, NDOFO, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'MOO', NTERM_MOO, NDOFO, NDOFO, SYM_MOO, I_MOO, J_MOO, MOO, MOO_FULL )

            CALL ALLOCATE_FULL_MAT ( 'GOA_FULL', NDOFO, NDOFA, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'GOA', NTERM_GOA, NDOFO, NDOFA, SYM_GOA, I_GOA, J_GOA, GOA, GOA_FULL )

            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFO, NDOFA, SUBR_NAME )

            CALL MATMULT_FFF ( MOO_FULL, GOA_FULL, NDOFO, NDOFO, NDOFA, DUM2 )

            CALL DEALLOCATE_FULL_MAT ( 'MOO_FULL' )

            CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFA, NDOFA, SUBR_NAME )

            CALL MATMULT_FFF_T (GOA_FULL, DUM2, NDOFO, NDOFA, NDOFA, DUM1 )

            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

            CALL DEALLOCATE_FULL_MAT ( 'GOA_FULL' )

            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFA, NDOFA, SUBR_NAME )

            ITRNSPB = 0
            CALL MATADD_FFF ( MAA_FULL, DUM1, NDOFA, NDOFA, ALPHA, BETA, ITRNSPB, DUM2 )
            DO I=1,NDOFA
               DO J=1,NDOFA
                  MAA_FULL(I,J) = DUM2(I,J)
               ENDDO
            ENDDO

            CALL DEALLOCATE_FULL_MAT ( 'DUM1' )
            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

         ENDIF

         CALL CNT_NONZ_IN_FULL_MAT ( 'MAA_FULL  ', MAA_FULL, NDOFA, NDOFA, SYM_MAA, NTERM_MAA, SMALL )

         CALL DEALLOCATE_SPARSE_MAT ( 'MAA' )
         CALL ALLOCATE_SPARSE_MAT ( 'MAA', NDOFA, NTERM_MAA, SUBR_NAME )

         IF (NTERM_MAA > 0) THEN                            ! Create new sparse arrays from MAA_FULL
!xx         CALL FULL_TO_SPARSE_CRS ( 'MAA_FULL  ', NDOFA, NDOFA, MAA_FULL, NTERM_MAA, SYM_MAA, I_MAA, J_MAA, MAA )
            CALL DEALLOCATE_FULL_MAT ( 'MAA_FULL' )
         ENDIF

      ELSE

         WRITE(ERR,911) SUBR_NAME,MATSPARS
         WRITE(F06,911) SUBR_NAME,MATSPARS
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
  911 FORMAT(' *ERROR   911: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER MATSPARS MUST BE EITHER ','Y',' OR ','N',' BUT VALUE IS ',A)

  932 FORMAT(' *ERROR   932: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER SPARSTOR MUST BE EITHER "SYM" OR "NONSYM" BUT VALUE IS ',A)

12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_MFF_TO_MAA


      SUBROUTINE REDUCE_PF_TO_PA ( PART_VEC_F_AO, PART_VEC_SUB )

! Call routines to reduce the PF grid point load matrix from the F-set to the A, O-sets. See Appendix B to the MYSTRAN User's
! Reference Manual for the derivation of the reduction equations.

! NOTE: This subr has code for sparse matrices as well as full matrices, (i.e. Bulk Data PARAM MATSPARS = 'Y' for sparse and 'N'
! for full). The code for full matrices was put in originally in order that the sparse code could be thoroughly checked. That task
! is complete and the remaining full matrix code has not been maintained since around 2005. In addition, new capability added to
! MYSTRAN since that approx time does not have full matrix code.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, KOO_SDIA, NDOFF, NDOFA, NDOFO, NSUB, NTERM_GOA, NTERM_PF,        &
                                         NTERM_PA, NTERM_PO
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE
      USE PARAMS, ONLY                :  EPSIL, MATSPARS
      USE SPARSE_MATRICES, ONLY       :  I_PF, J_PF, PF, I_PA, J_PA, PA, I_PO, J_PO, PO, I_GOA, J_GOA, GOA, I_GOAt, J_GOAt, GOAt
      USE SPARSE_MATRICES, ONLY       :  SYM_GOA, SYM_PF, SYM_PA, SYM_PO
      USE FULL_MATRICES, ONLY         :  PA_FULL, PO_FULL, GOA_FULL, DUM1, DUM2
      USE SCRATCH_MATRICES

      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM, MATMULT_SSS, MATMULT_SSS_NTERM, MATTRNSP_SS
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CCS_MAT, ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE SPARSE_FORMAT_CONVERSION, ONLY:  SPARSE_CRS_SPARSE_CCS, SPARSE_CRS_TO_FULL
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE FULL_MATRIX_LIFECYCLE, ONLY :  ALLOCATE_FULL_MAT, DEALLOCATE_FULL_MAT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATADD_FFF, MATMULT_FFF_T
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CNT_NONZ_IN_FULL_MAT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'REDUCE_PF_TO_PA'

      INTEGER(LONG), INTENT(IN)       :: PART_VEC_F_AO(NDOFF)! Partitioning vector (F set into A and O sets)
      INTEGER(LONG), INTENT(IN)       :: PART_VEC_SUB(NSUB)  ! Partitioning vector (1's for all subcases)
      INTEGER(LONG)                   :: AROW_MAX_TERMS      ! Output from MATMULT_SFS_NTERM and input to MATMULT_SFS
      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG), PARAMETER        :: ITRNSPB     = 0     ! Transpose indicator for matrix multiply routine
      INTEGER(LONG)                   :: NTERM_CRS1          ! Number of terms in matrix CRS1
      INTEGER(LONG)                   :: NTERM_CRS2          ! Number of terms in matrix CRS2
      INTEGER(LONG), PARAMETER        :: NUM1        = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2        = 2     ! Used in subr's that partition matrices
      INTEGER(LONG)                   :: PA_ROW_MAX_TERMS    ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: PO_ROW_MAX_TERMS    ! Output from subr PARTITION_SIZE (max terms in any row of matrix)


      REAL(DOUBLE)                    :: ALPHA = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: BETA  = ONE         ! Scalar multiplier for matrix
      REAL(DOUBLE)                    :: SMALL               ! A number used in filtering out small numbers from a full matrix

      INTRINSIC                       :: DABS




! **********************************************************************************************************************************
! Partition PA from PF (This is PA before reduction, or PA(bar) )

      IF (NDOFA > 0) THEN

         CALL PARTITION_SS_NTERM ( 'PF' , NTERM_PF , NDOFF, NSUB , SYM_PF , I_PF , J_PF      , PART_VEC_F_AO, PART_VEC_SUB,        &
                                    NUM1, NUM1, PA_ROW_MAX_TERMS, 'PA', NTERM_PA, SYM_PA )

         CALL ALLOCATE_SPARSE_MAT ( 'PA', NDOFA, NTERM_PA, SUBR_NAME )

         IF (NTERM_PA  > 0) THEN
            CALL PARTITION_SS ( 'PF' , NTERM_PF , NDOFF, NSUB , SYM_PF , I_PF , J_PF , PF , PART_VEC_F_AO, PART_VEC_SUB,           &
                                 NUM1, NUM1, PA_ROW_MAX_TERMS, 'PA', NTERM_PA , NDOFA, SYM_PA, I_PA , J_PA , PA  )
         ENDIF

      ENDIF

! Partition PO from PF

      IF (NDOFO > 0) THEN

         CALL PARTITION_SS_NTERM ( 'PF' , NTERM_PF , NDOFF, NSUB , SYM_PF , I_PF , J_PF ,      PART_VEC_F_AO, PART_VEC_SUB,        &
                                    NUM2, NUM1, PO_ROW_MAX_TERMS, 'PO', NTERM_PO, SYM_PO )

         CALL ALLOCATE_SPARSE_MAT ( 'PO', NDOFO, NTERM_PO, SUBR_NAME )

         IF (NTERM_PO  > 0) THEN
            CALL PARTITION_SS ( 'PF' , NTERM_PF , NDOFF, NSUB , SYM_PF , I_PF , J_PF , PF , PART_VEC_F_AO, PART_VEC_SUB,           &
                                 NUM2, NUM1, PO_ROW_MAX_TERMS, 'PO', NTERM_PO , NDOFO, SYM_PO, I_PO , J_PO , PO  )
         ENDIF

      ENDIF

! Reduce PF to PA = PA(bar) + GOA'*PO
! If PARAM MATSPARS = 'Y', then we use sparse matrix operations (multiply/add/transpose). If not, use full matrix operations

      IF (MATSPARS == 'Y') THEN


         IF (NTERM_PO > 0) THEN                            ! Calc GOAt*PO & add it orig PA

            CALL ALLOCATE_SPARSE_MAT ( 'GOAt', NDOFA, NTERM_GOA, SUBR_NAME )
            CALL MATTRNSP_SS ( NDOFO, NDOFA, NTERM_GOA, 'GOA', I_GOA, J_GOA, GOA, 'GOAt', I_GOAt, J_GOAt, GOAt )
                                                           ! CCS1 will be sparse CCS format version of sparse CRS matrix PO
!! *********************************************************************************************************************************
! ERROR: NDOFO should be NSUB for col size of PO
!!          CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NDOFO, NTERM_PO, SUBR_NAME )
!! *********************************************************************************************************************************
            CALL ALLOCATE_SCR_CCS_MAT ( 'CCS1', NSUB, NTERM_PO, SUBR_NAME )
            CALL SPARSE_CRS_SPARSE_CCS ( NDOFO, NSUB, NTERM_PO, 'PO', I_PO, J_PO, PO, 'CCS1', J_CCS1, I_CCS1, CCS1, 'Y' )

                                                           ! Sparse multiply to get CRS1 = GOAt*PO = GOAt*CCS1.
! *********************************************************************************************************************************
! ERROR: NDOFO should be NSUB for col size of PO
!!          CALL MATMULT_SSS_NTERM ( 'GOAt', NDOFA, NTERM_GOA, SYM_GOA, I_GOAt, J_GOAt,                                            &
!!                                   'PO'  , NDOFO, NTERM_PO , SYM_PO , J_CCS1, I_CCS1, AROW_MAX_TERMS,                            &
!!                                   'CRS1',        NTERM_CRS1 )
! *********************************************************************************************************************************
            CALL MATMULT_SSS_NTERM ( 'GOAt', NDOFA, NTERM_GOA, SYM_GOA, I_GOAt, J_GOAt,                                            &
                                     'PO'  , NSUB , NTERM_PO , SYM_PO , J_CCS1, I_CCS1, AROW_MAX_TERMS,                            &
                                     'CRS1',        NTERM_CRS1 )

            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFA, NTERM_CRS1, SUBR_NAME )

            CALL MATMULT_SSS ( 'GOAt', NDOFA, NTERM_GOA , SYM_GOA, I_GOAt, J_GOAt, GOAt,                                           &
                               'PO'  , NSUB , NTERM_PO  , SYM_PO , J_CCS1, I_CCS1, CCS1, AROW_MAX_TERMS,                           &
                               'CRS1', ONE  , NTERM_CRS1,          I_CRS1, J_CRS1, CRS1 )

            CALL DEALLOCATE_SCR_MAT ( 'CCS1' )

                                                           ! Sparse add to get CRS2 = PA-bar + CRS1
            CALL MATADD_SSS_NTERM ( NDOFA, 'PA-bar', NTERM_PA, I_PA, J_PA, SYM_PA, 'GOAt*PO', NTERM_CRS1, I_CRS1, J_CRS1, 'N',     &
                                         'CRS2', NTERM_CRS2 )
            CALL ALLOCATE_SCR_CRS_MAT ( 'CRS2', NDOFA, NTERM_CRS2, SUBR_NAME )
            CALL MATADD_SSS ( NDOFA, 'PA-bar', NTERM_PA, I_PA, J_PA, PA, ONE, 'GOAt*PO', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, ONE,    &
                                     'CRS2', NTERM_CRS2, I_CRS2, J_CRS2, CRS2 )

            CALL DEALLOCATE_SCR_MAT ( 'CRS1' )             ! Deallocate CRS1

            NTERM_PA = NTERM_CRS2                          ! Reallocate KAA to be size of CRS2
            WRITE(SC1, * ) '    Reallocate PA '
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PA ', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'PA' )
            WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PA ', CR13
            CALL ALLOCATE_SPARSE_MAT ( 'PA', NDOFA, NTERM_PA, SUBR_NAME )
                                                           ! Set KAA = CRS1
            DO I=1,NDOFA+1
               I_PA(I) = I_CRS2(I)
            ENDDO
            DO J=1,NTERM_PA
               J_PA(J) = J_CRS2(J)
                 PA(J) =   CRS2(J)
            ENDDO

            WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
            WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GOAt', CR13
            CALL DEALLOCATE_SPARSE_MAT ( 'GOAt' )
            CALL DEALLOCATE_SCR_MAT ( 'CRS2' )             ! Deallocate CRS2 and GOAt

         ENDIF

      ELSE IF (MATSPARS == 'N') THEN

         CALL ALLOCATE_FULL_MAT ( 'PA_FULL', NDOFA, NSUB, SUBR_NAME )  ! Calc reduced PA

         IF (NTERM_PA > 0) THEN
            CALL SPARSE_CRS_TO_FULL ( 'PA', NTERM_PA, NDOFA, NSUB, SYM_PA, I_PA, J_PA, PA, PA_FULL )
         ENDIF

         IF (NTERM_PO > 0) THEN                            ! Calc GOA(t)*PO and add to partitioned PA to get reduced PA

            CALL ALLOCATE_FULL_MAT ( 'PO_FULL', NDOFO, NSUB, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'PO', NTERM_PO, NDOFO, NSUB, SYM_PO, I_PO, J_PO, PO, PO_FULL )

            CALL ALLOCATE_FULL_MAT ( 'GOA_FULL', NDOFO, NDOFA, SUBR_NAME )
            CALL SPARSE_CRS_TO_FULL ( 'GOA', NTERM_GOA, NDOFO, NDOFA, SYM_GOA, I_GOA, J_GOA, GOA, GOA_FULL )

            CALL ALLOCATE_FULL_MAT ( 'DUM1', NDOFA, NSUB, SUBR_NAME )

            CALL MATMULT_FFF_T (GOA_FULL, PO_FULL, NDOFO, NDOFA, NSUB, DUM1 )

            CALL DEALLOCATE_FULL_MAT ( 'PO_FULL' )
            CALL DEALLOCATE_FULL_MAT ( 'GOA_FULL' )

            CALL ALLOCATE_FULL_MAT ( 'DUM2', NDOFA, NSUB, SUBR_NAME )
                                                           ! Final reduced PA
            CALL MATADD_FFF ( PA_FULL, DUM1, NDOFA, NSUB, ALPHA, BETA, ITRNSPB, DUM2 )
            DO I=1,NDOFA
               DO J=1,NSUB
                  PA_FULL(I,J) = DUM2(I,J)
               ENDDO
            ENDDO

            CALL DEALLOCATE_FULL_MAT ( 'DUM1' )
            CALL DEALLOCATE_FULL_MAT ( 'DUM2' )

         ENDIF

         CALL CNT_NONZ_IN_FULL_MAT ( 'PA_FULL  ', PA_FULL, NDOFA, NSUB, SYM_PA, NTERM_PA, SMALL )

         WRITE(SC1, * ) '    Reallocate PA '
   !xx   WRITE(SC1, * )                                    ! Advance 1 line for screen messages
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PA ', CR13  ;  CALL DEALLOCATE_SPARSE_MAT ( 'PA' )
         WRITE(SC1,12345,ADVANCE='NO') '       Allocate   PA ', CR13  ;  CALL ALLOCATE_SPARSE_MAT ('PA', NDOFA, NTERM_PA, SUBR_NAME)

         IF (NTERM_PA > 0) THEN                            ! Create new sparse arrays from PA_FULL
!xx         CALL FULL_TO_SPARSE_CRS ( 'PA_FULL   ', NDOFA, NSUB , PA_FULL,  NTERM_PA,  SYM_PA,  I_PA,  J_PA,  PA  )
            CALL DEALLOCATE_FULL_MAT ( 'PA_FULL' )
         ENDIF

      ELSE

         WRITE(ERR,911) SUBR_NAME,MATSPARS
         WRITE(F06,911) SUBR_NAME,MATSPARS
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ENDIF

! If there are applied loads on the O-set, solve for UO0 = KOO(-1) x PO which is part of the final solution for UO

      IF (NTERM_PO > 0) THEN
         CALL SOLVE_UO0
      ENDIF



      RETURN

! **********************************************************************************************************************************
  911 FORMAT(' *ERROR   911: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' PARAMETER MATSPARS MUST BE EITHER ','Y',' OR ','N',' BUT VALUE IS ',A)


12345 FORMAT(A,10X,A)

! **********************************************************************************************************************************

      END SUBROUTINE REDUCE_PF_TO_PA




      SUBROUTINE SOLVE_GOA

! Solves the sustem of equations: KOO*GOA = -KAO' for matrix GOA which is used in the reduction of the F set stiffness, mass and
! load matrices from the F-set to the A, O_sets

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, ERR, F06, SCR, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FACTORED_MATRIX, FATAL_ERR, KOO_SDIA, NDOFA, NDOFO, NTERM_GOA, NTERM_KOO,   &
                                         NTERM_KAO
      USE PARAMS, ONLY                :  EPSIL, PRTGOA
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  SOLLIB, SPARSE_FLAVOR
      USE SPARSE_MATRICES, ONLY       :  I2_GOA, I_GOA, J_GOA, GOA, I_KOO, J_KOO, KOO, I_KAO, J_KAO, KAO

! Interface module not needed for subr's DPBTRF and DPBTRS. These are "CONTAIN'ed" in module LAPACK_LIN_EQN_DPB, which
! is "USE'd" above

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE
      USE SPARSE_CRS_ACCESS, ONLY     :  GET_SPARSE_CRS_ROW
      USE LAPACK_ADAPTERS, ONLY       :  FBS_LAPACK
      USE SUPERLU_ADAPTERS, ONLY      :  FBS_SUPRLU
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  READ_MATRIX_1, READ_MATRIX_2, WRITE_SPARSE_CRS
      USE SORTING, ONLY               :  SORT_INT2_REAL1
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_GOA'
      CHARACTER(  1*BYTE)             :: CLOSE_IT          ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER(  8*BYTE)             :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER( 24*BYTE)             :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN
      CHARACTER( 22*BYTE)             :: MODNAM1           ! Name to write to screen to describe module being run
      CHARACTER(  1*BYTE)             :: READ_NTERM        ! 'Y' or 'N' Input to subr READ_MATRIX_1
      CHARACTER(  1*BYTE)             :: NULL_COL          ! 'Y' if a col of KAO(transpose) is null
      CHARACTER(  1*BYTE)             :: OPND              ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SCRFIL            ! File name

      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: INFO              ! Info on success of factorization or solve
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: GOA_COL(NDOFO)    ! A column of GOA solved for herein
      REAL(DOUBLE)                    :: INOUT_COL(NDOFO)  ! A column of KAO'
      REAL(DOUBLE)                    :: KOO_SCALE_FACS(NDOFO)
                                                           ! KOO scale facs. KOO will not be equilibrated so these are set to 1.0

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Make units for writing errors the screen and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Make sure that ABAND (KOO triangular factor) was successfully generated when SYM_MAT_DECOMP_LAPACK was run.

      IF (SOLLIB == 'BANDED  ') THEN
         IF (FACTORED_MATRIX(1:3) /= 'KOO') THEN
            WRITE(ERR,2504) SUBR_NAME, FACTORED_MATRIX
            WRITE(F06,2504) SUBR_NAME, FACTORED_MATRIX
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDIF

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Solve for GOA

! Open a scratch file that will be used to write GOA nonzero terms to as we solve for columns of GOA. After all col's
! of GOA have been solved for, and we have a count on NTERM_GOA, we will allocate memory to the GOA arrays and read
! the scratch file values into those arrays. Then, in the calling subroutine, we will write NTERM_GOA, followed by
! GOA  row/col/value to a permanent file

      SCRFIL(1:)  = ' '
      SCRFIL(1:9) = 'SCRATCH-991'
      OPEN (SCR(1),STATUS='SCRATCH',FORM='UNFORMATTED',ACTION='READWRITE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SCRFIL, OUNT )
         CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )
         CALL OUTA_HERE ( 'Y' )                            ! Can't open scratch file, so quit
      ENDIF

! Loop on columns of KAO(transpose)

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      NTERM_GOA = 0
      DO J=1,NDOFA

         DO I=1,NDOFO
            INOUT_COL(I)       = ZERO                      ! Initialize INOUT_COL since GET_SPARSE_CRS_COL won't
            KOO_SCALE_FACS(I)  = ONE                       ! Initialize these
         ENDDO
         NULL_COL = 'Y'
         CALL GET_SPARSE_CRS_ROW ( 'KAO',J, NTERM_KAO, NDOFA, NDOFO, I_KAO, J_KAO, KAO, -ONE, INOUT_COL, NULL_COL )

         IF (NULL_COL == 'N') THEN                         ! FBS will solve for GOA_COL & load it into GOA array

            CALL OURTIM
            MODNAM1 = '    Solve for GOA col '
            WRITE(SC1,22345) MODNAM1, J, NDOFA
                                                           ! FBS should not equilibrate since KOO was prevented from equilibrating
            IF      (SOLLIB == 'BANDED  ') THEN

               CALL FBS_LAPACK ( 'N', NDOFO, KOO_SDIA, KOO_SCALE_FACS, INOUT_COL )

            ELSE IF (SOLLIB == 'SPARSE  ') THEN

               IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

                  INFO = 0
                  CALL FBS_SUPRLU ( SUBR_NAME, 'KOO', NDOFO, NTERM_KOO, I_KOO, J_KOO, KOO, J, INOUT_COL, INFO )

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

            DO I=1,NDOFO
               GOA_COL(I) = INOUT_COL(I)
            ENDDO

            DO I=1,NDOFO                                   ! Count NTERM_GOA and write nonzero GOA to scratch file
               IF (DABS(GOA_COL(I)) > EPS1) THEN
                  NTERM_GOA = NTERM_GOA + 1
                  WRITE(SCR(1)) I,J,GOA_COL(I)
               ENDIF
            ENDDO

         ELSE

            DO I=1,NDOFO
               GOA_COL(I) = ZERO
            ENDDO

         ENDIF

      ENDDO

! The GOA data in SCRATCH-991 is written one col at a time. We need it to be  written for one row at a time with all
! rows in numerical order.

      CALL ALLOCATE_SPARSE_MAT ( 'GOA', NDOFO, NTERM_GOA, SUBR_NAME )
      CALL ALLOCATE_L2_GOA_2 ( SUBR_NAME )
      REWIND (SCR(1))
      MESSAG = 'SCRATCH: GOA ROW/COL/VAL'
      READ_NTERM = 'N'
      OPND       = 'Y'
      CLOSE_IT   = 'N'
      CLOSE_STAT = 'KEEP    '
      CALL READ_MATRIX_2 (SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG,'GOA',NDOFO, NTERM_GOA, READ_NTERM, I2_GOA, J_GOA, GOA)
      CALL SORT_INT2_REAL1 (SUBR_NAME, 'I2_GOA, J_GOA, GOA', NTERM_GOA, I2_GOA, J_GOA, GOA )
      REWIND (SCR(1))
      WRITE(SCR(1)) NTERM_GOA
      DO K=1,NTERM_GOA
         WRITE(SCR(1)) I2_GOA(K),J_GOA(K),GOA(K)
      ENDDO

! Reallocate memory to GOA based on the NTERM_GOA counted above and read values from scratch file into GOA arrays

      WRITE(SC1, * ) '    Reallocate GOA'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GOA'     , CR13   ;   CALL DEALLOCATE_SPARSE_MAT ( 'GOA' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate L2_GOA_2', CR13   ;   CALL DEALLOCATE_L2_GOA_2
      WRITE(SC1,12345,ADVANCE='NO') '       Allocate   GOA'     , CR13
      CALL ALLOCATE_SPARSE_MAT ('GOA', NDOFO, NTERM_GOA, SUBR_NAME)

      REWIND (SCR(1))
      SCRFIL(1:)  = ' '
      SCRFIL(1:9) = 'SCRATCH-991'
      MESSAG = 'SCRATCH: GOA ROW/COL/VAL'
      READ_NTERM = 'Y'
      OPND       = 'Y'
      CLOSE_IT   = 'Y'
      CLOSE_STAT = 'DELETE  '
      CALL READ_MATRIX_1 ( SCRFIL, SCR(1), OPND, CLOSE_IT, CLOSE_STAT, MESSAG, 'GOA', NTERM_GOA, READ_NTERM, NDOFO,                &
                           I_GOA, J_GOA, GOA)

      CALL FILE_CLOSE ( SCR(1), SCRFIL, 'DELETE' )

! Print out constraint matrix GOA, if requested

      IF ( PRTGOA == 1) THEN
         IF (NTERM_GOA > 0) THEN
            CALL WRITE_SPARSE_CRS ( 'CONSTRAINT MATRIX GOA', 'O ', 'A ', NTERM_GOA, NDOFO, I_GOA, J_GOA, GOA )
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 2504 FORMAT(' *ERROR  2504: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NAME OF MATRIX THAT HAS BEEN DECOMPOSED INTO TRIANGULAR FACTORS SHOULD BE "KOO".'                 &
                    ,/,14X,',HOWEVER, IT IS NAMED "',A,'". CANNOT CONTINUE')

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(A,10X,A)

22345 FORMAT(3X,A,I8,' of ',I8,A)

! **********************************************************************************************************************************

      END SUBROUTINE SOLVE_GOA


      SUBROUTINE SOLVE_UO0

! Solves KOO*UO0 = PO for matrix UO0

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L2F, LINK2F, L2F_MSG, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FACTORED_MATRIX, FATAL_ERR, KOO_SDIA, NDOFO, NSUB, NTERM_KOO, NTERM_PO
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  PRTUO0, SOLLIB, SPARSE_FLAVOR
      USE SPARSE_MATRICES, ONLY       :  I_PO, J_PO, PO, I_KOO, J_KOO, KOO
      USE COL_VECS, ONLY              :  UO0_COL

! Interface module not needed for subr DPBTRS. This is "CONTAIN'ed" in module LAPACK_LIN_EQN_DPB, which is "USE'd" above

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN
      USE COL_VEC_LIFECYCLE, ONLY     :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
      USE SPARSE_CRS_ACCESS, ONLY     :  GET_SPARSE_CRS_COL
      USE LAPACK_ADAPTERS, ONLY       :  FBS_LAPACK
      USE SUPERLU_ADAPTERS, ONLY      :  FBS_SUPRLU
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_VECTOR

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLVE_UO0'
      CHARACTER( 1*BYTE)              :: NULL_COL          ! 'Y' if a col of KAO(transpose) is null
      CHARACTER(22*BYTE)              :: MODNAM1           ! Name to write to screen to describe module being run

      INTEGER(LONG)                   :: I,J               ! DO loop indices or counters
      INTEGER(LONG)                   :: INFO              ! Info on success of SuperLU solve
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN


      REAL(DOUBLE)                    :: NULL_SCALE_FACS(NDOFO)
                                                           ! LAPACK_S values not used so null this vector
      REAL(DOUBLE)                    :: INOUT_COL(NDOFO)  ! Temp variable for one col of load matrix PO



! **********************************************************************************************************************************
! Make units for writing errors the screen and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! For SOLLIB = BANDED, Subr REDUCE_KFF_TO_KAA should have done the decomp of KOO, so here we run FBS over the columns of PO to get
! UO0. First make sure that ABAND (KOO triangular factor) was successfully generated when SOLVE_GOA was run.

      IF (SOLLIB == 'BANDED  ') THEN
         IF (FACTORED_MATRIX(1:3) /= 'KOO') THEN
            WRITE(ERR,2504) SUBR_NAME, FACTORED_MATRIX
            WRITE(F06,2504) SUBR_NAME, FACTORED_MATRIX
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ENDIF

! Open file for writing UO0
                                                           ! Write GOA matrix to file L2F
      CALL FILE_OPEN ( L2F, LINK2F, OUNT, 'REPLACE', L2F_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

! **********************************************************************************************************************************
! Solve for UO0 by looping on columns of PO ("loads") to get columns of UO0 ("displs")

!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      CALL ALLOCATE_COL_VEC ( 'UO0_COL', NDOFO, SUBR_NAME )
      DO J = 1,NSUB

         DO I=1,NDOFO
            INOUT_COL(I)       = ZERO                      ! Initialize INOUT_COL since GET_SPARSE_CRS_COL won't
            NULL_SCALE_FACS(I) = ZERO                      ! Initialize these (not used here snce KOO not equilibrated
         ENDDO
         NULL_COL = 'Y'
         CALL GET_SPARSE_CRS_COL ( 'PO',J, NTERM_PO, NDOFO, NSUB, I_PO, J_PO, PO, ONE, INOUT_COL, NULL_COL )

         IF (NULL_COL == 'N') THEN                         ! Solve for UO0_COL

            CALL OURTIM
            MODNAM1 = '    Solve for UO0 col '
            WRITE(SC1,12345,ADVANCE='NO') MODNAM1, J, NSUB, CR13
                                                           ! FBS should not equilibrate since KOO was prevented from equilibrating
            IF      (SOLLIB == 'BANDED  ') THEN

               CALL FBS_LAPACK ( 'N', NDOFO, KOO_SDIA, NULL_SCALE_FACS, INOUT_COL )

            ELSE IF (SOLLIB == 'SPARSE  ') THEN

               IF (SPARSE_FLAVOR(1:7) == 'SUPERLU') THEN

                  INFO = 0
                  CALL FBS_SUPRLU ( SUBR_NAME, 'KOO', NDOFO, NTERM_KOO, I_KOO, J_KOO, KOO, J, INOUT_COL, INFO )

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

            DO I=1,NDOFO
               UO0_COL(I) = INOUT_COL(I)
            ENDDO

         ELSE

            DO I=1,NDOFO
               UO0_COL(I) = ZERO
            ENDDO

         ENDIF

! Write UO0 for this S/C to L2F

         DO I=1,NDOFO
            WRITE(L2F) UO0_COL(I)
         ENDDO

         IF (PRTUO0 > 0) THEN
            WRITE(F06,201) J
            CALL WRITE_VECTOR ( 'UO0 = KOO(-1)*PO', 'DISPL', NDOFO, UO0_COL)
         ENDIF

      ENDDO

      WRITE(SC1,*) CR13

      CALL DEALLOCATE_COL_VEC ( 'UO0_COL' )

      CALL FILE_CLOSE ( L2F, LINK2F, 'KEEP' )



      RETURN

! **********************************************************************************************************************************
  201 FORMAT(59X,'COLUMN',I6)

 2504 FORMAT(' *ERROR  2504: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NAME OF MATRIX THAT HAS BEEN DECOMPOSED INTO TRIANGULAR FACTORS SHOULD BE "KOO".'                 &
                    ,/,14X,' HOWEVER, IT IS NAMED "',A,'". CANNOT CONTINUE')

 9991 FORMAT(' *ERROR  9991: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,A, ' = ',A,' NOT PROGRAMMED ',A)

12345 FORMAT(3X,A,I8,' of ',I8,A)

! **********************************************************************************************************************************

      END SUBROUTINE SOLVE_UO0


      SUBROUTINE ALLOCATE_L2_GOA_2 ( CALLING_SUBR )

! Allocate some arrays for use in LINK

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPP6
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NTERM_GOA, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE SPARSE_MATRICES, ONLY       :  I2_GOA

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_L2_GOA_2'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format
      CHARACTER(24*BYTE)              :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: NROWS             ! Number of rows in array
      INTEGER(LONG), PARAMETER        :: NCOLS     = 1     ! Number of cols in array


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
! Allocate array I2_GOA
! NOTE: I2_GOA has row no's for ALL terms in GOA - not the sparse format that is compressed row storage.
! This is a temporary allocate for subr SOLVE_GOA

      MB_ALLOCATED = ZERO
      NROWS = NTERM_GOA
      JERR = 0

! Allocate array I2_GOA

      NAME = 'I2_GOA                  '
      IF (ALLOCATED(I2_GOA)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (I2_GOA(NTERM_GOA),STAT=IERR)
         IF (IERR == 0) THEN
            DO I=1,NTERM_GOA
               I2_GOA(I) = 0
            ENDDO
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,1699) TRIM(SUBR_NAME), CALLING_SUBR
         WRITE(F06,1699) SUBR_NAME, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      MB_ALLOCATED = REAL(DOUBLE)*REAL(NROWS)/ONEPP6
      CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )

      RETURN

! **********************************************************************************************************************************
  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_L2_GOA_2


      SUBROUTINE DEALLOCATE_L2_GOA_2

! Deallocate some arrays used in LINK2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE SPARSE_MATRICES, ONLY       :  I2_GOA

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_L2_GOA_2'
      CHARACTER(24*BYTE)              :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

! Deallocate array I2_GOA

      IF (ALLOCATED(I2_GOA)) THEN
         DEALLOCATE (I2_GOA,STAT=IERR)
         NAME = 'I2_GOA                  '
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
         ENDIF
      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
      CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )


      RETURN

! **********************************************************************************************************************************
  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)


! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_L2_GOA_2


   END MODULE REDUCTION_F_TO_A
