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

   MODULE LINK6_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: LINK6

   CONTAINS

      SUBROUTINE LINK6

! LINK6 generates Craig-Bampton model matrices and outputs those matrices to unformatted files as requested in Executive Control
! via OUTPUT4 entries.

! For a description of Craig-Bamptom analyses, see Appendix D to the MYSTRAN User's Referance Manual

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, ERR, F06, ERRSTAT, MOU4, SC1,                             &
                                         L2I    , L2K    , L2L    , L2M    , L2N    , L3A    ,OU4,                                 &
                                         LINK2I , LINK2K , LINK2L , LINK2M , LINK2N , LINK3A ,OU4FIL,                              &
                                         L2I_MSG, L2K_MSG, L2L_MSG, L2M_MSG, L2N_MSG, L3A_MSG,                                     &
                                         L2ISTAT, L2KSTAT, L2LSTAT, L2MSTAT, L2NSTAT, OU4STAT, OU4_MSG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, COMM, FATAL_ERR, LINKNO, MBUG, NDOFA, NDOFL, NDOFR, NUM_CB_DOFS, NUM_EIGENS,&
                                         NVEC, WARN_ERR

      USE SCONTR, ONLY                :  NTERM_CG_LTM, NTERM_DLR   , NTERM_PHIZL , NTERM_IF_LTM,                                   &
                                         NTERM_IRR   , NTERM_KLL   , NTERM_KRL   , NTERM_KRR   , NTERM_KRRcb ,                     &
                                         NTERM_KXX   , NTERM_MLL   , NTERM_MRL   , NTERM_MRN   , NTERM_MRR   , NTERM_MRRcb  ,      &
                                         NTERM_MXX   , NTERM_PHIXA

      USE CONSTANTS_1, ONLY           :  ONE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  CUSERIN, CUSERIN_XSET, PRTPHIXA, SUPWARN
      USE MODEL_STUF, ONLY            :  MEFFMASS_CALC, MPFACTOR_CALC
      USE CB_EFFECTIVE_PROPERTIES, ONLY:  CALC_CB_MEFM_MPF
      USE CB_TRANSFORMATIONS, ONLY    :  CALC_PHIZL, MERGE_PHIXA
      USE CB_LOAD_TRANSFORMATIONS, ONLY: INTERFACE_FORCE_LTM, NET_CG_LOADS_LTM, MERGE_LTM, SOLVE_DLR
      USE CB_REDUCED_MATRICES, ONLY   :  CALC_KRRCB, CALC_MRRCB, CALC_MRN, MERGE_KXX, MERGE_MXX
      USE EIGEN_MATRICES_1, ONLY      :  GEN_MASS, EIGEN_VAL, EIGEN_VEC
      USE OUTPUT4_MATRICES
      USE LAPACK_DPB_MATRICES, ONLY   :  ABAND
      USE SPARSE_MATRICES, ONLY       :  I_CG_LTM, J_CG_LTM, CG_LTM, I_DLR   , J_DLR   , DLR   , I_PHIZL , J_PHIZL , PHIZL ,       &
                                         I_IF_LTM, J_IF_LTM, IF_LTM, I_IRR   , J_IRR   , IRR   , I_KLL   , J_KLL   , KLL   ,       &
                                         I_KRL   , J_KRL   , KRL   , I_KRR   , J_KRR   , KRR   ,                                   &
                                         I_KRRcb , J_KRRcb , KRRcb , I_KXX   , J_KXX   , KXX   ,                                   &
                                         I_MLL   , J_MLL   , MLL   , I_MRL   , J_MRL   , MRL   , I_MRN   , J_MRN   , MRN   ,       &
                                         I_MRRcb , J_MRRcb , MRRcb , I_MRR   , J_MRR   , MRR   , I_MXX   , J_MXX   , MXX   ,       &
                                         I_PHIXA , J_PHIXA , PHIXA

      USE CB_TRANSFORMATIONS, ONLY      :  CALC_PHIZL, MERGE_PHIXA
      USE CB_LOAD_TRANSFORMATIONS, ONLY:  INTERFACE_FORCE_LTM, NET_CG_LOADS_LTM, MERGE_LTM, SOLVE_DLR
      USE CB_REDUCED_MATRICES, ONLY   :  CALC_KRRCB, CALC_MRRCB, CALC_MRN, MERGE_KXX, MERGE_MXX
      USE CB_EFFECTIVE_PROPERTIES, ONLY:  CALC_CB_MEFM_MPF
      USE L6_WORKSPACE, ONLY          :  ALLOCATE_L6_2
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_ARRAY_ALLOC_STAT, LINK_MESSAGE, WRITE_ALLOC_MEM_TABLE

      USE SPARSE_CRS_ACCESS, ONLY     :  GET_SPARSE_CRS_COL
      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_VEC
      USE DATE_TIME_UTILS, ONLY       :  OURDAT, OURTIM, TIME_INIT
      USE EIGEN_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_EIGEN1_MAT, DEALLOCATE_EIGEN1_MAT
      USE LAPACK_MATRIX_LIFECYCLE, ONLY:  DEALLOCATE_LAPACK_MAT
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_INQUIRE, FILE_OPEN, READERR
      USE TEMP_FILE_READERS, ONLY     :  READ_L1A, READ_L1M
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_USERIN_BD_CARDS
      USE MATRIX_FILE_IO, ONLY        :  READ_MATRIX_1, WRITE_SPARSE_CRS
      USE OUTPUT4_FILE_IO, ONLY       :  OUTPUT4_PROC
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE, WRITE_L1A
      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK6'
      CHARACTER(  1*BYTE)             :: CLOSE_IT            ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER(  1*BYTE)             :: READ_NTERM          ! 'Y' or 'N' Input to subr READ_MATRIX_1
      CHARACTER(  1*BYTE)             :: NULL_COL            ! = 'Y' if col returned from subr GET_SPARSE_CRS_COL is null
      CHARACTER(  1*BYTE)             :: OPND                ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to open  a file or not

      INTEGER(LONG)                   :: I,J                 ! DO loop indices
      INTEGER(LONG)                   :: IERROR    = 0       ! Error count when reading records from a file.
      INTEGER(LONG)                   :: IOCHK               ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: NUM_SOLNS           ! No. of solutions to process (e.g. NSUB for STATICS)
      INTEGER(LONG)                   :: OUNT(2)             ! File units to write messages to. Input to subr UNFORMATTED_OPEN.
      INTEGER(LONG)                   :: PART_VEC_A_LR(NDOFA)! Partitioning vector (N set into F and S sets)
      INTEGER(LONG), PARAMETER        :: P_LINKNO  = 4       ! Prior LINK no's that should have run before this LINK can execute
      INTEGER(LONG)                   :: REC_NO              ! Record number when reading a file.

      REAL(DOUBLE)                    :: PHIZL_COL(NDOFL)    ! Full column from sparse PHIZL matrix

! **********************************************************************************************************************************
      LINKNO = 6

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
! Set NUM_CB_DOFS (since it was initialized as 0 in SCONTR and hasn't been calc'd yet, must do this AFTER we call READ_L1A)

      NUM_CB_DOFS = 2*NDOFR + NVEC

! Check COMM for successful completion of prior LINKs

      IF (COMM(P_LINKNO) /= 'C') THEN
         WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
         WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Prior LINK's didn't complete, so quit
      ENDIF

! Read eigenvalues. EIGEN_VAL was not deallocated in LINK4 (see LINK4 comment 01/11/19) so we do not allocate it here anymore


!xx   CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', NUM_EIGENS, 1, SUBR_NAME )
      CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM' , NUM_EIGENS, 1, SUBR_NAME )
      CALL ALLOCATE_EIGEN1_MAT ( 'GEN_MASS' , NUM_EIGENS, 1, SUBR_NAME )
      CALL READ_L1M ( IERROR )

!!Read KLL stiffness matrix
!xx
!xx   CALL ALLOCATE_SPARSE_MAT ( 'KLL', NDOFL, NTERM_KLL, SUBR_NAME )
!xx
!xx   IF (NTERM_KLL > 0) THEN                              ! Allocate and read KGG stiffness matrix
!xx
!xx      READ_NTERM = 'Y'
!xx      OPND       = 'N'
!xx      CLOSE_IT   = 'Y'
!xx      CALL LINK_MESSAGE('READ KLL STIFFNESS MATRIX                 ')
!xx      CALL READ_MATRIX_1 ( LINK2G, L2G, OPND, CLOSE_IT, L2GSTAT, L2G_MSG, 'KLL', NTERM_KLL, READ_NTERM, NDOFL                   &
!xx                         , I_KLL, J_KLL, KLL)
!xx   ENDIF

! Read KRL stiffness matrix

      CALL ALLOCATE_SPARSE_MAT ( 'KRL', NDOFR, NTERM_KRL, SUBR_NAME )

      IF (NTERM_KRL > 0) THEN                              ! Allocate and read KGG stiffness matrix

         READ_NTERM = 'Y'
         OPND       = 'N'
         CLOSE_IT   = 'Y'
         CALL LINK_MESSAGE('READ KRL STIFFNESS MATRIX                 ')
         CALL READ_MATRIX_1 ( LINK2K, L2K, OPND, CLOSE_IT, L2KSTAT, L2K_MSG, 'KRL', NTERM_KRL, READ_NTERM, NDOFR                   &
                            , I_KRL, J_KRL, KRL)
      ENDIF

! Read KRR stiffness matrix

      CALL ALLOCATE_SPARSE_MAT ( 'KRR', NDOFR, NTERM_KRR, SUBR_NAME )

      IF (NTERM_KRR > 0) THEN                              ! Allocate and read KGG stiffness matrix

         READ_NTERM = 'Y'
         OPND       = 'N'
         CLOSE_IT   = 'Y'
         CALL LINK_MESSAGE('READ KRR STIFFNESS MATRIX                 ')
         CALL READ_MATRIX_1 ( LINK2L, L2L, OPND, CLOSE_IT, L2LSTAT, L2L_MSG, 'KRR', NTERM_KRR, READ_NTERM, NDOFR                   &
                            , I_KRR, J_KRR, KRR)
      ENDIF

!xxRead MLL mass matrix
!xx
!xx   CALL ALLOCATE_SPARSE_MAT ( 'MLL', NDOFL, NTERM_MLL, SUBR_NAME )
!xx
!xx   IF (NTERM_MLL > 0) THEN                              ! Allocate and read KGG stiffness matrix
!xx
!xx      READ_NTERM = 'Y'
!xx      OPND       = 'N'
!xx      CLOSE_IT   = 'Y'
!xx      CALL LINK_MESSAGE('READ MLL MASS MATRIX                 ')
!xx      CALL READ_MATRIX_1 ( LINK2I, L2I, OPND, CLOSE_IT, 'KEEP', L2I_MSG, 'MLL', NTERM_MLL, READ_NTERM, NDOFL                   &
!xx                         , I_MLL, J_MLL, MLL)
!xx   ENDIF

! Read MRL mass matrix

      CALL ALLOCATE_SPARSE_MAT ( 'MRL', NDOFR, NTERM_MRL, SUBR_NAME )

      IF (NTERM_MRL > 0) THEN                              ! Allocate and read KGG stiffness matrix

         READ_NTERM = 'Y'
         OPND       = 'N'
         CLOSE_IT   = 'Y'
         CALL LINK_MESSAGE('READ MRL MASS MATRIX                 ')
         CALL READ_MATRIX_1 ( LINK2M, L2M, OPND, CLOSE_IT, L2MSTAT, L2M_MSG, 'MRL', NTERM_MRL, READ_NTERM, NDOFR                   &
                            , I_MRL, J_MRL, MRL)
      ENDIF

! Read MRR mass matrix

      CALL ALLOCATE_SPARSE_MAT ( 'MRR', NDOFR, NTERM_MRR, SUBR_NAME )

      IF (NTERM_MRR > 0) THEN                              ! Allocate and read KGG stiffness matrix

         READ_NTERM = 'Y'
         OPND       = 'N'
         CLOSE_IT   = 'Y'
         CALL LINK_MESSAGE('READ MRR MASS MATRIX                 ')
         CALL READ_MATRIX_1 ( LINK2N, L2N, OPND, CLOSE_IT, L2NSTAT, L2N_MSG, 'MRR', NTERM_MRR, READ_NTERM, NDOFR                   &
                            , I_MRR, J_MRR, MRR)
      ENDIF

! Open file that has L-set eigenvectors and read them

      CALL FILE_OPEN ( L3A, LINK3A, OUNT, 'OLD', L3A_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
      CALL LINK_MESSAGE('READ EIGENVECTORS FROM FILE')
      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC', NDOFL, NVEC, SUBR_NAME )
      DO J=1,NVEC
         DO I=1,NDOFL
            READ(L3A,IOSTAT=IOCHK) EIGEN_VEC(I,J)
            IF (IOCHK /= 0) THEN
               CALL READERR ( IOCHK, LINK3A, L3A_MSG, REC_NO, OUNT )
               IERROR = IERROR + 1
            ENDIF
         ENDDO
      ENDDO
      CALL FILE_CLOSE ( L3A, LINK3A, 'KEEP' )

      IF (IERROR > 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Solve for DLR

      CALL LINK_MESSAGE('SOLVE FOR DLR              ')
      IF ((NTERM_KLL > 0) .AND. (NTERM_KRL > 0)) THEN
         CALL SOLVE_DLR
      ELSE
         NTERM_DLR = 0
         CALL ALLOCATE_SPARSE_MAT ( 'DLR' , NDOFL, NTERM_DLR, SUBR_NAME )
         CALL ALLOCATE_SPARSE_MAT ( 'DLRt', NDOFR, NTERM_DLR, SUBR_NAME )
      ENDIF

! Merge an RxR identity matrix (IRR), an RxN null matrix, DLR and L-set eigenvector matrix to get CB transformation matrix PHIXA:

!                        |IRR  0RN |     IRR = R-set identity matrix, 0RN  = RxN null matrix (N = num modes)
!                PHIXA = |         |
!                        |DLR  PHIL|     DLR = C.B. boundary modes  , PHIL = EIGEN_VECS() LxN eigenvectors for the L-set x N modes

      NTERM_IRR = NDOFR                                    ! First, allocate and set identity matrix IRR
      CALL ALLOCATE_SPARSE_MAT ( 'IRR', NDOFR, NTERM_IRR, 'LINK6' )
      DO I=1,NDOFR+1
         I_IRR(I) = I
      ENDDO
      DO I=1,NDOFR
         J_IRR(I) = I
           IRR(I) = ONE
      ENDDO

      CALL PARTITION_VEC ( NDOFA, 'A ', 'L ', 'R ', PART_VEC_A_LR )
      NTERM_PHIXA = NDOFR + NTERM_DLR + NDOFL*NVEC
      CALL ALLOCATE_SPARSE_MAT ( 'PHIXA', NDOFA, NTERM_PHIXA, 'LINK6' )
      CALL LINK_MESSAGE('MERGE MATRICES INTO PHIXA  ')
      CALL MERGE_PHIXA ( PART_VEC_A_LR )
      IF (PRTPHIXA > 0) THEN
         CALL WRITE_SPARSE_CRS ( 'PHIXA','A ','  ', NTERM_PHIXA, NDOFA, I_PHIXA, J_PHIXA, PHIXA )
      ENDIF
!xx   WRITE(SC1, * )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate IRR      ', CR13
      CALL DEALLOCATE_SPARSE_MAT ( 'IRR')

! Calculate L-set PHIZL (cols are L-set CB vecs that will be written to L3A for LINK5 processing to G-set size)

      CALL LINK_MESSAGE('CALC DISPLACEMENT PHIZL    ')
      CALL CALC_PHIZL

! Calc KXX, the CB stiffness matrix (with DOF's: R-set displs and modal DOF's)

      CALL LINK_MESSAGE('CALC KRRcb                 ')
      CALL CALC_KRRcb

      CALL LINK_MESSAGE('MERGE MATRICES INTO KXX    ')
      CALL MERGE_KXX

! Calc MXX, the CB mass matrix (with DOF's: R-set displs and modal DOF's)

      CALL LINK_MESSAGE('CALC MRRcb                 ')
      CALL CALC_MRRcb

      CALL LINK_MESSAGE('CALC MRN                   ')
      CALL CALC_MRN

!xx   WRITE(SC1, * )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate DLRt     ', CR13
      CALL DEALLOCATE_SPARSE_MAT ( 'DLRt' )

      CALL LINK_MESSAGE('MERGE MATRICES INTO MXX    ')
      CALL MERGE_MXX

! Calculate LTM for interface forces

      CALL LINK_MESSAGE('CALC INTERFACE FORCES      ')
      CALL INTERFACE_FORCE_LTM

! Calculate LTM for net CG accels

      CALL LINK_MESSAGE('CALC NET CG LOADS          ')
      CALL NET_CG_LOADS_LTM

!xx   WRITE(SC1, * )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MRRcbn   ', CR13
      CALL DEALLOCATE_SPARSE_MAT ( 'MRRcbn' )

! Merge CG_LTM and IF_LTM into overall LTM

      CALL LINK_MESSAGE('MERGE LTM                  ')
      CALL MERGE_LTM

! Rewind L3A and write PHIZL to it and then close L3A

      CALL LINK_MESSAGE('WRITE L-set PHIZL to FILE')
      CALL FILE_OPEN ( L3A, LINK3A, OUNT, 'OLD', L3A_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
      NUM_SOLNS = NUM_CB_DOFS
      DO I=1,NUM_SOLNS
         CALL GET_SPARSE_CRS_COL ( 'PHIZL', I, NTERM_PHIZL, NDOFL, NUM_SOLNS, I_PHIZL, J_PHIZL, PHIZL, ONE, PHIZL_COL, NULL_COL )
         DO J=1,NDOFL
            WRITE(L3A) PHIZL_COL(J)
         ENDDO
      ENDDO
      CALL FILE_CLOSE ( L3A, LINK3A, 'KEEP' )

! Calc modal participation factors and modal mass

      IF ((MEFFMASS_CALC == 'Y') .OR. (MPFACTOR_CALC == 'Y')) THEN
         CALL LINK_MESSAGE('CALC MPF AND MEFFMASS')
         CALL CALC_CB_MEFM_MPF
      ENDIF

! Call OUTPUT4 processor to process output requests for OUTPUT4 matrices generated in this link

      IF (NUM_OU4_REQUESTS > 0) THEN
         CALL LINK_MESSAGE('WRITE OUTPUT4 NATRICES      ')
         WRITE(F06,*)
         CALL OUTPUT4_PROC ( SUBR_NAME )
      ENDIF

! If user has requested, write out B.D. CUSERIN card images for use in INPUT4 for this substructure

      IF (CUSERIN == 'Y') THEN
         IF (CUSERIN_XSET(1:2) /= 'R ') THEN
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,6888) CUSERIN_XSET
            IF (SUPWARN == 'N') THEN
               WRITE(F06,6888) CUSERIN_XSET
            ENDIF
            CUSERIN_XSET = 'R '                            ! Since we are in LINK6, CUSERIN_XSET must be the R-set
         ENDIF
         CALL WRITE_USERIN_BD_CARDS ( NDOFR, CUSERIN_XSET )
      ENDIF

! Deallocate arrays, but not PHIXA which may be needed in LINK5 to expand to G-set or IF_LTM needed in LINK9 for calc of QR_COL.
! Also do not deallocate EIGEN_VAL, it may be needed later

!xx   WRITE(SC1, * )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate CG_LTM   ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'CG_LTM'    )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate DLR      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'DLR'       )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PHIZL    ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'PHIZL'     )
!xx   WRITE(SC1,12345,ADVANCE='NO') '       Deallocate EIGEN_VAL', CR13  ;   CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate EIGEN_VEC', CR13  ;   CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VEC' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GEN_MASS ', CR13  ;   CALL DEALLOCATE_EIGEN1_MAT ( 'GEN_MASS'  )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MODE_NUM ', CR13  ;   CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM'  )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLL      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KLL'       )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KLLs     ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KLLs'      )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KRL      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KRL'       )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KRR      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KRR'       )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KRRcb    ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KRRcb'     )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KXX      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KXX'       )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate LTM      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'LTM'       )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MLL      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MLL'       )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MRL      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MRL'       )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MRN      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MRN'       )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MRR      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MRR'       )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MRRcb    ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MRRcb'     )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MXX      ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MXX'       )

      CALL DEALLOCATE_LAPACK_MAT ( 'ABAND')

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

! Write LINK6 end to F06

      CALL OURTIM
      WRITE(F06,151) LINKNO

! Close files

      IF (( DEBUG(193) == 6) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK6' )
      ENDIF

! Write LINK4 end to screen

      WRITE(SC1,153) LINKNO

! **********************************************************************************************************************************
  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END',/)

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  153 FORMAT(  ' >> LINK',I3,' END')

  290 FORMAT(23X,5A)

 6888 FORMAT(' *WARNING    : PARAMETER CUSERIN_XSET WAS READ FROM THE B.D. PARAM CUSERIN ENTRY AS "',A,'". IT HAS BEEN RESET TO'   &
                    ,/,14X,'"R " FOR THIS CB MODEL GENERATION RUN')

 9998 FORMAT(' *ERROR  9998: COMM ',I3,' INDICATES UNSUCCESSFUL LINK ',I2,' COMPLETION.'                                           &
                    ,/,14X,' FATAL ERROR - CANNOT START LINK ',I2)

12345 FORMAT(A,10X,A)






99887 format(32767(1es14.6))

! ##################################################################################################################################

      END SUBROUTINE LINK6

   END MODULE LINK6_MOD
