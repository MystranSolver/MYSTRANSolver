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

   MODULE LINK9_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: LINK9

   CONTAINS

      SUBROUTINE LINK9 ( LK9_PROC_NUM )

! Main driver for calculating outputs requested in Case Control once the G-set unknowns have been solved for in prior LINK's

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR

      USE IOUNT1, ONLY                :  ERR, F06, F25, L1E, L1M, L1R, L2A, L2B, L2C, L2D, L2I, L2J, L2R, L2S,                     &
                                         L5A, L5B, NEU, OT4, OU4, SC1

      USE IOUNT1, ONLY                :  F06FIL, F25FIL, LINK1B, LINK1E, LINK1M, LINK1R, LINK2A, LINK2B, LINK2C, LINK2D,           &
                                         LINK2I, LINK2J, LINK2R, LINK2S, LINK5A, LINK5B, MOT4  , MOU4  , NEUFIL, OT4FIL, OU4FIL

      USE IOUNT1, ONLY                :  L1ASTAT, L1ESTAT, L1MSTAT, L1RSTAT, L2ASTAT, L2BSTAT, L2CSTAT, L2ISTAT, L2JSTAT, L2RSTAT, &
                                         L2SSTAT, OT4STAT, OU4STAT

      USE IOUNT1, ONLY                :  F25_MSG, L1E_MSG, L1M_MSG, L1R_MSG, L2A_MSG, L2B_MSG, L2C_MSG, L2D_MSG, L2I_MSG,          &
                                         L2J_MSG, L2R_MSG, L2S_MSG, L5A_MSG, L5B_MSG, NEU_MSG,                                     &
                                         OT4_MSG, OU4_MSG, OT4_GRD_OTM, OT4_ELM_OTM, OU4_GRD_OTM, OU4_ELM_OTM

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_ENTRY_LEN, COMM, IBIT, INT_EIG_NUM, INT_SC_NUM, JTSUB, FATAL_ERR,        &
                                         FEMAP_VERSION, LINKNO, MBUG,                                                              &
                                         NDOFF, NDOFG, NDOFL, NDOFM, NDOFN, ndofo, NDOFR, NDOFS, NDOFSA, NGRID, NSUB, NVEC,        &
                                         NTERM_IF_LTM, NTERM_GMN, NTERM_HMN, NTERM_KFS, NTERM_KFSD, NTERM_LMN, NTERM_MFS,          &
                                         NTERM_MGG, NTERM_MLL,NTERM_PG, NTERM_PM, NTERM_PS, NTERM_QSYS,                            &
                                         NUM_BUCKLING_SUBS, NUM_CB_DOFS, NUM_EIGENS,                                               &
                                         NROWS_OTM_ACCE, NROWS_OTM_DISP, NROWS_OTM_MPCF, NROWS_OTM_SPCF,                           &
                                         NROWS_OTM_ELFE, NROWS_OTM_ELFN, NROWS_OTM_STRE, NROWS_OTM_STRN,                           &
                                         NROWS_TXT_ACCE, NROWS_TXT_DISP, NROWS_TXT_MPCF, NROWS_TXT_SPCF,                           &
                                         NROWS_TXT_ELFE, NROWS_TXT_ELFN, NROWS_TXT_STRE, NROWS_TXT_STRN, RESTART, SOL_NAME,        &
                                         WARN_ERR, MODE_SUBCASE

      USE SCONTR, ONLY                :  GROUT_ACCE_BIT, GROUT_DISP_BIT, GROUT_OLOA_BIT, GROUT_SPCF_BIT, GROUT_MPCF_BIT,           &
                                         GROUT_GPFO_BIT, ELOUT_ELFN_BIT, ELOUT_ELFE_BIT, ELOUT_STRE_BIT, ELOUT_STRN_BIT,           &
                                         ELDT_F25_U_P_BIT

      USE CC_OUTPUT_DESCRIBERS, ONLY  :  DISP_F06, ACCE_F06, OLOA_F06, SPCF_F06, MPCF_F06, FORC_F06, GPFO_F06, STRE_F06, STRN_F06
      USE TIMDAT, ONLY                :  STIME
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  EPSIL, MPFOUT, SUPINFO, SUPWARN, WTMASS, PRTNEU, POST
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE COL_VECS, ONLY              :  FG_COL, UG_COL, PG_COL, PM_COL, PS_COL, QSYS_COL, QGm_COL, QGr_COL, QGs_COL, QR_COL,      &
                                         PHIXG_COL, PHIXN_COL
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, GEN_MASS, MODE_NUM
      USE OUTPUT4_MATRICES, ONLY      :  NUM_OU4_REQUESTS, OU4_PART_MAT_NAMES, HAS_OU4_MAT_BEEN_PROCESSED, OU4_PART_MAT_NAMES
      USE OUTPUT4_MATRICES, ONLY      :  OTM_ACCE, OTM_DISP, OTM_MPCF, OTM_SPCF, OTM_ELFE, OTM_ELFN, OTM_STRE, OTM_STRN,           &
                                         TXT_ACCE, TXT_DISP, TXT_MPCF, TXT_SPCF, TXT_ELFE, TXT_ELFN, TXT_STRE, TXT_STRN

      USE SPARSE_MATRICES, ONLY       :  I_GMN , J_GMN , GMN , I_GMNt, J_GMNt, GMNt, I_HMN , J_HMN , HMN ,                         &
                                         I_KSF , J_KSF , KSF , I_KSFD, J_KSFD, KSFD, I_LMN , J_LMN , LMN ,                         &
                                         I_MGG , J_MGG , MGG , I_MLL , J_MLL , MLL , I_MSF , J_MSF , MSF ,                         &
                                         I_PG  , J_PG  , PG  , I_PM  , J_PM  , PM  , I_PS  , J_PS  , PS  , I_QSYS, J_QSYS, QSYS

      USE SPARSE_MATRICES, ONLY       :  I_IF_LTM, J_IF_LTM, IF_LTM, SYM_MGG, SYM_MSF, SYM_PG, SYM_PM

      USE DOF_TABLES, ONLY            :  TDOF

      USE MODEL_STUF, ONLY            :  ANY_ACCE_OUTPUT, ANY_DISP_OUTPUT, ANY_MPCF_OUTPUT, ANY_SPCF_OUTPUT, ANY_OLOA_OUTPUT,      &
                                         ANY_GPFO_OUTPUT, ANY_ELFE_OUTPUT, ANY_ELFN_OUTPUT, ANY_STRE_OUTPUT, ANY_STRN_OUTPUT,      &
                                         OELDT, OELOUT, OGROUT, GRID, GROUT, MEFFMASS_CALC, MPFACTOR_CALC, SCNUM, SUBLOD, TITLE,   &
                                         STITLE, LABEL
      USE LINK9_STUFF, ONLY           :  MAXREQ

      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURDAT, OURTIM, TIME_INIT
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_INQUIRE, FILE_OPEN, READERR
      USE TEMP_FILE_READERS, ONLY     :  READ_L1A, READ_L1M
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE, WRITE_L1A
      USE MODEL_STORAGE_DEALLOCATION, ONLY:  DEALLOCATE_MODEL_STUF
      USE MODEL_STORAGE_ALLOCATION, ONLY:  ALLOCATE_MODEL_STUF
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  READ_MATRIX_1
      USE MATRIX_PARTITIONING, ONLY   :  PARTITION_SS, PARTITION_SS_NTERM, PARTITION_VEC
      USE COL_VEC_LIFECYCLE, ONLY     :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
      USE SPARSE_CRS_ACCESS, ONLY     :  GET_SPARSE_CRS_COL
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATTRNSP_SS
      USE EIGEN_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_EIGEN1_MAT, DEALLOCATE_EIGEN1_MAT
      USE RIGID_BODY_STORAGE_LIFECYCLE, ONLY:  ALLOCATE_CB_ELM_OTM, ALLOCATE_CB_GRD_OTM
      USE GRID_RESULT_DISPATCH, ONLY  :  OFP1
      USE MODAL_RESULT_DISPATCH, ONLY :  OFP2
      USE SPARSE_FULL_MULTIPLICATION, ONLY:  MATMULT_SFF
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE GRID_POINT_FORCE_BALANCE, ONLY:  GP_FORCE_BALANCE_PROC
      USE ELEMENT_RESULT_DISPATCH, ONLY:  OFP3
      USE OUTPUT4_FILE_IO, ONLY       :  OUTPUT4_PROC, WRITE_OU4_FULL_MAT
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_ELM_OT4, WRITE_GRD_OT4
      USE MODAL_OUTPUT_WRITERS, ONLY  :  WRITE_MEFFMASS, WRITE_MPFACTOR
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE IN4_FILE_LIFECYCLE, ONLY    :  DEALLOCATE_IN4_FILES
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_ARRAY_ALLOC_STAT, LINK_MESSAGE, LINK_MESSAGE_I, WRITE_ALLOC_MEM_TABLE
      USE MATRIX_MERGING, ONLY        :  MERGE_COL_VECS

      USE DOF_TABLE_LIFECYCLE, ONLY   :  DEALLOCATE_DOF_TABLES
      IMPLICIT NONE

      LOGICAL                         :: WRITE_NEU         ! flag
      LOGICAL                         :: LEXIST            ! .TRUE. if a file exists
      LOGICAL                         :: LOPEN             ! .TRUE. if a file is opened

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK9'
      CHARACTER(LEN=3*CC_ENTRY_LEN+5) :: TSL               ! Concatenated TITLE, STITLE, LABEL for FEMAP block 450 in FEMAP NEU file
      CHARACTER( 1*BYTE)              :: CLOSE_IT          ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER( 8*BYTE)              :: CLOSE_STAT        ! What to do with file when it is closed
      CHARACTER(14*BYTE)              :: CTIME             ! A char variable to which STIME will be written (for use in NEU file)
      CHARACTER( 6*BYTE)              :: FEMAP_BLK='xxxxxx'! 3 digit number indicating the FEMAP data block
      CHARACTER( 1*BYTE)              :: NULL_ROW          ! 'Y'/'N' depending on whether a col in IF_LTM is null
      CHARACTER( 1*BYTE)              :: ZERO_GEN_STIFF    ! Indicator of whether there are zero gen stiffs (can't calc MEFFMASS)

      CHARACTER(24*BYTE)              :: MESSAG            ! File description. Input to subr UNFORMATTED_OPEN
      CHARACTER( 1*BYTE)              :: NULL_COL          ! An output from subr GET_SPARSE_CRS_COL
      CHARACTER( 1*BYTE)              :: PROC_PG_OUTPUT    ! 'Y' in general. However, for BUCKLING, set to 'N' for eigen subcase
      CHARACTER( 1*BYTE)              :: READ_SPCARRAYS    ! ='Y' if we need to read KSF, etc. See test below.

      INTEGER(LONG), INTENT(IN)       :: LK9_PROC_NUM      ! 2 if this is the LINK9 call for the linear buckling step of
!                                                            SOL_NAME = 'BUCKLING. Otherwise 1 to designate that, for BUCKLING,
!                                                            this call to LINK9 is for the linear statics (1st) portion of BUCKLING

      INTEGER(LONG)                   :: ANY_U_P_OUTPUT    ! > 0 if requests for output of elem loads/displs in a any S/C
      INTEGER(LONG)                   :: COL_NUM           ! Col number to get when subr GET_SPARSE_CRS_COL is called
      INTEGER(LONG)                   :: FEMAP_SET_ID  = 0 ! Set ID for FEMAP output
      INTEGER(LONG)                   :: FORM              ! Matrix format
      INTEGER(LONG)                   :: GDOF              ! G-set DOF number
      INTEGER(LONG)                   :: RDOF              ! R-set DOF number
      INTEGER(LONG)                   :: G_SET_COL         ! Col number in TDOF where the G-set DOF's exist
      INTEGER(LONG)                   :: R_SET_COL         ! Col number in TDOF where the R-set DOF's exist
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: ITE       = 0     ! Index (1 thru MOT4) for file unit num for elem OTM's text files
      INTEGER(LONG)                   :: ITG       = 0     ! Index (1 thru MOU4) for file unit num for grid OTM's text files
      INTEGER(LONG)                   :: IUE       = 0     ! Index (1 thru MOT4) for file unit num for elem OTM's unformatted files
      INTEGER(LONG)                   :: IUG       = 0     ! Index (1 thru MOU4) for file unit num for grid OTM's unformatted files
      INTEGER(LONG)                   :: IERROR            ! Error count
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: JVEC              ! DO loop index - output vector no. being processed (S/C or eigenvec no.)
      INTEGER(LONG), PARAMETER        :: NUM1      = 1     ! Used in subr's that partition matrices
      INTEGER(LONG), PARAMETER        :: NUM2      = 2     ! Used in subr's that partition matrices
      INTEGER(LONG)                   :: NUM_COLS          ! Number of cols to get when subr GET_SPARSE_CRS_COL is called
      INTEGER(LONG)                   :: NUM_SOLNS         ! No. of solutions to process (e.g. NSUB for STATICS)
      INTEGER(LONG)                   :: NUM_OU4_NOT_PART  ! Number of OU4 mats requested for partitioning that were not done
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: OT4_EROW  = 0     ! Row number in OT4 elem related files. Accumulated in OFP1,2 for OTM's
      INTEGER(LONG)                   :: OT4_GROW  = 0     ! Row number in OT4 grid related files. Accumulated in OFP1,2 for OTM's
      INTEGER(LONG)                   :: PART_G_NM(NDOFG)  ! Partitioning vector (G set into N and M sets)
      INTEGER(LONG)                   :: PART_SUB(NSUB)    ! Partitioning vector (1's for all subcases)
      INTEGER(LONG)                   :: P_LINKNO          ! Prior LINK no's that should have run before this LINK can execute
      INTEGER(LONG)                   :: PM_ROW_MAX_TERMS  ! Output from subr PARTITION_SIZE (max terms in any row of matrix)
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: SC_ACCE_OUTPUT    ! = 1 if requests for output of accels in a particular S/C
      INTEGER(LONG)                   :: SC_DISP_OUTPUT    ! = 1 if requests for output of displs in a particular S/C
      INTEGER(LONG)                   :: SC_OLOA_OUTPUT    ! = 1 if requests for output of applied loads in a particular S/C
      INTEGER(LONG)                   :: SC_SPCF_OUTPUT    ! = 1 if requests for output of SPC forces in a particular S/C
      INTEGER(LONG)                   :: SC_MPCF_OUTPUT    ! = 1 if requests for output of MPC forces in a particular S/C
      INTEGER(LONG)                   :: SC_GPFO_OUTPUT    ! = 1 if requests for output of G.P. force balance in a particular S/C
      INTEGER(LONG)                   :: SC_ELFN_OUTPUT    ! = 1 if requests for output of elem node forces in a particular S/C
      INTEGER(LONG)                   :: SC_ELFE_OUTPUT    ! = 1 if requests for output of elem engr forces in a particular S/C
      INTEGER(LONG)                   :: SC_STRE_OUTPUT    ! = 1 if requests for output of elem stresses in a particular S/C
      INTEGER(LONG)                   :: SC_STRN_OUTPUT    ! = 1 if requests for output of elem strains  in a particular S/C
      INTEGER(LONG)                   :: XTIME             ! Time stamp read from an unformatted file
      INTEGER(LONG)                   :: PREV_SC_NUM   = 0 ! Tracks previous INT_SC_NUM to detect subcase transitions for INT_EIG_NUM
      INTEGER(LONG)                   :: SC_VEC_COUNT  = 0 ! Running per-subcase eigenvector counter; reset on each subcase change


      REAL(DOUBLE)                    :: EPS1              ! Small number to compare against zero
      REAL(DOUBLE)                    :: UGV               ! A G-set vector read from file L5A
      REAL(DOUBLE)                    :: PHIXGV            ! A G-set vector read from file L5B
      INTEGER(LONG)                   :: ITABLE            !
      LOGICAL                         :: NEW_RESULT        ! Is this a new result

      INTRINSIC                       :: IAND

! **********************************************************************************************************************************
      LINKNO = 9

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

      EPS1 = EPSIL(1)

      WRITE_NEU = (PRTNEU == 'Y')

      ! Make units for writing errors the screen until we open output files
      OUNT(1) = SC1
      OUNT(2) = SC1

      ! Make units for writing errors the error file and output file
      OUNT(1) = ERR
      OUNT(2) = F06

      ! Write info to text files
      WRITE(ERR,150) LINKNO
      WRITE(F06,150) LINKNO

      ! Read LINK1A file
      CALL READ_L1A ( 'KEEP' )

      ! Check COMM for successful completion of prior LINKs
      IF (RESTART == 'Y') THEN
         P_LINKNO = 1
      ELSE
         P_LINKNO = 5
      ENDIF
      IF (COMM(P_LINKNO) /= 'C') THEN
         WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
         WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Prior LINK's didn't complete, so quit
      ENDIF

      ! Before reading file data in subr LINK9S, deallocate all of those arrays and then allocate them fresh
      CALL LINK_MESSAGE('DEALLOCATE ARRAYS BEFORE READING LINK9S')
                                                           ! Deallocate data in file LINK1D
      CALL DEALLOCATE_MODEL_STUF ( 'SCNUM' )
      CALL DEALLOCATE_MODEL_STUF ( 'TITLES' )
      CALL DEALLOCATE_MODEL_STUF ( 'SUBLOD' )
      CALL DEALLOCATE_MODEL_STUF ( 'GROUT, ELOUT' )
      CALL DEALLOCATE_MODEL_STUF ( 'ELDT' )
                                                           ! Deallocate data in file LINK1G
      CALL DEALLOCATE_MODEL_STUF ( 'ETYPE, EDAT, EPNT' )
      CALL DEALLOCATE_MODEL_STUF ( 'ESORT1' )
      CALL DEALLOCATE_MODEL_STUF ( 'ESORT2' )
      CALL DEALLOCATE_MODEL_STUF ( 'EOFF' )
      CALL DEALLOCATE_MODEL_STUF ( 'VVEC, OFFSETS, PLATE stuff' )
      CALL DEALLOCATE_MODEL_STUF ( 'ELEM PROPERTIES AND MATERIALS' )
                                                           ! Deallocate data in file LINK1K
      CALL DEALLOCATE_MODEL_STUF ( 'TPNT, TDATA' )
      CALL DEALLOCATE_MODEL_STUF ( 'GTEMP' )
                                                           ! Deallocate data in file LINK1Q
      CALL DEALLOCATE_MODEL_STUF ( 'PPNT, PDATA, PTYPE' )
      CALL DEALLOCATE_MODEL_STUF ( 'PLOAD4_3D_DATA' )

      CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR DATA READ IN LINK9S')
                                                           ! Allocate data to be read in LINK9S from file LINK1D
      CALL   ALLOCATE_MODEL_STUF ( 'SCNUM', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'TITLES', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'SUBLOD', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'GROUT, ELOUT', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'ELDT', SUBR_NAME )
                                                           ! Allocate data to be read in LINK9S from file LINK1G
      CALL   ALLOCATE_MODEL_STUF ( 'ETYPE, EDAT, EPNT', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'ESORT1', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'ESORT2', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'EOFF', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'VVEC, OFFSETS, PLATE stuff', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'ELEM PROPERTIES AND MATERIALS', SUBR_NAME )
                                                           ! Allocate data to be read in LINK9S from file LINK1K
      CALL   ALLOCATE_MODEL_STUF ( 'TPNT, TDATA', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'GTEMP', SUBR_NAME )
                                                           ! Allocate data to be read in LINK9S from file LINK1Q
      CALL   ALLOCATE_MODEL_STUF ( 'PPNT, PDATA, PTYPE', SUBR_NAME )
      CALL   ALLOCATE_MODEL_STUF ( 'PLOAD4_3D_DATA', SUBR_NAME )

      ! Read LINK9S data
      CALL LINK_MESSAGE('READ MODEL DATA ARRAYS')
      CALL LINK9S

      ! Determine MAXREQ (max number of output requests) so we can allocate memory to arrays below
      CALL MAXREQ_OGEL
      CALL DEALLOCATE_MODEL_STUF ( 'ESORT2' )
!-----------------------------------------------------------------------------------------------------------------------------------
! Read PG (G_set loads) if this is STATICS, NLSTATIC or BUCKLING (with LOAD_ISTEP = 1) and there are any output requests for applied
! loads or MPC forces or G.P. force balance.
! Note: need PG to partition PM if MPC force outout is requested and also need PG if G.P. force balance is requested.

      CALL ALLOCATE_SPARSE_MAT ( 'PG', NDOFG, NTERM_PG, SUBR_NAME )

      IF ((SOL_NAME(1:7)=='STATICS') .OR. (SOL_NAME(1:8)=='NLSTATIC') .OR. ((SOL_NAME(1:8)=='BUCKLING') .AND. (LOAD_ISTEP==1))) THEN

         IF ((ANY_OLOA_OUTPUT > 0) .OR. (ANY_MPCF_OUTPUT > 0) .OR. (ANY_GPFO_OUTPUT > 0) .OR. (WRITE_NEU)) THEN

            IF (NTERM_PG > 0) THEN

               CALL LINK_MESSAGE('ALLOCATING SPARSE ARRAYS FOR PG LOADS')

               CALL LINK_MESSAGE('READ PG LOADS')
               CLOSE_IT   = 'N'
               CALL READ_MATRIX_1 ( LINK1E, L1E, 'N', CLOSE_IT, 'KEEP', L1E_MSG, 'PG', NTERM_PG, 'Y', NDOFG,                       &
                                    I_PG, J_PG, PG)
               IF ((ANY_MPCF_OUTPUT > 0) .OR. (ANY_GPFO_OUTPUT > 0) .OR. (WRITE_NEU)) THEN
                  IF (NTERM_PM  > 0) THEN                  ! Partition PM from PG if there are any loads on the M-set
                     CALL PARTITION_VEC (NDOFG,'G ','N ','M ',PART_G_NM)
                     DO I=1,NSUB
                        PART_SUB = 1
                     ENDDO
                     CALL PARTITION_SS_NTERM ( 'PG' , NTERM_PG,  NDOFG, NSUB , SYM_PG , I_PG , J_PG ,      PART_G_NM, PART_SUB,    &
                                                NUM2, NUM1, PM_ROW_MAX_TERMS, 'PM', NTERM_PM, SYM_PM )
                     CALL ALLOCATE_SPARSE_MAT ( 'PM', NDOFM, NTERM_PM, SUBR_NAME )

                     CALL PARTITION_SS ( 'PG' , NTERM_PG , NDOFG, NSUB , SYM_PG , I_PG , J_PG , PG , PART_G_NM, PART_SUB,          &
                                          NUM2, NUM1, PM_ROW_MAX_TERMS, 'PM', NTERM_PM , NDOFM, SYM_PM, I_PM , J_PM , PM  )
                  ENDIF
               ENDIF

            ENDIF

         ENDIF

      ENDIF

      ! Read files with KSF, MSF, QSYS (used to calc SPC constraint forces, QS), but only if they will be needed.
      ! For any SOL_NAME they will be needed if any SPC constraint force output is requested or GP force balance or if WRITE_NEU.
      ! For non CB they will be needed also if MEFFMASS, MPFACTOR are to be calculated (done via SPC force total method)
      READ_SPCARRAYS = 'N'
      IF (SOL_NAME == 'GEN CB MODEL') THEN
         IF ((ANY_SPCF_OUTPUT > 0) .OR. (ANY_GPFO_OUTPUT > 0) .OR. (NDOFSA > 0) .OR. (WRITE_NEU)) THEN
            READ_SPCARRAYS = 'Y'
         ENDIF
      ELSE
         IF ((ANY_SPCF_OUTPUT > 0) .OR. (ANY_GPFO_OUTPUT > 0) .OR. (NDOFSA > 0) .OR. (WRITE_NEU) .OR.                          &
             (MEFFMASS_CALC == 'Y') .OR. (MPFACTOR_CALC == 'Y')) THEN
            READ_SPCARRAYS = 'Y'
         ENDIF
      ENDIF

      CALL ALLOCATE_SPARSE_MAT ( 'KSF' , NDOFS, NTERM_KFS , SUBR_NAME )
      CALL ALLOCATE_SPARSE_MAT ( 'KSFD', NDOFS, NTERM_KFSD, SUBR_NAME )
      CALL ALLOCATE_SPARSE_MAT ( 'MSF' , NDOFS, NTERM_MFS , SUBR_NAME )
      CALL ALLOCATE_SPARSE_MAT ( 'QSYS', NDOFS, NTERM_QSYS, SUBR_NAME )
      CALL ALLOCATE_SPARSE_MAT ( 'PS'  , NDOFS, NTERM_PS  , SUBR_NAME )
      CALL ALLOCATE_COL_VEC ('QSYS_COL',NDOFS,SUBR_NAME)! Alloc this here since OFP2 uses it (will be zero's if NTERM_QSYS = 0)

      IF (READ_SPCARRAYS == 'Y') THEN

         IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN

            IF (NTERM_KFSD > 0) THEN

               CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR, AND READ, KSFD')
               CALL LINK_MESSAGE('READ KSFD MATRIX')
               CLOSE_IT   = 'Y'
               CLOSE_STAT = 'KEEP'
               CALL READ_MATRIX_1 (LINK2B,L2B,'N',CLOSE_IT,CLOSE_STAT,L2B_MSG,'KSFD',NTERM_KFSD,'Y',NDOFS,I_KSFD,J_KSFD,KSFD)

            ENDIF

         ELSE

            IF (NTERM_KFS  > 0) THEN

               CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR, AND READ, KSF')
               CALL LINK_MESSAGE('READ KSF MATRIX')
               CLOSE_IT   = 'Y'
               CLOSE_STAT = 'KEEP'
               CALL READ_MATRIX_1 (LINK2B,L2B,'N',CLOSE_IT,CLOSE_STAT,L2B_MSG,'KSF ',NTERM_KFS ,'Y',NDOFS,I_KSF ,J_KSF ,KSF )

            ENDIF

            IF (NTERM_MFS > 0) THEN

               IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN

                                                              ! Allocate and read MSF
                  CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR, AND READ, MSF')
                  CALL LINK_MESSAGE('READ MSF MATRIX')
                  CLOSE_IT   = 'Y'
                  CALL READ_MATRIX_1 ( LINK2S, L2S, 'N', CLOSE_IT, L2SSTAT, L2S_MSG, 'MSF', NTERM_MFS , 'Y', NDOFS,                &
                                       I_MSF , J_MSF , MSF  )
               ENDIF

            ENDIF

            IF (NTERM_QSYS > 0) THEN                          ! Note this will be 0 unless this is STATICS

               CALL LINK_MESSAGE('ALLOCATE ARRAYS FOR, AND READ, QSYS')
               CALL LINK_MESSAGE('READ QSYS MATRIX')
               CLOSE_IT   = 'Y'
               CALL READ_MATRIX_1 ( LINK2C, L2C, 'N', CLOSE_IT, L2CSTAT, L2C_MSG, 'QSYS', NTERM_QSYS, 'Y', NDOFS,                  &
                                    I_QSYS, J_QSYS, QSYS )
               COL_NUM  = 1                                   ! Put QSYS nonzero terms into QSYS_COL.
               NUM_COLS = 1
               IF (NTERM_QSYS > 0) THEN
                  CALL GET_SPARSE_CRS_COL ('QSYS_COL  ', COL_NUM, NTERM_QSYS, NDOFS, NUM_COLS, I_QSYS, J_QSYS, QSYS, ONE,          &
                                            QSYS_COL, NULL_COL)
               ENDIF

            ENDIF

            IF (NTERM_PS > 0) THEN

               CALL LINK_MESSAGE('ALLOCATE SPARSE ARRAYS FOR PS LOADS')
               CALL LINK_MESSAGE('READ PS LOADS')
               CLOSE_IT   = 'N'
               CALL READ_MATRIX_1 ( LINK2D, L2D, 'N', CLOSE_IT, 'KEEP', L2D_MSG, 'PS', NTERM_PS, 'Y', NDOFS,                       &
                                    I_PS, J_PS, PS)
            ENDIF

         ENDIF

      ENDIF

      ! Read MPC constraint matrices
      IF ((ANY_MPCF_OUTPUT > 0) .OR. (ANY_GPFO_OUTPUT > 0) .OR. (WRITE_NEU)) THEN

         IF (NDOFM > 0) THEN

            IF (NTERM_GMN > 0) THEN

                                                           ! Allocate and read GMN and create GMNt
               CALL LINK_MESSAGE('READ GMN MATRIX')
               CLOSE_IT   = 'Y'
               CALL ALLOCATE_SPARSE_MAT ( 'GMN',  NDOFM, NTERM_GMN, SUBR_NAME )
               CALL READ_MATRIX_1 ( LINK2A, L2A, 'N', CLOSE_IT, 'KEEP', L2A_MSG, 'GMN', NTERM_GMN, 'Y', NDOFM                     &
                                  , I_GMN, J_GMN, GMN )
               CALL ALLOCATE_SPARSE_MAT ( 'GMNt', NDOFN, NTERM_GMN, SUBR_NAME )
               CALL MATTRNSP_SS ( NDOFM, NDOFN, NTERM_GMN, 'GMN', I_GMN, J_GMN, GMN, 'GMNt', I_GMNt, J_GMNt, GMNt )

            ENDIF

            IF (NTERM_HMN > 0) THEN                     ! Allocate and read HMN if there are any terms in it.

               CALL LINK_MESSAGE('READ HMN MATRIX')
               CLOSE_IT   = 'Y'
               CALL ALLOCATE_SPARSE_MAT ( 'HMN',  NDOFM, NTERM_HMN, SUBR_NAME )
               CALL READ_MATRIX_1 ( LINK2J, L2J, 'N', CLOSE_IT, L2JSTAT, L2J_MSG, 'HMN', NTERM_HMN, 'Y', NDOFM                     &
                                  , I_HMN, J_HMN, HMN )
            ENDIF

            IF (NTERM_LMN > 0) THEN                     ! Allocate and read LMN if there are any terms in it.

               CALL LINK_MESSAGE('READ LMN MATRIX')
               CLOSE_IT   = 'Y'
               CALL ALLOCATE_SPARSE_MAT ( 'LMN',  NDOFM, NTERM_LMN, SUBR_NAME )
               CALL READ_MATRIX_1 ( LINK2R, L2R, 'N', CLOSE_IT, 'KEEP', L2R_MSG, 'LMN', NTERM_LMN, 'Y', NDOFM                     &
                                  , I_LMN, J_LMN, LMN )
            ENDIF

         ENDIF

      ENDIF

      ! Read MGG mass matrix if this is a dynamics solution and GP force balance is requested
      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
         IF (ANY_GPFO_OUTPUT > 0) THEN
            CALL LINK_MESSAGE('ALLOCATE SPARSE ARRAYS FOR MGG MASS ARRAYS')
            CALL ALLOCATE_SPARSE_MAT ( 'MGG', NDOFG, NTERM_MGG, SUBR_NAME )
            IF (NTERM_MGG > 0) THEN
               CLOSE_IT   = 'Y'
               CALL OURTIM
               CALL READ_MATRIX_1 ( LINK1R, L1R, 'N', CLOSE_IT, L1RSTAT, L1R_MSG, 'MGG', NTERM_MGG, 'Y', NDOFG                     &
                                  , I_MGG, J_MGG, MGG)
            ENDIF
         ENDIF
      ENDIF

      ! Read MLL mass matrix if this is a dynamics solution and GP force balance is requested.
      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
         IF (ANY_GPFO_OUTPUT > 0) THEN
            CALL LINK_MESSAGE('ALLOCATE SPARSE ARRAYS FOR MLL MASS ARRAYS')
            CALL ALLOCATE_SPARSE_MAT ( 'MLL', NDOFL, NTERM_MLL, SUBR_NAME )
            IF (NTERM_MLL > 0) THEN
               CLOSE_IT   = 'Y'
               CALL OURTIM
               CALL READ_MATRIX_1 ( LINK2I, L2I, 'N', CLOSE_IT, L2ISTAT, L2I_MSG, 'MLL', NTERM_MLL, 'Y', NDOFL, I_MLL, J_MLL, MLL)
            ENDIF
         ENDIF
      ENDIF

! Determine if we need to open F25 to write element disp, loads to unformatted file

!xx   ANY_U_P_OUTPUT = IAND(OELDT,IBIT(ELDT_F25_U_P_BIT))
!xx   IF (ANY_U_P_OUTPUT > 0) THEN
!xx      CALL FILE_OPEN ( F25, F25FIL, OUNT, 'REPLACE', F25_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
!xx   ENDIF

      ! Open data files for reading displacements (will be read below in loop over number of subcases/vectors)
      CALL FILE_OPEN ( L5A, LINK5A, OUNT, 'OLD', L5A_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

      ! If this is an eigenvalue problem, determine if there are modes with zero gen stiffness. If so, cannot calc modal masses
      ! or modal participation factors (but only do this if not a CB soln since MPFACTOR and MEFFMASS were calc'd in LINK6 for CB)
      ! EIGEN_VAL was not deallocated in LINK4 (see LINK4 comment 01/11/19) so we do not allocate it here anymore
      ZERO_GEN_STIFF = 'N'
      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
                                                        ! MODE_NUM is not used to det gen stiff but it is read in subr READ_L1M
         CALL ALLOCATE_EIGEN1_MAT ( 'MODE_NUM' , NUM_EIGENS, 1, SUBR_NAME )
!xx      CALL ALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL', NUM_EIGENS, 1, SUBR_NAME )
         CALL ALLOCATE_EIGEN1_MAT ( 'GEN_MASS' , NUM_EIGENS, 1, SUBR_NAME )
         IERROR = 0
         CALL READ_L1M ( IERROR )
         CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM' )
         IF (IERROR /= 0) THEN
            WRITE(ERR,9995) LINKNO,IERROR
            WRITE(F06,9995) LINKNO,IERROR
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         IF (SOL_NAME /= 'GEN CB MODEL') THEN              ! Check for modes with zero generalized stiffness

            DO I=1,NVEC
               IF (DABS(EIGEN_VAL(I)*GEN_MASS(I)) < EPS1) THEN
                  ZERO_GEN_STIFF = 'Y'
                  EXIT
               ENDIF
            ENDDO

            IF (ZERO_GEN_STIFF == 'N') THEN                ! No zero gen stiff, so allocate arrays for eff mass, mpf if requested
               IF (MEFFMASS_CALC == 'Y') THEN
                  CALL ALLOCATE_EIGEN1_MAT ( 'MEFFMASS', NVEC, 6, SUBR_NAME )
               ENDIF
               IF (MPFACTOR_CALC  == 'Y') THEN
                  IF (MPFOUT == '6') THEN
                     CALL ALLOCATE_EIGEN1_MAT ( 'MPFACTOR_N6', NVEC, 6, SUBR_NAME )
                  ELSE
                     CALL ALLOCATE_EIGEN1_MAT ( 'MPFACTOR_NR', NVEC, 6, SUBR_NAME )
                  ENDIF
               ENDIF
            ELSE                                           ! There are 0 gen stiff's, so give msg if eff mass or mpf are requested
               IF ((MEFFMASS_CALC == 'Y') .OR. (MPFACTOR_CALC == 'Y')) THEN
                  WRITE(ERR,9990)
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,9990)
                  ENDIF
               ENDIF
            ENDIF

         ENDIF

      ENDIF

!      CALL WRITE_OP2_GEOM()

      ! Open FEMAP neutral file for writing, if WRITE_NEU, and write FEMAP data block 100
      IF (WRITE_NEU) THEN
         WRITE(CTIME,9000) STIME
         CALL FILE_OPEN ( NEU, NEUFIL, OUNT, 'REPLACE', NEU_MSG, 'WRITE_STIME', 'FORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         FEMAP_BLK = '   100'
         WRITE(NEU,9001)
         WRITE(NEU,9011) FEMAP_BLK
         WRITE(NEU,9012) STIME, F06FIL
         WRITE(NEU,9013) FEMAP_VERSION
         WRITE(NEU,9001)
      ENDIF

      CALL ALLOCATE_COL_VEC ( 'PG_COL', NDOFG, SUBR_NAME )

!-----------------------------------------------------------------------------------------------------------------------------------
      ! Allocate arrays particular to LINK9
      CALL ALLOCATE_LINK9_STUF ( SUBR_NAME )

      ! Initialize JTSUB which will become the col no in the elem thermal loads matrix corresponding to the subcases below.
      JTSUB = 0

      ! Set NUM_SOLNS for use in loop (below) to get outputs for each subcase/solution vector and size. Also, allocate memory for
      ! CB OTM matrices (if CB soln) and open CB OTM output files (OU4(8) for grid related OTM's and OU4(9) for elem related OTM's)
      PROC_PG_OUTPUT = 'Y'
      IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
         NUM_SOLNS = NSUB

      ELSE IF (SOL_NAME(1:8) == 'BUCKLING') THEN

         IF (LOAD_ISTEP == 1) THEN
            NUM_SOLNS = NSUB - NUM_BUCKLING_SUBS  ! all static preload subcases

         ELSE
            NUM_SOLNS = NVEC
            PROC_PG_OUTPUT = 'N'

         ENDIF

      ELSE IF  (SOL_NAME(1:5) == 'MODES') THEN
         NUM_SOLNS = NVEC

      ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
         NUM_SOLNS = NUM_CB_DOFS
         ! These will be allocated with zero size if no output ie requested
         CALL ALLOCATE_CB_GRD_OTM ( 'OTM_ACCE' )
         CALL ALLOCATE_CB_GRD_OTM ( 'OTM_DISP' )
         CALL ALLOCATE_CB_GRD_OTM ( 'OTM_MPCF' )
         CALL ALLOCATE_CB_GRD_OTM ( 'OTM_SPCF' )
         CALL ALLOCATE_CB_ELM_OTM ( 'OTM_ELFE' )
         CALL ALLOCATE_CB_ELM_OTM ( 'OTM_ELFN' )
         CALL ALLOCATE_CB_ELM_OTM ( 'OTM_STRE' )
         CALL ALLOCATE_CB_ELM_OTM ( 'OTM_STRN' )

         ! Get index for file unit nos for elem/grid related OTM unformatted files
         IUE = 0
         IUG = 0
         DO I=1,MOU4
            IF (OU4(I) == OU4_ELM_OTM) IUE = I
            IF (OU4(I) == OU4_GRD_OTM) IUG = I
         ENDDO
         IF ((IUE <= 0) .OR. (IUG <= 0)) THEN              ! Open files for unformatted OTM data
            WRITE(ERR,9901) SUBR_NAME, IUE, IUG, MOU4
            WRITE(F06,9901) SUBR_NAME, IUE, IUG, MOU4
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ELSE
            CALL FILE_OPEN (OU4(IUE),OU4FIL(IUE),OUNT,'REPLACE', OU4_MSG(IUE),'NEITHER','UNFORMATTED','WRITE','REWIND','Y','N')
            CALL FILE_OPEN (OU4(IUG),OU4FIL(IUG),OUNT,'REPLACE', OU4_MSG(IUG),'NEITHER','UNFORMATTED','WRITE','REWIND','Y','N')
         ENDIF

         ! Get index for file unit nos for elem/grid related OTM text files
         ITE = 0
         ITG = 0
         OT4_EROW = 0
         OT4_GROW = 0
         DO I=1,MOT4
            IF (OT4(I) == OT4_ELM_OTM) ITE = I
            IF (OT4(I) == OT4_GRD_OTM) ITG = I
         ENDDO
         IF ((ITE <= 0) .OR. (ITG <= 0)) THEN
            ! Open files for OTM text data to describe the unformatted files above
            WRITE(ERR,9901) SUBR_NAME, ITE, ITG, MOT4
            WRITE(F06,9901) SUBR_NAME, ITE, ITG, MOT4
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ELSE
            CALL FILE_OPEN (OT4(ITE), OT4FIL(ITE), OUNT,'REPLACE', OT4_MSG(ITE),'NEITHER','FORMATTED','WRITE','REWIND','Y','N')
            CALL FILE_OPEN (OT4(ITG), OT4FIL(ITG), OUNT,'REPLACE', OT4_MSG(ITG),'NEITHER','FORMATTED','WRITE','REWIND','Y','N')
         ENDIF

         IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN        ! We need cols of PHIXG to process NDOFR+NVEC cols of GPFO
            CALL FILE_OPEN ( L5B, LINK5B, OUNT, 'OLD', L5B_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
         ENDIF

      ENDIF

! Loop on the number of subcases, or eigenvectors or CB vecs (as the case may be) for all output (except CB accel - processed later)

j_do: DO JVEC=1,NUM_SOLNS
         IF     ((SOL_NAME(1: 7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
            INT_SC_NUM   = JVEC
            FEMAP_SET_ID = SCNUM(JVEC)

         ELSE IF (SOL_NAME(1: 8) == 'BUCKLING') THEN
            IF (LOAD_ISTEP == 2) THEN
               ! Each eigenvector is attributed to its owning buckling subcase via MODE_SUBCASE (same as MODES)
               INT_SC_NUM = 1
               IF (ALLOCATED(MODE_SUBCASE)) THEN
                  IF (JVEC <= SIZE(MODE_SUBCASE)) INT_SC_NUM = MODE_SUBCASE(JVEC)
               ENDIF
               FEMAP_SET_ID = JVEC
            ELSE
               ! Static preload pass: JVEC indexes the static subcases (1..NSUB-NUM_BUCKLING_SUBS)
               INT_SC_NUM   = JVEC
               FEMAP_SET_ID = SCNUM(JVEC)
            ENDIF

         ELSE IF (SOL_NAME(1: 5) == 'MODES') THEN
            ! Each mode is attributed to its owning subcase via MODE_SUBCASE (populated in LINK4). For legacy single-METHOD
            ! decks MODE_SUBCASE is uniformly the canonical subcase, so behaviour matches the original INT_SC_NUM=1 fallback.
            INT_SC_NUM = 1
            IF (ALLOCATED(MODE_SUBCASE)) THEN
               IF (JVEC <= SIZE(MODE_SUBCASE)) INT_SC_NUM = MODE_SUBCASE(JVEC)
            ENDIF
            FEMAP_SET_ID = JVEC

         ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
            INT_SC_NUM   = 1
            FEMAP_SET_ID = JVEC

         ENDIF

         ! Compute per-subcase local eigenvector index for eigen solutions (MODES or BUCKLING step 2).
         ! INT_EIG_NUM is reset to 1 whenever INT_SC_NUM changes (new subcase), and counts up within a subcase.
         ! For non-eigen solutions (STATICS, NLSTATIC, BUCKLING step 1, GEN CB MODEL) INT_EIG_NUM is set to 0.
         IF ((SOL_NAME(1:5) == 'MODES') .OR. ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2))) THEN
            IF (INT_SC_NUM /= PREV_SC_NUM) THEN
               SC_VEC_COUNT = 0
               PREV_SC_NUM  = INT_SC_NUM
            ENDIF
            SC_VEC_COUNT = SC_VEC_COUNT + 1
            INT_EIG_NUM  = SC_VEC_COUNT
         ELSE
            INT_EIG_NUM = 0
         ENDIF

         IF (WRITE_NEU) THEN
            FEMAP_BLK = '   450'
            CALL CONCATENATE_TITLES
            WRITE(NEU,9001)                                ! Write data block 450 to FEMAP NEU file
            WRITE(NEU,9011) FEMAP_BLK
            WRITE(NEU,9022) FEMAP_SET_ID
!           WRITE(NEU,9023) TITLE(JVEC), STITLE(JVEC), LABEL(JVEC)
            WRITE(NEU,9023) TSL
            WRITE(NEU,9024)
            WRITE(NEU,9025)
            WRITE(NEU,9026)
            WRITE(NEU,9001)

            FEMAP_BLK = '   451'                           ! Write header for FEMAP data block 451 (for output vectors)
            WRITE(NEU,9001)
            WRITE(NEU,9011) FEMAP_BLK
         ENDIF


         IF (SOL_NAME(1:8) == 'DIFFEREN') THEN
            JTSUB = 1
            INT_SC_NUM = 1
         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
            JTSUB = 1                                      ! eigenvectors: buckling subcases carry no thermal loads
            ! INT_SC_NUM was already set correctly in the j_do init block above; do not override
         ELSE
            IF (SUBLOD(INT_SC_NUM,2) > 0) THEN                ! JTSUB must only be used in the subrs called if this SUBLOD > 0
               JTSUB = JTSUB + 1
            ENDIF
         ENDIF

         ! Det if GP related outputs were requested in Case Control for this S/C
         SC_ACCE_OUTPUT = IAND(OGROUT(INT_SC_NUM),IBIT(GROUT_ACCE_BIT))
         SC_DISP_OUTPUT = IAND(OGROUT(INT_SC_NUM),IBIT(GROUT_DISP_BIT))
         SC_OLOA_OUTPUT = IAND(OGROUT(INT_SC_NUM),IBIT(GROUT_OLOA_BIT))
         SC_SPCF_OUTPUT = IAND(OGROUT(INT_SC_NUM),IBIT(GROUT_SPCF_BIT))
         SC_MPCF_OUTPUT = IAND(OGROUT(INT_SC_NUM),IBIT(GROUT_MPCF_BIT))
         SC_GPFO_OUTPUT = IAND(OGROUT(INT_SC_NUM),IBIT(GROUT_GPFO_BIT))

                                                           ! Write message to screen
         IF      ((SOL_NAME(1: 7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                            &
                 ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
            CALL LINK_MESSAGE_I('READ G-SET DISPLACEMENTS,                      Subcase', JVEC)

         ELSE IF ((SOL_NAME(1: 5) == 'MODES') .OR. ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2))) THEN
            CALL LINK_MESSAGE_I('READ G-SET EIGENVECTORS,                      Eigenvec', JVEC)

         ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
            CALL LINK_MESSAGE_I('READ G-SET CB VECTORS,                       CB vector', JVEC)

         ENDIF

                                                           ! Read the displ's for the DOF for this subcase/eigenvector
         CALL DEALLOCATE_COL_VEC ( 'UG_COL' )
         CALL ALLOCATE_COL_VEC ( 'UG_COL', NDOFG, SUBR_NAME )
         DO I=1,NDOFG
            READ(L5A,IOSTAT=IOCHK) UGV
            IF (IOCHK /=0) THEN
               REC_NO = I - 1
               CALL READERR (IOCHK, LINK5A, L5A_MSG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            UG_COL(I) = UGV
         ENDDO

         ! If this is a CB soln and JVEC <= NDOFR+NVEC, formulate a col of PHIXG from data in file L5B. Otherwise zero
         IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
            CALL ALLOCATE_COL_VEC ( 'PHIXG_COL', NDOFG, SUBR_NAME )
            IF (JVEC <= NDOFR+NVEC) THEN
               DO I=1,NDOFG
                  READ(L5B,IOSTAT=IOCHK) PHIXGV
                  IF (IOCHK /=0) THEN
                     REC_NO = I - 1
                     CALL READERR (IOCHK, LINK5B, L5B_MSG, REC_NO, OUNT )
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
                  PHIXG_COL(I) = PHIXGV
               ENDDO
            ELSE
               DO I=1,NDOFG
                  PHIXG_COL(I) = ZERO
               ENDDO
            ENDIF
         ENDIF

        ! Process acceleration output requests
        ! 10/01/14: Need BGRID for Femap displs to transform from global to basic
        CALL ALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS', SUBR_NAME )

        NEW_RESULT = .TRUE.
        ITABLE = -1
         IF ((SC_ACCE_OUTPUT > 0) .OR. (WRITE_NEU)) THEN
            IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
               CALL LINK_MESSAGE_I('PROCESS ACCEL OUTPUT REQUESTS,                    "',JVEC)
               CALL OFP1 ( JVEC, 'ACCE', SC_ACCE_OUTPUT, FEMAP_SET_ID, ITG, OT4_GROW, ITABLE, NEW_RESULT )
!              NEW_RESULT = .FALSE.
            ELSE
               WARN_ERR = WARN_ERR + 1
               WRITE(ERR,9453)
               IF (SUPWARN == 'N') THEN
                  WRITE(F06,9453)
               ENDIF
            ENDIF
         ENDIF

         ! Process displacement output requests
         IF ((SC_DISP_OUTPUT > 0) .OR. (WRITE_NEU)) THEN
            CALL LINK_MESSAGE_I('PROCESS DISPL OUTPUT REQUESTS,                    "',JVEC)
            CALL OFP1 ( JVEC, 'DISP', SC_DISP_OUTPUT, FEMAP_SET_ID, ITG, OT4_GROW, ITABLE, NEW_RESULT )
!           NEW_RESULT = .FALSE.
         ENDIF
         !CALL END_OP2_TABLE(ITABLE)

         ! Process applied load (OPG1) output requests
         NEW_RESULT = .TRUE.
         ITABLE = -1
         IF (PROC_PG_OUTPUT == 'Y') THEN
            IF ((SC_OLOA_OUTPUT > 0) .OR. (SC_GPFO_OUTPUT > 0) .OR. (WRITE_NEU)) THEN
               IF  ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'BUCKLING') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
                  CALL LINK_MESSAGE_I('PROCESS APPLIED LOAD OUTPUT REQS,                 "',JVEC)
                  CALL GET_SPARSE_CRS_COL ('PG_COL    ',JVEC      , NTERM_PG, NDOFG, NSUB, I_PG, J_PG, PG, ONE, PG_COL, NULL_COL)
                  CALL OFP1 ( JVEC, 'OLOAD', SC_OLOA_OUTPUT, FEMAP_SET_ID, ITG, OT4_GROW, ITABLE, NEW_RESULT )
!                 NEW_RESULT = .FALSE.
               ENDIF
            ENDIF
         ENDIF
         !CALL END_OP2_TABLE(ITABLE)

        ! Calc SPC forces and process SPC force output requests, if there are any or if GP force balance, modal effective mass and/or
        ! participation factor output is requested. Calc anyway if there are any DOF's in the SA (AUTOSPC) set
        NEW_RESULT = .TRUE.
        ITABLE = -1
         IF (SOL_NAME(1:5) == 'MODES') THEN
            IF (NDOFS == 0) THEN
               IF ((MEFFMASS_CALC == 'Y') .OR. (MPFACTOR_CALC == 'Y')) THEN
                  WRITE(ERR,111)
                  IF (SUPINFO == 'N') THEN
                     WRITE(F06,111)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF ((NDOFS > 0) .OR. (SC_SPCF_OUTPUT > 0) .OR. (SC_GPFO_OUTPUT > 0) .OR.                                                  &
             (MEFFMASS_CALC == 'Y') .OR. (MPFACTOR_CALC == 'Y') .OR. (WRITE_NEU)) THEN

            CALL ALLOCATE_COL_VEC ( 'PS_COL', NDOFS, SUBR_NAME )
            DO K=1,NDOFS
               PS_COL(K) = ZERO
            ENDDO

            IF      ((SOL_NAME(1: 7) == 'STATICS') .OR.  (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                        &
                    ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
               IF (NTERM_PS > 0) THEN
                  CALL GET_SPARSE_CRS_COL ('PS_COL', JVEC      , NTERM_PS, NDOFS, NSUB, I_PS, J_PS, PS, ONE, PS_COL, NULL_COL )
               ENDIF
            ELSE
               DO K=1,NDOFS
                   PS_COL(K) = ZERO
               ENDDO
            ENDIF

            CALL LINK_MESSAGE_I('PROCESS SPC FORCE OUTPUT REQUESTS,                "',JVEC)
            CALL ALLOCATE_COL_VEC ( 'QGs_COL', NDOFG, SUBR_NAME )
           CALL OFP2 ( JVEC, 'SPCF', SC_SPCF_OUTPUT, ZERO_GEN_STIFF, FEMAP_SET_ID, ITG, OT4_GROW, ITABLE, NEW_RESULT )
!           NEW_RESULT = .FALSE.
         ENDIF

         ! Process MPC force output requests, if there are any
         NEW_RESULT = .TRUE.
         IF (NDOFM > 0) THEN

            IF ((SC_MPCF_OUTPUT > 0) .OR. (SC_GPFO_OUTPUT > 0) .OR. (WRITE_NEU)) THEN

               CALL ALLOCATE_COL_VEC ( 'PM_COL', NDOFM, SUBR_NAME )
               DO K=1,NDOFM
                  PM_COL(K) = ZERO
               ENDDO
               IF  ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'BUCKLING') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
                  IF (NTERM_PM > 0) THEN
                     CALL GET_SPARSE_CRS_COL ('PM_COL', JVEC      , NTERM_PM, NDOFM, NSUB, I_PM, J_PM, PM, ONE, PM_COL, NULL_COL )
                  ENDIF
               ENDIF

               CALL LINK_MESSAGE_I('PROCESS MPC FORCE OUTPUT REQUESTS,                "',JVEC)
               CALL ALLOCATE_COL_VEC ( 'QGm_COL', NDOFG, SUBR_NAME )
               CALL OFP2 ( JVEC, 'MPCF', SC_MPCF_OUTPUT, ZERO_GEN_STIFF, FEMAP_SET_ID, ITG, OT4_GROW, ITABLE, NEW_RESULT )
!              NEW_RESULT = .FALSE.
            ENDIF
         ENDIF
         !CALL END_OP2_TABLE(ITABLE)

         ! Process grid point force balance requests
!zzzz    CALL ALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS', SUBR_NAME ) ! 10/01/14: Move to above CALL OFP1. Need for Femap disp
         NEW_RESULT = .TRUE.
         IF (SC_GPFO_OUTPUT > 0) THEN
            CALL ALLOCATE_COL_VEC ( 'FG_COL', NDOFG, SUBR_NAME )
                                                           ! Accel load is Mgg*Ug_ddot = -EIGEN_VAL*Mgg*Ug in eigen analyses.
            IF (SOL_NAME(1:5) == 'MODES') THEN
               CALL MATMULT_SFF ( 'MGG', NDOFG, NDOFG, NTERM_MGG, SYM_MGG, I_MGG, J_MGG, MGG, 'UG', NDOFG, 1, UG_COL, 'Y',         &
                                  'FG', -EIGEN_VAL(JVEC), FG_COL )
                                                           ! DEBUG(191): calc FG_COL as if all inertia force due to MAA*UA_DDOT
               IF ((NDOFO == 0) .AND. (DEBUG(191) == 2)) THEN
                  ! Disabled because it's nothing to do with DEBUG(191) (temperature averaging for thermal loads) and it crashes if used.
                  ! CALL GET_FG_INERTIA_FORCES
               ENDIF
                                                           ! Mult MGG*PHIXG for FG unless JVEC > NDOFR+NVEC, otherwise FG is null
            ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN! Get FG_COL from L5B for CB soln

               IF (JVEC <= NDOFR+NVEC) THEN
                  CALL MATMULT_SFF ( 'MGG', NDOFG, NDOFG, NTERM_MGG, SYM_MGG, I_MGG, J_MGG, MGG, 'PHIXG', NDOFG, 1, PHIXG_COL,     &
                                     'Y', 'FG', ONE, FG_COL )
               ELSE
                  DO I=1,NDOFG
                     FG_COL(I) = ZERO
                  ENDDO
               ENDIF
                                                           ! Get QR_COL, a col from IF_LTM, and put into QGr_COL, (G-set I/F forces)
               CALL ALLOCATE_COL_VEC ( 'QR_COL' , NDOFR, SUBR_NAME )
               CALL ALLOCATE_COL_VEC ( 'QGr_COL', NDOFG, SUBR_NAME )
               CALL GET_SPARSE_CRS_COL ( 'IF_LTM'  , JVEC   , NTERM_IF_LTM, NDOFR, NUM_CB_DOFS, I_IF_LTM, J_IF_LTM, IF_LTM, ONE ,  &
                                          QR_COL , NULL_ROW )
               IF (NULL_ROW == 'Y') THEN
                  DO I=1,NDOFR
                     QR_COL(I) = ZERO
                  ENDDO
               ENDIF
               DO I=1,NDOFG                                ! Calc SPC forces for all grids in requested output set (not only ones
                  CALL TDOF_COL_NUM ( 'R ', R_SET_COL )    ! that have a component in S-set)
                  CALL TDOF_COL_NUM ( 'G ', G_SET_COL )
                  RDOF = TDOF(I,R_SET_COL)
                  GDOF = TDOF(I,G_SET_COL)
                  IF (RDOF > 0) THEN
                     QGr_COL(GDOF) = QR_COL(RDOF)
                  ENDIF
               ENDDO
               CALL DEALLOCATE_COL_VEC ( 'QR_COL' )

            ENDIF

            CALL LINK_MESSAGE_I('PROCESS G.P. FORCE BALANCE REQUESTS,              "',JVEC)
            CALL GP_FORCE_BALANCE_PROC ( JVEC, 'Y' )
            CALL DEALLOCATE_COL_VEC ( 'FG_COL' )
            CALL DEALLOCATE_COL_VEC ( 'QGr_COL' )

         ENDIF

         CALL DEALLOCATE_COL_VEC ( 'PM_COL' )
         CALL DEALLOCATE_COL_VEC ( 'PS_COL' )
         CALL DEALLOCATE_COL_VEC ( 'QGs_COL' )
         CALL DEALLOCATE_COL_VEC ( 'QGm_COL' )
         WRITE(SC1,*) CR13

         ! Process element force/stress output requests
         SC_ELFE_OUTPUT = IAND(OELOUT(INT_SC_NUM),IBIT(ELOUT_ELFE_BIT))
         SC_ELFN_OUTPUT = IAND(OELOUT(INT_SC_NUM),IBIT(ELOUT_ELFN_BIT))
         SC_STRE_OUTPUT = IAND(OELOUT(INT_SC_NUM),IBIT(ELOUT_STRE_BIT))
         SC_STRN_OUTPUT = IAND(OELOUT(INT_SC_NUM),IBIT(ELOUT_STRN_BIT))
         IF((SC_ELFE_OUTPUT > 0) .OR. (SC_ELFN_OUTPUT > 0) .OR. (SC_STRE_OUTPUT > 0) .OR. (SC_STRN_OUTPUT > 0) .OR.                &
            ! (ANY_U_P_OUTPUT > 0) .OR.
            (WRITE_NEU)) THEN
            CALL LINK_MESSAGE_I('PROCESS ELEM FORCE/STRESS REQUESTS,               "',JVEC)
            IF ((DEBUG(176) == 0) .AND. (JVEC == 1)) THEN
               WRITE(ERR,98980)
               WRITE(ERR,98988) DEBUG(176)
               WRITE(ERR,98980)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,98980)
                  WRITE(F06,98988) DEBUG(176)
                  WRITE(F06,98980)
               ENDIF
            ELSE IF (JVEC == 1) THEN
               WRITE(ERR,98980)
               WRITE(ERR,98989) DEBUG(176)
               WRITE(ERR,98980)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,98980)
                  WRITE(F06,98989) DEBUG(176)
                  WRITE(F06,98980)
               ENDIF
            ENDIF
            CALL OFP3 ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )
!           NEW_RESULT = .FALSE.
         ENDIF
         CALL DEALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS' )

         ! Rewind files containing G-set & S-set loads and read STIME so we can read loads file again for next subcase
         INQUIRE ( FILE=LINK1E, EXIST=LEXIST, OPENED=LOPEN )
         IF (LOPEN) THEN
            REWIND (L1E)
            READ(L1E,IOSTAT=IOCHK) XTIME
            MESSAG = 'STIME                   '
            IF (IOCHK /= 0) THEN
               REC_NO = 1
               CALL READERR ( IOCHK, LINK1E, MESSAG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                      ! Can't read STIME from PG loads file
            ENDIF
         ENDIF

         INQUIRE ( FILE=LINK2D, EXIST=LEXIST, OPENED=LOPEN )
         IF (LOPEN) THEN
            REWIND (L2D)
            READ(L2D,IOSTAT=IOCHK) XTIME
            MESSAG = 'STIME                   '
            IF (IOCHK /= 0) THEN
               REC_NO = 1
               CALL READERR ( IOCHK, LINK2D, MESSAG, REC_NO, OUNT )
               CALL OUTA_HERE ( 'Y' )                      ! Can't read STIME from PS loads file
            ENDIF
         ENDIF
                                                           ! For BUCKLING we want to keep UG_COL from the linear statics portion of
         CALL DEALLOCATE_COL_VEC ( 'PHIXG_COL' )

         IF (WRITE_NEU) THEN
            WRITE(NEU,9001)                                ! End of FEMAP block 451 indicator
         ENDIF

      ENDDO j_do

      !IF (POST /= 0) THEN
      !ENDIF
      IF (WRITE_NEU) THEN
         WRITE(NEU,9001)                                   ! End of FEMAP block 451 indicator
         CALL FILE_CLOSE ( NEU, NEUFIL, 'KEEP' )
      ENDIF

      CALL DEALLOCATE_COL_VEC ( 'PG_COL' )

!-----------------------------------------------------------------------------------------------------------------------------------
      ! Write ACCE, DISP, MPCF and SPCF OTM's to OUTPUT4 file. Also write text file descriptions of the rows of the OTM's.
      IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN

         WRITE(F06,*)

         IF ((ANY_ACCE_OUTPUT > 0) .OR. (ANY_DISP_OUTPUT > 0) .OR. (ANY_MPCF_OUTPUT > 0) .OR. (ANY_SPCF_OUTPUT > 0)) THEN
            WRITE(OT4(ITG), 9040) OU4FIL(ITG), NDOFR, NVEC, OU4FIL(ITG)
            OU4STAT(IUG) = 'KEEP'
            OT4STAT(ITG) = 'KEEP'
         ENDIF

         IF (ANY_ACCE_OUTPUT > 0) THEN
            WRITE(OT4(ITG), 9041) NDOFR+NVEC
         ENDIF

         IF (ANY_DISP_OUTPUT > 0) THEN
            WRITE(OT4(ITG), 9042) 2*NDOFR+NVEC
         ENDIF

         IF (ANY_MPCF_OUTPUT > 0) THEN
            WRITE(OT4(ITG), 9043) 2*NDOFR+NVEC
         ENDIF

         IF (ANY_SPCF_OUTPUT > 0) THEN
            WRITE(OT4(ITG), 9044) 2*NDOFR+NVEC
         ENDIF

         IF ((ANY_ACCE_OUTPUT > 0) .OR. (ANY_DISP_OUTPUT > 0) .OR. (ANY_MPCF_OUTPUT > 0) .OR. (ANY_SPCF_OUTPUT > 0)) THEN
            WRITE(OT4(ITG), * )
         ENDIF

         IF (ANY_ACCE_OUTPUT > 0) THEN
            FORM = 2
            CALL WRITE_OU4_FULL_MAT ( 'OTM_ACCE', NROWS_OTM_ACCE, NDOFR+NVEC , FORM, 'N', OTM_ACCE, OU4(IUG) )
            CALL WRITE_GRD_OT4      ( 'OTM_ACCE', NROWS_OTM_ACCE, NROWS_TXT_ACCE, NDOFR+NVEC , TXT_ACCE, OT4(IUG) )
         ENDIF

         IF (ANY_DISP_OUTPUT > 0) THEN
            FORM = 2
            CALL WRITE_OU4_FULL_MAT ( 'OTM_DISP', NROWS_OTM_DISP, NUM_CB_DOFS, FORM, 'N', OTM_DISP, OU4(IUG) )
            CALL WRITE_GRD_OT4      ( 'OTM_DISP', NROWS_OTM_DISP, NROWS_TXT_DISP, NUM_CB_DOFS, TXT_DISP, OT4(IUG) )
         ENDIF

         IF (ANY_MPCF_OUTPUT > 0) THEN
            FORM = 2
            CALL WRITE_OU4_FULL_MAT ( 'OTM_MPCF', NROWS_OTM_MPCF, NUM_CB_DOFS, FORM, 'N', OTM_MPCF, OU4(IUG) )
            CALL WRITE_GRD_OT4      ( 'OTM_MPCF', NROWS_OTM_MPCF, NROWS_TXT_MPCF, NUM_CB_DOFS, TXT_MPCF, OT4(IUG) )
         ENDIF

         IF (ANY_SPCF_OUTPUT > 0) THEN
            FORM = 2
            CALL WRITE_OU4_FULL_MAT ( 'OTM_SPCF', NROWS_OTM_SPCF, NUM_CB_DOFS, FORM, 'N', OTM_SPCF, OU4(IUG) )
            CALL WRITE_GRD_OT4      ( 'OTM_SPCF', NROWS_OTM_SPCF, NROWS_TXT_SPCF, NUM_CB_DOFS, TXT_SPCF, OT4(IUG) )
         ENDIF

      ENDIF

! Write ELFE, ELFN and STRE OTM's to OUTPUT4 file. Also write text file descriptions of the rows of the OTM's.

      IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN

         WRITE(F06,*)

         IF ((ANY_ELFE_OUTPUT > 0) .OR. (ANY_ELFN_OUTPUT > 0) .OR. (ANY_STRE_OUTPUT > 0)) THEN
            WRITE(OT4(ITE), 9050) OU4FIL(ITE), NDOFR, NVEC, OU4FIL(ITE)
            OU4STAT(IUE) = 'KEEP'
            OT4STAT(ITE) = 'KEEP'
         ENDIF

         IF (ANY_ELFE_OUTPUT > 0) THEN
            WRITE(OT4(ITE), 9051) 2*NDOFR+NVEC
         ENDIF

         IF (ANY_ELFN_OUTPUT > 0) THEN
            WRITE(OT4(ITE), 9052) 2*NDOFR+NVEC
         ENDIF

         IF (ANY_STRE_OUTPUT > 0) THEN
            WRITE(OT4(ITE), 9053) 2*NDOFR+NVEC
         ENDIF

         IF (ANY_STRN_OUTPUT > 0) THEN
            WRITE(OT4(ITE), 9054) 2*NDOFR+NVEC
         ENDIF

         WRITE(OT4(ITE),*)
         WRITE(OT4(ITE), 9055)

         IF ((ANY_ELFE_OUTPUT > 0) .OR. (ANY_ELFN_OUTPUT > 0) .OR. (ANY_STRE_OUTPUT > 0)) THEN
            WRITE(OT4(ITE), * )
         ENDIF

         IF (ANY_ELFE_OUTPUT > 0) THEN
            FORM = 2
            CALL WRITE_OU4_FULL_MAT ( 'OTM_ELFE', NROWS_OTM_ELFE, NUM_CB_DOFS, FORM, 'N', OTM_ELFE, OU4(IUE) )
            CALL WRITE_ELM_OT4      ( 'OTM_ELFE', NROWS_OTM_ELFE, NROWS_TXT_ELFE, NUM_CB_DOFS, TXT_ELFE, OT4(IUE) )
         ENDIF

         IF (ANY_ELFN_OUTPUT > 0) THEN
            FORM = 2
            CALL WRITE_OU4_FULL_MAT ( 'OTM_ELFN', NROWS_OTM_ELFN, NUM_CB_DOFS, FORM, 'N', OTM_ELFN, OU4(IUE) )
            CALL WRITE_ELM_OT4      ( 'OTM_ELFN', NROWS_OTM_ELFN, NROWS_TXT_ELFN, NUM_CB_DOFS, TXT_ELFN, OT4(IUE) )
         ENDIF

         IF (ANY_STRE_OUTPUT > 0) THEN
            FORM = 2
            CALL WRITE_OU4_FULL_MAT ( 'OTM_STRE', NROWS_OTM_STRE, NUM_CB_DOFS, FORM, 'N', OTM_STRE, OU4(IUE) )
            CALL WRITE_ELM_OT4      ( 'OTM_STRE', NROWS_OTM_STRE, NROWS_TXT_STRE, NUM_CB_DOFS, TXT_STRE, OT4(IUE) )
         ENDIF

         IF (ANY_STRN_OUTPUT > 0) THEN
            FORM = 2
            CALL WRITE_OU4_FULL_MAT ( 'OTM_STRN', NROWS_OTM_STRN, NUM_CB_DOFS, FORM, 'N', OTM_STRN, OU4(IUE) )
            CALL WRITE_ELM_OT4      ( 'OTM_STRN', NROWS_OTM_STRN, NROWS_TXT_STRN, NUM_CB_DOFS, TXT_STRN, OT4(IUE) )
         ENDIF

      ENDIF

      ! Close OTM text files (unformatted OU4 files closed in subr CLOSE_LIJFILES)
      IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
         DO I=1,MOT4
            CALL FILE_CLOSE ( OT4(I), OT4FIL(I), OT4STAT(I) )
         ENDDO
      ENDIF

!-----------------------------------------------------------------------------------------------------------------------------------
! If sol is eigens (not CB) then MPFACTOR, MEFFMASS were calc'd in OFP2

      IF ((MPFACTOR_CALC  == 'Y') .AND. (ZERO_GEN_STIFF == 'N')) THEN
         CALL WRITE_MPFACTOR
      ENDIF

      IF ((MEFFMASS_CALC == 'Y') .AND. (ZERO_GEN_STIFF == 'N')) THEN
         IF (DABS(WTMASS) < EPS1) THEN
            WRITE(ERR,9991)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,9991)
            ENDIF
         ENDIF
         CALL WRITE_MEFFMASS
      ENDIF

      ! Call OUTPUT4 processor to process output requests for OUTPUT4 matrices generated in this link
      IF (NUM_OU4_REQUESTS > 0) THEN
         CALL LINK_MESSAGE('WRITE OUTPUT4 MATRICES      ')
         WRITE(F06,*)
         CALL OUTPUT4_PROC ( SUBR_NAME )
      ENDIF

      ! Deallocate MPFACTOR, MEFFMASS
      CALL DEALLOCATE_EIGEN1_MAT ( 'MPFACTOR_N6' )
      CALL DEALLOCATE_EIGEN1_MAT ( 'MPFACTOR_NR' )
      CALL DEALLOCATE_EIGEN1_MAT ( 'MEFFMASS' )

      IF ((DEBUG(195) > 0) .AND. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
         CALL WRITE_OTM_TO_F06
      ENDIF

      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
         CALL DEALLOCATE_EIGEN1_MAT ( 'EIGEN_VAL' )
         CALL DEALLOCATE_EIGEN1_MAT ( 'GEN_MASS' )
         CALL DEALLOCATE_EIGEN1_MAT ( 'MODE_NUM' )
      ENDIF

! Deallocate some arrays

!xx   WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KSF ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KSF' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate KSFD', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KSFD')
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MGG ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MGG' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PG  ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'PG' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PM  ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'PM' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate PS  ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'PS' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate QSYS', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'QSYS' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GMN ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'GMN' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate GMNt', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'GMNt' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate HMN ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'HMN' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MSF ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MSF' )

      ! save MLL from deallocation in case we need to use it for eigenvalue
      ! estimation in a next step of the eigen solution
      IF ((SOL_NAME(1:8) /= 'BUCKLING') .OR. (LOAD_ISTEP == 2)) THEN
         WRITE(SC1,12345,ADVANCE='NO') '       Deallocate MLL ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MLL' )
      END IF

      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate LMN ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'LMN' )
      WRITE(SC1,12345,ADVANCE='NO') '       Deallocate QSYS', CR13  ;   CALL DEALLOCATE_COL_VEC    ( 'QSYS_COL' )

      CALL DEALLOCATE_IN4_FILES  ( 'IN4FIL' )
                                                           ! Deallocate data in file LINK1D
      IF ((SOL_NAME(1:8) /= 'BUCKLING') .OR. (LOAD_ISTEP == 2)) THEN
         ! gotta make SCNUM survive past the 1st run because we use it in LINK4
         CALL DEALLOCATE_MODEL_STUF ( 'SCNUM' )
         ! titles too, for use in block labels during link4/5
         CALL DEALLOCATE_MODEL_STUF ( 'TITLES' )
      END IF
      CALL DEALLOCATE_MODEL_STUF ( 'GROUT, ELOUT' )
                                                           ! Deallocate data in file LINK1G (except ETYPE, EDAT, EPNT
      CALL DEALLOCATE_MODEL_STUF ( 'ESORT1' )
      CALL DEALLOCATE_MODEL_STUF ( 'ESORT2' )
                                                           ! Deallocate data in file LINK1K      ! Deallocate data in file LINK1Q
      CALL DEALLOCATE_MODEL_STUF ( 'PPNT, PDATA, PTYPE' )
      CALL DEALLOCATE_MODEL_STUF ( 'PLOAD4_3D_DATA' )

      CALL DEALLOCATE_MODEL_STUF ( 'GROUT, ELOUT' )

      IF (SOL_NAME(1:12) /= 'GEN CB MODEL') THEN
         CALL DEALLOCATE_IN4_FILES ( 'IN4FIL' )
      ENDIF

      CALL DEALLOCATE_LINK9_STUF

      IF ((SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC')) THEN
         CALL DEALLOCATE_MODEL_STUF ( 'ETYPE, EDAT, EPNT' )
         CALL DEALLOCATE_MODEL_STUF ( 'VVEC, OFFSETS, PLATE stuff' )
         CALL DEALLOCATE_MODEL_STUF ( 'ELEM PROPERTIES AND MATERIALS' )
         CALL DEALLOCATE_MODEL_STUF ( 'EOFF' )
         CALL DEALLOCATE_MODEL_STUF ( 'ELDT' )
         CALL DEALLOCATE_MODEL_STUF ( 'SUBLOD' )
         CALL DEALLOCATE_MODEL_STUF ( 'TPNT, TDATA' )
         CALL DEALLOCATE_MODEL_STUF ( 'GTEMP' )
         CALL DEALLOCATE_MODEL_STUF ( 'GRID_ELEM_CONN_ARRAY' )
         CALL DEALLOCATE_MODEL_STUF ( 'CORD, RCORD' )
         CALL DEALLOCATE_MODEL_STUF ( 'GRID_SNORM' )
!xx      CALL DEALLOCATE_MODEL_STUF ( 'GRID_ID' )
!xx      CALL DEALLOCATE_MODEL_STUF ( 'GRID, RGRID' )
!xx      CALL DEALLOCATE_MODEL_STUF ( 'ELEM PROPERTIES AND MATERIALS' )
!xx      CALL DEALLOCATE_MODEL_STUF ( 'EOFF' )
!xx      CALL DEALLOCATE_MODEL_STUF ( 'ELDT' )
         CALL DEALLOCATE_DOF_TABLES ( 'TSET' )             ! NB *** 10/05/21. Was commented out since ver 6.21 or 6.22
         CALL DEALLOCATE_DOF_TABLES ( 'TDOF' )             ! NB *** 10/05/21. Was commented out since ver 6.21 or 6.22
         CALL DEALLOCATE_DOF_TABLES ( 'TDOF_ROW_START' )   ! NB *** 10/05/21. Was commented out since ver 6.21 or 6.22
         CALL DEALLOCATE_DOF_TABLES ( 'TDOFI' )            ! NB *** 10/05/21. Was commented out since ver 6.21 or 6.22
!xx      CALL DEALLOCATE_MODEL_STUF ( 'GRID_SEQ, INV_GRID_SEQ' )
!xx      CALL DEALLOCATE_MODEL_STUF ( 'SUBLOD' )
!xx      CALL DEALLOCATE_MODEL_STUF ( 'VVEC, OFFSETS, PLATE stuff' )
      ENDIF

!-----------------------------------------------------------------------------------------------------------------------------------
! If there were matrices that were requested to be partitioned but were not (i.e. maybe can't be partitioned in this SOL_NAME), then
! write messages

      NUM_OU4_NOT_PART = 0
      DO I=1,NUM_OU4_REQUESTS
         IF ((OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') .AND. (HAS_OU4_MAT_BEEN_PROCESSED(I,1) == 'N')) THEN
            NUM_OU4_NOT_PART = NUM_OU4_NOT_PART + 1
         ENDIF
      ENDDO
      IF (NUM_OU4_NOT_PART > 0) THEN
         WRITE(ERR,101) NUM_OU4_NOT_PART, SOL_NAME
         WRITE(F06,101) NUM_OU4_NOT_PART, SOL_NAME
         K = 0
         DO I=1,NUM_OU4_REQUESTS
            IF ((OU4_PART_MAT_NAMES(I,1)(1:) /= ' ') .AND. (HAS_OU4_MAT_BEEN_PROCESSED(I,1) == 'N')) THEN
               K = K + 1
               WRITE(ERR,102) K, OU4_PART_MAT_NAMES(I,1)
               WRITE(F06,102) K, OU4_PART_MAT_NAMES(I,1)
            ENDIF
         ENDDO
      ENDIF
  101 FORMAT(' *INFORMATION: THE FOLLOWING ',I4,' MATRICES WERE REQUESTED TO BE PARTITIONED BUT WERE NOT AVAIL IN THIS SOL = ', A)

  102 FORMAT('               (',I2,')',1X,A)

! Process is now complete so set COMM(LINKNO)

      COMM(LINKNO) = 'C'

! Write data to L1A

      CALL WRITE_L1A ( L1ASTAT, 'Y' )

! Do file inquire, if requested

      IF (( DEBUG(193) == 9) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK9' )
      ENDIF

! Close F25

      INQUIRE ( FILE=F25FIL, EXIST=LEXIST, OPENED=LOPEN )
      IF (LOPEN) THEN
         CALL FILE_CLOSE ( F25, F25FIL, 'KEEP' )
      ELSE
         CALL FILE_CLOSE ( F25, F25FIL, 'DELETE' )
      ENDIF

      ! Check allocation status of allocatable arrays, if requested
      IF (DEBUG(100) > 0) THEN
         CALL CHK_ARRAY_ALLOC_STAT
         IF (DEBUG(100) > 1) THEN
            CALL WRITE_ALLOC_MEM_TABLE ( 'at the end of '//SUBR_NAME )
         ENDIF
      ENDIF

      ! Write LINK9 end to F06
      CALL OURTIM
      WRITE(F06,151) LINKNO

      ! Leave the closing of BUG, ERR, F06 files until after LINK9 returns to MYSTRAN.for

      ! Close some files
      IF ((SOL_NAME(1:8) == 'BUCKLING') .OR. (SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
         CALL FILE_CLOSE ( L1E, LINK1E, 'KEEP' )
      ELSE
         CALL FILE_CLOSE ( L1E, LINK1E, L1ESTAT )
      ENDIF

! Write LINK9 end to screen

      CALL OURTIM
      WRITE(SC1,153) LINKNO

! **********************************************************************************************************************************
  111 FORMAT(' *INFORMATION: CASE CONTROL REQUEST WAS MADE FOR MPFACTOR OR MEFFMASS BUT THESE CANNOT BE CALCULATED IN THIS SOL'    &
                    ,/,14X,' IF THERE ARE NO SINGLE POINT CONSTRAINTS TO GROUND THE STRUCTURE')

  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END')

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  153 FORMAT(  ' >> LINK',I3,' END')

 9000 FORMAT(I14)

 9001 FORMAT('   -1')

 9003 FORMAT('-------------------------------------------------------------------------------------------------------------------',&
             '-----------------')

 9011 FORMAT(A)

 9012 FORMAT(I14,', ',A,',')

 9013 FORMAT(F8.1,',')

 9022 FORMAT(I8,',')

 9023 FORMAT(A,1X,A,1X,A,',')

 9024 FORMAT('0,1,')

 9025 FORMAT('0.,')

 9026 FORMAT('0,')

 9040 FORMAT('This text file describes the rows of the grid related OTM matrices written to unformatted file: ',A,/,               &
             '-------------------------------------------------------------------------------------------------------------------',&
             '----------------'//,                                                                                                 &
             'The description for each of the matrices has the headers:',/,                                                        &
             '       ROW        : row number in the individual OTM described',/,                                                   &
             '       DESCRIPTION: what OTM is this',/,                                                                             &
             '       GRID       : grid number for this row of the OTM',/,                                                          &
             '       COMP       : displacement component number (1,2,3 translations and 4,5,6 rotations)',//,                      &
             'The number of rows for each OTM depends on the output requests, by the user, in Case Control',//,                    &
             'The number of cols for each OTM depends on the number of support DOFs (NDOFR) and the number of eigenvecors (NVEC)', &
             'where:',/,'       NDOFR = ',I8,/,'       NVEC  = ',I8,//,                                                            &
             'This text file has descriptions for the following grid relatad OTMs from ',A)

 9041 FORMAT('       Acceleration OTM (matrix OTM_ACCE) with   NDOFR + NVEC = ',I8,' cols')

 9042 FORMAT('       Displacement OTM (matrix OTM_DISP) with 2*NDOFR + NVEC = ',I8,' cols')

 9043 FORMAT('       MPC force    OTM (matrix OTM_MPCF) with 2*NDOFR + NVEC = ',I8,' cols')

 9044 FORMAT('       SPC force    OTM (matrix OTM_SPCF) with 2*NDOFR + NVEC = ',I8,' cols')

 9050 FORMAT('This text file describes the rows of the elem related OTM matrices written to unformatted file: ',A,/,               &
             '-------------------------------------------------------------------------------------------------------------------',&
             '----------------'//,                                                                                                 &
             'The description for each of the matrices has the headers:',/,                                                        &
             '       ROW        : row number in the individual OTM described',/,                                                   &
             '       DESCRIPTION: what OTM is this',/,                                                                             &
             '       TYPE       : element type',/,                                                                                 &
             '       EID        : element ID',/,                                                                                   &
             'Then, for the element nodal force OTM:',/,                                                                           &
             '       GRID       : grid number of the element that the OTM is for',/,                                               &
             '       COMP       : displacement component number (1,2,3 translations and 4,5,6 rotations)',/,                       &
             'and for element engineering force and element stress OTMs:',/,                                                       &
             '       ITEM       : element force or stress item (axial force, torque, etc)'//,                                      &
             'The number of rows for each OTM depends on the output requests, by the user, in Case Control',//,                    &
             'The number of cols for each OTM depends on the number of support DOFs (NDOFR) and the number of eigenvecors (NVEC)', &
             'where:',/,'       NDOFR = ',I8,/,'       NVEC  = ',I8,//,                                                            &
             'This text file has descriptions for the following element related OTMs from ',A)

 9051 FORMAT('       Element engr  force OTM (matrix OTM_ELFE) with 2*NDOFR + NVEC = ',I8,' cols')

 9052 FORMAT('       Element nodal force OTM (matrix OTM_ELFN) with 2*NDOFR + NVEC = ',I8,' cols')

 9053 FORMAT('       Element stress      OTM (matrix OTM_STRE) with 2*NDOFR + NVEC = ',I8,' cols')

 9054 FORMAT('       Element strain      OTM (matrix OTM_STRN) with 2*NDOFR + NVEC = ',I8,' cols')

 9055 FORMAT('The heading "LOCATION" for stresses and strains only has significance for the elements that allow output of these',/,&
             'quantities at specific locations as specified on the Case Control STRESS, STRAIN entries (see MYSTRAN Users Manual)')

 9095 FORMAT(1X,'********** Subcase No. ',I8,' **********')

 9096 FORMAT(1X,'********** Eigenvector No. ',I8,' **********')

 9097 FORMAT(1X,'********** CB vector No. ',I8,' **********')

 9453 FORMAT(' *WARNING    : ACCELERATION REQUESTS ONLY PROGRAMMED IN CB MODEL SOLUTION')

 9990 FORMAT(' *INFORMATION: CANNOT CALCULATE MODAL EFFECTIVE MASS OR MODAL PARTICIPATION FACTORS SINCE THERE ARE MODES WITH ZERO '&
                           ,'GENERALIZED STIFFNESS')

 9901 FORMAT(' *ERROR  9901: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VARIABLE IE = ',I8,' OR IG = ',I8,' WRONG. MUST BE IN RANGE 1 TO MOU4 = ',I8)

 9991 FORMAT(' *INFORMATION: CANNOT CONVERT UNITS OF MODAL EFFECTIVE MASS USING PARAM WTMASS SINCE IT IS ZERO')

 9995 FORMAT(/,' PROCESSING ENDED IN LINK ',I3,' DUE TO ABOVE ',I8,' ERRORS')

 9998 FORMAT(' *ERROR  9998: COMM ',I3,' INDICATES UNSUCCESSFUL LINK ',I2,' COMPLETION.'                                           &
                    ,/,14X,' FATAL ERROR - CANNOT START LINK ',I2)




98980 FORMAT(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

98988 FORMAT(' *INFORMATION: Due to DEBUG(176) = ',i3                                                                              &
                    ,/,14x,' Plate elem engr forces and stresses will be calculated by multiplying strains by the material matrix' &
                    ,/,14x,' Strains are calculated using the strain-displ matrices:'                                              &
                    ,/,14x,' BE1 (membrane), BE2 (bending), BE3 (transverse shear) times displacements')

98989 FORMAT(' *INFORMATION: Due to DEBUG(176) = ',i3                                                                              &
                    ,/,14x,' Plate elem engr forces and stresses will be calculated by multiplying stress-displ matrices:'         &
                    ,/,14x,' SE1 (membrane), SE2 (bending), SE3 (transv shear) times displacements')

12345 FORMAT(A,10X,A)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE CONCATENATE_TITLES

! Concatenate TITLE, STITLE, LABEL for FEMAP block 450 in FEMAP NEU file

      USE PENTIUM_II_KIND

      IMPLICIT NONE

      INTEGER(LONG)                        :: P1
      INTEGER(LONG)                        :: P2
      INTEGER(LONG)                        :: P3
      INTEGER(LONG)                        :: P4
      INTEGER(LONG)                        :: P5
      INTEGER(LONG)                        :: TITLE_LEN  = 0
      INTEGER(LONG)                        :: STITLE_LEN = 0
      INTEGER(LONG)                        :: LABEL_LEN  = 0

! **********************************************************************************************************************************
      DO I=CC_ENTRY_LEN,1,-1
         IF (TITLE(INT_SC_NUM)(I:I) /= ' ') THEN
            TITLE_LEN = I
            EXIT
         ENDIF
      ENDDO

      DO I=CC_ENTRY_LEN,1,-1
         IF (STITLE(INT_SC_NUM)(I:I) /= ' ') THEN
            STITLE_LEN = I
            EXIT
         ENDIF
      ENDDO

      DO I=CC_ENTRY_LEN,1,-1
         IF (LABEL(INT_SC_NUM)(I:I) /= ' ') THEN
            LABEL_LEN = I
            EXIT
         ENDIF
      ENDDO

      P1 = TITLE_LEN
      P2 = P1 + 2
      P3 = P2 + STITLE_LEN
      P4 = P3 + 2
      P5 = P4 + LABEL_LEN

      TSL(   1:P1) = TITLE(INT_SC_NUM)
      TSL(P1+1:P2) = '. '
      TSL(P2+1:P3) = STITLE(INT_SC_NUM)
      TSL(P3+1:P4) = '. '
      TSL(P4+1:P5) = LABEL(INT_SC_NUM)
      TSL(P5+1:  ) = ','

      END SUBROUTINE CONCATENATE_TITLES

! ##################################################################################################################################

      SUBROUTINE WRITE_OTM_TO_F06

      IMPLICIT NONE

! **********************************************************************************************************************************
! Write DISP, MPCF and SPCF OTM's to F06 file.

      IF ((SC_ACCE_OUTPUT > 0) .OR. (SC_DISP_OUTPUT > 0) .OR. (SC_MPCF_OUTPUT > 0) .OR. (SC_SPCF_OUTPUT > 0) .OR.                  &
          (SC_ELFE_OUTPUT > 0) .OR. (SC_ELFN_OUTPUT > 0) .OR. (SC_STRE_OUTPUT > 0)) THEN
         WRITE(F06,99770)
      ENDIF

      IF (SC_ACCE_OUTPUT > 0) THEN
         WRITE(F06,*) ' Output transformation matrix for acceleration: OTM_ACCE'
         WRITE(F06,*) ' -------------------------------------------------------'
         WRITE(F06,99771) (I,I=1,NDOFR+NVEC)
         DO I=1,NROWS_OTM_ACCE
            WRITE(F06,99778) I, (OTM_ACCE(I,J),J=1,NDOFR+NVEC)
         ENDDO
         WRITE(F06,*)
      ENDIF

      IF (SC_DISP_OUTPUT > 0) THEN
         WRITE(F06,*) ' Output transformation matrix for displacement: OTM_DISP'
         WRITE(F06,*) ' -------------------------------------------------------'
         WRITE(F06,99771) (I,I=1,NUM_CB_DOFS)
         DO I=1,NROWS_OTM_DISP
            WRITE(F06,99778) I, (OTM_DISP(I,J),J=1,NUM_CB_DOFS)
         ENDDO
         WRITE(F06,*)
      ENDIF

      IF (SC_MPCF_OUTPUT > 0) THEN
         WRITE(F06,*) ' Output transformation matrix for MPC forces: OTM_MPCF'
         WRITE(F06,*) ' -----------------------------------------------------'
         WRITE(F06,99771) (I,I=1,NUM_CB_DOFS)
         DO I=1,NROWS_OTM_MPCF
            WRITE(F06,99778) I, (OTM_MPCF(I,J),J=1,NUM_CB_DOFS)
         ENDDO
         WRITE(F06,*)
      ENDIF

      IF (SC_SPCF_OUTPUT > 0) THEN
         WRITE(F06,*) ' Output transformation matrix for SPC forces: OTM_SPCF'
         WRITE(F06,*) ' -----------------------------------------------------'
         WRITE(F06,99771) (I,I=1,NUM_CB_DOFS)
         DO I=1,NROWS_OTM_SPCF
            WRITE(F06,99778) I, (OTM_SPCF(I,J),J=1,NUM_CB_DOFS)
         ENDDO
         WRITE(F06,*)
      ENDIF

      IF (SC_ELFE_OUTPUT > 0) THEN
         WRITE(F06,*) ' Output transformation matrix for element engineering forces: OTM_ELFE'
         WRITE(F06,*) ' ---------------------------------------------------------------------'
         WRITE(F06,99771) (I,I=1,NUM_CB_DOFS)
         DO I=1,NROWS_OTM_ELFE
            WRITE(F06,99778) I, (OTM_ELFE(I,J),J=1,NUM_CB_DOFS)
         ENDDO
         WRITE(F06,*)
      ENDIF

      IF (SC_ELFN_OUTPUT > 0) THEN
         WRITE(F06,*) ' Output transformation matrix for element nodal forces: OTM_ELFN'
         WRITE(F06,*) ' ---------------------------------------------------------------'
         WRITE(F06,99771) (I,I=1,NUM_CB_DOFS)
         DO I=1,NROWS_OTM_ELFN
            WRITE(F06,99778) I, (OTM_ELFN(I,J),J=1,NUM_CB_DOFS)
         ENDDO
         WRITE(F06,*)
      ENDIF

      IF (SC_STRE_OUTPUT > 0) THEN
         WRITE(F06,*) ' Output transformation matrix for element stresses: OTM_STRE'
         WRITE(F06,*) ' -----------------------------------------------------------'
         WRITE(F06,99771) (I,I=1,NUM_CB_DOFS)
         DO I=1,NROWS_OTM_STRE
            WRITE(F06,99778) I, (OTM_STRE(I,J),J=1,NUM_CB_DOFS)
         ENDDO
         WRITE(F06,*)
      ENDIF

      RETURN

! **********************************************************************************************************************************
99770 FORMAT('*******************************************************************************************************************',&
             '****************',/,' DEBUG(195) output:',/,' -----------------',/)


99771 FORMAT(3X,32767(9X,I5))

99778 format(i8,32767(1es14.6))

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_OTM_TO_F06

! ##################################################################################################################################

      SUBROUTINE GET_FG_INERTIA_FORCES

      USE PENTIUM_II_KIND
      USE IOUNT1, ONLY                :  ERR, F06, LINK2I, L2I, L2I_MSG, L2ISTAT
      USE SCONTR, ONLY                :  NDOFA, NDOFF, NDOFG, NDOFL, NDOFM, NDOFN, NDOFO, NDOFS, NDOFR, NTERM_MLL
      USE SPARSE_MATRICES, ONLY       :  I_MLL, J_MLL, MLL, SYM_MLL
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL
      USE DOF_TABLES, ONLY            :  TDOFI
      USE COL_VECS, ONLY              :  FA_COL, FF_COL, FG_COL, FL_COL, FM_COL, FN_COL, FO_COL, FR_COL, FS_COL,                   &
                                         PHIXG_COL, PHIXL_COL


      IMPLICIT NONE

      INTEGER(LONG)                   :: K                 ! Counter
      INTEGER(LONG)                   :: A_SET_COL         ! Col no. in TDOF for A_SET displ set definition
      INTEGER(LONG)                   :: F_SET_COL         ! Col no. in TDOF for F_SET displ set definition
      INTEGER(LONG)                   :: G_SET_COL         ! Col no. in TDOF for G_SET displ set definition
      INTEGER(LONG)                   :: L_SET_COL         ! Col no. in TDOF for L_SET displ set definition
      INTEGER(LONG)                   :: M_SET_COL         ! Col no. in TDOF for M_SET displ set definition
      INTEGER(LONG)                   :: N_SET_COL         ! Col no. in TDOF for N_SET displ set definition
      INTEGER(LONG)                   :: O_SET_COL         ! Col no. in TDOF for O_SET displ set definition
      INTEGER(LONG)                   :: R_SET_COL         ! Col no. in TDOF for R_SET displ set definition
      INTEGER(LONG)                   :: S_SET_COL         ! Col no. in TDOF for S_SET displ set definition


! **********************************************************************************************************************************
! Partition PHIXG_COL to the L-set PHIXL_COL

      CALL ALLOCATE_COL_VEC ('PHIXL_COL', NDOFL, SUBR_NAME)
      CALL TDOF_COL_NUM ( 'L ', L_SET_COL )
      K = 0
      DO I=1,NDOFG
         IF (TDOFI(I,L_SET_COL) > 0) THEN
            K = K + 1
            PHIXL_COL(K) = PHIXG_COL(I)
         ENDIF
      ENDDO

! Multiply MLL times PHIXL_COL to get FL_COL

      CALL ALLOCATE_COL_VEC ('FL_COL', NDOFL, SUBR_NAME)
      CALL MATMULT_SFF ( 'MLL', NDOFL, NDOFL, NTERM_MLL, SYM_MLL, I_MLL, J_MLL, MLL, 'UL', NDOFL, 1, PHIXL_COL, 'Y',               &
                         'FL', EIGEN_VAL(JVEC), FL_COL )
      CALL DEALLOCATE_COL_VEC ( 'PHIXL_COL' )

! Build FL_COL to FG_COL by adding zeros
! (1) Build FA from FL and null FR

      CALL ALLOCATE_COL_VEC ( 'FA_COL', NDOFA, SUBR_NAME )
      CALL ALLOCATE_COL_VEC ( 'FR_COL', NDOFR, SUBR_NAME )
      IF (NDOFR > 0) THEN

         CALL TDOF_COL_NUM('A ', A_SET_COL)
         CALL TDOF_COL_NUM('L ', L_SET_COL)
         CALL TDOF_COL_NUM('R ', R_SET_COL)

         DO I=1,NDOFR
            FR_COL(I) = ZERO
         ENDDO

         CALL MERGE_COL_VECS ( L_SET_COL, NDOFL, FL_COL, R_SET_COL, NDOFR, FR_COL, A_SET_COL, NDOFA, FA_COL )

      ELSE

         DO I=1,NDOFA
            FA_COL(I) = FL_COL(I)
         ENDDO

      ENDIF

      CALL DEALLOCATE_COL_VEC ( 'FL_COL' )
      CALL DEALLOCATE_COL_VEC ( 'FR_COL' )

! (2) Build FF from FA and null FO

      CALL ALLOCATE_COL_VEC ( 'FF_COL' , NDOFF, SUBR_NAME )
      CALL ALLOCATE_COL_VEC ( 'FO_COL' , NDOFO, SUBR_NAME )
      IF (NDOFO > 0) THEN

         CALL TDOF_COL_NUM('A ', A_SET_COL)
         CALL TDOF_COL_NUM('O ', O_SET_COL)
         CALL TDOF_COL_NUM('F ', F_SET_COL)

         DO I=1,NDOFO
            FO_COL(I) = ZERO
         ENDDO

         CALL MERGE_COL_VECS ( A_SET_COL, NDOFA, FA_COL, O_SET_COL, NDOFO, FO_COL, F_SET_COL, NDOFF, FF_COL )

      ELSE

         DO I=1,NDOFF
            FF_COL(I) = FA_COL(I)
         ENDDO

      ENDIF

      CALL DEALLOCATE_COL_VEC ( 'FA_COL' )
      CALL DEALLOCATE_COL_VEC ( 'FO_COL' )

! (3) Build FN from FF and null FS

      CALL ALLOCATE_COL_VEC ( 'FN_COL' , NDOFN, SUBR_NAME )
      CALL ALLOCATE_COL_VEC ( 'FS_COL' , NDOFS, SUBR_NAME )
      IF (NDOFS > 0) THEN

         CALL TDOF_COL_NUM('N ', N_SET_COL)
         CALL TDOF_COL_NUM('F ', F_SET_COL)
         CALL TDOF_COL_NUM('S ', S_SET_COL)

         DO I=1,NDOFS
            FS_COL(I) = ZERO
         ENDDO

         CALL MERGE_COL_VECS ( F_SET_COL, NDOFF, FF_COL, S_SET_COL, NDOFS, FS_COL, N_SET_COL, NDOFN, FN_COL )

      ELSE

         DO I=1,NDOFN
            FN_COL(I) = FF_COL(I)
         ENDDO

      ENDIF

      CALL DEALLOCATE_COL_VEC ( 'FF_COL' )
      CALL DEALLOCATE_COL_VEC ( 'FS_COL' )

! (4) Build FG from FN and null FM (FG_COL already allocated in LINK9)

      CALL ALLOCATE_COL_VEC ( 'FM_COL' , NDOFM, SUBR_NAME )
      IF (NDOFM > 0) THEN

         CALL TDOF_COL_NUM('G ', G_SET_COL)
         CALL TDOF_COL_NUM('N ', N_SET_COL)
         CALL TDOF_COL_NUM('M ', M_SET_COL)

         DO I=1,NDOFM
            FM_COL(I) = ZERO
         ENDDO

         CALL MERGE_COL_VECS ( N_SET_COL, NDOFN, FN_COL, M_SET_COL, NDOFM, FM_COL, G_SET_COL, NDOFG, FG_COL )

      ELSE

         DO I=1,NDOFG
            FG_COL(I) = FN_COL(I)
         ENDDO

      ENDIF

      CALL DEALLOCATE_COL_VEC ( 'FN_COL' )
      CALL DEALLOCATE_COL_VEC ( 'FM_COL' )

      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE GET_FG_INERTIA_FORCES

      END SUBROUTINE LINK9


      SUBROUTINE LINK9S

! Reads data from files L1D, L1G, L1K and L1Q (created in LINK1) needed in LINK9

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_ERR, ERR, F06,                                         &
                                         L1D    , L1G    , L1K    , L1Q    ,                                                       &
                                         LINK1D , LINK1G , LINK1K , LINK1Q ,                                                       &
                                         L1D_MSG, L1G_MSG, L1K_MSG, L1Q_MSG,                                                       &
                                         L1DSTAT, L1GSTAT, L1KSTAT, L1QSTAT

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, DATA_NAM_LEN, MPCOMP0, MRPCOMP0, MPCOMP_PLIES,                              &
                                         MRPCOMP_PLIES, MRMATLC, MPBAR, MRPBAR, MPBEAM, MRPBEAM, MPBUSH, MRPBUSH, MPELAS, MRPELAS, &
                                         MPLOAD4_3D_DATA, MPROD, MRPROD, MPSHEAR, MRPSHEAR, MPSHEL, MRPSHEL, MPUSER1, MRPUSER1,    &
                                         MPUSERIN, MUSERIN_MAT_NAMES, MMATL, MPSOLID, NEDAT, NBAROFF, NBUSHOFF, NELE, NGRID,       &
                                         NMATANGLE, NMATL, NPBAR, NPBEAM, NPBUSH, NPCOMP, NPCARD, NPDAT, NPELAS, NPROD, NPSHEAR,   &
                                         NPSHEL, NPSOLID, NPLATEOFF, NPLATETHICK, NPLOAD4_3D, NPUSER1, NPUSERIN, NSEQ, NSUB,       &
                                         NTCARD, NTDAT, NTSUB, NVVEC, SOL_NAME

      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  CBMIN3, CBMIN4, IORQ1M, IORQ1S, IORQ1B, IORQ2B, IORQ2T

      USE MODEL_STUF, ONLY            :  BAROFF, BUSHOFF, EDAT, EOFF, EPNT, ESORT1, ESORT2, ETYPE, PLATEOFF, PLATETHICK, VVEC
      USE MODEL_STUF, ONLY            :  MATANGLE, MATL, RMATL, PBAR, RPBAR, PBEAM, RPBEAM, PBUSH, RPBUSH, PCOMP, RPCOMP, PELAS,   &
                                         RPELAS, PROD, RPROD, PSHEAR, RPSHEAR, PSHEL, PSOLID, RPSHEL, PUSER1, RPUSER1, PUSERIN,    &
                                         USERIN_MAT_NAMES
      USE MODEL_STUF, ONLY            :  ELDT, ELOUT, GROUT, OELDT, OELOUT, OGROUT, SCNUM, SUBLOD, TITLE, STITLE, LABEL
      USE MODEL_STUF, ONLY            :  GTEMP, TDATA, TPNT
      USE MODEL_STUF, ONLY            :  PDATA, PPNT, PLOAD4_3D_DATA, PTYPE
      USE MODEL_STUF, ONLY            :  ANY_ACCE_OUTPUT, ANY_DISP_OUTPUT, ANY_MPCF_OUTPUT, ANY_SPCF_OUTPUT, ANY_OLOA_OUTPUT,      &
                                         ANY_GPFO_OUTPUT, ANY_ELFE_OUTPUT, ANY_ELFN_OUTPUT, ANY_STRE_OUTPUT
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN
      USE TEMP_FILE_READERS, ONLY     :  READ_CHK
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  DATA_SET_NAME_ERROR, DATA_SET_SIZE_ERROR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK9S'
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: FILNAM            ! Name of a file that is to be opened for reading
      CHARACTER(132*BYTE)             :: MESSAG            ! Char message for file name
      CHARACTER(LEN=DATA_NAM_LEN)     :: NAME_Is           ! Name of data actually read from file
      CHARACTER(LEN=DATA_NAM_LEN)     :: NAME_ShouldBe     ! Name of data that should be read from file

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: INT2              ! An integer value read from a file from a file
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening or reading a file
      INTEGER(LONG)                   :: PCOMP_PLIES       ! Number of plies in 1 PCOMP entry incl sym plies not explicitly defined
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: REC_NO            ! Record number of a record read from a file
      INTEGER(LONG)                   :: UNT               ! Unit number of a file to be read




! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

!-----------------------------------------------------------------------------------------------------------------------------------
! Open L1D and read data

      FILNAM = LINK1D
      UNT    = L1D
      MESSAG = L1D_MSG

      CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

! Read subcase numbers, titling, load request data

      NAME_ShouldBe = 'S/C NUMBERS, TITLING, LOAD SET IDS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1D, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NSUB)     CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      DO I=1,NSUB

         READ(UNT,IOSTAT=IOCHK) SCNUM(I)                                       ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) TITLE(I)                                       ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) STITLE(I)                                      ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) LABEL(I)                                       ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) SUBLOD(I,1)                                    ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) SUBLOD(I,2)                                    ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      ENDDO

! Read OGROUT

      NAME_ShouldBe = 'OGROUT'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1D, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NSUB)     CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      DO I=1,NSUB
         READ(UNT,IOSTAT=IOCHK) OGROUT(I)                                      ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read OELOUT

      NAME_ShouldBe = 'OELOUT'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1D, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NSUB)     CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      DO I=1,NSUB
         READ(UNT,IOSTAT=IOCHK) OELOUT(I)                                      ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read OELDT

      NAME_ShouldBe = 'OELDT'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1D, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) OELDT                                             ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

! Read GROUT

      NAME_ShouldBe = 'GROUT'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1D, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NGRID)     CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NGRID', NGRID, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NSUB)     CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      DO I=1,NGRID
         DO J=1,NSUB
            READ(UNT,IOSTAT=IOCHK) GROUT(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read ELOUT

      NAME_ShouldBe = 'ELOUT'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1D, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NELE)     CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NELE', NELE, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NSUB)     CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      DO I=1,NELE
         DO J=1,NSUB
            READ(UNT,IOSTAT=IOCHK) ELOUT(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read ELDT

      NAME_ShouldBe = 'ELDT'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1D, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NELE)     CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NELE', NELE, INT2 )

      DO I=1,NELE
         READ(UNT,IOSTAT=IOCHK) ELDT(I)                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read ANY_xxxx_OUTPUT

      NAME_ShouldBe = 'ANY_xxxx_OUTPUT'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1D, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= ANY_ACCE_OUTPUT)      CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= ANY_DISP_OUTPUT)      CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= ANY_OLOA_OUTPUT)      CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= ANY_SPCF_OUTPUT)      CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= ANY_MPCF_OUTPUT)      CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= ANY_GPFO_OUTPUT)      CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= ANY_ELFN_OUTPUT)      CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= ANY_ELFE_OUTPUT)      CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= ANY_STRE_OUTPUT)      CALL DATA_SET_SIZE_ERROR ( LINK1D, NAME_Is, 'NSUB', NSUB, INT2 )

      IF ((SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC')) THEN
         CALL FILE_CLOSE ( L1D, LINK1D, L1DSTAT )
      ELSE
         CALL FILE_CLOSE ( L1D, LINK1D, 'KEEP' )
      ENDIF

! **********************************************************************************************************************************
! Open L1G

      FILNAM = LINK1G
      UNT    = L1G
      MESSAG = L1G_MSG

      CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

! Read ETYPE, EPNT, ESORT1 ESORT,2, EOFF

      NAME_ShouldBe = 'ETYPE, EPNT, ESORT1, ESORT2, EOFF'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NELE) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NELE', NELE, INT2 )
      DO I = 1,NELE

         READ(UNT,IOSTAT=IOCHK) ETYPE(I)                                       ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) EPNT(I)                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) ESORT1(I)                                      ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) ESORT2(I)                                      ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) EOFF(I)                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      ENDDO

! Read EDAT

      NAME_ShouldBe = 'EDAT'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NEDAT) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NEDAT', NEDAT, INT2 )
      DO I = 1,NEDAT
         READ(UNT,IOSTAT=IOCHK) EDAT(I)                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read element PARAMETERS

      NAME_ShouldBe = 'ELEM PARAMETERS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) IORQ1M                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) IORQ1S                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) IORQ1B                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) IORQ2B                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) IORQ2T                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) CBMIN3                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) CBMIN4                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

! Read BAR v vectors

      NAME_ShouldBe = 'V VECTORS IN GLOBAL COORDS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NVVEC) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NVVEC', NVVEC, INT2 )
      DO I=1,NVVEC
         READ(UNT,IOSTAT=IOCHK) (VVEC(I,J),J=1,3)                              ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read BAR offsets

      NAME_ShouldBe = 'BAR, BEAM OFFSETS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NBAROFF) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NBAROFF', NBAROFF, INT2 )
      DO I = 1,NBAROFF
         DO J = 1,6
            READ(UNT,IOSTAT=IOCHK) BAROFF(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
          ENDDO
      ENDDO

! Read BUSH offsets

      NAME_ShouldBe = 'BUSH OFFSETS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NBUSHOFF) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NBUSHOFF', NBUSHOFF, INT2 )
      DO I = 1,NBUSHOFF
         DO J = 1,6
            READ(UNT,IOSTAT=IOCHK) BUSHOFF(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
          ENDDO
      ENDDO

! Read plate offsets

      NAME_ShouldBe = 'PLATE OFFSETS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPLATEOFF) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPLATEOFF', NPLATEOFF, INT2 )
      DO I = 1,NPLATEOFF
         READ(UNT,IOSTAT=IOCHK) PLATEOFF(I)                                    ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read plate thicknesses from connection entries

      NAME_ShouldBe = 'PLATE THICKNESSES FROM CONNECTION ENTRIES'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPLATETHICK) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPLATETHICK', NPLATETHICK, INT2 )
      DO I = 1,NPLATETHICK
         READ(UNT,IOSTAT=IOCHK) PLATETHICK(I)                                  ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read PBAR, RPBAR

      NAME_ShouldBe = 'PBAR, RPBAR'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPBAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPBAR', NPBAR, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPBAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPBAR', MPBAR, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPBAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPBAR', MRPBAR, INT2 )

      DO I=1,NPBAR
         DO J=1,MPBAR
            READ(UNT,IOSTAT=IOCHK) PBAR(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPBAR
            READ(UNT,IOSTAT=IOCHK) RPBAR(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PBEAM, RPBEAM

      NAME_ShouldBe = 'PBEAM, RPBEAM'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPBEAM) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPBEAM', NPBEAM, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPBEAM) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPBEAM', MPBEAM, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPBEAM) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPBEAM', MRPBEAM, INT2 )

      DO I=1,NPBEAM
         DO J=1,MPBEAM
            READ(UNT,IOSTAT=IOCHK) PBEAM(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPBEAM
            READ(UNT,IOSTAT=IOCHK) RPBEAM(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PBUSH, RPBUSH

      NAME_ShouldBe = 'PBUSH, RPBUSH'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPBUSH) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPBUSH', NPBUSH, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPBUSH) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPBUSH', MPBUSH, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPBUSH) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPBUSH', MRPBUSH, INT2 )

      DO I=1,NPBUSH
         DO J=1,MPBUSH
            READ(UNT,IOSTAT=IOCHK) PBUSH(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPBUSH
            READ(UNT,IOSTAT=IOCHK) RPBUSH(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PROD, RPROD

      NAME_ShouldBe = 'PROD, RPROD'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPROD) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPROD', NPROD, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPROD) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPROD', MPROD, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPROD) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPROD', MRPROD, INT2 )

      DO I = 1,NPROD
         DO J=1,MPROD
            READ(UNT,IOSTAT=IOCHK) PROD(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPROD
            READ(UNT,IOSTAT=IOCHK) RPROD(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PELAS, RPELAS

      NAME_ShouldBe = 'PELAS, RPELAS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPELAS) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPELAS', NPELAS, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPELAS) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPELAS', MPELAS, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPELAS) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPELAS', MRPELAS, INT2 )

      DO I = 1,NPELAS
         DO J=1,MPELAS
            READ(UNT,IOSTAT=IOCHK) PELAS(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPELAS
            READ(UNT,IOSTAT=IOCHK) RPELAS(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PSHEAR, RPSHEAR

      NAME_ShouldBe = 'PSHEAR, RPSHEAR'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPSHEAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPSHEAR', NPSHEAR, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPSHEAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPSHEAR', MPSHEAR, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPSHEAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPSHEAR', MRPSHEAR, INT2 )

      DO I = 1,NPSHEAR
         DO J=1,MPSHEAR
            READ(UNT,IOSTAT=IOCHK) PSHEAR(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPSHEAR
            READ(UNT,IOSTAT=IOCHK) RPSHEAR(I,J)                                ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PSHEL, RPSHEL

      NAME_ShouldBe = 'PSHEL, RPSHEL'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPSHEL) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPSHEL', NPSHEL, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPSHEL) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPSHEL', MPSHEL, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPSHEL) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPSHEL', MRPSHEL, INT2 )

      DO I = 1,NPSHEL
         DO J=1,MPSHEL
            READ(UNT,IOSTAT=IOCHK) PSHEL(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPSHEL
            READ(UNT,IOSTAT=IOCHK) RPSHEL(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PCOMP, RPCOMP

      NAME_ShouldBe = 'PCOMP, RPCOMP'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPCOMP) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPCOMP', NPCOMP, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPCOMP0) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPCOMP0', MPCOMP0, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPCOMP_PLIES) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPCOMP_PLIES', MPCOMP_PLIES, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPCOMP0) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPCOMP0', MRPCOMP0, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPCOMP_PLIES) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPCOMP_PLIES', MRPCOMP_PLIES, INT2 )

      DO I = 1,NPCOMP

         READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         PCOMP_PLIES = INT2

         DO J=1,MPCOMP0+MPCOMP_PLIES*PCOMP_PLIES
            READ(UNT,IOSTAT=IOCHK) PCOMP(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO

         DO J=1,MRPCOMP0+MRPCOMP_PLIES*PCOMP_PLIES
            READ(UNT,IOSTAT=IOCHK) RPCOMP(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO

      ENDDO

! Read PSOLID

      NAME_ShouldBe = 'PSOLID'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPSOLID) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPSOLID', NPSOLID, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPSOLID) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPSOLID', MPSOLID, INT2 )

      DO I = 1,NPSOLID
         DO J=1,MPSOLID
            READ(UNT,IOSTAT=IOCHK) PSOLID(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PUSER1, RPUSER1

      NAME_ShouldBe = 'PUSER1, RPUSER1'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPUSER1) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPUSER1', NPUSER1, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPUSER1) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPUSER1', MPUSER1, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPUSER1) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPUSER1', MRPUSER1, INT2 )

      DO I = 1,NPUSER1
         DO J=1,MPUSER1
            READ(UNT,IOSTAT=IOCHK) PUSER1(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPUSER1
            READ(UNT,IOSTAT=IOCHK) RPUSER1(I,J)                                ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PUSERIN

      NAME_ShouldBe = 'PUSERIN'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPUSERIN) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPUSERIN', NPUSERIN, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPUSERIN) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPUSERIN', MPUSERIN, INT2 )

      DO I = 1,NPUSERIN
         DO J=1,MPUSERIN
            READ(UNT,IOSTAT=IOCHK) PUSERIN(I,J)                                ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read USERIN_MAT_NAMES

      NAME_ShouldBe = 'USERIN_MAT_NAMES'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPUSERIN) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPUSERIN', NPUSERIN, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MUSERIN_MAT_NAMES) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MUSERIN_MAT_NAMES', MUSERIN_MAT_NAMES, INT2 )

      DO I = 1,NPUSERIN
         DO J=1,MUSERIN_MAT_NAMES
            READ(UNT,IOSTAT=IOCHK) USERIN_MAT_NAMES(I,J)                       ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read material data from L1G

      NAME_ShouldBe = 'MATL, RMATL'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NMATL) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NMATL', NMATL, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MMATL) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MMATL', MMATL, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRMATLC) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRMATLC', MRMATLC, INT2 )

      DO I = 1,NMATL
         DO J=1,MMATL
            READ(UNT,IOSTAT=IOCHK) MATL(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRMATLC
            READ(UNT,IOSTAT=IOCHK) RMATL(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read material property angles

      NAME_ShouldBe = 'MATERIAL PROPERTY ANGLES'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NMATANGLE) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NMATANGLE', NMATANGLE, INT2 )
      DO I = 1,NMATANGLE
         READ(UNT,IOSTAT=IOCHK) MATANGLE(I)                                    ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

      IF ((SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC')) THEN
         CALL FILE_CLOSE ( L1G, LINK1G, L1GSTAT )
      ELSE
         CALL FILE_CLOSE ( L1G, LINK1G, 'KEEP' )
      ENDIF

! **********************************************************************************************************************************
! Open L1K and read data

      IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'BUCKLING') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN

         FILNAM = LINK1K
         UNT    = L1K
         MESSAG = L1K_MSG

         IF (NTCARD > 0) THEN

            CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

! Read TPNT

            NAME_ShouldBe = 'TPNT'
            REC_NO = 0

            READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1K, NAME_Is )

            READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (INT2 /= NELE) CALL DATA_SET_SIZE_ERROR ( LINK1K, NAME_Is, 'NELE', NELE, INT2 )

            READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (INT2 /= NTSUB) CALL DATA_SET_SIZE_ERROR ( LINK1K, NAME_Is, 'NTSUB', NTSUB, INT2 )

            DO I=1,NELE
               DO J=1,NTSUB
                  READ(UNT,IOSTAT=IOCHK) TPNT(I,J)                                ; REC_NO = REC_NO + 1
                  CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
               ENDDO
            ENDDO

! Read TDATA

            NAME_ShouldBe = 'TDATA'
            REC_NO = 0

            READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1K, NAME_Is )

            READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (INT2 /= NTDAT) CALL DATA_SET_SIZE_ERROR ( LINK1K, NAME_Is, 'NTDAT', NTDAT, INT2 )

            DO I=1,NTDAT
               READ(UNT,IOSTAT=IOCHK) TDATA(I)                                    ; REC_NO = REC_NO + 1
               CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            ENDDO

! Read GTEMP

            NAME_ShouldBe = 'GTEMP'
            REC_NO = 0

            READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1K, NAME_Is )

            READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (INT2 /= NGRID) CALL DATA_SET_SIZE_ERROR ( LINK1K, NAME_Is, 'NGRID', NGRID, INT2 )

            READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (INT2 /= NTSUB) CALL DATA_SET_SIZE_ERROR ( LINK1K, NAME_Is, 'NTSUB', NTSUB, INT2 )

            DO I=1,NGRID
               DO J=1,NTSUB
                  READ(UNT,IOSTAT=IOCHK) GTEMP(I,J)                               ; REC_NO = REC_NO + 1
                  CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
               ENDDO
            ENDDO

            IF ((SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC')) THEN
               CALL FILE_CLOSE ( L1K, LINK1K, L1KSTAT )
            ELSE
               CALL FILE_CLOSE ( L1K, LINK1K, 'KEEP' )
            ENDIF

         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Open L1Q and read data

      IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'BUCKLING') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN

         FILNAM = LINK1Q
         UNT    = L1Q
         MESSAG = L1Q_MSG

         IF (NPCARD > 0) THEN

            CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

! Read PPNT

            NAME_ShouldBe = 'PPNT'
            REC_NO = 0

            READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1Q, NAME_Is )

            READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (INT2 /= NELE) CALL DATA_SET_SIZE_ERROR ( LINK1Q, NAME_Is, 'NELE', NELE, INT2 )

            READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (INT2 /= NSUB) CALL DATA_SET_SIZE_ERROR ( LINK1Q, NAME_Is, 'NTSUB', NTSUB, INT2 )

            DO I=1,NELE
               DO J=1,NSUB
                  READ(UNT,IOSTAT=IOCHK) PPNT(I,J)                                ; REC_NO = REC_NO + 1
                  CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
               ENDDO
            ENDDO

! Read PDATA

            NAME_ShouldBe = 'PDATA'
            REC_NO = 0

            READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1Q, NAME_Is )

            READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (INT2 /= NPDAT) CALL DATA_SET_SIZE_ERROR ( LINK1Q, NAME_Is, 'NTDAT', NTDAT, INT2 )

            DO I=1,NPDAT
               READ(UNT,IOSTAT=IOCHK) PDATA(I)                                    ; REC_NO = REC_NO + 1
               CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            ENDDO

! Read PTYPE

            NAME_ShouldBe = 'PTYPE'
            REC_NO = 0

            READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1Q, NAME_Is )

            READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (INT2 /= NELE) CALL DATA_SET_SIZE_ERROR ( LINK1Q, NAME_Is, 'NELE' , NTDAT, INT2 )

            DO I=1,NELE
               READ(UNT,IOSTAT=IOCHK) PTYPE(I)                                    ; REC_NO = REC_NO + 1
               CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            ENDDO

! Read PLOAD4_3D_DATA

            NAME_ShouldBe = 'PLOAD4_3D_DATA'
            REC_NO = 0

            READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1Q, NAME_Is )

            READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            IF (INT2 /= NPLOAD4_3D) CALL DATA_SET_SIZE_ERROR ( LINK1Q, NAME_Is, 'NPLOAD4_3D', NPLOAD4_3D, INT2 )

            DO I=1,NPLOAD4_3D
               READ(UNT,IOSTAT=IOCHK) (PLOAD4_3D_DATA(I,J),J=1,MPLOAD4_3D_DATA)   ; REC_NO = REC_NO + 1
               CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            ENDDO

            IF ((SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC')) THEN
               CALL FILE_CLOSE ( L1Q, LINK1Q, L1QSTAT )
            ELSE
               CALL FILE_CLOSE ( L1Q, LINK1Q, 'KEEP' )
            ENDIF

         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE LINK9S


      SUBROUTINE MAXREQ_OGEL

! Count number of output requests to determine required leading dimension of array OGEL so memory can be allocated to it

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IBIT, LSUB, NDOFG, NELE, NGRID, METYPE, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  ELMTYP, ELOUT, ESORT2, ETYPE, GROUT, MEFFMASS_CALC, MPFACTOR_CALC, NELGP, NUM_PLIES,      &
                                         PCOMP_PROPS, SCNUM, TYPE
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_LOC, STRE_LOC, FORC_LOC
      USE LINK9_STUFF, ONLY           :  MAXREQ
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE COMPOSITE_SHELL_PREPARATION, ONLY:  GET_ELEM_NUM_PLIES, IS_ELEM_PCOMP_PROPS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MAXREQ_OGEL'
      CHARACTER(25*BYTE)              :: GROUT_NAME(0:15)
      CHARACTER(25*BYTE)              :: ELOUT_NAME(0:15)
      CHARACTER(46*BYTE)              :: MPF_MEFM_MSG         ! Message written if MAXREQ impacted by req for output of MPF/MEFM
      CHARACTER( 1*BYTE)              :: SKIP_ELEM_TYPE

      INTEGER(LONG)                   :: NUMBER_ROWS(0:15)    ! For elem output, this accounts for more than 1 row written to OGEL
      INTEGER(LONG)                   :: I,J,K,L              ! DO loop indices
      INTEGER(LONG)                   :: IB                   ! Result of IAND to determine if there is an output request
      INTEGER(LONG)                   :: INT_ELEM_ID          ! Internal element number
      INTEGER(LONG)                   :: MAXGROUT_SC          ! Max no. grid pt output requests in array GROUT for any 1 subcase
      INTEGER(LONG)                   :: MAXELOUT_SC          ! Max no. elem output requests in array ELOUT for any 1 subcase
      INTEGER(LONG)                   :: MAXGROUT             ! Max of MAXGROUT_SC for all subcases
      INTEGER(LONG)                   :: MAXELOUT             ! Max of MAXELOUT_SC for all subcases
      INTEGER(LONG)                   :: NREQ_EL(METYPE,0:15) ! No. of requests in ELOUT for each bit of ELOUT for 1 subcase
      INTEGER(LONG)                   :: NREQ_GR(0:15)        ! No. of requests in GROUT for each bit of GROUT for 1 subcase


      INTRINSIC                       :: IAND, MAX



! **********************************************************************************************************************************
! Initialize outputs

      MAXREQ = ZERO

      DO I=0,15
         GROUT_NAME(I) = '******** No Name ********'
         ELOUT_NAME(I) = '******** No Name ********'
      ENDDO

      GROUT_NAME(0) = 'grid displacement'
      GROUT_NAME(1) = 'grid applied load'
      GROUT_NAME(2) = 'grid SPC Force'
      GROUT_NAME(3) = 'grid MPC Force'
      GROUT_NAME(4) = 'grid point force balance'
      ELOUT_NAME(0) = 'element nodal force'
      ELOUT_NAME(1) = 'element engineering force'
      ELOUT_NAME(2) = 'element stress'
      ELOUT_NAME(3) = 'element strain'

! Count MAXGROUT, the max number of requests in GROUT.

      IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '11' )

      MAXGROUT = 0
      DO I=1,LSUB

         DO K=0,15
            NREQ_GR(K) = 0
         ENDDO

         DO J=1,NGRID
            DO K=0,15
               IB = IAND(GROUT(J,I),IBIT(K))
               IF (IB > 0) THEN
                  NREQ_GR(K) = NREQ_GR(K) + 1
               ENDIF
            ENDDO
         ENDDO

         MPF_MEFM_MSG(1:) = ' '
         IF ((MEFFMASS_CALC == 'Y') .OR. (MPFACTOR_CALC == 'Y')) THEN! Need to make sure MAXGROUT can cover MPF, MEFM
            IF (SOL_NAME /= 'GEN CB MODEL') THEN
               IF (NREQ_GR(2) < NDOFG) THEN
                  NREQ_GR(2) = NDOFG
                  MPF_MEFM_MSG = '(required for calculation of MPFACTOR/MEFMASS)'
               ENDIF
            ENDIF
         ENDIF

         IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '12' )

         MAXGROUT_SC = 0
         DO K=0,15
            IF (NREQ_GR(K) > 0) THEN
               IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '13' )
               IF (NREQ_GR(K) > MAXGROUT_SC) THEN
                  MAXGROUT_SC = NREQ_GR(K)
               ENDIF
               IF (NREQ_GR(K) > MAXGROUT) THEN
                  MAXGROUT = NREQ_GR(K)
               ENDIF
            ENDIF
         ENDDO
         IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '14' )
      ENDDO

      IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '15' )

! **********************************************************************************************************************************
! Count MAXELOUT, the max number of requests in ELOUT.

      DO I=0,15
         NUMBER_ROWS(I) = 0
      ENDDO

      MAXELOUT = 0
      DO I=1,LSUB

         DO L=1,METYPE
            DO K=0,15
               NREQ_EL(L,K) = 0
            ENDDO
         ENDDO

         DO L=1,METYPE
            DO J=1,NELE
               INT_ELEM_ID = ESORT2(J)
               DO K=0,15
                  TYPE = ETYPE(INT_ELEM_ID)
                  IF (TYPE == ELMTYP(L)) THEN
                     IB = IAND(ELOUT(INT_ELEM_ID,I),IBIT(K))
                     IF (IB > 0) THEN
                        CALL GET_NUM_ROWS_OUTPUT( L )      ! Determines the num rows of output for this element and request type
                        NREQ_EL(L,K) = NREQ_EL(L,K) + NUMBER_ROWS(K)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '21' )

         DO L=1,METYPE
            SKIP_ELEM_TYPE = 'Y'
            DO K=0,15
               IF (NREQ_EL(L,K) > 0) THEN
                  SKIP_ELEM_TYPE = 'N'
                  EXIT
               ENDIF
            ENDDO
            MAXELOUT_SC = 0
            IF (SKIP_ELEM_TYPE == 'N') THEN
               IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '22' )
               DO K=0,15
                  IF (NREQ_EL(L,K) > 0) THEN
                     IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '23' )
                     IF (NREQ_EL(L,K) > MAXELOUT_SC) THEN
                        MAXELOUT_SC = NREQ_EL(L,K)
                     ENDIF
                     IF (NREQ_EL(L,K) > MAXELOUT) THEN
                        MAXELOUT = NREQ_EL(L,K)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '24' )
         ENDDO

      ENDDO

      IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '25' )

      MAXREQ = MAX(MAXGROUT,MAXELOUT)

      IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '31' )



      RETURN

! **********************************************************************************************************************************




! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_NUM_ROWS_OUTPUT ( LETYPE )

      USE MODEL_STUF, ONLY            :  NUM_SEi

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: LETYPE            ! Index identifying which row to use from NUM_SEi array

! **********************************************************************************************************************************
      IF      (K == 0) THEN                                ! K = 0 is elem node force requests. (NELGP num of rows of output/elem)
!                                                            -----
         NUMBER_ROWS(K) = NELGP(L)

      ELSE IF (K == 1) THEN                                ! K = 1 is elem engr force output requests. (only 1 row of output/elem)
!                                                            -----
         NUMBER_ROWS(K) = 1

         IF (TYPE(1:5) == 'QUAD4') THEN
            IF (FORC_LOC == 'CENTER  ') THEN            !    PSHELL requires 2 rows of output/elem for FORC_LOC = 'CENTER'
               NUMBER_ROWS(K) = 1
            ELSE                                        !    PSHELL requires more lines of output for other FORC_LOC
               NUMBER_ROWS(K) = NUM_SEi(LETYPE)
            ENDIF
         ENDIF

      ELSE IF (K == 2) THEN                                ! K = 2 is elem stress output requests
!                                                            -----
         NUMBER_ROWS(K) = 1                                !    1 row of stress output/elem unless elem is BAR or shell

         CALL IS_ELEM_PCOMP_PROPS ( INT_ELEM_ID )          !    See if this is a PCOMP elem and get NUM_PLIES
         IF (PCOMP_PROPS == 'Y') THEN
            CALL GET_ELEM_NUM_PLIES ( INT_ELEM_ID )
         ENDIF

         IF       (TYPE(1:3) == 'BAR  ') THEN
               NUMBER_ROWS(K) = 2                          !    BAR stresses require 2 rows of output/elem
         ELSE IF ((TYPE(1:5) == 'TRIA3' ) .OR. (TYPE(1:5) == 'QUAD4')) THEN
            IF (PCOMP_PROPS == 'Y') THEN
               NUMBER_ROWS(K) = NUM_PLIES                  !    PCOMP requires NUM_PLIES rows of output/elem
            ELSE
               IF (STRE_LOC == 'CENTER  ') THEN            !    PSHELL requires 2 rows of output/elem for STRE_LOC = 'CENTER'
                  NUMBER_ROWS(K) = 2
               ELSE                                        !    PSHELL requires more lines of output for other STRE_LOC
                  NUMBER_ROWS(K) = 2*NUM_SEi(LETYPE)
               ENDIF
            ENDIF
         ELSE IF (TYPE(1:5) == 'QUAD8' ) THEN
            IF (PCOMP_PROPS == 'Y') THEN
               NUMBER_ROWS(K) = NUM_PLIES                  !    PCOMP requires NUM_PLIES rows of output/elem
            ELSE
               NUMBER_ROWS(K) = 2*NUM_SEi(LETYPE)          !    CQUAD8 stress output is CORNER even if CENTER is specified.
            ENDIF
         ELSE IF ((TYPE(1:4) == 'HEXA' ) .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
            NUMBER_ROWS(K) = NUM_SEi(LETYPE)
         ENDIF

      ELSE IF (K == 3) THEN                                ! K = 3 is elem strain output requests
!                                                            -----
         NUMBER_ROWS(K) = 1                                !    1 row of strain output/elem unless elem is BAR or shell

         CALL IS_ELEM_PCOMP_PROPS ( INT_ELEM_ID )          !    See if this is a PCOMP elem and get NUM_PLIES
         IF (PCOMP_PROPS == 'Y') THEN
            CALL GET_ELEM_NUM_PLIES ( INT_ELEM_ID )
         ENDIF

         IF ((TYPE(1:5) == 'TRIA3' ) .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'SHEAR')) THEN
            IF (PCOMP_PROPS == 'Y') THEN
               NUMBER_ROWS(K) = NUM_PLIES                  !    PCOMP requires NUM_PLIES rows of output/elem
            ELSE
               IF (STRN_LOC == 'CENTER  ') THEN            !    PSHELL requires 2 rows of output/elem for STRN_LOC = 'CENTER'
                  NUMBER_ROWS(K) = 2
               ELSE                                        !    PSHELL requires more lines of output for other STRN_LOC
                  NUMBER_ROWS(K) = 2*NUM_SEi(LETYPE)
               ENDIF
            ENDIF
         ELSE IF ((TYPE(1:4) == 'HEXA' ) .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
            NUMBER_ROWS(K) = NUM_SEi(LETYPE)
         ENDIF

      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE GET_NUM_ROWS_OUTPUT

! ##################################################################################################################################

      SUBROUTINE MAXREQ_OGEL_DEB ( WHICH )

      CHARACTER( 2*BYTE)              :: WHICH             ! Decides what to print out for this call to this subr

! **********************************************************************************************************************************
      IF      (WHICH == '11') THEN

         WRITE(F06,1101)
         WRITE(F06,1102)
         WRITE(F06,1103)

      ELSE IF (WHICH == '12') THEN

         IF  ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
            WRITE(F06,1201) SCNUM(I)
         ELSE
            WRITE(F06,1202)
         ENDIF

      ELSE IF (WHICH == '13') THEN

         WRITE(F06,1301) GROUT_NAME(K), NREQ_GR(K), MPF_MEFM_MSG

      ELSE IF (WHICH == '14') THEN

         IF (MAXGROUT_SC > 0) THEN
            IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
               WRITE(F06,1401) MAXGROUT_SC
            ENDIF
         ENDIF

      ELSE IF (WHICH == '15') THEN

         IF  ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
            WRITE(F06,1501) MAXGROUT
         ELSE
            WRITE(F06,1502) MAXGROUT
         ENDIF
         WRITE(F06,1103)

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (WHICH == '21') THEN

         IF  ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
            WRITE(F06,2101) SCNUM(I)
         ELSE
            WRITE(F06,2102)
         ENDIF

      ELSE IF (WHICH == '22') THEN

         WRITE(F06,2201) ELMTYP(L)

      ELSE IF (WHICH == '23') THEN

         WRITE(F06,1301) ELOUT_NAME(K), NREQ_EL(L,K)

      ELSE IF (WHICH == '24') THEN

         IF (MAXELOUT_SC > 0) THEN
            WRITE(F06,2401) ELMTYP(L), MAXELOUT_SC
         ENDIF

      ELSE IF (WHICH == '25') THEN

         WRITE(F06,2501) MAXELOUT
         WRITE(F06,1103)

      ELSE IF (WHICH == '31') THEN

         WRITE(F06,3101) MAXREQ
         WRITE(F06,3102)

      ENDIF

! **********************************************************************************************************************************
 1101 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' :::::::::::::::::::::::::::::::::::::::START DEBUG(91) OUTPUT FROM SUBROUTINE MAXREQ_OGEL:::::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1102 FORMAT(10X,'C A L C U L A T E   N U M B E R   O F   R O W S   O F   O U T P U T   N E E D E D   F O R   A R R A Y   O G E L'/)

 1103 FORMAT(1X,'****************************************************************************************************************',&
                '*******************')

 1201 FORMAT(1X,/,' Grid point related output requests for subcase ',I8                                                          ,/&
                  ' -------------------------------------------------------',/)

 1202 FORMAT(1X,/,' Grid point related output requests'                                                                          ,/&
                  ' ----------------------------------',/)

 1301 FORMAT('   Number of rows of output needed for ',A,':',I8,2X,A)

 1401 FORMAT(1X,/,'   The maximum number of rows needed for grid related outputs for this subcase is    : ',I8,/)

 1501 FORMAT(1X,/,'   The maximum number of rows needed for grid related outputs for all subcases is    : ',I8,/)

 1502 FORMAT(1X,/,'   The maximum number of rows needed for grid related outputs is                     : ',I8,/)

 2101 FORMAT(1X,/,' Element related output requests for subcase ',I8                                                             ,/&
             ' ----------------------------------------------------',/)

 2102 FORMAT(1X,/,' Element related output requests'                                                                             ,/&
                  ' -------------------------------',/)

 2201 FORMAT('   For element type: ',A                                                                                          ,/,&
             '   -------------------------')

 2401 FORMAT(1X,/,'   The maximum number of rows needed for ',A,' element outputs for this subcase is: ',I8,/)

 2402 FORMAT(1X,/,'   The maximum number of rows needed for ',A,' element outputs is                    : ',I8,/)

 2501 FORMAT(1X,/,'   The maximum number of rows needed for all element outputs for all subcases is     : ',I8,/)

 3101 FORMAT(1X,/,'   The maximum number of rows of output needed for any subcase is                    : ',I8,/)

 3102 FORMAT(' ::::::::::::::::::::::::::::::::::::::END DEBUG(91) OUTPUT FROM SUBROUTINE MAXREQ_OGEL::::::::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE MAXREQ_OGEL_DEB

      END SUBROUTINE MAXREQ_OGEL


      SUBROUTINE ALLOCATE_LINK9_STUF ( CALLING_SUBR )

! Allocate some arrays for use in LINK9

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, ONEPP6
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MELGP, MMSPRNT, MOGEL, TOT_MB_MEM_ALLOC
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE TIMDAT, ONLY                :  TSEC
      USE LINK9_STUFF, ONLY           :  GID_OUT_ARRAY, EID_OUT_ARRAY, FTNAME, MAXREQ, MSPRNT, OGEL, POLY_FIT_ERR,                 &
                                         POLY_FIT_ERR_INDEX

      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_LINK9_STUF'
      CHARACTER(24*BYTE)              :: NAME              ! Array name (used for output error message)
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! STAT from ALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: NROWS             ! Nunber of rows in array NAME being allocated
      INTEGER(LONG)                   :: NCOLS             ! Nunber of cols in array NAME being allocated


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called
      REAL(DOUBLE)                    :: RDOUBLE           ! Real value of DOUBLE
      REAL(DOUBLE)                    :: RLONG             ! Real value of LONG

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
      RDOUBLE = REAL(DOUBLE)
      RLONG   = REAL(LONG)

      MB_ALLOCATED = ZERO
      NROWS = MAXREQ
      JERR = 0

! Allocate array for GID_OUT_ARRAY

      NAME = 'GID_OUT_ARRAY'
      NCOLS = MELGP
      IF (ALLOCATED(GID_OUT_ARRAY)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (GID_OUT_ARRAY(MAXREQ,MELGP+1),STAT=IERR)
         MB_ALLOCATED = RLONG*REAL(MAXREQ)*REAL(MELGP)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            DO I=1,MAXREQ
               DO J=1,MELGP+1
                  GID_OUT_ARRAY(I,J) = 0
               ENDDO
            ENDDO
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Allocate arrays for EID_OUT_ARRAY

      NAME = 'EID_OUT_ARRAY'
      NCOLS = 1
      IF (ALLOCATED(EID_OUT_ARRAY)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (EID_OUT_ARRAY(MAXREQ,2),STAT=IERR)
         MB_ALLOCATED = RLONG*TWO*REAL(MAXREQ)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            DO I=1,MAXREQ
               EID_OUT_ARRAY(I,1) = 0
               EID_OUT_ARRAY(I,2) = 1                      ! Set number of plies to 1 so that we don't have to worry about it
            ENDDO                                          ! except for 2D elems (which will be set in another subr)
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Allocate arrays for FTNAME

      NAME = 'FTNAME                  '
      IF (ALLOCATED(FTNAME)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (FTNAME(MAXREQ),STAT=IERR)
         IF (IERR == 0) THEN
         MB_ALLOCATED =REAL(LEN(FTNAME)*MAXREQ)/ONEPP6
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            DO I=1,MAXREQ
               FTNAME(I) = 'none'
            ENDDO
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Allocate arrays for MSPRNT

      NAME = 'MSPRNT                  '
      NCOLS = MMSPRNT
      IF (ALLOCATED(MSPRNT)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (MSPRNT(MAXREQ,MMSPRNT),STAT=IERR)
         IF (IERR == 0) THEN
         MB_ALLOCATED =REAL(LEN(MSPRNT)*MAXREQ*MMSPRNT)/ONEPP6
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            DO I=1,MAXREQ
               DO J=1,MMSPRNT
                  MSPRNT(I,J)(1:) = ' '
               ENDDO
            ENDDO
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Allocate arrays for OGEL

      NAME = 'OGEL'
      NCOLS = MOGEL
      IF (ALLOCATED(OGEL)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (OGEL(MAXREQ,MOGEL),STAT=IERR)
         MB_ALLOCATED = RDOUBLE*REAL(MAXREQ)*REAL(MOGEL)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            DO I=1,MAXREQ
               DO J=1,MOGEL
                  OGEL(I,J) = ZERO
               ENDDO
            ENDDO
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Allocate arrays for POLY_FIT_ERR

      NAME = 'POLY_FIT_ERR'
      NCOLS = 1
      IF (ALLOCATED(POLY_FIT_ERR)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (POLY_FIT_ERR(MAXREQ),STAT=IERR)
         MB_ALLOCATED = RDOUBLE*REAL(MAXREQ)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            DO I=1,MAXREQ
               POLY_FIT_ERR(I) = ZERO
            ENDDO
         ELSE
            WRITE(ERR,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            WRITE(F06,991) MB_ALLOCATED,NAME,SUBR_NAME,IERR
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ENDIF
      ENDIF

! Allocate arrays for POLY_FIT_ERR_INDEX

      NAME = 'POLY_FIT_ERR_INDEX'
      NCOLS = 1
      IF (ALLOCATED(POLY_FIT_ERR_INDEX)) THEN
         WRITE(ERR,990) SUBR_NAME, NAME
         WRITE(F06,990) SUBR_NAME, NAME
         FATAL_ERR = FATAL_ERR + 1
         JERR = JERR + 1
      ELSE
         ALLOCATE (POLY_FIT_ERR_INDEX(MAXREQ),STAT=IERR)
         MB_ALLOCATED = BYTE*REAL(MAXREQ)/ONEPP6
         IF (IERR == 0) THEN
            CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            DO I=1,MAXREQ
               POLY_FIT_ERR_INDEX(I) = 0
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
         WRITE(ERR,1699) TRIM(SUBR_NAME), CALLING_SUBR
         WRITE(F06,1699) SUBR_NAME, CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_LINK9_STUF


      SUBROUTINE DEALLOCATE_LINK9_STUF

! Deallocate some arrays used in LINK9

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE LINK9_STUFF, ONLY           :  GID_OUT_ARRAY, EID_OUT_ARRAY, FTNAME, MSPRNT, OGEL, POLY_FIT_ERR,              &
                                         POLY_FIT_ERR_INDEX

      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  ALLOCATED_MEMORY
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_LINK9_STUF'
      CHARACTER(24*BYTE)              :: NAME              ! Array name (used for output error message)

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

! Deallocate array GID_OUT_ARRAY

      IF (ALLOCATED(GID_OUT_ARRAY)) THEN
         DEALLOCATE (GID_OUT_ARRAY,STAT=IERR)
         NAME = 'GID_OUT_ARRAY'
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
         ENDIF
      ENDIF

! Deallocate array EID_OUT_ARRAY

      IF (ALLOCATED(EID_OUT_ARRAY)) THEN
         DEALLOCATE (EID_OUT_ARRAY,STAT=IERR)
         NAME = 'EID_OUT_ARRAY'
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
         ENDIF
      ENDIF

! Deallocate array FTNAME

      IF (ALLOCATED(FTNAME)) THEN
         DEALLOCATE (FTNAME,STAT=IERR)
         NAME = 'FTNAME                  '
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
         ENDIF
      ENDIF

! Deallocate array MSPRNT

      IF (ALLOCATED(MSPRNT)) THEN
         DEALLOCATE (MSPRNT,STAT=IERR)
         NAME = 'MSPRNT                  '
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
         ENDIF
      ENDIF

! Deallocate array OGEL

      IF (ALLOCATED(OGEL)) THEN
         DEALLOCATE (OGEL,STAT=IERR)
         NAME = 'OGEL                  '
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
         ENDIF
      ENDIF

! Deallocate array POLY_FIT_ERR

      IF (ALLOCATED(POLY_FIT_ERR)) THEN
         DEALLOCATE (POLY_FIT_ERR,STAT=IERR)
         NAME = 'POLY_FIT_ERR                  '
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
         ENDIF
      ENDIF

! Deallocate array POLY_FIT_ERR_INDEX

      IF (ALLOCATED(POLY_FIT_ERR_INDEX)) THEN
         DEALLOCATE (POLY_FIT_ERR_INDEX,STAT=IERR)
         NAME = 'POLY_FIT_ERR_INDEX                  '
         CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
         IF (IERR /= 0) THEN
            WRITE(ERR,992) NAME,SUBR_NAME
            WRITE(F06,992) NAME,SUBR_NAME
            JERR = JERR + 1
         ENDIF
      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)

! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_LINK9_STUF

   END MODULE LINK9_MOD
