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

   MODULE LINK0_LINK1_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: LINK0, LINK1

   CONTAINS

      SUBROUTINE LINK0

! LINK0:

!  - reads in the data deck
!  - checks for restart
!  - processes grid and coord data
!  - sequences grids
!  - process Case Control output requests
!  - forms the DOF tables
!  - processes concentrated masses
!  - calcs rigid body mass props (GPWG)
!  - process temperature and pressure load input data to get arrays needed for element load calcs

      USE PENTIUM_II_KIND, ONLY       :  BYTE, SHORT, LONG, SINGLE, DOUBLE, QUAD

      USE IOUNT1, ONLY                :  MOU4, SC1, WRT_BUG
      USE IOUNT1, ONLY                :  BUG, ERR, F06, F21, F22, F23, F24, F25, IN1, L1B, L1C, L1D, L1F, L1G, L1H, L1I, L1K, &
                                         L1L, L1N, L1O, L1P, L1Q, L1S, L1T, L1U, L1V, L1W, L1X, L1Y, OP2, OU4, SEQ

      USE IOUNT1, ONLY                :  F21FIL, F22FIL, F23FIL, F24FIL, F25FIL, INFILE, LINK1B, LINK1C, LINK1D,              &
                                         LINK1F, LINK1H, LINK1I, LINK1K, LINK1L, LINK1N, LINK1O, LINK1P, LINK1Q, LINK1S, LINK1T,   &
                                         LINK1U, LINK1V, LINK1W, LINK1X, LINK1Y, OP2FIL, OU4FIL, SEQFIL

      USE IOUNT1, ONLY                :  L1LSTAT, L1NSTAT, L1OSTAT, L1QSTAT, L1YSTAT, OP2STAT

      USE IOUNT1, ONLY                :  F21_MSG, F22_MSG, F23_MSG, F24_MSG, F25_MSG, L1B_MSG, L1C_MSG, L1D_MSG, L1F_MSG, L1H_MSG, &
                                         L1I_MSG, L1K_MSG, L1L_MSG, L1N_MSG, L1O_MSG, L1P_MSG, L1Q_MSG, L1S_MSG, L1T_MSG, L1U_MSG, &
                                         L1V_MSG, L1W_MSG, L1X_MSG, L1Y_MSG, OP2_MSG, OU4_MSG, SEQ_MSG

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, BANDIT_ERR, CHKPNT, COMM, DEMO_GRID_LIMIT, ENFORCED, EPSIL1_SET, FATAL_ERR, &
                                         IBIT, KMAT_BW, KMAT_DEN, LGRID, LINKNO, MBUG, MELDTS, NAOCARD, NCUSERIN,                  &
                                         NDOFG, NDOFR, NDOFSE, NELE, NFORCE, NGRAV, NGRID, NMPC, NPCARD,                           &
                                         NPUSERIN, NRFORCE, NRIGEL, NSLOAD, NSNORM, NSPC, NSPC1, NTCARD, NTERM_KGG,                &
                                         NUM_PARTVEC_RECORDS, NUM_SUPT_CARDS, NUM_USET_RECORDS, PROG_NAME, RESTART, SOL_NAME,      &
                                         WARN_ERR

      USE SCONTR, ONLY                :  ELDT_BUG_DAT1_BIT, ELDT_BUG_DAT2_BIT, ELDT_BUG_ME_BIT  , ELDT_BUG_P_T_BIT ,               &
                                         ELDT_BUG_SE_BIT  , ELDT_BUG_KE_BIT  , ELDT_BUG_SHPJ_BIT, ELDT_BUG_BMAT_BIT,               &
                                         ELDT_BUG_BCHK_BIT, ELDT_BUG_U_P_BIT,  ELDT_F21_P_T_BIT , ELDT_F22_ME_BIT  ,               &
                                         ELDT_F23_KE_BIT  , ELDT_F24_SE_BIT  , ELDT_F25_U_P_BIT
      use scontr, only                :  ndofo
      USE DOF_TABLES, ONLY            :  TDOFI
      USE PARAMS, ONLY                :  CHKGRDS, EPSIL, EQCHK_OUTPUT, GRDPNT, GRDPNT_IN, GRIDSEQ, MEFMGRID, MEFMLOC, PRTCONN,     &
                                         PRTBASIC, PRTCORD, PRTDOF, PRTTSET, PRTSTIFD, PRTSTIFF, SETLKTK, SETLKTM, SUPINFO,        &
                                         SUPWARN, WTMASS
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MACHINE_PARAMS, ONLY        :  MACH_PREC
      USE MODEL_STUF, ONLY            :  ANY_GPFO_OUTPUT, EIG_METH, ELDT, ETYPE, MEFFMASS_CALC, NUM_EMG_FATAL_ERRS, PLY_NUM, OELDT,&
                                         GRID_ID, SNORM, RSNORM, GRID_SNORM
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE BANDIT_MODULE
      USE RIGID_BODY_DISP_MATS, ONLY  :  RBGLOBAL_GSET, TR6_CG, TR6_MEFM, TR6_0
      USE SPARSE_MATRICES, ONLY       :  SYM_KGG, I_KGG, J_KGG, KGG
      USE OUTPUT4_MATRICES, ONLY      :  NUM_OU4_REQUESTS, OU4_FILE_UNITS

      USE DATE_TIME_UTILS, ONLY       :  OURDAT, OURTIM, TIME_INIT
      USE INPUT_STAGE_DISPATCH, ONLY  :  LOADB, LOADB0, LOADB_RESTART, LOADC, LOADC0, LOADE
      USE IN4_FILE_LIFECYCLE, ONLY    :  ALLOCATE_IN4_FILES
      USE NONLINEAR_PARAMETER_LIFECYCLE, ONLY:  ALLOCATE_NL_PARAMS
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_INQUIRE, FILE_OPEN, SET_FILE_CLOSE_STAT
      USE MODEL_STORAGE_ALLOCATION, ONLY:  ALLOCATE_MODEL_STUF
      USE TEMP_FILE_READERS, ONLY     :  READ_DOF_TABLES, READ_L1Z
      USE LINK1_WORKFLOW_SUPPORT, ONLY:  LINK1_RESTART_DATA, PRINT_CONSTANTS_1, PRINT_ORDER, WRITE_ENF_TO_L1O
      USE SPARSE_FORMAT_CONVERSION, ONLY:  SET_SPARSE_MAT_SYM
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_ARRAY_ALLOC_STAT, GET_MACHINE_PARAMS, GET_MATRIX_DIAG_STATS, LINK_MESSAGE, WRITE_ALLOC_MEM_TABLE
      USE MODEL_STORAGE_DEALLOCATION, ONLY:  DEALLOCATE_MODEL_STUF
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_EDAT, WRITE_FIJFIL, WRITE_L1Z, WRITE_TDOF, WRITE_TSET, WRITE_USETSTR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE, WRITE_L1A
      USE ELEMENT_MODEL_INDEXING, ONLY:  ELESORT, ELEM_PROP_MATL_IIDS
      USE GRID_COORDINATE_PROCESSING, ONLY:  CORD_PROC, GRID_PROC, SEQ_PROC
      USE ELEMENT_LOOKUPS, ONLY       :  GET_ELEM_AGRID_BGRID, GRID_ELEM_CONN_TABLE
      USE SUBCASE_PROCESSING, ONLY    :  SUBCASE_PROC
      USE MATRIX_TEXT_OUTPUT, ONLY    :  WRITE_GRID_COORDS
      USE DOF_TABLE_LIFECYCLE, ONLY   :  ALLOCATE_DOF_TABLES
      USE DOF_ARRAY_INDEXING, ONLY    :  CALC_TDOF_ROW_START, GET_ARRAY_ROW_NUM
      USE EMG_MOD, ONLY               :  EMG
      USE DOF_NUMBERING, ONLY         :  DOF_PROC, TDOF_COL_NUM
      USE CONM2_PROCESSING, ONLY      :  CONM2_PROC_1, CONM2_PROC_2
      USE RIGID_BODY_STORAGE_LIFECYCLE, ONLY:  ALLOCATE_RBGLOBAL
      USE GRID_POINT_WEIGHT, ONLY     :  GPWG, GPWG_USERIN, RB_DISP_MATRIX_PROC
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  READ_MATRIX_1
      USE RESTART_FILE_IO, ONLY       :  PRT_MATS_ON_RESTART
      USE SPARSE_MATRIX_DEALLOCATION, ONLY:  DEALLOCATE_SPARSE_MAT
      USE REDUCTION_CHECKS, ONLY      :  STIFF_MAT_EQUIL_CHK
      USE LINK1_LOAD_SUPPORT, ONLY    :  PRESSURE_DATA_PROC, TEMPERATURE_DATA_PROC, YS_ARRAY
      USE COL_VEC_LIFECYCLE, ONLY     :  ALLOCATE_COL_VEC, DEALLOCATE_COL_VEC
      USE OUTPUT4_FILE_IO, ONLY       :  OUTPUT4_MATRIX_MSGS
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      USE OP2_GEOMETRY_OUTPUT, ONLY   :  WRITE_OP2_HEADER
      IMPLICIT NONE

      LOGICAL                         :: FILE_OPND         ! = .TRUE. if a file is opened

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK0'

      CHARACTER(1*BYTE)               :: OPEN_UNT(MOU4)    ! 'Y' if we need to open an OU4 file unit

                                                           ! Option indicator for subr EMG
      CHARACTER(1*BYTE)               :: OPT(6)    = (/'N', 'N', 'N', 'N', 'N', 'N'/)

      CHARACTER( 1*BYTE)              :: RBG_GSET_ALLOCATED! 'Y' or 'N' depending on whether array RBGLOBAL has been allocated
      CHARACTER(132*BYTE)             :: TDOF_MSG          ! Message to be printed out regarding at what point in the run the TDOF,I
!                                                            tables are printed out

      INTEGER(LONG)                   :: BANDIT_BW         ! Matrix bandwidth returned from subr BANDIT

                                                           ! Array used to tell subr ELDT_PROC_FOR_RESTART which ELDT to calc
      INTEGER(LONG)                   :: CHK_ELDT_BIT(0:MELDTS-1)

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: I1,I2             ! Intermediate integer variable
      INTEGER(LONG)                   :: JERROR    = 0     ! Local error count
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: PRINTIT   = 0     ! 'Y'/'N' depending on whether to tell subr GET_MATRIX_DIAG_STATS to
!                                                            print output
      INTEGER(LONG)                   :: R_SET_COL         ! Col number in TDOF where the R-set exists
      INTEGER(LONG)                   :: R_SET_DOF         ! Comp number in the R-set
      INTEGER(LONG)                   :: SETLKTK_DEF       ! Default value of SETLKTK
      INTEGER(LONG)                   :: SETLKTM_DEF       ! Default value of SETLKTM
      INTEGER(LONG) :: POST
      INTEGER(LONG)                   :: BGRID             ! Internal grid point number

      REAL(DOUBLE)                    :: KGG_DIAG(NDOFG)   ! Diagonal of KGG (needed for equil check on RESTART)
      REAL(DOUBLE)                    :: KGG_MAX_DIAG      ! Max diag term from KGG (needed for equil check on RESTART)
      !LOGICAL                        :: WRITE_F06  ! flag
      !LOGICAL                        :: WRITE_OP2  ! flag

      INTRINSIC                       :: IAND

! **********************************************************************************************************************************
      LINKNO = 0
      POST = -1

      ! Initialize WRT_BUG
      DO I=0,MBUG-1
         WRT_BUG(I) = 0
      ENDDO

      RBG_GSET_ALLOCATED = 'N'

      ! Set default values for SETLKT from values in module PARAMS
      SETLKTK_DEF = SETLKTK
      SETLKTM_DEF = SETLKTM

      ! Set time initializing parameters
      CALL TIME_INIT

      ! Get date and time, write to screen
      CALL OURDAT
      CALL OURTIM
      WRITE(SC1,152) LINKNO

      ! Write logo, date and copyright info to text files
      WRITE(F06,150) LINKNO
      WRITE(ERR,150) LINKNO

      ! Reset units for error output (were set to SC1 in MAIN1)
      OUNT(1) = ERR
      OUNT(2) = F06

      ! Read input data to count sizes of arrays (no. GRID's, elems, etc.)
      CALL LINK_MESSAGE('DETERMINE ARRAY SIZES - CASE CONTROL        ')
      CALL LOADC0

! Initial pass on Bulk Data when this is not a restart

res11:IF (RESTART == 'N') THEN
         CALL LINK_MESSAGE('DETERMINE ARRAY SIZES - BULK DATA           ')
         CALL LOADB0
      ENDIF res11

      ! Start back at beginnng of input file to read
      ! Exec, Case Control and Bulk Data decks
      REWIND (IN1)

      ! Processes the EXEC CONTROL DECK
      CALL LINK_MESSAGE('READ EXEC CONTROL DECK                      ')
      WRITE(F06,*)
      CALL ALLOCATE_IN4_FILES ( 'IN4FIL', 0, 0, SUBR_NAME )
      CALL LOADE

      IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
         CALL ALLOCATE_NL_PARAMS ( 'MYSTRAN' )
      ENDIF

! pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
! Code here is used to test whether user has permission to run problems.

! Code removed 10/02/15

! pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp

      ! If CHKPNT was in the Exec Control deck, reset file close status' to KEEP
      IF (CHKPNT == 'Y') THEN
         IF (SOL_NAME(1:8) /= 'NLSTATIC') THEN
            CALL SET_FILE_CLOSE_STAT ( 'KEEP' )
         ELSE
            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,9401) SOL_NAME
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9401) SOL_NAME
            ENDIF
         ENDIF
      ENDIF

      ! If nonlinear we need to keep files for iterations
      IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
         CALL SET_FILE_CLOSE_STAT ( 'KEEP' )
      ENDIF

      ! Processes the Case Control deck and write
      ! check-point data to L1Z if this is NOT a restart or
      ! verify that the check-pointed data
      ! agrees between the original run and this restart run
      CALL LINK_MESSAGE('ALLOCATE MEMORY FOR SOME MODEL DATA ARRAYS  ')
      CALL ALLOCATE_MODEL_STUF ( 'SETS ARRAYS', SUBR_NAME )
      CALL ALLOCATE_MODEL_STUF ( 'TITLES', SUBR_NAME )
      CALL ALLOCATE_MODEL_STUF ( 'SC_xxxx', SUBR_NAME )
      CALL ALLOCATE_MODEL_STUF ( 'SCNUM', SUBR_NAME )
      CALL ALLOCATE_MODEL_STUF ( 'SUBLOD', SUBR_NAME )
      CALL ALLOCATE_MODEL_STUF ( 'SPC_MPC_SET', SUBR_NAME )

      CALL LINK_MESSAGE('READ CASE CONTROL DECK                      ')
      CALL LOADC

! If this is a restart, do the following:

res12:IF (RESTART == 'Y') THEN

         CALL READ_L1Z

         CALL ALLOCATE_MODEL_STUF ( 'CONM2, RCONM2', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'CORD, RCORD', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'GRID, RGRID', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'GRID_ID', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'GRID_SEQ, INV_GRID_SEQ', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'ETYPE, EDAT, EPNT', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'ESORT1', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'ESORT2', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'EOFF', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'ELEM PROPERTIES AND MATERIALS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'VVEC, OFFSETS, PLATE stuff', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'PPNT, PDATA, PTYPE', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'PLOAD4_3D_DATA', SUBR_NAME )

         CALL LINK1_RESTART_DATA

         CALL LOADB_RESTART                                ! Read new Bulk Data entries allowed

         CALL SET_SPARSE_MAT_SYM                           ! Uses SPARSTOR which was read from L1A

         IF (DEBUG(1) == 1) THEN                           ! Print KIND parameters if DEBUG says to
            WRITE(F06,*)
            WRITE(F06,11) BYTE,SHORT,LONG,SINGLE,DOUBLE,QUAD
            WRITE(F06,*)
         ENDIF

         IF (DEBUG(2) == 1) THEN                           ! Print CONSTANTS_1 if DEBUG says to
            CALL PRINT_CONSTANTS_1
         ENDIF

         CALL GET_MACHINE_PARAMS
         IF (EPSIL1_SET == 'N') THEN
            EPSIL(1) = MACH_PREC
         ENDIF

         IF (DEBUG(5) == 1) THEN                           ! Print Gauss quadriture abscissas and weights if DEBUG says to
            CALL PRINT_ORDER
         ENDIF

      ENDIF res12

res13:IF (RESTART == 'N') THEN

! Open files needed for writing Bulk Data in LOADB.
! L1F: Rigid element data
! L1I: FORCE, MOMENT for loads requested in Case Control
! L1K: Temperature data for temp loads requested in Case Control
! L1M: Eigen extraction data for METHOD requested in Case Control
! L1N: ASET, OMIT data
! L1O: SPC data for SPC's requested in Case Control
! L1P: GRAV data for loads requested in Case Control
! L1Q: Element pressure load data for loads requested in Case Control
! L1S: Data from B.D. MPC cards
! L1T:
! L1U: Data from B.D. RFORCE cards
! L1W: Data from B.D. SLOAD  cards

         CALL FILE_OPEN ( L1F, LINK1F, OUNT, 'REPLACE', L1F_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1I, LINK1I, OUNT, 'REPLACE', L1I_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1K, LINK1K, OUNT, 'REPLACE', L1K_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1N, LINK1N, OUNT, 'REPLACE', L1N_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1O, LINK1O, OUNT, 'REPLACE', L1O_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1P, LINK1P, OUNT, 'REPLACE', L1P_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1Q, LINK1Q, OUNT, 'REPLACE', L1Q_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1S, LINK1S, OUNT, 'REPLACE', L1S_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1T, LINK1T, OUNT, 'REPLACE', L1T_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1U, LINK1U, OUNT, 'REPLACE', L1U_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1V, LINK1V, OUNT, 'REPLACE', L1V_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1W, LINK1W, OUNT, 'REPLACE', L1W_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL FILE_OPEN ( L1X, LINK1X, OUNT, 'REPLACE', L1X_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
!xx      IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
!xx         CALL FILE_OPEN ( L1M, LINK1M, OUNT, 'REPLACE', L1M_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
!xx      ENDIF

         ! Get machine parameters and set EPSIL(1).
         ! This may change later if user inputs a value on a PARAM EPSIL entry.
         CALL GET_MACHINE_PARAMS
         IF (EPSIL1_SET == 'N') THEN
            EPSIL(1) = MACH_PREC
         ENDIF

         ! Processes the Bulk Data deck
!         CALL LINK_MESSAGE('ALLOCATE MEMORY FOR SOME MODEL DATA ARRAYS  ')
         CALL ALLOCATE_MODEL_STUF ( 'SEQ1,2', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'FORMOM_SIDS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'PRESS_SIDS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'GRAV_SIDS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'LOAD_SIDS, LOAD_FACS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'MPC_SIDS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'MPCADD_SIDS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'SPC_SIDS, SPC1_SIDS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'SPCADD_SIDS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'ETYPE, EDAT, EPNT', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'ELEM PROPERTIES AND MATERIALS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'VVEC, OFFSETS, PLATE stuff', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'GRID, RGRID', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'CORD, RCORD', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'CONM2, RCONM2', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'CMASS, PMASS, RPMASS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'RFORCE_SIDS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'SLOAD_SIDS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'SNORM, RSNORM', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'RIGID_ELEM_IDS', SUBR_NAME )

         CALL LINK_MESSAGE('READ BULK DATA DECK                         ')
         CALL LOADB

         CALL DEALLOCATE_MODEL_STUF ( 'SPC_SIDS, SPC1_SIDS' )
         CALL DEALLOCATE_MODEL_STUF ( 'MPC_SIDS' )
         CALL DEALLOCATE_MODEL_STUF ( 'FORMOM_SIDS' )
         CALL DEALLOCATE_MODEL_STUF ( 'GRAV_SIDS' )
         CALL DEALLOCATE_MODEL_STUF ( 'MPCADD_SIDS' )
         CALL DEALLOCATE_MODEL_STUF ( 'RFORCE_SIDS' )
         CALL DEALLOCATE_MODEL_STUF ( 'SLOAD_SIDS' )
         CALL DEALLOCATE_MODEL_STUF ( 'SPC_MPC_SET' )
         WRITE(F06,*)

         ! Write EDAT table, if requested
         IF (DEBUG(108) > 0) THEN
            CALL WRITE_EDAT
         ENDIF

         ! Make sure we have a finite WTMASS
         ! (need to divide by it in several places later in the code)
         IF (DABS(WTMASS) < EPSIL(1)) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1804) WTMASS
            WRITE(F06,1804) WTMASS
         ENDIF

         ! Quit if FATAL_ERR > 0
         IF (FATAL_ERR > 0) THEN
            CALL OUTA_HERE ( 'Y' )
         ENDIF

! Reset MEFMGRID to GRDPNT if MEFMLOC was 'GRDPNT'. However, we need to make sure that GRDPNT is not the default (-1) otherwise
! we get problems in subr RB_DISP_MATRIX_PROC. Therefore reset GRDPNT if it is -1 to 0 (basic origin). This will mean that, if
! there is no GRDPNT in the Bulk data, we will still get results the same as if PARAM, GRDPNT, 0 was there.

         IF (MEFMLOC == 'GRDPNT') THEN
            IF (GRDPNT == -1) THEN                         !  We can't have GRDPNT = -1 if we want MEFMGRID = GRDPNT
               GRDPNT   = 0
            ENDIF
            MEFMGRID = GRDPNT
         ENDIF

         ! Allocate file names, unit numbers, etc for the
         ! USERIN elements whose matrices will be on In4 files.
         ! First check that there is a 1 to 1 correspondence of
         ! number of CUSERIN and PUSERIN B.D. entries
         IF (NPUSERIN /= NCUSERIN) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1802) NCUSERIN, NPUSERIN
            WRITE(F06,1802) NCUSERIN, NPUSERIN
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         ! If ENFORCED = 'Y' then all DOF's are SE and we will
         ! read their values from ENFFIL and use that to create LINK1O
         IF (ENFORCED == 'Y') THEN
            CALL WRITE_ENF_TO_L1O
         ENDIF

         ! Make sure that there is at least 1 grid
         IF (NGRID < 1) THEN
            WRITE(ERR,1827) NDOFG
            WRITE(F06,1827) NDOFG
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                         ! No G-set, so quit
         ENDIF

         ! Make sure that there is at least 1 element
         IF (NELE < 1) THEN
            WRITE(ERR,1820) NELE
            WRITE(F06,1820) NELE
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                         ! No G-set, so quit
         ENDIF

         ! Close files opened for writing Bulk data info to
         IF (NRIGEL > 0) THEN
            CALL FILE_CLOSE ( L1F, LINK1F, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1F, LINK1F, 'DELETE' )
         ENDIF

         IF (NFORCE > 0) THEN
            CALL FILE_CLOSE ( L1I, LINK1I, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1I, LINK1I, 'DELETE' )
         ENDIF

         IF (NTCARD > 0) THEN
            CALL FILE_CLOSE ( L1K, LINK1K, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1K, LINK1K, 'DELETE' )
         ENDIF

         IF (NAOCARD > 0) THEN
            CALL FILE_CLOSE ( L1N, LINK1N, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1N, LINK1N, 'DELETE' )
         ENDIF

         IF ((NSPC > 0) .OR. (NSPC1 > 0)) THEN
            CALL FILE_CLOSE ( L1O, LINK1O, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1O, LINK1O, 'DELETE' )
         ENDIF

         IF (NGRAV > 0) THEN
            CALL FILE_CLOSE ( L1P, LINK1P, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1P, LINK1P, 'DELETE' )
         ENDIF

         IF (NPCARD > 0) THEN
            CALL FILE_CLOSE ( L1Q, LINK1O, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1Q, LINK1O, 'DELETE' )
         ENDIF

         IF (NMPC > 0) THEN
            CALL FILE_CLOSE ( L1S, LINK1S, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1S, LINK1S, 'DELETE' )
         ENDIF

         IF (NUM_SUPT_CARDS > 0) THEN
            CALL FILE_CLOSE ( L1T, LINK1T, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1T, LINK1T, 'DELETE' )
         ENDIF

         IF (NRFORCE > 0) THEN
            CALL FILE_CLOSE ( L1U, LINK1U, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1U, LINK1U, 'DELETE' )
         ENDIF

         IF (NSLOAD > 0) THEN
            CALL FILE_CLOSE ( L1W, LINK1W, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1W, LINK1W, 'DELETE' )
         ENDIF

         IF (NUM_PARTVEC_RECORDS > 0) THEN
            CALL FILE_CLOSE ( L1V, LINK1V, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1V, LINK1V, 'DELETE' )
         ENDIF

         IF (NUM_USET_RECORDS > 0) THEN
            CALL FILE_CLOSE ( L1X, LINK1X, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( L1X, LINK1X, 'DELETE' )
         ENDIF

         CALL SET_SPARSE_MAT_SYM                           ! Set sparse matrix sym

      ENDIF res13

! Check error flag; if > 0 quit.

      IF (FATAL_ERR > 0)  THEN
          WRITE(ERR,9998) FATAL_ERR
          WRITE(F06,9998) FATAL_ERR
          CALL OUTA_HERE ( 'Y' )                           ! Finished processing data deck and there are input errors
      ENDIF

res14:IF (RESTART == 'N') THEN

!         CALL LINK_MESSAGE('ALLOCATE MEMORY FOR SOME MODEL DATA ARRAYS  ')
         CALL ALLOCATE_MODEL_STUF ( 'ESORT1', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'ESORT2', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'EOFF', SUBR_NAME )

         CALL LINK_MESSAGE('SORT AND CHECK ELEMENTS                     ')
         CALL ELESORT
         CALL DEALLOCATE_MODEL_STUF ( 'RIGID_ELEM_IDS' )

         ! Open L1B, OP2 for writing grid and coord data to.
         CALL FILE_OPEN ( L1B, LINK1B, OUNT, 'REPLACE', L1B_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

         CALL FILE_OPEN ( OP2, OP2FIL, OUNT, 'REPLACE', OP2_MSG, 'NEITHER', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         OP2STAT = 'KEEP    '
         CALL WRITE_OP2_HEADER(POST)


         ! Element processing to convert external PID's to internal.
         CALL LINK_MESSAGE('ELEM PROCESSOR                              ')
         CALL ELEM_PROP_MATL_IIDS

         ! Grid and coordinate system processing
!         CALL LINK_MESSAGE('ALLOCATE MEMORY FOR SOME MODEL DATA ARRAYS  ')
         CALL ALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'GRID_ID', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'GRID_SEQ, INV_GRID_SEQ', SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'TN', SUBR_NAME )

         CALL LINK_MESSAGE('GRID PROCESSOR                              ')
         CALL GRID_PROC
         CALL DEALLOCATE_MODEL_STUF ( 'TN' )

! Run through elements and make sure all grids on elem connection entries are defined. Need this so that if any are missing we
! won't send wrong NGRID value to Bandit. We don't need BGRID, but GET_ELEM_AGRID_BGRID checks that all grids on connection
! entries are defined

         IF (CHKGRDS == 'Y') THEN
            CALL LINK_MESSAGE('CHECK THAT ALL GRIDS FOR ELEMS EXIST        ')

            CALL COUNTER_INIT('       Process element', NELE)
            DO I=1,NELE
               CALL GET_ELEM_AGRID_BGRID ( I, 'Y' )
               CALL COUNTER_PROGRESS(I)
            ENDDO
         ENDIF
         WRITE(SC1,*) CR13

         ! Run Bandit, if requested
         REWIND (IN1)
         IF (GRIDSEQ(1:6) == 'BANDIT') THEN
            CALL LINK_MESSAGE('BANDIT - RESEQUENCE GRIDS                   ')
            CALL FILE_OPEN ( SEQ, SEQFIL, OUNT, 'REPLACE', SEQ_MSG, 'WRITE_STIME', 'FORMATTED', 'READWRITE', 'REWIND','Y','N' )
            CALL BANDIT ( NGRID, BANDIT_BW, KMAT_DEN, BANDIT_ERR )
            INQUIRE ( FILE=SEQFIL, OPENED=FILE_OPND )
            IF (FILE_OPND) THEN
               CALL FILE_CLOSE ( SEQ, SEQFIL, 'KEEP' )
            ENDIF
            WRITE(ERR,1030) BANDIT_ERR
            IF (SUPINFO == 'N') THEN
               WRITE(F06,1030) BANDIT_ERR
            ENDIF
            IF (BANDIT_ERR == 0) THEN
               KMAT_BW = 6*BANDIT_BW
            ELSE                                           ! If BANDIT had error can't use SETLKT = 1 (which depends on BANDIT_BW)
               IF (SETLKTK == 1) THEN
                  WRITE(F06,1001) 'SETLKTK', SETLKTK, SETLKTK_DEF
                  SETLKTK = SETLKTK_DEF
               ENDIF
               IF (SETLKTM == 1) THEN
                  WRITE(F06,1001) 'SETLKTM', SETLKTM, SETLKTM_DEF
                  SETLKTM = SETLKTM_DEF
               ENDIF
            ENDIF
         ENDIF
         CALL FILE_CLOSE ( IN1, INFILE, 'KEEP' )

         ! Grid point sequencing
         CALL LINK_MESSAGE('GRID SEQUENCE PROCESSOR                     ')
         CALL SEQ_PROC
         CALL FILE_CLOSE ( L1B, LINK1B, 'KEEP' )
         CALL DEALLOCATE_MODEL_STUF ( 'SEQ1,2' )

         ! Surface normal processing
         CALL ALLOCATE_MODEL_STUF ( 'GRID_SNORM', SUBR_NAME )
         DO I=1,NSNORM
            ! Find the internal grid point number for the actual grid point number specified in the SNORM entry.
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, SNORM(I,1), BGRID )

            IF (BGRID == -1) THEN
               WRITE(ERR,1828) SNORM(I,1)
               WRITE(F06,1828) SNORM(I,1)
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF

            ! When CID is supported, transform the vector here. See how it's done reading GRID entries.
            IF (SNORM(I,2) /= 0) THEN
               WRITE(ERR,1829) SNORM(I,1)
               WRITE(F06,1829) SNORM(I,1)
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF

            ! Normalize and store in the array indexed by internal grid point number.
            GRID_SNORM(BGRID,:) = RSNORM(I,:) / DSQRT(DOT_PRODUCT(RSNORM(I,:), RSNORM(I,:)))
         ENDDO
         CALL DEALLOCATE_MODEL_STUF ( 'SNORM, RSNORM' )


      ENDIF res14

      ! Subcase processing
      CALL ALLOCATE_MODEL_STUF ( 'GROUT, ELOUT', SUBR_NAME )
      CALL ALLOCATE_MODEL_STUF ( 'ELDT', SUBR_NAME )
      CALL FILE_OPEN ( L1D, LINK1D, OUNT, 'REPLACE', L1D_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
      CALL LINK_MESSAGE('SUBCASE PROCESSOR                           ')
      CALL SUBCASE_PROC
      CALL DEALLOCATE_MODEL_STUF ( 'GROUT, ELOUT' )
      CALL DEALLOCATE_MODEL_STUF ( 'SETS ARRAYS' )
      !CALL DEALLOCATE_MODEL_STUF ( 'TITLES' )
      CALL DEALLOCATE_MODEL_STUF ( 'SC_xxxx' )
      CALL FILE_CLOSE ( L1D, LINK1D, 'KEEP' )
                                                           ! If we will be writing data to the F2i disk files open them now
      IF ( IAND(OELDT,IBIT(ELDT_F21_P_T_BIT) ) > 0) THEN
         CALL FILE_OPEN ( F21, F21FIL, OUNT, 'REPLACE', F21_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
      ENDIF

      IF ( IAND(OELDT,IBIT(ELDT_F22_ME_BIT)  ) > 0) THEN
         CALL FILE_OPEN ( F22, F22FIL, OUNT, 'REPLACE', F22_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
      ENDIF

      IF ( IAND(OELDT,IBIT(ELDT_F23_KE_BIT)  ) > 0) THEN
         CALL FILE_OPEN ( F23, F23FIL, OUNT, 'REPLACE', F23_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
      ENDIF

      IF ( IAND(OELDT,IBIT(ELDT_F24_SE_BIT)  ) > 0) THEN
         CALL FILE_OPEN ( F24, F24FIL, OUNT, 'REPLACE', F24_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
      ENDIF

      IF ( IAND(OELDT,IBIT(ELDT_F25_U_P_BIT) ) > 0) THEN
         CALL FILE_OPEN ( F25, F25FIL, OUNT, 'REPLACE', F25_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
      ENDIF

      ! Print table showing the elements connected to each grid if requested
      IF (RESTART == 'N') THEN
         IF ((ANY_GPFO_OUTPUT > 0) .OR. (PRTCONN > 0)) THEN
            CALL LINK_MESSAGE('GRID/ELEMENT CONNECTION TABLE               ')
            CALL GRID_ELEM_CONN_TABLE
            IF (ANY_GPFO_OUTPUT == 0) THEN
               CALL DEALLOCATE_MODEL_STUF ( 'GRID_ELEM_CONN_ARRAY' )
            ENDIF
         ENDIF
      ENDIF

res15:IF (RESTART == 'Y') THEN

         IF (ANY_GPFO_OUTPUT > 0) THEN
            CALL LINK_MESSAGE('GRID/ELEMENT CONNECTION TABLE               ')
            CALL GRID_ELEM_CONN_TABLE
         ENDIF

         DO I=0,MELDTS-1
            CHK_ELDT_BIT(I) = 0
         ENDDO                                             ! Note below elem mass not done here since it would be done in subr GPWG
         IF ( IAND(OELDT,IBIT(ELDT_BUG_DAT1_BIT) ) > 0) THEN   ;   CHK_ELDT_BIT( 0) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_BUG_DAT2_BIT) ) > 0) THEN   ;   CHK_ELDT_BIT( 1) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_BUG_P_T_BIT)  ) > 0) THEN   ;   CHK_ELDT_BIT( 2) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_BUG_ME_BIT)   ) > 0) THEN   ;   CHK_ELDT_BIT( 3) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_BUG_KE_BIT)   ) > 0) THEN   ;   CHK_ELDT_BIT( 4) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_BUG_SE_BIT)   ) > 0) THEN   ;   CHK_ELDT_BIT( 5) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_BUG_U_P_BIT)  ) > 0) THEN   ;   CHK_ELDT_BIT( 6) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_BUG_SHPJ_BIT) ) > 0) THEN   ;   CHK_ELDT_BIT( 7) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_BUG_BMAT_BIT) ) > 0) THEN   ;   CHK_ELDT_BIT( 8) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_BUG_BCHK_BIT) ) > 0) THEN   ;   CHK_ELDT_BIT( 9) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_F21_P_T_BIT)  ) > 0) THEN   ;   CHK_ELDT_BIT(10) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_F22_ME_BIT)   ) > 0) THEN   ;   CHK_ELDT_BIT(11) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_F23_KE_BIT)   ) > 0) THEN   ;   CHK_ELDT_BIT(12) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_F24_SE_BIT)   ) > 0) THEN   ;   CHK_ELDT_BIT(13) = 1   ;   ENDIF
         IF ( IAND(OELDT,IBIT(ELDT_F25_U_P_BIT)  ) > 0) THEN   ;   CHK_ELDT_BIT(14) = 1   ;   ENDIF
         CALL ELDT_PROC_FOR_RESTART ( CHK_ELDT_BIT )
         CALL DEALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS' )

         IF (PRTCORD > 0) THEN
            CALL ALLOCATE_MODEL_STUF ( 'TN', SUBR_NAME )
            CALL CORD_PROC
            CALL DEALLOCATE_MODEL_STUF ( 'TN' )
         ENDIF

         IF (PRTBASIC > 0) THEN
            CALL WRITE_GRID_COORDS
         ENDIF

         CALL ALLOCATE_DOF_TABLES ( 'TSET' , SUBR_NAME )
         CALL ALLOCATE_DOF_TABLES ( 'TDOF' , SUBR_NAME )
         CALL ALLOCATE_DOF_TABLES ( 'TDOFI', SUBR_NAME )
         CALL ALLOCATE_DOF_TABLES ( 'TDOF_ROW_START' , SUBR_NAME )
         CALL READ_DOF_TABLES
         CALL CALC_TDOF_ROW_START ( 'N' )
!xx      CALL ALLOCATE_MODEL_STUF ( 'GRID_SEQ, INV_GRID_SEQ', SUBR_NAME )
         IF (PRTTSET > 0) THEN
            CALL WRITE_TSET
         ENDIF
         IF (PRTDOF > 0) THEN
            CALL WRITE_TDOF ( 'After all AUTOSPC' )
         ENDIF
         CALL DEALLOCATE_MODEL_STUF ( 'SCNUM' )
         CALL DEALLOCATE_MODEL_STUF ( 'SUBLOD' )
         CALL DEALLOCATE_MODEL_STUF ( 'ETYPE, EDAT, EPNT' )
         CALL DEALLOCATE_MODEL_STUF ( 'VVEC, OFFSETS, PLATE stuff' )
         CALL DEALLOCATE_MODEL_STUF ( 'ESORT1' )
         CALL DEALLOCATE_MODEL_STUF ( 'ESORT2' )
         CALL DEALLOCATE_MODEL_STUF ( 'EOFF' )
         CALL DEALLOCATE_MODEL_STUF ( 'GROUT, ELOUT' )
         CALL DEALLOCATE_MODEL_STUF ( 'ELDT' )
         CALL DEALLOCATE_MODEL_STUF ( 'GTEMP' )
         CALL DEALLOCATE_MODEL_STUF ( 'TPNT, TDATA' )
         CALL DEALLOCATE_MODEL_STUF ( 'PPNT, PDATA, PTYPE' )
         CALL DEALLOCATE_MODEL_STUF ( 'PLOAD4_3D_DATA' )

      ENDIF res15

res16:IF (RESTART == 'N') THEN

         NUM_EMG_FATAL_ERRS = 0                            ! Process ELDATA printed output requests for ELDATA item 0, 1
         JERROR  = 0
         DO I=1,NELE
            DO J=1,5                                       ! 05/27/07: Initialization of OPT added
               OPT(J) = 'N'
            ENDDO
            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_DAT1_BIT))
            I2 = IAND(ELDT(I),IBIT(ELDT_BUG_DAT2_BIT))
            IF ((I1 > 0) .OR. (I2 > 0)) THEN
               IF (I1 > 0) THEN
                  WRT_BUG(0) = 1
               ENDIF
               IF (I2 > 0) THEN
                  WRT_BUG(1) = 1
               ENDIF
               WRITE(BUG,10)
               PLY_NUM = 0
               CALL EMG ( I   , OPT, 'N', SUBR_NAME, 'Y' ) ! Calc no matrices (all OPT = 'N'). 'Y' means write to BUG file
               IF (NUM_EMG_FATAL_ERRS /=0) THEN
                  JERROR = JERROR + NUM_EMG_FATAL_ERRS
                  CYCLE
               ENDIF
            ENDIF
         ENDDO
         WRT_BUG(0) = 0
         WRT_BUG(1) = 0

         IF (JERROR > 0) THEN
            WRITE(ERR,9876) JERROR
            WRITE(F06,9876) JERROR
            CALL OUTA_HERE ( 'Y' )
         ENDIF

! Degree of Freedom processing. Also, enforced displs written to L1H.
!   Open L1F for reading rigid element data
!   Open L1S for reading MPC DOF's if NMPC > 0
!   Open L1O for reading SPC DOF's if NSPC or NSPC1 > 0
!   Open L1N for reading ASET, ASET1, OMIT, OMIT1 data if NAOCARD > 0
!   Open L1C for writing DOF tables to
!   Open L1H for writing enforced displ SPC's. Keep L1H for later processing in YS_ARRAY if SOL = STATICS

         IF (NRIGEL > 0) THEN
            CALL FILE_OPEN ( L1F, LINK1F, OUNT, 'OLD', L1F_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
         ENDIF

         IF (NMPC > 0) THEN
            CALL FILE_OPEN ( L1S, LINK1S, OUNT, 'OLD', L1S_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
         ENDIF

         IF ((NSPC > 0) .OR. (NSPC1 > 0)) THEN
            CALL FILE_OPEN ( L1O, LINK1O, OUNT, 'OLD', L1O_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
         ENDIF

         IF (NAOCARD > 0) THEN
            CALL FILE_OPEN ( L1N, LINK1N, OUNT, 'OLD', L1N_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
         ENDIF

         IF (NUM_SUPT_CARDS > 0) THEN
            CALL FILE_OPEN ( L1T, LINK1T, OUNT, 'OLD', L1T_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
         ENDIF

         CALL FILE_OPEN ( L1C, LINK1C, OUNT, 'REPLACE', L1C_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

         CALL FILE_OPEN ( L1H, LINK1H, OUNT, 'REPLACE', L1H_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

!         CALL LINK_MESSAGE('ALLOCATE MEMORY FOR SOME MODEL DATA ARRAYS  ')
         CALL ALLOCATE_DOF_TABLES ( 'TSET', SUBR_NAME )
         CALL ALLOCATE_DOF_TABLES ( 'TDOF', SUBR_NAME )
         CALL ALLOCATE_DOF_TABLES ( 'TDOFI', SUBR_NAME )
         CALL ALLOCATE_DOF_TABLES ( 'TDOF_ROW_START' , SUBR_NAME )
         CALL ALLOCATE_MODEL_STUF ( 'MPC_IND_GRIDS', SUBR_NAME )

         CALL LINK_MESSAGE('DOF PROCESSOR                               ')
         TDOF_MSG(1:)  = ' '
         TDOF_MSG(56:) = '(Before any AUTOSPC)'
         CALL DOF_PROC ( TDOF_MSG )

         ! Call subr to write any possible USETSTR output requests
         CALL WRITE_USETSTR

         CALL DEALLOCATE_MODEL_STUF ( 'SPCADD_SIDS' )

         IF (NMPC > 0) THEN
            CALL FILE_CLOSE ( L1S, LINK1S, 'KEEP' )
         ENDIF

         IF (NRIGEL > 0) THEN
            CALL FILE_CLOSE ( L1F, LINK1F, 'KEEP' )
         ENDIF

         IF ((NSPC > 0) .OR. (NSPC1 > 0)) THEN
            CALL FILE_CLOSE ( L1O, LINK1O, L1OSTAT )
         ENDIF

         IF (NDOFSE > 0) THEN
            IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                               &
               ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
               CALL FILE_CLOSE ( L1H, LINK1H, 'KEEP' )
            ELSE
               CALL FILE_CLOSE ( L1H, LINK1H, 'DELETE' )
            ENDIF
         ENDIF

         IF (NAOCARD > 0) THEN
            CALL FILE_CLOSE ( L1N, LINK1N, L1NSTAT )
         ENDIF

         ! CONM2 processing to get CONM2 data in basic coords at the mass
         CALL FILE_OPEN ( L1Y, LINK1Y, OUNT, 'REPLACE', L1Y_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
         CALL LINK_MESSAGE('CONM2 PROCESSOR #1                          ')
         CALL CONM2_PROC_1
         CALL FILE_CLOSE ( L1Y, LINK1Y, L1YSTAT )

         ! Grid point weight generator (model weight, c.g., etc.)
         CALL LINK_MESSAGE('ALLOCATE MEMORY FOR RBGLOBAL ARRAY          ')
         CALL ALLOCATE_RBGLOBAL ( 'G ', SUBR_NAME )
         RBG_GSET_ALLOCATED = 'Y'

         IF ((GRDPNT_IN >= 0) .OR. (MEFFMASS_CALC == 'Y')) THEN

            CALL LINK_MESSAGE('GRID POINT WEIGHT GENERATOR              ')

            IF (NCUSERIN > 0) THEN
               DO I=1,NELE
                  IF (ETYPE(I) == 'USERIN  ') THEN
                     CALL GPWG_USERIN ( I )
                  ENDIF
               ENDDO
            ENDIF

            IF (NCUSERIN > 0) THEN
               CALL GPWG ( 'RESIDUAL STR' )
            ENDIF

            CALL GPWG ( 'OA MODEL    ' )
            WRITE(SC1,*) CR13

         ENDIF

      ENDIF res16

res17:IF (RESTART == 'Y') THEN

         CALL DEALLOCATE_MODEL_STUF ( 'CMASS, PMASS, RPMASS' )
         CALL DEALLOCATE_MODEL_STUF ( 'ELEM PROPERTIES AND MATERIALS' )
         CLOSE ( L1B, STATUS='KEEP')
         CLOSE ( L1G, STATUS='KEEP')
         CLOSE ( L1K, STATUS='KEEP')
         CLOSE ( L1Y, STATUS='KEEP')

         CALL LINK_MESSAGE('CHECKING FOR MATRICES TO PRINT ON RESTART')
         IF ((PRTSTIFF(1) >= 1) .OR. (PRTSTIFD(1) >= 1) .OR. (EQCHK_OUTPUT(1) > 0)) THEN
            CALL ALLOCATE_SPARSE_MAT ( 'KGG', NDOFG, NTERM_KGG, SUBR_NAME )
            CALL READ_MATRIX_1 ( LINK1L, L1L, 'N', 'Y', L1LSTAT, L1L_MSG,'KGG', NTERM_KGG, 'Y', NDOFG, I_KGG, J_KGG, KGG)
         ENDIF
         CALL PRT_MATS_ON_RESTART
         WRITE(SC1, * ) '       Deallocate KGG          '
         CALL DEALLOCATE_SPARSE_MAT ( 'KGG' )

         IF (EQCHK_OUTPUT(1) > 0) THEN                     ! Do equilibrium check on the G-set stiffness matrix, if requested

            IF (RBG_GSET_ALLOCATED == 'N') THEN
               CALL ALLOCATE_RBGLOBAL ( 'G ', SUBR_NAME )
            ENDIF
            CALL RB_DISP_MATRIX_PROC ('EQCHK REF GRID', 0 )! Need rigid body displ matrix for the G-set if equil check requested
                                                           ! Get matrix diagonal and stats and print, if requested.
            IF ((PRTSTIFF(1) >= 1) .OR. (PRTSTIFD(1) >= 1)) THEN
               PRINTIT = 0                                 ! Diag stats were already printed in subr PRT_MATS_ON_RESTART
            ELSE
               PRINTIT = 1
            ENDIF
            CALL GET_MATRIX_DIAG_STATS ( 'KGG', 'G ', NDOFG, NTERM_KGG, I_KGG, J_KGG, KGG, PRINTIT, KGG_DIAG, KGG_MAX_DIAG  )

            CALL LINK_MESSAGE('EQUILIBRIUM CHECK ON KGG                  ')
            CALL STIFF_MAT_EQUIL_CHK ( EQCHK_OUTPUT(1), 'G ', SYM_KGG, NDOFG, NTERM_KGG, I_KGG, J_KGG, KGG, KGG_DIAG, KGG_MAX_DIAG,&
                                       RBGLOBAL_GSET )

         ENDIF

      ENDIF res17

res18:IF (RESTART == 'N') THEN

         ! CONM2 processing to get CONM2 data in global coords at the grid
         CALL LINK_MESSAGE('CONM2 PROCESSOR #2                          ')
         CALL CONM2_PROC_2

         ! Temperature data processing
         IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                                  &
            ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
            CALL LINK_MESSAGE('ALLOCATE MEM FOR ELEM AND GRID TEMP ARRAYS    ')
            CALL ALLOCATE_MODEL_STUF ( 'GTEMP', SUBR_NAME )
            CALL ALLOCATE_MODEL_STUF ( 'CGTEMP', SUBR_NAME )
            CALL ALLOCATE_MODEL_STUF ( 'ETEMP', SUBR_NAME )
            CALL ALLOCATE_MODEL_STUF ( 'CETEMP', SUBR_NAME )
            CALL ALLOCATE_MODEL_STUF ( 'TPNT, TDATA', SUBR_NAME )

            IF (NTCARD > 0) THEN
               OUNT(1) = ERR
               OUNT(2) = F06
               CALL FILE_OPEN ( L1K, LINK1K, OUNT, 'OLD', L1K_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
               CALL LINK_MESSAGE('GRID & ELEM TEMPERATURE DATA PROCESSOR      ')
               CALL TEMPERATURE_DATA_PROC
               CALL FILE_CLOSE ( L1K, LINK1K, 'KEEP' )
            ELSE
               CALL FILE_CLOSE ( L1K, LINK1K, 'DELETE' )
            ENDIF
         ENDIF
         CALL DEALLOCATE_MODEL_STUF ( 'ETEMP' )
         CALL DEALLOCATE_MODEL_STUF ( 'CGTEMP' )
         CALL DEALLOCATE_MODEL_STUF ( 'CETEMP' )

         ! Element pressure data processing.
         ! Open L1Q which contains element pressure Bulk Data
         IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                                  &
            ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
            CALL LINK_MESSAGE('ALLOCATE MEMORY FOR ELEM PRESSURE DATA ARRAYS ')
            CALL ALLOCATE_MODEL_STUF ( 'PPNT, PDATA, PTYPE', SUBR_NAME )
            CALL ALLOCATE_MODEL_STUF ( 'PLOAD4_3D_DATA', SUBR_NAME )

            IF (NPCARD > 0) THEN

               CALL FILE_OPEN ( L1Q, LINK1Q, OUNT, 'OLD', L1Q_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

               CALL LINK_MESSAGE('ELEMENT PRESSURE DATA PROCESSOR             ')
               CALL PRESSURE_DATA_PROC

               CALL FILE_CLOSE ( L1Q, LINK1Q, 'KEEP' )
            ENDIF
         ENDIF
         CALL DEALLOCATE_MODEL_STUF ( 'PRESS_SIDS' )

! Generate TR6_CG rigid body displ matrix for the R-set if this run has R-set DOF's. TR6_CG not needed until LINK6 but, by calc'ing
! it here we can avoid keeping CORD, RCORD, GRID, RGRID arrays from having to be allocated until then

         IF (NDOFR > 0) THEN

            CALL LINK_MESSAGE('GENERATE RIGID BODY DISPL MATRIX            ')

            CALL RB_DISP_MATRIX_PROC ( 'CG', 0 )
            CALL ALLOCATE_RBGLOBAL ( 'R ', SUBR_NAME )
            CALL TDOF_COL_NUM ( 'R ', R_SET_COL )

            DO I=1,NDOFG
               R_SET_DOF = TDOFI(I,R_SET_COL)
               IF (R_SET_DOF > 0) THEN
                  DO J=1,6
                     TR6_CG(R_SET_DOF,J) = RBGLOBAL_GSET(I,J)
                  ENDDO
               ENDIF
            ENDDO

            IF (MEFMLOC == 'CG  ') THEN
               DO I=1,NDOFG
                  R_SET_DOF = TDOFI(I,R_SET_COL)
                  IF (R_SET_DOF > 0) THEN
                     DO J=1,6
                        TR6_MEFM(R_SET_DOF,J) = TR6_CG(R_SET_DOF,J)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               IF (MEFMGRID == 0) THEN
                  CALL RB_DISP_MATRIX_PROC ( 'BASIC ORIGIN', 0 )
               ELSE
                  CALL RB_DISP_MATRIX_PROC ( 'GRID', MEFMGRID )
               ENDIF
               DO I=1,NDOFG
                  R_SET_DOF = TDOFI(I,R_SET_COL)
                  IF (R_SET_DOF > 0) THEN
                     DO J=1,6
                        TR6_MEFM(R_SET_DOF,J) = RBGLOBAL_GSET(I,J)
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF

            CALL RB_DISP_MATRIX_PROC ( 'BASIC ORIGIN', 0 )
            DO I=1,NDOFG
               R_SET_DOF = TDOFI(I,R_SET_COL)
               IF (R_SET_DOF > 0) THEN
                  DO J=1,6
                     TR6_0(R_SET_DOF,J) = RBGLOBAL_GSET(I,J)
                  ENDDO
               ENDIF
            ENDDO

         ENDIF

         ! Generate rigid body displ matrix for the G-set
         ! if requested by user (if it has been, some EQCHK_OUTPUT(i) will be > 0)
         IF((EQCHK_OUTPUT(1) > 0) .OR. (EQCHK_OUTPUT(2) > 0) .OR. (EQCHK_OUTPUT(3) > 0) .OR. (EQCHK_OUTPUT(4) > 0) .OR.            &
            (EQCHK_OUTPUT(5) > 0)) THEN
            IF (RBG_GSET_ALLOCATED == 'N') THEN
               CALL LINK_MESSAGE('ALLOCATE MEMORY FOR RBGLOBAL ARRAY          ')
               CALL ALLOCATE_RBGLOBAL ( 'G ', SUBR_NAME )
            ENDIF
                                                        ! already calculated above
            CALL LINK_MESSAGE('GENERATE RIGID BODY DISPL MATRIX            ')
            CALL RB_DISP_MATRIX_PROC ( 'EQCHK REF GRID', 0 )
         ENDIF

         ! Calculate YSe enforced displ array
         IF (NDOFSE > 0) THEN
            IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                               &
               ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN

               CALL LINK_MESSAGE('ALLOCATE MEMORY FOR Yse ARRAY               ')
               CALL ALLOCATE_COL_VEC ( 'YSe', NDOFSE, SUBR_NAME )

               OUNT(1) = ERR
               OUNT(2) = F06
               CALL FILE_OPEN ( L1H, LINK1H, OUNT, 'OLD', L1H_MSG, 'READ_STIME', 'UNFORMATTED', 'READWRITE', 'REWIND',             &
                               'Y', 'N' )
               CALL LINK_MESSAGE('CALCULATE YS ENFORCED DISPL ARRAY           ')
               CALL YS_ARRAY
               CALL DEALLOCATE_COL_VEC ( 'YSe' )
               CALL FILE_CLOSE ( L1H, LINK1H, 'KEEP' )

            ENDIF

         ENDIF

         ! Write problem size parameters
         WRITE(ERR,144) NGRID
         WRITE(ERR,145) NDOFG
         IF (SUPINFO == 'N') THEN
            WRITE(F06,144) NGRID
            WRITE(F06,145) NDOFG
         ENDIF
         IF (DEBUG(187) > 0) THEN
            CALL WRITE_ELEM_SUMMARY
         ENDIF
         WRITE(ERR,146) NELE
         IF (SUPINFO == 'N') THEN
            WRITE(F06,146) NELE
         ENDIF

         ! Write data to L1Z for restart
         !
         ! only need it if CHKPNT was in CC but write it anyways
         ! in case user changes their mind
         CALL WRITE_L1Z

      ENDIF res18

! Deallocate several items and we will re-allocate them at the beginning of LINK1. This is needed for BUCKLING and NLSTATIC
! since LINK0 does not get re-run at each step but LINK1 does (and we need these arrays allocated each time LINK1 starts.
! Earlier, we had 'SINGLE ELEMENT ARRAYS' allocated here (but not in LINK1).

      CALL DEALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS' )

      ! Write messages on what OUTPUT4 matrices were requested and open the OU4 files
      IF (NUM_OU4_REQUESTS > 0) THEN

         CALL OUTPUT4_MATRIX_MSGS ( OUNT )                 ! Write the messages on what OUTPUT4 matrices were requested

         DO I=1,MOU4                                       ! Determine which OU4 units will be needed
            OPEN_UNT(I) = 'N'
            DO J=1,NUM_OU4_REQUESTS
               IF (OU4_FILE_UNITS(J) == OU4(I)) THEN
                  OPEN_UNT(I) = 'Y'
               ENDIF
            ENDDO
         ENDDO

         DO I=1,MOU4                                       ! Open the OU4 units
            IF (OPEN_UNT(I) == 'Y') THEN
               CALL FILE_OPEN (OU4(I), OU4FIL(I), OUNT, 'REPLACE', OU4_MSG(I),'NEITHER','UNFORMATTED','WRITE','REWIND','Y','N')
            ENDIF
         ENDDO

      ENDIF

      ! Write data to L1A
      !
      ! rewrites L1A if this a restart - needed to tell LINK9 that this is a restart
      COMM(LINKNO) = 'C'
res20:IF (RESTART == 'N') THEN
         CALL WRITE_L1A ( 'KEEP', 'Y' )
      ENDIF res20

      ! Check allocation status of allocatable arrays, if requested
      IF (DEBUG(100) > 0) THEN
         CALL CHK_ARRAY_ALLOC_STAT
         IF (DEBUG(100) > 1) THEN
            CALL WRITE_ALLOC_MEM_TABLE ( 'at the end of '//SUBR_NAME )
         ENDIF
      ENDIF

      ! Write LINK0 end to F06
      CALL OURTIM
      WRITE(F06,151) LINKNO

      IF (( DEBUG(193) == 1) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'Near end of LINK0' )
      ENDIF

      ! Write LINK0 end to screen
      WRITE(SC1,154) LINKNO

! **********************************************************************************************************************************
   10 FORMAT(' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',&
             '>>>>>>>>>>>>>>>>>')

   11 FORMAT(/,' *INFORMATION: KIND PARAMETERS ARE: CHARACTER BYTE                                  = ',I4,' BYTES',/ &
               '                                    SHORT INTEGER                                   = ',I4,' BYTES',/ &
               '                                    LONG  INTEGER                                   = ',I4,' BYTES',/ &
               '                                    SINGLE PREC REAL                                = ',I4,' BYTES',/ &
               '                                    DOUBLE PREC REAL                                = ',I4,' BYTES',/ &
               '                                    QUAD   PREC REAL                                = ',I4,' BYTES')

  144 FORMAT(' *INFORMATION: NUMBER OF GRID POINTS                                                  = ',I12)

  145 FORMAT(' *INFORMATION: NUMBER OF G SET DEGREES OF FREEDOM (NDOFG)                             = ',I12,/)

  146 FORMAT(' *INFORMATION: TOTAL NUMBER OF ELEMENTS (NELE)                                        = ',I12,/)

  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END',/)

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  154 FORMAT(  ' >> LINK',I3,' END')

 1001 FORMAT(' *INFORMATION: BANDIT DID NOT RUN SUCCESSFULLY SO USER REQUEST FOR ',A,' = ',I2,' CANNOT BE USED.' ,                 &
                           ' IT IS RESET TO DEFAULT VALUE = ',I2)

 1030 FORMAT(' *INFORMATION: BANDIT WAS CALLED TO RESEQUENCE THE GRIDS AND HAS RETURNED WITH ERROR  = ',I8,/)

 1802 FORMAT(' *ERROR  1802: THERE MUST BE THE SAME NUMBER OF CUSERIN ELEM CONN ENTRIES, (NCUSERIN = ',I8,')'                      &
                    ,/,14X,'                    AND NUMBER OF PUSERIN PROPERTY  ENTRIES, (NPUSERIN = ',I8,')')

 1804 FORMAT(' *ERROR  1804: ABSOLUTE VALUE OF BULK DATA PARAM WTMASS MUST BE > 0 BUT IS = ',1ES14.6)

 1820 FORMAT(' *ERROR  1820: THERE MUST BE AT LEAST 1 ELEMENT IN THE BULK DATA. HOWEVER, NELE = ',I8)

 1827 FORMAT(' *ERROR  1827: THERE MUST BE AT LEAST 1 GRID IN THE BULK DATA. HOWEVER, NGRID = ',I8)

 1828 FORMAT(' *ERROR  1828: GRID POINT ',I8,' ON BULK DATA SNORM ENTRY DOES NOT EXIST')

 1829 FORMAT(' *ERROR  1829: CID IS NOT 0 OR BLANK ON BULK DATA SNORM ENTRY FOR GRID POINT ',I8)

 9401 FORMAT(' *WARNING    : REQUEST FOR CHKPNT IS NOT ALLOWED IN SOLUTION ',A)

 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')

 9998 FORMAT(/,' PROCESSING TERMINATED DUE TO ',I8,' INPUT ERRORS. CHECK OUTPUT FILE FOR ERROR MESSAGES')

12345 FORMAT(A, ', I = ', I8, ' of ', I8, A)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE ELDT_PROC_FOR_RESTART ( CHK_ELDT_BIT )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_FIJ
      USE SCONTR, ONLY                :  ELDT_BUG_DAT1_BIT, ELDT_BUG_DAT1_BIT, ELDT_BUG_ME_BIT  , ELDT_BUG_P_T_BIT ,               &
                                         ELDT_BUG_SE_BIT  , ELDT_BUG_KE_BIT  , ELDT_BUG_SHPJ_BIT, ELDT_BUG_BMAT_BIT,               &
                                         ELDT_BUG_BCHK_BIT, ELDT_F21_P_T_BIT , ELDT_F22_ME_BIT  , ELDT_F23_KE_BIT  ,               &
                                         ELDT_F24_SE_BIT  , MBUG, MFIJ, MELDTS

      IMPLICIT NONE

      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option flags for subr EMG (to tell it what to calc)

                                                           ! Array used to tell subr ELDT_PROC_FOR_RESTART which ELDT to calc
      INTEGER(LONG), INTENT(IN)       :: CHK_ELDT_BIT(0:MELDTS-1)

      INTEGER(LONG)                   :: II,jj             ! DO loop indices
      INTEGER(LONG)                   :: I00               ! Indicex
      INTEGER(LONG)                   :: IIERR             ! Error indicator

! **********************************************************************************************************************************
      DO II=1,6
         OPT(II) = 'N'
      ENDDO

      DO II=0,MBUG-1
         WRT_BUG(II) = 0
      ENDDO

      IF (CHK_ELDT_BIT( 0) > 0) THEN
         WRT_BUG(0) = 1
      ENDIF

      IF (CHK_ELDT_BIT( 1) > 0) THEN
         WRT_BUG(1) = 1
      ENDIF

      IF (CHK_ELDT_BIT( 2) > 0) THEN
         OPT(1)     = 'Y'
         WRT_BUG(2) = 1
      ENDIF

      IF (CHK_ELDT_BIT( 3) > 0) THEN
         OPT(2)     = 'Y'
         OPT(5)     = 'Y'
         WRT_BUG(3) = 1
      ENDIF

      IF (CHK_ELDT_BIT( 4) > 0) THEN
         OPT(3)     = 'Y'
         WRT_BUG(4) = 1
      ENDIF

      IF (CHK_ELDT_BIT( 5) > 0) THEN
         OPT(4)     = 'Y'
         WRT_BUG(5) = 1
      ENDIF

      IF (CHK_ELDT_BIT( 6) > 0) THEN
         OPT(5)     = 'Y'
         WRT_BUG(6) = 1
      ENDIF

      IF (CHK_ELDT_BIT( 7) > 0) THEN
         OPT(5)     = 'Y'
         WRT_BUG(7) = 1
      ENDIF

      IF (CHK_ELDT_BIT( 8) > 0) THEN
         OPT(5)     = 'Y'
         WRT_BUG(8) = 1
      ENDIF

      IF (CHK_ELDT_BIT( 9) > 0) THEN
         OPT(5)     = 'Y'
         WRT_BUG(9) = 1
      ENDIF

      NUM_EMG_FATAL_ERRS = 0
      IIERR = 0
      I00   = 0
      DO II=1,NELE

         I00 = 0

         DO JJ=0,MBUG-1
            WRT_BUG(JJ) = 0
         ENDDO

         DO JJ=1,MFIJ
            WRT_FIJ(JJ) = 0
         ENDDO

         IF ( IAND(ELDT(II),IBIT(ELDT_BUG_DAT1_BIT) ) > 0) THEN   ;   I00 = 1   ;   WRT_BUG(0) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_BUG_DAT2_BIT) ) > 0) THEN   ;   I00 = 1   ;   WRT_BUG(1) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_BUG_P_T_BIT)  ) > 0) THEN   ;   I00 = 1   ;   WRT_BUG(2) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_BUG_ME_BIT)   ) > 0) THEN   ;   I00 = 1   ;   WRT_BUG(3) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_BUG_KE_BIT)   ) > 0) THEN   ;   I00 = 1   ;   WRT_BUG(4) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_BUG_SE_BIT)   ) > 0) THEN   ;   I00 = 1   ;   WRT_BUG(5) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_BUG_U_P_BIT)  ) > 0) THEN   ;   I00 = 1   ;   WRT_BUG(6) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_BUG_SHPJ_BIT) ) > 0) THEN   ;   I00 = 1   ;   WRT_BUG(7) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_BUG_BMAT_BIT) ) > 0) THEN   ;   I00 = 1   ;   WRT_BUG(8) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_BUG_BCHK_BIT) ) > 0) THEN   ;   I00 = 1   ;   WRT_BUG(9) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_F21_P_T_BIT)  ) > 0) THEN   ;   I00 = 1   ;   WRT_FIJ(1) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_F22_ME_BIT)   ) > 0) THEN   ;   I00 = 1   ;   WRT_FIJ(2) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_F23_KE_BIT)   ) > 0) THEN   ;   I00 = 1   ;   WRT_FIJ(3) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_F24_SE_BIT)   ) > 0) THEN   ;   I00 = 1   ;   WRT_FIJ(4) = 1   ;   ENDIF
         IF ( IAND(ELDT(II),IBIT(ELDT_F25_U_P_BIT)  ) > 0) THEN   ;   I00 = 1   ;   WRT_FIJ(5) = 1   ;   ENDIF

         IF (I00 > 0) THEN

            WRITE(BUG,10)
            PLY_NUM = 0
            CALL EMG ( II  , OPT, 'N', SUBR_NAME, 'Y' )    ! 'Y' means write to BUG file

            IF (NUM_EMG_FATAL_ERRS /=0) THEN
               IIERR = IIERR + NUM_EMG_FATAL_ERRS
               CYCLE
            ENDIF

            IF ((WRT_FIJ(1) > 0) .OR. (WRT_FIJ(2) > 0) .OR. (WRT_FIJ(3) > 0) .OR. (WRT_FIJ(4) > 0) .OR. (WRT_FIJ(5) > 0)) THEN
               CALL WRITE_FIJFIL ( 9999, 0 )
            ENDIF

         ENDIF

      ENDDO

      IF (IIERR > 0) THEN
         WRITE(ERR,9876) IIERR
         WRITE(F06,9876) IIERR
         CALL OUTA_HERE ( 'Y' )                         ! Errors from subr EMG, so quit
      ENDIF

! **********************************************************************************************************************************
   10 FORMAT(' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',&
             '>>>>>>>>>>>>>>>>>')

 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')

! **********************************************************************************************************************************

      END SUBROUTINE ELDT_PROC_FOR_RESTART

! ##################################################################################################################################

      SUBROUTINE WRITE_ELEM_SUMMARY

      USE SCONTR, ONLY                :  NCBAR, NCBEAM, NCELAS1, NCELAS2, NCELAS3, NCELAS4, NCHEXA8, NCHEXA20, NCPENTA6, NCPENTA15,&
                                         NCQUAD4, NCQUAD4K, NCQUAD8, NCROD, NCTETRA4, NCTETRA10, NCTRIA3, NCTRIA3K, NCUSER1,       &
                                         NCUSERIN

      INTEGER(LONG)                   :: TOTAL = 0         ! Sum of number of elements written

! **********************************************************************************************************************************
      WRITE(F06,7776)

      IF (NCBAR     > 0) THEN
         TOTAL = TOTAL + NCBAR
         WRITE(F06,7777) 'BAR    ', NCBAR
      ENDIF

      IF (NCBEAM    > 0) THEN
         TOTAL = TOTAL + NCBEAM
         WRITE(F06,7777) 'BEAM   ', NCBEAM
      ENDIF

      IF (NCELAS1   > 0) THEN
         TOTAL = TOTAL + NCELAS1
         WRITE(F06,7777) 'CELAS1 ', NCELAS1
      ENDIF

      IF (NCELAS2   > 0) THEN
         TOTAL = TOTAL + NCELAS2
         WRITE(F06,7777) 'CELAS2 ', NCELAS2
      ENDIF

      IF (NCELAS3   > 0) THEN
         TOTAL = TOTAL + NCELAS3
         WRITE(F06,7777) 'CELAS3 ', NCELAS3
      ENDIF

      IF (NCELAS4   > 0) THEN
         TOTAL = TOTAL + NCELAS4
         WRITE(F06,7777) 'CELAS4 ', NCELAS4
      ENDIF

      IF (NCHEXA8   > 0) THEN
         TOTAL = TOTAL + NCHEXA8
         WRITE(F06,7777) 'HEXA8  ', NCHEXA8
      ENDIF

      IF (NCHEXA20  > 0) THEN
         TOTAL = TOTAL + NCHEXA20
         WRITE(F06,7777) 'HEXA20 ', NCHEXA20
      ENDIF

      IF (NCPENTA6  > 0) THEN
         TOTAL = TOTAL + NCPENTA6
         WRITE(F06,7777) 'PENTA6 ', NCPENTA6
      ENDIF

      IF (NCPENTA15 > 0) THEN
         TOTAL = TOTAL + NCPENTA15
         WRITE(F06,7777) 'PENTA15', NCPENTA15
      ENDIF

      IF (NCQUAD4   > 0) THEN
         TOTAL = TOTAL + NCQUAD4
         WRITE(F06,7777) 'QUAD4  ', NCQUAD4
      ENDIF

      IF (NCQUAD4K  > 0) THEN
         TOTAL = TOTAL + NCQUAD4K
         WRITE(F06,7777) 'QUAD4K ', NCQUAD4K
      ENDIF

      IF (NCQUAD8   > 0) THEN
         TOTAL = TOTAL + NCQUAD8
         WRITE(F06,7777) 'QUAD8  ', NCQUAD8
      ENDIF

      IF (NCROD     > 0) THEN
         TOTAL = TOTAL + NCROD
         WRITE(F06,7777) 'ROD    ', NCROD
      ENDIF

      IF (NCTETRA4  > 0) THEN
         TOTAL = TOTAL + NCTETRA4
         WRITE(F06,7777) 'TETRA4 ', NCTETRA4
      ENDIF

      IF (NCTETRA10 > 0) THEN
         TOTAL = TOTAL + NCTETRA10
         WRITE(F06,7777) 'TETRA10', NCTETRA10
      ENDIF

      IF (NCTRIA3   > 0) THEN
         TOTAL = TOTAL + NCTRIA3
         WRITE(F06,7777) 'TRIA3  ', NCTRIA3
      ENDIF

      IF (NCTRIA3K  > 0) THEN
         TOTAL = TOTAL + NCTRIA3K
         WRITE(F06,7777) 'TRIA3K ', NCTRIA3K
      ENDIF

      IF (NCUSER1   > 0) THEN
         TOTAL = TOTAL + NCUSER1
         WRITE(F06,7777) 'USER1  ', NCUSER1
      ENDIF

      IF (NCUSERIN  > 0) THEN
         TOTAL = TOTAL + NCUSERIN
         WRITE(F06,7777) 'USERIN ', NCUSERIN
      ENDIF

      IF (TOTAL /= NELE) THEN
         WRITE(F06,7779) TOTAL - NELE
      ENDIF

      WRITE(F06,7778)

! **********************************************************************************************************************************
 7776 FORMAT(' *INFORMATION: FOLLOWING IS A LIST OF ALL ELEMENTS IN THE INPUT DATA:')

 7777 FORMAT(15X,'Number of ',A,' elements                                             = ',I12)

 7778 FORMAT(15X,'                                                                         ------------')

 7779 FORMAT(15X,'Other                                                                  = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_ELEM_SUMMARY

      END SUBROUTINE LINK0


      SUBROUTINE LINK1

! LINK1 processes:

!  - MPC's and rigid elements resulting in a sparse RMG array of constraint coefficients
!  - applied forces (incl grid forces and moments, gravity, pressure, thermal, centrifugal, scalar). Result is a sparse load array
!  - G-set mass, stiffness and load matrices. Results are sparse stiffness, mass and load arrays
!  - process differential stiffness resulting in a sparse array of differential stiffness values

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG

      USE IOUNT1, ONLY                :  ERR, F06, F21, F22, F23, F24, L1C, L1F, L1I, L1G, L1J, L1P, L1S, L1U, L1W, SC1


      USE IOUNT1, ONLY                :  F21FIL, F22FIL, F23FIL, F24FIL, LINK1C, LINK1F, LINK1I, LINK1G, LINK1J, LINK1P, LINK1S,   &
                                         LINK1U, LINK1W

      USE IOUNT1, ONLY                :  L1FSTAT, L1ISTAT, L1PSTAT, L1SSTAT, L1USTAT, L1WSTAT

      USE IOUNT1, ONLY                :  F21_MSG, F22_MSG, F23_MSG, F24_MSG, L1F_MSG, L1G_MSG, L1I_MSG, L1J_MSG, L1P_MSG, L1S_MSG, &
                                         L1U_MSG, L1W_MSG

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, COMM, ELDT_F21_P_T_BIT, ELDT_F22_ME_BIT, ELDT_F23_KE_BIT, ELDT_F24_SE_BIT,  &
                                         FATAL_ERR, IBIT, LINKNO, LTERM_KGG, LTERM_KGGD, LTERM_MGGE, NDOFM, NFORCE,                &
                                         NGRAV, NMPC, NPLOAD, NRFORCE, NRIGEL, NSLOAD, NTERM_RMG, NTSUB, NUM_BUCKLING_SUBS,        &
                                         RESTART, SOL_NAME

      USE DOF_TABLES, ONLY            :  TDOFI

      USE PARAMS, ONLY                :  EMP0_PAUSE, ESP0_PAUSE, SETLKTK, SKIPMGG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  OELDT
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURDAT, OURTIM, TIME_INIT
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE, WRITE_L1A
      USE MODEL_STORAGE_ALLOCATION, ONLY:  ALLOCATE_MODEL_STUF
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_INQUIRE, FILE_OPEN
      USE MPC_PROCESSING, ONLY        :  MPC_PROC
      USE RIGID_ELEMENT_PROCESSING, ONLY:  RIGID_ELEM_PROC
      USE SPARSE_LOAD_CONSTRAINTS, ONLY:  SPARSE_PG, SPARSE_RMG
      USE FORCE_MOM_PROCESSING, ONLY  :  FORCE_MOM_PROC
      USE LINK1_LOAD_SUPPORT, ONLY    :  EPTL, SLOAD_PROC
      USE MASS_MATRIX_ASSEMBLY, ONLY  :  EMP0, EMP, MGGC_MASS_MATRIX, SPARSE_MGG
      USE LINK1_WORKSPACE_LIFECYCLE, ONLY:  ALLOCATE_EMS_ARRAYS, ALLOCATE_L1_MGG, ALLOCATE_STF_ARRAYS, DEALLOCATE_EMS_ARRAYS, DEALLOCATE_L1_MGG, DEALLOCATE_STF_ARRAYS
      USE MODEL_STORAGE_DEALLOCATION, ONLY:  DEALLOCATE_MODEL_STUF
      USE GRAV_PROCESSING, ONLY       :  GRAV_PROC
      USE RFORCE_PROCESSING, ONLY     :  RFORCE_PROC
      USE STIFFNESS_MATRIX_ASSEMBLY, ONLY:  ESP0, ESP, SPARSE_KGG, SPARSE_KGGD
      USE IN4_FILE_LIFECYCLE, ONLY    :  DEALLOCATE_IN4_FILES
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_DOF_TABLES
      USE ELEMENT_MODEL_INDEXING, ONLY:  ELSAVE
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_ARRAY_ALLOC_STAT, LINK_MESSAGE, WRITE_ALLOC_MEM_TABLE
      USE LINK1_WORKFLOW_SUPPORT, ONLY:  BUILD_KGGD_FROM_UG

      IMPLICIT NONE

      LOGICAL                         :: LEXIST            ! .TRUE. if a file exists
      LOGICAL                         :: LOPEN             ! .TRUE. if a file is opened

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK1'
      CHARACTER(10*BYTE)              :: LTERM_NAME        ! Name for an LTERM value
      CHARACTER( 1*BYTE)              :: RESPONSE          ! User response ('Y' or 'N') to a screen prompt

      INTEGER(LONG)                   :: BUCKLING_STEP     ! If SOL is BUCKLING then this is step 1 or 2 in the process, otherwise 0
      INTEGER(LONG)                   :: I1                ! Intermediate integer variable
      INTEGER(LONG)                   :: I2                ! Intermediate integer variable
      INTEGER(LONG)                   :: I3                ! Intermediate integer variable
      INTEGER(LONG)                   :: I4                ! Intermediate integer variable
      INTEGER(LONG)                   :: LTERM             ! Value for LTERM_KGGD or LTERM_KGG
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG), PARAMETER        :: P_LINKNO    = 0   ! Prior LINK no's that should have run before this LINK can execute

      INTRINSIC                       :: IAND

! **********************************************************************************************************************************
      LINKNO = 1

! Set time initializing parameters

      CALL TIME_INIT

! Get date and time, write to screen

      CALL OURDAT
      CALL OURTIM
      WRITE(SC1,152) LINKNO

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Write info to text files

      WRITE(F06,150) LINKNO
      WRITE(ERR,150) LINKNO

! Read LINK1A file

!xx   CALL READ_L1A ( 'KEEP' )
      CALL INIT_COUNTERS

! Check COMM for successful completion of prior LINKs

      IF (COMM(P_LINKNO) /= 'C') THEN
         WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
         WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Prior LINK's didn't complete, so quit
      ENDIF

! Allocate DOF tables

      CALL LINK_MESSAGE('ALLOCATE DOF TABLES')

! Set BUCKLING_STEP based on LOAD_ISTEP (see subr MYSTRAN.FOR)

      BUCKLING_STEP = 0
      IF (SOL_NAME(1:8) == 'BUCKLING') THEN
         BUCKLING_STEP = LOAD_ISTEP
      ENDIF

! Allocate several arrays since they are needed every time LINK1 gets re-run (for BUCKLING and NLSTATIC). They were used in
! LINK0 but deallocated there at the end of LINK0 so that we could allocate them again here

      CALL ALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS', SUBR_NAME )

res19:IF (RESTART == 'N') THEN

         CALL ALLOCATE_MODEL_STUF ( 'SYS_LOAD', SUBR_NAME )

! Process MPC's and rigid elements.

         IF ((NMPC > 0) .OR. (NRIGEL > 0)) THEN

            CALL FILE_OPEN ( L1J, LINK1J, OUNT, 'REPLACE', L1J_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

            IF (NMPC > 0) THEN
               CALL FILE_OPEN ( L1S, LINK1S, OUNT, 'OLD', L1S_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
               CALL LINK_MESSAGE('MPC PROCESSOR                               ')
               CALL MPC_PROC
               CALL FILE_CLOSE ( L1S, LINK1S, 'KEEP' )
            ENDIF

            IF (NRIGEL > 0) THEN                           ! Process rigid elements.
               CALL FILE_OPEN ( L1F, LINK1F, OUNT, 'OLD', L1F_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
               CALL LINK_MESSAGE('RIGID ELEMENT PROCESSOR                     ')
               CALL RIGID_ELEM_PROC
               CALL FILE_CLOSE ( L1F, LINK1F, 'KEEP' )
            ENDIF

            CALL FILE_CLOSE ( L1J, LINK1J, 'KEEP' )   ! Subr SPARSE_RMG will reopen LINK1S

            CALL LINK_MESSAGE('ALLOCATE MEMORY FOR RMG ARRAY               ')

                                                           ! Generate sparse RMG (constraint) matrix.
            CALL LINK_MESSAGE('SPARSE RMG PROCESSOR                        ')
            CALL SPARSE_RMG
            CALL FILE_CLOSE ( L1J, LINK1J, 'KEEP' )

         ENDIF

! FORCE/MOMENT processing. Open L1I which contains FORCE/MOMENT Bulk Data.

            IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                               &
               ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
            CALL LINK_MESSAGE('ALLOCATE MEMORY FOR SYS_LOAD ARRAY          ')
            IF (NFORCE > 0) THEN
               CALL FILE_OPEN ( L1I, LINK1I, OUNT, 'OLD', L1I_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
               CALL LINK_MESSAGE('FORCE/MOMENT PROCESSOR                      ')
               CALL FORCE_MOM_PROC
               CALL FILE_CLOSE ( L1I, LINK1I, L1ISTAT )
            ENDIF
         ENDIF

! Element thermal and pressure load processing

         IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR. (SOL_NAME(1:8) == 'BUCKLING')) THEN
!xx         ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
            CALL LINK_MESSAGE('ALLOCATE MEMORY FOR THERMAL LOAD ARRAYS       ')
            IF ((NPLOAD > 0) .OR. (NTSUB > 0)) THEN
               CALL LINK_MESSAGE('ELEMENT THERMAL AND PRESSURE LOAD PROCESSOR ')
               CALL EPTL
            ENDIF
         ENDIF
         INQUIRE ( FILE=F21FIL, EXIST=LEXIST, OPENED=LOPEN )
         IF (LOPEN) THEN
            CALL FILE_CLOSE ( F21, F21FIL, 'KEEP' )
         ELSE
            CALL FILE_CLOSE ( F21, F21FIL, 'DELETE' )
         ENDIF
! Generate G-set mass matrix, MGG
         IF (SKIPMGG == 'N') THEN

            CALL LINK_MESSAGE('CALCULATE ESTIMATE OF MASS MATRIX SIZE      ')
            CALL EMP0                                        ! Calcs estimate of LTERM_MGGE
            IF (EMP0_PAUSE == 'Y') THEN
               WRITE(SC1,'(A,I12)') 'From EMP0: LTERM_MGGE  = ',LTERM_MGGE
               WRITE(SC1,*) 'Do you want to change LTERM_MGGE estimate? (Y/N)'
               READ(*,*) RESPONSE
               IF ((RESPONSE == 'Y') .OR. (RESPONSE == 'y')) THEN
                  WRITE(SC1,*) 'Enter new LTERM_MGGE'
                  WRITE(SC1,*)
                  READ (*,*) LTERM_MGGE
                  WRITE(SC1,'(A,I12)') 'New LTERM_MGGE will be = ',LTERM_MGGE
               ENDIF
            ENDIF

            CALL LINK_MESSAGE('ALLOCATE MEM FOR EMSKEY, EMSCOL, EMSPNT, EMS')
            CALL ALLOCATE_EMS_ARRAYS ( SUBR_NAME )

            CALL LINK_MESSAGE('ELEMENT MASS MATRIX PROCESSOR               ')
            CALL EMP
            INQUIRE ( FILE=F22FIL, EXIST=LEXIST, OPENED=LOPEN )
            IF (LOPEN) THEN
               CALL FILE_CLOSE ( F22, F22FIL, 'KEEP' )
            ELSE
               CALL FILE_CLOSE ( F22, F22FIL, 'DELETE' )
            ENDIF

! Formulate MGGC mass matrix for concentrated masses

            CALL LINK_MESSAGE('CONCENTRATED MASS MATRIX PROCESSOR          ')
            CALL MGGC_MASS_MATRIX

            CALL LINK_MESSAGE('ALLOCATE MEMORY FOR ELEM MASS ARRAYS        ')
            CALL ALLOCATE_L1_MGG ( 'MGGE', SUBR_NAME )

! Convert system mass matrix from linked list format to sparse format

            CALL LINK_MESSAGE('SPARSE MGG PROCESSOR                        ')
            CALL SPARSE_MGG

            CALL DEALLOCATE_EMS_ARRAYS
            CALL DEALLOCATE_L1_MGG ( 'MGGE'   )
            CALL DEALLOCATE_L1_MGG ( 'MGGC'   )
            CALL DEALLOCATE_L1_MGG ( 'MGGS'   )
            IF ((SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC')) THEN
               CALL DEALLOCATE_MODEL_STUF ( 'CONM2, RCONM2' )
               CALL DEALLOCATE_MODEL_STUF ( 'CMASS, PMASS, RPMASS' )
            ENDIF

         ELSE

!           IF (OLD_NTERM_MGG > 0) THEN
!           ELSE
!           ENDIF

         ENDIF

! Gravity load processing. Open L1P which contains GRAV Bulk Data.

         IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                                  &
            ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
            IF (NGRAV > 0) THEN
               OUNT(1) = ERR
               OUNT(2) = F06
               CALL FILE_OPEN ( L1P, LINK1P, OUNT, 'OLD', L1P_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
               CALL LINK_MESSAGE('GRAV LOAD PROCESSOR                         ')
               CALL GRAV_PROC
               CALL FILE_CLOSE ( L1P, LINK1P, L1PSTAT )
            ENDIF
         ENDIF

! RFORCE load processing. Open L1U which contains RFORCE Bulk Data.

         IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                                  &
            ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
            IF (NRFORCE > 0) THEN
               OUNT(1) = ERR
               OUNT(2) = F06
               CALL FILE_OPEN ( L1U, LINK1U, OUNT, 'OLD', L1U_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
               CALL LINK_MESSAGE('RFORCE LOAD PROCESSOR                       ')
               CALL RFORCE_PROC
               CALL FILE_CLOSE ( L1U, LINK1U, L1USTAT )
            ENDIF
         ENDIF
         CALL DEALLOCATE_L1_MGG ( 'I2_MGG' )

! SLOAD load processing. Open L1W which contains SLOAD Bulk Data.

         IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                                  &
            ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
            IF (NSLOAD > 0) THEN
               OUNT(1) = ERR
               OUNT(2) = F06
               CALL FILE_OPEN ( L1W, LINK1W, OUNT, 'OLD', L1W_MSG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )
               CALL LINK_MESSAGE('SLOAD LOAD PROCESSOR                        ')
               CALL SLOAD_PROC
               CALL FILE_CLOSE ( L1W, LINK1W, L1WSTAT )
            ENDIF
         ENDIF

! Deallocate

         CALL DEALLOCATE_MODEL_STUF ( 'LOAD_SIDS, LOAD_FACS' )

! Convert loads to sparse form and write to file

         IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC') .OR.                                                  &
            ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1))) THEN
            CALL LINK_MESSAGE('CONVERT LOADS TO SPARSE MATRIX FORM         ')
            CALL SPARSE_PG
            CALL DEALLOCATE_MODEL_STUF ( 'SYS_LOAD' )
         ENDIF

! Estimate LTERM so arrays can be allocated for G-set stiffness matrix

         IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN

! KGGD assembly path: delegated to a re-entrant helper so the Phase 4 multi-buckling-subcase driver can rebuild KGGD
! once per buckling subcase (each with its own preload UG_COL) without re-entering all of LINK1. For the single-shot legacy
! path this is behaviorally identical to the previous inline ESP0 / ALLOCATE_STF_ARRAYS / ESP / SPARSE_KGGD sequence.

            CALL BUILD_KGGD_FROM_UG

            ! IF ((SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC') .AND. (SOL_NAME(1:8) /= 'DIFFEREN')) THEN
            !    CALL DEALLOCATE_MODEL_STUF ( 'SCNUM' )
            !    CALL DEALLOCATE_MODEL_STUF ( 'ELDT' )
            !    CALL DEALLOCATE_MODEL_STUF ( 'TPNT, TDATA' )
            !    CALL DEALLOCATE_MODEL_STUF ( 'PPNT, PDATA, PTYPE' )
            !    CALL DEALLOCATE_MODEL_STUF ( 'PLOAD4_3D_DATA' )
            !    CALL DEALLOCATE_MODEL_STUF ( 'GTEMP' )
            ! ENDIF

            INQUIRE ( FILE=F23FIL, EXIST=LEXIST, OPENED=LOPEN )
            IF (LOPEN) THEN
               CALL FILE_CLOSE ( F23, F23FIL, 'KEEP' )
            ELSE
               CALL FILE_CLOSE ( F23, F23FIL, 'DELETE' )
            ENDIF

            CALL DEALLOCATE_IN4_FILES ( 'IN4FIL' )

            INQUIRE ( FILE=F24FIL, EXIST=LEXIST, OPENED=LOPEN )
            IF (LOPEN) THEN
               CALL FILE_CLOSE ( F24, F24FIL, 'KEEP' )
            ELSE
               CALL FILE_CLOSE ( F24, F24FIL, 'DELETE' )
            ENDIF

            ! For SOL 105 multi-buckling decks LINK4 will iterate per buckling subcase, re-invoking BUILD_KGGD_FROM_UG
            ! (via REBUILD_KLLD_FROM_KGGD) which requires MPC_IND_GRIDS to stay alive. LINK4 frees it after its loop.
            IF (.NOT. ((SOL_NAME(1:8) == 'BUCKLING') .AND. (NUM_BUCKLING_SUBS > 1))) THEN
               CALL DEALLOCATE_MODEL_STUF ( 'MPC_IND_GRIDS' )
            ENDIF

         ELSE

            CALL ESP0
            CALL LINK_MESSAGE('CALCULATE ESTIMATE OF KGG MATRIX SIZE         ')
            LTERM_NAME = 'LTERM_KGG'
            LTERM      =  LTERM_KGG

            IF (ESP0_PAUSE == 'Y') THEN
               WRITE(SC1,'(A,A,A,I12)') ' From ESP0: ', LTERM_NAME,' = ',LTERM
               WRITE(SC1,'(A,A,A)') ' Do you want to change ',LTERM_NAME,' estimate? (Y/N)'
               READ(*,*) RESPONSE
               IF ((RESPONSE == 'Y') .OR. (RESPONSE == 'y')) THEN
                  WRITE(SC1,'(A,A)') 'Enter new ', LTERM_NAME
                  WRITE(SC1,*)
                  READ (*,*) LTERM
                  LTERM_KGG  = LTERM
                  WRITE(SC1,'(A,A,A,I12)') 'New ', LTERM_NAME,' will be = ',LTERM
               ENDIF
            ENDIF

            if (setlktk /= 3) then                            ! Subr ESP0 estimated LTERM conservatively. Now allocate this amount
               CALL LINK_MESSAGE('ALLOCATE MEM FOR STFKEY, STFCOL, STFPNT, STF')
               CALL ALLOCATE_STF_ARRAYS ( 'STFKEY', SUBR_NAME )
               CALL ALLOCATE_STF_ARRAYS ( 'STF3', SUBR_NAME )
            else
               Write(err,*) '*ERROR     : PROGRAMMING ERROR IN SUBR ',SUBR_NAME,' SETLKTK CANNOT = 3'
               Write(f06,*) '*ERROR     : PROGRAMMING ERROR IN SUBR ',SUBR_NAME,' SETLKTK CANNOT = 3'
            endif

! Compute element stiffness and merge into system stiffness matrix.

            CALL LINK_MESSAGE('G-SET STIFFNESS MATRIX PROCESSOR            ')
            CALL ESP

            IF ((SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC') .AND. (SOL_NAME(1:8) /= 'DIFFEREN')) THEN
               ! CALL DEALLOCATE_MODEL_STUF ( 'SCNUM' )
               CALL DEALLOCATE_MODEL_STUF ( 'ELDT' )
               CALL DEALLOCATE_MODEL_STUF ( 'TPNT, TDATA' )
               CALL DEALLOCATE_MODEL_STUF ( 'PPNT, PDATA, PTYPE' )
               CALL DEALLOCATE_MODEL_STUF ( 'PLOAD4_3D_DATA' )
               CALL DEALLOCATE_MODEL_STUF ( 'GTEMP' )
            ENDIF

            INQUIRE ( FILE=F23FIL, EXIST=LEXIST, OPENED=LOPEN )
            IF (LOPEN) THEN
               CALL FILE_CLOSE ( F23, F23FIL, 'KEEP' )
            ELSE
               CALL FILE_CLOSE ( F23, F23FIL, 'DELETE' )
            ENDIF

            CALL DEALLOCATE_IN4_FILES ( 'IN4FIL' )

            INQUIRE ( FILE=F24FIL, EXIST=LEXIST, OPENED=LOPEN )
            IF (LOPEN) THEN
               CALL FILE_CLOSE ( F24, F24FIL, 'KEEP' )
            ELSE
               CALL FILE_CLOSE ( F24, F24FIL, 'DELETE' )
            ENDIF

! Convert system stiff matrix from linked list format to sparse format (SPARSE_KGG calls grid singularity check subr)

            CALL LINK_MESSAGE('SPARSE KGG PROCESSOR                        ')
            CALL SPARSE_KGG
            CALL DEALLOCATE_MODEL_STUF ( 'MPC_IND_GRIDS' )
            CALL DEALLOCATE_STF_ARRAYS ( 'STFKEY' )
            CALL DEALLOCATE_STF_ARRAYS ( 'STF3' )

         ENDIF

! Write DOF tables and deallocate

         WRITE(SC1,*) CR13
         CALL LINK_MESSAGE('WRITE DOF TABLES TO FILE AND DEALLOCATE     ')
         CALL WRITE_DOF_TABLES
         CALL FILE_CLOSE ( L1C, LINK1C, 'KEEP' )

! Write element data to L1G. Save L1G for use in LINK9.
         IF (LOAD_ISTEP == 1) THEN

            CALL FILE_OPEN ( L1G, LINK1G, OUNT, 'REPLACE', L1G_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

            CALL LINK_MESSAGE('WRITE ELEMENT DATA TO FILE                  ')

            CALL ELSAVE

            IF ((SOL_NAME(1:8) /= 'BUCKLING') .AND. (SOL_NAME(1:8) /= 'NLSTATIC') .AND. (SOL_NAME(1:8) /= 'DIFFEREN')) THEN
               CALL DEALLOCATE_MODEL_STUF ( 'ETYPE, EDAT, EPNT' )
               CALL DEALLOCATE_MODEL_STUF ( 'VVEC, OFFSETS, PLATE stuff' )
               CALL DEALLOCATE_MODEL_STUF ( 'ELEM PROPERTIES AND MATERIALS' )
               CALL DEALLOCATE_MODEL_STUF ( 'EOFF' )
               CALL DEALLOCATE_MODEL_STUF ( 'ESORT1' )
               CALL DEALLOCATE_MODEL_STUF ( 'ESORT2' )
            ENDIF

            CALL FILE_CLOSE ( L1G, LINK1G, 'KEEP' )

         ENDIF

      ENDIF res19

! Deallocate
!
! For multi-buckling-subcase SOL 105 we defer these dealloc's so the LK4 multi-statsub driver can re-enter
! the LK1/LK2 reduction chain (via REBUILD_KLLD_FROM_KGGD -> BUILD_KGGD_FROM_UG -> ESP) without crashing on
! freed AGRID/BGRID/SUBLOD/etc. LK4 (or LK9) cleans them up after the buckling loop is done.

      IF (.NOT. ((SOL_NAME(1:8) == 'BUCKLING') .AND. (NUM_BUCKLING_SUBS > 1) .AND. (LOAD_ISTEP == 2))) THEN
         CALL DEALLOCATE_MODEL_STUF ( 'SINGLE ELEMENT ARRAYS' )
         IF ((SOL_NAME(1:8) /= 'NLSTATIC') .AND. (SOL_NAME(1:8) /= 'DIFFEREN')) THEN
            CALL DEALLOCATE_MODEL_STUF ( 'SUBLOD' )
         ENDIF
      ENDIF

! Check allocation status of allocatable arrays, if requested

      IF (DEBUG(100) > 0) THEN
         CALL CHK_ARRAY_ALLOC_STAT
         IF (DEBUG(100) > 1) THEN
            CALL WRITE_ALLOC_MEM_TABLE ( 'at the end of '//SUBR_NAME )
         ENDIF
      ENDIF

! Write data to L1A (this rewrites L1A if this a restart - needed to tell LINK9 that this is a restart)

      COMM(LINKNO) = 'C'
res20:IF (RESTART == 'N') THEN
         CALL WRITE_L1A ( 'KEEP', 'Y' )
      ENDIF res20

! Write LINK1 end to F06

      CALL OURTIM
      WRITE(F06,151) LINKNO

      IF (( DEBUG(193) == 1) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK1' )
      ENDIF

! Write LINK1 end to screen
      WRITE(SC1,154) LINKNO

! **********************************************************************************************************************************
  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END',/)

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  154 FORMAT(  ' >> LINK',I3,' END')

 9998 FORMAT(/,' PROCESSING TERMINATED DUE TO ',I8,' INPUT ERRORS. CHECK OUTPUT FILE FOR ERROR MESSAGES')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE INIT_COUNTERS

      USE PENTIUM_II_KIND
      USE SCONTR

      IMPLICIT NONE

! **********************************************************************************************************************************
      NTERM_KGG  = 0   ;   LTERM_KGG  = 0
      NTERM_KGGD = 0   ;   LTERM_KGGD = 0
      NTERM_MGG  = 0
      NTERM_MGGC = 0
      NTERM_MGGE = 0   ;   LTERM_MGGE = 0
      NTERM_MGGS = 0
      NTERM_PG   = 0
!     NTERM_RMG  = 0
!     NTERM_RMM  = 0
!     NTERM_RMN  = 0

! **********************************************************************************************************************************

      END SUBROUTINE INIT_COUNTERS

      END SUBROUTINE LINK1

   END MODULE LINK0_LINK1_MOD
