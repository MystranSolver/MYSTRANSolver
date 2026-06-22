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

      SUBROUTINE LINK9 ( LK9_PROC_NUM )

! Main driver for calculating outputs requested in Case Control once the G-set unknowns have been solved for in prior LINK's

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR

      USE IOUNT1, ONLY                :  ERR, F06, F25, L1E, L1M, L1R, L2A, L2B, L2C, L2D, L2I, L2J, L2R, L2S,                &
                                         L5A, L5B, NEU, OT4, OU4, PCH, SC1

      USE IOUNT1, ONLY                :  F06FIL, F25FIL, LINK1B, LINK1E, LINK1M, LINK1R, LINK2A, LINK2B, LINK2C, LINK2D,           &
                                         LINK2I, LINK2J, LINK2R, LINK2S, LINK5A, LINK5B, MOT4  , MOU4  , NEUFIL, OT4FIL, OU4FIL,   &
                                         PCHFIL

      USE IOUNT1, ONLY                :  L1ASTAT, L1ESTAT, L1MSTAT, L1RSTAT, L2ASTAT, L2BSTAT, L2CSTAT, L2ISTAT, L2JSTAT, L2RSTAT, &
                                         L2SSTAT, OT4STAT, OU4STAT, PCHSTAT

      USE IOUNT1, ONLY                :  F25_MSG, L1E_MSG, L1M_MSG, L1R_MSG, L2A_MSG, L2B_MSG, L2C_MSG, L2D_MSG, L2I_MSG,          &
                                         L2J_MSG, L2R_MSG, L2S_MSG, L5A_MSG, L5B_MSG, NEU_MSG, PCH_MSG,                            &
                                         OT4_MSG, OU4_MSG, OT4_GRD_OTM, OT4_ELM_OTM, OU4_GRD_OTM, OU4_ELM_OTM

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, CC_ENTRY_LEN, COMM, IBIT, INT_SC_NUM, JTSUB, FATAL_ERR,                     &
                                         FEMAP_VERSION, LINKNO, MBUG,                                                              &
                                         NDOFF, NDOFG, NDOFL, NDOFM, NDOFN, ndofo, NDOFR, NDOFS, NDOFSA, NGRID, NSUB, NVEC,        &
                                         NTERM_IF_LTM, NTERM_GMN, NTERM_HMN, NTERM_KFS, NTERM_KFSD, NTERM_LMN, NTERM_MFS,          &
                                         NTERM_MGG, NTERM_MLL,NTERM_PG, NTERM_PM, NTERM_PS, NTERM_QSYS,                            &
                                         NUM_CB_DOFS, NUM_EIGENS,                                                                  &
                                         NROWS_OTM_ACCE, NROWS_OTM_DISP, NROWS_OTM_MPCF, NROWS_OTM_SPCF,                           &
                                         NROWS_OTM_ELFE, NROWS_OTM_ELFN, NROWS_OTM_STRE, NROWS_OTM_STRN,                           &
                                         NROWS_TXT_ACCE, NROWS_TXT_DISP, NROWS_TXT_MPCF, NROWS_TXT_SPCF,                           &
                                         NROWS_TXT_ELFE, NROWS_TXT_ELFN, NROWS_TXT_STRE, NROWS_TXT_STRN, RESTART, SOL_NAME, WARN_ERR	
! --- cbeam_add begin --- !
      USE SCONTR, ONLY                :  INT_EIG_NUM, NUM_BUCKLING_SUBS,MODE_SUBCASE
! --- cbeam_add end --- !

      USE SCONTR, ONLY                :  GROUT_ACCE_BIT, GROUT_DISP_BIT, GROUT_OLOA_BIT, GROUT_SPCF_BIT, GROUT_MPCF_BIT,           &
                                         GROUT_GPFO_BIT, ELOUT_ELFN_BIT, ELOUT_ELFE_BIT, ELOUT_STRE_BIT, ELOUT_STRN_BIT,           &
                                         ELDT_F25_U_P_BIT

      USE CC_OUTPUT_DESCRIBERS, ONLY  :  DISP_OUT, ACCE_OUT, OLOA_OUT, SPCF_OUT, MPCF_OUT, FORC_OUT, GPFO_OUT, STRE_OUT, STRN_OUT
      USE TIMDAT, ONLY                :  STIME
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  EPSIL, MPFOUT, SUPINFO, SUPWARN, WTMASS, PRTF06, PRTOP2, PRTNEU
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

      USE LINK9_USE_IFs
      USE LINK_MESSAGE_Interface

      IMPLICIT NONE

      LOGICAL                         :: WRITE_F06, WRITE_OP2, WRITE_PCH, WRITE_NEU   ! flag
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
      ! setup PRTF06, PRTNEU, PRTOP2
      !IF (DEBUG(200) > 0) THEN
      !   PRTNEU = 'Y'
      !ENDIF
      !IF (PRTNEU == 'Y') THEN
      !   DEBUG(200) = 1
      !ENDIF
      IF (PRTF06 == 'Y') THEN
         DISP_OUT(1:1) = 'Y'  ! f06
         ACCE_OUT(2:2) = 'Y'
         OLOA_OUT(1:1) = 'Y'
         SPCF_OUT(1:1) = 'Y'
         MPCF_OUT(1:1) = 'Y'
         STRE_OUT(1:1) = 'Y'
         STRN_OUT(1:1) = 'Y'
         FORC_OUT(1:1) = 'Y'
         GPFO_OUT(1:1) = 'Y'
      ENDIF
      IF (PRTOP2 == 'Y') THEN
         DISP_OUT(2:2) = 'Y'  ! op2
         ACCE_OUT(2:2) = 'Y'
         OLOA_OUT(2:2) = 'Y'
         SPCF_OUT(2:2) = 'Y'
         MPCF_OUT(2:2) = 'Y'
         STRE_OUT(2:2) = 'Y'
         STRN_OUT(2:2) = 'Y'
         FORC_OUT(2:2) = 'Y'
         GPFO_OUT(2:2) = 'Y'
      ENDIF
      !DISP_OUT(3:3) = 'Y'  ! pch
      !------------------------------------------------------------------------

      ! Make units for writing errors the screen until we open output files
      OUNT(1) = SC1
      OUNT(2) = SC1

      ! Make units for writing errors the error file and output file
      OUNT(1) = ERR
      OUNT(2) = F06

      WRITE_F06 = (DISP_OUT(1:1) == 'Y')
      WRITE_OP2 = (DISP_OUT(2:2) == 'Y')
      WRITE_PCH = (DISP_OUT(3:3) == 'Y')
      IF (WRITE_PCH) THEN
         INQUIRE (FILE=PCHFIL, OPENED=LOPEN)
         IF (.NOT.LOPEN) THEN                          ! Otherwise we assume it is positioned at its end and ready for write
            CALL FILE_OPEN ( PCH, PCHFIL, OUNT, 'OLD', PCH_MSG, 'WRITE_STIME', 'FORMATTED', 'READWRITE', 'REWIND', 'Y', 'Y' )
         ENDIF
      ENDIF

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
! --- cbeam_add begin --- !
         !IF (LK9_PROC_NUM == 1) THEN	! original cbeam_add 
         !   NUM_SOLNS = 1              ! original cbeam_add
! --- cbeam_add end  --- !
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
! --- cbeam_add begin --- !
!            INT_SC_NUM   = LK9_PROC_NUM
!            FEMAP_SET_ID = LK9_PROC_NUM
! --- cbeam_add end --- !
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
            INT_SC_NUM   = 1
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
                  CALL GET_FG_INERTIA_FORCES
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
