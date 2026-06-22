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

      SUBROUTINE GP_FORCE_BALANCE_PROC ( JVEC, IHEADER )

      ! Processes grid point force balance output requests for one subcase.
      ! The effects of all forces on grids are included so totals
      ! should be zero

      USE PENTIUM_II_KIND, ONLY       :  BYTE, SHORT, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, OP2, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, GROUT_GPFO_BIT, IBIT, INT_SC_NUM, JTSUB, NDOFG, NDOFM, MELDOF, NDOFO, NDOFR,&
                                         NELE, NGRID, NUM_CB_DOFS, NVEC, SOL_NAME, MODE_SUBCASE
      USE CONSTANTS_1, ONLY           :  ZERO, ONE_HUNDRED
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  AGRID, EID, ELGP, ESORT1, ETYPE, GRID, GRID_ELEM_CONN_ARRAY, GRID_ID, &
                                         GROUT, LABEL, PLY_NUM, PEG, PTE, SCNUM, STITLE, SUBLOD, TITLE, TYPE
      USE COL_VECS, ONLY              :  FG_COL, PG_COL, QGm_COL, QGs_COL, QGr_COL, UG_COL
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  GPFO_OUT
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP

! --- cbeam_add begin --- !
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  NUM_EMG_FATAL_ERRS
      USE PARAMS, ONLY                :  EPSIL
      USE LINK9_STUFF, ONLY           :  GID_OUT_ARRAY
! --- cbeam_add end --- !
      USE GP_FORCE_BALANCE_PROC_USE_IFs

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GP_FORCE_BALANCE_PROC'
      CHARACTER(13*BYTE)              :: CHAR_PCT(6)       ! Character representation of MEFFMASS sum percents of total model mass
      CHARACTER(LEN=*) , INTENT(IN)   :: IHEADER           ! Indicator of whether to write an output header
      CHARACTER( 1*BYTE), PARAMETER   :: COORD_SYS   = 'G' ! Subr TRANSFORM_NODE_FORCES will calc elem node forces in global coords
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option flags for subr ELMTLB (to tell it what to transform)

      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG)                   :: AELEM             ! Actual element ID
      INTEGER(LONG)                   :: BNDY_COMP          ! Component (1-6) for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BNDY_GRID          ! Grid for a boundary DOF in CB analyses
      INTEGER(LONG)                   :: BNDY_DOF_NUM       ! DOF number for BNDY_GRID/BNDY_COMP
      INTEGER(LONG)                   :: BEG_ROW_NUM       ! Row num in PEG where we begin to get the 6 values of elem load KE*UE
      INTEGER(LONG)                   :: G_CID             ! Global coord system for GRID_NUM
      INTEGER(LONG)                   :: GRID_NUM          ! Actual grid point number
      INTEGER(LONG)                   :: G_SET_COL         ! Col number in TDOF where the G-set DOF's exist
      INTEGER(LONG)                   :: GDOF              ! G-set DOF number
      INTEGER(LONG)                   :: I,J,K,L           ! DO loop indices
      INTEGER(LONG)                   :: IB                ! Intermediate value used in determining NREQ
      INTEGER(LONG)                   :: IELE              ! Internal element ID
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: NREQ              ! Number of grids for which G.P. force balance is requested
      INTEGER(LONG)                   :: NUM_COMPS         ! Either 6 or 1 depending on whether grid is a physical grid or a SPOINT
      INTEGER(LONG)                   :: NUM_CONN_ELEMS    ! The number of elements that are connected to a specific grid
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: TDOF_ROW          ! Row no. in array TDOF to find GDOF DOF number


      INTEGER(SHORT), DIMENSION(1)    :: Udd = (/0220/)    ! 0220 is the 4 digit ASCII code for a capital U double-dot

      REAL(DOUBLE)                    :: DUM_KE(MELDOF,MELDOF)
      REAL(DOUBLE)                    :: FG1(6)            ! The 6 vals from FG_COL (grid inertia forces) for 1 grid point
      REAL(DOUBLE)                    :: PG1(6)            ! The 6 vals from PG_COL (grid applied loads) for 1 grid point
      REAL(DOUBLE)                    :: MAX_ABS(6)        ! The 6 abs largest values of any force item in the 6 components
      REAL(DOUBLE)                    :: MAX_ABS_PCT(6)    ! Modal mass as % of total mass
      REAL(DOUBLE)                    :: PEG1(6)           ! 6 vals of elem load for one grid point of an elem
      REAL(DOUBLE)                    :: PTE1(6)           ! 6 vals from thermal load, PTE, for one grid point, one elem
      REAL(DOUBLE)                    :: PTET(6)           ! 6 vals from thermal load, PTE, for one grid point, all elems
      REAL(DOUBLE)                    :: QGm1(6)           ! 6 vals from QGm_COL (G-set MPC constraint forces) for one grid point
      REAL(DOUBLE)                    :: QGr1(6)           ! 6 vals from QGr_COL (R-set I/F forces expanded to G-set) for 1 grid pt
      REAL(DOUBLE)                    :: QGs1(6)           ! 6 vals from QGs_COL (S-set SPC forces expanded to G-set) for 1 grid pt
      REAL(DOUBLE)                    :: TOTALS(6)         ! 6 totals of FG1, PG1, QGs1, QGm1, PTET, PEG1 (should be zero's)

      INTRINSIC IAND

      integer(long)                   :: max_abs_grid(6)    ! Grid where max-abs for GP force balance totals exists
      real(double)                    :: max_abs_all_grds(6)! The 6 max-abs from GP force balance totals

      INTEGER(LONG)                   :: NROWS, NCOLS, NNODE_GPFORCE, INODE_GPFORCE, IERR  ! GPFORCE table helper

      LOGICAL                         :: WRITE_F06, WRITE_OP2, IS_GPFORCE_SUMMARY_INFO  ! flags
      LOGICAL                         :: IS_MODES, IS_THERMAL, IS_APP, IS_SPC, IS_MPC

      INTEGER(LONG)                   :: ISUBCASE_INDEX         ! helper to get the title/subcase
      INTEGER(LONG)                   :: ISUBCASE               ! the correct subcase id

      ! op2 stuff
      ! table -3 info
      INTEGER(LONG)                   :: DEVICE_CODE        !
      INTEGER(LONG)                   :: ANALYSIS_CODE          ! static/modal/time/etc. flag
      INTEGER(LONG)                   :: ELEMENT_TYPE           ! the OP2 flag for the element
      LOGICAL                         :: FIELD_5_INT_FLAG       ! flag to trigger FIELD5_INT_MODE vs.
      INTEGER(LONG)                   :: FIELD5_INT_MODE        ! int value for field 5
      REAL(DOUBLE)                    :: FIELD6_EIGENVALUE      ! float value for field 6
      CHARACTER(LEN=128)              :: TITLEI                 ! the model TITLE
      CHARACTER(LEN=128)              :: SUBTITLEI              ! the subcase SUBTITLE
      CHARACTER(LEN=128)              :: LABELI                 ! the subcase LABEL
      INTEGER(LONG)                   :: ITABLE                 ! the table counter

      INTEGER, ALLOCATABLE            :: GPFORCE_NID_EID(:,:)    ! currently unused
      CHARACTER*8, ALLOCATABLE        :: GPFORCE_ETYPE(:)        ! currently unused
      REAL, ALLOCATABLE               :: GPFORCE_FXYZ_MXYZ(:,:)  ! currently unused



! **********************************************************************************************************************************

      ! GPFORCE is unsupported for buckling decks
      IF (SOL_NAME(1:8) == "BUCKLING") THEN
         RETURN
      ENDIF

      ! Print some summary info for max abs value of GP force balance for each solution vector
      IS_GPFORCE_SUMMARY_INFO = (DEBUG(192) > 0)

      IS_THERMAL = (SUBLOD(INT_SC_NUM,2) > 0)
      IS_MODES = ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL'))

      WRITE_F06 = (GPFO_OUT(1:1) == 'Y')
      INQUIRE ( UNIT=OP2, OPENED=WRITE_OP2 )      !WRITE_OP2 = (GPFO_OUT(2:2) == 'Y')

      ! OP2: Write output headers if this is not the first use of this subr.
      ANALYSIS_CODE = -1
      FIELD_5_INT_FLAG = .TRUE.
      FIELD5_INT_MODE = 0
      FIELD6_EIGENVALUE = 0.0

      INODE_GPFORCE = 1

      ! Initialize
      DO I=1,6
         FG1(I)         = ZERO
         PG1(I)         = ZERO
         MAX_ABS(I)     = ZERO
         MAX_ABS_PCT(I) = ZERO
         PEG1(I)        = ZERO
         PTE1(I)        = ZERO
         PTET(I)        = ZERO
         QGM1(I)        = ZERO
         QGR1(I)        = ZERO
         QGS1(I)        = ZERO
         TOTALS(I)      = ZERO
      ENDDO

      CALL CHK_COL_ARRAYS_ZEROS

      ! Write output headers.
      ANALYSIS_CODE = -1
      IF (IHEADER == 'Y') THEN   ! should we write the header?
         IF (WRITE_F06) THEN
            WRITE(F06,*)
            WRITE(F06,*)
         ENDIF
         !CALL WRITE_SUBCASE_EIGENVEC_HEADER(JVEC, WRITE_F06)
         ISUBCASE_INDEX = 0
         IF    (SOL_NAME(1:7) == 'STATICS') THEN
            ISUBCASE_INDEX = JVEC
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = SCNUM(JVEC)
            IF (WRITE_F06)  WRITE(F06,9101) SCNUM(JVEC)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN

            ISUBCASE_INDEX = JVEC
            ANALYSIS_CODE = 1
            FIELD5_INT_MODE = SCNUM(JVEC)

         ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN

            ISUBCASE_INDEX = 1  ! modes
            ANALYSIS_CODE = 2
            FIELD5_INT_MODE = JVEC

         ELSE IF (SOL_NAME(1:5) == 'MODES') THEN
            ! Per-mode subcase attribution (LINK4 populates MODE_SUBCASE). Falls back to 1 for legacy single-METHOD paths.
            ISUBCASE_INDEX = 1  ! modes
            ANALYSIS_CODE = 2
            FIELD5_INT_MODE = JVEC
            IF (WRITE_F06)  WRITE(F06,9102) JVEC

            IF (ALLOCATED(MODE_SUBCASE)) THEN
               IF (JVEC <= SIZE(MODE_SUBCASE)) ISUBCASE_INDEX = MODE_SUBCASE(JVEC)
            ENDIF



         ELSE IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN  ! Write info on what CB DOF the output is for
            ISUBCASE_INDEX = 1  ! modes

            IF ((JVEC <= NDOFR) .OR. (JVEC >= NDOFR+NVEC)) THEN
               IF (JVEC <= NDOFR) THEN
                  BNDY_DOF_NUM = JVEC
               ELSE
                  BNDY_DOF_NUM = JVEC-(NDOFR+NVEC)
               ENDIF
               CALL GET_GRID_AND_COMP ( 'R ', BNDY_DOF_NUM, BNDY_GRID, BNDY_COMP  )
            ENDIF

            IF (WRITE_F06) THEN
              IF       (JVEC <= NDOFR) THEN
                 WRITE(F06,9103) JVEC, NUM_CB_DOFS, 'acceleration', BNDY_GRID, BNDY_COMP
              ELSE IF ((JVEC > NDOFR) .AND. (JVEC <= NDOFR+NVEC)) THEN
                 WRITE(F06,9105) JVEC, NUM_CB_DOFS, JVEC-NDOFR
              ELSE
                 WRITE(F06,9103) JVEC, NUM_CB_DOFS, 'displacement', BNDY_GRID, BNDY_COMP
              ENDIF
            ENDIF
         ENDIF

         ISUBCASE = SCNUM(ISUBCASE_INDEX)

         ! -- F06 header for TITLE, SUBTITLE, LABEL (but only to F06)
         TITLEI = TITLE(INT_SC_NUM)
         SUBTITLEI = STITLE(INT_SC_NUM)
         LABELI = LABEL(INT_SC_NUM)

         IF (WRITE_F06) THEN
            IF (TITLEI(1:)  /= ' ') THEN
               WRITE(F06,9799) TITLE
            ENDIF

            IF (SUBTITLEI(1:) /= ' ') THEN
               WRITE(F06,9799) SUBTITLEI
            ENDIF

            IF (LABELI(1:)  /= ' ') THEN
               WRITE(F06,9799) LABELI
            ENDIF

            WRITE(F06,*)

            ! write f06 header
            IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
               WRITE(F06,8999)
            ELSE
               WRITE(F06,9200)
            ENDIF
         ENDIF

      ENDIF

      ! Process grid point force balance output requests.
      DO I=1,6
         FG1(I)     = ZERO
         PG1(I)     = ZERO
         PTET(I)    = ZERO
         QGs1(I)    = ZERO
         QGm1(I)    = ZERO
         PEG1(I)    = ZERO
         MAX_ABS(I) = ZERO
      ENDDO

      IF (IS_GPFORCE_SUMMARY_INFO) THEN                    ! If DEBUG 192 = 1, GPFO summary will be printed out to F06
         DO L=1,6                                          ! (1) Initialize the max valueof the 6 TOTAL components
            MAX_ABS_ALL_GRDS(L) = ZERO
         ENDDO
i_do1:   DO I=1,NGRID                                      ! (2) Set initial value for the grid at which the max occurs so that if
                                                           !     all TOTAL values are 0 we set the grid for max at the initial grid
            IB = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_GPFO_BIT))
            IF (IB > 0) THEN
               DO L=1,6
                  MAX_ABS_GRID(L) = GRID(I,1)
               ENDDO
               EXIT i_do1
            ENDIF
         ENDDO i_do1
      ENDIF

      NREQ = 0
      DO I=1,NGRID
         IB = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_GPFO_BIT))
         NREQ = NREQ + 1
      ENDDO

      DEVICE_CODE = 1

      NNODE_GPFORCE = NGRID
      DO I=1,NGRID
         IB = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_GPFO_BIT))
         GRID_NUM  = GRID(I,1)
         CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )

         IF ((IB > 0) .AND. (NUM_COMPS == 6)) THEN         ! Do not do force balance for SPOINT's
            G_CID = GRID(I,3)

            ! applied load, thermal, spc force, mpc force
            NNODE_GPFORCE = NNODE_GPFORCE + 4

            NUM_CONN_ELEMS = GRID_ELEM_CONN_ARRAY(I,2)
            NNODE_GPFORCE = NNODE_GPFORCE + NUM_CONN_ELEMS	

            !IF (IS_THERMAL) THEN
            !   NUM_CONN_ELEMS = GRID_ELEM_CONN_ARRAY(I,2)
            !
            !   DO J=1,NUM_CONN_ELEMS
            !      AELEM = GRID_ELEM_CONN_ARRAY(I,2+J)
            !      CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, AELEM, IELE )
            !      CALL GET_ELGP ( IELE )
            !      CALL GET_ELEM_AGRID_BGRID ( IELE, 'N' )
            !      DO K=1,ELGP
            !         IF (AGRID(K) == GRID_NUM) THEN
            !             NNODE_GPFORCE = NNODE_GPFORCE + 1
            !         ENDIF
            !      ENDDO
            !   ENDDO
            !ENDIF
         ENDIF
      ENDDO

      !------------------------------------------------------------------------
      ! ALLOCATE: GPFORCE_NID_EID, GPFORCE_ETYPE, GPFORCE_FXYZ_MXYZ
      !ref mystran SUB ALLOCATE_DOF_TABLES
      !
      !KTSTACK(5500,3)
      !
      IF(WRITE_OP2) THEN
          NROWS = NNODE_GPFORCE
          NCOLS = 2
          ALLOCATE (GPFORCE_NID_EID(NROWS,2),STAT=IERR)

          IF (IERR == 0) THEN
              DO I=1,NROWS
                  DO J=1,NCOLS
                      GPFORCE_NID_EID(I,J) = 0
                  ENDDO
              ENDDO
          ELSE
              WRITE(ERR,*) 'GPFORCE_NID_EID ALLOCATED error'
          ENDIF

          !----------------
          IF (ALLOCATED(GPFORCE_ETYPE)) THEN
              WRITE(ERR,*) 'ALLOCATED!'
          ELSE
              ALLOCATE (GPFORCE_ETYPE(NROWS),STAT=IERR)

              IF (IERR == 0) THEN
                  DO I=1,NROWS
                      GPFORCE_ETYPE(I) = "NA"
                  ENDDO
              ELSE
                  WRITE(ERR,*) 'GPFORCE_ETYPE ALLOCATED error'
              ENDIF
          ENDIF
          !----------------
          NCOLS = 6
          ALLOCATE (GPFORCE_FXYZ_MXYZ(NROWS,NCOLS),STAT=IERR)

          IF (IERR == 0) THEN
              DO I=1,NROWS
                  DO J=1,NCOLS
                      GPFORCE_FXYZ_MXYZ(I,J) = 0
                  ENDDO
              ENDDO
          ELSE
              WRITE(6,*) 'GPFORCE_FXYZ_MXYZ ALLOCATED error'
          ENDIF
      ELSE
         ! Dummy allocations to suppress warnings
         ALLOCATE (GPFORCE_NID_EID(0,0))
         ALLOCATE (GPFORCE_ETYPE(0))
         ALLOCATE (GPFORCE_FXYZ_MXYZ(0,0))

      ENDIF  ! write_op2 allocation
      !------------------------------------------------------------------------
      CALL COUNTER_INIT('Process grid ', NGRID)
      DO I=1,NGRID
         IB = IAND(GROUT(I,INT_SC_NUM),IBIT(GROUT_GPFO_BIT))
         GRID_NUM  = GRID(I,1)
         CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )

         IF ((IB > 0) .AND. (NUM_COMPS == 6)) THEN         ! Do not do force balance for SPOINT's
            G_CID = GRID(I,3)
            IF (WRITE_F06) THEN
               WRITE(F06,9201) GRID_NUM, G_CID
               WRITE(F06,9202)
            ENDIF

            ! Get element equiv thermal loads so we can sub them from applied loads
            OPT(1) = 'N'  ! OPT(1) is for calc of ME
            OPT(2) = 'Y'  ! OPT(2) is for calc of PTE
            OPT(3) = 'N'  ! OPT(3) is for calc of SEi, STEi
            OPT(4) = 'N'  ! OPT(4) is for calc of KE-linear
            OPT(5) = 'N'  ! OPT(5) is for calc of PPE
            OPT(6) = 'N'  ! OPT(6) is for calc of KE-diff stiff
            DO J=1,6
               PTET(J) = ZERO
            ENDDO
            IF (IS_THERMAL) THEN
               NUM_CONN_ELEMS = GRID_ELEM_CONN_ARRAY(I,2)
               DO J=1,NUM_CONN_ELEMS
                  AELEM = GRID_ELEM_CONN_ARRAY(I,2+J)
                  CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, AELEM, IELE )
                  CALL GET_ELGP ( IELE )
                  CALL GET_ELEM_AGRID_BGRID ( IELE, 'N' )
                  DO K=1,ELGP
                     IF (AGRID(K) == GRID_NUM) THEN
                        PLY_NUM = 0                        ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( IELE, OPT, 'N', SUBR_NAME, 'N' )
                        IF (ETYPE(IELE)(1:4) /= 'ELAS') THEN
                           CALL ELEM_TRANSFORM_LBG ( 'PTE', dum_ke, PTE )
                        ENDIF
                        DO L=1,6
                           BEG_ROW_NUM = 6*(K-1) + 1
                           PTE1(L) = PTE(BEG_ROW_NUM + L - 1,JTSUB)
                        ENDDO
                        DO L=1,6
                           PTET(L) = PTET(L) + PTE1(L)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

            ! get applied load, thermal load, SPC force, MPC force for a single node
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_NUM, IGRID )
            ROW_NUM_START = TDOF_ROW_START(IGRID)
            CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
            DO J=1,NUM_COMPS
               CALL TDOF_COL_NUM ( 'G ', G_SET_COL )
               TDOF_ROW = ROW_NUM_START + J - 1
               GDOF = TDOF(TDOF_ROW,G_SET_COL)
               IF (IS_MODES) THEN
                  FG1(J) = FG_COL(GDOF)
               ENDIF
               PG1(J)  = PG_COL(GDOF)
               QGs1(J) = QGs_COL(GDOF)
               IF (NDOFM > 0) THEN
                  QGm1(J) = QGm_COL(GDOF)
               ENDIF

               ! PTET = 0 unless INT_SC_NUM is a thermal subcase (see above)
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  QGr1(J) = QGr_COL(GDOF)
                  TOTALS(J) = -FG1(J) + PG1(J) - PTET(J) + QGs1(J) + QGm1(J) + QGr1(J)
               ELSE
                  TOTALS(J) = -FG1(J) + PG1(J) - PTET(J) + QGs1(J) + QGm1(J)
               ENDIF
            ENDDO

            IS_APP = IS_ABS_POSITIVE(PG1)
            IS_SPC = IS_ABS_POSITIVE(QGs1)
            IS_MPC = IS_ABS_POSITIVE(QGm1)

            IF(WRITE_OP2) THEN
                IF(IS_APP) THEN
                    GPFORCE_NID_EID(INODE_GPFORCE,1) = GRID_NUM
                    GPFORCE_NID_EID(INODE_GPFORCE,2) = 0
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,1) = PG1(1)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,2) = PG1(2)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,3) = PG1(3)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,4) = PG1(4)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,5) = PG1(5)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,6) = PG1(6)
                    GPFORCE_ETYPE(INODE_GPFORCE) = 'APP-LOAD'
                    INODE_GPFORCE = INODE_GPFORCE + 1
                ENDIF
                IF(IS_THERMAL) THEN
                    GPFORCE_NID_EID(INODE_GPFORCE,1) = GRID_NUM
                    GPFORCE_NID_EID(INODE_GPFORCE,2) = EID
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,1) = -PTET(1)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,2) = -PTET(2)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,3) = -PTET(3)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,4) = -PTET(4)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,5) = -PTET(5)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,6) = -PTET(6)
                    GPFORCE_ETYPE(INODE_GPFORCE) = 'THERMAL'
                    INODE_GPFORCE = INODE_GPFORCE + 1
                ENDIF
                IF(IS_SPC) THEN
                    GPFORCE_NID_EID(INODE_GPFORCE,1) = GRID_NUM
                    GPFORCE_NID_EID(INODE_GPFORCE,2) = EID
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,1) = QGs1(1)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,2) = QGs1(2)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,3) = QGs1(3)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,4) = QGs1(4)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,5) = QGs1(5)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,6) = QGs1(6)
                    GPFORCE_ETYPE(INODE_GPFORCE) = 'F-OF-SPC'
                    INODE_GPFORCE = INODE_GPFORCE + 1
                ENDIF
                IF(IS_MPC) THEN
                    GPFORCE_NID_EID(INODE_GPFORCE,1) = GRID_NUM
                    GPFORCE_NID_EID(INODE_GPFORCE,2) = EID
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,1) = QGm1(1)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,2) = QGm1(2)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,3) = QGm1(3)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,4) = QGm1(4)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,5) = QGm1(5)
                    GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,6) = QGm1(6)
                    GPFORCE_ETYPE(INODE_GPFORCE) = 'F-OF-MPC'
                    INODE_GPFORCE = INODE_GPFORCE + 1
                ENDIF
            ELSE
                IF(IS_APP) INODE_GPFORCE = INODE_GPFORCE + 1
                IF(IS_THERMAL) INODE_GPFORCE = INODE_GPFORCE + 1
                IF(IS_SPC) INODE_GPFORCE = INODE_GPFORCE + 1
                IF(IS_MPC) INODE_GPFORCE = INODE_GPFORCE + 1
            ENDIF

            IF (WRITE_F06) THEN
               IF (IS_APP) WRITE(F06,9203) (PG1(J),J=1,6)  ! applied load
               IF (IS_THERMAL) THEN
                  WRITE(F06,9204) (-PTET(J),J=1,6)  ! thermal
               ENDIF
               IF(IS_SPC) WRITE(F06,9205) (QGs1(J),J=1,6)  ! spc force
               IF(IS_MPC) WRITE(F06,9206) (QGm1(J),J=1,6)  ! mpc force

               IF (IS_MODES) THEN
                  DO J=1,6
                     IF(FG1(J) == ZERO) THEN;  FG1(J) = -ZERO;  ENDIF  ! Avoids writing -0.0 for -FG1(J) below
                  ENDDO
                  WRITE(F06,9207) (-FG1(J),J=1,6)  ! inertia force
               ENDIF
               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  WRITE(F06,9208) (QGr1(J),J=1,6)
               ENDIF
            ENDIF

            ! Calc elem forces
            OPT(1) = 'N'  ! OPT(1) is for calc of ME
            OPT(2) = 'Y'  ! OPT(2) is for calc of PTE
            OPT(3) = 'N'  ! OPT(3) is for calc of SEi, STEi
            OPT(4) = 'Y'  ! OPT(4) is for calc of KE
            OPT(5) = 'N'  ! OPT(5) is for calc of PPE
            OPT(6) = 'N'  ! OPT(6) is for calc of KE-diff stiff
            NUM_CONN_ELEMS = GRID_ELEM_CONN_ARRAY(I,2)
            DO J=1,NUM_CONN_ELEMS
               AELEM = GRID_ELEM_CONN_ARRAY(I,2+J)
               CALL GET_ARRAY_ROW_NUM ( 'ESORT1', SUBR_NAME, NELE, ESORT1, AELEM, IELE )
               CALL GET_ELGP ( IELE )
               CALL GET_ELEM_AGRID_BGRID ( IELE, 'N' )
               DO K=1,ELGP
                  IF (AGRID(K) == GRID_NUM) THEN
                     PLY_NUM = 0                           ! 'N' in call to EMG means do not write to BUG file
                     CALL EMG ( IELE, OPT, 'N', SUBR_NAME, 'N' )

                     CALL ELMDIS                           ! Get local displ, UEL, for this element
                                                           ! Calc elem loads, PEL.
                     CALL CALC_ELEM_NODE_FORCES            ! JTSUB will only be used if SUBLOD(INT_SC_NUM,2) > 0
                                                           ! Transform elem loads from local to global coords
                     CALL TRANSFORM_NODE_FORCES (COORD_SYS)! Now we have elem node forces in global coords for the elem that has
!                                                            internal grid K (actual grid GRID_BUM) for this elem
                     DO L=1,6
                        BEG_ROW_NUM = 6*(K-1) + 1
                        PEG1(L) = PEG(BEG_ROW_NUM + L - 1) ! PEG are elem forces at elem ends in global coords
                        IF(PEG1(L) == ZERO) THEN;  PEG1(L) = -ZERO;  ENDIF  ! Avoids writing -0.0 for -PEG1(J) below

                        TOTALS(L) = TOTALS(L) - PEG1(L)
                     ENDDO

                     IF(WRITE_OP2) THEN
                       GPFORCE_NID_EID(INODE_GPFORCE,1) = GRID_NUM
                       GPFORCE_NID_EID(INODE_GPFORCE,2) = EID
                       GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,1) = PEG1(1)
                       GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,2) = PEG1(2)
                       GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,3) = PEG1(3)
                       GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,4) = PEG1(4)
                       GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,5) = PEG1(5)
                       GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,6) = PEG1(6)
                       GPFORCE_ETYPE(INODE_GPFORCE) = TYPE
                     ENDIF
                     INODE_GPFORCE = INODE_GPFORCE + 1
                     IF (WRITE_F06) WRITE(F06,9209) TYPE, EID, (-PEG1(L),L=1,6)  ! element forces
                  ENDIF
               ENDDO
            ENDDO

            IF(WRITE_OP2) THEN
              GPFORCE_NID_EID(INODE_GPFORCE,1) = GRID_NUM
              GPFORCE_NID_EID(INODE_GPFORCE,2) = 0
              GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,1) = TOTALS(1)
              GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,2) = TOTALS(2)
              GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,3) = TOTALS(3)
              GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,4) = TOTALS(4)
              GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,5) = TOTALS(5)
              GPFORCE_FXYZ_MXYZ(INODE_GPFORCE,6) = TOTALS(6)
              GPFORCE_ETYPE(INODE_GPFORCE) = "*TOTALS*"
            ENDIF
            INODE_GPFORCE = INODE_GPFORCE + 1

            IF (WRITE_F06) THEN
              WRITE(F06,9210)
              IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
                 IF (NDOFO == 0) THEN
                    WRITE(F06,9211) (TOTALS(J),J=1,6)
                 ELSE
                    WRITE(F06,9310) (TOTALS(J),J=1,6)
                 ENDIF
              ELSE
                 WRITE(F06,9211) (TOTALS(J),J=1,6)     ! KEEP THIS
              ENDIF
            ENDIF

            IF (IS_GPFORCE_SUMMARY_INFO) THEN
               DO J=1,6
                  IF (DABS(TOTALS(J)) > MAX_ABS_ALL_GRDS(J)) THEN
                     MAX_ABS_ALL_GRDS(J) = DABS(TOTALS(J))
                     MAX_ABS_GRID(J) = GRID(I,1)
                  ENDIF
               ENDDO
            ENDIF

         ENDIF

         ! For each of the 6 components (J=1,6 for components T1, T2, T3, R1, R2, R3),
         ! calc % of grid force imbalance as a % of the largest
         ! force item in that component
         DO J=1,6
            IF (DABS( FG1(J)) > MAX_ABS(J)) MAX_ABS(J) = DABS( FG1(J))
            IF (DABS( PG1(J)) > MAX_ABS(J)) MAX_ABS(J) = DABS( PG1(J))
            IF (DABS(PTET(J)) > MAX_ABS(J)) MAX_ABS(J) = DABS(PTET(J))
            IF (DABS(QGs1(J)) > MAX_ABS(J)) MAX_ABS(J) = DABS(QGs1(J))
            IF (DABS(QGm1(J)) > MAX_ABS(J)) MAX_ABS(J) = DABS(QGm1(J))
            IF (DABS(PEG1(J)) > MAX_ABS(J)) MAX_ABS(J) = DABS(PEG1(J))
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO
      WRITE(SC1,*) CR13
      IF (WRITE_F06) THEN
          CALL CALCULATE_GPFB_IMBALANCE(CHAR_PCT, MAX_ABS, MAX_ABS_PCT, MAX_ABS_GRID, MAX_ABS_ALL_GRDS)
      ENDIF

      !----------------
      ! DEALLOCATE: GPFORCE_NID_EID, GPFORCE_ETYPE, GPFORCE_FXYZ_MXYZ
      !ref mystran SUB DEALLOCATE_DOF_TABLES
      !KTSTACK(5500,3)

      IF(WRITE_OP2) THEN
          ! write the data to the op2
          !(nid_device, eid, elem_name, f1, f2, f3, m1, m2, m3) = out
          !nid = nid_device // 10

          CALL OUTPUT2_WRITE_OGF(ISUBCASE, INODE_GPFORCE-1, &
                                 TITLEI, SUBTITLEI, LABELI, &
                                 ANALYSIS_CODE, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
          WRITE(OP2) (GPFORCE_NID_EID(I,1)*10+DEVICE_CODE, &
                      GPFORCE_NID_EID(I,2), &
                      GPFORCE_ETYPE(I), &
                      REAL(GPFORCE_FXYZ_MXYZ(I,1), 4), &
                      REAL(GPFORCE_FXYZ_MXYZ(I,2), 4), &
                      REAL(GPFORCE_FXYZ_MXYZ(I,3), 4), &
                      REAL(GPFORCE_FXYZ_MXYZ(I,4), 4), &
                      REAL(GPFORCE_FXYZ_MXYZ(I,5), 4), &
                      REAL(GPFORCE_FXYZ_MXYZ(I,6), 4), &
                      I=1,INODE_GPFORCE-1)

          ITABLE = -5
          CALL END_OP2_TABLE(ITABLE)   ! close the previous

          !---------------------
          ! deallocate the arrays

          DEALLOCATE(GPFORCE_NID_EID,STAT=IERR)
          IF (IERR /= 0) THEN
              WRITE(ERR,*) 'GPFORCE_NID_EID DEALLOCATE error'
          ENDIF

          ! QUAD4, TRIA3, BAR
          DEALLOCATE (GPFORCE_ETYPE,STAT=IERR)
          IF (IERR /= 0) THEN
              WRITE(ERR,*) 'GPFORCE_ETYPE DEALLOCATE error'
          ENDIF

          DEALLOCATE (GPFORCE_FXYZ_MXYZ,STAT=IERR)
          IF (IERR /= 0) THEN
              WRITE(ERR,*) 'GPFORCE_FXYZ_MXYZ DEALLOCATE error'
          ENDIF
      ELSE
          ! deallocate the dummy arrays
          DEALLOCATE (GPFORCE_NID_EID)
          DEALLOCATE (GPFORCE_ETYPE)
          DEALLOCATE (GPFORCE_FXYZ_MXYZ)
      ENDIF
      FLUSH(OP2)
      FLUSH(F06)
      FLUSH(ERR)



      RETURN

! **********************************************************************************************************************************
 9101 FORMAT(' OUTPUT FOR SUBCASE ',I8)

 9102 FORMAT(' OUTPUT FOR EIGENVECTOR ',I8)

 9103 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (boundary ',A,' for grid',I8,' component',I2,')')

 9105 FORMAT(' OUTPUT FOR CRAIG-BAMPTON DOF ',I8,' OF ',I8,' (modal acceleration for mode ',I8,')')

 9200 FORMAT(//,1X,'                                          G R I D   P O I N T   F O R C E   B A L A N C E',/,                  &
             1X,   '                                            (in global coordinate system at each grid)',/)

 8999 FORMAT(//,1X,'                                   C B   G R I D   P O I N T   F O R C E   B A L A N C E   O T M',/,           &
             1X,   '                                            (in global coordinate system at each grid)',/)

 9201 FORMAT(1X,   '                               FORCE BALANCE FOR GRID POINT ', I8, ' IN GLOBAL COORD SYSTEM ',I8,/)

 9202 FORMAT(1X,   '                              T1            T2            T3            R1            R2            R3',/)

 9203 FORMAT(1X,   'APPLIED FORCE          ',6(1ES14.6))

 9204 FORMAT(1X,   'MINUS EQUIV THERM LOAD ',6(1ES14.6))

 9205 FORMAT(1X,   'SPC FORCE              ',6(1ES14.6))

 9206 FORMAT(1X,   'MPC FORCE              ',6(1ES14.6))

 9207 FORMAT(1X,   'INERTIA FORCE          ',6(1ES14.6))

 9208 FORMAT(1X,   'SUBSTRUCTURE I/F FORCE ',6(1ES14.6))

 9209 FORMAT(1X,      A8  ,' ELEM ' ,I9     ,6(1ES14.6))

 9210 FORMAT(1X,   '                        ------------- ------------- ------------- ------------- ------------- -------------')

 9211 FORMAT(1X,   'TOTALS                :',6(1ES14.6),/,                                                                         &
             1X,   '(should all be 0)',//)

 9310 FORMAT(1X,   'TOTALS                :',6(1ES14.6),13X,/,                                                                     &
             1X,   '(may not be zero since there were OMIT''d DOF''s which can mean that the correct inertia forces at the G-set', &
                   ' are not -Mgg*a*g)',//)

 9799 FORMAT(1X,A)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE CHK_COL_ARRAYS_ZEROS

      IMPLICIT NONE

      INTEGER(LONG)                   :: II                ! DO loop index


      DO II=1,NDOFG

         IF (ALLOCATED( FG_COL)) THEN
            IF ( FG_COL(II) == -ZERO) THEN;    FG_COL(II) = ZERO;   ENDIF
         ENDIF

         IF (ALLOCATED( PG_COL)) THEN
            IF ( PG_COL(II) == -ZERO) THEN;    PG_COL(II) = ZERO;   ENDIF
         ENDIF

         IF (ALLOCATED(QGs_COL)) THEN
            IF (QGs_COL(II) == -ZERO) THEN;   QGs_COL(II) = ZERO;   ENDIF
         ENDIF

         IF (ALLOCATED(QGm_COL)) THEN
            IF (QGm_COL(II) == -ZERO) THEN;   QGm_COL(II) = ZERO;   ENDIF
         ENDIF

         IF (ALLOCATED(QGr_COL)) THEN
            IF (QGr_COL(II) == -ZERO) THEN;   QGr_COL(II) = ZERO;   ENDIF
         ENDIF

         IF (ALLOCATED( UG_COL)) THEN
            IF ( UG_COL(II) == -ZERO) THEN;    UG_COL(II) = ZERO;   ENDIF
         ENDIF

      ENDDO

      END SUBROUTINE CHK_COL_ARRAYS_ZEROS

      !------------------------------------------------------
      LOGICAL FUNCTION IS_ABS_POSITIVE(ARRAY)
      ! IS_ABS_POSITIVE has an implicit return that has the same name
      ! as the function

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      IMPLICIT NONE
      REAL(DOUBLE), INTENT(IN)        :: ARRAY(6)         ! The 6 values
      INTEGER(LONG)                   :: I                ! counter

      IS_ABS_POSITIVE = .FALSE.
      DO I=1,6
          IF (ABS(ARRAY(I)) > 0.0) THEN
              IS_ABS_POSITIVE = .TRUE.
          ENDIF
      ENDDO
      END FUNCTION IS_ABS_POSITIVE
      !------------------------------------------------------

      SUBROUTINE CALCULATE_GPFB_IMBALANCE(CHAR_PCT, MAX_ABS, MAX_ABS_PCT, MAX_ABS_GRID, MAX_ABS_ALL_GRDS)
      ! For each of the 6 components (J=1,6 for components T1, T2, T3, R1, R2, R3)
      !  - calc % of grid force imbalance as a % of the largest
      !    force item in that component
      USE PENTIUM_II_KIND, ONLY       :  BYTE, SHORT, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      IMPLICIT NONE
      CHARACTER(13*BYTE)              :: CHAR_PCT(6)        ! Character representation of MEFFMASS sum percents of total model mass
      INTEGER(LONG)                   :: I                  ! DO loop index

      REAL(DOUBLE), INTENT(IN)        :: MAX_ABS(6)         ! The 6 abs largest values of any force item in the 6 components
      REAL(DOUBLE), INTENT(INOUT)     :: MAX_ABS_PCT(6)     ! Modal mass as % of total mass

      !INTRINSIC IAND
      INTEGER(LONG), INTENT(IN)       :: max_abs_grid(6)    ! Grid where max-abs for GP force balance totals exists
      REAL(DOUBLE), INTENT(IN)        :: max_abs_all_grds(6)! The 6 max-abs from GP force balance totals

      IF (IS_GPFORCE_SUMMARY_INFO) THEN
         DO I=1,6
            IF (MAX_ABS(I) > ZERO) THEN
               MAX_ABS_PCT(I) = ONE_HUNDRED*MAX_ABS_ALL_GRDS(I)/MAX_ABS(I)
               WRITE(CHAR_PCT(I),9213) MAX_ABS_PCT(I)
            ELSE                                                 ! Denominator is zero so leave % blank
               CHAR_PCT(I)(1:) = ' '
            ENDIF
         ENDDO
         IF (WRITE_F06) THEN
            WRITE(F06,9214)
            WRITE(F06,9202)
            WRITE(F06,9215) (MAX_ABS_ALL_GRDS(I),I=1,6)
            IF (DEBUG(192) > 1) THEN
               WRITE(F06,9216) (MAX_ABS(I),I=1,6)
               WRITE(F06,9217) (CHAR_PCT(I),I=1,6)
            ENDIF
            WRITE(F06,9218) (MAX_ABS_GRID(I),I=1,6)
         ENDIF
      ENDIF

 9202 FORMAT(1X,   '                              T1            T2            T3            R1            R2            R3',/)

 9213 FORMAT(1ES12.2,'%')

 9214 FORMAT(1X,   '                                     Max abs values of force imbalance totals from above grids',/)

 9215 FORMAT(1X,'Max abs imbal any grid:',6(1ES14.6))

 9216 FORMAT(1X,'Max abs force any grid:',6(1ES14.6))

 9217 FORMAT(1X,'as % of max abs force : ',6(A14))

 9218 FORMAT(1X,   'Occurs at grid*       :',2X,I8,5(6X,I8),/,                                                                     &
             1X,   '(*for output set)')

      END SUBROUTINE CALCULATE_GPFB_IMBALANCE
!===================================================================================================================================

      SUBROUTINE OUTPUT2_WRITE_OGF(ISUBCASE, NUM, TITLE, SUBTITLE, LABEL, &
                                   ANALYSIS_CODE, FIELD5_INT_MODE, FIELD6_EIGENVALUE)
!     writes the grid point force results header.
!     Data is first written to character variables and then that character variable is output the F06.
!
!     Parameters
!     ==========
!
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  OP2, ERR
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_Value, IEEE_QUIET_NAN
      USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL32

      INTEGER(LONG), INTENT(IN)      :: ISUBCASE           ! subcase id
      INTEGER(LONG), INTENT(IN)      :: NUM                ! the number of "nodes" that will be printed
      CHARACTER(LEN=128), INTENT(IN) :: TITLE              ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN) :: SUBTITLE           ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN) :: LABEL              ! the subcase LABEL
      INTEGER(LONG), INTENT(IN) :: ANALYSIS_CODE           ! static, time, frequency, modal, etc. flag
      INTEGER(LONG), INTENT(IN) :: FIELD5_INT_MODE
      REAL(DOUBLE),  INTENT(IN) :: FIELD6_EIGENVALUE

      CHARACTER*8   :: TABLE_NAME
      INTEGER(LONG) :: ITABLE       ! the current subtable number
      INTEGER(LONG) :: NUM_WIDE     ! the number of "words" for a single element
      INTEGER(LONG) :: DEVICE_CODE  ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG) :: NVALUES      ! the number of "words" for the OP2 data
      INTEGER(LONG) :: NTOTAL       ! the number of bytes corresponding to nvalues
      REAL(REAL32)  :: NAN
      LOGICAL       :: IS_PRINT     ! is this a PRINT result -> F06

      ! TODO: assuming PLOT
      DEVICE_CODE = 1

      !-------
      ! we're already setup for subtable -3

      ! nid, eid, ETYPE (8-bytes), fx, fy, fz, mx, my, mz
      NUM_WIDE = 10

      ! dunno???
      TABLE_NAME = 'OGPFB1  '
      CALL WRITE_TABLE_HEADER(TABLE_NAME)
      ITABLE = -3

      !WRITE(ERR,*) "OUTPUT2_WRITE_OGF: ISUBCASE=",ISUBCASE
      !WRITE(ERR,*) "OUTPUT2_WRITE_OGF: ANALYSIS_CODE=",ANALYSIS_CODE
      CALL WRITE_OGF3(ITABLE, ISUBCASE, ANALYSIS_CODE, DEVICE_CODE, NUM_WIDE, &
                      TITLE, SUBTITLE, LABEL,                                          &
                      FIELD5_INT_MODE, FIELD6_EIGENVALUE)
      ! ITABLE = -4, -6, ...
      !NWORDS = NUM * NUM_WIDE
      !NTOTAL = NBYTES_PER_WORD * NWORDS

! 100  FORMAT("*DEBUG:    ITABLE=",I8, "; NUM=",I8,"; NVALUES=",I8,"; NTOTAL=",I8)
      NVALUES = NUM * NUM_WIDE
      NTOTAL = NVALUES * 4
      !WRITE(ERR,100) ITABLE,NUM,NVALUES,NTOTAL
      WRITE(OP2) NVALUES

! 102  FORMAT("*DEBUG: OUTPUT2_WRITE_OGF      ITABLE=",I8," (should be -5, -7,...)")
      !WRITE(ERR,102) ITABLE

      END SUBROUTINE OUTPUT2_WRITE_OGF
!==================================================================================================
      SUBROUTINE WRITE_OGF3(ITABLE, ISUBCASE, ANALYSIS_CODE, DEVICE_CODE, NUM_WIDE, &
                            TITLE, LABEL, SUBTITLE, &
                            FIELD5_INT_MODE, FIELD6_EIGENVALUE)
!      Parameters
!      ==========
!      analysis_code
!        the solution type flag
!      approach_code
!        ???
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR,OP2
      IMPLICIT NONE
      INTEGER(LONG), INTENT(INOUT) :: ITABLE               ! an OP2 subtable counter
      INTEGER(LONG), INTENT(IN) :: ANALYSIS_CODE           ! static, modal, time, freq, etc. flag
      INTEGER(LONG), INTENT(IN) :: ISUBCASE                ! the subcase ID
      INTEGER(LONG), INTENT(IN) :: DEVICE_CODE             ! PLOT, PRINT, PUNCH flag
      INTEGER(LONG), INTENT(IN) :: NUM_WIDE                ! the number of words per element

      CHARACTER(LEN=128), INTENT(IN) :: TITLE              ! the model TITLE
      CHARACTER(LEN=128), INTENT(IN) :: SUBTITLE           ! the subcase SUBTITLE
      CHARACTER(LEN=128), INTENT(IN) :: LABEL              ! the subcase LABEL

      CHARACTER(LEN=128) :: TITLE2                 ! the subcase TITLE
      CHARACTER(LEN=128) :: SUBTITLE2              ! the subcase SUBTITLE
      CHARACTER(LEN=128) :: LABEL2                 ! the subcase LABEL

      INTEGER(LONG), INTENT(IN) :: FIELD5_INT_MODE
      REAL(DOUBLE), INTENT(IN)  :: FIELD6_EIGENVALUE
      REAL(DOUBLE)              :: FIELD7
      INTEGER(LONG)             :: APPROACH_CODE
      INTEGER(LONG)             :: FORMAT_CODE        ! 1=real, 2=imag, 3=random
      INTEGER(LONG)             :: TABLE_CODE, LOAD_SET, THERMAL, ACOUSTIC_FLAG, STRESS_CODE

!     static is real
      FORMAT_CODE = 1

      TITLE2 = TITLE(1:100)
      SUBTITLE2 = SUBTITLE(1:67)
      LABEL2 = LABEL(1:100)

      CALL WRITE_ITABLE(ITABLE)  ! write the -3, -5, ... subtable header
! 1    FORMAT("WRITE_OGF3: ITABLE_START=",I8)
      !WRITE(ERR,1) ITABLE

      IF ((ANALYSIS_CODE == 1) .OR. (ANALYSIS_CODE == 10)) THEN
         ! statics
         FIELD7 = 0.0
      ELSE
         ! frequency in radians
         FIELD7 = SQRT(ABS(FIELD6_EIGENVALUE))
      ENDIF

      WRITE(OP2) 146
      ! grid point forces only
      TABLE_CODE = 19

      ! ???
      LOAD_SET = 1

      ! we're not doing acoustic/stress
      ACOUSTIC_FLAG = 0
      STRESS_CODE = 0

      ! not always 0 for stress, but for now
      THERMAL = 0

      APPROACH_CODE = ANALYSIS_CODE * 10 + DEVICE_CODE
!2     FORMAT(" APPROACH_CODE=",I4," TABLE_CODE=",I4," ISUBCASE=",I4)
      ! 584 bytes
      !WRITE(ERR,2) APPROACH_CODE, TABLE_CODE, ISUBCASE
      WRITE(OP2) APPROACH_CODE, TABLE_CODE, 0, ISUBCASE, FIELD5_INT_MODE, &
            REAL(FIELD6_EIGENVALUE, 4), REAL(FIELD7, 4),                          &
            LOAD_SET, FORMAT_CODE, NUM_WIDE, &
            STRESS_CODE, ACOUSTIC_FLAG, 0, 0, 0, &
            0, 0, 0, 0, 0, &
            0, 0, THERMAL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
            0, 0, 0, 0, &
            TITLE2, SUBTITLE2, LABEL2

      ITABLE = ITABLE - 1        ! flip it to -4, -6, ... so we don't have to do this later
! 3    FORMAT("WRITE_OGF3: ITABLE_END=",I8)
      !WRITE(ERR,3) ITABLE
      CALL WRITE_ITABLE(ITABLE)
      ITABLE = ITABLE - 1
      END SUBROUTINE WRITE_OGF3

      END SUBROUTINE GP_FORCE_BALANCE_PROC
