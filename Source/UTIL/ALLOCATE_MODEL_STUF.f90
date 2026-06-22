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

      SUBROUTINE ALLOCATE_MODEL_STUF ( NAME_IN, CALLING_SUBR )

! Allocates arrays defined in module MODEL_STUF (see comments there for definition of these matrices)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, THREE, SIX, ONEPP6
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE SCONTR, ONLY                :  LBAROFF, LBUSHOFF, LCMASS, LCONM2, LCORD, LEDAT, LELE, LFORCE, LGRAV, LGRID,              &
                                         LIND_GRDS_MPCS, LLOADC, LLOADR, LMATANGLE, LMATL, LMPC, LMPCADDC, LMPCADDR, LPBAR, LPBEAM,&
                                         LPBUSH, LPCOMP, LPCOMP_PLIES, LPDAT, LPELAS, LPLATEOFF, LPLATETHICK, LPLOAD, LPMASS,      &
                                         LPROD, LPSHEAR, LPSHEL, LPSOLID, LPUSER1, LPUSERIN, LRFORCE, LRIGEL, LSEQ, LSETLN, LSETS, &
                                         LSLOAD, LSNORM, LSPC, LSPC1, LSPCADDC, LSPCADDR, LSUB, LTDAT, LVVEC
      USE SCONTR, ONLY                :  MAX_ELEM_DEGREE, MAX_GAUSS_POINTS, MAX_STRESS_POINTS, MCMASS, MCONM2, MCORD, MDT, MELDOF, &
                                         MELGP, MGRID, MMATL, MOFFSET, MPBAR, MPBEAM, MPBUSH, MPCOMP_PLIES, MPCOMP0, MPELAS,       &
                                         MPMASS, MPRESS, MPLOAD4_3D_DATA, MPROD, MPSHEAR, MPSHEL, MPSOLID, MPUSER1, MPUSERIN,      &
                                         MRCONM2, MRCORD,                                                                          &
                                         MRGRID, MRMATLC, MRPBAR, MRPBEAM, MRPBUSH, MRPCOMP_PLIES, MRPCOMP0, MRPELAS, MRPMASS,     &
                                         MRPROD, MRPSHEAR, MRPSHEL, MRPUSER1, MRSNORM, MSNORM, MGRID_SNORM, MUSERIN_MAT_NAMES  
! --- cbeam_add begin --- !
      USE SCONTR, ONLY                :  MPBEAM_STATIONS
! --- cbeam_add end--- !                                         
      USE SCONTR, ONLY                :  NDOFG, NGRID, NMPC, NPCOMP, NPLOAD4_3D, NRBAR, NRBE1, NRBE2, NSPC, NTSUB, NUM_MPCSIDS,    &
                                         NUM_SPCSIDS
      USE TIMDAT, ONLY                :  TSEC

      USE MODEL_STUF, ONLY            :  AGRID, BE1, BE2, BE3, BGRID, DOFPIN, DT, ME, OFFDIS, OFFDIS_O, OFFDIS_B, OFFDIS_G, OFFSET,&
                                         KE, KEG, KED, KEM, PEB, PEG, PEL, PPE, PRESS, PTE, SE1, SE2, SE3, STE1, STE2, STE3,       &
                                         UEB, UEG, UEL, UGG, USERIN_NUM_ACT_GRDS, USERIN_ACT_COMPS, USERIN_ACT_GRIDS,              &
                                         USERIN_MAT_NAMES, XEB, XEL, XGL
      USE MODEL_STUF, ONLY            :  CMASS, CONM2, PMASS, RPMASS, RCONM2
      USE MODEL_STUF, ONLY            :  CORD, RCORD, TN
      USE MODEL_STUF, ONLY            :  SEQ1, SEQ2
      USE MODEL_STUF, ONLY            :  BAROFF, BUSHOFF, EDAT, EOFF, EPNT, ESORT1, ESORT2, ETYPE, PLATEOFF, PLATETHICK, VVEC
      USE MODEL_STUF, ONLY            :  PRESS_SIDS, PDATA, PPNT, PTYPE
      USE MODEL_STUF, ONLY            :  PLOAD4_3D_DATA
      USE MODEL_STUF, ONLY            :  FORMOM_SIDS
      USE MODEL_STUF, ONLY            :  GRAV_SIDS, RFORCE_SIDS, SLOAD_SIDS
      USE MODEL_STUF, ONLY            :  GRID, RGRID
      USE MODEL_STUF, ONLY            :  GRID_ID, GRID_SEQ, INV_GRID_SEQ, MPC_IND_GRIDS
      USE MODEL_STUF, ONLY            :  SNORM, RSNORM, GRID_SNORM
      USE MODEL_STUF, ONLY            :  MATL, RMATL, PBAR, RPBAR, PBEAM, RPBEAM, PBUSH, RPBUSH, PCOMP, RPCOMP, PELAS, RPELAS,     &
                                         PROD, RPROD, PSHEAR, RPSHEAR, PSHEL, RPSHEL, PSOLID, PUSER1, RPUSER1, PUSERIN,            &
                                         USERIN_ACT_COMPS, USERIN_ACT_GRIDS, USERIN_MAT_NAMES
! --- cbeam_add begin --- !
      USE MODEL_STUF, ONLY            :  PBEAM_NSTATIONS, PBEAM_XL, PBEAM_RPROPS
! --- cbeam_add end --- !                                                                                 
      USE MODEL_STUF, ONLY            :  MPC_SIDS, MPCSIDS, MPCADD_SIDS
      USE MODEL_STUF, ONLY            :  SPC_SIDS, SPC1_SIDS, SPCSIDS, SPCADD_SIDS
      USE MODEL_STUF, ONLY            :  ALL_SETS_ARRAY, ONE_SET_ARRAY, SETS_IDS, SC_ACCE, SC_DISP, SC_ELFN, SC_ELFE, SC_GPFO,     &
                                         SC_MPCF, SC_OLOA, SC_SPCF, SC_STRE, SC_STRN, LOAD_SIDS, LOAD_FACS
      USE MODEL_STUF, ONLY            :  ELDT, ELOUT, GROUT, OELOUT, OGROUT, LABEL, SCNUM, STITLE, SUBLOD, TITLE
      USE MODEL_STUF, ONLY            :  CC_EIGR_SID_SUB, CC_STATSUB_SUB, EIG_PARAMS, EIG_PARAMS_TYPE, IS_BUCKLING_SUBCASE,    &
                                         IS_MODES_SUBCASE, NUM_EIGENS_SUB
      USE MODEL_STUF, ONLY            :  SYS_LOAD
      USE MODEL_STUF, ONLY            :  CETEMP, CETEMP_ERR, CGTEMP, CGTEMP_ERR, ETEMP, ETEMP_INIT, GTEMP, GTEMP_INIT, TDATA, TPNT
      USE MODEL_STUF, ONLY            :  RIGID_ELEM_IDS
      USE MODEL_STUF, ONLY            :  MATANGLE
      USE MODEL_STUF, ONLY            :  SPCSETS, MPCSETS
      USE MODEL_STUF, ONLY            :  GRID_ELEM_CONN_ARRAY

      USE ALLOCATE_MODEL_STUF_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ALLOCATE_MODEL_STUF'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME_IN           ! Name of group of arrays to allocate
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Array name of the matrix to be allocated in sparse format
      CHARACTER(31*BYTE)              :: NAME              ! Specific array name used for output error message

      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator
      INTEGER(LONG)                   :: NCOLS             ! Number of cols allocated
      INTEGER(LONG)                   :: NROWS             ! Number of rows allocated


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM)
      REAL(DOUBLE)                    :: MB_ALLOCATED      ! Megabytes of mmemory allocated for the arrays to put into array
!                                                            ALLOCATED_ARRAY_MEM when subr ALLOCATED_MEMORY is called
      REAL(DOUBLE)                    :: RBYTE             ! Real value of BYTE
      REAL(DOUBLE)                    :: RDOUBLE           ! Real value of DOUBLE
      REAL(DOUBLE)                    :: RLONG             ! Real value of LONG

      INTRINSIC                       :: REAL



! **********************************************************************************************************************************
      RBYTE   = REAL(BYTE)
      RDOUBLE = REAL(DOUBLE)
      RLONG   = REAL(LONG)

      MB_ALLOCATED = ZERO
      JERR = 0

      IF (NAME_IN == 'SETS ARRAYS') THEN                      ! Allocate arrays for SET's

         NAME = 'SETS_IDS'
         IF (ALLOCATED(SETS_IDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SETS_IDS(LSETS),STAT=IERR)
            NROWS = LSETS
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSETS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSETS
                  SETS_IDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'ALL_SETS_ARRAY'
         IF (ALLOCATED(ALL_SETS_ARRAY)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (ALL_SETS_ARRAY(LSETLN),STAT=IERR)
            NROWS = LSETLN
            NCOLS = 1
            MB_ALLOCATED = RBYTE*REAL(LEN(ALL_SETS_ARRAY))*REAL(LSETLN)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSETLN
                  ALL_SETS_ARRAY(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'ONE_SET_ARRAY'
         IF (ALLOCATED(ONE_SET_ARRAY)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (ONE_SET_ARRAY(LSETLN),STAT=IERR)
            NROWS = LSETLN
            NCOLS = 1
            MB_ALLOCATED = RBYTE*REAL(LEN(ONE_SET_ARRAY)*LSETLN)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSETLN
                  ONE_SET_ARRAY(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'TITLES') THEN                      ! Allocate arrays for output title, subtitle, label

         NAME = 'TITLE'
         IF (ALLOCATED(TITLE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TITLE(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RBYTE*REAL(LEN(TITLE)*LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  TITLE(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'STITLE'
         IF (ALLOCATED(STITLE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (STITLE(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RBYTE*REAL(LEN(STITLE)*LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  STITLE(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'LABEL'
         IF (ALLOCATED(LABEL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (LABEL(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RBYTE*REAL(LEN(LABEL)*LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  LABEL(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SC_xxxx') THEN                     ! Allocate arrays for SC_DIDP, etc

         NAME = 'SC_ACCE'
         IF (ALLOCATED(SC_ACCE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SC_ACCE(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SC_ACCE = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_DISP'
         IF (ALLOCATED(SC_DISP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SC_DISP(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SC_DISP = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_ELFE'
         IF (ALLOCATED(SC_ELFE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SC_ELFE(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SC_ELFE = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_ELFN'
         IF (ALLOCATED(SC_ELFN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SC_ELFN(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SC_ELFN = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_GPFO'
         IF (ALLOCATED(SC_GPFO)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SC_GPFO(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SC_GPFO = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_MPCF'
         IF (ALLOCATED(SC_MPCF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SC_MPCF(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SC_MPCF = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_OLOA'
         IF (ALLOCATED(SC_OLOA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SC_OLOA(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SC_OLOA = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_SPCF'
         IF (ALLOCATED(SC_SPCF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SC_SPCF(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SC_SPCF = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_STRE'
         IF (ALLOCATED(SC_STRE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SC_STRE(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SC_STRE = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_STRN'
         IF (ALLOCATED(SC_STRN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SC_STRN(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SC_STRN = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SCNUM') THEN                       ! Allocate array SCNUM

         NAME = 'SCNUM'
         IF (ALLOCATED(SCNUM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SCNUM(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SCNUM(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         ! Piggy-back the per-subcase METHOD/eigen-parameter arrays on SCNUM. Their lifetime exactly matches SCNUM's,
         ! and any code path that allocates SCNUM needs these too. The per-mode result arrays (EIGEN_VAL_SUB, etc.)
         ! are allocated separately by LINK4 once NDOFL and the per-subcase mode counts are known.
         IF (.NOT. ALLOCATED(CC_EIGR_SID_SUB)) THEN
            ALLOCATE (CC_EIGR_SID_SUB(LSUB), STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,LSUB
                  CC_EIGR_SID_SUB(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) 0.0D0, 'CC_EIGR_SID_SUB', SUBR_NAME, IERR
               WRITE(F06,991) 0.0D0, 'CC_EIGR_SID_SUB', SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         IF (.NOT. ALLOCATED(IS_MODES_SUBCASE)) THEN
            ALLOCATE (IS_MODES_SUBCASE(LSUB), STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,LSUB
                  IS_MODES_SUBCASE(I) = 'N'
               ENDDO
            ELSE
               WRITE(ERR,991) 0.0D0, 'IS_MODES_SUBCASE', SUBR_NAME, IERR
               WRITE(F06,991) 0.0D0, 'IS_MODES_SUBCASE', SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         IF (.NOT. ALLOCATED(CC_STATSUB_SUB)) THEN
            ALLOCATE (CC_STATSUB_SUB(LSUB), STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,LSUB
                  CC_STATSUB_SUB(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) 0.0D0, 'CC_STATSUB_SUB', SUBR_NAME, IERR
               WRITE(F06,991) 0.0D0, 'CC_STATSUB_SUB', SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         IF (.NOT. ALLOCATED(IS_BUCKLING_SUBCASE)) THEN
            ALLOCATE (IS_BUCKLING_SUBCASE(LSUB), STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,LSUB
                  IS_BUCKLING_SUBCASE(I) = 'N'
               ENDDO
            ELSE
               WRITE(ERR,991) 0.0D0, 'IS_BUCKLING_SUBCASE', SUBR_NAME, IERR
               WRITE(F06,991) 0.0D0, 'IS_BUCKLING_SUBCASE', SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         IF (.NOT. ALLOCATED(NUM_EIGENS_SUB)) THEN
            ALLOCATE (NUM_EIGENS_SUB(LSUB), STAT=IERR)
            IF (IERR == 0) THEN
               DO I=1,LSUB
                  NUM_EIGENS_SUB(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) 0.0D0, 'NUM_EIGENS_SUB', SUBR_NAME, IERR
               WRITE(F06,991) 0.0D0, 'NUM_EIGENS_SUB', SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         IF (.NOT. ALLOCATED(EIG_PARAMS)) THEN
            ALLOCATE (EIG_PARAMS(LSUB), STAT=IERR)
            IF (IERR /= 0) THEN
               WRITE(ERR,991) 0.0D0, 'EIG_PARAMS', SUBR_NAME, IERR
               WRITE(F06,991) 0.0D0, 'EIG_PARAMS', SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
            ! Defaults from EIG_PARAMS_TYPE component initializers are applied automatically.
         ENDIF


      ELSE IF (NAME_IN == 'SUBLOD') THEN                      ! Allocate array SUBLOD

         NAME = 'SUBLOD'
         IF (ALLOCATED(SUBLOD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SUBLOD(LSUB,2),STAT=IERR)
            NROWS = LSUB
            NCOLS = 2
            MB_ALLOCATED = RLONG*REAL(LSUB)*TWO/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SUBLOD(I,1) = 0
                  SUBLOD(I,2) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SPC_MPC_SET') THEN                      ! Allocate array SUBLOD

         NAME = 'SPCSETS'
         IF (ALLOCATED(SPCSETS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SPCSETS(LSUB),STAT=IERR)
            NROWS = LSUB
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  SPCSETS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MPCSETS'
         IF (ALLOCATED(MPCSETS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MPCSETS(LSUB),STAT=IERR)
            NROWS = LSUB
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  MPCSETS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SEQ1,2') THEN                      ! Allocate arrays for SEQ1,2

         NAME = 'SEQ1'
         IF (ALLOCATED(SEQ1)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SEQ1(LSEQ),STAT=IERR)
            NROWS = LSEQ
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSEQ)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSEQ
                  SEQ1(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SEQ2'
         IF (ALLOCATED(SEQ2)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SEQ2(LSEQ),STAT=IERR)
            NCOLS = 1
            NROWS = LSEQ
            MB_ALLOCATED = RLONG*REAL(LSEQ)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSEQ
                  SEQ2(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'FORMOM_SIDS') THEN                 ! Allocate arrays for force/moment set ID's

         NAME = 'FORMOM_SIDS'
         IF (ALLOCATED(FORMOM_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (FORMOM_SIDS(LFORCE),STAT=IERR)
            NROWS = LFORCE
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LFORCE)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LFORCE
                  FORMOM_SIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PRESS_SIDS') THEN                  ! Allocate arrays for pressure load set ID's

         NAME = 'PRESS_SIDS'
         IF (ALLOCATED(PRESS_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PRESS_SIDS(LPLOAD),STAT=IERR)
            NROWS = LPLOAD
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LPLOAD)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPLOAD
                  PRESS_SIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRAV_SIDS') THEN                   ! Allocate arrays for grav load set ID's

         NAME = 'GRAV_SIDS'
         IF (ALLOCATED(GRAV_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GRAV_SIDS(LGRAV),STAT=IERR)
            NROWS = LGRAV
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LGRAV)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LGRAV
                  GRAV_SIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'LOAD_SIDS, LOAD_FACS') THEN        ! Allocate arrays for load set ID's and factors

         NAME = 'LOAD_SIDS'
         IF (ALLOCATED(LOAD_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (LOAD_SIDS(LLOADR,LLOADC),STAT=IERR)
            NROWS = LLOADR
            NCOLS = LLOADC
            MB_ALLOCATED = RLONG*REAL(LLOADR)*REAL(LLOADC)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LLOADR
                  DO J=1,LLOADC
                     LOAD_SIDS(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'LOAD_FACS'
         IF (ALLOCATED(LOAD_FACS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (LOAD_FACS(LLOADR,LLOADC),STAT=IERR)
            NROWS = LLOADR
            NCOLS = LLOADC
            MB_ALLOCATED = RDOUBLE*REAL(LLOADR)*REAL(LLOADC)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LLOADR
                  DO J=1,LLOADC
                     LOAD_FACS(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MPC_SIDS') THEN                    ! Allocate arrays for MPC set ID's

         NAME = 'MPC_SIDS'
         IF (ALLOCATED(MPC_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MPC_SIDS(LMPC),STAT=IERR)
            NROWS = LMPC
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LMPC)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LMPC
                  MPC_SIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MPCSIDS') THEN                     ! Allocate arrays for MPC set ID's used in 1 execution

         NAME = 'MPCSIDS'
         IF (ALLOCATED(MPCSIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MPCSIDS(NUM_MPCSIDS),STAT=IERR)
            NROWS = NUM_MPCSIDS
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(NUM_MPCSIDS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NUM_MPCSIDS
                  MPCSIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MPCADD_SIDS') THEN                 ! Allocate arrays for MPCADD set ID's

         NAME = 'MPCADD_SIDS'
         IF (ALLOCATED(MPCADD_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MPCADD_SIDS(LMPCADDR,LMPCADDC),STAT=IERR)
            NROWS = LMPCADDR
            NCOLS = LMPCADDC
            MB_ALLOCATED = RLONG*REAL(LMPCADDR)*REAL(LMPCADDC)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LMPCADDR
                  DO J=1,LMPCADDC
                     MPCADD_SIDS(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'RFORCE_SIDS') THEN                 ! Allocate arrays for grav load set ID's

         NAME = 'RFORCE_SIDS'
         IF (ALLOCATED(RFORCE_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RFORCE_SIDS(LRFORCE),STAT=IERR)
            NROWS = LRFORCE
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LRFORCE)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LRFORCE
                  RFORCE_SIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SLOAD_SIDS') THEN                  ! Allocate arrays for grav load set ID's

         NAME = 'SLOAD_SIDS'
         IF (ALLOCATED(SLOAD_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SLOAD_SIDS(LSLOAD),STAT=IERR)
            NROWS = LSLOAD
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSLOAD)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSLOAD
                  SLOAD_SIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SPC_SIDS, SPC1_SIDS') THEN         ! Allocate arrays for SPC, SPC1 set ID's

         NAME = 'SPC_SIDS'
         IF (ALLOCATED(SPC_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SPC_SIDS(LSPC),STAT=IERR)
            NROWS = LSPC
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSPC)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSPC
                  SPC_SIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SPC1_SIDS'
         IF (ALLOCATED(SPC1_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SPC1_SIDS(LSPC1),STAT=IERR)
            NROWS = LSPC1
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSPC1)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSPC1
                  SPC1_SIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SPCSIDS') THEN                     ! Allocate arrays for SPC set ID's used in 1 execution

         NAME = 'SPCSIDS'
         IF (ALLOCATED(SPCSIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SPCSIDS(NUM_SPCSIDS),STAT=IERR)
            NROWS = NUM_SPCSIDS
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(NUM_SPCSIDS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NUM_SPCSIDS
                  SPCSIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SPCADD_SIDS') THEN                 ! Allocate arrays for SPCADD set ID's

         NAME = 'SPCADD_SIDS'
         IF (ALLOCATED(SPCADD_SIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SPCADD_SIDS(LSPCADDR,LSPCADDC),STAT=IERR)
            NROWS = LSPCADDR
            NCOLS = LSPCADDC
            MB_ALLOCATED = RLONG*REAL(LSPCADDR)*REAL(LSPCADDC)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSPCADDR
                  DO J=1,LSPCADDC
                     SPCADD_SIDS(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ETYPE, EDAT, EPNT') THEN           ! Allocate arrays for ETYPE, EDAT, EPNT

         NAME = 'ETYPE'
         IF (ALLOCATED(ETYPE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (ETYPE(LELE),STAT=IERR)
            NROWS = LELE
            NCOLS = 1
            MB_ALLOCATED = RBYTE*REAL(LEN(ETYPE))*REAL(LELE)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  ETYPE(I)(1:) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'EDAT'
         IF (ALLOCATED(EDAT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (EDAT(LEDAT),STAT=IERR)
            NROWS = LEDAT
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LEDAT)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LEDAT
                  EDAT(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'EPNT'
         IF (ALLOCATED(EPNT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (EPNT(LELE+1),STAT=IERR)
            NROWS = LELE+1
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL((LELE+1))/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE+1
                  EPNT(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'EOFF') THEN                        ! Allocate array for EOFF

         NAME = 'EOFF'
         IF (ALLOCATED(EOFF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (EOFF(LELE),STAT=IERR)
            NROWS = LELE
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LELE)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  EOFF(I) = 'N'
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ELEM PROPERTIES AND MATERIALS') THEN

         NAME = 'MATL'
         IF (ALLOCATED(MATL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MATL(LMATL,MMATL),STAT=IERR)
            NROWS = LMATL
            NCOLS = MMATL
            MB_ALLOCATED = RLONG*REAL(LMATL)*REAL(MMATL)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LMATL
                  DO J=1,MMATL
                     MATL(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RMATL'
         IF (ALLOCATED(RMATL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RMATL(LMATL,MRMATLC),STAT=IERR)
            NROWS = LMATL
            NCOLS = MRMATLC
            MB_ALLOCATED = RDOUBLE*REAL(LMATL)*REAL(MRMATLC)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LMATL
                  DO J=1,MRMATLC
                     RMATL(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PBAR'
         IF (ALLOCATED(PBAR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PBAR(LPBAR,MPBAR),STAT=IERR)
            NROWS = LPBAR
            NCOLS = MPBAR
            MB_ALLOCATED = RLONG*REAL(LPBAR)*REAL(MPBAR)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPBAR
                  DO J=1,MPBAR
                     PBAR(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPBAR'
         IF (ALLOCATED(RPBAR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RPBAR(LPBAR,MRPBAR),STAT=IERR)
            NROWS = LPBAR
            NCOLS = MRPBAR
            MB_ALLOCATED = RDOUBLE*REAL(LPBAR)*REAL(MRPBAR)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPBAR
                  DO J=1,MRPBAR
                     RPBAR(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PBEAM'
         IF (ALLOCATED(PBEAM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PBEAM(LPBEAM,MPBEAM),STAT=IERR)
            NROWS = LPBEAM
            NCOLS = MPBEAM
            MB_ALLOCATED = RLONG*REAL(LPBEAM)*REAL(MPBEAM)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPBEAM
                  DO J=1,MPBEAM
                     PBEAM(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPBEAM'
         IF (ALLOCATED(RPBEAM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RPBEAM(LPBEAM,MRPBEAM),STAT=IERR)
            NROWS = LPBEAM
            NCOLS = MRPBEAM
            MB_ALLOCATED = RDOUBLE*REAL(LPBEAM)*REAL(MRPBEAM)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPBEAM
                  DO J=1,MRPBEAM
                     RPBEAM(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

! --- CBEAM_add begin --- !
         NAME = 'PBEAM_NSTATIONS'
         IF (ALLOCATED(PBEAM_NSTATIONS)) THEN
            CONTINUE
         ELSE
            ALLOCATE (PBEAM_NSTATIONS(LPBEAM),STAT=IERR)
            NROWS = LPBEAM
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LPBEAM)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPBEAM
                  PBEAM_NSTATIONS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PBEAM_XL'
         IF (ALLOCATED(PBEAM_XL)) THEN
            CONTINUE
         ELSE
            ALLOCATE (PBEAM_XL(LPBEAM,MPBEAM_STATIONS),STAT=IERR)
            NROWS = LPBEAM
            NCOLS = MPBEAM_STATIONS
            MB_ALLOCATED = RDOUBLE*REAL(LPBEAM)*REAL(MPBEAM_STATIONS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPBEAM
                  DO J=1,MPBEAM_STATIONS
                     PBEAM_XL(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PBEAM_RPROPS'
         IF (ALLOCATED(PBEAM_RPROPS)) THEN
            CONTINUE
         ELSE
            ALLOCATE (PBEAM_RPROPS(LPBEAM,MPBEAM_STATIONS,6),STAT=IERR)
            NROWS = LPBEAM
            NCOLS = MPBEAM_STATIONS*6
            MB_ALLOCATED = RDOUBLE*REAL(LPBEAM)*REAL(MPBEAM_STATIONS)*6.0/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPBEAM
                  DO J=1,MPBEAM_STATIONS
                     PBEAM_RPROPS(I,J,1) = ZERO
                     PBEAM_RPROPS(I,J,2) = ZERO
                     PBEAM_RPROPS(I,J,3) = ZERO
                     PBEAM_RPROPS(I,J,4) = ZERO
                     PBEAM_RPROPS(I,J,5) = ZERO
                     PBEAM_RPROPS(I,J,6) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF
! --- CBEAM_add end --- !

         NAME = 'PBUSH'
         IF (ALLOCATED(PBUSH)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PBUSH(LPBUSH,MPBUSH),STAT=IERR)
            NROWS = LPBUSH
            NCOLS = MPBUSH
            MB_ALLOCATED = RLONG*REAL(LPBUSH)*REAL(MPBUSH)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPBUSH
                  DO J=1,MPBUSH
                     PBUSH(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPBUSH'
         IF (ALLOCATED(RPBUSH)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RPBUSH(LPBUSH,MRPBUSH),STAT=IERR)
            NROWS = LPBUSH
            NCOLS = MRPBUSH
            MB_ALLOCATED = RDOUBLE*REAL(LPBUSH)*REAL(MRPBUSH)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPBUSH
                  DO J=1,MRPBUSH
                     RPBUSH(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PCOMP'
         IF (ALLOCATED(PCOMP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PCOMP(LPCOMP,MPCOMP0+MPCOMP_PLIES*LPCOMP_PLIES),STAT=IERR)
            NROWS = LPCOMP
            NCOLS = MPCOMP0+MPCOMP_PLIES*LPCOMP_PLIES
            MB_ALLOCATED = RLONG*REAL(LPCOMP)*REAL(MPCOMP0+MPCOMP_PLIES*LPCOMP_PLIES)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPCOMP
                  DO J=1,MPCOMP0+MPCOMP_PLIES*LPCOMP_PLIES
                     PCOMP(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPCOMP'
         IF (ALLOCATED(RPCOMP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RPCOMP(LPCOMP,MRPCOMP0+MRPCOMP_PLIES*LPCOMP_PLIES),STAT=IERR)
            NROWS = LPCOMP
            NCOLS = MRPCOMP0+MRPCOMP_PLIES*LPCOMP_PLIES
            MB_ALLOCATED = RDOUBLE*REAL(LPCOMP)*REAL(MRPCOMP0+MRPCOMP_PLIES*LPCOMP_PLIES)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPCOMP
                  DO J=1,MRPCOMP0+MRPCOMP_PLIES*LPCOMP_PLIES
                     RPCOMP(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PELAS'
         IF (ALLOCATED(PELAS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PELAS(LPELAS,MPELAS),STAT=IERR)
            NROWS = LPELAS
            NCOLS = MPELAS
            MB_ALLOCATED = RLONG*REAL(LPELAS)*REAL(MPELAS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPELAS
                  DO J=1,MPELAS
                     PELAS(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPELAS'
         IF (ALLOCATED(RPELAS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RPELAS(LPELAS,MRPELAS),STAT=IERR)
            NROWS = LPELAS
            NCOLS = MRPELAS
            MB_ALLOCATED = RDOUBLE*REAL(LPELAS)*REAL(MRPELAS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPELAS
                  DO J=1,MRPELAS
                     RPELAS(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PROD'
         IF (ALLOCATED(PROD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PROD(LPROD,MPROD),STAT=IERR)
            NROWS = LPROD
            NCOLS = MPROD
            MB_ALLOCATED = RLONG*REAL(LPROD)*REAL(MPROD)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPROD
                  DO J=1,MPROD
                     PROD(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPROD'
         IF (ALLOCATED(RPROD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RPROD(LPROD,MRPROD),STAT=IERR)
            NROWS = LPROD
            NCOLS = MRPROD
            MB_ALLOCATED = RDOUBLE*REAL(LPROD)*REAL(MRPROD)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPROD
                  DO J=1,MRPROD
                     RPROD(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PSHEAR'
         IF (ALLOCATED(PSHEAR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PSHEAR(LPSHEAR,MPSHEAR),STAT=IERR)
            NROWS = LPSHEAR
            NCOLS = MPSHEAR
            MB_ALLOCATED = RLONG*REAL(LPSHEAR)*REAL(MPSHEAR)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPSHEAR
                  DO J=1,MPSHEAR
                     PSHEAR(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPSHEAR'
         IF (ALLOCATED(RPSHEAR)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RPSHEAR(LPSHEAR,MRPSHEAR),STAT=IERR)
            NROWS = LPSHEAR
            NCOLS = MRPSHEAR
            MB_ALLOCATED = RDOUBLE*REAL(LPSHEAR)*REAL(MRPSHEAR)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPSHEAR
                  DO J=1,MRPSHEAR
                     RPSHEAR(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PSHEL'
         IF (ALLOCATED(PSHEL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PSHEL(LPSHEL,MPSHEL),STAT=IERR)
            NROWS = LPSHEL
            NCOLS = MPSHEL
            MB_ALLOCATED = RLONG*REAL(LPSHEL)*REAL(MPSHEL)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPSHEL
                  DO J=1,MPSHEL
                     PSHEL(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPSHEL'
         IF (ALLOCATED(RPSHEL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RPSHEL(LPSHEL,MRPSHEL),STAT=IERR)
            NROWS = LPSHEL
            NCOLS = MRPSHEL
            MB_ALLOCATED = RDOUBLE*REAL(LPSHEL)*REAL(MRPSHEL)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPSHEL
                  DO J=1,MRPSHEL
                     RPSHEL(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PSOLID'
         IF (ALLOCATED(PSOLID)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PSOLID(LPSOLID,MPSOLID),STAT=IERR)
            NROWS = LPSOLID
            NCOLS = MPSOLID
            MB_ALLOCATED = RLONG*REAL(LPSOLID)*REAL(MPSOLID)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPSOLID
                  DO J=1,MPSOLID
                     PSOLID(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PUSER1'
         IF (ALLOCATED(PUSER1)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PUSER1(LPUSER1,MPUSER1),STAT=IERR)
            NROWS = LPUSER1
            NCOLS = MPUSER1
            MB_ALLOCATED = RLONG*REAL(LPUSER1)*REAL(MPUSER1)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPUSER1
                  DO J=1,MPUSER1
                     PUSER1(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPUSER1'
         IF (ALLOCATED(RPUSER1)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RPUSER1(LPUSER1,MRPUSER1),STAT=IERR)
            NROWS = LPUSER1
            NCOLS = MRPUSER1
            MB_ALLOCATED = RDOUBLE*REAL(LPUSER1)*REAL(MRPUSER1)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPUSER1
                  DO J=1,MRPUSER1
                     RPUSER1(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PUSERIN'
         IF (ALLOCATED(PUSERIN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PUSERIN(LPUSERIN,MPUSERIN),STAT=IERR)
            NROWS = LPUSERIN
            NCOLS = MPUSERIN
            MB_ALLOCATED = RLONG*REAL(LPUSERIN)*REAL(MPUSERIN)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPUSERIN
                  DO J=1,MPUSERIN
                     PUSERIN(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'USERIN_MAT_NAMES'
         IF (ALLOCATED(USERIN_MAT_NAMES)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (USERIN_MAT_NAMES(LPUSERIN,MUSERIN_MAT_NAMES),STAT=IERR)
            NROWS = LPUSERIN
            NCOLS = MUSERIN_MAT_NAMES
            MB_ALLOCATED = RLONG*REAL(LPUSERIN)*REAL(MUSERIN_MAT_NAMES)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPUSERIN
                  DO J=1,MUSERIN_MAT_NAMES
                     USERIN_MAT_NAMES(I,J)(1:) = ' '
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'USERIN_ACT_GRDS, USERIN_ACT_COMPS') THEN

         NAME = 'USERIN_ACT_GRIDS'
         IF (ALLOCATED(USERIN_ACT_GRIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (USERIN_ACT_GRIDS(USERIN_NUM_ACT_GRDS),STAT=IERR)
            NROWS = USERIN_NUM_ACT_GRDS
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(USERIN_NUM_ACT_GRDS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,USERIN_NUM_ACT_GRDS
                  USERIN_ACT_GRIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'USERIN_ACT_COMPS'
         IF (ALLOCATED(USERIN_ACT_COMPS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (USERIN_ACT_COMPS(USERIN_NUM_ACT_GRDS),STAT=IERR)
            NROWS = USERIN_NUM_ACT_GRDS
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(USERIN_NUM_ACT_GRDS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,USERIN_NUM_ACT_GRDS
                  USERIN_ACT_COMPS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'VVEC, OFFSETS, PLATE stuff') THEN   ! Allocate arrays for bar v vectors and bar and plate elem offsets

         NAME = 'VVEC'
         IF (ALLOCATED(VVEC)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (VVEC(LVVEC,3),STAT=IERR)
            NROWS = LVVEC
            NCOLS = 3
            MB_ALLOCATED = RDOUBLE*REAL(LVVEC)*THREE/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LVVEC
                  DO J=1,3
                     VVEC(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BAROFF'
         IF (ALLOCATED(BAROFF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (BAROFF(LBAROFF,6),STAT=IERR)
            NROWS = LBAROFF
            NCOLS = 6
            MB_ALLOCATED = RDOUBLE*REAL(LBAROFF)*SIX/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LBAROFF
                  DO J=1,6
                     BAROFF(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BUSHOFF'
         IF (ALLOCATED(BUSHOFF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (BUSHOFF(LBUSHOFF,6),STAT=IERR)
            NROWS = LBUSHOFF
            NCOLS = 6
            MB_ALLOCATED = RDOUBLE*REAL(LBUSHOFF)*SIX/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LBUSHOFF
                  DO J=1,6
                     BUSHOFF(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PLATEOFF'
         IF (ALLOCATED(PLATEOFF)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PLATEOFF(LPLATEOFF),STAT=IERR)
            NROWS = LPLATEOFF
            NCOLS = 1
            MB_ALLOCATED = RDOUBLE*REAL(LPLATEOFF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPLATEOFF
                  PLATEOFF(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PLATETHICK'
         IF (ALLOCATED(PLATETHICK)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PLATETHICK(LPLATETHICK),STAT=IERR)
            NROWS = LPLATETHICK
            NCOLS = 1
            MB_ALLOCATED = RDOUBLE*REAL(LPLATETHICK)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPLATETHICK
                  PLATETHICK(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MATANGLE'
         IF (ALLOCATED(MATANGLE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MATANGLE(LMATANGLE),STAT=IERR)
            NROWS = LMATANGLE
            NCOLS = 1
            MB_ALLOCATED = RDOUBLE*REAL(LMATANGLE)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LMATANGLE
                  MATANGLE(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRID, RGRID') THEN                 ! Allocate arrays for GRID, RGRID

         NAME = 'GRID'
         IF (ALLOCATED(GRID)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GRID(LGRID,MGRID),STAT=IERR)
            NROWS = LGRID
            NCOLS = MGRID
            MB_ALLOCATED = RLONG*REAL(LGRID)*REAL(MGRID)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LGRID
                  DO J=1,MGRID
                     GRID(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RGRID'
         IF (ALLOCATED(RGRID)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RGRID(LGRID,MRGRID),STAT=IERR)
            NROWS = LGRID
            NCOLS = MRGRID
            MB_ALLOCATED = RDOUBLE*REAL(LGRID)*REAL(MRGRID)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LGRID
                  DO J=1,MRGRID
                     RGRID(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CORD, RCORD') THEN                 ! Allocate arrays for CORD, RCORD

         NAME = 'CORD'
         IF (ALLOCATED(CORD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (CORD(LCORD,MCORD),STAT=IERR)
            NROWS = LCORD
            NCOLS = MCORD
            MB_ALLOCATED = RLONG*REAL(LCORD)*REAL(MCORD)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LCORD
                  DO J=1,MCORD
                     CORD(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RCORD'
         IF (ALLOCATED(RCORD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RCORD(LCORD,MRCORD),STAT=IERR)
            NROWS = LCORD
            NCOLS = MRCORD
            MB_ALLOCATED = RDOUBLE*REAL(LCORD)*REAL(MRCORD)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LCORD
                  DO J=1,MRCORD
                     RCORD(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CMASS, PMASS, RPMASS') THEN        ! Allocate arrays for CMASS, PMASS, RPMASS

         NAME = 'CMASS'
         IF (ALLOCATED(CMASS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (CMASS(LCMASS,MCMASS),STAT=IERR)
            NROWS = LCMASS
            NCOLS = MCMASS
            MB_ALLOCATED = RLONG*REAL(LCMASS)*REAL(MCMASS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LCMASS
                  DO J=1,MCMASS
                     CMASS(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PMASS'
         IF (ALLOCATED(PMASS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PMASS(LPMASS,MPMASS),STAT=IERR)
            NROWS = LPMASS
            NCOLS = MPMASS
            MB_ALLOCATED = RLONG*REAL(LPMASS)*REAL(MPMASS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPMASS
                  DO J=1,MPMASS
                     PMASS(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPMASS'
         IF (ALLOCATED(RPMASS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RPMASS(LPMASS,MRPMASS),STAT=IERR)
            NROWS = LPMASS
            NCOLS = MRPMASS
            MB_ALLOCATED = RDOUBLE*REAL(LCMASS)*REAL(MRPMASS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPMASS
                  DO J=1,MRPMASS
                     RPMASS(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CONM2, RCONM2') THEN                 ! Allocate arrays for CONM2, RCONM2

         NAME = 'CONM2'
         IF (ALLOCATED(CONM2)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (CONM2(LCONM2,MCONM2),STAT=IERR)
            NROWS = LCONM2
            NCOLS = MCONM2
            MB_ALLOCATED = RLONG*REAL(LCONM2)*REAL(MCONM2)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LCONM2
                  DO J=1,MCONM2
                     CONM2(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RCONM2'
         IF (ALLOCATED(RCONM2)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RCONM2(LCONM2,MRCONM2),STAT=IERR)
            NROWS = LCONM2
            NCOLS = MRCONM2
            MB_ALLOCATED = RDOUBLE*REAL(LCONM2)*REAL(MRCONM2)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LCONM2
                  DO J=1,MRCONM2
                     RCONM2(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ESORT1') THEN                      ! Allocate arrays for ESORT1

         NAME = 'ESORT1'
         IF (ALLOCATED(ESORT1)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (ESORT1(LELE),STAT=IERR)
            NROWS = LELE
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LELE)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  ESORT1(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ESORT2') THEN                      ! Allocate arrays for ESORT1

         NAME = 'ESORT2'
         IF (ALLOCATED(ESORT2)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (ESORT2(LELE),STAT=IERR)
            NROWS = LELE
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LELE)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  ESORT2(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRID_ID') THEN                     ! Allocate arrays for GRID_ID

         NAME = 'GRID_ID'
         IF (ALLOCATED(GRID_ID)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GRID_ID(LGRID),STAT=IERR)
            NROWS = LGRID
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LGRID)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LGRID
                  GRID_ID(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRID_SEQ, INV_GRID_SEQ') THEN      ! Allocate arrays for grid sequence

         NAME = 'GRID_SEQ'
         IF (ALLOCATED(GRID_SEQ)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GRID_SEQ(LGRID),STAT=IERR)
            NROWS = LGRID
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LGRID)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LGRID
                  GRID_SEQ(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'INV_GRID_SEQ'
         IF (ALLOCATED(INV_GRID_SEQ)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (INV_GRID_SEQ(LGRID),STAT=IERR)
            NROWS = LGRID
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LGRID)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LGRID
                  INV_GRID_SEQ(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF


      ELSE IF (NAME_IN == 'SNORM, RSNORM') THEN            ! Allocate arrays for SNORM, RSNORM

         NAME = 'SNORM'
         IF (ALLOCATED(SNORM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SNORM(LSNORM,MSNORM),STAT=IERR)
            NROWS = LSNORM
            NCOLS = MSNORM
            MB_ALLOCATED = RLONG*REAL(LSNORM)*REAL(MSNORM)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSNORM
                  DO J=1,MSNORM
                     SNORM(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RSNORM'
         IF (ALLOCATED(RSNORM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RSNORM(LSNORM,MRSNORM),STAT=IERR)
            NROWS = LSNORM
            NCOLS = MRSNORM
            MB_ALLOCATED = RDOUBLE*REAL(LSNORM)*REAL(MRSNORM)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSNORM
                  DO J=1,MRSNORM
                     RSNORM(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRID_SNORM') THEN            ! Allocate array for GRID_SNORM

         NAME = 'GRID_SNORM'
         IF (ALLOCATED(GRID_SNORM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GRID_SNORM(NGRID,MGRID_SNORM),STAT=IERR)
            NROWS = NGRID
            NCOLS = MGRID_SNORM
            MB_ALLOCATED = RLONG*REAL(NGRID)*REAL(MGRID_SNORM)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NGRID
                  DO J=1,MGRID_SNORM
                     GRID_SNORM(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'TN') THEN                          ! Allocate array TN

         NAME = 'TN'
         IF (ALLOCATED(TN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TN(3,3,LCORD),STAT=IERR)
            NROWS = 3
            NCOLS = 3*LCORD
            MB_ALLOCATED = RDOUBLE*THREE*THREE*REAL(LCORD)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,3
                  DO J=1,3
                     DO K=1,LCORD
                        TN(I,J,K) = ZERO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GROUT, ELOUT') THEN                ! Allocate arrays for OGROUT, OELOUT, etc

         NAME = 'OGROUT'
         IF (ALLOCATED(OGROUT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OGROUT(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  OGROUT(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'GROUT'
         IF (ALLOCATED(GROUT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GROUT(LGRID,LSUB),STAT=IERR)
            NROWS = LGRID
            NCOLS = LSUB
            MB_ALLOCATED = RLONG*REAL(LGRID)*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LGRID
                  DO J=1,LSUB
                     GROUT(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OELOUT'
         IF (ALLOCATED(OELOUT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OELOUT(LSUB),STAT=IERR)
            NROWS = LSUB
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LSUB
                  OELOUT(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'ELOUT'
         IF (ALLOCATED(ELOUT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (ELOUT(LELE,LSUB),STAT=IERR)
            NROWS = LELE
            NCOLS = LSUB
            MB_ALLOCATED = RLONG*REAL(LELE)*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  DO J=1,LSUB
                     ELOUT(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ELDT') THEN                        ! Allocate array ELDT

         NAME = 'ELDT'
         IF (ALLOCATED(ELDT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (ELDT(LELE),STAT=IERR)
            NROWS = LELE
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LELE)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  ELDT(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SYS_LOAD') THEN                    ! Allocate array SYS_LOAD

         NAME = 'SYS_LOAD'
         IF (ALLOCATED(SYS_LOAD)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SYS_LOAD(NDOFG,LSUB),STAT=IERR)
            NROWS = NDOFG
            NCOLS = LSUB
            MB_ALLOCATED = RDOUBLE*REAL(NDOFG)*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NDOFG
                  DO J=1,LSUB
                     SYS_LOAD(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GTEMP') THEN                       ! Allocate arrays for GTEMP

         NAME = 'GTEMP'
         IF (ALLOCATED(GTEMP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GTEMP(LGRID,NTSUB),STAT=IERR)
            NROWS = LGRID
            NCOLS = NTSUB
            MB_ALLOCATED = RDOUBLE*REAL(LGRID)*REAL(NTSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LGRID
                  DO J=1,NTSUB
                     GTEMP(I,J) = GTEMP_INIT
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CGTEMP') THEN                      ! Allocate arrays for CGTEMP

         NAME = 'CGTEMP'
         IF (ALLOCATED(CGTEMP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (CGTEMP(LGRID,NTSUB),STAT=IERR)
            NROWS = LGRID
            NCOLS = NTSUB
            MB_ALLOCATED = RBYTE*REAL(LEN(CGTEMP))*REAL(LGRID)*REAL(NTSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LGRID
                  DO J=1,NTSUB
                     CGTEMP(I,J) = CGTEMP_ERR
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ETEMP') THEN                       ! Allocate arrays for ETEMP

         NAME = 'ETEMP'
         IF (ALLOCATED(ETEMP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (ETEMP(LELE,NTSUB),STAT=IERR)
            NROWS = LELE
            NCOLS = NTSUB
            MB_ALLOCATED = RDOUBLE*REAL(LELE)*REAL(NTSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  DO J=1,NTSUB
                     ETEMP(I,J) = ETEMP_INIT
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CETEMP') THEN                      ! Allocate arrays for ETEMP

         NAME = 'CETEMP'
         IF (ALLOCATED(CETEMP)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (CETEMP(LELE,NTSUB),STAT=IERR)
            NROWS = LELE
            NCOLS = NTSUB
            MB_ALLOCATED = RBYTE*REAL(LEN(CETEMP))*REAL(LELE)*REAL(NTSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  DO J=1,NTSUB
                     CETEMP(I,J) = CETEMP_ERR
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'TPNT, TDATA') THEN                 ! Allocate arrays for TPNT, TDATA

         NAME = 'TPNT'
         IF (ALLOCATED(TPNT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TPNT(LELE,NTSUB),STAT=IERR)
            NROWS = LELE
            NCOLS = NTSUB
            MB_ALLOCATED = RLONG*REAL(LELE)*REAL(NTSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  DO J=1,NTSUB
                     TPNT(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TDATA'
         IF (ALLOCATED(TDATA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (TDATA(LTDAT),STAT=IERR)
            NROWS = LTDAT
            NCOLS = 1
            MB_ALLOCATED = RDOUBLE*REAL(LTDAT)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LTDAT
                  TDATA(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PPNT, PDATA, PTYPE') THEN       ! Allocate arrays for PPNT, PDATA

         NAME = 'PPNT'
         IF (ALLOCATED(PPNT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PPNT(LELE,LSUB),STAT=IERR)
            NROWS = LELE
            NCOLS = LSUB
            MB_ALLOCATED = RLONG*REAL(LELE)*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  DO J=1,LSUB
                     PPNT(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PDATA'
         IF (ALLOCATED(PDATA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PDATA(LPDAT),STAT=IERR)
            NROWS = LPDAT
            NCOLS = 1
            MB_ALLOCATED = RDOUBLE*REAL(LPDAT)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LPDAT
                  PDATA(I) = ZERO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PTYPE'
         IF (ALLOCATED(PTYPE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PTYPE(LELE),STAT=IERR)
            NROWS = LELE
            NCOLS = 1
            MB_ALLOCATED = REAL(BYTE)*REAL(LPDAT)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LELE
                  PTYPE(I) = ' '
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PLOAD4_3D_DATA') THEN  ! Allocate arrays for PLOAD4_3D_DATA

         NAME = 'PLOAD4_3D_DATA'
         IF (ALLOCATED(PLOAD4_3D_DATA)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PLOAD4_3D_DATA(NPLOAD4_3D,MPLOAD4_3D_DATA),STAT=IERR)
            NROWS = NPLOAD4_3D
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(NPLOAD4_3D)*REAL(MPLOAD4_3D_DATA)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NPLOAD4_3D
                  DO J=1,MPLOAD4_3D_DATA
                     PLOAD4_3D_DATA(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SINGLE ELEMENT ARRAYS') THEN    ! Allocate arrays for elem thermal and pressure loads

         NAME = 'AGRID'
         IF (ALLOCATED(AGRID)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (AGRID(MELGP+1),STAT=IERR)
            NROWS = MELGP+1
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(MELGP+1)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELGP+1
                  AGRID(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BE1'
         IF (ALLOCATED(BE1)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (BE1(3,MELDOF,MAX_STRESS_POINTS+1),STAT=IERR)
            NROWS = 3
            NCOLS = MELDOF*MAX_STRESS_POINTS+1
            MB_ALLOCATED = RDOUBLE*REAL(3)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,3
                  DO J=1,MELDOF
                     DO K=1,MAX_STRESS_POINTS+1
                        BE1(I,J,K) = ZERO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BE2'
         IF (ALLOCATED(BE2)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (BE2(3,MELDOF,MAX_STRESS_POINTS+1),STAT=IERR)
            NROWS = 3
            NCOLS = MELDOF*MAX_STRESS_POINTS+1
            MB_ALLOCATED = RDOUBLE*REAL(3)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,3
                  DO J=1,MELDOF
                     DO K=1,MAX_STRESS_POINTS+1
                        BE2(I,J,K) = ZERO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BE3'
         IF (ALLOCATED(BE3)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (BE3(3,MELDOF,MAX_STRESS_POINTS+1),STAT=IERR)
            NROWS = 3
            NCOLS = MELDOF*MAX_STRESS_POINTS+1
            MB_ALLOCATED = RDOUBLE*REAL(3)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,3
                  DO J=1,MELDOF
                     DO K=1,MAX_STRESS_POINTS+1
                        BE3(I,J,K) = ZERO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BGRID'
         IF (ALLOCATED(BGRID)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (BGRID(MELGP+1),STAT=IERR)
            NROWS = MELGP+1
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(MELGP+1)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELGP+1
                  BGRID(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'DOFPIN'
         IF (ALLOCATED(DOFPIN)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (DOFPIN(MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  DOFPIN(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'DT'
         IF (ALLOCATED(DT)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (DT(MDT,LSUB),STAT=IERR)
            NROWS = MDT
            NCOLS = LSUB
            MB_ALLOCATED = RDOUBLE*REAL(MDT)*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MDT
                  DO J=1,LSUB
                     DT(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KE'
         IF (ALLOCATED(KE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KE(MELDOF,MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = MELDOF
            MB_ALLOCATED = RDOUBLE*REAL(MELDOF)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  DO J=1,MELDOF
                     KE(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KEG'
         IF (ALLOCATED(KEG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KEG(MELDOF,MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = MELDOF
            MB_ALLOCATED = RDOUBLE*REAL(MELDOF)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  DO J=1,MELDOF
                     KEG(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KED'
         IF (ALLOCATED(KED)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KED(MELDOF,MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = MELDOF
            MB_ALLOCATED = RDOUBLE*REAL(MELDOF)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  DO J=1,MELDOF
                     KED(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KEM'
         IF (ALLOCATED(KEM)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (KEM(MELDOF,MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = MELDOF
            MB_ALLOCATED = RDOUBLE*REAL(MELDOF)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  DO J=1,MELDOF
                     KEM(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'ME'
         IF (ALLOCATED(ME)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (ME(MELDOF,MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = MELDOF
            MB_ALLOCATED = RDOUBLE*REAL(MELDOF)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  DO J=1,MELDOF
                     ME(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OFFDIS'
         IF (ALLOCATED(OFFDIS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OFFDIS(MOFFSET,3),STAT=IERR)
            NROWS = MOFFSET
            NCOLS = 3
            MB_ALLOCATED = RDOUBLE*REAL(MOFFSET)*REAL(3)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MOFFSET
                  DO J=1,3
                     OFFDIS(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OFFDIS_O'
         IF (ALLOCATED(OFFDIS_O)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OFFDIS_O(MOFFSET,3),STAT=IERR)
            NROWS = MOFFSET
            NCOLS = 3
            MB_ALLOCATED = RDOUBLE*REAL(MOFFSET)*REAL(3)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MOFFSET
                  DO J=1,3
                     OFFDIS_O(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OFFDIS_B'
         IF (ALLOCATED(OFFDIS_B)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OFFDIS_B(MOFFSET,3),STAT=IERR)
            NROWS = MOFFSET
            NCOLS = 3
            MB_ALLOCATED = RDOUBLE*REAL(MOFFSET)*REAL(3)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MOFFSET
                  DO J=1,3
                     OFFDIS_B(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OFFDIS_G'
         IF (ALLOCATED(OFFDIS_G)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OFFDIS_G(MOFFSET,3),STAT=IERR)
            NROWS = MOFFSET
            NCOLS = 3
            MB_ALLOCATED = RDOUBLE*REAL(MOFFSET)*REAL(3)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MOFFSET
                  DO J=1,3
                     OFFDIS_G(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OFFSET'
         IF (ALLOCATED(OFFSET)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (OFFSET(MOFFSET),STAT=IERR)
            NROWS = MOFFSET
            NCOLS = 1
            MB_ALLOCATED = RDOUBLE*REAL(MOFFSET)*REAL(3)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MOFFSET
                  OFFSET(I)(1:) = 'N'
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PEB'
         IF (ALLOCATED(PEB)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PEB(MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  PEB(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PEG'
         IF (ALLOCATED(PEG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PEG(MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  PEG(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PEL'
         IF (ALLOCATED(PEL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PEL(MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  PEL(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PPE'
         IF (ALLOCATED(PPE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PPE(MELDOF,LSUB),STAT=IERR)
            NROWS = MELDOF
            NCOLS = LSUB
            MB_ALLOCATED = RDOUBLE*REAL(MELDOF)*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  DO J=1,LSUB
                     PPE(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PRESS'
         IF (ALLOCATED(PRESS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PRESS(MPRESS,LSUB),STAT=IERR)
            NROWS = MPRESS
            NCOLS = LSUB
            MB_ALLOCATED = RDOUBLE*REAL(MPRESS)*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MPRESS
                  DO J=1,LSUB
                     PRESS(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PTE'
         IF (ALLOCATED(PTE)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (PTE(MELDOF,LSUB),STAT=IERR)
            NROWS = MELDOF
            NCOLS = LSUB
            MB_ALLOCATED = RDOUBLE*REAL(MELDOF)*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  DO J=1,LSUB
                     PTE(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SE1'
         IF (ALLOCATED(SE1)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SE1(3,MELDOF,MAX_STRESS_POINTS+1),STAT=IERR)
            NROWS = 3
            NCOLS = MELDOF*MAX_STRESS_POINTS+1
            MB_ALLOCATED = RDOUBLE*REAL(3)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,3
                  DO J=1,MELDOF
                     DO K=1,MAX_STRESS_POINTS+1
                        SE1(I,J,K) = ZERO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SE2'
         IF (ALLOCATED(SE2)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SE2(3,MELDOF,MAX_STRESS_POINTS+1),STAT=IERR)
            NROWS = 3
            NCOLS = MELDOF*MAX_STRESS_POINTS+1
            MB_ALLOCATED = RDOUBLE*REAL(3)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,3
                  DO J=1,MELDOF
                     DO K=1,MAX_STRESS_POINTS+1
                        SE2(I,J,K) = ZERO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SE3'
         IF (ALLOCATED(SE3)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (SE3(3,MELDOF,MAX_STRESS_POINTS+1),STAT=IERR)
            NROWS = 3
            NCOLS = MELDOF*MAX_STRESS_POINTS+1
            MB_ALLOCATED = RDOUBLE*REAL(3)*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,3
                  DO J=1,MELDOF
                     DO K=1,MAX_STRESS_POINTS+1
                        SE3(I,J,K) = ZERO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'STE1'
         IF (ALLOCATED(STE1)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (STE1(3,LSUB,MAX_STRESS_POINTS+1),STAT=IERR)
            NROWS = 3
            NCOLS = LSUB*MAX_STRESS_POINTS+1
            MB_ALLOCATED = RDOUBLE*THREE*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,3
                  DO J=1,LSUB
                     DO K=1,MAX_STRESS_POINTS+1
                        STE1(I,J,K)   = ZERO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'STE2'
         IF (ALLOCATED(STE2)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (STE2(3,LSUB,MAX_STRESS_POINTS+1),STAT=IERR)
            NROWS = 3
            NCOLS = LSUB*MAX_STRESS_POINTS+1
            MB_ALLOCATED = RDOUBLE*THREE*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,3
                  DO J=1,LSUB
                     DO K=1,MAX_STRESS_POINTS+1
                        STE2(I,J,K)   = ZERO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'STE3'
         IF (ALLOCATED(STE3)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (STE3(3,LSUB,MAX_STRESS_POINTS+1),STAT=IERR)
            NROWS = 3
            NCOLS = LSUB*MAX_STRESS_POINTS+1
            MB_ALLOCATED = RDOUBLE*THREE*REAL(LSUB)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,3
                  DO J=1,LSUB
                     DO K=1,MAX_STRESS_POINTS+1
                        STE3(I,J,K)   = ZERO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'UEB'
         IF (ALLOCATED(UEB)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (UEB(MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  UEB(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'UEG'
         IF (ALLOCATED(UEG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (UEG(MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  UEG(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'UEL'
         IF (ALLOCATED(UEL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (UEL(MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  UEL(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'UGG'
         IF (ALLOCATED(UGG)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (UGG(MELDOF),STAT=IERR)
            NROWS = MELDOF
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(MELDOF)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELDOF
                  UGG(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'XEB'
         IF (ALLOCATED(XEB)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (XEB(MELGP+1,3),STAT=IERR)
            NROWS = MELGP+1
            NCOLS = 3
            MB_ALLOCATED = RDOUBLE*REAL(MELGP+1)*REAL(3)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELGP+1
                  DO J=1,3
                     XEB(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'XEL'
         IF (ALLOCATED(XEL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (XEL(MELGP,3),STAT=IERR)
            NROWS = MELGP
            NCOLS = 3
            MB_ALLOCATED = RDOUBLE*REAL(MELGP)*REAL(3)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MELGP
                  DO J=1,3
                     XEL(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'XGL'
         IF (ALLOCATED(XGL)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (XGL(MAX_GAUSS_POINTS*MAX_GAUSS_POINTS,3),STAT=IERR)
            NROWS = MAX_GAUSS_POINTS*MAX_GAUSS_POINTS
            NCOLS = 3
            MB_ALLOCATED = RDOUBLE*REAL(MELGP)*REAL(3)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,MAX_GAUSS_POINTS
                  DO J=1,2
                     XGL(I,J) = ZERO
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'RIGID_ELEM_IDS') THEN              ! Allocate arrays for rigid element ID's

         NAME = 'RIGID_ELEM_IDS'
         IF (ALLOCATED(RIGID_ELEM_IDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (RIGID_ELEM_IDS(LRIGEL),STAT=IERR)
            NROWS = LRIGEL
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LRIGEL)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LRIGEL
                  RIGID_ELEM_IDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MPC_IND_GRIDS') THEN               ! Allocate array for MPC_IND_GRIDS

         NAME = 'MPC_IND_GRIDS'
         IF (ALLOCATED(MPC_IND_GRIDS)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (MPC_IND_GRIDS(LIND_GRDS_MPCS),STAT=IERR)
            NROWS = LIND_GRDS_MPCS
            NCOLS = 1
            MB_ALLOCATED = RLONG*REAL(LIND_GRDS_MPCS)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,LIND_GRDS_MPCS
                  MPC_IND_GRIDS(I) = 0
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRID_ELEM_CONN_ARRAY') THEN        ! Allocate array for GRID_ELEM_CONN_ARRAY

         NAME = 'GRID_ELEM_CONN_ARRAY'
         IF (ALLOCATED(GRID_ELEM_CONN_ARRAY)) THEN
            WRITE(ERR,990) SUBR_NAME, NAME
            WRITE(F06,990) SUBR_NAME, NAME
            FATAL_ERR = FATAL_ERR + 1
            JERR = JERR + 1
         ELSE
            ALLOCATE (GRID_ELEM_CONN_ARRAY(NGRID,MAX_ELEM_DEGREE+2),STAT=IERR)
            NROWS = NGRID
            NCOLS = MAX_ELEM_DEGREE + 2
            MB_ALLOCATED = RLONG*REAL(NGRID)*REAL(MAX_ELEM_DEGREE+2)/ONEPP6
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, MB_ALLOCATED, 'ALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
               DO I=1,NGRID
                  DO J=1,MAX_ELEM_DEGREE+2
                     GRID_ELEM_CONN_ARRAY(I,J) = 0
                  ENDDO
               ENDDO
            ELSE
               WRITE(ERR,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               WRITE(F06,991) MB_ALLOCATED, NAME,SUBR_NAME, IERR
               FATAL_ERR = FATAL_ERR + 1
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,915) SUBR_NAME, 'ALLOCATED', NAME_IN
         WRITE(F06,915) SUBR_NAME, 'ALLOCATED', NAME_IN
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         WRITE(ERR,1699) TRIM(SUBR_NAME), CALLING_SUBR
         WRITE(F06,1699) TRIM(SUBR_NAME), CALLING_SUBR
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  990 FORMAT(' *ERROR   990: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT ALLOCATE MEMORY TO ARRAY ',A,'. IT IS ALREADY ALLOCATED')

  991 FORMAT(' *ERROR   991: CANNOT ALLOCATE ',F10.3,' MB OF MEMORY TO ARRAY ',A,' IN SUBROUTINE ',A                               &
                    ,/,14X,' ALLOCATION STAT = ',I8)

 1699 FORMAT('               THE SUBR IN WHICH THESE ERRORS WERE FOUND (',A,') WAS CALLED BY SUBR ',A)

! **********************************************************************************************************************************

      END SUBROUTINE ALLOCATE_MODEL_STUF
