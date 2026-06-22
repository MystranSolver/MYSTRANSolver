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

      SUBROUTINE DEALLOCATE_MODEL_STUF ( NAME_IN )

! DEallocates arrays defined in module MODEL_STUF (see comments there for definition of these matrices)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, TOT_MB_MEM_ALLOC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      USE MODEL_STUF, ONLY            :  AGRID, BE1, BE2, BE3, BGRID, DOFPIN, DT, KE, KEG, KED, KEM, ME,                           &
                                         OFFDIS, OFFDIS_B, OFFDIS_G, OFFDIS_O, OFFSET, PEB, PEG, PEL, PPE, PRESS, PTE,             &
                                         SE1, SE2, SE3, STE1, STE2, STE3, UEB, UEG, UEL, UGG, XEB, XEL, XGL
      USE MODEL_STUF, ONLY            :  CMASS, CONM2, PMASS, RPMASS, RCONM2
      USE MODEL_STUF, ONLY            :  CORD, RCORD, TN
      USE MODEL_STUF, ONLY            :  SEQ1, SEQ2
      USE MODEL_STUF, ONLY            :  BAROFF, BUSHOFF, EDAT, EOFF, EPNT, ESORT1, ESORT2, ETYPE, VVEC
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
      USE MODEL_STUF, ONLY            :  PBEAM_NSTATIONS, PBEAM_XL, PBEAM_RPROPS ! 
! --- cbeam_add emd --- !
      USE MODEL_STUF, ONLY            :  MPC_SIDS, MPCSIDS, MPCADD_SIDS
      USE MODEL_STUF, ONLY            :  SPC_SIDS, SPC1_SIDS, SPCSIDS, SPCADD_SIDS
      USE MODEL_STUF, ONLY            :  ALL_SETS_ARRAY, ONE_SET_ARRAY, SETS_IDS, SC_ACCE, SC_DISP, SC_ELFN, SC_ELFE, SC_GPFO,     &
                                         SC_MPCF, SC_OLOA, SC_SPCF, SC_STRE, SC_STRN, LOAD_SIDS, LOAD_FACS
      USE MODEL_STUF, ONLY            :  ELDT, ELOUT, GROUT, OELOUT, OGROUT, LABEL, SCNUM, STITLE, SUBLOD, TITLE
      USE MODEL_STUF, ONLY            :  CC_EIGR_SID_SUB, CC_STATSUB_SUB, EIG_PARAMS, IS_BUCKLING_SUBCASE, IS_MODES_SUBCASE,    &
                                         NUM_EIGENS_SUB
      USE MODEL_STUF, ONLY            :  SYS_LOAD
      USE MODEL_STUF, ONLY            :  CETEMP, CETEMP_ERR, CGTEMP, CGTEMP_ERR, ETEMP, ETEMP_INIT, GTEMP, GTEMP_INIT, TDATA, TPNT
      USE MODEL_STUF, ONLY            :  RIGID_ELEM_IDS
      USE MODEL_STUF, ONLY            :  MATANGLE, PLATEOFF, PLATETHICK
      USE MODEL_STUF, ONLY            :  GRID_ELEM_CONN_ARRAY
      USE MODEL_STUF, ONLY            :  SPCSETS, MPCSETS

      USE DEALLOCATE_MODEL_STUF_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DEALLOCATE_MODEL_STUF'
      CHARACTER(LEN=*), INTENT(IN)    :: NAME_IN           ! Name of group of arrays to allocate
      CHARACTER(31*BYTE)              :: NAME              ! Specific array name used for output error message

      INTEGER(LONG)                   :: IERR              ! STAT from DEALLOCATE
      INTEGER(LONG)                   :: JERR              ! Local error indicator


      REAL(DOUBLE)                    :: CUR_MB_ALLOCATED  ! MB of memory that is currently allocated to ARRAY_NAME when subr
!                                                            ALLOCATED_MEMORY is called (before entering MB_ALLOCATED into array
!                                                            ALLOCATED_ARRAY_MEM



! **********************************************************************************************************************************
      JERR = 0

      IF (NAME_IN == 'SETS ARRAYS') THEN                   ! Deallocate arrays for SET's

         NAME = 'SETS_IDS'
         IF (ALLOCATED(SETS_IDS)) THEN
            DEALLOCATE (SETS_IDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'ALL_SETS_ARRAY'
         IF (ALLOCATED(ALL_SETS_ARRAY)) THEN
            DEALLOCATE (ALL_SETS_ARRAY,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'ONE_SET_ARRAY'
         IF (ALLOCATED(ONE_SET_ARRAY)) THEN
            DEALLOCATE (ONE_SET_ARRAY,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'TITLES') THEN                   ! Deallocate arrays for output title, subtitle, label

         NAME = 'TITLE'
         IF (ALLOCATED(TITLE)) THEN
            DEALLOCATE (TITLE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'STITLE'
         IF (ALLOCATED(STITLE)) THEN
            DEALLOCATE (STITLE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'LABEL'
         IF (ALLOCATED(LABEL)) THEN
            DEALLOCATE (LABEL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SC_xxxx') THEN                  ! Deallocate arrays for SC_DIDP, etc

         NAME = 'SC_ACCE'
         IF (ALLOCATED(SC_ACCE)) THEN
            DEALLOCATE (SC_ACCE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_DISP'
         IF (ALLOCATED(SC_DISP)) THEN
            DEALLOCATE (SC_DISP,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_ELFE'
         IF (ALLOCATED(SC_ELFE)) THEN
            DEALLOCATE (SC_ELFE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_ELFN'
         IF (ALLOCATED(SC_ELFN)) THEN
            DEALLOCATE (SC_ELFN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_GPFO'
         IF (ALLOCATED(SC_GPFO)) THEN
            DEALLOCATE (SC_GPFO,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_MPCF'
         IF (ALLOCATED(SC_MPCF)) THEN
            DEALLOCATE (SC_MPCF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_OLOA'
         IF (ALLOCATED(SC_OLOA)) THEN
            DEALLOCATE (SC_OLOA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_SPCF'
         IF (ALLOCATED(SC_SPCF)) THEN
            DEALLOCATE (SC_SPCF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_STRE'
         IF (ALLOCATED(SC_STRE)) THEN
            DEALLOCATE (SC_STRE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SC_STRN'
         IF (ALLOCATED(SC_STRN)) THEN
            DEALLOCATE (SC_STRN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SCNUM') THEN                    ! Deallocate array SC_NUM

         NAME = 'SCNUM'
         IF (ALLOCATED(SCNUM)) THEN
            DEALLOCATE (SCNUM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         ! NOTE: CC_EIGR_SID_SUB, CC_STATSUB_SUB, IS_MODES_SUBCASE, IS_BUCKLING_SUBCASE, NUM_EIGENS_SUB, and EIG_PARAMS
         ! are deliberately *not* deallocated alongside SCNUM. LINK9 deallocates+reallocates SCNUM at startup, and we
         ! need the per-subcase eigen/buckling state populated by LK1/LK4 to survive that cycle so LINK9 can iterate
         ! modes per their owning subcase. The matching ALLOCATE_MODEL_STUF('SCNUM') skips them when they are already
         ! allocated, so they remain populated.

      ELSE IF (NAME_IN == 'SUBLOD') THEN                   ! Deallocate array SUBLOD

         NAME = 'SUBLOD'
         IF (ALLOCATED(SUBLOD)) THEN
            DEALLOCATE (SUBLOD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SPC_MPC_SET') THEN                   ! Deallocate array SUBLOD

         NAME = 'SPCSETS'
         IF (ALLOCATED(SPCSETS)) THEN
            DEALLOCATE (SPCSETS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MPCSETS'
         IF (ALLOCATED(MPCSETS)) THEN
            DEALLOCATE (MPCSETS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SEQ1,2') THEN                   ! Deallocate arrays for SEQ1,2

         NAME = 'SEQ1'
         IF (ALLOCATED(SEQ1)) THEN
            DEALLOCATE (SEQ1,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SEQ2'
         IF (ALLOCATED(SEQ2)) THEN
            DEALLOCATE (SEQ2,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'FORMOM_SIDS') THEN              ! Deallocate arrays for force/moment set ID's

         NAME = 'FORMOM_SIDS'
         IF (ALLOCATED(FORMOM_SIDS)) THEN
            DEALLOCATE (FORMOM_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PRESS_SIDS') THEN               ! Deallocate arrays for pressure load set ID's

         NAME = 'PRESS_SIDS'
         IF (ALLOCATED(PRESS_SIDS)) THEN
            DEALLOCATE (PRESS_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRAV_SIDS') THEN                ! Deallocate arrays for grav load set ID's

         NAME = 'GRAV_SIDS'
         IF (ALLOCATED(GRAV_SIDS)) THEN
            DEALLOCATE (GRAV_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'LOAD_SIDS, LOAD_FACS') THEN     ! Deallocate arrays for load set ID's and factors

         NAME = 'LOAD_SIDS'
         IF (ALLOCATED(LOAD_SIDS)) THEN
            DEALLOCATE (LOAD_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'LOAD_FACS'
         IF (ALLOCATED(LOAD_FACS)) THEN
            DEALLOCATE (LOAD_FACS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MPC_SIDS') THEN                 ! Deallocate arrays for MPC set ID's

         NAME = 'MPC_SIDS'
         IF (ALLOCATED(MPC_SIDS)) THEN
            DEALLOCATE (MPC_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MPCSIDS') THEN                  ! Deallocate arrays for MPC set ID's used in 1 execution

         NAME = 'MPCSIDS'
         IF (ALLOCATED(MPCSIDS)) THEN
            DEALLOCATE (MPCSIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MPCADD_SIDS') THEN              ! Deallocate arrays for MPCADD set ID's

         NAME = 'MPCADD_SIDS'
         IF (ALLOCATED(MPCADD_SIDS)) THEN
            DEALLOCATE (MPCADD_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'RFORCE_SIDS') THEN              ! Deallocate arrays for grav load set ID's

         NAME = 'RFORCE_SIDS'
         IF (ALLOCATED(RFORCE_SIDS)) THEN
            DEALLOCATE (RFORCE_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SLOAD_SIDS') THEN               ! Deallocate arrays for grav load set ID's

         NAME = 'SLOAD_SIDS'
         IF (ALLOCATED(SLOAD_SIDS)) THEN
            DEALLOCATE (SLOAD_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SPC_SIDS, SPC1_SIDS') THEN      ! Deallocate arrays for SPC, SPC1 set ID's

         NAME = 'SPC_SIDS'
         IF (ALLOCATED(SPC_SIDS)) THEN
            DEALLOCATE (SPC_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SPC1_SIDS'
         IF (ALLOCATED(SPC1_SIDS)) THEN
            DEALLOCATE (SPC1_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SPCSIDS') THEN                  ! Deallocate arrays for SPC set ID's used in 1 execution

         NAME = 'SPCSIDS'
         IF (ALLOCATED(SPCSIDS)) THEN
            DEALLOCATE (SPCSIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SPCADD_SIDS') THEN              ! Deallocate arrays for SPCADD set ID's

         NAME = 'SPCADD_SIDS'
         IF (ALLOCATED(SPCADD_SIDS)) THEN
            DEALLOCATE (SPCADD_SIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ETYPE, EDAT, EPNT') THEN        ! Deallocate arrays for ETYPE, EDAT, EPNT

         NAME = 'ETYPE'
         IF (ALLOCATED(ETYPE)) THEN
            DEALLOCATE (ETYPE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'EDAT'
         IF (ALLOCATED(EDAT)) THEN
            DEALLOCATE (EDAT,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'EPNT'
         IF (ALLOCATED(EPNT)) THEN
            DEALLOCATE (EPNT,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'EOFF') THEN                     ! Deallocate array for EOFF

         NAME = 'EOFF'
         IF (ALLOCATED(EOFF)) THEN
            DEALLOCATE (EOFF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ELEM PROPERTIES AND MATERIALS') THEN ! Deallocate arrays for element property and materials

         NAME = 'MATL'
         IF (ALLOCATED(MATL)) THEN
            DEALLOCATE (MATL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RMATL'
         IF (ALLOCATED(RMATL)) THEN
            DEALLOCATE (RMATL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PBAR'
         IF (ALLOCATED(PBAR)) THEN
            DEALLOCATE (PBAR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPBAR'
         IF (ALLOCATED(RPBAR)) THEN
            DEALLOCATE (RPBAR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PBEAM'
         IF (ALLOCATED(PBEAM)) THEN
            DEALLOCATE (PBEAM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPBEAM'
         IF (ALLOCATED(RPBEAM)) THEN
            DEALLOCATE (RPBEAM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

! --- cbeam_add begin ---!
         NAME = 'PBEAM_NSTATIONS'
         IF (ALLOCATED(PBEAM_NSTATIONS)) THEN
            DEALLOCATE (PBEAM_NSTATIONS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PBEAM_XL'
         IF (ALLOCATED(PBEAM_XL)) THEN
            DEALLOCATE (PBEAM_XL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PBEAM_RPROPS'
         IF (ALLOCATED(PBEAM_RPROPS)) THEN
            DEALLOCATE (PBEAM_RPROPS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF
! --- cbeam_add end ---!
         NAME = 'PBUSH'
         IF (ALLOCATED(PBUSH)) THEN
            DEALLOCATE (PBUSH,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPBUSH'
         IF (ALLOCATED(RPBUSH)) THEN
            DEALLOCATE (RPBUSH,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PCOMP'
         IF (ALLOCATED(PCOMP)) THEN
            DEALLOCATE (PCOMP,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPCOMP'
         IF (ALLOCATED(RPCOMP)) THEN
            DEALLOCATE (RPCOMP,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PELAS'
         IF (ALLOCATED(PELAS)) THEN
            DEALLOCATE (PELAS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPELAS'
         IF (ALLOCATED(RPELAS)) THEN
            DEALLOCATE (RPELAS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PROD'
         IF (ALLOCATED(PROD)) THEN
            DEALLOCATE (PROD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPROD'
         IF (ALLOCATED(RPROD)) THEN
            DEALLOCATE (RPROD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PSHEAR'
         IF (ALLOCATED(PSHEAR)) THEN
            DEALLOCATE (PSHEAR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPSHEAR'
         IF (ALLOCATED(RPSHEAR)) THEN
            DEALLOCATE (RPSHEAR,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PSHEL'
         IF (ALLOCATED(PSHEL)) THEN
            DEALLOCATE (PSHEL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPSHEL'
         IF (ALLOCATED(RPSHEL)) THEN
            DEALLOCATE (RPSHEL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PSOLID'
         IF (ALLOCATED(PSOLID)) THEN
            DEALLOCATE (PSOLID,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PUSER1'
         IF (ALLOCATED(PUSER1)) THEN
            DEALLOCATE (PUSER1,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPUSER1'
         IF (ALLOCATED(RPUSER1)) THEN
            DEALLOCATE (RPUSER1,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PUSERIN'
         IF (ALLOCATED(PUSERIN)) THEN
            DEALLOCATE (PUSERIN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'USERIN_MAT_NAMES'
         IF (ALLOCATED(USERIN_MAT_NAMES)) THEN
            DEALLOCATE (USERIN_MAT_NAMES,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'USERIN_ACT_GRDS, USERIN_ACT_COMPS') THEN

         NAME = 'USERIN_ACT_COMPS'
         IF (ALLOCATED(USERIN_ACT_COMPS)) THEN
            DEALLOCATE (USERIN_ACT_COMPS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'USERIN_ACT_GRIDS'
         IF (ALLOCATED(USERIN_ACT_GRIDS)) THEN
            DEALLOCATE (USERIN_ACT_GRIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'VVEC, OFFSETS, PLATE stuff') THEN! Deallocate arrays for bar v vectors and offsets

         NAME = 'VVEC'
         IF (ALLOCATED(VVEC)) THEN
            DEALLOCATE (VVEC,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BAROFF'
         IF (ALLOCATED(BAROFF)) THEN
            DEALLOCATE (BAROFF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BUSHOFF'
         IF (ALLOCATED(BUSHOFF)) THEN
            DEALLOCATE (BUSHOFF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PLATEOFF'
         IF (ALLOCATED(PLATEOFF)) THEN
            DEALLOCATE (PLATEOFF,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PLATETHICK'
         IF (ALLOCATED(PLATETHICK)) THEN
            DEALLOCATE (PLATETHICK,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'MATANGLE'
         IF (ALLOCATED(MATANGLE)) THEN
            DEALLOCATE (MATANGLE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRID, RGRID') THEN              ! Deallocate arrays for GRID, RGRID

         NAME = 'GRID'
         IF (ALLOCATED(GRID)) THEN
            DEALLOCATE (GRID,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RGRID'
         IF (ALLOCATED(RGRID)) THEN
            DEALLOCATE (RGRID,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CORD, RCORD') THEN              ! Deallocate arrays for coord systems

         NAME = 'CORD'
         IF (ALLOCATED(CORD)) THEN
            DEALLOCATE (CORD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RCORD'
         IF (ALLOCATED(RCORD)) THEN
            DEALLOCATE (RCORD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CMASS, PMASS, RPMASS') THEN     ! Deallocate arrays for CMASS's

         NAME = 'CMASS'
         IF (ALLOCATED(CMASS)) THEN
            DEALLOCATE (CMASS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PMASS'
         IF (ALLOCATED(PMASS)) THEN
            DEALLOCATE (PMASS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RPMASS'
         IF (ALLOCATED(RPMASS)) THEN
            DEALLOCATE (RPMASS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CONM2, RCONM2') THEN            ! Deallocate arrays for CONM2's

         NAME = 'CONM2'
         IF (ALLOCATED(CONM2)) THEN
            DEALLOCATE (CONM2,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RCONM2'
         IF (ALLOCATED(RCONM2)) THEN
            DEALLOCATE (RCONM2,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ESORT1') THEN                   ! Deallocate arrays for ESORT1

         NAME = 'ESORT1'
         IF (ALLOCATED(ESORT1)) THEN
            DEALLOCATE (ESORT1,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ESORT2') THEN                   ! Deallocate arrays for ESORT2

         NAME = 'ESORT2'
         IF (ALLOCATED(ESORT2)) THEN
            DEALLOCATE (ESORT2,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRID_ID') THEN                  ! Deallocate arrays for GRID_ID

         NAME = 'GRID_ID'
         IF (ALLOCATED(GRID_ID)) THEN
            DEALLOCATE (GRID_ID,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRID_SEQ, INV_GRID_SEQ') THEN   ! Deallocate arrays for grid sequence

         NAME = 'GRID_SEQ'
         IF (ALLOCATED(GRID_SEQ)) THEN
            DEALLOCATE (GRID_SEQ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'INV_GRID_SEQ'
         IF (ALLOCATED(INV_GRID_SEQ)) THEN
            DEALLOCATE (INV_GRID_SEQ,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SNORM, RSNORM') THEN            ! Deallocate arrays for SNORM, RSNORM

         NAME = 'SNORM'
         IF (ALLOCATED(SNORM)) THEN
            DEALLOCATE (SNORM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'RSNORM'
         IF (ALLOCATED(RSNORM)) THEN
            DEALLOCATE (RSNORM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRID_SNORM') THEN            ! Deallocate arrays for GRID_SNORM

         NAME = 'GRID_SNORM'
         IF (ALLOCATED(GRID_SNORM)) THEN
            DEALLOCATE (GRID_SNORM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'TN') THEN                       ! Deallocate array TN

         NAME = 'TN'
         IF (ALLOCATED(TN)) THEN
            DEALLOCATE (TN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GROUT, ELOUT') THEN             ! Deallocate arrays for OGROUT, OELOUT, etc

         NAME = 'OGROUT'
         IF (ALLOCATED(OGROUT)) THEN
            DEALLOCATE (OGROUT,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'GROUT'
         IF (ALLOCATED(GROUT)) THEN
            DEALLOCATE (GROUT,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OELOUT'
         IF (ALLOCATED(OELOUT)) THEN
            DEALLOCATE (OELOUT,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'ELOUT'
         IF (ALLOCATED(ELOUT)) THEN
            DEALLOCATE (ELOUT,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ELDT') THEN                     ! Deallocate array ELDT

         NAME = 'ELDT'
         IF (ALLOCATED(ELDT)) THEN
            DEALLOCATE (ELDT,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SYS_LOAD') THEN                 ! Deallocate array SYS_LOAD

         NAME = 'SYS_LOAD'
         IF (ALLOCATED(SYS_LOAD)) THEN
            DEALLOCATE (SYS_LOAD,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GTEMP') THEN                    ! Deallocate arrays for GTEMP

         NAME = 'GTEMP'
         IF (ALLOCATED(GTEMP)) THEN
            DEALLOCATE (GTEMP,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CGTEMP') THEN                   ! Deallocate arrays for CGTEMP

         NAME = 'CGTEMP'
         IF (ALLOCATED(CGTEMP)) THEN
            DEALLOCATE (CGTEMP,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'ETEMP') THEN                    ! Deallocate arrays for ETEMP

         NAME = 'ETEMP'
         IF (ALLOCATED(ETEMP)) THEN
            DEALLOCATE (ETEMP,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'CETEMP') THEN                   ! Deallocate arrays for CETEMP

         NAME = 'CETEMP'
         IF (ALLOCATED(CETEMP)) THEN
            DEALLOCATE (CETEMP,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'TPNT, TDATA') THEN              ! Deallocate arrays for TPNT, TDATA

         NAME = 'TPNT'
         IF (ALLOCATED(TPNT)) THEN
            DEALLOCATE (TPNT,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'TDATA'
         IF (ALLOCATED(TDATA)) THEN
            DEALLOCATE (TDATA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PPNT, PDATA, PTYPE') THEN       ! Deallocate arrays for PPNT, PDATA

         NAME = 'PPNT'
         IF (ALLOCATED(PPNT)) THEN
            DEALLOCATE (PPNT,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PDATA'
         IF (ALLOCATED(PDATA)) THEN
            DEALLOCATE (PDATA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PTYPE'
         IF (ALLOCATED(PTYPE)) THEN
            DEALLOCATE (PTYPE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'PLOAD4_3D_DATA') THEN  ! Deallocate arrays for PLOAD4_3D_DATA

         NAME = 'PLOAD4_3D_DATA'
         IF (ALLOCATED(PLOAD4_3D_DATA)) THEN
            DEALLOCATE (PLOAD4_3D_DATA,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'SINGLE ELEMENT ARRAYS') THEN! Deallocate arrays for element subcase arrays

         NAME = 'AGRID'
         IF (ALLOCATED(AGRID)) THEN
            DEALLOCATE (AGRID,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BE1'
         IF (ALLOCATED(BE1)) THEN
            DEALLOCATE (BE1,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BE2'
         IF (ALLOCATED(BE2)) THEN
            DEALLOCATE (BE2,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BE3'
         IF (ALLOCATED(BE3)) THEN
            DEALLOCATE (BE3,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'BGRID'
         IF (ALLOCATED(BGRID)) THEN
            DEALLOCATE (BGRID,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'DOFPIN'
         IF (ALLOCATED(DOFPIN)) THEN
            DEALLOCATE (DOFPIN,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'DT'
         IF (ALLOCATED(DT)) THEN
            DEALLOCATE (DT,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KE'
         IF (ALLOCATED(KE)) THEN
            DEALLOCATE (KE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KEG'
         IF (ALLOCATED(KEG)) THEN
            DEALLOCATE (KEG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KED'
         IF (ALLOCATED(KED)) THEN
            DEALLOCATE (KED,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'KEM'
         IF (ALLOCATED(KEM)) THEN
            DEALLOCATE (KEM,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'ME'
         IF (ALLOCATED(ME)) THEN
            DEALLOCATE (ME,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OFFDIS'
         IF (ALLOCATED(OFFDIS)) THEN
            DEALLOCATE (OFFDIS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OFFDIS_O'
         IF (ALLOCATED(OFFDIS_O)) THEN
            DEALLOCATE (OFFDIS_O,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OFFDIS_B'
         IF (ALLOCATED(OFFDIS_B)) THEN
            DEALLOCATE (OFFDIS_B,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OFFDIS_G'
         IF (ALLOCATED(OFFDIS_G)) THEN
            DEALLOCATE (OFFDIS_G,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'OFFSET'
         IF (ALLOCATED(OFFSET)) THEN
            DEALLOCATE (OFFSET,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PEB'
         IF (ALLOCATED(PEB)) THEN
            DEALLOCATE (PEB,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PEG'
         IF (ALLOCATED(PEG)) THEN
            DEALLOCATE (PEG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PEL'
         IF (ALLOCATED(PEL)) THEN
            DEALLOCATE (PEL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PPE'
         IF (ALLOCATED(PPE)) THEN
            DEALLOCATE (PPE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PRESS'
         IF (ALLOCATED(PRESS)) THEN
            DEALLOCATE (PRESS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'PTE'
         IF (ALLOCATED(PTE)) THEN
            DEALLOCATE (PTE,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SE1'
         IF (ALLOCATED(SE1)) THEN
            DEALLOCATE (SE1,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SE2'
         IF (ALLOCATED(SE2)) THEN
            DEALLOCATE (SE2,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'SE3'
         IF (ALLOCATED(SE3)) THEN
            DEALLOCATE (SE3,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'STE1'
         IF (ALLOCATED(STE1)) THEN
            DEALLOCATE (STE1,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'STE2'
         IF (ALLOCATED(STE2)) THEN
            DEALLOCATE (STE2,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'STE3'
         IF (ALLOCATED(STE3)) THEN
            DEALLOCATE (STE3,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'UEB'
         IF (ALLOCATED(UEB)) THEN
            DEALLOCATE (UEB,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'UEG'
         IF (ALLOCATED(UEG)) THEN
            DEALLOCATE (UEG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'UEL'
         IF (ALLOCATED(UEL)) THEN
            DEALLOCATE (UEL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'UGG'
         IF (ALLOCATED(UGG)) THEN
            DEALLOCATE (UGG,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'XEB'
         IF (ALLOCATED(XEB)) THEN
            DEALLOCATE (XEB,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'XEL'
         IF (ALLOCATED(XEL)) THEN
            DEALLOCATE (XEL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

         NAME = 'XGL'
         IF (ALLOCATED(XGL)) THEN
            DEALLOCATE (XGL,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'RIGID_ELEM_IDS') THEN           ! Deallocate array for rigid element ID's

         NAME = 'RIGID_ELEM_IDS'
         IF (ALLOCATED(RIGID_ELEM_IDS)) THEN
            DEALLOCATE (RIGID_ELEM_IDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'MPC_IND_GRIDS') THEN            ! Allocate array for MPC_IND_GRIDS

         NAME = 'MPC_IND_GRIDS'
         IF (ALLOCATED(MPC_IND_GRIDS)) THEN
            DEALLOCATE (MPC_IND_GRIDS,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE IF (NAME_IN == 'GRID_ELEM_CONN_ARRAY') THEN     ! Allocate array for GRID_ELEM_CONN_ARRAY

         NAME = 'GRID_ELEM_CONN_ARRAY'
         IF (ALLOCATED(GRID_ELEM_CONN_ARRAY)) THEN
            DEALLOCATE (GRID_ELEM_CONN_ARRAY,STAT=IERR)
            IF (IERR == 0) THEN
               CALL ALLOCATED_MEMORY ( NAME, ZERO, 'DEALLOC', 'Y', CUR_MB_ALLOCATED, SUBR_NAME )
            ELSE
               WRITE(ERR,992) NAME, SUBR_NAME
               WRITE(F06,992) NAME, SUBR_NAME
               JERR = JERR + 1
            ENDIF
         ENDIF

      ELSE                                                 ! NAME not recognized, so coding error

         WRITE(ERR,915) SUBR_NAME, 'DEALLOCATED', NAME_IN
         WRITE(F06,915) SUBR_NAME, 'DEALLOCATED', NAME_IN
         FATAL_ERR = FATAL_ERR + JERR
         JERR = JERR + 1

      ENDIF

! Quit if there were errors

      IF (JERR /= 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
  915 FORMAT(' *ERROR   915: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' NAME OF ARRAY TO BE ',A,' IS INCORRECT. INPUT NAME WAS ',A)

  992 FORMAT(' *ERROR   992: CANNOT DEALLOCATE MEMORY FROM ARRAY ',A,' IN SUBROUTINE ',A)

! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_MODEL_STUF
