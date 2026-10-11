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

   MODULE MASS_MATRIX_ASSEMBLY

   USE ELEMENT_TRANSFORMATIONS, ONLY:  ELEM_TRANSFORM_LBG

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: EMP0, EMP, MGGC_MASS_MATRIX, SPARSE_MGG

   CONTAINS

      SUBROUTINE EMP0

! Provides an estimate of LTERM_MGGE, which is the size that is required for the G-set mass array for element (not grid point) mass,
! MGGE. The estimate is needed so that MGGE can be allocated.

! If field 3 of PARAM SETLKTK is 0: then the estimate of LTERM_MGGE is based
!    on full elem ME matrices unconnected (i.e. connected DOF's recounted). If the value in field 3 of the PARAM SETLKTM
!    card is PAUSE then LINK1 is PAUSE'd after subr EMP0 so user can change the estimate of LTERM_MGGE

! If field 3 of PARAM SETLKTK is 1: then the estimate of LTERM_MGGE is based
!    on Bandit bandwidth of MGG times NDOFG. If the value in field 4 of the B.D. PARAM SETLKTM card is PAUSE then LINK1 is
!    PAUSE'd after subr EMP0 so user can change the estimate of LTERM_MGGE

! If field 3 of PARAM SETLKTK is 2: then the estimate of LTERM_MGGE is based on actual elem ME matrices unconnected.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  LTERM_MGGE, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  GRIDSEQ, SETLKTM, USR_LTERM_MGG

      USE ELEMENT_LOOKUPS, ONLY       :  GET_ELGP
      USE EMG_MOD, ONLY               :  EMG
      USE ELEMENT_TRANSFORMATIONS, ONLY:  ELEM_TRANSFORM_LBG
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EMP0'





! **********************************************************************************************************************************
      IF      (SETLKTM == 0) THEN                          ! LTERM_MGG based on full elem mass matrices not connected

          CALL EMP0_0

      ELSE IF (SETLKTM == 3) THEN                          ! LTERM_MGG based on actual elem ME matrices unconnected

         CALL EMP0_3

      ELSE IF (SETLKTM == 4) THEN                          ! Use estimate from user

         LTERM_MGGE = USR_LTERM_MGG

      ENDIF



      RETURN

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE EMP0_0

! Estimates LTERM_MGGE based on full elem mass matrices in an unassembled state (not connected)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  LTERM_MGGE, NELE, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SPARSTOR, SUPINFO
      USE MODEL_STUF, ONLY            :  EDAT, EID, EPNT, ETYPE, ELGP, TYPE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EMP0_0'

      INTEGER(LONG)                   :: DELTA_LTERM_MGGE  ! Increment of LTERM_MGGE for one element
      INTEGER(LONG)                   :: I                 ! DO loop index





! **********************************************************************************************************************************
! Process the elements: Assume mass has no coupling from one grid to another


      LTERM_MGGE = 0
      DO I = 1,NELE

         CALL GET_ELGP ( I )

         IF (SPARSTOR == 'SYM   ') THEN
            DELTA_LTERM_MGGE = 3*ELGP*(6*ELGP + 1)
            LTERM_MGGE = LTERM_MGGE + DELTA_LTERM_MGGE     ! SYM has terms only on diag and above
         ELSE
            DELTA_LTERM_MGGE = (6*ELGP)*(6*ELGP)
            LTERM_MGGE = LTERM_MGGE + DELTA_LTERM_MGGE     ! NONSYM can have full matrix
         ENDIF


      ENDDO

      WRITE(ERR,4321) LTERM_MGGE, SETLKTM
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4321) LTERM_MGGE, SETLKTM
      ENDIF



      RETURN

! **********************************************************************************************************************************
 4321 FORMAT(' *INFORMATION: IN EMPO_0: ESTIMATE OF NUMBER OF NONZEROS IN MASS MATRIX MGGE IS       = ',I12,                       &
                           ' BASED ON PARAM SETLKTM = ',I3)

! **********************************************************************************************************************************

      END SUBROUTINE EMP0_0

! ##################################################################################################################################

      SUBROUTINE EMP0_3

! Estimates LTERM_MGG based on actual element mass matrices unconnected.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, LTERM_MGGE, MELDOF, NELE, NSUB
      USE PARAMS, ONLY                :  EPSIL, SETLKTM, SPARSTOR, SUPINFO
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  ELDOF, NUM_EMG_FATAL_ERRS, ME, PLY_NUM, TYPE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EMP0_3'
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option flags for subr EMG (to tell it what to calc)

      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG)                   :: KSTART            ! Index


      REAL(DOUBLE)                    :: DQE(MELDOF,NSUB)  ! Dummy array in call to ELEM_TRANSFORM_LBG
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages

      EPS1 = EPSIL(1)

! Null dummy array DQE used in call to ELEM_TRANSFORM_LBG

      DO I=1,MELDOF
         DO J=1,NSUB
            DQE(I,J) = ZERO
         ENDDO
      ENDDO

! Set up the option flags for EMG:

      OPT(1) = 'Y'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

! Process the elements:

      IERROR = 0
      LTERM_MGGE = 0
      CALL COUNTER_INIT('     Estimate size of MGG: process elem  ', NELE)
elems:DO I=1,NELE

         PLY_NUM = 0
         CALL EMG ( I   , OPT, 'N', SUBR_NAME, 'N' )       ! 'N' means do not write to BUG file

         IF (NUM_EMG_FATAL_ERRS /=0) THEN
            IERROR = IERROR + NUM_EMG_FATAL_ERRS
            CYCLE elems
         ENDIF

! Transform ME from local at the elem ends to basic at elem ends to global at elem ends to global at grids.
                                                           ! Transform PTE from local-basic-global
         IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  '))THEN
            CALL ELEM_TRANSFORM_LBG ( 'ME', ME, DQE )
         ENDIF

! Count nonzero terms in transformed ME

         DO J=1,ELDOF
            IF (SPARSTOR == 'SYM') THEN
               KSTART = J
            ELSE
               KSTART = 1
            ENDIF
            DO K=KSTART,ELDOF
               IF (DABS(ME(J,K)) < EPS1) THEN
                  CYCLE
               ELSE
                  LTERM_MGGE = LTERM_MGGE + 2
               ENDIF
            ENDDO
         ENDDO
         CALL COUNTER_PROGRESS(I)
      ENDDO elems

      WRITE(SC1,*) CR13

! Quit if IERROR > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,9876) IERROR
         WRITE(F06,9876) IERROR
         CALL OUTA_HERE ( 'Y' )                            ! IERROR is count of all subr EMG errors, so quit
      ENDIF

      WRITE(ERR,4321) LTERM_MGGE, SETLKTM
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4321) LTERM_MGGE, SETLKTM
      ENDIF



      RETURN

! **********************************************************************************************************************************
 4321 FORMAT(' *INFORMATION: IN EMPO_3: ESTIMATE OF NUMBER OF NONZEROS IN MASS MATRIX MGGE IS       = ',I12,                       &
                           ' BASED ON PARAM SETLKTM = ',I3)

 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')

! **********************************************************************************************************************************

      END SUBROUTINE EMP0_3

      END SUBROUTINE EMP0


      SUBROUTINE EMP

! Element mass processor

! EMP generates the portion of the G-set mass matrix due to element mass and puts it into the 1D array EMS of nonzero mass terms
! above the diagonal. Integer arrays EMSKEY, EMSPNT and EMSCOL are generated to form a linked list for the mass terms.


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, F22, F22FIL, F22_MSG, SC1, WRT_BUG, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_ME_BIT, ELDT_F22_ME_BIT, FATAL_ERR, IBIT, LINKNO, LTERM_MGGE,   &
                                         MBUG, MELDOF, NDOFG, NELE, NGRID, NTERM_MGGE, NSUB
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL, SPARSTOR
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  AGRID, ELDT, ELDOF, ELGP, GRID_ID, NUM_EMG_FATAL_ERRS, ME, OELDT, PLY_NUM, TYPE
      USE EMS_ARRAYS, ONLY            :  EMS, EMSCOL, EMSKEY, EMSPNT
      USE ELEMENT_TRANSFORMATIONS, ONLY:  ELEM_TRANSFORM_LBG

      USE EMG_MOD, ONLY               :  EMG
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_FIJFIL
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EMP'
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option flags for subr EMG (to tell it what to calc)

      INTEGER(LONG)                   :: EDOF(MELDOF)      ! A list of the G-set DOF's for an elem
      INTEGER(LONG)                   :: EDOF_ROW_NUM      ! Row number in array EDOF
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no. in array TDOF where G-set DOF's are kept
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: I1                ! Intermediate variable resulting from an IAND operation
      INTEGER(LONG)                   :: IDUM              ! Dummy variable used when flipping DOF's
      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IS                ! A pointer into arrays EMSKEY and EMSPNT
      INTEGER(LONG)                   :: ISS               ! A particular value of IS
      INTEGER(LONG)                   :: KSTART            ! Used in deciding whether to process all elem mass terms or only
!                                                            the ones on and above the diagonal (controlled by param SPARSTOR)
      INTEGER(LONG)                   :: MAX_NUM           ! MAX of NTERM_MGGE/NDOFG (used for DEBUG printout)
      INTEGER(LONG)                   :: MGG_ROW           ! A row no. in MGG
      INTEGER(LONG)                   :: MGG_ROWJ          ! Another row no. in EMS
      INTEGER(LONG)                   :: MGG_COL           ! A col no. in MGG
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID is a physical grid, 1 if a scalar point
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr READERR
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: TDOF_ROW_NUM      ! Row number in array TDOF
                                                           ! Indicator for output of elem data to BUG file


      REAL(DOUBLE)                    :: DQE(MELDOF,NSUB)  ! Dummy array in call to ELEM_TRANSFORM_LBG
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      INTRINSIC                       :: DABS, IAND



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

!! Null dummy array DQE used in call to ELEM_TRANSFORM_LBG

      DO I=1,MELDOF
         DO J=1,NSUB
            DQE(I,J) = ZERO
         ENDDO
      ENDDO

! Set up the option flags for EMG:

      OPT(1) = 'Y'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

! Process the elements:

      IS  = 0
      ISS = IS
      IF ((DEBUG(10) == 22) .OR. (DEBUG(10) == 23) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
         CALL DUMPEMS ( '0', 0, 0, 0, 0, 0, 0 )
      ENDIF

      IERROR = 0
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      CALL COUNTER_INIT('     Calculating mass matrix. Process elem   ', NELE)
      elems:DO I=1,NELE

         DO J=0,MBUG-1
            WRT_BUG(J) = 0
         ENDDO

         IF (LINKNO == 1) THEN                             ! Only want element mass matrix, ME, written to BUG file in LINK 1
            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_ME_BIT))       ! WRT_BUG(3): printed output of ME
            IF (I1 > 0) THEN
               WRT_BUG(3) = 1
            ENDIF
         ENDIF

         IF ((DEBUG(10) == 22) .OR. (DEBUG(10) == 23) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
            WRITE(F06,14001)
         ENDIF

         NUM_EMG_FATAL_ERRS = 0
         PLY_NUM = 0
         CALL EMG ( I   , OPT, 'N', SUBR_NAME, 'Y' )       ! 'Y' means write to BUG file
         IF (NUM_EMG_FATAL_ERRS /=0) THEN
            IERROR = IERROR + NUM_EMG_FATAL_ERRS
            CYCLE elems
         ENDIF

         I1 = IAND(OELDT,IBIT(ELDT_F22_ME_BIT))            ! Do we need to write elem mass matrices to F22 files
         IF (I1 > 0) THEN
            CALL WRITE_FIJFIL ( 2, 0 )
         ENDIF

         EDOF_ROW_NUM = 0                                  ! Generate element DOF'S
         DO J = 1,ELGP
!           CALL CALC_TDOF_ROW_NUM ( AGRID(J), ROW_NUM_START, 'N' )
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(J), IGRID )
            ROW_NUM_START = TDOF_ROW_START(IGRID)
            CALL GET_GRID_NUM_COMPS ( IGRID, NUM_COMPS, SUBR_NAME )
            DO K = 1,NUM_COMPS
               CALL TDOF_COL_NUM ( 'G ',  G_SET_COL_NUM )
               TDOF_ROW_NUM       = ROW_NUM_START + K - 1
               EDOF_ROW_NUM       = EDOF_ROW_NUM + 1
               EDOF(EDOF_ROW_NUM) = TDOF(TDOF_ROW_NUM, G_SET_COL_NUM)
            ENDDO
         ENDDO

! Transform ME from local at the elem ends to basic at elem ends to global at elem ends to global at grids.

                                                           ! Transform PTE from local-basic-global
         IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  '))THEN

            CALL ELEM_TRANSFORM_LBG ( 'ME', ME, DQE )
24357 format(6(1es14.6))

         ENDIF

! Put the element mass matrix, ME, into EMS array. J ranges over rows, K over cols of elem mass matrix, ME

mgg_rows:DO J = 1,ELDOF
            MGG_ROWJ  = EDOF(J)
            IF ((DEBUG(10) == 22) .OR. (DEBUG(10) == 23) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
               WRITE(F06,*)
            ENDIF

            IF (SPARSTOR == 'SYM') THEN                    ! Set KSTART depending on SPARSTOR
               KSTART = J                                  ! Process only upper right portion of ME
            ELSE
               KSTART = 1                                  ! Process all of ME
            ENDIF

mgg_cols:   DO K = KSTART,ELDOF
               MGG_ROW  = MGG_ROWJ                         ! Make sure we have correct row num. It may have been flipped w/ col
               MGG_COL  = EDOF(K)
               IF (DABS(ME(J,K)) < EPS1) THEN
                  CYCLE mgg_cols
               ENDIF

               IF (SPARSTOR == 'SYM') THEN                 ! If 'SYM', Flip MGG_COL,MGG_ROW if MGG_COL < MGG_ROW
                  IF (MGG_COL < MGG_ROW) THEN
                     IDUM    = MGG_ROW
                     MGG_ROW = MGG_COL
                     MGG_COL = IDUM
                  ENDIF
               ENDIF

               IS = EMSKEY(MGG_ROW)                        ! Get pointer to first term in row MGG_ROW of global mass matrix
               IF (IS == 0) THEN                           ! EMSKEY(MGG_ROW)=0 means no current terms in global mass matrix at row
                                                           ! MGG_ROW update NTERM_MGGE and reset EMSKEY, EMSCOL, EMSPNT, EMS arrays
                  NTERM_MGGE = NTERM_MGGE + 1

                  IF (NTERM_MGGE > LTERM_MGGE) THEN
                     WRITE(ERR,1624) SUBR_NAME, 'MASS    ', 'LTERM_MGGE', LTERM_MGGE
                     WRITE(F06,1624) SUBR_NAME, 'MASS    ', 'LTERM_MGGE', LTERM_MGGE
                     CALL OUTA_HERE ( 'Y' )                        ! MYSTRAN limitation, so quit
                  ENDIF

                  EMSKEY(MGG_ROW) = NTERM_MGGE
                  EMSCOL(NTERM_MGGE) = MGG_COL
                  EMSPNT(NTERM_MGGE) = 0
                  EMS(NTERM_MGGE) = ME(J,K)

                  IF ((DEBUG(10) == 22) .OR. (DEBUG(10) == 23) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
                     CALL DUMPEMS ( 'A', J, K, MGG_ROW, MGG_COL, IS, ISS )
                  ENDIF

               ELSE                                        ! EMSKEY(MGG_ROW) /= 0 means there are already some terms in row MGG_ROW

emspnt0:          DO                                       ! so, run this loop until we find a place to put ME(J,K). If there is
                                                           ! already a term in this row w/ same DOF's as ME(J,K), loop runs once.
                                                           ! If not, then this loop runs until it finds EMSPNT=0, and insetrs term.

                     IF (MGG_COL == EMSCOL(IS)) THEN       ! There is a term that exists with same DOF'S as ME(J,K) so add terms

                        EMS(IS) = EMS(IS) + ME(J,K)
                        IF ((DEBUG(10) == 22) .OR. (DEBUG(10) == 23) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
                           CALL DUMPEMS ( 'B', J, K, MGG_ROW, MGG_COL, IS, ISS )
                        ENDIF

                        CYCLE mgg_cols                     ! We have added a term to EMS so exit this loop and do next col of MGG

                     ELSE                                  ! This is a new term for row J. Need to cycle until we find EMSPNT = 0.
                                                           ! Then we can put ME(J,K) in EMS
                        ISS = IS
                        IS  = EMSPNT(IS)
                        IF (IS == 0) THEN                  ! We are at end of where terms are in this row, so ME(J,K) goes here
                           IF (NTERM_MGGE+1 > LTERM_MGGE) THEN
                              WRITE(ERR,1624) SUBR_NAME, 'MASS', 'LTERM_MGGE', LTERM_MGGE
                              WRITE(F06,1624) SUBR_NAME, 'MASS', 'LTERM_MGGE', LTERM_MGGE
                              CALL OUTA_HERE ( 'Y' )       ! MYSTRAN limitation, so quit
                           ENDIF
                           NTERM_MGGE        = NTERM_MGGE+1! Increment NTERM_MGGE
                           EMSPNT(ISS)       = NTERM_MGGE  ! EMSPNT for the current ME(J,K) term
                           EMSPNT(NTERM_MGGE) = 0          ! Latest EMSPNT is set to 0 so we will know when to insert next ME(J,K)
                           EMSCOL(NTERM_MGGE) = MGG_COL    ! EMSCOL always is MGG_COL
                           EMS   (NTERM_MGGE) = ME(J,K)
                           IF ((DEBUG(10) == 22) .OR. (DEBUG(10) == 23) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
                              CALL DUMPEMS ( 'C', J, K, MGG_ROW, MGG_COL, IS, ISS )
                           ENDIF

                           CYCLE mgg_cols                  ! We put ME(J,K) into EMS so exit this loop and do next col of MGG
                        ELSE                               ! EMSPNT /= 0 so cycle this loop until we get it = 0
                           CYCLE emspnt0
                        ENDIF

                     ENDIF

                  ENDDO emspnt0

               ENDIF

            ENDDO mgg_cols

         ENDDO mgg_rows
         CALL COUNTER_PROGRESS(I)

      ENDDO elems

      WRITE(SC1,*) CR13

! Debug output:

      IF((DEBUG(10) == 21) .OR. (DEBUG(10) == 22) .OR. (DEBUG(10) == 23) .OR.                                                     &
         (DEBUG(10) == 31) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
         WRITE(F06,1260)
         MAX_NUM = MAX(NTERM_MGGE,NDOFG)
         DO I=1,MAX_NUM
            IF      (MAX_NUM == NTERM_MGGE) THEN
               IF (NDOFG >= I) THEN
                  WRITE(F06,1261) I,EMSKEY(I),EMSCOL(I),EMSPNT(I),EMS(I)
               ELSE
                  WRITE(F06,1262) I,          EMSCOL(I),EMSPNT(I),EMS(I)
               ENDIF
            ELSE IF (MAX_NUM == NDOFG) THEN
               IF (NTERM_MGGE >= I) THEN
                  WRITE(F06,1261) I,EMSKEY(I),EMSCOL(I),EMSPNT(I),EMS(I)
               ELSE
                  WRITE(F06,1263) I,EMSKEY(I)
               ENDIF
            ENDIF
         ENDDO
         WRITE(F06,*)
      ENDIF

! Reset subr EMG option flags:

      OPT(3) = 'N'
      OPT(4) = 'N'

! Quit if IERROR > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,9876) IERROR
         WRITE(F06,9876) IERROR
         CALL OUTA_HERE ( 'Y' )                                    ! IERROR is count of all subr EMG errors, so quit
      ENDIF


      RETURN

! **********************************************************************************************************************************
 1260 FORMAT(/,'            I   EMSKEY(I)   EMSCOL(I)   EMSPNT(I)           EMS(I)')

 1261 FORMAT(1X,I12,I12,I12,I12,3X,1ES21.14)

 1262 FORMAT(1X,I12,12X,I12,I12,3X,1ES21.14)

 1263 FORMAT(1X,I12,I12)

 1624 FORMAT(' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY NON-ZERO TERMS IN THE ',A,' MATRIX. LIMIT IS ',A,' = ',I12)

 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')

14001 FORMAT(' ******************************************************************************************************************&
&******')

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DUMPEMS ( WHAT, J, K, MGG_ROW, MGG_COL, IS, ISS )

! Prints out info on the formulation of stiffness arrays for subr ESP, which generates the arrays

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  NTERM_MGGE
      USE MODEL_STUF, ONLY            :  EID
      USE EMS_ARRAYS, ONLY            :  EMS, EMSCOL, EMSKEY, EMSPNT

      IMPLICIT NONE

      CHARACTER(1*BYTE), INTENT(IN)   :: WHAT              ! Indicator of where this subr was called from in subr EMP

      INTEGER(LONG)                   :: IS                ! A pointer into arrays EMSKEY and EMSPNT
      INTEGER(LONG)                   :: ISS               ! A particular value of IS
      INTEGER(LONG)    , INTENT(IN)   :: J                 ! Row number of elem mass matrix term, ME(J,K)
      INTEGER(LONG)    , INTENT(IN)   :: K                 ! Col number of elem mass matrix term, ME(J,K)
      INTEGER(LONG)    , INTENT(IN)   :: MGG_COL           ! Row number of MGG matrix where ME(J,K) goes
      INTEGER(LONG)    , INTENT(IN)   :: MGG_ROW           ! Col number of MGG matrix where ME(J,K) goes

! **********************************************************************************************************************************
      IF      (WHAT == '0') THEN

         WRITE(F06,8910)

      ELSE IF (WHAT == 'A') THEN

         WRITE(F06,8930) EID, J, K, MGG_ROW, MGG_COL, EMSKEY(MGG_ROW), NTERM_MGGE, EMSCOL(NTERM_MGGE), EMSPNT(NTERM_MGGE),         &
                         EMS(NTERM_MGGE), IS, ISS

      ELSE IF (WHAT == 'B') THEN

         IF (ISS /= 0) THEN
         WRITE(F06,8940) EID, J, K, MGG_ROW, MGG_COL, EMSKEY(MGG_ROW), NTERM_MGGE, EMSCOL(NTERM_MGGE), EMSPNT(NTERM_MGGE),        &
                         EMS(NTERM_MGGE), IS, ISS, EMSPNT(ISS), EMS(IS)

         ELSE
         WRITE(F06,8941) EID, J, K, MGG_ROW, MGG_COL, EMSKEY(MGG_ROW), NTERM_MGGE, EMSCOL(NTERM_MGGE), EMSPNT(NTERM_MGGE),        &
                         EMS(NTERM_MGGE), IS, ISS,              EMS(IS)
         ENDIF

      ELSE IF (WHAT == 'C') THEN

         WRITE(F06,8950) EID, J, K, MGG_ROW, MGG_COL, EMSKEY(MGG_ROW), NTERM_MGGE, EMSCOL(NTERM_MGGE), EMSPNT(NTERM_MGGE),         &
                         EMS(NTERM_MGGE), IS, ISS

      ENDIF

      RETURN

! **********************************************************************************************************************************
 8910 FORMAT(4X,'ELEM       J       K MGG_ROW MGG_COL  EMSKEY  NKTERM  EMSCOL  EMSPNT      EMS         IS     ISS  EMSPNT     EMS'&
          ,/,4X,' ID                                 (MGG_ROW)       (NKTERM) (NKTERM)   (NKTERM)                   (ISS)     (IS)')

 8920 FORMAT(1X,'---------------------------------------------------------------------------------------------------------------', &
'--------')

 8930 FORMAT(1X,'A',I6,I8,I8,I8,I8,I8,I8,I8,I8,1ES12.3,I8,I8)

 8940 FORMAT(1X,'B',I6,I8,I8,I8,I8,I8,I8,I8,I8,1ES12.3,I8,I8,I8,1ES12.3)

 8941 FORMAT(1X,'B',I6,I8,I8,I8,I8,I8,I8,I8,I8,1ES12.3,I8,I8,' -------',1ES12.3)

 8950 FORMAT(1X,'C',I6,I8,I8,I8,I8,I8,I8,I8,I8,1ES12.3,I8,I8)

! **********************************************************************************************************************************
      END SUBROUTINE DUMPEMS

      END SUBROUTINE EMP


      SUBROUTINE MGGC_MASS_MATRIX

! Forms the mass matrix, MGGC, for concentrated masses by calling subr MGG_CONM2_PROC to process the concentrated masses

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  NGRID, NTERM_MGGC, BLNK_SUB_NAM
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  EPSIL, SPARSTOR, WTMASS
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  AGRID, GRID_ID, INV_GRID_SEQ
      USE SPARSE_MATRICES, ONLY       :  I_MGGC, J_MGGC, MGGC

      USE LINK1_WORKSPACE_LIFECYCLE, ONLY:  ALLOCATE_L1_MGG
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1, GET_GRID_NUM_COMPS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MGGC_MASS_MATRIX'
      CHARACTER( 1*BYTE)              :: MGG_CONM2_NONZERO ! 'Y'/'N' indicator if a nonzero MGG_CONM2 6 x 6 matrix was created

      INTEGER(LONG)                   :: GRID_NUM          ! The actual grid number for which we create a 6 x 6 mass matrix for
!                                                            one CONM2 or 1 element (if one exists for this grid)
      INTEGER(LONG)                   :: DELTA_KTERM_MGGC  ! Coumt of nonzero terms in MGGC array for 1 grid
      INTEGER(LONG)                   :: KTERM_MGGC        ! Coumt of nonzero terms in MGGC array
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: IJ                ! Index
      INTEGER(LONG)                   :: IROW_START        ! Row number in TDOF where data begins for IGRID
      INTEGER(LONG)                   :: NUM_COMPS         ! Number of displ components (1 for SPOINT, 6 for physical grid)
      INTEGER(LONG)                   :: KSTART            ! Used in deciding whether to process all elem mass terms or only
!                                                            the ones on and above the diagonal (controlled by param SPARSTOR)
      INTEGER(LONG)                   :: MGGC_COL_NUM      ! A calculated col number for a nonzero term in MGG arrays


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: MGG_CONM2(6,6)    ! 6 X 6 mass matrix in global coords for one CONM2

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Initialize

      NTERM_MGGC = 36*NGRID                                ! Max possible NTERM_MGGC. Will be set to actual later by count nonzeros
      CALL ALLOCATE_L1_MGG ( 'MGGC', SUBR_NAME )

      KTERM_MGGC = 0
      I_MGGC(1)  = 1
      IROW_START = 1
i_do1:DO I=1,NGRID

         GRID_NUM = GRID_ID(INV_GRID_SEQ(I))               ! GRID_NUM's are in TDOFI order (internal DOF order)
         CALL MGG_CONM2_PROC ( I, GRID_NUM, MGG_CONM2, MGG_CONM2_NONZERO )
         CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )

         IF (MGG_CONM2_NONZERO == 'Y') THEN
            DO J=1,NUM_COMPS
               DELTA_KTERM_MGGC = 0

               IF (SPARSTOR == 'SYM') THEN                 ! Set KSTART depending on SPARSTOR
                  KSTART = J                               ! Process only upper right portion of MGG_CONM2
               ELSE
                  KSTART = 1                               ! Process all of MGG_CONM2
               ENDIF
               DO K=KSTART,NUM_COMPS
                  IF (DABS(MGG_CONM2(J,K)) >= EPS1) THEN
                     MGGC_COL_NUM       = IROW_START + K - 1
                     KTERM_MGGC         = KTERM_MGGC + 1
                     DELTA_KTERM_MGGC   = DELTA_KTERM_MGGC + 1
                     IF (KTERM_MGGC > NTERM_MGGC) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, KTERM_MGGC, 'MGGC' )
                      J_MGGC(KTERM_MGGC) = MGGC_COL_NUM
                        MGGC(KTERM_MGGC) = MGG_CONM2(J,K)
                  ENDIF
               ENDDO
!xx            IJ = 6*(I-1) + J
               IJ = NUM_COMPS*(I-1) + J
               I_MGGC(IJ+1) = I_MGGC(IJ) + DELTA_KTERM_MGGC
            ENDDO

         ELSE

            DO J=1,NUM_COMPS
               IJ = IROW_START + J - 1
               I_MGGC(IJ+1) = I_MGGC(IJ)
            ENDDO

         ENDIF

         IROW_START = IROW_START + NUM_COMPS

      ENDDO i_do1

      NTERM_MGGC = KTERM_MGGC

! Do not deallocate MGGC arrays - they are needed later in subr SPARSE_MGG



      RETURN

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE MGG_CONM2_PROC ( INT_GRID_ID, GRID_NUM, MGG_CONM2, MGG_CONM2_NONZERO )

! Generates 6 x 6 mass matrix, MGG_CONM2, for one CONM2 for grid GRID_NUM (if there is any CONM2 connected to this grid)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  NCONM2, NGRID, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  CONM2, RCONM2
      USE PARAMS, ONLY                :  ART_MASS, ART_ROT_MASS, ART_TRAN_MASS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MGG_CONM2_PROC'
      CHARACTER( 1*BYTE), INTENT(OUT) :: MGG_CONM2_NONZERO ! 'Y'/'N' indicator if a nonzero MGG_CONM2 6 x 6 matrix was created

      INTEGER(LONG), INTENT(IN)       :: INT_GRID_ID       ! The internal grid number for which we create a 6 x 6 mass matrix for
!                                                            one CONM2 (if one exists for this grid)
      INTEGER(LONG), INTENT(IN)       :: GRID_NUM          ! The actual grid number for internal grid ID INT_GRID_ID
      INTEGER(LONG)                   :: I,J,L             ! DO loop indices or counters


      REAL(DOUBLE) , INTENT(OUT)      :: MGG_CONM2(6,6)    ! 6 X 6 mass matrix in global coords for one CONM2



! **********************************************************************************************************************************
! Initialize

      MGG_CONM2_NONZERO = 'N'

      DO I=1,6
         DO J=1,6
            MGG_CONM2(I,J) = ZERO
         ENDDO
      ENDDO

! Add artificial mass terms to the grid, if requested

      IF (ART_MASS == 'Y') THEN

         DO I=1,3
            MGG_CONM2(I,I) = ART_TRAN_MASS
         ENDDO
         MGG_CONM2_NONZERO = 'Y'

         DO I=4,6
            MGG_CONM2(I,I) = ART_ROT_MASS
         ENDDO
         MGG_CONM2_NONZERO = 'Y'

      ENDIF

! Process CONM2's for this grid

      DO L=1,NCONM2

         IF (CONM2(L,2) == GRID_NUM) THEN

            WRITE(SC1,12345,ADVANCE='NO') CONM2(L,1), INT_GRID_ID, NGRID, CR13
            WRITE(SC1,*) CR13

            MGG_CONM2_NONZERO = 'Y'

            MGG_CONM2(1,1) =  MGG_CONM2(1,1) + RCONM2(L, 1)
            MGG_CONM2(2,2) =  MGG_CONM2(2,2) + RCONM2(L, 1)
            MGG_CONM2(3,3) =  MGG_CONM2(3,3) + RCONM2(L, 1)

            MGG_CONM2(1,5) =  MGG_CONM2(1,5) + RCONM2(L, 1)*RCONM2(L,4)
            MGG_CONM2(1,6) =  MGG_CONM2(1,6) - RCONM2(L, 1)*RCONM2(L,3)

            MGG_CONM2(2,4) =  MGG_CONM2(2,4) - RCONM2(L, 1)*RCONM2(L,4)
            MGG_CONM2(2,6) =  MGG_CONM2(2,6) + RCONM2(L, 1)*RCONM2(L,2)

            MGG_CONM2(3,4) =  MGG_CONM2(3,4) + RCONM2(L, 1)*RCONM2(L,3)
            MGG_CONM2(3,5) =  MGG_CONM2(3,5) - RCONM2(L, 1)*RCONM2(L,2)

            MGG_CONM2(4,4) =  MGG_CONM2(4,4) + RCONM2(L, 5)
            MGG_CONM2(4,5) =  MGG_CONM2(4,5) - RCONM2(L, 6)
            MGG_CONM2(4,6) =  MGG_CONM2(4,6) - RCONM2(L, 8)

            MGG_CONM2(5,5) =  MGG_CONM2(5,5) + RCONM2(L, 7)
            MGG_CONM2(5,6) =  MGG_CONM2(5,6) - RCONM2(L, 9)

            MGG_CONM2(6,6) =  MGG_CONM2(6,6) + RCONM2(L,10)

            DO I=1,6
               DO J=1,I-1
                  MGG_CONM2(I,J) = MGG_CONM2(J,I)
               ENDDO
            ENDDO

         ENDIF

      ENDDO



      RETURN

! **********************************************************************************************************************************
12345 format(5X,'Process mass for CONM2 ',I8,' int G.P. ',I8,' of ',I8,'           ', A)

! **********************************************************************************************************************************

      END SUBROUTINE MGG_CONM2_PROC

      END SUBROUTINE MGGC_MASS_MATRIX


      SUBROUTINE MGGS_MASS_MATRIX

! Forms the sparse scalar mass matrix, MGGS, (for masses defined on Bulk Data CMASS)
! Each mass is at the G-set DOF of its point and component (CMASS column 5 or 7; 0 or 1 for a scalar point). Before, every mass
! was put at the DOF of the first row whose component label was 1 counted in GRID_ID order: the component of a CMASS1/CMASS2 on a
! grid was ignored (a mass on T3 went to T1, Agent 2 a2_m13: MLL singular), and with the TDOF labels wrong for points sequenced out
! of ID order (fixed in TDOF_PROC) the masses of scalar points went to other DOFs. Masses at the same DOF are now added (before,
! only the first one counted).

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, NCMASS, NDOFG, NGRID, NPMASS, NTERM_MGGS, BLNK_SUB_NAM
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  SPARSTOR, WTMASS
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  CMASS, GRID_ID, PMASS, RPMASS
      USE SPARSE_MATRICES, ONLY       :  I_MGGS, J_MGGS, MGGS

      USE LINK1_WORKSPACE_LIFECYCLE, ONLY:  ALLOCATE_L1_MGG
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1, GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SORTING, ONLY               :  SORT_INT2_REAL1

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MGGS_MASS_MATRIX'
      CHARACTER( 1*BYTE)              :: FOUND             ! 'Y'/'N' indicator of whether we found something

      INTEGER(LONG)                   :: COMP              ! Component of a scalar mass at its point
      INTEGER(LONG)                   :: NCOMPS            ! Number of components of that point (1 scalar point, 6 grid)
      INTEGER(LONG)                   :: G_SET_COL         ! Col in TDOF where G-set exists
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: IERROR            ! Local error count
      INTEGER(LONG)                   :: IDOF(NCMASS)      ! G-set DOF number
      INTEGER(LONG)                   :: KTERM_MGGS        ! Count of number of terma going into MGGS
      INTEGER(LONG)                   :: ROW_NUM           ! Row number in TDOF where data begins for IGRID
      INTEGER(LONG)                   :: SGRID(NCMASS)     ! Grid number for a scalar mass (from array CMASS)
      INTEGER(LONG)                   :: PMASS_ID(NCMASS)  ! Prop ID for the CMASS that is attached to SGRID(I) (from array PMASS)


      REAL(DOUBLE)                    :: PMASS_VAL(NCMASS) ! Value for the mass attached to SGRID(I)

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
! Initialize

      IERROR = 0

      NTERM_MGGS = NDOFG+1                                 ! Max possible NTERM_MGGS. Will be set to actual later by nonzero count
      CALL ALLOCATE_L1_MGG ( 'MGGS', SUBR_NAME )

      CALL TDOF_COL_NUM ( 'G ', G_SET_COL )

      DO I=1,NCMASS

         IF (CMASS(I,4) /= 0) THEN                         ! The scalar point is in either col 4 or 6 in CMASS (checked in BD read)
            SGRID(I) = CMASS(I,4)
            COMP     = CMASS(I,5)
         ELSE
            SGRID(I) = CMASS(I,6)
            COMP     = CMASS(I,7)
         ENDIF
         PMASS_ID(I) = CMASS(I,3)

         ROW_NUM = -1
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, SGRID(I), ROW_NUM )
         IF (ROW_NUM /= -1) THEN
            CALL GET_GRID_NUM_COMPS ( ROW_NUM, NCOMPS, SUBR_NAME )
            IF ((NCOMPS == 1) .AND. (COMP == 0)) COMP = 1
            IF ((COMP < 1) .OR. (COMP > NCOMPS)) THEN
               IERROR    = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,1362) CMASS(I,1), SGRID(I), COMP, NCOMPS
               WRITE(F06,1362) CMASS(I,1), SGRID(I), COMP, NCOMPS
               IDOF(I) = 0
            ELSE
               IDOF(I) = TDOF(TDOF_ROW_START(ROW_NUM)+COMP-1,G_SET_COL)
            ENDIF
         ELSE
            IERROR    = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1361) 'GRID OR SPOINT', SGRID, 'BULK DATA CMASS ENTRY'
            WRITE(F06,1361) 'GRID OR SPOINT', SGRID, 'BULK DATA CMASS ENTRY'
         ENDIF

      ENDDO

! Get mass value at each SGRID (PMASS_VAL array)

i_do1:DO I=1,NCMASS
         FOUND = 'N'
j_do1:   DO J=1,NPMASS
            IF (PMASS_ID(I) == PMASS(J,1)) THEN
               PMASS_VAL(I) = RPMASS(J,1)
               FOUND = 'Y'
               EXIT j_do1
            ENDIF
         ENDDO j_do1
         IF (FOUND == 'N') THEN
            IERROR    = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1601) PMASS_ID(I), CMASS(I,1)
            WRITE(ERR,1601) PMASS_ID(I), CMASS(I,1)
         ENDIF
      ENDDO i_do1

      IF (DEBUG(182) > 0) CALL DEB_MGGS ( 1 )

! Quit if IERROR > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,9999) IERROR
         WRITE(F06,9999) IERROR
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Sort arrays IDOF, SGRID, PMASS_VAL so that IDOF is in numerical order (G-set DOF order)

      CALL SORT_INT2_REAL1 ( SUBR_NAME, 'IDOF, SGRI, PMASS_VAL', NCMASS, IDOF, SGRID, PMASS_VAL )
      IF (DEBUG(182) > 0) CALL DEB_MGGS ( 2 )

! Formulate sparse matrix MGGS. Note that the MGGS matrix is diagonal since CMASS is attached to only 1 grid/scalar point

      KTERM_MGGS = 0                                       ! IDOF is sorted: walk it once, adding the masses at the same DOF
      I_MGGS(1) = 1
      J = 1
i_do2:DO I=1,NDOFG
         I_MGGS(I+1) = I_MGGS(I)
         DO WHILE (J <= NCMASS)
            IF (IDOF(J) > I) EXIT
            IF (IDOF(J) == I) THEN
               IF (I_MGGS(I+1) == I_MGGS(I)) THEN          ! The first mass at this DOF
                  I_MGGS(I+1)        = I_MGGS(I) + 1
                  KTERM_MGGS         = KTERM_MGGS + 1
                  IF (KTERM_MGGS > NTERM_MGGS) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, NTERM_MGGS, 'MGGS' )
                  J_MGGS(KTERM_MGGS) = I                   ! Since MGGS is diagonal
                  MGGS(KTERM_MGGS)   = PMASS_VAL(J)
               ELSE
                  MGGS(KTERM_MGGS)   = MGGS(KTERM_MGGS) + PMASS_VAL(J)
               ENDIF
            ENDIF
            J = J + 1
         ENDDO
      ENDDO i_do2

      IF (DEBUG(182) > 0) CALL DEB_MGGS ( 3 )

! Reset NTERM_MGGS to what was counted above

      NTERM_MGGS = KTERM_MGGS



      RETURN

! **********************************************************************************************************************************
 1361 FORMAT(' *ERROR  1361: UNDEFINED ',A,I8,' ON ',A)

 1362 FORMAT(' *ERROR  1362: SCALAR MASS ',I8,' IS ON POINT ',I8,' COMPONENT ',I2,', BUT THE POINT HAS ',I2,' COMPONENT(S)',    &
             ' (A SCALAR POINT TAKES 0 OR BLANK, A GRID 1-6)')

 1601 FORMAT(' *ERROR  1601: UNDEFINED SCALAR MASS PROPERTY ID = ',I8,' ON CMASS ID ',I8)

 9999 FORMAT(' PROCESSING TERMINATED DUE TO ABOVE ',I8,' ERRORS')


! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEB_MGGS ( WHAT )

      USE PENTIUM_II_KIND, ONLY       :  LONG

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: WHAT              ! What to print on this call
      INTEGER(LONG)                   :: II                ! DO loop index

! **********************************************************************************************************************************
      IF (WHAT == 1) THEN
         WRITE(F06,1001)
         WRITE(F06,*) 'I, IDOF, SGRID, PMASS_ID before sorting IDOF, SGRID, PMASS_VAL on IDOF'
         DO II=1,NCMASS
            WRITE(F06,1002) II, IDOF(II), SGRID(II), PMASS_ID(II)
         ENDDO
         WRITE(F06,*)
      ENDIF

      IF (WHAT == 2) THEN
         WRITE(F06,*) 'I, IDOF(I), SGRID(I), PMASS_VAL(I) after  sort on IDOF'
         DO II=1,NCMASS
            WRITE(F06,2001) II, IDOF(II), SGRID(II), PMASS_VAL(II)
         ENDDO
         WRITE(F06,*)
      ENDIF

      IF (WHAT == 3) THEN
         WRITE(F06,*) ' Sparse arrays for MGGS'
         DO II=1,NDOFG+1
            WRITE(F06,3001) II, I_MGGS(II)
         ENDDO
         WRITE(F06,*)
         DO II=1,KTERM_MGGS
            WRITE(F06,3002) II, J_MGGS(II), MGGS(II)
         ENDDO
         WRITE(F06,*)
         WRITE(F06,3999)
      ENDIF

! **********************************************************************************************************************************
 1001 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::START DEBUG(182) OUTPUT FROM SUBROUTINE MGGS_MASS_MATRIX::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1002 format(' In MGGS_MASS_MATRIX: I, IDOF(I), SGRID(I), PMASS_ID(I)  = ',4i8)

 2001 format(' In MGGS_MASS_MATRIX: I, IDOF(I), SGRID(I), PMASS_VAL(I) = ',3i8,1es14.6)

 3001 format(' In MGGS_MASS_MATRIX: I, I_MGGS(I)                       = ',2i8)

 3002 format(' In MGGS_MASS_MATRIX: I, J_MGGS(I), MGGS(I)              = ',2i8,1es14.6)

 3999 FORMAT(' ::::::::::::::::::::::::::::::::::::::::END DEBUG(15) OUTPUT FROM SUBROUTINE CONM2_PROC_1:::::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE DEB_MGGS

      END SUBROUTINE MGGS_MASS_MATRIX


      SUBROUTINE SPARSE_MGG

! Add sparse arrays for concentrated masses (array MGGC), scalar masses (array MGGS) and element mass (array EMS) to get the final
! sparse G-set mass matrix, MGG. Rows are sorted to be in numerical G-set DOF order and the final MGG is written to file LINK1R

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1R, L1R_MSG, LINK1R, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NCMASS, NDOFG, NGRID, NTERM_MGG, NTERM_MGGC, NTERM_MGGE,         &
                                         NTERM_MGGS, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE DOF_TABLES,ONLY             :  TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  GRID_ID
      USE PARAMS, ONLY                :  EPSIL, PRTMASS, SUPINFO, WTMASS
      USE EMS_ARRAYS, ONLY            :  EMS, EMSCOL, EMSKEY, EMSPNT
      USE SPARSE_MATRICES, ONLY       :  I2_MGG, I_MGG, J_MGG, MGG, I_MGGC, J_MGGC, MGGC, I_MGGE, J_MGGE, MGGE,                    &
                                         I_MGGS, J_MGGS, MGGS,  SYM_MGGC, SYM_MGGE, SYM_MGGS
      USE SCRATCH_MATRICES, ONLY      :  I_CRS1, J_CRS1, CRS1

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, FILE_OPEN
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1, GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE SORTING, ONLY               :  SORT_INT1_REAL1
      USE SPARSE_MATRIX_ALGEBRA, ONLY :  MATADD_SSS, MATADD_SSS_NTERM
      USE SCRATCH_MATRIX_LIFECYCLE, ONLY:  ALLOCATE_SCR_CRS_MAT, DEALLOCATE_SCR_MAT
      USE LINK1_WORKSPACE_LIFECYCLE, ONLY:  ALLOCATE_L1_MGG
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS
      USE LINK1_LOAD_SUPPORT, ONLY    :  GET_GRID_6X6_MASS
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_MGG'
      CHARACTER(  1*BYTE)             :: FOUND             ! 'Y' if there is a mass matrix for this grid and 'N' otherwise
      CHARACTER(LEN=LEN(SYM_MGGE))    :: SYM_CRS1          ! 'Y'/'N' Symmetry indicator for scratch matrix CRS1

      INTEGER(LONG)                   :: GRID_NUM          ! An actual grid ID
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IERR              ! Local error count
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IK                ! Index for array I_MGGE
      INTEGER(LONG)                   :: IS                ! Index into arrays EMSPNT, EMSLIS, EMSCOL
      INTEGER(LONG)                   :: KTERM_MGGE        ! Count of terms written to MGG file LINK1R to compare with NTERM_MGGE
      INTEGER(LONG)                   :: MAX_NUM_IN_ROW    ! largest number of terms in any row of the MGG mass matrix
      INTEGER(LONG)                   :: NTERM_CRS1        ! Count of nonzero terms in matrix CRS1
      INTEGER(LONG)                   :: NUM               ! Count of the actual number of nonzero terms in a row of MGG
      INTEGER(LONG)                   :: NUM_IN_ROW_I      ! Number of nonzero terms in a row of MGG
      INTEGER(LONG)                   :: NUM_COMPS         ! Number of displ components (1 for SPOINT, 6 for physical grid)
      INTEGER(LONG)                   :: NZERO   = 0       ! Count on zero terms in array EMS
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: RJ(NDOFG)         ! Column numbers corresponding to the terms in REMS(I).
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero
      REAL(DOUBLE)                    :: GRID_MGG(6,6)     ! 6 x 6 mass matrix for a grid
      REAL(DOUBLE)                    :: REMS(NDOFG)       ! 1D array of the terms from EMS(I) pertaining to one row of the G-set
!                                                            mass matrix. Initially, the cols are not in increasing global DOF
!                                                            order. REMS is sorted, prior to writing the G-set mass matrix
!                                                            to file LINK1R, so that the cols are in increasing DOF order.

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)
! Pass # 1: Determine final NTERM_MGGE (may be less due to zero terms)

      NZERO = 0
i_do0:DO I = 1,NDOFG                                       ! Start conversion.

         IS = EMSKEY(I)
         IF (IS == 0) CYCLE i_do0                          ! Check for null row in mass matrix and CYCLE if it is

         NUM   = 0                                         ! Count zero terms so we can debit NTERM_MGGE before writing it to file
j_do0:   DO J = 1,NDOFG
            IF (DABS(EMS(IS)) < EPS1) THEN
               NZERO = NZERO + 1
            ELSE
               NUM = NUM + 1
            ENDIF
            IS = EMSPNT(IS)
            IF (IS == 0) THEN
               EXIT j_do0
            ENDIF
         ENDDO j_do0
         IF (IS /= 0) THEN
            WRITE(ERR,1626) SUBR_NAME,I
            WRITE(F06,1626) SUBR_NAME,I
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                         ! Coding error, so quit
         ENDIF

      ENDDO i_do0

      NTERM_MGGE = NTERM_MGGE - NZERO

      WRITE(ERR,146) NTERM_MGGE
      IF (SUPINFO == 'N') THEN
         WRITE(F06,146) NTERM_MGGE
      ENDIF

! **********************************************************************************************************************************
! Pass # 2: Reformulate rows and write to file LINK1R

! Open L1R to write mass.

      OUNT(1) = ERR
      OUNT(2) = F06
      CALL FILE_OPEN ( L1R, LINK1R, OUNT, 'REPLACE', L1R_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )

      KTERM_MGGE = 0
      I_MGGE(1) = 1
      WRITE(SC1, * )
      CALL COUNTER_INIT('     Working on grid ', NGRID)
i_do: DO I = 1,NGRID

         GRID_NUM = GRID_ID(I)

!xx      CALL CALC_TDOF_ROW_NUM ( GRID_NUM, IROW_START, 'N' )
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_NUM, IGRID )
         ROW_NUM_START = TDOF_ROW_START(IGRID)
         CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
k_do:    DO K=1,NUM_COMPS

            IK = ROW_NUM_START + K - 1
            IS = EMSKEY(IK)

            IF (IS == 0) THEN                              ! Check for null row in mass matrix
               I_MGGE(IK+1) = I_MGGE(IK)
               CYCLE k_do
            ENDIF

            NUM = 0
j_do1:      DO J=1,NDOFG
               IF (DABS(EMS(IS)) >= EPS1) THEN
                  NUM = NUM + 1
                  REMS(NUM) = EMS(IS)
                  RJ(NUM)   = EMSCOL(IS)
               ENDIF
               IS = EMSPNT(IS)
               IF (IS == 0) THEN
                  EXIT j_do1
               ENDIF
            ENDDO j_do1

            I_MGGE(IK+1) = I_MGGE(IK) + NUM

            IF (IS /= 0) THEN
               WRITE(ERR,1626) SUBR_NAME,I
               WRITE(F06,1626) SUBR_NAME,I
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Coding error, so quit
            ENDIF

            IF (NUM /= 1) THEN                             ! Sort row by the shell method so that RJ is in numerical order
               CALL SORT_INT1_REAL1 ( SUBR_NAME, 'RJ, REMS', NUM, RJ, REMS )
            ENDIF

j_do3:      DO J = 1,NUM
               KTERM_MGGE = KTERM_MGGE + 1                 ! KTERM_MGGE is a count on the number of records
               IF (KTERM_MGGE > NTERM_MGGE) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, NTERM_MGGE, 'MGGE' )
                  J_MGGE(KTERM_MGGE) = RJ(J)
                    MGGE(KTERM_MGGE) = REMS(J)
            ENDDO j_do3

         ENDDO k_do
         CALL COUNTER_PROGRESS(I)
      ENDDO i_do

      WRITE(SC1,*) CR13


      IF (KTERM_MGGE /= NTERM_MGGE) THEN                   ! Check KTERM_MGGE = NTERM_MGGE
         WRITE(ERR,1614) SUBR_NAME,LINK1R,KTERM_MGGE,NTERM_MGGE
         WRITE(F06,1614) SUBR_NAME,LINK1R,KTERM_MGGE,NTERM_MGGE
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
      ENDIF

! *********************************************************************************************************************************
! Call subr to calc MGGS matrix of scalar masses

      IF (NCMASS > 0) THEN
         CALL MGGS_MASS_MATRIX
      ENDIF

! Add MGGC, MGGE and MGGS to get MGG. This is done in 2 steps: add MGGC and MGGE to get temporary CRS1 then add CRS1 to MGGS

!  (1) add MGGC and MGGE to get CRS1 (do not mult by WTMASS here)
!  --------------------------------------------------------------

      CALL MATADD_SSS_NTERM ( NDOFG, 'MGGC', NTERM_MGGC, I_MGGC, J_MGGC, SYM_MGGC, 'MGGE', NTERM_MGGE, I_MGGE, J_MGGE, SYM_MGGE,&
                                     'CRS1' , NTERM_CRS1 )

      CALL ALLOCATE_SCR_CRS_MAT ( 'CRS1', NDOFG, NTERM_CRS1, SUBR_NAME )

      CALL MATADD_SSS ( NDOFG, 'MGGC', NTERM_MGGC, I_MGGC, J_MGGC, MGGC, ONE, 'MGGE', NTERM_MGGE, I_MGGE, J_MGGE, MGGE,            &
                        ONE, 'CRS1', NTERM_CRS1, I_CRS1, J_CRS1, CRS1 )


!  (2) add CRS1 = MGGC + MGGE and MGGS to get CMGG (mult by WTMASS here)
!  ---------------------------------------------------------------------

      IF (NTERM_MGGS > 0) THEN

         SYM_CRS1 = SYM_MGGS
         CALL MATADD_SSS_NTERM ( NDOFG, 'CRS1', NTERM_CRS1, I_CRS1, J_CRS1, SYM_CRS1, 'MGGS', NTERM_MGGS, I_MGGS, J_MGGS, SYM_MGGS,&
                                        'MGG' , NTERM_MGG )

         CALL ALLOCATE_L1_MGG ( 'I2_MGG', SUBR_NAME )
         CALL ALLOCATE_SPARSE_MAT ( 'MGG', NDOFG, NDOFG, SUBR_NAME )

         CALL MATADD_SSS ( NDOFG, 'CRS1', NTERM_CRS1, I_CRS1, J_CRS1, CRS1, WTMASS, 'MGGS', NTERM_MGGS, I_MGGS, J_MGGS, MGGS,      &
                           WTMASS, 'MGG', NTERM_MGG, I_MGG, J_MGG, MGG )

      ELSE

         NTERM_MGG = NTERM_CRS1
         CALL ALLOCATE_L1_MGG ( 'I2_MGG', SUBR_NAME )
         CALL ALLOCATE_SPARSE_MAT ( 'MGG', NDOFG, NTERM_MGG, SUBR_NAME )
         DO I=1,NDOFG+1
            I_MGG(I) = I_CRS1(I)
         ENDDO
         DO I=1,NTERM_MGG
            J_MGG(I) = J_CRS1(I)
              MGG(I) = WTMASS*CRS1(I)
         ENDDO


      ENDIF

! Deallocate CRS1

      CALL DEALLOCATE_SCR_MAT ( 'CRS1' )

      IF (PRTMASS(1) >= 2) THEN
         IF (NTERM_MGG > 0) THEN
            IF (ALLOCATED(MGGC)) THEN
               CALL WRITE_SPARSE_CRS (  'Conc mass matrix, MGGC', 'G ', 'G ', NTERM_MGGC, NDOFG, I_MGGC, J_MGGC, MGGC )
            ENDIF
            IF (ALLOCATED(MGGE)) THEN
               CALL WRITE_SPARSE_CRS (  'Elem mass matrix, MGGE', 'G ', 'G ', NTERM_MGGE, NDOFG, I_MGGE, J_MGGE, MGGE )
            ENDIF
            IF (ALLOCATED(MGGS)) THEN
               CALL WRITE_SPARSE_CRS ('Scalar mass matrix, MGGS', 'G ', 'G ', NTERM_MGGS, NDOFG, I_MGGS, J_MGGS, MGGS )
            ENDIF
         ENDIF
      ENDIF

! *********************************************************************************************************************************
! Write row, col, value to L1R for matrix MGG

      K = 0
      WRITE(L1R) NTERM_MGG
      IF (NTERM_MGG > 0) THEN
         DO I=1,NDOFG
            NUM_IN_ROW_I = I_MGG(I+1) - I_MGG(I)
            DO J=1,NUM_IN_ROW_I
               K = K + 1
               IF (K > NTERM_MGG)  CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, K, 'MGG' )
               I2_MGG(K) = I
               WRITE(L1R) I2_MGG(K), J_MGG(K), MGG(K)
            ENDDO
         ENDDO
      ENDIF

      CALL FILE_CLOSE ( L1R, LINK1R, 'KEEP' )

! Get stats on MGG to write to F06

      IF (NTERM_MGG > 0) THEN

         MAX_NUM_IN_ROW = 0
         DO I=1,NDOFG
            IK = I_MGG(I+1) - I_MGG(I)
            IF (IK > MAX_NUM_IN_ROW) THEN
               MAX_NUM_IN_ROW = IK
            ENDIF
         ENDDO

         WRITE(ERR,147) NTERM_MGG
         WRITE(ERR,101) MAX_NUM_IN_ROW
         IF (SUPINFO == 'N') THEN
            WRITE(F06,147) NTERM_MGG
            WRITE(F06,101) MAX_NUM_IN_ROW
         ENDIF

      ENDIF

! Debug output (print grid 6x6 mass for every grid)

      IERR = 0
      IF (DEBUG(36) > 0) THEN
         WRITE(F06,1101)
         DO K=1,NGRID
            CALL GET_GRID_NUM_COMPS ( K, NUM_COMPS, SUBR_NAME )
            IF (NUM_COMPS == 6) THEN                       ! Only do output for actual grids, not SPOINT's
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_ID(K), IGRID )
               IF (IGRID == -1) THEN
                  IERR      = IERR + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1314) 'GRID ', GRID_ID(K), ' CANNOT PROCESS GRID FOR DEBUG(36) OUTPUT'
                  WRITE(F06,1314) 'GRID ', GRID_ID(K), ' CANNOT PROCESS GRID FOR DEBUG(36) OUTPUT'
               ENDIF
               CALL GET_GRID_6X6_MASS (  GRID_ID(K), IGRID, FOUND, GRID_MGG )
               WRITE(F06,1102) GRID_ID(K)
               DO I=1,3
                  WRITE(F06,1103) (GRID_MGG(I,J),J=1,6)
               ENDDO
               WRITE(F06,*)
               DO I=4,6
                  WRITE(F06,1103) (GRID_MGG(I,J),J=1,6)
               ENDDO
               WRITE(F06,*)
            ENDIF
         ENDDO
         WRITE(F06,1104)
      ENDIF



      RETURN

! **********************************************************************************************************************************
  146 FORMAT(' *INFORMATION: NUMBER OF NONZERO TERMS IN THE MGGE MASS MATRIX (ELEMS) IS             = ',I12,/)

  147 FORMAT(' *INFORMATION: NUMBER OF NONZERO TERMS IN THE MGG MASS MATRIX (ELEMS + CONM) IS       = ',I12,/)

  101 FORMAT(' *INFORMATION: MAX NUMBER OF NONZERO TERMS IN A ROW OF THE G-SET MASS MATRIX          = ',I12,/)

 1101 FORMAT(' ___________________________________________________________________________________________________________________'&
            ,'________________'                                                                                                ,//,&
             ' ::::::::::::::::::::::::::::::::::::::::START DEBUG(36) OUTPUT FROM SUBROUTINE SPARSE_MGG:::::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1102 FORMAT('6 x 6 mass matrix for grid ',I8,/,'-----------------------------------')

 1103 FORMAT(3(1ES14.6),2X,3(1ES14.6))

 1104 FORMAT(' :::::::::::::::::::::::::::::::::::::::::END DEBUG(36) OUTPUT FROM SUBROUTINE SPARSE_MGG::::::::::::::::::::::::::',&
             ':::::::::::::::::'                                                                                                ,/,&
             ' ___________________________________________________________________________________________________________________'&
            ,'________________',/)

 1314 FORMAT(' *ERROR  1314: UNDEFINED ',A,I8,A)

 1626 FORMAT(' *ERROR  1626: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' EMSPNT ARRAY INDICATES THERE IS MORE DATA IN ARRAY EMS FOR ROW ',I12,' OF THE MGG STIFF MATRIX.'      &
                    ,/,14X,' ALTHOUGH THE DOF COUNT IS AT THE END OF THE ROW')

 1614 FORMAT(' *ERROR  1614: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF G-SET MASS MATRIX RECORDS WRITTEN TO FILE:'                                             &
                    ,/,15X,A                                                                                                       &
                    ,/,14X,' WAS KTERM_MGGE = ',I12,'. IT SHOULD HAVE BEEN NTERM_MGGE = ',I12)







! **********************************************************************************************************************************

      END SUBROUTINE SPARSE_MGG

   END MODULE MASS_MATRIX_ASSEMBLY
