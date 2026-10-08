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

   MODULE STIFFNESS_MATRIX_ASSEMBLY

   USE ELEMENT_TRANSFORMATIONS, ONLY:  ELEM_TRANSFORM_LBG

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ESP0, ESP, SPARSE_KGG, SPARSE_KGGD

   CONTAINS

      SUBROUTINE ESP0

! Provides an estimate of the size that is required for array KGG or KGGD. The estimate is needed so that the G-set stiffness
! matrix can be allocated. There are 4 possible means of providing the estimate of this size, and they all use Bulk Data
! PARAM SETLKTK. The variable used here to estimate the size of either the KGG or KGGD stiffness matrix is LTERM. After all
! processing in subr ESP0, LTERM will be set to be either LTERM_KGG or LTERM_LGGD, depending on which matrix is being processed

! Bulk Data PARAM card SETLKTK is used to control how LTERM is estimated. Field 3 of the B.D. PARAM SETLKTK card can
! be either 0, 1, 2 or 3:

! If field 3 of PARAM SETLKTK is 0: then the estimate of LTERM is based
!    on full elem KE matrices unconnected (i.e. connected DOF's recounted). If the value in field 3 of the PARAM SETLKTK
!    card is PAUSE then LINK1 is PAUSE'd after subr ESP0 so user can change the estimate of LTERM

! If field 3 of PARAM SETLKTK is 1: then the estimate of LTERM is based
!    on Bandit bandwidth of KGG times NDOFG. If the value in field 4 of the B.D. PARAM SETLKTK card is PAUSE then LINK1 is
!    PAUSE'd after subr ESP0 so user can change the estimate of LTERM

! If field 3 of PARAM SETLKTK is 2: then the estimate of LTERM is based on actual elem KE matrices unconnected.

! If field 3 of PARAM SETLKTK is 3: then the value in field 4 is used as the estimate for LTERM.


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, KMAT_BW, KMAT_DEN, LTERM_KGG, LTERM_KGGD, SOL_NAME
      USE PARAMS, ONLY                :  GRIDSEQ, SETLKTK, SUPINFO, USR_LTERM_KGG
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO

      USE ELEMENT_LOOKUPS, ONLY       :  GET_ELGP
      USE EMG_MOD, ONLY               :  EMG
      USE ELEMENT_TRANSFORMATIONS, ONLY:  ELEM_TRANSFORM_LBG
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ESP0'

      INTEGER(LONG)                   :: LTERM             ! Count of number of estimated terms in KGG or KGGD




! **********************************************************************************************************************************
      IF      (SETLKTK == 0) THEN                          ! LTERM based on full elem stiffness matrices unconnected

          CALL ESP0_0 ( LTERM )

      ELSE IF (SETLKTK == 1) THEN                          ! LTERM based on BW returned from subr BANDIT (if BANDIT was run).

         IF (KMAT_BW > 0) THEN
            CALL ESP0_1 ( LTERM )
         ELSE
            CALL ESP0_0 ( LTERM )
         ENDIF

      ELSE IF (SETLKTK == 2) THEN                          ! LTERM based on matrix density returned from subr BANDIT

         IF (KMAT_DEN > ZERO) THEN
            CALL ESP0_2 ( LTERM )
         ELSE
            CALL ESP0_0 ( LTERM )
         ENDIF

      ELSE IF (SETLKTK == 3) THEN                          ! LTERM based on actual elem KE matrices unconnected

         CALL ESP0_3 ( LTERM )

      ELSE IF (SETLKTK == 4) THEN                          ! Use estimate from user

         LTERM = USR_LTERM_KGG

      ENDIF

! Now use LTERM to be the estimate for the appropriate G-set stiffness:

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         LTERM_KGGD = LTERM
      ELSE
         LTERM_KGG  = LTERM
      ENDIF



      RETURN

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE ESP0_0 ( LTERM )

! Estimates LTERM based on full elem stiffness matrices in an unassembled state (not connected)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NELE, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  SPARSTOR
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  EDAT, EPNT, ETYPE, ELGP, TYPE
      use model_stuf, only            :  eid

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ESP0_0'

      INTEGER(LONG), INTENT(OUT)      :: LTERM             ! Count of number of estimated terms in KGG or KGGD
      INTEGER(LONG)                   :: DELTA_LTERM       ! Increment of LTERM for one element
      INTEGER(LONG)                   :: I                 ! DO loop index





! **********************************************************************************************************************************
! Process the elements: Asume each is has a stiffness matrix that is completely full


      LTERM = 0
      DO I = 1,NELE

         CALL GET_ELGP ( I )

         IF (SPARSTOR == 'SYM   ') THEN
            DELTA_LTERM = 3*ELGP*(6*ELGP + 1)
            LTERM = LTERM + DELTA_LTERM        ! SYM has terms only on diag and above
         ELSE
            DELTA_LTERM = (6*ELGP)*(6*ELGP)
            LTERM = LTERM + DELTA_LTERM        ! NONSYM can have full matrix
         ENDIF


      ENDDO

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         WRITE(ERR,4321) LTERM, SETLKTK
         IF (SUPINFO == 'N') THEN
            WRITE(F06,4321) LTERM, SETLKTK
         ENDIF
      ELSE
         WRITE(ERR,4322) LTERM, SETLKTK
         IF (SUPINFO == 'N') THEN
            WRITE(F06,4322) LTERM, SETLKTK
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 4321 FORMAT(' *INFORMATION: IN ESPO_0: ESTIMATE OF NUMBER OF NONZEROS IN STIFF MATRIX KGGD IS      = ',I12,                       &
                           ' BASED ON PARAM SETLKTK = ',I3)

 4322 FORMAT(' *INFORMATION: IN ESPO_0: ESTIMATE OF NUMBER OF NONZEROS IN STIFF MATRIX KGG IS       = ',I12,                       &
                           ' BASED ON PARAM SETLKTK = ',I3)

! **********************************************************************************************************************************

      END SUBROUTINE ESP0_0

! ##################################################################################################################################

      SUBROUTINE ESP0_1 ( LTERM )

! Estimates LTERM based on number of rows in the stiff matrix times the stiffness matrix bandwidth from BANDIT.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  KMAT_BW, KMAT_DEN, NDOFG, BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ESP0_1'

      INTEGER(LONG), INTENT(OUT)      :: LTERM             ! Count of number of estimated terms in KGG or KGGD





! **********************************************************************************************************************************
! Estimate number of nonzero terms as the number of rows in the stiff matrix times the stiff matrix bandwidth:

      LTERM = NDOFG*KMAT_BW

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         WRITE(ERR,4321) LTERM, SETLKTK, NDOFG, KMAT_BW
         IF (SUPINFO == 'N') THEN
            WRITE(F06,4321) LTERM, SETLKTK, NDOFG, KMAT_BW
         ENDIF
      ELSE
         WRITE(ERR,4322) LTERM, SETLKTK, NDOFG, KMAT_BW
         IF (SUPINFO == 'N') THEN
            WRITE(F06,4322) LTERM, SETLKTK, NDOFG, KMAT_BW
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 4321 FORMAT(' *INFORMATION: IN ESP0_1: ESTIMATE OF NUMBER OF NONZEROS IN STIFF MATRIX KGGD IS      = ',I12,                       &
                           ' BASED ON PARAM SETLKTK = ',I3                                                                         &
                   ,/,100X,' NDOFG         = ',I12,/,100X,' BANDIT BW     = ',I12)

 4322 FORMAT(' *INFORMATION: IN ESP0_1: ESTIMATE OF NUMBER OF NONZEROS IN STIFF MATRIX KGG IS       = ',I12,                       &
                           ' BASED ON PARAM SETLKTK = ',I3                                                                         &
                   ,/,100X,' NDOFG         = ',I12,/,100X,' BANDIT BW     = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE ESP0_1

! ##################################################################################################################################

      SUBROUTINE ESP0_2 ( LTERM )

! Estimates LTERM based on the full size of the stiffness matrix times the density returned from subr BANDIT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, KMAT_BW, KMAT_DEN, NDOFG, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE_HUNDRED
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ESP0_2'

      INTEGER(LONG), INTENT(OUT)      :: LTERM             ! Count of number of estimated terms in KGG or KGGD





! **********************************************************************************************************************************
! Estimate number of nonzero terms as the number of rows in the stiff matrix times the stiff matrix bandwidth:

      LTERM = NINT((KMAT_DEN/ONE_HUNDRED)*NDOFG*NDOFG)  ! KMAT_DEN is in %

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         WRITE(ERR,4321) LTERM, SETLKTK, NDOFG, KMAT_DEN
         IF (SUPINFO == 'N') THEN
            WRITE(F06,4321) LTERM, SETLKTK, NDOFG, KMAT_DEN
         ENDIF
      ELSE
         WRITE(ERR,4322) LTERM, SETLKTK, NDOFG, KMAT_DEN
         IF (SUPINFO == 'N') THEN
            WRITE(F06,4322) LTERM, SETLKTK, NDOFG, KMAT_DEN
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 4321 FORMAT(' *INFORMATION: IN ESP0_2: ESTIMATE OF NUMBER OF NONZEROS IN STIFF MATRIX KGGD IS      = ',I12,                       &
                           ' BASED ON PARAM SETLKTK = ',I3                                                                         &
                   ,/,100X,' NDOFG         = ',I12,/,100X,' MATRIX DENSITY= ',1ES11.3,'%')

 4322 FORMAT(' *INFORMATION: IN ESP0_2: ESTIMATE OF NUMBER OF NONZEROS IN STIFF MATRIX KGG IS       = ',I12,                       &
                           ' BASED ON PARAM SETLKTK = ',I3                                                                         &
                   ,/,100X,' NDOFG         = ',I12,/,100X,' MATRIX DENSITY= ',1ES11.3,'%')

! **********************************************************************************************************************************

      END SUBROUTINE ESP0_2

! ##################################################################################################################################

      SUBROUTINE ESP0_3 ( LTERM )

! Estimates LTERM based on actual element stiffness matrices unconnected.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MELDOF, NELE, NSUB, SOL_NAME
      USE PARAMS, ONLY                :  EPSIL, SETLKTK, SPARSTOR
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  ELDOF, NUM_EMG_FATAL_ERRS, PLY_NUM, KE, TYPE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ESP0_3'
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option flags for subr EMG (to tell it what to calc)

      INTEGER(LONG), INTENT(OUT)      :: LTERM             ! Count of number of estimated terms in KGG or KGGD
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

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE

      IF      ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         OPT(4) = 'N'                                      ! OPT(4) is for calc of KE-linear
         OPT(6) = 'Y'                                      ! OPT(6) is for calc of KE-nonlinear
      ELSE IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
         OPT(4) = 'Y'                                      ! OPT(4) is for calc of KE-linear
         OPT(6) = 'Y'                                      ! OPT(6) is for calc of KE-nonlinear
      ELSE
         OPT(4) = 'Y'                                      ! OPT(4) is for calc of KE-linear
         OPT(6) = 'N'                                      ! OPT(6) is for calc of KE-nonlinear
      ENDIF

! Process the elements:

      IERROR = 0
      LTERM  = 0
      CALL COUNTER_INIT('Estimate size of KGG: process elem  ', NELE)
elems:DO I=1,NELE

         PLY_NUM = 0
         CALL EMG ( I   , OPT, 'N', SUBR_NAME, 'N' )       ! 'N' means do not write to BUG file

         IF (NUM_EMG_FATAL_ERRS /=0) THEN
            IERROR = IERROR + NUM_EMG_FATAL_ERRS
            CYCLE elems
         ENDIF

! Transform KE from local at the elem ends to basic at elem ends to global at elem ends to global at grids.
                                                           ! Transform PTE from local-basic-global
         IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  '))THEN
            CALL ELEM_TRANSFORM_LBG ( 'KE', KE, DQE )
         ENDIF

! Count nonzero terms in transformed KE

kgg_rows:DO J=1,ELDOF
            IF (SPARSTOR == 'SYM') THEN
               KSTART = J
            ELSE
               KSTART = 1
            ENDIF
kgg_cols:   DO K=KSTART,ELDOF
               IF (DABS(KE(J,K)) < EPS1) THEN
                  CYCLE kgg_cols
               ELSE
                  LTERM = LTERM + 1
               ENDIF
            ENDDO kgg_cols
         ENDDO kgg_rows
         CALL COUNTER_PROGRESS(I)
      ENDDO elems
      WRITE(SC1,*) CR13

! Quit if IERROR > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,9876) IERROR
         WRITE(F06,9876) IERROR
         CALL OUTA_HERE ( 'Y' )                            ! IERROR is count of all subr EMG errors, so quit
      ENDIF

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         WRITE(ERR,4321) LTERM, SETLKTK
         IF (SUPINFO == 'N') THEN
            WRITE(F06,4321) LTERM, SETLKTK
         ENDIF
      ELSE
         WRITE(ERR,4322) LTERM, SETLKTK
         IF (SUPINFO == 'N') THEN
            WRITE(F06,4322) LTERM, SETLKTK
         ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 4321 FORMAT(' *INFORMATION: IN ESP0_3: INITIAL EST. OF NUMBER OF NONZEROS IN STIFF MATRIX KGGD IS  = ',I12,                       &
                           ' BASED ON PARAM SETLKTK = ',I3)

 4322 FORMAT(' *INFORMATION: IN ESP0_3: INITIAL EST. OF NUMBER OF NONZEROS IN STIFF MATRIX KGG IS   = ',I12,                       &
                           ' BASED ON PARAM SETLKTK = ',I3)

 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')

! **********************************************************************************************************************************

      END SUBROUTINE ESP0_3

! ##################################################################################################################################

      SUBROUTINE DUMPSTF0 ( WHAT, J, K, KGG_ROW, KGG_COL )

! Prints out info on the formulation of stiffness arrays for subr ESP0_3 which estimates LTERM for subr ESP

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE MODEL_STUF, ONLY            :  EID

      IMPLICIT NONE

      CHARACTER(1*BYTE), INTENT(IN)   :: WHAT              ! Indicator of where this subr was called from in subr ESP0_3

      INTEGER(LONG)    , INTENT(IN)   :: J                 ! Row number of elem stiff matrix term, KE(J,K)
      INTEGER(LONG)    , INTENT(IN)   :: K                 ! Col number of elem stiff matrix term, KE(J,K)
      INTEGER(LONG)    , INTENT(IN)   :: KGG_COL           ! Row number of KGG matrix where KE(J,K) goes
      INTEGER(LONG)    , INTENT(IN)   :: KGG_ROW           ! Col number of KGG matrix where KE(J,K) goes

! **********************************************************************************************************************************
      IF      (WHAT == '0') THEN

         WRITE(F06,8910)

      ELSE IF (WHAT == 'A') THEN

         WRITE(F06,8930) EID, J, K, KGG_ROW, KGG_COL

      ENDIF

      RETURN

! **********************************************************************************************************************************
 8910 FORMAT(1X,'     ELEM        J        K  KGG_ROW  KGG_COL')

 8930 FORMAT(1X,'A',I8,I8,I8,I8,I8)

! **********************************************************************************************************************************

      END SUBROUTINE DUMPSTF0

      END SUBROUTINE ESP0


      SUBROUTINE ESP0_FINAL

! Estimate number of terms in stiffness matrix by going through the complete process of generating it except that terms are not put
! into array STF. This way, we only need to allocate arrays STFCOL, STFPNT with the conservative estimate of LTERM_KGG. Once this
! subr finishes, we will hae an exact value of LTERM_KGG at which time we deallocate STFCOL, STFPNT and allocate them plus STF with
! the exact LTERM_KGG and then calculate them in subr ESP.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, IBIT, LTERM_KGG, MELDOF, NELE, NGRID, NTERM_KGG, NSUB
      USE PARAMS, ONLY                :  EPSIL, SPARSTOR, SUPINFO
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE MODEL_STUF, ONLY            :  AGRID, ELDT, ELDOF, ELGP, GRID_ID, NUM_EMG_FATAL_ERRS, PLY_NUM, KE, TYPE
      USE STF_ARRAYS, ONLY            :  STFKEY, STF3
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE ELEMENT_TRANSFORMATIONS, ONLY:  ELEM_TRANSFORM_LBG

      USE EMG_MOD, ONLY               :  EMG
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ESP0_FINAL'
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option flags for subr EMG (to tell it what to calc)

      INTEGER(LONG)                   :: EDOF(MELDOF)      ! A list of the G-set DOF's for an elem
      INTEGER(LONG)                   :: EDOF_ROW_NUM      ! Row number in array EDOF
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no. in array TDOF where G-set DOF's are kept
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IDUM              ! Dummy variable used when flipping DOF's
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IS                ! A pointer into array STF3
      INTEGER(LONG)                   :: ISS               ! A particular value of IS
      INTEGER(LONG)                   :: KGG_ROW           ! A row no. in KGG
      INTEGER(LONG)                   :: KGG_ROWJ          ! Another row no. in KGG
      INTEGER(LONG)                   :: KGG_COL           ! A col no. in KGG
      INTEGER(LONG)                   :: KSTART            ! Used in deciding whether to process all elem stiffness terms or only
!                                                            the ones on and above the diagonal (controlled by param SPARSTOR)
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID is a physical grid, 1 if a scalar point
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: TDOF_ROW_NUM      ! Row number in array TDOF
                                                           ! Indicator for output of elem data to BUG file


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

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'Y'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

! Process the elements:

      IS  = 0
      ISS = IS
      LTERM_KGG = 0
      CALL COUNTER_INIT('Estimate size of KGG: process elem  ', NELE)
elems:DO I=1,NELE

         PLY_NUM = 0
         CALL EMG ( I   , OPT, 'N', SUBR_NAME, 'N' )       ! 'N' means do not write to BUG file

         EDOF_ROW_NUM = 0                                  ! Generate element DOF'S
         DO J = 1,ELGP
!xx         CALL CALC_TDOF_ROW_NUM ( AGRID(J), ROW_NUM_START, 'N' )
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

! Transform KE from local at the elem ends to basic at elem ends to global at elem ends to global at grids.
                                                           ! Transform PTE from local-basic-global
         IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  '))THEN
            CALL ELEM_TRANSFORM_LBG ( 'KE', KE, DQE )
         ENDIF

! Put the element stiff matrix, KE (now in global coords), into STF array. J ranges over rows, K over cols of elem stiff matrix, KE

kgg_rows:DO J = 1,ELDOF
            KGG_ROWJ  = EDOF(J)

            IF (SPARSTOR == 'SYM') THEN                    ! Set KSTART depending on SPARSTOR
               KSTART = J                                  ! Process only upper right portion of ME
            ELSE
               KSTART = 1                                  ! Process all of ME
            ENDIF

kgg_cols:   DO K = KSTART,ELDOF
               KGG_ROW  = KGG_ROWJ                         ! Make sure we have correct row num. It may have been flipped w/ col
               KGG_COL  = EDOF(K)
               IF (DABS(KE(J,K)) < EPS1) THEN
                  CYCLE kgg_cols
               ENDIF

               IF (SPARSTOR == 'SYM') THEN                 ! If 'SYM', Flip KGG_COL,KGG_ROW if KGG_COL < KGG_ROW
                  IF (KGG_COL < KGG_ROW) THEN
                     IDUM    = KGG_ROW
                     KGG_ROW = KGG_COL
                     KGG_COL = IDUM
                  ENDIF
               ENDIF

               IS = STFKEY(KGG_ROW)                        ! Get pointer to first term in row KGG_ROW of global stiff matrix

               IF (IS == 0) THEN                           ! STFKEY(KGG_ROW)=0 means no current terms in global stiff matrix at row
                                                           ! KGG_ROW update NTERM_KGG and reset STFKEY, STFCOL, STFPNT, STF arrays
                  LTERM_KGG = LTERM_KGG + 1

                  STFKEY(KGG_ROW)       = LTERM_KGG
                  STF3(LTERM_KGG)%Col_1 = KGG_COL
                  STF3(LTERM_KGG)%Col_2 = 0

               ELSE                                        ! STFKEY(KGG_ROW) /= 0 means there are already some terms in row KGG_ROW

stfpnt0:          DO                                       ! so, run this loop until we find a place to put KE(J,K). If there is
                                                           ! already a term in this row w/ same DOF's as KE(J,K), loop runs once.
                                                           ! If not, then this loop runs until it finds STFPNT=0, and inserts term.

                     IF (KGG_COL == STF3(IS)%Col_1) THEN   ! There is a term that exists with same DOF's as KE(J,K) so add terms

                        CYCLE kgg_cols                     ! We have added a term to STF so exit this loop and go to next col of KGG

                     ELSE                                  ! This is a new term for row J. Need to cycle until we find STFPNT = 0.
                                                           ! Then we can put KE(J,K) in STF
                        ISS = IS
                        IS  = STF3(IS)%Col_2
                        IF (IS == 0) THEN                  ! We are at end of where terms are in this row, so KE(J,K) goes here
                           LTERM_KGG         = LTERM_KGG+1 ! Increment LTERM_KGG
                           STF3(ISS)%Col_2       = LTERM_KGG! STFPNT for the current KE(J,K) term
                           STF3(LTERM_KGG)%Col_2 = 0       ! Latest STFPNT is set to 0 so we will know when to insert next KE(J,K)
                           STF3(LTERM_KGG)%Col_1 = KGG_COL ! STFCOL always is KGG_COL

                           CYCLE kgg_cols                  ! We put KE(J,K) into STF so exit this loop and go to next col of KGG
                        ELSE                               ! STFPNT /= 0 so cycle this loop until we get it = 0
                           CYCLE stfpnt0
                        ENDIF

                     ENDIF

                  ENDDO stfpnt0

               ENDIF

            ENDDO kgg_cols

         ENDDO kgg_rows
         CALL COUNTER_PROGRESS(I)
      ENDDO elems
      WRITE(SC1,*) CR13

! Reset subr EMG option flags:

      OPT(3) = 'N'
      OPT(4) = 'N'

      WRITE(ERR,4321) LTERM_KGG
      IF (SUPINFO == 'N') THEN
         WRITE(F06,4321) LTERM_KGG
      ENDIF



      RETURN

! **********************************************************************************************************************************
 4321 FORMAT(' *INFORMATION: FINAL   LTERM_KGG EST OF THE NUMBER OF NONZEROS IN STIFF MATRIX KGG IS = ',I12)

! **********************************************************************************************************************************

      END SUBROUTINE ESP0_FINAL


      SUBROUTINE ESP

! Element stiffness processor

! ESP generates the G-set stiffness matrix and puts it into the 1D array STF of nonzero stiffness terms above the
! diagonal.

! ESP processes the elements sequentially to generate element KE matrix using the EMG set of routines. The element
! stiffness are transformed from local to basic to global coords for each grid and then merged into the system
! stiffness, STF, array. See explanation, with an example, in module STF_ARRAYS


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, F23, F23FIL, F23_MSG, F24, F24FIL, F24_MSG, FILE_NAM_MAXLEN, SC1, SCR,     &
                                         WRT_BUG, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_KE_BIT, ELDT_BUG_SE_BIT,                                           &
                                         ELDT_F23_KE_BIT, ELDT_F24_SE_BIT, ELDT_BUG_BCHK_BIT, ELDT_BUG_BMAT_BIT, ELDT_BUG_SHPJ_BIT,&
                                         FATAL_ERR, IBIT, LINKNO, LTERM_KGG, LTERM_KGGD, MBUG, MELDOF, NDOFG, NELE, NGRID,         &
                                         NTERM_KGG, NTERM_KGGD, NSUB, SOL_NAME
      USE PARAMS, ONLY                :  EPSIL, SPARSTOR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  AGRID, ELDT, ELDOF, ELGP, GRID_ID, NUM_EMG_FATAL_ERRS, PLY_NUM, OELDT, KE, KED, TYPE
      USE STF_ARRAYS, ONLY            :  STFKEY, STF3
      USE DERIVED_DATA_TYPES, ONLY    :  INT2_REAL1
      USE STF_TEMPLATE_ARRAYS, ONLY   :  CROW, TEMPLATE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE LINK1_WORKSPACE_LIFECYCLE, ONLY:  ALLOCATE_STF_ARRAYS, ALLOCATE_TEMPLATE, DEALLOCATE_STF_ARRAYS, DEALLOCATE_TEMPLATE
      USE EMG_MOD, ONLY               :  EMG
      USE TEMP_FILE_WRITERS, ONLY     :  WRITE_FIJFIL
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE FILE_LIFECYCLE, ONLY        :  FILE_CLOSE, READERR
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      TYPE(INT2_REAL1), ALLOCATABLE   :: STF3_KEEP(:)      ! In-memory copy of STF3(1:NTERM) while STF3 is reallocated

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ESP'
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option flags for subr EMG (to tell it what to calc)
      CHARACTER(24*BYTE)              :: NAME              ! Name for output error purposes
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: SCRFIL            ! File name

      INTEGER(LONG), PARAMETER        :: DEB_NUM   = 46    ! Debug number for output error message
      INTEGER(LONG)                   :: EDOF(MELDOF)      ! A list of the G-set DOF's for an elem
      INTEGER(LONG)                   :: EDOF_ROW_NUM      ! Row number in array EDOF
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no. in array TDOF where G-set DOF's are kept
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: I1                ! Intermediate variable resulting from an IAND operation
      INTEGER(LONG)                   :: IDUM              ! Dummy variable used when flipping DOF's
      INTEGER(LONG)                   :: IERROR            ! Local error indicator
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening a file
      INTEGER(LONG)                   :: IS                ! A pointer into arrays STFKEY and STFPNT
      INTEGER(LONG)                   :: ISS               ! A particular value of IS
      INTEGER(LONG)                   :: KGG_ROW           ! A row no. in KGG or KGGD
      INTEGER(LONG)                   :: KGG_ROWJ          ! Another row no. in KGG or KGGD
      INTEGER(LONG)                   :: KGG_COL           ! A col no. in KGG or KGGD
      INTEGER(LONG)                   :: KSTART            ! Used in deciding whether to process all elem stiffness terms or only
!                                                            the ones on and above the diagonal (controlled by param SPARSTOR)
      INTEGER(LONG)                   :: MAX_NUM           ! MAX of NTERM_KGG/NDOFG (used for DEBUG printout)
      INTEGER(LONG)                   :: NTERM             ! Either NTERM_KGGD (BUCKLING) or NTERM_KGG otherwise
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID is a physical grid, 1 if a scalar point
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to.
      INTEGER(LONG)                   :: PKTERM            ! Count of the terms in TEMPLATE for nonzero stiffness terms
      INTEGER(LONG)                   :: REC_NO            ! Record number when reading a file
      INTEGER(LONG)                   :: ROW_NUM_START     ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: TDOF_ROW_NUM      ! Row number in array TDOF
                                                           ! Indicator for output of elem data to BUG file
      INTEGER(LONG)                   :: LTERM             ! Either LTERM_KGGD (BUCKLING) or LTERM_KGG otherwise


      REAL(DOUBLE)                    :: DQE(MELDOF,NSUB)  ! Dummy array in call to ELEM_TRANSFORM_LBG
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare real zero

      INTRINSIC                       :: DABS
      INTRINSIC                       :: IAND
      INTRINSIC                       :: MAX



! **********************************************************************************************************************************
      EPS1  = EPSIL(1)
      NTERM = 0

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

! Null dummy array DQE used in call to ELEM_TRANSFORM_LBG

      DO I=1,MELDOF
         DO J=1,NSUB
            DQE(I,J) = ZERO
         ENDDO
      ENDDO

! LTERM is used in this subr to make sure we do not try to use more array dimension than was allocated. So here, set LTERM to be
! either LTERM_KGG or LTERM_KGGD depending on BUCKLING

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         LTERM = LTERM_KGGD
      ELSE
         LTERM = LTERM_KGG
      ENDIF


! DEBUG(10) = 13 or 33 requests that array TEMPLATE be printed

      IF ((DEBUG(10) == 13) .OR. (DEBUG(10) == 33)) THEN
         CALL ALLOCATE_TEMPLATE ( SUBR_NAME )
      ENDIF

! Set up the option flags for EMG:

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'N'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE

      IF      ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         OPT(4) = 'N'                                      ! OPT(4) is for calc of KE-linear
         OPT(6) = 'Y'                                      ! OPT(6) is for calc of KE-nonlinear
      ELSE IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
         OPT(4) = 'Y'                                      ! OPT(4) is for calc of KE-linear
         OPT(6) = 'Y'                                      ! OPT(6) is for calc of KE-nonlinear
      ELSE
         OPT(4) = 'Y'                                      ! OPT(4) is for calc of KE-linear
         OPT(6) = 'N'                                      ! OPT(6) is for calc of KE-nonlinear
      ENDIF

! Process the elements:

      IS  = 0
      ISS = IS
      IF ((DEBUG(10) == 12) .OR. (DEBUG(10) == 13) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
         CALL DUMPSTF ( '0', 0, 0, 0, 0, 0, 0 )
      ENDIF

      IERROR = 0
!xx   WRITE(SC1, * )                                       ! Advance 1 line for screen messages
      CALL COUNTER_INIT('     Calculating stiff matrix. Process elem  ', NELE)
      elems:DO I=1,NELE

         IF ((DEBUG(10) == 12) .OR. (DEBUG(10) == 13) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
            WRITE(F06,14001)
         ENDIF

         WRT_BUG = 0

         IF (LINKNO == 1) THEN                             ! Only want element stiff matrix, KE, written to BUG file in LINK 1

            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_KE_BIT))       ! WRT_BUG(4): printed output of KE
            IF (I1 > 0) THEN
               WRT_BUG(4) = 1
            ENDIF

            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_SE_BIT))       ! WRT_BUG(5): printed output of SEi, STEi
            IF (I1 > 0) THEN
               WRT_BUG(5) = 1
            ENDIF

            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_SHPJ_BIT))     ! WRT_BUG(7): printed output of shape fcns and Jacobians for some elems
            IF (I1 > 0) THEN
               WRT_BUG(7) = 1
            ENDIF

            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_BMAT_BIT))     ! WRT_BUG(8): printed output of strain-displ matrices for some elements
            IF (I1 > 0) THEN
               WRT_BUG(8) = 1
            ENDIF

            I1 = IAND(ELDT(I),IBIT(ELDT_BUG_BCHK_BIT))     ! WRT_BUG(9): printed output of R.B., const strain checks for some elems
            IF (I1 > 0) THEN
               WRT_BUG(9) = 1
            ENDIF

         ENDIF

         OPT(3) = 'N'                                      ! OPT(3) is for calc of SEi, STEi
         IF ((WRT_BUG(4) == 1) .OR. (WRT_BUG(5) == 1)) THEN
            OPT(3) = 'Y'
         ENDIF

         PLY_NUM = 0
         CALL EMG ( I   , OPT, 'Y', SUBR_NAME, 'Y' )       ! 'N' means do not write to BUG file

         IF (NUM_EMG_FATAL_ERRS /=0) THEN
            IERROR = IERROR + NUM_EMG_FATAL_ERRS
            CYCLE elems
         ENDIF

         I1 = IAND(OELDT,IBIT(ELDT_F23_KE_BIT))           ! Do we need to write elem stiff matrices to F23 files?
         IF (I1 > 0) THEN
            CALL WRITE_FIJFIL ( 3, 0 )
         ENDIF

         I1 = IAND(OELDT,IBIT(ELDT_F24_SE_BIT))           ! Do we need to write elem stress recovery matrices to F24 files?
         IF (I1 > 0) THEN
            CALL WRITE_FIJFIL ( 4, 0 )
         ENDIF

         EDOF_ROW_NUM = 0                                  ! Generate element DOF'S
         DO J = 1,ELGP
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

! Write diagonostics on negative diag stiffness before transformation to global

         IF ((DEBUG(189) == 1) .OR. (DEBUG(189) == 3)) THEN
            CALL WRITE_NEG_DIAG_STIFFNESS ( 1 )
         ENDIF

! Transform KE from local at the elem ends to basic at elem ends to global at elem ends to global at grids.
                                                           ! Transform PTE from local-basic-global
         IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  ')) THEN
            IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
               CALL ELEM_TRANSFORM_LBG ( 'KED', KED, DQE )
            ELSE
               CALL ELEM_TRANSFORM_LBG ( 'KE' , KE , DQE )
            ENDIF
         ENDIF

! Write diagonostics on negative diag stiffness after  transformation to global

         IF ((DEBUG(189) == 2) .OR. (DEBUG(189) == 3)) THEN
            CALL WRITE_NEG_DIAG_STIFFNESS ( 2 )
         ENDIF

! Put the elem stiff matrix, KE, or KED (now in global coords), into STF array. J ranges over rows, K over cols of elem stiff mat

kgg_rows:DO J = 1,ELDOF
            KGG_ROWJ  = EDOF(J)
            IF ((DEBUG(10) == 12) .OR. (DEBUG(10) == 13) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
               WRITE(F06,*)
            ENDIF

            IF (SPARSTOR == 'SYM') THEN                    ! Set KSTART depending on SPARSTOR
               KSTART = J                                  ! Process only upper right portion of ME
            ELSE
               KSTART = 1                                  ! Process all of ME
            ENDIF

kgg_cols:   DO K = KSTART,ELDOF
               KGG_ROW  = KGG_ROWJ                         ! Make sure we have correct row num. It may have been flipped w/ col
               KGG_COL  = EDOF(K)
               IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
                  IF (DABS(KED(J,K)) < EPS1) THEN
                     CYCLE kgg_cols
                  ENDIF
               ELSE
                  IF (DABS(KE(J,K)) < EPS1) THEN
                     CYCLE kgg_cols
                  ENDIF
               ENDIF

               IF (SPARSTOR == 'SYM') THEN                 ! If 'SYM', Flip KGG_COL,KGG_ROW if KGG_COL < KGG_ROW
                  IF (KGG_COL < KGG_ROW) THEN
                     IDUM    = KGG_ROW
                     KGG_ROW = KGG_COL
                     KGG_COL = IDUM
                  ENDIF
               ENDIF

               IS = STFKEY(KGG_ROW)                        ! Get pointer to first term in row KGG_ROW of global stiff matrix

               IF (IS == 0) THEN                           ! STFKEY(KGG_ROW)=0 means no current terms in global stiff matrix at row
                                                           ! KGG_ROW so update NTERM & reset STF arrays
                  NTERM = NTERM + 1

                  IF ((DEBUG(10) == 13) .OR. (DEBUG(10) == 33)) THEN
                     IF (ALLOCATED(TEMPLATE)) THEN
                        TEMPLATE(KGG_ROW,KGG_COL) = .TRUE.
                     ELSE
                        NAME = 'TEMPLATE                '
                        WRITE(ERR,1628) SUBR_NAME,DEB_NUM,NAME
                        WRITE(F06,1628) SUBR_NAME,DEB_NUM,NAME
                        FATAL_ERR = FATAL_ERR + 1
                        CALL OUTA_HERE ( 'Y' )             ! Coding error (TEMPLATE should be allocated), so quit
                     ENDIF
                  ENDIF

                  IF (NTERM > LTERM) THEN
                     WRITE(ERR,1624) SUBR_NAME, 'STIFFNESS','LTERM', LTERM
                     WRITE(F06,1624) SUBR_NAME, 'STIFFNESS','LTERM', LTERM
                     FATAL_ERR = FATAL_ERR + 1
                     CALL OUTA_HERE ( 'Y' )                ! MYSTRAN limitation, so quit
                  ENDIF

                  STFKEY(KGG_ROW)   = NTERM
                  STF3(NTERM)%Col_1 = KGG_COL
                  STF3(NTERM)%Col_2 = 0
                  IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
                     STF3(NTERM)%Col_3 = KED(J,K)
                  ELSE
                     STF3(NTERM)%Col_3 = KE(J,K)
                  ENDIF

                  IF ((DEBUG(10) == 12) .OR. (DEBUG(10) == 13) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
                     CALL DUMPSTF ( 'A', J, K, KGG_ROW, KGG_COL, IS, ISS )
                  ENDIF

               ELSE                                        ! STFKEY(KGG_ROW) /= 0 means there are already some terms in row KGG_ROW

stfpnt0:          DO                                       ! so, run this loop until we find a place to put stiff(J,K). If there is
                                                           ! already a term in this row w/ same DOF's as stiff(J,K), loop runs once.
                                                           ! If not, then this loop runs until it finds STFPNT=0, and inserts term.

                     IF (KGG_COL == STF3(IS)%Col_1) THEN   ! There is a term that exists with same DOF'S as stiff(J,K) so add terms

                        IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
                           STF3(IS)%Col_3 = STF3(IS)%Col_3 + KED(J,K)
                        ELSE
                           STF3(IS)%Col_3 = STF3(IS)%Col_3 + KE(J,K)
                        ENDIF
                        IF ((DEBUG(10) == 13) .OR. (DEBUG(10) == 33)) THEN
                           IF (ALLOCATED(TEMPLATE)) THEN
                              TEMPLATE(KGG_ROW,KGG_COL) = .TRUE.
                           ELSE
                              NAME = 'TEMPLATE                '
                              WRITE(ERR,1628) SUBR_NAME,DEB_NUM,NAME
                              WRITE(F06,1628) SUBR_NAME,DEB_NUM,NAME
                              FATAL_ERR = FATAL_ERR + 1
                              CALL OUTA_HERE ( 'Y' )       ! Coding error (TEMPLATE should be allocated), so quit
                           ENDIF
                        ENDIF
                        IF ((DEBUG(10) == 12) .OR. (DEBUG(10) == 13) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
                           CALL DUMPSTF ( 'B', J, K, KGG_ROW, KGG_COL, IS, ISS )
                        ENDIF

                        CYCLE kgg_cols                     ! We have added a term to STF so exit this loop and go to next col of KGG

                     ELSE                                  ! This is a new term for row J. Need to cycle until we find STFPNT = 0.
                                                           ! Then we can put KE(J,K), or KED(J,K) in STF
                        ISS = IS
                        IS  = STF3(IS)%Col_2
                        IF (IS == 0) THEN                  ! We are at end of where terms are in this row, so stiff(J,K) goes here
                           IF ((DEBUG(10) == 13) .OR. (DEBUG(10) == 33)) THEN
                              IF (ALLOCATED(TEMPLATE)) THEN
                                 TEMPLATE(KGG_ROW,KGG_COL) = .TRUE.
                              ELSE
                                 NAME = 'TEMPLATE                '
                                 WRITE(ERR,1628) SUBR_NAME,DEB_NUM,NAME
                                 WRITE(F06,1628) SUBR_NAME,DEB_NUM,NAME
                                 FATAL_ERR = FATAL_ERR + 1
                                 CALL OUTA_HERE ( 'Y' )    ! Coding error (TEMPLATE should be allocated), so quit
                              ENDIF
                           ENDIF
                           IF (NTERM+1 > LTERM) THEN
                              WRITE(ERR,1624) SUBR_NAME, 'STIFFNESS','LTERM', LTERM
                              WRITE(F06,1624) SUBR_NAME, 'STIFFNESS','LTERM', LTERM
                              FATAL_ERR = FATAL_ERR + 1
                              CALL OUTA_HERE ( 'Y' )       ! MYSTRAN limitation, so quit
                           ENDIF
                           NTERM = NTERM+1                 ! Increment NTERM
                           STF3(ISS)%Col_2   = NTERM       ! STFPNT for the current stiff(J,K) term
                           STF3(NTERM)%Col_2 = 0           ! Latest STFPNT is set to 0 so we will know when to insert next term
                           STF3(NTERM)%Col_1 = KGG_COL     ! STFCOL always is KGG_COL
                           IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
                              STF3(NTERM)%Col_3 = KED(J,K)
                           ELSE
                              STF3(NTERM)%Col_3 = KE(J,K)
                           ENDIF
                           IF ((DEBUG(10) == 12) .OR. (DEBUG(10) == 13) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
                              CALL DUMPSTF ( 'C', J, K, KGG_ROW, KGG_COL, IS, ISS )
                           ENDIF

                           CYCLE kgg_cols                  ! We put stiff(J,K) into STF so exit this loop and go to next col of KGG
                        ELSE                               ! STFPNT /= 0 so cycle this loop until we get it = 0
                           CYCLE stfpnt0
                        ENDIF

                     ENDIF

                  ENDDO stfpnt0

               ENDIF

            ENDDO kgg_cols

         ENDDO kgg_rows
         CALL COUNTER_PROGRESS(I)
      ENDDO elems
      WRITE(SC1,*) CR13

! Reset subr EMG option flags:

      OPT(3) = 'N'
      OPT(4) = 'N'
      OPT(6) = 'N'

! Quit if IERROR > 0

      IF (IERROR > 0) THEN
         WRITE(ERR,9876) IERROR
         WRITE(F06,9876) IERROR
         CALL OUTA_HERE ( 'Y' )                            ! IERROR is count of all subr EMG errors, so quit
      ENDIF

! Print out TEMPLATE which shows where the nonzero values are in the upper triangle of the stiffness matrix

      IF ((DEBUG(10) == 13) .OR. (DEBUG(10) == 33)) THEN

         IF (ALLOCATED(TEMPLATE)) THEN

            PKTERM = 0                                     ! Count nonzero terms in K based on TEMPLATE array.
            DO I=1,NDOFG                                   ! Call this PKTERM and it should be same as LTERM
               DO J=I,NDOFG
                  IF (TEMPLATE(I,J)) THEN
                     PKTERM = PKTERM + 1
                  ENDIF
               ENDDO
            ENDDO

            WRITE(F06,14002) PKTERM
            WRITE(F06,*)
            DO I=1,NDOFG
               DO J=1,NDOFG
                  CROW(J) = ' '
               ENDDO
               DO J=I,NDOFG
                  IF (TEMPLATE(I,J)) THEN
                     CROW(J) = 'K'
                  ELSE
                     CROW(J) = '_'
                  ENDIF
               ENDDO
               WRITE(F06,*) (CROW(J),J=1,NDOFG)
            ENDDO
            WRITE(F06,*)

         ELSE

            WRITE(ERR,1628) SUBR_NAME,DEB_NUM,NAME
            WRITE(F06,1628) SUBR_NAME,DEB_NUM,NAME
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                         ! Coding error (TEMPLATE should be allocated), so quit
            ENDIF

      ENDIF

! Deallocate TEMPLATE and CROW arrays

      IF ((ALLOCATED(TEMPLATE)) .OR. (ALLOCATED(CROW))) THEN
         CALL DEALLOCATE_TEMPLATE
      ENDIF

! Reallocate STF3 with the NTERM terms it holds, so that SPARSE_KGG does not carry the estimated size. The terms go through an
! in-memory copy. (They used to go through a scratch file, one record per term, which took seconds on large models; and LTERM was
! reset only after the reallocation, so STF3 came back at the estimated size.) Peak memory during the copy is LTERM + NTERM terms.

      ALLOCATE ( STF3_KEEP(NTERM) )
      STF3_KEEP(1:NTERM) = STF3(1:NTERM)
      CALL DEALLOCATE_STF_ARRAYS ( 'STF3' )

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         LTERM_KGGD = MAX(NTERM,1)                         ! Set LTERM before reallocating, so STF3 is reallocated with NTERM terms
      ELSE
         LTERM_KGG  = MAX(NTERM,1)
      ENDIF
      CALL ALLOCATE_STF_ARRAYS ( 'STF3', SUBR_NAME )

      STF3(1:NTERM) = STF3_KEEP(1:NTERM)
      DEALLOCATE ( STF3_KEEP )

! Reset LTERM and NTERM to appropriate values

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         NTERM_KGGD = NTERM
         LTERM_KGGD = NTERM_KGGD                           ! reset LTERM now that we have det. actual number of terms in KGGD
      ELSE
         NTERM_KGG  = NTERM
         LTERM_KGG  = NTERM_KGG                            ! reset LTERM now that we have det. actual number of terms in KGG
      ENDIF


! **********************************************************************************************************************************
! Debug output:

      IF((DEBUG(10) == 11) .OR. (DEBUG(10) == 12) .OR. (DEBUG(10) == 13) .OR.                                                     &
         (DEBUG(10) == 31) .OR. (DEBUG(10) == 32) .OR. (DEBUG(10) == 33)) THEN
         WRITE(F06,1260)
         MAX_NUM = MAX(NTERM,NDOFG)
         DO I=1,MAX_NUM
            IF      (MAX_NUM == NTERM) THEN
               IF (NDOFG >= I) THEN
                  WRITE(F06,1261) I,STFKEY(I),STF3(I)
               ELSE
                  WRITE(F06,1262) I, STF3(I)
               ENDIF
            ELSE IF (MAX_NUM == NDOFG) THEN
               IF (NTERM >= I) THEN
                  WRITE(F06,1261) I,STFKEY(I),STF3(I)
               ELSE
                  WRITE(F06,1263) I,STFKEY(I)
               ENDIF
            ENDIF
         ENDDO
         WRITE(F06,*)
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1260 FORMAT(/,'            I   STFKEY(I)   STFCOL(I)   STFPNT(I)           STF(I)')

 1261 FORMAT(1X,I12,I12,I12,I12,3X,1ES21.14)

 1262 FORMAT(1X,I12,12X,I12,I12,3X,1ES21.14)

 1263 FORMAT(1X,I12,I12)

 1624 FORMAT(' *ERROR  1624: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' TOO MANY NON-ZERO TERMS IN THE ',A,' MATRIX. LIMIT IS ',A,'    = ',I12)

 1628 FORMAT(' *ERROR  1628: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' BASED ON DEBUG ',I3,' VALUE, ARRAY ',A,' SHOULD BE ALLOCATED BUT IT IS NOT')

 9876 FORMAT(/,' PROCESSING ABORTED DUE TO ABOVE ',I8,' ELEMENT GENERATION ERRORS')

14001 FORMAT(' ********************************************************************************************************************&
&******')

14002 FORMAT('From subr ESP: TEMPLATE array showing ',I12,' actual nonzero terms in upper triangle of K')

35791 format(' In ESP: I, EID, J, K, BGRID(J), TDOF_ROW_NUM, 6*(BGRID(J)-1)+K, Diff = ',8i8)

88770 format(' In ESP:                         KGG_MAX_DIAG_TERM  = ',47X,1ES15.6)

88771 format(' In ESP #1: J, K, IS, NTERM_KGG, KGG(diagonal term) = ',4i8,1es15.6)

88772 format(' In ESP #2: J, K, IS, NTERM_KGG, KGG(diagonal term) = ',4i8,1es15.6)

88773 format(' In ESP #3: J, K, IS, NTERM_KGG, KGG(diagonal term) = ',4i8,1es15.6)

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DUMPSTF ( WHAT, J, K, KGG_ROW, KGG_COL, IS, ISS )

! Prints out info on the formulation of stiffness arrays for subr ESP, which generates the arrays

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, ERR, F06
      USE MODEL_STUF, ONLY            :  EID
      USE STF_ARRAYS, ONLY            :  STF3

      IMPLICIT NONE

      CHARACTER(1*BYTE), INTENT(IN)   :: WHAT              ! Indicator of where this subr was called from in subr ESP

      INTEGER(LONG)                   :: IS                ! A pointer into arrays STFKEY and STFPNT
      INTEGER(LONG)                   :: ISS               ! A particular value of IS
      INTEGER(LONG)    , INTENT(IN)   :: J                 ! Row number of elem stiff matrix term, KE(J,K), or KED(J,K)
      INTEGER(LONG)    , INTENT(IN)   :: K                 ! Col number of elem stiff matrix term, KE(J,K), or KED(J,K)
      INTEGER(LONG)    , INTENT(IN)   :: KGG_COL           ! Row number of KGG matrix where KE(J,K), or KED(J,K), goes
      INTEGER(LONG)    , INTENT(IN)   :: KGG_ROW           ! Col number of KGG matrix where KE(J,K), or KED(J,K), goes

! **********************************************************************************************************************************
      IF      (WHAT == '0') THEN

         WRITE(F06,8910)

      ELSE IF (WHAT == 'A') THEN

         WRITE(F06,8930) EID,J,K,KGG_ROW,KGG_COL,STFKEY(KGG_ROW),NTERM,STF3(NTERM)%Col_1,STF3(NTERM)%Col_2,                        &
                         STF3(NTERM)%Col_3,IS,ISS

      ELSE IF (WHAT == 'B') THEN

         WRITE(F06,8940) EID,J,K,KGG_ROW,KGG_COL,STFKEY(KGG_ROW),NTERM,STF3(NTERM)%Col_1,STF3(NTERM)%Col_2,                        &
                         STF3(NTERM)%Col_3,IS ,ISS,STF3(ISS)%Col_2,STF3(IS)%Col_3

      ELSE IF (WHAT == 'C') THEN

         WRITE(F06,8950) EID,J,K,KGG_ROW,KGG_COL,STFKEY(KGG_ROW),NTERM,STF3(NTERM)%Col_1,STF3(NTERM)%Col_2,                        &
                         STF3(NTERM)%Col_3,IS,ISS

      ENDIF

      RETURN

! **********************************************************************************************************************************
 8910 FORMAT(4X,'ELEM       J       K KGG_ROW KGG_COL  STFKEY  NKTERM  STFCOL  STFPNT      STF         IS     ISS  STFPNT     STF'&
          ,/,4X,' ID                                 (KGG_ROW)       (NKTERM) (NKTERM)   (NKTERM)                   (ISS)     (IS)')

 8920 FORMAT(1X,'---------------------------------------------------------------------------------------------------------------', &
'--------')

 8930 FORMAT(1X,'A',I6,I8,I8,I8,I8,I8,I8,I8,I8,1ES12.3,I8,I8)

 8940 FORMAT(1X,'B',I6,I8,I8,I8,I8,I8,I8,I8,I8,1ES12.3,I8,I8,I8,1ES12.3)

 8950 FORMAT(1X,'C',I6,I8,I8,I8,I8,I8,I8,I8,I8,1ES12.3,I8,I8)

! **********************************************************************************************************************************
      END SUBROUTINE DUMPSTF

! ##################################################################################################################################

      SUBROUTINE WRITE_NEG_DIAG_STIFFNESS ( WHAT )

      USE CONSTANTS_1,ONLY            :  ZERO
      USE MODEL_STUF, ONLY            :  AGRID, BGRID, EID, ELGP, ELDOF, TYPE

      IMPLICIT NONE

      INTEGER(LONG)                   :: II,JJ,KK,LL,MM    ! DO loop indices or counters
      INTEGER(LONG)                   :: NUM_DIAG_NEGS     ! Number of neg values on the diag of the quad element stiffness matrix
      INTEGER(LONG)                   :: WHAT              ! What header to write out

      REAL(DOUBLE)                    :: MAX_ABS_DIAG      ! Max absolute diagonal term
      REAL(DOUBLE)                    :: RATIO             ! Ratio of diagonal term to MAX_ABS_DIAG
      REAL(DOUBLE)                    :: ZE(ELDOF,ELDOF)   ! Either KE or KED (KED for BUCKLING)

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         DO II=1,ELDOF
            DO JJ=1,ELDOF
               ZE(II,JJ) = KED(II,JJ)
            ENDDO
         ENDDO
      ELSE
         DO II=1,ELDOF
            DO JJ=1,ELDOF
               ZE(II,JJ) = KE(II,JJ)
            ENDDO
         ENDDO
      ENDIF


      MAX_ABS_DIAG = ZERO
      NUM_DIAG_NEGS = 0
      DO II=1,ELDOF
         IF (DABS(ZE(II,II)) > MAX_ABS_DIAG) THEN
            MAX_ABS_DIAG = ZE(II,II)
         ENDIF
         IF (ZE(II,II) < 0.D0) THEN
            NUM_DIAG_NEGS = NUM_DIAG_NEGS + 1
         ENDIF
      ENDDO

      IF (NUM_DIAG_NEGS > 0) THEN

         IF (WHAT == 1) THEN
            WRITE(F06,97531) TYPE, EID, NUM_DIAG_NEGS
         ELSE IF (WHAT == 2) THEN
            WRITE(F06,97532) TYPE, EID, NUM_DIAG_NEGS
         ENDIF

         WRITE(F06,97533)
         KK=0
         DO LL=1,ELGP
            CALL GET_GRID_NUM_COMPS ( BGRID(LL), NUM_COMPS, SUBR_NAME )
            DO MM=1,NUM_COMPS
               KK = KK + 1
               RATIO = ZERO
               IF (MAX_ABS_DIAG > ZERO) THEN
                  RATIO = ZE(KK,KK)/MAX_ABS_DIAG
               ENDIF
               IF (ZE(KK,KK) >= ZERO) THEN
                  IF (MM == 1) THEN
                     WRITE(F06,90869) AGRID(LL), MM, ZE(KK,KK), RATIO
                  ELSE
                     WRITE(F06,90870) MM, ZE(KK,KK), RATIO
                  ENDIF
               ELSE
                  IF (MM == 1) THEN
                     WRITE(F06,90871) AGRID(LL), MM, ZE(KK,KK), RATIO
                  ELSE
                     WRITE(F06,90872) MM, ZE(KK,KK), RATIO
                  ENDIF
               ENDIF
            ENDDO
            WRITE(F06,*)
         ENDDO
         WRITE(F06,*)
      ENDIF

97531 FORMAT(' Diagonal stiffnesses for ',A,' element ',I8, ' in local elem coords. Element has ', I3,' negative diagonal terms')

97532 FORMAT(' Diagonal stiffnesses for ',A,' element ',I8, ' in global     coords. Element has ', I3,' negative diagonal terms')

97533 FORMAT( '    Grid  Comp  Diagonal stiffness      Ratio')

90869 FORMAT( I8,I6,1ES17.6,1ES17.6)

90870 FORMAT( 8X,I6,1ES17.6,1ES17.6)

90871 FORMAT( I8,I6,1ES17.6,1ES17.6,' *****NEGATIVE STIFFNESS*****')

90872 FORMAT( 8X,I6,1ES17.6,1ES17.6,' *****NEGATIVE STIFFNESS*****')

      END SUBROUTINE WRITE_NEG_DIAG_STIFFNESS

      END SUBROUTINE ESP


      SUBROUTINE KGG_SINGULARITY_PROC ( AGRID, KGRD, NUM_ASPC_BY_COMP )

! Grid point singularity processor. The algorithm is based on input matrix KGRD that is the 6x6 matrix from the diagonal of KGG for
! one grid point. The 2 3x3 diagonal partitions from KGRD are checked to see if there are any singularities based on:

!     1) Get eigenvalues and eigenvectors of each of the 3x3 matrices (one for translation and 1 for rotation). Calc the ratios of
!        the 3 eigenvales to the max value (among the 3) and, if the ratio is less than AUTOSPC_RAT, mark the DOF for AUTOSPC.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, SPC
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFSA, NGRID, NUM_PCHD_SPC1
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE PARAMS, ONLY                :  AUTOSPC, AUTOSPC_INFO, AUTOSPC_RAT, EPSIL, PCHSPC1, SPC1SID, SUPINFO
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START, TDOFI, TSET
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  GRID_ID

      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE TEXT_FIELD_UTILS, ONLY      :  CONVERT_INT_TO_CHAR
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'KGG_SINGULARITY_PROC'
      CHARACTER( 1*BYTE)              :: AUTOSPC_SOME_COMP  ! 'Y'/'N' indicator if some component of a grid is to be AUTOSPC'd
      CHARACTER( 1*BYTE)              :: CONSTR_COMP(6)     ! If DOF i is in the S or M sets, CONSTR_COMP(i) = 'Y'. Otherwise 'N'
      CHARACTER( 1*BYTE)              :: SINGLR_COMP(6)     ! If DOF i is singular, SINGLR_COMP(i) = 'Y'. Otherwise 'N'
      CHARACTER( 6*BYTE)              :: SINGLR_COMP_CHAR   ! Array of indicators of whether a displ comp in KGRD is singular. If
!                                                             DOF's 2, 4 and 5 are each singular then SINGLR_COMP_CHAR = ' 2 45 '
      CHARACTER( 1*BYTE)              :: WRITE_HEADER       ! 'Y'/'N' indicator for writing header to F06 for debug output

      INTEGER(LONG), INTENT(IN)       :: AGRID              ! Actual grid ID for IGRID
      INTEGER(LONG), INTENT(INOUT)    :: NUM_ASPC_BY_COMP(6)! The number of DOF's AUTOSPC'd for each displ component
      INTEGER(LONG)                   :: EIGENVAL_NUM(6)    ! Array to hold the eigenvalue number used in finding a SINGLR_COMP
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM    ! Row number in array GRID_ID where AGRID is found
      INTEGER(LONG)                   :: IDOF               ! Internal DOF number
      INTEGER(LONG)                   :: I,J,K              ! DO loop indices
      INTEGER(LONG)                   :: IGRID              ! Internal grid ID
      INTEGER(LONG)                   :: INFO               ! See subr K33_EIGENS (CONTAIN'ed herein)
      INTEGER(LONG)                   :: I2,J2              ! Index into an array
      INTEGER(LONG)                   :: NUM_COMPS          ! 6 if a physical grid,  if a scalar point (SPOINT)
      INTEGER(LONG)                   :: M_SET_COL          ! Col no. in array TDOF where the  M-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: S_SET_COL          ! Col no. in array TDOF where the  S-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: O_SET_COL          ! Col no. in array TDOF where the  O-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: R_SET_COL          ! Col no. in array TDOF where the  R-set is (from subr TDOF_COL_NUM)
      INTEGER(LONG)                   :: ROW_NUM_START      ! DOF number where TDOF data begins for a grid


      REAL(DOUBLE) , INTENT(IN)       :: KGRD(6,6)          ! 6x6 diagonal stiffness matrix for grid point AGRID
      REAL(DOUBLE)                    :: FAC                ! Multipling factor used in an intermediate calc
      REAL(DOUBLE)                    :: EPS1               ! Small value used in comparison to determine a real zero
      REAL(DOUBLE)                    :: K33(3,3)           ! Partition of the 6x6 stiff matrix for rotation    DOF's for a grid
      REAL(DOUBLE)                    :: K33_LAMBDAS(3)     ! Eigenvalues  of K33_ROT
      REAL(DOUBLE)                    :: K33_LAMBDA_MAX     ! Max value from  K33_ROT K33_LAMBDAS(i)
      REAL(DOUBLE)                    :: K33_VECS(3,3)      ! Eigenvectors of K33_ROT
      REAL(DOUBLE)                    :: K33_VEC_MAX        ! Max value from one eigenvector (column) of K33_VECS

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      SINGLR_COMP_CHAR(1:6) = ' '
      DO I=1,6
         SINGLR_COMP(I) = 'N'
         CONSTR_COMP(I) = 'N'
         EIGENVAL_NUM(I) = 0
      ENDDO

      CALL TDOF_COL_NUM ( 'M ',  M_SET_COL )
      CALL TDOF_COL_NUM ( 'S ',  S_SET_COL )
      CALL TDOF_COL_NUM ( 'O ',  O_SET_COL )
      CALL TDOF_COL_NUM ( 'R ',  R_SET_COL )

      WRITE_HEADER = 'Y'

! Check individual DOF's for singularities, and set SINGLR_COMP and CONSTR_COMP for this grid

!xx   CALL CALC_TDOF_ROW_NUM  ( AGRID, ROW_NUM_START, 'N' )! Det where in TDOF (not TDOFI) the DOF data begins for AGRID
      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, IGRID )
      ROW_NUM_START = TDOF_ROW_START(IGRID)
      CALL GET_GRID_NUM_COMPS ( IGRID, NUM_COMPS, 'N' )
comps:IF (NUM_COMPS == 6) THEN                             ! Physical grid with 6 components represented in KGRD

comps6:  DO K=1,2                                          ! K=1 is for translational DOF's and K=2 is for rotational DOF's

            DO I=1,3                                       ! K33 is the 3x3 partition of KGRD for translation or rotation
               I2 = I + 3*(K - 1)
               DO J=1,3
                  J2 = J + 3*(K - 1)
                  K33(I,J) = KGRD(I2,J2)
               ENDDO
            ENDDO
                                                           ! Get the 3 eigenvalues and eigenvectors of K33
            CALL K33_EIGENS ( K33, K33_LAMBDAS, K33_VECS, INFO )

eigs_ok:    IF (INFO == 0) THEN                            ! K33_EIGENS returned with no error from LAPACK routine

               K33_LAMBDA_MAX = ZERO                       ! Find max eigenvalue (max K33_LAMBDA_MAX)
               DO I=1,3
                  IF (DABS(K33_LAMBDAS(I)) > DABS(K33_LAMBDA_MAX)) THEN
                     K33_LAMBDA_MAX = K33_LAMBDAS(I)
                  ENDIF
               ENDDO

               IF (DABS(K33_LAMBDA_MAX) > EPS1) THEN       ! If max eigenvalue > 0, use it to  normalize all eigenvalues
                  FAC = ONE/K33_LAMBDA_MAX
               ELSE                                        ! If max eigenvalue = 0, use 1.0 to normalize all eigenvalues
                  FAC = ONE
               ENDIF

               DO J=1,3                                    ! For the J-th eigenvalue (1 to 3) get any comp that needs to be SPC'd
                  IF(DABS(FAC*K33_LAMBDAS(J))<=AUTOSPC_RAT) THEN
                     K33_VEC_MAX = ZERO
                     I2 = 0
                     DO I=1,3                              ! Scan the J-th eigenvector for comp with largest absolute value
                        IF (DABS(K33_VECS(I,J)) > K33_VEC_MAX) THEN
                           K33_VEC_MAX = DABS(K33_VECS(I,J))
                           I2 = I + 3*(K - 1)              ! After this loop, I2 will be the comp number (1-6) where the eigenvector
                        ENDIF                              ! is the max absolute value
                     ENDDO
                     SINGLR_COMP(I2) = 'Y'
                     EIGENVAL_NUM(I2) = J
                     CALL CONVERT_INT_TO_CHAR ( I2, SINGLR_COMP_CHAR(I2:I2) )
                  ENDIF
               ENDDO

               DO I=1,3                                    ! Check to see if any of the 3 comps (for K=1,2) are in M, S, O or R set
                  I2 = I + 3*(K - 1)                       ! (i.e., do not AUTOSPC DOF's that are members of the M, S, O or R-sets)
                  IDOF = ROW_NUM_START + I2 - 1
                  IF ((TDOF(IDOF,M_SET_COL) /= 0) .OR. (TDOF(IDOF,S_SET_COL) /= 0) .OR.                                            &
                      (TDOF(IDOF,O_SET_COL) /= 0) .OR. (TDOF(IDOF,R_SET_COL) /= 0)) THEN
                     CONSTR_COMP(I2) = 'Y'                 ! If there are, set CONSTR_COMP to 'Y'
                  ENDIF
               ENDDO

               AUTOSPC_SOME_COMP = 'N'
               DO I=1,6
                  IF ((SINGLR_COMP(I) == 'Y') .AND. (CONSTR_COMP(I) == 'N')) THEN
                     AUTOSPC_SOME_COMP = 'Y'
                  ENDIF
               ENDDO

deb_17:        IF (DEBUG(17) > 0) THEN

                  IF ((AUTOSPC_SOME_COMP == 'Y') .OR. (DEBUG(17) > 1)) THEN

                     IF (WRITE_HEADER == 'Y') THEN
                        CALL KGG_SING_PROC_DEBUG ( 1 )
                        WRITE_HEADER = 'N'
                     ENDIF

                     CALL KGG_SING_PROC_DEBUG ( 2 )

                     DO J=1,3
                        I2 = J + 3*(K - 1)
                        IF (SINGLR_COMP(I2) == 'Y') THEN
                           CALL KGG_SING_PROC_DEBUG ( 3 )
                        ENDIF
                     ENDDO

                     WRITE(F06,*)

                  ENDIF

               ENDIF deb_17

            ENDIF eigs_ok

         ENDDO comps6

      ELSE                                                 ! Scalar point - only need to check KGRD(1,1)

         IDOF = ROW_NUM_START
         IF ((TDOF(IDOF,M_SET_COL) /= 0) .OR. (TDOF(IDOF,S_SET_COL) /= 0) .OR.                                                     &
             (TDOF(IDOF,O_SET_COL) /= 0) .OR. (TDOF(IDOF,R_SET_COL) /= 0)) THEN
            CONSTR_COMP(1) = 'Y'                           ! If this is M, S, O, R set, set CONSTR_COMP to 'Y'
         ENDIF

         IF (KGRD(1,1) < EPS1) THEN
            SINGLR_COMP(1) = 'Y'
            CALL CONVERT_INT_TO_CHAR ( 1, SINGLR_COMP_CHAR(1:1) )
         ENDIF

         AUTOSPC_SOME_COMP = 'N'
         DO I=1,6
            IF ((SINGLR_COMP(I) == 'Y') .AND. (CONSTR_COMP(I) == 'N')) THEN
               AUTOSPC_SOME_COMP = 'Y'
            ENDIF
         ENDDO

      ENDIF comps


      DO I=1,6                                             ! Change SINGLR_COMP_CHAR if a component belongs to S or M set
         IF (CONSTR_COMP(I) == 'Y') THEN
            SINGLR_COMP_CHAR(I:I) = ' '
         ENDIF
      ENDDO

      DO I=1,6                                             ! Reset table TSET for the DOF's AUTOSPC'd and increment NUM_ASPC_BY_COMP
         IF (SINGLR_COMP_CHAR(I:I) /= ' ') THEN
            IF (AUTOSPC == 'Y') THEN
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID, GRID_ID_ROW_NUM )
               TSET ( GRID_ID_ROW_NUM, I ) = 'SA'
               NDOFSA = NDOFSA + 1
               NUM_ASPC_BY_COMP(I) = NUM_ASPC_BY_COMP(I) + 1
            ENDIF
         ENDIF
      ENDDO

      IF (AUTOSPC_INFO == 'Y') THEN                        ! Write singularity messages to output file if requested
         IF (SINGLR_COMP_CHAR /= '      ') THEN
            IF (AUTOSPC == 'Y') THEN
               WRITE(ERR,101) AGRID,(SINGLR_COMP_CHAR(I:I),I=1,6),AUTOSPC
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,101) AGRID,(SINGLR_COMP_CHAR(I:I),I=1,6),AUTOSPC
               ENDIF
            ELSE
               WRITE(ERR,102) AGRID,(SINGLR_COMP_CHAR(I:I),I=1,6)
               IF (SUPINFO == 'N') THEN
                  WRITE(F06,102) AGRID,(SINGLR_COMP_CHAR(I:I),I=1,6)
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      DO I=1,6                                             ! If requested, write SPC1 card images to text file for singulaqr DOF's
         IF (SINGLR_COMP_CHAR(I:I) /= ' ') THEN
            IF (PCHSPC1 == 'Y') THEN
               WRITE(SPC,109) SPC1SID, I, AGRID
               NUM_PCHD_SPC1 = NUM_PCHD_SPC1 + 1
            ENDIF
         ENDIF
      ENDDO



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *INFORMATION: GRID POINT ',I8,' HAS SINGULARITY FOR DISPL COMPONENT(S) ',6A1,'. SINCE PARAM AUTOSPC = ',A,          &
                          ', THESE WILL BE AUTOSPC''d')

  102 FORMAT(' *INFORMATION: GRID POINT ',I8,' HAS SINGULARITY FOR DISPL COMPONENT(S) ',6A1)

  109 FORMAT('SPC1    ',3I8)


! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE K33_EIGENS (K33, K33_LAMBDAS, K33_VECS, INFO )

! Jacobi solution for 3x3 eigenvalue problem used in finding the eigenvalues of a 3x3 diag partition of a 6x6 grid stiffness matrix

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE CONSTANTS_1, ONLY           :  ZERO

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: CALLED_SUBR = ' ' ! Name of a called subr (for output error purposes)
      CHARACTER( 1*BYTE), PARAMETER   :: JOBZ      = 'V'   ! Indicates to solve for eigenvalues and vectors in LAPACK subr DSYEV
      CHARACTER( 1*BYTE), PARAMETER   :: UPLO      = 'U'   ! Indicates array A is the upper triangular part of K33

      INTEGER(LONG), INTENT(OUT)      :: INFO              ! = 0:  successful exit from subr DSYEV
!                                                            < 0:  if INFO = -i, the i-th argument had an illegal value
!                                                            > 0:  if INFO = i, the algorithm failed to converge; i off_diag
!                                                            elems of an intermediate tridiag form did not converge to zero
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG), PARAMETER        :: N         = 3     ! Order of matrix K33
      INTEGER(LONG), PARAMETER        :: LWORK     = 3*N-1 ! Size of array WORK

      REAL(DOUBLE) , INTENT(IN)       :: K33(N,N)          ! A 3x3 diag partition of the 6x6 stiff matrix for a grid
      REAL(DOUBLE) , INTENT(OUT)      :: K33_LAMBDAS(N)    ! The eigenvalues of input K33 (if INFO = 0)
      REAL(DOUBLE) , INTENT(OUT)      :: K33_VECS(N,N)     ! Prior to entry to DSYEV, K33_VECS is set = K33.
!                                                            On exit, K33_VECS contains the eigenvectors of K33
      REAL(DOUBLE)                    :: WORK(LWORK)       ! Workspace for subr DSYEV

! **********************************************************************************************************************************
! Initialize outputs

      INFO = 0

      DO I=1,N
         K33_LAMBDAS(I) = ZERO
      ENDDO

      DO I=1,N
         DO J=1,N
            K33_VECS(I,J) = ZERO
         ENDDO
      ENDDO

! Use LAPACK driver DSYEV to get all eigenvalues of K33

      DO I=1,N
         DO J=1,N
            K33_VECS(I,J) = K33(I,J)
         ENDDO
      ENDDO

      CALL DSYEV ( JOBZ, UPLO, N, K33_VECS, N, K33_LAMBDAS, WORK, LWORK, INFO )
      CALLED_SUBR = 'DSYEV'

      IF      (INFO < 0) THEN                              ! LAPACK subr XERBLA should have reported error on an illegal argument
                                                           ! in a call to a LAPACK subr, which would only occur if the call to DSYEV
         WRITE(ERR,993) SUBR_NAME, CALLED_SUBR             ! was incorrect.
         WRITE(F06,993) SUBR_NAME, CALLED_SUBR
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )

      ELSE IF (INFO > 0) THEN                              ! No convergence in subr DSYEV. Return INFO > 0 for calling routine to
                                                           ! deal with
         WRITE(ERR,1612) CALLED_SUBR, SUBR_NAME
         WRITE(F06,1612) CALLED_SUBR, SUBR_NAME
         FATAL_ERR = FATAL_ERR + 1

      ENDIF

      RETURN

! **********************************************************************************************************************************
  993 FORMAT(' *ERROR   993: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' LAPACK SUBR XERBLA SHOULD HAVE REPORTED AN ERROR ON AN ILLEGAL ARGUMENT IN A CALL TO LAPACK SUBR '    &
                    ,/,15X,A,' (OR A SUBR CALLED BY IT) AND THEN ABORTED')

 1612 FORMAT(' *ERROR  1612: LAPACK DRIVER ',A8,' CALLED BY SUBROUTINE ',A                                                         &
                    ,/,14X,' CANNOT CONVERGE IN ATTEMPTING TO FIND K33 EIGENVALUES AND EIGENVECTORS'                               &
                    ,/,14X,' THE ALGORITHM HAS FAILED TO FIND ALL THE EIGENVALUES (PRINCIPAL MOIs) IN 90 ITERATIONS')

! **********************************************************************************************************************************

      END SUBROUTINE K33_EIGENS

! ##################################################################################################################################

      SUBROUTINE KGG_SING_PROC_DEBUG ( WHAT )

! Debug output for KGG singularity calcs

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR

      IMPLICIT NONE

      INTEGER(LONG)                   :: COMP_NUM          ! Displ component number (1 - 6)
      INTEGER(LONG)                   :: II,JJ             ! DO loop indices
      INTEGER(LONG)                   :: WHAT              ! What to print out

! **********************************************************************************************************************************
      IF     ( WHAT == 1) THEN

         WRITE(F06,201) NDOFSA+1, AGRID                    ! Use NDOFSA+1 since NDOFSA is not incremented until later

      ELSE IF (WHAT == 2) THEN

         IF      (K == 1) THEN
            WRITE(F06,202)
            WRITE(F06,204)
         ELSE IF (K == 2) THEN
            WRITE(F06,203)
            WRITE(F06,204)
         ENDIF

         DO II=1,3
            COMP_NUM = II + 3*(K - 1)
            WRITE(F06,205) COMP_NUM,(K33(II,JJ),JJ=1,3),II,K33_LAMBDAS(II),(K33_VECS(II,JJ),JJ=1,3)
         ENDDO
         WRITE(F06,*)

      ELSE IF (WHAT == 3) THEN

         WRITE(F06,101) I2, EIGENVAL_NUM(I2), I2, EIGENVAL_NUM(I2), I2

      ENDIF

! **********************************************************************************************************************************
  101 FORMAT(' AUTOSPC comp',I2,' since (eigenvalue',I2,')/max eigenvalue < AUTOSPC_RAT (comp',I2,' is AUTOSPC''d since eigenvec', &
             I2,' is max abs for comp',I2,')')

  201 FORMAT(' ******************************************************************************************************************',&
             '*****************'                                                                                                   &
          ,/,I8, ')', ' Results from KGG_SINGULARITY_PROC for grid ',I8                                                            &
          ,/,' ------------------------------------------------------------')

  202 FORMAT(50X,'Results for translation 3x3 partitions')

  203 FORMAT(50X,' Results for rotation 3x3 partitions')

  204 FORMAT(4X,'Comp',16X,'G-set stiffness',24X,'Eigenvalues',28X,'Eigenvectors',/,93X,'1             2             3')

  205 FORMAT(5X,I2,3X,3(1ES14.6),8X,I2,1ES14.6,10X,3(1ES14.6))

! **********************************************************************************************************************************

      END SUBROUTINE KGG_SING_PROC_DEBUG

      END SUBROUTINE KGG_SINGULARITY_PROC


      SUBROUTINE SPARSE_KGG

! (1) Converts the system KGG matrix from a sparse linked list format to a row, col, val format. It sorts each row to be
!     in G-set DOF numerical order. The transformed matrix is written out to file LINK1L row by row.

! (2) Call KGG_SINGULARITY_PROC to check 6x6 diagonal stiffness partitions for singularities (which resets TSET table if there are
!     any grid point singularities)

! (3) Call TDOF_PROC to regenerate TDOF, TDOFI tables if KGG_SINGULARITY_PROC found singularities

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1L, L1L_MSG, LINK1L, SC1, SPCFIL, SPC, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NGRID, NIND_GRDS_MPCS,                                    &
                                         NTERM_KGG, NUM_PCHD_SPC1, SOL_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  AUTOSPC, AUTOSPC_RAT, EPSIL, PRTTSET, PRTSTIFF, SPC1QUIT, SUPINFO, SUPWARN
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID, GRID_SEQ, MPC_IND_GRIDS, INV_GRID_SEQ
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START, TDOFI, TSET
      USE STF_ARRAYS, ONLY            :  STFKEY, STF3
      USE SPARSE_MATRICES, ONLY       :  I_KGG, J_KGG, KGG
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE FILE_LIFECYCLE, ONLY   :  OPNERR, OUTA_HERE
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE FILE_LIFECYCLE, ONLY        :  FILERR, FILE_CLOSE, FILE_OPEN
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM, TDOF_PROC
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE SORTING, ONLY               :  SORT_INT1_REAL1
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  AUTOSPC_SUMMARY_MSGS
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_KGG'
      CHARACTER(  7*BYTE)             :: ASPC_SUM_MSG1      ! Message to be printed out in the AUTOSPC summary table
      CHARACTER(100*BYTE)             :: ASPC_SUM_MSG2      ! Message to be printed out in the AUTOSPC summary table
      CHARACTER( 13*BYTE)             :: ASPC_SUM_MSG3      ! Message to be printed out in the AUTOSPC summary table
      CHARACTER(132*BYTE)             :: TDOF_MSG           ! Message to be printed out regarding at what pt in the run the TDOF,I
!                                                             tables are printed out
      CHARACTER( 1*BYTE)              :: SKIPIT      = 'N'  ! Indicator of whether to skip KGG sing proc (for MPC indep grids)

      INTEGER(LONG)                   :: AGRIDI             ! Actual grid ID
      INTEGER(LONG)                   :: G_SET_COL          ! Col in TDOF where G-set DOF's are
      INTEGER(LONG)                   :: I,J,K,L,N          ! DO loop indices
      INTEGER(LONG)                   :: IGRID              ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK              ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: IS                 ! Index into array STF3
      INTEGER(LONG)                   :: KGG_COL_NUM        ! The col num in G-set stiff matrix where stiff for DOF I begins
      INTEGER(LONG)                   :: KGG_ROW_NUM        ! The row num in G-set stiff matrix where stiff for DOF I begins
      INTEGER(LONG)                   :: KGG_II_COL_NUM     ! Col number in the 6x6 stiff matrix for 1 grid
      INTEGER(LONG)                   :: KGG_NUM_ASPC       ! Sum of NUM_ASPC_BY_COMP(6) (this is also NDOFSA but we need to test
      INTEGER(LONG)                   :: KTERM_KGG          ! Count of terms written to KGG file LINK1L to compare with NTERM_KGG
      INTEGER(LONG)                   :: NUM_NONZERO_IN_ROW ! Count of the actual number of nonzero terms in a row of KGG
      INTEGER(LONG)                   :: NUM_ASPC_BY_COMP(6)! The number of SPC1's for each displ component
      INTEGER(LONG)                   :: NUM_MAX = 0        ! largest number of terms in any row of the KGG stiffness matrix
      INTEGER(LONG)                   :: NUM_COMPS          ! Number of displ components (1 for SPOINT, 6 for physical grid)
      INTEGER(LONG)                   :: NZERO   = 0        ! Count on zero terms in array STF
      INTEGER(LONG)                   :: OUNT(2)            ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: ROW_NUM_START      ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: RJ(NDOFG)          ! Column numbers corresponding to the terms in RSTF(I).


      REAL(DOUBLE)                    :: EPS1               ! A small number to compare real zero
      REAL(DOUBLE)                    :: KGG_II(6,6)        ! 6 x 6 diagonal stiffness matrices for 1 grid
      REAL(DOUBLE)                    :: RSTF(NDOFG)        ! 1D array of terms from STF(I) pertaining to one row of the G-set
!                                                             stiffness matrix. Initially, the cols are not in increasing global
!                                                             DOF order. RSTF is sorted, prior to writing the G-set stiff matrix
!                                                             to file LINK1L, so that the cols are in increasing DOF order.


      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Pass # 1: Determine final NTERM_KGG (may be less due to terms stripped)

      NZERO = 0
      DO I = 1,NDOFG                                        ! Start conversion.
         IS = STFKEY(I)

         IF (IS == 0) CYCLE                                 ! Check for null row in stiffness matrix and CYCLE if it is

         NUM_NONZERO_IN_ROW = 0                             ! Count zero terms so we can debit NTERM_KGG before writing it to file
         DO J = 1,NDOFG
            IF (DABS(STF3(IS)%Col_3) < EPS1) THEN
               NZERO = NZERO + 1
            ELSE
               NUM_NONZERO_IN_ROW = NUM_NONZERO_IN_ROW + 1
            ENDIF
            IS = STF3(IS)%Col_2
            IF (IS == 0) THEN
               EXIT
            ENDIF
         ENDDO
         IF (NUM_NONZERO_IN_ROW > NUM_MAX) THEN
            NUM_MAX = NUM_NONZERO_IN_ROW
         ENDIF
         IF (IS /= 0) THEN
            WRITE(ERR,1625) SUBR_NAME,I
            WRITE(F06,1625) SUBR_NAME,I
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                          ! Coding error, so quit
         ENDIF

      ENDDO


      NTERM_KGG = NTERM_KGG - NZERO

      WRITE(ERR,146) NTERM_KGG
      IF (SUPINFO == 'N') THEN
         WRITE(F06,146) NTERM_KGG
      ENDIF

      CALL ALLOCATE_SPARSE_MAT ( 'KGG', NDOFG, NTERM_KGG, SUBR_NAME )

! **********************************************************************************************************************************
! Pass # 2: Reformulate rows and write to file LINK1L. Call KGG_SINGULARITY_PROC


      IF (NTERM_KGG <= 0) THEN
         WRITE(ERR,1611) NTERM_KGG
         WRITE(F06,1611) NTERM_KGG
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                             ! Quit if the no. nonzero terms in KGG is <= 0
      ENDIF

! Open L1L to write stiffness.

      OUNT(1) = ERR
      OUNT(2) = F06
      CALL FILE_OPEN ( L1L, LINK1L, OUNT, 'REPLACE', L1L_MSG, 'WRITE_STIME', 'UNFORMATTED', 'WRITE', 'REWIND', 'Y', 'N' )
      WRITE(L1L) NTERM_KGG

! Open SPC to write SPC1 records if KGG_SINGULARITY_PROC finds singularities

      OPEN (SPC,FILE=SPCFIL,STATUS='REPLACE',IOSTAT=IOCHK)
      IF (IOCHK /= 0) THEN
         CALL OPNERR ( IOCHK, SPCFIL, OUNT )
         CALL FILERR ( OUNT )
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      DO I=1,6                                             ! Initialize NUM_ASPC_BY_COMP. It will get accumulated in KGG_SING_PROC
         NUM_ASPC_BY_COMP(I) = 0
      ENDDO

      IF (DEBUG(17) > 0) THEN                              ! Write leading seperator for DEBUG output
         WRITE(F06,9901) AUTOSPC_RAT
      ENDIF
!xx   WRITE(SC1,*)
      KTERM_KGG = 0                                        ! I runs over the number of rows (or grids)
      CALL TDOF_COL_NUM ( 'G ', G_SET_COL )
      KGG_ROW_NUM = 0
      I_KGG(1) = 1
      CALL COUNTER_INIT('     Working on grid ', NGRID)
i_do: DO I = 1,NGRID
         SKIPIT = 'N'

         DO K=1,6                                          ! Make KGG_II 6x6 even though for SPOINT's we only use 1-1 term
            DO L=1,6
               KGG_II = ZERO
            ENDDO
         ENDDO

!xx      CALL CALC_TDOF_ROW_NUM ( GRID_ID(INV_GRID_SEQ(I)), IROW_START, 'N' )
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_ID(INV_GRID_SEQ(I)), IGRID )
         ROW_NUM_START = TDOF_ROW_START(IGRID)
         KGG_COL_NUM = TDOF(ROW_NUM_START,G_SET_COL)
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
k_do:    DO K=1,NUM_COMPS

            KGG_ROW_NUM = KGG_ROW_NUM + 1
            IS = STFKEY(KGG_ROW_NUM)

            IF (IS == 0) THEN                              ! Check for null row in stiffness matrix
               I_KGG(KGG_ROW_NUM+1) = I_KGG(KGG_ROW_NUM)
               CYCLE k_do
            ENDIF

            NUM_NONZERO_IN_ROW = 0                         ! Form row of non-zero's in arrays RJ, RSTF
j_do1:      DO J=1,NDOFG
               IF (DABS(STF3(IS)%Col_3) >= EPS1) THEN
                  NUM_NONZERO_IN_ROW = NUM_NONZERO_IN_ROW + 1
                  RSTF(NUM_NONZERO_IN_ROW) = STF3(IS)%Col_3
                  RJ(NUM_NONZERO_IN_ROW)   = STF3(IS)%Col_1
               ENDIF
               IS = STF3(IS)%Col_2
               IF (IS == 0) THEN
                  EXIT j_do1
               ENDIF
            ENDDO j_do1

            IF (NUM_NONZERO_IN_ROW > NUM_MAX) THEN
               NUM_MAX = NUM_NONZERO_IN_ROW
            ENDIF
            IF (IS /= 0) THEN
               WRITE(ERR,1625) SUBR_NAME,I
               WRITE(F06,1625) SUBR_NAME,I
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                       ! Coding error, so quit
            ENDIF

            IF (NUM_NONZERO_IN_ROW /= 1) THEN               ! Sort row by the shell method so that RJ is in numerical order
               CALL SORT_INT1_REAL1 ( SUBR_NAME, 'RJ, RSTF', NUM_NONZERO_IN_ROW, RJ, RSTF )
            ENDIF


n_do:       DO N=1,NUM_NONZERO_IN_ROW                      ! Formulate the K-th row of KGG_II
               IF ((RJ(N) >= KGG_COL_NUM) .AND. (RJ(N) <= KGG_COL_NUM+NUM_COMPS-1)) THEN
                  KGG_II_COL_NUM = RJ(N) - (KGG_COL_NUM - 1)
                  KGG_II(K,KGG_II_COL_NUM) = RSTF(N)
               ENDIF
            ENDDO n_do

j_do3:      DO J=1,NUM_NONZERO_IN_ROW
               KTERM_KGG = KTERM_KGG + 1                    ! KTERM_KGG is a count on the no. records written
               WRITE(L1L) KGG_ROW_NUM, RJ(J), RSTF(J)
               J_KGG(KTERM_KGG) = RJ(J)
                 KGG(KTERM_KGG) = RSTF(J)
            ENDDO j_do3

            I_KGG(KGG_ROW_NUM+1) = I_KGG(KGG_ROW_NUM) + NUM_NONZERO_IN_ROW

         ENDDO k_do

         DO K=1,6                                           ! Set lower portion of KGG_II to be symmetric
            DO J=1,K-1
               KGG_II(K,J) = KGG_II(J,K)
            ENDDO
         ENDDO


         AGRIDI = GRID_ID(INV_GRID_SEQ(I))               ! We don't want to call KGG_SINGULARITY_PROC for grids that are indep
j_do4:   DO J=1,NIND_GRDS_MPCS                           ! on MPC's since they may have zero stiff at this point and get stiff
            IF (AGRIDI == MPC_IND_GRIDS(J)) THEN         ! later when MPC's are eliminated. These grids will get SKIPIT = 'Y'
               SKIPIT = 'Y'
               EXIT j_do4
            ENDIF
         ENDDO j_do4
         IF (SKIPIT == 'N') THEN
            CALL KGG_SINGULARITY_PROC ( AGRIDI, KGG_II, NUM_ASPC_BY_COMP )
         ENDIF
         CALL COUNTER_PROGRESS(I)
      ENDDO i_do
      IF (DEBUG(17) > 0) THEN                              ! Write trailing seperator for DEBUG output
         WRITE(F06,9902)
      ENDIF
      WRITE(F06,*)

      WRITE(SC1,*) CR13

      IF (PRTSTIFF(1) >= 1) THEN
         CALL WRITE_SPARSE_CRS ( 'STIFFNESS MATRIX KGG' , 'G ', 'G ', NTERM_KGG, NDOFG, I_KGG, J_KGG, KGG )
      ENDIF

! If AUTOSPC = Y reprint TSET table, and reset and reprint TDOF tables if we have AUTOSPC'd any DOF's

      IF (AUTOSPC == 'Y') THEN

         KGG_NUM_ASPC = 0
         DO I=1,6
            KGG_NUM_ASPC = KGG_NUM_ASPC + NUM_ASPC_BY_COMP(I)
         ENDDO

         IF (KGG_NUM_ASPC > 0) THEN

            IF (PRTTSET > 0) THEN
               WRITE(F06,56)
               WRITE(F06,57)
               DO J = 1,NGRID
                  WRITE(F06,58) GRID(J,1), GRID_SEQ(J), (TSET(J,K),K = 1,6)
               ENDDO
               WRITE(F06,'(//)')
            ENDIF

            ASPC_SUM_MSG1(1:) = 'Stage 1:'
            ASPC_SUM_MSG2(1:) = 'after identification of AUTOSPC''s at the grid level'
            ASPC_SUM_MSG3(1:) = 'in this stage'
            CALL AUTOSPC_SUMMARY_MSGS ( ASPC_SUM_MSG1, ASPC_SUM_MSG2, ASPC_SUM_MSG3, 'Y', NUM_ASPC_BY_COMP )

            TDOF_MSG(1:)  = ' '
            TDOF_MSG(39:) = ASPC_SUM_MSG2(1:)
            CALL TDOF_PROC ( TDOF_MSG )

         ENDIF

      ENDIF

      IF (NUM_PCHD_SPC1 > 0) THEN                       ! Close SPC file and, if any records were written to it, save it
         CALL FILE_CLOSE ( SPC, SPCFIL, 'KEEP' )
         IF (SPC1QUIT == 'Y') THEN
            WRITE(ERR,9991) SUBR_NAME, SPC1QUIT
            WRITE(F06,9991) SUBR_NAME, SPC1QUIT
            CALL OUTA_HERE ( 'Y' )
         ENDIF
      ELSE
         CALL FILE_CLOSE ( SPC, SPCFIL, 'DELETE' )
      ENDIF

      IF (KTERM_KGG /= NTERM_KGG) THEN                      ! Check KTERM_KGG = NTERM_KGG
         WRITE(ERR,1623) SUBR_NAME, LINK1L, KTERM_KGG, NTERM_KGG
         WRITE(F06,1623) SUBR_NAME, LINK1L, KTERM_KGG, NTERM_KGG
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                             ! Coding error, so quit
      ENDIF

      CALL FILE_CLOSE ( L1L, LINK1L, 'KEEP' )
      WRITE(ERR,101) NUM_MAX
      IF (SUPINFO == 'N') THEN
         WRITE(F06,101) NUM_MAX
      ENDIF



      RETURN

! **********************************************************************************************************************************
   56 FORMAT(64X,'DEGREE OF FREEDOM SET TABLE (TSET)')

   57 FORMAT(33x,'     GRID SEQUENCE       T1       T2       T3       R1       R2       R3',/)

   58 FORMAT(33x,2(1X,I8),6(7X,A2))

  146 FORMAT(' *INFORMATION: NUMBER OF NONZERO TERMS IN THE KGG STIFFNESS MATRIX IS                 = ',I12,/)

  101 FORMAT(' *INFORMATION: MAX NUMBER OF NONZERO TERMS IN A ROW OF THE G-SET STIFFNESS MATRIX     = ',I12,/)

 1625 FORMAT(' *ERROR  1625: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' 1ST COL OF ARRAY STF3 INDICATES THERE IS MORE DATA IN ARRAY STF3 FOR ROW ',I12,' OF THE KGG STIFF'    &
                    ,/,14X,' MATRIX ALTHOUGH THE DOF COUNT IS AT THE END OF THE ROW')

 1611 FORMAT(' *ERROR  1611: THE G-SET STIFFNESS MATRIX, KGG, MUST HAVE SOME NONZERO TERMS. HOWEVER IT HAS ',I12,' TERMS')

 1623 FORMAT(' *ERROR  1623: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF G-SET STIFFNESS MATRIX RECORDS WRITTEN TO FILE:'                                        &
                    ,/,15X,A                                                                                                       &
                    ,/,14X,' WAS KTERM_KGG = ',I12,'. IT SHOULD HAVE BEEN NTERM_KGG = ',I12)

 9901 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::START DEBUG(17) OUTPUT FROM SUBROUTINE KGG_SINGULARITY_PROC:::::::::::::::::::',&
             ':::::::::::::::::',//,54X,'AUTOSPC_RAT = ',1ES13.6,/)

 9902 FORMAT(' :::::::::::::::::::::::::::::::::::::END DEBUG(17) OUTPUT FROM SUBROUTINE KGG_SINGULARITY_PROC::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

 9991 FORMAT(' PROCESSING ABORTED IN SUBR ',A,' BASED ON PARAMETER SPC1QUIT = ',A)




99111 format('         I      GRID     COMPS   KGG_COL    II_ROW     II_COL        N           KGG_II',/,                          &
             '         -      ----     -----   -------    ------     ------        -        ------------')
99121 format(7i10,1es20.6)

! **********************************************************************************************************************************

      END SUBROUTINE SPARSE_KGG


      SUBROUTINE SPARSE_KGGD

! Converts the system KGGD differential stiff matrix from a sparse linked list format to a sparse row, col, val format. It sorts
! each row to be in G-set DOF numerical order.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1, SPCFIL, SPC, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NGRID, NIND_GRDS_MPCS,                                    &
                                         NTERM_KGGD, NUM_PCHD_SPC1, SOL_NAME, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE PARAMS, ONLY                :  AUTOSPC, AUTOSPC_RAT, EPSIL, PRTSTIFF, SPC1QUIT, SUPINFO, SUPWARN
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  GRID, GRID_ID, GRID_SEQ, MPC_IND_GRIDS, INV_GRID_SEQ
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START, TDOFI, TSET
      USE STF_ARRAYS, ONLY            :  STFKEY, STF3
      USE SPARSE_MATRICES, ONLY       :  I_KGGD, J_KGGD, KGGD

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SPARSE_MATRIX_ALLOCATION, ONLY:  ALLOCATE_SPARSE_MAT
      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS
      USE SORTING, ONLY               :  SORT_INT1_REAL1
      USE MATRIX_FILE_IO, ONLY        :  WRITE_SPARSE_CRS
      USE PROGRESS_COUNTERS, ONLY     :  COUNTER_INIT, COUNTER_PROGRESS

      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SPARSE_KGGD'

      INTEGER(LONG)                   :: G_SET_COL          ! Col in TDOF where G-set DOF's are
      INTEGER(LONG)                   :: I,J,K,L,N          ! DO loop indices
      INTEGER(LONG)                   :: IGRID              ! Internal grid ID
      INTEGER(LONG)                   :: IS                 ! Index into array STF3
      INTEGER(LONG)                   :: KGGD_COL_NUM       ! The col num in G-set stiff matrix where stiff for DOF I begins
      INTEGER(LONG)                   :: KGGD_ROW_NUM       ! The row num in G-set stiff matrix where stiff for DOF I begins
      INTEGER(LONG)                   :: KGGD_II_COL_NUM    ! Col number in the 6x6 stiff matrix for 1 grid
      INTEGER(LONG)                   :: KTERM_KGGD         ! Count of terms written to KGGD to compare with NTERM_KGGD
      INTEGER(LONG)                   :: NUM_NONZERO_IN_ROW ! Count of the actual number of nonzero terms in a row of KGGD
      INTEGER(LONG)                   :: NUM_MAX = 0        ! largest number of terms in any row of the KGGD stiffness matrix
      INTEGER(LONG)                   :: NUM_COMPS          ! Number of displ components (1 for SPOINT, 6 for physical grid)
      INTEGER(LONG)                   :: NZERO   = 0        ! Count on zero terms in array STF
      INTEGER(LONG)                   :: ROW_NUM_START      ! DOF number where TDOF data begins for a grid
      INTEGER(LONG)                   :: RJ(NDOFG)          ! Column numbers corresponding to the terms in RSTF(I).


      REAL(DOUBLE)                    :: EPS1               ! A small number to compare real zero
      REAL(DOUBLE)                    :: KGGD_II(6,6)       ! 6 x 6 diagonal stiffness matrices for 1 grid
      REAL(DOUBLE)                    :: RSTF(NDOFG)        ! 1D array of terms from STF(I) pertaining to one row of the G-set
!                                                             stiffness matrix. Initially, the cols are not in increasing global
!                                                             DOF order. RSTF is sorted so that the cols are in incr DOF order.

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Pass # 1: Determine final NTERM_KGGD (may be less due to terms stripped)

      NZERO = 0
      DO I = 1,NDOFG                                        ! Start conversion.
         IS = STFKEY(I)

         IF (IS == 0) CYCLE                                 ! Check for null row in stiffness matrix and CYCLE if it is

         NUM_NONZERO_IN_ROW = 0                             ! Count zero terms so we can debit NTERM_KGGD before writing it to file
         DO J = 1,NDOFG
            IF (DABS(STF3(IS)%Col_3) < EPS1) THEN
               NZERO = NZERO + 1
            ELSE
               NUM_NONZERO_IN_ROW = NUM_NONZERO_IN_ROW + 1
            ENDIF
            IS = STF3(IS)%Col_2
            IF (IS == 0) THEN
               EXIT
            ENDIF
         ENDDO
         IF (NUM_NONZERO_IN_ROW > NUM_MAX) THEN
            NUM_MAX = NUM_NONZERO_IN_ROW
         ENDIF
         IF (IS /= 0) THEN
            WRITE(ERR,1625) SUBR_NAME,I
            WRITE(F06,1625) SUBR_NAME,I
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                          ! Coding error, so quit
         ENDIF

      ENDDO


      NTERM_KGGD = NTERM_KGGD - NZERO

      WRITE(ERR,146) NTERM_KGGD
      IF (SUPINFO == 'N') THEN
         WRITE(F06,146) NTERM_KGGD
      ENDIF

      CALL ALLOCATE_SPARSE_MAT ( 'KGGD', NDOFG, NTERM_KGGD, SUBR_NAME )

! **********************************************************************************************************************************
! Pass # 2: Reformulate rows


      IF (NTERM_KGGD <= 0) THEN
         WRITE(ERR,1611) NTERM_KGGD
         WRITE(F06,1611) NTERM_KGGD
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                             ! Quit if the no. nonzero terms in KGGD is <= 0
      ENDIF

      KTERM_KGGD = 0                                        ! I runs over the number of rows (or grids)
      CALL TDOF_COL_NUM ( 'G ', G_SET_COL )
      KGGD_ROW_NUM = 0
      I_KGGD(1) = 1

!xx   WRITE(SC1, * )
      CALL COUNTER_INIT('     Working on grid ', NGRID)
i_do: DO I = 1,NGRID

         DO K=1,6                                          ! Make KGGD_II 6x6 even though for SPOINT's we only use 1-1 term
            DO L=1,6
               KGGD_II = ZERO
            ENDDO
         ENDDO

!xx      CALL CALC_TDOF_ROW_NUM ( GRID_ID(INV_GRID_SEQ(I)), IROW_START, 'N' )
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, GRID_ID(INV_GRID_SEQ(I)), IGRID )
         ROW_NUM_START = TDOF_ROW_START(IGRID)
         KGGD_COL_NUM = TDOF(ROW_NUM_START,G_SET_COL)
         CALL GET_GRID_NUM_COMPS ( INV_GRID_SEQ(I), NUM_COMPS, SUBR_NAME )
k_do:    DO K=1,NUM_COMPS

            KGGD_ROW_NUM = KGGD_ROW_NUM + 1
            IS = STFKEY(KGGD_ROW_NUM)

            IF (IS == 0) THEN                              ! Check for null row in stiffness matrix
               I_KGGD(KGGD_ROW_NUM+1) = I_KGGD(KGGD_ROW_NUM)
               CYCLE k_do
            ENDIF

            NUM_NONZERO_IN_ROW = 0                         ! Form row of non-zero's in arrays RJ, RSTF
j_do1:      DO J=1,NDOFG
               IF (DABS(STF3(IS)%Col_3) >= EPS1) THEN
                  NUM_NONZERO_IN_ROW = NUM_NONZERO_IN_ROW + 1
                  RSTF(NUM_NONZERO_IN_ROW) = STF3(IS)%Col_3
                  RJ(NUM_NONZERO_IN_ROW)   = STF3(IS)%Col_1
               ENDIF
               IS = STF3(IS)%Col_2
               IF (IS == 0) THEN
                  EXIT j_do1
               ENDIF
            ENDDO j_do1

            IF (NUM_NONZERO_IN_ROW > NUM_MAX) THEN
               NUM_MAX = NUM_NONZERO_IN_ROW
            ENDIF
            IF (IS /= 0) THEN
               WRITE(ERR,1625) SUBR_NAME,I
               WRITE(F06,1625) SUBR_NAME,I
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                       ! Coding error, so quit
            ENDIF

            IF (NUM_NONZERO_IN_ROW /= 1) THEN               ! Sort row by the shell method so that RJ is in numerical order
               CALL SORT_INT1_REAL1 ( SUBR_NAME, 'RJ, RSTF', NUM_NONZERO_IN_ROW, RJ, RSTF )
            ENDIF


n_do:       DO N=1,NUM_NONZERO_IN_ROW                      ! Formulate the K-th row of KGGD_II
               IF ((RJ(N) >= KGGD_COL_NUM) .AND. (RJ(N) <= KGGD_COL_NUM+NUM_COMPS-1)) THEN
                  KGGD_II_COL_NUM = RJ(N) - (KGGD_COL_NUM - 1)
                  KGGD_II(K,KGGD_II_COL_NUM) = RSTF(N)
               ENDIF
            ENDDO n_do

j_do3:      DO J=1,NUM_NONZERO_IN_ROW
               KTERM_KGGD = KTERM_KGGD + 1                    ! KTERM_KGGD is a count on the no. records written
               J_KGGD(KTERM_KGGD) = RJ(J)
                 KGGD(KTERM_KGGD) = RSTF(J)
            ENDDO j_do3

            I_KGGD(KGGD_ROW_NUM+1) = I_KGGD(KGGD_ROW_NUM) + NUM_NONZERO_IN_ROW

         ENDDO k_do

         DO K=1,6                                           ! Set lower portion of KGGD_II to be symmetric
            DO J=1,K-1
               KGGD_II(K,J) = KGGD_II(J,K)
            ENDDO
         ENDDO

         CALL COUNTER_PROGRESS(I)
      ENDDO i_do

      WRITE(SC1,*) CR13

      IF (PRTSTIFF(1) >= 1) THEN
         CALL WRITE_SPARSE_CRS ( 'STIFFNESS MATRIX KGGD', 'G ', 'G ', NTERM_KGGD, NDOFG, I_KGGD, J_KGGD, KGGD )
      ENDIF

      WRITE(ERR,101) NUM_MAX
      IF (SUPINFO == 'N') THEN
         WRITE(F06,101) NUM_MAX
      ENDIF



      RETURN

! **********************************************************************************************************************************
  101 FORMAT(' *INFORMATION: MAX NUMBER OF NONZERO TERMS IN A ROW OF THE G-SET STIFFNESS MATRIX     = ',I12,/)

  146 FORMAT(' *INFORMATION: NUMBER OF NONZERO TERMS IN THE KGGD STIFFNESS MATRIX IS                 = ',I12,/)

 1611 FORMAT(' *ERROR  1611: THE G-SET DIFFERENTIAL STIFF MATRIX, KGGD, MUST HAVE SOME NONZERO TERMS. HOWEVER IT HAS ',I12,' TERMS')

 1625 FORMAT(' *ERROR  1625: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' 1ST COL OF ARRAY STF3 INDICATES THERE IS MORE DATA IN ARRAY STF3 FOR ROW ',I12,' OF THE KGGD STIFF'   &
                    ,/,14X,' MATRIX ALTHOUGH THE DOF COUNT IS AT THE END OF THE ROW')

! **********************************************************************************************************************************

      END SUBROUTINE SPARSE_KGGD

   END MODULE STIFFNESS_MATRIX_ASSEMBLY
